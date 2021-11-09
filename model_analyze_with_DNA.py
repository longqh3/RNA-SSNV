# 2021.10.8 #
# 测试完成，一切正常#
# 测试命令 #
# python /home/lqh/Codes/Python/RNA-SSNV/model_analyze_with_DNA.py \
# --step 1 \
# --cancer_type BLCA \
# --DNA_info /home/lqh/Codes/Data/TCGA_maf_files/TCGA-BLCA \
# --RNA_info /home/lqh/Codes/Python/RNA-SSNV/output/BLCA.table \
# --WXS_target_interval /home/lqh/resources/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals_add_chr_to_hg38_rm_alt.bed \
# --exon_interval /home/lqh/resources/database/gencode/GRCh38_GENCODE_v22_exon_rm_alt.bed \
# --RNA_calling_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/tables/info/BLCA_RNA_somatic_calling_info.tsv \
# --RNA_bam_folder /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/BLCA/RNA/apply_BQSR \
# --Mutect2_target_detected_sites /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/BLCA/RNA/RNA_somatic_mutation/VariantsToTable/SNP_WES_Interval_exon.table \
# --project_folder /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results \
# --num_threads 40 \
# --output_file_path /home/lqh/Codes/Python/RNA-SSNV/output/BLCA_DNA_step_1.class

# --template_vcf_file /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer_.training_data_col \

# 测试完成，一切正常#
# python /home/lqh/Codes/Python/RNA-SSNV/model_analyze_with_DNA.py \
# --step 2 \
# --force_call_RNA_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/GBM/RNA/RNA_somatic_mutation/VcfAssembly_new/Mutect2_force_call.txt \
# --instance_path /home/lqh/Codes/Python/RNA-SSNV/output/GBM_DNA_step_1.class \
# --model_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer.model \
# --one_hot_encoder_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer.one_hot_encoder \
# --training_columns_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer.training_data_col \
# --output_file_path /home/lqh/Codes/Python/RNA-SSNV/output/GBM.final.table

# 导入相关需求包
# 基础包
import os
import numpy as np    #导入Python科学计算的基础软件包numpy
import pandas as pd     #导入python的一个数据分析包pandas
# 数据分析包
import matplotlib.pyplot as plt    #导入Python可视化Matplotlib模块的绘图pyplot函数
# 其他包
import pickle # 保存相应模型
import joblib # 其他方式保存模型
import warnings    #导入Python中的warnings模块
import pysam
warnings.filterwarnings('ignore')    #warnings模块利用过滤器来实现忽略告警
# 导入模型解释包
import lime
import lime.lime_tabular
# 导入数据后处理相关包
import os
from multiprocessing import Pool
from multiprocessing import Manager
import portion as P
import vcf

import argparse

# description参数可以用于描述脚本的参数作用，默认为空
parser=argparse.ArgumentParser(description="A discriminate model construction pipeline for RNA-SSNV.")
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.')
# Generic parameter
parser.add_argument('--step', type=int, help='Step indicator of model utilization with DNA.')
parser.add_argument("--output_file_path", help="Path for temp/final output file.")
# Step 1
parser.add_argument('--cancer_type', help='Cancer type for current study.')
parser.add_argument('--DNA_info', help='DNA somatic mutations used to validate.')
parser.add_argument('--RNA_info', help='RNA somatic mutations with model predicting status.')
parser.add_argument('--WXS_target_interval', help='WXS target interval bed file.')
parser.add_argument('--exon_interval', help='GENCODE v22 exon interval bed file.')
parser.add_argument('--RNA_calling_info', help='Tabular info used in RNA somatic mutation calling.')
parser.add_argument('--RNA_bam_folder', help='RNA bam file folder path.')
parser.add_argument('--Mutect2_target_detected_sites', help='Table for Mutect2 detected sites.')
parser.add_argument('--template_vcf_file', default=os.path.join(os.path.split(os.path.realpath(__file__))[0], "resources/template.vcf"), help='Template vcf file to conduct force-calling.')
parser.add_argument('--project_folder', help='Project folder path.')
parser.add_argument('--num_threads', '-nt', type=int, help='Number of threads allowed.')
# Step 2
parser.add_argument('--force_call_RNA_info', help='Force called DNA only mutations in RNA.')
parser.add_argument('--instance_path', help='Data analysis instance path from step 1.')
parser.add_argument("--model_path", help="Path for constructed model.")
parser.add_argument("--one_hot_encoder_path", help="Path for one-hot encoder.")
parser.add_argument("--training_columns_path", help="Path for constructed model.")

args=parser.parse_args()

# Assign allele change rules
ALLELE_CHANGE_DICT = {
    "T>C":"A>G",
    "C>T":"G>A",
    "T>G":"A>C",
    "G>T":"C>A",
    "T>A":"A>T",
    "G>C":"C>G"
}

def RNA_EDIT_process(REDIprotal, DARNED):
    # 获取REDIportal和DARNED数据库中RNA编辑位点信息
    REDIprotal_info = pd.read_table(REDIprotal, header=None)
    DARNED_info = pd.read_table(DARNED, header=None)
    ## 对数据库信息进行预处理，提取所需信息
    REDIprotal_info = REDIprotal_info[[0, 2]]
    DARNED_info = DARNED_info[[0, 2]]
    print(f"REDIportal和DARNED数据库中的RNA编辑位点数目分别为{len(REDIprotal_info)}和{len(DARNED_info)}")
    ## 对REDIpotal和DARNED数据库信息进行合并
    RNA_EDIT_INFO = pd.concat([REDIprotal_info, DARNED_info], ignore_index=True)
    RNA_EDIT_INFO.columns = ["Chromosome", "Start_Position"]
    ## 合并所得数据库信息需要去重
    RNA_EDIT_INFO = RNA_EDIT_INFO.drop_duplicates(keep="first")
    print(f"合并数据库信息并去重后，所得的RNA编辑位点总数为{len(RNA_EDIT_INFO)}")

    return RNA_EDIT_INFO

# 工具函数1：提取pysam对象列表中相应信息
# 根据vcf文件中单个record的突变信息，获取所有特定case的DNA肿瘤样本中var最大AD及其他相应信息，其返回值纳入其他信息的dataframe中
def tumor_bam_record_retrive(df_record, samfile_list):
    try:
        # 获取基本信息
        ref = df_record.Reference_Allele
        alt = df_record.Tumor_Allele1
        chrom = df_record.Chromosome
        position = df_record.Start_Position
        # 获取每个bam文件中对应的碱基coverage信息
        other_coverage_list = []
        # 获取每个bam文件中对应的平均mapping质量信息
        other_median_mapping_quality_list = []

        # 应用首个bam文件中coverage信息完成初始化(设置base的最低质量值为20，与GATK要求保持一致)
        initial_count_coverage_info = count_coverage_decompose(
            samfile_list[0].count_coverage(contig=chrom, start=position - 1, stop=position, quality_threshold=20))
        initial_median_mapping_quality = np.median([x.mapping_quality for x in samfile_list[0].fetch(contig=chrom, start=position - 1, stop=position)])
        # 什么乱七八糟的数据结构啊= =我佛了，还需要强行转为str后才能作为key去访问字典中元素
        optimal_ref_coverage = initial_count_coverage_info[str(ref)]
        optimal_alt_coverage = initial_count_coverage_info[str(alt)]
        optimal_median_mapping_quality = initial_median_mapping_quality
        # 遍历所有bam文件，并进行条件判断
        for samfile in samfile_list[1:]:
            # pysam apply 0-based coordinate system(设置base的最低质量值为20，与GATK要求保持一致)
            count_coverage_info = count_coverage_decompose(
                samfile.count_coverage(contig=chrom, start=position - 1, stop=position, quality_threshold=20))
            # 当前覆盖度信息
            current_ref_coverage = count_coverage_info[str(ref)]
            current_alt_coverage = count_coverage_info[str(alt)]
            current_median_mapping_quality = np.median([x.mapping_quality for x in samfile.fetch(contig=chrom, start=position - 1, stop=position)])
            # 当前覆盖度信息进行对比
            # 比较alt的coverage信息
            if current_alt_coverage > optimal_alt_coverage:
                # 保存当前值，并将更好的值重新赋值
                other_coverage_list.append(str(optimal_ref_coverage) + "/" + str(optimal_alt_coverage))
                other_median_mapping_quality_list.append(str(current_median_mapping_quality))
                optimal_ref_coverage, optimal_alt_coverage, optimal_median_mapping_quality = current_ref_coverage, current_alt_coverage, current_median_mapping_quality
            elif current_alt_coverage == optimal_alt_coverage:
                # 进一步比较ref的coverage信息
                if current_ref_coverage > optimal_ref_coverage:
                    # 保存当前值，并将更好的值重新赋值
                    other_coverage_list.append(str(optimal_ref_coverage) + "/" + str(optimal_alt_coverage))
                    other_median_mapping_quality_list.append(str(current_median_mapping_quality))
                    optimal_ref_coverage, optimal_alt_coverage, optimal_median_mapping_quality = current_ref_coverage, current_alt_coverage, current_median_mapping_quality
                else:
                    # 保存当前值，不进行更新
                    other_coverage_list.append(str(current_ref_coverage) + "/" + str(current_alt_coverage))
                    other_median_mapping_quality_list.append(str(current_median_mapping_quality))
            else:
                # 保存当前值，不进行更新
                other_coverage_list.append(str(current_ref_coverage) + "/" + str(current_alt_coverage))
                other_median_mapping_quality_list.append(str(current_median_mapping_quality))

        return optimal_ref_coverage, optimal_alt_coverage, ";".join(other_coverage_list), optimal_median_mapping_quality, ";".join(other_median_mapping_quality_list)
    except Exception as ex:
        print(ex)
        print("tumor_bam_record_retrive错误！！！")

# 工具函数2：将pysam中count_coverage函数的返回值转换为字典形式
def count_coverage_decompose(count_coverage_info):
    try:
        count_coverage_dict = {}
        # 函数返回值为ACGT的顺序情况
        count_coverage_dict['A'] = count_coverage_info[0][0]
        count_coverage_dict['C'] = count_coverage_info[1][0]
        count_coverage_dict['G'] = count_coverage_info[2][0]
        count_coverage_dict['T'] = count_coverage_info[3][0]
        return count_coverage_dict
    except Exception as ex:
        print(ex)
        print("count_coverage_decompose错误！！！")

# 工具函数10：自定义的error_callback函数
def print_error(self, value):
    print("error: ", value)

# 工具函数10：添加RNA bam信息进入突变位点中
def RNA_bam_info_adder(case_id, gdc_validate_all_info_DNA_only_RNA_missing_case_list,
                       gdc_validate_all_info_DNA_only_RNA_missing_case_only, RNA_calling_info, RNA_bam_folder_loc,
                       tumor_bam_record_retrive):
    print(f"开始根据{case_id}中的突变记录获取对应的RNA测序情况")
    # 获取所有RNA肿瘤aliquots_id信息
    RNA_tumor_aliquots_id = RNA_calling_info.loc[(RNA_calling_info['case_id'] == case_id) & (
            RNA_calling_info['sample_type'] == "Primary Tumor"), 'aliquots_id']
    print(f"{case_id}对应肿瘤RNA测序样本的aliquots_id为{RNA_tumor_aliquots_id}")
    # 获取所有RNA tumor bam文件路径
    RNA_tumor_case_file_paths = RNA_bam_folder_loc + "/" + RNA_tumor_aliquots_id + ".bam"
    # 应用pysam打开所有RNA tumor bam文件，并将对应的pysam对象存放于列表中
    samfile_list = [pysam.AlignmentFile(RNA_tumor_case_file_path, "rb") for RNA_tumor_case_file_path in
                    RNA_tumor_case_file_paths]
    # 对于case内所有突变记录进行分析
    # 新建列以保存相关信息
    gdc_validate_all_info_DNA_only_RNA_missing_case_only['ref_AD_tumor_bam_RNA'] = ""
    gdc_validate_all_info_DNA_only_RNA_missing_case_only['alt_AD_tumor_bam_RNA'] = ""
    gdc_validate_all_info_DNA_only_RNA_missing_case_only['other_AD_tumor_bam_RNA'] = ""
    gdc_validate_all_info_DNA_only_RNA_missing_case_only['median_MQ_tumor_bam_RNA'] = ""
    gdc_validate_all_info_DNA_only_RNA_missing_case_only['other_median_MQ_tumo_bam_RNA'] = ""
    for i in gdc_validate_all_info_DNA_only_RNA_missing_case_only.index:
        # 获得DNA tumor的最佳ref、alt和其他DNA tumor样本的碱基覆盖度信息
        # 将其保存于DataFrame中
        gdc_validate_all_info_DNA_only_RNA_missing_case_only.loc[i, 'ref_AD_tumor_bam_RNA'], \
        gdc_validate_all_info_DNA_only_RNA_missing_case_only.loc[i, 'alt_AD_tumor_bam_RNA'], \
        gdc_validate_all_info_DNA_only_RNA_missing_case_only.loc[i, 'other_AD_tumor_bam_RNA'], \
        gdc_validate_all_info_DNA_only_RNA_missing_case_only.loc[i, 'median_MQ_tumor_bam_RNA'], \
        gdc_validate_all_info_DNA_only_RNA_missing_case_only.loc[
            i, 'other_median_MQ_tumo_bam_RNA'] = tumor_bam_record_retrive(
            gdc_validate_all_info_DNA_only_RNA_missing_case_only.loc[i,], samfile_list)
    # 关闭所有RNA tumor bam文件对应的pysam对象
    [samfile.close() for samfile in samfile_list]
    # 将修改后DataFrame保存于list中
    gdc_validate_all_info_DNA_only_RNA_missing_case_list.append(
        gdc_validate_all_info_DNA_only_RNA_missing_case_only)

# 为联合DNA突变进行分析提供代码支持
class data_prepare():

    # 提取GDC中所有突变的位点信息——Tumor_Sample_UUID Chromosome Start_Position Reference_Allele Tumor_Allele1 Tumor_Allele2
    # 位点层面GDC信息GDC_SNP_info（不包含相应注释，无法用以maftools分析）
    def GDC_site_info_retrieve(self, GDC_info):
        ## 整理GDC相应信息
        print("GDC项目中的突变位点情况为：\n", GDC_info["Variant_Type"].value_counts())
        ### 仅选择SNP突变信息进行分析
        self.GDC_SNP_info = GDC_info[GDC_info['Variant_Type'] == "SNP"]
        ### 重新编制index信息
        self.GDC_SNP_info.reset_index(drop=True, inplace=True)
        ### 将“Tumor_Sample_UUID”列中信息进行切片&提取，使之与RNA体细胞突变中的case_id信息保持一致
        del self.GDC_SNP_info['Tumor_Sample_UUID']
        self.GDC_SNP_info["Tumor_Sample_UUID"] = pd.Series(["-".join(GDC_sample_info.split("-")[0:3])
                                                       for GDC_sample_info in
                                                       self.GDC_SNP_info["Tumor_Sample_Barcode"]])
        ### 仅选取染色体、碱基和case_id信息来进行合并
        self.GDC_SNP_info = self.GDC_SNP_info.loc[:,
                       ["Chromosome", "Start_Position", "Reference_Allele", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_UUID"]]
        ### 重新命名，便于进行合并(Tumor_Allele1对应突变碱基、Tumor_Allele2对应参考碱基)
        self.GDC_SNP_info.columns = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele2", "Tumor_Allele1", "Tumor_Sample_UUID"]
        ### GDC突变去重（由于多个肿瘤样本数据随机组合的缘故，体细胞突变会存在重复），在LUAD项目中，去重前后比例为3:1
        self.GDC_SNP_info = self.GDC_SNP_info.drop_duplicates(keep="first")
        print(f"整理后，GDC项目中的SNP突变数目共计{len(self.GDC_SNP_info)}，对象名称为GDC_SNP_info")

        return self.GDC_SNP_info

    # 获取WXS靶向区域和exon所对应的区间信息（bed格式）
    def WXS_exon_region_dict_generate(self, WXS_target_interval_path, exon_interval_path, num_threads):
        print("开始获取WXS靶向区域对应的区间信息...")
        self.WXS_target_interval = dict(self.IntervalSet_generate_control(
            pd.read_table(WXS_target_interval_path, names=["chr", "start", "end", "target_info", "strand"]), num_threads))
        print("WXS靶向区域对应的区间信息获取完成，对象名为WXS_target_interval")

        print("="*100)

        print("开始获取exon对应的区间信息...")
        self.exon_interval = dict(self.IntervalSet_generate_control(
            pd.read_table(exon_interval_path, names=["chr", "start", "end", "exon_info", "count", "strand"]), num_threads))
        print("exon对应的区间信息获取完成，对象名为exon_interval")

    # 对于target区间相应信息 #
    # 工具函数
    def IntervalSet_generate(self, region_df):
        """Create corresponding IntervalSet for each chromosome.

        Generate dict containing IntervalSet for each chromosome

        Args:
            region_df: a dataframe containing chr, start, end and other info

        Returns:
            A dict containing Interval objects whose key was chromosome name.
        """
        region_interval_dict = {}
        for single_chr in set(region_df["chr"]):
            single_chr_region_df = region_df.loc[region_df["chr"] == single_chr,]
            single_chr_interval_set = P.empty()
            # generate Interval object for each region
            for i in single_chr_region_df.index:
                single_chr_interval_set = single_chr_interval_set | P.closed(single_chr_region_df['start'][i] + 1,
                                                                             single_chr_region_df['end'][i])
                # avoid possible '-' strand start>end issue
                if single_chr_region_df['start'][i] > single_chr_region_df['end'][i]:
                    print(
                        f"'{single_chr_region_df['exon_info'][i]}'s end({single_chr_region_df['end'][i]}) was less than start({single_chr_region_df['start'][i]}) position, need to switch order.")
            # generate IntervalSet object for this chromosome
            region_interval_dict[single_chr] = single_chr_interval_set
        return region_interval_dict

    # 工具函数2：便于进行多进程处理
    def IntervalSet_generate_subprocess(self, region_df, single_chr, region_interval_dict):
        print(f"开始对{single_chr}上的区间信息进行处理......")
        single_chr_region_df = region_df.loc[region_df["chr"] == single_chr,]
        single_chr_interval_set = P.empty()
        # generate Interval object for each region
        for i in single_chr_region_df.index:
            single_chr_interval_set = single_chr_interval_set | P.closed(single_chr_region_df['start'][i] + 1,
                                                                         single_chr_region_df['end'][i])
            # avoid possible '-' strand start>end issue
            if single_chr_region_df['start'][i] > single_chr_region_df['end'][i]:
                print(
                    f"'{single_chr_region_df['exon_info'][i]}'s end({single_chr_region_df['end'][i]}) was less than start({single_chr_region_df['start'][i]}) position, need to switch order.")
        # generate IntervalSet object for this chromosome
        region_interval_dict[single_chr] = single_chr_interval_set
        print(f"{single_chr}上的区间信息处理完成")

    def IntervalSet_generate_control(self, region_df, num_threads):
        print('Parent process %s.' % os.getpid())
        p = Pool(num_threads)
        manager = Manager()
        region_interval_dict = manager.dict()
        for single_chr in set(region_df["chr"]):
            p.apply_async(self.IntervalSet_generate_subprocess, args=(region_df, single_chr, region_interval_dict,))
        print('Waiting for all subprocesses done...')
        p.close()
        p.join()
        print('All subprocesses done.')
        return (region_interval_dict)

    # 获取驱动基因相关信息
    def driver_gene_info_extract(self, driver_gene_loc, cancer_type):
        print("正在读取全部驱动基因信息")
        driver_genes = pd.read_table(driver_gene_loc)
        print(f"开始获取{cancer_type}所对应的驱动基因集合")
        self.cancer_spec_driver_genes = set(driver_genes[driver_genes.CANCER_TYPE==cancer_type]['SYMBOL'])
        print(f"{cancer_type}所对应的驱动基因集合获取完成，共计{len(self.cancer_spec_driver_genes)}个驱动基因，名称为cancer_spec_driver_genes")

    # 获取预测模型相关所有信息，并储存于实例对象当中
    def model_predict_interpret_prepare(self, RF_model_loc, one_hot_encoder_loc, training_col_loc):
        print("开始读取预测模型相关所有信息...")
        # 读取one-hot编码模型、机器学习模型、训练数据列名
        f1, f2 = open(one_hot_encoder_loc, 'rb'), open(RF_model_loc, 'rb')
        self.enc, self.rf_gcv = joblib.load(f1), joblib.load(f2)
        f1.close()
        f2.close()
        self.training_col = pd.read_table(training_col_loc)["0"]

        # print("开始构建模型解释相关信息...")
        # self.shap_explainer = shap.TreeExplainer(self.rf_gcv)

        print("全部构建完成")

    # 有待优化：将实例自带属性进行修改
    # 根据给定one-hot-encoder、机器学习模型，对给定数据进行预测，获得预测信息
    # 输入：待预测数据信息
    # 输出：预测概率值、预测标签
    def model_predict(self, data_info):
        # 编制修改规则
        ALLELE_CHANGE_DICT = {
            "T>C": "A>G",
            "C>T": "G>A",
            "T>G": "A>C",
            "G>T": "C>A",
            "T>A": "A>T",
            "G>C": "C>G"
        }

        # 将index重置为连续状态（避免添加碱基改变标记时，index不一致出现问题）
        self.data_info = data_info.reset_index(drop=True)

        ### 2020.12.7 考虑添加碱基改变信息后，构建模型判断效果改善与否
        data_info_allele_change = self.data_info['Reference_Allele'] + '>' + self.data_info['Tumor_Allele1']
        #### 2021.5.16——修改突变碱基改变信息
        data_info_allele_change = pd.Series([ALLELE_CHANGE_DICT[allele_change] if ALLELE_CHANGE_DICT.keys().__contains__(allele_change) else allele_change for allele_change in data_info_allele_change])
        data_info_allele_change_df = pd.DataFrame(
            self.enc.transform(data_info_allele_change.values.reshape(-1, 1)).toarray(),
            columns=self.enc.categories_[0])
        # 当data_info已有相应碱基改变信息时，删除相应列
        if set(self.enc.categories_[0]).issubset(set(self.data_info.columns)):
            self.data_info.drop(list(self.enc.categories_[0]), axis=1, inplace=True)
        # 新建碱基改变对应列
        self.data_info = pd.concat([self.data_info, data_info_allele_change_df], axis=1)

        # 数据预处理
        # 将类别列和其他无关列（CHROM、POS、REF、ALT、FILTER、CONTQ）从数据集删除，并指定列名顺序（避免预测错误）
        self.data_info_training = self.data_info[self.training_col]

        # 数据概览
        # 检验测试集概率分布情况
        print("\n独立测试数据集概率分布情况如下所示：")
        self.data_info_pred = self.rf_gcv.predict_proba(self.data_info_training)  # 概率形式预测测试集的类别
        plt.hist(self.data_info_pred[:, 1])

        print("---------------------------------------------")
        print("data_info行数为%d" % (len(self.data_info)))
        # 数据预测
        self.data_info_pred_thre = self.rf_gcv.predict(self.data_info_training)  # 概率形式预测测试集的类别
        print(f"预测完成，所得预测结果data_info_pred_thre行数为{len(self.data_info_pred_thre)}")
        # 获取预测阳性、阴性的结果，报告并导出
        self.data_info_positive_thre = self.data_info[self.data_info_pred_thre == 1]
        self.data_info_negative_thre = self.data_info[self.data_info_pred_thre == 0]
        print("总突变数目为%d，其中通过模型判别为阳性的突变位点数目为%d，占比为1:%d" % (
            len(self.data_info), len(self.data_info_positive_thre),
            len(self.data_info) / len(self.data_info_positive_thre)))

        return(self.data_info_pred, self.data_info_pred_thre)

    # 工具函数6：根据给定机器学习模型，对给定数据行进行预测，获得预测信息
    def data_row_interpret(self, data_row):
        shap_explainer = shap.TreeExplainer(self.rf_gcv)
        shap_value = shap_explainer.shap_values(data_row[self.training_col])
        shap.force_plot(shap_explainer.expected_value[1], shap_value[1], features=data_row[self.training_col],
                        feature_names=self.training_col, matplotlib=True)

    # 工具函数7：根据给定机器学习模型，对给定数据集合进行预测，获得预测信息
    def dataset_interpret(self, dataset):
        shap_explainer = shap.TreeExplainer(self.rf_gcv)
        # 计算数据集所对应的所有SHAP值
        self.shap_values = shap_explainer.shap_values(dataset[self.training_col])
        # 为数据集中每个样本绘制其对应特征的SHAP值，这可以更好地理解整体模式，并允许发现预测异常值
        shap.summary_plot(self.shap_values[1], dataset[self.training_col])
        # 绘制简易版本的特征信息汇总图（柱状图）
        shap.summary_plot(self.shap_values[1], dataset[self.training_col], plot_type="bar")

    # 工具函数8：恢复本地数据分析类
    def unpickle(self, instance_path):
        with open(instance_path, 'rb') as f:
            return joblib.load(f)

# 根据所提供的信息来完成RNA与DNA信息的完全整合
# RNA信息：预测为positive、negative的相应RNA体细胞突变，包含需要训练的特征
# DNA信息：输入仅有6列的信息——Tumor_Sample_UUID Chromosome Start_Position Reference_Allele Tumor_Allele1 Tumor_Allele2
# DNA信息需要仅纳入上述列，以保证软件流程的可泛化性
class model_analyze_with_DNA(object):

    # 暂时删除Variant_Classification列，因其在RNA和GDC的注释情况中有极少数不同，但可能会引起合并过程中出现差异
    # 因此只保留一定共有的列 + 各自独有的列信息
    predict_required_columns = ['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Allele1', "Hugo_Symbol",
                        'ref_AD_tumor_RNA', 'alt_AD_tumor_RNA', 'ref_AD_normal', 'alt_AD_normal', 'ref_AD_tumor_DNA', 'alt_AD_tumor_DNA',
                        "Tumor_Sample_UUID", 'record_filter', 'Protein_Change']
    gdc_required_columns = ['Tumor_Sample_UUID', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Allele1', 'Tumor_Allele2']

    # 初始化类实例，positive_info、negative_info、REDIportal_df、gdc_validate_df分别为输入的预测为阳性和阴性的dataframe数据集、REDIportal数据集和相应项目的GDC突变数据集
    def __init__(self, gdc_validate_all_info, RNA_info, WXS_target_interval, exon_interval, num_threads):
        self.gdc_validate_all_info = gdc_validate_all_info[['Tumor_Sample_UUID', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Allele1', 'Tumor_Allele2']]
        self.RNA_info = RNA_info
        self.WXS_target_interval = WXS_target_interval
        self.exon_interval = exon_interval
        self.num_threads = num_threads

    # 关键函数1：根据所得RNA突变来对GDC中突变进行分割与详细解释——已检查，代码无问题，所计算所得P-R与之前一致
    # 输入信息
    # 输入1：具有所有信息的GDC突变——gdc_validate_all_info（已下载的TCGA项目突变数据情况）
    # 输入2：模型预测为positive的RNA突变——positive_RNA_info
    # 输入3：模型预测为negative的RNA突变——negative_RNA_info
    # 输入4：根据区间所对应的IntervalSet_1——WXS靶向区间
    # 输入5：根据区间所对应的IntervalSet_2——外显子区间
    # 输出信息
    # 输出1：DNA与RNA重叠的突变信息（真实突变）——解决
    # 输出2：仅在DNA中出现的部分（不在RNA检测区间内）——解决
    # 输出3：仅在DNA中出现的部分（在RNA检测区间内但被RNA判别为negative，需解释模型判别为negative的原因所在，应用Shap库进行说明）——解决
    # 输出4：仅在DNA中出现的部分（在RNA检测区间内但未被检测到，需给出RNA测序文件中的测序深度等信息，暂时搁置，需要单个函数来说明）——解决，采用关键函数5来说明
    # 输出5：仅在RNA中出现的部分（RNA中新增加的信息量，目前仅将其聚焦于驱动基因上）
    # 分析依据：肿瘤case id和突变位置等信息作为判定是否交叉的信息来源
    def GDC_RNA_discrepancy_analysis_total(self):
        print(f"开始对DNA和RNA-seq中相应突变的交叉、独立情况进行分析")

        print("="*100)

        # 首先对RNA、GDC交叉的case情况作分析，找出是否存在RNA所独有的case，并输出其结果
        GDC_set = set(self.gdc_validate_all_info.Tumor_Sample_UUID)
        RNA_set = set(self.RNA_info.Tumor_Sample_UUID)
        GDC_RNA_intersection_set = set(self.gdc_validate_all_info.Tumor_Sample_UUID).intersection(set(self.RNA_info.Tumor_Sample_UUID))
        print(f"DNA对应case数为{len(GDC_set)}；RNA对应case数为{len(RNA_set)}；DNA、RNA间存在交叉的case数为{len(GDC_RNA_intersection_set)}")

        print(f"其中，缺少DNA证据支持的case信息为{RNA_set - GDC_RNA_intersection_set}，在后续分析中剔除，以避免影响P-R值计算，但需要重点关注")
        self.gdc_validate_all_info = self.gdc_validate_all_info[self.gdc_validate_all_info['Tumor_Sample_UUID'].isin(GDC_RNA_intersection_set)]
        self.RNA_info = self.RNA_info[self.RNA_info['Tumor_Sample_UUID'].isin(GDC_RNA_intersection_set)]
        print(f"最终所得用以和DNA联合分析的RNA体细胞突变数目为{len(self.RNA_info)}")

        # 应用DNA数据对其进行初步评估
        self.RNA_info_TP = pd.merge(self.RNA_info, self.gdc_validate_all_info, on=list(self.gdc_validate_all_info.columns))
        print(f"RNA_info_TP类获取完成，数目为{len(self.RNA_info_TP)}（对象名称为RNA_info_TP）")
        self.RNA_info_TN = self.RNA_info.append(self.RNA_info_TP)
        self.RNA_info_TN = self.RNA_info_TN.drop_duplicates(keep=False)
        print(f"RNA_info_TN类获取完成，数目为{len(self.RNA_info_TN)}（对象名称为RNA_info_TN）")
        # 分别为TP和TN添加标签列
        self.RNA_info_TP['DNA_label'] = 1
        self.RNA_info_TN['DNA_label'] = 0

        print(f"其中，有DNA证据支持的RNA突变数为{len(self.RNA_info_TP)}，缺少DNA证据支持的RNA突变数为{len(self.RNA_info_TN)}")

        # 将TP、TN重新合并（注意此时已经将index重置）
        self.RNA_info = self.RNA_info_TP.append(self.RNA_info_TN, ignore_index=True)
        self.positive_RNA_info = self.RNA_info[self.RNA_info['pred_label'] == 1]
        self.negative_RNA_info = self.RNA_info[self.RNA_info['pred_label'] == 0]


        print("=" * 100)

        self.gdc_validate_all_info = self.gdc_validate_all_info[self.gdc_validate_all_info.Tumor_Sample_UUID.isin(GDC_RNA_intersection_set)]
        # 接下来确认最终纳入比较分析的GDC突变信息
        print(f"纳入比较分析的GDC突变位点数为{len(self.gdc_validate_all_info)}")
        self.gdc_validate_all_info_outside_region = self.gdc_validate_all_info[self.gdc_validate_all_info.apply(lambda x : False if (self.WXS_target_interval[x['Chromosome']].contains(x['Start_Position'])) & (self.exon_interval[x['Chromosome']].contains(x['Start_Position'])) else True, axis=1)]
        self.gdc_validate_all_info = self.gdc_validate_all_info.append(self.gdc_validate_all_info_outside_region)
        self.gdc_validate_all_info = self.gdc_validate_all_info.drop_duplicates(keep=False)
        print(f"不位于RNA研究区间内的的GDC突变集gdc_validate_all_info_outside_region的数目为{len(self.gdc_validate_all_info_outside_region)}，最终纳入比较分析的位于RNA研究区间内的GDC突变位点数为{len(self.gdc_validate_all_info)}")

        print("=" * 100)

        # 最后确认最终纳入比较分析的RNA突变信息
        print(f"纳入比较分析的RNA positive突变位点数为{len(self.positive_RNA_info)}、negative突变位点数为{len(self.negative_RNA_info)}")
        RNA_only_set = set(self.positive_RNA_info.Tumor_Sample_UUID) - GDC_RNA_intersection_set
        print(f"RNA突变集case数共为{len(set(self.positive_RNA_info.Tumor_Sample_UUID))}，其中无GDC突变证据支持的RNA突变集case为{RNA_only_set}")
        self.positive_RNA_info_RNA_case_only = self.positive_RNA_info[self.positive_RNA_info.Tumor_Sample_UUID.isin(RNA_only_set)]
        self.positive_RNA_info = self.positive_RNA_info[~self.positive_RNA_info.Tumor_Sample_UUID.isin(RNA_only_set)]
        self.negative_RNA_info = self.negative_RNA_info[~self.negative_RNA_info.Tumor_Sample_UUID.isin(RNA_only_set)]
        print(f"无GDC证据支持的RNA positive突变集positive_RNA_info_RNA_case_only的数目为{len(self.positive_RNA_info_RNA_case_only)}，因此，最终纳入比较分析的RNA positive突变位点数为{len(self.positive_RNA_info)}、negative突变位点数为{len(self.negative_RNA_info)}")

        print("=" * 100)

        # 首先求DNA与RNA重叠的突变信息
        # 不采用"Hugo_Symbol"求交集，避免可能出现的遗漏/不一致
        self.positive_RNA_info_GDC_intersect = self.positive_RNA_info.merge(self.gdc_validate_all_info, how = 'inner')
        positive_RNA_info_GDC_intersect_case_set = set(self.positive_RNA_info_GDC_intersect.Tumor_Sample_UUID.value_counts().index)
        print(f"RNA、GDC中共有的突变positive_RNA_info_GDC_intersect数目为{len(self.positive_RNA_info_GDC_intersect)}，分布在{len(positive_RNA_info_GDC_intersect_case_set)}个case中")

        print("="*100)

        # 接下来求仅在DNA中出现的突变部分
        self.gdc_validate_all_info_DNA_only = self.gdc_validate_all_info.append(self.positive_RNA_info_GDC_intersect)
        self.gdc_validate_all_info_DNA_only = self.gdc_validate_all_info_DNA_only.drop_duplicates(subset=['Tumor_Sample_UUID', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Allele1', 'Tumor_Allele2'], keep=False)
        self.gdc_validate_all_info_DNA_only = self.gdc_validate_all_info_DNA_only.dropna(axis='columns', how="all") # 去除所有值均为空的列（RNA中所有注释列与特征列）
        print(f"仅在DNA中出现的突变gdc_validate_all_info_DNA_only数目为{len(self.gdc_validate_all_info_DNA_only)}，分布在{len(set(self.gdc_validate_all_info_DNA_only.Tumor_Sample_UUID))}个case中")
        # 进一步对仅在DNA中出现的突变部分进行细分
        print("对仅在DNA中出现的突变进一步细分")
        # 与RNA negative突变集部分求交集，以确认其在RNA中误判为negative的情况
        self.gdc_validate_all_info_DNA_only_RNA_negative = self.gdc_validate_all_info_DNA_only.merge(self.negative_RNA_info, how = 'inner')
        print(f"在RNA中误判为negative的突变gdc_validate_all_info_DNA_only_RNA_negative数目为{len(self.gdc_validate_all_info_DNA_only_RNA_negative)}，分布在{len(set(self.gdc_validate_all_info_DNA_only_RNA_negative.Tumor_Sample_UUID))}个case中")
        # 剩余部分GDC突变则属于既不在RNA positive，也不在negative部分的突变，自然，其属于未被RNA检测到的突变（需要引入其他工具函数来对相应位点在RNA内测序深度进行探索、汇报）
        self.gdc_validate_all_info_DNA_only_RNA_missing = self.gdc_validate_all_info_DNA_only.append(self.gdc_validate_all_info_DNA_only_RNA_negative)
        self.gdc_validate_all_info_DNA_only_RNA_missing = self.gdc_validate_all_info_DNA_only_RNA_missing.drop_duplicates(subset=['Tumor_Sample_UUID', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Allele1', 'Tumor_Allele2'], keep=False)
        self.gdc_validate_all_info_DNA_only_RNA_missing = self.gdc_validate_all_info_DNA_only_RNA_missing.dropna(axis='columns', how="all")  # 去除所有值均为空的列（RNA中所有注释列与特征列）
        print(f"未被RNA检测到的突变gdc_validate_all_info_DNA_only_RNA_missing数目为{len(self.gdc_validate_all_info_DNA_only_RNA_missing)}，分布在{len(set(self.gdc_validate_all_info_DNA_only_RNA_missing.Tumor_Sample_UUID))}个case中")

        print("=" * 100)

        # 最后求仅在RNA中出现的突变部分
        ## 阳性突变
        self.positive_RNA_info_RNA_only = self.positive_RNA_info.append(self.positive_RNA_info_GDC_intersect)
        self.positive_RNA_info_RNA_only = self.positive_RNA_info_RNA_only.drop_duplicates(keep=False)
        ## 阴性突变
        self.negative_RNA_info_RNA_only = self.negative_RNA_info.append(self.gdc_validate_all_info_DNA_only_RNA_negative)
        self.negative_RNA_info_RNA_only = self.negative_RNA_info_RNA_only.drop_duplicates(keep=False)

        print(f"仅在RNA中出现的阳性突变positive_RNA_info_RNA_only数目为{len(self.positive_RNA_info_RNA_only)}，分布在{len(set(self.positive_RNA_info_RNA_only.Tumor_Sample_UUID))}个case中；"
              f"阴性突变negative_RNA_info_RNA_only数目为{len(self.negative_RNA_info_RNA_only)}，分布在{len(set(self.negative_RNA_info_RNA_only.Tumor_Sample_UUID))}个case中")

    # 关键函数2：检查GDC中未在RNA内检测到突变的情况，确定其在RNA中是否有表达/表达情况如何（为相应GDC突变添加bam文件信息）
    # GDC突变信息内被添加了ref_AD_tumor_bam_RNA、alt_AD_tumor_bam_RNA、other_AD_tumor_bam_RNA、median_MQ_tumor_bam_RNA、other_median_MQ_tumo_bam_RNA相应列信息——gdc_validate_all_info_DNA_only_RNA_missing_RNA_info
    def GDC_RNA_missing_check(self, RNA_calling_info, RNA_bam_folder_loc, Mutect2_target_detected_sites):
        print(f"未被RNA检测到的GDC突变数为{len(self.gdc_validate_all_info_DNA_only_RNA_missing)}")
        print("开始对相应突变在原始RNA测序文件中的覆盖度/测序深度进行分析......")

        print("="*100)

        # 新建列表以保存获取信息后的DataFrame
        # 此处转入多进程并行处理（2021.10.13）
        manager = Manager()
        gdc_validate_all_info_DNA_only_RNA_missing_case_list = manager.list()

        print('Parent process %s.' % os.getpid())
        p = Pool(self.num_threads)
        for case_id in self.gdc_validate_all_info_DNA_only_RNA_missing.Tumor_Sample_UUID.value_counts().index:
            p.apply_async(func=RNA_bam_info_adder, args=(case_id, gdc_validate_all_info_DNA_only_RNA_missing_case_list, self.gdc_validate_all_info_DNA_only_RNA_missing[self.gdc_validate_all_info_DNA_only_RNA_missing.Tumor_Sample_UUID == case_id], RNA_calling_info, RNA_bam_folder_loc, tumor_bam_record_retrive,), error_callback=print_error)
        print('Waiting for all subprocesses done...')
        p.close()
        p.join()
        print('All subprocesses done.')

        gdc_validate_all_info_DNA_only_RNA_missing_case_list = list(gdc_validate_all_info_DNA_only_RNA_missing_case_list)
        self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info = pd.concat(gdc_validate_all_info_DNA_only_RNA_missing_case_list)
        print(f"在RNA中缺少突变信息的GDC突变所对应的RNA测序深度信息获取完成，对象名为gdc_validate_all_info_DNA_only_RNA_missing_RNA_info")

        # 归因0：Mutect2中已有记录，但未出现在最终positive或negative突变集中（位点被上述方式滤过）
        self.Mutect2_target_detected_sites = Mutect2_target_detected_sites
        gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_temp = self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info.copy(deep=True)
        self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_multi_filtered_sites = pd.merge(gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_temp, self.Mutect2_target_detected_sites)
        print(f"经RNA预处理后所缺失的突变数为{len(self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_multi_filtered_sites)}，"
              f"对象名为gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_multi_filtered_sites")

        # 归因1：Mutect2中无相应记录，需要对其进行force-call以获取相应信息并进行判别
        print(f"添加原始RNA测序文件中的覆盖度/测序深度信息后，需要提取位点相应信息，并应用Mutect2 force-calling来进行检查的突变数为{len(self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info)}")

    # 暂时搁置？主要是仅由RNA bam文件的coverage信息来判定原因过于单薄
    # 关键函数3：进一步分析强制call突变后所得信息的基本情况
    def GDC_Mutect2_missing_analysis(self, GDC_Mutect2_missing_force_call_RNA_info):
        # 去除Funcotator中因重复注释而引入的重复记录信息
        self.GDC_Mutect2_missing_force_call_Mutect2_info = GDC_Mutect2_missing_force_call_RNA_info.drop_duplicates(subset=['Chromosome', 'Start_Position', 'Tumor_Sample_UUID'], keep='first', inplace=False)
        # 检查突变集数目基本情况
        print(f"经过Mutect2 Force-calling检查的、已提取信息的突变集数目为{len(self.GDC_Mutect2_missing_force_call_Mutect2_info)}")
        self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_force_call_failed = self.GDC_Mutect2_missing_force_call_Mutect2_info.append(self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info)
        self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_force_call_failed = self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_force_call_failed.drop_duplicates(
            subset=['Tumor_Sample_UUID', 'Chromosome', 'Start_Position'], keep=False)
        print(f"其中仍未被force-call检查的突变位点gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_force_call_failed数目为{len(self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_force_call_failed)}个")
        # 将force call后包含注释+特征信息的内容补充添加至原有GDC列（添加RNA bam信息）中
        ## 删除相应ref、alt碱基信息以完成合并
        ## 自然因merge的特性，而排除掉了failed force call的三个位点对应情况
        self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info.drop(['Reference_Allele', 'Tumor_Allele1', 'Tumor_Allele2'], axis=1, inplace=True)
        self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info = self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info.merge(self.GDC_Mutect2_missing_force_call_Mutect2_info)

        print(f"最终，添加force call内注释+特征信息后的原有GDC列（添加RNA bam信息）构建完成，对象名为gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info，数目为{len(self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info)}")

        print("="*100)

        print("开始进行归因分析，并准备对概率值进行预测")

        temp_gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info = self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info.copy(deep=True)

        # 归因0：无法获取相应信息（non-SNP情况）
        self.GDC_Mutect2_missing_force_call_RNA_info_non_SNP = temp_gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info[temp_gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info.apply(lambda x : True if (len(x['Reference_Allele'])>1)or(len(x['Tumor_Allele1'])>1) else False, axis=1)]
        print(f"{len(self.GDC_Mutect2_missing_force_call_RNA_info_non_SNP)}个突变在Mutect2中被识别为non-SNP突变，无法提取相应信息，对象名称为GDC_Mutect2_missing_force_call_RNA_info_non_SNP")
        temp_gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info = temp_gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info[temp_gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info.apply(lambda x : False if (len(x['Reference_Allele'])>1)or(len(x['Tumor_Allele1'])>1) else True, axis=1)]

        # 归因1：预测概率值过低（成为体细胞突变的证据不足）（确认其判别后是否全部为阴性位点）
        self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info_final = temp_gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info
        print(f"已排除所有无法进行模型预测的突变，最终剩余的突变gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info_final数目为{len(self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info_final)}")

    # 关键函数4：总结函数，为上述所有突变类型分类添加注释信息（Tag、Sub_Tag、Predict_Prob、Predict_Tag）后完成合并
    # 主标签：Tag——RNA_DNA_overlap, RNA_only, DNA_only
    # 子标签：Sub_Tag——non_SNP, force_called, force_call_failed
    # 概率值：Predicted_prob；预测标签：Predicted_label
    # 注意：检查并确认相应数目是否一致，保证不遗漏任何一个突变位点
    def GDC_RNA_info_summary(self, data_prepare_demo):
        ### 检查各部分基本信息，添加tag和sub_tag注释 ###
        print("="*100)
        # 对于DNA和RNA重叠部分
        print(f"首先完成DNA和RNA重叠部分的合并，其中positive对象为positive_RNA_info_GDC_intersect，数目为{len(self.positive_RNA_info_GDC_intersect)}；\n"
              f"negative对象为gdc_validate_all_info_DNA_only_RNA_negative，数目为{len(self.gdc_validate_all_info_DNA_only_RNA_negative)}")
        RNA_DNA_overlap = pd.concat([self.positive_RNA_info_GDC_intersect, self.gdc_validate_all_info_DNA_only_RNA_negative])
        RNA_DNA_overlap['Tag'] = "RNA_DNA_overlap"

        # 对于仅在RNA中出现部分
        print(f"\n接下来完成仅在RNA中出现部分的合并，其中positive对象为positive_RNA_info_RNA_only，数目为{len(self.positive_RNA_info_RNA_only)}；\n"
              f"negative对象为negative_RNA_info_RNA_only，数目为{len(self.negative_RNA_info_RNA_only)}")
        RNA_only = pd.concat([self.positive_RNA_info_RNA_only, self.negative_RNA_info_RNA_only])
        RNA_only['Tag'] = "RNA_only"

        # 上述两部分均具有突变注释 + 模型预测相关特征信息，可直接应用于模型预测过程中
        # 下述部分中仅有force_called子类型可直接应用于模型预测过程中，其他子类型无法进行预测

        # 对于仅在DNA中出现部分
        print(f"\n最后完成仅在DNA中出现部分的合并，可分为多个子部分，包括\n"
              f"non-SNP对象为GDC_Mutect2_missing_force_call_RNA_info_non_SNP，数目为{len(self.GDC_Mutect2_missing_force_call_RNA_info_non_SNP)}\n"
              f"force_call_failed对象为gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_force_call_failed，数目为{len(self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_force_call_failed)}\n"
              f"force_called对象为gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info_final，数目为{len(self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info_final)}\n")
        self.GDC_Mutect2_missing_force_call_RNA_info_non_SNP['Sub_Tag'] = "non_SNP"
        self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_force_call_failed['Sub_Tag'] = "force_call_failed"
        self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info_final['Sub_Tag'] = "force_called"
        DNA_only = pd.concat([self.GDC_Mutect2_missing_force_call_RNA_info_non_SNP, self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_force_call_failed, self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info_final])
        DNA_only['Tag'] = "DNA_only"

        ### 增加Predicted_prob和Predicted_label列信息 ###
        # 首先取得可进行模型预测的突变信息并获取预测概率值、预测标签
        predictable_part = pd.concat([RNA_DNA_overlap, RNA_only, DNA_only[DNA_only['Sub_Tag']=="force_called"]])
        pred_prob, pred_label = data_prepare_demo.model_predict(predictable_part)
        predictable_part['pred_prob'] = pred_prob[:, 1]
        predictable_part['pred_label'] = pred_label
        # 接下来将其他突变合并进来
        self.final_info = pd.concat([predictable_part, DNA_only[DNA_only['Sub_Tag']!="force_called"]])
        self.final_info.reset_index(drop=True, inplace=True)

        print(f"评估结束，对于所有RNA突变，其precision为{len(self.positive_RNA_info_GDC_intersect)/(len(self.positive_RNA_info_GDC_intersect)+len(self.positive_RNA_info_RNA_only))}，recall为{len(self.positive_RNA_info_GDC_intersect)/(len(self.positive_RNA_info_GDC_intersect)+len(self.gdc_validate_all_info_DNA_only_RNA_negative))}")

        print(f"\n汇总DNA、RNA所有突变信息的汇总表final_info构建完成！包含{len(self.final_info)}个突变对应信息")
        print("主标签：Tag——RNA_DNA_overlap, RNA_only, DNA_only")
        print("子标签：Sub_Tag——non_SNP, force_called, force_call_failed")
        print("概率值：pred_prob；预测标签：pred_label")

    # 工具函数1：提取pysam对象列表中相应信息
    # 根据vcf文件中单个record的突变信息，获取所有特定case的DNA肿瘤样本中var最大AD及其他相应信息，其返回值纳入其他信息的dataframe中
    def tumor_bam_record_retrive(self, df_record, samfile_list):
        try:
            # 获取基本信息
            ref = df_record.Reference_Allele
            alt = df_record.Tumor_Allele1
            chrom = df_record.Chromosome
            position = df_record.Start_Position
            # 获取每个bam文件中对应的碱基coverage信息
            other_coverage_list = []
            # 获取每个bam文件中对应的平均mapping质量信息
            other_median_mapping_quality_list = []

            # 应用首个bam文件中coverage信息完成初始化(设置base的最低质量值为20，与GATK要求保持一致)
            initial_count_coverage_info = self.count_coverage_decompose(
                samfile_list[0].count_coverage(contig=chrom, start=position - 1, stop=position, quality_threshold=20))
            initial_median_mapping_quality = np.median([x.mapping_quality for x in samfile_list[0].fetch(contig=chrom, start=position - 1, stop=position)])
            # 什么乱七八糟的数据结构啊= =我佛了，还需要强行转为str后才能作为key去访问字典中元素
            optimal_ref_coverage = initial_count_coverage_info[str(ref)]
            optimal_alt_coverage = initial_count_coverage_info[str(alt)]
            optimal_median_mapping_quality = initial_median_mapping_quality
            # 遍历所有bam文件，并进行条件判断
            for samfile in samfile_list[1:]:
                # pysam apply 0-based coordinate system(设置base的最低质量值为20，与GATK要求保持一致)
                count_coverage_info = self.count_coverage_decompose(
                    samfile.count_coverage(contig=chrom, start=position - 1, stop=position, quality_threshold=20))
                # 当前覆盖度信息
                current_ref_coverage = count_coverage_info[str(ref)]
                current_alt_coverage = count_coverage_info[str(alt)]
                current_median_mapping_quality = np.median([x.mapping_quality for x in samfile.fetch(contig=chrom, start=position - 1, stop=position)])
                # 当前覆盖度信息进行对比
                # 比较alt的coverage信息
                if current_alt_coverage > optimal_alt_coverage:
                    # 保存当前值，并将更好的值重新赋值
                    other_coverage_list.append(str(optimal_ref_coverage) + "/" + str(optimal_alt_coverage))
                    other_median_mapping_quality_list.append(str(current_median_mapping_quality))
                    optimal_ref_coverage, optimal_alt_coverage, optimal_median_mapping_quality = current_ref_coverage, current_alt_coverage, current_median_mapping_quality
                elif current_alt_coverage == optimal_alt_coverage:
                    # 进一步比较ref的coverage信息
                    if current_ref_coverage > optimal_ref_coverage:
                        # 保存当前值，并将更好的值重新赋值
                        other_coverage_list.append(str(optimal_ref_coverage) + "/" + str(optimal_alt_coverage))
                        other_median_mapping_quality_list.append(str(current_median_mapping_quality))
                        optimal_ref_coverage, optimal_alt_coverage, optimal_median_mapping_quality = current_ref_coverage, current_alt_coverage, current_median_mapping_quality
                    else:
                        # 保存当前值，不进行更新
                        other_coverage_list.append(str(current_ref_coverage) + "/" + str(current_alt_coverage))
                        other_median_mapping_quality_list.append(str(current_median_mapping_quality))
                else:
                    # 保存当前值，不进行更新
                    other_coverage_list.append(str(current_ref_coverage) + "/" + str(current_alt_coverage))
                    other_median_mapping_quality_list.append(str(current_median_mapping_quality))

            return optimal_ref_coverage, optimal_alt_coverage, ";".join(other_coverage_list), optimal_median_mapping_quality, ";".join(other_median_mapping_quality_list)
        except Exception as ex:
            print(ex)
            print("tumor_bam_record_retrive错误！！！")

    # 工具函数2：将pysam中count_coverage函数的返回值转换为字典形式
    def count_coverage_decompose(self, count_coverage_info):
        try:
            count_coverage_dict = {}
            # 函数返回值为ACGT的顺序情况
            count_coverage_dict['A'] = count_coverage_info[0][0]
            count_coverage_dict['C'] = count_coverage_info[1][0]
            count_coverage_dict['G'] = count_coverage_info[2][0]
            count_coverage_dict['T'] = count_coverage_info[3][0]
            return count_coverage_dict
        except Exception as ex:
            print(ex)
            print("count_coverage_decompose错误！！！")

    # 工具函数3：针对性转换整体突变集（不在乎分类信息，仅提取chr、start、end信息，用以Mutect2强制call突变以检查原因）
    # 输入：整体突变集情况（maf格式等）
    # 输入2：指定输出文件夹
    # 输出1：各个case所对应待检查的突变数目信息（并将其保存于输出文件夹中）
    def interval_info_extractor(self, total_mutation_set, output_folder):
        print(f"开始将给定突变集({len(total_mutation_set)}个突变)中信息提取为单个case的interval信息")
        # 初始化保存case id的list
        result_case_list = []
        # 遍历每个case，获取信息
        for case_id in total_mutation_set.Tumor_Sample_UUID.value_counts().index:
            # 求单个case的子集
            case_mutation_set = total_mutation_set[total_mutation_set.Tumor_Sample_UUID==case_id]
            # 提取相应信息并保存于输出文件夹中
            print(f"{case_id}中有{len(case_mutation_set)}个突变信息待提取...")
            case_mutation_set[['Chromosome', 'Start_Position', 'Start_Position']].to_csv(os.path.join(output_folder, f"{case_id}.interval"), index=False, sep="\t", header=False)
            print(f"信息提取完成，interval信息保存于{os.path.join(output_folder)}")
        # 将所有提取后的case_id信息同样保存于相同文件夹内，便于进行case确定
        pd.Series(list(total_mutation_set.Tumor_Sample_UUID.value_counts().index)).to_csv(os.path.join(output_folder, f"case_info.table"), index=False)

    # 工具函数4：针对性转换整体突变集（不在乎分类信息，仅提取chr、start、end信息，用以Mutect2强制call突变以检查原因）
    # 此处输出信息为模板化的vcf文件
    # 输入：整体突变集情况（maf格式等）
    # 输入2：指定输出文件夹
    # 输出1：各个case所对应待检查的突变数目信息（并将其保存于输出文件夹中）
    def vcf_info_extractor(self, total_mutation_set, template_vcf_file, output_folder):
        # 检查文件夹是否存在
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        # 获取vcf读取后模板情况
        vcf_reader = vcf.Reader(open(template_vcf_file, 'r'))
        record_list = [record for record in vcf_reader]
        # 开始信息提取过程
        print(f"开始将给定突变集({len(total_mutation_set)}个突变)中信息提取为单个case的vcf信息")
        # 遍历每个case，获取信息
        for case_id in total_mutation_set.Tumor_Sample_UUID.value_counts().index:
            # 求单个case的子集
            case_mutation_set = total_mutation_set[total_mutation_set.Tumor_Sample_UUID==case_id]
            # 提取相应信息并保存于输出文件夹中
            print(f"{case_id}中有{len(case_mutation_set)}个突变信息待提取...")
            # 新建vcf文件以保存相应信息
            vcf_writer = vcf.Writer(open(os.path.join(output_folder, f"{case_id}.vcf"), 'w'), vcf_reader)
            for i in case_mutation_set.index:
                record_list[0].CHROM = case_mutation_set.loc[i, "Chromosome"]
                record_list[0].POS = case_mutation_set.loc[i, "Start_Position"]
                record_list[0].REF = case_mutation_set.loc[i, "Reference_Allele"]
                record_list[0].ALT = list(case_mutation_set.loc[i, "Tumor_Allele1"])
                vcf_writer.write_record(record_list[0])
            vcf_writer.close()
            print(f"信息提取完成，将vcf信息保存于{os.path.join(output_folder)}")
        # 将所有提取后的case_id信息同样保存于相同文件夹内，便于进行case确定
        pd.Series(list(total_mutation_set.Tumor_Sample_UUID.value_counts().index)).to_csv(os.path.join(output_folder, f"case_info.table"), index=False)

    # 似乎也没有必要自己来实现excel的功能，暂且搁置
    # 工具函数7：根据给定条件（Hugo_Symbol、Tumor_Sample_UUID、Start_Position等）来进行筛选
    # 目前暂时通过输入字典的方式，来进行maf信息的筛选分析
    # 应用map方法，针对不同提示信息完成筛选
    # Chromosome、Start_Position、Hugo_Symbol、Tumor_Sample_UUID、Tag、Sub_Tag、Predicted_prob、Predicted_label均为潜在筛选标准
    def multiple_query(self, data_info, query_terms):
        data_bool = True
        for col in query_terms.keys():
            data_bool = data_info[col].map(lambda x: query_terms[col].__contains__(x)) & data_bool
        return data_info[data_bool]

    # 工具函数8：将所得结果导出成为可为maftools处理的maf格式文件
    def result_2_standard_maf(self, result_info, output_maf_loc):
        print("开始对结果文件进行处理使之成为标准的maf格式，添加End_Position、Variant_Type、Tumor_Sample_Barcode列信息")
        # result_info['End_Position'] = result_info['Start_Position']
        result_info['Variant_Type'] = "SNP"
        result_info['Tumor_Sample_Barcode'] = result_info['Tumor_Sample_UUID']
        result_info.insert(2, 'End_Position', result_info['Start_Position']) # 避免KGGSeq因空值而出现读取偏离的情况
        print("将标准maf格式文件导出...")
        result_info.to_csv(output_maf_loc, sep="\t", index=False)
        print("标准maf格式文件导出成功！")
        result_info.drop(['End_Position', 'Variant_Type', 'Tumor_Sample_Barcode'],axis=1,inplace=True)

    # 工具函数9：将数据分析类保存至本地文件
    def pickle(self, output_file_path):
        with open(output_file_path, 'wb') as f:
            joblib.dump(self, f, compress=3)

    # 工具函数10：自定义的error_callback函数
    def print_error(self, value):
        print("error: ", value)

    # 工具函数10：添加RNA bam信息进入突变位点中
    def RNA_bam_info_adder(self, case_id, gdc_validate_all_info_DNA_only_RNA_missing_case_list,
                           gdc_validate_all_info_DNA_only_RNA_missing_case_only, RNA_calling_info, RNA_bam_folder_loc,
                           tumor_bam_record_retrive):
        print(f"开始根据{case_id}中的突变记录获取对应的RNA测序情况")
        # 获取所有RNA肿瘤aliquots_id信息
        RNA_tumor_aliquots_id = RNA_calling_info.loc[(RNA_calling_info['case_id'] == case_id) & (
                RNA_calling_info['sample_type'] == "Primary Tumor"), 'aliquots_id']
        print(f"{case_id}对应肿瘤RNA测序样本的aliquots_id为{RNA_tumor_aliquots_id}")
        # 获取所有RNA tumor bam文件路径
        RNA_tumor_case_file_paths = RNA_bam_folder_loc + "/" + RNA_tumor_aliquots_id + ".bam"
        # 应用pysam打开所有RNA tumor bam文件，并将对应的pysam对象存放于列表中
        samfile_list = [pysam.AlignmentFile(RNA_tumor_case_file_path, "rb") for RNA_tumor_case_file_path in
                        RNA_tumor_case_file_paths]
        # 对于case内所有突变记录进行分析
        # 新建列以保存相关信息
        gdc_validate_all_info_DNA_only_RNA_missing_case_only['ref_AD_tumor_bam_RNA'] = ""
        gdc_validate_all_info_DNA_only_RNA_missing_case_only['alt_AD_tumor_bam_RNA'] = ""
        gdc_validate_all_info_DNA_only_RNA_missing_case_only['other_AD_tumor_bam_RNA'] = ""
        gdc_validate_all_info_DNA_only_RNA_missing_case_only['median_MQ_tumor_bam_RNA'] = ""
        gdc_validate_all_info_DNA_only_RNA_missing_case_only['other_median_MQ_tumo_bam_RNA'] = ""
        for i in gdc_validate_all_info_DNA_only_RNA_missing_case_only.index:
            # 获得DNA tumor的最佳ref、alt和其他DNA tumor样本的碱基覆盖度信息
            # 将其保存于DataFrame中
            gdc_validate_all_info_DNA_only_RNA_missing_case_only.loc[i, 'ref_AD_tumor_bam_RNA'], \
            gdc_validate_all_info_DNA_only_RNA_missing_case_only.loc[i, 'alt_AD_tumor_bam_RNA'], \
            gdc_validate_all_info_DNA_only_RNA_missing_case_only.loc[i, 'other_AD_tumor_bam_RNA'], \
            gdc_validate_all_info_DNA_only_RNA_missing_case_only.loc[i, 'median_MQ_tumor_bam_RNA'], \
            gdc_validate_all_info_DNA_only_RNA_missing_case_only.loc[
                i, 'other_median_MQ_tumo_bam_RNA'] = tumor_bam_record_retrive(
                gdc_validate_all_info_DNA_only_RNA_missing_case_only.loc[i,], samfile_list)
        # 关闭所有RNA tumor bam文件对应的pysam对象
        [samfile.close() for samfile in samfile_list]
        # 将修改后DataFrame保存于list中
        gdc_validate_all_info_DNA_only_RNA_missing_case_list.append(
            gdc_validate_all_info_DNA_only_RNA_missing_case_only)


if __name__ == '__main__':

    if args.step == 1:
        # 预处理 #
        # 读取DNA位点信息
        DNA_info = pd.read_table(args.DNA_info)
        # 读取经过模型预测的RNA位点信息
        RNA_info = pd.read_table(args.RNA_info)
        # 读取portion相关区间信息
        data_prepare_demo = data_prepare()
        WXS_target_interval_path = args.WXS_target_interval
        exon_interval_path = args.exon_interval
        data_prepare_demo.WXS_exon_region_dict_generate(WXS_target_interval_path, exon_interval_path, args.num_threads)

        # 额外处理（仅对GDC对应maf文件进行）
        DNA_info = data_prepare_demo.GDC_site_info_retrieve(DNA_info)

        # 正式分析 #
        # 初始化分析类
        result_analysis_demo = model_analyze_with_DNA(DNA_info,
                                                   RNA_info,
                                                   data_prepare_demo.WXS_target_interval,
                                                   data_prepare_demo.exon_interval,
                                                   args.num_threads)

        # 对DNA和RNA中相应突变的交叉、独立情况进行分析
        result_analysis_demo.GDC_RNA_discrepancy_analysis_total()

        # 检查RNA完全未检测到的GDC突变
        # 确认其在RNA中表达/覆盖度情况
        RNA_calling_info = pd.read_table(args.RNA_calling_info)
        RNA_bam_folder_loc = args.RNA_bam_folder
        Mutect2_target_detected_sites = pd.read_table(args.Mutect2_target_detected_sites)
        result_analysis_demo.GDC_RNA_missing_check(RNA_calling_info, RNA_bam_folder_loc, Mutect2_target_detected_sites)
        # 进一步导出需要force-call的vcf文件信息
        result_analysis_demo.vcf_info_extractor(result_analysis_demo.gdc_validate_all_info_DNA_only_RNA_missing,
                                                args.template_vcf_file,
                                                os.path.join(args.project_folder, args.cancer_type, "RNA/RNA_somatic_mutation/MAFToVCF/DNA_only_RNA_missing_Mutect2_check"))
        # 最后将其信息一次导出
        # result_analysis_demo.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info.to_csv(args.output_table_path, sep="\t", index=False)
        result_analysis_demo.pickle(args.output_file_path)

    if args.step == 2:
        # 预处理 #
        # 读取force-call后数据文件
        DNA_only_RNA_missing_RNA_info_Mutect2_force_call = pd.read_table(args.force_call_RNA_info)
        # 读取数据分析类
        data_prepare_demo = data_prepare()
        result_analysis_demo = data_prepare_demo.unpickle(args.instance_path)
        # 读取实际预测模型
        data_prepare_demo.model_predict_interpret_prepare(args.model_path, args.one_hot_encoder_path, args.training_columns_path)

        # 正式分析 #
        # 添加force-call后信息
        result_analysis_demo.GDC_Mutect2_missing_analysis(DNA_only_RNA_missing_RNA_info_Mutect2_force_call)
        # 获取最终的大table信息！
        result_analysis_demo.GDC_RNA_info_summary(data_prepare_demo)
        # 导出最终的大table信息
        result_analysis_demo.final_info.to_csv(args.output_file_path, sep="\t", index=False)