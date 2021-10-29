# python /home/lqh/Codes/Python/RNA-SSNV/lib/own_data_vcf_info_retriver.py \
# --cancer_type BLCA \
# --RNA_calling_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/tables/info/BLCA_RNA_somatic_calling_info.tsv \
# --project_folder /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results \
# --exon_interval /home/lqh/resources/database/gencode/GRCh38_GENCODE_v22_exon_rm_alt.bed \
# --output_table_path /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/BLCA/RNA/RNA_somatic_mutation/VcfAssembly_new/SNP_WES_Interval_exon.txt \
# --num_threads 60

# 特征提取核心代码
# 提取待分析的文件夹内所有突变vcf文件中的所有位点
# 以文件为单位，获取每个case中的所有位点信息 + 30 features取值
# 将所获取的多个特征文件保存于指定文件夹中
# 最后将文件夹内所有特征文件（所有case）纵向合并，得到最终的训练集、测试集、验证集

# 2021.1.4更新：从特征矩阵中删除表达量（TPM）相关特征信息，因并非所有数据都有输入的TPM信息
# 2021.4.8更新：添加record_filter列，保存突变位点相关filter信息，用以评估有多少filter过滤的record被rescue回来

import vcf
import pandas as pd
import os
from multiprocessing import Pool
import pysam
import traceback
import collections

import argparse

# description参数可以用于描述脚本的参数作用，默认为空
parser=argparse.ArgumentParser(description="A discriminate model construction pipeline for RNA-SSNV.")
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.')
# Generic parameter
parser.add_argument('--cancer_type', help='Cancer type for current study')
parser.add_argument("--RNA_calling_info", help="Tabular info for RNA somatic mutation calling.")
parser.add_argument('--project_folder', help='Project folder path.')
parser.add_argument('--exon_interval', help='GENCODE v22 exon interval bed file.')
parser.add_argument("--output_table_path", help="Path for final output table.")
parser.add_argument('--num_threads', '-nt', type=int, help='Number of threads allowed.')

args=parser.parse_args()

# 主要处理函数
# 输入信息
# tumor_tsv：用以获取当前的case对应的其他相关信息（如配对正常样本信息等）
# 对于每个case id对应的vcf文件，进行文件信息提取，并将获取的相应feature信息保存至指定文件夹内
def case_vcf_info_retrive(tumor_tsv, case_id, vcf_folder_path, table_folder_path, Exon_loc, funcotator_folder_path):
    """Retrieve all features based on given case_id.

    Retrieve specific case_id's vcf-related, tumor-DNA-bam-related, annotation-related and other-related (exon distance, etc) features.
    Store features above into provided table folder in order to facilitate multi-process.

    Args:
        tumor_tsv: A tsv file containing RNA somatic calling info.
        case_id: Specific case id.
        vcf_folder_path: A path to the vcf folder containing case's vcf file (format: vcf.gz).
        table_folder_path: A path to the table folder which will store case's feature-info table.
        DNA_tumor_tsv: A tsv file containing DNA somatic calling info.
        DNA_tumor_folder_path: A path to the bam folder containing case's DNA tumor bam.
        Exon_loc: A path to the bed file containing all GENCODE v22's exon info.
        funcotator_folder_path: A path to the maf folder containing case's annotation maf info.

    Raises:
        IOError: An error occurred accessing the tumor_tsv, DNA_tumor_tsv or Exon_loc. Table object.
    """
    try:
        # 读取对应case的vcf文件，并保存相应vcf单位点信息列表
        vcf_reader = vcf.Reader(filename=os.path.join(vcf_folder_path, case_id+".vcf.gz"))

        # 2021.6.11
        # 为使特征提取代码适应于包含"AS_UNIQ_ALT_READ_COUNT"信息的vcf文件，添加下述代码段
        # 直接修改INFO中类型值会报错，因此需要进行根据相应INFO格式进行转换
        pre_info = vcf_reader.infos['AS_UNIQ_ALT_READ_COUNT']
        after_info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc', 'source', 'version'])
        after_info = after_info(pre_info.id, pre_info.num,
                                "String", pre_info.desc,
                                pre_info.source, pre_info.version)
        vcf_reader.infos['AS_UNIQ_ALT_READ_COUNT'] = after_info
        # 2021.6.12
        # 为使特征提取代码适应于包含"AS_MQ"信息的vcf文件，添加下述代码段
        # 直接修改INFO中类型值会报错，因此需要进行根据相应INFO格式进行转换
        pre_info = vcf_reader.infos['AS_MQ']
        after_info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc', 'source', 'version'])
        after_info = after_info(pre_info.id, pre_info.num,
                                "String", pre_info.desc,
                                pre_info.source, pre_info.version)
        vcf_reader.infos['AS_MQ'] = after_info

        vcf_record_list = [record for record in vcf_reader]
        # 获取RNA tumor、DNA normal相应的所有aliquots_id（RNA somatic calling info中单个case内所包含的所有aliquots信息）
        RNA_tumor_aliquots_id = tumor_tsv.loc[(tumor_tsv['case_id']==case_id)&(tumor_tsv['sample_type']=="Primary Tumor"),'aliquots_id']
        DNA_normal_aliquots_id = tumor_tsv.loc[(tumor_tsv['case_id']==case_id)&(tumor_tsv['sample_type'].isin(["Blood Derived Normal", "Solid Tissue Normal"])),'aliquots_id']

        # 将经过人工筛选的vcf文件相应特征，引入DataFrame中便于分析——首先指定最后将要构建的dataframe列名——亦即位点基本信息+位点特征
        # 新增加column names
        # "Hugo_Symbol", "Variant_Classification", "Gencode_28_secondaryVariantClassification",
        # "HGNC_Ensembl_Gene_ID", "gc_content", "COSMIC_total_alterations_in_gene", "Expression_TPM"
        # 新增加column names（2021.6.11）
        # "ClippingRankSum",
        # "FS", "LikelihoodRankSum",
        # "SOR", "ROQ"
        # 部分column会出现缺失值状况，选择删除
        # "ReadPosRankSum", "BaseQRankSum"
        # "ClippingRankSum", "LikelihoodRankSum"
        col_names = ['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Allele1',
                     'Tumor_Allele2', 'Tumor_Sample_UUID',
                     'CONTQ', 'DP', 'ECNT', 'GERMQ',
                     'MBQ_ref', 'MBQ_alt',
                     'MFRL_ref', 'MFRL_alt',
                     'MMQ_ref', 'MMQ_alt',
                     'MPOS', 'NALOD', 'NLOD', 'POPAF', 'SEQQ', 'STRANDQ', 'TLOD',
                     'ref_AD_tumor_RNA', 'alt_AD_tumor_RNA',
                     'ref_AD_normal', 'alt_AD_normal',
                     'AF_tumor', 'AF_normal',
                     'DP_tumor', 'DP_normal',
                     'F1R2_tumor', 'F1R2_normal',
                     'F2R1_tumor', 'F2R1_normal',
                     'AD_other_RNA_tumor', 'AD_other_normal',
                     'transcript_ID', 'exon_distance',
                     'record_filter',
                     "FS",
                     "SOR", "ROQ"]
        # 添加GENCODE v22版本对应exon结构信息 + Funcotator注释信息

        # 读取GENCODE v22所对应的exon区间信息，建立以chromosome为key的exon dict对象存储exon区间信息
        exon_region_dict = Exon_region_create(Exon_loc)
        # 对case内的每个record进行分析，获取其基本信息+特征
        all_record_list = [record_select(single_record, case_id, RNA_tumor_aliquots_id, DNA_normal_aliquots_id, exon_region_dict) for single_record in vcf_record_list]
        # 汇总分离后的record
        all_record_df = pd.DataFrame(all_record_list, columns=col_names)

        # 信息补充
        # 添加Funcotator注释信息（读取相应funcotator注释信息后，再根据位置+碱基改变信息进行合并）
        # 2021.1.15 引入Strand、Transcript_Strand、Codon_Change、Protein_Change、DrugBank信息
        case_funcotator_df = pd.read_table(os.path.join(funcotator_folder_path, case_id+".maf"), comment="#")[["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2",
                                                                                                               "Hugo_Symbol", "Variant_Classification", "Gencode_28_secondaryVariantClassification",
                                                                                                               "HGNC_Ensembl_Gene_ID", "gc_content", "COSMIC_total_alterations_in_gene", "ref_context",
                                                                                                               "Strand", "Transcript_Strand", "Codon_Change", "Protein_Change", "DrugBank"]]
        # 2021.9.2 基于位置信息进行去重，保证multi-allelic处理过程中不出错（暂时注释掉）
        case_funcotator_df = case_funcotator_df.drop_duplicates(subset=["Chromosome", "Start_Position", "Reference_Allele"], keep='first', inplace=False)
        all_record_df = all_record_df.merge(case_funcotator_df, on=["Chromosome", "Start_Position", "Reference_Allele"], how="left")

        # 导出case table文件至对应文件夹中
        all_record_df.to_csv(os.path.join(table_folder_path, case_id+".table"), sep="\t", index=False)
        return all_record_df
    except Exception as ex:
        print(ex.args)
        print(traceback.format_exc())
        print("case_vcf_info_retrive function had failed!!!")
        print(f"其对应case信息为{case_id}")


# 其实这里也是偷懒，为了避免不同特征间合并起来的困难而在一个函数中进行
# 选择vcf文件中单个record内RNA肿瘤样本var最大AD、正常样本ref最大AD
# 同时选择DNA肿瘤样本var、ref最大AD和exon-distance来作为feature，共同组成dataframe
def record_select(record, case_id, RNA_tumor_aliquots_id, DNA_normal_aliquots_id, exon_region_dict):
    """Retrieve all features based on given case_id.

    Retrieve specific vcf record's vcf-related, tumor-DNA-bam-related and other-related (exon distance, etc) features.
    Store features above into a list with its order consistent with col_names in the case_vcf_info_retrive function.

    Args:
        record: A record within specific case_id's vcf info.
        case_id: Specific case id.
        RNA_tumor_aliquots_id: All aliquots_id of specific case_id's RNA tumor, help to identify corresponding vcf info.
        DNA_normal_aliquots_id: All aliquots_id of specific case_id's DNA normal, help to identify corresponding vcf info.
        samfile_list: A list containing all DNA tumors' pysam-objects which help to retrieve corresponding AD info.
        exon_region_dict: A dict whose keys were chromosome and values were regions' interval information.

    Raises:
        IOError: An error occurred accessing the tumor_tsv, DNA_tumor_tsv or Exon_loc. Table object.
    """
    try:
        # vcf-contained multiple cases' info retrieve
        # 初始化RNA tumor / DNA normal aliquots_id
        high_tumor_AD_aliquots_id = RNA_tumor_aliquots_id.values[0]
        high_normal_AD_aliquots_id = DNA_normal_aliquots_id.values[0]
        # 找出最大的RNA tumor aliquots_id——以便于获取相应tumor AD信息(最大的定义为最大的alt AD值)
        for single_tumor_aliquots_id in RNA_tumor_aliquots_id:
            if record.genotype(single_tumor_aliquots_id).data.AD[1] > record.genotype(high_tumor_AD_aliquots_id).data.AD[1]:
                high_tumor_AD_aliquots_id = single_tumor_aliquots_id
        # 找出最大的DNA normal aliquots_id——便于获取相应normal AD信息
        for single_normal_aliquots_id in DNA_normal_aliquots_id:
            if record.genotype(single_normal_aliquots_id).data.AD[0] > record.genotype(high_normal_AD_aliquots_id).data.AD[0]:
                high_normal_AD_aliquots_id = single_normal_aliquots_id
        # 获得其他aliquots的AD信息以全面获取features
        AD_other_tumor = ";".join(["/".join([str(record.genotype(single_tumor_aliquots_id).data.AD[0]), str(record.genotype(single_tumor_aliquots_id).data.AD[1])])
                               for single_tumor_aliquots_id in RNA_tumor_aliquots_id.drop(high_tumor_AD_aliquots_id)])
        AD_other_normal = ";".join(["/".join([str(record.genotype(single_normal_aliquots_id).data.AD[0]), str(record.genotype(single_normal_aliquots_id).data.AD[1])])
                               for single_normal_aliquots_id in DNA_normal_aliquots_id.drop(high_normal_AD_aliquots_id)])

        # other related info: exon distance
        # 获得GENCODE v22所对应的exon区间位置信息
        transcript_ID, exon_distance = Exon_region_record_retrive(record, exon_region_dict)

        # other maybe useful info: record_filter
        record_filter = "PASS" if record.FILTER==[] else ",".join(record.FILTER)

        # 根据获得的最大AD对应aliquots_id信息，来完成相应dataframe组成元素的构建
        # related column names
        # col_names = ['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Allele1',
        #              'Tumor_Allele2', 'Tumor_Sample_UUID',
        #              'CONTQ', 'DP', 'ECNT', 'GERMQ',
        #              'MBQ_ref', 'MBQ_alt',
        #              'MFRL_ref', 'MFRL_alt',
        #              'MMQ_ref', 'MMQ_alt',
        #              'MPOS', 'NALOD', 'NLOD', 'POPAF', 'SEQQ', 'STRANDQ', 'TLOD',
        #              'ref_AD_tumor_RNA', 'alt_AD_tumor_RNA',
        #              'ref_AD_normal', 'alt_AD_normal',
        #              'AF_tumor', 'AF_normal',
        #              'DP_tumor', 'DP_normal',
        #              'F1R2_tumor', 'F1R2_normal',
        #              'F2R1_tumor', 'F2R1_normal',
        #              'AD_other_RNA_tumor', 'AD_other_normal',
        #              'transcript_ID', 'exon_distance',
        #              'record_filter',
        #              "FS",
        #              "SOR", "ROQ"]
        tmp_vcf_info = [record.CHROM, record.POS, record.REF, ",".join([str(alt) for alt in record.ALT]),
                        record.REF, case_id,
                        record.INFO['CONTQ'], record.INFO['DP'], record.INFO['ECNT'], record.INFO['GERMQ'],
                        record.INFO['MBQ'][0], record.INFO['MBQ'][1],
                        record.INFO['MFRL'][0], record.INFO['MFRL'][1],
                        record.INFO['MMQ'][0], record.INFO['MMQ'][1],
                        record.INFO['MPOS'][0], record.INFO['NALOD'][0], record.INFO['NLOD'][0],
                        record.INFO['POPAF'][0], record.INFO['SEQQ'], record.INFO['STRANDQ'], record.INFO['TLOD'][0],
                        record.genotype(high_tumor_AD_aliquots_id).data.AD[0], record.genotype(high_tumor_AD_aliquots_id).data.AD[1],
                        record.genotype(high_normal_AD_aliquots_id).data.AD[0], record.genotype(high_normal_AD_aliquots_id).data.AD[1],
                        record.genotype(high_tumor_AD_aliquots_id).data.AF, record.genotype(high_normal_AD_aliquots_id).data.AF, # 存在形如[0.015, 0.005935]的AF值，而不报错的情况，通过强制类型转换来避免出错
                        record.genotype(high_tumor_AD_aliquots_id).data.DP, record.genotype(high_normal_AD_aliquots_id).data.DP,
                        record.genotype(high_tumor_AD_aliquots_id).data.F1R2[1], record.genotype(high_normal_AD_aliquots_id).data.F1R2[0],
                        record.genotype(high_tumor_AD_aliquots_id).data.F2R1[1], record.genotype(high_normal_AD_aliquots_id).data.F2R1[0],
                        AD_other_tumor, AD_other_normal,
                        transcript_ID, exon_distance,
                        record_filter,
                        record.INFO['FS'],
                        record.INFO['SOR'], record.INFO['ROQ']
                        ]
        return tmp_vcf_info
    except Exception as ex:
        # print("record_select错误！！！对应Record信息如下所示")
        # print(record)
        # print(record.INFO)
        # print(f"其对应case信息为{case_id}")

        tmp_vcf_info = [record.CHROM, record.POS, record.REF, ",".join([str(alt) for alt in record.ALT]),
                        record.REF, case_id,
                        record.INFO['CONTQ'], 0, record.INFO['ECNT'], record.INFO['GERMQ'],
                        record.INFO['MBQ'][0], record.INFO['MBQ'][1],
                        record.INFO['MFRL'][0], record.INFO['MFRL'][1],
                        record.INFO['MMQ'][0], record.INFO['MMQ'][1],
                        record.INFO['MPOS'][0], record.INFO['NALOD'][0], record.INFO['NLOD'][0],
                        record.INFO['POPAF'][0], record.INFO['SEQQ'], record.INFO['STRANDQ'], record.INFO['TLOD'][0],
                        0, 0,
                        0, 0,
                        0, 0,
                        0, 0,
                        0, 0,
                        0, 0,
                        0, 0,
                        0, 0,
                        "PASS" if record.FILTER==[] else ",".join(record.FILTER),
                        record.INFO['FS'],
                        record.INFO['SOR'], record.INFO['ROQ']
                        ]
        return tmp_vcf_info


# 根据vcf文件中单个record的突变信息，获取所有特定case的DNA肿瘤样本中var最大AD及其他相应信息，其返回值纳入其他信息的dataframe中
def DNA_tumor_bam_record_retrive(record, samfile_list):
    try:
        # 获取基本信息
        ref = record.REF
        alt = record.ALT[0]
        chrom = record.CHROM
        position = record.POS
        # 获取每个bam文件中对应的碱基coverage信息
        other_coverage_list = []
        # 应用首个bam文件中coverage信息完成初始化(设置base的最低质量值为30)
        initial_count_coverage_info = count_coverage_decompose(samfile_list[0].count_coverage(contig=chrom, start=position-1, stop=position, quality_threshold=30))
        # 什么乱七八糟的数据结构啊= =我佛了，还需要强行转为str后才能作为key去访问字典中元素
        optimal_ref_coverage = initial_count_coverage_info[str(ref)]
        optimal_alt_coverage = initial_count_coverage_info[str(alt)]
        # 遍历所有bam文件，并进行条件判断
        for samfile in samfile_list[1:]:
            # pysam apply 0-based coordinate system(设置base的最低质量值为30)
            count_coverage_info = count_coverage_decompose(samfile.count_coverage(contig=chrom, start=position-1, stop=position, quality_threshold=30))
            # 当前覆盖度信息
            current_ref_coverage = count_coverage_info[str(ref)]
            current_alt_coverage = count_coverage_info[str(alt)]
            # 当前覆盖度信息进行对比
            # 比较alt的coverage信息
            if current_alt_coverage > optimal_alt_coverage:
                # 保存当前值，并将更好的值重新赋值
                other_coverage_list.append(str(optimal_ref_coverage)+"/"+str(optimal_alt_coverage))
                optimal_ref_coverage, optimal_alt_coverage = current_ref_coverage, current_alt_coverage
            elif current_alt_coverage == optimal_alt_coverage:
                # 进一步比较ref的coverage信息
                if current_ref_coverage > optimal_ref_coverage:
                    # 保存当前值，并将更好的值重新赋值
                    other_coverage_list.append(str(optimal_ref_coverage) + "/" + str(optimal_alt_coverage))
                    optimal_ref_coverage, optimal_alt_coverage = current_ref_coverage, current_alt_coverage
                else:
                    # 保存当前值，不进行更新
                    other_coverage_list.append(str(current_ref_coverage) + "/" + str(current_alt_coverage))
            else:
                # 保存当前值，不进行更新
                other_coverage_list.append(str(current_ref_coverage) + "/" + str(current_alt_coverage))
        return optimal_ref_coverage, optimal_alt_coverage, ";".join(other_coverage_list)
    except Exception as ex:
        print(ex)
        print("DNA_tumor_bam_record_retrive错误！！！")


# 将pysam中count_coverage函数的返回值转换为字典形式
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


# 根据record中的chrom pos来确定exon位置
def Exon_region_create(Exon_loc):
    exon_region_df = pd.read_table(Exon_loc, sep="\t")
    exon_region_df.columns=["chr", "start", "end", "exon_info", "count", "strand"]
    # 将bed格式转为position格式
    exon_region_df["start"] = exon_region_df["start"] + 1
    # 新建dict保存各染色体对应dataframe信息
    exon_region_dict = {}
    # 遍历每个染色体上的exon相关信息
    for single_chr in set(exon_region_df["chr"]):
        exon_region_df_single_chr = pd.DataFrame(exon_region_df.loc[exon_region_df["chr"] == single_chr, ])
        # 对具体坐标进行排序——根据start从小到大进行排序
        exon_region_df_single_chr.sort_values("start", inplace=True)
        # 添加排序后dataframe进入dict内
        exon_region_dict[single_chr] = exon_region_df_single_chr
    # 返回保存各染色体对应dataframe信息的dict
    return exon_region_dict


def Exon_region_record_retrive(record, exon_region_dict):
    try:
        # 获取基本信息
        chrom = str(record.CHROM)
        position = int(record.POS)
        # 进入染色体对应dataframe
        current_chrom_df = exon_region_dict[chrom]
        # 找出距离最近的exon区间——判断位于exon内或exon外
        # 分别找出离start、end位置的距离
        start_pos_distance = position - current_chrom_df["start"]
        end_pos_distance = current_chrom_df["end"] - position
        # 若其离start position更近
        if start_pos_distance.abs().min() <= end_pos_distance.abs().min():
            # 找出其所对应的index
            min_dis_exon_index = start_pos_distance.abs().idxmin()
            # 确定返回值
            transcript_ID = current_chrom_df.loc[min_dis_exon_index, "exon_info"].split(".")[0]
            exon_distance = start_pos_distance.abs().min()
            return (transcript_ID, exon_distance)
        # 若其离end position更近
        elif end_pos_distance.abs().min() < start_pos_distance.abs().min():
            # 找出其所对应的index
            min_dis_exon_index = end_pos_distance.abs().idxmin()
            # 确定返回值
            transcript_ID = current_chrom_df.loc[min_dis_exon_index, "exon_info"].split(".")[0]
            exon_distance = end_pos_distance.abs().min()
            return (transcript_ID, exon_distance)
    except Exception as ex:
        print(ex)
        print("Exon_region_record_retrive错误！！！")


if __name__ == '__main__':
    CANCER_TYPE = args.cancer_type

    # Input（仅需要在此修改相关路径信息即可）
    # 包含有RNA tumor相关注释信息的tsv文件，通常为RNA体细胞突变检测对应信息列表文件
    tumor_tsv = pd.read_table(args.RNA_calling_info)
    tumor_tsv['case_id'] = tumor_tsv['case_id'].astype(str)
    # 存储该项目所有case的突变vcf文件的文件夹路径
    vcf_folder_path = os.path.join(args.project_folder, args.cancer_type, "RNA/RNA_somatic_mutation/SelectVariants_new/SNP_WES_Interval_exon")
    # TODO: 更新不同基因对应的exon structure信息
    # 根据UCSC数据库中GENCODE v22中的exon信息所构造的结构信息文件为/home/lqh/resources/database/gencode/GRCh38_GENCODE_v22_exon_rm_alt.bed
    Exon_loc = args.exon_interval
    # 存储Funcotator注释后信息文件（完成对所有突变位点的注释）的文件夹路径
    funcotator_folder_path = os.path.join(args.project_folder, args.cancer_type, "RNA/RNA_somatic_mutation/Funcotator_new/SNP")

    # Output（修改相应路径信息即可）
    table_folder_path = os.path.join(args.project_folder, args.cancer_type, "RNA/RNA_somatic_mutation/VcfAssembly_new/SNP_WES_Interval_exon")
    if not os.path.exists(table_folder_path):
        os.makedirs(table_folder_path)
    # 修改tumor_tsv行名
    tumor_tsv.index = tumor_tsv["aliquots_id"]
    # 对于指定文件夹内的所有记录在案vcf文件进行整合
    all_vcf_info = pd.DataFrame()
    # 开始多进程处理
    print('Parent process %s.' % os.getpid())
    p = Pool(args.num_threads)
    # 遍历每一个case(直接用set()也可以= =)
    for single_case_id in tumor_tsv["case_id"].value_counts().index:
        p.apply_async(case_vcf_info_retrive, args=(tumor_tsv, single_case_id, vcf_folder_path, table_folder_path,
                                                   Exon_loc,
                                                   funcotator_folder_path))
    print('Waiting for all subprocesses done...')
    p.close()
    p.join()
    print('All subprocesses done.')
    # 汇总多进程处理后结果——本质上是对所提供vcf文件对应突变位点信息的集合汇总（因为本代码是基于每个突变vcf文件内的所有位点进行分析）
    single_vcf_info_list = [pd.read_table(os.path.join(os.path.join(table_folder_path, single_case_id+".table"))) for single_case_id in tumor_tsv["case_id"].value_counts().index]
    all_vcf_info = pd.concat(single_vcf_info_list, ignore_index=True)
    all_vcf_info.to_csv(f"/home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/{CANCER_TYPE}/RNA/RNA_somatic_mutation/VcfAssembly_new/SNP_WES_Interval_exon.txt", sep="\t", index=False)