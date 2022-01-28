# test success
# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_miscellaneous.py \
# --analyze_type preprocess_assessment \
# --data_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/GBM/RNA/RNA_somatic_mutation/VcfAssembly_new/SNP_WES_Interval_exon.txt \
# --GDC_info /home/lqh/Codes/Data/TCGA_maf_files/TCGA-GBM \
# --REDIportal_info /home/lqh/resources/database/RNA_edit/REDIportal/REDIportal_main_table.hg38.bed \
# --DARNED_info /home/lqh/resources/database/RNA_edit/DARNED_hg19_to_bed_to_hg38_rm_alt.bed \
# --metric_info /home/lqh/Codes/Python/RNA-SSNV/results/GBM.multi_filtering.metric

import pandas as pd

import argparse

parser=argparse.ArgumentParser(description="Miscellaneous tasks for results generation.")
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.')
# General parameters
parser.add_argument('--analyze_type', help='Type of analysis.')
parser.add_argument('--data_info', help='Location of RNA somatic mutations info.')
# Multi-filtering strategy related parameters
parser.add_argument('--GDC_info', help='GDC maf file location.')
parser.add_argument('--REDIportal_info', help='REDIportal database file location.')
parser.add_argument('--DARNED_info', help='DARNED database file location.')
# output info
parser.add_argument('--metric_info', help='Location of performance metric output.')

#
args = parser.parse_args()

# 整理GDC maf数据集相应信息（取SNP+去重）
def GDC_trim(GDC_info):
    ## 整理GDC相应信息
    print("GDC项目中的突变位点情况为：\n", GDC_info["Variant_Type"].value_counts())
    ### 仅选择SNP突变信息进行分析
    GDC_SNP_info = GDC_info[GDC_info['Variant_Type'] == "SNP"]
    ### 重新编制index信息
    GDC_SNP_info.reset_index(drop=True, inplace=True)
    ### 将“Tumor_Sample_UUID”列中信息进行切片&提取，使之与RNA体细胞突变中的case_id信息保持一致
    del GDC_SNP_info['Tumor_Sample_UUID']
    GDC_SNP_info["Tumor_Sample_UUID"] = pd.Series(["-".join(GDC_sample_info.split("-")[0:3])
                                                        for GDC_sample_info in
                                                        GDC_SNP_info["Tumor_Sample_Barcode"]])
    ### 仅选取染色体、碱基和case_id信息来进行合并
    GDC_SNP_info = GDC_SNP_info.loc[:,
                        ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_UUID"]]
    ### 重新命名，便于进行合并(Tumor_Allele1对应突变碱基、Tumor_Allele2对应参考碱基)
    GDC_SNP_info.columns = ["Chromosome", "Start_Position", "Tumor_Allele2", "Tumor_Allele1", "Tumor_Sample_UUID"]
    ### GDC突变去重（由于多个肿瘤样本数据随机组合的缘故，体细胞突变会存在重复），在LUAD项目中，去重前后比例为3:1
    GDC_SNP_info = GDC_SNP_info.drop_duplicates(keep="first")
    print(f"整理后，GDC项目中的突变数目共计{len(GDC_SNP_info)}")

    return GDC_SNP_info

# 整理RNA编辑位点相关信息
def RNA_edit_combine(REDIportal_info, DARNED_info):
    ## 对数据库信息进行预处理，提取所需信息
    REDIportal_info = REDIportal_info[[0, 2]]
    DARNED_info = DARNED_info[[0, 2]]
    print(f"REDIportal和DARNED数据库中的RNA编辑位点数目分别为{len(REDIportal_info)}和{len(DARNED_info)}")
    ## 对REDIpotal和DARNED数据库信息进行合并
    RNA_EDIT_INFO = pd.concat([REDIportal_info, DARNED_info], ignore_index=True)
    RNA_EDIT_INFO.columns = ["Chromosome", "Start_Position"]
    ## 合并所得数据库信息需要去重
    RNA_EDIT_INFO = RNA_EDIT_INFO.drop_duplicates(keep="first")
    print(f"合并数据库信息并去重后，所得的RNA编辑位点总数为{len(RNA_EDIT_INFO)}")

    return RNA_EDIT_INFO

# 获取FilterMutectCalls表现矩阵
def RNA_GDC_analysis(RNA_info, GDC_info):
    print(f"RNA突变数据集和GDC突变数据集对应突变位点数目为{len(RNA_info)}、{len(GDC_info)}")
    # split RNA_info into two parts: RNA_info_TP, RNA_info_TN
    RNA_info_TP = pd.merge(RNA_info, GDC_info, on=list(GDC_info.columns))
    print(f"RNA_info_TP类获取完成，数目为{len(RNA_info_TP)}")
    RNA_info = RNA_info.append(RNA_info_TP)
    RNA_info_TN = RNA_info.drop_duplicates(keep=False)
    print(f"RNA_info_TN类获取完成，数目为{len(RNA_info_TN)}")
    # retrieve confusion matrix for RNA_info
    print("开始获取基于FilterMutectCalls的混淆矩阵")
    FilterMutectCalls_TP = RNA_info_TP[RNA_info_TP['record_filter']=="PASS"]
    FilterMutectCalls_FN = RNA_info_TP[RNA_info_TP['record_filter'] != "PASS"]
    FilterMutectCalls_TN = RNA_info_TN[RNA_info_TN['record_filter'] != "PASS"]
    FilterMutectCalls_FP = RNA_info_TN[RNA_info_TN['record_filter'] == "PASS"]
    confusion_matrix = pd.DataFrame([[len(FilterMutectCalls_TP), len(FilterMutectCalls_FN)], [len(FilterMutectCalls_FP), len(FilterMutectCalls_TN)]],
                                    index=['Real_Positive', 'Real_Negative'],
                                    columns=['Predict_Positive', 'Predict_Negative'])
    print(confusion_matrix)
    print(f"Precision为{len(FilterMutectCalls_TP)/(len(FilterMutectCalls_TP) + len(FilterMutectCalls_FP))}")
    print(f"Recall为{len(FilterMutectCalls_TP) / (len(FilterMutectCalls_TP) + len(FilterMutectCalls_FN))}")

# 检验多重过滤过程中GDC突变变化情况
def multi_filtering_stats(all_info, GDC_SNP_info, RNA_EDIT_INFO):
    # 构建list来存储区间多层过滤后所得结果
    # 存储结构为(filter_stage, total_RNA_num, GDC_num, non_GDC_num)
    result_list = []
    # 开始对所有RNA突变进行筛选，筛选掉基本无法处理&没必要处理的multiallelic位点——self.all_info
    print(f"所有RNA突变在进行multi-allelic筛选前，对应的突变数目为{len(all_info)}")

    all_info_GDC = pd.merge(all_info, GDC_SNP_info, on=list(['Chromosome', 'Start_Position', 'Tumor_Sample_UUID']))
    print(f"其中，位于GDC数据库中的突变数为{len(all_info_GDC)}，不位于GDC数据库中突变数为{len(all_info) - len(all_info_GDC)}")
    result_list.append(['raw', len(all_info), len(all_info_GDC), len(all_info) - len(all_info_GDC)])
    print("=" * 100)

    all_info = all_info[
        all_info['record_filter'].apply(lambda x: False if x.__contains__("multiallelic") else True)]
    print(f"所有RNA突变在经过multi-allelic筛选后，剩余的突变数目为{len(all_info)}")

    all_info_GDC = pd.merge(all_info, GDC_SNP_info, on=list(GDC_SNP_info.columns))
    all_info_non_GDC = all_info.append(all_info_GDC)
    all_info_non_GDC = all_info_non_GDC.drop_duplicates(keep=False)

    print(f"其中，位于GDC数据库中的突变数为{len(all_info_GDC)}，不位于GDC数据库中突变数为{len(all_info_non_GDC)}")
    result_list.append(['multi-allelic', len(all_info), len(all_info_GDC), len(all_info_non_GDC)])
    print("=" * 100)

    # 进一步对RNA编辑位点数据库进行分析，判断其是否位于数据库内
    ## 找出位于RNA编辑数据库中的相应RNA突变位点——保存于self.all_info_RNA_edit中
    all_info_RNA_edit = pd.merge(all_info, RNA_EDIT_INFO)
    all_info = all_info.append(all_info_RNA_edit)
    all_info = all_info.drop_duplicates(keep=False)

    print(f"位于RNA编辑位点上的突变数为{len(all_info_RNA_edit)}（对象名称为all_info_RNA_edit），"
          f"不属于RNA编辑位点的突变数为{len(all_info)}")

    all_info_GDC = pd.merge(all_info, GDC_SNP_info, on=list(GDC_SNP_info.columns))
    all_info_non_GDC = all_info.append(all_info_GDC)
    all_info_non_GDC = all_info_non_GDC.drop_duplicates(keep=False)

    print(f"其中，位于GDC数据库中的突变数为{len(all_info_GDC)}，不位于GDC数据库中突变数为{len(all_info_non_GDC)}")
    result_list.append(['RNA-edit', len(all_info), len(all_info_GDC), len(all_info_non_GDC)])
    print("=" * 100)

    # 进一步排除位于免疫球蛋白基因内突变
    all_info_immunoglobulin = all_info[
        all_info.apply(lambda x: True if x['Hugo_Symbol'].startswith(("IGH", 'IGK', 'IGL')) else False, axis=1)]
    all_info = all_info[
        all_info.apply(lambda x: False if x['Hugo_Symbol'].startswith(("IGH", 'IGK', 'IGL')) else True, axis=1)]

    print(f"位于免疫球蛋白基因内的突变数为{len(all_info_immunoglobulin)}（对象名称为all_info_immunoglobulin），"
          f"不位于免疫球蛋白基因内的突变数为{len(all_info)}")

    all_info_GDC = pd.merge(all_info, GDC_SNP_info, on=list(GDC_SNP_info.columns))
    all_info_non_GDC = all_info.append(all_info_GDC)
    all_info_non_GDC = all_info_non_GDC.drop_duplicates(keep=False)

    print(f"其中，位于GDC数据库中的突变数为{len(all_info_GDC)}，不位于GDC数据库中突变数为{len(all_info_non_GDC)}")
    result_list.append(['immunoglobulin', len(all_info), len(all_info_GDC), len(all_info_non_GDC)])
    print("=" * 100)

    # 进一步排除HLA基因内突变
    all_info_HLA = all_info[
        all_info.apply(lambda x: True if x['Hugo_Symbol'].startswith("HLA") else False, axis=1)]
    all_info = all_info[
        all_info.apply(lambda x: False if x['Hugo_Symbol'].startswith("HLA") else True, axis=1)]

    print(f"位于HLA基因内的突变数为{len(all_info_HLA)}（对象名称为all_info_HLA），"
          f"不位于HLA基因内的突变数为{len(all_info)}")

    all_info_GDC = pd.merge(all_info, GDC_SNP_info, on=list(GDC_SNP_info.columns))
    all_info_non_GDC = all_info.append(all_info_GDC)
    all_info_non_GDC = all_info_non_GDC.drop_duplicates(keep=False)

    print(f"其中，位于GDC数据库中的突变数为{len(all_info_GDC)}，不位于GDC数据库中突变数为{len(all_info_non_GDC)}")
    result_list.append(['HLA', len(all_info), len(all_info_GDC), len(all_info_non_GDC)])

    print("=" * 100)

    return(pd.DataFrame(result_list, columns=['filter_stage', 'total_RNA_num', 'GDC_num', 'non_GDC_num']), all_info)

if __name__ == '__main__':

    # 读取RNA体细胞突变信息
    RNA_info = pd.read_table(args.data_info)

    # 读取gdc数据库突变信息
    GDC_info = pd.read_table(args.GDC_info)
    # 整理gdc数据库信息
    GDC_SNP_info = GDC_trim(GDC_info)

    # 读取RNA编辑数据库信息
    REDIportal_info = pd.read_table(args.REDIportal_info, header=None)
    DARNED_info = pd.read_table(args.DARNED_info, header=None)
    RNA_EDIT_INFO = RNA_edit_combine(REDIportal_info, DARNED_info)

    # 根据case情况求得最终纳入分析的所有RNA、DNA突变
    GDC_case_info_set = set(GDC_SNP_info['Tumor_Sample_UUID'].value_counts().index)
    RNA_case_info_set = set(RNA_info['Tumor_Sample_UUID'].value_counts().index)
    # 找到RNA中缺少GDC数据支持的case
    RNA_only_case_info_set = RNA_case_info_set - GDC_case_info_set
    # 找到两者交叉的case信息
    RNA_DNA_overlap_case_info_set = set.intersection(RNA_case_info_set, GDC_case_info_set)
    print(f"RNA突变集包含的case数为{len(RNA_case_info_set)}，GDC DNA突变集包含的case数为{len(GDC_case_info_set)}"
          f"其中，RNA中缺少GDC数据支持的case信息为{RNA_only_case_info_set}，将其排除后得到最终RNA、DNA突变集的case数为{len(RNA_DNA_overlap_case_info_set)}")
    # 获取最终用以评估的GDC DNA突变集（具有证据支持）
    # RNA突变集并不求子集，保留仅在RNA中出现的case中假阴性位点用以训练
    RNA_info = RNA_info[RNA_info['Tumor_Sample_UUID'].isin(RNA_DNA_overlap_case_info_set)]
    GDC_SNP_info = GDC_SNP_info[GDC_SNP_info['Tumor_Sample_UUID'].isin(RNA_DNA_overlap_case_info_set)]

    # 根据分析需求的区别，确定分析路线
    if args.analyze_type == "preprocess_assessment":
        print("="*100)
        print("开始进行多重过滤并评估其过滤效果")
        multi_filtering_metric, RNA_info = multi_filtering_stats(RNA_info, GDC_SNP_info, RNA_EDIT_INFO)
        multi_filtering_metric.to_csv(args.metric_info, sep="\t", index=False)
        print("最终，多重过滤结果矩阵为")
        print(multi_filtering_metric)

        print("="*100)
        print("应用多重过滤后结果，评估FilterMutectCalls效果")
        RNA_GDC_analysis(RNA_info, GDC_SNP_info)