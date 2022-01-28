# test successful
# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_utility.py \
# --cancer_type LUSC \
# --Intogen_database_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/tables/gene/Compendium_Cancer_Genes.tsv \
# --result_folder /home/lqh/Codes/Python/RNA-SSNV/output \
# --output_info /home/lqh/Codes/Python/RNA-SSNV/results/LUSC.RNA_DNA_comparison_matrix.driver_gene.csv \
# --utilize_type driver_gene
#
python /home/lqh/Codes/Python/RNA-SSNV/lib/result_utility.py \
--cancer_type GBM \
--Intogen_database_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/tables/gene/Compendium_Cancer_Genes.tsv \
--result_folder /home/lqh/Codes/Python/RNA-SSNV/output \
--output_info /home/lqh/Codes/Python/RNA-SSNV/results/GBM.RNA_DNA_comparison_matrix.csv \
--utilize_type all

import pandas as pd
import os

import argparse

parser=argparse.ArgumentParser(description="Extract required info from given DNA and RNA datasets.")
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.')
# Generic parameter
parser.add_argument('--cancer_type', help='Type of cancer.')
parser.add_argument('--Intogen_database_info', help='Intogen database info location.')
parser.add_argument('--result_folder', help='Folder containing results.')
parser.add_argument('--output_info', help='Extracted info from DNA and RNA datasets.')
# Special parameter
parser.add_argument('--utilize_type', help='Type of result utilization.')

args = parser.parse_args()

def value_counts_statistic(df_list, col_name, columns_list, calculate_type=None, mutation_type=None, cancer_driver_gene_list=None):
    result_list = []

    if calculate_type=="median":
        for i in range(len(df_list)):
            median_value_dict = {}
            for case in df_list[i]['Tumor_Sample_UUID'].value_counts().index:
                median_value_dict[case] = df_list[i][df_list[i]['Tumor_Sample_UUID'] == case][col_name].median()
            median_value_df = pd.DataFrame(pd.Series(median_value_dict))
            median_value_df.columns = [columns_list[i]]
            result_list.append(median_value_df)
    else:
        for i in range(len(df_list)):
            current_df = df_list[i]
            if mutation_type=="non_synonymous":
                current_df = current_df[current_df.Variant_Classification != "Silent"]
            elif mutation_type=="non_synonymous_non_driver":
                current_df = current_df[current_df.Variant_Classification != "Silent"]
                current_df = current_df[~current_df.Hugo_Symbol.isin(cancer_driver_gene_list)]
            column_counts = pd.DataFrame(current_df[col_name].value_counts())
            column_counts.columns = [columns_list[i]]
            result_list.append(column_counts)

    return pd.concat(result_list, axis=1)

def precision_recall_statistic(assessment_df):

    result_list = []
    case_list = list(assessment_df['Tumor_Sample_UUID'].value_counts().index)

    for case in case_list:
        case_assessment_df = assessment_df[assessment_df.Tumor_Sample_UUID == case]

        raw_TP_count = len(case_assessment_df[case_assessment_df.DNA_label == 1])
        positive_TP_count = len(case_assessment_df[(case_assessment_df.DNA_label == 1)&(case_assessment_df.pred_label == 1)])
        positive_count = len(case_assessment_df[case_assessment_df.pred_label == 1])

        if positive_count == 0:
            print(f"{case}不存在positive位点，precision为0")
            precision = 0
        else:
            precision = round((positive_TP_count / positive_count)*100, 2)

        if raw_TP_count == 0:
            print(f"{case}不存在TP位点，precision&recall均为0")
            recall = 0
        else:
            recall = round((positive_TP_count / raw_TP_count)*100, 2)

        result_list.append([positive_TP_count, precision, recall])

    return pd.DataFrame(result_list, columns=["true_RNA_mutation_counts", "precision", 'recall'], index=case_list)

# 获取驱动基因相关信息
def driver_gene_info_extract(driver_gene_loc, cancer_type):
    print("正在读取全部驱动基因信息")
    driver_genes = pd.read_table(driver_gene_loc)
    print(f"开始获取{cancer_type}所对应的驱动基因集合")
    cancer_spec_driver_genes = set(driver_genes[driver_genes.CANCER_TYPE == cancer_type]['SYMBOL'])
    print(f"{cancer_type}所对应的驱动基因集合获取完成，共计{len(cancer_spec_driver_genes)}个驱动基因，名称为cancer_spec_driver_genes")
    return(cancer_spec_driver_genes)

if __name__ == '__main__':
    cancer_type = args.cancer_type
    cancer_spec_driver_genes = driver_gene_info_extract(args.Intogen_database_info, args.cancer_type)

    # 读取相关信息
    if args.utilize_type == "all":
        RNA_total = pd.read_table(os.path.join(args.result_folder, f"{cancer_type}.RNA_total.table"))
        RNA_total_all = pd.read_table(os.path.join(args.result_folder, f"{cancer_type}.RNA_total.all.table"))
        DNA_total = pd.read_table(os.path.join(args.result_folder, f"{cancer_type}.DNA_total.table"))
    elif args.utilize_type == "driver_gene":
        RNA_total = pd.read_table(os.path.join(args.result_folder, f"{cancer_type}.RNA_total.driver_gene.table"))
        RNA_total_all = pd.read_table(os.path.join(args.result_folder, f"{cancer_type}.RNA_total.driver_gene.all.table"))
        DNA_total = pd.read_table(os.path.join(args.result_folder, f"{cancer_type}.DNA_total.driver_gene.table"))

    # 报告整体信息
    raw_TP_count = len(RNA_total_all[RNA_total_all.DNA_label == 1])
    positive_TP_count = len(RNA_total_all[(RNA_total_all.DNA_label == 1) & (RNA_total_all.pred_label == 1)])
    positive_count = len(RNA_total_all[RNA_total_all.pred_label == 1])
    print(f"整体precision为{positive_TP_count/positive_count},"
          f"整体recall为{positive_TP_count/raw_TP_count}")

    # 添加突变counts信息（整体、非同义突变两类）
    mutation_counts_df = value_counts_statistic([DNA_total, RNA_total, RNA_total_all], 'Tumor_Sample_UUID',
                                                ['DNA_mutation_counts', 'predicted_RNA_mutation_counts',
                                                 'all_RNA_mutation_counts'])
    non_syn_mutation_counts_df = value_counts_statistic([DNA_total, RNA_total, RNA_total_all], 'Tumor_Sample_UUID',
                                                ['non_syn_DNA_mutation_counts', 'non_syn_predicted_RNA_mutation_counts',
                                                 'non_syn_all_RNA_mutation_counts'], mutation_type="non_synonymous")
    non_syn_non_driver_mutation_counts_df = value_counts_statistic([DNA_total, RNA_total, RNA_total_all], 'Tumor_Sample_UUID',
                                                ['non_syn_non_driver_DNA_mutation_counts', 'non_syn_non_driver_predicted_RNA_mutation_counts',
                                                 'non_syn_non_driver_all_RNA_mutation_counts'], mutation_type="non_synonymous_non_driver", cancer_driver_gene_list=cancer_spec_driver_genes)

    # 添加中位数coverage
    coverage_df = value_counts_statistic([DNA_total, RNA_total, RNA_total_all], 'DP_tumor',
                                         ['DNA_mutation_median_coverage', 'predicted_RNA_mutation_median_coverage',
                                          'all_RNA_mutation_median_coverage'],
                                         'median')
    # 添加P-R信息
    PR_df = precision_recall_statistic(RNA_total_all)

    # 汇总RNA、DNA层面突变信息（各类可用信息进行比较+离群点探究）
    # 需注意存在为0或缺失的情况（因为所取DNA位点并非均位于编码区内，存在为0的情况）
    RNA_DNA_comparison_matrix = pd.concat([non_syn_non_driver_mutation_counts_df, non_syn_mutation_counts_df, mutation_counts_df, coverage_df, PR_df], axis=1)
    RNA_DNA_comparison_matrix = RNA_DNA_comparison_matrix.fillna(0)

    RNA_DNA_comparison_matrix.to_csv(args.output_info, index=True)