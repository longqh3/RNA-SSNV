# test success
# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_driver_extractor.py \
# --result_info /home/lqh/Codes/Python/RNA-SSNV/output/GBM.RNA_only.table \
# --cancer_type GBM \
# --driver_gene_db /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/tables/gene/Compendium_Cancer_Genes.tsv \
# --output_info /home/lqh/Codes/Python/RNA-SSNV/output/GBM.RNA_only.driver.table

import pandas as pd

import argparse

# description参数可以用于描述脚本的参数作用，默认为空
parser=argparse.ArgumentParser(description="A discriminate model construction pipeline for RNA-SSNV.")
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.')
# Generic parameter
parser.add_argument('--result_info', help='Result info location')
parser.add_argument('--cancer_type', help='Cancer type')
parser.add_argument('--driver_gene_db', help='Database for cancer driver gene location')
parser.add_argument('--output_info', help='Output info location')

args = parser.parse_args()

# 获取驱动基因相关信息
def driver_gene_info_extract(driver_gene_loc, cancer_type):
    print("正在读取全部驱动基因信息")
    driver_genes = pd.read_table(driver_gene_loc)
    print(f"开始获取{cancer_type}所对应的驱动基因集合")
    cancer_spec_driver_genes = set(driver_genes[driver_genes.CANCER_TYPE == cancer_type]['SYMBOL'])
    print(f"{cancer_type}所对应的驱动基因集合获取完成，共计{len(cancer_spec_driver_genes)}个驱动基因，名称为cancer_spec_driver_genes")
    return cancer_spec_driver_genes

if __name__ == '__main__':
    cancer_spec_driver_genes = driver_gene_info_extract(args.driver_gene_db, args.cancer_type)
    # 读取突变文件
    result_info = pd.read_table(args.result_info)
    # 求驱动基因所对应子集
    output_info = result_info[result_info['Hugo_Symbol'].isin(cancer_spec_driver_genes)]
    # 仅输出指定列，保持前后一致性
    output_info = output_info[['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Allele1', "Hugo_Symbol",
                               'ref_AD_tumor_RNA', 'alt_AD_tumor_RNA', 'ref_AD_normal', 'alt_AD_normal',
                               'ref_AD_tumor_DNA', 'alt_AD_tumor_DNA',
                               "Tumor_Sample_UUID", 'Protein_Change', 'pred_prob', 'pred_label', 'Tag']]
    print(f"{len(result_info)}个突变中有{len(output_info)}个突变位于驱动基因内")
    output_info.to_csv(args.output_info, sep="\t", index=False)