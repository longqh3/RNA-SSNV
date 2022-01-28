# test success
# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_piler.py \
# --result_info /home/lqh/Codes/Python/RNA-SSNV/output/GBM.DNA_only.all.table \
# --result_info /home/lqh/Codes/Python/RNA-SSNV/output/GBM.RNA_DNA_overlap.all.table \
# --output_info /home/lqh/Codes/Python/RNA-SSNV/output/GBM.DNA_total.table

# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_piler.py \
# --result_info /home/lqh/Codes/Python/RNA-SSNV/output/GBM.DNA_only.all.table \
# --result_info /home/lqh/Codes/Python/RNA-SSNV/output/GBM.RNA_DNA_overlap.all.table \
# --output_info /home/lqh/Codes/Python/RNA-SSNV/output/GBM.DNA_total.driver_gene.table \
# --Intogen_database_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/tables/gene/Compendium_Cancer_Genes.tsv \
# --cancer_type GBM

# --pred_label 1 \

import pandas as pd

import argparse

parser=argparse.ArgumentParser(description="Pile up different results of RNA-SSNV into a total table "
                                           "and add other informations.")
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.')
# Generic parameter
parser.add_argument('--result_info', action='append', help='Result info location')
parser.add_argument('--output_info', help='Tags required')
# Special parameter
parser.add_argument('--Intogen_database_info', help='Intogen database info location.')
parser.add_argument('--cancer_type', help='Type of cancer.')

args = parser.parse_args()

# 获取驱动基因相关信息
def driver_gene_info_extract(driver_gene_loc, cancer_type):
    print("正在读取全部驱动基因信息")
    driver_genes = pd.read_table(driver_gene_loc)
    print(f"开始获取{cancer_type}所对应的驱动基因集合")
    cancer_spec_driver_genes = set(driver_genes[driver_genes.CANCER_TYPE == cancer_type]['SYMBOL'])
    print(f"{cancer_type}所对应的驱动基因集合获取完成，共计{len(cancer_spec_driver_genes)}个驱动基因，名称为cancer_spec_driver_genes")
    return(cancer_spec_driver_genes)

if __name__ == '__main__':

    print(f"待合并结果为{args.result_info}")

    output_info = pd.concat([pd.read_table(result_info) for result_info in args.result_info], axis=0)

    if args.cancer_type is not None:
        print(f"取驱动基因内突变前，突变数目为{len(output_info)}")
        cancer_spec_driver_genes = driver_gene_info_extract(args.Intogen_database_info, args.cancer_type)
        print(f"添加/修改Hugo_Symbol筛选标准为{args.cancer_type}对应的{len(cancer_spec_driver_genes)}个驱动基因")
        output_info = output_info[output_info.Hugo_Symbol.isin(cancer_spec_driver_genes)]
        print(f"取驱动基因内突变后，突变数目为{len(output_info)}")

    output_info.to_csv(args.output_info, sep="\t", index=False)