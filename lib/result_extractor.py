# 结果正确性验证完成
# Tag: RNA_DNA_overlap DNA_only RNA_only
# pred_label: 1 0
# output_type: simple_maf total_maf
# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_extractor.py \
# --result_info /home/lqh/Codes/Python/RNA-SSNV/output/BLCA.final.table \
# --Tag RNA_DNA_overlap \
# --output_type total_maf \
# --pred_label 1 \
# --Hugo_Symbol TP53 \
# --output_info /home/lqh/Codes/Python/RNA-SSNV/output/BLCA.RNA_DNA_overlap.TP53.table

# --pred_label 1 \
# --Hugo_Symbol KMT2C \

# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_extractor.py \
# --result_info /home/lqh/Codes/Python/RNA-SSNV/output/LUSC.final.table \
# --Tag DNA_only \
# --Sub_Tag force_called \
# --Variant_Classification Missense_Mutation \
# --output_type total_maf \
# --output_info /home/lqh/Codes/Python/RNA-SSNV/output/LUSC.DNA_only.driver_gene.table \
# --Intogen_database_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/tables/gene/Compendium_Cancer_Genes.tsv \
# --cancer_type LUSC

# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_extractor.py \
# --result_info /home/lqh/Codes/Python/RNA-SSNV/output/LUSC.final.table \
# --Tag RNA_DNA_overlap \
# --Variant_Classification Missense_Mutation \
# --pred_label 1 \
# --output_type total_maf \
# --output_info /home/lqh/Codes/Python/RNA-SSNV/output/LUSC.RNA_DNA_overlap.driver_gene.table \
# --Intogen_database_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/tables/gene/Compendium_Cancer_Genes.tsv \
# --cancer_type LUSC

import pandas as pd

import argparse

# description参数可以用于描述脚本的参数作用，默认为空
parser=argparse.ArgumentParser(description="A discriminate model construction pipeline for RNA-SSNV.")
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.')
# Generic parameter
parser.add_argument('--result_info', help='Result info location')
parser.add_argument('--Tag', help='Tags required', nargs='+')
parser.add_argument("--Sub_Tag", help="Sub_Tag required.", nargs='+')
parser.add_argument("--pred_label", help="Predicted label required.", nargs='+', type=int)
parser.add_argument('--Hugo_Symbol', help='Hugo Symbols required.', nargs='+')
parser.add_argument('--Variant_Classification', help='Variant Classification required.', nargs='+')
parser.add_argument('--output_type', help='Output type: simple_maf, total_maf')
parser.add_argument('--output_info', help='Tags required')
# Specific parameters
parser.add_argument('--Intogen_database_info', help='Intogen database info location.')
parser.add_argument('--cancer_type', help='Intogen database info location.')

args = parser.parse_args()

# 获取驱动基因相关信息
def driver_gene_info_extract(driver_gene_loc, cancer_type):
    print("正在读取全部驱动基因信息")
    driver_genes = pd.read_table(driver_gene_loc)
    print(f"开始获取{cancer_type}所对应的驱动基因集合")
    cancer_spec_driver_genes = set(driver_genes[driver_genes.CANCER_TYPE == cancer_type]['SYMBOL'])
    print(f"{cancer_type}所对应的驱动基因集合获取完成，共计{len(cancer_spec_driver_genes)}个驱动基因，名称为cancer_spec_driver_genes")
    return(cancer_spec_driver_genes)

# 似乎也没有必要自己来实现excel的功能，暂且搁置
# 工具函数7：根据给定条件（Hugo_Symbol、Tumor_Sample_UUID、Start_Position等）来进行筛选
# 目前暂时通过输入字典的方式，来进行maf信息的筛选分析
# 应用map方法，针对不同提示信息完成筛选
# Chromosome、Start_Position、Hugo_Symbol、Tumor_Sample_UUID、Tag、Sub_Tag、Predicted_prob、Predicted_label均为潜在筛选标准
def multiple_query(data_info, query_terms, output_type):
    data_bool = True
    for col in query_terms.keys():
        data_bool = data_info[col].map(lambda x: query_terms[col].__contains__(x)) & data_bool
    if output_type=="simple_maf":
        # 提取maftools分析所需的相关列信息
        data_info = data_info[['Hugo_Symbol', "Chromosome", "Start_Position", "Start_Position",
             'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Classification',
             'Tumor_Sample_UUID', 'Protein_Change']]
        # 修改列名，使之适配相应临床数据中命名
        data_info.columns = ['Hugo_Symbol', "Chromosome", "Start_Position", "End_Position",
                            'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Classification',
                            'Tumor_Sample_Barcode', 'Protein_Change']
        # 添加新列
        data_info['Variant_Type'] = "SNP"
    result = data_info[data_bool]
    print(f"输出突变数为{len(result)}")
    return result

# 获取所有参数作为筛选标准
arg_terms = args.__dict__
query_terms = {}
# 挑选出有值的筛选标准来进行筛选
for key in ["Tag", "Sub_Tag", "pred_label", "Hugo_Symbol"]:
    if not (arg_terms[key] is None):
        query_terms[key] = arg_terms[key]
print(f"Select criteria: {query_terms}")
# 读取待筛选dataframe
data_info = pd.read_table(args.result_info)
# 判断是否需要求突变子集
if args.cancer_type is not None:
    cancer_spec_driver_genes = driver_gene_info_extract(args.Intogen_database_info, args.cancer_type)
    print(f"添加/修改Hugo_Symbol筛选标准为{args.cancer_type}对应的{len(cancer_spec_driver_genes)}个驱动基因")
    query_terms["Hugo_Symbol"] = cancer_spec_driver_genes

# 导出筛选后dataframe
multiple_query(data_info, query_terms, args.output_type).to_csv(args.output_info, sep="\t", index=False)
