# python /home/lqh/Codes/Python/RNA-SSNV/lib/Mutect2_calls_prepare_to_table.py \
# --cancer_type LUAD \
# --project_folder /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results \
# --RNA_calling_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/tables/info/LUAD_RNA_somatic_calling_info.tsv \
# --output_file_path /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/LUAD/RNA/RNA_somatic_mutation/VariantsToTable/SNP_WES_Interval_exon.table

import pandas as pd
import os

import argparse

# description参数可以用于描述脚本的参数作用，默认为空
parser=argparse.ArgumentParser(description="A discriminate model construction pipeline for RNA-SSNV.")
# Generic parameter
parser.add_argument('--cancer_type', type=str, help='Cancer type of current project.')
parser.add_argument("--project_folder", help="Path for current project.")
parser.add_argument("--RNA_calling_info", help="Tabular info for RNA somatic calling.")
parser.add_argument("--output_file_path", help="Path for output table.")

args=parser.parse_args()

cancer_type = args.cancer_type

# 获取项目case_id相应信息
tumor_tsv = pd.read_table(args.RNA_calling_info)
# 新建dataframe来保存Mutect2 WES突变位点信息
column_names = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele1", "Tumor_Sample_UUID"]
tumor_table_info = pd.DataFrame(columns=column_names)

# 遍历每一个case(直接用set()也可以= =)
for single_case_id in tumor_tsv["case_id"].value_counts().index:
    print(f"开始合并{single_case_id}对应信息")
    # 读取单case内突变信息
    case_table_info = pd.read_table(os.path.join(args.project_folder, args.cancer_type, f"RNA/RNA_somatic_mutation/VariantsToTable/SNP_WES_Interval_exon/{single_case_id}.table"))
    # 删除奇怪的列
    del case_table_info['Unnamed: 4']
    # 新建case列
    case_table_info['Tumor_Sample_UUID'] = single_case_id
    # 指定列名
    case_table_info.columns = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele1", "Tumor_Sample_UUID"]
    # 完成合并
    tumor_table_info = pd.concat([tumor_table_info, case_table_info], axis=0, ignore_index=True)

# 输出合并后结果
tumor_table_info.to_csv(args.output_file_path, sep="\t", index=False,columns=column_names)
