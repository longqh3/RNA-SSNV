import pandas as pd

cancer_type = "BLCA"

# 获取项目case_id相应信息
tumor_tsv = pd.read_table(f"/home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/tables/info/{cancer_type}_RNA_somatic_calling_info.tsv")
# 新建dataframe来保存Mutect2 WES突变位点信息
column_names = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele1", "Tumor_Sample_UUID"]
tumor_wes_info = pd.DataFrame(columns=column_names)

# 遍历每一个case(直接用set()也可以= =)
for single_case_id in tumor_tsv["case_id"].value_counts().index:
    print(f"开始合并{single_case_id}对应信息")
    # 读取单case内突变信息
    case_wes_info = pd.read_table(f"/home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/{cancer_type}/RNA/RNA_somatic_mutation/VariantsToTable/SNP_WES_Interval_exon/{single_case_id}.table")
    # 删除奇怪的列
    del case_wes_info['Unnamed: 4']
    # 新建case列
    case_wes_info['Tumor_Sample_UUID'] = single_case_id
    # 指定列名
    case_wes_info.columns = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele1", "Tumor_Sample_UUID"]
    # 完成合并
    tumor_wes_info = pd.concat([tumor_wes_info, case_wes_info], axis=0, ignore_index=True)

# 输出合并后结果
tumor_wes_info.to_csv(f"/home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/{cancer_type}/RNA/RNA_somatic_mutation/VariantsToTable/SNP_WES_Interval_exon.table", sep="\t", index=False,columns=column_names)
