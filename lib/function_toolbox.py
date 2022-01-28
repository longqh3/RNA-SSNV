import pandas as pd

def multi_query(query_terms, data_info):
    data_bool = True
    for col in query_terms.keys():
        data_bool = data_info[col].map(lambda x: query_terms[col].__contains__(x)) & data_bool
    return data_info[data_bool]

# 预先提取相应ENSG信息列以转换为Hugo_Symbol
# UCSC_tpm = pd.read_table("/home/lqh/resources/database/UCSC_Xena/GBM/GDC/expression/TCGA-GBM.tpm.tsv")
# UCSC_tpm[['Ensembl_ID']].to_csv("/home/lqh/resources/database/UCSC_Xena/LUSC/GDC/expression/TCGA-LUSC.ENSG.tsv", header=False, index=False)
# UCSC_tpm_converter("/home/lqh/resources/database/UCSC_Xena/GBM/GDC/expression/TCGA-GBM.tpm.tsv", "/home/lqh/resources/database/UCSC_Xena/GBM/GDC/expression/TCGA-GBM.ENSG.tsv.trans.txt",
#                    "/home/lqh/resources/database/UCSC_Xena/GBM/GDC/expression/TCGA-GBM.tpm.Hugo_Symbol.tsv")
def UCSC_tpm_converter(UCSC_tpm_loc, transformed_symbol_loc, new_UCSC_tpm_loc):
    # 读取数据
    UCSC_tpm = pd.read_table(UCSC_tpm_loc)
    transformed_symbol = pd.read_table(transformed_symbol_loc)
    # 添加symbol列
    UCSC_tpm['Hugo_Symbol'] = transformed_symbol['Symbol']
    # 提取sample ID相关信息
    columns = list(UCSC_tpm.columns)
    columns.remove('Ensembl_ID')
    columns.remove('Hugo_Symbol')
    # 将sample ID转换为case ID相关的信息
    final_columns = []
    for sample_id in columns:
        sample_id_split = sample_id.split(".")
        if sample_id_split[-1] == ("01A") or sample_id_split[-1] == ("01B"):
            final_columns.append("-".join(sample_id_split[0:-1] + ["01A"]))
        elif sample_id_split[-1] == ("11A") or sample_id_split[-1] == ("11B"):
            final_columns.append("-".join(sample_id_split[0:-1] + ["11A"]))
    # 完整修改column信息
    UCSC_tpm.columns = final_columns+['Ensembl_ID', 'Hugo_Symbol']
    # 导出添加信息后表达量矩阵
    UCSC_tpm.to_csv(new_UCSC_tpm_loc, index=False, sep="\t")