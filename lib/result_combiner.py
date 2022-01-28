# 与DoCM数据库数据结合分析
# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_combiner.py \
# --combine_type DoCM \
# --result_info /home/lqh/Codes/Python/RNA-SSNV/output/GBM.DNA_only.driver_gene.table \
# --DoCM_info /home/lqh/resources/database/DoCM/variants_hg38.tsv \
# --cancer_type GBM \
# --output_info /home/lqh/Codes/Python/RNA-SSNV/output/GBM.DNA_only.driver.DoCM.table

# 与TMB数据结合分析
# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_combiner.py \
# --combine_type TMB \
# --mutation_table_info /home/lqh/Codes/Python/RNA-SSNV/output/BLCA.DNA_only.LRP1B.table \
# --TMB_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/tables/cBioportal/BLCA/data_clinical_sample.txt \
# --mutation_TMB_info /home/lqh/Codes/Python/RNA-SSNV/results/fig4.BLCA.LRP1B.mutation.TMB.table \
# --no_mutation_TMB_info /home/lqh/Codes/Python/RNA-SSNV/results/fig4.BLCA.no.LRP1B.mutation.TMB.table

import pandas as pd

import argparse

parser=argparse.ArgumentParser(description="Integrate other database's info into our results.")
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.')
# Generic parameter
parser.add_argument('--combine_type', help='Type of data to combine')
# DoCM data type
parser.add_argument('--result_info', help='Result info location')
parser.add_argument('--DoCM_info', help='DoCM table location')
parser.add_argument('--output_info', help='Output info location')
parser.add_argument('--cancer_type', help='Type of cancer to query DoCM records')
# TMB data type
parser.add_argument('--mutation_table_info', help='Gene corresponding mutation table location')
parser.add_argument('--TMB_info', help='TMB table location')
parser.add_argument('--mutation_TMB_info', help='Output TMB info for cases harboring mutations')
parser.add_argument('--no_mutation_TMB_info', help='Output TMB info for cases without mutation')


args = parser.parse_args()

if __name__ == '__main__':

    if args.combine_type=="DoCM":
        # 读取待处理的突变集
        result_info = pd.read_table(args.result_info)
        print(f"DoCM突变取子集前，位于驱动基因内的突变数为{len(result_info)}")
        # 读取并处理DoCM数据集
        DoCM_info = pd.read_table(args.DoCM_info)
        # 根据指定字典内信息来求DoCM数据集子集
        query_dict = {
            "LUSC": ['cancer', 'lung squamous cell carcinoma', 'non-small cell lung carcinoma'],
            "BLCA": ['cancer', 'urinary bladder urothelial carcinoma'],
            "GBM": ['cancer', 'glioblastoma multiforme']
        }
        cancer_requried = DoCM_info['Cancer_Type'].map(lambda x: len(set(x).intersection(set(query_dict[args.cancer_type])))>0)
        cancer_DoCM_info = DoCM_info[cancer_requried]
        cancer_DoCM_info.drop_duplicates(subset=["Chromosome", 'Start_Position', 'Reference_Allele', "Tumor_Allele1"], keep='first', inplace=True)
        # 完成待处理的突变集和DoCM数据集的合并
        output_info = pd.merge(result_info, DoCM_info)
        print(f"DoCM突变取子集后，属于DoCM数据库的突变数为{len(output_info)}")
        # 仅输出指定列信息，便于导出至excel
        output_info = output_info[['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Allele1', "Hugo_Symbol",
                'ref_AD_tumor_RNA', 'alt_AD_tumor_RNA', 'ref_AD_normal', 'alt_AD_normal', 'ref_AD_tumor_DNA', 'alt_AD_tumor_DNA',
                "Tumor_Sample_UUID", 'Protein_Change', 'pred_prob', 'pred_label', 'Tag', 'Cancer_Type', 'Pubmed_Sources']]
        output_info.to_csv(args.output_info, sep="\t", index=False)
    elif args.combine_type=="TMB":
        mutation_table_info = pd.read_table(args.mutation_table_info)
        TMB_info = pd.read_table(args.TMB_info, comment="#")
        # 去除na
        print(f"去除NA值前，共有{len(TMB_info)}个case信息")
        TMB_info = TMB_info.dropna(subset=['TMB_NONSYNONYMOUS'])
        print(f"去除NA值后，共有{len(TMB_info)}个case信息")
        # 挑出突变基因对应所有突变信息
        mutation_case = mutation_table_info[mutation_table_info.Variant_Classification.isin(
            ["Missense_Mutation", "Nonsense_Mutation"])].Tumor_Sample_UUID.value_counts().index
        # 挑出含有突变基因的case
        mutation_TMB_info = TMB_info[TMB_info.PATIENT_ID.isin(mutation_case)]
        no_mutation_TMB_info = TMB_info[~TMB_info.PATIENT_ID.isin(mutation_case)]
        # 输出基本信息
        print(f"携带有该突变的case对应的TMB中位数为{mutation_TMB_info.TMB_NONSYNONYMOUS.median()}、均值为{mutation_TMB_info.TMB_NONSYNONYMOUS.mean()}")
        print(f"不携带该突变的case对应的TMB中位数为{no_mutation_TMB_info.TMB_NONSYNONYMOUS.median()}、均值为{no_mutation_TMB_info.TMB_NONSYNONYMOUS.mean()}")
        # 导出相应信息
        mutation_TMB_info[['PATIENT_ID', 'TMB_NONSYNONYMOUS']].to_csv(args.mutation_TMB_info, sep="\t", index=False)
        no_mutation_TMB_info[['PATIENT_ID', 'TMB_NONSYNONYMOUS']].to_csv(args.no_mutation_TMB_info, sep="\t", index=False)