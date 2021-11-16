# test successful
# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_editor.py \
# --result_info /home/lqh/Codes/Python/RNA-SSNV/output/GBM.RNA_DNA_overlap.driver_gene.table \
# --edit_type oncoKB \
# --output_info /home/lqh/Codes/Python/RNA-SSNV/data/therapeutic/GBM.RNA_DNA_overlap.driver_gene.maf
#
# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_editor.py \
# --result_info /home/lqh/Codes/Python/RNA-SSNV/output/GBM.DNA_only.driver_gene.table \
# --edit_type oncoKB \
# --output_info /home/lqh/Codes/Python/RNA-SSNV/data/therapeutic/GBM.DNA_only.driver_gene.maf
#
# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_editor.py \
# --result_info /home/lqh/Codes/Python/RNA-SSNV/output/LUSC.RNA_only.driver_gene.table \
# --edit_type oncoKB \
# --output_info /home/lqh/Codes/Python/RNA-SSNV/data/therapeutic/LUSC.RNA_only.driver_gene.maf

import pandas as pd
from scipy.stats import chi2_contingency
import numpy as np
import statsmodels.api as sm

import argparse

parser=argparse.ArgumentParser(description="Edit result info for special needs.")
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.')
# Generic parameter
parser.add_argument('--result_info', help='Specified result info location')
parser.add_argument('--edit_type', help='Type of data to be edited into.')
parser.add_argument('--output_info', help='Specified output location')

args = parser.parse_args()

if __name__ == '__main__':
    # read in data info
    result_info = pd.read_table(args.result_info)
    # edit result accordingly
    if args.edit_type == "oncoKB":
        # new column
        result_info['End_Position'] = result_info['Start_Position']
        result_info['Tumor_Allele2'] = result_info['Tumor_Allele1']
        # modify column names
        result_info.rename(columns={'Chromosome': 'CHROMOSOME', 'Start_Position': 'START_POSITION', 'End_Position': 'END_POSITION',
                                    'Reference_Allele': 'REFERENCE_ALLELE', 'Tumor_Allele1':'TUMOR_SEQ_ALLELE1',
                                    'Tumor_Allele2':'TUMOR_SEQ_ALLELE2'}, inplace=True)
    # export result
    result_info.to_csv(args.output_info, sep="\t", index=False)