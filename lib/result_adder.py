# test success
# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_adder.py \
# --result_info /home/lqh/Codes/Python/RNA-SSNV/output/GBM.DNA_only.table \
# --result_info /home/lqh/Codes/Python/RNA-SSNV/output/GBM.RNA_DNA_overlap.table \
# --output_info /home/lqh/Codes/Python/RNA-SSNV/output/GBM.DNA_total.table

# --pred_label 1 \

import pandas as pd

import argparse

parser=argparse.ArgumentParser(description="Pile up different results of RNA-SSNV into a total table "
                                           "and add other informations.")
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.')
# Generic parameter
parser.add_argument('--result_info', action='append', help='Result info location')
parser.add_argument('--output_info', help='Tags required')
#
args = parser.parse_args()

print(f"待合并结果为{args.result_info}")

output_info = pd.concat([pd.read_table(result_info) for result_info in args.result_info], axis=0)
output_info.to_csv(args.output_info, sep="\t", index=False)