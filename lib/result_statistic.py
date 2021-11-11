# test successful
# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_statistic.py \
# --RNA_DNA_overlap_info /home/lqh/Codes/Python/RNA-SSNV/data/scores/LUSC.RNA_DNA_overlap.all.txt \
# --DNA_only_info /home/lqh/Codes/Python/RNA-SSNV/data/scores/LUSC.DNA_only.all.txt \
# --test_type chisq \
# --test_column_name SIFT_converted_rankscore \
# --column_threshold 0.39575
#
# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_statistic.py \
# --RNA_DNA_overlap_info /home/lqh/Codes/Python/RNA-SSNV/data/scores/GBM.RNA_DNA_overlap.all.txt \
# --DNA_only_info /home/lqh/Codes/Python/RNA-SSNV/data/scores/GBM.DNA_only.all.txt \
# --test_type chisq \
# --test_column_name FATHMM_converted_rankscore \
# --column_threshold 0.81332
#
# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_statistic.py \
# --RNA_DNA_overlap_info /home/lqh/Codes/Python/RNA-SSNV/data/scores/GBM.RNA_DNA_overlap.all.txt \
# --DNA_only_info /home/lqh/Codes/Python/RNA-SSNV/data/scores/GBM.DNA_only.all.txt \
# --test_type chisq \
# --test_column_name CADD_raw_rankscore \
# --column_threshold 0.5

import pandas as pd
from scipy.stats import chi2_contingency
import numpy as np
import statsmodels.api as sm

import argparse

parser=argparse.ArgumentParser(description="Assessment of pathogenic scores for RNA-SSNV's missense results.")
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.')
# Generic parameter
parser.add_argument('--RNA_DNA_overlap_info', help='RNA_DNA_overlap info location')
parser.add_argument('--DNA_only_info', help='DNA_only info location')
parser.add_argument('--test_type', help='Statistical test type.')
parser.add_argument('--test_column_name', help='Column name for conducting statistical test')
parser.add_argument('--column_threshold', type=float, help='Threshold for specified column')

args = parser.parse_args()

if __name__ == '__main__':
    # 读取突变文件
    RNA_DNA_overlap_info = pd.read_table(args.RNA_DNA_overlap_info)
    DNA_only_info = pd.read_table(args.DNA_only_info)
    # 根据检验类型进行后续操作
    if args.test_type == "chisq":
        # 根据指定阈值转换指定列值
        RNA_DNA_overlap_status = RNA_DNA_overlap_info[args.test_column_name].apply(lambda x : 1 if x >= args.column_threshold else 0).value_counts()
        DNA_only_status = DNA_only_info[args.test_column_name].apply(lambda x : 1 if x >= args.column_threshold else 0).value_counts()
        print("Distribution for RNA_DNA_overlap_status: ")
        print(RNA_DNA_overlap_status)
        print("Distribution for DNA_only_status: ")
        print(DNA_only_status)
        # 对指定列进行卡方检验
        kf_data = np.array([[RNA_DNA_overlap_status[1], RNA_DNA_overlap_status[0]],
                                [DNA_only_status[1], DNA_only_status[0]]])
        kf = chi2_contingency(kf_data)
        print('chisq-statistic=%.4f, p-value=%.4f, df=%i expected_frep=%s' % kf)
        # OR值计算
        contingency_table = sm.stats.Table2x2(kf_data)
        print("Summary information is")
        print(contingency_table.summary())
        print(f"OR value is {contingency_table.oddsratio}, p value is {contingency_table.oddsratio_pvalue()}")