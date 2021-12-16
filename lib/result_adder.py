# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_adder.py \
# --result_info /home/lqh/Codes/Python/RNA-SSNV/output/LUSC.final.table \
# --output_info /home/lqh/Codes/Python/RNA-SSNV/output/LUSC.final.DNA_coverage.table \
# --add_type DNA \
# --DNA_calling_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/tables/info/LUSC_WXS_somatic_calling_info.tsv \
# --DNA_tumor_folder /public1/data/projects/tumor/multi/TCGA/raw/WXS/LUSC \
# --num_threads 80


import pandas as pd
import numpy as np
import os
from multiprocessing import Pool
from multiprocessing import Manager
import pysam

import argparse

parser=argparse.ArgumentParser(description="Add specified information into given result table.")
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.')
# Generic parameter
parser.add_argument('--result_info', help='Result info location')
parser.add_argument('--output_info', help='Result info after adding infos')
# Specific parameter
parser.add_argument('--add_type', help='Type of info to add into result table.')
parser.add_argument("--DNA_calling_info", help="Tabular info for DNA somatic mutation calling.")
parser.add_argument('--DNA_tumor_folder', help='DNA tumor bam folder path.')
parser.add_argument('--num_threads', '-nt', type=int, help='Number of threads allowed.')

args = parser.parse_args()

# tool functions
def DNA_tumor_bam_record_retrive(case_result_row, samfile_list):
    # 根据vcf文件中单个record的突变信息，获取所有特定case的DNA肿瘤样本中var最大AD及其他相应信息，其返回值纳入其他信息的dataframe中
    try:
        # 获取基本信息
        ref = case_result_row['Reference_Allele']
        alt = case_result_row['Tumor_Allele1']
        chrom = case_result_row['Chromosome']
        position = case_result_row['Start_Position']
        # 获取每个bam文件中对应的碱基coverage信息
        other_coverage_list = []
        # 应用首个bam文件中coverage信息完成初始化(设置base的最低质量值为30)
        initial_count_coverage_info = count_coverage_decompose(samfile_list[0].count_coverage(contig=chrom, start=position-1, stop=position, quality_threshold=30))
        # 什么乱七八糟的数据结构啊= =我佛了，还需要强行转为str后才能作为key去访问字典中元素
        optimal_ref_coverage = initial_count_coverage_info[str(ref)]
        optimal_alt_coverage = initial_count_coverage_info[str(alt)]
        # 遍历所有bam文件，并进行条件判断
        for samfile in samfile_list[1:]:
            # pysam apply 0-based coordinate system(设置base的最低质量值为30)
            count_coverage_info = count_coverage_decompose(samfile.count_coverage(contig=chrom, start=position-1, stop=position, quality_threshold=30))
            # 当前覆盖度信息
            current_ref_coverage = count_coverage_info[str(ref)]
            current_alt_coverage = count_coverage_info[str(alt)]
            # 当前覆盖度信息进行对比
            # 比较alt的coverage信息
            if current_alt_coverage > optimal_alt_coverage:
                # 保存当前值，并将更好的值重新赋值
                other_coverage_list.append(str(optimal_ref_coverage)+"/"+str(optimal_alt_coverage))
                optimal_ref_coverage, optimal_alt_coverage = current_ref_coverage, current_alt_coverage
            elif current_alt_coverage == optimal_alt_coverage:
                # 进一步比较ref的coverage信息
                if current_ref_coverage > optimal_ref_coverage:
                    # 保存当前值，并将更好的值重新赋值
                    other_coverage_list.append(str(optimal_ref_coverage) + "/" + str(optimal_alt_coverage))
                    optimal_ref_coverage, optimal_alt_coverage = current_ref_coverage, current_alt_coverage
                else:
                    # 保存当前值，不进行更新
                    other_coverage_list.append(str(current_ref_coverage) + "/" + str(current_alt_coverage))
            else:
                # 保存当前值，不进行更新
                other_coverage_list.append(str(current_ref_coverage) + "/" + str(current_alt_coverage))
        return optimal_ref_coverage, optimal_alt_coverage, ";".join(other_coverage_list)
    except Exception as ex:
        print(f"For {case_result_row['Reference_Allele']}>{case_result_row['Tumor_Allele1']}, DNA_tumor_bam_record_retrive failed and emitted null values.")

        return np.NaN, np.NaN, ""

# 将pysam中count_coverage函数的返回值转换为字典形式
def count_coverage_decompose(count_coverage_info):
    try:
        count_coverage_dict = {}
        # 函数返回值为ACGT的顺序情况
        count_coverage_dict['A'] = count_coverage_info[0][0]
        count_coverage_dict['C'] = count_coverage_info[1][0]
        count_coverage_dict['G'] = count_coverage_info[2][0]
        count_coverage_dict['T'] = count_coverage_info[3][0]
        return count_coverage_dict
    except Exception as ex:
        print(ex)
        print("count_coverage_decompose错误！！！")

def print_error(value):
    print("error: ", value)

# main functions
def add_DNA_coverage_info(case_result_info, case_id, DNA_calling_info, DNA_tumor_folder_path, result_list):
    print(f"Start to retrieve DNA coverage info from {case_id} case...")
    # determine aliquots ids for tumor DNA sequence data
    DNA_tumor_case_df = DNA_calling_info.loc[
        (DNA_calling_info['case_id'] == case_id) & (DNA_calling_info['sample_type'] == "Primary Tumor"), ['file_id', 'file_name', 'aliquots_id']]
    DNA_tumor_case_file_paths = DNA_tumor_folder_path + "/" + DNA_tumor_case_df['file_id'] + "/" + DNA_tumor_case_df['file_name']
    # read in all case_id's corresponding DNA tumor bam files and store objects into a list
    samfile_list = [pysam.AlignmentFile(DNA_tumor_case_file_path, "rb") for DNA_tumor_case_file_path in DNA_tumor_case_file_paths]
    # retrieve all mutations' optimal DNA coverage info
    DNA_coverage_info = list()
    for i in case_result_info.index:
        DNA_coverage_info.append(list(DNA_tumor_bam_record_retrive(case_result_info.loc[i, ], samfile_list)))
    DNA_coverage_info = pd.DataFrame(DNA_coverage_info, columns=['ref_AD_tumor_DNA', 'alt_AD_tumor_DNA', 'other_AD_tumor_DNA'])
    # concat corresponding info back into result info
    case_result_info['ref_AD_tumor_DNA'] = list(DNA_coverage_info['ref_AD_tumor_DNA'])
    case_result_info['alt_AD_tumor_DNA'] = list(DNA_coverage_info['alt_AD_tumor_DNA'])
    case_result_info['other_AD_tumor_DNA'] = list(DNA_coverage_info['other_AD_tumor_DNA'])

    result_list.append(case_result_info)


if __name__ == '__main__':
    result_info = pd.read_table(args.result_info)
    num_threads = args.num_threads

    if args.add_type == "DNA":
        DNA_calling_info = pd.read_table(args.DNA_calling_info)
        DNA_tumor_folder_path = args.DNA_tumor_folder

        # 开始多进程处理
        print('Parent process %s.' % os.getpid())
        p = Pool(num_threads)
        manager = Manager()
        result_list = manager.list()
        # 遍历每一个case
        for case_id in result_info["Tumor_Sample_UUID"].value_counts().index:
            case_result_info = result_info.loc[result_info['Tumor_Sample_UUID']==case_id, ]
            p.apply_async(add_DNA_coverage_info, args=(
            case_result_info, case_id, DNA_calling_info, DNA_tumor_folder_path, result_list), error_callback=print_error)
        print('Waiting for all subprocesses done...')
        p.close()
        p.join()
        print('All subprocesses done.')

        result_list = list(result_list)
        # concat corresponding info back into result info
        output_info = pd.concat(result_list)

        output_info.to_csv(args.output_info, sep="\t", index=False)