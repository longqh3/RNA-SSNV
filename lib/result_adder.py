# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_adder.py \
# --result_info /home/lqh/Codes/Python/RNA-SSNV/output/GBM.DNA_only.table \
# --output_info /home/lqh/Codes/Python/RNA-SSNV/output/GBM.DNA_total.table

# --pred_label 1 \

import pandas as pd
import os
from multiprocessing import Pool
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

args = parser.parse_args()

# tool functions
# 根据vcf文件中单个record的突变信息，获取所有特定case的DNA肿瘤样本中var最大AD及其他相应信息，其返回值纳入其他信息的dataframe中
def DNA_tumor_bam_record_retrive(record, samfile_list):
    try:
        # 获取基本信息
        ref = record.REF
        alt = record.ALT[0]
        chrom = record.CHROM
        position = record.POS
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
        print(ex)
        print("DNA_tumor_bam_record_retrive错误！！！")

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

# main functions


if __name__ == '__main__':
    result_info = pd.read_table(args.result_info)

    if args.add_type == "DNA":
        # 开始多进程处理
        print('Parent process %s.' % os.getpid())
        p = Pool(args.num_threads)
        # 遍历每一个case(直接用set()也可以= =)
        for case_id in result_info["Tumor_Sample_UUID"].value_counts().index:
            case_result_info = result_info.loc[result_info['Tumor_Sample_UUID']==case_id, ]
            p.apply_async(case_vcf_info_retrive, args=(
            tumor_tsv, single_case_id, vcf_folder_path, table_folder_path, DNA_tumor_tsv, DNA_tumor_folder,
            Exon_loc,
            funcotator_folder_path))
        print('Waiting for all subprocesses done...')
        p.close()
        p.join()
        print('All subprocesses done.')

output_info.to_csv(args.output_info, sep="\t", index=False)