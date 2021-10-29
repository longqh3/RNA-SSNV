# 2021.10.9 #
# 测试完成，一切正常#
# 测试命令 #
python /home/lqh/Codes/Python/RNA-SSNV/model_utilize.py \
--REDIportal /home/lqh/resources/database/RNA_edit/REDIportal/REDIportal_main_table.hg38.bed \
--DARNED /home/lqh/resources/database/RNA_edit/DARNED_hg19_to_bed_to_hg38_rm_alt.bed \
--raw_RNA_mutations /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/LUAD/RNA/RNA_somatic_mutation/VcfAssembly_new/SNP_WES_Interval_exon.txt \
--model_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer.model \
--one_hot_encoder_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer.one_hot_encoder \
--training_columns_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer.training_data_col \
--output_table_path /home/lqh/Codes/Python/RNA-SSNV/output/LUAD.table

# 导入相关需求包
# 基础包
import os
import math
import numpy as np    #导入Python科学计算的基础软件包numpy
import pandas as pd     #导入python的一个数据分析包pandas
# 数据分析包
## 模型包
from sklearn.feature_selection import RFECV
from sklearn import decomposition    #导入数据降维包decomposition，以便后面导入PCA包
from sklearn.ensemble import RandomForestClassifier  # 导入随机森林算法
## 模型相关包——数据预处理
from sklearn import preprocessing  # 导入数据预处理包
from sklearn.model_selection import train_test_split  # 导入训练集和测试集划分函数tain_test_split
## 模型相关包——导入模型训练相关优化函数
from sklearn.model_selection import StratifiedKFold     #导入将数据划分函数StratifiedKFold
from sklearn.model_selection import GridSearchCV    #导入网格搜索自动调参函数GridSearchCV
## 模型相关包——模型评估
from sklearn.metrics import *    #导入metrics模块的所有函数，metrics模块实现了一些函数，用来评估预测误差。已使用：precision_recall_curve
# 绘图包
import scikitplot as skplt # 绘制模型P-R曲线
import matplotlib.pyplot as plt    #导入Python可视化Matplotlib模块的绘图pyplot函数
import seaborn as sns    #导入Python数据可视化matplotlib模块扩展版seaborn模块
# 其他包
import pickle # 保存相应模型
import joblib # 其他方式保存模型
import warnings    #导入Python中的warnings模块
warnings.filterwarnings('ignore')    #warnings模块利用过滤器来实现忽略告警
# 导入模型解释包
import lime
import lime.lime_tabular
# 导入数据后处理相关包
import os
from multiprocessing import Pool
from multiprocessing import Manager
import portion as P
import vcf

import argparse

# description参数可以用于描述脚本的参数作用，默认为空
parser=argparse.ArgumentParser(description="A discriminate model construction pipeline for RNA-SSNV.")
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.')

parser.add_argument('--REDIportal', help='hg38 REDIportal bed file.')
parser.add_argument('--DARNED', help='hg38 DARNED bed file.')
parser.add_argument('--raw_RNA_mutations', '-r', help='raw RNA somatic single nucleotide variants.')
parser.add_argument("--model_path", help="Path for constructed model.")
parser.add_argument("--one_hot_encoder_path", help="Path for one-hot encoder.")
parser.add_argument("--training_columns_path", help="Path for constructed model.")
parser.add_argument("--output_table_path", help="Path for final output table.")

args=parser.parse_args()

# 编制碱基改变对应规则
ALLELE_CHANGE_DICT = {
    "T>C":"A>G",
    "C>T":"G>A",
    "T>G":"A>C",
    "G>T":"C>A",
    "T>A":"A>T",
    "G>C":"C>G"
}

def RNA_EDIT_process(REDIprotal, DARNED):
    # 获取REDIportal和DARNED数据库中RNA编辑位点信息
    REDIprotal_info = pd.read_table(REDIprotal, header=None)
    DARNED_info = pd.read_table(DARNED, header=None)
    ## 对数据库信息进行预处理，提取所需信息
    REDIprotal_info = REDIprotal_info[[0, 2]]
    DARNED_info = DARNED_info[[0, 2]]
    print(f"REDIportal和DARNED数据库中的RNA编辑位点数目分别为{len(REDIprotal_info)}和{len(DARNED_info)}")
    ## 对REDIpotal和DARNED数据库信息进行合并
    RNA_EDIT_INFO = pd.concat([REDIprotal_info, DARNED_info], ignore_index=True)
    RNA_EDIT_INFO.columns = ["Chromosome", "Start_Position"]
    ## 合并所得数据库信息需要去重
    RNA_EDIT_INFO = RNA_EDIT_INFO.drop_duplicates(keep="first")
    print(f"合并数据库信息并去重后，所得的RNA编辑位点总数为{len(RNA_EDIT_INFO)}")

    return RNA_EDIT_INFO

class model_utilize_no_DNA(object):

    # 初始化类实例
    # 输入1：所有原始RNA SSNV突变
    # 输入2：模型文件存储路径
    # 输入3：RNA编辑位点
    def __init__(self, independent_info, RNA_edit_info, model_path, one_hot_encoder_path, training_columns_path, output_table_path):
        # 初始化数据
        self.independent_info = independent_info
        self.RNA_edit_info = RNA_edit_info
        # 初始化模型
        self.rf_gcv = joblib.load(model_path)
        self.enc = joblib.load(one_hot_encoder_path)
        self.training_columns = pd.Series(pd.read_table(training_columns_path)["0"])
        self.output_table_path = output_table_path

    # 2021.07.05 根据上述最佳模型，绘制独立验证集中的PR曲线，并根据使F1最大化的阈值来输出相应positive、negative突变集
    # 应用目前效果最好的模型进行测试
    # 输入：实例新建的属性rf_gcv、scaler，代表经过网格搜索所得到的最优RF模型、标准化工具
    # 输入2：实例新建的属性independent_info,为独立测试数据集对应的dataframe
    # 输出：实例新建的属性independent_info_training、independent_info_training_scaler,为独立测试数据集对应的特征信息、标准化后的特征信息
    # 输出2：实例新建的属性independent_positive_thre，为独立测试数据集根据最大化F1值所得阈值对应的阳性位点（negative为阴性位点）
    # 输出3：实例新建的属性gdc_validate_df_SNP，为经过突变去重+SNP选取的gdc突变数据集
    # 输出4：实例新建的属性independent_pred、independent_pred_thre、independent_label，为根据相应模型和gdc突变数据集结果，所得独立验证集的预测概率值、预测标签和实际标签
    # 输出5：实例新建的属性all_info_final、GDC_SNP_info_final，最终确定用以构建标签的RNA突变数据集和GDC突变数据集
    def common_RF_utilize(self):

        # 开始对所有RNA突变进行筛选，筛选掉基本无法处理&没必要处理的multiallelic位点——self.all_info
        print(f"所有RNA突变在进行multi-allelic筛选前，对应的突变数目为{len(self.independent_info)}")
        self.independent_info = self.independent_info[self.independent_info['record_filter'].apply(lambda x: False if x.__contains__("multiallelic") else True)]
        print(f"所有RNA突变在经过multi-allelic筛选后，剩余的突变数目为{len(self.independent_info)}")

        print("="*100)

        # 进一步对RNA编辑位点数据库进行分析，判断其是否位于数据库内
        ## 找出位于RNA编辑数据库中的相应RNA突变位点——保存于self.all_info_RNA_edit中
        self.independent_info_RNA_edit = pd.merge(self.independent_info, self.RNA_edit_info)
        self.independent_info = self.independent_info.append(self.independent_info_RNA_edit)
        self.independent_info = self.independent_info.drop_duplicates(keep=False)

        print(f"位于RNA编辑位点上的突变数为{len(self.independent_info_RNA_edit)}（对象名称为independent_info_RNA_edit），"
              f"不属于RNA编辑位点的突变数为{len(self.independent_info)}")

        print("="*100)

        # 进一步排除位于免疫球蛋白基因内突变
        self.independent_info_immunoglobulin = self.independent_info[self.independent_info.apply(lambda x: True if x['Hugo_Symbol'].startswith(("IGH", 'IGK', 'IGL')) else False, axis=1)]
        self.independent_info = self.independent_info[self.independent_info.apply(lambda x: False if x['Hugo_Symbol'].startswith(("IGH", 'IGK', 'IGL')) else True, axis=1)]

        print(f"位于免疫球蛋白基因内的突变数为{len(self.independent_info_immunoglobulin)}（对象名称为independent_info_immunoglobulin），"
              f"不位于免疫球蛋白基因内的突变数为{len(self.independent_info)}")

        print("="*100)

        # 进一步排除HLA基因内突变
        self.independent_info_HLA = self.independent_info[self.independent_info.apply(lambda x: True if x['Hugo_Symbol'].startswith("HLA") else False, axis=1)]
        self.independent_info = self.independent_info[self.independent_info.apply(lambda x: False if x['Hugo_Symbol'].startswith("HLA") else True, axis=1)]

        print(f"位于HLA基因内的突变数为{len(self.independent_info_HLA)}（对象名称为independent_info_HLA），"
              f"不位于HLA基因内的突变数为{len(self.independent_info)}")

        print("=" * 100)

        ### 将index重置为连续状态（避免添加碱基改变标记时，index不一致出现问题）
        self.independent_info = self.independent_info.reset_index(drop=True)

        ### 2020.12.7 考虑添加碱基改变信息后，构建模型判断效果改善与否
        independent_info_allele_change = self.independent_info['Reference_Allele'] + '>' + self.independent_info['Tumor_Allele1']
        #### 2021.5.16——修改突变碱基改变信息
        independent_info_allele_change = pd.Series([ALLELE_CHANGE_DICT[allele_change] if ALLELE_CHANGE_DICT.keys().__contains__(allele_change) else allele_change for allele_change in independent_info_allele_change])

        independent_info_allele_change_df = pd.DataFrame(self.enc.transform(independent_info_allele_change.values.reshape(-1, 1)).toarray(), columns=self.enc.categories_[0])
        self.independent_info = pd.concat([self.independent_info, independent_info_allele_change_df], axis=1)

        # 数据预处理
        # 将类别列和其他无关列（CHROM、POS、REF、ALT、FILTER、CONTQ）从数据集删除，并指定列名顺序（避免预测错误）
        # self.independent_info_training = self.independent_info.drop(DROP_COLUMNS, axis=1)
        self.independent_info_training = self.independent_info[self.training_columns]
        self.independent_info_training = self.independent_info_training.dropna()
        # 缺值删除——仅对需要纳入模型的列进行处理——删除包含缺失值的记录，以保持后续切片过程中，行数的一致
        self.independent_info = self.independent_info.dropna(subset=self.independent_info_training.columns)

        # 数据导出
        ## 概率、标签共同导出
        self.independent_pred = self.rf_gcv.predict_proba(self.independent_info_training)[:, 1]  # 概率形式预测测试集的类别
        self.independent_info['pred_prob'] = self.independent_pred
        self.independent_pred_thre = self.rf_gcv.predict(self.independent_info_training)  # 标签形式预测测试集的类别
        self.independent_info['pred_label'] = self.independent_pred_thre
        print(f"预测完成，positive数目为{self.independent_info['pred_label'].value_counts()[1]}，negative数目为{self.independent_info['pred_label'].value_counts()[0]}")
        ## 导出至指定路径
        self.independent_info.to_csv(self.output_table_path, sep="\t", index=False)

if __name__ == '__main__':
    independent_info = pd.read_table(args.raw_RNA_mutations)
    RNA_EDIT_INFO = RNA_EDIT_process(args.REDIportal, args.DARNED)

    model_utilize = model_utilize_no_DNA(independent_info, RNA_EDIT_INFO, args.model_path, args.one_hot_encoder_path, args.training_columns_path, args.output_table_path)
    model_utilize.common_RF_utilize()