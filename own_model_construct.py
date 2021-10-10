# 2021.10.8 #
# 测试完成，一切正常#
# 测试命令 #
# python /home/lqh/Codes/Python/RNA-SSNV/own_model_construct.py \
# --REDIportal /home/lqh/resources/database/RNA_edit/REDIportal/REDIportal_main_table.hg38.bed \
# --DARNED /home/lqh/resources/database/RNA_edit/DARNED_hg19_to_bed_to_hg38_rm_alt.bed \
# --raw_RNA_mutations /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/LUAD/RNA/RNA_somatic_mutation/VcfAssembly_new/SNP_WES_Interval_exon.txt \
# --DNA_mutations /home/lqh/Codes/Data/TCGA_maf_files/TCGA-LUAD \
# --model_folder_path /home/lqh/Codes/Python/RNA-SSNV/model

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


import argparse

# description参数可以用于描述脚本的参数作用，默认为空
parser=argparse.ArgumentParser(description="A discriminate model construction pipeline for RNA-SSNV.")
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.')

parser.add_argument('--REDIportal', help='hg38 REDIportal bed file.')
parser.add_argument('--DARNED', help='hg38 DARNED bed file.')
parser.add_argument('--raw_RNA_mutations', '-r', help='raw RNA somatic single nucleotide variants.')
parser.add_argument("--DNA_mutations", help="GDC mutations.")
parser.add_argument("--model_folder_path", help="Folder path for constructed model.")

args=parser.parse_args()

# 固定相关需从训练数据中drop的变量
# 删除"Expression_TPM"特征，因其在RNA测序数据中难以获得（需要额外进行软件运算），且事实上影响仅为0.1%
# 保存"ref_AD_normal", "alt_AD_normal"（在"ref_AD_tumor_DNA", "alt_AD_tumor_DNA"之前）
# 删除"ref_context"列，避免干扰结果（2021.1.11）
# 删除多个后添加的列，避免干扰结果构建
# 纳入"exon_distance"列，增强其对于不恰当的位点识别的敏感性（2021.02.20）
# 删除"exon_distace"、"target_symbol"列，以保证仅有"target_distance"列纳入模型中（2021.03.07）
# 进一步纳入UCSC GENCODE v22中更新后的"exon_distance"列，并排除其他相关列"transcript_ID"（2021.03.09）如"target_symbol"、"gene_symbol", "gene_strand", "gene_exon", "mutation_region",
# 2021.5.7 尝试删除"AF_tumor"这一特征（避免部分alt AD很高的位点被误认为是FN）
# 2021.5.9 尝试删除"DP"这一特征（避免DNA和RNA对应DP的无理由组合）——尝试失败，效果更差了
# 2021.5.10 尝试删除"COSMIC"这一特征（避免泛化性差这一可能后果）
# 2021.6.13 尝试删除"AS_UNIQ_ALT_READ_COUNT"这一特征（根据GATK官方，其仅在罕见场合可应用）
# 2021.8.17
DROP_COLUMNS = ['Attribute', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Allele1', 'Tumor_Allele2',
             'Tumor_Sample_UUID', "ref_AD_tumor_DNA", "alt_AD_tumor_DNA",
             "AD_other_RNA_tumor", "AD_other_normal", "AD_other_DNA_tumor",
             "transcript_ID",
             "Tumor_Seq_Allele2", "Hugo_Symbol", "Variant_Classification",
             "Gencode_28_secondaryVariantClassification", "HGNC_Ensembl_Gene_ID",
             "Reference_Allele", "ref_context",
             'Strand', 'Transcript_Strand', 'Codon_Change', 'Protein_Change', 'DrugBank',
             'record_filter',
             'COSMIC_total_alterations_in_gene',
             'AS_UNIQ_ALT_READ_COUNT']

# 编制碱基改变对应规则
ALLELE_CHANGE_DICT = {
    "T>C":"A>G",
    "C>T":"G>A",
    "T>G":"A>C",
    "G>T":"C>A",
    "T>A":"A>T",
    "G>C":"C>G"
}


# 提取GDC中所有突变的位点信息——Tumor_Sample_UUID Chromosome Start_Position Reference_Allele Tumor_Allele1 Tumor_Allele2
# 位点层面GDC信息GDC_SNP_info（不包含相应注释，无法用以maftools分析）
def GDC_site_info_retrieve(GDC_info):
    ## 整理GDC相应信息
    print("GDC项目中的突变位点情况为：\n", GDC_info["Variant_Type"].value_counts())
    ### 仅选择SNP突变信息进行分析
    GDC_SNP_info = GDC_info[GDC_info['Variant_Type'] == "SNP"]
    ### 重新编制index信息
    GDC_SNP_info.reset_index(drop=True, inplace=True)
    ### 将“Tumor_Sample_UUID”列中信息进行切片&提取，使之与RNA体细胞突变中的case_id信息保持一致
    del GDC_SNP_info['Tumor_Sample_UUID']
    GDC_SNP_info["Tumor_Sample_UUID"] = pd.Series(["-".join(GDC_sample_info.split("-")[0:3])
                                                        for GDC_sample_info in
                                                        GDC_SNP_info["Tumor_Sample_Barcode"]])
    ### 仅选取染色体、碱基和case_id信息来进行合并
    GDC_SNP_info = GDC_SNP_info.loc[:,
                        ["Chromosome", "Start_Position", "Reference_Allele", "Reference_Allele", "Tumor_Seq_Allele2",
                         "Tumor_Sample_UUID"]]
    ### 重新命名，便于进行合并(Tumor_Allele1对应突变碱基、Tumor_Allele2对应参考碱基)
    GDC_SNP_info.columns = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele2", "Tumor_Allele1", "Tumor_Sample_UUID"]
    ### GDC突变去重（由于多个肿瘤样本数据随机组合的缘故，体细胞突变会存在重复），在LUAD项目中，去重前后比例为3:1
    GDC_SNP_info = GDC_SNP_info.drop_duplicates(keep="first")
    print(f"整理后，GDC项目中的SNP突变数目共计{len(GDC_SNP_info)}，对象名称为GDC_SNP_info")

    return GDC_SNP_info

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

class exon_RNA_analysis_newer(object):

    # 初始化类实例
    # all_info为输入的所有突变数据集的dataframe
    # DNA_info为对应癌症项目的DNA突变数据集的dataframe
    # RNA_edit_info为REDIportal数据库中RNA编辑位点数据集的dataframe
    # 输出：实例自带的all_info_TP和all_info_TN
    # 输出2：最终确定用以构建三类标签的RNA突变数据集和GDC突变数据集，分别为all_info_final和GDC_SNP_info_final
    def __init__(self, all_info, DNA_info, RNA_edit_info):
        self.all_info = all_info
        self.DNA_info = DNA_info[["Chromosome", "Start_Position", "Tumor_Allele2", "Tumor_Allele1", "Tumor_Sample_UUID"]] # 仅选取染色体、碱基和case_id信息来进行合并
        self.RNA_edit_info = RNA_edit_info

        # 开始对所有RNA突变进行筛选，筛选掉基本无法处理&没必要处理的multiallelic位点——self.all_info
        print(f"所有RNA突变在进行multi-allelic筛选前，对应的突变数目为{len(self.all_info)}")
        self.all_info = self.all_info[self.all_info['record_filter'].apply(lambda x: False if x.__contains__("multiallelic") else True)]
        print(f"所有RNA突变在经过multi-allelic筛选后，剩余的突变数目为{len(self.all_info)}")

        print("="*100)

        # 进一步对RNA编辑位点数据库进行分析，判断其是否位于数据库内
        ## 找出位于RNA编辑数据库中的相应RNA突变位点——保存于self.all_info_RNA_edit中
        self.all_info_RNA_edit = pd.merge(self.all_info, self.RNA_edit_info)
        self.all_info = self.all_info.append(self.all_info_RNA_edit)
        self.all_info = self.all_info.drop_duplicates(keep=False)

        print(f"位于RNA编辑位点上的突变数为{len(self.all_info_RNA_edit)}（对象名称为all_info_RNA_edit），"
              f"不属于RNA编辑位点的突变数为{len(self.all_info)}")

        print("="*100)

        # 进一步排除位于免疫球蛋白基因内突变
        self.all_info_immunoglobulin = self.all_info[self.all_info.apply(lambda x: True if x['Hugo_Symbol'].startswith(("IGH", 'IGK', 'IGL')) else False, axis=1)]
        self.all_info = self.all_info[self.all_info.apply(lambda x: False if x['Hugo_Symbol'].startswith(("IGH", 'IGK', 'IGL')) else True, axis=1)]

        print(f"位于免疫球蛋白基因内的突变数为{len(self.all_info_immunoglobulin)}（对象名称为all_info_immunoglobulin），"
              f"不位于免疫球蛋白基因内的突变数为{len(self.all_info)}")

        print("="*100)

        # 进一步排除HLA基因内突变
        self.all_info_HLA = self.all_info[self.all_info.apply(lambda x: True if x['Hugo_Symbol'].startswith("HLA") else False, axis=1)]
        self.all_info = self.all_info[self.all_info.apply(lambda x: False if x['Hugo_Symbol'].startswith("HLA") else True, axis=1)]

        print(f"位于HLA基因内的突变数为{len(self.all_info_HLA)}（对象名称为all_info_HLA），"
              f"不位于HLA基因内的突变数为{len(self.all_info)}")

        print("="*100)

        # 开始预备TP_info, TN_info获取工作
        print(f"DNA突变数目共计{len(self.DNA_info)}")

        print("="*100)

        print("以GDC突变数据集为依据，检查RNA突变训练数据集中缺少GDC注释信息的case情况")
        # 获取对应的case_id信息
        DNA_case_info_set = set(self.DNA_info['Tumor_Sample_UUID'].value_counts().index)
        RNA_case_info_set = set(self.all_info['Tumor_Sample_UUID'].value_counts().index)
        # 找到RNA中缺少GDC数据支持的case
        RNA_only_case_info_set = RNA_case_info_set - DNA_case_info_set
        print(f"RNA中缺少GDC数据支持的case信息为{RNA_only_case_info_set}，需要对其进行重点关注")

        print(f"GDC数据集在根据RNA突变数据集信息取子集前，包含的case数目为{len(self.DNA_info['Tumor_Sample_UUID'].value_counts().index)}")
        self.DNA_info = self.DNA_info[self.DNA_info['Tumor_Sample_UUID'].isin(RNA_case_info_set)]
        print(f"GDC数据集在根据RNA突变数据集信息取子集后，包含的case数目为{len(self.DNA_info['Tumor_Sample_UUID'].value_counts().index)}")
        self.all_info_final = self.all_info.copy(deep=True)
        self.DNA_info_final = self.DNA_info.copy(deep=True)
        print("至此，纳入训练的RNA突变数据集和GDC突变数据集已确定，分别命名为all_info_final和DNA_info_final")

        print("="*100)

        print(f"开始将all_info分为TP_info, TN_info两类，初始数目为{len(self.all_info)}")
        # 获取符合GDC数据集内case_specific突变信息的RNA体细胞突变——self.TP_info
        self.TP_info = pd.merge(self.all_info, self.DNA_info_final, on=list(self.DNA_info_final.columns))
        self.TN_info = self.all_info.append(self.TP_info)
        self.TN_info = self.TN_info.drop_duplicates(keep=False)
        print(f"TP_info类获取完成，数目为{len(self.TP_info)}（对象名称为TP_info）")
        print(f"TN_info类获取完成，数目为{len(self.TN_info)}（对象名称为TN_info）")

    # 检查训练数据TP、TN相关信息
    # 输入：实例自带的TP_info和TN_info
    # 输出：打印TP、TN突变相关组成信息
    def data_check(self):
        ## 检查筛选后的TP、TN基本信息
        ### 注意，样本集中突变总数为1615931——其中位于exon区域内的TP、TN数目分别为：34187、98205
        ### 2021.1.5更新——注意，样本集中突变总数为1615931——其中位于coding区域内的TP、TN数目分别为：34435、101942
        print("TP、TN的突变数分别为：%d、%d" % (len(self.TP_info), len(self.TN_info)))
        print("TP、TN的比例数为：1:%d \n" % (len(self.TN_info) / len(self.TP_info)))

        print("TP（%d个突变）的组成信息为：" % (len(self.TP_info)))
        print(self.TP_info.Variant_Classification.value_counts())
        print("其中，Splice_Site的类型组成情况为：位点总数%d" % (len(self.TP_info.loc[self.TP_info.Variant_Classification=='Splice_Site', ])))
        print(self.TP_info.loc[self.TP_info.Variant_Classification=='Splice_Site', ].Gencode_28_secondaryVariantClassification.value_counts())

        print("\nTN（%d个突变）的组成信息为：" % (len(self.TN_info)))
        print(self.TN_info.Variant_Classification.value_counts())
        print("其中，Splice_Site的类型组成情况为：位点总数%d" % (len(self.TN_info.loc[self.TN_info.Variant_Classification=='Splice_Site', ])))
        print(self.TN_info.loc[self.TN_info.Variant_Classification=='Splice_Site', ].Gencode_28_secondaryVariantClassification.value_counts())

        print("\nTP、TN列名为：")
        print(self.TP_info.columns)

    # 对TP和TN进行合并，并进行特征处理
    # 输入：实例自带的TP_info和TN_info
    # 输出：实例新建的属性all_info、training_data和y，分别代表TP、TN合并后的所有信息、可训练的信息和类别信息
    # 输出2：实例新建的属性enc，为碱基改变情况转为one-hot变量的处理模型
    def data_preprocess(self):
        ## 2020.12.7——尝试将碱基改变情况转为one-hot变量，并纳入模型中进行判别
        ### 获取碱基改变情况信息
        TP_info_allele_change = self.TP_info['Reference_Allele'] + '>' + self.TP_info['Tumor_Allele1']
        #### 2021.5.16——修改突变碱基改变信息
        TP_info_allele_change = pd.Series([ALLELE_CHANGE_DICT[allele_change] if ALLELE_CHANGE_DICT.keys().__contains__(allele_change) else allele_change for allele_change in TP_info_allele_change])

        self.enc = preprocessing.OneHotEncoder()  # 尝试将碱基改变情况转为one-hot变量，并纳入模型中进行判别
        self.enc.fit(TP_info_allele_change.values.reshape(-1, 1))  # 根据TP中碱基改变情况分布，构建one-hot变量

        ## 训练数据合并
        ### 分别为TP和TN添加属性列
        self.TN_info['Attribute'] = "TN"
        self.TP_info['Attribute'] = "TP"
        ### 将TP、TN合并（注意此时已经将index重置）
        self.all_info = self.TP_info.append(self.TN_info, ignore_index=True)
        ### 2020.11.30 考虑添加碱基改变信息后，再构建模型判断效果改善与否
        all_info_allele_change = self.all_info['Reference_Allele'] + '>' + self.all_info['Tumor_Allele1']
        #### 2021.5.16——修改突变碱基改变信息
        all_info_allele_change = pd.Series([ALLELE_CHANGE_DICT[allele_change] if ALLELE_CHANGE_DICT.keys().__contains__(allele_change) else allele_change for allele_change in all_info_allele_change])

        all_info_allele_change_df = pd.DataFrame(self.enc.transform(all_info_allele_change.values.reshape(-1, 1)).toarray(), columns=self.enc.categories_[0])
        self.all_info = pd.concat([self.all_info, all_info_allele_change_df], axis=1)

        ## 数据进行初步处理
        ### 将training_data与Attribute列分开
        ### 将类别值转化数值，便于后面损失函数的计算
        self.all_info['Attribute'] = self.all_info['Attribute'].map({'TP': 1, 'TN': 0})
        ### 强制将‘Attribute’数据类型转化为"category"类型
        self.all_info['Attribute'] = self.all_info['Attribute'].astype('category')
        # 零值补全（TPM、COSMIC_total_alterations_in_gene）
        self.all_info = self.all_info.fillna({"Expression_TPM": 0, "COSMIC_total_alterations_in_gene": 0})
        # 将类别列等其他无关列（CHROM、POS、REF、ALT、FILTER、CONTQ）从数据集删除，并将新数据集保存为training_data
        # , "Reference_Allele_x", "Reference_Allele_y"
        # 2020.12.8 我为啥把ref_AD_normal和alt_AD_normal这两个关键特征给删掉了？？我佛了，赶紧弄回来
        self.training_data = self.all_info.drop(DROP_COLUMNS, axis=1)
        ### 按列统计na值数目，并去除所有包含na的行
        print("检查各列na值数目，情况如下所示：\n")
        print(self.training_data.isna().sum())
        self.training_data = self.training_data.dropna()
        # 缺值删除——仅对需要纳入模型的列进行处理——即删除包含缺失值的记录行
        self.all_info = self.all_info.dropna(subset=self.training_data.columns)
        # 切片，得到标签y
        self.y = self.all_info['Attribute']

    # 训练前数据相关性展示
    # 输入：实例新建的属性all_info、training_data，分别代表TP、TN合并后的所有信息、可训练的信息
    # 输出：所有特征间的相关性分析图，PCA分析图等
    # 输出2：pca计算对应变量变换pca_scaler
    def data_description(self):
        # 特征相关性检验
        corr = self.all_info.corr()  # 求数据集两两特征之间的相关性
        # sns为seaborn（基于matplotlib），diverging_palette为生成调色板对象，赋给cmap对象
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        # 创建分散颜色：h_neg, h_pos ：起始/终止颜色值
        # s ---> 值区间0-100，饱和度
        # l ---> 值区间0-100，亮度
        # n ---> 颜色个数
        # center ---> 中心颜色为浅色还是深色'light', 'dark', 默认为light
        f, ax = plt.subplots(figsize=(21, 19))  # 创建画布，定义画布大小
        sns.heatmap(corr, cmap=cmap, center=0, annot=True,
                    square=True, linewidths=.5, cbar_kws={"shrink": .5});  # 画出corr数据的热图：
        # cmap:从数字到色彩空间的映射，取值是matplotlib包里的colormap名称或颜色对象，或者表示颜色的列表；改参数默认值：根据center参数设定
        # center:数据表取值有差异时，设置热力图的色彩中心对齐值；annot(annotate的缩写):默认取值False；如果是True，在热力图每个方格写入数据；
        # square:设置热力图矩阵小块形状，默认值是False
        # linewidths:定义热力图里“表示两两特征关系的矩阵小块”之间的间隔大小
        # linecolor:切分热力图上每个矩阵小块的线的颜色，默认值是’white’
        # cbar_kws:热力图侧边绘制颜色刻度条时，相关字体设置，默认值是None

        # PCA分析——不知道是否有必要进行
        from sklearn.preprocessing import StandardScaler  # 导入数据预处理模块StandardScaler
        self.pca_scaler = StandardScaler()
        X = self.training_data
        X_scaled = self.pca_scaler.fit_transform(X)  # 将数据进行归一化处理

        self.pca = decomposition.PCA(n_components=2)  # 定义pca降维函数,保证降维后的数据只有二维
        X_pca_scaled = self.pca.fit_transform(X_scaled)  # 使用PCA方法对数据进行降维处理

        print('Projecting %d-dimensional data to 2D' % X_scaled.shape[
            1])  # 输出说明数据的维度转化：“Projecting 30-dimensional data to 2D”

        plt.figure(figsize=(12, 10))  # 定义画图板的宽和高
        plt.scatter(X_pca_scaled[:, 0], X_pca_scaled[:, 1], c=self.all_info['Attribute'].apply(lambda item: "red" if item else "blue"), alpha=1, s=1);  # 画出散点图
        # plt.colorbar()  # 画出颜色柱,先决条件是c=df['diagnosis']
        plt.title('PCA projection')  # 定义画板的名称

    # 输入：实例新建的属性training_data和y，分别代表TP、TN合并后可训练的信息和类别信息
    # 输出：实例新建的属性X_train、X_holdout、y_train、y_holdout，分别代表训练集、测试集的相关特征、实际分类（9:1进行分割）
    def data_prepare(self):
        ## 数据进行最终处理
        ### 进行训练数据标准化与训练集、测试集分割

        ### 训练集、测试集分割
        self.X_train, self.X_holdout, self.y_train, self.y_holdout = train_test_split(self.training_data, self.y, test_size=0.1, random_state=17)    #划分原数据：分为训练数据，测试数据，训练集标签和测试集标签

        ### 检查其特征是否服从正态分布
        print("如下所示为特征所对应的分布情况：")
        self.X_train.hist(figsize=(20, 15), color='c')

    # 注意：本函数仅针对常规随机森林模型进行构建，对于平衡随机森林模型，需要其他函数来处理
    # 输入：实例新建的属性X_train、X_holdout，分别代表经过标准化后的训练集、测试集的相关特征
    # 输出：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    def common_RF_build(self):
        # 考虑到为二次测试，直接进行参数选择，而且经杨哥建议，直接增加树的颗树，效果会更好，故目前仅对树的颗树进行网格搜索
        # Stratified split for the validation process
        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=17)  # 设置划分数据集的参数——10 fold
        # 字典内均为参数，需要进行网格搜索
        # initialize the set of parameters for exhaustive search and fit to find out the optimal parameters
        rfc_params = {'n_estimators': [i * 100 for i in range(13, 14)],
                      'criterion': ['gini'],
                      'max_depth': range(30, 31),
                      'min_samples_split': [2],
                      'min_samples_leaf': range(1, 2),
                      'max_features': [i / 10 for i in range(4, 5)],
                      'class_weight': ["balanced"]}  # 已确定随机森林参数组合
        # rfc_params = {'n_estimators': [i * 100 for i in range(12, 15)],
        #               'criterion': ['gini'],
        #               'max_depth': range(30, 31),
        #               'min_samples_split': [2],
        #               'min_samples_leaf': range(1, 2),
        #               'max_features': [i / 10 for i in range(4, 5)],
        #               'class_weight': ["balanced"]}  # 定义随机森林参数组合
        # n_jobs=-1设置线程数
        rfc = RandomForestClassifier(random_state=17, n_jobs=-1, oob_score=True)  # 定义随机森林算法
        # cv表示交叉验证的参数，可设置为10（10-fold）或者直接传入交叉验证对象
        self.rf_gcv = GridSearchCV(rfc, rfc_params, n_jobs=-1, cv=skf, scoring='f1', verbose=1)  # 定义网络搜索算法，优化目标指定为f1
        # 开始进行网格搜索
        print("开始进行网格搜索......")
        print("目标参数为：")
        print(rfc_params)
        self.rf_gcv.fit(self.X_train, self.y_train)  # 使用不同的参数组合训练拟合训练集
        # 网格搜索结束
        print("\n网格搜索结束......")
        print("最优参数、训练数据得分（F1）、袋外得分(oob)分别为：")
        print((self.rf_gcv.best_params_, self.rf_gcv.best_score_, self.rf_gcv.best_estimator_.oob_score_))

        # 模型初步评估
        # 对不平衡测试集的评估
        rf_pred = self.rf_gcv.predict(self.X_train)  # 预测训练集的类别
        print("\n不调整阈值，对训练集的预测得分如下所示：")
        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_train, rf_pred))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ", recall_score(self.y_train, rf_pred))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ", precision_score(self.y_train, rf_pred))  # 打印精度

        # 模型进一步评估
        # 对不平衡测试集的评估
        rf_pred = self.rf_gcv.predict(self.X_holdout)  # 预测测试集的类别
        print("\n不调整阈值，对测试集的预测得分如下所示：")
        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_pred))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ", recall_score(self.y_holdout, rf_pred))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ", precision_score(self.y_holdout, rf_pred))  # 打印精度

    # 输入：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：特征选择(feature selection)相关指标信息 + 经过特征选择后的rf_gcv_select
    def common_RF_feature_selection(self):

        # rfc_params = {'n_estimators': [i * 100 for i in range(13, 14)],
        #               'criterion': ['gini'],
        #               'max_depth': range(30, 31),
        #               'min_samples_split': [2],
        #               'min_samples_leaf': range(1, 2),
        #               'max_features': [i / 10 for i in range(4, 5)],
        #               'class_weight': ["balanced"]}  # 已确定随机森林参数组合
        # rfc = RandomForestClassifier(random_state=17, n_jobs=-1, oob_score=True)  # 定义随机森林算法

        print(f"特征筛选前，训练数据中特征数为{len(self.X_train.columns)}")
        print(f"具体特征名为{list(self.X_train.columns)}")
        print("开始进行特征筛选......")
        self.rf_gcv_select = RFECV(self.rf_gcv.best_estimator_, step=1, cv=StratifiedKFold(3),
                                   scoring="f1", verbose=1, n_jobs=-1, min_features_to_select=1)
        self.rf_gcv_select = self.rf_gcv_select.fit(self.X_train, self.y_train)
        self.selected_columns = [self.X_train.columns[i] for i in range(len((self.rf_gcv_select.support_))) if self.rf_gcv_select.support_[i]]
        print(f"特征筛选结束，筛选后所得最优特征数为{self.rf_gcv_select.n_features_}")
        print(f"具体特征名为{self.selected_columns}")
        print(f"变量筛选后的F1值为{self.rf_gcv_select.score(self.X_train, self.y_train)}")

        print("="*100)

        # 模型进一步评估
        # 对不平衡测试集的评估
        rf_pred = self.rf_gcv_select.predict(self.X_holdout)  # 预测测试集的类别
        print("\n不调整阈值，经变量筛选后所得模型对测试集的预测得分如下所示：")
        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_pred))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ", recall_score(self.y_holdout, rf_pred))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ", precision_score(self.y_holdout, rf_pred))  # 打印精度

        # Plot number of features VS. cross-validation scores
        plt.figure()
        plt.xlabel("Number of features selected")
        plt.ylabel("Cross validation score (F1 score)")
        plt.plot(range(1, len(self.rf_gcv_select.grid_scores_) + 1), self.rf_gcv_select.grid_scores_)
        plt.show()

    # 输入：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：实例新建的属性common_RF_ROC_thre、common_RF_ROC_thre_2，经由ROC曲线不同标准（tpr-fpr最大、距离(0,1)最远）评估后所获得的阈值
    def common_RF_tuning_ROC(self):
        # 应用ROC曲线进行优化
        print("\n应用ROC曲线完成对一般随机森林的优化")
        # 首先检验测试集概率分布情况
        rf_pred = self.rf_gcv.predict_proba(self.X_holdout)  # 以概率形式预测测试集的类别
        plt.hist(rf_pred[:, 1])
        # 绘制模型ROC曲线，得到ROC曲线所对应的最优界值（使tpr-fpr最大的阈值，存在不合理的可能，因为存在样本不平衡的情况）
        fpr, tpr, thresholds = roc_curve(self.y_holdout, rf_pred[:, 1])
        skplt.metrics.plot_roc(self.y_holdout, rf_pred, title="") # figsize=(6,6)可对图形尺寸进行调整
        print("模型对应ROC曲线下面积为", roc_auc_score(self.y_holdout, rf_pred[:, 1]))
        optimal_idx = np.argmax(tpr - fpr)
        self.common_RF_ROC_thre = thresholds[optimal_idx]
        print("ROC曲线上使tpr-fpr最大的阈值为：", self.common_RF_ROC_thre)

        # 找出距离(0,1)最近的阈值点（第二种思路）
        distance = pow((fpr - 0), 2) + pow((tpr - 1), 2)
        distance_sqrt = [math.sqrt(single_distance) for single_distance in distance]
        optimal_distance_idx = np.argmin(distance_sqrt)
        self.common_RF_ROC_thre_2 = thresholds[optimal_distance_idx]
        print("ROC曲线上距离(0,1)最近的阈值为：", self.common_RF_ROC_thre_2)

        # 根据对应的不同阈值来评估相应效果
        print("\n应用ROC曲线阈值一来完成测试集评估：")
        rf_predicted = (rf_pred[:, 1] >= self.common_RF_ROC_thre).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, rf_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, rf_predicted))  # 打印精度

        # 根据对应的不同阈值来评估相应效果
        print("\n应用ROC曲线阈值二来完成测试集评估：")
        rf_predicted = (rf_pred[:, 1] >= self.common_RF_ROC_thre_2).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, rf_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, rf_predicted))  # 打印精度

    # 输入：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：实例新建的属性common_RF_PR_thre、common_RF_PR_thre_2，经由PR曲线不同标准（f1最大、距离(0,1)最远）评估后所获得的阈值
    def common_RF_tuning_PR(self):
        # 应用PR曲线进行优化
        print("\n应用PR曲线完成对一般随机森林的优化")
        # 首先检验测试集概率分布情况
        rf_pred = self.rf_gcv.predict_proba(self.X_holdout)  # 以概率形式预测测试集的类别
        plt.hist(rf_pred[:, 1])
        # 绘制模型PR曲线，得到ROC曲线所对应的最优界值（使f1最大的阈值，相对比较合理）
        skplt.metrics.plot_precision_recall_curve(self.y_holdout, rf_pred, title="") # figsize=(6,6)可对图形尺寸进行调整
        precisions, recalls, thresholds = precision_recall_curve(self.y_holdout, rf_pred[:, 1])
        print("模型对应PR曲线下面积为", auc(recalls, precisions))
        optimal_idx = np.argmax((2 * precisions * recalls) / (precisions + recalls))
        self.common_RF_PR_thre = thresholds[optimal_idx]
        print("PR曲线上使F1最大的阈值为：", self.common_RF_PR_thre)

        # 找出距离(0,1)最近的阈值点（第二种思路）
        distance = pow((precisions - 1), 2) + pow((recalls - 1.0), 2)
        distance_sqrt = [math.sqrt(single_distance) for single_distance in distance]
        optimal_distance_idx = np.argmin(distance_sqrt)
        self.common_RF_PR_thre_2 = thresholds[optimal_distance_idx]
        print("PR曲线上距离(0,1)最近的阈值为：", self.common_RF_PR_thre_2)

        # 根据对应的不同阈值来评估相应效果
        print("\n应用PR曲线阈值一来完成测试集评估：")
        rf_predicted = (rf_pred[:, 1] >= self.common_RF_PR_thre).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, rf_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, rf_predicted))  # 打印精度

        # 根据对应的不同阈值来评估相应效果
        print("\n应用PR曲线阈值二来完成测试集评估：")
        rf_predicted = (rf_pred[:, 1] >= self.common_RF_PR_thre_2).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, rf_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, rf_predicted))  # 打印精度

    # 输入：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：打印特征权重对应图像
    def common_RF_assessment(self):
        # Create Data frame of Regression coefficients
        feature_importances = pd.DataFrame(self.rf_gcv.best_estimator_.feature_importances_)  # 创建新的数据表，存储回归系数
        # Merge Regression coefficients with feature names
        df_columns = pd.DataFrame(self.training_data.columns)  # 抽取特征名
        importance_and_feat = pd.merge(feature_importances, df_columns, left_index=True, right_index=True,
                                       how="left")  # 将特征和回归系数对应整合
        importance_and_feat.columns = ["feature_importance", "features"]  # 命名列名
        importance_and_feat = importance_and_feat.sort_values(by="feature_importance",
                                                              ascending=False)  # 根据特征重要性对数据进行排序

        # Set up the matplotlib figure
        plt.rcParams['figure.figsize'] = (10, 8)  # 定义画布大小
        # Let's draw top 10 important features
        sns.barplot(x='features', y='feature_importance', data=importance_and_feat).set_title('Feature importance plot')  # 画出直方图
        plt.xticks(rotation=90)  # 将x坐标上对应名称旋转90度

    # 输入：实例新建的属性rf_gcv、scaler，代表经过网格搜索所得到的最优RF模型、标准化工具
    # 输出：保存相应模型+工具
    def common_RF_save(self, model_folder_path):
        # 判断模型文件夹是否存在
        if not os.path.exists(model_folder_path):
            os.mkdir(model_folder_path)
        # 保存当前最佳模型
        with open(os.path.join(model_folder_path, self.__class__.__name__ + '.own.model'), 'wb') as f:
            joblib.dump(self.rf_gcv.best_estimator_, f, compress=3)
        # 保存当前one-hot-encoder工具
        with open(os.path.join(model_folder_path, self.__class__.__name__ + '.own.one_hot_encoder'), 'wb') as f:
            joblib.dump(self.enc, f)
        # 保存当前用以训练的列名信息
        pd.DataFrame(pd.Series(self.training_data.columns)).to_csv(os.path.join(model_folder_path, self.__class__.__name__ + '.own.training_data_col'), index=False)

    # 2021.07.05 根据上述最佳模型，绘制独立验证集中的PR曲线，并根据使F1最大化的阈值来输出相应positive、negative突变集
    # 应用目前效果最好的模型进行测试
    # 输入：实例新建的属性rf_gcv、scaler，代表经过网格搜索所得到的最优RF模型、标准化工具
    # 输入2：实例新建的属性independent_info,为独立测试数据集对应的dataframe
    # 输出：实例新建的属性independent_info_training、independent_info_training_scaler,为独立测试数据集对应的特征信息、标准化后的特征信息
    # 输出2：实例新建的属性independent_positive_thre，为独立测试数据集根据最大化F1值所得阈值对应的阳性位点（negative为阴性位点）
    # 输出3：实例新建的属性gdc_validate_df_SNP，为经过突变去重+SNP选取的gdc突变数据集
    # 输出4：实例新建的属性independent_pred、independent_pred_thre、independent_label，为根据相应模型和gdc突变数据集结果，所得独立验证集的预测概率值、预测标签和实际标签
    # 输出5：实例新建的属性all_info_final、GDC_SNP_info_final，最终确定用以构建标签的RNA突变数据集和GDC突变数据集
    def common_RF_utilize(self, independent_info, gdc_validate_df):

        self.independent_info = independent_info

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

        # print("="*100)
        #
        # print("位点筛选完成，进一步根据已有信息建立新特征")
        # self.independent_info = self.new_feature_constructor(self.independent_info)

        # self.independent_info_low_exp = self.independent_info.loc[self.independent_info['DP_tumor'] <= 10, ]
        # self.independent_info = self.independent_info.loc[self.independent_info['DP_tumor'] > 10,]
        # print(f"independent_info类中低表达量部分和正常表达量部分获取完成，数目分别为{len(self.independent_info_low_exp)}（对象名称为independent_info_low_exp）和{len(self.independent_info)}（对象名称为independent_info）")

        print("---------------------------------------------")
        print("对gdc数据库进行预处理，所得去重后的gdc数据库SNP位点信息为gdc_validate_df_SNP")
        # gdc数据预处理
        # 对gdc突变数据集进行预处理
        # 仅选择SNP突变信息进行分析
        self.gdc_validate_df_SNP = gdc_validate_df[gdc_validate_df['Variant_Type'] == "SNP"]
        # 重新编制index信息
        self.gdc_validate_df_SNP.reset_index(drop=True, inplace=True)
        # 添加“Tumor_Sample_UUID”列，使之与RNA体细胞突变中信息保持一致（经过Funcotator重新构建列名）
        del self.gdc_validate_df_SNP['Tumor_Sample_UUID']
        self.gdc_validate_df_SNP["Tumor_Sample_UUID"] = pd.Series(["-".join(GDC_sample_info.split("-")[0:3]) for GDC_sample_info in
             self.gdc_validate_df_SNP["Tumor_Sample_Barcode"]])
        # 仅选取染色体、碱基信息和case_id信息，便于GDC相关信息和待处理所有突变信息进行合并和分离
        # 2021.1.13 删除"Variant_Classification"列选取
        self.gdc_validate_df_SNP = self.gdc_validate_df_SNP.loc[:,["Chromosome", "Start_Position", "Reference_Allele", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_UUID"]]
        # 重新命名，便于进行信息合并(注意在Funcotator的注释信息中，Tumor_Allele1对应突变碱基、Tumor_Allele2对应参考碱基)
        self.gdc_validate_df_SNP.columns = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele2", "Tumor_Allele1", "Tumor_Sample_UUID"]
        # 去除GDC数据集中重复的突变，避免重复合并
        print("用以验证的GDC项目中原有%d个体细胞突变" % (len(self.gdc_validate_df_SNP)))
        print("开始去重.......................")
        self.gdc_validate_df_SNP = self.gdc_validate_df_SNP.drop_duplicates(keep="first")
        print("经过去重后剩余%d个突变" % (len(self.gdc_validate_df_SNP)))

        print("="*100)

        print("以GDC突变数据集为依据，检查RNA突变训练数据集中缺少GDC注释信息的case情况")
        # 获取对应的case_id信息
        GDC_case_info_set = set(self.gdc_validate_df_SNP['Tumor_Sample_UUID'].value_counts().index)
        RNA_case_info_set = set(self.independent_info['Tumor_Sample_UUID'].value_counts().index)
        # 找到RNA中缺少GDC数据支持的case
        RNA_only_case_info_set = RNA_case_info_set - GDC_case_info_set
        print(f"RNA中case总数为{len(RNA_case_info_set)}，其中缺少GDC数据支持的case信息为{RNA_only_case_info_set}，需要对其进行重点关注")

        print(f"GDC数据集在根据RNA突变数据集信息取子集前，包含的case数目为{len(self.gdc_validate_df_SNP['Tumor_Sample_UUID'].value_counts().index)}")
        self.gdc_validate_df_SNP = self.gdc_validate_df_SNP[self.gdc_validate_df_SNP['Tumor_Sample_UUID'].isin(RNA_case_info_set)]
        print(f"GDC数据集在根据RNA突变数据集信息取子集后，包含的case数目为{len(self.gdc_validate_df_SNP['Tumor_Sample_UUID'].value_counts().index)}")

        self.independent_info_final = self.independent_info.copy(deep=True)
        self.gdc_validate_df_SNP_final = self.gdc_validate_df_SNP.copy(deep=True)
        print("至此，RNA突变数据集和GDC突变数据集已确定，分别命名为independent_info_final和gdc_validate_df_SNP_final")

        print("=" * 100)

        # 检查independent_info和GDC数据间的重叠/差异情况
        self.independent_info_TP = pd.merge(self.independent_info, self.gdc_validate_df_SNP, on=list(self.gdc_validate_df_SNP.columns))
        print(f"independent_info_TP类获取完成，数目为{len(self.independent_info_TP)}（对象名称为independent_info_TP）")
        self.independent_info = self.independent_info.append(self.independent_info_TP)
        self.independent_info_TN = self.independent_info.drop_duplicates(keep=False)
        print(f"independent_info_TN类获取完成，数目为{len(self.independent_info_TN)}（对象名称为independent_info_TN）")
        # 分别为TP和TN添加标签列
        self.independent_info_TN['DNA_label'] = 0
        self.independent_info_TP['DNA_label'] = 1
        # 将TP、TN合并（注意此时已经将index重置）
        self.independent_info = self.independent_info_TP.append(self.independent_info_TN, ignore_index=True)

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
        self.independent_info_training = self.independent_info[self.training_data.columns]
        self.independent_info_training = self.independent_info_training.dropna()
        # 缺值删除——仅对需要纳入模型的列进行处理——删除包含缺失值的记录，以保持后续切片过程中，行数的一致
        self.independent_info = self.independent_info.dropna(subset=self.independent_info_training.columns)

        # 切片，得到标签independent_label
        self.independent_label = self.independent_info['DNA_label']

        # 数据概览
        # 检验测试集概率分布情况
        print("\n独立测试数据集概率分布情况如下所示：")
        self.independent_pred = self.rf_gcv.predict_proba(self.independent_info_training)  # 概率形式预测测试集的类别
        plt.hist(self.independent_pred[:, 1])

        print("---------------------------------------------")
        print("开始获取independent_info所对应的gdc突变注释信息independent_label")
        # 数据验证
        # 应用gdc数据来完成P-R曲线绘制+最优参数选取工作
        # 绘制模型PR曲线，得到ROC曲线所对应的最优界值（使f1最大的阈值，相对比较合理）
        skplt.metrics.plot_precision_recall_curve(self.independent_label, self.independent_pred)
        precisions, recalls, thresholds = precision_recall_curve(self.independent_label, self.independent_pred[:, 1])
        print("模型对应PR曲线下面积为", auc(recalls, precisions))

        print("---------------------------------------------")
        print("应用模型对独立验证集进行预测")
        print("independent_info行数为%d" % (len(self.independent_info)))
        # 数据预测
        self.independent_pred_thre = self.rf_gcv.predict(self.independent_info_training)  # 概率形式预测测试集的类别
        print(f"预测完成，所得预测结果independent_pred_thre行数为{len(self.independent_pred_thre)}")
        # 获取预测阳性、阴性的结果，报告并导出
        self.independent_positive_thre = self.independent_info[self.independent_pred_thre == 1]
        self.independent_negative_thre = self.independent_info[self.independent_pred_thre == 0]
        print("总SNP数目为%d，其中通过模型判别为阳性的突变位点independent_positive_thre数目为%d，占比为1:%d" % (
        len(self.independent_info), len(self.independent_positive_thre), len(self.independent_info) / len(self.independent_positive_thre)))

        positive_RNA_case_info_set = set(self.independent_positive_thre['Tumor_Sample_UUID'].value_counts().index)
        negative_RNA_case_info_set = set(self.independent_negative_thre['Tumor_Sample_UUID'].value_counts().index)
        print(f"阴性突变位点independent_negative_thre所对应的case数为{len(negative_RNA_case_info_set)}，阳性突变位点所对应的case数为{len(positive_RNA_case_info_set)}，两者差集为{negative_RNA_case_info_set - positive_RNA_case_info_set}，需要重点关注")

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.independent_label, self.independent_pred_thre))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.independent_label, self.independent_pred_thre))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.independent_label, self.independent_pred_thre))  # 打印精度

    # 暂时搁置（并非很好的阳性结果）
    # 应用PCA检查相应independent_info的分布情况——观察独立验证集特征分布是否与训练集（LUAD）一致
    def common_RF_independent_check(self, gdc_validate_df_SNP):
        # PCA分析
        # X = self.independent_info_training
        # X_scaled = self.pca_scaler.transform(X)  # 使用训练集中对应的标准化规则将数据进行归一化处理
        #
        # X_pca_scaled = self.pca.transform(X_scaled)  # 使用训练集中对应的PCA方法对数据进行降维处理
        #
        # print('Projecting %d-dimensional data to 2D' % X_scaled.shape[1])  # 输出说明数据的维度转化：“Projecting 30-dimensional data to 2D”
        #
        # plt.figure(figsize=(12, 10))  # 定义画图板的宽和高
        # plt.scatter(X_pca_scaled[:, 0], X_pca_scaled[:, 1], alpha=0.7, s=10);  # 画出散点图
        # plt.colorbar()  # 画出颜色柱,先决条件是c=df['diagnosis']
        # plt.title('PCA projection')  # 定义画板的名称

        # 新PCA分析（引入第三方验证系统）
        X = self.independent_info
        X_TP = pd.merge(X, gdc_validate_df_SNP)
        X_TN = X.append(X_TP)
        X_TN = X_TN.drop_duplicates(keep=False)
        X_TN['Attribute'] = 0
        X_TP['Attribute'] = 1
        X = X_TP.append(X_TN, ignore_index=True)
        X['Attribute'] = X['Attribute'].astype('category')
        X_training = X[self.training_data.columns]

        X_scaled = self.pca_scaler.transform(X_training)  # 使用训练集中对应的标准化规则将数据进行归一化处理

        X_pca_scaled = self.pca.transform(X_scaled)  # 使用训练集中对应的PCA方法对数据进行降维处理

        print('Projecting %d-dimensional data to 2D' % X_scaled.shape[1])  # 输出说明数据的维度转化：“Projecting 30-dimensional data to 2D”

        plt.figure(figsize=(12, 10))  # 定义画图板的宽和高
        plt.scatter(X_pca_scaled[:, 0], X_pca_scaled[:, 1], c=X['Attribute'].apply(lambda item: "red" if item else "blue"), alpha=1, s=1);  # 画出散点图
        # plt.colorbar()  # 画出颜色柱,先决条件是c=df['diagnosis']
        plt.title('PCA projection')  # 定义画板的名称

    # 对单个case（单个突变位点）的判别理由进行说明
    # 输入：实例对应训练模型rf_gcv，代表对应的训练模型
    # 输入2：实例新建属性independent_info，代表输入的所有独立验证集信息
    def common_RF_interpret(self):
        explainer = lime.lime_tabular.LimeTabularExplainer(self.X_train, feature_names=self.X_train.columns, class_names=self.y_train, discretize_continuous=True)

    # 工具函数1：对训练数据进行处理，基于本身数据建立新特征
    def new_feature_constructor(self, data_info):
        # 2021.08.17
        # 建立Read Orientation Artifact相关的新特征（F1R2或F2R1的最大占比）并加以检查
        # 工具函数
        def F1R2_F2R1_tumor_alt_ratio_construct(x):
            # F1R2和F2R1之和为0时，比例直接为1（经过随机筛选验证所得）
            if (x['F1R2_tumor'] + x['F2R1_tumor']==0):
                return 1
            else:
                if x['F1R2_tumor'] > x['F2R1_tumor']:
                    return(x['F1R2_tumor'] / (x['F1R2_tumor'] + x['F2R1_tumor']))
                else:
                    return(x['F2R1_tumor'] / (x['F1R2_tumor'] + x['F2R1_tumor']))

        # 正式运行
        print("根据F1R2、F2R1信息建立Read Orientation Artifact相关的新特征——F1R2_F2R1_tumor_alt_ratio")
        data_info['F1R2_F2R1_tumor_alt_ratio'] = data_info.apply(F1R2_F2R1_tumor_alt_ratio_construct, axis=1)

        return data_info

if __name__ == '__main__':
    # 获取模型训练所需数据文件
    all_info = pd.read_table(args.raw_RNA_mutations)
    DNA_info = pd.read_table(args.DNA_mutations)

    ## 附注：仅对GDC需要做此处理
    DNA_info = GDC_site_info_retrieve(DNA_info)

    RNA_EDIT_INFO = RNA_EDIT_process(args.REDIportal, args.DARNED)

    # 初始化所有训练数据
    common_RF_exon = exon_RNA_analysis_newer(all_info, DNA_info, RNA_EDIT_INFO)

    # 检查TP、TN分布情况
    common_RF_exon.data_check()

    # 检查是否有NaN值存在干扰模型构建
    common_RF_exon.data_preprocess()

    # 进行数据分割——train_test_split
    common_RF_exon.data_prepare()

    # 构建加权随机森林模型WRF
    common_RF_exon.common_RF_build()

    # 获取对应P-R信息
    common_RF_exon.common_RF_tuning_PR()

    # 保存所得模型
    common_RF_exon.common_RF_save(args.model_folder_path)