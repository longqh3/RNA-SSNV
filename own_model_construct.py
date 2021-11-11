# 2021.10.8 #
# test successful
# python /home/lqh/Codes/Python/RNA-SSNV/own_model_construct.py \
# --REDIportal /home/lqh/resources/database/RNA_edit/REDIportal/REDIportal_main_table.hg38.bed \
# --DARNED /home/lqh/resources/database/RNA_edit/DARNED_hg19_to_bed_to_hg38_rm_alt.bed \
# --raw_RNA_mutations /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/LUAD/RNA/RNA_somatic_mutation/VcfAssembly_new/SNP_WES_Interval_exon.txt \
# --DNA_mutations /home/lqh/Codes/Data/TCGA_maf_files/TCGA-LUAD \
# --model_folder_path /home/lqh/Codes/Python/RNA-SSNV/model

# Basic packages
import os
import math
import numpy as np    #导入Python科学计算的基础软件包numpy
import pandas as pd     #导入python的一个数据分析包pandas
# Data analysis packages
## Models
from sklearn.ensemble import RandomForestClassifier  # 导入随机森林算法
## Data pre-process
from sklearn import preprocessing  # 导入数据预处理包
from sklearn.model_selection import train_test_split  # 导入训练集和测试集划分函数tain_test_split
## Data training
from sklearn.model_selection import StratifiedKFold     #导入将数据划分函数StratifiedKFold
from sklearn.model_selection import GridSearchCV    #导入网格搜索自动调参函数GridSearchCV
## Model assessing
from sklearn.metrics import *    #导入metrics模块的所有函数，metrics模块实现了一些函数，用来评估预测误差。已使用：precision_recall_curve
# Visualize
import matplotlib.pyplot as plt    #导入Python可视化Matplotlib模块的绘图pyplot函数
# Others
import joblib # 其他方式保存模型
import warnings    #导入Python中的warnings模块
warnings.filterwarnings('ignore')


import argparse


parser=argparse.ArgumentParser(description="Own model construction for RNA-SSNV.")
# demo option
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.')

parser.add_argument('--REDIportal', help='hg38 REDIportal bed file.')
parser.add_argument('--DARNED', help='hg38 DARNED bed file.')
parser.add_argument('--raw_RNA_mutations', '-r', help='raw RNA somatic single nucleotide variants.')
parser.add_argument("--DNA_mutations", help="GDC mutations.")
parser.add_argument("--model_folder_path", help="Folder path for constructed model.")
parser.add_argument("--num_threads", type=int, help="Number of threads.")

args=parser.parse_args()

# Specify variables to be dropped before building own training dataset
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

# Hard-code rules for base changes
ALLELE_CHANGE_DICT = {
    "T>C":"A>G",
    "C>T":"G>A",
    "T>G":"A>C",
    "G>T":"C>A",
    "T>A":"A>T",
    "G>C":"C>G"
}


# For GDC maf file specifically
# De-duplicate and extract sites info from all GDC SNP records which included columns:
# Tumor_Sample_UUID Chromosome Start_Position Reference_Allele Tumor_Allele1 Tumor_Allele2
def GDC_site_info_retrieve(GDC_info):
    print("Mutational status for GDC project: \n", GDC_info["Variant_Type"].value_counts())
    # select SNP
    GDC_SNP_info = GDC_info[GDC_info['Variant_Type'] == "SNP"]
    GDC_SNP_info.reset_index(drop=True, inplace=True)
    # reformat "Tumor_Sample_UUID" column to accommodate case ID within RNA mutations
    del GDC_SNP_info['Tumor_Sample_UUID']
    GDC_SNP_info["Tumor_Sample_UUID"] = pd.Series(["-".join(GDC_sample_info.split("-")[0:3])
                                                        for GDC_sample_info in
                                                        GDC_SNP_info["Tumor_Sample_Barcode"]])
    # select specified columns, rename columns and de-duplicate
    GDC_SNP_info = GDC_SNP_info.loc[:,
                        ["Chromosome", "Start_Position", "Reference_Allele", "Reference_Allele", "Tumor_Seq_Allele2",
                         "Tumor_Sample_UUID"]]
    GDC_SNP_info.columns = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele2", "Tumor_Allele1", "Tumor_Sample_UUID"]
    GDC_SNP_info = GDC_SNP_info.drop_duplicates(keep="first")
    print(f"After sorting, counts for SNPs within GDC project: {len(GDC_SNP_info)}, Attribute: GDC_SNP_info")

    return GDC_SNP_info

def RNA_EDIT_process(REDIportal, DARNED):
    """Fetches RNA editing sites info from REDIportal and DARNED databases.

    Args:
        REDIportal: Location of the bed file downloaded from REDIportal database.
        DARNED: Location of the bed file downloaded from DARNED database.

    Returns:
        Merged sites info extracted from REDIportal and DARNED databases.
        Only two columns retained: "Chromosome", "Start_Position"
    """
    # Read in sites info as DataFrames
    REDIportal_info = pd.read_table(REDIportal, header=None)
    DARNED_info = pd.read_table(DARNED, header=None)
    # Pre-process database info to extract required sites info
    REDIportal_info = REDIportal_info[[0, 2]]
    DARNED_info = DARNED_info[[0, 2]]
    print(f"REDIportal and DARNED database contained {len(REDIportal_info)}和{len(DARNED_info)} RNA editing sites respectively")
    # Merge sites info  from two databases
    RNA_EDIT_INFO = pd.concat([REDIportal_info, DARNED_info], ignore_index=True)
    RNA_EDIT_INFO.columns = ["Chromosome", "Start_Position"]
    # De-duplicate sites
    RNA_EDIT_INFO = RNA_EDIT_INFO.drop_duplicates(keep="first")
    print(f"After de-duplication，we got {len(RNA_EDIT_INFO)} RNA editing sites totally.")

    return RNA_EDIT_INFO

class exon_RNA_analysis_newer(object):
    """Training + testing datasets preparation,
    model construction and assessment,
    model save.

    Attributes:
        all_info: The DataFrame containing all features and annotations from "own_data_vcf_info_retriver.py".
        DNA_info: The DataFrame containing sites, allele change and case info from self-called DNA mutations (GDC project maf file).
        RNA_edit_info: The DataFrame containing RNA editing sites info from REDIportal and DARNED databases.
        num_threads: Number of threads available for model training.
    """

    def __init__(self, all_info, DNA_info, RNA_edit_info, num_threads):
        """Inits exon_RNA_analysis with two new attributes categories:
        TP&TN-all_info_TP, all_info_TP
        final utilized datasets-all_info_final, GDC_SNP_info_final"""

        self.all_info = all_info
        self.DNA_info = DNA_info[["Chromosome", "Start_Position", "Tumor_Allele2", "Tumor_Allele1", "Tumor_Sample_UUID"]] # 仅选取染色体、碱基和case_id信息来进行合并
        self.RNA_edit_info = RNA_edit_info
        self.num_threads = num_threads

        # Primary filter: remove multi-allelic sites hard to deal with
        # update "all_info" attribute
        print(f"Before multi-allelic filtering, counts for all RNA somatic mutations was {len(self.all_info)}")
        self.all_info = self.all_info[self.all_info['record_filter'].apply(lambda x: False if x.__contains__("multiallelic") else True)]
        print(f"After multi-allelic filtering, counts for all RNA somatic mutations was {len(self.all_info)}")

        print("="*100)

        # Second filter: remove RNA editing sites
        # update attribute: "all_info"
        # add attribute: "all_info_RNA_edit"
        self.all_info_RNA_edit = pd.merge(self.all_info, self.RNA_edit_info)
        self.all_info = self.all_info.append(self.all_info_RNA_edit)
        self.all_info = self.all_info.drop_duplicates(keep=False)

        print(f"Counts for mutations located within RNA editing sites: {len(self.all_info_RNA_edit)}（attribute: all_info_RNA_edit），"
              f"Counts for mutations not located within RNA editing sites: {len(self.all_info)}")

        print("="*100)

        # Third filter: remove mutations within immunoglobin genes
        # update attribute: "all_info"
        # add attribute: "all_info_immunoglobulin"
        self.all_info_immunoglobulin = self.all_info[self.all_info.apply(lambda x: True if x['Hugo_Symbol'].startswith(("IGH", 'IGK', 'IGL')) else False, axis=1)]
        self.all_info = self.all_info[self.all_info.apply(lambda x: False if x['Hugo_Symbol'].startswith(("IGH", 'IGK', 'IGL')) else True, axis=1)]

        print(f"Counts for mutations located within immunoglobin genes: {len(self.all_info_immunoglobulin)}（attribute: all_info_immunoglobulin），"
              f"Counts for mutations not located within immunoglobin genes: {len(self.all_info)}")

        print("="*100)

        # Forth filter: remove mutations within HLA genes
        # update attribute: "all_info"
        # add attribute: "all_info_HLA"
        self.all_info_HLA = self.all_info[self.all_info.apply(lambda x: True if x['Hugo_Symbol'].startswith("HLA") else False, axis=1)]
        self.all_info = self.all_info[self.all_info.apply(lambda x: False if x['Hugo_Symbol'].startswith("HLA") else True, axis=1)]

        print(f"Counts for mutations located within HLA genes: {len(self.all_info_HLA)}（attribute: all_info_HLA），"
              f"Counts for mutations not located within HLA genes: {len(self.all_info)}")

        print("="*100)

        # Determine cases used in model training

        # We use overlap cases of DNA and RNA mutation sets as training cases.
        # We also pay special attention to cases lacking DNA mutation evidence.
        # add attribute: "all_info_final", "DNA_info_final"

        print(f"DNA mutations count: {len(self.DNA_info)}")

        print("="*100)

        print("Use DNA mutation set as reference, check the case status within RNA mutation set.")
        DNA_case_info_set = set(self.DNA_info['Tumor_Sample_UUID'].value_counts().index)
        RNA_case_info_set = set(self.all_info['Tumor_Sample_UUID'].value_counts().index)
        RNA_only_case_info_set = RNA_case_info_set - DNA_case_info_set
        print(f"The RNA cases lacking DNA support were {RNA_only_case_info_set}, special attention was required.")

        print(f"Before subsetting, DNA mutation set contained {len(self.DNA_info['Tumor_Sample_UUID'].value_counts().index)} cases")
        self.DNA_info = self.DNA_info[self.DNA_info['Tumor_Sample_UUID'].isin(RNA_case_info_set)]
        print(f"After subsetting, DNA mutation set contained {len(self.DNA_info['Tumor_Sample_UUID'].value_counts().index)} cases")
        self.all_info_final = self.all_info.copy(deep=True)
        self.DNA_info_final = self.DNA_info.copy(deep=True)
        print(f"Finally, we determined final RNA and GDC mutation sets: all_info_final and DNA_info_final, "
              f"which contained {len(self.all_info_final)} and {len(self.DNA_info_final)} mutations respectively")

        print("="*100)

        # Split RNA mutations into two categories

        # We use DNA mutations as references to split RNA mutations.
        # Finally, we constructed training data with two categories (TP_info and TN_info)
        # add attribute: "TP_info" and "TN_info"

        print(f"Start to split RNA mutations into two categories, initial mutation count was {len(self.all_info)}")
        # First, RNA somatic mutations concordant with DNA evidence were retrieved——self.TP_info
        self.TP_info = pd.merge(self.all_info, self.DNA_info_final, on=list(self.DNA_info_final.columns))
        self.TN_info = self.all_info.append(self.TP_info)
        self.TN_info = self.TN_info.drop_duplicates(keep=False)
        print(f"TP_info category: {len(self.TP_info)}（attribute: TP_info）")
        print(f"TN_info category: {len(self.TN_info)}（attribute: TN_info）")

    def data_check(self):
        """Check TP, TN basic composition
        """

        print("TP、TN mutation counts：%d、%d" % (len(self.TP_info), len(self.TN_info)))
        print("TP:TN proportion：1:%d \n" % (len(self.TN_info) / len(self.TP_info)))

        print(f"TP ({len(self.TP_info)}) mutations: ")  # TP
        print(self.TP_info.Variant_Classification.value_counts())
        print("Splice_Site mutations：%d" % (len(self.TP_info.loc[self.TP_info.Variant_Classification=='Splice_Site', ])))
        print(self.TP_info.loc[self.TP_info.Variant_Classification=='Splice_Site', ].Gencode_28_secondaryVariantClassification.value_counts())

        print("\nTN (%d mutations) mutations：" % (len(self.TN_info)))  # TN
        print(self.TN_info.Variant_Classification.value_counts())
        print("Splice_Site mutations：%d" % (len(self.TN_info.loc[self.TN_info.Variant_Classification=='Splice_Site', ])))
        print(self.TN_info.loc[self.TN_info.Variant_Classification=='Splice_Site', ].Gencode_28_secondaryVariantClassification.value_counts())

        print("\nColumn names for TP、TN categories: ")
        print(self.TP_info.columns)

    def data_preprocess(self):
        """Build new features and combine TP, TN info

        Transfer allele changes into one-hot variables and included into model training.
        Construct standard training dataset and one-hot encoder.

        """

        # change allele change status into one-hot variables accordingly

        TP_info_allele_change = self.TP_info['Reference_Allele'] + '>' + self.TP_info['Tumor_Allele1']
        TP_info_allele_change = pd.Series([ALLELE_CHANGE_DICT[allele_change] if ALLELE_CHANGE_DICT.keys().__contains__(allele_change) else allele_change for allele_change in TP_info_allele_change])

        self.enc = preprocessing.OneHotEncoder()
        self.enc.fit(TP_info_allele_change.values.reshape(-1, 1))  # train OneHotEncoder using TP's distribution

        # construct training dataset by combining TP&TN and build one-hot variables

        self.TN_info['Attribute'] = "TN"
        self.TP_info['Attribute'] = "TP"
        self.all_info = self.TP_info.append(self.TN_info, ignore_index=True)
        all_info_allele_change = self.all_info['Reference_Allele'] + '>' + self.all_info['Tumor_Allele1']
        all_info_allele_change = pd.Series([ALLELE_CHANGE_DICT[allele_change] if ALLELE_CHANGE_DICT.keys().__contains__(allele_change) else allele_change for allele_change in all_info_allele_change])
        all_info_allele_change_df = pd.DataFrame(self.enc.transform(all_info_allele_change.values.reshape(-1, 1)).toarray(), columns=self.enc.categories_[0])
        self.all_info = pd.concat([self.all_info, all_info_allele_change_df], axis=1)

        # split X and y from training dataset
        # add attribute: "training_data", "y"

        self.all_info['Attribute'] = self.all_info['Attribute'].map({'TP': 1, 'TN': 0})
        self.all_info['Attribute'] = self.all_info['Attribute'].astype('category')
        # self.all_info = self.all_info.fillna({"Expression_TPM": 0, "COSMIC_total_alterations_in_gene": 0})  # fill na
        self.training_data = self.all_info.drop(DROP_COLUMNS, axis=1)
        print("Check NA count for each column, result listed below: \n")
        print(self.training_data.isna().sum())
        self.training_data = self.training_data.dropna()
        self.all_info = self.all_info.dropna(subset=self.training_data.columns)   # drop na
        self.y = self.all_info['Attribute']  # get y

    def data_prepare(self):
        """Prepare for final training and testing data (9:1 split).
        Check for alll features' distribution using histogram.
        add attribute: "X_train", "X_holdout", "y_train", "y_holdout"
        """

        self.X_train, self.X_holdout, self.y_train, self.y_holdout = train_test_split(self.training_data, self.y, test_size=0.1, random_state=17)    #划分原数据：分为训练数据，测试数据，训练集标签和测试集标签

        # self.X_train.hist(figsize=(20, 15), color='c')

    def common_RF_build(self):
        """Build a weighted random forest model tuned by Grid-Search through 10x cross-validation.
        Primary model assessment using testing dataset
        add attribute: "rf_gcv"
        """
        # model training

        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=17)  # set up 10x cross-validation
        rfc_params = {'n_estimators': [i * 100 for i in range(13, 14)],
                      'criterion': ['gini'],
                      'max_depth': range(30, 31),
                      'min_samples_split': [2],
                      'min_samples_leaf': range(1, 2),
                      'max_features': [i / 10 for i in range(4, 5)],
                      'class_weight': ["balanced"]}  # parameter combination determined by Grid-Search
        rfc = RandomForestClassifier(random_state=17, n_jobs=self.num_threads, oob_score=True)  # init weighted random forest model
        self.rf_gcv = GridSearchCV(rfc, rfc_params, n_jobs=-1, cv=skf, scoring='f1', verbose=1)  # init Grid-Search
        print("Start Grid-Searching......")
        print("Target parameter matrix was: ")
        print(rfc_params)
        self.rf_gcv.fit(self.X_train, self.y_train)  # fit model using training data by Grid-Searching
        print("\nGrid-Search finished......")
        print("Best parameters, cross-validated F1 score and OOB(out of bag) score were: ")
        print((self.rf_gcv.best_params_, self.rf_gcv.best_score_, self.rf_gcv.best_estimator_.oob_score_))

        # model assessment within training and testing dataset

        rf_pred = self.rf_gcv.predict(self.X_train)  # training dataset
        print("\nModel performance within training dataset: ")
        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_train, rf_pred))  # accuracy
        print("Recall Score (how much of TP were predicted correctly) : ", recall_score(self.y_train, rf_pred))  # recall
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ", precision_score(self.y_train, rf_pred))  # precision

        rf_pred = self.rf_gcv.predict(self.X_holdout)  # testing dataset
        print("\nModel performance within testing dataset: ")
        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_pred))  # accuracy
        print("Recall Score (how much of TP were predicted correctly) : ", recall_score(self.y_holdout, rf_pred))  # recall
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ", precision_score(self.y_holdout, rf_pred))  # precision

    def common_RF_tuning_PR(self):
        """Use P-R curve to assess fitted model.
        Find two optimal thresholds and assess their corresponding performances.
        add attribute: "common_RF_PR_thre", "common_RF_PR_thre_2"
        """

        print("\nUsing P-R curve to assess trained model's performance within testing dataset...")
        rf_pred = self.rf_gcv.predict_proba(self.X_holdout)  # predict probabilities
        plt.hist(rf_pred[:, 1])  # plot prob distribution
        precisions, recalls, thresholds = precision_recall_curve(self.y_holdout, rf_pred[:, 1])  # P-R values
        print("AUC for PR curve: ", auc(recalls, precisions))
        optimal_idx = np.argmax((2 * precisions * recalls) / (precisions + recalls))
        self.common_RF_PR_thre = thresholds[optimal_idx]
        print("Optimal threshold 1 for maximizing F1 score: ", self.common_RF_PR_thre)

        distance = pow((precisions - 1), 2) + pow((recalls - 1.0), 2)
        distance_sqrt = [math.sqrt(single_distance) for single_distance in distance]
        optimal_distance_idx = np.argmin(distance_sqrt)
        self.common_RF_PR_thre_2 = thresholds[optimal_distance_idx]
        print("Optimal threshold 2 for the shortest distance with (0,1) in PR curve：", self.common_RF_PR_thre_2)

        print("\nUsing threshold 1 to assess performance：")
        rf_predicted = (rf_pred[:, 1] >= self.common_RF_PR_thre).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_predicted))  # accuracy
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, rf_predicted))  # recall
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, rf_predicted))  # precision

        print("\nUsing threshold 2 to assess performance：")
        rf_predicted = (rf_pred[:, 1] >= self.common_RF_PR_thre_2).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_predicted))  # accuracy
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, rf_predicted))  # recall
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, rf_predicted))  # precision

    def common_RF_save(self, model_folder_path):
        """Save trained model, one-hot encoder and features used in training.

        :param model_folder_path: Folder path to store saved model.
        """

        if not os.path.exists(model_folder_path):
            os.mkdir(model_folder_path)  # new folder
        with open(os.path.join(model_folder_path, self.__class__.__name__ + '.own.model'), 'wb') as f:
            joblib.dump(self.rf_gcv.best_estimator_, f, compress=3)  # save model
        with open(os.path.join(model_folder_path, self.__class__.__name__ + '.own.one_hot_encoder'), 'wb') as f:
            joblib.dump(self.enc, f)  # save encoder
        pd.DataFrame(pd.Series(self.training_data.columns)).to_csv(os.path.join(model_folder_path, self.__class__.__name__ + '.own.training_data_col'), index=False)  # save features

if __name__ == '__main__':
    # read in all required files
    all_info = pd.read_table(args.raw_RNA_mutations)
    DNA_info = pd.read_table(args.DNA_mutations)
    num_threads = args.num_threads

    # Be advised, only GDC maf file required this process
    DNA_info = GDC_site_info_retrieve(DNA_info)

    RNA_EDIT_INFO = RNA_EDIT_process(args.REDIportal, args.DARNED)

    # init model training process
    common_RF_exon = exon_RNA_analysis_newer(all_info, DNA_info, RNA_EDIT_INFO, num_threads)

    # check TP, TN distribution
    common_RF_exon.data_check()

    # check na value and prepare data
    common_RF_exon.data_preprocess()

    # split data into training and testing datasets
    common_RF_exon.data_prepare()

    # train weighted random forest model
    common_RF_exon.common_RF_build()

    # retrieve P-R info
    common_RF_exon.common_RF_tuning_PR()

    # save trained model
    common_RF_exon.common_RF_save(args.model_folder_path)