# 2021.10.8 #
# Finished testing#
# python /home/lqh/Codes/Python/RNA-SSNV/model_construct.py \
# --REDIportal /home/lqh/resources/database/RNA_edit/REDIportal/REDIportal_main_table.hg38.bed \
# --DARNED /home/lqh/resources/database/RNA_edit/DARNED_hg19_to_bed_to_hg38_rm_alt.bed \
# --raw_RNA_mutations /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/LUAD/RNA/RNA_somatic_mutation/VcfAssembly_new/SNP_WES_Interval_exon.txt \
# --GDC_mutations /home/lqh/Codes/Data/TCGA_maf_files/TCGA-LUAD \
# --WES_mutations /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/LUAD/WXS/WXS_somatic_mutation/VariantsToTable/PASS_SNP_WES_Interval_exon.table \
# --model_folder_path /home/lqh/Codes/Python/RNA-SSNV/model \
# --num_threads 20

# Basic packages
import os
import math
import numpy as np
import pandas as pd
# Data analysis packages
## Models
from sklearn.feature_selection import RFECV
from sklearn import decomposition
from sklearn.ensemble import RandomForestClassifier
## Data pre-process
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
## Data training
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV
## Model assessing
from sklearn.metrics import *
# Visualize
import matplotlib.pyplot as plt
import seaborn as sns
# Others
import joblib
import warnings
# Model explanation
import lime
import lime.lime_tabular

import argparse
warnings.filterwarnings('ignore')

parser=argparse.ArgumentParser(description="Model construction for RNA-SSNV.")
# demo option
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.')
parser.add_argument('--REDIportal', help='hg38 REDIportal bed file.')
parser.add_argument('--DARNED', help='hg38 DARNED bed file.')
parser.add_argument('--raw_RNA_mutations', '-r', help='raw RNA somatic single nucleotide variants.')
parser.add_argument("--GDC_mutations", help="GDC mutations.")
parser.add_argument("--WES_mutations", help="WES mutations.")
parser.add_argument("--model_folder_path", help="Folder path for constructed model.")
parser.add_argument("--num_threads", type=int, help="Number of threads.")

args=parser.parse_args()

# Specify variables to be dropped before building training dataset
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

class exon_RNA_analysis(object):
    """Training + testing datasets preparation,
    model construction and assessment,
    model save.

    Attributes:
        all_info: The DataFrame containing all features and annotations from "own_data_vcf_info_retriver.py".
        GDC_info: The DataFrame containing info from TCGA project maf file downloaded from "gdc-maf-tool.py".
        WES_info: The DataFrame containing sites, allele change and case info from self-called WES somatic mutations.
        RNA_edit_info: The DataFrame containing RNA editing sites info from REDIportal and DARNED databases.
        num_threads: Number of threads available for model training.
    """

    def __init__(self, all_info, GDC_info, WES_info, RNA_edit_info, num_threads):
        """Inits exon_RNA_analysis with two new attributes categories:
        TP&TN-all_info_TP, all_info_TP
        final utilized datasets-all_info_final, GDC_SNP_info_final"""

        self.all_info = all_info
        self.GDC_info = GDC_info
        self.WES_info = WES_info
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

        # Prepare to construct training dataset (three classes): Counts for mutations not located within HLA genes:

        # Sort GDC related info
        # We select SNP sites info from GDC mutations (columns renamed) and
        # re-format "Tumor_Sample_UUID" column to help conduct de-duplication
        # update attribute: "GDC_SNP_info"

        print("Mutation count for current GDC cancer project was: \n", self.GDC_info["Variant_Type"].value_counts())
        self.GDC_SNP_info = self.GDC_info[self.GDC_info['Variant_Type'] == "SNP"] # SNP
        self.GDC_SNP_info.reset_index(drop=True, inplace=True)
        del self.GDC_SNP_info['Tumor_Sample_UUID'] # re-format Tumor_Sample_UUID" column
        self.GDC_SNP_info["Tumor_Sample_UUID"] = pd.Series(["-".join(GDC_sample_info.split("-")[0:3])
                                                            for GDC_sample_info in self.GDC_SNP_info["Tumor_Sample_Barcode"]])
        self.GDC_SNP_info = self.GDC_SNP_info.loc[:, ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_UUID"]]
        self.GDC_SNP_info.columns = ["Chromosome", "Start_Position", "Tumor_Allele2", "Tumor_Allele1", "Tumor_Sample_UUID"] # rename columns
        self.GDC_SNP_info = self.GDC_SNP_info.drop_duplicates(keep="first") # de-duplication
        print(f"After sorting, mutation count reduced into {len(self.GDC_SNP_info)}")

        # Combine WES mutations with GDC

        # Through overlapping GDC mutations with WES somatic mutations,
        # we introduce two parts (GDC_SNP_WES_info, WES_info_GDC_trimmed) to help
        # construct TP, TN and Ambiguous parts
        # add attribute: "GDC_SNP_WES_info", "WES_info_GDC_trimmed"

        print(f"For WES, Mutect2 called {len(self.WES_info)} mutations")
        print("Next, we combined WES mutations with GDC")
        self.GDC_SNP_WES_info = pd.merge(self.GDC_SNP_info, self.WES_info)  # overlapping
        print(f"Their intersection contained {len(self.GDC_SNP_WES_info)} mutations (attribute: GDC_SNP_WES_info)")
        WES_info_temp = self.WES_info.append(self.GDC_SNP_WES_info)
        self.WES_info_GDC_trimmed = WES_info_temp.drop_duplicates(keep=False)  # WES only
        print(f"Their difference set contained {len(self.WES_info_GDC_trimmed)} mutations (attribute: WES_info_GDC_trimmed)")

        print("="*100)

        # Determine cases used in model training

        # We use overlap cases of GDC and RNA mutation sets as training cases.
        # We also pay special attention to cases lacking GDC mutation evidence.
        # add attribute: "all_info_final", "GDC_SNP_info_final"

        print("Use GDC mutation set as reference, check the case status within RNA mutation set.")
        GDC_case_info_set = set(self.GDC_SNP_info['Tumor_Sample_UUID'].value_counts().index)
        RNA_case_info_set = set(self.all_info['Tumor_Sample_UUID'].value_counts().index)  # case sets
        RNA_only_case_info_set = RNA_case_info_set - GDC_case_info_set  # RNA only cases
        print(f"The RNA cases lacking GDC support were {RNA_only_case_info_set}, special attention was required.")

        print(f"Before subsetting, GDC mutation set contained {len(self.GDC_SNP_info['Tumor_Sample_UUID'].value_counts().index)} cases.")
        self.GDC_SNP_info = self.GDC_SNP_info[self.GDC_SNP_info['Tumor_Sample_UUID'].isin(RNA_case_info_set)]
        print(f"After subsetting, GDC mutation set contained {len(self.GDC_SNP_info['Tumor_Sample_UUID'].value_counts().index)} cases")
        self.all_info_final = self.all_info.copy(deep=True)
        self.GDC_SNP_info_final = self.GDC_SNP_info.copy(deep=True)  # final preparation training data construction
        print(f"Finally, we determined final RNA and GDC mutation sets: all_info_final and GDC_SNP_info_final, "
              f"which contained {len(self.all_info_final)} and {len(self.GDC_SNP_info_final)} mutations respectively")

        print("="*100)

        # Split RNA mutations into three categories

        # We use GDC and WES mutations as references to split RNA mutations.
        # Finally, we constructed training data with three categories (TP_info, TN_info and Ambiguity_info)
        # add attribute: "TP_info", "Ambiguity_info_Mutect2" and "TN_info"

        print(f"Start to split RNA mutations into three categories, initial mutation count was {len(self.all_info)}")
        # First, RNA somatic mutations concordant with GDC evidence were retrieved——self.TP_info
        self.TP_info = pd.merge(self.all_info, self.GDC_SNP_info, on=list(self.GDC_SNP_info.columns))
        self.all_info = self.all_info.append(self.TP_info)
        self.all_info = self.all_info.drop_duplicates(keep=False)
        print(f"TP_info category: {len(self.TP_info)}（attribute: TP_info）")

        # Subsequently, RNA somatic mutations concordant with Mutect2 evidence were retrieved——self.Ambiguity_info_Mutect2
        self.Ambiguity_info_Mutect2 = pd.merge(self.all_info, self.WES_info_GDC_trimmed, on=list(self.WES_info_GDC_trimmed.columns))
        self.TN_info = self.all_info.append(self.Ambiguity_info_Mutect2)
        self.TN_info = self.TN_info.drop_duplicates(keep=False)
        print(f"Ambiguity_info category: {len(self.Ambiguity_info_Mutect2)}（attribute: Ambiguity_info_Mutect2）")
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
        self.training_data = self.all_info.drop(DROP_COLUMNS, axis=1)  # drop un-necessary columns and get X
        print("Check NA count for each column, result listed below: \n")
        print(self.training_data.isna().sum())
        self.training_data = self.training_data.dropna()
        self.all_info = self.all_info.dropna(subset=self.training_data.columns)  # drop na
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
        self.rf_gcv = GridSearchCV(rfc, rfc_params, n_jobs=self.num_threads, cv=skf, scoring='f1', verbose=1)  # init Grid-Search
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
        print("Accuracy Score : (how much of variants type was predicted correctly): ",
              accuracy_score(self.y_holdout, rf_pred))  # accuracy
        print("Recall Score (how much of TP were predicted correctly): ", recall_score(self.y_holdout, rf_pred))  # recall
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
        with open(os.path.join(model_folder_path, self.__class__.__name__ + '.model'), 'wb') as f:
            joblib.dump(self.rf_gcv.best_estimator_, f, compress=3)  # save model
        with open(os.path.join(model_folder_path, self.__class__.__name__ + '.one_hot_encoder'), 'wb') as f:
            joblib.dump(self.enc, f)  # save encoder
        pd.DataFrame(pd.Series(self.training_data.columns)).to_csv(os.path.join(model_folder_path, self.__class__.__name__ + '.training_data_col'), index=False)  # save features

if __name__ == '__main__':
    # read in all required files
    all_info = pd.read_table(args.raw_RNA_mutations)
    GDC_info = pd.read_table(args.GDC_mutations)
    WES_info = pd.read_table(args.WES_mutations)
    RNA_EDIT_INFO = RNA_EDIT_process(args.REDIportal, args.DARNED)
    num_threads = args.num_threads

    # init model training process
    common_RF_exon = exon_RNA_analysis(all_info, GDC_info, WES_info, RNA_EDIT_INFO, num_threads)

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