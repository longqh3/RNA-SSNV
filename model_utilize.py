# 2021.10.9 #
# finished test
# python /home/lqh/Codes/Python/RNA-SSNV/model_utilize.py \
# --REDIportal /home/lqh/resources/database/RNA_edit/REDIportal/REDIportal_main_table.hg38.bed \
# --DARNED /home/lqh/resources/database/RNA_edit/DARNED_hg19_to_bed_to_hg38_rm_alt.bed \
# --raw_RNA_mutations /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/LUAD/RNA/RNA_somatic_mutation/VcfAssembly_new/SNP_WES_Interval_exon.txt \
# --model_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer.model \
# --one_hot_encoder_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer.one_hot_encoder \
# --training_columns_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer.training_data_col \
# --output_table_path /home/lqh/Codes/Python/RNA-SSNV/output/LUAD.table

# Basic packages
import pandas as pd     #导入python的一个数据分析包pandas
import joblib # 其他方式保存模型
import warnings    #导入Python中的warnings模块
warnings.filterwarnings('ignore')

import argparse

parser=argparse.ArgumentParser(description="Model utilization for RNA-SSNV.")
# demo option
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.')

parser.add_argument('--REDIportal', help='hg38 REDIportal bed file.')
parser.add_argument('--DARNED', help='hg38 DARNED bed file.')
parser.add_argument('--raw_RNA_mutations', '-r', help='raw RNA somatic single nucleotide variants.')
parser.add_argument("--model_path", help="Path for constructed model.")
parser.add_argument("--one_hot_encoder_path", help="Path for one-hot encoder.")
parser.add_argument("--training_columns_path", help="Path for constructed model.")
parser.add_argument("--output_table_path", help="Path for final output table.")

args=parser.parse_args()

# Assign allele change rules
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

class model_utilize_no_DNA(object):
    """Read in model and utilize it into called mutations.

    Attributes:
        independent_info: The DataFrame containing all features and annotations from "own_data_vcf_info_retriver.py".
        RNA_edit_info: The DataFrame containing RNA editing sites info from REDIportal and DARNED databases.
        model_path: Path of machine learning model.
        one_hot_encoder_path: Path of one-hot encoder.
        training_columns_path: Path of training data's features.
        output_table_path: Path of output table (added predicted probabilities and labels).
    """

    def __init__(self, independent_info, RNA_edit_info, model_path, one_hot_encoder_path, training_columns_path, output_table_path):
        """Inits model utilization with data and model"""

        # data
        self.independent_info = independent_info
        self.RNA_edit_info = RNA_edit_info
        # model
        self.rf_gcv = joblib.load(model_path)
        self.enc = joblib.load(one_hot_encoder_path)
        self.training_columns = pd.Series(pd.read_table(training_columns_path)["0"])
        self.output_table_path = output_table_path

    def common_RF_utilize(self):
        """Utilize multi-filtering strategy and model discrimination
        """


        # Primary filter: remove multi-allelic sites hard to deal with
        # update "independent_info" attribute
        print(f"Before multi-allelic filtering, counts for all RNA somatic mutations was {len(self.independent_info)}")
        self.independent_info = self.independent_info[self.independent_info['record_filter'].apply(lambda x: False if x.__contains__("multiallelic") else True)]
        print(f"After multi-allelic filtering, counts for all RNA somatic mutations was {len(self.independent_info)}")

        print("="*100)

        # Second filter: remove RNA editing sites
        # update attribute: "independent_info"
        # add attribute: "independent_info_RNA_edit"
        self.independent_info_RNA_edit = pd.merge(self.independent_info, self.RNA_edit_info)
        self.independent_info = self.independent_info.append(self.independent_info_RNA_edit)
        self.independent_info = self.independent_info.drop_duplicates(keep=False)

        print(f"Counts for mutations located within RNA editing sites: {len(self.independent_info_RNA_edit)}（attribute: independent_info_RNA_edit），"
              f"Counts for mutations not located within RNA editing sites: {len(self.independent_info)}")

        print("="*100)

        # Third filter: remove mutations within immunoglobin genes
        # update attribute: "independent_info"
        # add attribute: "independent_info_immunoglobulin"
        self.independent_info_immunoglobulin = self.independent_info[self.independent_info.apply(lambda x: True if x['Hugo_Symbol'].startswith(("IGH", 'IGK', 'IGL')) else False, axis=1)]
        self.independent_info = self.independent_info[self.independent_info.apply(lambda x: False if x['Hugo_Symbol'].startswith(("IGH", 'IGK', 'IGL')) else True, axis=1)]

        print(f"Counts for mutations located within immunoglobin genes: {len(self.independent_info_immunoglobulin)}（attribute: independent_info_immunoglobulin），"
              f"Counts for mutations not located within immunoglobin genes: {len(self.independent_info)}")

        print("="*100)

        # Forth filter: remove mutations within HLA genes
        # update attribute: "independent_info"
        # add attribute: "independent_info_HLA"
        self.independent_info_HLA = self.independent_info[self.independent_info.apply(lambda x: True if x['Hugo_Symbol'].startswith("HLA") else False, axis=1)]
        self.independent_info = self.independent_info[self.independent_info.apply(lambda x: False if x['Hugo_Symbol'].startswith("HLA") else True, axis=1)]

        print(f"Counts for mutations located within HLA genes: {len(self.independent_info_HLA)}（attribute: independent_info_HLA），"
              f"Counts for mutations not located within HLA genes: {len(self.independent_info)}")

        print("=" * 100)

        # Prepare to construct utilization dataset: Counts for mutations not located within HLA genes:

        # reset index for dataset
        self.independent_info = self.independent_info.reset_index(drop=True)

        # add allelic change info
        independent_info_allele_change = self.independent_info['Reference_Allele'] + '>' + self.independent_info['Tumor_Allele1']
        independent_info_allele_change = pd.Series([ALLELE_CHANGE_DICT[allele_change] if ALLELE_CHANGE_DICT.keys().__contains__(allele_change) else allele_change for allele_change in independent_info_allele_change])
        independent_info_allele_change_df = pd.DataFrame(self.enc.transform(independent_info_allele_change.values.reshape(-1, 1)).toarray(), columns=self.enc.categories_[0])
        self.independent_info = pd.concat([self.independent_info, independent_info_allele_change_df], axis=1)

        # retrieve only training columns
        self.independent_info_training = self.independent_info[self.training_columns]
        # drop rows containing NA
        self.independent_info_training = self.independent_info_training.dropna()
        self.independent_info = self.independent_info.dropna(subset=self.independent_info_training.columns)

        # add prediction info
        self.independent_pred = self.rf_gcv.predict_proba(self.independent_info_training)[:, 1]  # pred prob
        self.independent_info['pred_prob'] = self.independent_pred
        self.independent_pred_thre = self.rf_gcv.predict(self.independent_info_training)  # pred label
        self.independent_info['pred_label'] = self.independent_pred_thre
        print(f"Prediction finished, positive count: {self.independent_info['pred_label'].value_counts()[1]}, negative count: {self.independent_info['pred_label'].value_counts()[0]}")

        # export dataset
        self.independent_info.to_csv(self.output_table_path, sep="\t", index=False)

if __name__ == '__main__':
    independent_info = pd.read_table(args.raw_RNA_mutations)
    RNA_EDIT_INFO = RNA_EDIT_process(args.REDIportal, args.DARNED)

    model_utilize = model_utilize_no_DNA(independent_info, RNA_EDIT_INFO, args.model_path, args.one_hot_encoder_path, args.training_columns_path, args.output_table_path)
    model_utilize.common_RF_utilize()