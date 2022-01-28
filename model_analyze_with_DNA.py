# 2021.10.8 #
# 测试完成，一切正常#
# 测试命令 #
# python /home/lqh/Codes/Python/RNA-SSNV/model_analyze_with_DNA.py \
# --step 1 \
# --cancer_type LUSC \
# --DNA_info /home/lqh/Codes/Data/TCGA_maf_files/TCGA-LUSC \
# --RNA_info /home/lqh/Codes/Python/RNA-SSNV/output/LUSC.table \
# --WXS_target_interval /home/lqh/resources/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals_add_chr_to_hg38_rm_alt.bed \
# --exon_interval /home/lqh/resources/database/gencode/GRCh38_GENCODE_v22_exon_rm_alt.bed \
# --RNA_calling_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/tables/info/LUSC_RNA_somatic_calling_info.tsv \
# --RNA_bam_folder /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/LUSC/RNA/apply_BQSR \
# --Mutect2_target_detected_sites /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/LUSC/RNA/RNA_somatic_mutation/VariantsToTable/SNP_WES_Interval_exon.table \
# --project_folder /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results \
# --num_threads 64 \
# --output_file_path /home/lqh/Codes/Python/RNA-SSNV/output/test_DNA_step_1.class

# --template_vcf_file /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer_.training_data_col \

# 测试完成，一切正常#
# python /home/lqh/Codes/Python/RNA-SSNV/model_analyze_with_DNA.py \
# --step 2 \
# --force_call_RNA_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/BLCA/RNA/RNA_somatic_mutation/VcfAssembly_new/Mutect2_force_call.txt \
# --instance_path /home/lqh/Codes/Python/RNA-SSNV/output/BLCA_DNA_step_1.class \
# --model_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis.model \
# --one_hot_encoder_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis.one_hot_encoder \
# --training_columns_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis.training_data_col \
# --output_file_path /home/lqh/Codes/Python/RNA-SSNV/output/BLCA.final.table

# 导入相关需求包
# basic
import os
import numpy as np    #导入Python科学计算的基础软件包numpy
import pandas as pd     #导入python的一个数据分析包pandas
# data prepare
import portion as P
import vcf
import pysam
# multi-process
from multiprocessing import Pool
from multiprocessing import Manager
# other
import joblib # 其他方式保存模型
import warnings    #导入Python中的warnings模块
warnings.filterwarnings('ignore')

import argparse

parser=argparse.ArgumentParser(description="A discriminate model construction pipeline for RNA-SSNV.")
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.') # demo
# Generic parameter
parser.add_argument('--step', type=int, help='Step indicator of model utilization with DNA.')
parser.add_argument("--output_file_path", help="Path for temp/final output file.")
# Step 1
parser.add_argument('--cancer_type', help='Cancer type for current study.')
parser.add_argument('--DNA_info', help='DNA somatic mutations used to validate.')
parser.add_argument('--RNA_info', help='RNA somatic mutations with model predicting status.')
parser.add_argument('--WXS_target_interval', help='WXS target interval bed file.')
parser.add_argument('--exon_interval', help='GENCODE v22 exon interval bed file.')
parser.add_argument('--RNA_calling_info', help='Tabular info used in RNA somatic mutation calling.')
parser.add_argument('--RNA_bam_folder', help='RNA bam file folder path.')
parser.add_argument('--Mutect2_target_detected_sites', help='Table for Mutect2 detected sites.')
parser.add_argument('--template_vcf_file', default=os.path.join(os.path.split(os.path.realpath(__file__))[0], "resources/template.vcf"), help='Template vcf file to conduct force-calling.')
parser.add_argument('--project_folder', help='Project folder path.')
parser.add_argument('--num_threads', '-nt', type=int, help='Number of threads allowed.')
# Step 2
parser.add_argument('--force_call_RNA_info', help='Force called DNA only mutations in RNA.')
parser.add_argument('--instance_path', help='Data analysis instance path from step 1.')
parser.add_argument("--model_path", help="Path for constructed model.")
parser.add_argument("--one_hot_encoder_path", help="Path for one-hot encoder.")
parser.add_argument("--training_columns_path", help="Path for constructed model.")

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

# toolkit: 0. Custom function for "error_callback" parameter
def print_error(self, value):
    print("error: ", value)

# toolkit: 1. Extract information from samfile objects
def tumor_bam_record_retrive(df_record, samfile_list):
    """Retrieve best and other allelic coverages and mapping qualities info
       from bam objects listed in samfile_list.

    Args:
        df_record: Single record from pyvcf retrieved record list.
        samfile_list: List of pysam retrieved samfile objects.

    Returns:
        Coverage and mapping quality info
    """

    try:
        # retrieve basic allelic info
        ref = df_record.Reference_Allele
        alt = df_record.Tumor_Allele1
        chrom = df_record.Chromosome
        position = df_record.Start_Position
        # init other info
        other_coverage_list = []
        other_median_mapping_quality_list = []

        # init coverage&mapping quality info using first samfile object (minimum base quality: 20, concordant with GATK)
        initial_count_coverage_info = count_coverage_decompose(
            samfile_list[0].count_coverage(contig=chrom, start=position - 1, stop=position, quality_threshold=20))
        initial_median_mapping_quality = np.median([x.mapping_quality for x in samfile_list[0].fetch(contig=chrom, start=position - 1, stop=position)])
        optimal_ref_coverage = initial_count_coverage_info[str(ref)]
        optimal_alt_coverage = initial_count_coverage_info[str(alt)]
        optimal_median_mapping_quality = initial_median_mapping_quality
        # traverse other samfile objects and conduct judgement to update
        for samfile in samfile_list[1:]:
            # pysam apply 0-based coordinate system (minimum base quality: 20, concordant with GATK)
            count_coverage_info = count_coverage_decompose(
                samfile.count_coverage(contig=chrom, start=position - 1, stop=position, quality_threshold=20))
            current_ref_coverage = count_coverage_info[str(ref)]
            current_alt_coverage = count_coverage_info[str(alt)]
            current_median_mapping_quality = np.median([x.mapping_quality for x in samfile.fetch(contig=chrom, start=position - 1, stop=position)])
            # compare with optimal values from alt and ref allelic changes
            if current_alt_coverage > optimal_alt_coverage:
                # save previous optimal values and update
                other_coverage_list.append(str(optimal_ref_coverage) + "/" + str(optimal_alt_coverage))
                other_median_mapping_quality_list.append(str(current_median_mapping_quality))
                optimal_ref_coverage, optimal_alt_coverage, optimal_median_mapping_quality = current_ref_coverage, current_alt_coverage, current_median_mapping_quality
            elif current_alt_coverage == optimal_alt_coverage:
                if current_ref_coverage > optimal_ref_coverage:
                    # save previous optimal values and update
                    other_coverage_list.append(str(optimal_ref_coverage) + "/" + str(optimal_alt_coverage))
                    other_median_mapping_quality_list.append(str(current_median_mapping_quality))
                    optimal_ref_coverage, optimal_alt_coverage, optimal_median_mapping_quality = current_ref_coverage, current_alt_coverage, current_median_mapping_quality
                else:
                    # save current values
                    other_coverage_list.append(str(current_ref_coverage) + "/" + str(current_alt_coverage))
                    other_median_mapping_quality_list.append(str(current_median_mapping_quality))
            else:
                # save current values
                other_coverage_list.append(str(current_ref_coverage) + "/" + str(current_alt_coverage))
                other_median_mapping_quality_list.append(str(current_median_mapping_quality))

        return optimal_ref_coverage, optimal_alt_coverage, ";".join(other_coverage_list), optimal_median_mapping_quality, ";".join(other_median_mapping_quality_list)
    except Exception as ex:
        print(ex)
        print("tumor_bam_record_retrive error......")

# toolkit: 1.1 Convert returned value of "count_coverage function" within pysam library into dict
def count_coverage_decompose(count_coverage_info):
    try:
        count_coverage_dict = {}
        # returned value had ACGT order
        count_coverage_dict['A'] = count_coverage_info[0][0]
        count_coverage_dict['C'] = count_coverage_info[1][0]
        count_coverage_dict['G'] = count_coverage_info[2][0]
        count_coverage_dict['T'] = count_coverage_info[3][0]
        return count_coverage_dict
    except Exception as ex:
        print(ex)
        print("count_coverage_decompose error......")

# toolkit: 2. Add RNA bam info into mutations info
def RNA_bam_info_adder(case_id, gdc_validate_all_info_DNA_only_RNA_missing_case_list,
                       gdc_validate_all_info_DNA_only_RNA_missing_case_only, RNA_calling_info, RNA_bam_folder_loc,
                       tumor_bam_record_retrive):

    print(f"Retrieve RNA sequencing info according to mutational records within {case_id}")
    # retrieve tumor aliquots_id info from calling table
    RNA_tumor_aliquots_id = RNA_calling_info.loc[(RNA_calling_info['case_id'] == case_id) & (
            RNA_calling_info['sample_type'] == "Primary Tumor"), 'aliquots_id']
    print(f"aliquots_id for tumor RNA sequencing within {case_id} were {RNA_tumor_aliquots_id}")
    # retrieve RNA tumor bam file locations and read in as samfile objects
    RNA_tumor_case_file_paths = RNA_bam_folder_loc + "/" + RNA_tumor_aliquots_id + ".bam"
    samfile_list = [pysam.AlignmentFile(RNA_tumor_case_file_path, "rb") for RNA_tumor_case_file_path in
                    RNA_tumor_case_file_paths]

    # init info
    gdc_validate_all_info_DNA_only_RNA_missing_case_only['ref_AD_tumor_bam_RNA'] = ""
    gdc_validate_all_info_DNA_only_RNA_missing_case_only['alt_AD_tumor_bam_RNA'] = ""
    gdc_validate_all_info_DNA_only_RNA_missing_case_only['other_AD_tumor_bam_RNA'] = ""
    gdc_validate_all_info_DNA_only_RNA_missing_case_only['median_MQ_tumor_bam_RNA'] = ""
    gdc_validate_all_info_DNA_only_RNA_missing_case_only['other_median_MQ_tumo_bam_RNA'] = ""
    for i in gdc_validate_all_info_DNA_only_RNA_missing_case_only.index:
        # retrieve corresponding info and store
        gdc_validate_all_info_DNA_only_RNA_missing_case_only.loc[i, 'ref_AD_tumor_bam_RNA'], \
        gdc_validate_all_info_DNA_only_RNA_missing_case_only.loc[i, 'alt_AD_tumor_bam_RNA'], \
        gdc_validate_all_info_DNA_only_RNA_missing_case_only.loc[i, 'other_AD_tumor_bam_RNA'], \
        gdc_validate_all_info_DNA_only_RNA_missing_case_only.loc[i, 'median_MQ_tumor_bam_RNA'], \
        gdc_validate_all_info_DNA_only_RNA_missing_case_only.loc[i, 'other_median_MQ_tumo_bam_RNA'] = \
            tumor_bam_record_retrive(gdc_validate_all_info_DNA_only_RNA_missing_case_only.loc[i,], samfile_list)

    [samfile.close() for samfile in samfile_list]
    # add info
    gdc_validate_all_info_DNA_only_RNA_missing_case_list.append(gdc_validate_all_info_DNA_only_RNA_missing_case_only)

class data_prepare():
    """Provide function codes for analysis."""

    def GDC_site_info_retrieve(self, GDC_info):
        """De-duplicate and extract sites info from all GDC SNP records which included columns:
           Tumor_Sample_UUID Chromosome Start_Position Reference_Allele Tumor_Allele1 Tumor_Allele2
        """
        
        print("Mutational status for GDC project: \n", GDC_info["Variant_Type"].value_counts())
        # select SNP
        self.GDC_SNP_info = GDC_info[GDC_info['Variant_Type'] == "SNP"]
        self.GDC_SNP_info.reset_index(drop=True, inplace=True)
        # reformat "Tumor_Sample_UUID" column to accommodate case ID within RNA mutations
        del self.GDC_SNP_info['Tumor_Sample_UUID']
        self.GDC_SNP_info["Tumor_Sample_UUID"] = pd.Series(["-".join(GDC_sample_info.split("-")[0:3])
                                                       for GDC_sample_info in
                                                       self.GDC_SNP_info["Tumor_Sample_Barcode"]])
        # select specified columns, rename columns and de-duplicate
        self.GDC_SNP_info = self.GDC_SNP_info.loc[:,
                       ["Chromosome", "Start_Position", "Reference_Allele", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_UUID"]]
        self.GDC_SNP_info.columns = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele2", "Tumor_Allele1", "Tumor_Sample_UUID"]
        self.GDC_SNP_info = self.GDC_SNP_info.drop_duplicates(keep="first")
        print(f"After sorting, counts for SNPs within GDC project: {len(self.GDC_SNP_info)}, Attribute: GDC_SNP_info")

        return self.GDC_SNP_info

    def WXS_exon_region_dict_generate(self, WXS_target_interval_path, exon_interval_path, num_threads):
        """Fetch WXS target regions and exon regions (bed foramt)
        """

        print("Start to fetch interval info for WXS target regions...")
        self.WXS_target_interval = dict(self.IntervalSet_generate_control(
            pd.read_table(WXS_target_interval_path, names=["chr", "start", "end", "target_info", "strand"]), num_threads))
        print("Fetch complete, Attribute: WXS_target_interval")

        print("="*100)

        print("Start to fetch interval info for exon regions...")
        self.exon_interval = dict(self.IntervalSet_generate_control(
            pd.read_table(exon_interval_path, names=["chr", "start", "end", "exon_info", "count", "strand"]), num_threads))
        print("Fetch complete, Attribute: exon_interval")

    # auxiliary function
    def IntervalSet_generate(self, region_df):
        """Create corresponding IntervalSet for each chromosome.

        Generate dict containing IntervalSet for each chromosome

        Args:
            region_df: a dataframe containing chr, start, end and other info

        Returns:
            A dict containing Interval objects whose key was chromosome name.
        """

        region_interval_dict = {}
        for single_chr in set(region_df["chr"]):
            single_chr_region_df = region_df.loc[region_df["chr"] == single_chr,]
            single_chr_interval_set = P.empty()
            # generate Interval object for each region
            for i in single_chr_region_df.index:
                single_chr_interval_set = single_chr_interval_set | P.closed(single_chr_region_df['start'][i] + 1,
                                                                             single_chr_region_df['end'][i])
                # avoid possible '-' strand start>end issue
                if single_chr_region_df['start'][i] > single_chr_region_df['end'][i]:
                    print(
                        f"'{single_chr_region_df['exon_info'][i]}'s end({single_chr_region_df['end'][i]}) was less than start({single_chr_region_df['start'][i]}) position, need to switch order.")
            # generate IntervalSet object for this chromosome
            region_interval_dict[single_chr] = single_chr_interval_set
        return region_interval_dict

    # auxiliary function
    def IntervalSet_generate_subprocess(self, region_df, single_chr, region_interval_dict):
        """Subprocess for multi-process function
        """

        print(f"Start to process interval info located within {single_chr}......")
        single_chr_region_df = region_df.loc[region_df["chr"] == single_chr,]
        single_chr_interval_set = P.empty()
        # generate Interval object for each region
        for i in single_chr_region_df.index:
            single_chr_interval_set = single_chr_interval_set | P.closed(single_chr_region_df['start'][i] + 1,
                                                                         single_chr_region_df['end'][i])
            # avoid possible '-' strand start>end issue
            if single_chr_region_df['start'][i] > single_chr_region_df['end'][i]:
                print(
                    f"'{single_chr_region_df['exon_info'][i]}'s end({single_chr_region_df['end'][i]}) was less than start({single_chr_region_df['start'][i]}) position, need to switch order.")
        # generate IntervalSet object for this chromosome
        region_interval_dict[single_chr] = single_chr_interval_set
        print(f"Finished processing interval info located within {single_chr}.")

    # auxiliary function
    def IntervalSet_generate_control(self, region_df, num_threads):
        """Multi-process mode to generate interval set
        """

        print('Parent process %s.' % os.getpid())
        p = Pool(num_threads)
        manager = Manager()
        region_interval_dict = manager.dict()
        for single_chr in set(region_df["chr"]):
            p.apply_async(self.IntervalSet_generate_subprocess, args=(region_df, single_chr, region_interval_dict,))
        print('Waiting for all subprocesses done...')
        p.close()
        p.join()
        print('All subprocesses done.')
        return (region_interval_dict)

    def driver_gene_info_extract(self, driver_gene_loc, cancer_type):
        """Retrieve specified cancer driver gene set from IntoGen database
        """

        driver_genes = pd.read_table(driver_gene_loc)
        print(f"Fetch {cancer_type}'s corresponding driver genes.")
        self.cancer_spec_driver_genes = set(driver_genes[driver_genes.CANCER_TYPE==cancer_type]['SYMBOL'])
        print(f"Totally, we got {len(self.cancer_spec_driver_genes)} driver genes for {cancer_type}, Attribute: cancer_spec_driver_genes")

    def model_predict_interpret_prepare(self, RF_model_loc, one_hot_encoder_loc, training_col_loc):
        """Fetch all info correlated with prediction model, one-hot encoder and training columns
        """

        print("Start to read all info correlated with prediction model, one-hot encoder and training columns...")
        f1, f2 = open(one_hot_encoder_loc, 'rb'), open(RF_model_loc, 'rb')
        self.enc, self.rf_gcv = joblib.load(f1), joblib.load(f2)
        f1.close()
        f2.close()
        self.training_col = pd.read_table(training_col_loc)["0"]

        # print("开始构建模型解释相关信息...")
        # self.shap_explainer = shap.TreeExplainer(self.rf_gcv)

        print("Read in complete.")

    def model_predict(self, data_info):
        """Predict prob and labels for given data

        :return: predict prob, labels
        """

        # Assign allele change rules
        ALLELE_CHANGE_DICT = {
            "T>C": "A>G",
            "C>T": "G>A",
            "T>G": "A>C",
            "G>T": "C>A",
            "T>A": "A>T",
            "G>C": "C>G"
        }

        self.data_info = data_info.reset_index(drop=True)

        # add allelic change info
        data_info_allele_change = self.data_info['Reference_Allele'] + '>' + self.data_info['Tumor_Allele1']
        data_info_allele_change = pd.Series([ALLELE_CHANGE_DICT[allele_change] if ALLELE_CHANGE_DICT.keys().__contains__(allele_change) else allele_change for allele_change in data_info_allele_change])
        data_info_allele_change_df = pd.DataFrame(self.enc.transform(data_info_allele_change.values.reshape(-1, 1)).toarray(), columns=self.enc.categories_[0])
        # avoid redundant columns
        if set(self.enc.categories_[0]).issubset(set(self.data_info.columns)):
            self.data_info.drop(list(self.enc.categories_[0]), axis=1, inplace=True)
        self.data_info = pd.concat([self.data_info, data_info_allele_change_df], axis=1)

        # retrieve only training columns
        self.data_info_training = self.data_info[self.training_col]

        # predict prob
        # print("\nProbability distribution exhibited in plot：")
        self.data_info_pred = self.rf_gcv.predict_proba(self.data_info_training)
        # plt.hist(self.data_info_pred[:, 1])

        print("---------------------------------------------")
        print("data_info contained %d records" % (len(self.data_info)))
        # predict label
        self.data_info_pred_thre = self.rf_gcv.predict(self.data_info_training)  # 概率形式预测测试集的类别
        print(f"Prediction complete, predicting result 'data_info_pred_thre' contained {len(self.data_info_pred_thre)} records")
        self.data_info_positive_thre = self.data_info[self.data_info_pred_thre == 1]
        self.data_info_negative_thre = self.data_info[self.data_info_pred_thre == 0]
        print("Total mutations count: %d, positive mutations count: %d, proportion for positive records: 1:%d" % (
            len(self.data_info), len(self.data_info_positive_thre),
            len(self.data_info) / len(self.data_info_positive_thre)))

        return(self.data_info_pred, self.data_info_pred_thre)

    def unpickle(self, instance_path):
        """load local model
        """
        with open(instance_path, 'rb') as f:
            return joblib.load(f)

class model_analyze_with_DNA(object):
    """Using provided information to complete the integration of RNA and DNA mutations.

    Attributes:
        gdc_validate_all_info: The DataFrame containing required six columns (Tumor_Sample_UUID Chromosome Start_Position Reference_Allele Tumor_Allele1 Tumor_Allele2) for DNA evidence
        RNA_info: The DataFrame containing all features and annotations for RNA somatic mutations from "own_data_vcf_info_retriver.py".
        WXS_target_interval: Interval set for WXS target regions.
        exon_interval: Interval set for exon regions.
        num_threads: Number of available threads.
    """

    # required columns
    predict_required_columns = ['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Allele1', "Hugo_Symbol",
                        'ref_AD_tumor_RNA', 'alt_AD_tumor_RNA', 'ref_AD_normal', 'alt_AD_normal', 'ref_AD_tumor_DNA', 'alt_AD_tumor_DNA',
                        "Tumor_Sample_UUID", 'record_filter', 'Protein_Change']
    gdc_required_columns = ['Tumor_Sample_UUID', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Allele1', 'Tumor_Allele2']

    def __init__(self, gdc_validate_all_info, RNA_info, WXS_target_interval, exon_interval, num_threads):
        """Init data and target intervals
        """

        self.gdc_validate_all_info = gdc_validate_all_info[['Tumor_Sample_UUID', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Allele1', 'Tumor_Allele2']]
        self.RNA_info = RNA_info
        self.WXS_target_interval = WXS_target_interval
        self.exon_interval = exon_interval
        self.num_threads = num_threads

    def GDC_RNA_discrepancy_analysis_total(self):
        """Using RNA mutations to split and elaborately interpret DNA mutations.

        :param gdc_validate_all_info: DNA mutations with essential info
        :param positive_RNA_info: predicted positive RNA mutations with all info
        :param negative_RNA_info: predicted negative RNA mutations with all info

        :return positive_RNA_info_GDC_intersect: positive RNA-DNA overlap part -- best part
        :return gdc_validate_all_info_DNA_only_RNA_negative: negative RNA-DNA overlap part -- our model's defect
        :return gdc_validate_all_info_DNA_only_RNA_missing: DNA only (RNA missing) part -- RNA not-captured part, required explanation
        :return positive_RNA_info_RNA_only: positive RNA only part -- RNA rescued part
        :return negative_RNA_info_RNA_only: negative RNA only part -- inherent erroneous part
        """

        print(f"Start to analyze the intersection and difference set for DNA and RNA mutations.")

        print("="*100)

        # analyze case status for DNA and RNA respectively
        GDC_set = set(self.gdc_validate_all_info.Tumor_Sample_UUID)
        RNA_set = set(self.RNA_info.Tumor_Sample_UUID)
        GDC_RNA_intersection_set = set(self.gdc_validate_all_info.Tumor_Sample_UUID).intersection(set(self.RNA_info.Tumor_Sample_UUID))
        print(f"DNA contained {len(GDC_set)} cases；RNA contained {len(RNA_set)} cases；DNA、RNA had {len(GDC_RNA_intersection_set)} intersected cases")

        print(f"P.S: {RNA_set - GDC_RNA_intersection_set} lacked DNA evidence support "
              f"which were excluded in subsequent analysis to avoid impact on Precision-Recall calculation, but special attention was required.")
        self.gdc_validate_all_info = self.gdc_validate_all_info[self.gdc_validate_all_info['Tumor_Sample_UUID'].isin(GDC_RNA_intersection_set)]
        self.RNA_info = self.RNA_info[self.RNA_info['Tumor_Sample_UUID'].isin(GDC_RNA_intersection_set)]
        print(f"Finally, we had {len(self.RNA_info)} RNA somatic mutations for integration analysis.")

        # primary assessment using overlapping status
        self.RNA_info_TP = pd.merge(self.RNA_info, self.gdc_validate_all_info, on=list(self.gdc_validate_all_info.columns))
        print(f"Attribute: RNA_info_TP (with DNA evidence) fetched, mutations count: {len(self.RNA_info_TP)}")
        self.RNA_info_TN = self.RNA_info.append(self.RNA_info_TP)
        self.RNA_info_TN = self.RNA_info_TN.drop_duplicates(keep=False)
        print(f"Attribute: RNA_info_TN (without DNA evidence) fetched, mutations count: {len(self.RNA_info_TN)}")
        self.RNA_info_TP['DNA_label'] = 1
        self.RNA_info_TN['DNA_label'] = 0
        # re-combine TP and TN records
        self.RNA_info = self.RNA_info_TP.append(self.RNA_info_TN, ignore_index=True)
        self.positive_RNA_info = self.RNA_info[self.RNA_info['pred_label'] == 1]
        self.negative_RNA_info = self.RNA_info[self.RNA_info['pred_label'] == 0]

        print("=" * 100)

        # RNA mutations confirmed, determine DNA mutations subsequently
        self.gdc_validate_all_info = self.gdc_validate_all_info[self.gdc_validate_all_info.Tumor_Sample_UUID.isin(GDC_RNA_intersection_set)]
        print(f"DNA mutations included within our integrative analysis: {len(self.gdc_validate_all_info)}")
        self.gdc_validate_all_info_outside_region = self.gdc_validate_all_info[self.gdc_validate_all_info.apply(lambda x : False if (self.WXS_target_interval[x['Chromosome']].contains(x['Start_Position'])) & (self.exon_interval[x['Chromosome']].contains(x['Start_Position'])) else True, axis=1)]
        self.gdc_validate_all_info = self.gdc_validate_all_info.append(self.gdc_validate_all_info_outside_region)
        self.gdc_validate_all_info = self.gdc_validate_all_info.drop_duplicates(keep=False)
        print(f"P.S: DNA mutations not located within our research regions (gdc_validate_all_info_outside_region): {len(self.gdc_validate_all_info_outside_region)}")
        print(f"Finally, we had {len(self.gdc_validate_all_info)} DNA somatic mutations for integration analysis.")

        print("=" * 100)

        # RNA and DNA mutations confirmed, start to exclude special mutations
        print(f"Out of {len(self.positive_RNA_info)} positive and {len(self.negative_RNA_info)} negative RNA somatic mutations. ")
        RNA_only_set = set(self.positive_RNA_info.Tumor_Sample_UUID) - GDC_RNA_intersection_set
        print(f"Positive RNA mutations contained {len(set(self.positive_RNA_info.Tumor_Sample_UUID))} cases totally, but {RNA_only_set} had no DNA evidence (P-R cannot be assessed, ignore)")
        self.positive_RNA_info_RNA_case_only = self.positive_RNA_info[self.positive_RNA_info.Tumor_Sample_UUID.isin(RNA_only_set)]
        self.positive_RNA_info = self.positive_RNA_info[~self.positive_RNA_info.Tumor_Sample_UUID.isin(RNA_only_set)]
        self.negative_RNA_info = self.negative_RNA_info[~self.negative_RNA_info.Tumor_Sample_UUID.isin(RNA_only_set)]
        print(f"As a sub-part of positive RNA mutations, these without DNA evidence (positive_RNA_info_RNA_case_only) contained {len(self.positive_RNA_info_RNA_case_only)} mutations. ")
        print(f"Finally, we had {len(self.positive_RNA_info)} positive and {len(self.negative_RNA_info)} negative RNA somatic mutations for integration analysis.")

        print("=" * 100)

        # First, we studied positive RNA-DNA overlapping part: positive_RNA_info_GDC_intersect
        self.positive_RNA_info_GDC_intersect = self.positive_RNA_info.merge(self.gdc_validate_all_info, how = 'inner')
        positive_RNA_info_GDC_intersect_case_set = set(self.positive_RNA_info_GDC_intersect.Tumor_Sample_UUID.value_counts().index)
        print(f"Positive RNA-DNA overlapping part: positive_RNA_info_GDC_intersect contained {len(self.positive_RNA_info_GDC_intersect)} mutations and spread over {len(positive_RNA_info_GDC_intersect_case_set)} cases")

        print("="*100)

        # Second, we studied mutations theoretically within DNA only (due to labeled negative or missing)
        # gdc_validate_all_info_DNA_only_RNA_negative, gdc_validate_all_info_DNA_only_RNA_missing
        self.gdc_validate_all_info_DNA_only = self.gdc_validate_all_info.append(self.positive_RNA_info_GDC_intersect)
        self.gdc_validate_all_info_DNA_only = self.gdc_validate_all_info_DNA_only.drop_duplicates(subset=['Tumor_Sample_UUID', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Allele1', 'Tumor_Allele2'], keep=False)
        self.gdc_validate_all_info_DNA_only = self.gdc_validate_all_info_DNA_only.dropna(axis='columns', how="all") # drop NA rows
        print(f"Theoretically DNA only mutations (mislabeled negative or missing): gdc_validate_all_info_DNA_only contained {len(self.gdc_validate_all_info_DNA_only)} mutations and spread over {len(set(self.gdc_validate_all_info_DNA_only.Tumor_Sample_UUID))} cases")
        print("Further split into detailed sub-part.")
        # labeled negative
        self.gdc_validate_all_info_DNA_only_RNA_negative = self.gdc_validate_all_info_DNA_only.merge(self.negative_RNA_info, how = 'inner')
        print(f"These mislabeled negative RNA mutations: gdc_validate_all_info_DNA_only_RNA_negative contained {len(self.gdc_validate_all_info_DNA_only_RNA_negative)} mutations and spread over {len(set(self.gdc_validate_all_info_DNA_only_RNA_negative.Tumor_Sample_UUID))} cases")
        # missing in RNA
        self.gdc_validate_all_info_DNA_only_RNA_missing = self.gdc_validate_all_info_DNA_only.append(self.gdc_validate_all_info_DNA_only_RNA_negative)
        self.gdc_validate_all_info_DNA_only_RNA_missing = self.gdc_validate_all_info_DNA_only_RNA_missing.drop_duplicates(subset=['Tumor_Sample_UUID', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Allele1', 'Tumor_Allele2'], keep=False)
        self.gdc_validate_all_info_DNA_only_RNA_missing = self.gdc_validate_all_info_DNA_only_RNA_missing.dropna(axis='columns', how="all")  # 去除所有值均为空的列（RNA中所有注释列与特征列）
        print(f"These missing DNA mutations: gdc_validate_all_info_DNA_only_RNA_missing contained {len(self.gdc_validate_all_info_DNA_only_RNA_missing)} mutations and spread over {len(set(self.gdc_validate_all_info_DNA_only_RNA_missing.Tumor_Sample_UUID))} cases")

        print("=" * 100)

        # Finally, we studied RNA mutations without any DNA evidence
        # positive RNA_only
        self.positive_RNA_info_RNA_only = self.positive_RNA_info.append(self.positive_RNA_info_GDC_intersect)
        self.positive_RNA_info_RNA_only = self.positive_RNA_info_RNA_only.drop_duplicates(keep=False)
        # negative RNA_only
        self.negative_RNA_info_RNA_only = self.negative_RNA_info.append(self.gdc_validate_all_info_DNA_only_RNA_negative)
        self.negative_RNA_info_RNA_only = self.negative_RNA_info_RNA_only.drop_duplicates(keep=False)

        print(f"Positive RNA only mutations: positive_RNA_info_RNA_only contained {len(self.positive_RNA_info_RNA_only)} mutations and spread over {len(set(self.positive_RNA_info_RNA_only.Tumor_Sample_UUID))} cases; "
              f"And negative RNA only mutations: negative_RNA_info_RNA_only contained {len(self.negative_RNA_info_RNA_only)} mutations and spread over {len(set(self.negative_RNA_info_RNA_only.Tumor_Sample_UUID))} cases")

        print("DNA_RNA_discrepancy_analysis finished.\n")

    def GDC_RNA_missing_check(self, RNA_calling_info, RNA_bam_folder_loc, Mutect2_target_detected_sites):
        """Check DNA only part (gdc_validate_all_info_DNA_only_RNA_missing) and annotate them with RNA coverage&mapping quality

        :param RNA_calling_info: Location for RNA somatic mutation calling info
        :param RNA_bam_folder_loc: Location for RNA bam file storage
        :param Mutect2_target_detected_sites: Location for Mutect2 called RNA somatic mutations (generated along with framework)
        :return gdc_validate_all_info_DNA_only_RNA_missing_RNA_info: DNA only mutation info with annotation of RNA sites
        (ref_AD_tumor_bam_RNA、alt_AD_tumor_bam_RNA、other_AD_tumor_bam_RNA、median_MQ_tumor_bam_RNA、other_median_MQ_tumo_bam_RNA)
        """

        print(f"For {len(self.gdc_validate_all_info_DNA_only_RNA_missing)} missing DNA mutations")
        print("We chose to analyze these mutations' coverage and other info within raw RNA sequencing data......")

        print("="*100)

        # start multi-process to accelerate the retrieval of annotation info
        manager = Manager()
        gdc_validate_all_info_DNA_only_RNA_missing_case_list = manager.list()

        print('Parent process %s.' % os.getpid())
        p = Pool(self.num_threads)
        for case_id in self.gdc_validate_all_info_DNA_only_RNA_missing.Tumor_Sample_UUID.value_counts().index:
            p.apply_async(func=RNA_bam_info_adder, args=(case_id, gdc_validate_all_info_DNA_only_RNA_missing_case_list, self.gdc_validate_all_info_DNA_only_RNA_missing[self.gdc_validate_all_info_DNA_only_RNA_missing.Tumor_Sample_UUID == case_id], RNA_calling_info, RNA_bam_folder_loc, tumor_bam_record_retrive,), error_callback=print_error)
        print('Waiting for all subprocesses done...')
        p.close()
        p.join()
        print('All subprocesses done.')

        gdc_validate_all_info_DNA_only_RNA_missing_case_list = list(gdc_validate_all_info_DNA_only_RNA_missing_case_list)
        self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info = pd.concat(gdc_validate_all_info_DNA_only_RNA_missing_case_list)
        print(f"Finished the retrieval of annotation for DNA only mutations: gdc_validate_all_info_DNA_only_RNA_missing_RNA_info")

        # sub-part 0: Mutect2 called, but not present in positive or negative mutation sets (falsly filtered by multi-filtering strategy)
        self.Mutect2_target_detected_sites = Mutect2_target_detected_sites
        gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_temp = self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info.copy(deep=True)
        self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_multi_filtered_sites = pd.merge(gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_temp, self.Mutect2_target_detected_sites)
        print(f"Multi-filtering strategy falsly filtered {len(self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_multi_filtered_sites)} mutations，"
              f"Attribute: gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_multi_filtered_sites")

        # sub-part 1: Mutect2 cannot call, missing in RNA
        print(f"Mutations which Mutect2 cannot call and required force-calling to check their status: {len(self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info)}; "
              f"Attribute: gdc_validate_all_info_DNA_only_RNA_missing_RNA_info\n")

    def GDC_Mutect2_missing_analysis(self, GDC_Mutect2_missing_force_call_RNA_info):
        """Analyze status for force-called mutations and prepare for prediction

        :param GDC_Mutect2_missing_force_call_RNA_info: Mutect2 force-called mutations
        :return gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info_final
        """

        # remove duplicate records due to Funcotator's annotation
        self.GDC_Mutect2_missing_force_call_Mutect2_info = GDC_Mutect2_missing_force_call_RNA_info.drop_duplicates(subset=['Chromosome', 'Start_Position', 'Tumor_Sample_UUID'], keep='first', inplace=False)
        # check basic status
        print(f"Mutect2 force-called {len(self.GDC_Mutect2_missing_force_call_Mutect2_info)} mutations")
        self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_force_call_failed = self.GDC_Mutect2_missing_force_call_Mutect2_info.append(self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info)
        self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_force_call_failed = self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_force_call_failed.drop_duplicates(
            subset=['Tumor_Sample_UUID', 'Chromosome', 'Start_Position'], keep=False)
        print(f"However, there're still {len(self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_force_call_failed)} mutations had not been force-called. Attribute: gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_force_call_failed")
        # add force-called sites' info back into initial DNA only part
        self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info.drop(['Reference_Allele', 'Tumor_Allele1', 'Tumor_Allele2'], axis=1, inplace=True)
        self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info = self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info.merge(self.GDC_Mutect2_missing_force_call_Mutect2_info)

        print(f"Finally, we added force-called mutations' features back into DNA only part: {len(self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info)}. Attribute: gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info")

        print("="*100)

        print("Start to predict these force-called mutations...")

        temp_gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info = self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info.copy(deep=True)

        # sub-part 0: non-SNP
        self.GDC_Mutect2_missing_force_call_RNA_info_non_SNP = temp_gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info[temp_gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info.apply(lambda x : True if (len(x['Reference_Allele'])>1)or(len(x['Tumor_Allele1'])>1) else False, axis=1)]
        print(f"Be advised, {len(self.GDC_Mutect2_missing_force_call_RNA_info_non_SNP)} mutations were recognized as non-SNP (indel) within Mutect2 force-called results, cannot be predicted. Attribute: GDC_Mutect2_missing_force_call_RNA_info_non_SNP")
        temp_gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info = temp_gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info[temp_gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info.apply(lambda x : False if (len(x['Reference_Allele'])>1)or(len(x['Tumor_Allele1'])>1) else True, axis=1)]

        # sub-part 1: SNP
        self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info_final = temp_gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info
        print(f"All non-SNP mutations had been excluded for prediction, finally we had {len(self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info_final)} mutations for prediction. Attribute: gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info_final")

    def GDC_RNA_info_summary(self, data_prepare_demo):
        """Add tag&sub_tag info for all integrated RNA&DNA mutations and finish aggregating

        :param data_prepare_demo: Data preparation class.
        :return final_info: Final table containing all available information
        """

        print("="*100)

        # For RNA-DNA overlap part
        print(f"\nFirst, we conducted the combination for RNA-DNA overlap part and assigned 'RNA_DNA_overlap' tag. \n"
              f"Positive: positive_RNA_info_GDC_intersect, {len(self.positive_RNA_info_GDC_intersect)}; \n"
              f"Negative: gdc_validate_all_info_DNA_only_RNA_negative, {len(self.gdc_validate_all_info_DNA_only_RNA_negative)}.")
        RNA_DNA_overlap = pd.concat([self.positive_RNA_info_GDC_intersect, self.gdc_validate_all_info_DNA_only_RNA_negative])
        RNA_DNA_overlap['Tag'] = "RNA_DNA_overlap"

        # For RNA only part
        print(f"\nSecond, we conducted the combination for RNA only part and assigned 'RNA_only' tag."
              f"Positive: positive_RNA_info_RNA_only, {len(self.positive_RNA_info_RNA_only)}; \n"
              f"Negative: negative_RNA_info_RNA_only, {len(self.negative_RNA_info_RNA_only)}.")
        RNA_only = pd.concat([self.positive_RNA_info_RNA_only, self.negative_RNA_info_RNA_only])
        RNA_only['Tag'] = "RNA_only"

        # For DNA only part
        # only force_called mutations can be predicted
        print(f"\nFinally, we conducted the combination for DNA only part and assigned 'DNA_only' tag. However, there're different categories which required sub-tags. \n"
              f"non-SNP(Sub_Tag): GDC_Mutect2_missing_force_call_RNA_info_non_SNP, {len(self.GDC_Mutect2_missing_force_call_RNA_info_non_SNP)}\n"
              f"force_call_failed(Sub_Tag): gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_force_call_failed, {len(self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_force_call_failed)}\n"
              f"force_called(Sub_Tag): gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info_final, {len(self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info_final)}\n")
        self.GDC_Mutect2_missing_force_call_RNA_info_non_SNP['Sub_Tag'] = "non_SNP"
        self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_force_call_failed['Sub_Tag'] = "force_call_failed"
        self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info_final['Sub_Tag'] = "force_called"
        DNA_only = pd.concat([self.GDC_Mutect2_missing_force_call_RNA_info_non_SNP, self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_force_call_failed, self.gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_info_final])
        DNA_only['Tag'] = "DNA_only"

        print("=" * 100)

        # Add prediction probability and label columns

        # get prob and label
        predictable_part = pd.concat([RNA_DNA_overlap, RNA_only, DNA_only[DNA_only['Sub_Tag']=="force_called"]])
        pred_prob, pred_label = data_prepare_demo.model_predict(predictable_part)
        predictable_part['pred_prob'] = pred_prob[:, 1]
        predictable_part['pred_label'] = pred_label
        # merge into final table
        self.final_info = pd.concat([predictable_part, DNA_only[DNA_only['Sub_Tag']!="force_called"]])
        self.final_info.reset_index(drop=True, inplace=True)

        print(f"Integrative analysis with DNA evidence finished. For RNA somatic mutations, "
              f"their precision was{len(self.positive_RNA_info_GDC_intersect)/(len(self.positive_RNA_info_GDC_intersect)+len(self.positive_RNA_info_RNA_only))}, "
              f"recall was{len(self.positive_RNA_info_GDC_intersect)/(len(self.positive_RNA_info_GDC_intersect)+len(self.gdc_validate_all_info_DNA_only_RNA_negative))}")

        print(f"\nFinal table (final_info) aggregating all RNA DNA somatic mutation info has been constructed, contained {len(self.final_info)} records.")
        print("Tag——RNA_DNA_overlap, RNA_only, DNA_only")
        print("Sub_Tag——non_SNP, force_called, force_call_failed")
        print("Predicted probability column：pred_prob")
        print("Predicted label column：pred_label")

    def vcf_info_extractor(self, total_mutation_set, template_vcf_file, output_folder):
        """Extracted sites for Mutect2 force-calling

        :param total_mutation_set: mutation dataset
        :param template_vcf_file: vcf file as a template
        :param output_folder: folder path for vcf files
        :return: vcf files only contained sites and variants (CHROM、POS、REF、ALT)
        """

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)  # new folder
        # read in template vcf file
        vcf_reader = vcf.Reader(open(template_vcf_file, 'r'))
        record_list = [record for record in vcf_reader]
        # extract case's sites info respectively
        print(f"Start to extract single cases' sites information into vcf files from {len(total_mutation_set)} mutations to be force-called...")
        # traverse each case
        for case_id in total_mutation_set.Tumor_Sample_UUID.value_counts().index:
            case_mutation_set = total_mutation_set[total_mutation_set.Tumor_Sample_UUID==case_id]
            print(f"{case_id} had {len(case_mutation_set)} mutations waiting to be extracted...")
            # new vcf file
            vcf_writer = vcf.Writer(open(os.path.join(output_folder, f"{case_id}.vcf"), 'w'), vcf_reader)
            for i in case_mutation_set.index:
                record_list[0].CHROM = case_mutation_set.loc[i, "Chromosome"]
                record_list[0].POS = case_mutation_set.loc[i, "Start_Position"]
                record_list[0].REF = case_mutation_set.loc[i, "Reference_Allele"]
                record_list[0].ALT = list(case_mutation_set.loc[i, "Tumor_Allele1"])
                vcf_writer.write_record(record_list[0])
            vcf_writer.close()
            print(f"Sites extraction complete, corresponding vcf info had been stored within {os.path.join(output_folder)}")
        # store case IDs into same folder
        pd.Series(list(total_mutation_set.Tumor_Sample_UUID.value_counts().index)).to_csv(os.path.join(output_folder, f"case_info.table"), index=False)

    def pickle(self, output_file_path):
        """Save model or instances into given file path
        """
        with open(output_file_path, 'wb') as f:
            joblib.dump(self, f, compress=3)

    def print_error(self, value):
        """Self-specified action for multi-process's error_callback
        """
        print("error: ", value)

if __name__ == '__main__':

    if args.step == 1:
        # pre-process #
        # read in mutation info
        DNA_info = pd.read_table(args.DNA_info)
        RNA_info = pd.read_table(args.RNA_info)
        # read in interval info for WES and exon
        data_prepare_demo = data_prepare()
        WXS_target_interval_path = args.WXS_target_interval
        exon_interval_path = args.exon_interval
        data_prepare_demo.WXS_exon_region_dict_generate(WXS_target_interval_path, exon_interval_path, args.num_threads)

        # extra process (for GDC maf file input only)
        # DNA_info = data_prepare_demo.GDC_site_info_retrieve(DNA_info)

        # analysis #
        # init
        result_analysis_demo = model_analyze_with_DNA(DNA_info,
                                                   RNA_info,
                                                   data_prepare_demo.WXS_target_interval,
                                                   data_prepare_demo.exon_interval,
                                                   args.num_threads)

        # analyze the overlap and difference status for RNA and DNA mutations
        result_analysis_demo.GDC_RNA_discrepancy_analysis_total()

        # check the RNA coverage info for DNA only mutations
        RNA_calling_info = pd.read_table(args.RNA_calling_info)
        RNA_bam_folder_loc = args.RNA_bam_folder
        Mutect2_target_detected_sites = pd.read_table(args.Mutect2_target_detected_sites)
        result_analysis_demo.GDC_RNA_missing_check(RNA_calling_info, RNA_bam_folder_loc, Mutect2_target_detected_sites)

        # export the vcf files for Mutect2 force-calling
        result_analysis_demo.vcf_info_extractor(result_analysis_demo.gdc_validate_all_info_DNA_only_RNA_missing,
                                                args.template_vcf_file,
                                                os.path.join(args.project_folder, args.cancer_type, "RNA/RNA_somatic_mutation/MAFToVCF/DNA_only_RNA_missing_Mutect2_check"))

        # save the analysis class for step 1
        result_analysis_demo.pickle(args.output_file_path)

        print("Step 1 completed, please proceed to prepare for step 2. ")

    if args.step == 2:
        # pre-process #
        # read in mutations info after force-calling
        DNA_only_RNA_missing_RNA_info_Mutect2_force_call = pd.read_table(args.force_call_RNA_info)

        # read in the analysis class for step 1
        data_prepare_demo = data_prepare()
        result_analysis_demo = data_prepare_demo.unpickle(args.instance_path)

        # read in prediction model
        data_prepare_demo.model_predict_interpret_prepare(args.model_path, args.one_hot_encoder_path, args.training_columns_path)

        # analysis #
        # add force-calling info into our table
        result_analysis_demo.GDC_Mutect2_missing_analysis(DNA_only_RNA_missing_RNA_info_Mutect2_force_call)

        # generate and export the final table
        result_analysis_demo.GDC_RNA_info_summary(data_prepare_demo)
        result_analysis_demo.final_info.to_csv(args.output_file_path, sep="\t", index=False)