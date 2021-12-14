# python /home/lqh/Codes/Python/RNA-SSNV/lib/tools/maf2vcf.py \
# --total_mutation_set /home/lqh/Codes/Python/RNA-SSNV/output/RNA_overlooked_cancer_driver_genes_mutations.table \
# --template_vcf_file /home/lqh/Codes/Python/RNA-SSNV/resources/template_simple.vcf \
# --output_vcf_path /home/lqh/Codes/Python/RNA-SSNV/output/RNA_overlooked_cancer_driver_genes_mutations.vcf

import os
import vcf
import pandas as pd

import argparse

parser=argparse.ArgumentParser(description="A discriminate model construction pipeline for RNA-SSNV.")
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.') # demo
# Generic parameter
parser.add_argument('--total_mutation_set', help='Total mutation set.')
parser.add_argument("--template_vcf_file", help="Path for template vcf file.")
parser.add_argument("--output_vcf_path", help="Path for final output file.")

args=parser.parse_args()

def total_vcf_info_extractor(total_mutation_set, template_vcf_file, output_vcf_path):
    """Extracted sites for Mutect2 force-calling

    :param total_mutation_set: mutation dataset
    :param template_vcf_file: vcf file as a template
    :param output_vcf_path: path for vcf file
    :return: vcf files only contained sites and variants (CHROM、POS、REF、ALT)
    """
    # read in template vcf file
    vcf_reader = vcf.Reader(open(template_vcf_file, 'r'))
    record_list = [record for record in vcf_reader]
    # new vcf file
    vcf_writer = vcf.Writer(open(output_vcf_path, 'w'), vcf_reader)
    for i in total_mutation_set.index:
        record_list[0].CHROM = total_mutation_set.loc[i, "Chromosome"]
        record_list[0].POS = total_mutation_set.loc[i, "Start_Position"]
        record_list[0].REF = total_mutation_set.loc[i, "Reference_Allele"]
        record_list[0].ALT = list(total_mutation_set.loc[i, "Tumor_Allele1"])
        vcf_writer.write_record(record_list[0])
    vcf_writer.close()

def vcf_info_extractor(total_mutation_set, template_vcf_file, output_folder):
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
    print(
        f"Start to extract single cases' sites information into vcf files from {len(total_mutation_set)} mutations to be force-called...")
    # traverse each case
    for case_id in total_mutation_set.Tumor_Sample_UUID.value_counts().index:
        case_mutation_set = total_mutation_set[total_mutation_set.Tumor_Sample_UUID == case_id]
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
    pd.Series(list(total_mutation_set.Tumor_Sample_UUID.value_counts().index)).to_csv(
        os.path.join(output_folder, f"case_info.table"), index=False)

if __name__ == '__main__':
    total_mutation_set = pd.read_table(args.total_mutation_set)

    total_vcf_info_extractor(total_mutation_set, args.template_vcf_file, args.output_vcf_path)