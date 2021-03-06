"""
WXS short variant discovery (SNPs + Indels)
https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-WXSseq-short-variant-discovery-SNPs-Indels-
"""
import pandas as pd
import os
from snakemake.utils import min_version

min_version("3.2")

workdir: config["workDir"]+"/"+config["cancerType"]+"/"+"WXS"

# Prepare corresponding sample info
samples = pd.read_table(config["samples"]).set_index("aliquots_id", drop=False)
tumor_samples = samples.loc[samples["sample_type"]=="Primary Tumor", ]
normal_samples = samples.loc[samples["sample_type"].isin(["Solid Tissue Normal", "Blood Derived Normal"]), ]
best_normal_samples = samples.loc[samples["sample_type"]=="Best Normal", ].set_index("case_id", drop=False)

# Prepare lists to iterate over
WXS_primary_tumor_aliquots_id = tumor_samples.index.values
case_ids = set(samples["case_id"])

rule all:
    input:
         # "WXS_germline_mutation/KGG_annotate/"+config["sampleType"]+"_SNP.flt.xlsx",
         # expand("apply_BQSR/{aliquots_id}.bam", aliquots_id = WXS_primary_tumor_aliquots_id),
         expand("WXS_somatic_mutation/FilterMutectCalls/{case_id}.vcf.gz", case_id = case_ids),
         expand("WXS_somatic_mutation/SelectVariants/PASS_SNP/{case_id}.vcf.gz", case_id = case_ids),
         expand("WXS_somatic_mutation/SelectVariants/PASS_SNP_WES_Interval/{case_id}.vcf.gz", case_id = case_ids),
         expand("WXS_somatic_mutation/SelectVariants/PASS_SNP_WES_Interval_exon/{case_id}.vcf.gz", case_id = case_ids),
         expand("WXS_somatic_mutation/VariantsToTable/PASS_SNP_WES_Interval_exon/{case_id}.table", case_id = case_ids)

# takes 22h
# cores 2
rule Mutect2:
    input:
        ref = config["ref"]["genome"],
        gnomad = config["ref"]["variant"]["gnomad"],
        PoN = config["ref"]["variant"]["PoN"],
        tumor_bams = lambda wildcards : [os.path.join(config["sampleDir"],tumor_samples["file_id"][tumor_aliquots_id],tumor_samples["file_name"][tumor_aliquots_id])
                                        for tumor_aliquots_id in tumor_samples.loc[tumor_samples["case_id"]==wildcards.case_id, ].index],
        normal_bams = lambda wildcards : [os.path.join(config["normalSampleDir"],normal_samples["file_id"][normal_aliquots_id],normal_samples["file_name"][normal_aliquots_id])
                                        for normal_aliquots_id in normal_samples.loc[normal_samples["case_id"]==wildcards.case_id, ].index]
    output:
        vcf = protected("WXS_somatic_mutation/Mutect2/{case_id}.vcf.gz"),
        f1r2 = protected("WXS_somatic_mutation/Mutect2/{case_id}.f1r2.tar.gz")
    params:
        tumor_bams = lambda wildcards : ["-I "+os.path.join(config["sampleDir"],tumor_samples["file_id"][tumor_aliquots_id],tumor_samples["file_name"][tumor_aliquots_id])
                                        for tumor_aliquots_id in tumor_samples.loc[tumor_samples["case_id"]==wildcards.case_id, ].index],
        normal_bams = lambda wildcards : ["-I "+os.path.join(config["normalSampleDir"],normal_samples["file_id"][normal_aliquots_id],normal_samples["file_name"][normal_aliquots_id])
                                        for normal_aliquots_id in normal_samples.loc[normal_samples["case_id"]==wildcards.case_id, ].index],
        normal_sample_name = lambda wildcards : expand("-normal {normal_sample_names}",
                                                       normal_sample_names = normal_samples.loc[normal_samples["case_id"]==wildcards.case_id, "aliquots_id"]),
        threads = 2
    threads: 2
    log:
        "logs/Mutect2/{case_id}.log"
    shell:
        """
         gatk Mutect2 \
         -R {input.ref} \
         {params.tumor_bams} \
         {params.normal_bams} \
         {params.normal_sample_name} \
         --germline-resource {input.gnomad} \
         --panel-of-normals {input.PoN} \
         --native-pair-hmm-threads {params.threads} \
         --f1r2-tar-gz {output.f1r2} \
         -O {output.vcf}
        """

# takes
# cores 2
rule FilterMutectCalls_combined:
    input:
        ref = config["ref"]["genome"],
        tumor_bams = lambda wildcards : [os.path.join(config["sampleDir"],tumor_samples["file_id"][tumor_aliquots_id],tumor_samples["file_name"][tumor_aliquots_id])
                                        for tumor_aliquots_id in tumor_samples.loc[tumor_samples["case_id"]==wildcards.case_id, ].index],
        best_normal_bam = lambda wildcards : os.path.join(config["normalSampleDir"],best_normal_samples["file_id"][wildcards.case_id],best_normal_samples["file_name"][wildcards.case_id]),
        vcf = "WXS_somatic_mutation/Mutect2/{case_id}.vcf.gz",
        f1r2 = "WXS_somatic_mutation/Mutect2/{case_id}.f1r2.tar.gz",
        gnomad_SNP = config["ref"]["variant"]["gnomad_SNP"]
    output:
        best_normal_pileups_table = "WXS_somatic_mutation/GetPileupSummaries/{case_id}-best-normal-pileups.table",
        read_orientation_model = "WXS_somatic_mutation/Mutect2/{case_id}.read-orientation-model.tar.gz",
        vcf = protected("WXS_somatic_mutation/FilterMutectCalls/{case_id}.vcf.gz")
    params:
        tumor_bams = lambda wildcards : [os.path.join(config["sampleDir"],tumor_samples["file_id"][tumor_aliquots_id],tumor_samples["file_name"][tumor_aliquots_id])
                                        for tumor_aliquots_id in tumor_samples.loc[tumor_samples["case_id"]==wildcards.case_id, ].index]
    threads: 2
    log:
        "logs/FilterMutectCalls_combined/{case_id}.log"
    run:
        commands = []
        # For Best Normal GetPileupSummaries
        normal_GetPileupSummaries_command = f"gatk GetPileupSummaries -I {input.best_normal_bam} -L {input.gnomad_SNP} -V {input.gnomad_SNP} -O {output.best_normal_pileups_table}"
        commands.append(normal_GetPileupSummaries_command)
        # For --contamination-table info
        contamination_table = ""
        # For --tumor-segmentation info
        tumor_segmentation = ""
        for tumor_bam in list(input.tumor_bams):
            print("\n"+tumor_bam)
            tumor_bam_file_name = os.path.split(tumor_bam)[-1]
            # For Tumor GetPileupSummaries
            GetPileupSummaries_command = f"gatk GetPileupSummaries " \
                                         f"-I {tumor_bam} " \
                                         f"-L {input.gnomad_SNP} " \
                                         f"-V {input.gnomad_SNP} " \
                                         f"-O WXS_somatic_mutation/GetPileupSummaries/{tumor_bam_file_name}-pileups.table"
            print(GetPileupSummaries_command)
            commands.append(GetPileupSummaries_command)
            # For CalculateContaminations
            CalculateContaminations_command = f"gatk CalculateContamination " \
                                              f"-I WXS_somatic_mutation/GetPileupSummaries/{tumor_bam_file_name}-pileups.table " \
                                              f"-matched {output.best_normal_pileups_table} " \
                                              f"-O WXS_somatic_mutation/GetPileupSummaries/{tumor_bam_file_name}-contamination.table " \
                                              f"--tumor-segmentation WXS_somatic_mutation/GetPileupSummaries/{tumor_bam_file_name}-segments.table"
            print(CalculateContaminations_command)
            commands.append(CalculateContaminations_command)
            # Add corresponding table info
            contamination_table = contamination_table + f" --contamination-table WXS_somatic_mutation/GetPileupSummaries/{tumor_bam_file_name}-contamination.table"
            tumor_segmentation = tumor_segmentation + f" --tumor-segmentation WXS_somatic_mutation/GetPileupSummaries/{tumor_bam_file_name}-segments.table"
        print("--------------------------------------------------")
        # For LearnReadOrientationModel
        LearnReadOrientationModel_command = f"gatk LearnReadOrientationModel -I {input.f1r2} -O {output.read_orientation_model}"
        print(LearnReadOrientationModel_command)
        commands.append(LearnReadOrientationModel_command)
        # For FilterMutectCalls
        FilterMutectCalls_command = f"gatk FilterMutectCalls " \
                                    f"-R {input.ref} " \
                                    f"-V {input.vcf} " \
                                    f"{contamination_table} " \
                                    f"{tumor_segmentation} " \
                                    f"--orientation-bias-artifact-priors {output.read_orientation_model} " \
                                    f"-O {output.vcf}"
        print(FilterMutectCalls_command)
        commands.append(FilterMutectCalls_command)

        for c in commands:
            shell(c)

rule SelectVariants:
    input:
        vcf = "WXS_somatic_mutation/FilterMutectCalls/{case_id}.vcf.gz"
    output:
        vcf = "WXS_somatic_mutation/SelectVariants/PASS/{case_id}.vcf.gz"
    threads: 2
    log:
        "logs/SelectVariants/{case_id}.log"
    shell:
        """
        gatk SelectVariants \
        -V {input.vcf} \
        -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY \
        --exclude-filtered true \
        -O {output.vcf}
        """

rule SelectVariants_SNP:
    input:
        vcf = "WXS_somatic_mutation/SelectVariants/PASS/{case_id}.vcf.gz"
    output:
        vcf = "WXS_somatic_mutation/SelectVariants/PASS_SNP/{case_id}.vcf.gz"
    threads: 2
    log:
        "logs/SelectVariants/{case_id}.log"
    shell:
        """
        gatk SelectVariants \
        -V {input.vcf} \
        --select-type-to-include SNP \
        -O {output.vcf}
        """

rule SelectVariants_PASS_SNP_WES_Interval:
    input:
        vcf = "WXS_somatic_mutation/SelectVariants/PASS_SNP/{case_id}.vcf.gz",
        target_interval = config["ref"]["targetInterval"]
    output:
        vcf = "WXS_somatic_mutation/SelectVariants/PASS_SNP_WES_Interval/{case_id}.vcf.gz"
    threads: 2
    log:
        "logs/SelectVariants/{case_id}.log"
    shell:
        """
        gatk SelectVariants \
        -V {input.vcf} \
        -L {input.target_interval} \
        -O {output.vcf}
        """

rule SelectVariants_PASS_SNP_WES_Interval_exon:
    input:
        vcf = "WXS_somatic_mutation/SelectVariants/PASS_SNP_WES_Interval/{case_id}.vcf.gz",
        exon_interval = config["ref"]["exonInterval"]
    output:
        vcf = "WXS_somatic_mutation/SelectVariants/PASS_SNP_WES_Interval_exon/{case_id}.vcf.gz"
    threads: 2
    log:
        "logs/SelectVariants/{case_id}.log"
    shell:
        """
        gatk SelectVariants \
        -V {input.vcf} \
        -L {input.exon_interval} \
        -O {output.vcf}
        """

rule VariantsToTable_PASS_SNP_WES_Interval_exon:
    input:
        vcf = "WXS_somatic_mutation/SelectVariants/PASS_SNP_WES_Interval_exon/{case_id}.vcf.gz"
    output:
        table = "WXS_somatic_mutation/VariantsToTable/PASS_SNP_WES_Interval_exon/{case_id}.table"
    threads: 2
    shell:
        """
        gatk VariantsToTable \
        -V {input.vcf} \
        -F CHROM -F POS -F REF -F ALT \
        -O {output.table}
        """