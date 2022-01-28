import pandas as pd
import os

from snakemake.utils import min_version

min_version("3.2")

workdir: config["workDir"]+"/"+config["cancerType"]+"/"+"RNA"

# Prepare corresponding sample info
samples = pd.read_table(config["samples"]).set_index("aliquots_id", drop=False)

# To prevent starnge problem happen (case_id like 12327882 format shall be parsed into string format)
samples['case_id'] = samples['case_id'].astype(str)

tumor_samples = samples.loc[samples["sample_type"]=="Primary Tumor", ]
normal_samples = samples.loc[samples["sample_type"].isin(["Solid Tissue Normal", "Blood Derived Normal"]), ]
best_normal_samples = samples.loc[samples["sample_type"]=="Best Normal", ].set_index("case_id", drop=False)

# Prepare lists to iterate over
RNA_primary_tumor_aliquots_id = tumor_samples.index.values
case_ids = set(samples["case_id"])

rule all:
    input:
         # "RNA_germline_mutation/KGG_annotate/"+config["sampleType"]+"_SNP.flt.xlsx",
         expand("apply_BQSR/{aliquots_id}.bam", aliquots_id = RNA_primary_tumor_aliquots_id),
         expand("RNA_somatic_mutation/FilterMutectCalls_new/{case_id}.vcf.gz", case_id = case_ids),
         # SelectVariants_SNP
         expand("RNA_somatic_mutation/SelectVariants_new/SNP/{case_id}.vcf.gz", case_id = case_ids),
         # Funcotator_SNP
         expand("RNA_somatic_mutation/Funcotator_new/SNP/{case_id}.maf", case_id = case_ids),
         # SelectVariants_SNP_WES_interval_exon
         expand("RNA_somatic_mutation/SelectVariants_new/SNP_WES_Interval_exon/{case_id}.vcf.gz", case_id = case_ids),
         # VariantsToTable_SNP_WES_interval_exon
         expand("RNA_somatic_mutation/VariantsToTable/SNP_WES_Interval_exon/{case_id}.table", case_id = case_ids)

# time 30min+
# cores 10
rule mark_duplicates:
    input:
        bam = lambda wildcards: os.path.join(config["sampleDir"],tumor_samples["file_id"][wildcards.aliquots_id],tumor_samples["file_name"][wildcards.aliquots_id]),
        ref = config["ref"]["genome"]
    output:
        bam = temp("marked_duplicates/{aliquots_id}.bam"),
        metric = temp("marked_duplicates/{aliquots_id}_marked_dup_metrics.txt")
    params:
        gatk = config["software"]['gatk']
    threads: 8
    log:
        "logs/mark_duplicates/{aliquots_id}.log"
    shell:
        """
        {params.gatk} MarkDuplicates \
        -I {input.bam} \
        -O {output.bam} \
        -M {output.metric} \
        -R {input.ref} \
        --VALIDATION_STRINGENCY SILENT \
        --CREATE_INDEX true
        """

# time 6h+
# cores 4
rule split_N_cigar_reads:
    input:
        bam = "marked_duplicates/{aliquots_id}.bam",
        ref = config["ref"]["genome"]
    output:
        bam = temp("splited_n_cigar_reads/{aliquots_id}.bam")
    params:
        gatk = config["software"]['gatk']
    threads: 4
    log:
        "logs/split_N_cigar_reads/{aliquots_id}.log"
    shell:
        """
        {params.gatk} SplitNCigarReads \
        -I {input.bam} \
        -O {output.bam} \
        -R {input.ref}
        """

# time 10min+
# cores 4
rule base_reclibrate:
    input:
        bam = "splited_n_cigar_reads/{aliquots_id}.bam",
        ref = config["ref"]["genome"],
        kgSNP = config["ref"]["variant"]["kgSNP"],
        kgINDEL = config["ref"]["variant"]["kgINDEL"]
    output:
        table = temp("base_reclibrate/{aliquots_id}.table")
    params:
        gatk = config["software"]['gatk']
    threads: 4
    log:
        "logs/base_reclibrate/{aliquots_id}.log"
    shell:
        """
        {params.gatk} BaseRecalibrator \
        -I {input.bam} \
        -R {input.ref} \
        --known-sites {input.kgSNP} \
        --known-sites {input.kgINDEL} \
        -O {output.table}
        """

# time 30min+
# cores 4
rule apply_BQSR:
    input:
        ref = config["ref"]["genome"],
        bam = "splited_n_cigar_reads/{aliquots_id}.bam",
        table = "base_reclibrate/{aliquots_id}.table"
    output:
        bam = protected("apply_BQSR/{aliquots_id}.bam")
    threads: 4
    params:
        gatk = config["software"]['gatk']
    log:
        "logs/apply_BQSR/{aliquots_id}.log"
    shell:
        """
        {params.gatk} ApplyBQSR \
        -R {input.ref} \
        -I {input.bam} \
        --bqsr-recal-file {input.table} \
        -O {output.bam}
        """

# takes 22h
# cores 2
# deleted          --panel-of-normals {input.PoN} \ and         PoN = config["ref"]["variant"]["PoN"],
# because that was for consistency of LUAD WXS and RNA-seq processing and apparently RNA-seq cannot find proper PoN.
# deleted         --enable-all-annotations true \
# because FilterMutectCalls within 4.2.0.0 cannot accept some annotations
# 2021.6.6 recover the usage of "--enable-all-annotations true" in order to maximize features
# 2021.6.8 delete "         --bam-output {output.bam} \" in order to increase performance
rule Mutect2:
    input:
        ref = config["ref"]["genome"],
        gnomad = config["ref"]["variant"]["gnomad"],
        tumor_bams = lambda wildcards : expand("apply_BQSR/{tumor_aliquots_id}.bam",
                                        tumor_aliquots_id = tumor_samples.loc[tumor_samples["case_id"]==wildcards.case_id, ].index),
        normal_bams = lambda wildcards : [os.path.join(config["normalSampleDir"],normal_samples["file_id"][normal_aliquots_id],normal_samples["file_name"][normal_aliquots_id])
                                        for normal_aliquots_id in normal_samples.loc[normal_samples["case_id"]==wildcards.case_id, ].index]
    output:
        vcf = protected("RNA_somatic_mutation/Mutect2_new/{case_id}.vcf.gz"),
        f1r2 = protected("RNA_somatic_mutation/Mutect2_new/{case_id}.f1r2.tar.gz"),
    params:
        gatk = config["software"]['gatk_new'],
        tumor_bams = lambda wildcards : expand("-I apply_BQSR/{tumor_aliquots_id}.bam",
                                        tumor_aliquots_id = tumor_samples.loc[tumor_samples["case_id"]==wildcards.case_id, ].index),
        normal_bams = lambda wildcards : ["-I "+os.path.join(config["normalSampleDir"],normal_samples["file_id"][normal_aliquots_id],normal_samples["file_name"][normal_aliquots_id])
                                        for normal_aliquots_id in normal_samples.loc[normal_samples["case_id"]==wildcards.case_id, ].index],
        normal_sample_name = lambda wildcards : expand("-normal {normal_sample_names}",
                                                       normal_sample_names = normal_samples.loc[normal_samples["case_id"]==wildcards.case_id, "aliquots_id"]),
        threads = 2
    threads: 2
    log:
        "logs/Mutect2_new/{case_id}.log"
    shell:
        """
         {params.gatk} Mutect2 \
         -R {input.ref} \
         {params.tumor_bams} \
         {params.normal_bams} \
         {params.normal_sample_name} \
         --germline-resource {input.gnomad} \
         --native-pair-hmm-threads {params.threads} \
         --f1r2-tar-gz {output.f1r2} \
         --dont-use-soft-clipped-bases true \
         --enable-all-annotations true \
         -O {output.vcf}
        """

# extremely fuzzy task
# 首先需要执行GetPileupSummaries来获取best normal测序数据所对应的PileupSummary信息
# 接着（多次）执行GetPileupSummaries来获取所有tumor测序数据所对应的PileupSummary信息，同时应用CalculateContamination来获取所有tumor-best normal对的contamination.table、segments.table信息
# 并将相应contamination-table、tumor-segmentation信息纳入变量中保存
# 而后执行LearnReadOrientationModel来获取read_orientation_model信息
# 最后应用上述所有信息，共同完成单case的FilterMutectCalls操作，得到经过过滤的vcf文件信息
# takes
# cores 2
rule FilterMutectCalls_combined:
    input:
        ref = config["ref"]["genome"],
        tumor_bams = lambda wildcards : expand("apply_BQSR/{tumor_aliquots_id}.bam",
                                        tumor_aliquots_id = tumor_samples.loc[tumor_samples["case_id"]==wildcards.case_id, ].index),
        best_normal_bam = lambda wildcards : os.path.join(config["normalSampleDir"],best_normal_samples["file_id"][wildcards.case_id],best_normal_samples["file_name"][wildcards.case_id]),
        vcf = "RNA_somatic_mutation/Mutect2_new/{case_id}.vcf.gz",
        f1r2 = "RNA_somatic_mutation/Mutect2_new/{case_id}.f1r2.tar.gz",
        gnomad_SNP = config["ref"]["variant"]["gnomad_SNP"]
    output:
        best_normal_pileups_table = "RNA_somatic_mutation/GetPileupSummaries_new/{case_id}-best-normal-pileups.table",
        read_orientation_model = "RNA_somatic_mutation/Mutect2_new/{case_id}.read-orientation-model.tar.gz",
        vcf = protected("RNA_somatic_mutation/FilterMutectCalls_new/{case_id}.vcf.gz")
    params:
        gatk = config["software"]['gatk'],
        tumor_bams = lambda wildcards : ["apply_BQSR/%s.bam" % (tumor_aliquots_id) for tumor_aliquots_id in tumor_samples.loc[tumor_samples["case_id"]==wildcards.case_id, ].index]
    threads: 2
    log:
        "logs/FilterMutectCalls_combined/{case_id}.log"
    run:
        commands = []
        # For Best Normal GetPileupSummaries
        normal_GetPileupSummaries_command = f"{params.gatk} GetPileupSummaries -I {input.best_normal_bam} -L {input.gnomad_SNP} -V {input.gnomad_SNP} -O {output.best_normal_pileups_table}"
        commands.append(normal_GetPileupSummaries_command)
        # For --contamination-table info
        contamination_table = ""
        # For --tumor-segmentation info
        tumor_segmentation = ""
        for tumor_bam in list(input.tumor_bams):
            print("\n"+tumor_bam)
            # For Tumor GetPileupSummaries
            GetPileupSummaries_command = f"{params.gatk} GetPileupSummaries -I {tumor_bam} -L {input.gnomad_SNP} -V {input.gnomad_SNP} -O {tumor_bam}-pileups.table"
            print(GetPileupSummaries_command)
            commands.append(GetPileupSummaries_command)
            # For CalculateContaminations
            CalculateContaminations_command = f"{params.gatk} CalculateContamination -I {tumor_bam}-pileups.table -matched {output.best_normal_pileups_table} -O {tumor_bam}-contamination.table --tumor-segmentation {tumor_bam}-segments.table"
            print(CalculateContaminations_command)
            commands.append(CalculateContaminations_command)
            # Add corresponding table info
            contamination_table = contamination_table + f" --contamination-table {tumor_bam}-contamination.table"
            tumor_segmentation = tumor_segmentation + f" --tumor-segmentation {tumor_bam}-segments.table"
        print("--------------------------------------------------")
        # For LearnReadOrientationModel
        LearnReadOrientationModel_command = f"{params.gatk} LearnReadOrientationModel -I {input.f1r2} -O {output.read_orientation_model}"
        print(LearnReadOrientationModel_command)
        commands.append(LearnReadOrientationModel_command)
        # For FilterMutectCalls
        FilterMutectCalls_command = f"{params.gatk} FilterMutectCalls -R {input.ref} -V {input.vcf} {contamination_table} {tumor_segmentation} --orientation-bias-artifact-priors {output.read_orientation_model} -O {output.vcf}"
        print(FilterMutectCalls_command)
        commands.append(FilterMutectCalls_command)

        for c in commands:
            shell(c)

rule SelectVariants_SNP:
    input:
        vcf = "RNA_somatic_mutation/FilterMutectCalls_new/{case_id}.vcf.gz"
    output:
        vcf = "RNA_somatic_mutation/SelectVariants_new/SNP/{case_id}.vcf.gz"
    params:
        gatk = config["software"]['gatk']
    threads: 2
    shell:
        """
        {params.gatk} SelectVariants \
        --select-type-to-include SNP \
        -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY \
        -V {input.vcf} \
        -O {output.vcf}
        """

rule Funcotator_SNP:
    input:
        vcf = "RNA_somatic_mutation/SelectVariants_new/SNP/{case_id}.vcf.gz",
        ref = config["ref"]["genome"]
    output:
        maf = "RNA_somatic_mutation/Funcotator_new/SNP/{case_id}.maf"
    params:
        data_source = config["ref"]["annotation"]["gatk_funcotator"],
        case_barcode = "{case_id}",
        gatk = config["software"]['gatk']
    threads: 2
    shell:
        """
        {params.gatk} Funcotator \
        -R {input.ref} \
        -V {input.vcf} \
        -O {output.maf} \
        --output-file-format MAF \
        --data-sources-path {params.data_source} \
        --annotation-default tumor_barcode:{params.case_barcode} \
        --annotation-default normal_barcode:{params.case_barcode} \
        --ref-version hg38
        """

rule SelectVariants_SNP_WES_interval:
    input:
        vcf = "RNA_somatic_mutation/SelectVariants_new/SNP/{case_id}.vcf.gz",
        target_interval = config["ref"]["targetInterval"]
    output:
        vcf = "RNA_somatic_mutation/SelectVariants_new/SNP_WES_Interval/{case_id}.vcf.gz"
    params:
        gatk = config["software"]['gatk']
    threads: 2
    shell:
        """
        {params.gatk} SelectVariants \
        -L {input.target_interval} \
        -V {input.vcf} \
        -O {output.vcf}
        """

rule SelectVariants_SNP_WES_interval_exon:
    input:
        vcf = "RNA_somatic_mutation/SelectVariants_new/SNP_WES_Interval/{case_id}.vcf.gz",
        exon_interval = config["ref"]["exonInterval"]
    output:
        vcf = "RNA_somatic_mutation/SelectVariants_new/SNP_WES_Interval_exon/{case_id}.vcf.gz"
    params:
        gatk = config["software"]['gatk']
    threads: 2
    shell:
        """
        {params.gatk} SelectVariants \
        -L {input.exon_interval} \
        -V {input.vcf} \
        -O {output.vcf}
        """

rule VariantsToTable_SNP_WES_interval_exon:
    input:
        vcf = "RNA_somatic_mutation/SelectVariants_new/SNP_WES_Interval_exon/{case_id}.vcf.gz"
    output:
        table = "RNA_somatic_mutation/VariantsToTable/SNP_WES_Interval_exon/{case_id}.table"
    params:
        gatk = config["software"]['gatk']
    threads: 1
    shell:
         """
         {params.gatk} VariantsToTable \
         -V {input.vcf} \
         -F CHROM -F POS -F REF -F ALT \
         --show-filtered true \
         -O {output.table}
         """