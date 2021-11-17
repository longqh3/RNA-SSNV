# WGS相关
# 从原始WGS fastq.gz文件起始，到可用以突变检测的WGS bam文件为终止
# 应用BWA完成alignment和mapping操作，应用GATK相关工具完成对RNA reads的前处理事宜

"""
MarkDuplicate过程中出现问题"Mapped mate should have mate reference name"
https://github.com/cbrueffer/tophat-recondition/issues/1
"""
import pandas as pd
import os

from snakemake.utils import min_version

min_version("3.2")

workdir: config["workDir"]+"/"+config["diseaseType"]

# Prepare corresponding sample info dataframe
samples = pd.read_table(config["samples"]).set_index("file_name", drop=False)
# remove duplicate case fastq files
# samples.drop_duplicates(subset=["aliquots_id"], keep='first', inplace=True)

# Prepare lists to iterate over
aliquots_ids = set(samples["aliquots_id"])

rule all:
    input:
         expand("apply_BQSR/{aliquots_id}.bam", aliquots_id = aliquots_ids)

# time unknown
# cores 10
rule bwa:
    input:
        fastqs = lambda wildcards : [os.path.join(config["sampleDir"], samples["file_id"][file_name], samples["file_name"][file_name])
                                for file_name in samples.loc[samples["aliquots_id"]==wildcards.aliquots_id, ].index],
        # fq1 = lambda wildcards: os.path.join(config["sampleDir"], samples["Sample"][wildcards.sample_id]+"_"+samples["Library"][wildcards.sample_id], samples["Sample"][wildcards.sample_id]+"_"+samples["Library"][wildcards.sample_id]+"_1.clean.fq.gz"),
        # fq2 = lambda wildcards: os.path.join(config["sampleDir"], samples["Sample"][wildcards.sample_id]+"_"+samples["Library"][wildcards.sample_id], samples["Sample"][wildcards.sample_id]+"_"+samples["Library"][wildcards.sample_id]+"_2.clean.fq.gz"),
        ref = config["ref"]["genome"]
    output:
        bam = "bwa/{aliquots_id}.bam"
    params:
        ID = "{aliquots_id}",
        LB = "{aliquots_id}",
        SM = "{aliquots_id}",
        threads = 10
    threads: 10
    log:
        "logs/bwa/{aliquots_id}.log"
    shell:
        """
        bwa mem -t {params.threads} -R '@RG\\tID:{params.ID}\\tPL:illumina\\tLB:{params.LB}\\tSM:{params.SM}' {input.ref} {input.fastqs} | samtools view -S -b - > {output.bam}
        """

# time 6h+
# cores 4
rule sort_sam:
    input:
        bam = "bwa/{aliquots_id}.bam"
    output:
        bam = "sorted_sam/{aliquots_id}.bam"
    params:
        tmp = config["tmpDir"]
    threads: 4
    log:
        "logs/sort_sam/{aliquots_id}.log"
    shell:
        """
        gatk SortSam \
        -I {input.bam} \
        -O {output.bam} \
        --SORT_ORDER coordinate \
        --TMP_DIR {params.tmp}
        """

# time 30min+
# cores 10
rule mark_duplicates:
    input:
        bam = "sorted_sam/{aliquots_id}.bam",
        ref = config["ref"]["genome"]
    output:
        bam = "marked_duplicates/{aliquots_id}.bam",
        metric = "marked_duplicates/{aliquots_id}.txt"
    threads: 4
    log:
        "logs/mark_duplicates/{aliquots_id}.log"
    shell:
        """
        gatk MarkDuplicates \
        -I {input.bam} \
        -O {output.bam} \
        -M {output.metric} \
        -R {input.ref} \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT 
        """

# time 10min+
# cores 4
rule base_reclibrate:
    input:
        bam = "marked_duplicates/{aliquots_id}.bam",
        ref = config["ref"]["genome"],
        kgSNP = config["ref"]["variant"]["kgSNP"],
        kgINDEL = config["ref"]["variant"]["kgINDEL"]
    output:
        table = "base_reclibrate/{aliquots_id}.table"
    threads: 4
    log:
        "logs/base_reclibrate/{aliquots_id}.log"
    shell:
        """
        gatk BaseRecalibrator \
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
        bam = "marked_duplicates/{aliquots_id}.bam",
        table = "base_reclibrate/{aliquots_id}.table"
    output:
        bam = protected("apply_BQSR/{aliquots_id}.bam")
    threads: 4
    log:
        "logs/apply_BQSR/{aliquots_id}.log"
    shell:
        """
        gatk ApplyBQSR \
        -R {input.ref} \
        -I {input.bam} \
        --bqsr-recal-file {input.table} \
        -O {output.bam}
        """