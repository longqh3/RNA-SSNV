# RNA相关
# 从原始RNA fastq.gz文件起始，到可用以突变检测的RNA bam文件为终止
# 应用STAR完成alignment和mapping操作，应用GATK相关工具完成对RNA reads的前处理事宜

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
         expand("STAR-2-Pass/{aliquots_id}.bam", aliquots_id = aliquots_ids)

# specific for STAR-2.6.0c
rule STAR_2_pass:
    input:
        fastqs = lambda wildcards : [os.path.join(config["sampleDir"], samples["file_id"][file_name], samples["file_name"][file_name])
                                        for file_name in samples.loc[samples["aliquots_id"]==wildcards.aliquots_id, ].index]
        # fastq1 = lambda wildcards: os.path.join(config["sampleDir"], wildcards.aliquots_id+"m", wildcards.aliquots_id+"_1.fq.gz"),
        # fastq2 = lambda wildcards: os.path.join(config["sampleDir"], wildcards.aliquots_id+"m", wildcards.aliquots_id+"_2.fq.gz")
    params:
        STAR_genome_dir = config["ref"]["genome_dir"],
        RG_ID = "{aliquots_id}",
        RG_SM = "{aliquots_id}",
        RG_PL = "ILLUMINA",
        output_prefix = "STAR-2-Pass/{aliquots_id}-"
    output:
        bam = "STAR-2-Pass/{aliquots_id}.bam"
    threads: 20
    shell:
        """
        STAR \
        --readFilesIn {input.fastqs} \
        --outSAMattrRGline ID:{params.RG_ID} SM:{params.RG_SM} PL:{params.RG_PL} \
        --alignIntronMax 1000000 \
        --alignIntronMin 20 \
        --alignMatesGapMax 1000000 \
        --alignSJDBoverhangMin 1 \
        --alignSJoverhangMin 8 \
        --alignSoftClipAtReferenceEnds Yes \
        --chimJunctionOverhangMin 15 \
        --chimMainSegmentMultNmax 1 \
        --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
        --chimSegmentMin 15 \
        --genomeDir {params.STAR_genome_dir} \
        --genomeLoad NoSharedMemory \
        --limitSjdbInsertNsj 1200000 \
        --outFileNamePrefix {params.output_prefix} \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.1 \
        --outFilterMultimapNmax 20 \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterType BySJout \
        --outSAMattributes NH HI AS nM NM ch \
        --outSAMstrandField intronMotif \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --quantMode TranscriptomeSAM GeneCounts \
        --readFilesCommand zcat \
        --runThreadN 20 \
        --twopassMode Basic && mv {params.output_prefix}Aligned.sortedByCoord.out.bam {output.bam}
        """