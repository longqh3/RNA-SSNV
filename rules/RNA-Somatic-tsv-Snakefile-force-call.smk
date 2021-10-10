# 核心代码
# 根据指定Interval信息（未被Mutect2 call出来的突变）来检查原因

"""
RNAseq short variant discovery (SNPs + Indels)
https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-
GATK Best Practices — step4 转录组SNP/INDEL（RNAseq SNPs + Indels）
https://www.jianshu.com/p/a6852891d1f4
最新版针对RNA-seq数据的GATK找变异流程
https://cloud.tencent.com/developer/article/1536221
Variant calling from RNA-seq data using STAR, Picard and GATK tools
http://rpubs.com/bhagirathi_das/494122
MarkDuplicate过程中出现问题"Mapped mate should have mate reference name"
https://github.com/cbrueffer/tophat-recondition/issues/1
"""
import pandas as pd
import os

from snakemake.utils import min_version

min_version("3.2")

workdir: config["workDir"]+"/"+config["cancerType"]+"/"+"RNA"

# Prepare corresponding sample info
samples = pd.read_table(config["samples"]).set_index("aliquots_id", drop=False)
tumor_samples = samples.loc[samples["sample_type"]=="Primary Tumor", ]
normal_samples = samples.loc[samples["sample_type"].isin(["Solid Tissue Normal", "Blood Derived Normal"]), ]
best_normal_samples = samples.loc[samples["sample_type"]=="Best Normal", ].set_index("case_id", drop=False)

# Prepare lists to iterate over
RNA_primary_tumor_aliquots_id = tumor_samples.index.values
case_ids = set(samples["case_id"])

# use force-call case list to exclude
force_call_case_ids = set(pd.read_table(config["force_call_samples"])["0"])
case_ids = case_ids.intersection(force_call_case_ids)

rule all:
    input:
         # expand("RNA_somatic_mutation/MAFToVCF/gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_check/{case_id}.vcf.idx", case_id = case_ids),
         expand("RNA_somatic_mutation/Funcotator_new_force_call/{case_id}.maf", case_id = case_ids)

rule SortVcf:
    input:
        force_call_vcf = "RNA_somatic_mutation/MAFToVCF/gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_check/{case_id}.vcf"
    output:
        force_call_vcf_sorted = "RNA_somatic_mutation/MAFToVCF/gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_check_sorted/{case_id}.vcf"
    shell:
        """
         gatk SortVcf \
         -I {input.force_call_vcf} \
         -O {output.force_call_vcf_sorted}
        """

# rule IndexFeatureFile:
#     input:
#         force_call_vcf = "RNA_somatic_mutation/MAFToVCF/gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_check/{case_id}.vcf"
#     output:
#         force_call_vcf_index = "RNA_somatic_mutation/MAFToVCF/gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_check/{case_id}.vcf.idx"
#     shell:
#         """
#          gatk IndexFeatureFile \
#          -I {input.force_call_vcf}
#         """

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
                                        for normal_aliquots_id in normal_samples.loc[normal_samples["case_id"]==wildcards.case_id, ].index],
        force_call_vcf_sorted = "RNA_somatic_mutation/MAFToVCF/gdc_validate_all_info_DNA_only_RNA_missing_RNA_info_Mutect2_check_sorted/{case_id}.vcf"
    output:
        vcf = protected("RNA_somatic_mutation/Mutect2_new_force_call/{case_id}.vcf.gz"),
        f1r2 = protected("RNA_somatic_mutation/Mutect2_new_force_call/{case_id}.f1r2.tar.gz"),
        bam_out = protected("RNA_somatic_mutation/Mutect2_new_force_call/{case_id}.bamout.bam")
    params:
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
         /home/lqh/software/GATK-4.2.0.0/gatk Mutect2 \
         -R {input.ref} \
         {params.tumor_bams} \
         {params.normal_bams} \
         {params.normal_sample_name} \
         --force-call-filtered-alleles true \
         -L {input.force_call_vcf_sorted} \
         --alleles {input.force_call_vcf_sorted} \
         --germline-resource {input.gnomad} \
         --native-pair-hmm-threads {params.threads} \
         --f1r2-tar-gz {output.f1r2} \
         --dont-use-soft-clipped-bases true \
         --enable-all-annotations true \
         -bamout {output.bam_out} \
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
        vcf = "RNA_somatic_mutation/Mutect2_new_force_call/{case_id}.vcf.gz",
        f1r2 = "RNA_somatic_mutation/Mutect2_new_force_call/{case_id}.f1r2.tar.gz",
        gnomad_SNP = config["ref"]["variant"]["gnomad_SNP"]
    output:
        best_normal_pileups_table = "RNA_somatic_mutation/GetPileupSummaries_new_force_call/{case_id}-best-normal-pileups.table",
        read_orientation_model = "RNA_somatic_mutation/Mutect2_new_force_call/{case_id}.read-orientation-model.tar.gz",
        vcf = protected("RNA_somatic_mutation/FilterMutectCalls_new_force_call/{case_id}.vcf.gz")
    params:
        tumor_bams = lambda wildcards : ["apply_BQSR/%s.bam" % (tumor_aliquots_id) for tumor_aliquots_id in tumor_samples.loc[tumor_samples["case_id"]==wildcards.case_id, ].index]
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
            # For Tumor GetPileupSummaries
            GetPileupSummaries_command = f"gatk GetPileupSummaries -I {tumor_bam} -L {input.gnomad_SNP} -V {input.gnomad_SNP} -O {tumor_bam}-pileups.table"
            print(GetPileupSummaries_command)
            commands.append(GetPileupSummaries_command)
            # For CalculateContaminations
            CalculateContaminations_command = f"gatk CalculateContamination -I {tumor_bam}-pileups.table -matched {output.best_normal_pileups_table} -O {tumor_bam}-contamination.table --tumor-segmentation {tumor_bam}-segments.table"
            print(CalculateContaminations_command)
            commands.append(CalculateContaminations_command)
            # Add corresponding table info
            contamination_table = contamination_table + f" --contamination-table {tumor_bam}-contamination.table"
            tumor_segmentation = tumor_segmentation + f" --tumor-segmentation {tumor_bam}-segments.table"
        print("--------------------------------------------------")
        # For LearnReadOrientationModel
        LearnReadOrientationModel_command = f"gatk LearnReadOrientationModel -I {input.f1r2} -O {output.read_orientation_model}"
        print(LearnReadOrientationModel_command)
        commands.append(LearnReadOrientationModel_command)
        # For FilterMutectCalls
        FilterMutectCalls_command = f"gatk FilterMutectCalls -R {input.ref} -V {input.vcf} {contamination_table} {tumor_segmentation} --orientation-bias-artifact-priors {output.read_orientation_model} -O {output.vcf}"
        print(FilterMutectCalls_command)
        commands.append(FilterMutectCalls_command)

        for c in commands:
            shell(c)

rule Funcotator:
    input:
        vcf = "RNA_somatic_mutation/FilterMutectCalls_new_force_call/{case_id}.vcf.gz",
        ref = config["ref"]["genome"]
    output:
        maf = "RNA_somatic_mutation/Funcotator_new_force_call/{case_id}.maf"
    params:
        data_source = config["ref"]["annotation"]["gatk_funcotator"],
        case_barcode = "{case_id}"
    threads: 2
    shell:
        """
        gatk Funcotator \
        -R {input.ref} \
        -V {input.vcf} \
        -O {output.maf} \
        --output-file-format MAF \
        --data-sources-path {params.data_source} \
        --annotation-default tumor_barcode:{params.case_barcode} \
        --annotation-default normal_barcode:{params.case_barcode} \
        --ref-version hg38
        """