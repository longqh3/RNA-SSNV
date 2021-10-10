#!/bin/bash
#PBS -N LUAD_RNA_Somatic@pan02
#PBS -l nodes=pan02:ppn=80
#PBS -l mem=800G
#PBS -q default
#PBS -o /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/logs/qsub/LUAD_RNA_somatic_calling_output.log
#PBS -e /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/logs/qsub/LUAD_RNA_somatic_calling_error.log

# 一般一个原始测序数据的bam文件
snakemake --cores 120 \
-s /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/rules/RNA-Somatic-tsv-Snakefile-new-option.smk \
--configfile /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/configs/LUAD_RNA_Somatic_config.yaml \
--rerun-incomplete