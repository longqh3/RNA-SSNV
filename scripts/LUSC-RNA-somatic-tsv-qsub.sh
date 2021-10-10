#!/bin/bash
#PBS -N LUSC_RNA_Somatic@pan02_7days
#PBS -l nodes=pan02:ppn=1
#PBS -l mem=500G
#PBS -q default
#PBS -o /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/logs/qsub/LUSC_RNA_Somatic_validation_output.log
#PBS -e /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/logs/qsub/LUSC_RNA_Somatic_validation_error.log

# 一般一个原始测序数据的bam文件
snakemake --cores 60 \
-ns /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/rules/RNA-Somatic-tsv-Snakefile-new-option.smk \
--configfile /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/configs/LUSC_RNA_Somatic_config.yaml \
--rerun-incomplete