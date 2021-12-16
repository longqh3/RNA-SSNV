#!/bin/bash
#PBS -N GBM_RNA_force_call@pan01_1day
#PBS -l nodes=pan01:ppn=60
#PBS -l mem=500G
#PBS -q default
#PBS -o /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/logs/qsub/GBM_RNA_Somatic_force_call_output.log
#PBS -e /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/logs/qsub/GBM_RNA_Somatic_force_call_error.log

# 一般一个原始测序数据的bam文件
snakemake --cores 60 \
-s /home/lqh/Codes/Python/RNA-SSNV/rules/RNA-Somatic-tsv-Snakefile-force-call.smk \
--configfile /home/lqh/Codes/Python/RNA-SSNV/configs/GBM_RNA_Somatic_config_force_call.yaml \
--rerun-incomplete