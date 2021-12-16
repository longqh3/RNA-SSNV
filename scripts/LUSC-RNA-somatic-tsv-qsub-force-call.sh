#!/bin/bash
#PBS -N LUSC_RNA_force_call@pan02_1day
#PBS -l nodes=pan02:ppn=80
#PBS -l mem=400G
#PBS -q default
#PBS -o /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/logs/qsub/LUSC_RNA_Somatic_force_call_output.log
#PBS -e /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/logs/qsub/LUSC_RNA_Somatic_force_call_error.log

# 一般一个原始测序数据的bam文件
snakemake --cores 80 \
-s /home/lqh/Codes/Python/RNA-SSNV/rules/RNA-Somatic-tsv-Snakefile-force-call.smk \
--configfile /home/lqh/Codes/Python/RNA-SSNV/configs/LUSC_RNA_Somatic_config_force_call.yaml \
--rerun-incomplete