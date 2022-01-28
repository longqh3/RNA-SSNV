#!/bin/bash
#PBS -N RNA_force_call@pan02_1day
#PBS -l nodes=pan02:ppn=80
#PBS -l mem=800G
#PBS -q default
#PBS -o /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/logs/qsub/RNA_Somatic_force_call_output.log
#PBS -e /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/logs/qsub/RNA_Somatic_force_call_error.log

# 一般一个原始测序数据的bam文件
snakemake --cores 80 \
-ns rules/RNA-Somatic-tsv-Snakefile-force-call.smk \
--configfile configs/project_force_call_config.yaml \
--rerun-incomplete