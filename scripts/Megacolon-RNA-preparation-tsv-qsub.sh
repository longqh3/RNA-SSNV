#!/bin/bash
#PBS -N Megacolon_RNA_prepare
#PBS -l nodes=pan01:ppn=60
#PBS -l mem=600G
#PBS -q default
#PBS -o /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/logs/qsub/Megacolon_RNA_preparation_output.log
#PBS -e /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/logs/qsub/Megacolon_RNA_preparation_error.log

# 一般一个原始测序数据的bam文件
snakemake --cores 64 \
-s /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/rules/RNA-Preparation-tsv-Snakefile.smk \
--configfile /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/configs/Megacolon_RNA_preparation_config.yaml