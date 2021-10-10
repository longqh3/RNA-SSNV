#!/bin/bash
#PBS -N Megacolon_WGS_preparation
#PBS -l nodes=pan02:ppn=70
#PBS -l mem=700G
#PBS -q default
#PBS -o /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/logs/qsub/Megacolon_WGS_preparation_output.log
#PBS -e /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/logs/qsub/Megacolon_WGS_preparation_error.log

# 一般一个原始测序数据的bam文件
snakemake --cores 80 \
-s /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/rules/WGS-Preparation-tsv-Snakefile.smk \
--configfile /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/configs/Megacolon_WGS_preparation_config.yaml