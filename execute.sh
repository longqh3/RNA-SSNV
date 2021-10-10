python /home/lqh/Codes/Python/RNA-SSNV/model_construct.py \
--REDIportal /home/lqh/resources/database/RNA_edit/REDIportal/REDIportal_main_table.hg38.bed \
--DARNED /home/lqh/resources/database/RNA_edit/DARNED_hg19_to_bed_to_hg38_rm_alt.bed \
--raw_RNA_mutations /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/LUAD/RNA/RNA_somatic_mutation/VcfAssembly_new/SNP_WES_Interval_exon.txt \
--GDC_mutations /home/lqh/Codes/Data/TCGA_maf_files/TCGA-LUAD \
--WES_mutations /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/LUAD/WXS/WXS_somatic_mutation/VariantsToTable/PASS_SNP_WES_Interval_exon.table \
--model_folder_path /home/lqh/Codes/Python/RNA-SSNV/model