RNA-SSNV
=======================================================

The RNA-SSNV is a scalable and efficient analysis method for RNA somatic mutation detection from RNA-WES(tumor-normal) paired sequencing data which utilized Mutect2 as core-caller and multi-filtering strategy & Machine-learning based model to maximize precision & recall performance. It runs highly automated once snakemake and related configs & infos get configurated properly. It reports an aggregated mutation file (maf format) to facilitate downstream analysis and clinical decision. 

Pre-requirements
~~~~~~~~~~~~~~~~~

Download codes
----------------------

.. code:: sh
    
  git clone https://github.com/longqh3/RNA-SSNV

Softwares
----------

Beware, all the following softwares should be added into **environmental variables** (*~/.bashrc*, *~/.bash_profile*, etc)

- GATK-v4.1.6.0
  
  `Official download link <https://github.com/broadinstitute/gatk/releases/download/4.1.6.0/gatk-4.1.6.0.zip>`_ & `Back-up download link <http://link>`_ 

- Snakemake-v5.10.0

  .. code:: sh

    pip install snakemake==5.10.0

- Python packages: 
    
    - scikit-learn: v4.1


  .. code:: sh

    pip install -r requirements.txt


Modify config & table file
---------------------------

- *configs/project_config.yaml*: pipeline-related configurations, modify it accordingly. 
- *tables/project_RNA_somatic_calling_info.tsv*: project-related sequencing data information, modify under following instruction.

.. list-table:: Project RNA somatic calling  sample info
    :widths: auto
    :header-rows: 1
    :align: center

    * - file_name
      - file_id
      - aliquots_id
      - case_id
      - sample_type
    * - Name of bam file
      - Name of folder which contains bam file
      - ID of aliquot sequenced
      - ID of corresponding case
      - Type of aliquot's origin

Modify entry shell scripts
--------------------------

- *scripts/project_RNA_somatic-tsv-qsub.sh*: shell script as entry command for whole project, modify it accordingly (support PBS task management system).

Select interval files
---------------------

- *resources/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals_add_chr_to_hg38_rm_alt.bed*: bed-format interval file for paired-normal Whole Exome Sequence(WES) targets, canonical for TCGA projects. (change it with your own WES target interval file)
- *resources/GRCh38_GENCODE_v22_exon_rm_alt.bed*: bed-format interval file for GENCODE v22 exon regions. (change it with your own exon regions of interest)

  All interval files should be *bed* format and contain column names for "*chr*  *start* *end*". 

Run Pipeline
~~~~~~~~~~~~~~~

Once configurated correctly, our pipeline is ready to go. Please execute the following commands step by step, make sure everything works normally before moving forward. 

Call and annotate raw RNA somatic mutations
-----------------------------------------------

.. code:: sh
    
    # dry run to see if the mutation calling pipeline works
    snakemake --cores {num_of_cores} \
    -ns rules/RNA-Somatic-tsv-Snakefile.smk \
    --configfile configs/project_RNA_Somatic_config.yaml

    # run the pipeline
    snakemake --cores {num_of_cores} \
    -s rules/RNA-Somatic-tsv-Snakefile.smk \
    --configfile configs/project_RNA_Somatic_config.yaml

Beware, owing to the breakpoint-run feature of snakemake, our pipeline also supports taking any final files (listed below) as starting point. 

Prepare features for raw RNA somatic mutations
-----------------------------------------------

.. code:: sh

    # run feature-extraction codes
    python lib/own_data_vcf_info_retriver.py \
    --cancer_type BLCA \
    --RNA_calling_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/tables/info/BLCA_RNA_somatic_calling_info.tsv \
    --project_folder /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results \
    --exon_interval /home/lqh/resources/database/gencode/GRCh38_GENCODE_v22_exon_rm_alt.bed \
    --output_table_path /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/BLCA/RNA/RNA_somatic_mutation/VcfAssembly_new/SNP_WES_Interval_exon.txt \
    --num_threads 60

Predict reliable RNA somatic mutations
------------------------------------------

For the generated result, the records with *pred_label* being 1 should be considered as reliable RNA somatic mutations. 

.. code:: sh

    # run model predicting codes
    python /home/lqh/Codes/Python/RNA-SSNV/model_utilize.py \
    --REDIportal /home/lqh/resources/database/RNA_edit/REDIportal/REDIportal_main_table.hg38.bed \
    --DARNED /home/lqh/resources/database/RNA_edit/DARNED_hg19_to_bed_to_hg38_rm_alt.bed \
    --raw_RNA_mutations /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/GBM/RNA/RNA_somatic_mutation/VcfAssembly_new/SNP_WES_Interval_exon.txt \
    --model_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer.model \
    --one_hot_encoder_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer.one_hot_encoder \
    --training_columns_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer.training_data_col \
    --output_table_path /home/lqh/Codes/Python/RNA-SSNV/output/GBM.table

Pairwise analysis for DNA and RNA somatic mutations (only do it with DNA evidence)
----------------------------------------------------------------------------------------

Step 1: Generate RNA-omitted DNA mutations to force-call
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: sh

    python /home/lqh/Codes/Python/RNA-SSNV/model_analyze_with_DNA.py \
    --step 1 \
    --cancer_type BLCA \
    --DNA_info /home/lqh/Codes/Data/TCGA_maf_files/TCGA-BLCA \
    --RNA_info /home/lqh/Codes/Python/RNA-SSNV/output/BLCA.table \
    --WXS_target_interval /home/lqh/resources/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals_add_chr_to_hg38_rm_alt.bed \
    --exon_interval /home/lqh/resources/database/gencode/GRCh38_GENCODE_v22_exon_rm_alt.bed \
    --RNA_calling_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/tables/info/BLCA_RNA_somatic_calling_info.tsv \
    --RNA_bam_folder /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/BLCA/RNA/apply_BQSR \
    --Mutect2_target_detected_sites /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/BLCA/RNA/RNA_somatic_mutation/VariantsToTable/SNP_WES_Interval_exon.table \
    --project_folder /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results \
    --num_threads 40 \
    --output_file_path /home/lqh/Codes/Python/RNA-SSNV/output/BLCA_DNA_step_1.class

Step 1.1: Force calling all DNA only mutations and extract features
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Modify config file for force-calling process

- *configs/project_force_call_config.yaml*: pipeline-related configurations, modify it accordingly. 

Run commands sequencially.

.. code:: sh
    
    # dry run to see if the mutation calling pipeline works
    snakemake --cores {num_of_cores} \
    -ns rules/RNA-Somatic-tsv-Snakefile-force-call.smk \
    --configfile configs/project_RNA_Somatic_config_force_call.yaml \
    --rerun-incomplete

    # run formally
    snakemake --cores {num_of_cores} \
    -s rules/RNA-Somatic-tsv-Snakefile.smk \
    --configfile configs/project_RNA_Somatic_config.yaml

    # run feature extraction codes after force-calling
    python force_call_data_vcf_info_retriver.py \
    --cancer_type GBM \
    --RNA_calling_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/tables/info/GBM_RNA_somatic_calling_info.tsv \
    --project_folder /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results \
    --exon_interval /home/lqh/resources/database/gencode/GRCh38_GENCODE_v22_exon_rm_alt.bed \
    --output_table_path /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/GBM/RNA/RNA_somatic_mutation/VcfAssembly_new/Mutect2_force_call.txt \
    --num_threads 80


Step 2: Combine force-called results with RNA somatic mutations to finish RNA-DNA integrative analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: py

    python /home/lqh/Codes/Python/RNA-SSNV/model_analyze_with_DNA.py \
    --step 2 \
    --force_call_RNA_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/GBM/RNA/RNA_somatic_mutation/VcfAssembly_new/Mutect2_force_call.txt \
    --instance_path /home/lqh/Codes/Python/RNA-SSNV/output/GBM_DNA_step_1.class \
    --model_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer.model \
    --one_hot_encoder_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer.one_hot_encoder \
    --training_columns_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer.training_data_col \
    --output_file_path /home/lqh/Codes/Python/RNA-SSNV/output/GBM.final.table

Output folders & files
~~~~~~~~~~~~~~~~~~~~~~~~~~

The pipeline outputs several folders containing intermediate files and final project-level mutations annotation file (maf format). Here, we describe the `results/` folder's schema. 

Sequencing data pre-process
------------------------------

- *results/project_name/RNA/marked_duplicates*: temporary folder containing MarkDuplicates tool's output.
- *results/project_name/RNA/splited_n_cigar_reads*: temporary folder containing SplitNCigarReads tool's output.
- `results/project_name/RNA/base_reclibrate`: temporary folder containing BaseRecalibrate tool's output.
- *results/project_name/RNA/apply_BQSR*: permanent folder containing ApplyBQSR tool's output, **final** files (bam format) used to call RNA somatic mutations, **applicable** for other analysis.

Calling process - RNA somatic mutation
-----------------------------------------

- *results/project_name/RNA/RNA_somatic_mutation/Mutect2*: permanent folder containing Mutect2 tool's output. 
- *results/project_name/RNA/RNA_somatic_mutation/GetPileupSummaries*: permanent folder containing GetPileupSummaries tool's output (best normal sample's pileup summary info).
- *results/project_name/RNA/RNA_somatic_mutation/FilterMutectCalls*: permanent folder containing FilterMutectCalls tool's output, **final** files (vcf format) used to discriminate true RNA somatic mutations, applicable for other filtering strategy. 

Model prediction process - RNA somatic mutation
---------------------------------------------------------

- *results/project_name/RNA/RNA_somatic_mutation/Funcotator/SNP*: permanent folder containing Funcotator's annnotation info for raw RNA SNP calls. 
- *results/project_name/RNA/RNA_somatic_mutation/SelectVariants/SNP_WES_interval*: permanent folder containing raw RNA SNP calls subsetted via given WES target intervals. 
- *results/project_name/RNA/RNA_somatic_mutation/SelectVariants/SNP_WES_interval_exon*: permanent folder containing **final** raw RNA SNP calls subsetted via given WES target intervals and exon regions.

Pair-wise analysis with DNA process - RNA-DNA somatic mutation
-----------------------------------------------------------------------

- *results/project_name/RNA/RNA_somatic_mutation/VcfAssembly/SNP_WES_interval_exon*: permanent folder containing extracted features and other info per case. 
- *results/project_name/RNA/RNA_somatic_mutation/VcfAssembly/SNP_WES_interval_exon_positive.maf*: **final result** file for whole project - total project's Mutect2 calls marked as **positive** by our discriminant model and default threshold.

Pipeline explaination
~~~~~~~~~~~~~~~~~~~~~~~~~

Essential codes
------------------

- *rules/RNA_Somatic-tsv-Snakefile.smk*: snakemake-style codes to describe whole pipeline (modify at your own risk!!!). 
- *codes/vcf_info_retriver_tsv.py*: python codes to extract features (variant, genotype and annotation level) from different sources. 
- *codes/function_based_RNA_somatic_random_forest_prediction.py*: python codes to predict the probability of given Mutect2 calls being true RNA somatic mutations. 

Pre-trained models
----------------------

- *models/data_ormalization_model.model*: data normalization model which adapted to following model.
- *models/classic_random_forest_model.model*: random forest discriminant model trained using whole TCGA LUAD project data.

Resource files
------------------

- *resources/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals_add_chr_to_hg38_rm_alt.bed*: bed-format interval file for paired-normal Whole Exome Sequence(WES) targets. (canonical for TCGA projects)
- *resources/GRCh38_GENCODE_v22_exon_rm_alt.bed*: bed-format interval file for GENCODE v22 exon regions. 


P.S. Train your own discriminant model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Although we used 511 cases of TCGA LUAD RNA-WES paired data to train our discriminant model, other non-cancerous RNA somatic mutations or non-bulk RNA-Seq data may exhibit **different patterns of FP calls**. In that case, our model may not served as expected, and a customized model was required to be trained on your own. 

Data-preparation
--------------------

- Gold-standard TP mutations for given project (maf-format) with required columns: "Chromosome", "Start_Position", "Tumor_Allele2", "Tumor_Allele1", "Tumor_Sample_UUID"

Train customized model
-----------------------

- Using gold-standard TP mutations with their corresponding RNA somatic mutations to train customized model. The performance matrix for model training will be generated. 

  .. code:: sh
    
    # run feature-extraction codes
    python lib/own_data_vcf_info_retriver.py \
    --cancer_type BLCA \
    --RNA_calling_info /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/tables/info/BLCA_RNA_somatic_calling_info.tsv \
    --project_folder /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results \
    --exon_interval /home/lqh/resources/database/gencode/GRCh38_GENCODE_v22_exon_rm_alt.bed \
    --output_table_path /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/BLCA/RNA/RNA_somatic_mutation/VcfAssembly_new/SNP_WES_Interval_exon.txt \
    --num_threads 60

    # train your own model
    python /home/lqh/Codes/Python/RNA-SSNV/own_model_construct.py \
    --REDIportal /home/lqh/resources/database/RNA_edit/REDIportal/REDIportal_main_table.hg38.bed \
    --DARNED /home/lqh/resources/database/RNA_edit/DARNED_hg19_to_bed_to_hg38_rm_alt.bed \
    --raw_RNA_mutations /home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/results/LUAD/RNA/RNA_somatic_mutation/VcfAssembly_new/SNP_WES_Interval_exon.txt \
    --DNA_mutations /home/lqh/Codes/Data/TCGA_maf_files/TCGA-LUAD \
    --model_folder_path /home/lqh/Codes/Python/RNA-SSNV/model

Utilize customized model
-------------------------

- Back to the beginning of our pipeline, edit the **model** path within config file, start our pipeline and good to go!

Q & A
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Process failed
--------------------

Check your log file with `grep -C 10 your_log_file.log` 