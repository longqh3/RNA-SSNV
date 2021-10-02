RNA-Mutect2
=======================================================

The RNA-Mutect2 pipeline is a scalable and efficient analysis pipeline for RNA somatic mutation detection from RNA-WES(tumor-normal) paired sequencing data which utilized Mutect2 as core-caller and ML-filter model to maximize precision & recall. It runs highly automated once snakemake and related-codes get configurated properly. It reports aggregated mutation file (maf format) for given cases to facilitate downstream analysis. 

Pre-requirements
~~~~~~~~~~~~~~~~~

Download codes
----------------------

.. code:: sh
    
  git clone https://github.com/longqh3/Machine-learning-based-RNA-somatic-mutation-detection

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

Change interval files
---------------------

- *resources/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals_add_chr_to_hg38_rm_alt.bed*: bed-format interval file for paired-normal Whole Exome Sequence(WES) targets, canonical for TCGA projects. (change it with your own WES target interval file)
- *resources/GRCh38_GENCODE_v22_exon_rm_alt.bed*: bed-format interval file for GENCODE v22 exon regions. (change it with your own exon regions of interest)

Run Pipeline
~~~~~~~~~~~~~~~

Once configurated correctly, our pipeline is ready to go. Please execute the following commands line by line, make sure everything works fine before moving forward. 

.. code:: sh
    
    # dry run to see if everything works
    snakemake --cores num_of_cores \
    -s rules/RNA-Somatic-tsv-Snakefile.smk \
    --configfile configs/LUAD_RNA_Somatic_config.yaml \
    -n
    # run the pipeline
    snakemake --cores num_of_cores \
    -s rules/RNA-Somatic-tsv-Snakefile.smk \
    --configfile configs/LUAD_RNA_Somatic_config.yaml
    # run feature-extraction codes
    python codes/4_10_newest_validation_vcf_info_retriver_tsv.py
    # run model-discriminant codes
    python codes/new_function_based_RNA_somatic_random_forest_exon_analysis.py
    # run downstream analysis
    R cancer_survival_analysis.R

Beware, owing to the breakpoint-run feature of snakemake, our pipeline also supports taking any final files (listed below) as starting point. 

Output folders & files
~~~~~~~~~~~~~~~~~~~~~~~~~~

The pipeline outputs several folders containing intermediate files and final project-level mutations annotation file (maf format). Here, we describe the `results/` folder's schema. 

Sequencing data pre-process
------------------------------

- *results/project_name/RNA/marked_duplicates*: temporary folder containing MarkDuplicates tool's output.
- *results/project_name/RNA/splited_n_cigar_reads*: temporary folder containing SplitNCigarReads tool's output.
- `results/project_name/RNA/base_reclibrate`: temporary folder containing BaseRecalibrate tool's output.
- *results/project_name/RNA/apply_BQSR*: permanent folder containing ApplyBQSR tool's output, **final** files (bam format) used to call RNA somatic mutations, applicable for other analysis.

Calling process - RNA somatic mutation
-----------------------------------------

- *results/project_name/RNA/RNA_somatic_mutation/Mutect2*: permanent folder containing Mutect2 tool's output. 
- *results/project_name/RNA/RNA_somatic_mutation/GetPileupSummaries*: permanent folder containing GetPileupSummaries tool's output (best normal sample's pileup summary info).
- *results/project_name/RNA/RNA_somatic_mutation/FilterMutectCalls*: permanent folder containing FilterMutectCalls tool's output, **final** files (vcf format) used to discriminate true RNA somatic mutations, applicable for other filtering strategy. 

Refining process - RNA somatic mutation
----------------------------------------

- *results/project_name/RNA/RNA_somatic_mutation/Funcotator/SNP*: permanent folder containing Funcotator's annnotation info for raw RNA SNP calls. 
- *results/project_name/RNA/RNA_somatic_mutation/SelectVariants/SNP_WES_interval*: permanent folder containing raw RNA SNP calls subsetted via given WES target intervals. 
- *results/project_name/RNA/RNA_somatic_mutation/SelectVariants/SNP_WES_interval_exon*: permanent folder containing **final** raw RNA SNP calls subsetted via given WES target intervals and exon regions.

Filtering process - RNA somatic mutation
----------------------------------------

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

Although we used 511 cases of TCGA LUAD RNA-WES paired data to train our discriminant model, other non-cancerous RNA somatic mutations may exhibit **different patterns of FP calls**. In that case, our model may not served as expected, and a customized model was required to be trained on your own.

Data-preparation
--------------------

- Gold-standard TP mutations for given project (maf-format)

Train customized model
-----------------------

- Using TP mutations to assess FilterMutectCalls tool's performance - output a performance matrix and UpSet plot to visualize FN's patterns. 

  .. code:: sh
      
    python custom_model_train_preparation.py

- Specify FN patterns (list of filters within FilterMutectCalls outputs) and corresponding folder paths - output a performance matrix for model training. 

  .. code:: sh
        
    python custom_model_train_process.py

Utilize customized model
-------------------------

- Back to the beginning of our pipeline, edit the model path within config file, start our pipeline and good to go!

Q & A
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Process failed
--------------------

Check your log file with `grep -C 10 your_log_file.log`