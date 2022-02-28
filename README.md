# CReSIL - Construction based Rolling-circle-amplification for eccDNA Sequence Identification and Location 

Python scripts and pipeline for detecting eccDNA from Nanopore reads



* Creating an environment with commands:
    ```bash
    ## Install Conda environment
    cd cresil
    conda env create -f environment.yml

    ## Activate CReSIL environment
    conda activate cresil

    ## Install CReSIL
    pip install .

    ## Run CReSIL
    cresil -h
    cresil trim -h
    cresil identify -h
    cresil annotate -h
    cresil visualize -h
    ```
    
* Run CReSIL:
    ```bash
    ## Go to an example folder
    cd example
    
    ## Run trim 
    cresil trim -t 4 -i exp_reads.fastq -r hg19.25chr.mmi -o cresil_result

    ## Run eccDNA identification [enrichment dataset]
    cresil identify -t 4 -fa reference.fa -fai reference.fa.fai -fq exp_reads.fastq -trim cresil_result/trim.txt

    ## Run eccDNA identification [whole-genome long-read sequencing dataset]
    cresil identify_wgls -t 4 -r hs19.25chr.mmi -fa reference.fa -fai reference.fa.fai -fq exp_reads.fastq -trim cresil_result/trim.txt

    ## Run eccDNA annotation
    cresil annotate -t 4 -rp reference.rmsk.bed -cg reference.cpg.bed -gb reference.gene.bed -identify cresil_result/eccDNA_final.txt

    ## Run visualize eccDNA
    cresil visualize -t 4 -c ec1 -identify cresil_result/eccDNA_final.txt

    ## Run Circos
    circos -noparanoid -conf circos
    ```
