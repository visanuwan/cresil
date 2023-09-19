# CReSIL: Accurate Identification of Extrachromosomal Circular DNA from Long-read Sequences

A tool for detecting eccDNA from Nanopore reads



## Installation instructions
```
## Install Conda environment
cd cresil
conda env create -f environment.yml

## Activate CReSIL environment
conda activate cresil

## Install CReSIL
pip install .
```

## Testing (Human eccDNA)
```
## Go to an example folder
cd example

## Run trim 
cresil trim -t 4 -fq exp_reads.fastq -r reference.mmi -o cresil_result

## Run eccDNA identification for enriched data
cresil identify -t 4 -fa reference.fa -fai reference.fa.fai -fq exp_reads.fastq -trim cresil_result/trim.txt

## Run eccDNA annotation
cresil annotate -t 4 -rp reference.rmsk.bed -cg reference.cpg.bed -gb reference.gene.bed -identify cresil_result/eccDNA_final.txt

## Run visualize eccDNA (ec1)
cresil visualize -t 4 -c ec1 -identify cresil_result/eccDNA_final.txt

## Run Circos (ec1)
cd cresil_result/for_Circos/ec1
circos -noparanoid -conf circos.conf
```

The structure of files and directories generated from CReSIL is as follows:

```
cresil_result/
├── trim.txt
├── eccDNA_final.txt
├── cresil_run/
│   ├── subGraphs.summary.txt
│   ├── tmp/
│   └── assemGraph/
│       ├── ec1/
│       ├── ec2/
│       └── ..
├── cresil_gAnnotation/
│   ├── gene.annotate.txt
│   ├── CpG.annotate.txt
│   ├── repeat.annotate.txt
│   └── variant.annotate.txt
└── for_Circos/
    └── ec1/
        ├── circos.conf
        ├── ticks.conf
        ├── ideogram.conf
        ├── karyotype.conf
        ├── data/
        └── tmp/
```
For more details about mapped reads for each eccDNA, please see the eccDNA folder in **assemGraph**.

## Additional useful commands

To run CReSIL to identify eccDNA in a relaxed mode (from the trimming step)

```
## Run eccDNA identification for with minimizing parameters (minimum size of eccDNA, average depth, number of supported breakpoints)
cresil identify -t 4 -minrsize 40 -depth 1 -break 1 -fa reference.fa -fai reference.fa.fai -fq exp_reads.fastq -trim cresil_result/trim.txt
```

To use CReSIL to identify eccDNA in whole-genome long-read (WGLS) sequencing data (from the trimming step)

```
## Run eccDNA identification for whole-genome long-read (WGLS) sequencing data
cresil identify_wgls -t 4 -r reference.mmi -fa reference.fa -fai reference.fa.fai -fq exp_reads.fastq -trim cresil_result/trim.txt
```

## Interface
 - `cresil trim` - find and trim potential eccDNA regions from ONT reads
 - `cresil identify` - identify and verify eccDNA from enriched data
 - `cresil identify_wgls` - identify and verify eccDNA from whole-genome long-read (WGLS) sequencing data
 - `cresil annotate` - annotate identified eccDNA
 - `cresil visualize` - generate Circos configuration files for specified eccDNA

## Sample data
[References and simulated data](https://app.box.com/s/5leixacmp1xx8qs7qtcuyz93lgqpzgkn) for the user to go through the example of CReSIL pipeline and used in the CReSIL benchmarks.
- References
- Simulated data

## Citation
Please cite the following article if you use CReSIL in your research
> Visanu Wanchai, Piroon Jenjareonpun, Thongpun Leangapichat, Gerard Arrey, Charles M Burnham, Maria C Tümmle, Jesus Delgado-Calle, Birgitte Regenberg, Intawat Nookaew, CReSIL: Accurate Identification of Extrachromosomal Circular DNA from Long-read Sequences, *Briefings in Bioinformatics.* 2022;, bbac422, https://doi.org/10.1093/bib/bbac422.<br>

## License and Copyright

CReSIL is distributed under the terms of the MIT License
