# CD40 signalling in the carcinoma cells
This was a 12-month research placement investigating the mechanisms of transcriptional and post-translational regulation of gene expression following CD40 receptor ligation in carcinoma cells.

Required packages: FastQC, FastQ Screen, HISAT2, SAMtools, StringTie.

**Overview**
Step |Procedure | Outputs
------------ | ------------ | -------------
1 | paired-end RNA-seq | *H. sapiens* FASTQ files
2 | FastQC | N/A
3 | FastQ Screen | N/A
4 | HISAT2 | aligned SAM files
5 | SAMtools | processed and sorted BAM files
6 | FastQ Screen | N/A
7 | StringTie | assembled GTF files, merged GTF file, GTF files with estimated transcript abundances
8 | prepDE.py | CSV read counts

Fist step before starting the analysis is to install [conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html) package manager to be able to install the required packages.

## Quality control

