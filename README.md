# CD40 signalling in the carcinoma cells
This was a 12-month research placement investigating the mechanisms of transcriptional and post-translational regulation of gene expression following CD40 receptor ligation in carcinoma cells.

Required packages: FastQC, FastQ Screen, HISAT2, SAMtools, StringTie.

**Overview**
Step | Procedure | Outputs
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

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a quality control package which produces a sequence quality report, which is useful in determining whether or not your sequences need to be trimmed before the next step. 

[FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) is a package which allows you to screen your sequence library against a set of sequence databases to see the composition of your sequence library. This is an optional step but I did it out of curiosity, because the CD40 ligand (CD40L) used for receptor ligation in human carcinoma cells, was membrane-bound on a murine fibroblast, so it was expected that there would be a mixture of murine and human genome in the sequenced reads. 

The latest version of both [FastQC](https://bioconda.github.io/recipes/fastqc/README.html) and [FastQ Screen](https://bioconda.github.io/recipes/fastq-screen/README.html?) can be installed using conda:
```
conda install fastqc | conda install fastq-screen
```
## Genome alignment
