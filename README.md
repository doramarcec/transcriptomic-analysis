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
conda install fastqc
conda install fastq-screen
```

### RNA-strandedness

Before moving onto the genome alignment step, I needed to figure out the RNA strandedness of my reads to make sure that I use the right strandedness setting when running HISAT2. To do this, I used a [custom python script](https://github.com/signalbash/how_are_we_stranded_here) and [kallisto](https://pachterlab.github.io/kallisto/manual).
```
pip install how_are_we_stranded_here
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
kallisto index -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx Homo_sapiens.GRCh38.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz
gunzip Homo_sapiens.GRCh38.108.gtf.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
cd ..
```
GTF transcript annotation file, transcripts FASTA file, and FASTQ read files from one sample are needed to run the script.
```
check_strandedness --gtf kallisto/Homo_sapiens.GRCh38.108.gtf --transcripts kallisto/Homo_sapiens.GRCh38.cdna.all.fa --reads_1 data/untrimmed_fastq/LE2_6h_1.fq --reads_2 data/untrimmed_fastq/LE2_6h_2.fq -k kallisto/Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx
```

## Genome alignment

HISAT2 is an alignment program which maps NGS reads (both RNA and DNA) to population of human genomes as well as to a single reference genome.

The latest version of HISAT2 can be installed using [conda](https://bioconda.github.io/recipes/hisat2/README.html?):
```
conda install hisat2
```

To automate the genome alignment across all of my samples I used the *foor* loop as below:
```
for infile in LE1_6h LE2_6h NE1_6h NE2_6h
do
hisat2 -p 20 -q --rna-strandness RF --dta -x HISAT2/grch38/genome -1 data/untrimmed_fastq/${infile}_1.fq -2 data/untrimmed_fastq/${infile}_2.fq -S results/4.2.aligned_dta/${infile}.aligned.sam | echo "Alignment for ${infile} has started."
done
```

**Explanation of used settings**
Argument | Explanation
--- | ---
-p | Used for performance tuning - if the computer has multiple processors/cores, -p instructs HISAT2 on how many to use.
-q | Specifies that the input reads are FASTQ files.
--RNA-strandedness | Takes strand-specific information, e.g. unstranded (the default), F (forward-stranded) or R (reverse-stranded) single-end reads, RF (reverse-stranded) or FR (direct-stranded) paired-end reads.
--dta / --downstream-transcriptome-assembly | Reports alignments tailored for transcript assemblers such as StringTie.
-x | The path to and the basename of the reference genome indices.
-1 | Files containing mate 1s (filename usually includes _1). If FR, this will be the reads from the sense (forward) strand. If RF, this will be reads from the antisense (reverse) strand.
-2 | Files containing mate 2s (filename usually includes _2). If FR, this will be the reads from the antisense (reverse) strand. If RF, this will be reads from the sense (forward) strand.
-S | Path to and the basename of the output SAM alignment file. 

Find out more in the [HISAT2 manual.](https://daehwankimlab.github.io/hisat2/manual/) 

### SAMtools




