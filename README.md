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

```
fastqc data/untrimmed_fastq/*.fq -o results/untrimmed_fastqc
```

[FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) is a package which allows you to screen your sequence library against a set of sequence databases to see the composition of your sequence library. This is an optional step but I did it out of curiosity, because the CD40 ligand (CD40L) was presented to the CD40 receptor, in human carcinoma cells, in a membrane bound form on a murine fibroblast, so it was expected that there would be a mixture of murine and human genome in the sequenced reads (there was a 50:50 mixture due to a co-culture). 

```
fastq_screen ~/analysis/data/untrimmed_fastq/*.fq.gz
```

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
-p | Used for performance tuning - if the computer has multiple CPUs, -p specifies how many to use.
-q | Specifies that the input reads are FASTQ files.
--RNA-strandedness | Takes strand-specific information, e.g. unstranded (the default), F (forward-stranded) or R (reverse-stranded) single-end reads, RF (reverse-stranded) or FR (direct-stranded) paired-end reads.
--dta / --downstream-transcriptome-assembly | Reports alignments tailored for transcript assemblers such as StringTie.
-x | The path to and the basename of the reference genome indices.
-1 | Files containing mate 1s (filename usually includes _1). If FR, this will be the reads from the sense (forward) strand. If RF, this will be reads from the antisense (reverse) strand.
-2 | Files containing mate 2s (filename usually includes _2). If FR, this will be the reads from the antisense (reverse) strand. If RF, this will be reads from the sense (forward) strand.
-S | Path to and the basename of the output SAM alignment file. 

Find out more in the [HISAT2 manual.](https://daehwankimlab.github.io/hisat2/manual/) 

### Post-alignment processing

Post-alignment processing involved sorting the aligned reads by their genomic location, and (if needed) converting the file to a BAM format, which contains binary, compressed version of the SAM files, accounting for a better storage format and faster processing. To do this, I used SAMtools, a program commonly used to manipulate aligned reads in the SAM, BAM, or CRAM formats. 

**SAM to BAM conversion**
```
for infile in LE1_6h LE2_6h NE1_6h NE2_6h
do
samtools view -S -b ~/analysis/results/5.2.processed_aligned_dta/${infile}.aligned.sam > ~/analysis/results/5.2.processed_aligned_dta/${infile}.aligned.bam
done
```

Before sorting the BAM files and moving onto transcript assembly, I needed to make sure that all of the murine genome is removed, and I did that by removing singletons and reads with a mate mapped to a different chromosome.

**Removal of reads with a mate mapped to a different chromosome**
```
for infile in LE1_6h LE2_6h NE1_6h NE2_6h
do
samtools view -H ${infile}.aligned.bam > ${infile}.aligned.sam
samtools view ${infile}.aligned.bam | awk '$7 == "=" {print}' >> ${infile}.aligned.sam
samtools view -bS ${infile}.aligned.sam > ${infile}.chr_removed.bam
samtools flagstat ${infile}.chr_removed.bam
done
```

**Singleton removal**
```
for infile in LE1_6h LE2_6h NE1_6h NE2_6h
do
samtools view -@ 8 -F 0x08 -b ${infile}.aligned.sam > ${infile}.chr_singleton_removed.bam
samtools flagstat ${infile}.chr_singleton_removed.bam
done
```

**Sorting BAM files**
```
for infile in LE1_6h LE2_6h NE1_6h NE2_6h
do
samtools sort -o ${infile}.chr_singleton_removed.sorted.bam ${infile}.chr_singleton_removed.bam
done
```

Detailed instructions and explanations of used settings can be found in [SAMtools manual](http://www.htslib.org/doc/samtools.html).

Before moving onto transcript assembly, I re-run FastQ Screen to check if all of the murine genome was removed.
```
for infile in LE1_6h LE2_6h NE1_6h NE2_6h
do
samtools bam2fq ${infile}.chr_singleton_removed.sorted.bam > ${infile}.chr_singleton_removed.sorted.fq
fastq_screen ${infile}.chr_singleton_removed.sorted.fq
done
```

## Transcript assembly

StringTie is a transcript assembler and quantification tool for RNA-seq alignments. It is compatible with the HISAT2 output, and its' output is compatible with DESeq2, used in the next step, making it ideal for my workflow. 

**Assembling the reads**
```
for infile in LE1_6h LE2_6h NE1_6h NE2_6h
do
stringtie -p 20 --rf -o results/6.assembled/${infile}.assembled.gtf -G kallisto/Homo_sapiens.GRCh38.108.gtf results/5.2.processed_aligned_dta/${infile}.chr_singleton_removed.sorted.bam
done
```

**Merging the reads**
```
for infile in LE1_6h LE2_6h NE1_6h NE2_6h
do
stringtie --merge -G kallisto/Homo_sapiens.GRCh38.108.gtf -o results/6.assembled/merged.gtf results/6.assembled/${infile}.assembled.gtf
done
```

**Using a [perl script](https://gist.github.com/gpertea/b83f1b32435e166afa92a2d388527f4b) to add gene IDs to MSTRG gene names from the output**
```
chmod 777 mstrg_prep.pl # permission to execute the script
mstrg_prep.pl results/6.assembled/merged.gtf > results/6.assembled/merged_prep.gtf
```

**Estimation of transcript abundances**
```
for infile in LE1 LE2 NE1 NE2
do
stringtie -e -B -G results/6.assembled/merged_prep.gtf -o results/7.expression_estimation/${infile}.re-estimated_prep.gtf results/5.2.processed_aligned_dta/${infile}_6h.chr_singleton_removed.sorted.bam
done
```

**Explanation of used settings**
Argument | Explanation
--- | ---
-p | Used for performance tuning - if the computer has multiple CPUs, -p specifies how many to use.
--rf | Specifies RF RNA-strandedness.
-o | Output path and basename.
-G | Reference genome annotation file.
--merge | Transcript merge mode. 
-e | Instructs StringTie to operate in expression estimation mode. Used in combination with -G and the output includes the expressed reference transcripts as well as any novel transcripts that were assembled in the GTF file format.
-B | Outputs CTAB file for Ballgown input (was not necessary in the end).

More information can be found in the [StringTie manual](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual).

Finally, to generate read counts in the CSV file format from the estimated transcript abundances, [prepDE.py](https://github.com/gpertea/stringtie/blob/master/prepDE.py) script was used. 



