#!/bin/bash

# STEP 1: Quality control (fastqc)
fastqc data/untrimmed_fastq/*.fq -o results/untrimmed_fastqc

# Run FastQ Screen to check genome contents
fastq_screen ~/analysis/data/untrimmed_fastq/*.fq.gz

# STEP 2: Check RNA-strandedness
pip install how_are_we_stranded_here
mkdir kallisto
cd kallisto
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
kallisto index -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx Homo_sapiens.GRCh38.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz
gunzip Homo_sapiens.GRCh38.108.gtf.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
cd ..
check_strandedness --gtf kallisto/Homo_sapiens.GRCh38.108.gtf --transcripts kallisto/Homo_sapiens.GRCh38.cdna.all.fa --reads_1 data/untrimmed_fastq/LE2_6h_1.fq --reads_2 data/untrimmed_fastq/LE2_6h_2.fq -k kallisto/Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx

# STEP 3: Genome alignment (HISAT2)
for infile in LE1_6h LE2_6h NE1_6h NE2_6h
do
hisat2 -p 20 -q --rna-strandness RF --dta -x HISAT2/grch38/genome -1 data/untrimmed_fastq/${infile}_1.fq -2 data/untrimmed_fastq/${infile}_2.fq -S results/4.2.aligned_dta/${infile}.aligned.sam | echo "Alignment for ${infile} has started."
done

cd ~/analysis/results/5.2.processed_aligned_dta/

# 3.1 Convert SAM to BAM
for infile in LE1_6h LE2_6h NE1_6h NE2_6h
do
samtools view -S -b ~/analysis/results/5.2.processed_aligned_dta/${infile}.aligned.sam > ~/analysis/results/5.2.processed_aligned_dta/${infile}.aligned.bam
done

# 3.2 Remove reads with a mate mapped to a different chromosome
for infile in LE1_6h LE2_6h NE1_6h NE2_6h
do
samtools view -H ${infile}.aligned.bam > ${infile}.aligned.sam
samtools view ${infile}.aligned.bam | awk '$7 == "=" {print}' >> ${infile}.aligned.sam
samtools view -bS ${infile}.aligned.sam > ${infile}.chr_removed.bam
samtools flagstat ${infile}.chr_removed.bam
done

# 3.3 Remove singletons
for infile in LE1_6h LE2_6h NE1_6h NE2_6h
do
samtools view -@ 8 -F 0x08 -b ${infile}.aligned.sam > ${infile}.chr_singleton_removed.bam
samtools flagstat ${infile}.chr_singleton_removed.bam
done

# 3.4 Sort BAM files
for infile in LE1_6h LE2_6h NE1_6h NE2_6h
do
samtools sort -o ${infile}.chr_singleton_removed.sorted.bam ${infile}.chr_singleton_removed.bam
done

# 3.5 Re-run fastq screen to check if there is only human genome left
for infile in LE1_6h LE2_6h NE1_6h NE2_6h
do
samtools bam2fq ${infile}.chr_singleton_removed.sorted.bam > ${infile}.chr_singleton_removed.sorted.fq
fastq_screen ${infile}.chr_singleton_removed.sorted.fq
done

# STEP 4: Transcript assembly (StringTie)

# 4.1 Assemble the reads
for infile in LE1_6h LE2_6h NE1_6h NE2_6h
do
stringtie -p 20 --rf -o results/6.assembled/${infile}.assembled.gtf -G kallisto/Homo_sapiens.GRCh38.108.gtf results/5.2.processed_aligned_dta/${infile}.chr_singleton_removed.sorted.bam
done

# 4.2 Merge the reads
for infile in LE1_6h LE2_6h NE1_6h NE2_6h
do
stringtie --merge -G kallisto/Homo_sapiens.GRCh38.108.gtf -o results/6.assembled/merged.gtf results/6.assembled/${infile}.assembled.gtf
done

# 4.2.1 Add gene IDs to MSTRG gene names
chmod 777 mstrg_prep.pl
mstrg_prep.pl results/6.assembled/merged.gtf > results/6.assembled/merged_prep.gtf

# 4.3 Estimate transcript abundances and generate read coverage tables
for infile in LE1 LE2 NE1 NE2
do
stringtie -e -B -G results/6.assembled/merged_prep.gtf -o results/7.expression_estimation/${infile}.re-estimated_prep.gtf results/5.2.processed_aligned_dta/${infile}_6h.chr_singleton_removed.sorted.bam
done

# 4.4 Run prepDE.py to generate .csv read counts
cd 7.expression_estimation/
prepDE.py
