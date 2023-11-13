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

Finally, to generate read counts in the CSV file format from the estimated transcript abundances, [prepDE.py](https://github.com/gpertea/stringtie/blob/master/prepDE.py) script was used, and generated read count matrix was exported as a CSV file and imported into R. 

## Differential gene expression analysis in R

Before importing the read count matrix and the metadata, required packages were installed and loaded from the library:

```
install.packages(c("tidyverse", "here", "dplyr"))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "Glimma"))
```
Data import, tidying and the analysis itself were performed following the [Analyzing RNA-seq data with DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) vignette by the respective authors of the DESeq2 package.

In order for DESeq2 to be able to perform the differential expression analysis, gene count matrix file has to be accompanied by the file containing metadata. Moreover, the data needs to be tidied and the row names in the metadata file (colData) need to match the column names in the gene count matrix (gene_counts).

```
gene_counts <- read.csv(here("data", "gene_count_matrix.csv"))
colData <- read.csv(here("data", "sample_info.csv"))

# Are the row names in colData matching the column names in gene_counts?
all(rownames(colData) %in% colnames(gene_counts)) # FALSE

gene_counts <- column_to_rownames(gene_counts, var = "gene_id")
gene_counts

colData <- column_to_rownames(colData, var = "sample_id")
colData

all(rownames(colData) %in% colnames(gene_counts)) #TRUE

# Are the row names in colData and column names in gene_counts in the same order?
all(rownames(colData) == colnames(gene_counts)) #TRUE
```

After the data is tidied, a DESeqDataSet is created from the count matrix:

```
dds <- DESeqDataSetFromMatrix(countData = gene_counts, colData = colData, design = ~ treatment)
dds
```

As I have used StringTie for transcript abundance estimation, the gene names in my dataset involved MSTRG IDs, which were removed, along with the duplicate genes names.

```
# Remove MSTRG IDs from rownames
newrownames <- str_replace(rownames(dds), "MSTRG\\.[0-9]+\\|", "")
rownames(dds) <- newrownames

# Remove duplicate gene names and symbols from rownames
newrownames2 <- str_replace(rownames(dds), "[|].*", "")
rownames(dds) <- newrownames2
```

Now that I've tidied up my rownames, I removed rows with low gene counts (< 10) across all the samples, and set the factor level, which is essentially a "reference" level against which the expression in treated samples is going to be compared.

```
# Pre-filter to remove rows with low gene counts across all the samples
keep <- rowSums(counts(dds)) >= 10 # keeps only rows that have at least 10 reads total
dds <- dds[keep,]
dds

# Set the factor level
# Use control as a reference level to compare against the treated
dds$treatment <- relevel(dds$treatment, ref = "control")
```

After this, the data is ready and the differential expression analysis can be performed by DESeq2:

```
dds <- DESeq(dds)
res <- results(dds)
res
```

I filtered the results to the genes that had a p adjusted value < 0.01 and a log2 fold change in the expression levels > 2, in either direction. 

```
# Filtering of the results
sig_res <- res_df %>%
  dplyr::filter(padj < 0.01) %>%
  dplyr::filter(log2FoldChange > 2 | log2FoldChange < -2)

# Results tables
res_df <- as.data.frame(res)
sig_res_df <- as.data.frame(sig_res)
```

The exploratory analysis as well as all the visualisations I did can be found in the RNA-seq-analysis-workflow.Rmd file, but as visualisation is my favourite part of any type of data analysis, here are some of my highlights :)

**plotCounts for genes of interest**
```
# CD40
# Find the Ensembl Id of CD40
grch38[grch38$symbol == "CD40", "ensgene"] # ENSG00000101017
# Plot CD40 expression
plotCounts(dds, gene = "ENSG00000101017", intgroup = "treatment", main = "CD40")

# CD40L
# Find the Ensembl Id of CD40L
grch38[grch38$symbol == "CD40LG", "ensgene"] # ENSG00000102245
# Plot CD40L expression
plotCounts(dds, gene = "ENSG00000102245", intgroup = "treatment", main = "CD40L")
```

**PCA plot**
```
plotPCA(vsd, "treatment") + theme(
  panel.background = element_rect(fill = "white", colour = "grey50"), 
  panel.border = element_rect(linetype = "solid", fill = NA, colour = "grey50"),
  axis.ticks = element_line(colour = "grey50"))
```

**Volcano plot with labels of top 10 significantly DE genes**
```
# Additional packages required to execute the code below
library(annotables)
library(ggrepel)
library(cowplot)

# Label top 10 genes (lowest padj) on this plot
# Add all gene symbols as a column from grch38
res_table1 <- bind_cols(res_table1, symbol=grch38$symbol[match(res_table1$gene, grch38$ensgene)])

# Create an empty column to indicate which genes to label
res_table1 <- res_table1 %>% mutate(genelabels = "")

# Sort by padj
res_table1 <- res_table1 %>% arrange(padj)

# Populate the genelabels column with contents of the gene symbols column for the first 10 rows (top 10 most significantly dif expressed genes)
res_table1$genelabels[1:10] <- as.character(res_table1$symbol[1:10])
res_table1

ggplot(res_table1) + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) + 
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = genelabels)) +
  ggtitle("Gene expression following CD40 ligation") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") + theme(
  panel.background = element_rect(fill = "white", colour = "grey80"), 
  panel.border = element_rect(linetype = "solid", fill = NA, colour = "grey80"),
  axis.ticks = element_line(colour = "grey80"),
  axis.title = element_text(size = rel(1.15)),
  axis.text = element_text(size = rel(1.05)),
  plot.title = element_text(size = rel(1.5), hjust = 0.5),
  legend.position = "top",
  legend.key = element_rect(fill = "white", color = NA),
  legend.title = element_blank()) +
  scale_color_manual(values = c("#1c9e77", "#d95f02"), na.translate = F)
```

Before moving onto the functional analyses, the DE results have to be sorted as below:

```
# Sort out the DE results for Functional Analysis

# Return the IDs for the gene symbols in the DE results
idx <- grch38$ensgene %in% rownames(res_df)
ids <- grch38[idx, ]

# Remove duplicate IDs
non_duplicates <- which(duplicated(ids$ensgene) == FALSE)
ids <- ids[non_duplicates, ]

# Merge the IDs with the results
res_ids <- inner_join(res_table, ids, by=c("gene"="ensgene")) 

# Repeat the above steps for significant genes
a <- grch38$ensgene %in% rownames(sig_res_df)
b <- grch38[a, ]
non_duplicates2 <- which(duplicated(b$ensgene) == FALSE)
b <- b[non_duplicates2, ]
sig_res_ids <- inner_join(sig_res_table, b, by=c("gene"="ensgene"))

sig_res_ids %>%
  arrange(padj)

# Create a background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- as.character(res_ids$gene)

# Significant results
sig_genes <- as.character(sig_res_ids$gene)
```

## Functional profiling in R

For functional profiling, I performed the GO enrichment analysis, KEGG pathway annotation and Gene Set Enrichment Analysis (GSEA). The required packages are outlined below:

```
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(pathview)
library(AnnotationHub)
library(ensembldb)
```

**GO enrichment**
```
## Biological processes
ego <- enrichGO(gene = sig_genes,
                universe = all_genes, 
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

BP <- data.frame(ego) %>%
  dplyr::filter(p.adjust < 0.01)
write.csv(BP, here("output/GO/biological_processes.csv"))
```

This is the code exerpt from the GO enrichment to obtain the enriched biological processes, but the same applies for molecular functions and cellular components, only the "BP" in *ont* parameter has to be changed to "MF" or "CC".

For the KEGG and GSEA analyses, the data had to be slightly manipulated, as below.
```
# Remove NAs
res_entrez <- sig_res_ids %>%
  dplyr::filter(entrez != "NA")

# Remove Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrez) == F), ]

# Extract the fold changes
l2f <- res_entrez$log2FoldChange

# Name each fold change with the corresponding Entrez ID
names(l2f) <- res_entrez$entrez

# Sort fold changes in descending order
l2f <- sort(l2f, decreasing = TRUE)
head(l2f)
```

**GSEA and KEGG pathway annotation**
```
# GSEA using gene sets from KEGG pathways
set.seed(123456) # this is to use the same permutations every time below code is run
gseaKEGG <- gseKEGG(geneList = l2f, 
                     organism = "hsa", 
                     nPerm = 1000, 
                     minGSSize = 20, 
                     pvalueCutoff = 0.05,
                     verbose = FALSE)

# Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result
gseaKEGG_results <- setReadable(gseaKEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

# Save results to csv file
view(gseaKEGG_results)
write.csv(gseaKEGG_results, here("output/KEGG/GSEA_KEGG_results.csv"), quote=F)

# Output images for a single significant KEGG pathway
detach("package:dplyr", unload=TRUE) # first unload dplyr to avoid conflicts
pathview(gene.data = l2f,
         pathway.id = "hsa04657",
         species = "hsa",
         limit = list(gene = 2, # value gives the max/min limit for fold changes
                      cpd = 1))

# Output image for all significant pathways
kegg_plots <- function(x) {
  pathview(gene.data = l2f, pathway.id = gseaKEGG_results$ID[x], species = "hsa", limit = list(gene = 2, cpd = 1))
}
purrr::map(1:length(gseaKEGG_results$ID), kegg_plots)
```

**GSEA plots**

When it comes to visualising the enriched gene sets, each gene set can be plotted invididually, like in 1) and 2), or together with the other gene sets in a single plot, like in 3). 

```{r}
# 1)
gseaplot(gseaKEGG, geneSetID = 'hsa04060', title = gseaKEGG$Description[1])
gseaplot(gseaKEGG, geneSetID = 'hsa04668', title = gseaKEGG$Description[2])
gseaplot(gseaKEGG, geneSetID = 'hsa04064', title = gseaKEGG$Description[3])

# 2)
gseaplot2(gseaKEGG, geneSetID = 'hsa04060', title = gseaKEGG$Description[1])
gseaplot2(gseaKEGG, geneSetID = 'hsa04668', title = gseaKEGG$Description[2])
gseaplot2(gseaKEGG, geneSetID = 'hsa04064', title = gseaKEGG$Description[3])

# 3)
gp1 <- gseaplot2(gseaKEGG, geneSetID = 1:3)
ggsave(here("output/Combined_GSEA_plot1.png"), width = 12, height = 10)
```

