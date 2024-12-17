# Differential expression analysis

```{r}
# Install packages

install.packages(c("tidyverse", "here", "dplyr"))

# Install Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")

# Install DESeq2 package -> performs differential expression analysis
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2"))

# Install Glimma package -> for interactive plots
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Glimma")
```

```{r}
# Load package libraries

library(tidyverse)
library(here)
library(DESeq2)
library(dplyr)
library(Glimma)
```

```{r}
# Import the data

gene_counts <- read.csv(here("data", "gene_count_matrix.csv"))
colData <- read.csv(here("data", "sample_info.csv"))
```

```{r}
# Tidy the data
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

```{r}
# Create a DESeqDataSet from count matrix and labels
dds <- DESeqDataSetFromMatrix(countData = gene_counts, colData = colData, design = ~ treatment)
dds

# Remove MSTRG IDs from rownames
newrownames <- str_replace(rownames(dds), "MSTRG\\.[0-9]+\\|", "")
rownames(dds) <- newrownames

# Remove duplicate gene names and symbols from rownames
newrownames2 <- str_replace(rownames(dds), "[|].*", "")
rownames(dds) <- newrownames2

# Pre-filter to remove rows with low gene counts across all the samples
keep <- rowSums(counts(dds)) >= 10 # keeps only rows that have at least 10 reads total
dds <- dds[keep,]
dds

# Set the factor level
# Use control as a reference level to compare against the treated
dds$treatment <- relevel(dds$treatment, ref = "control")
```

```{r}
# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)
res

# Specify the coefficient/contrast you want to build a results table for
res <- results(dds, name="treatment_mCD40L_vs_control")
res <- results(dds, contrast=c("treatment","mCD40L","control"))
res
```

```{r}
# p-values and adjusted p-values
resOrdered <- res[order(res$padj),]
summary(resOrdered)

# p < 0.05
res05 <- results(dds, alpha=0.05)
summary(res05)

# How many adjusted p-values were less than 0.05?
sum(res05$padj < 0.05, na.rm = TRUE)

# p < 0.01
res01 <- results(dds, alpha=0.01)
summary(res01)

# How many adjusted p-values were less than 0.01?
sum(res01$padj < 0.01, na.rm = TRUE)

# p < 0.001
res001 <- results(dds, alpha=0.001)
summary(res001)

# How many adjusted p-values were less than 0.001?
sum(res001$padj < 0.001, na.rm = TRUE)
```

## Exploring results

```{r}
# MA plots
plotMA(res, ylim=c(-2,2))

# Interactive Glimma MD-plot
status <- as.numeric(res$padj < .01)
anno <- data.frame(GeneID=rownames(res))
glMDPlot(res, status=status, counts=counts(dds,normalized=TRUE),
         groups=dds$treatment, transform=FALSE,
         samples=colnames(dds), anno=anno, folder="../output/interactive_MD_01", launch=FALSE)

status1 <- as.numeric(res$padj < .001)
glMDPlot(res, status=status1, counts=counts(dds,normalized=TRUE),
         groups=dds$treatment, transform=FALSE,
         samples=colnames(dds), anno=anno, folder="../output/interactive_MD_001", launch=FALSE)
```


```{r}
# Filtering of the results
sig_res <- res_df %>%
  dplyr::filter(padj < 0.01) %>%
  dplyr::filter(log2FoldChange > 2 | log2FoldChange < -2)

# Results tables
res_df <- as.data.frame(res)
sig_res_df <- as.data.frame(sig_res)

res_table <- rownames_to_column(res_df, var = "gene")
sig_res_table <- rownames_to_column(sig_res_df, var = "gene")

sig_res_table %>%
  arrange(padj)
```


```{r}
# countPlots

plotCounts(dds, gene = which.min(res$padj), intgroup = "treatment")
# plots the count for the gene which had the smallest p value

# Customized plotting
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="treatment", 
                returnData=TRUE)
ggplot(d, aes(x=treatment, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0), shape=21) +
  theme_bw() + 
  ggtitle("ENSG00000185215|TNFAIP2") +
  xlab("Treatment") + ylab("Normalized count") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_line(colour = "white")) +
  scale_x_discrete(labels = c("Control", "mCD40L")) +
  scale_y_continuous(trans = "log10")

sig_res_ids <- sig_res_ids %>%
  arrange(padj)

# Most significantly DE protein-coding genes
png(file = "top_9_most_sig_proteincoding_genes.png", width = 900, height = 900, units = "px", pointsize = 24)
head(sig_res_ids, n = 10)
par(mfrow=c(3,3))
plotCounts(dds, gene="ENSG00000185215", intgroup = "treatment", main = "TNFAIP2")
plotCounts(dds, gene="ENSG00000090339", intgroup = "treatment", main = "ICAM1")
plotCounts(dds, gene="ENSG00000056558", intgroup = "treatment", main = "TRAF1")
plotCounts(dds, gene="ENSG00000125347", intgroup = "treatment", main = "IRF1")
plotCounts(dds, gene="ENSG00000077150", intgroup = "treatment", main = "NFKB2")
plotCounts(dds, gene="ENSG00000100906", intgroup = "treatment", main = "NFKBIA")
plotCounts(dds, gene="ENSG00000118503", intgroup = "treatment", main = "TNFAIP3")
plotCounts(dds, gene="ENSG00000110047", intgroup = "treatment", main = "EHD1")
plotCounts(dds, gene="ENSG00000145901", intgroup = "treatment", main = "TNIP1")
dev.off()

# Most significantly DE non-protein coding genes
sig_res_non_prot_cod <- sig_res_ids %>%
dplyr::filter(biotype != "protein_coding")
sig_res_non_prot_cod %>%
  arrange(padj)

png(file = "top_9_most_sig_non-proteincoding_genes.png", width = 900, height = 900, units = "px", pointsize = 24)
head(sig_res_non_prot_cod, n = 10)
par(mfrow=c(3,3))
plotCounts(dds, gene="ENSG00000285427", intgroup = "treatment", main = "SOD2-OT1")
plotCounts(dds, gene="ENSG00000228432", intgroup = "treatment", main = "DHFRP2")
plotCounts(dds, gene="ENSG00000262198", intgroup = "treatment", main = "Novel transcript")
plotCounts(dds, gene="ENSG00000277117", intgroup = "treatment", main = "Artifact")
plotCounts(dds, gene="ENSG00000261618", intgroup = "treatment", main = "LINC02605")
plotCounts(dds, gene="ENSG00000234685", intgroup = "treatment", main = "NUS1P2")
plotCounts(dds, gene="ENSG00000285744", intgroup = "treatment", main = "Novel transcript")
plotCounts(dds, gene="ENSG00000267607", intgroup = "treatment", main = "ICAM4-AS1")
plotCounts(dds, gene="ENSG00000280303", intgroup = "treatment", main = "ERICD")
dev.off()

sig_res_ids <- sig_res_ids %>%
  arrange(desc(log2FoldChange))

# Top upregulated sig DE protein-coding genes
png(file = "top_9_protein_coding_by_l2f1.png", width = 900, height = 900, units = "px", pointsize = 24)
head(sig_res_ids, n = 10)
par(mfrow=c(3,3))
plotCounts(dds, gene="ENSG00000102245", intgroup = "treatment", main = "CD40LG")
plotCounts(dds, gene="ENSG00000162654", intgroup = "treatment", main = "GBP4")
plotCounts(dds, gene="ENSG00000271503", intgroup = "treatment", main = "CCL5")
plotCounts(dds, gene="ENSG00000164400", intgroup = "treatment", main = "CSF2")
plotCounts(dds, gene="ENSG00000130477", intgroup = "treatment", main = "UNC13A")
plotCounts(dds, gene="ENSG00000105371", intgroup = "treatment", main = "ICAM4")
plotCounts(dds, gene="ENSG00000163121", intgroup = "treatment", main = "NEURL3")
plotCounts(dds, gene="ENSG00000137571", intgroup = "treatment", main = "SLCO5A1")
plotCounts(dds, gene="ENSG00000108342", intgroup = "treatment", main = "CSF3")
dev.off()

sig_res_non_prot_cod <- sig_res_non_prot_cod %>%
  arrange(desc(log2FoldChange))
```

```{r}
# plotCounts for genes of interest
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

# NLRP3
# Find the Ensembl Id of NLRP3
grch38[grch38$symbol == "NLRP3", "ensgene"] # ENSG00000162711
# Plot NLRP3 expression
plotCounts(dds, gene = "ENSG00000162711", intgroup = "treatment", main = "NLRP3")

# Export
png(file = "genes_of_interest.png", width = 2700, height = 900, units = "px", pointsize = 24)
head(sig_res_non_prot_cod, n = 10)
par(mfrow=c(1,3))
plotCounts(dds, gene = "ENSG00000101017", intgroup = "treatment", main = "CD40")
plotCounts(dds, gene = "ENSG00000102245", intgroup = "treatment", main = "CD40L")
plotCounts(dds, gene = "ENSG00000162711", intgroup = "treatment", main = "NLRP3")
dev.off()
```

```{r}
# Principal components plot

plotPCA(vsd, "treatment") + theme(
  panel.background = element_rect(fill = "white", colour = "grey50"), 
  panel.border = element_rect(linetype = "solid", fill = NA, colour = "grey50"),
  axis.ticks = element_line(colour = "grey50"))

ggsave(here("output/DE/pca_plot.png"), width = 10, height = 7)
```

```{r}
# Volcano plot
res_table1 <- res_table %>%
  mutate(threshold = padj < 0.001 & abs(log2FoldChange) >= 2) %>%
  mutate(threshold = case_when(
    log2FoldChange >= 2 & padj <= 0.001 ~ "upregulated",
    log2FoldChange <= -2 & padj <= 0.001 ~ "downregulated"
  ))
 
ggplot(res_table1) + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
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

ggsave(here("output/DE/volcano_plot1.png"), width = 10, height = 7)

# Label top 10 genes (lowest padj) on this plot
# Add all gene symbols as a column from grch38
library(annotables)
library(ggrepel)
library(cowplot)

res_table1 <- bind_cols(res_table1, symbol=grch38$symbol[match(res_table1$gene, grch38$ensgene)])

# Create an empty columb to indicate which genes to label
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

ggsave(here("output/DE/labeled_volcano_plot1.png"), width = 10, height = 7)
```

```{r}
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

# Functional analysis

```{r}
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(pathview)
library(AnnotationHub)
library(ensembldb)
```

## GO Enrichment

```{r}
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

## Molecular functions
ego1 <- enrichGO(gene = sig_genes,
                universe = all_genes, 
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db,
                ont = "MF",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

MF <- data.frame(ego1) %>%
  dplyr::filter(p.adjust < 0.01)
write.csv(MF, here("output/GO/molecular_functions.csv"))

## Cellular components
ego2 <- enrichGO(gene = sig_genes,
                universe = all_genes, 
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db,
                ont = "CC",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

CC <- data.frame(ego2) %>%
  dplyr::filter(p.adjust < 0.01)
write.csv(CC, here("output/GO/cellular_components.csv"))
```

```{r}
# Dotplots
dotplot(ego) + scale_color_gradient(high = "#F5A200", low = "#CC5800")
ggsave(here("output/bio_processes_dotplot.png"), width = 10, height = 7)

dotplot(ego, showCategory=10, split=".sign") + facet_grid(.~.sign)

dotplot(ego1)
ggsave(here("output/mol_functions_dotplot.png"), width = 10, height = 7)

dotplot(ego2)
ggsave(here("output/cell_compartments_dotplot.png"), width = 10, height = 7)
```

```{r}
# Enrichment GO plots
library(enrichplot)

ego <- enrichplot::pairwise_termsim(ego)
ego1 <- enrichplot::pairwise_termsim(ego1)
ego2 <- enrichplot::pairwise_termsim(ego2)

emapplot(ego, showCategory = 15, shadowtext = FALSE, repel = TRUE)
ggsave(here("output/bio_processes_emapplot.png"), width = 10, height = 7)

emapplot(ego1, shadowtext = FALSE, repel = TRUE)
ggsave(here("output/mol_functions_emapplot.png"), width = 10, height = 7)

emapplot(ego2, shadowtext = FALSE, repel = TRUE)
ggsave(here("output/cell_compartments_emapplot.png"), width = 10, height = 7)
```

```{r}
# Category netplots

# To colour the genes by log2Fold changes, log2Fold change column needs to be extracted from the results table, creating a vector
DE_foldchanges <- sig_res_ids$log2FoldChange
names(DE_foldchanges) <- sig_res_ids$gene

# Cnetplot details the genes associated with one or more terms
## By default it gives top 5 significant terms (by padj)
DE_foldchanges1 <- ifelse(DE_foldchanges > 2, 2, DE_foldchanges)
DE_foldchanges1 <- ifelse(DE_foldchanges < -2, 2, DE_foldchanges)

cnetplot(ego,
         categorySize="pvalue",
         showCategory = 5,
         foldChange=DE_foldchanges1,
         vertex.label.font=2) + scale_colour_gradient(low = "#1c9e77", high = "#BFFFEC")
ggsave(here("output/cnetplot.png"), width = 10, height = 7)

cnetplot(ego1,
         categorySize="pvalue",
         showCategory = 5,
         foldChange=DE_foldchanges1,
         vertex.label.font=2) 
ggsave(here("output/mol_functions_cnetplot.png"), width = 10, height = 7)

cnetplot(ego2,
         categorySize="pvalue",
         showCategory = 5,
         foldChange=DE_foldchanges1,
         vertex.label.font=2) 
ggsave(here("output/cell_compartments_cnetplot.png"), width = 10, height = 7)
```

## Gene Set Enrichment Analysis (GSEA)

```{r}
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

```{r}
# GSEA using gene sets from KEGG pathways
set.seed(123456) # this is to use the same permutations every time below code is run
gseaKEGG <- gseKEGG(geneList = l2f, # ordered named vector of fold changes (entrez IDs are the associated names)
                     organism = "hsa", 
                     nPerm = 1000, # default number permutations
                     minGSSize = 20, # minimum gene set size (# genes in a set)
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

```{r}
# GSEA plots
# enrichment score peaks at the red dotted line, which is in the positive log2 fold changes, suggesting the up-regulation of this pathway

# hsa04060: Cytokine-cytokine receptor interaction
gseaplot(gseaKEGG, geneSetID = 'hsa04060', title = gseaKEGG$Description[1])
gseaplot2(gseaKEGG, geneSetID = 'hsa04060', title = gseaKEGG$Description[1])

# hsa04668: TNF signaling pathway
gseaplot(gseaKEGG, geneSetID = 'hsa04668', title = gseaKEGG$Description[2])
gseaplot2(gseaKEGG, geneSetID = 'hsa04668', title = gseaKEGG$Description[2])

# hsa04064: NF-kappa B signaling pathway
gseaplot(gseaKEGG, geneSetID = 'hsa04064', title = gseaKEGG$Description[3])
gseaplot2(gseaKEGG, geneSetID = 'hsa04064', title = gseaKEGG$Description[3])

# Multiple gene sets
gp1 <- gseaplot2(gseaKEGG, geneSetID = 1:3)
ggsave(here("output/Combined_GSEA_plot1.png"), width = 12, height = 10)
```

```{r}
# GSEA using gene sets associated with BP Gene Ontology terms
gseaGO <- gseGO(geneList = l2f, 
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                keyType = "ENTREZID",
                minGSSize = 20,
                pvalueCutoff = 0.05,
                verbose = FALSE)
gseaGO_res <- as.data.frame(gseaGO)
gseaplot(gseaGO, geneSetID = 'GO:0006954')
```

```{r}
# GMT file of gene sets
library(GSEABase)

# Load in GMT file of gene sets 
# downloaded from the Broad Institute [website](http://software.broadinstitute.org/gsea/msigdb/collections.jsp) for MSigDB)

# GSEA using MSigDB gene sets
c2 <- read.gmt(here("data", "c2.cp.v2023.1.Hs.entrez.gmt"))
msig <- GSEA(l2f, TERM2GENE=c2, verbose=FALSE)
msig_df <- data.frame(msig)
msig_df
```

## Visualizations

```{r}
# Convert gene ID to symbol
egox <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')
ego1x <- setReadable(ego1, 'org.Hs.eg.db', 'ENTREZID')
ego2x <- setReadable(ego2, 'org.Hs.eg.db', 'ENTREZID')

# Heatmap-like functional classification
sig_l2f_df <- sig_res_ids %>%
  dplyr::filter(log2FoldChange > 2 | log2FoldChange < -2)
sig_l2f <- sig_l2f_df$log2FoldChange
names(sig_l2f) <- sig_l2f_df$symbol

p1 <- heatplot(egox, foldChange=sig_l2f, showCategory=16)
p1
ggsave(here("output/bio_process_hmfc.png"), width = 14, height = 7)

p2 <- heatplot(ego1x, foldChange=sig_l2f)
p2
ggsave(here("output/mol_functions_hmfc.png"), width = 14, height = 7)

p3 <- heatplot(ego2x, foldChange=sig_l2f)
p3
ggsave(here("output/cell_compartment_hmfc.png"), width = 14, height = 7)
```

```{r}
# Tree plot

p4 <- treeplot(egox)
p4
ggsave(here("output/bio_process_tree.png"), width = 14, height = 7)

p5 <- treeplot(ego1x)
p5
ggsave(here("output/mol_functions_tree.png"), width = 14, height = 7)

p6 <- treeplot(ego2x)
p6
ggsave(here("output/cel_compartments_tree.png"), width = 14, height = 7)
```

























