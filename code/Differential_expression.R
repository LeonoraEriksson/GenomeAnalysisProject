
# Install correct packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library("DESeq2")

install.packages("pheatmap")
library ("pheatmap")

# Read files
read_counts <- read.table("C:/Users/Admin/Documents/GenomeAnalysis/Read_counts.txt", sep = '\t', header = FALSE)
data_types <- read.table("C:/Users/Admin/Documents/GenomeAnalysis/Sample_types.txt", sep = '\t', header = FALSE)
gene_annotation = read.table("C:/Users/Admin/Documents/GenomeAnalysis/Genes_ID_Description.txt",sep = '\t', header = FALSE)

# Create appropriate dataframes
read_counts[,1] = as.data.frame(gsub('_g','', read_counts[,1])) # Remove _g from gene id/quiery id for merging
gene_annotation[,1] <- as.data.frame(gsub('_t','', gene_annotation[,1])) # Remove _t from gene id/quiery id for merging
gene_annotation[,2] <- as.data.frame(gsub(' ', '_', gene_annotation[,2])) # Make gene description readable for DESeq
gene_annotation[,2] <- as.data.frame(gsub('-', '_', gene_annotation[,2]))

colnames(gene_annotation) = c('ID','Name')
colnames(read_counts) = c('ID',"Aril_1", "Aril_2", "Aril_3", "Leaf", "Root", "Stem" )

genes = merge(read_counts, gene_annotation, by='ID', all=TRUE) # Merge read count data and gene descriptions based on id
genes[is.na(genes)] <- "" # Make gene description readable for DESeq

# Create dataframe for DESeq
countData <- as.matrix(genes[ ,c(2,3,4,5,6,7)])
rownames(countData) <- genes[ , 8]
colnames(countData) <- c("Aril_1", "Aril_2", "Aril_3", "Leaf", "Root", "Stem")
colData <- as.matrix(data_types[ , 2:3])
colnames(colData) <- c("type", "organ")
rownames(colData) <- c("Aril_1", "Aril_2", "Aril_3", "Leaf", "Root", "Stem")

# Create DESeq dataset
dds = DESeqDataSetFromMatrix(countData=countData, colData = colData, design= ~ organ)
dds_type =  DESeqDataSetFromMatrix(countData=countData, colData = colData, design= ~ type)

# Do differential expression analysis, comparing fruit organs vs non-fruit organs
dds <- DESeq(dds)
dds$type <- factor(dds$type, levels = c("Aril", "Leaf", "Root", "Stem"))
dds$organ <- factor(dds$organ, levels = c("Fruit", "No_fruit"))
res <- results(dds, contrast=c("organ", "Fruit", "No_fruit"))

# Order results by smallest adjusted p-value
resOrdered <- res[order(res$padj),]
summary(res)

# How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)

# Change cut-off adj. p-value from 0.1 to 0.05
res05 <- results(dds, alpha=0.05)
summary(res05)
# How many adjusted p-values were less than 0.05?
sum(res05$padj < 0.05, na.rm=TRUE)

# Create subset of the results with only the significant genes
resSig <- subset(resOrdered, padj < 0.05)
resSig

# Examine counts for gene with the smallest pvalue
plotCounts(dds, gene=which.min(res$padj), intgroup="organ")      

 #What tests were done?
mcols(res)$description

#Create heatmap gene count top 20 genes
rld <- rlog(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("organ", "type")])
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df)

#PCA plot
plotPCA(rld, intgroup=c("type"))
