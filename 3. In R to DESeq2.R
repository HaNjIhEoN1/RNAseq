library(dplyr) # data wrangling
library(ggplot2) # plotting
library(DESeq2) # rna-seq
library(edgeR) # rna-seq
library(tximport) # importing kalisto transcript counts to geneLevels
library(readr) # Fast readr of files.
library(rhdf5) # read/convert kalisto output files.  

# If tximport and rhdf5 are not installed you can install them now.
# Prepare about packages's version
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("tximport")
BiocManager::install("rhdf5")

# metadata parsing 
## case 1
mr <- read.csv("data/metadata_raw.csv",header=TRUE,stringsAsFactors=F,row.names=1)
## case 2
mr<- read.csv("data/metadata_raw.txt",sep = "\t")

# making abundance.tsv files's path
files <- paste("dir",
               list.files(path = "dir",pattern = "abundance.tsv / quant.sf", recursive = TRUE),
               sep = "/")

# designate names for each file
names(files) <- mr$SampleName # SampleName or Run

# ref cdna file parsing
tx2gene <- read.csv("ref.cdna.csv")

# creating kallisto.tsv
txi.kallisto.tsv <- tximport(files, type = "program name", tx2gene = tx2gene, ignoreAfterBar = TRUE)

# making DESeq file 
mr = mr %>% mutate($factor = as.factor($factor))
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, mr, ~$factor) # >>>>> to normalieze true
dds$'factor' <- relevel(dds$'factor', ref = '$ref')

keep <- rowSums(counts(dds)) >= n # n = filtering reads count
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds)

# filtering
up <- filter(as.data.frame(res),padj<0.05,log2FoldChange>0.1)
down <- filter(as.data.frame(res),padj<0.05,abs(log2FoldChange)<1)
 resFilt <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
res <- as.data.frame(resfilt)
write.csv(res,file="path/result.csv")

#----------------------------------normalized = true version-------------------------------------------------
head(counts(dds)
norm <- counts(dds, normalized = TRUE)
write.table(as.data.frame(norm),file = "./norm.txt") ## norm.txt == filename 
aa <- read.table("norm.txt",header = TRUE,row.names = 1)
aa <- round(aa, 0)
head(aa)
deseq2.coldata <- data.frame(condition=factor(c(rep("JUL",3),rep("WT",3))))
deseq2.dds <- DESeqDataSetFromMatrix(aa,deseq2.coldata,design = ~ condition)
deseq2.dds$condition <- relevel(deseq2.dds$condition, ref = 'WT')
deseq2.dds <- DESeq(deseq2.dds)
deseq2.res <- results(deseq2.dds)
resultsNames(deseq2.res)

deseq2.df<-as.data.frame(deseq2.res)
deseq2.df005 <- deseq2.df[which(deseq2.df$padj < 0.05 & abs(deseq2.df$log2FoldChange) >= 2),]
deseq2.df001 <- deseq2.df[which(deseq2.df$padj < 0.01 & abs(deseq2.df$log2FoldChange) >= 2),]
write.csv(deseq2.df005,file="deseq2df005.csv")
