library(dplyr) # data wrangling
library(ggplot2) # plotting
library(DESeq2) # rna-seq
library(tximport) # importing kalisto transcript counts to geneLevels
library(readr) # Fast readr of files.
library(rhdf5)

setwd("C:/Users/jtige/OneDrive/R")
mr <- read.csv('rt/metadata.txt', sep = '\t')
files <- paste('rt/kallisto', list.files(path = 'rt/kallisto',pattern = 'abundance.h5',recursive = TRUE),sep = '/')
names(files) <- mr$Run
tx2gene <- read.csv('rt/rice.cdna.csv')
txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
mr = mr %>% mutate(Library.Name = as.factor(Library.Name))
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, mr, ~Library.Name)
dds$Library.Name <- relevel(dds$Library.Name, ref = 'CK')
dds <- DESeq(dds)
res <- results(dds)
results(dds)

as <- results(dds, name="Library.Name_As_vs_CK")
cd <- results(dds, name="Library.Name_Cd_vs_CK")
asd <- as.data.frame(as)
cdd <- as.data.frame(cd)
asf <- asd[which(asd$padj < 0.01 & abs(asd$log2FoldChange) >= 2),]
cdf <- cdd[which(cdd$padj < 0.01 & abs(cdd$log2FoldChange) >= 2),]

write.csv(asd,file="rt/deseqas.csv")
write.csv(asd,file="rt/deseqcd.csv")
