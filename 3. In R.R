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
files <- paste("kallisto_quantification",
               list.files(path = "kallisto_quantification",pattern = "abundance.tsv", recursive = TRUE),
               sep = "/")

# designate names for each file
names(files) <- mr$SampleName # SampleName or Run


tx2gene <- read_csv("ref.cdna.csv")
