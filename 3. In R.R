library(dplyr) # data wrangling
library(ggplot2) # plotting
library(DESeq2) # rna-seq
library(edgeR) # rna-seq
library(tximport) # importing kalisto transcript counts to geneLevels
library(readr) # Fast readr of files.
library(rhdf5) # read/convert kalisto output files.  

#If tximport and rhdf5 are not installed you can install them now.
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("tximport")
BiocManager::install("rhdf5")
