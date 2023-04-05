getwd()
setwd("C:/Users/jtige/OneDrive/R")
library(dplyr) # data wrangling
library(ggplot2) # plotting
library(edgeR) # rna-seq
library(tximport) # importing kalisto transcript counts to geneLevels
library(readr) # Fast readr of files.
library(rhdf5)
library(statmod)

mr <- read.csv('rt/metadata.txt', sep = '\t')
files <- paste('rt/kallisto', list.files(path = 'rt/kallisto',pattern = 'abundance.h5',recursive = TRUE),sep = '/')
names(files) <- mr$Run
txi.kallisto <- tximport(files, type = 'kallisto', txOut = TRUE)
head(txi.kallisto$counts)
cts <- txi.kallisto$counts
normMat <- txi.kallisto$length
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

eff.lib <- calcNormFactors(normCts)*colSums(normCts)

normMat <- sweep(normMat,2,eff.lib,"*")
normMat <- log(normMat)

group <- mr$Library.Name
y <- DGEList(cts, group = factor(group))
keep <- filterByExpr(y = y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
design <- model.matrix(~0+group, data = y$samples)
colnames(design) <- levels(y$samples$group)
fit <- glmQLFit(y,design)

As <- glmQLFTest(fit, contrast = c(1,0,-1))
topTags(As)

Cd <- glmQLFTest(fit, contrast = c(0,1,-1))
topTags(Cd)

Asd <- as.data.frame(As)
Asf <- Asd[which(abs(Asd$logFC)>=2),]
write.csv(Asf,file="rt/edgeras.csv")

Cdd <- as.data.frame(Cd)
Cdf <- Cdd[which(abs(Cdd$logFC)>=2),]
write.csv(Cdf,file="rt/edgercd.csv")

as2 <- topTags(As,n = Inf)
keep <- as2$table$FDR < 0.01& abs(as2$table$logFC) >=2
asdeg <- as2$table[keep,]
asdd <- as.data.frame(asdeg)
write.csv(asdd,file="rt/asddedger.csv")

cd2 <- topTags(Cd,n = Inf)
keep <- cd2$table$FDR < 0.01& abs(cd2$table$logFC) >=2
cddeg <- cd2$table[keep,]
cddd <- as.data.frame(cddeg)
write.csv(cddd,file="rt/cdddedger.csv")
