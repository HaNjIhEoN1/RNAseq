library(dplyr) # data wrangling
library(ggplot2) # plotting
library(edgeR) # rna-seq
library(tximport) # importing kalisto transcript counts to geneLevels
library(readr) # Fast readr of files.
library(rhdf5)
install.packages("statmod")
library(statmod)

mr <- read.csv('rt/CKvsAS/metadata.txt', sep = '\t')
files <- paste('rt/CKvsAS', list.files(path = 'rt/CKvsAS',pattern = 'abundance.h5',recursive = TRUE),sep = '/')
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

y <- edgeR::DGEList(cts)
y <- scaleOffset(y, normMat)

design <- model.matrix(~Library.Name, data = mr)
keep <- filterByExpr(y,design)
y <- y[keep,]

cpms <- cpm(y, offset = y$offset, log = FALSE)

plotMDS(y)

y1 <- estimateDisp(y, design = design, robust = TRUE)
y1$common.dispersion

plotBCV(y1)

fit <- glmFit(y1, design = design)
lrt <- glmLRT(fit)
topTags(lrt)


fit1 <- glmQLFit(y1, design = design, robust = TRUE)
plotQLDisp(fit1)

qlf <- glmQLFTest(fit1, coef = 2)
topTags(qlf)
FDR <- p.adjust(qlf$table$PValue, method = 'BH')
sum(FDR < 0.05)
qlf <- glmQLFTest(fit1)
topTags(qlf)

ql <- as.data.frame(qlf)
qlfilter <- ql[which(abs(ql$logFC) >= 2),]


plotMD(lrt)
abline(h=c(-1,1),col='blue')

summary(decideTests(lrt))

plotMD(qlf)
abline(h=c(-1,1),col='blue')

summary(decideTests(qlf))

-----------------------------------------------------------------------------------------------------

getwd()
setwd("C:/Users/jtige/OneDrive/R")
library(dplyr) # data wrangling
library(ggplot2) # plotting
library(edgeR) # rna-seq
library(tximport) # importing kalisto transcript counts to geneLevels
library(readr) # Fast readr of files.
library(rhdf5)
library(statmod)

mr <- read.csv('rt/CKvsAS/metadata.txt', sep = '\t')
files <- paste('rt/CKvsAS', list.files(path = 'rt/CKvsAS',pattern = 'abundance.h5',recursive = TRUE),sep = '/')
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

group = mr$Library.Name
y <- DGEList(cts, group = group)
y <- scaleOffset(y, normMat)
y$samples$group <- relevel(y$samples$group, ref = 'CK')
design <- model.matrix(~group)
keep <- filterByExpr(y,design)

y <- y[keep,]
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
y$common.dispersion

et<- exactTest(y)
topTags(et)
summary(decideTests(et))
etd <- as.data.frame(et)
etf <- etd[which(abs(etd$logFC)>=2),]
write.csv(etf,file="rt/CKvsAS/REALAS.csv")


plotMD(et)
abline(h=c(-1,1),col='blue')

