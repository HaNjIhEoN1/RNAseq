#Differntial expression analysis
dds <- DESeq(dds)
res <- results(dds)

#specified the coefficient or contrast 
##case 1
res <- results(dds, name="$factor_$factor1_vs_$factor2")
##case 2
res <- results(dds, contrast=c('$factor', '$factor1', '$factor2'))

#Expoloring and exporting reults
##MA-plot
plotMA(res,ylim=c(-2,2))

##Plot counts
plotCounts(dds, gene = '$genename', intgroup='$factor')

##Heatmap
ntd <- normTransform(dds)
library('pheatmap')
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c('$factor1','$factor2')])
pheatmap(assay(ntd)[select,],cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE, annotation_col=df)

vsd <- vst(dds, blind = FALSE)
library('RColorBrewer')
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Sample_Type)
colnames(sampleDistMatrix)  <- NULL
pheatmap::pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists)

##compare edgeR, DESeq2 degs

edgeR <- read.csv('rt/CKvsAS/edgerresult.csv')
deseq2 <- read.csv('rt/CKvsAS/as.csv')
install.packages('VennDiagram')

install.packages('gplots')
library(gplots)
e <- as.data.frame(c(edgeR['gene']))
d <- as.data.frame(c(deseq2['gene']))
data <- list(edgeR = e, DESeq2 = d)
venn(data)

library(VennDiagram)

data<- list(edgeR = edgeR$gene, DESeq2 = deseq2$gene)
venn.diagram(data, fill = c(3,2),
             alpha = c(0.5, 0.5),
             lty = c(1,1), filename = 'test.tiff'
             )
