# kegg

setwd("C:/Users/jtige/OneDrive/R")
a <- read.csv('rt/cdddedger.csv')
a$gene

library(riceidconverter)

msu <- list()
for(i in a$gene){
  b <- RiceIDConvert(i, 'RAP', toType = 'MSU')
  msu <- append(msu,b)
}

repeat{
  msu['RAP'] = NULL
  if(length(msu['RAP'])==0) break
}

capture.output(msu, file = "my_list.txt")

BiocManager::install("clusterProfiler", version = "3.16", force = TRUE)

# KEGG in R
library(clusterProfiler)
search_kegg_organism('dosa', by='kegg_code')
rice <- search_kegg_organism('Oryza sativa japonica', by = 'scientific_name')
dim(rice)
head(rice)
genelist = a$gene
R.utils::setOption("clusterProfiler.download.method","auto")

mkk2 <- enrichKEGG(genelist, organism = 'dosa', keyType = 'kegg' ,pvalueCutoff = 0.05)
head(mkk2)
browseKEGG(mkk2, "dosa00480")

library(DOSE)
library(enrichplot)

dotplot(mkk2, showCategory = 30, font.size = 10 )

b <- read.csv('rt/asddedger.csv')
b$gene
genelistb = b$gene
mkk2b <- enrichKEGG(genelistb, organism = 'dosa', keyType = 'kegg' ,pvalueCutoff = 0.05)
dotplot(mkk2b, showCategory = 30, font.size = 10 )

