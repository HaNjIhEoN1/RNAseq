BiocManager::install('biomaRt')
library(biomaRt)
ensembl = useMart(biomart = "plants_mart",host = "https://plants.ensembl.org",dataset="osativa_eg_gene")
get_go <- getBM(attributes = c("ensembl_gene_id","go_id"), mart = ensembl)
get_go <- get_go[get_go$go_id != '',]
geneID2GO <- by(get_go$go_id, get_go$ensembl_gene_id,function(x) as.character(x))
all.genes <- sort(unique(as.character(get_go$ensembl_gene_id)))




