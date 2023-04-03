setwd("C:/Users/jtige/OneDrive/R")
a <- read.csv('rt/CKvsAS/convertid.csv')
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
