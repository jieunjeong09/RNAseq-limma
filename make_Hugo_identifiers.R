if (file.exists("HUGO.RDS")) quit() # avoid unnecessary downloads
library(readr)
TT <- readLines('https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt')
Symbol <- character() # approved symbols
Entrez <- character() # Entrez gene id
Ensembl <- character() # Ensembl gene idT
Type <- character() # gene type
n = 0
for (L in TT[2:length(TT)]) {
  n = n+1
  E <- strsplit(L,'\t')[[1]]
  Symbol[n] <- E[2]
  Type[n] <- E[5]
  Entrez[n] <- E[19]
  Ensembl[n] <- E[20]
}
H <- data.frame(Symbol,Entrez,Ensembl,Type)
saveRDS(H,"HUGO.RDS")

