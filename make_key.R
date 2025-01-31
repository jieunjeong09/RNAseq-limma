library(stringr)
library(tidyverse)
if (file.exists('GSE263611_key.RDS')) # avoid unnecessary downloads
  quit()
GEOpref = file.path('https://ftp.ncbi.nlm.nih.gov/geo/series',
	'GSE263nnn/GSE263611/matrix')
Url = file.path(GEOpref,'GSE263611_series_matrix.txt.gz')
download.file(Url, destfile = 'matrix.gz', mode = "wb")
con <- gzfile('matrix.gz', "rt")
Line <- readLines(con)
close(con)
file.remove('matrix.gz')

Line <- gsub('"','',Line)
Accession <- character()
Sample <- character()
Treatment <- character()
Assay <- character()
CellLine <- character()
B <- character() # for extracting Assay and CellLine
E1 <- str_split_1(Line[30],'\t')
E2 <- str_split_1(Line[39],'\t')
E3 <- str_split_1(Line[47],'\t')

for (i in 2:17) {
  Accession[i-1] <- E1[i]
  Sample[i-1] <- paste0("S",substr(E1[i],9,11))
  t <- sub('treatment: ','',E2[i])
  Treatment[i-1] <- sub('Untreated ','',t)
  B <- str_split_1(E3[i], ' ')
  Assay[i-1] <- B[1]
  CellLine[i-1] <- B[3]
}
DF <- data.frame(Sample,Accession,Treatment,Assay,CellLine)
column_to_rownames(DF, "Sample")
saveRDS(DF,'GSE263611_key.RDS')