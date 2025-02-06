setwd('/Users/jieun/Work/Git_test/RNAseq_limma')
library(dplyr)
library(tidyverse)

readRDS('Ri1p.RDS',Ri1p)
df <- data.frame(a = Ri1p$core_enrichment)
Ri1p_gene_counts <- df %>% separate_rows(a,sep='/') %>%
  count(a,name='occurences') %>% arrange(desc(occurences))
  
readRDS('Ri2p.RDS',Ri2p)
df <- data.frame(a = Ri2p$core_enrichment)
Ri2p_gene_counts <- df %>% separate_rows(a,sep='/') %>%
  count(a,name='occurences') %>% arrange(desc(occurences))
# the second design identifies many more significant pathways
dim(Ri1p)
dim(Ri2p)
# the second design identifies many more differential
# genes in significant pathways
dim(Ri1p_gene_counts)
dim(Ri2p_gene_counts)
# WNT5A is prominent in the second design, undedected in the first
filter(Ri2p_gene_counts,a=="WNT5A")
filter(Ri2p_gene_counts,a=="WNT5A")