setwd("/Users/jieun/Work/Git_test/RNAseq_limma")
library(dplyr,quietly = T)
library(tidyr,quietly = T)
library(tidyverse,quietly = T)
L1 = readRDS("Limma_1.RDS")
L2 = readRDS("Limma_2.RDS")
R1 = readRDS("GSEA_1.RDS")
R2 = readRDS("GSEA_2.RDS")

Li1 <- L1 %>% rownames_to_column() %>% select(-t,-B,-P.Value) %>%
  mutate(L10.adj.P = -log10(adj.P.Val)) %>% select(-adj.P.Val)
Li2 <- L2 %>% rownames_to_column() %>% select(-t,-B,-P.Value) %>%
  mutate(L10.adj.P = -log10(adj.P.Val)) %>% select(-adj.P.Val)
LiComb <- Li1 %>% inner_join(Li2,by="rowname") %>% mutate(AveExpr=AveExpr.x) %>%
  select(-AveExpr.x,-AveExpr.y) %>% column_to_rownames("rowname")
saveRDS(LiComb,'LiComb.RDS')
# while the formula for fold change is different in two designs, the 
# results is minor
ggplot(LiComb,aes(x=logFC.x,logFC.y)) + geom_point(size=.5) +
  geom_abline(intercept=0,slope=1,color="orange")
# the adjusted p-values are dramatically different for many genes
ggplot(LiComb,aes(x=L10.adj.P.x,y=L10.adj.P.y)) + geom_point(size=.5) +
  geom_abline(intercept=0,slope=1,color="orange")

#### GSEA
library(clusterProfiler,quietly = T)
library(org.Hs.eg.db,quietly = T)

geneL <- LiComb$L10.adj.P.x %>% setNames(rownames(LiComb))
geneL <- sort(geneL,decreasing=T)
Ri1p <- gseGO(geneList = geneL, OrgDb = org.Hs.eg.db,
                               key = "SYMBOL", ont = "BP", scoreType = "pos")
saveRDS(Ri1p,"Ri1p.RDS")

geneL <- LiComb$L10.adj.P.y %>% setNames(rownames(LiComb))
geneL <- sort(geneL,decreasing=T)
Ri2p <- gseGO(geneList = geneL, OrgDb = org.Hs.eg.db,
              key = "SYMBOL", ont = "BP", scoreType = "pos")
saveRDS(Ri1p,"Ri2p.RDS")



