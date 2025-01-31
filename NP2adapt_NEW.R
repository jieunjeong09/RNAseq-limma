setwd("/Users/jieun/Work/R_Jan_25")
library(tidyverse, quietly=T)
library(janitor, quietly=T)
library(limma, quietly=T)
library(edgeR, quietly=T)
library(readr, quietly=T)

# collecting input tables #

# table of Symbol, Entrez, Ensembl and locus type from hgnc_complete_set.txt
hugo <- readRDS('HUGO.RDS')
# on GSE263611:
#   samples of with Sample (short name), Accession, Treatment, Assay, # CellType
key <- readRDS('GSE263611_key.RDS') %>% filter(Assay == 'RNA-seq')
#   fragment counts for each sample and gene (Sample are column names)
counts <- readRDS('GSE263611_counts.RDS')

types_of_interest <- c("gene with protein product", "pseudogene", 
  "RNA, long non-coding")
counts <- counts %>% left_join(hugo, by="Ensembl") %>%
  filter(Type %in% types_of_interest) %>% select(-Type)

countsS <- counts %>% filter(Symbol != "") %>% distinct(Symbol, .keep_all = T) %>%
  select(-Entrez, -Ensembl)
countsS <- as.matrix(column_to_rownames(countsS, "Symbol"))

countsM <- counts %>% filter(Entrez != "") %>% distinct(Entrez, .keep_all = T) %>%
 select(-Symbol, -Ensembl) %>% column_to_rownames("Entrez")
countsM <- as.matrix(countsM)

### Define the condition or group
treatment <- factor(key[ , 'Treatment'])
cell_line <- factor(key[ , 'CellLine'])

design <- model.matrix(~ 0 + treatment)
# making DGEList object, TMM normalisation
deg <- DGEList(countsS)
#  Voom normalization
v <- voom(countsS, design)

### Linear model fitting
fit <- lmFit(v, design)

### Control group setup
contrast <- makeContrasts(treatmentIFNg - treatmentControl, levels = design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)

### Check the result
# https://support.bioconductor.org/p/53177/
results <- topTable(fit, adjust.method = "BH", coef = 1, number = Inf)

head(results)

### Visualization of DEG

## Volcano plot

# Volcano plot using volcanoplot()
volcanoplot(fit, highlight = 10)

fit$genes <- rownames(countsN)

volcanoplot(fit, highlight = 10, names = fit$genes)


## Volcano plot using ggplot2
## https://sbc.shef.ac.uk/rnaseq-r-online/session3.nb.html
## max.overlaps: https://ggrepel.slowkow.com/articles/examples.html
top_genes <- rownames(results[-log10(results$adj.P.Val) > 3 & abs(results$logFC) > 2,])

library(ggrepel)

results %>%
  mutate(Significant = -log10(adj.P.Val) > 3 & abs(logFC) > 2) %>% 
  mutate(Label = 
           ifelse(rownames(results) %in% top_genes, rownames(results), "")) %>%
  ggplot(aes(x = logFC, y = -log10(adj.P.Val),
             col = Significant,
             label = Label)) + 
  geom_point() + 
  geom_text_repel(col = "blue", max.overlaps = 15)

#### EXTRA 9.GSEA 
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

geneList <- results[order(results$logFC,decreasing = T),'logFC']

names(geneList) <- 
  rownames(results[order(results$logFC,decreasing = T),])

gsea_result <- gseGO(geneList = geneList,
                     OrgDb = org.Hs.eg.db,
                     key = "SYMBOL",
                     ont = "BP",
)

gsea_result %>% as.data.frame %>% head()
str(gsea_result)


### GSEA visualization: 
# https://sbc.shef.ac.uk/rnaseq-r-online/session3.nb.html
# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html


# Ridgeplot
ridgeplot(gsea_result)


# Upsetplot
enrichplot::upsetplot(gsea_result)


# GSEA plot
gseaplot(gsea_result, geneSetID = "GO:0140546")
gseaplot(gsea_result, geneSetID = "GO:0051607")
gseaplot(gsea_result, geneSetID = "GO:0009615")
gseaplot(gsea_result, geneSetID = "GO:0019079")
gseaplot(gsea_result, geneSetID = "GO:0048525")
gseaplot(gsea_result, geneSetID = "GO:0045069")


# Dotplot
dotplot(gsea_result, showCategory = 10)


### Reactome
geneListN <- geneList

names(geneListN) <- genenames$NCBI_Gene_ID[match(names(geneListN),
                                             genenames$Approved_symbol)]

reactome_result <- gsePathway(geneList = geneListN, 
                              pvalueCutoff = 0.05)

# View the results
head(as.data.frame(reactome_result))

# Ridgeplot
ridgeplot(reactome_result)


# Upsetplot
enrichplot::upsetplot(reactome_result)


# GSEA plot
gseaplot(reactome_result, geneSetID = "R-HSA-913531")
gseaplot(reactome_result, geneSetID = "R-HSA-909733")
gseaplot(reactome_result, geneSetID = "R-HSA-877300")
gseaplot(reactome_result, geneSetID = "R-HSA-1169410")
gseaplot(reactome_result, geneSetID = "R-HSA-1236974")
gseaplot(reactome_result, geneSetID = "R-HSA-983170")


# Dotplot
dotplot(reactome_result, showCategory = 10)

