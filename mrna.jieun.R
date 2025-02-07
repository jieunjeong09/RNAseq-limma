setwd('/Users/jieun/Work/Git_test/RNAseq_limma')
library(tidyverse, quietly=T)
library(dplyr)
library(janitor, quietly=T)

# making input tables
hugo <- readRDS('HUGO.RDS')
samples <- readRDS('GSE263611_key.RDS') %>% filter(Assay == "RNA-seq")
types_of_interest <- c("gene with protein product", "pseudogene", "RNA,
                       long non-coding")
counts <- readRDS('GSE263611_counts.RDS') %>% left_join(hugo, by="Ensembl") %>%
  filter(Type %in% types_of_interest) %>% select(-Type)

countsAS <- counts %>% filter(Symbol != "") %>% distinct(Symbol, .keep_all = T) %>%
  select(-Entrez, -Ensembl)
countsAS <- as.matrix(column_to_rownames(countsAS, "Symbol"))

countsN <- counts %>% filter(Entrez != "") %>% distinct(Entrez, .keep_all = T) %>%
  select(-Symbol, -Ensembl)
countsN <- as.matrix(column_to_rownames(countsN, "Entrez"))

### Definitions for limma
treatment <- as.factor(samples$Treatment)
cell_line <- as.factor(samples$CellLine)
library(edgeR, quietly=T)
design <- model.matrix(~ 0 + treatment)
contrast <- makeContrasts(treatmentIFNg - treatmentControl, levels = design)

# Normalization for limma
dge <- DGEList(counts = countsAS) # TMM normalization, reformmating
dge <- calcNormFactors(dge)
v <- voom(countsAS, design, plot = T) # log normalization and adjustment
# checking the normalized distribution of gene expression
A <- as.data.frame(v$E)
A <- rownames_to_column(A,var='gene')
AL <- pivot_longer(A,-gene,names_to = "sample", values_to = "logExpr")
AL$sample <- as.factor(AL$sample)
p <- ggplot(AL,aes(x=sample,y=logExpr)) + geom_violin(trim=F,width=0.8)
p + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red",fun.args=list(mult=1))
### Clustering on v$E
library(SummarizedExperiment)

sample_info <- samples %>% select(-Assay, -Accession)
se <- SummarizedExperiment(assays = list(counts = v$E),
                           colData = sample_info)
colData(se)
colData = as.data.frame(colData(se))

library(pheatmap)

# Select the top 100 most variable genes among the samples
VG <- apply(v$E, 1, var)
selectedGenes <- names(VG[order(VG, decreasing = T)])[1:100]

p <- pheatmap(v$E[selectedGenes, ], 
         scale = 'row',
         show_rownames = FALSE,
         annotation_col = colData)
ggsave("Heat_heat_map_top_100.pdf",plot=p)
### PCA on v$E
library(factoextra) # for PCA

df_t <- t(v$E)

# Using only most variable genes
pr1 <- prcomp(df_t[,selectedGenes], scale. = TRUE)

p <- fviz_pca_ind(pr1,
             geom.ind = "point",
             col.ind = colData$Treatment,
             palette = "jco",
             addEllipses = TRUE,
             legene.title = "Treatment")
p
ggsave("PCA_of_100.pdf",plot=p)

pca_coords <- as.data.frame(pr1$x)
pca_coords$Sample <- rownames(pca_coords)
pca_coords$CellLine <- colData$CellLine
pca_coords$Treatment <- colData$Treatment  # Add treatment groups
p <- ggplot(pca_coords, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 4) +
  # stat_ellipse(type = "norm", level = 0.95) +  # 95% confidence ellipses
  geom_line(aes(group = CellLine), color = "gray") +
  xlim(-10, 10) + ylim(-10, 10) +
  theme_minimal()
p
ggsave("PCA_vectors.pdf",plot=p)
# Fitting by limma
fit <- lmFit(v, design) # linear

### Control treatment setup
fit2 <- contrasts.fit(fit, contrast)
fit3 <- eBayes(fit2)

### Check the result
results <- topTable(fit3, adjust.method = "BH", coef = 1, number = Inf)

head(results)
saveRDS(results,"Limma_1.RDS")

### Visualization of DEG

## Volcano plot

# Volcano plot using volcanoplot()
volcanoplot(fit3, highlight = 10, names = rownames(fit3$coefficients))


## Volcano plot using ggplot2
top_genes <-
  rownames(results[-log10(results$adj.P.Val) > 3 & abs(results$logFC) > 2,])

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
### limma with confounding factor: cell line
#### GSEA
library(clusterProfiler)

geneList <- results[order(results$logFC,decreasing = T),'logFC']

names(geneList) <-
  rownames(results[order(results$logFC,decreasing = T),])

gsea_result <- gseGO(geneList = geneList,
                     OrgDb = org.Hs.eg.db,
                     key = "SYMBOL",
                     ont = "BP")


gsea_result %>% as.data.frame %>% head()
saveRDS(gsea_result,"GSEA_1.RDS")


### GSEA visualization:
# Ridgeplot for GO-BP
ridgeplot(gsea_result, showCategory = 15)


# Upsetplot for GO-BP
enrichplot::upsetplot(gsea_result)

# GSEA plot for top GP-BP
gseaplot(gsea_result, geneSetID = "GO:0051607")
gseaplot(gsea_result, geneSetID = "GO:0140546")
gseaplot(gsea_result, geneSetID = "GO:0009615")
gseaplot(gsea_result, geneSetID = "GO:0050792")
gseaplot(gsea_result, geneSetID = "GO:0019079")
gseaplot(gsea_result, geneSetID = "GO:0048525")

# Dotplot for GO-BP
dotplot(gsea_result, showCategory = 15)


### Reactome
geneListN <- geneList

names(geneListN) <- hugo$Entrez[match(names(geneListN), 
                                              hugo$Symbol)]

reactome_result <- gsePathway(geneList = geneListN, pvalueCutoff = 0.05)

# View the results for Reactome
head(as.data.frame(reactome_result))

# Ridgeplot for Reactome
ridgeplot(reactome_result, showCategory = 15)

# Upsetplot for Reactome
enrichplot::upsetplot(reactome_result)

# GSEA plot for Top Reactome
gseaplot(reactome_result, geneSetID = "R-HSA-913531")
gseaplot(reactome_result, geneSetID = "R-HSA-877300")
gseaplot(reactome_result, geneSetID = "R-HSA-909733")
gseaplot(reactome_result, geneSetID = "R-HSA-1169410")
gseaplot(reactome_result, geneSetID = "R-HSA-9705671")
gseaplot(reactome_result, geneSetID = "R-HSA-449147")

# Dotplot for Reactome
dotplot(reactome_result, showCategory = 15)

### limma with a potentially better design


# alternative:
design <- model.matrix(~0 + treatment + cell_line)
contrast <- makeContrasts(treatmentIFNg - treatmentControl, levels = design)

# Normalization for limma
dge <- DGEList(counts = countsAS) # TMM normalization, reformmating
dge <- calcNormFactors(dge)
v <- voom(countsAS, design, plot = T)

# Fitting by limma
fit <- lmFit(v, design) # linear

### Control treatment setup
fit2 <- contrasts.fit(fit, contrast)
fit3 <- eBayes(fit2)

### Check the result
results <- topTable(fit3, adjust.method = "BH", coef = 1, number = Inf)

head(results)
saveRDS(results,"Limma_2.RDS")

geneList <- results[order(results$logFC,decreasing = T),'logFC']

names(geneList) <-
  rownames(results[order(results$logFC,decreasing = T),])

gsea_result <- gseGO(geneList = geneList,
                     OrgDb = org.Hs.eg.db,
                     key = "SYMBOL",
                     ont = "BP")


gsea_result %>% as.data.frame %>% head()
saveRDS(gsea_result,"GSEA_2.RDS")