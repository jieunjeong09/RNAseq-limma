library(GEOquery, quietly = TRUE)
library(readr, quietly = TRUE)
library(dplyr, quietly = TRUE)
Gse <- getGEO("GSE263611", GSEMatrix = FALSE)
Sample <- GSMList(Gse)
Url <- sapply(Sample, function(x) Meta(x)$supplementary_file)
Url <- Url[grep("genes",Url)] # change if you want ATAC-seq
n_samples <- length(Url)
samples <- sapply(Url, function(x) sub("^.*GSM","GSM",x))
coltypes = c(rep("c",2),rep("n",5))
for (i in 1:n_samples) {
  download.file(Url[[i]], destfile = 'sample.gz', mode = "wb")
  df <- read_tsv('sample.gz', col_types = coltypes, progress = FALSE)
  file.remove('sample.gz')
  df <- df %>% select(gene_id,expected_count) %>%
    mutate(gene_id = substr(gene_id,1,15)) %>%
    distinct(gene_id,.keep_all = TRUE)
  if (i == 1) {
    Ensembl <- df$gene_id
    counts <- data.frame(df$expected_count)
  } else
    counts[[i]] <- df$expected_count
}
counts[[n_samples+1]] <- Ensembl
colnames(counts) <- c(paste0("S",substr(samples,9,10)),"Ensembl")
counts <- counts[apply(counts[,1:8], 1, sd) != 0, ] # remove genes with 0 std
saveRDS(counts,'GSE263611_counts.RDS')
