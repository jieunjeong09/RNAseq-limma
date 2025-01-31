# RNAseq workflow with download from NCBI GEO, limma, PCA etc.

## File downloads and making data tables

### HGNC table for identifying gene type and identifiers: Ensembl, Symbol and Entrez

make_Hugo_identifiers.R downloads HGNC data file and makes HUGO.RDS of the table

### meta data table for samples from the study

make_key.R downloads the series matrix and makes GSE263611_key.RDS, a table with rownames that identify samples and columns that provide GEO access ID, treatment, cell line and assay.

### table of counts for each Ensembl gene identifier (rows) and each sample (columns)

make_counts.R identifies and downloads 8 sample files and creates GSE263611_counts.RDS.

## The analytical workflow

