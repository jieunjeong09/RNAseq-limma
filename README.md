# RNAseq workflow with download from NCBI GEO, limma, PCA etc.

**1** File downloads and making data tables

`make_Hugo_identifiers.R` downloads HGNC data file and makes `HUGO.RDS` of the table of gene `Type` and identifiers: `Ensembl`, `Symbol` and `Entrez`.

Different programs require different identifiers, mapping program RSEM used Ensembl that annotates more genes, and there is a potential of finding roles of genes without identifiers in Entrez or HGNC approved symbols.

Processing requires meta data key that can use information from soft file or series matrix file, the latter has a format that is easier to use.  In either case, one has to modify this script for other data series.

`make_key.R` downloads the series matrix and makes `GSE263611_key.RDS` of table with rownames that identify samples and columns that provide GEO `Accession`, `Treatment`, `CellLine` and `Assay`.  For nicer tables, sample with accessions GSM8195226 to  GSM8195233 get sample names `S26` to `S33`.  Such short names are practical for GEO data samples, but when there are more samples, we may use more terminal digits from Accessions

Some data series include count data in data series matrix, but not in this series.

`make_counts.R` identifies and downloads 8 sample files and creates `GSE263611_counts.RDS` with columns `S26` to `S33` and `Ensembl` with raw counts estimates made by RSEM.

**2** The analytical workflow

'mrna.kieun.R' finds direrential genes using edgeR, voom and limma, finds significant biological processes and pathways and other findings.

A. Reads RDS files

B. Makes two version of counts file by changing identifiers to Symbol and to Entrez (required by ReactomePA package for biological pathways.

C. Shows the heatmap for 100 most variable genes.

D. Shows results of PCA analysis in two ways: standard and vectors from Control to IFNg treatment in the same cell line.

E. Finds differential genes using edgeR, voom and limma workflow.

F. Shows gene set analysis for GO Biological Precess and Reactome pathways.

**3** Plan to add on February 2

Comparison with potentially better identification of differential genes using CellLine as a confounding factor in "design",
it does not change most significant genes sets but it changes the set of differential genes etc.



