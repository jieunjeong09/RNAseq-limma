# RNAseq workflow with download from NCBI GEO, limma, PCA etc.

**1** File downloads and making data tables

`make_Hugo_identifiers.R` downloads HGNC data file and makes `HUGO.RDS` of the table of gene `Type` and identifiers: `Ensembl`, `Symbol` and `Entrez`.

Different programs require different identifiers, mapping program RSEM used Ensembl that annotates more genes, and there is a potential of finding roles of genes without identifiers in Entrez or HGNC approved symbols.

Processing requires meta data key that can use information from soft file or series matrix file, the latter has a format that is easier to use.  In either case, one has to modify this script for other data series.

`make_key.R` downloads the series matrix and makes `GSE263611_key.RDS` of table with rownames that identify samples and columns that provide GEO `Accession`, `Treatment`, `CellLine` and `Assay`.  For nicer tables, sample with accessions GSM8195226 to  GSM8195233 get sample names `S26` to `S33`.  Such short names are practical for GEO data samples, but when there are more samples, we may use more terminal digits from Accessions

Some data series include count data in data series matrix, but not in this series.

`make_counts.R` identifies and downloads 8 sample files and creates `GSE263611_counts.RDS` with columns `S26` to `S33` and `Ensembl` with raw counts estimates made by RSEM.

**2** The analytical workflow

`mrna.jieun.R` finds direrential genes using edgeR, voom and limma, finds significant biological processes and pathways and other findings.

A. Reads RDS files

B. Makes two version of counts file by changing identifiers to Symbol and to Entrez (required by ReactomePA package for biological pathways.

C. Shows the heatmap for 100 most variable genes.

D. Shows results of PCA analysis in two ways: standard and vectors from Control to IFNg treatment in the same cell line.

E. Finds differential genes using edgeR, voom and limma workflow.

F. Shows gene set analysis for GO Biological Precess and Reactome pathways.

**3** Comparing two designs and two ranking methods

`limma` workflow may use simpler approach `design <- model.matrix(~ 0 + treatment)` or include the confounding factor `design <- model.matrix(~ 0 + treatment + cell_line)`.  Subsequently, when we form gene list for gsea analysis, we can rank it by fold change or by adjusted p-value (of course, other approaches exist too).  

To simplify R code, we run two scripts:

A. `compare_limma_1.R` reads RDS objects saved by `mrna.jieun.R` and saves two gsea objects.  It also visualizes differences between fold changes between two designs, which are small, and adjusted p-values of genes which are large.  Thus later we use adjusted p-values.

B. `compare_limma_2.R reads RDS objects save by `compare_limma_1.R and compare the results.
Comparison with potentially better identification of differential genes using `cell_line` as a confounding factor in "design", ie. `design <- model.matrix(~ 0 + treatment)` 
and `design <- model.matrix(~ 0 + treatment + cell_line)`

If we rank genes according to fold change, the differences in ranking of values were small and the differences in GO Biological Process.

But the values of adjusted p-values change dramatically, and so the results for GO Biological Process when we use this ranking.  Including `cell_line` nearly doubles the number of identified processes, and identifies roles of almost 5 times more genes.  The most dramatic example is WNT5A that does not occur in any process when `cell_line` is excluded from the design, and in 210 if it is included.

This shows importance of including the confounding factor in the analysis.



