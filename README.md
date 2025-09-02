# RNA-seq workflow: Automated GEO Download, limma Analysis, and PCA

**1** File download and data table creation

`make_Hugo_identifiers.R` downloads HGNC data file and makes `HUGO.RDS` of the table of gene `Type` and identifiers: `Ensembl`, `Symbol` and `Entrez`.

Different programs require different identifiers.  For example, RSEM uses Ensembl which annotates more genes.  This increases the chance of identifying genes that lack Entrez IDs or HGNC-approved symbols.

Processing requires meta data key that can use information from soft file or series matrix file, the latter has a format that is easier to use.  In either case, one has to modify this script for other data series.

`make_key.R` downloads the series matrix and makes `GSE263611_key.RDS` of table with rownames that identify samples and columns that provide GEO `Accession`, `Treatment`, `CellLine` and `Assay`.  For nicer tables, sample with accessions GSM8195226 to  GSM8195233 get sample names `S26` to `S33`.  Such short names are practical for GEO data samples, but when there are more samples, we may use more terminal digits from Accessions

Some data series include count data in data series matrix, but not in this series.

`make_counts.R` identifies and downloads 8 sample files and creates `GSE263611_counts.RDS` with columns `S26` to `S33` and `Ensembl` with raw counts estimates made by RSEM.

**2** The analytical workflow

`mrna.jieun.R` finds direrential genes using edgeR, voom and limma, finds significant biological processes.
This vignette is meant to present completely initial processing, and thus it has chunks with downloads and making data tables.


A. Reads RDS files (if not computed already)

B. Makes two version of counts file by changing identifiers to Symbol and to Entrez (required by ReactomePA package for biological pathways.

C. Shows the heatmap for 100 most variable genes.

D. Shows results of PCA analysis in two ways: standard and vectors from Control to IFNg treatment in the same cell line.

E. Finds differential genes using edgeR, voom and limma workflow.

F. Shows gene set analysis for GO Biological Process.

**3** Comparing two designs

The design that includes CellLine leads to identification of several times more significant biological processes,
it produces much smaller adjusted p-values and thus more genes sets with adjusted p-value below 0.05.
This shows importance of including the confounding factor in the analysis.

**4** Further plans: more steps in gene set analysis



