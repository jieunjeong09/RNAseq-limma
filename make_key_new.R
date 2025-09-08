#!/usr/bin/env Rscript
# Purpose: build a sample key for GSE263611 from the series matrix
# Output: GSE263611_key.RDS with columns: Sample, Accession, Treatment, Assay, CellLine
# this version makes full use of vector operation and identifies needed lines by patterns

suppressPackageStartupMessages({
  library(stringr)
  library(tidyverse)
})

# Cache: skip if the key already exists
if (file.exists("GSE263611_key.RDS")) quit(save = "no")  # end quietly in script mode

# Download series matrix to a temp file (cleaned up automatically)
tmp <- tempfile(fileext = ".gz")
on.exit(unlink(tmp), add = TRUE)

GEOpref <- file.path(
  "https://ftp.ncbi.nlm.nih.gov/geo/series",
  "GSE263nnn/GSE263611/matrix"
)
url <- file.path(GEOpref, "GSE263611_series_matrix.txt.gz")
download.file(url, destfile = tmp, mode = "wb", quiet = TRUE)

# Read all lines; remove quotes to simplify parsing
con <- gzfile(tmp, open = "rt")
Line <- readLines(con, warn = FALSE)
close(con)
Line <- gsub('"', "", Line, fixed = TRUE)

# Helper: extract a tab-delimited GEO "row" by tag, return vector of fields
get_values <- function(regEx) {
  x <- grep(regEx, Line, value = TRUE)
  if (length(x) != 1) stop("regEx not found or duplicated: ", regEx)
  vals <- str_split_1(x, "\t")
  return(vals[-1])
}

# 1) Accession / Sample short ID from accession
#    accession have the form: GSMxxxxxxx 
Accession <- get_values("!Sample_geo_accession")
# acc[1] is the tag itself; 2..N are per-sample values
Sample <- paste0("S", substr(Accession, 9, 11))     # short ID rule for THIS project

# 2) Treatment from only characteristics line containing "treatment:"
# entry examples: "treatment: IFNg"       "treatment: Untreated Control"
Treatment <- get_values("^!Sample_characteristics_ch1\t.*treatment")
Treatment <- sub("^treatment: ?", "", Treatment) |> sub("^Untreated ?", "", x = _)

# 3) Assay and CellLine from entries like  "ATAC-seq for UDID148"  "RNA-seq for UDID006"
Assay <- get_values("^!Sample_description")
CellLine <- sub("^.* for ","",Assay)
Assay <- sub(" for .*$","",Assay)

# Build data frame (no rownames; keep Sample as a proper column)
DF <- tibble( Sample, Accession, Treatment, Assay, CellLine)

# Save RDS (single-object cache)
saveRDS(DF, "GSE263611_key.RDS")
