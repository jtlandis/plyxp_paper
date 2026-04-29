# plyxp_paper

The files relating to plyxp application note

## Getting Started

This repository requires a some packages from Bioconductor and CRAN to get started.
Below is an example of how you can install the required packages to run either [code_chunks.R](./code_chunks.R) or [comparison.R](./comparison.R)

```r

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Bioconductor packages
if (!require("plyxp", quietly = TRUE)) {
   BiocManager::install("plyxp")
}
if (!require("airway", quietly = TRUE)) {
   BiocManager::install("airway")
}
if (!require("tidySummarizedExperiment", quietly = TRUE)) {
   BiocManager::install("tidySummarizedExperiment")
}

# CRAN package for profiling
if (!require("bench", quietly = TRUE)) {
  install.packages("bench")
}
if (!require("purrr", quietly = TRUE)) {
  install.packages("purrr")
}
if (!require("vctrs", quietly = TRUE)) {
  install.packages("vctrs")
}
```
