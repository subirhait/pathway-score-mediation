#!/usr/bin/env Rscript
## Install packages required by the pathway-score reproducibility workflow.
cran <- c("data.table", "survival", "pseudo", "glmnet", "dplyr", "UCSCXenaTools")
need <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("limma", quietly = TRUE))
  BiocManager::install("limma", ask = FALSE, update = FALSE)

cat("Package check complete.\n")
