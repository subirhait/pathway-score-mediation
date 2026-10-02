#!/usr/bin/env Rscript
## Validate the processed TCGA LUAD object used by all empirical analyses.
## This does NOT recreate the object from raw/public downloads.
## Expected path: tcga_lung/luad_stage1.rds

ROOT <- Sys.getenv("PATHWAY_SCORE_ROOT", unset = getwd())
setwd(ROOT)
file <- "tcga_lung/luad_stage1.rds"
if (!file.exists(file)) stop("Missing ", file,
  ". Copy the existing processed stage-1 object to this location.")

x <- readRDS(file)
req <- c("Mval", "ph", "ann", "alpha")
miss <- setdiff(req, names(x))
if (length(miss)) stop("Stage-1 object is missing: ", paste(miss, collapse = ", "))

Mval <- x$Mval; ph <- x$ph; ann <- x$ann; alpha <- x$alpha
if (!is.matrix(Mval) && !is.data.frame(Mval)) stop("Mval must be matrix/data.frame-like")
if (ncol(Mval) != nrow(ph)) stop("Mval columns do not match ph rows")
need_ph <- c("A", "age", "sex", "stage4", "purity_cpe", "OS.time", "OS")
miss_ph <- setdiff(need_ph, names(ph))
if (length(miss_ph)) stop("ph is missing: ", paste(miss_ph, collapse = ", "))
need_ann <- c("UCSC_RefGene_Name", "UCSC_RefGene_Group")
miss_ann <- setdiff(need_ann, names(ann))
if (length(miss_ann)) stop("ann is missing: ", paste(miss_ann, collapse = ", "))
if (is.null(rownames(Mval))) stop("Mval needs CpG row names")
if (is.null(rownames(ann))) stop("ann needs CpG row names")
if (is.null(rownames(alpha))) stop("alpha needs CpG row names")
if (!all(c("P.Value", "adj.P.Val") %in% names(alpha)))
  stop("alpha must contain P.Value and adj.P.Val")

cat("Stage-1 input validated successfully.\n")
cat("Samples:", nrow(ph), "\n")
cat("CpGs in Mval:", nrow(Mval), "\n")
cat("Smoking-responsive CpGs at adj.P.Val < 0.05:", sum(alpha$adj.P.Val < 0.05, na.rm = TRUE), "\n")
