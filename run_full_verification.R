#!/usr/bin/env Rscript
## One-command clean verification run after tcga_lung/luad_stage1.rds exists.
## The code can live outside PATHWAY_SCORE_ROOT; scripts are sourced from CODEDIR.

this_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
  error = function(e) NA_character_
)
CODEDIR <- if (!is.na(this_file)) dirname(this_file) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
ROOT <- Sys.getenv("PATHWAY_SCORE_ROOT", unset = getwd())
ROOT <- normalizePath(ROOT, winslash = "/", mustWork = TRUE)
Sys.setenv(PATHWAY_SCORE_CODE_DIR = CODEDIR)

source_code <- function(fname) {
  f <- file.path(CODEDIR, fname)
  if (!file.exists(f)) stop("Missing verification script: ", f)
  source(f, chdir = FALSE)
}

setwd(ROOT)
dir.create("verification", showWarnings = FALSE, recursive = TRUE)
logfile <- file.path("verification", "full_rerun_console.log")
zz <- file(logfile, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")
on.exit({
  while (sink.number(type = "message") > 2) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(zz)
}, add = TRUE)

cat("FULL R VERIFICATION RUN\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Project root:", ROOT, "\n")
cat("Code folder:", CODEDIR, "\n\n")

source_code("01_validate_stage1_input.R")
source_code("02_simulation_and_application.R")
source_code("03_gap_inference_simulation.R")

## The frozen v8.4 analyses used 5,000 gene-bootstrap replicates for the
## full-data genome-wide benchmark and 2,000 for the cross-screen audit.
## Keep these two Monte Carlo layers separate so a single environment setting
## cannot silently change one of the archived confidence intervals.
genome_boot <- Sys.getenv("GENOME_AUDIT_GENE_BOOT", "5000")
cross_boot  <- Sys.getenv("CROSS_AUDIT_GENE_BOOT", Sys.getenv("AUDIT_GENE_BOOT", "2000"))

Sys.setenv(
  AUDIT_SPLITS = Sys.getenv("AUDIT_SPLITS", "25"),
  AUDIT_GENE_BOOT = genome_boot
)
source_code("05_genome_wide_audit.R")

Sys.setenv(
  AUDIT_INNER_SPLITS = Sys.getenv("AUDIT_INNER_SPLITS", "25"),
  AUDIT_GENE_BOOT = cross_boot,
  AUDIT_TOPK = Sys.getenv("AUDIT_TOPK", "25,50,100,250")
)
source_code("06_cross_screened_audit.R")
source_code("04_make_figures.R")
source_code("99_session_info.R")
source_code("07_verify_against_frozen_v8_4.R")
cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
