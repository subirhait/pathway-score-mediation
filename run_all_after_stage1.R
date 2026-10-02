#!/usr/bin/env Rscript
## Run the full R analysis after tcga_lung/luad_stage1.rds has been created.
this_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
  error = function(e) NA_character_
)
CODEDIR <- if (!is.na(this_file)) dirname(this_file) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
ROOT <- normalizePath(Sys.getenv("PATHWAY_SCORE_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
Sys.setenv(PATHWAY_SCORE_CODE_DIR = CODEDIR)
source_code <- function(fname) source(file.path(CODEDIR, fname), chdir = FALSE)
setwd(ROOT)
source_code("01_validate_stage1_input.R")
source_code("02_simulation_and_application.R")
source_code("03_gap_inference_simulation.R")
source_code("05_genome_wide_audit.R")
source_code("06_cross_screened_audit.R")
source_code("04_make_figures.R")
source_code("99_session_info.R")
