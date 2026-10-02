#!/usr/bin/env Rscript
ROOT <- Sys.getenv("PATHWAY_SCORE_ROOT", unset = getwd())
setwd(ROOT)
capture.output(sessionInfo(), file = "sessionInfo_current.txt")
cat("Saved sessionInfo_current.txt\n")
