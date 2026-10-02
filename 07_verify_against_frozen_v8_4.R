#!/usr/bin/env Rscript
## Verify a fresh R rerun against the frozen v8.4 manuscript results.
## This compares the full CSV result tables used in the final manuscript and
## checks representative outputs from the simulation/application RDS objects.

ROOT <- Sys.getenv("PATHWAY_SCORE_ROOT", unset = getwd())
setwd(ROOT)

CODEDIR <- Sys.getenv("PATHWAY_SCORE_CODE_DIR", unset = getwd())
REFDIR <- file.path(CODEDIR, "reference_results_v8_4")
OUTDIR <- "verification"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(REFDIR)) stop("Missing reference directory: ", REFDIR)

ABS_TOL <- as.numeric(Sys.getenv("VERIFY_ABS_TOL", "1e-8"))
REL_TOL <- as.numeric(Sys.getenv("VERIFY_REL_TOL", "1e-7"))

report <- list()
add_result <- function(group, item, status, max_abs_diff = NA_real_, note = "") {
  report[[length(report) + 1L]] <<- data.frame(
    group = group, item = item, status = status,
    max_abs_diff = max_abs_diff, note = note,
    stringsAsFactors = FALSE
  )
}

compare_csv <- function(current, reference, keys = character()) {
  if (!file.exists(current)) {
    add_result("CSV", current, "FAIL", note = "fresh-run file missing")
    return(invisible(FALSE))
  }
  if (!file.exists(reference)) {
    add_result("CSV", current, "FAIL", note = "reference file missing")
    return(invisible(FALSE))
  }
  a <- read.csv(current, stringsAsFactors = FALSE, check.names = FALSE)
  b <- read.csv(reference, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(names(a), names(b))) {
    add_result("CSV", current, "FAIL", note = "column names/order differ")
    return(invisible(FALSE))
  }
  if (nrow(a) != nrow(b)) {
    add_result("CSV", current, "FAIL", note = sprintf("row count %d vs %d", nrow(a), nrow(b)))
    return(invisible(FALSE))
  }
  if (length(keys)) {
    miss <- setdiff(keys, names(a))
    if (length(miss)) stop("Bad verification key(s) for ", current, ": ", paste(miss, collapse=", "))
    oa <- do.call(order, c(a[keys], list(na.last = TRUE)))
    ob <- do.call(order, c(b[keys], list(na.last = TRUE)))
    a <- a[oa, , drop = FALSE]
    b <- b[ob, , drop = FALSE]
    rownames(a) <- rownames(b) <- NULL
  }
  maxdiff <- 0
  problems <- character()
  for (nm in names(a)) {
    xa <- a[[nm]]; xb <- b[[nm]]
    if (is.numeric(xa) && is.numeric(xb)) {
      if (!identical(is.na(xa), is.na(xb))) {
        problems <- c(problems, paste0(nm, ": NA pattern differs")); next
      }
      ok <- is.finite(xa) & is.finite(xb)
      if (any(xor(is.finite(xa), is.finite(xb)), na.rm = TRUE)) {
        problems <- c(problems, paste0(nm, ": finite/nonfinite pattern differs")); next
      }
      if (any(ok)) {
        d <- abs(xa[ok] - xb[ok])
        scale <- pmax(1, abs(xb[ok]))
        tol <- ABS_TOL + REL_TOL * scale
        ## Bootstrap quantile endpoints are Monte Carlo summaries.  Across R builds/
        ## RNG implementations, a fixed seed can still yield tiny endpoint changes.
        ## Keep strict tolerance for all deterministic quantities, but allow 0.002
        ## absolute tolerance for the two descriptive median-bootstrap CI endpoints.
        if (basename(current) == "genome_wide_audit_benchmark.csv" &&
            nm %in% c("gene_boot_median_lo", "gene_boot_median_hi")) {
          tol <- pmax(tol, 0.002)
        }
        maxdiff <- max(maxdiff, max(d), na.rm = TRUE)
        if (any(d > tol)) problems <- c(problems, paste0(nm, ": numeric difference above tolerance"))
      }
    } else {
      ca <- ifelse(is.na(xa), "<NA>", as.character(xa))
      cb <- ifelse(is.na(xb), "<NA>", as.character(xb))
      if (!identical(ca, cb)) problems <- c(problems, paste0(nm, ": values differ"))
    }
  }
  if (length(problems)) {
    add_result("CSV", current, "FAIL", maxdiff, paste(unique(problems), collapse = "; "))
    return(invisible(FALSE))
  }
  add_result("CSV", current, "PASS", maxdiff, "matches frozen v8.4 reference")
  invisible(TRUE)
}

checks <- list(
  list("gap_inference_results_R.csv", "gap_inference_results_R.csv", c("n","p","scenario")),
  list("tcga_lung/genome_wide_audit.csv", "genome_wide_audit.csv", "gene"),
  list("tcga_lung/genome_wide_audit_summary.csv", "genome_wide_audit_summary.csv", "score"),
  list("tcga_lung/genome_wide_audit_sensitivity.csv", "genome_wide_audit_sensitivity.csv", "score"),
  list("tcga_lung/genome_wide_audit_diagnostics.csv", "genome_wide_audit_diagnostics.csv", "metric"),
  list("tcga_lung/genome_wide_audit_benchmark.csv", "genome_wide_audit_benchmark.csv", "score"),
  list("tcga_lung/genome_wide_audit_benchmark_by_p.csv", "genome_wide_audit_benchmark_by_p.csv", "p_bin"),
  list("tcga_lung/genome_wide_pooled_summary.csv", "genome_wide_pooled_summary.csv", "score"),
  list("tcga_lung/genome_wide_pooled_by_p.csv", "genome_wide_pooled_by_p.csv", "p_bin"),
  list("tcga_lung/cross_screened_audit.csv", "cross_screened_audit.csv", c("direction","gene")),
  list("tcga_lung/cross_screened_audit_summary.csv", "cross_screened_audit_summary.csv", c("direction","selection_rule","score")),
  list("tcga_lung/cross_screened_audit_mechanism.csv", "cross_screened_audit_mechanism.csv", c("direction","selection_rule")),
  list("tcga_lung/cross_screened_audit_sign_summary.csv", "cross_screened_audit_sign_summary.csv", c("direction","selection_rule","min_same_sign_fraction")),
  list("tcga_lung/cross_screened_pooled_summary.csv", "cross_screened_pooled_summary.csv", c("direction","selection_rule","score")),
  list("tcga_lung/cross_screened_screening_counts.csv", "cross_screened_screening_counts.csv", "direction"),
  list("tcga_lung/cross_screened_screening_overlap.csv", "cross_screened_screening_overlap.csv", "top_k")
)
for (z in checks) compare_csv(z[[1]], file.path(REFDIR, z[[2]]), z[[3]])

## ---- Representative checks for 02_simulation_and_application.R -------------
check_scalar <- function(group, item, got, expected, tol) {
  if (!is.finite(got) || abs(got - expected) > tol) {
    add_result(group, item, "FAIL", abs(got - expected),
               sprintf("got %.10g; expected %.10g; tolerance %.3g", got, expected, tol))
  } else {
    add_result(group, item, "PASS", abs(got - expected),
               sprintf("got %.10g; expected %.10g", got, expected))
  }
}

simfile <- "tcga_lung/simulation_results.rds"
if (!file.exists(simfile)) {
  add_result("RDS", simfile, "FAIL", note = "missing")
} else {
  sim <- readRDS(simfile)
  req <- c("truth_pop","summp","chk","spec","sens","e2","p4")
  if (!all(req %in% names(sim))) add_result("RDS", simfile, "FAIL", note="missing expected objects")
  else {
    add_result("RDS", simfile, "PASS", note="all expected objects present")
    chk <- as.data.frame(sim$chk)
    getcg <- function(rep, col) chk[chk$rep == rep, col][1]
    check_scalar("Simulation profile", "PC1 m", getcg("PC1","m"), 0.283, 0.0015)
    check_scalar("Simulation profile", "PC1 gamma", getcg("PC1","gamma"), 0.481, 0.0015)
    check_scalar("Simulation profile", "PC20 gamma", getcg("PC20","gamma"), 0.236, 0.0015)
    check_scalar("Simulation profile", "PC50 gamma", getcg("PC50","gamma"), 0.218, 0.0015)
    check_scalar("Simulation profile", "PC100 gamma", getcg("PC100","gamma"), 0.208, 0.0015)
    check_scalar("Simulation profile", "sign gamma", getcg("sign","gamma"), 0.412, 0.0015)
    check_scalar("Simulation profile", "means gamma", getcg("means","gamma"), 0.327, 0.0015)
    check_scalar("Simulation profile", "alpha_wt gamma", getcg("alpha_wt","gamma"), 0.503, 0.0015)

    summp <- as.data.frame(sim$summp)
    gets <- function(rep, col) summp[summp$rep == rep, col][1]
    check_scalar("Sparse simulation", "sign median relative gap", gets("sign","med_rel_gap"), 0.249, 0.0015)
    check_scalar("Sparse simulation", "PC1 pct >50%", gets("PC1","pct_gt50"), 34.4, 0.051)
    check_scalar("Sparse simulation", "PC10 pct >50%", gets("PC10","pct_gt50"), 5.4, 0.051)
    check_scalar("Sparse simulation", "PC100 false mediation", gets("PC100","false_med"), 0.5, 0.051)

    e2 <- as.data.frame(sim$e2)
    e2s <- aggregate(cbind(m4, g2_4) ~ rep, e2, function(x) c(mean=mean(x), neg=mean(x<0)))
    ## Directly recreate the manuscript summaries from raw e2.
    for (rp in c("PC20","PC50","PC100")) {
      d <- e2[e2$rep == rp,]
      g4 <- sqrt(max(mean(d$g2_4),0))
      expg <- c(PC20=0.220, PC50=0.156, PC100=0.000)[rp]
      check_scalar("Four-fold estimator", paste(rp,"gamma_hat"), g4, expg, 0.0015)
    }
    d100 <- e2[e2$rep == "PC100",]
    check_scalar("Four-fold estimator", "PC100 negative gamma2 percent", 100*mean(d100$g2_4<0), 65, 0.051)
  }
}

appfile <- "tcga_lung/application_results.rds"
if (!file.exists(appfile)) {
  add_result("RDS", appfile, "FAIL", note = "missing")
} else {
  app <- readRDS(appfile)
  req <- c("te_tab","prof_app","rob","nqo1")
  if (!all(req %in% names(app))) add_result("RDS", appfile, "FAIL", note="missing expected objects")
  else {
    add_result("RDS", appfile, "PASS", note="all expected objects present")
    te <- as.data.frame(app$te_tab)
    gette <- function(g, col) te[te$gene == g,col][1]
    check_scalar("Application", "AHRR total effect", gette("AHRR","TE"), 0.62495, 0.00006)
    check_scalar("Application", "UCHL1 total effect", gette("UCHL1","TE"), 0.85375, 0.00006)
    check_scalar("Application", "NQO1 total effect", gette("NQO1","TE"), 0.50972, 0.00006)

    rob <- as.data.frame(app$rob)
    getrob <- function(g, spec, rep, col) rob[rob$gene==g & rob$spec==spec & rob$rep==rep,col][1]
    check_scalar("Application", "AHRR full indirect effect", getrob("AHRR","main","IE_M","est"), 0.3996, 0.00015)
    check_scalar("Application", "AHRR gene-mean indirect effect", getrob("AHRR","main","gene_mean","est"), -0.0408, 0.00015)
    check_scalar("Application", "AHRR gene-mean gap", getrob("AHRR","main","gene_mean","gap"), 0.44042, 0.00006)
    check_scalar("Application", "UCHL1 full indirect effect", getrob("UCHL1","main","IE_M","est"), 0.6840, 0.00015)
    check_scalar("Application", "UCHL1 gene-mean indirect effect", getrob("UCHL1","main","gene_mean","est"), 0.5137, 0.00015)

    nq <- as.data.frame(app$nqo1)
    getnq <- function(rep, col) nq[nq$rep==rep,col][1]
    check_scalar("Application", "NQO1 full indirect effect", getnq("IE_M","est"), 0.07489, 0.00006)
    check_scalar("Application", "NQO1 gene-mean indirect effect", getnq("gene_mean","est"), -0.00859, 0.00006)
    check_scalar("Application", "NQO1 PC1 refit indirect effect", getnq("PC1_refit","est"), -0.01871, 0.00006)
  }
}

## Figures should exist and be nonempty; PDF bytes can differ because of metadata.
for (f in c("figures/figure1_profiles.pdf","figures/figure2_application.pdf","figures/figure3_pooled_benchmark.pdf")) {
  if (file.exists(f) && file.info(f)$size > 1000) add_result("Figure", f, "PASS", note="exists and nonempty")
  else add_result("Figure", f, "FAIL", note="missing or unexpectedly small")
}

repdf <- do.call(rbind, report)
write.csv(repdf, file.path(OUTDIR,"verification_report.csv"), row.names = FALSE)

nfail <- sum(repdf$status == "FAIL")
lines <- c(
  "Pathway-score full R rerun verification",
  paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("R:", R.version.string),
  paste("Checks:", nrow(repdf)),
  paste("PASS:", sum(repdf$status == "PASS")),
  paste("FAIL:", nfail),
  "",
  capture.output(print(repdf, row.names = FALSE))
)
writeLines(lines, file.path(OUTDIR,"verification_report.txt"))
cat(paste(lines[1:6], collapse="\n"), "\n")
if (nfail > 0) stop("Verification FAILED: ", nfail, " checks failed. See verification/verification_report.txt")
cat("All verification checks passed.\n")
