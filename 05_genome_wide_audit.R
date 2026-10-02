#!/usr/bin/env Rscript

## ============================================================================
## Genome-wide outcome-free aggregation audit
## "When Does a Pathway Score Preserve a Mediated Effect?"
##
## Purpose
##   For smoking-responsive genes in TCGA LUAD, estimate the exposure-information
##   loss omega_B^2 = 1 - kappa_B/kappa for common annotation-based summaries,
##   without using gene expression or any downstream outcome.
##
## Required input
##   tcga_lung/luad_stage1.rds
## created by the data-preparation workflow used by 02_simulation_and_application.R.
##
## Primary output
##   tcga_lung/genome_wide_audit.csv
##   tcga_lung/genome_wide_audit.rds
##   tcga_lung/genome_wide_audit_summary.csv
##
## Estimator
##   Two-fold cross-fitted Wishart-corrected estimators of kappa_B and kappa.
##   Repeated splits are averaged at the kappa level before forming omega_B^2.
##   omega_B^2 is intentionally NOT truncated to [0,1]. Small excursions outside
##   that range are possible in finite samples, especially when kappa is weak.
##
## Notes
##   * Primary scores are means of the original methylation M-values, matching the
##     application analysis. Set STANDARDIZED_MEANS=1 for a sensitivity analysis
##     based on means of standardized CpGs.
##   * No outcome is used anywhere in this script.
##   * Base R + data.table only.
## ============================================================================

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages(library(data.table))

SEED <- as.integer(Sys.getenv("AUDIT_SEED", "2026"))
NSPLIT <- as.integer(Sys.getenv("AUDIT_SPLITS", "25"))  # main manuscript run
MIN_CPG <- as.integer(Sys.getenv("AUDIT_MIN_CPG", "5"))
MAX_CPG <- as.integer(Sys.getenv("AUDIT_MAX_CPG", "150"))
FDR_CUT <- as.numeric(Sys.getenv("AUDIT_FDR", "0.05"))
KAPPA_MIN <- as.numeric(Sys.getenv("AUDIT_KAPPA_MIN", "0"))
GENE_BOOT <- as.integer(Sys.getenv("AUDIT_GENE_BOOT", "5000"))
STANDARDIZED_MEANS <- as.integer(Sys.getenv("STANDARDIZED_MEANS", "0")) == 1L

infile <- "tcga_lung/luad_stage1.rds"
if (!file.exists(infile)) {
  stop("Missing ", infile,
       ". Run the stage-1 data-preparation workflow first, then rerun this script.")
}

s2 <- readRDS(infile)
Mval <- s2$Mval
ph <- s2$ph
ann <- s2$ann
tt <- s2$alpha
rm(s2)

if (!all(c("A", "age", "sex", "stage4", "purity_cpe") %in% names(ph)))
  stop("ph is missing one or more required columns: A, age, sex, stage4, purity_cpe")
if (!"P.Value" %in% names(tt)) stop("s2$alpha must contain P.Value")
if (is.null(rownames(tt))) stop("s2$alpha must have CpG row names")
if (is.null(rownames(ann))) stop("s2$ann must have CpG row names")
if (!all(c("UCSC_RefGene_Name", "UCSC_RefGene_Group") %in% names(ann)))
  stop("Annotation must contain UCSC_RefGene_Name and UCSC_RefGene_Group")

## Align samples explicitly when identifiers are available.
if ("sample" %in% names(ph) && !is.null(colnames(Mval))) {
  jj <- match(ph$sample, colnames(Mval))
  if (anyNA(jj)) stop("Some ph$sample IDs are absent from Mval columns")
  Mval <- Mval[, jj, drop = FALSE]
}
if (ncol(Mval) != nrow(ph)) stop("Mval columns and ph rows are not aligned")

ph$stage4 <- droplevels(ph$stage4)
A <- ph$A
if (anyNA(A)) stop("Exposure A contains missing values")
D0 <- model.matrix(~ A + age + sex + stage4 + purity_cpe, data = ph)
acoef <- which(colnames(D0) == "A")
if (length(acoef) != 1L) stop("Could not identify the A column in the design matrix")

## ---- CpG-to-gene annotation -------------------------------------------------
clean_split <- function(x) {
  x <- ifelse(is.na(x), "", x)
  strsplit(x, ";", fixed = TRUE)
}
gene_list <- clean_split(ann$UCSC_RefGene_Name)
grp_list <- clean_split(ann$UCSC_RefGene_Group)
if (!all(lengths(gene_list) == lengths(grp_list))) {
  stop("Gene-name and gene-group annotation lengths differ for at least one CpG")
}

gl <- data.table(
  cpg  = rep(rownames(ann), lengths(gene_list)),
  gene = unlist(gene_list, use.names = FALSE),
  grp  = unlist(grp_list,  use.names = FALSE)
)
gl <- gl[gene != "" & !is.na(gene)]
gl[, prom := grepl("TSS1500|TSS200|5'UTR|1stExon", grp)]
## One CpG can carry multiple annotation strings for the same gene.
gl <- gl[, .(prom = any(prom)), by = .(gene, cpg)]

pmap <- data.table(cpg = rownames(tt), pval = tt$P.Value)
gl <- pmap[gl, on = "cpg"]

simes <- function(p) {
  p <- p[is.finite(p)]
  if (!length(p)) return(NA_real_)
  p <- sort(p)
  min(1, min(length(p) * p / seq_along(p)))
}

gs <- gl[, .(N = .N, p_simes = simes(pval)), by = gene]
gs <- gs[N >= MIN_CPG & N <= MAX_CPG & is.finite(p_simes)]
gs[, fdr := p.adjust(p_simes, "BH")]
genes <- gs[fdr < FDR_CUT, gene]
cat("Eligible smoking-responsive genes:", length(genes), "\n")
if (!length(genes)) stop("No genes pass the requested CpG-count/FDR filters")

## ---- Repeated stratified two-fold splits -----------------------------------
make_split <- function(seed) {
  set.seed(seed)
  fold <- integer(length(A))
  for (lev in unique(A)) {
    w <- which(A == lev)
    w <- sample(w)
    fold[w] <- rep(1:2, length.out = length(w))
  }
  fold
}
folds <- lapply(seq_len(NSPLIT), function(s) make_split(SEED + s - 1L))

## Wishart-corrected cross-fitted kappa estimator for a fixed B.
kappa_cf <- function(a1, a2, S, df, B) {
  B <- as.matrix(B)
  k <- ncol(B)
  if (df <= k + 1L) return(NA_real_)
  G <- crossprod(B, S %*% B)
  z1 <- crossprod(B, a1)
  z2 <- crossprod(B, a2)
  ans <- tryCatch(
    (df - k - 1) / df * drop(crossprod(z1, solve(G, z2))),
    error = function(e) NA_real_
  )
  ans
}

fit_gene_split <- function(Mg, fold, Blist) {
  p <- ncol(Mg)
  fit <- vector("list", 2L)
  RSS <- matrix(0, p, p)
  df <- 0L
  for (h in 1:2) {
    ii <- which(fold == h)
    Xh <- D0[ii, , drop = FALSE]
    Mh <- Mg[ii, , drop = FALSE]
    Q <- qr(Xh)
    cf <- qr.coef(Q, Mh)
    ah <- cf[acoef, ]
    Rh <- qr.resid(Q, Mh)
    dh <- length(ii) - Q$rank
    if (dh <= 0) return(NULL)
    RSS <- RSS + crossprod(Rh)
    df <- df + dh
    fit[[h]] <- ah
  }
  if (df <= p + 1L) return(NULL)
  S <- RSS / df
  kF <- kappa_cf(fit[[1]], fit[[2]], S, df, diag(p))
  vals <- vapply(Blist, function(B)
    kappa_cf(fit[[1]], fit[[2]], S, df, B), numeric(1))
  c(kappa = kF, vals)
}

mean_abs_cor <- function(M) {
  if (ncol(M) < 2L) return(NA_real_)
  C <- suppressWarnings(cor(M, use = "pairwise.complete.obs"))
  mean(abs(C[upper.tri(C)]), na.rm = TRUE)
}

audit_gene <- function(g) {
  cp <- gl[gene == g]
  cp <- cp[match(unique(cpg), cpg)]
  jj <- match(cp$cpg, rownames(Mval))
  ok <- !is.na(jj)
  cp <- cp[ok]
  jj <- jj[ok]
  if (length(jj) < MIN_CPG || length(jj) > MAX_CPG) return(NULL)

  Mg <- t(Mval[jj, , drop = FALSE])
  keep <- apply(Mg, 2, function(x) all(is.finite(x)) && sd(x) > 0)
  Mg <- Mg[, keep, drop = FALSE]
  cp <- cp[keep]
  p <- ncol(Mg)
  if (p < MIN_CPG) return(NULL)

  ## Primary: fixed averages on the original M-value scale.
  ## Optional sensitivity: average standardized CpGs.
  if (STANDARDIZED_MEANS) {
    s <- apply(Mg, 2, sd)
    Mg <- sweep(Mg, 2, colMeans(Mg), "-")
    Mg <- sweep(Mg, 2, s, "/")
  }

  Blist <- list(mean = matrix(1 / p, p, 1L))
  if (any(cp$prom))
    Blist$promoter <- matrix(as.numeric(cp$prom) / sum(cp$prom), p, 1L)
  if (any(cp$prom) && any(!cp$prom))
    Blist$prom_nonprom <- cbind(
      as.numeric(cp$prom) / sum(cp$prom),
      as.numeric(!cp$prom) / sum(!cp$prom)
    )

  ee <- lapply(folds, function(f) fit_gene_split(Mg, f, Blist))
  ee <- ee[!vapply(ee, is.null, logical(1))]
  if (!length(ee)) return(NULL)
  EE <- do.call(rbind, ee)
  kF <- mean(EE[, "kappa"], na.rm = TRUE)
  if (!is.finite(kF)) return(NULL)

  getkb <- function(nm) mean(EE[, nm], na.rm = TRUE)
  kb_mean <- getkb("mean")
  kb_prom <- if ("promoter" %in% colnames(EE)) getkb("promoter") else NA_real_
  kb_pnp <- if ("prom_nonprom" %in% colnames(EE)) getkb("prom_nonprom") else NA_real_

  data.table(
    gene = g,
    p = p,
    n_prom = sum(cp$prom),
    p_simes = gs[gene == g, p_simes][1],
    fdr = gs[gene == g, fdr][1],
    kappa = kF,
    kappa_split_sd = sd(EE[, "kappa"], na.rm = TRUE),
    mean_abs_cor = mean_abs_cor(Mg),
    kappaB_mean = kb_mean,
    omega2_mean = 1 - kb_mean / kF,
    kappaB_prom = kb_prom,
    omega2_prom = 1 - kb_prom / kF,
    kappaB_prom_nonprom = kb_pnp,
    omega2_prom_nonprom = 1 - kb_pnp / kF,
    n_splits = nrow(EE),
    standardized_means = STANDARDIZED_MEANS
  )
}

cat("Running outcome-free audit...\n")
audit <- rbindlist(lapply(genes, audit_gene), fill = TRUE)
setorder(audit, -omega2_mean)

if (!nrow(audit)) stop("No gene produced a usable audit estimate")

## Main summary: retain untruncated omega^2 values.
summary_one <- function(x) c(
  median = median(x, na.rm = TRUE),
  q25 = unname(quantile(x, 0.25, na.rm = TRUE)),
  q75 = unname(quantile(x, 0.75, na.rm = TRUE)),
  pct_gt_0.5 = 100 * mean(x > 0.5, na.rm = TRUE),
  pct_gt_0.8 = 100 * mean(x > 0.8, na.rm = TRUE)
)

Smean <- summary_one(audit$omega2_mean)
Sprom <- summary_one(audit$omega2_prom)
Spnp  <- summary_one(audit$omega2_prom_nonprom)

summ <- data.table(
  score = c("Gene mean", "Promoter mean", "Promoter + non-promoter"),
  genes = c(sum(is.finite(audit$omega2_mean)),
            sum(is.finite(audit$omega2_prom)),
            sum(is.finite(audit$omega2_prom_nonprom))),
  median_omega2 = c(Smean["median"], Sprom["median"], Spnp["median"]),
  q25_omega2 = c(Smean["q25"], Sprom["q25"], Spnp["q25"]),
  q75_omega2 = c(Smean["q75"], Sprom["q75"], Spnp["q75"]),
  pct_gt50 = c(Smean["pct_gt_0.5"], Sprom["pct_gt_0.5"], Spnp["pct_gt_0.5"]),
  pct_gt80 = c(Smean["pct_gt_0.8"], Sprom["pct_gt_0.8"], Spnp["pct_gt_0.8"])
)

cat("\nGenome-wide outcome-free audit summary:\n")
print(summ, digits = 4)

## Fixed-threshold sensitivity summaries used in the manuscript.
sens_thresholds <- c(0.2, 0.5)
sens_list <- lapply(sens_thresholds, function(thr) {
  aa <- audit[kappa >= thr]
  if (!nrow(aa)) return(NULL)
  sm <- summary_one(aa$omega2_mean)
  sp <- summary_one(aa$omega2_prom)
  sn <- summary_one(aa$omega2_prom_nonprom)
  data.table(
    kappa_min = thr,
    score = c("Gene mean", "Promoter mean", "Promoter + non-promoter"),
    genes = c(sum(is.finite(aa$omega2_mean)), sum(is.finite(aa$omega2_prom)),
              sum(is.finite(aa$omega2_prom_nonprom))),
    median_omega2 = c(sm["median"], sp["median"], sn["median"]),
    q25_omega2 = c(sm["q25"], sp["q25"], sn["q25"]),
    q75_omega2 = c(sm["q75"], sp["q75"], sn["q75"]),
    pct_gt50 = c(sm["pct_gt_0.5"], sp["pct_gt_0.5"], sn["pct_gt_0.5"]),
    pct_gt80 = c(sm["pct_gt_0.8"], sp["pct_gt_0.8"], sn["pct_gt_0.8"])
  )
})
sens_tab <- rbindlist(sens_list, fill = TRUE)
## Preserve the frozen-v8.4 public column order. The earlier verification
## failure was only a column-order mismatch, not a numerical discrepancy.
setcolorder(sens_tab, c("score", "genes", "median_omega2", "q25_omega2",
                        "q75_omega2", "pct_gt50", "pct_gt80", "kappa_min"))

## Sensitivity restricted to stronger exposure information.
if (KAPPA_MIN > 0) {
  aa <- audit[kappa >= KAPPA_MIN]
  cat("\nSensitivity subset with kappa >=", KAPPA_MIN, ":", nrow(aa), "genes\n")
  if (nrow(aa)) {
    print(aa[, .(
      median_omega2_mean = median(omega2_mean, na.rm = TRUE),
      pct_mean_gt50 = 100 * mean(omega2_mean > 0.5, na.rm = TRUE),
      median_omega2_prom = median(omega2_prom, na.rm = TRUE),
      median_omega2_pnp = median(omega2_prom_nonprom, na.rm = TRUE)
    )], digits = 4)
  }
}

## Descriptive relation between mean-score loss, dimension, and correlation.
regdat <- audit[is.finite(omega2_mean) & is.finite(mean_abs_cor) & p > 0]
if (nrow(regdat) >= 10L) {
  cat("\nDescriptive regression: omega2_mean ~ log(p) + mean_abs_cor\n")
  print(summary(lm(omega2_mean ~ log(p) + mean_abs_cor, data = regdat)))
}

## Dimension-only random-subspace benchmark and descriptive gene bootstrap.
boot_med <- function(x, seed) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(c(lo = NA_real_, hi = NA_real_))
  set.seed(seed)
  z <- replicate(GENE_BOOT, median(sample(x, length(x), replace = TRUE)))
  unname(quantile(z, c(0.025, 0.975), na.rm = TRUE))
}

bench_def <- list(
  list(score = "Gene mean", omega = "omega2_mean", k = 1L),
  list(score = "Promoter mean", omega = "omega2_prom", k = 1L),
  list(score = "Promoter + non-promoter", omega = "omega2_prom_nonprom", k = 2L)
)
bench_tab <- rbindlist(lapply(seq_along(bench_def), function(j) {
  d <- bench_def[[j]]
  x <- audit[[d$omega]]
  rb <- 1 - d$k / audit$p
  ok <- is.finite(x) & is.finite(rb)
  rr <- (1 - x[ok]) / (1 - rb[ok])
  ci <- boot_med(x[ok], SEED + 50000L + j)
  data.table(
    score = d$score,
    genes = sum(ok),
    median_omega2 = median(x[ok]),
    gene_boot_median_lo = ci[1],
    gene_boot_median_hi = ci[2],
    median_random_subspace_loss = median(rb[ok]),
    pct_better_than_random = 100 * mean(x[ok] < rb[ok]),
    pct_worse_than_random = 100 * mean(x[ok] > rb[ok]),
    median_retained_information_ratio = median(rr),
    spearman_omega2_p = if (d$score == "Gene mean") cor(x[ok], audit$p[ok], method = "spearman") else NA_real_,
    spearman_omega2_mean_abs_cor = if (d$score == "Gene mean") cor(x[ok], audit$mean_abs_cor[ok], method = "spearman") else NA_real_
  )
}))

bins <- list(c(5L,10L), c(11L,25L), c(26L,60L), c(61L,150L))
bench_by_p <- rbindlist(lapply(bins, function(z) {
  aa <- audit[p >= z[1] & p <= z[2] & is.finite(omega2_mean)]
  rb <- 1 - 1 / aa$p
  data.table(
    p_bin = paste0(z[1], "-", z[2]),
    genes = nrow(aa),
    median_gene_mean_loss = median(aa$omega2_mean),
    median_random_subspace_loss = median(rb),
    pct_better_than_random = 100 * mean(aa$omega2_mean < rb),
    median_retained_ratio = median((1 - aa$omega2_mean) / (1 - rb))
  )
}))

## Finite-sample diagnostics: do not hide estimates outside [0,1].
cat("\nFinite-sample diagnostics:\n")
cat("omega2_mean outside [0,1]:",
    sum(is.finite(audit$omega2_mean) & (audit$omega2_mean < 0 | audit$omega2_mean > 1)), "\n")
cat("kappa <= 0:", sum(is.finite(audit$kappa) & audit$kappa <= 0), "\n")

outdir <- "tcga_lung"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
fwrite(audit, file.path(outdir, "genome_wide_audit.csv"))
saveRDS(audit, file.path(outdir, "genome_wide_audit.rds"))
fwrite(summ, file.path(outdir, "genome_wide_audit_summary.csv"))
fwrite(sens_tab, file.path(outdir, "genome_wide_audit_sensitivity.csv"))

diag <- data.table(
  metric = c("eligible_genes", "kappa_nonpositive", "mean_omega2_below0",
             "mean_omega2_above1", "mean_omega2_outside01",
             "prom_omega2_outside01", "pnp_omega2_outside01"),
  value = c(nrow(audit), sum(is.finite(audit$kappa) & audit$kappa <= 0),
            sum(is.finite(audit$omega2_mean) & audit$omega2_mean < 0),
            sum(is.finite(audit$omega2_mean) & audit$omega2_mean > 1),
            sum(is.finite(audit$omega2_mean) & (audit$omega2_mean < 0 | audit$omega2_mean > 1)),
            sum(is.finite(audit$omega2_prom) & (audit$omega2_prom < 0 | audit$omega2_prom > 1)),
            sum(is.finite(audit$omega2_prom_nonprom) &
                (audit$omega2_prom_nonprom < 0 | audit$omega2_prom_nonprom > 1)))
)
fwrite(diag, file.path(outdir, "genome_wide_audit_diagnostics.csv"))
fwrite(bench_tab, file.path(outdir, "genome_wide_audit_benchmark.csv"))
fwrite(bench_by_p, file.path(outdir, "genome_wide_audit_benchmark_by_p.csv"))

capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo_genome_wide_audit.txt"))

cat("\nSaved:\n",
    file.path(outdir, "genome_wide_audit.csv"), "\n",
    file.path(outdir, "genome_wide_audit.rds"), "\n",
    file.path(outdir, "genome_wide_audit_summary.csv"), "\n",
    file.path(outdir, "sessionInfo_genome_wide_audit.txt"), "\n", sep = "")
