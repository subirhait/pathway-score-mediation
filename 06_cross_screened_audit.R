#!/usr/bin/env Rscript

## ============================================================================
## Independent-screening sensitivity audit for pathway-score adequacy
## "When Does a Pathway Score Preserve a Mediated Effect?"
##
## Purpose
##   Address selection-on-the-same-sample in the genome-wide outcome-free audit.
##   One exposure-stratified half screens smoking-responsive genes; the other half
##   estimates kappa, kappa_B, and omega_B^2. The roles are then reversed.
##
##   The script additionally reports:
##     * the whitened-isotropic random-subspace benchmark 1 - k/p;
##     * a descriptive gene bootstrap interval for the median loss;
##     * Spearman associations with mediator dimension and within-gene correlation;
##     * sign consistency of exposure coefficients measured in the screening half;
##     * pooled ratio-of-sums loss with a gene-level bootstrap interval.
##
## Required input
##   tcga_lung/luad_stage1.rds
##
## Outputs
##   tcga_lung/cross_screened_audit.csv
##   tcga_lung/cross_screened_audit_summary.csv
##   tcga_lung/cross_screened_audit_mechanism.csv
##   tcga_lung/cross_screened_audit_sign_summary.csv
##   tcga_lung/cross_screened_pooled_summary.csv
##   tcga_lung/genome_wide_pooled_summary.csv  (if genome_wide_audit.csv is present)
##   tcga_lung/genome_wide_pooled_by_p.csv     (if genome_wide_audit.csv is present)
##   tcga_lung/sessionInfo_cross_screened_audit.txt
##
## Notes
##   * Screening uses limma on the screening half with the same adjustment set.
##   * Estimation uses repeated exposure-stratified two-fold splits inside the
##     independent estimation half and the Wishart-corrected cross-fitted kappa
##     estimator from Proposition 6.
##   * omega_B^2 is intentionally not truncated to [0,1].
##   * The 1-k/p benchmark is an expected loss for a uniformly random k-dimensional
##     subspace in whitened coordinates; it is a geometric reference, not a null
##     distribution for raw random weights under correlated mediators.
##   * The gene bootstrap treats the audited genes as the empirical units. It does
##     not replace participant-level resampling and is labelled descriptive.
## ============================================================================

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

SEED <- as.integer(Sys.getenv("AUDIT_SEED", "2026"))
INNER_SPLITS <- as.integer(Sys.getenv("AUDIT_INNER_SPLITS", "25"))
BOOT_REPS <- as.integer(Sys.getenv("AUDIT_GENE_BOOT", "2000"))
MIN_CPG <- as.integer(Sys.getenv("AUDIT_MIN_CPG", "5"))
MAX_CPG <- as.integer(Sys.getenv("AUDIT_MAX_CPG", "150"))
FDR_CUT <- as.numeric(Sys.getenv("AUDIT_FDR", "0.05"))
TOPK <- as.integer(strsplit(Sys.getenv("AUDIT_TOPK", "25,50,100,250"), ",", fixed = TRUE)[[1]])
TOPK <- sort(unique(TOPK[is.finite(TOPK) & TOPK > 0L]))
if (!length(TOPK)) stop("AUDIT_TOPK must contain at least one positive integer")
KMAX <- max(TOPK)
STANDARDIZED_MEANS <- as.integer(Sys.getenv("STANDARDIZED_MEANS", "0")) == 1L

infile <- "tcga_lung/luad_stage1.rds"
if (!file.exists(infile)) {
  stop("Missing ", infile, ". Run the stage-1 data-preparation workflow first.")
}

s2 <- readRDS(infile)
Mval <- s2$Mval
ph <- s2$ph
ann <- s2$ann
rm(s2)

if (!all(c("A", "age", "sex", "stage4", "purity_cpe") %in% names(ph))) {
  stop("ph is missing one or more required columns: A, age, sex, stage4, purity_cpe")
}
if (is.null(rownames(ann))) stop("Annotation must have CpG row names")
if (!all(c("UCSC_RefGene_Name", "UCSC_RefGene_Group") %in% names(ann))) {
  stop("Annotation must contain UCSC_RefGene_Name and UCSC_RefGene_Group")
}

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

full_rank_design <- function(X) {
  q <- qr(X)
  keep <- sort(q$pivot[seq_len(q$rank)])
  if (!(acoef %in% keep)) stop("Exposure coefficient A is not estimable in a split")
  list(X = X[, keep, drop = FALSE], acoef = match(acoef, keep))
}

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
  grp  = unlist(grp_list, use.names = FALSE)
)
gl <- gl[gene != "" & !is.na(gene)]
gl[, prom := grepl("TSS1500|TSS200|5'UTR|1stExon", grp)]
gl <- gl[, .(prom = any(prom)), by = .(gene, cpg)]

## Keep only annotated CpGs that are present in Mval.
annot_cpg <- intersect(unique(gl$cpg), rownames(Mval))
gl <- gl[cpg %in% annot_cpg]

simes <- function(p) {
  p <- p[is.finite(p)]
  if (!length(p)) return(NA_real_)
  p <- sort(p)
  min(1, min(length(p) * p / seq_along(p)))
}

## ---- Outer exposure-stratified split ---------------------------------------
make_stratified_half <- function(seed, idx = seq_along(A)) {
  set.seed(seed)
  half <- integer(length(idx))
  Ai <- A[idx]
  for (lev in unique(Ai)) {
    wloc <- which(Ai == lev)
    wloc <- sample(wloc)
    half[wloc] <- rep(1:2, length.out = length(wloc))
  }
  out <- integer(length(A))
  out[idx] <- half
  out
}

outer <- make_stratified_half(SEED)
if (!all(outer %in% c(1L, 2L))) stop("Outer split failed")

## ---- Independent screening --------------------------------------------------
screen_genes <- function(screen_idx) {
  ds <- full_rank_design(D0[screen_idx, , drop = FALSE])

  ## limma works directly with CpGs x samples. Nuisance columns that are
  ## unidentifiable in a half-sample are dropped, while A must remain estimable.
  fit <- lmFit(Mval[annot_cpg, screen_idx, drop = FALSE], design = ds$X)
  fit <- eBayes(fit)
  pval <- fit$p.value[, ds$acoef]

  pmap <- data.table(cpg = annot_cpg, pval = as.numeric(pval))
  gg <- pmap[gl, on = "cpg"]
  gs <- gg[, .(N = .N, p_simes = simes(pval)), by = gene]
  gs <- gs[N >= MIN_CPG & N <= MAX_CPG & is.finite(p_simes)]
  gs[, fdr := p.adjust(p_simes, "BH")]
  setorder(gs, p_simes)
  gs[, screen_rank := .I]
  gs[]
}

## ---- Wishart-corrected estimation on an independent half --------------------
kappa_cf <- function(a1, a2, S, df, B) {
  B <- as.matrix(B)
  k <- ncol(B)
  if (df <= k + 1L) return(NA_real_)
  G <- crossprod(B, S %*% B)
  z1 <- crossprod(B, a1)
  z2 <- crossprod(B, a2)
  tryCatch(
    (df - k - 1) / df * drop(crossprod(z1, solve(G, z2))),
    error = function(e) NA_real_
  )
}

make_inner_split <- function(est_idx, seed) {
  set.seed(seed)
  fold <- integer(length(est_idx))
  Ai <- A[est_idx]
  for (lev in unique(Ai)) {
    w <- which(Ai == lev)
    w <- sample(w)
    fold[w] <- rep(1:2, length.out = length(w))
  }
  fold
}

fit_gene_split <- function(Mg, D_est, fold, Blist) {
  p <- ncol(Mg)
  fit <- vector("list", 2L)
  RSS <- matrix(0, p, p)
  df <- 0L
  for (h in 1:2) {
    ii <- which(fold == h)
    ds <- full_rank_design(D_est[ii, , drop = FALSE])
    Xh <- ds$X
    Mh <- Mg[ii, , drop = FALSE]
    Q <- qr(Xh)
    cf <- qr.coef(Q, Mh)
    ah <- cf[ds$acoef, ]
    Rh <- qr.resid(Q, Mh)
    dh <- length(ii) - Q$rank
    if (dh <= 0L) return(NULL)
    RSS <- RSS + crossprod(Rh)
    df <- df + dh
    fit[[h]] <- ah
  }
  if (df <= p + 1L) return(NULL)
  S <- RSS / df
  kF <- kappa_cf(fit[[1]], fit[[2]], S, df, diag(p))
  vals <- vapply(Blist, function(B) kappa_cf(fit[[1]], fit[[2]], S, df, B), numeric(1))
  c(kappa = kF, vals)
}

mean_abs_cor <- function(M) {
  if (ncol(M) < 2L) return(NA_real_)
  C <- suppressWarnings(cor(M, use = "pairwise.complete.obs"))
  mean(abs(C[upper.tri(C)]), na.rm = TRUE)
}

estimate_gene <- function(g, gs_screen, screen_idx, est_idx, direction_id) {
  cp <- gl[gene == g]
  cp <- cp[match(unique(cpg), cpg)]
  jj <- match(cp$cpg, rownames(Mval))
  ok <- !is.na(jj)
  cp <- cp[ok]
  jj <- jj[ok]
  if (length(jj) < MIN_CPG || length(jj) > MAX_CPG) return(NULL)

  Mg <- t(Mval[jj, est_idx, drop = FALSE])
  Mg_screen <- t(Mval[jj, screen_idx, drop = FALSE])
  keep <- apply(Mg, 2, function(x) all(is.finite(x)) && sd(x) > 0)
  Mg <- Mg[, keep, drop = FALSE]
  Mg_screen <- Mg_screen[, keep, drop = FALSE]
  cp <- cp[keep]
  p <- ncol(Mg)
  if (p < MIN_CPG || p > MAX_CPG) return(NULL)

  if (STANDARDIZED_MEANS) {
    s <- apply(Mg, 2, sd)
    Mg <- sweep(Mg, 2, colMeans(Mg), "-")
    Mg <- sweep(Mg, 2, s, "/")
  }

  Blist <- list(mean = matrix(1 / p, p, 1L))
  if (any(cp$prom)) {
    Blist$promoter <- matrix(as.numeric(cp$prom) / sum(cp$prom), p, 1L)
  }
  if (any(cp$prom) && any(!cp$prom)) {
    Blist$prom_nonprom <- cbind(
      as.numeric(cp$prom) / sum(cp$prom),
      as.numeric(!cp$prom) / sum(!cp$prom)
    )
  }

  D_est <- D0[est_idx, , drop = FALSE]

  ## Sign consistency is defined in the SCREENING half only. It is therefore
  ## independent of the kappa/kappa_B estimates formed in the estimation half.
  D_screen <- D0[screen_idx, , drop = FALSE]
  ds_screen <- full_rank_design(D_screen)
  Q_screen <- qr(ds_screen$X)
  a_screen <- qr.coef(Q_screen, Mg_screen)[ds_screen$acoef, ]
  nz <- is.finite(a_screen) & a_screen != 0
  frac_same_sign <- if (any(nz)) {
    max(mean(a_screen[nz] > 0), mean(a_screen[nz] < 0))
  } else NA_real_

  inner <- lapply(seq_len(INNER_SPLITS), function(s) {
    make_inner_split(est_idx, SEED + 10000L * direction_id + s - 1L)
  })
  ee <- lapply(inner, function(f) fit_gene_split(Mg, D_est, f, Blist))
  ee <- ee[!vapply(ee, is.null, logical(1))]
  if (!length(ee)) return(NULL)
  EE <- do.call(rbind, ee)

  kF <- mean(EE[, "kappa"], na.rm = TRUE)
  if (!is.finite(kF)) return(NULL)
  getkb <- function(nm) if (nm %in% colnames(EE)) mean(EE[, nm], na.rm = TRUE) else NA_real_
  kb_mean <- getkb("mean")
  kb_prom <- getkb("promoter")
  kb_pnp <- getkb("prom_nonprom")

  om_mean <- 1 - kb_mean / kF
  om_prom <- 1 - kb_prom / kF
  om_pnp <- 1 - kb_pnp / kF

  random_mean <- 1 - 1 / p
  random_prom <- random_mean
  random_pnp <- if (p >= 2L) 1 - 2 / p else NA_real_

  gsrow <- gs_screen[gene == g][1]
  data.table(
    gene = g,
    p = p,
    n_prom = sum(cp$prom),
    screen_p_simes = gsrow$p_simes,
    screen_fdr = gsrow$fdr,
    screen_rank = gsrow$screen_rank,
    strict_fdr05 = is.finite(gsrow$fdr) && gsrow$fdr < FDR_CUT,
    kappa = kF,
    kappa_split_sd = sd(EE[, "kappa"], na.rm = TRUE),
    mean_abs_cor = mean_abs_cor(Mg),
    frac_same_sign_alpha_screen = frac_same_sign,
    kappaB_mean = kb_mean,
    omega2_mean = om_mean,
    random_loss_mean = random_mean,
    retained_ratio_mean_vs_random = (1 - om_mean) / (1 - random_mean),
    kappaB_prom = kb_prom,
    omega2_prom = om_prom,
    random_loss_prom = random_prom,
    retained_ratio_prom_vs_random = (1 - om_prom) / (1 - random_prom),
    kappaB_prom_nonprom = kb_pnp,
    omega2_prom_nonprom = om_pnp,
    random_loss_prom_nonprom = random_pnp,
    retained_ratio_pnp_vs_random = (1 - om_pnp) / (1 - random_pnp),
    n_splits = nrow(EE),
    standardized_means = STANDARDIZED_MEANS
  )
}

boot_median_ci <- function(x, seed) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(c(lo = NA_real_, hi = NA_real_))
  set.seed(seed)
  bm <- replicate(BOOT_REPS, median(sample(x, length(x), replace = TRUE)))
  unname(quantile(bm, c(0.025, 0.975), na.rm = TRUE))
}


pooled_loss_one <- function(audit, kb_col, random_k, seed) {
  ok <- is.finite(audit$kappa) & is.finite(audit[[kb_col]]) & is.finite(audit$p)
  aa <- audit[ok]
  if (!nrow(aa)) {
    return(c(genes = 0, pooled_loss = NA, boot_lo = NA, boot_hi = NA,
             random_loss = NA, median_kappa = NA, pct_kappa_le0 = NA,
             pct_kappa_lt01 = NA, pct_omega_outside = NA))
  }
  den <- sum(aa$kappa)
  pooled <- if (abs(den) > .Machine$double.eps)
    1 - sum(aa[[kb_col]]) / den else NA_real_
  rb <- if (abs(den) > .Machine$double.eps)
    sum(aa$kappa * (1 - random_k / aa$p)) / den else NA_real_

  set.seed(seed)
  n <- nrow(aa)
  boot <- replicate(BOOT_REPS, {
    ii <- sample.int(n, n, replace = TRUE)
    dd <- sum(aa$kappa[ii])
    if (abs(dd) <= .Machine$double.eps) NA_real_
    else 1 - sum(aa[[kb_col]][ii]) / dd
  })
  ci <- unname(quantile(boot, c(.025, .975), na.rm = TRUE))

  omega_col <- if (kb_col == "kappaB_mean") "omega2_mean" else
    if (kb_col == "kappaB_prom") "omega2_prom" else "omega2_prom_nonprom"
  om <- aa[[omega_col]]

  c(genes = n, pooled_loss = pooled, boot_lo = ci[1], boot_hi = ci[2],
    random_loss = rb, median_kappa = median(aa$kappa, na.rm = TRUE),
    pct_kappa_le0 = 100 * mean(aa$kappa <= 0, na.rm = TRUE),
    pct_kappa_lt01 = 100 * mean(aa$kappa < 0.1, na.rm = TRUE),
    pct_omega_outside = 100 * mean(om < 0 | om > 1, na.rm = TRUE))
}

pooled_summary_rule <- function(audit, direction_label, rule, seed_offset = 0L) {
  aa <- subset_rule(audit, rule)
  defs <- list(
    list(score = "Gene mean", kb = "kappaB_mean", k = 1L),
    list(score = "Promoter mean", kb = "kappaB_prom", k = 1L),
    list(score = "Promoter + non-promoter", kb = "kappaB_prom_nonprom", k = 2L)
  )
  if (is.null(aa) || !nrow(aa)) {
    return(rbindlist(lapply(defs, function(d) data.table(
      direction = direction_label, selection_rule = rule, score = d$score,
      genes = 0L, pooled_loss = NA_real_, boot_lo = NA_real_,
      boot_hi = NA_real_, kappa_weighted_random_loss = NA_real_,
      median_kappa = NA_real_, pct_kappa_le0 = NA_real_,
      pct_kappa_lt01 = NA_real_, pct_omega_outside = NA_real_
    ))))
  }
  rbindlist(lapply(seq_along(defs), function(j) {
    d <- defs[[j]]
    z <- pooled_loss_one(aa, d$kb, d$k,
                         SEED + 60000L + seed_offset + 100L * j)
    data.table(
      direction = direction_label, selection_rule = rule, score = d$score,
      genes = as.integer(z["genes"]), pooled_loss = z["pooled_loss"],
      boot_lo = z["boot_lo"], boot_hi = z["boot_hi"],
      kappa_weighted_random_loss = z["random_loss"],
      median_kappa = z["median_kappa"],
      pct_kappa_le0 = z["pct_kappa_le0"],
      pct_kappa_lt01 = z["pct_kappa_lt01"],
      pct_omega_outside = z["pct_omega_outside"]
    )
  }))
}

summarize_direction <- function(audit, direction) {
  defs <- list(
    list(score = "Gene mean", omega = "omega2_mean", rand = "random_loss_mean", ratio = "retained_ratio_mean_vs_random"),
    list(score = "Promoter mean", omega = "omega2_prom", rand = "random_loss_prom", ratio = "retained_ratio_prom_vs_random"),
    list(score = "Promoter + non-promoter", omega = "omega2_prom_nonprom", rand = "random_loss_prom_nonprom", ratio = "retained_ratio_pnp_vs_random")
  )
  rbindlist(lapply(seq_along(defs), function(j) {
    d <- defs[[j]]
    x <- audit[[d$omega]]
    b <- audit[[d$rand]]
    rr <- audit[[d$ratio]]
    ok <- is.finite(x) & is.finite(b)
    ci <- boot_median_ci(x[ok], SEED + 30000L + 100L * direction + j)
    data.table(
      direction = direction,
      score = d$score,
      genes = sum(ok),
      median_omega2 = median(x[ok], na.rm = TRUE),
      q25_omega2 = unname(quantile(x[ok], 0.25, na.rm = TRUE)),
      q75_omega2 = unname(quantile(x[ok], 0.75, na.rm = TRUE)),
      gene_boot_median_lo = ci[1],
      gene_boot_median_hi = ci[2],
      pct_gt50 = 100 * mean(x[ok] > 0.5),
      pct_gt80 = 100 * mean(x[ok] > 0.8),
      median_random_loss = median(b[ok]),
      pct_better_than_random = 100 * mean(x[ok] < b[ok]),
      pct_worse_than_random = 100 * mean(x[ok] > b[ok]),
      median_retained_ratio_vs_random = median(rr[ok], na.rm = TRUE)
    )
  }))
}

summarize_mechanism <- function(audit, direction) {
  aa <- audit[is.finite(omega2_mean) & is.finite(p) & is.finite(mean_abs_cor)]
  data.table(
    direction = direction,
    genes = nrow(aa),
    spearman_omega2_p = if (nrow(aa) > 2) cor(aa$omega2_mean, aa$p, method = "spearman") else NA_real_,
    spearman_omega2_mean_abs_cor = if (nrow(aa) > 2) cor(aa$omega2_mean, aa$mean_abs_cor, method = "spearman") else NA_real_
  )
}

summarize_sign <- function(audit, direction) {
  rbindlist(lapply(c(0.80, 0.90, 1.00), function(thr) {
    aa <- audit[is.finite(frac_same_sign_alpha_screen) & frac_same_sign_alpha_screen >= thr & is.finite(omega2_mean)]
    data.table(
      direction = direction,
      min_same_sign_fraction = thr,
      genes = nrow(aa),
      median_omega2_mean = if (nrow(aa)) median(aa$omega2_mean) else NA_real_,
      pct_mean_gt50 = if (nrow(aa)) 100 * mean(aa$omega2_mean > 0.5) else NA_real_,
      median_random_loss = if (nrow(aa)) median(aa$random_loss_mean) else NA_real_,
      pct_better_than_random = if (nrow(aa)) 100 * mean(aa$omega2_mean < aa$random_loss_mean) else NA_real_
    )
  }))
}

run_direction <- function(screen_half, estimate_half, direction_id) {
  screen_idx <- which(outer == screen_half)
  est_idx <- which(outer == estimate_half)
  label <- paste0("screen", screen_half, "_estimate", estimate_half)
  cat("\n", label, ": screening on n=", length(screen_idx),
      ", estimating on n=", length(est_idx), "\n", sep = "")

  gs_all <- screen_genes(screen_idx)
  n_strict <- sum(gs_all$fdr < FDR_CUT, na.rm = TRUE)
  cat("Strict FDR ", FDR_CUT, " discoveries: ", n_strict, "\n", sep = "")
  cat("Rank-based sensitivity: estimating the top ",
      paste(TOPK, collapse = ", "), " screened genes independently\n", sep = "")

  ## Evaluate the union of all strict discoveries and the largest requested top-K
  ## set. Smaller top-K analyses are nested subsets and therefore do not require
  ## re-estimation.
  eval_genes <- unique(c(
    gs_all[fdr < FDR_CUT, gene],
    gs_all[screen_rank <= min(KMAX, .N), gene]
  ))
  if (!length(eval_genes)) return(list(audit = NULL, screen = gs_all))

  out <- rbindlist(
    lapply(eval_genes, function(g)
      estimate_gene(g, gs_all, screen_idx, est_idx, direction_id)),
    fill = TRUE
  )
  if (!nrow(out)) return(list(audit = NULL, screen = gs_all))
  out[, direction := label]
  list(audit = out[], screen = gs_all[])
}

subset_rule <- function(audit, rule) {
  if (is.null(audit) || !nrow(audit)) return(NULL)
  if (rule == "strict_FDR05") return(audit[strict_fdr05 %in% TRUE])
  if (grepl("^top[0-9]+$", rule)) {
    kk <- as.integer(sub("^top", "", rule))
    return(audit[screen_rank <= kk])
  }
  stop("Unknown selection rule: ", rule)
}

summarize_rule <- function(audit, direction_label, direction_id, rule) {
  aa <- subset_rule(audit, rule)
  if (is.null(aa) || !nrow(aa)) {
    return(data.table(
      direction = direction_label, selection_rule = rule, score = "NONE",
      genes = 0L, median_omega2 = NA_real_, q25_omega2 = NA_real_,
      q75_omega2 = NA_real_, gene_boot_median_lo = NA_real_,
      gene_boot_median_hi = NA_real_, pct_gt50 = NA_real_, pct_gt80 = NA_real_,
      median_random_loss = NA_real_, pct_better_than_random = NA_real_,
      pct_worse_than_random = NA_real_, median_retained_ratio_vs_random = NA_real_
    ))
  }
  ## summarize_direction uses its direction argument only to construct a
  ## reproducible bootstrap seed, so pass a numeric id and relabel afterward.
  z <- summarize_direction(aa, direction_id)
  ## summarize_direction creates an integer seed-id column. Remove it before
  ## attaching the human-readable direction label to avoid data.table coercion.
  z[, direction := NULL]
  z[, direction := direction_label]
  z[, selection_rule := rule]
  setcolorder(z, c("direction", "selection_rule",
                   setdiff(names(z), c("direction", "selection_rule"))))
  z[]
}

mechanism_rule <- function(audit, direction_label, rule) {
  aa <- subset_rule(audit, rule)
  if (is.null(aa) || !nrow(aa)) {
    return(data.table(direction = direction_label, selection_rule = rule, genes = 0L,
                      spearman_omega2_p = NA_real_,
                      spearman_omega2_mean_abs_cor = NA_real_))
  }
  z <- summarize_mechanism(aa, direction_label)
  z[, selection_rule := rule]
  setcolorder(z, c("direction", "selection_rule",
                   setdiff(names(z), c("direction", "selection_rule"))))
  z[]
}

sign_rule <- function(audit, direction_label, rule) {
  aa <- subset_rule(audit, rule)
  if (is.null(aa) || !nrow(aa)) {
    return(data.table(direction = direction_label, selection_rule = rule,
                      min_same_sign_fraction = c(0.80, 0.90, 1.00),
                      genes = 0L, median_omega2_mean = NA_real_,
                      pct_mean_gt50 = NA_real_, median_random_loss = NA_real_,
                      pct_better_than_random = NA_real_))
  }
  z <- summarize_sign(aa, direction_label)
  z[, selection_rule := rule]
  setcolorder(z, c("direction", "selection_rule",
                   setdiff(names(z), c("direction", "selection_rule"))))
  z[]
}

cat("Running independent-screening sensitivity audit...\n")
r12 <- run_direction(1L, 2L, 1L)
r21 <- run_direction(2L, 1L, 2L)
a12 <- r12$audit
a21 <- r21$audit

## The strict-FDR result is retained even when one direction has zero discoveries.
## That lack of half-sample discoveries is itself a screening-stability result;
## the rank-based analyses below provide estimable independent-screen sensitivities.
screen_counts <- rbindlist(list(
  data.table(
    direction = "screen1_estimate2",
    n_screen = sum(outer == 1L),
    n_estimate = sum(outer == 2L),
    eligible_genes = nrow(r12$screen),
    strict_fdr05_genes = sum(r12$screen$fdr < FDR_CUT, na.rm = TRUE)
  ),
  data.table(
    direction = "screen2_estimate1",
    n_screen = sum(outer == 2L),
    n_estimate = sum(outer == 1L),
    eligible_genes = nrow(r21$screen),
    strict_fdr05_genes = sum(r21$screen$fdr < FDR_CUT, na.rm = TRUE)
  )
))

audit_list <- Filter(function(x) !is.null(x) && nrow(x), list(a12, a21))
if (!length(audit_list)) stop("No genes were estimable in either direction.")
audit <- rbindlist(audit_list, fill = TRUE)
setcolorder(audit, c("direction", setdiff(names(audit), "direction")))

rules <- c("strict_FDR05", paste0("top", TOPK))
summ <- rbindlist(c(
  lapply(rules, function(rr) summarize_rule(a12, "screen1_estimate2", 1L, rr)),
  lapply(rules, function(rr) summarize_rule(a21, "screen2_estimate1", 2L, rr))
), fill = TRUE)

mech <- rbindlist(c(
  lapply(rules, function(rr) mechanism_rule(a12, "screen1_estimate2", rr)),
  lapply(rules, function(rr) mechanism_rule(a21, "screen2_estimate1", rr))
), fill = TRUE)

signsum <- rbindlist(c(
  lapply(rules, function(rr) sign_rule(a12, "screen1_estimate2", rr)),
  lapply(rules, function(rr) sign_rule(a21, "screen2_estimate1", rr))
), fill = TRUE)

pooled_cross <- rbindlist(c(
  lapply(seq_along(rules), function(j)
    pooled_summary_rule(a12, "screen1_estimate2", rules[j], 1000L + 10L * j)),
  lapply(seq_along(rules), function(j)
    pooled_summary_rule(a21, "screen2_estimate1", rules[j], 2000L + 10L * j))
), fill = TRUE)

## Screening-list overlap is descriptive only; each direction remains independent
## because adequacy is always estimated on the opposite half.
overlap <- rbindlist(lapply(TOPK, function(kk) {
  g1 <- r12$screen[screen_rank <= kk, gene]
  g2 <- r21$screen[screen_rank <= kk, gene]
  data.table(
    top_k = kk,
    overlap_n = length(intersect(g1, g2)),
    union_n = length(union(g1, g2)),
    jaccard = length(intersect(g1, g2)) / length(union(g1, g2))
  )
}))


full_pooled <- NULL
full_by_p <- NULL
full_file <- file.path("tcga_lung", "genome_wide_audit.csv")
if (file.exists(full_file)) {
  full_audit <- fread(full_file)
  zz <- list(
    pooled_loss_one(full_audit, "kappaB_mean", 1L, SEED + 71001L),
    pooled_loss_one(full_audit, "kappaB_prom", 1L, SEED + 71002L),
    pooled_loss_one(full_audit, "kappaB_prom_nonprom", 2L, SEED + 71003L)
  )
  labs_score <- c("Gene mean", "Promoter mean", "Promoter + non-promoter")
  full_pooled <- rbindlist(lapply(seq_along(zz), function(j) {
    z <- zz[[j]]
    data.table(score = labs_score[j], genes = as.integer(z["genes"]),
               pooled_loss = z["pooled_loss"], boot_lo = z["boot_lo"],
               boot_hi = z["boot_hi"],
               kappa_weighted_random_loss = z["random_loss"],
               median_kappa = z["median_kappa"],
               pct_kappa_le0 = z["pct_kappa_le0"],
               pct_kappa_lt01 = z["pct_kappa_lt01"],
               pct_omega_outside = z["pct_omega_outside"])
  }))

  brks <- c(4, 10, 25, 60, 150)
  labs <- c("5-10", "11-25", "26-60", "61-150")
  full_audit[, p_bin := cut(p, breaks = brks, labels = labs, include.lowest = TRUE)]
  full_by_p <- rbindlist(lapply(labs, function(lb) {
    dd <- full_audit[p_bin == lb]
    z1 <- pooled_loss_one(dd, "kappaB_mean", 1L,
                          SEED + 72000L + match(lb, labs))
    z2 <- pooled_loss_one(dd, "kappaB_prom_nonprom", 2L,
                          SEED + 73000L + match(lb, labs))
    data.table(
      p_bin = lb, genes = nrow(dd),
      gene_mean_pooled_loss = z1["pooled_loss"],
      gene_mean_random_loss = z1["random_loss"],
      two_component_genes = as.integer(z2["genes"]),
      two_component_pooled_loss = z2["pooled_loss"],
      two_component_random_loss = z2["random_loss"]
    )
  }))
}

cat("\nStrict half-sample screening counts:\n")
print(screen_counts)
cat("\nCross-screened adequacy summaries:\n")
print(summ, digits = 4)
cat("\nMechanism summaries:\n")
print(mech, digits = 4)
cat("\nSign-consistency summaries (signs estimated in screening half):\n")
print(signsum, digits = 4)
cat("\nPooled ratio-of-sums summaries:\n")
print(pooled_cross, digits = 4)
cat("\nTop-K screening overlap:\n")
print(overlap, digits = 4)

outdir <- "tcga_lung"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
fwrite(audit, file.path(outdir, "cross_screened_audit.csv"))
fwrite(screen_counts, file.path(outdir, "cross_screened_screening_counts.csv"))
fwrite(summ, file.path(outdir, "cross_screened_audit_summary.csv"))
fwrite(mech, file.path(outdir, "cross_screened_audit_mechanism.csv"))
fwrite(signsum, file.path(outdir, "cross_screened_audit_sign_summary.csv"))
fwrite(pooled_cross, file.path(outdir, "cross_screened_pooled_summary.csv"))
fwrite(overlap, file.path(outdir, "cross_screened_screening_overlap.csv"))
if (!is.null(full_pooled))
  fwrite(full_pooled, file.path(outdir, "genome_wide_pooled_summary.csv"))
if (!is.null(full_by_p))
  fwrite(full_by_p, file.path(outdir, "genome_wide_pooled_by_p.csv"))
capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo_cross_screened_audit.txt"))

cat("\nSaved:\n",
    file.path(outdir, "cross_screened_audit.csv"), "\n",
    file.path(outdir, "cross_screened_screening_counts.csv"), "\n",
    file.path(outdir, "cross_screened_audit_summary.csv"), "\n",
    file.path(outdir, "cross_screened_audit_mechanism.csv"), "\n",
    file.path(outdir, "cross_screened_audit_sign_summary.csv"), "\n",
    file.path(outdir, "cross_screened_pooled_summary.csv"), "\n",
    file.path(outdir, "cross_screened_screening_overlap.csv"), "\n",
    if (!is.null(full_pooled)) file.path(outdir, "genome_wide_pooled_summary.csv") else "", "\n",
    if (!is.null(full_by_p)) file.path(outdir, "genome_wide_pooled_by_p.csv") else "", "\n",
    file.path(outdir, "sessionInfo_cross_screened_audit.txt"), "\n", sep = "")
