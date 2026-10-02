## =====================================================================
## Reproducibility script: simulations and application
## "When Does a Pathway Score Preserve a Mediated Effect?"
##
## Requires tcga_lung/luad_stage1.rds, created by 01_data_prep.R
## (the earlier tcga_luad_aggregation.R, stages 0-3).
## Every block below is code that was run for the manuscript.
## =====================================================================
## Run from the repository/project root. Optionally set PATHWAY_SCORE_ROOT.
ROOT <- Sys.getenv("PATHWAY_SCORE_ROOT", unset = getwd())
setwd(ROOT)
library(data.table); library(survival); library(pseudo); library(glmnet)
library(dplyr); library(UCSCXenaTools)

s2 <- readRDS("tcga_lung/luad_stage1.rds")
Mval <- s2$Mval; ph <- s2$ph; ann <- s2$ann; tt <- s2$alpha; rm(s2); gc()
ph[, stage4 := droplevels(stage4)]

## ---- 1. Screened mediators and the cautionary survival analysis -------
S  <- rownames(tt)[tt$adj.P.Val < 0.05]                  # 4,625 CpGs
MS <- t(Mval[S, ])
Xc <- model.matrix(~ age + sex + stage4 + purity_cpe, data = ph)[, -1]
A  <- ph$A
tau_rmst <- 1095
ph[, Y := pseudomean(OS.time, OS, tmax = tau_rmst) / 365.25]
Y <- ph$Y
print(summary(lm(Y ~ A + Xc))$coefficients["A", ])      # total effect

set.seed(2026)
H <- cbind(Xc, MS); pf <- c(rep(0, ncol(Xc)), rep(1, ncol(MS)))
fold <- sample(rep(1:5, length.out = nrow(H))); rY <- rA <- numeric(nrow(H))
for (k in 1:5) {
  tr <- fold != k; te <- fold == k
  fy <- cv.glmnet(H[tr, ], Y[tr], penalty.factor = pf)
  fa <- cv.glmnet(H[tr, ], A[tr], penalty.factor = pf)
  rY[te] <- Y[te] - predict(fy, H[te, ], s = "lambda.min")
  rA[te] <- A[te] - predict(fa, H[te, ], s = "lambda.min")
}
DE_M <- sum(rA * rY) / sum(rA^2)
se_M <- sqrt(sum(rA^2 * (rY - DE_M * rA)^2)) / sum(rA^2)
cat("DE_M =", DE_M, "SE =", se_M, "\n"); print(summary(A - rA))

## ---- 2. Population: real alpha, shrunk real Sigma --------------------
MSs   <- scale(MS)
a_hat <- coef(lm(MSs ~ A + Xc))["A", ]
p     <- ncol(MSs); up_i <- a_hat > 0; dn_i <- a_hat < 0
ord   <- order(-abs(a_hat))
s <- 50; size <- 0.02

R  <- resid(lm(MSs ~ A + Xc))
n  <- nrow(R)
sv <- svd(R / sqrt(n)); V <- sv$v; d2 <- sv$d^2
mu <- sum(d2) / p
delta <- 0.2
ev <- (1 - delta) * d2 + delta * mu
Sig_x  <- function(X) { X <- as.matrix(X); VX <- crossprod(V, X)
  V %*% (VX * ((1 - delta) * d2)) + delta * mu * X }
Sinv_x <- function(X) { X <- as.matrix(X); VX <- crossprod(V, X)
  V %*% (VX / ev) + (X - V %*% VX) / (delta * mu) }

make_beta <- function(scen, s = 50, size = 0.02) {
  b <- numeric(p)
  if (scen == "aligned")    { j <- ord[1:s];     b[j] <- size * sign(a_hat[j]) }
  if (scen == "cancel")     { j <- ord[1:(2*s)]; b[j] <- size * sign(a_hat[j]) * rep(c(1, -1), s) }
  if (scen == "hyper_only") { j <- which(up_i)[1:s]; b[j] <- size }
  if (scen == "weak_alpha") { j <- tail(ord, s); b[j] <- size * 5 }
  b
}

pca <- prcomp(MSs, rank. = 100)
Bs <- list(sign        = ifelse(up_i, 1 / sum(up_i), -1 / sum(dn_i)),
           means       = cbind(up_i / sum(up_i), dn_i / sum(dn_i)),
           alpha_wt    = a_hat,
           disc_oracle = drop(Sinv_x(a_hat)))
for (k in c(1, 2, 5, 10, 20, 50, 100))
  Bs[[paste0("PC", k)]] <- pca$rotation[, 1:k, drop = FALSE]

prep <- function(Blist, Sx) lapply(Blist, function(B) { B <- as.matrix(B); SB <- Sx(B)
  list(aB = crossprod(a_hat, B), K = solve(crossprod(B, SB)), SB = SB) })
pre <- prep(Bs, Sig_x)
ie_pop <- function(b, pr = pre) sapply(pr, function(P)
  drop(P$aB %*% P$K %*% crossprod(P$SB, b)))

## constructed scenarios (Supplementary Table)
truth_pop <- rbindlist(lapply(c("aligned", "cancel", "hyper_only", "weak_alpha"),
  function(sc) { b <- make_beta(sc)
    data.table(scenario = sc, rep = names(Bs), IE_M = sum(a_hat * b), IE_Z = ie_pop(b)) }))
print(dcast(truth_pop, rep ~ scenario, value.var = "IE_Z"), digits = 3)

## random outcome models (Table 1)
set.seed(2)
simp <- rbindlist(lapply(1:1000, function(r) {
  b <- numeric(p); j <- sample(p, s); b[j] <- rnorm(s, 0, size)
  data.table(r = r, rep = names(Bs), IE_M = sum(a_hat * b), IE_Z = ie_pop(b)) }))
simp[, D_B := IE_M - IE_Z]
ref <- abs(simp[rep == "sign", IE_M])
cut_hi <- quantile(ref, 0.5); cut_lo <- quantile(ref, 0.2)
summp <- simp[, .(
  med_rel_gap = median(abs(D_B / IE_M)[abs(IE_M) > cut_hi]),
  pct_gt50    = 100 * mean((abs(D_B / IE_M) > 0.5)[abs(IE_M) > cut_hi]),
  flip_pct    = 100 * mean((sign(IE_Z) != sign(IE_M))[abs(IE_M) > cut_hi]),
  false_med   = 100 * mean((abs(IE_Z) > cut_hi)[abs(IE_M) < cut_lo])), by = rep]
print(summp, digits = 3)

## ---- 3. Cauchy profile: theory vs simulation, spectrum shape ---------
cauchy_pars <- function(pr, a = a_hat) {
  rbindlist(lapply(names(pr), function(nm) { P <- pr[[nm]]
    v  <- drop(a - P$SB %*% P$K %*% t(P$aB))
    aa <- sum(a^2); va <- sum(v * a)
    data.table(rep = nm, m = va / aa,
               gamma = sqrt(max(sum(v^2) * aa - va^2, 0)) / aa) }))
}
th <- cauchy_pars(pre)

set.seed(4)
Bd <- matrix(rnorm(p * 2000), p, 2000)
dense <- rbindlist(lapply(names(pre), function(nm) { P <- pre[[nm]]
  IEM <- drop(crossprod(a_hat, Bd))
  IEZ <- drop(P$aB %*% P$K %*% crossprod(P$SB, Bd))
  Rr <- (IEM - IEZ) / IEM
  data.table(rep = nm, sim_median = median(Rr), sim_halfIQR = IQR(Rr) / 2) }))
sparse <- simp[, .(sp_median = median(D_B / IE_M),
                   sp_halfIQR = IQR(D_B / IE_M) / 2), by = rep]
chk <- Reduce(function(x, y) merge(x, y, by = "rep"), list(th, dense, sparse))
print(chk[order(-gamma)], digits = 3)

make_pop <- function(dd, delta) {
  ev <- (1 - delta) * dd + delta * mu
  list(Sx  = function(X) { X <- as.matrix(X); VX <- crossprod(V, X)
                           V %*% (VX * ((1 - delta) * dd)) + delta * mu * X },
       Six = function(X) { X <- as.matrix(X); VX <- crossprod(V, X)
                           V %*% (VX / ev) + (X - V %*% VX) / (delta * mu) })
}
spec <- rbindlist(lapply(c(0.5, 0.75, 1, 1.25, 1.5), function(g) {
  d2g <- d2^g; d2g <- d2g * sum(d2) / sum(d2g)
  pp <- make_pop(d2g, 0.2)
  B0 <- Bs; B0$disc_oracle <- drop(pp$Six(a_hat))
  cauchy_pars(prep(B0, pp$Sx))[, g := g]
}))
print(dcast(spec, rep ~ g, value.var = "gamma"), digits = 3)

## sensitivity grid: shrinkage x sparsity
run_pop <- function(delta, s, n_draw = 1000, seed = 2) {
  pp <- make_pop(d2, delta)
  B0 <- Bs; B0$disc_oracle <- drop(pp$Six(a_hat))
  pr <- prep(B0, pp$Sx)
  set.seed(seed)
  out <- rbindlist(lapply(1:n_draw, function(r) {
    b <- numeric(p); j <- sample(p, s); b[j] <- rnorm(s)
    data.table(r = r, rep = names(pr), IE_M = sum(a_hat * b), IE_Z = ie_pop(b, pr)) }))
  refm <- abs(out[rep == "sign", IE_M])
  hi <- quantile(refm, 0.5); lo <- quantile(refm, 0.2)
  out[, .(pct_gt50  = 100 * mean((abs((IE_M - IE_Z) / IE_M) > 0.5)[abs(IE_M) > hi]),
          false_med = 100 * mean((abs(IE_Z) > hi)[abs(IE_M) < lo])),
      by = rep][, `:=`(delta = delta, s = s)]
}
grid <- CJ(delta = c(0.1, 0.2, 0.5), s = c(10, 50, 200))
sens <- rbindlist(Map(run_pop, grid$delta, grid$s))
print(dcast(sens[rep %in% c("sign", "alpha_wt", "PC1", "PC10", "PC100")],
            delta + s ~ rep, value.var = "pct_gt50"), digits = 3)

## ---- 4. Estimation study ---------------------------------------------
Dm <- cbind(1, A, Xc)
Cf <- coef(lm(MSs ~ A + Xc))
r_ <- length(d2); sdev <- sqrt((1 - delta) * d2)
gen_M <- function()
  Dm %*% Cf + matrix(rnorm(n * r_), n, r_) %*% (sdev * t(V)) +
    sqrt(delta * mu) * matrix(rnorm(n * p), n, p)
fitfold <- function(Mg, idx) {
  Q <- qr(Dm[idx, ])
  list(a = qr.coef(Q, Mg[idx, ])[2, ], R = qr.resid(Q, Mg[idx, ]),
       df = length(idx) - ncol(Dm))
}
SX <- function(f, X) crossprod(f$R, f$R %*% X) / f$df
split_k <- function(Kf = 4) { f <- integer(n)
  for (g0 in 0:1) { w <- which(A == g0); f[w] <- sample(rep(1:Kf, length.out = length(w))) }
  lapply(1:Kf, function(k) which(f == k)) }
perms <- as.matrix(expand.grid(1:4, 1:4, 1:4, 1:4))
perms <- perms[apply(perms, 1, function(x) length(unique(x)) == 4), ]

est_pars4 <- function(B, fl, K) {
  SB <- lapply(fl, function(f) SX(f, B))
  BA <- lapply(fl, function(f) crossprod(B, f$a))
  q  <- function(x, y, m) drop(crossprod(BA[[x]], K %*% crossprod(SB[[m]], fl[[y]]$a)))
  acc <- c(aa = 0, va = 0, num = 0, aa4 = 0)
  for (r0 in seq_len(nrow(perms))) {
    i <- perms[r0, 1]; j <- perms[r0, 2]; k <- perms[r0, 3]; l <- perms[r0, 4]
    aij <- sum(fl[[i]]$a * fl[[j]]$a); akl <- sum(fl[[k]]$a * fl[[l]]$a)
    va1 <- aij - q(i, j, k); va2 <- akl - q(k, l, i)
    vv  <- akl - 2 * q(k, l, i) +
      drop(crossprod(BA[[k]], K %*% crossprod(SB[[i]], SB[[j]]) %*% K %*% BA[[l]]))
    acc <- acc + c(aij, va1, vv * aij - va1 * va2, aij * akl)
  }
  c(m = unname(acc["va"] / acc["aa"]), gamma2 = unname(acc["num"] / acc["aa4"]))
}
pop_eval <- function(b, nm) {
  pr <- setNames(prep(list(b), Sig_x), nm); cauchy_pars(pr)
}
keepB <- lapply(Bs[c("sign", "means", "alpha_wt", "PC1", "PC10", "PC20",
                     "PC50", "PC100", "disc_oracle")], as.matrix)

set.seed(6)
e2 <- rbindlist(lapply(1:100, function(it) {
  Mg <- gen_M(); ff <- fitfold(Mg, 1:n)
  fl <- lapply(split_k(4), function(ix) fitfold(Mg, ix))
  rbindlist(lapply(names(keepB), function(nm) { B <- keepB[[nm]]
    K <- solve(crossprod(B, SX(ff, B)))
    e <- est_pars4(B, fl, K)
    data.table(it = it, rep = nm, m4 = e[["m"]], g2_4 = e[["gamma2"]]) }))
}))
print(merge(e2[, .(m4 = mean(m4), g4 = sqrt(max(mean(g2_4), 0)),
                   pct_g2_neg = 100 * mean(g2_4 < 0)), by = rep],
            th[, .(rep, m_true = m, g_true = gamma)], by = "rep"), digits = 3)

pf_M <- c(rep(0, ncol(Xc)), rep(1, p))
set.seed(9)
p4 <- rbindlist(lapply(1:30, function(it) {
  Mg <- gen_M(); ff <- fitfold(Mg, 1:n)
  up_h <- ff$a > 0; dn_h <- ff$a < 0
  W <- list(sign_hat     = ifelse(up_h, 1 / sum(up_h), -1 / sum(dn_h)),
            means_hat    = cbind(up_h / sum(up_h), dn_h / sum(dn_h)),
            alpha_wt_hat = ff$a)
  pcs <- prcomp(Mg, rank. = 50)$rotation
  for (k in c(1, 10, 20, 50)) W[[paste0("PC", k, "_hat")]] <- pcs[, 1:k, drop = FALSE]
  Hs <- cbind(Xc, Mg)
  for (al in c(0.5, 1)) {
    fit <- cv.glmnet(Hs, A, alpha = al, penalty.factor = pf_M, nfolds = 5)
    w <- as.numeric(coef(fit, s = "lambda.min"))[-(1:(ncol(Xc) + 1))]
    if (any(w != 0)) W[[paste0("AonM_alpha", al)]] <- w
  }
  rbindlist(lapply(names(W), function(nm) pop_eval(W[[nm]], nm)))[, it := it]
}))
print(p4[, .(m = mean(m), sd_m = sd(m), gamma = mean(gamma)), by = rep][order(gamma)],
      digits = 3)

saveRDS(list(truth_pop = truth_pop, summp = summp, chk = chk, spec = spec,
             sens = sens, e2 = e2, p4 = p4), "tcga_lung/simulation_results.rds")

## ---- 5. Application: smoking -> gene CpGs -> gene expression ---------
pth <- "tcga_lung/TCGA.LUNG.sampleMap/"
fx  <- list.files(pth, pattern = "^HiSeqV2(\\.gz)?$", full.names = TRUE)
if (!length(fx)) {
  XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
    XenaFilter(filterCohorts  = "TCGA Lung Cancer \\(LUNG\\)") %>%
    XenaFilter(filterDatasets = "sampleMap/HiSeqV2$") %>%
    XenaQuery() %>% XenaDownload(destdir = "tcga_lung")
  fx <- list.files(pth, pattern = "^HiSeqV2(\\.gz)?$", full.names = TRUE)
}
ex <- fread(fx[1]); setnames(ex, 1, "gene")
common <- intersect(ph$sample, names(ex))
phx <- ph[match(common, ph$sample)]; phx[, stage4 := droplevels(stage4)]
Ax  <- phx$A
getY <- function(g) as.numeric(unlist(ex[gene == g, ..common]))
Xx  <- model.matrix(~ age + sex + stage4 + purity_cpe, data = phx)[, -1]

cand <- c("CYP1B1", "CYP1A1", "AHRR", "ALDH3A1", "NQO1", "AKR1B10", "SLC7A11", "UCHL1")
te_tab <- rbindlist(lapply(intersect(cand, ex$gene), function(g) {
  cf <- summary(lm(getY(g) ~ Ax + Xx))$coefficients["Ax", ]
  data.table(gene = g, TE = cf[1], se = cf[2], t = cf[3])
}))[order(-abs(t))]
print(te_tab, digits = 3)

gene_list <- strsplit(ann$UCSC_RefGene_Name, ";")
grp_list  <- strsplit(ann$UCSC_RefGene_Group, ";")
gene_cpgs <- function(g0) {
  hit <- vapply(gene_list, function(z) g0 %in% z, logical(1))
  grp <- mapply(function(nm, gr) paste(unique(gr[nm == g0]), collapse = ";"),
                gene_list[hit], grp_list[hit])
  list(Mm = t(Mval[rownames(ann)[hit], match(common, colnames(Mval)), drop = FALSE]),
       prom = grepl("TSS1500|TSS200|5'UTR|1stExon", grp))
}

## plug-in adequacy profile for fixed aggregations (detailed application table: columns m, gamma)
profile_gene <- function(g0) {
  gc_ <- gene_cpgs(g0); Mm <- gc_$Mm; prom <- gc_$prom; pg <- ncol(Mm)
  D0 <- cbind(1, Ax, Xx); sds <- apply(Mm, 2, sd)
  up <- qr.coef(qr(D0), Mm)[2, ] > 0; pc <- prcomp(Mm, scale. = TRUE)
  W <- list(gene_mean = matrix(1 / pg, pg, 1))
  if (any(prom))  W$promoter_mean <- matrix(prom / sum(prom), pg, 1)
  if (any(!prom)) W$nonpromoter_mean <- matrix((!prom) / sum(!prom), pg, 1)
  if (any(prom) && any(!prom)) W$prom_plus_nonprom <- cbind(prom / sum(prom), (!prom) / sum(!prom))
  if (any(up) && any(!up)) W$sign_score <- matrix(ifelse(up, 1 / sum(up), -1 / sum(!up)), pg, 1)
  Bl <- lapply(W, function(w) w * sds)
  for (k in c(1, 3)) if (k < pg) Bl[[paste0("PC", k)]] <- pc$rotation[, 1:k, drop = FALSE]
  Ms <- scale(Mm); QR <- qr(D0)
  a_s <- qr.coef(QR, Ms)[2, ]
  Sg  <- crossprod(qr.resid(QR, Ms)) / (nrow(Ms) - ncol(D0))
  rbindlist(lapply(names(Bl), function(nm) { B <- Bl[[nm]]; SB <- Sg %*% B
    K <- solve(crossprod(B, SB))
    v <- a_s - drop(SB %*% K %*% crossprod(B, a_s))
    aa <- sum(a_s^2); va <- sum(v * a_s)
    data.table(gene = g0, rep = nm, m = va / aa,
               gamma = sqrt(max(sum(v^2) * aa - va^2, 0)) / aa) }))
}
prof_app <- rbindlist(lapply(te_tab$gene[1:3], profile_gene))
print(prof_app, digits = 3)

## bootstrap with data-derived weights refit in each resample (Tables 3, S2)
app_rob <- function(g0, fml, label, nboot = 2000, seed = 12) {
  gc_ <- gene_cpgs(g0); Mm <- gc_$Mm; prom <- gc_$prom
  y <- getY(g0); pg <- ncol(Mm)
  D0 <- cbind(1, Ax, model.matrix(fml, data = phx)[, -1, drop = FALSE])
  fixW <- list(gene_mean = rep(1 / pg, pg))
  if (any(prom))  fixW$promoter_mean <- prom / sum(prom)
  if (any(!prom)) fixW$nonpromoter_mean <- (!prom) / sum(!prom)
  a0 <- lm.fit(D0, Mm)$coefficients[2, ]
  use_sign <- any(a0 > 0) && any(a0 < 0)
  est <- function(id) {
    D <- D0[id, , drop = FALSE]; yy <- y[id]; Mi <- Mm[id, , drop = FALSE]
    bA <- function(DD) unname(lm.fit(DD, yy)$coefficients[2])
    TE <- bA(D); ie <- function(z) TE - bA(cbind(D, z))
    out <- c(TE = TE, IE_M = ie(Mi), sapply(fixW, function(w) ie(Mi %*% w)))
    if (use_sign) {
      up <- lm.fit(D, Mi)$coefficients[2, ] > 0
      out["sign_refit"] <- if (any(up) && any(!up))
        ie(Mi %*% ifelse(up, 1 / sum(up), -1 / sum(!up))) else NA
    }
    pcx <- prcomp(Mi, scale. = TRUE)$x
    out["PC1_refit"] <- ie(pcx[, 1, drop = FALSE])
    if (pg > 3) out["PC3_refit"] <- ie(pcx[, 1:3])
    out
  }
  boot_id <- function() unlist(lapply(split(seq_along(y), Ax),
                                      function(i) sample(i, replace = TRUE)))
  set.seed(seed)
  pt <- est(seq_along(y))
  bt <- t(replicate(nboot, est(boot_id())))
  zn <- setdiff(names(pt), c("TE", "IE_M")); nm <- c("TE", "IE_M", zn)
  q  <- function(x) quantile(x, c(.025, .975), na.rm = TRUE)
  gq <- apply(bt[, "IE_M"] - bt[, zn, drop = FALSE], 2, q)
  data.table(gene = g0, spec = label, rep = nm, est = pt[nm],
             lo = apply(bt[, nm], 2, q)[1, ], hi = apply(bt[, nm], 2, q)[2, ],
             gap = c(NA, NA, pt["IE_M"] - pt[zn]),
             gap_lo = c(NA, NA, gq[1, ]), gap_hi = c(NA, NA, gq[2, ]))
}
specs <- list(main        = ~ age + sex + stage4 + purity_cpe,
              no_stage    = ~ age + sex + purity_cpe,
              no_purity   = ~ age + sex + stage4,
              demographic = ~ age + sex)
rob <- rbindlist(lapply(c("AHRR", "UCHL1"), function(g)
  rbindlist(lapply(names(specs), function(s0) app_rob(g, specs[[s0]], s0)))))
## NQO1 provenance note: The detailed application table retains the earlier fixed-score
## run used when that table was assembled. The following call refits learned
## PC/sign weights inside each resample (seed 11) and is a robustness analysis;
## its learned-score intervals can therefore differ slightly from the detailed application table.
nqo1 <- app_rob("NQO1", specs$main, "main", nboot = 1000, seed = 11)
print(rob, digits = 3); print(nqo1, digits = 3)

saveRDS(list(te_tab = te_tab, prof_app = prof_app, rob = rob, nqo1 = nqo1),
        "tcga_lung/application_results.rds")
sessionInfo()
