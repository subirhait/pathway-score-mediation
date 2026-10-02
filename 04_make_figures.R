#!/usr/bin/env Rscript
## R-only figure regeneration for the pathway-score manuscript.
## Requires outputs from 02_simulation_and_application.R and 06_cross_screened_audit.R.

ROOT <- Sys.getenv("PATHWAY_SCORE_ROOT", unset = getwd())
setwd(ROOT)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

## ---------------- Figure 1: adequacy profiles ----------------
simfile <- "tcga_lung/simulation_results.rds"
if (!file.exists(simfile)) stop("Missing ", simfile, ". Run 02_simulation_and_application.R first.")
sim <- readRDS(simfile)
chk <- as.data.frame(sim$chk)
spec <- as.data.frame(sim$spec)

pc_names <- c("PC1","PC2","PC5","PC10","PC20","PC50","PC100")
keep_names <- c("sign","means","alpha_wt","disc_oracle",pc_names)
prof <- chk[match(keep_names, chk$rep), ]

pdf("figures/figure1_profiles.pdf", width = 8.2, height = 4.2, useDingbats = FALSE)
par(mfrow = c(1,2), mar = c(4.2,4.4,2.3,0.8), mgp = c(2.4,0.7,0))
pc <- prof[prof$rep %in% pc_names, ]
plot(pc$m, pc$gamma, type = "b", pch = 16,
     xlab = expression("Systematic distortion "*m[B]),
     ylab = expression("Unpredictable distortion "*gamma[B]),
     main = "(a) Population profiles", ylim = c(0,0.60))
text(pc$m, pc$gamma, labels = sub("PC","",pc$rep), pos = 4, cex = 0.72)
for (nm in c("sign","means","alpha_wt","disc_oracle")) {
  z <- prof[prof$rep == nm, ]
  points(z$m, z$gamma, pch = switch(nm, sign=15, means=17, alpha_wt=18, disc_oracle=8), cex = 1.15)
}
legend("bottomright", legend = c("PCs","Sign score","Hyper/hypo means",expression(alpha*"-weighted"),expression(Sigma^{-1}*alpha*" (exact)")),
       pch = c(16,15,17,18,8), bty = "n", cex = 0.72)
abline(v = 0, lty = 3)

sel <- c("sign","alpha_wt","means","PC1","PC10","PC100")
plot(NULL, xlim = range(spec$g), ylim = c(0.15,0.62),
     xlab = "Spectrum exponent g (larger = stronger correlation)",
     ylab = expression(gamma[B]), main = "(b) Dependence on correlation strength")
pchs <- c(15,18,17,16,1,0)
for (j in seq_along(sel)) {
  d <- spec[spec$rep == sel[j], ]
  lines(d$g, d$gamma, type = "b", pch = pchs[j])
}
legend("bottomright", legend = c("Sign score",expression(alpha*"-weighted"),"Hyper/hypo means","PC1","PC10","PC100"),
       pch = pchs, lty = 1, bty = "n", cex = 0.72, ncol = 2)
dev.off()

## ---------------- Figure 2: application forest plot ----------------
appfile <- "tcga_lung/application_results.rds"
if (!file.exists(appfile)) stop("Missing ", appfile, ". Run 02_simulation_and_application.R first.")
app <- readRDS(appfile)
rob <- as.data.frame(app$rob)
nqo1 <- as.data.frame(app$nqo1)
get_rows <- function(g) {
  if (g == "NQO1") d <- nqo1[nqo1$gene == g & nqo1$spec == "main", ]
  else d <- rob[rob$gene == g & rob$spec == "main", ]
  wanted <- c("IE_M","gene_mean","promoter_mean","nonpromoter_mean","PC1_refit","PC3_refit","sign_refit")
  d <- d[d$rep %in% wanted, ]
  lab <- c(IE_M="Full CpG set", gene_mean="Gene mean", promoter_mean="Promoter mean",
           nonpromoter_mean="Non-promoter mean", PC1_refit="PC1 (refit)",
           PC3_refit="PC3 (refit)", sign_refit="Sign score (refit)")
  d$label <- unname(lab[d$rep])
  d
}

genes <- c("AHRR","UCHL1","NQO1")
pdf("figures/figure2_application.pdf", width = 9.0, height = 4.1, useDingbats = FALSE)
par(mfrow = c(1,3), mar = c(4.2,7.0,2.2,1.0), mgp = c(2.3,0.7,0))
for (g in genes) {
  d <- get_rows(g)
  y <- rev(seq_len(nrow(d)))
  xr <- range(c(d$lo,d$hi), na.rm = TRUE)
  pad <- 0.08 * diff(xr); if (!is.finite(pad) || pad == 0) pad <- 0.1
  plot(d$est, y, xlim = xr + c(-pad,pad), ylim = c(0.5,nrow(d)+0.5), pch = 19,
       yaxt = "n", ylab = "", xlab = "Indirect effect (log2 expression)", main = g)
  segments(d$lo, y, d$hi, y, lwd = 1.5)
  axis(2, at = y, labels = d$label, las = 1, cex.axis = 0.72)
  abline(v = 0, lty = 2)
}
dev.off()

## ---------------- Figure 3: pooled audit benchmark ----------------
poolfile <- "tcga_lung/genome_wide_pooled_by_p.csv"
if (!file.exists(poolfile)) stop("Missing ", poolfile, ". Run 06_cross_screened_audit.R first.")
pp <- read.csv(poolfile, check.names = FALSE)
if ("p_bin" %in% names(pp)) {
  xlab <- pp$p_bin
  gm <- pp$gene_mean_pooled_loss
  gr <- pp$gene_mean_random_loss
  pm <- pp$two_component_pooled_loss
  pr <- pp$two_component_random_loss
} else {
  xlab <- pp[["CpGs per gene"]]
  gm <- pp[["Gene mean pooled loss"]]
  gr <- pp[["Gene mean random benchmark"]]
  pm <- pp[["Two-component pooled loss"]]
  pr <- pp[["Two-component random benchmark"]]
}
x <- seq_along(xlab)
pdf("figures/figure3_pooled_benchmark.pdf", width = 8.2, height = 5.0, useDingbats = FALSE)
par(mar = c(4.3,4.5,1.2,0.8), mgp = c(2.5,0.7,0))
plot(x, gm, type = "b", pch = 16, ylim = c(0,1.02), xaxt = "n",
     xlab = "CpGs per gene", ylab = "Exposure-information loss")
axis(1, at = x, labels = xlab)
lines(x, gr, type = "b", pch = 1, lty = 2)
lines(x, pm, type = "b", pch = 15)
lines(x, pr, type = "b", pch = 0, lty = 2)
legend("bottomright", legend = c("Gene mean: pooled loss","Gene mean: random benchmark",
                                 "Promoter + non-promoter: pooled loss","Two-component: random benchmark"),
       pch = c(16,1,15,0), lty = c(1,2,1,2), bty = "n", cex = 0.82)
dev.off()

cat("Wrote figures/figure1_profiles.pdf\n")
cat("Wrote figures/figure2_application.pdf\n")
cat("Wrote figures/figure3_pooled_benchmark.pdf\n")
