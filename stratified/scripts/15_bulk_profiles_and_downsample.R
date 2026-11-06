#!/usr/bin/env Rscript
# 15_bulk_profiles_and_downsample.R — two per-gene summaries of the isolated-hepatocyte bulk cohort, for the
# concordance readouts against hepatocyte pseudobulk (Spearman, significant -> not significant flips):
#   bulk_group_logcpm.csv   mean log2 CPM (TMM) per gene over the normal donors and over the in-range ESLD preps
#   bulk_normal_halves_logcpm.csv  the same over two random halves of the normal donors (same-platform ceiling)
#   bulk_downsample_t.csv.gz  the cohort's reference model (~ group + logit(MT) + age + sex + log10(assigned)) refitted on
#                           N_DRAW random 6-ESLD vs 4-donor subsets, the size of the PSC/PBC vs healthy pseudobulk
#                           design: how many full-cohort calls survive at that n with the platform held fixed
# Only per-gene statistics are written (no per-patient values).
# usage: Rscript 15_bulk_profiles_and_downsample.R <dir with meta.csv, counts_noMT_noNUMT.csv.gz> <out>
#        [<etiology_group> <n_disease> <n_healthy> <draws file name>]
# With the optional arguments only the draws are written: disease preps restricted to one etiology group, drawn at
# the size of another pseudobulk design (e.g. MASH 9 vs 4 for the Gribben et al. cirrhotic-MASH comparison).
suppressPackageStartupMessages({ library(edgeR); library(limma) })
a <- commandArgs(trailingOnly = TRUE)
inp <- a[1]; out <- a[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
STRATUM <- if (length(a) >= 3) a[3] else NA
N_DRAW <- 50
N_DIS <- if (length(a) >= 4) as.integer(a[4]) else 6
N_HEALTHY <- if (length(a) >= 5) as.integer(a[5]) else 4
DRAWS_FILE <- if (length(a) >= 6) a[6] else "bulk_downsample_t.csv.gz"
cnt  <- as.matrix(read.csv(gzfile(file.path(inp, "counts_noMT_noNUMT.csv.gz")), row.names = 1, check.names = FALSE))
meta <- read.csv(file.path(inp, "meta.csv"), row.names = 1, stringsAsFactors = FALSE)[colnames(cnt), ]
normal <- rownames(meta)[meta$role == "normal_donor"]
esld   <- rownames(meta)[meta$role == "esld" & meta$esld_in_mt_range == "True"]
if (!is.na(STRATUM)) esld <- esld[meta[esld, "etiology_group"] %in% STRATUM]

if (is.na(STRATUM)) {
  y <- calcNormFactors(DGEList(cnt[, c(normal, esld)]), method = "TMM")
  lc <- cpm(y, log = TRUE, prior.count = 1)
  write.csv(data.frame(normal_mean_log2cpm = rowMeans(lc[, normal]), esld_mean_log2cpm = rowMeans(lc[, esld])),
            file.path(out, "bulk_group_logcpm.csv"))
}

fit_t <- function(ids) {
  grp <- factor(ifelse(ids %in% normal, "normal", "disease"), levels = c("normal", "disease"))
  m <- meta[ids, ]
  mt <- pmin(pmax(m$mt_pct / 100, 1e-4), 1 - 1e-4)
  cv <- data.frame(group = grp, logit_mt = log(mt / (1 - mt)), age = m$age, sex = factor(m$sex),
                   log10_assigned = log10(m$assigned))
  des <- model.matrix(~ group + logit_mt + age + sex + log10_assigned, cv)
  if (qr(des)$rank < ncol(des)) des <- model.matrix(~ group + logit_mt + log10_assigned, cv)  # e.g. one sex drawn
  yy <- DGEList(cnt[, ids])
  yy <- calcNormFactors(yy[filterByExpr(yy, group = grp), , keep.lib.sizes = FALSE], method = "TMM")
  f <- eBayes(lmFit(voomWithQualityWeights(yy, des), des), robust = TRUE)
  tt <- topTable(f, coef = "groupdisease", number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), t = signif(tt$t, 5), logFC = signif(tt$logFC, 5), adj.P.Val = signif(tt$adj.P.Val, 5))
}
set.seed(20261104)
res <- lapply(seq_len(N_DRAW), function(d) {
  r <- fit_t(c(sample(normal, N_HEALTHY), sample(esld, N_DIS)))
  r$draw <- d
  r
})
write.csv(do.call(rbind, res), gzfile(file.path(out, DRAWS_FILE)), row.names = FALSE)

# same-platform ceiling for the expression Spearman: two random halves of the normal donors (own seed, after the
# draws so they are unchanged)
if (is.na(STRATUM)) {
  set.seed(1)
  h <- sample(normal)
  half <- floor(length(h) / 2)
  write.csv(data.frame(normal_half1_log2cpm = rowMeans(lc[, h[1:half]]), normal_half2_log2cpm = rowMeans(lc[, h[(half + 1):length(h)]])),
            file.path(out, "bulk_normal_halves_logcpm.csv"))
}
cat("draws:", N_DRAW, "| median BH<0.05 genes per draw:",
    median(sapply(res, function(r) sum(r$adj.P.Val < 0.05))), "\n")
