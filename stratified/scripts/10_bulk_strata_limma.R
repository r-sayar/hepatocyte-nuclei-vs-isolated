#!/usr/bin/env Rscript
# 10_bulk_strata_limma.R — private bulk cohort (isolated normal vs ESLD hepatocytes; not distributed
# "20251006_Bulk RNAseq_Isolated Normal and ESLD Hepatocytes"): per-stratum disease vs normal-donor
# moderated t, for the stratified bulk / scRNA / snRNA inflammation-GSEA comparison.
#
# Same model as the cohort's reference limma model:
#   counts  : htseq, GENCODE v44, MT-* and NUMT pseudogenes removed (precomputed counts_noMT_noNUMT)
#   set     : adult normal donors + ESLD preps inside the donor MT range (precomputed roles; see stratified/README.md)
#   model   : voomWithQualityWeights, robust eBayes, ~ group + logit(MT) + age + sex + log10(assigned)
#   universe: edgeR::filterByExpr per stratum
# Strata (each vs the same normal donors): all in-range ESLD, ALD, MASH, Child-Pugh B, Child-Pugh C.
# Output is per-gene statistics only (no per-patient data): <out>/bulk_<stratum>_limma.csv
#
# usage: Rscript 10_bulk_strata_limma.R <dir with meta.csv, counts_noMT_noNUMT.csv.gz, gencode_v44_genes.csv> <out>
suppressPackageStartupMessages({ library(edgeR); library(limma) })
a <- commandArgs(trailingOnly = TRUE)
inp <- a[1]; out <- a[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

cnt  <- as.matrix(read.csv(gzfile(file.path(inp, "counts_noMT_noNUMT.csv.gz")), row.names = 1, check.names = FALSE))
meta <- read.csv(file.path(inp, "meta.csv"), row.names = 1, stringsAsFactors = FALSE)
gt   <- read.csv(file.path(inp, "gencode_v44_genes.csv"), stringsAsFactors = FALSE)
meta <- meta[colnames(cnt), ]

normal <- rownames(meta)[meta$role == "normal_donor"]
esld   <- rownames(meta)[meta$role == "esld" & meta$esld_in_mt_range == "True"]
cp     <- as.character(meta$child_pugh)
strata <- list(
  esld_all = esld,
  ALD      = esld[meta[esld, "etiology_group"] %in% "ALD"],
  MASH     = esld[meta[esld, "etiology_group"] %in% "MASH"],
  CP_B     = esld[cp[match(esld, rownames(meta))] %in% "B"],
  CP_C     = esld[cp[match(esld, rownames(meta))] %in% "C"])

pc <- unique(gt$gene_name[gt$gene_type == "protein_coding"])
summ <- list()
for (s in names(strata)) {
  ids <- c(normal, strata[[s]])
  grp <- factor(ifelse(ids %in% normal, "normal", "disease"), levels = c("normal", "disease"))
  m   <- meta[ids, ]
  mt  <- pmin(pmax(m$mt_pct / 100, 1e-4), 1 - 1e-4)
  cv  <- data.frame(group = grp, logit_mt = log(mt / (1 - mt)), age = m$age, sex = factor(m$sex),
                    log10_assigned = log10(m$assigned))
  des <- model.matrix(~ group + logit_mt + age + sex + log10_assigned, cv)
  y   <- DGEList(cnt[, ids])
  y   <- calcNormFactors(y[filterByExpr(y, group = grp), , keep.lib.sizes = FALSE], method = "TMM")
  f   <- eBayes(lmFit(voomWithQualityWeights(y, des), des), robust = TRUE)
  tt  <- topTable(f, coef = "groupdisease", number = Inf, sort.by = "none")
  tt$protein_coding <- rownames(tt) %in% pc
  write.csv(tt, file.path(out, paste0("bulk_", s, "_limma.csv")))
  summ[[s]] <- data.frame(stratum = s, n_disease = length(strata[[s]]), n_normal = length(normal),
                          n_genes = nrow(tt), n_protein_coding = sum(tt$protein_coding),
                          bh05 = sum(tt$adj.P.Val < 0.05),
                          median_mt_disease = median(meta[strata[[s]], "mt_pct"]),
                          median_mt_normal = median(meta[normal, "mt_pct"]))
  cat(s, ":", length(strata[[s]]), "vs", length(normal), "|", nrow(tt), "genes |", sum(tt$adj.P.Val < 0.05), "BH<0.05\n")
}
write.csv(do.call(rbind, summ), file.path(out, "bulk_strata_summary.csv"), row.names = FALSE)
