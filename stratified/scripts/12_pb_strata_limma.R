#!/usr/bin/env Rscript
# 12_pb_strata_limma.R — PSC/PBC (6) vs healthy (4) livers on scRNA and snRNA pseudobulks, whole-tissue and
# hepatocyte-only, with the same limma engine as the bulk strata (10_bulk_strata_limma.R):
# filterByExpr, TMM, voomWithQualityWeights, robust eBayes. Design ~ group (no clinical covariates are
# available for these livers). Genes restricted to ENSG annotated in every library (11b). Hepatocyte-only
# fits use libraries with >= MIN_HEP labelled hepatocytes (11c).
# Output: results/pb_<modality>_<material>_limma.csv (per-gene statistics, ENSG rows, symbol + type columns)
suppressPackageStartupMessages({ library(edgeR); library(limma) })
D <- "/srv/public/sayar99/liver_strat_gsea/data"
R <- "/srv/public/sayar99/liver_strat_gsea/results"
dir.create(R, showWarnings = FALSE)
lib <- read.delim(file.path(D, "libraries.tsv"), stringsAsFactors = FALSE)
shared <- readLines(file.path(D, "shared_ensg.txt"))
ann <- read.delim(file.path(D, "ensg2symbol_v46.tsv"), header = FALSE, col.names = c("ensg", "symbol", "type"))
ann <- ann[!duplicated(ann$ensg), ]
rownames(ann) <- ann$ensg

MIN_HEP <- 50  # libraries with fewer labelled hepatocytes leave the hepatocyte-only fit
summ <- list()
for (material in c("whole", "hep", "hepdecont")) {  # hepdecont: etiology_matched/scripts/23_decontaminate.py
  X <- as.matrix(read.csv(file.path(D, paste0("pb_", material, "_ensg.csv")), row.names = 1, check.names = FALSE))
  X <- X[intersect(rownames(X), shared), ]
  for (mod in c("sc", "sn")) {
    L <- lib[lib$modality == mod, ]
    if (material %in% c("hep", "hepdecont")) {
      low <- L$library[L$n_hep < MIN_HEP]
      if (length(low)) cat(mod, "hepatocyte arm: left out (<", MIN_HEP, "hepatocytes):", low, "\n")
      L <- L[L$n_hep >= MIN_HEP, ]
    }
    grp <- factor(ifelse(L$group == "healthy", "healthy", "PSC_PBC"), levels = c("healthy", "PSC_PBC"))
    des <- model.matrix(~ grp)
    y <- DGEList(X[, L$library])
    y <- calcNormFactors(y[filterByExpr(y, group = grp), , keep.lib.sizes = FALSE], method = "TMM")
    f <- eBayes(lmFit(voomWithQualityWeights(y, des), des), robust = TRUE)
    tt <- topTable(f, coef = "grpPSC_PBC", number = Inf, sort.by = "none")
    tt$symbol <- ann[rownames(tt), "symbol"]
    tt$protein_coding <- ann[rownames(tt), "type"] %in% "protein_coding"
    key <- paste0(mod, "_", material)
    write.csv(tt, file.path(R, paste0("pb_", key, "_limma.csv")))
    summ[[key]] <- data.frame(set = key, n_psc_pbc = sum(grp == "PSC_PBC"), n_healthy = sum(grp == "healthy"),
                              n_genes = nrow(tt), n_protein_coding = sum(tt$protein_coding),
                              bh05 = sum(tt$adj.P.Val < 0.05))
    cat(key, ":", nrow(tt), "genes,", sum(tt$adj.P.Val < 0.05), "BH<0.05\n")
  }
}
write.csv(do.call(rbind, summ), file.path(R, "pb_strata_summary.csv"), row.names = FALSE)
