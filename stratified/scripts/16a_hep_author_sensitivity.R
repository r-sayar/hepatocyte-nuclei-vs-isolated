#!/usr/bin/env Rscript
# 16a_hep_author_sensitivity.R — hepatocyte-only PSC/PBC vs healthy fits with the HEALTHY side built from the
# authors' own QC'd hepatocytes instead of rule-called, cluster-labelled cells (11/11c).
# Why: the >= 500 UMI rule calls 39k / 53k barcodes in Healthy3_sn / Healthy4_sn against ~9k author cells; the extra
# barcodes are ambient-heavy droplets, and liver ambient RNA is mostly hepatocyte RNA, so they can leak into a
# "hepatocyte" pseudobulk. The disease side is unchanged (no author cell lists exist for GSE247128).
# Author-cell pseudobulks per type come from the paired analysis (liver_pb_gsea/results/CT, read only).
# Same engine and filters as 12_pb_strata_limma.R. Output: results/pb_<sc|sn>_hepauthor_limma.csv
suppressPackageStartupMessages({ library(edgeR); library(limma) })
D <- "/srv/public/sayar99/liver_strat_gsea/data"
R <- "/srv/public/sayar99/liver_strat_gsea/results"
AU <- "/srv/public/sayar99/liver_pb_gsea/results/CT/celltype_pseudobulk_ensg.csv"
MIN_HEP <- 50
lib <- read.delim(file.path(D, "libraries.tsv"), stringsAsFactors = FALSE)
shared <- readLines(file.path(D, "shared_ensg.txt"))
ann <- read.delim(file.path(D, "ensg2symbol_v46.tsv"), header = FALSE, col.names = c("ensg", "symbol", "type"))
ann <- ann[!duplicated(ann$ensg), ]; rownames(ann) <- ann$ensg
H <- as.matrix(read.csv(file.path(D, "pb_hep_ensg.csv"), row.names = 1, check.names = FALSE))
A <- as.matrix(read.csv(AU, row.names = 1, check.names = FALSE))
g <- intersect(intersect(rownames(H), rownames(A)), shared)
for (mod in c("sc", "sn")) {
  dis <- lib$library[lib$modality == mod & lib$group != "healthy" & lib$n_hep >= MIN_HEP]
  hea <- paste0(sub(paste0("_", mod, "$"), "", lib$library[lib$modality == mod & lib$group == "healthy"]), "_", mod)
  X <- cbind(A[g, paste0(hea, "__Hepatocyte"), drop = FALSE], H[g, dis, drop = FALSE])
  colnames(X) <- c(hea, dis)
  grp <- factor(c(rep("healthy", length(hea)), rep("PSC_PBC", length(dis))), levels = c("healthy", "PSC_PBC"))
  des <- model.matrix(~ grp)
  y <- DGEList(X)
  y <- calcNormFactors(y[filterByExpr(y, group = grp), , keep.lib.sizes = FALSE], method = "TMM")
  f <- eBayes(lmFit(voomWithQualityWeights(y, des), des), robust = TRUE)
  tt <- topTable(f, coef = "grpPSC_PBC", number = Inf, sort.by = "none")
  tt$symbol <- ann[rownames(tt), "symbol"]
  tt$protein_coding <- ann[rownames(tt), "type"] %in% "protein_coding"
  write.csv(tt, file.path(R, paste0("pb_", mod, "_hepauthor_limma.csv")))
  write.csv(X, gzfile(file.path(D, paste0("pb_", mod, "_hepauthor_counts.csv.gz"))))
  cat(mod, "hepatocyte (author healthy):", length(hea), "healthy vs", length(dis), "PSC/PBC |", nrow(tt), "genes |",
      sum(tt$adj.P.Val < 0.05), "BH<0.05\n")
}
