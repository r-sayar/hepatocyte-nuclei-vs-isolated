#!/usr/bin/env Rscript
# 21_gribben_limma.R — disease vs healthy in Gribben et al. 2024 (GSE202379) patient pseudobulks, hepatocyte nuclei
# and all nuclei, with the engine used everywhere else in this repo: filterByExpr, TMM, voomWithQualityWeights,
# robust eBayes, design ~ group. Contrasts vs the 4 healthy controls:
#   cirrMASH  = NASH with cirrhosis + end stage   (the match for end-stage MASH isolated-hepatocyte preps)
#   endstage  = end stage only (explants)
#   NASH      = NASH w/o cirrhosis
#   NAFLD     = NAFLD
# Output: <out>/gribben_<hep|all>_<contrast>_limma.csv (symbol rows, protein_coding flag from GENCODE v46)
# usage: Rscript 21_gribben_limma.R <dir with pb_hep.csv, pb_all.csv, patients.csv> <ensg2symbol_v46.tsv> <out>
suppressPackageStartupMessages({ library(edgeR); library(limma) })
a <- commandArgs(trailingOnly = TRUE)
D <- a[1]; out <- a[3]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
ann <- read.delim(a[2], header = FALSE, col.names = c("ensg", "symbol", "type"))
pc <- unique(ann$symbol[ann$type == "protein_coding"])
pt <- read.csv(file.path(D, "patients.csv"), stringsAsFactors = FALSE)
pt$Patient.ID <- as.character(pt$Patient.ID)
CONTRASTS <- list(cirrMASH = c("NASH with cirrhosis", "end stage"), endstage = "end stage",
                  NASH = "NASH w/o cirrhosis", NAFLD = "NAFLD")
# Sample source: patients with HL*/CL* IDs (all end-stage explants and 2 of the 4 healthy controls) separate from
# the numeric-ID patients (all NAFLD / NASH and the other 2 healthy controls) on PC2 of the hepatocyte
# pseudobulks. Every disease group sits in one source only, so an unadjusted disease-vs-healthy contrast partly
# compares sources. The primary fit therefore uses all patients with source as a covariate (healthy controls exist
# in both sources, which identifies it); the unadjusted per-contrast fits are kept as *_unadjusted_limma.csv.
pt$source <- factor(ifelse(grepl("^(HL|CL)", pt$Patient.ID), "HL_CL", "numeric"))
pt$grp <- factor(c("Healthy control" = "healthy", "NAFLD" = "NAFLD", "NASH w/o cirrhosis" = "NASH",
                   "NASH with cirrhosis" = "NASHcirr", "end stage" = "endstage")[pt$Disease.status],
                 levels = c("healthy", "NAFLD", "NASH", "NASHcirr", "endstage"))
fit_and_write <- function(y, des, coef_or_contrast, file, label, n_dis, n_h) {
  f <- lmFit(voomWithQualityWeights(y, des), des)
  if (is.matrix(coef_or_contrast)) f <- contrasts.fit(f, coef_or_contrast)
  f <- eBayes(f, robust = TRUE)
  tt <- topTable(f, coef = if (is.matrix(coef_or_contrast)) 1 else coef_or_contrast, number = Inf, sort.by = "none")
  tt$symbol <- rownames(tt)
  tt$protein_coding <- rownames(tt) %in% pc
  write.csv(tt, file.path(out, file))
  data.frame(fit = label, n_disease = n_dis, n_healthy = n_h, n_genes = nrow(tt),
             n_protein_coding = sum(tt$protein_coding), bh05 = sum(tt$adj.P.Val < 0.05))
}
summ <- list()
for (material in c("hep", "all", "hep_decont")) {  # hep_decont: 23_decontaminate.py
  X <- as.matrix(read.csv(file.path(D, paste0("pb_", material, ".csv")), row.names = 1, check.names = FALSE))
  # primary: all patients, ~ 0 + group + source
  ids <- pt$Patient.ID
  des <- model.matrix(~ 0 + grp + source, pt)
  colnames(des) <- sub("^grp", "", colnames(des))
  y <- DGEList(X[, ids])
  y <- calcNormFactors(y[filterByExpr(y, group = pt$grp), , keep.lib.sizes = FALSE], method = "TMM")
  for (k in names(CONTRASTS)) {
    g <- c(cirrMASH = "(NASHcirr+endstage)/2", endstage = "endstage", NASH = "NASH", NAFLD = "NAFLD")[[k]]
    cm <- makeContrasts(contrasts = paste0(g, "-healthy"), levels = des)
    summ[[paste(material, k, "adj")]] <- cbind(material = material, contrast = k,
      fit_and_write(y, des, cm, paste0("gribben_", material, "_", k, "_limma.csv"), "source-adjusted, all patients",
                    sum(pt$Disease.status %in% CONTRASTS[[k]]), sum(pt$grp == "healthy")))
  }
  # fibrosis trend inside the numeric-ID (biopsy) source only: one sample source, no dependence on the healthy
  # controls; slope per fibrosis stage F0-F4 (NASH with cirrhosis = F4 belongs to this source)
  ft <- pt[pt$source == "numeric" & !is.na(pt$Fibrosis.score..F0.4.), ]
  dft <- model.matrix(~ fib, data.frame(fib = as.numeric(ft$Fibrosis.score..F0.4.)))
  yft <- DGEList(X[, ft$Patient.ID])
  yft <- calcNormFactors(yft[filterByExpr(yft, design = dft), , keep.lib.sizes = FALSE], method = "TMM")
  summ[[paste(material, "fibrosis_trend")]] <- cbind(material = material, contrast = "fibrosis_trend",
    fit_and_write(yft, dft, "fib", paste0("gribben_", material, "_fibrosis_trend_limma.csv"),
                  "slope per fibrosis stage, numeric-ID source", nrow(ft), NA))
  # unadjusted, two-group fits (kept for the record)
  for (k in names(CONTRASTS)) {
    ids <- pt$Patient.ID[pt$Disease.status %in% c("Healthy control", CONTRASTS[[k]])]
    grp <- factor(ifelse(pt$Disease.status[match(ids, pt$Patient.ID)] == "Healthy control", "healthy", "disease"),
                  levels = c("healthy", "disease"))
    d2 <- model.matrix(~ grp)
    y2 <- DGEList(X[, ids])
    y2 <- calcNormFactors(y2[filterByExpr(y2, group = grp), , keep.lib.sizes = FALSE], method = "TMM")
    summ[[paste(material, k, "unadj")]] <- cbind(material = material, contrast = k,
      fit_and_write(y2, d2, "grpdisease", paste0("gribben_", material, "_", k, "_unadjusted_limma.csv"), "unadjusted",
                    sum(grp == "disease"), sum(grp == "healthy")))
  }
}
S <- do.call(rbind, summ)
write.csv(S, file.path(out, "gribben_limma_summary.csv"), row.names = FALSE)
print(S, row.names = FALSE)
