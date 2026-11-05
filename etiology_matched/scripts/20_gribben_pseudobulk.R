#!/usr/bin/env Rscript
# 20_gribben_pseudobulk.R — patient-level pseudobulks from the MASLD-spectrum snRNA-seq atlas of Gribben et al.
# 2024, Nature (PMID 38778114), GSE202379: the etiology-matched (MASH) comparator for isolated-hepatocyte bulk.
#
# Input: GSE202379_SeuratObject_AllCells.rds.gz from the GEO supplement. It is gzipped twice: read it through
# gzcon(gzfile(...)). The authors' labels are used as deposited: cell.annotation (cell type), Disease.status
# (Healthy control, NAFLD, NASH w/o cirrhosis, NASH with cirrhosis, end stage), Patient.ID.
# Output (same directory): pb_hep.csv and pb_all.csv (genes x patients, RNA-assay counts summed over hepatocyte
# nuclei / all nuclei; patients with several libraries are summed), patients.csv, celltype_counts.csv
suppressPackageStartupMessages({ library(SeuratObject); library(Matrix) })
D <- commandArgs(trailingOnly = TRUE)[1]
o <- readRDS(gzcon(gzfile(file.path(D, "GSE202379_SeuratObject_AllCells.rds.gz"), "rb")))
md <- o@meta.data
cnt <- tryCatch(LayerData(o, assay = "RNA", layer = "counts"),
                error = function(e) GetAssayData(o, assay = "RNA", slot = "counts"))
stopifnot(ncol(cnt) == nrow(md))
cat("nuclei:", ncol(cnt), "| genes:", nrow(cnt), "\n")
print(table(md$cell.annotation))
print(table(md$Disease.status[!duplicated(md$Patient.ID)]))

sum_by <- function(cols, groups) {
  g <- factor(groups)
  G <- sparseMatrix(i = seq_along(g), j = as.integer(g), x = 1, dims = c(length(g), nlevels(g)), dimnames = list(NULL, levels(g)))
  as.matrix(cnt[, cols, drop = FALSE] %*% G)  # groups are already the subset's labels
}
hep_label <- grep("^hepatocyte", levels(factor(md$cell.annotation)), ignore.case = TRUE, value = TRUE)
cat("hepatocyte label(s):", hep_label, "\n")
is_hep <- md$cell.annotation %in% hep_label
pid <- as.character(md$Patient.ID)
write.csv(sum_by(which(is_hep), pid[is_hep]), file.path(D, "pb_hep.csv"))
write.csv(sum_by(seq_len(ncol(cnt)), pid), file.path(D, "pb_all.csv"))

keep <- intersect(c("Patient.ID", "Disease.status", "Fibrosis.score..F0.4.", "SAF.Score", "Sex", "sex", "Age", "age"), colnames(md))
pt <- md[!duplicated(md$Patient.ID), keep]
pt$n_nuclei <- as.integer(table(pid)[as.character(pt$Patient.ID)])
pt$n_hep <- as.integer(table(factor(pid[is_hep], levels = unique(pid)))[as.character(pt$Patient.ID)])
pt$n_libraries <- as.integer(tapply(md$orig.ident, pid, function(x) length(unique(x)))[as.character(pt$Patient.ID)])
write.csv(pt, file.path(D, "patients.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(table(pid, md$cell.annotation)), file.path(D, "celltype_counts.csv"))
cat("patients:", nrow(pt), "| hepatocyte nuclei:", sum(is_hep), "\n")
