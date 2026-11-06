# stratified

Hepatocyte-only pseudobulks from scRNA-seq and snRNA-seq of the same livers (Andrews et al.: 4 healthy donor livers, GSE185477; 6 PSC/PBC explants, GSE247128), compared disease-vs-healthy against a bulk RNA-seq cohort of isolated hepatocytes analysed in its own disease strata. The bulk cohort is private and not distributed; only its code is here.

| step | script | data | what |
|---|---|---|---|
| 1 | `scripts/10_bulk_strata_limma.R` | private | bulk limma-voom per disease stratum vs normal donors, `~ group + logit(MT) + age + sex + log10(assigned reads)`, quality weights, robust eBayes |
| 2 | `scripts/15_bulk_profiles_and_downsample.R` | private | group mean expression, two donor halves, and 50 random 6-vs-4 refits (the sample-size yardstick) |
| 3 | `scripts/run_strat.sh` → `11_build_sc_sn_strata.py`, `11b_shared_genes.py`, `11c_label_pseudobulk.py`, `12_pb_strata_limma.R`, `13_gsea_strata.py`, `14_figures.py` | public (+ step 1 output) | cell calling (≥ 500 UMIs, ≥ 200 genes, < 50% MT), Leiden per library, reference-based cluster labels (author healthy cell types, leave-one-liver-out), whole-liver and hepatocyte-only pseudobulks, limma-voom PSC/PBC vs healthy per modality, GSEA |
| 4 | `scripts/16a_hep_author_sensitivity.R` | public | healthy side rebuilt from the authors' QC'd hepatocytes |
| 5 | `scripts/16_concordance.py`, `scripts/17_concordance_figures.py` | public + private | pathway-rank Spearman, pathways significant → not significant both ways, top-25 lists, gene-level support; all against the yardstick |

`results/` holds only public-data outputs: hepatocyte label check against the authors' labels (`label_check.tsv`: agreement 83–98%, hepatocyte precision 0.85–1.00, recall 0.84–1.00), per-library cell counts (`libraries.tsv`), per-cluster labels (`cluster_labels.tsv`) and the sc/sn limma summary (`pb_strata_summary.csv`).
