# etiology_matched

MASH hepatocyte nuclei from Gribben et al. 2024, Nature, "Acquisition of epithelial plasticity in human chronic liver disease" (PMID 38778114, GSE202379): 99,809 nuclei, 47 patients, author labels (69,426 hepatocyte nuclei). Used as the etiology-matched comparator for isolated-hepatocyte bulk, next to the PSC/PBC nuclei of Andrews et al.

| script | what |
|---|---|
| `scripts/20_gribben_pseudobulk.R` | patient pseudobulks (hepatocyte nuclei, all nuclei) from the deposited Seurat object (gzipped twice: read through `gzcon(gzfile())`); `SeuratObject` is enough |
| `scripts/21_gribben_limma.R` | limma-voom, all patients with sample source as covariate (primary), the fibrosis-stage trend inside the biopsy source, and unadjusted two-group fits for the record |
| `scripts/22_etiology_matched_concordance.py` | pathway-level comparison with an isolated-hepatocyte bulk cohort (private; code only) and between the two nuclei datasets |
| `scripts/23_decontaminate.py` | SoupX-style ambient correction at pseudobulk level for both nuclei datasets |
| `scripts/24a_define_nuclei_programme.py` | the hepatocyte progression programme, defined from decontaminated nuclei only |
| `scripts/24_programme_in_bulk.py` | tests the programme on an isolated-hepatocyte bulk cohort (private; code only) |
| `scripts/25_programme_in_mash_cohorts.py` | replication along fibrosis stage in public biopsy cohorts GSE162694 and GSE130970, with an expression-matched null |

`results/` (public data only): `gribben_patients.csv` (per-patient stage, fibrosis, nuclei counts), `gribben_limma_summary.csv`, `decont_ecm_check.csv` (per-sample contamination fraction and how much of the ECM signal it explains), `public_hallmark_hepatocyte_nuclei.csv` (Hallmark NES / p / FDR for every nuclei and cell contrast, before and after decontamination), `nuclei_hepatocyte_programme.csv`, `programme_in_mash_cohorts.csv`, `programme_by_stage.csv`, `programme_genes_in_mash_cohorts.csv`. Findings: `docs/findings-public.md`.
