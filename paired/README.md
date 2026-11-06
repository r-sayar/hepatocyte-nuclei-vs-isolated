# liver_pb_gsea

Does inflammation GSEA change when the same livers are bulk-sequenced, or pseudobulked from scRNA-seq or snRNA-seq? Write-up: `docs/Liver inflammation GSEA — bulk vs scRNA vs snRNA pseudobulk.md`. Figures: `figures/liver_pb_gsea/`.

Runs on FU Berlin IMP compute04 in `/srv/public/sayar99/liver_pb_gsea/` (node-local; raw data and per-analysis DESeq2/GSEA tables stay there). This folder holds the code and the summary tables.

| file | what it does |
|---|---|
| `encode_manifest.tsv`, `andrews_urls.txt` | inputs: ENCODE4 liver bulk / snRNA files (4 donors, 5 pairs); GEO URLs for GSE185477 + GSE247128 |
| `scripts/00_gene_lengths.py` | exonic length, span, intron fraction per gene (GENCODE v46) |
| `scripts/01_build_encode.py` | ENCODE bulk (RSEM expected counts) + snRNA pseudobulk, ENSG-keyed |
| `scripts/02_build_andrews.py` | sc + sn pseudobulk per liver; one cell-calling rule plus the authors' cells for the healthy livers |
| `scripts/liver_common.py` | shared helpers; imports `ref/common.py` (verbatim copy of `bulktech-assimilation/scripts/pathway_bulk_vs_pb/common.py`) |
| `scripts/03_gsea_liver.py` | all-genes run: platform tests, PSC/PBC vs healthy on sc and on sn, leave-one-liver-out, author-cell check |
| `scripts/04_protein_coding_length.py` | the same on protein-coding genes, plus the length-adjusted ranking |
| `scripts/05_celltype.py` | sc vs sn within cell types and composition-matched (healthy livers) |
| `scripts/06_figures.py` | `fig1_platform_heatmap.png`, `fig2_contrast_sc_vs_sn.png` |
| `scripts/launch.sh` | starts 03, 04 and 05 under nohup on compute04 |
| `results/` | summary CSVs: `pc_*` = protein-coding (primary), unprefixed = all genes, `ct_*` = cell-type, `fig*_values*` = figure data |
