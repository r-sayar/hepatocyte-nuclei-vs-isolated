# Liver inflammation GSEA — bulk vs scRNA vs snRNA pseudobulk of the same livers

*2026-11-04. Question: if the same liver samples are measured by bulk RNA-seq, or pseudobulked from scRNA-seq or snRNA-seq, does GSEA give different answers, specifically for inflammation? Computed on FU Berlin IMP compute04 (`/srv/public/sayar99/liver_pb_gsea/`). Code in `liver_pb_gsea/scripts/`, summary tables in `liver_pb_gsea/results/`, figures in `figures/liver_pb_gsea/`.*

---

## 0. The answer

**Yes. The library type changes the inflammation result, and scRNA vs snRNA is the bigger split.**

1. **Same livers, scRNA vs snRNA pseudobulk (10 livers).** All 8 Hallmark inflammation sets are higher in scRNA pseudobulk (protein-coding NES +1.5 to +2.9, all p < 0.05). Allograft rejection and IFN-α response are higher in sc in 10 of 10 livers individually. The gap survives every control we tried. It holds within a single cell type (hepatocytes, LSECs, macrophages, NK/T, B cells), after re-weighting sc to the same cell-type mix as sn, with the authors' own QC'd cells instead of our cell calling (NES r 0.998), and after adjusting for gene length.
2. **Only TNFα/NF-κB is dissociation stress.** Removing the warm-dissociation stress genes from the ranking abolishes the TNFα/NF-κB gap in every cell type, e.g. hepatocytes +1.63 → +0.86 (p 0.76). IFN-α, IFN-γ, allograft rejection and complement do not move: hepatocytes IFN-γ +1.77 → +1.74. Their gap is a nuclear-vs-cytoplasmic capture difference. The clearest case is MHC-II. CD74 is 3,800–10,200 CPM in PSC/PBC scRNA pseudobulk but 7–66 CPM in snRNA of the same liver, a 58–538× gap. PTPRC (CD45), which marks the immune cells themselves, differs only 1–9× (median about 6×).
3. **The same disease contrast gives a different inflammation call on sc vs sn.** For PSC/PBC vs healthy on the same 10 livers:
   - **sc pseudobulk** calls all 8 inflammation sets up (p < 0.001).
   - **sn pseudobulk** calls 5 of 8.
   - **IFN-α response** flips from called to not called: sc +2.37 vs sn +0.77 (p 0.89). Dropping any one liver leaves the two apart (sc +2.19 to +2.49, sn −0.88 to +0.98).
   - **IFN-γ response and complement** are also lost in sn.
   - **GO BP:** 120 vs 58 inflammation-related terms are called up.
4. **snRNA vs bulk (ENCODE, 4 healthy donors, 5 donor × lobe pairs).** Inflammation sets come out lower in snRNA pseudobulk (NES −1.1 to −1.8), but so do 34 of 50 Hallmark sets. Inflammation is not singled out: no inflammation set is in the reported bottom-15. Per-pair signs are mixed: TNFα is higher in sn in 4 of 5 pairs. With all annotated genes in the ranking (the default pipeline), 46 of 50 Hallmark sets look "bulk-high", because nuclear libraries load intronic, long and non-coding genes. **Restrict to protein-coding genes before comparing sn with bulk.**

**What this means for us.** Comparing anything sc-derived with anything sn- or bulk-derived will produce inflammation differences, especially IFN and antigen presentation, that come from the library rather than the biology. Within one platform a contrast is safer. But in the one dataset where the same disease contrast was measured both ways, the IFN conclusion depended on whether the livers were profiled as cells or as nuclei.

---

## 1. Data: what "same samples" means here

No public human-liver dataset has bulk RNA-seq and sc/snRNA-seq from the same samples together with an inflammatory contrast. GEO, ArrayExpress/BioStudies and ENCODE were searched (2026-11-04). The two closest designs fail for different reasons. The Innsbruck machine-perfusion livers (bulk GSE263614, sc GSE216584) cannot be linked sample-to-sample from GEO. The Gubra NASH biopsies with bulk and snRNA were never deposited. The internal collections checked hold liver bulk and spatial data only, with no liver snRNA-seq. So the question is answered with two paired designs that each cover one side.

| set | material | libraries | pairing | contrast |
|---|---|---|---|---|
| **ENCODE4** (Snyder snRNA-seq, Wold total RNA-seq) | healthy donor liver, right lobe (+ one left lobe) | bulk total RNA-seq (RSEM expected counts) vs 10x snRNA-seq (STARsolo GeneFull_Ex50pAS, cell-called) | same donor and lobe, separate tissue pieces; 4 donors (16–47 y), 5 pairs | none (platform only) |
| **Andrews 2022** GSE185477 | healthy donor livers C41, C58, C70, C72 | 10x 3′ sc (v2) vs 10x 3′ sn (TST nuclei; v3 by barcode whitelist) | same liver | — |
| **Andrews 2024** GSE247128 | PSC (5) and PBC (1) explants | 10x 5′ sc vs 10x 3′ sn | same explant | PSC/PBC vs healthy (with the 2022 livers) |

Pseudobulk = sum over called cells. The GEO matrices are unfiltered, so every library gets one rule: ≥ 500 UMIs, ≥ 200 genes, < 50% MT, the MT ceiling being the authors'. Per-library counts are in `liver_pb_gsea/results/andrews_libraries.tsv` and `encode_samples.tsv`.

## 2. Method

Same convention as the lab's standard pathway figures and an earlier pathway note (`ref/common.py` is a verbatim copy of `bulktech-assimilation/scripts/pathway_bulk_vs_pb/common.py`):

- **DE:** DESeq2 (PyDESeq2) with the rowSums ≥ 10 prefilter. Platform tests are paired: `~sample + platform`. The disease contrast is `~group`, run separately on sc and on sn pseudobulk.
- **Ranking:** genes ranked by the Wald statistic, MT-* dropped, symbols from GENCODE v46 (joined on ENSG), de-duplicated by baseMean.
- **GSEA:** prerank with 1000 permutations, set size 10–500, seed 1234, MSigDB 2026.1. Hallmark and KEGG legacy use gseapy, which keeps the leading edges. Reactome and GO BP use the project's numpy re-implementation, which reproduces gseapy NES at r = 0.9999.
- **Gene universe:** the primary analysis uses **protein-coding genes**. The all-genes version (script 03) is kept as "what the default pipeline does".
- **Inflammation sets:**
  - Hallmark: TNFA_SIGNALING_VIA_NFKB, INFLAMMATORY_RESPONSE, INTERFERON_ALPHA/GAMMA_RESPONSE, IL6_JAK_STAT3, IL2_STAT5, COMPLEMENT and ALLOGRAFT_REJECTION.
  - Other libraries: any term matching inflammation, interferon, cytokine, chemokine, NF-κB, TNF, interleukin, Toll, complement or acute phase.
- **Controls:**
  - **Length-adjusted ranking:** Wald stat minus the median stat of genes in the same exonic-length × span decile cell.
  - **Stress genes removed:** rankings re-scored without the 133-gene warm-dissociation stress/IEG union (van den Brink 2017 + O'Flanagan 2019).
  - **Per-cell-type pseudobulk:** using the authors' cell labels, healthy livers.
  - **Composition-matched pseudobulk:** sc re-weighted to the same liver's sn cell-type proportions.
  - **Author cells vs our cell calling.**
  - **Leave-one-liver-out** for the disease contrast.

## 3. Results

![Hallmark inflammation sets across platform comparisons](figures/liver_pb_gsea/fig1_platform_heatmap.png)

Columns 1–2 are ENCODE snRNA vs bulk. Columns 3–6 are scRNA vs snRNA of the same livers. Columns 7–8 are the composition controls (healthy livers, author cells).

**Where the inflammation sets rank among all 50 Hallmark sets** (1 = most higher in the first-named library):

| set | sn vs bulk | sc vs sn, 10 livers | sc vs sn, length-adj. | sc vs sn, hepatocytes only |
|---|---|---|---|---|
| Allograft rejection | 17 | **3** | **2** | 13 |
| IFN-γ response | 32 | **4** | **4** | 18 |
| IFN-α response | 35 | **8** | **9** | 19 |
| Complement | 29 | 14 | **5** | 15 |
| TNFα via NF-κB | 20 | 15 | 31 | 23 |
| Inflammatory response | 12 | 16 | 13 | 36 |

In sn vs bulk, the inflammation sets sit mid-pack in a general protein-coding shift. In sc vs sn they are among the most platform-shifted sets in Hallmark.

**Dissociation-stress genes removed** (Hallmark NES, sc vs sn; bold = p < 0.05):

| cell type | TNFα via NF-κB, full ranking → stress genes removed | IFN-γ response, full → removed | allograft rejection, full → removed |
|---|---|---|---|
| hepatocytes | **+1.63** → +0.86 | **+1.77** → **+1.74** | **+1.95** → **+1.97** |
| LSECs | **+1.80** → +1.12 | **+1.79** → **+1.76** | **+1.64** → **+1.66** |
| macrophages | **+1.94** → +1.24 | **+1.84** → **+1.82** | **+2.14** → **+2.18** |
| NK/T cells | **+1.67** → +0.87 | **+1.88** → **+1.85** | **+2.33** → **+2.37** |
| same cell mix, whole liver | **+1.84** → +1.31 | **+1.98** → **+1.97** | **+2.18** → **+2.22** |

![PSC/PBC vs healthy, sc vs sn](figures/liver_pb_gsea/fig2_contrast_sc_vs_sn.png)

**PSC/PBC vs healthy, same 10 livers, protein-coding:**

| Hallmark set | sc pseudobulk NES (p) | sn pseudobulk NES (p) | leave-one-out, sc | leave-one-out, sn |
|---|---|---|---|---|
| Allograft rejection | **+3.04** (<0.001) | **+1.87** (<0.001) | 2.79 – 3.23 | 1.74 – 2.21 |
| IFN-γ response | **+2.78** (<0.001) | +1.27 (0.061) | 2.57 – 2.85 | 1.11 – 1.68 |
| Inflammatory response | **+2.61** (<0.001) | **+2.10** (<0.001) | 2.30 – 2.71 | 2.00 – 2.45 |
| IFN-α response | **+2.37** (<0.001) | +0.77 (0.89) | 2.19 – 2.49 | −0.88 – 0.98 |
| TNFα via NF-κB | **+2.16** (<0.001) | **+1.72** (<0.001) | 1.76 – 2.40 | 1.47 – 2.71 |
| IL6–JAK–STAT3 | **+2.06** (<0.001) | **+1.35** (0.047) | 1.76 – 2.21 | 1.30 – 1.87 |
| IL2–STAT5 | **+1.99** (<0.001) | **+1.68** (<0.001) | 1.80 – 2.23 | 1.50 – 1.96 |
| Complement | **+1.77** (<0.001) | +1.12 (0.22) | 1.50 – 2.00 | −1.19 – 1.39 |

Across all 50 Hallmark sets the two runs agree at NES r = 0.75; 8 of 15 top-15 and 9 of 15 bottom-15 entries are shared. Across the 8 inflammation sets alone, r = 0.34. The gene-level Wald statistics correlate at r = 0.49. Inflammation-related sets called up: Hallmark 8 vs 5, Reactome 27 vs 15, GO BP 120 vs 58.

## 4. Why

- **Nuclei hold far less mature cytoplasmic mRNA of short, highly expressed transcripts.** MHC class I/II, B2M and interferon-stimulated genes are this kind of transcript, and they fill the IFN-γ, IFN-α and allograft sets. The CD74 / HLA-DRA numbers above are the clearest instance. This is a per-gene effect, not a per-cell-type one, so re-weighting composition cannot remove it.
- **Gene length is not the explanation by itself.** The platform statistic does track gene span strongly (Spearman +0.54 for sn vs bulk, −0.65 for sc vs sn: long genes are higher in nuclei). But removing the length trend leaves the inflammation gap intact, and after adjustment 34 of 50 Hallmark sets are higher in sc rather than 26. Transcript localisation and stability beyond length were not tested here; this is the likely remaining factor.
- **TNFα/NF-κB is the one set that is dissociation.** Warm enzymatic dissociation induces FOS/JUN/EGR1/IER/NFKBIA in cells, never in nuclei. The stress list is an assembled union that includes real NF-κB targets (NFKBIA, TNFAIP3, CXCL1/2, SOCS3). So "removing it abolishes TNFα" cannot separate dissociation-induced NF-κB from biological NF-κB activity. It only shows that the sc-vs-sn TNFα gap lives in those genes.
- **snRNA vs bulk is mostly a gene-class shift.** Nuclear counts include introns (GeneFull) while bulk RSEM is exonic. With all genes ranked, every protein-coding pathway drifts to the bulk side. On protein-coding genes the drift remains (34 of 50 Hallmark sets lower in sn), but inflammation is not specifically affected.

## 5. Caveats

- **The disease contrast spans two studies.** Healthy livers are from 2022 (3′ sc, TST nuclei) and PSC/PBC explants from 2024 (5′ sc, 3′ sn), so absolute disease effects are confounded with study on both platforms. The disease sc libraries are also much more immune-rich than the healthy sc libraries: ALB mostly 2–5k vs 12–35k CPM, CD74 3.8–10.2k vs 57–493. PSC019's sc library has no hepatocyte reads at all. This inflates the sc contrast. The claim here is only that **the same livers give a different inflammation answer depending on cells vs nuclei**, not what PSC biology is.
- **Small n.** 4 healthy and 6 disease livers; ENCODE 4 donors (5 pairs). The ENCODE bulk and sn come from separate pieces of the same lobe.
- **Cell calling.** Our rule over-calls nuclei in two healthy sn libraries (C70 TST 39k vs 9k author cells, C72 TST 53k vs 9k), probably ambient-heavy droplets. Using the authors' cells instead gives NES r = 0.998 (Hallmark), 0.98–0.99 (others).
- **No direct sc vs bulk in liver.** That needs a paired liver set that does not exist publicly. The thesis note on non-liver paired sets (bulktech-assimilation, 2026-10-13) found allograft rejection, IFN-α and IFN-γ significantly pseudobulk-high vs bulk in 4, 2 and 1 of 5 datasets, and TNFα and inflammatory response 1 of 5 each way (`scripts/pathway_bulk_vs_pb/results/C1_consensus.csv`).

## 6. Reproduce

On compute04 (`ssh compute04-fu`), project dir `/srv/public/sayar99/liver_pb_gsea/` (node-local disk; `/storage/mi/sayar99` is over quota):

```
venv/bin/python scripts/00_gene_lengths.py            # GENCODE v46 exonic length / span
venv/bin/python scripts/01_build_encode.py            # ENCODE bulk + sn pseudobulk (inputs: encode_manifest.tsv)
venv/bin/python scripts/02_build_andrews.py 8         # Andrews sc + sn pseudobulk (inputs: andrews_urls.txt)
bash scripts/launch.sh                                # 03 all genes, 04 protein-coding + length-adj, 05 cell types
venv/bin/python scripts/06_figures.py
```

Raw downloads: md5-verified against ENCODE, and size-verified against the GEO file lists. Every step caches its DESeq2 and GSEA tables under `results/`.

## References

- Andrews TS et al. 2022, Hepatol Commun — [PubMed 34792289](https://pubmed.ncbi.nlm.nih.gov/34792289/)
- Andrews TS et al. 2024, J Hepatol — [PubMed 38199298](https://pubmed.ncbi.nlm.nih.gov/38199298/)
- MacParland SA et al. 2018, Nat Commun (C41 sc library) — [PubMed 30348985](https://pubmed.ncbi.nlm.nih.gov/30348985/)
- ENCODE Project Consortium 2012, Nature — [PubMed 22955616](https://pubmed.ncbi.nlm.nih.gov/22955616/); Moore JE et al. 2020, Nature (ENCODE phase III) — [PubMed 32728249](https://pubmed.ncbi.nlm.nih.gov/32728249/)
- van den Brink SC et al. 2017, Nat Methods — [PubMed 28960196](https://pubmed.ncbi.nlm.nih.gov/28960196/)
- O'Flanagan CH et al. 2019, Genome Biol — [PubMed 31623682](https://pubmed.ncbi.nlm.nih.gov/31623682/)
