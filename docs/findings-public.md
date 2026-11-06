# Findings from public data (2026-11-05)

Data: Gribben et al. 2024 (MASLD spectrum, GSE202379, PMID 38778114) and Andrews et al. 2022/2024 (healthy and PSC/PBC livers, GSE185477 / GSE247128, PMIDs 34792289 / 38199298); for replication, two MASLD liver-biopsy RNA-seq cohorts (GSE162694; GSE130970). Hepatocyte pseudobulks per patient, limma-voom, GSEA as in the rest of the repo. Numbers: `etiology_matched/results/`.

**The main result is section 4.** A hepatocyte programme of fibrosis progression was discovered in contamination-corrected hepatocyte nuclei of two diseases (MASH and PSC/PBC). It replicates, with the same slope, in two independent biopsy cohorts (190 patients). Sections 1–3 are the groundwork that made the nuclei trustworthy enough to discover it.

## 4. A hepatocyte programme of fibrosis progression: discovered in nuclei, replicated in two biopsy cohorts

**Discovery, nuclei only.** Genes moving the same way (t > 2.5) in both decontaminated hepatocyte-nucleus contrasts:

- the MASH fibrosis trend within one sample source (GSE202379, 40 patients, F0–F4);
- PSC/PBC vs healthy (GSE247128 / GSE185477, 6 vs 4 livers).

That gives 76 genes up and 64 down, against 28 and 38 expected by chance. Restricted to genes hepatocytes make (MacParland-atlas hepatocyte fraction ≥ 0.5), the programme is **17 genes up and 47 down** (`results/nuclei_hepatocyte_programme.csv`; `scripts/24a_define_nuclei_programme.py`). No bulk data was used to choose them.

**Replication in whole-biopsy RNA-seq.** Independent patients. Score = mean z of the up genes minus the down genes, modelled as score ~ fibrosis stage + hepatocyte content + age + sex. The hepatocyte-content term is a 10-gene hepatocyte marker score, so the result is not "fewer hepatocytes". The null is 1,000 random sets of hepatocyte-made genes matched on expression (`scripts/25_programme_in_mash_cohorts.py`):

| cohort | NAFLD patients | slope per fibrosis stage | t | p vs matched random sets |
|---|---|---|---|---|
| GSE162694 | 112 (+ 31 normal) | 0.454 | 8.6 | < 0.001 |
| GSE130970 | 78 | 0.454 | 5.9 | 0.003 |

The two slopes agree to three decimals. Dose-response, mean score:

| | normal | F0 | F1 | F2 | F3 | F4 |
|---|---|---|---|---|---|---|
| GSE162694 | −0.45 (31) | −0.39 (35) | −0.26 (30) | +0.25 (27) | +0.60 (8) | +2.01 (12) |
| GSE130970 | — | −0.51 (25) | −0.18 (28) | +0.29 (9) | +0.83 (14) | +1.82 (2) |

It is flat through F0–F1, rises from F2 and jumps at cirrhosis.

**Per gene:** 40 of the 64 genes replicate on their own (t > 2 in the programme direction in both cohorts; `results/programme_genes_in_mash_cohorts.csv`).

- **Up (9 of 17):** SPSB1, PHLDA1, AKR1B10, FAT1, TNFRSF11B, GDF15, TMEM92, CLDN7, NEDD4L.
- **Down (31 of 47),** strongest first: ADCY1, KCNN2, PPP1R1A, COQ10A, ENPP3, ALDH2, CTNNA3, SULT1E1, AMOTL1, TCEA3, MME, CNDP1.

**What is new.** PubMed title/abstract searches, run 2026-11-05, for the gene plus "liver fibrosis / hepatic fibrosis / cirrhosis / NAFLD / NASH / MASLD / MASH", and for the gene plus "hepatocyte":

| gene | with liver fibrosis / cirrhosis / MASLD | with hepatocyte |
|---|---|---|
| GDF15 | 123 | 85 |
| AKR1B10 | 44 | 25 |
| NEDD4L | 9 | 8 |
| PHLDA1 | 6 | 3 |
| CLDN7 | 1 | 4 |
| **SPSB1** | **0** | **0** |
| **TMEM92** | **0** | **0** |
| ADCY1 | 3 | 0 |
| PPP1R1A | 1 | 0 |
| KCNN2 | 2 | 0 |
| CTNNA3 / TCEA3 / KANK4 | 0 | 0 |

GDF15 and AKR1B10 are known MASH hepatocyte markers and act as positive controls. Most of the programme has no prior link to fibrosis or to hepatocytes. Two themes stand out; both are hypotheses for follow-up, not established mechanisms:

- **Hepatocytes switch on brakes on TGF-β signalling as fibrosis advances.** SPSB1 targets the TGF-β type II receptor for degradation, and NEDD4L ubiquitinates SMAD2/3 and the type I receptor. With them rise the decoy receptor osteoprotegerin (TNFRSF11B), the stress hepatokine GDF15 and PHLDA1: a graded self-protection programme.
- **Hepatocytes lose their cAMP-signalling apparatus.** ADCY1 (adenylyl cyclase) and PPP1R1A (inhibitor-1, the cAMP/PKA amplifier) are the strongest falling genes in both cohorts, with the K⁺ channels KCNN2 and ABCC9. Glucagon acts through cAMP, so this is a candidate cell-intrinsic route to the hepatic glucagon resistance described in MASLD (to be checked against the literature with verified citations; research-plan task T12).

**Blood-measurable members:** osteoprotegerin and GDF15 rise; butyrylcholinesterase (BCHE) and carnosinase (CNDP1) fall. Together they could read out the programme without a biopsy (task T11).

**Caveats.**

- **Biopsy bulk is not hepatocytes.** Adjusting for hepatocyte content removes the main composition effect but not all of it. The atlas restriction and the matched null guard against generic hepatocyte loss.
- **Cross-sectional design**, so the programme marks progression; it is not shown to drive it.
- **Discovery thresholds** (t > 2.5 in each nuclei contrast) are nominal; replication, not discovery, carries the inference.
- **Few F4 patients** in GSE130970 (2).

## 1. In GSE202379 the healthy controls come from two sample sources, and the end-stage contrast inherits it

Patients with HL*/CL* IDs (all 5 end-stage explants and 2 of the 4 healthy controls) separate from numeric-ID patients (all NAFLD and NASH, the other 2 healthy controls) on PC2 of the hepatocyte pseudobulks (18% of variance; group means −25 and −40 for HL/CL healthy and end stage, −1 to +25 for the numeric-ID groups). Each disease group sits in one source only.

The consequence shows up in the most basic hepatocyte readout. End-stage MASH hepatocyte nuclei vs healthy controls give *higher* xenobiotic metabolism, fatty-acid metabolism, OXPHOS and cholesterol homeostasis (Hallmark NES +1.5 to +1.9, p ≤ 0.004, with the source covariate in the model). That is the opposite of the metabolic decline seen in PSC/PBC hepatocyte nuclei and cells, and expected in end-stage disease. A source covariate cannot rescue it, because only two healthy controls exist per source.

A within-source design does rescue it. The fibrosis-stage trend inside the biopsy source (40 patients, F0–F4, no dependence on the healthy controls) puts bile-acid, xenobiotic and fatty-acid metabolism back on the decline side (NES −1.0 to −1.2, not significant), while cholesterol homeostasis rises with stage (+1.4, p = 0.026). Recommendation for users of this atlas: use stage trends within a source, not healthy-vs-disease contrasts.

## 2. Ambient RNA is 12–33% of a "hepatocyte" nucleus pseudobulk and rises with fibrosis

The contamination fraction is estimated per patient from immune and endothelial genes that hepatocytes do not express, with the all-nuclei pseudobulk as the ambient pool (SoupX logic at pseudobulk level). It is a median 17% in healthy controls, 21% in NAFLD, 25% in NASH, 33% in NASH with cirrhosis and 28% in end-stage livers (Gribben et al.), and 12–16% in the Andrews livers. The ambient pool explains a median 62–90% of the collagen/ECM signal (COL1A1, COL1A2, COL3A1, DCN, LUM, ACTA2, PDGFRB, SPARC, COL6A3, FBN1) found in hepatocyte nuclei.

Because the fraction grows with fibrosis, it manufactures apparent hepatocyte induction of stromal programmes along disease progression. After subtraction, the MASH fibrosis trend loses Hallmark myogenesis (+1.6 → +1.0), angiogenesis (+1.8 → +1.2) and hypoxia (+1.4 → +1.2). In PSC/PBC nuclei the immune side moves instead: allograft rejection +2.0 → +1.2, complement +1.0 → −1.1.

## 3. What is left of "hepatocyte EMT" is a hepatocyte injury programme, not a mesenchymal transition

Hallmark EMT survives decontamination in both nuclei datasets (MASH fibrosis trend +2.2 → +1.7; PSC/PBC +2.3 → +2.0), but its content changes.

- **Before correction,** the MASH trend is led by stromal genes (COL1A1, COL3A1, DCN, VIM, FBLN5, IGFBP3; MacParland-atlas hepatocyte fraction ≤ 0.17).
- **After correction,** the leaders are genes hepatocytes make themselves: SERPINE1 (atlas hepatocyte fraction 0.90), PLOD2 (0.90), COL5A3 (0.97), HTRA1 (0.85), CAP2 (0.82), TNFRSF12A/Fn14 (0.60), PVR (0.82). They come with TNF/NF-κB targets (IL32, CXCL8, TNFAIP3, JUN).
- **No canonical EMT transcription factor or mesenchymal marker** (SNAI1/2, ZEB1/2, TWIST1, CDH2) is among the leaders.

The "EMT" gene-set call in human hepatocyte nuclei is therefore ambient stromal RNA plus a hepatocyte wound / secretory response that shares genes with the EMT set. This is consistent with mouse lineage tracing that finds no hepatocyte EMT in fibrosis. Calling it EMT from GSEA alone is not supported.

## Caveats

- The ambient pool is approximated by the all-nuclei pseudobulk because no empty droplets are deposited for GSE202379. An empty-droplet profile (Andrews raw matrices) would sharpen it (research-plan task T3).
- The contamination fraction treats doublets of hepatocyte and non-hepatocyte nuclei as contamination. That is right for pseudobulk purposes but means it is an upper bound on true ambient RNA.
- The healthy controls are four per dataset.
