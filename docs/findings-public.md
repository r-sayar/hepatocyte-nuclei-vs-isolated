# Findings from public data (2026-11-05)

Three findings about hepatocyte-nucleus pseudobulks that matter for anyone comparing them with bulk RNA-seq, or reading disease pathways off them. Data: Gribben et al. 2024 (MASLD spectrum, GSE202379, PMID 38778114) and Andrews et al. 2022/2024 (healthy and PSC/PBC livers, GSE185477 / GSE247128, PMIDs 34792289 / 38199298). Hepatocyte pseudobulks per patient, limma-voom, GSEA as in the rest of the repo. Numbers: `etiology_matched/results/`.

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
