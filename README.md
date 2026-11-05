# Hepatocyte nuclei vs isolated hepatocytes

Do single-nucleus RNA-seq pseudobulks of hepatocytes and bulk RNA-seq of freshly isolated hepatocytes tell the same story about liver disease? And where they disagree, is one of them showing an artifact — or is there biology that only one of them can see?

The aim is to find something new that the usual readout hides, or to show that an accepted result comes from an artifact (isolation stress, carry-over of other cell types, ambient RNA in droplets, nuclear vs cytoplasmic RNA, sample size).

## Why this is hard: no paired data

No public human-liver dataset has bulk RNA-seq of isolated hepatocytes and sn/scRNA-seq from the same livers (searched GEO, ArrayExpress/BioStudies and ENCODE, 2026-11-04; see `docs/paired-analysis-note.md`). So the comparison is **stratified**: each platform is analysed within its own disease-vs-healthy strata, and only the answers are compared. Four things make that trustworthy:

1. **Same material.** Hepatocytes are taken from the sn/sc data by clustering and reference-based labelling (checked against the original authors' labels), so hepatocytes are compared with hepatocytes.
2. **Same statistics.** The same limma-voom engine and the same GSEA settings on every side.
3. **A sample-size yardstick.** The bulk is cut down to the pseudobulk's size (random 6-vs-4 draws) and scored identically. A difference only counts if it exceeds what an equally small bulk shows.
4. **A within-dataset control.** scRNA and snRNA of the *same* livers: what cells and nuclei share is cohort or disease, what nuclei do alone is the nucleus.

## What is here

| folder | what | data |
|---|---|---|
| `paired/` | Two paired designs that each cover one side: ENCODE4 liver (bulk vs snRNA, same 4 healthy donors) and Andrews et al. livers (scRNA vs snRNA of the same 4 healthy livers and 6 PSC/PBC explants). Write-up: `docs/paired-analysis-note.md` | public |
| `stratified/` | Hepatocyte-only pseudobulks from the Andrews livers, compared by disease stratum with a bulk RNA-seq cohort of isolated hepatocytes | sn/sc public; the bulk cohort is private (see below) |
| `docs/` | notes, research plan, open questions for the original authors | — |

**The isolated-hepatocyte bulk cohort is unpublished and not in this repo.** The code that analyses it is included (`stratified/scripts/10_*.R`, `15_*.R`) but needs local inputs; nothing derived from it (numbers, tables, figures) is committed here.

## Public data used

- Andrews et al. 2022, Hepatol Commun, PMID 34792289 — GSE185477 (healthy livers, sc + nuclei, author cell labels).
- Andrews et al. 2024, J Hepatol, PMID 38199298 — GSE247128 (PSC/PBC explants, sc 5′ + sn 3′, raw matrices).
- ENCODE4 liver (Snyder / Wold labs): snRNA-seq and total RNA-seq of the same donors (accessions in `paired/encode_manifest.tsv`).
- MacParland et al. 2018, Nat Commun, PMID 30348985 — healthy-liver atlas (hepatocyte fraction per gene).
- MSigDB 2026.1 (Hallmark, KEGG legacy, Reactome, GO BP).

Download lists: `paired/andrews_urls.txt`, `paired/encode_manifest.tsv`.

## Running

Python 3.13 with scanpy, gseapy 1.3.2, pydeseq2 0.5.4, leidenalg; R with limma/edgeR. Scripts were run on a 56-core Linux box; paths (`/srv/public/...`) are constants at the top of each script, and `GMT_DIR` points to the MSigDB GMT files. Order: `paired/scripts/00 → 06`; `stratified/scripts/run_strat.sh` (11 → 14), then `16a`, `16`, `17`.

## Status (2026-11-05)

- Paired analyses done (see `docs/paired-analysis-note.md`).
- Stratified pipeline done; its comparisons with the private cohort are kept locally.
- Next steps and tasks that can run on public data alone: `docs/research-plan.md`.
- Questions for the original authors: `docs/open-questions-for-researchers.md`.
