# Research plan

Goal: find a biologically interesting difference between hepatocyte nuclei (snRNA pseudobulk) and isolated hepatocytes (bulk RNA-seq), or show that an accepted result is an artifact. Each hypothesis says what would count as evidence and which data it needs. Mark a task `[in progress: <who>, <date>]` or `[done: <date>, see <file>]` when you take or finish it.

## Hypotheses

**H1. Stromal carry-over makes "hepatocyte" disease signals in isolated preps.** Isolated-hepatocyte preparations carry a few percent of stellate, endothelial and other non-parenchymal cells, more so from cirrhotic livers. Their genes (collagens, ECM, angiogenesis, EMT) can then look like a hepatocyte response, which matters because hepatocyte EMT in human cirrhosis is still claimed while lineage tracing in mice argues against it. Evidence: the signal shrinks when the bulk model includes a non-parenchymal marker score, its leading edge is dominated by genes with low atlas hepatocyte fraction, and it is absent from ambient-corrected hepatocyte nuclei.

**H2. Ambient immune RNA inflates immune pathways in hepatocyte nuclei from inflamed livers.** In PSC/PBC explants, T/B-cell activation and inflammasome terms rise in hepatocyte-nucleus pseudobulks. Nuclei-suspension ambient RNA from infiltrating immune cells is the obvious suspect. Evidence: the terms shrink after ambient correction (empty-droplet profile per library, SoupX-style subtraction at pseudobulk level), and their leading-edge genes are immune-restricted (PTPRC, CD3E, CD247, LCK...).

**H3. Coagulation-factor output in disease.** Isolated-hepatocyte bulk and hepatocyte nuclei disagree on the direction of the Hallmark coagulation set in disease. Decompose the set into hepatocyte-made factors (FGA/FGB/FGG, F2, F9, F10, SERPINC1, PLG...) and endothelial ones (VWF, F8, THBD...). If the hepatocyte factors fall on both platforms and only endothelial carry-over rises in the bulk, the bulk "up" is an artifact and the nuclei show the real biology (falling synthesis, which matches the coagulopathy of end-stage disease).

**H4. Ribosome / translation.** Reversed between isolated bulk and both cells and nuclei of the same livers. Either an isolation response that differs by disease state, or etiology. Needs H5 to separate.

**H5. Etiology-matched nuclei.** The public sn/sc comparator is PSC/PBC (cholestatic). A MASH or alcohol-related cirrhosis snRNA dataset with healthy controls would match the most common etiologies of isolated-hepatocyte cohorts and is the single biggest trust upgrade. Candidate to verify first: Gribben et al. 2024 (MASLD spectrum snRNA-seq, GSE202379). Check that hepatocyte nuclei are labelled or labellable, which fibrosis stages exist and how many healthy controls there are.

**H6. Re-test isolated-hepatocyte findings in nuclei.** The repo owner has a local list of findings made in isolated-hepatocyte bulk. Check each in hepatocyte nuclei (local only, because the findings are unpublished).

## Tasks

Public data only (cloud sessions can do these):

- [ ] T1 (H5) Find and vet etiology-matched liver snRNA datasets: MASH and/or ALD cirrhosis plus healthy controls, with raw counts and ideally cell labels. Record accession, PMID (verified), n per group, chemistry and labels in `docs/datasets.md`.
- [ ] T2 (H5) For the best T1 dataset: hepatocyte-nucleus pseudobulk per donor using the same cell rule and reference labelling as `stratified/scripts/11*`. Then limma-voom disease vs healthy and Hallmark/KEGG/Reactome/GO BP GSEA with the same settings. Commit the tables.
- [ ] T3 (H2) Ambient correction for the Andrews libraries: per-library ambient profile from barcodes with 10–100 UMIs, a contamination fraction per library (from genes absent in hepatocytes), and corrected hepatocyte pseudobulks. Then re-run the PSC/PBC vs healthy hepatocyte GSEA and report which pathways move.
- [ ] T4 (H3) Coagulation decomposition on the public data: per-gene log fold change of hepatocyte-made vs endothelial coagulation genes in hepatocyte nuclei and cells, PSC/PBC vs healthy, and in T2.
- [ ] T5 Literature check, with every PMID verified: direct comparisons of isolated-hepatocyte bulk with hepatocyte nuclei or cells; non-parenchymal contamination of hepatocyte isolations; human hepatocyte EMT evidence; hepatocyte coagulation-factor transcription in cirrhosis.
- [ ] T6 Draft (do not send) the emails in `docs/open-questions-for-researchers.md` into `docs/outreach/`.

Local only (repo owner's machine, private cohort):

- [ ] L1 (H1) Refit the bulk disease model with a non-parenchymal marker covariate; list pathways that lose significance.
- [ ] L2 (H3) Coagulation split by cell of origin in the bulk.
- [ ] L3 (H6) Re-test the local list of findings in hepatocyte nuclei.
