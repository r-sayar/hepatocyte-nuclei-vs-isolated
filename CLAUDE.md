# Instructions for Claude sessions (local and cloud)

Read `README.md`, then `docs/research-plan.md`. Pick an open task there, mark it in progress in that file, and record what you did and found.

## Data rules

- **Never commit anything derived from the private isolated-hepatocyte bulk cohort** (numbers, tables, figures, gene lists, or prose results). It is unpublished. Cloud sessions do not have it; only the local machine does. Public-data results are fine to commit.
- Raw and processed data do not go into git (`.gitignore` covers `data/`). Commit code, small summary tables and figures from public data.
- New public datasets: record accession, paper, PMID and what was downloaded in the relevant README.

## Scientific rules

- **Never write a PMID, DOI or author email from memory.** Look up every PMID with NCBI E-utilities (`esummary.fcgi?db=pubmed&id=...`) and keep the title next to it. A citation must directly support the sentence it is attached to.
- Compare like with like. Hepatocytes vs hepatocytes; same limma-voom engine; same GSEA settings (gseapy prerank 1000 permutations, set size 10–500, seed 1234, or the validated `fast_gsea` in `*/scripts/ref/common.py`).
- Every "significant here, not there" claim needs the sample-size yardstick: the larger dataset cut to the smaller one's design and scored identically.
- Before believing a pathway hit, read its leading edge: which genes drive it, and from which cell type they normally come (MacParland atlas hepatocyte fraction).
- Do not call a sample "excluded" or "failed"; report the numbers and leave the decision to the repo owner.

## Contacting researchers

Cloud sessions may **draft** emails to original authors (questions in `docs/open-questions-for-researchers.md`). Do not send anything. Put drafts in `docs/outreach/` as markdown, with the recipient as listed on the paper's corresponding-author line (verify it on the publisher's page, never guess an address), a short subject, and one specific question per email. The repo owner reviews and sends.

## Budget

Keep runs small: prefer pseudobulk-level work over re-processing raw FASTQ, and avoid downloads above a few GB without a note in the research plan explaining why.

## Writing

Markdown prose is one line per paragraph (no hard wrapping); the notes are read in Obsidian.
