#!/usr/bin/env python3
"""00_gene_lengths.py — per-gene union exonic length, genomic span and intron fraction (GENCODE v46)."""
import re, sys
import numpy as np, pandas as pd
GTF = "/storage/mi/sayar99/scrna_realign/genome/gencode.v46.annotation.gtf"
span, exons = {}, {}
rx = re.compile(r'gene_id "(ENSG\d+)')
with open(GTF) as fh:
    for line in fh:
        if line[0] == "#":
            continue
        f = line.split("\t", 9)
        if f[2] not in ("gene", "exon"):
            continue
        g = rx.search(f[8]).group(1)
        s, e = int(f[3]), int(f[4])
        if f[2] == "gene":
            span[g] = max(span.get(g, 0), e - s + 1)
        else:
            exons.setdefault(g, []).append((s, e))
rows = []
for g, iv in exons.items():
    iv.sort()
    tot, cs, ce = 0, iv[0][0], iv[0][1]
    for s, e in iv[1:]:
        if s > ce + 1:
            tot += ce - cs + 1; cs, ce = s, e
        else:
            ce = max(ce, e)
    tot += ce - cs + 1
    rows.append((g, tot, span.get(g, tot)))
d = pd.DataFrame(rows, columns=["ensg", "exonic_len", "span"]).set_index("ensg")
d["intron_frac"] = 1 - d.exonic_len / d.span.clip(lower=1)
d.to_csv("/srv/public/sayar99/liver_pb_gsea/data/gene_lengths_v46.tsv", sep="\t")
print(d.describe().round(2).to_string())
