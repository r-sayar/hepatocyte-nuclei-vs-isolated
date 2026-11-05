#!/usr/bin/env python3
"""11b_shared_genes.py — ENSG present in EVERY Andrews library's annotation (the 2018 C41 sc library uses an
older Ensembl release; a gene missing from one annotation would otherwise enter the pseudobulk as a real zero).
Also copies the GENCODE v46 ENSG -> symbol/type table used for protein-coding filtering and GSEA symbols.
Output: data/shared_ensg.txt, data/ensg2symbol_v46.tsv
"""
import gzip, io, os, shutil, zipfile
import pandas as pd

SRC = "/srv/public/sayar99/liver_pb_gsea/data"
OUT = "/srv/public/sayar99/liver_strat_gsea/data"
A = os.path.join(SRC, "andrews")
sets = []
for f in sorted(os.listdir(A)):
    if f.endswith("_features.tsv.gz"):
        sets.append(set(pd.read_csv(os.path.join(A, f), sep="\t", header=None)[0].str.split(".").str[0]))
    elif f.endswith("_raw_counts.zip"):
        with zipfile.ZipFile(os.path.join(A, f)) as z:
            n = next(n for n in z.namelist() if os.path.basename(n).split(".")[0] in ("features", "genes") and n.endswith(".gz"))
            sets.append(set(pd.read_csv(io.BytesIO(gzip.decompress(z.read(n))), sep="\t", header=None)[0].str.split(".").str[0]))
shared = sorted(set.intersection(*sets))
with open(os.path.join(OUT, "shared_ensg.txt"), "w") as fh:
    fh.write("\n".join(shared) + "\n")
shutil.copy(os.path.join(SRC, "ensg2symbol_v46.tsv"), os.path.join(OUT, "ensg2symbol_v46.tsv"))
print(f"{len(sets)} library annotations, {len(shared)} shared ENSG")
