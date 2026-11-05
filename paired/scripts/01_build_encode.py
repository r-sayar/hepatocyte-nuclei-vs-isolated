#!/usr/bin/env python3
"""01_build_encode.py — ENCODE4 human liver: bulk total RNA-seq vs snRNA-seq pseudobulk, same donor and lobe.

Pairs (4 donors, 5 donor x lobe pairs; manifest in data/encode_manifest.tsv):
  bulk  Wold lab total RNA-seq, RSEM gene quantifications, expected_count (GRCh38, GENCODE V29)
  sn    Snyder lab 10x snRNA-seq, STARsolo GeneFull_Ex50pAS *filtered* matrix of unique reads
        (cell-called barcodes; exons + introns), summed over all nuclei = pseudobulk

Both are keyed on unversioned ENSG. Output (data/encode/):
  bulk_counts_ensg.csv, pseudobulk_counts_ensg.csv   genes x pair_id (integer)
  samples.tsv                                         pair_id, donor, lobe, n_nuclei, sn_umis, bulk_reads
"""
import io, os, tarfile
import numpy as np, pandas as pd
from scipy.io import mmread

P = "/srv/public/sayar99/liver_pb_gsea/data"
E = os.path.join(P, "encode")
man = pd.read_csv(os.path.join(P, "encode_manifest.tsv"), sep="\t")
man["pair_id"] = man["donor"] + "_" + man["lobe"]

bulk, pb, meta = {}, {}, {}
for _, r in man.iterrows():
    if r.modality == "bulk":
        q = pd.read_csv(os.path.join(E, f"{r.file}.tsv"), sep="\t", usecols=["gene_id", "expected_count"])
        q = q[q.gene_id.str.startswith("ENSG")]
        q["ensg"] = q.gene_id.str.split(".").str[0]
        # PAR_Y duplicates collapse onto the same ENSG; sum them
        bulk[r.pair_id] = q.groupby("ensg")["expected_count"].sum()
        meta.setdefault(r.pair_id, {})["bulk_reads"] = float(bulk[r.pair_id].sum())
    else:
        with tarfile.open(os.path.join(E, f"{r.file}.tar.gz")) as t:
            names = {os.path.basename(m.name): m for m in t.getmembers() if "/filtered/" in m.name}
            feat = pd.read_csv(t.extractfile(names["features.tsv"]), sep="\t", header=None)
            bc = pd.read_csv(t.extractfile(names["barcodes.tsv"]), sep="\t", header=None)
            M = mmread(io.BytesIO(t.extractfile(names["matrix.mtx"]).read())).tocsr()
        assert M.shape == (len(feat), len(bc)), (r.file, M.shape, len(feat), len(bc))
        s = pd.Series(np.asarray(M.sum(axis=1)).ravel(), index=feat[0].str.split(".").str[0])
        pb[r.pair_id] = s.groupby(level=0).sum()
        mt = feat[1].str.upper().str.startswith("MT-").values
        meta.setdefault(r.pair_id, {}).update(
            donor=r.donor, lobe=r.lobe, n_nuclei=int(M.shape[1]), sn_umis=float(s.sum()),
            sn_mt_frac=float(s[mt].sum() / s.sum()),
            sn_median_umi=float(np.median(np.asarray(M.sum(axis=0)).ravel())))
        print(r.pair_id, "nuclei", M.shape[1], "UMIs", int(s.sum()), flush=True)

B = pd.DataFrame(bulk).fillna(0)
S = pd.DataFrame(pb).fillna(0)
cols = sorted(B.columns)
assert cols == sorted(S.columns)
B, S = B[cols].round().astype(int), S[cols].round().astype(int)
B.to_csv(os.path.join(E, "bulk_counts_ensg.csv"))
S.to_csv(os.path.join(E, "pseudobulk_counts_ensg.csv"))
m = pd.DataFrame(meta).T.loc[cols]
m.index.name = "pair_id"
m.to_csv(os.path.join(E, "samples.tsv"), sep="\t")
print(m.to_string())
print("bulk genes", B.shape, "sn genes", S.shape, "shared ENSG", len(B.index.intersection(S.index)))
