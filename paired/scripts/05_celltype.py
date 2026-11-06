#!/usr/bin/env python3
"""05_celltype.py — is the sc-vs-sn inflammation shift cell composition, or a per-cell effect?

Healthy livers only (GSE185477), where the authors' QC'd cells carry cell-type labels on BOTH the
sc and the sn (TST) library of each liver (GSE185477_Final_Metadata.txt.gz, Subcluster_Group).
  per-type   pseudobulk each cell type separately; sc vs sn within the type (~liver + modality,
             protein-coding). Types kept only if every one of the 8 libraries has >= MIN_CELLS.
  matched    composition-matched whole-liver pseudobulk: every library re-weighted to the same
             cell-type mix (the per-liver sn proportions), so only per-cell differences remain
  markers    CPM of composition markers in every pseudobulk / bulk (Andrews 20 libs + ENCODE)
Outputs: results/CT/..., results/ct_*.csv
"""
import gzip, io, os, zipfile
import numpy as np, pandas as pd
from scipy.io import mmread
from scipy import sparse
import liver_common as L
from liver_common import C

OUT = os.path.join(L.RES, "CT")
os.makedirs(OUT, exist_ok=True)
MIN_CELLS = 30
LIBS = {("Healthy1", "sc"): ("GSE185477_GSM3178784_C41_SC_raw_counts.zip", "C41"),
        ("Healthy1", "sn"): ("GSM5616001_C41_TST_raw_counts.zip", "C41_TST"),
        ("Healthy2", "sc"): ("GSM5616002_C58_SC_raw_counts.zip", "C58_RESEQ"),
        ("Healthy2", "sn"): ("GSM5616003_C58_TST_raw_counts.zip", "C58_TST"),
        ("Healthy3", "sc"): ("GSM5616004_C70_SC_raw_counts.zip", "C70_RESEQ"),
        ("Healthy3", "sn"): ("GSM5616005_C70_TST_raw_counts.zip", "C70_TST"),
        ("Healthy4", "sc"): ("GSM5616006_C72_SC_raw_counts.zip", "C72_RESEQ"),
        ("Healthy4", "sn"): ("GSM5616007_C72_TST_raw_counts.zip", "C72_TST")}

pb_file = os.path.join(OUT, "celltype_pseudobulk_ensg.csv")
if not os.path.exists(pb_file):
    meta = pd.read_csv(os.path.join(L.A, "GSE185477_Final_Metadata.txt.gz"), sep="\t",
                       usecols=["sample", "cell_barcode", "Subcluster_Group"])
    cols, ncell = {}, []
    for (liver, mod), (zname, samp) in LIBS.items():
        with zipfile.ZipFile(os.path.join(L.A, zname)) as z:
            names = z.namelist()
            pick = lambda key: next(n for n in names if os.path.basename(n).startswith(key) and n.endswith(".gz"))
            rd = lambda n: io.BytesIO(gzip.decompress(z.read(n)))
            bc = pd.read_csv(rd(pick("barcodes")), sep="\t", header=None)[0].astype(str).str.split("-").str[0].to_numpy(dtype=object)
            ftn = next(n for n in names if os.path.basename(n).split(".")[0] in ("features", "genes") and n.endswith(".gz"))
            ensg = pd.read_csv(rd(ftn), sep="\t", header=None)[0].astype(str).str.split(".").str[0].to_numpy(dtype=object)
            M = sparse.csc_matrix(mmread(rd(pick("matrix"))))
        lab = meta[meta["sample"] == samp].set_index("cell_barcode")["Subcluster_Group"]
        lab = lab[~lab.index.duplicated()]
        ct = pd.Series(bc).map(lab).to_numpy(dtype=object)
        for t in pd.unique(lab.values):
            sel = ct == t
            cols[f"{liver}_{mod}__{t}"] = pd.Series(np.asarray(M[:, sel].sum(axis=1)).ravel(), index=ensg).groupby(level=0).sum()
            ncell.append(dict(liver=liver, modality=mod, celltype=t, n_cells=int(sel.sum())))
        print(liver, mod, "labelled cells matched", int(pd.notna(ct).sum()), "of", len(lab), flush=True)
    pd.DataFrame(cols).fillna(0).round().astype(int).to_csv(pb_file)
    pd.DataFrame(ncell).to_csv(os.path.join(OUT, "celltype_ncells.csv"), index=False)

PB = pd.read_csv(pb_file, index_col=0)
PB = PB.reindex(sorted(L.andrews_shared_ensg())).fillna(0).astype(int)
NC = pd.read_csv(os.path.join(OUT, "celltype_ncells.csv"))
wide = NC.pivot_table(index="celltype", columns=["liver", "modality"], values="n_cells").fillna(0)
types = [t for t in wide.index if (wide.loc[t] >= MIN_CELLS).all()]
print("cell types with >=", MIN_CELLS, "cells in all 8 libraries:", types, flush=True)
livers = sorted(NC.liver.unique())

INFL, SUMM, NODIS = [], [], []
for t in types:
    cols = [f"{l}_{m}__{t}" for l in livers for m in ("sc", "sn")]
    X = PB[cols].rename(columns=lambda c: c.split("__")[0])
    md = L.paired_md(list(X.columns), "liver", "modality")
    res, tabs = L.run(L.protein_coding(X), md, "~liver + modality", ("modality", "sc", "sn"),
                      f"{OUT}/sc_vs_sn_{t}", libs=["hallmark", "kegg", "reactome", "gobp"] if t == "Hepatocyte" else ["hallmark"])
    lab = f"P5 healthy sc vs sn within {t}"
    INFL.append(L.infl_table(tabs, lab))
    SUMM.append(L.infl_summary(tabs, lab))
    NODIS.append(L.without_dissoc(res, lab))

# composition-matched whole-liver pseudobulk: sc re-weighted to the sn cell-type mix of the same liver
allt = sorted(NC.celltype.unique())
matched = {}
for l in livers:
    frac_sn = wide[(l, "sn")] / wide[(l, "sn")].sum()
    for m in ("sc", "sn"):
        prof = {}
        for t in allt:
            c = f"{l}_{m}__{t}"
            n = wide.at[t, (l, m)]
            if c in PB and n >= 1:
                prof[t] = PB[c] / n  # mean UMIs per cell of type t
        tot_cells = wide[(l, m)].sum()
        s = sum(prof[t] * frac_sn.get(t, 0) * tot_cells for t in prof)
        matched[f"{l}_{m}"] = s
MT = pd.DataFrame(matched).round().astype(int)
md = L.paired_md(list(MT.columns), "liver", "modality")
res_m, tabs_m = L.run(L.protein_coding(MT), md, "~liver + modality", ("modality", "sc", "sn"), f"{OUT}/sc_vs_sn_composition_matched")
INFL.append(L.infl_table(tabs_m, "P5 healthy sc vs sn, composition-matched"))
SUMM.append(L.infl_summary(tabs_m, "P5 healthy sc vs sn, composition-matched"))
NODIS.append(L.without_dissoc(res_m, "P5 healthy sc vs sn, composition-matched"))
# same 4 livers, same author cells, NOT matched (reference for the matched run)
UN = pd.DataFrame({f"{l}_{m}": PB[[c for c in PB.columns if c.startswith(f"{l}_{m}__")]].sum(axis=1) for l in livers for m in ("sc", "sn")})
res_u, tabs_u = L.run(L.protein_coding(UN), md, "~liver + modality", ("modality", "sc", "sn"), f"{OUT}/sc_vs_sn_author_unmatched")
INFL.append(L.infl_table(tabs_u, "P5 healthy sc vs sn, author cells unmatched"))
SUMM.append(L.infl_summary(tabs_u, "P5 healthy sc vs sn, author cells unmatched"))

# composition markers, CPM
MARK = ["ALB", "APOA1", "CYP3A4", "KRT19", "EPCAM", "PECAM1", "STAB2", "DCN", "COL1A1", "PTPRC", "CD74",
        "HLA-DRA", "CD68", "MARCO", "CD3E", "NKG7", "MS4A1", "HBB", "FOS", "JUN", "EGR1", "HSPA1A"]
inv = {s: e for e, s in L.SYM.items() if s in MARK}
frames = {}
RB, _ = L.load_andrews("rule")
for c in RB.columns:
    frames[f"Andrews {c}"] = RB[c]
Bk, Sn = L.load_encode()
for c in Bk.columns:
    frames[f"ENCODE {c} bulk"] = Bk[c]
    frames[f"ENCODE {c} sn"] = Sn[c]
cpm = pd.DataFrame({k: (v / v.sum() * 1e6).reindex([inv[s] for s in MARK if s in inv]) for k, v in frames.items()})
cpm.index = [L.SYM[e] for e in cpm.index]
cpm.T.round(1).to_csv(f"{L.RES}/ct_marker_cpm.csv")

pd.concat(INFL).to_csv(f"{L.RES}/ct_inflammation_sets.csv", index=False)
pd.concat(SUMM).to_csv(f"{L.RES}/ct_inflammation_summary.csv", index=False)
pd.concat(NODIS).to_csv(f"{L.RES}/ct_inflammation_without_dissoc_genes.csv", index=False)
pd.set_option("display.width", 250)
H = pd.concat(INFL)
print(H[H.lib == "hallmark"].pivot_table(index="term", columns="analysis", values="NES").round(2).T.to_string())
print(pd.concat(NODIS).round(3).to_string(index=False))
print(cpm.T.round(0).to_string())
print(wide.astype(int).to_string())
