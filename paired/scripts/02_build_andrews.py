#!/usr/bin/env python3
"""02_build_andrews.py — same-liver scRNA-seq vs snRNA-seq pseudobulk (MacParland lab, Toronto).

Healthy donor livers, GSE185477 (Andrews 2022, Hepatol Commun): sc vs sn (TST nuclei protocol)
  Donor 1 = C41 (sc is MacParland 2018 GSM3178784), Donor 2 = C58, Donor 3 = C70, Donor 4 = C72
PSC / PBC explants, GSE247128 (Andrews 2024, J Hepatol): sc 5' vs sn 3' of the same explant
  PBC001, PSC005, PSC012, PSC014, PSC018, PSC019

GEO holds raw (unfiltered) barcode matrices, so one cell-calling rule is applied to every library:
  UMIs >= 500, detected genes >= 200, MT fraction < 0.50 (the authors' own MT ceiling).
For the four healthy pairs the authors' published QC'd barcodes (GSE185477_Final_Metadata) give a
second pseudobulk ("author" cells), used to check that the rule does not drive the result.

Output (data/andrews/):
  pseudobulk_rule_ensg.csv, pseudobulk_author_ensg.csv   genes x library (integer)
  libraries.tsv                                           library, liver, group, modality, n_cells, ...
"""
import gzip, io, os, sys, zipfile
from multiprocessing import Pool
import numpy as np, pandas as pd
from scipy.io import mmread
from scipy import sparse

A = "/srv/public/sayar99/liver_pb_gsea/data/andrews"
MIN_UMI, MIN_GENES, MAX_MT = 500, 200, 0.50

LIBS = []  # (library, liver, group, modality, source)
for d, cid, sc_zip, sn_zip in [
        ("Healthy1", "C41", "GSE185477_GSM3178784_C41_SC_raw_counts.zip", "GSM5616001_C41_TST_raw_counts.zip"),
        ("Healthy2", "C58", "GSM5616002_C58_SC_raw_counts.zip", "GSM5616003_C58_TST_raw_counts.zip"),
        ("Healthy3", "C70", "GSM5616004_C70_SC_raw_counts.zip", "GSM5616005_C70_TST_raw_counts.zip"),
        ("Healthy4", "C72", "GSM5616006_C72_SC_raw_counts.zip", "GSM5616007_C72_TST_raw_counts.zip")]:
    LIBS.append((f"{d}_sc", d, "healthy", "sc", ("zip", sc_zip, cid, "sc")))
    LIBS.append((f"{d}_sn", d, "healthy", "sn", ("zip", sn_zip, cid, "sn")))
GSM24 = {"PBC001": ("GSM7883930", "GSM7883929"), "PSC005": ("GSM7883933", "GSM7883934"),
         "PSC012": ("GSM7883937", "GSM7883938"), "PSC014": ("GSM7883940", "GSM7883939"),
         "PSC018": ("GSM7883942", "GSM7883943"), "PSC019": ("GSM7883944", "GSM7883945")}
for liv, (gsc, gsn) in GSM24.items():
    grp = liv[:3]
    LIBS.append((f"{liv}_sc", liv, grp, "sc", ("triplet", f"{gsc}_{liv}_SC_5pr")))
    LIBS.append((f"{liv}_sn", liv, grp, "sn", ("triplet", f"{gsn}_{liv}_SN_3pr")))


def _read(fh_bc, fh_ft, fh_mx):
    bc = pd.read_csv(fh_bc, sep="\t", header=None)[0].astype(str).to_numpy(dtype=object)
    ft = pd.read_csv(fh_ft, sep="\t", header=None)
    M = sparse.csc_matrix(mmread(fh_mx))  # genes x barcodes
    assert M.shape == (len(ft), len(bc)), (M.shape, len(ft), len(bc))
    return bc, ft, M


def load(src):
    if src[0] == "triplet":
        p = os.path.join(A, src[1])
        return _read(p + "_barcodes.tsv.gz", p + "_features.tsv.gz", p + "_matrix.mtx.gz")
    _, zname, _, _ = src
    with zipfile.ZipFile(os.path.join(A, zname)) as z:
        names = z.namelist()
        pick = lambda key: next(n for n in names if os.path.basename(n).startswith(key) and n.endswith(".gz"))
        bcn, mxn = pick("barcodes"), pick("matrix")
        ftn = next(n for n in names if os.path.basename(n).split(".")[0] in ("features", "genes") and n.endswith(".gz"))
        rd = lambda n: io.BytesIO(gzip.decompress(z.read(n)))
        return _read(rd(bcn), rd(ftn), rd(mxn))


def build(lib):
    name, liver, group, modality, src = lib
    bc, ft, M = load(src)
    bc = np.array([b.split("-")[0] for b in bc])
    ensg = ft[0].astype(str).str.split(".").str[0].to_numpy(dtype=object)
    mt = ft[1].astype(str).str.upper().str.startswith("MT-").to_numpy(dtype=bool)
    umi = np.asarray(M.sum(axis=0)).ravel()
    ngene = np.diff(M.indptr)
    mtf = np.asarray(M[mt].sum(axis=0)).ravel() / np.maximum(umi, 1)
    keep = (umi >= MIN_UMI) & (ngene >= MIN_GENES) & (mtf < MAX_MT)
    rule = pd.Series(np.asarray(M[:, keep].sum(axis=1)).ravel(), index=ensg).groupby(level=0).sum()
    info = dict(library=name, liver=liver, group=group, modality=modality, n_barcodes=len(bc),
                n_cells_rule=int(keep.sum()), umis_rule=float(rule.sum()),
                median_umi_rule=float(np.median(umi[keep])) if keep.any() else 0.0,
                mt_frac_rule=float(np.asarray(M[mt][:, keep].sum()).item() / max(rule.sum(), 1)))
    author = None
    if src[0] == "zip":
        cid, mod = src[2], src[3]
        meta = pd.read_csv(os.path.join(A, "GSE185477_Final_Metadata.txt.gz"), sep="\t",
                           usecols=["sample", "cell_barcode"])
        samp = cid if mod == "sc" and cid == "C41" else (f"{cid}_RESEQ" if mod == "sc" else f"{cid}_TST")
        want = set(meta.loc[meta["sample"] == samp, "cell_barcode"].astype(str))
        sel = np.isin(bc, list(want))
        author = pd.Series(np.asarray(M[:, sel].sum(axis=1)).ravel(), index=ensg).groupby(level=0).sum()
        info.update(author_sample=samp, n_cells_author=int(sel.sum()), n_author_listed=len(want),
                    overlap_rule_author=int((sel & keep).sum()))
    print(name, {k: v for k, v in info.items() if k.startswith("n_")}, flush=True)
    return name, rule, author, info


if __name__ == "__main__":
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 8) as pool:
        res = pool.map(build, LIBS)
    R = pd.DataFrame({n: r for n, r, _, _ in res}).fillna(0).round().astype(int)
    Au = pd.DataFrame({n: a for n, _, a, _ in res if a is not None}).fillna(0).round().astype(int)
    R.to_csv(os.path.join(A, "pseudobulk_rule_ensg.csv"))
    Au.to_csv(os.path.join(A, "pseudobulk_author_ensg.csv"))
    info = pd.DataFrame([i for *_, i in res]).set_index("library")
    info.to_csv(os.path.join(A, "libraries.tsv"), sep="\t")
    pd.set_option("display.width", 250)
    print(info.to_string())
