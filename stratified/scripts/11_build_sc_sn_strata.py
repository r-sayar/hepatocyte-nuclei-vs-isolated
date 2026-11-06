#!/usr/bin/env python3
"""11_build_sc_sn_strata.py — step 1 of 2: cluster every scRNA / snRNA library of the same livers (Andrews /
MacParland lab) and save what the labelling step needs.

Libraries: healthy donor livers (GSE185477; sc vs TST nuclei, 4 livers) and PSC/PBC explants (GSE247128;
sc 5' vs sn 3', 6 livers). Raw matrices are read from the paired analysis' download directory (read only).
Cells: UMIs >= 500, genes >= 200, MT < 0.50 (same rule as the paired analysis).
Per library: scanpy log-normalise (1e4) -> 2000 HVG -> PCA 30 -> kNN 15 -> Leiden 1.0, then save
  - per-cluster summed counts (pseudobulk building blocks) and mean log-expression
  - for the healthy libraries, mean log-expression of the authors' own cell types
    (GSE185477_Final_Metadata, Subcluster_Group), which become the labelling reference in step 2
Output: data/lib/<library>.npz, data/pb_whole_ensg.csv, data/genes.tsv
"""
import gzip, io, os, sys, zipfile
from multiprocessing import Pool
import numpy as np, pandas as pd
import scanpy as sc
import anndata as ad
from scipy.io import mmread
from scipy import sparse

SRC = "/srv/public/sayar99/liver_pb_gsea/data/andrews"   # downloads of the paired analysis (read only)
OUT = "/srv/public/sayar99/liver_strat_gsea/data"
LIBDIR = os.path.join(OUT, "lib")
MIN_UMI, MIN_GENES, MAX_MT = 500, 200, 0.50
LIBS = []  # (library, liver, group, modality, source, author sample id)
for d, cid, sc_zip, sn_zip in [
        ("Healthy1", "C41", "GSE185477_GSM3178784_C41_SC_raw_counts.zip", "GSM5616001_C41_TST_raw_counts.zip"),
        ("Healthy2", "C58", "GSM5616002_C58_SC_raw_counts.zip", "GSM5616003_C58_TST_raw_counts.zip"),
        ("Healthy3", "C70", "GSM5616004_C70_SC_raw_counts.zip", "GSM5616005_C70_TST_raw_counts.zip"),
        ("Healthy4", "C72", "GSM5616006_C72_SC_raw_counts.zip", "GSM5616007_C72_TST_raw_counts.zip")]:
    LIBS.append((f"{d}_sc", d, "healthy", "sc", ("zip", sc_zip), cid if cid == "C41" else f"{cid}_RESEQ"))
    LIBS.append((f"{d}_sn", d, "healthy", "sn", ("zip", sn_zip), f"{cid}_TST"))
for liv, (gsc, gsn) in {"PBC001": ("GSM7883930", "GSM7883929"), "PSC005": ("GSM7883933", "GSM7883934"),
                        "PSC012": ("GSM7883937", "GSM7883938"), "PSC014": ("GSM7883940", "GSM7883939"),
                        "PSC018": ("GSM7883942", "GSM7883943"), "PSC019": ("GSM7883944", "GSM7883945")}.items():
    LIBS.append((f"{liv}_sc", liv, "PSC_PBC", "sc", ("triplet", f"{gsc}_{liv}_SC_5pr"), None))
    LIBS.append((f"{liv}_sn", liv, "PSC_PBC", "sn", ("triplet", f"{gsn}_{liv}_SN_3pr"), None))


def load(src):
    """-> barcodes (no -1 suffix), features DataFrame, genes x barcodes CSC."""
    if src[0] == "triplet":
        p = os.path.join(SRC, src[1])
        fb, ff, fm = p + "_barcodes.tsv.gz", p + "_features.tsv.gz", p + "_matrix.mtx.gz"
    else:
        z = zipfile.ZipFile(os.path.join(SRC, src[1]))
        names = z.namelist()
        pick = lambda key: next(n for n in names if os.path.basename(n).startswith(key) and n.endswith(".gz"))
        rd = lambda n: io.BytesIO(gzip.decompress(z.read(n)))
        fb, fm = rd(pick("barcodes")), rd(pick("matrix"))
        ff = rd(next(n for n in names if os.path.basename(n).split(".")[0] in ("features", "genes") and n.endswith(".gz")))
    bc = pd.read_csv(fb, sep="\t", header=None)[0].astype(str).str.split("-").str[0].to_numpy(dtype=object)
    ft = pd.read_csv(ff, sep="\t", header=None)
    M = sparse.csc_matrix(mmread(fm))
    assert M.shape == (len(ft), len(bc))
    return bc, ft, M


def group_sum(X, codes, k):
    """rows of X (cells x genes, CSR) summed per group code -> k x genes dense."""
    G = sparse.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(k, len(codes)))
    return np.asarray((G @ X).todense())


def build(lib):
    name, liver, group, modality, src, author = lib
    bc, ft, M = load(src)
    ensg = ft[0].astype(str).str.split(".").str[0].to_numpy(dtype=object)
    sym = ft[1].astype(str).to_numpy(dtype=object)
    mt = pd.Series(sym).str.upper().str.startswith("MT-").to_numpy(dtype=bool)
    umi = np.asarray(M.sum(axis=0)).ravel()
    keep = (umi >= MIN_UMI) & (np.diff(M.indptr) >= MIN_GENES) & (np.asarray(M[mt].sum(axis=0)).ravel() / np.maximum(umi, 1) < MAX_MT)
    raw = M[:, keep].T.tocsr().astype(np.float32)
    a = ad.AnnData(X=raw.copy(), obs=pd.DataFrame(index=bc[keep]))
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    logx = a.X.copy()
    sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor="seurat")
    b = a[:, a.var.highly_variable].copy()
    sc.pp.scale(b, max_value=10)
    sc.tl.pca(b, n_comps=30, random_state=0)
    sc.pp.neighbors(b, n_neighbors=15, random_state=0)
    sc.tl.leiden(b, resolution=1.0, flavor="igraph", n_iterations=2, random_state=0)
    cl = b.obs["leiden"].astype(int).to_numpy()
    k = int(cl.max()) + 1
    ncl = np.bincount(cl, minlength=k)
    out = dict(ensg=ensg, symbol=sym, cluster=cl, cluster_n=ncl, barcodes=a.obs_names.to_numpy(dtype=str),
               cluster_counts=group_sum(raw, cl, k), cluster_logmean=group_sum(logx, cl, k) / ncl[:, None])
    if author is not None:
        meta = pd.read_csv(os.path.join(SRC, "GSE185477_Final_Metadata.txt.gz"), sep="\t",
                           usecols=["sample", "cell_barcode", "Subcluster_Group"])
        au = meta[meta["sample"] == author].drop_duplicates("cell_barcode").set_index("cell_barcode")["Subcluster_Group"]
        lab = pd.Series(a.obs_names).map(au).to_numpy(dtype=object)
        types = sorted(pd.unique(au.values))
        has = pd.notna(lab)
        codes = np.array([types.index(x) for x in lab[has]])
        out.update(author_labels=lab.astype(str), author_types=np.array(types),
                   author_n=np.bincount(codes, minlength=len(types)),
                   author_logsum=group_sum(logx[has], codes, len(types)))
    np.savez_compressed(os.path.join(LIBDIR, f"{name}.npz"), **out)
    whole = pd.Series(np.asarray(raw.sum(axis=0)).ravel(), index=ensg).groupby(level=0).sum()
    print(name, "cells", int(keep.sum()), "clusters", k, flush=True)
    return name, whole


if __name__ == "__main__":
    os.makedirs(LIBDIR, exist_ok=True)
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 4) as pool:
        res = pool.map(build, LIBS)
    pd.DataFrame(dict(res)).fillna(0).round().astype(int).to_csv(os.path.join(OUT, "pb_whole_ensg.csv"))
    pd.DataFrame([dict(library=l[0], liver=l[1], group=l[2], modality=l[3]) for l in LIBS]).to_csv(
        os.path.join(OUT, "libraries_design.tsv"), sep="\t", index=False)
