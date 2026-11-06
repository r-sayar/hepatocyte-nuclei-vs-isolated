#!/usr/bin/env python3
"""11c_label_pseudobulk.py — step 2 of 2: label every cluster against the authors' own healthy-liver cell types
and build hepatocyte-only pseudobulks, the material that matches the isolated-hepatocyte bulk cohort.

Labelling is reference-based at cluster level (the SingleR idea, Aran et al. 2019): for each modality (sc, sn)
the reference is the mean log-expression of each author cell type (Subcluster_Group) in the healthy libraries
of that modality; a cluster takes the type whose reference correlates best (Spearman) with the cluster's mean
log-expression over a marker panel (per type: top 50 genes by margin over the next-highest type, reference
mean > 0.3). Healthy libraries are labelled with references built WITHOUT their own liver (leave-one-liver-out)
and scored against the authors' per-cell labels; PSC/PBC libraries use all four healthy livers.
Output (data/): pb_hep_ensg.csv, libraries.tsv, label_check.tsv, cluster_labels.tsv
"""
import os
import numpy as np, pandas as pd

D = "/srv/public/sayar99/liver_strat_gsea/data"
LIBDIR = os.path.join(D, "lib")
MIN_REF_CELLS, TOP, MIN_MEAN = 10, 50, 0.3
design = pd.read_csv(os.path.join(D, "libraries_design.tsv"), sep="\t")
shared = pd.Index(open(os.path.join(D, "shared_ensg.txt")).read().split())


def load(name):
    z = np.load(os.path.join(LIBDIR, f"{name}.npz"), allow_pickle=True)
    ensg = pd.Index(z["ensg"])
    first = ~ensg.duplicated()
    pos = pd.Series(np.flatnonzero(first), index=ensg[first]).reindex(shared).to_numpy()
    on_shared = lambda m: m[:, pos]  # every shared ENSG is present in every library by construction (11b)
    out = dict(z=z, logmean=on_shared(z["cluster_logmean"]))
    if "author_types" in z.files:
        n = z["author_n"].astype(float)
        out["ref"] = {t: (on_shared(z["author_logsum"])[i] / n[i], n[i]) for i, t in enumerate(z["author_types"])}
    return out


L = {n: load(n) for n in design.library}


def reference(modality, exclude_liver=None):
    """per-type centroid = cell-weighted mean of the healthy libraries' type means; plus the marker panel."""
    libs = design[(design.modality == modality) & (design.group == "healthy") & (design.liver != exclude_liver)].library
    acc = {}
    for n in libs:
        for t, (mu, k) in L[n]["ref"].items():
            if k >= MIN_REF_CELLS:
                s, w = acc.get(t, (0, 0))
                acc[t] = (s + mu * k, w + k)
    cent = pd.DataFrame({t: s / w for t, (s, w) in acc.items()}, index=shared)
    panel = set()
    for t in cent.columns:
        margin = cent[t] - cent.drop(columns=t).max(axis=1)
        panel |= set(margin[cent[t] > MIN_MEAN].nlargest(TOP).index)
    return cent.loc[sorted(panel)]


def label_clusters(name, cent):
    X = pd.DataFrame(L[name]["logmean"].T, index=shared).loc[cent.index]
    rho = pd.concat([X, cent], axis=1).rank().corr().loc[X.columns, cent.columns]  # Spearman
    return rho.idxmax(axis=1).to_numpy(), rho.max(axis=1).to_numpy()


info, checks, rows, hep = [], [], [], {}
for _, d in design.iterrows():
    z = L[d.library]["z"]
    cent = reference(d.modality, exclude_liver=d.liver if d.group == "healthy" else None)
    lab, r = label_clusters(d.library, cent)
    ncl = z["cluster_n"]
    is_hep = lab == "Hepatocyte"
    cc = pd.DataFrame(z["cluster_counts"][is_hep].sum(axis=0) if is_hep.any() else np.zeros(len(z["ensg"])),
                      index=z["ensg"])[0]
    hep[d.library] = cc.groupby(level=0).sum()
    comp = pd.Series(ncl).groupby(lab).sum()
    info.append(dict(library=d.library, liver=d.liver, group=d.group, modality=d.modality, n_cells=int(ncl.sum()),
                     n_clusters=len(ncl), n_hep=int(ncl[is_hep].sum()), hep_frac=float(ncl[is_hep].sum() / ncl.sum()),
                     hep_umis=float(cc.sum()), **{f"n_{t}": int(comp.get(t, 0)) for t in cent.columns}))
    rows.append(pd.DataFrame({"library": d.library, "cluster": np.arange(len(ncl)), "n": ncl, "label": lab, "spearman": r.round(3)}))
    if "author_labels" in z.files:
        cell_lab = lab[z["cluster"]]
        au = z["author_labels"]
        has = au != "nan"
        ours, theirs = cell_lab[has] == "Hepatocyte", au[has] == "Hepatocyte"
        checks.append(dict(library=d.library, n_author_cells_matched=int(has.sum()),
                           agree_all_types=float((cell_lab[has] == au[has]).mean()),
                           hep_recall=float((ours & theirs).sum() / max(theirs.sum(), 1)),
                           hep_precision=float((ours & theirs).sum() / max(ours.sum(), 1))))
    print(d.library, info[-1]["n_hep"], "hepatocytes of", info[-1]["n_cells"], flush=True)

pd.DataFrame(hep).fillna(0).round().astype(int).to_csv(os.path.join(D, "pb_hep_ensg.csv"))
I = pd.DataFrame(info)
I.to_csv(os.path.join(D, "libraries.tsv"), sep="\t", index=False)
K = pd.DataFrame(checks)
K.to_csv(os.path.join(D, "label_check.tsv"), sep="\t", index=False)
pd.concat(rows).to_csv(os.path.join(D, "cluster_labels.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250)
print(I.round(3).to_string(index=False))
print(K.round(3).to_string(index=False))
