"""Shared helpers: Pittsburgh-convention pathway enrichment on paired bulk / pseudobulk data.

Question: when we run pathway enrichment the way the lab's standard pipeline does, are the
top-15 / bottom-15 pathways actually different between bulk and pseudobulk of the same samples?

Pittsburgh convention, the lab's standard pathway pipeline:
  DE      DESeq2, rowSums(counts) >= 10 prefilter, donor in the design when paired
          
  rank    Wald stat; MT-* genes dropped; symbol-deduplicated keeping the highest baseMean
          (scripts/143_gsea_full_sweep.py)
  GSEA    gseapy prerank, 1000 permutations, set size 10-500, seed 1234
  report  sets with nominal p < 0.05, top 15 each way by NES
          (figures/gsea_gobp/*_p05.png, figures/deck_2026_10_03/gsea_kegg/*_p05.png)
Gene sets are the exact GMTs cached by the Pittsburgh pipeline (MSigDB 2026.1.Hs).
"""
from __future__ import annotations

import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

import gseapy as gp  # noqa: E402
from pydeseq2.dds import DeseqDataSet  # noqa: E402
from pydeseq2.default_inference import DefaultInference  # noqa: E402
from pydeseq2.ds import DeseqStats  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PL = os.path.join(ROOT, "data", "paired_local")
GMT_DIR = os.environ.get("GMT_DIR", "genesets")  # MSigDB 2026.1.Hs GMTs; set per machine
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

LIB_FILES = {
    "hallmark": "h.all.2026.1.Hs.gmt",
    "kegg": "c2.cp.kegg_legacy.2026.1.Hs.gmt",
    "reactome": "c2.cp.reactome.2026.1.Hs.gmt",
    "gobp": "c5.go.bp.2026.1.Hs.gmt",
}
NPERM, MIN_SIZE, MAX_SIZE, SEED = 1000, 10, 500, 1234
TOP_N, P_GATE = 15, 0.05
N_CPUS = int(os.environ.get("PBV_CPUS", "8"))

_GS = None


def gene_sets() -> dict[str, dict[str, list[str]]]:
    global _GS
    if _GS is None:
        _GS = {k: gp.read_gmt(os.path.join(GMT_DIR, f)) for k, f in LIB_FILES.items()}
    return _GS


# ---- datasets ---------------------------------------------------------------------------
# Each contrast dataset: metadata builder, design, (factor, test level, reference level).
CONTRASTS = {
    "GSE226163": "iVC vs iEC — iPSC smooth-muscle vs endothelial (3 v 3)",
    "GSE323136": "CF vs healthy — nasal curettage (5 v 5)",
    "GSE296381": "AML relapse vs diagnosis — marrow, paired (10 pairs)",
    "GSE295599": "tumour day-10 vs baseline — radiotherapy biopsy, paired (7 pairs)",
}
PLATFORM_ONLY = ["GSE217517"]  # no biological contrast; used for the platform-gap ranking only

CF = {"subj_001", "subj_003", "subj_004", "subj_011", "subj_015"}  # GEO GSE323136 'disease state'


def load_pair(ds: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    b = pd.read_csv(os.path.join(PL, ds, "bulk.csv"), index_col=0)
    p = pd.read_csv(os.path.join(PL, ds, "pseudobulk.csv"), index_col=0)
    b = b[~b.index.duplicated()]
    p = p[~p.index.duplicated()]
    g = b.index.intersection(p.index)
    cols = list(b.columns)
    assert cols == list(p.columns), ds
    return b.loc[g, cols].round().astype(int), p.loc[g, cols].round().astype(int)


def design_for(ds: str, cols: list[str]):
    """-> (metadata, design, contrast, pair_column or None). Paired sets keep complete pairs only."""
    if ds == "GSE226163":
        md = pd.DataFrame({"cond": [c[:3] for c in cols]}, index=cols)
        return md, "~cond", ("cond", "iVC", "iEC"), None
    if ds == "GSE323136":
        md = pd.DataFrame({"cond": ["CF" if c in CF else "HV" for c in cols]}, index=cols)
        return md, "~cond", ("cond", "CF", "HV"), None
    if ds == "GSE296381":
        md = pd.DataFrame({"pat": [c.split("_")[0] for c in cols],
                           "cond": [c.split("_")[1] for c in cols]}, index=cols)
        return md, "~pat + cond", ("cond", "RL", "DX"), "pat"
    if ds == "GSE295599":
        md = pd.DataFrame({"pat": [c[:-3] if c.endswith("d10") else c[:-1] for c in cols],
                           "cond": ["d10" if c.endswith("d10") else "b" for c in cols]}, index=cols)
        n = md.groupby("pat").cond.nunique()
        md = md[md.pat.isin(n[n == 2].index)]
        return md, "~pat + cond", ("cond", "d10", "b"), "pat"
    raise KeyError(ds)


# ---- DE ---------------------------------------------------------------------------------
def deseq(counts: pd.DataFrame, md: pd.DataFrame, design: str, contrast) -> pd.DataFrame:
    """counts genes x samples (int). Pittsburgh prefilter rowSums >= 10. Returns results_df."""
    counts = counts.loc[:, md.index]
    counts = counts[counts.sum(axis=1) >= 10]
    inf = DefaultInference(n_cpus=N_CPUS)
    dds = DeseqDataSet(counts=counts.T, metadata=md, design=design, inference=inf, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=list(contrast), inference=inf, quiet=True)
    st.summary()
    res = st.results_df.copy()
    res.index.name = "gene"
    return res


def ranks(res: pd.DataFrame) -> pd.Series:
    """Pittsburgh ranking: Wald stat, MT-* dropped, symbol-deduplicated by highest baseMean."""
    d = res.dropna(subset=["stat"]).copy()
    d["sym"] = d.index.astype(str).str.upper()
    d = d[~d.sym.str.startswith("MT-")]
    d = d.sort_values("baseMean", ascending=False).drop_duplicates("sym")
    return d.set_index("sym")["stat"].sort_values(ascending=False).astype(float)


def gsea(rnk: pd.Series, libs=None, nperm: int = NPERM, threads: int = N_CPUS) -> dict[str, pd.DataFrame]:
    out = {}
    for name, gs in gene_sets().items():
        if libs is not None and name not in libs:
            continue
        pre = gp.prerank(rnk=rnk, gene_sets=gs, permutation_num=nperm, min_size=MIN_SIZE,
                         max_size=MAX_SIZE, outdir=None, seed=SEED, threads=threads, verbose=False)
        d = pre.res2d.rename(columns={"NOM p-val": "p", "FDR q-val": "fdr", "Lead_genes": "lead"})
        d = d[["Term", "ES", "NES", "p", "fdr", "lead"]].copy()
        for c in ("ES", "NES", "p", "fdr"):
            d[c] = pd.to_numeric(d[c], errors="coerce")
        out[name] = d.set_index("Term")
    return out


def _running_extremes(P: np.ndarray, w: np.ndarray, N: int) -> np.ndarray:
    """ES of gene sets given sorted hit positions P (m x k) in a ranked list of N genes.

    Same statistic as gseapy (weighted KS, weight=1): the running sum peaks at a hit and bottoms
    out just before one, and it always ends at exactly 0, so only those 2k points are needed.
    """
    k = P.shape[1]
    ww = w[P]
    cw = np.cumsum(ww, axis=1)
    nr = cw[:, -1:]
    miss = (P - np.arange(k)[None, :]) / (N - k)  # misses before (and at) hit j
    hi = np.maximum((cw / nr - miss).max(axis=1), 0.0)
    lo = np.minimum(((cw - ww) / nr - miss).min(axis=1), 0.0)
    return np.where(np.abs(hi) > np.abs(lo), hi, lo)


def fast_gsea(rnk: pd.Series, libs=None, nperm: int = NPERM, seed: int = SEED) -> dict[str, pd.DataFrame]:
    """numpy preranked GSEA reproducing gseapy's ES / NES / nominal p / FDR formulas.

    Difference from gseapy: the gene-set-permutation null is shared by all sets of the same
    size (fgsea-style), which is statistically equivalent and ~20x faster on GO BP.
    """
    rnk = rnk[~rnk.index.duplicated()].sort_values(ascending=False)
    N = len(rnk)
    w = np.abs(rnk.to_numpy(float))
    pos_of = pd.Series(np.arange(N), index=rnk.index)
    out = {}
    for name, gs in gene_sets().items():
        if libs is not None and name not in libs:
            continue
        terms, poss = [], []
        for term, genes in gs.items():
            idx = pos_of.reindex(pd.unique(np.asarray(genes, dtype=object))).dropna().to_numpy(int)
            if MIN_SIZE <= len(idx) <= MAX_SIZE:
                terms.append(term)
                poss.append(np.sort(idx))
        sizes = np.array([len(x) for x in poss])
        rng = np.random.default_rng(seed)
        kmax = sizes.max()
        R = np.stack([rng.choice(N, size=kmax, replace=False) for _ in range(nperm)])
        es = np.empty(len(terms))
        nes = np.empty(len(terms))
        pv = np.empty(len(terms))
        null_norm, null_owner = [], []
        for k in np.unique(sizes):
            members = np.flatnonzero(sizes == k)
            e = _running_extremes(np.stack([poss[i] for i in members]), w, N)
            null = _running_extremes(np.sort(R[:, :k], axis=1), w, N)
            pos_mean = null[null >= 0].mean() if (null >= 0).any() else np.nan
            neg_mean = null[null < 0].mean() if (null < 0).any() else np.nan
            es[members] = e
            nes[members] = np.where(e >= 0, e / pos_mean, -e / neg_mean)
            pv[members] = np.where(e >= 0, (null[None, :] >= e[:, None]).sum(1) / max((null >= 0).sum(), 1),
                                   (null[None, :] < e[:, None]).sum(1) / max((null < 0).sum(), 1))
            nn = np.where(null >= 0, null / pos_mean, -null / neg_mean)
            null_norm.append(nn)
            null_owner.append(len(members))
        # gseapy FDR: pooled normalised nulls of every tested set vs the observed NES
        nvals = np.sort(np.concatenate([np.repeat(nn[None, :], c, axis=0).ravel()
                                        for nn, c in zip(null_norm, null_owner)]))
        nnes = np.sort(nes)
        fdr = np.empty(len(terms))
        for i, x in enumerate(nes):
            if x >= 0:
                all_pos = len(nvals) - np.searchsorted(nvals, 0, "left")
                all_hi = len(nvals) - np.searchsorted(nvals, x, "left")
                obs_pos = len(nnes) - np.searchsorted(nnes, 0, "left")
                obs_hi = len(nnes) - np.searchsorted(nnes, x, "left")
            else:
                all_pos = np.searchsorted(nvals, 0, "left")
                all_hi = np.searchsorted(nvals, x, "right")
                obs_pos = np.searchsorted(nnes, 0, "left")
                obs_hi = np.searchsorted(nnes, x, "right")
            with np.errstate(divide="ignore", invalid="ignore"):
                f = (all_hi / all_pos) / (obs_hi / obs_pos)
            fdr[i] = min(f, 1.0) if np.isfinite(f) else 1.0
        out[name] = pd.DataFrame({"ES": es, "NES": nes, "p": pv, "fdr": fdr, "size": sizes},
                                 index=pd.Index(terms, name="Term"))
    return out


def top_lists(t: pd.DataFrame, n: int = TOP_N, gate: float | None = P_GATE) -> dict[str, list[str]]:
    """Pittsburgh figure rule: nominal p < gate, then top n each way by NES."""
    d = t if gate is None else t[t.p < gate]
    up = d[d.NES > 0].sort_values("NES", ascending=False).index[:n]
    dn = d[d.NES < 0].sort_values("NES", ascending=True).index[:n]
    return {"top15_up": list(up), "bottom15_down": list(dn)}


def compare(ref: pd.DataFrame, alt: pd.DataFrame, gate: float | None = P_GATE) -> list[dict]:
    """Compare two GSEA tables (same library). ref = the analysis you would have reported."""
    common = ref.index.intersection(alt.index)
    R, A = ref.loc[common], alt.loc[common]
    nes_r = R.NES.corr(A.NES)
    nes_rho = R.NES.corr(A.NES, method="spearman")
    lr, la = top_lists(R, gate=gate), top_lists(A, gate=gate)
    rows = []
    for direction, sign in (("top15_up", 1), ("bottom15_down", -1)):
        lref, lalt = lr[direction], la[direction]
        ov = len(set(lref) & set(lalt))
        # where do the ref's listed sets land in alt's ranking of the same direction?
        arank = (sign * A.NES).rank(ascending=False)
        held = [s for s in lref if np.sign(A.at[s, "NES"]) == sign and A.at[s, "p"] < P_GATE]
        rows.append(dict(
            direction=direction, n_tested=len(common), nes_r=nes_r, nes_rho=nes_rho,
            n_ref=len(lref), n_alt=len(lalt), overlap=ov,
            chance=len(lref) * len(lalt) / max(len(common), 1),
            same_sign=int(sum(np.sign(A.at[s, "NES"]) == sign for s in lref)),
            still_called=len(held),  # same direction AND nominal p < 0.05 in alt
            in_alt_top50=int(sum(arank.at[s] <= 50 for s in lref)),
            median_alt_rank=float(np.median([arank.at[s] for s in lref])) if lref else np.nan,
        ))
    return rows


def save_gsea(tabs: dict[str, pd.DataFrame], prefix: str) -> None:
    for lib, t in tabs.items():
        t.to_csv(f"{prefix}_{lib}_gsea.csv")


def load_gsea(prefix: str, libs=None) -> dict[str, pd.DataFrame]:
    out = {}
    for lib in LIB_FILES:
        if libs is not None and lib not in libs:
            continue
        f = f"{prefix}_{lib}_gsea.csv"
        if os.path.exists(f):
            out[lib] = pd.read_csv(f, index_col=0)
    return out
