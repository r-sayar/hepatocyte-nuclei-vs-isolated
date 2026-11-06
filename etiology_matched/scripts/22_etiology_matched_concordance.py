#!/usr/bin/env python3
"""22_etiology_matched_concordance.py — isolated-hepatocyte bulk (MASH stratum) vs MASH hepatocyte nuclei.

The PSC/PBC nuclei of 16_concordance.py differ from the bulk cohort in etiology. Gribben et al. 2024 (GSE202379)
gives hepatocyte nuclei from end-stage MASH explants vs healthy livers, the etiology-matched comparator. Two
independent nuclei datasets with different etiologies then classify every bulk-vs-nuclei disagreement:
  platform   the bulk disagrees with BOTH nuclei datasets in the same way  -> isolated cells vs nuclei
  etiology   only the PSC/PBC nuclei disagree                                -> cholestatic vs MASH disease
Yardstick: the bulk MASH stratum refitted on 50 random 5-vs-4 draws (15_bulk_profiles_and_downsample.R), the size of
the end-stage vs healthy nuclei design. Scoring: fast_gsea on protein-coding limma t (as in 16_concordance.py).
Outputs: results/em_*.csv
"""
import os, sys
from multiprocessing import Pool
import numpy as np, pandas as pd

ROOT = "/srv/public/sayar99/liver_strat_gsea"
sys.path.insert(0, os.path.join(ROOT, "scripts", "ref"))
import common as C  # noqa: E402

C.GMT_DIR = os.path.join(ROOT, "genesets")
D, R = os.path.join(ROOT, "data"), os.path.join(ROOT, "results")
LIBS = list(C.LIB_FILES)


def rank_from(d, sym_col=None):
    d = d[d.protein_coding.astype(bool)].copy()
    d["sym"] = d[sym_col] if sym_col else d.index
    d = d[d.sym.notna() & ~d.sym.astype(str).str.upper().str.startswith("MT-")]
    if "AveExpr" in d:
        d = d.sort_values("AveExpr", ascending=False)
    d = d.drop_duplicates("sym")
    return d.set_index("sym")["t"].astype(float).sort_values(ascending=False)


rk = {"bulk_MASH": rank_from(pd.read_csv(os.path.join(D, "bulk_MASH_limma.csv"), index_col=0)),
      "bulk_esld_all": rank_from(pd.read_csv(os.path.join(D, "bulk_esld_all_limma.csv"), index_col=0))}
for k in ("endstage", "cirrMASH", "NASH", "fibrosis_trend"):
    rk[f"mash_decont_{k}"] = rank_from(pd.read_csv(os.path.join(R, "gribben", f"gribben_hep_decont_{k}_limma.csv"), index_col=0))
    rk[f"mash_nuclei_{k}"] = rank_from(pd.read_csv(os.path.join(R, "gribben", f"gribben_hep_{k}_limma.csv"), index_col=0))
rk["mash_allnuclei_endstage"] = rank_from(pd.read_csv(os.path.join(R, "gribben", "gribben_all_endstage_limma.csv"), index_col=0))
for k in ("sn_hep", "sc_hep", "sn_hepdecont", "sc_hepdecont"):
    rk[f"pscpbc_{k}"] = rank_from(pd.read_csv(os.path.join(R, f"pb_{k}_limma.csv"), index_col=0), "symbol")
ds = pd.read_csv(os.path.join(D, "bulk_downsample_MASH_5v4_t.csv.gz"))
ds["protein_coding"] = ds.gene.isin(rk["bulk_MASH"].index)
draw_rk = [rank_from(g.set_index("gene")) for _, g in ds.groupby("draw")]

with Pool(8) as pool:
    names = list(rk)
    tabs = dict(zip(names, pool.map(C.fast_gsea, [rk[n] for n in names])))
    dtabs = pool.map(C.fast_gsea, draw_rk)
os.makedirs(os.path.join(R, "pathway_nes"), exist_ok=True)
for n, t in tabs.items():
    for lib, tt in t.items():
        tt.to_csv(os.path.join(R, "pathway_nes", f"em_{n}_{lib}.csv"))


def calls(t):
    return pd.Series(np.where(t.p < 0.05, np.sign(t.NES), 0), index=t.index)


def rho(a, b, terms=None):
    g = a.index.intersection(b.index)
    if terms is not None:
        g = g.intersection(terms)
    return float(a.NES[g].corr(b.NES[g], method="spearman"))


def lost(ref, alt):
    """fraction of pathways called in ref that are not called the same way in alt."""
    r, a = calls(ref), calls(alt)
    hit = r[r != 0]
    aa = a.reindex(hit.index).dropna()
    return float((aa != hit.loc[aa.index]).mean()), int(len(aa)), int(((aa != 0) & (aa != hit.loc[aa.index])).sum())


def top_overlap(ref, alt, which):
    return len(set(C.top_lists(ref, n=25)[which]) & set(C.top_lists(alt, n=25)[which]))


B, N = "bulk_MASH", "mash_nuclei_endstage"
rows = []
for lib in LIBS:
    ref = tabs[B][lib]
    called = ref.index[ref.p < 0.05]
    for other in [n for n in tabs if n != B]:
        o = tabs[other][lib]
        lb, nb, ob = lost(ref, o)
        lr, nr, orr = lost(o, ref)
        rows.append(dict(lib=lib, other=other, rho_all=rho(ref, o), rho_bulk_called=rho(ref, o, called),
                         bulk_to_other_lost=lb, bulk_called_tested=nb, bulk_to_other_opposite=ob,
                         other_to_bulk_lost=lr, other_called_tested=nr, other_to_bulk_opposite=orr,
                         top25_up_overlap=top_overlap(ref, o, "top15_up"), top25_down_overlap=top_overlap(ref, o, "bottom15_down")))
    X = pd.DataFrame([dict(rho_all=rho(ref, t[lib]), rho_bulk_called=rho(ref, t[lib], called),
                           bulk_to_other_lost=lost(ref, t[lib])[0], other_to_bulk_lost=lost(t[lib], ref)[0],
                           top25_up_overlap=top_overlap(ref, t[lib], "top15_up"), top25_down_overlap=top_overlap(ref, t[lib], "bottom15_down"))
                      for t in dtabs])
    rows.append(dict(lib=lib, other="bulk MASH 5v4 draws (median)", **{c: float(X[c].median()) for c in X.columns},
                     **{f"{c}_p5": float(X[c].quantile(.05)) for c in X.columns}, **{f"{c}_p95": float(X[c].quantile(.95)) for c in X.columns}))
    # nuclei vs nuclei across etiologies (same platform, different disease)
    rows.append(dict(lib=lib, other="MASH nuclei vs PSC/PBC nuclei", rho_all=rho(tabs[N][lib], tabs["pscpbc_sn_hep"][lib]),
                     top25_up_overlap=top_overlap(tabs[N][lib], tabs["pscpbc_sn_hep"][lib], "top15_up"),
                     top25_down_overlap=top_overlap(tabs[N][lib], tabs["pscpbc_sn_hep"][lib], "bottom15_down")))
SUM = pd.DataFrame(rows)

# every pathway called by the bulk MASH stratum or by MASH hepatocyte nuclei, classified with the PSC/PBC nuclei
cls = []
for lib in LIBS:
    b, n, p = calls(tabs[B][lib]), calls(tabs[N][lib]), calls(tabs["pscpbc_sn_hep"][lib])
    dsame = pd.concat([calls(t[lib]) for t in dtabs], axis=1)
    for term in b.index[b != 0].union(n.index[n != 0]):
        bb, nn, pp = b.get(term, np.nan), n.get(term, np.nan), p.get(term, np.nan)
        if pd.isna(bb) or pd.isna(nn):
            continue
        agree = bb == nn
        both_nuclei_disagree = (not agree) and pd.notna(pp) and pp == nn
        cls.append(dict(lib=lib, term=term, bulk_MASH_call=bb, mash_nuclei_call=nn, pscpbc_nuclei_call=pp,
                        status="agree" if agree else ("platform (both nuclei sets)" if both_nuclei_disagree else "MASH nuclei only"),
                        small_bulk_calls_bulk_direction=float((dsame.loc[term] == bb).mean()) if (bb != 0 and term in dsame.index) else np.nan,
                        bulk_NES=tabs[B][lib].NES.get(term), mash_nuclei_NES=tabs[N][lib].NES.get(term),
                        pscpbc_nuclei_NES=tabs["pscpbc_sn_hep"][lib].NES.get(term), pscpbc_cells_NES=tabs["pscpbc_sc_hep"][lib].NES.get(term),
                        mash_allnuclei_NES=tabs["mash_allnuclei_endstage"][lib].NES.get(term)))
CLS = pd.DataFrame(cls)
SUM.to_csv(os.path.join(R, "em_summary.csv"), index=False)
CLS.to_csv(os.path.join(R, "em_pathway_classes.csv"), index=False)

pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(SUM[["lib", "other", "rho_all", "rho_bulk_called", "bulk_to_other_lost", "other_to_bulk_lost", "top25_up_overlap", "top25_down_overlap"]].round(2).to_string(index=False))
print(SUM[SUM.other.str.contains("draws")].filter(regex="lib|_p5|_p95").round(2).to_string(index=False))
print(CLS.groupby(["lib", "status"]).size().unstack(fill_value=0).to_string())
