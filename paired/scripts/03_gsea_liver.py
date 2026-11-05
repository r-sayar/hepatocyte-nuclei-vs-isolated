#!/usr/bin/env python3
"""03_gsea_liver.py — does GSEA of inflammation change when liver samples are pseudobulked from
scRNA-seq or snRNA-seq instead of bulk-sequenced? All annotated genes (the default pipeline).

Parts (outputs in results/<part>/, summaries in results/*.csv):
  P1  ENCODE4 liver, bulk vs sn-pseudobulk, same donor+lobe, ~pair + platform   (+ = higher in sn)
  P2  Andrews livers, sc vs sn pseudobulk, same liver, ~liver + modality         (+ = higher in sc)
      all 10 livers / healthy 4 (both 3') / PSC-PBC 6 (sc 5' vs sn 3')
  P3  PSC/PBC (6) vs healthy (4): the same contrast run on sc-pseudobulk and on sn-pseudobulk
  P4  checks: leave-one-liver-out noise floor for P3; inflammation sets re-scored with the
      dissociation-stress genes removed; per-pair sign consistency; author vs rule cells
04 repeats P1-P3 on protein-coding genes and with a length-adjusted ranking; 05 splits P2 by cell type.
"""
import os
import numpy as np, pandas as pd
import liver_common as L
from liver_common import C

RES = L.RES
ALL_INFL, ALL_SUM, ALL_NODIS, ALL_PAIR, CMP = [], [], [], [], []


def record(tabs, res, label):
    ALL_INFL.append(L.infl_table(tabs, label))
    ALL_SUM.append(L.infl_summary(tabs, label))
    ALL_NODIS.append(L.without_dissoc(res, label))


def pair_signs(X, Y, units, label):
    """Per-unit Hallmark NES on the log2 CPM ratio X/Y (genes with CPM>=1 in both)."""
    rows = []
    for u in units:
        a, b = X[u] / X[u].sum() * 1e6, Y[u] / Y[u].sum() * 1e6
        ok = (a >= 1) & (b >= 1)
        lr = np.log2(a[ok] + 0.5) - np.log2(b[ok] + 0.5)
        lr.index = [L.SYM.get(e, e) for e in lr.index]
        lr = lr[~lr.index.str.upper().str.startswith("MT-")].groupby(level=0).mean().sort_values(ascending=False)
        t = C.fast_gsea(lr, libs=["hallmark"])["hallmark"]
        rows += [dict(analysis=label, unit=u, term=x, NES=t.at[x, "NES"], p=t.at[x, "p"]) for x in L.INFL_H if x in t.index]
    return pd.DataFrame(rows)


# P1 ENCODE bulk vs sn
Bk, Sn = L.load_encode()
pairs = list(Bk.columns)
M = pd.concat([Bk[pairs].add_suffix("__bulk"), Sn[pairs].add_suffix("__sn")], axis=1)
res, tabs = L.run(M, L.paired_md(list(M.columns), "pair", "platform"), "~pair + platform",
                  ("platform", "sn", "bulk"), f"{RES}/P1/encode_sn_vs_bulk")
record(tabs, res, "P1 ENCODE sn vs bulk")
ALL_PAIR.append(pair_signs(Sn, Bk, pairs, "P1 ENCODE sn vs bulk"))

# P2 Andrews sc vs sn
PB, LIB = L.load_andrews("rule")
livers = sorted(LIB.liver.unique())
healthy = [l for l in livers if l.startswith("Healthy")]
disease = [l for l in livers if l not in healthy]
print(f"Andrews shared ENSG across all libraries: {len(PB)}", flush=True)
for name, sub in (("all", livers), ("healthy", healthy), ("disease", disease)):
    cols = [f"{l}_{m}" for l in sub for m in ("sc", "sn")]
    res, tabs = L.run(PB[cols], L.paired_md(cols, "liver", "modality"), "~liver + modality",
                      ("modality", "sc", "sn"), f"{RES}/P2/andrews_sc_vs_sn_{name}")
    record(tabs, res, f"P2 Andrews sc vs sn ({name}, n={len(sub)})")
SC = PB[[f"{l}_sc" for l in livers]].rename(columns=lambda c: c.rsplit("_", 1)[0])
SN = PB[[f"{l}_sn" for l in livers]].rename(columns=lambda c: c.rsplit("_", 1)[0])
ALL_PAIR.append(pair_signs(SC, SN, livers, "P2 Andrews sc vs sn"))

# P3 PSC/PBC vs healthy, on sc and on sn
md3 = pd.DataFrame({"group": ["disease" if l in disease else "healthy" for l in livers]}, index=livers)
res_sc, tabs_sc = L.run(SC, md3, "~group", ("group", "disease", "healthy"), f"{RES}/P3/disease_vs_healthy_sc")
res_sn, tabs_sn = L.run(SN, md3, "~group", ("group", "disease", "healthy"), f"{RES}/P3/disease_vs_healthy_sn")
record(tabs_sc, res_sc, "P3 disease vs healthy, sc-pseudobulk")
record(tabs_sn, res_sn, "P3 disease vs healthy, sn-pseudobulk")
for lib in C.LIB_FILES:
    for ref, alt, tr, ta in (("sn", "sc", tabs_sn, tabs_sc), ("sc", "sn", tabs_sc, tabs_sn)):
        for row in C.compare(tr[lib], ta[lib]):
            it = [x for x in tr[lib].index.intersection(ta[lib].index) if L.is_infl(lib, x)]
            row.update(lib=lib, ref=ref, alt=alt, n_infl=len(it),
                       infl_nes_r=float(tr[lib].loc[it, "NES"].corr(ta[lib].loc[it, "NES"])) if len(it) > 2 else np.nan)
            CMP.append(row)
gene_r = float(res_sc.set_index("ensg").stat.corr(res_sn.set_index("ensg").stat))

# P4 leave-one-liver-out noise floor for P3 (Hallmark)
jk = []
for drop in livers:
    keep = [l for l in livers if l != drop]
    for plat, X in (("sc", SC), ("sn", SN)):
        _, t = L.run(X[keep], md3.loc[keep], "~group", ("group", "disease", "healthy"),
                     f"{RES}/P4/jk_{plat}_minus_{drop}", libs=["hallmark"])
        h = t["hallmark"]
        jk += [dict(platform=plat, dropped=drop, term=x, NES=h.at[x, "NES"], p=h.at[x, "p"]) for x in L.INFL_H if x in h.index]
JK = pd.DataFrame(jk)

# P4 author cells vs rule cells (healthy sc vs sn)
AU, _ = L.load_andrews("author")
cols = [f"{l}_{m}" for l in healthy for m in ("sc", "sn")]
md = L.paired_md(cols, "liver", "modality")
_, tabs_au = L.run(AU[cols], md, "~liver + modality", ("modality", "sc", "sn"), f"{RES}/P4/andrews_sc_vs_sn_healthy_authorcells")
_, tabs_ru = L.run(PB[cols], md, "~liver + modality", ("modality", "sc", "sn"), f"{RES}/P2/andrews_sc_vs_sn_healthy")
cell_check = {lib: float(tabs_au[lib].NES.corr(tabs_ru[lib].NES.reindex(tabs_au[lib].index))) for lib in C.LIB_FILES}
ALL_INFL.append(L.infl_table(tabs_au, "P4 healthy sc vs sn, author cells"))

pd.concat(ALL_INFL).to_csv(f"{RES}/inflammation_sets.csv", index=False)
pd.concat(ALL_SUM).to_csv(f"{RES}/inflammation_summary.csv", index=False)
pd.concat(ALL_NODIS).to_csv(f"{RES}/inflammation_without_dissoc_genes.csv", index=False)
pd.concat(ALL_PAIR).to_csv(f"{RES}/per_pair_signs.csv", index=False)
pd.DataFrame(CMP).to_csv(f"{RES}/P3_sc_vs_sn_list_comparison.csv", index=False)
JK.to_csv(f"{RES}/P4_jackknife_hallmark.csv", index=False)
pd.Series(cell_check, name="NES_r_author_vs_rule").to_csv(f"{RES}/P4_author_vs_rule_cells.csv")

pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
H = pd.concat(ALL_INFL)
print("\n== Hallmark inflammation NES by analysis\n", H[H.lib == "hallmark"].pivot_table(index="term", columns="analysis", values="NES").round(2).T.to_string())
print("\n== inflammation-set summary\n", pd.concat(ALL_SUM).round(3).to_string(index=False))
print("\n== without dissociation genes\n", pd.concat(ALL_NODIS).round(3).to_string(index=False))
print("\n== P3 sc vs sn list comparison (ref=sn)\n", pd.DataFrame(CMP).query("ref=='sn'").round(2).to_string(index=False))
print("\n== P3 gene-level Wald stat r (sc vs sn):", round(gene_r, 3))
print("\n== jackknife NES range\n", JK.groupby(["term", "platform"]).NES.agg(["min", "median", "max"]).round(2).to_string())
print("\n== author vs rule cells NES r:", {k: round(v, 3) for k, v in cell_check.items()})
ps = pd.concat(ALL_PAIR)
print("\n== per-unit sign consistency (NES>0 / total)\n", ps.assign(pos=ps.NES > 0).groupby(["analysis", "term"]).pos.agg(["sum", "count"]).to_string())
