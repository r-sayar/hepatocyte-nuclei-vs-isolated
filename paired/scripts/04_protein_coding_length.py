#!/usr/bin/env python3
"""04_protein_coding_length.py — is the inflammation shift specific, or just gene class / gene length?

03 ranks every annotated gene. Nuclear libraries load intronic, long and non-coding genes, so in
03 nearly every protein-coding pathway lands on the bulk side. Here:
  pc      DESeq2 on protein-coding genes only (the universe Hallmark sets are drawn from)
  lenadj  pc Wald stat minus the median stat of genes with similar exonic length and genomic span
          (10 x 10 quantile grid); asks whether inflammation genes move more than their length predicts
for P1 (ENCODE sn vs bulk) and P2 (Andrews sc vs sn: all / healthy / disease), and the pc version of
P3 (PSC/PBC vs healthy on sc vs on sn) with its leave-one-liver-out floor.
Outputs: results/PC/..., results/pc_*.csv
"""
import os
import numpy as np, pandas as pd
import liver_common as L
from liver_common import C

OUT = os.path.join(L.RES, "PC")
INFL, SUMM, NODIS, TILT, CMP, JK = [], [], [], [], [], []


def gsea_cached(rnk, pre, libs=None):
    want = libs or list(C.LIB_FILES)
    tabs = C.load_gsea(pre, libs=want)
    if len(tabs) < len(want):
        tabs = L.gsea(rnk, libs=want)
        C.save_gsea(tabs, pre)
        print(f"[gsea] {os.path.basename(pre)}", flush=True)
    return tabs


def both(counts, md, design, contrast, pre, label):
    res, tabs = L.run(L.protein_coding(counts), md, design, contrast, pre)
    tabs_la = gsea_cached(L.length_adjusted_ranks(res), pre + "_lenadj")
    for lab, t in ((f"{label} | pc", tabs), (f"{label} | pc+length-adj", tabs_la)):
        INFL.append(L.infl_table(t, lab))
        SUMM.append(L.infl_summary(t, lab))
    NODIS.append(L.without_dissoc(res, f"{label} | pc"))
    NODIS.append(L.without_dissoc(res, f"{label} | pc+length-adj", rank_fn=L.length_adjusted_ranks))
    TILT.append(dict(analysis=label, **L.length_tilt(res)))
    return res, tabs


# P1 ENCODE sn vs bulk
Bk, Sn = L.load_encode()
pairs = list(Bk.columns)
M = pd.concat([Bk[pairs].add_suffix("__bulk"), Sn[pairs].add_suffix("__sn")], axis=1)
both(M, L.paired_md(list(M.columns), "pair", "platform"), "~pair + platform", ("platform", "sn", "bulk"),
     f"{OUT}/P1_encode_sn_vs_bulk", "P1 ENCODE sn vs bulk (5 pairs)")

# P2 Andrews sc vs sn
PB, LIB = L.load_andrews("rule")
livers = sorted(LIB.liver.unique())
healthy = [l for l in livers if l.startswith("Healthy")]
disease = [l for l in livers if l not in healthy]
for name, sub in (("all", livers), ("healthy", healthy), ("disease", disease)):
    cols = [f"{l}_{m}" for l in sub for m in ("sc", "sn")]
    both(PB[cols], L.paired_md(cols, "liver", "modality"), "~liver + modality", ("modality", "sc", "sn"),
         f"{OUT}/P2_andrews_sc_vs_sn_{name}", f"P2 Andrews sc vs sn, {name} (n={len(sub)})")

# P3 PSC/PBC vs healthy, sc vs sn, protein-coding
SC = PB[[f"{l}_sc" for l in livers]].rename(columns=lambda c: c.rsplit("_", 1)[0])
SN = PB[[f"{l}_sn" for l in livers]].rename(columns=lambda c: c.rsplit("_", 1)[0])
md3 = pd.DataFrame({"group": ["disease" if l in disease else "healthy" for l in livers]}, index=livers)
r_sc, t_sc = L.run(L.protein_coding(SC), md3, "~group", ("group", "disease", "healthy"), f"{OUT}/P3_disease_vs_healthy_sc")
r_sn, t_sn = L.run(L.protein_coding(SN), md3, "~group", ("group", "disease", "healthy"), f"{OUT}/P3_disease_vs_healthy_sn")
for lab, t, r in (("P3 disease vs healthy, sc-pseudobulk | pc", t_sc, r_sc), ("P3 disease vs healthy, sn-pseudobulk | pc", t_sn, r_sn)):
    INFL.append(L.infl_table(t, lab))
    SUMM.append(L.infl_summary(t, lab))
    NODIS.append(L.without_dissoc(r, lab))
for lib in C.LIB_FILES:
    for ref, alt, tr, ta in (("sn", "sc", t_sn, t_sc), ("sc", "sn", t_sc, t_sn)):
        for row in C.compare(tr[lib], ta[lib]):
            it = [x for x in tr[lib].index.intersection(ta[lib].index) if L.is_infl(lib, x)]
            row.update(lib=lib, ref=ref, alt=alt, n_infl=len(it),
                       infl_nes_r=float(tr[lib].loc[it, "NES"].corr(ta[lib].loc[it, "NES"])) if len(it) > 2 else np.nan)
            CMP.append(row)
gene_r = float(r_sc.set_index("ensg").stat.corr(r_sn.set_index("ensg").stat))

# leave-one-liver-out, both platforms, Hallmark
for drop in livers:
    keep = [l for l in livers if l != drop]
    for plat, X in (("sc", SC), ("sn", SN)):
        _, t = L.run(L.protein_coding(X[keep]), md3.loc[keep], "~group", ("group", "disease", "healthy"),
                     f"{OUT}/P4_jk/{plat}_minus_{drop}", libs=["hallmark"])
        h = t["hallmark"]
        for term in L.INFL_H:
            if term in h.index:
                JK.append(dict(platform=plat, dropped=drop, term=term, NES=h.at[term, "NES"], p=h.at[term, "p"]))

pd.concat(INFL).to_csv(f"{L.RES}/pc_inflammation_sets.csv", index=False)
pd.concat(SUMM).to_csv(f"{L.RES}/pc_inflammation_summary.csv", index=False)
pd.concat(NODIS).to_csv(f"{L.RES}/pc_inflammation_without_dissoc_genes.csv", index=False)
pd.DataFrame(TILT).to_csv(f"{L.RES}/pc_length_tilt.csv", index=False)
pd.DataFrame(CMP).to_csv(f"{L.RES}/pc_P3_sc_vs_sn_list_comparison.csv", index=False)
pd.DataFrame(JK).to_csv(f"{L.RES}/pc_P4_jackknife_hallmark.csv", index=False)

pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
H = pd.concat(INFL)
print(H[H.lib == "hallmark"].pivot_table(index="term", columns="analysis", values="NES").round(2).T.to_string())
print(pd.concat(SUMM).round(3).to_string(index=False))
print(pd.concat(NODIS).round(3).to_string(index=False))
print(pd.DataFrame(TILT).round(3).to_string(index=False))
print(pd.DataFrame(CMP).query("ref=='sn'").round(2).to_string(index=False))
print("P3 pc gene-level Wald r sc vs sn:", round(gene_r, 3))
J = pd.DataFrame(JK)
print(J.groupby(["term", "platform"]).NES.agg(["min", "median", "max"]).round(2).to_string())
