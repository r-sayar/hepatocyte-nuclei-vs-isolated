"""Hepatocyte disease programme defined from nuclei only (public data): genes moving the same way, t > 2.5, in BOTH
decontaminated hepatocyte-nucleus contrasts (MASH fibrosis trend within the biopsy source, Gribben et al.;
PSC/PBC vs healthy, Andrews et al.), restricted to hepatocyte-made genes (MacParland atlas hepatocyte fraction >= 0.5)."""
import pandas as pd, numpy as np
co = pd.read_csv("data/cell_of_origin.csv").drop_duplicates("symbol").set_index("symbol")["hep_fraction"]


def load(f, symcol=None):
    d = pd.read_csv(f, index_col=0)
    d = d[d.protein_coding.astype(bool)]
    if symcol:
        d = d.set_index(symcol)
    return d[~d.index.duplicated()]["t"]


a = load("results/gribben/gribben_hep_decont_fibrosis_trend_limma.csv")
b = load("results/pb_sn_hepdecont_limma.csv", "symbol")
g = a.index.intersection(b.index)
T = pd.DataFrame({"mash_nuc_trend_t": a[g], "pscpbc_nuc_t": b[g], "hep_frac": co.reindex(g)})
TH = 2.5
out = []
for sign, lab in ((1, "up"), (-1, "down")):
    hit = (sign * T.mash_nuc_trend_t > TH) & (sign * T.pscpbc_nuc_t > TH)
    exp = ((sign * T.mash_nuc_trend_t > TH).mean() * (sign * T.pscpbc_nuc_t > TH).mean()) * len(T)
    sel = T[hit & (T.hep_frac >= 0.5)].assign(direction=lab)
    print(f"{lab}: both nuclei sets {hit.sum()} (chance {exp:.1f}); hepatocyte-made {len(sel)}")
    out.append(sel)
S = pd.concat(out)
S.to_csv("results/nuclei_hepatocyte_programme.csv")
print(S.sort_values(["direction", "mash_nuc_trend_t"]).round(2).to_string())
