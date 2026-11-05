#!/usr/bin/env python3
"""23_decontaminate.py — is the ECM / EMT / immune signal in hepatocyte-nucleus pseudobulks contamination?

SoupX's logic (Young & Behjati 2020) at pseudobulk level. Each sample's hepatocyte pseudobulk is taken as true
hepatocyte expression plus a fraction alpha of the sample's ambient pool. The ambient pool is approximated by the
sample's all-nuclei pseudobulk (no empty droplets are deposited for GSE202379). Alpha is estimated per sample from
genes hepatocytes do not express, deliberately excluding ECM / stellate genes so the test is not circular:
    alpha = sum(hep CPM of G0) / sum(all-nuclei CPM of G0),  G0 = immune + endothelial markers
The predicted contamination (alpha x all-nuclei counts, scaled to the hepatocyte library) is subtracted, floored at
zero. Outputs feed the same limma / GSEA steps:
  Gribben et al. 2024 (GSE202379): data/gse202379/pb_hep_decont.csv
  Andrews et al. PSC/PBC + healthy:  data/pb_hepdecont_ensg.csv
  results/decont_alpha.csv (alpha per sample), results/decont_ecm_check.csv (observed vs predicted for ECM genes)
"""
import os
import numpy as np, pandas as pd

ROOT = "/srv/public/sayar99/liver_strat_gsea"
D, R = os.path.join(ROOT, "data"), os.path.join(ROOT, "results")
G0 = ["PTPRC", "CD3E", "CD2", "NKG7", "CD79A", "MS4A1", "LYZ", "CD68", "CD163", "MARCO",
      "PECAM1", "VWF", "STAB2", "CLEC4G", "FCN3"]
ECM = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "ACTA2", "PDGFRB", "SPARC", "COL6A3", "FBN1"]
ann = pd.read_csv(os.path.join(D, "ensg2symbol_v46.tsv"), sep="\t", header=None, names=["ensg", "symbol", "type"]).drop_duplicates("ensg").set_index("ensg")


def decontaminate(hep, allp, sym, label):
    """hep, allp: genes x samples counts (same index); sym: index -> symbol Series."""
    hep, allp = hep.align(allp, join="inner", axis=0)
    cols = hep.columns.intersection(allp.columns)
    hep, allp = hep[cols], allp[cols]
    s = sym.reindex(hep.index)
    hc, ac = hep / hep.sum() * 1e6, allp / allp.sum() * 1e6
    g0 = s.isin(G0).to_numpy()
    alpha = hc[g0].sum() / ac[g0].sum()
    pred = ac.mul(alpha, axis=1) * hep.sum() / 1e6       # predicted contaminating counts in the hepatocyte library
    clean = (hep - pred).clip(lower=0).round().astype(int)
    ecm = s.isin(ECM).to_numpy()
    chk = pd.DataFrame({"dataset": label, "sample": cols, "alpha": alpha.values,
                        "ecm_observed_cpm": hc[ecm].sum().values, "ecm_predicted_cpm": (ac[ecm].mul(alpha, axis=1)).sum().values})
    chk["ecm_explained_frac"] = chk.ecm_predicted_cpm / chk.ecm_observed_cpm.replace(0, np.nan)
    return clean, chk


# Gribben (symbol rows)
gh = pd.read_csv(os.path.join(D, "gse202379", "pb_hep.csv"), index_col=0)
ga = pd.read_csv(os.path.join(D, "gse202379", "pb_all.csv"), index_col=0)
g_clean, g_chk = decontaminate(gh, ga, pd.Series(gh.index, index=gh.index), "Gribben MASLD")
g_clean.to_csv(os.path.join(D, "gse202379", "pb_hep_decont.csv"))
# Andrews (ENSG rows)
ah = pd.read_csv(os.path.join(D, "pb_hep_ensg.csv"), index_col=0)
aa = pd.read_csv(os.path.join(D, "pb_whole_ensg.csv"), index_col=0)
a_clean, a_chk = decontaminate(ah, aa, ann.symbol, "Andrews healthy + PSC/PBC")
a_clean.to_csv(os.path.join(D, "pb_hepdecont_ensg.csv"))

CHK = pd.concat([g_chk, a_chk])
CHK.to_csv(os.path.join(R, "decont_ecm_check.csv"), index=False)
pt = pd.read_csv(os.path.join(D, "gse202379", "patients.csv"), dtype={"Patient.ID": str}).set_index("Patient.ID")
g_chk["group"] = pt.reindex(g_chk["sample"])["Disease.status"].values
pd.set_option("display.width", 200)
print(g_chk.groupby("group")[["alpha", "ecm_observed_cpm", "ecm_predicted_cpm", "ecm_explained_frac"]].median().round(3).to_string())
print(a_chk.assign(group=a_chk["sample"].str.split("_").str[0].str[:3])
      .groupby("group")[["alpha", "ecm_observed_cpm", "ecm_predicted_cpm", "ecm_explained_frac"]].median().round(3).to_string())
