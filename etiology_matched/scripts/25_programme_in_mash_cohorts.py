#!/usr/bin/env python3
"""25_programme_in_mash_cohorts.py — does the nuclei-defined hepatocyte programme track fibrosis stage in public MASLD
liver-biopsy RNA-seq cohorts?

Programme: results/nuclei_hepatocyte_programme.csv (hepatocyte-made genes moving the same way in decontaminated
hepatocyte nuclei of MASH and PSC/PBC livers; no bulk data used to choose them). Cohorts (whole-biopsy bulk):
  GSE162694  normal histology + NAFLD F0-F4 (raw counts, Ensembl IDs)
  GSE130970  NAFLD F0-F4 (Hoang et al. 2019; salmon/tximport counts, Entrez IDs)
Per sample: score = mean z (log2 CPM) of the up genes minus the down genes. Because biopsies lose hepatocytes as
fibrosis advances, every model adjusts for hepatocyte content (mean z of ALB, APOA1, APOC3, APOH, FGB, HP, SERPINA1,
TF, TTR, RBP4): score ~ fibrosis stage (0-4, numeric) + hepatocyte content + age + sex, NAFLD samples only.
Calibration: 1,000 random sets of hepatocyte-made genes (MacParland atlas >= 0.5) matched on expression decile.
Per gene: the same model with that gene as the outcome; "replicated" = t > 2 in the programme direction in both.
Output: results/programme_in_mash_cohorts.csv, results/programme_by_stage.csv, results/programme_genes_in_mash_cohorts.csv
"""
import gzip, os, re
import numpy as np, pandas as pd
import statsmodels.formula.api as smf

ROOT = "/srv/public/sayar99/liver_strat_gsea"
D, R = os.path.join(ROOT, "data", "mash_bulk"), os.path.join(ROOT, "results")
HEP = ["ALB", "APOA1", "APOC3", "APOH", "FGB", "HP", "SERPINA1", "TF", "TTR", "RBP4"]
N_NULL, SEED = 1000, 20261106
prog = pd.read_csv(os.path.join(R, "nuclei_hepatocyte_programme.csv"), index_col=0)
co = pd.read_csv(os.path.join(ROOT, "data", "cell_of_origin.csv")).drop_duplicates("symbol").set_index("symbol")["hep_fraction"]


def series_matrix(path):
    """sample title + characteristics_ch1 'key: value' fields -> DataFrame indexed by GSM order."""
    rows, titles = {}, None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                key = vals[0].split(":")[0].strip().lower()
                rows[key] = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
    return pd.DataFrame(rows, index=titles)


def load_162694():
    c = pd.read_csv(os.path.join(D, "GSE162694_raw_counts.csv.gz"), index_col=0)
    ann = pd.read_csv(os.path.join(ROOT, "data", "ensg2symbol_v46.tsv"), sep="\t", header=None, names=["ensg", "symbol", "type"]).drop_duplicates("ensg").set_index("ensg")
    c.index = ann.symbol.reindex(c.index.str.split(".").str[0]).to_numpy()
    c = c[pd.notna(c.index)].groupby(level=0).sum()
    m = series_matrix(os.path.join(D, "GSE162694_series_matrix.txt.gz"))
    m.index = [t.split()[-1] for t in m.index]
    m["fibrosis"] = pd.to_numeric(m["fibrosis stage"], errors="coerce")  # "normal liver histology" -> NaN
    m["normal"] = m["fibrosis stage"].str.contains("normal")
    return c, m.rename(columns={"sex": "sex"})


def load_130970():
    c = pd.read_csv(os.path.join(D, "GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz"), index_col=0)
    h = pd.read_csv(os.path.join(D, "hgnc_complete_set.txt"), sep="\t", usecols=["symbol", "entrez_id"], low_memory=False).dropna()
    h = h.assign(entrez_id=h.entrez_id.astype(int)).drop_duplicates("entrez_id").set_index("entrez_id").symbol
    c.index = h.reindex(c.index.astype(int)).to_numpy()
    c = c[pd.notna(c.index)].groupby(level=0).sum().round()
    m = series_matrix(os.path.join(D, "GSE130970_series_matrix.txt.gz"))
    m["fibrosis"] = pd.to_numeric(m["fibrosis stage"], errors="coerce")
    m["normal"] = False
    m["age"] = m["age at biopsy"]
    return c, m


rng = np.random.default_rng(SEED)
rows, by_stage, per_gene = [], [], {}
for name, loader in (("GSE162694", load_162694), ("GSE130970", load_130970)):
    c, m = loader()
    m = m.reindex(c.columns)
    c = c[(c >= 5).sum(axis=1) >= 0.25 * c.shape[1]]
    lc = np.log2(c / c.sum() * 1e6 + 1)
    Z = lc.sub(lc.mean(axis=1), axis=0).div(lc.std(axis=1).replace(0, np.nan), axis=0)
    up = [g for g in prog.index[prog.direction == "up"] if g in Z.index]
    dn = [g for g in prog.index[prog.direction == "down"] if g in Z.index]
    d = m.assign(hep=Z.reindex(HEP).mean(), age=pd.to_numeric(m["age"], errors="coerce"),
                 sex=m["sex"].str[0].str.upper())
    nafld = d.index[d.fibrosis.notna()]

    def fit(u, dd):
        s = Z.loc[u].mean() - Z.loc[dd].mean()
        f = smf.ols("score ~ fibrosis + hep + age + C(sex)", data=d.loc[nafld].assign(score=s[nafld])).fit()
        return f.params["fibrosis"], f.tvalues["fibrosis"], f.pvalues["fibrosis"], s

    est, t, p, s = fit(up, dn)
    pool = pd.Index([g for g in Z.index if co.get(g, 0) >= 0.5 and g not in prog.index])
    dec = pd.qcut(lc.mean(axis=1), 10, labels=False)
    draw = lambda genes: [rng.choice(pool[dec.reindex(pool).to_numpy() == k]) for k in dec.reindex(genes).to_numpy()]
    nt = np.array([fit(draw(up), draw(dn))[1] for _ in range(N_NULL)])
    # per gene: same model with one gene's z as the outcome
    per_gene[name] = pd.Series({g: smf.ols("y ~ fibrosis + hep + age + C(sex)", data=d.loc[nafld].assign(y=Z.loc[g, nafld])).fit().tvalues["fibrosis"]
                                for g in up + dn})
    rows.append(dict(cohort=name, n_nafld=len(nafld), n_normal=int(d.normal.sum()), up_genes=len(up), down_genes=len(dn),
                     slope_per_stage=est, t=t, p_ols=p, p_matched_null=float((np.abs(nt) >= abs(t)).mean()),
                     null_abs_t_p95=float(np.percentile(np.abs(nt), 95))))
    st = d.assign(score=s).assign(stage=lambda x: np.where(x.normal, "normal", x.fibrosis.astype("Int64").astype(str)))
    by_stage.append(st.groupby("stage").score.agg(["mean", "count"]).assign(cohort=name))
R1 = pd.DataFrame(rows)
R1.to_csv(os.path.join(R, "programme_in_mash_cohorts.csv"), index=False)
B = pd.concat(by_stage)
B.to_csv(os.path.join(R, "programme_by_stage.csv"))
G = pd.DataFrame(per_gene).join(prog[["direction", "hep_frac"]])
sgn = np.where(G.direction == "up", 1, -1)
G["replicated_both_cohorts"] = (G.GSE162694 * sgn > 2) & (G.GSE130970 * sgn > 2)
G.round(3).to_csv(os.path.join(R, "programme_genes_in_mash_cohorts.csv"))
pd.set_option("display.width", 200)
print(R1.round(4).to_string(index=False))
print(B.round(3).to_string())
print(G.groupby("direction").replicated_both_cohorts.agg(["sum", "count"]).to_string())
