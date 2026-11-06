#!/usr/bin/env python3
"""24_programme_in_bulk.py — test a hepatocyte programme defined in nuclei on an isolated-hepatocyte bulk cohort.

The programme (results/nuclei_hepatocyte_programme.csv) is defined from public nuclei data only: hepatocyte-made
genes (MacParland atlas hepatocyte fraction >= 0.5) that move the same way, t > 2.5, in both decontaminated
hepatocyte-nucleus contrasts (MASH fibrosis trend within one sample source; PSC/PBC vs healthy). The bulk cohort is
never used to choose genes. Per sample: score = mean per-gene z (log2 CPM) of the up genes minus that of the down
genes. Models (OLS) use the cohort's covariates: logit(MT fraction), age, sex, log10(assigned reads).
  A  disease vs normal donors          B  Child-Pugh C vs B (disease only)          C  MASH vs ALD (disease only)
Calibration: 1,000 random gene sets of hepatocyte-made genes, matched to the programme's up and down genes on
mean-expression decile, scored and modelled identically; empirical two-sided p from |t|.
The bulk cohort is private: inputs are read from a local directory and outputs written to a private directory.
usage: python 24_programme_in_bulk.py <dir: meta.csv, counts_noMT_noNUMT.csv.gz, cell_of_origin.csv> <programme csv> <out dir>
"""
import os, sys
import numpy as np, pandas as pd
import statsmodels.formula.api as smf

inp, prog_file, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
N_NULL, SEED = 1000, 20261105
cnt = pd.read_csv(os.path.join(inp, "counts_noMT_noNUMT.csv.gz"), index_col=0)
meta = pd.read_csv(os.path.join(inp, "meta.csv"), index_col=0).reindex(cnt.columns)
co = pd.read_csv(os.path.join(inp, "cell_of_origin.csv")).drop_duplicates("symbol").set_index("symbol")["hep_fraction"]
prog = pd.read_csv(prog_file, index_col=0)

ids = meta.index[(meta.role == "normal_donor") | ((meta.role == "esld") & (meta.esld_in_mt_range.astype(str) == "True"))]
m = meta.loc[ids].copy()
m["group"] = np.where(m.role == "normal_donor", "normal", "disease")
mt = (m.mt_pct / 100).clip(1e-4, 1 - 1e-4)
m["logit_mt"] = np.log(mt / (1 - mt))
m["log10_assigned"] = np.log10(m.assigned)
m["cp"] = m.child_pugh.astype(str)
X = cnt[ids]
X = X[(X >= 5).sum(axis=1) >= 3]
lc = np.log2(X / X.sum() * 1e6 + 1)
Z = lc.sub(lc.mean(axis=1), axis=0).div(lc.std(axis=1).replace(0, np.nan), axis=0)

up = [g for g in prog.index[prog.direction == "up"] if g in Z.index]
dn = [g for g in prog.index[prog.direction == "down"] if g in Z.index]


def score(u, d):
    return Z.loc[u].mean() - Z.loc[d].mean()


COV = "logit_mt + age + C(sex) + log10_assigned"
MODELS = {"A_disease_vs_normal": (m.index, "group", "disease", "normal"),
          "B_ChildPugh_C_vs_B": (m.index[(m.group == "disease") & m.cp.isin(["B", "C"])], "cp", "C", "B"),
          "C_MASH_vs_ALD": (m.index[(m.group == "disease") & m.etiology_group.isin(["MASH", "ALD"])], "etiology_group", "MASH", "ALD")}


def fit(s):
    res = {}
    for name, (rows, var, lvl, ref) in MODELS.items():
        d = m.loc[rows].assign(score=s[rows])
        f = smf.ols(f"score ~ C({var}, Treatment('{ref}')) + {COV}", data=d).fit()
        term = f"C({var}, Treatment('{ref}'))[T.{lvl}]"
        res[name] = (f.params[term], f.tvalues[term], f.pvalues[term], int((d[var] == lvl).sum()), int((d[var] == ref).sum()))
    return res


obs = fit(score(up, dn))
# expression-matched null from hepatocyte-made genes outside the programme
pool = pd.Index([g for g in Z.index if co.get(g, 0) >= 0.5 and g not in prog.index])
dec = pd.qcut(lc.mean(axis=1), 10, labels=False)
rng = np.random.default_rng(SEED)


def draw(genes):
    return [rng.choice(pool[dec.reindex(pool).to_numpy() == k]) for k in dec.reindex(genes).to_numpy()]


null = [fit(score(draw(up), draw(dn))) for _ in range(N_NULL)]
rows = []
for name in MODELS:
    est, t, p, n1, n0 = obs[name]
    nt = np.array([r[name][1] for r in null])
    rows.append(dict(model=name, n_level=n1, n_ref=n0, estimate=est, t=t, p_ols=p,
                     p_matched_null=float((np.abs(nt) >= abs(t)).mean()), null_t_median=float(np.median(nt)),
                     null_t_p95=float(np.percentile(np.abs(nt), 95))))
R = pd.DataFrame(rows)
R.to_csv(os.path.join(out, "programme_in_bulk.csv"), index=False)
# per-gene direction in Child-Pugh C vs B
d = m.loc[MODELS["B_ChildPugh_C_vs_B"][0]]
diff = lc.loc[up + dn, d.index[d.cp == "C"]].mean(axis=1) - lc.loc[up + dn, d.index[d.cp == "B"]].mean(axis=1)
agree = int((diff[up] > 0).sum() + (diff[dn] < 0).sum())
print(f"programme genes in bulk: up {len(up)}, down {len(dn)}; genes moving in the programme direction, CP C vs B: {agree}/{len(up) + len(dn)}")
pd.set_option("display.width", 200)
print(R.round(4).to_string(index=False))
