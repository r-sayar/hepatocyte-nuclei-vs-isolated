#!/usr/bin/env python3
"""16_concordance.py — hepatocyte nuclei (and cells) pseudobulk vs the isolated-hepatocyte bulk cohort.

Readouts asked for, at PATHWAY level: (a) Spearman correlation of the pathway ranks (NES) between them, (b) which
pathways turn from significant to not significant (and the reverse); and (f) the same for the reported lists only,
the top 25 pathways each way. The same at gene level and the expression-profile Spearman are kept as support. No matched samples, so every comparison is disease-vs-healthy within each
modality, and every flip count is set against the same bulk cut down to the pseudobulk's size
(15_bulk_profiles_and_downsample.R: 50 random 6-ESLD vs 4-donor draws, same model). A flip only says something
about the platform if it is more common than in the n-matched bulk draws.

Genes are joined on ENSG (bulk is keyed by GENCODE v44 gene name, mapped through unambiguous names), restricted
to protein-coding genes. Significance: limma BH < 0.05 ("sig"), and BH < 0.05 with |log2FC| > 1 ("DEG").
Pathways: fast_gsea (the validated numpy re-implementation of gseapy prerank; scripts/ref/common.py) on every
ranking, so full bulk, strata, draws and pseudobulks are scored identically; called at nominal p < 0.05 (the
project's figure gate) and at FDR < 0.25.
Gene flips are also split by cell of origin (MacParland-atlas classes), because isolated-hepatocyte preps
carry some non-parenchymal cells whose genes a hepatocyte-only pseudobulk should not reproduce.
Outputs: results/conc_*.csv
"""
import os, sys
from multiprocessing import Pool
import numpy as np, pandas as pd

ROOT = "/srv/public/sayar99/liver_strat_gsea"
sys.path.insert(0, os.path.join(ROOT, "scripts", "ref"))
import common as C  # noqa: E402

C.GMT_DIR = os.path.join(ROOT, "genesets")
D, R = os.path.join(ROOT, "data"), os.path.join(ROOT, "results")
MIN_HEP = 50
ann = pd.read_csv(os.path.join(D, "ensg2symbol_v46.tsv"), sep="\t", header=None, names=["ensg", "symbol", "type"])
ann = ann.drop_duplicates("ensg").set_index("ensg")
PC = set(ann.index[ann.type == "protein_coding"])
n2e = pd.read_csv(os.path.join(D, "v44_name2ensg.csv")).set_index("gene_name")["ensg"]
lib = pd.read_csv(os.path.join(D, "libraries.tsv"), sep="\t")


def by_ensg(d, key):
    """index -> ENSG (bulk: via v44 gene name), protein-coding only, one row per ENSG."""
    d = d.copy()
    d["ensg"] = n2e.reindex(d.index).to_numpy() if key == "name" else d.index
    d = d.dropna(subset=["ensg"]).drop_duplicates("ensg").set_index("ensg")
    return d[d.index.isin(PC)]


bulk = {s: by_ensg(pd.read_csv(os.path.join(D, f"bulk_{s}_limma.csv"), index_col=0), "name")
        for s in ("esld_all", "ALD", "MASH", "CP_B", "CP_C")}
PBK = ["sn_hep", "sc_hep", "sn_hepauthor", "sc_hepauthor", "sn_whole", "sc_whole"]
pb = {k: by_ensg(pd.read_csv(os.path.join(R, f"pb_{k}_limma.csv"), index_col=0), "ensg") for k in PBK}
ds = pd.read_csv(os.path.join(D, "bulk_downsample_t.csv.gz")).set_index("gene")
draws = {k: by_ensg(g.drop(columns="draw"), "name") for k, g in ds.groupby("draw")}
FULL = bulk["esld_all"]

# ---------------- (a) expression Spearman ---------------------------------------------------------
bm = by_ensg(pd.read_csv(os.path.join(D, "bulk_group_logcpm.csv"), index_col=0), "name")
cnt = {"hep": pd.read_csv(os.path.join(D, "pb_hep_ensg.csv"), index_col=0),
       "whole": pd.read_csv(os.path.join(D, "pb_whole_ensg.csv"), index_col=0)}
auth = {m: pd.read_csv(os.path.join(D, f"pb_{m}_hepauthor_counts.csv.gz"), index_col=0) for m in ("sc", "sn")}
genes = sorted(set(bm.index) & set(cnt["hep"].index) & PC)


def logcpm(v):
    v = v.reindex(genes).fillna(0)
    return np.log2(v / v.sum() * 1e6 + 1)


profiles = {}
for mat in ("hep", "whole"):
    for mod in ("sc", "sn"):
        for grp, key in (("healthy", "healthy"), ("PSC_PBC", "disease")):
            L = lib[(lib.modality == mod) & (lib.group == grp)]
            if mat == "hep":
                L = L[L.n_hep >= MIN_HEP]
            profiles[f"{mod}_{mat}_{key}"] = logcpm(cnt[mat][L.library].sum(axis=1))
for mod in ("sc", "sn"):
    hcols = [c for c in auth[mod].columns if c.startswith("Healthy")]
    profiles[f"{mod}_hepauthor_healthy"] = logcpm(auth[mod][hcols].sum(axis=1))
profiles["bulk_normal"] = bm.reindex(genes)["normal_mean_log2cpm"]
profiles["bulk_esld"] = bm.reindex(genes)["esld_mean_log2cpm"]
halves = by_ensg(pd.read_csv(os.path.join(D, "bulk_normal_halves_logcpm.csv"), index_col=0), "name").reindex(genes)
profiles["bulk_normal_half1"], profiles["bulk_normal_half2"] = halves.iloc[:, 0], halves.iloc[:, 1]
P = pd.DataFrame(profiles)
expr_bulk = 2 ** P["bulk_normal"] - 1 >= 1
rows = []
pairs = [("sn_hep_healthy", "bulk_normal"), ("sn_hepauthor_healthy", "bulk_normal"), ("sc_hep_healthy", "bulk_normal"),
         ("sc_hepauthor_healthy", "bulk_normal"), ("sn_whole_healthy", "bulk_normal"), ("sc_whole_healthy", "bulk_normal"),
         ("sn_hep_disease", "bulk_esld"), ("sc_hep_disease", "bulk_esld"), ("sn_whole_disease", "bulk_esld"),
         ("sc_whole_disease", "bulk_esld"),
         ("bulk_normal_half1", "bulk_normal_half2"), ("bulk_normal", "bulk_esld"),
         ("sn_hep_healthy", "sc_hep_healthy"), ("sn_hep_healthy", "sn_hep_disease")]
for a, b in pairs:
    both = (2 ** P[a] - 1 >= 1) & (2 ** P[b] - 1 >= 1)
    rows.append(dict(a=a, b=b, rho_all_pc=P[a].corr(P[b], method="spearman"), n_all=len(P),
                     rho_expressed_both=P.loc[both, a].corr(P.loc[both, b], method="spearman"), n_expressed=int(both.sum()),
                     pearson_log_expressed=P.loc[both, a].corr(P.loc[both, b])))
EXPR = pd.DataFrame(rows)
per_lib = []
for mat in ("hep", "whole"):
    for _, l in lib.iterrows():
        if mat == "hep" and l.n_hep < MIN_HEP:
            continue
        v = logcpm(cnt[mat][l.library])
        ref = P["bulk_normal"] if l.group == "healthy" else P["bulk_esld"]
        per_lib.append(dict(library=l.library, material=mat, modality=l.modality, group=l.group,
                            rho_all_pc=v.corr(ref, method="spearman"),
                            rho_bulk_expressed=v[expr_bulk].corr(ref[expr_bulk], method="spearman")))
PERLIB = pd.DataFrame(per_lib)

# ---------------- (b) DE t Spearman ---------------------------------------------------------------
def t_rho(x, y):
    g = x.index.intersection(y.index)
    return float(x.loc[g, "t"].corr(y.loc[g, "t"], method="spearman")), len(g)


de_rows = []
for s, B in bulk.items():
    for k, Q in pb.items():
        r, n = t_rho(B, Q)
        de_rows.append(dict(bulk=s, alt=k, rho_t=r, n_genes=n))
ceiling = [t_rho(FULL, d)[0] for d in draws.values()]
draw_vs_sn = [t_rho(d, pb["sn_hep"])[0] for d in draws.values()]
de_rows += [dict(bulk="esld_all", alt="bulk 6v4 draws (median)", rho_t=float(np.median(ceiling)), n_genes=np.nan,
                 lo5=float(np.percentile(ceiling, 5)), hi95=float(np.percentile(ceiling, 95))),
            dict(bulk="bulk 6v4 draws (median)", alt="sn_hep", rho_t=float(np.median(draw_vs_sn)), n_genes=np.nan,
                 lo5=float(np.percentile(draw_vs_sn, 5)), hi95=float(np.percentile(draw_vs_sn, 95)))]
DE = pd.DataFrame(de_rows)

# ---------------- (c) gene-level flips ------------------------------------------------------------
def calls(d, deg):
    s = d["adj.P.Val"] < 0.05
    if deg:
        s &= d["logFC"].abs() > 1
    return pd.Series(np.where(s, np.sign(d["t"]), 0), index=d.index)


def fate(ref, alt, deg):
    """genes called in ref -> same / not significant / opposite in alt (fractions of those tested in alt)."""
    cr, ca = calls(ref, deg), calls(alt, deg)
    hit = cr[cr != 0]
    tested = hit.index.intersection(ca.index)
    a = ca.loc[tested]
    same, opp = int((a == hit.loc[tested]).sum()), int(((a != 0) & (a != hit.loc[tested])).sum())
    ns = len(tested) - same - opp
    return dict(n_called_ref=len(hit), n_tested_alt=len(tested), same=same, not_sig=ns, opposite=opp,
                frac_same=same / max(len(tested), 1), frac_not_sig=ns / max(len(tested), 1),
                frac_opposite=opp / max(len(tested), 1))


flip_rows = []
for deg in (False, True):
    tag = "DEG (BH<0.05 & |log2FC|>1)" if deg else "sig (BH<0.05)"
    for k in PBK:
        flip_rows.append(dict(call=tag, direction="bulk -> pseudobulk", ref="bulk esld_all", alt=k, **fate(FULL, pb[k], deg)))
        flip_rows.append(dict(call=tag, direction="pseudobulk -> bulk", ref=k, alt="bulk esld_all", **fate(pb[k], FULL, deg)))
    fw = pd.DataFrame([fate(FULL, d, deg) for d in draws.values()])
    rv = pd.DataFrame([fate(d, FULL, deg) for d in draws.values()])
    for direction, X in (("bulk -> pseudobulk", fw), ("pseudobulk -> bulk", rv)):
        med = X.median(numeric_only=True)
        flip_rows.append(dict(call=tag, direction=direction, ref="bulk esld_all" if direction.startswith("bulk") else "bulk 6v4 draw",
                              alt="bulk 6v4 draw" if direction.startswith("bulk") else "bulk esld_all",
                              **{c: float(med[c]) for c in X.columns},
                              frac_not_sig_p5=float(X.frac_not_sig.quantile(.05)), frac_not_sig_p95=float(X.frac_not_sig.quantile(.95)),
                              frac_opposite_p95=float(X.frac_opposite.quantile(.95))))
FLIP = pd.DataFrame(flip_rows)

# the same gene flips split by cell of origin: isolated-hepatocyte preps carry some non-parenchymal cells, so a
# bulk "disease" gene can be carry-over that a hepatocyte-only pseudobulk rightly lacks. Classes are the project's
# (MacParland healthy-liver atlas hepatocyte fraction: hepatocyte >= 0.5, mixed 0.2-0.5, non-parenchymal < 0.2)
coo = pd.read_csv(os.path.join(D, "cell_of_origin.csv")).dropna(subset=["symbol"]).set_index("symbol")
ORIGIN = pd.Series(coo["cell_of_origin"].to_numpy(), index=n2e.reindex(coo.index).to_numpy())
ORIGIN = ORIGIN[pd.notna(ORIGIN.index) & ~ORIGIN.index.duplicated()]
orow = []
for cls in ("hepatocyte", "mixed", "non-parenchymal"):
    ids = ORIGIN.index[ORIGIN == cls]
    sub = lambda d: d[d.index.isin(ids)]
    for k in ("sn_hep", "sn_hepauthor", "sc_hep"):
        orow.append(dict(origin=cls, direction="bulk -> pseudobulk", alt=k, **fate(sub(FULL), sub(pb[k]), False)))
        orow.append(dict(origin=cls, direction="pseudobulk -> bulk", alt=k, **fate(sub(pb[k]), sub(FULL), False)))
    for direction, X in (("bulk -> pseudobulk", pd.DataFrame([fate(sub(FULL), sub(d), False) for d in draws.values()])),
                         ("pseudobulk -> bulk", pd.DataFrame([fate(sub(d), sub(FULL), False) for d in draws.values()]))):
        orow.append(dict(origin=cls, direction=direction, alt="bulk 6v4 draws (median)",
                         **{c: float(X[c].median()) for c in X.columns},
                         frac_not_sig_p5=float(X.frac_not_sig.quantile(.05)), frac_not_sig_p95=float(X.frac_not_sig.quantile(.95)),
                         frac_opposite_p95=float(X.frac_opposite.quantile(.95))))
FLIP_ORIGIN = pd.DataFrame(orow)

# named genes: robust in small bulk, lost in sn hepatocytes; and sn-only calls the far larger bulk does not support.
# No gene is significant in >= 80% of the 6v4 draws (max 72%), so "robust at this n" = significant, same direction,
# in at least half of them.
ROBUST = 0.5
cf = calls(FULL, False)
dcalls = pd.DataFrame({k: calls(d, False).reindex(FULL.index) for k, d in draws.items()})
keep_frac = (dcalls.eq(cf, axis=0) & (cf != 0).to_numpy()[:, None]).mean(axis=1)
sn = pb["sn_hep"]
csn = calls(sn, False)
G = pd.DataFrame({"symbol": ann.reindex(FULL.index).symbol, "bulk_logFC": FULL.logFC, "bulk_t": FULL.t,
                  "bulk_adjP": FULL["adj.P.Val"], "draws_same_sig_frac": keep_frac,
                  "sn_hep_logFC": sn.logFC.reindex(FULL.index), "sn_hep_t": sn.t.reindex(FULL.index),
                  "sn_hep_adjP": sn["adj.P.Val"].reindex(FULL.index), "sn_call": csn.reindex(FULL.index)})
robust = G[(cf != 0) & (G.draws_same_sig_frac >= ROBUST)]
robust.assign(origin=ORIGIN.reindex(robust.index).fillna("unknown").to_numpy(), sn_status=np.where(robust.sn_call.isna(), "not tested", np.where(robust.sn_call == np.sign(robust.bulk_t), "same",
              np.where(robust.sn_call == 0, "not sig", "opposite")))).to_csv(os.path.join(R, "conc_genes_robust_in_small_bulk.csv"))
lost = G[(cf != 0) & (G.draws_same_sig_frac >= ROBUST) & (G.sn_call.fillna(0) != np.sign(G.bulk_t))]
lost = lost.assign(origin=ORIGIN.reindex(lost.index).fillna("unknown").to_numpy(), sn_status=np.where(lost.sn_call.isna(), "not tested", np.where(lost.sn_call == 0, "not sig", "opposite")))
gained = pd.DataFrame({"symbol": ann.reindex(sn.index).symbol, "sn_hep_logFC": sn.logFC, "sn_hep_t": sn.t,
                       "sn_hep_adjP": sn["adj.P.Val"], "bulk_logFC": FULL.logFC.reindex(sn.index),
                       "bulk_t": FULL.t.reindex(sn.index), "bulk_adjP": FULL["adj.P.Val"].reindex(sn.index),
                       "bulk_call": cf.reindex(sn.index)})
gained = gained[(csn != 0) & (gained.bulk_call.fillna(0) != np.sign(gained.sn_hep_t))]
gained = gained.assign(origin=ORIGIN.reindex(gained.index).fillna("unknown").to_numpy(), bulk_status=np.where(gained.bulk_call.isna(), "not tested", np.where(gained.bulk_call == 0, "not sig", "opposite")))

# ---------------- (d) pathways: significant -> not significant, both directions --------------------
# Every ranking is scored with fast_gsea on v46 symbols so full bulk, strata, draws and pseudobulks are comparable.
GATES = {"p<0.05": ("p", 0.05), "FDR<0.25": ("fdr", 0.25)}  # project figure gate; the usual GSEA FDR gate


def rnk(d):
    s = d["t"].copy()
    s.index = ann.reindex(s.index).symbol.to_numpy()
    s = s[pd.notna(s.index)]
    return s[~s.index.duplicated()].sort_values(ascending=False)


def fg(d):
    return C.fast_gsea(rnk(d))


with Pool(8) as pool:
    dtabs = pool.map(fg, list(draws.values()))
ftab = fg(FULL)
stab = {s: fg(bulk[s]) for s in bulk if s != "esld_all"}
ptab = {k: fg(pb[k]) for k in ("sn_hep", "sc_hep", "sn_hepauthor", "sn_whole", "sc_whole")}


def pcalls(t, gate):
    col, thr = GATES[gate]
    return pd.Series(np.where(t[col] < thr, np.sign(t.NES), 0), index=t.index)


def pfate(ref_c, alt_c):
    """pathways called in ref -> same / not significant / opposite in alt."""
    hit = ref_c[ref_c != 0]
    a = alt_c.reindex(hit.index).dropna()
    same, opp = int((a == hit.loc[a.index]).sum()), int(((a != 0) & (a != hit.loc[a.index])).sum())
    n = max(len(a), 1)
    return dict(n_called=len(hit), n_tested=len(a), same=same, not_sig=len(a) - same - opp, opposite=opp,
                frac_lost=(len(a) - same) / n, frac_not_sig=(len(a) - same - opp) / n, frac_opposite=opp / n)


prow, pstat = [], []
for gate in GATES:
    for lib in C.LIB_FILES:
        fc = pcalls(ftab[lib], gate)
        dc = [pcalls(t[lib], gate) for t in dtabs]
        for k, t in ptab.items():
            pc = pcalls(t[lib], gate)
            prow.append(dict(gate=gate, lib=lib, direction="bulk -> pseudobulk", other=k, **pfate(fc, pc)))
            prow.append(dict(gate=gate, lib=lib, direction="pseudobulk -> bulk", other=k, **pfate(pc, fc)))
        for direction, X in (("bulk -> pseudobulk", pd.DataFrame([pfate(fc, d) for d in dc])),
                             ("pseudobulk -> bulk", pd.DataFrame([pfate(d, fc) for d in dc]))):
            prow.append(dict(gate=gate, lib=lib, direction=direction, other="bulk 6v4 draws (median)",
                             **{c: float(X[c].median()) for c in X.columns},
                             frac_lost_p5=float(X.frac_lost.quantile(.05)), frac_lost_p95=float(X.frac_lost.quantile(.95))))
        # named: every pathway called by the full bulk or by hepatocyte nuclei
        Dm = pd.concat(dc, axis=1)
        nt, ct = ptab["sn_hep"][lib], ptab["sc_hep"][lib]
        nc, cc = pcalls(nt, gate), pcalls(ct, gate)
        for term in fc.index[fc != 0].union(nc.index[nc != 0]):
            b, n = fc.get(term, np.nan), nc.get(term, np.nan)
            status = ("not tested in bulk" if pd.isna(b) else "not tested in nuclei" if pd.isna(n) else
                      "same" if b == n else "opposite" if (b != 0 and n != 0) else "bulk only" if b != 0 else "nuclei only")
            want = b if (pd.notna(b) and b != 0) else n
            pstat.append(dict(gate=gate, lib=lib, term=term, status=status,
                              bulk_NES=ftab[lib].NES.get(term, np.nan), bulk_p=ftab[lib].p.get(term, np.nan),
                              bulk_fdr=ftab[lib].fdr.get(term, np.nan),
                              nuclei_NES=nt.NES.get(term, np.nan), nuclei_p=nt.p.get(term, np.nan), nuclei_fdr=nt.fdr.get(term, np.nan),
                              cells_NES=ct.NES.get(term, np.nan), cells_p=ct.p.get(term, np.nan),
                              cells_call=cc.get(term, np.nan),
                              small_bulk_calls_it=float((Dm.loc[term] == want).mean()) if term in Dm.index else np.nan))
PFLIP, PSTAT = pd.DataFrame(prow), pd.DataFrame(pstat)

# ---------------- (e) pathway-rank Spearman ---------------------------------------------------------
# The rank of the pathways, not of the genes: Spearman of NES between two disease-vs-healthy GSEA results, per
# library, over all pathways scored in both and over the pathways the full bulk calls (p < 0.05).
TABS = {"bulk esld_all": ftab, **{f"bulk {s}": t for s, t in stab.items()}, **ptab}
RPAIRS = ([("bulk esld_all", k) for k in ptab] + [(f"bulk {s}", "sn_hep") for s in stab]
          + [("bulk ALD", "bulk MASH"), ("bulk CP_B", "bulk CP_C"), ("sc_hep", "sn_hep"), ("sc_whole", "sn_whole")])


def nes_rho(x, y, terms=None):
    g = x.index.intersection(y.index)
    if terms is not None:
        g = g.intersection(terms)
    return float(x.NES[g].corr(y.NES[g], method="spearman")), len(g)


rrow = []
for lib in C.LIB_FILES:
    called = ftab[lib].index[ftab[lib].p < 0.05]
    for a, b in RPAIRS:
        r_all, n_all = nes_rho(TABS[a][lib], TABS[b][lib])
        r_c, n_c = nes_rho(TABS[a][lib], TABS[b][lib], called)
        rrow.append(dict(lib=lib, a=a, b=b, rho_all=r_all, n_all=n_all, rho_bulk_called=r_c, n_bulk_called=n_c))
    for a, b, vals in (("bulk esld_all", "bulk 6v4 draws", [(nes_rho(ftab[lib], t[lib]), nes_rho(ftab[lib], t[lib], called)) for t in dtabs]),
                       ("bulk 6v4 draws", "sn_hep", [(nes_rho(t[lib], ptab["sn_hep"][lib]), nes_rho(t[lib], ptab["sn_hep"][lib], called)) for t in dtabs])):
        ra = np.array([v[0][0] for v in vals]); rc = np.array([v[1][0] for v in vals])
        rrow.append(dict(lib=lib, a=a, b=b, rho_all=float(np.median(ra)), rho_all_p5=float(np.percentile(ra, 5)),
                         rho_all_p95=float(np.percentile(ra, 95)), n_all=vals[0][0][1], rho_bulk_called=float(np.median(rc)),
                         rho_bulk_called_p5=float(np.percentile(rc, 5)), rho_bulk_called_p95=float(np.percentile(rc, 95)),
                         n_bulk_called=vals[0][1][1]))
PRANK = pd.DataFrame(rrow)

# ---------------- (f) the reported lists: top 25 pathways each way ----------------------------------
# The project's figure rule (nominal p < 0.05, then the top N each way by NES) with N = 25. For each library and
# direction: how many of the bulk's top 25 are in the other's top 25, how many stay significant in the same
# direction, and Spearman of NES over the bulk's top 50 (25 up + 25 down); the 6v4 draws are the yardstick.
TOPN = 25
DIRS = (("up", "top15_up", 1), ("down", "bottom15_down", -1))  # keys as returned by C.top_lists


def top_compare(ref, alt, which, sign):
    lr, la = C.top_lists(ref, n=TOPN)[which], C.top_lists(alt, n=TOPN)[which]
    a = alt.reindex(lr)
    still = int(((np.sign(a.NES) == sign) & (a.p < 0.05)).sum())
    opp = int(((np.sign(a.NES) == -sign) & (a.p < 0.05)).sum())
    return dict(n_bulk_top=len(lr), n_other_top=len(la), overlap=len(set(lr) & set(la)), still_sig=still,
                opposite=opp, not_sig=int(a.NES.notna().sum()) - still - opp, not_tested=int(a.NES.isna().sum()))


def top50_rho(ref, alt):
    tl = C.top_lists(ref, n=TOPN)
    g = pd.Index(tl["top15_up"] + tl["bottom15_down"]).intersection(alt.index)
    return float(ref.NES[g].corr(alt.NES[g], method="spearman"))


trow, tnamed = [], []
for lib in C.LIB_FILES:
    ref = ftab[lib]
    for dname, which, sign in DIRS:
        for k in ("sn_hep", "sn_hepauthor", "sc_hep"):
            trow.append(dict(lib=lib, direction=dname, other=k, **top_compare(ref, ptab[k][lib], which, sign)))
        X = pd.DataFrame([top_compare(ref, t[lib], which, sign) for t in dtabs])
        trow.append(dict(lib=lib, direction=dname, other="bulk 6v4 draws (median)", **{c: float(X[c].median()) for c in X.columns},
                         overlap_p5=float(X.overlap.quantile(.05)), overlap_p95=float(X.overlap.quantile(.95)),
                         still_sig_p5=float(X.still_sig.quantile(.05)), still_sig_p95=float(X.still_sig.quantile(.95))))
        # named: the bulk's top 25 in this direction, with the nuclei rank / status and how stable it is in small bulk
        lr = C.top_lists(ref, n=TOPN)[which]
        nt = ptab["sn_hep"][lib]
        nrank = (sign * nt.NES).rank(ascending=False)
        in_top = pd.DataFrame({i: pd.Series(1, index=C.top_lists(t[lib], n=TOPN)[which]).reindex(lr).fillna(0) for i, t in enumerate(dtabs)}).mean(axis=1)
        for i, term in enumerate(lr, 1):
            n_nes, n_p = nt.NES.get(term, np.nan), nt.p.get(term, np.nan)
            status = ("not tested" if pd.isna(n_nes) else "significant, same" if (np.sign(n_nes) == sign and n_p < 0.05)
                      else "significant, opposite" if n_p < 0.05 else "not significant")
            tnamed.append(dict(lib=lib, direction=dname, bulk_rank=i, term=term, bulk_NES=ref.at[term, "NES"], bulk_p=ref.at[term, "p"],
                               nuclei_NES=n_nes, nuclei_p=n_p, nuclei_rank_same_direction=nrank.get(term, np.nan),
                               nuclei_status=status, in_nuclei_top25=term in C.top_lists(nt, n=TOPN)[which],
                               cells_NES=ptab["sc_hep"][lib].NES.get(term, np.nan), cells_p=ptab["sc_hep"][lib].p.get(term, np.nan),
                               small_bulk_top25_frac=float(in_top[term])))
    rv = [top50_rho(ref, t[lib]) for t in dtabs]
    trow.append(dict(lib=lib, direction="top 50 (25 up + 25 down)", other="bulk 6v4 draws (median)", rho_top50=float(np.median(rv)),
                     rho_top50_p5=float(np.percentile(rv, 5)), rho_top50_p95=float(np.percentile(rv, 95))))
    for k in ("sn_hep", "sn_hepauthor", "sc_hep"):
        trow.append(dict(lib=lib, direction="top 50 (25 up + 25 down)", other=k, rho_top50=top50_rho(ref, ptab[k][lib])))
TOP, TOPNAMED = pd.DataFrame(trow), pd.DataFrame(tnamed)
NESDIR =os.path.join(R, "pathway_nes")  # per-ranking fast_gsea tables, used by the figures
os.makedirs(NESDIR, exist_ok=True)
for k, t in {"bulk_esld_all": ftab, **{f"bulk_{s}": v for s, v in stab.items()}, **{f"pb_{k}": v for k, v in ptab.items()}}.items():
    for lib, tt in t.items():
        tt.to_csv(os.path.join(NESDIR, f"{k}_{lib}.csv"))

for name, df in (("conc_expression_spearman", EXPR), ("conc_expression_spearman_per_library", PERLIB),
                 ("conc_de_t_spearman", DE), ("conc_gene_flips", FLIP), ("conc_gene_flips_by_origin", FLIP_ORIGIN),
                 ("conc_genes_lost_in_sn_hep", lost), ("conc_genes_sn_hep_only", gained),
                 ("conc_pathway_flips", PFLIP), ("conc_pathway_status", PSTAT), ("conc_pathway_rank_spearman", PRANK),
                 ("conc_top25_summary", TOP), ("conc_top25_pathways", TOPNAMED)):
    df.to_csv(os.path.join(R, f"{name}.csv"), index=name in ("conc_genes_lost_in_sn_hep", "conc_genes_sn_hep_only"))

pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 45)
print(EXPR.round(3).to_string(index=False))
print(DE.round(3).to_string(index=False))
print(FLIP.round(3).to_string(index=False))
print(FLIP_ORIGIN.round(3).to_string(index=False))
print(PRANK.round(3).to_string(index=False))
print(TOP.round(2).to_string(index=False))
print(PFLIP.round(3).to_string(index=False))
print(PSTAT.groupby(["gate", "lib", "status"]).size().unstack(fill_value=0).to_string())
