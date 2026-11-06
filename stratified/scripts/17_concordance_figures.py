#!/usr/bin/env python3
"""17_concordance_figures.py — headline figures: hepatocyte nuclei pseudobulk vs isolated-hepatocyte bulk.

figC_expression_spearman.png  left: healthy-liver hepatocyte nuclei pseudobulk vs normal-donor isolated hepatocytes
                              (mean log2 CPM, protein-coding genes); right: Spearman rho of every profile pair
                              against the same-platform references
figD_flips.png                significant -> same / not significant / opposite, both directions, each pseudobulk set
                              against the bulk cut to the same size (median of 50 6-vs-4 draws); all genes and
                              hepatocyte-origin genes
Style as 14_figures.py (reference diverging pair, neutral midpoint, text in ink colours).
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/srv/public/sayar99/liver_strat_gsea"
D, R = os.path.join(ROOT, "data"), os.path.join(ROOT, "results")
FIG = os.path.join(R, "figures")
os.makedirs(FIG, exist_ok=True)
INK, MUTED, GRID, SURF = "#1f1f1e", "#6b6b67", "#d9d8d4", "#fcfcfb"
BLUE, RED, NEUTRAL = "#2a78d6", "#e34948", "#c9c8c3"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": GRID,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": INK,
                     "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF})

# ---- fig C ------------------------------------------------------------------------------------------
ann = pd.read_csv(os.path.join(D, "ensg2symbol_v46.tsv"), sep="\t", header=None, names=["ensg", "symbol", "type"]).drop_duplicates("ensg").set_index("ensg")
PC = set(ann.index[ann.type == "protein_coding"])
n2e = pd.read_csv(os.path.join(D, "v44_name2ensg.csv")).set_index("gene_name")["ensg"]
bm = pd.read_csv(os.path.join(D, "bulk_group_logcpm.csv"), index_col=0)
bm.index = n2e.reindex(bm.index).to_numpy()
bm = bm[pd.notna(bm.index) & ~bm.index.duplicated()]
lib = pd.read_csv(os.path.join(D, "libraries.tsv"), sep="\t")
hep = pd.read_csv(os.path.join(D, "pb_hep_ensg.csv"), index_col=0)
genes = sorted(set(bm.index) & set(hep.index) & PC)
v = hep[lib[(lib.modality == "sn") & (lib.group == "healthy")].library].sum(axis=1).reindex(genes)
x = np.log2(v / v.sum() * 1e6 + 1)
y = bm.loc[genes, "normal_mean_log2cpm"]
E = pd.read_csv(os.path.join(R, "conc_expression_spearman.csv"))

fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.6, 5.2), gridspec_kw={"width_ratios": [1, 1.0]})
hb = a1.hexbin(y, x, gridsize=70, bins="log", cmap="Blues", mincnt=1, linewidths=0)
lim = (min(x.min(), y.min()) - .3, max(x.max(), y.max()) + .3)
a1.plot(lim, lim, color=GRID, lw=1)
a1.set_xlim(lim); a1.set_ylim(lim)
a1.set_xlabel("Isolated hepatocytes, bulk (normal donors), mean log2 CPM")
a1.set_ylabel("Hepatocyte nuclei pseudobulk (4 healthy livers), log2 CPM")
r = E[(E.a == "sn_hep_healthy") & (E.b == "bulk_normal")].iloc[0]
a1.text(0.03, 0.97, f"Spearman ρ = {r.rho_all_pc:.2f}, all {len(genes):,} protein-coding genes\n"
        f"ρ = {r.rho_expressed_both:.2f}, the {int(r.n_expressed):,} expressed in both (CPM ≥ 1)",
        transform=a1.transAxes, va="top", fontsize=8, color=MUTED)
a1.set_title("Healthy liver: hepatocyte nuclei vs isolated hepatocytes", loc="left", fontsize=10, color=INK)
for s in ("top", "right"):
    a1.spines[s].set_visible(False)

rows = [("bulk_normal_half1", "bulk_normal_half2", "Bulk vs bulk: two halves of the normal donors"),
        ("bulk_normal", "bulk_esld", "Bulk vs bulk: normal donors vs ESLD"),
        ("sc_hep_healthy", "bulk_normal", "Hepatocyte CELLS (scRNA) vs bulk, healthy"),
        ("sc_hep_disease", "bulk_esld", "Hepatocyte CELLS (scRNA) vs bulk, disease"),
        ("sn_hep_healthy", "bulk_normal", "Hepatocyte NUCLEI (snRNA) vs bulk, healthy"),
        ("sn_hepauthor_healthy", "bulk_normal", "  same, authors' QC'd nuclei"),
        ("sn_hep_disease", "bulk_esld", "Hepatocyte NUCLEI (snRNA) vs bulk, disease"),
        ("sn_hep_healthy", "sc_hep_healthy", "Hepatocyte nuclei vs cells, same 4 livers")]
vals = [E[(E.a == a) & (E.b == b)].iloc[0] for a, b, _ in rows]
yy = np.arange(len(rows))[::-1]
for yi, rr, (a, b, lab) in zip(yy, vals, rows):
    col = INK if a.startswith("bulk") else (RED if "sn_" in a else BLUE)
    a2.plot([rr.rho_expressed_both, rr.rho_all_pc], [yi, yi], color=GRID, lw=2, zorder=1)
    a2.scatter(rr.rho_all_pc, yi, s=46, facecolor=SURF, edgecolor=col, lw=1.6, zorder=3)
    a2.scatter(rr.rho_expressed_both, yi, s=46, color=col, edgecolor=SURF, lw=1.2, zorder=3)
    a2.text(1.025, yi, f"{rr.rho_expressed_both:.2f}", va="center", fontsize=8, color=INK, clip_on=False)
a2.set_yticks(yy, [l for *_, l in rows], fontsize=8.5)
a2.set_xlim(0.35, 1.0)
a2.set_xlabel("Spearman ρ\nfilled: genes expressed in both (CPM ≥ 1), value printed;  open: all protein-coding genes")
a2.set_title("Profile agreement against the same-platform references", loc="left", fontsize=10, color=INK)
a2.tick_params(axis="y", length=0)
a2.grid(axis="x", color=GRID, lw=.6)
for s in ("top", "right", "left"):
    a2.spines[s].set_visible(False)
fig.tight_layout()
fig.subplots_adjust(wspace=0.75)
fig.savefig(os.path.join(FIG, "figC_expression_spearman.png"), dpi=200)

# ---- fig D ------------------------------------------------------------------------------------------
F = pd.read_csv(os.path.join(R, "conc_gene_flips.csv"))
FO = pd.read_csv(os.path.join(R, "conc_gene_flips_by_origin.csv"))
LAB = {"sn_hep": "Hepatocyte nuclei (snRNA)", "sn_hepauthor": "  same, authors' nuclei", "sc_hep": "Hepatocyte cells (scRNA)",
       "bulk 6v4 draw": "Bulk cut to 6 vs 4 (median of 50)", "bulk 6v4 draws (median)": "Bulk cut to 6 vs 4 (median of 50)"}
fig, axes = plt.subplots(2, 2, figsize=(12.4, 5.6), sharex=True)
for ci, direction in enumerate(("bulk -> pseudobulk", "pseudobulk -> bulk")):
    for ri, (title, T, draw) in enumerate((("All protein-coding genes", F[F.call == "sig (BH<0.05)"], "bulk 6v4 draw"),
                                           ("Hepatocyte-origin genes only (atlas hepatocyte fraction ≥ 0.5)",
                                            FO[FO.origin == "hepatocyte"], "bulk 6v4 draws (median)"))):
        ax = axes[ri, ci]
        T = T[T.direction == direction]
        # conc_gene_flips names the pseudobulk in `alt` (bulk -> pb) or `ref` (pb -> bulk); the by-origin table always in `alt`
        col = "alt" if (direction == "bulk -> pseudobulk" or ri == 1) else "ref"
        order = ["sn_hep", "sn_hepauthor", "sc_hep", draw]
        sel = [T[T[col] == o].iloc[0] for o in order]
        yy = np.arange(len(sel))[::-1]
        for yi, s, o in zip(yy, sel, order):
            left = 0
            for frac, col in ((s.frac_same, BLUE), (s.frac_not_sig, NEUTRAL), (s.frac_opposite, RED)):
                ax.barh(yi, frac, left=left, color=col, height=.62, edgecolor=SURF, lw=2)
                left += frac
            n = int(s.n_tested_alt)
            ax.text(1.01, yi, f"{s.frac_not_sig + s.frac_opposite:.0%} lost  (n = {n:,})", va="center", fontsize=8, color=INK)
        ax.set_yticks(yy, [LAB[o] for o in order], fontsize=8.5)
        ax.set_xlim(0, 1)
        ax.tick_params(axis="y", length=0)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        head = ("Significant in full bulk (ESLD vs donors) → status in …" if direction == "bulk -> pseudobulk"
                else "Significant in … → status in full bulk")
        ax.set_title(f"{head}\n{title}", loc="left", fontsize=9, color=INK)
axes[1, 0].set_xlabel("fraction of genes"); axes[1, 1].set_xlabel("fraction of genes")
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (BLUE, NEUTRAL, RED)]
fig.legend(handles, ["significant, same direction", "not significant", "significant, opposite direction"],
           loc="lower center", ncol=3, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.01))
fig.suptitle("Significant → not significant: hepatocyte pseudobulk (PSC/PBC vs healthy) against isolated-hepatocyte bulk "
             "(ESLD vs donors)\nyardstick: the same bulk cut to the pseudobulk's size, 6 vs 4", x=0.01, ha="left",
             fontsize=10.5, color=INK)
fig.tight_layout(rect=(0, 0.04, 0.93, 0.94))
fig.savefig(os.path.join(FIG, "figD_flips.png"), dpi=200)

# ---- fig E: pathway-rank Spearman ---------------------------------------------------------------------
LIBS = [("hallmark", "Hallmark"), ("kegg", "KEGG legacy"), ("reactome", "Reactome"), ("gobp", "GO BP")]
PR = pd.read_csv(os.path.join(R, "conc_pathway_rank_spearman.csv"))
NES = os.path.join(R, "pathway_nes")
hb_ = pd.read_csv(os.path.join(NES, "bulk_esld_all_hallmark.csv"), index_col=0)
hn_ = pd.read_csv(os.path.join(NES, "pb_sn_hep_hallmark.csv"), index_col=0)
g = hb_.index.intersection(hn_.index)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.2, 5.4), gridspec_kw={"width_ratios": [1, 1.05]})
sig_b, sig_n = hb_.loc[g, "p"] < .05, hn_.loc[g, "p"] < .05
same_dir = np.sign(hb_.loc[g, "NES"]) == np.sign(hn_.loc[g, "NES"])
col = np.where(sig_b & sig_n & same_dir, BLUE, np.where(sig_b & sig_n, RED, np.where(sig_b | sig_n, "#8a8984", NEUTRAL)))
a1.scatter(hb_.loc[g, "NES"], hn_.loc[g, "NES"], s=30, c=col, edgecolor=SURF, lw=1, zorder=3)
lim = (min(hb_.NES.min(), hn_.NES.min()) - .3, max(hb_.NES.max(), hn_.NES.max()) + .3)
a1.plot(lim, lim, color=GRID, lw=1, zorder=1)
a1.axhline(0, color=GRID, lw=.8); a1.axvline(0, color=GRID, lw=.8)
a1.set_xlim(lim); a1.set_ylim(lim)
for t in g[(sig_b != sig_n) | (sig_b & sig_n & ~same_dir)]:
    a1.annotate(t.replace("HALLMARK_", "").replace("_", " ").lower(), (hb_.at[t, "NES"], hn_.at[t, "NES"]),
                xytext=(4, 3), textcoords="offset points", fontsize=6.8, color=INK)
r0 = PR[(PR.lib == "hallmark") & (PR.a == "bulk esld_all") & (PR.b == "sn_hep")].iloc[0]
a1.text(0.02, 0.98, f"50 Hallmark sets, Spearman ρ = {r0.rho_all:.2f}\nblue: significant in both, same direction; red: opposite;\n"
        "dark grey: significant in one only (labelled)", transform=a1.transAxes, va="top", fontsize=8, color=MUTED)
a1.set_xlabel("NES, isolated-hepatocyte bulk (ESLD vs donors)")
a1.set_ylabel("NES, hepatocyte nuclei (PSC/PBC vs healthy)")
a1.set_title("Where each Hallmark pathway ranks", loc="left", fontsize=10, color=INK)
for s in ("top", "right"):
    a1.spines[s].set_visible(False)

ROWS = [("sn_hep", "nuclei", RED), ("sn_hepauthor", "nuclei, authors' cells", "#f0a3a2"), ("sc_hep", "cells (scRNA)", BLUE)]
ytick, ylab, y = [], [], 0
for lib, lname in LIBS[::-1]:
    d = PR[(PR.lib == lib) & (PR.a == "bulk esld_all") & (PR.b == "bulk 6v4 draws")].iloc[0]
    a2.plot([d.rho_all_p5, d.rho_all_p95], [y, y], color=GRID, lw=6, solid_capstyle="butt", zorder=1)
    a2.scatter(d.rho_all, y, s=40, color=INK, zorder=3)
    ytick.append(y); ylab.append(f"{lname}: bulk cut to 6 vs 4")
    for k, lab, c in ROWS:
        y += 1
        r = PR[(PR.lib == lib) & (PR.a == "bulk esld_all") & (PR.b == k)].iloc[0]
        a2.scatter(r.rho_all, y, s=40, color=c, edgecolor=SURF, lw=1, zorder=3)
        a2.text(r.rho_all + .012, y, f"{r.rho_all:.2f}", va="center", fontsize=7.5, color=INK)
        ytick.append(y); ylab.append(f"{lname}: {lab}")
    y += 1.8
a2.set_yticks(ytick, ylab, fontsize=8)
a2.set_xlim(0.35, 1.0)
a2.tick_params(axis="y", length=0)
a2.grid(axis="x", color=GRID, lw=.6)
for s in ("top", "right", "left"):
    a2.spines[s].set_visible(False)
a2.set_xlabel("Spearman ρ of pathway NES with the full isolated-hepatocyte bulk\n(grey band: 5–95% of 50 bulk 6-vs-4 draws, black dot their median)")
a2.set_title("Pathway-rank agreement, all pathways scored in both", loc="left", fontsize=10, color=INK)
fig.tight_layout()
fig.subplots_adjust(wspace=0.55)
fig.savefig(os.path.join(FIG, "figE_pathway_rank_spearman.png"), dpi=200)

# ---- fig F: pathways significant -> not significant ----------------------------------------------------
PF = pd.read_csv(os.path.join(R, "conc_pathway_flips.csv"))
PF = PF[PF.gate == "p<0.05"]
fig, axes = plt.subplots(1, 2, figsize=(12.6, 6.0), sharey=True)
for ax, direction, head in ((axes[0], "bulk -> pseudobulk", "Called in the bulk (p < 0.05) → status in …"),
                            (axes[1], "pseudobulk -> bulk", "Called in … → status in the full bulk")):
    ytick, ylab, y = [], [], 0
    for lib, lname in LIBS[::-1]:
        for k, lab in (("bulk 6v4 draws (median)", "bulk cut to 6 vs 4 (median)"), ("sc_hep", "cells (scRNA)"),
                       ("sn_hepauthor", "nuclei, authors' cells"), ("sn_hep", "nuclei")):
            s = PF[(PF.lib == lib) & (PF.direction == direction) & (PF.other == k)].iloc[0]
            n = max(s.n_tested, 1)
            left = 0
            for v, c in ((s.same / n, BLUE), (s.not_sig / n, NEUTRAL), (s.opposite / n, RED)):
                ax.barh(y, v, left=left, color=c, height=.66, edgecolor=SURF, lw=1.5)
                left += v
            if k.startswith("bulk"):
                ax.plot([1 - s.frac_lost_p95, 1 - s.frac_lost_p5], [y, y], color=INK, lw=1.2)
            ax.text(1.01, y, f"{s.frac_lost:.0%} lost (n = {int(s.n_tested):,})", va="center", fontsize=7.5, color=INK)
            ytick.append(y); ylab.append(f"{lname}: {lab}")
            y += 1
        y += 0.8
    ax.set_yticks(ytick, ylab, fontsize=8)
    ax.set_xlim(0, 1)
    ax.tick_params(axis="y", length=0)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.set_title(head, loc="left", fontsize=9.5, color=INK)
    ax.set_xlabel("fraction of pathways\n(black line on bulk rows: 5–95% of the draws' 'same' share)")
fig.legend([plt.Rectangle((0, 0), 1, 1, color=c) for c in (BLUE, NEUTRAL, RED)],
           ["significant, same direction", "not significant", "significant, opposite direction"],
           loc="lower center", ncol=3, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.005))
fig.suptitle("Pathways that turn significant → not significant: hepatocyte pseudobulk (PSC/PBC vs healthy) vs isolated-hepatocyte "
             "bulk (ESLD vs donors)\nyardstick: the same bulk cut to the pseudobulk's size, 6 vs 4", x=0.01, ha="left",
             fontsize=10.5, color=INK)
fig.tight_layout(rect=(0, 0.05, 0.9, 0.92))
fig.savefig(os.path.join(FIG, "figF_pathway_flips.png"), dpi=200)
print("written:", sorted(os.listdir(FIG)))
