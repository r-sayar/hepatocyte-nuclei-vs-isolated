#!/usr/bin/env python3
"""06_figures.py — figures for the liver bulk / sc / sn inflammation-GSEA note (protein-coding universe).

fig1_platform_heatmap.png   Hallmark inflammation NES per platform comparison (rows = sets)
fig2_contrast_sc_vs_sn.png  PSC/PBC vs healthy: NES from sc-pseudobulk vs sn-pseudobulk, 50 Hallmark sets,
                            inflammation sets with leave-one-liver-out ranges
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

RES = "/srv/public/sayar99/liver_pb_gsea/results"
FIG = os.path.join(RES, "figures")
os.makedirs(FIG, exist_ok=True)
INK, MUTED, GRID, SURF = "#1f1f1e", "#6b6b67", "#d9d8d4", "#fcfcfb"
BLUE, MID, RED = "#2a78d6", "#f0efec", "#e34948"  # reference diverging pair + neutral midpoint
CMAP = LinearSegmentedColormap.from_list("div", [BLUE, MID, RED])
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": GRID,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": INK,
                     "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF})

INFL = ["HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_IL2_STAT5_SIGNALING",
        "HALLMARK_COMPLEMENT", "HALLMARK_ALLOGRAFT_REJECTION"]
short = lambda t: t.replace("HALLMARK_", "").replace("_", " ").replace("SIGNALING", "sig.").title() \
    .replace("Tnfa", "TNFα").replace("Nfkb", "NF-κB").replace("Il6 Jak Stat3", "IL6–JAK–STAT3") \
    .replace("Il2 Stat5", "IL2–STAT5").replace("Interferon Alpha", "IFN-α").replace("Interferon Gamma", "IFN-γ") \
    .replace("Via", "via")

COLS = [  # (label, file, reference platform)
    ("snRNA vs bulk\nENCODE, 5 pairs", "PC/P1_encode_sn_vs_bulk_hallmark_gsea.csv"),
    ("snRNA vs bulk\nlength-adjusted", "PC/P1_encode_sn_vs_bulk_lenadj_hallmark_gsea.csv"),
    ("scRNA vs snRNA\nall 10 livers", "PC/P2_andrews_sc_vs_sn_all_hallmark_gsea.csv"),
    ("scRNA vs snRNA\nhealthy 4", "PC/P2_andrews_sc_vs_sn_healthy_hallmark_gsea.csv"),
    ("scRNA vs snRNA\nPSC/PBC 6", "PC/P2_andrews_sc_vs_sn_disease_hallmark_gsea.csv"),
    ("scRNA vs snRNA\nlength-adjusted", "PC/P2_andrews_sc_vs_sn_all_lenadj_hallmark_gsea.csv"),
    ("scRNA vs snRNA\nhepatocytes only", "CT/sc_vs_sn_Hepatocyte_hallmark_gsea.csv"),
    ("scRNA vs snRNA\nsame cell mix", "CT/sc_vs_sn_composition_matched_hallmark_gsea.csv"),
]
N, P = pd.DataFrame(index=INFL), pd.DataFrame(index=INFL)
RANK = pd.DataFrame(index=INFL)
for lab, f in COLS:
    t = pd.read_csv(os.path.join(RES, f), index_col=0)
    N[lab], P[lab] = t.NES.reindex(INFL), t.p.reindex(INFL)
    RANK[lab] = t.NES.rank(ascending=False).reindex(INFL)
N.to_csv(os.path.join(FIG, "fig1_values_NES.csv")); P.to_csv(os.path.join(FIG, "fig1_values_p.csv"))

fig, ax = plt.subplots(figsize=(11.2, 4.6))
lim = float(np.nanmax(np.abs(N.values)))
ax.imshow(N.values, cmap=CMAP, norm=TwoSlopeNorm(0, -lim, lim), aspect="auto")
for i in range(N.shape[0]):
    for j in range(N.shape[1]):
        v, p = N.iat[i, j], P.iat[i, j]
        ax.text(j, i, f"{v:+.2f}" + ("" if p < .05 else "\nn.s."), ha="center", va="center", fontsize=8.5,
                color=INK, fontweight="bold" if p < .05 else "normal")
for j in (1.5, 5.5):
    ax.axvline(j, color=SURF, lw=4)
ax.set_xticks(range(N.shape[1]), N.columns, fontsize=8.5)
ax.set_yticks(range(N.shape[0]), [short(t) for t in INFL])
ax.tick_params(length=0)
ax.xaxis.tick_top()
for s in ax.spines.values():
    s.set_visible(False)
fig.suptitle("Same liver, different library: Hallmark inflammation sets (protein-coding genes, GSEA NES)",
             x=0.01, ha="left", fontsize=11, color=INK, y=0.995)
fig.text(0.01, 0.012, "Red = higher in the first-named library (snRNA pseudobulk vs bulk; scRNA vs snRNA pseudobulk). "
         "Bold = nominal p < 0.05.\nPaired DESeq2 (~sample + platform), Wald-stat ranking, gseapy prerank 1000 permutations. "
         "Length-adjusted = Wald stat minus the median of genes with similar exonic length and span.",
         fontsize=8, color=MUTED)
fig.tight_layout(rect=(0, 0.07, 1, 0.97))
fig.savefig(os.path.join(FIG, "fig1_platform_heatmap.png"), dpi=200)

# ---- fig 2: PSC/PBC vs healthy, sc vs sn -------------------------------------------------------
sc = pd.read_csv(os.path.join(RES, "PC/P3_disease_vs_healthy_sc_hallmark_gsea.csv"), index_col=0)
sn = pd.read_csv(os.path.join(RES, "PC/P3_disease_vs_healthy_sn_hallmark_gsea.csv"), index_col=0)
jk = pd.read_csv(os.path.join(RES, "pc_P4_jackknife_hallmark.csv"))
terms = sc.index.intersection(sn.index)
LABEL_OFFSET = {  # (dx, dy, ha) in points, to keep the crowded upper-right labels apart
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": (10, 6, "left"), "HALLMARK_IL6_JAK_STAT3_SIGNALING": (-8, 6, "right"),
    "HALLMARK_IL2_STAT5_SIGNALING": (10, -10, "left"), "HALLMARK_COMPLEMENT": (-8, -12, "right"),
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": (-8, 6, "right"), "HALLMARK_INTERFERON_GAMMA_RESPONSE": (-8, 6, "right")}
fig, ax = plt.subplots(figsize=(6.6, 6.2))
other = [t for t in terms if t not in INFL]
ax.scatter(sn.loc[other, "NES"], sc.loc[other, "NES"], s=26, color="#b9b8b3", edgecolor=SURF, lw=1, zorder=2,
           label="other Hallmark sets")
lo = min(sn.NES.min(), sc.NES.min()) - .3
hi = max(sn.NES.max(), sc.NES.max()) + .3
ax.plot([lo, hi], [lo, hi], color=GRID, lw=1, zorder=1)
ax.axhline(0, color=GRID, lw=.8, zorder=1); ax.axvline(0, color=GRID, lw=.8, zorder=1)
for t in INFL:
    if t not in terms:
        continue
    x, y = sn.at[t, "NES"], sc.at[t, "NES"]
    jx = jk[(jk.term == t) & (jk.platform == "sn")].NES
    jy = jk[(jk.term == t) & (jk.platform == "sc")].NES
    ax.errorbar(x, y, xerr=[[x - jx.min()], [jx.max() - x]], yerr=[[y - jy.min()], [jy.max() - y]],
                fmt="o", ms=7, color=RED, mec=SURF, mew=1.2, ecolor="#f0a3a2", elinewidth=1.2, capsize=0, zorder=3)
    dx, dy, ha = LABEL_OFFSET.get(t, (6, 4, "left"))
    ax.annotate(short(t), (x, y), xytext=(dx, dy), textcoords="offset points", fontsize=8, color=INK, ha=ha)
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
ax.set_xlabel("NES from snRNA pseudobulk  (PSC/PBC vs healthy)")
ax.set_ylabel("NES from scRNA pseudobulk  (same 10 livers)")
ax.set_title("PSC/PBC vs healthy liver: does sc or sn change the inflammation call?", loc="left", fontsize=10.5, color=INK)
r = sc.loc[terms, "NES"].corr(sn.loc[terms, "NES"])
ax.text(0.98, 0.02, f"50 Hallmark sets, NES r = {r:.2f}\nred = inflammation sets\nbars = leave-one-liver-out range",
        transform=ax.transAxes, va="bottom", ha="right", fontsize=8, color=MUTED)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig2_contrast_sc_vs_sn.png"), dpi=200)
pd.DataFrame({"NES_sc": sc.NES, "p_sc": sc.p, "NES_sn": sn.NES, "p_sn": sn.p}).to_csv(os.path.join(FIG, "fig2_values.csv"))
print("figures written:", os.listdir(FIG))
print(N.round(2).to_string()); print(RANK.astype(int).to_string())
