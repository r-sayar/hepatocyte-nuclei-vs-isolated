#!/usr/bin/env python3
"""14_figures.py — figures for the stratified bulk / scRNA / snRNA inflammation-GSEA comparison.

figA_strata_inflammation_heatmap.png  Hallmark inflammation NES per ranking (disease vs healthy in every column)
figB_strata_agreement_matrix.png      Hallmark NES correlation between every pair of rankings
Style follows the paired analysis' fig1 (reference diverging pair, values printed, bold = nominal p < 0.05).
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

R = "/srv/public/sayar99/liver_strat_gsea/results"
FIG = os.path.join(R, "figures")
os.makedirs(FIG, exist_ok=True)
INK, MUTED, GRID, SURF = "#1f1f1e", "#6b6b67", "#d9d8d4", "#fcfcfb"
CMAP = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#f0efec", "#e34948"])
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": GRID,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": INK,
                     "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF})
INFL = ["HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_IL2_STAT5_SIGNALING",
        "HALLMARK_COMPLEMENT", "HALLMARK_ALLOGRAFT_REJECTION"]
SHORT = {"HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNFα via NF-κB", "HALLMARK_INFLAMMATORY_RESPONSE": "Inflammatory response",
         "HALLMARK_INTERFERON_ALPHA_RESPONSE": "IFN-α response", "HALLMARK_INTERFERON_GAMMA_RESPONSE": "IFN-γ response",
         "HALLMARK_IL6_JAK_STAT3_SIGNALING": "IL6–JAK–STAT3", "HALLMARK_IL2_STAT5_SIGNALING": "IL2–STAT5",
         "HALLMARK_COMPLEMENT": "Complement", "HALLMARK_ALLOGRAFT_REJECTION": "Allograft rejection"}
COLS = [("bulk_esld_all", "Bulk\nall ESLD"), ("bulk_ALD", "Bulk\nALD"), ("bulk_MASH", "Bulk\nMASH"),
        ("bulk_CP_B", "Bulk\nChild-Pugh B"), ("bulk_CP_C", "Bulk\nChild-Pugh C"),
        ("pb_sc_whole", "scRNA\nwhole liver"), ("pb_sn_whole", "snRNA\nwhole liver"),
        ("pb_sc_hep", "scRNA\nhepatocytes"), ("pb_sn_hep", "snRNA\nhepatocytes")]

inf = pd.read_csv(os.path.join(R, "strata_hallmark_inflammation.csv"))
N = inf.pivot_table(index="term", columns="ranking", values="NES").reindex(index=INFL, columns=[c for c, _ in COLS])
P = inf.pivot_table(index="term", columns="ranking", values="p").reindex(index=INFL, columns=[c for c, _ in COLS])
N.to_csv(os.path.join(FIG, "figA_values_NES.csv")); P.to_csv(os.path.join(FIG, "figA_values_p.csv"))

fig, ax = plt.subplots(figsize=(11.6, 4.8))
lim = float(np.nanmax(np.abs(N.values)))
ax.imshow(N.values, cmap=CMAP, norm=TwoSlopeNorm(0, -lim, lim), aspect="auto")
for i in range(N.shape[0]):
    for j in range(N.shape[1]):
        v, p = N.iat[i, j], P.iat[i, j]
        ax.text(j, i, f"{v:+.2f}" + ("" if p < .05 else "\nn.s."), ha="center", va="center", fontsize=8.5,
                color=INK, fontweight="bold" if p < .05 else "normal")
for j in (4.5, 6.5):
    ax.axvline(j, color=SURF, lw=4)
ax.set_xticks(range(len(COLS)), [l for _, l in COLS], fontsize=8.5)
ax.set_yticks(range(len(INFL)), [SHORT[t] for t in INFL])
ax.tick_params(length=0)
ax.xaxis.tick_top()
for s in ax.spines.values():
    s.set_visible(False)
ax.text(2, 8.0, "Bulk: isolated hepatocytes, disease vs normal donors", ha="center", va="top", fontsize=8, color=MUTED)
ax.text(6.5, 8.0, "PSC/PBC (6) vs healthy (4) livers, same livers for sc and sn", ha="center", va="top", fontsize=8, color=MUTED)
fig.suptitle("No matched samples, so each modality is compared by disease stratum: Hallmark inflammation sets (GSEA NES)",
             x=0.01, ha="left", fontsize=11, color=INK, y=0.995)
fig.text(0.01, 0.012, "Red = higher in disease. Bold = nominal p < 0.05. limma-voom moderated t on protein-coding genes "
         "(bulk adjusted for MT%, age, sex, depth), gseapy prerank 1000 permutations.", fontsize=8, color=MUTED)
fig.tight_layout(rect=(0, 0.04, 1, 0.97))
fig.savefig(os.path.join(FIG, "figA_strata_inflammation_heatmap.png"), dpi=200)

# ---- agreement matrix ------------------------------------------------------------------------------
pa = pd.read_csv(os.path.join(R, "strata_pairwise_agreement.csv"))
keys = [c for c, _ in COLS]
M = pd.DataFrame(np.eye(len(keys)), index=keys, columns=keys)
for _, r in pa.iterrows():
    M.at[r.a, r.b] = M.at[r.b, r.a] = r.hallmark_nes_r
M.to_csv(os.path.join(FIG, "figB_values_hallmark_nes_r.csv"))
fig, ax = plt.subplots(figsize=(7.4, 6.4))
ax.imshow(M.values, cmap=CMAP, norm=TwoSlopeNorm(0, -1, 1))
for i in range(len(keys)):
    for j in range(len(keys)):
        if i != j:
            ax.text(j, i, f"{M.iat[i, j]:.2f}", ha="center", va="center", fontsize=8, color=INK)
lab = [l.replace("\n", " ") for _, l in COLS]
ax.set_xticks(range(len(keys)), lab, rotation=40, ha="right", fontsize=8)
ax.set_yticks(range(len(keys)), lab, fontsize=8)
for v in (4.5, 6.5):
    ax.axvline(v, color=SURF, lw=3); ax.axhline(v, color=SURF, lw=3)
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
ax.set_title("How alike are the GSEA answers? Hallmark NES correlation, all 50 sets", loc="left", fontsize=10.5, color=INK)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "figB_strata_agreement_matrix.png"), dpi=200)
print("written:", sorted(os.listdir(FIG)))
