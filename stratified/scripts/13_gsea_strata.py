#!/usr/bin/env python3
"""13_gsea_strata.py — stratified inflammation GSEA: isolated-hepatocyte bulk ( disease strata vs normal
donors) against scRNA and snRNA pseudobulk (PSC/PBC vs healthy livers, whole tissue and hepatocytes only).

There are no matched samples, so each modality is analysed within its own strata and only the GSEA answers
are compared. Every ranking: limma moderated t (10_bulk_strata_limma.R, 12_pb_strata_limma.R), protein-coding
genes, MT-* dropped, one row per symbol (highest AveExpr). GSEA: the Pittsburgh convention imported unchanged
from scripts/ref/common.py (gseapy prerank 1000 perms, size 10-500, seed 1234, MSigDB 2026.1) — gseapy for
Hallmark / KEGG (leading edges), the validated numpy fast_gsea for Reactome / GO BP.
Outputs: results/gsea/<ranking>_<lib>_gsea.csv and results/strata_*.csv
"""
import os, re, sys
import numpy as np, pandas as pd

ROOT = "/srv/public/sayar99/liver_strat_gsea"
sys.path.insert(0, os.path.join(ROOT, "scripts", "ref"))
import common as C  # noqa: E402

C.GMT_DIR = os.path.join(ROOT, "genesets")
D, R = os.path.join(ROOT, "data"), os.path.join(ROOT, "results")
G = os.path.join(R, "gsea")
os.makedirs(G, exist_ok=True)

INFL_H = ["HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE",
          "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
          "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_IL2_STAT5_SIGNALING",
          "HALLMARK_COMPLEMENT", "HALLMARK_ALLOGRAFT_REJECTION"]
INFL_RE = re.compile(r"INFLAM|INTERFERON|CYTOKINE|CHEMOKINE|NF_?KAPPA|NFKB|TNF|INTERLEUKIN|"
                     r"\bIL\d|_IL\d|TOLL|COMPLEMENT|ACUTE_PHASE|ALLOGRAFT", re.I)
# warm-dissociation stress / IEG union (bulktech-assimilation gap_explain::_DISSOC; van den Brink 2017 +
# O'Flanagan 2019, assembled union)
DISSOC = set("""ADM ANKRD1 ARL4D ATF3 ATF4 BAG3 BHLHE40 BRD2 BTG1 BTG2 CCN1 CYR61 CCN2 CTGF CCNL1 CDKN1A CEBPB
CEBPD CSRNP1 CXCL1 CXCL2 DDIT3 DDIT4 DNAJA1 DNAJB1 DNAJB4 DNAJB6 DNAJB9 DUSP1 DUSP5 DUSP8 EGR1 EGR2 EGR3 EGR4
EIF1 ERRFI1 FOS FOSB FOSL1 FOSL2 GADD45B GADD45G GDF15 GEM HERPUD1 HES1 HIF1A HSP90AA1 HSP90AB1 HSPA1A HSPA1B
HSPA5 HSPA8 HSPB1 HSPE1 HSPH1 ID1 ID2 ID3 IER2 IER3 IER5 IER5L INHBA IRF1 JMJD6 JUN JUNB JUND KLF2 KLF4 KLF6 KLF9
KLF10 LMNA MAFF MAP3K8 MCL1 MIDN MT1X MT2A MYC NAMPT NCOA7 NFE2L2 NFKBIA NFKBIE NFKBIZ NR4A1 NR4A2 NR4A3 OSGIN1
PER1 PIM1 PIM3 PLK2 PLK3 PNRC1 PPP1R15A PPP1R15B PTGS2 RASD1 RGS1 RGS2 RHOB SAT1 SERTAD1 SGK1 SIK1 SKIL SLC20A1
SLC25A25 SLC2A3 SLC3A2 SNAI1 SOCS3 SOX9 SQSTM1 THBS1 TIPARP TNFAIP3 TNFAIP6 TRA2A TRA2B TRIB1 TSC22D1 UBC USP2
VEGFA XBP1 ZFP36 ZFP36L1 ZFP36L2 ZYX""".split())
FAST = {"reactome", "gobp"}
BULK = {"bulk_esld_all": "Bulk hepatocytes: all ESLD vs donors",
        "bulk_ALD": "Bulk hepatocytes: ALD vs donors", "bulk_MASH": "Bulk hepatocytes: MASH vs donors",
        "bulk_CP_B": "Bulk hepatocytes: Child-Pugh B vs donors",
        "bulk_CP_C": "Bulk hepatocytes: Child-Pugh C vs donors"}
PB = {"pb_sc_whole": "scRNA pseudobulk, whole liver: PSC/PBC (6) vs healthy (4)",
      "pb_sn_whole": "snRNA pseudobulk, whole liver: PSC/PBC (6) vs healthy (4)",
      "pb_sc_hep": "scRNA pseudobulk, hepatocytes: PSC/PBC vs healthy",
      "pb_sn_hep": "snRNA pseudobulk, hepatocytes: PSC/PBC vs healthy"}


def ranking(name):
    if name.startswith("bulk_"):
        d = pd.read_csv(os.path.join(D, f"{name}_limma.csv"), index_col=0)
        d["symbol"] = d.index
    else:
        d = pd.read_csv(os.path.join(R, f"{name}_limma.csv"), index_col=0)
    d = d[d.protein_coding & d.symbol.notna()]
    d = d[~d.symbol.str.upper().str.startswith("MT-")]
    d = d.sort_values("AveExpr", ascending=False).drop_duplicates("symbol")
    return d.set_index("symbol")["t"].astype(float).sort_values(ascending=False)


def gsea(rnk, pre):
    tabs = C.load_gsea(pre)
    if len(tabs) == len(C.LIB_FILES):
        return tabs
    slow = [l for l in C.LIB_FILES if l not in FAST]
    tabs = C.gsea(rnk, libs=slow)
    for k, t in C.fast_gsea(rnk, libs=sorted(FAST)).items():
        tabs[k] = t.assign(lead=np.nan)[["ES", "NES", "p", "fdr", "lead"]]
    C.save_gsea(tabs, pre)
    print(f"[gsea] {os.path.basename(pre)}", flush=True)
    return tabs


is_infl = lambda lib, x: x in INFL_H if lib == "hallmark" else INFL_RE.search(x) is not None
names = list(BULK) + list(PB)
RNK = {n: ranking(n) for n in names}
TAB = {n: gsea(RNK[n], os.path.join(G, n)) for n in names}

infl, summ, nodis = [], [], []
for n in names:
    for lib, t in TAB[n].items():
        s = t[[is_infl(lib, x) for x in t.index]]
        tl = C.top_lists(t)
        summ.append(dict(ranking=n, lib=lib, n_genes=len(RNK[n]), n_sets=len(t), n_infl=len(s),
                         infl_sig_up=int(((s.NES > 0) & (s.p < .05)).sum()),
                         infl_sig_down=int(((s.NES < 0) & (s.p < .05)).sum()),
                         infl_in_top15_up=int(sum(is_infl(lib, x) for x in tl["top15_up"])),
                         infl_in_bottom15=int(sum(is_infl(lib, x) for x in tl["bottom15_down"]))))
        if lib == "hallmark":
            for x in INFL_H:
                lead = [g for g in str(t.at[x, "lead"]).split(";") if g and g != "nan"]
                infl.append(dict(ranking=n, term=x, NES=t.at[x, "NES"], p=t.at[x, "p"], fdr=t.at[x, "fdr"],
                                 rank_of_50=int(t.NES.rank(ascending=False)[x]), n_lead=len(lead),
                                 lead_dissoc=sum(g in DISSOC for g in lead),
                                 lead_dissoc_genes=";".join(g for g in lead if g in DISSOC)))
    nd = C.gsea(RNK[n][~RNK[n].index.isin(DISSOC)], libs=["hallmark"])["hallmark"]
    for x in INFL_H:
        nodis.append(dict(ranking=n, term=x, NES_all=TAB[n]["hallmark"].at[x, "NES"], p_all=TAB[n]["hallmark"].at[x, "p"],
                          NES_no_dissoc=nd.at[x, "NES"], p_no_dissoc=nd.at[x, "p"]))
INF, SUM, NOD = pd.DataFrame(infl), pd.DataFrame(summ), pd.DataFrame(nodis)

# agreement between every pair of rankings: Hallmark NES r (all 50 / 8 inflammation), gene-level t r
pairs = []
for i, a in enumerate(names):
    for b in names[i + 1:]:
        row = dict(a=a, b=b, gene_t_r=float(RNK[a].corr(RNK[b].reindex(RNK[a].index))),
                   n_shared_genes=int(len(RNK[a].index.intersection(RNK[b].index))))
        for lib in C.LIB_FILES:
            ta, tb = TAB[a][lib], TAB[b][lib]
            com = ta.index.intersection(tb.index)
            it = [x for x in com if is_infl(lib, x)]
            row[f"{lib}_nes_r"] = float(ta.loc[com, "NES"].corr(tb.loc[com, "NES"]))
            row[f"{lib}_infl_nes_r"] = float(ta.loc[it, "NES"].corr(tb.loc[it, "NES"])) if len(it) > 2 else np.nan
        pairs.append(row)
PAIR = pd.DataFrame(pairs)

# list-level replication: each bulk stratum's reported Hallmark/KEGG/Reactome/GO BP 15-lists in each pseudobulk
rep = []
for a in BULK:
    for b in PB:
        for lib in C.LIB_FILES:
            for r in C.compare(TAB[a][lib], TAB[b][lib]):
                r.update(ref=a, alt=b, lib=lib)
                rep.append(r)
REP = pd.DataFrame(rep)

INF.to_csv(os.path.join(R, "strata_hallmark_inflammation.csv"), index=False)
SUM.to_csv(os.path.join(R, "strata_inflammation_summary.csv"), index=False)
NOD.to_csv(os.path.join(R, "strata_inflammation_without_dissoc.csv"), index=False)
PAIR.to_csv(os.path.join(R, "strata_pairwise_agreement.csv"), index=False)
REP.to_csv(os.path.join(R, "strata_bulk_list_replication.csv"), index=False)
pd.Series({**BULK, **PB}, name="label").to_csv(os.path.join(R, "strata_labels.csv"))

pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 50)
print(INF.pivot_table(index="term", columns="ranking", values="NES").reindex(columns=names).round(2).to_string())
print(INF.pivot_table(index="term", columns="ranking", values="p").reindex(columns=names).round(3).to_string())
print(SUM.round(3).to_string(index=False))
print(NOD.round(2).to_string(index=False))
print(PAIR[["a", "b", "gene_t_r", "hallmark_nes_r", "hallmark_infl_nes_r", "kegg_nes_r", "gobp_nes_r", "gobp_infl_nes_r"]].round(2).to_string(index=False))
print(REP[REP.lib == "hallmark"][["ref", "alt", "direction", "n_ref", "overlap", "still_called"]].to_string(index=False))
