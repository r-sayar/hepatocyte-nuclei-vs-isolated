"""liver_common.py — shared helpers for the liver bulk / sc / sn pseudobulk GSEA comparison.

Pathway convention is imported unchanged from bulktech-assimilation/scripts/pathway_bulk_vs_pb/common.py
(copied to scripts/ref/): DESeq2 rowSums>=10, Wald-stat ranking, MT-* dropped, symbol de-dup by
baseMean, gseapy prerank 1000 perms, size 10-500, seed 1234, MSigDB 2026.1.
"""
import gzip, io, os, re, sys, zipfile
import numpy as np, pandas as pd

ROOT = "/srv/public/sayar99/liver_pb_gsea"
sys.path.insert(0, os.path.join(ROOT, "scripts", "ref"))
import common as C  # noqa: E402

C.GMT_DIR = os.path.join(ROOT, "genesets")
RES = os.path.join(ROOT, "results")
D = os.path.join(ROOT, "data")
A = os.path.join(D, "andrews")
E = os.path.join(D, "encode")

INFL_H = ["HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE",
          "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
          "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_IL2_STAT5_SIGNALING",
          "HALLMARK_COMPLEMENT", "HALLMARK_ALLOGRAFT_REJECTION"]
INFL_RE = re.compile(r"INFLAM|INTERFERON|CYTOKINE|CHEMOKINE|NF_?KAPPA|NFKB|TNF|INTERLEUKIN|"
                     r"\bIL\d|_IL\d|TOLL|COMPLEMENT|ACUTE_PHASE|ALLOGRAFT", re.I)

# Dissociation-stress / IEG union used across the thesis project
# (bulktech-assimilation/scripts/gap_explain/build_gene_features_rich.py::_DISSOC;
#  van den Brink 2017 Nat Methods + O'Flanagan 2019 Genome Biol, assembled union)
DISSOC = {
    "ADM", "ANKRD1", "ARL4D", "ATF3", "ATF4", "BAG3", "BHLHE40", "BRD2", "BTG1",
    "BTG2", "CCN1", "CYR61", "CCN2", "CTGF", "CCNL1", "CDKN1A", "CEBPB",
    "CEBPD", "CSRNP1", "CXCL1", "CXCL2", "DDIT3", "DDIT4", "DNAJA1", "DNAJB1",
    "DNAJB4", "DNAJB6", "DNAJB9", "DUSP1", "DUSP5", "DUSP8", "EGR1", "EGR2",
    "EGR3", "EGR4", "EIF1", "ERRFI1", "FOS", "FOSB", "FOSL1", "FOSL2",
    "GADD45B", "GADD45G", "GDF15", "GEM", "HERPUD1", "HES1", "HIF1A",
    "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA5", "HSPA8", "HSPB1",
    "HSPE1", "HSPH1", "ID1", "ID2", "ID3", "IER2", "IER3", "IER5", "IER5L",
    "INHBA", "IRF1", "JMJD6", "JUN", "JUNB", "JUND", "KLF2", "KLF4", "KLF6",
    "KLF9", "KLF10", "LMNA", "MAFF", "MAP3K8", "MCL1", "MIDN", "MT1X", "MT2A",
    "MYC", "NAMPT", "NCOA7", "NFE2L2", "NFKBIA", "NFKBIE", "NFKBIZ", "NR4A1",
    "NR4A2", "NR4A3", "OSGIN1", "PER1", "PIM1", "PIM3", "PLK2", "PLK3",
    "PNRC1", "PPP1R15A", "PPP1R15B", "PTGS2", "RASD1", "RGS1", "RGS2", "RHOB",
    "SAT1", "SERTAD1", "SGK1", "SIK1", "SKIL", "SLC20A1", "SLC25A25",
    "SLC2A3", "SLC3A2", "SNAI1", "SOCS3", "SOX9", "SQSTM1", "THBS1",
    "TIPARP", "TNFAIP3", "TNFAIP6", "TRA2A", "TRA2B", "TRIB1", "TSC22D1",
    "UBC", "USP2", "VEGFA", "XBP1", "ZFP36", "ZFP36L1", "ZFP36L2", "ZYX",
}

_ann = pd.read_csv(os.path.join(D, "ensg2symbol_v46.tsv"), sep="\t", header=None,
                   names=["ensg", "symbol", "type"]).drop_duplicates("ensg").set_index("ensg")
SYM, GTYPE = _ann["symbol"], _ann["type"]
LEN = pd.read_csv(os.path.join(D, "gene_lengths_v46.tsv"), sep="\t", index_col=0)


def is_infl(lib, term):
    return term in INFL_H if lib == "hallmark" else INFL_RE.search(term) is not None


def symbolise(res):
    r = res.copy()
    r["ensg"] = r.index
    r.index = [SYM.get(e, e) for e in r.index]
    return r


def protein_coding(counts):
    return counts.loc[[g for g in counts.index if GTYPE.get(g) == "protein_coding"]]


FAST_LIBS = {"reactome", "gobp"}


def gsea(rnk, libs=None):
    """gseapy prerank for Hallmark / KEGG (keeps leading edges); the project's numpy fast_gsea for
    Reactome / GO BP, which reproduces gseapy NES at r = 0.9999 (pathway_bulk_vs_pb/03_validate_fast_gsea.py)."""
    want = libs or list(C.LIB_FILES)
    slow = [l for l in want if l not in FAST_LIBS]
    out = C.gsea(rnk, libs=slow) if slow else {}
    fast = [l for l in want if l in FAST_LIBS]
    if fast:
        for k, t in C.fast_gsea(rnk, libs=fast).items():
            out[k] = t.assign(lead=np.nan)[["ES", "NES", "p", "fdr", "lead"]]
    return {l: out[l] for l in want}


def run(counts, md, design, contrast, pre, libs=None, rank_fn=None):
    """DESeq2 (cached) + GSEA (cached). rank_fn(res) -> ranked Series; default Pittsburgh C.ranks."""
    os.makedirs(os.path.dirname(pre), exist_ok=True)
    want = libs or list(C.LIB_FILES)
    if os.path.exists(f"{pre}_deseq2.csv"):
        res = pd.read_csv(f"{pre}_deseq2.csv", index_col=0)
    else:
        res = symbolise(C.deseq(counts, md, design, contrast))
        res.to_csv(f"{pre}_deseq2.csv")
    tabs = C.load_gsea(pre, libs=want)
    if len(tabs) < len(want):
        tabs = gsea((rank_fn or C.ranks)(res), libs=want)
        C.save_gsea(tabs, pre)
        print(f"[run] {os.path.basename(pre)}: {counts.shape[1]} libs, {res.stat.notna().sum()} genes", flush=True)
    return res, tabs


def length_adjusted_ranks(res, nbins=10):
    """Wald stat minus the median stat of genes with similar exonic length and genomic span
    (nbins x nbins quantile grid on log10 lengths). Same MT / de-dup rules as C.ranks."""
    d = res.dropna(subset=["stat"]).copy()
    L = LEN.reindex(d["ensg"].values)
    d["lx"] = np.log10(L.exonic_len.values)
    d["ls"] = np.log10(L.span.values)
    d = d.dropna(subset=["lx", "ls"])
    bx = pd.qcut(d.lx, nbins, labels=False, duplicates="drop")
    bs = pd.qcut(d.ls, nbins, labels=False, duplicates="drop")
    d["stat"] = d.stat - d.groupby([bx, bs]).stat.transform("median")
    return C.ranks(d)


def length_tilt(res):
    d = res.dropna(subset=["stat"])
    L = LEN.reindex(d["ensg"].values)
    return dict(rho_exonic_len=float(pd.Series(d.stat.values).corr(pd.Series(np.log10(L.exonic_len.values)), method="spearman")),
                rho_span=float(pd.Series(d.stat.values).corr(pd.Series(np.log10(L.span.values)), method="spearman")),
                rho_intron_frac=float(pd.Series(d.stat.values).corr(pd.Series(L.intron_frac.values), method="spearman")))


def infl_table(tabs, label):
    rows = []
    for lib, t in tabs.items():
        sel = t[[is_infl(lib, x) for x in t.index]]
        for term, r in sel.iterrows():
            lead = [g for g in str(r.get("lead", "")).split(";") if g and g != "nan"]
            rows.append(dict(analysis=label, lib=lib, term=term, NES=r.NES, p=r.p, fdr=r.fdr,
                             n_lead=len(lead), lead_dissoc=sum(g in DISSOC for g in lead),
                             lead_dissoc_genes=";".join(g for g in lead if g in DISSOC)))
    return pd.DataFrame(rows)


def infl_summary(tabs, label):
    """Inflammation-related sets called (nominal p<.05) each way, and how many reach the top-15 lists."""
    out = []
    for lib, t in tabs.items():
        s = t[[is_infl(lib, x) for x in t.index]]
        tl = C.top_lists(t)
        out.append(dict(analysis=label, lib=lib, n_sets=len(t), n_infl_sets=len(s),
                        frac_all_sig=float((t.p < .05).mean()),
                        infl_sig_up=int(((s.NES > 0) & (s.p < .05)).sum()),
                        infl_sig_down=int(((s.NES < 0) & (s.p < .05)).sum()),
                        infl_in_top15_up=int(sum(is_infl(lib, x) for x in tl["top15_up"])),
                        infl_in_bottom15=int(sum(is_infl(lib, x) for x in tl["bottom15_down"]))))
    return pd.DataFrame(out)


def without_dissoc(res, label, rank_fn=None):
    """Hallmark inflammation sets re-scored with the dissociation-stress genes removed from the ranking."""
    rnk = (rank_fn or C.ranks)(res)
    a = C.gsea(rnk, libs=["hallmark"])["hallmark"]
    b = C.gsea(rnk[~rnk.index.isin(DISSOC)], libs=["hallmark"])["hallmark"]
    keep = [x for x in INFL_H if x in a.index]
    return pd.DataFrame({"analysis": label, "term": keep, "NES_all": a.loc[keep, "NES"].values,
                         "p_all": a.loc[keep, "p"].values,
                         "NES_no_dissoc": b.reindex(keep)["NES"].values, "p_no_dissoc": b.reindex(keep)["p"].values})


def andrews_shared_ensg():
    """ENSG present in EVERY Andrews library's annotation (pseudobulk fillna(0) would fake zeros otherwise)."""
    sets = []
    for f in sorted(os.listdir(A)):
        if f.endswith("_features.tsv.gz"):
            sets.append(set(pd.read_csv(os.path.join(A, f), sep="\t", header=None)[0].str.split(".").str[0]))
        elif f.endswith("_raw_counts.zip"):
            with zipfile.ZipFile(os.path.join(A, f)) as z:
                n = next(n for n in z.namelist() if os.path.basename(n).split(".")[0] in ("features", "genes") and n.endswith(".gz"))
                sets.append(set(pd.read_csv(io.BytesIO(gzip.decompress(z.read(n))), sep="\t", header=None)[0].str.split(".").str[0]))
    return set.intersection(*sets)


def load_encode():
    Bk = pd.read_csv(os.path.join(E, "bulk_counts_ensg.csv"), index_col=0)
    Sn = pd.read_csv(os.path.join(E, "pseudobulk_counts_ensg.csv"), index_col=0)
    g = Bk.index.intersection(Sn.index)
    return Bk.loc[g], Sn.loc[g]


def load_andrews(which="rule"):
    PB = pd.read_csv(os.path.join(A, f"pseudobulk_{which}_ensg.csv"), index_col=0)
    LIB = pd.read_csv(os.path.join(A, "libraries.tsv"), sep="\t", index_col=0)
    shared = sorted(andrews_shared_ensg())
    return PB.reindex(shared).fillna(0).astype(int), LIB


def paired_md(cols, unit, factor):
    return pd.DataFrame({unit: [c.rsplit("_", 1)[0] if "__" not in c else c.split("__")[0] for c in cols],
                         factor: [c.rsplit("_", 1)[1] if "__" not in c else c.split("__")[1] for c in cols]}, index=cols)
