#!/bin/bash
# pipeline on compute04: 11 cluster -> 11b shared genes -> 11c label + hepatocyte pseudobulk -> 12 limma (R)
# -> 13 GSEA -> 14 figures. Bulk strata (10_bulk_strata_limma.R) run locally beforehand; only their per-gene
# tables are copied to data/.
cd /srv/public/sayar99/liver_strat_gsea/scripts
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMBA_NUM_THREADS=2
{ ../venv/bin/python -u 11_build_sc_sn_strata.py 4 > ../logs/11_build.log 2>&1 \
  && ../venv/bin/python -u 11b_shared_genes.py > ../logs/11b.log 2>&1 \
  && ../venv/bin/python -u 11c_label_pseudobulk.py > ../logs/11c_label.log 2>&1 \
  && Rscript 12_pb_strata_limma.R > ../logs/12_limma.log 2>&1 \
  && PBV_CPUS=8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ../venv/bin/python -u 13_gsea_strata.py > ../logs/13_gsea.log 2>&1 \
  && ../venv/bin/python -u 14_figures.py > ../logs/14_figures.log 2>&1; } ; echo "PIPE_DONE rc=$?" >> ../logs/pipeline.log
