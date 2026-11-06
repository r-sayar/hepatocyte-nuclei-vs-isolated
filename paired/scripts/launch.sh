#!/bin/bash
# kill leftover joblib/loky workers from this project's venv, then launch 03/04/05
pkill -u sayar99 -f "liver_pb_gsea/(scripts/\.\./)?venv/bin/python" ; sleep 2
cd /srv/public/sayar99/liver_pb_gsea/scripts
for s in 03_gsea_liver 04_protein_coding_length 05_celltype; do
  log=../logs/$s.log
  nohup bash -c "PBV_CPUS=8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ../venv/bin/python -u $s.py > $log 2>&1; echo CHAIN_DONE >> $log" >/dev/null 2>&1 &
done
sleep 3
pgrep -u sayar99 -f "venv/bin/python -u" -a | cut -c1-120
