#!/bin/bash
# e26 Gate B low-defect negative-order sweep — detached run A (n=10).
# Strict defect caps below the tool's 2/5 baseline, tight against available
# full-set-containing seeds; the min-a walk then has room to descend below
# the 14/45 record (70/225).  Any finalist at defect <= 69/225 beats the
# record.  Deterministic: seed 20260830 e26-namespace, restart-indexed
# append-only checkpoints, resume keys include score/steps/orders/seed/cap.
#
# Usage: bash e26gateb_runA_n10.sh
set -u
cd "$(dirname "$0")/../../.."   # math/

SIZES_A="13,15,20"
CAP_A="73/225"
# cap 73/225 has limit 72 at m=15 (admits the 8/25 seed), 60 at m=13, 126 at
# m=20 (admits the 17/50 seed 136? 137/400 needed for that).
# m=20 needs cap 137/400; run it as a second pass.
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -I -B \
  uc/gate_b/search_subtwofifths_q.py \
  --dimension 10 --sizes 15 \
  --restarts 12 --steps 7500 --screen-orders 48 \
  --seed 20260830 --score min-a --defect-cap 73/225 \
  --defect-checkpoint uc/gate_b/experiments/local_defect_n10_checkpoint.jsonl \
  --checkpoint uc/gate_b/experiments/e26gateb_n10_m15_sub32_checkpoint.jsonl \
  --report uc/gate_b/candidates/e26gateb_n10_m15_sub32_report.json \
  > uc/gate_b/experiments/e26gateb_n10_m15_sub32.log 2>&1
echo "pass m=15 done rc=$?" >> uc/gate_b/experiments/e26gateb_runA_status.log

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -I -B \
  uc/gate_b/search_subtwofifths_q.py \
  --dimension 10 --sizes 20 \
  --restarts 12 --steps 7500 --screen-orders 48 \
  --seed 20260830 --score min-a --defect-cap 137/400 \
  --defect-checkpoint uc/gate_b/experiments/local_defect_n10_checkpoint.jsonl \
  --checkpoint uc/gate_b/experiments/e26gateb_n10_m20_sub34_checkpoint.jsonl \
  --report uc/gate_b/candidates/e26gateb_n10_m20_sub34_report.json \
  > uc/gate_b/experiments/e26gateb_n10_m20_sub34.log 2>&1
echo "pass m=20 done rc=$?" >> uc/gate_b/experiments/e26gateb_runA_status.log

# m=13 pass waits for the pinned m=13 seed prep; cap 2/5 uses limit 67 which
# the prep must reach.  If the prep's checkpoint has no qualifying family the
# pass is skipped (the tool prints "no dominant seed below" and continues).
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -I -B \
  uc/gate_b/search_subtwofifths_q.py \
  --dimension 10 --sizes 13 \
  --restarts 12 --steps 7500 --screen-orders 48 \
  --seed 20260830 --score min-a --defect-cap 2/5 \
  --defect-checkpoint uc/gate_b/experiments/e26gateb_n10_m13_fullsetseeds_checkpoint.jsonl \
  --extra-seed-checkpoint uc/gate_b/experiments/local_defect_n10_checkpoint.jsonl \
  --checkpoint uc/gate_b/experiments/e26gateb_n10_m13_sub40_checkpoint.jsonl \
  --report uc/gate_b/candidates/e26gateb_n10_m13_sub40_report.json \
  > uc/gate_b/experiments/e26gateb_n10_m13_sub40.log 2>&1
echo "pass m=13 done rc=$?" >> uc/gate_b/experiments/e26gateb_runA_status.log
echo "runA complete" >> uc/gate_b/experiments/e26gateb_runA_status.log
