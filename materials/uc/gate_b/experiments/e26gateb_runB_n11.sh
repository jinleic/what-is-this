#!/bin/bash
# e26 Gate B low-defect negative-order sweep — detached run B (n=11).
# Seeds manufactured by e26gateb_seed_prep.py (full-set-pinned annealing);
# caps at 45/113 for m=13 (limit 67) and 73/225 for m=15 (limit 72).  Each
# pass is skipped naturally if its seed pool is empty.  Resume keys include
# score/steps/orders/seed/cap; checkpoints append-only with e26 namespace.
#
# Usage: bash e26gateb_runB_n11.sh
set -u
cd "$(dirname "$0")/../../.."   # math/

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -I -B \
  uc/gate_b/search_subtwofifths_q.py \
  --dimension 11 --sizes 13 \
  --restarts 12 --steps 7500 --screen-orders 48 \
  --seed 20260830 --score min-a --defect-cap 45/113 \
  --defect-checkpoint uc/gate_b/experiments/e26gateb_n11_fullsetseeds_checkpoint.jsonl \
  --extra-seed-checkpoint uc/gate_b/experiments/local_defect_n10_checkpoint.jsonl \
  --checkpoint uc/gate_b/experiments/e26gateb_n11_m13_sub40_checkpoint.jsonl \
  --report uc/gate_b/candidates/e26gateb_n11_m13_sub40_report.json \
  > uc/gate_b/experiments/e26gateb_n11_m13_sub40.log 2>&1
echo "n11 pass m=13 done rc=$?" >> uc/gate_b/experiments/e26gateb_runB_status.log

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -I -B \
  uc/gate_b/search_subtwofifths_q.py \
  --dimension 11 --sizes 15 \
  --restarts 12 --steps 7500 --screen-orders 48 \
  --seed 20260830 --score min-a --defect-cap 73/225 \
  --defect-checkpoint uc/gate_b/experiments/e26gateb_n11_fullsetseeds_checkpoint.jsonl \
  --extra-seed-checkpoint uc/gate_b/experiments/local_defect_n10_checkpoint.jsonl \
  --checkpoint uc/gate_b/experiments/e26gateb_n11_m15_sub32_checkpoint.jsonl \
  --report uc/gate_b/candidates/e26gateb_n11_m15_sub32_report.json \
  > uc/gate_b/experiments/e26gateb_n11_m15_sub32.log 2>&1
echo "n11 pass m=15 done rc=$?" >> uc/gate_b/experiments/e26gateb_runB_status.log
echo "runB complete" >> uc/gate_b/experiments/e26gateb_runB_status.log
