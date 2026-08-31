#!/bin/bash
# e26 Gate B second-wave n=11 m=13 sweep, seeded by the defect-66/169
# full-set family found by e26gateb_seed_prep.py (checkpoint C).  Cap 2/5
# gives strict limit 67 at m=13, so every walked family has defect <= 67/169
# = 0.39645 < 2/5, and any finalist <= 68/169 already beats nothing new but
# ties the sub-14/45 question only through the frontier exact enclosures.
# Deterministic resume keys (score/steps/orders/seed/cap) with e26 namespace.
set -u
cd "$(dirname "$0")/../../.."   # math/

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -I -B \
  uc/gate_b/search_subtwofifths_q.py \
  --dimension 11 --sizes 13 \
  --restarts 12 --steps 7500 --screen-orders 48 \
  --seed 20260830 --score min-a --defect-cap 2/5 \
  --defect-checkpoint uc/gate_b/experiments/e26gateb_n11_m13_fullsetseeds_c_checkpoint.jsonl \
  --extra-seed-checkpoint uc/gate_b/experiments/e26gateb_n11_m13_fullsetseeds_b_checkpoint.jsonl \
  --checkpoint uc/gate_b/experiments/e26gateb_n11_m13_sub40b_checkpoint.jsonl \
  --report uc/gate_b/candidates/e26gateb_n11_m13_sub40b_report.json \
  > uc/gate_b/experiments/e26gateb_n11_m13_sub40b.log 2>&1
echo "n11 m13 sub40b done rc=$?" >> uc/gate_b/experiments/e26gateb_runB_status.log
