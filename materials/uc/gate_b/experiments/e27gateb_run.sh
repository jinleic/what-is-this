#!/bin/bash
# e27 Gate B neighborhood descent around the lowest-defect e26 families.
# Strict caps in the tool grammar admit the e27 seeds and force every walked
# family to defect <= the seed floor; any finalist below it advances the
# frontier.  All runs deterministic (seed 20260831, e27 namespace,
# append-only checkpoints, resume keys include cap).  Max 2 concurrent with
# the n11 seed prep; each pass runs sequentially inside this script.
#
# Usage: bash e27gateb_run.sh
set -u
cd "$(dirname "$0")/../../.."   # math/

# Pass 1: n10 m15, cap 63/225 (limit 62; seeds 62/60/56).
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -I -B \
  uc/gate_b/search_subtwofifths_q.py \
  --dimension 10 --sizes 15 \
  --restarts 12 --steps 7500 --screen-orders 48 \
  --seed 20260831 --score min-a --defect-cap 63/225 \
  --defect-checkpoint uc/gate_b/experiments/e27gateb_seeds_checkpoint.jsonl \
  --extra-seed-checkpoint uc/gate_b/experiments/e26gateb_n10_m15_sub32_checkpoint.jsonl \
  --checkpoint uc/gate_b/experiments/e27gateb_n10_m15_desc_checkpoint.jsonl \
  --report uc/gate_b/candidates/e27gateb_n10_m15_desc_report.json \
  > uc/gate_b/experiments/e27gateb_n10_m15_desc.log 2>&1
echo "e27 n10 m15 done rc=$?" >> uc/gate_b/experiments/e27gateb_status.log

# Pass 2: n10 m20, cap 131/400 (limit 130 = the 13/40 seed).
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -I -B \
  uc/gate_b/search_subtwofifths_q.py \
  --dimension 10 --sizes 20 \
  --restarts 12 --steps 7500 --screen-orders 48 \
  --seed 20260831 --score min-a --defect-cap 131/400 \
  --defect-checkpoint uc/gate_b/experiments/e26gateb_n10_m20_sub34_checkpoint.jsonl \
  --extra-seed-checkpoint uc/gate_b/experiments/e27gateb_seeds_checkpoint.jsonl \
  --checkpoint uc/gate_b/experiments/e27gateb_n10_m20_desc_checkpoint.jsonl \
  --report uc/gate_b/candidates/e27gateb_n10_m20_desc_report.json \
  > uc/gate_b/experiments/e27gateb_n10_m20_desc.log 2>&1
echo "e27 n10 m20 done rc=$?" >> uc/gate_b/experiments/e27gateb_status.log

# Pass 3: n10 m13, cap 65/169 (limit 64 = the seed).
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -I -B \
  uc/gate_b/search_subtwofifths_q.py \
  --dimension 10 --sizes 13 \
  --restarts 12 --steps 7500 --screen-orders 48 \
  --seed 20260831 --score min-a --defect-cap 65/169 \
  --defect-checkpoint uc/gate_b/experiments/e26gateb_n10_m13_sub40_checkpoint.jsonl \
  --checkpoint uc/gate_b/experiments/e27gateb_n10_m13_desc_checkpoint.jsonl \
  --report uc/gate_b/candidates/e27gateb_n10_m13_desc_report.json \
  > uc/gate_b/experiments/e27gateb_n10_m13_desc.log 2>&1
echo "e27 n10 m13 done rc=$?" >> uc/gate_b/experiments/e27gateb_status.log

# Pass 4 (runs last; n11 seeds may still be manufacturing — the pass is
# naturally skipped if the seed pool is empty at that point):
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -I -B \
  uc/gate_b/search_subtwofifths_q.py \
  --dimension 11 --sizes 15 \
  --restarts 12 --steps 7500 --screen-orders 48 \
  --seed 20260831 --score min-a --defect-cap 17/45 \
  --defect-checkpoint uc/gate_b/experiments/e27gateb_n11_m15_fullsetseeds_checkpoint.jsonl \
  --extra-seed-checkpoint uc/gate_b/experiments/e27gateb_n11_m15_hotseed_checkpoint.jsonl \
  --checkpoint uc/gate_b/experiments/e27gateb_n11_m15_desc_checkpoint.jsonl \
  --report uc/gate_b/candidates/e27gateb_n11_m15_desc_report.json \
  > uc/gate_b/experiments/e27gateb_n11_m15_desc.log 2>&1
echo "e27 n11 m15 done rc=$?" >> uc/gate_b/experiments/e27gateb_status.log
echo "e27 complete" >> uc/gate_b/experiments/e27gateb_status.log
