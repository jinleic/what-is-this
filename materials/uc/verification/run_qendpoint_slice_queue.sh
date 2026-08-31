#!/bin/bash
# Detached scheduler entry point for qendpoint-lift slices 02..16.
#
# Launch:   setsid nohup uc/verification/run_qendpoint_slice_queue.sh \
#               >> uc/verification/results/queue_logs/queue.log 2>&1 &
# Status:   python3 -c 'import json;print(json.load(open(
#               "uc/verification/results/queue_state/queue.status")))'
# Tail:     tail -f uc/verification/results/queue_logs/sliceKK.log
#
# Idempotent: finished slices (report + lineage close on disk, queue.status
# stage=done) are skipped; partial slices rerun only their missing stage.
set -u
cd "$(dirname "$0")/../.." || exit 1
exec /usr/bin/nice -n 19 /usr/bin/env -i \
  HOME="$HOME" PATH="/usr/bin:/bin:/usr/local/bin" \
  PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  ./.venv/bin/python -I -B -u \
    uc/verification/run_qendpoint_slice_queue.py "$@"
