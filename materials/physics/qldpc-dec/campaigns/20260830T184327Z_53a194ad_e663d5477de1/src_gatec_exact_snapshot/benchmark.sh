#!/bin/bash
# Gate C benchmark: exact degenerate-ML reference vs the harness decoders
# on BB[[36,4,4]], code-capacity depolarizing Z-memory.
# Reproduction: both noise points, 2000 shots each. Wall < 5 min on one core
# (nice -n 10). Reports land in qldpc-dec/scratch/gatec/.
set -euo pipefail
cd "$(dirname "$0")/.."   # qldpc-dec/src
PY="${PY:-../../.venv/bin/python}"

nice -n 10 "$PY" -m gatec_exact.run_gatec --p 1e-2 --shots 2000 --csv
echo
nice -n 10 "$PY" -m gatec_exact.run_gatec --p 3e-2 --shots 2000 --csv
echo
echo "reports: $(cd ../scratch/gatec && pwd)"
