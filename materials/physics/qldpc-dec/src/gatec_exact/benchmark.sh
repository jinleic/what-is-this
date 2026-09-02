#!/bin/bash
# Gate C benchmark: exact degenerate ML vs harness decoders on BB[[36,4,4]]
# under code-capacity depolarizing noise. Set BASIS=X or BASIS=Z.
# Reproduction: both noise points, 2000 shots each; reports land in scratch.
set -euo pipefail
cd "$(dirname "$0")/.."   # qldpc-dec/src
PY="${PY:-../../.venv/bin/python}"
BASIS="${BASIS:-Z}"

nice -n 10 "$PY" -m gatec_exact.run_gatec --basis "$BASIS" --p 1e-2 --shots 2000 --csv
echo
nice -n 10 "$PY" -m gatec_exact.run_gatec --basis "$BASIS" --p 3e-2 --shots 2000 --csv
echo
echo "reports: $(cd ../scratch/gatec && pwd)"
