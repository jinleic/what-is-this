#!/bin/sh
# Run one census stratum with implementation A and diff canonically against
# the published file. Usage: run_stratum.sh <n> <e> <expected_iso_classes> [workers]
set -eu
n="$1"; e="$2"; expect="$3"; workers="${4:-8}"
cd "$(dirname "$0")/.."
PY=/Users/jinleic/jinleic-workspace/math/.venv/bin/python
$PY src/glue_census.py --n "$n" --e "$e" --workers "$workers"
labelg "outA/r45_${n}_${e}.g6" 2>/dev/null | sort -u > "outA/r45_${n}_${e}.canon"
got=$(wc -l < "outA/r45_${n}_${e}.canon" | tr -d ' ')
pub="data/r45extreme/r45${n}.${e}.g6"
if [ -f "$pub" ]; then
  labelg "$pub" 2>/dev/null | sort -u > "outA/r45_${n}_${e}.pub.canon"
  if diff -q "outA/r45_${n}_${e}.canon" "outA/r45_${n}_${e}.pub.canon" >/dev/null; then
    echo "== n=$n e=$e: $got iso classes (expect $expect) — canonical sets EQUAL to published"
  else
    echo "!! n=$n e=$e: $got iso classes (expect $expect) — canonical sets DIFFER"
    exit 1
  fi
else
  echo "?? n=$n e=$e: $got iso classes (expect $expect) — no published file to diff"
fi
