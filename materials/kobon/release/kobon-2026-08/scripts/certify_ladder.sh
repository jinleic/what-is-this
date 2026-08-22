#!/bin/bash
# Certify a broad-convention ladder rung end to end:
#   build CNF -> kissat --unsat with DRAT proof -> drat-trim -w
# Appends one JSON line per rung to scratch/kobon/ladder_certificates.jsonl.
# Usage: certify_ladder.sh N T
set -u
W=/Users/jinleic/jinleic-workspace
PY=$W/scratch/kobon-audit/venv/bin/python
KISSAT=$W/scratch/kobon-audit/tools/kissat/build/kissat
DRAT=$W/scratch/kobon-audit/tools/drat-trim/drat-trim
N=$1; T=$2
CNF=$W/scratch/kobon/ladder_n${N}_t${T}.cnf
PROOF=$W/scratch/kobon/ladder_n${N}_t${T}.drat
LOG=$W/scratch/kobon/ladder_n${N}_t${T}.dratcheck.log
LEDGER=$W/scratch/kobon/ladder_certificates.jsonl

"$PY" - "$N" "$T" "$CNF" <<'EOF'
import sys
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/math/kobon')
import engine
n, T, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
cnf, pool = engine.build_model(n, T)
tc = engine.add_triangle_crossing_indicators(cnf, pool, n)
engine.add_face_bound(cnf, pool, n, T, crossing_lits=tc, per_line=True)
engine.add_exact_selection(cnf, pool, n, T)
cnf.to_file(out)
print(f"{cnf.nv} {len(cnf.clauses)}")
EOF
DIMS=$(head -1 "$CNF" | awk '{print $3" "$4}')

t0=$(date +%s)
"$KISSAT" --unsat --seed=0 "$CNF" "$PROOF" > /dev/null 2>&1
rc=$?
t1=$(date +%s)
if [ "$rc" -ne 20 ]; then
  echo "{\"n\":$N,\"T\":$T,\"solver_rc\":$rc,\"verdict\":\"NOT-UNSAT\"}" >> "$LEDGER"
  echo "rc=$rc (not UNSAT) — no certificate"
  exit 1
fi
"$DRAT" "$CNF" "$PROOF" -w > "$LOG" 2>&1
t2=$(date +%s)
status=$(grep -c 's VERIFIED' "$LOG")
psize=$(stat -f%z "$PROOF")
echo "{\"n\":$N,\"T\":$T,\"vars_clauses\":\"$DIMS\",\"solve_s\":$((t1-t0)),\"verify_s\":$((t2-t1)),\"proof_bytes\":$psize,\"drat_verified\":$status}" >> "$LEDGER"
echo "n=$N T=$T solve=$((t1-t0))s verify=$((t2-t1))s proof=${psize}B VERIFIED=$status"
