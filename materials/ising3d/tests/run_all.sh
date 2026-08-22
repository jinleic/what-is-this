#!/usr/bin/env bash
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="$ROOT/.venv/bin/python"

if [[ ! -x "$PYTHON" ]]; then
  printf 'ERROR: interpreter is not executable: %s\n' "$PYTHON" >&2
  exit 2
fi

shopt -s nullglob
tests=("$ROOT"/tests/test_*.py)

passed=0
failed=0
for test_path in "${tests[@]}"; do
  test_name="tests/${test_path##*/}"
  printf '\n=== %s ===\n' "$test_name"
  if (cd "$ROOT" && "$PYTHON" "$test_path"); then
    printf 'SUMMARY %s: PASS\n' "$test_name"
    ((passed += 1))
  else
    exit_code=$?
    printf 'SUMMARY %s: FAIL (exit %d)\n' "$test_name" "$exit_code"
    ((failed += 1))
  fi
done

total=${#tests[@]}
printf '\nFINAL: %d total, %d passed, %d failed\n' "$total" "$passed" "$failed"
((failed == 0))
