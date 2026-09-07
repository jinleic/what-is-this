#!/bin/bash
# Kobon n=12, target 39: launch the supervised 15-cube cover run (v2).
#
# Why this exists: the 2026-08-16 attempt ran every cube with proof logging
# enabled, filled the disk with 13-14 GB partial DRATs, and was killed.  The
# 2026-09-03 pass separates DECIDING from CERTIFYING (proof-free verdicts
# first).  This file is now only an exec wrapper; every safety property lives
# in supervise.py, which is hash-pinned by manifest.json (schema /2):
#   * permanent nonblocking fcntl campaign lock (.campaign.lock, never
#     deleted) held for the supervisor's full lifetime;
#   * byte-for-byte manifest pinning of the 15 CNFs, engine, this wrapper,
#     the supervisor implementation, and the exact resolved Kissat binary;
#     This binds existing bytes together; it does not prove the engine
#     emitted those CNFs.  Byte-regeneration/provenance audit remains required;
#   * a prelaunch manifest digest record appended to launches.jsonl;
#   * manifest-keyed v2 ledger + log namespace (the legacy five-column
#     verdicts.tsv and logs/ stay untouched, discovery-only);
#   * one process group per Kissat, terminated and reaped on INT/TERM/exit;
#   * bounded concurrency, resumable exact-manifest SAT/UNSAT rows, retries
#     for nonterminal rows, refusal on unknown/duplicate terminal rows;
#   * all-UNSAT declared only with exactly one terminal row per named cube.
# Evidence boundary: SAT is a Boolean candidate only; proofless all-UNSAT is
# discovery evidence only.  DRAT emission is a separate pass.  See
# supervise.py, make_manifest.py, and manifest.json.
set -eu
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec python3 "$HERE/supervise.py" "$@"
