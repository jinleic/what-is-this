# Kobon K_gen(10) certification audit — 2026-08-15

Question: does the ~67 GB DRAT family in
`/Users/jinleic/Downloads/files-to-transfer-math-phys/workspace-files/scratch/kobon/` verify?

## Verdict

**No — the family as shipped cannot certify the upper bound.** 2 of the 11 cube proofs are
truncated mid-write (SIGHUP killed the solver fleet before completion) and contain no empty-clause
addition, so they are unverifiable in principle, not merely unchecked:

| proof | size | records | literals | empty clause |
|---|---|---|---|---|
| `proof_n10_t26_c0.drat` | 8.83 GB | — | — | **ABSENT** (ends mid-clause, `… 80 08 82 07`) |
| `proof_n10_t26_par_conc.drat` | 10.08 GB | 203 M | 4.71 B | **ABSENT** (ends mid-varint, `… 05 93`) |
| other 9 n=10 proofs + n=9 | 1.8–7.9 GB | 44–161 M | 0.8–3.7 B | present, exactly at EOF |

Method: binary-DRAT byte census. In the format, 0x00 occurs only as clause terminator; the pattern
`00 61 00` (or leading `61 00`) detects every empty-clause addition. Truncation therefore existed
at the source machine (bundle used APFS clonefile, size-verified copies).

Per `CubeCoverAudit.md` (cover correct, all 11 cubes required), the missing cubes leave open:
- **c0**: all pairs cross ∧ concurrency in the D10 orbit of triple (0,1,2);
- **par_conc**: parallel pair {0,1} ∧ some concurrency.

## What does hold

- **Lower bound `K_gen(10) ≥ 25`: re-verified** here via `engine.verify_n10_lower_bound()` —
  exact rational separating-axis check, 10 lines / 25 distinct triangles → `True, ok`.
- **CNF provenance: sealed.** All 12 shipped `kobon_n10_t26_proof*.cnf` regenerate
  **byte-identically** from `engine.py::dump_instances(10, 26)` (fresh venv, python-sat 1.9.dev14).
- **Positive control: `proof_n9_t22.drat` vs `kobon_n9_t22_exact.cnf` → `s VERIFIED`**
  (drat-trim, 386 s, 121,436,534 resolution steps — matches the session's recorded ~120 M / 444 s).
  So K(9) ≤ 21 is certified, and the toolchain works end-to-end.
- `report.md` §4.2 now records this incomplete certificate state explicitly;
  no upper-bound conclusion is claimed while the two replacement proofs run.

All relative paths below resolve against the working directory
`~/jinleic-workspace/scratch/kobon-audit/`; this document is canonical at
`math/kobon/AUDIT.md`. The canonical engine and report live beside it in
`math/kobon/`.

## In flight (supervised, survive this session)

- `kobon-drat-fleet` (hub process): drat-trim over the 9 complete n=10 proofs, 4-wide,
  biggest-first. Verdicts append to `results.tsv`; per-cube logs in `logs/<cube>.log`
  (`s VERIFIED` / `s NOT VERIFIED`, wall time, max RSS). Expected wall ≈ 4–5 h total.
- `kobon-solve-c0`, `kobon-solve-parconc` (hub processes): kissat 4.0.4 re-solving the two
  truncated cubes from the byte-identical CNFs; binary DRAT streaming to `resolve/`.
  Completion time unbounded (originals died at 8.8/10.1 GB of proof). On UNSAT exit (code 20),
  check with: `tools/drat-trim/drat-trim <cnf> resolve/<proof> -w`.
- Seed racers `kobon-{c0,parconc}-s{1..4}` were stopped 2026-08-15 with user
  approval (memory: each kissat historically peaks at 60–95 GB RSS; the machine
  has 96 GB, no swap). Their partial `resolve/*_sN.drat` files are inert.

## Split family (primary path for the two missing cubes)

`split/kobon_n10_t26_split_{c0,par_conc}_{0..7}.cnf` — verified byte-for-byte
equal to the corresponding shipped cube CNF plus exactly **3 unit clauses**
(header clause count +3):

- `c0` splits on X-vars **424, 586, 740** = X(3,4,5), X(5,6,7), X(7,8,9);
- `par_conc` splits on **340, 514, 740** = X(2,3,4), X(4,6,7), X(7,8,9);
- subcube `k` uses polarity `+v` iff bit `b` of `k` is set (bit order as listed).

Soundness: the 8 polarity patterns over 3 Boolean variables form a tautological
case split — any model of the cube CNF satisfies exactly one pattern — so
`s VERIFIED` UNSAT proofs for all 8 subcubes certify UNSAT of the shipped cube
CNF. No symmetry argument is used or needed.

Runs: `ksplit-{c0,parconc}-{k}` (kissat `--unsat`, detached), proofs to
`resolve/split_{cube}_{k}.drat`; launched in waves as memory frees
(drat-trim fleet peaks ~51 GB). Verification: `drat-trim <split cnf> <drat> -w`
per subcube. `kobon-memwatch` logs solver RSS to `memwatch.log` every 120 s.

## Progress log (2026-08-15, evening)

- **`c0` cube: CLOSED.** All 8 split subcubes solved UNSAT by kissat
  (exit 20; 21–35 min each, proofs 0.9–1.4 GB) and all 8 DRATs independently
  verified by drat-trim backward checking (`s VERIFIED`; 18–49 min each;
  ledger `splits.tsv`). With the byte-verified tautological split, this
  certifies UNSAT of the shipped `kobon_n10_t26_proof_c0.cnf`. The redundant
  whole-cube re-solve `kobon-solve-c0` was stopped after certification
  (5h54m, superseded).
- **Fleet ledger repair.** drat-trim emits `\n\r` line endings; the fleet's
  `grep -E '^s '` therefore recorded `NO-VERDICT` for finished cubes. The
  logs contain real verdicts; CR-tolerant extraction
  (`grep -aE '^\r?s '`) shows `s VERIFIED` for **c2, c5, c6, c7, simple**
  (5/9). `run_checks.sh` is patched for future runs; `results.tsv` rows
  remain raw, the corrected verdicts live in the per-cube logs and below.
- `par_conc` cube: 6/8 subcubes UNSAT (0,1,2,3,5,7), 2 solving (4,6);
  verified so far: 1 (`s VERIFIED`, 92.5 min); 0,2,3,5,7 verifying.
- **Fleet: COMPLETE — closure criterion 1 met.** All 9 shipped cube proofs
  independently re-verified: `s VERIFIED` for `simple`, `par_noc`, `c1`–`c7`
  (0.6–4.0 h each, ≈23.4 CPU-h total; peak RSS 2.8–9.6 GB). Corrected ledger:
  `results_verified.tsv` (regenerated from the per-cube logs; the raw
  `results.tsv` predates the CR fix).
- **`par_conc` cube: CLOSED.** All 8 split subcubes solved UNSAT by kissat
  (exit 20; 58 min–5 h each, proofs 1.9–6.3 GB) and all 8 DRATs independently
  verified (`s VERIFIED`; 42–110 min each). This certifies UNSAT of the
  shipped `kobon_n10_t26_proof_par_conc.cnf`. The backup whole-cube re-solve
  `kobon-solve-parconc` was stopped after certification (8h04m, superseded);
  `kobon-memwatch` retired.

## VERDICT — 2026-08-15, night: ALL CLOSURE CRITERIA MET

1. 9/9 shipped cube proofs `s VERIFIED` (`results_verified.tsv`).
2. `c0` and `par_conc` both closed by 8/8 `s VERIFIED` split subproofs each
   (`splits.tsv`; splits byte-verified as parent CNF + 3 units, a
   tautological case split — no symmetry argument used).
3. Lower bound `K_gen(10) ≥ 25` re-verified from the canonical engine.

**Hence `K_gen(10) = 25` stands on a fully machine-checked certificate
chain**: 25 is achieved by an explicit integer arrangement, and no 10-line
configuration under the broad convention (parallels, multiple points,
crossed triangles allowed) admits 26 pairwise interior-disjoint triangles.
Total independent verification ≈ 23.4 CPU-h (shipped proofs) + 16.4 CPU-h
(split proofs); replacement proof mass 44 GB across 16 DRATs.

## Closure criteria for K_gen(10) = 25

1. All 9 fleet rows in `results.tsv` read `s VERIFIED`;
2. for each of `c0` and `par_conc`, EITHER its whole-cube re-solve exits 20
   (UNSAT) with a `drat-trim`-verified fresh DRAT, OR all 8 of its split
   subcubes are UNSAT with `drat-trim`-verified DRATs (tautological case
   split — see above);
3. then upper bound ≤ 25 is certified over the audited cover, and with the verified lower
   certificate, `K_gen(10) = 25` under the report's broad convention.

Tooling: drat-trim @ github master (built 2026-08-15), kissat 4.0.4, both in `tools/`.
