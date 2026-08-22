# Experiment ledger

Snapshot time: **2026-08-22T01:45:58Z** (2026-08-21 local campaign day).
Statuses below come from the files and live process command lines observed at
that snapshot.  A live status can change after this document is written.

Evidence labels have strict meanings:

- **MACHINE-CHECKED** — an exact verifier completed, or a fixed-arrangement SAT
  optimum has `optimality_proved: true`.
- **DISCOVERY-ONLY** — a solver/model/probe result without an independently
  checked proof; useful for search, not a theorem.
- **PENDING** — in flight, queued, killed without a verdict, or no run found.
- A solver's `UNSATISFIABLE`/exit `20` is still discovery-only until its complete
  proof is independently verified.

DIMACS sizes below were read from the actual `p cnf` lines, not inferred from
filenames.  Byte counts are filesystem sizes at the snapshot.

## 1. Completed exact verification and packing experiments

| Run | Instance size / target | Settled result | Evidence | Transcript or artifact |
|---|---:|---|---|---|
| Maiorana fail-closed verifier | 14 lines, claimed 54 | All seven stages passed; final `K_gen_14_lower_bound=54` | MACHINE-CHECKED | `scratch/kobon/n14/maiorana_verify.log`; copied to `math/kobon/verification/maiorana_verify.log` |
| Maiorana all-solution sweep | 15 independent 14-line files, 54 each | 15/15 hashes, distinctness checks, and exact censuses passed | MACHINE-CHECKED | `scratch/kobon/n14/sweep_all_verify.log`; copied to `math/kobon/verification/sweep_all_verify.log` |
| Theorem-M finite checker, n=14 | 165 paired placements + 12 triad placements + controls | All obstruction legs pass | MACHINE-CHECKED finite leg | `scratch/kobon/n14/escape/appendix_independent_check.py`; no standalone `.log` exists on disk |
| Theorem-M finite checker, n=18 | 455 paired placements, 19 dihedral orbits, 16 triad placements + controls | All obstruction legs pass | MACHINE-CHECKED finite leg | `scratch/kobon/n18/appendix_check_n18.py`; no standalone `.log` exists on disk |
| Maiorana max-packing census | 15 fixed arrangements, n=14 | Every optimum is exactly 54; UNSAT at 55; all 15 have `optimality_proved=true` | MACHINE-CHECKED | `scratch/kobon/n14/maiorana/sol{1..15}_lines_rational.maxpacking.json` |
| Bader 53 max packing | fixed n=14 arrangement | optimum 53, proved | MACHINE-CHECKED | `scratch/kobon/discoveries/n14_bader53_verified.maxpacking.json` |
| Other exact 53 max packing | fixed n=14 arrangement | optimum 53, proved | MACHINE-CHECKED | `scratch/kobon/discoveries/n14_best_exact53.maxpacking.json` |
| Known n=11 packing | fixed n=11 arrangement | optimum 32, proved | MACHINE-CHECKED | `scratch/kobon/discoveries/n11_k32.maxpacking.json` |
| Known n=18 packing | fixed n=18 arrangement | optimum 93, proved | MACHINE-CHECKED | `scratch/kobon/discoveries/n18_known93.maxpacking.json` |
| Known n=20 packing | fixed n=20 arrangement | optimum 116, proved | MACHINE-CHECKED | `scratch/kobon/discoveries/n20_known116.maxpacking.json` |
| n=20 reconstruction | 20 lines; overlap graph reported 1,122 triples / 294,134 clashes | 116 selected and exact `verify_selection=True`; later max-packing run proves optimum 116 | MACHINE-CHECKED after later prober | `scratch/kobon/n20/known116_reconstruct.log`, `scratch/kobon/discoveries/n20_known116.json` |
| n=9 SAT-pipeline rehearsal | `p cnf 62517 379298`, target 21 | Kissat SAT; persisted model/rehearsal certificate | DISCOVERY-ONLY SAT model; exact witness artifact retained | `scratch/kobon/n9_rehearsal/n9_t21.kissat.stdout`, `scratch/kobon/n9_rehearsal/rehearsal_report.json`, `scratch/kobon/discoveries/n9_t21_satpipe_rehearsal.json` |

The obstruction scripts check the finite cycle enumerations.  Their hand-proof
interfaces and hypotheses are separately stated in
`scratch/kobon/n14/escape_obstruction.md` and
`scratch/kobon/n18/obstruction_n18.md`.

## 2. n=14 production CNFs and live runs

### Broad-pattern cubes and monoliths

| Instance | Actual DIMACS header | Bytes | Snapshot status | Log / proof artifact |
|---|---:|---:|---|---|
| `scratch/kobon/n14_monolith_t54.cnf` | `p cnf 1357615 7841316` | 240,994,798 | two Kissat runs (seeds 0,1) in flight; no verdict | `scratch/kobon/n14_monolith_t54.kissat0.log`, `scratch/kobon/n14_monolith_t54.kissat1.log`; no DRAT |
| `scratch/kobon/n14_cube_q0_t54.cnf` | `p cnf 350209 5824853` | 198,349,702 | Kissat seed 0 in flight; no verdict | `scratch/kobon/n14_cube_q0_t54.kissat0.log`; no DRAT |
| `scratch/kobon/n14_cube_q1_t54.cnf` | `p cnf 341539 5807521` | 197,994,509 | Kissat seed 0 and a CaDiCaL run in flight; no verdict | `scratch/kobon/n14_cube_q1_t54.kissat0.log`; CaDiCaL had no separate transcript found |
| `scratch/kobon/n14_cube_q2_t54.cnf` | `p cnf 332861 5790173` | 197,638,988 | Kissat seed 0 in flight; no verdict | `scratch/kobon/n14_cube_q2_t54.kissat0.log`; no DRAT |
| `scratch/kobon/n14_cube_q3paired_t54.cnf` | `p cnf 324175 5772809` | 197,283,139 | Kissat seed 0 in flight; no verdict | `scratch/kobon/n14_cube_q3paired_t54.kissat0.log`; no DRAT |
| `scratch/kobon/n14_cube_q3triad_t54.cnf` | `p cnf 324169 5772797` | 197,282,896 | Kissat seed 0 in flight; no verdict | `scratch/kobon/n14_cube_q3triad_t54.kissat0.log`; no DRAT |
| `scratch/kobon/n14/n14_monolith_t55.cnf` | `p cnf 1349454 7824997` | 240,627,409 | Kissat seed 0 in flight; no verdict | `scratch/kobon/n14_monolith_t55.kissat0.log`; no DRAT |
| `scratch/kobon/n14/n14_monolith_t56.cnf` | `p cnf 1341289 7808671` | 240,259,861 | Kissat seed 0 in flight; no verdict | `scratch/kobon/n14_monolith_t56.kissat0.log`; no DRAT |

Every row in this subtable is **PENDING / DISCOVERY-ONLY**.  In particular, the
analytically closed Q=3 branch does not turn an in-flight solver log into a
certificate, and the live monoliths have not settled.

### Planted-concurrency and position cubes

| Instance | Actual DIMACS header | Bytes | Snapshot status | Log / proof artifact |
|---|---:|---:|---|---|
| `scratch/kobon/n14/n14c-q0tp1.cnf` | `p cnf 1357615 7841407` | 240,995,268 | proof-producing Kissat in flight | `.cnf.kissat.log`; incomplete `.cnf.drat` |
| `scratch/kobon/n14/n14c-q0tp2s.cnf` | `p cnf 1357615 7841408` | 240,995,273 | discovery Kissat in flight | `scratch/kobon/n14/n14c-q0tp2s.kissat.log`; no DRAT |
| `scratch/kobon/n14/n14c-q1tp1d.cnf` | `p cnf 1357615 7841407` | 240,995,271 | proof-producing Kissat in flight | `.cnf.kissat.log`; incomplete `.cnf.drat` |
| `scratch/kobon/n14/n14c-q1tp1s.cnf` | `p cnf 1357615 7841407` | 240,995,270 | proof-producing Kissat in flight | `.cnf.kissat.log`; incomplete `.cnf.drat` |
| `scratch/kobon/n14/q2pos/q2pos-001.cnf` | `p cnf 332861 5790173` | 197,638,988 | Kissat seed 0 in flight | sibling `q2pos-001.kissat.log`; no DRAT |
| `scratch/kobon/n14/q2pos/q2pos-002.cnf` | `p cnf 332861 5790173` | 197,638,988 | Kissat seed 0 in flight | sibling `q2pos-002.kissat.log`; no DRAT |
| `scratch/kobon/n14/q3pos/q3pos-000.cnf` | `p cnf 324175 5772809` | 197,283,139 | Kissat seed 0 in flight | sibling `q3pos-000.kissat.log`; no DRAT |
| `scratch/kobon/n14/q3pos/q3pos-001.cnf` | `p cnf 324175 5772809` | 197,283,139 | Kissat seed 0 in flight | sibling `q3pos-001.kissat.log`; no DRAT |

All eight are **PENDING**.  A growing `.drat` is explicitly not a settled proof.
The remaining generated `q2pos-*.cnf` and `q3pos-*.cnf` files have no matching
run transcript at this snapshot unless listed above.

## 3. n=18 target-94 cube batch

| Instance | Actual DIMACS header | Bytes | Snapshot status | Log / artifact |
|---|---:|---:|---|---|
| `scratch/kobon/n18/n18_cube_q0_t94.cnf` | `p cnf 1708348 29932035` | 1,063,540,279 | Kissat in flight; no verdict | `scratch/kobon/n18/n18_cube_q0_t94.kissat0.log`; no DRAT |
| `scratch/kobon/n18/n18_cube_q1_t94.cnf` | `p cnf 1682318 29879983` | 1,062,343,488 | Kissat and CaDiCaL in flight; no verdict | `scratch/kobon/n18/q1_cadical.log`; no DRAT |
| `scratch/kobon/n18/n18_cube_q2_t94.cnf` | `p cnf 1656280 29827915` | 1,061,146,329 | Kissat seed 0 in flight; no verdict | `scratch/kobon/n18/n18_cube_q2_t94.kissat0.log`; no DRAT |
| `scratch/kobon/n18/n18_cube_q3paired_t94.cnf` | `p cnf 1630234 29775831` | 1,059,948,802 | held, not launched in batch 2; analytically obstructed branch | `scratch/kobon/n18/batch2_launch.md`; no run log |
| `scratch/kobon/n18/n18_cube_q3triad_t94.cnf` | `p cnf 1630228 29775819` | 1,059,948,532 | held, not launched in batch 2; analytically obstructed branch | `scratch/kobon/n18/batch2_launch.md`; no run log |

The first three are **PENDING / DISCOVERY-ONLY**; the last two are **PENDING as
solver experiments**, even though the corresponding mathematical Q=3 branch is
closed by the separately documented obstruction.

## 4. n=20 target-117 cube batch

| Instance | Actual DIMACS header | Bytes | Snapshot status | Log / artifact |
|---|---:|---:|---|---|
| `scratch/kobon/n20/n20_cube_q0_t117.cnf` | `p cnf 3455871 59195137` | 2,131,260,256 | no active run or transcript found | CNF only |
| `scratch/kobon/n20/n20_cube_q1_t117.cnf` | `p cnf 3414933 59113269` | 2,129,377,829 | Kissat and CaDiCaL in flight; no verdict | process command lines observed; no on-disk solver transcript found |
| `scratch/kobon/n20/n20_cube_q2_t117.cnf` | `p cnf 3373987 59031385` | 2,127,495,034 | no active run or transcript found | CNF only |
| `scratch/kobon/n20/n20_cube_q3paired_t117.cnf` | `p cnf 3333033 58949485` | 2,125,611,871 | no active run or transcript found | CNF only |
| `scratch/kobon/n20/n20_cube_q3triad_t117.cnf` | `p cnf 3333027 58949473` | 2,125,611,601 | no active run or transcript found | CNF only |
| `scratch/kobon/n20/n20_cube_q4four_t117.cnf` | `p cnf 3292071 58867569` | 2,123,728,340 | no active run or transcript found | CNF only |
| `scratch/kobon/n20/n20_cube_q4pairtriad_t117.cnf` | `p cnf 3292065 58867557` | 2,123,728,070 | no active run or transcript found | CNF only |

All seven are **PENDING**.  The n=20 Q=3 note is conditional and explicitly does
not prove the branch at target 117; `defect_lemma.md` records exact cycle-level
escape families.  The exploratory diagnostics `scratch/kobon/n20/explore1.log`
through `explore4.log` are settled **DISCOVERY-ONLY** geometry diagnostics, not
solver verdicts.

## 5. n=11 target-33 monolith

| Instance | Actual DIMACS header | Bytes | Snapshot status | Log / proof artifact |
|---|---:|---:|---|---|
| `scratch/kobon/n11_monolith_t33.cnf` | `p cnf 252622 1521724` | 44,481,369 | proof-producing Kissat seed 0, a second Kissat seed, and a PySAT/CaDiCaL solve were in flight; no verdict | `scratch/kobon/n11_monolith_t33.kissat.log`; growing `scratch/kobon/n11_monolith_t33.drat` |

This is **PENDING**.  The separate no-concurrency target-33 solve reported
UNSAT after 12,018,922 conflicts / 1,600.7 s but emitted no DRAT, so it remains
**DISCOVERY-ONLY** (`scratch/kobon/n11/section_draft.md`).

## 6. Savchuk-tooling probe experiments

### Three interrupted global runs

| Run | Instance size | Settled status | Log |
|---|---:|---|---|
| global missing-6 | `p cnf 98426 378658` | SIGTERM after about 2,523 s; no SAT/UNSAT verdict — PENDING | `scratch/kobon/n14_savchuk_global6.kissat.log` |
| exact-sym missing-6 | `p cnf 177014 522743` | SIGTERM after about 2,412 s; no verdict — PENDING | `scratch/kobon/n14_savchuk_global6_exact_sym.kissat.log` |
| seed/unit variant | no separate CNF header artifact found | SIGTERM after about 671 s; no verdict — PENDING | `scratch/kobon/n14_savchuk_global6-s1.kissat.log` |

### Bounded probe batch

All rows below are in `scratch/kobon/n14/work/probe_results.log`.  SAT controls
were subsequently straightened and exact-checked in `stack_validation.json`;
probe UNSAT rows have no DRAT and remain discovery-only.

| Probe | Header from log | Solver result |
|---|---:|---|
| `sanity-n5-perfect` | `p cnf 200 2200` | SAT, 0 s; exact stack control passed |
| `sanity-n7-parity` | `p cnf 672 10080` | UNSAT, 0 s; DISCOVERY-ONLY |
| `sanity-n6-m3` | `p cnf 450 6271` | UNSAT, 0 s; DISCOVERY-ONLY |
| `sanity-n9-perfect` | `p cnf 1584 30096` | SAT, 0 s; exact stack control passed |
| `n14-m6-line01` | `p cnf 7670 341303` | UNSAT, 89 s; DISCOVERY-ONLY |
| `n14-m33-lines1-2` | `p cnf 7670 341304` | UNSAT, 420 s; DISCOVERY-ONLY |
| `n14-m33-opposite` | `p cnf 7670 341304` | UNSAT, 251 s; DISCOVERY-ONLY |
| `n14-m222-lines1-3` | `p cnf 7670 341305` | UNSAT, 507 s; DISCOVERY-ONLY |
| `n14-m222-spread` | `p cnf 7670 341305` | UNSAT, 510 s; DISCOVERY-ONLY |
| `n14-m2211` | `p cnf 8138 414319` | UNSAT, 417 s; DISCOVERY-ONLY |
| `n14-m111111-lines1-6` | `p cnf 7670 341308` | UNKNOWN at 900 s |
| `n14-m111111-spread` | `p cnf 7670 341308` | UNKNOWN at 901 s |
| `n14-m8-lines1-8` | header not retained in the combined log | UNKNOWN at 3,600 s |
| `n14-m9-lines1-9` | header not retained in the combined log | UNKNOWN at 3,600 s |
| `n14-m9-spread` | header not retained in the combined log | UNKNOWN at 3,600 s |

## 7. Historical n=12 verification attempts retained on disk

These are included to prevent old files from being mistaken for certified
closure:

| Run family | Count | Settled status | Logs |
|---|---:|---|---|
| first DRAT verification attempt: `c0..c11`, `simple`, `par_noc`, `split_par_conc_0..7`, and even/decorated deep tiles | 28 | `s NOT VERIFIED` / `ERROR: no conflict`; not evidence | `scratch/kobon/n12/verify_attempt1_20260816/*.log` |
| odd deep `par_conc` tiles `0_1,0_3,0_5,0_7` | 4 | `s VERIFIED` for those four subinstances only | same directory |
| xz smoke controls + `parnoc_{1,3,5,7}_a2` | 6 | `s VERIFIED` | `scratch/kobon/n12/verify_xz/*.log` |
| current deep n=12 solver fleet | multiple CNFs under `deeper_c1_7/`, `deep_parnoc/`, `deep_simple/`, `deep_c0/`, `deep0/` | in flight at snapshot; no consolidated final verdict | live command lines; see `scratch/kobon/n12/closeout_draft.md` and `path_arbiter_rec.md` for the earlier ledger |

A verified subcube does not close `K_gen(12)`: the cover and all remaining
branches still need proof-checked closure.  The canonical `math/kobon/report.md`
and `AUDIT.md` separately record the already certified n=10 campaign and its
complete DRAT cover.
