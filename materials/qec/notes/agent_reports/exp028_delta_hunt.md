NEGATIVE — no reversal among **247 perturbations across 3 parents** on the
$(6,6)$ lattice; the 90-minute wall budget reached neither $(6,4)$ nor $(4,6)$.

## Evidence

- The full pre-registered budgeted run sampled 300 translation-canonical, connected weight-3/weight-3 parents (100 per target lattice), then processed the adversarial delta>0 stratum first. The fixed wall budget reached 3 parents on (6,6): 247 PBB perturbations received both exact distance solves, and all 247 parent/PBB status pairs were `OPTIMAL`. No candidate had `d(PBB) > d(parent)`. There were 0 preliminary reversals and therefore no hit to double-verify.
- Exact per-lattice solver coverage: (6,6): 3 parents, 247 perturbations, 247/247 both statuses `OPTIMAL`; (6,4): 0 parents, 0 perturbations; (4,6): 0 parents, 0 perturbations. The zero coverage on the latter two lattices is explicit and must not be read as evidence about them.
- Enumeration accounting on the three reached (6,6) parents: 2,324 commuting supports of total perturbation weight <=3, 844 translation-unique perturbations outside the catalogue orbit, 454 retained under per-parent caps, and 247 delta>0 positive-k candidates solver-tested before the budget-reserve stop.
- The observed parent distances for the three reached parents were exact(OPTIMAL): 2, 6, and 2. Every tested PBB distance was exact(OPTIMAL), with complete witnesses and per-sector statuses stored in the JSON.
- This is a budgeted search, not an exhaustive lattice-wide result.

## Catalogue perturbation-weight audit

The belief that the catalogue used perturbation weight at most 2 is false. Across all 368 rows, `len(C_terms)+len(D_terms)` has histogram: weight 2: 62; 3: 45; 4: 212; 5: 23; 6: 9; 7: 2; 8: 5; 9: 5; 10: 3; 11: 1; 13: 1. Maximum stored C weight is 6, maximum D weight is 7, and maximum combined weight is 13. Thus weight 3 itself is not new to the catalogue. The search remains outside the catalogue because each tested `(A,B,C,D)` safe translation orbit was absent from all catalogue rows; additionally, (6,4) and (4,6) have no catalogue rows, though the wall budget did not reach them.

## Protocol

- Seed 20260812; double-verification seed reserved as 20260813.
- Lattices `(ell,m)=(6,6),(6,4),(4,6)`; general translation-canonical weight-3 supports for A and B; connected Tanner graph; `4 <= k <= 16`.
- Perturbations satisfy `commutation_defect(A,B,C,D)==0`, have total weight <=3, and have at least one stabilizer row with nonzero X-sector and nonzero Z-sector support. Up to 200 per parent were selected by deterministic reservoir sampling after translation and catalogue-orbit deduplication.
- Corollary 1 rules out delta=0 reversals, so solver budget was assigned to delta>0, positive-k candidates first.
- `exact_distance_css(parent)` and `exact_distance_symplectic(PBB)`, 240 s limits, one exact solve at a time, 8 CP-SAT workers. `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=VECLIB_MAXIMUM_THREADS=NUMEXPR_NUM_THREADS=1` before spawning.
- A hit required both statuses `OPTIMAL`, strict `d_PBB>d_parent`, then a fresh rebuild, independent rank/witness checks, and seed+1 solver reruns. No preliminary hit triggered that path.

## Reproduce

```bash
cd /Users/jinleic/jinleic-workspace/qec-codesign
PYTHONPATH=src .venv/bin/python experiments/exp028_delta_hunt.py
PYTHONPATH=src .venv/bin/python -m pytest -q tests/test_exp028_delta_hunt.py
```

Observed search wall time: 4,718.159 s. The experiment reserved 480 s inside the 90-minute budget for mandatory hit verification and stopped when insufficient budget remained for another 240 s candidate solve plus that reserve. External interactive desktop workload dominated the run queue and QoS-deprioritized the experiment; this reduced wall-budget coverage but does not affect exact CP-SAT statuses, witnesses, GF(2) ranks, or candidate counts.

Focused verification: 3 tests passed in 113.85 s. The artifact-routing guard was exercised twice with forced partial writes; both routed to `results/partial_runs/`, and the canonical artifact hash/state was unchanged. The clean completed budgeted run is stored at `results/processed/exp028_delta_hunt.json`.

## Caveats

- Coverage is only 3 of the 300 sampled parents and only on (6,6); it is not evidence of absence on (6,4) or (4,6).
- Parent processing was descending k within each uniformly sampled lattice set to prioritize the only strata observed to admit delta>0 low-weight children. This affects budgeted prefix coverage, not the validity of the 247 exact comparisons.
- The search uses the finite stated translation equivalence for deduplication; it does not claim classification under arbitrary code equivalence.
