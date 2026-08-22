BREAKTHROUGH_CANDIDATE — Theorem C6 basis-independence proved computationally for all 14 catalogue PBB members.

# EXP-023 light generating sets

## What was decided

For each code, `V_w` is the GF(2) span of all stabilizers of symplectic weight at most `w`, and `w*` is the least `w` with `rank(V_w) = n - k = 132`. The decisive question per PBB member is whether any stabilizer of weight at most 7 has nonzero X-part. It is settled exhaustively, not by a single opaque solver call:

- Every stabilizer of a PBB code is `(u H_X | u P + s)` with `s` ranging over the pure-Z subgroup `S_Z`, so a weight-<=7 element with nonzero X-part exists **iff** some nonzero codeword `x = u H_X` of `rowspace(H_X)` with `|x| <= 7` admits an `s` leaving at most `7 - |x|` qubits of Z-support outside `supp(x)`.
- Stage 1 enumerates **all** such codewords `x` with CP-SAT; the `OPTIMAL` status after exhaustive enumeration is the completeness certificate. Stage 2 decides each candidate by exact GF(2) coset arithmetic, with no solver involved.
- `r_Zsub`, the exact rank of the full pure-Z subgroup `rowspace(H_Z) + {u[C D] : u H_X = 0}`, is computed by plain linear algebra and cross-checked against the direct definition `{c H : (c H)_X = 0}`.

The two proof routes recorded per member are:

- `gate_plus_structural`: no mixed light element exists, hence `V_7` lies inside the pure-Z subgroup and `rank(V_7) <= r_Zsub < 132`. The published weight-6 Z rows span exactly `r_Zsub`, so this upper bound meets a matching lower bound and `rank(V_6) = rank(V_7) = r_Zsub` is exact.
- `exact_closure`: the iterative closure loop terminated on its own, either spanning the full rowspace or returning CP-SAT `INFEASIBLE` against a complete basis of physical annihilator functionals.

Both routes prove the same statement: every generating set of the stabilizer group contains a generator of symplectic weight at least 8, so by Proposition C2 every one-ancilla syndrome-extraction round has depth at least 8 for that code, independently of the measured basis.

| code | canonical max | gate w=7 | r_Zsub | rank V6 | rank V7 | proof route | w* | depth >= 8 basis-independent? |
|---|---:|---|---:|---|---|---|---|---|
| CSS Gross [[144,12,12]] (control) | 6 | n/a | n/a | w=5: 0 exact | w=6: 132 exact | exact_closure | 6 | no (w*=6) |
| 12_6_0187 | 10 | INFEASIBLE | 66 | 66 exact | 66 exact | gate_plus_structural | in [8, 10] | yes |
| 12_6_0188 | 9 | INFEASIBLE | 66 | 66 exact | 66 exact | gate_plus_structural | in [8, 9] | yes |
| 12_6_0189 | 10 | INFEASIBLE | 66 | 66 exact | 66 exact | gate_plus_structural | in [8, 10] | yes |
| 12_6_0190 | 9 | INFEASIBLE | 66 | 66 exact | 66 exact | gate_plus_structural | in [8, 9] | yes |
| 12_6_0191 | 10 | INFEASIBLE | 66 | 66 exact | 66 exact | gate_plus_structural | in [8, 10] | yes |
| 12_6_0192 | 10 | INFEASIBLE | 66 | 66 exact | 66 exact | gate_plus_structural | in [8, 10] | yes |
| 12_6_0193 | 8 | INFEASIBLE | 66 | 66 exact | 66 exact | gate_plus_structural | 8 | yes |
| 12_6_0194 | 8 | INFEASIBLE | 66 | 66 exact | 66 exact | gate_plus_structural | 8 | yes |
| catalogue_row_0036_b261ff5e77ae456d | 10 | INFEASIBLE | 66 | 66 exact | 66 exact | gate_plus_structural | in [8, 10] | yes |
| catalogue_row_0037_40bae617111e7148 | 10 | INFEASIBLE | 66 | 66 exact | 66 exact | gate_plus_structural | in [8, 10] | yes |
| catalogue_row_0039_1219913c49eb1daf | 8 | INFEASIBLE | 66 | 66 exact | 66 exact | gate_plus_structural | 8 | yes |
| catalogue_row_0040_4a77a996000a5d40 | 10 | INFEASIBLE | 66 | 66 exact | 66 exact | gate_plus_structural | in [8, 10] | yes |
| catalogue_row_0041_08166cccb951f700 | 9 | INFEASIBLE | 66 | 66 exact | 66 exact | gate_plus_structural | in [8, 9] | yes |
| catalogue_row_0043_393b2d4a60427543 | 8 | INFEASIBLE | 66 | 66 exact | 66 exact | gate_plus_structural | 8 | yes |

`rank V6`/`rank V7` are the secondary rank-bookkeeping tier; an `UNDECIDED` there never blocks the theorem, which rests on the decisive gate plus `r_Zsub`. The Gross control columns are `w=5` and `w=6`.

`w*` is exact where a single value is shown; where an interval is shown, the lower bound is the proved one (no generating set can stay under weight 8) and the upper bound is the published generating set's own maximum weight, which this experiment did not try to improve.

## Independent verification and controls

Coefficients live on an explicitly recorded 132-row independent subset `H_B` of the published rows, so `c -> c H_B` is a bijection onto the stabilizer rowspace and no rank or closure statement can be corrupted by the 12-dimensional kernel of the full 144-row map. Every closure functional is physical: an `h` in `F_2^288` annihilating the already-found stabilizers, pulled back as `f = H_B h`. Every witness was recomputed as `v = c H_B` with NumPy GF(2) arithmetic, weighed with both a local support count and the project's `symplectic_weight`, confirmed to lie in `rowspace(H)` by rank comparison, and ranked in 288-bit physical space by three independent implementations.

Because every scientific claim here rests on a negative (infeasible) answer, each code also runs controls that must come out positive, and a failure aborts the run before any canonical artifact is written:

- CP-SAT positive control: a nonzero stabilizer of weight at most the canonical max exists, with a physically verified witness.
- CP-SAT positive control (PBB): a stabilizer of weight at most the canonical max with nonzero X-part exists, with a physically verified witness. This exercises exactly the query shape that returns `INFEASIBLE` at w=7.
- Decomposition controls: every published weight-<=7 X-check row appears in the stage-1 enumeration; the stage-2 coset test recovers a published mixed generator at its own budget; an impossible (negative) budget is rejected.
- Cross-encoding check on the benchmark member: the monolithic single-query CP-SAT formulation was run independently and its answer must agree with the decomposition.

A separate cross-check (`experiments/exp023_verify_artifact.py`, output `results/raw/exp023_independent_verification.json`) does not import the producing script: it rebuilds every code from the catalogue, recomputes `rank(H)` and `r_Zsub`, re-checks every stored coefficient witness, brute-forces all published-row combinations of coefficient weight at most 3, and independently repeats the exact GF(2) coset test for every published light X row. This is a bounded independent cross-check; the all-coefficient-space completeness proof remains the stage-1 CP-SAT `OPTIMAL` enumeration.

`independent_verification_flags_all_true`: True.

## Timeouts and statistics

Timeouts: none.

Decisive enumeration wall time across the 14 PBB members: min 33.1 s, median 109.8 s, max 226.3 s.
Cross-encoding check on `12_6_0193`: the monolithic single-query CP-SAT formulation returned `INFEASIBLE` in 558.4 s against the decomposition's `INFEASIBLE` in 167.4 s; methods agree: True.

This is an exact finite feasibility experiment, not Monte Carlo: shots, failures, and Clopper-Pearson intervals do not apply. Every reported outcome is `OPTIMAL`, `INFEASIBLE`, or an explicitly labelled `UNDECIDED`.

## Reproduction

```bash
cd /Users/jinleic/jinleic-workspace/qec-codesign
PYTHONPATH=src .venv/bin/python experiments/exp023_light_gensets.py
PYTHONPATH=src .venv/bin/python -m pytest -q tests/test_exp023_light_gensets.py
PYTHONPATH=src .venv/bin/python experiments/exp023_verify_artifact.py
```

Observed end-to-end wall time: 2263.442 s. The decisive enumeration was capped at 900 s per code (single worker, as CP-SAT enumeration requires) and secondary rank-closure calls at 600 s with 8 workers. Base seed 20260812; every call seed is recorded in the JSON artifacts. Timing was measured on a heavily loaded shared machine (system load average above 100 during the run, dominated by an external desktop workload) and is not a benchmark; the SAT/UNSAT outcomes themselves are timing-independent.

## Artifacts and caveats

- Processed results: `results/processed/exp023_light_gensets.json`.
- Coefficient witnesses, physical closure functionals, pure-Z subgroup bases, and any explicit light generating set: `results/raw/exp023_light_elements.json`.
- A non-null `w_star_if_determined` is exact; otherwise `w_star_interval` is a rigorous bracket and this pre-registered w=6,7 threshold experiment claims nothing sharper.
- Catalogue distance metadata is used only for input selection; EXP-023 makes no distance claim.
- Independent verification pass: `results/raw/exp023_independent_verification.json`.
- No PBB member reached `rank(V_7) = 132`, so no depth-7 PBB generating set exists to save. The one explicit light generating set stored is the Gross control's own 132 weight-6 rows, which is what `w* = 6` means for it.
- The proved statement is scoped to one-ancilla syndrome extraction via Proposition C2; it says nothing about multi-ancilla or measurement-gadget schemes.
- This supersedes the partial attempt in `experiments/exp022_basis_independence.py`, which measured a weaker quantity (`w_mix` on a single element, and an X-part-only rank) and left the question open with unproven-optimal solver values.
- Write-integrity incident (`notes/failed_routes.md` FR-012 pattern): an early `--only gross` smoke run of this script wrote the canonical paths. Output routing now goes through `qec_research.artifacts.canonical_route`, so only a run covering the full declared scope can touch them; partial runs land in `results/partial_runs/` with a timestamped suffix. The guard was verified by re-running `--only gross` and confirming the canonical files' checksums were unchanged, and the canonical artifacts here were regenerated from scratch by the full 15-code run.
