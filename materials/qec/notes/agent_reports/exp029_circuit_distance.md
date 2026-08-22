POSITIVE — both circuits have certified **undecomposed-DEM-mechanism** lower
bounds of at least 4 and validated upper bounds of at most 12.  These are
mechanism-model circuit-distance bounds, not an exact physical-location circuit
distance.

## Evidence

The fault model is **flattened undecomposed DEM error mechanisms**: one positive-probability `error` instruction from `circuit.detector_error_model(decompose_errors=False)` counts as one fault. Probabilities are ignored. The nonzero value `p=0.001` only materializes the DEM structure, so these distance bounds do not depend on p.

| Circuit | schedule | mechanisms | max detector degree | certified LB | best UB | graphlike |
|---|---:|---:|---:|---:|---:|---:|
| CSS-BB [[144,12,12]] Gross | 7 | 67032 | 9 | >= 4 CERTIFIED | <= 12 (certified mechanism-level upper bound) | <= not returned |
| nonCSS-PBB [[144,12,12]] 12_6_0193 | 8 | 82800 | 14 | >= 4 CERTIFIED | <= 12 (certified mechanism-level upper bound) | <= not returned |

### Certified lower bounds

- **CSS-BB [[144,12,12]] Gross**: all 67032 mechanisms passed the weight-1 condition (a detector fires or no observable flips), and hashing exact detector bitsets found no equal-detector/different-observable weight-2 pair. This certifies distance >= 3 in the stated mechanism model.
  The streamed singles-vs-pairs XOR search also completed: 36235296 unique shared-detector pairs checked with no weight-3 logical, certifying distance >= 4.
- **nonCSS-PBB [[144,12,12]] 12_6_0193**: all 82800 mechanisms passed the weight-1 condition (a detector fires or no observable flips), and hashing exact detector bitsets found no equal-detector/different-observable weight-2 pair. This certifies distance >= 3 in the stated mechanism model.
  The streamed singles-vs-pairs XOR search also completed: 85714200 unique shared-detector pairs checked with no weight-3 logical, certifying distance >= 4.

### Upper bounds

- **CSS-BB [[144,12,12]] Gross**: best validated upper bound is **distance <= 12**, via explicit weight-12 data-logical DEM witness. Witness: `results/raw/exp029_witnesses.json#circuits/gross/explicit_data_logical`. Stim heuristic attempts returning a witness: 0/4; heuristic-search wall time 2702.557 s.
  Graphlike search returned no bound: `ValueError: Failed to find any graphlike logical errors.` after 0.039 s on the undecomposed DEM.
  Accepted witness: explicit initial-data X_ERROR mechanisms on a pure-X logical support; provenance: existing exact Gross code-distance certificate witness_X; observables flipped: [0, 1, 2, 3, 4]; detector XOR is empty and all twelve signatures were found in the flattened undecomposed DEM.
- **nonCSS-PBB [[144,12,12]] 12_6_0193**: best validated upper bound is **distance <= 12**, via explicit weight-12 data-logical DEM witness. Witness: `results/raw/exp029_witnesses.json#circuits/pbb/explicit_data_logical`. Stim heuristic attempts returning a witness: 0/4; heuristic-search wall time 2701.921 s.
  Graphlike search returned no bound: `ValueError: Failed to find any graphlike logical errors.` after 0.029 s on the undecomposed DEM.
  Accepted witness: explicit initial-data X_ERROR mechanisms on a pure-X logical support; provenance: EXP-029 exact pure-X sector CP-SAT search: weight 12, all nontrivial sectors OPTIMAL or INFEASIBLE; observables flipped: [6, 9]; detector XOR is empty and all twelve signatures were found in the flattened undecomposed DEM.

### Gross versus PBB and hook implication

The best upper bounds are tied at 12. Gross has depth 7, weight-6 pure checks, while PBB has depth 8, weight-8 checks and 72 mixed checks. Thus the PBB circuit admits longer, mixed-sector single-mechanism hook propagation, but the reported comparison only shows whether the available witnesses exploit that structure; an upper-bound tie or gap is not by itself a proof of exact circuit distance. The certified lower bounds rule out one- and two-mechanism logical hooks in both circuits.

## Reproduction

From the repository root:

```bash
PYTHONPATH=src .venv/bin/python experiments/exp029_circuit_distance.py --phase all --circuit all
```

This command enforces a 45-minute total budget for the heuristic Stim search on each circuit, a separate 15-minute guard for graphlike search, a 45-minute weight-3 budget, and an 8 GiB RSS guard. It writes canonical artifacts only after one clean run has covered both circuits and all stages; subset runs go under `results/partial_runs/`.

## Runtime and caveats

- Recorded summed per-circuit wall time: **5515.039 s**.
- Stim version: `1.16.0`; Python: `3.13.14`.
- `decompose_errors=False` is used for every certified lower-bound signature and for the graphlike search. The graphlike result explicitly skips mechanisms with detector degree >2; it is an upper bound from that subset, not a lower bound and not an exact distance.
- `Circuit.search_for_undetectable_logical_errors` is heuristic under the recorded caps. Every returned entry was mapped back to an undecomposed DEM signature and XOR-validated before being accepted as a mechanism-level upper bound.
- No Monte Carlo was performed (shots=0, fails=0, no rate and no Clopper-Pearson interval).
- A budget-limited search is not exhaustive unless its record explicitly says completed/certified.
- Shared external desktop load materially inflated solver wall times. Exhaustive counts and logical-witness validity are unaffected; timeouts are wall-budget outcomes.
- Partial-run guard verified empirically: `--phase lower --circuit gross --w3-budget-s 0.001` wrote `results/partial_runs/exp029_lower_gross_state.json` and left the canonical raw state/witness SHA-256 digests unchanged.
- Recorded invocations for this artifact: finalize/all (canonical, 0.0 s); finalize/all (canonical, 0.0 s).
- Execution history of the recorded numbers: the first `--phase all --circuit all` invocation computed the Gross lower and upper stages plus the PBB lower stage, then aborted because the PBB pure-X CP-SAT witness returned a different (equally valid) translate of the weight-12 orbit than the hard-coded support; the code now adopts the freshly certified support, and `--phase finalize --circuit all` completed the PBB upper stage and wrote every canonical artifact. Both circuits were rebuilt from scratch in each invocation and their schedule/circuit/DEM SHA-256 digests matched.
