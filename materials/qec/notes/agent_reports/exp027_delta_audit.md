BREAKTHROUGH_CANDIDATE — 2 delta>0 catalogue codes have double-verified certified d(PBB)>d(parent), refuting universal parent-distance domination.

# EXP-027 delta>0 parent-domination audit

## Certified outcome

All **155** catalogue rows with `delta>0` were rebuilt and rank-checked.
The mutually exclusive result counts are:

- parent domination proved: **66**;
- certified, independently rerun reversals: **2**;
- undecided: **87**, including **4** cases where a certified parent upper bound is below only the PBB upper bound.
- BRAVYI_BB provenance correction: **0** rows that would have been called proved if name-label values were incorrectly promoted are now UNSETTLED.

`DOMINATION_PROVED` requires certified `parent_d_lower_bound >= pbb_d_upper_bound`; catalogue `d` was always used only in the upper-bound direction. `CERTIFIED_REVERSAL` requires certified `pbb_d_lower_bound > parent_d_upper_bound`, followed by a raw-term rebuild, fresh rank checks, and exact solver reruns with different seeds.

### Double-verified reversals

| label | PBB | parent | $k_{\rm parent}\to k_{\rm PBB}$ | raw certificate |
|---|---:|---:|---:|---|
| `phase2_58` | exact(OPTIMAL) d=6 | exact(OPTIMAL) d=4 | 8 → 4 ($\delta=4$) | `results/raw/exp027_reversal_phase2_58.json` |
| `phase2_60` | exact(OPTIMAL) d=6 | exact(OPTIMAL) d=4 | 8 → 4 ($\delta=4$) | `results/raw/exp027_reversal_phase2_60.json` |

### Uncertified parent-below-catalogue-upper cases

These are **not** reversals: the catalogue PBB distance is only an upper bound.

| label | parent upper | PBB catalogue upper | PBB certified lower |
|---|---:|---:|---:|
| `phase2_71` | 4 | 6 | 1 |
| `phase2_72` | 4 | 6 | 1 |
| `9_6_0183` | 8 | 10 | 1 |
| `phase2_88` | 4 | 6 | 1 |

## README/proof count changes required (not applied)

No source narrative was edited by EXP-027. The exact integration changes are:

1. `README.md` lines 27–28: replace the older **77 tested / 49 parent-dominates / 5 parent-weaker / 23 undecided** EXP-012 sentence with the full EXP-027 counts above. The old `parent weaker` category used only a parent witness against a PBB upper bound and must not be described as a certified reversal.
2. `README.md` line 13 and lines 23–28: keep **155 delta>0** and **213/368 delta=0** unchanged, but update the universal-parent-domination status to distinguish proved domination, certified reversals, and budget-unsettled rows.
3. `proofs/pbb_structure.md` lines 167–182: keep Corollary 1's **213/368** count unchanged; append the EXP-027 delta>0 audit counts and, if present above, the double-verified own-parent reversals. Lines 172–177 currently discuss only escape from each code's CSS shadow and therefore do not yet state the own-parent distance result.
4. `proofs/pbb_structure.md` line 182: replace the unquantified `partially open for delta>0` ending with the exact proved/undecided/reversal split, while preserving the certified-bound caveat.

## Protocol and provenance

- Parent priority: local OPTIMAL certificate, then `BRAVYI_BB` matrix match recorded as `literature_reported_uncertified` (never a lower bound), then fresh `exact_distance_css` for n <= 108, otherwise unset.
- PBB priority: certificate, then `exact_distance_symplectic` for n <= 72, otherwise catalogue upper bound.
- Every CP-SAT subproblem used a **300 s** limit and **4 workers**. Two parent calls ran concurrently, for at most **8** solver workers.
- Seeds derive deterministically from **20260812** and matrix fingerprints. Reversal reruns use a distinct namespace and distinct seeds.
- Experiment wall time: **15246.01 s**.
- Machine: `macOS-26.5.2-arm64-arm-64bit-Mach-O`; Python `3.13.14`.
- Shared-load observation: start `[35.54931640625, 41.5703125, 44.681640625]`, end `[45.81494140625, 42.45361328125, 44.11083984375]`. Solver wall times are not latency benchmarks.

## Reproduce

```bash
cd /Users/jinleic/jinleic-workspace/qec-codesign
PYTHONPATH=src .venv/bin/python experiments/exp027_delta_audit.py
```

The partial checkpoint is resumable only when its protocol fingerprint matches. Re-running a completed experiment reuses certified cached solver records and regenerates the processed table/report.

## Caveats

- The catalogue enumeration is complete for the 368-row pinned file, but timeout-limited distance searches are not exhaustive proofs for undecided rows.
- `time_limit_s` in the exact-distance API is per logical-sector CP-SAT problem, not a whole-code wall-clock limit.
- Any parent call not returning complete `d_exact=true`/OPTIMAL is `UNSETTLED`; partial CP-SAT bounds are retained only as diagnostics and never affect classification.
- `BRAVYI_BB` numeric name/comment labels are `literature_reported_uncertified`. They never supply a certified parent lower bound; the four matching local JSON certificates take precedence.
- No Monte Carlo sampling is involved, so shots, failures, and Clopper–Pearson intervals are not applicable.
