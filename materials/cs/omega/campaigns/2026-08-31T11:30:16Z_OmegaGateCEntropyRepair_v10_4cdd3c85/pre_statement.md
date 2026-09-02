# PRE-STATEMENT — OmegaGateCEntropyRepair v10 producer

Written 2026-08-31T11:30:16Z before any v10 evaluation. PARENT:
`2026-08-31T09:51:45Z_OmegaGateCEntropyRepair_1ec82a2` (v9, owner-accepted
after frozen-artifact audit).

## Parent abort (honest, untouched)

The Main-released one-shot v9 run ABORTED at the BASE point in 0.98 s
(EXIT 1): `entropy_certificate` line 604 asserted
`len(p) == matrix_full.shape[1]` and the real Part's `split_dist` length
(num_split) did not match what the old `exact_feature_matrix` built. v9 is
preserved UNTOUCHED as an honest ABORT record; this campaign does not edit
it.

## Root cause and v10 interface fix (owner diagnosis + refinement)

`JointToMargin.mats[t]` has shape `(num_distribution_entries, margin_width)`
because `apply()` computes `dist @ m`. The certificate needs
rows = margin coordinates and cols = distribution entries. Two further
facts drive the final construction:

1. TRANSPOSE each block (`m.T`), stacking t=0,1,2. (v9 omitted this.)
2. IMPORTANT REFINEMENT (Main, pre-freeze): stored Part lambdas are
   COMPACT — `lam_margin[(r,t)]` has length `highs[t]-lows[t]+1`, while
   `j2m.mats[t]` has full width `sum_half+1`. Blindly stacking all m.T
   fixes columns but misaligns `branch_lambdas[selected_rows]` and
   includes structural-zero rows absent from the stored multipliers.

Therefore v10 builds B in EXACT lambda coordinate order:

    for each dim t in 0,1,2:
        Part (has lam_low/lam_high): m_t[:, lo:hi+1].T
        GlobalStage (no lam_low/lam_high): m_t.T (full width)
    vstack t=0,1,2

with these asserted guarantees:

- every omitted Part margin column is structurally all zero (hard assert);
- each retained slice width equals that lambda vector length;
- total B rows equals len(branch_lambdas(branch, region)) for each of the
  3 regions — verified structurally against ParamManager pm.size table
  (no live distribution access);
- B columns equal the registered distribution length
  (num_split for Part, n_shape for glob);
- entries exactly 0/1.

Verified against the REAL Workspace in prepare-only: all 63 real Parts and
the GlobalStage pass B_rows==lam width (e.g. Part (1,1,6) has B 7x4 with
lam widths 2+2+3=7; glob has B 27x45 with lam widths 9*3=27).

## Non-square + nonzero-lo mock plant (new, additive)

The non-square mock uses NONZERO `lo` offsets per dim (los = 1, 0, 3, with
widths 2, 3, 2), omitted columns structurally zero. Old no-slice stacking
gives rows = 3*9 = 27, which does NOT match total lambda width 7 — old
construction is structurally misaligned. New lambda-slice stacking gives
B shape (7, 5) — exactly the lambda-aligned form. This catches compact-row
misalignment that uniform-lo mocks hide.

## NO other changes

- NO change to dual/projection/decision logic (decision semantics carry
  the STRICT DESCENT wording from the parent).
- All v9 plants retained; the non-square + nonzero-lo plant is additive.
- All pre-run checksums from parent campaign remain valid history.

## ZERO real point compute on this campaign

- `py_compile` and `--prepare-only` smoke only (constructs the real
  Workspace, asserts the structural guard, records it, stops; NO entropy
  solve, NO evaluate_point, NO objective evaluation, no registered point
  evaluation).
- Mock self-test runs through `--mock-self-test` flag only (no real
  Workspace evaluation; mock objects only).
- NO heavy run and no point evaluation until Main releases the next CPU
  boundary. At release the frozen v10 artifact would execute the six-point
  ladder (base, R/{16,8,4,2,1}) exactly like the v9 release command with
  `results_v10.json` and `run_v10_stdout.log`.

## Rule 7 scope sentence (unchanged from parent, restated)

The swept set remains exactly base plus five fixed steps R/{16,8,4,2,1} on
one exact integer direction d (A*d=0, sum(d)=0, entries in {-1,0,1}),
changing only the 45 region-0 glob dist[0] coordinates of the frozen VXXZ24
`K100_2.37155181` vector under the transcribed single-p_comp max_level=3
q=5 program, preserving observed mass via sum(new)-sum(base)=0 with the
generalized-entropy/absorbed-defect objective at that mass m (never
normalized to 1). Nothing else searched: no other direction, step, radius,
box, region or parameter block, multi-block/cross-block class, q/max-level
regime, unpublished omega<2.371177 parameters, or other construction. Even
a strict result would establish only STRICT DESCENT of that objective on
the named fixed-m affine ray at the released binary64 center — not a
feasible construction, not a counterexample to the constrained program's
local optimality, and not an exponent improvement; OPEN remains OPEN. No
root ledgers touched (README feasible-ray sentence remains flagged for
owner post-run correction); no paper or record claim.
