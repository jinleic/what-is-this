# Gate C diagnostic/repair report

## Headline and verdict

**V1 is repaired and machine-equal to the frozen endpoint, but the honest
branch-safe full-box enclosure is 73.85977508485647 times the inherited
signal. Verdict: FAILURE TO CERTIFY; the 21-dimensional question remains
OPEN.** A numerical candidate direction exists, but this campaign establishes
neither an improving witness nor local optimality/no improvement.

Campaign: `2026-08-31T07:57:00Z_OmegaGateCDiag_4859327`.

Freeze/source commits:

- `4859327`: original pre-statement at `omega/pre_statement.md`, before every
  computation in this campaign.
- `256e22f`: campaign-directory byte-copy of that pre-statement, before the
  final authoritative replay (but after the first DIAGNOSTIC ONLY pass; this
  timing deviation is disclosed rather than backdated).
- `5b1e327`: replay, quotient-rule, exact-basis, LP, and endpoint repairs.
- `e9c6414`: branch-safe full-box source executed for the final replay.

## The requested first diagnostic

The first attempted point replay ran for 122.5 s and then crashed in the
new diagnostic harness because a scalar `Slope` was indexed as a list; it
produced no verdict. The first successful pass was explicitly **DIAGNOSTIC
ONLY** and overbroad: it zeroed every Part and glob hash-penalty RHS, not just
the `lam_sum`/Lemma-1 residual. That scope error is retracted in
`diagnostic_protocol_correction.md`; the result is retained as a
counterfactual, never as a certified path.

- Real replay raw endpoint: `2.371553835835081`
  (**COMPUTATIONAL-EVIDENCE**).
- All hash penalties zeroed: `2.3715518061842102`
  (**COMPUTATIONAL-EVIDENCE**).
- Shift: `-2.029650870927213e-6` (**COMPUTATIONAL-EVIDENCE**).

The authoritative final run performs the prescribed narrower test: zero only
the `2*eps` Lemma-1 contribution while retaining entropy mismatch `hm` and
`p_comp`. All seven resulting R branches, R, M, and raw endpoint are
machine-equal to the frozen `include_lemma1=false` record:

- real raw minus Lemma-1-zero raw = `2.029648699330977e-6`
  (**MACHINE-VERIFIED**);
- Lemma-1-zero raw = `2.3715518061863814`
  (**MACHINE-VERIFIED**).

Thus there is no approximately `3e5` downstream amplification. The earlier
`+0.9008715` is not a penalty effect.

## V1 and the stale `+0.9008715` provenance

The repaired endpoint replay machine-matches all seven frozen branches:

| branch | value (**MACHINE-VERIFIED**) |
|---|---:|
| R_comp[0,3] | 0.24530661807075219 |
| R_comp[1,3] | 0.24715974210930786 |
| R_comp[2,3] | 0.2456557897939231 |
| R_comp[0,2] | 0.5883309785276519 |
| R_glob[0] | 0.49684932826618977 |
| R_glob[1] | 0.4968521822710635 |
| R_glob[2] | 0.49684892842208755 |

`R = 2.8170035674609757`, `M = 2.0942543887102634`, and raw endpoint
`2.3715538358350807` are **MACHINE-VERIFIED** against the frozen record. The
V1 residual against the registered raw anchor `2.3715538358350803` is
`4.440892098500626e-16`, so V1 passes.

The prior manifest's raw `3.272425321778391` and residual
`+0.9008714859433109` are **REPORTED**, not reproduced. The advertised stale
source has SHA-256
`4504365a7917f80362424a7252203cd770b7c51392061e517f8c10e62fc946fa` and,
when executed, raises `NameError: name 's_const' is not defined` at line 191
before its first aggregate checkpoint (**COMPUTATIONAL-EVIDENCE**). Therefore
that number has no reproducible producer among the advertised frozen bytes;
it is not retrofitted to repaired code. If one algebraically combines that
reported raw value with the correct M, the implied stale R
`0.9303495043602614` and R deficit `1.8866540631007143` are only
**INFERENCE**.

Classification of the requested alternatives: **(b), a fourth independent
provenance/aggregation defect**, not amplified epsilon and not (c). More
precisely, the advertised failing path is non-executable and the repaired
path already reproduces the frozen endpoint before any penalty toggle.

## Repairs and planted old-path failures

1. The endpoint slope now differentiates `Omega=(T-R)/M` directly with the
   interval quotient rule; the undefined `s_const` and incorrect M-chain term
   are gone.
2. V1 is a hard gate: no LP executes on a V1 failure.
3. The actual kernel basis comes from the exact rational JSON (all
   denominators 1, entries `-1/0/1`, exact integer `A V = 0`). It is not equal
   to the old dense SVD array; the latter has max float residual
   `4.127839899618667e-15` (**COMPUTATIONAL-EVIDENCE**).
4. LP units use `|V u|_inf <= 1` and `delta=r V u`. The old extra `/r`
   counterfactual changes the effective radius from `1e-7` to
   `9.999999999999998e-15` and the midpoint objective from
   `-1.577332901516852e-7` to `-1.577332255043032e-14`
   (**COMPUTATIONAL-EVIDENCE**).
5. The correct dual has 90 inequality-row marginals (45 upper plus 45 lower).
   The old 21 variable-bound marginals cannot multiply `V^T` of shape 21x45;
   the counterfactual records the shape failure.
6. Slope leaves now form endpoints in Arb as exact
   `IC.const(center) +/- IC.const(radius)`. The old binary64 subtraction made
   24 lower endpoints and 25 upper endpoints inward
   (**MACHINE-VERIFIED**).
7. Every R/M minimum hulls every interval-overlapping possible minimizer.
   Endpoint-argmin branch pinning is not consumed: several branches overlap
   even in the nominal point pass, so pinning one would be unsound in the
   permissive direction.
8. The ANY-y control is executed on the original Arb gradient intervals over
   the full named box. `y=0` gives swing upper
   `1.3787196625196844e-5`; the fitted arbitrary y strictly tightens it to
   `5.841014555619999e-6`, so the instrument is demonstrably capable of a
   nontrivial response (**MACHINE-VERIFIED**).

The inherited V2/V3 float arrays disagree with the repaired branch-safe
expression (max coordinate-reference deviation `1.2898328055870454`, max
kernel-reference deviation `0.44761912584538077`; both
**COMPUTATIONAL-EVIDENCE**). They are retained only as provenance diagnostics
and are not certificate gates.

## LP candidate versus rigorous interval result

The exact-basis midpoint HiGHS LP returns `-1.577332901516852e-7`
(**COMPUTATIONAL-EVIDENCE**). This is a numerical candidate direction only.
HiGHS feasibility/KKT residuals print as zero in binary64, but no solver
verdict is consumed by the certificate.

For arbitrary fixed dyadic y and every gradient g in the branch-safe Arb
bundle and every delta in D,

`g·delta = (g-A^T y)·delta >= -r ||g-A^T y||_1`.

Direct Arb evaluation gives the enclosure
`[-5.841014555619999e-6, +5.841014555619999e-6]`, width
`1.168202911124e-5` (**MACHINE-VERIFIED**). The predecessor signal
`1.5816497000997742e-7` is a **CITED-DEPENDENCY**; the width/signal ratio
`73.85977508485647` is **INFERENCE**. Because zero lies in the enclosure, the
rigorous verdict is **FAILURE TO CERTIFY / OPEN**. It is not “no improvement,”
not local optimality, and not an improving witness.

## Mandatory rule-7 scope sentence

The swept set was exactly `D = {delta in R^45 : A delta = 0 exactly,
|delta_i| <= r for every i=0,...,44}`, where A is the frozen exact 27x45 0/1
margin matrix (rank 24 over Q), `r` is the binary64 number
`0x1.ad7f29abcaf48p-24` (display `1e-7`), and the center is all 45 binary64
coordinates of the region-0 glob `dist[0]` block of the frozen VXXZ24
`K100_2.37155181` MATLAB parameter vector; the run used one point pass at
radius 0 and one un-subdivided full box whose exact Arb endpoints were
`IC.const(center_i) +/- IC.const(r)`, with every other parameter held at its
frozen center, under the transcribed 3-region single-`p_comp` program at
`max_level=3`, `q=5`. It did **not** search region-1/2 glob blocks, any
part/split/lambda/region_prop or other parameter block, multi-block or
cross-block directions, any other radius, subdivision, Monte Carlo, any
point outside D, any other q/max-level regime, the unpublished parameters for
`omega < 2.371177`, or any construction differing from the enclosed
laser-method program (including asymptotic-rank/centroid-based improvements).
Nothing here bears on the published or record exponents, and no record claim
is made.

## Premises disproved and next campaign

Disproved premises:

- The frozen source did not reproducibly produce `+0.9008715`; it crashes
  before that checkpoint.
- The penalty/epsilon path cannot carry that residual: the exact narrower
  toggle moves raw Omega by only `2.029648699330977e-6`.
- Single-argmin branch pinning and binary64-rounded box endpoints were not
  certificate-safe; both could be permissive and were replaced before the
  final result.
- The anticipated approximately 1.81x width depended on unsafe branch
  pinning. The honest hull is 73.85977508485647x the imported signal.

Named next campaign: **`OmegaGateCCandidateWitness`**. Lift the midpoint LP
solution through the exact integer kernel basis into an exact dyadic feasible
delta (scale inward if an exact coordinate reaches beyond r), then directly
Arb-certify the moved endpoint against the frozen center. This asks the
cheapest decisive question—whether the numerical candidate itself is a real
improving witness—without claiming anything about all of D. Estimated cost:
two certified point evaluations plus exact rational feasibility checks,
approximately **2–5 CPU-minutes**, with no full-box sweep.
