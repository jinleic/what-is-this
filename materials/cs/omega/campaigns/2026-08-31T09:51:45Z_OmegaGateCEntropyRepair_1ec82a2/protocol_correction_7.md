# Protocol correction 7 — round-5: mass-aware generalized entropy (2026-08-31)

Status: withdrawn static acceptance honored. ZERO compute — only
`py_compile`, `--prepare-only`, and the extended mock self-test ran. No real
point evaluation, no entropy solve, no `results.json`, no decision.

## Owner-found defect, independently verified

The frozen `protocol_static_v3.json` records that the first moved
45-coordinate global distribution has EXACT dyadic mass
`37778931862949057407573 / 37778931862957161709568 = 0.9999999999997855`,
NOT 1. `vxxz24_float.entropy_vec` is generalized `-sum x log x` with no
renormalization. Therefore the pre-statement/simplex model was wrong for the
actual finite point: the old `exact_ok` rejected live p; the old projection
forced mass 1; raw pmax validation forced total 1; and `logsumexp - y*b`
bounded a normalized simplex problem that is infeasible for b of mass
m != 1. The KKT mass defect 2.1455e-13 < 2^-32 explains the false comfort.

## Repair (all ZERO-compute, mock-verified)

1. `mass_aware_feasible(p_intervals, matrix, pmax_intervals)` — exact mass
   m = sum(p) (Fraction, m > 0), positivity, pmax validated against
   total == m and rows == B p, all exact Fractions.
2. `exact_rational_primal(..., target_mass=None)` — normalization RHS is m
   (uniform fallback m/n, not 1/n); validation requires `sum(q) == m`.
3. `dual_float(..., mass)` — float objective/gradient consume
   `q(y) = m*softmax(B^T y)`, matching the generalized dual.
4. `interval_dual(..., mass)` — Arb-certificate uses
   `U_m(y) = m*logsumexp(B^T y) - y^T b - m*log(m)` (Fenchel-equivalent)
   with `q(y) = m*softmax(B^T y)`; KKT residual compares B q to b.
5. `m` is treated as an exact Fraction and lifted to Arb via
   `_fraction_point`. Primal H(p) is genuinely feasible without
   pretending m = 1.
6. Full-row/nonnegativity/mass exact checks are mandatory everywhere.

## m != 1 plants (added to the mock)

- mass-3/4 vector `[1/4, 1/4, 1/4]` FAILS the old simplex total==1 check
  and PASSES `mass_aware_feasible` (pmax `[1/4,1/4,1/4]` shares mass and
  marginals exactly).
- projection with `target_mass=3/4` succeeds with `sum(q) == 3/4 != 1`
  and BOTH full-row marginal equalities exact.
- `q(0) = m*softmax(0) = m/n` per coordinate sums to m = 3/4, NOT 1.
- the generalized dual upper at y=0 contains the brute exact feasible
  entropy (`m3q4_generalized_dual_contains_brute_entropy: true`).
- the observed base mass record is verified != 1:
  `37778931862949057407573/37778931862957161709568`.

## Interpretation correction (Main directive, both directives)

- The guard proves `sum(new) - sum(base) = 0` — mass is PRESERVED from the
  base to each rung — NOT `sum(new) = sum(base) = 0`; both masses equal the
  positive observed m ~= 0.9999999999997855.
- Therefore a strict result on these frozen points cannot be called a
  feasible witness or a counterexample to the constrained program's local
  optimality. The registered points/domain are preserved; the center is
  NOT normalized or moved post hoc.
- The only possible strict result on this campaign is:
  **"STRICT DESCENT of the repaired generalized-entropy/absorbed-defect
  objective on the named fixed-m affine ray at the released binary64
  center."** It explicitly does NOT produce a feasible construction, a
  local-optimality counterexample, or an exponent improvement. OPEN remains
  OPEN either way. This wording is now in the source `decision` dict.
- The earlier README sentence claiming a normalization/feasible ray is due
  for owner correction AFTER the run; no README edit has been made.
