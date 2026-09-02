# Protocol correction 6 — fourth re-audit response (2026-08-31, static only)

Status: ZERO compute. `py_compile`, `--prepare-only`, and the mock self-test
only. No real point evaluation, no entropy solve, no `results.json`.

## Owner-reproduced defects, both confirmed by direct probe

Probe (`.venv`, flint arb):

```
Fraction(arb) RAISES: TypeError argument should be a string or a Rational
                      instance or have the as_integer_ratio() method
Fraction.from_float(float(mid)) -> 3/8
tiny = 1/2^1200 : float(tiny) == 0.0 -> True ; tiny == 0 -> False
```

1. **Every pmax seed was silently None.** `Fraction(mid)` on a `flint.arb`
   raises TypeError, so the previous `try/except` converted EVERY seed to
   `None` and the projection fell back to the uniform `1/n_active` value —
   contrary to my report that pmax midpoints were used. FIXED: the seed is
   now `Fraction.from_float(float(mid))` after an `np.isfinite` assertion
   (the binary64 midpoint is itself an exact dyadic, and it is only a
   seed), plus a hard assertion that every seed is a `Fraction`.

2. **Float-based forced-zero detection was unsound.** `float(marginal.lo)
   == 0.0 and float(marginal.hi) == 0.0` underflows a tiny positive Arb
   interval to zero and would delete live support from the dual. FIXED:
   the rule is now the exact Arb test `marginal.lo == 0 and marginal.hi
   == 0`.

## New plants (both through the REAL code paths)

- **Arb -> seed conversion plant.** Seeds are built from actual Arb
  intervals (one with genuine width, midpoint 7/16), not prebuilt
  Fractions. Asserts: every seed is a Fraction; the seed vector equals
  `[1/16, 7/16, 3/8, 1/8]`; the vector is NONUNIFORM; no seed equals the
  uniform `1/4` (i.e. no silent uniform fallback); the projection
  succeeds; and the nonuniform Arb-derived values land in the FREE slots
  `{1,3}` as `q[1]=7/16`, `q[3]=1/8`, with exact normalization and both
  full rows exact.
- **Tiny positive Arb plant.** `tiny = 1/2^1200` has
  `float(lo)==float(hi)==0.0` yet `lo != 0` exactly. Asserts the old float
  rule WOULD force the column and the exact Arb rule does NOT.

## Mock self-test (verbatim, `mock_self_test_v4.json`)

```
arb_to_seed_conversion_plant_ok: true
arb_derived_seeds: [[1,16],[7,16],[3,8],[1,8]]
arb_seeds_reach_free_slots: [[7,16],[1,8]]
tiny_positive_arb_float_underflows: true
tiny_positive_arb_not_forced_by_exact_rule: true
honest_fallback_used: false        <-- still NOT exercised; reported as-is
```

The prior rounds' plants (exp endpoint bounds, wrapper invocation for both
Part and GlobalStage, overlap/negative rejection, noncontiguous-free and
forced-column projections) continue to pass unchanged.
