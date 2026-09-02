# Protocol correction 8 — round-6 response + release-integrity (2026-08-31)

Status: ZERO real compute (py_compile + prepare-only + mock self-test only).

## Fixes to the three round-6 defects

1. **m*log(m) algebra (was bare log m).** `interval_dual` now computes
   `upper = m_point * (ln Z + shift) - y^T b - m_point * ln(m)`, matching
   the Fenchel generalized dual `U_m(y)=m*logsumexp(B^T y)-y^T b-m*log(m)`.
   Replaced one-sided containment with a DIRECT two-sided identity plant:
   at y=0, |U_certified - (m*log(n) - m*log(m))|_upper = 9.03e-17 < 2^-40
   (tight), while the wrong (bare log m) formula is looser by
   m|log m| = 0.21576 (> 2^-8, recorded as
   wrong_formula_loosens_by_m_logm_value).

2. **Bad-multiplier plant call signature.** The base plant now calls
   `interval_dual(matrix, b_intervals, planted, target_mass,
   residual_matrix=residual_matrix, residual_b=residual_b)` with keywords —
   previously residual_matrix was passed as mass and residual_b as
   residual_matrix, which would error the production base plant before any
   decision.

3. **Fabricated mass eliminated.** A non-point/non-exact/non-positive p is a
   HARD abort (RuntimeError), and `mass_aware_feasible` recovering no
   positive mass also hard-aborts. The old silent `target_mass = ... else
   Fraction(1)` fallback is removed (`target_mass = exact_mass` after the
   positive-mass proof). The duplicate `stored dist_max REJECTED` predicate
   is replaced with a single check plus an unexpected-status abort.

## Round-trip hardening (Main's robustness directive)

Centralized `exact_dyadic_fraction(interval)` (module level): the value is
an exact dyadic point ONLY if a single binary64 round-trips to BOTH
original Arb endpoints exactly; any mismatch is a HARD abort (never a
silent b/m change). Plants: a 71-bit exact dyadic `(2^70+1)/2^70` —
old float conversion silently changes it to 1 (over53bit_changed_by_old_float
= [1,1]), converter rejects via round-trip mismatch
(roundtrip_converter_rejects_over53bit = true), while real 0.25 passes
(roundtrip_passes_real_dyadics = true). Production p conversion now goes
through this converter.

## Release-integrity (Main directive, provenance-only)

`save()` is atomic: write to a fixed same-directory temp (OUT + '.tmp'),
flush + os.fsync, then `os.replace(temp, OUT)`. A 900s SIGALRM can no
longer destroy the last complete JSON. No retries/reruns; no compute.

## Mock self-test result (verbatim, mock_self_test_v6.json)

```
dual_identity_y0_tight_enclosure: true
dual_identity_abs_separation_upper: 9.029685343895233e-17
wrong_formula_loosens_by_m_logm_value: 0.2157615543388357
over53bit_changed_by_old_float: [1, 1]
over53bit_true_value: [1180591620717411303425, 1180591620717411303424]
roundtrip_converter_rejects_over53bit: true
roundtrip_passes_real_dyadics: true
honest_fallback_used: false   <-- fallback still NOT exercised; as-is
```
