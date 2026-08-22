# Analysis of the extended simple-cubic series

## Exact input and scope

`experiments/e18_series_extend.py` derives the reduced free-energy series through
`v^20` at high temperature and `x^28` at low temperature.  Every coefficient is
an exact `Fraction`.  The high-temperature extension requires five free boxes
above the old 62-site limit, including `4x4x4` and `4x4x5` (80 sites); modular
propagation uses three 31-bit primes for each of those boxes.  A separate
`4x5x2` calculation exercised a 20-spin cross-section and agreed exactly with
the int64 engine.  Direct GF(2) cycle-space enumeration, with no spin transfer
matrix, agrees with the connected-cluster/free-energy coefficients through
`v^10`.

The singularity analysis uses

\[
F(z)=\phi(v)-\log 2,\qquad z=v^2,
\]

so the ten nonzero HT coefficients become ten consecutive powers of `z`.  Root
finding is performed at 80 decimal digits, but the uncertainty below is
set by series length, not numerical precision.

## Critical singularity

The combined, untuned estimate is

\[
v_c=0.2193 ± 0.0048,\qquad
K_c=\operatorname{atanh}v_c=0.2229 ± 0.0050.
\]

The central value is the median of the highest-complexity nondefective
near-diagonal Dlog-Pade poles and all balanced first-order inhomogeneous
differential-approximant poles at order `v^20`.  The quoted uncertainty is the
full stability envelope: the largest displacement from that central value
among every retained approximant at truncations `v^16`, `v^18`, and `v^20`.
This deliberately conservative rule was fixed without reference to the
published benchmark.

For comparison only, `0.221654626(5)` differs from our central `K_c` by
`0.00128`, which is
`0.25` of our quoted
uncertainty.  The benchmark was not used to select approximants, filters, or an
error bar.

### Exponent

Dlog-Pade was applied to `F''(z)`.  Its two highest-complexity nondefective
approximants give a median specific-heat exponent
`alpha = 0.51489236822484049102106684153289762`.  The differential approximants fit
`Q1(z) F'(z) + Q0(z) F(z) = P(z)` and give a median
`alpha = 0.32295249787210382439050649353962534`.  Their disagreement and strong
order drift are material: the combined value is only

\[
\alpha=0.35 ± 0.34,\qquad
2-\alpha=1.65 ± 0.34.
\]

The exponent uncertainty is the analogous envelope over all retained
approximants at `v^16`, `v^18`, and `v^20`.  This short series locates the
singularity much more stably than it determines the weak specific-heat
exponent; the exponent estimate should not be treated as precision physics.

## Finite D-finiteness search

The exact Euler-operator ansatz was

\[
\sum_{j=0}^r Q_j(z)\,(z\,d/dz)^j F(z)=0,
\qquad \deg Q_j\le d.
\]

All `20` pairs `(r,d)` satisfying
`(r+1)(d+1) <= 11 - 3` were tested.
The first `8` coefficients fit the
candidate; the last `3` were held out.  No candidate passed all held-out
coefficients.  The failed training relations are recorded in
`results/series/extended_series_analysis.json`, including their nonzero exact
residuals.

This is **not a proof of non-D-finiteness**.  It excludes only homogeneous
Euler-form ODEs inside the tested order-degree budget on the available prefix.
It says nothing about larger `(r,d)`, and a finite prefix can always support
accidental relations at sufficient complexity.  We would require at least
`10` newly derived coefficients,
not used in fitting, before calling any discovered ODE credible; this
conservative threshold is at least ten and at least the number of fitted ODE
parameters.

As a positive control, the identical code was applied to 19 exact coefficients
of `exp(z)`.  It found the order-one, degree-one equation
`z F - (z d/dz)F = 0` from the fitting prefix and obtained exactly zero on all
three held-out coefficients.  Thus the negative Ising result is not caused by a
search implementation that cannot recognize a simple holonomic function.

## Closed-form and integer-relation searches

No integer-relation or closed-form search was run.  The series-derived
uncertainty supplies far too few reliable digits for such a search to have
non-accidental predictive power.  Consequently there are no closed-form
candidates to report, and no exactness claim is made.
