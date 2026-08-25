# Wave 26 research plan

**Status:** frozen before implementation, with the normalization correction below recorded before any artifact was written.
**Date:** 2026-08-24.

## Goal

Advance the ranked open `2x4` spectral front from a length-96 coefficient quotient to an exact,
nontrivial trace-nine projection theorem. The breakthrough target remains the complete primitive
`F9` norm, exact isolation of every `q>1` root, and shifted-Stieltjes disposition of every branch.
The required bounded fallback is stronger than the existing graph-uniform finiteness theorem:
certify at one physical rational `q0>1` that multiplication by `F9` is invertible in the
96-dimensional `F4,...,F8` quotient, then derive an explicit degree bound for one nonzero
univariate projection polynomial.

A result lands only with an exact producer, a clean-room verifier, a scoped proof, a machine
artifact, current source hashes, an adversarial review, ledger/checkpoint integration, and a
validated evidence bundle. No finite-layer statement is promoted to a thermodynamic solution or a
critical-point bound.

## Evidence before the freeze

- e248 exactly supplies `H1,...,H9`, the normalized eight-mode Lucas incidence, equation degrees
  `(2,2,2,3,4,3)`, and the parameter-independent 35-relation leading ideal. Its 96 standard
  monomials have maximum degree eight. It does not materialize the `F9` norm.
- The complete-intersection Hilbert series
  `(1+t)^3(1+t+t^2)(1+t+t^2+t^3)` predicts the degree vector
  `(1,5,12,19,22,19,12,5,1)`. e251 must recompute this from the literal standard monomials.
- The strongest local implementation pattern is e238's fixed-leading-form quotient reduction.
  Generic `QQ(q)` F5B is closed by e248's 900-second non-result. e95 supplies the applicable
  modular determinant and CRT discipline, but its quadratic interpolation does not scale to a
  six-figure degree envelope.
- Current primary CAS sources were checked for modular Groebner, multiplication-matrix, and
  polynomial-matrix determinant support. Singular, Sage, msolve, and python-flint are not installed
  in this project environment. No dependency will be added for this wave.
- **Pre-implementation correction.** The first exploratory modular preflight used `q^(+10k)`
  instead of e248's required `q^(-10k)` normalization. Its determinant residue was discarded
  before producer implementation and supplies no evidence. Producer and verifier must independently
  rebuild every witness with the correct negative exponent.
- Pre-freeze host load was `5.65 / 4.22 / 3.67` on the shared 28-core, 96-GB workstation. Every
  heavy stage remains single-process, below 900 process CPU seconds and 2 GiB task RSS.

## Theorem and evidence contract

### A. Exact modular nonvanishing witness

1. Consume the landed e248 trace artifact after checking every inherited semantic and provenance
   hash. Do not recompute the 1.24-GiB trace-nine recurrence inside the elimination process.
2. Reconstruct the six abstract Lucas-resultant rows independently of their stored hashes.
3. At deterministic physical rational points, evaluate every trace target modulo safe primes,
   reject every zero denominator, and build the specialized `F4,...,F8` Groebner basis.
4. Require the exact generic leading-monomial signature, 35 relations, the 96 literal standard
   monomials, and Hilbert vector `(1,5,12,19,22,19,12,5,1)`.
5. Reduce all 96 products `F9*b`, form the `96x96` multiplication matrix, and take its determinant
   exactly in the finite field. At least one nonzero determinant residue is required.
6. Include a synthetic positive eight-mode control whose known mode polynomial annihilates all six
   rows and forces determinant zero.

A good-prime nonzero residue proves the characteristic-zero quotient norm is nonzero at that exact
`q0`; it is not a numerical inference.

### B. Explicit projection-degree theorem

Let `D` be the product of the three denominators in `Q(0),Q(2),Q(3)`. From the e248 target records,
prove `deg D=64` and that every coefficient of `D Q(u)` has `q`-degree at most 74. Clear the six
resultant rows and verify the safe coefficient-degree bounds

`(192,200,224,296,384,392)`.

For equation degrees `(2,2,2,3,4,3)`, verify the projective-resultant multidegrees

`(144,144,144,96,72,96)`

and their dot product `182400`. The fixed leading system has no projective zero, so the modular
nonvanishing witness proves that the cleared resultant is a nonzero univariate polynomial. The
permitted theorem is: the trace-four-through-nine incidence, and therefore every positive full
eight-mode spectrum of the finite open `2x4` transfer matrix, can occur at no more than 182,400
complex `q` values at which the normalized targets are defined. The physical `q>1` set is a finite
subset. This does not locate a root or prove emptiness.

### C. Close cheap Lucas-product shortcuts

Derive the shifted factors at `x=u-4>=0` and certify the seven necessary target inequalities

- `r4 > r6/r2`,
- `r7/r1 > r5/r1`,
- `r8 > r9/r1`,
- `r4^2 > r8`,
- `r4*(r6/r2) > (r5/r1)^2`,
- `(r5/r1)*(r7/r1) >= (r6/r2)^2`,
- `(r6/r2)*r8 > (r7/r1)^2`.

For the literal open-`2x4` targets, clear positive denominators and require every numerator shifted
by `q=1+x` to have coefficients of the claimed sign. These checks are exact filters only. Their
success explains why the quotient projection remains necessary; it is not evidence for a mode
representation.

## Aggressive route and stop rules

- After the witness, pilot a fixed-leading-form or modular interpolation route only if a proved
  denominator, degree, and coefficient-height envelope projects below the per-stage wall.
- Never infer degree from apparent stabilization, reconstruct integer coefficients without a
  `modulus > 2*height_bound` certificate, or discard chart-denominator roots.
- The coarse bound needs 182,401 interpolation nodes. A generic specialized Groebner node measured
  in preflight is already tens of process-CPU seconds; therefore that direct interpolation is
  forbidden. A compiled lifted reducer may be pursued only with an exact template and a fresh
  projected-cost gate.
- Stop before 840 process CPU seconds or 1.75 GiB observed task RSS to preserve write/cleanup
  margin. Record the first blocked gate rather than silently narrowing a full-norm claim.
- No emptiness claim without the primitive norm, complete exact `q>1` root isolation, and exact
  disposition of every branch. Repeated mode roots require shifted `H0/H1` moment-localizer checks;
  Hermite distinct-root signature alone is insufficient.

## Verification contract

- Producer and verifier share no implementation code. One uses heap-based finite-field reduction
  and vectorized modular elimination; the other uses a separately written reducer and scalar
  Gaussian elimination.
- The verifier rebuilds the abstract incidence, every finite-field matrix, determinant residue,
  Hilbert count, degree-bound arithmetic, Lucas inequalities, theorem text, scope, and source hashes.
- Universal tables fail on any uncovered row. Resource measurements use process CPU and Darwin
  `mach_task_basic_info.resident_size_max`.
- Final review attacks reduction order, good-prime logic, projective/resultant degree accounting,
  denominator positivity, synthetic-control validity, and the boundary between finiteness and
  emptiness.
