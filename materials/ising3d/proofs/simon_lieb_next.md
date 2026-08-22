# Simon--Lieb next step: an exact finite-family obstruction, no improvement

Artifacts: `experiments/e66_simon_lieb_next.py`, `results/bounds/simon_lieb_next.json`,
`tests/test_simon_lieb_next.py`.

Provenance note: the producing agent was killed by a provider quota after writing the
experiment and artifact; the standalone verifier and this note were completed by the lead
session from the on-disk state, and every certificate below was independently recomputed
before acceptance. An earlier in-flight run that appeared to show an improvement was an
error in that run's certificate points and exited nonzero; the final exact result is the
negative one recorded here.

## 1. Question and criterion

Can the already-audited Simon--Lieb finite-set criterion, evaluated exactly on larger boxes,
push the certified lower bound `K_c >= atanh(v*) = 0.21221190116616783933...` (the exact
`c_36` submultiplicativity endpoint) any higher?

For a finite box `S` with free boundary and the high-temperature parameter `t = tanh K`, the
criterion coefficient is the rational function

    phi_S(t) = t * Q_S(t) / P_S(t),

and `phi_S(t) < 1` implies exponential decay of correlations at `t`, hence `atanh(t) < K_c`.
The sign of `phi_S(t) - 1` at rational `t = p/q` is the sign of the exact integer
`p*Q - q*P` evaluated at `p/q` in scaled form (`tQ - P` below), so every certificate is one
integer sign — no floating point.

## 2. The exact family and the challenge point

**[COMPUTATION]** Twenty open boxes were evaluated exactly: all sorted triples
`2 <= a <= b <= 4, b <= c <= 5` with transfer cross-section `a*b <= 16`, plus the prisms
`4x4x8`, `4x4x12`, `4x4x16`, `4x4x24` — a fixed size-then-lexicographic order chosen before
any evaluation, with no benchmark input. For each shape, a modest-denominator safe point
`t = p/10000` has certified `tQ - P < 0` (exact integer sign), giving the per-shape
lower-bound endpoints recorded in the artifact. The best completed endpoint is

    atanh(2023/10000) >= 0.20512954053188980029...   (shapes 4x4x16 and 4x4x24)

which is BELOW the incumbent `0.21221190116616783933...`.

**[COMPUTATION — the obstruction]** The fixed challenge point `t = 23/110` satisfies
`atanh(23/110) >= 0.21222050478358502326... >` incumbent, so ANY shape certifying it would
strictly improve the interval. Every one of the twenty shapes REJECTS it: the exact integer
`tQ - P` is nonnegative (`phi_S(23/110) >= 1`) on all of them. Both sign families carry
sha256 fingerprints of the exact integers, and the standalone test recomputes four of them
from the library alone (never importing the experiment).

## 3. Controls

* Complete spin enumeration on `2x2x2` produces the full integer polynomials `P, Q`; the
  independent parity-transfer route reproduces both at `t = 1/5` exactly.
* A second control shape `2x3x2` at `t = 2/9` checks the direct residual stream against the
  full `P, Q` values exactly.
* All `atanh` endpoints use mpmath interval arithmetic at 90 dps with outward rounding and
  40-place floor reporting; the test reproduces both decisive decimals independently.

## 4. Verdict and scope

**[UNRESOLVED — honest negative]** The scalar Simon--Lieb criterion on this completed exact
family cannot reach the incumbent: the saturating cross-section `<= 16` family tops out near
`0.2051` and rejects the first useful challenge point outright. This is an exact obstruction
for THIS route, not a global no-go: larger cross-sections (`2^{a*b}` parity-transfer states),
diamonds/cylinders, vector-valued or multi-point refinements of the criterion, and any
sharper exact contraction all remain open. The certified interval endpoint is unchanged:

    K_c >= 0.2122119011661678393310862783954278184914   (c_36 route, unchanged)

This result is consistent with, and complementary to, the wave-8 infrared barrier theorems:
both endpoints of the current interval now carry explicit method-class obstructions.
