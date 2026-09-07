# N6 precompute amendment — control and evidence-label clarifications

Run `20260904T052006Z_e5f2adce_0e7798f95678`; gate
`n6-qi-real13-descent`; preregistration commit `f8187e1`, SHA-256
`f6c7648a1a3a39ae84ded246b5dc45bbc4f7891890cc3821909d70e46b455564`.
Bound before any N6 scientific computation.

## 1. Krawczyk negative control

The preregistered exact-substitution rejection of a corrupted plant tests only
the tensor/witness serializer and is insufficient to test the certificate
instrument. Strengthen the fixed plant control as follows. At the accepted
plant seed, use the exact same dyadic center, QRCP chart-selection rule, and
complete radius ladder in the copied Krawczyk core on:

1. the exact planted tensor, which must produce strict containment with the
   outward-Arb check; and
2. a box-local negative formed by increasing only planted target entry
   `[0,0,0]` by exactly `1000`, which must produce no strict containment at
   any registered rung.

The Jacobian is target-independent, so the two calls must serialize identical
solved/free column lists and the same chart rule. Any containment for the
corrupted box, chart mismatch, or failure to contain the true plant is INVALID
INSTRUMENT and stops before TF. “No containment” is only rejection of these
fixed boxes around this fixed center; it is not global infeasibility and makes
no tensor-rank claim. The original same-witness exact substitution mismatch is
also retained.

## 2. Composite success label

A successful exact/Krawczyk real 13-term construction gives the upper bound
`rank_R(TF)<=13` labeled MACHINE-VERIFIED. Equality with 13 also consumes the
already frozen lower-13 chain with its existing cited-dependency component.
Therefore every equality restatement must use the composite label
**upper bound MACHINE-VERIFIED / lower bound FROZEN-CHAIN
(CITED-DEPENDENCY component retained)**. The whole equality must never be
called simply MACHINE-VERIFIED.

## 3. Continuation and certification charts

The continuation chart is the N6 QRCP proposal verified nonsingular by exact
rational realification. The copied Krawczyk core independently QRCP-selects its
certification chart at the admitted dyadic center; these charts need not be
the same. Candidate/results must serialize both chart kinds and both
solved/free column lists. For the certification chart, serialize every exact
dyadic center value at its 55 held-free coordinates. The independent replay
must verify (a) its reselected solved/free lists equal the serialized
certification lists and (b) every held-free coordinate equals the serialized
exact value before accepting containment. Exact-substitution mode records that
no separate certification chart was used and independently rechecks all 247
coordinates through the tensor equations.

No target coefficient, continuation schedule, scientific threshold, seed,
positive acceptance rule, or frontier inference changes.
