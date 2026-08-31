# Pre-statement — new-decomposition off-diagonal census (`mws59`, `stapleton60`)

**Created 2026-08-31T083323Z by the producer and committed BEFORE any target computation.**
Campaign: `/Users/jinleic/jinleic-workspace/cs/mm3/campaigns/2026-08-31T083323Z_55447c11-5887-4fdd-aa2f-30e31f2332e8_b3046e6a6c50`.  Producer SHA-256 prefix: `b3046e6a6c50`.

**Rule-5 correction / prior abort:** producer hash `4dc928ac6dbe` at campaign
`2026-08-31T082412Z_a2d87381-157e-45bb-ad3a-519ddf899bc1_4dc928ac6dbe`
aborted in 0.26 s before any anchor or target count: it attempted to iterate the
intentional `None` in `META["perminov58"]["split"]`.  That failed source and log are
frozen, with no target claim.  This code path instead decodes the pinned Perminov JSON,
asserts zero factor mismatches and `complexity.reduced == 58`, and receives this new hash
and separate pre-registration before rerun.

## 1. Fixed question and falsifiable outcomes

The completed off-diagonal programme covers only `paper55` and `sun56`; its undecided
remainder is empty.  This campaign therefore takes the highest-value listed extension:
a **different decomposition**.  It fixes exactly the two remaining on-disk, previously
729/729-verified public rank-23 decompositions: `mws59` (published 59 additions) and
`stapleton60` (published 60 additions).  `perminov58` is not a new decomposition here:
its factor blocks are the already-covered `paper55` blocks.

Exactly one verdict is admissible per fully decided data triple:

* **NOGO:** exact lower bound at least 55, via the exact d-counter, complete subset-DFS
  floor census, and the +14 transposition gap with its rank/activity precondition checked.
* **LIVE:** lower bound at most 54; halt before recording the row and escalate to Main.
  A low lower bound is not a scheme.  It then enters a separately frozen exact scheduling
  adjudication; any SAT witness is replayed over Z, and any reported UNSAT must have DRAT
  converted to LRAT and accepted by both `drat-trim` and `lrat-check`.
* **RECORD:** an explicit <=54 circuit, only after exact gate-list replay and all 729 Brent
  identities over Z; escalate to Main before writing anywhere.

The primary output is the exact survivor count and exact certified-lower-bound histogram
for each decomposition.  If a wall cap intervenes, only the lexicographic prefix is
reported and the exact remainder named; no total claim is made for it.

## 2. Domain fixed before results (rules 14-16)

Let T be the **set** (not a group) of all 6960 ternary 3x3 integer matrices with
determinant +/-1, exactly re-enumerated from all 3^9=19,683 ternary matrices.  For each
D in {`mws59`, `stapleton60`}, enumerate every independent triple
(P,Q,R) in T^3 and then every post-sandwich sigma power k in {0,1,2}, where
sigma(U,V,W)=(V,W^T,U^T).  Full raw domain: 2,022,921,216,000 scheme-instances
(1,011,460,608,000 per decomposition).  The 6960 set's failure to close under
multiplication is irrelevant: each P,Q,R is used independently and no product is assumed
to remain in T.

Each factor row is a row-major 3x3 block.  In **standard matrix semantics** the fixed
family is

    U' = P^-1 U Q^-T,   V' = Q^T V R^-T,   W' = P^T W R.

This is re-derived from F(A,B,C)=tr((AB)^T C), not tr(ABC).  Convention warning: frozen
`gatec_sweep.sandwich` uses transposed-L letters; on monomials frozen (X,Y,Z) equals the
standard family at (P,Q,R)=(X^T,Y^T,Z^T), checked bit-for-bit.

An instance is admitted only when all 23 rows of all three transformed blocks remain in
the alphabet {-1,0,1}^9.  Brent validity follows from the proved automorphism and is
also spot-checked 729/729 in two exact paths.  P^-1,Q^-1,R^-1 need not themselves be
ternary.  Non-ternary factor images are counted as excluded from the fixed alphabet, not
silently dropped.

## 3. Exact factorization and completeness

T is partitioned by free right multiplication by the 48 signed monomials M:
T=A*M, |A|=145, every orbit size 48, representatives chosen lexicographically.  If
P=a1 p1, Q=a2 p2, R=a3 p3, the outer p_i induce a common signed permutation of the
nine input coordinates on each side.  Ternarity, d, floor existence, activity, and
circuit cost are invariant.  Thus every data triple (a1,a2,a3) represents exactly
48^3 sandwich triples, and sigma contributes a further factor 3.  The campaign decides
exactly 6,097,250 (decomposition, data-triple) points across the two decompositions
before the proven fiber multiplier.
The factorization is accepted only after: all 6960 elements partition into 145 size-48
orbits with unique factorization; a full anchored 145^2 direct-vs-table subcube for each
decomposition has zero mismatches; 200 deterministic random direct/table triples per
decomposition agree; 24 numpy-vs-pure-Python entry-level checks agree; and 25 random
right-monomial fiber predicates per decomposition agree.  Any mismatch aborts before a
count is reported.

Counting uses exact int64 only after the hard entry bounds |U'|<=36, |V'|<=18,
|W'|<=9 are proved from |G|<=1, |G^-1|<=2; counts are accumulated in Python integers
and independently as a bounded int64 matrix contraction (<145^3<2^63).  No float,
interval arithmetic, HiGHS, SAT, or tolerance enters the main verdict.

## 4. Bound semantics and direction

For one side F, d(F) is the number of distinct non-input sign classes among its 23
needed forms.  Each addition creates at most one new class, so C(F)>=d(F).  At exactly d
gates every gate must be a needed class; complete memoized subset search enumerates every
ordering and every representation t=+/-a+/-b (including a=b).  Failure is therefore an
exact finite census proving C(F)>=d+1, not a solver UNSAT label.  A found floor schedule
is replayed gate by gate.

For the output map Wfac^T:23->9, reverse-mode transposition of an active A-gate circuit
uses A+9-23=A-14 additions for Wfac:9->23, hence C(output)>=C(Wfac)+14.  The permissive
direction is explicitly excluded.  Startup asserts every slot has 23 nonzero rows and
exact Q-rank 9; the honest action is invertible on the nine coordinates, so the
precondition holds on the full domain.  Total lower bound is

    U_lb + V_lb + Wfac_lb + 14.

Sigma only cycles the three side costs and transposes blocks (a coordinate permutation),
so the total is invariant; it multiplies the orientation count by exactly 3.  Twelve
survivor totals per decomposition are replayed through all three sigma powers and through
independent right-monomial fibers.

## 5. Mandatory anchors and counterfactual controls — before new counts

1. Reproduce 55/58/56/59/60 published totals and 729/729 Brent over Z in independent
   Python-int, fmpz, and frozen paths.
2. Assert the **all-monomial data triple is exactly 55** on `paper55` and `sun56`
   before any new histogram.  Also assert the inherited all-monomial certified totals
   56 for `mws59` and 58 for `stapleton60`, across their sigma classes.
3. Reproduce the frozen d/floor tuples for paper55, sun56, mws59, and stapleton60 with
   both this independent DFS and the frozen DFS.
4. Fixed nonmonomial G=[[-1,-1,-1],[0,1,1],[0,0,1]]: the final family must be ternary
   and Brent 729/729 on sun56; the wrong W'=G W G^T and transpose-as-inverse plants
   must fail.  A deliberate nonzero W-coefficient flip must produce Brent failures;
   det-2 must be rejected.  These are the controls that give a zero-event result meaning.

## 6. Caps, order, escalation, and non-goals

Lexicographic order in (decomposition,p,q,r), with mws59 before stapleton60.  Fixed seed
20260831 only for controls.  Caps: anchors/controls 2 h, tables 8 h, decisions 24 h,
total 48 h; unfinished work is an exact named remainder.  Any <=54 row halts before the
row is written and is escalated to Main; no computation is adjusted toward a published
value.

Not searched: `paper55`/`sun56` beyond startup controls (their full domain is already
closed); decompositions outside the five named public ones (Laderman, other Smirnov,
Schwartz-Vaknin, unpublished, F2-only); non-ternary factor alphabets; GL(3,Q) or GL(3,Z)
beyond T; actions outside the proved family; anti-cyclic pair swaps (compute B*A); and
upper-bound circuit synthesis unless a <=54 adjudication triggers.  No radii exist: all
objects are exact integers.
