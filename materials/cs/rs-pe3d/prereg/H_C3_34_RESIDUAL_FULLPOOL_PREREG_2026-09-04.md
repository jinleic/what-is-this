# PREREG-RSPE3D-21-C3-34-RESIDUAL-FULLPOOL — gate **H-C3-34-RESIDUAL-FULLPOOL**

Written 2026-09-04 UTC by agent `RsSparkClose`, immediately after
H-C3-34-RESIDUAL-POINT (run `20260904T051205Z_d03471d8_0604b708af38`) was
frozen and closed **FROZEN-INCONCLUSIVE**. That run proved the universal
Bezout residual-length lemma (Lemma R1) as a separate subresult and completed
a registered one-million-pair GF(17) sample with outcome E-empty. Because the
sample covered only 1,000,000 of the 41,587,200 compatible unordered pairs of
the pinned pool, its emptiness was inconclusive for the family. This gate
removes exactly that incompleteness: it exhaustively classifies **all
41,587,200 compatible unordered equation pairs** of the same pinned pool. It
is still a finite GF(17) statement about one registered factor family; no
universal, extension-field, or characteristic-zero claim may be inferred from
it, with or without a positive witness.

Binding order: this preregistration -> path-scoped commit -> `campaign.py
init` -> byte-identical in-run prereg copy with source commit and SHA-256 ->
standalone instrument (vendored from the frozen POINT driver; no import of
frozen code) -> controls before the sweep -> fixed-chunk checkpoints ->
classification record -> freeze -> exactly one close verdict. No frozen run
and no root `cs/RESULTS.md` is edited.

## 1. Registered population (no sampling)

The pools are reconstructed bit-identically from the frozen POINT
preregistration: pool seed `0xC33417`, 96 canonical PGL $(1,1)$ factors, 96
canonical coprime degree-two $(1,2)$ factors over GF(17), combined pool SHA-256

`0b6d92cf60e051e5fad09f4101da320b8c75430b976bab095734fc3b0be5b120`.

Equations are the 9,216 products $L_MQ_R$. The sweep enumerates every
unordered pair of distinct equation indices exactly once, in plain
lexicographic order $(first, second)$ with $0\le first<second<9216$; there
are

$$\binom{9216}{2}=42,462,720$$

such pairs. A pair is **compatible** iff its two equations share neither
factor: $l_1\ne l_2$ and $q_1\ne q_2$. Because equation indices biject with
$(l,q)$, two distinct indices can share at most one factor class, so the
incompatible pairs are exactly the $2\cdot 96\binom{96}{2}=875{,}520$
shared-factor pairs and the compatible population is exactly

$$42,462,720-875,520=41,587,200=2\binom{96}{2}^2.$$

The instrument asserts this count both by exhaustive enumeration and against
the closed form. No seeds, no permutation, no cap, no adaptive extension:
every compatible pair is checked exactly once.

## 2. Per-pair classification (identical semantics to the frozen POINT run)

For each compatible pair, the affine common-zero set on the full $17\times17$
grid is computed by two exact routes — direct evaluation of $F=L_MQ_R$,
$G=L_NQ_T$, and the union of the four factor-pair intersections — which must
agree. Per-pair record: intersection size $s\in\{0,\dots,12\}$ (a
size histogram over all 41,587,200 pairs is a registered output). A pair is
**qualifying** iff $s=12$ and the 12 points are pairwise distinct with 12
distinct $x$- and 12 distinct $y$-values (coprimality and Bezout length
$L\cdot Q+Q\cdot L=12$ then make these 12 points the reduced projective
intersection, so nothing is unseen).

Every qualifying base set $Z_0$ is classified exactly as in the frozen POINT
preregistration, with rank computed before any deletion work:

1. base rank of the 12 evaluation columns in the
   $H^0(\mathcal O(2,3))\cong F^{12}$ tensor-product Vandermonde basis
   (dual independent implementations must agree); recorded for
   base-rank-10 vs below-10;
2. for each of the 12 residual deletions $S_q=Z_0\setminus\{q\}$: rank
   (must be exactly 10, else recorded as rank-deficient with witness), then
   relation-space dimension (must be one) and full support (all eleven
   relation coefficients nonzero); if rank-10 but the relation is not
   full-support, the deletion is recorded **nonminimal** with its full
   witness (ranks of all eleven 10-point sub-deletions included); if
   full-support, the deletion is recorded **circuit**.

Thus each qualifying base receives 12 rank-before-deletion checks and up to
132 ten-deletion checks, with no early exit that could suppress a later
failure. First witnesses of each negative class are pinned. The four
size-12 non-injective pairs already found in the POINT audit
(`[1663,6459]`, `[4534,5905]`, `[4791,5611]`, `[5861,6545]`) are registered
sanity anchors: the sweep must rediscover exactly these four size-12
intersections unless a hash-collision-defying difference appears, in which
case the discrepancy is itself a recorded outcome.

## 3. Controls before the sweep (same battery, seeds removed)

1. Dual rank implementations agree on the 12-column tensor-product
   Vandermonde basis (rank 12) and on all later promoted supports.
2. Circuit ACCEPT: ten independent ambient basis vectors plus their
   full-support sum — rank 10, nullity one, full-support relation, all
   eleven deletions rank 10.
3. Full-rank promotion REJECT: eleven independent vectors, rank 11, all
   deletions rank 10 — must be rejected by the rank-before-deletion gate.
4. Rank-10 nonminimal REJECT: $e_0..e_9$ plus duplicate $e_0$ — rank 10,
   nullity one, deletion-rank vector $(10,9,9,9,9,9,9,9,9,9,10)$, non-full
   support — must be rejected as nonminimal.
5. Geometric ACCEPT: the frozen GF(13) $(3,3)$ complete-intersection anchor,
   standalone (no frozen import), rank 7 with every deletion rank 7.
6. Factor accepts/rejects: valid PGL and degree-two accepted; singular PGL
   and common-factor degree-two rejected.
7. Shared-factor rejects: planted shared-$L$ and shared-$Q$ pairs rejected;
   distinct-factor pair accepted as no-common-component.
8. Route checks: direct and four-factor routes agree on a plant; a corrupted
   factor evaluation breaks the agreement; a repeated-projection 12-set
   fails the all-distinct gate.
9. Generator control: exhaustive enumeration of compatible pairs reaches
   exactly 41,587,200 and equals the closed form $2\binom{96}{2}^2$; first
   compatible pair is $(1,97)$-indexed rank $(0,96)$ semantics — pinned by
   direct computation.

## 4. Registered outcomes and verdict rule

- **E-positive:** at least one qualifying GF(17) 12-set exists and all its
  residual deletions are classified. Verdict **FROZEN-CERTIFIED** for the
  completed exhaustive pool classification (the classification is complete;
  its scope remains the finite pool).
- **E-negative:** some qualifying rank-10 11-set is nonminimal (fails a
  10-deletion rank). This falsifies no-fixed-component sufficiency on the
  full registered pool. Verdict **FROZEN-NEGATIVE** with the full witness;
  all other classifications still recorded.
- **E-empty:** no qualifying 12-set occurs anywhere in the full compatible
  population. Because the sweep is exhaustive over the registered family,
  this is a certified finite emptiness. Verdict **FROZEN-CERTIFIED**,
  scoped: "the registered 96x96 factored GF(17) pool admits no qualifying
  all-distinct residual configuration". This is not universal C3
  nonexistence; Lemma R1 stays a separate analytic subresult; any claim
  beyond the pool requires a new preregistration.
- Incomplete, failed, or budget-exceeded runs: **FROZEN-INCONCLUSIVE**.
- Instrument defect: uniquely named broken artifact; nothing deleted; the
  affected stage reruns from the beginning with assertions intact; the
  defect is logged before any rerun.

## 5. Runtime and supervision

`run_fullpool.py` is standalone: exact Python integers modulo 17, no floats,
no NumPy/SymPy, no import of frozen campaign code (vendored and modified from
the frozen POINT driver, whose SHA-256
`1aee81d5bd840c4332564e7307e857113d53a47f6ab600d0004b9dfe4aa5090a` is
recorded in provenance). `sys.dont_write_bytecode=True` before other imports;
launched with `PYTHONDONTWRITEBYTECODE=1`; all five thread caps equal 1;
`nice -n 10`; `RLIMIT_CPU=(5400,5400)` asserted in-run. Hard caps: 5400 CPU
seconds, 6000 wall seconds. Fixed-chunk checkpoints every 4,000,000 checked
pairs and on every qualifying base set; cumulative stats and the size
histogram are persisted at every checkpoint so a partial run has an exact
record. Expected cost from the frozen POINT run is about 12.4 CPU
seconds per million pairs, i.e. roughly 9-17 CPU minutes for the full sweep
— well inside the caps, but the run is launched under process supervision
and is resumable only by a fresh full rerun (no mid-sweep shortcut is
permitted).

## 6. Independent replay

Before close, `independent_residual_audit.py`-style bitset replay is
re-run against the full sweep: the pool SHA, the per-size histogram over all
41,587,200 pairs, the four size-12 anchors, and the classification of every
qualifying set must reproduce from the preregistered seeds alone, driver-free.
The replay artifacts are hashed into the freeze ledger.
