# Pre-statement addendum S3 — three-slice L-family bound (agent OctRankS3) [DECIDED AHEAD, 2026-08-30]

**(Committed BEFORE any S3 computation ran. This is the pre-registration the
assignment requires for item 1; thresholds, routes, budgets and the
adjudication rule below were fixed after reading cs/README.md rules 1-17c,
the oct-rank README (full), pre_statement.md clauses C3 and the C4/C3
addendum, the C4 VERDICT, and the Lean chain Oct18/Peel/Pencil/Tau7 — and
BEFORE running any new script. No thresholds may be adjusted after seeing
results; any re-scoping gets an appended "Re-scoping log" entry keeping the
original visible. Rule 16: fixed domains, no domain-shopping. Rule 5:
retractions inline. Rule 14-17c: every route below must be able to ESTABLISH
its conclusion, or it is not evidence and is labelled as such.)

## The object and the question (rule 17c: provenance fixed up front)

Object fixed for this entire campaign: T_O = 8x8x8 real structure tensor of
octonion multiplication built by Cayley-Dickson doubling
H(+) H l with basis (1, i, j, k, l, il, jl, kl),
(a,b)(c,d) = (ac - conj(d)*b, da + b*conj(c)) — the SAME convention as
gate_a_verify.py A1 (entries in {-1,0,1}; cross-checked entrywise against
the upstream table by gate A). Slice p of T_O is the 8x8 left-multiplication
matrix L_{e_p}. A "3-slice L-family" is (L_u, L_v, L_w) for linearly
independent u, v, w in R^8 (independence in the vector sense).

S3 question, verbatim from the target README honest-gaps section: is
rank(L_u, L_v, L_w) >= 14 for independent u, v, w in R^8? The chain gives
5 + 13 = 18 via the k-family induction (Oct18Peel.slicesRank_L_family_ge:
k-slice L-family of independent vectors costs >= n + n/2 + k = 13 at k = 3, n
= 8); a 3-family floor of 14 would improve 18 to 19.

Published context [CITED-DEPENDENCY, owner-read per target README]:
arXiv:2608.16649 (Hardik Jain) proves 18 <= R_R(T_O) <= 25 via peeling plus
the n + n/2 pencil floor; sharpening the three-slice floor 13 to 14 would
lift the published 18 to 19. Escalation rule (assignment + pre-statement
C3): any candidate >= 14 is NOT written into any authoritative file; freeze
evidence, hub-message Main, stop.

## What a failure to reach 14 certifies (fixed IN ADVANCE)

If no route below establishes rank(L_u, L_v, L_w) >= 14 for the declared
class, the honest outcome is: "best certified lower bound on
rank(L_u, L_v, L_w) over independent u, v, w is B with exact route R", plus
an explicit statement of which proofs were attempted and what each gives.
A certified 12 or 13 with a named obstruction is a legitimate, reportable
closing of the item. It is NOT dressable as anything else: no claim about
R_R(T_O) >= 18+1/2 or similar is made; the only statement is about the
worst case of the 3-slice floor over the declared domain
(independent u, v, w in R^8).

## Exactness discipline (rule: no floats in load-bearing steps)

- Every load-bearing rank/linear-algebra computation runs in sympy over
  Q or python-flint fmpq. Floats only for pre-computational reconnaissance,
  and every float-derived claim is re-established exactly before use, or
  it is labelled COMPUTATIONAL-EVIDENCE and does not enter any bound.
- Arb appears nowhere in this campaign unless an interval refinement of an
  exact quantity is explicitly needed (not planned).
- Solver-tolerance trap (README constraint): no LP/MIP solver output is
  used as a bound directly; if an LP relaxation is ever used to HUNT for
  structure, its output is re-derived exactly in closed form; reported
  tolerances are compared against the magnitude sought.

## Anchors to establish FIRST (before any new claim)

Anchor-1 (n=2 => 3, complex): the machinery reproduces R_R(T_C) = 3.
  Already machine-verified in gate C (both directions, self-contained).
  Re-run here on MY exact scripts (cheap) as a harness check — the anchor
  validates the harness, not the published fact.
  Pass: s3_anchor_n2.py exits 0 with exact identity checks.
Anchor-2 (n=4 => 8, quaternion): the machinery reproduces R_R(T_H) = 8.
  Upper 8 = hadamard+sparse witness (exact, prior campaign); lower
  peel(2)+pencil(6) [CITED-DEPENDENCY Lean chain].
  Re-run the exact witness + a fresh exact-restatement of the chain bookkeeping.
  Pass: s3_anchor_n4.py exits 0 with witness equality exact.
Anchor-3 (NEW, the load-bearing one for S3, rule 14: test what the instrument
  CANNOT do): the 3-slice floor at n = 4 must NOT be 7 in a fixed witness —
  the THREE-slice quaternion tensor tau = (L_1, L_i, L_j) on H has
  rank EXACTLY 7 [MACHINE-VERIFIED earlier this campaign, both directions].
  This is the CONTROL for any claim of the form "the 3-family floor is
  n + n/2 + k + ceil(n/k)-ish": any general bounding scheme that would
  promise >= 8 at k = 3, n = 4 is REFUTED BY ANCHOR-3 and cannot be the S3
  route. Any proposed S3 route is first checked against anchor-3; a route
  that would also lift tau to 8 is dead on arrival and labeled as such.
  (Anchors 1-2 pin the pencil; anchor 3 pins the 3-family floor at n = 4.)

## Fixed domain for all routes (rule 15/16)

Domain: all (L_u, L_v, L_w) with u, v, w linearly independent in R^8.
Worst-case floor meaning: the largest B such that EVERY independent triple
has rank >= B. This is the object the chain consumes. NOT in scope for any
bound claim in this campaign: rank conditionals on non-generic triples
(that would be a different, specific-triple statement, of no use to the
chain), upper bounds except to orient, and the n = 8 8-slice object T_O
itself beyond what the 3-family statement feeds it.

## Routes (pre-registered; every one reported numerically, misses prominent)

Route A — Strassen commutator: with A = L_u (u != 0, invertible),
  R >= n + (1/2)*rank([A^-1 B, A^-1 C]) = 8 + (1/2)*rank(commutator of
  L_{u-bar v}, L_{u-bar w}) computed EXACTLY as a matrix rank over Q.
  Exact-control for this route: on tau (n = 4), the commutator of
  L_{u-bar v}, L_{u-bar w} has rank 2 = n/2 for u = 1, v = i, w = j and the
  route gives 4 + 1 = 5 — a WEAKER number than the true 7; recorded so the
  route's n = 8 value is never misread as tight-by-default.
  [Concretely: for A = L_1 = I, B = L_i, C = L_j: [B, C] = L_{ij} - L_{ji} =
  2 L_k, rank 4 on H; formula normalized gives 4 + 2 = 6 on tau... the
  EXACT tau number is fixed by the script; the pre-committed expectation
  recorded here is that the route on tau yields EXACTLY 6 < 7, i.e. the
  route is lossy. If the script's tau-control value differs from this
  derivation, the discrepancy is reported inline (rule 5) — the pre-statement
  guess is a hypothesis, per the owner-calibration directive.]
  n = 8 analogue: u, v, w the darkest independent triples the chain allows.
  The commutator rank is computed exactly for a spanning family of triples
  (u, v, w); the REPORTED number is the minimum over the spanning family.

Route B — Blaeser/Lickteig-type improvements over Strassen for 3-slice
  tensors: the known theorem R >= n + rank([A^-1 B, A^-1 C]) for
  A, B, C with [B, C] having rank > ... (Lickteig 1984; Blaeser's formula
  for (n, m, 2) three-slice tensors). Implementation: exact matrix-rank
  computations; the route is expected to give n + rank([B, C]) = 8 + 4 = 12
  worst-case on triples where the commutator has rank 4 in a nested
  algebra — i.e. <= 12, ABOVE the machinery's 13, dead end, but computed
  and reported with labels.

Route C — Substitution-method (AFT) floor on the 3-family directly: instead
  of peel-1 + pencil, analyze rank((L_u, L_v, L_w)) >= n + pencilFloor of
  the residual pair (L_v - lam L_u, L_w - mu L_u) — the same value 1 + 12
  = 13 the machinery gives; the route tests whether a residual pencil can
  EXCEED the n + n/2 floor when its normalized complex structure is
  restricted to the octonion L-family. That is the S1/C4 question in
  pencil form; C4 proved the consumed machine gives 12 exactly on
  (I_8, J_8); the route's only NEW content would be a J outside the
  conjugacy orbit of J_8 — but every J with J^2 = -1 on R^8 is GL-conjugate
  to J_8 (standard canonical form; [INFERENCE] flagged in C4), so the
  route is expected to close at 13 with a machine-checked reason: the
  pencil-floor-12 is TIGHT on the consumed class. Report: the exact
  arithmetic showing (a) C^-1 D has irreducible quadratic (chain
  identities), (b) the floor-12 witness of C4 extends to a rank-(13)-1
  ... [concretely: 1 + 12], (c) hence the substitution route cannot exceed
  13 as a universal statement.

Route D — Laderman/Griesser-type equation counting, border-rank
  obstructions, flattening ranks: flattening ranks are capped by
  min(64, 64, 24) = 24 for (L_u, L_v, L_w) as an 8x8x3 tensor — wait,
  rule 5 retraction for the pre-statement: the tensor (L_u, L_v, L_w) is
  8 x (8x3) = 8 x 24 under (u) vs (v, w) partition — flattening rank
  <= 8; the (L_u L_v, L_w) partition is not a flattening of a 3-slice in
  slice-rank sense; the useful flattening is the (u) x (v, w) one and it
  caps at 8, i.e. NO, flattening gives at most 8 << 14; recorded as the
  expected miss with the exact number 8. Laderman-type lower bounds
  (equation counting on 3x3x3-like cells) do not apply at n = 8 3-slice
  with these dimensions; border-rank obstructions over C (e.g. commutator
  with a nilpotent) — the ./scripts compute the organic candidates
  (rank of the 2x(8x8) + 1x(8x8) flattening = 8; commutator rank facts)
  and every miss is recorded.

Route E — exact ideal-theoretic infeasibility: decide RANK >= 14 by
  establishing NO rank-13 slice decomposition exists for SOME specific
  (L_u, L_v, L_w). Feasible only if the variety is zero-dimensional or
  has small Gröbner basis; pre-committed budget: 30 min wall-clock on ONE
  chosen triple (u, v, w) = (e0, e1, e2) [the darkest structural triple:
  1, i, j, whose residual pencil is the standard complex structure J_8 on
  H-blocks]; Gröbner basis over Q for 24 unknowns is UNDECIDABLE in
  general in that budget; pre-commit the expected verdict: BUDGET-FAIL
  (= not attempted to completion, NOT a certified negative, recorded as
  such with the exact sympy/groebner call and failure mode).

## Budget (fixed)

- Wall-clock cap: 4 h total across all routes (this session).
- Anchor runs: not counted (seconds each, on prior campaign artifacts).
- Route E sub-budget: 30 min on the single pre-declared triple, then halt.
- No solver starts unless its tolerance is below the quantity sought by
  >= 100x; no solver output counts as evidence without exact re-derivation.

## Adjudication rule (fixed)

The campaign outcome is exactly one of:
  (i) CANDIDATE: some route establishes rank(L_u, L_v, L_w) >= 14 with an
      exact-arithmetic certificate whose argument survives self-audit.
      Action: freeze evidence, hub-message Main, STOP (no README claim).
  (ii) NEGATIVE: no route reaches 14. Action: report best certified B with
      exact route; classify each route's numeric outcome incl. misses;
      append README section with rule-7 scope; file in campaigns/ +
      scratch/ per C3 grade of the artifacts.
A route's output that refines an upper bound or gives a specific-triple
bound (generic-only) does NOT count as (i) or close the item: the chain
consumes the WORST-CASE floor over independent triples only.

## Rule-7 scope sentence (pre-drafted skeleton, final numbers filled at close)

"Covered: [exact routes A-E with their certified numeric outputs on the
declared domain, ALL independent (u, v, w) triples]. Not covered:
[enumeration of what was not searched — non-linear Phi, other n, upper
side beyond orientation, specific-triple rank-14 statements,
non-L-family 3-slice residuals arising from non-L peeling — and no claim
about R_R(T_O) itself beyond the existing 18]."
