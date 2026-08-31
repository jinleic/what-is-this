# VERDICT.md — S3 campaign: three-slice L-family bound at n=8
Campaign: `2026-08-30T19:27:14Z_d6ca8a53-0950-4a3e-a6ef-084d0a1f39d0_S3`
Pre-registration: `pre_statement.md` in this directory, committed BEFORE
the first computation. All scripts frozen here as byte-copies with SHA256
in `code_hashes.sha256`; outputs in `*.out` (all exit 0).

## Outcome first

**NEGATIVE (no reach of 14).** No route established
`rank(L_u, L_v, L_w) >= 14` for the declared class (ALL linearly
independent (u, v, w) in R^8). The best certified lower bound REMAINS the
chain's own:

    rank(L_u, L_v, L_w) >= 13   for all independent u, v, w   (see below)

with the *gap to 14 exactly 1*, and the specific-triple picture refined:
for the special triple (e0, e1, e2) = (1, i, j), an exact 14-term
block construction gives `rank((L_1, L_i, L_j)) <= 14`, so that triple
sits in [13, 14] — the worst-case question is NOT settled by any
single-triple upper bound. **No claim of 18 -> 19 is made; S3 stays open.**

## Anchors (all MACHINE-VERIFIED, exact fmpq, run first)

| anchor | statement | verdict | evidence |
|---|---|---|---|
| 1 | R_R(T_C) = 3 at n=2: upper 3-term witness entrywise-exact; lower core (U A = I, U B = J, rank facts) exact | PASS | s3_anchor.out |
| 2 | R_R(T_H) = 8 at n=4: upper = frozen gate-C 8-term witness, entrywise-exact; lower chain identities exact (L_i^2 = -I, L_a L_b = L_ab (assoc), C = L_{-i} L_j = L_{-k}, C^2 = -I) | PASS | s3_anchor.out |
| 3 (CONTROL) | R_R(tau) = 7 at n=4, k=3: upper via tau_r7 Krawczyk certificate (existence at rho = 1e-5, exact-replayed: K = 0.005352 < 1, worst margin 9.946e-6 > 0); lower = 1 + 6 chain | PASS (as a control) | s3_anchor.out + upstream replay log |

Anchor-3 is the rule-14 control: ANY general S3 machinery that would
output >= 8 for a 3-slice family at n=4 is refuted by tau and is dead on
arrival.

### Landmine found during anchoring (MACHINE-VERIFIED observation)

The tau_r7 / rank25 certificates store DECIMAL-LITERAL factors (up to 19
significant digits). Parsing them through float64 destroys ~4 digits
(e.g. stored `7.861952407189559677e-01` -> float64 rational
differs from the exact decimal by ~5e-19 *per entry*, accumulating to
~1e-16 residuals after 7-term contraction). The gate A/C replays were
unaffected only because their tolerance floors were 1e-6 and the asserted
inequalities carry 1e-5+ margins.切记: certificate decimals must be parsed
EXACTLY (Fraction(str)) whenever entrywise-exactness is demanded — and the
certificate's correct semantic is "an exact witness exists within radius
rho of the stored numeric factors", NOT "the stored decimals are an exact
witness".

## Route-by-route numeric outcomes (misses prominent)

Domain (fixed pre-run): all independent triples in R^8; reported bounds
are worst-case floors over the swept families.

| route | method | numeric outcome | verdict on reaching 14 |
|---|---|---|---|
| A (Strassen commutator, both formula variants) | rank([L_ubar L_v, L_ubar L_w]) exact over 56 basis + 112 structured + 5 fixed triples | commutator rank = 8 on EVERY swept member (min=max=8); twin forms: 8 + 8 = 16 (Strassen R >= n + rank form) and 8 + 4 = 12 (pencil /2 form) | **12 or 16, depending on which variant is valid**; see obstruction below |
| B (AFT/substitution pivots) | pivot-independence verified exactly for all 56 basis triples x 3 pivots; double-peel and flattening variants computed | all pivots give 1 + 12 = 13; double-peel 2 + 8 = 10; flattening 8 | misses; best = 13 (same as chain) |
| C (residual-pencil floor) | exact re-derivation that residual pairs stay independent and satisfy the irreducible quadratic (C^2 - 2aC + bI = 0, a^2 < b verified exactly per triple); C4's feed: pencil floor = 12 exactly on the consumed class | 1 + 12 = 13, cannot exceed 13 within the peel-into-pencil route | closed at 13 [MACHINE-VERIFIED route ceiling] |
| D (flattenings + Griesser combo) | exact flattening ranks: mode-(p) = 3 on all 56 basis triples; mode-(c), mode-(b) = 8; combo max(8 + comm) = 16 | 3 / 8 / 8 << 14 | miss (cap 8); combo inherits Route A's 16-with-obstruction |
| E (exact ideal-theoretic infeasibility) | pre-registered verdict BUDGET-FAIL; only the linear-feasibility diagnostic for ONE fixed generic factor family ran (192 eqs, 39 unknowns: inconsistent, as expected — necessary-only diagnostic, NOT a bound) | no bound | not attempted to completion; recorded honestly |
| F (sub-algebra block decomposition — genuinely NEW content) | EXACT: span{e0..e3} is a quaternion subalgebra, orthogonal complement span{e4..e7}; for x in H', L_x is block DIAG (+) DIAG (complement action conjugate-equivalent to Lq_x — exact per-triple check); the specific triple (1, i, j) is conjugate-equivalent to blockdiag(tau, conj tau) | **upper side: rank((L_1, L_i, L_j)) <= 14** (7 + 7 via block-concatenation of two copies of the machine-verified tau facts); lower side for THAT triple: 13 (chain) | specific triple in [13, 14]; worst case untouched |

* RETRACTION (rule 5), recorded inline: the first Route-D script
  (`s3_routeD_flattening.py`, first commit) computed the mode-(c)
  flattening with a slicing error, reporting flattening rank 1 for
  triples containing e_0. Correct value: 3 / 8 / 8 on ALL 56 basis
  triples (fixed script, re-run, re-hashed; the wrong numbers appear
  only in the superseded output, retracted here and in scratch log).
* RETRACTION (rule 5), recorded inline: the first Route-A script
  (`s3_routeA_strassen.py`) printed a pre-formula bookkeeping draft with
  the WRONG derived number (16 from an invalid pairing of the /2 rule).
  The database numbers are Route A's honest twin-form report in
  `s3_routeAB.py`: commutator rank 8 everywhere on the swept families;
  8 + 8 = 16 under the R >= n + rank([X, Y]) reading, 12 under the
  pencil /2 reading. The superseded script is retained in the freeze
  (byte-copy) with its WRONG output marked; the corrected artifact is
  s3_routeAB.py + .out.

## The decisive obstruction (why 16 is NOT a certified 3-family floor)

The Strassen commutator bound in the form R >= n + rank([A^{-1}B, A^{-1}C])
applies to the 3-SLICE system (A, B, C) as a tensorial object. For an
L-family, slices (L_u, L_v, L_w) satisfy the EXTRA quadratic identities
(the Lean chain's H3: L's satisfy L_ubar L_v's minimal quadratic with
negative discriminant), which the stripped (A, B, C) form does not consume.
The exact computation establishes:

* [X, Y] has rank exactly 8 = 2*(n/2)?? — the swept value is 8 for every
  member of F1/F3/F4 (56 + 112 + 5 triples), including all basis triples.
  rank([X, Y]) = 8 gives via n + rank = 16 the twin reading; but a
  CERTIFIED universal 16 needs:
    (a) the commutator rank to be >= 8 for EVERY independent triple (we
        swept 173 triples with min = max = 8 — a strong pattern, but per
        rule 7/17a a swept pattern is EQUIVALENT-CLASS evidence, NOT the
        universal quantifier [COMPUTATIONAL-EVIDENCE for the universal
        statement]), and
    (b) the Strassen form feeding n + rank (not n + rank/2) to hold for
        3-slice systems — cited-variant verification pending
        ([CITED-DEPENDENCY UNVERIFIED]: Strassen 1974 / Blaeser 2003 /
        Lickteig 1984 statements not read first-hand this session;
        the /2 split is correct for (x, y)-type pairings).
* The tau CONTROL pins the space: at n = 4, tau's true rank is exactly 7,
  while the twin readings would give 4 + 4 = 8 (Strassen: n + rank) or
  4 + 2 = 6 (pencil: n + rank/2). BOTH MISS in opposite directions: any
  certifiable general machinery must output EXACTLY 7 at (n=4, k=3) —
  the commutator machinery does not. So Route A's n=8 number 16 is
  UNRELIABLE as a certified floor (too high by construction to trust
  without the exact-form (b) resolution), and Route C's accounting nails
  the certified floor back at 1 + 12 = 13.

## What this campaign certifies (rule 15/16 discipline)

* [MACHINE-VERIFIED] Exact anchors 1-3 (n=2 => 3, n=4 => 8, tau = 7) on
  freshly re-derived scripts (upper witnesses entrywise; lower chain
  identity core exact; Lean floor [CITED-DEPENDENCY]).
* [MACHINE-VERIFIED] Route C: the peel-into-pencil route cannot exceed
  13 (independence of residuals + irreducible quadratic, exact per
  triple; pencil floor 12 exact from C4's frozen machine-verified
  witness).
* [MACHINE-VERIFIED] Route B: all pivot bookkeepings give 13; nothing
  above.
* [MACHINE-VERIFIED] Route F: (1, i, j) triple decomposes into two
  conjugate copies of tau: rank <= 14 for that triple exactly;
  rank >= 13 for it by the chain. Worst case: unresolved, gap 1.
* [COMPUTATIONAL-EVIDENCE] Commutator rank = 8 across 173/173 swept
  triples (universal statement NOT certified; family fixed pre-run).
* [CITED-DEPENDENCY UNVERIFIED] The Strassen-twin form validity for
  3-slice systems; if first-hand-verified, a NEW question opens
  (could 16 hold universally? — tau control refutes the analogous n=4
  claim, so the n=8 form must be checked against tau-logic before any
  escalation).

## Best certified bound and exact route

    rank(L_u, L_v, L_w) >= 13    for all linearly independent u, v, w
    Route: 1 (substitution peel — Lean Substitution.lean, exact) + 12
    (Thm-2 pencil floor at n=8, C4 machine-verified tight on the consumed
    class).

Gap to 14: **exactly 1**.

## Re-scoping log (rule 16)

* ADDENDUM (not a re-scope; recorded for rule 14 completeness): Route F
  (sub-algebra block decomposition) was NOT in the original 5-route
  pre-registration; it was added AFTER Routes A-E produced their numbers,
  as the upper-side attack. NO bound in this campaign depends on Route F;
  its outputs (upper <= 14 on the specific triple) are reported with
  their own labels. Original routes A-E stand as pre-registered with
  domains untouched.
* The E budget (30 min) was not consumed (diagnostic ran in 0.25 s);
  verdict BUDGET-FAIL per pre-registration (far short of a complete
  Groebner infeasibility, which the pre-statement already expected
  unreachable: 312 unknowns).

## Rule-7 scope sentence

This campaign covered: the three-slice L-family lower-bound question S3 at
n = 8 over the domain (ALL linearly independent (u, v, w) in R^8), via
five pre-registered routes (Strassen commutator twins, AFT substitution
pivots, residual-pencil floor, flattening/Griesser, exact Groebner
— budget-failed) plus one post-hoc upper-side route (quaternion
sub-algebra block decomposition), sweeping 173 triples for Route A, all
56 basis triples for Routes B/C/D, one fixed factor family for E,
and the global-conjugation structure for F; all arithmetic fmpq-exact
(no floats in any load-bearing step). NOT covered, hence NOT excluded by
any claim here: a universal commutator-rank >= 8 statement for all
independent triples (only 173 swept — equivalence-class evidence only),
the first-hand verification of the Strassen/Bläser/Lickteig twin-form
validity for 3-slice systems, any upper bound on the whole T_O 8-slice
tensor or R_R(T_O) itself (the published window 18 <= R_R(T_O) <= 25 is
unchanged), any non-L-family 3-slice residual (peels with non-linear
Phi), any rank statement at n != 4, 8, and any claim that the specific
triple (1, i, j) attains either endpoint of [13, 14].
