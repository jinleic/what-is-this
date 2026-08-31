# Pre-statement AMENDMENT 1 — honest sandwich map corrected by the A2 control
**Timestamp 2026-08-30T20:0xZ. Committed BEFORE any census compute (pass B/C untouched;
zero census numbers exist). Amends pre_statement.md §2/§3; original left intact per rule 16.**

## What happened
Pre-statement §3 declared the honest action as `u' = Pᵀ·U·Q⁻ᵀ, v' = Qᵀ·V·R⁻ᵀ, w' = Rᵀ·W·P⁻ᵀ`
(the naive cyclic-conjugation family, derived from `tr(ABC)` invariance). The A2 control
(§4, 10 seeded non-monomial triples) REFUTED it: Brent failures up to 728, and monomial
mismatch with the frozen code. Per the pre-registered rule, compute was halted and the
derivation redone.

## Root cause
`tr(ABC)` is the wrong invariant. The mm3 tensor's trilinear form is
`Phi(U,V,W-decomposition) = Σ_{i,k,j} u[(i,k)] v[(k,j)] w[(i,j)]` with delta structure
`[i==ip][k==kp][j==jp]` tying A-rows to OUTPUT rows (not to a third matrix-space),
A-cols to B-rows, B-cols to output cols.

## The corrected honest action (delta-collapse solved exactly)
For ANY invertible 3×3 P,Q,R over ℤ, the factor map

    u2_blockwise = P^{-T} · U_blk · Q^{-1}
    v2_blockwise = Q      · V_blk · R^{-1}
    w2_blockwise = P      · W_blk · R^{T}

preserves the tensor EXACTLY. Evidence, in order of increasing strength:
1. **Symbolic delta-collapse:** the three pairwise contractions reduce to
   `Pinvᵀ·Pᵀ = I`, `Q·Qinv = I`, `R·Rinv = I` — verified by sympy with fully symbolic
   P,Q,R entries (identity for ALL invertible matrices, no det/orthogonality condition).
2. **Frozen-limit agreement:** at (P,Q,R)=(X,Y,Z) arbitrary MONOMIAL, this map equals
   the frozen `gatec_sweep.sandwich` output bit-for-bit on all three 23×9 blocks
   (measured, MONO[5],MONO[6],MONO[7] trial).
3. **Empirical: 45/45 random non-monomial (P,Q,R) triples across BOTH `paper55` and
   `sun56`: 729/729 Brent over ℤ, zero failures** (seed 20260901 battery + earlier
   20260830 battery on paper55 post-correction).

Interpretive note (scope, not a paper defect — nothing escalated): the *published*
action sentence `(X U Y⁻¹, Y V Z⁻¹, Z W X⁻¹)` matches this family only when read with
contragredient letter pairing; read with three independent letters it is an automorphism
exactly on the orthogonal subgroup — the same shape as the frozen-code defect documented
in pre_statement.md §0. All published results remain true under their own (monomial/orbit)
conventions.

## Consequences for the census design (all положения unchanged except the letters)
- Diagonal slice under the corrected map: `u2 = P^{-T} U P^{-1}, v2 = P V P^{-1},
  w2 = P W Pᵀ`. The A1 dual-map anchor will be RE-RUN with this map as the third variant;
  the 6912/0/48 target must reproduce (it holds for both prior maps by relabeling).
- Right-monomial factorization under the corrected letters: `P = a·p` (data×monomial):
  `u2 = p^{-T}(a^{-T} U a'^{-1}) q'` — permutations factor OUT exactly; the alphabet
  predicate and bound layer depend only on DATA-triples. Per-decomposition arithmetic:
  A = 1160 data nodes (determinant ±1 split recorded), pair tables on A², triangle count
  on A³, ×(monomial multiplicity)×3σ. The exact multiplicative bridge must reproduce
  BOTH frozen touchstones before any off-diagonal number is released:
    (i) 48³ = 110,592 per (D,σ) monomial-only orientations
    (ii) 48³×3 = 331,776 per D across σ
- A2 battery EXTENDED (power upgrade per Main): non-monomial triples now sampled across
  both decompositions and both det branches; acceptance = 45/45 plus the symbolic proof.
  This supersedes §3's 10-triple spec (which caught the bug at 10/10 fail — its power
  is proven by the catch).
- Wall-caps, escalation, adjudication (§1, §5, §7, §9): UNCHANGED.
- Falsifiable outcome unchanged (§1): exact survivor count per (D,σ) as primary
  deliverable; ≤54 escalation rule unchanged.

## Retraction recorded (rule 5)
The §3 sentence "honest cyclic action `u' = Pᵀ·U·Q⁻ᵀ …`" is RETRACTED and replaced by
the corrected family above. The §0 statement that the frozen diagonal formula equals the
honest map at `(P,Q,R)=(G⁻¹,G⁻¹,G⁻¹)` is likewise RETRACTED as stated (it held for the
retracted family); the corrected diagonal claim — frozen@G equals honest@G on all 48
monomials and the census verdicts 6912/0/48 agree under both maps on all 6960 — is the
A1-measured statement and was measured BEFORE this amendment (pass A part 1), so the
anchor layers of the campaign were NOT invalidated by the derivation error: the frozen
code was used as the frozen code, and the corrected map is new work.
