# Theorem (corrected form) — the dual product-code mechanism of the
# weight-≤3 window emptiness, first-hand derivation + machine anchors

Campaign `2026-09-01T02-33-35Z_mech2_mechanismCorrected` (Campaign M). This
file is the §4 deliverable of the pre-statement (git `29b0dfe`, committed
pre-compute; runner frozen as-committed `e117ae2`, sha256 `8c1e9251…`,
working copy sha256 at freeze `2919905c…` — the two documented rule-5 fixes
below are in the working copy; see §F). Provenance: the quotient skeleton —
tensor code as a sum of C_i-slices, coset independence, separator,
sharpness — **originates with the owner** (frozen campaign
`2026-08-31T08-26-54Z_hgateMechanism`, credited; owner review corrections 1-7
are the reason this campaign exists). The derivation below is RE-DERIVED
first-hand in this campaign's own words in the owner-supplied CORRECTED form
(no per-support biconditional, C-2 relation-vector form, C-4 dimension form,
C-5 set-only triple-sum). No mechanism novelty is claimed by this agent.

Label: **HUMAN-AUDITED derivation (independently re-derived in corrected
form; owner skeleton credited) + MACHINE-VERIFIED controls (exact F_q
arithmetic end to end; 16/16 controls PASS, zero deviation).**

------------------------------------------------------------
## (a) Hypotheses and scope

- q prime (all checks run on prime fields; nothing in L1-L6 uses primality,
  only field structure — prime-power fields are in scope by the same proofs
  but were NOT swept by machine here).
- s = (s0,s1,s2), s_i ≥ 1; t = (t0,t1,t2) with **1 ≤ t_i ≤ s_i** (H1; t_i = 0
  is out of scope — the zero code, d = ∞ by convention, never silently
  treated as s_i + 1).
- Evaluation sets S_i ⊂ F_q, |S_i| = s_i, **pairwise distinct points within
  each S_i** (H2; subgroup elements qualify). NO coprimality condition on
  the s_i enters anywhere in L1-L6.
- Λ_i = (λ_i(a))_{a∈S_i} with all entries nonzero (H3) — the diagonal
  scaling is support-preserving per coordinate.
- A_i = F_q^{s_i}; C_i = { (λ_i(a)·p(a))_{a∈S_i} : p ∈ F_q[x], deg p < t_i };
  N = s0 s1 s2; points indexed by idx = (a0 s1 + a1) s2 + a2.
- d_i := d(C_i) = s_i − t_i + 1 (L4); d := min_i d_i.
- V ⊂ F_q^N := Σ_i A_0⊗…⊗C_i⊗…⊗A_2 (L1), the code instrumented by
  `DeltaEngine.V_basis` (rref of the stacked line classes).
- For a point support S ⊆ {0..N−1}: V ∩ F^S = {v ∈ V : supp(v) ⊆ S}.

## (b) Statement (correction-integrated, C-1 form)

**Theorem T.** Under (a):
1. For EVERY point support S with |S| < d: V ∩ F^S = {0}.
2. There EXISTS a support S with |S| = d and V ∩ F^S ≠ {0}: take
   i* ∈ argmin d_i; a minimum-distance word of C_i* (L4) embedded on one
   direction-i* line is a V-codeword supported on d_i* = d points.
3. No claim is made for an individual S with |S| ≥ d beyond 2 (C-1: the
   per-support biconditional is FALSE in general).
Census semantics (the union window): {S : |S| ≤ w, V ∩ F^S ≠ {0}} is EMPTY
iff w < d; its first nonempty weight is exactly d; for larger w the window
is nonempty (monotone in w by inclusion) — nonemptiness for w ≥ d follows
from 2 alone.
**Corollary (t = (1,1,1) slice; d_i = s_i, d = min_i s_i):** the weight-≤3
window is empty iff min_i s_i ≥ 4; nonempty iff min_i s_i ≤ 3 — H-GATE's
min-form, as a theorem on this slice (including the s_i = 1 corners: d = 1
= that s_i, C_i = F_q^1 = A_i, and V = A_0⊗A_1⊗A_2 there; corner behavior
machine-anchored by E1corners). The previously-prospective (4,5,7) class is
inside the theorem: d = 4 ⟹ w=3 empty, w=4 nonempty — both sides first-hand
machine-anchored here (C-PN4: 457,450 supports EMPTY; E2p: 35/35 direction-0
4-lines carry codewords).

## (c) The derivation, line by line

Notation: for v ∈ A_i, supp(v) = {a : v_a ≠ 0}; wt(v) = |supp(v)|. π_i :
A_i → Q_i := A_i/C_i. "Line" = All axis-parallel fiber; direction-i line L
fixes the other two coordinates; the lifted class of a word c ∈ C_i on line
L is the vector w ∈ F_q^N with w on L = c and 0 elsewhere.

- **L1 (V as lifted slices).** The direction-i line space on one line is
  {c embedded on L : c ∈ C_i} (at t_i = 1 this is span(all-ones on L);
  in general the Λ-scaled GRS words on L). V := span over ALL lines of ALL
  directions of these classes. Under idx, embedding on a direction-i line is
  A_i ⊗ (coordinate basis of the other axes), so
  V = Σ_i A_0 ⊗ … ⊗ C_i ⊗ … ⊗ A_2 as subspaces of A_0⊗A_1⊗A_2 = F_q^N.
  This matches the implemented generator stack (lift_rref(i) = rref of all
  direction-i line classes).
- **L2 (kernel identity).** The tensor product is right-exact: with π =
  π_0 ⊗ π_1 ⊗ π_2 : A_0⊗A_1⊗A_2 → Q_0⊗Q_1⊗Q_2,
  ker π = Σ_i A_0⊗…⊗ker π_i⊗…⊗A_2 = Σ_i A_0⊗…⊗C_i⊗…⊗A_2 = V.
  (Right-exactness of ⊗ concrete: LHS ⊇ each summand since π kills C_i in
  slot i; equality by dimension count below.) By rank-nullity and the
  pure-tensor span of the images, rank π = Π_i rank(π_i) = Π_i (s_i − t_i)
  (rank π_i = s_i − dim C_i = s_i − t_i: the evaluation map p ↦ (λ p(a)) has
  kernel {deg < t_i polys vanishing on all s_i distinct points} = 0, so
  dim C_i = t_i by H1/H2/H3). Hence
  **dim V = dim ker π = N − Π_i (s_i − t_i)** [C-4 form]. At t = (1,1,1):
  dim V = N − Π(s_i − 1).
- **L3 (coset independence; C-2 form).** Claim: for T ⊆ {0..s_i−1} with
  |T| < d_i, the cosets {π_i(e_a) : a ∈ T} are F_q-independent in Q_i.
  A relation Σ_{a∈T} c_a π_i(e_a) = 0 means c := Σ_{a∈T} c_a e_a ∈ C_i
  (preimage nonempty ⟹ its π-image zero... precisely: c ∈ ker π_i = C_i).
  Now supp(c) ⊆ T (the c_a are the only nonzero slots — λ_i-scaling is
  per-coordinate nonzero multiplication and does not spread support), so
  wt(c) ≤ |T| < d_i = d(C_i), forcing c = 0 by the minimum distance, so all
  c_a = 0. ∎ [C-2 honored: NO use is made of c having equal/indicator
  entries; the argument never claims 1_T ∈ C_i.]
- **L4 (MDS distance).** d(C_i) = s_i − t_i + 1.
  (≤) Fix U ⊆ {0..s_i−1} with |U| = t_i − 1. The evaluation-conditions on U
  are t_i − 1 independent linear conditions on the t_i-dimensional space
  P_{<t_i} of polynomials of degree < t_i (evaluations at distinct points are
  independent functionals on P_{<t_i}), so W := {p ∈ P_{<t_i} : p|_U ≡ 0}
  has dim ≥ t_i − (t_i − 1) = 1; pick 0 ≠ p ∈ W. Since deg p < t_i and p
  already vanishes at the t_i − 1 distinct points of U, p can have NO
  further roots among S_i (a nonzero deg-<(t_i) poly has ≤ t_i − 1 roots
  over a field — already attained), so the word λ_i·(p(a))_{a∈S_i} is
  nonzero on exactly s_i − (t_i − 1) = s_i − t_i + 1 coordinates, and it is
  a member of C_i. Hence d(C_i) ≤ s_i − t_i + 1, witnessed by a word
  vanishing at exactly t_i − 1 points [C-3 wording: "vanishing at t_i − 1
  points" — not s_i − 1; at t_i = 1, U = ∅ and the word is the constant,
  vanishing nowhere: the all-ones (Λ-scaled) vector].
  (≥) Conversely if c = λ·(p(a))_{a∈S_i} ≠ 0 has wt(c) ≤ s_i − t_i, then p
  has ≥ t_i roots among the s_i distinct points of S_i, so p ≡ 0 (deg < t_i),
  i.e. c = 0 — contradiction. So every nonzero word has wt ≥ s_i − t_i + 1.
  Together: d(C_i) = s_i − t_i + 1. ∎ (K3 machine-pins both the count and the
  value at t_i = 1 and t_i = 2.)
- **L5 (separator; the induction-free core of 1).** Let S be a point support
  with |S| = w < d; suppose v ∈ V ∩ F^S, v ≠ 0; pick x ∈ S with v_x ≠ 0.
  For each axis i define T_i := {x_i} ∪ { z_i : the point z ∈ S∖{x} is
  ASSIGNED to axis i, where the assignment sends each z ≠ x to some axis
  with z_i ≠ x_i (exists: z ≠ x differs in ≥ 1 coordinate); T_i is a SET
  (duplicates collapsed) so |T_i| ≤ 1 + |S∖{x}| = w < d ≤ d_i. By L3 the
  cosets {π_i(e_a) : a ∈ T_i} are independent, and e_{x_i} ∉
  span{e_a : a ∈ T_i∖{x_i}} mod C_i (else a relation c ∈ C_i with
  supp(c) ⊆ {x_i} ∪ (T_i∖{x_i}) and c_{x_i} ≠ 0, wt ≤ w < d_i — same
  contradiction). Hence the (t_i + |T_i|)×s_i system
      ⟨φ, λ_i·(deg-<t_i monomials)⟩ = 0 (t_i rows) and φ(x_i-coord) = 1,
      φ(a-coord) = 0 for a ∈ T_i∖{x_i}
  is consistent (its non-solvability would exhibit exactly the forbidden
  relation); pick φ_i := ℓ_i ∘ π_i ∈ A_i^* a solution — a functional that
  kills C_i (hence descends to Q_i) and separates e_{x_i}'s coset from the
  T_i∖{x_i} cosets. Set y := φ_0 ⊗ φ_1 ⊗ φ_2 ∈ (F_q^N)^*. Then:
    (i) ⟨y, u⟩ = 0 for every generator A⊗C_i⊗A of V (the i-factor kills
    C_i; product functionals vanish on pure tensors with one dead factor) —
    by linearity ⟨y, u⟩ = 0 ∀ u ∈ V, i.e. y ⊥ V.
    (ii) For z ∈ S: ⟨y, e_z⟩ = Π_i φ_i(z_i-coord). If z = x all three
    factors hit their "1" coordinate ⟹ 1. If z ≠ x, its assigned axis i has
    φ_i(z_i) = 0 ⟹ 0.
  Now 0 = ⟨y, v⟩ (by (i), v ∈ V) = Σ_{z∈S} v_z ⟨y, e_z⟩ (v supported in S)
  = v_x (by (ii)) ≠ 0 — contradiction. Hence V ∩ F^S = {0}. ∎
  [Machine closure of the mechanism: K4 verifies (i)+(ii) for ALL 1,953
  w=2 supports through the fixed point at (17,(4,4,4)).]
- **L6 (sharpness; part 2 + census equivalence).** The L4 witness word of
  C_i*, embedded on any direction-i* line, lies in V with support exactly
  that line's d_i* = d points (at t = (1,1,1) it is the all-ones vector of a
  minimum-size line). So V ∩ F^{line} ≠ {0} with |line| = d, proving 2. With
  1, the union window W_w := {S : |S| ≤ w, V∩F^S ≠ {0}} is empty for w < d
  and nonempty for w ≥ d (2 gives |S| = d; superset supports stay nonempty).
  So the union-window statement "empty iff w < d" holds with first-nonzero
  weight exactly d. ∎ **No per-support claim is made at |S| ≥ d** [C-1]:
  indeed individual supports of size ≥ d may be empty (E2's lines are
  merely EXISTENCE witnesses; nothing asserts every |S| = d support
  carries a codeword — and at w ≥ d the census counts supports, not
  codewords).
- **L7 (Λ-generality).** The pointwise map M: F_q^N → F_q^N,
  M(v)_z = v_z / (λ_0(z_0)λ_1(z_1)λ_2(z_2)) is a linear bijection (nonzero
  scale per coordinate). For a direction-i line class c ⊗ (embedding), M
  maps the Λ-instance class to the Id-instance class scaled per coordinate —
  the class of (p(a)) in place of (λ_i p(a)): since {λ_i(a)p(a) : deg < t_i}
  = {p'(a) : deg < t_i} as SETS (λ_i multiplication is a bijection of the
  polynomial-evaluation space of fixed degree), M(V_Λ) = V_Id exactly.
  M preserves supports pointwise and d(C_i) is Λ-free. So Theorem T at Id
  implies T for every admissible Λ. [L7' machine-anchors the dim + the
  pointwise-division Bijection at (13,(2,2,4)) on 3 random nonzero draws.]
- **L8 (triple-sum, C-5 scope, structure only).** At t = (1,1,1), a
  direction-0 line class is a function of (a1,a2) alone, so V ⊆
  {x0(a1,a2) + x1(a0,a2) + x2(a0,a1)}; conversely each x0-term is constant
  on direction-0 lines, hence is a combination of direction-0 line classes
  (idempotent family spanning all (a1,a2)-profiles), so the set of such
  triple-sums ⊆ V. Equality of SETS. Representations are NOT unique
  (constants are representable three ways) — no uniqueness claim is made
  [C-5]. T1' machine-verifies set equality at (31,(2,3,5)) (and the frozen
  campaign verified it at (13,(2,2,4))). L8 is sanity structure; no
  theorem weight rests on it.

## (d) Line → control map

| line | claim | machine control (exact F_q) | result |
|------|-------|------------------------------|--------|
| L1 | V = Σ slices = instrument's generator stack | K1 rank(generator stack) == dim V == N − Π(s_i−t_i) at 3 instances incl. general-t | PASS |
| L2 | dim V = N − Π(s_i−t_i) | A1 (9 instances, incl. two general-t slices) + K1 | PASS (all equalities exact) |
| L3 | coset independence below d_i; dependence attainable at d_i | K2 at C_0 of (31,(2,3,5)): all 30 point-cosets nonzero (independence at |T|=1<2) AND a dependent pair exists at |T|=2=d_0 (2-point axis-0 line word ∈ C_0) | PASS |
| L4 | d(C_i) = s_i − t_i + 1; word vanishes at t_i − 1 points | K3(i): at t_i=1 all line-class weights = 4 = s_i (weights of the actual class basis, not a hardcoded count); K3(ii): exhaustive 13² polys at (s,t)=(4,2): min weight = 3 exactly, exhibited | PASS |
| L5 | separator mechanism: y ⊥ V, δ-pattern on S | K4: ALL 1,953 w=2 supports through x* at (17,(4,4,4)); y⊥V-basis and ⟨y,e_z⟩ = 1 iff z=x* for z ∈ S, every support | PASS |
| L6 | existence at d on minimum lines | E2 (16 lines at (4,4,4) + 20 at (3,4,5)); E2p (35 lines at PN4); E1corners (14/14, 41/41 — every support nonempty at d=1) | PASS |
| 1 (statement) | V∩F^S = {0} ∀|S|<d | E1union (43,744 supports at (17,(4,4,4)), d=4); E1union2 (30 singles at (2,3,5) d=2; 1,830 at (3,4,5) d=3) | PASS (zero nonzero) |
| (b)/Corollary | w=3 window at min(s)≥4 empty; (4,5,7) inside theorem | C-PN4: 457,450 supports, 0 nonzero; C-P1 (43,744, 0); C-P2 (85,400, 0) | PASS |
| instrument validity | counts match frozen instrument | A2: PN6 → 3,486 = 3,486; PP2 → 917 = 917 (byte-equal to frozen tallies under this campaign's fresh code path) | PASS |
| L7 | Λ-generality | L7' (3 random nonzero-Λ draws: dim = 13 and pointwise-division lands in V_Id) | PASS |
| L8 | triple-sum set equality | T1' at (31,(2,3,5)): dims equal, mutual membership | PASS |
| A3 | startup count asserts | every census block asserted its pre-registered total before enumerating | PASS (asserts live in blocks) |

## (e) Rule-5 watch-points that fired (all fixed pre-verdict; none touched the
theorem's mathematical content)

1. **K1 middle term (block K1, 1 fail-fix).** First draft asserted
   r_q := Π(s_i−t_i) == dim V; wrong — Π(s_i−t_i) = 27 is the QUOTIENT-image
   dim at (17,(4,4,4)), never dim V = 37. The correct machine content of L2
   is the triple rank(gen) == dim V == N − Π(s_i−t_i). Root-cause note frozen
   in `scratch_k1_rootcause.py` (run: 37 / 37 / 27 as documented).
2. **K4 off-S assert (block K4, 1 fail-fix).** First draft also asserted
   ⟨y, e_z⟩ = 0 off S; L5 makes no such claim (y's off-S values are
   irrelevant since v vanishes off S). Removed; the L5-exact checks (y⊥V;
   δ-pattern ON S) PASS on all 1,953 supports. Recorded: every K4 separator
   has nonzero off-S values (1,953/1,953 supports) — consistent with L5, and
   a useful honesty datum: the mechanism works WITHOUT any off-S vanishing.
3. **L7' wrong Λ in the division (block L7p, 1 fail-fix).** First draft
   divided by the Id-instance's Λ (all ones) — a no-op. Fixed to divide by
   the draw's own Λ. (This bug would have weakened L7' to a dim-only check;
   the fixed version re-proves the pointwise-division bijection.)
4. **K2 dead code (pre-run, caught in self-review before first launch):**
   a leftover no-op loop; removed before any run of K2.
5. **T1' generator construction (block T1p, 1 fail-fix):** first draft
   iterated TWO free axes per generator family (an IndexError); the
   correct generator of the x_{missing axis} family is ONE (d0,d1)-fiber
   profile (constant on the remaining axis). Fixed; T1' PASSes.
No other block required a fix. Total fail-fix cycles across all blocks: 4
(plus 1 pre-run cleanup) — below the pre-statement's 5-per-block escalation
threshold on every block.

## (f) Frozen label

**HUMAN-AUDITED derivation + MACHINE-VERIFIED controls.** Theorem T stands
as stated in (b), with the (4,5,7) class (and every min(s) ≥ 4 shape at
t = (1,1,1)) inside its grip: at those instances the weight-≤3 window is
provably empty and the ratio-visible window begins at support weight ≥ 4.
H-GATE's min-form is a THEOREM on the t = (1,1,1) slice; the "20 consistent
rows" are explained, not merely consistent. This file's claims are subject
to owner/Main review, which this campaign does not preempt.
