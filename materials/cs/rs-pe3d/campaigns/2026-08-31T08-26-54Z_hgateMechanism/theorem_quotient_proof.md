# The quotient proof of H-GATE at t=(1,1,1) — audited statement + machine checks

Campaign `2026-08-31T08-26-54Z_hgateMechanism`. This file records (1) the
theorem with its exact hypothesis class, every link audited from scratch;
(2) the machine controls that anchor each algebraic claim; (3) the honest
scope of what a MACHINE-VERIFIED label covers here.

## 1. Setup and hypotheses (stated precisely)

- Indices/conventions: `idx = (a0·s1 + a1)·s2 + a2`, `a_i` indexes S_i ⊂ F_q^×
  in ascending order (src/rs.py `Inst`).
- A_i = F_q^{s_i} with coordinate-per-point evaluation at S_i.
- C_i ⊂ A_i is the Λ_i-diagonal GRS code: C_i = { (λ_i(a)·p(a))_{a∈S_i} :
  p ∈ F_q[x], deg p < t_i }, with hypotheses:
  - **(H1) 1 ≤ t_i ≤ s_i** — needed for d(C_i) = s_i − t_i + 1 (MDS distance)
    and for the minimum-word upper bound. If the paper ever allows t_i = 0,
    that is the zero code C_i = {0}, d(C_i) = ∞, handled separately (never
    silently as s_i + 1). All frozen rows have positive t_i.
  - **(H2) S_i = distinct evaluation points** (they are, being subgroup
    elements; pairwise coprimality of the s_i is IRRELEVANT to this theorem);
  - **(H3) all λ_i(a) ≠ 0** (enforced by `Inst`'s assert) — the diagonal
    scaling is weight-preserving, so d(C_i) is Λ-free.
- π_i : A_i → Q_i := A_i / C_i. V := span of the lifted line bases =
  Σ_i A_0⊗…⊗C_i(correct factor)…⊗A_2 (the implemented `DeltaEngine` lift
  stack; `lift_rref(i)` = rref of all i-lines' bases).

## 2. The five links, each audited (audit done independently before compute)

1. **Kernel identity.** Tensor is right-exact: ker(π_0⊗π_1⊗π_2) =
   Σ_i A_0⊗…⊗ker π_i⊗…⊗A_2 = V. Exactness of this needs no coprimality and
   no Λ condition (H3 only ensures C_i has the stated GRS form).
2. **Coset independence.** {π_i(e_a) : a ∈ T_i} is F_q-independent whenever
   |T_i| < d(C_i): a relation Σ c_a π_i(e_a) = 0 ⟺ 1_T ∈ C_i (mod scaling by
   λ_i — nonzero on every coordinate, so support unchanged), and
   wt(1_T) = |T_i| < d(C_i) contradicts the minimum distance. [MDS; the
   Λ-diagonal is support-preserving.]
3. **Point isolation (the owner's assignment step, corrected).** Given a
   support S of w points, w < d := min_i d(C_i), and a target x ∈ S: for each
   i, define the SET T_i := {x_i} ∪ { z_i : the point z ≠ x is assigned to
   coordinate i, chosen with z_i ≠ x_i }. Assignment always exists (any z ≠ x
   differs in some coordinate). Duplicates in T_i are collapsed (it is a SET);
   |T_i| ≤ w < d. So {π_i(e_a) : a ∈ T_i} is independent by link 2.
4. **Separator construction.** Since e_{x_i} ∉ span{ e_a : a ∈ T_i \ {x_i} }
   mod C_i (same weight argument), there is a functional
   ℓ_i ∈ Q_i^* with ℓ_i(π_i(e_{x_i})) = 1 and ℓ_i(π_i(e_a)) = 0 for all
   a ∈ T_i \ {x_i}. Concretely ℓ_i∘π_i is an A_i^* vector solving an exact
   (t_i + |T_i|)×s_i F_q system (t_i "polynomial-moment" rows + |T_i| point
   rows; solvable because the T_i ¬ {x_i} cosets are independent). The
   coordinate functionals do NOT descend to Q_i by themselves — ℓ_i is
   constructed modulo C_i; the pullback ℓ_i∘π_i is what multiplies.
   Then y := ⊗_i (ℓ_i∘π_i) ∈ A_0^*⊗A_1^*⊗A_2^* satisfies: each factor kills
   its own C_i, so ⟨y, u⟩ = 0 for every u ∈ Σ_i A⊗C_i(correct slot)⊗A = V,
   i.e. y ⊥ V. And for z ∈ S: ⟨y, e_z⟩ = Π_i [ℓ_i∘π_i](e_{z_i}) equals 1 if
   z = x (each ℓ_i∘π_i pairs 1 with e_{x_i}) and 0 otherwise (z ≠ x is
   assigned to some coordinate i with z_i ∈ T_i\{x_i}, where ℓ_i∘π_i
   vanishes). Now suppose v ∈ V ∩ F^S, v ≠ 0, and pick x with v_x ≠ 0.
   Then 0 = ⟨y, v⟩ = Σ_{z∈S} v_z⟨y, e_z⟩ = v_x ≠ 0 — contradiction.
   Hence V ∩ F^S = {0}. ∎
5. **Sharpness (both directions at once).** A minimum-distance word of C_i
   (for the i-line: the polynomial of degree t_i−1 vanishing at s_i−1 of the
   s_i points and nonzero at the last — exists since t_i ≤ s_i; at t_i = 1 it
   is literally the all-ones line indicator) lifts to a V-codeword supported
   on that i-line. So V∩F^S ≠ {0} for S = that line, |S| = s_i = d at t=1.
   Therefore: **window emptiness ⟺ |S| < min_i d(C_i)** — the biconditional,
   with d(V ∩ F^{·} first-nonzero-support-weight) = min_i (s_i − t_i + 1)
   exactly.

**At t = (1,1,1):** d(C_i) = s_i, so the weight-≤3 window is empty ⟺
min_i s_i ≥ 4 — H-GATE's min-form, now a THEOREM on this slice. Hypothesis
class, exactly: prime power q (we run prime fields), pairwise-coprime-or-not
irrelevant, any s = (s_0,s_1,s_2) with s_i ≥ 1 (s_i = 1 fine: C_i = F_q^1,
d = 1, and V = A there when s_i = t_i), any Λ_i with all entries nonzero,
any t with 1 ≤ t_i ≤ s_i (general-t predicate: empty ⟺ |S| < min_i(s_i−t_i+1)).

**Λ-uniformity made explicit:** the quotient argument never uses Λ except
through d(C_i), which is Λ-free. Alternative Λ=Id route (this campaign's
original derivation): pointwise division by Λ_0(a_0)Λ_1(a_1)Λ_2(a_2) ≠ 0 maps
V_Λ bijectively to V_Id preserving supports. Machine-checked: dim V at
(13,(2,2,4)) is 13 for three random nonzero Λ draws and for Λ=Id.

## 3. Machine checks actually run (all exact F_q, no floats)

| id | check | instance(s) | result |
|----|-------|-------------|--------|
| A1a | dim V == N − Π(s_i−1) at t=(1,1,1) | (13,(2,2,4)) | 13 = 13 PASS |
| A1b | dim V == N − Π(s_i−1) | (31,(2,3,5)) | 22 = 22 PASS |
| A1c | dim V == N − (s_0−1)^3 | (17,(4,4,4)) | 37 = 37 PASS |
| K1 | kernel identity: rank(Σ_i A⊗C_i⊗A) == dim V == dim ker(π_0⊗π_1⊗π_2) with dim ker = N − Π t_i | (17,(4,4,4)) | 37 = 37 PASS |
| E1 | V ∩ F^S = 0 for ALL supports w < d (complete enumeration below the theorem's threshold) | (13,(2,2,4)) d=2 w<2: 16 checks; (31,(2,3,5)) d=2 w<2: 30; (61,(3,4,5)) d=3 w≤2: 1,830 | all EMPTY PASS |
| E2 | boundary: every i-line (w = d) carries a codeword | (61,(3,4,5)) dir-0 lines; (41,(4,4,5)) all axes | all nonempty PASS |
| E3 | w = 3 complete census on a d = 4 probe shape | (41,(4,4,5)) all C(80,3) = 82,160 supports | 0 nonzero PASS (52.9 s) |
| E4 | w ≤ 2 on second probe shape | (61,(4,5,5)) 5,050 supports | 0 nonzero PASS |
| C1 | s_i = 1 corner: EVERY support nonempty (⇐ arm incl. corner) | (5,(1,2,2)) N=4: 14/14 | PASS |
| C2 | s_i = 1 corner | (7,(1,2,3)) N=6: 41/41 | PASS |
| A0 | instrument re-validation vs frozen counts | PN6 (43,(2,6,7)) → 3486; PP2 (43,(2,3,7)) → 917 | 3486=3486, 917=917 PASS |
| L1 | Λ-uniformity | (13,(2,2,4)), 3 random nonzero Λ | dim V = 13 each PASS |

Failed intermediates are NOT hidden (rule 5).
Control bugs, root-caused and rerun: (i) a `kernel_mod` arity mistake
(`field.kernel_mod(rows, q, ncol)`; my first control passed the arguments in
the wrong order), (ii) a wrong self-assertion in the first K1 draft
(I asserted `rank(Σ_i A⊗C_i⊗A) == N − dim V`; the correct assertion is
`rank == dim V` because V IS that kernel space — fixed, K1 passes).
Instrument-truth discovery: the difference-operator identity "V = ker
D0D1D2" is WRONG as an equality of spaces. The machine check at
(13,(2,2,4)) found codewords of V not killed by the D_i (28 of the stacked
difference rows pair nontrivially with V-basis vectors; also
dim ker(D0D1D2) = 1 there while dim V = 13 — at (2,2,4) the three summands
overlap so heavily across degenerate directions that the naive
difference-operator identity overcounts). The pre-statement §1 sentence
"Equivalently V = ker D0D1D2" is RETRACTED as an equality-of-spaces claim
(recorded here and in the report). What SURVIVES, verified on every instance
checked (A1a/A1b/A1c): the DIMENSION identity
dim V = N − Π_i (s_i − t_i) at t=(1,1,1), pre-registered as checkable in the
pre-statement and anchored by A1.

**Label discipline.** E1/E2/E3/E4/C1/C2/A0: MACHINE-VERIFIED (exact F_q
arithmetic end to end; semantics audited). A1/L1/K1: MACHINE-VERIFIED. The
THEOREM itself: proven by the audited argument above (steps 1-5, each link
mathematically derived, with links 2-4 numerically anchored by E1/E2 and K1);
the theorem is DERIVED + machine-anchored, not machine-enumerated (except on
the enumerated slices listed). Not swept, plainly: full P1/P3 censuses beyond
the partial controls above; the general-t boundary at t_i = s_i exactly
(corner where C_i = A_i, d = 1, handled by theorem but not enumerated);
t_i = 0 (out of scope per Main: d = ∞ convention, not applied anywhere).

## 4. What remains for the heavy slot (pre-registered, held)

P1 (17,(4,4,4)) full w≤3 census (43,744 supports), P2 (41,(4,4,5)) already
w=3-complete (82,160, empty — E3), P3 (61,(4,5,5)) w=3 remainder
(161,700 + 4,950 + 100). These are now CONFIRMATORY for the theorem
(a theorem no longer needs them, but they were pre-registered before the
proof existed and Main's instruction is to keep them held and preserve
order). They are also the census-level check that the implemented V-basis
matches the theorem's hypothesis class (C_i exactly the GRS codes).

## 5. Consequence for the H-GATE ledger

- H-GATE at t=(1,1,1) is a THEOREM with hypothesis class {all q, all
  s_i ≥ 1, all nonzero Λ, d-hypotheses (H1)(H2)(H3)}; the empirical "20
  consistent rows" are now explained, not merely consistent. The
  paper-conjecture content of Conjecture 4.2 (rho(r, η, β) form) is NOT
  touched by this — the gate program's effective domain becomes exactly
  {instances with an order ≤ 3 at t=(1,1,1)}, and **ρ^window is undefined
  (empty window) at min(s) ≥ 4** there; measuring rho at such instances
  requires the weight-≥4 window, still OPEN everywhere.
- The general-t form: window emptiness threshold = min_i (s_i − t_i + 1).
  The three frozen q=13 t-profiles (1/2, 3/8, 1/4 at wt 2/3/1... README rows)
  are consistent with this reading at d = min(2,2,4−t_2+1): t=(1,1,4) gives
  d = 2 (still wt-1 witness... [[INFERENCE — the frozen rows' wt/delta
  pattern is quoted from README lines 172-174 and was NOT re-derived here;
  treat the t>1 slice as COROLLARY-CONSISTENT, not machine-verified in this
  campaign]]).
