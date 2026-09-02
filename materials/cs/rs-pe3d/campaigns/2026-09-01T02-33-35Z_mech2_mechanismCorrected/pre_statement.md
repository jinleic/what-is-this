# Pre-statement — campaign `2026-09-01T02-33-35Z_mech2_mechanismCorrected` (Campaign M, mechanism)

Committed as the FIRST file of the campaign, before any control or census
compute. Everything below is pre-registered. Agent: RsPe3dMW. The campaign dir
name was minted before this file existed; nothing was computed under it before
this commit (the only prior touch is `mkdir`).

## 0. Background, inheritance, why this campaign exists

H-GATE (frozen P2, README 2026-08-30): the weight-≤3 window is non-empty iff
min_i s_i ≤ 3. The mechanism campaign `2026-08-31T08-26-54Z_hgateMechanism`
claimed a theorem. The owner review (README section "Owner correction —
2026-08-31") REJECTED five textual claims in that campaign plus two provenance
items (corrections 1-7), while stating "**The theorem survives; the wording
did not.**" The rejected wording includes correction 1: the frozen "V ∩ F^S =
{0} iff |S| < min_i d(C_i)" — read as a per-support biconditional — is FALSE;
the correct theorem is an exists-threshold statement about the union window,
with NO per-support claim at |S| ≥ d.

This campaign's job (assignment, Main dispatch 2026-09-01): derive the dual
product-code mechanism **in corrected form, first-hand** — every step either
proved with the argument written out line by line (HUMAN-AUDITED, this
campaign's own derivation, not a copy of the frozen file's prose), or
machine-checked where it reduces to exact arithmetic (MACHINE-VERIFIED). The
theorem was previously prospective for the (4,5,7) class; at t = (1,1,1) it
replaces H-GATE's 20-row evidence with a theorem on that slice.

Inheritance liabilities, accepted up front:
- C-1 (per-support biconditional false): my derivation must output the
  two-part form (union-window emptiness below d; existence at a line support
  of size d = d_min; no per-support claim at |S| ≥ d), or be marked
  OBSTRUCTION.
- C-2 (link-2 wording): a relation Σc_a π_i(e_a) = 0 yields c = Σc_a e_a ∈ C_i
  with supp(c) ⊆ T; it does NOT in general yield 1_T ∈ C_i. Independence
  follows from wt(c) ≤ |T| < d(C_i).
- C-3 (minimum GRS word): the minimum line word is the deg-(t_i−1) polynomial
  vanishing at t_i−1 of the s_i evaluation points; at t_i = 1 both readings
  coincide (the degree-0 constant polynomial, vanishing at t_i−1 = 0 points),
  which is why the t = (1,1,1) slice has the all-ones line indicator. The MDS
  fact d(C_i) = s_i − t_i + 1 is unaffected.
- C-4 (kernel dimension): dim ker(π_0 ⊗ π_1 ⊗ π_2) = N − Π_i (s_i − t_i),
  which at t = (1,1,1) equals N − Π_i (s_i − 1). The frozen text's printed
  "N − Π t_i" is wrong as printed.
- C-5 (triple-sum non-unique): V = {x0(a1,a2) + x1(a0,a2) + x2(a0,a1)} only as
  an equality of SETS, never with unique representations.
- C-6 (provenance): pre-registration claims are commit-anchored in THIS
  campaign (this file is git-committed before any compute; the frozen
  mech campaign's "first file" claim is REPORTED-grade only).
- C-7 (tuple typo): the A0 anchor instance is PP2 = (43,(2,3,7)), tally 917.
- Also inherited as established fact (machine-falsified in the frozen
  campaign, accepted here BEFORE re-derivation): "V = ker D0D1D2" is not an
  equality of spaces (falsified at (13,(2,2,4))); no difference-operator
  characterization of V is claimed in this campaign. The surviving structural
  identities are the quotient kernel identity and the triple-sum set equality.
- The K^-1623 exponent question (Gate C / AI-proof form) is OUT of scope.

## 1. Question under test, in the exact target form

**Theorem target T (union-window form), to be re-derived first-hand:**
Fix q prime (we run prime fields), s = (s0,s1,s2) with s_i ≥ 1, t = (t0,t1,t2)
with 1 ≤ t_i ≤ s_i (H1), evaluation sets S_i = pairwise-distinct points of
F_q (here: multiplicative subgroups; H2), and diagonal scalings λ_i with all
entries nonzero (H3). Let N = s0 s1 s2, C_i = { (λ_i(a) p(a))_{a ∈ S_i} :
deg p < t_i } ⊂ A_i = F_q^{s_i}, d_i = d(C_i) = s_i − t_i + 1, d = min_i d_i
(if any t_i = s_i then that d_i = 1; t_i = 0 is out of scope — the zero code,
d = ∞ by convention, never silently as s_i + 1). Let V ⊂ F_q^N be the span of
the lifted line-direction classes (L1 below), and for a point support S write
V ∩ F^S for the elements of V supported inside S. Then:
(a) for EVERY point support S with |S| < d: V ∩ F^S = {0};
(b) there EXISTS a support S with |S| = d and V ∩ F^S ≠ {0} (at
    i = argmin d_i, a direction-i line support of size d_i = d carries a
    minimum-word lift; concretely at t = (1,1,1) the all-ones vector of a
    minimum-size line);
(c) NO claim either way for an individual S with |S| ≥ d beyond the existence
    in (b).
Equivalently (census semantics): the union window {S : |S| ≤ w,
V ∩ F^S ≠ {0}} is empty iff w < d. At t = (1,1,1): d = min_i s_i, giving
H-GATE's min-form on that slice — weight-≤3 window non-empty iff
min_i s_i ≤ 3 — with the (4,5,7) class inside the theorem's grip (d = 4:
w = 3 empty; w = 4 non-empty by (b)).

Λ-generality: the quotient argument touches Λ only through d(C_i) (Λ-free by
H3), so the theorem holds for every admissible Λ; the machine checks run at
Λ = Id, with the pointwise-division bijection as the explicit reduction,
re-derived (not copy-asserted) in L7.

## 2. Derivation plan (the links L1-L8, to be audited in corrected form)

Each link is re-derived in this campaign's own words into
`theorem_corrected.md`, then anchored by the controls of §3:

- L1 (setup): V = Σ_i A_0 ⊗ … ⊗ C_i ⊗ … ⊗ A_2 under idx = (a0 s1 + a1) s2 +
  a2; the implemented `DeltaEngine.lift_rref(i)` spans A⊗C_i⊗A exactly.
- L2 (kernel identity): tensor right-exactness gives
  ker(π_0 ⊗ π_1 ⊗ π_2) = Σ_i A_0 ⊗ … ⊗ ker π_i ⊗ … ⊗ A_2 = V; rank-nullity
  gives dim ker = N − rank(π_0⊗π_1⊗π_2) and rank(π_0⊗π_1⊗π_2) = Π_i (s_i −
  t_i) [the C-4 form; NOT "N − Π t_i"]. The full-direction rank identity is
  derived via the pure-tensor span (A_0/C_0)⊗(A_1/C_1)⊗(A_2/C_2) of dimension
  Π_i(s_i − t_i).
- L3 (coset independence): for any T ⊆ S_i-index-set with |T| < d_i, the
  cosets {π_i(e_a) : a ∈ T} are F_q-linearly independent in Q_i = A_i/C_i.
  Corrected argument [C-2]: a relation Σ_{a∈T} c_a π_i(e_a) = 0 lifts to
  c = Σ c_a e_a ∈ C_i with supp(c) ⊆ T (the λ_i-diagonal preserves supports),
  so wt(c) ≤ |T| < d_i contradicts the minimum distance. NO claim that c has
  equal or all-ones entries.
- L4 (MDS distance): d(C_i) = s_i − t_i + 1. Upper bound: the deg-(t_i−1)
  polynomial vanishing at the s_i − t_i points S_i \ U for a fixed (t_i − 1)-
  set U, nonzero on U, has support s_i − t_i + 1 [C-3 wording], so d(C_i) ≤
  s_i − t_i + 1. Lower bound: a nonzero word λ·(p(a)) has support size ≥ the
  number of non-roots of p among S_i; if wt ≤ s_i − t_i then p has ≥ t_i
  roots among the distinct points S_i, so p ≡ 0 (deg p < t_i) — contradiction.
  The exact vanishing count (s_i − t_i zeros, not t_i − 1) is spelled in the
  proof file and pinned by the worked t_i = 1 and t_i = 2 examples in K3.
- L5 (separator): given a support S, |S| = w < d, and x ∈ S with v_x ≠ 0 for
  the alleged v: for each i define the SET T_i := {x_i} ∪ { z_i : the point
  z ∈ S \ {x} is assigned to i, chosen with z_i ≠ x_i }, duplicates collapsed
  (assignment exists: any z ≠ x differs in some coordinate), so |T_i| ≤ w <
  d ≤ d_i. By L3 the cosets {π_i(e_a) : a ∈ T_i} are independent, and
  π_i(e_{x_i}) ∉ span{π_i(e_a) : a ∈ T_i \ {x_i}} (same weight argument);
  hence there is a functional ℓ_i ∈ Q_i^* with ℓ_i(π_i(e_{x_i})) = 1 and
  ℓ_i(π_i(e_a)) = 0 for a ∈ T_i \ {x_i} (linear-algebra dual basis of an
  independent set); its pullback ℓ_i ∘ π_i ∈ A_i^* is what multiplies (raw
  coordinate functionals need not descend). Then y := ⊗_i (ℓ_i ∘ π_i)
  satisfies ⟨y, u⟩ = 0 for all u ∈ V (each factor kills its own C_i and the
  other factors are unrestricted, so y kills each summand A⊗C_i⊗A, hence
  their sum V), and ⟨y, e_z⟩ = Π_i (ℓ_i∘π_i)(e_{z_i}) = 1 if z = x, 0
  otherwise (for z ≠ x some coordinate has z_i ∈ T_i \ {x_i}). Hence
  0 = ⟨y, v⟩ = v_x ≠ 0 — contradiction; V ∩ F^S = {0}.
- L6 (sharpness / existence at d): the L4 minimum word of C_i, embedded as
  the line word on any direction-i line, lies in V with support exactly that
  line of size s_i... at general t the support is s_i − t_i + 1 points on
  that line; at t = (1,1,1) it is the whole line of size s_i. Combined with
  (a): the union window's first nonempty weight is exactly d = min_i(s_i −
  t_i + 1); union-window emptiness ⟺ w < d. No per-support converse at
  |S| ≥ d [C-1].
- L7 (Λ-reduction): the pointwise multiplication by the nowhere-zero function
  1/(Λ_0(a0)Λ_1(a1)Λ_2(a2)) maps V_Λ bijectively onto V_Id preserving
  supports pointwise; d(C_i) is Λ-free; so Λ = Id results hold for all
  admissible Λ.
- L8 (triple-sum, corrected scope [C-5]): at t = (1,1,1),
  V = { x0(a1,a2) + x1(a0,a2) + x2(a0,a1) } as SETS (each lifted line class
  has this form; conversely each such x0 is constant along a0-lines so it is
  a combination of direction-0 line classes, etc.). Uniqueness NOT claimed.
  Structure sanity only; carries no theorem weight.

## 3. Pre-registered controls (exact F_q; expected values fixed HERE, before compute)

Semantics notes fixed in advance: (i) "V-basis" = rref of the stacked
`DeltaEngine.lift_rref(i)` rows (the frozen instrument's convention; the PN6
3486 / PP2 917 anchors are under this convention). (ii) V ∩ F^S computed by
left-kernel of the V-basis' off-S columns (`DeltaEngine.V_inter_support`),
each returned vector re-verified to lie in V (membership vs the V-basis) and
to vanish off S before it counts as nonzero.

Anchors (must pass before census-class controls are reported):
- A1 `dim V == N − Π(s_i − t_i)` at t=(1,1,1): (13,(2,2,4))→13,
  (31,(2,3,5))→22, (17,(4,4,4))→37, (421,(4,5,7))→68, (41,(4,4,5))→44,
  (5,(1,2,2))→4, (7,(1,2,3))→6; general-t: (13,(2,2,4),(1,1,2))→14,
  (13,(2,2,4),(1,1,4))→16.
- A2 instrument cross-validation on the frozen non-empty anchors under THIS
  campaign's code path: PN6 (43,(2,6,7)) w≤3 census → 3,486 nonzero supports;
  PP2 (43,(2,3,7)) w≤3 → 917. Mismatch ⟹ HALT, hub Main, fix instrument,
  re-run A2 before anything else is reported.
- A3 count-assertions at startup of every census block (rule 17b): the block
  computes C(N,1)+C(N,2)+C(N,3) (or its w-prefix) itself and asserts the
  pre-registered total before enumerating.

Proof-machine controls (per link):
- K1 kernel identity: rank of the stacked A⊗C_i⊗A generator rows == dim V ==
  N − Π(s_i−t_i) == dim ker(π_0⊗π_1⊗π_2) (the last as N minus the exact rank
  of the quotient-map matrix) at (17,(4,4,4)) = 37, (13,(2,2,4)) = 13, and
  general-t (13,(2,2,4),(1,1,2)) = 14.
- K2 coset independence AND its boundary: at C_0 of (31,(2,3,5)) (t_0 = 1,
  d_0 = 2): (i) π_0(e_a) ≠ 0 for every point a (1-point support cannot carry
  a nonzero constant polynomial — exact membership check of each e_a against
  C_0's span); (ii) a DEPENDENT pair exists at the boundary: construct the
  C_0 minimum word (constant nonzero on 2 points, zero on the rest —
  d_0 = 2), embed as c = e_a − μ e_{a'}; then exact check: c ∈ C_0, so
  π_0(e_{a'}) = μ^{-1} π_0(e_a) and the pair is dependent. Both directions
  machine-verified: independence STRICTLY below d, dependence attainable AT d.
- K3 MDS both bounds machine-checked exactly: (i) at (17,(4,4,4)) (t_i = 1):
  every nonzero C_i word is λ·(constant), weight exactly s_i = 4 = s_i − t_i
  + 1 — min and max weight both 4 over nonzero words (constants enumerated);
  (ii) at (13,(2,2,4),(1,1,2)) axis 2 (s_2 = 4, t_2 = 2, d_2 = 3): enumerate
  ALL 13² coefficient pairs (a,b), evaluate 4-point λ-scaled words, exact min
  weight over nonzero words = 3, and exhibit a weight-3 word (deg-1 poly with
  exactly one root in S_2). Pins both the count (s_i − t_i zeros) and the
  bound value.
- K4 separator end-to-end: implement the L5 constructor (dual-basis solve for
  ℓ_i as a length-s_i vector solving the exact F_q system: t_i moment rows
  AND the point rows — concretely ℓ∘π_i must (a) kill every C_i word [t_i
  independent rows] and (b) hit the prescribed 0/1 pattern on T_i) and verify
  on ALL supports S containing the fixed point x* = (S_0[0], S_1[0], S_2[0])
  with |S| = d−1 at (17,(4,4,4)): all C(63,2) = 1,953 of them: for each,
  y = ⊗_i(ℓ_i∘π_i) must satisfy ⟨y, basis_j⟩ = 0 ∀j (V-basis, exact F_q dot
  products) AND ⟨y, e_z⟩ = 1 iff z = x* else 0 on S. Failure on any ⟹ HALT
  (the separator mechanism is load-bearing for (a)).
- T1' triple-sum set equality at a second instance: (31,(2,3,5)): dims equal
  AND every V-basis vector is an exact triple-sum (solve for the x_i) AND
  every triple-sum generator class is in V (exact membership). No uniqueness
  claim [C-5].
- L7' Λ-check: dim V at (13,(2,2,4)) under 3 fresh random nonzero-Λ draws ==
  13 each, and the pointwise-division map sends the Λ-V basis INTO V_Id
  exactly (membership), verifying the bijection machine-side.

Census-class controls (each asserts its count at startup per A3):
- E1union: V ∩ F^S = {0} for EVERY support |S| < d at (17,(4,4,4)): all
  64 + 2,016 + 41,664 = 43,744 supports of w ≤ 3 (d = 4). Expected: 43,744
  EMPTY.
- E1union' at (31,(2,3,5)) (d = 2): 30 singletons EMPTY; (61,(3,4,5))
  (d = 3): 60 + 1,770 = 1,830 supports of w ≤ 2 EMPTY.
- E1corners: (5,(1,2,2)) [N = 4, all s_i ≥ 1, d = 1있다... d = min(1,2,2) =
  1], (7,(1,2,3)) [d = 1]: theorem (a) is vacuous and (b) asserts existence;
  instrument sanity expects EVERY support non-empty: 14/14 and 41/41.
- E2 boundary existence at d: at (17,(4,4,4)): all 16 direction-0 lines
  (|line| = 4 = d) verified to carry a V-codeword (a line class supported
  exactly on the line); at (61,(3,4,5)): all 20 direction-0 lines (|line| = 3
  = d); corners covered by E1corners.
- E2' boundary existence at PN4 (421,(4,5,7)), d = 4: the all-ones
  direction-0 line word (4 points) is verified to lie in the dim-68 V-basis
  span (exact membership). Machine-side closure of theorem (b) at the W
  campaign's instance.
- C-PN4 (confirmatory): full w ≤ 3 census at (421,(4,5,7)): 140 + 9,730 +
  447,580 = 457,450 supports — expected 457,450 EMPTY.
- C-P1 (confirmatory, the frozen campaign's HELD probe): complete w ≤ 3
  census at (17,(4,4,4)) — same 43,744 enumeration as E1union.
- C-P2 (confirmatory): complete w ≤ 3 census at (41,(4,4,5)): 80 + 3,160 +
  82,160 = 85,400 supports — expected fully EMPTY.

Adjudication per control: PASS = measured value equals the pre-registered
value EXACTLY (exact integers/rationals; no tolerance). Any FAIL: halt,
root-cause, record in the report; if the failure indicates the THEOREM form
is wrong (not an instrument bug), the campaign re-scopes to the OBSTRUCTION
statement (where the proof breaks, what invariant is missing), per the
assignment's fallback.

## 4. Theorem-deliverable format and label

`theorem_corrected.md` will contain, in order: (a) hypotheses H1-H3 + scope;
(b) the correction-integrated statement (§1 form); (c) the full line-by-line
derivation (L1-L8) in this campaign's own words; (d) a table mapping each
numbered step of (c) to the §3 controls bearing on it; (e) rule-5 notes of
every watch-point that fired during the audit (unknown in advance — that is
what the table is for); (f) the final label.

Label policy (fixed now): if all §3 controls PASS, the theorem is output as
**HUMAN-AUDITED derivation (independently re-derived in corrected form;
owner skeleton credited — the quotient construction originates with the owner
per the frozen campaign; no mechanism novelty claimed here) with
MACHINE-VERIFIED controls (exact F_q arithmetic end to end)**. The theorem is
not machine-enumerated beyond the listed census slices. If any census-class
control FAILS while the proof-machine controls PASS, the label drops to
OBSTRUCTION with the precise break. If K-controls fail, no theorem is claimed.

## 5. Compute budget, resources, order

ONE background process at a time, `nice -n 10`, single slot. Pricing (from
frozen measured rates): A/K/T/L7 blocks ≈ minutes; the three big censuses ≈
583,594 supports total ≈ 6-10 min at vectorized 1,500-2,500 supports/s.
Execution order FIXED: A1 → K1 → K2 → K3 → K4 → T1' → L7' → A2 → E1union →
E1union' → E1corners → E2 → E2' → C-P1 → C-P2 → C-PN4. Every block asserts
its expectations BEFORE its loop (rule 17b) and appends a JSON line (with its
own canonical sha256) to `controls_results.jsonl`. A crash mid-census resumes
by re-running that block whole; resume = skip blocks whose line is
byte-present and checksum-valid (the delcap `validate_resume` pattern,
adapted; a mid-block crash leaves no line, so block-level re-run is clean).

## 6. What this campaign does not claim (rule 7 scope, fixed now)

Swept: exact-F_q checks, t = (1,1,1) and the two general-t slices, Λ = Id
unless stated (L7'), on EXACTLY these (q,s): (13,(2,2,4)), (31,(2,3,5)),
(17,(4,4,4)), (41,(4,4,5)), (61,(3,4,5)), (421,(4,5,7)), (43,(2,6,7)),
(43,(2,3,7)), (5,(1,2,2)), (7,(1,2,3)). Completeness limited to: all supports
w ≤ 3 at (17,(4,4,4)), (41,(4,4,5)), (421,(4,5,7)); all supports w ≤ 2 at
(61,(3,4,5)); w ≤ 1 at (31,(2,3,5)); all 14/41 supports at the corners; the
boundary line checks listed. NOT swept: any other (q,s,t,Λ); every t beyond
the named slices; nonprime fields; supports of weight ≥ 4 anywhere (Campaign
W's separate pre-statement owns the PN4 w = 4 slice — no claim about it is
made here); per-support behavior at |S| ≥ d beyond the listed line checks
[C-1]; δ/∆ and all delta-exactness sweeps (none run this campaign); ρ_inst
everywhere; Conjecture 4.2's ρ(r,η,β) content; the AI-proof exponent
question. This campaign re-DERIVES and re-ANCHORS a theorem; it does not
re-run the frozen gate-B table and makes no new ratio claims.

## 7. Adjudication of the campaign itself (pre-committed)

- Verdict THEOREM-CORRECTED-FORM: iff all §3 controls PASS and
  `theorem_corrected.md` (c) covers L1-L8 with the (d) map leaving no
  load-bearing line uncovered. Then H-GATE's min-form at t = (1,1,1) is
  recorded as THEOREM (owner-skeleton credit + correction lineage C-1..C-7
  cited), replacing the "20-row evidence" status on that slice — subject to
  owner/Main review, which this campaign does not preempt.
- Verdict OBSTRUCTION: if a control fails in a way that implicates the
  mathematical form (not the instrument) after documented root-cause. The
  obstruction statement must name the exact step that breaks, the missing
  invariant, and the smallest instance exhibiting the break (or "none found"
  if formal).
- No middle outcomes: an instrument-bug FAIL, fixed and re-run cleanly,
  counts as its control's PASS with the bug noted in the report; 5+ FAIL-fix
  cycles on one block ⟹ escalate to Main before the next attempt.
- Any H-GATE-contradicting observation (a nonempty support at w < d, or an
  all-empty E2/E1corners) is reported to Main within the same session it
  lands, BEFORE any further compute — falsifying either the theorem form or
  the instrument outweighs the schedule.
