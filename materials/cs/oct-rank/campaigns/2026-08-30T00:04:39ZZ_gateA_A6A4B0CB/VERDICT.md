# Gate A verdict — re-verification of the published rank-25 Krawczyk certificate

**VERDICT: PASS (MACHINE-VERIFIED, as a re-verification).** Timestamp of
freeze: 2026-08-30T00:0xZ. Campaign directory immutable after this write.

## Inputs
* Paper: arXiv:2608.16649v1 (Hardik Jain, 2026-08-17), read first-hand in full.
* Upstream code: github.com/hxrdxkxvxd/octonion-rank @ 816a01e16acb5849ba5ad9858e80cbb5a39b967d
  (2026-08-19), snapshot at scratch/upstream_ref with SHA256 inventory below.

## What was re-verified — exact rational inequalities (gate_a_verify.py,
## python-flint fmpq + arb, all arithmetic EXACT, final comparisons also
## re-evaluated with OUTWARD-ROUNDED arb balls):
1. A1 tensor: my independent Cayley–Dickson construction of T_O equals the
   upstream table entrywise; unit + norm multiplicativity checked.
2. A2 preconditioner consistency: ||Y·J_S − I||_inf (exact rational)
   = 4352331780372686777606956781732931893/1461501637330902918203684832716283019655932542976
   = 2.978×10⁻¹² < 10⁻⁶.
3. A3 contraction: K = 934065388564663273070817924426182673785089593718087090786657
   /2923003274661805836407369665432566039311865085952000000000000
   = 0.3195567369566958… < 1 (exact; matches paper's K ≈ 0.32 / repo 0.319557).
4. A4 Banach margin: (1−K)·ρ − max_e |(YF(x₀))_e| = +6.804×10⁻⁷ > 0 (exact;
   off_max = 1.314×10⁻¹⁴, matching paper's ≈1.3×10⁻¹⁴).
5. A5 inclusion radii: every one of the 512 per-equation margins
   ρ − (off_e + ρ·rs0_e + ρ²·rs1_e + ρ³·rs2_e) is a strictly positive exact
   rational; worst = same +6.804×10⁻⁷ (matches paper's 6.8×10⁻⁷).
6. Outward-rounded re-evaluation of 3–5 with arb balls (all strict):
   K_arb = [0.319556736956696 ± 2.14e-16] upper < 1;
   banach margin ball lb ≥ 0; worst margin ball [6.80443262388140e-7 ± 9.89e-22] lb > 0.

## Independent Krawczyk radius re-derivation (gate_a_radius.py) [DERIVED]:
Bisection (60 iters, exact + outward) gives: hypotheses hold for ALL
ρ < ρ* > 3.1293×10⁻⁶ (ρ* ∈ (3.129332e-06 bisect interval), exact dyadic
lower bound 8400235957/2684354560000000), K(ρ*) → 0.9999997545967169.
The paper's ρ = 10⁻⁶ therefore has ≈3.1× slack inside the certified radius.
All three flags fail together just above ρ* — the radius is
contraction-limited (hypothesis (i)), not off-vector-limited.

## Lean 4 replay (lean_replay/):
* Toolchain leanprover/lean4:v4.29.1 (repo-pinned), mathlib v4.29.1
  (rev 5e932f97dd25535344f80f9dd8da3aab83df0fe6).
* `lake exe cache get` (8232 files) + `lake build`: **completed
  successfully, 2398 jobs.**
* Axiom audit (`#print axioms`) on oct_rank_ge_eighteen,
  eighteen_le_slicesRank_T_O, tau_rank_ge_seven, seven_le_slicesRank_tau,
  Krawczyk.krawczyk_exists, Krawczyk.krawczyk_exists_of_bounds,
  Krawczyk.slicesRank_le_of_residual_zero(+_tau), isSlicesOfRank_transpose,
  Lmat8_conj_mul_self, Lmat8_polar, inn8_sq_lt_of_linearIndependent,
  pencil_rank_ge, slicesRank_L_family_ge, substitution_pivot0,
  substitution_pivot_last: ALL depend exactly on
  [propext, Classical.choice, Quot.sound]. No sorry, no additional axioms.
* Kernel-accepted statements recorded via #check in lean_replay.log; they
  match the paper's prose claims (18 ≤ rank(T_O); 7 ≤ rank(τ); Proposition-6
  Krawczyk existence in both Banach and fderiv-bounded forms; bridge lemmas
  residual-zero ⇒ slice decomposition for r = 25 and r = 7). NO statement
  found materially weaker than the paper's prose.
* Note: the certificate inequalities themselves are NOT inside the kernel —
  evaluated by scripts (exact rational), exactly as the paper states
  ("the boundary between the two is stated where it arises").

## Provenance tags
* Numeric certificate values: [REPRODUCED] (independent implementation).
* Radius ρ* ≈ 3.1293×10⁻⁶: [DERIVED] (my bisection; not in the paper).
* Lean theorems: [REPRODUCED] replay of upstream artifact.
* Fiduccia–Zalcstein ≥ 15 / Cariow–Cariowa ≤ 30: [REPORTED] (not fetched);
  no gate depends on them.

## Evidence labels
* 18 ≤ R_R(T_O) ≤ 25 re-verification: MACHINE-VERIFIED (this artifact:
  gate_a_verify.py PASS + gate_a_radius.py + lean_replay/lean_replay.log).
* The underlying mathematics: CITED-DEPENDENCY (arXiv:2608.16649 + its Lean
  artifact). Nothing here independently re-proves the peeling/pencil math;
  the algebra was AUDITED by reading the Lean sources line by line
  (Pencil.lean key_bound; Substitution/Oct18Peel; identity theorems in
  Oct18Identities.lean) — kernel-side statements as above.
