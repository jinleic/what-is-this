# Route AF verdict — 2026-08-31

Campaign: `campaigns/2026-08-31T08:02:18Z_routeAF/` (pre-statement committed
first, before computation).  Rules 1–17c of the repo README and the
oct-rank README gate ladder were read before any code.

## Headline

**ROUTE F: OPEN at 13-or-14 — neither certificate obtained; verdict remains
`[13,14]` for the specific triple (1,i,j).  ROUTE A: CLOSED NEGATIVE —
the proposed rank-16 floor is refuted; the exact Strassen-form import gives
only 12.  RESOURCE BREACH: two numerical searches ran during KgBandClose's
exclusive heavy slot; all timing interpretation is forbidden.**

## Ordered exact statements

Conventions are the pre-statement's: Cayley--Dickson octonions in basis
(1,i,j,k,l,il,jl,kl), `L_x[c,b]=(x*e_b)_c`.  Let tau=(L_1,L_i,L_j) on the
quaternion subalgebra (rank exactly 7), and let `s = e_0 ⊗ (f_1⊗f_1 + f_2⊗f_2)`.

1. The conjugated Route-F target for the triple (1,i,j) is the
   **shared-first-factor block duplication**
   `T_F = tau ⊠ s`, equivalent to `sym(P diag(d) Q⁻¹)` for
   `P=[I₄,0;0,I₄]` (8x8) and even multiset `d ∈ ℝ¹⁰`;
   the slices are the three blockdiag(tau_slice, tau_slice) 8x8 matrices.
   [MACHINE-VERIFIED — exact `fmpq` identities, `upper14_replay.out`,
   `routeF_block_replay.out`.]

2. Upper 14 on the specific triple (1,i,j): certified 14-term CP witness
   obtained by concatenating two frozen rank-7 tau decompositions blockwise.
   [MACHINE-VERIFIED modulo frozen tau r7 certificate exact facts; the
   certificate decimals are approximate, the tau-facts layer is exact —
   `anchor_replay.out`, inherited `tau_r7` artifacts.]

3. Lower 13 on the specific triple: inherited chain
   `1 + pencil_floor(12)` reproduced exactly after the unimodular shift
   identity (det = 1, rank-preserving).  [MACHINE-VERIFIED —
   `lower13_chain_replay.out` replaying S3 chain artifacts.]

4. The two facts together give triple rank ∈ {13,14}; deciding needs an
   exact rank-13 witness (down) or rank ≥ 14 proof (up).  Not obtained:
   the numerical candidates best relative residuals are 1.33e-3 (direct
   CP, 5 seeds, 5000 evals each, all failing the 1e-11 trigger) and 1.00e-4
   (commuting-extension model, 3 seeds); least-squares nonconvergence is
   **not** an infeasibility proof.  [COMPUTATIONAL-EVIDENCE only.]

## Route A closure (negative)

Human algebraic proof (frozen in `commutator_derivation.txt`), anchored by
exact coefficient checks (`route_af_commutator.py.asrun` →
`commutator_proof.out`): for every linearly independent real octonion
triple (u,v,w), with U=N(u)>0,

`[L_ūL_v, L_ūL_w]² = −4·U·det Gram(u,v,w)·I`, det Gram > 0,

so the commutator is invertible and its rank is **exactly 8** on every
independent triple; on dependent triples the square vanishes and rank ≤ 4,
so dependence cannot improve any universal floor.  This is a proof plus
exact coefficient anchors — explicitly **not wholly machine-proved**.

Imported theorem (CITED-DEPENDENCY; Strassen-primary provenance
UNESTABLISHED): Landsberg, *Geometry and the Complexity of Matrix
Multiplication*, §6.1, Theorem 6.1.1 — exact transcription in
`source_audit.txt` — states, for a full-rank slice,
`Rank[T_{α,α₁},T_{α,α₂}] ≤ 2(r−b)`.  Order/conjugacy is pinned exactly:
`L_vL_u^{-1}=L_u X L_u^{-1}`, so the survey-order commutator and the
X,Y-order commutator have equal rank.  Substituting b=8, commutator
rank 8:

`8 ≤ 2(r−8)`  ⟹  `r ≥ 12` (**not 16**).

The retracted reading `r ≥ b + rank(comm)` (which would give 16) is refuted
by the mandatory tau control: it would assert `4+4=8` for a tau whose exact
rank is 7; the correct rearrangement gives `4+4/2=6 ≤ 7`.  Koiran
(arXiv:2006.02374, Theorem 1 / Lemma 6) confirms the ½-factor form.
The 3n/2-vertical bound: 12 ≤ 3·8/2 = 12 — exactly attained but NOT tight;
floor 13 from the chain stands, so Route A offers no path to 14.

**Retraction:** the earlier phrase "trace formula: rank ≥ a+rank(commutator)"
is retracted as an unfaithful paraphrase of the source; the correct form has
the ½ factor.  This retraction is recorded here and in the target README.

## Route F theorem audit (negative)

Frozen in `routeF_theorem_audit.txt`: no source-backed additivity /
multiplicativity theorem covers `tau ⊠ s`.  Checked and rejected:
Christandl–Jensen–Zuiddam Proposition 22 (multiplicativity requires the
nonmatrix factor to be a matrix tensor, i.e. a = 2 over ℂ; here a=3, and
rank is wanted over ℝ — the authors call their case "essentially minimal").
Strassen's direct-sum additivity applies to full direct sums
(⊥-doubling of all three factors), not to the shared-first-factor block
duplication.  Without an exact 13-witness the model miss stays
COMPUTATIONAL-EVIDENCE.

## Resource-contract breach

Frozen in `resource_contract_breach.txt`.  Two commands
(`route_af_search.py`, then `route_af_extension.py`) started at
08:09:46Z and 08:13:26Z inside KgBandClose's exclusive heavy slot.
Harness-recorded wall times: 105.30 s and 32.67 s.  CPU times were NOT
logged.  **All timing/performance/cost interpretation of these runs is
forbidden.**  The residual rows are usable only as
COMPUTATIONAL-EVIDENCE about the mathematical miss.  No heavy script was
run after this finding; the finding was self-reported to Main at 09:07Z
and is recorded in the target README.

## Rule-7 scope sentence

This campaign swept exactly the fixed conjugated Route-F tensor
`tau ⊠ s` (slices = blockdiag copies of the three tau slices on the fixed
Cayley–Dickson basis) at target rank 13, with 5 random+ALS full-CP seeds
(5000 evals each) and 3 real-diagonal commuting-extension seeds
(5000 evals each, n=8, r=13, k=5), plus the two frozen inherited replay
chains and the exact commutator coefficient program (octonion: 8+64+1+28
checks, 56 independent rank controls; quaternion control analogous); it did
NOT sweep other triples of the 8×8×8 octonion multiplication tensor, other
ranks, complex decompositions, border rank, non-CP constructions, the full
direct-sum `tau ⊥ tau` (a different tensor), finite-field/rational
enumeration lower bounds, and it found no rank-14 obstruction among
Strassen-form theorems available in this campaign's sources — "no improving
direction found among those tried" is not "none exists".

## Premises of this prompt disproved

1. "Route A (Strassen commutator) gives 16" (S3 row) — refuted: the
   rearranged Strassen form gives 12; the naked twin overshoots the tau
   control (8 vs true 7) and misses the ½ factor.
2. "The (1,i,j) verdict may be decided by the least-squares slot I was
   promised" — disproved: numerical refinement cannot certify either side;
   only exact/certified artifacts can, per the pre-registered adjudication.

## Next campaign (named)

**`af-triple13-existence`:** decide 13-vs-14 for the specific triple by
attacking the exact 10-variable symmetric-real realification of the rank-13
CP system with exact-rational Krawczyk/interval Newton on a least-squares
critical point (positivedefiniteness of the 20x20 symmetrized Jacobian is the
hypothesis to check — ≤ 30 min pinned on the free slot), with fallback
Gröbner elimination on the 1+20-variable exact system (risky, un-budgeted;
estimate 2–6 h).  Downward success ⇒ rank exactly 13 for (1,i,j); a genuine
first.

## Freeze inventory

`pre_statement.md` (first file), `source_audit.txt`,
`commutator_derivation.txt`, `routeF_theorem_audit.txt`,
`resource_contract_breach.txt`, `.asrun` byte-copies of all executed
sources (`route_af_search.py.asrun`, `route_af_extension.py.asrun`,
`route_af_commutator.py.asrun`, inherited-replay sources
`s3_anchor.py.asrun`, `s3_routeC_substitution.py.asrun`,
`s3_routeF2_blocks.py.asrun`, `s3_routeF3_global.py.asrun`), `.out`
artifacts, both candidate `.npz` files, `checksums.sha256`, `manifest.json`.
Unexecuted drafts (`route_af_refine.py`, `route_af_structured.py`) were
deleted rather than mislabeled as-run; their edit history is in the session
transcript.
