# Pre-statement — Gate A (waterfall census), mceliece target

**Committed 2026-08-29T17:00Z, BEFORE any computation run.** This file fixes the
quantity computed, the conventions, the budgets, the certificate format, and the
pass/fail threshold. Thresholds may not be adjusted after seeing results; any
re-scoping is recorded below with the original left visible.

## 1. Primary-source basis (all read first-hand this session)

- Apon, ePrint **2026/1810** — Theorem 1 (lower bound, eqs. (10)/(11) hypotheses),
  Definition 4 (auxiliary map T(A) = Q_A + ρ_A·Z·F), Lemmas 2/3 (binary
  differential identity ΠF′ + Π′F = κG²F⁽²⁾), eq. (16) (F′ ≡ R_G·F mod G²),
  Lemma 6/7 (flag coset law [T(A)]_{a,r} = A(a)·[T(1)]_{a,r}), Remark 8
  (rank exactly 1 iff F(a), F′(a) E-independent; T(1)(a) = G(a)⁻²F′(a) + σ(a)F(a)),
  Lemma 9 (injectivity via Δ_{p,q} = f_p f_q′ − f_q f_p′, eq. (25)),
  eqs. (3)/(6)/(8) (Sol(S) definition with true flags), §2 protocol ("perfect
  oracle"), Table 1 (12 small instances, c_need = 2t+3 in every row).
- Vedenev, ePrint **2026/1747 v2** — eq. (13) flag conditions, Remark 6
  (per-position count ≥ (k choose 2) + r_i − k + 1), Conjecture 1, §8.1
  c* = 2δ+3 for binary Goppa, rejection criterion (24)
  dim W(s,a) < 2δ+4−c, Conjecture 5, F4 experiments (dim W = 60, 28 at c = 4).
- Saarinen, ePrint **2026/1786** — Eqs. (29)/(30) reading of Apon (his
  cross-check), interpolation-lattice eq. (27), strict-profile count B = (k_ℓ
  choose 2) eq. (28)/(31)/(32), Proposition 3 (unindexed complete-flag count is
  exact), bundle criterion (36). Relation-generation layer only; not used for
  Step-3 assembly.
- Ghoshal–Ishai–Jain–Sun (GIJS), ePrint **2026/1630** — Fact 2.3 (multipliers
  λ_i = Γ(α_i)/Π′(α_i), λ′_i = Γ(α_i)²/Π′(α_i) for square-free case; binary
  Goppa with irreducible G is square-free), Definition 2.7 (Hasse derivatives),
  Thm 4.1 (all-kernel distinguisher), §7/Thm 7.2 + App. C (conditional key
  recovery; nullity excess ms + m − 1 prediction).

## 2. The exact quantity computed (Gate A)

For each of the nine ladder instances (m, n, t; k = n − mt, D = n − 2t − 1):

Build the **true** instance: Γ (called G) a seeded random monic irreducible
degree-t polynomial over E = F_{2^m}; support L = (a_1..a_n) distinct elements of
E (seeded prefix of one seeded shuffle of a fixed enumeration); Π = ∏(Z − a_i);
multipliers λ_i = G(a_i)²/Π′(a_i) (Apon eq. (1), κ = 1; GIJS λ′ — square-free
case applies since irreducible ⇒ square-free); the public code generator rows
and the hidden vector polynomial F = (f_1..f_k) interpolating them:
**f_j = the unique polynomial of degree ≤ D with λ_i·f_j(a_i) = (G_pub)_{j,i} ∈ F₂^k**,
realized concretely by choosing k linearly independent rows and requiring
y_i = λ_i F(a_i) ∈ F_2^k (i.e. F(a_i) = λ_i⁻¹ y_i, which must itself lie in the
subfield — checked exactly).

Measured quantities, per instance, per held count c ∈ {2t+1, 2t+2, 2t+3, 2t+4},
held set S_c = prefix of the (single, seeded, shuffled) support ordering:

1. **Family nullity** Nfam(c) = dim_E { A ∈ E[Z]≤2t+2 : T(A) satisfies
   A(a)=0 ∀ a ∈ S_c }... precisely: by Apon Lemma 7 the flag conditions on the
   family T(A) hold iff A vanishes on S_c (one scalar condition per point, rank
   ≤ 1 per point). We measure this directly as
   Nfam(c) = dim_E { A ∈ E[Z]_{≤2t+2} : A(a) = 0 ∀a ∈ S_c } = 2t+3 − c for
   c ≤ 2t+3, computed by exact echelon on the (c × (2t+3)) evaluation matrix
   **and verified against the nontrivial sufficient direction**: for the
   vanishing polynomials A ∈ P_S·E[Z]_{≤2t+2−c}, T(A) is machine-checked to
   satisfy the actual flag conditions (∂[j]T(A))(a) ∈ V_a^[j]) up to depth R at
   every held point. This certifies EF + T(A_S) ⊆ Sol(S_c) not just at c = 2t+3
   but at the exact rank level asserted by Lemma 7.
2. **Per-point rank certificates.** For each held point a ∈ S_c and each
   derivative order j < R (= 4, see §3): form the flag condition
   (∂[j]H)(a) ∈ V_a^[j] restricted to the (2t+3)-dimensional A-family, i.e. the
   E-linear map A ↦ (∂[j]T(A))(a) mod V_a^[j]; compute its rank r_{a,j} by
   exact echelon. **Apon's Lemma 7 predicts r_{a,j} ≤ 1 for every a, j.**
   Additionally the **order-0 nondegeneracy witness**: verify
   T(1)(a) ∉ E·F(a) at every held point (Remark 8), so the rank-1 is attained,
   i.e. **exactly** 1, not 0.
3. **New-equation increments.** Δ(c) = Nfam(c−1) − Nfam(c) (=1 for all
   c ≤ 2t+3 under the theorem). A refutation event is any measured deviation
   (see §5).

Guards machine-checked at build (all exact, per Apon §2 / README): (α)
differential identity ΠF′ + Π′F = κG²F⁽²⁾ coefficient-wise; (β) Δ_{p,q} ≠ 0 for
an exhibited pair p < q; (γ) k = n − mt; (δ) λ_i F(a_i) ∈ F_2^k exactly for all
i; (ε) gcd of F's coordinates = 1 and max coordinate degree exactly D.

## 3. Fixed conventions

- **Convention (declared):** Apon's Step-3 system with full oracle — Sol(S) per
  eqs. (6)/(8) with the TRUE support, TRUE derivative flags V_a^[j] =
  span_E{F(a),…,(∂[j]F)(a)}, TRUE Frobenius alignment. This is Apon's own §2
  experiment ("grants the attack a perfect oracle … and isolates only the rank
  of Step 3"). Rationale: both disputants' claims are statements about this
  system's rank behavior ( Vedenev v2 Theorem 1 quote and criterion (24) are
  about exactly this W(s,a)). A refutation found under the perfect-oracle
  convention cannot be dismissed as a convention mismatch.
- **Flag depth bound R = 4** for the per-point rank certificates: orders
  j ∈ {0,1,2,3}. Justification: (i) j = 0 already exercises the identity
  [T(A)]_{a,r} = A(a)[T(1)]_{a,r} at its sharpest (Remark 8's rank-exactly-1);
  (ii) j = 1,2,3 exercise Lemma 6's mechanism (orders j − q + 1 ≤ j of F at a,
  parity kills even scalars) at increasing depth; (iii) depth R is a free fixed
  choice made here, before results, and the theorem's claim is depth-uniform —
  deeper-r generality at ladder scale is recorded as CITED-DEPENDENCY (Apon
  Lemma 7) and is machine-verified at small t in Gate C. R = 4 keeps each
  per-point certificate ≈ 4 rank computations of width 2t+3.
- **Held-set convention:** one seeded shuffle of a fixed enumeration of E per
  instance (seed recorded); S_c = first c distinct support elements in that
  order; c runs 2t+1 … 2t+4 (window from the README). No adaptivity: S_c is a
  prefix chain.
- **T(1) normalization:** T(1) = Q_1 + ρ_1·Z·F with G = the instance Goppa
  polynomial, computed by Definition 4 exactly (B_1 ≡ R_G mod G², deg B_1 < 2t).
- **Precision:** exact finite-field arithmetic throughout (integer Zech-log
  tables over uint16 representation, validated against python-flint fq_default
  on random pairs; one echelon cross-checked against flint if the API allows).
  No floating point anywhere. Rounding mode: N/A (exact).

## 4. Budget & certificate format

- One process at a time, `nice -n 10`, all BLAS/OMP threads pinned to 1.
- Per instance: 1 build + guards, 4 echelons of size ≤ (2t+4)×(2t+3) (family
  nullity), per held point 4 echelons of size (2t+3)×k... restricted rank maps
  (∂[j]T(A))(a) mod V_a^[j] are k-vectors linear in A ∈ E^{2t+3}; the quotient
  map has matrix (k − dim V_a^[j]) rows × (2t+3) columns, ≤ 4 echelons of size
  ≤ (c)·k columns... bounded by (2t+3) × (k) per certificate — all trivially
  feasible per the README runtime basis.
- **Certificate artifacts (frozen per campaign dir):** `instance.json` (m, n,
  t, k, D, seed, shelf: G coefficients, full support in order, guard-pair
  (p,q) with Δ_{p,q} ≠ 0, code and data hashes); `results.json` (per
  instance: table of c, Nfam(c), per-(point,order) ranks r_{a,j}, order-0
  witness booleans, guard check booleans); `echelon/` — per instance and c, the
  row-echelon form (snapshots) of (i) the family evaluation matrix at c = 2t+3,
  (ii) two sampled per-point rank maps (largest and smallest k); a SHA-256 of
  every artifact; the exact git-less code hash (sha256 of src/ files recorded
  in instance.json).
- Rebuildability: the polynomial F itself at t = 96 is k×(D+1) ≈ 7.7M field
  elements — NOT stored; stored instead are seed, G, and support (n field
  elements each direction) which suffice to re-derive F bit-exactly by the
  recorded deterministic construction.

## 5. Pass/fail (verbatim from README — fixed before any run)

- **Pass (Apon):** in all nine instances each additional held position in the
  critical window contributes **exactly one** new independent equation, and the
  count stabilizes exactly at 2t+3. Operationally in quantities (1)–(3):
  Nfam(c) = 2t+3−c for c ≤ 2t+3 and 0 at c = 2t+4 (read with the direction:
  the increment Nfam(c−1) − Nfam(c) = 1 for c ≤ 2t+3), all per-point ranks
  r_{a,j} = 1 at order 0 (nondegenerate) and ≤ 1 at orders 1..3.
- **Refutation (Vedenev):** any instance where a held position contributes
  ≥ 2 new independent equations below the 2t+3 threshold → exhibit the
  instance, its Goppa polynomial, its support ordering, and the echelon
  certificate. **Escalation (hard, from repo constraints):** any such event is
  messaged to Main via hub BEFORE being written anywhere as established.
- Evidence labels: what the echelon computations show is MACHINE-VERIFIED;
  the depth-uniform Lemma 7 claim at ladder scale is CITED-DEPENDENCY;
  extrapolation to standardized parameters would be COMPUTATIONAL-EVIDENCE
  only (a finite-n theorem with N, never asymptotic).

## 6. Scope notes fixed in advance

- Gate A's nine-instance census is a finite theorem with its N (nine instances
  × four held-counts): it decides the *measured collapse behavior* on these
  instances; it does not prove c_need = 2t+3 in general (Apon explicitly does
  not claim that upper bound either) and does not prove c_need = 2t+4 for any
  standardized parameter set.
- If any guard fails to hold on a generated instance (e.g. irreducibility
  rejected, binary-valued check fails), that instance is regenerated with the
  next seed in the recorded stream and the substitution is logged; guards are
  not adjusted.

— Mceliece agent, 2026-08-29.

---

# ADDENDUM 2 — guard chain revision (committed 2026-08-30, BEFORE the re-run)

Committed before implementing or running anything below. Fixes the guard
chain, the per-row evidence labels, and the abort conditions in advance.

## (i) delta by the exact Lagrange route, with ABORT on degree failure

F is constructed as f_j = SUM_i Y[j,i] * lam_i^{-1} * L_i, with
L_i = Pi/(Z - a_i) / Pi'(a_i).  Off-diagonal vanishing L_i(a_l) = 0 (l != i)
holds BY CONSTRUCTION: Pi/(Z-a_i) retains the factor (Z - a_l), and the
division is exact (zero remainder, asserted per i).  Therefore

    lam_l * f_j(a_l) = Y[j,l]   for ALL (j,l)

follows EXACTLY from the n scalar checks L_i(a_i) = 1, at cost O(n D)
instead of O(k n D).  This replaces both the earlier full k x n grid check
(m = 6,7,8) and the row-wise pivot+8-sample check (m = 9 onward).

LOAD-BEARING DEGREE CONDITION, checked per instance, ABORT (not warn) on
failure: max_j deg(f_j) <= D.  This is what licenses the Lemma-2 degree
argument (both sides of Apon (12) have degree <= n-1: LHS terms
Pi/(Z-a_i) have degree n-1; RHS G^2 f_j has degree <= 2t + D = n-1).  The
bound is essential, not decorative: f and f + Pi agree at every support
point yet differ, so without deg <= D the agreement-at-n-points argument
collapses.  A run whose degree check fails produces NO delta verdict.

## (ii) uniform re-run across ALL THIRTEEN instances, ONE code path

The delta guard is re-run for every one of the 13 certified instances
(m = 6,7,8,9,10,11) through a single implementation, so the ladder rests
on one standard verified by one code path.  Re-running only m = 9 (the row
whose standard I silently changed) would leave two implementations behind
one claim, which is not acceptable.  Per-instance results are reported.
If any of m = 6,7,8 fails to reproduce under the new route, that is a
FINDING: escalate to Main BEFORE writing anything.

## (iii) alpha stays MEASURED on the certified ladder; split by row

alpha (the differential identity Pi F' + Pi' F = G^2 F^(2)) is the
INDEPENDENT CHECK ON APON LEMMA 3 ITSELF — it is what would catch the
lemma being wrong, and it is the evidence behind the headline phrase
"confirmed by measurement".  Therefore:

  * m <= 11 (the 13 certified instances): alpha remains a per-instance
    MEASUREMENT, unchanged.  It is already done and affordable at these
    sizes.  Its label does not change.
  * m = 12, IF AND ONLY IF measuring alpha there proves unaffordable:
    alpha may be DERIVED from the exact delta plus Apon Lemma 3, and that
    row's alpha is then labelled CITED-DEPENDENCY.  This is a different
    KIND of evidence, not a different confidence level, and the artifact
    must state it in one sentence: at m <= 11 the mechanism is MEASURED;
    at m = 12 it is INFERRED from a cited lemma given a measured delta.

No global relabelling of alpha.  A reader must be able to see which rows
test Apon and which assume him.

## (iv) everything else unchanged

The m <= 11 labels, thresholds, conventions, flag depth R = 4, pass/fail
criteria and the enumerated claim of record are unchanged by this
addendum.  No universal ("m <= N complete") phrasing is used for any N.

-- Mceliece, 2026-08-30, before the re-run.
