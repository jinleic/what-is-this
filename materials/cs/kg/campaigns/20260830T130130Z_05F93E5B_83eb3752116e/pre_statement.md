# Gate C pre-statement addendum (agent KgGateC-2, 2026-08-30T~12:00Z, committed BEFORE any septic computation)

Scope: upper side (a) septic extension + lower side (b) degree-5 analogue of the affine
coefficient inequality. Anchor, tripwires, and labels per Main's assignment. This addendum
supersedes the Gate C section above ONLY where it restates box/budget; everything else stands.

## Literature scan loaded first (Main's directive, 2026-08-30; REPORTED unless noted)

- **Hei26b — arXiv:2606.00247v1 (owner-read in full this session).** Heilman, "An Upper
  Bound on Grothendieck's Constant". THE COLLISION CANDIDATE. Planar Krivine schemes
  f(x)=sgn(x2−η h_d(x1)) (d=5 symmetric; d=3 anti-pair), one random plane, NO mixing
  coordinate. Perturbative in q=η²: H_q(z)=arcsin z + q·(2(σz^d−z)/(π√π√(1−z²))) + O(q²),
  inverse B_q(z)=sin z + qΨ_d(z)+O(q²), A_q(L)=1−(24 or 800)√2 q²/π+o(q²) at
  L=log(1+√2). Certified result (§16, Sage/interval): KG < Krivine − 1.013e−5, at the
  SINGLE point η=0.04249900400783211, d=3 anti-pair. Their theorem-1.9 machinery proves
  A_q(L)<1−7.25e−6 for that one η — a perturbative neighborhood certificate, NOT an
  optimizer search, NOT a K_G lower bound, and NOT applicable at fixed non-small η.
- **JM26 — arXiv:2603.30039v1 (owner-read in full).** Lower bound: Hermite projection
  games A = Π1 − λ*I − εΠ3, KG >= KDR + 1e−12. Mechanism: every DR near-extremizer has
  Ω(1) weight on degree-3 Hermite coefficients (Lemma 3.3: E(Π3 f)(Π3 g) >= 0.046), so a
  small cubic perturbation raises val (Theorem 4.1). HOLDING a lower bound via the SDP
  side; no statement about admissible γ of explicit rounding schemes — the paper does not
  use the inverse-majorant chain at all.
- **Hei26a — arXiv:2603.22616** (abstract; Main's scan): KG >= KDR + 1e−26, same
  perturbation strategy, different analysis. Lower side only.
- **XZ — arXiv:2608.14817** (abstract; Main retracted his earlier cap-hypothesis about it
  after owner-reading; I verified the abstract's load-bearing sentences myself:
  "One question in [BMMN] attempts to determine the Grothendieck constant through
  alternating Krivine rounding schemes arising from König's bilinear form in high
  dimension" and "this gives a negative answer to the high-dimensional aspect of the
  question in [BMMN]"). NOT load-bearing for Gate C: it does not cap finite-dimensional
  families (König-value 1 in the limit is an upper-side improvement direction per
  Corollary 1.3 as Main quoted: "There is a randomized Krivine scheme such that
  K_G < pi/(2 log(1+sqrt2))"), and our scheme is a 2-coordinate tensor scheme, not a
  König-form instance. Recorded as adjacent, non-load-bearing.

**Verdict on duplication (answers Main's three questions):**
(a) Heilman's family is NOT 2608.11158's family: his f,g live on ONE random plane with
threshold structure sgn(x2 ∓ η h_d(x1)); ours mixes a ρ-correlated coordinate w with
sgn-partitions in the transverse Gaussian — his expansion point is arcsin (κ=1/2, i.e.
b1 = π/2·E[U V] with E U² = 1/2 conditioning), ours is the two-coordinate ρ-scheme with
A-grid coefficients. Both are Krivine rounding schemes; the parameterization, the
sufficiency chain (his A_q(γ)=1 exact vs our γ+Δ_H<b₁ head-only), and the exact constants
differ. HOWEVER his §16 certificate (Krivine −1.0e−5) is the SAME PROGRAM (Hermite
thresholding + interval arithmetic) and strictly weaker than our headline 1.7818413e+0
(Krivine −3.7265e−4 at Gate A's γ*); no conflict, but our Gate C must cite him as the
independent-interval-arithmetic precedent at d∈{3,5}.
(b) Heilman's interval computation does NOT cover our septic box: it is q∈{e^−450}
perturbative (his §8) plus ONE exact-η interval run at η=0.0425, d=3. Our box is
(s3,s5,s7)∈[0,0.6]³ to depth 7 with certified head b1 − Δ_H; his machinery cannot reach
non-small η and only handles planar threshold signs, not ρ-power-correlated coordinates.
(c) Gate C(b) is NOT answered by him: he has no b5-constraint analogue (that is a
universal FCC (finitely constrained-coefficient) statement in OUR chain, cf. §10.2
Proposition 10.2 which uses only b3); JM26 works on the SDP-value side of Π3 perturbation
and never writes admissible-γ or inverse-majorant constraints. Both are ADJACENT, not
duplicating. C(b) stays open as assigned.

## (a) Septic family — DERIVATION anchored in the paper's own text (quoted, not assumed)

Eq (6) at paper line 672 (raw dump): "ρD(t) = (t + Σ_{d∈D} σd s_d² t^d)/V_D" — the general-D
correlation map. Line 668: "σ_d = (−1)^{(d−1)/2}, V_D = 1 + Σ_{d∈D} s_d²". Line 685: 
"V = 1 + s₃² + s₅², ϑ = η/√V". Line 687 (raw): "ρ(t) = (t − s₃² t³ + s₅² t⁵)/V" — matches
(6) with σ_3 = −1, σ_5 = +1 for D={3,5}. Line 1721 repeats (8) with the explicit 1/(1+s₃²+s₅²).
Line 857: "D := t d/dt"; line 1826 rho_j := D_t^j ρ; line 1816 Φ_{k+1} = D_t Φ_k − t ∂x∂y Φ_k.
Septic instance adopted (D = {3,5,7}): σ_7 = (−1)^3 = −1 ⇒
**ρ(t) = (t − s₃² t³ + s₅² t⁵ − s₇² t⁷)/V, V = 1 + s₃² + s₅² + s₇²**.
Correspondingly rho_j(t) = (t − 3^j s₃² t³ + 5^j s₅² t⁵ − 7^j s₇² t⁷)/V — the pattern
(m^j coefficient on t^m, sign σ_m alternates with odd m) follows from D_t^j t^m = m^j t^m
(line 858: "the Euler operator … multiplies the coefficient of t^m by m") applied to (6).
Constraints on (s3,s5,s7): same odd-Hermite bracket as D→∞ in line 668 — the paper restricts
to D ⊆ {3,5,7,…} i.e. allows {3,5,7} directly; my box keeps s_d ≤ 0.6 which keeps
Σ_d s_d² ≤ 1.08 bounded, and V ≥ 1 ⇒ ρ allowable for every box point (|coefficients|
sum = 1 exactly by construction regardless of s values).

**Septic c-table:** the (p,a) ∈ ℐ rows and the (1,2),(2,1),(1,1),(1,0),(2,0),(3,0)
involving rho_j generalize by substitution 5^j → 5^j, 3^j → 3^j, 7^j → 7^j inside rho_j
only — c_{1,0} = rho_3(new), c_{2,0} = 3 rho_1 rho_2, c_{3,0} = rho_1³,
c_{1,1} = −3t(rho_1 + rho_2), c_{2,1} = −3t rho_1², c_{1,2} = 3t² rho_1 — exactly as
printed in paper lines 1836-1852 with rho_j = D_t^j ρ (line 1824), only the rho-polynomial
underneath changes. Rows (0,*) (t, 3t², −t³) don't involve rho. G_{p,a} machinery from
d3h_kernels/d3h_dxdy is scalar in r = rho(t) — unchanged. The leg2 dxdy shift uses
G_{p,a+1}: unchanged.

## (a) Septic sweep protocol [committed BEFORE the box is spent]

- **Pre-statement anchor (Gate A-style, MUST pass first):** specialize s₇ = 0 and
  parameters (s₃, s₅, ϑ) to the paper's exact decimals 0.34101124, 0.05276111,
  0.136419125/√V. The 49-point θ-sweep via d3h_sweep49_asm.py machinery with my
  independently-implemented septic c-table must reproduce **avg(S) = 123.377609** and
  **B3 = 14.2459830** ⇒ **B3² = 202.948032** (from B3² = (π²/6)·avg(S), paper line 1866
  chain), **avg(leg1) = 1.939156**, **avg(leg2) = 30.359613** to 6 digits. NO extension
  runs before this reproduces. If it misses, STOP, do not proceed.
- **Tripwires (asserted at runtime, residuals reported, abort on failure above 1e−10):**
  (T1) pi-periodicity. All rho_j have odd exponents only ⇒ rho_j(−t) = −rho_j(t) for
  s₇ ≠ 0 as well, so S(θ+π) = S(θ); on the uniform 49-node grid we check S(0) = S(π)
  (measured target 9.93e−16 relative on the quintic reference).
  (T2) reflection. S even + π-periodic ⇒ S(θ) = S(π − θ); 24 independent node-pair
  equalities on the 49-pt grid (measured target 7.94e−15 relative on the quintic
  reference). Reflection is the cheaper catch-all tripwire.
  (T3) tail-shape guard: min S at θ = π/2 exactly (grid midpoint) as an additional
  invariant of the septic family (the reflection critical point on an odd grid).
- **Parameter box (committed, no adjustment after the run):** (s3, s5, s7) ∈
  [0, 0.6]³. Depth 7 per axis ⇒ 2^21 = 2,097,152 leaf cells; NOT feasible in the
  ≤ 200,000-box budget, so the committed protocol is depth-7 GOAL with a two-phase
plan: (i) coarse pass 6³ = 216 uniform cells, retain the top-K (K=64) by certified head
  upper estimator; (ii) local subdivision to depth 7 (7 refinement rounds on the
  survivors), total live boxes capped at 200,000. Budget overrun ⇒ freeze the reachable
  verdict and ship it; no post-hoc box change to chase a number.
- **Objective (interim, CONDITIONAL):** maximize head-only γ*_head(s3,s5,s7) =
  b1 − head(3..251) over the box, computed at prec 256 in Arb from the frozen
  A-grid machinery (DEPENDS ON ϑ = η/√V_with_s7 — same formal path). The Δ_H tail
  needs a certified septic ‖D³H‖ — my float-scale machinery is the sweep above; the
  certified part stays CONDITIONAL until its tail boundary pass succeeds. Verdict per
  the committed outcome taxonomy: certified-head improvement ⇒ CONDITIONAL-HEAD with
  an explicit open tail; certified cap ⇒ close with MACHINE-VERIFIED label on the
  swept set ONLY (rule 7: exact swept set stated, not universal).
- **Reporting:** every evaluated box value reported including misses; no
  near-target-coincidence explanations; retraction of any false start uploaded
  to this file as an inline correction, not a silent removal.

## (b) Degree-5 analogue of the affine constraint [open question as of this addendum]

Status: Formulation. The paper's Proposition 10.2 chain uses only b3 ≥ 2b1 − 11/6; a
b5-tightening would need a UNIVERSAL pair-inequality b5 >= α b1 + β b3 + δ (affine in
coefficients, mixture-stable, closed under coefficientwise limits) certified by the SAME
1-D Gaussian-interval program. There is no such statement in the paper, none in Hei26b
(SDP side + threshold schemes with sgn structure,19§16.1), none in JM26 (Hermite
projection games, Π3 perturbation, no inverse-majorant chain). Protocol: (i) derive the
candidate b5 analogue from the ternary u machinery extended to degree 5
(V5(u) = Σ_{j≤5}(E[u ψ_j])², β same, threshold-form optimizers); (ii) interval-certify
the resulting 2-parameter (a_0..a_5 unit norm, c) family analogous to Gate B's (p,c)
cover at degree 3; (iii) combine with Proposition 10.2's M_H inequality to see if
γ <= 11/12 can drop. HONEST EXPECTATION (recorded now, refutable): the degree-5 term is
NOT functionally independent — heuristically b5 ~ (small, sign-symmetric in b3) for the
paper's sign families, so it likely does NOT tighten γ < 11/12 below 11/12; the
pre-registered outcome covers both branches.

## One-way doors (cannot adjust after results)

- Box = [0,0.6]³, budget ≤ 200k live boxes, depth-7 target, 6^3 coarse pre-pass,
  top-64 survivors, 49-pt θ/grid sweeps: FIXED as of this file's commit.
- Anchors: s₇=0 must give avg(S) ∈ [123.37760][6-digit window], B3 = 14.2459830,
  leg1 = 1.939156, leg2 = 30.359613 (6 digits) else STOP.
- Tripwires T1/T2/T3 with abovementioned measured targets; abort above 1e−10.
- The A-grid regeneration for septic ϑ (needed for certified head) is priced before
  launch: reuse the ρ-power machinery at new ϑ ∈ ϑ(s7), head-only, 251 cutoff
  (paper line 1895). Head-only budget: ≤ 8 wall-hours; else freeze with CONDITIONAL.

## ADDENDUM 2 (agent KgGateC-2, appended BEFORE the corresponding runs) — Gate C(a) refinement stages

Committed after the anchor passed and after the certified head machinery was validated
(gate_c_fast.py reproduces b1 = 0.8815738220495995485818877 and
head = 1.132885992768967e-5 at the paper point, matching Gate A and paper Prop 6.3),
and before any search over the box. Reason for the addendum: the coarse-lattice-only
protocol of Addendum 1 is provably blind to the structure just measured, and saying so
in advance is required rather than silently changing method.

**Measured structure that forces the refinement (all certified, prec 256):**
- head(s) := sum_{3<=m<=251 odd} |b_m| is a razor valley, not a smooth bowl. At the paper
  point b_3 = +3.453873e-09 and b_5 = -9.043554e-11 -- i.e. the paper's TWO free
  s-parameters annihilate b_3 and b_5 to their decimal rounding; the surviving head
  1.132886e-05 is carried by b_11 = -5.743e-06, b_13 = -1.449e-06, b_15..b_23.
- b1 depends on the parameters ONLY through V: b1(V) = (pi/2)(A_{1,0}^2/V - A_{0,1}^2),
  strictly decreasing, with certified A_{1,0}^2 = 0.628071295048328 and
  A_{0,1}^2 = 1.54451017356478e-5, so db1/dV = -(pi/2)A_{1,0}^2/V^2 = -0.7878 at V_paper.
- Consequently ANY septic gain is bounded: since head >= 0, a septic point can beat the
  quintic point only if 0.7878*(V_new - V_paper) < head_paper = 1.13e-5, i.e. only inside
  the thin shell V_new - V_paper < 1.4375e-5. A 6^3 or even 128^3 uniform lattice on
  [0,0.6]^3 cannot resolve a shell of thickness 1.4e-5; the lattice pass therefore
  measures the landscape but CANNOT decide the gate.

**Stages added (budgets fixed here):**
- (S1) CERTIFIED CAP by branch-and-bound on the pre-registered box with the bound
  gamma*(s) <= b1(V_min(B)) - head_lo(B) for all s in B (valid because tail >= 0 and
  head >= 0). head_lo(B) computed from interval rho coefficients; boxes with bound
  <= gamma_ref are certified incapable of beating gamma_ref. Depth 7, <= 200,000 boxes,
  8-way bisection, gamma_ref in {gamma_paper = 0.881545409,
  gamma_gateA = 0.881557917504162}. Cheap m-truncated variant of head_lo permitted
  (Sum_{3<=m<=M} |b_m| with M = 23 is still a certified LOWER bound on head since every
  term is non-negative; measured to carry 97.8%-99.99% of the head at probe points).
- (S2) VALLEY CONTINUATION: on the b_3 = b_5 = 0 manifold (the paper's own construction,
  continued in s7), solve for (s3, s5) by Newton at each fixed s7 in the committed grid
  s7 in {0, 0.0002, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01, 0.015,
  0.02, 0.03, 0.05}, and report gamma_head at EVERY node including misses. Also test the
  septic-natural triple-annihilation b_3 = b_5 = b_7 = 0 (three equations, three
  unknowns) and report its (V, head, gamma_head) whether or not it wins.
- (S3) FREE LOCAL SEARCH: Nelder-Mead on (s3,s5,s7) maximizing the certified
  gamma_head midpoint, started from the paper point and from the S2 winner, <= 400
  evaluations each, reported with every improvement step.
- (S4) The winner (if gamma_head exceeds gamma_gateA) gets a float 49-pt S-sweep for its
  tail (COMPUTATIONAL-EVIDENCE only) and a certified head restatement; if no point
  exceeds gamma_gateA, the verdict is the certified cap from (S1) plus the S2/S3 tables,
  stated as a cap over the swept set with the swept set named exactly (rule 7).

**One-way door:** the objective is gamma_head = b1 - head (certified) with tail >= 0
handled by inequality; the reference values gamma_paper and gamma_gateA above are fixed;
no post-hoc reference change. If the septic family loses, that is the result.

## ADDENDUM 3 (agent KgGateC-2) — CONVENTION CORRECTION recorded before the runs it affects

**Defect found in Addendum 1/2 as written, and its fix.** A_{a,b} = E[psi_b(X) q_a(vartheta
psi_3(X))] (paper lines 1717-1719) depends on vartheta, and the paper sets
vartheta = eta/sqrt(V) (line 685). The certified A-grid (grid_M251.json) is computed at the
paper's vartheta_paper = 0.136419125/sqrt(V_paper) = 0.128957369921412. Therefore a septic
sweep that holds eta fixed moves vartheta (since V moves), and the imported A-grid would no
longer be the right input — the head/b1 machinery would be evaluating a scheme whose
A-coefficients it does not have. My gate_c_fast/gate_c_bnb/gate_c_valley code uses the
A-grid as CONSTANTS, i.e. it silently holds vartheta FIXED.

**Resolution (chosen deliberately, not to fit a number):** the swept family is
    f(w,x) = sgn(w + vartheta psi_3(x)),  g(w,x) = sgn(w - vartheta psi_3(x)),
    rho(t) = (t - s3^2 t^3 + s5^2 t^5 - s7^2 t^7)/V,  V = 1 + s3^2 + s5^2 + s7^2,
    **vartheta FIXED at vartheta_paper = 0.128957369921412** (equivalently eta(s) =
    vartheta_paper*sqrt(V(s)) — legitimate because eta is a free parameter of the paper's
    own construction, line 675 "for any eta in R", so (s3,s5,s7,vartheta) and
    (s3,s5,s7,eta) parameterize the same family bijectively).
At s7 = 0 with the paper's (s3,s5) this is EXACTLY the paper's scheme, so the anchor is
unaffected (vartheta_paper = eta_paper/sqrt(V_paper) identically).

**Consequence for the claim, stated now:** every Gate C(a) certified statement is over the
3-D slice {(s3,s5,s7) in [0,0.6]^3, vartheta = vartheta_paper}. The vartheta-axis is NOT
swept. Sweeping it requires regenerating the certified A-grid at each new vartheta
(16,002 parity-allowed 1-D Gaussian integrals with a+b <= 251; priced in Addendum 1 at
<= 8 wall-hours, NOT run in this campaign). This is a Rule-7 scope statement, not a
claim of universality.

**Also declared before running (supplementary, same bound, same box, same budget cap):**
depth-12 and depth-16 refinements of the identical certified branch-and-bound. Rationale:
the pre-registered depth-7 run finished at 2,409 evaluations, i.e. 1.2% of the 200,000-box
budget, and depth 7 provably cannot resolve the 1.4e-5-thick valley (one depth-7 cell in
s3 moves b1 by 2.5e-3, ~200x the entire head budget). The deeper runs are the SAME test at
finer resolution, reported alongside — not a replacement for the pre-registered one.

**Float cross-check declared:** an independent float A-grid (my own Gauss-Hermite
quadrature of A_{a,b} = E[psi_b(X) q_a(vartheta psi_3(X))]) will be compared against the
imported certified grid entries (second-route requirement), and used ONLY to map the
vartheta-sensitivity as COMPUTATIONAL-EVIDENCE.

## ADDENDUM 4 (agent KgGateC-2) — residual-slab refinement, declared before running

The pre-registered depth-7 run and the supplementary depth-10 run both terminate with a
survivor hull contained in s3 in [0.33223, 0.34160], s7 in [0, 0.019922]. Depth-12/14 runs
on the FULL box exhaust the 200,000-box budget mid-tree (depth-first, 8-way splits) and
therefore prune less volume than depth 10; that is the honest budget-capped outcome and is
reported as such. To shrink the residual instead of re-spending the budget on already-pruned
outer regions, the SAME certified bound is now applied to the survivor region as its own
root box:
    R := [0.325, 0.350] x [0.00, 0.080] x [0.000, 0.020]
at depths up to 16, budget <= 200,000 boxes, references unchanged (gamma_gateA,
gamma_qhead). This is subdivision of survivors — step (ii) of Addendum 1's two-phase
protocol — not a new box: R is a superset of both survivor hulls above, so
"pruned outside R" + "pruned inside R" composes into a statement about the whole
pre-registered box. Any part of R that survives at budget is reported with its exact
hull and volume.
