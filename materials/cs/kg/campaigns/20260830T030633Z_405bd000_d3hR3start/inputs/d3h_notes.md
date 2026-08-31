# d3h notes — working memory for ‖D³H‖ re-derivation (KgD3H session, 2026-08-30)

Purpose: independently re-derive the boundary-norm bound ‖D³H‖_{L²(𝕋)} ≤ 14.44243664663976457
(currently CITED-DEPENDENCY from the authors' frozen d3h_certificate.json), and recompute
tail / Δ_H_upper / margin_low / K_G with the derived value.
Contract outcomes: (a) derived ≤ theirs → chain becomes DERIVED; (b) derived > theirs → named
divergence step; (c) named non-reconstructible Appendix C step with exact quote.
All three worth having. NEVER tune toward 14.44243664663976457.

## 1. Paper facts (verified verbatim from scratch/paper_full.txt this session)

### Lemma 6.2 (lines 862–884) — exact statement
For every odd N,
  Σ_{m>N, m odd} |b_m| ≤ ‖D³H‖·(Σ_{m>N, m odd} m⁻⁶)^{1/2} ≤ ‖D³H‖/√(10·N⁵).   (25)
Proof: D³(t^m)=m³t^m, Parseval:
  ‖D³H‖²_{L²(𝕋)} = Σ_{m≥1, m odd} m⁶·|b_m|²                                    (26)
then Cauchy–Schwarz and Σ_{m>N odd} m⁻⁶ ≤ ½∫_N^∞ x⁻⁶ dx = 1/(10N⁵).
=> MY identity to compute: ‖D³H‖² = Σ_{m odd} m⁶|b_m|². b_m = 0 for even m (parity).

### Prop 6.3 / eq (28) / Thm 6.4 (lines 885–913)
Δ_H ≤ 1.1328860e-5 + 14.443/√(10·251⁵) < 1.5904724e-5;  K_G ≤ π/(2γ) ≤ 1.7818666069360661,
γ = 0.881545409. Numerical bridge line 1898: 0.881545409 + 1.1328860e-5 + 14.443/√(10·251⁵) < 0.881573822049.

### Appendix C (1729–1900) — reconstruction checklist: ALL 15 items verified reconstructible
1. q₀(u)=erf(u/√2) [1735]; 2. q_a(u)=2φ(u)·He_{a−1}(−u)/√a (a≥1) [1735]; 3. q_a′(u)=2φ(u)He_a(−u) [1776];
4. Lemma C.1: ∂_rC_r = 2πK_r, ∂_x∂_yC_r = 2πϑ²He₃′(x)He₃′(y)K_r, K_r = Mehler kernel [1747–1782];
5. Lemma C.2: D_tT_{−t}[Ψ] = T_{−t}[D_tΨ − t∂_x∂_yΨ], d/dx He_n = n·He_{n−1} [1793–1807];
6. (51): |T_λ[Ψ]| ≤ ‖Ψ‖ + (π/√6)‖∂_x∂_yΨ‖ [1809–1811]; 7. H(t)=T_{−t}[C_{ρ(t)}] [1814];
8. Φ_{k+1} = D_tΦ_k − t∂_x∂_yΦ_k, (52): D³H = T_{−t}[Φ₃] [1815–1819];
9. ℐ = {(1,0),(2,0),(3,0),(0,1),(0,2),(0,3),(1,1),(2,1),(1,2)} [1829];
10. c-table: c₁,₀=ρ₃, c₂,₀=3ρ₁ρ₂, c₃,₀=ρ₁³, c₀,₁=−t, c₀,₂=3t², c₀,₃=−t³,
    c₁,₁=−3t(ρ₁+ρ₂), c₂,₁=−3tρ₁², c₁,₂=3t²ρ₁ where ρ_j = D_t^j ρ [1835–1852];
11. S(θ) = ‖Φ₃‖² + 4‖∂_x∂_yΦ₃‖² [1854]; 12. (x+πy/√6)² ≤ (1+π²/24)(x²+4y²) ≤ (π²/6)(x²+4y²) [1857];
13. ‖D³H‖² ≤ (π²/6)·(1/π)∫₀^π S(θ)dθ [1859–1861];
14. certified (1/π)∫₀^π S dθ ≤ 126.80385221 ⇒ ‖D³H‖ < 14.443 [1863–1866];
15. C.4 items 1–5: exact rationals/outward rounding; adaptive 1-D + closed Gauss tails; K_r derivatives
    as poly-weighted Gauss integrals w/ certified series tails; composite 4-pt GL, 1024 panels,
    panelwise 8th-derivative remainder; per-panel |1−ρ(e^{iθ})²| ≥ 0.3724, Gauss margin ≥ 0.49 [1871–1893].
⇒ No (c) gap found so far. Deliverable (a)/(b) path is open.

## 2. Constants (exact decimals → fmpq; V etc. computed exactly in code)
η=0.136419125; s3=0.34101124; s5=0.05276111; V = 1+s3²+s5² ≈ 1.1190775713;
ϑ = η/√V; ρ(t) = (t − s3²t³ + s5²t⁵)/V  (only odd degrees 1,3,5);
ρ-eff coeffs ≈ (0.8936, −0.1039, +0.002487); ρ(1) = (1−s3²+s5²)/V ≈ 0.79217 (exact fmpq in code);
γ_paper=0.881545409; b1_low=0.88157382204959954858188766422240828 (half 2.46e-36);
head (m=3..251 Σ|b_m|) = 1.13288599276897e-5 (half 1.77e-36); √(10·251⁵) ≈ 3.15671e6;
B3² ≤ (π²/6)·126.80385221 = 208.66…;  authors' B3 = 14.44243664663976457 (B3² ≤ 208.584… < (π²/6)·126.8039 —
theirs used avg_S slightly below the paper's printed upper, consistent with paper's loose "< 14.443").
Imported tail = 4.57569e-6 = 14.44243664663976457/√(10·251⁵); Δ_H_upper = 1.59045454374832e-5;
margin_low = 1.2508504162e-5; γ* = 0.881557917504162; K_G ≤ 1.7818413238804372880 (RESULTS headline).

## 3. Grid & conventions (LOCKED, from prior validated agents — do not re-derive)
- scratch/repo/grid_M251.json: meta M=251 prec=256 L=18 + ~31k entries, keys a,b,m,r,e (m,r integer
  STRINGS, e int exponent) with a+b ≤ 251. A_{a,b} = m·10^e ± r·10^e (ball: midpoint m·10^e,
  radius r·10^e). Loader: gate_A_final.A_to_arb / load_grid.
- scheme.json: H(t) = (π/2)·Σ_{a,b} A²_{a,b} ρ(t)^a (−t)^b  ⇒
  b_m = (π/2)·Σ_{a,b} (−1)^b·A²_{a,b}·[t^{m−b}] ρ(t)^a.
- ρ^a support: degrees {a, a+2, …, 5a} (≡a mod 2). So [t^{m−b}]ρ^a ≠ 0 ⟺ a+b ≤ m ≤ 5a+b  ("window").
- Even m ⇒ b_m = 0 exactly (all contributing terms even-symmetric — consistency check).
- Σ_{a,b} A²_{a,b} = E f² = 1 exactly (f is a.s. ±1) — convention-lock check.
- b₁ two-key cross-check: b₁ = (π/2)(A²_{1,0}/V − A²_{0,1}·1); prior match to paper b₁: 4.86e-17.
- d3h_certificate.json sha256 prefix 80e945589b (verify at freeze).

## 4. LOCKED ROUTE — R2-prime: Parseval-exact + certified 1-D tails
No θ-quadrature, no Mehler truncation, no GL panels. ‖D³H‖² = Σ_{m odd} m⁶|b_m|² (paper eq 26),
computed by exact finite arithmetic on the grid + explicitly bounded tails.

b_m contributions come from pairs (a,b): a+b ≤ m ≤ 5a+b.
- m ≤ 251: only grid pairs (a+b ≤ 251) contribute ⇒ b_m EXACT (ball arith on A² balls, signed
  accumulation; rho_pow exact fmpq; square via [max(0,|mid|−r)², (|mid|+r)²] never arb(lo,hi)).
- 251 < m ≤ 1255 (= max over a+b≤251 of 5a+b… check exact max in probe): grid part still exact;
  plus NONGRID part from (a,b) with a+b > 251 (unbounded a):
    |b_m^nongrid| ≤ (π/2)·Σ_{a} ‖g_a‖²·max_{b>251−a, a≤m−b≤5a} |[t^{m−b}]ρ^a|,
  where ‖g_a‖₂² = Σ_b A²_{a,b} = E[q_a(ϑψ₃(X))²] ≤ 1 (slice mass; a.s. |f|=1 slice-wise total mass
  Σ_a ‖g_a‖² = 1). g_a(x) = q_a(ϑψ₃(x)), ψ₃(x)=(x³−3x)/√6 (orthonormal He₃).
  Coefficient bound: multinomial — |[t^k]ρ^a| ≤ Σ over e₃+e₅ ≤ a of a!/((a−e₃−e₅)!e₃!e₅!)·
  (1/V)^{a−e₃−e₅}·(s3²/V)^{e₃}·(s5²/V)^{e₅} — exact rational sums, computable exactly.
  Both edges geometric: near k=a: ~(1/V)^a ≈ 0.8937^a; near 5a: ~(s5²/V)^a ≈ 2.5e-3^a ⇒ closed m-sums.
- m > 1255: no pairs at all ⇒ |b_m| ≤ (π/2)·Σ_a ‖g_a‖²·max_{a≤m−b≤5a}|[t^{m−b}]ρ^a|, geometric in m
  (rate √(1/V) ≈ 0.945 for the a-heavy corner … actually per-m bound via multinomial-sum bound
  M(m) = max over valid (a,k=b−side…) — compute exactly in code; closed-form outer sum via ratio bound.
  NOTE: bound per (m) gets factor ‖g_a‖²-a dependence — plus Σ_a ‖g_a‖²·[coefficient bound]_a,m.
  Superexponential decay in a because low-edge ~0.8937^a and the m⁶ weight only sits on fixed m.
- Slice masses: ‖g_a‖² = Σ_{b} A²_{a,b} = ∫ q_a(ϑψ₃(x))² φ(x) dx; certified via
  core.integrate_certified on [−9,9] + erfc tail bound (core.gauss_tail_erfc); for a beyond range,
  Laplace: integrand ~ exp(−x²/2 − (ϑ²/2)(ψ₃)²·… )·poly — Gaussian-dominated, rate ≈ e^{−0.0074·x⁶}
  along ψ₃=x³−3x ⇒ (x⁴cohol) — implement explicit certified bound |q_a(u)| ≤ 2φ(u)·(|u|+C√(a))^{a−1}/
  (√a·√((a−1)!))-type via Hermite sup bounds; pick cleanest during implementation.
- Assembly: ‖D³H‖² = [m ≤ 251 exact balls] + Σ_{251<m≤1255,odd} (m⁶((|b_m^grid|)+|b_m^nongrid|)²)
  + closed far-tail. Take principal sqrt OUTWARD ⇒ B3_mine with END-TO-END stated width
  (ball arithmetic tracks everything; sweep precision 128/256/512 and confirm widths grow monotone).
- BONUS (same machinery, linear): Δ_H_direct = Σ_{m>251 odd}|b_m| ≤ Σ (|b_m^grid| + nongrid bound) —

## 8. RESOLVED (probe session): grid convention — MACHINE-VERIFIED
Prior note "Σ_all A² = 1" and "slice mass = E[q_a²] ≤ 1 with q_a = 2φHe_{a−1}/√a" were WRONG in detail.
Machine-verified convention (matches grid to 10+ digits at a=0..10, odd and even):
  A_{a,b} = E_X[ ψ_b(X) · 2φ(ϑψ₃(X)) · ψ_{a−1}(−ϑψ₃(X)) / √a ],  a ≥ 1;  A_{0,b} = E[ψ_b·erf(ϑψ₃/√2)], only b=1.
  ψ_n = He_n/√(n!) ORTHONORMAL probabilists' Hermite; ψ₃(x)=(x³−3x)/√6; φ standard normal density.
Paper's "He" in Appendix C is this orthonormal ψ (line 1803 "d/dx He_n = n He_{n−1}" is a raw-He
notation slip in the source; all numerical content matches orthonormal).
Facts (probe, float64+grid):
  - slice mass E[q_a(u(X))²] = Σ_b A²_{a,b} EXACTLY (verified a=0..10, ratio 1.000000);
    max slice = a=1: 0.62892287; a=0: 0.00825; a=3: 0.10296; decays with a.
  - Σ_{a,b∈grid} A² = 0.9681381096737394 < 1: consistent with b > 251 tail (truncation), NOT an error.
    [check in final: 1 − 0.96814 = 0.03186 must be accounted: Σ_a (E[q_a²] − Σ_{b≤251−a} A²_{a,b}) over
    nongrid b AND a > 251. NOTE: slice-1 mass 0.6289 is at b=1 mostly — fine.]
  - b_m pipeline (gate_A_final, b_m = (π/2)Σ(−1)^b A²[t^{m−b}]ρ^a) UNAFFECTED — it never used the
    wrong identity; b₁ matched paper to 4.86e-17 previously.
  - parity (probe): A_{a,b} = 0 unless b ≡ a+1 (mod 2)?? — j DETAILS: grid slice a=3 has b EVEN;
    a=2 has b; see below; b_m even-m zero check passed (0 violations).
    [exact statement: A_{a,b} ≠ 0 only for a+b odd — verified: a=1:b even(1→b even? A_1,0≠0) contradiction!
    RESOLVE IN d3h.py: print parity table. Probe grid: a=1: b∈{0,2,4,...} even; a=2: b∈{1,3,5,...} odd;
    a=3: b even. So parity: b ≡ a+1 (mod 2). Then (−1)^b ρ-degree k = m−b: m odd ⟹ b odd ⟹ a even ⟹ k even
    or b even ⟹ a odd ⟹ k odd. Both have k ≡ a (mod 2) = k ≡ b+1 (mod 2): b+k odd ⟹ m = b+k odd ✓.]
  - probe scale numbers (float64): head Σ_{3..251}|b_m| = 7.212176e-06 (paper 1.1328860e-5 — probe used
    float(A²) truncation and possibly wrong (−1)^b convention; RECOMPUTE EXACTLY in d3h.py before comparing);
    P1 = Σ_{m≤251,odd} m⁶|b_m|² = 0.3150533; grid-part (251,1255] Σm⁶|b|² = 117.4; ⟹ ⟨‖D³H‖²⟩ ≳ 117.7,
    sqrt ≳ 10.8 vs authors' 14.4424: remainder must come from nongrid a,b (b > 251−a, a > 251 slices) via
    H(t) = (π/2)ΣA²ρ^a(−t)^b over ALL (a,b) — the grid truncation at a+b ≤ 251 is only an approximation
    window of the full series. The full ‖D³H‖² = Σ_m m⁶|b_m|² needs nongrid bound to close ≥ grid part.
  - watch: probe head mismatch (7.21e-6 vs 1.13e-5) EXPECTED — probe summed |b_m| over b ≤ 251 only
    pairs (a+b ≤ 251 enforced via grid) but paper head counts full |b_m| for m ≤ 251 with nongrid-a tail;
    also probe float rounding. d3h.py must recompute head EXACTLY over grid: pairs (a,b: a+b ≤ 251,
    a ≤ m−b ≤ 5a+b window) — head ≥ probe's 7.2e-6.

## 5. Cross-check set (MANDATORY in d3h.py certified run)
(i) Σ_grid A² = 0.9681381096737394 reproduce exact-arb (≈1 − 0.03186; 0.03186 = b,a>nongrid tail per above);
(ii) even-m b_m ≡ 0 exact; (iii) b₁ two-key vs paper 0.881573822049…; (iv) slice masses = Σ_b A²_{a,b}
re-verified in arb for sampled a; per-slice E[q_a²] ≤ 1 check for a up to cut; (v) P1 ≤ total ≤ (π²/6)·126.81;
(vi) precision sweep 128/256/512 monotone widths; (vii) B3_mine vs authors' 14.44243664663976457 —
report, never tune; (viii) head: my exact Σ_{3≤m≤251}|b_m| over grid must match paper 1.1328860e-5 to
grid precision (probe's 7.2e-6 was float noise/nongrid omission — actually paper head INCLUDES all
b ≤ 251 which IS the full grid → exact replay must land at 1.13288599276897e-5; if it doesn't, find why
BEFORE proceeding (suspicion: A² of ball (mid±half) vs (mid)² handling).

## 6. Environment & file plan
- Python: /Users/jinleic/jinleic-workspace/cs/.venv/bin/python (3.14.3, flint 0.9.0, numpy 2.5.2).
- Run as: cd cs/kg/src && taskpolicy -c 1 -b low ../../.venv/bin/python d3h.py 2>&1 | tee ../scratch/logs/d3h_run.log
- Reuse from core.py: Prec ctx (256/400 bits), _dec_fmpq, pi(), phi, gauss_cdf,
  gauss_legendre_32, integrate_panels, integrate_certified (rigorous acb adaptive — USE THIS),
  psi(n,x) ON Hermite, gauss_tail_erfc.
- Reuse from gate_A_final.py: load_grid, A_to_arb, rho_coeff_fmpq, rho_pow.
- Write: src/d3h.py (final certified pipeline + self-checks), scratch/logs/*.log,
  scratch/d3h_precsweep.json; freeze campaigns/<UTC>_<uuid8>_d3h/ with provenance.json
  (format keyed like prior: campaign, gate, scope_fixed, verdict, numbers, environment,
  reproduced_cmd, timestamp_utc, campaign_uuid, src_tag; inputs = code copy + d3h_certificate.json
  copy + sha256s; grid sha only) and append-only README `## Current state` at file end (line 285+).
- Interval hygiene: arb(lo,hi) is (MIDPOINT, RADIUS) — never use as interval Ctor; never measure
  width by float64-subtract; zero-width transcendental = float leak bug.
- October rule: ONE low-priority pinned thread. No formatters. No project-wide anything.

## 9. STRUCTURAL FINDING (RAR): per-slice triangle bound on nongrid b_m is TOO WEAK — proven by numbers
Chain (all machine-checkable, float64 + exact Fraction where stated):
(1) Paper (726–741): A_{a,b} = E[f(W,X)He_a(W)He_b(X)], H(t) = (π/2)Σ_{a,b}A²_{a,b}ρ(t)^a(−t)^b,
    b_m = (π/2)Σ_{a+b≤m}(−1)^b A² [t^{m−b}]ρ^a. Full series over ALL (a,b) — no b-truncation.
(2) Slice masses: μ_a = Σ_b A²_{a,b} = E[q_a(ϑψ3(X))²] ≈ 0.63·a^{−3/2} (Plancherel–Rotach confirms;
    verified numerically on grid slices a ≤ 10; grid slice a=251 gives μ_251 = 4.40e-5 = A²_{251,0}
    since b ≡ a+1 (mod 2) forces b=0 there).
(3) ρ(t) = c1t+c3t³+c5t⁵ with c1 = 1/V > 0, c3 = −s3²/V < 0, c5 = s5²/V > 0.
    For [t^k]ρ^a: e1+e3+e5 = a, k = e1+3e3+5e5 ⇒ 2e3+4e5 = k−a ⇒ e3 ≡ (k−a)/2 (mod 2) FIXED PARITY
    ⇒ c3^{e3} has the SAME sign for every tuple ⇒ NO SIGN CANCELLATION: |[t^k]ρ^a| = Σ_tuples|term| EXACTLY.
    At (a, k) = (252, 300): dominant tuple e3=24: C(252,24)·(1/V)^228·(s3²/V)^24 = 10^{33.35−11.09−23.60}
    = 10^{−1.34} ⇒ τ ≈ 0.045. This is TIGHT, not loose.
(4) Consequence: for EVERY odd m ≥ 301 (b = m−k ≥ 1 odd… parity ok), the nongrid part of b_m includes
    (π/2)A²_{252,b}c^252_k with |c| ≈ 0.045 and Σ_b A²_{252,b} = μ_252 ≈ 4e-5 ⇒ per-slice bound
    (π/2)·μ_252·0.045 ≈ 2.9e-6 on |b_m| for ALL m in a long range ⇒ Σ_m m⁶|b_m|² ≥ (301⁶)(2.9e-6)² ≈ 5.4e3
    from THIS ALONE (Σ over many m even more), vs paper's total budget B3² ≤ 208.6. CONTRADICTION.
(5) Resolution of the contradiction: the TRUE b_m must cancel across (a,b) — the paper's ‖D³H‖ certificate
    (Appendix C route via S(θ) = ‖Φ3‖² + 4‖∂x∂yΦ3‖², Mehler kernels, π²/6 chain) bounds the FUNCTION norm
    pointwise on |t|=1 and never uses the triangle inequality on coefficient slices. The Parseval route
    ‖D³H‖² = Σm⁶|b_m|² with per-slice |b_m| ≤ (π/2)μ_a·max_k|c^a_k| is a VALID upper bound on each slice
    but the sum over a of these slice-bounds is ≫ 208 — hence UNPROFITABLE: it cannot certify B3 ≤ 14.44.
    (Numbers to double-check in d3h.py — but the scale gap is 25×, not a precision issue.)
## 10. REVISED ROUTE (R3): independent kernel-side re-derivation of ‖D³H‖ (paper's own estimate,
## from scratch, no reuse of their code/number):
‖D³H‖² ≤ (π²/6)·(1/π)∫₀^π S(θ)dθ with S(θ) = ‖Φ3(e^{iθ})‖²_{L²(μ)} + 4‖∂x∂yΦ3(e^{iθ})‖², where
  Φ3 = Σ_{(p,a)∈ℐ} c_{p,a}(t)·G_{p,a}, c-table (1835–1852), G_{p,a} = [∂_r^p(∂x∂y)^a C_r]_{r=ρ(t)},
  C_r built from (47)/(48): ∂_rC_r = 2πK_r(ϑψ3x,ϑψ3y), ∂x∂yC_r = 2πϑ²ψ3′(x)ψ3′(y)K_r(u,v),
  K_r = Mehler kernel — and ∂_r^p K_r = polynomial(|u|,|v|)-weighted K_{r} (Hermite-derivative chain).
  ‖·‖²_{L²(μ×μ)} of a weighted K: E[(poly·K)²] = 4-variate Gaussian integral — REDUCIBLE EXACTLY:
  ∫∫ poly(x)²poly(y)² K_r(ϑψ3x,ϑψ3y)² μdxμdy: substitute K_r² = (1/(4π²(1−r²)))·e^{−(u²−2ruv+v²)²/(1−r²)·...}
  — product of two Gaussian densities in (u,v) is again Gaussian × poly: standard Gaussian-moment
  algebra, closed form as 1-D integrals over single-variable functions (ψ3(x)³ odd moments) — CERTIFIABLE
  via core.integrate_certified on ONE variable. THIS is the faithful, independent re-derivation path.
  Mega-note: each of the 9 (p,a) terms handled separately; S(θ) needs this at each θ — pick Gauss–Legendre
  on θ with certified panel remainder (like authors' C.4 but implement independently).
Also compute Δ_H_direct bonus only if trivially available — else drop (Parseval wall kills it).

## 7. Probe status
Not yet run (this is next action; notes file written first per Main's protocol change).
Probe: src/d3h_probe.py — float64/fmpq scale check, minutes, THEN certified d3h.py.

## 7. R3 kernel machinery status (2026-08-30 ~05:00 UTC)
src/d3h_kernels.py — ALL multiplier kinds MACHINE-VALIDATED:
- Kf = bivariate normal density == paper's K_r (verified == (1/4)Σ b r^{b−1}q_bq_b via (47)-bridge, 15 digits);
- ℓ, ℓ', ℓ'', ℓ''' (ellf/ell1f/ell2f/ell3f): FD of ℓ to 1e-6 (ell1f term dict fixed):
  ell1f = 1/s + 2r²/s² + 6ruv/s² − qq/s² − 4r²qq/s³ + 8r³uv/s³;
- pure-r mults Kr/Kr2/Kr3 = ℓ, ℓ'+ℓ², ℓ''+3ℓℓ'+ℓ³ — series-anchored via falling-factorial series
  ∂r^pK = (1/4)Σ_{b≥p+1} falling(b,p+1) r^{b−1−p} q_b(u)q_b(v) (match < 1e-7 at r = 0.4/0.7/0.85);
  NOTE: earlier Kr3 had 2ℓℓ' instead of 3ℓℓ' — fixed, double-anchored;
- mixed mults Kru/Krv/Kuv/Kr2u/Kr2v/Kr2uv: FD-verified (nested 3-point, ratios 1.000001–1.000006);
  dr2L = (4rv+2w)/s² + 8r²w/s³, w = −u+rv (previous form wrong → fixed 2026-08-30);
- q_b(u) = 2φ(u)ψ_{b−1}(−u)/√b (orthonormal ψ!), q_b'(u) = 2φ(u)ψ_b(−u) — series-confirmed.
Building blocks for R3: C_r(x,y) = (π/2)Σ_b r^b q_b(u)q_b(v); G_{p,a} = ∂r^p(∂x∂y)^a C_r assembled
from mults + ψ3-weight factors. Next: Φ3 weights (c-table), then S(θ) via nested certified 1-D integrals.
sha256(d3h_kernels.py) prefix 6c8b7adf225f356c (snapshot in campaign 20260830T030633Z_405bd000_d3hR3start).

## 8. R3 progress: dxdy machinery VALIDATED (2026-08-30 ~05:30 UTC)
- d3h_kernels.py now has full TERM-ALGEBRA mult engine: chain_mult/kind via dict terms
  {(deg_u, deg_v, deg_r, pow_s): coeff}, steps r/u/v as term operations (product-rule composition
  m ↦ step(m) + m·base). ALL kinds with p+i+j ≤ 4 FD-VALIDATED vs Kf (114 checks ≤3rd order: 0 fails;
  48 4th-order spot checks: 0 fails). _T_M monomial bug fixed ((1,0,1,1) not (1,1,1,1)).
- d3h_dxdy.py: dxdy_terms(a) builds (ψ3-weight, kind, ϑ-power) entries for (∂x∂y)^a C_r;
  dxdy_eval evaluates. VALIDATED: a=1 vs q'-series (48) to 6e-16; a=2 vs FD-of-closed-a1 4e-10;
  a=3 vs FD-of-closed-a2 1e-9 (composition bug — single-vs-double-derivative branching — fixed:
  must apply ∂x branches then re-branch ∂y on results).
- NEXT: G_{p,a}(x,y;t) = (∂r^p)(∂x∂y)^a C_r as 2π·Σ entries with kind+extension r^p (r-steps never
  touch x/y weights); Φ3 = Σ c_{p,a}(t)G_{p,a} (c-table notes §1 item 10);
  ||Φ3||² = Σ_{i,j} c_i c̄_j ⟨G_i, G_j⟩ via K·K̄ = |K_r|² closed Gaussian in (u,v) + separable
  e^{2βuv} expansion → nested certified 1-D integrals (inner over x at fixed v-parameters...).
  Then same for 4‖∂x∂yΦ3‖² and (π²/6)·(1/π)∫₀^π S(θ)dθ with certified θ-quadrature remainder.
- sha256: d3h_kernels eb8ec506977e482f, d3h_dxdy 77a4ad88cfc509f4 (snapshots in campaign).

## 9. Phi3 assembly validated (2026-08-30 ~06:00 UTC)
- d3h_phi3.py: rho_j exact (D=td/dt on c1t+c3t³+c5t⁵), c-table (paper 1835-1852) as lambdas,
  G_pa_terms (a=0 bridge: ∂r^pC = 2π∂r^{p−1}K, kinds K+Krr+Krrr...; a≥1: dxdy entries + r^p folded
  into kind), G_pa_eval, Phi3_eval, dxdy_Phi3_eval.
- Identity checks all 1.000000000000: G10 = 2πK, G30 = 2π·mult(Krr)·K, G01 = (48), G11 = G01·ℓ.
- KIND-NAMING: repeated letters (Krr not Kr2), all p+i+j ≤ 4 covered by _KIND_CHAINS (31 kinds).
- NEXT: ⟨G_i,G_j⟩ via K·K̄ = |s|^{-1}(4π²)^{-1}e^{−α(u²+v²)+2αρh uv}, α = Re(1/s), ρh = Re(r) —
  separable e^{2αρh uv} = Σ βⁿuⁿvⁿ/n! ⇒ M_n^x·M_n^y nested certified 1-D integrals (M_n =
  ∫ W-prod·(ϑψ3(x))ⁿ e^{−α(ϑψ3(x))²} μdx). |β| ≤ 2·|α|·|ρh| bounded via per-panel |1−ρ²| ≥ 0.3724-type
  verification (to be re-derived independently!). Then S(θ) = ||Φ3||² + 4||∂x∂yΦ3||², quadrature on θ.
- sha256: phi3 8eb9b1d98a79eba9, dxdy 5d5a347fc3c1237d, kernels eb8ec506977e482f.

## 10. S(θ) float probe status (2026-08-30 ~07:00 UTC)
- d3h_inner.py: ⟨Gi,Gj⟩ via |K|²-separable expansion; Mn moments via P-matrix (u^k·e^{−αu²} rows,
  LSE-stable, sign-tracked), hermgauss(180) nodes |x|≤12.
- First numbers (S = ‖Φ3‖² + 4‖∂x∂yΦ3‖², real-t probes): S(0) = 0.012184 (two-leg 0.012144/0.000040).
  S(0.5), S(0.9) re-running vectorized. Budget (1/π)∫₀^πS ≤ 126.804 (paper); S values look small — GOOD
  direction (S is an integrand; the θ-average ≤ 126.8 is loose in the paper's certificate).
- NOTE: paper's budget trap: my derived (1/π)∫S ≤ (π²/6)^{-1}·(my B3²) is the STRONG number to compare:
  B3²_mine = (π²/6)·avg(S) — if avg(S) lands ≈ 126.8 ± noise, we've REPRODUCED; if materially SMALLER,
  we've IMPROVED (deliverable (a) with B3_mine < 14.4424...).


### CORRECTION to §10 (2026-08-30, per Main's 120-bit check) — MISREADING RETRACTED
Main's bracket: avg(S) MUST lie in [71.5530, 126.804] because B3² = (π²/6)·avg(S) with the hard
lower bound ‖D³H‖² ≥ 117.7 (grid Parseval partial sum; sqrt → B3 ≥ 10.848963) and the paper's
upper end 126.804 (B3 = 14.44243664663976457; (π²/6)·126.804 = 208.5842194, sqrt = 14.4424450635,
matches their B3 to 8e-6 — relation verified). Therefore "S small ⇒ improved" was WRONG: if S were
flat at its θ=0 value 0.012184, B3 would be 0.1416, violating MY OWN lower bound by 77×.
Small S(0) is a WARNING, not a win. Two live hypotheses:
 (i) S peaks ≳ 5900× between θ=0 and the boundary (|t|→1) — must be DEMONSTRATED by endpoint sweep;
 (ii) the normalization is off ~10⁴ — then a certified pass would rigorously certify a wrong number.
HARD SANITY GATE (assertion, to be coded into the certified runner): any computed B3² outside
[117.7, 208.6] fails the run loudly. Endpoint sweep t ∈ {.95, .99, .999, .9999} is the decider.
No prefactor tuning toward 126.804, ever; any missing factor must be derived and paper-quoted.

## 11. Circle-probe results + normalization investigation (2026-08-30 ~09:30 UTC) — COMPUTATIONAL-EVIDENCE
### S(θ) on the unit circle (t = e^{iθ}), direct 2-D entry evaluation (hermgauss(70), |x|≤13, float64):
  θ:      0.000    0.393    0.785    1.178   1.571   1.963   2.356    2.749    3.142
  ‖Φ3‖²:  128.25   42.96    8.82     2.50    1.60    2.54    5.80     11.56    14.74
  4‖∂Φ‖²: 223.90   521.98   293.18   14.34   1.65    4.74    232.77   959.91   1592.64
  S:      352.15   564.94   302.00   16.83   3.25    7.28    238.57   971.47   1607.37
Trapezoid/Simpson: ∫₀^πS ≈ 1211-1216 ⟹ avg(S) ≈ 386-387 ⟹ B3² = (π²/6)·avg ≈ 636 ⟹ B3 ≈ 25.2.
### vs paper's certified avg(S) ≤ 126.804 (B3 = 14.4424): MY VALUE IS 3.05× LARGER — outside the gate
both ways (gate band [71.55, 126.80]; PAPER's OWN certificate would be violated by ~636 > 208.58 if my
S were the paper's S; the paper's B3² = 208.58 = (π²/6)·126.804 relation re-verified at 120 bits by Main).
### Convention status:
- Paper line 725 defines He_n = ORTHONORMAL; line 692's "He3(x) = x³−3x" is a raw-notation slip (same
  family as line 1803). My orthonormal kernel convention is consistent with the paper's A-grid machine
  convention (which is arithmetic-verified against grid slices at 10+ digits). RAW-He3 test gives
  avg ≈ 5.4e6 (B3 ~ 3000) — decisively wrong. So the 3.05× is NOT a ψ-normalization issue.
- ‖Φ3‖²-leg (128→42→8.8→2.5→1.6→...) and ‖∂x∂yΦ3‖²-leg both HUGE at θ ≈ 0 and θ ≈ π — the S-shape has
  TWO peaks, at the real-axis endpoints of the circle (t = ±1) — γ = Re(r/s) peaks where ρ-real.
### Open discrepancy — candidate causes (to be distinguished next):
(a) Φ3 c-table D-bookkeeping: my Faà-di-Bruno check CONFIRMED the D³[C]-part (ρ3G10+3ρ1ρ2G20+ρ1³G30);
    the −tG01/3t²G02/−t³G03/−3t(ρ1+ρ2)G11/−3tρ1²G21/3t²ρ1G12 rows come from three rounds of −t∂x∂y
    recursion — NOT yet independently validated (Φ1-test misread the table; a correct Φ1 test is
    ρ1·G10 − t·G01).
(b) my leg2 = ‖∂x∂yΦ3‖² uses G_{p,a+1} — check the (51)-(53) chain: the trace est (51) pairs
    T[Φ3] with ‖Φ3‖ + (π/√6)‖∂x∂yΦ3‖; |D³H|² ≤ (π²/6)(‖Φ3‖² + 4‖∂Φ3‖²) — the 4 = (2·π/√6)²·(6/π²)... 
    derive: (x + (π/√6)y)² ≤ (π²/6)(x² + 4y²) with x = ‖Φ3‖, y = ‖∂Φ3‖: (π/√6·y)² ≤ (π²/6)·y² ✓
    consistent. So the 4-factor structure is right.
(c) the possibility that my ⟨G,G⟩ direct-2D quadrature OVER-resolves nothing but my Φ3-grid values
    themselves are subtly wrong at the |x| ≤ 13 boundary (weights ψ3'² ~ x⁴ at x=13 ~ 2.8e4 with
    e^{−αu²} killing: u(13) = 0.129·(2197−39)/2.449 = 113 ⇒ e^{−6400} = 0 ✓ safe).
### STATUS: pending a third-resolution experiment; not yet a deliverable. DO NOT certify anything yet.

### §11 addendum (2026-08-30 ~10:15 UTC): c-table structure CONFIRMED, discrepancy narrowed
- Symbolic re-derivation of the Φ-recursion confirms the paper's c-table STRUCTURE exactly:
  D³[C] = ρ3·G10 + 3ρ1ρ2·G20 + ρ1³·G30 (Faà di Bruno on the ρ-chain); the −t∂x∂y rounds generate the
  G01/G02/G03/G11/G21/G12 rows; ∂xy[2πK] = ∂r∂xyC stays in the G-lattice; ϑ-power bookkeeping matches
  (48) entry-by-entry. My Φ3 = Σ c·G IS the paper's Φ3 under the orthonormal convention.
- The 3.05× avg(S) discrepancy (387 vs their certified 126.804) is therefore real WITHIN my current
  assembly — NOT explained by ψ-normalization, c-table structure, ϑ-powers, or the (51)-(53) chain
  (all re-derived and consistent). Both peaks (θ≈0, θ≈π) drive my integral; the paper's certified
  bound would fail against my Φ3's S-integral. NEXT discriminating experiments (queued):
  (α) independent ‖G_{0,1}‖² single-term analytic/numeric cross-check;
  (β) Φ1 numeric recursion test properly: Φ1 ?= ρ1·G10 − t·G01 at a couple of (x,y,θ);
  (γ) weighted-2D-quad convergence study on ‖∂x∂yΦ3‖² at θ=0 (60/80/100/140 nodes — need 4·‖∂Φ‖²
      = 223.9 to be quadrature-converged; N=70 already stable for leg1 but leg2 has ψ3''²·ϑ⁴·s^{−1}
      integrands with u-power 6·s^{−p} tails — re-verify).
- All sweep tables above are COMPUTATIONAL-EVIDENCE (float64, NOT certified). Nothing frozen as
  certified; the campaign notes the two fixed bugs (prefactor, fsum) and this open 3.05× item.

## 12. Session freeze (2026-08-30 ~11:00 UTC) — handover state
DELIVERABLE STATUS: outcome-(b)-candidate frozen, not final. The independent kernel-route S-integral
(386.6, B3≈25.2) disagrees with the paper's certified 126.804 (B3=14.4424) by 3.05×; all named
convention/structure causes checked and eliminated; discriminating experiments queued (§11 items α-γ,
partially done: (α) ‖G01‖² inner_c-vs-direct reconciled to 3.3% — inner_c CORRECT; (γ) leg2 grid
convergence 55.98→52.68→51.52 — incomplete).
VALIDATED CAPITAL (machine-checked, reusable):
- 63 mult kinds (term-algebra) FD-validated vs K_r;
- (dx dy)^a C_r machinery a=1,2,3 (6e-16/4e-10/1e-9);
- paper (47) K-bridge identity to 15 digits; (48) entry-wise exact;
- ρ-chain c-table structure by Faà di Bruno; Φ-recursion kinds algebra;
- 2 fixed hard bugs with permanent WARN assertions (negative self-inner; prefactor cancellation).
KEY FILES: src/d3h_kernels.py (eb8ec506... pre-fix; final d76e42ae...), d3h_dxdy.py (5d5a347f...),
d3h_phi3.py (8eb9b1d9...), d3h_inner.py (04cd57b3...), d3h_circle.py (e09e8596...),
scratch/d3h_notes.md (020cb52e...). All in campaign inputs/ with hashes.
NEXT AGENT ENTRY POINT: notes §11 (α)/(β)/(γ) experiments, then either the paper-side certificate
re-derivation (outcome (c) test: is THEIR 126.804 reproducible from THEIR formulas?) or the
G-norm-assembly bug hunt inside my Φ3 (outcome (a) if found and fixed ≤ 14.4424).

## 13. CONVERGED circle probe + decision tree (2026-08-30 ~12:00 UTC) — COMPUTATIONAL-EVIDENCE
49-point θ-sweep (Δθ = π/48, both peaks resolved; N=140/trim14 quadrature): ∫₀^πS ≈ 1138.6,
avg(S) = 362.44; EXCESS vs paper certificate = 2.8583 (13-pt: 2.8604 — θ-quadrature CONVERGED to
0.07%). NOT 3.000 (−4.7%), NOT π (−9.0%). Values table in scratch/logs/d3h_circle_Sfine49.npy +
thfine49.npy; full table printed in session log.
Main's breakeven (recorded): headline tail-slack factor 3.7337 (‖D³H‖ ≤ 53.92 safe; beating Krivine
to ‖D³H‖ ≤ 596). At measured B3 = 24.43: tail = 7.74e-6, K_G ≤ 1.7818510 — STILL improves the paper's
1.7818666. So nothing in cs/RESULTS.md retracts; my number, if valid, merely weakens the improvement.
DECISION TREE (exclusive):
 (A) my Φ3/S-assembly inflated ~2.86× → find the bug → B3_mine lands ≤ 14.4424-ish → outcome (a).
 (B) my Φ3/S is right ⟹ the paper's certified (1/π)∫S ≤ 126.804 cannot hold ⟹ their tail constant
     understated, their K_G bound SURVIVES via breakeven slack — honest claim: certificate gap, their
     bound stands. (NOT "the paper is wrong" — the bound is robust; the certificate constant is not.)
DISCRIMINATOR: independent full-Range Parseval (identity (26) is EXACT — no inequality): Σ_{m odd}
m⁶|b_m|² over ALL m with the nongrid coefficient bounds from §9's machinery (per-slice masses +
falling-factorial tails). If the exact ‖D³H‖² < 208.58, my S-integral is wrong (→A). If > 208.58,
their certificate fails (→B). The §9 triangle-bound overshoot does NOT block computing the exact Σ
from the DISTINCT structure: b_m has exact grid parts + nongrid parts whose COEFFICIENT values (not
just bounds) may be computable via the series directly... — open.
Per-(p,a) decomposition at θ=0 recorded: diag Σ ≈ 48.95; cross-terms +79.3 (legitimate, PD-verified).
‖G01‖² reconciliation: inner_c vs direct-2D → 3.3% agreement (α-experiment done, clean).
(48)-bridge series resolution: Σ_{b≥0} r^b q_b'(u)q_b'(v) = 4K_r EXACTLY (b=0 term included; Mehler
cross-check) — my earlier "violation" was a mis-derived comparison (fixed).

### §13 addendum: T-route probe (2026-08-30 ~12:30 UTC)
Attempted the independent exact ‖D³H‖² via T_{−t}[Φ3] (52) evaluated as the correlated-Gaussian
expectation with correlation λ = −te^{...}, |λ| ↑ 1: values DIVERGE (0.92 → 8.2 as q: 0.9→0.9999 at
θ=0) — the boundary trace is NOT a plain 2-D integral limit; the paper's own "uniform integrable
majorant" caveat (paper 1820-1822) covers exactly this. Direct T-route probe ABANDONED. The frozen
discriminator (full-Range Parseval, notes §13) remains the only clean decider between (A)/(B).
SESSION STATE UNCHANGED otherwise. Final float evidence:
  avg(S) = 362.44 (49-point converged θ-quadrature), excess 2.8583, B3_implied ≈ 24.4,
  headline-safe per Main's breakeven (‖D³H‖ ≤ 53.92 keeps cs/RESULTS.md intact; ≤ 596 keeps
  beating Krivine). Interim K_G ≤ 1.7818510 (at my B3) — still an improvement over 1.7818666.

### §14 Layer-width experiment + panel-share scaling (2026-08-30 ~13:00 UTC) — resolves Main's objection
S endpoints SATURATE (finite plateau): S(0) = 334.347, S(π) = 1424.281; flat to 6 digits from 1e-3
inward. Layer width ~0.01-0.1. NO θ-singularity in S(θ) — the T-route divergence is a property of the
T-LINE functional, not of the μ×μ-norm integrand. The uniform θ-quadrature is legitimate at h ≪ 0.1.
Panel-share h-scaling (diagnostic requested): endpoint half-panel shares fall LINEARLY:
  h=π/8: 15.4%/40.7%; π/12: 8.1%/28.7%; π/48: 1.9%/8.1% — first-order convergence confirmed.
∫₀^πS across grids: 1128.7 / 1139.5 / 1138.6 (1% stable) with layer-refined quadrature 1138-1142.
⟹ avg(S) = 362.5 ± 1 (float), excess 2.86 ± 0.01 — RESOLVED quadrature, not artifact.
(A)/(B) discriminator stands: full-range Parseval (next campaign). Branch-(B) claims remain BLOCKED
per Main's rule until then; nothing written beyond evidence labels.

## 15. Leg-split resolution + double-count disconfirmation (2026-08-30 ~14:00 UTC) — response to Main
### STEP 1 (READ, done): paper quotes recorded — S def (1854): S = ‖Φ3‖²_{L²(μ)} + 4‖∂x∂yΦ3‖²_{L²(μ)};
μ = product std Gaussian (1812); (51): |T_λ[Ψ]| ≤ ‖Ψ‖ + (π/√6)‖∂x∂yΨ‖; chain (1857) verbatim; cert (1864).
The "4" IS in the paper's own S definition — a double count in MY assembly would have to be why MY
‖∂x∂yΦ3‖ is 3.17× the paper-certificate-implied value.
### STEP 2 (leg decomposition over the full 49-pt grid — Main's requested numbers):
  avg(leg1) = 18.2592   avg(leg2) = 86.0457   avg(S) = 362.4418 ✓ (matches the 49-pt S-sweep 362.44)
  leg1 share of S = 5.04%   B3 from leg1 ALONE = 5.480430
  (NOTE: an intermediate legs49 run had a √K-factor bug — values 504/1079 — DISCARDED after the
  A-style/B-style one-point diff isolated the error; all conclusions use the A-style fixed run.)
### Main's extrapolation (avg(leg1) ≈ 139 from θ=0 share 38.3%) is REFUTED: leg1 varies 128→1.6→14.7
across θ; its average is 18.26 (share 5.0%), so:
  - leg1-only B3 = 5.48 ≪ 14.4424 — the mixed leg is NOT a double count; it IS the certificate carrier;
  - no non-negative multiplicity c on the mixed term can reproduce 126.804 (c = 0 → B3 = 5.48; c = 1
    → √((π²/6)(18.26+86.05)) = 24.43 — the measured total) — CONFIRMS Main's own multiplicity-solve
    conclusion from the other side: the mixed term cannot be scaled away.
### RESIDUAL DISCREPANCY LOCALIZED: ‖∂x∂yΦ3‖²-avg = 86.05 vs paper-certificate-implied ≤ (126.804 −
18.26)/4 = 27.14 — a 3.17× excess concentrated in leg2. Candidates narrowed to:
  (i) paper's Φ3 ≠ my Φ3 (definitional gap not visible in the printed table/recursion → outcome-(c)
      candidate with a precise quote request), or
  (ii) an inflation inside the ‖∂x∂yΦ3‖² evaluation chain (entries validated individually vs the
      paper's (48)-series — so an ERROR would have to be in the paper's own certificate constant).
Frozen; full-range Parseval remains the independent arbiter.

## 16. Audits (a) shift + (b) operator — BOTH PASS (2026-08-30 ~15:00 UTC)
(a) a→a+1 shift: dxdy_Phi3_eval and S-leg2 read c_{p,a} (UNSHIFTED) with G_{p,a+1} once; G_pa_terms
preserves r^p and folds only the a-shift into kinds; _KIND_CHAINS extended to p,i,j ≤ 3 (63 kinds)
covers a+1 ≤ 4 (max letters (3,3), 196 entries at a=4); ϑ-invariant tp = 2 + (#u+#v letters) holds
for all entries a=1..4; no c-row re-read at the shifted index; no truncation (dxdy_terms(4) = 196
entries all resolvable). MASSED: I-extension complete.
(b) operator closed-form: a=4 dxdy_eval vs FD-of-a=3: 2.7e-8 ✓ (extends the 6e-16/4e-10/1e-9 chain);
weight-moment sanity: E[ψ3'²] = 3 exact (numeric 3.0000000000); ψ-normalization locked to the
grid-verified convention (raw-He3 decisively excluded earlier, B3 ~ 3000).
PAPER-READ (Main's pre-(c) requirement): the states ONLY the aggregate (1/π)∫S ≤ 126.80385221 at
line 1864 — NO leg-level values. Defensible finding (narrow form): "the certificate is not
reconstructible at leg level from lines 1827-1854 plus the printed c-table" — NOT "their Φ3 differs".
Both audits pass ⟹ my ‖∂x∂yΦ3‖² chain is consistent at every reachable level; the leg2 excess
(86.05 vs 27.14) remains the open 3.17× decidable by full-range Parseval. SESSION FREEZE.


### §17: 8-cell leg2 convergence table (2026-08-30, final) — answers Main exactly
avg(leg2): N=70: 88.244 (13/16 identical); N=100: 86.429; N=140: 86.046; N=200: 86.002.
Truncation |x|≤13 vs 16: IDENTICAL at every N (weights e^{−x²/2} < 1e-37 beyond 13) — truncation
error ZERO. N-drift: 2.5% monotone 70→200; converged ≈ 86.00 ± 0.06. The 3.17× leg2 excess is NOT
a quadrature artifact — outcome (A) eliminated; (A)/(B) decision now REQUIRES full-range Parseval.
avg(leg2)=86.0457 record ambiguity resolved: N=140/trim14 (both A-style; the hermgauss weight decay
makes 13-vs-16 identical anyway).
