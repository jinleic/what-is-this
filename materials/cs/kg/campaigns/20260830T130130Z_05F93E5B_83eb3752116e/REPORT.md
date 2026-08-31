# Gate C campaign — septic extension (a) + degree-5 affine analogue (b)

Agent `KgGateC-2`, 2026-08-30. Campaign frozen; nothing in this directory is edited after
writing. Pre-statement (with Addenda 1-4, each committed BEFORE the runs it governs) is
copied here as `pre_statement.md`.

Toolchain: `cs/.venv` Python 3.14.3, python-flint 0.9.0 (Arb), numpy 2.5.2. Working
precision 256 bits for the coefficient chain, 300 bits for the b5 arithmetic, 400 bits for
the K_G restatements. All parameters are exact decimals -> FMPQ; every interval is an Arb
ball built in midpoint/radius form (never `arb(lo,hi)`).

---

## 0. ANCHOR (mandatory, passed before any extension)

`code/septic_legs.py` — the septic S(theta) evaluator, specialised to s7 = 0 at the paper's
decimals (s3, s5) = (0.34101124, 0.05276111), vartheta = 0.128957369921412:

| quantity | reproduced | target (frozen corrected run) | rel |
|---|---|---|---|
| avg(leg1) | 1.9391561733360843 | 1.939156 | 8.94e-08 |
| avg(leg2) | 30.359613148308657 | 30.359613 | 4.89e-09 |
| **avg(S)** | **123.37760876657072** | **123.377609** | **1.89e-09** |
| **B3** | **14.245983003864794** | **14.2459830** | **2.71e-10** |

Tripwire residuals (measured, not asserted-as-passed):
* **T1 pi-periodicity** S(0) vs S(pi): rel **9.93e-16** (S(0) = 171.6916408202522,
  S(pi) = 171.69164082025236).
* **T2 reflection** S(theta) = S(pi - theta), 24 independent node pairs on the 49-node
  grid: worst rel **7.94e-15**.
* **T3** min S attained exactly at theta = pi/2: **1.266545**; max S = 310.80695994744156.
Both T1 and T2 match the frozen cubic-quintic reference (9.93e-16 / 7.936e-15) digit for
digit, so the septic evaluator inherits the corrected Euler-rho_j structure.

Certified head-side anchor (`code/gate_c_fast.py`, prec 256, from the certified A-grid):
* b1 = **0.8815738220495995485818877** +/- 3.58e-26 (paper Prop 6.3: b1 >= 0.881573822049)
* head = sum_{3<=m<=251 odd}|b_m| = **1.132885992768967e-5** +/- 1.92e-23
  (paper Prop 6.3: <= 1.1328860e-5; Gate A: 1.13288599276897e-5)
* gamma_head = b1 - head = **0.8815624931896718589117068** +/- 2.49e-26

---

## 1. SEPTIC DERIVATION (verified against the paper, not assumed)

Paper lines (scratch/paper_full.txt, raw dump):
* **668**: `sigma_d = (-1)^{(d-1)/2},  V_D = 1 + sum_{d in D} s_d^2`
* **672 (eq 6)**: `rho_D(t) = (t + sum_{d in D} sigma_d s_d^2 t^d)/V_D`
* **666**: `D subseteq {3,5,7,...}` — so D = {3,5,7} is inside the paper's own family
* **687 (eq 8)** / **1721**: `rho(t) = (t - s3^2 t^3 + s5^2 t^5)/(1+s3^2+s5^2)` — the
  D = {3,5} instance, matching (6) with sigma_3 = -1, sigma_5 = +1
* **857**: `D := t d/dt` (Euler); **858**: "multiplies the coefficient of t^m by m"
* **1824**: `rho_j = D_t^j rho`; **1816**: `Phi_{k+1} = D_t Phi_k - t d_x d_y Phi_k`
* **1836-1852**: the nine printed c_{p,a} rows in terms of rho_1, rho_2, rho_3 and t

Hence, with **sigma_7 = (-1)^3 = -1** read off line 668 (NOT assumed):

    rho(t)   = (t - s3^2 t^3 + s5^2 t^5 - s7^2 t^7) / V,   V = 1 + s3^2 + s5^2 + s7^2
    rho_j(t) = (t - 3^j s3^2 t^3 + 5^j s5^2 t^5 - 7^j s7^2 t^7) / V

The c-table rows are unchanged in form (they are universal in rho_j and t); the
G_{p,a} kernel machinery depends on the parameters only through r = rho(t) and vartheta.
The oddness of every exponent survives (1,3,5,7 all odd), so rho_j(-t) = -rho_j(t) and both
tripwires remain valid identities for the septic family — that is why T1/T2 are asserted in
the septic sweep, not only at s7 = 0.

Convention (Addendum 3, recorded before the dependent runs): **vartheta is held FIXED** at
vartheta_paper = 0.128957369921412, equivalently eta(s) = vartheta_paper*sqrt(V(s)). This is
required for the imported certified A-grid to be the right input (A_{a,b} = E[psi_b(X)
q_a(vartheta psi_3(X))], paper 1717-1719, depends on vartheta), and is legitimate because
eta is free in the paper's construction (line 675 "for any eta in R"). **The vartheta axis
is therefore NOT swept** — Rule-7 scope statement, see section 5.

---

## 2. (a) UPPER SIDE — VERDICT: CERTIFIED CAP, no septic improvement

### 2.1 Certified branch-and-bound over the pre-registered box

Bound (needs no tail input; valid because tail >= 0 and every |b_m| >= 0):

    for all s in B:  gamma*(s) = b1(V(s)) - head(s) - tail(s) <= b1(V_min(B)) - head_lo(B)

with b1(V) = (pi/2)(A_{1,0}^2/V - A_{0,1}^2) strictly decreasing in V (certified
A_{1,0}^2 = 0.628071295048328, A_{0,1}^2 = 1.54451017356478e-5) and head_lo(B) the sum of
certified inf|b_m| over the box for 3 <= m <= 23 (a certified LOWER bound on the full head
since the discarded terms are non-negative).

| run | ref | evals | boxes pruned | box volume certified-pruned | survivor hull |
|---|---|---|---|---|---|
| depth 7 (pre-registered) | gamma_gateA | 2,409 | 1,649 | **99.97811%** | s3 [0.32344,0.34219], s5 [0,0.07031], s7 [0,0.05625] |
| depth 7 | gamma_qhead | 2,409 | 1,650 | 99.97816% | same hull |
| depth 10 (supplementary) | gamma_gateA | 59,233 | 38,290 | **99.998739%** | s3 [0.332227,0.341602], s5 [0,0.055078], s7 [0,0.019922] |
| depth 10 | gamma_qhead | 59,057 | 38,273 | 99.998752% | same hull |
| depth 12 / 14 | both | 200,000 (budget hit) | ~143k | 49.98% / 49.85% | budget-capped, depth-first order |

The depth-12/14 full-box runs exhaust the pre-registered 200,000-box budget mid-tree
(depth-first, 8-way splits) and therefore prune *less* volume than depth 10; reported as the
budget-capped outcome it is, not hidden.

### 2.2 Residual-slab refinement (Addendum 4; same bound, survivors as root box)

Root box R = [0.325,0.350] x [0,0.080] x [0,0.020], which **contains the depth-10 survivor
hull** (verified: `logs/gc_bnb_d10_hull.json`, `contained_in_R: true`), so the two statements
compose into a statement about the whole pre-registered box.

Anisotropic bisection (split the axis maximising hi^2 - lo^2, i.e. the axis dominating the
b1 slack), FIFO order, 200,000 boxes:

| ref | pruned volume of R | survivors | survivor volume | survivor hull |
|---|---|---|---|---|
| gamma_gateA = 0.881557917504162 | **99.9999%** | 69,169 cells | 3.088e-11 | s3 [0.340803,0.341016], s5 [0.052061,0.052808], **s7 [0, 0.002344]** |
| gamma_qhead = 0.881562493211119 | **100.0000%** | 40,713 cells | 1.232e-13 | s3 [0.340987,0.341012], s5 [0.052681,0.052766], **s7 [0, 0.000781]** |

**Composed certified statement (MACHINE-VERIFIED, prec 256, exact-rational parameters,
A-grid CITED-DEPENDENCY):** for every (s3,s5,s7) in [0,0.6]^3 at vartheta = vartheta_paper,
except a residual set of volume <= 3.088e-11 (1.43e-8 of the box) contained in
s3 in [0.340803,0.341016], s5 in [0.052061,0.052808], **s7 <= 0.002344**,

    gamma*(s3,s5,s7)  <=  gamma_gateA = 0.881557917504162,

i.e. no such scheme improves on the repo's Gate A certificate. Against the sharper
reference gamma_qhead the residual shrinks to s7 <= 0.000781 and volume 1.232e-13.

### 2.3 Certified point verdicts inside the residual (where the box bound abstains)

`code/gate_c_certify.py`, prec 256. The septic-natural candidate is the **triple
annihilation** b3 = b5 = b7 = 0 (three parameters, three conditions — the exact analogue of
the paper's two-parameter b3 = b5 = 0), found by Newton at
s = (0.34101121, 0.05276111, **0.00036449**) — which lies INSIDE the residual slab:

| quantity | polished quintic (s7 = 0) | septic triple annihilation |
|---|---|---|
| V | 1.1190724020120809760 | 1.1190725129270563000 |
| b1 | 0.8815738208857831710 | 0.8815737335076835567 |
| head | 1.13298064851813860e-5 | 1.12762110436504438e-5 |
| gamma_head | 0.8815624910792979896 | 0.8815624572966399063 |

* **DELTA(gamma_head) = -3.37826580833e-8 +/- 3.99e-20 — CERTIFIED NEGATIVE.**
* Decomposition: head gain = +5.35954415309e-8, b1 cost = +8.73780996143e-8,
  net = -3.37826580833e-8. The septic direction pays more in b1 than it recovers in head.
* Certified sensitivities: **db1/dV = -0.787793605140965** +/- 5.97e-17;
  b7(quintic) = 1.17615565426e-7, annihilated at u7 = s7^2 = 1.32852960100e-7, i.e. an
  |b7| gain rate of **0.885306321646 per unit u7**. The gain rate *exceeds* the cost rate
  0.7878 — but it saturates when b7 crosses zero, and the m >= 9 coefficients grow, so only
  5.36e-8 of the available 1.176e-7 survives: less than the linear 8.74e-8 cost. That is
  the mechanism of the no-go, quantified.

Certified continuation table on the paper's own manifold b3 = b5 = 0 (Newton in (s3,s5) at
each fixed s7); every node reported, including misses:

| s7 | gamma_head (certified) | delta vs quintic | head | beats gamma_gateA |
|---|---|---|---|---|
| 0 | 0.88156249107929798962 | 0 | 1.132980649e-5 | yes |
| 0.0002 | 0.88156248264862742758 | -8.430670562e-9 | 1.131234773e-5 | yes |
| 0.0005 | 0.88156222944241032589 | -2.616368877e-7 | 1.143426688e-5 | yes |
| 0.001 | 0.88156073809921689517 | -1.752980081e-6 | 1.245457964e-5 | yes |
| 0.002 | 0.88155477430218404691 | -7.716777114e-6 | 1.654675063e-5 | no |
| 0.003 | 0.88154483522661939692 | -1.765585268e-5 | 2.335821967e-5 | no |
| 0.005 | 0.88151302705209462958 | -4.946402720e-5 | 4.516615860e-5 | no |
| 0.01 | 0.88136310052048284231 | -1.993905588e-4 | 1.482122640e-4 | no |
| 0.02 | 0.88075742577032384299 | -8.050653090e-4 | 5.664104871e-4 | no |
| 0.05 | 0.87651867877491502025 | -5.043812304e-3 | 3.494359608e-3 | no |

Every delta is certified negative. Free Nelder-Mead on gamma_head (`code/gate_c_valley.py`,
250 steps, two independent starts: the paper point and the triple-annihilation point)
converged to **s7 = 0** from both — best float value 0.881562493211119 at
(0.341011244, 0.052761125, 0.000000000).

**RETRACTION (Rule 5, recorded inline):** the float Nelder-Mead value above appeared to beat
the paper's decimals by +2.14e-11. Certified re-evaluation at the rounded optimiser digits
(0.34101124, 0.052761124, 0) gives gamma_head = 0.8815624910792979896, which is
**2.11e-9 BELOW** the paper decimals' certified 0.8815624931896718589. The apparent
improvement does not survive certification and is withdrawn; the best certified point in the
entire sweep remains the paper's own (s3,s5) at s7 = 0.

### 2.4 Tail-side arithmetic (mixed labels, stated explicitly)

At the certified quintic head value, with the tail from Lemma 6.2, tail = B3/sqrt(10*251^5):

| B3 input | label | tail | gamma* | K_G <= | improvement over Krivine |
|---|---|---|---|---|---|
| 14.44243664663976457 (paper's certificate) | CITED-DEPENDENCY | 4.57568550979349e-6 | 0.88155791539378819614 | 1.781841328136936456209 | 3.72650054433e-4 |
| 14.245983003864794 (our corrected 49-pt sweep) | COMPUTATIONAL-EVIDENCE (float quadrature) | 4.51344462145968e-6 | 0.88155797763467652994 | 1.781841202333086968398 | 3.72775858282e-4 |

The second row is **not** a certified improvement: its B3 comes from float quadrature, so it
carries COMPUTATIONAL-EVIDENCE and nothing in the repo's headline moves on it. It is
recorded because it is the honest consequence of our own corrected B3 if that B3 were ever
certified.

---

## 3. (b) LOWER SIDE — degree-5 analogue: reduction + certified ceiling, program OPEN

`code/gate_c_b5.py`, prec 300.

### 3.1 Reduction (derived; the lambda = 0 case reproduces the paper exactly)

For lambda >= 0 define the affine functional A_lambda(f,g) := 2 b1 - b3 - lambda b5, i.e.
the strip b3 + lambda b5 >= 2 b1 - C_lambda with C_lambda := sup over odd sign pairs of
A_lambda. Affine in the coefficient vector, hence stable under mixtures and coefficientwise
limits (the property the Naor-Regev transfer of paper Section 11 requires).

Reversion of H gives a1 = 1/b1, a3 = -b3/b1^4, **a5 = (3 b3^2 - b1 b5)/b1^7** (paper line
1276 for a1, a3; a5 by the same coefficient comparison, cross-checked against Heilman
2606.00247 section 15.5, which prints a5 = (3 beta_3^2 - beta_1 beta_5)/beta_1^7).
The adversary minimising M_H(c) sets b3 = max(0, 2b1 - 11/6) and b5 = max(0, (2b1 -
C_lambda)/lambda); with b3 = b5 = 0 one has M_H(c) >= c/b1 > 1 only for b1 < c, so the
largest b1 with both coefficients annihilated is min(11/12, C_lambda/2). Hence

    gamma <= min(11/12, C_lambda/2),
    the b5 route tightens the barrier  <=>  C_lambda < 11/6,  and then K_G >= pi/C_lambda.

At lambda = 0 this is exactly the paper's Prop 10.2 (C_0 = 11/6 -> gamma <= 11/12 ->
K_G >= 6pi/11 = 1.713595992867160).

### 3.2 Certified anchors

The hyperplane pair f = g = sgn(X_1), computed from exact Gaussian tail moments in Arb:

    b1 = 1,   b3 = 1/6,   b5 = 3/40,   A_0 = 2b1 - b3 = 11/6

so the paper's constant is **sharp** at lambda = 0 (paper lines 969-973 say exactly this),
and **A_lambda(hyperplane) = 11/6 - (3/40) lambda < 11/6 for every lambda > 0**: the
hyperplane stops saturating as soon as the b5 term is switched on.

### 3.3 Certified ceiling on the whole b5 program

Any admissible gamma obeys the barrier, and Naor-Regev (paper Lemma 11.2) supplies
admissible gammas converging to pi/(2 K_G); hence for **every** lambda >= 0

    C_lambda >= pi / K_G >= pi / K_G^upper = 1.763115834999343717805  (+/- 3.18e-22)

using our own certified K_G^upper = 1.7818413238804372880. Combined with
A_lambda(hyperplane) <= C_lambda this gives the **certified critical value**

    lambda* = (11/6 - pi/K_G^upper) * 40/3 = **0.93623331111986** (+/- 3.52e-15),

with the consequence table (certified):

| lambda | hyperplane constant 11/6 - 3L/40 | implied gamma bound | implied K_G >= | consistent with our certified upper bound |
|---|---|---|---|---|
| 0 | 1.833333333333333 | 0.916666666666667 | 1.713595992867160 | yes (this is the paper) |
| 0.5 | 1.795833333333333 | 0.897916666666667 | 1.749378739817982 | yes |
| 0.9347 | 1.763230833333333 | 0.881615416666667 | 1.781725111765831 | yes |
| 1.0 | 1.758333333333333 | 0.879166666666667 | 1.786687765074764 | **no** |
| 2.0 | 1.683333333333333 | 0.841666666666667 | 1.866292665498887 | **no** |

So: **for lambda > lambda* = 0.9362, sign pairs with A_lambda > 11/6 - (3/40)lambda MUST
exist** — a hyperplane-attained b5-strip is impossible there, certified, without any new
variational work. And at lambda = lambda* the barrier would give exactly
gamma <= 0.881557917499672, K_G >= 1.781841323880437 — numerically our own upper bound.
**The degree-5 affine program is therefore not a small refinement: at its ceiling it closes
the entire remaining K_G gap, and is exactly as hard as determining K_G.**

### 3.4 Swept family evidence (exact swept set stated)

Two explicit odd sign-pair families, all quantities in closed form (exact Gaussian tail
moments) in Arb; every member is a genuine pair, so each value is a valid LOWER bound on
C_lambda:
* **F1** (1-D strip): h = sgn(S)1_{|S|>=c}, k = sgn(S)1_{|S|<c}; c in {0.02, 0.05, 0.1,
  0.2, 0.25573, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5} (11 members).
* **F3** (2-D Davie-Reeds shape): h = sgn(S)1_{|S|>=c}, k = sgn(T)1_{|S|<c}; c in {0.02,
  0.05, 0.1, 0.15, 0.2, 0.25573, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0} (13 members).
Full tables in `logs/gc_b5.json`. Findings:
* b5 changes sign inside both families (F1 at c ~ 0.5, F3 at c ~ 0.4).
* Best A_0 among members with b5 <= 0: **A_0 = 1.47033139820178559** (F3, c = 0.4,
  b5 = -0.0133514242499), i.e. 11/6 - A_0 = 0.363001935132 — far from saturating.
* Swept-family lambda_max = min over b5<=0 members of (11/6 - A_0)/|b5| =
  **9.14497596448671** (binding member F1, c = 1.0), i.e. no swept member obstructs any
  lambda <= 9.14.
* Since lambda_max(swept) = 9.14 >> lambda* = 0.936, **the swept families provably do not
  contain the degree-5 extremisers** (for lambda in (lambda*, 9.14] consistency forces a
  pair with A_lambda above the hyperplane value, and no swept member provides one). This is
  consistent with the paper's own remark at line 1013 that stronger strips require
  controlling ||P_1 k||^2 as well.

### 3.5 (b) verdict

* DERIVED + MACHINE-VERIFIED: the reduction, the sharpness of 11/6 at lambda = 0, the
  hyperplane's (1, 1/6, 3/40), the certified ceiling lambda* = 0.93623331111986, and the
  structural no-go for hyperplane-attained strips at lambda > lambda*.
* **NOT settled:** the value of C_lambda for 0 < lambda <= lambda*, which is what would
  actually move gamma <= 11/12. Determining it requires the degree-5 fiber inequality — the
  analogue of the paper's Theorem 9.3 (V(u) <= 3 nu beta - beta^2) with
  V_5(u) = sum_{j<=5}(E[u psi_j])^2 over the threshold family u = sgn(p_a)1_{|p_a| >= tau|z|}
  with p_a in span{psi_0..psi_5} — a 6-parameter certified cover, versus the 2-parameter
  cover of Gate B, which is itself only PARTIAL at budget in this repo. Priced, not run.
* Duplication check (both papers read in full this session): Heilman arXiv:2606.00247 has no
  affine coefficient inequality at all — his "degree five" is a threshold polynomial inside
  f, not a coefficient constraint in the inverse-majorant chain; Jones-Malavolta
  arXiv:2603.30039 perturbs the Davie-Reeds *operator* by -eps*Pi_3 on the SDP side and
  never writes an admissible-gamma statement. C(b) is adjacent to both and duplicated by
  neither.

---

## 4. What the anchor/tripwires bought

The reflection tripwire T2 (24 independent equalities) and T1 were asserted in every septic
sweep, not just at the anchor. Residuals stayed at 7.94e-15 / 9.93e-16 — the same order as
the frozen cubic-quintic reference — which is the evidence that the septic rho_j retains the
Euler-operator parity structure. A defect of the kind that produced the retracted
avg(S) = 379.72 era would have fired T2 immediately.

## 5. Exact scope (Rule 7)

* Swept set for (a): {(s3,s5,s7) in [0,0.6]^3} at **fixed vartheta = 0.128957369921412**,
  D = {3,5,7}, sigma_7 = -1 from paper line 668. Depth 7 (pre-registered) and 10
  (supplementary) on the full box; anisotropic refinement inside
  R = [0.325,0.350]x[0,0.080]x[0,0.020] which contains the depth-10 survivors.
* **NOT swept:** the vartheta axis (needs certified A-grid regeneration at each vartheta:
  16,002 parity-allowed 1-D Gaussian integrals, priced at <= 8 wall-hours, not run); the
  "second Hermite mixing parameter" variant (sign-function perturbations beyond psi_3 —
  that is Heilman's d = 5 territory, uncertified there too, and it changes the A-grid);
  D beyond {3,5,7}; negative s_d (equivalent by s_d^2); the residual cells listed in 2.2.
* For (b): the two named families F1/F3 only, 24 members total. No claim about the full
  supremum C_lambda.
* No paper-side disagreement anywhere in this campaign: every certified quantity that the
  paper also states (b1, head, avg(S) budget, 11/6, 6pi/11) agrees with it.
