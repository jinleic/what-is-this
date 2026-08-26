# H10 over Q — session results (2026-08-24, rev 21)

## Deliverables
- `NOTES.md` — problem map with exact theorem statements and arXiv ids.
- `THEOREMS.md` — precise frontier statements through L34: published and
  audited anchors; the audited six-unknown L6 architecture; the unchanged
  L19–L22 proof that classical Schinzel H implies the conditional record;
  L23's unconditional half-sieve and exact analytic/algebraic barriers;
  L24's theorem that both unscaled self-coupled diagonal images are empty
  on $\Phi$ over $\mathbb Q_2$; L25's scaled dyadic escape; L26's
  reciprocal quartic-trace tie; L27's one-piece shear with complete
  dyadic/aligned-target local points and three section no-gos; L28's
  uniform trace-field and bad-sign-cover irreducibility; L29's complete
  reciprocal-lift squareclass classification and fixed-field reduction;
  L30's lower-degree constant-two quartic cover, square-branch rigidity,
  and exact trace-base/lift separation; L31's first proved fixed-field
  reciprocal member, exact off-cube quartic point, dispersion correction,
  and AP1 equivalence; L32's automatic-\(\Phi\) all-target maps plus
  two section-menu obstructions; L33's positive unweighted two-large
  root-pair asymptotic; and L34's sharp relative pair-moment criterion
  for AP1.  Schinzel H remains the sole conjectural input
  (`THEOREMS.md`, L6 and L19–L34).
- `h10q.py` — stdlib-only executable model of both pillars, running on a
  **proven-primality arithmetic engine**: deterministic Miller–Rabin below the exact A014233
  13-base bound, generalized Pocklington $n{-}1$ certificates above it, explicit
  `PrimalityBound`/`FactorBudget` refusal otherwise — no probabilistic acceptance anywhere.
  Its `_L6_WITNESSES` and `_L9_STEERED` tables are the **evidence
  authorities**; every run regenerates each authority's canonical
  serialization and requires `data/l6_witnesses.jsonl` and
  `data/l9_steered.jsonl` to match **byte-for-byte** (regenerate via
  `--export-evidence`), failing on drift or absence.
  `python3 h10q.py` reproduces every in-suite claim below (~25 s);
  `--extended` additionally re-verifies both full 61-prime witness tables and all
  6,494 L5 pairs.
- `l6_search.py` — standalone tied-certificate searcher. `--canonical` enforces
  $\chi_w(A)=-1$ with $A$ a $w$-unit before the $b$ loop, tries only W0's branch,
  and refuses uncertified factorization per candidate.
- `l9_steer.py` — reciprocity-steered assembly searcher (certificate logic
  lives in `h10q._l9_steered_solvable`; the tool only supplies pools,
  budgets, and a deep-rescue mode). Run artifacts under `data/l9_*.jsonl`
  are search-side provenance; the frozen table is the citable object.


## Engine hardening (advisory 2026-08-12, upheld and repaid)
- **Bug class confirmed live.** The prior `factorint` trusted 12-base Miller–Rabin as exact.
  Two defects: (a) the stated 3.3e24 bound belongs to 13 bases (2..41) — the 12-base bound is
  $\psi_{12} = 318665857834031151167461 \approx 3.2\cdot10^{23}$ (A014233); (b) nothing
  enforced any bound. During the hardened re-run the guard **fired in practice**: the $w=61$
  witness search hit an MR-passing 31-digit cofactor ($\approx 3.75\cdot10^{30}$) that the old
  engine would have silently certified.
- **Consequence for prior claims.** The earlier interactive "witnesses at every $w < 120$" run
  used the unbounded engine; it was retracted and rebuilt from scratch under the proven engine
  (outcome below — the rebuilt result is stronger, $w < 300$).
- **Regression tests** (`_verify_factoring`, in-suite): $\psi_{12}$ must split (it does:
  $399165290221 \cdot 798330580441$), $\psi_{13}$ must never be certified prime, $2^{89}-1$
  (prime, above the MR range) must be *proven* prime via Pocklington, classical
  spsp/Carmichael/perfect-power/semiprime cases, and 300 random round-trips with independent
  trial-division certification of every returned factor.
- **Independent review pass** (subagent, 23 min): verdict "correct", no soundness defects; three
  findings applied — (i) the refuse-never-guess paths are now exercised in-suite (a
  Pocklington-certified prime $2^{83}3^{10}{+}1$ with 3-smooth $n{-}1$; a rho-hard-$p{-}1$
  prime $2^{90}{+}4449$ that must refuse, with `factorint(3p)` propagating the refusal; a
  balanced 122-bit semiprime raising `FactorBudget`); (ii) the `represents_q` place-restriction
  argument corrected for arity $\le 2$ (even valuations + positivity $\Rightarrow$ $-\det$ a
  perfect square); (iii) `_mr_screen` guarded against $n < 2$ (latent non-termination).

## What the artifact verifies, with exact scopes
All rows asserted on every run (`_selftest` + `_verify_results` + the audit sections); positive
membership checks are sampled or windowed, never claimed beyond their stated scope.

| Claim | Source | Scope actually checked |
|---|---|---|
| Hilbert symbol: symmetry, bimultiplicativity, square-class invariance | classical | 200 random triples × 6 places |
| Product formula $\prod_v (a,b)_v = 1$; even ramification | classical | 300 random pairs |
| Serre isotropy / Hasse–Minkowski representation | classical | fixed families + 8 sum-of-squares cases |
| Zero-representation = isotropy convention | classical | 5 definite/indefinite regression forms |
| Proven factoring: bound, Pocklington, adversarial cases | A014233, Crandall–Pomerance 4.1.3 | $\psi_{12}$, $\psi_{13}$, $2^{89}{-}1$, 11 spsp/Carmichael, powers, 300 fuzz |
| Poonen $\forall\exists$ definition of $\mathbb{Z}$ in $\mathbb{Q}$ | math/0703907 | certificates on 4 integers + 4 non-integers; $\forall^2$ coverage all $p < 50$ |
| L2′ structure theorem: exact APs, congruence description, emptiness | proved in THEOREMS.md | primes $q < 50$, two Weierstrass models (minimal + 5-rescaled), $n \le 41$ |
| Sun identities (2.2)–(2.4) | 2607.28606 | 150 random exact instances |
| Sun block soundness (Prop 2.1 route) | 2607.28606 | ≈3000 instances, 0 violations |
| Sun Lemma 2.2 | 2607.28606 | exhaustive over $\mathbb{F}_q$, $q \in \{2\} \cup$ odd $q < 50$ |
| **L4: Lemma 5.1 for small fields** | **proved here** | **exhaustive, all odd prime powers $q \le 25$ (+ sanity $q \le 49$)** |
| **L5: Lemma 5.1 = Hasse for Legendre family** | **proved here (full proof + machine check)** | **identity, Hasse, $4 \mid \#E$, steps (iii)–(vi): all odd prime powers $q \le 49$ in-suite, $q \le 250$ extended (6494 pairs)** |
| Bridge certificates at input $z=w$ (bridge evaluated at $z^3=w^3$), $\Delta = \{2,w\}$ | 2607.28606 completeness, corroboration | every odd prime $w < 100$ in-suite; $w < 300$ (61 targets) via `--extended` |
| Ternary ($s{=}0$) restriction: sound but incomplete | proved here (certificate) | 224 integral instances, 24 completeness losses |
| **L6/W0–W2: canonical selector, tied local criterion, exact norm obstruction** | **proved here** | **character identity on all odd prime powers $q\le49$; selector all odd primes $w<300$; 2000 exact local instances (0 mismatches); 200 exact norm identities; $M_\tau{=}0$ zero case** |
| **L6/E6: canonical tied global certificates (bounded evidence)** | **new architecture; assembly still open** | **all 24 odd primes $w<100$ in-suite; all 61 odd primes $w<300$ extended; $A$ nonsquare $w$-unit, $\Delta=\{2,w\}$ exactly; 32 $\tau=0$, 29 $\tau=2a/A$; export `data/l6_witnesses.jsonl` byte-identical to the authority's serialization on every run** |
| Equidistribution / apparition (pillar 2) | math/0306277 | 35 primes < 200; bins for $n \le 4000$, primes < 40000 |
| **L8a: absorption (even valuations ⟹ local solubility at odd $v$)** | **proved here** | **150 exact local instances in-suite (400 extended), plus the odd-valuation counterexample re-checked** |
| **L8b: 2-adic parity wall $v_2(-2AbP)=1$** | **proved here** | **600 random admissible witnesses in-suite (2000 extended); exact valuation arithmetic** |
| **L8c: exact-matching gauge is empty** | **proved here (corollary of L8b + $s\ne0$ Pell parametrization; $s=0$ slice forced $v_2(b)$ odd)** | **≈2300 gauge points (both Pell roots) in-suite, none 2-adically admissible** |
| **L8d: rescues are alignment events** | **evidence, scoped** | **10/10 frozen rescues: wild odd-mult primes $1.6\times10^5$–$4.5\times10^{21}$, all symbols $+1$** |
| **L9: steering lemma (one unchecked place is forced)** | **proved here** | **200 random product-over-$T$ instances in-suite (600 extended)** |
| **L9a: prime-$b$ symbol $(x,d)_{q_1}=(A\vert q_1)$, $v_{q_1}(x)=-4$, $A\equiv5\bmod8$** | **proved here (unit hypotheses factor by factor)** | **100 exact instances in-suite (300 extended); the refuted product-form guard and the refuted class-modulus recipe are both re-refuted every run** |
| **L9-T: steered assembly witnesses** | **bounded evidence; assembly still open** | **345/353 cells (odd prime $w<100$; fixed 15-value $u$-pool filtered by $v_w(u)\ge0$ — 11/14/13 retained at $w=3/5/7$, 15 elsewhere) frozen in `_L9_STEERED`; full table replayed every run (steered forcing genuinely exercised, wild symbol cross-replayed); export `data/l9_steered.jsonl` byte-identical; 8 named cells remain open (search-exhaustion, not certified failure)** |
| **L10-0: free-place lemma $v_p(d)=0$ and $v_p(x)$ even $\Rightarrow$ symbol $+1$** | **the implication is proved; the SET EQUALITY once claimed here is RETRACTED** (the frozen set is *not* $\{2,\infty\}\cup\operatorname{supp}(\alpha)\cup\{q_1\}$ — $\delta_\tau$ carries $A$, so $p\mid A$ has odd $v_p(x)$; complete set in L11f) | **$v_p(d)=0$ and $v_p(x)$ even $\Rightarrow$ symbol $+1$: 150 instances in-suite (400 extended); exact Taylor shift re-verified on 40 random polynomials $\times$ 4 offsets** |
| **L10a: $\alpha=-1$ exactly on the canonical branch $\tau=2a/A$** | **identity proved; its "no place dividing $A$" corollary RETRACTED** (L11f) | **100 random admissible $s$ in-suite (300 extended); identity $\delta_\tau=1/A$ exact** |
| **L10b: $\tau=0$ wall, $\prod_{p\mid A}(x,d)_p\cdot(A\vert q_1)=(-2\varepsilon\vert A)=-1$** | **proved here for $A$ prime integral with $a$ odd (admissible), $A\nmid2zD_z$; evidence otherwise** | **80 in-suite instances (250 extended) with hypotheses asserted; 200 in-suite Jacobi-form instances per run (600 extended) for rational/non-squarefree $A$ (NOT proved); the non-admissible case $a=2$, $A=17$ is frozen as a regression showing the wall fails there; no quantification over the infinite $s$-range is claimed** |
| **L10c: aligned classes exhibited** | **SUPERSEDED — 130 of 164 rows invalid** | the certificate omitted the controlled places $p\mid A$; 130 rows have symbol $-1$ there, the other 34 carry stale moduli under L11f. Table retired; 8 rows of its shape frozen as regressions that must stay refused |
| **L11a: branch completion is free** | **proved here** (soundness is $\tau$-uniform; a $k$-fold disjunction is one polynomial in the same $(y,r)$) | **120 exact instances in-suite (400 extended)**; hypothesis: denominator of $\tau$ proved nonvanishing on $\Phi$ |
| **L11b: the four square classes** | **proved here** (exact identity per branch) | $\tau\in\{0,2a/A,1,(1+2a^2)/(1+4a^2)\}\to\alpha\equiv-A,-1,A,1$ |
| **L11c: square-branch family $\tau_d$** | **proved here** | $\alpha=((A-d^2)/(2d))^2$ for every rational $d\ne0$; only finite subfamilies are adjoinable |
| **L11d: target-place character** | **proved here as a character identity only** | $\chi_w(\alpha)=+1$ for all odd $w\nmid2a$ (122 pairs in-suite); **NOT W1**, which also needs $v_w(A)=v_w(\delta)=0$ — false at $a=1,w=5$ |
| **L11e: aligned class at every cell** | **RETRACTED — the claim was FALSE** | claimed 353/353; the certificate omitted controlled places with odd $v_p(x)$. Regression frozen: cell $(3,1)$, $a=1$, $q_1=41$ has $(x,d)_5=-1$ while $2,41,\infty$ are $+1$. 293/353 rows were bad |
| **L11f: the controlled set is complete** | **proved here** | $S=\{2,3,5,7\}\cup\operatorname{supp}(\alpha)\cup\operatorname{supp}(\delta)$; even $p$-content of $P$ + $\deg P=8<p-1$ for $p\ge11$; the $p=3$ fixed divisor exhibited (60 instances in-suite, 200 extended) |
| **L11g: the wall is branch-INDEPENDENT** | **proved here** ($A$ prime, admissible $a$, $v_A(z)\le0$) | $(x,d)_A(x,d)_q=-1$ on **608 (cell, branch, $q$) instances in-suite, 2136 extended, across all four branches**; subsumes L10b as its $\tau=0$ case |
| **L11h: characterization of reachable cells** | **sufficiency proved (constructive); necessity proved for $A$ prime, EVIDENCE for composite/rational $A$** | reachable **iff $z$ has a numerator prime $\equiv1\bmod4$**: **190 of 353** rows replayed with escape prime and constructed $a$; for the other **163**, necessity is a theorem for all admissible $a$ via L13e (was: L11g for prime $A$ + recorded 878,400-certificate search) |
| **L12a: the wall generalizes to every coprime $b$** | **proved here** ($A$ prime, $v_A(z)\le0$, odd places of $b$ coprime to the cell data) | $(x,d)_p=(A\vert p)$ at every odd place of $b$, so $\prod_{\{A\}\cup\operatorname{oddsupp}(b)}(x,d)_p=-1$: **231 instances in-suite over prime, semiprime, three-prime, squarefull and rational $b$, all four branches**; subsumes L11g and refutes the $(m\vert A)=-1$ cofactor escape |
| **L12b: the non-coprimality door** | **constructive + class alignment VERIFIED (L13a)** | $b=\varepsilon fq$ with $f\mid z$ gives every controlled symbol $+1$ at **103 of the 163 L11h-walled cells**; minimal set $S_{\min}=\{2,3,5,7,f\}$, median modulus $2.4\times10^7$, **883 prime members sampled, 0 violations** (`l12_class.py`); emergent places stay step (ii) — 467 of 959 spot-probed carry $-1$, count always even |
| **L13a: verified class alignment of the 103 escapes** | **proved here (frozen side)** | 12 rows × 2 members re-certified in-suite per default run (40 rows extended); full 883-member sweep in `data/l12b_class_sample.json`; excluded-member set (174 primes) named per row |
| **L13b: grid witness coverage 353/353 — COMPLETE** | **proved per row, suite-asserted** | 60 walled-no-escape cells fully soluble (55 Hasse-Minkowski + 5 steered certs, all lead-replayed) + 7 frozen residual-cell witnesses ($\_L13\_RESIDUAL$); the last cell $(89,(-2,1))$ closed 2026-08-18 by L13f (witness $(3,\,89/367)$, $\tau=19/37$, tied True, lead-replayed through the five-rung cofactor ladder) |
| **L13c: no square-class shortcut exists** | **proved per row (706/706 sympy; 24/96 in-suite Frobenius mod $p<1000$)** | $P$ irreducible deg 8 over $\mathbb Q$ on every (cell, branch); content sqf $\in\{\pm1,\pm5\}$; 137.5M square-class tests + 40 designed families: **0 hits**; Galois $D_8\wr C_2$ is EVIDENCE (28,240 certified unramified Frobenius samples; $S_8$ excluded) |
| **L13d: emergent-free members exist** | **certified instances, ~25% empirical density** | **24 zero-bad members, 15 distinct cells**, 16 frozen rows replayed every run; 15/57 fully factored members in the dedicated sweep, they NEVER vanish everywhere simultaneously — step (ii) is a density statement |
| **L13e: the wall extends to composite squarefree $A$** | **proved here** (derivation audited; $A\equiv5\bmod8$ composite or rational included) | 40 in-suite instances per run across all four branches ($a\le77$ composite pools), lead-replayed; full 316-instance sweep + derivation in `data/l13_compositewall.{md,json}`; consequence: L11h necessity theorem for all admissible $a$ |
| **L22a: fixed-$Z$ factor elimination on $\tau^\dagger=(1+2a^2)/A$, $\delta=-4a^4/A$** | **PROVED for every fixed $Z\in\mathbb Q^\times$** | $H=(A/4)P$ has equal endpoints $a^8A^2Z^4$; the $A$-adic slopes force even factor degree, the $a=0$ and infinity endpoints exclude degree $2$, and degree $4$ forces $Z^4=1024(1-Z)^2$ with nonsquare discriminants $1152$ and $896$.  $Z=0$ is explicitly reducible and excluded; $Z=1$ has an exact $a=1$, mod-$11$ certificate.  Symbolic proof: `/tmp/l22_elimination.md`; replay: `l22_elimination.py`, `data/l22_elimination.jsonl`. |
| **L22b: reciprocal all-$Z$ theorem on the same named branch** | **PROVED independently** | At $a=1,A=5,s=0,\tau=3/5=\tau^\dagger,\delta=-4/5$, $P=b^4T(b+b^{-1})$; trace reducibility forces $125D^2+64Z^4$ square, whose reduced numerator is $5\bmod8$.  The lift obstruction maps to $E:y^2=x(x+5)(x+80)$; complete $2$-isogeny descent gives rank $0$ and only its four $2$-torsion points, excluding positive $x$.  Proof: `/tmp/l22_reciprocal_cube.md`; replay: `l22_reciprocal_cube.py`, `data/l22_reciprocal_cube.jsonl`. |
| **L22 infinity clusters** | **PROVED local structure; NOT the closure** | Four irreducible quadratic Laurent clusters permit only proper degrees $2,4,6$; the natural $a=0$ place usually narrows to $4+4$ but has an explicit exceptional locus.  L22a, not this local method, eliminates the survivors.  Proof: `/tmp/l22_infinity.md`; replay: `l22_infinity.py`, `data/l22_infinity.jsonl`. |
| **L22 factor-tuple route** | **PROVED criterion/no-go; CONDITIONAL prime values** | The odd-exponent factor tuple plus one all-plus residue in its resultant-character period is sufficient under classical Schinzel H.  The admissible tuple $(8t+5,8t+7,8t+23)$ has signs $(-1,-1)$ while their product is $+1$, proving reciprocity controls only the product.  This route relocates rather than removes the obligation and is not the closure.  Proof: `/tmp/l22_factor_tuple.md`; replay: `l22_factor_tuple.py`, `data/l22_factor_tuple.jsonl`. |
| **L22 free-$\lambda$ square branch** | **PROVED off-chain algebra/class theorem; NOT APPLICABLE to the six-count** | The $\lambda$-pencil is irreducible and yields aligned classes after HIT, but the selected $\lambda$ ranges through the infinite L11c family after the cell is fixed; making it variable costs an additional witness.  It must not be cited as the record closure.  Proof: `/tmp/l22_square_branch.md`; replay: `l22_square_branch.py`, `data/l22_square_branch.jsonl`; scope audit: `agent://BranchAudit`. |
| **L22 geometric Noether route** | **PROVED reduction; OPEN auxiliary nonvanishing only** | The norm pencil is generically geometrically integral and fixed-fibre reducibility is detected by the gcd of the $483$-minors of an explicit $903\times483$ matrix.  Nonvanishing at every cell parameter remains OPEN in this side route, but L22a makes it unnecessary for the chain.  Proof: `/tmp/l22_fiber_geometry.md`; replay: `l22_fiber_geometry.py`, `data/l22_fiber_geometry.jsonl`. |
| **L22d: classical Schinzel H $\Rightarrow$ intermediate H $\Rightarrow$ Theorem C** | **PROVED implication; conclusion CONDITIONAL on classical Schinzel H** | L22 plus quantitative HIT selects an irreducible $a$ in L20's odd character-admissible progression; L20 supplies the aligned $f=w$ class and conditions (b)–(d), while L22 completes condition (a); L19 supplies a member from classical Schinzel H.  All four conditions are therefore PROVED.  The finite branch menu and exact six-count are unchanged (`THEOREMS.md`, L19–L22; `CONDITIONAL.md`, §§1–2). |
| **L23a: half-dimensional bad-root sieve** | **PROVED unconditional small-prime-clean theorem; NOT a member theorem** | For a selected aligned linear/octic pair with $u(\theta)$ nonsquare, $r_-(p)$ has Chebotarev mean $1/2$, $\kappa=1/2$, and the mod-$p$ root oversieve has BV level $D=X^{1/2}/(\log X)^B$.  The semilinear lower sieve with $z=X^{49/100}$ gives $\gg X/(\log X)^{3/2}$ prime members clean at every bad root prime below $z$.  Exact odd valuation would require $p^2$ moduli; finite local replays are not the analytic proof (`l23_half_sieve.py`, `data/l23_half_sieve.jsonl`; `agent://HalfSieveAudit`). |
| **L23b: large-bad-divisor barrier** | **PROVED barrier; unconditional global member existence remains OPEN** | Since $|G(t)|=X^{8+o(1)}$, as many as $16$ factors above $X^{0.49}$ survive, and reciprocity permits $0,2,\ldots,16$ bad factors.  The first missing estimate is a saving bound for the sector $p_1,p_2\ge z$, $p_1p_2>D$, both bad and of odd valuation.  Even EH leaves eight; Kao gives only total $P_{12}$ after an AP adaptation, whereas reciprocity needs the exact threshold $R_{\rm bad}\le1$. |
| **L23c: Capell square-in-$K$ route** | **PROVED empty on the selected refined protocol; strict residual scope retained** | For irreducible $P$, $2\theta\in K^2$ iff $P(u^2/2)$ is reducible; $N(2\theta)=256$ is not sufficient.  Odd $v_w(z)$ gives an odd left-edge valuation for every admissible $a$; otherwise the free refinement $a\not\equiv1\pmod w$ gives the simple root $\beta=(2As^2)^{-1}$ with $2\beta=1/(As^2)$ nonsquare.  A pre-existing even-$v_w(z)$ class with $w\mid s$ remains OPEN uniformly, though all $353$ stored rows are proved individually.  No universal $2$-adic Capell claim is made. |
| **L23d: fixed-cell conic bundle** | **PROVED rank/Brauer calculation; no unconditional point theorem applies** | $X_{a,Z}:U^2-2bV^2=P(b)W^2$ has ten geometric degenerate fibres in closed degrees $1+8+1$, splitting classes $(A,2\theta,A)$, and $K_X^2=-2$.  On the proved nonsquare scopes its non-split rank is $10$ and $\operatorname{Br}(X)/\operatorname{Br}(\mathbb Q)=\mathbb Z/2$, generated by $(A,b)$; the hypothetical square case would have rank $2$ and trivial quotient.  HWW first fails at rank $\le2$; the non-split degree-$8$ closed point defeats the other checked low-rank hypotheses (`l23_fibration.py`, `data/l23_fibration.jsonl`; `agent://FibrationAudit`). |
| **L23e: norm and squareclass walls** | **PROVED in the stated ansatzes; residual curves OPEN** | $H=X^2+ALY^2$.  Exact matching $2b=-ALq^2$ plus $(A,L)=1$ is $\Phi$-empty by $v_2$; no polynomial-in-$b$ norm section exists because $P(0)=(2a^4Z^2)^2A$.  On $\Phi$, $v_2(P)=4+4\min(v_2(Z),0)$, so $P=\pm2b\,y^2$ is empty, while $P=A\,\mathrm{Norm}$ forces $(P,2b)_2=-1$.  The fixed-$P$ square locus is genus $1$.  The first admissible factor ansatz $L=-A\rho^2$ factors $H$ but reduces to scoped bad symbols and genus-$4$ curves, not a section (`l23_norm_section.py`, `l23_rational_section.py`, `l23_squareclass.py`; paired data). |
| **L23f: moving-$a$ complexity** | **PROVED for the tested diagonals; no universal minimum claim** | On the $a(0)=0$ branch, exactly $Q^4\Vert H$; removing it gives a squarefree degree-$18$ representative with factor degrees $2,16$.  The separate $a(0)\ne0$ sharp diagonal has raw degree $21$ and factors $2,19$.  This corrects the earlier apparent minimum-$21$ claim: degree $18$ is the smallest certified example in this audit only (`l23_multivar.py`, `data/l23_multivar.jsonl`). |
| **L24a: self-coupled diagonal covers** | **PROVED finite flat degree $8$; nominal $\le5$ count is vacuous** | The two exact systems eliminate to $Q_I=Au^2+(A^2-c^2-16AB)u+16A^2Bs^2-16A$ and $Q_{II}=Au^2-(c^2+16AB)u+16A^2B(s^2-1)-16A$, with $u$ and $A+u$ both squares.  Clearing by $E^2=(Ab^2D)^2$ introduces no inverse witness.  Each guarded coordinate algebra is free of rank $8$, so tied rank is at most $1$ before local admissibility (`l23_absorption.py`; L24 producers below). |
| **L24b: both diagonal orientations are $\mathbb Q_2$-empty** | **PROVED uniformly on all of $\Phi$; diagonal five-count CLOSED/FALSE** | $v_2(c)\ge4$.  Orientation I has normalized coefficient valuations $(4,0,0)$ and possible $v_2(u)=0,4$, incompatible with simultaneous squareness of $u$ and $A+u$.  Orientation II is $u^2-32hu+16e$ with unit $h,e$; $v_2(u)=2$, and $u=4t^2$ gives $0\equiv8+2Ab(s^2-1)\pmod{16}$, which would require the impossible $v_2(s^2-1)=2$.  Both images are empty; this does not touch the fixed canonical six-count (`l24_diagonal_geometry.py`, `l24_diagonal_arithmetic.py`; paired data; `agent://DiagonalAudit`). |
| **L24c: local geometry and bounded search** | **PROVED local/geometry statements; zero-hit scan EVIDENCE only** | Odd target places are locally soluble on orientation II's standard $k=1<v_w(Z)$ stratum, so the global obstruction is genuinely dyadic.  The generic cover is integral degree $8$ with Galois closure $V_4\wr C_2$ of order $32$, no base-rational section, genus-$35$ actual-cell slices and a genus-$53$ actual-$z$ slice; absolute surface rationality remains OPEN.  The scan covered $370$ cells, $2{,}013{,}141$ bases and $4{,}026{,}282$ orientation attempts with zero hits; every-cell emptiness comes from the mod-$16$ proof, never the scan (`l24_diagonal_{local,geometry,arithmetic,search}.py`; paired data). |
| **L25: scaled non-diagonal couplings** | **PROVED local escape; global five-shape route remains OPEN** | With $u=\rho^2$, fixed scalings give exact eliminants.  The type-I scalings $(y,r)=(2X,\rho)$ and $(2X,\rho/2)$ are $2$-adically admissible exactly on even-$s$ and odd-$s$ strata, uniformly in unit classes and every $v_2(c)\ge4$; their disjunction covers all of $\Phi$, so L24's obstruction is not universal over coupled formulas.  A freely chosen unit-$b$ unramified target stratum has an exact smooth construction at every odd $w\ge7$ (Jacobian $-2400$); the streamlined proof does not lift at $w=3,5$, and compatibility with preselected L20 target residues remains separate.  Two guarded bridge samples provide parity-wise real certificates.  No rational coupled point, five-count improvement, or member theorem is claimed (`l25_scaled_coupling.py`, `data/l25_scaled_coupling.jsonl`; `THEOREMS.md`, L25). |
| **L26: reciprocal even-pullback tie** | **PROVED structural/local theorem; uniform member theorem remains OPEN** | With $Z=z^2$ and $\eta=Z(b+1)/(Db)$, $P_{\rm rec}=b^4T(b+b^{-1})$ for an exact quartic $T$; the canonical tie is split at $2$, standard guarded W1 targets lift, and $(2b|p)=(2(u+2)|p)$.  L28 proves the trace field and bad-sign cover; L29 classifies the distinct lift; L31 proves one exact member at \(w=13\) but not a per-target theorem (`l26_reciprocal_tie.py`, `data/l26_reciprocal_tie.jsonl`; `THEOREMS.md`, L26–L31). |
| **L27: one-piece triangular shear** | **PROVED complete local theorem + three scoped section no-gos; global cover remains OPEN** | The fixed shear $(2X+28\rho,sX+\rho)$ covers every dyadic $s$-parity and every standard aligned odd target, including $3,5$; one guarded real sample works.  The linear-$B$ coefficient factorization excludes its natural cancellation section; $C_2\equiv128\bmod256$ excludes $\lambda=\pm1,\pm A$ uniformly; and $y=0$ forces a conic with no $\mathbb Q_3$-point.  None is a no-point theorem for the full degree-$8$ cover (`l27_triangular_shear.py`, `data/l27_triangular_shear.jsonl`; `THEOREMS.md`, L27). |
| **L28: uniform reciprocal trace field** | **PROVED over $\mathbb Q_2$ and $\mathbb Q$; lift question superseded by L29** | The shifted trace quartic has a one-segment Newton polygon, and its Ferrari resolvent has no square root in any of the five possible dyadic valuation cases.  Thus $T$ is irreducible uniformly.  Moreover $N(2(\theta+2))\in A\mathbb Q_2^{\times2}$ with $A\equiv5\bmod8$, so the bad-sign extension is nontrivial and $T(v^2/2-2)$ is irreducible.  L29 classifies the distinct reciprocal lift; L31 proves one exact \(w=13\) member, but no uniform member theorem follows (`l28_trace_field.py`, `data/l28_trace_field.jsonl`; `THEOREMS.md`, L28–L31). |
| **L29: reciprocal lift and dyadic compressor** | **PROVED local classification + exact reductions; uniform member/parity step remains OPEN** | For square $Z$ with even $t=v_2(Z)$, $\theta^2-4$ is square in the trace field exactly for $t=-2$ or $t\ge2$, and nonsquare for $t=0$ or $t\le-4$.  The negative cases exhaust every possible Eisenstein-quartic factor congruence; the exceptional split case uses a two-slope ordinary resultant, independently validated after replacing an under-specified higher-polygon presentation.  The pullback $Z=8z^2/(1+z^2)$ preserves every positive odd valuation and forces $t\ge2$.  The slice $b=w\rho^2$ freezes the norm field but has no generic section.  L31 promotes the \(w=13\) row to a proved global member; the remaining finite rows stay EVIDENCE and no per-target theorem follows (`l29_reciprocal_frontier.py`, `data/l29_reciprocal_frontier.jsonl`; `THEOREMS.md`, L29–L31). |
| **L30: constant-two quartic frontier** | **PROVED degree drop + dyadic/all-target selected-fibre theorem; global rational point remains OPEN** | The specialization $(y,r)=(2,sX+\rho)$ gives the exact quadratic $H_2(m)=0$ in $m=\lambda^2$, hence a genuine even quartic on every $\Phi$ base rather than L27's octic.  A normalized Hensel argument solves every dyadic stratum; a two-quadratic character count gives fibre-regular aligned rows for every odd $w\ge5$; and the exact fibre $(a,b,z)=(5,3,3)$ is solved by strong Hensel at $w=3$.  The control $(1,3,3)$ is $\mathbb Q_3$-empty, so the theorem selects $a$ rather than claiming fixed-$a$ uniformity.  One guarded real stratum works.  The producer also proves generic square-branch/shear rigidity and exact trace-base/lift separation.  The $13{,}224$-row zero-hit global scan is EVIDENCE only (`l30_quartic_frontier.py`, `data/l30_quartic_frontier.jsonl`; `THEOREMS.md`, L30). |
| **L31: one exact reciprocal member and sharpened remaining walls** | **PROVED exact member/reductions; uniform closure remains OPEN** | On the mandatory fixed-field slice, \(w=13,a=3,Z=169,q_0=1,b=13\) gives \(M=2^4Q/(3^2 37^3 83^2 1033^2)\), \(Q=14082426920623718389\) prime by a recursive Pocklington certificate, and \((M,26)_v=+1\) at every place.  The reciprocal lift is explicit: \(\rho=14/13,\lambda=12\).  Separately, \(a=1,c=-64/25,\lambda=3,b=-3253/3125,Z=2033125/6101423\) is an exact rational \(H_2\)-point before the cube condition; \(Z\) is not a cube and the shared prime \(3253\) leaves \(c\) a unit, so it is not an L30 target point.  Cube compatibility reduces to two genus-\(2\) curves.  The two-large sector is corrected to a positive-main Buchstab/Hilbert-detector problem beyond BV, not a signed-error-only estimate.  AP1 (\(R_{\rm bad}\le1\) in one selected class per cell) is equivalent to intermediate H by even Hilbert parity but remains OPEN (`l31_frontier_push.py`, `data/l31_frontier_push.jsonl`; `THEOREMS.md`, L31). |
| **L32: automatic-\(\Phi\) target maps and section-menu obstructions** | **PROVED local parameterization; global rational point remains OPEN** | The maps \(\sigma(t)=t/(1+2t^2)\), \(a(t)=1+2\sigma(t)\), and \(b_\kappa(u)=(u^2+\kappa)/(5u^2+\kappa)\), \(\kappa\in\{1,2,-2\}\), lie in the dyadic \(\Phi\)-conditions for every rational parameter.  Quartics \(F_5=20t^4+32t^3+36t^2+16t+5\) and \(F_9=36t^4+32t^3+52t^2+16t+9\), with resultant \(2^{28}\), give \(N_w\ge(w-13\sqrt w)/4-4>0\) for \(w\ge211\); exact exhaustion closes \(5\le w<211\).  The identity \((-1|w)(-2|w)(2|w)=1\) lets the \(b\)-menu achieve \(v_w(b)=1\) at every odd target.  Thus one finite rational-map menu is locally soluble at both \(2\) and every \(w\ge5\).  Conversely, \((b-1)^2/b=4ra^2\) with fixed \(r\) can meet a unit-\(a\) target only when \(v_w(r)=-1\), so every finite fixed \(r\)-menu is non-uniform; \(b=1\pm2a\) covers only \(w\equiv5,19\bmod24\).  Emergent numerator/denominator primes remain uncontrolled (`l32_local_parameter_frontier.py`, `data/l32_local_parameter_frontier.jsonl`; `THEOREMS.md`, L32). |
| **L33: unweighted two-large root-pair main** | **PROVED positive asymptotic; prime-weighted exact-odd transfer remains OPEN** | For every \(0<\alpha<\beta<1/2\), L23 Chebotarev gives \(\sum_{X^\alpha\le p\le X^\beta}r_-(p)/p=\frac12\log(\beta/\alpha)+o(1)\).  CRT supplies \(r_-(p)r_-(q)\) classes and \(X/(pq)+O(1)\) integers per class, while the total error is \(O(X^{2\beta})=o(X)\).  Hence \(\mathcal U_2(X;\alpha,\beta)\sim\frac18\log^2(\beta/\alpha)X\).  With \(\alpha=0.49\), this is a proved positive main in the first beyond-BV range.  The theorem is unweighted and uses the mod-\(p\) root oversieve; \(\Lambda(Q(t))\), exact odd valuations \(p^2q^2\), and small-prime cleanliness remain the pointwise fixed-family detector problem (`l33_unweighted_pair_main.py`, `data/l33_unweighted_pair_main.jsonl`; `THEOREMS.md`, L33). |
| **L34: sharp pair-moment threshold for AP1** | **PROVED exact sufficient criterion; prime-weighted estimate remains OPEN** | Even Hilbert parity gives \(\mathbf1_{\{R_{\rm bad}=0\}}\ge1-\binom{R_{\rm bad}}2\).  Since \(R_{\rm bad}\le R_{\rm root}\), \(\mathcal N_{\rm good}\ge\#\mathcal A_z-\mathcal P_2^{\rm root}\); exact odd valuations are unnecessary for this upper-bound route.  With \(\mu(\vartheta)=\frac12\log(8/\vartheta)\), the unrestricted pair envelope \(\mu(\vartheta)^2/2\) is below \(1\) for \(\vartheta>8e^{-2\sqrt2}\); the divisor-product constraint can only lower its region.  At \(\vartheta=0.49\) the envelope is \(0.9749604961\ldots\).  Thus \(\mathcal P_2^{\rm root}\le(0.9749604961\ldots+o(1))\#\mathcal A_z\) would imply AP1.  Current sieve constants and beyond-BV prime weighting after small-prime conditioning do not prove this safe \(2.5\%\)-margin estimate (`l34_ap1_pair_threshold.py`, `data/l34_ap1_pair_threshold.jsonl`; `THEOREMS.md`, L34). |
| **L21d/e: precursor on the constructed branch** | **PROVED generically + PROVED per cell on 353/353; SUPERSEDED AS THE FRONTIER by L22** | Generic: $H=(A/4)P$ has $L=a^8A^2Z^4$ as both end coefficients; the $(a,Z)=(1,27)$ specialization is irreducible mod $17$, proving irreducibility over $\mathbb Q(a,Z)$.  The former vertical-line repair used 353/353 exact fibre certificates plus Cohen–Serre quantitative HIT; L22 now proves every fixed nonzero rational $Z$ directly (`data/l21_irred_generic.jsonl`; `/tmp/l22_elimination.md`). |
| **L21: the reducible locus, and its avoidance** | **PROVED (both lemmas)** | (a) $s=0$ and $\delta_\tau=\sigma^2$ a square $\Rightarrow$ $P=(4DAb^2-\sigma a^2Z^2N_g)(4DAb^2+\sigma a^2Z^2N_g)$ — an explicit $4\times4$ factorization, so blanket irreducibility is **false**; 56/56 exact, 8/8 controls inapplicable. (b) On the constructed branch $\delta_{\tau^\dagger}=-4a^4/A<0$, never a square, so the degeneration **cannot occur there**, for any $a$ or cell (9/9). Square class of $\delta_{\tau^\dagger}$ is that of $-A$. (c) $P+32A^3s^2D^2b^5$ is palindromic, $P_0=P_8$ always (288/288) $\Rightarrow$ Galois group in $C_2\wr S_4$ (`data/l21_reducible_locus.jsonl`) |
| **L20: uniform class existence** | **PROVED** | For every cell take $f=w$ (available since $v_w(z)\ge1$); pick odd $a$ with $(A\mid w)=-1$.  The character sum is $-1$, with exactly $\bigl(w-\left(\frac{-1}{w}\right)\bigr)/2$ qualifying residues: $(w+1)/2$ for $w\equiv3\bmod4$, $(w-1)/2$ for $w\equiv1\bmod4$.  Thus the set is always nonempty.  The explicit residue system for $q_1$ mod $M=4A\prod_Sp$ has $\varphi(M)/2^{\#\{p\}}$ classes, and Dirichlet finishes. **Clause (ii) of H is removed for every cell.** 103/103 canonical rows reproduced (1,113,000 residues enumerated), 353/353 grid certificates clean, composite $A$ replayed, 0 refusals (`data/l19_classexist.jsonl`). Fixing $a=1$ collides exactly at $w=11,19,31,59,71,79$ — the independently derived 5-wall set |
| **L19: Schinzel H $\Rightarrow$ per-class H** | **PROVED implication (conclusion CONDITIONAL on Schinzel)** | For a fixed verified aligned class with $F(t)=c\,G(t)$ (the **square class** of $c$ is $S$-supported; $G$ is primitive with positive lead), Schinzel H for $\{q_1+Nt,G(t)\}$ gives infinitely many ladder-zero members: the $S$-strip leaves one odd place $R=G(t)$, $v_Q(P(b))=0$, all other outside symbols are $+1$, and reciprocity forces $(x_0,d_0)_R=+1$.  The recorded audit covers 293/293 canonical pairs and a 24-row prime-rung cross-check (`data/l18_schinzel_implies_h.jsonl`); L22 supplies irreducibility uniformly. |
| **L19: the canonical 5-wall and its break** | **PROVED (identity + criterion, canonical protocol); PROVED break** | $\mathrm{Hilb}_5=-(w|5)(q_1|5)$ vs wild $(5|q_1)$ $\Rightarrow$ canonical ESC 5-walled iff $(w|5)=+1$; with L11 iff $w\equiv1\bmod4$, canonical obstruction is exactly $w\equiv11,19\bmod20$. Audit: 53 primes in $[101,400]$, 26 L11 / 13 ESC / 14 obstructed, **0 mismatches** (`data/l18_fivewall.jsonl`). Break: $a{=}7,\tau{=}0,f{=}w$ closes $131,139,151$ (L12b non-coprimality door; no conflict with the L10 $\tau=0$ coprime wall) |
| **L19: the 5-wall is always breakable (class side)** | **PROVED** | Factor-by-factor: $\prod_{p\mid A}(x,d)_p(A\mid q_1)=(-2\varepsilon\mid A)(f\mid A)$, so with $A\equiv5\bmod8$ the L10b anti-correlation breaks **iff $(f\mid A)=-1$**. For the canonically obstructed primes $w\equiv11,19\bmod20$ — hence explicitly $w\equiv3\bmod4$ — the character sum gives exactly $(w+1)/2$ good residues; CRT fixes the unit hypotheses and Dirichlet supplies $q_1$. Thus every such $w$ has an aligned $\tau=0$, $f=w$ class. Audit 16,744 candidates, 908 aligned, **0 mismatches** (`data/l19_tauzero.jsonl`). Scope: class side only; member existence is **CONDITIONAL** on classical Schinzel H through L19. Fixed $a=7$ is breakable iff $(w\mid197)=-1$. |
| **L19: off-grid horizon** | **PROVED per row (lead-replayed)** | Closed $w\in\{101,\dots,179\}$ — all sixteen tried cells at $u=-1$ plus $[131,[5,1]]$; **no off-grid cell remains open**. $w=179$: $a=7$, $\varepsilon=-1$, $f=179$, $q_1=Q=251$, $k=0$, decided fraction 94.7% (10,282/10,858), 222-digit max cofactor (`data/l19_cell179.jsonl`). Divergence model: L18 tail closes 51% of the excess multiplier / 78% of the log gap vs the 15.67$\times$ Bateman–Horn ledger gap; stationary model misses the selected $k=0$ spike (`data/l18_divergence_model.jsonl`) |
| **L18: step-(ii) density law** | **PROVED (local bound; no-uniform-mass for the 293 current classes); EVIDENCE (measured exponents)** | $m_3=m_5=m_7=0$ and $m_p\le 8/(p+1)$ for $p\ge11$ (singular roots incl.); lead-verified on 102,439 local pairs, 0 violations. Chebotarev exponent $c=1/2$ in every class (293 replayed irreducibility certificates + 293 replayed nonsquare witnesses) $\Rightarrow \prod_{p\le X}(1-m_p)\asymp(\log X)^{-1/2}$: **no cutoff-uniform positive clean mass**. Measured $c$: 0.4778 / 0.4801 / 0.4590 (`data/l17_stepii.jsonl`) |
| **L18-horizon wave 2 (historical frontier, SUPERSEDED by L19)** | **PROVED per row; one then-open protocol obstruction** | The recorded closures were $w\in\{113,127,137\}$ at $u=-1$, while $w=131$ was obstructed inside that protocol after 336 attempts and 326 eligible collisions (`data/l17_horizon{113,127,131,137}.jsonl`). L19 subsequently closed $w=131$ through the non-coprime $f=w$ route (`data/l18_route131.jsonl`); it is not open now. |
| **L17-horizon: off-grid closures $w\in\{101,103,107,109\}$** | **PROVED per row (independently replayed by lead)** | Verified horizon beyond $w\le97$: $[101,[-1,1]]$ L11 $a{=}5$ member $k{=}16$ factorint rung; $[103,[-1,1]]$ ESC $f{=}103$ $k{=}0$; $[107,[-1,1]]$ ESC $f{=}107$ $k{=}0$; $[109,[-1,1]]$ L11 $a{=}71$ $k{=}0$ (auxiliary tied/ramified budget-refused — accepted L14 closure form); off-grid machinery: fresh-L11 roots + `_l10_class_cert` / temporary `_L12_ESCAPE` + L13 ladder (`data/l17_horizon{101,103,107,109}.jsonl`) |
| **L17: mechanistic p-adic sieve — hit-conditioned matched fit exact; cofactor = reciprocity-level proved; rate band-only** | **per-row PROVED; matched model validation EVIDENCE (marginal 0.9993; union z=-0.09, empirical tails 0.61/0.53 under marginal-preserving permutation)** | (i) Residue-determination at full grid: 293/293 classes, 0 counterexamples (`data/l17_sieve.jsonl`, 8,760 members, $k\le119$). (ii) Exact bad masses $m_p$ (simple $1/(p{+}1)$; step recursively lifted/certified; 56 partial classes labeled), closure-authority bad-root residues reconciled 5,375/5,375 (`data/l17_badroots_closures.jsonl`): matched statistics in `data/l17_matched.jsonl` (regenerated by `l17_matched.py` on the v2 resolved-mass artifact) — marginal $4481.3$ vs $4478$ (0.9993); hit-conditioned union $3477.9$ vs $3475$; per-class residuals mean $-0.00039$/SD $0.0121$; permutation null (2000) $z=-0.09$, empirical tails $0.61/0.53$ — no dependence. The class-weighted Haar-vs-window gap $\approx+0.0195$ is characterized as window-vs-$\mathbb Z_p$ mass, not a model failure. (iii) Cofactor layer proved = reciprocity-level: 197/197 singleton auto-$+1$ (given frozen/small $+1$); parity product $+1$ on 293/293; factorint individual signs closure-conditioned only (`data/l17_cofactor.jsonl`). (iv) Two-layer rate, predeclared cohort, unknowns both ways: family $c\in[0.017,0.80]$; wave-1+2 predicted $[114.3,293.0]$ vs 117; band brackets $p_0=0.0400$ (`data/l17_ratemodel_censored.jsonl`) |
| **L16: emergent-symbol structure — small layer is sieve-shaped** | **measured, per-row checks PROVED; layering EVIDENCE** | Lead probe, 7 classes × prime $k\le600$, sign-decorated emergent lists (`data/l16_char.jsonl`): the product character is $+1$ on all aligned members (parity), so content is per-prime; for every small emergent $p$ recurring $\ge4$ times (all 7 classes, zero counterexamples) the symbol $(x_0,d_0)_p=(2\alpha b/p)$ is a function of $k\bmod p$ — the small-emergent obstruction is a $p$-adic sieve: bad-signature root classes of $k\bmod p$ lifted by valuation parity (mass $1/(p{+}1)$ per simple root; singular roots require recursive lift certification or exclusion-labeled products). Two-layer split: clean-small rate 35–58%/class vs $\approx4\%$ closure → the big-cofactor layer (parity-forced sparse: one emergent prime, or a small product) carries its own $\sim20\times$ obstruction share |
| **L15: remainder law + counting law for H on the grid** | **PROVED per row; statistics exact (deterministic playback)** | Independent exact audit of all 293 closure members (`data/l15_remainders.jsonl`; 0 refusals, 0 alarms): stripped remainder $R$ **squarefree in every record**; **$R=1$ never occurs**; exactly two shapes — prime $R=p$ (197/293, never $pC^2$; median 48 digits, max 119) and factorint (96/293: squarefree products of $e\in\{2,3,4\}$ emergent primes, 76/15/5; 15 odd-$|E|$ rows exactly on the parity law). Counting law (`data/l15_density.json`): first member at $k=0$ in 207/293 (70.65%; L11 65.3% / ESC 80.6%); median 0, tail median 28, max 1694; per-prime-member emergent-free rate $\approx0.040$ (≈15.7× slower than pure Bateman–Horn first-prime wait); resistant subclass $a=1,w\equiv3\bmod4$ (30.8% at $k{=}0$; both max outliers). The uniform all-$w$ H now has its quantitative benchmark: "$p_0$ bounded below and the tail geometric on every grid class" |
| **L14: hypothesis H instance-verified on the whole grid — 293/293 classes** | **PROVED per row, suite-asserted (default 60-seed / extended all)** | Every aligned class (190 L11 + 103 ESC) carries a certified emergent-free member, frozen in `h10q.py::_L13_H_CLASSES` and replayed every run: ladder verdict *zero* on every row (197 prime-rung + 96 factorint), `ramified(x_0,d_0)` empty on 271 (extended: 242; 22/51 auxiliary refusals logged, never evidence). Three waves: frozen classes $k\le200$ (97) → BLS-prover rerun (+20) → alternate classes (+164): Hensel-root $a$-lifts, $f/\varepsilon$ variants, $w\equiv3\bmod4$ via the $a=1$ path; last cells: $(89,(-1,1))$ at $q=1097\ k{=}0$, $(89,(7,3))$ $q=2347\ k{=}0$ factorint, $(67,(-1,1))$ $f{=}67\ \varepsilon{=}-1\ q{=}811\ k{=}0$. Lead-replayed 293/293 (`data/l13h_replay.jsonl`); zero alarms. **What remains of H: the uniform statement over all $w$** (analytic, blocked at degree 2 unconditionally per LitScout-pinning) |
| **L13f: cofactor ladder; cell $(89,(-2,1))$ CLOSED** | **search-side discovery, lead-replayed, suite-frozen witness; method PROVED per row** | Zero-bad ⟸ (remainder of $P(b)$ after stripping $S\cup\operatorname{supp}(b)$) is a perfect square or a single proved prime with odd valuation (parity law then forces its symbol $+1$) — no full factorization needed. **Not a biconditional:** the ladder also returns *zero* after exact factorization whenever every proved odd-exponent factor has symbol $+1$, which is exactly how the 103 factorization-rung rows close. Ground truth: ladder-zero on 16/16 frozen zero-bad rows; INCONSISTENT = 0 across 5,767 decided rows. Cell box: 309,392 structured candidates → 10,515 aligned → 4,481 cofactor-big → **236 certified soluble** (133 proved-prime rung + 103 factorization rung) + 89 in the dense $|b|<2\cdot10^4$ box; 516 unproved-Jacobi rows evidence-tier only; 1,896 + 1,637 refusals recorded with reasons (never evidence). Frozen witness $(89,(-2,1))$: $(a,b)=(3,89/367)$, $\tau=19/37$ — emergent places $\{43,90947,204917,471137\}$ (all $+1$) and 29-digit proved-prime remainder $R$ ($v_R=1$, symbol $+1$); `_l7_tied_status` True; frozen in `_L13_RESIDUAL` (7 rows). The square-subfamily phase (290,928 decided, 0 hits, prior $\approx10^{-11}$/candidate) is recorded as non-evidence, per scope discipline |

QS2/L13-audit crash note (2026-08-18): the hi-q leg of the same driver ($a{=}17$, sq branch, $q\in3001..8000$) died on a benign `NameError` in the loop body — fix pending, no result claimed or counted.

## Findings (rev 21)
1. **L1 is published.** Cornelissen–Zahidi math/0006140, Ex. 2.2(a) + Rem. 2.4 ⟹ no one-witness
   definition of $\mathbb{Z}$ in $\mathbb{Q}$.
2. **L2 repaired into L2′ — a structure theorem** (advisory upheld; written proof in
   THEOREMS.md, computationally corroborated on the stated finite scopes). Cofinite-$S$
   integrality predicates are decidable boolean combinations of congruences; $D$ is
   model-relative, the dichotomy is not.
3. **A1 audit of Daans' 10** (unchanged): $12 - 2 = 10$ is the fixed point of that proof shape;
   three slack channels named.
4. **Sun's claimed 7 realizes A1 channel 2** — and his paper credits "AI's analysis of Daans's
   paper" for the approach. Audit A2, two passes:
   - **Corroborated mechanically**: identities; soundness (0 violations); Lemma 2.2 exhaustively;
     count bookkeeping $2+(3+3-1)=7$; witnesses at every odd prime $w < 300$.
   - **Refuted reading**: (9.5)∨(9.6) is NOT pointwise complete — deterministic counterexample
     $(2,-7/3,-9)$; independently, 8 of 61 witness targets *required* $\tau \notin \{0,\tau_1\}$.
     The Λ-family freezing is load-bearing, twice over.
   - **Editorial defects, then chain traced**: Prop 9.1 cites a phantom "Section 11" and
     self-referential "Sections 7 and 8" — stale draft numbering; the actual chain (§3 freezing →
     §4 Λ → §5 target-place Hensel → §7 global parameters → §8 Hasse–Minkowski + Lemma 8.1) is
     present and was traced by hand. **The local–global completeness proof is hand-traced, not
     machine-verified** — no suite can currently verify it.
   - **Verdict: 7 = established modulo refereeing** (unrefereed; editorially defective; no
     substantive gap found). Refereed anchor: 10 (Daans, JLMS 2024). Q5.7: $2 \le m \le 7$.
5. **The quantifier race = a fibre-dimension problem** (THEOREMS A3, DDF 2102.06941):
   minimal unknowns $= \mathrm{efd}_\mathbb{Q}(\mathbb{Q}\setminus\mathbb{Z}) + 1$; the race and
   the lower-bound wall are one geometric question, $\mathrm{efd} \in [1,6]$ ($[1,9]$ refereed).
6. **L4 (new, proved).** Sun's Lemma 5.1 holds for **every** odd prime power $5 \le q \le 25$
   and every admissible $\tau$ (exhaustion = proof at each $q$); $q = 3$ is vacuous
   ($\mathbb{F}_3$ has no admissible $\tau$). Hence his hypothesis $|k| > 25$ can be replaced by
   $|k| > 3$, and over $\mathbb{Q}$ the $|\kappa_w| \le 25$ patch component of
   $E_{\mathrm{exc}}$ shrinks from 8 places to $\{3\}$. Effectivity refinement; count unchanged.
7. **Channel-2 barrier, certified.** The $s = 0$ ternary restriction of $\Psi_\tau$ inherits
   soundness but provably loses completeness (24/224 integral instances on the witness grid,
   first at $(1,-3,2,0)$): a 2-witness certificate cannot be a specialization of $\Psi$ — the
   quaternion trace mechanism genuinely uses its 4th coordinate.
8. **Witness table.** Explicit frozen $(a, b, \tau)$ with $\Delta(Q_{a,b}) = \{2, w\}$ for all
   61 odd primes $w < 300$, found and re-verified exclusively under the proven engine
   (candidates whose Hasse decisions were uncertifiable were *skipped*, never trusted).
9. **L6 witness-tie breakthrough (conditional, not a theorem of completeness).**
   Set $a=1+2s$, identifying Sun's third block witness with Daans' $\Phi$ congruence
   parameter. The count is then $2+(3+2-1)=6$. Direct source verification of DDF v5
   Thm 1.4 and an independent adversarial audit found no counting defect:
   polynomial pullback preserves the $\exists_3$ upper bound; DDF fusion works for
   correlated subsets of the same base; projection adds $(b,s)$; denominators and
   $\delta_\tau\ne0$ are free because $\Phi$ forces
   $A=1+4a^2\equiv5\pmod8$, a nonsquare in $\mathbb{Q}_2$.
   - **Proved (W0–W1):** at the target $w$, the tied conic is soluble iff
     $\chi_w(-\delta_\tau A)=1$. A direct character-sum argument selects a
     nonsquare $A=1+4a^2$ and one of the two formula branches
     $\tau=0$ or $2a/A$ over **every** odd residue field, including $\mathbb F_3$.
   - **Proved (W2):** with $E_\tau=\mathbb Q(\sqrt{-\delta_\tau AB})$ and
     $M_\tau=16-\delta_\tau c^2-16ABs^2$, global solvability is exactly
     $$M_\tau=0\quad\text{or}\quad-\delta_\tau A M_\tau\in
       N_{E_\tau/\mathbb Q}(E_\tau^\times)$$
     ($M_\tau=0$ is solved by $(y,r)=(0,0)$, while $0\notin N(E_\tau^\times)$,
     so the disjunct cannot be absorbed into the norm group).
   - **Verified evidence:** fresh guarded searches give canonical certificates for
     inputs $z=w$ (the formula evaluates $h(a,b,z^3)$) at all 61 odd primes
     $w<300$ (32 $\tau=0$, 29 $\tau=2a/A$), with $A$ a nonsquare $w$-unit and
     $\Delta=\{2,w\}$. SSOT: `_L6_WITNESSES` in `h10q.py` is the authority;
     `data/l6_witnesses.jsonl` is its canonical serialization, byte-checked on
     every suite run — no export field exists without being derived from the
     authority. Search provenance (timings, found flags, raw buckets, the old
     rescue log, and the earlier unguarded run that exposed nonunit $A=5$ at
     $w=5$) lives in `data/superseded/` only.
   - **Still open:** force the W2 criterion — $M_\tau=0$ or the norm
     membership — while preserving $\Phi$, $v_w(b)=1$, and the W0 residue
     class. This is the exact global assembly lemma.

10. **L7 structure of the assembly problem (corollaries + bounded probe; status unchanged).**
    THEOREMS.md L7, machine checks `_verify_L7`:
    - **L7a (corollary of W2):** at odd places with $v(\delta_\tau A)=v(B)=0$
      the tied conic fails iff $v(M_\tau)$ odd and $\chi_v(-\delta_\tau
      AB)=-1$; where **additionally** $v(A)=0$ (all ternary coefficients
      units) the *untied* quaternary block is soluble regardless of the RHS
      (Chevalley–Warning + Hensel; isotropic regular ⇒ universal). The unit
      hypothesis is **necessary** and coefficient places carry real
      obstructions — all machine-witnessed in-suite: ternary anisotropic
      over $\mathbb{Q}_5$ for $a{=}1,b{=}11$ ($5\mid A$); a guarded
      $w{=}5$ shape insoluble at $v{=}3\mid B$ with $v_3(M_\tau)=0$; and a
      cancellation case $b{=}1/5$ with $v_5(AB)=0$ yet $v{=}5$ obstructed
      (so the bad set is *not* "$v\mid AB$"). **Obstruction support:
      $\{2,\infty\}\cup\mathrm{wild}(M_\tau)\cup\{v:v(A)\ne0\text{ or
      }v(B)\ne0\}$.** Tied dictionary: 232 places in-suite (337 extended); untied
      universality: 193 (283 extended) fully-unit places.
    - **L7b (corollary, conditional on L6 soundness):** no fixed
      $(s,b,\tau)$ has a rational-function witness $Y,R\in\mathbb{Q}(z)$ — it
      would specialize to a solution at some non-target $z=\pm2^k$. Probe: 8
      canonical witnesses ($w\le23$) × $|k|\le3$ × both branches = 0 solutions
      found on 224 cells (205 certified failures, 19 budget refusals unknown).
      This kills only uniform sections; separately, certified target-cell
      failures (e.g. $w{=}5$, $z{=}-5$, both branches — 7 such cells
      re-certified in-suite every run, covering **all five** probed pairs
      $w=3,5,7,11,13$) show each probed canonical pair fails to cover
      $\mathfrak m_w$ pointwise — whether some other fixed pair could
      remains open.
    - **L7c (corollary of W2 + Hilbert reciprocity):** a failing tied conic is
      obstructed at an **even** number of places, never one; repairs flip in
      pairs. Verified on 8 disproved cells in-suite per run; 42 more in the
      recorded session sweep.
    - **Bounded completeness probe (evidence only):** fixed canonical witnesses
      cover a minority of $z\in\mathfrak m_w$ (6/14, 6/16, 2/19, 2/22, 6/22 for
      $w=3..13$); free-$\Phi$ witness switching rescued **71/71** uncovered
      cells — 93/93 target cells positive; non-target probe: 0 solutions
      found on 224 cells (205 certified failures, 19 refusals). First-hit
      witnesses scatter (26/71 use $s=0$); ten spot-checks frozen as
      `_L7_RESCUES` and re-verified every run. **The assembly lemma itself
      remains open; no count or record claim changes.**
11. **L8 — the alignment wall (unconditional structure; status unchanged).**
    Exact value identities: $(Ab^2D_z)^2u_\tau = X^2 - D_\tau Y^2$ with
    $X=a^2z^6N_g$, $Y=4Ab^2D_z$, $D_0=P$, $D_1=AP$, $P=1-ABs^2$; both
    branches satisfy $D_\tau\cdot\operatorname{disc}(E_\tau)\equiv-2AbP$.
    - **L8a (absorption, proved):** odd $v$ with $v(\delta_\tau A)$, $v(B)$,
      $v(M_\tau)$ all even is soluble — the obstruction support is contained
      in $\{2,\infty\}\cup\{v:\ v(A),v(B),\text{ or }v(M_\tau)\text{ odd}\}$.
    - **L8b (parity wall, proved):** admissibility forces $v_2(P)=0$, hence
      $v_2(-2AbP)=1$: the source and target square classes $D_\tau$,
      $\operatorname{disc}E_\tau$ are distinct for **every** admissible
      witness (the source algebra can be split when $P\in\mathbb Q^{\times2}$),
      and $E_\tau$ is always a field (soundness protection re-proved).
    - **L8c (exact-matching emptiness, proved):** the two square-class
      identities that would kill wild value-primes identically (source $=$
      target; split target) have no admissible points; the first is also
      Pell-parametrized for $s\ne0$ and machine-swept (the $s=0$ slice
      forces $v_2(b)$ odd). General
      fixed-class matching is *not* ruled out (value-primes of a thin family
      need not equidistribute) — scope frozen in THEOREMS L8c.
    - **L8d (evidence):** all 10 frozen rescues succeed by wild-prime
      alignment (largest odd-mult primes $1.6\times10^5$–$4.5\times10^{21}$,
      symbols all $+1$); recorded sweep: 13/161 admissible pairs soluble
      ($\approx8\%$; grid in NOTES, not persisted in-suite), observed
      obstruction sets all even. Consistent with — not proving — the
      absence of any uniform mechanism beyond L8c's two.
    - **Computed residue data (recorded):** over $\mathbb Q(b)$ the class
      $(Au_1,-2b)$ ramifies at $b=0$ with residue $A$, and $-2b$ is
      **non-square** in the octic residue field of the degree-8 value
      divisor at all sampled $(s,Z)$ — finite-field splitting witnesses in
      NOTES. Any CTS/Schinzel-style attack must engage this residue.
12. **L9 — reciprocity steering (proved lemmas + the strongest assembly
    evidence so far; status unchanged).**
    - **L9 (proved):** with $T=\{2,\infty\}\cup\operatorname{supp}(x)\cup
      \operatorname{supp}(d)$, the product of $(x,d)_v$ over $T$ is $1$;
      one unchecked place is always forced. A single wild odd-multiplicity
      prime of the norm value is therefore *never* an obstruction once the
      controlled places align.
    - **L9a (proved, corrected hypotheses):** on the prime-$b$ family
      $b=\varepsilon q_1$ with $v_{q_1}=0$ imposed **factor by factor** on
      $\delta_\tau,A,a,z,D_z$, the new place contributes
      $(x,d)_{q_1}=(A|q_1)$ exactly ($v_{q_1}(x)=-4$, unit part $\equiv A$
      mod squares), controllable by Dirichlet classes; admissibility forces
      $A\equiv5\bmod8$, so $A$ is never a square and $+1$-classes exist.
      Audit P1: the earlier product-form guard $q_1\nmid2\delta_\tau A$ is
      **insufficient** — on $\tau=2a/A$, $\delta_\tau A=1$, and at
      $w=3,z=6,q_1=5=A$ one gets $v_5(x)=-5$ (odd) with symbol $-1$;
      frozen as an in-suite boundary regression. This still repairs the
      naive reduction (an arbitrary $b$-progression lets uncontrolled
      primes into $d$ — earlier advisory upheld).
    - **L9-T (evidence):** steered search covered **345/353 cells**
      (odd prime $w<100$; fixed 15-value $u$-pool filtered by
      $v_w(u)\ge0$: 11/14/13 retained at $w=3/5/7$, 15 elsewhere) with PROVED-soluble admissible
      witnesses — 344 steered (one wild prime, symbol concluded from
      $T\setminus\{q_0\}$ by L9, cross-replayed; largest forced prime 45
      digits), 1 smooth; 51 single-path (independent replay refused), 294
      double-proof. The 8 open cells are search-exhaustions with spread
      bad-place profiles (no local wall). Fixed-pair baseline was 13/161
      pairs: per-cell witness steering is the difference.
    - **Reduction (now split by L10):** per cell, assembly reduces to (i) an
      aligned class on the prime-$b$ family, and (ii) a Schinzel-type single
      prime value of the wild cofactor. L10 settles (i): the modulus freezes
      only the fixed places $S_0=\{2\}\cup\operatorname{supp}(\alpha)$, while
      the symbols that must all be $+1$ are
      $\Sigma=\{\infty\}\cup S_0\cup\{q_1\}$ — the moving prime $q_1$ is held
      by L9a and its residue class, never inserted into $N$ (that would
      contradict $\gcd(q_1,N)=1$). This supersedes the earlier
      $\{p\le H\}\cup$ data-support recipe: big fixed primes such as
      $40077920921911\mid D_z$ at $(73,5)$ are wild
      candidates, not controlled places; on $\tau=2a/A$ one has $\alpha=-1$
      and aligned classes are exhibited for all 164 $w\equiv1\bmod4$ cells;
      on $\tau=0$ the class is walled for prime $A$. (ii) remains
      conjectural. Frozen in THEOREMS L9/L10.
15. **L12 — the coprime-$b$ direction is closed, and the door is named.**
    The obvious repair of L11g — a cofactor $m$ with $(m\mid A)=-1$ — fails for a
    reason that generalizes: at every odd place of $b$ coprime to the cell data the
    unit part of $x$ is $\equiv A$, so $(x,d)_p=(p\mid A)$, and the product over
    $\{A\}\cup\operatorname{oddsupp}(b)$ telescopes to $-1$ **for every admissible
    $b$** — prime, semiprime, squarefull, rational alike. The $-1$ relocates, never
    cancels. What survives is exactly one hypothesis: coprimality. Breaking it with
    $b=\varepsilon fq$, $f\mid z$, lights up **103 of the 163 walled cells**
    (base-point alignment; class constancy not sampled, moduli $\sim10^{13}$).
14. **L11 — branches are free, the wall is not; and the reachable cells are now characterized.**
    Adjoining a branch costs nothing (soundness is $\tau$-uniform; $k$ conics in the
    same $(y,r)$ multiply to one polynomial), and the four square classes of
    $\langle-1,A\rangle$ are hit by explicit rational $\tau$, with $\tau^\dagger$
    making $\alpha$ a perfect square. **But the premise that the branch set was the
    deficient object is refuted by L11g:** the pairing $(x,d)_A(x,d)_q=-1$ holds on
    *every* branch whenever $v_A(z)\le0$, because admissibility forces $A\equiv5\bmod8$
    and hence $(2\mid A)=-1$. L10b is its $\tau=0$ instance. The escape is to make $A$
    share a prime with $z$, which is possible exactly when $z$ has a numerator prime
    $\equiv1\bmod4$ — giving **L11h: reachable iff such a prime exists, 190 of 353
    cells.** For the other 163, impossibility is *proved* for prime $A$ and is
    *evidence* for composite or rational $A$ (878,400 certificates, 0 hits) —
    the same scope split L10b already carries.
    An intermediate claim of 353/353 (L11e) was **false and is retracted**; the same
    defect invalidated 130 of L10c's 164 rows. Both were found by the adversarial
    audit of this layer, and both are frozen as in-suite regressions.

13. **L10 — the class side splits by branch (two proofs and a wall).**
    - **L10-0 (proved):** at $p\notin\operatorname{supp}(d)=S_0\cup\{q_1\}$,
      $v_p(d)=0$ so $(x,d)_p=(d|p)^{v_p(x)}$, which is $+1$ for even
      $v_p(x)$. Only $\Sigma=\{\infty\}\cup S_0\cup\{q_1\}$ matters;
      explicit exponents come from the exact Taylor shift of $P$ at $b_0$.
    - **L10a (proved):** $\tau=2a/A\Rightarrow\delta_\tau=1/A\Rightarrow
      \alpha=-1$ exactly, so no place divides $\alpha$ and only three symbols
      must be aligned.
    - **L10b (proved for prime $A$; evidence otherwise):** on $\tau=0$,
      $(x,d)_A\cdot(A|q_1)=(-2\varepsilon|A)=-1$ for both signs because
      $A\equiv5\bmod8$ (which needs $a$ odd) — two frozen symbols permanently anti-correlated, so
      **no aligned class exists** there. Not claimed for rational or
      non-squarefree $A$ (200 in-suite instances per run, unproved), and not claimed
      uniformly over the infinite $s$-range.
    - **L10c (proved per cell):** 164 exhibited classes = exactly the canonical $w\equiv1\bmod4$ cells,
      with $\gcd(q_1,N)=1$ so Dirichlet is non-vacuous; $N\in\{40,640\}$.
      Step (i) of the reduction is therefore **done** on those cells.
    - **Consequence — SUPERSEDED by L11g/L11h:** the claim that the reduction
      is established on 164 of 353 cells is **wrong** (130 rows invalid); the
      correct statement is 190 cells, characterized by a numerator prime
      $\equiv1\bmod4$, not by $w\bmod4$. Historical text: On the other
      189 it is proved unavailable **for prime $A$ with $a$ odd** and unobserved
      otherwise (finite scan: 11 $s$-values $\times$ 2 signs $\times$ 300
      primes per cell); universal impossibility is **not** claimed. Their
      L9-E pointwise certificates are
      unaffected — any wall is in the reduction, not in solubility.


16. **L23 — the unconditional sieve is real, but stops at the
    fixed-family detector problem.**  The semilinear dimension-\(1/2\)
    sieve supplies \(\gg X/(\log X)^{3/2}\) prime members clean below
    \(X^{0.49}\).  Degree \(8\) still permits sixteen large factors and
    reciprocity only forces an even bad count.  L31 corrects the
    analytic target: the two-large sector is the first term beyond BV,
    but its unsigned root-class expansion has a positive main term.
    One needs a Buchstab/Hilbert-detector asymptotic with positive
    zero-bad constant, or directly AP1; no unconditional member theorem
    follows.
17. **L24 — the apparent unscaled diagonal five-count is vacuous.**  Both
    exact self-coupled covers are finite flat of degree $8$, but
    orientation I is killed by its $(4,0,0)$ dyadic Newton polygon and
    orientation II by the complete mod-$16$ contradiction.  Their union is
    empty on $\Phi$.  The $4{,}026{,}282$ zero-hit search is
    corroboration only; the uniform theorem is local at $2$.
18. **L25 — fixed scaling changes the local obstruction into parity
    control.**  The type-I couplings $(2X,\rho)$ and $(2X,\rho/2)$ are
    $2$-adically admissible exactly on even- and odd-$s$ strata,
    respectively.  Their disjunction therefore has no uniform dyadic
    obstruction; global rational coupled points and any count improvement
    remain OPEN.
19. **L26 — an even pullback and reciprocal tie expose a quartic trace
    character.**  With $Z=z^2$ and
    $\eta=Z(b+1)/(Db)$,
    $P_{\rm rec}(b)=b^4T(b+b^{-1})$ for an exact quartic $T$.  The
    canonical square branch is automatically split at $2$, guarded W1
    targets lift, and $(2b|p)=(2(u+2)|p)$.  L28 proves the trace
    field and bad-sign extension uniformly; L31 proves one exact member,
    while the uniform member/parity step remains.
20. **L27 — one genuinely non-diagonal shear closes the displayed local
    coupling gaps.**  The fixed matrix
    $\left(\begin{smallmatrix}2&28\\s&1\end{smallmatrix}\right)$ covers
    every $2$-adic $s$-parity.  A Weil character-sum bound, patched by
    $43$ exact small-prime rows, gives a smooth point on the standard
    aligned target stratum for every odd $w$; one guarded real sample
    works.  The exact linear-$B$ factorization excludes the natural
    cancellation section, but the full bridge-specialized degree-$8$
    cover has no proved global rational point, so no five-count follows.
21. **L28 — the quartic trace field closes uniformly at \(2\).**  After
    shifting \(w=u-2\), the quartic Newton polygon excludes linear
    factors.  A Ferrari resolvent and five exact mod-\(8\) valuation
    cases exclude quadratic factors.  Finally
    \(N(2(\theta+2))\in A\mathbb Q_2^{\times2}\) is nonsquare, so the
    bad-sign cover has degree \(8\).  This is not a global member
    theorem and does not improve the count.

22. **L29 — the distinct reciprocal lift is now classified exactly.**
    For even \(t=v_2(Z)\), \(\theta^2-4\) is square in the trace field
    exactly for \(t=-2\) or \(t\ge2\).  The \(t=-2\) proof uses an
    ordinary resultant with two Newton slopes; the other strata use
    strong Hensel or exhaustive necessary factor congruences.  The map
    \(Z=8z^2/(1+z^2)\) preserves every positive odd valuation and forces
    the split local stratum.  Freezing \(b=w\rho^2\) exposes a fixed
    quadratic norm problem but retains obstruction pairs.  L31 proves
    one \(w=13\) member, not a uniform theorem or count improvement.
23. **L30 — fixing \(y=2\) halves the coupled-cover equation.**
    The exact equation is quadratic in \(m=\lambda^2\), hence even
    quartic in \(\lambda\).  A normalized Hensel argument solves every
    dyadic stratum; fibre-regular rows exist for every aligned
    \(w\ge5\), and \((a,b,z)=(5,3,3)\) gives an exact \(w=3\)
    strong-Hensel fibre.  A nearby empty control shows the selected
    \(a\) is load-bearing.  The same layer proves generic rigidity for
    every rational square-branch tie and separates the trace-base bundle
    from its mandatory reciprocal lift.  A global rational quartic root
    remains open.
24. **L31 — one reciprocal member is proved and the remaining inputs
    are exact.**  The fixed-field row
    \((w,a,Z,q_0,b,\rho,\lambda)=(13,3,169,1,13,14/13,12)\) has a
    recursive Pocklington certificate and every Hilbert symbol \(+1\).
    The constant-\(c=-64/25\) section gives a rational \(H_2\)-point,
    but its \(Z\) is not a cube; the cube step becomes two squarefree
    genus-\(2\) curves.  An external complete Magma calculation finds
    only \(z=0\) on both, closing this unique section.  The unsigned
    two-large-divisor expansion has a positive main term beyond BV, and
    AP1 is equivalent to intermediate H.  Neither the uniform member
    theorem, another controlled quartic point, fixed-family detector
    asymptotic, nor AP1 is proved.
25. **L32 — the dyadic and target base now has an explicit rational
    parameterization.**  Automatic-\(\Phi\) maps select smooth local
    L30 bases at \(2\) and every odd target.  Fixed proportional-trace
    menus and \(b=1\pm2a\) are proved non-uniform.  Emergent map primes
    remain uncontrolled.
26. **L33 — the unweighted pair main is positive.**  For
    \(0<\alpha<\beta<1/2\),
    \(\mathcal U_2\sim\frac18\log^2(\beta/\alpha)X\).  This puts a
    proved order-\(X\) unsigned main in the first beyond-BV window.
27. **L34 — a sharp root-pair upper bound below \(1\) would prove
    AP1.**  The root oversieve avoids exact odd valuations, and the
    unrestricted \(\vartheta=0.49\) pair envelope is
    \(0.9749604961\ldots\).  The conditioned prime-weighted estimate is
    unproved.


## The frontier, stratified (rev 21)
- **F1** ($\exists^k$ definability of $\mathbb{Z}$): $k=1$ is closed
  [CZ 2000], while $k=2$ is the first open case (sources and bounds:
  `THEOREMS.md`, L1 and consequences table).
- **F2** (quantifier records): the established unrefereed range is
  $2\le m\le7$ and the refereed upper bound is $10$ (`THEOREMS.md`, A1–A2).
  **Classical Schinzel H conditionally narrows the range to
  $2\le m\le6$**, equivalently
  $\operatorname{efd}\in[1,5]$ (`THEOREMS.md`, L22d and L24e–L34;
  `CONDITIONAL.md`, §1).  L23–L34 neither remove Schinzel H nor improve
  this record.
- **F3** (cofinite subrings): L2′ closes the elliptic-integrality route
  decidably; infinite excluded prime sets are necessary
  (`THEOREMS.md`, L2′).
- **F4** (low-variable H10/$\mathbb Q$): the one-variable case is
  decidable and the two-variable case is open; the
  $\forall_9\exists_7$ result is unrefereed and the
  $\forall_9\exists_{10}$ anchor is refereed (`THEOREMS.md`,
  consequences table).
- **F5** (unconditional member frontier): small-prime cleanliness for
  the selected linear/octic sequence is a theorem, and L31 proves one
  globally good reciprocal member at \(w=13\), but no per-target member
  theorem.  L30 supplies the degree-\(4\) constant-two cover, while L32
  parameterizes bases locally at \(2\) and every odd target; other map
  primes remain uncontrolled.  L33 proves that the unweighted
  root-pair sector has a positive order-\(X\) main.  L34 shows that the
  conditioned prime-weighted root-oversieve bound
  \[
  \mathcal P_2^{\rm root}
  \le(0.9749604961\ldots+o(1))\#\mathcal A_z
  \]
  would imply AP1 without exact odd valuations.  That beyond-BV sharp
  estimate is not proved.  AP1 remains equivalent to intermediate H
  (`THEOREMS.md`, L23–L34).


## Honest assessment
No Hilbert problem fell.  The record remains **CONDITIONAL**, not
unconditional: classical Schinzel H is unproved, and H10/$\mathbb Q$
remains open (`THEOREMS.md`, L22d and L24e–L34;
`CONDITIONAL.md`, §1).  L31's \(w=13\) row is a genuine global member,
but one row does not prove the uniform per-cell hypothesis.  Its exact
quartic point fails the cube/target condition.  AP1, the fixed-family
detector asymptotic, and every count improvement remain open.
Every finite scan retains its stated EVIDENCE or per-row scope.

Produced in this arc, to our knowledge new:
(1) L2′ — a written structure theorem replacing a flawed step;
(2) the two-pass audit of Sun's unrefereed seven-quantifier claim;
(3) L4 — the exact small-field sharpening of Sun's finite-field lemma;
(4) the certified obstruction to the naive two-witness specialization;
(5) the guarded witness tables and proven-primality engine;
(6) L6 — the audited six-unknown architecture;
(7) L8–L18 — the alignment, wall, steering, class, grid and density
structure, with every bounded statement kept at its recorded scope;
(8) L19–L20 — classical Schinzel H implies the member clause and an
aligned $f=w$ class exists for every cell;
(9) L21 — the reducible off-branch locus and the generic/per-grid
precursor;
(10) L22 — two independent fixed-branch irreducibility proofs for every
nonzero rational $Z$, closing the last non-Schinzel premise;
(11) L23 — the unconditional half-sieve, its exact two-large-divisor
barrier, and the Capell/fibration/norm/squareclass/multivariable
no-go reductions at their stated scopes; and
(12) L24 — exact finite-flat unscaled diagonal covers and the uniform
theorem that both their images are empty on $\Phi$ over $\mathbb Q_2$;
and
(13) L25 — the exact scaled-coupling parity theorem showing that
$(2X,\rho)$ and $(2X,\rho/2)$ jointly have no uniform $2$-adic
obstruction, with global rationality and target-compatibility left OPEN;
(14) L26 — the even-pullback reciprocal tie, quartic trace model,
automatic dyadic splitting and trace-character descent;
(15) L27 — the fixed $28$-shear, one-piece dyadic Hensel theorem and
all-odd-prime aligned-target character-sum theorem, with global
rationality left OPEN;
(16) L28 — the uniform dyadic irreducibility theorem for the L26 trace
quartic and its bad-sign Capell cover;
(17) L29 — the complete dyadic classification of the reciprocal lift,
the exact target-preserving compressor, and the fixed-field slice with
its section obstruction and scoped parity evidence
(`THEOREMS.md`, L23–L29; `agent://HalfSieveAudit`;
`agent://FibrationAudit`; `agent://DiagonalAudit`;
`agent://ReciprocalTieAudit`; `agent://TriangularShearAudit`); and
(18) L30 — the constant-two even quartic, its uniform dyadic and
all-target selected-fibre theorem, its fixed-\(a\) \(w=3\) control, the
general square-branch/shear rigidity theorems, and the exact trace-base
reciprocal-lift separation
(`THEOREMS.md`, L30; `l30_quartic_frontier.py`); and
(19) L31 — the exact \(w=13\) fixed-field reciprocal member, the unique
constant-\(c\) bridge section and its off-cube rational quartic point,
the externally complete genus-\(2\) closure of that section, the
positive-main correction to the two-large-divisor dispersion target,
and the AP1 equivalence
(`THEOREMS.md`, L31; `l31_frontier_push.py`;
`data/l31_magma_genus2.json`).

### Next

1. **Primary reciprocal-tie target:** extend the exact \(w=13\)
   globally good member to a per-target theorem.  The member itself is
   now PROVED, including the mandatory lift, Pocklington chain, and every
   Hilbert symbol.  The general slice is a sign-decorated degree-\(16\)
   binary-form norm problem; one row does not settle its density.
2. **Primary coupled-cover target:** leave the now-closed
   constant-\(c=-64/25\) section.  It gives a rational
   \(H_2(\lambda^2)=0\) point before the cube, but its \(Z\) is not a
   cube; Magma's complete rational-point calculation proves that both
   resulting genus-\(2\) cube curves have only \(z=0\).  A different
   section producing a target-specific rational point at all controlled
   places remains OPEN.
3. **Primary analytic target:** prove a pointwise fixed-family
   Buchstab/analytic-Hilbert-detector asymptotic with positive zero-bad
   constant.  The first unavailable term has
   \(p_1,p_2\ge X^{0.49}\) and
   \(p_1p_2>D_{\rm BV}\).  Its unsigned root-class expansion has a
   positive main term; signed character cancellation alone cannot
   discard it.
4. **Exact Schinzel replacement:** prove AP1 for every cell: one
   selected aligned class with \(Q(t)\) prime and
   \(R_{\rm bad}\le1\).  Even Hilbert parity makes AP1 equivalent to
   intermediate H.  L23 permits \(R_{\rm bad}=0,2,\ldots,16\), so it
   does not prove AP1.
5. **Scope guard:** the single \(w=13\) construction is not a uniform
   member theorem.  Classical Schinzel H remains the sole conjectural
   input and the six-count remains CONDITIONAL.  No unconditional count
   or H10/\(\mathbb Q\) conclusion changes.
6. **Optional side geometry:** the L22 Noether polynomial nonvanishing
   lemma and absolute rationality of the L24 total surfaces remain OPEN
   but are not chain obligations.  Continue the source watch for Sun's
   unrefereed paper and the priority search for L2′; neither affects the
   unchanged L22 implication (`THEOREMS.md`, A2 and L2′).
