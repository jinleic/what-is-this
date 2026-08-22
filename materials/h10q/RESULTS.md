# H10 over Q — session results (2026-08-17, rev 10)

## Deliverables
- `NOTES.md` — problem map with exact theorem statements and arXiv ids.
- `THEOREMS.md` — precise frontier statements: L1 (published, attributed), L2′ (proved here,
  full written proof), L4 (finite-field exhaustion), L5 (proved identity:
  Sun's Lemma 5.1 count $=\#E_\lambda(\mathbb{F}_q)/4$), **L6 (new: audited
  6-unknown witness-tie architecture; canonical two-branch target theorem and exact norm
  obstruction proved; global assembly explicitly open)**, A1/A2 (audits of Daans' 10 and
  Sun's 7), A3 (efd reframing).
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
| **L19: Schinzel H $\Rightarrow$ per-class H** | **PROVED implication (conclusion CONDITIONAL on Schinzel)** | For a fixed verified aligned class with $F(t)=c\,G(t)$ ($c$ $S$-supported, $G$ primitive, positive lead), Schinzel H for $\{q_1+Nt,G(t)\}$ gives infinitely many ladder-zero members: the $S$-strip leaves one odd place $R=G(t)$, $v_Q(P(b))=0$, all other outside symbols $+1$, and reciprocity forces $(x_0,d_0)_R=+1$. Schinzel hypotheses verified exactly on 293/293 canonical pairs; 24-row prime-rung cross-check (`data/l18_schinzel_implies_h.jsonl`). **The bespoke analytic hypothesis is replaced by a named classical conjecture** |
| **L19: the canonical 5-wall and its break** | **PROVED (identity + criterion, canonical protocol); PROVED break** | $\mathrm{Hilb}_5=-(w|5)(q_1|5)$ vs wild $(5|q_1)$ $\Rightarrow$ canonical ESC 5-walled iff $(w|5)=+1$; with L11 iff $w\equiv1\bmod4$, canonical obstruction is exactly $w\equiv11,19\bmod20$. Audit: 53 primes in $[101,400]$, 26 L11 / 13 ESC / 14 obstructed, **0 mismatches** (`data/l18_fivewall.jsonl`). Break: $a{=}7,\tau{=}0,f{=}w$ closes $131,139,151$ (L12b non-coprimality door; no conflict with the L10 $\tau=0$ coprime wall) |
| **L19: the 5-wall is always breakable (class side)** | **PROVED** | Factor-by-factor: $\prod_{p\mid A}(x,d)_p(A\mid q_1)=(-2\varepsilon\mid A)(f\mid A)$, so with $A\equiv5\bmod8$ the L10b anti-correlation breaks **iff $(f\mid A)=-1$**. Uniform existence: $\sum_r(1+4r^2\mid w)=-1$ gives $(w+1)/2$ good residues, CRT fixes the unit hypotheses, Dirichlet supplies $q_1$ — so every canonically obstructed $w\equiv11,19\bmod20$ has an aligned $\tau=0$, $f=w$ class. Audit 16,744 candidates, 908 aligned, **0 mismatches** (`data/l19_tauzero.jsonl`). Scope: class side only; members still need Schinzel (OPEN). Fixed $a=7$ is breakable iff $(w\mid197)=-1$ |
| **L19: off-grid horizon** | **PROVED per row (lead-replayed)** | Closed $w\in\{101,\dots,179\}$ — all sixteen tried cells at $u=-1$ plus $[131,[5,1]]$; **no off-grid cell remains open**. $w=179$: $a=7$, $\varepsilon=-1$, $f=179$, $q_1=Q=251$, $k=0$, decided fraction 94.7% (10,282/10,858), 222-digit max cofactor (`data/l19_cell179.jsonl`). Divergence model: L18 tail closes 51% of the excess multiplier / 78% of the log gap vs the 15.67$\times$ Bateman–Horn ledger gap; stationary model misses the selected $k=0$ spike (`data/l18_divergence_model.jsonl`) |
| **L18: step-(ii) density law** | **PROVED (local bound; no-uniform-mass for the 293 current classes); EVIDENCE (measured exponents)** | $m_3=m_5=m_7=0$ and $m_p\le 8/(p+1)$ for $p\ge11$ (singular roots incl.); lead-verified on 102,439 local pairs, 0 violations. Chebotarev exponent $c=1/2$ in every class (293 replayed irreducibility certificates + 293 replayed nonsquare witnesses) $\Rightarrow \prod_{p\le X}(1-m_p)\asymp(\log X)^{-1/2}$: **no cutoff-uniform positive clean mass**. Measured $c$: 0.4778 / 0.4801 / 0.4590 (`data/l17_stepii.jsonl`) |
| **L18-horizon wave 2** | **PROVED per row (lead-replayed); one PROVED protocol obstruction** | Closures $w\in\{113,127,137\}$, $u=-1$ (ESC $f{=}w$ for 113/127; L11 $a{=}87$ for 137) extend the off-grid set to $\{101,103,107,109,113,127,137\}$. $w=131$: L11 empty ($(-1|131)=-1$) and ESC forces $\mathrm{Hilb}_5=-(5|q_1)$ — 336 attempts, 326 obstructed, 0 aligned, cell OPEN. Composite $w$ out of scope (place constructor); pseudo-cell $[21,[-1,1]]$ closes at honest prime $f=3$ (`data/l17_horizon{113,127,131,137}.jsonl`, `data/l17_composite_w.jsonl`) |
| **L17-horizon: off-grid closures $w\in\{101,103,107,109\}$** | **PROVED per row (independently replayed by lead)** | Verified horizon beyond $w\le97$: $[101,[-1,1]]$ L11 $a{=}5$ member $k{=}16$ factorint rung; $[103,[-1,1]]$ ESC $f{=}103$ $k{=}0$; $[107,[-1,1]]$ ESC $f{=}107$ $k{=}0$; $[109,[-1,1]]$ L11 $a{=}71$ $k{=}0$ (auxiliary tied/ramified budget-refused — accepted L14 closure form); off-grid machinery: fresh-L11 roots + `_l10_class_cert` / temporary `_L12_ESCAPE` + L13 ladder (`data/l17_horizon{101,103,107,109}.jsonl`) |
| **L17: mechanistic p-adic sieve — hit-conditioned matched fit exact; cofactor = reciprocity-level proved; rate band-only** | **per-row PROVED; matched model validation EVIDENCE (marginal 0.9993; union z=-0.09, empirical tails 0.61/0.53 under marginal-preserving permutation)** | (i) Residue-determination at full grid: 293/293 classes, 0 counterexamples (`data/l17_sieve.jsonl`, 8,760 members, $k\le119$). (ii) Exact bad masses $m_p$ (simple $1/(p{+}1)$; step recursively lifted/certified; 56 partial classes labeled), closure-authority bad-root residues reconciled 5,375/5,375 (`data/l17_badroots_closures.jsonl`): matched statistics in `data/l17_matched.jsonl` (regenerated by `l17_matched.py` on the v2 resolved-mass artifact) — marginal $4481.3$ vs $4478$ (0.9993); hit-conditioned union $3477.9$ vs $3475$; per-class residuals mean $-0.00039$/SD $0.0121$; permutation null (2000) $z=-0.09$, empirical tails $0.61/0.53$ — no dependence. The class-weighted Haar-vs-window gap $\approx+0.0195$ is characterized as window-vs-$\mathbb Z_p$ mass, not a model failure. (iii) Cofactor layer proved = reciprocity-level: 197/197 singleton auto-$+1$ (given frozen/small $+1$); parity product $+1$ on 293/293; factorint individual signs closure-conditioned only (`data/l17_cofactor.jsonl`). (iv) Two-layer rate, predeclared cohort, unknowns both ways: family $c\in[0.017,0.80]$; wave-1+2 predicted $[114.3,293.0]$ vs 117; band brackets $p_0=0.0400$ (`data/l17_ratemodel_censored.jsonl`) |
| **L16: emergent-symbol structure — small layer is sieve-shaped** | **measured, per-row checks PROVED; layering EVIDENCE** | Lead probe, 7 classes × prime $k\le600$, sign-decorated emergent lists (`data/l16_char.jsonl`): the product character is $+1$ on all aligned members (parity), so content is per-prime; for every small emergent $p$ recurring $\ge4$ times (all 7 classes, zero counterexamples) the symbol $(x_0,d_0)_p=(2\alpha b/p)$ is a function of $k\bmod p$ — the small-emergent obstruction is a $p$-adic sieve: bad-signature root classes of $k\bmod p$ lifted by valuation parity (mass $1/(p{+}1)$ per simple root; singular roots require recursive lift certification or exclusion-labeled products). Two-layer split: clean-small rate 35–58%/class vs $\approx4\%$ closure → the big-cofactor layer (parity-forced sparse: one emergent prime, or a small product) carries its own $\sim20\times$ obstruction share |
| **L15: remainder law + counting law for H on the grid** | **PROVED per row; statistics exact (deterministic playback)** | Independent exact audit of all 293 closure members (`data/l15_remainders.jsonl`; 0 refusals, 0 alarms): stripped remainder $R$ **squarefree in every record**; **$R=1$ never occurs**; exactly two shapes — prime $R=p$ (197/293, never $pC^2$; median 48 digits, max 119) and factorint (96/293: squarefree products of $e\in\{2,3,4\}$ emergent primes, 76/15/5; 15 odd-$|E|$ rows exactly on the parity law). Counting law (`data/l15_density.json`): first member at $k=0$ in 207/293 (70.65%; L11 65.3% / ESC 80.6%); median 0, tail median 28, max 1694; per-prime-member emergent-free rate $\approx0.040$ (≈15.7× slower than pure Bateman–Horn first-prime wait); resistant subclass $a=1,w\equiv3\bmod4$ (30.8% at $k{=}0$; both max outliers). The uniform all-$w$ H now has its quantitative benchmark: "$p_0$ bounded below and the tail geometric on every grid class" |
| **L14: hypothesis H instance-verified on the whole grid — 293/293 classes** | **PROVED per row, suite-asserted (default 60-seed / extended all)** | Every aligned class (190 L11 + 103 ESC) carries a certified emergent-free member, frozen in `h10q.py::_L13_H_CLASSES` and replayed every run: ladder verdict *zero* on every row (197 prime-rung + 96 factorint), `ramified(x_0,d_0)` empty on 271 (extended: 242; 22/51 auxiliary refusals logged, never evidence). Three waves: frozen classes $k\le200$ (97) → BLS-prover rerun (+20) → alternate classes (+164): Hensel-root $a$-lifts, $f/\varepsilon$ variants, $w\equiv3\bmod4$ via the $a=1$ path; last cells: $(89,(-1,1))$ at $q=1097\ k{=}0$, $(89,(7,3))$ $q=2347\ k{=}0$ factorint, $(67,(-1,1))$ $f{=}67\ \varepsilon{=}-1\ q{=}811\ k{=}0$. Lead-replayed 293/293 (`data/l13h_replay.jsonl`); zero alarms. **What remains of H: the uniform statement over all $w$** (analytic, blocked at degree 2 unconditionally per LitScout-pinning) |
| **L13f: cofactor ladder; cell $(89,(-2,1))$ CLOSED** | **search-side discovery, lead-replayed, suite-frozen witness; method PROVED per row** | Zero-bad ⟺ (remainder of $P(b)$ after stripping $S\cup\operatorname{supp}(b)$) is a perfect square OR a single proved prime with odd valuation (parity law forces its symbol $+1$) — no full factorization needed. Ground truth: ladder-zero on 16/16 frozen zero-bad rows; INCONSISTENT = 0 across 5,767 decided rows. Cell box: 309,392 structured candidates → 10,515 aligned → 4,481 cofactor-big → **236 certified soluble** (133 proved-prime rung + 103 factorization rung) + 89 in the dense $|b|<2\cdot10^4$ box; 516 unproved-Jacobi rows evidence-tier only; 1,896 + 1,637 refusals recorded with reasons (never evidence). Frozen witness $(89,(-2,1))$: $(a,b)=(3,89/367)$, $\tau=19/37$ — emergent places $\{43,90947,204917,471137\}$ (all $+1$) and 29-digit proved-prime remainder $R$ ($v_R=1$, symbol $+1$); `_l7_tied_status` True; frozen in `_L13_RESIDUAL` (7 rows). The square-subfamily phase (290,928 decided, 0 hits, prior $\approx10^{-11}$/candidate) is recorded as non-evidence, per scope discipline |

QS2/L13-audit crash note (2026-08-18): the hi-q leg of the same driver ($a{=}17$, sq branch, $q\in3001..8000$) died on a benign `NameError` in the loop body — fix pending, no result claimed or counted.

## Findings (rev 9)
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


## The frontier, stratified (rev 9)
- **F1** ($\exists^k$ definability of $\mathbb{Z}$): $k=1$ closed [CZ 2000]; $k=2$ first open
  case, guarded only by conjectures (Mazur for surfaces; strong BL per Koenigsmann Cor 23).
- **F2** (quantifier records): established/unrefereed $\forall$-side
  $2\le m\le7$ (Sun; $\le10$ refereed). **L6 conditionally narrows this to
  $2\le m\le6$**; equivalently $\mathrm{efd}\in[1,5]$, but only after its open
  assembly lemma. The naive witness-side specialization is certified dead; L6
  instead removes a parameter direction.
- **F3** (cofinite subrings): L2′ closes the elliptic-integrality route *decidably*; infinite
  excluded prime sets are necessary, qualitatively.
- **F4** (low-variable H10/$\mathbb{Q}$): $n=1$ decidable; $n=2$ open (effective Faltings/BSD).
  $\forall_9\exists_7$-theory of $\mathbb{Q}$ undecidable (Sun, unrefereed;
  $\forall_9\exists_{10}$ refereed via Daans Cor 6.2).


## Honest assessment
No Hilbert problem fell; every claim above carries its exact scope. Produced this session, to
our knowledge new:
(1) L2′ — written proof + bounded computational corroboration, replacing a flawed step;
(2) an independent two-pass audit of the 14-day-old claimed record 7: statement-level trace of
the completeness chain (hand-traced, not machine-verified) + mechanical corroboration at every
accessible point + the fixed-$\tau$ refutation + identification of the stale cross-references
[INFERENCE: no other public independent verification of 2607.28606 found];
(3) **L4** — a proved sharpening of his Lemma 5.1 constant ($25 \to 3$), shrinking
$E_{\mathrm{exc}}$ over $\mathbb{Q}$ from 8 patched places to one;
(4) a certified barrier: 2-witness specialization of $\Psi$ is impossible (soundness inherited,
completeness provably lost);
(5) the witness table to $w < 300$ with 8 targets exercising the Λ-license beyond
$\{0, \tau_1\}$;
(6) a proven-primality arithmetic engine with adversarial regression tests, after the advisory's
bug class was confirmed live at $w = 61$.
(7) **L6** — an audited six-unknown architecture, a proved canonical two-branch
target selector, the exact quadratic norm obstruction, and 61 guarded canonical
global examples. This is genuine progress on Daans Question 5.7, but **not** a
new quantifier record until the global assembly lemma is proved.
(8) **L8** — the alignment wall: an exact-value-identity framework
($u_\tau$ as binary-form values $X^2-D_\tau Y^2$; the invariant
$D_\tau\operatorname{disc}E_\tau\equiv-2AbP$ for both branches), the proved
2-adic parity wall $v_2(-2AbP)=1$ separating the source and target square
classes on every admissible witness, the proved emptiness of both exact
square-class matching gauges, an absorption lemma bounding the obstruction
support by the odd-valuation places, and scoped evidence that the ten
frozen rescue witnesses succeed by wild-prime alignment. This localizes
*why* assembly is hard without changing its status.
(9) **L9** — reciprocity steering: the forcing lemma and the prime-$b$
symbol law (both proved after audit repair, machine-checked with both
refuted drafts frozen as regressions), an engine that certifies per-cell
solubility from trial division plus a bounded number of proven-primality
checks (no hard factoring, no wild-prime luck),
and the largest assembly evidence to date:
345/353 cells covered with in-suite-replayed certificates, 8 named
residual cells. Assembly itself remains OPEN; count 6 stays conditional.
(10) **L10–L13 — the class side, closed and verified.** L11h's exact
characterization (190 reachable: $z$ has a numerator prime $\equiv1\bmod4$),
the L12a generalized wall killing every coprime $b$ on every branch, L13a's
1. **L6 assembly, primary target — now exactly localized.** The class side is
   solved and verified (L11h + L12a + L13a): every grid cell either has an
   aligned class (190 prime-$b$ + 103 non-coprime escapes) or a directly
   certified witness (60 escape2 + 7 residuals — all replayable in-suite;
   total coverage **353/353, complete since 2026-08-18** when L13f's
   cofactor ladder closed $(89,(-2,1))$). What remains is a single
   analytic statement: **some member of each aligned class is emergent-free**
   (all odd-valuation primes of $P(b)$ outside $S\cup\{Q\}$ carry symbol
   $+1$; the bad count is automatically even, so the condition is "zero
   bad"). Certify it per-class where factoring permits (cofactor-ladder
   instances now number in the hundreds for the last cell alone),
   then engage the literature input (LitScout's verified report,
   `data/litscout_h10q.md`: square-class prescribed values unconditional
   only to degree 2, degree 3 conditional on elliptic Parity, degree 8
   fully open) for the Bunyakovsky/Schinzel form.
2. Watch arXiv 2607.28606 for v2; remaining referee items in Sun's own 7-proof are §3 dyadic
   freezing and §7 simultaneous norm approximation. Offer L5's exact identity for a corrected
   exceptional-set statement.
3. Deeper priority search for L2′ (Shlapentokh's book; MathOverflow); if clean, write a short
   expository note.
4. Implement Koenigsmann Cor 2 ($\mathbb{Q}\setminus\mathbb{Z}$ Diophantine) atop `h10q.py`.
5. arXiv 1309.0441, “Undecidability in number theory,” is Koenigsmann's survey
   (often misattributed to Poonen, who has a 2008 Notices article with the same title).
