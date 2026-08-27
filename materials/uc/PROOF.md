# The union-closed inequality at $t = 0.3820660112501052$: candidate proof and evidence chain

This document separates ordinary mathematics, cited standard results,
machine-verified finite arithmetic, and numerical evidence.  The Arb
certificate is replayable; the support reduction, Margin Lemma, and
set-family bridge are human-audited rather than formalized.  The overall claim
must not be described as end-to-end machine checked.

**Candidate claim.** Set
$$\psi = \tfrac{3-\sqrt5}{2} = 0.3819660112501051\ldots,\qquad
t_{\mathrm{cert}} = 0.3820660112501052,\qquad \alpha = 0.0356069.$$
The human reduction and replayed certificate imply that every union-closed
family $\mathcal{F} \neq \{\varnothing\}$ has an element in at least a
$t_{\mathrm{cert}}$ fraction of its sets.  The exact comparison
$t_{\mathrm{cert}}\ge\psi+10^{-4}$ is machine checked.  Novelty and priority
are addressed separately in `LITERATURE_ORIGINALITY.md`; external mathematical
review remains outstanding.

The margin is *rational-exact*: for $r \in [0, 3/2]$, $r \ge \psi
\iff r^2 - 3r + 1 \le 0$ — a rational check by
`cert3_par.target_relation_holds`.

## 1. From the certificate to union-closed sets (bridge)

Cambie's Question 2 (**source, arXiv:2212.12500v2, Q2 + Section 4**) asks for
the following inequality when the common expectation is strictly below \(c\).
For iid $p,q\sim\mu$ on $[0,1]$ and $(p,r)$ any self-coupling,
$$(1-\alpha)\,\mathbb{E}h(p+q-pq) + \alpha\,\mathbb{E}h\!\big(\max(p,r,\min(p+r,\tfrac12))\big) \;\ge\; \mathbb{E}h(p).$$
Our certificate covers the stronger closed domain
\(\mathbb E p\le t_{\mathrm{cert}}\).  The revised paper and `AUDIT.md`
reconstruct Section 4's sequential Bernoulli coupling and prove the required
strict first-coordinate margin, so the implication is no longer an opaque
black-box citation.

Two exact human bridge steps, with finite symbolic controls in `bridge_uc.py`:

1. $s^\star(p,r) := \min(\max(\tfrac12,p,r),\ \min(p+r,1))
  = \max(p,r,\min(p+r,\tfrac12))$: three exhaustive cases on $s = p+r$ and
  $m = \max(p,r)$, closed by inspection ($s\le\tfrac12$;
  $m \le \tfrac12 \le s$; $m \ge \tfrac12$), guarded by exact rational controls.
2. Any coupling $\pi$ of $\mu$ with itself symmetrizes
  $(\pi+\pi^{\mathsf T})/2$ and folds to $\nu$ on
  $\Delta = \{0\le p\le q\le 1\}$ through
  $(p,r)\mapsto(\min(p,r),\max(p,r))$, preserving both marginals and the cost
  $\int h(s^\star)\,dM$: an identity of integrals, checked term-by-term by the
  `functional_identity` exact rational schema (all three $s^\star$ branches).

So the certified functional — the exact quantity minimized in
`thmB3_proof.py` —
$$F(\mu) = (1-\alpha)\,Q + \alpha\,C - L,\qquad
Q = \iint h(p+q-pq)\,d\mu d\mu,\quad L = \int h\,d\mu,\quad
C = \inf_{M \in \mathrm{SymCpl}(\mu)} \int h(s^\star)\,dM,$$
satisfies: $F(\mu)\ge0$ on every feasible $\mu$ implies the union-closed bound
at $t$ by the human entropy-chain proof in the revised manuscript.

**Exhaustive finite control (2026-08-27).**
`verification/entropy_bridge_exhaustive.py` rebuilds that entropy-chain proof's
construction from the corollary statement alone — importing neither
`bridge_uc.py` nor any campaign or certificate module — and runs it on every
nonempty family of subsets of \([n]\) for \(n\le4\): 65,808 families,
1,631,880 coupled prefix states, 48.2 s on one core.  Exact rational checks
cover the four coupling masses, both Bernoulli marginals, the \(s^\star\) OR
identity, prefix-by-prefix uniformity of \(A\) and \(C\),
\(\mathcal L(P_i)=\mathcal L(R_i)\) with mean the frequency of element \(i\),
and the chain rule; 256-bit Arb certifies both conditioning inequalities, each
enclosure being provably nonnegative or a structural tie inside
\(\pm2^{-200}\).  All 5,096 enumerated union-closed families satisfy the
\(t_{\mathrm{cert}}\) conclusion, the extreme case being exactly \(1/2\).  The
universal statement is still the human induction; this is a finite control,
and its fail-closed tests reject a broken \(s^\star\) clip, a constant prefix
probability, and an inflated union bit.

## 2. Exact scaffold: Theorems A and A′

Natural-log entropy $H = \ln2\cdot h$ inside proofs.  **Lemma 1 — PROVED,
`decomposition.py`**: $T(u) := u + (1-u)\ln(1-u) = \sum_{n\ge2}
\frac{u^n}{n(n-1)}$ (displayed series rearrangement), so $T(xy)$ is PSD:
$\iint T(xy)\,d\nu d\nu = \sum_{n\ge2} M_n^2/(n(n-1))$.

**Theorem A — PROVED, `decomposition.py`.**  With $x = 1-p$, $A = \int (-x\ln x)d\nu$, $B = \int x\,d\nu$, $G = A+B = \int x(1-\ln x)d\nu$:
$$\iint H(xy)\,d\nu d\nu = G(\nu)^2 - A(\nu)^2 - \iint T(xy)\,d\nu d\nu,$$
i.e. $H(xy) = g\otimes g - a\otimes a - T(xy)$: SymPy simplifies LHS − RHS to
$0$ identically for a symbolic 3-atom measure — an exact algebraic tautology.
The OR-entropy kernel has exactly one positive square, in the direction $G$.
(An earlier stronger claim is *retracted* in the same file: the coupling
minimum $C$ is convex — it is a supremum of linear functionals by
Monge–Kantorovich duality — with an explicit witness.)

**Theorem A′ — PROVED, `reduction.py`.**  Let
$$W(x,y) = xy - xT(y) - yT(x) + T(xy) = \sum_{n\ge2}\frac{(x-x^n)(y-y^n)}{n(n-1)},$$
a PSD kernel, and $E(\mu) = \frac1{\ln2}\iint W(1-x,1-y)\,d\mu d\mu \ge 0$.
Then $\ Q(\mu) = 2B(\mu)L(\mu) - E(\mu)$, $B = 1 - \mathbb{E}p$;
equivalently $\ln2\,E = \sum_{n\ge2}(M_1-M_n)^2/(n(n-1))$, a sum of squared
moment defects $M_1 - M_n \ge 0$ (rearrangement termwise, controls at
$10^{-30}$).  Recorded consequences:

* **Corollary 1**: the certified form —
  $F = [2(1-\alpha)B - 1]\,L + \alpha C - (1-\alpha)E$.
* **Diagonal**: $W(x,x)/\ln2 = 2x\,h(x) - h(x^2)$ (checked at the golden
  point to $10^{-30}$).
* **Sharpness (Corollary 2)**: at $\alpha = 0$,
  $F(\delta_p) = h((1-p)^2) - h(p)$, so $F(\delta_p) \ge 0 \iff
  (1-p)^2 \ge p \iff p \le \psi$: $\psi$ is exactly where the PSD defect
  $E$ first overtakes $(2B-1)L$, the golden ratio entering only via
  $x^2 = 1 - x$.  A uniqueness lemma ($dF/dp < 0$ strictly on
  $[p_0,\tfrac12)$, $p_0 = 1-1/\sqrt2$, every factor's sign argued casewise
  — proved, not sampled) makes that root unique.
* **Pointwise route insufficient**: dropping the variance gives only
  $0.3475266\ldots < \psi$ (maximiser at $x \approx 0.8835$); the variance
  term is load-bearing, hence a *support reduction* is needed.  (All four
  items: controls in `reduction.py`.)

## 3. Support reduction: Theorem B′′′

**Theorem B′′′ — HUMAN-AUDITED, `thmB3_proof.py`** (the text proof uses
explicitly identified Riesz/Banach--Alaoglu/Stone--Weierstrass and Bauer
inputs; its executable controls are sampled):
$$\min\{\,\Phi_{\mathrm{exact}}(\nu) : \nu\in P(\Delta),\;
\mathbb{E}_{\mu_\nu}p \le t\,\},\qquad
\Phi_{\mathrm{exact}}(\nu) = (1-\alpha)Q(\mu_\nu) + \alpha\int c\,d\nu - L(\mu_\nu),$$
with $\mu_\nu = \int\frac{\delta_p+\delta_q}{2}d\nu$, $c(p,q) =
h(s^\star(p,q))$, exists and equals $\inf\{\,F(\mu) : \mathbb{E}_\mu p \le
t\,\}$; it is attained at a $\nu$ with **at most 2 atoms** (marginal ≤ 4
atoms).  Proof skeleton:

1. *Compactness*: $\Delta$ compact metric; $P(\Delta)$ weak-* closed in the
   unit ball of $C(\Delta)^\ast$ (Riesz representation and Banach–Alaoglu,
   quoted literally from van Neerven, arXiv:2112.11166v7, Thms 4.2, 4.50).
2. *Continuity*: $h$, $s^\star$ (min/max compositions), $T$, $W$ continuous
   on closed domains (from $r\ln r\to0$); $E$ weak-* continuous via
   Stone–Weierstrass; so $\Phi_{\mathrm{exact}}$ attains its minimum on the
   weak-* compact feasible set $H_t$.
3. *Slice concavity*: fix a minimizer $\nu^\star$, $\beta^\star =
   B(\nu^\star)$; on $K_{\beta^\star} = \{B = \beta^\star\}$ the bilinear
   $BL$ is a constant by A′, so $\Phi_{\mathrm{exact}} =
   [2(1-\alpha)\beta^\star - 1]L + \alpha\int c\,d\nu - (1-\alpha)E$ is
   **concave** — a PSD quadratic form is convex (S3), and $\alpha \le 1$ is
   the only place $\alpha\le1$ is used.  The file's control detects failure
   when the fixed-$B$ premise is removed.
4. *Bauer*: a minimizer that is extreme in $K_{\beta^\star}$ exists (Bauer
   1958, DOI 10.1007/BF01898615, quoted via Bru & de Siqueira Pedra,
   arXiv:1610.03411, Lemma 3.3).
5. *Extreme points ≤ 2 atoms*: $K_{\beta^\star}$ is a Choquet slice (Dirac
   extremes) cut by one moment; Winkler, Math. Oper. Res. 13(4) (1988)
   581–587 Thm 2.1(b) — exact statement quoted via Pinelis, arXiv:1204.0249,
   Thm 12 — gives ≤ 2; an independent three-point perturbation argument in
   the same file proves the bound without the citation.

(The superseded B′, B″ — three orbits from freezing $G$ or $B+L$ — are kept
as documented history in `decomposition.py`, `margin_lemma.py`; B′′′ is
tight: the obstruction uses exactly two orbits.)

## 4. The Margin Lemma and the certified functional

Write $W(x,y)=\langle\psi(x),\psi(y)\rangle\ln2$; then
$E=\|\int\psi\,d\mu\|^2$ and
\[
\|\psi(x)\|^2=\frac{W(x,x)}{\ln2}
=\operatorname{rh}(p)=\rho(p)h(p),
\]
where
$$\operatorname{rh}(p):=2(1-p)h(p)-h\big((1-p)^2\big)
=2(1-p)h(p)-h(2p-p^2),\qquad \rho:=\operatorname{rh}/h.$$

**Margin Lemma — HUMAN-AUDITED, `margin_lemma.py`.**  Put
$\sigma(\mu) = \int\sqrt{\operatorname{rh}}\,d\mu$ and, for $L>0$,
$$\Lambda(\mu) := 2(1-\alpha)B - 1 + \alpha\,\frac{C}{L} - (1-\alpha)\frac{\sigma^2}{L},\qquad\text{then}\quad
F(\mu) \ge L(\mu)\,\Lambda(\mu),$$
by the Bochner triangle inequality $E = \|\int\psi\,d\mu\|^2 \le
(\int\|\psi\|d\mu)^2 = \sigma^2$ and Corollary 1.  Equality whenever ≤ 1
atom lies outside $\{0,1\}$ — in particular at the binding obstruction
$\mu^\star = a\delta_1 + (1-a)\delta_b$ ($a = 0.078877292705923173412$,
$b = 0.32945473850303697239$; the $p=1$ atom is a sink,
$\operatorname{rh}(1) = 0$), where $\Lambda = F/L$ exactly ($10^{-12}$
assert).  Two decisive properties: **scale-free** — as $\mu\to$ the sink set,
$\sigma^2$ is quadratic in the non-sink mass while $L$ is linear, so
$\Lambda \to 2(1-\alpha)(1-t) - 1 > 0$, and always
$\sigma^2/L \le \max_p\rho$ so $\Lambda$ is never singular; and **tight** —
the Jensen-form bound loses $1.94\times10^{-2}$ there, recorded as a failure
in the same file.  ($\max_p\rho = 0.3049467$ at $p\approx0.1165$: NUMERICAL,
sampled.)

So the certificate targets the **5-parameter** measure
$\nu = w\,\delta_{(p_1,q_1)} + (1-w)\,\delta_{(p_2,q_2)}$, $p_i\le q_i$:
$$\Phi(\nu) = \big[2(1-\alpha)(1-M) - 1\big]\,L + \alpha\,C - (1-\alpha)\,S^2,\qquad
M = \textstyle\int\frac{p+q}{2}d\nu,$$
$$L = \int\frac{h(p)+h(q)}{2}d\nu,\qquad
C = \int h(s^\star(p,q))d\nu,\qquad
S = \int\frac{\sqrt{\operatorname{rh}(p)}+\sqrt{\operatorname{rh}(q)}}{2}d\nu,$$
all linear in $\nu$.  Since $\Phi_{\mathrm{exact}} - \Phi =
(1-\alpha)(S^2 - E) \ge 0$, certifying $\Phi \ge 0$ on the 2-orbit family
certifies $\min\Phi_{\mathrm{exact}} = \inf F \ge 0$.

## 5. The $\operatorname{rh}$ lemmas and endpoint bounds

* **$\operatorname{rh}(x) \le 2x$ on $[0,\frac12]$ — PROVED,
  `lemma_rh_proof.py`.**  With $f := 2x - \operatorname{rh}$, exact identity
  (SymPy-asserted):
  $\ln2\cdot f(x) = x^2\ln(2/x) - x(2-x)\ln(1-x/2)$ — both terms strictly
  positive on $(0,\tfrac12]$.  Near 0: with
  $\rho(u) = -\ln(1-u)-u = \int_0^u \frac{t}{1-t}dt \le u^2/[2(1-U)]$
  (every algebraic step machine-checked),
  $\ln2\,f \ge x^2[\ln(2/\varepsilon) + 1 - \varepsilon/2]$,
  $\varepsilon = \frac1{16}$ — a positive multiple of $x^2$
  ($\ge 4.434\,x^2$ as printed).  On $[\frac1{16},\frac12]$: adaptive Arb
  cover (44 accepted cells, min certified bound $2.4\times10^{-3}$);
  $f(0)=0$ at the exact endpoint.
  The right half, $\operatorname{rh}(x)\le 2(1-x)$ on $[\frac12,1]$, is
  immediate from $h\le1$ — PROVED in the same file's docstring.
* **Endpoint monotonicity: $h(z) - zh'(z) = -\log_2(1-z) > 0$ — PROVED,
  `cert2.py` (`prove_rho_endpoint_monotonicity`, exact SymPy identity
  asserted at import).**  Hence $z/h(z)$ increasing, and by symmetry
  $(1-z)/h(z)$ decreasing: the $\operatorname{rh}\le 2\min$ collars give
  endpoint-attained $\rho$ caps ($\rho \le 2v/h(v)$ on
  $[u,v]\subset(0,\tfrac12]$; $\rho \le 2(1-u)/h(u)$ on the right half).
* **Global caps — PROVED:** $\rho$: `cert2._rho_global`
  ($n = 20000$; per-cell the tighter of a direct dependency-safe bound and
  the proved collar cap; observed value $\approx 0.3054$ in every campaign
  run, encoding $\sigma^2/L \le \mathrm{RHO\_GMAX}$);
  $\operatorname{rh}$: `arbcore`/`diag_exhaust` ($n = 20000$), with the
  pinned interior assertion $0.2342294 \le \mathrm{RH\_GMAX} \le 0.2350$.
* **Endpoint derivative lemma — HUMAN-AUDITED, with symbolic subchecks in `cert3.py`.**  With
  $x = 1-z$: $\ln2\,\operatorname{rh}(1-x) = (1-x)^2\ln(1-x) +
  (1-x^2)\ln(1+x)$ (import-time asserted identity), so
  $\operatorname{rh} \ge c\,x^2/\ln2$, $c = (2-3X+X^2)/2 > 0$ on $x\le X$,
  and $|d\sqrt{\operatorname{rh}}/dx| \le 1/((1-X)\sqrt{c\ln2}) =$
  `DQ_ENDPOINT` $\in (1.2, 1.4)$ at $X = 1/16$ — the $q_2\to1$ derivative
  of $\sqrt{\operatorname{rh}}$ is finite and certified.

## 6. The interval certificate (5-D branch-and-bound)

**Machine-checked reduced arithmetic.**  All certificate discharge inequalities are Arb ball
arithmetic (`flint`): ranges of $h$; $s^\star$ (monotone corners);
$\operatorname{rh}$; $\sqrt{\operatorname{rh}}$; sound min/max extensions;
exact shortest-decimal parsing of float endpoints into Arb/`Fraction`.
Every rule below is an independently sound lower bound or infeasibility
proof: skipping one can only send a box to the split branch.

**$\lambda$-shifts.**  For $\lambda\ge0$:
$\Phi_\lambda = \Phi + \lambda(M-t) \le \Phi$ on feasible ($M\le t$) points,
so any certified lower bound on $\Phi_\lambda$ over a box lower-bounds $\Phi$
on its feasible part — always sound.  The family
$\{0,1,1.626,2.2\}\cup\{0.98\lambda^\star,1.02\lambda^\star,\lambda^\star\}$
with $\lambda^\star = 1.3682158100251758$ comes from a numerical face-KKT
solve (residual $<10^{-35}$; **NUMERICAL** seed: correctness never depends
on it).

**The rule chain (`cert3.certify`), cheap-first** (trace codes in §7):
1. **Mean contractor** (`diag_exhaust.mean_contract`, PROVED): with exact
   rationals $l_i = (p_i^{\rm lo}+q_i^{\rm lo})/2$ from shortest-decimal
   endpoints, every box point has $M \ge A(w) := l_2 + w(l_1-l_2)$; solve
   $A(w)\le t$ exactly in `Fraction` arithmetic — contract $w$ or prove
   infeasibility.
2. **Corner bound** (`diag_exhaust.phi_corner`, PROVED): with
   $m_\star = \min(M_{\rm hi}, t)$, requiring
   $\kappa_\star = 2(1-\alpha)(1-m_\star)-1>0$:
   $\Phi \ge \kappa_\star L_{\rm lo} + \alpha C_{\rm lo} -
   (1-\alpha)S_{\rm hi}^2$.
3. **Ratio rule** (`cert2.ratio_rule`, PROVED): the Margin-Lemma
   Cauchy–Schwarz, $\sigma^2 = (\int\sqrt{\rho}\sqrt h\,d\mu)^2
   \le L\int\rho\,d\mu$, so
   $$\Phi \ge L\big[\kappa_t - (1-\alpha)\rho_{\rm hi}\big],\qquad
   \kappa_t = 2(1-\alpha)(1-t)-1 = 0.1918\ldots,$$
   $\rho_{\rm hi}$ the endpoint max over the box's *own* $w$-window of the
   convex combination of the two orbits' certified $\rho$ caps (the
   combination is linear in $w$; the unsound slope-form co-variation was
   caught and replaced with an endpoint max).  This handles the degenerate
   sink collars where $\Phi\to0$ and no plain sign test terminates.
4. **KKT-shifted center rules** (gated by width $\le 1/32$ and corner
   $\ge -0.02$), mean-value form over the *prepared* ordered rectangle with
   SymPy-derived interval gradients, $\Phi_\lambda \ge
   \Phi_\lambda(\mathrm{mid})^{\rm lo} - \sum_i \sup|\partial_i\Phi_\lambda
   + \lambda s_i|\,\mathrm{rad}_i$: `centered5_best`; then, when the full
   gradient is unusable at a sink, **block-mixed** `centered5_mixed` (MVT in
   $(p_1,q_1,w)$, natural interval enclosure in $(p_2,q_2)$); for
   collar-touching boxes the same after the exact orbit swap
   (`centered5_mixed_swap`), and the derivative-free-in-atoms **weight-only**
   MVT bound `centered_w_best` ($\partial_w\Phi_\lambda$ given by an
   asserted exact SymPy identity).
5. **$q_2$-pin face rule** (`cert3.q2_pin` + `face_bb`, PROVED): on boxes
   with $q_2^{\rm lo}\ge 15/16$ and $p_2^{\rm hi}\le q_2^{\rm lo}$
   (so $s^\star(p_2,q_2) = q_2$ on the whole strip), prove
   $\partial\Phi_\lambda/\partial q_2 \le 0$ on the **strip**: the negative
   terms $(\kappa v/2 + \alpha v)h'(q_2)\to-\infty$ ($\kappa>0$ certified)
   dominate positives containing only bounded aggregates and `DQ_ENDPOINT`;
   the omitted $-(1-\alpha)vL$ is $\le 0$.  Integration along $q_2\to1$
   yields $\Phi_\lambda \ge \Phi_\lambda|_{q_2=1}$; a 4-D face B&B
   (`face_bb`) then certifies $(p_1,q_1,p_2,w)$ with *full-face* bounds
   only — a feasibility-restricted candidate after pinning is **unsound**
   (a pin can map a feasible point to an infeasible one; regression frozen).
6. **Split**: influence-weighted widest coordinate
   (`cert3.split_by_influence`: atom coordinates of orbit $i$ carry at most
   its weight cap $w^{\rm hi}$ resp. $1-w^{\rm lo}$).  **Any split is
   sound** (children cover the parent exactly); the choice affects only cost.
7. **Residual**: all coordinates at the floor ($10^{-3}$ ordinary,
   $1.25\times10^{-4}$ collar-touching) — a COMPLETE slice has *zero*.

**Orbit swap — HUMAN-AUDITED** (direct termwise invariance, plus an exact
decimal-endpoint enclosure check and sampled value regression in `verify_orbit_swap`):
$(p_1,q_1,p_2,q_2,w)\mapsto(p_2,q_2,p_1,q_1,1-w)$ preserves $M,L,C,S,\Phi$,
so the root $[0,1]^4\times[\frac12,1]$ covers the full cube; the root is
split into 8 exact dyadic $w$-slices. Campaign I used min width $10^{-3}$,
collar floor $1.25\times10^{-4}$, face floor $10^{-3}$, face budget $10^5$,
box budget $2\times10^9$, 72 h/slice, 160-bit covers, and 80-bit work.

## 7. Proof trace and independent replay

**Trace — PROVED machinery, `cert3.py`:** one byte per processed DFS node
(`cert3-trace-v1`): codes 0–9 terminal outcomes (mean-infeasible,
corner-infeasible, corner, ratio, centered5, mixed, mixed-swap,
weight-only, face, residual), code $16+j$ = split of coordinate $j$.
**Replay — PROVED machinery, `cert3_replay.py`**: from the *exact dyadic
root*, it re-runs `mean_contract` itself, re-proves each claimed rule in Arb
(corner/ratio/each centered variant/pin + face sub-minimization without
deadline), follows recorded splits (midpoint strict and representable), and
requires an empty stack, no residual events, and tallies equal to the
committed record.

**Secure source-independent replay — MACHINE-VERIFIED.**
`verification/independent_arithmetic_replay_secure.py` imports no frozen
campaign arithmetic, derives centered gradients by interval AD, and
rediscovers nested face proofs.  Eight no-resume workers started at trace byte
zero and reproduced all 488,465,854 events, 244,232,923 splits, 244,232,931
leaves, exact per-rule tallies, and zero residuals.  The accepted trust boundary
still includes the trace partition, Arb primitives, published report-lock
custody, and a nonmalicious system stack.  Trusted report-lock raw SHA-256:
`af86480901c2c497739916bb410f1e117e4eaf3eefb82575ce494b7d8c04e730`;
secure composite canonical SHA-256:
`4acbd3b935bda7e51ed387e42e0598debf22f4a85b7a24ad976c3eaca4f23243`.

* **Adversarial suite — all 8 rejections VERIFIED**: flipped rule byte;
  injected residual; clear-replaced-by-split; truncation; extension; unknown
  code; false infeasibility; wrong expected tallies.  An INCOMPLETE slice
  trace replays and is rejected for exactly its pending stack.
* **Determinism**: campaigns E and F ran slice 0 five hours apart on
  different snapshots: 16,231,621 nodes each, **traces byte-identical**
  (SHA-256 `9eed0cc1382f5485…`).

## 8. Immutable campaign protocol (v2)

`cert3_par.py create` writes a fresh UUID-named campaign: the 15 pinned
sources (7 executable: `cert3.py`, `cert2.py`, `cert3_par.py`,
`diag_exhaust.py`, `bound_kkt.py`, `arbcore.py`, `entropy.py`; 8 review:
`bridge_uc.py`, `cert3_collect.py`, `cert3_replay.py`, `decomposition.py`,
`lemma_rh_proof.py`, `margin_lemma.py`, `reduction.py`, `thmB3_proof.py`)
into `snapshot/`, SHA-256-hashed, with an atomically written `launch.json`
(exact dyadic roots, UUID4 run ids, exact parameters).  Workers run from the
snapshot under `python -B`; trace → hidden temp → fsync → SHA-256+size
recorded → supervisor re-hashes → atomic install via Darwin `renamex_np`
RENAME_EXCL; nothing is ever overwritten.

`cert3_collect.py` accepts a campaign only by **re-proving it**: hardcoded
inventories (the manifest's own inventory is not trusted), pristine snapshot
listing, exact roots, unique runs, exit code 0, COMPLETE verdicts, positive
work, `stack`=`residual`=`budget_boxes`=`budget_time`$=0$, trace
digest/size consistent — then it replays **all eight traces** and only then
prints `COMPOSITE CERTIFICATE`.

## 9. Campaign record

| campaign | launched (UTC) | budget/slice | outcome |
|---|---|---|---|
| E (`…5bde6bea1eca`) | 2026-08-15 21:36 | 8 h | 6/8 COMPLETE replay-verified; slices 6,7 stopped on the wall by **12 and 18 boxes**, zero residuals |
| F (`…eb00229a822b`) | 2026-08-16 05:38 | 12 h | 7/8 COMPLETE replay-verified; slice 7: wall with **residual 1,447,470** — a heuristic *gate*, diagnosed below |
| G (`…4c2b69644e5f`) | 2026-08-16 17:56 | 24 h | **7/8 COMPLETE, zero residuals, replay-verified**; slice 7 pure compute |
| H (`…2f23a58ebdb8`) | 2026-08-17 21:20 | 24 h | **5/8 COMPLETE (0–3, 7) replay-verified**; slices 4–6 wall-only (stacks 25/27/30), zero residuals — see **A1** |

**G per slice** (committed artifacts; residual = 0 in every box of all eight):

| slice | processed | ratio clears | verdict |
|---|---|---|---|
| 0 | 15,528,149 | 192,277 | COMPLETE, replayed |
| 1 | 15,226,435 | 228,386 | COMPLETE, replayed |
| 2 | 15,381,807 | 295,807 | COMPLETE, replayed |
| 3 | 17,831,917 | 281,389 | COMPLETE, replayed |
| 4 | 19,030,155 | 82,944 | COMPLETE, replayed |
| 5 | 22,640,283 | 78,913 | COMPLETE, replayed |
| 6 | 29,889,721 | 271,760 | COMPLETE, replayed |
| 7 | 52,086,331 | 1,193,438 | INCOMPLETE — 24 h wall, `stack` 16, residual 0 |

**Diagnoses.**  *E*: cost curve: slices 6,7 reached 21.5 M and 26.1 M boxes
by the wall, both residual-free.  *F*: all 1.45 M residuals classify `sink`,
pinned to the near-total-sink slab $p_1,q_1\to(0,0)$, $q_2\to1$, $w\to1$,
$p_2 \in [0.406,0.420]$ — the degenerate family the Margin Lemma was built
for; the ratio rule cleared 400/400 dumped leaves (NUMERICAL sampling), but
the hand-tuned `ratio_plausible` gate (every atom $\le0.05$ or $\ge0.42$)
*hid the rule*: $p_2$ missed $0.42$ by $1.4\times10^{-2}$.  Gate removed
(the rule is its own cheapest test: 35.6 µs vs the 27.5 µs corner bound —
1.29×; soundness unaffected): the frozen slab certifies COMPLETE in 36,511
boxes, residual 0.  *G*: gate removed; 3.7× more ratio clears than F on
slice 0 (1.9× in the all-slice aggregate) with *smaller* trees; slice 7
passed F's residual onset at zero.  The eight G
traces total **187,614,798 processed nodes**, every box `residual`-free.
Split cost is distributed (NUMERICAL probe: 48.5 % `other`,
35.1 % mean-boundary, 16.4 % sink) — compute, not mathematics.  *H*:
splitter fix = width × influence; F's former slab: 36,511 →
**131 boxes (279×)**; obstruction root 231 → 109 (2.1×); throughput
$w\in[0.99,1]$: 2.8× (NUMERICAL).

**H final numbers** — `certH0–7`, launched 2026-08-17 21:20 UTC, 24 h/slice,
campaign dir `uc/campaigns/cert3_20260817T212001Z_4e6268eb7ce44aa3883ef15ea401c7c7_2f23a58ebdb8/`,
code SHA-256 `2f23a58ebdb8…78ce05`, launch SHA-256 `3a79d116fa924f8e…48add3`:
slices 0–3 and 7 COMPLETE — 21,135,611 / 20,827,065 / 23,486,265 /
34,186,855 / 78,257,861 processed boxes, stacks and budget flags zero,
residual 0 — with traces independently replayed (slice 7's 78,257,861-event
trace in 2 h 11 m, PASS, every tally matching the committed records);
slices 4–6 INCOMPLETE at the 24 h wall only — 54,161,482 / 59,867,104 /
74,343,239 processed, stacks 25/27/30, residual 0. All eight slices sum to
**366,265,482 processed boxes with the residual count identically zero**.
Per-slice table and evidence status: **A1**. The chain of §§1–6 is complete;
the frontier is hours, not mathematics.

## 10. Trust base

The verdict relies **explicitly** on: CPython semantics + startup hooks and
pinned-source execution (`python -B`); **flint/Arb** = the hashed interval
arithmetic itself; **mpmath, numpy, scipy, sympy** at the versions pinned in
`launch.json`; Darwin **`renamex_np`** (RENAME_EXCL) no-replace atomicity;
**SHA-256**; the hashed Python sources; and the published citations: Bauer
(1958), Winkler (1988) via Pinelis (2016), van Neerven arXiv:2112.11166v7
(Riesz / Banach–Alaoglu / Stone–Weierstrass), Bru & de Siqueira Pedra (Lemma
3.3), and standard finite Shannon entropy facts. Cambie
arXiv:2212.12500v2 supplies the source construction for the locally
reconstructed set-family bridge.  The analytic chain remains human-audited.

## References

* S. Cambie, *Better bounds for the union-closed sets conjecture using the
  entropy approach*, arXiv:2212.12500v2 — Question 2, Section 4.
* W. Sawin, arXiv:2211.11504v3, Theorem 1 (the $\psi$ sketch); Z. Chase and
  S. Lovett, arXiv:2211.10563; Alweiss–Huang–Sellke; the $\psi$ record is as
  summarized in `PROGRESS.md` Session 5 (literature list verbatim).
* Yu, arXiv:2212.00658v2 (≈ 0.38234, numerical); Liu, arXiv:2306.08824v1
  (non-explicit $c' > c^\star$, conditional); Cambie's Thm 3
  $c \approx 0.3823455$ carries the paper's "slightly less rigorous" caveat.
* J. van Neerven, *Functional Analysis*, arXiv:2112.11166v7 — Thms 2.5, 4.2, 4.50.
* H. Bauer, Arch. Math. 9 (1958) 389–393, DOI 10.1007/BF01898615.
* J.-B. Bru and W. de Siqueira Pedra, arXiv:1610.03411, Lemma 3.3.
* G. Winkler, Math. Oper. Res. 13(4) (1988) 581–587, DOI 10.1287/moor.13.4.581.
* I. Pinelis, Math. Methods Oper. Res. 83 (2016) 325–349, arXiv:1204.0249, Thm 12.

## A1. Campaign H/I evidence status

This appendix is the evidence ledger for the two campaigns sharing code hash
`2f23a58ebdb8…78ce05` (the influence-weighted splitter build).

**H's eight committed per-slice results are the comparison corpus** for the
campaign-record statements of §9 and the paper's campaign table
(`tab:campaigns`). Campaign I now also has eight supervisor-committed result
records:

| slice | H processed | H verdict | I processed (committed) | I note |
|---|---:|---|---:|---|
| 0 | 21,135,611 | COMPLETE, replayed | 21,135,611<sup>‡</sup> | trace SHA-256 identical to H |
| 1 | 20,827,065 | COMPLETE, replayed | 20,827,065<sup>‡</sup> | trace SHA-256 identical to H |
| 2 | 23,486,265 | COMPLETE, replayed | 23,486,265<sup>‡</sup> | trace SHA-256 identical to H |
| 3 | 34,186,855 | COMPLETE, replayed | 34,186,855<sup>‡</sup> | trace SHA-256 identical to H |
| 4 | 54,161,482 | INCOMPLETE at 24 h (stack 25, residual 0) | 71,405,987<sup>§</sup> | I total exceeds H wall |
| 5 | 59,867,104 | INCOMPLETE at 24 h (stack 27, residual 0) | 117,709,435<sup>§</sup> | I total exceeds H wall |
| 6 | 74,343,239 | INCOMPLETE at 24 h (stack 30, residual 0) | 121,456,775<sup>§</sup> | I total exceeds H wall |
| 7 | 78,257,861 | COMPLETE, replayed | 78,257,861<sup>‡</sup> | trace SHA-256 identical to H |
| **total** | **366,265,482** | — | **488,465,854** | **CERTIFIED 2026-08-21 — `COMPOSITE CERTIFICATE`** |

<sup>‡</sup> The H/I trace files are byte-identical by SHA-256 on slices
0–3 and 7. <sup>§</sup> The committed I total exceeds the corresponding H
wall total on slices 4–6. Campaign I's committed row therefore totals
**488,465,854 processed boxes** across all eight slices; **CERTIFIED
2026-08-21** — the frozen collector printed `COMPOSITE CERTIFICATE` after
replaying all eight traces inside its own gate.

**Determinism (SHA-256).** Campaign I
(`uc/campaigns/cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8/`,
launched 2026-08-18 21:26 UTC, 72 h/slice, workers `certI0–7`) re-ran the
identical code hash and the same exact rational target. Its slices 0–3 and 7
committed traces **byte-identical to H by SHA-256**; slices 4–6 committed
71,405,987, 117,709,435, and 121,456,775 boxes, respectively, each exceeding
the corresponding H wall total. The frozen collector replayed all eight
committed traces in its own gate and printed **`COMPOSITE CERTIFICATE`**:
**Φ ≥ 0 on the feasible 5-parameter family at certified rational
$t = 0.3820660112501052 \geq \psi + 10^{-4}$ over $w \in [1/2,1]$** (orbit-swap
symmetry covers $w \in [0,1]$). Certificate verbatim:
[`campaigns/CERTIFICATE.txt`](campaigns/CERTIFICATE.txt).

**Acceptance path.** H's replay evidence (independent PASS replays of the
five COMPLETE traces, tallies matching the committed records) stands on its
own; the §8 frozen-collector pass replayed all eight Campaign I traces on the
same bytes and printed `COMPOSITE CERTIFICATE` (2026-08-21), the sole
acceptance event for the certificate at $t = \psi + 10^{-4}$.
