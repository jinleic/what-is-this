# Perturbing a bivariate-bicycle code cannot buy distance for free

**An exact, perturbation-independent rate–distance obstruction for perturbed bivariate-bicycle (PBB) codes, with machine-checked certificates**

*Status: publication-ready draft. Every numerical claim in this file is produced by a
committed experiment in this repository and is reproducible from a clean checkout; the
provenance table in §8 gives the exact artifact path for each one.*

---

## Abstract

Perturbed bivariate-bicycle (PBB) codes are non-CSS stabilizer codes obtained from a CSS
bivariate-bicycle (BB) parent $P$ by adding a $Z$-type perturbation $[C\;D]$ to the
$X$-check block. They have been proposed as a route past the CSS BB rate–distance
envelope, on the strength of a published catalogue of $368$ codes including a non-CSS
$[[144,12,12]]$ and a $[[360,12,\le 24]]$ candidate.

We prove that the perturbation degree of freedom is not free. Let
$\Delta=\{\lambda[C\;D] : \lambda[A;B]=0\}$ be the *dressing space* of the perturbation.
We show (i) $\Delta$ is an $R$-submodule of the $Z$-sector for
$R=\mathbb{F}_2[x,y]/(x^{\ell}-1,y^{m}-1)$, (ii) the PBB dimension obeys the exact
identity $k_Q = k_P - \dim\bar\Delta$, and (iii) any minimum-weight $Z$-logical of the
parent that is *not* absorbed by $\Delta$ survives as a logical of $Q$, certifying
$d_Q \le d_Z(P)$.

Because $\Delta$ is a submodule, absorbing one minimum-weight parent logical absorbs its
entire translation orbit. This upgrades the per-perturbation criterion to a
**parent-level no-go theorem**: with
$T(P)=\dim\big(M(P)+S_Z\big)/S_Z$, where $M(P)$ is the span of the translation orbits of
*all* minimum-weight $Z$-logicals of $P$,

$$\boxed{\;d_Q > d_Z(P)\ \Longrightarrow\ k_Q \le k_P - T(P)\quad\text{for \emph{every} perturbation } [C\;D]\;}$$

$T(P)$ depends only on the parent and is computed *exactly* — not bounded — by a SAT
enumeration that terminates in a certified UNSAT. When $T(P)=k_P$ the entire family is
closed: no PBB over $P$ retains a single logical qubit while exceeding $d_Z(P)$.

We compute $T$ exactly for $139$ of the $202$ distinct parents in the published
catalogue — exhaustively for every parent at $n\le 144$. **$63$ parents are family-closed**, including the Gross code
$A=x^3+y+y^2,\;B=y^3+x+x^2$ itself ($T=12=k_P$) and all $11$ distinct catalogue
$[[144,12,12]]$ parents. This settles $253$ of the $368$ catalogue rows with no per-row
search at all, and it forecloses the headline direction of the construction over the
most important parent in the family.

The bound is not merely valid but *exact*.  A second elementary identity — the
dressing space is the image of the left kernel of $[A\;B]$, whose dimension is
$k_P/2$ on every BB parent — caps $\dim\bar\Delta\le k_P/2$.  Whenever
$T(P)\ge k_P/2$, which holds on $138$ of the $139$ certified parents, the two
bounds **sandwich**:

$$\boxed{\;d_Q > d_Z(P)\ \text{and}\ T(P)\ge k_P/2\quad\Longrightarrow\quad
k_Q = k_P - T(P)\ \text{exactly}\;}$$

Perturbation does not just pay at least $T$ logical qubits for distance; on
every parent in the theorem's class it pays *precisely* $T$, no more and no
less.  All $7$ independently certified distance-increasing perturbations sit on
the law exactly (slack $0$, hypotheses replay-certified), and a controlled
small-lattice probe — $64$ sampled parents across four lattices, $26{,}898$
valid $\delta>0$ perturbations of $Z$-support $\le 4$, distances decided by
exact symplectic meet-in-the-middle — found $496$ further
strict increases, **every one on the law, zero above it**: instances we now
know were theorem-forced rather than lucky.

On the $X$ side the companion collapse mechanism is classified completely. The
set of perturbation directions that *cannot* demote a parent $X$-stabilizer is
an invariant of the single ideal $I=\operatorname{Ann}_R(A,B)$: writing
$I^\infty$ for its stable power, $\dim S=2\dim I^\infty$, so a parent is
demote-full iff $I$ is nilpotent, demote-immune iff $I$ is idempotent, and
mixed iff neither (Theorem J-G) — reproducing the $202$-parent catalogue
classification in $9$ s with no module computation, and sharpening it: $I^2=0$
on all $192$ demoting parents, $I^2=I$ on all $10$ immune ones. Two corollaries
are structural rather than computational. If $\ell$ and $m$ are both odd, $R$ is
semisimple and *no* parent can demote at all — of the seven published BB
instances exactly one, the $[[90,8,10]]$ code on the $(15,3)$ lattice, is
protected this way. And splitting $G=\mathbb Z_\ell\times\mathbb Z_m$ into its
$2$-part and odd part shows that exact vanishing at a local factor costs two
support points per occupied $2$-coset, whence **no pair with
$\operatorname{wt}(A),\operatorname{wt}(B)\le3$ — the entire family used in the
literature — is ever mixed** (Theorem J-I), while weight $4$ already realizes
it: $A=(1+x)(y+y^2)$, $B=Ay$ on $(2,3)$ is mixed with $240$ of $255$ classes
demoting (Theorem J-J, verified by three independent routes). An exhaustive
census of $653{,}022{,}021$ weight-$\le3$ pairs over $18$ lattices confirms the
weight law with zero exceptions.

---

## 1. Setting and conventions

Let $R=\mathbb{F}_2[x,y]/(x^{\ell}-1,\,y^{m}-1)$ and let $n=2\ell m$. A polynomial
$f\in R$ acts as an $\ell m\times \ell m$ binary matrix via the commuting cyclic shifts
$x,y$. For $A,B\in R$ the CSS **bivariate-bicycle** parent is

$$H_X=[\,A\;\;B\,],\qquad H_Z=[\,B^{\mathsf T}\;\;A^{\mathsf T}\,],$$

on $n$ qubits split into two blocks of $\ell m$. Commutation is
$H_XH_Z^{\mathsf T}=AB+BA=0$, automatic since $R$ is commutative.

A **perturbed bivariate-bicycle** (PBB) code adds a $Z$-type perturbation $[C\;D]$,
$C,D\in R$, to the $X$-check block:

$$H_Q=\begin{pmatrix} A & B & \big|\ & C & D\\[2pt] 0 & 0 & \big|\ & B^{\mathsf T} & A^{\mathsf T}\end{pmatrix},$$

in the binary symplectic representation $(x\,|\,z)$ with symplectic product
$\langle (x|z),(x'|z')\rangle_s = xz'^{\mathsf T}+zx'^{\mathsf T}$. Validity requires
$M+M^{\mathsf T}=0$ with $M=AC^{\mathsf T}+BD^{\mathsf T}$, i.e. **symmetry of $M$**, not
its vanishing. (This is stated as $AC^{\mathsf T}+BD^{\mathsf T}=0$ in some summaries; the
weaker symmetric condition is the correct one and is what we verify.)

Throughout, $S_Z=\operatorname{rowspace}(H_Z)\subseteq\mathbb{F}_2^{\,n}$ is the parent's
pure-$Z$ stabilizer space, $k_P=n-\operatorname{rank}H_X-\operatorname{rank}H_Z$ is the
parent dimension, $k_Q$ the PBB dimension, and

$$d_Z(P)=\min\{\,\mathrm{wt}(z) : z\in\ker[A\;B]\setminus S_Z\,\}$$

is the parent's $Z$-distance. We write $d_Q$ for the (symplectic) distance of $Q$ and
reserve *code distance* strictly for the algebraic quantity — never a decoder estimate,
never a circuit or detector-graph distance.

**Terminology.** $C=D=0$ recovers the CSS parent. We call a catalogue row a *reversal*
when $d_Q > d_Z(P)$ is certified, and say the perturbation is *dominated* when
$d_Q \le d_Z(P)$ is certified.

---

## 2. The dressing space and the dimension identity

Define the **dressing space**

$$\Delta \;=\; \big\{\, \lambda[C\;D] \;:\; \lambda\in\mathbb{F}_2^{\,\ell m},\ \lambda[A;B]=0 \,\big\} \;\subseteq\; \mathbb{F}_2^{\,n},$$

the set of $Z$-parts of those first-block row combinations whose $X$-part cancels, and
write $\bar\Delta=(\Delta+S_Z)/S_Z$.

**Theorem G.** *For every parent $P$ and every valid perturbation $[C\;D]$:*

1. *(Centralizer invariance.)* The pure-$Z$ centralizer of $Q$ equals that of $P$, namely
   $\ker[A\;B]$.
2. *(Dimension identity.)* The pure-$Z$ stabilizer group of $Q$ is $S_Z+\Delta$, and
   $$k_Q \;=\; k_P - \dim\bar\Delta .$$
3. *(Survival criterion.)* If some minimum-weight $z\in\ker[A\;B]\setminus S_Z$ satisfies
   $z\notin S_Z+\Delta$, then $(0|z)$ is a nontrivial logical of $Q$ of weight $d_Z(P)$,
   hence $d_Q\le d_Z(P)$.

*Proof.* (i) A vector $(0|z)$ commutes with every second-block row automatically ($Z$
against $Z$). Against a first-block row $(a\,|\,c)$ the symplectic product is $a\cdot z$.
Hence $(0|z)$ is in the centralizer iff $[A\;B]z=0$ — no dependence on $[C\;D]$.

(ii) A combination $\lambda$ of first-block rows is pure $Z$ iff its $X$-part
$\lambda[A;B]$ vanishes, in which case it contributes exactly its $Z$-part
$\lambda[C\;D]$. Adding the second block gives pure-$Z$ stabilizer space $S_Z+\Delta$.
For the dimension, project the row space onto the $X$-coordinate: the image is
$\operatorname{rowspace}[A\;B]$ for both codes, while the kernel is $\{0\}\times S_Z$
for $H_P$ and $\{0\}\times(S_Z+\Delta)$ for $H_Q$ — the pure-$Z$ stabilizer space,
by the first part of (ii). Hence
$\operatorname{rank}H_Q=\operatorname{rank}H_P+\dim\bar\Delta$, and $k=n-\operatorname{rank}H$
gives $k_Q=k_P-\dim\bar\Delta$. (Note that $\operatorname{rowspace}H_P$ is **not**
contained in $\operatorname{rowspace}H_Q$ — parent $X$-stabilizers are demoted —
so the rank identity, not a containment, is what makes the count work.)
Equivalently, via the sector identity valid for every
stabilizer code, $\dim(\mathrm{Xcen}/S_X)=\dim(\mathrm{Zcen}/S_Z)=k$: the unchanged
pure-$Z$ centralizer (i) minus the enlarged stabilizers $S_Z+\Delta$ is directly
$k_P-\dim\bar\Delta$. (The shortcut "the $X$-sector of the quotient is unchanged, so only
the pure-$Z$ quotient loses dimension" is *false as an argument*: the pure-$X$
centralizer also shrinks, $\mathrm{Xcen}(Q)=\mathrm{Xcen}(P)\cap\ker[C\;D]$ — e.g. from
dimension $40$ to $12$ on `phase2_58` — and both the rank identity and the sector
identity are what make the count work; `notes/theorem_j_xsector.md`, J.1–J.3.)

(iii) Immediate from (i) and (ii): $z\in\ker[A\;B]$ by (i), and $z\notin S_Z+\Delta$
means its class is nontrivial in the quotient of (ii). $\square$

Part (iii) is the operational content: **a single unabsorbed minimum-weight parent
logical caps the PBB distance**, with no search over PBB codewords. Machine checks of
(i)–(iii) on catalogue instances are in `tests/test_pbb_survival.py`; the dimension
identity (ii) holds exactly on all tested catalogue rows (§8).

---

## 3. From per-perturbation to parent-level: the module structure

The step that makes a family-wide theorem possible is the following.

**Lemma 1 (submodule).** *$\Delta$ is an $R$-submodule of $\mathbb{F}_2^{\,n}$; so is
$S_Z$, and hence so is $S_Z+\Delta$.*

*Proof.* The coefficient set $L=\{\lambda : \lambda[A;B]=0\}$ is an ideal of $R$: if
$\lambda A=\lambda B=0$ then $(x^ay^b\lambda)A = x^ay^b(\lambda A)=0$ and likewise for
$B$, using commutativity of $R$. Since $\Delta=L\cdot[C\;D]$ and multiplication by
$[C\;D]$ is $R$-linear, $\Delta$ is closed under the translations $x^ay^b$, i.e. is a
submodule. $S_Z=\operatorname{rowspace}[B^{\mathsf T}\;A^{\mathsf T}]$ is a submodule for
the same reason. $\square$

**Corollary 2 (orbit absorption).** If $\Delta$ absorbs a minimum-weight $Z$-logical $z$
(i.e. $z\in S_Z+\Delta$) then it absorbs the entire translation orbit
$\{x^ay^b z\}$ of $z$. Translations preserve weight and map $\ker[A\;B]$ to itself, so
every orbit element is again a minimum-weight $Z$-logical.

Now define, **from the parent alone**,

$$M(P)=\operatorname{span}\Big(\bigcup\nolimits_{z}\ \mathrm{orbit}(z)\ :\ z \text{ a minimum-weight } Z\text{-logical of } P\Big),\qquad T(P)=\dim\frac{M(P)+S_Z}{S_Z}.$$

**Theorem H (parent-level no-go).** *For every parent $P$ and every valid perturbation
$[C\;D]$,*

$$d_Q > d_Z(P)\quad\Longrightarrow\quad k_Q \;\le\; k_P - T(P).$$

*In particular, if $T(P)=k_P$ then no perturbation of $P$ retains a logical qubit while
exceeding $d_Z(P)$: the whole family over $P$ is closed.*

*Proof.* Suppose $d_Q>d_Z(P)$. By Theorem G(iii), contrapositively, **every**
minimum-weight $Z$-logical $z$ of $P$ must be absorbed: $z\in S_Z+\Delta$. By
Corollary 2 the whole orbit of each such $z$ is absorbed, so
$M(P)\subseteq S_Z+\Delta$, whence $\bar M\subseteq\bar\Delta$ and
$\dim\bar\Delta\ge T(P)$. Theorem G(ii) gives
$k_Q=k_P-\dim\bar\Delta\le k_P-T(P)$. The final sentence is the case $T=k_P$, forcing
$k_Q\le 0$. $\square$

Three features are worth emphasising. First, the hypothesis is *only* a distance
increase — no assumption on $C,D$, on weight, or on locality. Second, $T(P)$ is a
property of the parent, so one computation covers the infinite family of perturbations
over $P$, including every $[C\;D]$ nobody has enumerated. Third, the conclusion bounds
$k_Q$, not $d_Q$: the theorem does not forbid distance increase, it *prices* it.

### 3.1 Computing $T$ exactly

$T$ is obtained by a monotone SAT enumeration rather than estimated. Maintain
$W=S_Z+M$, initially $S_Z$, and repeatedly ask a solver for
$z\in\ker[A\;B]$ with $\mathrm{wt}(z)\le d_Z(P)$ and $z\notin W$. Nontriviality with
respect to $W$ is encoded exactly as ordinary logical nontriviality — pair $z$ against a
basis of $W^{\perp}$, since $\exists f\in W^{\perp}: f\cdot z=1$ iff $z\notin W$. Each
SAT hit is necessarily a *minimum-weight* logical (nothing outside $S_Z$ has weight below
$d_Z$), and its full orbit is adjoined to $M$. The enumeration terminates in UNSAT, which
certifies that $W$ contains every minimum-weight $Z$-logical; therefore the resulting $T$
is **exact**, not a lower bound. If the budget is exhausted first, the partial $T$ is a
certified *lower* bound — the sound direction, since Theorem H only needs
$\dim\bar\Delta\ge T$.

Two soundness details are enforced in code. (a) Every SAT witness is re-verified through
two independent GF(2) paths and its weight is checked to equal the supplied
$d_Z(P)$ exactly; a mismatch is raised, never absorbed. (b) A translation symmetry break
(anchor support at index $0$ of one block) is applied to the decision queries. This is
satisfiability-preserving precisely because $\ker[A\;B]$, $S_Z$ and every $M$ built from
full orbits are translation invariant — the solution set of each query is translation
invariant, so any solution can be translated to the anchored form. The premise is
machine-validated in `tests/test_translation_symmetry_break.py` (automorphism action,
transitivity, and $30$ differential SAT/UNSAT agreements); it buys $6.2\times$ on the
decisive UNSAT proofs ($77.8$ s $\to$ $12.5$ s on the $[[144,12,12]]$ cap-$11$ instance).

### 3.2 The trade is exact: forced saturation

Theorem H bounds $k_Q$ from above. A second, independent observation bounds
$\dim\bar\Delta$ from above — and for most parents the two meet.

**Lemma 2 (left-kernel identity).** *For every CSS BB parent, with
$L=\{\lambda:\lambda A=\lambda B=0\}$ the left kernel of $[A\;B]$,
$k_P = 2\dim L$.*

*Proof.* $\operatorname{rank}H_X=\ell m-\dim L$ ($L$ is the left nullspace of
$H_X$). The left nullspace of $H_Z=[B^{\mathsf T}\,A^{\mathsf T}]$ is the
transpose of $K=\ker A\cap\ker B$, so $\operatorname{rank}H_Z=\ell m-\dim K$,
and $k_P=\dim L+\dim K$. The coordinate-reversal map conjugates each cyclic
shift to its transpose, so $\dim L=\dim K$. $\square$

*(Verified on all $202$ distinct parents, `tests/test_pbb_theorems.py`.)*

**Lemma 3 (dressing ceiling).** *$\Delta=L\cdot[C\;D]$ is the image of a
linear map on $L$, so $\dim\bar\Delta\le\dim L=k_P/2$ for every valid
perturbation.* *(Verified on all $368$ catalogue rows via
$\dim\bar\Delta=k_P-k_Q$.)*

> **Theorem I (forced saturation).** *If $T(P)\ge k_P/2$ and $d_Q>d_Z(P)$, then
> $\dim\bar\Delta=T(P)$ and $k_Q=k_P-T(P)$ **exactly**.*

*Proof.* Theorem H gives $\dim\bar\Delta\ge T$; Lemma 3 gives
$\dim\bar\Delta\le k_P/2\le T$. Squeeze, then apply Theorem G(ii). $\square$

**Scope.** $T\ge k_P/2$ holds on $138$ of $139$ certified parents. The single
exception is `9a7638586033` ($n=144$, $k_P=12$, $T=4$, not family-closed),
where increase is only sandwiched to $4\le\dim\bar\Delta\le 6$. On every other
certified parent the exact trade law is a theorem: a distance increase pays
precisely $T(P)$ logical qubits. This subsumes the empirical saturation of §5.2
— all $7$ reversals have $T=k_P/2$ — and, retroactively, explains the
small-lattice probe of §7: its $64$ parents all have $T\in\{k_P/2,k_P\}$, so
its $496$ strict increases could not have landed anywhere but on the law.

---

## 4. The Gross code: its entire perturbation family is closed

The canonical Gross code (Bravyi et al., arXiv:2308.07915) is the BB code with
$\ell=12$, $m=6$, $A=x^3+y+y^2$, $B=y^3+x+x^2$, giving $[[144,12,12]]$. Running the
enumeration of §3.1 on it:

| parent | catalogue row | $\ell,m$ | $A$ | $B$ | $k_P$ | $d_Z(P)$ | $T$ | $T$ exact | family closed |
|---|---|---|---|---|---|---|---|---|---|
| `aebb649586c1f5ab` | `12_6_0194` | $12,6$ | $x^3+y+y^2$ | $y^3+x+x^2$ | $12$ | $12$ | $\mathbf{12}$ | yes (UNSAT) | **yes** |

$T=k_P=12$, so by Theorem H:

> **Corollary 3 (Gross closure).** No perturbed bivariate-bicycle code built on the Gross
> code retains a single logical qubit while exceeding $Z$-distance $12$. The entire
> perturbation family over the Gross code — not merely the catalogued members — is closed.

The same holds for **all $11$ distinct $[[144,12,12]]$ parents present in the published
catalogue**, each with $T=12=k_P$ certified by a terminating UNSAT, including
`4c7eb964a6a3e38c`, the parent of the benchmark non-CSS PBB $[[144,12,12]]$
(`12_6_0193`, $220.4$ s, $3$ witnesses). That benchmark code is independently interesting
— we certify it is *genuinely* non-CSS (not CSS after row operations, any qubit
permutation, or any local Hadamard; a $19$-parity-mask plus one-hot certificate is
CP-SAT infeasible) and that its exact distance is $d=12$, two-sided: weights $\le 5$
excluded by exact meet-in-the-middle, weight $6$ by a complete classification (all
$72/72$ zero-syndrome vectors coincide with the stored pure-$Z$ checks, so no weight-$6$
logical exists), weights $7$–$11$ by four CaDiCaL UNSAT sector proofs
($502/674/422/576$ s), and a weight-$12$ witness verified through two independent GF(2)
paths. It matches the Gross code's $[[144,12,12]]$
parameters; Corollary 3 shows that no sibling perturbation over the same parent can
improve on them without giving up all $12$ logical qubits.

---

## 5. Catalogue-scale results

The published catalogue (arXiv:2606.02418) contains $368$ PBB rows over $202$ distinct
CSS BB parents. Applying §3.1 to the parents:

| quantity | value |
|---|---|
| distinct parents | $202$ |
| parents with **exact** $T$ (terminating UNSAT) | $139$ |
| parents **family-closed** ($T=k_P$) | $\mathbf{63}$ |
| catalogue rows capped *a priori* ($k_Q>k_P-T$, no per-row work) | $\mathbf{253}$ of $368$ |
| catalogue rows Theorem H permits, of the $284$ classified (candidates for a genuine increase) | $31$ |

The per-length breakdown matters more than the totals, because the enumeration is
**complete at every length up to $n=144$**:

| $n$ | parents | certified | family-closed | rows | capped *a priori* |
|---|---|---|---|---|---|
| $36$ | $3$ | $3/3$ | $3$ | $47$ | $47$ |
| $72$ | $8$ | $8/8$ | $5$ | $34$ | $32$ |
| $108$ | $68$ | $68/68$ | $28$ | $107$ | $80$ |
| $144$ | $38$ | $38/38$ | $21$ | $66$ | $65$ |
| **$\le 144$** | **$117$** | **$117/117$** | **$57$** | **$254$** | **$224$** |
| $180$ | $46$ | $20/46$ | $5$ | $64$ | $27$ |
| $360$ | $39$ | $2/39$ | $1$ | $50$ | $2$ |

So over the entire catalogue at $n\le 144$ — which contains every $[[144,12,12]]$
candidate and both $[[72,\cdot,\cdot]]$ reversals — the classification is exhaustive: all
$117$ parents have an exact $T$, $57$ of them are family-closed, and $224$ of $254$ rows
are settled with no per-row search. The $63$ uncertified parents are exactly the $n=180$
and $n=360$ ones, where a single terminating UNSAT can cost hours (the slowest completed
parent took $16{,}493$ s); that sweep continues, and because $T$ is monotone under further
enumeration these figures are lower bounds on the final closure count, never upper bounds.

### 5.1 A structural confinement

Cross-tabulating the $139$ certified parents against their $284$ catalogue rows exposes a
sharp pattern. Only four ratios $k_Q/k_P$ occur, and only three ratios $T/k_P$:

| $k_Q/k_P$ | $T/k_P=1$ | $T/k_P=1/2$ | $T/k_P=1/3$ | capped *a priori* |
|---|---|---|---|---|
| $1$ | $109$ | $69$ | $1$ | yes |
| $5/6$ | $11$ | $10$ | — | yes |
| $3/4$ | $25$ | $13$ | — | yes |
| $1/2$ | $15$ | $\mathbf{31}$ | — | $15$ yes / $\mathbf{31}$ **no** |

Every row that keeps more than half of the parent's logical qubits is capped
*a priori*: $d_Q\le d_Z(P)$ with no search. The only rows Theorem H permits to increase
distance are the $k$-halving ones with $T/k_P=1/2$ — exactly $31$ of $284$. The observed
distance-increasing perturbations are confined to a $k$-halving class, and this is now a
proved confinement rather than an empirical regularity.

### 5.2 The bound is saturated, not merely valid

The catalogue contains $7$ perturbations for which we hold a certified reversal
$d_Q>d_Z(P)$. Theorem H must permit every one of them or it is refuted. A refutation
requires *both* halves — hypothesis and violated conclusion — so we verify the hypothesis
independently by a replayed UNSAT lower bound on $d_Q$, and count any check whose
hypothesis is not independently established as vacuous rather than as support.

| row | $n$ | $k_P$ | $k_Q$ | $d_Z(P)$ | certified $d_Q>$ | $T$ | ceiling $k_P-T$ | slack | hypothesis |
|---|---|---|---|---|---|---|---|---|---|
| `12_6_0217` | $144$ | $8$ | $4$ | $8$ | $9$ | $4$ | $4$ | $\mathbf{0}$ | certified |
| `phase2_58` | $72$ | $8$ | $4$ | $4$ | $4$ | $4$ | $4$ | $\mathbf{0}$ | certified |
| `phase2_60` | $72$ | $8$ | $4$ | $4$ | $4$ | $4$ | $4$ | $\mathbf{0}$ | certified |
| `phase2_71` | $108$ | $12$ | $6$ | $4$ | $5$ | $6$ | $6$ | $\mathbf{0}$ | certified |
| `phase2_72` | $108$ | $12$ | $6$ | $4$ | $5$ | $6$ | $6$ | $\mathbf{0}$ | certified |
| `phase2_88` | $108$ | $4$ | $2$ | $4$ | $5$ | $2$ | $2$ | $\mathbf{0}$ | certified |
| `9_6_0183` | $108$ | $4$ | $2$ | $8$ | $9$ | $2$ | $2$ | $\mathbf{0}$ | certified |

$7/7$ checked, $7/7$ non-vacuous, $0$ violations, and **slack exactly $0$ in every
case**: $k_Q = k_P - T(P)$ identically. Every certified distance increase pays exactly
the price Theorem H charges, never less and never more — and Theorem I (§3.2)
proves they never could: all seven have $T=k_P/2$, inside the theorem's class, so
equality was forced before it was observed.

> **Observation 4 (saturation — now a theorem for $T\ge k_P/2$).** On every
> certified reversal in the published catalogue, $k_Q = k_P - T(P)$. Theorem I
> (§3.2) proves this is forced for every parent with $T\ge k_P/2$ — all
> certified parents except `9a7638586033`; §7 carries the residual.

The gate is executable: `experiments/exp039_nogo_module.py gate` re-derives the table and
exits nonzero on any violation, and is locked by `tests/test_pbb_nogo.py`.

---

## 6. What this does and does not settle

**Settled.** Over the $63$ family-closed parents — including the Gross code and every
catalogue $[[144,12,12]]$ parent — the perturbation degree of freedom cannot produce a
code that beats the parent's $Z$-distance while retaining a logical qubit. This is a
theorem about all $[C\;D]$, not a search result, and it is certified by terminating UNSAT
proofs with hash-bound, replayable CNFs.

**Priced, not forbidden.** For the remaining parents, distance increase is possible but
pays at least $T(P)$ logical qubits, and — by Theorem I — pays *exactly* $T(P)$ on
every parent with $T\ge k_P/2$, which is every certified parent but one. It would be
wrong to claim the trade is always unprofitable *relative to the parent*: on five of the
seven certified reversals it is mildly profitable in $kd^2/n$, because $+2$ distance at
$d_Z(P)=4$ outweighs halving $k$ (e.g. `phase2_71`: parent $12\cdot 4^2/108=1.78$ versus
PBB $6\cdot 6^2/108=2.00$). The correct comparison is not the parent but the CSS envelope
at the same length, and there the reversals lose decisively — every one of the seven is
dominated at equal $n$ by a CSS BB code with $k_C\ge k_Q$ and certified exact
$d_C\ge$ the PBB's upper bound (all seven machine-checked in `results/processed/exp036_envelope_check.json`, schema
`exp036-envelope-check-v2`, against $162$ certified-exact CSS candidates; the dominator
shown is the strongest at that length by $kd^2$):

| reversal | $[[n,k_Q,d_Q{\le}]]$ | $kd^2/n$ | dominating CSS | $kd^2/n$ | factor |
|---|---|---|---|---|---|
| `12_6_0217` | $[[144,4,10]]$ | $2.78$ | $[[144,12,12]]$ | $12.00$ | $4.3\times$ |
| `9_6_0183` | $[[108,2,10]]$ | $1.85$ | $[[108,8,10]]$ | $7.41$ | $4.0\times$ |
| `phase2_71`, `phase2_72` | $[[108,6,6]]$ | $2.00$ | $[[108,8,10]]$ | $7.41$ | $3.7\times$ |
| `phase2_88` | $[[108,2,6]]$ | $0.67$ | $[[108,8,10]]$ | $7.41$ | $11.1\times$ |
| `phase2_58`, `phase2_60` | $[[72,4,6]]$ | $2.00$ | $[[72,12,6]]$ | $6.00$ | $3.0\times$ |

So the perturbation can beat *its own parent*, but never the best CSS code of the same
length: the reversals trade a $3.0$–$11.1\times$ deficit against the envelope for a $\le 1.13\times$
gain against the parent.

**The $X$ sector.** The full analysis is in `notes/theorem_j_xsector.md` (J.0–J.7) and
`notes/j5_attempt.md`; what the main results need is: $d_X(P)=d_Z(P)$ for every BB parent
(reversal symmetry, J.0), so the $Z$-side framing loses no generality; the pure-$X$
centralizer of $Q$ shrinks by
$\rho_X=\dim((R_{CD}+S_Z)/S_Z)\ge\dim\bar\Delta$ (J.1); and $Q$'s pure-$X$ logicals are
either survivors of weight $\ge d_X(P)$ or demoted parent stabilizers paired with
parent $Z$-logicals (J.4). Two further statements are now *decided*, one in each
direction:

- **$Z$-distance is monotone (L1, proved).** The pure-$Z$ centralizer is unchanged
  (Theorem G(i)) and $S_Z(Q)=S_Z+\Delta\supseteq S_Z$ grows, so every parent $Z$-logical
  class either survives or is absorbed — but any *surviving representative* is a parent
  $Z$-logical itself. Hence $d_Z(Q)\ge d_Z(P)$ for **every** valid perturbation. This is
  the dual counterpart of the caps and requires no hypotheses.
- **$X$-distance can collapse (J.5 refuted, machine-verified).** The candidate
  monotonicity $d_X(Q)\ge d_X(P)$ is **false**: there exist valid perturbations that
  demote a low-weight parent $X$-stabilizer to a genuine $X$-logical of $Q$ with *no
  loss of logical qubits* ($\dim\bar\Delta=0$). Verified witnesses on four parents with
  independently certified exact parent distances, each confirmed on two independent
  GF(2) code paths ($x_0\in\ker H_Z\cap\ker[C\;D]\cap\operatorname{rowspace}H_X$,
  $x_0\notin S_X(Q)$; $z_0=e_0[C\;D]\in\ker H_X$, $z_0\notin S_Z+\Delta$):
  all four with $x_0$ of weight $6$ —
  `6_6_0099` ($n{=}72$: $\le6<8$), `9_6_0136` ($n{=}108$: $\le6<12$),
  `12_6_0193` ($n{=}144$: $\le6<12$, *the same parent as the headline non-CSS code*),
  `15_6_0219` ($n{=}180$: $\le6<14$).
  The mechanism is the **syzygy family** $\{(C,D):AC^{\mathsf T}+BD^{\mathsf T}=0\}$:
  its relative density inside the validity space is exponentially small
  ($\approx2^{-34}$ on the flagship parent — the reason the earlier uniform hunts
  returned zero violations), a single-row demotion exists **only if** $M=0$ (Lemma L3),
  its nontrivial excess over the row-removable subfamily has dimension exactly $k_P$
  on every parent tested (Lemma L4), and no $Z$-side mirror argument exists (L5).
  Earlier evidence stands with corrected interpretation: $0$ collapses in $6.47$M
  *uniformly sampled* valid perturbations is a concentration statement about the
  complement of the syzygy class, not about the class itself.
  **The collapse is ubiquitous, not exceptional**: a syzygy-targeted census over $102$
  parsed parents ($n\in\{144,180,360\}$) finds a collapse in **every bounded parent
  checked ($49/49$: all $28$ at $n=144$, all $20$ bounded at $n=180$, the one bounded
  at $n=360$)**; every demoted logical found has weight exactly $6$; the $53$
  violation-free rows all lack a distance bound and score nothing
  (`exp044_n180_hunt.json`).
  **The mechanism is computable in polynomial time** (EXP-046, Theorem J-C of
  `notes/theorem_jc_syzygy_constructor.md`): the collapse family is the annihilator
  $c\in\ker A$ with $D=0$ — validity $M=AC^{\mathsf T}=0$ and partner membership
  $z_0=(c,0)\in\ker H_X$ are both *automatic*, so the demotion decision is a single
  GF(2) rank test with $\dim\ker A\ge k_P/2$; the deterministic nullspace constructor
  demotes on **192 of all 202 catalogue parents** (100 % from single basis vectors)
  and collapses **every bounded open-window parent but one (99/100: $9_6_0175$,
  certified channel-free below its distance by EXP-051: the one-shot CP-SAT query
  — rowspace membership at weight $\le d_P-1$ beside every original translated row
  ruled out — is INFEASIBLE, proving multi-row light channels below $8$ absent,
  and its single-row syzygy channel is immune by Theorem J-E exactly)**;
  **EXP-048 upgrades those collapse witnesses to exact child distances**:
  for every one of the $133$ distance-bounded demotion parents it rebuilds the
  exact $D=0$ sibling, re-verifies the parent fingerprint, $M=AC^{\mathsf T}=0$,
  $c\in\ker A$, and the pure-$X$ weight-$6$ witness, and decides the full
  symplectic distance below $6$ with \textsc{cadical195}; $133/133$ end
  terminal-UNSAT after bounded descent (five require a second call at cap $3$).
  The exact child distances are $128$ at $6$ and $5$ at $4$; the strict exact
  parent drops are $99$, and **$98/99$ are exact strict drops at preserved
  dimension** ($\bar\Delta=0$), the sole exception \texttt{30\_6\_0287}
  ($k=12\to10$, exact child distance $6$ against parent bound $12$).
  **Theorem J-E upgrades the mechanism to an exact characterization**
  (`notes/theorem_je_exact_decision.md`, proved three ways incl.\ Nakayama's lemma,
  machine-confirmed 202/202): with $M:=\ker H_X/S_Z$ as an $R$-module and
  $L_{\rm pre}\subset R$ the left null space of $[A\;B]$,
  *no single-row syzygy demotion exists $\iff 1\in L_{\rm pre}+\operatorname{Ann}_R(M)$*
  — one GF(2) rank test decides collapse-hazard per parent; exhaustive class
  enumeration agrees everywhere (zero upgrades/regressions); the demote fixed
  set is structurally the stable image $I^\infty M$ (Theorem J-E′′, proved:
  fixed classes are exactly the idempotent part on $I$-full local factors;
  immune $\iff IM=M$, ubiquitous $\iff I$ nilpotent on $M$, mixed otherwise;
  for $k_P>0$ the three cases are mutually exclusive), computed in EXP-047's
  quotient coordinates by iterating $F_{r+1}=IF_r$ with descent-asserted
  termination (EXP-052): **on all $202$
  catalogue parents the demote-fraction is exactly $1$ or exactly $0$
  ($192/10$, zero mixed), with chain length $\le 3$, $I^2M=0$ on every
  demoting parent, and the immune set reproducing EXP-047's ten exactly**;
  an independent exact enumeration over the $196$ parents with $k_P\le 20$
  (2{,}693{,}100 nonzero classes, EXP-050) agrees with the algebraic
  predictions to the integer on every parent, and the six
  $k_P\in\{24,40,60\}$ parents are all demote-full through the algebraic
  route. The trichotomy is now *fully classified*: it is an invariant of the
  single ideal $I=\operatorname{Ann}_R(A,B)$, mixed parents exist but never at
  weight $\le3$, and problem J-A is resolved — see §6.1,
  while on the immune side
  **`9_6_0175` is CERTIFIED X-monotone over its entire valid family**
  (syzygy channel exactly empty — Theorem J-E; multi-row light exclusion below
  $8$ proved by an exhaustive one-shot CP-SAT channel query, EXP-051:
  dual-kernel membership encoding, weight $\le 7$, every original row forbidden
  by an exclusion clause, INFEASIBLE; the same certificate holds for
  `phase2_75/76/77/83/84/87` and `phase2_109` below their respective
  distances $6$). The remaining two
  immune parents (`15_6_0256`, `30_6_0289`) each carry an independently
  verified weight-$8$ multi-row light channel, so their $X$-monotonicity is
  organically undecided rather than merely uncertified.
  The X-side of the catalogue is
  thereby classified by mechanism: 192 demotion-realized (99 with exact
  strict-collapse certificates from EXP-048; the rest at certified parity,
  vacuous UB, or awaiting the parent's own distance) / 8 **certified** X-monotone
  (`9_6_0175` and `phase2_75/76/77/83/84/87/109`, exhaustive EXP-051 certificates)
  / 2 X-U with verified weight-8 light channels (`15_6_0256`, `30_6_0289`,
  monotonicity not decidable by exclusion),
  for $192+8+2=202$. Practical reading: a PBB candidate can be *silently* low-distance
  on the $X$ side (the flagship's sibling is exactly $d_Q=d_X=6$ at full $k$,
  certified by EXP-048) — any candidate
  used downstream needs its own two-sided distance certificate in the fashion of
  Theorem F; and the audit is one GF(2) rank test away (Theorems J-C/J-E).

The caps in Theorems G/H/I do not use any of this: they cap
$d_Q$ through pure-$Z$ survivors only, and L1 is consistent with them.

**The residual class is small.** The known violation of $T\ge k_P/2$ might still be the
tip of a large uncharted class. It is not, at least at small lattice size: among
$10{,}645$ distinct connected $k_P\ge2$ parents with $A,B$ each of weight $\le3$
(support families $2{\times}2$, $2{\times}3$, $3{\times}3$), exhaustive over all $43$
ordered lattices with $\ell\cdot m\le 40$ and quotiented by independent block
translations with $\{A,B\}$ unordered ($4{,}621$ parents for $\ell\cdot m\le24$ plus
$6{,}024$ for $24<\ell\cdot m\le40$; lattices $6{\times}7$ and larger skipped by an
extrapolated budget model), **none** has $T<k_P/2$ — in every case computed,
$T\in\{k_P/2,\,k_P\}$ exactly. The only known exception (`9a7638586033`, $k_P{=}12$,
$T{=}4$) sits at $\ell\cdot m=72$, outside this scope. Whether the residual class
(parents with $T<k_P/2$) is empty beyond `9a7638586033` at larger lattices is under
computation; for the known witness the residual question is *decided* (§7, EXP-045):
no perturbation of it can increase distance at all.

**Not settled by this paper.** (i) Whether *every* $\delta>0$ PBB is dominated by *some*
CSS code — the universal envelope statement remains open (EXP-037 classifies
only the 368-row published catalogue). **Within the catalogue, every row through
$n=180$ is dominated**: $254/254$ dominated at $n\le144$ and $64/64$ dominated
at $n=180$; only the $n=360$ catalogue rows remain unresolved ($28$ candidate
escapees with $22$ pending). The former $n=180$
holdouts are all decided adversely to PBB: `15_6_0224` and `15_6_0229` are
dominated by the same certified CSS $[[180,8,16]]$
(\texttt{855e8bf135ba4e05}); the replacement candidates
(`bliss_e02255ad9ff31483`, `bliss_b2215e94b5b2b50f`, $k=6$) were re-tightened
from UB $18$ to $14$ and $16$, hence also dominated; and the final candidate
`phase2_107` — an exact PBB $[[180,20,6]]$ — is matched \emph{with equality}
by its own CSS parent certified exact $[[180,20,6]]$
(\texttt{15654a1cebe3e639}), so it too is dominated. The $28$ candidates plus
the $22$ pending $n=360$ rows are coverage artifacts of the incomplete $n=360$
CSS census (the running EXP-039 sweep). (ii) The $n=360$ parents remain open
($38$ of $39$ not yet certified). (iii) Everything here is
a *code-level* statement. Circuit-level cost is a separate axis: our detector-error-model
mechanism-distance intervals for the two $[[144,12,12]]$ circuits are tied at the
certified interval $[5,12]$ on each — weight-$4$ exclusion completed on both circuits
(perfect-matching class by meet-in-the-middle, star class by exact cover; no logical
failure from any set of $\le4$ DEM mechanisms; structural-probability DEM), with
weight-$12$ data-logical witnesses as upper bounds (EXP-029/EXP-040) — the PBB
schedule is $\approx 6$–$17\times$ slower on an idle machine ($p_{50}$ $10.3\times$,
$p_{95}$ $16.6\times$, $p_{99}$ $6.3\times$), and no circuit-level result has favoured the
PBB candidate. **Provisional matched-cell Monte Carlo (EXP-040, 2026-08-19):** the
first completed matched schedule pair at $p=0.002$ circuit-level noise separates
clearly — CSS Gross sched#4382: block LER $6.33\times10^{-3}$ (95\% simultaneous CI
$(5.39,7.37)\times10^{-3}$, $52{,}000$ shots) vs PBB sched#4107: $2.21\times10^{-2}$
(CI $(2.112,2.311)\times10^{-2}$, $172{,}000$ shots) — disjoint intervals, the PBB
$\approx3.5\times$ higher; cumulative decode-core time per shot is $\approx10\times$
higher for the PBB (throughput proxy only; pooled latency quantiles are *not*
claimed). This is **one finished cell of a 10-cell grid** — a provisional
single-pair data point, not a verdict; the grid-complete matched comparison is
recorded only in the assembled EXP-040 artefact once all cells certify.
A verified one-ancilla mixed-stabilizer syndrome-extraction circuit for
the non-CSS PBB, benchmarked against the Gross code under identical noise, remains the
decisive end-to-end experiment.

**Honest limitations.** $T$ is exact only where the enumeration terminated in UNSAT
($139/202$ parents; all $117$ at $n\le 144$); elsewhere it is a certified lower bound, which is the sound
direction. $d_Z(P)$ is used as supplied and is independently certified for the parents
in the table. Catalogue-reported distances are treated as untrusted throughout; where we
quote a distance we quote our own two-sided certificate, and one catalogue heuristic
bound was found to be $\ge 2.5\times$ loose ($d\le 40$ advertised, $\le 16$ re-verified).

### 6.1 The demote law is an ideal invariant, and no trinomial pair is mixed

The single-row demote mechanism above is decided by a *module* computation in
§6 ($S=I^\infty M$ inside $M=\ker H_X/S_Z$). It is in fact governed by the
ideal alone. Write $R=\prod_\beta R_\beta$ into local Artinian factors and
recall $I=L_{\rm pre}=\operatorname{Ann}_R((a,b))$. This is *not* the module
annihilator $\operatorname{Ann}_R(M)$ that also appears in Theorem J-E, and the
two cannot be substituted for one another: $I$ is full exactly on the factors
where $a_\beta=b_\beta=0$ and zero where either polynomial is a unit, whereas
$\operatorname{Ann}_R(M)=(\bar a,\bar b)$ is full exactly on the unit factors —
dual to $I$, with $\dim I+\dim\operatorname{Ann}_R(M)=\ell m$.

**Theorem J-G (ideal invariance; proved,
`notes/theorem_jg_ideal_invariant.md`).** $\dim S=2\dim_{\mathbb F_2}I^\infty$,
and for $k_P>0$
$$\textbf{demote-full}\iff I\ \text{nilpotent},\qquad
\textbf{immune}\iff I^2=I,\qquad
\textbf{mixed}\iff 0\neq I^\infty\neq I,$$
with $M=0\iff I=0$. *Proof sketch.* $I=\prod I_\beta$ and
$M=\bigoplus M_\beta$ split componentwise; $I_\beta=R_\beta$ exactly on the
factors where $A_\beta=B_\beta=0$, and there the localized $H_X$ is zero, so
$M_\beta=R_\beta^2$; $I^\infty=e_IR$ because proper $I_\beta\subseteq\mathfrak
m_\beta$ is nilpotent. Theorem J-E″ then gives
$S=e_IM=\bigoplus_{\rm full}M_\beta$ of dimension $2\dim e_IR$. The immune case
uses that $R=\mathbb F_2[G]$ and each of its blocks is Frobenius, whence
$\dim I_\beta=\dim R_\beta/(A_\beta,B_\beta)=\tfrac12\dim M_\beta$, so
$M_\beta=0\iff I_\beta=0$. $\square$

Consequences, all machine-checked (EXP-053 classifies the whole catalogue in
$8.9$ s — the module chains of EXP-052 are no longer needed; the local-factor
audit below adds $\approx4$ minutes):

* the ideal route reproduces EXP-052 on **all $202$** parents ($k_P=2\dim I$,
  $\dim S=2\dim I^\infty$, identical cases, zero mismatches);
* **sharper than §6**: on all $192$ demoting parents $I^2=0$ (nilpotency index
  exactly $2$) and on all $10$ immune ones $I^2=I$ — the catalogue realizes only
  the two extremes, never an intermediate ideal;
* **Corollary J-G1 (odd lattices).** If $\ell,m$ are both odd then $R$ is
  semisimple, every ideal is idempotent, and *no* parent admits a single-row
  syzygy demotion: the collapse mechanism needs an even lattice dimension.
  Machine: $3{,}600$ parents over nine odd lattices, zero violations. Among the
  seven published BB instances exactly one is protected this way — the
  $[[90,8,10]]$ code on $(\ell,m)=(15,3)$ is **immune**, while
  $[[72,12,6]]$, $[[108,8,10]]$, the Gross $[[144,12,12]]$, $[[288,12,18]]$,
  $[[360,12,\le24]]$ and $[[756,16,\le34]]$ all have $I^2=0$, i.e. *every*
  nonzero class demotes.
* **the dual ideal, audited independently.** The module annihilator points the
  other way and carries a bar: since $H_X$ acts through transposed circulants,
  $\operatorname{Ann}_R(M)=(\bar a,\bar b)$, verified on $202/202$ parents by two
  independent computations (the quotient *action* on $M$ versus the ideal
  generated by the two $H_Z$ blocks), together with the Frobenius duality
  $\dim I+\dim\operatorname{Ann}_R(M)=\ell m$. The bar is not cosmetic: $(a,b)$
  is bar-invariant on only $156$ of the $202$ parents. Consequently the two
  annihilators are never interchangeable — substituting
  $\operatorname{Ann}_R(M)$ into the formula returns $2(\ell m-k_P/2)$ — that is
  $100$, $172$ and $344$ across the ten immune parents — instead of the correct
  $\dim S=k_P\in\{8,16\}$, and it misclassifies the mixed witness as
  demote-full. The
  zero/deficient factor split is tracked by $\dim I^\infty$ and
  $\dim I-\dim I^\infty$: immune $\iff$ the nilpotent part vanishes,
  demote-full $\iff$ the zero part does, and the Theorem J-J witness is the only
  instance with both nonzero;

The weight dependence is now also exact. Split $G=\mathbb Z_\ell\times\mathbb
Z_m$ into its $2$-part and odd part, $G=G_2\times G_{\rm odd}$. Then
$\mathbb F_2[G_{\rm odd}]=\prod_\chi F_\chi$ is semisimple and $\mathbb
F_2[G_2]$ is local with nilpotent maximal ideal, so the local factors of $R$ are
$R_\chi=F_\chi[G_2]$, in which $G_2$ is an $F_\chi$-basis.

**Theorem J-H (exact-vanishing coset law).** For $a\in R$ partition its support
by $2$-part, $S_h=\{u\in G_{\rm odd}:hu\in\operatorname{supp}a\}$. Then
$a_\chi=\sum_{h\in G_2}\bigl(\sum_{u\in S_h}\chi(u)\bigr)h$; hence $a_\chi=0$
iff every coset sum vanishes, $a_\chi$ is a non-unit iff the total sum
vanishes, and $a_\chi=0$ forces $|S_h|\ge2$ on every occupied coset. In
particular a polynomial of weight $\le3$ that vanishes exactly at some factor
occupies a **single** coset, and is then a *scalar* — zero or unit — at every
factor.

**Theorem J-I (weight law).** If $\operatorname{wt}(A)\le3$ and
$\operatorname{wt}(B)\le3$ then $(A,B)$ is **not mixed**, on any lattice: mixed
needs a factor where both vanish (forcing single cosets, hence zero-or-unit
everywhere) and a factor where both are non-units but not both zero — which is
then impossible. Equivalently, a mixed parent requires weight $\ge4$ spread
over $\ge2$ cosets in one of $A,B$.

This *explains* the all-or-nothing law of §6 without any enumeration: the whole
published catalogue is weight-$(3,3)$. Two machine confirmations: an exhaustive
census of **653{,}022{,}021** weight-$\le3$ pairs (translation-normalised) over
$18$ lattices — every catalogue lattice, the previously underexplored $(12,12)$
and $(15,12)$, and a parity grid — finds **zero** mixed pairs (EXP-054, $450$
sampled pairs re-decided by an independent route, zero disagreements); and the
single-coset predicate separates the catalogue **exactly**: all $10$ immune
parents are single-coset, all $192$ demoting parents are not (EXP-053). For
$m=6$ this reads off by eye: immunity $\iff$ all $A$-monomials share one
$y$-parity and all $B$-monomials share one $y$-parity.

**Theorem J-J (mixed exists; the bound $4$ is attained).** On $(\ell,m)=(2,3)$,
$A=(1+x)(y+y^2)$ and $B=Ay$ give $\dim I=4$, $\dim I^2=\dim I^\infty=2$: a
**mixed** parent with $k_P=8$, $\dim S=4$, and exactly $2^8-2^4=240$ of $255$
nonzero classes demoting. Verified three ways — ideal chain, the EXP-052 module
chain $[8,4,4]$, and brute-force evaluation of all $255$ classes through the
EXP-047 rank test. Problem J-A ("does a mixed parent exist?") is therefore
**resolved affirmatively**, while Theorem J-I shows it can never happen inside
the trinomial family used in the literature. What remains open is the sharp
weight/lattice map at weight $\ge4$ and whether any mixed parent has useful
$[[n,k,d]]$ (the witness above is a $[[12,8]]$ toy).


### 6.2 A certified solver-free distance ceiling on odd lattices

Corollary J-G1 says the collapse channel is empty when $\ell,m$ are both odd.
That raises the obvious design question — *what actually lives in the
collapse-free region?* — and answering it required one new tool and one
correction to our own bookkeeping.

**What is already known (reproduced, not claimed).** For $\ell,m$ odd, $|G|$ is
odd, so $R=\mathbb F_2[G]$ is semisimple by Maschke and splits as
$\prod_\chi\mathbb F_\chi$ over the Frobenius orbits of characters.
Multiplication by $a$ is then diagonal, so
$\operatorname{Ann}(a)=\bigoplus_{a(\chi)=0}\mathbb F_\chi$ and
$$k=2\dim_{\mathbb F_2}I
=2\sum_{\chi:\,a(\chi)=b(\chi)=0}[\mathbb F_\chi:\mathbb F_2],$$
the degree-weighted odd-lattice common-root law.
Panteleev–Kalachev (arXiv:1904.02703, Prop. 1) as
$k=2\deg\gcd(a,b,x^\ell-1)$ in the cyclic case, Lin–Pryadko
(arXiv:2306.16400, Eq. 47), Wang–Mueller (arXiv:2408.10001v4, Eq. 11) for
coprime lattices via $\pi=xy$, Postema–Kokkelmans (arXiv:2502.17052v4,
Thm. 2.6), and Eberhardt–Steffan (arXiv:2407.03973v1, Cor. 2.11–2.12: "if
$\ell$ and $m$ are odd, all BB codes are principal"). We use it, we do not
claim it.

**Theorem K (reciprocal-pole logical isomorphism; bar-aware form of the
published principal-code structure).** The coefficient ideal
$I=\operatorname{Ann}_{\rm left}(a,b)$ is what `nullspace(HX.T)` computes. It
is **not** itself the physical right kernel. Put $J=\bar I$, where
$\bar{\cdot}$ is $x\mapsto x^{-1},y\mapsto y^{-1}$, and define
$\mathcal P=J\oplus J$. Then
$$\mathcal P\subseteq\ker H_X,\qquad
  \mathcal P\cap S_Z=0,\qquad
  \dim\mathcal P=k,$$
so the quotient map is an isomorphism
$$\boxed{\mathcal P\cong\ker H_X/S_Z.}$$
*Proof.* In the character decomposition $I$ is supported on the common-zero set
$Z$ and $J=\bar I$ on $Z^{-1}$. The physical column action of $H_X$ evaluates
the reciprocal character, so $\mathcal P\subseteq\ker H_X$. Meanwhile
$S_Z=\{(\lambda\bar b,\lambda\bar a)\}$, and at
$\chi\in Z^{-1}$ its generator is
$(b(\chi^{-1}),a(\chi^{-1}))=(0,0)$; hence
$\mathcal P\cap S_Z=0$. Finally $\dim\mathcal P=2\dim I=k$ by the published rate
law. $\square$

This exact logical transversal is the principal-code isomorphism of
Eberhardt–Steffan in the repository's matrix convention. The solver-free
corollaries are
$$d(P)=\min_{0\ne p\in\mathcal P}\min_{s\in S_Z}\operatorname{wt}(p+s),
\qquad
d(P)\le d_{\rm pole}:=\min_{0\ne u\in J}\operatorname{wt}(u).$$
The raw pole ceiling costs $O(2^{k/2}\ell m)$ bit operations — microseconds for
$k\le32$, no solver, decoder, or sampling. Reducing selected pole vectors
modulo $S_Z$ under several information-set orders gives much tighter
**self-certifying** logical witnesses: randomisation changes only tightness,
never soundness. On the published $[[90,8,10]]$, for example, the raw ceiling
is $20$ and the reduced witness has weight exactly $10$.

The bar is load-bearing. On the published $(7,7)\,[[98,6,12]]$ code every basis
vector of raw $I$ fails $H_X(u,0)^T=0$, whereas every basis vector of
$J=\bar I$ passes. An intermediate EXP-055 draft used raw $I$ in a fallback;
the regression caught it before the final PDF and the route is retracted. The
earlier restricted $I\cap\bar I$ statement happened to be sound but was
unnecessarily weak. Full statement, proof, and retraction:
`notes/theorem_k_certified_ceiling.md`.

**Validation against the literature, and a correction.** A read-only novelty
check found the pole isomorphism already published and did not find the
explicit general-BB minimum-pole-weight ceiling or its exhaustive use as a
rejection oracle; the latter is nevertheless a short corollary and is not
oversold as a deep new theorem. The same check refuted a claim we had been
carrying, that $[[90,8,10]]$ on $(15,3)$ is the only odd$\times$odd BB instance
in print. It is not: we sourced $27$ instances from arXiv:2308.07915,
arXiv:2407.03973v1, arXiv:2408.10001v4 and arXiv:2502.17052v4, including
$(9,9)\,[[162,8,12]]$, $(9,15)\,[[270,8,18]]$,
$(7,7)\,[[98,6,12]]$ and $(3,27)\,[[162,8,14]]$. Our earlier statement was
scoped to the seven instances of our own baseline table, where it is correct;
stated of the literature it was wrong, and it is withdrawn.

On those $27$ instances (EXP-055), the reciprocal-pole isomorphism passes
**27/27** independent rank audits. We recomputed $k$ by three independent
routes — annihilator dimension, $k=n-\operatorname{rank}H_X-
\operatorname{rank}H_Z$, and the published $\gcd$ formula on coprime lattices.
All three agree internally on all $27$ and with printed $k$ on $25$. Two
transcribed Wang–Mueller App. C rows, $(5,9)$ and $(7,11)$, give internal $k=0$
against printed $4$ and $6$; because we did not read the appendix ourselves,
we record them as unreproduced and exclude them. All 25 reproduced rows satisfy
`pole ceiling >= reported d`, but this is only a sanity check:
Wang–Mueller's distances are BP-OSD `distance_upperbound` outputs and
Postema–Kokkelmans labels its table values Monte-Carlo estimates. Six rows are
independently exact-certified here (Bravyi $[[90,8,10]]$, three small Postema
rows, Wang–Mueller's $[[126,12,10]]$ exactified by EXP-055, and its
$[[162,8,14]]$ exactified by EXP-056): zero ceiling violations; slack min/median/max $2/20/22$. The screen admits **only independent two-sided
certificates** as domination thresholds; no decoder estimate or un-replayed
source distance can reject a candidate.

**Exhaustive census of the region.** Over all $65$ odd lattices with
$\ell m\le180$ ($n\le360$) we enumerated every weight-$\le3$ pair —
$4.23\times10^{9}$ pairs — computing $k$ by the ideal route and cross-checking
against the matrix route with **zero mismatches**, and re-verifying the J-G1
mechanism directly ($273$ idempotence tests, zero violations, on lattices
outside our catalogue's $m\in\{3,6\}$). Grouping polynomials by
$\operatorname{Ann}$ makes this exact rather than sampled: $k$ and the ceiling
depend only on the pair of annihilators.

**Exact discoveries before fixed-point closure.** The first exact screen wave
produced two connected, row-space-indecomposable BB codes:
$$
\begin{array}{c|c|c|c}
(\ell,m)&A&B&[[n,k,d]]\\\hline
(5,3)&1+x+x^3y&1+xy+x^4y&[[30,8,4]]\\
(9,3)&1+y+x^3y^2&1+x+x^2&[[54,8,6]]\\
(7,9)&1+\pi+\pi^{58}&\pi^3+\pi^{16}+\pi^{44}&[[126,12,10]]
\end{array}
$$
All sector distances are independently CP-SAT exact and carry explicit
witnesses. A primary-table search found no matching $[[30,8,4]]$ BB
constructor (Wang--Mueller's only $n=30$ row is $[[30,4,6]]$), but Grassl's
QECC table contains an explicit general stabilizer $[[30,8,7]]$: the BB
constructor is new within the checked BB corpus but **globally dominated**, not
a code-parameter or end-to-end Pareto improvement. The $[[54,8,6]]$ and
$[[126,12,10]]$ results independently exactify parameters Wang--Mueller
reported through BP-OSD (the latter costs 306s). All three become hash-bound
local references for the fixed-point screen; their certificates are
`results/certificates/exp055_discovered_references.json`.

**Orbit-class exactification at $n=162$ (EXP-056).** For the Wang–Mueller
$(3,27)$ code
$A=1+y^{10}+y^{14}$, $B=y^{12}+x+x^2$, the $255$ nonzero logical classes form
$20$ orbits under the verified order-$162$ translation/reflection group.
Each representative is encoded as an affine coset of the 77-dimensional
$Z$-stabilizer rowspace. Every $H_X$ column has odd degree three, so all kernel
words have even weight and excluding weight $\le12$ also excludes weight
$13$. Kissat proves all 20 class CNFs UNSAT and repeats all 20 on fresh,
digest-bound replay; an explicit weight-14 word passes independent NumPy and
bitset checks. The exact BB duality permutation gives
$d_X=d_Z=14$, promoting the source's BP-OSD estimate to a local two-sided
certificate (`exp056_wm_162_8_14_distance.json`). Generic fixed-functional
sectors do **not** inherit the monolithic origin anchor; a synthetic
counterexample caught and repaired that latent alternative-mode bug (FR-027).

**Fixed-point Pareto screen (complete through $n=162$).** We promoted the three
discoveries and the EXP-056 exactification into the hash-bound reference set
and reran every weight-3 pair on all 13 odd lattices with a nonempty frontier
through $n=162$, $8\le k\le24$. Exact
translation/unit/block-swap/$x\leftrightarrow y$ quotienting leaves **2,132**
classes representing **51,769** translation-normalised pairs. Of these, 1,928
have an admissible independently exact reference: **all 1,928 are dominated** —
1,804 by reduced-pole logical witnesses, 119 by bounded CDCL witnesses, and
five by the exact CP-SAT fallback. All 1,923 persisted witnesses are checked in
$\ker H_X\setminus S_Z$. The other 204 high-$k$ classes have no certified
reference and are labelled `no_reference`, not dominated. There are **zero survivors and zero undecided** among referenced classes. Every shard is bound
to the census SHA-256, certificate-hashed reference set, pure-validator version
and decision protocol; the assembled screen is
`results/processed/exp055_odd_lattice_screen.json`.

Scope is load-bearing: the algebraic $k$ census covers all 65 odd lattices
through $n=360$, but the exact Pareto screen stops at $n=162$. No
$n>162$ distance-closure claim is made.

---

## 7. The exact trade law, its residual, and remaining conjectures

**Theorem I (restated; proved in §3.2).** If $T(P)\ge k_P/2$ and $d_Q>d_Z(P)$,
then $k_Q=k_P-T(P)$ exactly. This covers $138$ of $139$ certified parents and
all $57$ family-closed ones at $n\le 144$; the saturation observed in §5.2 and
the $496/496$ small-lattice increases are instances of it.

Three layers of machine evidence now sit beneath the theorem and delimit its
residual:

1. *Catalogue upper law.* On all $284$ classified rows,
   $\dim\bar\Delta\le T(P)$ with zero exceptions (spread of
   $\dim\bar\Delta-T$ only non-positive, down to $-24$). For $T\ge k_P/2$
   parents this is Lemma 3; the law holding also on the one $T<k_P/2$ parent's
   rows is evidence the residual class is benign.
2. *Small-lattice probe (EXP-040).* $64$ parents across four lattices,
   $26{,}898$ unseen $\delta>0$ perturbations: $\dim\bar\Delta<T$ on $23{,}634$,
   $=T$ on $3{,}264$, $>T$ on **zero**; all $496$ strict increases at $=T$, with
   actual containment $\bar M=\bar\Delta$ verified, not just dimensions. All
   probed parents had $T\in\{k_P/2,k_P\}$, so these are theorem-forced.
3. *Reversals.* All $7$ catalogue reversals at $\dim\bar\Delta=T=k_P/2$ with
   $\bar M=\bar\Delta$.

**Conjecture B′ (residual saturation) — decided on the known witness.** For parents
with $T<k_P/2$ (one known: `9a7638586033`, where increase sandwiches to
$4\le\dim\bar\Delta\le 6$),
distance increase still forces $\dim\bar\Delta=T$. On the only known residual parent
the question is now *decided*, and the answer is stronger than saturation: increase
over `9a7638586033` is **impossible**, because the dressing rows of every valid
perturbation live in a fixed two-dimensional space (EXP-045): over its entire valid
family ($\dim V=112$, $2{,}533{,}006{,}645$ instances enumerated exhaustively through
weight $6$, per-weight counts matching $\binom{112}{k}$ exactly),
$\dim\bar\Delta\in\{0,2\}\subset[0,2]<4=T$, so the Theorem-H premise is never met and
$d_Q\le d_Z(P)=6$ for *every* valid $[C\;D]$ of that parent. The sandwich
$[4,6]$ is never attained. Whether the residual class contains any other parent at
larger lattices ($\ell\cdot m>40$) remains open (§6); on the known witness the
$k$-halving confinement picture is therefore complete with no exceptional direction.

**Conjecture C ($k$-halving confinement).** Any perturbation with $d_Q>d_Z(P)$
has $k_Q/k_P=1/2$. On every $T=k_P/2$ parent this is immediate from
Theorem I; the conjecture is residual exactly on $T<k_P/2$ parents and on the
$63$ uncertified ones ($n\in\{180,360\}$): 26 at $n=180$ and 37 at $n=360$. Every one of the $284$ classified
rows with $k_Q/k_P\in\{1,3/4,5/6\}$ is capped *a priori* (§5.1, zero
exceptions).

**Conjecture D (distance parity).** Certified reversals gain exactly $+2$ in
distance. *Evidence:* the uniform gap $d_Q-d_Z(P)=2$ across all $7$ reversals,
whose exact $d_Q$ values are certified; consistent with $d_Z$ even throughout.
(The probe stores only the increase decision, not the achieved $d_Q$, so its
$496$ increases are not parity evidence.)

The open problem that remains is not *whether* perturbation pays — it pays
exactly $T$ — but *whether anything on the CSS side of the envelope can ever be
beaten this way*: §6's domination table says no for all seven known increases.

---

## 8. Reproduction and provenance

Environment: macOS Darwin $25.5.0$ arm64, Apple M3 Ultra, $28$ logical CPUs, $96$ GiB;
Python $3.13.9$ in `.venv`; `numpy 2.4.6`, `scipy 1.18.0`, `stim 1.16.0`,
`pymatching 2.4.0`, `sinter 1.16.0`, `ldpc 2.4.1`, `galois 0.4.11`, `ortools 9.15.6755`,
`python-sat 1.9.dev13`; solvers CaDiCaL $1.9.5$ and Kissat $4.0.4$ via PySAT. Full suite: $972$ passing tests, $1$ skipped ($973$ collected).

| claim | artifact | experiment |
|---|---|---|
| Theorem G, machine-checked on catalogue rows | `src/qec_research/codes/pbb_survival.py`, `tests/test_pbb_survival.py` ($15$ checks) | EXP-038 |
| Theorem H, module + exact $T$ enumeration | `src/qec_research/codes/pbb_nogo.py`, `tests/test_pbb_nogo.py` ($17$ checks) | EXP-039 |
| Gross closure, $63$ closed parents, $253$ rows capped, complete at $n\le 144$ | `results/partial_runs/exp039_nogo_module.json` | EXP-039 |
| Saturation table, $7/7$ non-vacuous, gate | `exp039_nogo_module.py gate` | EXP-039 |
| Certified $d_Q$ lower bounds for reversals | `results/partial_runs/exp037/row_*.json` | EXP-037 |
| Theorem I legs: $k_P=2\dim L$ (202 parents), $\dim\bar\Delta\le k_P/2$ (368 rows), scope exception | `tests/test_pbb_theorems.py`, `results/processed/exp040_saturation_probe.json` | EXP-039/040 |
| All $7$ reversals CSS-dominated at equal $n$ ($162$ exact candidates) | `results/processed/exp036_envelope_check.json` | EXP-036/037 |
| $5$ certified reversals, $41$ dominations of $87$ | `results/partial_runs/exp036_delta_closure.json` | EXP-036 |
| Exact $d=12$ for `12_6_0193`, two-sided | `results/certificates/pbb_12_6_0193_distance.json` | EXP-035 |
| Genuinely non-CSS (`12_6_0193`) | `notes/novelty_matrix.md`, EXP-034 records | EXP-034 |
| Catalogue-wide CSS-envelope classification | `results/partial_runs/exp037_envelope_classification.json` | EXP-037 |
| Strict-provenance $\delta>0$ audit ($155$ rows) | `results/processed/exp027_delta_audit.json` | EXP-027 |
| Translation symmetry break soundness, $6.2\times$ | `tests/test_translation_symmetry_break.py` | FR-022 |
| DEM mechanism-distance intervals, LER, latency | `reports/technical_report.md` §5–6 | EXP-016/029 |
| Theorem J-C: $D{=}0$ collapse family decidable, $192/202$ parents | `results/processed/exp046_generic_syzygy.json` | EXP-046 |
| X-U distance certificates $d_P=10$ + family-closure ($T=8=k_P$ respectively $T=16=k_P$), `phase2_109` $d_P=6$, $T=4$ | `results/partial_runs/exp039/parent_e15b3d6480981d1e.json`, `results/partial_runs/exp039/parent_765498fc67b2b4b9.json`, `results/partial_runs/exp039/parent_dc98a5e4c350a656.json` | EXP-039/049 |
| X-monotone channel exclusion on the $10$ immune parents: $8$ INFEASIBLE certificates + $2$ weight-$8$ witnesses | `results/processed/exp051_xmonotone_lightscan.json` + `tests/test_exp051_xmonotone_lightscan.py` ($3$ checks) | EXP-051 |
| Final X-side partition $192/8/2$ ($=202$) | jointly derived: EXP-047 demote/immune decision $+$ EXP-049 dedicated certificates $+$ EXP-051 exclusions | EXP-047/049/051 |
| Theorem J-E immunity decision, $10$ immune exact | `results/processed/exp047_exact_demotion_decision.json` | EXP-047 |
| Exact collapse census: $133/133$ terminal-UNSAT, $99$ strict drops ($98$ full-$k$) | `results/processed/exp048_exact_collapse_sweep.json` | EXP-048 |
| Demote-fraction enumeration, $196$ parents ($k_P\le 20$, 2{,}693{,}100 classes) | `results/processed/exp050_demote_fractions.json` + `results/partial_runs/exp050/fraction_*.json` | EXP-050 |
| Theorem J-E′′: demote fixed set $= I^\infty M$; $202/202$ fractions $\in\{0,1\}$ ($192/10$, $0$ mixed), chains $\le 3$, $196$ integer crosschecks | `results/processed/exp052_ideal_power_trichotomy.json` + `tests/test_exp052_ideal_power_trichotomy.py` ($11$ checks) | EXP-052 |
| Theorems J-G/J-H/J-I/J-J: ideal-invariance ($\dim S=2\dim I^\infty$, $202/202$ vs EXP-052), $I^2=0$ on all $192$ demoting / $I^2=I$ on all $10$ immune, odd-lattice corollary ($3{,}600$ parents), baseline table ($[[90,8,10]]$ immune), coset criterion exact ($10$ vs $192$), lemma battery, mixed witness (3 routes) | `results/processed/exp053_ideal_classification.json` + `tests/test_exp053_ideal_invariant.py` ($13$ checks) | EXP-053 |
| Weight-$\le3$ mixed census: $653{,}022{,}021$ pairs over $18$ lattices, zero mixed, $450$ independent cross-checks | `results/processed/exp054_mixed_census.json` | EXP-054 |
| Odd-lattice algebraic census: $65$ lattices, $4{,}229{,}823{,}962$ pairs, zero $k$ mismatches, $273$ idempotence checks | `results/processed/exp055_odd_lattice_sweep.json` + config-bound `results/partial_runs/exp055/*.json` | EXP-055 |
| Reciprocal-pole isomorphism $27/27$; six locally exact ceiling checks; source-estimate audit | `results/processed/exp055_literature_validation.json` + `notes/theorem_k_certified_ceiling.md` | EXP-055/056 |
| Exact $[[162,8,14]]$: 20/20 logical-class orbits UNSAT and replayed, weight-14 witness, $d_X=d_Z$ duality | `results/certificates/exp056_wm_162_8_14_distance.json` + `tests/test_exp056_odd_distance.py` | EXP-056 |
| Fixed-point screen through $n=162$: 2,132 classes / 51,769 pairs, $1,928/1,928$ referenced dominated, 204 no-reference, zero survivor/undecided; 1,923 explicit witnesses | `results/processed/exp055_odd_lattice_screen.json` + `results/certificates/exp055_odd_lattice_survivors.json` + `tests/test_exp055_odd_lattice.py` | EXP-055/056 |

Every SAT decision records its canonical CNF SHA-256 and encoding version; verdicts are
re-derived from rebuilt matrices on replay, and stamps that fail to hash-bind are
rejected rather than trusted. Forgery regression tests confirm that crude tampering dies
arithmetically or under replay.

**Reproduce the headline result:**

```bash
cd math/qec
PYTHONPATH=src .venv/bin/python experiments/exp039_nogo_module.py run --ns 144
PYTHONPATH=src .venv/bin/python experiments/exp039_nogo_module.py gate
PYTHONPATH=src .venv/bin/python experiments/exp055_odd_lattice_sweep.py run \
  --max-dim 180 --workers 24 --force
PYTHONPATH=src .venv/bin/python experiments/exp056_odd_distance.py \
  run-classes --indexes all --force
PYTHONPATH=src .venv/bin/python experiments/exp056_odd_distance.py \
  run-classes --indexes all --replay
PYTHONPATH=src .venv/bin/python experiments/exp056_odd_distance.py assemble
PYTHONPATH=src .venv/bin/python experiments/exp055_odd_lattice_sweep.py \
  literature
PYTHONPATH=src .venv/bin/python experiments/exp055_odd_lattice_sweep.py screen \
  --lattices 3x3,5x3,7x3,9x3,9x5,15x3,7x7,21x3,9x7,15x5,25x3,9x9,27x3 \
  --workers 8 --time-limit 120 --force
PYTHONPATH=src .venv/bin/python experiments/exp055_odd_lattice_sweep.py \
  screen-certify
PYTHONPATH=src .venv/bin/python -m pytest \
  tests/test_exp055_odd_lattice.py tests/test_exp056_odd_distance.py -q
```

---

## 9. Related work and novelty

Bivariate-bicycle codes and the Gross code are due to Bravyi et al.
(arXiv:2308.07915). The PBB construction and the $368$-row catalogue are from
arXiv:2606.02418. Cruz-Benito et al. (arXiv:2606.02418 v1, Sec. III.2, Lemma 0) give the
validity condition; we derive it independently from the binary symplectic Gram matrix and
confirm the correct condition is *symmetry* of $AC^{\mathsf T}+BD^{\mathsf T}$, not its
vanishing. Multivariate bicycle generalisations are in arXiv:2406.19151; existence and
characterisation results for BB codes in arXiv:2502.17052. ASC (arXiv:2603.21499) certifies
no-depth-$6$ for IBM BB codes; our depth criterion is independent and basis-independent.

To our knowledge the following are new here: (a) the identification of $\Delta$ as an
$R$-submodule and the resulting orbit-absorption argument; (b) the exact dimension
identity $k_Q=k_P-\dim\bar\Delta$ as a certified structural fact rather than an observed
coincidence; (c) the perturbation-independent no-go bound (Theorem H) and its exact
$T$ enumeration with a terminating UNSAT certificate; (d) the closure of the entire
perturbation family over the Gross code; (e) the saturation observation $k_Q=k_P-T$ on
all certified reversals; and (f) the proved confinement of distance increase to a
$k$-halving class.

---

## 10. Conclusion

The PBB construction's extra freedom is real but not free. It is governed by a single
submodule $\Delta$ whose dimension is exactly the dimension deficit, and whose module
structure forces any distance increase to absorb entire translation orbits of the
parent's minimum-weight logicals. This yields a computable, perturbation-independent
obstruction that closes $63$ of $202$ published parents outright — the Gross code among
them — caps $253$ of $368$ catalogue rows with no search, confines distance increase to a
$k$-halving class, and is saturated with slack exactly zero on all $7$ certified
distance-increasing perturbations known to us.

The abstract $[[n,k,d]]$ route past the CSS bivariate-bicycle envelope is, over these
parents, closed by theorem. What remains genuinely open is the circuit-level question:
whether a verified fault-tolerant syndrome-extraction circuit for a non-CSS PBB can beat
the Gross code end-to-end under identical noise. Nothing in our circuit-level evidence so
far favours it.
