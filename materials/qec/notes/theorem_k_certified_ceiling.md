# Theorem K — the reciprocal-pole logical isomorphism and certified BB distance bounds

**Status (2026-08-22):** proved in the repository's matrix convention and
machine-verified on 27 sourced odd-lattice BB instances. EXP-056 adds a
replayed exact-distance closure for the published $(3,27)$ $[[162,8,14]]$
instance. The underlying
principal-code logical-space isomorphism is **published** (Eberhardt–Steffan,
arXiv:2407.03973v1, Corollaries 2.11–2.12); do not claim it as new. A novelty
search did not find the explicit general-BB *minimum-pole-weight ceiling* or its
use as a solver-free exhaustive rejection oracle, but both are short corollaries
of the published isomorphism and should be presented as such, not oversold as a
deep new theorem.

## 0. Setting and the load-bearing bar

Let $\ell,m$ be odd, $G=\mathbb Z_\ell\times\mathbb Z_m$, and
$$R=\mathbb F_2[G]
  =\mathbb F_2[x,y]/(x^\ell-1,y^m-1).$$
For $a,b\in R$, the CSS bivariate-bicycle code is
$$H_X=[A\ B],\qquad H_Z=[B^{\mathsf T}\ A^{\mathsf T}],\qquad n=2\ell m,$$
where $A,B$ are the binary circulants of $a,b$. Put
$S_Z=\operatorname{rowspace}H_Z$ and let $\bar{\cdot}$ be the involution
$x\mapsto x^{-1},y\mapsto y^{-1}$.

The code stores polynomial coefficient vectors as **rows**. Consequently
```python
I = nullspace(HX.T)
```
is a *left*-annihilator coefficient ideal. A physical column vector in the
right kernel of $H_X$ is obtained only after applying the reciprocal map. This
is the convention that must never be elided:
$$I=\operatorname{Ann}_{\rm left}(a,b),\qquad J=\bar I
  =\operatorname{Ann}_{\rm right}(a,b).$$
The distinction is invisible when $I=\bar I$ (for example the published
$[[90,8,10]]$ code) and load-bearing otherwise. On the published
$(7,7)\ [[98,6,12]]$ code every basis vector of raw $I$ fails the equation
$H_X(u,0)^{\mathsf T}=0$, while every basis vector of $J=\bar I$ satisfies it.
`test_bar_convention_on_noninvariant_code` locks this regression.

## 1. Published rate and logical-space input

Since $|G|$ is odd, Maschke's theorem makes $R$ semisimple:
$$R\cong\prod_{\chi\in X}\mathbb F_\chi,$$
with $X$ the Frobenius orbits of characters. Multiplication is diagonal, so if
$Z=\{\chi:a(\chi)=b(\chi)=0\}$, then
$$I\cong\bigoplus_{\chi\in Z}\mathbb F_\chi,\qquad
  k=2\dim_{\mathbb F_2}I
   =2\sum_{\chi\in Z}[\mathbb F_\chi:\mathbb F_2].$$
This is published in equivalent forms by:

* Panteleev–Kalachev, arXiv:1904.02703, Proposition 1
  ($k=2\deg\gcd(a,b,x^\ell-1)$ in the cyclic case);
* Lin–Pryadko, arXiv:2306.16400, Eq. (47);
* Wang–Mueller, arXiv:2408.10001v4, Eq. (11), for coprime lattices via
  $\pi=xy$;
* Postema–Kokkelmans, arXiv:2502.17052v4, Theorem 2.6;
* Eberhardt–Steffan, arXiv:2407.03973v1, Corollaries 2.11–2.12
  (if $\ell,m$ are odd, all BB codes are principal).

EXP-055 reproduces $k$ independently by annihilator dimension, direct matrix
ranks, and the published gcd formula where $\gcd(\ell,m)=1$.

## 2. The exact reciprocal-pole isomorphism

**Theorem K (matrix-convention form of the published principal-code
isomorphism).** Define
$$\mathcal P=J\oplus J\subseteq\mathbb F_2^{2\ell m},
  \qquad J=\bar I.$$
Then
$$\mathcal P\subseteq\ker H_X,\qquad
  \mathcal P\cap S_Z=0,\qquad
  \dim\mathcal P=k.$$
Therefore the quotient map restricts to an isomorphism
$$\boxed{\mathcal P\;\cong\;\ker H_X/S_Z.}$$
Every logical class has a unique reciprocal-pole representative $(u,v)$ with
$u,v\in J$.

**Proof.** In the character decomposition, $I$ is supported on $Z$ and
$J=\bar I$ on $Z^{-1}$. Physical column multiplication by $H_X$ evaluates the
reciprocal character, so both blocks of every $(u,v)\in J\oplus J$ are killed:
$\mathcal P\subseteq\ker H_X$. Meanwhile
$$S_Z=\{(\lambda\bar b,\lambda\bar a):\lambda\in R\}.$$
At a character $\chi\in Z^{-1}$,
$(\bar b(\chi),\bar a(\chi))=(b(\chi^{-1}),a(\chi^{-1}))=(0,0)$; hence $S_Z$
has zero component on every block that supports $\mathcal P$, and
$\mathcal P\cap S_Z=0$. Finally
$\dim\mathcal P=2\dim J=2\dim I=k=\dim(\ker H_X/S_Z)$ by the published rate
law. The injective quotient map is therefore surjective. $\square$

This proof also states precisely why using raw $I$ was wrong: $I$ is supported
on $Z$, while the physical pole and the zero blocks of $S_Z$ are supported on
$Z^{-1}$.

## 3. Exact variational formula and solver-free bounds

**Corollary K1 (exact pole-coset formula).**
$$\boxed{
 d(P)=\min_{0\ne p\in\mathcal P}\ \min_{s\in S_Z}\operatorname{wt}(p+s).
}$$
This is just the definition of distance after replacing arbitrary logical
representatives by the exact pole transversal $\mathcal P$.

**Corollary K2 (pole ceiling).** Taking $s=0$ and one pole block zero gives
$$\boxed{
 d(P)\le d_{\rm pole}:=
 \min\{\operatorname{wt}(u):0\ne u\in J\}.
}$$
Bar preserves Hamming weight, so the number equals
$\min\operatorname{wt}(I\setminus\{0\})$, but the *physical witness* is in
$J=\bar I$. For $k\le32$, exhaustive enumeration costs
$O(2^{k/2}\ell m)$ bit operations — microseconds, no solver, decoder, or
sampling.

**Corollary K3 (certified coset reduction).** Pick any nonzero pole
$p\in\mathcal P$ and any stabilizer $s\in S_Z$. Since
$\mathcal P\cap S_Z=0$, $p+s$ is automatically a nontrivial logical. Row
reducing selected pole vectors modulo $S_Z$, under several information-set
orders, therefore returns explicit self-certifying logical witnesses. Every
returned weight is a rigorous upper bound; randomisation changes only tightness,
never soundness.

This second route is much tighter in practice. On the published
$[[90,8,10]]$, the raw pole ceiling is $20$ while the reduced witness has weight
exactly $10$. In the Pareto screen, it turns many 30–90 second CP-SAT calls into
$\approx0.01$ s witness checks.

## 4. Validation and limitations

`experiments/exp055_odd_lattice_sweep.py literature` audits 29 sourced
odd-lattice instances:

* the reciprocal-pole isomorphism holds **29/29**: $\mathcal P\subseteq\ker H_X$,
  $\mathcal P\cap S_Z=0$, and $\dim\mathcal P=k$ by independent rank tests;
* 27 rows reproduce their printed $k$; two transcribed Wang–Mueller App. C rows
  ($(5,9)$ and $(7,11)$) give internal $k=0$ by both our matrix and gcd routes
  against printed 4 and 6, so they are excluded rather than reconciled;
* all 27 reproduced rows satisfy `pole ceiling >= reported d`, but this is only
  a **sanity check** for decoder-estimated rows: Wang–Mueller uses BP-OSD
  `distance_upperbound`, and Postema–Kokkelmans labels its table distances
  Monte-Carlo estimates;
* eight rows are independently exact-certified here: the previous six plus
  Liang et al.'s two $[[234,8,18]]$ rows exactified by EXP-067. There are zero
  pole-ceiling violations, with slack min/median/max $2/22/34$.

The pole ceiling is sound but loose because setting $s=0$ ignores cancellations
inside a logical coset. The reduced-witness route repairs much of that looseness
but remains an upper-bound search, not an exact distance algorithm. Exactness
still requires exhaustive coset search, CP-SAT/SAT, matching bounds, or a proof.
No BP-OSD estimate is admitted as a domination threshold in EXP-055.

The semisimple character proof requires both $\ell,m$ odd. On even lattices the
ring has nilpotents and this pole decomposition is not valid without a local
Artinian refinement.

## 5. Orbit-generated exact distance and the $[[162,8,14]]$ closure

The pole transversal also reduces exact lower-bound work. Let
$L_X=I\oplus I\cong\ker H_Z/S_X$. If vectors $u_1,\ldots,u_r\in L_X$ have
translation orbits spanning $L_X$, then for every weight cap $c$,
$$
\exists z\in\ker H_X\setminus S_Z,\ \operatorname{wt}(z)\le c
\iff
\exists i,\ z'\in\ker H_X:\ \langle z',u_i\rangle=1,\
\operatorname{wt}(z')\le c.
$$
Indeed a nontrivial $z$ pairs with some translated $g u_i$; translating $z$
by $g^{-1}$ preserves its class status and weight and moves the pairing to
$u_i$. Since every ideal in the odd semisimple group algebra is principal,
one pole generator per physical block suffices. EXP-056 machine-checks the
two orbit ranks, their union quotient rank $k$, and the opposite-quotient
pairing rank $k$ on the target.

A global translation-origin clause is **not sound after fixing one arbitrary
functional sector**: translation can move that functional to another sector.
`sector_instance` therefore drops it. EXP-056 uses only the subgroup fixing
the selected functional modulo stabilizers and inserts one support coordinate
from each subgroup orbit. The hypotheses—code invariance, functional
invariance, and coordinate-orbit coverage—are all replayed before UNSAT may
count. A two-variable synthetic regression demonstrates that retaining the
global clause can change SAT to UNSAT.
The two-sector route is retained as a diagnostic speed probe only: it has no
digest-bound replay writer and `canonical_route_enabled` is false. The exact
certificate below rests exclusively on the fully replayed fixed-class route.

For the Wang–Mueller code
$$
(\ell,m)=(3,27),\quad
a=1+y^{10}+y^{14},\quad b=y^{12}+x+x^2,
$$
EXP-056 uses the stronger fixed-class form. The $2^8-1=255$ nonzero logical
classes form **20** orbits under the verified order-$162$ automorphism group
(translations and $x$-reflection), with sizes $6,9,18$. Each representative
is the affine coset $p+S_Z$ and is encoded by 77 independent sparse
$Z$-stabilizer generators. Every $H_X$ column has odd degree three, so summing
all $H_X$ equations proves that every kernel word has even weight; excluding
weight at most 12 therefore excludes weight at most 13.

Kissat returns UNSAT on all 20 hash-bound class CNFs and repeats all 20 UNSATs
on fresh replay. Initial per-class wall times are 39.79–89.19 s (median
62.00 s; 1,312.75 s serial total); replay total is 1,287.67 s. An explicit
weight-14 reciprocal-pole coset representative passes independent NumPy and
bitset checks. Finally the exact coordinate permutation
$QH_XP=H_Z,\ QH_ZP=H_X$ proves $d_X=d_Z$. Hence
$$\boxed{[[162,8,14]]\ \text{with}\ d_X=d_Z=14\ \text{exactly}.}$$
This promotes Wang–Mueller's BP-OSD `distance_upperbound` value to a local
two-sided certificate; the source value remains correctly labelled an
estimate.

## 6. Adaptive fixed-point extension through $n=210$

EXP-057--065 retain the same bar-aware physical convention but add three
independent mechanisms. First, all-odd check-column degree forces even kernel
weight and exact BB duality reduces a distance proof to one side. An adaptive
raw-kissat climb exactifies $[[170,16,10]]$, seven $[[186,10,14]]$ classes,
thirteen $[[210,18,8]]$ classes, $[[210,24,4]]$, $[[210,14,12]]$, and two
$[[210,10,16]]$ classes. A separate $k=8$ row has independent weight-16
witnesses and is dominated, but has no certified lower bound and is not
claimed exact.

Second, for $\gcd(\ell,m)=1$ the map $\pi\mapsto(x,y)$ gives an explicit
coordinate permutation between every factorization of $N=\ell m$. EXP-059
checks every one of the 105 monomials, so the $(15,7)$, $(21,5)$ and $(35,3)$
screens are one exact $N=105$ problem. Third, EXP-063 proves reference
thresholds are monotone: old domination proofs survive a reference extension,
while carried witnesses are rechecked and prior shards are archived.

The final fixed-point screen is complete through $n=210$: 20
nonempty-frontier lattices, 4,020 classes / 122,833 represented pairs, all
3,816 referenced classes dominated, 204 no-reference, zero
survivors/undecided. Of the dominations, 3,364 use reduced-pole witnesses, 395
use bounded/adaptive CDCL witnesses, and 57 use exact CP-SAT fallback records.
All 3,759 explicit witnesses are rechecked in $\ker H_X\setminus S_Z$.

Two new alternatives fail closed rather than become claims: exact affine
trellis widths 48--99 exceed the state gate, and the exact $H=7$ quotient
projection gives physical lower bound only 4.

## 7. Full-automorphism bundling at $n=234$ (EXP-066 bounded phase)

For $G=\mathbb Z_{39}\times\mathbb Z_3\cong
\mathbb Z_{13}\times\mathbb F_3^2$,
$$
\operatorname{Aut}(G)\cong
\mathbb F_{13}^{\times}\times\operatorname{GL}(2,3),
\qquad |\operatorname{Aut}(G)|=12\cdot48=576.
$$
EXP-066 enumerates all 576 coordinate maps, verifies that each is a bijective
homomorphism, and checks both CSS rowspaces before transporting a witness.
The 182 $(39,3)$ classes not closed by the initial pole reduction partition
into 30 exact code-equivalence bundles. A source witness closes 28 bundles /
158 classes; every transported word is then independently checked in the
target's rebuilt $\ker H_X\setminus S_Z$.

The two remaining bundles contain 24 classes. Each has a verified weight-18
logical, but no complete cap-16 exclusion. Therefore EXP-066 proves only
$$
d\le18
$$
for these classes and records them as undecided. It does **not** prove
$[[234,8,18]]$. The expanded comparison surface has 22 lattices, 4,862
classes / 150,581 pairs, 4,634/4,658 referenced classes dominated, 204
no-reference, zero survivors, and 24 undecided. Exact fixed-point closure
remains $n\le210$.

## 8. Exact $n=234$ closure by rooted connected clusters (EXP-067)

The novelty audit missed Liang--Liu--Song--Chen
(arXiv:2503.03827v3, Table III), which already prints two
$[[234,8,18]]$ generalized-toric/BB rows and states that table distances
through 20 were computed exactly by integer programming. EXP-067 enumerates
all quotient-lattice generator images and finds one map for each row into
$\mathbb Z_{39}\times\mathbb Z_3$. The mapped CSS matrices have exactly the
two open EXP-066 representatives: both $H_X$ and $H_Z$ rowspaces agree after
an explicit qubit permutation. Thus $[[234,8,18]]$ is **published prior art**;
the result here is an independent machine certificate and fixed-point closure,
not a new code.

The lower certificate uses the connected-cluster algorithm of
Webster--Jacob--Higgott (arXiv:2603.22532) through the pinned
QEC-pages/dist-m4ri implementation. A minimum nontrivial kernel word has no
proper zero-syndrome subset: a logical subset is already a lighter logical,
while a stabilizer subset can be removed to leave a lighter representative of
the same class. Therefore the first-unsatisfied-check recursion reaches every
minimum logical from its least support coordinate. Translation sends any
supported qubit in the first 117-qubit block to coordinate 0; an exact
block-column-swapped presentation covers words supported only in the second
block. Two rooted searches therefore replace 234 unrooted searches.

For each open representative, EXP-067 exhausts weights 1--16 in the original
and block-swapped presentations and repeats both runs. The accepted engine is
the race-free legacy `STANDALONE` entrypoint; the multithreaded coordinator is
rejected by FR-033. Every run terminates at `-16`, all check columns have odd
degree and hence every kernel word has even weight, and independent NumPy and
bitset paths verify a weight-18 logical plus the exact BB $X/Z$ isometry.
Consequently both representatives have
$$
d_X=d_Z=d=18.
$$

Promoting either exact $[[234,8,18]]$ reference raises the $k=8$ threshold
from 16 to 18. The 24 former undecided classes then carry rechecked
weight-18 witnesses and become dominated. The exact fixed-point surface is
closed through $n=234$: 22 lattices, 4,862 classes / 150,581 represented
pairs, all 4,658/4,658 referenced classes dominated, 204 high-$k$
`no_reference`, and zero survivors or undecided.

EXP-068 closes the global proof-surface gap exposed by final review. The 57
legacy CP-SAT fallback rows now bind explicit physical $Z$-logicals (29 solved
independently, 28 obtained by exact CRT transport). Validator v11 checks every
comparison-class identity or transport cover, current threshold/source and
physical witness on all 22 shards; monotone rebinding validates its input,
output and archive before the aggregate can move.

The trusted-computing boundary is explicit: EXP-067 hardcodes the audited
dist-m4ri binary, source, license and M4RI source digests rather than trusting
the mutable build manifest. The native binary does not emit a proof trace.
Initial/replay agreement, independent algebraic checks, and the published
exact result are redundant evidence, not formal verification of the native
executable.

---

## 9. Retraction record



An intermediate EXP-055 draft used raw $I$ as a physical kernel and proposed a
member fallback when $I\cap\bar I=0$. That fallback was **wrong**; the new
regression on $(7,7)$ caught it before a final PDF or ledger entry was shipped.
The earlier restricted statement using $I\cap\bar I$ happened to be sound
because that subspace is bar-invariant, but it was unnecessarily weak. The
surviving result is the exact bar-aware isomorphism with $J=\bar I$ above.

Machine sources:

* `results/processed/exp055_literature_validation.json`
* `results/processed/exp055_odd_lattice_sweep.json`
* `results/processed/exp055_odd_lattice_screen.json`
* `results/certificates/exp055_odd_lattice_survivors.json`
* `results/certificates/exp056_wm_162_8_14_distance.json`
* `experiments/exp056_odd_distance.py`
* `results/certificates/exp057_170_16_10_distance.json`
* `results/certificates/exp058_186_10_14_distance.json`
* `results/certificates/exp060_210_18_8_distance.json`
* `results/certificates/exp064_210_*_distance.json`
* `results/processed/exp059_coprime_transport.json`
* `results/processed/exp060_n210_ratchet.json`
* `results/processed/exp061_affine_trellis_profile.json`
* `results/processed/exp062_n105_pole_stabilizers.json`
* `results/processed/exp065_n105_k8_projection.json`
* `experiments/exp057_odd_frontier.py` through `exp065_quotient_projection.py`
* `results/processed/exp066_n234_frontier.json`
* `experiments/exp066_n234_frontier.py`
* `tests/test_exp066_n234_frontier.py`
* `results/certificates/exp067_234_8_18_*_distance.json`
* `results/partial_runs/exp067_n234_cluster/`
* `experiments/exp067_n234_connected_cluster.py`
* `experiments/exp068_screen_proof_repair.py`
* `results/partial_runs/exp068_screen_witnesses/`
* `tests/test_exp067_n234_connected_cluster.py`
* `tests/test_exp056_odd_distance.py`
* `tests/test_exp055_odd_lattice.py`
