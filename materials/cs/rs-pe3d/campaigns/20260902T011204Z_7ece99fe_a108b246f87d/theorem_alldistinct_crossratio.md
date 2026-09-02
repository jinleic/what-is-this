# The all-distinct size-four channel is cross-ratio preservation
## gate H-ALLDISTINCT-CROSSRATIO, run `20260902T011204Z_7ece99fe_a108b246f87d`

Agent `RsPe3dH2`, 2026-09-01/02 UTC. Prereg:
`prereg/H_ALLDISTINCT_CROSSRATIO_PREREG_2026-09-01.md`, path-scoped commit
`44a491993b229ecf1789a9f7d64b4b8130e2dcb9`, sha256
`b2b7bfc454abc07047fce8b91a8d2d5dda359464a836a16851c1067e141cf635`,
byte-identical run copy `pre_statement.md` (cmp-verified). No
parameter-dependent compute preceded that commit or `campaign.py init`; no
frozen campaign directory was modified. The characterization supplied by Main
was treated as a candidate and is derived independently in §§1–5. Exact
machine record: `controls_results.json` (176 aggregate asserts, 0 failures;
CPU 0.466/600 s, wall 0.491/900 s; `nice` 15; four thread caps 1; exact
integer GF($p$) arithmetic; bytecode disabled).

## 0. Setting and theorem

Let $A$ have columns $a_u=(1,x_u)^{\mathsf T}$ and $B$ have columns
$b_v=(1,y_v)^{\mathsf T}$, with the $x_u$ pairwise distinct and the $y_v$
pairwise distinct. Then $d_A=d_B=3$. The row-major tensor convention gives

$$h(u,v)=a_u\otimes b_v=(1,y_v,x_u,x_u y_v)^{\mathsf T}\in F^4.$$

An all-distinct size-four support is
$S=\{(u_k,v_k):1\le k\le4\}$ with all $u_k$ distinct and all $v_k$ distinct.
Write $x_k=x_{u_k}$ and $y_k=y_{v_k}$. For four distinct elements define

$$\operatorname{CR}(a,b,c,d)=
\frac{(a-c)(b-d)}{(a-d)(b-c)}.$$

**Theorem AD-CR.** For an all-distinct size-four support the following are
equivalent:

1. $S$ is a circuit of $A\otimes B$;
2. $\det [h(u_1,v_1)\ h(u_2,v_2)\ h(u_3,v_3)\ h(u_4,v_4)]=0$;
3. a nonzero bilinear form $c_0+c_1y+c_2x+c_3xy$ vanishes on all four paired
   points;
4. the four pairs lie on the graph of one nondegenerate Möbius map
   $y=-(c_2x+c_0)/(c_3x+c_1)$, with no selected pole;
5. $\operatorname{CR}(x_1,x_2,x_3,x_4)=
   \operatorname{CR}(y_1,y_2,y_3,y_4)$ under the simultaneous pairing.

Consequently, after ordering each four-subset $U=\{u_1<\cdots<u_4\}$,

$$N_{\rm all}=
\#\left\{(U,\sigma): |U|=4,\ \sigma:U\hookrightarrow[n_B],\
\operatorname{CR}(x_{u_1},\ldots,x_{u_4})=
\operatorname{CR}(y_{\sigma(u_1)},\ldots,y_{\sigma(u_4)})\right\}.$$

This is the exact source of the field dependence left open by the crossing
run: the count changes exactly when the finite-field cross-ratio coincidence
set changes.

## 1. P1 — tensor coordinates and why determinant zero is a circuit

The tensor-coordinate identity is immediate but load-bearing:

$$
(1,x_u)^{\mathsf T}\otimes(1,y_v)^{\mathsf T}
=(1\cdot1,1\cdot y_v,x_u\cdot1,x_u y_v)^{\mathsf T}
=(1,y_v,x_u,x_u y_v)^{\mathsf T}.
$$

It remains to show that a dependent all-distinct four-set is minimal. More
generally, any three columns $a_i\otimes b_i$ with pairwise nonparallel
$a_i\in F^2$ and pairwise nonparallel $b_i\in F^2$ are independent. Suppose

$$\lambda_1a_1\otimes b_1+\lambda_2a_2\otimes b_2+
\lambda_3a_3\otimes b_3=0.$$

Choose a nonzero functional $\phi$ on $F^2$ with $\phi(a_1)=0$. Its kernel is
the line $\langle a_1\rangle$. Since $a_2,a_3$ are not parallel to $a_1$,
$\phi(a_2),\phi(a_3)\ne0$. Applying $\phi\otimes\mathrm{id}$ gives

$$\lambda_2\phi(a_2)b_2+\lambda_3\phi(a_3)b_3=0.$$

$b_2,b_3$ are independent, so $\lambda_2=\lambda_3=0$; then
$\lambda_1a_1\otimes b_1=0$, and both factors are nonzero, so
$\lambda_1=0$. In this GRS setting distinct $x$'s make the $a_i$ pairwise
nonparallel, and distinct $y$'s make the $b_i$ pairwise nonparallel.
Therefore every three-subset of an all-distinct four-support is independent.

The four columns live in $F^4$. Their $4\times4$ determinant vanishes iff their
rank is at most three. Any three already have rank three, so determinant zero
forces rank exactly three and every proper subset independent: a genuine
circuit. Conversely a circuit is dependent, hence its determinant vanishes.
This proves 1 $\Longleftrightarrow$ 2 and prevents the common overclaim
“dependent” $\Rightarrow$ “circuit” without a minimality check.

## 2. P2 — determinant zero iff a bilinear form vanishes

Let $M$ be the $4\times4$ matrix whose $k$th column is
$(1,y_k,x_k,x_ky_k)^{\mathsf T}$. Since $M$ is square,
$\det M=0$ iff its left kernel is nonzero. A left-kernel vector
$c=(c_0,c_1,c_2,c_3)$ satisfies, for every column,

$$cM_{\cdot k}=c_0+c_1y_k+c_2x_k+c_3x_ky_k=0.$$

Thus a nonzero left-kernel vector is exactly a nonzero bilinear form vanishing
on all four points; conversely such a form is a nonzero left-kernel vector.
This proves 2 $\Longleftrightarrow$ 3 without identifying the left-kernel
vector with the (different) circuit relation among columns.

## 3. P3 — the bilinear curve is exactly a nondegenerate Möbius graph

Write

$$F(x,y)=c_0+c_1y+c_2x+c_3xy=D(x)y+N(x),\qquad
D(x)=c_3x+c_1,\quad N(x)=c_2x+c_0.$$

Assume $F\ne0$ vanishes on four pairs with distinct $x_k$ and distinct $y_k$.
There are three possible degeneracies to exclude.

### 3.1 $D$ is not the zero polynomial

If $D\equiv0$, then $c_1=c_3=0$ and $F=N(x)$. A nonzero linear polynomial
cannot vanish at four distinct $x_k$; if it does, $N\equiv0$, contradicting
$F\ne0$. Hence $D\not\equiv0$.

### 3.2 The Möbius determinant cannot vanish

Put

$$\Delta=c_0c_3-c_1c_2.
$$

The degenerate condition named in the prereg is
$c_1c_2=c_0c_3$, i.e. $\Delta=0$. It says the coefficient pairs of $N$ and
$D$ are proportional: because $D\ne0$, $N=\lambda D$ for some $\lambda$.
At every selected $x_k$ with $D(x_k)\ne0$, the equation becomes
$y_k=-\lambda$. A nonzero linear $D$ has at most one root, so at least three
of the four distinct $x_k$ have $D(x_k)\ne0$ and would share the same $y$,
contradicting pairwise-distinct $y_k$. Therefore $\Delta\ne0$.

Equivalently, when $c_3\ne0$ the degenerate curve factors with zero right-hand
side:
$(x+c_1/c_3)(y+c_2/c_3)=0$, a union of one vertical and one horizontal line.
With distinct $x$'s and distinct $y$'s those two lines cannot contain four
selected pairs. This is the geometric form of the same exclusion.

### 3.3 No selected point is a pole

If $D(x_k)=0$, $F(x_k,y_k)=0$ also forces $N(x_k)=0$. But two linear
polynomials $D,N$ have a common root exactly when their coefficient
determinant is zero; that determinant is $c_1c_2-c_0c_3=-\Delta$. Since
$\Delta\ne0$, no selected $x_k$ is a root of $D$. Thus

$$y_k=-\frac{N(x_k)}{D(x_k)}
=-\frac{c_2x_k+c_0}{c_3x_k+c_1}.$$

The matrix

$$\begin{pmatrix}-c_2&-c_0\\c_3&c_1\end{pmatrix}$$

has determinant $\Delta\ne0$, so this is a nondegenerate Möbius map.

The two affine forms requested in the prereg now follow exactly:

- If $c_3=0$, $\Delta=-c_1c_2\ne0$, hence $c_1,c_2\ne0$ and
  $y=-(c_2/c_1)x-c_0/c_1$ is a non-axis-parallel affine line.
- If $c_3\ne0$, completing the product gives
  $$
  (x+c_1/c_3)(y+c_2/c_3)
  =\frac{c_1c_2-c_0c_3}{c_3^2}=:r.
  $$
  $r=-\Delta/c_3^2\ne0$. The potential pole $x=-c_1/c_3$ makes the left side
  zero and therefore cannot be a point on a curve with $r\ne0$.

Conversely, for a nondegenerate Möbius map
$y=(ax+b)/(cx+d)$ with $ad-bc\ne0$ and no selected denominator zero, set
$c_2=-a$, $c_0=-b$, $c_3=c$, $c_1=d$. Then
$F(x,y)=0$ is exactly its graph and
$c_0c_3-c_1c_2=ad-bc\ne0$. Hence 3 $\Longleftrightarrow$ 4, with both the
factor-degenerate case and the pole handled explicitly.

## 4. P4 — Möbius graph iff simultaneous cross-ratio equality

Let $M(z)=(az+b)/(cz+d)$, $ad-bc\ne0$, and suppose its denominators are
nonzero at the selected arguments. Direct subtraction gives

$$M(s)-M(t)=\frac{(ad-bc)(s-t)}{(cs+d)(ct+d)}.$$

Substituting this identity into the four differences defining cross-ratio
cancels the determinant factors and every denominator factor, yielding

$$\operatorname{CR}(M(x_1),M(x_2),M(x_3),M(x_4))
=\operatorname{CR}(x_1,x_2,x_3,x_4).$$

Thus a Möbius graph implies equality of the paired cross-ratios.

Conversely, three distinct projective points determine a unique Möbius map.
Explicitly,

$$T_{r,s,t}(z)=\frac{(s-t)(z-r)}{(s-r)(z-t)}$$

maps $(r,s,t)$ to $(0,1,\infty)$, so
$M=T_{y_1,y_2,y_3}^{-1}\circ T_{x_1,x_2,x_3}$ is the unique map sending the
first three $x$'s to the first three $y$'s. Möbius invariance gives

$$\operatorname{CR}(y_1,y_2,y_3,M(x_4))
=\operatorname{CR}(x_1,x_2,x_3,x_4).
$$

If the right side equals $\operatorname{CR}(y_1,y_2,y_3,y_4)$, then
$M(x_4)=y_4$: with fixed distinct $y_1,y_2,y_3$, the function
$z\mapsto\operatorname{CR}(y_1,y_2,y_3,z)$ is itself a nonconstant Möbius
function and hence injective on the projective line. In particular $M(x_4)$
is the finite value $y_4$, so $x_4$ is not a pole. This proves
5 $\Longrightarrow$ 4 and completes 4 $\Longleftrightarrow$ 5.

The support has no preferred ordering. Applying the same permutation to both
quadruples either preserves the cross-ratio or applies the same one of its six
fractional-linear transforms to both values. More fundamentally, the ordering-
free statement “one Möbius map sends every paired $x$ to its paired $y$” has
just been proved equivalent to equality. Therefore simultaneous reordering
does not change the support predicate.

## 5. P5 — count, field dependence, and the exact higher-rank boundary

Fix the increasing order $U=\{u_1<\cdots<u_4\}$. An all-distinct support with
row projection $U$ contains exactly one cell in each row, so it defines a
unique injection $\sigma:U\hookrightarrow[n_B]$ by
$(u,\sigma(u))\in S$; injectivity is exactly the distinct-column condition.
Conversely every such $(U,\sigma)$ defines one all-distinct support. This is a
bijection, not a quotient by $4!$, because the order of $U$ is fixed. Applying
§4 to each support proves the displayed formula for $N_{\rm all}$.

For fixed reduced evaluation sets, every part of the formula except the field
equality is combinatorial. Thus the field dependence of this channel is
precisely the dependence of the cross-ratio coincidence set on field
arithmetic. The in-run symmetric example makes this visible while its crossing
population stays fixed: all-distinct counts $8,4,12,4,4$ over
$p=7,11,13,17,31$, but 144 crossing supports at every prime.

For general GRS row counts, write (including nonzero GRS column multipliers)

$$a_u=\alpha_u(1,x_u,\ldots,x_u^{r_A-1})^{\mathsf T},\qquad
b_v=\beta_v(1,y_v,\ldots,y_v^{r_B-1})^{\mathsf T}.$$

Then $a_u\otimes b_v$ is the nonzero scalar $\alpha_u\beta_v$ times the
bi-Vandermonde evaluation vector

$$\bigl(x_u^iy_v^j\bigr)_{0\le i<r_A,\,0\le j<r_B}.$$

Column scalars do not change rank. For a selected set of $m$ pairs,
dependence is exactly column rank $<m$ of this $r_Ar_B\times m$ matrix;
circuitness requires that rank failure plus full column rank on every proper
selected submatrix. Only when the relevant matrix is square can dependence be
expressed by one determinant. In particular, if $r_Ar_B>d+1$, the
size-$(d+1)$ condition is a rank condition (equivalently, simultaneous
vanishing of all maximal minors), not the single $4\times4$ determinant above.
This run proves no higher-rank cross-ratio criterion and no general
higher-rank circuit classification.

## 6. Exact in-run evidence

Every accepted support was checked three independent ways: determinant zero,
cross-ratio equality, and membership in the profile-$(4,4)$ population of a
full exact circuit census. On every determinant-zero support the run also
solved a nonzero bilinear left-null vector, asserted exact vanishing,
$\Delta\ne0$, no selected pole, exact Möbius values, and either the line or
nonzero-$r$ hyperbola normal form.

### T1 — symmetric four-point pair, five primes

For $x=y=(1,2,3,4)$ all 24 bijections were swept at each of
$p\in\{7,11,13,17,31\}$ (120 determinant/CR comparisons), while each full
product census scanned all $\binom{16}{4}=1820$ supports. Set equality held at
every prime:

| $p$ | det-zero | CR-equal | census $(4,4)$ | census $(3,3)$ | total circuits | line / hyperbola |
|---:|---:|---:|---:|---:|---:|---:|
| 7 | 8 | 8 | 8 | 144 | 152 | 2 / 6 |
| 11 | 4 | 4 | 4 | 144 | 148 | 2 / 2 |
| 13 | 12 | 12 | 12 | 144 | 156 | 2 / 10 |
| 17 | 4 | 4 | 4 | 144 | 148 | 2 / 2 |
| 31 | 4 | 4 | 4 | 144 | 148 | 2 / 2 |

This reproduces the frozen $8/4/12/4/4$ all-distinct counts and in particular
$148=144+4$ at GF(31), with zero other circuit profiles. The identity and
reversal pairings were pure-tensor circuit ACCEPTs at every prime, both affine
lines. At GF(13), the supplied pairing
$(1,1),(2,3),(3,4),(4,2)$ is a genuine hyperbola circuit: the run found
$(c_0,c_1,c_2,c_3)=(5,1,6,1)$,
$\Delta=12\ne0$, and
$(x+1)(y+6)=1$; all four denominators $(2,3,4,5)$ are nonzero.

### T2 — asymmetric five-point pair, two primes

For $x=(1,2,3,4,5)$ and $y=(2,3,5,7,11)$, every
$\binom54\binom54 4!=600$ all-distinct pairing was checked at each prime
(1,200 total); each full product census scanned $\binom{25}{4}=12,650$
supports. Again the three sets were exactly equal:

| $p$ | det-zero = CR-equal = census $(4,4)$ | census $(3,3)$ | total | line / hyperbola |
|---:|---:|---:|---:|---:|
| 13 | 84 | 900 | 984 | 6 / 78 |
| 17 | 40 | 900 | 940 | 4 / 36 |

The change $84\to40$ at fixed point labels is direct finite evidence for the
cross-ratio field-dependence mechanism; it is not used to prove the algebra.

### T3 — controls in both directions

- **Cross-ratio mismatch REJECT (GF(13)).** The pairing
  $x=(1,2,3,4)$, $y=(1,2,4,3)$ has
  $\operatorname{CR}_x=10\ne4=\operatorname{CR}_y$ and determinant $6\ne0$;
  exact rank rejects it as a circuit.
- **Non-tensor plant REJECT.** The pristine identity support has rank 3 and
  equal cross-ratios. Replacing actual $h(0,0)$ by
  $h(0,0)+h(0,1)$ while retaining its labels raises actual rank to 4; the
  label criterion still says ACCEPT, so the fail-loud disagreement rejects
  the non-tensor instrument.
- **Duplicated-column REJECT.** Cloning one factor column drops measured spark
  from 3 to 2 and exposes dependent pair $(0,4)$.
- **Concat/Kronecker guard.** Independently generated 2-row by 3-row pure
  tensors have exact multiplicative coordinates and dimension 6; planted
  block concatenation has dimension 5 and different coordinates.

Runtime controls asserted `PYTHONDONTWRITEBYTECODE=1`,
`sys.dont_write_bytecode`, all thread caps 1, and nice value 15. The exact run
used CPU 0.466 s and wall 0.491 s, below the preregistered 600/900 s caps.

## 7. Scope and nonclaims

**Certified:** for two 2-row GRS factors with distinct evaluation points, the
all-distinct size-four circuit predicate is exactly determinant zero, exactly
one nondegenerate bilinear/Möbius graph, and exactly simultaneous cross-ratio
equality; the injection count formula is exact; field dependence is precisely
cross-ratio-coincidence dependence. The five-prime and two-factor-pair sweeps
are exact finite corroboration, not the source of the proof.

**Not certified:** a simpler closed formula evaluating the coincidence count
for arbitrary point sets; spark-2 branches; repeated-coordinate crossing or
fiber populations (separate theorems); a cross-ratio rule when $r_Ar_B>4$;
general higher-rank circuit classification beyond the bi-Vandermonde rank
statement; three or more factors; asymptotics or characteristic-specific
uniformity beyond the algebraic theorem.

## 8. Defect disclosure and verdict inputs

**D1 (launch wrapper only, before Python).** The first launch command used the
shell builtin `exec`, unavailable in the harness command wrapper, and exited
127 with `error: command not found: exec` before Python started. It created no
run output and consumed no parameter family. The command was changed only by
removing `exec`; the instrument and parameters were unchanged, and the entire
T1–T3 matrix then ran from the beginning. This failure and remediation are
preserved in `defect_log.json` and the session transcript. No mathematical or
instrument assertion failed; no pin or parameter was changed.

P1–P5 are discharged above, including three-column minimality, the degenerate
condition $c_1c_2=c_0c_3$, selected poles, the converse Möbius construction,
and the higher-rank boundary. All three set populations are equal at all seven
configurations; all supplied anchors reproduce; known-true ACCEPTs pass; every
planted REJECT fires; runtime and budget controls pass. Under the preregistered
promotion rule this supports exactly one verdict: **FROZEN-CERTIFIED**.
