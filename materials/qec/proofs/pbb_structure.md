# Structure theorems for perturbed bivariate-bicycle (PBB) codes

Every statement below is labelled with its epistemic status:
**[proved]** analytic proof given here · **[exact computation]** certified by
exact GF(2) arithmetic or a solver that reported proven optimality/infeasibility ·
**[statistical]** Monte Carlo with a stated confidence interval.

---

## 0. Setting and notation

Let

$$R \;=\; \mathbb{F}_2[x,y]\big/\big(x^{\ell}-1,\;y^{m}-1\big),\qquad
\dim = \ell m,\qquad n = 2\,\dim .$$

Identify $f\in R$ with the $\dim\times\dim$ matrix $M(f)$ of the regular
representation, indexed by $G=\mathbb Z_\ell\times\mathbb Z_m$ with
$M(x^ay^b)_{g,\,g+(a,b)}=1$.  Two facts are used constantly:

* **(R1)** $M(f)M(h)=M(h)M(f)$ — $R$ is commutative.
* **(R2)** $M(f)^{\mathsf T}=M(f^{*})$ where $f^{*}(x,y)=f(x^{-1},y^{-1})$;
  transposition is the ring involution.

For $A,B,C,D\in R$ write

$$H_X=[\,A\;\;B\,],\qquad H_Z=[\,B^{\mathsf T}\;\;A^{\mathsf T}\,],\qquad P=[\,C\;\;D\,],$$

each of shape $\dim\times n$.  The **PBB code** $Q=\mathrm{PBB}(A,B,C,D)$ has
symplectic check matrix

$$H \;=\; \begin{pmatrix} H_X & P\\[2pt] 0 & H_Z\end{pmatrix},$$

with the convention $v=(x\,|\,z)$ and
$\langle (x|z),(x'|z')\rangle_s = x\!\cdot\!z' + z\!\cdot\!x' \bmod 2$.
Setting $C=D=0$ recovers the CSS bivariate-bicycle code
$\mathrm{BB}(A,B)$, henceforth the **parent**.

---

## 1. Commutation condition

> **Lemma 0.** **[proved]** The rows of $H$ pairwise commute if and only if
> $$M \;:=\; AC^{\mathsf T}+BD^{\mathsf T}\quad\text{is \emph{symmetric}.}$$

*Proof.*  Write $H=(H_x\,|\,H_z)$ with
$H_x=\begin{pmatrix}A&B\\0&0\end{pmatrix}$,
$H_z=\begin{pmatrix}C&D\\B^{\mathsf T}&A^{\mathsf T}\end{pmatrix}$.
Then
$$H_xH_z^{\mathsf T}=\begin{pmatrix}A&B\\0&0\end{pmatrix}
\begin{pmatrix}C^{\mathsf T}&B\\ D^{\mathsf T}&A\end{pmatrix}
=\begin{pmatrix}AC^{\mathsf T}+BD^{\mathsf T} & AB+BA\\ 0&0\end{pmatrix}
=\begin{pmatrix}M&0\\0&0\end{pmatrix},$$
using $AB=BA$ from **(R1)**.  Since
$H_zH_x^{\mathsf T}=(H_xH_z^{\mathsf T})^{\mathsf T}$, the validity condition
$H_xH_z^{\mathsf T}+H_zH_x^{\mathsf T}=0$ reads $M+M^{\mathsf T}=0$, i.e.
$M=M^{\mathsf T}$ over $\mathbb F_2$. $\blacksquare$

**Remark (agreement with the primary source).**  In polynomial form, with
$M=M(ac^{*}+bd^{*})$ and **(R2)**, the condition says $ac^{*}+bd^{*}$ is
invariant under $f\mapsto f^{*}$ — *not* that it vanishes.  arXiv:2606.02418v1
§III.2 states this correctly and verbatim: "All rows commute if and only if
$(AC^{\top}+BD^{\top})\bmod 2$ is symmetric over $\mathbb{F}_2$ (this is the
only nontrivial commutativity condition; commutativity between block 1 and
block 2 is automatic from BB ring commutativity)."  The reference
implementation `qcode-discovery/evaluation/pbb_code.py:71-73` tests the same
predicate.  Our derivation reproduces it independently.  **[exact computation]**
all 368 catalogue codes satisfy symmetry and our independent construction
rebuilds every one with zero commutation failures.

---

## 2. The pure-$Z$ sector is inherited from the parent

> **Proposition 1.** **[proved]** For $w\in\mathbb F_2^{\,n}$, the operator
> $(0\,|\,w)$ commutes with every row of $H$ **iff** $H_X w^{\mathsf T}=0$.
> Hence the pure-$Z$ centraliser of $Q$ equals $\ker H_X$ — *identical to that
> of the parent* $\mathrm{BB}(A,B)$, and independent of $C,D$.

*Proof.*  $\langle (0|w),(h_x|h_z)\rangle_s = w\cdot h_x$.  Block-1 rows have
$h_x = (H_X)_i$; block-2 rows have $h_x=0$. $\blacksquare$

> **Proposition 2.** **[proved]** The pure-$Z$ elements of $\operatorname{rowspace}(H)$ are
> $$S_Z \;=\; \operatorname{rowspace}(H_Z)\;+\;\{\,uP \;:\; u\in \mathrm{LN}(H_X)\,\},
> \qquad \mathrm{LN}(H_X)=\{u: uH_X=0\}.$$

*Proof.*  A general element is $(u,w)H=(uH_X,\;uP+wH_Z)$; its $x$-part vanishes
iff $uH_X=0$. $\blacksquare$

> **Proposition 3 (dimension bookkeeping).** **[proved]** Put
> $r_X=\operatorname{rank}H_X$, $r_Z=\operatorname{rank}H_Z$ and
> $$\boxed{\;\delta \;:=\; \dim S_Z - r_Z \;\in\; [\,0,\;\dim - r_X\,].\;}$$
> Then $\operatorname{rank}H = r_X+r_Z+\delta$ and
> $$k(Q) \;=\; k\big(\mathrm{BB}(A,B)\big)\;-\;\delta .$$

*Proof.*  The map $\pi:\operatorname{rowspace}(H)\to\mathbb F_2^{\,n}$,
$(u,w)H\mapsto uH_X$, has image $\operatorname{rowspace}(H_X)$ of dimension
$r_X$ and kernel exactly $S_Z$ of dimension $r_Z+\delta$.  So
$\operatorname{rank}H=r_X+r_Z+\delta$ and $k=n-\operatorname{rank}H$.  The
parent has $\operatorname{rank}=r_X+r_Z$. $\blacksquare$

**[exact computation]** Verified on all 368 catalogue codes: the catalogue's
$k$, the formula $k_{\mathrm{BB}}-\delta$, and a direct symplectic rank
computation agree in **368/368** cases.  Observed distribution
$\delta\in\{0:213,\;2:115,\;4:34,\;6:2,\;8:2,\;10:2\}$.

---

## 3. The distance ceiling

> **Theorem 1 (PBB distance ceiling).** **[proved]** If $\delta=0$ — equivalently
> $k(Q)=k(\mathrm{BB}(A,B))$ — then the nontrivial pure-$Z$ logical operators of
> $Q$ are **exactly** those of the parent $\mathrm{BB}(A,B)$, and therefore
> $$d(Q)\;\le\;d_Z\big(\mathrm{BB}(A,B)\big).$$

*Proof.*  $\delta=0$ means $S_Z=\operatorname{rowspace}(H_Z)$.  By Prop. 1 the
two codes have the same pure-$Z$ centraliser $\ker H_X$; by Prop. 2 they have
the same pure-$Z$ stabiliser subgroup.  The quotients therefore coincide.  Let
$w$ realise $d_Z(\mathrm{BB}(A,B))$; then $(0|w)$ is a nontrivial logical of
$Q$ of symplectic weight $|w|$, so $d(Q)\le|w|$. $\blacksquare$

> **Proposition 4 (BB self-duality).** **[proved]** Let $\sigma$ permute the
> $n=2\dim$ qubits by $(\mathrm L,g)\mapsto(\mathrm R,-g)$,
> $(\mathrm R,g)\mapsto(\mathrm L,-g)$.  Then $\sigma$ carries
> $\operatorname{rowspace}(H_X)$ onto $\operatorname{rowspace}(H_Z)$ and back.
> Consequently **every** bivariate-bicycle code satisfies $d_X=d_Z$.

*Proof.*  Swapping the two halves turns $[A\;B]$ into $[B\;A]$; the lattice
inversion $g\mapsto-g$ then replaces each circulant by its transpose via
**(R2)**, giving $[B^{\mathsf T}\;A^{\mathsf T}]=H_Z$. $\blacksquare$

**[exact computation]** Checked on all 7 Bravyi instances and all 368
catalogue parents: **375/375**.

> **Corollary 1 (parent domination).** **[proved]** If $\delta=0$ then the
> parent $P=\mathrm{BB}(A,B)$ is a CSS bivariate-bicycle code with
> $$n(P)=n(Q),\qquad k(P)=k(Q),\qquad d(P)=d_Z(P)\;\ge\;d(Q).$$
> $P$ therefore **weakly dominates $Q$ in $[[n,k,d]]$**, while having maximum
> check weight $|A|+|B|$ and admitting a purely CSS syndrome circuit.

*Proof.*  $n,k$ equal by Prop. 3 with $\delta=0$; $d(P)=\min(d_X,d_Z)=d_Z(P)$
by Prop. 4; then Theorem 1. $\blacksquare$

**[exact computation]** Corollary 1 applies to **213 of the 368** published PBB
codes (57.9 %), *including all 14 members of the headline
$[[144,12,12]]$ family*.  For the representative `12_6_0187` we certified by
CP-SAT (status `OPTIMAL`) that the PBB minimum pure-$Z$ logical weight is
$12$ and the parent's is also $12$ — Theorem 1's equality prediction, confirmed.
Over those 213 codes the parent's maximum check weight is $6$ in every case,
versus $7$–$16$ for the PBB code.

> **Theorem 2 (CSS shadow; arbitrary $\delta$).** **[proved]** For any
> $Q=\mathrm{PBB}(A,B,C,D)$ define the CSS code $Q'$ with $X$-checks $H_X$ and
> $Z$-checks $S_Z$.  Then $H_XS_Z^{\mathsf T}=0$, and
> $$n(Q')=n(Q),\qquad k(Q')=k(Q),\qquad d(Q)\;\le\;d_Z(Q').$$

*Proof.*  $S_Z\subseteq\ker H_X$ by Prop. 1, so $Q'$ is a legitimate CSS code.
$k(Q')=n-r_X-\dim S_Z=k(Q)$ by Prop. 3.  The nontrivial pure-$Z$ logicals of
$Q$ are the $Z$-logicals of $Q'$ by Props. 1–2, giving the inequality exactly
as in Theorem 1. $\blacksquare$

**Interpretation, stated at exactly the strength proved.**  Theorem 2 gives a
*pure-$Z$ ceiling*: $d(Q)\le d_Z(Q')$.  It does **not** give
$d(Q)\le d(Q')=\min(d_X(Q'),d_Z(Q'))$, because nothing here bounds $d_X(Q')$
from below.  So Theorem 2 alone does **not** say the shadow dominates $Q$.

Full $[[n,k,d]]$ domination is proved only in the $\delta=0$ case
(Corollary 1), where $Q'$ is the parent BB code and Proposition 4 supplies the
missing half via $d_X=d_Z$.  That covers **213 of the 368** catalogue codes,
including all 14 of the $[[144,12,12]]$ family.

For $\delta>0$ the gap is real, and EXP-027's strict-provenance audit now
gives it a certified two-sided shape.  `phase2_58` and `phase2_60`
($[[72,4,6]]$, $\delta=4$) beat their own CSS shadow ($d_X(Q')=4$), **and**
beat their parent outright: the parent is a $[[72,8,4]]$ BB code with
$d=4$ **exact** (fresh `exact_distance_css`, OPTIMAL in every sector,
witness verified), while the perturbed code has $d=6$ **exact**
(`exact_distance_symplectic`, OPTIMAL, witness verified; both rebuilt from
raw terms and re-solved with distinct seeds — `results/raw/exp027_reversal_*.json`).
So a $\delta=4$ perturbation **strictly increased the distance** while halving
$k$: parent and child are incomparable in $[[n,k,d]]$.  (Both remain dominated
by the certified Bravyi $[[72,12,6]]$ — $k=12$ against 4 at equal $n$ and $d$ —
so the CSS *envelope* still stands; what falls is only the idea that
perturbations never buy distance.)

The honest global statement, with the full EXP-027 ledger over all 155
$\delta>0$ rows (66 DOMINATION\_PROVED, 2 CERTIFIED\_REVERSAL, 87 UNDECIDED —
literature values never accepted as certificates): the perturbation can never
*raise* $k$ above the parent's (Prop. 3, all $\delta$), and $d$ is capped by a
same-$(n,k)$ CSS code's $Z$-distance (Theorem 2, all $\delta$); full
$[[n,k,d]]$ **domination holds for $\delta=0$ (Corollary 1) and is FALSE in
general for $\delta>0$**.  EXP-036 (2026-08-16) closes part of the undecided
tail with the Theorem-F CDCL method: of the 87, all 29 at $n=108/144$ plus 5
at $n=180/360$ are decided — 29 further dominations and **5 further certified
reversals**
(`phase2_71`, `phase2_72`, `phase2_88`, `9_6_0183`, `12_6_0217`), each a
verified witness upper bound plus a replayed UNSAT lower bound.  Running
totals: **107 dominations / 7 reversals / 41 undecided** ($n\in\{180,360\}$,
in flight).  All seven certified reversals are CSS-dominated at equal $n$ — EXP-036's five
and EXP-027's two, each against the strongest certified-exact CSS code at that
length (machine-checked, `results/processed/exp036_envelope_check.json`) — so the
envelope question below stays open only in its universal form.

## Theorem G (survival criterion, EXP-038): why perturbation almost never buys distance

EXP-036 decided its rows one expensive UNSAT at a time and said nothing about
*why* reversals are rare.  Theorem G answers that structurally, and turns the
answer into a domination test that needs no search on the PBB at all.

Fix a parent CSS BB code $P$ with $H_X=[A\;B]$, $H_Z=[B^{T}\;A^{T}]$ on
$n=2\ell m$ qubits, and the PBB code $Q$ built on it,
$$H_Q=\begin{pmatrix}A&B&C&D\\0&0&B^{T}&A^{T}\end{pmatrix}.$$
Write $S_Z=\operatorname{rowspace}(H_Z)\subseteq\mathrm{GF}(2)^{n}$ for the
parent's $Z$-stabilizers and define the **dressing space**
$$\Delta=\bigl\{\lambda[C\;D]\;:\;\lambda[A\;B]=0\bigr\}\subseteq\mathrm{GF}(2)^{n},$$
the $Z$-parts of those first-block combinations whose $X$-part cancels.

**(i) The pure-$Z$ centralizer is unchanged.**  The symplectic product of
$(0|z)$ with a second-block row $(0|b)$ vanishes identically, and with a
first-block row $(a|c)$ it equals $a\cdot z$.  Hence
$\{z:(0|z)\in C(Q)\}=\ker[A\;B]=\{z:(0|z)\in C(P)\}$.  The perturbation does
not touch this space.

**(ii) The pure-$Z$ stabilizers grow by exactly $\Delta$, and
$k_Q=k_P-\dim(\Delta+S_Z)/S_Z$.**  A combination with first-block coefficients
$\lambda$ and second-block coefficients $\mu$ has $X$-part $\lambda[A\;B]$, so
it is pure $Z$ iff $\lambda\in\ker_L[A\;B]$, and then its $Z$-part is
$\lambda[C\;D]+\mu H_Z$.  So the pure-$Z$ stabilizer group is $S_Z+\Delta$.
For the dimension count, the left kernel of $H_Q$ is
$$\{(\lambda,\mu):\lambda[A\;B]=0,\ \lambda[C\;D]=\mu H_Z\},$$
whose dimension is
$\dim\{\lambda\in\ker_L[A\;B]:\lambda[C\;D]\in S_Z\}+\dim\ker_L H_Z$
$=\dim\ker_L[A\;B]-\dim(\Delta+S_Z)/S_Z+\dim\ker_L H_Z$.
Since $\operatorname{rank}H_P=\operatorname{rank}H_X+\operatorname{rank}H_Z$
(disjoint column blocks), this gives
$\operatorname{rank}H_Q=\operatorname{rank}H_P+\dim(\Delta+S_Z)/S_Z$, and
$k=n-\operatorname{rank}H$ yields the claim. $\square$

This is a *proved* identity, and it is independently machine-verified on all
368 catalogue rows (`tests/test_dimension_identity_368.py`): $k_P-k_Q=\dim\Delta$ with
zero mismatches, including the 213 $\delta=0$ rows where $\Delta=0$ and $k$ is
unchanged.  It also gives a free cross-check on both dimensions, which is how a
stacked-rank error that would have *inflated* $k_C$ was caught.

**(iii) Survival bounds the distance.**  Let $z\in\ker[A\;B]$ with
$\operatorname{wt}(z)=d_Z(P)$ and $z\notin S_Z+\Delta$.  By (i) $z$ is in
$C(Q)$; by (ii) it is not in $Q$'s stabilizer group; so it is a nontrivial
logical of $Q$ of weight $d_Z(P)$, whence
$$d_Q\le d_Z(P).$$
The survivor *is* the certificate — an explicit vector, re-checked through two
independent GF(2) paths.  No weight-bounded search on $Q$ is performed.

**(iv) Counting form.**  Contrapositive of (iii): if $d_Q>d_Z(P)$ then $\Delta$
must absorb *every* minimum-weight $Z$-logical of $P$, so it contains their
whole span modulo $S_Z$.  If that span has rank $t$ then
$$d_Q>d_Z(P)\ \Longrightarrow\ k_Q\le k_P-t .$$
Equivalently: exhibiting $\dim\Delta+1$ independent minimum-weight $Z$-logicals
of the parent certifies $d_Q\le d_Z(P)$ outright.  BB translations are code
automorphisms preserving weight and the logical quotient, so a single witness
generates its entire $\ell m$-element orbit for free; $t$ is then the orbit's
rank modulo $S_Z$.

### What this explains

The empirical pattern that EXP-036 found but could not account for.  Across the
155 $\delta>0$ catalogue rows the ratio $k_Q/k_P$ takes exactly three values —
$1/2$ (58 rows), $3/4$ (54 rows), $5/6$ (43 rows) — and **all seven certified
reversals lie in the $k$-halving class**, none among the 97 rows above it.
Theorem G (iv) is the reason: a reversal needs $\dim\Delta=k_P-k_Q\ge t$, and
for these codes the minimum-weight orbit span $t$ is large, so only the maximal
dressing spaces can absorb it.

Sharper still, and machine-checked on every certified reversal: the margin
$t-\dim\Delta$ is **exactly zero**.  The dressing space is precisely large
enough to absorb the minimum-weight orbit span and no larger — reversals sit on
the knife edge of the criterion rather than comfortably inside it.  That is why
they exist at all, and why they are worth so little: each pays a factor-two loss
in $k$ for the smallest distance step the parity of these codes permits, $+2$.

**Scope.**  (i)-(iv) are proved for every parent/PBB pair of the stated shape,
with no restriction on $C,D$.  (iii) concludes $d_Q\le d_Z(P)$; converting that
into domination of the parent's *code* distance $d_P=\min(d_X,d_Z)$ additionally
needs $d_Z(P)\le d_X(P)$, which is why the experiment reports the certified
inequality $d_Q\le d_Z(P)$ rather than silently comparing against $d_P$.

---

## Theorem H (parent-level no-go, EXP-039): one certificate closes an entire family

Theorem G (iv) decides one perturbation at a time: $t$ was the orbit rank of the
logicals a *given* $\Delta$ absorbs.  The next observation removes $[C\;D]$ from
the statement entirely.

> **Lemma 1 (module structure).** **[proved]** Let
> $R=\mathbb F_2[x,y]/(x^{\ell}-1,y^{m}-1)$ act on $\mathbb F_2^{\,\ell m}$ and
> on the two-block $Z$-sector $\mathbb F_2^{\,n}=R\oplus R$ by simultaneous
> translation.  Then:
> * $S_Z=\operatorname{rowspace}[B^{\mathsf T}\;A^{\mathsf T}]$ is an
>   $R$-submodule;
> * $L:=\{\lambda : \lambda A=\lambda B=0\}$ is an **ideal** of $R$, hence an
>   $R$-submodule of $\mathbb F_2^{\,\ell m}$;
> * $\Delta=L\cdot[C\;D]$ is an $R$-submodule of the $Z$-sector, for every
>   $[C\;D]$.

>*Proof.*  $A$ and $B$ are $R$-linear maps, so their rowspaces, kernels and
images are $R$-submodules; in particular translating a row of $H_Z$ by
$x^{a}y^{b}$ yields another row.  If $\lambda A=\lambda B=0$ then
$(x^{a}y^{b}\lambda)A=x^{a}y^{b}(\lambda A)=0$ and likewise for $B$, using
commutativity of $R$; so $L$ is closed under multiplication by $R$.  Finally
$\Delta=\{\lambda[C\;D]:\lambda\in L\}$ is the image of the $R$-module
homomorphism $L\to R\oplus R,\ \lambda\mapsto\lambda[C\;D]$ (here $C,D$ act by
multiplication, which commutes with translation).  $\square$

Consequently, since translations preserve weight, commute with $[A\;B]$, and
stabilise $S_Z$: **absorbing one minimum-weight $Z$-logical of the parent
absorbs its entire $\ell m$-element translation orbit.**  Define, from the
parent alone,
$$M(P)=\operatorname{span}\{\text{translation orbits of all minimum-weight }
Z\text{-logicals of }P\},\qquad T(P)=\dim\frac{M(P)+S_Z}{S_Z}.$$

> **Theorem H.** **[proved]** For every parent $P$ and every valid perturbation
> $[C\;D]$,
> $$d_Q>d_Z(P)\ \Longrightarrow\ k_Q\le k_P-T(P).$$
> In particular $T(P)=k_P$ closes the family: no PBB over $P$ keeps a logical
> qubit while exceeding $d_Z(P)$.

>*Proof.*  By Theorem G (iii), $d_Q>d_Z(P)$ forces *every* minimum-weight
$Z$-logical $z$ of $P$ into $S_Z+\Delta$.  By Lemma 1 the absorbed set is a
submodule containing each orbit, so $M(P)\subseteq S_Z+\Delta$, i.e.
$\bar M\subseteq\bar\Delta$; hence $\dim\bar\Delta\ge T(P)$ and Theorem G (ii)
gives $k_Q=k_P-\dim\bar\Delta\le k_P-T(P)$.  $\square$

$T(P)$ is computed *exactly* — not bounded — by a SAT enumeration that
terminates in a certified UNSAT (`src/qec_research/codes/pbb_nogo.py`,
EXP-039): grow $W=S_Z+M$ by asking for weight-$\le d_Z$ elements of
$\ker[A\;B]$ outside $W$ (nontriviality vs $W$ encoded through $W^{\perp}$);
each hit is necessarily a minimum-weight logical and donates its full orbit;
UNSAT certifies $W$ contains them all.  Independently cross-checked by a
brute-force MITM implementation with no SAT solver (EXP-041: 12/12 in-budget
parents agree exactly, including three $T=k_P$ closures).

Machine results: $T$ exact for 134/202 catalogue parents — exhaustively all
117 at $n\le 144$ — with **61 family-closed**, including the Gross code
($A=x^3+y+y^2$, $B=y^3+x+x^2$, $T=12=k_P$) and all eleven catalogue
$[[144,12,12]]$ parents; 249/368 rows capped a priori with no per-row search.

---

## Theorem I (forced saturation, EXP-039/040): the trade is exact

All seven certified reversals sit at $k_Q=k_P-T$ with margin exactly zero
(§Theorem G, "sharper still").  That is no coincidence, and the small-lattice
probe of EXP-040 (26,898 unseen $\delta>0$ perturbations over 64 parents:
$\dim\bar\Delta<T$ 23,634 times, $=T$ 3,264 times, $>T$ **zero** times; all
496 strict distance increases at $=T$) demanded the reason.  It is two ranks.

> **Lemma 2 (dimension identity for the left kernel).** **[proved]** For every
> CSS BB parent, with $L=\{\lambda:\lambda A=\lambda B=0\}$,
> $$k_P \;=\; 2\dim L .$$

>*Proof.*  Everything is rank-nullity.  $\operatorname{rank}H_X=\ell m-\dim L$,
> since $L$ is exactly the left nullspace of $H_X=[A\;B]$.  Dually, the left
> nullspace of $H_Z=[B^{\mathsf T}\;A^{\mathsf T}]$ is
> $\{c : cB^{\mathsf T}=0,\ cA^{\mathsf T}=0\}$, i.e. the transpose of
> $K=\ker A\cap\ker B$ (right kernels), so
> $\operatorname{rank}H_Z=\ell m-\dim K$.  Therefore
$$k_P = n-\operatorname{rank}H_X-\operatorname{rank}H_Z
      = 2\ell m-(\ell m-\dim L)-(\ell m-\dim K)
      = \dim L+\dim K .$$
Finally $\dim L=\dim K$: the coordinate-reversal map $J$ (which conjugates each
cyclic shift $S$ to $S^{\mathsf T}=S^{-1}$) is a linear isomorphism carrying
$\ker A$ onto $\ker A^{\mathsf T}$ and $\ker B$ onto $\ker B^{\mathsf T}$, so
the two intersections have equal dimension.  Hence $k_P=2\dim L$.  $\square$

*(Machine-checked: $\dim L=k_P/2$ on all 134 certified parents, zero
exceptions; `tests/test_pbb_theorems.py`.)*

> **Lemma 3 (dressing ceiling).** **[proved]** For every valid perturbation,
> $$\dim\bar\Delta \;\le\; \dim L \;=\; k_P/2 .$$

>*Proof.*  $\Delta$ is the image of the linear map $L\to\mathbb F_2^{\,n}$,
$\lambda\mapsto\lambda[C\;D]$, so $\dim\Delta\le\dim L$; quotienting by $S_Z$
cannot increase dimension.  Apply Lemma 2.  $\square$

> **Theorem I (forced saturation).** **[proved]** If $T(P)\ge k_P/2$ and
> $d_Q>d_Z(P)$, then
> $$\dim\bar\Delta = T(P) \quad\text{and hence}\quad k_Q = k_P - T(P)
> \quad\text{exactly}.$$

>*Proof.*  Theorem H gives $\dim\bar\Delta\ge T$; Lemma 3 gives
$\dim\bar\Delta\le k_P/2\le T$.  Squeeze.  Theorem G (ii) converts dimensions
into $k_Q$.  $\square$

**Scope on the catalogue.**  $T\ge k_P/2$ holds on 133 of the 134 certified
parents — the single exception is `9a7638586033` ($n=144$, $k_P=12$, $T=4$,
not family-closed), where Theorem I degrades to the two-sided bound
$4\le\dim\bar\Delta\le 6$ under increase.  On every other certified parent,
**the exact trade law is now a theorem**: a distance increase pays precisely
$T(P)$ logical qubits, never fewer, never more.  The EXP-040 probe's parents
all had $T\in\{k_P/2,\,k_P\}$, so its 496 strict increases at $\dim\bar\Delta=T$
were theorem-forced, not lucky.

The residual conjecture — saturation when $T<k_P/2$ — survives only on the
exceptional class, pending the $n=180/360$ sweep.

---

## 4. Circuit-level consequences

> **Lemma C1 (measurement-correctness criterion).** **[proved]** In a
> one-ancilla-per-check circuit — ancilla $a$ prepared in $|+\rangle$, gates
> $\mathrm{C}\text{-}P^{a}_j$ applied for $j\in\operatorname{supp}(g_a)$, ancilla
> measured in the $X$ basis — let
> $J(a,b)=\{\,j:\ P^{a}_j \text{ and } P^{b}_j \text{ anticommute}\,\}$.
> The circuit measures the intended generators **iff** for every pair $(a,b)$
> $$\#\{\,j\in J(a,b) : a \text{ acts before } b\,\}\ \text{ is even.}$$

*Proof.*  $\mathrm{C}\text{-}P_{a,j}$ and $\mathrm{C}\text{-}Q_{b,j}$ commute
when $[P,Q]=0$; when they anticommute,
$\mathrm{C}\text{-}P_{a,j}\,\mathrm{C}\text{-}Q_{b,j}
=\mathrm{C}\text{-}Q_{b,j}\,\mathrm{C}\text{-}P_{a,j}\cdot \mathrm{CZ}_{a,b}$,
because the two orderings differ by $PQ$ versus $QP=-PQ$ on the
$|1_a1_b\rangle$ component.  Reordering the schedule into the canonical
"all of $a$ first" order therefore emits one $\mathrm{CZ}_{a,b}$ per element of
$J(a,b)$ whose order is inverted.  Since $\mathrm{CZ}^2=I$, the accumulated
operator is trivial iff that count is even.  $|J(a,b)|$ is even because $g_a$
and $g_b$ commute, so "$b$ first is even" $\iff$ "$a$ first is even".
$\blacksquare$

**[exact computation]** A König edge-colouring of the Steane code violates the
criterion on 5 pairs and fires 1402 detectors in 200 noiseless shots; the
CP-SAT schedule that enforces it fires **0** detectors in 4000 noiseless shots
on Steane $[[7,1,3]]$, Shor $[[9,1,3]]$ and the non-CSS $[[5,1,3]]$ code.

> **Proposition C2 (depth lower bound).** **[proved]** Any such circuit needs at
> least $T\ \ge\ \max\big(\max_a|\operatorname{supp}(g_a)|,\ \max_j \deg j\big)$
> two-qubit layers.  For a PBB code with $C,D$ not both zero, some block-1 check
> has weight $\ge |A|+|B|+1$, so $T\ge 7$ whenever $|A|=|B|=3$.

*Proof.*  An ancilla, and a data qubit, can each take part in at most one
two-qubit gate per layer. $\blacksquare$

> **Theorem C3 (depth-7 is forced within the translation-invariant class, and
> achieved, for weight-6 BB codes).**
> **[exact computation, TI schedule class]** For $\mathrm{BB}$ codes
> $[[72,12,6]]$, $[[108,8,10]]$, $[[144,12,12]]$ and $[[288,12,18]]$, CP-SAT
> proves **INFEASIBLE** at $T=6$ under the criterion of Lemma C1 **within the
> translation-invariant (orbit) schedule class** (`exp004` passes
> `symmetry_orbits`; every entry of `results/raw/exp004_schedules.json` records
> `"mode": "orbit"`), and returns a valid schedule at $T=7$ (verified
> independently, 0 noiseless detector firings).  A depth-7 schedule is a
> schedule, so $T=7$ is achieved unconditionally; Proposition C2 gives the
> unconditional bound $T\ge6$.  Our artifacts do **not** exclude an
> unrestricted depth-6 schedule.  ASC (arXiv:2603.21499) independently reports
> certifying that no depth-6 syndrome circuit exists for the IBM BB instances
> — an external result we did not reproduce.

Within the TI class this *derives* the depth-7 syndrome cycle of
arXiv:2308.07915 from the correctness criterion rather than assuming it; the
unrestricted 6-vs-7 question is closed only externally (ASC).  This scoping was
a silent-restriction defect of the same class as FR-005/FR-014, caught by the
2026-08-13 adversarial audit.

> **Theorem C4 (depth separation for the $[[144,12,12]]$ family).**
> **[proved + exact computation]** Every one of the 14 catalogue PBB
> $[[144,12,12]]$ codes has maximum check weight $\ge 8$ in its **published
> generating set**.  By Proposition C2 any one-ancilla schedule **that measures
> those generators** — translation-invariant or not — needs $T\ge 8$ two-qubit
> layers, whereas the CSS Gross code achieves $T=7$ (Theorem C3).

> **Theorem C6 (basis-independent generator-weight bound).**
> **[proved + exact computation, EXP-023]**  For each of the 14 catalogue PBB
> $[[144,12,12]]$ codes: (i) the pure-$Z$ subgroup
> $S_Z=\{v\in S:\ X\text{-part}(v)=0\}$ has GF(2) rank **exactly 66**
> (three independent rank computations per code); (ii) **no element of $S$
> with nonzero $X$-part has symplectic weight $\le 7$** — decided by a
> certified-complete CP-SAT enumeration of all weight-$\le7$ $X$-codewords
> plus exact GF(2) coset-budget tests and positive SAT controls.  An independent
> artifact rebuilt every code, repeated NumPy/bitset GF(2) ranks, and reproduced
> the bounded row-combination and coset checks in physical 288-bit space.  (A
> second monolithic CP-SAT encoding was run only on benchmark member
> `12_6_0193`, where it agreed; it was not run on the other 13.)
> $V_7\subseteq S_Z$, so $\operatorname{rank}(V_7)=66<132=\operatorname{rank}(S)$,
> and **every generating set contains a generator of weight $\ge 8$**.
> By Proposition C2, every one-ancilla syndrome round — for **any** choice of
> measured generating set — needs $T\ge8$.  Theorem C4's published-generator
> caveat is removed.

In the terminology of arXiv:2601.19848 this determines the *optimal code
weight*: $W\ge 8$ for all 14 members (exactly $8$ for `12_6_0193`,
`12_6_0194` and two further members; three others have interval $[8,9]$ and
the remaining seven have interval $[8,10]$), versus $W=6$ exactly for the CSS
Gross code ($\operatorname{rank}(V_5)=0$; the published weight-6 rows span all
132 dimensions).  These upper-bound differences do not affect the depth bound.
Evidence: `results/processed/exp023_light_gensets.json`,
independent re-verification in `results/raw/exp023_independent_verification.json`,
`tests/test_exp023_light_gensets.py`.  The earlier decision procedure sketch
(`experiments/exp022_basis_independence.py`) is superseded by EXP-023.

**Remaining scope note.**  The depth consequence stays within the one-ancilla
model of Proposition C2 (one two-qubit gate per check-qubit incidence,
no ancilla-ancilla forwarding), exactly as before.  EXP-024 addresses the
multi-ancilla question separately, and only structurally so far: its specified
two-ancilla 4+4 cat construction has a verified depth-9 TI circuit, exact
total-two-qubit-depth-7 infeasibility, and verdict `NEGATIVE_STRUCTURAL` in
`results/processed/exp024_cat_extraction_structural.json`; its **canonical
Monte Carlo benchmark was never executed**, and the production CLI fails
closed until the frozen supervised v2 driver exists (FR-020).

Concretely, per syndrome round at $n=144$: the Gross code uses **864**
two-qubit gates in **7** layers with 0 mixed checks; the PBB member with the
smallest combinatorial depth bound (`12_6_0193`) uses **1008** gates in **8**
layers with **72** mixed checks.  Both use 144 data + 144 ancilla qubits.
("Smallest depth bound" is not "best circuit": schedule choice alone moves the
logical error rate by 2.8x within a code, and we did not optimise either.)

> **Proposition C5 (translation-invariant schedules are strictly weaker).**
> **[exact computation]** Over the 318 catalogue PBB codes with $n\le180$, 235
> admit a translation-invariant schedule (depths realised
> $\{8:122,\ 10:101,\ 12:4,\ 14:6,\ 16:2\}$) and 83 do not, with CP-SAT
> returning `INFEASIBLE` for every depth up to (number of monomial directions)
> $+2$.  However, on the three smallest such codes the **unrestricted** model
> finds valid schedules — `phase2_26` at depth 7, `phase2_21` and `phase2_44`
> at depth 9 — each verified with 0 detector firings in 4000 noiseless shots.

**Correction notice.**  An earlier draft of Theorem C4 asserted that those 83
codes were "unschedulable".  That claim was **falsified by our own
experiment** (`exp009`): the obstruction belongs to the translation-invariant
*schedule class*, not to the codes.  The statement above is the corrected,
narrower one, and the depth separation for $[[144,12,12]]$ has been re-derived
from Proposition C2 so that it no longer depends on any schedule-class
restriction.  See `notes/failed_routes.md` FR-005.

> **Proposition C7 (the tested translation-invariant penalty is not a code
> penalty).** **[exact computation]** EXP-031 found three weight-9
> $[[144,12,12]]$ PBB members (`12_6_0188`, `12_6_0190`, catalogue line 41)
> for which the translation-invariant model is `INFEASIBLE` at every depth
> $9,\ldots,13$, refuting the proposed two-value law *in that restricted
> schedule class*.  EXP-033 removes the restriction: the unrestricted CP-SAT
> model returns an `OPTIMAL` schedule at depth **9** for all three, each with
> a persisted slot-map witness that is machine-checked offline — structurally
> (edge cover, ancilla/qubit slot-disjointness, parity criterion) and by Stim
> (0 detector firings and 0 observable flips in 4000 noiseless 12-round shots;
> `results/processed/exp033_noiseless_check.json`)
> (wall times: first run 466/117/1093 s at a 1500 s cap; witness-persisting
> rerun re-certified all three `OPTIMAL` at depth 9).  Since check weight
> lower-bounds depth by 9, these minima are exact.  Thus the apparent
> $>4$-layer penalty was caused entirely by imposing translation invariance:
> for these three instances
> $$T_{\rm unrestricted}=9,\qquad T_{\rm TI}>13.$$
> EXP-031 supports the two-value law on 150 unrestricted random commuting
> sets; EXP-033 supports it on these three adversarial PBB instances.  These
> finite data do **not** prove the law universally.

The EXP-031 exact-CSP predictor was exact on all 208 previously decided cases
(150 random commuting sets in the unrestricted class; 47 BB and 11 PBB
instances in the translation-invariant class).  Its weak closed-form proxy
(“some anticommuting overlap exists”) is not useful in general: over all 211
attempted instances the confusion matrix is 51 true positives, 157 false
positives, 3 true negatives and no false negatives.  The exact parity CSP, not
that proxy, is the candidate characterization still to prove or refute in the
unrestricted class.

---

## 5. Open sub-questions left by these results

1. Is there a PBB code with $|C|+|D|$ small enough (after $A\cap C$, $B\cap D$
   overlaps produce $Y$s) that max check weight reaches $7$, matching the BB
   depth?  Two catalogue codes have PBB weight 7; `phase2_26` was scheduled at
   depth 7.  Even then the check count and two-qubit gate count exceed the
   parent's, and Corollary 1 still caps the distance.
2. **ANSWERED, yes.**  `phase2_58` and `phase2_60` ($[[72,4,6]]$, $\delta=4$)
   have $d_X(Q')=4<6=d(Q)=d_Z(Q')$, so they beat their own CSS shadow.  Theorem
   2's ceiling still holds; Corollary 1 does not apply ($\delta\ne0$).  Both are
   nevertheless dominated by the certified Bravyi $[[72,12,6]]$.  **Still open:**
   whether *every* $\delta>0$ PBB code is dominated by *some* CSS code.  The
   strict-provenance audit (EXP-027, all 155 $\delta>0$ rows) proves 66
   dominations, certifies 2 reversals (`phase2_58`, `phase2_60` — two of the
   legacy EXP-012 "parent weaker" five), and leaves 87 undecided at budget;
   EXP-036 decides 46 of those (41 dominations, 5 new reversals; totals
   107/7/41); the remaining legacy tallies (49/5/23 on $n\le108$) used
   uncertified labels and are superseded.
3. Does the even-depth pattern among *translation-invariant* schedules
   ($8,10,12,14,16$ only) follow from a parity invariant of Lemma C1 on
   group-structured codes?
4. Do flag qubits or multi-ancilla (cat-state) extraction recover the lost
   depth for mixed checks, and at what ancilla cost?  **Still open, but one
   natural route to closing it is now ruled out.**

   *Failed proof route (recorded so it is not retried).*  It is tempting to
   make Theorem C4 unconditional by counting incidences: every check-qubit
   incidence $(g,j)$ *seems* to need its own two-qubit gate touching $j$, and
   gates touching $j$ must occupy distinct layers, which would give
   $T \ge \max_j \deg(j)$ for *any* ancilla scheme.

   **Two things are wrong with this, and neither is that the inequality is
   false.**

   *(a) The arithmetic does not reach 8 family-wide.*  Incidence counts differ
   across the 14 members.  Ten of them have $1080$ or $1152$ incidences over
   $144$ qubits (average degree $7.5$ or $8$), giving $T\ge8$; but **four** —
   including `12_6_0193`, the member used in every circuit benchmark here — have
   $\sum_g|\mathrm{supp}(g)| = 1008$, an average degree of **exactly 7**, so
   counting yields only $T\ge7$ and does not separate them from the Gross code.
   A family-wide counting proof is therefore impossible.  Theorem C4's $T\ge8$
   instead comes from the *maximum check weight* $\ge8$, which does hold for all
   14 (**[exact]** min 8, max 10) — but only for the published generators.  Theorem C4's $T\ge 8$ comes from
   the *maximum check weight* $\ge 8$ via Proposition C2, a different and valid
   argument that is restricted to one-ancilla schedules.

   *(b) The premise fails for multi-ancilla schemes.*  Ancilla-ancilla CNOTs
   forward parity, so one gate at a shared data qubit can serve two checks.
   For $g_1 = Z_0Z_1$, $g_2 = Z_0Z_2$ (so $\deg(0)=2$):

   ```
   R a1 a2 ; CX q0->a1 ; CX a1->a2 ; CX q1->a1 ; CX q2->a2 ; M a1 a2
   ```

   measures $g_1$ on `a1` and $g_2$ on `a2` (verified in
   `tests/test_conventions.py`) using **one** gate at $q_0$, not $\deg(0)=2$,
   and three data-touching gates rather than $|g_1|+|g_2| = 4$.

   Note carefully what this circuit does and does not show.  It uses **three**
   two-qubit layers while $\max_j\deg(j)=2$, so it *satisfies*
   $T\ge\max_j\deg(j)$ — it is **not** a counterexample to the inequality.
   What it refutes is the **premise** "one direct data gate per check
   incidence", and therefore the incidence-counting *derivation*.  Whether
   $T\ge\max_j\deg(j)$ holds for arbitrary multi-ancilla schemes is
   **open**; we have neither a proof nor a counterexample.  Theorem C4 must
   accordingly stay scoped to one-ancilla schedules.  What the trick costs is exactly what this programme measures:
   the forwarding edge makes a single ancilla fault corrupt *both* checks, so
   it buys gate count by manufacturing precisely the correlated hook errors
   that already make the PBB circuit worse.  Quantifying that trade for PBB
   requires building a multi-ancilla compiler and is genuinely not attempted
   here.

5. **Target equivalence (EXP-034) — now fully decided.**  For the benchmark
   member `12_6_0193`: the stabilizer row space satisfies
   $\dim S_X + \dim S_Z = 20 + 66 = 86 < 132 = \dim S$, so the code is not CSS
   under arbitrary row operations; these three ranks are invariant under every
   qubit permutation; an exact affine GF(2) feasibility system proves no
   local-Hadamard assignment makes it CSS (contradiction $0=1$ persisted); and
   the grouped X/Z column matroid is one circuit-connected component, so the
   row space is **indecomposable** under row operations, qubit permutations,
   and arbitrary independent local Cliffords (each local Clifford acts
   invertibly inside one qubit's two-column block and preserves every
   partition rank).

   > **Theorem E (genuine non-CSS-ness of the target). [proved + exact
   > computation]**  No transformation in the group generated by arbitrary
   > stabilizer row operations, arbitrary qubit permutations, and arbitrary
   > independent single-qubit Clifford transformations maps the stabilizer
   > group of `12_6_0193` to a CSS stabilizer group.
   >
   > *Proof.*  CSS-ness after a local Clifford assignment $Q$ is equivalent
   > to invariance of $S$ under the per-qubit conjugated rank-one projection
   > $Q P_X Q^{-1}$, giving, over the 864 one-hot selectors $y[j,g]$, the
   > persisted homogeneous parity masks; exactly-one per qubit linearly
   > implies $\bigoplus_g y[j,g] = 1$.  Nineteen parity masks plus the
   > qubit-0 one-hot row XOR to $0=1$ over GF(2), so the system has no
   > solution of any kind — in particular no one-hot solution.  Permutations
   > reduce to this case because CSS-ness is permutation-invariant and a
   > permutation conjugate of a local Clifford is again a local Clifford.
   > The 20-row certificate is persisted and machine-rechecked in
   > `results/processed/exp034_target_equivalence.json` and
   > `results/certificates/exp034_target_full_lc_decision.json`; OR-Tools
   > CP-SAT independently returns INFEASIBLE on the full one-hot model.
   > The $[[5,1,3]]$ perfect code is a non-vacuousness control: it is
   > LC-inequivalent to CSS (exhaustive over $6^5$) yet its linear
   > relaxation is feasible, so the linear route does not trivially refute
   > every code. $\square$

6. **Independent target distance (EXP-035) — now fully decided.**

   > **Theorem F (exact distance of the target). [proved + exact]**
   > $d(\texttt{12\_6\_0193}) = 12$.
   >
   > *Proof.*  Lower bound $d \ge 12$ by exhaustion over weight strata:
   > (i) weights $\le 5$: exact meet-in-the-middle exclusion of every
   > weight-$\le5$ centralizer element; (ii) weight 6: the disjoint
   > triple-triple syndrome MITM enumerates **exactly 72** weight-6
   > zero-syndrome vectors, each re-verified and shown to lie in $S$ (the
   > set equals the 72 stored pure-$Z$ checks), so no weight-6 logical
   > exists; (iii) weights 7–11: the machine-checked translation-orbit
   > reduction (72 verified lattice translations; four representatives
   > whose transported functionals span the full 24-dimensional quotient
   > dual, rank 24 on two independent GF(2) paths) reduces the stratum to
   > four weight-$\le11$ feasibility sectors, and all four are **UNSAT** —
   > proved by CaDiCaL (PySAT, equisatisfiable CNF: totalizer cardinality
   > + chained XOR) in 502/674/422/576 s under the frozen identity/trace
   > binding.  Upper bound $d \le 12$ by an independently re-verified
   > weight-12 witness (nontrivial on both GF(2) paths).  The canonical
   > certificate `results/certificates/pbb_12_6_0193_distance.json` binds
   > the exact catalogue row, matrix and basis hashes, all four sector
   > records, and the witness. $\square$

   The racing OR-Tools CP-SAT backend (7 workers, multi-hour budgets)
   produced no proof on any sector; the catalogue's `deep_milp` $d=12$ is
   hereby independently confirmed rather than trusted.
