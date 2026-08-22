# Finite-torus infrared/GKS LP floors evaporate as $L^{-1}$

Artifacts: `experiments/e93_mag_floor.py`, `results/bounds/mag_floor.json`, and
`tests/test_mag_floor.py`.

## 0. Result and scope

Let $T_L=(\mathbb Z/L\mathbb Z)^3$, $N=L^3$, and

\[
 \lambda(k)=3-\sum_{i=1}^3\cos k_i,
 \qquad
 C_L(z)=\frac1N\sum_{k\ne0}\frac{\cos(k\cdot z)}{\lambda(k)}.
\]

Write

\[
 D_L=C_L(0)-\min_z C_L(z),\qquad
 \kappa_L=\frac{D_L}{2},\qquad
 u_L=\frac{-\min_z C_L(z)}{D_L}.
\tag{1}
\]

As in `proofs/kc_interval2.md`, let $S_*(K,L)$ be the maximum of the exact
finite-torus Parseval/infrared/GKS/cap/energy linear program and set

\[
 F(K,L)=1-S_*(K,L).
\tag{2}
\]

Thus the audited finite-volume argument gives the valid physical lower bound
$M_L^2(K)\ge F(K,L)$.

> **[THEOREM, exact LP obstruction]** For every finite even $L\ge4$ and
> $K>0$, $F(K,L)>0$.  Nevertheless, for every $0<K<I_3/2$,
> \[
> F(K,L)\le\frac{5168}{525L}
> \quad\text{for all sufficiently large even }L,
> \tag{3}
> \]
> with an explicit sufficient cutoff below.  At the endpoint $K=I_3/2$ the
> same bound holds for every even $L\ge4$.
>
> Consequently there are no $0<K\le I_3/2$, $m>0$, and $L_0$ for which the
> **LP floor** satisfies $F(K,L)\ge m$ for every even $L\ge L_0$.

This is a negative result about the **audited LP certificate route**.  It
**does not** claim that the true Ising quantity $M_L^2(K)$ tends to zero.  In
particular, it does not contradict the physical possibility of ordering below
$I_3/2$.  It says only that this finite-torus constraint system cannot certify
the missing uniform floor.

The exact finite computation also extends the stored rational sides from
$L=4,6$ to $L=8$ (over $\mathbb Q(\sqrt2)$ for the LP), and obtains exact
Green-function break thresholds through $L=12$.

## 1. The exact LP and two finite-volume facts

For a full (not yet symmetrised) nonzero-mode profile $x(k)=\widehat G(k)$,
put

\[
 S=\frac1N\sum_{k\ne0}x(k),\qquad M^2=1-S,
\qquad G_x(z)=M^2+\frac1N\sum_{k\ne0}x(k)\cos(k\cdot z).
\tag{4}
\]

The inherited finite-torus constraints are

\[
\begin{aligned}
0&\le x(k)\le\frac1{2K\lambda(k)} &&(k\ne0),\\
0&\le G_x(z)\le1 &&(z\ne0),\\
S&\le1,\\
G_x(e)&\ge1-\frac{1-N^{-1}}{6K}.
\end{aligned}
\tag{5}
\]

The last row follows already by summing the pointwise ceiling against
$\lambda(k)$, but is retained in the stored LP to match the audited system.
Cubic signed-permutation averaging maps every full-LP feasible profile to an
orbit-constant feasible profile with the same $S$.  Thus the orbit LP stored in
the JSON is an **exact reduction of (5)**.  The LP itself is, of course, a
relaxation of the set of actual Ising two-point functions.

### 1.1 Positivity at every finite side

> **[LEMMA]** $F(K,L)>0$ for every nonempty finite LP in (5).

*Proof.*  If $S=1$, then $M^2=0$, $G_x(0)=1$, and the absent zero Fourier mode
implies

\[
\sum_{z\in T_L}G_x(z)=0.
\]

Hence $\sum_{z\ne0}G_x(z)=-1$, contradicting $G_x(z)\ge0$ for all $z\ne0$.
The feasible polytope is compact because every $x(k)$ has a finite ceiling, so
its maximum $S_*$ is attained and is strictly below one. $\square$

Physical finite-torus correlations supply a feasible point for the inherited
system, so the nonemptiness premise is available at every $K>0$ in the audited
setting.

### 1.2 Exact all-ceiling break and the compensating Green profile

> **[THEOREM]** The all-ceiling profile
> \[
> x_{\rm ceil}(k)=\frac1{2K\lambda(k)}
> \tag{6}
> \]
> is feasible if and only if $2K\ge D_L$.  In that regime it is componentwise
> maximal, and therefore
> \[
> F(K,L)=1-\frac{C_L(0)}{2K}.
> \tag{7}
> \]
> If $2K\le D_L$, the profile
> \[
> x_D(k)=\frac1{D_L\lambda(k)}
> \tag{8}
> \]
> is feasible and gives
> \[
> F(K,L)\le u_L.
> \tag{9}

*Proof.*  With (6), equation (4) becomes

\[
G_{\rm ceil}(z)=1-\frac{C_L(0)-C_L(z)}{2K}.
\tag{10}
\]

The GKS lower rows are therefore equivalent to $2K\ge D_L$.  The upper rows
hold because $C_L(z)\le C_L(0)$ termwise, the cap holds because
$D_L>C_L(0)$, and the energy row follows from the ceilings.  Since every
objective coefficient is positive and (6) saturates every individual upper
bound, it is optimal.

For (8),

\[
G_D(z)=1-\frac{C_L(0)-C_L(z)}{D_L}\in[0,1].
\tag{11}
\]

The ceiling holds precisely when $2K\le D_L$; the cap and energy row again
follow.  Its value is $S=C_L(0)/D_L$, yielding (9). $\square$

Thus the exact finite-L point at which the ordinary all-ceiling constraint
chain breaks is not heuristic: it is exactly $K=\kappa_L$.  Below that value,
the most negative Green site makes (10) negative.  At the minimizer,

\[
G_{\rm ceil}(z_{\min})=1-\frac{D_L}{2K}<0\quad(K<\kappa_L).
\tag{12}
\]

## 2. A uniform exact $O(L^{-1})$ upper bound for the LP floor

The following three lemmas give the all-sizes part of the result.  No
floating-point comparison decides any displayed inequality.

**[EXTERNAL]** The only nonlocal numerical input in the all-sizes comparison is
the certified rational lower endpoint for $I_3$ recorded in
`results/bounds/upper_infrared.json`; it is used only in the indicated
one-sided comparisons.  The finite-torus identities and all LP decisions in
this note are recomputed with exact rational or quadratic-field arithmetic.

### 2.1 Controlling the negative torus-Green defect

Define

\[
J_2(L)=\frac1N\sum_{k\ne0}\lambda(k)^{-2}.
\tag{13}
\]

> **[LEMMA]** For every $\varepsilon>0$ and every $z$,
> \[
> C_L(z)\ge-\frac1{\varepsilon N}-\varepsilon J_2(L).
> \tag{14}
> \]
> In particular,
> \[
> -\min_zC_L(z)\le2\sqrt{\frac{J_2(L)}N}.
> \tag{15}

*Proof.*  Let

\[
H_{L,\varepsilon}(z)=\frac1N\sum_k
 \frac{\cos(k\cdot z)}{\lambda(k)+\varepsilon}.
\]

The massive geometric expansion

\[
\frac1{3+\varepsilon-\gamma(k)}
=\frac1{3+\varepsilon}\sum_{n\ge0}
 \left(\frac{\gamma(k)}{3+\varepsilon}\right)^n
\]

is uniformly convergent because $|\gamma|\le3<3+\varepsilon$.  Its real-space
coefficients are nonnegative torus walk weights, hence
$H_{L,\varepsilon}(z)\ge0$.  After removing the $k=0$ term and writing

\[
\frac1\lambda=\frac1{\lambda+\varepsilon}
 +\frac{\varepsilon}{\lambda(\lambda+\varepsilon)},
\]

the second term is bounded below by $-\varepsilon J_2(L)$, proving (14).
Optimising at $\varepsilon=(NJ_2(L))^{-1/2}$ gives (15). $\square$

### 2.2 An elementary lattice-sum bound

> **[LEMMA]** For every $L\ge2$,
> \[
> J_2(L)\le\frac{L\Sigma_4}{64},\qquad
> \Sigma_4=\sum_{j\in\mathbb Z^3\setminus\{0\}}\frac1{|j|^4}
> \le4\pi^2+\frac{\pi^4}{45}
> \le\left(\frac{323}{50}\right)^2.
> \tag{16}

*Proof.*  Choose each momentum representative with coordinates in
$(-L/2,L/2]$.  The elementary bound

\[
1-\cos t\ge\frac{2t^2}{\pi^2}\qquad(|t|\le\pi)
\]

follows because $(1-\cos t)/t^2$ decreases on $(0,\pi]$.  Hence
$\lambda(2\pi j/L)\ge8|j|^2/L^2$, which proves the first inequality in
(16).  Grouping the integer lattice by $\lVert j\rVert_\infty=m$ gives a
shell of $24m^2+2$ points, so

\[
\Sigma_4\le\sum_{m\ge1}\frac{24m^2+2}{m^4}
=24\zeta(2)+2\zeta(4)=4\pi^2+\frac{\pi^4}{45}.
\]

Finally use the exact classical rational bound $\pi<355/113$ and square
comparison with $323/50$. $\square$

### 2.3 A lower Riemann-sum estimate for $C_L(0)$

Let $C_\mu(0)=\langle(\lambda+\mu^2)^{-1}\rangle$ denote the infinite-lattice
massive Green value.

> **[LEMMA]** With
> \[
> c_\Delta=
> \frac{99}{70}\,\frac{(355/113)^2+2(355/113)}{16}
> =\frac{583407}{408608},
> \qquad A=c_\Delta+1=\frac{992015}{408608},
> \tag{17}
> \]
> one has
> \[
> C_L(0)\ge I_3-\frac{A}{L}.
> \tag{18}

*Proof.*  Set $\mu=L^{-1}$.  Positivity of the massive torus walk expansion
shows

\[
\begin{aligned}
C_L(0)
&\ge \frac1N\sum_{k\ne0}\frac1{\lambda(k)+\mu^2}\\
&=H_{L,\mu}(0)-\frac1{N\mu^2}
 \ge C_\mu(0)-\frac1L.
\end{aligned}
\tag{19}
\]

For the continuum difference,

\[
I_3-C_\mu(0)=\mu^2\left\langle
\frac1{\lambda(\lambda+\mu^2)}\right\rangle.
\tag{20}
\]

On $\{\lambda\ge\mu^2\}$, the bound $\lambda\le|k|^2/2$ forces
$|k|\ge\sqrt2\,\mu$, while $\lambda\ge2|k|^2/\pi^2$ bounds the integrand by
$\mu^2\pi^4/(4|k|^4)$.  Integrating over the containing exterior ball gives
$\pi^2\mu/(8\sqrt2)$.  On $\{\lambda<\mu^2\}$, the lower dispersion bound
forces $|k|<\pi\mu/\sqrt2$, and the integrand is at most
$\pi^2/(2|k|^2)$; its ball integral is $\pi\mu/(4\sqrt2)$.  Thus

\[
I_3-C_\mu(0)\le
\frac{\sqrt2(\pi^2+2\pi)}{16}\,\mu\le c_\Delta\mu.
\tag{21}
\]

Combining (19) and (21) proves (18). $\square$

### 2.4 Assembly and the endpoint

> **[THEOREM]** Let
> \[
> c_* = \frac{5168}{525}.
> \tag{22}
> \]
> For every $0<K<I_3/2$ and every even
> \[
> L\ge
> \max\left\{4,\left\lceil\frac{A}{I_3-2K}\right\rceil\right\},
> \tag{23}
> \]
> one has $F(K,L)\le c_*/L$.  At $K=I_3/2$, the same conclusion holds for
> every even $L\ge4$.

*Proof.*  Since $\lambda(k)\le6$,

\[
D_L\ge C_L(0)\ge\frac{N-1}{6N}.
\tag{24}
\]

Equations (15)–(16) give

\[
-\min C_L\le\frac{\sqrt{\Sigma_4}}{4L}.
\]

Therefore

\[
 u_L\le\frac{6N}{N-1}\frac{\sqrt{\Sigma_4}}{4L}
 \le\frac{32}{21}\frac{323}{50L}
 =\frac{5168}{525L}.
\tag{25}
\]

For (23), (18) yields $D_L\ge C_L(0)\ge2K$, so (9) and (25) apply.

At the endpoint take the single profile

\[
 x(k)=\frac{c}{\lambda(k)},\qquad c=\min\{I_3^{-1},D_L^{-1}\}.
\tag{26}
\]

It is feasible by the same proof as (8).  If $D_L\ge I_3$, it gives the
$u_L$ bound.  If $D_L\le I_3$, it gives

\[
F(I_3/2,L)\le1-\frac{C_L(0)}{I_3}
\le\frac{A}{I_3L}
\le\frac{c_*}{L},
\]

where the last comparison is an exact rational check using the certified lower
endpoint for $I_3$. $\square$

> **[FALSIFIED] LP-uniform-floor premise.** The quantified premise required by
> the route in `proofs/kc_interval2.md`, namely a positive $m$ with
> $F(K,L)\ge m$ at all large even $L$, fails for every $0<K\le I_3/2$.  This is a
> statement about $F$, not about the true $M_L^2$.

## 3. Exact finite certificates

**[COMPUTATION]** The exact Green data are as follows.  Every minimum occurs at
the antipodal site $(L/2,L/2,L/2)$; the standalone verifier checks the Green
defect equation and zero spatial sum at every site.

| $L$ | $C_L(0)$ | $\min C_L$ | $\kappa_L$ | $u_L$ | $L u_L$ |
|---:|---|---|---|---|---|
| 4 | $1517/3840$ | $-49/1280$ | $13/60$ | $147/1664$ | $147/416$ |
| 6 | $1289503/2993760$ | $-13861/598752$ | $5147/22680$ | $69305/1358808$ | $69305/226468$ |
| 8 | $118783817/264427520$ | $-13281829/793282560$ | $36097/154938$ | $13281829/369633280$ | $13281829/46204160$ |
| 10 | $852594421766635003/1851857825844180000$ | $-24381639636139997/1851857825844180000$ | $313487315/1323944776$ | $24381639636139997/876976061402775000$ | $24381639636139997/87697606140277500$ |
| 12 | $733533368777249/1567766752231680$ | $-39905507821/3671584899840$ | $422619944041/1765503099360$ | $17039651839567/750573020616816$ | $17039651839567/62547751718068$ |

The exact inequalities

\[
\frac{13}{60}<\frac{5147}{22680}<\frac{36097}{154938}
<\frac{313487315}{1323944776}
<\frac{422619944041}{1765503099360}<\frac{I_3}{2}
\tag{27}
\]

are checked with the certified lower endpoint of $I_3$.  They quantify the
finite-L all-ceiling failure requested in the target: its threshold rises
through all five checked sides toward the incumbent infrared endpoint.

**[COMPUTATION]** The $L=8$ LP is solved over the ordered real field
$\mathbb Q(\sqrt2)$; a value `Q[a|b|2]` in the JSON means $a+b\sqrt2$.  Exact
primal and dual certificates agree at all stored couplings:

| $K$ | $F(K,4)$ | $F(K,6)$ | $F(K,8)$ |
|---:|---|---|---|
| $1/5$ | $251/4608$ | $253/15552$ | $1813/208896-\sqrt2/816$ |
| $21/100$ | $607/8064$ | $133/5832$ | $13171/1354752$ |
| $11/50$ | $863/8448$ | $34663/940896$ | $5829/382976$ |
| $6/25$ | $1631/9216$ | $737509/7185024$ | $40706963/634626048$ |

For each displayed $K$, the three floors strictly decrease from $L=4$ to
$L=6$ to $L=8$, an exact finite computation.  The $L=4,6$ rows reproduce all
eight stored certificates from `kc_interval2.json`; the verifier rebuilds all
12 (including $L=8$) matrices, primal vectors, dual vectors, objectives, and
bases without importing this producer.

> **[FALSIFIED] Candidate dyadic LP-floor monotonicity.** The natural
> nondecreasing relation $F(K,2L)\ge F(K,L)$ already fails exactly at
> $(K,L)=(1/5,4)$:
> \[
> F(1/5,4)-F(1/5,8)
> =\frac{28697}{626688}+\frac{\sqrt2}{816}>0.
> \tag{27a}
> \]
> This is a certificate against that LP-floor monotonicity, not a statement
> about a possible reflection-positivity relation for the true Ising
> magnetization with different hypotheses or error terms.

## 4. Modified constraint systems

> **[THEOREM, limited meta-negative]** Any construction obtained solely by
> dropping a row, averaging a pointwise constraint over momenta/orbits, or
> replacing a constraint by a convex consequence is a relaxation of (5).  Its
> feasible set is larger, so its floor is at most $F(K,L)$ and therefore also
> has the $O(L^{-1})$ obstruction above.

This rules out the suggested momentum-averaging and paired-momentum operations
when they are merely repackagings/weakenings of the existing ceiling and GKS
rows.  A genuinely stronger paired-momentum inequality could only help if it
is a **new valid physical input** not implied by the audited system (for
example, a random-current, Simon--Lieb, GKS-II four-point, or ABF-type input).

> **[UNRESOLVED]** No such new stronger constraint is proved here.  This note
> does not rule out an upper-endpoint refinement by a method outside the audited
> two-point finite-torus LP class.

## 5. A sharp finitely-checkable conditional and its countercertificate

For fixed rational $K,m>0$, define the exact predicate

\[
H_{K,m}(L):\quad F(K,L)\ge m.
\tag{28}
\]

Each instance is decidable by the stored Fraction/quadratic-field primal-dual
certificate.

> **[CONDITIONAL]** If $H_{K,m}(L)$ holds for every sufficiently large even
> $L$, then the finite-torus infrared argument supplies a uniform physical
> magnetization floor $M_L^2(K)\ge m$, and the standard thermodynamic-limit
> implication used in `kc_interval2.md` yields $K_c\le K$.

The pre-fixed e85 grid point $K=11/50$ and $m=1/100$ give a concrete test of
this conditional.  Exact checks are

\[
\frac{863}{8448}>\frac1{100},\qquad
\frac{34663}{940896}>\frac1{100},\qquad
\frac{5829}{382976}>\frac1{100},
\tag{29}
\]

so $H_{11/50,1/100}(4)$, $H_{11/50,1/100}(6)$, and
$H_{11/50,1/100}(8)$ all hold.

> **[FALSIFIED]** The all-L premise of this particular conditional is false:
> the exact rate certificate gives $2K\le I_3-A/L$ for every $L\ge38$, and
> $c_*/L<1/100$ for every $L\ge985$.  Therefore, for every even
> \[
> L\ge986,
> \qquad F(11/50,L)\le\frac{c_*}{L}<\frac1{100},
> \tag{30}
> \]
> so $H_{11/50,1/100}(L)$ fails.  The `986` cutoff is an exact rational
> comparison, not a float locator.

Thus the requested $L=4,6,8$ conditional checks are available and their
failure at all sufficiently large sides is now certified, rather than merely
unresolved.

## 6. Reproduction and classification

**[COMPUTATION]** From the repository root, run

```sh
timeout 3600 .venv/bin/python experiments/e93_mag_floor.py
timeout 3600 .venv/bin/python tests/test_mag_floor.py
```

The producer took about `67` seconds in the recorded artifact; the clean-room
verifier took about `11` seconds in this run.  The 3600-second timeout
exceeds twice either wall time.

| item | status |
|---|---|
| finite-L strict positivity of the audited LP floor | [THEOREM] |
| exact all-ceiling threshold and low-K Green-profile certificate | [THEOREM] |
| all-even-L $O(L^{-1}) LP-floor collapse | [THEOREM] |
| uniform LP floor needed by the infrared upper refinement | [FALSIFIED] |
| $L=8$ quadratic-field LP certificates; $L=10,12$ Green data | [COMPUTATION] |
| the $H_{11/50,1/100}$ finite checks | [CONDITIONAL], then [FALSIFIED] globally |
| strengthened constraints outside the audited two-point system | [UNRESOLVED] |
| true physical $M_L^2$ uniformity below $I_3/2$ | [UNRESOLVED; not claimed false] |
