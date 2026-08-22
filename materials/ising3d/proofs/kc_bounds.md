# Certified elementary bounds for the simple-cubic Ising critical coupling

## Result and conventions

For the nearest-neighbor ferromagnetic Ising model

\[
  Z_\Lambda(K)=\sum_{\sigma\in\{-1,1\}^{\Lambda}}
  \exp\!\left(K\sum_{\langle x,y\rangle}\sigma_x\sigma_y\right),
  \qquad v=\tanh K,
\]

the computation in `experiments/e05_kc_bounds.py` gives the **certified** interval

\[
\boxed{
  0.2074277114992039908436804465100887760027
  \;\le K_c^{\rm sc}\le\;
  0.2527310098586630030260020266135701299926 .
}
\]

The displayed endpoints are rounded outward.  The exact constants behind them are

\[
 \operatorname{atanh}\!\left(4468911678^{-1/14}\right)
 \le K_c^{\rm sc}
 \le \frac{W_{\rm sc}}6,
\]

where

\[
 W_{\rm sc}=\frac{\sqrt6}{32\pi^3}
 \Gamma\!\left(\frac1{24}\right)
 \Gamma\!\left(\frac5{24}\right)
 \Gamma\!\left(\frac7{24}\right)
 \Gamma\!\left(\frac{11}{24}\right).
\]

The benchmark `0.221654626(5)` is inside the interval.  It was used only as a
falsification check and nowhere in deriving either endpoint.

## Dependency ledger

Every constant in the final interval has the following dependency chain.

| Constant | Logical dependencies | Status |
|---|---|---|
| `4468911678 = c_14` | exhaustive visited-set backtracking on `Z^3`; cubic-symmetry orbit weights; independent comparison with [OEIS A001412](https://oeis.org/A001412/b001412.txt) | exact, computer-assisted |
| `4468911678^(1/14)` | `c_(m+n) <= c_m c_n` and Fekete's lemma | rigorous theorem plus interval evaluation |
| lower endpoint | high-temperature SAW domination, then monotonicity of `atanh` | rigorous theorem plus exact count and interval evaluation |
| `log(1+sqrt(2))/2` | exact square-lattice critical point; GKS-II monotonicity when vertical bonds are added | rigorous theorem |
| factor `1/2` in `I_3/2` | FSS Gaussian domination with the physical bond normalization derived below; Parseval spin sum rule | rigorous theorem |
| `W_sc/3 = I_3` | exact change between Watson-integral conventions | algebraic identity |
| gamma product with denominator `32*pi^3` | Glasser--Zucker exact evaluation | rigorous published identity |
| decimal endpoints | 80-decimal-digit directed-rounding interval arithmetic, then outward decimal rounding | computer-assisted rigorous |

The `d=4,5` values and the independent Bessel quadrature are checks only and do
not enter the certified interval.

---

## A. Lower bound from self-avoiding-walk domination

### A.1 Precise correlation theorem

**Finite-graph SAW domination.**  Let `G=(V,E)` be a finite graph, let every
edge coupling satisfy `K_e >= 0`, put `v_e=tanh(K_e)`, and impose zero external
field and free boundary conditions.  For distinct vertices `a,b`,

\[
  \langle\sigma_a\sigma_b\rangle_G
  \le
  \sum_{\gamma:a\to b\ {
m self\!\!-avoiding}}
  \prod_{e\in\gamma}v_e.                                      \tag{A.1}
\]

The same inequality holds in an infinite locally finite graph after taking a
finite-volume exhaustion.  This is the path form of the high-temperature
Simon--Lieb type correlation bound.  The more general separator inequalities
are due to B. Simon, *Correlation inequalities and the decay of correlations
in ferromagnets*, Commun. Math. Phys. **77** (1980), 111--126,
[doi:10.1007/BF01982711](https://doi.org/10.1007/BF01982711), and E. H. Lieb,
*A refinement of Simon's correlation inequality*, Commun. Math. Phys. **77**
(1980), 127--135,
[doi:10.1007/BF01982712](https://doi.org/10.1007/BF01982712).  Equation (A.1)
is proved directly here, so no unquoted normalization from either paper is
being assumed.

**Proof of (A.1).**  The high-temperature expansion gives

\[
 \langle\sigma_a\sigma_b\rangle_G=
 \frac{\displaystyle\sum_{A\subseteq E:\,\partial A=\{a,b\}}
                    \prod_{e\in A}v_e}
      {\displaystyle\sum_{B\subseteq E:\,\partial B=\varnothing}
                    \prod_{e\in B}v_e},                         \tag{A.2}
\]

where `partial A` is the set of odd-degree vertices of the subgraph `A`.
Every `A` in the numerator contains an `a`-to-`b` trail; erase loops from a
fixed deterministic choice of such a trail to obtain a self-avoiding path
`gamma(A) subset A`.  Then

\[
 B=A\setminus\gamma(A)=A\mathbin\triangle\gamma(A)
\]

has even degree at every vertex, and
`prod_(e in A) v_e = (prod_(e in gamma) v_e)(prod_(e in B) v_e)`.
The map `A -> (gamma(A),B)` is injective because `A=gamma union B`.  Dropping
the restrictions that `B` be disjoint from `gamma` and that `gamma` be the
deterministically selected path can only increase the positive sum.  The
numerator of (A.2) is therefore at most the denominator times the right side
of (A.1).  This proves the claim.  Positivity and monotone finite-volume
limits give the infinite-volume version.  QED.

For uniform coupling on `Z^3`, sum (A.1) over the endpoint `x`.  Each rooted
`n`-step SAW has one endpoint, so

\[
 \chi(K)=\sum_x\langle\sigma_0\sigma_x\rangle
 \le \sum_{n\ge0}c_n v^n.                                      \tag{A.3}
\]
The same inequality directly rules out spontaneous magnetization.  Represent
plus boundary conditions in a box by a fixed ghost spin joined to every
boundary site with coupling `+infinity`, hence ghost-edge weight
`tanh(+infinity)=1`.  Applying (A.1) from the origin to the ghost bounds the
finite-volume plus magnetization by the tail
`sum_(n >= dist(0,boundary)) c_n v^n`.  If `v < 1/mu`, the SAW generating
function converges by the Cauchy--Hadamard formula, and this tail tends to zero
as the box grows.  Thus the infinite-volume plus magnetization vanishes.
Equivalently, (A.3) also shows finite susceptibility.  Therefore a point with
spontaneous order (or divergent susceptibility) must have `v >= 1/mu`, where
`mu` is the connective constant derived next.  Thus

\[
 v_c\ge\frac1\mu.                                                \tag{A.4}
\]

### A.2 Why a finite count gives an upper, not a lower, bound on `mu`

Let `c_n` be the number of rooted `n`-step SAWs on `Z^3`.  Splitting an
`(m+n)`-step walk after step `m`, translating the suffix to start at the
origin, and retaining the prefix defines an injection into a pair consisting
of an `m`-step SAW and an `n`-step SAW.  Not every such pair concatenates to a
SAW, hence

\[
 c_{m+n}\le c_m c_n.                                             \tag{A.5}
\]

Fekete's lemma applied to `log c_n` yields

\[
 \mu=\lim_{n\to\infty}c_n^{1/n}=\inf_{n\ge1}c_n^{1/n}.           \tag{A.6}
\]

This classical connective-constant construction goes back to J. M.
Hammersley, *Percolation processes II. The connective constant*, Proc.
Cambridge Philos. Soc. **53** (1957), 642--645.  The direction needed here is
therefore

\[
 \mu\le c_n^{1/n}\quad\text{for every }n,                        \tag{A.7}
\]

not the reverse.

Combining (A.4) and (A.7) gives, for every exactly counted `c_n`,

\[
 K_c=\operatorname{atanh}(v_c)
 \ge \operatorname{atanh}(c_n^{-1/n}).                           \tag{A.8}
\]

### A.3 Independent exact enumeration certificate

The program performs ordinary backtracking.  The state is the current lattice
vertex and an exact set of visited integer-coordinate vertices.  Every one of
the six nearest-neighbor moves is tried, and a branch is rejected exactly when
its next vertex is already visited.  No tabulated count is used.

For speed only, all length-six prefixes are grouped under the 48 signed
coordinate permutations of the cubic lattice.  The program first enumerates
every actual prefix and records how many map to each canonical representative;
that enumerated number, not an assumed value of 48, is the orbit weight.  It
then backtracks every representative to length 14 and multiplies by the exact
orbit weight.  Different prefix orbits are processed in independent worker
processes.  Python arbitrary-precision integers are used throughout.

| `n` | exact `c_n` |
|---:|---:|
| 0 | 1 |
| 1 | 6 |
| 2 | 30 |
| 3 | 150 |
| 4 | 726 |
| 5 | 3534 |
| 6 | 16926 |
| 7 | 81390 |
| 8 | 387966 |
| 9 | 1853886 |
| 10 | 8809878 |
| 11 | 41934150 |
| 12 | 198842742 |
| 13 | 943974510 |
| 14 | **4468911678** |

All fifteen values were compared at run time with the live OEIS b-file for
[A001412](https://oeis.org/A001412/b001412.txt); all agree.  The finite set of
instances of (A.5) with `m+n <= 14` was also checked exactly.  The OEIS check is
an independent bug detector, not an input to the enumeration.

At 80 decimal digits, directed-rounding interval arithmetic gives

\[
\begin{aligned}
  4.889901675184127756927959967629262085367217791730769405727852339786059179862818607633
  &\le c_{14}^{1/14}\\
  &\le
  4.889901675184127756927959967629262085367217791730769405727852339786059179862818632934,
\end{aligned}
\]

and

\[
\begin{aligned}
0.2074277114992039908436804465100887760027696059720007957779765967418812997681495035536
&\le \operatorname{atanh}(c_{14}^{-1/14})\\
&\le
0.2074277114992039908436804465100887760027696059720007957779765967418812997681495067163.
\end{aligned}
\]

Taking the lower endpoint and rounding downward produces the certified lower
endpoint displayed at the start.

---

## B1. Upper bound from square-lattice layers

The GKS-II inequality for a finite ferromagnetic Ising model says

\[
 \langle\sigma_A\sigma_B\rangle
 -\langle\sigma_A\rangle\langle\sigma_B\rangle\ge0.              \tag{B.1}
\]

Therefore differentiation with respect to any ferromagnetic bond coupling
`K_e` gives

\[
 \frac{\partial}{\partial K_e}\langle\sigma_A\rangle
 =\langle\sigma_A\sigma_e\rangle
  -\langle\sigma_A\rangle\langle\sigma_e\rangle\ge0.             \tag{B.2}
\]

The same conclusion holds with plus boundary conditions (equivalently, a
fixed plus ghost spin), and survives the thermodynamic limit.  Thus adding
ferromagnetic bonds cannot decrease correlations or plus-state
magnetization.  Sources are R. B. Griffiths, *Correlations in Ising
Ferromagnets. I*, J. Math. Phys. **8** (1967), 478--483,
[doi:10.1063/1.1705219](https://doi.org/10.1063/1.1705219), and D. G. Kelly and
S. Sherman, *General Griffiths' Inequalities on Correlations in Ising
Ferromagnets*, J. Math. Phys. **9** (1968), 466--484,
[doi:10.1063/1.1664600](https://doi.org/10.1063/1.1664600).

Set every vertical simple-cubic bond to zero.  The model is then a disjoint
stack of square-lattice Ising models.  For any
`K > K_c^square`, a layer has positive plus-state magnetization.  Turning the
vertical couplings from zero up to `K` cannot destroy it by (B.2).  Taking the
infimum over ordered `K` gives

\[
 K_c^{\rm sc}\le K_c^{\rm square}
 =\frac12\log(1+\sqrt2)
 =0.440686793509771512616304662489896\ldots .                    \tag{B.3}
\]

The exact square-lattice value is from L. Onsager, *Crystal Statistics. I. A
Two-Dimensional Model with an Order-Disorder Transition*, Phys. Rev. **65**
(1944), 117--149,
[doi:10.1103/PhysRev.65.117](https://doi.org/10.1103/PhysRev.65.117); the
spontaneous-magnetization statement is C. N. Yang, Phys. Rev. **85** (1952),
808--816, [doi:10.1103/PhysRev.85.808](https://doi.org/10.1103/PhysRev.85.808).
This route is rigorous but weaker than B2.

---

## B2. Reflection positivity and the infrared bound

### B2.1 The theorem with its normalization fixed

The Gaussian-domination/infrared theorem of J. Fröhlich, B. Simon, and T.
Spencer applies to reflection-positive ferromagnetic interactions on even
tori.  In the convention

\[
 H=-\sum_{x,y}j(y-x)\sigma_x\sigma_y,
 \qquad j(z)=j(-z)\ge0,
\]

define

\[
 E(k)=\sum_z j(z)\,[1-\cos(k\cdot z)].                            \tag{B.4}
\]

For every nonzero torus momentum,

\[
 \widehat G_L(k)\le\frac1{2\beta E(k)},
 \qquad
 \widehat G_L(k)=\sum_x e^{ik\cdot x}\langle\sigma_0\sigma_x\rangle_L.
                                                                    \tag{B.5}
\]

The original source is J. Fröhlich, B. Simon, and T. Spencer, *Infrared
bounds, phase transitions and continuous symmetry breaking*, Commun. Math.
Phys. **50** (1976), 79--95,
[doi:10.1007/BF01608557](https://doi.org/10.1007/BF01608557).  An explicit
Ising statement with `E(k)` is Proposition 1.3, equations (1.15)--(1.17), of
M. Aizenman, H. Duminil-Copin, and V. Sidoravicius, *Random Currents and
Continuity of Ising Model's Spontaneous Magnetization*, Commun. Math. Phys.
**334** (2015), 719--742,
[doi:10.1007/s00220-014-2093-y](https://doi.org/10.1007/s00220-014-2093-y).
The nearest-neighbor interaction is reflection positive, satisfying the
hypotheses of the theorem.

Here the ordered sum counts each physical undirected bond twice.  To reproduce
our Hamiltonian `-J sum_<x,y> sigma_x sigma_y`, one must set

\[
 j(+e_i)=j(-e_i)=\frac J2.                                       \tag{B.6}
\]

Consequently

\[
 E(k)=J\,[3-\cos k_1-\cos k_2-\cos k_3],                         \tag{B.7}
\]

and, because `K=beta J`, (B.5) becomes

\[
 \widehat G_L(k)\le
 \frac1{2K(3-\cos k_1-\cos k_2-\cos k_3)}.                       \tag{B.8}
\]

This line is the source of the crucial factor two.  Putting `j(+/-e_i)=J`
inside the ordered-pair Hamiltonian would make each physical bond have
strength `2J`; calling `beta J` the physical `K` after that substitution is a
normalization error.  It would yield the false candidate `I_3/4 =
0.1263655049293315...`, below the falsification target even though the real
model is disordered there.  Equation (B.6) removes that double counting.

### B2.2 Sum rule and phase-transition bound

Parseval and `sigma_0^2=1` give on a torus of `N` sites

\[
 1=\frac1N\sum_k\widehat G_L(k)
  =\underbrace{\frac{\widehat G_L(0)}N}_{M_L^2}
   +\frac1N\sum_{k\ne0}\widehat G_L(k).                           \tag{B.9}
\]

Apply (B.8), then take even-torus thermodynamic limits.  In three dimensions
the singularity is integrable, so

\[
 \liminf_{L\to\infty}M_L^2
 \ge 1-\frac{I_3}{2K},                                           \tag{B.10}
\]

where

\[
 I_3=\frac1{(2\pi)^3}\int_{[-\pi,\pi]^3}
 \frac{d^3k}{3-\cos k_1-\cos k_2-\cos k_3}.                      \tag{B.11}
\]

If `K>I_3/2`, (B.10) proves long-range order.  The usual FSS thermodynamic
limit then gives nonzero spontaneous magnetization.  Therefore

\[
 K_c^{\rm sc}\le\frac{I_3}{2}.                                  \tag{B.12}
\]

### B2.3 Independent numerical integral and the exact Watson identity

For a numerical cross-check, the script uses

\[
 \frac1a=\int_0^\infty e^{-at}\,dt,
 \qquad
 \frac1{2\pi}\int_{-\pi}^{\pi}e^{t\cos k}\,dk=I_0(t),
\]

to obtain

\[
 I_d=\int_0^\infty e^{-dt}I_0(t)^d\,dt.                           \tag{B.13}
\]

It evaluates (B.13) after `t=(s/(1-s))^2`; the identity
`e^(-t) I_0(t) = 1F1(1/2;1;-2t)` prevents overflow.  At 80 decimal digits the
numerical quadrature gives

\[
 I_3=0.5054620197173260060520040532271402599851290148174208921889934878860288\ldots .
\]

This quadrature is **numerical only**.  The certificate instead uses the exact
closed form.  With the standard Watson convention,

\[
\begin{aligned}
 W_{\rm sc}
 &=\frac1{\pi^3}\int_{[0,\pi]^3}
   \frac{d^3k}{1-(\cos k_1+\cos k_2+\cos k_3)/3}\\
 &=3I_3.
\end{aligned}                                                     \tag{B.14}
\]

M. L. Glasser and I. J. Zucker, *Extended Watson integrals for the cubic
lattices*, Proc. Natl. Acad. Sci. USA **74** (1977), 1800--1801,
[doi:10.1073/pnas.74.5.1800](https://doi.org/10.1073/pnas.74.5.1800), give

\[
 W_{\rm sc}=\frac{\sqrt6}{32\pi^3}
 \Gamma\!\left(\frac1{24}\right)
 \Gamma\!\left(\frac5{24}\right)
 \Gamma\!\left(\frac7{24}\right)
 \Gamma\!\left(\frac{11}{24}\right)
 =1.516386059151978018156012159681420779955\ldots .               \tag{B.15}
\]

The expression in the task used `1/(4*pi^3)`, not `1/(32*pi^3)`.  Direct
high-precision evaluation gives

\[
 \frac{\text{task-stated expression}}{W_{\rm sc}}=8
\]

to more than 70 digits; algebraically the ratio is already exactly eight.
Thus the supplied formula does **not** match the integral.  It is not forced
to agree: the denominator `32` in the cited identity is used, and the
factor-eight discrepancy is recorded as a passing detection check in the
JSON artifact.

Evaluating (B.15)/3 with `mpmath.iv` at `iv.dps=80` gives the outward interval

\[
\begin{aligned}
0.5054620197173260060520040532271402599851290148174208921889934878860287734511738095325
\le I_3\le{}\\
0.5054620197173260060520040532271402599851290148174208921889934878860287734511738232374.
\end{aligned}
\]

Therefore

\[
\begin{aligned}
0.2527310098586630030260020266135701299925645074087104460944967439430143867255869042392
\le \frac{I_3}{2}\le{}\\
0.2527310098586630030260020266135701299925645074087104460944967439430143867255869121458.
\end{aligned}
\]

The upper endpoint, rounded upward to 40 places, is the certified final upper
endpoint.

### B2.4 Dimension and normalization cross-checks

Near `k=0`, `d-sum_i cos(k_i) = |k|^2/2 + O(|k|^4)`, so the integral behaves
locally like `integral r^(d-3) dr`.  Hence `I_2` diverges logarithmically.  The
infrared criterion is therefore vacuous in two dimensions, as it must be; it
does not incorrectly prove 2D long-range order from (B.10).

The same independent quadrature gives the following **numerical-only** trend:

| `d` | `I_d` | `I_d/2` |
|---:|---:|---:|
| 2 | diverges | vacuous |
| 3 | 0.505462019717326006052004053227140260 | 0.252731009858663003026002026613570130 |
| 4 | 0.309866780462120428169674416214750178 | 0.154933390231060214084837208107375089 |
| 5 | 0.231261624968046235741427024387713397 | 0.115630812484023117870713512193856699 |

The monotone decrease is consistent with increasing coordination.  These
`d=4,5` point values are not interval certificates and are not used in the
three-dimensional conclusion.

---

## C. What is and is not rigorous

| Component | Classification | Reason |
|---|---|---|
| SAW correlation domination | rigorous theorem | finite-graph proof above; infinite-volume monotone limit |
| `mu = inf c_n^(1/n)` | rigorous theorem | injection (A.5) plus Fekete |
| `c_0,...,c_14` | computer-assisted rigorous | exhaustive exact-integer backtracking; symmetry weights enumerated exactly |
| OEIS comparison | independent validation | all available terms match the fetched b-file; OEIS is not a proof input |
| square-layer bound | rigorous theorem | GKS-II plus exact 2D critical point |
| infrared inequality and sum rule | rigorous theorem | reflection positivity/Gaussian domination plus Parseval |
| Watson gamma identity | rigorous published identity | Glasser--Zucker, with the corrected denominator 32 |
| interval endpoints | computer-assisted rigorous | directed rounding at 80 decimal digits and outward decimal rounding |
| Bessel quadrature | numerical only | independent 60-digit agreement check |
| `d=4,5` trend | numerical only | quadrature check, not used in the bound |
| `0.221654626(5)` comparison | numerical benchmark only | no fit and no logical role in either endpoint |

The benchmark source is A. M. Ferrenberg, J. Xu, and D. P. Landau, *Pushing
the Limits of Monte Carlo Simulations for the Three-Dimensional Ising Model*,
Phys. Rev. E **97** (2018), 043301,
[doi:10.1103/PhysRevE.97.043301](https://doi.org/10.1103/PhysRevE.97.043301).

## Reproduction

From the repository root:

```text
.venv/bin/python experiments/e05_kc_bounds.py
.venv/bin/python tests/test_bounds.py
```

The experiment writes `results/bounds/kc_bounds.json`.  Its provenance records
`precision: 80`, the exact SAW table, live OEIS comparison, interval endpoints,
all sanity checks, and the factor-eight Watson-formula discrepancy.
