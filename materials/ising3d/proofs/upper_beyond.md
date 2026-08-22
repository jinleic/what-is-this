# Beyond the two-point LP: the Peierls contour route to the upper endpoint is
# certified impossible at the incumbent (limitation theorem)

**File classification.** [COMPUTATION] exact enumeration data (fixed polycube
census through \(n\le10\), contour-area census through \(A\le28\), bond strata,
tree strata, cavity certificate).  [THEOREM] the plain Peierls contour-sum
certificate of `proofs/kc_upper_peierls.md`, eq. (2.5)–(2.6), **fails for
every** \(K\le K^{*}\) with a certified \(K^{*}>\bar K=I_3/2\), decided by an
exact rational inequality at the incumbent coupling; hence the exact-count
head-extension program can neither reach nor beat the certified incumbent
upper endpoint \(\bar K=I_3/2\) along this route.  [COMPUTATION] the wave-6
exponential tail \(N(A)\le(78/625)\lambda^{A}\) is machine-reproved and its
head-plus-tail certified endpoint improved, with the certified conclusion
that this specific tail can never approach the incumbent (its convergence
threshold \((1/2)\log\lambda\) is a certified \(1.42\ldots\) above \(\bar K\)).

Producer: `experiments/e140_upper_beyond.py`.  Standalone auditor:
`tests/test_upper_beyond.py` (regenerates the complete census through
\(n=10\) with its own engine, plus every threshold/residual value, and
re-decides every inequality on exact `Fraction` endpoints produced by a
*different* transcendental enclosure).  Artifact: `results/bounds/upper_beyond.json`.

> **Constraint-class statement.**  The audited reflection-positivity/infrared
> two-point constraint class is exactly optimal at \(I_3/2\)
> (`proofs/upper_infrared.md`, method-optimality theorem) and its
> zero-magnetization profile saturates the constraint with torus LP floor
> \(O(1/L)\) (`proofs/mag_floor.md`).  It is therefore *ruled out* for any
> improvement, and this front offers **no two-point-only variant**: all work
> below is strictly in the Peierls outer-contour counting class.

---

## 1. The contour model (recall)

For the upper-bound direction, wave 6 (`proofs/kc_upper_peierls.md` §§2–3)
proved the complete Peierls inequality chain: if

\[
  \sum_{A\ge6} N(A)\,e^{-2KA}<\tfrac12, \tag{1.1}
\]

then the plus state has positive spontaneous magnetisation at \(K\), hence
\(K_c\le K\).  Here \(\gamma\) runs over outer contours,
\(A=|\gamma|\) is the contour area, and

\[
  N(A)=\sum_{P} n(P), \tag{1.2}
\]

the sum over translation classes \(P\) of finite face-connected cell sets
with face-connected complement and boundary area \(A\) (rotations/reflections
distinct), \(n(P)=|P|\): a class with \(n\) cells contributes exactly
\(n\) cell-rootings at the fixed origin (§3.1 of the wave-6 note).  Hence the
*unweighted* class count is \(c(n,b)\) with \(b\) the number of internal bonds,
and

\[
  N(A)=\sum_{6n-2b=A} n\,c(n,b). \tag{1.3}
\]

## 2. The exact census through area 28

### 2.1 Data

All counts are translation classes on \(\mathbb Z^3\) (orientations distinct).
The census is a breadth-first canonical growth (translation-minimum
representative packed as a single integer), evaluated through \(n=10\) cells:

\[
\begin{array}{c|rrrrrrrrrr}
n&1&2&3&4&5&6&7&8&9&10\\\hline
\#\text{classes}&1&3&15&86&534&3481&23502&162913&1152870&8294738
\end{array}
\]

matching OEIS A001931 through \(n=10\).  The full bond strata \(c(n,b)\) are
stored in the artifact (`census.bond_strata`, 32 rows).

**Area completeness (discrete Loomis–Whitney).**  Every \(n\)-cell class has
boundary area \(A\) with \(A^3\ge216n^2\) (wave-6 note, §3.2, eq. (3.1)).
Since \(216\cdot10^2=21600\le28^3=21952<26136=216\cdot11^2\), the enumeration
through \(n=10\) is **complete for every area \(A\le28\)**; similarly
\(216\cdot9^2\le26^3<216\cdot10^2\) certifies completeness through \(n=9\)
for \(A\le26\).

**Cavity certificate.**  A class contributes to \(N(A)\) *iff* its complement
is face-connected (wave-6 note §3.1).  The producer ran an exact flood fill of
the complement of **every one of the 9 638 143 level-\(n\) classes,
\(n\le10\), finding \(0\) cavity-bearing classes**, so
\(c(n,b)=c^{\mathrm{cf}}(n,b)\) for all \(n\le10\) and the unfiltered table is
the contour table.  This is sharper than needed for the head, and is also the
load-bearing hypothesis for the tail components below (a cavity-bearing class
would inflate the lower bounds).  The result is independently sharpened:

**[LEMMA: no cavity below 11 cells; a tree witness at 11].**  Let \(P\) be a
connected cell set whose complement has a bounded component \(D\).
\(N(D)\subseteq P\) (the neighbourhood of \(D\) must be occupied), and

* \(|D|=1\): the 6 spokes of the enclosed cell lie in \(P\); the exact frontier
  enumeration from the six-spoke base (adding one adjacent cell at a time,
  complete by the reverse-spanning-tree argument) finds *no connected*
  spoke superset on \(\le10\) cells at all, and an 11-cell connected one —
  the smallest single-cell cavity needs exactly 11 cells;
* \(|D|=2\): \(P\supseteq N(\text{domino})\), \(|N(\text{domino})|=10\), and
  \(N(\text{domino})\) is disconnected, so \(|P|\ge11\);
* \(3\le|D|\le5\): \(D\) contains a connected 3-subset \(D_3\);
  \(\min|N(D_3)|=13\) over the 15 three-cell classes (exhaustive), so
  \(|P|\ge|N(D)|\ge13-(|D|-3)\ge11\);
* \(|D|\ge6\): excluded by the direct computation above (zero cavity-bearing
  classes with \(\le10\) cells).

Moreover the 11-cell tree \(P_0\) of the artifact (`CAVITY_WITNESS`) — the six
spokes plus five corner connectors, \(b=10=n-1\) — has a disconnected
complement, so the bound is attained. \(\square\)

**Tree strata and their closed form.**  The bond strata with \(b=n-1\) are
the tree classes

\[
T_n=1,\ 3,\ 15,\ 83,\ 486,\ 2967,\ 18748,\ 121725,\ 807381,\ 5447203
\qquad(n=1,\dots,10), \tag{2.2}
\]

counted twice: as strata of the census and by an **independent leaf-extension
enumeration** (grow by adjoining a cell adjacent to *exactly one* existing
cell — the candidate has exactly one new bond, so the result is precisely the
tree classes).  The two enumerations agree at every \(n\le10\).  Every tree
class has area

\[
  A=6n-2(n-1)=4n+2, \tag{2.3}
\]

an exact closed form for the tree contribution: the tree part of \(N(A)\)
is \(n\,T_n\) with \(n=(A-2)/4\) whenever \(A\equiv2\pmod4\).

**Area census.**  The exact rooted counts:

\[
\begin{array}{c|rrrrrrrrrrr}
A&6&8&10&12&14&16&18&20&22&24&26\\\hline
N(A)&1&0&6&0&45&12&332&240&2538&3040&20448
\end{array}
\]

reproduces the wave-6 artifact `results/bounds/peierls_upper.json`
(`exact_head.N_of_A`) **verbatim** through \(A=24\), and adds

\[
  N(26)=6\cdot2967+7\cdot378=20448, \qquad
  N(28)=7\cdot4368+8\cdot306+9\cdot24=33240. \tag{2.4}
\]

(\(N(26)\) has an \(n=6\) tree part \(6T_6\) and a one-cycle part
\(7\,c(7,8)\); \(N(28)=7\,c(7,7)+8\,c(8,10)+9\,c(9,13)\), with
\(c(10,16)=0\) — the census rows whose Loomis–Whitney completeness at
\(A\le28\) is certified above.)

### 2.2 The certified tail components (exact lower bounds for the tail)

For the areas \(A=28,30,\dots,42\) the enumeration is *not* complete
(\(n\ge11\) classes contribute), but the \(n\le10\) components give exact
lower bounds

\[
  M_A:=\sum_{n\le10,\,6n-2b=A} n\,c(n,b)\ \le\ N(A), \tag{2.5}
\]

with the full per-\(n\) decompositions stored in the artifact
(`tail_lower_bound_components`); the values are

\[
\begin{array}{c|rrrrrrrr}
A&28&30&32&34&36&38&40&42\\\hline
M_A&33240&171739&342072&1484100&3321028&12156189&22807920&54472030
\end{array}
\]

(\(M_{28}=N(28)\) exactly; \(M_{42}=10T_{10}\) is the pure tree part).

## 3. Exact transcendental arithmetic

Every proof step is exact rational arithmetic.  The only transcendental input
is \(x=e^{-2K}\) at rational \(K\) (and \(x^{*}=e^{-I_3}\)), enclosed as
follows.

**[LEMMA: alternating-enclosure correctness]** For \(0<u\le1/2\) the
alternating Taylor series of \(e^{-u}\) has strictly decreasing absolute
terms, so every odd partial sum is a strict lower bound and every even
partial sum a strict upper bound.  Subdividing any rational \(t>0\) as
\(t=Nu\) with \(N=\lceil 2t\rceil\) and \(u\le1/2\), Taylor order 44, and
raising to the \(N\)-th power (monotone on positives) encloses
\(e^{-t}\) between two exact `Fraction`s. \(\square\)

The \(I_3\) anchor is the repo-certified interval of
`results/bounds/upper_infrared.json` (`certified_constants.I3`), intersected
with a fresh `mpmath` interval-arithmetic evaluation at `iv.dps=100` of

\[
  I_3=\frac{\sqrt6\;\Gamma(\tfrac1{24})\Gamma(\tfrac5{24})
            \Gamma(\tfrac7{24})\Gamma(\tfrac{11}{24})}{96\pi^3}, \tag{3.1}
\]

the fresh interval nesting inside the repo interval widened by \(10^{-90}\).
The certified enclosure used for proof steps is

\[
  x^{*}\in[x_{-},x_{+}],\qquad
  x_{-}=0.603226808354577933805005459382\ldots,\quad
  x_{+}=0.603226808354577933805005459383\ldots \tag{3.2}
\]

with width \(<10^{-60}\), and it overlaps the independent `mpmath` enclosure
of the same quantity.

**Standalone auditor.**  The test makes its own enclosures by a *different*
formula: \(e^{-t}=(e^{u})^{-N}\) with the all-positive series for \(e^{u}\)
truncated at order 50 and the exact geometric remainder bound
\(\mathrm{rem}\le\frac{u^{51}}{51!}(1-\tfrac{u}{52})^{-1}\);
every threshold below is re-decided on the test's own endpoints, and each
certificate interval produced by the producer is asserted to overlap the
test's own.  Both routes are pure exact `Fraction` arithmetic (no floating
point in any proof step; `mpmath` appears only as a cross-check and for
anchoring, never as a decider).

## 4. The limitation theorem

**[THEOREM P — certified failure of the Peierls certificate at the incumbent]**
Let \(\bar K=I_3/2\) (the certified incumbent upper endpoint) and
\(x^{*}=e^{-2\bar K}=e^{-I_3}\).  Define the certified partial contour sum

\[
  V(x):=\sum_{A\le26}N(A)\,x^{A}+\sum_{A\in\{28,30,\dots,42\}}M_A\,x^{A},
  \qquad V(x)\le\sum_{A\ge6}N(A)x^{A}, \tag{4.1}
\]

the inequality by (2.5) and exactness of the head (§2.1) with every summand
nonnegative.  Then on the exact `Fraction` endpoints of (3.2),

\[
  V(x^{*})\ \ge\ V(x_{-})
  \ =\ \underbrace{\sum_{A\le26}N(A)x_{-}^{A}}_{=\,0.2691403125\ldots}
  \ +\ \underbrace{\sum_{28\le A\le42}M_Ax_{-}^{A}}_{=\,0.3191792150\ldots}
  \ =\ 0.5883195276\ldots\ >\ \tfrac12, \tag{4.2}
\]

the last inequality a certified exact rational comparison with margin

\[
  V(x_{-})-\tfrac12
  =0.08831952760029326734\ldots\;>\;0. \tag{4.3}
\]

In particular the tail mass provably exceeds the head slack:
\(\sum_{A\ge28}M_Ax_{-}^{A}=0.3191792150\ldots>
\tfrac12-\sum_{A\le26}N(A)x_{-}^{A}=0.2308596874\ldots\).

**Consequences.**

1. **The certificate fails at the incumbent.** By (2.5),
   \(\sum_AN(A)e^{-2\bar KA}\ge V(x_{-})>1/2\); condition (1.1) is false at
   \(K=\bar K\).  No sharpening of the *counts* can repair it: the exact head
   values through \(A=26\) are final (completeness), and every additional
   exact or estimated term only *increases* the sum.
2. **It fails on a certified whole interval above the incumbent.**  Let
   \(K^{*}:=\sup\{K:\text{the certified partial sum at }K\text{ exceeds }1/2\}\).
   \(V(e^{-2K})\) is decreasing in \(K\) (positive polynomial in \(x=e^{-2K}\)),
   and the exact 3-way bisection certifies
   \[
     K^{*}\in[\,0.255832613172027451\ldots,\ 0.255832613172027452\ldots\,],
     \qquad V(e^{-2K})>\tfrac12\quad(K\le K^{*}_{-}), \tag{4.4}
   \]
   with \(V(e^{-2K^{*}_{+}})<1/2\) certified on the upper end, and a certified
   positive gap above the incumbent:
   \[
     K^{*}_{-}-\bar K_{\max}
     =0.0031016033133644\ldots\;>\;0. \tag{4.5}
   \]
   So the plain certificate fails for EVERY \(K\le K^{*}_{-}\), hence on a
   certified positive-measure interval strictly above the incumbent.
3. **The head-truncation crossing is spurious.**  Solving head-only
   \(\sum_{A\le26}N(A)x^{A}=1/2\) produces a crossing
   \(K_{26}\in[\,0.234693300322076899\ldots,\ 0.234693300322076900\ldots\,]\),
   *well below* \(\bar K_{\min}\) — a feigned gain of
   \(\bar K_{\min}-K_{26}=0.018037709536\ldots\).
   At \(K=K_{26}\) the certified *full* partial sum still satisfies
   \(V(e^{-2K_{26}})\ge1.6528764620\ldots>\tfrac12\): the proved tail mass
   annihilates the feigned gain.  This is the precise sense in which **the
   proved tail error exceeds any possible gain at \(n=26\)**.

**[THEOREM: quantitative shape of the wave-6 tail route]**  The wave-6
exponential tail is machine-reproved here (Lagrange-inversion identity
\([z^{A-1}]z(1+T)^{12}=\frac{12}{A-1}\binom{11A}{A-2}\) verified by iterative
series reversion for \(A=2..33\), and the exponential majorants
\(100\binom{11A}{A-2}\le\lambda^{A}\),
\(625AB_A\le78\lambda^{A}\), \(\lambda=11^{11}/10^{10}\), verified over
exact-integers rows; cf. wave-6 §4 for the all-\(A\) proofs via the binomial
theorem).  With the extended exact head through \(A=28\),

* the least \(10^{-6}\)-grid coupling \(K_g\) with
  \(\sum_{A\le28}N(A)e^{-2KA}+\frac{78}{625}\frac{(\lambda e^{-2K})^{30}}
  {1-(\lambda e^{-2K})^2}<\frac12\) is
  \[
    K_g=\frac{423823}{250000}=1.695292,
    \qquad K_0-K_g=\frac{803}{500000}=0.001606, \tag{4.6}
  \]
  improving the wave-6 certified endpoint \(K_0=848449/500000=1.696898\)
  (the grid point one step below fails certifiedly);
* nevertheless the majorant converges only when
  \(\lambda e^{-2K}<1\iff K>\tfrac12\log\lambda\), with
  \[
    \tfrac12\log\lambda\in[\,1.675498535420\ldots,1.675498535421\ldots\,],
    \qquad \tfrac12\log\lambda-\bar K_{\max}=1.422767525562\ldots, \tag{4.7}
  \]
  a certified \(1.42\ldots\) *above* the incumbent \(\bar K\): **no finite
  exact head can make this exponential tail argument competitive** (wave-6
  §5, restated with quantitative gap).

## 5. What is *not* claimed

* No improved upper endpoint: this is (b) of the assignment — the limitation
  theorem with exact numbers, not (a).
* The true series \(\sum_AN(A)(x^{*})^A\) need not diverge; the theorem
  certifies only that a certified finite partial sum of it already exceeds
  \(1/2\).  Whether the full series converges at \(x^{*}\) is **open**.  The
  tree subseries alone converges iff
  \(\lim T_{n+1}/T_n<(x^{*})^{-4}\approx7.552\); the successive tree ratios
  through \(n=10\) rise to \(T_{10}/T_9=5447203/807381\approx6.747\) — below,
  but with no proved limiting behaviour.
* The contour-sum threshold \(K^{*}\) is a property of the *certificate*, not
  of \(K_c\): no statement about the location of \(K_c\) relative to
  \(K^{*}\) is made.

## 6. Reproduction

```
PYTHONPATH=src .venv/bin/python experiments/e140_upper_beyond.py   # regenerate artifact
PYTHONPATH=src .venv/bin/python tests/test_upper_beyond.py         # independent audit
```

The audit regenerates the entire census through \(n=10\) (own canonical form,
own growth, own flood fill), the tree census through \(n=10\) (own leaf
extension), all threshold values (own transcendental enclosures), and the
wave-6 table verbatim.
