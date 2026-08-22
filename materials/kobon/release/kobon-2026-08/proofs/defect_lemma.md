# n=20, Q=3, T=117: defect-tolerant direction-cycle check

**Verdict:** `PENDING-VERIFICATION: defect-aware Q=3 branch`. The saturated
three-pair and triad cycles are impossible, but the three-unit capacity slack
does not presently imply the saturation hypotheses. At the cycle level there
are exact M2-feasible affine escapes: a paired `b=2` family (including the
`(6,0,8)` witness from `obstruction_n20.md`) and a triad `b=4` witness. No
geometric realization of either affine witness is claimed.

**Data hygiene / no compute.** This is an analytic derivation only. No solver,
interpreter, checker, arrangement construction, or local computation was run
for this artifact. Any proposed finite defect enumeration below is marked
`PENDING-VERIFICATION`; the displayed recurrence calculations are exact hand
algebra.

## 1. What is proved by the capacity ledger

In the no-concurrency branch there are no finite points of multiplicity at
least three. For `n=20`, `Q=3`,

\[
 n(n-2)-2Q=20\cdot18-6=354.
\]

A family of `T=117` triangles has `3T=351` selected-side incidences. With
`C` the number of selected-triangle/interior outside-line incidences and
`e_\ell` the linewise residual,

\[
 e_\ell=m_\ell-\sigma_\ell-\gamma_\ell,\qquad
 E=\sum_\ell e_\ell,
\]

satisfy the **only currently established target-level identity**

\[
 \boxed{E+C=354-351=3},
 \qquad (C,E)\in\{(0,3),(1,2),(2,1),(3,0)\}.
 \tag{1}
\]

For `[2,2,2]+14` the line capacities are `17` on the six paired lines and
`18` on the fourteen singleton lines. For `[3]+17` they are `16` on the
three triad lines and `18` on the seventeen singleton lines. Equation (1) is
a global line/gap ledger; it does not say that every gap is used.

### The symbol `d`

The n14/n18 obstruction documents do not define a separate `d`. If `d` is
intended to count direction-cycle defects, the frequently written expression

\[
 d+E+C=3
 \tag{2?}
\]

is **not a consequence of the capacity theorem**. One must first prove a
charging map from mixed-extreme or endpoint defects to distinct linewise
residual/chord tokens. In the notation below the natural cycle defect is
`b=\sum_i b_i`; no argument currently proves either `b+E+C=3`, `b/2+E+C=3`,
or even `b\le E+C`. Treating (2?) as an identity would double-count (or
reclassify without proof) the same slack already present in (1).

This point matters: `C\ge1` is not itself incompatible with no concurrency.
An outside line may cross the interior of a selected triangle at ordinary
simple crossings. Such a crossing need not be an extreme vertex and need not
change the number of double-extreme crossings between adjacent direction
classes. Likewise, (M2) only bounds double-extreme vertices; it does not bound
`C`. Thus the proposed implication “`C\ge1` contradicts no concurrency +
(M2)” is unavailable without an additional geometric lemma.

## 2. Defect-aware extremity equations

Let `C_0,\ldots,C_{R-1}` be the direction classes in circular order, with
sizes `s_i`. Let `m_i` count double-extreme vertices at
`C_i\cap C_{i+1}`. A mixed-extreme incidence is an extreme incidence that is
extreme on exactly one of its two support lines. Let `b_i` count such
incidences in class `C_i`, and put `b=\sum_i b_i`.

Direction adjacency remains valid for every double-extreme vertex in the
simple (no-three-concurrent) arrangement: it is the same wedge argument as
L4 in the n14/n18 documents. Therefore the exact defect-aware replacement for
(M1), conditional only on retaining adjacency for double extremes, is

\[
 m_{i-1}+m_i=2s_i-b_i,
 \qquad 0\le m_i\le s_i s_{i+1},
 \qquad m_i\in\mathbb Z.
 \tag{M1b--M2}
\]

Summing gives

\[
 2\sum_i m_i=40-b,
 \qquad b\equiv0\pmod2,
 \qquad c:=\sum_i m_i=20-b/2.
 \tag{3}
\]

The saturated Theorem-M case is `b=0`. It is not forced by (1).

For convenient comparison with the prior appendix, write

\[
 y_i=m_i-1,\qquad x_i=s_i-1.
\]

Then

\[
 y_{i-1}+y_i=2x_i-b_i.\tag{4}
\]

The edge bounds in the paired layout are

\[
\begin{array}{c|c|c}
\text{edge}&s_i s_{i+1}&y_i=m_i-1\\ \hline
1\text{--}1&1&-1\le y_i\le0\\
2\text{--}1\text{ or }1\text{--}2&2&-1\le y_i\le1\\
2\text{--}2&4&-1\le y_i\le3,
\end{array}\tag{5}
\]

and in the triad layout the triad-single edge has cap `3`, so
`-1\le y_i\le2`, while singleton-single edges still have
`-1\le y_i\le0`.

## 3. Paired layout `[2,2,2]+14`, `R=17`

There are three doubled direction classes and fourteen singleton classes.
If `r_0,r_1,r_2` are the singleton-run lengths between successive doubled
classes, then

\[
 r_0+r_1+r_2=14.
\tag{6}
\]

### 3.1 The zero-defect row is still closed

For `b=0`, (4) is the M1 system in `obstruction_n20.md`. The hand recurrence
across a run is

\[
 a_{j+1}=2-(-1)^{r_j}a_j,
\tag{7}
\]

where `a_j` is the `y` value on the edge leaving the `j`th doubled class.
Because the total run length is even, either all three runs are even, giving
`a_0=a_1=a_2=1`, or exactly two runs are odd, giving (up to cyclic relabeling)
`(a_0,a_1,a_2)=(-1,1,3)`. In the all-even case, a run long enough to contain
a single-single `y=1` violates (5); keeping every run short gives total at
most `6<14`. In the two-odd case, a positive even run puts `m=4` on a
 doubled-single edge (cap `2`) and then produces `m=-2`; a zero even run still
forces both odd runs to have length `1`, total `2<14`. Hence

\[
 b=0\quad\Longrightarrow\quad\text{paired M2 infeasibility}.
\tag{8}
\]

This is exactly the saturated result, and it applies separately to every
capacity partition in (1) **only if** a separate argument has first shown
`b=0`.

### 3.2 Exact paired `b=2` affine escape family

The danger case cannot be removed by (M1b--M2). The following gives not just
one vector but a hand-parametrized family. Label doubled classes
`p_0,p_1,p_2` cyclically, and let `r_j` be the number of singleton classes
from `p_j` to `p_{j+1}`. Put

\[
 b_{p_0}=2,\qquad b_i=0\ (i\ne p_0).
\tag{9}
\]

Let `a_j` be the departure value `y` on the edge leaving `p_j` and set
`\varepsilon_j=(-1)^{r_j}`. The two ordinary doubled equations and the
one defective equation give

\[
 a_1=2-\varepsilon_0a_0,\qquad
 a_2=2-\varepsilon_1a_1,\qquad
 a_0=-\varepsilon_2a_2,
 \qquad \varepsilon_0\varepsilon_1\varepsilon_2=1.
\tag{10}
\]

Solving (10) exactly:

* if `r_1` is even, `(a_0,a_1,a_2)=(0,2,0)`;
* if `r_1` is odd, `(a_0,a_1,a_2)=(2\varepsilon_0,0,2)`.

In the first case, (5) is satisfied precisely when `r_1=0`: the value
`y=2` then lies on the doubled-doubled edge `p_1p_2`, whose allowance is
`-1\le y\le3`; if `r_1>0`, the first edge is doubled-single and `y=2`
violates its cap.

In the second case, `a_0=2\varepsilon_0` is admissible only if `r_0=0`
(otherwise it is a doubled-single departure of magnitude `2`). Then
`r_0=0` forces `\varepsilon_0=1`, so `a_0=2`; since `r_1` is odd and
`r_0+r_1+r_2=14`, `r_2` is also odd and hence positive. The departure
`a_2=2` is therefore on a doubled-single edge and violates its cap. Thus the
M2-feasible members of this one-class-defect family are exactly

\[
 (r_0,r_1,r_2)=(u,0,14-u),\qquad 0\le u\le14,
\tag{11}
\]

with the defect at the remaining doubled class `p_0`.

For every member of (11), all edges have `m_i=1` except the adjacent
doubled-doubled edge `p_1p_2`, where `m_i=3`. At the three doubled classes,

\[
 1+1=4-2,
 \qquad 1+3=4,
 \qquad 3+1=4,
\tag{12}
\]

and every singleton equation has `1+1=2`. The bounds are exact:
`1\le1` on singleton-single edges, `1\le2` on doubled-single edges, and
`3\le4` on the doubled-doubled edge. In particular, the concrete placement

\[
 (r_0,r_1,r_2)=(6,0,8),\qquad
 (b_{p_0},b_{p_1},b_{p_2})=(2,0,0)
\tag{13}
\]

has departure values `(0,2,0)` and is the witness recorded in the n20
obstruction note.

This is an **exact affine M2 counter-model**, not a geometric arrangement:
it says that endpoint/extreme incidence constraints beyond (M1b--M2) are
needed. It is compatible with the cycle equations and with the numerical
budget `(E,C)=(3,0)` at the level of separate ledgers. Whether a real
117-triangle arrangement can realize the endpoint pattern (13), for any of
the four `(E,C)` rows, is `PENDING-VERIFICATION`.

## 4. Triad layout `[3]+17`, `R=18`

Let `t` be the position of the size-three class. Alternating (M1b) around
the even cycle gives the exact necessary condition

\[
 \sum_{i=0}^{17}(-1)^i b_i=4(-1)^t.
\tag{14}
\]

Indeed, the all-single baseline cancels and replacing `s_t=1` by `s_t=3`
contributes `4(-1)^t` to the alternating sum. Since `b_i\ge0`, (14) implies

\[
 b\ge4.
\tag{15}
\]

Thus a conditional geometric estimate `b\le3` would close the triad branch
immediately. That estimate is not proved by the capacity theorem.

There is also an exact minimum affine escape at `b=4`: put

\[
 b_t=4,\qquad b_i=0\ (i\ne t),\qquad m_i=1\quad(0\le i<18).
\tag{16}
\]

At every singleton class the equation is `1+1=2`; at the triad class it is
`1+1=6-4=2`. The M2 bounds hold (`1\le1` on singleton-singleton edges and
`1\le3` on triad-single edges), and (14) is satisfied for either parity of
`t`. Therefore the first M-vector obstruction to the triad is exactly the
unproved geometric exclusion of four mixed-extreme incidences, not the
alternating equation itself.

Again (16) is only a direction-cycle counter-model. It does not assert that
four mixed extremes can be embedded in a real no-concurrency arrangement with
117 selected triangles.

## 5. The four capacity partitions and what the cycle sees

The exact target partitions are:

\[
\begin{array}{c|c|c|l}
 C&E& E+C&\text{cycle-level consequence}\\
0&3&3&b=0\text{ gives (8); paired }b=2\text{ witness (13) remains}\\
1&2&3&b=0\text{ gives (8)/triad parity; }C\text{ does not enter (M1b--M2)}\\
2&1&3&\text{same; no established map from the two chord tokens to }b_i\\
3&0&3&\text{same; all three units are crossing incidences, but no M2 contradiction follows}\\
\end{array}\tag{17}
\]

For the paired layout, if a future lemma proves the candidate bound
`b\le E+C=3`, parity leaves only `b=0` or `b=2`; (8) closes the first and
(13) is the precise remaining family to exclude geometrically. For the triad
layout, the same candidate bound would give `b\le3`, contradicting (15), so
triad closure is conditional on that one bound.

If the requested notation `d+E+C=3` means that every mixed-extreme defect is
charged injectively (`d=b`), then this is exactly the missing geometric
charging lemma. If `d` is intended to count mixed-extreme *vertices* (`d=b/2`)
or edge corrections instead, no identity has been supplied and the bound must
be re-derived from local gap/XOR geometry. In no interpretation does the
existing M2 inequality alone provide it.

## 6. Exact residual open statements

A complete analytic closure needs both of the following (or an equivalent
stronger statement):

1. **Defect charging:** prove a local injection from the mixed-extreme
   incidences to unused-gap/chord defects, strong enough to imply `b\le3`
   from (1). This would eliminate the triad minimum `b=4`.
2. **Paired endpoint/location lemma:** exclude the paired `b=2` pattern in
   (11), in particular (13). M2 cannot do this: the edge vector is integral,
   obeys every class-product cap, and has the correct total `c=19`.

The smallest exact algebraic escapes that a future check must test are therefore

\[
\begin{array}{c|c|c|c}
\text{layout}&\text{minimum }b&\text{mixed-extreme vector}&\text{M2-feasible }m\\ \hline
[2,2,2]+14&2&b_{p_0}=2,\ (r_0,r_1,r_2)=(u,0,14-u)&m=1\text{ except one adjacent }2\text{--}2\text{ edge }m=3\\
[3]+17&4&b_t=4& m_i=1\text{ for all }i
\end{array}\tag{18}
\]

The geometric status of every row of (18), and of its realizability under each
of the four `(C,E)` partitions in (17), is `PENDING-VERIFICATION`. A future
exact defect-aware appendix must retain linewise capacities `17/18` (paired)
and `16/18` (triad), rather than checking only (M1b--M2).

## 7. Verdict for Main

> **n=20, Q=3, no-concurrency, T=117:** the saturated `b=0` paired
> `R=17` cycle is M2-infeasible and the saturated triad `R=18` cycle fails
> alternating parity, but `E+C=3` does not imply saturation. The exact
> cycle-level residuals are the paired `b=2` family (11), including the
> `(6,0,8)` witness (13), and the triad `b=4` vector (16). No-concurrency
> does not by itself contradict `C\ge1`; a geometric defect-charging lemma
> and a paired endpoint/location lemma remain `PENDING-VERIFICATION`.
