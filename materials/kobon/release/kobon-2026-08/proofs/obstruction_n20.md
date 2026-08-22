# The Q=3 no-concurrency branch of K(20), target T=117

## Status and data-hygiene correction

This is a hand derivation. **No interpreter, checker, solver, build, or other
local computation was run for this artifact.** The exact checks that remain are
listed in §8 as `PENDING-VERIFICATION` items.

There are two different statements to keep separate:

1. The literal target branch is a **near-saturation** branch. Its capacity
   slack is three units, so the saturation hypothesis used by Theorem M at
   `n=14` and `n=18` is not forced here. The honest status of the full target
   branch is therefore `PENDING-VERIFICATION` until a defect lemma (or an exact
   defect-aware check) is supplied.
2. The saturated Theorem-M subcases are nevertheless completely determined by
   the direction-cycle arithmetic. The three-pair subcase is
   `CONDITIONAL-CLOSED-mfeasibility`; the triad subcase is
   `CONDITIONAL-CLOSED-parity`. Those are conditional statements about the
   formal saturated leg, not a silent upgrade of the target's three-unit
   budget to equality.

**Prompt/table correction.** The current `scratch/kobon/n20/signature_algebra.md`
and `scratch/kobon/n20/cube_build.py` define `Q` as the number of parallel
pairs. Under that definition, `n=20,Q=3` has exactly these two class patterns:

* `[2,2,2]+14` singleton classes: three disjoint parallel pairs, hence
  `R=3+14=17` direction classes, an odd cycle;
* `[3]+17` singleton classes: one parallel triad, hence `R=1+17=18`
  direction classes, an even cycle.

A `[3,2]+15` pattern is **pair plus triad**, but it has
`Q=\binom32+\binom22=3+1=4` and is the `q4pairtriad` case, outside this
Q=3 branch. The value `R=15` for a three-pair layout belongs to `n=18`
(`[2,2,2]+12`), not to `n=20`. No `R=15` or pair-plus-triad Q=3 leg is
silently inserted below.

---

## 1. Capacity arithmetic and the three-unit defect budget

The crossing-refined capacity theorem is

$$
3|S|+C\le n(n-2)-2Q-
\sum_{p}k_p(k_p-4),
$$

where `C` counts incidences of selected-triangle interiors with non-supporting
arrangement lines. In the no-concurrency branch there are no finite
multipoints `k_p\ge3`, so the right side is

$$
20\cdot18-2\cdot3=360-6=354.
$$

For a putative `T=117` family,

$$
3T=3\cdot117=351,
\qquad 351+C\le354,
\qquad 0\le C\le3.
\tag{1.1}
$$

The face budget is also harmless in this branch:

$$
F_b=\binom{19}{2}-Q=171-3=168,
\qquad 117\le168.
\tag{1.2}
$$

The three integer-budget rows from `signature_algebra.md` reduce to

$$
C+2Q=C+6\le9,
\qquad Q=3\le54,
\qquad 0\le190-Q=187.
\tag{1.3}
$$

Thus the only active budget is exactly `C\le3`.

For each arrangement line `\ell`, let `q_\ell` be the number of other lines
parallel to it. In the simple finite-point regime, the per-line capacity is

$$
\sigma_\ell+\gamma_\ell\le m_\ell:=20-2-q_\ell,
\tag{1.4}
$$

where `\sigma_\ell` is the number of selected sides on `\ell` and
`\gamma_\ell` is the number of selected-triangle interiors crossed by `\ell`.
Summing (1.4) gives 354. Define the nonnegative linewise deficit

$$
 e_\ell:=m_\ell-\sigma_\ell-\gamma_\ell,
 \qquad E:=\sum_\ell e_\ell.
$$

Since `\sum_\ell\sigma_\ell=3T=351` and
`\sum_\ell\gamma_\ell=C`, we have the exact ledger

$$
\boxed{E+C=354-351=3.}
\tag{1.5}
$$

The four possible total cases are

$$
(C,E)=(0,3),(1,2),(2,1),(3,0).
\tag{1.6}
$$

This is the obstruction's critical difference from `n=18,T=94`, where the
corresponding sum was `282=3T` and forced both `C=0` and `E=0`.

### 1A. Three disjoint parallel pairs: `[2,2,2]+14`

There are six lines in the three size-two parallel classes and fourteen
singleton lines. Therefore

$$
\begin{array}{c|c|c|c}
\text{line type}&\#\text{ lines}&q_\ell&m_\ell\\ \hline
\text{paired}&6&1&17\\
\text{singleton}&14&0&18
\end{array}
$$

and the exact total is

$$
6\cdot17+14\cdot18=102+252=354.
\tag{1.7}
$$

Consequently each paired line obeys `\sigma_\ell+\gamma_\ell\le17`, each
singleton line obeys `\sigma_\ell+\gamma_\ell\le18`, and their deficits obey
(1.5).

### 1B. One parallel triad: `[3]+17`

The three triad lines have `q_\ell=2`; the seventeen singleton lines have
`q_\ell=0`. Hence

$$
\begin{array}{c|c|c|c}
\text{line type}&\#\text{ lines}&q_\ell&m_\ell\\ \hline
\text{triad}&3&2&16\\
\text{singleton}&17&0&18
\end{array}
$$

and

$$
3\cdot16+17\cdot18=48+306=354.
\tag{1.8}
$$

Thus each triad line has cap 16 and each singleton line cap 18, again with
(1.5) as the complete target-level defect ledger.

---

## 2. What saturation would imply (and why it is not forced at T=117)

The gap-and-chord injection gives (1.4). If `C=0`, a selected side cannot run
through an interior vertex: the line through that vertex would cross the
triangle interior. Hence each selected side is an elementary gap and each
selected triangle is an empty cell, exactly as in Lemma E0 of
`scratch/kobon/n18/obstruction_n18.md`.

At `n=18,T=94`, the total side count equaled the total gap capacity, so every
gap was used and the local four-quadrant XOR forced the strict zigzag and
extremity lemmas. Here, however, even the `C=0` row of (1.6) leaves

$$
E=3
$$

unused capacity units. There are three unassigned gap positions after the
side map (counted linewise); they cannot be discarded as if they were zero.
For `C>0`, the chord map consumes some unused gaps and the remaining side-map
defect is `E`, still with `E+C=3`. Therefore the hypotheses

$$
C=0,\qquad E=0,\qquad\text{every bounded gap used}
\tag{2.1}
$$

of the saturated Theorem-M proof are incompatible with `T=117` itself:
`E=0` would require 354 side incidences, i.e. 118 triangles, before any
crossing incidence is counted.

Sections 3--5 record the exact saturated cycle obstruction, because it is the
leg a future defect lemma must recover. Section 6 gives the defect-aware
cycle equations and shows the precise unresolved cases rather than claiming
that (2.1) follows from (1.1).

---

## 3. The saturated direction-cycle equations

Assume temporarily the saturated hypotheses used by Theorem M: every bounded
gap is assigned exactly once, `C=0`, and all local gap/XOR consequences hold.
Let `C_1,\ldots,C_R` be the direction classes in circular order, with class
sizes `s_i`. A double-extreme vertex is an intersection at which both lines
are extreme on their respective lines. The extremity and direction-adjacency
lemmas give

$$
m_i:=\#\{\text{double-extreme vertices at }C_i\cap C_{i+1}\},
$$

with

$$
\tag{M1}m_{i-1}+m_i=2s_i,
\qquad
\tag{M2}0\le m_i\le s_i s_{i+1},
\qquad m_i\in\mathbb Z.
$$

The upper bound in (M2) is the number of available crossings between two
neighboring direction classes; distinct double-extreme vertices use distinct
line pairs.

### 3A. Three-pair layout, `R=17` (odd): the exact M1/M2 leg

There are three doubled positions and fourteen singleton positions. If the
numbers of singleton classes strictly between successive doubled positions
are `r_1,r_2,r_3`, then

$$
r_1+r_2+r_3=14.
\tag{3.1}
$$

Put `x_i=s_i-1\in\{0,1\}` and `y_i=m_i-1`. Then (M1) is

$$
y_{i-1}+y_i=2x_i.
\tag{3.2}
$$

Let `a_j` be the `y` value on the edge leaving the `j`th doubled position.
Across a run of `r_j` singleton classes,

$$
a_{j+1}=2-(-1)^{r_j}a_j.
\tag{3.3}
$$

Because the total in (3.1) is even, there are two cases.

* **All three runs even.** The odd number of doubled positions makes the
  recurrence have the unique solution `a_1=a_2=a_3=1`. Along an even run the
  `y` values alternate `1,-1,1,-1,\ldots`; a `y=1` on a singleton-singleton
  edge means `m=2`, violating the upper bound `m\le1` in (M2). Avoiding such
  an edge forces every even run to have length at most 2, so their total is at
  most 6, contradicting (3.1).

* **Exactly two runs odd and one run even.** Up to cyclic relabeling, (3.3)
  gives `(a_1,a_2,a_3)=(-1,1,3)`, hence edge multiplicities `(m_1,m_2,m_3)
  =(0,2,4)` on the doubled-class departures. If the even run has length
  `e>0`, its first edge is doubled-single and has cap 2, whereas `m=4`;
  the following alternating value is `m=-2`. If `e=0`, the `m=4` edge is
  doubled-doubled and its cap 4 is not by itself a violation, but each odd
  run must then have length 1 to avoid a singleton-singleton edge with `m=2`.
  This gives total run length at most 2, again contradicting (3.1).

Thus the saturated M1/M2 system is infeasible for every placement. The raw
placement ledger is hand-checkable:

$$
\#\text{placements}=\binom{17}{3}=680.
\tag{3.4}
$$

The ordered rooted gap compositions have `\binom{16}{2}=120` members; dividing
`17\cdot120` by the three choices of distinguished doubled class gives 680.
The parity subtotals are

$$
\begin{array}{c|c|c|c}
\text{run type}&\#\text{ordered compositions}&\text{raw placements}&\text{leg}\ \hline
\text{all even}&\binom92=36&17\cdot36/3=204&\text{M2: }m=2>1\\
\text{two odd, even run }e=0&21&17\cdot21/3=119&\text{M2: }m=2>1\\
\text{two odd, }e>0&63&17\cdot63/3=357&\text{M2: }m=4>2\text{ and/or }m<0\\ \hline
\text{total}&120&680&\text{all infeasible}
\end{array}
\tag{3.5}
$$

This is the corrected `n=20` analogue of the odd-cycle leg; `R=17`, not
`R=15`. It is a `CONDITIONAL-CLOSED-mfeasibility` result **only under
the explicitly stated saturation hypothesis**.

### 3B. Triad layout, `R=18` (even): the exact saturated parity leg

There is one size-three class and seventeen singleton classes. Place the triad
at position `t\in\{0,\ldots,17\}`. The necessary alternating sum of (M1) is

$$
0=\sum_{i=0}^{17}(-1)^i(m_{i-1}+m_i)
 =\sum_{i=0}^{17}(-1)^i2s_i.
\tag{3.6}
$$

The all-single baseline cancels on an even cycle. Replacing one `s_t=1` by
`s_t=3` adds 4 at position `t`, so

$$
\sum_{i=0}^{17}(-1)^i2s_i=4(-1)^t\in\{+4,-4\}\ne0.
\tag{3.7}
$$

Thus all 18 triad placements have no rational M1 solution. The saturated
triad leg is therefore `CONDITIONAL-CLOSED-parity`, only under (2.1).

---

## 4. Defect-aware extremity ledger: the missing implication

For the actual target, let `b_i` count extreme incidences in direction class
`C_i` that are mixed-extreme (extreme on exactly one support line), rather than
double-extreme. If all double-extreme vertices still obey direction adjacency,
the natural defect-aware replacement of (M1) is

$$
m_{i-1}+m_i=2s_i-b_i,
\qquad \sum_i b_i=b,
\tag{4.1}
$$

where `b` is the total number of mixed-extreme incidences. The extremity count
also gives `2c+b=40`, with `c` the number of double-extreme vertices, hence
`b` is even.

Equation (4.1) makes the unresolved slack visible.

### 4A. Paired layout with `R=17`

Reducing (4.1) modulo 2 around the odd cycle gives

$$
0\equiv\sum_i b_i=b\pmod2.
\tag{4.2}
$$

The budget `E+C=3` does **not**, by itself, prove `b=0`; no lemma in the
current n20 materials bounds mixed-extreme incidences by the three defect
units. If one could prove `b\le3`, only `b=0` or `b=2` would remain. The
`b=0` row is exactly the infeasible M1/M2 analysis of §3A.

The `b=2` row cannot be erased by M-vector arithmetic alone. For an explicit
hand-checkable feasibility ledger, place doubled classes at positions
`0,7,8`, so the singleton runs have lengths `(6,0,8)`. Set

$$
(b_0,b_7,b_8)=(2,0,0),\qquad b_i=0\text{ elsewhere}.
$$

Writing `y_i=m_i-1`, the recurrence gives departure values

$$
(a_0,a_7,a_8)=(0,2,0).
$$

Therefore `m_i=1` on every edge except the edge between the adjacent doubled
classes at positions 7 and 8, where `m_i=3`. This obeys every M2 bound:

* singleton-singleton edges have `m=1\le1`;
* doubled-single edges have `m=1\le2`;
* the doubled-doubled edge has `m=3\le4`;
* at the three doubled classes, the sums are respectively
  `1+1=4-2`, `1+3=4`, and `3+1=4`.

This is not a geometric existence claim; it is an exact demonstration that
an M-vector check cannot close the `b=2` defect row without an additional
geometric bound on where mixed extremes may occur. Consequently the actual
q3paired target branch is currently

$$
\boxed{\texttt{PENDING-VERIFICATION: defect-aware paired leg}}
$$
while its saturated subcase is `CONDITIONAL-CLOSED-mfeasibility` by §3A.

### 4B. Triad layout with `R=18`

Applying the alternating sum to (4.1) yields

$$
0=4(-1)^t-\sum_{i=0}^{17}(-1)^i b_i,
\qquad
\sum_{i=0}^{17}(-1)^i b_i=4(-1)^t.
\tag{4.3}
$$

Thus at least four mixed-extreme incidences are needed by the defect-aware
cycle equation. If a future geometric defect lemma proves `b\le3` from
`E+C=3`, the triad branch closes immediately on the parity leg. That
inequality is not proved in the current templates: a mixed-extreme vertex can
be made possible precisely when the all-gaps-used argument has a defect at a
flanking gap, and the present capacity theorem only counts defects globally.

The exact current status is therefore

$$
\boxed{\texttt{PENDING-VERIFICATION: defect-aware triad leg}}
$$
with the saturated subcase `CONDITIONAL-CLOSED-parity` by §3B.

---

## 5. Branch verdict table

| branch | parallel classes | `R` | per-line total | saturated cycle result | actual `T=117` status |
|---|---|---:|---:|---|---|
| `q3paired` | `[2,2,2]+1^{14}` | 17 odd | `6·17+14·18=354`; `E+C=3` | `CONDITIONAL-CLOSED-mfeasibility` (680 placements, M2 failure) | `PENDING-VERIFICATION` (defect rows `b=2` not ruled out) |
| `q3triad` | `[3]+1^{17}` | 18 even | `3·16+17·18=354`; `E+C=3` | `CONDITIONAL-CLOSED-parity` (18 placements, alternating sum `±4`) | `PENDING-VERIFICATION` (need a defect lemma, e.g. `b\le3`) |
| `q4pairtriad` (out of scope) | `[3,2]+1^{15}` | 17 | `q=2,1,0` capacities; Q=4 | not part of Q=3 proof | not covered here |

The precise verdict line for Main is:

> **n=20, Q=3, no-concurrency, T=117:** the saturated three-pair leg dies by
> odd-cycle M2 infeasibility (`R=17`, all 680 placements), and the saturated
> triad leg dies by alternating-sum parity (`R=18`, all 18 placements). The
> actual target branch has three units of capacity slack (`E+C=3`), so both
> full branches remain `PENDING-VERIFICATION` until the defect-aware extension
> is proved; pair-plus-triad is Q=4 and out of scope.

---

## 6. What would constitute a complete defect closure

A valid unconditional closure must not assume `E=0`. It must do one of the
following, separately or together:

1. prove a geometric inequality converting the global defect ledger `E+C=3`
   into a bound on mixed-extreme incidences (`b`, and the locations of the
   associated non-adjacent extremes); or
2. perform an exact finite defect-aware cycle check whose input ledger is
   explicitly all four cases `(C,E)=(0,3),(1,2),(2,1),(3,0)`, with linewise
   capacities 17/18 in the paired layout and 16/18 in the triad layout.

A check that tests only (M1)+(M2), or only the all-gaps-used case, does not
cover the target branch and must not be reported as a Q=3/T117 obstruction.

---

## 7. Hand-arithmetic audit

All numerical identities used above are:

$$
20(20-2)=360,\quad 2Q=6,\quad 360-6=354,\quad 3T=351,
$$

$$
354-351=3,\quad \binom{19}{2}-3=168,
$$

$$
6(20-2-1)+14(20-2)=6\cdot17+14\cdot18=354,
$$

$$
3(20-2-2)+17(20-2)=3\cdot16+17\cdot18=354,
$$

$$
R_{\rm paired}=3+14=17,\qquad R_{\rm triad}=1+17=18,
$$

$$
\binom{17}{3}=680,
\quad
204+119+357=680,
\quad
4(-1)^t\ne0.
$$

No machine-generated count is used in these rows.

---

## 8. `PENDING-VERIFICATION` queue (nothing run here)

1. **Defect-aware paired cycle check.** Enumerate the four exact defect totals
   `(C,E)=(0,3),(1,2),(2,1),(3,0)`. For each, enumerate all admissible
   distributions of linewise deficits over the six cap-17 and fourteen cap-18
   lines, all mixed-extreme incidence vectors `b_i` allowed by the local
   geometry, and all direction placements of the three doubled classes. Check
   (4.1), `0\le m_i\le s_i s_{i+1}`, and the endpoint/extreme incidence
   constraints. The check must reproduce the saturated subtotals
   `204+119+357=680` before adding defects. In particular it must report
   whether the hand-feasible arithmetic `b_0=2` witness in §4A survives the
   geometric endpoint constraints.

2. **Defect-aware triad cycle check.** For each triad position
   `t=0,\ldots,17`, enumerate the same `(C,E)` partitions and require the
   signed equation (4.3), together with M2 and the linewise 16/18 capacities.
   If the geometry proves `b\le3`, report the immediate parity closure from
   (4.3); otherwise retain any surviving row as a genuine pending branch.

3. **Proof-interface check.** Independently verify which of E0/E1/L2/L3/L4
   survives with `E+C=3`; in particular, do not import the n14/n18 assertion
   “every bounded gap is used” into n20/T117. Controls should include the
   all-single odd/even cycles and the feasible affine-hexagon M-vector control.

No item above was executed or silently approximated in preparing this file.
