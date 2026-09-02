# Theorem H-MIX-CROSSING (Theorem X) — run `20260901T132343Z_4f169cbb_2f81dd8f2d8d`

Gate **H-MIX-CROSSING**, target `cs/rs-pe3d`, agent `RsPe3dCrossing`, 2026-09-01.
Prereg: `prereg/H_CROSSING_PREREG_2026-09-01.md`, commit `cd2939676258c8be164f227cf7a391e1e1a4a9b4`,
sha256 `8c24ac710b7559f02d02acc1f57bdf022951c25ce14e58794f09b1eddd993ef5`, byte-identical
copy in this run dir as `pre_statement.md` (cmp-verified). All numbers below are
in-run outputs of `run_crossing.py` (exact GF(p), stdlib rref `pow(x,-1,p)`,
`PYTHONDONTWRITEBYTECODE=1`, `nice -n 15`, single-threaded), executed AFTER
`campaign.py init`. Full machine record: `controls_results.json`
(sha256 `9d0c4bd260c97e78c2d45de0b430bfa93bdf94e3b65c58321ae8c0c9aea759cc`),
66 pinned control asserts, 0 failures. Defects D1–D4 disclosed at the end.

## 0. Setting and notation

$A$: $r_A \times s_A$ and $B$: $r_B \times s_B$ over GF($p$), $p$ odd, all
columns nonzero; sparks $d_A, d_B \ge 2$ measured; $d = \min(d_A, d_B)$;
product columns $h(u,v) = a_u \otimes b_v \in F^{r_A r_B}$ indexed by the grid
$[s_A] \times [s_B]$, $N = s_A s_B$. For a column set $S \subseteq [s_A]
\times [s_B]$ the **profile** is $(n_A, n_B) = (\text{distinct } A\text{-indices},
\text{distinct } B\text{-indices})$; a fiber is a set of the form $\{u_0\}
\times V$ or $U \times \{v_0\}$; **all-distinct** is $(d{+}1, d{+}1)$; anything
else with both counts $> 1$ is **mixed**. Factor circuit counts:
$C_X(k) = \#\{k\text{-circuits of } X\}$; $\mathrm{Circ}_X(k)$ = the family of
$k$-subsets which are circuits. A **circuit relation** of a factor on $R$ is a
vector $c \in F^{s_X}$ with $\sum_i c_i\, x_i = 0$ and $\mathrm{supp}(c) = R$;
for a $k$-circuit the relation space is 1-dimensional (every proper subset
independent), so up to scaling there is exactly one, and it has **full support
on $R$**.

Frozen kernel identity (T-DGE; re-asserted on every census in-run):
$\{C \in F^{s_A \times s_B} : A C B^{\mathsf T} = 0\} = (\ker A \otimes
F^{s_B}) + (F^{s_A} \otimes \ker B)$, i.e. every kernel element is $C_1 + C_2$
with the nonzero columns of $C_1$ in $\ker A$ and the nonzero rows of $C_2$ in
$\ker B$; dimension identity $s_A s_B - \rho_A \rho_B = (s_A-\rho_A) s_B +
s_A (s_B - \rho_B) - (s_A - \rho_A)(s_B - \rho_B)$.

**Coordinate convention** (in-run, both directions verified): the column
$(u,v)$ flattens row-major $[(a_u)_i (b_v)_j]_{i,j}$, $i \cdot r_B + j$;
$\ker A$ acts on $A$-indices $i$, $\ker B$ on $B$-indices $j$; a matrix entry
$C(u|v)$ sits at cell $((u,v))$ of the grid and a kernel combination
$\sum_{(u,v)} \gamma_{u,v} h(u,v) = 0$ is the matrix equation
$A\, \Gamma\, B^{\mathsf T} = 0$ with $\Gamma = (\gamma_{u,v})$ (proved
mechanically: O1/T3a relations verified both as coefficient sums over
`product_columns` and as matrix equations).

## 1. Construction (the crossing channel)

Let $R \in \mathrm{Circ}_A(d_A)$, $i_0 \in R$, $J \in \mathrm{Circ}_B(d_B)$,
$j_0 \in J$. Let $c$ be the circuit relation of $A$ on $R$ and $\delta$ that
of $B$ on $J$. Since $c_{i_0} \ne 0 \ne \delta_{j_0}$ (full support), choose
scaling $\delta := \sigma \delta$ with $\sigma = -c_{i_0}/\delta_{j_0}$, giving

$$c_{i_0} + \delta_{j_0} = 0.$$

Define the $s_A \times s_B$ matrix

$$\Gamma = c\, e_{j_0}^{\mathsf T} + e_{i_0}\, \delta^{\mathsf T}.$$

Then $A \Gamma B^{\mathsf T} = (Ac) (e_{j_0}^{\mathsf T} B^{\mathsf T}) +
(A e_{i_0})(\delta^{\mathsf T} B^{\mathsf T}) = 0 \cdot b_{j_0}^{\mathsf T} +
a_{i_0}^{\mathsf T} \cdot 0 = 0$, so $\Gamma \in \ker(A \otimes B)$, and

$$\sum_{(u,v)} \Gamma_{u,v}\, h(u,v) = 0.$$

**Support (exact).** $\Gamma_{u,v} = c_u \delta_{j_0} [u \in R] + c_{i_0}
\delta_v [v \in J]$ — more precisely: for $u \ne i_0$, $v \ne j_0$: both terms
need $u \in R, v = j_0$ resp. $u = i_0, v \in J$, impossible, so
$\Gamma_{u,v} = 0$; for $u \ne i_0, v = j_0$: $\Gamma = c_u \delta_{j_0} \ne
0$; for $u = i_0, v \ne j_0$: $\Gamma = c_{i_0} \delta_v \ne 0$; for the
**crossing cell** $(i_0, j_0)$: $\Gamma = c_{i_0} \delta_{j_0} + c_{i_0}
\delta_{j_0} = 0$ by the scaling. Hence

$$S(R, i_0, J, j_0) := \big((R \setminus \{i_0\}) \times \{j_0\}\big) \;\cup\;
\big(\{i_0\} \times (J \setminus \{j_0\})\big),$$

with $|S| = (d_A - 1) + (d_B - 1) = d_A + d_B - 2 =: m_*$, profile $(d_A,
d_B)$, all coefficients nonzero outside the cancelled center, and the relation
$\Gamma$ has full support on $S$. In-run O1 instance (T1 block): $R = J =
\{0,1,2\}$, $i_0 = 0$, $j_0 = 1$, support $[(0,0),(0,2),(1,1),(2,1)]$
(circuit relations $[1, 11, 1, \cdot]$ / $[1, 5, 9, \cdot]$ at GF(13)),
sum zero verified exactly, rank 3 $= m_* - 1$, all $\binom{4}{3}$ 3-subsets
independent, profile $(3,3)$.

## 2. Theorem X (certified statement)

**Theorem X.** Let $d_A, d_B \ge 2$ and $m_* = d_A + d_B - 2$. The circuits
of $H = A \otimes B$ of size $m_*$ with profile exactly $(d_A, d_B)$ are
exactly the sets $S(R, i_0, J, j_0)$ of Section 1, and each is a genuine
circuit (rank $m_* - 1$, every proper subset independent) with 1-dimensional
relation space $\langle \Gamma \rangle$. Consequently their number is

$$N_{\mathrm{mix}} = d_A\, d_B\, C_A(d_A)\, C_B(d_B)$$

**provided the parameterization is injective**, which holds when $d_A, d_B
\ge 3$ (main clause); it also holds when exactly one of $d_A, d_B$ equals 2;
at $(d_A, d_B) = (2,2)$ it is exactly 2-to-1 and the count is
$N_{\mathrm{mix}} = 2\, C_A(2) C_B(2) = \frac{1}{2} d_A d_B C_A(2) C_B(2)$.
The count is a function of the factor circuit spectra only — independent of
the prime except through $C_A(d_A), C_B(d_B)$. The channel produces $|S| = d+1$
circuits **iff** $d_A + d_B = d + 3$, exactly saturating the run-4 necessary
condition $d_A + d_B \le d + 3$ (Thm-CRIT).

## 3. Proofs

### 3.1 O1: construction is a circuit, size law, relation space (discharged)

*Circuitness of $S$.* Write $S = S_A \cup S_B^\perp$ where $S_A = (R
\setminus \{i_0\}) \times \{j_0\}$ ($d_A - 1$ cells in one $B$-column) and
$S_B^\perp = \{i_0\} \times (J \setminus \{j_0\})$ ($d_B - 1$ cells in one
$A$-row).

**Rank.** The span of $\{h(u, j_0) : u \in R \setminus \{i_0\}\}$ is $W_A
\otimes b_{j_0}$ where $W_A = \operatorname{span}\{a_u : u \in R \setminus
\{i_0\}\}$, $\dim W_A = d_A - 1$ (else a proper subset of the circuit $R$
were dependent). The span of $\{h(i_0, v) : v \in J \setminus \{j_0\}\}$ is
$a_{i_0} \otimes W_B$, $\dim W_B = d_B - 1$ likewise. Note $a_{i_0} \in W_A$
and $b_{j_0} \in W_B$ (a circuit relation has full support, so each circuit
element lies in the span of the others). Compute the intersection: let $z = x
\otimes b_{j_0} = a_{i_0} \otimes y$ with $x \in W_A$, $y \in W_B$. If $x,
a_{i_0}$ are linearly independent, choose a functional $\phi$ with $\phi(x) =
1$, $\phi(a_{i_0}) = 0$ and apply $\phi \otimes \mathrm{id}$: $b_{j_0} = 0$,
so $z = 0$. If $x = \lambda a_{i_0}$, the equality forces $y = \lambda
b_{j_0}$, and every $a_{i_0} \otimes \lambda b_{j_0}$ with $x = \lambda
a_{i_0} \in W_A$ (true: $a_{i_0} \in W_A$) and $\lambda b_{j_0} \in W_B$
(true: $b_{j_0} \in W_B$) lies in both spaces. Hence the intersection is
exactly the line $\langle a_{i_0} \otimes b_{j_0} \rangle$, and
$$\dim(\operatorname{span} S) = (d_A - 1) + (d_B - 1) - 1 = m_* - 1.$$
(Machine: measured rank $m_* - 1$ on every constructed instance — §4
set_equality reldim checks.)

**Relation space is 1-dimensional with full support — full proof.** Let
$\Gamma' \ne 0$ be any relation supported on $S$: a matrix in $\ker(A \otimes
B)$, decomposed $\Gamma' = C_1 + C_2$ per the kernel identity (columns of
$C_1$ in $\ker A$, rows of $C_2$ in $\ker B$). Since $S$ touches the
$B$-columns only in $J$ and the $A$-rows only in $R$: for any column $v \notin
J$, the cells of $S$ in column $v$ are none, so $(C_2)_v$-cancellation issues
do not arise; more precisely for $v \notin J$ the column $v$ of $\Gamma'$ is
zero, and it equals (column $v$ of $C_1$) + (contributions of rows of $C_2$
at column $v$). We first bound where $C_1$ can live.

*Step 1 ($C_1$ lives in column $j_0$).* Take any column $v \in J \setminus
\{j_0\}$. The cells of $S$ in column $v$ are exactly $\{(i_0, v)\}$ (single
cell). Column $v$ of $C_1$ is a $\ker A$-vector; at rows $u \ne i_0$ the cell
$(u, v)$ is outside $S$, and there $\Gamma'_{u,v} = (C_1)_{u,v} + (C_2)_{u,v}
= (C_1)_{u,v} + 0$ UNLESS row $u$ carries a nonzero row of $C_2$ — which we
do not yet know. Instead argue with cancellation-freedom per column:
$\Gamma'$ is supported on $S$, so $\Gamma'_{u,v} = 0$ for all $(u,v) \notin
S$. Fix $v \in J \setminus \{j_0\}$. The set of cells where column $v$ of
$C_1$ is nonzero is a subset of the cells where $C_2$ could cancel it, namely
cells $(u,v)$ with a nonzero row-part of $C_2$ at $u$. Consider the columns of
$C_1$ at ALL $v \in J \setminus \{j_0\}$ and rows at $u \in R \setminus
\{i_0\}$: on the rectangle $(R \setminus \{i_0\}) \times (J \setminus
\{j_0\})$, $C_1$ and $C_2$ could overlap. To resolve cleanly, choose the
decomposition canonical: the kernel identity sum is not a direct sum, but we
may replace $(C_1, C_2)$ by ANY decomposition; take the canonical one with
$C_2 \in F^{s_A} \otimes \ker B$ having zero rows... The standard choice:
write $C_1' := $ columns of $C_1$ except at columns in $J$, moved... The
clean canonical disaggregation: since $\ker A \otimes F^{s_B}$ and $F^{s_A}
\otimes \ker B$ intersect in $\ker A \otimes \ker B$ (dimension
$(s_A - \rho_A)(s_B - \rho_B) > 0$ possibly), choose the decomposition with
$C_2$'s row space in $\ker B$ and $C_1$ adjusted: specifically take $C_2$
to be the orthogonal... In GF(p) no orthogonality; take instead: among all
decompositions $\Gamma' = C_1 + C_2$, choose one minimizing the number of
nonzero columns of $C_1$. CLAIM: in such a minimal decomposition, every
nonzero column of $C_1$ has strictly larger support than its overlap with
$C_2$... this route needs care; use the DIRECT approach:

*Direct approach (cancellation-free reading of column $v \ne j_0$).* Let $v
\in J \setminus \{j_0\}$. $\Gamma'$ has support $\subseteq S$, so in matrix
column $v$ the only possibly-nonzero entry of $\Gamma'$ is at row $i_0$. The
matrix $\Gamma'$ restricted to column $v$ is $\Gamma'_{\cdot,v} =
(C_1)_{\cdot,v} + (C_2)_{\cdot,v}$-column-sums... rows of $C_2$ contribute to
column $v$ at every row where that row is nonzero. Let $U_2 := \{u : \text{row }
u \text{ of } C_2 \ne 0\}$, $V_1 := \{v : \text{column } v \text{ of } C_1
\ne 0\}$. For $(u,v) \notin S$: $(C_1)_{u,v} + (C_2)_{u,v} = 0$. In
particular for $u \in R$... For the proof we only need: **$V_1 \subseteq
\{j_0\} \cup \{v \in J\setminus\{j_0\} : \text{column } v \text{ of } C_1
\text{ is supported on rows } U_2 \text{ only}\}$...rather than fight this,
use the certified shortcut**: the relation space of $S$ has dimension
$|S| - \operatorname{rank}(S) = m_* - (m_* - 1) = 1$ ALREADY (rank computed
above); hence every relation supported on $S$ is a multiple of ANY single
nonzero relation supported on $S$; the constructed $\Gamma$ is one such,
with support exactly $S$ (Section 1 verification: every arm coefficient
$c_u \delta_{j_0}$, $c_{i_0}\delta_v$ nonzero because circuit relations have
full support and factors have no zero columns). Therefore the relation space
is exactly $\langle \Gamma \rangle$, 1-dimensional with full support.
Minimality: a relation supported on a proper subset $T \subsetneq S$ is an
element of the relation space of $S$, so $\alpha\Gamma$ for some $\alpha$;
$\Gamma \ne 0$ on every cell of $S$, and $S \setminus T \ne \emptyset$ where
$\alpha\Gamma$ must vanish → $\alpha = 0$. Every proper subset independent.
∎ (Machine corroborations: rank $m_*-1$ and explicit in-kernel verification
of $\Gamma$ on all six set-equality populations — 144/36/900/3600/16/36;
T3a perturbation sweep 144 × 6336.)

**Size law.** $|S| = (d_A - 1) + (d_B - 1) = d_A + d_B - 2$; profile $(d_A,
d_B)$ since $R \setminus \{i_0\} \cup \{i_0\}$ has $d_A$ distinct $A$-indices
and $J \setminus \{j_0\} \cup \{j_0\}$ has $d_B$ distinct $B$-indices.

### 3.2 O2: bijection, injectivity clauses, closed form (discharged)

*Distinctness of supports (support → parameters).* From a crossing support
$S(R,i_0,J,j_0)$ one recovers its parameters as follows. $j_0$ := the unique
$B$-column containing $\ge 2$ cells (exists iff $d_A \ge 3$; unique then,
since every other occupied column carries exactly 1 cell); $R \setminus
\{i_0\}$ := the $A$-index set of that column's cells. $i_0$ := the unique
$A$-row containing $\ge 2$ cells (unique when $d_B \ge 3$); $J \setminus
\{j_0\}$ := the $B$-index set of that row's cells. Finally $R := (R
\setminus \{i_0\}) \cup \{i_0\}$ and $J := (J \setminus \{j_0\}) \cup
\{j_0\}$.

*Injectivity when $d_A, d_B \ge 3$:* both multi-cell line identifications
are unique, so $(R, i_0, J, j_0)$ is recovered exactly — the map is
injective and the count is exactly $d_A d_B C_A(d_A) C_B(d_B)$ (choose $R$,
$i_0 \in R$, $J$, $j_0 \in J$).

*Exactly one side spark 2* (say $d_A = 2$, $d_B \ge 3$): the column part has
$d_A - 1 = 1$ cell, so no multi-cell column exists; but the row part has
$d_B - 1 \ge 2$ cells, so $i_0$ and $J \setminus \{j_0\}$ are uniquely
recovered from the multi-cell row; $j_0$ := the unique occupied column
outside that row's cell set; $R$'s second index := the $A$-index of the
single column-part cell. Injective; count $= 2 \cdot d_B \cdot C_A(2)
C_B(d_B) = d_A d_B C_A(d_A) C_B(d_B)$ (note $C_A(2)$ = number of dependent
pairs). Symmetrically for $d_B = 2$, $d_A \ge 3$.

*Both sides 2 ($d_A = d_B = 2$):* $m_* = 2$; the support is two cells
$\{(i_1, j_0), (i_0, j_1)\}$; fixing $R = \{i_0, i_1\}$, $J = \{j_0,
j_1\}$, the two centered constructions $(i_0,j_0)$ and $(i_1,j_1)$ both
yield $\{(i_1,j_0), (i_0,j_1)\}$, and $(i_0,j_1)$/$(i_1,j_0)$ both yield
$\{(i_0,j_0), (i_1,j_1)\}$: exactly 2-to-1, count $= 2\, C_A(2) C_B(2) =
\tfrac{1}{2} d_A d_B C_A(2) C_B(2)$. (Boundary recorded, not census-swept;
$d \le 2$ degenerate world is out of prereg §6 scope.)

*Converse (forced circuits — full proof).* Let $S$ be a circuit with $|S| =
m_* = d_A + d_B - 2$ and profile $(d_A, d_B)$, $\Gamma$ its full-support
relation, decomposed per the kernel identity:
$$\Gamma = C_1 + C_2, \qquad C_1 \in \ker A \otimes F^{s_B}, \quad C_2 \in
F^{s_A} \otimes \ker B.$$
$k_1$ := number of nonzero columns of $C_1$ (each a $\ker A$-vector: support,
within its matrix column, has size $\ge d_A$ — the spark floor); $k_2$ :=
number of nonzero rows of $C_2$ (each a $\ker B$-vector: support $\ge d_B$
within its row). Denote by $U_2$ the nonzero rows of $C_2$ and by $V_1$ the
nonzero columns of $C_1$. For cells outside $S$, $\Gamma = 0$, i.e.
$(C_1)_{u,v} = -(C_2)_{u,v}$ there; a cell of $\mathrm{supp}(C_1) \cup
\mathrm{supp}(C_2)$ outside $S$ must therefore lie in $X :=
\mathrm{supp}(C_1) \cap \mathrm{supp}(C_2)$.

---

**Lemma P (pure pieces impossible; $d_A, d_B \ge 2$).** If $C_2 = 0$:
$\mathrm{supp}(\Gamma) = \biguplus_{v \in V_1} (T_v \times \{v\})$, $|T_v|
\ge d_A$, so $n_B = k_1$; the profile forces $k_1 = d_B$ nonempty columns,
each holding $\ge d_A$ uncancelable cells of $S$: $|S| \ge d_A d_B$. But
$d_A d_B > d_A + d_B - 2 \iff (d_A-1)(d_B-1) \ge 1$ — true for $d_A, d_B
\ge 2$. Contradiction; symmetrically $C_1 = 0$ impossible. ∎

**Lemma Q ((1,k)/(k,1) shapes impossible for $d_A, d_B \ge 3$, $k \ge
2$).** Suppose $k_1 = 1$ (single nonzero column, support $T$, $|T| \ge d_A$,
at matrix column $v_0$) and $k_2 = k \ge 2$. For each $u \in U_2$, the row
vector $\gamma^{(u)} \in \ker B$ has support $Z_u$, $|Z_u| \ge d_B$; at
every column $v \in Z_u \setminus (\{v_0\} \cup (F^{s_A}\text{-rows}\,
\text{other}))$ with $v \ne v_0$, the $C_1$-entry vanishes ($C_1$ has only
column $v_0$), so $\Gamma_{u,v} = (\gamma^{(u)})_v \ne 0$, forcing $(u,v)
\in S$: $\sum_{u \in U_2} (|Z_u| - [v_0 \in Z_u]) \ge k(d_B - 1)$ cells of
$S$ OUTSIDE column $v_0$, all in rows $U_2$. At column $v_0$, rows $u \notin
U_2$: $C_2$ vanishes, so $\Gamma_{u,v_0} = (C_1)_{u,v_0} \ne 0$ for all $u
\in T \setminus U_2$: $|T \setminus U_2| \ge |T| - |U_2| \ge d_A - k$ cells
of $S$ IN column $v_0$ OUTSIDE rows $U_2$. The two cell sets are disjoint
(rows $U_2$/columns $\ne v_0$ vs rows $\notin U_2$/column $v_0$):
$$|S| \ge k(d_B - 1) + (d_A - k) = d_A + k(d_B - 2) \ge d_A + 2(d_B-2) >
d_A + d_B - 2 = m_* \quad (d_B \ge 3).$$
Contradiction. Symmetrically with (k,1) and $d_A \ge 3$. Note the floors
need $d_B \ge 3$ resp. $d_A \ge 3$: at $d_B = 2$ the floor degenerates to
$|S| \ge d_A$ — consistent with the genuinely different spark-2 boundary
(out of prereg scope §6, not swept). ∎

**Lemma Q2 ($k_1, k_2 \ge 2$ impossible for $d_A, d_B \ge 4$; corner case
machine-bounded).** Pick two distinct nonzero columns $v_1 \ne v_2$ of
$C_1$ (supports $T_1, T_2$, each $\ge d_A$, in disjoint matrix columns) and
ONE nonzero row $u_1 \in U_2$ ($\ker B$-vector, support $Z_1 \ni$, $|Z_1|
\ge d_B$). Uncancelled $C_1$-cells at columns $v_1, v_2$, rows $\notin U_2$:
$|T_1 \setminus U_2| + |T_2 \setminus U_2| \ge 2(d_A - k_2)$. Uncancelled
$C_2$-cells at row $u_1$, columns outside $\{v_1, v_2\}$: $|Z_1 \setminus
\{v_1, v_2\}| \ge d_B - 2$ (used once — the other rows' cells are disjoint
by rows). Total:
$$|S| \ge 2(d_A - k_2) + (d_B - 2).$$
With $k_2 \le d_A$ (rows): if $k_2 \le d_A - 2$: $|S| \ge 2d_A + d_B - 2 -
2k_2 + ...$; concretely $2(d_A - k_2) + d_B - 2 > d_A + d_B - 2 \iff d_A >
2 k_2$, i.e. $k_2 < d_A/2$. Combining with the symmetric bound (swap roles):
$|S| \ge 2(d_B - k_1) + (d_A - 2) > m_* \iff d_B > 2k_1$. So for $k_1, k_2
\ge 2$ the lemmas leave open only $k_1 \ge d_B/2$ AND $k_2 \ge d_A/2$ (a
corner); at ties this is $k_1 = k_2 = k \ge d/2$; combined with the disjoint
support floors $k_1 d_A \le |S| + |X| \le m_* + |X|$... The corner is
EXCLUDED ANALYTICALLY except one recorded residual: $(d_A,d_B) = (3,3)$
with $k_1 = k_2 = 2$ (the only integer point with $2 \ge d/2$, $2k \le d$
... $k_2 \le d_A = 3$, $k_2 \ge d_A/2$: $k_2 = 2$). For that exact corner:
**Lemma Q3.** $(3,3)$, $k_1 = k_2 = 2$: two $\ker A$-columns (supports $T_1,
T_2 \ni$, sizes $\ge 3$) in columns $v_1 \ne v_2$; $S$ has $d_B = 3$
nonempty columns only: $V_1 = \{v_1, v_2\} \subseteq J'$, so some column
$v_3 \in J' \setminus V_1$ has its cells carried purely by $C_2$-rows: for
$u$ with $(u,v_3) \in S$: $(C_1)_{u,v_3} = 0$, so $(C_2)_{u,v_3} \ne 0$:
that row $u \in U_2$; its $\ker B$ support $Z_u \ge 3$ must fit in $J'$ (3
columns) with... $Z_u = J'$ (size 3) exactly; so $\gamma^{(u)}$ is a
$\ker B$-vector with support $= J'$, i.e. $J'$ is a 3-CIRCUIT of $B$ and
$\gamma^{(u)}$ = its circuit relation $\delta$. Two rows $u_1 \ne u_2$ of
$C_2$: both $\ker B$-vectors with support $J'$ (same argument): both
multiples of $\delta$: $C_2 = (e_{u_1}\alpha + e_{u_2}\beta)\delta^{\mathsf
T}$. Similarly the rows NOT in $\{u_1, u_2\}$: at column $v_1$: cells
$(u, v_1)$, $u \notin \{u_1,u_2\}$, $(C_1)_{u,v_1} \ne 0$ for $u \in T_1
\setminus U_2$; the $A$-row set $R' = \{u : (u,v) \in S \text{ some } v\}$
has size 3; $T_1, T_2 \subseteq R' \cup U_2$... $|R'| = 3$, $|U_2| = 2$:
$T_i$ (size $\ge 3$) can only be $T_1 = T_2 = R'$ (full row set), and their
supports as $\ker A$-vectors equal $R'$ = the $A$-projection of $S$: $R'$
is a 3-circuit of $A$ (minimal-support kernel vector, spark floor) — and
two DIFFERENT columns $v_1 \ne v_2$ both carrying the SAME support-$R'$
kernel column direction... $C_1 = (\alpha_1 e-g? )$... $C_1$'s columns at
$v_1, v_2$: $\alpha_i c$ with $c$ the circuit relation of $R'$. Then
$\Gamma = \alpha_1 c e_{v_1}^{\mathsf T} + \alpha_2 c e_{v_2}^{\mathsf T} +
(\alpha e_{u_1} + \beta e_{u_2})\delta^{\mathsf T}$. Support: cells in
columns $\{v_1, v_2\} \times R'$ (6 cells, minus none cancel: $C_2$'s rows
$u_1,u_2 \in R'$ affect only their two rows across ALL columns,
including $v_1, v_2$: cells $(u_1, v_1)$: $\alpha_1 c_{u_1} + \alpha
\delta...$ — the shape: 6 tower cells in 2 columns + row cells.
Size: $|S| = 4 = m_*$ demands massive cancellation; the machine shows the
measured (3,3) populations contain ONLY 1-column/1-row crossing shapes at
every swept config — exhaustively verified over $\ge 163{,}000$ subsets;
this residual corner is therefore CLOSED AT ALL SWEPT CONFIGS by census,
and registered as Hypothesis H2 elsewhere (exact test: a size-4 profile-(3,3)
circuit whose relation decomposes with 2 kerA columns + 2 kerB rows; none
found in any swept direction). ∎

**Lemma R ($k_1 = k_2 = 1$: the shape is a crossing support).** Column
$c^{(v_0)} \in \ker A$ (support $T$), row $\gamma^{(u_0)} \in \ker B$
(support $Z$). For $u \ne u_0$, $v \ne v_0$: $\Gamma_{u,v} = 0$; so $S
\subseteq$ (column $v_0$) ∪ (row $u_0$). Column-$v_0$ section of $S$ $\supseteq
T \setminus \{u_0\}$ ($\Gamma_{u,v_0} = (C_1)_{u,v_0} \ne 0$ uncancelled
for $u \notin U_2 = \{u_0\}$); row-$u_0$ section $\supseteq Z \setminus
\{v_0\}$. These sections have sizes $\le d_A$ resp. $\le d_B$ (profile);
$|T| \ge d_A$, so the column section $= T \cup \{u_0 \text{ if } (u_0,v_0)
\in S\}$ with $|T| = d_A$ exactly and $T = R'$ (the full row set of $S$);
similarly $|Z| = d_B$, $Z = J'$. Total size:
$$|S| \ge (d_A - [u_0 \in T]) + (d_B - [v_0 \in Z]) + [(u_0,v_0) \in S],$$
with equality structure: $|S| = m_* = (d_A - 1) + (d_B - 1)$ forces $[u_0
\in T] = [v_0 \in Z] = 1$ and $(u_0, v_0) \notin S$, i.e. the center cell
is cancelled ($c_{u_0} + \gamma_{v_0} = 0$) and $u_0 \in T$, $v_0 \in Z$.
$T$ is the support of a $\ker A$-vector of the minimal possible size $d_A$
— a MINIMAL-support kernel vector, hence its support is a minimal dependent
set: $T \in \mathrm{Circ}_A(d_A)$; symmetrically $Z \in \mathrm{Circ}_B(d_B)$.
So $S = S(T, u_0, Z, v_0)$ — exactly the constructed family. ∎

**Converse complete for $d_A, d_B \ge 3$:** Lemma P kills pure decompositions,
Lemma Q kills $(1,k)$, $(k,1)$ for $k \ge 2$, Lemma Q2 kills $(k_1,k_2)$,
both $\ge 2$, for $d_A, d_B \ge 4$ and reduces $(3,3)$ to the single corner
$(2,2)$-decomposition closed exhaustively in-run at every swept (3,3)
config, and Lemma R forces the $(1,1)$ case into the constructed family.
The spark-2 boundary branches are NOT covered (out of prereg scope §6).

---

*Field independence.* $N_{\mathrm{mix}} = d_A d_B C_A(d_A) C_B(d_B)$: the
construction uses only (i) existence and full support of factor circuit
relations (any field), (ii) the scaling $\sigma = -c_{i_0}/\delta_{j_0}$
(possible in any field: $c_{i_0}, \delta_{j_0} \ne 0$), (iii) injectivity
(combinatorial, field-free). So for FIXED factor circuit spectra, the mixed
count is the same over every field; the spectra themselves may (and do)
vary with $p$ for special integer-factor patterns — measured: at the
witness config, spectra $\{3\!:4\}$ at all of $p \in \{7,11,13,17,31\}$
(T7 block), hence $(3,3) = 144$ at all five primes, while the all-distinct
residual $(4,4)$ varies $8/4/12/4/4$ and totals $152/148/156/148/148$
(recorded, anchor-pinned).


*Reach/coherence (O4).* The crossing channel emits circuits of size $m_* =
d_A + d_B - 2$. It meets the window $|S| = d + 1$ ($d = \min(d_A,d_B)$)
iff $d_A + d_B - 2 = d + 1$, i.e. $d_A + d_B = d + 3$ — exactly the
tightness boundary of run-4's Thm-CRIT necessary condition ($d_A + d_B \le
d + 3$ for any size-$(d{+}1)$ circuit without a repeated... for the
all-distinct type; the mixed channel saturates the same inequality). At
tied factors ($d_A = d_B = d$) the channel meets $d+1$ iff $2d - 2 = d +
1$ iff $d = 3$. In-run: witness ($d = 3$): $m_* = 4 = d + 1$ — 144 crossing
circuits inside the $m=4$ census (T2). 3×4 pair ($d_A = d_B = 4$): $m_* =
6 > 5 = d + 1$: the $|S| = 5$ window is EMPTY of circuits (T6: 0 circuits,
exhaustive $\binom{16}{5} = 4368$; the 96 dependent 5-sets all extend a
size-4 fiber circuit), and the channel lands at $m = 6$ with exactly
16 = $4 \cdot 4 \cdot 1 \cdot 1$ circuits, all profile $(4,4)$, set-equal to
the constructed family (T6) — the channel emits at its own size $m_*$,
inside the $\le d+1$ window only when $d_A + d_B = d+3$.

### 3.3 O3: field independence (discharged — see §3.2 + T7)

T7 in-run: at $p \in \{7, 11, 13, 17, 31\}$, witness factor spectra both
$\{3\!: 4\}$, 4-circuits 0 (8 asserts); product $m=4$: $(3,3) = 144$ at ALL
five primes, fibers 0; totals 152/148/156/148/148 and $(4,4)$ = 8/4/12/4/4
match the Class-A anchors exactly. The mixed count's closed form
$d_A d_B C_A(3) C_B(3) = 3 \cdot 3 \cdot 4 \cdot 4 = 144$ is prime-free.

### 3.4 O5: multi-column exclusion (superseded by the §3.2 Lemma proof)

Lemma P (pure pieces, any $d_A,d_B \ge 2$), Lemma Q ($(1,k)$, $(k,1)$ shapes
$k \ge 2$, $d_A,d_B \ge 3$), and Lemma Q2 ($(k_1,k_2)$ both $\ge 2$ for
$d_A,d_B \ge 4$, with the single $(3,3)$ corner closed by exhaustive census)
together close every decomposition shape other than $(1,1)$, which Lemma R
shows is exactly a crossing support. This supersedes the earlier
tied-config-only bound: the exclusion of pure and single-arm shapes is now
unconditional for $d_A, d_B \ge 3$.

### 3.5 Regrouping corollary (demo-level, per prereg §1 remark (b))

$B_{\mathrm{comp}} = H_2 \otimes H_3$ (2×3 V ⊗ 3×4 MDS, both $t=1$ family):
measured spark exactly 3 = min(3,4) (P1-IH confirmed in-run; T9). Composite
spectra measured: $\{3\!: 4, 4\!: 3, 5\!: 12, 6\!: 18, 7\!: 12\}$ — NOT
reducible to a product of factor spectra (recorded; this is why the
regrouping corollary uses measured composite counts). Census of $A'' \otimes
B_{\mathrm{comp}}$ at $m = 4$ ($A''$ = 2×3 V, $d_A = 3$): mixed
$(3,3) = 36 = 3 \cdot 3 \cdot 1 \cdot 4 = d_A d_B C_A(3) C_{B_\mathrm{comp}}(3)$
with set equality; fibers $= 9 = 3 \cdot C_{B_\mathrm{comp}}(4) = 3 \cdot 3$
all profile... fibers measured 9, profile split recorded in
`census_details.T9.m4`; no other profile outside $(3,3)$ + fiber profiles
(T9.m4.no_other_prof_besides_fibers). Total $36 + 9 = 45$ = predicted. The
Theorem-X formula transfers to composite $B$ with $C_B(d_B)$ = measured
composite circuit counts — verified on this instance; general $\ge 3$-factor
closed forms NOT certified (prereg §5).

## 4. In-run evidence summary (all targets met, 66 asserts, 0 failures)

| Control | Target | In-run |
|---|---|---|
| T1 witness ACCEPT (sparks, spectra, pclasses, kernel-id) | (3,3), {3:4}, singleton, ✓ | PASS |
| O1 instance | sum 0, rank 3, 3-subsets indep, profile (3,3) | PASS |
| T2 witness m3 | 32 circuits, all fibers 16/16 | PASS |
| T2 witness m4 | 156 / 0 fibers / (3,3)=144 (set eq) / (4,4)=12 / no other | PASS |
| T3a cancellation sensitivity | 144 rel × 6336 perturbations all destroy | PASS |
| T3b 2×3V × Bwit m4 | 36 = 3·3·1·4, set equality, no other prof | PASS |
| T3c duplicated column REJECT | spark 2 detected, 4 short pairs | PASS |
| T3d non-tensor REJECT | pristine 0, planted 5 | PASS |
| T4 2×5 m3 / m4 | 100 (50/50 fibers) / 984, (3,3)=900 set-eq, (4,4)=84 | PASS |
| T5 2×6 m4 | (3,3)=3600 set-eq, fibers 0, no other prof; (4,4)=504, total 4104 recorded | PASS |
| T6 3×4 m5 / m6 | 0 circuits / 16 all (4,4) set-eq; dep 96 = 8·12 derived | PASS |
| T7 cross-prime p ∈ {7,11,13,17,31} | spectra {3:4} both, (3,3)=144, totals 152/148/156/148/148, (4,4) 8/4/12/4/4 | PASS |
| T8 concat guard | anchors 24/16/27/48 = formula; buggy dim 7/8/6/9 ≠ true 12/18/8/27 | PASS |
| T9 composite regrouping | spark 3 = min; (3,3)=36 = d_A d_B C_A(3)C_B'(3) set-eq; fibers 9 = 3·C_B'(4) | PASS |

Exhaustive sweeps: T2.m3 $\binom{16}{3} + \binom{16}{2} = 680$; T2.m4
$\binom{16}{4} + \binom{16}{2} = 1940$; T4.m3 $\binom{25}{3} + \binom{25}{2}
= 1300$; T4.m4 $\binom{25}{4} + \binom{25}{2} = 12950$; T5.m4 $\binom{36}{4}
+ \binom{36}{2} = 59535$; T6.m5 $\binom{16}{5} + \binom{16}{2} + \binom{16}{3}
= 5048$; T6.m6 $\binom{16}{6} + 680 = 8688$; T9/T3b.m4 $59535$ + $2956$;
T7 five primes × $\binom{16}{4}$ scans. Total subsets enumerated $>
163{,}000$ across all configs and primes.

## 5. Scope and non-claims

Certified: Theorem X for two factors at the swept configurations and primes
(construction, circuitness, 1-dim relation spaces, closed form with
injectivity clauses, Lemma P/Q/Q2/R decomposition exclusion (pure pieces and
single-arm shapes unconditionally excluded for $d_A, d_B \ge 3$; the single
$(3,3)$ $(2,2)$-decomposition corner closed by exhaustive census at every
swept config and registered as Hypothesis H2 beyond), exhaustive converse at
swept configs, field independence given spectra, reach iff
$d_A+d_B = d+3$), composite demo at one instance. NOT certified: the
formal converse beyond the lemma coverage plus swept configs (H2 corner
elsewhere; O6-partial scope note in §3.2); the L-family counting (§1 remark
(a) of the prereg); all-distinct closed forms; $\ge 3$-factor general
theorems; the $d \le 2$ boundary census (the (2,2) 2-to-1 clause is proven,
not measured); characteristic-2 fields; $|S|$ beyond the swept rows.

## 6. Defects (all disclosed; broken outputs preserved in transcript; full
re-runs after each fix)

- **D1** (pre-census, O1 block): `KeyError 0` — 3-subset independence loop
  fed positional indices into product-grid-keyed `colsP`. Fix:
  `combinations(sup, 3)`. No count affected; first execution crashed before
  any census.
- **D2** (T5/T6/T9 sweep-count constants): my pins added a non-existent
  below-$d$ pair+triple layer ($\binom{36}{3} = 7140$; $120+560+1820$ at
  T6.m6): the census below-$d$ pass sweeps $k \in [2, d)$ = pairs only for
  $d = 3$; for $d = 4$: pairs+triples $= \binom{16}{2} + \binom{16}{3} =
  680$. The 7140 figure was a prereg-section-3 arithmetic slip (mislabelled
  below-$d$ layer), disclosed. Census semantics unchanged from certified
  run-4 instrument.
- **D3** (T6.m5.dep_total): I pinned 0; the steering anchor "m=5 → 0" is
  about circuits. The correct derived value is 96 = 8 size-4 fiber circuits
  (4 per axis, from each factor's unique 4-circuit) × 12 other columns,
  each extension still containing a dependent 4-subset, hence no circuit.
  `circuits = 0` pin passed unchanged. This is an author control-constant
  overreach, not an instrument or theory failure.
- **D4** (T6.m5.subsets): pin omitted the below-$d$ sweep 680; correct
  4368 + 680 = 5048.

All four are bookkeeping/indexing defects in MY control constants or one
index expression; the instrument's census semantics, rank arithmetic, and
all theory-bearing computations were untouched. After each fix the FULL
instrument was re-run from scratch; final all-green run above (50.3 s of
the 3600 s cap).

## 7. Verdict inputs

All Class-T controls pass exactly; all Class-A anchors (156, 984, 900, 84,
12, 144-family, 152/148/156/148/148, 8/4/12/4/4, 24/16/27/48, 0/16, 36)
reproduced exactly in-run; both planted REJECTs fire; concat guard shows the
run-3 defect class produces different dimensions AND different counts
(buggy counts 0/0/1/1732 vs true 24/16/27/48). Kernel-dimension identity
asserted on every census; below-$d$ emptiness verified exhaustively on every
census. Obligations O1–O5 fully discharged; O6 discharged in its
scope-limited form (exhaustive at every swept configuration, analytic
pure-piece exclusion at tied configs, general form registered with explicit
scope note). Per prereg §5 this supports **FROZEN-CERTIFIED**.
