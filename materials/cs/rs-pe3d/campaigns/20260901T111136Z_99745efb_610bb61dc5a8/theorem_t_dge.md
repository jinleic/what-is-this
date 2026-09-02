# THEOREM T-DGE — complete proof, discharging P1-P5 — gate H-MINLINE-DGE3-PROOF

Run: `20260901T111136Z_99745efb_610bb61dc5a8`. This artifact discharges the
five obligations P1-P5 registered in the committed prereg
(`pre_statement.md`, sha256 `558ba6fe...d84c`, commit `f3d1b03`). Every step
is written out; the machine controls M1-M6 of the prereg (run separately in
this run dir) anchor the statement's finite content. Verdict language per
the prereg's fixed rule.

Notation. $F$ is an arbitrary field; $F^m$ means the $m$-dimensional
$F$-vector space of column vectors. For a matrix $M$ over $F$, $\mathrm{col}(M)$
is its column space; the **spark** $d(M)$ is the least size of a linearly
DEPENDENT set of columns ($= \infty$ if all columns are independent, which
cannot happen here when the row count is finite and columns exceed it, and
is irrelevant otherwise). A dependent set of columns is a **circuit** if all
its proper subsets are independent. "Parallel" means linearly dependent as
a 2-element set, i.e. one vector is a scalar multiple of the other; this
relation is written $u \parallel v$. For vectors $u_i \in F^{r_i}$ and
coordinates $i$, $\bigotimes_{i=1}^n u_i$ is the Kronecker product column,
an element of $F^{r_1 \cdots r_n}$; it is **simple** (a pure tensor) by
construction. $\mathrm{id}$ is the identity functional/map and for a
functional $\varphi \in (F^{r})^*$, $\varphi \otimes \mathrm{id}$ acts on
$F^{r} \otimes F^{s} \cong F^{rs}$ (row-major flattening: index
$(a,b) \mapsto a \cdot s + b$) by
$(\varphi \otimes \mathrm{id})(u \otimes v) = \varphi(u)\, v$ and linear
extension; the dual action on general tensors is the linear map induced by
$u \otimes v \mapsto \varphi(u) v$, which is well defined because
$u \otimes v \mapsto \varphi(u) v$ is bilinear.

Standing hypothesis throughout: each factor $H_i$ has **all columns
nonzero**, spark $d_i$, and $d := \min_i d_i \ge 2$. (Nonzero columns are
load-bearing: they make $d_i \ge 2$ and every tensor column nonzero.)

-----------------------------------------------------------------------

## P3 (discharged first — used by P1's converse step and P2)

**P3 NO-PARALLEL PRODUCT LEMMA.** Let $u_i, v_i \in F^{r_i}$ with
$u_i \ne 0$ and $v_i \ne 0$ for all $i$. Then
$\bigotimes_i u_i \parallel \bigotimes_i v_i$ iff $u_i \parallel v_i$ for
every $i$.

*Proof.* ($\Leftarrow$) If $v_i = \lambda_i u_i$ with all $\lambda_i \ne 0$
(none can be zero because $v_i \ne 0$, and a parallel pair with a zero
member would force the other to be zero), then
$\bigotimes_i v_i = (\prod_i \lambda_i) \bigotimes_i u_i$, and
$\prod_i \lambda_i \ne 0$ in a field. So the tensors are parallel.

($\Rightarrow$) Suppose $\bigotimes_i u_i = c \bigotimes_i v_i$ for some
scalar $c \ne 0$ (parallel nonzero tensors: both sides are nonzero —
Kronecker products of nonzero vectors are nonzero, since each coordinate
$(a_1,\dots,a_n)$ of the product is $\prod_i u_i(a_i) \ne 0$ — so neither
side is zero and the proportionality scalar is nonzero). Fix $i$. We must
show $u_i \parallel v_i$. For each $j \ne i$, choose a linear functional
$\varphi_j \in (F^{r_j})^*$ with $\varphi_j(v_j) \ne 0$; this is possible
over ANY field because $v_j \ne 0$ is a nonzero vector in a finite-
dimensional space: pick a basis of $F^{r_j}$ containing $v_j$ (any linearly
independent set extends to a basis, by the Steinitz exchange theorem, which
holds over every field) and let $\varphi_j$ read the $v_j$-coordinate,
i.e. the coordinate functional dual to the basis element $v_j$; then
$\varphi_j(v_j) = 1 \ne 0$. Now apply the linear map
$\bigl(\bigotimes_{j \ne i} \varphi_j\bigr) \otimes \mathrm{id}$ to both
sides: this map sends $\bigotimes_k w_k$ to
$\bigl(\prod_{j \ne i} \varphi_j(w_j)\bigr) w_i$ (bilinearity of the
extraction in each slot). We get
$$\Bigl(\prod_{j \ne i} \varphi_j(u_j)\Bigr) u_i \;=\; c \Bigl(\prod_{j
\ne i} \varphi_j(v_j)\Bigr) v_i.$$
The right-hand side is nonzero, because $c \ne 0$, every
$\varphi_j(v_j) \ne 0$ by choice, and $v_i \ne 0$ (nonzero scalar times
nonzero vector). Hence the left scalar $\prod_{j \ne i} \varphi_j(u_j)$ is
nonzero and $u_i = c' v_i$ with
$c' = c \prod_{j\ne i} \varphi_j(v_j) / \prod_{j \ne i} \varphi_j(u_j) \ne
0$, i.e. $u_i \parallel v_i$. Since $i$ was arbitrary, this holds for every
coordinate. $\blacksquare$

**Corollary P3'.** If no factor $H_i$ has two distinct parallel columns,
then no product column pair is parallel: for distinct tuples $a \ne b$,
$h(a) \parallel h(b)$ would force $h_i(a_i) \parallel h_i(b_i)$ for every
$i$ (P3, converse direction), and since $a \ne b$ some coordinate has
$a_i \ne b_i$, where $H_i$ has no such pair — contradiction. In
particular **every factor with $d_i \ge 3$ has no parallel pair of
distinct columns** (a parallel pair is a dependent 2-set), and so does
their product. Conversely, $d_i = 2$ means SOME 2-subset of columns is
dependent, i.e. some distinct parallel pair exists (a dependent 2-set of
nonzero vectors is exactly a parallel pair). So, given all columns
nonzero: "$d_i \ge 3$" $\iff$ "no parallel pair of distinct columns."
$\square$

-----------------------------------------------------------------------

## P1 (spark of the product)

**P1 THEOREM.** $\mathrm{spark}(H_1 \otimes \cdots \otimes H_n) = d =
\min_i d_i$.

*Proof, by induction on $n$.*

*Base $n = 1$:* trivial, $d = d_1$.

*Step.* Assume the theorem for $n - 1$ factors ($n \ge 2$). Regroup
$$H = A \otimes B, \qquad A := H_1, \quad B := H_2 \otimes \cdots \otimes
H_n.$$
$B$ is a legitimate factor for the induction: (i) all its columns are
nonzero — its columns are $\bigotimes_{i\ge2} h_i(a_i)$, each a product of
nonzero vectors, hence nonzero in every coordinate; (ii) by the induction
hypothesis $\mathrm{spark}(B) = \min_{i \ge 2} d_i =: d_B$, and $d_A :=
d_1$, so $d = \min(d_A, d_B)$.

We show $\mathrm{spark}(A \otimes B) = \min(d_A, d_B)$.

**(1) $\mathrm{spark}(A \otimes B) \ge \min(d_A, d_B)$: no dependent set of
size $k < \min(d_A, d_B)$.** Suppose
$\sum_{(u,v) \in S} c_{uv}\, a_u \otimes b_v = 0$ with $|S| = k <
\min(d_A, d_B)$ and not all $c_{uv}$ zero; discard zero-coefficient
terms, so WLOG all $c_{uv} \ne 0$. Let $U \subseteq [s_A]$ be the set of
DISTINCT first indices, $r := |U| \le k < \min(d_A, d_B) \le d_A$. Hence
the distinct columns $\{a_u : u \in U\}$ are linearly independent (any
dependent subfamily would be a dependent set of $\le r \le k < d_A$
columns). Choose $u_0 \in U$ and, by independence, a functional $\varphi$
on $F^{r_A}$ with $\varphi(a_{u_0}) = 1$ and $\varphi(a_u) = 0$ for all
$u \in U \setminus \{u_0\}$ (basis extension as in P3: the independent
list $a_{u_0}, a_{u_1}, \ldots$ extends to a basis; $\varphi$ reads off
the $a_{u_0}$-coordinate). Apply $\varphi \otimes \mathrm{id}$:
$$\sum_{v \in B_{u_0}} c_{u_0 v}\, b_v \;=\; 0, \qquad B_{u_0} := \{v :
(u_0, v) \in S\},$$
where all terms with $u \ne u_0$ vanish because $\varphi(a_u) = 0$. The
indices $v \in B_{u_0}$ are DISTINCT (distinct elements $(u_0, v) \ne
(u_0, v')$ of $S$ have $v \ne v'$), their count satisfies $|B_{u_0}| \le
|S| = k < \min(d_A, d_B) \le d_B$, every $c_{u_0 v} \ne 0$ and every
$b_v \ne 0$. The dependence is NONTRIVIAL with $\ge 2$ terms: if
$|B_{u_0}| = 1$ then $c\, b_v = 0$ with $b_v \ne 0$ forces $c = 0$,
excluded. So $B_{u_0}$ is a nontrivially dependent set of distinct
nonzero columns of $B$ of size $|B_{u_0}| \le k < d_B$, contradicting
$\mathrm{spark}(B) = d_B$. Hence no such $S$ exists, and
$\mathrm{spark}(A \otimes B) \ge \min(d_A, d_B) = d$.

**(2) $\mathrm{spark}(A \otimes B) \le d$.** Take $U_A$ a dependent set of
$d_A$ columns of $A$ and fix any $v_0$ (exists, $B$ has at least one
column, all nonzero); then $\{a_u \otimes b_{v_0}\}$ is dependent
(hyperplane argument below), unless $d_A$ is the min, in which case we are
done; if instead $d_B < d_A$, take $V_B$ dependent of size $d_B$ in $B$
and fix $u_0$: $\{a_{u_0} \otimes b_v\}_{v \in V_B}$ is dependent by the
same argument. In either case a dependent set of size
$\min(d_A, d_B) = d$ exists, so $\mathrm{spark}(A \otimes B) \le d$.
*Hyperplane argument.* For a dependent set $\{x_1, \dots, x_m\}$ in a
vector space, $y \in V$ nonzero implies $\{x_i \otimes y\}$ is dependent:
if $\sum_i \gamma_i x_i = 0$ is a nontrivial dependence, then
$\sum_i \gamma_i (x_i \otimes y) = \bigl(\sum_i \gamma_i x_i\bigr) \otimes
y = 0 \otimes y = 0$ is a nontrivial dependence (same coefficients, which
are not all zero). The members $x_i \otimes y$ are DISTINCT: the map $x
\mapsto x \otimes y$ is injective for $y \ne 0$, since for any functional
$\psi$ with $\psi(y) = 1$ (basis extension), the linear map
$\mathrm{id} \otimes \psi$ satisfies $(\mathrm{id} \otimes \psi)(x \otimes
y) = \psi(y)\, x = x$, recovering $x$ from $x \otimes y$. A nontrivial
dependence among distinct vectors transports, coefficient for coefficient,
to a nontrivial dependence among their images.

By (1) and (2), $\mathrm{spark}(A \otimes B) = d$, closing the induction.
$\blacksquare$

*Remark on the $r = d$ edge:* in step (1) we only ever have $r \le k < d
\le d_A$ STRICTLY, so no circuit/degenerate case arises below $d$; the
size-exactly-$d$ analysis is P2's business. The base case $n = 1$ needs no
isolation. Regrouping for $n \ge 3$ uses the canonical associativity
isomorphism $(F^{r_1} \otimes F^{r_2}) \otimes F^{r_3} \cong F^{r_1}
\otimes (F^{r_2} \otimes F^{r_3})$, under which all arguments above are
basis-free; the row-major flattening pinned in the instrument is pure
index bookkeeping and carries no mathematical content.

-----------------------------------------------------------------------

## P2 (classification of size-$d$ dependent sets)

By P1, spark$(H) = d$, so any dependent set of exactly $d$ columns is a
circuit of $H$ (a proper dependent subset would have size $< d$,
impossible). Consequently all coefficients in any dependence among its
members are nonzero, and the space of dependences on the $d$ involved
columns is 1-DIMENSIONAL: if two linearly independent dependence
coefficient vectors existed, a nontrivial linear combination of them
would vanish at one member while remaining a dependence, producing a
dependence of a proper subset — contradiction. Every nonzero dependence
vector on a circuit is therefore everywhere nonzero (any zero coordinate
would itself be a dependence of a proper subset).

Let $S = \{(u_1, v_1), \dots, (u_d, v_d)\}$ (2-factor notation, $A$-index
$u$, $B$-index $v$; for $n \ge 3$ regroup as in P1) be dependent:
$$\sum_{k=1}^{d} c_k\, a_{u_k} \otimes b_{v_k} = 0, \qquad c_k \ne 0
\ \forall k.$$
Let $r = |\{u_1, \dots, u_d\}|$ (distinct A-indices). Two cases.
**Case 1: $r < d_A$.** The distinct columns $\{a_u : u \in U\}$ are
independent ($|U| = r < d_A$). Dual-isolate: for each $u_0 \in U$ pick a
functional $\varphi$ with $\varphi(a_{u_0}) = 1$ and $\varphi(a_u) = 0$
for $u \in U$, $u \ne u_0$ (basis extension). Applying $\varphi \otimes
\mathrm{id}$ collapses the dependence into the $u_0$-group:
$\sum_{v \in G_{u_0}} c_{u_0 v}\, b_v = 0$ where $G_{u_0} := \{v :
(u_0, v) \in S\}$, a NONTRIVIAL dependence (all $c$'s nonzero, all
$b$'s nonzero) among the distinct columns $\{b_v : v \in G_{u_0}\}$.
A nontrivial dependence among distinct columns needs at least $d_B$
columns, and $d_B \ge d$; so $|G_{u_0}| \ge d_B \ge d$ for EVERY $u_0 \in
U$. The groups $G_{u_0}$ partition $S$ ($\sum_{u_0} |G_{u_0}| = d$), so
there is exactly ONE group, containing all $d$ elements, and $d_B = d$.
All $u_k$ are equal — $S$ agrees on the A-axis — and the $v$-indices form
a dependent set of exactly $d = d_B$ distinct B-columns: a $B$-circuit
(minimal: any proper subset would be a dependent set of $< d$ B-columns).
This is an axis fiber along $B$ (A-side fixed) whose varying side is a
size-$d$ circuit of $B$.

**Case 2: $r = d$** (all A-indices distinct; note $r \le d \le d_A$
always, and $d = d_A$ will be forced inside the case). The tuples are
$(u_k, v_k)$ with the $u_k$ pairwise distinct. For EVERY functional $\psi
\in (F^{r_B})^*$, apply $\mathrm{id} \otimes \psi$ to the dependence:
$$\sum_{k=1}^{d} c_k\, \psi(b_{v_k})\, a_{u_k} = 0 \qquad (\*)$$
a dependence among the $r = d$ distinct columns $\{a_{u_k}\}$. Whether
these are independent depends on $d$ versus $d_A$:

*Subcase 2a: $d < d_A$.* Then $\{a_{u_k}\}$ is independent (any dependent
$d$-set would witness spark $\le d < d_A$, contradiction). So in $(*)$
all coefficients vanish: $c_k \psi(b_{v_k}) = 0$ for every $\psi$ and
every $k$. Since $c_k \ne 0$: $\psi(b_{v_k}) = 0$ for all $\psi$; the
functionals separate points of the finite-dimensional space $F^{r_B}$
(basis-coordinate functionals: if $w \ne 0$, the functional dual to a
basis element containing $w$ reads a nonzero coordinate — over any field,
by basis extension), so $b_{v_k} = 0$ — contradicting all columns nonzero.
Hence Subcase 2a is impossible.

*Subcase 2b: $d = d_A$.* If $\{a_{u_k}\}$ were independent, the same
argument as 2a applies ($c_k \psi(b_{v_k}) = 0$ for all $\psi$ gives $b=0$,
contradiction), so $\{a_{u_k}\}$ is DEPENDENT; being $d = d_A$-sized in a
matrix of spark $d_A$, it is a circuit with 1-dimensional dependence space
spanned by some $\alpha \in F^d$ with ALL $\alpha_k \ne 0$ (everywhere
nonzero: if $\alpha_j = 0$ then dropping $u_j$ leaves a dependent set of
size $d-1$, contradicting circuit minimality). Every dependence on
$\{a_{u_k}\}$ is a multiple of $\alpha$. Now $(*)$ says exactly that
$\bigl(c_k \psi(b_{v_k})\bigr)_k$ is a dependence on the circuit, so
$$c_k\, \psi(b_{v_k}) = \lambda(\psi)\, \alpha_k \quad \text{for all } k,
\text{ for some scalar } \lambda(\psi) \text{ depending on } \psi.$$
**No zero-coefficient case split is needed:** $\lambda(\psi) = 0$ is
simply the member $0 \in \mathrm{span}(\alpha)$; the conclusion that all
$k$ share one value holds uniformly. Divide by $c_k \alpha_k \ne 0$
(field):
$$\psi\!\left( \frac{c_k b_{v_k}}{\alpha_k} \right) = \lambda(\psi)
\qquad \text{for every } k \text{ and every } \psi.$$
For fixed $\psi$ the value $\lambda(\psi)$ is the same for all $k$; hence
for every $j, k$ and every $\psi$:
$\psi\bigl( c_j b_{v_j}/\alpha_j - c_k b_{v_k}/\alpha_k \bigr) = 0$.
Functionals separate points, so $c_j b_{v_j}/\alpha_j = c_k b_{v_k}/
\alpha_k$ for all $j,k$: all the vectors $c_k b_{v_k}/\alpha_k$ equal a
SINGLE vector $w$. The common value is nonzero: $c_k \ne 0$ (circuit, all
coefficients nonzero) and $b_{v_k} \ne 0$, so $w = c_k b_{v_k}/\alpha_k
\ne 0$. Therefore $b_{v_j} = w \alpha_j / c_j$: **all $b_{v_j}$ are
parallel** (multiples of $w$; or, if one prefers, of $b_{v_1}$, since
$b_{v_1} \ne 0$).

Now invoke the parallel structure of $B$:
- If $B$ has NO two distinct parallel columns (in particular whenever
  $d_B \ge 3$, by P3'), then $b_{v_j} \parallel b_{v_k}$ with $b$'s
  nonzero forces $v_j = v_k$. All B-indices coincide: $S$ is an axis fiber
  along $A$ whose varying A-columns form a size-$d$ circuit of $A$ (the
  fiber set is $\{(u_k, v_0)\}_{k}$ — distinct tuples, agreeing on the
  B-axis).
- If $B$ DOES have distinct parallel columns ($d_B = 2$), Case 2 can end
  in a "diagonal": the tuple set keeps $u_k$ distinct and has $v_j \ne
  v_k$ parallel partners. This is exactly the regime-2 phenomenon; it
  cannot occur when $d = \min(d_A, d_B) \ge 3$. When exactly one factor
  has $d_i = 2$ (the unique degenerate axis, $d = 2$), a size-2 dependent
  set cannot land in Case 2 at all: $r = d = 2$ distinct A-indices with a
  dependent A-pair would force $d_A = 2$ (subcase 2b), contradicting the
  uniqueness of the degenerate axis, and subcase 2a is impossible anyway;
  so Case 1 governs every dependent pair — see P4's regime lemma for the
  direct P3 route. This completes the structural dichotomy; the regime
  bookkeeping is finished in P4.

**Summary of P2:** every size-$d$ dependent set of product columns is an
axis fiber over a size-$d$ circuit of one factor, EXCEPT when the
other-side columns include a distinct parallel pair — which requires
$d_B = 2$ on that other side. $\blacksquare$

*Audit points re-verified:* (a) basis extension over an arbitrary field
gives the isolation/separation functionals — Steinitz exchange holds over
every field (it uses only the field axioms and finiteness of dimension);
(b) the circuit kernel is 1-dim and everywhere nonzero — a zero coordinate
would yield a smaller dependent subset, contradicting minimality;
(c) $0 \in \mathrm{span}(\alpha)$ makes the $\lambda(\psi) = 0$ branch a
non-case: dividing by $\alpha_k \ne 0$ is legitimate in both branches and
the all-equal conclusion is uniform, with nonzeroness of $w$ coming
separately from $c_k, b_{v_k} \ne 0$; (d) functionals separate points in
finite dimension over any field — the coordinate functionals of a basis
are defined for every field (no ordering, no topology used).

-----------------------------------------------------------------------

## P4 (regime classification + closed-form counts)

**P4 REGIME LEMMA.** For a 2-factor product $A \otimes B$ with $d =
\min(d_A, d_B)$:
(i) if $d_A, d_B \ge 3$ (equiv. no factor has a distinct parallel pair):
every size-$d$ dependent set is a fiber (P2 Case 1 or Case 2b with the
no-parallel branch; both directions covered — Case 1 gives a fiber along
B, Case 2b a fiber along A);
(ii) if $d_A = 2 < d_B$ (exactly one degenerate factor, $d = 2$): every
dependent pair $h(a) = a_u \otimes b_v,\ h(a') = a_{u'} \otimes b_{v'}$
forces each coordinate pair parallel — $a_u \parallel a_{u'}$ and $b_v
\parallel b_{v'}$ — because a dependent pair of product columns means one
product column is a scalar multiple of the other, which is precisely P3's
converse hypothesis. Since $d_B \ge 3$, distinct $b_v \parallel b_{v'}$ is
impossible, so $v = v'$: the pair agrees off the degenerate axis, with
$u \ne u'$ (the tuples are distinct) forming a dependent pair of
A-columns, i.e. an axis-A fiber whose varying pair is a 2-circuit of $A$.
Conversely every such fiber is a dependent pair (the A-pair dependence
transported by the hyperplane argument of P1(2)). If additionally $d_A =
d_B = 2$ (both factors degenerate): pairs agreeing up to parallel class on
BOTH axes with distinct tuples on both are dependent (P3 converse) but NOT
fibers: these are the diagonals.
(iii) For $n \ge 3$ the classification follows by regrouping (P1's
associativity induction): a size-$d$ dependent set that is not a fiber
requires, by P2, a distinct parallel pair on some axis outside the fiber
witness — which forces that axis to have $d_i = 2$; with at most one
spark-2 factor this cannot occur (P3' + P2). Completing (B): if $d \ge 3$
then no factor has spark 2 (each $d_i \ge 3$), and P2 closes both
groupings; if $d = 2$ with exactly one factor of spark 2, P3 alone
finishes, as in (ii).

**P4 COUNTS.**
(1) Dependent-pair count at $d = 2$ (any regime): for factor $i$, let the
parallel relation (proportionality) partition the column set $[s_i]$ into
maximal classes of sizes $m_{i,1}, \dots, m_{i,k_i}$ (singletons included —
parallelism is an equivalence relation on nonzero columns), and set $A_i
:= \sum_k m_{i,k}^2$; let $N := \prod_i s_i$. By P3's converse, a pair of
DISTINCT tuples $\{p, q\}$ is a dependent pair iff for EVERY $i$ the
factor columns are parallel, i.e. $q_i$ lies in $p_i$'s parallel class
for each $i$: such $\{p, q\}$ are exactly the class-aligned distinct
tuples. Counting ordered class-aligned pairs: $\prod_i A_i$; of these $N$
are the self-pairs $p = q$; every remaining one is a distinct pair counted
twice (once in each order), so
$$\#\{\text{dependent unordered distinct pairs}\} = \frac{\prod_i A_i -
N}{2}.$$
Note the count is $0$ exactly when $\prod_i A_i = N$, i.e. every class of
every factor is a singleton (no factor has a parallel pair): the
parallel-free regime, consistent with regime 1 (no dependent pairs at
all).
(2) Fiber subcount at $d = 2$: a fiber pair DISAGREES only on one axis
$j$ (all other coordinates are equal) with the varying pair inside one
parallel class of factor $j$:
$$\#\text{fibers} = \sum_j \Bigl(\prod_{i \ne j} s_i\Bigr) \sum_k
\binom{m_{j,k}}{2},$$
and the remainder $\frac{\prod A_i - N}{2} - \#\text{fibers}$ are the
off-fiber (diagonal) pairs. These are consistent by the combinatorial
identity $A_i = s_i + \sum_k (m_{i,k}^2 - m_{i,k}) = s_i + 2\sum_k
\binom{m_{i,k}}{2}$: expanding $\prod_i A_i = \prod_i \bigl(s_i + 2\sum_k
\binom{m_{i,k}}{2}\bigr)$, each factor chooses per axis between "equal"
($s_i$ ways) and "distinct-within-class" ($2\binom{m}{2}$ ordered ways);
subtracting the all-"equal" term $N$ and halving maps each axis-choice
combination to exactly one unordered pair, with at least one
"distinct-within-class" choice: fibers are those with exactly one such
axis, diagonals those with several.
(3) Fiber-regime carrier count ($d \ge 2$, at most one $d_i = 2$):
$$\#\text{carriers} = \sum_{i : d_i = d} \Bigl(\prod_{j \ne i}
s_j\Bigr) \cdot \#\{\text{size-}d\text{ circuits of } H_i\}.$$
Only axes with $d_i = d$ contribute: an axis fiber over factor $i$ needs a
size-$d$ dependent set of $H_i$, so $d_i \le d$, i.e. $d_i = d$.
Conversely, each size-$d$ circuit of $H_i$ combined with each fixed tuple
on the other axes gives a dependent $d$-set (P1's hyperplane argument).
Distinct (axis, circuit, other-tuple) triples give distinct column sets,
and no carrier is counted twice: for $n \ge 2$ a $d$-set ($d \ge 2$,
distinct tuples) that agrees off axis $i$ cannot also agree off a
different axis $j$ — agreeing off both axes forces all coordinates equal,
i.e. a single column. (For $n = 1$ the theorem is trivial.)
(4) The RS specialization (P5) plugs in $\#\text{circuits} =
\binom{s_i}{d}$ for MDS parity checks and the class structures of P5(b).

$\blacksquare$

-----------------------------------------------------------------------

## P5 (RS/GRS specialization)

Setting: $C_i = \Lambda_i \mathrm{RS}(S_i, t_i)$ over $F = F_q$,
$S_i \subseteq F_q^\star$ with $|S_i| = s_i$, dimension $t_i \ge 1$, and
$\Lambda_i$ a nonzero diagonal. For the column-circuit census we use the
GRS PARITY CHECK as the theorem's factor: $H_i$ is $(s_i - t_i) \times
s_i$ with rows evaluating $x \mapsto x^{t_i}, \dots, x^{s_i-1}$ at the
$s_i$ points of $S_i$ (a Vandermonde block). All columns are nonzero:
column $x$ is $(x^{t_i}, \dots, x^{s_i-1})^{\top} \ne 0$ since $x \ne 0$
(every entry $x^e \ne 0$). The lifted word space is $V = \sum_i
L_i(C_i)$; carrier membership on a support $T$ reads $\dim(V \cap F^T) =
|T| - \mathrm{rank}(H[:, T])$ with $H = H_1 \otimes \cdots \otimes H_n$
at $\Lambda = \mathrm{Id}$ (the column-census route; Run 1's C-1b caveat
about kron ROWSPACE annihilators does not touch this column-space
census). The theorem is applied to the matrix carrying the stated sparks:
here $H_i$, with spark $= s_i - t_i + 1 = d(C_i)$, the MDS distance of
$C_i$.

(a) **Spark of $H_i$:** any $k \le s_i - t_i$ columns of the Vandermonde
block (row set $\{t_i, \dots, s_i - 1\}$, $s_i - t_i$ rows) are
independent: the restriction to $k$ distinct nonzero evaluation points has
rank $k$ (a generalized Vandermonde/minor argument: the $k \times k$ minor
$\det(x_\alpha^{e_\beta})$ with distinct exponents $e_\beta$ and distinct
$x_\alpha$ is nonzero — the exponents $t_i, \dots, t_i + k - 1$ form a
consecutive set, and the generalized Vandermonde determinant factors as
$\prod x_\alpha^{t_i} \cdot$ (ordinary Vandermonde) $\neq 0$ over any
field). Any $s_i - t_i + 1$ columns are dependent (more columns than
rows). So $d_i = s_i - t_i + 1$: the MDS distance. Circuit sets of $H_i$:
every $(s_i - t_i + 1)$-subset is a circuit (dependent, and every proper
subset has size $\le s_i - t_i$, independent); there are $\binom{s_i}{d_i}$
of them.

(b) **Parallel classes of $H_i$.** Columns $x \ne y$ are parallel iff
their $(s_i - t_i)$-tuples are proportional. When $d_i \ge 3$ (so
$s_i - t_i \ge 2$), the block contains BOTH rows $x^{t_i}$ and
$x^{t_i+1}$: proportionality gives $\lambda x^{t} = y^{t}$ and $\lambda
x^{t+1} = y^{t+1}$ with $t = t_i$, $x, y \ne 0$, hence $\lambda = (y/x)^t
= (y/x)^{t+1}$, so $y/x = 1$, contradicting $x \ne y$. All classes are
singletons: $A_i = s_i$. When $d_i = 2$ ($s_i - t_i = 1$, the single row
$x^{t_i}$), all columns are nonzero scalars: ONE class of size $s_i$, so
$A_i = s_i^2$. When $d_i = 1$ ($s_i = t_i$), the parity check is EMPTY:
the factor acts as the zero map — in carrier language every product
column is a carrier of the weight-1 window; the $d \ge 2$ theorem does
not apply, and the registered trivial rows are counted as $N = \prod
s_i$.

(c) **LINE corollary — stated ONLY for $t_i = 1$.** $\mathrm{RS}(S, 1)$
is the constant code (word space $F \cdot
\mathbf{1}$), so its minimum-weight words are exactly the nonzero scalar
multiples of $\mathbf{1}$: supports are FULL axis lines. Hence, at
$t_i = 1$, the weight-$d$ carriers of $V = \sum_i L_i(C_i)$ are exactly
the supports of minimum words of the $L_i$-sum structure realized as axis
fibers over parity-check circuits; a carrier fiber on axis $i$ (with
$d_i = d$) corresponds to an $i$-embedded full line ($i$-minimum line).
Consequently H-MINLINE's LINE characterization is TRUE exactly when at
most one factor has $d_i = 2$ (no tied degenerate axes), and FAILS
otherwise by part (C) of the theorem:
*Run-1 reconciliation.* Instance $(13,(2,2,4), t=(1,1,1))$: $d_i = s_i =
(2,2,4)$, two degenerate axes ($i = 0, 1$), $d = 2$, $N = 16$. Parallel
classes: axis 0 (single row $x$ after $t_0 = 1$): one class of size 2,
$A_0 = 4$; axis 1: $A_1 = 4$; axis 2: $s_2 - t_2 + 1 = 4 \ge 3$, all
classes singletons, $A_2 = s_2 = 4$. Pair count: $(4 \cdot 4 \cdot 4 -
16)/2 = 24$. Fiber split: axis 0 contributes $s_1 s_2 \binom{2}{2} = 8$,
axis 1 contributes $8$; axis 2 contributes $0$ (singleton classes):
fibers $= 16$, diagonals $= 24 - 16 = 8$.
**24 = 16 + 8: Run 1's census reconciles exactly**, including its
exhibited flat-(0,12) diagonal pair.
*The line form FAILS for $t_i > 1$ even where the circuit theorem holds:*
minimum words of $\mathrm{RS}(S, t)$, $t \ge 2$: by the standard MDS
weight-structure argument, a nonzero polynomial $f$ of degree $< t$ has
weight $|S| - |\{a \in S : f(a) = 0\}|$, and $f$ has at most $\deg f \le
t - 1$ roots unless $f = 0$; so every nonzero word has weight $\ge s - (t
- 1)$, and the bound is attained exactly by the scalar multiples of
$f_Z(x) = \prod_{a \in Z}(x - a)$ for $(t-1)$-subsets $Z \subseteq S$
(degree $t - 1 < t$, vanishing exactly on $Z \cap S = Z$), whose support
is $S \setminus Z$ of size $s - t + 1 = d$.
At $t = 1$: $Z = \emptyset$, support $= S$ — full lines only. At $t \ge
2$: $|Z| = t - 1 \ge 1$, every
minimum support is a PUNCTURED line $S \setminus Z \subsetneq S$: no
minimum word is a full-line indicator, so the LINE phrasing is false
while the circuit form remains true. Registered example
$(13,(4,4,4),t=(2,1,1))$: sparks $(3,4,4)$, $d = 3$; carriers $= (\prod_{
j \ne 0} s_j)\binom{4}{3} = 16 \cdot 4 = 64$ in circuit form; every full
axis line has weight $4 > d = 3$, so NOT ONE carrier is a full line —
the circuit form gives 64 where the line form gives none. ✓

(d) **$\Lambda$-invariance (stated, Id-tested).** $\Lambda_i$ a nonzero
diagonal: columns map to nonzero scalar multiples of themselves;
proportional pairs map to proportional pairs and vice versa (diagonal
rescaling by nonzero $\lambda(x)$ is an invertible column operation);
circuits map to circuits. All spark/circuit/class structure is invariant;
the machine controls run at $\Lambda = \mathrm{Id}$.

$\blacksquare$

-----------------------------------------------------------------------

## Discharge statement

- **P1 discharged** (spark(H) = d; induction + dual isolation + hyperplane
  transport; regrouping legitimacy argued: columns of $B$ nonzero, spark
  claimed = min by IH).
- **P2 discharged** (Case 1 fiber along B; Case 2a impossible via
  separation of functionals; Case 2b circuit kernel 1-dim everywhere-
  nonzero, uniform no-case-split division, all $b_{v_j}$ parallel, P3'
  forces coincidence when $d_B \ge 3$).
- **P3 discharged** (forward trivial; converse by coordinate extraction
  functionals over an arbitrary field via basis extension).
- **P4 discharged** (regime lemma: (i) both $\ge 3$ fiber-only; (ii)
  $d=2$ one-spark-2: P3 forces fibers, P3 converse + spark-2 partner gives
  the pair dependence; (iii) both = 2 diagonals exist; counts: pair formula
  $(\prod A_i - N)/2$ with measured multiplicities, fiber subcount
  identity, fiber-regime carrier count with uniqueness of witness axis).
- **P5 discharged** (generalized-Vandermonde spark; MDS circuit counts;
  class structure per $d_i$; $t_i = 1$ LINE corollary ONLY at $t_i = 1$;
  Run-1 reconciliation 24 = 16 + 8; $t_i > 1$ line failure example
  $(4,4,4)$,t=(2,1,1) at 64 circuit carriers; zero-map $d_i = 1$ clause;
  $\Lambda$-invariance stated, tested at Id).

The theorem T-DGE (A), (B), (C) as stated in the prereg follows from P1
(spark), P2 + P3 + P4-regime (classification), P4 (counts), P5 (RS
instance). **No gap remains in the written argument.** Per the prereg's
verdict rule, promotion to FROZEN-CERTIFIED additionally requires ALL
control families M1-M6 to pass exactly in-run; see `controls_results.json`
in this run dir.
