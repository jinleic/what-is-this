# Complete circuit classification for two 2-row GRS factors
## gate `H-2ROW-COMPLETE-3ROW-OPEN`, run `20260902T023120Z_8ac826db_96f766bfe6cc`

Agent `RsPe3dH2`, 2026-09-02 UTC. Preregistration:
`prereg/H_2ROW_COMPLETE_3ROW_OPEN_PREREG_2026-09-01.md`, path-scoped
commit `4486a168f8cc11a07f3e1585531e768502ca3408`, SHA-256
`18f60f53b292e873f31d01231f8634d24250c37e237eabb76c2859c402ee60b1`,
with a byte-identical in-run copy at `pre_statement.md`. There was no amendment
and no parameter-dependent computation before the preregistration commit.

## 1. Statement

Let

\[
a_x=(1,x)^{\mathsf T}\quad(x\in X),\qquad
b_y=(1,y)^{\mathsf T}\quad(y\in Y)
\]

be two 2-row GRS/Vandermonde factors over a field, with distinct evaluation
points within each factor and arbitrary nonzero column multipliers. The
multipliers do not change the represented matroid, so suppress them. The
product columns are

\[
h(x,y)=a_x\otimes b_y=(1,y,x,xy)^{\mathsf T}\in F^4.
\]

Then every circuit of the product has size 3, 4, or 5, and the following list
is exhaustive.

1. **Size 3:** exactly an axis fiber supported by a factor triple.
2. **Size 4:** exactly either
   - a profile-$(3,3)$ crossing support, or
   - a profile-$(4,4)$ all-distinct matching whose four paired points preserve
     cross-ratio, equivalently lie on one nondegenerate bilinear/Möbius graph.
3. **Size 5:** exactly the graph/minor families in the table below.
4. **All other sizes:** none. In particular, no circuit has size at least 6.

For a five-edge bipartite support, write $P_k$ for a path on $k$ vertices and
$K_2$ for an isolated edge.

| profile | graph | exact remaining predicate |
|---|---|---|
| $(3,3)$ | $C_4+K_2$ | automatic |
| $(3,4)$ | $P_5+K_2$ | automatic |
| $(3,5)$ | two row-centered $P_3$'s $+K_2$ | automatic |
| $(4,3),(5,3)$ | transposes of the preceding two | automatic |
| $(4,4)$ | $P_4+2K_2$ | its unique all-distinct 4-minor is nonzero |
| $(4,5)$ | one row-centered $P_3+3K_2$ | both all-distinct 4-minors are nonzero |
| $(5,4)$ | transpose | both all-distinct 4-minors are nonzero |
| $(5,5)$ | $5K_2$ | all five 4-minors are nonzero |

Equivalently, a five-set is a circuit iff its maximum row and column degree is
at most two, it contains no size-four crossing, and every all-distinct
four-subset has nonzero determinant.

## 2. Ambient-rank bound

The following observation does not require the GRS hypothesis. If a matrix
$H$ has column rank at most $R$ and $S$ is a circuit of $H$, then

\[
\operatorname{rank}H_S=|S|-1\le R,
\]

and consequently

\[
|S|\le R+1. \tag{1}
\]

For two factor matrices of row ranks at most $r_A,r_B$, the tensor product has
rank at most $r_A r_B$. Thus every two-factor product circuit obeys

\[
\boxed{|S|\le r_A r_B+1}. \tag{2}
\]

This bound is attainable, not merely formal: the size-five families in §1
have rank four and attain (2) for $r_A=r_B=2$. This is an existence claim;
it does not say every pair of factors attains equality.

In the present 2-row-by-2-row setting, (2) gives $|S|\le5$. More concretely,
every five columns in $F^4$ are dependent, so every six-set contains a
proper dependent five-subset and cannot be a circuit. The same proper subset
excludes every larger set.

## 3. Why the list is exhaustive

Both factors are MDS with spark three. The certified minimum-circuit theorem
therefore gives exactly the factor triples on one axis at size three.

At size four, the certified near-minimum and crossing theorems leave two and
only two possibilities. A repeated-coordinate circuit has profile $(3,3)$
and, by the complete Theorem-X converse, is a crossing. An all-distinct
support has profile $(4,4)$, and the certified bilinear determinant theorem
says it is a circuit exactly when

\[
\det[1,y,x,xy]_{(x,y)\in S}=0,
\]

which is equivalent to simultaneous cross-ratio equality and to containment
in the graph of a unique nondegenerate Möbius transformation. The intermediate
profiles $(3,4)$ and $(4,3)$ are independent: if, for example, two cells share
$x=x_0$, a bilinear polynomial vanishing at their two distinct $y$ values has
both coefficients at $1,y$ vanish at $x_0$; its two remaining zeros at distinct
$x$ and distinct $y$ force the polynomial to be zero. Hence no further
size-four profile occurs.

At size five, dependence is automatic in $F^4$, and circuitness is exactly
independence of all five four-deletions. A row or column of degree three
contains a proper factor-triple circuit, so maximum degree is two. A graph of
maximum degree two is a disjoint union of paths and even cycles. Exhausting
those five-edge graphs, and excluding precisely the already-classified
size-four circuits among their deletions, gives the table in §1. This is the
frozen `H-DP2-SIZE5` theorem; its proof is graph-theoretic and uses the exact
four-minors, not a finite-field fit. Combining it with (1) completes the list.

The field-free size-five counts for $|X|=|Y|=n$ are

\[
N_{33}=9\binom n3^2,
\]

\[
N_{34}=N_{43}=72\binom n3\binom n4,
\]

\[
N_{35}=N_{53}=18n\binom n3\binom{n-1}4.
\]

Writing $M_4=24\binom{|X|}4\binom{|Y|}4$ and
$N_{\rm all}^{(4)}(X,Y)$ for the all-distinct size-four circuit count,

\[
N_{44}^{(5)}=12\bigl(M_4-N_{\rm all}^{(4)}(X,Y)\bigr).
\]

The $(4,5),(5,4),(5,5)$ populations remain exact finite determinant sums and
are genuinely field-dependent; the theorem does not replace them by a false
field-free polynomial.

## 4. Correct degree-bound scope

For tied spark-$d$ MDS factors, any circuit $S$ with $|S|>d$ satisfies

\[
\deg_{\rm row},\deg_{\rm column}\le d-1. \tag{3}
\]

Indeed, $d$ cells on one axis contain a factor $d$-circuit, which would be a
proper dependent subset of $S$.

For unequal MDS sparks $d_A,d_B$, with A indices drawn as row vertices, the
correct asymmetric statement is

\[
\deg_A\le d_B-1,\qquad \deg_B\le d_A-1 \tag{4}
\]

for a nonfiber circuit larger than both minimum fiber sizes. The factor named
on the right supplies the forbidden fiber. A universal bound by
$\min(d_A,d_B)-1$ is false: when $d_B>d_A$, a minimum B-fiber has size $d_B$
and is itself a circuit, with A-row degree $d_B$. This is the required explicit
exception and fixes the earlier convention hazard.

For the tied 2-row case, (3) specializes to maximum degree two for every
circuit of size four or five.

## 5. In-run exact controls

`run_frontier.py` used Python integer arithmetic modulo $p$ and true row-major
Kronecker columns. It passed all 121 aggregate assertions. For each of
$n=4,5$ and $p=7,13$ it exhausted every six-subset:

| $n$ | six-sets swept per prime | six-circuits | every six-set has dependent five-subset |
|---:|---:|---:|---:|
| 4 | $\binom{16}{6}=8{,}008$ | 0 | yes |
| 5 | $\binom{25}{6}=177{,}100$ | 0 | yes |

A known size-five support was accepted with rank four and all four-deletions
independent, while a planted six-set was rejected because it contains a
dependent five-subset. These controls corroborate both directions of the
ambient edge; the proof of (1), not the census, supplies the general result.

The run consumed 40.594 CPU seconds and 42.775 wall seconds under the
preregistered caps of 600 and 900 seconds, with `nice` level 15, all declared
thread counts pinned to one, `PYTHONDONTWRITEBYTECODE=1`, and
`sys.dont_write_bytecode=True`.

## 6. Status

The theorem in §1 and the general bound (2) are **proved and certified**.
Their proof dependencies are the frozen size-three theorem, complete
Theorem-X crossing converse, all-distinct cross-ratio theorem, and complete
size-five graph/minor theorem. The companion higher-row investigation in this
same run is deliberately separate and remains incomplete; its inconclusive
status does not weaken this two-row capstone.
