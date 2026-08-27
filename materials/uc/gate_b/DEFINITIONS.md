# Gate B: authoritative definitions

This file freezes the definitions used by the Gate B theorem. They are the
current definitions in `math/uc`, not a replacement objective.

## Admissible families

Let \(\mathcal F\subseteq 2^{[n]}\) be a nonempty finite simple family and put
\(m=|\mathcal F|\). All random rows below are uniform on \(\mathcal F\).
The repository's Gate B admissibility conditions are

\[
\max_{i\in[n]}|\{A\in\mathcal F:i\in A\}|
\le \left\lfloor\frac{2m}{5}\right\rfloor
\tag{cap}
\]

and

\[
\sum_{A\in\mathcal F}|A|\ge \frac{m\log_2m}{2}.
\tag{Reimer incidence}
\]

The search scripts exclude the trivial sizes \(m<3\). This convention is
immaterial for the construction, whose sizes are \(45^k\).

The exact integer Reimer threshold is

\[
R_m=\left\lceil\frac{m\log_2m}{2}\right\rceil.
\]

It is computed without floating point: \(R_m\) is the least integer \(r\)
for which \(2^{2r}\ge m^m\). The previous `n=7` producer used
`ceil(m*log2(m)/2 - 1e-12)`. An exact audit found the same values for every
\(3\le m\le128\), so the prior census is unchanged; the producer and all new
Gate B code now use the integer comparison.

A family is **active** if no coordinate column is zero and **separating** if
its coordinate columns are pairwise distinct. The repository calls active and
separating families **normalized**. Normalization was a search reduction and
reporting condition, not an extra condition in the displayed definition of
\(c_{\rm cl}^\star\). This distinction does not affect the theorem: the base
family and every Cartesian power used below are normalized.

## Fixed constant

The repository fixes

\[
\alpha=0.0356069=\frac{356069}{10^7}.
\]

Gate B does not optimize or reinterpret \(\alpha\).

## Shapley iid cost

Fix a coordinate order \(\pi\). At a coordinate \(i\), after a predecessor
prefix \(u\) of a uniform row, put

\[
p_i(u)=\Pr[X_i=1\mid X_{\operatorname{pred}_\pi(i)}=u].
\]

For two independent uniform prefixes \(u,v\), the conditional OR probability
is \(p_i(u)+p_i(v)-p_i(u)p_i(v)\). With binary entropy \(h\) in bits, define

\[
Q_\pi(\mathcal F)
=\sum_i \mathbb E_{u,v}
 h\!\left(p_i(u)+p_i(v)-p_i(u)p_i(v)\right).
\]

Then \(Q(\mathcal F)=\mathbb E_\pi Q_\pi(\mathcal F)\), where \(\pi\) is
uniform over all \(n!\) coordinate orders. The predecessor-subset Shapley
formula in the code is exactly this order average.

## One-sided causal Bellman cost

At a pair of row prefixes \((u,v)\), write
\(p=p_i(u)\), \(r=p_i(v)\), and

\[
s^*(p,r)=\max\{p,r,\min(p+r,1/2)\},
\qquad U(p,r)=\min(1,p+r).
\]

A one-sided causal action chooses the conditional OR probability
\(s\in[s^*(p,r),U(p,r)]\). Its four next-bit probabilities are

\[
P_{00}=1-s,\quad P_{10}=s-r,\quad
P_{01}=s-p,\quad P_{11}=p+r-s.
\]

For a fixed order, \(C_{+,\pi}(\mathcal F)\) is the maximum, over all such
causal policies, of the expected sum of \(h(s)\). Equivalently, with
continuation values \(V_{ac}\), its Bellman recurrence maximizes

\[
h(s)+(1-s)V_{00}+(s-r)V_{10}+(s-p)V_{01}+(p+r-s)V_{11}
\]

over the same interval. If
\(D=-V_{00}+V_{10}+V_{01}-V_{11}\), the unconstrained maximizer is
\(1/(1+2^{-D})\), clamped to \([s^*,U]\). Finally,

\[
C_+(\mathcal F)=\mathbb E_\pi C_{+,\pi}(\mathcal F).
\]

## Gate B objective and closure defect

The objective is

\[
A_+(\mathcal F)
=(1-\alpha)Q(\mathcal F)+\alpha C_+(\mathcal F)-\log_2m.
\]

The closure defect uses **ordered pairs with replacement**:

\[
\varepsilon_\vee(\mathcal F)
=\Pr_{X,Y\sim\operatorname{Unif}(\mathcal F)}
 [X\cup Y\notin\mathcal F]
=\frac{|\{(X,Y)\in\mathcal F^2:X\cup Y\notin\mathcal F\}|}{m^2}.
\]

The exact Gate B quantity recorded before this round was

\[
c_{\rm cl}^\star=
\sup_{\substack{\mathcal F\text{ cap/Reimer}\\A_+(\mathcal F)<0}}
\frac{-A_+(\mathcal F)}{\varepsilon_\vee(\mathcal F)}.
\]
The quotient is interpreted in the extended nonnegative reals. If an
admissible family has \(A_+(\mathcal F)<0\) and
\(\varepsilon_\vee(\mathcal F)=0\), its value is \(+\infty\); otherwise the
displayed quotient is used with \(\varepsilon_\vee(\mathcal F)>0\). Equivalently,
one may first ask whether such a zero-defect negative family exists and, if
not, take the supremum only over positive defects. The construction below has
\(0<\varepsilon_\vee(\mathcal F_k)<1\), so its conclusion does not depend on
this convention.

There is no small-defect restriction in this supremum. The phrase
"near-union-closed" in the search plan described a proposed discovery method,
not a changed admissible class. A separate local-stability question with
\(\varepsilon_\vee\to0\) remains meaningful, but it is not Gate B as defined
above.

## Cartesian products

For \(\mathcal F\subseteq2^{[n]}\) and
\(\mathcal G\subseteq2^{[d]}\), place the ground sets in disjoint blocks and
define

\[
\mathcal F\boxtimes\mathcal G
=\{A\cup(n+B):A\in\mathcal F,\ B\in\mathcal G\}.
\]

The \(k\)-fold power is \(\mathcal F^{\boxtimes k}\). This is the product used
in the Gate B theorem; it is not a dummy-coordinate lift and does not alter any
normalization or probability convention.
