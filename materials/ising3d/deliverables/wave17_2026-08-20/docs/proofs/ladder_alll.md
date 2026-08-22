# Two-leg ladder all-length search: exact multi-word obstruction

## 1. Statement and outcome

For the open two-leg ladder \(\Lambda_L\), put

\[
 A_L=\sum_{v\in\Lambda_L}X_v,
 \qquad
 B_L=\sum_{(u,v)\in E(\Lambda_L)}Z_uZ_v,
\]

and define the depth filtration

\[
 F_1=\operatorname{span}\{A_L,B_L\},\qquad
 F_{k+1}=F_k+[A_L,F_k]+[B_L,F_k],\qquad
 D_k=\dim F_k.
\]

The target is \(D_{2L}\geq 2^L\) for every \(L\geq2\).

**[UNRESOLVED].** This investigation does not prove or refute that all-length inequality. It does
verify the existing binary multi-word family one length farther, through \(L=7\), and it gives an
exact obstruction to the requested triangular induction. The obstruction is already present in the
base lift \(L=2\to3\): after exact rational echelonization, eight independent lifted elements have
only seven leading Pauli words. The two colliding elements are not proportional; an exact
\(2\times2\) coefficient minor has determinant \(6\).

More strongly, **every** order in the natural class of 96 product-lexicographic Pauli orders fails
at this same base lift. The maximum is seven distinct leaders, whereas eight are required. This is a
finite obstruction to that order class and to this particular pair of extension moves, not an
all-orders no-go theorem and not evidence against the inequality itself.

## 2. Exact conventions

Write

\[
 Q_{(a\mid b)}=X^aZ^b,
 \qquad (a\mid b)\in\mathbb F_2^{4L}.
\]

The integer encoding is `a | (b << 2L)`. The experiment divides every commutator by two:

\[
 \frac12[Q_g,Q_v]=
 \begin{cases}
  0,&\langle g,v\rangle=0,\\
  \epsilon(g,v)Q_{g\oplus v},&\langle g,v\rangle=1,
 \end{cases}
 \qquad \epsilon(g,v)\in\{+1,-1\}.
\]

All Pauli coefficients and all collision minors are computed with exact integers or
`fractions.Fraction`. Modular echelonization uses

\[
 p_1=2147483647,\qquad p_2=2147483629.
\]

**[LEMMA — modular lower bound].** A rank \(r\) modulo either prime is a rigorous lower bound on
the rank over \(\mathbb Q\).

**Proof.** Modular rank \(r\) exhibits an \(r\times r\) integer minor whose determinant is nonzero
modulo the prime. That integer determinant is nonzero, so the same rows have rational rank at least
\(r\). \(\square\)

Agreement at two primes is only an arithmetic cross-check; it is not claimed to prove rational rank
equality.

## 3. Recomputed depth bases and their leading words

Fix the product order that compares rungs from right to left, top before bottom on each rung, with
local descending order

\[
 X>Y>Z>I.
\]

`experiments/e53_ladder_alll.py` expands every nonzero generator-left-nested word through depth
\(2L\) for \(2\leq L\leq5\), and performs sparse modular echelonization in that fixed monomial
order. Each accepted row is a genuine bracket word in \(A_L,B_L\); its stored pivot is therefore the
leading Pauli word of a concrete multi-word basis element after triangular reduction.

**[COMPUTATION].** The recomputed ranks and distinct echelon leaders are

| \(L\) | depth | rank modulo \(p_1\) | distinct echelon leaders |
|---:|---:|---:|---:|
| 2 | 4 | 6 | 6 |
| 3 | 6 | 16 | 16 |
| 4 | 8 | 38 | 38 |
| 5 | 10 | 85 | 85 |

The JSON artifact stores every basis word, its depth, leading Pauli code and decoded rung string,
raw support size, and SHA-256 digest of all exact `(code, coefficient)` pairs. These data demonstrate
the intended escape from the wave-6 single-word obstruction: an element need not be a single Pauli
word to participate in a triangular certificate.

As an independent implementation control, the repository
`ising.clifford.fast_lie.GeneratedAlgebra` engine recomputes the complete algebra dimensions
\(11\) and \(263\), with support-closure sizes \(56\) and \(1056\), for \(L=2,3\), respectively.
These complete-algebra dimensions are controls, not claims for arbitrary \(L\).

## 4. The tested multi-word induction

The binary family begins at \(L=2\) with

\[
\begin{array}{c|c}
00&A,\\
01&[B,A],\\
10&[A,[B,A]],\\
11&[B,[A,[B,A]]].
\end{array}
\]

For every later bit, define the two local-in-word extension moves

\[
 T_0(W)=[A,[B,W]],
 \qquad
 T_1(W)=[B,[B,W]]. \tag{1}
\]

Thus \(w_{L+1}(s0)=T_0(w_L(s))\) and
\(w_{L+1}(s1)=T_1(w_L(s))\). Every step adds exactly two brackets, so every member has depth at
most \(2L\) and lies in \(F_{2L}\).

**[COMPUTATION — finite rank certificate].** Exact expansion followed by independent reduction at
both primes gives

| \(L\) | family size | rank mod \(p_1\) | rank mod \(p_2\) |
|---:|---:|---:|---:|
| 2 | 4 | 4 | 4 |
| 3 | 8 | 8 | 8 |
| 4 | 16 | 16 | 16 |
| 5 | 32 | 32 | 32 |
| 6 | 64 | 64 | 64 |
| 7 | 128 | 128 | 128 |

Consequently \(D_{2L}\geq2^L\) is rigorously certified for \(2\leq L\leq7\). This is a finite
computation, not an all-length theorem.

## 5. Exact failure of triangular leading-term propagation

A legitimate leading-term induction must work on basis elements, not merely on raw bracket words.
At length \(L\), the experiment therefore echelonizes the \(2^L\) family exactly over
\(\mathbb Q\), retaining each basis element as a rational linear combination of the original
bracket words. It then applies each move in (1) linearly to that same combination evaluated on the
\((L+1)\)-rung ladder. Every lifted element remains an explicit rational linear combination of
bracket words and hence lies in \(F_{2L+2}\).

For the fixed order, the numbers of lifted rows and distinct leaders are

| step | lifted elements | distinct leaders |
|---:|---:|---:|
| \(2\to3\) | 8 | 7 |
| \(3\to4\) | 16 | 13 |
| \(4\to5\) | 32 | 16 |
| \(5\to6\) | 64 | 27 |
| \(6\to7\) | 128 | 42 |

Thus the candidate loses leading-term control at the first possible doubling step.

**[COMPUTATION — minimal exact counterexample].** Number the four rational-echelon parent rows at
\(L=2\) in insertion order. The \(T_0\) child of parent 1 and the \(T_1\) child of parent 3 share
leader

\[
 \mathtt{II\ ZI\ XY}
\]

(code 2212), with coefficients \(-2\) and \(4\), respectively. On the two Pauli supports

\[
 q=\mathtt{YI\ XI\ ZI}\quad(323),
 \qquad
 r=\mathtt{XY\ ZI\ II}\quad(649),
\]

their coefficient minor has determinant

\[
 \det\begin{pmatrix}
 c_1(q)&c_1(r)\\
 c_2(q)&c_2(r)
 \end{pmatrix}=6\ne0. \tag{2}
\]

Equation (2) proves that the collision is specifically a failure of distinct leading terms; it is
not proportionality or linear dependence of the two lifted elements.

## 6. Exhaustive finite locality/order check

The requested locality argument would need a position-independent local leading rule. Before trying
to promote such a rule, the experiment exhausts the standard bounded local order class:

\[
 2\text{ rung directions}\times
 2\text{ within-rung directions}\times
 24\text{ strict local orders on }\{I,X,Y,Z\}=96.
\]

For each order it repeats the exact rational base echelonization and the eight lifts from
\(L=2\to3\). The histogram of distinct lifted leaders is

| distinct leaders | number of orders |
|---:|---:|
| 5 | 16 |
| 6 | 36 |
| 7 | 44 |

**[COMPUTATION — product-order obstruction].** No order in this class reaches the required eight
leaders. Because the failure occurs on the three-rung base window itself, translation covariance or
checking the same local rule at additional positions cannot repair it. Any successful induction
must change at least one essential ingredient: the extension maps, the parent family, or the
monomial order beyond this product-lex class.

This is the precise point where a proposed all-length proof stops. A finite local calculation can
support an induction only after its statement is true on the complete local window. Here the
candidate local statement is exactly false, so no locality argument is written and no finite
verification is promoted to an all-length theorem.

## 7. Reproducibility and scope

Run from the repository root:

```sh
.venv/bin/python experiments/e53_ladder_alll.py
.venv/bin/python tests/test_ladder_alll.py
```

The experiment and standalone test each print final `PASS`. The standalone test independently
reimplements the Pauli brackets, recomputes the complete \(L=6,7\) families and both modular ranks,
checks every stored leader and exact-vector digest at those two consecutive lengths, reconstructs
the determinant-6 lift collision, exhausts the 96 orders, and validates the result envelope.

The machine-readable certificate is
`results/algebra_growth/ladder_alll.json`. Its status and target are explicitly tagged
`UNRESOLVED`; no asymptotic extrapolation or benchmark fitting is used.
