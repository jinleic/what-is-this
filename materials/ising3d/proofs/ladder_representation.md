# Ladder AB/BB family: exact representation obstruction at eight rungs

## 1. Statement and outcome

For the open two-leg ladder \(\Lambda_L\), let

\[
 A_L=\sum_{v\in\Lambda_L}X_v,
 \qquad
 B_L=\sum_{(u,v)\in E(\Lambda_L)}Z_uZ_v,
\]

and use the depth filtration

\[
 F_1=\operatorname{span}\{A_L,B_L\},\qquad
 F_{k+1}=F_k+[A_L,F_k]+[B_L,F_k],\qquad
 D_k=\dim F_k.
\]

The intended theorem is \(D_{2L}\ge 2^L\) for every \(L\ge2\).

**[UNRESOLVED].** The all-length inequality is not proved here. The representation-evaluation route
instead gives an exact negative certificate for the proposed AB/BB family. That family is independent
through \(L=7\), but its 256 elements have rank exactly 248 over \(\mathbb Q\) at \(L=8\). Eight
explicit primitive integer relations are stored in the artifact and vanish coefficientwise in the
full Pauli basis. Thus no nonzero-determinant recurrence, tensor induction, Walsh minor, or
computational-basis restriction can prove all-length independence of this precise family: its full
operator Gram determinant is already zero at \(L=8\).

This is a family-specific obstruction, not a bound on \(D_{16}\). It neither proves nor refutes the
all-length target.

## 2. Exact AB/BB family and depth

Write a word such as `ABBA` for the left-nested adjoint
\(\operatorname{ad}_{A}\operatorname{ad}_{B}\operatorname{ad}_{B}(A)\), dividing every bracket by
two. Begin at two rungs with

\[
 W_2(00)=A,\quad W_2(01)=BA,\quad W_2(10)=ABA,\quad W_2(11)=BABA.
\]

For each later bit set

\[
 W_{L+1}(s0)=ABW_L(s),\qquad W_{L+1}(s1)=BBW_L(s). \tag{1}
\]

Each step adds two brackets, so every \(W_L(s)\) has depth at most \(2L\) and belongs to
\(F_{2L}\). If all \(2^L\) elements were independent for every \(L\), they would prove the target.

The exact Pauli convention is

\[
 Q_{(a\mid b)}=X^aZ^b,
 \qquad
 \frac12[Q_g,Q_v]=
 \begin{cases}
 0,&\langle g,v\rangle=0,\\
 \epsilon(g,v)Q_{g+v},&\langle g,v\rangle=1,
 \end{cases}
\]

with \(\epsilon(g,v)\in\{+1,-1\}\). Hence every family vector has integer coefficients.

## 3. Representation evaluation is injective

Let \(n=2L\), let \(\{|x\rangle:x\in\mathbb F_2^n\}\) be the computational basis, and define
matrix-unit functionals

\[
 \phi_{x,y}(T)=\langle y|T|x\rangle.
\]

For a Pauli basis vector one has exactly

\[
 \phi_{x,y}(Q_{(a\mid b)})
   =\mathbf 1_{y=x+a}(-1)^{b\cdot x}. \tag{2}
\]

**[LEMMA — Walsh matrix-unit transform].** The map
\(T\mapsto(\phi_{x,y}(T))_{x,y}\) is injective. If \(C\) is a row-by-Pauli coefficient matrix and
\(E\) is the corresponding row-by-matrix-unit evaluation matrix, then

\[
 EE^{\mathsf T}=2^n CC^{\mathsf T}. \tag{3}
\]

**Proof.** For fixed flip \(a=y-x\), equation (2) is the Walsh transform in the phase mask \(b\).
The Walsh characters obey

\[
 \sum_x(-1)^{b\cdot x}(-1)^{b'\cdot x}=2^n\mathbf 1_{b=b'}.
\]

Different flips have disjoint matrix-unit supports. This proves (3), and Walsh inversion proves
injectivity. \(\square\)

Thus representation evaluation genuinely escapes leading-Pauli-word arguments, but it cannot hide a
linear relation: the evaluation Gram determinant is nonzero exactly when the coefficient Gram
determinant is nonzero.

**[COMPUTATION — exact base determinants].** Direct integer construction of both sides of (3) gives
nonzero coefficient Gram determinants at \(L=2,3,4\):
\[
\det G_2=2^{15},
\]
\[
\det G_3=2^{30}3^6 5^4\cdot13\cdot3881,
\]
\[
\det G_4=
2^{76}3^{13}\cdot7\cdot11\cdot157^2\cdot2999\cdot175481
\cdot750229\cdot1368803\cdot30026837.
\]
The artifact stores each determinant, its integer factorization, the evaluation determinant, and the
exact identity

\[
 \det(EE^{\mathsf T})=(2^{2L})^{2^L}\det(CC^{\mathsf T}). \tag{4}
\]

These are finite base cases, not an inferred recurrence.

## 4. Exact finite ranks and first failure

Sparse elimination at the two recorded 31-bit primes gives

| \(L\) | family size | rank mod \(2147483647\) | rank mod \(2147483629\) |
|---:|---:|---:|---:|
| 2 | 4 | 4 | 4 |
| 3 | 8 | 8 | 8 |
| 4 | 16 | 16 | 16 |
| 5 | 32 | 32 | 32 |
| 6 | 64 | 64 | 64 |
| 7 | 128 | 128 | 128 |
| 8 | 256 | 248 | 248 |

**[LEMMA — modular lower bound].** Rank 248 modulo either prime proves rational rank at least 248,
because it exhibits an integer \(248\times248\) minor nonzero modulo that prime and hence nonzero
over \(\mathbb Q\). Agreement of the two primes is only a cross-check.

The upper bound is characteristic zero, not modular. Eight independent dependency rows were first
recovered at eight good primes. Their coefficients were combined by the Chinese remainder theorem,
rationally reconstructed, cleared to primitive integers, and then substituted back into the full
integer Pauli expansions. Every residual dictionary is empty. The eight relations have distinct
latest nonzero indices in insertion order, so their coefficient vectors are linearly independent.
Consequently

\[
 \operatorname{rank}_{\mathbb Q}\{W_8(s):s\in\{0,1\}^8\}
 \le 256-8=248.
\]

Together with the modular lower bound:

\[
 \boxed{\operatorname{rank}_{\mathbb Q}\{W_8(s)\}=248}. \tag{5}
\]

**[COMPUTATION — sparse exact certificates].** Two of the eight relations use only six and seven
family members. For example, the six-term relation is

\[
\begin{aligned}
0={}&-1105002572961\,W_8(00011110)
 +1784495478700\,W_8(00101110)\\
 &-823704922418\,W_8(00110110)
 +151372442200\,W_8(00111010)\\
 &-8978829696\,W_8(00111100)
 +336529200\,W_8(10111110).
\end{aligned} \tag{6}
\]

The seven-term relation is

\[
\begin{aligned}
0={}&-3215763564\,W_8(01011111)
 +9056323925\,W_8(01101111)\\
 &-9030348822\,W_8(01110111)
 +3800131005\,W_8(01111011)\\
 &-668129264\,W_8(01111101)
 +38382720\,W_8(01111110)\\
 &+4851000\,W_8(11111111).
\end{aligned} \tag{7}
\]

Equations (6) and (7) are exact identities in the complete Pauli coordinate space. All eight
relations, including the larger ones, appear in
`results/algebra_growth/ladder_representation.json`.

## 5. Boundary computational block

A more economical representation restriction keeps only Pauli terms that are identity on the first
\(L-3\) rungs. Equivalently, it evaluates the operator on the last-three-rung matrix block after
projecting all earlier rungs to the identity coefficient. This linear projection has exact ranks

| \(L\) | full family rank mod both primes | last-three-rung projection rank mod both primes |
|---:|---:|---:|
| 7 | 128 | 128 |
| 8 | 248 | 248 |

**[COMPUTATION].** The boundary block already sees the complete rank at these two consecutive
lengths. In particular, the eight-rung defect is already visible after discarding remote-support
coordinates: the small boundary projection has exactly the same rank deficiency. This observation
suggests a finite-state transfer description, but no all-length rank recurrence is claimed.

## 6. Precise negative conclusion

**[THEOREM — AB/BB family obstruction].** The recursive family (1) cannot prove
\(D_{2L}\ge2^L\) for every \(L\), by representation evaluation or by any other linear functional,
because its 256 elements are linearly dependent at \(L=8\). In particular, there is no recurrence
for full-family Gram determinants with nonzero value at every \(L\).

**Proof.** Equation (6), or any one of the stored exact residual-zero relations, is a nontrivial
integer linear dependence. Applying any collection of linear functionals preserves it, so every
256-column evaluation matrix has zero determinant. \(\square\)

The theorem rules out only the specified AB/BB family. It does not rule out replacing eight members,
using other depth-\(\le16\) words, constructing a different quotient representation, or proving the
target by another family.

## 7. Reproducibility

Run from the repository root:

```sh
.venv/bin/python experiments/e65_ladder_representation.py
.venv/bin/python tests/test_ladder_representation.py
```

Both commands print final `PASS` only if every embedded check succeeds. The standalone test
independently reconstructs the Pauli brackets, verifies the exact computational-evaluation Gram
identity at independent base cases, recomputes the consecutive \(L=7,8\) full and boundary ranks,
and substitutes all eight stored integer relations into fresh \(L=8\) vectors.
