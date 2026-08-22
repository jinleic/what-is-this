# Compact local commutant of the nondegenerate cubic transverse-field Ising generator

Status labels in this note are literal: **[THEOREM]** and **[LEMMA]** are
proved below; **[COMPUTATION]** is an exact finite check only; **[UNRESOLVED]**
marks scope deliberately not decided here.

---

## 1. Algebra, support, and the meaning of the infinite sums

Let \(\Lambda=\mathbb Z^3\), with nearest-neighbour bonds \(E\), and let
\(e_1,e_2,e_3\) be the coordinate unit vectors. For a finite
\(F\Subset\Lambda\), write

\[
\mathcal A_F=\bigotimes_{x\in F}M_2(\mathbb C)_x
\]

embedded in the spin algebra by identities outside \(F\). The **compact local
algebra** is

\[
\mathcal A_{\rm loc}=\bigcup_{F\Subset\Lambda}\mathcal A_F.
\]

An operator in this note always means an element of \(\mathcal A_{\rm loc}\),
unless stated otherwise. Its support is the least finite set
\(\operatorname{supp}(Q)\) on which it acts nontrivially. Equivalently, in the
repository's ordered-Pauli basis

\[
Q_{a,b}=X^aZ^b,\qquad Q=\sum_{a,b}q_{a,b}Q_{a,b},
\]

where only finitely many masks are nonzero,

\[
\operatorname{supp}(Q)=
\bigcup_{q_{a,b}\ne0}\{x:(a_x,b_x)\ne(0,0)\}.
\tag{1}
\]

This is the same canonical convention used by `e52` and `e61`; at one site
\(XZ=-iY\). Thus referring below to a \(Y\)-coefficient introduces no new
Pauli convention: it is the \(XZ\)-slot written in the usual Hermitian basis.

For \(a,b\in\mathbb C\), define the *derivation notation*

\[
H_{a,b}=a\sum_{x\in\Lambda}X_x+
        b\sum_{\{x,y\}\in E}Z_xZ_y
\tag{2}
\]

only through its commutator on local operators:

\[
\delta_{a,b}(Q):=[Q,H_{a,b}]
 :=a\sum_{x\in\Lambda}[Q,X_x]
   +b\sum_{\{x,y\}\in E}[Q,Z_xZ_y].
\tag{3}
\]

There is no assertion that either bare infinite sum in (2) is an element of
\(\mathcal A_{\rm loc}\) or of its norm closure.

**Lemma 1 [LEMMA — finite derivations].** If
\(S=\operatorname{supp}(Q)\), then (3) is the finite expression

\[
\delta_{a,b}(Q)=
 a\sum_{x\in S}[Q,X_x]
 +b\sum_{\substack{\{x,y\}\in E\\\{x,y\}\cap S\ne\varnothing}}
      [Q,Z_xZ_y].
\tag{4}
\]

Moreover,
\(\operatorname{supp}(\delta_A(Q))\subseteq S\) and
\(\operatorname{supp}(\delta_B(Q))\subseteq S\cup N(S)\), where \(N(S)\)
is the nearest-neighbour set of \(S\).

*Proof.* A generator supported disjointly from \(S\) commutes with
\(\mathcal A_S\). There are at most \(|S|\) field terms and at most
\(6|S|\) bonds meeting \(S\). The support bounds follow from tensor-factor
locality of a commutator. \(\square\)

For unambiguous coefficient extraction, fix the two-factor order with the
named site first:

\[
\mathcal A_{\rm loc}
 \cong M_2(\mathbb C)_u\otimes
        \mathcal A_{\rm loc}(\Lambda\setminus\{u\}).
\]

Every local operator has the unique expansion

\[
O=I_u\otimes O_{00}+X_u\otimes O_{10}
  +Z_u\otimes O_{01}+X_uZ_u\otimes O_{11}.
\tag{5}
\]

Define \(\pi^{(u)}_{\alpha\beta}(O)=O_{\alpha\beta}\). Thus the first
factor of every displayed tensor product is always the named site \(u\), and
the coefficient is an operator on its complement. Uniqueness follows by
grouping the ordered-Pauli basis in (1) by the two bits at \(u\).

---

## 2. The all-size geometric step

**Lemma 2 [LEMMA — private outside neighbour].** Let
\(\varnothing\ne S\Subset\mathbb Z^3\). Choose \(v\in S\) maximizing its
first coordinate, and put \(w=v+e_1\). Then

\[
w\notin S,\qquad S\cap N(w)=\{v\}.
\tag{6}
\]

*Proof.* The point \(w\) has first coordinate \(v_1+1\), so it is not in
\(S\). Its six neighbours are

\[
N(w)=\{v,\ w+e_1,\ w+e_2,\ w-e_2,\ w+e_3,\ w-e_3\}.
\]

Every displayed neighbour other than \(v\) has first coordinate at least
\(v_1+1\), strictly above the maximum attained in \(S\). Hence none lies in
\(S\), while \(v\in S\). \(\square\)

The point of (6) is not merely that \(w\) is outside the support: it gives a
single, uniquely supported boundary-bond contribution to a tensor slot.

---

## 3. The two coefficient isolations

**Lemma 3 [LEMMA — the \(Z_w\)-slot isolates \([Q,Z_v]\)].** Let
\(Q\in\mathcal A_{\rm loc}\), let \(S=\operatorname{supp}(Q)\ne\varnothing\),
and choose \(v,w\) as in Lemma 2. Then

\[
\pi^{(w)}_{01}\bigl(\delta_{a,b}(Q)\bigr)
 =b[Q,Z_v].
\tag{7}
\]

In particular, if \(b\ne0\) and \(\delta_{a,b}(Q)=0\), then
\([Q,Z_v]=0\).

*Proof.* Since \(w\notin S\), \([Q,Z_w]=0\). Every field term in (4)
is supported on \(S\), and every bond term not incident to \(w\) has the
identity in the \(w\)-factor. None contributes to
\(\pi^{(w)}_{01}\).

For a bond incident to \(w\), with other endpoint \(u\in N(w)\), the
Leibniz rule gives

\[
[Q,Z_wZ_u]=[Q,Z_w]Z_u+Z_w[Q,Z_u]=Z_w[Q,Z_u].
\tag{8}
\]

Here the tensor order is exactly the one fixed in (5): the factor \(Z_w\) is
first and \([Q,Z_u]\) acts on the complement of \(w\). Thus its
\(Z_w\)-coefficient is \([Q,Z_u]\).

It is important to account explicitly for the other five bonds at \(w\).
For \(u\notin S\), \([Q,Z_u]=0\), since \(Q\) acts trivially at \(u\).
Those formal bond terms may be reinserted into the finite sum in (4) as zero;
they do not hide a boundary contribution. Lemma 2 says that the only
\(u\in N(w)\cap S\) is \(v\). Therefore the sum of all w-incident
contributions is exactly \(b[Q,Z_v]\), proving (7). \(\square\)

**Lemma 4 [LEMMA — \(Z_v\)-centrality removes the \(X_v,Y_v\) slots].**
If \([Q,Z_v]=0\), then uniquely

\[
Q=I_v\otimes Q_I+Z_v\otimes Q_Z.
\tag{9}
\]

*Proof.* In the usual Hermitian spelling of the same local basis, write

\[
Q=I_v\otimes Q_I+X_v\otimes Q_X+Y_v\otimes Q_Y+Z_v\otimes Q_Z.
\]

Using \([X,Z]=-2iY\) and \([Y,Z]=2iX\),

\[
[Q,Z_v]=-2iY_v\otimes Q_X+2iX_v\otimes Q_Y.
\]

Coefficient uniqueness forces \(Q_X=Q_Y=0\), giving (9). \(\square\)

In the ordered convention used by the exact calculation, this is the
injective map on the \(a_v=1\) sector

\[
\frac12[Z_v,Q_{a,b}]
 =-\mathbf1_{a_v=1}Q_{a,b+e_v}.
\tag{10}
\]

The map \(b\mapsto b+e_v\) is a permutation, so there can be no cancellation
among those coefficients. Equation (10) also shows why the proof is clean
over \(\mathbb C\) while the finite matrices can be recorded with integer
entries. The same ordered-basis proof works over any field of characteristic
not two; characteristic two is intentionally excluded.

**Lemma 5 [LEMMA — the \(Y_v\)-slot peels the remaining \(Z_v\) slot].**
Assume (9). Then

\[
\pi^{(v)}_Y\bigl(\delta_{a,b}(Q)\bigr)=2ia\,Q_Z,
\tag{11}
\]

where \(\pi^{(v)}_Y\) denotes the coefficient of \(Y_v\) in the Hermitian
form of (5). Consequently, if \(a\ne0\) and \(\delta_{a,b}(Q)=0\), then
\(Q_Z=0\).

*Proof.* The \(x=v\) field term is

\[
[Q,aX_v]=a[Z_v,X_v]\otimes Q_Z
         =2iaY_v\otimes Q_Z.
\tag{12}
\]

For \(x\ne v\), \([Q,X_x]\) leaves the \(v\)-factor in the span of
\(I_v,Z_v\), so it has no \(Y_v\)-component.

The bond terms have the same property. If \(v\notin\{x,y\}\), then
\(Z_xZ_y\) acts entirely on the second factor in (9), so its commutator
preserves the local span \(\operatorname{span}\{I_v,Z_v\}\). If the bond is
\(\{v,y\}\), use Lemma 4 in the same Leibniz identity as before:

\[
[Q,Z_vZ_y]=[Q,Z_v]Z_y+Z_v[Q,Z_y]=Z_v[Q,Z_y].
\tag{13}
\]

Again this has only \(I_v\) and \(Z_v\) factors. Hence (12) is the unique
\(Y_v\)-contribution, proving (11). Since \(2ia\ne0\), it forces
\(Q_Z=0\). \(\square\)

For comparison with the exact JSON maps, let
\(\operatorname{ad}_h=\frac12[h,\cdot]\), the generator-left orientation
used in `e52` and `e61`. On the subspace selected by Lemma 4 (the
\(a_v=0\) columns), the same statement is the integer coefficient identity

\[
\pi^{(v)}_{11}\bigl(a\operatorname{ad}_A+b\operatorname{ad}_B\bigr)
 =a\,\pi^{(v)}_{01}.
\tag{14}
\]

Here both sides are coefficients on the complement of \(v\). Without the
\(a_v=0\) restriction, bond and other-field terms can have a general
\((1,1)_v\) component; the stored maps retain those entries rather than
silently discarding them. This is exactly why Lemma 4 is a logically prior
step.

---

## 4. Main theorem and simultaneous-centralizer corollary

**Theorem 6 [THEOREM — compact local commutant for one nondegenerate
\(H_g\)].** Let \(a,b\in\mathbb C\) satisfy \(ab\ne0\). Then

\[
\boxed{
\{Q\in\mathcal A_{\rm loc}:\delta_{a,b}(Q)=0\}=\mathbb C I.
}
\tag{15}
\]

Equivalently, every compact finite-support operator commuting with the formal
transverse-field Ising generator

\[
H_g=a\sum_{x\in\mathbb Z^3}X_x+
    b\sum_{\langle x,y\rangle}Z_xZ_y
\]

in the finite-commutator-derivation sense is scalar.

*Proof by induction on \(m=|\operatorname{supp}(Q)|\).* For \(m=0\), a
local operator supported on the empty set is a scalar multiple of \(I\).

Let \(m>0\), assume the statement for all supports of size at most \(m-1\),
and take \(Q\) with support size \(m\) and \(\delta_{a,b}(Q)=0\). Choose
\(v\) with maximal first coordinate. Lemma 3 gives \([Q,Z_v]=0\), because
\(b\ne0\). Lemma 4 gives the form (9). Lemma 5 gives \(Q_Z=0\), because
\(a\ne0\). Hence

\[
Q=I_v\otimes Q_I,
\qquad
\operatorname{supp}(Q)=\operatorname{supp}(Q_I)\subseteq S\setminus\{v\}.
\tag{16}
\]

This is the same global operator \(Q\), now with support at most \(m-1\);
its already-assumed identity \(\delta_{a,b}(Q)=0\) is unchanged. The
induction hypothesis therefore makes it scalar. This completes the
induction. \(\square\)

Equivalently, (16) contradicts the minimal-support characterization (1) at
the first nonempty peeling step. The induction presentation makes explicit
that no fixed finite volume or finite-box base case is being used.

**Corollary 7 [THEOREM — simultaneous centralizer].** With

\[
A=\sum_xX_x,\qquad B=\sum_{\langle x,y\rangle}Z_xZ_y
\]

understood through Lemma 1,

\[
\{Q\in\mathcal A_{\rm loc}:[Q,A]=[Q,B]=0\}=\mathbb C I.
\tag{17}
\]

*Proof.* Simultaneous commutation implies
\([Q,A+B]=0\). Apply Theorem 6 with \(a=b=1\). \(\square\)

The main theorem is strictly stronger than this corollary: it already rules
out compact \(Q\) commuting with any *one* nondegenerate linear combination
\(A+gB\), not only with the pair separately.

---

## 5. Exact finite checks (illustrations, not premises)

All calculations in this section are generated by
`experiments/e115_local_commutant.py`, recorded in
`results/integrability/local_commutant.json`, and independently reconstructed
by `tests/test_local_commutant.py`. They use the canonical ordered-Pauli basis
and exact integers/Fractions.

**[COMPUTATION — exhaustive private-neighbour geometry].** Every nonempty
subset of each listed box was tested using the full \(\mathbb Z^3\)
nearest-neighbour relation, not a truncated-box relation:

| box | sites | nonempty subsets | maximal-site checks | failures |
|---|---:|---:|---:|---:|
| \(2\times2\times2\) | 8 | 255 | 544 | 0 |
| \(3\times2\times2\) | 12 | 4,095 | 8,736 | 0 |
| \(2\times2\times3\) | 12 | 4,095 | 12,480 | 0 |
| \(3\times3\times1\) | 9 | 511 | 876 | 0 |
| \(3\times3\times2\) | 18 | 262,143 | 798,912 | 0 |
| **total** | — | **271,099** | **821,548** | **0** |

The producer also stores explicit nonmaximal choices where the conclusion
fails, so the maximal-first-coordinate condition is actually exercised. This
finite census illustrates Lemma 2; it is **not** its proof.

**[COMPUTATION — exact local coefficient maps].** For a one-site support and
for a bond support, the JSON stores every nonzero integer entry of the two
pencil maps

\[
\pi^{(w)}_{01}(a\operatorname{ad}_A+b\operatorname{ad}_B),
\qquad
\pi^{(v)}_{11}(a\operatorname{ad}_A+b\operatorname{ad}_B),
\]

including the full finite site legend used for each target coefficient. It
checks exactly that

\[
\pi^{(w)}_{01}\operatorname{ad}_A=0,
\qquad
\pi^{(w)}_{01}\operatorname{ad}_B=\operatorname{ad}_{Z_v},
\tag{18}
\]

and that the \(b\)-part of the second map vanishes on all \(a_v=0\) columns,
whereas the \(a\)-part gives (14). The verifier deliberately drops all
boundary bonds as a sabotage control; then (18) fails. It also checks that a
nonzero \(Z_v\)-slotted column has a nonzero \((1,1)_v\) image, excluding a
spurious local-basis cancellation.

**[COMPUTATION — representative exact kernels].** For every row below, the
full \(\mathbb Z^3\) bond set meeting the support was retained. The table
lists the exact-\(\mathbb Q\) rank for each of the six nonzero coefficient
pairs

\[
(a,b)\in\{(1,1),(1,2),(2,1),(3,5),(1,-1),(2,3)\};
\]

all six ranks agree within a row. The tagged simultaneous map
\((\operatorname{ad}_A,\operatorname{ad}_B)\) has the same displayed rank.
Thus every listed nullity is exactly one, with identity spanning the kernel.
Two prime-field reconstructions at
\(2\,147\,483\,647\) and \(2\,147\,483\,629\) agree with the exact ranks.

| support | \(|S|\) | columns \(4^{|S|}\) | bonds meeting \(S\) (internal/boundary) | exact rank | nullity |
|---|---:|---:|---:|---:|---:|
| one site | 1 | 4 | 6 (0/6) | 3 | 1 |
| bond | 2 | 16 | 11 (1/10) | 15 | 1 |
| length-3 path | 3 | 64 | 16 (2/14) | 63 | 1 |
| L tromino | 3 | 64 | 16 (2/14) | 63 | 1 |
| 3D claw | 4 | 256 | 21 (3/18) | 255 | 1 |
| plaquette | 4 | 256 | 20 (4/16) | 255 | 1 |
| degree-4 planar star | 5 | 1,024 | 26 (4/22) | 1,023 | 1 |
| 3D corner | 5 | 1,024 | 25 (5/20) | 1,023 | 1 |
| \(2\times3\) slab | 6 | 4,096 | 29 (7/22) | 4,095 | 1 |

**[COMPUTATION — why the boundary terms matter].** With only internal bonds,
the finite \(2\times2\) and \(2\times3\) graph pairs have the exact
simultaneous-commutant nullities 27 and 12, respectively, as certified in
`e52_full_commutant.py`. Restoring the \(\mathbb Z^3\) boundary bonds gives
nullity 1 for both embedded supports. The e115 producer cross-checks the
finite values modulo both primes and the standalone verifier rederives the
full-boundary kernels. This is an explanatory contrast, not an all-size
finite-volume classification.

None of the finite computations is used in the proof of Theorem 6. Lemmas
1–5, valid for an arbitrary finite support, are the complete proof.

---

## 6. Exact scope and unresolved routes

**[UNRESOLVED — excluded operator classes].** The theorem has the following
sharp boundaries.

1. It is only about **compact finite-support** \(Q\). Extensive translation
   sums of local densities, including the familiar type of charges in the 1D
   integrable chain, have infinite support and are not elements of
   \(\mathcal A_{\rm loc}\).
2. Quasilocal, weakly local, and other infinite-support operators are not
   classified. In particular, no statement is made about a limiting charge in
   a larger operator topology.
3. On a finite torus, \(P=\prod_xX_x\) has full torus support and commutes
   with both finite \(A\) and \(B\). It is not a compact observable on
   \(\mathbb Z^3\), so it is not a counterexample to (15) or (17).
4. Both coefficients matter. For \(b=0\), a single-site \(X_v\) commutes
   with the pure field; for \(a=0\), a single-site \(Z_v\) commutes with the
   pure diagonal interaction. These degenerate one-generator centralizers
   are intentionally not ruled out.
5. The theorem does not exclude nonlocal transformations, enlarged Hilbert
   spaces, or a charge commuting with a nonlinear transfer operator/product
   while failing to commute with every nondegenerate linear generator
   \(H_{a,b}\).
6. It makes no claim about transfer-operator spectra, thermodynamic
   integrability, or an exact solution of the 3D Ising model.

Thus Theorem 6 closes precisely the compact-local conserved-observable route
for every nonzero transverse field and Ising coupling; it does not convert
that local obstruction into a general non-integrability theorem.
