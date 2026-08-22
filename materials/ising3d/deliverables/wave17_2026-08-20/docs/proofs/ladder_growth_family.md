# Two-leg ladder family: exact obstructions to single-string induction

## 1. Result and scope

For the open two-leg ladder \(\Lambda_L\), let

\[
 A_L=\sum_{v\in\Lambda_L}X_v,
 \qquad
 B_L=\sum_{(u,v)\in E(\Lambda_L)}Z_uZ_v,
 \qquad
 \mathfrak g_L=\langle A_L,B_L\rangle_{\mathrm{Lie}}.
\]

Use the depth filtration from `proofs/ladder_growth.md`, so that \(D_{2L}\) is the dimension of the
span of generator-left-nested words of depth at most \(2L\).

**[THEOREM — single-Pauli membership obstruction].** For every \(L\geq2\), any two single Pauli
words that individually belong to \(\mathfrak g_L\) commute. Consequently, no family obtained by
iterated commutators of a single-Pauli seed with single-Pauli elements of \(\mathfrak g_L\) can be
the desired \(2^L\)-element injection: every nonempty iterated commutator in that class is zero.
This conclusion is independent of the bracket depth and of the total order used to select leading
monomials.

**[LEMMA — fixed-support seed obstruction].** More generally, even without imposing Lie-algebra
membership, a single Pauli seed supported on a fixed number \(b\) of rungs and translated
single-rung Pauli controls can produce at most \(2^b\) nonzero nested commutators. A corner seed
therefore cannot encode all subsets of arbitrarily many rungs.

**[LEMMA — orbit-symmetrized reflection obstruction].** The obvious attempt to repair membership
by replacing each rung control with its row/reflection orbit sum also fails. For commuting depth-one
orbit controls, reflection gives the exact vector identity

\[
 w_S=w_{\rho(S)}.
\]

Thus there are at most

\[
 \frac{2^L+2^{\lceil L/2\rceil}}2<2^L
\]

distinct outputs for \(L\geq2\), regardless of the monomial order.

**[COMPUTATION]** An exhaustive exact census of 96 natural product-lex orders finds a repaired order
that has enough raw leaders through \(L=4\), but it has only 26 leaders at \(L=5\), where 32 are
required. The maximum over all 96 orders at \(L=5\) is 27. A separate exhaustive search of 3,825
fixed two-rung-seed/rung-control families at each \(L=2,3,4,5\) finds no family that survives past
\(L=2\). Exact certificates are in
`results/algebra_growth/ladder_induction.json`.

**[CONJECTURE / UNRESOLVED]** This work does **not** prove or refute
\(D_{2L}\geq2^L\) for all \(L\), and it does not prove or refute the corresponding lower bound for
the full algebra \(\mathfrak g_L\). It rules out the requested single-string construction and its
most direct symmetry-orbit repair. A proof based on nonlocal, multi-string invariant elements remains
possible. No finite verification is promoted to an all-length claim.

## 2. Exact Pauli convention

Write

\[
 Q_{(a\mid b)}=X^aZ^b,
 \qquad (a\mid b)\in\mathbb F_2^{4L}.
\]

For Pauli codes \(p,q\),

\[
 \frac12[Q_p,Q_q]=
 \begin{cases}
 0,&\langle p,q\rangle=0,\\
 \epsilon(p,q)Q_{p+q},&\langle p,q\rangle=1,
 \end{cases}
 \qquad \epsilon(p,q)\in\{+1,-1\}.
\]

All computations below retain the exact integer coefficients in this basis. Dividing every
commutator by 2 changes only nonzero scalar factors, never support, membership, or linear
independence.

## 3. Attempt 1: the Pauli signature works, but not inside the algebra

There is a particularly simple all-length Pauli injection if membership is ignored. On rung \(i\),
put

\[
 P_L=\bigotimes_{i=1}^L (ZI)_i,
 \qquad C_i=(XX)_i.
\]

The \(C_i\) commute with one another and each anticommutes with \(P_L\). Therefore, for
\(S\subseteq\{1,\ldots,L\}\),

\[
 w_S=\operatorname{ad}_{C_{s_1}}\cdots
     \operatorname{ad}_{C_{s_k}}(P_L)
\]

is a single nonzero Pauli word. In the normalized bracket convention its rung signature is

\[
 (w_S)_i=
 \begin{cases}
 ZI,&i\notin S,\\
 YX,&i\in S.
 \end{cases}
\]

The unnormalized coefficient has magnitude \(2^{|S|}\). Hence \(S\mapsto w_S\) is injective for
all \(L\). `experiments/e46_ladder_induction.py` verifies all \(2^L\) words exactly for
\(2\leq L\leq7\).

This is **not** a theorem about \(\mathfrak g_L\): \(P_L\notin\mathfrak g_L\), and the required
individual off-centre \(C_i\notin\mathfrak g_L\). The failure is proved, not inferred from a finite
search, in Section 5.

## 4. Locality goes the wrong way for a corner seed

**[LEMMA 1 — fixed-support lemma].** Let \(P\) be a single Pauli word, and let \(C_i\) be a single
Pauli word supported only on rung \(i\). The \(C_i\) commute pairwise. If \(C_j\) commutes with
\(P\), then every nested commutator using \(C_j\) once and any collection of the other \(C_i\)'s is
zero. In particular, if rung \(j\) is outside the support of \(P\), every such word is zero.

**Proof.** Disjoint rung controls commute, so

\[
 [\operatorname{ad}_{C_i},\operatorname{ad}_{C_j}]
 =\operatorname{ad}_{[C_i,C_j]}=0.
\]

Move \(\operatorname{ad}_{C_j}\) through the other adjoints until it acts directly on \(P\). It
then gives \([C_j,P]=0\). If rung \(j\) is outside the support of \(P\), disjointness already implies
this commutation. \(\square\)

**Corollary.** If the seed touches \(b\) rungs, only subsets of those \(b\) rungs can have nonzero
outputs, so there are at most \(2^b\) nonzero outputs. In particular a fixed corner seed cannot
produce \(2^L\) nonzero strings for unbounded \(L\).

This corrects the tempting but incomplete statement that “a commutator of Pauli strings is another
Pauli string.” It is another Pauli string only in the anticommuting case; for disjoint support it is
zero. A rung-local commutator does not attach a new remote rung to the seed.

**[COMPUTATION]** The exact search used every nonidentity seed on the first two rungs
(\(16^2-1=255\) choices) and every nonidentity translated rung control (15 choices), for 3,825
families at each length:

| \(L\) | fully injective families | maximum distinct nonzero outputs | required |
|---:|---:|---:|---:|
| 2 | 960 | 4 | 4 |
| 3 | 0 | 4 | 8 |
| 4 | 0 | 4 | 16 |
| 5 | 0 | 4 | 32 |

The canonical edge seed \(P=(ZI)_1(ZI)_2\), which is a top-row horizontal \(ZZ\) edge, with
\(C_i=(XX)_i\), realizes all four outputs at \(L=2\). For every \(L\geq3\), the labels
\(S=\{3\}\) and \(S=\{1,3\}\) both give the zero vector. This is an exact collision, not a
leading-term collision.

## 5. Why isolated single strings cannot be elements of the ladder algebra

Let \(\tau\) exchange the two rows of the ladder, and let \(\rho\) reflect its columns. Both are
graph automorphisms. Their permutation operators fix the two generators:

\[
 \tau(A_L)=A_L,\quad \tau(B_L)=B_L,
 \qquad
 \rho(A_L)=A_L,\quad \rho(B_L)=B_L.
\]

Hence they fix every Lie polynomial in \(A_L,B_L\), and therefore every element of
\(\mathfrak g_L\).

**[THEOREM 2 — single-Pauli no-go].** If a nonzero scalar multiple of a single Pauli word \(Q\)
lies in \(\mathfrak g_L\), then \(Q\) has identical top and bottom labels on every rung. Any two
Pauli words with this property commute.

**Proof.** Pauli words form a linearly independent basis and row permutation introduces no phase.
Since every algebra element is fixed by \(\tau\), a single basis word occurring alone must satisfy
\(\tau(Q)=Q\). Thus its labels on rung \(i\) have the form \(q_iq_i\), with
\(q_i\in\{I,X,Y,Z\}\).

Take two such words \(Q\) and \(R\). At a given rung, the binary symplectic contribution from the
top site is repeated identically at the bottom site. Their total contribution is therefore twice a
bit and vanishes modulo 2. Summing over rungs gives \(\langle Q,R\rangle=0\), so \([Q,R]=0\).
\(\square\)

This theorem directly blocks the proposed “isolate single-string elements first” strategy. If a
single-string seed and every single-string \(C_i\) were actually isolated by nested commutators and
linear combinations, they would be elements of \(\mathfrak g_L\), hence they would commute and the
first inner bracket would vanish.

Longitudinal reflection supplies a second, independent membership check for the formal injection in
Section 3. The seed \((ZI)^L\) is not fixed by row exchange. Although \((XX)_i\) is row-fixed,
reflection sends it to \((XX)_{L+1-i}\); an individual off-centre rung word is not fixed and cannot
belong to \(\mathfrak g_L\).

**[COMPUTATION]** Exhaustive enumeration gives:

| \(L\) | row-fixed Pauli words | row-and-reflection-fixed words | anticommuting pairs among row-fixed words |
|---:|---:|---:|---:|
| 2 | 16 | 4 | 0 |
| 3 | 64 | 16 | 0 |
| 4 | 256 | 16 | 0 |
| 5 | 1024 | 64 | 0 |

The counts are the proved formulas \(4^L\) and \(4^{\lceil L/2\rceil}\); the enumeration is a finite
implementation check of the symmetry maps and symplectic convention.

A separate search over all 225 repeated-local-pattern pairs found 120 algebra-free Pauli injections
(the 120 anticommuting ordered pairs in the nonzero two-qubit Pauli space), valid for every tested
\(L=2,3,4,5\). Exactly zero of those pairs pass even the necessary row-swap condition for both the
seed pattern and the control pattern. Thus the Pauli bookkeeping is not the missing ingredient;
membership is.

## 6. Attempt 2: orbit symmetrization creates intrinsic reflection collisions

One may replace a forbidden local word by a symmetry orbit sum. Let a translated depth-one rung
pattern be symmetrized under row exchange and longitudinal reflection, and call the resulting
control \(C_i\). Then

\[
 C_i=C_{L+1-i}.
\]

For the rung Pauli orbit sums tested here, all \(C_i\)'s commute: summands on different rung orbits
are disjoint, while row-swapped summands on one rung differ on zero or two sites and also commute.
Therefore the adjoints commute.

**[LEMMA 3 — reflection collision].** For any fixed seed \(P\) and any such commuting control family,

\[
 w_S=\prod_{i\in S}\operatorname{ad}_{C_i}(P)
     =\prod_{i\in S}\operatorname{ad}_{C_{L+1-i}}(P)
     =w_{\rho(S)}.
\]

The number of reflection orbits of subsets is, by the two-element Burnside count,

\[
 N_L=\frac{2^L+2^{\lceil L/2\rceil}}2.
\]

For \(L\geq2\), \(N_L<2^L\), so this family cannot be injective. The equality is between complete
vectors before any leading term is selected, and hence no total order can repair it. \(\square\)

**[COMPUTATION]** The script enumerated three corner-seed orbit types and all 15 nonidentity
single-rung control patterns. It checked \(w_S=w_{\rho(S)}\) for every subset at \(L=2,3,4,5\), a
total of 180, 360, 720, and 1,440 exact identities respectively. A nonzero collision occurs uniformly:
use the orbit sum of \(X\) on the four corners as seed and the row/reflection orbit sum of \(IY\) as
control. Then

\[
 w_{\{1\}}=w_{\{L\}}\ne0,
\]

with support equal to the four boundary strings having \(ZI\) or \(IZ\) on one end rung and identity
elsewhere (each coefficient is \(-1\) in the normalized convention).

The lemma covers the natural translation-covariant depth-one orbit repair. It does not cover
noncommuting, nonlocal, multi-string controls.

## 7. Repaired leading orders: a finite improvement and its exact failure

The previous order was right-boundary-first, top-before-bottom, with local order
\(X>Y>I>Z\). The repaired order changes only the last two local priorities:

\[
 X>Y>Z>I.
\]

The experiment enumerated every nonzero generator-left-nested \(A/B\) word of depth at most \(2L\)
and every order in the class

\[
 2\text{ longitudinal directions}\times
 2\text{ within-rung directions}\times
 24\text{ local Pauli orders}=96.
\]

**[COMPUTATION]** The exact capacities are:

| \(L\) | nonzero eligible words | old order | repaired order | maximum over 96 | required \(2^L\) |
|---:|---:|---:|---:|---:|---:|
| 2 | 16 | 5 | 6 | 6 | 4 |
| 3 | 64 | 9 | 9 | 10 | 8 |
| 4 | 256 | 15 | 16 | 16 | 16 |
| 5 | 1024 | 27 | 26 | 27 | 32 |

Thus the changed order repairs the *capacity* failure at \(L=4\): the artifact stores one actual
eligible bracket word for each of its 16 leaders. It does not repair the old recursive family. Under
both the old and repaired orders,

\[
 0101\mapsto\mathtt{BBABBA},
 \qquad
 0110\mapsto\mathtt{ABBBBA}
\]

still share leader

\[
 \mathtt{II\ II\ ZI\ XY}
\]

(code 33928, coefficient 8 in each vector). The vectors are nevertheless independent: the exact
coefficient minor on supports `YI XI ZI II` and `ZI XI YI II` is

\[
 \det\begin{pmatrix}14&20\\20&20\end{pmatrix}=-120\ne0.
\]

At \(L=5\), the repaired-order recursive labels 01001 and 01010, words
`BBABABBA` and `ABBBABBA`, share leader `II II II YZ XX` (code 270872), with coefficients
\(-28\) and \(-16\). An exact two-support minor has determinant \(-2268\), again proving that this is
a leader collision rather than proportionality.

More strongly, the repaired order has only 26 leaders among **all** eligible depth-10 words, so no
choice of 32 words can have distinct leaders under it. None of the 96 product-lex orders passes the
capacity requirement simultaneously for \(L=2,3,4,5\). This is a finite obstruction to this order
class, not a theorem about every conceivable total order.

## 8. Reproducibility and remaining route

Run from the repository root:

```sh
.venv/bin/python experiments/e46_ladder_induction.py
.venv/bin/python tests/test_ladder_induction.py
```

The standalone test reimplements the Pauli encoding, commutators, 96-order census, fixed-seed search,
symmetry census, and orbit collision; it does not obtain the claims by merely reading the JSON.
Both programs print a final `PASS`.

The exact machine-readable certificate is
`results/algebra_growth/ladder_induction.json`. Its top-level status is `OBSTRUCTION_LEMMA` and the
target all-length bound is explicitly `UNRESOLVED`.

A future positive proof must leave the class ruled out here. In particular it must use multi-string
invariant elements (not isolated single Pauli words), and it must control cancellation and membership
inside \(\mathfrak g_L\) symbolically. Finite modular independence alone remains insufficient.
