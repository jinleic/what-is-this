# Two-slice observability for the ladder cyclic space and the first bonding-truncation kernel

Artifacts:

- `experiments/e230_wlaw_rank_mechanism.py`
- `results/algebra_growth/wlaw_rank_mechanism.json`
- `tests/test_wlaw_rank_mechanism.py` (independent; no producer import)

Producer command, from the repository root:

```sh
.venv/bin/python experiments/e230_wlaw_rank_mechanism.py
```

The independent verifier is intentionally separate:

```sh
.venv/bin/python tests/test_wlaw_rank_mechanism.py
```

## 0. Verdict and scope

Let

\[
 U_L=\operatorname{alg}(A_L,B_L)\psi\subset K_L,
 \qquad W_L=U_L^\perp\subset K_L,
\]

in the open \(2\times L\) ladder conventions of `proofs/wlaw.md`.  In particle
sector \(2m\), write

\[
 \Delta(x)=\bigl|n_{\rm top}(x)-n_{\rm bottom}(x)\bigr|,
\]

and split the invariant configuration-orbit space orthogonally as

\[
 K_L^{(2m)}=E_{L,m}\oplus H_{L,m},
 \quad E_{L,m}=\operatorname{span}\{\Delta=0,2\},
 \quad H_{L,m}=\operatorname{span}\{\Delta\ge4\}.
\tag{1}
\]

The finite exact result is:

> **[COMPUTATION] — exact over \(\mathbb Q\).** For every \(L=3,\ldots,8\)
> and every \(m=0,\ldots,L\), coordinate restriction
> \(R_{\le2}:U_L^{(2m)}\to E_{L,m}\) is injective.  Consequently its exact
> rank is the already certified exact rank of \(U_L^{(2m)}\).  In particular,
> the exceptional rows have ranks \(20\) at \((L,m)=(4,2)\) and \(2484\) at
> \((8,4)\).

This is a new direct rank mechanism beyond the wave-20 natural-map theorem:
only the balanced and minimally unbalanced physical slices are needed to
observe every vector in the exact finite cyclic spaces.  It is not an all-\(L\)
rank theorem.

A second result closes the most economical rungwise bonding/antibonding
truncation:

> **[THEOREM].** In the unnormalised rung basis
> \(b=|10\rangle+|01\rangle\), \(a=|10\rangle-|01\rangle\), the projection
> retaining at most two \(a\)-factors first fails at \((L,m)=(4,2)\).  Its
> exact kernel on \(U_4^{(4)}\) is the line spanned by
> \(a_0a_1a_2a_3\), and its exact rank is \(19\), not \(20\).  Hence every
> linear evaluation that factors through this truncation is non-injective.

No sequence was fitted, no signed-monomial reflection map was introduced, no
new proof of the existing invariant-theory dimension identity is restated, no
\(L=9\) or \(L=10\) closure was run, and \(K_c\) is absent.

---

## 1. [LEMMA] The omitted-slice duality criterion

The mechanism is a finite-dimensional duality identity rather than a proposed
intertwiner.

**Lemma 1 (observability by the annihilator).**  Let \(K=E\oplus H\) be an
orthogonal coordinate decomposition over \(\mathbb Q\), let \(U\le K\), and
let \(W=U^\perp\).  Define

\[
 \Theta_H:W\longrightarrow H^*,\qquad
 \Theta_H(w)(h)=\langle w,h\rangle.
\tag{2}
\]

If \(R_E:K\to E\) is coordinate restriction, then

\[
 \dim\ker(R_E|_U)=\dim H-\operatorname{rank}_{\mathbb Q}\Theta_H,
\tag{3}
\]

and therefore

\[
 R_E|_U\text{ is injective}
 \quad\Longleftrightarrow\quad
 \operatorname{rank}_{\mathbb Q}\Theta_H=\dim H.
\tag{4}
\]

Equivalently,

\[
 \operatorname{rank}(R_E|_U)
 =\dim U-\dim H+\operatorname{rank}\Theta_H.
\tag{5}
\]

*Proof.*  The annihilator in \(H\) of \(\Theta_H(W)\) is

\[
 \{h\in H:\langle w,h\rangle=0\ \forall w\in W\}
 =H\cap W^\perp=H\cap U.
\]

The pairing on the coordinate space \(H\) is non-degenerate, so
rank-nullity for the annihilator gives
\(\dim(H\cap U)=\dim H-\operatorname{rank}\Theta_H\).  But
\(H\cap U=\ker(R_E|_U)\), proving (3)--(5).  ∎

For the ladder orbit-sum basis \(e_o=\sum_{x\in o}|x\rangle\), the
configuration pairing is diagonal with nonzero weights \(|o|\in\{1,2,4\}\).
Thus the matrix of (2) is obtained by restricting every exact annihilator
basis vector to omitted orbit representatives and multiplying column \(o\) by
\(|o|\).

This lemma is general; the full-column-rank statement used below is finite and
computed only for \(L\le8\).

---

## 2. [COMPUTATION] Exact two-sided ranks for all rows \(L=3,\ldots,8\)

The source `results/ladder/l8_saturation.json` contains complete integer bases
of \(W_L=U_L^\perp\) for \(L=3,\ldots,8\).  Those bases were already validated
over \(\mathbb Q\): homogeneous, independent, vacuum-orthogonal, and invariant
under \(A_L\) and \(B_L\); their dimensions close the exact cyclic-space
sandwiches.  The producer pins the source by SHA-256 and rechecks all of those
contracts before using it.

For each row, the producer independently enumerates the
\(\langle\tau,\rho\rangle\)-orbits, forms the weighted restriction matrix
\(M_{L,m}\) of \(W_L^{(2m)}\) on \(H_{L,m}\), and ranks it at the two primes
\(999983\) and \(1000003\).  The logical sandwich is always

\[
 \operatorname{rank}_{\mathbb F_p}M_{L,m}
 \le \operatorname{rank}_{\mathbb Q}M_{L,m}
 \le \dim H_{L,m}.
\tag{6}
\]

Every modular rank in the artifact equals the right-hand coordinate ceiling.
Thus every displayed rank is exact over \(\mathbb Q\); agreement of the two
primes is a control, not a substitute for the ceiling.

The following table gives one representative of each particle-hole pair.
Here \(K=\dim K_L^{(2m)}\), \(E=\dim E_{L,m}\), \(H=\dim H_{L,m}\),
\(W=\dim W_L^{(2m)}\), and the column `rank` is the exact rank of
\(\Theta_H\).  Rows beyond the middle are exact mirrors and are present in the
machine artifact, so all 39 sectors—not just this compressed table—are
covered.

| \(L\) | \(m\) | \(K\) | \(E\) | \(H\) | \(W\) | \(\operatorname{rank}\Theta_H\) | \(\dim U=\operatorname{rank}R_{\le2}(U)\) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 0 | 1 | 1 | 0 | 0 | 0 | 1 |
| 3 | 1 | 6 | 6 | 0 | 0 | 0 | 6 |
| 4 | 0 | 1 | 1 | 0 | 0 | 0 | 1 |
| 4 | 1 | 10 | 10 | 0 | 0 | 0 | 10 |
| 4 | 2 | 22 | 21 | 1 | 2 | 1 | 20 |
| 5 | 0 | 1 | 1 | 0 | 0 | 0 | 1 |
| 5 | 1 | 15 | 15 | 0 | 0 | 0 | 15 |
| 5 | 2 | 60 | 57 | 3 | 5 | 3 | 55 |
| 6 | 0 | 1 | 1 | 0 | 0 | 0 | 1 |
| 6 | 1 | 21 | 21 | 0 | 0 | 0 | 21 |
| 6 | 2 | 135 | 126 | 9 | 15 | 9 | 120 |
| 6 | 3 | 246 | 227 | 19 | 36 | 19 | 210 |
| 7 | 0 | 1 | 1 | 0 | 0 | 0 | 1 |
| 7 | 1 | 28 | 28 | 0 | 0 | 0 | 28 |
| 7 | 2 | 266 | 247 | 19 | 35 | 19 | 231 |
| 7 | 3 | 777 | 698 | 79 | 147 | 79 | 630 |
| 8 | 0 | 1 | 1 | 0 | 0 | 0 | 1 |
| 8 | 1 | 36 | 36 | 0 | 0 | 0 | 36 |
| 8 | 2 | 476 | 438 | 38 | 70 | 38 | 406 |
| 8 | 3 | 2044 | 1804 | 240 | 448 | 240 | 1596 |
| 8 | 4 | 3270 | 2837 | 433 | 786 | 433 | 2484 |

Since the rank column equals \(H\) everywhere, Lemma 1 proves injectivity.
Both rank directions are then exact on the audited range:

1. \(\operatorname{rank}R_{\le2}(U)\le\dim U\) because restriction cannot
   increase rank;
2. \(\ker(R_{\le2}|_U)=0\) by (3) and the exact full omitted-column rank, so
   \(\operatorname{rank}R_{\le2}(U)\ge\dim U\).

The value of \(\dim U\) in the last column is the prior characteristic-zero
sandwich from the exact annihilator source.  The comparison with
\(\binom{\binom Lm+1}{2}-\delta\) is recorded only after these ranks are
obtained and selects no certificate.

### Exceptional rows

- At \((4,2)\), \(H\) is one-dimensional and the weighted annihilator
  restriction has exact rank \(1\), so the two-slice restriction has exact
  rank \(20\).
- At \((8,4)\), \(H\) has dimension \(433\).  Both finite-field ranks are
  \(433\), hence the rational rank is exactly \(433\), and the retained
  restriction has exact rank \(2484\).

These are the two defect sectors requested explicitly; neither is tuned or
selected using any numerical critical coupling.

---

## 3. The one missing all-\(L\) lemma—stated, not assumed

The finite mechanism suggests a precise next theorem rather than an
extrapolated sequence.

> **[UNRESOLVED] — two-slice exact-sequence lemma.**  For every \(L,m\), the
> weighted map
> \[
>   \Theta_H:W_L^{(2m)}\longrightarrow H_{L,m}^*
> \]
> is onto, and
> \[
>   \dim\ker\Theta_H
>   =\dim E_{L,m}
>    -\left(\binom{\binom Lm+1}{2}-\delta(L,m)\right).
> \tag{7}
> \]

This single exact-sequence lemma would give both rank inclusions.  Indeed,
(7) and surjectivity would imply

\[
 \dim W=\dim H+\dim E-
 \left(\binom{\binom Lm+1}{2}-\delta\right)
 =\dim K-
 \left(\binom{\binom Lm+1}{2}-\delta\right),
\]

hence the desired value of \(\dim U=\dim K-\dim W\); surjectivity and Lemma 1
would simultaneously make \(R_{\le2}|_U\) injective.  The artifact verifies
this exact sequence only for \(L=3,\ldots,8\).  It is **not** asserted for
\(L\ge9\), and the finite calculation is not an all-\(L\) proof.

---

## 4. [LEMMA] What the rungwise bonding/antibonding transform actually does

On one rung use the unnormalised basis

\[
 |0\rangle=|00\rangle,\qquad
 |b\rangle=|10\rangle+|01\rangle,\qquad
 |a\rangle=|10\rangle-|01\rangle,\qquad
 |d\rangle=|11\rangle.
\tag{8}
\]

Leg swap is diagonal: it fixes \(0,b,d\) and negates \(a\).  Therefore the
\(\tau\)-invariant space has an even number of \(a\)-factors.

Let \(R=X_{\rm top}X_{\rm bottom}\) be a rung-edge toggle and let
\(H=X_{{\rm top},r}X_{{\rm top},r+1}+
X_{{\rm bottom},r}X_{{\rm bottom},r+1}\) be the sum of the two horizontal
edge toggles between adjacent rungs.  Exact change of basis gives

\[
\begin{array}{c|rrrr}
 v&0&b&a&d\\ \hline
 2Rv&2d&2b&-2a&2\,0,
\end{array}
\tag{9}
\]

and, among the cases that change the number of \(a\)-factors,

\[
\begin{aligned}
 2H|00\rangle&=|bb\rangle+|aa\rangle,\\
 2H|0d\rangle&=|bb\rangle-|aa\rangle,\\
 2H|d0\rangle&=|bb\rangle-|aa\rangle,\\
 2H|dd\rangle&=|bb\rangle+|aa\rangle,\\
 2H|bb\rangle&=4(|00\rangle+|0d\rangle+|d0\rangle+|dd\rangle),\\
 2H|aa\rangle&=4(|00\rangle-|0d\rangle-|d0\rangle+|dd\rangle),\\
 2H|ba\rangle&=2H|ab\rangle=0.
\end{aligned}
\tag{10}
\]

The remaining eight local sources preserve the number of \(a\)-factors.
The machine artifact contains and checks all 16 horizontal sources and all
four rung sources, not merely the displayed representatives.

> **Lemma 2 (integral block tridiagonality).**  In basis (8), \(2B_L\) has
> integer entries and changes the number of \(a\)-factors only by
> \(0,\pm2\).  After splitting by particle grade, the same is true separately
> for \(2D,2F,2D^+\).  Thus the even-\(a\) block is block-tridiagonal in
> \(a\)-count.

*Proof.*  Equation (9) preserves the \(a\)-count.  The exhaustive two-rung
table underlying (10) has changes only \(0,\pm2\).  Particle number in every
entry changes by \(0,\pm2\), so selecting the creation, hopping, or annihilation
grade cannot introduce another \(a\)-change.  Scaling by two clears the only
basis-change denominators and is invertible over \(\mathbb Q\).  ∎

This is a genuine non-monomial structural simplification, but it does not by
itself prove the rank law: the next section shows that its smallest tempting
truncation loses a real cyclic vector.

---

## 5. [THEOREM] Exact first kernel of the \(\#a\le2\) truncation

Let \(Q_{\le2}\) discard all transformed coordinates containing four or more
\(a\)-factors.  Call an evaluation a **factor-through-\(Q_{\le2}\) map** if it
has the form \(\Phi=T\circ Q_{\le2}\) for an arbitrary linear map \(T\).
This names the entire map class being closed, not one selected matrix.

For \(L=3\), no rung word can contain four \(a\)-factors, so the truncation is
the identity in every sector.  For \(L=4\), it remains the identity in particle
sectors \(0\) and \(2\).  The lexicographically first possible omitted
coordinate is therefore \((L,m)=(4,2)\), where the omitted
\(\tau\)-plus, \(\rho\)-invariant subspace is the single line
\(\mathbb Q\,a_0a_1a_2a_3\).

Define its physical expansion

\[
 v=\bigotimes_{r=0}^{3}
 \left(|10\rangle_r-|01\rangle_r\right).
\tag{11}
\]

It has 16 nonzero coefficients, all \(\pm1\), and exactly four particles.
Moreover \(\tau v=(-1)^4v=v\), while rung reversal merely permutes its four
identical factors, so \(v\in K_4^{(4)}\).

The exact annihilator certificate has two basis vectors in
\(W_4^{(4)}\).  Direct configuration pairing, including orbit sizes, gives

\[
 \langle w_1,v\rangle=0,
 \qquad
 \langle w_2,v\rangle=0.
\tag{12}
\]

The complete basis and all 16 coefficients of (11) are stored in the artifact;
the independent verifier rebuilds them in a different site encoding.  Since
the cited characteristic-zero certificate proves
\(W_4^{(4)}=(U_4^{(4)})^\perp\), (12) implies
\(v\in U_4^{(4)}\).  It is nonzero and spans the entire omitted line.  Hence

\[
 \ker(Q_{\le2}|_{U_4^{(4)}})=\mathbb Qv,
 \qquad
 \operatorname{rank}_{\mathbb Q}Q_{\le2}(U_4^{(4)})=20-1=19.
\tag{13}
\]

Finally, every \(\Phi=T\circ Q_{\le2}\) kills \(v\), independently of the
choice of \(T\).  Thus the whole named truncation class fails, and (11) is the
lexicographically first exact kernel witness.  This obstruction concerns the
route, not the W-law itself.

---

## 6. Claim ledger and limitations

| claim | status |
|---|---|
| omitted-slice duality criterion (3)--(5) | **[LEMMA]**, all finite-dimensional coordinate splittings |
| \(R_{\le2}:U_L^{(2m)}\to E_{L,m}\) injective for every sector \(L=3,\ldots,8\) | **[COMPUTATION]**, exact over \(\mathbb Q\) by modular lower bound = omitted-column ceiling |
| exact exceptional retained ranks \(20\) at \((4,2)\), \(2484\) at \((8,4)\) | **[COMPUTATION]**, exact finite rows |
| integral, even-\(a\), block-tridiagonal local action of \(D,F,D^+\) | **[LEMMA]**, exhaustive local table |
| first failure of \(Q_{\le2}\), kernel \(\mathbb Q a_0a_1a_2a_3\), every factor-through map non-injective | **[THEOREM]** |
| two-slice exact-sequence lemma (7), hence the rank law for all \(L\) | **[UNRESOLVED]**; verified only for \(L=3,\ldots,8\) |
| any new rank at \(L\ge9\), any \(L=10\) falsification, or any use of \(K_c\) | **[UNRESOLVED]**; not computed and not used |

The result is deliberately a structural finite closure plus one exact failed
map class.  It neither promotes modular evidence without a rational ceiling
nor converts the exact \(L\le8\) mechanism into an all-size theorem.

---

## 7. [COMPUTATION] Integral exact-sequence certificate through \(L=9\)

Wave 22 adds:

- `experiments/e237_wlaw_exact_sequence.py`;
- `results/algebra_growth/wlaw_exact_sequence.json`;
- `tests/test_wlaw_exact_sequence.py` (clean-room; it does not import the
  producer).

Run, from the repository root,

```sh
.venv/bin/python experiments/e237_wlaw_exact_sequence.py
.venv/bin/python tests/test_wlaw_exact_sequence.py
```

This update uses the pre-existing exact annihilator bases through \(L=9\);
it does not run a ladder closure.  It strengthens the finite rank certificate
from a modular sandwich to an explicit nonzero integral minor and settles the
first exact row beyond the range of Section 2.

### 7.1 [LEMMA] Leaf-plus-core integral minors

Let \(M\) be an integer matrix with column set \(H\).  Starting with all
columns active, repeatedly choose an unused source row having exactly one
active nonzero coordinate, record that coordinate as its pivot, and remove
the pivot column.  Suppose this selects \(s\) rows and leaves \(c\) columns.
In selection order, the chosen rows have a triangular \(s\times s\) block
with nonzero integer diagonal, and they vanish on the final \(c\) columns.
If the unused rows on those \(c\) columns contain a \(c\times c\) minor
\(C\) with \(\det C\ne0\), then \(M\) has an \((s+c)\)-minor of the form

\[
 \begin{pmatrix}
  T&0\\
  *&C
 \end{pmatrix},
 \qquad
 \det T=\prod_{j=1}^{s}d_j\ne0.
\tag{14}
\]

Thus \(\operatorname{rank}_{\mathbb Q}M\ge s+c\).  When \(s+c=|H|\),
the coordinate count gives the reverse inequality and \(M\) has full column
rank over \(\mathbb Q\).  This is field-independent: reduction modulo primes
is unnecessary for the implication.

The producer applies this lemma to the orbit-size-weighted matrix of
\(\Theta_H:W_L^{(2m)}\to H_{L,m}^*\).  Every leaf step is stored explicitly
as a source-basis index, pivot representative, and nonzero integer diagonal.
The residual minors are evaluated by fraction-free Bareiss elimination.
The verifier rebuilds the map in separate top-mask/bottom-mask coordinates,
replays every stored leaf, repeats peeling with the opposite source-row
order, and evaluates the core determinants by exact
`fractions.Fraction` Gaussian elimination.

### 7.2 [COMPUTATION] The first exact row beyond \(L=8\) is positive

For \(L=9\), one representative of each particle-hole pair is:

| \(m\) | \(\dim K_9^{(2m)}\) | \(\dim E_{9,m}\) | \(\dim H_{9,m}\) | \(\dim W_9^{(2m)}\) | leaf pivots | core | \(\operatorname{rank}_{\mathbb Q}\Theta_H\) | \(\dim U_9^{(2m)}\) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 |
| 1 | 45 | 45 | 0 | 0 | 0 | 0 | 0 | 45 |
| 2 | 792 | 726 | 66 | 126 | 66 | 0 | 66 | 666 |
| 3 | 4704 | 4090 | 614 | 1134 | 614 | 0 | 614 | 3570 |
| 4 | 11034 | 9345 | 1689 | 3033 | 1556 | 133 | 1689 | 8001 |

The rows \(m=5,\ldots,9\) are the exact particle-hole mirrors and are present
separately in the artifact.  At \(m=4\) and \(m=5\), the unpeeled matrix is
already square, \(133\times133\), and both exact determinants are

\[
 \det C_{9,4}=\det C_{9,5}
 =-2^{288}3^5
 =-120849546447578154043757886299677464144410844643755495338032746385037757774045610254532608.
\tag{15}
\]

Consequently (14) gives rank \(1689\), exactly the omitted-coordinate
ceiling.  Lemma 1 then gives

\[
 \ker\!\left(R_{\le2}:U_9^{(2m)}\to E_{9,m}\right)=0
 \quad(0\le m\le9).
\tag{16}
\]

Hence there is no two-slice counterexample at the first exact size beyond
\(L=8\).  The exact retained ranks are

\[
 1,45,666,3570,8001,8001,3570,666,45,1.
\tag{17}
\]

These values come from the cited characteristic-zero annihilator sandwich;
they are not fitted or used to select a certificate.  At each of the two
control primes \(999983\) and \(1000003\), the independently reconstructed
rank also meets the same rational coordinate ceiling.  Those modular ranks
are controls only: the integral minor is already a rational lower bound.

### 7.3 [COMPUTATION][UNRESOLVED] Exact frontier of the raw leaf route

The leaf construction tests a precise direct-triangular method class: it may
use the stored exact annihilator source vectors in any successive singleton
order, but it may not first take rational row combinations.  That method does
not close every row.

Its first residual core is at \((L,m)=(8,4)\): \(416\) leaf pivots leave
\(17\) columns.  There are \(18\) active source rows, of residual degrees
\(2^3,3^{10},4^5\) (degree \(d\) with multiplicity shown as an exponent);
a selected \(17\times17\) minor has

\[
 \det C_{8,4}=-3\cdot2^{32}.
\tag{18}
\]

At \((9,4)\) and \((9,5)\), \(1556\) leaves leave \(133\) columns, and all
\(133\) surviving active rows have degree at least two (degrees
\(2^5,3^{36},4^{20},5^{10},6^{27},7^{30},8^5\)).  Thus no raw singleton
pivot remains; (15), a genuinely coupled block determinant, is required.
The opposite-order verifier obtains exactly the same cores.

This is an exact limitation of the named source-vector leaf method, not a
counterexample to row-combined triangular families and not a failure of
two-slice injectivity.  It narrows a successor proof to a block-core theorem
(or a structural row-combination/Koszul construction); merely asserting one
raw pivot per omitted coordinate is already false on the certified range.

### 7.4 Updated scope ledger

| claim | status |
|---|---|
| leaf-plus-core criterion (14) | **[LEMMA]**, every finite integer matrix |
| \(R_{\le2}:U_L^{(2m)}\to E_{L,m}\) injective in all sectors \(L=3,\ldots,9\) | **[COMPUTATION]**, exact over \(\mathbb Q\) by explicit integral minors meeting the coordinate ceiling |
| \(L=9\) core determinants (15) and retained ranks (17) | **[COMPUTATION]**, exact integers/rational ranks |
| raw uncombined-source singleton peeling closes every omitted coordinate | **[REFUTED METHOD]**; exact residual cores begin at \((8,4)\) |
| all-\(L\) two-slice exact-sequence lemma | **[UNRESOLVED]**; no exact annihilator source or rank is claimed for \(L\ge10\) |
| first counterexample beyond \(L=9\) | **[UNRESOLVED]**; no larger closure was run |

This appended section supersedes only the old finite frontier \(L\le8\) in
the ledger above.  It does not promote the finite integral certificates to an
all-size theorem, and it uses neither a sequence fit nor \(K_c\).

**Confluence note.**  The residual raw leaf core is independent of the order
of singleton choices.  If a row has active support \(\{c\}\), deleting other
columns either leaves \(\{c\}\) or deletes \(c\) first; it can never make a
different column the unique survivor of that same row while \(c\) remains.
Thus singleton deletions commute, and every fair order reaches the same
maximal fixed core.  The producer's increasing-order and verifier's
decreasing-order reconstructions are exact controls of this general argument.
