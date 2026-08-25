# W-law support defects: the finite `17=16+1` and `133=128+5` cores

Artifacts:

- `experiments/e249_wlaw_support_defect.py`
- `results/algebra_growth/wlaw_support_defect.json`
- `tests/test_wlaw_support_defect.py` (clean-room; it imports neither producer)
- exact predecessors `e230_wlaw_rank_mechanism.py` and `e237_wlaw_exact_sequence.py`

The producer and verifier use only exact integer arithmetic and Smith normal forms over
\(\mathbb Z\). They load the stored exact \(W_8\) and \(W_9\) bases; they do not launch an
\(L=10\) closure. Both enforce a process-CPU budget of 120 seconds and a 2 GiB RSS wall. On Darwin,
peak RSS is `mach_task_basic_info.resident_size_max`, not the platform-dependent `ru_maxrss`
shortcut.

## 1. Verdict and scope

There are two logically different results.

> **[LEMMA — all \(L\), elementary].** For any rung word of length \(L\) with \(2m\) occupied
> sites, let \(T,B,E,D\) be its top-only, bottom-only, empty, and double-occupied rung sets. Then
> \[
> |E|-|D|=L-2m,
> \qquad q:=|E|+|D|\ge |L-2m|. \tag{1}
> \]
> Equality holds exactly when \(\min(|E|,|D|)=0\). Equivalently, there is no cancelling
> empty--double pair: after the \(|L-2m|\) forced empty or double defects, every remaining rung is
> singly occupied.

> **[THEOREM + EXACT FINITE COMPUTATION].** For the deterministic raw-leaf residual minors stored by
> e237, and only for \((L,m)=(8,4),(9,4),(9,5)\), every residual coordinate is on the equality
> boundary in (1). The \((8,4)\) core consists of 16 \((6,2)\) no-empty/no-double orbits and one
> \((8,0)\) orbit. The \((9,4)\) core consists of 128 \((6,2)\) one-empty/no-double orbits and five
> \((8,0)\) one-empty orbits. The \((9,5)\) core is its particle-hole mirror, with one double rung
> in place of one empty rung. Exact support-position blocks and their Smith forms explain
> \(-3\,2^{32}\) and \(-3^5 2^{288}\).

Only the first statement is an all-\(L\) theorem. The residual classification is a theorem about
three literal, basis- and leaf-certificate-bound finite matrices. It proves no recurrence, no
\(L\ge10\) rank, no arbitrary-basis obstruction, no Pfaffian-number lower bound, and no
thermodynamic or exact three-dimensional Ising solution.

## 2. Rung coordinates and the elementary defect lemma

Number rungs by \(0,\ldots,L-1\). A two-bit rung has one of four states:

\[
E=00,\qquad T=10,\qquad B=01,\qquad D=11.
\]

The four corresponding rung sets partition all rungs, so

\[
L=|T|+|B|+|E|+|D|. \tag{2}
\]

A top-only or bottom-only rung contributes one particle and a double rung contributes two. In
particle sector \(2m\),

\[
2m=|T|+|B|+2|D|. \tag{3}
\]

Subtracting (3) from (2) proves the first identity in (1). For nonnegative integers \(E=|E|\) and
\(D=|D|\),

\[
(E+D)-|E-D|=2\min(E,D)\ge0. \tag{4}
\]

This proves the inequality and its equality criterion. Notice what the criterion does and does not
say. It says that no empty defect is cancelled by a double defect. It does not assert that every
coordinate on this boundary survives e237's leaf peeling.

The producer exhaustively checks (1)--(4) on every even-particle rung word through \(L=9\). That
finite check validates the decoder; the proof above, not the enumeration, gives the all-\(L\)
lemma.

## 3. The explicit Klein-four action

As in e230 and e237, coordinates are orbits under

\[
K=\langle\tau,\rho\rangle\cong C_2\times C_2.
\]

Here \(\tau\) swaps top and bottom on every rung and \(\rho\) reverses rung order:

\[
\tau:T\leftrightarrow B,\ E\mapsto E,\ D\mapsto D,
\qquad
\rho:r\mapsto L-1-r. \tag{5}
\]

The implementation applies all four maps \(1,\tau,\rho,\tau\rho\) to literal integer-encoded
configurations. For every claimed family it forms the complete configuration closure, partitions it
directly into canonical orbits, and independently evaluates Burnside's formula

\[
\#(X/K)=\frac{|\operatorname{Fix}(1)|+|\operatorname{Fix}(\tau)|
 +|\operatorname{Fix}(\rho)|+|\operatorname{Fix}(\tau\rho)|}{4}. \tag{6}
\]

Thus the counts below are outputs of explicit actions, not inserted Burnside constants.

## 4. Complete finite residual classification

Write \((a,b)\), with \(a\ge b\), for the unordered pair
\((|T|,|B|)\). The reconstructed residual-coordinate table is

| core | \((|E|,|D|)\) | \(q\) | \((6,2)\) orbits | \((8,0)\) orbits | total |
|---|---:|---:|---:|---:|---:|
| \((8,4)\) | \((0,0)\) | 0 | 16 | 1 | 17 |
| \((9,4)\) | \((1,0)\) | 1 | 128 | 5 | 133 |
| \((9,5)\) | \((0,1)\) | 1 | 128 | 5 | 133 |

Every one of the 283 representatives is decoded in the machine artifact with its full
\((T,B,E,D)\) rung sets. The producer independently reproduces e237's deterministic singleton
leaf sequence, remaining representatives, active source rows, selected minor rows, and determinant.
It also enumerates every configuration with the displayed defect and single-occupancy family. The
union of the reconstructed core orbits equals that exhaustive family set, with no missing
configuration and no extra coordinate.

This last sentence is deliberately family-specific. For example, a \((7,1)\) word can also satisfy
\(q=|L-2m|\) and have leg imbalance at least four, but those coordinates are removed by the earlier
leaf block and are not part of these residual cores.

### 4.1 Burnside calculation at \(L=8\)

For the \((6,2)\) family, the orbit closure contains
\(2\binom82=56\) configurations. The fixed-point counts in the order
\((1,\tau,\rho,\tau\rho)\) are

\[
(56,0,8,0),
\]

so (6) gives \((56+8)/4=16\) orbits. The two all-top/all-bottom configurations form the
\((8,0)\) closure. Their fixed-point counts are \((2,0,2,0)\), giving one orbit. Hence
\(17=16+1\).

### 4.2 Burnside calculation at \(L=9\)

For \(m=4\), choose one empty rung, then choose the minority single leg. The \((6,2)\) closure has

\[
9\cdot 2\binom82=504
\]

configurations. Its fixed-point vector is \((504,0,8,0)\), so it has 128 orbits. The \((8,0)\)
closure has \(9\cdot2=18\) configurations and fixed-point vector \((18,0,2,0)\), hence five
orbits. The \(m=5\) fixed-point data are identical after exchanging empty and double rungs.

A finer reflection-position count explains both summands. Let

\[
p=\min(r,8-r)\in\{0,1,2,3,4\}
\]

for the unique empty or double rung at position \(r\). Then

| \(p\) | \((6,2)\) configurations | nonidentity fixed points | \((6,2)\) orbits | \((8,0)\) configurations | nonidentity fixed points | \((8,0)\) orbits |
|---:|---:|---:|---:|---:|---:|---:|
| 0,1,2,3 (each) | 112 | 0 | 28 | 4 | 0 | 1 |
| 4 | 56 | \(\operatorname{Fix}(\rho)=8\) | 16 | 2 | \(\operatorname{Fix}(\rho)=2\) | 1 |

Therefore

\[
128=4\cdot28+16,
\qquad
5=4\cdot1+1. \tag{7}
\]

### 4.3 Particle-hole bijection

Let \(c(x)\) complement all \(2L\) occupation bits. On each rung it exchanges

\[
E\leftrightarrow D,
\qquad T\leftrightarrow B. \tag{8}
\]

Bit complement commutes with both maps in (5), so it descends to a bijection of canonical
\(K\)-orbits. At \(L=9\) it sends sector eight to sector ten, hence \(m=4\) to \(m=5\), and
preserves the reflection class \(p\). The certificate records the complete 133-column permutation.
After that permutation, the selected \((9,5)\) core matrix equals the selected \((9,4)\) matrix
entry by entry, with the same row order.

## 5. Reconstructing the literal e237 matrices

Let \(o_j\) be a residual configuration orbit and let

\[
w_i=\sum_o c_{i,o}e_o
\]

be a selected exact annihilator basis vector, where \(e_o\) is the orbit sum. As in e230/e237, the
literal weighted restriction matrix is

\[
M_{ij}=|o_j|c_{i,o_j}. \tag{9}
\]

The orbit sizes are recomputed from (5), not copied. The clean-room verifier independently rebuilds
all high-imbalance coordinates and all sector-homogeneous rows from the stored \(W_8,W_9\) bases,
replays the singleton queue, selects the first nonzero residual minor in e237's bounded order, and
checks (9). It does not import e237 or e249.

The exact Smith forms of the three literal matrices are

| core | Smith diagonal multiplicities | determinant |
|---|---|---:|
| \((8,4)\), \(17\times17\) | \(2^4,4^{11},8,24\) | \(-3\,2^{32}\) |
| \((9,4)\), \(133\times133\) | \(2^4,4^{103},8^{21},24^5\) | \(-3^5 2^{288}\) |
| \((9,5)\), \(133\times133\) | \(2^4,4^{103},8^{21},24^5\) | \(-3^5 2^{288}\) |

Here exponents on diagonal values denote multiplicity, not prime valuation. Every invariant factor is
nonzero, so every matrix has full rational rank. The product of the invariant factors gives the
absolute determinant; recorded row/column permutation signs recover each displayed negative literal
determinant.

## 6. Packing/polar blocks

Order columns first by the \((6,2)\) family and then by \((8,0)\). Order rows first by those having
zero support on all polar columns. In all three cases the number of first rows equals the number of
packing columns, and the number of remaining rows equals the number of polar columns. The matrix is
therefore

\[
P_r M P_c^T=
\begin{pmatrix}
A&0\\
C&D
\end{pmatrix}. \tag{10}
\]

No Schur division is needed. Exact block data are

| core | Smith diagonal of \(A\) | \(D\) |
|---|---|---|
| \((8,4)\) | \(2^3,4^{11},8,24\) | \([2]\) |
| \((9,4)\), \((9,5)\) | \(2^3,4^{99},8^{21},24^5\) | \(\operatorname{diag}(4,4,4,4,2)\) |

The lower-left block \(C\) is not discarded: it has exact ranks 1 at \(L=8\) and 5 in each \(L=9\)
core. It does not affect the determinant in (10), and the full Smith forms above are computed from
the literal matrices, not inferred from the diagonal blocks.

## 7. Defect-position direct sum and the determinant mechanism

For \(L=9\), every literal row as well as every column is supported on one reflection class
\(p\). Simultaneous row and column permutations therefore split the full core into four
\(29\times29\) off-centre blocks and one central \(17\times17\) block:

\[
M_9\sim M_{9,0}\oplus M_{9,1}\oplus M_{9,2}
          \oplus M_{9,3}\oplus M_{9,4}. \tag{11}
\]

For \(p=0,1,2,3\), each block has Smith diagonal

\[
4^{23},8^5,24,
\qquad |\det M_{9,p}|=3\,2^{64}. \tag{12}
\]

For \(p=4\), delete the central empty rung. For \(m=5\), first use particle-hole after deleting the
central double rung. The resulting 17 columns are exactly the \(L=8\) columns. A recorded row
permutation with one row sign change makes the central block equal to the literal \(L=8\) core.
Thus

\[
\operatorname{SNF}(M_{9,4})=(2^4,4^{11},8,24),
\qquad |\det M_{9,4}|=3\,2^{32}. \tag{13}
\]

Equations (11)--(13) give the complete prime-power mechanism:

\[
|\det M_9|=(3\,2^{64})^4(3\,2^{32})
          =3^5 2^{288}. \tag{14}
\]

Likewise (10) at \(L=8\) gives

\[
|\det M_8|=(3\,2^{31})\cdot2=3\,2^{32}. \tag{15}
\]

All signs in (14)--(15) are recovered by the explicit permutation signs stored in the artifact. No
modular rank, floating-point determinant, or fitted formula is promoted.

## 8. The first exact obstruction to the count-driven recurrence

The finite identity (11) contains one inherited \(L=8\) block, but it does not give a uniform
recurrence. A particularly tempting count-only proposal is

\[
A_9\stackrel{?}{\sim}A_8^{\oplus8}, \tag{16}
\]

because the packing dimensions satisfy \(128=8\cdot16\). If `~` means unimodular row/column
equivalence, determinant magnitude is an invariant. But

\[
|\det(A_8^{\oplus8})|=(3\,2^{31})^8=3^8 2^{248},
\qquad
|\det A_9|=3^5 2^{279}. \tag{17}
\]

They differ, so (16) is impossible. The first off-centre shell already pinpoints the failure. Every
entry of a \(29\times29\) shell is divisible by four, but after dividing by four its Smith diagonal
is

\[
1^{23},2^5,6. \tag{18}
\]

The invariant 6 contributes a factor of 3. Thus a proposed dyadic/unimodular new-shell recurrence
fails at the first sorted defect class \(p=0\), and each of the four off-centre shells contributes
one of the four new ternary factors in (14).

This is an exact obstruction only to the named eight-copy/dyadic-shell proposal. It is not a no-go
for every possible row combination, filtration, or block recurrence.

## 9. Why no all-\(L\) recurrence is promoted

There are only two residual sizes here, and the second uses the stored \(W_9\) basis. No exact
\(W_{10}\) source is available in this front. In particular:

1. the equality boundary in (1) is elementary, but equality between that boundary family and a
   future raw-leaf residual core is not;
2. the row support by a single defect-position class is checked only in the two \(L=9\) matrices;
3. central-rung deletion is an exact finite identity from \(L=9\) to \(L=8\), not an induction
   hypothesis;
4. the four new blocks in (11) have nontrivial ternary Smith torsion, so even the simplest proposed
   recurrence fails; and
5. fitting determinant exponents or a recurrence from \(L=8,9\) is explicitly forbidden and was not
   done.

Accordingly all residual-core and two-slice rank questions for \(L\ge10\) remain
**[UNRESOLVED]**.
