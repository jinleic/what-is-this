# Open \(3\times L\) two-sum DLA: finite exponential family and exact obstruction

## 1. Result and scope

Let \(G_L\) be the open \(3\times L\) grid, with sites numbered row-major, and define

\[
 A_L=\sum_{v\in G_L}X_v,
 \qquad
 B_L=\sum_{(u,v)\in E(G_L)}Z_uZ_v,
 \qquad
 \mathfrak g_L=\operatorname{Lie}_{\mathbb Q}\langle A_L,B_L\rangle.
\]

This front sought an all-length bound \(\dim\mathfrak g_L\ge c^L\), with \(c>1\), or a rigorous
finite advance beyond the existing two-leg work.

**[UNRESOLVED — all lengths].** No all-\(L\) exponential theorem is proved. The positive result is a
new finite family certificate:

\[
 \dim_{\mathbb Q}\mathfrak g_L\ge 2^L
 \qquad (2\le L\le6). \tag{1}
\]

For each length in (1), the certificate consists of \(2^L\) explicit integer left-nested Lie words
and a nonzero modular minor at each of two primes. For \(2\le L\le5\), the selected square
coefficient matrix was also evaluated by fraction-free integer elimination and its nonzero
integer determinant is stored. Thus the family is an actual collection of Lie elements, not a
support count or a fit.

The requested induction was attempted and stopped honestly. The chosen minors do not obey a proved
factor recurrence, no representation-evaluation intertwiner was established, and every order in a
natural 96-element product-lex class fails to assign distinct raw leaders by \(L=4\). Finite ranks
are not promoted to an all-size theorem.

A second finite advance extends the open \(3\times4\) depth profile by one complete depth:

\[
 D_{14}=741,\qquad D_{15}=1242 \pmod {2147483647}, \tag{2}
\]

so \(\dim_{\mathbb Q}\mathfrak g_4\ge1242\). This is not the full \(3\times4\) dimension.

## 2. Exact Pauli convention

Encode

\[
 Q_{(a\mid b)}=X^aZ^b
\]

by the integer `a | (b << n)`, where \(n=3L\). The experiment divides every commutator by two:

\[
 \frac12[Q_g,Q_v]=
 \begin{cases}
 0,&\langle g,v\rangle=0,\\
 \epsilon(g,v)Q_{g\oplus v},&\langle g,v\rangle=1,
 \end{cases}
 \qquad \epsilon(g,v)\in\{+1,-1\}. \tag{3}
\]

All word coefficients and all reported integer determinants are evaluated exactly over
\(\mathbb Z\). Dividing by two changes only nonzero scalar factors and hence changes neither the
Lie span nor linear independence.

The modular eliminations use

\[
 p_1=2147483647,
 \qquad
 p_2=2147483629.
\]

**[LEMMA — modular independence].** If an integer coefficient matrix has rank \(r\) modulo a prime
\(p\), then the represented Lie elements have rank at least \(r\) over \(\mathbb Q\).

**Proof.** Modular rank \(r\) supplies an \(r\times r\) integer minor whose determinant is nonzero
modulo \(p\). That integer determinant is nonzero, so the same minor is nonzero over
\(\mathbb Q\). \(\square\)

Agreement at two primes is only an arithmetic cross-check. It is not a rational upper-bound
certificate.

## 3. Modular closure and depth profiles

For

\[
 F_1=\operatorname{span}\{A_L,B_L\},
 \qquad
 F_{k+1}=F_k+[A_L,F_k]+[B_L,F_k],
 \qquad D_k=\dim F_k,
\]

`experiments/e74_dla_3xl.py` records the following complete-depth modular profiles:

| shape | complete modular profile \(D_1,D_2,\ldots\) | scope |
|---|---|---|
| \(3\times2\) | 2, 3, 5, 7, 11, 16, 26, 38, 57, 85, 129, 182, 245, 262, 263, 263 | saturated mod \(p_1,p_2\), paired with an exact-Q certificate; \(\dim_{\mathbb Q}=263\) |
| \(3\times3\) | 2, 3, 5, 7, 11, 16, 26, 41, 65, 98, 153, 236, 375, 592, 947, 1504, 2391, 3735, 5805, 7473, 8034, 8034 | saturated mod \(p_1,p_2\) in the prior certificate; \(\dim_{\mathbb Q}\ge8034\) |
| \(3\times4\) | 2, 3, 5, 7, 11, 16, 26, 41, 67, 105, 169, 271, 448, 741, 1242 | complete only through depth 15 mod \(p_1\); \(\dim_{\mathbb Q}\ge1242\) |

The \(3\times2\) graph is isomorphic to the previously studied \(2\times3\) graph by transposition.
The new script nevertheless recomputes its orbit closure at both primes. It obtains 304 reachable
Pauli-symmetry orbits, modular row count 263, grading-block column counts
\((92,60,76,76)\), and row counts \((78,51,67,67)\).

**[THEOREM — exact finite dimension].**
\(\dim_{\mathbb Q}\mathfrak g_{3\times2}=263\).

**Proof.** Grid transposition \((r,c)\mapsto(c,r)\) is a graph isomorphism from the open
\(3\times2\) grid to the open \(2\times3\) grid. It maps \(A=\sum_vX_v\) and
\(B=\sum_{uv\in E}Z_uZ_v\) term-by-term, so it induces an isomorphism of the generated rational Lie
algebras. The existing `results/algebra_structure/char0_levi.json` certificate independently gives
both a 263-word exact Pauli closure and faithful exact-\(\mathbb Q\) matrix rank 263 for the
\(2\times3\) algebra. Hence the isomorphic \(3\times2\) algebra has exact dimension 263.
\(\square\)

The full \(3\times3\) modular row count 8034 was already present. This front verifies its scope:
it is a completed finite-field closure, and therefore a rigorous lower bound over \(\mathbb Q\),
not an exact rational dimension equality. No characteristic-zero upper certificate is supplied.

For \(3\times4\), the orbit engine uses row and column reflections, the global-spin-flip restriction
that every reached Pauli word has even \(Z\)-parity, and the exact four-block grading by
`(X-count parity, transpose parity)`. The depth-15 run used 773,276,448 bytes of sparse-basis
accounting, observed peak RSS 989,380,608 bytes, and 456.747006 seconds. It stopped at the requested
depth boundary, not at saturation. Consequently 1242 in (2) is a lower bound only.

## 4. Explicit binary Lie-word family

A word string is outermost-first: `ABBA` means
\([A,[B,[B,A]]]\). At \(L=2\), define

\[
\begin{array}{c|c}
00&A,\\
01&BA,\\
10&ABA,\\
11&BABA.
\end{array}
\]

For every added bit use

\[
 w_{L+1}(s0)=AB\,w_L(s),
 \qquad
 w_{L+1}(s1)=BB\,w_L(s). \tag{4}
\]

Identify a label \(s=s_1\cdots s_L\) with the subset
\(S(s)=\{i:s_i=1\}\subseteq\{1,\ldots,L\}\). Thus (4) is an explicit family indexed by every
subset of columns; this indexing does not assert that the resulting Pauli support is column-local.

Each step adds two brackets, so every word in (4) has depth at most \(2L\) and belongs to
\(F_{2L}\subseteq\mathfrak g_L\).

**[COMPUTATION — finite family].** Exact integer expansion followed by modular elimination gives

| \(L\) | words | maximum depth | rank mod \(p_1\) | rank mod \(p_2\) |
|---:|---:|---:|---:|---:|
| 2 | 4 | 4 | 4 | 4 |
| 3 | 8 | 6 | 8 | 8 |
| 4 | 16 | 8 | 16 | 16 |
| 5 | 32 | 10 | 32 | 32 |
| 6 | 64 | 12 | 64 | 64 |

This proves (1) by the modular-independence lemma. The JSON certificate stores every binary label,
its corresponding subset of columns, concrete word, exact support size, and SHA-256 digest of the
full `(Pauli code, integer coefficient)` vector.
The standalone test independently expands and eliminates the three consecutive sizes
\(L=4,5,6\).

### Exact integer minors

Sparse elimination modulo \(p_1\) selects an explicit Pauli column for each row. Re-evaluating those
same square matrices over \(\mathbb Z\) by Bareiss elimination gives

| \(L\) | minor size | exact determinant |
|---:|---:|---:|
| 2 | 4 | \(4\) |
| 3 | 8 | \(-41472\) |
| 4 | 16 | \(621154997895168000\) |
| 5 | 32 | \(821526282062651783776259662615193653063452217835520000\) |

Every determinant is nonzero modulo both recorded primes as well as nonzero in \(\mathbb Z\). The
artifact stores the actual selected Pauli columns and decoded three-letter rung strings, so these
are reproducible explicit certificates rather than rank-only claims.

The determinant sequence is **not** claimed to obey a recurrence. The columns were selected by
global finite-size elimination, and neither their boundary behavior nor a block factorization was
proved.

## 5. Exact obstruction to the natural leading-term induction

A natural route to promote (4) would assign each raw word a distinct leading Pauli monomial under a
position-independent product order. The experiment exhausts

\[
2\text{ column directions}
\times2\text{ row directions}
\times24\text{ strict local orders on }\{I,X,Y,Z\}
=96
\]

orders.

The exact histograms of distinct raw leaders are

| \(L\) | histogram `leaders: orders` | maximum | required |
|---:|---|---:|---:|
| 4 | 5:12, 6:12, 9:24, 10:4, 11:16, 12:16, 13:12 | 13 | 16 |
| 5 | 5:12, 6:12, 10:12, 12:4, 13:16, 15:12, 16:4, 18:24 | 18 | 32 |

Thus no order in this class proves independence of the family even at \(L=4\). This is stronger
than exhibiting one collision: the entire tested order class has insufficient leader capacity.
It does not prove that another monomial order, an echelon multiword family, or a representation
certificate cannot work.

A second possible promotion would factor the selected evaluation matrix under the lift
\(L\mapsto L+1\). The finite minors above do not do so: their pivot columns change globally, and no
local block recurrence with nonzero determinant was found. Since a recurrence statement is the
induction hypothesis—not a consequence of finitely many ranks—none is asserted.

## 6. What is proved and what remains open

The exact conclusions are:

1. \(3\times2\) recomputes to 263 modular closure rows at both primes. Graph transposition paired
   with the existing exact-\(\mathbb Q\) \(2\times3\) certificate proves
   \(\dim_{\mathbb Q}\mathfrak g_2=263\).
2. The existing 8034 count for \(3\times3\) is a saturated finite-field count. The rational
   conclusion justified by that certificate is \(\dim\mathfrak g_3\ge8034\), not equality.
3. The \(3\times4\) complete depth profile is extended to \(D_{15}=1242\), proving
   \(\dim\mathfrak g_4\ge1242\).
4. The explicit family (4) proves \(\dim\mathfrak g_L\ge2^L\) for exactly \(2\le L\le6\).
5. The 96 natural product-lex orders cannot turn that family into a raw-leading-term induction from
   \(L=4\) onward.

**[UNRESOLVED].** No \(c>1\) is certified for every \(L\). A successful next step must supply an
actual all-size mechanism: for example, a representation evaluation whose matrix has a proved local
block recurrence, or a different multiword family with a stable boundary signature. More finite
sizes alone cannot close that gap.

## 7. Reproducibility

From the repository root run

```sh
.venv/bin/python experiments/e74_dla_3xl.py
.venv/bin/python tests/test_dla_3xl.py
```

Both programs print final `PASS`. The experiment artifact is
`results/algebra_growth/dla_3xl.json` and has the required `provenance`, `data`, and `checks`
envelope. The standalone test does not import the experiment: it independently rebuilds the tiny
\(3\times2\) orbit closure at \(p_2\), the family ranks at \(L=4,5,6\), the exact integer minors at
\(L=4,5\), both 96-order histograms, and all scope guards.
