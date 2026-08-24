# Boundary-row restriction for the two-generator rectangle algebra

## 1. Scope and convention

**[EXTERNAL].** For an open rectangle, write

\[
 A_{r,L}=\sum_v X_v,\qquad B_{r,L}=\sum_{uv\in E(r\times L)}Z_uZ_v,
 \qquad \mathfrak g_{r,L}=\operatorname{Lie}_{\mathbb Q}\langle iA_{r,L},iB_{r,L}\rangle .
\]

**[LEMMA].** It is convenient to use the Hermitian integral bracket

\[
 \{S,T\}:=\frac{[S,T]}{2i}.
\]

Every word in this bracket is a nonzero scalar multiple of the corresponding anti-Hermitian Lie word, so membership and linear independence are unchanged.

**[EXTERNAL].** `proofs/ladder_w8.md` proves the desired dimension inequality for every open `2xL` ladder with `L>=3`; `proofs/dla_3xl.md` supplies only finite `3xL` anchors and explicitly leaves an all-length induction unresolved; `proofs/twogen_allsize.md` shows that `ad_A` eigenvalue separation cannot isolate path bonds from chords.

**[THEOREM — headline route closure].** Deleting one boundary row does not give a tensor restriction or a Lie-homomorphic `|+>` compression from the `(r+1)xL` two-sum algebra to the `rxL` two-sum algebra.  The tensor restriction is impossible for every fixed rank-one row state, and the `|+>` compression first fails at the depth-three word `\{B,\{A,B\}\}` by an explicit noncentral boundary field for every `r>=2` and `L>=1`.

## 2. No rank-one frozen-row restriction

Let `C` be the `rxL` core and let the removed boundary row have sites `e_1,...,e_L`.  Denote the adjacent core boundary sites by `s_1,...,s_L`.  The generators decompose as

\[
\begin{aligned}
 A_{r+1,L}&=A_C\otimes I+I\otimes A_E,\\
 B_{r+1,L}&=B_C\otimes I+I\otimes H_E+\sum_{j=1}^L Z_{s_j}\otimes Z_{e_j},
\end{aligned}
\]

where `A_E=sum_j X_{e_j}` and `H_E=sum_{j=1}^{L-1}Z_{e_j}Z_{e_{j+1}}`.

**[LEMMA — operator-Schmidt separation].** If the full tensor copy

\[
 \mathcal H_C\otimes\operatorname{span}\{\Omega\}
\]

is invariant under `B_(r+1,L)`, then `Z_(e_j) Omega` belongs to `span{Omega}` for every `j`.

**[THEOREM — proof of the lemma].** Let `q` be the projection of the extra-row space onto the quotient by `span{Omega}`.  Invariance under `B` gives, for every core vector `v`,

\[
 v\otimes qH_E\Omega+
 \sum_{j=1}^L Z_{s_j}v\otimes qZ_{e_j}\Omega=0.
\]

The core Pauli operators `I,Z_(s_1),...,Z_(s_L)` are linearly independent.  Taking their Hilbert--Schmidt coefficient functionals forces `qH_E Omega=0` and `qZ_(e_j)Omega=0` separately.  This proves the lemma.

**[THEOREM — no frozen row].** No nonzero `Omega` makes the same tensor copy invariant under both `A_(r+1,L)` and `B_(r+1,L)`.

**[THEOREM — proof].** The preceding lemma makes `Omega` a simultaneous eigenvector of every commuting involution `Z_(e_j)`, hence a computational-basis vector.  Invariance under `A` would additionally require `A_E Omega` to lie in `span{Omega}`.  Instead,

\[
 A_E\Omega=\sum_{j=1}^L X_{e_j}\Omega
\]

is a sum of `L` distinct computational-basis vectors, each orthogonal to `Omega`; it is nonzero.  This contradiction holds for every positive `L`.

**[UNRESOLVED].** This theorem excludes a tensor copy of the entire core with a one-dimensional frozen ancillary row.  It does not exclude a higher-rank ancillary code, a smaller specially chosen invariant subspace, or a non-tensor/nonlocal representation embedding.

## 3. Exact failure of `|+>` row compression

Put `Omega_+=|+>^(tensor L)` on the removed row and define the partial matrix element

\[
 \Phi(T)=(I\otimes\langle\Omega_+|)T(I\otimes|\Omega_+\rangle).
\]

Let `S_partial=sum_(j=1)^L X_(s_j)` be the field on the core boundary adjacent to the removed row.

**[LEMMA — agreement through the first bracket].** Exact Pauli multiplication gives

\[
\begin{aligned}
 \Phi(A_{r+1,L})&=A_{r,L}+L I,\\
 \Phi(B_{r+1,L})&=B_{r,L},\\
 \Phi(\{A_{r+1,L},B_{r+1,L}\})
   &=\{A_{r,L},B_{r,L}\}.
\end{aligned}
\]

**[THEOREM — depth-three anomaly].** For every `r>=1` and `L>=1`,

\[
\boxed{
 \Phi\!\left(\{B_{r+1,L},\{A_{r+1,L},B_{r+1,L}\}\}\right)
 =\{B_{r,L},\{A_{r,L},B_{r,L}\}\}
  +S_{\partial}+(3L-2)I .
}
\]

**[LEMMA — one-edge identity].** For an edge `e={u,v}` with `b_e=Z_uZ_v`,

\[
 \{b_e,\{A,b_e\}\}=X_u+X_v.
\]

**[THEOREM — proof of the anomaly].** Expand the depth-three word over ordered pairs of edges.  The contribution from two core edges is exactly the core depth-three word.  A pair of distinct edges contributes only if the edges meet; its Pauli word has `X` on the common endpoint and `Z` on both unmatched endpoints.  Every such mixed or removed-row cross term retains a `Y` or `Z` on the removed row and therefore has zero `|+>` matrix element.  Only equal outside edges survive.  Each of the `L` seam edges contributes `X_(s_j)+I`.  Each of the `L-1` horizontal edges in the removed row contributes `2I`.  Their sum is

\[
 S_{\partial}+L I+2(L-1)I=S_{\partial}+(3L-2)I,
\]

which proves the formula.

**[THEOREM — noncentral obstruction].** For `r>=2`, the compressed depth-three word is not in `mathfrak g_(r,L)+Q I`.

**[THEOREM — proof].** Reflection of the core across its rows fixes `A_(r,L)` and `B_(r,L)` pointwise, hence fixes every element of `mathfrak g_(r,L)` and also fixes `I`.  It sends `S_partial` to the field on the opposite boundary row.  These are distinct Pauli sums when `r>=2`.  The boxed expression is therefore not reflection-fixed and cannot belong to `mathfrak g_(r,L)+Q I`.

**[THEOREM — exact counterexample to induction].** The map `Phi` agrees with the desired generator images and even with their first bracket, but it is not a Lie homomorphism on the generated word algebra.  Thus applying the same ladder Lie words upstairs and compressing them cannot be used as a boundary-row induction without additional correction operators; the first missing correction is exactly `S_partial+(3L-2)I`.

**[THEOREM].** This obstruction is independent of the `ad_A` grade-separation wall: it comes from excursions into and back out of the removed row, equivalently from the failure of partial matrix elements to preserve products and brackets.

## 4. Bounded exact anchors

**[COMPUTATION].** `experiments/e197_rectangle_words.py` independently expands integral Hermitian Pauli sums for core shapes `2x1`, `2x2`, `2x3`, `2x4`, `3x2`, and `3x3`.  In every case it verifies the two generator compressions, first-bracket agreement, the boxed depth-three defect, and the nonzero row-reflection residual.

**[COMPUTATION].** `experiments/e198_rectangle_induction.py` uses exact rational elimination on the evaluation columns of `I,Z_1,...,Z_L` and exact integer computational-basis actions for `1<=L<=6`.  It obtains rank `L+1`, transverse-field support size and squared norm `L`, zero overlap with the starting joint-`Z` eigenvector, and scalar accounting `L+2(L-1)=3L-2`.

**[COMPUTATION].** `experiments/e199_rectangle_certificate.py` writes `results/algebra_growth/rectangle_twogen.json` with top-level `meta`, `data`, and computed `checks`.  `tests/test_rectangle_twogen.py` is a clean-room verifier using a bit-mask symplectic Pauli representation and imports none of the producers.

**[COMPUTATION].** Reproduce from the repository root with

```sh
nice -n 10 .venv/bin/python experiments/e197_rectangle_words.py
nice -n 10 .venv/bin/python experiments/e198_rectangle_induction.py
nice -n 10 .venv/bin/python experiments/e199_rectangle_certificate.py
nice -n 10 .venv/bin/python tests/test_rectangle_twogen.py
```

## 5. Result boundary

**[THEOREM].** The natural ladder-to-grid route by freezing one added row is closed for every width, and the natural `|+>` compression route is closed by one explicit all-size depth-three counterexample rather than by finite trend extrapolation.

**[UNRESOLVED].** No all-`3xL` or all-rectangle lower bound \(\dim\mathfrak g_{r,L}>rL(2rL-1)\) is proved here.  In particular, the exact route closure neither proves nor disproves that inequality for `3xL` with `L>=5` or for rectangles of height at least four.

**[UNRESOLVED].** Direct Lie-word families with stable boundary signatures, higher-rank invariant codes, and nonlocal representation restrictions remain possible.  Any future row induction must explicitly cancel or generate the boundary field in the boxed anomaly before it can inherit a ladder certificate.

**[THEOREM].** No finite anchor in this note is promoted to an all-size dimension statement, and no claim is made about a thermodynamic free energy or critical coupling.
