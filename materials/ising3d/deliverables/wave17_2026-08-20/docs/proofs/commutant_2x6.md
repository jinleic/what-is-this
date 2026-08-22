# Exact joint commutant of the open `2 x 6` Ising ladder

## Statement and scope

Let `G` be the open two-row, six-column ladder, with its twelve vertices in
row-major order, and define

\[
 A=\sum_{v\in G}X_v,\qquad B=\sum_{\{u,v\}\in E(G)}Z_uZ_v.
\]

**[THEOREM — finite `2 x 6` result].** Over the complex numbers,

\[
 \dim_{\mathbb C}\operatorname{Comm}(A,B)=71,
 \qquad
 \operatorname{Comm}(A,B)\cong
 M_7(\mathbb C)\oplus M_2(\mathbb C)^2\oplus\mathbb C^{14}.
\]

The integer-coordinate commutator systems below show that the rational and
complex dimensions agree, so the dimension is also exactly `71` over
\(\mathbb Q\).  The split-type assertion in this note is deliberately stated
over \(\mathbb C\).

This is a theorem about this one finite open ladder.  It neither gives an
all-`L` formula nor supplies a local or extensive conserved charge.

The producer is `experiments/e89_commutant_2x6.py`, the certificate is
`results/integrability/commutant_2x6.json`, and the independent rebuild is
`tests/test_commutant_2x6.py`.

## Exact character-Hom decomposition

Let `R`, `C`, and `P` respectively denote row reflection, column reflection,
and global spin flip.  They commute with both `A` and `B` and generate
\(H\cong C_2^3\).  For
\(\chi=(\chi_R,\chi_C,\chi_P)\in\mathbb F_2^3\), write

\[
 V_\chi=\{v:Rv=(-1)^{\chi_R}v,\ Cv=(-1)^{\chi_C}v,
 Pv=(-1)^{\chi_P}v\}.
\]

**[LEMMA — exact orbit basis].** For an `H`-orbit representative `r`, the
nonzero vectors

\[
 s_{\chi,r}=\sum_{g\in H}(-1)^{\chi\cdot g}|gr\rangle
\]

supply a rational basis of `V_chi`.  In this basis `B` is diagonal with
integer entries, `A` has integer entries, and the positive diagonal orbit
Gram matrix makes both restrictions self-adjoint.

*Justification.* Every group element preserves both the ladder edge set and
the transverse-field sum.  Thus character projection produces invariant
summands.  All states in an orbit have the same `B` energy, and flipping one
spin gives integral orbit-basis coefficients for `A`.  Direct exact
verification of the orbit-metric Hermiticity is included in the producer and
certificate.  ∎

The sector dimensions, in the order

```text
000, 001, 010, 011, 100, 101, 110, 111
```

are exactly

\[
(560,512,496,512,496,512,496,512).
\]

Since `A` and `B` are block diagonal for this splitting,

\[
 \operatorname{End}(V)=\bigoplus_{\chi,\psi}
 \operatorname{Hom}(V_\psi,V_\chi),
\]

and the joint commutant is the direct sum of the joint Hom kernels.  This is
an equality of vector spaces, not a symmetry-based lower bound.  Orbit-metric
adjunction gives an exact bijection between opposite Hom blocks, so only the
36 upper-triangle blocks need be eliminated.

For a target/source Hom matrix `X`, first impose `XB=BX`.  Because `B` is
diagonal, this retains only entries joining equal `B` energies.  The remaining
integer equations are

\[
 XA_{\rm source}-A_{\rm target}X=0. \tag{1}
\]

The largest resulting system has `43,664` variables, compared with
\(4096^2=16,777,216\) variables in the raw full-endomorphism calculation.

## Characteristic-zero certificate

**[LEMMA — same-basis modular/exact closure].** Let `M` be one of the integer
matrices for (1), with `n` columns.  Suppose its rank over
\(\mathbb F_p\) is `r`, and there are `n-r` rational vectors in its exact
kernel that are independent over \(\mathbb Q\).  Then

\[
 \operatorname{rank}_{\mathbb Q}M=r,
 \qquad \dim_{\mathbb Q}\ker M=n-r.
\]

*Proof.* Reduction modulo `p` gives
\(\operatorname{rank}_{\mathbb F_p}M\leq\operatorname{rank}_{\mathbb Q}M\).
The displayed exact kernel basis gives
\(\operatorname{rank}_{\mathbb Q}M\leq n-(n-r)=r\).  ∎

**[COMPUTATION — rational kernel construction].** For every upper-triangle
block, the producer performs exact sparse RREF at

\[
 p_1=1{,}000{,}000{,}007,\qquad p_2=2{,}147{,}483{,}647.
\]

The two ranks agree blockwise.  For every nonzero modular nullspace, its two
bases are normalized at the same independent coordinate anchors, combined by
CRT, and rationally reconstructed.  Every reconstructed vector is then
substituted into every integer equation (1), and its anchor identity matrix
proves rational independence.  The JSON stores the two ranks, variable count,
anchor minor residues, coefficient-bit bound, and a SHA-256 digest of each
basis.

Agreement of the two modular ranks is only a cross-check: it is **not** used
as an equality claim by itself.  The exact reconstructed vectors supply the
other inequality in the lemma.

The exact Hom-dimension table, with rows indexed by target and columns by
source, is

\[
H=\begin{pmatrix}
3&0&3&0&2&0&2&0\\
0&3&0&0&0&0&0&0\\
3&0&6&0&4&0&4&0\\
0&0&0&3&0&0&0&0\\
2&0&4&0&5&0&5&0\\
0&0&0&0&0&3&0&0\\
2&0&4&0&5&0&5&0\\
0&0&0&0&0&0&0&3
\end{pmatrix}.
\]

Its entries sum to

\[
\boxed{\sum_{\chi,\psi}H_{\chi\psi}=71}.
\]

Because the systems have integer coordinates in a full matrix basis, tensoring
the exact rational kernels with \(\mathbb C\) gives the complex kernels with
the same dimensions.

## Common-zero block and complex split type

**[COMPUTATION — exact common-zero space].** Exact rational elimination of
\(W=\ker A\cap\ker B\) inside the eight character sectors gives

\[
w=(\dim(W\cap V_\chi))_\chi=(1,0,2,0,2,0,2,0),
\qquad \dim W=7.
\]

`W` is reducing: `A` and `B` are self-adjoint and vanish on it.  Thus
\(\operatorname{End}(W)\cong M_7(\mathbb C)\) is a direct commutant summand.
There are no intertwiners from `W` to `W^\perp`: their images would be common
zero vectors in `W^\perp`; adjunction eliminates the reverse direction as
well.  Consequently the exact Hom table on `W^\perp` is
\(H-ww^{\mathsf T}\), namely

\[
H_{W^\perp}=\begin{pmatrix}
2&0&1&0&0&0&0&0\\
0&3&0&0&0&0&0&0\\
1&0&2&0&0&0&0&0\\
0&0&0&3&0&0&0&0\\
0&0&0&0&1&0&1&0\\
0&0&0&0&0&3&0&0\\
0&0&0&0&1&0&1&0\\
0&0&0&0&0&0&0&3
\end{pmatrix}. \tag{2}
\]

**[LEMMA — reading the residual Hom table].** The complex commutant of a
finite Hermitian pair is a finite-dimensional complex `*`-algebra and hence
semisimple.  A diagonal endomorphism algebra of dimension two or three can
therefore only be \(\mathbb C^2\) or \(\mathbb C^3\): a non-scalar matrix
block has dimension at least four.

In (2), the `[000,010]` subtable `[[2,1],[1,2]]` therefore comprises one
shared multiplicity-two irreducible and two unmatched scalar components.  It
contributes \(M_2(\mathbb C)\oplus\mathbb C^2\).  The `[100,110]` subtable
`[[1,1],[1,1]]` comprises one shared multiplicity-two irreducible and
contributes a second \(M_2(\mathbb C)\).  The four isolated diagonal blocks
of dimension three contribute \(\mathbb C^{12}\).  These exact nonzero Hom
spaces are the intertwiner witnesses; zero entries exclude any further
identifications.

Combining this with the `M_7` common-zero block gives

\[
M_7(\mathbb C)\oplus M_2(\mathbb C)^2\oplus\mathbb C^{14},
\]

whose dimension is `49 + 4 + 4 + 14 = 71` and whose centre has dimension
`17`.  This proves the theorem.

## Controls and previous bounds

**[COMPUTATION — same-pipeline `2 x 5` control].** Re-running the same
character split, two-prime reconstruction, and exact substitution on `2 x 5`
reproduces all 36 stored Hom dimensions, their sum `20`, and

\[
M_2(\mathbb C)^2\oplus\mathbb C^{12}.
\]

The final producer records `11.760111` seconds for this control and
`274.407279` seconds for the `2 x 6` closure; its complete run took
`286.207903` seconds on the stated workstation.  A fresh producer-plus-
standalone-verifier command completed in `585.80` seconds and printed

```text
clean-room 2x5=20; 2x6=71; C-type=M_7(C)+M_2(C)^2+C^14
PASS
```

**[COMPUTATION — old gap cross-check].** The prior family certificate recorded
an exact lower bound `57` from the seven-dimensional common-zero space and an
upper bound `71` from a single modular scout.  The present exact Hom closure
reproduces the same common-zero data and modular table, then proves

\[
57\leq 71=\dim\operatorname{Comm}(A,B)\leq71.
\]

**[UNRESOLVED].** This finite theorem does not establish a recurrence in ladder
length, the locality of any non-geometric matrix unit, integrability of the
finite ladder, or an exact solution of the three-dimensional Ising model.
