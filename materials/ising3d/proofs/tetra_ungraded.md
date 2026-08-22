# Arbitrary-$8\times8$ auxiliary-$R$ tetrahedron/RLLL obstruction

[THEOREM] This note removes the auxiliary $\mathbb Z_2$ grading assumption from the identical-$L$, two-dimensional-leg tetrahedron ansatz. It does **not** solve the three-dimensional Ising model.

## Exact setup

[DEFINITION] For binary triples $o,i\in\{0,1\}^3$, let $|o|$ denote Hamming weight and set

$$
L(q)_{o,i}=
\begin{cases}
q^{(|o|+|i|)/2},& |o|+|i|\equiv0\pmod 2,\\
0,& |o|+|i|\equiv1\pmod 2.
\end{cases}
$$

This is the verified six-leg Ising tensor in the convention $L[\mathrm{output},\mathrm{input}]=T[\mathrm{input},\mathrm{output}]/2$. The omitted factor two occurs three times on both sides and therefore does not change the equation.

[DEFINITION] On six binary spaces, in zero-based bit positions, use

$$
(123),(145),(246),(356)=(0,1,2),(0,3,4),(1,3,5),(2,4,5).
$$

Write

$$
F(q)=L_{145}(q)L_{246}(q)L_{356}(q),\qquad
G(q)=L_{356}(q)L_{246}(q)L_{145}(q).
$$

The equation under study is

$$
R_{123}F(q)=G(q)R_{123}. \tag{1}
$$

There is no condition on $R$: its 64 matrix entries are independent. Expanding all $64^2$ components of (1) gives the exact homogeneous system

$$
B(q)\,\operatorname{vec}(R)=0,\qquad B(q)\in\operatorname{Mat}_{4096\times64}(\mathbb Z[q]). \tag{2}
$$

The implementation is `experiments/e70_tetra_ungraded.py`; its artifact is `results/integrability/tetra_ungraded.json`.

## Why the block ordering is not a grading assumption

[LEMMA] The 64 columns of $B$ split into 32 parity-preserving and 32 parity-reversing columns, but **both** sets are retained.

Each local $L$ preserves the parity of the three spaces on which it acts. Hence $F$ and $G$ preserve total six-bit parity. A matrix unit of $R$ with auxiliary parity change $\epsilon$ contributes only to component rows with total parity change $\epsilon$. Consequently, after ordering the 32 parity-preserving coordinates first and the 32 parity-reversing coordinates second, (2) is a direct sum of two 4096-by-32 column systems supported on disjoint sets of 2048 rows.

[CERTIFICATE] At the exact specialization $q=1/4$, all 4096 component rows are nonzero in the full system. The test compares all $4096\times64=262144$ derived coefficients with direct 64-by-64 embedded-operator products. Both off-sector blocks contain exactly zero nonzero entries. Thus the direct-sum organization is only a proof device; no entry of $R$ has been discarded.

## Complementary maximal minors

[CERTIFICATE] Number a component row by

$$
\rho=64s_{\rm out}+s_{\rm in},\qquad 0\le s_{\rm out},s_{\rm in}<64.
$$

For the parity-preserving columns, select the rows

```text
0, 3, 5, 6, 9, 10, 12, 15, 17, 18, 33, 34, 36, 39, 40, 43,
45, 46, 65, 66, 68, 71, 192, 195, 257, 263, 264, 267, 269, 270, 272, 275.
```

Fraction-free exact elimination gives the 32-by-32 determinant

$$
d_{\rm even}(q)=-2q^{56}(q-1)^{45}(q+1)^{39}. \tag{3}
$$

For the parity-reversing columns, the primary row set is

```text
1, 2, 4, 7, 8, 11, 13, 14, 16, 19, 32, 35, 37, 38, 41, 42,
44, 47, 64, 67, 69, 70, 193, 194, 256, 259, 261, 262, 265, 266, 268, 271.
```

Its determinant is

$$
d_{\rm odd,0}(q)=-q^{53}(q-1)^{47}(q+1)^{38}(q^3+2q-1). \tag{4}
$$

Replace row 262 by row 274 to obtain the alternate set

```text
1, 2, 4, 7, 8, 11, 13, 14, 16, 19, 32, 35, 37, 38, 41, 42,
44, 47, 64, 67, 69, 70, 193, 194, 256, 259, 261, 265, 266, 268, 271, 274.
```

Its determinant is

$$
d_{\rm odd,1}(q)=-q^{54}(q-1)^{47}(q+1)^{38}(q^2+1). \tag{5}
$$

Because the complementary row/column blocks vanish exactly, concatenating the even rows with either odd row set gives literal 64-by-64 minors of the full matrix $B(q)$. Their determinants are

$$
\begin{aligned}
D_0(q)&=2q^{109}(q-1)^{92}(q+1)^{77}(q^3+2q-1),\\
D_1(q)&=2q^{110}(q-1)^{92}(q+1)^{77}(q^2+1).
\end{aligned} \tag{6}
$$

Let

$$
C(q)=q^{109}(q-1)^{92}(q+1)^{77},\quad
f(q)=q^3+2q-1,\quad g(q)=q(q^2+1).
$$

The two extra factors are coprime over $\mathbb Q$. The stored exact Bezout certificate is

$$
-\frac{q^2+q+2}{2}f(q)+\frac{q^2+q+3}{2}g(q)=1. \tag{7}
$$

Equivalently, the monic polynomial gcd of the two full-minor determinants over $\mathbb Q[q]$ is exactly $C(q)$.

## All-parameter conclusion

[THEOREM] Let $k$ be a characteristic-zero field and $q\in k$ with $q\notin\{0,1,-1\}$. Then the only arbitrary 8-by-8 matrix $R$ over $k$ satisfying (1) is $R=0$.

Indeed, $C(q)\ne0$. If both full minors in (6) vanished, then $f(q)=g(q)=0$, contradicting (7). Thus at least one 64-by-64 minor is nonsingular. Since $B(q)$ has exactly 64 columns,

$$
\operatorname{rank}B(q)=64,\qquad \ker B(q)=\{0\}.
$$

This is a symbolic all-$q$ statement, not an extrapolation from finitely many rational samples.

[COROLLARY] For every finite ferromagnetic Ising weight $0<q<1$, no nonzero, and hence no invertible, arbitrary 8-by-8 auxiliary $R$ satisfies the identical-$L$ relation (1).

## Candidate algebraic roots and true exceptions

[CERTIFICATE] The cubic in (4) introduces no true exception. In $\mathbb Q[q]/(q^3+2q-1)$,

$$
\left(\frac{q^2+q+3}{2}\right)q(q^2+1)=1.
$$

Therefore the alternate full minor is nonzero at every one of the three conjugate cubic roots. The exact rank is 64 and the nullspace is $\{0\}$ at each root.

[CERTIFICATE] Likewise, the quadratic in (5) introduces no true exception. In $\mathbb Q[q]/(q^2+1)$,

$$
-\frac{q+1}{2}(q^3+2q-1)=1.
$$

Thus the primary full minor is nonzero at both roots of $q^2+1$, again giving rank 64 and nullspace $\{0\}$.

[EXACT EXCEPTIONS] Exact RREF over $\mathbb Q$ gives the following complete specialization table. The two middle columns show the nullities in the parity-preserving and parity-reversing coordinate sectors; their sum is the full nullity.

| $q$ | rank $B(q)$ | even nullity | odd nullity | full nullity | invertible solution |
|---:|---:|---:|---:|---:|:---|
| $0$ | 14 | 26 | 24 | 50 | $R=I_8$, determinant 1 |
| $1$ | 26 | 19 | 19 | 38 | $R=I_8$, determinant 1 |
| $-1$ | 0 | 32 | 32 | 64 | $R=I_8$, determinant 1 |

For each row, the artifact stores every exact nullspace basis vector as a sparse list of coordinate indices and integer coefficients. It also stores the RREF pivot and free columns. The standalone test reconstructs all basis matrices, verifies $B(q)N=0$, verifies independence and rank-nullity 64, rebuilds $I_8$ from the unrestricted coordinate order, and checks both its RLLL residual and determinant exactly.

## Reproduction and scope

[CHECK] Run only the standalone certificate test:

```bash
.venv/bin/python tests/test_tetra_ungraded.py
```

It reconstructs the full exact component matrix, the three sector minors, the two full-minor products, the polynomial gcd and Bezout identity, the algebraic-root quotient certificates, and all three exceptional nullspaces.

[SCOPE] The theorem concerns three identical copies of this fixed two-dimensional-leg Ising $L(q)$ and an arbitrary 8-by-8 $R$ in the displayed vertex-type RLLL equation. It does not classify higher auxiliary dimensions, nonidentical or spectral local factors, IRF/dynamical relations, or singular dimension-changing projections. No claim of solving the 3D Ising model is made.
