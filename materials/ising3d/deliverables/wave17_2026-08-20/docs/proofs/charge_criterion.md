# Graph criterion for the exceptional proper-support charges

## 1. Scope and conventions

**[THEOREM — finite geometry, not a solution of the 3D model].** This note
classifies the four exceptional support-orbit representatives found at size
seven on the open `3x3` layer, gives conceptual formulas for every element of
their exact kernels, and explains why the same graph mechanism cannot occur on
a size-seven support of the open `3x4` layer. It does not claim an
all-rectangle classification, a thermodynamic conserved charge, or a solution
of the three-dimensional Ising model.

For a finite simple graph $G=(V,E)$, write

$$
 A=\sum_{v\in V}X_v,\qquad B=\sum_{\{u,v\}\in E}Z_uZ_v.
$$

An operator supported in $S\subseteq V$ means $O_S\otimes I_{V\setminus S}$.
The exact Pauli convention used by the calculation is
$Q_{a,b}=X^aZ^b$. All identities below are identities over $\mathbb Q$.

## 2. The pointwise-stabilizer theorem

**Theorem 1 [THEOREM].** Let $\sigma\in\operatorname{Aut}(G)$ fix every
vertex of $V\setminus S$. Let $U_\sigma$ be the unitary that permutes tensor
factors according to $\sigma$. Then $U_\sigma$ is supported in $S$ and

$$
 [U_\sigma,A]=[U_\sigma,B]=0.
$$

If $\sigma$ is nonidentity, then $U_\sigma$ is non-scalar.

**Proof.** Tensor-factor permutation gives

$$
 U_\sigma X_vU_\sigma^{-1}=X_{\sigma(v)},\qquad
 U_\sigma Z_vU_\sigma^{-1}=Z_{\sigma(v)}.
$$

The first sum is invariant because $\sigma$ permutes $V$; the second is
invariant because it permutes $E$. Pointwise fixation of $V\setminus S$
means the tensor permutation factors as an operator on $S$ tensored with the
identity outside. A nonidentity vertex permutation moves at least one
computational-basis vector, so its tensor-permutation matrix is not scalar.
$\square$

The theorem yields an exact graph-theoretic criterion for this mechanism:

$$
 \boxed{\operatorname{Aut}(G)_{(V\setminus S)}\ne 1.}
 \tag{1}
$$

Here the subscript denotes the pointwise stabilizer. This criterion is
necessary and sufficient for a nonidentity *site-permutation* charge; it is
only conjectured to be necessary for an arbitrary local charge outside the
finite budgets stated below.

**Lemma 2 [LEMMA — boundary-neighborhood form].** A permutation
$\pi\in\operatorname{Sym}(S)$ extends by the identity on $V\setminus S$ to an
automorphism of $G$ if and only if

1. $\pi$ is an automorphism of the induced graph $G[S]$, and
2. for every $s\in S$,
   $$N(s)\cap(V\setminus S)=N(\pi(s))\cap(V\setminus S).$$

**Proof.** Internal edges are preserved exactly by condition 1. Edges with
one endpoint outside are preserved exactly by condition 2. Edges with both
endpoints outside are unchanged pointwise. These exhaust all edges. $\square$

Thus equal external neighborhoods are the relevant equitable cells, but an
individual pair of vertices need not be twins in the full graph: the
permutation must also preserve all internal edges simultaneously.

## 3. The four `3x3` support orbits

Use row-major labels

```text
0 1 2
3 4 5
6 7 8
```

The exact `flux_support` kernels and their pointwise stabilizers are:

| support mask | omitted vertices | reflection cycles | supported fixed connector $f$ | orbit size | exact kernel dimension |
|---:|---|---|---:|---:|---:|
| 239 | $\{4,8\}$ | $(1\,3)(2\,6)(5\,7)$ | 0 | 4 | 4 |
| 254 | $\{0,8\}$ | $(1\,3)(2\,6)(5\,7)$ | 4 | 2 | 4 |
| 367 | $\{4,7\}$ | $(0\,2)(3\,5)(6\,8)$ | 1 | 4 | 2 |
| 381 | $\{1,7\}$ | $(0\,2)(3\,5)(6\,8)$ | 4 | 2 | 2 |

Each pointwise stabilizer has order two. The first two rows use reflection in
a diagonal; the last two use reflection in an axial line. In every row the
six moved vertices alone induce two disconnected components. Adding one fixed
vertex joins them, so seven is the first connected support size that can
contain the permutation. The four representatives expand to twelve actual
size-seven supports under the full order-eight square automorphism group.

For every reflection pair, the two vertices have equal neighborhoods outside
the displayed support, as required by Lemma 2. None of the three individual
transpositions is a graph automorphism. The conserved permutation is the
simultaneous product of all three swaps. This rules out the naive explanation
as three independent twin swaps.

## 4. Alternating Bell-cell construction

The diagonal-reflection kernels have two additional generators. They come
from a general quotient-graph construction.

For a two-cycle $C_i=\{u_i,v_i\}$ of an involutive graph automorphism, define

$$
 |\beta_{-}\rangle=\frac{|01\rangle-|10\rangle}{\sqrt2},\qquad
 |\beta_{+}\rangle=\frac{|00\rangle-|11\rangle}{\sqrt2}.
$$

They obey the two exact identities

$$
 (X_{u_i}+X_{v_i})|\beta_{\eta_i}\rangle=0,
 \qquad
 Z_{v_i}|\beta_{\eta_i}\rangle=\eta_i Z_{u_i}|\beta_{\eta_i}\rangle,
 \quad \eta_i\in\{-1,+1\}.
 \tag{2}
$$

**Theorem 3 [THEOREM — alternating Bell cells].** Let $\sigma$ be an
involutive graph automorphism, with fixed vertices $F$ and two-cycle cells
$C_i$. Suppose signs $\eta_i\in\{-1,+1\}$ can be chosen so that

1. $\eta_i=-1$ whenever a fixed vertex is adjacent to $C_i$; and
2. $\eta_i\eta_j=-1$ whenever an edge orbit joins distinct cells $C_i,C_j$.

Let

$$
 |\Omega\rangle=\bigotimes_i|\beta_{\eta_i}\rangle,
 \qquad P=|\Omega\rangle\langle\Omega|
$$

on the moved vertices, with identity on all fixed vertices. Then
$[P,A]=[P,B]=0$. Moreover, if a fixed vertex $f$ has no fixed-fixed neighbor
and all cells adjacent to $f$ have sign $-1$, then

$$
 [X_fP,A]=[X_fP,B]=0.
$$

**Proof.** Equation (2) makes the transverse-field sum on every cell
annihilate $|\Omega\rangle$. A fixed-to-cell edge orbit contributes
$Z_w(Z_{u_i}+Z_{v_i})$, which annihilates the Bell factor when
$\eta_i=-1$. A matching intercell edge orbit contributes

$$
 Z_{u_i}Z_{u_j}+Z_{v_i}Z_{v_j},
$$

which acts on $|\Omega\rangle$ as
$(1+\eta_i\eta_j)Z_{u_i}Z_{u_j}=0$. A crossed orbit instead gives the factor
$\eta_i+\eta_j=0$. An edge inside one cell has a fixed Bell eigenvalue, and
edges among fixed vertices act only on identity tensor factors. Consequently
the range of $P$ is reducing for both $A$ and $B$, proving the first claim.

$X_f$ commutes with $A$. Every $B$ term not incident to $f$ also commutes
with it. By the extra assumptions, the terms incident to $f$ occur only as
$Z_f(Z_u+Z_v)$ on negative cells and annihilate $P$. This proves the second
claim. $\square$

### Why diagonal and axial reflections differ

For diagonal reflection the quotient of the three two-cycle cells is a path

$$
 \{1,3\}\;--\;\{2,6\}\;--\;\{5,7\}.
$$

Only the two endpoint cells touch fixed vertices. The unique assignment is

$$
 (\eta_{13},\eta_{26},\eta_{57})=(-1,+1,-1).
 \tag{3}
$$

Every diagonal fixed vertex has no fixed-fixed neighbor, so the supported
connector $f=0$ or $f=4$ satisfies Theorem 3.

For axial reflection the quotient cells are
$\{0,2\}--\{3,5\}--\{6,8\}$. Each cell touches one of the fixed vertices
$1,4,7$, so condition 1 forces all three signs to $-1$, while each quotient
edge forces opposite signs. The system is inconsistent. This is the precise
boundary obstruction to the extra Bell projector on masks 367 and 381.

## 5. Minimal exact operator forms

For a pair $u,v$, let

$$
\begin{aligned}
 W_{uv}&=\tfrac12\bigl(I+X_uX_v+Z_uZ_v-(X_uZ_u)(X_vZ_v)\bigr),\\
 S_{uv}&=\tfrac14\bigl(I-X_uX_v-Z_uZ_v+(X_uZ_u)(X_vZ_v)\bigr),\\
 F_{uv}&=\tfrac14\bigl(I-X_uX_v+Z_uZ_v-(X_uZ_u)(X_vZ_v)\bigr).
\end{aligned}
 \tag{4}
$$

In the ordered $X^aZ^b$ convention, $W$ is the swap, $S$ projects onto
$|\beta_-\rangle$, and $F$ projects onto $|\beta_+\rangle$.

### Diagonal supports: masks 239 and 254

Set

$$
 U=W_{13}W_{26}W_{57},\qquad
 P=S_{13}F_{26}S_{57},
 \tag{5}
$$

and take $f=0$ for mask 239 or $f=4$ for mask 254. The four basis vectors
stored by `flux_support` are exactly, coefficient by coefficient,

$$
 \boxed{I,\quad 4U-32P,\quad I-4U-32P,\quad -64X_fP.}
 \tag{6}
$$

Their ordered-Pauli term counts are respectively $1,28,35,64$, exactly the
stored counts. Direct Fraction substitution in the symbolic formulas for
$\operatorname{ad}_A$ and $\operatorname{ad}_B$ gives empty residuals.

A more structural basis consists of four mutually orthogonal projectors:

$$
\begin{aligned}
 E_-&=\frac{I-U}{2},\\
 E_{+,0}&=\frac{I+U}{2}-P,\\
 E_{X,+}&=\frac{P+X_fP}{2},\\
 E_{X,-}&=\frac{P-X_fP}{2}.
\end{aligned}
 \tag{7}
$$

Their ranks on the seven-site Hilbert space are

$$
 56,\quad 70,\quad 1,\quad 1.
$$

They are symmetric idempotents, are pairwise orthogonal, sum to identity, and
commute with $A,B$. After tensoring identity on the two omitted vertices they
are therefore actual projectors in the full `3x3` joint commutant, with global
ranks $224,280,4,4$. Thus the exceptional basis is not merely reminiscent of
global commutant projectors: it is an exact linear combination of explicit
such projectors. In projector order $(E_-,E_{+,0},E_{X,+},E_{X,-})$, the four
rows of (6) have eigenvalue vectors

$$
 (1,1,1,1),\quad(-4,4,-28,-28),\quad(5,-3,-35,-35),\quad(0,0,-64,64).
$$

The local kernel algebra is consequently $\mathbb C^4$. The four independent
operators in (6), together with the prior exact kernel dimension four, exhaust
the kernel.

### Axial supports: masks 367 and 381

Here

$$
 U=W_{02}W_{35}W_{68}.
$$

The exact stored basis is simply

$$
 \boxed{I,\quad I-8U.}
 \tag{8}
$$

It has term counts $1,63$. The projectors $(I+U)/2$ and $(I-U)/2$ have local
ranks $72$ and $56$ (global ranks $288$ and $224$), and (8) has eigenvalues
$(-7,9)$ on these two blocks. Hence this kernel algebra is $\mathbb C^2$.
The prior exact dimension two proves that the inconsistent Bell assignment did
not hide another local charge.

## 6. Exact rejection of an induced-$C_4$ explanation

**[THEOREM — finite exact statements].** An induced square is neither
necessary nor sufficient for these charges.

- Masks 239, 367, and 381 contain no induced $C_4$ at all, yet have exact
  non-scalar kernels.
- Mask 27, the single induced square on vertices $\{0,1,3,4\}$, has exact
  kernel dimension one.
- Mask 254 contains induced squares $\{1,2,4,5\}$ and $\{3,4,6,7\}$. For each
  square, exact rational elimination of coefficients outside the square shows
  that the intersection of its operator space with the four-dimensional
  kernel has dimension one: identity only.

Thus none of the nonidentity support-seven charges is an operator localized to
an induced square.

## 7. Why `3x4` has no size-seven instance

The open `3x4` rectangle has automorphism group $C_2\times C_2$. Its three
nonidentity elements move respectively $8,12,12$ vertices. The element moving
eight vertices exchanges the top and bottom rows and fixes the middle row.
Its moved set is the disjoint union of two four-vertex paths, so any connected
support containing it must include at least one fixed middle-row connector.
Therefore the smallest connected support satisfying (1) has size nine, not
seven. There are exactly four such size-nine supports.

**[COMPUTATION — exhaustive graph search].** Among all 234 connected
size-seven subsets of the `3x4` grid, none has a nontrivial complement
pointwise stabilizer. The stored table quotients them to 62 orbits, and its
exact-Q kernel dimension is one for every orbit. The first statement excludes
the proved automorphism mechanism; the second exact computation also excludes
any different local mechanism on those supports.

The same automorphism-only search, without constructing any $4^{|S|}$ kernel,
gives:

| grid | nonidentity moved-set sizes | smallest connected criterion support | criterion-positive size-seven supports |
|---|---|---:|---:|
| `2x4` | $8,8,8$ | 8 | 0 |
| `3x3` | $6,6,6,6,8,8,8$ | 7 | 12 |
| `3x4` | $8,12,12$ | 9 | 0 |
| `4x4` | $12,12,16,16,16,16,16$ | 13 | 0 |
| `3x5` | $10,12,14$ | 11 | 0 |

There are eight minimum supports on `4x4` and five on `3x5`; Theorem 1 gives
a permutation charge on each. No claim about their complete local kernels is
made.

## 8. Finite converse budget on general graphs

**[COMPUTATION].** The experiment enumerated every connected simple graph on
$2\le n\le6$ up to isomorphism and every nonempty proper connected support up
to its graph automorphism group. The budget is:

| $n$ | connected graph classes | proper connected support orbits | nontrivial pointwise stabilizer | trivial pointwise stabilizer |
|---:|---:|---:|---:|---:|
| 2 | 1 | 1 | 0 | 1 |
| 3 | 2 | 5 | 1 | 4 |
| 4 | 6 | 29 | 9 | 20 |
| 5 | 21 | 220 | 74 | 146 |
| 6 | 112 | 2427 | 774 | 1653 |
| **total** | **142** | **2682** | **858** | **1824** |

For each support, the integer ordered-Pauli commutator matrix was reduced at
primes $2147483647$ and $2147483629$. The nullities agreed at both primes. In
every trivial-stabilizer case the modular nullity was one. A nonzero
$(4^{|S|}-1)$-minor modulo either prime is a nonzero integer minor, so the
rational nullity is at most one; explicit identity makes it exactly one. In
every positive case Theorem 1 supplies two rationally independent kernel
elements, $I$ and $U_\sigma$. Hence, on this entire finite budget,

$$
 \ker_{\mathbb Q}(\operatorname{ad}_A\oplus\operatorname{ad}_B)
 \text{ is non-scalar}
 \iff
 \operatorname{Aut}(G)_{(V\setminus S)}\ne1.
$$

This is an exact finite conclusion, not an extrapolation from floating point.
The complete ordered case stream has SHA-256
`7cba88cf79625ade4fc992d5f1ba26a8d176290049c1f726cbc5cc6300c0fb19`.

**[CONJECTURE].** The equivalence above holds for every finite connected
simple graph and every nonempty proper connected support. Theorem 1 proves the
forward construction, but necessity beyond the explicit finite graph and grid
budgets remains **[UNRESOLVED]**.

## 9. Reproducibility

- Experiment: `experiments/e61_charge_criterion.py`
- Exact artifact: `results/integrability/charge_criterion.json`
- Independent test: `tests/test_charge_criterion.py`
- Source exact bases: `results/integrability/flux_support.json`

The experiment finishes with `PASS` only when every embedded exact check
passes. The independent test reconstructs the conceptual operators without
importing the experiment, substitutes all of them in both commutators, checks
the projector algebras, verifies criterion equivalence on all 250 stored orbit
rows, repeats every larger-grid automorphism search, and independently repeats
the general-graph scan through five vertices.
