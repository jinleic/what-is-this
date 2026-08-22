# The Dolan–Grady defect of the Ising transfer operator

**Main result of Track A.** A closed-form, exactly proved expression for the failure of the
Onsager mechanism on an arbitrary layer graph, together with the identification of the failure as
a *quartic* operator whose multiplicity is exactly the number of claws in the frustration graph.

Verification: `experiments/e16_dg_theorem_check.py` — term-by-term integer agreement on 16 graphs
(paths, cycles, stars, complete graphs, trees, square grids, a torus layer, a cubic graph).
Supporting data: `experiments/e14_dolan_grady.py`, `experiments/e15_dg_local_structure.py`.

---

## 1. Setting

For a finite simple graph `Gamma = (V,E)` with `|V| = n`, put

    A = sum_{v in V} X_v ,        B = sum_{(ij) in E} Z_i Z_j .

The layer-to-layer transfer operator of the classical Ising model whose layer graph is `Gamma` is
`V = e^{K^* A} e^{K B}` up to a scalar (see `proofs/algebraic_obstruction.md` §0).  For a 1D layer
`Gamma` is a path or a cycle and the classical model is the 2D Ising model; for a 2D layer
`Gamma` is a square grid and the classical model is the 3D Ising model.

**Dolan–Grady criterion** (Dolan & Grady, Phys. Rev. D **25**, 1587 (1982)).  If two operators
`A, B` satisfy

    (DG1)   [A,[A,[A,B]]] = 16 [A,B],
    (DG2)   [B,[B,[B,A]]] = 16 [B,A],

then they generate an Onsager algebra, the family `A + k B` has an infinite tower of mutually
commuting charges, and the associated transfer matrix is diagonalisable by Onsager's method.
This is the mechanism of Onsager's 1944 solution.

---

## 2. The theorem

> **Theorem DG.** For every finite simple graph `Gamma`:
>
> **(i)** `[A,[A,[A,B]]] - 16[A,B] = 0` identically — (DG1) holds for *every* layer graph.
>
> **(ii)** `[B,[B,[B,A]]] - 16[B,A] = R(Gamma)` with
>
>     R(Gamma) = 24 i * sum_{v in V} Y_v [ (deg v - 2) * sum_{u ~ v} Z_u
>                                        + 2 * sum_{{u,u',u''} subset N(v)} Z_u Z_{u'} Z_{u''} ]
>
> **(iii)** Hence `R(Gamma) = 0` **iff every vertex has degree 0 or 2**, i.e. iff `Gamma` is a
> disjoint union of cycles and isolated vertices; and the **quartic part** of `R(Gamma)`
> vanishes **iff** `Delta(Gamma) <= 2`.
>
> *(The `deg = 0` case is real but physically empty: an isolated site carries a spin decoupled
> from all others, `X_v` commutes with `B`, and `R_v = 0` because both `S_v` and `T_v` are empty
> sums.  It is included here because `experiments/e16_dg_theorem_check.py` exercises edgeless
> vertices and would otherwise contradict a bare "iff 2-regular" claim.  Every layer graph of
> interest has no isolated vertices, so "2-regular" is the operative condition throughout.)*
>
> **(iv)** The number of 4-body Pauli terms in `R(Gamma)` is exactly
> `sum_v C(deg v, 3)`, which is exactly the number of induced claws `K_{1,3}` in the frustration
> graph of the standard generator set.

### 2.1 Proof of (i)

`[A,B] = sum_{(ij) in E} [A, Z_iZ_j]`, and only `X_i, X_j` fail to commute with `Z_iZ_j`, so the
computation is edge-local.  With `[X,Z] = -2iY`, `[X,Y] = 2iZ`:

    [A, Z_iZ_j] = [X_i,Z_i]Z_j + Z_i[X_j,Z_j] = -2i ( Y_i Z_j + Z_i Y_j )
    [A, Y_iZ_j] = [X_i,Y_i]Z_j + Y_i[X_j,Z_j] = 2i Z_iZ_j - 2i Y_iY_j
    [A, Z_iY_j] = -2i Y_iY_j + 2i Z_iZ_j
    => [A,[A,Z_iZ_j]] = -2i( 4i Z_iZ_j - 4i Y_iY_j ) = 8 Z_iZ_j - 8 Y_iY_j
    [A, Y_iY_j] = 2i Z_iY_j + 2i Y_iZ_j
    => [A,[A,[A,Z_iZ_j]]] = 8(-2i)(Y_iZ_j+Z_iY_j) - 8(2i)(Z_iY_j+Y_iZ_j)
                          = -32i (Y_iZ_j + Z_iY_j)
                          = 16 * [A, Z_iZ_j].

Summing over edges gives (i). `[]`

*(Interpretation: (DG1) is automatic because `A` is a sum of commuting single-site terms and `B`
is a sum of two-site terms; it carries no information about `Gamma`.)*

### 2.2 Proof of (ii)

By linearity in `A`, `R = sum_v R_v` with `R_v = [B,[B,[B,X_v]]] - 16[B,X_v]`.  Fix `v`, let
`d = deg v` and `N(v) = {u_1,...,u_d}`.  Write `S = sum_k Z_{u_k}`,
`T = sum_{k<l<m} Z_{u_k}Z_{u_l}Z_{u_m}`.

**Step 1.** Only the `d` edges at `v` fail to commute with `X_v`, and
`[Z_vZ_u, X_v] = (Z_vX_v - X_vZ_v)Z_u = -2 X_vZ_v Z_u = 2i Y_v Z_u`, so

    [B, X_v] = 2i Y_v S.                                                     (1)

**Step 2.** Compute `[B, Y_v Z_{u_k}]`.  An edge `Z_aZ_b` with `a,b != v` commutes with
`Y_vZ_{u_k}` (it is diagonal, and it carries `I` at `v`).  An edge `Z_{u_k}Z_w`, `w != v`,
also commutes.  Only the edges `Z_vZ_{u_l}` contribute, and `Z_v` anticommutes with `Y_v` while
`Z_{u_l}` commutes with `Z_{u_k}`, so the two operators anticommute and the commutator is twice
the product.  Using `Z_vY_v = -i X_v` and `Z_{u_k}^2 = 1`:

    [B, Y_vZ_{u_k}] = -2i X_v ( 1 + sum_{l != k} Z_{u_l} Z_{u_k} ).

Hence, with (1),

    [B,[B,X_v]] = 2i * (-2i) X_v sum_k ( 1 + sum_{l!=k} Z_{u_l}Z_{u_k} )
                = 4 X_v ( d + 2 sum_{k<l} Z_{u_k}Z_{u_l} ).                  (2)

**Step 3.** Compute `[B, X_v Z_{u_k} Z_{u_l}]` (`k<l`).  Again only the edges `Z_vZ_{u_m}`
contribute; `Z_v` anticommutes with `X_v`, `Z_{u_m}` commutes with the `Z`'s.  Using
`Z_vX_v = i Y_v`:

    [B, X_vZ_{u_k}Z_{u_l}] = 2i Y_v ( Z_{u_k} + Z_{u_l} + sum_{m notin {k,l}} Z_{u_m}Z_{u_k}Z_{u_l} ).

Applying `[B, .]` to (2) and using
`sum_{k<l}(Z_{u_k}+Z_{u_l}) = (d-1) S` and
`sum_{k<l} sum_{m notin{k,l}} Z_{u_m}Z_{u_k}Z_{u_l} = 3 T`:

    [B,[B,[B,X_v]]] = 4d * (2i Y_v S) + 8 * 2i Y_v ( (d-1) S + 3 T )
                    = 8i Y_v ( (3d - 2) S + 6 T ).

**Step 4.** Subtract `16[B,X_v] = 32 i Y_v S`:

    R_v = 8i Y_v ( (3d-2-4) S + 6T ) = 24 i Y_v ( (d-2) S + 2 T ).           `[]`

### 2.3 Proof of (iii) and (iv)

The Pauli strings `Y_v Z_u` (over `u ~ v`) and `Y_v Z_uZ_{u'}Z_{u''}` (over triples in `N(v)`) are
pairwise distinct for distinct `v` as well, because the unique `Y` locates `v`.  Hence there is no
cancellation between different `v`, and

    #terms(R_v) = d * [d != 2] + C(d,3),        #4-body terms(R_v) = C(d,3).

So `R = 0` iff every vertex has `d*[d!=2] + C(d,3) = 0`, i.e. `d in {0,2}`; the quartic part
vanishes iff every `C(d,3) = 0`, i.e. `d <= 2`.
By Theorem 5 of `proofs/algebraic_obstruction.md` the frustration graph of the standard generator
set is the subdivision `S(Gamma)`, whose induced claws are exactly the triples of edges meeting at
a common vertex; their number is `sum_v C(deg v, 3)`. `[]`

*(Global phase: the machine check uses the non-Hermitian basis `Q_v = X^aZ^b`, in which the whole
expression carries an extra factor `-1` relative to the Hermitian convention above.  This is
recorded by `e16_dg_theorem_check.py`, which prints the phase and then demands exact integer
agreement of every coefficient.)*

---

## 3. What this says about 2D versus 3D

Define the **quartic defect density**

    q(Gamma) := (1/n) * sum_{v} C(deg v, 3) .

| classical model | layer graph | `Delta` | `q` (bulk) | `R = 0`? |
|---|---|---|---|---|
| 1D Ising | single site, no bonds | 0 | 0 | yes (trivially) |
| **2D Ising**, cylinder | cycle `C_n` | 2 | **0** | **yes, exactly** |
| **2D Ising**, strip | path `P_n` | 2 | **0** | no, but only 2 boundary terms, no quartic part |
| **3D Ising**, `L x L` layer | square grid | 4 | **`4 (L-2)^2/L^2 + O(1/L) -> 4`** | no |
| **3D Ising**, `L x L` torus layer | 4-regular grid | 4 | **exactly 4** | no |

* For a **1D layer the quartic defect is identically zero**, for every boundary condition and
  every length.  The only failure of (DG2) is the two `(d-2) = -1` boundary terms of an open
  chain: an `O(1)` surface effect that is removed by closing the chain, and that does not obstruct
  free-fermionisation (the frustration graph of a path is a path, a line graph).
* For a **2D layer the quartic defect is extensive**: exactly `4` per bulk site — one for each of
  the `C(4,3) = 4` triples of neighbours.  It cannot be removed by any boundary condition, and it
  survives every thermodynamic limit.

**This is the precise, quantitative answer to "what fails in three dimensions".**  Not the
dimension of the lattice as such, and not "the algebra becomes complicated": the obstruction is
the coordination number of the *layer*, entering through a single binomial coefficient
`C(Delta, 3)`.  The obstruction first appears — and can be seen in a system of four spins — the
moment one vertex acquires a third neighbour.

Three logically independent characterisations coincide **exactly**:

    Gamma has a vertex of degree >= 3
      <=>  the frustration graph of {X_i} u {Z_iZ_j} contains an induced claw
           (Theorem 5, proofs/algebraic_obstruction.md)
      <=>  no injective term-wise map to Majorana bilinears exists
           (Theorem 4 here + Chapman-Flammia Thm 1)
      <=>  the Dolan-Grady residual has a nonzero QUARTIC part
           (Theorem DG (iii))

and the multiplicities agree: *number of quartic obstruction terms = number of claws*.

---

## 4. Locality of the obstruction in the Pauli presentation

Theorem DG (ii) identifies the *first* obstruction in the Dolan-Grady residual: it is **4-local
in the Pauli presentation**, supported on a site and three of its neighbours, and there are
exactly `sum_v C(deg v,3)` such distinct 4-site Pauli monomials `Y_v Z_{u}Z_{u'}Z_{u''}`.  No
cubic or 5-local-or-higher Pauli term appears in the residual `R` at this order, and the count is
exact.

**What this does and does not establish -- read carefully, this was an overclaim in an earlier
draft.**

* It DOES establish that the DG2 residual is 4-local in the Pauli basis, with the stated count,
  for the natural generators `A = sum X_v`, `B = sum Z_iZ_j`.
* It does NOT establish that the "minimal fermionic extension" is quartic.  A 4-site Pauli string
  is not generally a four-Majorana term: under a Jordan-Wigner map the string operator attached to
  distant sites can involve many Majoranas, and a different Majorana choice changes the count
  again.  Pauli 4-locality and fermionic 4-body are distinct notions.
* It does NOT establish that adding quartics *closes* the algebra, nor that the `sum C(deg,3)`
  monomials are independent as generators, nor that they are minimal.  Proving minimality would
  require a closure construction showing that the bilinear-plus-quartic span is closed under
  commutation and reproduces the transfer operator, or a lower bound on the number of
  generators of any extension.  Neither is provided here; the spectral result of
  `proofs/spectral_gaussianity.md` is the independent, basis-independent obstruction.

**Caveat, stated explicitly.**  Theorem DG is a statement about the *Onsager/Dolan–Grady route*.
Failure of (DG2) proves that this particular mechanism does not apply; it does not by itself prove
that no exact solution exists.  Two independent facts limit how far it can be pushed:

1. Open 1D chains also violate (DG2) (by boundary terms) yet are exactly solvable — so a nonzero
   `R` is not sufficient for non-solvability.  This is why the *quartic* part, which is exactly
   zero for every 1D layer including open ones, is the meaningful invariant.
2. Highly symmetric graphs with claws can still have small algebras: `experiments/e13` finds
   `dim <A,B>_Lie = 24` for the star `K_{1,3}` and `15` for `K_4`, versus `2952` for the `2x4`
   grid.  A large quartic defect density is therefore evidence, not proof, of an exponentially
   large algebra; whether that growth is exponential is an open question (sec. 7.5 of the
   final report).
