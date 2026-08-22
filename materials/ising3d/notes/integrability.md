# Integrability structures tested for the 3-D Ising layer

## Reproduction and scope

Run from the repository root:

```text
.venv/bin/python experiments/e10_conserved_charges.py
.venv/bin/python experiments/e11_tetrahedron.py
.venv/bin/python tests/test_integrability.py
```

The machine-readable records are
`results/integrability/conserved_charges.json` and
`results/integrability/tetrahedron.json`.  The first experiment uses exact finite-field or rational
arithmetic except for the explicitly labelled large transfer-matrix spectra.  Those spectra are
computed independently at 50 and 80 decimal digits.  The tetrahedron specialization is exact over
`Z[w]`.

These calculations test three finite, precisely defined structures.  They do **not** prove that the
3-D Ising model has no exact solution, no useful nonlocal conserved operators, or no other
higher-dimensional integrability representation.

## 1. Translation-invariant local conserved charges

### Ansatz and exact linear system

For a periodic layer, let `O(P)` be the full sum of all translations of a Pauli string `P`.  At
range `r` the ansatz contains one coefficient for every translation orbit `O(P)` for which the
nonidentity sites of `P` are contained in at least one connected lattice animal with at most `r`
vertices.  This includes identities inside the connected cover.  Rotations and reflections are
**not** quotiented, so differently oriented densities have independent coefficients.

The Hamiltonian used in the calculation is

```text
H = -g sum_i X_i - sum_<ij> Z_i Z_j.
```

A packed Pauli string is `Q_(a|b)=X^a Z^b`.  For every ansatz column and every Hamiltonian term the
code uses exactly

```text
[Q_v,Q_w] = ((-1)^(b_v.a_w) - (-1)^(b_w.a_v)) Q_(v+w),
```

whose coefficients are `0` or `+-2`.  Equal output translation orbits are collected, giving a
sparse integer matrix for `Q -> [Q,H]`.  Its rank is computed independently in
`F_2147483647` and `F_2147483629`.  The reported nullities agree in both fields in every case.
This is exact finite-field linear algebra; two-field agreement is a guard against an unlucky
characteristic, not a formal proof that every reported rank equals the rank over `Q`.
For the 2-D rows there is also an exact rational conclusion.  An integer matrix has
`rank_Fp <= rank_Q`, so its rational nullity cannot exceed its finite-field nullity.  The explicit
rational vectors `I` and `H` give rational nullity at least two, while either prime calculation
gives nullity at most two.  Therefore the 2-D kernel is exactly `span_Q{I,H}` over `Q`.  The larger
1-D dimensions in the table are stated as exact dimensions in each listed prime field; agreement
at two primes is the reported lift check.

The containing rings/tori have side `L >= 2r+1`: `L=11` for the range-five chain and `L=9` for the
range-four square layer.  Thus the local supports and the supports created by one commutator do not
alias through the periodic boundary.  The same dimensions were obtained at the two generic fields
`g=2` and `g=3`.

### What “trivial” means here

The trivial local subspace is defined to be

```text
span{ full translation sum of the identity density, H }.
```

It has dimension two for every `r>=2`.  Higher powers of the extensive `H` contain products of
terms at arbitrarily separated sites and are not bounded-support extensive charges in this
ansatz; the local part of the polynomial algebra is therefore the affine span above.  The global
spin-flip operator `prod_i X_i` and lattice translation/reflection permutation operators are
obvious exact symmetries, but they are nonlocal on the chosen `L>r` containers and are absent from
the local Pauli-density ansatz.  Consequently “nontrivial dimension” below means exactly
`kernel dimension - 2`, not quotienting by any unlisted object.

### Exact dimensions

| layer | `r` | ansatz dimension | rank of commutator | kernel dimension | trivial | nontrivial |
|---|---:|---:|---:|---:|---:|---:|
| 1-D periodic chain | 2 | 13 | 10 | 3 | 2 | 1 |
| 1-D periodic chain | 3 | 49 | 44 | 5 | 2 | 3 |
| 1-D periodic chain | 4 | 193 | 186 | 7 | 2 | 5 |
| 1-D periodic chain | 5 | 769 | 760 | 9 | 2 | 7 |
| 2-D periodic square layer | 2 | 22 | 20 | 2 | 2 | 0 |
| 2-D periodic square layer | 3 | 220 | 218 | 2 | 2 | 0 |
| 2-D periodic square layer | 4 | 2461 | 2459 | 2 | 2 | 0 |

The 1-D control is the expected free-fermion tower: enlarging the allowed connected cover by one
site adds two independent conserved charges.  Reproducing the total dimensions `3,5,7` at
`r=2,3,4` (and `9` at `r=5`) validates that the ansatz and commutator equations do not erase the
known integrable structure.

For the 2-D layer, all `2461` translation-orbit coefficients allowed through `r=4` leave only the
two defined trivial directions.  This is a strong local obstruction with a precise cutoff: at
`g=2` and `g=3`, on the stated translation-invariant Pauli ansatz, the rational kernel is exactly
`span_Q{I,H}` through connected-cover range four.  It follows already from either prime rank plus
the two explicit rational kernel vectors, and the second prime independently reproduces the rank.
This says nothing about range five or larger, non-translation-invariant charges, quasilocal
or nonlocal charges, special isolated values of `g`, or charges outside the Pauli/operator algebra
being tested.

## 2. Full commutant of the finite transfer matrix

### Exact matrix construction

For an open layer with `n` sites and internal bond set `E`, the tested transfer operator is

```text
V = exp(K* sum_i X_i) exp(K sum_(ij in E) Z_i Z_j),
tanh(K*) = exp(-2K) = x.
```

The two generic values were `x=1/2` and `x=1/3`, corresponding to
`K=(log 2)/2` and `K=(log 3)/2`.  Although `K` is transcendental rather than rational, the matrix
is made exactly integral up to one irrelevant positive scalar.  If `x=p/q` and
`v=tanh K=(q-p)/(q+p)`, the code replaces each transverse factor by `q I+p X` and each interaction
factor by `v_den I+v_num Z_i Z_j`.  No floating-point matrix entries enter this construction.

Write the resulting matrix as `A D`, where `A` is symmetric positive definite and `D` is positive
diagonal.  It is similar to the real symmetric positive-definite matrix
`sqrt(D) A sqrt(D)`, hence it is diagonalizable.  Therefore

```text
dim Comm(V) = sum_lambda multiplicity(lambda)^2.
```

For dimensions at most 64, SymPy computes the integer characteristic polynomial and its exact
square-free decomposition over `Q`; those commutant dimensions are exact.  For the `n=8` chain and
`3x3` layer, the symmetric matrix is first split by global spin flip and coordinate reflections.
Every symmetry block is diagonalized with `mpmath` twice: 50 digits with relative clustering
threshold `1e-30`, and 80 digits with threshold `1e-45`.  Both precision runs give identical
multiplicity histograms.  These largest two rows are precision-stable numerical multiplicity
classifications, not interval proofs that no closer pair was missed.

### Dimensions at both generic couplings

| open layer | matrix dimension `D` | distinct eigenvalues | multiplicities | `dim Comm(V)` | method |
|---|---:|---:|---|---:|---|
| chain `n=4` | 16 | 16 | `16 x 1` | 16 | exact characteristic polynomial |
| chain `n=6` | 64 | 64 | `64 x 1` | 64 | exact characteristic polynomial |
| chain `n=8` | 256 | 256 | `256 x 1` | 256 | 50/80-digit stable |
| grid `2x2` | 16 | 11 | `8 x 1, 2 x 2, 1 x 4` | 32 | exact characteristic polynomial |
| grid `2x3` | 64 | 58 | `52 x 1, 6 x 2` | 76 | exact characteristic polynomial |
| grid `3x3` | 512 | 378 | `260 x 1, 116 x 2, 2 x 10` | 924 | 50/80-digit stable |

Every row is unchanged between `x=1/2` and `x=1/3`.

The result corrects a tempting but invalid expectation about this test.  The free-fermion chain
has many **local** commuting charges, yet its generic open-chain transfer matrix has simple
spectrum.  The unrestricted commutant then has the minimal possible dimension `D`; every commuting
matrix on that finite space can be interpolated as a polynomial in `V`.  Locality and the behavior
as `n` grows are lost in that interpolation.  Conversely, the square grids have larger full
commutants because of eigenvalue degeneracies (including spatial-symmetry degeneracies), but that
does not establish integrability.  Thus the unrestricted finite-matrix commutant is reported as
requested, but the range-bounded local commutant in Section 1 is the discriminating obstruction.

## 3. Bond-dimension-two tensor and tetrahedron equation

### Local tensor for the partition function

Set `v=tanh K` and introduce `w` with `w^2=v`.  For `sigma=+-1`, define the two-component vector

```text
u_sigma = (1, sigma w).
```

Then an edge contraction is

```text
sum_a u_sigma[a] u_sigma'[a] = 1 + v sigma sigma'
                             = exp(K sigma sigma') / cosh K.
```

Summing the physical spin at each cubic-lattice vertex gives one identical rank-six tensor

```text
T[a,b,c,d,e,f]
  = sum_(sigma=+-1) u_sigma[a]u_sigma[b]u_sigma[c]u_sigma[d]u_sigma[e]u_sigma[f]
  = 2 w^(a+b+c+d+e+f)  if a+b+c+d+e+f is even,
    0                    otherwise.
```

With ordinary Kronecker-delta contractions on every bond,

```text
Z = (cosh K)^number_of_bonds * contraction(product_vertices T).
```

This is an exact local tensor representation of the zero-field isotropic 3-D Ising partition
function with bond dimension two.

The tetrahedron residual code drops the uniform factor `2` from `T`; both sides are homogeneous
of degree four in `R`, so this leaves every equation and every zero unchanged.

### Polynomial tetrahedron system and exact specialization

The chosen operator reshaping is

```text
R[out_x,out_y,out_z ; in_x,in_y,in_z]
  = T[in_x,in_y,in_z,out_x,out_y,out_z].
```

For a general `8x8` operator `R`, the constant Zamolodchikov equation is

```text
R_123 R_145 R_246 R_356 = R_356 R_246 R_145 R_123
```

For six-bit strings `a,b`, define the embedded entries by

```text
(R_ijk)[a,b] = R[a_i,a_j,a_k ; b_i,b_j,b_k]
               * product_(ell not in {i,j,k}) delta[a_ell,b_ell].
```

The polynomial attached to each ordered pair `a,b in {0,1}^6` is explicitly

```text
F[a,b] = sum_(s,t,u in {0,1}^6)
           (R_123)[a,s] (R_145)[s,t] (R_246)[t,u] (R_356)[u,b]
       - sum_(s,t,u in {0,1}^6)
           (R_356)[a,s] (R_246)[s,t] (R_145)[t,u] (R_123)[u,b] = 0.
```

This acts on `(C^2)^tensor_6` and gives `64^2=4096` homogeneous quartic equations in the 64
tensor entries.  In components each side has 64 terms before identical monomials are collected.

After substituting the isotropic Ising tensor, exact integer coefficient propagation leaves 832
nonzero component equations and 18 distinct nonzero univariate polynomials.  Their monic common
generator—the gcd, equivalently the reduced one-variable Groebner-basis generator—is

```text
w^2 (w-1)^3 (w+1)^3 (w^2+1)^2.
```

Hence all 4096 equations hold simultaneously only at
`w in {0, 1, -1, i, -i}`.  A finite ferromagnetic coupling has `0<w<1`, so the raw isotropic tensor
in this reshaping does **not** satisfy the constant tetrahedron equation.  The endpoints include
degenerate infinite- or zero-temperature limits and are not a generic integrable family.

### What happens under gauge transformations

There are two materially different meanings of “gauge” here.

1. **Arbitrary covariant `GL(2)` action on all six legs.**  For `w!=0`,

   ```text
   G = [[1/2,  1/(2w)],
        [1/2, -1/(2w)]]
   ```

   maps `u_+` to `(1,0)` and `u_-` to `(0,1)`.  Thus `G^tensor_6 T` is the six-leg GHZ tensor
   `|000000>+|111111>`.  As an `8x8` operator it is the diagonal projector onto
   `|000>` and `|111>`, so all four embedded copies commute and the tetrahedron equation holds.
   The code verifies its 4096 residual components exactly.

2. **The same gauge on every half-edge while retaining the standard delta bond.**  This requires
   `G^T G=I` up to a scalar.  In the operator reshaping it is a similarity transformation, under
   which satisfaction of the tetrahedron equation is invariant.  Since the raw tensor fails for
   every `0<w<1`, no gauge in this delta-preserving identical-tensor class can repair it.

The nonorthogonal GHZ gauge in item 1 changes the bond metric if applied identically at both ends.
The partition function is preserved only if that inverse metric is kept, or if inverse gauges and
therefore alternating tensors are used on the bipartite cubic lattice.  Consequently the GHZ
normal form is a genuine answer to the broad six-leg `GL(2)`-orbit question, but it is not a
commuting family of the original identical Ising transfer tensors with ordinary delta bonds.

No calculation here classifies alternating bipartite tensors, nontrivial bond metrics, other
input/output pairings, IRF constructions, bond dimension greater than two, or
spectral-parameter-dependent tetrahedron equations.  In particular, bond dimension two is not
proved impossible: its unrestricted gauge orbit even contains the degenerate GHZ solution.  What
is ruled out is only the raw finite-temperature isotropic tensor and every compatible
TE-invariant/delta-preserving identical-tensor gauge of that stated reshaping.

## Bottom line

The 1-D control passes decisively: a range-growing free-fermion charge tower is recovered by exact
Pauli linear algebra.  The corresponding 2-D Ising layer has no charge beyond `I` and `H` through
connected-cover range four in a 2461-dimensional ansatz, at two generic fields and in two prime
characteristics.  The unrestricted finite transfer commutant does not provide an integrability
diagnostic, and the canonical bond-dimension-two site tensor fails the constant tetrahedron
equation throughout the physical finite-temperature interval.  Together these are quantified
obstructions to the two standard 2-D solution mechanisms, not a general no-go theorem for an
exact 3-D solution.
