# Spectral tetrahedron and commuting-transfer tests

## Scope and reproduction

```text
.venv/bin/python experiments/e35_tetrahedron_spectral.py
.venv/bin/python tests/test_tetrahedron_spectral.py
```

The machine record is `results/integrability/tetrahedron_spectral.json`.  All arithmetic is exact: Python integers and `Fraction`, and SymPy polynomial algebra over `QQ`.  No external value of `K_c` enters any construction or selection.

[COMPUTATION] The frozen Wave-2 calculation in `experiments/e11_tetrahedron.py`, `notes/integrability.md`, and `results/integrability/tetrahedron.json` tested the **constant**, spectral-parameter-free equation

```text
R_123 R_145 R_246 R_356 = R_356 R_246 R_145 R_123
```

for one reshaping of the one-parameter isotropic bond-dimension-two Ising site tensor.  It found 832 nonzero components among 4096; their raw common factor was

```text
w^2 (w-1)^3 (w+1)^3 (w^2+1)^2,
```

so no raw tensor with `0<w<1` passed.  It also showed that a covariant, nonorthogonal `GL(2)` action maps a fixed tensor to GHZ but changes the delta bond metric, while a delta-preserving identical-tensor gauge cannot alter TE satisfaction.  It explicitly did **not** classify spectral dependence, alternating tensors, changed metrics, auxiliary bond dimensions, IRF representations, or other input/output pairings.  The work below is a strict extension only in the explicitly stated six-rapidity parity-tensor and two-weight transfer families; it is not a classification of general 64-entry `R` matrices.

## Verified vertex/tensor representation

[LEMMA] Let `v=tanh K`, choose `w` with `w^2=v`, and define

```text
u_sigma = (1, sigma w),                  sigma in {-1,+1}.
```

Then

```text
sum_a u_sigma[a] u_tau[a]
  = 1 + v sigma tau
  = exp(K sigma tau) / cosh K.
```

Summing a physical spin at every cubic-lattice site produces the explicit rank-six tensor

```text
T[a,b,c,d,e,f]
  = sum_{sigma=+-1} u_sigma[a]u_sigma[b]u_sigma[c]
                         u_sigma[d]u_sigma[e]u_sigma[f]
  = 2 w^(a+b+c+d+e+f),  if a+b+c+d+e+f is even,
    0,                   otherwise.
```

Consequently, with an ordinary Kronecker delta on every bond,

```text
Z(K) = (cosh K)^E contraction(product_sites T).
```

[COMPUTATION] Before using this representation in any integrability test, the experiment contracted its binary legs exactly.  The contraction dynamic program retains the vertex-parity mask and selected-edge degree after every edge.  Thus it computes the tensor polynomial independently of spin enumeration.  It agrees coefficient by coefficient with `ising.exact_enumeration.even_subgraph_polynomial` on three open boxes:

| box | sites | bonds | exact tensor/enumerator polynomial `P(v)` | common scaled partition integer at `v=1/3` |
|---|---:|---:|---|---:|
| `2x2x2` | 8 | 12 | `[1,0,0,0,6,0,16,0,9,0,0,0,0]` | `36450` |
| `2x2x3` | 12 | 20 | `[1,0,0,0,11,0,36,0,99,0,184,0,157,0,20,0,4,0,0,0,0]` | `16394562` |
| `2x3x3` | 18 | 33 | `[1,0,0,0,20,0,78,0,402,0,1640,0,4909,0,11636,0,17685,0,16840,0,9374,0,2542,0,376,0,32,0,1,0,0,0,0,0]` | `246853161090` |

The integer column compares the density-of-states polynomial at `exp(2K)=2` with the tensor expression `2^N P(1/3)` multiplied by `(3/2)^E`; both sides are exact integers.  This control fixes the tensor normalization and lattice bond conventions before the TE calculation.

## Six-space edge-rapidity tetrahedron equation

### Ansatz and equations

[COMPUTATION] Give each of the six shared TE vector spaces an independent multiplicative edge rapidity:

```text
T_[w1,...,w6](a1,...,a6)
  = 2 product_i w_i^a_i,  if sum_i a_i is even,
    0,                     otherwise,

R[out1,out2,out3 ; in1,in2,in3]
  = T[in1,in2,in3,out1,out2,out3].
```

The common factor 2 may be dropped because each TE side is homogeneous of degree four.  Put `x_i=w_i^2`.  The embedded `R_ijk` uses the variables attached to spaces `i,j,k` on both its input and output legs.  This is a strict spectral extension of the constant tensor and spans all six shared-space edge rapidities while retaining the verified parity support.  It is **not** the general equation with four independently parametrized `R_123,R_145,R_246,R_356`, nor a classification of arbitrary 64-entry `R` matrices.

For each `alpha,beta in {0,1}^6`, the exact component equation is

```text
F[alpha,beta]
 = sum_{s,t,u}
     (R_123)[alpha,s](R_145)[s,t](R_246)[t,u](R_356)[u,beta]
   - sum_{s,t,u}
     (R_356)[alpha,s](R_246)[s,t](R_145)[t,u](R_123)[u,beta]
 = 0.
```

[COMPUTATION] Exact propagation gives 1792 nonzero component polynomials.  On the open torus where every `w_i` is nonzero, each equation's monomial factor is invertible; removing it and replacing even powers `w_i^(2k)` by `x_i^k` leaves 28 distinct polynomials supported on 64 monomials.  Exact rational row reduction has rank seven.  The seven retained independent equations are recorded verbatim in the JSON artifact; no modular-rank inference is used.

### Elimination certificate and the manifold

[COMPUTATION] Let `I` be the ideal generated by those seven equations.  To ask for a branch on which no leg is exceptional, the experiment adds

```text
q product_i x_i(x_i-1)(x_i+1) - 1 = 0
```

and computes a Gröbner basis over `QQ[q,x1,...,x6]` in graded reverse lexicographic order.  The exact basis is

```text
[1].
```

Therefore no solution exists with all six `x_i` outside `{0,+1,-1}`.  Every solution in this six-rapidity ansatz has at least one zero-weight or zero-/infinite-temperature/complex exceptional leg.  This is the precise nonexceptional-manifold no-go; zero branches themselves are not classified by the monomial-stripped ideal.

[COMPUTATION] The variety is not globally empty.  For example, direct substitution into all seven independent equations verifies the three-dimensional degenerate component

```text
x4=1, x5=1, x6=-1,        x1,x2,x3 arbitrary.
```

This fixes three tensor legs at zero/infinite-temperature or complex exceptional values.  It is an algebraic TE family but not a finite ferromagnetic Ising spectral family.

On the physical diagonal `x_1=...=x_6=t=tanh K`, the reduced polynomial gcd is

```text
(t-1)^3 (t+1)^2,
```

and the raw residual gcd in `w` reproduces the frozen constant result exactly:

```text
w^2 (w-1)^3 (w+1)^3 (w^2+1)^2.
```

[COMPUTATION] Hence no finite ferromagnetic physical weight, `0<t<1`, lies on this manifold.

## The three proposed escapes

### (a) Tensor-leg gauge transformation

[LEMMA] A partition-preserving delta-bond gauge uses paired changes of basis `G` and `G^{-T}` on the two half-edges of every contracted bond.  In the displayed operator reshaping, the four embedded `R` operators undergo a simultaneous global tensor-product similarity.  The TE residual is conjugated by an invertible matrix, so its vanishing is invariant.

[COMPUTATION] Since the physical diagonal residual is nonzero for every `0<t<1`, this gauge class **fails** to move the physical tensor onto the TE variety.

[COMPUTATION] The broader covariant action already exhibited in Wave 2 does map each fixed `w != 0` tensor to GHZ.  It is not a counterexample to the preceding result: applied identically on both ends it changes the delta bond metric, and retaining the partition function requires inverse metrics or alternating inverse gauges.  It therefore does not produce an identical-tensor spectral family for the original vertex model.

### (b) Auxiliary-state extension

[COMPUTATION] The minimal explicit extension tested enlarges each local leg from dimension two to three and defines

```text
R_3 = R_Ising on {0,1}^3,
R_3 = identity on a local triple containing auxiliary state 2,
no mixing between the two blocks.
```

At the exact physical point `w=1/2` (`t=1/4`), its TE component with output state 20 and input state 0 is

```text
-675/4096.
```

The physical two-state subspace is invariant, so the old nonzero residual is a principal block of the extended residual.  This direct-sum auxiliary construction **fails**.

[UNDECIDED] General interacting bond-dimension-three or larger extensions, in which auxiliary states mix with the physical block and the physical partition function is recovered by a projection or trace, are not classified.  The explicit direct-sum failure is not promoted to a universal auxiliary-state no-go.

### (c) Controlled limiting procedure

[LEMMA] A polynomial zero set is closed under coefficientwise finite limits.  Thus, if exact TE tensors in this six-rapidity parity ansatz converged coefficientwise to a finite isotropic tensor with `0<t<1`, that limiting tensor would itself obey all TE component equations.

[COMPUTATION] The physical-diagonal gcd excludes such a tensor, so every regular finite coefficientwise limiting procedure in this ansatz **fails**.  Exact solutions occur only on components with some `x_i=+1` or `-1`; these exceptional degenerations do not approach an interior finite-ferromagnetic point regularly.

[UNDECIDED] Singular limits using divergent gauges, dimension-changing projections, or cancellations between divergent auxiliary blocks are outside the polynomial-closedness argument.

## Commuting layer-to-layer transfer matrices

[COMPUTATION] On a periodic `L_x x L_y` layer, including both parallel bonds when a side is two, define the exact two-weight matrix (up to a scalar)

```text
T(a,b)[s,t] = a^Hamming(s,t) b^broken_periodic_layer(t).
```

For the Ising transfer operator, `a=e^{-2K_z}` and `b=e^{-2K_parallel}`.  The isotropic physical curve is

```text
b = (1-a)/(1+a),       0<a<1,
```

when `a=tanh K*` and `b=tanh K` in the conventional factorized normalization.

[COMPUTATION] Two exact untuned points on this curve, `(a,b)=(1/2,1/3)` and `(1/3,1/2)`, do not commute.  The first nonzero commutator components are

| periodic layer | component | exact value |
|---|---|---:|
| `2x2` | `(0,0)` | `-1213745/30233088` |
| `2x3` | `(0,0)` | `-127393295/2176782336` |

Thus the natural physical one-parameter family already fails exactly on both requested small layers.

A stronger local test treats `a,b` as independent coordinates.  At an isotropic interior point `T(x,x)`, a differentiable commuting curve with velocity `(da,db)` must solve `[T,dT]=0`.  On the periodic `2x2` layer, the rows from matrix components `(1,0)` and `(1,15)` give the exact coefficient matrix stored in the artifact.  Its determinant factors as

```text
4 x^6 (x-1)^4 (x+1)^4 (x^2+1)^2 (x^4+3).
```

[COMPUTATION] This is nonzero for every `0<x<1`; hence `da=db=0`.  No nonconstant differentiable commuting curve in the displayed two-weight ansatz passes through an interior isotropic point.  Coincident matrices commute trivially, and exceptional/decoupled boundary branches are not claimed as physical spectral families.

[UNDECIDED] This is not a classification of arbitrary layer-weight parametrizations, singular exceptional branches, or arbitrary `2^(L_xL_y)`-square matrices.  The finite exact result is a no-go for the explicit natural two-weight Ising parametrization, with an all-interior local tangent certificate on `2x2` and direct pair witnesses on both `2x2` and `2x3`.

## Decision

[COMPUTATION] In the six-space edge-rapidity ansatz, the verified Ising parity tensor admits degenerate spectral TE components, but every solution on its nonzero-weight torus has an exceptional leg `x_i in {+1,-1}`; the saturated calculation more generally excludes solutions with all six `x_i` outside `{0,+1,-1}`.  Physical isotropic weights `0<t<1` are excluded exactly.  Partition-preserving leg gauges, the tested direct-sum auxiliary state, and regular finite limits do not repair this; broader interacting auxiliary extensions and singular dimension-changing limits remain explicitly undecided.  The natural physical layer-transfer curve also fails exact commutativity on periodic `2x2` and `2x3`, and has no nonzero differentiable commuting tangent through any interior isotropic point in its two-weight ansatz.

[UNDECIDED] The general four-independent-`R` spectral equation and other rapidity representations remain open; the exact no-go must not be extended beyond the displayed edge-rapidity parametrization.
