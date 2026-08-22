# Exact Dolan–Grady deformation search: a cancellation and bounded no-go results

**Status.** Exact symbolic computation over `Q`, produced by
`experiments/e17_dg_deformation.py` and recorded in
`results/deformation/dg_deformation_search.json`.

The most important outcome is **positive but not an integrable solution**:
within family F3 there is a translation- and point-group-invariant deformation of `A` that cancels
the complete second Dolan–Grady residual on the `3 x 3` torus.  The explicit deformed pair fails
the first Dolan–Grady relation for every floating normalization.  The computation therefore rules
out several natural deformation subfamilies, and answers the quartic-cancellation question, but it
does **not** prove a no-go theorem for all deformations or produce a Dolan–Grady integrable 3D
Ising model.

## 1. Basis, graphs, and symmetry

All expansion uses the repository convention

    Q_(a|b) = X^a Z^b,
    Q_v Q_w = (-1)^(b_v . a_w) Q_(v+w).

A displayed `Y` denotes the canonical phase-free word `XZ`.  In this basis every commutator
coefficient is an integer and every solved linear or polynomial system is over `Q`.  The explicit
positive F3 solution below contains only commuting-site products of `X` and `Z`, so it is also a
real Hermitian operator in the ordinary Pauli convention; it does not depend on the phase-free
`Y` convention.

The searched graphs and imposed spatial groups are:

| layer | sites | bonds | group | order |
|---|---:|---:|---|---:|
| open `2 x 3` grid | 6 | 7 | rectangle point group `D2` | 4 |
| open `3 x 3` grid | 9 | 12 | square point group `D4` | 8 |
| periodic `3 x 3` grid | 9 | 18 | `(Z_3 x Z_3) semidirect D4` | 72 |

The group permutes sites and does not rotate the Pauli labels `X,Y,Z`.  Every parameter multiplies
a full orbit sum.  On the open rectangles there is no nonidentity graph translation; on the torus
all 9 translations are imposed.

## 2. Literal local families and parameter counts

Words are deduplicated before taking symmetry orbits.

* **F1:** every nonidentity Z-only word of weight `1..4` whose support lies in one elementary
  plaquette.  This includes one-site Z fields, edge and diagonal ZZ words, plaquette triples, and
  the four-site plaquette word.
* **F2:** every nonidentity word of weight 1 or 2 supported in one nearest-neighbour edge, with
  independent labels in `{X,Y,Z}`.  Thus it includes all one-site words and all 9 two-site label
  choices such as `XX`, `XY`, ..., `ZZ`.
* **F3:** every nonidentity word of weight `1..4` whose support lies in the closed neighbourhood
  `N[v]` of at least one site.
* **F4:** every phase-free word `Y_v product_(u in T) Z_u`, where `T` is a three-element subset of
  the open neighbourhood of `v`.  These are the quartic words occurring in the undeformed DG
  residual.
* **F_ALL:** the deduplicated union of F1–F4.

The table gives `number of literal words / number of invariant orbit parameters`:

| layer | F1 | F2 | F3 | F4 | F_ALL |
|---|---:|---:|---:|---:|---:|
| `2 x 3` open | `27 / 9` | `81 / 27` | `621 / 186` | `2 / 1` | `623 / 187` |
| `3 x 3` open | `49 / 11` | `135 / 27` | `1800 / 318` | `8 / 2` | `1804 / 319` |
| `3 x 3` torus | `90 / 5` | `189 / 9` | `6102 / 136` | `36 / 1` | `6111 / 137` |

These are literal enumerated counts, not estimates.

Overall rescaling is handled by an explicit projective chart.  For a `B` deformation, the
coefficient of the lexicographically first invariant nearest-neighbour `ZZ` orbit is fixed to its
Ising value.  For an `A` deformation, the analogous first one-site `X` orbit is fixed.  On the
homogeneous torus there is one such orbit, so the charts contain every deformation with nonzero
coefficient along the original `B` or `A`, respectively.  Solutions on the hyperplanes where that
coefficient vanishes are outside these charts unless they are separately classified below.
The floating `lambda` absorbs the corresponding nonzero rescaling of `B`; the DG2 equation is
homogeneous in `A`.

## 3. Controls

For fixed `A,B`, floating `lambda` is a one-column exact linear problem:

    ad_B^3(A) = 16 lambda ad_B(A).

The same exact Pauli machinery gives:

| control | Pauli equations | unknowns | result |
|---|---:|---:|---|
| cycle `C6` | 12 | 1 | `lambda = 1`, residual has 0 terms |
| path `P6` plus bond `(1,4)` | 14 | 1 | no solution; 2 quartic equations are nonzero |

As a stronger pipeline check, the general three-parameter tridiagonal relation from
`experiments/e21_tridiagonal.py` was solved with the same multiplication and rational elimination.
For `C6`, 96 equations in `(beta,gamma,rho)` give `(2,0,16)`.  For the path plus one bond, 84
equations have no solution.  Thus the pipeline finds the known ring relation and detects the first
degree-three obstruction after one bond is added.

## 4. Complete fixed-B search over each A family

### 4.1 Why this nonlinear-looking search is finite and exact

Fix the standard Ising `B` and write

    A' = A + sum_j s_j P_j.

For fixed `lambda`,

    ad_B (ad_B^2 - 16 lambda) A' = 0

is linear in all `s_j`.  Completeness in `lambda` does not require a Gröbner basis.  In the common
Z eigenbasis, the eigenvalues of `B = sum_e Z_i Z_j` are integers and every nonzero matrix element
of `A'` has an energy difference

    Delta = E(z) - E(z') = 2m,       |m| <= number of bonds.

For a noncommuting matrix element the DG2 polynomial requires

    Delta (Delta^2 - 16 lambda) = 0,

so necessarily `lambda = m^2/4`.  The script enumerates every `m=1,...,number of bonds`, and also
solves separately for an affine `A'` commuting with `B` (which would satisfy DG2 for arbitrary
`lambda`).  Each candidate is one exact rational rank computation.  Hence this is a complete
search in the stated fixed-B projective chart, over `Q` and also over `C` as far as consistency and
dimension are concerned.

On the clean `3 x 3` torus the exact polynomial-system sizes and results are:

| family | raw orbits | scale-fixed `s` | unknowns incl. `lambda` | Pauli equations (orbit equations) | result |
|---|---:|---:|---:|---:|---|
| F1 | 5 | 5 | 6 | `72 (2)` | no commuting or noncommuting solution |
| F2 | 9 | 8 | 9 | `1368 (33)` | no commuting or noncommuting solution |
| **F3** | 136 | 135 | 136 | **`53424 (899)`** | **one candidate, `lambda=1`; affine dimension 33** |
| F4 | 1 | 1 | 2 | `144 (6)` | no commuting or noncommuting solution |

There are 18 possible nonzero `lambda` candidates on this graph.  F3 is consistent only at
`lambda=1`: the exact matrix has rank 102 in 135 variables, hence a 33-dimensional affine solution
space.  One sparse particular point is

    A' = sum_v X_v [ 1 - ( Z_(v+x) Z_(v-x) + Z_(v+y) Z_(v-y) ) / 2 ],
    B  = sum_<uv> Z_u Z_v,
    lambda = 1.                                                     (1)

It uses one F3 orbit, represented on the `3 x 3` torus by `X0 Z1 Z2`, of size 18, with coefficient
`-1/2`.  Its support has weight 3 and lies in `N[v]`, exactly as required by F3.

### 4.2 Direct proof of the cancellation

For a Z-basis configuration, put

    S_v = z_(v+x) + z_(v-x) + z_(v+y) + z_(v-y),
    F_v = 1 - [z_(v+x)z_(v-x) + z_(v+y)z_(v-y)]/2.

Flipping `v` changes the eigenvalue of `B` by `Delta_v = +/- 2 S_v`.  There are three cases:

* `|S_v|=4`: both opposite-neighbour products are `+1`, so `F_v=0`; the unwanted
  `|Delta_v|=8` matrix element is removed.
* `|S_v|=2`: one opposite-neighbour product is `+1` and the other is `-1`, so `F_v=1` and
  `Delta_v^2=16`.
* `S_v=0`: `Delta_v=0`, so every odd commutator with `B` vanishes regardless of `F_v`.

Therefore every matrix element of `A'=sum_v X_v F_v` obeys
`Delta_v^3=16 Delta_v`, which proves exactly

    [B,[B,[B,A']]] = 16 [B,A'].

The sparse Pauli expansion independently gives zero residual term by term.  Thus the answer to the
weaker quartic question is **yes**, and in fact (1) cancels the entire DG2 residual, not merely its
quartic part.

### 4.3 Why this is not yet an integrable model

After deforming `A`, DG1 is no longer automatic.  The script floats a second constant `mu` in

    [A',[A',[A',B]]] = 16 mu [A',B].

For the particular point (1), this is an exact system of 1872 Pauli equations (32 symmetry-orbit
equations) in one unknown.  Its coefficient rank is 1 and augmented rank is 2, so no `mu` exists.
An explicit obstruction is the phase-free word

    Y0 Z1 X2 X3:
       coefficient in ad_A'^3(B) = 3,
       coefficient in ad_A'(B)   = 0.

Consequently (1) is **not** a Dolan–Grady pair and does not establish an integrable 3D-Ising-family
model.  The other points in the 33-dimensional DG2 solution space were not subjected to a complete
simultaneous DG1 polynomial solve, so this calculation also does not prove that every point in that
space fails DG1.

### 4.4 Open-layer A-deformation results

The same complete fixed-B calculation was run, not merely counted, on both open grids.  Boundary
site orbits make these finite systems substantially less rigid than the homogeneous torus:

| layer | family | commuting affine space | noncommuting `lambda : affine dimension` |
|---|---|---:|---|
| `2 x 3` | F1 | none | none |
|  | F2 | none | `1 : 13` |
|  | F3 | dimension 22 | `1/4 : 64`, `1 : 65`, `9/4 : 26`, `4 : 29`, and `25/4, 9, 49/4 : 22` each |
|  | F4 | none | none |
| `3 x 3` open | F1 | none | none |
|  | F2 | none | `1 : 8` |
|  | F3 | dimension 23 | `1/4 : 49`, `1 : 73`, `9/4 : 27`, `4 : 28`, and `25/4, 9, 49/4, 16, 81/4, 25, 121/4, 36 : 23` each |
|  | F4 | none | none |

Here “commuting” means an affine `A'` with `[B,A']=0`, so DG2 is vacuous for every `lambda`;
these are classified as degenerate DG2 solutions, not integrable models.  The noncommuting rows are
genuine DG2 solution spaces in the stated scale chart.  No representative from these open-layer
spaces was tested against DG1.  They are therefore recorded as additional finite-boundary DG2
solutions, with **no integrability claim**.  The periodic torus is singled out above because it is
the boundaryless, translation-invariant layer relevant to the bulk model.

## 5. B deformations: exact first-order quartic systems

With `A` fixed, write `B' = B + sum_j t_j Q_j`.  The full residual is cubic in `t`.  Before attempting
large Gröbner systems, the script forms exactly the Newton/first-order quartic equation

    R_4(B) + D R_4(B)[delta B] = 0.                                (2)

The derivative includes all three placements of `delta B` in `ad_B^3(A)` and the derivative of
`-16 ad_B(A)`.  The overall-B direction is removed by the scale chart described above.  There is
one equation per Pauli word; symmetry-related equal rows are compressed only after all word
coefficients have been generated and checked.

The table gives `scale-fixed unknowns`, `full Pauli equations (orbit equations)`, and
`rank / augmented rank`:

| layer | family | unknowns | equations | rank / augmented | (2) solvable? |
|---|---|---:|---:|---:|---|
| `2 x 3` | F1 | 8 | `18 (5)` | `3 / 3` | yes, trivial bond deletion |
|  | F2 | 26 | `286 (78)` | `19 / 19` | yes, same deletion |
|  | F3 | 185 | `1044 (282)` | `88 / 88` | yes, same deletion |
|  | F4 | 1 | `22 (8)` | `1 / 2` | no |
|  | F_ALL | 186 | `1044 (282)` | `89 / 89` | yes, same deletion |
| `3 x 3` open | F1 | 10 | `60 (10)` | `4 / 5` | no |
|  | F2 | 26 | `884 (126)` | `17 / 18` | no |
|  | F3 | 317 | `5113 (711)` | `172 / 173` | no |
|  | F4 | 2 | `80 (13)` | `2 / 3` | no |
|  | F_ALL | 318 | `5113 (711)` | `173 / 174` | no |
| `3 x 3` torus | F1 | 4 | `252 (5)` | `2 / 3` | no |
|  | F2 | 8 | `2304 (40)` | `5 / 6` | no |
|  | F3 | 135 | `9909 (188)` | `86 / 87` | no |
|  | F4 | 1 | `432 (10)` | `1 / 2` | no |
|  | F_ALL | 136 | `9909 (188)` | `87 / 88` | no |

On the torus the undeformed quartic residual is one invariant orbit of 36 words, represented by
`Y0 Z1 Z2 Z3`, with phase-free coefficient `-48`.  In every scale-fixed family its invariant
orbit-sum lies outside the image of the linearized map; the rank increase in the last column is an
exact certificate.

On the `2 x 3` grid the particular solution in F1/F2/F3 changes only the central edge orbit
`Z1 Z4` by `-1`.  It deletes that bond, leaving the perimeter cycle `C6`.  The full nonlinear
weighted-edge solve below confirms that this is a trivial 2-regular solution, not a surviving
2D interaction.

Equation (2) is only the exact affine first-order model.  Its inconsistency is **not** asserted to
be a global theorem about roots of the full cubic system.

## 6. Full nonlinear B subfamilies

### 6.1 All invariant nearest-neighbour ZZ couplings

This reduced family is contained in F1, F2, and F3 and allows edge deletion, anisotropy, and overall
rescaling.  Full DG2 Gröbner systems over `Q` give:

* **`2 x 3`:** 16 Pauli equations, 5 orbit equations, 4 unknowns
  `(c0,c1,c2,lambda)`.  The complete variety consists of:
  * `c2=0`, `c1^2=c0^2`, `lambda=c0^2`: the equal-magnitude perimeter `C6`;
  * `c0=0`, with `c1(c1^2-4lambda)=0` and `c2(c2^2-4lambda)=0`: disconnected
    vertical dimers or `B'=0`.
* **`3 x 3` open:** 32 Pauli equations, 6 orbit equations, 3 unknowns.  The Gröbner basis forces
  `c1=0`; either `lambda=c0^2`, giving the outer `C8` plus an isolated centre, or `B'=0`.
* **`3 x 3` torus:** 72 Pauli equations, 2 orbit equations, 2 unknowns.  The basis contains
  `c0^3` and forces `c0=0`; only `B'=0` remains.

Thus this complete weighted-edge subfamily has no nontrivial connected two-dimensional solution.
Every solution is a path/cycle-type active component, isolated dimers/sites, or the zero generator.

### 6.2 Full torus F1 deformation of B

The five F1 orbit representatives are ordered as

    Z0,
    Z0 Z1                 (nearest-neighbour edges),
    Z0 Z1 Z3,
    Z0 Z1 Z3 Z4           (plaquette),
    Z0 Z4                 (plaquette diagonal).

For `B'=B+sum_(j=0)^4 t_j Q_j`, the full cubic system has 2160 Pauli equations, 45 symmetry-orbit
equations, and 6 unknowns `(t0,...,t4,lambda)`.  A 53-element Gröbner basis over `Q`, followed by
radical case analysis, gives exactly:

1. `t1=-1`, `t2=t3=t4=0`, `lambda=t0^2/4`, so
   `B'=t0 sum_v Z_v`, a set of decoupled one-site fields;
2. its `t0=0` component, `B'=0`, for arbitrary `lambda`.

There is no connected solution in the full invariant F1 B-deformation family.

### 6.3 Full F4 deformation of B

The complete nonlinear F4 systems have unit Gröbner ideal:

| layer | unknowns incl. `lambda` | Pauli equations | orbit equations | Gröbner result |
|---|---:|---:|---:|---|
| `2 x 3` | 2 | 57 | 20 | `[1]`, empty variety |
| `3 x 3` open | 3 | 424 | 63 | `[1]`, empty variety |
| `3 x 3` torus | 2 | 8601 | 156 | `[1]`, empty variety |

Thus adding the residual quartic words directly to `B` cannot satisfy DG2, even before excluding
trivial/decoupled solutions.

### 6.4 Expanded but unsolved full B systems

For reproducibility, the larger cubic systems that were expanded but not globally solved are:

| layer | family | unknowns incl. `lambda` | Pauli equations | orbit equations | status |
|---|---|---:|---:|---:|---|
| `2 x 3` | F1 | 10 | 96 | 28 | expanded, not globally solved |
|  | F2 | 28 | 1213 | 338 | expanded, not globally solved |
| `3 x 3` open | F1 | 12 | 400 | 71 | expanded, not globally solved |
|  | F2 | 28 | 3763 | 556 | expanded, not globally solved |
| `3 x 3` torus | F1 | 6 | 2160 | 45 | solved as above |
|  | F2 | 10 | 8037 | 162 | expanded, not globally solved |

The full B-deformed F3 cubic systems were deliberately not formed: they have respectively 187,
319, and 137 unknowns including `lambda`, and naive cubic expansion over thousands of Pauli words
is not a credible Gröbner calculation.  No conclusion about their global nonlinear varieties is
drawn.  Their exact fixed-B A systems and exact first-order B systems are the ones reported above.

## 7. Precise conclusions

The computation establishes all of the following, with exact arithmetic:

1. The controls behave correctly: `C6` solves DG2 at `lambda=1`; one added bond creates an
   uncancellable quartic equation for the undeformed pair.
2. The quartic defect **can** be cancelled within F3.  Equation (1) cancels the entire DG2
   residual on the `3 x 3` torus and belongs to a 33-dimensional affine DG2 solution space.
3. That particular positive solution fails DG1 for every floating normalization, so it is not a
   Dolan–Grady integrable pair.
4. There is an exact scale-normalized first-order B-deformation obstruction for F1, F2, F3, F4,
   and their union on both the open `3 x 3` grid and the `3 x 3` torus.
5. Full nonlinear no-go results hold for the invariant weighted-edge family, the torus F1
   B-deformation family (apart from decoupled fields), and every F4 B-deformation system listed
   above.

## 8. What is **not** proved

This search does **not** cover:

* simultaneous independent nonlinear deformations of both `A` and `B` over the union F1–F4;
* the complete nonlinear B-deformed F2 or F3 varieties;
* points outside the chosen nonzero-original-generator projective charts;
* translation-breaking or point-group-breaking coefficients;
* Pauli words of weight greater than 4 or support outside the stated local families;
* arbitrary non-Pauli or nonlocal deformations;
* the simultaneous DG1/DG2 variety inside the 33-dimensional F3 DG2 solution space.

Accordingly, this is **not** a no-go theorem for all deformations, all integrability mechanisms, or
an exact solution of the 3D Ising model.  It is an exact classification of the subfamilies stated
above and an explicit demonstration that DG2 quartic cancellation alone is possible but
insufficient.
