# Kac--Ward in 2D, and the first obstruction in 3D

Status of each statement is marked: **[THEOREM]** (external theorem, used as such),
**[COMPUTATION]** (exact machine computation in this repository, script named),
**[NUMERICAL]** (high-precision floating computation, not a proof), **[NOT COVERED]**.

Conventions: `v = tanh K`, `P(v) = sum_{E even} v^|E| = 2^{-N} Z / (cosh K)^{n_b}`,
`Z = 2^N (cosh K)^{n_b} P(v)`.

Scripts: `experiments/e24_kac_ward.py`, `tests/test_kac_ward.py`.
Library: `src/ising/fermions/kac_ward.py`.
Data: `results/kac_ward/kac_ward.json`.

---

## 1. The 2D control identity and its normalisation

Let `G` be a finite planar graph drawn in the plane and `D(G)` its set of directed edges
(`2|E|` of them). For a directed edge `e = (tail(e), head(e))` and a successor `e'` with
`head(e) = tail(e')` and `e' != reverse(e)`, let `theta(e, e') in (-pi, pi)` be the signed
turning angle from the direction of `e` to the direction of `e'`. The Kac--Ward matrix is the
`2|E| x 2|E|` matrix

    Lambda_{e, e'}(v) = v * exp(i theta(e, e') / 2)            when head(e) = tail(e'), e' != reverse(e),
    Lambda_{e, e'}(v) = 0                                       otherwise.

**[THEOREM — Kac & Ward 1952; Whitney 1931]**. For a finite *planar* graph `G` with free
boundary conditions

    det(I - Lambda(v)) = P_G(v)^2 .

Equivalently, with `Z = 2^N (cosh K)^{n_b} P_G(v)`,

    Z^2 = 4^N (cosh K)^{2 n_b} det(I - Lambda).

This is the *squared* normalisation. A common alternative statement is `det(I - Lambda) = Z^2` up
to the prefactor `4^N (cosh K)^{2 n_b}`; we use the polynomial form above because it is the one
that admits an exact integer check.

**[COMPUTATION — `exact_kac_ward_polynomial`]**. We verified the identity *exactly*, coefficient
by coefficient in `Z[v]`, for the free square lattices `3x3`, `3x4`, `4x4`, `4x5`. Method:

* Build `Lambda(v)` over the ring `Z[zeta_8]` with `zeta_8 = exp(i pi/4)`. The turning angle
  on the square lattice is always a multiple of `pi/2`, so `exp(i theta/2) in {1, zeta_8, zeta_8^{-1}}`.
* Compute `det(I - v Lambda)` modulo several primes `p = 1 (mod 8)`, in each of the four
  complex embeddings `zeta_8 -> zeta_8^{1,3,5,7}`. A `4 x 4` inverse Vandermonde recovers the
  `1, zeta_8, zeta_8^2, zeta_8^3` coefficients of every characteristic-polynomial coefficient
  modulo `p`.
* CRT-reconstruct over enough primes that the product modulus exceeds twice a Hadamard bound on
  every cyclotomic-basis coefficient, then centre on `(-M/2, M/2]`.

**The normalisation is fixed empirically.** Every reconstructed coefficient has zero
`zeta_8, zeta_8^2, zeta_8^3` part, and the rational part equals the corresponding coefficient of
`P(v)^2` exactly, with `P(v)` from `ising.exact_enumeration.even_subgraph_polynomial`. The
maximum coefficient discrepancy across all four lattices is **0**. A simultaneous 90-digit
`mpmath` evaluation of `det(I - v Lambda)` agrees with `P(v)^2` to `< 1.5e-90`.

| lattice | `2|E|` | max exact coeff discrepancy | 90-dps determinant error |
|---|---|---|---|
| 3x3 free | 24 | 0 | 2.4e-91 |
| 3x4 free | 34 | 0 | 0 |
| 4x4 free | 48 | 0 | 1.5e-90 |
| 4x5 free | 62 | 0 | 9.8e-91 |

This control is the prerequisite for the 3D question: we have a *certified* 2D identity, against
which any 3D construction must be measured.

---

## 2. Why 2D works: Whitney's theorem and the loop expansion

**[THEOREM — Whitney 1931; the half-angle loop argument]**. The identity is a consequence of two
facts:

**(a) Loop expansion of the determinant.** `det(I - v Lambda) = exp(-sum_{k>=1} (v^k/k) Tr(Lambda^k))`,
and `Tr(Lambda^k)` is the sum over closed non-backtracking walks of length `k` of the product of
the half-angle weights `exp(i theta/2)` along the walk. Equivalently, the expansion of the
determinant in principal minors enumerates collections of edge-disjoint closed non-backtracking
walks ("polygon decompositions"), each weighted by the product of its half-angle phases.

**(b) Whitney's turning-number theorem.** For a closed immersed planar curve with `c`
self-crossings, the total signed turning angle is `2 pi (1 + (c mod 2))` modulo `4 pi`; in the
half-angle convention the phase accumulated around a single closed non-backtracking walk is

    prod exp(i theta/2) = (-1) * (-1)^c .

(For a single simple closed curve `c = 0` and the phase is `-1` — the familiar `-1` per loop.
Each added self-crossing contributes an extra factor `-1`.)

The decisive point is that in *2D* the signs from multiple loops combine correctly:

* a collection of `L` edge-disjoint closed curves, with `C` total crossings between them, carries
  total phase `(-1)^{L + C}`,
* whereas the same even subgraph, viewed as a union of `L` vertex-disjoint cycles, contributes
  `+1` to `P(v)`,
* and the standard Whitney sign reconciliation `(-1)^{L+C} -> +1` per polygon decomposition is
  *exactly* what converts the determinant's loop sum into `P(v)^2`.

The squaring `P^2` is itself a Whitney phenomenon: `P(v)^2` counts ordered pairs of even
subgraphs, and the two "slots" correspond to the two spinor sectors of the planar Kac--Ward
determinant (equivalently, to the two chiralities of the polygon decomposition).

**Why this is planar-specific.** The argument leans on (i) a *signed* turning angle `theta`,
which requires a global orientation of the ambient plane, and (ii) the Whitney turning-number
theorem, which is a statement about curves in the plane. In 3D there is no turning number: the
sign of `theta` between two directions is not defined, and the self-crossing count is no longer
related to any topological invariant of the curve. Section 3 asks whether *any* weighting can
nevertheless reproduce `P(v)^2`.

---

## 3. The 3D question, posed sharply

**Question.** *Does there exist an assignment of complex weights `U_{e,e'}` to consecutive
directed-edge pairs of the simple cubic lattice such that, for every finite free-boundary box,
`det(I - v U) = P(v)^2` coefficient by coefficient in `v`?*

We attack this in four steps: the trace equations (3a), the symmetry reduction (3b), the exact
no-go on the cubic-covariant slice (3c), the independent numerical reproduction (3d), and the
scope of the full translation-invariant ansatz (3e).

### 3a. The trace equations

**[THEOREM — MacMahon / Newton]**. With `U` any directed-edge weight matrix,

    log det(I - v U) = - sum_{k>=1} (v^k / k) Tr(U^k) ,

and `Tr(U^k) = sum_{gamma in CNBW(k)} prod_{t=1}^{k} U_{e_t, e_{t+1}}` is a sum over closed
non-backtracking walks of length `k`. Matching `log det(I - v U) = 2 log P(v)` coefficient by
coefficient therefore requires

    Tr(U^k) = -2 k [v^k] log P(v)                      (***)

for every `k`. This is an infinite sequence of polynomial equations in the unknowns `U_{e,e'}`.
On the simple cubic lattice the first nontrivial orders are `k = 4, 6, 8, ...` (all odd orders
vanish on the bipartite lattice, and `k = 2` vanishes because the only length-2 closed
non-backtracking walk would require an immediate reversal).

### 3b. Symmetry reduction

**Translation invariance.** We restrict to weights shared by every vertex, so that
`U_{e,e'} = U(d(e), d(e'))` depends only on the ordered pair of lattice directions
`(d, d')`. The non-backtracking condition `d' != -d` leaves `6 * 5 = 30` unknowns
`U(d, d')`.

**Gauge freedom.** The transformation `U(d, d') -> g(d)^{-1} U(d, d') g(d')` for nonzero
`g(d)` leaves `Tr(U^k)` — hence every equation (***) — invariant. We quotient this gauge by fixing
the five entries on a directed spanning tree of the direction-state graph to `1`:

    U(+x,+y) = U(+x,-y) = U(+x,+z) = U(+x,-z) = U(+y,-x) = 1.

This is a valid gauge slice **on the chart where these five entries are nonzero**; it leaves
`30 - 5 = 25` unknowns. **[NOT COVERED]** The zero branches of this gauge — points of the variety
where a fixed tree weight vanishes — are *not* covered by our computation. A complete
classification would have to cover those branches separately (e.g. by a different spanning tree).

**Strict scalar cubic covariance (the slice we can decide).** If we further impose that the
weights respect the *proper* cubic rotation group `O` (order 24, the orientation-preserving
octahedral group; the full group `O_h` of order 48 adds reflections, which are not needed here)
— i.e. that the weight depends only on the geometric relation between `d` and `d'`, not on the
particular directions — the 30 unknowns collapse to two: `a = U(d, d)` for a straight step and
`b = U(d, d')` for an orthogonal turn. `O` is transitive on the 6 straight ordered pairs `(d, d)`
and on the 24 ordered perpendicular pairs `(d, d')` with `d' perpendicular to `d`: for any
ordered perpendicular pair `(a, b)` the triple `(a, b, a x b)` is right-handed, hence the image
of `(+x, +y, +z)` under a unique proper rotation. Two orbits give two scalar weights, so this is
a genuine 2-parameter slice (not merely a gauge choice). It is the natural 3D analogue of the 2D
Kac--Ward weight, where the single turn weight `exp(i theta/2)` is fixed by planar isotropy.

### 3c. Exact no-go on the cubic-covariant slice

**[COMPUTATION — `bulk_cubic_trace_monomials`, `strict_cubic_equations`]**. We enumerate the
closed non-backtracking walks of length `k = 4, 6, 8` on the infinite cubic lattice, classified by
the number of straight transitions, and read off the exact trace as a polynomial in `a, b`. The
required right-hand side comes from the independently established simple-cubic high-temperature
series (`results/series/sc_ht_free_energy.json`):

| `k` | walks/site (straight, turn) | `Tr(U^k)` (per site) | required `-2k [v^k]log P` |
|---|---|---|---|
| 4 | 24 (all turn) | `24 b^4` | `-24` |
| 6 | 192 turn + 72 straight-turn | `72 a^2 b^4 + 192 b^6` | `-264` |
| 8 | 1080 turn, 2304 one-straight, 1056 two-straight, 144 four-straight | `144 a^4 b^4 + 1056 a^2 b^6 + 2304 a b^7 + 1080 b^8` | `-3000` |

(The `k = 8` trace contains an `a b^7` term because a length-8 closed walk can have exactly one
straight transition followed by seven turns; this odd monomial is what makes the slice
inconsistent.)

Setting each `Tr(U^k) = required` gives three polynomial equations. **`groebner([f4, f6], a, b) =
{2 polynomials}`** — the first two equations are consistent and leave a 0-dimensional variety of
candidates. **`groebner([f4, f6, f8], a, b) = {1}`** — adding the `k = 8` equation makes the
ideal the whole ring. **There is no solution `(a, b)` to all three equations.**

This is an exact, rational-arithmetic proof that **no strict cubic-covariant scalar turn rule
reproduces `P(v)^2` past order `v^6`**; the first obstruction is at order `v^8`.

### 3d. Independent numerical reproduction

**[COMPUTATION — `fit_strict_cubic_candidate`, `numeric_determinant_prefix`]**. As an independent
check that does not rely on the Groebner computation, we fit `(a, b)` from the `k = 4, 6`
equations on the *finite* `3x3x2` free box (whose `P(v)` is exactly known from the repo), then
ask where the resulting `det(I - v U)` first departs from `P(v)^2`.

On `3x3x2` free, `P(v) = 1 + 20 v^4 + 78 v^6 + 402 v^8 + ...`, so `P(v)^2 = 1 + 40 v^4 + 156 v^6
+ 1204 v^8 + ...`. The `k = 4, 6` fit gives `b = zeta_8 = (1+i)/sqrt2`, `a^2 = (39 - 32 i)/7`.
With these weights the determinant prefix is `1 + 40 v^4 + 156 v^6 + (8676/7 + ...)v^8`, which
**disagrees with the target `1204` at order `v^8`** — exactly the order predicted by the bulk
Groebner computation. A dense 90-digit `mpmath` matrix-power computation of `det(I - v U)` agrees:
first (and only, through `v^8`) failure at `v^8`, coefficient `~798.6 + 1575.0 i` against target
`1204`.

The finite-box first-failure order (`v^8`) coincides with the bulk first-failure order (`v^8`),
which is the expected consistency check: the `k = 8` trace is the first one to see the
odd-straight monomial `a b^7`, and that monomial has no counterpart in `log P`.

### 3e. The natural 3D guess, falsified

**[COMPUTATION — `natural_three_dimensional_guess` in `e24_kac_ward.py`]**. The most natural 3D
analogue of the 2D Kac--Ward weight assigns the straight step weight `1` and every orthogonal turn
the unsigned half-angle `exp(i pi/4)` (unsigned because 3D supplies no left/right turning sign).
On the `2x2x2` free box this guess reproduces `det(I - v U)` through `v^4` (`12 v^4`, matching
`P^2`) but **fails already at `v^6`**: the determinant coefficient is `32 i` while `P(v)^2`
requires `32`. The guess is purely imaginary where the target is real — a direct, numerical,
Groebner-independent falsification. The earlier failure order (6 vs 8) reflects that this
particular guess does not even satisfy the `k = 6` trace equation, whereas the fitted candidate of
3d does.

---

## 4. Quantifying the obstruction: walks versus even subgraphs

**[COMPUTATION]**. The following counts make the obstruction quantitative. In 2D the Kac--Ward
identity *forces* a precise relationship between the weighted trace sum and the even-subgraph
count; in 3D no such relationship holds for any strict cubic-covariant scalar turn
rule.  (Scope: the full 30-weight scalar family is NOT disproved here -- see sec. 5.  The
quantifier was narrowed after an adversarial audit finding.)

### 4a. 2D (the identity forces the relation)

| lattice | `k` | closed NBW | even subgraphs of size `k` | required weighted trace `-2k [v^k] log P` |
|---|---|---|---|---|
| 3x3 free | 4 | 32 | 4 | `-32` |
| 3x3 free | 6 | 48 | 4 | `-48` |
| 3x3 free | 8 | 240 | 7 | `+16` |
| 4x4 free | 4 | 72 | 9 | `-72` |
| 4x4 free | 6 | 144 | 12 | `-144` |
| 4x4 free | 8 | 936 | 50 | `-152` |
| 4x5 free | 4 | 96 | 12 | `-96` |
| 4x5 free | 6 | 204 | 17 | `-204` |
| 4x5 free | 8 | 1392 | 89 | `-272` |

At each order, the Kac--Ward weighted trace `Tr(Lambda^k)` (in the `1, zeta_8, zeta_8^2,
zeta_8^3` basis) has *exactly* the required rational value `-2k [v^k] log P` and zero
non-rational components — this is the content of the certified identity, verified per order in
`kac_ward_trace_basis`. The unweighted NBW count and the even-subgraph count are *not* directly
equal (e.g. 32 vs 4 at `k = 4`); the half-angle phases perform a nontrivial signed cancellation
that reduces the former to the latter.

### 4b. 3D (the relation cannot hold for a strict cubic-covariant scalar turn rule)

On the `3x3x2` free cubic box:

| `k` | closed NBW | even subgraphs of size `k` | required weighted trace `-2k [v^k] log P` |
|---|---|---|---|
| 4 | 160 | 20 | `-160` |
| 6 | 936 | 78 | `-936` |
| 8 | 6688 | 402 | `-3232` |

The required weighted trace is a *specific* integer at each order. For the scalar cubic-covariant
rule the trace is `sum_s N_{k,s} a^s b^{k-s}` where `N_{k,s}` is the number of length-`k` closed
NBWs with `s` straight transitions. At `k = 4` (all-turn, `N_{4,0} = 160`) the equation
`160 b^4 = -160` is solvable (`b^4 = -1`). At `k = 6` the two-monomial equation is still
solvable. **At `k = 8` the four-monomial trace `144 a^4 b^4 + 1056 a^2 b^6 + 2304 a b^7 + 1080 b^8`
cannot equal the integer `-3232`** for any `(a, b)` that solved the lower orders — this is the
content of the `groebner = {1}` verdict. The analogous relation in 3D therefore *fails at order
`v^8`*, the first order at which the closed-walk structure of the cubic lattice is rich enough to
produce a monomial (`a b^7`) incompatible with the even-subgraph count.

---

## 5. Scope: what is proved, what is computed, what is not covered

### 5a. Theorem (external, used as such)

* The 2D Kac--Ward identity `det(I - Lambda) = P^2` for finite planar graphs (Kac & Ward 1952;
  Whitney 1931), together with the Whitney turning-number theorem that underlies it.

### 5b. Our computation (exact, in this repository)

* Certified exact verification of `det(I - Lambda) = P^2` on `3x3, 3x4, 4x4, 4x5` free square
  lattices, by `Z[zeta_8]` modular reconstruction with a Hadamard certificate
  (`exact_kac_ward_polynomial`); max coefficient discrepancy 0.
* Exact enumeration of closed non-backtracking cubic walks by straight-transition count
  (`bulk_cubic_trace_monomials`, `trace_type_counts`).
* Exact Groebner proof that the strict cubic-covariant 2-parameter slice `(a, b)` is
  inconsistent at order `k = 8`: `groebner([f4, f6, f8], a, b) = {1}`.
* Independent dense-matrix numerical reproduction on `3x3x2` free: first determinant-coefficient
  failure at `v^8` (`numeric_determinant_prefix`).
* Numerical falsification of the natural unsigned-half-angle guess on `2x2x2` free: first failure
  at `v^6`, determinant coefficient `32 i` vs target `32`.

### 5c. The full translation-invariant 30-weight family — explicitly unresolved

We also built the exact finite-box trace system for the full translation-invariant,
gauge-quotiented 30-weight ansatz (25 unknowns after fixing the five-entry spanning-tree gauge on
its nonzero chart), across ten free boxes `{2,3,4} x {2,3,4} x {2,3}` at orders `k = 4, 6, 8`
(`full_weight_finite_system`, 30 equations, 25 unknowns).

* Through order `k = 6`: the Groebner basis has 13 elements and does *not* contain `1` — the
  system is **consistent**, no obstruction.
* Through order `k = 8`: the Groebner computation did not terminate within the allotted budget
  (20 s); the system is **unresolved**. We do not claim a no-go theorem for this family.

The scalar-covariant obstruction of 5b is therefore a *proper subspace* result: it rules out the
2-parameter isotropic slice but does not rule out the full 25-parameter (or 30-parameter,
pre-gauge) family. The two results are reported separately and not conflated.

### 5d. What is NOT covered

The following families of constructions are **not** addressed by anything in this note and must
not be read as excluded:

1. **The full 30-weight translation-invariant family** at order `k >= 8` (unresolved, 5c).
2. **Zero-entry gauge branches** — points where a spanning-tree weight fixed to `1` in our gauge
   slice actually vanishes; a different tree would be needed there.
3. **Projective / spin-structure covariance** — weights that transform projectively under the
   cubic group, or that depend on a choice of spin structure (Kasteleyn orientation analogue).
4. **Edge- or vertex-dependent weights** — dropping translation invariance, so that `U` depends
   on the actual edges, not just their directions.
5. **Sums over spin structures / Kasteleyn orientations** — a determinant that is itself a sum of
   determinants over several signed orientations.
6. **Larger auxiliary space / higher bond dimension** — replacing each scalar weight `U(d, d')`
   by an `m x m` matrix (bond dimension `m >= 2`). The obstruction here is proved only for `m = 1`.

### 5e. Bond dimension 2 — status

The acceptance criterion asks, if possible, for a test of bond dimension 2. We did not implement
it. With `m = 2`, each of the 30 direction-pair weights becomes a `2 x 2` complex matrix (8 real
parameters per pair, 240 real parameters pre-gauge), the trace equations become matrix-polynomial
equations in these parameters, and even the `k = 4` equation is a system of 4 matrix equations in
240 unknowns with no obvious symmetry reduction. The Groebner approach that decides the scalar
slice is not feasible at this size, and we have no independent certificate. **Bond dimension 2 is
therefore an open problem; it is recorded as `bond_dimension_2_tested: false` in the result JSON,
and no claim is made about it.**

---

## 6. Reproduction

```text
.venv/bin/python tests/test_kac_ward.py     # prints PASS
.venv/bin/python experiments/e24_kac_ward.py  # writes results/kac_ward/kac_ward.json, prints PASS
```

The standalone test checks the 2D exact coefficient identity on four free square lattices, the 2D
90-digit determinant agreement, the strict-cubic `groebner([f4,f6,f8])={1}` verdict, and the
independent finite-matrix first-failure at `v^8`, then prints `PASS`.

## 7. Dependency ledger

| Claim | Logical dependency | Status |
|---|---|---|
| `det(I - Lambda) = P^2` (2D) | Kac--Ward 1952 / Whitney 1931 | external theorem |
| exact coefficient reconstruction | `Z[zeta_8]` modular charpoly + CRT + Hadamard bound | exact, this repo |
| simple-cubic `[v^k] log P` | `results/series/sc_ht_free_energy.json` (FLM series, independently exact) | exact, this repo |
| closed NBW counts on cubic lattice | exhaustive backtracking enumeration | exact, this repo |
| strict-cubic `groebner = {1}` at `k = 8` | sympy `groebner` over `Q[a, b]` | exact, this repo |
| finite-box `v^8` first failure | dense `mpmath` matrix powers, 90 dps | numerical corroboration |
| full 30-weight family, `k <= 6` | sympy `groebner`, 25 unknowns | exact, consistent |
| full 30-weight family, `k = 8` | sympy `groebner`, 25 unknowns | **unresolved** (timeout) |
| natural guess `v^6` failure | dense `mpmath` matrix powers | numerical falsification |

---

## 8. Flat `Z2` twists of the phase-free Hashimoto determinant

**[DEFINITION — graph 1-complex]** For an undirected edge signing
`s:E(G)->{+1,-1}`, let the phase-free twisted Hashimoto matrix be

    B_s[(u,v),(v,w)] = s_{vw}     if w != u,
    B_s[(u,v),(x,w)] = 0          if v != x or w = u.

Vertex switching is gauge.  For a connected graph, fixing all signs on a
spanning tree to `+1` leaves one sign on each chord and therefore exactly
`2^b1` representatives, where `b1=|E|-|V|+1`.  Here “flat” refers to the graph
as a 1-complex: it has no 2-cells, so every `Z2` edge cochain is a cocycle.

**[LEMMA — all finite connected graphs]** Uniform averaging over
`H^1(G;Z2)` projects the determinant circuit-collection expansion onto total
mod-2 edge class zero.  Indeed, the signing product of a collection with
class `h in H_1(G;Z2)` is the character value `chi(h)`, and

    2^(-b1) sum_{chi in H^1(G;Z2)} chi(h) = 1 if h=0, and 0 otherwise.

This statement concerns the total class of a determinant collection; it does
not assert that each individual circuit in a surviving collection is
homology-trivial.

**[THEOREM — exact acyclic/cyclic boundary].** Let \(G\) be a finite
connected simple graph, and average the phase-free determinant uniformly over
all graph-\(H^1(G;\mathbb Z_2)\) twists.  This average equals \(P_G(v)^2\) for
a forest and fails for every cyclic \(G\).  More precisely, if \(g\) is the
girth and \(c_g>0\) the number of simple \(g\)-cycles, then the first mismatch is

```text
[v^g] uniform twist average = 0,
[v^g] P_G(v)^2             = 2 c_g.
```

Every shortest closed nonbacktracking walk is an oriented simple \(g\)-cycle,
whose nonempty edge set is a nonzero class in the homology of the graph
one-complex.  Character orthogonality kills it.  A determinant collection with
two positive-length circuits has degree at least \(2g\), so nothing else
contributes at degree \(g\).  On the Ising side, the smallest nonempty even edge
sets are exactly the simple \(g\)-cycles, hence
\([v^g]P_G=c_g\) and \([v^g]P_G^2=2c_g\).  For a forest the nonbacktracking
matrix is nilpotent and \(P_G=1\).

This theorem is only about the **uniform phase-free average over every graph
one-complex character**.  It does not exclude a fixed twist, nonuniform
combination, a surface spin-structure formula with 2-cells, or another local
turning-phase rule.

### 8a. Exact 2D positive and negative controls

**[COMPUTATION — finite path control, `e239_twist_hashimoto.py`]** On the
three-vertex open path, `P(v)=1`, the phase-free Hashimoto matrix is nilpotent,
and both its determinant and the planar Kac--Ward determinant are exactly
`1`.  This is a positive construction control, but an acyclic graph cannot
distinguish local turning phases because it has no closed nonbacktracking
circuit.

**[COMPUTATION — finite `C4` controls, `e239_twist_hashimoto.py`]** For the
embedded square cycle,

    P(v)^2 = (1+v^4)^2 = 1+2v^4+v^8.

The exact directed-edge determinants are:

| local/twist data | exact determinant |
|---|---|
| phase-free, trivial holonomy | `(1-v^4)^2 = 1-2v^4+v^8` |
| phase-free, fixed nontrivial holonomy | `(1+v^4)^2` |
| planar signed half-angle phases, trivial holonomy | `(1+v^4)^2` |
| uniform average of the two phase-free twists | `1+v^8` |
| uniform average of the two phased twists | `1+v^8` |

**[COMPUTATION — finite control conclusion]** The planar signed half-angle
turns supply loop phase `-1`.  A fixed nontrivial `Z2` holonomy can mimic that
single loop sign on `C4`, but uniform twist averaging deletes the `v^4` term
instead of recreating it.  Thus the path checks implementation and the square
identifies the indispensable signed-turn/loop phase data.

### 8b. Complete open-cube twist calculation

**[COMPUTATION — finite open `2x2x2` cube only]** The cube graph has
`n=8`, `m=12`, and `b1=5`.  In the root-0 breadth-first spanning-tree gauge,
`e239_twist_hashimoto.py` enumerates all 32 chord-sign classes and computes
every polynomial exactly with the signed Bass--Ihara identity

    det(I-v B_s) = (1-v^2)^(m-n) det(I-v A_s+v^2(D-I)).

There are six distinct determinant polynomials, with multiplicities
`1,1,3,3,12,12`; the artifact stores all 32 representatives and all 25
coefficients of each polynomial.

**[COMPUTATION — finite even-subgraph target]** Direct enumeration of the
32 even edge sets gives

    P(v) = 1 + 6v^4 + 16v^6 + 9v^8,

and hence

    P(v)^2 = 1 + 12v^4 + 32v^6 + 54v^8 + 192v^10
             + 364v^12 + 288v^14 + 81v^16.

**[COMPUTATION — finite negative certificate]** The exact uniform twist
average is

    2^(-5) sum_s det(I-v B_s)
      = 1 + 6v^8 + 16v^12 - 48v^14 - 87v^16 - 48v^18
        + 672v^20 - 768v^22 + 256v^24.

The first coefficient mismatch is therefore already

    [v^4] average = 0,        [v^4] P(v)^2 = 12.

This is an honest finite negative result on the open cube, not an all-size
no-go theorem.

### 8c. Homology-trivial circuit certificate at the first obstruction

**[COMPUTATION — finite cube, lengths 4, 6, and 8]** Exact rooted
closed-nonbacktracking-walk enumeration gives

| length | all rooted closed NB walks | homology-trivial rooted walks |
|---|---:|---:|
| 4 | 48 | 0 |
| 6 | 192 | 0 |
| 8 | 336 | 48 |

**[COMPUTATION — finite length-4 classification]** The 48 length-4 walks are
exactly six faces times two orientations times four roots.  Every face has a
nonzero mod-2 edge class, so character averaging retains none.  In contrast,
the target coefficient 12 consists of the ordered pairs `(empty,face)` and
`(face,empty)` for the six faces.  This explains the degree-4 mismatch rather
than merely locating it.

**[COMPUTATION — finite first-survivor classification]** All 48
homology-trivial rooted length-8 walks go twice around one face; there is no
second type.  In the determinant circuit-collection expansion, the first
nonconstant survivors are the six products of the two oppositely oriented
primitive factors of the same face.  Each contributes `+v^8`, explaining the
average coefficient `+6` at degree 8.

**[COMPUTATION — finite fixed-twist diagnostic]** There is a unique gauge
class whose six face holonomies are all `-1` (class `11110` in the stored
chord order).  Its determinant matches the target coefficient `12` at
degree 4 but first fails at degree 6, where it has `-16` instead of `32`.
Thus assigning the right sign to every square face is not, by itself, the
missing three-dimensional turning rule.

### 8d. A coordinate-phased ansatz, separately scoped

**[DEFINITION — frame-dependent ansatz]** Fix the standard oriented
coordinate frame and `n=(1,1,1)`.  Give a straight transition phase `1`; for
an orthogonal coordinate turn `d->d'`, give phase
`zeta_8^sigma`, where `sigma=sign((d cross d') dot n)`.

**[COMPUTATION — finite open cube through degree 8]** Exact
`Z[zeta_8]` traces and Newton identities give the untwisted determinant prefix

    1 + 12v^4 + 0v^6 + 30v^8,

against target prefix `1+12v^4+32v^6+54v^8`.  Its first mismatch is degree 6.

**[NOT COVERED — coordinate ansatz scope]** The rule selects a frame and
diagonal, is not invariant under the full cubic group, and was tested only as
an ansatz on this one finite graph.  It is not asserted to be a universal 3D
Kac--Ward formula.

### 8e. Artifact, digests, verifier, and scope

**[COMPUTATION — exact artifact]** The complete certificate is
`results/kac_ward/twist_hashimoto.json`, produced by
`experiments/e239_twist_hashimoto.py`.  Its canonical content digests are:

| content | SHA-256 |
|---|---|
| 2D controls | `ef744a87f47c5b08d67c40cb456ac1bbfe5262e884c4a4e6b0cd77ec177f2ecf` |
| all 32 twist determinants | `85da746f0c38d1df412f593f8ec3d5a64a533a23ece2432584048fb14e60b1c5` |
| distinct-polynomial inventory | `02bd32c7e834395b80a928266e128efb2903b83630f0df1b753a7f983f4c5d21` |
| cube comparison | `b6d860dc2adae612620ebbcb81f7b1c68b3e2da0122d800044b5f33fa902975d` |
| homology circuit certificate | `8fe3021ac762552b050636f6db0358314a8bf1591548bc9a8156e8b5304d099b` |
| coordinate ansatz | `451ce532390c02ba58a586cd9fb9050007f42cf14d87f77ecf79a7dd0520912c` |

**[REPRODUCTION — independent exact verifier]**
`tests/test_twist_hashimoto.py` does not import the producer.  It uses a
different spanning tree, builds every `24x24` Hashimoto matrix directly, and
recovers its determinant from exact traces and Newton identities before
comparing the complete polynomial multiset and content digests.

**[SCOPE — explicit nonclaims]** No larger graph or periodic box was launched.
The finite determinant representation is not claimed to provide a tractable
thermodynamic limit, an infinite-volume free energy, a critical point, or a
critical exponent.  The benchmark value of `K_c` played no role.

## 9. Orientable-genus bound for open simple-cubic boxes

Write

\[
G(a,b,c)=P_a\mathbin{\square}P_b\mathbin{\square}P_c
\]

for positive integers \(a,b,c\), with open boundary conditions.  Thus the
vertices are triples
\((x,y,z)\in\{0,\ldots,a-1\}\times\{0,\ldots,b-1\}
\times\{0,\ldots,c-1\}\), and an edge changes exactly one coordinate by one.
Let \(\gamma(G)\) denote the minimum genus of a closed orientable surface in
which \(G\) embeds.

**[THEOREM — all open boxes].** If \(\min(a,b,c)\geq 2\), then
\(G(a,b,c)\) is finite, simple, connected, bipartite, and bridgeless, and

\[
\boxed{\quad
\gamma(G(a,b,c))\ \geq\
\max\!\left(0,\left\lceil
\frac{abc-ab-bc-ca+4}{4}
\right\rceil\right).
\quad}
\]

For arbitrary positive side lengths, the exact planarity classification is

\[
G(a,b,c)\text{ is planar}
\quad\Longleftrightarrow\quad
\min(a,b,c)=1
\ \text{ or, up to permutation, }\ (a,b,c)=(2,2,L).
\]

This is a graph-topological theorem.  In particular, the positive values on
the right are lower bounds; they are not asserted to be exact genera.

### 9a. Counts and elementary graph properties

There are plainly

\[
V=abc.
\]

The edges in the three positive coordinate directions form disjoint classes
of sizes

\[
(a-1)bc,\qquad a(b-1)c,\qquad ab(c-1).
\]

Consequently

\[
E=(a-1)bc+a(b-1)c+ab(c-1)
  =3abc-(ab+bc+ca).
\]

Every edge has distinct endpoints, and its two endpoints determine its
coordinate direction, so the graph is simple.  Moving one coordinate at a
time joins any two vertices, so it is connected.  The parity of \(x+y+z\)
gives a bipartition because every edge flips that parity.

Now suppose \(\min(a,b,c)\geq2\), and take an edge in one coordinate
direction.  Choose either transverse coordinate.  Since its side length is
at least two, the current transverse coordinate has an in-range neighbour:
move by \(+1\) when possible and otherwise by \(-1\).  Shifting both endpoints
of the original edge by that move produces a unit square containing the
edge.  The other three sides of this plaquette still join the endpoints when
the original edge is deleted.  Thus every edge lies on a unit plaquette and
no edge is a bridge.

### 9b. Why a minimum embedding may be taken cellular

**[LEMMA — self-contained regular-neighbourhood argument].** Every finite
connected graph has a cellular embedding whose orientable genus equals its
minimum embedding genus.

Take a minimum-genus embedding of a connected graph \(G\) in a closed
orientable surface \(S_g\), and take a sufficiently small connected regular
neighbourhood \(N\) of the embedded graph.  Concretely, \(N\) is made from a
small disk at every vertex and a narrow rectangular band along every edge.
Its boundary circles are the boundary walks of the resulting ribbon graph.
Cap every boundary circle with a disk.  Shrinking the vertex disks and edge
bands back to their cores carries the regions adjacent to the caps to open
disks.  Hence the capped surface \(S_h\) contains the same embedded graph and
every complementary component is a disk: the new embedding is cellular.

It remains to check that capping did not hide a genus change.  Suppose \(N\)
has genus \(h\) and \(q\) boundary components.  Let the closures of the
\(r\) components of \(S_g\setminus N\) have genera \(k_i\) and \(q_i\)
boundary components.  Every such component has a boundary and
\(\sum_iq_i=q\), so \(q\geq r\).  Euler-characteristic additivity along the
boundary circles gives

\[
\begin{aligned}
2-2g
 &= (2-2h-q)+\sum_{i=1}^{r}(2-2k_i-q_i),\\
g
 &= h+\sum_{i=1}^{r}k_i+q-r
 \ \geq\ h.
\end{aligned}
\]

Capping \(N\) produces a closed surface of genus \(h\) in which \(G\)
embeds, while \(g\) was minimal; hence \(g\leq h\).  Therefore \(g=h\), and
the cellular embedding has minimum genus.  No external cellularity theorem
is needed.

### 9c. Face counting and the genus bound

Apply the lemma to \(G(a,b,c)\) with all sides at least two.  Every facial
boundary in a cellular minimum embedding has length at least four.  Indeed,
a length-one boundary would be a loop.  A length-two boundary using two
different edges would be a parallel-edge digon; if it uses the same edge in
both directions, the local cyclic order forces the edge to be traversed out
and immediately back on both ends, which is the bridge case.  Simplicity and
bridgelessness exclude both possibilities.  A length-three facial boundary
would be an odd closed walk, excluded by bipartiteness.

Each edge has two sides in the cellular embedding, counted with
multiplicity, so

\[
2E=\sum_f |\partial f|\geq4F,
\qquad\text{and therefore}\qquad F\leq\frac E2.
\]

Euler's formula on the genus-\(\gamma\) surface now gives

\[
2-2\gamma=V-E+F
\leq V-\frac E2.
\]

Equivalently,

\[
\gamma\geq\frac{E-2V+4}{4}
=\frac{abc-ab-bc-ca+4}{4}.
\]

Since genus is a nonnegative integer, taking the ceiling and then the maximum
with zero proves the displayed bound.

### 9d. Exact planarity classification

If one side equals one, the graph is a rectangular path grid (with paths and
the single vertex included as boundary cases), so its coordinate-grid drawing
is planar.

For \(2\times2\times L\), regard every \(2\times2\) layer as a 4-cycle.
Draw these \(L\) cycles as successively nested squares.  Join corresponding
corners of consecutive squares by the four radial segments in the annulus
between them.  Distinct annuli and distinct corner rays are disjoint, so this
is a crossing-free drawing for every \(L\).

Conversely, order the dimensions as \(2\leq a\leq b\leq c\).  Put

\[
N(a,b,c)=abc-ab-bc-ca+4.
\]

If \(a=2\), then

\[
N(2,b,c)=(b-2)(c-2),
\]

which vanishes when \(b=2\) and is positive when \(b,c\geq3\).  If
\(a\geq3\), the numerator is coordinatewise nondecreasing on side lengths at
least two, since

\[
N(a+1,b,c)-N(a,b,c)=bc-b-c=(b-1)(c-1)-1\geq0,
\]

and cyclically.  In this case
\(N(a,b,c)\geq N(3,3,3)=4\).  Thus every box not covered by the two planar
constructions has \(\gamma\geq1\), completing the classification.

### 9e. Cube and elongated-family corollaries

For open cubes,

\[
N(L,L,L)=L^3-3L^2+4=(L-2)^2(L+1),
\]

so, for \(L\geq2\),

\[
\gamma(G(L,L,L))
\geq\left\lceil\frac{(L-2)^2(L+1)}4\right\rceil.
\]

For \(L=2,\ldots,8\), these lower bounds are respectively

\[
0,\ 1,\ 5,\ 14,\ 28,\ 50,\ 81.
\]

For a fixed \(a\times b\) cross-section and \(L\geq2\),

\[
\gamma(G(a,b,L))
\geq
\max\!\left(0,\left\lceil
\frac{(ab-a-b)L-ab+4}{4}
\right\rceil\right).
\]

The numerator has positive slope \(ab-a-b\) for every \(a,b\geq2\) except
\((a,b)=(2,2)\).  Hence the lower bound, and therefore the genus, is unbounded
along every other fixed open cross-section.  Useful exact specialisations are

\[
\begin{aligned}
\gamma(G(2,b,L))
 &\geq \left\lceil\frac{(b-2)(L-2)}4\right\rceil,\\
\gamma(G(2,3,L))
 &\geq \left\lceil\frac{L-2}{4}\right\rceil,\\
\gamma(G(2,4,L))
 &\geq \left\lceil\frac{L-2}{2}\right\rceil,\\
\gamma(G(3,3,L))
 &\geq \left\lceil\frac{3L-5}{4}\right\rceil.
\end{aligned}
\]

These are lower-bound families only.  No nonplanar row here is labelled with
an exact genus.

### 9f. Exact artifact and independent verifier

**[COMPUTATION — exact integer controls,
`e241_box_genus_bound.py`]** The producer rebuilds all \(8^3=512\) ordered
boxes with \(1\leq a,b,c\leq8\).  On every row it verifies the vertex and edge
counts, simplicity, connectivity, the parity bipartition, all unit-plaquette
counts and edge incidences, the exact quotient/remainder ceiling arithmetic,
and the planarity-classification arithmetic.  Every edge in every
\(\min(a,b,c)\geq2\) row has a validated plaquette witness.  Twelve boundary
and representative rows are pinned literally.  Complete selected witnesses
are stored for \(2\times2\times2\) and \(2\times3\times3\).

The artifact is `results/topology/box_genus_bound.json`.  Its canonical
content digests are:

| content | SHA-256 |
|---|---|
| theorem and proof data | `bea49ac1a6776c5378e2abe0e452e4f2493114fc833b643262907d1583ddd16a` |
| 512-row finite table | `d0fcbebd0f07622db20c24d6bc797337633e878db2310370b0e62ad7c1a71df6` |
| pinned rows | `a8f893521a1c48828e892ccb4803d3e52fb25cbfbd65fd6a46706b41a89552f7` |
| ceiling controls | `0bfbca7704609fd0134a1ce12d71a45bb8280fa54ee22518e0144806d7f68238` |
| cube controls | `4272245f842f2bcb98d1e22f404269378bdaa3180bc1744a320fd64cfca9ae1d` |
| elongated controls | `228f1d420d551bfb3eb08f7953aa5eb4858a9225dddba3b5d28a27d880a0cae5` |
| explicit plaquette witnesses | `203e6c97e5b162f45a7a87c20a1abc505460366b0e9d6d30bc38850a52647991` |
| scope | `84dcc32d347c2db571966554ffec72f3648f722c23c80e8adc3daf7cf1fef09c` |

**[REPRODUCTION — standalone exact verifier]**
`tests/test_box_genus_bound.py` does not import the producer.  It independently
builds each graph from six-neighbour adjacency, enumerates the `xy`, `xz`, and
`yz` plaquettes in separate loops, reconstructs every selected witness and
every table row, checks the literal pins and corollary rows, and verifies all
content and source digests.

### 9g. Scoped surface relevance and exclusions

The classical genus-\(g\) spin-structure construction, when invoked, has
\(2^{2g}=4^g\) terms.  This is only a term count within that construction.  It
is not used in the theorem above and is not a lower bound on the number of
Pfaffians required by an arbitrary formula.

This section excludes periodic length-two multigraphs and all other periodic
identifications.  It does not determine exact genera beyond the planar cases,
does not prove an arbitrary-formula or arbitrary-Pfaffian lower bound, and
does not provide a thermodynamic solution, infinite-volume free energy,
critical point, or critical exponent.  No benchmark value is used.
