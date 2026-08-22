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
