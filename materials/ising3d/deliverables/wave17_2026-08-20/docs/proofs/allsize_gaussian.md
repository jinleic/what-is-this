# An all-size Gaussian no-go for Ising layers with a branching site

Artifacts: `experiments/e128_allsize_gaussian.py`,
`tests/test_allsize_gaussian.py`, `results/spectral/allsize_gaussian.json`.

Reproduce:

```
PYTHONPATH=src .venv/bin/python experiments/e128_allsize_gaussian.py
PYTHONPATH=src .venv/bin/python tests/test_allsize_gaussian.py
```

**Status convention.** Every mathematical statement is tagged **[THEOREM]**,
**[LEMMA]**, **[COMPUTATION]**, **[EXTERNAL]**, or **[UNRESOLVED]**.
**[COMPUTATION]** means an exact integer / finite-field / exact-rational
calculation recorded in the artifact and independently rebuilt by the standalone
test.  Modular gcd degrees are used only in the one-sided direction justified in
§3.3, where they are *upper* bounds for the characteristic-zero gcd degree and
hence *lower* bounds for the collision count — the direction a no-go needs.

The predecessors `proofs/pair_product_obstruction.md` (one layer, one coupling),
`proofs/pair_product_scale.md` (one layer, generic coupling) and
`proofs/parity_pair_product.md` are *finite* certificates: each fixes a layer
size.  The point of this note is to remove the size quantifier.  It succeeds
outright for the inhomogeneous (per-site field, per-bond coupling) layer family
(Theorem F, §7), partially for the uniform-field family (Theorems U and F2,
§10), and fails, in an exactly locatable way, on the fully isotropic curve; §9
states that failure as a theorem about the method, not as a gap in the
write-up.  Along the way §10b closes the `3x3` gap that
`proofs/pair_product_scale.md` §6 left open at a resource wall.

---

## 1. Results

Throughout, `G = (V,E)` is a finite simple graph ("layer graph"), `n = |V|`, and
a *full n-mode Gaussian subset-product multiset* is a multiset of `2^n` numbers
of the form

```
    Lambda(a, u) = { a * prod_{i in S} u_i : S subset of {1..n} },     a, u_i in C,
```

listed with multiplicity `1` for each of the `2^n` subsets.  This is the
spectrum of a free-fermion (Bogoliubov) transfer operator with `n` modes; it is
the object that all three predecessor certificates rule out for a single layer.

### [THEOREM F — family no-go, all sizes]

*Let `G` be any finite simple graph with at least one vertex of degree `>= 3`
(equivalently: `G` is not a disjoint union of paths and cycles; equivalently,
by `proofs/algebraic_obstruction.md` Thm. 5, the frustration graph of the
standard generator set contains an induced claw).  Let `n = |V| >= 4` and let*

```
    R_G(t, w)  =  P(t) D(w) P(t),
    P(t)_{k,l} = prod_{v : k_v != l_v} t_v ,
    D(w)_{k,k} = prod_{e = (uv) : sigma_u(k) = sigma_v(k)} w_e ,
```

*be the inhomogeneous layer operator of §2 (a nonzero scalar multiple of the
physical transfer operator `V` at transverse fields `h_v = 2 artanh t_v` and
bond couplings `J_e = (1/2) log w_e`).  Put `S = C(2^n, 2)`.  Then there is a
nonzero polynomial `rho_G` in `Q[t_v, w_e]` with total degree*

```
    deg rho_G  <=  4 (2n + |E|) * S * (S - 1)                                (1)
```

*such that for every parameter point with `rho_G(t,w) != 0` the `2^n`
eigenvalues of `R_G(t,w)`, listed by algebraic multiplicity, have at least*

```
    3^n - 2^n + 48 * 3^(n-4)                                                 (2)
```

*distinct products over unordered pairs of distinct eigenvalue slots.  Since a
full `n`-mode Gaussian subset-product multiset has at most `3^n - 2^n` such
products (Lemma 1), the spectrum of `R_G(t,w)` — hence of the physical layer
transfer operator — is **not** a full `n`-mode Gaussian subset-product multiset,
for arbitrary complex `a` and `u_1..u_n`.*

The excess `48 * 3^(n-4) = (16/27) * 3^n` is not a boundary effect: it is a
*constant fraction* `16/27 ≈ 59.3%` of `3^n`.  Relative to the ceiling
`3^n - 2^n` itself the ratio is `(16/27)/(1 - (2/3)^n)`, which is `64.96%` at
`n = 6`, `61.67%` at `n = 8` and decreases monotonically to `16/27 = 59.26%` as
`n -> infinity`.

Two immediate corollaries, both uniform in `n`:

* **[COROLLARY F1]** The exceptional set `{rho_G = 0}` is a proper algebraic
  hypersurface, so it has Lebesgue measure zero in `R^{|V|+|E|}`.  For
  Lebesgue-almost-every inhomogeneous (random-field, random-bond) assignment of
  couplings, a layer with a branching site is not spectrally Gaussian.
* **[COROLLARY F2]** By Schwartz–Zippel, for any finite set `T` of rational
  values with `|T| > 4 (2n+|E|) S (S-1)` there is a coupling vector in
  `T^{|V|+|E|}` at which the conclusion holds; in fact all but a fraction
  `deg rho_G / |T|` of them work.

**[SCOPE — read this as part of Theorem F, not as a caveat].**  The quantifier
is *every graph with a branching site and every `n`*, at **anisotropic**
parameters: each site carries its own `t_v` and each bond its own `w_e`, and the
conclusion is asserted only off the hypersurface `{rho_G = 0}`.  The **physical
isotropic layer** — all `t_v` equal to one `t` and all `w_e` equal to
`q(t) = (1+t^2)/(2t)` — is a curve in that parameter space and is **NOT covered**
by Theorem F, at any `n`.  Two reasons, both proved below and neither cosmetic:

1. The specialization lemma (Lemma 7) transfers an **upper** bound on the
   gcd degree from special parameters to the generic point, i.e. a **lower**
   bound on the collision count from special to generic — never the reverse.  A
   generic theorem therefore cannot by itself decide any single named point, and
   the isotropic curve is exactly such a named locus.
2. The obstruction is real and not merely a limitation of the transfer
   direction: **[COMPUTATION]** the very core of the proof, the claw at
   *isotropic* parameters, sits exactly **at** the Gaussian ceiling
   (`r = 65 = 3^4 - 2^4`, `r_all = 73 < 3^4`) at all seven declared couplings,
   so it carries no excess to propagate (§9.3), and equal-field decoupling makes
   the count only linear in the number of free sites (Lemma 9, §9.2).

§10 states the strongest *uniform-field* statement we could prove and the exact
hypothesis that separates it from Theorem F.

### [THEOREM C — the mandatory 1D control, all sizes]

*If `G` is a disjoint union of paths then, for every `t in (0,1)^V` and every
`w in (0,inf)^E`, the spectrum of `R_G(t,w)` **is** a full `n`-mode Gaussian
subset-product multiset.  Consequently its pair-product count is at most
`3^n - 2^n` and the invariant of Theorem F never flags a chain, at any size and
at any coupling.*

[COMPUTATION] The count is *exactly* `3^n - 2^n`, and the "with repetition"
count is exactly `3^n`, at the declared anisotropic rational parameters for
`n = 3,...,8`; see the table in §8.  In particular the invariant sits precisely
on, and never above, the Gaussian ceiling for the 1D control.

### [THEOREM R — bipartite reciprocality, all sizes]

*If `G` is bipartite (every grid layer is), then for all parameters the spectrum
of `R_G(t,w)` is closed under `lambda -> c/lambda` with*

```
    c = prod_v (1 - t_v^2)^2 * prod_e w_e .
```

The underlying matrix identity `(W R_G W^T) R_G = c I` is polynomial and holds at
every parameter point; the *spectral* conclusion needs `R_G` invertible, i.e.
`t_v != +-1` and `w_e != 0`.  Every Gaussian subset-product multiset has this
property (`lambda_S lambda_{S^c} = a^2 prod_i u_i`), so Theorem R says the
reciprocality necessary condition is *automatically satisfied* by every bipartite
layer at every admissible coupling.  Any invariant built only from it is
therefore worthless here.  This is why §3 uses the strictly finer pair-product
count.

### [THEOREM U — isotropic no-go without any anisotropy, `n <= 9`]

*Let `T` be the 5-vertex tree with a degree-3 vertex whose branches have lengths
`2,1,1`.  For `m = 0,1,2,3,4` the layer graph `T u P_m` at **fully uniform
physical couplings** (`t_v = t`, `w_e = q(t)`) is not a full `(5+m)`-mode
Gaussian subset-product multiset, for all `t` outside a finite set of degree
bounded by `4(2n+2|E|) S (S-1)` — in particular at `t = 1/3` and, for
`m <= 3`, also at the named couplings `1/4, 2/5, 1/5`.*

Row `m = 4` was computed at `t = 1/3` only; the `t`-quantifier for every row comes
from Theorem G of §10c, not from the grid.  The proof uses **only** certified
lower bounds and Lemma 1.  It is an isotropic no-go for branching layers at
`n = 5,...,9` with no anisotropy at all.  What is **not** claimed: that the
counts `3893` and `35597` are exact, or that any multiplicative-independence /
injectivity hypothesis is certified — see the warning in §10.

### [THEOREM F2 — uniform-field family theorem, doubly conditional]

*Assume the two hypotheses CZ and INJ(n-5) of §10 (an exact characteristic-zero
value for the core counts, and injectivity of the Lemma 5 labelling map).  Then
every layer graph that decouples into `T u P_{n-5}` when one designated bond class
is switched off satisfies `r >= 3^n - 2^n + 202 * 3^{n-5}` for all `(t, w, y)`
outside an explicit hypersurface, where `t` is the common site field, `w` the
common weight of the surviving bonds and `y` the weight of the cut class.*  In
particular, at fixed generic `(t,w)`, for all but finitely many `y`.  The
degeneration point has **uniform fields and uniform surviving bonds** — only the
cut class is special — so it lives on a physical three-parameter family, not on a
per-site/per-bond one.  Neither hypothesis is verified, and the theorem asserts
nothing about whether the physical value `y = w` is exceptional.  §10.

### [THEOREM 3x3 — a gap in the repository closed]

*The open `3x3` Ising layer at the isotropic coupling `t = 1/3`, `e^{2K} = 5/3`
is **not** a full nine-mode Gaussian subset-product multiset:*
`deg gcd(C_2, C_2') = 59641` *at `p = 1000003`, hence at least
`130816 - 59641 = 71175` distinct pair products against a Gaussian ceiling of
`3^9 - 2^9 = 19171`.
`proofs/pair_product_scale.md` §6 records a resource wall at exactly this
computation and states that no `3x3` conclusion is claimed; §10b completes it.
This is the first layer with a degree-4 vertex certified in the repository.*

### [UNRESOLVED — the isotropic curve]

Theorem F says nothing on the physical curve `t_v = t`, `w_e = q(t) =
(1+t^2)/(2t)`.  §9 proves that this is a genuine obstruction of the *method*,
not an accident of the write-up, and isolates exactly what would repair it.

---

## 2. The inhomogeneous layer model

Let `X_v, Z_v` be Pauli operators on `(C^2)^{tensor n}` and

```
    V(h, J) = exp( (1/2) sum_v h_v X_v ) exp( sum_{e=(uv)} J_e Z_u Z_v )
              exp( (1/2) sum_v h_v X_v ).
```

With `t_v = tanh(h_v/2)` and `w_e = exp(2 J_e)`,

```
    exp( (1/2) sum_v h_v X_v )  =  ( prod_v cosh(h_v/2) ) * P(t),
    exp( sum_e J_e Z_u Z_v )    =  ( prod_e w_e^{-1/2} ) * D(w),
```

because `exp(theta X) = cosh(theta) (1 + tanh(theta) X)` and
`exp(J s s') = e^{-J} w^{(1+ss')/2}`.  Hence

```
    V(h, J) = kappa * R_G(t, w),   kappa = prod_v cosh^2(h_v/2) * prod_e w_e^{-1/2} != 0.   (3)
```

A nonzero global scalar multiplies every pair product by `kappa^2` and therefore
does not change any collision count; similarity does not either.  So it suffices
to work with `R_G`.

**Isotropic (physical 3D) curve.**  Setting all `t_v = t` and all `w_e = q(t) =
(1+t^2)/(2t)` reproduces the canonical construction of
`experiments/e38_gaussianity_certificate.py`:

**[COMPUTATION]** For the open `2x3` and `2x2` layers at `t = 1/3` the producer
verifies entrywise that `R_G(t, q) = q^{(|E|+eps)/2} * R_{e38}`, with the
observed ratios `625/81 = (5/3)^4` and `25/9 = (5/3)^2`.  So the family below
literally contains the operator all previous certificates used.

**Integral model.**  Writing `t_v = a_v/b_v` and `w_e = c_e/d_e`, the matrix

```
    M = ( prod_v b_v )^2 ( prod_e d_e ) * R_G(t,w)
      = ( kron_v [[b_v, a_v],[a_v, b_v]] ) * diag(D_int) * ( kron_v [[b_v, a_v],[a_v, b_v]] )
```

is integral **by construction** (no denominator search is needed), where
`D_int(k) = prod_e (c_e if e aligned in k else d_e)`.  The producer builds `M`
this way in `O(n 4^n)` integer operations.  The standalone test never imports the
producer: it rebuilds every case by an independent triple-loop rational
`P D P` product followed by an actual-LCM denominator clearing, recomputes the
invariants, and only then compares the recomputed numbers with the stored
artifact.

---

## 3. The invariant

For a multiset `Lambda` of `N` eigenvalue **slots** (repeated by algebraic
multiplicity), set

```
    r(Lambda)     = #{ distinct values lambda_i lambda_j , i < j },
    r_all(Lambda) = #{ distinct values lambda_i lambda_j , i <= j }.
```

Slots, not values: two numerically equal eigenvalues still give two slots.  This
is what makes the invariant ordering-free and degeneracy-proof.

### [LEMMA 1 — Gaussian ceilings]

If `Lambda` is a full `n`-mode subset-product multiset then

```
    r(Lambda) <= 3^n - 2^n        and        r_all(Lambda) <= 3^n .
```

*Proof.* Match the `2^n` slots bijectively with the subsets `S`, i.e. with bit
vectors `eps in {0,1}^n`, even when displayed values coincide.  A pair of slots
`(eps, delta)` has product `a^2 prod_i u_i^{eps_i + delta_i}`, whose exponent
vector `f = eps + delta` lies in `{0,1,2}^n`; there are `3^n` such vectors.  For
two *distinct* slots, `eps != delta`, so some coordinate of `f` equals `1`, and
the number of such `f` is `3^n - 2^n`.  Equal or zero parameters, repeated
eigenvalues and accidental multiplicative relations can only merge displayed
values, never create a new exponent class. `[]`

### [LEMMA 2 — derivative gcd counts distinct roots]

Over a field of characteristic zero, a monic `P` of degree `N` with `r` distinct
roots in an algebraic closure satisfies `deg gcd(P, P') = N - r`.

*Proof.* As in `proofs/pair_product_obstruction.md` Lemma 2 (the common-root
multiplicity at a root of multiplicity `m` is exactly `m-1` because `m != 0` in
characteristic zero, including at the root `0`; gcd degree is stable under field
extension). `[]`

### 3.3 [LEMMA 3 — exact reconstruction and one-sided modular reduction]

For an integer matrix `M` of size `d`, let

```
    C_2(z)   = charpoly( Lambda^2 M )   (degree C(d,2),   roots rho_i rho_j, i<j),
    C_all(z) = charpoly( Sym^2 M )      (degree C(d+1,2), roots rho_i rho_j, i<=j).
```

Both are monic in `Z[z]`.  For every prime `p`,

```
    deg gcd_{F_p}( C, C' )  >=  deg gcd_{Q}( C, C' )                        (4)
```

for `C` either of them, so a modular computation yields the *lower* bounds
`r >= C(d,2) - deg gcd_{F_p}(C_2, C_2')` and
`r_all >= C(d+1,2) - deg gcd_{F_p}(C_all, C_all')`.

*Proof.* Integrality and monicity: `Lambda^2 M` and `Sym^2 M` are integer
matrices and triangularizing `M` over an algebraic closure exhibits the stated
eigenvalues with algebraic multiplicity.  (4): the monic rational gcd `g` and
both cofactors are integral by Gauss's lemma, so reducing the two factorizations
modulo `p` keeps `bar g` a common monic divisor of the same degree. `[]`

**[COMPUTATION — the reconstruction pipeline].**  The producer never forms
`Lambda^2 M` or `Sym^2 M` explicitly for `d > 16`.  It computes `s_k = tr(M^k)`
for `k <= d` by exact modular matrix products (float64 BLAS with the explicit
guard `d (p-1)^2 < 2^53`, so the products are exactly representable), gets
`charpoly(M)` by Newton's identities, extends `s_k` through `2 C(d+1,2)` by the
Cayley–Hamilton trace recurrence, and then uses

```
    tr( (Lambda^2 M)^m ) = (s_m^2 - s_{2m})/2 ,
    tr( (Sym^2 M)^m )    = (s_m^2 + s_{2m})/2 ,
```

with Newton's identities again, followed by an explicit Euclidean gcd in
`F_p[z]`.  Cayley–Hamilton is verified by matrix Horner evaluation, the
recurrence is spot-checked against binary matrix powers, and the returned gcd is
verified to divide both inputs with zero remainder.  The standalone test
recomputes the decisive `16`-dimensional case *twice*: once through this
pipeline and once through the characteristic polynomials of the **explicitly
assembled** `120 x 120` exterior square and `136 x 136` symmetric square, and
requires the two polynomials to be equal coefficient by coefficient.

---

## 4. The family engine: exact tensor decoupling

Fix a vertex `v*` of degree `>= 3` and three of its edges; let
`C = {v*, x, y, z}` be the four endpoints and `E_C` the three chosen edges, so
`(C, E_C)` is a claw `K_{1,3}`.  Write `m = n - 4`.

### [LEMMA 4 — the decoupled operator is a tensor product]

Specialize `w_e -> 1` for every `e` not in `E_C`.  Then

```
    R_G  =  R_C(t_C, w_C)  tensor  ( kron_{v not in C} [[1, t_v],[t_v, 1]]^2 ),
```

and `[[1,t],[t,1]]^2` has eigenvalues `(1+t)^2` and `(1-t)^2`.

*Proof.* `w_e = 1` removes `e` from the product defining `D(w)`, so `D` depends
only on the spins in `C`: `D = D_C tensor I`.  `P(t)` is a Kronecker product by
definition.  Hence `R_G = P D P = (P_C D_C P_C) tensor (kron_{v not in C} P_v^2)`. `[]`

### [LEMMA 5 — exact pair-product count of a decoupled configuration]

Let the `m` free parameters `t_v` (`v not in C`) be independent indeterminates
over `Qbar`, let `Lambda_C` be the `16` eigenvalue slots of `R_C` (all nonzero;
`R_C` is positive definite for `0 < t < 1`, `w > 0`), and let `Lambda` be the
`2^n` eigenvalue slots of the decoupled operator of Lemma 4.  Then, exactly,

```
    r(Lambda)     = (3^m - 2^m) * r_all(Lambda_C)  +  2^m * r(Lambda_C),      (5)
    r_all(Lambda) = 3^m * r_all(Lambda_C).                                    (6)
```

*Proof.* Slots of `Lambda` are pairs `(s, delta)` with `s` one of the 16 claw
slots and `delta in {0,1}^m`; the eigenvalue is `alpha_s * prod_v c_{v,delta_v}`
with `c_{v,0} = (1+t_v)^2`, `c_{v,1} = (1-t_v)^2`.  A pair of slots
`(s,delta) != (s',delta')` has product

```
    alpha_s alpha_{s'} * prod_v (1+t_v)^{p_v} (1-t_v)^{q_v},
    (p_v, q_v) = (4,0), (2,2), (0,4)  according as  f_v := delta_v + delta'_v = 0, 1, 2.
```

In `Qbar[t_v : v not in C]` the linear forms `1+t_v`, `1-t_v` are pairwise
non-associate primes and the `alpha`'s are units, so unique factorization
recovers `f` from the value, and then recovers `alpha_s alpha_{s'}`.  Hence the
value determines the pair `(f, alpha_s alpha_{s'})` and conversely.  Now split:

* `f` has at least one coordinate `1` (there are `3^m - 2^m` such `f`).  Then
  `delta != delta'` already, so `s` and `s'` range over *all* pairs including
  `s = s'`; the number of values contributed is `r_all(Lambda_C)`.
* `f in {0,2}^m` (there are `2^m` such `f`).  Then `delta = delta'` is forced,
  so `s != s'`; the number of values contributed is `r(Lambda_C)`.

Summing gives (5).  Allowing also `(s,delta) = (s',delta')` adds nothing new to
the first bullet and upgrades the second to "all pairs", giving (6). `[]`

Lemma 5 is validated three ways in the artifact: (i) on synthetic multisets
built from pairwise coprime integers, where the count can be done by brute force
and unique factorization guarantees the hypothesis (`m = 1,2,3`); (ii) against
*direct* modular computation of the full decoupled layer at rational (not
indeterminate) free fields, for `n = 5,6,7,8,9`, which reproduces (5) and (6)
*exactly* (`355, 1097, 3355, 10193, 30835`), so the rational specializations
already realize the generic count; (iii) by the internal consistency check
`r_all - r <= 16`.

Dropping the independence hypothesis turns (5) into an inequality that is valid
unconditionally and is used in §10: **Lemma 5'** says
`r <= (3^m - 2^m) r_all(Lambda_H) + 2^m r(Lambda_H)` always, with equality iff
the hypothesis holds.  Together with a modular lower bound this is a
*two-sided* tool: whenever the two meet, `r` is pinned exactly and independence
is certified rather than assumed.

### [LEMMA 6 — excess propagation, general core]

Let `H` be any core on `n_H` sites, `E := r_all(Lambda_H) - 3^{n_H}`.  Then for
the decoupled configuration of Lemma 4 with `m = n - n_H` free sites,

```
    r(Lambda) - (3^n - 2^n)  =  3^m * E  -  2^m * ( (r_all - r)(Lambda_H) - 2^{n_H} )
                             >=  3^m * E ,                                    (7)
```

because `r_all(Lambda_H) - r(Lambda_H) <= 2^{n_H}`: the extra values counted by
`r_all` all lie in the set of at most `2^{n_H}` squares `alpha_s^2`.

*Proof.* Substitute (5) and `3^n - 2^n = 3^{n_H} 3^m - 2^{n_H} 2^m` and
rearrange. `[]`

So **one** core whose "with repetition" count strictly exceeds `3^{n_H}` forces
an exponentially growing excess at every larger size.

---

## 5. From one point to the generic point

### [LEMMA 7 — monic specialization]

Let `A(t,w;z) = charpoly_z( Lambda^2 M(t,w) )`, monic in `z` of degree `S` with
coefficients in `Z[t,w]`, let `G` be the monic
`gcd_{Q(t,w)[z]}(A, d A/dz)` and `g = deg_z G`.  For every ring homomorphism
`psi : Q[t,w] -> R` into an integral domain of characteristic zero,

```
    deg_z gcd( psi A, psi (dA/dz) )  >=  g .
```

*Proof.* `A` is monic in `z`, so `G`, `A/G` and `(dA/dz)/G` all lie in
`Q[t,w][z]` by Gauss's lemma; applying `psi` keeps `psi G` monic of degree `g`
and a common divisor. `[]`

Consequently, with `r_gen := S - g` the number of distinct pair products at the
generic point,

```
    r_gen  >=  r(Lambda_specialized)   for every specialization.              (8)
```

### [LEMMA 8 — an explicit exceptional set]

Put `A_0 = A/G` and `B_0 = (dA/dz)/(S G)`.  Since `A` is monic of degree `S`, the
derivative `dA/dz` has leading coefficient `S`, so this normalisation (and not
`(dA/dz)/G`) is what makes `B_0` monic of degree `n_0 - 1`; dividing by the
nonzero integer `S` changes no `t`- or `w`-degree, so every bound below is
unaffected.  `A_0` and `B_0` are coprime over `Q(t,w)`; set
`rho_G = Res_z(A_0, B_0) in Q[t,w]`.  Then `rho_G != 0`; for every parameter
point `theta` with `rho_G(theta) != 0` one has
`deg gcd(A(theta,z), dA/dz(theta,z)) = g`, hence exactly `r(theta) = r_gen`;
and

```
    deg rho_G  <=  4 (2n + |E|) * S (S-1).
```

*Proof.* Non-vanishing and the equality clause are the standard resultant
argument of `proofs/pair_product_scale.md` Lemma 6.  For the degree, restrict to
a generic affine line `theta(tau) = theta_0 + tau theta_1`: since `resultants of
monic polynomials commute with specialization`, `rho_G` restricted to the line is
the resultant of the restrictions, and a generic line realizes the total degree.
Along the line, entries of `P` have `tau`-degree `<= n`, entries of `D` have
`tau`-degree `<= |E|`, so entries of `M` have degree `<= d_R := 2n + |E|` and
entries of `Lambda^2 M` degree `<= 2 d_R`; the `k`-th coefficient of `A` is a sum
of `k x k` minors, of degree `<= 2 d_R k`.  Substituting `z = tau^{2 d_R} y` and
dividing by `tau^{2 d_R S}` makes `A` monic in `y` with coefficients in the UFD
`Q[tau^{-1}]`; Gauss's lemma gives the same indexed bound `2 d_R i` for the
`i`-th coefficient of every monic factor, and monic long division propagates it
to `B_0`.  The Sylvester determinant of `A_0` (degree `n_0 <= S`) and `B_0`
(degree `n_0 - 1`) is homogeneous of degree `n_0 - 1` in the first set of
coefficients and `n_0` in the second, so its degree is at most
`2 * 2 d_R * n_0 (n_0 - 1) <= 4 d_R S (S-1)`. `[]`

---

## 6. The claw certificate — the only decisive finite input

**[COMPUTATION — claw].**  At the predeclared rational point

```
    (t_{v*}, t_x, t_y, t_z) = (1/3, 1/5, 1/7, 2/11),
    (w_{v*x}, w_{v*y}, w_{v*z}) = (5/3, 7/2, 11/5),
```

the `16 x 16` operator `R_{K_{1,3}}` gives, at **both** primes
`p_1 = 1000003` and `p_2 = 2000003`,

| core (4 sites, 3 edges) | pair slots | `deg gcd(C_2, C_2')` | `r >=` | sym slots | `deg gcd(C_all, C_all')` | `r_all >=` |
|---|---:|---:|---:|---:|---:|---:|
| claw `K_{1,3}` | 120 | 7 | **113** | 136 | 7 | **129** |
| path `P_4` (same parameters) | 120 | 55 | 65 | 136 | 55 | 81 |

The Gaussian ceilings at `n = 4` are `3^4 - 2^4 = 65` and `3^4 = 81`.  So:

* the claw **exceeds** both ceilings: `113 > 65` and `129 > 81`;
* the path `P_4`, at the *same* four fields and the *same* three bond weights,
  attains both ceilings **exactly** and never exceeds them.

The two graphs have the same number of vertices, the same number of edges and
the same parameters.  The single difference is the degree-3 vertex.  This is
precisely the claw of `proofs/algebraic_obstruction.md` Theorem 5(b): the
frustration graph of `{X_i} u {Z_iZ_j}` is the subdivision `S(Gamma)`, which
contains an induced claw iff `Gamma` has a vertex of degree `>= 3`.

The excess constant is

```
    E = r_all(claw) - 3^4 = 129 - 81 = 48 .
```

Also recorded: `r_all - r = 129 - 113 = 16 = 2^4`, saturating the inequality used
in Lemma 6.

---

## 7. Proof of Theorem F

Let `G` have a vertex `v*` of degree `>= 3`; choose three incident edges and let
`C` be as in §4, `m = n - 4`.  Define the specialization
`psi : Q[t,w] -> Q[t_v : v not in C]`:

```
    w_e -> 1                    for every e not in E_C,
    (t_C, w_C) -> the rational values of §6,
    t_v -> t_v                  for v not in C   (kept as indeterminates).
```

By Lemma 4 the specialized operator is `R_C tensor (kron_v P_v^2)`.  By Lemma 5
and the certificate of §6,

```
    r( psi Lambda ) = (3^m - 2^m) r_all(Lambda_C) + 2^m r(Lambda_C)
                    >= (3^m - 2^m) * 129 + 2^m * 113
                    =  129 * 3^m - 16 * 2^m
                    =  (3^n - 2^n) + 48 * 3^{n-4} .                            (9)
```

(The last line is the exact identity `129 * 3^m - 16 * 2^m = 81 * 3^m - 16 * 2^m
+ 48 * 3^m`.)  By Lemma 2 applied over the characteristic-zero field
`Q(t_v : v not in C)`, this is the statement
`deg gcd(psi A, psi A') = S - r(psi Lambda) <= S - (3^n - 2^n) - 48 * 3^{n-4}`.
Lemma 7 transports it to the generic point:

```
    g  <=  S - (3^n - 2^n) - 48 * 3^{n-4},     hence   r_gen >= (3^n-2^n) + 48*3^{n-4}.
```

Lemma 8 produces `rho_G` with the degree bound (1) such that `r(theta) = r_gen`
whenever `rho_G(theta) != 0`.  Finally Lemma 1 says a full `n`-mode Gaussian
subset-product multiset can have at most `3^n - 2^n` distinct pair products, and
`48 * 3^{n-4} > 0`.  By (3) the conclusion transfers from `R_G` to the physical
`V(h,J)`. `[]`

**Remark (general cores).**  Nothing used the claw beyond `r_all > 3^{n_H}`.  By
Lemma 6, any core `H` on `n_H` sites, occurring in `G` as a chosen vertex subset
plus a chosen edge subset, with certified `r_all(Lambda_H) = 3^{n_H} + E`, yields
the same theorem with excess `3^{n - n_H} E`.  The claw is used because it is the
*smallest* core that occurs in **every** graph with a branching vertex.

**Remark (why `r_all`, not `r`, is the right core quantity).**  Lemma 6 shows
that the propagated excess is governed by `r_all - 3^{n_H}` and not by
`r - (3^{n_H} - 2^{n_H})`.  A core with excess only in `r` would be diluted away.
The claw has excess in both (`113 > 65`, `129 > 81`), and it is the `r_all`
excess `48` that survives to all sizes.

---

## 8. Verification tables

All rows below are in `results/spectral/allsize_gaussian.json` and are
recomputed by the standalone test.

### 8.1 The family bound at every requested size

`r_family(n) = (3^{n-4} - 2^{n-4}) * 129 + 2^{n-4} * 113`; excess `= 48 * 3^{n-4}`.

| layer | `n` | `|E|` | max deg | Gaussian ceiling `3^n - 2^n` | family lower bound | excess |
|---|---:|---:|---:|---:|---:|---:|
| `2x3` | 6 | 7 | 3 | 665 | 1097 | 432 |
| `2x4` | 8 | 10 | 3 | 6305 | 10193 | 3888 |
| `3x3` | 9 | 12 | 4 | 19171 | 30835 | 11664 |
| `2x5` | 10 | 13 | 3 | 58025 | 93017 | 34992 |
| `3x4` | 12 | 17 | 4 | 527345 | 842273 | 314928 |
| `4x4` | 16 | 24 | 4 | 42981185 | 68490353 | 25509168 |
| `5x5` | 25 | 40 | 4 | 847255055011 | 1349352008755 | 502096953744 |
| `10x10` | 100 | 180 | 4 | `3^100 - 2^100` (48 digits) | `129*3^96 - 16*2^96` | `48 * 3^96` |

The excess column is always exactly `48 * 3^{n-4}`.  The artifact stores every
entry as an exact integer, including the 48-digit `10x10` row and the
exceptional-set degree bound `4 (2n + |E|) S (S-1)` per row (for example
`308,730,240` for `2x3`), and the standalone test recomputes every row
independently.  `2x3`, `2x4` and `3x3` are the layers this repository had
attacked with the pair-product invariant (`2x3` and `2x4` certified at single
couplings; `3x3` walled out, and now certified in §10b).  `3x4`, `4x4`, `5x5`,
`10x10` are sizes no pair-product certificate in the repository reaches, and
`2x5` was previously reached only by the different, ordering-dependent
Theorem-S route (`results/spectral/spectral_2x5.json`).  The family theorem
covers all of them, and every other size, because the only hypothesis is
`max degree >= 3`.

### 8.2 Direct verification of the specialized count (illustration of Lemma 5)

| `n` | free sites `m` | direct `r` | Lemma 5 | direct `r_all` | Lemma 5 | `3^n - 2^n` | excess |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 1 | 355 | 355 | 387 | 387 | 211 | 144 |
| 6 | 2 | 1097 | 1097 | 1161 | 1161 | 665 | 432 |
| 7 | 3 | 3355 | 3355 | 3483 | 3483 | 2059 | 1296 |
| 8 | 4 | 10193 | 10193 | 10449 | 10449 | 6305 | 3888 |
| 9 | 5 | 30835 | 30835 | — | — | 19171 | 11664 |

These are *not* the theorem; they are direct exact recomputations of one
specialization, at rational free fields, confirming that Lemma 5's count is
attained already over `Q` and not merely over the function field.

### 8.3 The 1D control (Theorem C)

Open chain `P_n` at the declared anisotropic rational parameters, prime `p_1`:

| `n` | `r` | `3^n - 2^n` | `r_all` | `3^n` | flagged? |
|---:|---:|---:|---:|---:|---|
| 3 | 19 | 19 | 27 | 27 | no |
| 4 | 65 | 65 | 81 | 81 | no |
| 5 | 211 | 211 | 243 | 243 | no |
| 6 | 665 | 665 | 729 | 729 | no |
| 7 | 2059 | 2059 | 2187 | 2187 | no |
| 8 | 6305 | 6305 | 6561 | 6561 | no |

The chain attains the Gaussian ceiling exactly at every tested size and never
exceeds it, as Theorem C requires.  (Equality also shows that at these rational
points the chain's `n` modes are already multiplicatively independent, so the
test is not passing for a degenerate reason.)

**Proof of Theorem C.**  Order the vertices along each path and use the
Jordan–Wigner map `X_v = -i gamma_{2v-1} gamma_{2v}`,
`Z_v Z_{v+1} = -i gamma_{2v} gamma_{2v+1}` (the strings cancel for nearest
neighbours in the path order).  Both `sum_v h_v X_v` and `sum_e J_e Z_uZ_v` are
then quadratic forms in the `2n` Majoranas, so by §2 `R_G` is a product
`e^{Q_1} e^{Q_2} e^{Q_1}` of exponentials of Majorana bilinears and therefore
lies in the spin group; `Ad_{R_G}` is the corresponding element of the **complex**
orthogonal group `O(2n,C)`, whose eigenvalues occur in pairs
`{mu_k, mu_k^{-1}}` because it preserves the symmetric bilinear form.
**[EXTERNAL, standard]** For the open transverse-field Ising chain this element
is diagonalizable and a (complex) Bogoliubov transformation brings `R_G` to
`a * exp( sum_k (log mu_k) d_k^dagger d_k )`, whose spectrum is
`{ a prod_k mu_k^{eps_k} }`: this is the classical free-fermion diagonalization
of the open chain (Lieb–Schultz–Mattis, *Ann. Phys.* **16** (1961) 407;
Schultz–Mattis–Lieb, *Rev. Mod. Phys.* **36** (1964) 856), which we cite rather
than re-derive.  (An earlier draft of this note asserted that `Ad_{R_G}` is a
symmetric positive-definite element of `SO(2n,R)`; that is false — a real matrix
that is both symmetric and orthogonal has eigenvalues `+-1` — and the real form
plays no role in the argument.)  Lemma 1 then gives `r <= 3^n - 2^n`.  The table
above is the exact confirmation at `n = 3..8`. `[]`

### 8.4 Bipartite reciprocality (Theorem R)

*Proof.* Let `(P, V\P)` be a bipartition and
`W = ( prod_{v in P} X_v ) ( prod_{v in V} Z_v )`.  Conjugation by
`prod_v Z_v` sends `X_v -> -X_v` and fixes `Z_v`; conjugation by
`prod_{v in P} X_v` fixes `X_v` and sends `Z_v -> -Z_v` exactly for `v in P`, so
it flips the sign of every `Z_uZ_v` with `(uv) in E` (each edge has exactly one
endpoint in `P`).  Hence `W` sends `A -> -A` and `B -> -B`, and therefore

```
    (W R_G W^T) R_G  =  c I ,        c = prod_v (1-t_v^2)^2 * prod_e w_e ,     (R)
```

an identity between matrices whose entries are polynomials in `(t, w)`, hence
valid at **every** parameter point (including `t_v = +-1` and `w_e = 0`).  The
*spectral* consequence needs `R_G` invertible: for `t_v != +-1` and `w_e != 0`
the factors `P(t)` and `D(w)` are invertible, `c != 0`, and (R) says
`W R_G W^{-1} = c R_G^{-1}`, so `spec R_G = c / spec R_G` as multisets.  At
`t_v = +-1` the block `[[1,t_v],[t_v,1]]` is singular and at `w_e = 0` the
diagonal is singular; there (R) still holds but carries no reciprocality
statement. `[]`

**[COMPUTATION]** The producer verifies the *exact rational matrix identity*
`(W R W^T) R = c I` for the claw, for `P_4`, and for the `2x3` layer both at
anisotropic and at isotropic parameters; the observed centres agree with the
closed form (e.g. `c = 2392326144/285333125` for the claw parameters of §6, and
`c = 2736521232777216/5102703580975` for the anisotropic `2x3`).  A
non-bipartite graph (triangle) is rejected by the bipartiteness test, as the
lemma requires.

---

## 9. The isolated obstruction: the isotropic curve

Theorem F is a statement about generic *inhomogeneous* couplings.  The physical
3D Ising layer lives on the curve `t_v = t`, `w_e = q(t) = (1+t^2)/(2t)`.  The
following three facts, all exact, show that the obstruction is a property of the
localization method itself.

### 9.1 [COMPUTATION] Switching off the transverse bonds restores Gaussianity

For the open `2x3` layer with the three rung bonds set to `w = 1` (the layer
splits into two open 3-chains), at prime `p_1` and `p_2`:

| `2x3` configuration | `deg gcd(C_2,C_2')` | certified `r >=` | Gaussian ceiling | flagged? |
|---|---:|---:|---:|---|
| isotropic `t=1/3`, all bonds on | 385 | 1631 | 665 | **yes** |
| isotropic `t=1/3`, rungs off | 1899 | 117 | 665 | no |
| anisotropic, rungs off | **1351** | **665** | 665 | no |
| anisotropic, all bonds on | 31 | 1985 | 665 | **yes** |

Row 1 independently reproduces the repository value
`deg gcd = 385`, `r >= 1631` of `proofs/pair_product_obstruction.md`, from a
different construction route.  Row 3 shows the decoupled layer sitting *exactly*
on the Gaussian floor `C(64,2) - (3^6-2^6) = 2016 - 665 = 1351` — exact because
the two decoupled chains are Gaussian, so `665` is simultaneously a certified
lower bound and the Lemma 1 upper bound; the six modes are therefore
multiplicatively independent at that point.  Row 2 shows that with *isotropic*
parameters the two decoupled chains are identical, their modes coincide
pairwise, and the certified count drops far *below* the Gaussian ceiling
(`117 < 665`, gcd `1899 > 1351`).

**Conclusion.** No no-go can be uniform in the transverse coupling: at zero
transverse coupling the spectrum *is* Gaussian.  Any family theorem must exclude
a neighbourhood of, or at least the locus containing, the decoupled point — as
Theorem F does through `rho_G`.

### 9.2 [LEMMA 9 — equal-field decoupling can never be flagged]

*Let `Lambda = Lambda_H x Lambda_free` where the `m` free sites all carry the
**same** field parameter, so `Lambda_free = { c u^j } with multiplicity C(m,j)`.
Then*

```
    r(Lambda)  <=  (2m + 1) * r_all(Lambda_H)  <=  (2m+1) * 2^{n_H - 1} (2^{n_H} + 1),
```

*which is linear in `m`, whereas the Gaussian ceiling `3^{n_H + m} - 2^{n_H+m}`
is exponential.  Hence for all `m` beyond an explicit threshold the invariant
provably does **not** flag an equal-field decoupled configuration.*

*Proof.* Every pair product is `alpha_s alpha_{s'} u^{j+j'}` with
`j + j' in {0,...,2m}`. `[]`

For `n_H = 4` the *unconditional* bound `136 (2m+1)` (using
`r_all <= C(2^4+1,2) = 136`) is compared with the ceiling below; it drops below
it exactly from `m = 3` on.  **[COMPUTATION]** with the isotropic claw as core,
`m = 1,2,3,4` extra equal-field sites give

| `m` | `n` | certified `r >=` | proved upper bound `136 (2m+1)` | ceiling `3^n - 2^n` | decisive? |
|---:|---:|---:|---:|---:|---|
| 1 | 5 | 203 | 408 | 211 | no |
| 2 | 6 | 349 | 680 | 665 | no |
| 3 | 7 | 495 | 952 | 2059 | **yes** |
| 4 | 8 | 641 | 1224 | 6305 | **yes** |

From `m = 3` on the *proved* upper bound is below the ceiling, so Lemma 9
**proves** that no equal-field decoupled configuration with a 4-site core can be
flagged, however large the layer.  (The certified lower bounds in column 3 grow
by exactly `146` per extra site, consistent with the linear behaviour, but they
are lower bounds and are not used for the impossibility statement.)

**Consequence.**  A localization proof must (i) delete edges, which forces
`w_e = 1 != q(t)`, and (ii) keep the free part's Gaussian modes distinct, which
equal fields cannot do — and by Lemma 9 the failure of (ii) is not a technical
nuisance but a proved impossibility.  Both requirements leave the isotropic
curve.  §10 removes obstruction (ii) and shows precisely how much of the
isotropic statement that buys.

### 9.3 [COMPUTATION] The isotropic claw yields no certificate

At every coupling of the declared grid `t in {4/20, 5/20, ..., 10/20}` the
isotropic claw gives

```
    r >= 65 = 3^4 - 2^4   and   r_all >= 73 < 81 .
```

Both are *lower* bounds (Lemma 3 is one-sided), so the correct reading is: **the
invariant returns no certified excess for the isotropic claw.**  It does *not*
say the true count is `65`, and it does *not* say the isotropic claw has no
excess — that would require an upper bound, which we do not have.  What it does
say is decisive for the method: the claw-localized proof of §7 needs a certified
`r_all(core) > 3^{n_core}` to propagate, and isotropically there is none.

**[COMPUTATION] Where the isotropic obstruction does start.**  At `t = 1/3`, on
5 sites:

| isotropic 5-site graph | max deg | certified `r >=` | `3^5-2^5 = 211` | certified flagged? |
|---|---:|---:|---:|---|
| `K_{1,4}` star | 4 | 147 | 211 | no (non-certificate) |
| `T` = claw with one branch of length 2 | 3 | **417** | 211 | **yes** |
| claw `+` isolated site | 3 | 203 | 211 | no (non-certificate) |
| paw (triangle `+` pendant) | 3 | **373** | 211 | **yes** |
| path `P_5` | 2 | 211 | 211 | no — and here `r = 211` **exactly**, because a path is Gaussian so Lemma 1 supplies the matching upper bound |
| cycle `C_5` | 2 | 141 | 211 | no (non-certificate) |

The smallest graph with a branching site for which the isotropic invariant
*returns a certificate* is the 5-vertex tree `T` (a degree-3 vertex whose
branches have lengths `2,1,1`).  Stars return none; a plausible reason is their
large automorphism group, but that is an explanation, not a proof of absence.
This already makes the isotropic situation different from Theorem F: the core
that works changes.

### 9.4 What would repair it

Two concrete, precisely stated repairs:

* **(R1) A non-localization lower bound on the isotropic curve.**  One needs
  `r >= 3^n - 2^n + 1` at one point of the curve *for every `n`*, obtained
  without decoupling.  The two natural degenerations of the curve are
  `t -> 0` (where `V` becomes the diagonal classical Ising layer and the
  eigenvalue Puiseux data are governed by degenerate perturbation theory on
  equal-energy spin configurations) and `t -> 1` (where `q -> 1`, `V` becomes
  `exp(K^* A)` and the first-order splitting inside each `X`-basis Hamming-weight
  sector is the hard-core-boson hopping operator of `G` restricted to `k`
  particles).  Since distinct Puiseux germs of `lambda_i lambda_j` give distinct
  *functions*, a lower bound on the number of distinct germs at either endpoint
  is a valid lower bound for `r` at generic `t` on the curve.  We did not carry
  this out; it is a genuinely different technique.
* **(R2) One nonvanishing certificate per size.**  Theorem F plus
  `rho_G(t, q(t)) != 0` for one `t` gives the isotropic statement for that `G`.
  That is available for `n = 6` (`r >= 1631 > 665`, §9.1 and
  `proofs/pair_product_obstruction.md`) and `n = 8` (`r >= 23311 > 6305`,
  `proofs/pair_product_scale.md`), which is exactly the per-size character we set
  out to remove.  See §10 for what the present run adds.

A third repair, suggested during the run, *does* work partially and is developed
in §10: keep the fields uniform but replace the decoupled free **sites** by a
decoupled open **chain**, whose Gaussian modes are pairwise distinct.

---

## 10. The uniform-field route: core plus open chain

Lemma 9 fails because `m` free sites at one field give `m` equal modes.  An open
chain `P_m` at the same uniform field is Gaussian (Theorem C) with `m` pairwise
distinct modes, so Lemma 5 might apply verbatim.  This section carries that out
exactly.  Two ingredients are new.

### [LEMMA 5' — the tensor count is always an UPPER bound]

*For any `Lambda = Lambda_H x Lambda_free` with `Lambda_free` a full `m`-mode
Gaussian,*

```
    r(Lambda)  <=  (3^m - 2^m) r_all(Lambda_H) + 2^m r(Lambda_H),
```

*with equality if and only if the map `(f, alpha_s alpha_{s'}) -> value` of Lemma
5 is injective.  (Equality is exactly injectivity of that finite labelling map on
`{0,1,2}^m x` core pair products; it is **not** the same as multiplicative
independence of the free modes modulo the core pair-product group, and we do not
claim the latter anywhere.)*

*Proof.* Every pair product is `alpha_s alpha_{s'} prod_k u_k^{f_k}` and the
counting of Lemma 5 enumerates the possible pairs `(f, alpha_s alpha_{s'})`;
collisions can only merge them. `[]`

**[WARNING — the direction that does NOT work].**  Lemma 5' bounds `r(Lambda)` by
an expression in the **characteristic-zero** core counts `r_all(Lambda_H)` and
`r(Lambda_H)`.  Below we only have *lower* bounds for those (`>= 445` and
`>= 417`).  Substituting lower bounds into that expression produces a lower
estimate of an unknown upper bound, which bounds nothing.  Consequently a
coincidence between the certified composite count and that substituted value
proves **nothing** — in particular it does not pin `r`, and it does not certify
the injectivity hypothesis.  An earlier draft of this note asserted the
opposite; the resulting "Theorem U2" is **withdrawn** and appears below only as a
conditional statement.  Everything else in §10 uses certified lower bounds only.

### [COMPUTATION] Core `T` and chain at uniform physical couplings

`T` is the 5-vertex tree of §9.3, with certified `r >= 417 > 211` and
`r_all >= 445 > 243`, so `E_T := r_all(T) - 3^5 >= 202`.  The following are
**named-point** certificates: rows `m = 0,1,2,3` were computed at
`t in {1/3, 1/4, 2/5, 1/5}` and are identical at all four of those couplings;
row `m = 4` was computed at `t = 1/3` **only**.  The "reference" column
substitutes the certified core lower bounds into the tensor expression and is,
per the warning above, neither an upper nor a lower bound — it is recorded only
because the certified count happens to equal it at `m = 2, 4`.

| `m` | `n = 5+m` | certified `r >=` | certified excess | `3^n - 2^n` | reference value | equal to reference? | couplings computed |
|---:|---:|---:|---:|---:|---:|---|---|
| 0 | 5 | 417 | 206 | 211 | — | — | 1/3, 1/4, 2/5, 1/5 |
| 1 | 6 | 1219 | 554 | 665 | 1279 | no | 1/3, 1/4, 2/5, 1/5 |
| 2 | 7 | 3893 | 1834 | 2059 | 3893 | yes | 1/3, 1/4, 2/5, 1/5 |
| 3 | 8 | 7885 | 1580 | 6305 | 11791 | no | 1/3, 1/4, 2/5, 1/5 |
| 4 | 9 | 35597 | 16426 | 19171 | 35597 | yes | 1/3 only |

The one rigorous conclusion:

* **[THEOREM U — isotropic no-go without any anisotropy]** At the named uniform
  couplings above, the layer graph `T u P_m` (`m = 0,1,2,3,4`, so `n = 5,...,9`)
  is **not** a full `n`-mode Gaussian subset-product multiset.  Only the
  certified lower bounds and Lemma 1 are used.  §10c upgrades the `t`-quantifier
  from these named points to *all but finitely many* `t`.

* **[CONDITIONAL — what would have been Theorem U2]** *If* the
  characteristic-zero core counts equal the certified values, `r_all(T) = 445`
  and `r(T) = 417`, *and* the Lemma 5 labelling map is injective at `m = 2, 4`,
  *then* Lemma 5 gives `r = 3893` and `r = 35597` exactly and Lemma 6 gives the
  excess `3^m * 202 - 2^m * (28 - 32) = 202 * 3^m + 4 * 2^m`, i.e. `1834` and
  `16426`, matching the observed values.  Note `r_all(T) - r(T) = 445 - 417 = 28`
  and **not** `2^5 = 32`, which is why the excess carries the extra `4 * 2^m`;
  the unconditional half of Lemma 6 only gives `>= 202 * 3^m`, and even that is
  conditional on the injectivity hypothesis.  Neither hypothesis is certified.

### [COMPUTATION] Where the shortfall at odd `m` comes from: a parity pattern

Disjoint unions of open chains at a *shared* coupling, `t = 1/3`.  Here both
components are Gaussian, so Lemma 1 supplies the **true** upper bound
`r <= 3^n - 2^n`; a modular lower bound meeting it therefore pins `r` exactly and
certifies injectivity of the pair-exponent map.  A shortfall certifies nothing,
except in the equal-length rows, where the two components have literally
identical mode multisets and dependence is immediate:

| components | `n` | certified `r >=` | `3^n - 2^n` | injectivity certified? |
|---|---:|---:|---:|---|
| `P_1 u P_1` | 2 | 3 | 5 | no — and here dependence is immediate (equal modes) |
| `P_1 u P_2` | 3 | 19 | 19 | **yes** (`r = 19` exactly) |
| `P_1 u P_3` | 4 | 51 | 65 | no (non-certificate) |
| `P_1 u P_4` | 5 | 211 | 211 | **yes** (`r = 211` exactly) |
| `P_1 u P_5` | 6 | 603 | 665 | no (non-certificate) |
| `P_2 u P_2` | 4 | 21 | 65 | no — dependence immediate (equal modes) |
| `P_2 u P_3` | 5 | 211 | 211 | **yes** |
| `P_2 u P_4` | 6 | 603 | 665 | no (non-certificate) |
| `P_2 u P_5` | 7 | 2059 | 2059 | **yes** |
| `P_3 u P_3` | 6 | 117 | 665 | no — dependence immediate (equal modes) |
| `P_3 u P_4` | 7 | 2059 | 2059 | **yes** |

The pattern is exact on this table: **injectivity is certified precisely when the
two vertex counts have opposite parity** (equivalently when `n` is odd).  Only
the equal-length rows are proved dependent; the unequal same-parity rows are
non-certificates.  The `T` core behaves like an odd component, which is why the
`T u P_m` count equals the reference value exactly for `m` even.

**[COMPUTATION — floating point diagnostic, not used in any proof]** The
mechanism is an alternating-product invariant.  Sorting the modes of `P_m`
ascending and putting `S(m) = sum_k (-1)^{k+1} log u_k`, the extracted modes at
`t = 1/3` satisfy, to `< 1e-9`,

```
    S(m) = 2 log( (1-t)/(1+t) )      for m odd,
    S(m) = - log q(t)                for m even,
```

i.e. `prod_k u_k^{(-1)^{k+1}}` depends only on the parity of `m`.  Two
same-parity chains therefore satisfy the multiplicative relation
`prod u^{(-1)^{k+1}} / prod v^{(-1)^{k+1}} = 1`, whose exponent vector lies in
`{-1,0,1}^n` and hence merges pair products.  This is a diagnostic: the
*rigorous* statements are the certified equalities above and the modular lower
bounds, never this float computation.

### [THEOREM F2 — uniform-field family theorem, conditional]

*Let `G` be a layer graph on `n` sites with a designated edge set `X` whose
removal leaves `T u P_{n-5}` (a spanning subgraph consisting of the 5-vertex tree
`T` and an open chain on the other `n-5` vertices).  Give every site the same
field `t`, every edge of `G \ X` the same weight `w`, and every edge of `X` the
weight `y`.  Assume both of*

> **CZ:** *the characteristic-zero core counts are the certified values,
> `r_all(T) = 445` and `r(T) = 417`, at the relevant `(t,w)`;*
>
> **INJ(n-5):** *the labelling map `(f, alpha_s alpha_{s'}) -> value` of Lemma 5
> is injective for the decoupled operator `R_T(t,w) tensor R_{P_{n-5}}(t,w)`.*

*Then there is a nonzero `rho` in `Q[t,w,y]`, of degree at most
`4(2n+|E|) S (S-1)`, such that for all `(t,w,y)` with `rho != 0`,*

```
    r  >=  3^n - 2^n + 202 * 3^{n-5} ,
```

*so `G` is not a full `n`-mode Gaussian subset-product multiset.*

*Proof.* Specialize `y -> 1`; by Lemma 4 the operator becomes
`R_T(t,w) tensor R_{P_{n-5}}(t,w)`.  INJ makes Lemma 5 an equality, CZ makes
`E_T = 445 - 3^5 = 202`, and Lemma 6 gives the excess (at least `202 * 3^{n-5}`,
in fact `202 * 3^{n-5} + 4 * 2^{n-5}` because `r_all(T) - r(T) = 28 < 32`).
Lemmas 7 and 8 transport the lower bound to the generic point of the
three-parameter family `(t,w,y)`. `[]`

**What this buys, stated exactly.**  The degeneration point has **uniform fields
and uniform surviving bonds**; only the cut bonds are special.  So the conclusion
holds on the whole three-parameter *physical* family (transverse field, in-layer
bond weight, one distinguished bond class) off a hypersurface — and in
particular, at fixed generic `(t,w)`, for **all but finitely many** cut-bond
weights `y`.  That is strictly closer to physics than Theorem F, which needed a
different parameter for every site and every bond.

**What it does not buy.**  The fully isotropic point is `y = w`.  Theorem F2 does
not decide it: `y = w` is a single *unevaluated* value of the parameter, and the
theorem asserts nothing about whether it lies in `{rho = 0}` — it may or may not.
Both hypotheses CZ and INJ are unverified: CZ would need an exact
characteristic-zero gcd degree for a degree-`496` and a degree-`528` integer
polynomial, and INJ is only *consistent with* the certified numbers at `m = 2, 4`
(§10, warning).  Theorem U is independent of both and stands unconditionally.

---

## 10b. Bounded direct isotropic attempts — including a new `3x3` certificate

The producer also runs, under a predeclared wall, the *direct* isotropic
pair-product computation at `t = 1/3` for the `2x4` and `3x3` layers, to extend
the per-size evidence of (R2).  Both completed:

| isotropic layer, `t = 1/3` | dim | pair slots | `deg gcd(C_2,C_2')` | `r >=` | `3^n - 2^n` | verdict |
|---|---:|---:|---:|---:|---:|---|
| `2x3` | 64 | 2016 | 385 | 1631 | 665 | not Gaussian |
| `2x4` | 256 | 32640 | 9329 | 23311 | 6305 | not Gaussian |
| `3x3` | 512 | 130816 | 59641 | **71175** | 19171 | **not Gaussian** |

**[THEOREM 3x3 — new].**  The open `3x3` Ising layer at `t = tanh(K*/2) = 1/3`,
`e^{2K} = 5/3`, isotropic, is **not** a full nine-mode Gaussian subset-product
multiset.  `proofs/pair_product_scale.md` §6 records a 90-second resource wall
at exactly this computation and explicitly claims no `3x3` conclusion; the
optimized pipeline here (float64-exact BLAS modular products, contiguous
reversed Newton buffer, chunked exact dot products) completes it.  This is the
first layer in the repository with a degree-4 vertex to be certified.

The `2x3` row independently reproduces `proofs/pair_product_obstruction.md`
(`deg gcd = 385`) and the `2x4` row reproduces `proofs/pair_product_scale.md`
(`deg gcd = 9329`), both from a different construction route; the `385`
reproduction is kept as a regression assertion in the standalone test.  A wall,
had one fired, would have been reported as a wall and never as a result.

---

## 10c. Isotropic-curve genericity: from named points to all but finitely many `t`

### [LEMMA 10 — the isotropic curve is a polynomial family]

*Fix `G` and restrict to the isotropic curve `w_e = q(t) = (1+t^2)/(2t)`.  Then*

```
    M(t) := (2t)^{|E|} P_t D(q(t)) P_t   in   Z[t]^{2^n x 2^n},
```

*with every entry of `t`-degree at most `d_R := 2n + 2|E|`, and `M(t)` is a
nonzero scalar multiple of the layer operator at each `t != 0`.*

*Proof.* `(2t)^{|E|} q^{a} = (2t)^{|E|-a}(1+t^2)^a` is an integer polynomial of
degree `|E| + a <= 2|E|` for `0 <= a <= |E|`, and the two factors `P_t` contribute
`t`-degree at most `n` each. `[]`

### [THEOREM G — genericity along the isotropic curve]

*Let `G` be a layer graph for which a certified lower bound
`r(t_0) >= r_0 > 3^n - 2^n` is available at one rational `t_0` on the isotropic
curve.  Then there is a nonzero `rho in Q[t]` with*

```
    deg rho  <=  4 (2n + 2|E|) * S (S - 1),        S = C(2^n, 2),
```

*such that for **every** `t` with `rho(t) != 0` the isotropic layer at coupling
`t` has `r(t) >= r_0 > 3^n - 2^n`, hence is not a full `n`-mode Gaussian
subset-product multiset.*

*Proof.* Lemma 10 puts the family in `Z[t]`.  Lemma 7 applied to the
specialization `t -> t_0` gives `g <= S - r_0` for the generic gcd degree `g` over
`Q(t)`, so `r_gen = S - g >= r_0`.  Lemma 8, with `d_R = 2n + 2|E|` and `B_0`
normalized as in the corrected statement, produces `rho` and the equality
`r(t) = r_gen` off its roots. `[]`

**[COMPUTATION]** Applying Theorem G to every isotropic certificate in this note
upgrades all of them from single named couplings to all-but-finitely-many `t`:

| isotropic graph | `n` | `|E|` | certified `r >=` at `t = 1/3` | `3^n - 2^n` | `deg rho <=` |
|---|---:|---:|---:|---:|---:|
| `T` | 5 | 4 | 417 | 211 | 17,677,440 |
| `T u P_1` | 6 | 4 | 1219 | 665 | 324,979,200 |
| `T u P_2` | 7 | 5 | 3893 | 2059 | 6,341,400,576 |
| `T u P_3` | 8 | 6 | 7885 | 6305 | 119,317,739,520 |
| `T u P_4` | 9 | 7 | 35597 | 19171 | 2,190,424,965,120 |
| `2x3` | 6 | 7 | 1631 | 665 | 422,472,960 |
| `2x4` | 8 | 10 | 23311 | 6305 | 153,408,522,240 |
| `3x3` | 9 | 12 | 71175 | 19171 | 2,874,932,766,720 |

The `2x3` row reproduces `proofs/pair_product_scale.md`'s exceptional-set bound
`422,472,960` exactly, which is an independent check of the degree bookkeeping
(their `52` is our `2 d_R = 2(2 * 6 + 2 * 7)`).  So **Theorem U holds for all but
finitely many `t`** on the isotropic curve, for all five graphs `T u P_m`, and the
same upgrade applies to `2x3`, `2x4` and the new `3x3` certificate.  The
exceptional sets are bounded but not located.

### Resource budgets (observed, not preflight estimates)

Predeclared walls: `900 s` per ordinary case, `1500 s` for the four large cases
(`n = 9` decoupled specialization, uniform `T u P_4`, isotropic `2x4`, isotropic
`3x3`), RSS cap `7,500,000,000` bytes.  **No wall fired in the recorded run.**
Observed: 98 checks, 0 failures, `899.8 s` total wall, peak process RSS
`162,004,992` bytes (2.2% of the cap).  The slowest single case is the isotropic
`3x3` gcd.  The standalone test, which independently rebuilds the isotropic
`2x4`, isotropic `3x3` and uniform `T u P_4` certificates including their
coefficient digests, runs `53` assertions in `155 s`.  All timings are in the
artifact.

---

## 11. What is **not** established

* **[UNRESOLVED]** Theorem F does not cover the isotropic curve `t_v = t`,
  `w_e = q(t)` for any `n` beyond the graphs with an explicit certificate
  (`2x3`, `2x4`, `3x3` and `T u P_m` for `m <= 4`, each now for all but finitely
  many `t` via §10c).  §9 shows this is a limitation of localization and
  §9.4/§10 state the repairs.
* **[UNRESOLVED — the largest single gap]** Every core count in this note
  (`113`, `129` for the claw; `417`, `445` for `T`) is a proved *lower* bound
  only: Lemma 3 is one-sided and agreement at two primes is an implementation
  cross-check, not a proof.  Consequently the Lemma 5' *upper* bound cannot be
  evaluated, so no composite count in §10 is pinned, and the injectivity
  hypothesis of Lemma 5 is nowhere certified for a shared coupling.  Pinning the
  characteristic-zero values would need an exact `Q`-gcd degree for integer
  polynomials of degree `120`/`136` (claw) and `496`/`528` (`T`); we did not
  attempt it.  A larger true value would only strengthen Theorem F's bound (2).
* **[UNRESOLVED]** At `m = 1, 3` the certified count is below the lower-bound
  tensor reference (`1219` vs `1279`, `7885` vs `11791`).  This is a
  non-certificate in both directions: it neither pins `r` nor proves that any
  injectivity or independence hypothesis fails.
* **[UNRESOLVED]** The parity closed form
  `prod_k u_k^{(-1)^{k+1}} = ((1-t)/(1+t))^2` (odd `m`), `= 1/q` (even `m`) for
  the open chain is a floating-point diagnostic verified to residual `< 3e-15`
  for `m = 1..7`; it is not proved and no statement depends on it.  Among the
  disjoint-chain rows only the equal-length ones (`P_k u P_k`) are proved
  dependent, and that for the trivial reason of identical mode multisets.
* **[UNRESOLVED]** All exceptional sets (`{rho_G = 0}` in Theorem F,
  `{rho = 0}` in Theorem G) are bounded in degree but neither located nor
  evaluated.  No statement here is an interval or neighbourhood statement, and in
  particular Theorem F2 asserts nothing about the physical value `y = w`.
* **[UNRESOLVED]** "Not certified flagged" is never an assertion that a layer has
  no excess.  The isotropic claw, `K_{1,4}`, `C_5` and the isotropic rungs-off
  `2x3` are all non-certificates; only where a Gaussian upper bound (Lemma 1) or
  Lemma 9's proved bound applies is a count pinned or an impossibility proved.
* **[UNRESOLVED]** Nothing here decides parity-projected sectors, restrictions of
  Gaussian spectra on larger Hilbert spaces, or interacting-fermion
  representations; those are different quantifiers (cf.
  `proofs/parity_pair_product.md`).
* **[UNRESOLVED]** Cycles.  A layer graph of maximum degree 2 containing a cycle
  is free-fermionic only sector-by-sector, so it is neither covered by Theorem F
  (no branching site) nor by Theorem C (not a union of paths).  The `C_5`
  computation in §9.3 (`r = 141 < 211`) is consistent with the sector picture but
  proves nothing about it.
* **[UNRESOLVED]** No thermodynamic-limit statement is made.  Theorem F is
  uniform in `n` but is a statement about each finite layer.
* **[EXTERNAL]** The Bogoliubov step of Theorem C is standard free-fermion
  theory (Lieb–Schultz–Mattis 1961; Schultz–Mattis–Lieb 1964) and is used as
  cited, not re-derived.  It is only used for the *control*, never for
  Theorem F.
