# Isotropic all-size spectral route: surviving cores and an exact decoupling obstruction

Artifacts: `results/spectral/isotropic_allsize.json` (integrated certificate),
`results/spectral/isotropic_core_search.json`, and
`results/spectral/isotropic_decoupling.json`.  Producers:
`experiments/e163_isotropic_core_search.py`,
`experiments/e164_isotropic_decoupling.py`, and
`experiments/e165_isotropic_allsize.py`.  Independent verifier:
`tests/test_isotropic_allsize.py`.

## 1. Outcome and exact scope

**[UNRESOLVED — requested all-size spectral theorem].**  This run does **not**
prove that every finite graph with a vertex of degree at least three is
spectrally non-Gaussian at all but finitely many isotropic couplings.  Therefore
there is no proved `n_0` to report.  The statement remains open even though the
finite-core half of the proposed repair succeeds.

**[COMPUTATION — the piece that survives].**  At the uniform physical point
`t=1/3`, `q(t)=(1+t^2)/(2t)=5/3`, three small cores have exact modular
pair-polynomial derivative-gcd degrees, identical at `p=1000003` and
`p=2000003`, which give characteristic-zero pair-product lower bounds strictly
above their own Gaussian ceilings.  They are the paw, the five-site tree `T`
with branch lengths `2,1,1`, and the open `2x2` grid with one pendant vertex.

**[THEOREM — the piece that fails].**  The graph-piece tensor localization of
Lemma 5 in `proofs/allsize_gaussian.md` cannot be made into an identity on the
isotropic curve across a nonempty edge cut: for a cut containing `c>=1` uniform
edges, a necessary rank-one minor is `q(t)^(2c)-1`, which is a nonzero rational
function.  Moreover, the valid all-size branching family

```
G_m = K_{1,3} disjoint-union m K_1
```

has only identical free-site modes at uniform field.  For every `m>=3`, its
number of distinct unordered-slot pair products is at most
`136(2m+1)`, strictly below the Gaussian ceiling
`3^(m+4)-2^(m+4)`.  Thus the proposed Lemma-5 pair-product localization route
cannot prove the requested every-graph theorem.

**[SCOPE — part of the preceding theorem].**  This is an exact obstruction to
that proof method, not a proof that `G_m` is Gaussian, and not a disproof of an
isotropic all-size no-go obtainable from another invariant or a non-localizing
argument.  The obstruction family starts at `m=3`, hence at `n=7`; this `7` is
**not** an `n_0` for a spectral no-go.

## 2. Imported exact lemmas and direction discipline

**[LEMMA — Gaussian ceiling, imported].**  By Lemma 1 of
`proofs/allsize_gaussian.md`, a full `n`-mode Gaussian subset-product multiset
has at most

```
3^n - 2^n
```

distinct products over unordered pairs of distinct eigenvalue slots.

**[LEMMA — monic modular direction, imported].**  For the actual-denominator
integer model `M=D R`, the exterior-pair characteristic polynomial

```
A_M(z) = charpoly_z(wedge^2 M)
```

is monic in `Z[z]`.  Consequently, for every prime used here,

```
deg gcd_Q(A_M,A_M') <= deg gcd_Fp(A_M mod p, A_M' mod p),
```

so a modular degree `g_p` proves only the one-sided statement

```
r_Q >= C(2^n,2) - g_p.
```

This is Lemma 2 of `proofs/pair_product_scale.md`; agreement at two primes is an
independent arithmetic cross-check and is never promoted to equality in
characteristic zero.

**[LEMMA — monic-in-parameter upgrade, imported].**  Along the isotropic curve,
Lemma 10 and Theorem G of `proofs/allsize_gaussian.md` use

```
M_H(t) = (2t)^|E(H)| P_t D_H(q(t)) P_t in Z[t]^(2^k x 2^k)
```

and the monicity of `charpoly(wedge^2 M_H(t))` in `z`.  One rational point with
`r(t_0)>=r_0` bounds the generic derivative-gcd degree in the helpful direction;
a nonzero resultant then confines exceptional couplings to the roots of a
nonzero polynomial of degree at most

```
4 (2k + 2|E(H)|) S(S-1),       S=C(2^k,2).
```

## 3. Exact isotropic core search

**[COMPUTATION — raw certificate table].**  Every row below uses the uniform
point `t=1/3`, `q=5/3`.  The producer constructs the actual-LCM monic integer
model and reconstructs `charpoly(wedge^2 M)` by exact modular power traces,
Newton identities, the Cayley--Hamilton recurrence, and Euclidean gcd.  The
standalone verifier independently rebuilds `R=PDP` by rational triple loops and
repeats both prime computations.

| graph | `k` | edges | pair slots | gcd degree at `p1` | gcd degree at `p2` | proved `r_Q >=` | `3^k-2^k` | result |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| claw `K_{1,3}` | 4 | 3 | 120 | 55 | 55 | 65 | 65 | non-certificate |
| `T`: claw with one branch extended | 5 | 4 | 496 | 79 | 79 | **417** | 211 | certified |
| paw: triangle with one pendant | 4 | 4 | 120 | 1 | 1 | **119** | 65 | certified |
| open `2x2` grid plus pendant | 5 | 5 | 496 | 21 | 21 | **475** | 211 | certified |
| star `K_{1,4}` | 5 | 4 | 496 | 349 | 349 | 147 | 211 | non-certificate |

**[THEOREM — three named isotropic core families].**  For each of the following
three named finite graphs, there is a nonzero polynomial `rho_H(t)` such that,
whenever `rho_H(t)!=0`, the isotropic layer spectrum is not a full `k`-mode
Gaussian subset-product multiset:

| graph `H` | proved generic `r >=` | Gaussian ceiling | `deg rho_H <=` |
|---|---:|---:|---:|
| paw | 119 | 65 | 913,920 |
| five-site tree `T` | 417 | 211 | 17,677,440 |
| open `2x2` grid plus pendant | 475 | 211 | 19,641,600 |

**[THEOREM — proof].**  The modular direction lemma gives the displayed
characteristic-zero lower bound at `t_0=1/3`; each lower bound is strictly above
`3^k-2^k`.  The Gaussian ceiling excludes a full subset-product spectrum at that
point.  The monic-in-parameter lemma upgrades the same lower bound to every
`t` outside the roots of a nonzero `rho_H`, and substitution into
`4(2k+2|E|)S(S-1)` gives the three displayed degree bounds. `[]`

**[SCOPE — finite cores only].**  No monotonicity of pair-product counts under
adding edges or vertices is known or used.  In particular, these three named
core theorems do not imply the same conclusion for a graph merely containing
one of the cores as a subgraph.

## 4. Why an isotropic graph-piece decoupling identity fails

**[LEMMA — nonzero cut minor].**  Partition the sites as `A disjoint-union B`
and suppose `c>=1` graph edges cross the cut.  For `t^2!=1`, `P_t` is invertible.
A Lemma-5 graph-piece factorization across this site partition would therefore
force the diagonal Boltzmann weight to be a product of a function of the spins
in `A` and a function of the spins in `B`.

**[LEMMA — proof].**  Restrict the crossing Boltzmann factor to the four
configurations in which all spins in each part are separately constant.  After
removing the internal-edge row and column factors, its `2x2` matrix is

```
[ q^c   1  ]
[  1   q^c ]
```

and has determinant `q^(2c)-1`.  On substituting
`q(t)=(1+t^2)/(2t)` and clearing denominators, its numerator is

```
(1+t^2)^(2c) - (2t)^(2c).
```

This polynomial has constant coefficient `1` and leading coefficient `1`, so it
is nonzero for every integer `c>=1`.  Hence the required rank-one factorization
is not an identity in `Q(t)`. `[]`

**[THEOREM — physical interval].**  For `0<t<1`,

```
q(t)-1 = (t-1)^2/(2t) > 0,
```

so `q(t)^(2c)-1>0` for every `c>=1`; graph-piece tensor decoupling across a
nonempty cut is impossible at every such physical isotropic coupling.

**[LEMMA — the only positive endpoint collapses].**  The positive solution of
`q(t)=1` is `t=1`.  At `t=1`, every one-site factor
`[[1,t],[t,1]]` has rank one, hence every `P_t` and every layer operator
`P_t D P_t` has rank at most one.  Its unordered distinct-slot pair products
therefore have only the value zero.  The endpoint cannot supply distinct free
modes or a pair-product seed.

**[SCOPE — what the cut lemma excludes].**  The lemma excludes precisely the
site-partition graph-component identity required by Lemma 5.  It does not claim
that an arbitrary matrix could not admit an unrelated tensor factorization
after a nonlocal change of basis.

## 5. An all-size exact counterfamily to the localization method

**[LEMMA — exact uniform tensor identity].**  Let
`G_m=K_{1,3} disjoint-union m K_1`.  It is a finite simple graph on `m+4` sites
with maximum degree three.  At uniform parameters its operator factorizes
exactly as

```
R_Gm(t,q(t)) = R_K1,3(t,q(t)) tensor S(t)^(tensor m),
S(t) = [[1,t],[t,1]]^2.
```

The two eigenvalues of `S(t)` are `(1+t)^2` and `(1-t)^2`.

**[LEMMA — identical-mode count].**  The free eigenvalue slots are indexed by
`j=0,...,m` and have values proportional to
`(1+t)^(2(m-j))(1-t)^(2j)`.  A product of two such slots depends on the free
indices only through `j+j' in {0,...,2m}`, so the free part contributes at most
`2m+1` pair-product classes.  The 16-slot claw has at most
`C(16+1,2)=136` pair products when repetitions are allowed.  Therefore

```
r(G_m) <= U_m := 136(2m+1).
```

**[THEOREM — strict collapse for every `m>=3`].**  Put
`C_m=3^(m+4)-2^(m+4)`.  At `m=3`,

```
U_3 = 952 < 2059 = C_3.
```

For every `m>=3`, the exact identities

```
2U_m - U_(m+1) = 136(2m-1) > 0,
C_(m+1) - 2C_m = 81*3^m > 0
```

propagate `U_m<C_m` by induction.  Thus the pair-product invariant is provably
below its Gaussian ceiling on every `G_m` with `m>=3`.

**[THEOREM — consequence for the requested repair].**  The remainder of `G_m`
after selecting its claw consists only of isolated vertices; edge deletion
cannot create an open chain.  Uniform fields make all of their one-site modes
identical, and the preceding theorem proves the resulting linear collapse.
Thus replacing free sites by a non-isomorphic or opposite-parity open chain is
not available uniformly over the requested class of **every** branching graph.
Together with the nonzero-cut-minor lemma, this exactly blocks both ways of
obtaining the required isotropic Lemma-5 identity.

**[SCOPE — no converse].**  The inequality `r(G_m)<3^(m+4)-2^(m+4)` is only an
upper bound on this obstruction invariant.  It does not certify a Gaussian
representation of the spectrum.

## 6. Final ledger statement

**[COMPUTATION — two-piece verdict].**  Piece 1 survives: three uniform cores
have decisive exact modular certificates, and Theorem G makes each named core
non-Gaussian for all but finitely many isotropic couplings.  Piece 2 fails: an
isotropic graph-piece decoupling across a nonempty cut has a nonzero separating
minor, while the valid no-cut family `G_m` has identical modes and an exact
linear pair-count ceiling.

**[UNRESOLVED — remaining route].**  Closing the isotropic all-size spectral
gap now requires a non-localizing invariant or a direct isotropic certificate
for every graph/size.  Neither the finite core table nor generic anisotropic
Theorem F supplies that missing quantifier.

## 7. Reproduction

**[COMPUTATION — producer].**  From the repository root:

```bash
.venv/bin/python experiments/e165_isotropic_allsize.py
```

The command writes `results/spectral/isotropic_allsize.json` and prints a final
`PASS e165 isotropic all-size obstruction` only when every stored check passes.

**[COMPUTATION — independent verifier].**  From the repository root:

```bash
.venv/bin/python tests/test_isotropic_allsize.py
```

The verifier imports none of `e163`, `e164`, or `e165`; it rebuilds every raw
operator and both modular gcd degrees, checks every finite audit row, and checks
the closed-form base and induction identities for the universal `m>=3` claim.
Its final line is `PASS test_isotropic_allsize (58 checks)`.

**[COMPUTATION — component producers].**  The component artifacts can be
regenerated independently with

```bash
.venv/bin/python experiments/e163_isotropic_core_search.py
.venv/bin/python experiments/e164_isotropic_decoupling.py
```
