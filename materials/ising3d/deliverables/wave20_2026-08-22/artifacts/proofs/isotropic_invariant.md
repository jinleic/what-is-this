# A direct isotropic trace invariant: odd-cycle no-go at every physical coupling

Artifacts: `experiments/e200_isotropic_puiseux.py`,
`experiments/e201_isotropic_invariant.py`,
`experiments/e202_isotropic_certificate.py`,
`results/spectral/isotropic_invariant.json`, and the independent verifier
`tests/test_isotropic_invariant.py`.

## 1. Outcome and exact scope

**[THEOREM — all-size, every-coupling result].** Let `G` be any finite simple
**non-bipartite** graph on `n` vertices. For every real physical isotropic
parameter `0<t<1`, the `2^n` eigenvalue slots of its Ising layer transfer
operator are not a full `n`-mode subset-product multiset, even with arbitrary
complex nonzero subset parameters. The exceptional set inside the physical
interval is empty.

**[THEOREM — new connected branching family].** In particular, let `L_r` be a
triangle with a path of `r>=1` edges attached at one triangle vertex. Then
`L_r` is connected, has a degree-three vertex, has `n=r+3`, and is spectrally
non-Gaussian at every `0<t<1`. Thus one connected branching example is proved
at every size `n>=4` by one exact recurrence, not by a finite census.

**[THEOREM — method boundary].** The same invariant vanishes identically for
every bipartite graph at every coupling. Hence this route exactly separates
bipartite from non-bipartite graphs, but it does not decide open rectangular
grids or the requested theorem for **every** branching graph.

**[UNRESOLVED].** No conclusion is obtained here for bipartite branching
layers, no thermodynamic-limit identity is asserted, and no part of this note
solves the three-dimensional Ising model. The path and square controls are
counterexamples to universal efficacy of this trace invariant, not
counterexamples to the non-bipartite theorem.

**[THEOREM — route separation].** The proof works directly with a
characteristic-coefficient trace invariant on
`(y-1)x^2-(y+1)=0`. It neither switches off graph edges to tensor-localize a
finite core nor differentiates a family of varying conjugators.

## 2. A determinant-center-free subset-product invariant

**[LEMMA — invariant definition].** Let `W` be an invertible `N x N` matrix,
where `N=2^n`, and put `M=N/2`. Write `e_k(W)` for the `k`-th elementary
symmetric polynomial of its eigenvalue slots, including algebraic
multiplicity. Define

```text
Phi_N(W) = e_N(W)^(M-1) e_1(W)^M - e_(N-1)(W)^M.             (1)
```

This is a conjugacy-invariant polynomial in the characteristic coefficients.
It is homogeneous under nonzero scalar multiplication of `W`, so the equation
`Phi_N(W)=0` is unaffected by the scalar normalization of a transfer matrix.

**[LEMMA — homogeneity proof].** Under `W -> sW`, the first term in (1) has
`s`-degree

```text
N(M-1)+M = M(N-1),
```

where the equality uses `N=2M`; the second term has the same degree because
`e_(N-1)` has degree `N-1`. `[]`

**[LEMMA — subset-product necessity].** Every invertible full `n`-mode
subset-product multiset

```text
lambda_S = a product_{i in S} u_i,       S subseteq {1,...,n},
a,u_i in C^times,
```

obeys `Phi_N=0`.

**[LEMMA — proof].** Put `c=a^2 product_i u_i`. Complementary slots obey
`lambda_S lambda_(S^c)=c`. There are `M` complementary pairs, so

```text
e_N = product_S lambda_S = c^M.
```

Complementation also permutes all slots, giving

```text
e_1 = sum_S lambda_S = c sum_S lambda_S^(-1),
e_(N-1) = e_N sum_S lambda_S^(-1).
```

Raising the last two identities to the `M`-th power gives
`e_N^(M-1)e_1^M=e_(N-1)^M`, which is (1). No ordering, distinctness, or
diagonalizability assumption is used. `[]`

**[LEMMA — determinant-one form].** If `det W=1`, then
`e_(N-1)(W)=tr(W^(-1))`, so

```text
Phi_N(W) = tr(W)^M - tr(W^(-1))^M.                         (2)
```

**[COMPUTATION — nontriviality controls].** The producer constructs exact
rational subset-product slots for `n=1,2,3,4` and obtains `Phi_N=0` in every
case. For the positive determinant-one four-slot multiset
`{2,3,1/5,5/6}`, it obtains

```text
e_1=181/30,       e_3=211/30,       Phi_4=-196/15 != 0.
```

These are algebra audits only; Lemma 2 is the all-`n` proof.

## 3. Exact reduction to a graph polynomial

**[LEMMA — physical determinant-one representative].** Put

```text
A_G = sum_{v in V} X_v,
B_G = sum_{uv in E} Z_u Z_v,
V_G(a,b) = exp(a A_G) exp(b B_G).
```

Every Pauli word in `A_G` and `B_G` is traceless, hence `det V_G=1`.
Furthermore `V_G` is similar to the symmetric physical sandwich

```text
exp(a A_G/2) exp(b B_G) exp(a A_G/2),
```

so it has the physical layer spectrum. A different conventional nonzero
scalar normalization does not change the zero or nonzero status of (1).

**[LEMMA — alignment polynomials].** Let `m=|E|`, let `a_G(sigma)` be the
number of aligned edges in a spin configuration `sigma in {+-1}^V`, and define

```text
Q_G(y)   = sum_sigma y^a_G(sigma),
Q_G^*(y) = y^m Q_G(1/y)
         = sum_sigma y^(m-a_G(sigma)),
I_G(y)   = Q_G(y)-Q_G^*(y).                                  (3)
```

All three are integer polynomials of degree at most `m`.

**[LEMMA — trace reduction].** In multiplicative coordinates
`x=exp(a)`, `y=exp(2b)`,

```text
tr V_G      = cosh(a)^n y^(-m/2) Q_G(y),
tr V_G^(-1) = cosh(a)^n y^(-m/2) Q_G^*(y).                    (4)
```

Consequently, with `M=2^(n-1)`, (2) is a nonzero common factor times

```text
Psi_G(y) = Q_G(y)^M - Q_G^*(y)^M
         = I_G(y) sum_{j=0}^{M-1} Q_G(y)^(M-1-j) Q_G^*(y)^j.  (5)
```

**[LEMMA — proof].** In the `Z` basis the diagonal entries of
`exp(aA_G)=tensor_v exp(aX_v)` all equal `cosh(a)^n`, while `exp(bB_G)` is
diagonal. Its spin weight is
`exp(b(2a_G(sigma)-m))=y^(a_G(sigma)-m/2)`, proving the first formula.
Cyclicity gives
`tr V_G^(-1)=tr(exp(-aA_G)exp(-bB_G))`; `cosh` is even and `b -> -b`
replaces aligned by disagreeing exponents, proving the second. Equation (5)
is the difference-of-powers factorization. `[]`

## 4. The high-temperature leading term

**[LEMMA — Eulerian expansion].** Let `C(G)` be the set of edge subsets in
which every vertex has even selected degree. Then

```text
I_G(y) = 2^(n+1-m)
         sum_{F in C(G), |F| odd}
             (y-1)^|F| (y+1)^(m-|F|).                         (6)
```

**[LEMMA — proof].** Put `u=tanh b=(y-1)/(y+1)`. Expanding each edge factor
as

```text
exp(b sigma_u sigma_v)=cosh(b)(1+u sigma_u sigma_v)
```

and summing over spins kills exactly the edge subsets with a nonempty mod-two
boundary. Therefore

```text
Z_G(b)=2^n cosh(b)^m sum_{F in C(G)} u^|F|.
```

Subtracting `Z_G(-b)` retains twice the odd-cardinality terms. Finally,
`I_G(y)=y^(m/2)(Z_G(b)-Z_G(-b))` and
`y^(1/2)cosh(b)=(y+1)/2`, giving (6). `[]`

**[THEOREM — exact graph characterization].** `I_G` is the zero polynomial if
and only if `G` is bipartite.

**[THEOREM — proof].** Every Eulerian edge set decomposes into edge-disjoint
cycles. If `G` is bipartite, every cycle is even, so (6) has no term. If `G`
is non-bipartite, an odd cycle is itself an odd-cardinality Eulerian edge set;
all coefficients in (6) have the same positive sign for real `y>1`, so the sum
is nonzero. Equivalently, the usual two-colouring spin flip bijects aligned and
disagreeing edges exactly in the bipartite case. `[]`

**[LEMMA — leading Taylor/Puiseux germ].** Suppose `G` has odd girth `g`, and
let `c_g` be its number of shortest odd cycles. With `s=y-1`,

```text
I_G(1+s) = 2^(n+1-g) c_g s^g + O(s^(g+1)).                    (7)
```

**[LEMMA — proof].** An odd Eulerian edge set contains an odd cycle, so it has
at least `g` edges. At cardinality `g` it can only be one shortest odd cycle.
Taking the lowest `s`-degree in (6) gives
`2^(n+1-m)c_g 2^(m-g)=2^(n+1-g)c_g`. `[]`

**[THEOREM — meaning of the Puiseux route].** Individual eigenvalues may split
into genuine Puiseux branches at the endpoint, but the symmetric invariant
already lies in `Z[y]`; its first nonzero Puiseux exponent is therefore the
integer exponent `g` in (7). No matching or ordering of eigenvalue branches is
needed.

## 5. Nonzero modulo the physical curve and finite exceptions

**[LEMMA — physical curve].** With

```text
t=tanh(a/2),
x=(1+t)/(1-t),
y=q(t)=(1+t^2)/(2t),
```

the multiplicative coordinates satisfy

```text
F_phys(x,y)=(y-1)x^2-(y+1)=0.                                (8)
```

**[LEMMA — exact nonzero point for every non-bipartite graph].** At `t=1/3`,
`(x,y)=(2,5/3)` lies on (8). For a non-bipartite `G`, (6) gives

```text
I_G(5/3)
 = 2^(n+1-m)(8/3)^m
   sum_{F Eulerian, |F| odd} (1/4)^|F|  > 0.                  (9)
```

Since `Q_G(5/3)>Q_G^*(5/3)>0`, equation (5) gives
`Psi_G(5/3)>0`.

**[THEOREM — nonzero modulo the curve].** For every non-bipartite `G`, the
class of `Psi_G(y)` modulo `F_phys(x,y)` is nonzero. Indeed, the point in (9)
is on the curve and evaluates `Psi_G` to a nonzero exact rational. Thus this is
curve nonvanishing, not ambient two-parameter noncontainment.

**[THEOREM — finite complex exceptional set].** Fix a non-bipartite `G` with
`m` edges and put `M=2^(n-1)`. Then

```text
K_G(t)=(2t)^(mM) Psi_G((1+t^2)/(2t)) in Z[t]                 (10)
```

is nonzero and has degree exactly `2mM`. On the nonsingular multiplicative part
of the complex physical curve, every subset-product point is a root of `K_G`.
Hence there are at most `2mM` complex exceptional parameters, counting roots
without claiming that every root is Gaussian.

**[THEOREM — proof].** In a non-bipartite graph no spin assignment disagrees
on every edge, whereas assignments constant on each connected component align
every edge. Thus `Q_G` has nonzero `y^m` coefficient and `Q_G^*` does not;
`Psi_G` has degree exactly `mM`. Clearing the denominator after `y=q(t)` gives
(10); its `t^(2mM)` coefficient comes only from the leading `y^(mM)` term and
is nonzero. Away from the finitely many parameter points where the
multiplicative coordinates degenerate, the common factor in (4) is nonzero,
so subset-product necessity `Phi_N=0` forces `Psi_G=0`, hence `K_G=0`. `[]`

**[THEOREM — quantifier discipline].** Equations (9) and (10) prove two
different statements: (9) proves nonzero at one exact point and therefore
nonzero modulo the curve; (10) then confines possible complex exceptions to a
finite root set. Neither statement alone says that there are no exceptional
real physical couplings.

## 6. Every real physical coupling

**[THEOREM — no physical exceptions].** If `0<t<1`, then

```text
q(t)-1=(t-1)^2/(2t)>0.
```

For non-bipartite `G`, every summand in (6) is positive at `y=q(t)>1`, so
`I_G(q(t))>0`. Also `Q_G(q(t))` and `Q_G^*(q(t))` are positive; hence (5) gives
`Psi_G(q(t))>0`, and (4) gives `Phi_N(V_G)>0`. This contradicts the necessary
subset-product identity (1) at every `0<t<1`. `[]`

**[THEOREM — independence from numerical benchmarks].** The proof uses only
the symbolic inequality `q(t)>1` and exact polynomial identities. No critical
coupling, numerical fit, or finite coupling grid selects any part of the
argument; `t=1/3` is used only as the exact algebraic witness in (9).

## 7. Graph-extension recurrence and an infinite family

**[LEMMA — one-leaf recurrence].** If `G^+` is obtained from `G` by attaching a
new leaf with one new edge, summing over the leaf spin gives

```text
Q_(G^+)(y)   =(y+1)Q_G(y),
Q_(G^+)^*(y) =(y+1)Q_G^*(y),
I_(G^+)(y)   =(y+1)I_G(y).                                  (11)
```

If `M=2^(n-1)`, the new half-dimension is `2M`, and

```text
Psi_(G^+)(y)
 =(y+1)^(2M) Psi_G(y) (Q_G(y)^M+Q_G^*(y)^M).                 (12)
```

Thus nonvanishing propagates under every leaf attachment in the integral domain
`Z[y]`, without a tensor-product spectral localization.

**[LEMMA — proof].** For each old spin configuration, the two choices of the
new spin contribute one aligned and one disagreeing edge, hence the factor
`y+1` in (11). Equation (12) follows from
`A^(2M)-B^(2M)=(A^M-B^M)(A^M+B^M)`. `[]`

**[THEOREM — explicit lollipop family].** Starting with the triangle and
iterating (11), the triangle-with-tail graph `L_r` satisfies

```text
I_(L_r)(y)=2(y-1)^3(y+1)^r,
(2t)^(r+3) I_(L_r)(q(t))=2(t-1)^6(t+1)^(2r).                 (13)
```

For `r>=1`, `L_r` is connected and branching. Formula (13), or the general
non-bipartite theorem, excludes a subset-product spectrum at every physical
coupling and every size `n=r+3>=4`.

**[THEOREM — finite-exception versus every-coupling scope for the family].**
The second polynomial in (13) proves a particularly simple nonzero trace-skew
germ, while the determinant-center-free `Psi_(L_r)` supplies the complex
finite-exception theorem. Extra complex roots of `Psi` are not claimed absent;
positivity is what removes every root on `0<t<1`.

## 8. Exact finite audits and route counterexamples

**[COMPUTATION — lemma-test table].** The producer independently constructs
`I_G` from all spins and from all Eulerian edge subsets for the following small
graphs. The two coefficient arrays agree in every row.

| graph | bipartite | odd girth | shortest odd cycles | leading coefficient in `I(1+s)` | `I(5/3)` | `deg Psi` | complex root bound |
|---|---|---:|---:|---:|---:|---:|---:|
| path `P_4` | yes | — | 0 | 0 | 0 | — | — |
| square `C_4` | yes | — | 0 | 0 | 0 | — | — |
| triangle | no | 3 | 1 | 2 | `16/27` | 12 | 24 |
| paw | no | 3 | 1 | 4 | `128/81` | 32 | 64 |
| triangle with two-edge tail | no | 3 | 1 | 8 | `1024/243` | 80 | 160 |
| `C_5` | no | 5 | 1 | 2 | `64/243` | 80 | 160 |
| `K_4` | no | 3 | 4 | 16 | `8192/729` | 48 | 96 |

**[COMPUTATION — recurrence tests].** Exact coefficient arithmetic checks
triangle `->` paw `->` two-edge tail under (11) and checks (13) for
`r=0,1,2,3,4,5`. These rows test the implementation of the uniform recurrence;
they do not replace its one-line all-`r` proof.

**[THEOREM — exact blind locus].** If `G` is bipartite, flipping every spin in
one colour class is a bijection taking `a_G(sigma)` to `m-a_G(sigma)`.
Therefore `Q_G=Q_G^*`, `I_G=0`, and `Psi_G=0` identically. The path and square
rows are exact counterexamples to any attempt to promote this invariant to all
branching or all graph classes. Every open rectangular grid is bipartite, so a
different invariant is still required there.

**[UNRESOLVED].** Vanishing of `Phi_N` is only a necessary subset-product
condition; it does not prove that a bipartite layer is Gaussian. In particular,
the blind-locus theorem is a limitation theorem for this direct trace route,
not a positive construction of fermionic modes.

## 9. Reproduction

**[COMPUTATION — producer].** From the repository root:

```bash
.venv/bin/python experiments/e202_isotropic_certificate.py
```

The command writes `results/spectral/isotropic_invariant.json` and prints
`PASS e202 isotropic invariant certificate` only if every stored exact check
passes.

**[COMPUTATION — independent verifier].** From the repository root:

```bash
.venv/bin/python tests/test_isotropic_invariant.py
```

The verifier imports none of `e200`, `e201`, or `e202`. It rebuilds the spin and
Eulerian enumerators using different iteration schemes, reconstructs every
stored graph polynomial and digest, checks the subset-product invariant on a
separate rational parameter family, and rechecks the leaf recurrence. The
finite audits remain lemma tests only.
