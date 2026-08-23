# All-L saturation and classification of the open `2 x L` layer algebra

**[THEOREM]** For every integer `L >= 2`, with `n = 2L` sites and
`m = 2^(n-1)`, the characteristic-zero local-term Lie algebra of the open
`2 x L` Ising layer is

```text
g_{2xL}(Q) = so_m(Q) (+) so_m(Q)     when L is even,
             sp_m(Q) (+) sp_m(Q)     when L is odd,

dim_Q g_{2xL} = 2^(n-1) (2^(n-1) - (-1)^L).
```

**[THEOREM]** In particular, the all-`L` saturation step left open in
`proofs/alternation_law.md` is now closed: the reachable Pauli labels equal the
entire alternation-law target set for every `L >= 2`.  This is an all-size
proof, not an extrapolation from `L <= 7`.

## 1. Exact production system and target

**[LEMMA]** Let `V_L = {(r,c): r in {0,1}, 1 <= c <= L}` and let `E_L`
contain the `L` rung edges and the `2(L-1)` horizontal edges.  Put
`J = sum_{v in V_L} e_v`, and let `d` be either checkerboard color class.
Every edge `uv` satisfies `(e_u+e_v).d = 1` over `F_2`.

**[LEMMA]** Use the rational Pauli monomial
`M(a,b) = X^a Z^b`, in a fixed site order, for `(a|b) in F_2^(2n)`.  Define

```text
<(a|b),(c|e)> = a.e + b.c  in F_2,
(a|b) + (c|e) = (a+c | b+e).
```

**[LEMMA]** The exact generator set and exact production rule are

```text
G_L = {(e_v|0): v in V_L} union {(0|e_u+e_v): uv in E_L},

(P)  p,q already reached and <p,q>=1  ==>  p+q is reached.
```

**[LEMMA]** Rule `(P)` is sound over `Q`.  Indeed,
`M(a,b)M(c,e)=(-1)^(b.c)M(a+c,b+e)`, so an odd symplectic pairing gives
`[M(a,b),M(c,e)] = +/- 2 M(a+c,b+e)`.  Thus every produced label lies in
`g_L`.  Conversely, every nonzero nested commutator of generators is a
nonzero rational multiple of one label produced by `(P)`, and the Pauli
monomials are linearly independent.  Therefore the span of the reached-label
set is exactly `g_L` and its cardinality is `dim_Q g_L`.

**[LEMMA]** The target predicate is

```text
S_L = {(a|b): J.b=0,
       [b != J and (J+b).a = 1 + b.d]
       or [b = J and L is odd]}.
```

**[LEMMA]** The `b=J` clause has no restriction on `a`.

**[LEMMA]** Every generator lies in `S_L`, and generator-monotone production
cannot leave `S_L`.  For `X_i`, a nonzero bracket requires `b_i=1`, while the
coefficient `(J+b)_i` is then zero, so toggling `a_i` preserves the predicate.
For a bond `h=e_u+e_v`, a nonzero bracket requires `h.a=1`; when `b+h != J`,

```text
(J+b+h).a = (J+b).a + h.a
            = 1+b.d+1 = 1+(b+h).d,
```

because `h.d=1`.  If `b+h=J`, the source law forces
`h.a = 1+(J+h).d = L (mod 2)`: the bracket vanishes for even `L` and lands in
the allowed full `J` fiber for odd `L`.  Starting from the `J` fiber when `L`
is odd, a nonzero bond bracket has `h.a=1`, exactly the ordinary-fiber law for
`J+h`.  Hence the already-proved containment `g_L subseteq span_Q(S_L)` is
recovered directly.

**[LEMMA]** Consequently, the arbitrary reached-pair rule `(P)` used below
also cannot escape `S_L`: its output is a single Pauli monomial in `g_L`, while
`g_L subseteq span_Q(S_L)` and Pauli monomials are linearly independent.

## 2. The bounded production schemas

**[LEMMA] (edge fiber)** For every edge `uv`, all four labels

```text
(alpha e_u + beta e_v | e_u+e_v),  alpha,beta in F_2,
```

are reached.  Start from `(0|e_u+e_v)` and apply `X_u` and `X_v`; each endpoint
application has odd pairing.

**[LEMMA] (honest telescope)** If `u-v-w` is a length-two path, then all four
labels

```text
(alpha e_u + e_v + beta e_w | e_u+e_w),  alpha,beta in F_2,
```

are reached.  They are the products of the two edge-fiber labels
`(alpha e_u|e_u+e_v)` and `(e_v+beta e_w|e_v+e_w)`, whose pairing is one.
The interior bit is forced to one.  Thus this lemma deliberately supplies
only the four endpoint choices; it does not repeat the wave-17 false claim of
arbitrary path-interior freedom.

**[LEMMA] (second independent local move)** Let `p-r-s-q-p` be one ladder
square.  Every pure-`X` three-corner label is reached:

```text
(p+q+r|0), (q+r+s|0), (p+r+s|0), (p+q+s|0).
```

For the diagonal `{p,s}`, the two length-two routes give
`(e_r|p+s)` and `(e_p+e_q|p+s)`; their pairing is one and their product is
`(p+q+r|0)`.  Replacing the second parent by `(e_q+e_s|p+s)` gives
`(q+r+s|0)`.  The two routes across diagonal `{r,q}` give the other two
triples in the same way.  This bounded square move is the missing second move:
it creates directions not available from a single forced-interior telescope.

## 3. No-rank-loss transport of a complete pair fiber

**[LEMMA] (transport theorem)** For a two-site support `z`, write

```text
F_L(z) = {a: (a|z) in S_L},
D_L(z) = {delta: (J+z).delta=0}.
```

For every pair `z`, `F_L(z)` is an affine coset with tangent `D_L(z)` of
dimension `n-1`.  Suppose `z'={p,r}` and `z={p,s}`, where `r-s` is an actual
ladder bond and `p,r,s` are distinct.  If the whole source fiber `F_L(z')` is
reached, then the whole target fiber `F_L(z)` is reached.

**[LEMMA] (proof of transport theorem)** Put `h=e_r+e_s`.  Bracketing the
source with bond `(0|h)` transfers exactly the affine half satisfying
`a.h=1`.  This half is nonempty because `e_r in D_L(z')` and `e_r.h=1`.  Its
difference space is

```text
E = D_L(z') intersect ker(h.),       dim E = n-2.
```

Moreover, `J+z = (J+z')+h`, so `E subset D_L(z)`.  The newly occupied endpoint
satisfies `e_s in D_L(z)` but `e_s.h=1`, hence `e_s notin E`.  Applying the
single-site generator `X_s` to every transferred target label is legal because
`s in z`, and it supplies the missing direction.  The reached target set has
difference space

```text
E (+) span{e_s} = D_L(z)
```

and therefore has all `2^(n-1)` target points.  No interior-degree assumption
was used: at the boundary columns `c=1` and `c=L`, the argument still uses only
the actual moved bond `r-s` and the always-present local generator `X_s`.

**[LEMMA] (two-token connectivity)** If a finite graph `G` is connected and
has at least three vertices, then the graph whose vertices are two-element
subsets of `V(G)` and whose edges slide one element across an edge of `G` into
an unoccupied vertex is connected.

**[LEMMA] (proof of two-token connectivity)** It is enough to use a spanning
tree `T`.  Induct on `|V(T)|`.  For a leaf `v` with neighbor `u`, configurations
not containing `v` are connected inside the two-token graph of `T-v`.  A
configuration `{v,w}` with `w != u` slides `v` to `u`.  For `{v,u}`, when the
tree has at least three vertices, slide `u` to another neighbor `t` and then
slide `v` to `u`; the result avoids `v`.  The two-vertex base has one
configuration.  Thus every configuration connects to the inductive subgraph,
so the two-token graph is connected for every `L`, not merely for tested
widths.

## 4. Uniform induction for all weight-two fibers

**[THEOREM] (inductive invariant)** For every `L >= 2`, the invariant `I(L)`
holds:

```text
I(L): for every two-element z subset V_L, every label (a|z) in S_L is reached.
```

**[LEMMA] (base `L=2`)** On the single square `p-r-s-q-p`, take the top edge
`z={p,r}`.  Its edge fiber gives the two independent endpoint directions
`e_p,e_r`.  The square triple `e_p+e_q+e_s` intersects `z` once, so it legally
toggles the edge fiber and supplies a third independent direction.  The target
tangent has dimension three; hence this fiber is complete.  The transport
theorem and two-token connectivity complete all six pair fibers, proving
`I(2)`.

**[LEMMA] (step `I(L) => I(L+1)`)** Embed the width-`L` ladder as the first
`L` columns of width `L+1`.  Let `t,p` be the last old top horizontal edge,
let `q` be the bottom site under `p`, and let `r,s` be the new top and bottom
sites.  By `I(L)`, the old fiber at `z={t,p}` is complete and remains reached
after appending two zero coordinates.  The new square `p-r-s-q-p` supplies

```text
tau_1=e_p+e_q+e_r,       tau_2=e_p+e_q+e_s.
```

Each intersects `z` only at `p`, so both are legal toggles.  Their projections
to the two new coordinates are `(1,0)` and `(0,1)`, so they are independent
modulo the embedded old tangent.  The reached affine set therefore has size

```text
2^(2L-1) * 4 = 2^(2(L+1)-1),
```

which is the full new fiber size.  The transport theorem and two-token
connectivity then fill every pair fiber at width `L+1`.  This proves `I(L+1)`
with no circular use of weight-zero or higher-weight fibers.

**[THEOREM]** The base and step prove `I(L)` for every `L >= 2`.

## 5. From pair fibers to every target fiber

**[LEMMA] (all odd pure-X strings)** Every `(a|0)` with `J.a=1` is reached.
Singletons are generators.  Otherwise, because `n` is even, the odd support of
`a` is a nonempty proper subset of the connected ladder, so a boundary edge
`z` satisfies `z.a=1`.  In the complete pair fiber choose any `x`; since
`(J+z).a=1+1=0`, both `(x|z)` and `(x+a|z)` lie in that fiber.  Their pairing is
`z.a=1`, and their product is `(a|0)`.

**[LEMMA] (one seed for every ordinary even support)** Let `b` be nonzero,
even, and different from `J`.  Partition its support into disjoint pairs
`z_1,...,z_w`.  Pair fibers are complete.  Suppose a reached partial product
has support `C=z_1+...+z_{j-1}`.  On `F_L(z_j)`, the required pairing depends
affinely on `C.x`.  This functional is nonconstant on `D_L(z_j)` unless
`C=J+z_j`; that exceptional equality would force the final support to be `J`
and is excluded.  Choose the required value and multiply.  Induction on `j`
produces one reached seed in `F_L(b)`.

**[LEMMA] (ordinary fiber completion)** For nonzero even `b != J`, let

```text
Q_b = {q: J.q=1 and b.q=1}.
```

Every `(q|0)` for `q in Q_b` is reached and legally toggles a label with support
`b`.  The set `Q_b` is an affine coset of
`ker(J.) intersect ker(b.)`, of dimension `n-2`; its linear span has dimension
`n-1` and equals `D_L(b)=ker((J+b).)`.  Therefore the seed from the preceding
lemma expands to the whole fiber.  The case `b=0` is exactly the odd pure-`X`
fiber already proved.

**[LEMMA] (the exceptional `J` fiber)** Let `h` be any ladder bond and put
`b'=J+h`.  For odd `L`, the ordinary-fiber law forces

```text
h.a = 1+b'.d = 1+J.d+h.d = L = 1 (mod 2).
```

Thus any reached seed at `b'` brackets with `(0|h)` to produce a `J`-support
label.  The singleton `X_i` generators then toggle all `n` directions, so the
entire `J` fiber is reached.  For even `L` the same forced value is zero, so the
last bracket vanishes and `J` remains excluded, exactly as the target law
requires.

**[THEOREM] (saturation)** Every label in `S_L` is reached for every `L >= 2`.
Together with the forward containment, the reached set equals `S_L`.

## 6. Dimension, form, type, and structural explanation

**[THEOREM] (dimension)** There are `2^(n-1)-1` ordinary even supports
`b != J`, each with `2^(n-1)` labels.  For odd `L`, the exceptional `J` fiber
adds `2^n` labels.  Hence

```text
|S_L| = 2^(n-1)(2^(n-1)-1)  for even L,
        2^(n-1)(2^(n-1)+1)  for odd L.
```

Pauli independence gives the stated dimension of `g_L`.

**[THEOREM] (form and type)** The all-`L` invariance identity already proved in
`proofs/alternation_law.md` places each product-`X` half in the isometry algebra
of the nondegenerate signed-checkerboard matching form `Omega_L`.  It is
symmetric for even `L` and alternating for odd `L`.  The saturation dimension
is exactly twice `dim so_m = m(m-1)/2` or twice `dim sp_m = m(m+1)/2`, so the
containments are equalities.  On each half the invariant form is unique up to
scale.  In the symmetric case, checkerboard translation pairs the class basis
without fixed points, and every signed `2 x 2` block is hyperbolic; hence the
orthogonal form is split over `Q`.

**[LEMMA] (independent Clifford cross-check, restricted to `2 x L`)** The
Hamiltonian-path/Jordan-Wigner construction in this paragraph was supplied by
the wave-18 `AllLQuadraticNoGo` sibling front; its later general-graph
extension belongs to the separate `CliffordGradeClass` front and is not claimed
here.  Order the sites around the perimeter Hamiltonian path and use the
rational Majoranas

```text
Gamma_(2j-1) = X_1...X_(j-1) Z_j,
Gamma_(2j)   = X_1...X_(j-1) (X_j Z_j).
```

The fields and path bonds are consecutive bilinears, so their brackets produce
all grade-two monomials.  The remaining rung at column `c` is one monomial of
grade `4L-4c+2`; together the seeded grades are exactly
`2,6,10,...,4L-2`.

**[LEMMA] (Clifford orbit coverage)** Bracketing a grade-`k` basis monomial
with a grade-two monomial that has exactly one index in its support replaces
that index by the other index.  The one-index-swap graph on `k`-subsets is
connected, so one rung spans every monomial of its grade, including the middle
grade `k=2L` when `L` is odd.  Thus no middle-grade half-orbit assumption is
used by the primary induction; even this independent route covers the middle
directly over `Q`.  Equivalently, after base change to an algebraic closure, a
middle basis monomial has nonzero projections to both Hodge-star summands; the
orbit span is `Q`-rational, so equality after base change descends to `Q`.

**[LEMMA] (Clifford dimension identity)** Grades congruent to `2 mod 4` are
closed under commutator: two even monomials anticommute only when their support
intersection has odd size, and the output grade remains `2 mod 4`.  The
root-of-unity filter gives

```text
sum_(k = 2 mod 4) C(4L,k)
 = (2^(4L) - (1+i)^(4L) - (1-i)^(4L))/4
 = 2^(4L-2) - (-1)^L 2^(2L-1),
```

which is exactly the saturation dimension.  Reversal acts on grade `k` by
`(-1)^(k(k-1)/2)`, so this is its `-1` eigenspace.  On the half-spin module for
`so(2n)`, the invariant form is symmetric for `n=0 mod 4` and alternating for
`n=2 mod 4`; since `n=2L`, this structurally forces the even-`L` orthogonal /
odd-`L` symplectic alternation.  This Clifford route is an independent
cross-check, not a premise of Sections 2--5.

## 7. Exact finite verification and resource record

**[COMPUTATION]** The producer retained exact Pauli labels in every pair fiber,
checked every ordinary even-support row without sampling, and performed the
Clifford integer census at `L=2,...,7`.  It additionally enumerated the full
generator-monotone closure and checked target-law membership for every reached
label at `L=2,...,6`.  All timings below use `time.process_time()`; full-closure
RSS is measured in an isolated child process.

| `L` | exact dimension | full closure enumeration and membership | exact construction/rank check | measured CPU / peak RSS |
|---:|---:|---|---|---|
| 2 | 56 | **[COMPUTATION]** 56 reached, 0 violations | **[COMPUTATION]** all 7 ordinary rows; 6 pair fibers rank 3 | closure 0.000090 s / 28,770,304 B; construction 0.000466 s / 28,377,088 B cumulative |
| 3 | 1,056 | **[COMPUTATION]** 1,056 reached, 0 violations | **[COMPUTATION]** all 31 ordinary rows; 15 pair fibers rank 5 | closure 0.001822 s / 28,393,472 B; construction 0.001768 s / 28,377,088 B cumulative |
| 4 | 16,256 | **[COMPUTATION]** 16,256 reached, 0 violations | **[COMPUTATION]** all 127 ordinary rows; 28 pair fibers rank 7 | closure 0.036629 s / 29,720,576 B; construction 0.005739 s / 28,573,696 B cumulative |
| 5 | 262,656 | **[COMPUTATION]** 262,656 reached, 0 violations | **[COMPUTATION]** all 511 ordinary rows; 45 pair fibers rank 9 | closure 0.938014 s / 52,494,336 B; construction 0.047740 s / 29,704,192 B cumulative |
| 6 | 4,192,256 | **[COMPUTATION]** 4,192,256 reached, 0 violations; new full-closure confirmation | **[COMPUTATION]** all 2,047 ordinary rows; 66 pair fibers rank 11 | closure 19.089286 s / 436,371,456 B; construction 0.302009 s / 35,880,960 B cumulative |
| 7 | 67,117,056 | **[COMPUTATION]** **not a full closure enumeration** | **[COMPUTATION]** all 8,191 ordinary rows; 91 pair fibers rank 13; independent Clifford sum | construction 2.288197 s / 67,862,528 B cumulative |

**[COMPUTATION]** The exact producer artifact is
`results/algebra_growth/alll_saturation.json`.  Its `data.checks` list contains
only passed checks and records the exact environment, input hashes, per-point
evidence strength, process time, and RSS.  The standalone verifier independently
reconstructs every decisive quantity and repeats the full closure through
`L=6` without importing producer code or conclusions.

**[COMPUTATION]** Reproduction commands from the repository root are:

```bash
.venv/bin/python experiments/e157_saturation_local_move.py
.venv/bin/python experiments/e158_pair_fiber_bfs.py
.venv/bin/python experiments/e159_alll_saturation.py
.venv/bin/python tests/test_alll_saturation.py
```

## 8. Status by range and honest scope

**[THEOREM]** The status ledger after this proof is:

| range | saturation | equality/dimension | type | invariant form | finite evidence in this artifact |
|---|---|---|---|---|---|
| `L=2,...,5` | **[THEOREM]** | **[THEOREM]** | **[THEOREM]** even orthogonal / odd symplectic | **[THEOREM]** explicit, nondegenerate, unique per half | **[COMPUTATION]** full closure plus construction |
| `L=6` | **[THEOREM]** | **[THEOREM]** | **[THEOREM]** orthogonal | **[THEOREM]** symmetric split matching, unique per half | **[COMPUTATION]** new 4,192,256-string full closure plus construction |
| `L=7` | **[THEOREM]** | **[THEOREM]** | **[THEOREM]** symplectic | **[THEOREM]** alternating matching, unique per half | **[COMPUTATION]** exhaustive construction/ranks and Clifford census; no full closure enumeration |
| every `L>=8` | **[THEOREM]** | **[THEOREM]** | **[THEOREM]** parity-alternating classical type | **[THEOREM]** explicit matching, unique per half | **[UNRESOLVED]** no finite computation claimed here |

**[THEOREM] (Gaussian consequence)** Every injective,
Lie-bracket-preserving realization of these local terms by quadratic operators
on `M` fermionic modes has image dimension at most
`dim so(2M)=M(2M-1)`.  Therefore

```text
M >= ceil((1 + sqrt(1+8 dim g_L))/4),
```

which grows exponentially in `L`.  On the physical `2^n`-dimensional spin
space one has `M=n=2L`; for `n>=4`, `2^(n-1)>=2n`, so
`dim g_L >= 2^(n-1)(2^(n-1)-1) > n(2n-1)`.  Hence no change of basis on the
physical space makes all open-`2 x L` local terms Majorana-quadratic, for any
`L>=2`.

**[UNRESOLVED]** The preceding mode bound does not cover a non-invariant code
compression `U^dagger q U`, because compression need not preserve commutators.

**[UNRESOLVED]** This theorem is restricted to open `2 x L` ladders and
characteristic zero (in particular `Q`).  It does not classify periodic
ladders, higher-row layers, or arbitrary graphs; it does not solve the 3D
Ising model or determine its free energy.  The general Hamiltonian-path graph
classification is deliberately outside this file's ownership and scope.
