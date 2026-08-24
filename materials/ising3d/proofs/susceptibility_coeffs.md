# Exact susceptibility coefficients: concatenation inequalities, the smallest counterexample, and the endpoint-route no-go

Artifacts: `experiments/e218_susceptibility_coeffs.py`,
`experiments/e219_coeff_inequality.py`, `experiments/e220_coeff_bound.py`,
`results/bounds/susceptibility_coeffs.json`,
`tests/test_susceptibility_coeffs.py`.

## Scope

[COMPUTATION]/[THEOREM] This note derives the simple-cubic zero-field
high-temperature susceptibility coefficients exactly from finite graphs,
settles the direction question for coefficient-level concatenation
inequalities, and audits every Fekete / radius-of-convergence route to a
critical-coupling endpoint.  It does **not** solve the three-dimensional
Ising model, it proves **no** new endpoint, and it promotes **no** finite
computation to an all-size statement.  The two headline theorems are a
smallest exact counterexample (Theorem 1) and a finite-verification no-go
for the proposed endpoint route (Theorem 4).  The benchmark value
`0.221654626(5)` is used nowhere.

## 1. Conventions and the exact objects

For a finite simple graph `G = (V, E)` with uniform coupling, `v = tanh K`,
free boundary conditions and zero field, the high-temperature representation
used throughout the repository (`proofs/kc_bounds.md`, eq. (A.2)) is

    <s_a s_b>_G = N_G(a,b) / Z_G,
    N_G(a,b) = sum_{A subset E, dA = {a,b}} v^|A|,
    Z_G      = sum_{A subset E, dA = empty} v^|A|,

where `dA` is the set of odd-degree vertices of `A`.  For a box
`Lambda subset Z^3` write

    chi_Lambda(v) = sum_{x in Lambda} <s_0 s_x>_Lambda .

Both `N` and `Z` are polynomials with nonnegative integer coefficients and
`Z(0) = 1`, so `<s_0 s_x>_Lambda` and `chi_Lambda` are well-defined formal
power series with integer coefficients.

**Definition (susceptibility coefficients).**  `a_n` is the common value of
`[v^n] chi_Lambda` over all boxes `Lambda` containing the ball
`B(0, 2n) = {x : ||x||_1 <= 2n}` centred at the source.  Lemma 3 proves the
common value exists.  This is the standard formal HT susceptibility series;
its prefix `1, 6, 30, 150, 726, ...` coincides with the published series
(OEIS A002913, [EXTERNAL], comparison only).

## 2. The cluster identity

**[LEMMA 1] (decomposition).**  On every finite graph, for distinct `a, b`:

    <s_a s_b>_G = sum_C v^|C| * Z_G^{V \ V(C)} / Z_G ,

where `C` runs over connected edge sets with `dC = {a, b}` exactly, `V(C)`
is the vertex set of `C`, and `Z_G^{W}` counts fully even edge subsets of
the subgraph induced on `W`.

*Proof.*  Let `A` satisfy `dA = {a,b}`.  The connected components of `A`
partition its edges and are pairwise vertex-disjoint.  All edges of `A`
incident to `a` lie in one component (they meet at `a`), so `a`'s full
degree is realized in a single component; a component has an even number of
odd vertices, so the component containing `a` has odd set exactly `{a,b}`
and every other component is even.  Write `A = C ⊔ B` with `C` that
component; then `B` is an arbitrary fully even subset of the graph induced
on `V \ V(C)` (its components can share no vertex with `C`), and
`|A| = |C| + |B|`.  The map `A -> (C, B)` is a bijection onto pairs
(connected two-odd `C` with `dC = {a,b}`) x (even `B` on `V \ V(C)`).
Summing `v^|A|` gives `N_G(a,b) = sum_C v^|C| Z_G^{V \ V(C)}`; divide by
`Z_G`.  QED.

**[LEMMA 2] (dressing recursion).**  For vertex sets `U subset X` put
`g_X(U) = Z^{X\U} / Z^{X}`.  Then

    1 / g_X(U) = 1 + sum_{B1} v^|B1| * g_{X\U}(V(B1) \ U),

where `B1` runs over nonempty even edge subsets of `E(X)` all of whose
components have a vertex in `U`.

*Proof.*  Split any even `B subset E(X)` as `B = B1 ⊔ B2`, `B1` = union of
components meeting `U`.  Then `B2` is an arbitrary even subset of the graph
on `X \ U \ V(B1)`: its components avoid `U` by definition of the split and
avoid `V(B1)` by disjointness of components.  The split is bijective, so
`Z^X = sum_{B1} v^|B1| Z^{(X\U) \ (V(B1)\U)}`.  Divide by `Z^{X\U}`.  QED.

Every even connected subgraph of a bipartite graph has even size (Veblen
decomposition into edge-disjoint cycles; all cycles of `Z^3` are even), so
each `B1` has size in `{4, 6, 8, ...}` and the recursion terminates modulo
`v^{K+1}` after at most `K/4` levels.  Series inversion of `1 + q` with
`q(0) = 0` keeps integer coefficients, so all `g` coefficients are integers.

**[LEMMA 3] (locality and stabilization).**  Modulo `v^{n+1}`:
(i) `g_X(U)` depends only on the isomorphism type of the ball of graph
radius `n` around `U` inside `X`; (ii) `[v^n] <s_0 s_x>_Lambda` is the same
for every box containing `B(0, 2n)`; (iii) `[v^n] chi_Lambda` is the same
for every box containing `B(0, 2n)`, which defines `a_n`, and

    a_n = sum_{[C]} 2 * [v^{n-|C|}] g(V(C)),

where `[C]` runs over translation classes of connected two-odd clusters
with `|C| <= n` and `g` is the dressing on `Z^3`.

*Proof.*  (i) Fully expanding Lemma 2 to depth `<= n/4` writes
`g_X(U) mod v^{n+1}` as a finite integer combination of products over chains
`(B_1, ..., B_k)` in which `B_1` meets `U`, each `B_{j+1}` meets
`V(B_j) \ (previous removals)`, and `sum |B_j| <= n`.  Every edge of every
chain therefore lies within graph distance `n` of `U`, and each summation
range is determined by that ball.  (ii) By Lemma 1 with `a = 0`, `b = x`,
the order-`n` coefficient involves clusters `C` with `|C| <= n` (hence
inside `B(0, n)`, since `C` is connected and contains `0`) dressed by
`g(V(C))` truncated at `n - |C|`; by (i) these depend only on the ball of
radius `n` around `V(C)`, which lies in `B(0, 2n)`.  (iii) Only `x` with
`||x||_1 <= n` contribute at order `n` (a contributing `C` joins `0` to `x`
with at most `n` edges).  Summing (ii) over `x`: the pairs `(x, C)` with
`dC = {0, x}` are exactly the clusters with `0` in their odd pair, and each
translation class contains precisely two placements of its odd pair on the
origin.  QED.

Lemma 3 also yields the crude but fully explicit growth envelope used by
Theorem 4:

**[LEMMA 6] (envelope).**  Without assuming any coefficient sign,

    |a_n| <= (1000 (n+1)^3)^{n+2}.

*Proof.*  Take `Lambda ⊇ B(0,2n)` inside the box of side `4n+1`, with
`E <= 3 (4n+1)^3` edges.  `|[v^k] N| <= C(E,k) <= E^k` and, writing
`Z = 1 + p` with `|[v^t] p| <= E^t`, the composition expansion of
`1/(1+p)` gives `|[v^j] 1/Z| <= sum_{k>=1} C(j-1, k-1) E^j = 2^{j-1} E^j`.
Hence `|[v^n] N/Z| <= n 2^n E^n`, and summing over at most `(4n+1)^3`
targets, `|a_n| <= (4n+1)^3 n 2^n E^n <= (125(n+1)^3)(n+1) (750 (n+1)^3)^n
<= (1000 (n+1)^3)^{n+2}`.  QED.

## 3. The exact coefficients

[COMPUTATION] `experiments/e218_susceptibility_coeffs.py` enumerates every
translation class of connected clusters with at most `N = 10` edges
(anchored at the lexicographically least edge; each class generated exactly
once), classifies odd sets, and dresses every two-odd cluster through
Lemma 2.  The derivation uses **no published input**.  Result:

    n   : 0  1  2   3    4    5     6      7      8       9        10
    a_n : 1  6  30  150  726  3510  16710  79494  375174  1769686  8306862

with census `two_odd = [., 3, 15, 75, 363, 1815, 8775, 43053, 207639,
1010405, 4863651]`.  Controls, all recorded as computed checks in the
artifact:

* the complete brute-force parity enumeration on the `2x2x2` box equals the
  cluster identity for all seven targets through `v^8`, and on `3x2x2` for
  two targets (this validates Lemmas 1-2 and the code on graphs where
  complete enumeration is possible);
* hand-countable census anchors (3 edge classes, `15 = C(6,2)` wedges,
  `20 = C(6,3)` stars, 3 plaquettes, 20 four-cycles meeting an edge);
* the `N = 8` census pass is a prefix of the final census;
* an independent SAW enumeration gives `a_n = c_n` for `n <= 4` and
  `a_n < c_n` for `5 <= n <= 10` with exact gaps
  `24, 216, 1896, 12792, 84200, 503016`;
* the order-5 identity `a_5 = 2*1815 - 2*3*20 = 3510` and the order-6
  identity `a_6 = 2*8775 - 2*15*28 = 16710` (every wedge is dressed by
  exactly 28 four-cycles);
* [EXTERNAL cross-check] the derived prefix equals OEIS A002913 (Butera--
  Comi PRB 62 (2000) 14837 through `v^23`; Campostrini--Pelissetto--Rossi--
  Vicari PRE 65 (2002) 066127 through `v^25`; Fujiwara--Arisue workshop
  slides for `n >= 26` -- the last is a weaker provenance grade and enters
  no theorem-grade claim).

Telemetry: 2808 CPU s, peak RSS 289 MiB, single process.

## 4. Theorem 1: the smallest exact counterexample

**[THEOREM 1] (supermultiplicativity is false).**  The susceptibility
coefficients are **not** supermultiplicative: the smallest possible pair is
already a counterexample,

    a_2 = 30 < 36 = a_1 * a_1 ,

and both sides are theorem-grade finite-graph counts: order `<= 2`
coefficients receive no dressing (the smallest even cluster has four
edges), so `a_1 = 2 * #{edge classes} = 2 * 3` and
`a_2 = 2 * #{two-edge path classes} = 2 * C(6,2) = 30`.

Moreover, every anchored variant fails at its own smallest self-pair:
`a_{2k} < a_k^2` for every `k <= 5` (theorem-grade) and every `k <= 16`
([EXTERNAL]-grade), and no pair `(m, n)` with `m + n <= 32` satisfies
`a_{m+n} >= a_m a_n`.

*Consequence.*  A Fekete argument in the sup direction (supermultiplicative
`a_n` would give `lim a_n^{1/n} = sup_n a_n^{1/n}`, hence per-`n` **upper**
endpoints `K_c <= atanh(a_n^{-1/n})`) is unavailable for this sequence: its
hypothesis fails at `(1,1)`.  This is the exact reason the proposed
upper-endpoint route through susceptibility coefficients dies.  [THEOREM]

## 5. Theorem 2: the direction audit

**[THEOREM 2] (what each inequality direction would yield).**  Let
`t_n = atanh(a_n^{-1/n})`.

| statement (all `m, n >= 1`) | Fekete consequence | endpoint direction | status |
|---|---|---|---|
| `a_{m+n} <= a_m a_n` together with `a_n > 0` | `lim a_n^{1/n} = inf_n a_n^{1/n} <= a_n^{1/n}` | lower: `K_c >= t_n` for every `n`, **conditional also on P1'** (Theorem 5) | both all-`n` positivity (G0) and the inequality (G1) are [CONJECTURE]: true throughout the checked range, unprovable from finite data (Theorem 4) |
| `a_{m+n} >= a_m a_n` | `lim = sup` | upper: `K_c <= t_n` | **[THEOREM] false at (1,1)** |
| `a_n^2 >= a_{n-1} a_{n+1}` (log-concavity) | implies the sub row (Lemma 5) | lower | [CONJECTURE]: non-strict throughout the checked range, with equality at `n=2` |
| `a_{n+1} >= rho a_n` (ratio floor) | `a_n >= c rho^n`, radius `<= 1/rho` | upper: `K_c <= atanh(1/rho)`, conditional on complex-disk analyticity of `chi` up to the radius (open) | [UNRESOLVED]: no finite prefix certifies any floor -- the observed ratios are nonincreasing, with the unique checked plateau `a_2/a_1=a_3/a_2=5` |
| `a_{n+1} <= 6 a_n` (ratio ceiling) | radius `>= 1/6` | lower: `K_c >= atanh(1/6)` | holds at every checkable index, but **valueless**: `atanh(1/6) <= 0.16823611831060646526 < 0.2122159753...` |

**[LEMMA 5].**  If `a_0 = 1`, `a_n > 0`, and `a_n^2 >= a_{n-1} a_{n+1}` for
all `n >= 1`, then `a_{m+n} <= a_m a_n` for all `m, n >= 0`.

*Proof.*  Log-concavity says `r_k = a_{k+1}/a_k` is nonincreasing.  Then
`a_{m+n}/a_n = prod_{k=n}^{n+m-1} r_k <= prod_{k=0}^{m-1} r_k = a_m / a_0`.
QED.  (So proving log-concavity for all `n` would suffice for the sub row;
it is exactly as open.)

[COMPUTATION] All finite checks: submultiplicativity holds at every pair
with `m+n <= 10` (theorem-grade) and `m+n <= 32` ([EXTERNAL]); log-concavity
holds non-strictly on both ranges, with equality at `n=2`, and the coefficient
ratios are otherwise strictly decreasing from `6 = a_1/a_0`.

## 6. Theorem 3: the natural bridge premises, and their exact failures

Even granting the sub row's hypothesis, deducing `K_c >= t_n` needs a bridge
from coefficient growth to the model.  The naive bridge premises are
**false**, by exact finite counterexamples:

**[COMPUTATION] (negative box coefficients).**  On the `2x2x1` box with a
corner source,

    chi_box = 1 + 2v + 2v^2 + 2v^3 + 0v^4 - 2v^5 - 2v^6 - 2v^7 + 0v^8 + ...

(hand check: `chi = 1 + (2v + 2v^2 + 2v^3)/(1 + v^4)`).  So
`[v^5] chi_box = -2 < 0`: box susceptibility coefficients are not
nonnegative, and no monotone-convergence bridge exists.

**[COMPUTATION] (coefficient-wise volume monotonicity fails).**
`[v^5] <s_0 s_{e_1}>` drops from `0` on `2x1x1` to `-1` on `2x2x1` although
the first box is a subgraph of the second containing both sites: the
coefficient-wise analogue of the GKS volume monotonicity is false.

**[CONJECTURE G0] (all-order positivity).**  `a_n > 0` for every `n`.
Every derived/published coefficient through `n = 32` is positive; no
all-order proof is supplied here.  (P1' below would imply nonnegativity, but
strict positivity is stated separately so the ordinary logarithmic form of
Fekete is never applied silently.)

**[CONJECTURE P1'] (absolute domination).**  `|[v^n] chi_Lambda| <= a_n`
for every finite box, every source, and every `n`.  Finite evidence: true
on every tested box (`2x1x1` through `3x3x3`, corner and centre sources,
`n <= 8`), with equality attained.  This is the weakest natural premise
that still bridges (Theorem 5); it is open.

**[LEMMA 4] (conditional bridge).**  Assume P1'.  If `sum_n a_n v^n`
converges at some `v` with `0 <= v < 1`, then
`sup_Lambda chi_Lambda(v) <= sum_n a_n v^n < infinity`, and consequently
there is no spontaneous magnetization at `v`; under the repository's
critical-coupling conventions (`proofs/kc_bounds.md`, section A.1) every
such `v` satisfies `v <= v_c`, i.e. `K_c >= atanh(v)` for every
`v < radius`.

*Proof.*  `chi_Lambda` is a rational function of `v`; P1' bounds its Taylor
coefficients at `0` by `a_n`, so its Taylor radius is at least the radius
`rho` of `sum a_n v^n`; a rational function's Taylor radius is the distance
to its nearest pole, so `chi_Lambda` is analytic on `|v| < rho` and equals
its Taylor sum there; term-by-term `chi_Lambda(v) <= sum |[v^n]chi_Lambda|
v^n <= sum a_n v^n` for `0 <= v < rho`.  For the magnetization: by GKS the
free-boundary `<s_0 s_x>_Lambda` increases in `Lambda`, so partial sums
over any fixed finite set of targets are bounded by
`sup_Lambda chi_Lambda(v)`; hence `chi(v) = sum_x <s_0 s_x> < infinity` and
the shell sums `S_r = sum_{x in dB_r} <s_0 s_x>` satisfy
`liminf_r S_r = 0`.  Lieb's inequality ([EXTERNAL THEOREM], E. H. Lieb,
Commun. Math. Phys. 77 (1980) 127, doi:10.1007/BF01982712, already in this
repository's dependency ledger; the separator form with free boundary
conditions on the inner box) applied to the ghost-spin representation of
plus boundary conditions gives, for every outer box and every `r`,
`m_+ = <s_0 s_ghost> <= sum_{x in dB_r} <s_0 s_x>_{B_r, free} * 1 <= S_r`,
using `<s_0 s_x>_{B_r,free} <= <s_0 s_x>` (GKS volume monotonicity, at
positive `v`, function level).  Letting `r` run along the vanishing
subsequence, `m_+ = 0` uniformly in the outer box.  QED.

**[THEOREM 5] (the conditional route, fully assembled).**  IF
`a_n > 0` for all `n` (G0), `a_{m+n} <= a_m a_n` for all `m, n >= 1`
(G1), AND P1' holds, THEN for every exactly known `a_n`,

    K_c >= atanh(a_n^{-1/n}) .

*Proof.*  G0 and G1 allow Fekete's lemma for the subadditive sequence
`log a_n`, hence
`lim a_m^{1/m} = inf_m a_m^{1/m} <= a_n^{1/n}`.  By
Cauchy--Hadamard the series radius is therefore at least
`a_n^{-1/n}`, and Lemma 4 applies to every `v < a_n^{-1/n}`.  QED.

[COMPUTATION] Certified 40-decimal directed enclosures (exact integer
`n`-th-root brackets with denominator `10^40`, then rational `atanh` Taylor
enclosures with exact tail bounds; no floating point):

    t_10 in [0.2061328462722200968354785911248920631301, ...303]   (theorem-grade a_10)
    t_25 in [0.2133819845800037925035678837191924967804, ...806]   ([EXTERNAL] peer-reviewed grade)
    t_32 in [0.2147544289685755010198749449249222101981, ...983]   ([EXTERNAL] slide grade)

Against the incumbent certified lower endpoint
`0.2122159753270231627267517174278577806020` (wave 14, `proofs/saw_union4.md`):
the conditional route would first beat the incumbent at order `n = 21`, and
`t_25`, `t_32` would improve it by `1.2e-3` resp. `2.5e-3` -- **if** G0,
G1 and P1' were theorems and the external integers were audited at their
primary tables.  None of that is available, so **no bound is claimed**
(`data.directed.claimed_new_bound = null`).  Within the theorem-grade range
the route is worthless even conditionally: `t_10 < incumbent`.

## 7. Theorem 4: finite verification cannot carry the route

**[THEOREM 4] (no-go for finite evidence).**  Let `F` be the full verified
constraint set: the exact values `a_0..a_10` (theorem-grade), the
[EXTERNAL] values `a_11..a_32`, every inequality among them checked in this
front (all pairs `m+n <= 32`), and the provable envelope of Lemma 6.  Then
`F` entails no upper bound on `limsup b_n^{1/n}` below `6` and no lower
bound above `inf_{n<=32} a_n^{1/n}`; in particular `F` does not decide
whether `K_c >= atanh(a_n^{-1/n})` holds for any `n`, and no finite
extension of the checked range changes this conclusion structurally.

*Proof.*  Two explicit completions agree with every element of `F`:

(i) the min-plus closure `s_n = a_n` for `n <= 32`,
`s_n = min_{1<=m<n} s_m s_{n-m}` for `n > 32`.  By induction every pair is
submultiplicative: for `m+n <= 32` it is the verified data; for `m+n > 32`,
`s_{m+n}` is by construction `<=` every split product, in particular
`s_m s_n`.  So `F` is consistent with full submultiplicativity, hence with
`limsup s_n^{1/n} <= s_32^{1/32} < 4.7279` -- the scenario in which the
route's arithmetic works.  ([COMPUTATION]: verified exhaustively for all
pairs with `m+n <= 64`; `s_33 = 23301326963237632125732 = a_1 a_32`.)

(ii) the divergent completion `d_n = a_n` for `n <= 32`, `d_n = 6^n` for
`n > 32`.  It matches every verified value, satisfies every checked pair
inequality (each involves only indices `<= 32`), respects the envelope
`6^n <= (1000(n+1)^3)^{n+2}`, and has `limsup d_n^{1/n} = 6`, so the series
radius would be `1/6` and the best conclusion `K_c >= atanh(1/6)`, weaker
than the incumbent.  Its first submultiplicativity violation sits at
`(m,n) = (1,32)`, i.e. at order `33`, invisible to every check through
order `32`.

Since both completions are consistent with `F` and give contradictory
endpoint conclusions, `F` entails neither.  Checking through any order `N`
merely moves the first invisible violation to `N+1`.  QED.

This is the precise sense in which the assigned rule "finite evidence must
never be promoted to an all-size theorem" is *forced* here, not merely
prudent: Theorem 5's conclusion genuinely needs all pairs, and the data can
never supply them.

## 8. Why the SAW concatenation proof does not transfer

[COMPUTATION]/[LEMMA 7] The classical proof of `c_{m+n} <= c_m c_n` splits
a walk at a vertex and translates.  At coefficient level the objects are
clusters *dressed by `g`*, and the dressing is not multiplicative across a
split:

* **[LEMMA 7] (factorization at long distance).**  If
  `dist(U_1, U_2) > K` then `g(U_1 ∪ U_2) = g(U_1) g(U_2) mod v^{K+1}`.
  *Proof sketch (fully local):* in Lemma 2 every `B_1` has at most `K`
  edges, so each component meets exactly one `U_i`, no component can meet
  both, and no chain built from `U_1` within total size `K` can reach
  `U_2`'s removals; the sum over `B_1` factorizes as a product of the two
  one-sided sums, and induction over the truncation order factorizes the
  inner dressings.  QED.
* **[COMPUTATION] (failure at contact).**  At order 4 the dressing of a
  single edge is `-20` (twenty 4-cycles meet it); the dressing of the
  concatenated collinear pair `{0, e_1} ∪ {e_1, 2e_1}` is `-28`, not
  `-40 = -20 - 20`; a distance-10 pair is dressed by exactly `-40`.

So any coefficient-level concatenation injection must repair `O(1)`-sized
boundary terms at the splice, and those terms carry *both signs* (`g` is
signed).  This blocks the natural transfer of the (A.1)-style loop-erasure
injection, which controls functions at fixed `v >= 0` but not coefficients.
[UNRESOLVED] This is an obstruction to the natural proof strategy, not an
impossibility theorem for all conceivable proofs of G1.

## 9. What is and is not established

| claim | classification |
|---|---|
| Lemmas 1-3 (cluster identity, recursion, locality) with finite-graph brute-force validation | [THEOREM]/[LEMMA] with [COMPUTATION] controls |
| `a_0..a_10` exact from finite graphs only | [COMPUTATION], theorem-grade; equals the published prefix ([EXTERNAL] cross-check) |
| supermultiplicativity false; smallest counterexample `(1,1)`, `30 < 36`; all anchored variants fail at checked anchors | [THEOREM] |
| Fekete direction table; ratio-ceiling route valueless; Lemma 5 | [THEOREM] |
| box coefficient negativity and non-monotonicity counterexamples | [COMPUTATION], exact, hand-checkable |
| conditional bridge Lemma 4 and route Theorem 5 | [LEMMA]/[THEOREM], conditional on G0, G1 and P1'; uses Lieb (1980) [EXTERNAL THEOREM] |
| G0 (all-order positivity), G1 (submultiplicativity), P1' (absolute domination), log-concavity | [CONJECTURE]: every applicable finite check passes; nothing promoted |
| finite-verification no-go (two completions), envelope Lemma 6 | [THEOREM] |
| conditional targets `t_10, t_25, t_32` and crossing order 21 | [COMPUTATION], certified enclosures; **no bound claimed** |
| `a_n <= c_n` for all `n`, exact asymptotics | [UNRESOLVED] (finite range only) |

## Reproduction

From the repository root:

    .venv/bin/python experiments/e218_susceptibility_coeffs.py   # ~47 min elapsed on a loaded host
    .venv/bin/python experiments/e219_coeff_inequality.py        # seconds
    .venv/bin/python experiments/e220_coeff_bound.py             # seconds
    .venv/bin/python tests/test_susceptibility_coeffs.py         # standalone verifier, a few minutes
