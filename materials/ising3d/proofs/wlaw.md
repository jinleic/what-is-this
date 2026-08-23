# The symmetric-square law: an exact traceless-exterior model for the defect, the balanced polarization's rank and kernel, and closure of the natural-map route

Artifacts: `experiments/e191_wlaw_structure.py`,
`experiments/e192_wlaw_module.py`,
`experiments/e193_wlaw_certificate.py`,
`tests/test_wlaw.py`,
`results/algebra_growth/wlaw.json`.

Reproduce (from the repository root):

```sh
.venv/bin/python experiments/e191_wlaw_structure.py
.venv/bin/python experiments/e192_wlaw_module.py
.venv/bin/python experiments/e193_wlaw_certificate.py   # writes results/algebra_growth/wlaw.json
.venv/bin/python tests/test_wlaw.py                     # independent verifier, no producer imports
```

## 0. Setting, target, and verdict

Notation is that of `proofs/ladder_w8.md` and `proofs/ladder_alll_proof.md`:
open `2xL` ladder, vacuum `psi`, generators `A_L`, `B_L`, cyclic space
`U_L = alg(A_L,B_L) psi`, invariant ceiling `K_L`, annihilator
`W_L = U_L^perp` inside `K_L`, hard-core pieces `D` (pair creation), `F`
(hop), `D+` (pair annihilation), all in `alg(A_L,B_L)`; `tau` = leg swap,
`rho` = rung reversal.  The wave-18 conjecture (H459) is

```
rank U_L^(2m) = T(L,m) - delta(L,m),
T(L,m)     = C(L,m)(C(L,m)+1)/2,
delta(L,m) = 1 iff L = 0 mod 4 and m = L/2, else 0.
```

This note attacks the conjecture structurally.  Its verdict, in one line:
the conjectured right-hand side is exactly the dimension of the
top-wedge-traceless symmetric square `Sym^2_0(Lambda^m Q^L)`, with the
`delta` defect forced by invariant theory (Theorems 1-2); but the two
*natural* ways to realize that model inside the ladder — the balanced-leg
polarization, and any reflection-equivariant map in a signed subset basis —
are both closed off by exact counter-certificates (Theorems 3-5).  The law
itself remains `[UNRESOLVED]`; nothing here is an all-`L` rank claim.

Throughout, `V = Q^L` is the rung space with basis `e_0..e_{L-1}`,
`n = C(L,m)`, and `P_m` is the *permutation* module with basis
`{u_S : S an m-subset of rungs}`.  `Lambda^m V` has basis
`{e_S = e_{s_1} ^ ... ^ e_{s_m}, s_1 < ... < s_m}`.

---

## 1. [THEOREM] Dimension identity: the candidate space has exactly the conjectured dimension

Define the **top-wedge trace** `tr : Sym^2(Lambda^m V) -> Q` by

```
tr(x . y) = coefficient of e_0 ^ e_1 ^ ... ^ e_{L-1} in x ^ y + y ^ x
```

(the symmetrization of `b(x,y) := [x ^ y]_top`, which makes `tr` well
defined on the symmetric product for every `m`; identically zero unless
`2m = L`), and the **candidate space** `C(L,m) := ker(tr) <= Sym^2(Lambda^m V)`.

**Theorem 1.** For every `L >= 1` and `0 <= m <= L`,

```
dim C(L,m) = T(L,m) - delta(L,m).
```

*Proof.*  `dim Sym^2(Lambda^m V) = T(L,m)`.  If `2m != L` then `x ^ y` has
degree `2m != L` and `tr = 0`, so the corank is `0 = delta`.  If `2m = L`,
then for `x, y in Lambda^m V`,

```
y ^ x = (-1)^(m^2) x ^ y = (-1)^m x ^ y,
```

so `b` is symmetric for `m` even and antisymmetric for `m` odd.  For `m`
odd the symmetrized functional `tr(x.y) = b(x,y) + b(y,x) = 0`, corank
`0 = delta`.  For `m` even, `tr(e_S . e_{S^c}) = 2 b(e_S, e_{S^c}) != 0`,
so `tr != 0` and the corank is `1 = delta`.  (Here `4 | L` and `m = L/2`
is exactly `2m = L` with `m` even.)  ∎

**Lemma 1.1 (totals).** Summing over sectors,
`sum_m (T(L,m) - delta(L,m)) = (C(2L,L) + 2^L)/2 - [4|L]`, by the
Vandermonde identity `sum_m C(L,m)^2 = C(2L,L)` and
`sum_m C(L,m) = 2^L`.  So the sector law and the total law
(`dim U_L`, hence `W_L = dim K_L - dim U_L`) are equivalent packages.
Verified exactly for `L = 2..12` (`e193`/scratch; controls below).

## 2. [THEOREM] Defect dichotomy: `delta` is forced by invariant theory

**Theorem 2.** Let `0 < m < L`.  The space of `SL_L(Q)`-invariant linear
functionals on `Sym^2(Lambda^m V)` has dimension exactly
`delta(L,m) = [2m = L and m even]`, and when it is nonzero it is spanned
by `tr`.  Consequently, for interior sectors an `SL_L`-canonical
corank-one cut of the symmetric square exists **iff** `L = 0 mod 4` and
`m = L/2`, and then `ker(tr)` is the unique such cut.

*Proof.*  Functionals on `Sym^2 X` are symmetric bilinear forms on `X`
(char 0).  Let `beta` be an invariant bilinear form on `X = Lambda^m V`.

*Torus support.*  For `t = diag(t_0..t_{L-1})` with `prod t_i = 1`,
invariance at the pair `(e_S, e_T)` forces
`prod_{i in S} t_i * prod_{i in T} t_i = 1` whenever `beta(e_S,e_T) != 0`.
Taking `t = (s, 1/s, 1, ..., 1)` and permutations, this forces the
multiplicity vector `1_S + 1_T` to be a constant vector `c * (1,..,1)`
with `c in {0,1,2}`.  For `0 < m < L`: `c = 0` and `c = 2` are impossible
(they force `m = 0` or `m = L`), so `c = 1`, i.e. `T = S^c` and `2m = L`.
Hence `beta = 0` whenever `2m != L`, and otherwise `beta` is supported on
complementary pairs: `beta(e_S, e_{S^c}) = lambda_S`, all other entries
zero.

*Unipotent transport.*  Fix `i != j`, `g = I + s E_{ji} in SL_L(Q)`.  On
`Lambda^m`, `g e_S = e_S + s * sigma(S,i,j) e_{S'}` when `i in S`,
`j notin S` (`S' = S \ i + j`, `sigma` the resorting sign), and
`g e_S = e_S` otherwise.  Apply invariance to the pair
`(e_S, e_{S'^c})` where `i in S`, `j notin S`.  The order-`s^0` and
order-`s^2` terms pair non-complementary subsets and vanish by the
support result; the order-`s` term gives

```
sigma(S,i,j) * lambda_{S'} + sigma(S'^c,i,j) * lambda_S = 0.
```

Since any `m`-subset is reachable from any other by such single-element
moves, all `lambda_S` are determined by one of them: the invariant-form
space has dimension at most one.  Finally `b` itself is invariant
(`g` acts on `Lambda^L` by `det g = 1`), and `b` is symmetric iff `m`
is even (Theorem 1's sign computation), so the symmetric invariant space
has dimension `[2m = L][m even] = delta`.  ∎

**Endpoint remark (scope).**  At `m in {0, L}` the `SL_L`-module
`Lambda^m V` is trivial, every functional is invariant, and Theorem 2's
dichotomy does not apply; there `tr = 0`, `T(L,m) = 1 = T - delta`, and
the certified endpoint sector ranks are trivially `1`.  The dichotomy is
an interior-sector statement and is stated as such.

**What this explains.**  The long-noted "exceptional" corrections
`W_4 = 2` and `rank U_8^(8) = T(8,4) - 1` sit exactly at the only
interior sectors where an invariant hyperplane exists at all.  Under the
(unproved) identification `U_L^(2m) ~ C(L,m)`, the defect is not an
anomaly: it is the removal of the unique invariant trace direction, and
it can only fire at `4 | L`, `m = L/2` — nowhere else is a canonical
corank-one cut even available.  `[CONJECTURE]` for the identification;
`[THEOREM]` for everything about the model itself.

---

## 3. [THEOREM] The canonical bilinear map: exact image, kernel, and rank

The canonical symmetric-square map suggested by the sector data is the
**balanced polarization**: for `m`-subsets `S, T` let `x_{S,T}` be the
configuration with top-leg occupation `S` and bottom-leg occupation `T`,
and define

```
pi : Sym^2(P_m) -> K_L,     pi(u_S . u_T) = sum of the <tau,rho>-orbit of x_{S,T}.
```

(`tau x_{S,T} = x_{T,S}` makes this well defined on unordered pairs; the
orbit sum lands in `K_L^(2m)`.)  Let `f = f(L,m)` be the number of
reversal-fixed `m`-subsets; from the cycle structure of `rho` on rungs
(`floor(L/2)` two-cycles plus `L mod 2` fixed points),

```
f(L,m) = [x^m] (1+x)^(L mod 2) (1+x^2)^floor(L/2)
       = C(floor(L/2), floor(m/2))   (L odd),
         C(L/2, m/2) [m even]        (L even).
```

**Theorem 3.**  For every `L >= 2`, `0 <= m <= L`:

```
ker pi  = the rho-odd part of Sym^2(P_m),   dim ker pi  = (n^2 - f^2)/4,
rank pi = (n^2 + 2n + f^2)/4,
```

where `rho` acts on unordered pairs by `{S,T} -> {rho S, rho T}`.

*Proof.*  `pi({S,T}) = pi({rho S, rho T})` because the two configurations
lie in the same group orbit; hence the odd part is killed and
`rank pi <= dim` of the even part.  Burnside on the two-element group
`{1, rho}` acting on unordered pairs gives even-part dimension
`(T(L,m) + (f^2+n)/2)/2 = (n^2+2n+f^2)/4` (the number of `rho`-fixed
unordered pairs is `(f^2+n)/2`: `f(f+1)/2` pointwise-fixed pairs plus
`(n-f)/2` swapped pairs `{S, rho S}`).  Conversely the map
`(S,T) -> {S,T}` induces a bijection between `<tau,rho>`-orbits of
balanced `(m,m)` configurations and `rho`-orbits of unordered pairs, and
orbit sums of distinct configuration orbits have disjoint supports, hence
are linearly independent; so the image has dimension equal to the orbit
count `(n^2+2n+f^2)/4`, and the kernel is exactly the odd part.  ∎

`[COMPUTATION]` All three descriptions (pair enumeration, direct
configuration-orbit census over all `4^L` configurations, closed forms)
agree on all 52 rows `L = 2..9`, `0 <= m <= L`
(`wlaw.json: structure_records_exact_L2_9`, checks C1; independently
re-derived in `tests/test_wlaw.py` V2/V6 with a generating-function `f`
and a fresh census).

---

## 4. [THEOREM] The balanced image is not a `D/F/D+`-submodule, for any `L >= 2`

Let `B_L = span of <tau,rho>-orbit sums of leg-balanced configurations`
(`= sum_m image(pi_m)`).  All three hard-core pieces leave it:

**Theorem 4.**  For every `L >= 2`, none of `D B_L`, `F B_L`, `D+ B_L`
is contained in `B_L + (unbalanced components zero)`; concretely each of
`D`, `F`, `D+` maps an explicit vector of `B_L` to a vector with a
nonzero leg-unbalanced component:

* `D`: the vacuum `psi in B_L^(0)`; `D psi` contains the top-leg dimer
  configuration (top rungs `{0,1}`) with coefficient `1`.
* `F`: the orbit sum of `x` (top `{0}`, bottom `{1}`) lies in `B_L^(2)`;
  `F` applied to it contains the configuration top `{0,1}`, bottom empty,
  with coefficient `>= 2` (exactly `2` on the certified range `L = 2..9`).
* `D+`: the orbit sum of the `2x2` block (top `{0,1}` = bottom) lies in
  `B_L^(4)`; `D+` applied to it contains the configuration bottom `{0,1}`,
  top empty, with coefficient `>= 1`.

*Proof.*  Each of `D`, `F`, `D+` is a sum over ladder edges of matrices
with entries in `{0,1}` in the configuration basis, and each witness
vector has all coefficients `+1`; therefore no cancellation can occur and
a single explicit contribution proves the component nonzero.  For `D`:
the top edge `(0,1)` creates the top dimer from the vacuum.  For `F`: the
rung-1 hop moves the bottom particle of `x` up, and the rung-0 hop does
the same in `tau x`; both target top `{0,1}`.  For `D+`: annihilating the
top edge `(0,1)` of the block leaves bottom `{0,1}`.  Every target has
leg counts `(2,0)` or `(0,2)`.  ∎

`[COMPUTATION]` All 24 rows `(L, piece)`, `L = 2..9`, recomputed exactly
with full sparse vectors: sources balanced and `<tau,rho>`-invariant,
outputs graded exactly by `+2 / 0 / -2` particles, leakage nonzero and
itself `<tau,rho>`-invariant (`wlaw.json: module_records_exact_L2_9`,
checks C5; verifier V3 re-derives every row in a different site
encoding).

**Exceptional-row refinement.**  At the defect row `(L,m) = (4,2)` one
could hope the leak lives only in the trace direction that the candidate
`C(4,2)` removes.  It does not: the `D+` witness corresponds to the
symmetric-square element `u_S.u_S + u_{rho S}.u_{rho S}` (`S = {0,1}`),
whose top-wedge trace is `b(e_S,e_S) + b(e_{rho S},e_{rho S}) = 0`; the
witness lies in the traceless candidate and still leaks (check C6).

**Consequence.**  `pi` is not a module map onto its image for the
`{D,F,D+}` action in any sense, so the balanced polarization cannot by
itself transport `Sym^2`-structure to `U_L^(2m)`.  This is a rigorous
closure of the most obvious route to the law.

---

## 5. [THEOREM] No natural reflection-equivariant map can realize the law — with an unconditional counter-certificate at `(3,1)`

Call a linear map `Phi : Sym^2(P_m) -> K_L^(2m)` **natural** if it
intertwines rung reversal, where reversal acts on `P_m` by *any* signed
subset permutation `rho^ u_S = eps_S u_{rho S}` (`eps_S in {+1,-1}`,
`(rho^)^2 = 1`) and acts on `K_L` as the identity (it does: `K_L` is the
`<tau,rho>`-invariant subspace).  Global sign twists, and in particular
both the permutation module `P_m` and the exterior module `Lambda^m V`
with its signed reversal, are instances.

**Theorem 5.**  For every gauge `eps`:

```
rank Phi <= dim (Sym^2 P_m)^(Sym^2 rho^) <= (n^2 + 2n + f^2)/4,
```

and for every interior sector `0 < m < L`,

```
(n^2 + 2n + f^2)/4  <  T(L,m) - delta(L,m).
```

*Proof.*  Naturality forces `Phi o Sym^2(rho^) = Phi`, so `Phi` kills the
`Sym^2(rho^)`-odd part and `rank Phi <= dim` even part.  On a two-element
pair-orbit `{{S,T},{rho S,rho T}}` the even part contributes one
dimension regardless of signs.  On a fixed pair of swap type `{S, rho S}`
the diagonal coefficient is `eps_S eps_{rho S} = +1` (forced by
`(rho^)^2 = 1`).  On a fixed pair of pointwise type (`rho S = S`,
`rho T = T`) it is `eps_S eps_T in {+1,-1}`, contributing at most one.
Hence the even part is at most the orbit count `(n^2+2n+f^2)/4`, with
equality at `eps = 1`.  (A projective lift `(rho^)^2 = -1` needs
`eps_S eps_{rho S} = -1` for all `S`, impossible when `f > 0`, and gives
a strictly smaller even part when `f = 0`.)

For the gap: interior `m` gives a non-fixed subset (`{0..m-1}` is moved
by `rho`), and non-fixed subsets pair up, so `f <= n-2`; hence

```
T - delta - (n^2+2n+f^2)/4 = (n^2 - f^2)/4 - delta >= n - 1 - delta > 0,
```

using `n >= 2` when `delta = 0` and `n = C(L, L/2) >= 6` on defect rows.  ∎

**Corollary 5.1 (route closure, unconditional at the earliest row).**
The certified rank `rank U_3^(2) = 6 = T(3,1)` (modular rank 6 at
`q = 999983` is a rational lower bound, and `dim K_3^(2) = 6` is the
ceiling, so the rank is exactly 6 over `Q`).  Every natural map at
`(L,m) = (3,1)` has rank at most `(9+6+1)/4 = 4 < 6`.  Therefore **no
natural reflection-equivariant symmetric-square map can produce the
certified sector space even at the earliest controlled row**, and by
Theorem 5's strict gap the same failure recurs at every interior sector
of every `L` *if* the law's value is the true rank there.  This closes
the natural symmetric-square route to the law; any proof must use a
non-monomial intertwiner (reversal acting on the model by something
other than a signed subset permutation) or avoid an explicit intertwiner
altogether.

`[COMPUTATION]` checks C3 (all 36 interior rows: nonzero forced kernel
and strict ceiling), C4 (the `(3,1)` counter-certificate), verifier V5
(re-derivation plus the formula-level strict gap for all `L <= 30`).

---

## 6. [COMPUTATION] Held-out controls: seven totals, 49 sector values, zero fit inputs

The candidate `T(L,m) - delta` is parameter-free and was constructed by
`e193` *before* loading any ladder artifact.  The prior certified rows
are used only as controls:

| L | certified sector ranks (`k = 0,2,...,2L`) | total | W_L | all `= T(L,m) - delta`? |
|---|---|---|---|---|
| 3 | 1, 6, 6, 1 | 14 | 0 | yes |
| 4 | 1, 10, 20, 10, 1 | 42 | 2 | yes (`delta` at `k=4`) |
| 5 | 1, 15, 55, 55, 15, 1 | 142 | 10 | yes |
| 6 | 1, 21, 120, 210, 120, 21, 1 | 494 | 66 | yes |
| 7 | 1, 28, 231, 630, 630, 231, 28, 1 | 1780 | 364 | yes |
| 8 | 1, 36, 406, 1596, 2484, 1596, 406, 36, 1 | 6562 | 1822 | yes (`delta` at `k=8`) |
| 9 | 1, 45, 666, 3570, 8001, 8001, 3570, 666, 45, 1 | 24566 | 8586 | yes |

Sources: `results/algebra_growth/ladder_w8_regression.json` (fresh
mirror closures `L = 3..8` at `q = 999983`) and
`results/ladder/w9_saturation.json` (wave-16 exact-Q sandwich), pinned by
sha256 in `wlaw.json: meta.source_artifacts`.  Sector `K`-dimensions are
additionally re-derived by a direct `4^L` orbit census in the verifier.
Zero of these 49 sector values or 7 totals entered the construction
(checks C7, C8; protocol block `control_protocol`).

---

## 7. Scope, honest negatives, counterexample discipline

* The rank law itself remains `[UNRESOLVED]`.  Neither inclusion of
  `rank U_L^(2m) = T(L,m) - delta` is proved for any `L >= 10`, and no
  such claim is made.  `W_10` was not computed (standing lead directive);
  the finite controls stop at `L = 9`.
* Theorems 1-2 are statements about the *model* `Sym^2_0(Lambda^m Q^L)`;
  they prove that the conjectured dimensions are canonical, not that
  `U_L^(2m)` realizes them.  If a future exact closure finds
  `rank U_L^(2m) != T(L,m) - delta` at some row, Theorems 1-5 all
  survive; only the conjectural identification dies.
* Theorem 5's no-go is scoped to *natural* maps (signed subset
  permutation lifts of reversal).  A non-monomial reversal action on the
  model, or a rank argument with no explicit intertwiner, is not
  excluded — that is exactly the remaining open route, and the theorem
  tells any successor front not to spend effort on monomial gauges.
* The endpoint sectors `m in {0, L}` are outside Theorem 2's dichotomy
  (trivial `SL_L`-modules); the dimension identity (Theorem 1) covers
  them, and their certified ranks are trivially 1.
* No counterexample to the law was found: the earliest-row certificate of
  Corollary 5.1 refutes the *route*, not the *law* (the certified rank 6
  equals `T(3,1)`; it is the natural maps that cannot reach it).
* Modular ranks cited from source artifacts are rational lower bounds;
  at `(3,1)` the bound is saturated by the ceiling `dim K_3^(2) = 6`, so
  the counter-certificate needs no unproved exactness.
* Resource honesty: producer `e193` runs in `< 1 s` CPU, peak RSS
  `~26 MB`; the verifier's heaviest step is the `4^9` census.  No
  closures were launched and nothing here touched `K_c`.

## 8. Claim ledger

| claim | status |
|---|---|
| `dim ker(tr) = T(L,m) - delta(L,m)` for all `0 <= m <= L` | `[THEOREM]` (Sec. 1) |
| interior invariant-functional dichotomy; `ker(tr)` unique canonical cut; defect location `4|L, m = L/2` forced | `[THEOREM]` (Sec. 2) |
| `sum_m (T - delta) = (C(2L,L)+2^L)/2 - [4|L]` | `[LEMMA]` (Sec. 1) |
| `ker pi` = reflection-odd part; `rank pi = (n^2+2n+f^2)/4` | `[THEOREM]` (Sec. 3) |
| balanced span not invariant under any of `D`, `F`, `D+`, every `L >= 2`; traceless refinement at `(4,2)` | `[THEOREM]` (Sec. 4) |
| natural-map ceiling `(n^2+2n+f^2)/4`, strict interior gap, unconditional `(3,1)` counter-certificate | `[THEOREM]` (Sec. 5) |
| 52 structure rows, 24 module rows, exact `L = 2..9` | `[COMPUTATION]` (checks C1-C6) |
| 7 totals + 49 sector controls match, zero fit inputs | `[COMPUTATION]` (checks C7-C8) |
| the symmetric-square rank law; `U_L^(2m) ~ Sym^2_0(Lambda^m Q^L)` | `[CONJECTURE]`, identification `[UNRESOLVED]` |
| any all-`L` value of `rank U_L^(2m)` or `W_L` for `L >= 10` | `[UNRESOLVED]`, not claimed |
