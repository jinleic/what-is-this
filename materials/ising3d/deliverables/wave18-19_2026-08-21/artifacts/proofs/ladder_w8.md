# A partition-family growth theorem for the two-generator ladder algebra, the all-L free-fermion no-go, and the symmetric-square annihilator law

Artifacts: `experiments/e172_ladder_w8_regression.py`,
`experiments/e173_ladder_nogo_family.py`,
`experiments/e174_ladder_gap_certificates.py`,
`tests/test_ladder_w8.py`,
`results/algebra_growth/ladder_w8_regression.json`,
`results/algebra_growth/ladder_w8_family.json`,
`results/algebra_growth/ladder_w8_certificates.json`,
`results/algebra_growth/ladder_w8.json` (consolidated envelope).

Reproduce:

```sh
PYTHONPATH=src .venv/bin/python experiments/e172_ladder_w8_regression.py
PYTHONPATH=src .venv/bin/python experiments/e173_ladder_nogo_family.py
PYTHONPATH=src .venv/bin/python experiments/e174_ladder_gap_certificates.py
PYTHONPATH=src .venv/bin/python experiments/e174_ladder_gap_certificates.py --assemble
PYTHONPATH=src .venv/bin/python tests/test_ladder_w8.py
```

**Status convention.** `[THEOREM]` is proved here or in a cited repository
note; `[COMPUTATION]` is an exact finite calculation (modular ranks are
one-sided: characteristic-zero lower bounds only); `[CONJECTURE]` is
unproved; `[UNRESOLVED]` is deliberately undecided.

**Naming note.** This front was dispatched under the file family
`ladder_w8.*`.  The datum `W_8 = 1822` was already certified in wave 15
(`results/ladder/l8_saturation.json`, H393) and `W_9 = 8586` in wave 16
(`proofs/ladder_l9.md`, H409); the wave-18 lead retargeted the front, and
`W_10` was deliberately not attempted (projected `>= 15x` the L=9 CPU with
peak RSS near or over the 4 GB cap; a wall datum decides nothing).

---

## 0. Setting

Open `2xL` ladder, rungs `0..L-1`, `X`-eigenbasis occupation
configurations, vacuum `psi = |+>^(2L)` (empty configuration);

```
A_L = sum_v X_v = 2L - 2N,     B_L = sum_{(u,v) in E(2xL)} Z_u Z_v,
g_L = <A_L, B_L>_Lie  (over Q),      U_L = alg(A_L, B_L) psi,
K_L = even-particle <tau,rho>-invariant space,  dim K_L = 2^{2L-3} + 3*2^{L-2},
W_L = (U_L)^perp inside K_L   (pairing kernel / annihilator).
```

`A_L` has distinct eigenvalues `2L - 2k` on distinct particle sectors, so
every sector projector `P_k` is a Lagrange polynomial in `A_L`, and the
hard-core pieces `D = P_{k+2} B P_k` (pair creation), `F = P_k B P_k`
(hop), `D+ = P_{k-2} B P_k` (pair annihilation) lie in the unital
associative algebra `alg(A_L, B_L)`
(`proofs/ladder_alll_proof.md` Theorem 1; `proofs/ladder_l9.md` sec. 2).
The evaluation bound `dim g_L >= dim U_L` is Theorem 2 of
`proofs/ladder_alll_proof.md`.

Certified input data (waves 14-16, two-sided exact-Q sandwich):

```
L               3     4     5     6     7      8      9
dim U_L         14    42    142   494   1780   6562   24566
dim K_L         14    44    152   560   2144   8384   33152
W_L             0     2     10    66    364    1822   8586
```

---

## 1. [THEOREM] The partition-family lower bound (headline result)

Index set (stated twice to remove any ambiguity about zero parts):

* **Sequence form.**  `I(L)` is the set of weakly decreasing integer
  sequences `lambda = (a_1 >= a_2 >= ... >= a_m >= 0)` with `m >= 1` and
  `|lambda| + 2m <= L`, where `|lambda| = sum a_i` and **`m` counts every
  entry, zero entries included**.  Sequences of different lengths are
  different indices even when they agree after deleting zeros: a zero
  entry is an *unspread dimer* — one extra particle pair in a higher
  sector — not a formal padding symbol.  So `(2)` and `(2, 0)` are
  distinct members of `I(6)`.
* **Bijective partition form.**  Mapping `lambda` to
  `nu := (a_1 + 2, a_2 + 2, ..., a_m + 2)` is a bijection from `I(L)` onto
  the set of nonempty partitions with **all parts `>= 2` and total
  `<= L`** (with `m = len(nu)` recovered).  Hence, unambiguously,

```
Q(L) := #I(L) = #{ nonempty partitions nu, every part >= 2, |nu| <= L }.
```

First values (`L = 2..12`): `1, 2, 4, 6, 10, 14, 21, 29, 41, 55, 76`.
These are strictly larger than the count of standard partitions
satisfying `|mu| + 2 len(mu) <= L` with positive parts only
(`1, 2, 3, 4, 6, 8, 11, 15, 20, 26, 35`): the machine-verified family
ranks below (`21, 29, 41` at `L = 8, 9, 10`) match the first count and
exceed the second, settling which object the theorem is about.

For `lambda in I(L)` define the walk vector (read right to left; the
`j`-th created dimer is followed by the phase `F^{a_j}`, so the first
dimer receives the largest phase):

```
w(lambda) = F^{a_m} D  F^{a_{m-1}} D  ...  F^{a_2} D  F^{a_1} D  psi   in U_L.
```

**Theorem 1.1.**  For every `L >= 2`:

```
dim_Q g_L  >=  dim_Q U_L  >=  Q(L).
```

*Proof.*

*(P) Positivity.*  All matrix entries of `D` and `F` in the configuration
basis are 0 or 1, so the coefficient of a configuration `C` in
`w(lambda') psi` equals the number of *histories* — sequences of elementary
events (a `D`-event creates a pair on an empty edge, an `F`-event moves one
particle across an edge to an empty site) — ending at `C`.  Worldlines
through a history are well defined (each `F`-event moves one identified
particle), each particle descends from exactly one `D`-event, and a
particle created at site `s` ending at site `t` consumes at least
`graphdist(s, t)` `F`-events.

*(W) Witnesses.*  For `lambda in I(L)` put `d_j := a_j + 1` and place `m`
particle pairs on the TOP leg, left to right, pair `j` occupying rungs
`x_j` and `x_j + d_j`, with `x_1 = 0` and `x_{j+1} = x_j + d_j + 1`
(particle gap exactly one rung).  The rightmost rung used is
`|lambda| + 2m - 1 <= L - 1`, so `C(lambda)` exists; `lambda` is
recoverable from `C(lambda)` (consecutive-pair distances), so
`lambda -> C(lambda)` is injective.

*(C) Cost lemma.*  Let `h` be a history of any word with `m'` `D`-events
and `T'` `F`-events ending at `C(lambda)`, whose particles sit at top-leg
rungs `p_1 < ... < p_{2m}`.  Particle number forces `m' = m`.  Let `Pi(h)`
be the worldline pairing.  A pair covering top-leg targets `(q, q')`,
`q < q'`, costs at least `q' - q - 1` if its `D`-event used a top
leg-edge (minimize `|q - y| + |q' - (y+1)|` over the creation rung `y`),
and at least `(q' - q - 1) + 2` if it used a bottom leg-edge (two rung
crossings) or a rung edge (`|y - q| + |y - q'| + 1 >= q' - q + 1`).  Hence

```
T'  >=  sum_{pairs in Pi(h)} (dist - 1)  +  2 #(non-top-created pairs).
```

*(M) Pairing minimality.*  For points `p_1 < ... < p_{2m}` on a line,
represent a pairing by its interval set and count, for each elementary gap
`[p_i, p_{i+1}]`, the covering number `c_i`.  Parity forces
`c_i ≡ i (mod 2)`.  The consecutive pairing `(p_1p_2)(p_3p_4)...` realizes
the minimum `c_i = (i mod 2)`; any other pairing has some even gap with
`c_i >= 2` or some odd gap with `c_i >= 3`, an excess of at least
`2 * (gap length) >= 2`.  The consecutive value on `C(lambda)` is
`sum_j d_j`.

*(B) Budget necessity.*  The witness cost is
`T := sum_j (d_j - 1) = |lambda|`.  If `T' < T` the coefficient of
`C(lambda)` is zero.  If `T' = T`, lemmas (C) + (M) force: every dimer
top-leg-created, the worldline pairing consecutive, and every pair's moves
exactly `dist - 1` (no wasted moves).

*(S) Scheduling (Hall) necessity.*  In `w(lambda')` the pair created by
the `j`-th `D`-event moves only during phases `j..m`, of total budget
`sum_{i >= j} a'_i`.  Under budget-exactness the geometric pairs' costs
are the multiset `lambda` (sorted decreasingly), and *whatever* the
assignment of geometric pairs to creation slots, the pairs in slots
`j..m` cost at least the sum of the `m - j + 1` smallest elements, i.e.
`sum_{i >= j} lambda_i`.  So a nonzero coefficient of `C(lambda)` in
`w(lambda') psi` with `|lambda'| = |lambda|` requires
`sum_{i >= j} lambda_i <= sum_{i >= j} lambda'_i` for every `j`, which at
equal totals is exactly dominance `lambda |> lambda'`.

*(D) Diagonal.*  For `lambda' = lambda`, the history that creates dimer
`j` at rungs `(x_j, x_j + 1)` at the `j`-th `D`-event and spends the whole
phase `F^{a_j}` moving its right particle to rung `x_j + d_j` is
admissible (territories disjoint, creation edges empty when used), so the
diagonal coefficient is a positive integer.

*(T) Triangularity.*  Order `I(L)` by `(m, |lambda|,` any linear extension
of dominance`)`.  A witness `C(lambda)` receives zero coefficient from
every strictly earlier word: a different `m` is impossible (particle
number), a smaller total is impossible (budget), and at equal totals the
dominance necessity applies.  The witness-coefficient matrix is
block-lower-triangular with strictly positive diagonal over `Z`, so the
`Q(L)` integral vectors `w(lambda) psi` are linearly independent over `Q`.
Both inequalities follow.  `QED`

**[COMPUTATION] Exhaustive instance verification.**  At `L = 8, 9, 10` the
full witness-coefficient matrices were computed exactly
(`ladder_w8_family.json` and independently in the standalone test):
diagonal positive, sector/budget block structure exact, dominance
necessity exact, and the family rank over `F_p` equals
`Q(L) = 21, 29, 41` — the constructed vectors realize the full rank, not
merely the proof's lower bound.

**Remark (why this evades the word-family obstructions).**  The certified
obstructions of `proofs/ladder_all_l.md` and `proofs/ladder_nonstat.md`
concern literal bracket-word families and leading-term transfers.  The
present family lives in the state space: walk positivity plus
transport-cost accounting replaces leading-term bookkeeping.  The price is
a subexponential count (partition growth), far below the certified
`~3.7^L` evaluation dimensions — but far above every polynomial threshold,
which is all the no-go needs.

### 1.2 [THEOREM] Threshold arithmetic

Let `tau(L) := 2L(4L-1) = 8L^2 - 2L`.

* Exact DP (`ladder_w8_family.json`; the DP is over the sequence form and
  agrees with direct enumeration of the parts `>= 2` partition form):
  `Q(L) > tau(L) + 1` for every `32 <= L <= 300`.  Crossover detail:

  ```
  L      30     31     32     33     35     40
  Q      5603   6841   8348   10142  14882  37337
  tau+1  7141   7627   8129   8647   9731   12721
  ```

  `Q(31) = 6841 <= 7627` and `Q(32) = 8348 > 8129`: the family alone
  certifies exactly `L >= 32`.  (Counting only positive-part partitions
  with `|mu| + 2 len(mu) <= L` would instead cross near `L = 41`; that is
  a strictly smaller index set than `I(L)` — see the two counts in
  Section 1 — and is *not* the object the theorem uses.)
* Quartic tail: partitions into at most 4 positive parts of total
  `<= L - 8`, zero-padded to `m = 4` exactly (equivalently: `nu` with
  exactly 4 parts, all `>= 2`, `|nu| <= L`), are valid indices, so
  `Q(L) >= C(L-4, 4)/24 >= (L-7)^4/576`.  For `L >= 301`,
  `(L-7)^4/576 > 8L^2` (machine-checked for `301 <= L <= 5000`; for
  `L > 5000`, `(L-7)^4/576 - 8L^2` is a quartic-dominant polynomial with
  positive leading coefficient, increasing there).

**Corollary 1.3.**  `dim_Q g_L >= Q(L) > 8L^2 - 2L + 1` for every
`L >= 32`, unconditionally.

---

## 2. The all-`L` free-fermion no-go for the two-generator pair

**Lemma 2.1 (quadratic ceiling).**  Suppose there exist an invertible `V`
on `(C^2)^{tensor 2L}` and Majorana operators `c_1, ..., c_{4L}`
(`{c_i, c_j} = 2 delta_{ij}`) with `V A_L V^{-1}` and `V B_L V^{-1}` in
`span_C{ i c_j c_k : j < k } + C 1`.  Then `dim_Q g_L <= 2L(4L-1)`.

*Proof.*  `A_L`, `B_L` are traceless (sums of non-identity Pauli strings)
while `tr(c_j c_k) = 0` for `j != k` and `tr 1 = 2^{2L} != 0`; taking
traces kills the identity components, so both conjugated generators lie in
`Q_span := span{i c_j c_k}`, a Lie algebra isomorphic to `so(4L)` of
dimension `2L(4L-1)`.  The generated Lie algebra stays inside, and a
`Q`-basis of the integral algebra `g_L` remains `C`-independent.  `QED`

**Theorem 2.2 (no-go, assembled; seam table).**  For every `L >= 3`,
`dim_Q g_L > 2L(4L-1)`; hence the pair `(A_L, B_L)` on the open `2xL`
ladder is **not** simultaneously Majorana-quadratic in any basis of the
physical `2^{2L}`-dimensional space.  The assembly, with one row per
range and its certificate type (no gap between 3 and 32):

| `L` | bound on `dim g_L` | vs `2L(4L-1)` | certificate type |
|---|---|---|---|
| 3 | `263` | `> 66` | exact closure dimension `[THEOREM]` (`notes/algebra_growth.md`, H070) |
| 4 | `2952` | `> 120` | exact closure dimension `[THEOREM]` (H255, `proofs/char0_complete_2x4.md`) |
| 5 | `191` | `> 190` | Pauli-word one-prime certificate `[COMPUTATION]` (`e174`, `p = 2147483647`) |
| 6 | `494` | `> 276` | `K_L - W_L` exact-Q sandwich (wave 14) |
| 7 | `1780` | `> 378` | `K_L - W_L` exact-Q sandwich (wave 14) |
| 8 | `6562` | `> 496` | `K_L - W_L` exact-Q sandwich (wave 15, H393) |
| 9 | `24566` | `> 630` | `K_L - W_L` exact-Q sandwich (wave 16, H409) |
| 10..31 | `>= 8L^2 - 2L + 1` per row | `>` by 1 | sector-capped state one-prime certificates `[COMPUTATION]` (`e174`, `q = 999983`), per-`L` rows in `ladder_w8_certificates.json` |
| `>= 32` | `Q(L)` | `>` (Cor. 1.3) | partition-family `[THEOREM]` (Sec. 1) |

The `L = 10..31` certificates store only exact integral images of `psi`
under sector-filtered words `P_k B P_{k'}` (all in `alg(A_L, B_L)`)
restricted to particle sectors `0, 2, 4`; a modular rank of an integral
family is a rational lower rank, and the evaluation bound finishes the
row.  If any per-`L` run had ended at its process-time budget the row
would be recorded NON-DECISIVE and excluded from the claim; the shipped
artifact contains no such row (the standalone test re-derives every row
independently and fails loudly otherwise).

**Honest exception `L = 2`.**  `dim g_2 = 11 <= 28 = dim so(8)`: the
dimension criterion is inconclusive, and no no-go is claimed at `L = 2`.
The `2x2` ladder is the 4-cycle transverse-field Ising chain (ring DLA
dimension `3n - 1 = 11`, machine-verified in
`results/algebra_growth/profiles.json`), i.e. the known integrable case —
consistent with the criterion's silence.

**Contrast with the local-term algebra (different object, different
verdict).**  The wave-18 sibling front treats the LOCAL-term algebra
`<X_i, Z_iZ_j>` of the same lattice, which already at `L = 2` has
dimension `56 > 28`, so its no-go includes `L = 2`.  No contradiction:
`dim <A_2, B_2> = 11` (two global sums) versus `dim <X_i, Z_iZ_j> = 56`
(all local terms); at `2x3` the pair is `263` versus `1056`.  This note
claims nothing about the local-term algebra.

**Scope.**  Theorem 2.2 excludes simultaneous quadraticity of the two
generators in `4L` Majoranas on the physical space.  It does not exclude
free-fermion descriptions on enlarged spaces or after nonlinear changes of
variables, and at `L = 2` it excludes nothing.

---

## 3. [COMPUTATION] Regression control

A fresh, self-contained sector-pure mirror closure (independent code, same
proved theorems: sector-pure closure and particle-hole mirror,
`proofs/ladder_l9.md` sec. 2-3) at `q = 999983` reproduces every certified
row at `L = 3..8`: cyclic `14, 42, 142, 494, 1780, 6562`, annihilator
`0, 2, 10, 66, 364, 1822`, and every per-sector deficit split, in 67 s CPU
(`ladder_w8_regression.json`).  The standalone test repeats this with a
different prime (`999979`) and a different pivot convention.  For `L = 9`
the wave-16 artifacts are re-checked (basis sha256, sandwich
`24566 + 8586 = 33152`, low-sector ranks `{1, 45, 666, 3570, 8001}`,
sector palindrome); the 2.5 ks closure was not re-run (lead directive).

---

## 4. [CONJECTURE] The symmetric-square annihilator law

**Identification.**  The stored `L = 9` low-sector cyclic ranks are
exactly the symmetric-square triangle
`T(L,m) := C(L,m)(C(L,m)+1)/2` (A107105):

```
rank U^{(2m)}_9 = 1, 45, 666, 3570, 8001 = T(9, m),   m = 0..4.
```

**Law.**  For all `L >= 3`, `0 <= m <= L`:

```
rank U^{(2m)}_L = T(L, m) - delta,    delta = 1  iff  L ≡ 0 (mod 4) and m = L/2,
dim U_L = (C(2L,L) + 2^L)/2 - [4 | L],
W_L     = 2^{2L-3} + 3*2^{L-2} - (C(2L,L) + 2^L)/2 + [4 | L].
```

**Epistemic position.**  The law is `[CONJECTURE]`, but it is not a curve
fit.  It was *identified* from the `L = 9` sector vector alone (a
structural match to the symmetric-square triangle), consuming **zero** of
the seven certified totals as fit input, and then *tested* against
everything else:

* all seven certified totals `dim U_L` and `W_L` at `L = 3..9` — exact;
* all certified per-sector deficit splits at `L = 4..9` (28 sector
  values), via `W^{(2m)} = dim K^{(2m)} - T(L,m) + delta` — exact;
* the `delta`-correction fires exactly at `L = 4` (explaining the
  long-noted "exceptional" `W_4 = 2`) and at `L = 8` (middle sector
  `2484 = T(8,4) - 1`), and nowhere else — exact.

Zero failures, zero free parameters.  Still unproved: the verdict is
*unfalsified*, never *confirmed*.

**Falsifiable prediction (next wave).**
`cyclic_10 = 92890`, `W_10 = 131840 - 92890 = 38950`.

**[THEOREM, conditional on the law] the `2^L` conjecture.**  If the law
holds for a given `L`, then
`dim g_L >= K_L - W_L = (C(2L,L)+2^L)/2 - [4|L] >= 2^L`
(since `C(2L,L) >= 2^L + 2` for `L >= 2`).  If it holds for every `L`, the
standing `D_{2L} >= 2^L` conjecture follows at every `L`, with the much
stronger asymptotics `dim U_L ~ 4^L / (2 sqrt(pi L))`.

**[REFUTED APPROACH — do not retry.]**  The wave-18 lead's suggested
route — prove `W_L <= (1 - c) K_L` for a constant `c > 0` — is dead if the
law holds, and the certified data already trend that way:
`dim U_L / K_L ~ (4^L / (2 sqrt(pi L))) / (4^L / 8) = 4 / sqrt(pi L) -> 0`,
i.e. `W_L / K_L -> 1`.  The certified fractions
`0, .045, .066, .118, .170, .217, .259` (`L = 3..9`) are monotonically
climbing.  No constant-`c` annihilator bound can hold under the law; the
correct all-`L` object is a direct lower bound on `dim U_L` (Section 1),
not an upper bound on `W_L` relative to `K_L`.

**Direction discipline** (lead directive): a *lower* bound on `W_L` would
strengthen the cyclicity refutation and weaken nothing; the derived bound
`dim g_L >= K_L - W_L` needs an *upper* bound on `W_L`.  These are
different statements and are kept separate throughout this note.

---

## 5. Observed annihilator shapes `[COMPUTATION]` + `[UNRESOLVED]`

From the stored exact bases (`results/ladder/w9_basis.json`,
`results/ladder/l8_saturation.json`): every sector-4 kernel vector at
`L <= 9` is an *alternating four-single-rung* vector: for a rung subset
`S = {r_1 < r_2 < r_3 < r_4}` and leg pattern `s in {T,B}^S`, the (orbit)
coefficient is `(-1)^{s_{r_2} + s_{r_4}}`, split into two leg-parity
halves; the count is `C(L,4)` for `5 <= L <= 9`, with the `L = 4`
exception explained by the law's `delta`.  Sector-6 kernels contain the
`d`-decorated copies (one doubly-occupied rung tensored with the same
alternating pattern: 630 of 1134 at `L = 9`) plus genuinely new shapes of
support 6, 8, 11.  The guessed all-`L` family closes under the rung part
of `B` (alternating-sign cancellation) and under hops, and pairwise-cancels
fusions, but pair annihilation adjacent to a decoration generates the
support-6+ shapes; a proved all-`L` invariant family was not completed.
All-`L` statements about `W_L` in either direction remain `[UNRESOLVED]`;
certified values stop at `L = 9`.

---

## 6. [COMPUTATION] Recurrence audit of `W_3..W_9` (holdout discipline)

Exact `Fraction` fits on `(0, 2, 10, 66, 364, 1822, 8586)`
(`ladder_w8_regression.json`):

| candidate | free params | fit pts | holdout pts | verdict |
|---|---|---|---|---|
| `W_{L+1} = c W_L` | 1 | 1 | 5 | no exact fit (degenerate at `W_3 = 0`) |
| order-2 constant | 2 | 2 | 3 | fits `(5, 8)`; REFUTED at `L = 7` |
| order-3 constant | 3 | 3 | 1 | fits `(263/23, -556/23, -1713/23)`; REFUTED at `L = 9` |
| order-1 P-recursive, linear coeffs | 3 | 3 | 2 | REFUTED at `L = 8` |
| symmetric-square closed form (Sec. 4) | 0 | 0 | 7 | UNFALSIFIED on all 7; `[CONJECTURE]` |

Seven terms do not determine a law; every ansatz that consumes points dies
on the remaining ones.  The only survivor consumed none.

---

## 7. Scope, honest negatives, claim ledger

* `W_10` was not computed (lead directive; resource projection).  The
  law's `W_10 = 38950` is a prediction, not a datum.
* The symmetric-square law is `[CONJECTURE]`; both inclusions are open.
  Sector-2 saturation for all `L` is open (certified `L <= 9`).
* Theorem 1.1's family is subexponential; it does not prove `2^L` growth.
  The `2^L` conjecture remains open unconditionally (it follows from the
  law; Section 4).
* The `L = 5` and `L = 10..31` rows of Theorem 2.2 are one-prime modular
  lower-bound certificates, not exact dimensions.
* Nothing here claims an exact `dim g_L` for any `L >= 5`, and nothing is
  claimed about `W_L` for `L >= 10`.
* Resource honesty: the initial certificate sweep and its verifier
  (single processes covering all of `L = 10..31`) peaked at 5.53 GB and
  5.75 GB RSS on the tail rows — above the wave's ~4 GB budget.  The
  per-row peaks are recorded inside `ladder_w8_certificates.json`
  (`L = 30`: 4.6 GB, `L = 31`: 5.5 GB cumulative).  RSS is operational,
  not mathematical: the ranks are unaffected.  The engines now store
  echelon rows as `int32` pairs; the fixed engine re-verifies `L = 20`
  (rank 3161, peak 0.28 GB isolated) and the worst row `L = 31`
  identically, within budget (measurement recorded in the consolidated
  envelope's `rss_remediation` field).

| statement | status |
|---|---|
| `dim g_L >= Q(L)` for all `L >= 2` | `[THEOREM]` (Sec. 1) |
| `Q(L) > 8L^2 - 2L + 1` for all `L >= 32` | `[THEOREM]` (Sec. 1.2) |
| `dim g_L > 2L(4L-1)` at `L = 5, 10..31` | `[COMPUTATION]` one-prime certificates |
| no simultaneous Majorana-quadratic rep of `(A_L, B_L)`, all `L >= 3` | `[THEOREM]` assembled per the seam table (Sec. 2) |
| `L = 2` | `[UNRESOLVED]` by this route; known integrable ring |
| symmetric-square law; `W_10 = 38950` | `[CONJECTURE]` with recorded holdout; falsifiable |
| law implies `K_L - W_L >= 2^L` for all `L` | `[THEOREM]` conditional on the law |
| `W_L <= (1-c)K_L` route | REFUTED as an approach (under the law `W_L/K_L -> 1`) |
| regression `L = 3..8` + `L = 9` artifacts | `[COMPUTATION]`, PASS |
| all-`L` behaviour of `W_L` | `[UNRESOLVED]` |
