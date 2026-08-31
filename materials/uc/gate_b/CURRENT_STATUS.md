# Gate B: current status

Verified current facts only, newest results first. Every row names the artifact
that proves it. Labels are literal: **PROVED** means proved or machine-checked
here, **CERTIFIED** means an exact two-sided rational enclosure, **NUMERICAL**
means searched or float64 and never a proof step, **OPEN** means undecided.

Chronology lives in [EXPERIMENTS.md](EXPERIMENTS.md); the definitions this page
depends on are frozen in [DEFINITIONS.md](DEFINITIONS.md); the proofs are in
[PROOF.md](PROOF.md).

## The theorem

**PROVED.** With the admissible class and objective of `DEFINITIONS.md`,

```text
c_cl^star = sup { -A_+(F) / eps_vee(F) : F cap/Reimer-admissible, A_+(F) < 0 } = +infinity.
```

No fixed scalar multiple of the ordered-pair closure defect repairs `A_+`.
Three arithmetically disjoint routes certify it: a 256-bit Arb computation with
clamp classification, a standard-library dyadic computation that makes no clamp
decision, and an exact rational computation with no interval library, no clamp
classification and no relaxation of the action set.

## Certified records

Ten bases are registered and exactly certified; one additional symbolic clone
family is certified by exact fixed-order enclosures plus Proposition 25's
reweighting identity. Every number below is the upper endpoint of an exact
rational two-sided enclosure. Normalization (active *and* separating) is a
search and reporting condition in `DEFINITIONS.md`, not a condition in the
displayed supremum, so both conventions are reported.

| record | value | base | separating | certificate |
|---|---|---|---|---|
| **lowest defect with `A_+ < 0`** | **`14/45 = 0.311111`** | `n8m15_clone_reachable`, ratio-41 symbolic clone | no | `gate_b_subtwofifths_clone_rational_v1.json` |
| lowest defect, separating | `1144/1875 = 0.610133` | `n8lo` | yes | `gate_b_n8_rational_v1.json` |
| highest repair ratio | `147/2536 = 0.0579653` | `n8clone_hi` | no | `gate_b_n8_clone_rational_v1.json` |
| highest repair ratio, separating | `531/11360 = 0.0467430` | `n8best` | yes | `gate_b_n8_k4_rational_v1.json` |
| **best asymptotic slope** | **`3/640 = 0.0046875`** | `n8clone_hi` | no | `gate_b_n8_clone_rational_v1.json` |
| **best asymptotic slope, separating** | **`177/40000 = 0.004425`** | `n8best` | yes | `gate_b_n8_k4_rational_v1.json` |

Published values before this round: defect `1336/2025 = 0.659753`, ratio
`0.0362591`, slope `1/250 = 0.004`. All three are improved, and the slope is
improved twice over -- once with a cloned base and once, by `10.6%`, without
leaving the separating class.

The new defect record starts from a normalized fifteen-row core on eight
coordinates, with incidence `45 >= R_15 = 30` and defect `14/45`. Its uniform
objective is positive, but one fixed order is exactly negative. Proposition 25
turns that order into a finite ratio-41 clone family with exact
`A_+ in [-0.000084288238,-0.000084288237]` at unchanged defect. The clone
family is active but not separating, so the normalized frontier remains
`1144/1875`.

## Growth of the numerator: settled, and the constant improved twice

**PROVED, two-sided.** Write `Lambda(n)` for the supremum of `-A_+` over
admissible families on at most `n` coordinates. Then

```text
(3/640) n - 3/80  <=  Lambda(n)  <=  4n/5.
```

The upper bound is unconditional: incidence is at most `n*floor(2m/5)` while
Reimer demands at least `m log2 m / 2`, so **every** admissible family obeys the
size ceiling `log2 m <= 4n/5`, tested exactly as `m^5 <= 2^(4n)` (`PROOF.md`
Lemma 9, Corollary 10). It is close to sharp: the largest admissible sizes are
45 at `n=7` and 445 at `n=11`, against ceilings `2^5.6` and `2^8.8`.

The lower bound improved twice this round. With the cloned-coordinate base the
slope rises from `1/250` to `3/640`, a factor `75/64 = 1.1719` (`PROOF.md`
Proposition 16). Independently, `n8best` from the `S_4 x S_4` class reaches
slope `177/40000 = 0.004425` while remaining **separating**, so the growth
constant improves by `10.6%` even for a reader who insists on normalized
families. So `Lambda(n) = Theta(n)` and only the *denominator* of the repair
ratio remains open.

## Defect-free coordinate extensions: classified, and the old monotonicity claim retracted

**PROVED.** Adding a coordinate with one-rows
\(\mathcal U\subseteq\mathcal F\) never repairs a failed join and preserves the
defect exactly when

```text
1_U(A union B) = 1_U(A) or 1_U(B)
```

for every successful join. Equivalently, `U` is an up-set inside the family and
every successful join landing in `U` has a part in `U`. If `|U| <= floor(2m/5)`,
admissibility is automatic: old degrees stay fixed, the new degree meets the
cap, and incidence rises against a fixed `R_m` (`PROOF.md`, Lemma 19).

The class is larger than clones and OR-columns, but the record core is completely
enumerated: 57 join-consistent row subsets, six usable under cap 6. They are its
five coordinate columns and one genuinely new column. Exact objectives:

| extension of the `n8tiny` core | `A_+ <=` | effect |
|---|---:|---|
| none | `+0.010695694102` | baseline |
| clone coordinate 2 | `+0.003478204249` | improves |
| clone coordinate 0, 1, 3, or 4 | `+0.012500066566` | **worsens** |
| unique non-OR column | `+0.114769081463` | worsens sharply |

So the previous statement that cloning always lowers `Q`, `C_+`, and `A_+` was
**false**. Structural preservation is universal; objective direction is not.
Concentrating three clones on coordinate 2 reaches `-0.000478464900`, while
spreading one over coordinates 0, 1, and 2 leaves `+0.006494611034`, all at
defect `4/9` (`classify_defect_free_extensions.py`).

**PROVED, sharp exact closure.** Clones only reweight original fixed-order
objectives. For multiplicities `r_i`, the collapsed first-appearance order has
the exact Plackett--Luce law

```text
P_r(pi) = product_k r[pi[k]] / sum_{j>=k} r[pi[j]].
```

Geometric multiplicities concentrate on any chosen order, so the infimum and
supremum over all finite clone multisets equal the minimum and maximum original
fixed-order values. In particular, a finite clone multiset is negative **iff**
some original fixed order is negative.

For one coordinate of total multiplicity `r`, the earlier rank law remains

```text
P(first clone has rank j) = C(n+r-j-2,r-1) / C(n+r-1,r).
```

Its convergence is `O(1/r)`, not geometric. For `n8tiny` coordinate 2 the exact
limit is `-0.002338401508...`; total multiplicity 4 first crosses zero. For the
new `14/45` core, no one-coordinate limit is negative, but multicoordinate
ratio `41` is. `audit_clone_limit.py` exactly reweights all 40,320 orders.

## Why the cloned bases are sound

**No lemma uses separation.** Lemma 1 needs only that the two block join events
are independent under independent uniform product rows; Lemmas 2 and 3 need only
block independence and that the four transition probabilities sum to one;
admissibility of the powers is computed from degrees and incidence. Separation
appears once in `PROOF.md`, where it is *concluded* for the powers of a
normalized base, never assumed. Furthermore

```text
Reimer for F^k  <=>  k I m^(k-1) >= ceil(k m^k log2(m)/2)  <=  2I >= m log2 m  <=>  m^m <= 2^(2I),
```

which is exactly the base condition `I >= R_m`: **one base-level Reimer check
certifies every power.** Audited by direct instantiation in
`audit_clone_power.py` (verdict `EVERY_POWER_ADMISSIBLE`; at `k=2` the 4,900-row
product recomputed from its rows has defect `210171/240100 = 1-(173/490)^2`
exactly).

## The local regime: both defect and negativity crossed `2/5`

The original \(n=9,10\) pipeline covered 13 and 15 sizes respectively,
including maxima 145 and 255. A corrected minimum-fixed-order search then used
Proposition 25's exact reachability criterion. Exact and discovery labels remain
separate:

| quantity | value | convention | label | artifact |
|---|---:|---|---|---|
| lowest searched admissible defect | `14/45 = 0.311111...` | normalized *and* displayed | exact witness / heuristic coverage | `local_defect_frontier_n9.json` |
| **lowest certified negative defect** | **`14/45 = 0.311111...`** | displayed, non-separating | exact rational two-sided certificate | `gate_b_subtwofifths_clone_rational_v1.json` |
| lowest certified negative, separating | `1144/1875 = 0.610133...` | normalized | exact rational | `gate_b_n8_rational_v1.json` |

Two normalized fifteen-row cores now explain the mechanism precisely. The
earlier seven-coordinate core
`(0,2,4,6,8,9,16,22,31,32,96,105,112,125,127)` has defect `14/45`,
`A_+ in [0.305904110097,0.305904110098]`, and every fixed order at least
`0.145914214376`; every defect-free extension stays positive. The new
eight-coordinate row-changing core
`(0,1,2,4,5,8,10,43,64,190,192,193,245,254,255)` has the same defect and
positive average `A_+ in [0.230029313333,0.230029313334]`, but fixed order
`(6,1,2,0,3,4,5,7)` lies in
`[-0.003540262562,-0.003540262561]`.

Geometric ratio `41` yields dimension `199623130728`, incidence
`1197738780965`, the same `m=15` and defect, and exact
`A_+ in [-0.000084288238,-0.000084288237]`. The row/order search is heuristic;
the core facts, selected order, all-order reweighting, and final sign are exact.

## Structural obstructions proved this round

* **Lemma 11, union-growth floor.** `eps_vee >= ((8/5) sbar - M) / (n - M)` with
  `M` the largest set size and `sbar` the mean, because the cap forces
  `E|X or Y| >= (8/5) sbar` while a successful union has size at most `M`.
* **Corollary 12, no dominant set means bounded ratio.** If
  `M <= theta (8/5) sbar` with `theta < 1` then `-A_+/eps_vee <= 5n/(4(1-theta))`.
  So any divergence faster than `Theta(n)`, and any approach to the local
  regime, **requires** a set of size at least `(8/5) sbar >= (4/5) log2 m`. The
  `n=6` base sits below its threshold and its floor `7/25` is active; both `n=7`
  bases and every `n=8` witness is above it.
* **Lemma 13, downset capacity.** `eps_vee >= 1 - m^-2 sum_A N(A)^2` with
  `N(A) = #{B in F : B subseteq A}`. Attained with equality: the exhaustive
  audit over all 366 admissible families with `n <= 4` reports tightest slack 0.
* **Proposition 14, the certified bases are maximal.** 25 and 45 are the largest
  admissible sizes at `n=6` and `n=7`, all three certified bases attain them
  with every degree exactly at the cap, and **no set can be added to any of
  them** -- verified over all `2^n` candidates. The one amplification route that
  lowers `eps_vee` and raises `log2 m` at once, adding the missing unions, is
  therefore empty at every certified base.

* **Lemma 19, defect-free extensions.** Join-consistency is necessary and
  sufficient; the defect can never decrease by appending a coordinate. The
  `n8tiny` core has 57 join-consistent subsets, six cap-usable, with exactly one
  improving direction.
* **Proposition 21, clone convex hull.** Arbitrary clone multiplicities only
  reweight original fixed-order objectives. This corrects two earlier claims:
  objective direction is coordinate-dependent, and one-coordinate convergence
  is `O(1/r)`, not geometric.
* **Proposition 22, cap-only floor is false [REPORTED source + in-repo
  derivation].** Chase-Lovett Example 1.4 is active, separating, contains `[n]`,
  has `eps_vee=o(1)`, and every frequency is `psi+o(1)<2/5`; it fails Reimer by
  `h(psi)/2 - psi > 0.0977433528` per coordinate.
* **Proposition 23, zero-defect endpoint.** Reimer's cited average-set-size
  theorem and the dominant set are automatic for a cap union-closed family.
  Deleting zero/duplicate columns normalizes it without changing the endpoint.
* **Proposition 24, no robust Reimer bootstrap.** The Chase-Lovett sequence
  refutes every inequality
  `log2(m)/2-sbar <= n*g(eps_vee)` with `g(t)->0`, even with cap,
  normalization, and a full set.
* **Proposition 25, sharp clone reachability.** The exact Plackett--Luce order
  law makes the clone infimum equal the fixed-order minimum; negative sign is
  finitely reachable iff one fixed order is negative.
* **Proposition 26, exact sub-`2/5` witness.** Ratio-41 multicoordinate cloning
  of the new normalized core gives defect `14/45` and a strictly negative exact
  two-sided `A_+` enclosure under the displayed convention.
* **Exact iid-chain barrier.** Six of ten registered bases satisfy
  `max_order Q_order < log2 m`; the chain fails at every order, with widest
  exact margin `0.0393134578`. All ten bases have maximum frequency exactly
  `2/5`. Cambie's `c* = 0.3823455333667027...` is used only as a
  human-audited upper obstruction; its matching lower verification remains
  **OPEN / COMPUTATIONAL-EVIDENCE** in the broader `math/uc` audit.

## The census objective is a screen, and now we know exactly why

Earlier notes in this file attributed the census/exact discrepancy to float64
imprecision. **That was wrong, and the correction matters.** On the same order
set the float evaluator agrees with the exact rational one to `1e-14`, and the
exact orbit route agrees with all `8! = 40320` orders to the last bit. The fault
is the *order set*:

`block_order_representatives(k)` returns one order per `S_k x S_(8-k)` pattern,
which is a valid representative system for the objective only when the family's
automorphism group **is** that block group. At `k=1` it is, and the `k=1`
census value is exactly right (`-0.027762696467`, ratio `0.034335490331`,
re-verified exactly). At `k=2` it is not: `n8lo` has `|Aut| = 1440`, the same as
`|S_2 x S_6|`, but only **120** of those elements preserve the block partition,
so the 28 block-pattern orders hit just **12 of the 28** true orbits,
double-counting six (multiplicity up to six) and missing sixteen. The resulting
`C_+` is off by `2.3e-2` and `A_+` by `alpha` times that, `8.3e-4`.

Consequences, measured across all four re-ranked classes: of 150 screen-reported
negatives, 147 are certified; four are exactly non-negative; and one family the
screen reported as non-negative is exactly **negative** (`k=3`, defect
`94/125`). The screen errs in both directions, so a threshold of `0` would lose
witnesses.

## Class coverage at eight coordinates: complete

`S_k x S_(8-k)` pairs `k` with `8-k`, so `k in {1,2,3,4}` exhausts every
block-symmetric class at `n = 8`. All four are now enumerated and every
candidate with screen value below `+0.01` has been exactly re-ranked.

| class | cells | raw masks | canonical families | screen-negative | certified negative | false neg | missed |
|---|---|---|---|---|---|---|---|
| `k=1` | 16 | 65,535 | 136 | 2 | 2 | 0 | 0 |
| `k=2` | 21 | 2,097,151 | 2,272 | 21 | 20 | 1 | 0 |
| `k=3` | 24 | 16,777,215 | 13,470 | 54 | 52 | 3 | 1 |
| `k=4` | 25 | 33,554,431 | 11,553 | 73 | 73 | 0 | 0 |
| total | | 52,494,332 | 27,431 | 150 | 147 | 4 | 1 |

The screen is exact at `k=1` and `k=4` and wrong at `k=2` and `k=3`, which is
consistent with the diagnosis above: it is correct exactly when the family's
automorphism group is the block group.

The direct block-symmetric records other than the new symbolic clone family
come from `k=2` or `k=4`. `k=3`'s best, `n8max`, sits at the maximum admissible
size `m = 80` with incidence `256 = 8*32` and held the separating ratio record
until `k=4`.


## Open

* **`c_loc` is undecided.** The certified negative frontier is now `14/45` in
  the displayed cap/Reimer class and `1144/1875` under active separation, but
  neither gives a sequence with defect tending to zero.
* **Further defect reduction has an exact search criterion.** A row-changing
  core can be rescued by clones iff its minimum fixed-order objective is
  negative. The strict search reached `14/45`; no lower-defect seed is currently
  known. Search coverage is heuristic, not an optimality theorem.
* **A Gate B floor must be Reimer-essential.** Proposition 22 kills cap +
  normalization + dominant-set floors, while Proposition 24 kills every
  pair-defect-continuous attempt to recover Reimer. A useful floor must assume
  Reimer as a hard hypothesis or exploit stronger union-image structure.
* **The zero endpoint is still the `2/5` frequency problem.** Proposition 23
  shows Reimer and the dominant set add nothing there. Excluding zero defect in
  all dimensions would beat the rigorously audited `psi = 0.381966...` frontier
  by `0.0180`. Frankl is verified for `n <= 12` [REPORTED], settling only those
  dimensions.
* **The separating frontier remains open.** The core of the new witness is
  normalized, but nontrivial cloning duplicates columns. Reaching below
  `1144/1875` while preserving separation needs a distinct amplification or a
  direct negative family.
* **Proof-assistant formalization.** None of the Gate B certificates is a
  proof-assistant artifact.

## Next frontier

Block-symmetric `n=8` coverage, the requested `n=9,10` low-defect sweep,
full-order `Q` gating, defect-free extension classification, and clone
reachability are closed. Remaining work:

1. **Search below `14/45` by the correct objective.** First produce
   cap/Reimer/dominant-set seeds with smaller defect, then rank by sampled
   minimum fixed-order `A_+` and certify only the selected order. Full-order
   `Q` remains a rejection gate, not the ranking objective.
2. **Exploit Reimer as a hard hypothesis.** Chase-Lovett misses it by
   `0.097743... n`, but that gap cannot be bounded continuously by pair defect.
   A useful theorem must use density/support information unavailable to their
   construction and degrade at the exact zero-defect endpoint.
3. **Preserve separation.** Seek row-changing cores that are already negative
   on the uniform order average, or a non-clone reweighting that keeps columns
   distinct; ordinary cloning cannot move the normalized frontier.
4. **Consider a stronger union-image invariant.** Chase-Lovett is sharp for the
   pair-weighted defect, not for the number or entropy of distinct successful
   union values [REPORTED open direction].
