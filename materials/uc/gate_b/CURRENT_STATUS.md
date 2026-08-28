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

Ten bases are registered and exactly certified. Every number below is the upper
endpoint of an exact rational two-sided enclosure. Normalization (active *and*
separating) is a search and reporting condition in `DEFINITIONS.md`, not a
condition in the displayed supremum, so both conventions are reported and every
row says which it uses.

| record | value | base | separating | certificate |
|---|---|---|---|---|
| **lowest defect with `A_+ < 0`** | **`4/9 = 0.444444`** | `n8tiny` | no | `gate_b_n8_k4_rational_v1.json` |
| lowest defect, separating | `1144/1875 = 0.610133` | `n8lo` | yes | `gate_b_n8_rational_v1.json` |
| highest repair ratio | `147/2536 = 0.0579653` | `n8clone_hi` | no | `gate_b_n8_clone_rational_v1.json` |
| highest repair ratio, separating | `531/11360 = 0.0467430` | `n8best` | yes | `gate_b_n8_k4_rational_v1.json` |
| **best asymptotic slope** | **`3/640 = 0.0046875`** | `n8clone_hi` | no | `gate_b_n8_clone_rational_v1.json` |
| **best asymptotic slope, separating** | **`177/40000 = 0.004425`** | `n8best` | yes | `gate_b_n8_k4_rational_v1.json` |

Published values before this round: defect `1336/2025 = 0.659753`, ratio
`0.0362591`, slope `1/250 = 0.004`. All three are improved, and the slope is
improved twice over -- once with a cloned base and once, by `10.6%`, without
leaving the separating class.

`n8tiny` is worth a second look: it reaches the lowest defect anywhere with only
**fifteen rows**, four of whose eight coordinates coincide. Low defect did not
require a large family; it required a dominant set, which it has (`[8]` itself,
against the Corollary 12 threshold `128/25 = 5.12`).

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

## What a cloned coordinate buys — two different things

Duplicating a coordinate is a bijection commuting with union, so it preserves
`m`, `eps_vee`, every degree (hence the cap), and it only raises the incidence
while `R_m` stays put. **Cloning therefore preserves admissibility
unconditionally and cannot move the defect.** It is free in both quantities the
local question measures. What it moves is `A_+`, through `Q` and `C_+`, and it
moves it *down*.

That gives two distinct uses, and the registered bases use both:

* **Feasibility.** `n8clone_lo` and `n8clone_hi` collapse to *inadmissible*
  seven-coordinate families: incidence 196 against `R_70 = 215`, with the cap
  satisfied. Size 70 simply does not exist at `n=7`, and the eighth coordinate
  is spent entirely on clearing Reimer. The reward is a base whose `-A_+` per
  coordinate, `0.037548/8 = 0.004694`, beats anything available at seven.
* **Negativity.** `n8tiny` is the opposite case. Its four identical coordinates
  collapse to an admissible **five**-coordinate family of the same size 15 and
  the same defect `4/9` — so feasibility was never the issue — but that core has
  `A_+ = +0.010695694`, *positive*. Three clones drive it to `-0.000478465`. The
  clones, not the ground set, are what make it a witness.

**The gain saturates.** Successive clones of one coordinate give deltas
`-7.22e-3, -2.71e-3, -1.25e-3, -6.59e-4`, with ratios `0.375, 0.462, 0.527`: the
improvement is geometric and the budget is bounded, roughly `0.013` for this
family. So cloning flips a family that is positive but *small*, and cannot
rescue one that is far positive — the lowest-defect admissible family known at
`n=7, m=45` has `A_+ = +0.0847`, an order of magnitude outside the budget.

This is the amplification question of the original plan answered in its sharpest
form: the numerator *can* be amplified at zero defect cost, and the amplification
is finite. Audited in `audit_clone_saturation.py`, verdict
`CLONING_IS_DEFECT_FREE_AND_SATURATES`.

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

## The local regime

The only escape left by the theorem is `eps_vee -> 0`. Beyond the frontier rows
above:

| quantity | value | label | artifact |
|---|---|---|---|
| lowest defect of any admissible family, `n=7, m=45` | `692/2025 = 0.341728...` | **NUMERICAL** | `experiments/local_defect_checkpoint.jsonl` |

**The combinatorial floor is not the obstruction.** Admissible families at
`n=7, m=45` reach defect `0.3418`, while the best *negative* one there sits at
`0.6598`. What fails as the defect falls is negativity itself: at fixed `m=45`,
`Q` rises through `log2 45 = 5.4919` (`5.4474` at defect `0.79`, `5.4758` at
`0.66`, `5.5608` at `0.34`), so `A_+` turns positive. Negativity is won only at
the maximum admissible size -- and indeed `n8max`, the best separating ratio,
sits at `m=80` with incidence `256 = 8*32`, the largest any `n=8` family can
have.

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

Every record in the table at the top of this file comes from `k=2` or `k=4`.
`k=3`'s best, `n8max`, sits at the maximum admissible size `m = 80` with
incidence `256 = 8*32` and held the separating ratio record until `k=4`.


## Open

* **`c_loc`, the local ratio, is undecided.** Its upper end is certified at
  defect `4/9` (any admissible) and `1144/1875` (separating); nothing negative
  has been found below `4/9`.
* **Finiteness is Frankl-equivalent.** `c_cl^star` restricted to at most `n`
  coordinates is finite exactly when no admissible family on `n` coordinates has
  `eps_vee = 0` and `A_+ < 0`. A zero-defect admissible family is union closed
  with every frequency at most `2/5`, so proving finiteness for all `n` gives a
  `2/5` frequency bound, about `0.018` beyond the published `0.381966`. Frankl
  is verified for `n <= 12` [REPORTED], which settles those dimensions only.
* **The denominator's growth is open.** The numerator is `Theta(n)`; whether
  `c_cl^star(n)` is `Theta(n)` depends on a positive defect floor, which is the
  Frankl-equivalent statement above.
* **Proof-assistant formalization.** None of the certificates is a
  proof-assistant artifact.

## Next frontier

Block-symmetric coverage at `n = 8` is complete, so the cheap seam is exhausted.

1. `n = 9, 10`, where the maximum admissible sizes are 145 and 255. Three
   lessons from this round change how to search there. Corollary 12 says a
   dominant set is *necessary*, so condition on containing `[n]` or a set of size
   at least `(4/5) log2 m`. `n8tiny` shows the winner may be **small** -- it took
   the defect record on fifteen rows -- so do not restrict to near-maximal sizes.
   And by Proposition 18 the right pipeline is: enumerate low-defect admissible
   families, rank them by `A_+`, and clone only those already within about `0.01`
   of zero. Cloning a cheap positive is the highest-yield move available;
   cloning an expensive one is wasted work, because the budget is finite.
2. Measure the cloning budget per family instead of assuming `n8tiny`'s. The
   saturation ratios `0.375 -> 0.577` are one family's; the limiting sum of
   deltas is what decides which positives are reachable, and it is cheap to
   compute for any candidate.
3. Strengthen the defect lower bound for families that *do* have a dominant set.
   Corollary 12's floor is vacuous exactly when one is present, which is exactly
   the regime every record lives in. Proposition 18 now prices the remaining gap:
   closing the distance from the certified `4/9` to the searched combinatorial
   floor `0.3418` needs `0.0847` of objective, against a cloning budget of order
   `0.013`, so a different mechanism is required and not just more cloning.
4. Fix the screen rather than working around it. The two-stage architecture is
   sound and every label already comes from the exact route, but the screen
   mis-ranks by up to `3.24e-3` whenever the automorphism group is not the block
   group. A per-family automorphism computation costs 40320 permutation tests; a
   cheaper invariant detecting "Aut is larger than the block group" would let the
   screen fall back to all orders only where it must.
5. Whether `4/9` can be pushed below the cap value `2/5`. Nothing yet rules out a
   certified negative family there, and the pair floor `2/m^2` is orders of
   magnitude lower, so the obstruction is not combinatorial.
