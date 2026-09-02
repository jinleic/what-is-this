# Report — new-decomposition off-diagonal census

## Headline

**No new lower-bound-55 data triple enters from `mws59` or `stapleton60`; `paper55` remains the only verified 55-addition circuit in the expanded swept set, while `sun56` remains LB 55 / UB 56.** No lower bound at most 54 occurs. `[MACHINE-VERIFIED]`

Frozen campaign: `cs/mm3/campaigns/2026-08-31T083323Z_55447c11-5887-4fdd-aa2f-30e31f2332e8_b3046e6a6c50/`.

## Exact verdicts

| decomposition | survivor data triples | valid parameterized orientations | certified LB range | upper median | minimizing triples | nonmonomial-involving minimum | LB at most 54 |
|---|---:|---:|---:|---:|---:|---:|---:|
| `mws59` | 5,796 | 1,922,973,696 | 56–74 | 67 | 2 | 56 | 0 |
| `stapleton60` | 4,800 | 1,592,524,800 | 58–76 | 67 | 9 | 58 | 0 |
| aggregate | 10,596 | 3,515,498,496 | 56–76 | 67 | 2 | 56 | 0 |

Every count and bound in the table is an exact finite census or complete finite subset-DFS result. `[MACHINE-VERIFIED]` The 215.63-second total is operational metadata only. `[COMPUTATIONAL-EVIDENCE]`

These are certified lower bounds for the registered data triples and parameterized orientations. They are neither synthesized upper-bound circuits nor counts of distinct canonical tensor arrays.

## Known witnesses and remaining gaps

The startup controls reproduce all five published decompositions over the full 729-entry Brent tensor: `paper55`, Perminov-58, `sun56`, `mws59`, and Stapleton-60 all pass 729/729 under exact integer and `fmpz` evaluation. Their published totals replay as 55, 58, 56, 59, and 60. `[MACHINE-VERIFIED]`

Both `paper55` and `sun56` have a certified lower-bound-55 all-monomial data triple. Only `paper55` has an exact verified 55-addition upper-bound circuit. `sun56` therefore remains LB 55 / UB 56. The new `mws59` and `stapleton60` rows acquire no upper-bound witness from this campaign; their all-monomial lower bounds are 56 and 58, respectively. Published circuit provenance is `[CITED-DEPENDENCY]`; the local exact replays are `[MACHINE-VERIFIED]`.

## Instrument adequacy

The run exactly enumerated all 19,683 ternary `3 x 3` matrices, including 6,960 unimodular matrices, and proved the unique `6,960 = 145 x 48` right-monomial factorization. It independently checked dense and factored survivor counts, exact pair-table bounds, a 21,025-triple anchored subcube, 200 random concordance cases per decomposition, 24 NumPy/Python entries, and 25 right-monomial fibers. Honest, wrong-`W`, transpose-as-inverse, direct-corruption, and determinant-2 plants all gave their registered verdicts. All twelve activity rows have three sides of rank 9 with 23 nonzero rows and nine nonzero columns. The `+14` transposition direction is rank-audited. Every admitted pair uses a complete subset-DFS floor census; no float, SAT, or UNSAT claim enters the adjudication. `[MACHINE-VERIFIED]`

An independent owner validator regenerated survivor tables, side lower bounds, totals, histograms, medians, minimizers, crosschecks, aggregate counts, and all escalation predicates. It passed 476 assertions. Its first attempt expected `fiber_total_lb` in all 24 fresh checks, although the frozen schema provides that field only in the first 12; it aborted before a verdict and left an empty write-once output. The corrected retry preserves that failed attempt and passes. This audit-script incident does not alter any frozen campaign output.

## Rule-5 incident

The predecessor campaign `2026-08-31T082412Z_a2d87381-157e-45bb-ad3a-519ddf899bc1_4dc928ac6dbe/` remains frozen. It aborted in 0.26 seconds before its first anchor by iterating the intentional null Perminov split. It produced no target claim. The elapsed time is `[COMPUTATIONAL-EVIDENCE]` only.

## Rule-7 scope

This campaign swept exactly the two decompositions `mws59` and `stapleton60` under every independent `(P,Q,R)` in `T^3`, where `T` is the exactly enumerated 6,960-element set of ternary `3 x 3` integer determinant-`+/-1` matrices, followed by all three post-sandwich sigma powers, for 2,022,921,216,000 raw scheme-instances represented without sampling by the proved `145^3` data-node factorization per decomposition; it admitted exactly 10,596 ternary data-triples (3,515,498,496 orientations after the `48^3` and three-sigma multiplicities), and decided every admitted data-triple with exact integer `d`-counting, complete subset-DFS floor census, and the rank-audited `+14` transposition gap under the standard action `U'=P^-1 U Q^-T`, `V'=Q^T V R^-T`, `W'=P^T W R`; the alphabet was exactly `{-1,0,1}`, parameters were the two fixed decompositions, all 6,960 choices for each of `P,Q,R`, and `sigma=0,1,2`, radii were none, and no instance in that named domain was sampled or left undecided; it did not sweep `paper55`/`sun56` beyond startup controls, any decomposition outside the five named public ones (including Laderman, other Smirnov schemes, Schwartz–Vaknin, unpublished, or F2-only constructions), non-ternary factor alphabets, `GL(3,Q)` or `GL(3,Z)` sandwiches beyond `T`, actions outside the proved family, anti-cyclic swaps that compute `B*A`, or upper-bound circuit synthesis, so it establishes no universal no-54 claim outside the named domain.

## Prompt premises disproved

No published-paper claim is disproved. The bounded possibility that this registered ternary standard-action family introduces a lower-bound-55-or-less data triple from `mws59` or `stapleton60` is closed negatively: their exact minima are 56 and 58. This says nothing about other actions, alphabets, decompositions, or upper-bound synthesis.

## Named next campaign

`laderman23_source_lock`: first lock primary-source bytes, factor conventions, and a verified Laderman-23 decomposition, then apply the same census only if the action and transposition accounting remain source-compatible. Expected cost is 2–6 hours of provenance/convention work plus 15–60 minutes of census at the observed state scale, longer if its DFS state spaces grow. `[INFERENCE]`
