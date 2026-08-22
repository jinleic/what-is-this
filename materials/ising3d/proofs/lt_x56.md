# Exact low-temperature simple-cubic series through `x^56`, and the re-audit of the seven dormant legacy LT truncation-unobservable survivors

Producer: `experiments/e153_lt_x56.py`.  Artifact: `results/series/lt_x34.json`.
Standalone independent test: `tests/test_lt_x34.py`.

## 0. Conventions

`K=beta J`, `x=e^{-2K}`, `u=x^2=e^{-4K}`, and the reduced low-temperature free energy is normalised as in `problem_specification.md` and `notes/flm_derivation.md`:

    phi(K) = 3*K + sum_{n>=1} c_n x^n,   f_LT(u) = sum_n c_{2n} u^n.

`Xi_A(x)` denotes the fixed-plus-boundary broken-bond polynomial of the `a x b x c` box `A`, i.e. `sum_sigma x^{#unsatisfied bonds}` with the exterior frozen to `+1`.

## 1. Order bound and completeness of the box inventory

[THEOREM] (bounding-box surface bound; proved in `notes/flm_derivation.md`, Eq. (10), and not reproved here).  A connected set `S` of flipped sites spanning `(a,b,c)` has broken-bond surface `|dS| >= 2[(a+b-1)+(a+c-1)+(b+c-1)] = 4*(a+b+c)-6`, and a connected Mayer cluster of droplets inherits the same lower bound through the exterior surface of its union.  Consequently the finite-lattice weight `W_A` of a box `A` satisfies `[x^m] W_A = 0` for all `m < 4*(a+b+c)-6`.

[THEOREM] (inventory completeness at order 56).  Every box whose weight can contribute at or below `x^56` satisfies `4*(a+b+c)-6 <= 56`, i.e. `a+b+c <= 15`; a box with side sum 16 has surface at least `4*16-6 = 58 > 56`.  The computation therefore uses exactly the boxes with `a+b+c <= 15`: 455 ordered boxes (`= C(15,3)`, the number of positive integer triples with that side-sum bound) grouped into 102 permutation classes, one cached polynomial per class.  This inventory is complete for the achieved order and, by the same inequality, not complete for `x^{TARGET_ORDER + 2}`.

[COMPUTATION] The bound is also verified empirically, not merely assumed: every one of the 455 ordered box weights vanishes identically below `4*(a+b+c)-6`, and every canonical box has its first nonzero weight coefficient exactly at `4*(a+b+c)-6` (artifact key `canonical_weight_support`).  A missing or misplaced box would break one of those two statements.

## 2. Exact arithmetic and the CRT uniqueness certificate

All box polynomials, logarithms, bounding-box inversions and survivor residuals are computed in exact Python `int`/`Fraction` arithmetic.  Boxes with at most 62 spins use the `int64` transfer engine; the 23 boxes with more than 62 spins are reconstructed by the repository's deterministic 31-bit-prime CRT engine, truncated at the target order.

[THEOREM] (coefficient bound).  `[x^q] Xi_A` counts spin configurations of the `N`-site box with exactly `q` unsatisfied bonds, hence is a non-negative integer at most `2^N`, and the coefficients sum to exactly `2^N`.

[COMPUTATION] For every CRT box the recorded modulus product strictly exceeds `2 * 2^N`, i.e. strictly more than twice the proved coefficient bound, so the reconstruction of each coefficient is unique.  This is a two-sided certificate: the residues fix the coefficient modulo the product, and the proved bound confines it to an interval shorter than half the modulus.

Because the CRT boxes are truncated at the target order, the global normalisation check `sum_q [x^q] Xi_A = 2^N` is not available for them, and the artifact records that status explicitly rather than silently claiming it.  Uniqueness does not need it: the per-coefficient bound `0 <= [x^q] Xi_A <= 2^N` is a combinatorial fact about each individual coefficient, and only the low-degree prefix enters the logarithm and the bounding-box inversion at this order.  The truncation is what keeps the run inside the memory budget: for the largest box `5x5x5` the unlaunched input-size estimate is 45902462976 bytes truncated against 363193171968 bytes at full degree, the latter exceeding the 51539607552-byte budget.

| box | sites | cross-section | primes | product bits | `2*2^N` | product `> 2*2^N` |
| --- | --- | --- | --- | --- | --- | --- |
| `2x4x8` | 64 | 8 | `2147483647, 2147483629, 2147483587` | 93 | `36893488147419103232` | `True` |
| `2x4x9` | 72 | 8 | `2147483647, 2147483629, 2147483587` | 93 | `9444732965739290427392` | `True` |
| `2x5x7` | 70 | 10 | `2147483647, 2147483629, 2147483587` | 93 | `2361183241434822606848` | `True` |
| `2x5x8` | 80 | 10 | `2147483647, 2147483629, 2147483587` | 93 | `2417851639229258349412352` | `True` |
| `2x6x6` | 72 | 12 | `2147483647, 2147483629, 2147483587` | 93 | `9444732965739290427392` | `True` |
| `2x6x7` | 84 | 12 | `2147483647, 2147483629, 2147483587` | 93 | `38685626227668133590597632` | `True` |
| `3x3x7` | 63 | 9 | `2147483647, 2147483629, 2147483587` | 93 | `18446744073709551616` | `True` |
| `3x3x8` | 72 | 9 | `2147483647, 2147483629, 2147483587` | 93 | `9444732965739290427392` | `True` |
| `3x3x9` | 81 | 9 | `2147483647, 2147483629, 2147483587` | 93 | `4835703278458516698824704` | `True` |
| `3x4x6` | 72 | 12 | `2147483647, 2147483629, 2147483587` | 93 | `9444732965739290427392` | `True` |
| `3x4x7` | 84 | 12 | `2147483647, 2147483629, 2147483587` | 93 | `38685626227668133590597632` | `True` |
| `3x4x8` | 96 | 12 | `2147483647, 2147483629, 2147483587, 2147483579` | 124 | `158456325028528675187087900672` | `True` |
| `3x5x5` | 75 | 15 | `2147483647, 2147483629, 2147483587` | 93 | `75557863725914323419136` | `True` |
| `3x5x6` | 90 | 15 | `2147483647, 2147483629, 2147483587` | 93 | `2475880078570760549798248448` | `True` |
| `3x5x7` | 105 | 15 | `2147483647, 2147483629, 2147483587, 2147483579` | 124 | `81129638414606681695789005144064` | `True` |
| `3x6x6` | 108 | 18 | `2147483647, 2147483629, 2147483587, 2147483579` | 124 | `649037107316853453566312041152512` | `True` |
| `4x4x4` | 64 | 16 | `2147483647, 2147483629, 2147483587` | 93 | `36893488147419103232` | `True` |
| `4x4x5` | 80 | 16 | `2147483647, 2147483629, 2147483587` | 93 | `2417851639229258349412352` | `True` |
| `4x4x6` | 96 | 16 | `2147483647, 2147483629, 2147483587, 2147483579` | 124 | `158456325028528675187087900672` | `True` |
| `4x4x7` | 112 | 16 | `2147483647, 2147483629, 2147483587, 2147483579` | 124 | `10384593717069655257060992658440192` | `True` |
| `4x5x5` | 100 | 20 | `2147483647, 2147483629, 2147483587, 2147483579` | 124 | `2535301200456458802993406410752` | `True` |
| `4x5x6` | 120 | 20 | `2147483647, 2147483629, 2147483587, 2147483579` | 124 | `2658455991569831745807614120560689152` | `True` |
| `5x5x5` | 125 | 25 | `2147483647, 2147483629, 2147483587, 2147483579, 2147483563` | 155 | `85070591730234615865843651857942052864` | `True` |

## 3. The new exact coefficients

[COMPUTATION] The stored predecessor prefix is reproduced exactly before any new coefficient is claimed: all 33 coefficients `x^0,...,x^32` of `results/series/extended2_sc_lt_free_energy.json` agree.  Every odd-order coefficient through `x^56` is exactly zero.  The 12 new coefficients are

    x^34 = 666750
    x^36 = -12520405/6
    x^38 = 6583341
    x^40 = -83409453/4
    x^42 = 464980286/7
    x^44 = -424819905/2
    x^46 = 682202205
    x^48 = -17588511087/8
    x^50 = 35552605353/5
    x^52 = -46131904167/2
    x^54 = 675410878105/9
    x^56 = -979227570369/4

In the `u=x^2` convention this extends the reduced series from 17 to 29 coefficients, `u^0,...,u^28`.

## 4. Post-computation external witness

[EXTERNAL] Guttmann and Enting, *Series studies of the Potts model. I. The simple cubic Ising model*, J. Phys. A 26 (1993) 807--821, Table 2 (local copy `sources/fulltext/num_guttmann_enting1993.pdf`, sha256 `b866f359ca4bc0f75d534f8561cc2a286edd0680a0b37e3f7d8a741830609bff`).  The stored PDF was opened only after every local box polynomial, logarithm, inversion, coefficient check and survivor residual above had been computed; no external number is an input to any local computation, a fit target, or a selection criterion.

Exponentiating the local reduced series, `Lambda_0(u) = exp(f_LT(u))`, reproduces all 25 printed `lambda_n` values of that table, including every new order `x^34`, `x^36`, `x^38`, `x^40`, `x^42`, `x^44`, `x^46`, `x^48`, `x^50`, `x^52`.

This comparison IS a validation gate, and the artifact says so: it is recorded in `checks` as `post_computation_external_table_witness`, so a mismatch fails the producer.  The precise split is therefore: no external value is a computational input and none was used for fitting, selection or tuning; the hashed table is used as a post-computation validation witness, opened after the computation.

| `n` | `x` order | table cell | published `lambda_n` | local `exp` value | match |
| --- | --- | --- | --- | --- | --- |
| 0 | `x^0` | row_aligned | `1` | `1` | `True` |
| 3 | `x^6` | row_aligned | `1` | `1` | `True` |
| 4 | `x^8` | row_aligned | `0` | `0` | `True` |
| 5 | `x^10` | row_aligned | `3` | `3` | `True` |
| 6 | `x^12` | row_aligned | `-3` | `-3` | `True` |
| 7 | `x^14` | row_aligned | `15` | `15` | `True` |
| 8 | `x^16` | row_aligned | `-30` | `-30` | `True` |
| 9 | `x^18` | row_aligned | `101` | `101` | `True` |
| 10 | `x^20` | row_aligned | `-261` | `-261` | `True` |
| 11 | `x^22` | row_aligned | `807` | `807` | `True` |
| 12 | `x^24` | row_aligned | `-2308` | `-2308` | `True` |
| 13 | `x^26` | row_aligned | `7065` | `7065` | `True` |
| 14 | `x^28` | row_aligned | `-21171` | `-21171` | `True` |
| 15 | `x^30` | row_aligned | `65337` | `65337` | `True` |
| 16 | `x^32` | row_aligned | `-200934` | `-200934` | `True` |
| 17 | `x^34` | row_aligned | `627249` | `627249` | `True` |
| 18 | `x^36` | row_aligned | `-1962034` | `-1962034` | `True` |
| 19 | `x^38` | row_aligned | `6192066` | `6192066` | `True` |
| 20 | `x^40` | merged_column_block_token | `-19610346` | `-19610346` | `True` |
| 21 | `x^42` | row_aligned | `62482527` | `62482527` | `True` |
| 22 | `x^44` | merged_column_block_token | `-199807110` | `-199807110` | `True` |
| 23 | `x^46` | merged_column_block_token | `641837193` | `641837193` | `True` |
| 24 | `x^48` | merged_column_block_token | `-2068695927` | `-2068695927` | `True` |
| 25 | `x^50` | merged_column_block_token | `6691611633` | `6691611633` | `True` |
| 26 | `x^52` | merged_column_block_token | `-21710041944` | `-21710041944` | `True` |

The paper's Table 1 independently states that `4x5` cross-sections (`s=14`) give the low-temperature Ising series to `u^26`, which is exactly `x^52` and exactly the local budget `a+b+c <= 15`; its `5x5` row (`s=16`, `2^25` vector elements) is the next step, matching the local wall below.

The six merged-block values were located as verbatim tokens appearing in increasing-n order inside the single merged Table-2 column dump of the stored PDF text layer; they are a weaker witness than the row-aligned cells.

## 5. Re-audit of the eighteen legacy LT survivors

[COMPUTATION] The eighteen rows are the stored LT `CANDIDATE_SURVIVES(!)` truncation-unobservable cells of `results/series/structure_certificates.json` (produced by `experiments/e43_series_structure.py`).  The frozen protocol is reused verbatim: for `U` unknown monomial coefficients the training prefix is the first `U+2` equations, every later equation is a holdout, algebraic rows use equation orders `0..N-1` and first-order differential-algebraic rows only the derivative-safe orders `0..N-2`.  Nothing is refitted or reselected: the training prefix depends only on `U`, and for all eighteen rows the recomputed training kernel basis and training equation count equal the stored ones.

[COMPUTATION] Result: 13 of 18 rows are `CANDIDATE_REFUTED_BY_HOLDOUT` at the newly observable orders and 5 remain truncation-unobservable.

| # | class | degrees | frozen zero columns | verdict | first refuting / next testable | residual on frozen training kernel |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | algebraic | P(x,f), deg_x=0, deg_f=6 | [6] | CANDIDATE_REFUTED_BY_HOLDOUT | refuted at u^18 | `-155793/8, 6833/6, -35/2, 1` |
| 2 | algebraic | P(x,f), deg_x=0, deg_f=7 | [6, 7] | CANDIDATE_REFUTED_BY_HOLDOUT | refuted at u^21 | `-76321/2, 10885/6, -21, 1` |
| 3 | algebraic | P(x,f), deg_x=0, deg_f=8 | [6, 7, 8] | CANDIDATE_REFUTED_BY_HOLDOUT | refuted at u^24 | `337673527/240, -402125/6, 10799/4, -49/2, 1` |
| 4 | algebraic | P(x,f), deg_x=0, deg_f=9 | [6, 7, 8, 9] | CANDIDATE_REFUTED_BY_HOLDOUT | refuted at u^27 | `-1006003633/20, 389007109/144, -109068, 45871/12, -28, 1` |
| 5 | algebraic | P(x,f), deg_x=0, deg_f=10 | [6, 7, 8, 9, 10] | STILL_TRUNCATION_UNOBSERVABLE | next testable u^30 (requires x^60) | `0 at every available holdout` |
| 6 | algebraic | P(x,f), deg_x=0, deg_f=11 | [6, 7, 8, 9, 10, 11] | STILL_TRUNCATION_UNOBSERVABLE | next testable u^30 (requires x^60) | `0 at every available holdout` |
| 7 | algebraic | P(x,f), deg_x=0, deg_f=12 | [6, 7, 8, 9, 10, 11, 12] | STILL_TRUNCATION_UNOBSERVABLE | next testable u^30 (requires x^60) | `0 at every available holdout` |
| 8 | algebraic | P(x,f), deg_x=0, deg_f=13 | [6, 7, 8, 9, 10, 11, 12, 13] | STILL_TRUNCATION_UNOBSERVABLE | next testable u^30 (requires x^60) | `0 at every available holdout` |
| 9 | algebraic | P(x,f), deg_x=1, deg_f=6 | [12, 13] | CANDIDATE_REFUTED_BY_HOLDOUT | refuted at u^19 | `-35/2, 0, 1` |
| 10 | algebraic | P(x,f), deg_x=0, deg_f=14 | [6, 7, 8, 9, 10, 11, 12, 13, 14] | STILL_TRUNCATION_UNOBSERVABLE | next testable u^30 (requires x^60) | `0 at every available holdout` |
| 11 | differential_algebraic_order_1 | Q(x,f,f'), deg_x=0, deg_f=0, deg_f'=8 | [8] | CANDIDATE_REFUTED_BY_HOLDOUT | refuted at u^16 | `426465, 76545, 6561` |
| 12 | differential_algebraic_order_1 | Q(x,f,f'), deg_x=0, deg_f=0, deg_f'=9 | [8, 9] | CANDIDATE_REFUTED_BY_HOLDOUT | refuted at u^18 | `7554627, 1683990, 262440, 19683` |
| 13 | differential_algebraic_order_1 | Q(x,f,f'), deg_x=0, deg_f=0, deg_f'=10 | [8, 9, 10] | CANDIDATE_REFUTED_BY_HOLDOUT | refuted at u^20 | `32684715, 6429780, 885735, 59049` |
| 14 | differential_algebraic_order_1 | Q(x,f,f'), deg_x=0, deg_f=0, deg_f'=11 | [8, 9, 10, 11] | CANDIDATE_REFUTED_BY_HOLDOUT | refuted at u^22 | `609282891, 135655236, 23914845, 2952450, 177147` |
| 15 | differential_algebraic_order_1 | Q(x,f,f'), deg_x=0, deg_f=5, deg_f'=1 | [11] | CANDIDATE_REFUTED_BY_HOLDOUT | refuted at u^17 | `15, -63, 3` |
| 16 | differential_algebraic_order_1 | Q(x,f,f'), deg_x=0, deg_f=0, deg_f'=12 | [8, 9, 10, 11, 12] | CANDIDATE_REFUTED_BY_HOLDOUT | refuted at u^24 | `2711333250, 544845123, 87097275, 9743085, 531441` |
| 17 | differential_algebraic_order_1 | Q(x,f,f'), deg_x=0, deg_f=0, deg_f'=13 | [8, 9, 10, 11, 12, 13] | CANDIDATE_REFUTED_BY_HOLDOUT | refuted at u^26 | `51891290172, 11615705937, 2130783165, 311778720, 31886460, 1594323` |
| 18 | differential_algebraic_order_1 | Q(x,f,f'), deg_x=0, deg_f=6, deg_f'=1 | [6, 12, 13] | CANDIDATE_REFUTED_BY_HOLDOUT | refuted at u^20 | `18, -147/2, 3` |

[COMPUTATION] Why `x^34` is the first order that can test any of these rows: a surviving kernel direction is a literal zero column, i.e. a monomial whose formal first term lies beyond the available prefix.  With the 17 stored `u`-coefficients the derivative-safe equation orders stop at 15, so the `(f')^8` column of the `(0,0,8)` differential-algebraic row (formal valuation 16) had no admissible equation.  The new `x^34 = u^17` coefficient makes equation order 16 derivative-safe and that column acquires the entry `(3*f_3)^8 = 3^8 = 6561`, refuting the row. The standalone test asserts both halves of this: with the old 17-coefficient series the row is `STILL_TRUNCATION_UNOBSERVABLE` with zero column `[8]`, and with the new series it is refuted at holdout order 16 with that exact residual entry.

[UNRESOLVED] The 5 remaining rows are not evidence for a relation.  Each still has a literal zero column and is refutable only once the series reaches the coefficient order listed above: every one of them becomes testable at exactly `x^60`, i.e. precisely the order blocked by the cross-section guard of Section 6, so this re-audit closes at a clean boundary rather than at an arbitrary cut.  A truncation-unobservable survivor is a finite-prefix bookkeeping fact, never a claim that an algebraic or differential-algebraic relation exists.

## 6. Resource walls

[COMPUTATION] OBSERVED for the achieved `x^56` run: wall 2190.220183 s against a 21600 s budget, process peak RSS 47028420608 bytes against a 51539607552 byte budget.  These are measured values (`resource.getrusage` and `time.perf_counter`), not estimates.  The largest box is `5x5x5` (125 spins, 25-spin cross-section) with an unlaunched input-size estimate of 45902462976 bytes at the target truncation.

[COMPUTATION] The next order `x^58` reuses the SAME box inventory (side-sum budget 16) with deeper truncation only.  The first newly walled order is `x^62`: budget 17 requires the canonical `5x6x6` box, whose 30-spin cross-section exceeds the (post-lift) transfer guard 25. No x^56 or higher computation was started, so no larger-order timing or memory figure is claimed.

## 7. Scope

[UNRESOLVED] These are finite-prefix facts about the simple-cubic low-temperature expansion.  Eleven new exact coefficients and the survivor-verdict counts recorded above do not prove non-algebraicity, non-differential-algebraicity or non-D-finiteness of the true low-temperature free energy, do not constrain relations outside the displayed finite ansatz spaces, do not yield a closed form, and do not bound `K_c`.

[EXTERNAL] External-value role, stated without overclaiming in either direction: the critical-coupling benchmark `K_c` is used nowhere here, for fitting, selection, tuning or validation.  The hashed Guttmann--Enting table is likewise not a computational input and was not used for fitting, selection or tuning, but it IS used as a post-computation validation witness, and that comparison is a producer gate (Section 4).  The artifact records exactly this split under `kc_benchmark_used_for_fit_selection_or_validation` and `external_value_role`.
