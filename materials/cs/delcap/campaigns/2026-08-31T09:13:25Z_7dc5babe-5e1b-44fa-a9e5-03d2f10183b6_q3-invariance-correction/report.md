# q=3 output-orbit-mass correction report

## Headline and verdict

**MACHINE-VERIFIED:** all 20 prior q=3 finite-n rows were re-certified with exact total-output-orbit-mass dual references. Exact audit found that **16/20** frozen accepted `ba_word` vectors split an output orbit; those 16 representative-compressed upper endpoints are retracted and replaced below. The other 4 legacy vectors were invariant. Every corrected row still strictly improves both the exact Tavakoli-Nguyen-Bose `LB+` and `UB`; no external theorem or printed sandwich value is refuted.

**FAILURE TO CERTIFY / NOT RUN for q=4:** the parent q=4 campaign stopped at its startup gate. Exactly zero q=4 certificate rows ran. Owner integration of this correction is required before a q=4 restart.

Correction campaign: `campaigns/2026-08-31T09:13:25Z_7dc5babe-5e1b-44fa-a9e5-03d2f10183b6_q3-invariance-correction/`.

Frozen sources:

- correction runner SHA-256 `ed84c1a3982f243e49e8f8f35367ee9e8788bb3c8a47fed16b883bd9da9efa21`;
- localized v3 dependency SHA-256 `6d12cc77cefe4ff7cce1efa7770e69cb5c470ca6797f66446756bd4d2b141cdd`;
- conservative renderer SHA-256 `bc8901f1038d09650fae9dea1caaaab25861f5102ae595efa3e17601782e1b15`.

## Retraction and source fix

The frozen source comment claiming that equal invariant float values necessarily receive equal per-word largest-remainder counts is false: stable tie-breaking can split equal values at the rounding boundary. The first in-grid witness was `(3,6,1/5)`, flat output orbit 51 of size 12, split between exact counts 578028 and 578029. The exact 20-row audit found 16 affected rows.

Primal and dual roles are kept separate in the correction:

- **Primal:** snap total input-orbit masses to `2^30`, set `p_x=mass(O)/(2^30|O|)`, induce the output marginal exactly, and require direct-full-alphabet mutual information to overlap the orbit-weighted value.
- **Dual:** snap total output-orbit masses to `2^30`, apply the full-support bump in orbit coordinates, set `D_y=mass(T)/(S|T|)`, and require representative KL maxima to overlap direct maxima over every one of the `3^n` inputs. The only fallback is exact uniform `D`; candidate choice is the smaller direct-full-alphabet Arb upper endpoint in fixed order.
- The legacy per-word vector is audit-only and is rejected before representative compression whenever split.

Every chosen dual was `ba_total_orbit_mass`. Across the 20 replacements, the runner performed **1,058,508 direct full-input KL evaluations**, **5,821,704 exact group-generator checks**, and **237,744,816 positive-W/positive-D term checks**; zero positive-mass input words were skipped. These counts and all certificate statements are `MACHINE-VERIFIED`. Runtime `841.039 s` wall / `838.785 s` CPU is `COMPUTATIONAL-EVIDENCE` only.

## Corrected 20-row table

The `new certified interval` endpoints below are conservative outward decimals from the full 400-bit per-symbol Arb balls. Each decimal's exact binary rational was asserted outside the corresponding Arb endpoint. The binary64 `cert_lo_per_symbol`, `cert_hi_per_symbol`, and width fields in `corrected_q3_rows.jsonl` are **REPORT-ONLY**, not certificate endpoints. Every row's published-margin verdict is `CERT_LOWER_BEATS_LBplus + CERT_UPPER_BEATS_UB` (`MACHINE-VERIFIED`).

| n | d | frozen report interval | new certified interval (outward) | width upper (outward) | legacy D split? |
|---:|:---:|:---|:---|---:|:---:|
| 6 | 1/2 | `[0.48612690873525904, 0.4861310395537632]` | `[0.48612690873525893, 0.48612726836569464]` | `3.5963043559742424e-07` | no |
| 6 | 1/5 | `[0.9992805994038764, 0.9992807615729179]` | `[0.9992805994038761, 0.9992806185950381]` | `1.919116167554878e-08` | yes |
| 6 | 1/10 | `[1.2624515385268185, 1.2624516450349905]` | `[1.262451538526818, 1.2624515509732321]` | `1.2446413657816448e-08` | yes |
| 6 | 1/20 | `[1.4158072326903548, 1.4158073313105772]` | `[1.4158072326903546, 1.4158072449168306]` | `1.2226475645090283e-08` | yes |
| 7 | 1/2 | `[0.4657259682042402, 0.46575328765404533]` | `[0.46572596820424017, 0.46573204857053757]` | `6.080366297333852e-06` | yes |
| 7 | 1/5 | `[0.9769505177740003, 0.9769508993374995]` | `[0.9769505177740003, 0.976950563790745]` | `4.6016744441483404e-08` | yes |
| 7 | 1/10 | `[1.2473597256274753, 1.2473599906772777]` | `[1.2473597256274749, 1.247359775851469]` | `5.022399355970927e-08` | yes |
| 7 | 1/20 | `[1.4070914100656966, 1.40709164898683]` | `[1.4070914100656966, 1.4070914501217517]` | `4.005605471757835e-08` | no |
| 8 | 1/2 | `[0.4491602653164888, 0.4492922601763715]` | `[0.4491602653164887, 0.44918246544877494]` | `2.220013228610702e-05` | yes |
| 8 | 1/5 | `[0.9582453246473194, 0.9582467564944788]` | `[0.9582453246473193, 0.9582455337027981]` | `2.09055478620426e-07` | yes |
| 8 | 1/10 | `[1.2344520895003903, 1.2344530557652682]` | `[1.23445208950039, 1.2344521805124733]` | `9.101208270304075e-08` | yes |
| 8 | 1/20 | `[1.3995415062373293, 1.3995424502513716]` | `[1.399541506237329, 1.399541571855359]` | `6.56180294580417e-08` | yes |
| 9 | 1/2 | `[0.4353814110192391, 0.4353857919473774]` | `[0.43538141101923905, 0.4353817875980781]` | `3.765788389715567e-07` | yes |
| 9 | 1/5 | `[0.9423194816428749, 0.942324555218804]` | `[0.9423194816428748, 0.9423200290509085]` | `5.474080334643068e-07` | no |
| 9 | 1/10 | `[1.2232499579738394, 1.2232523623416272]` | `[1.2232499579738394, 1.2232502176397095]` | `2.5966586972603587e-07` | yes |
| 9 | 1/20 | `[1.3929068203510149, 1.39290943510016]` | `[1.3929068203510147, 1.3929070242522732]` | `2.03901257961709e-07` | yes |
| 10 | 1/2 | `[0.42369975643641544, 0.4237122878711516]` | `[0.4236997564364154, 0.42370155403328696]` | `1.7975968714839528e-06` | no |
| 10 | 1/5 | `[0.9285797540357409, 0.9285961038799806]` | `[0.9285797540357409, 0.9285822403778092]` | `2.4863420681301732e-06` | yes |
| 10 | 1/10 | `[1.213412964521387, 1.2134231653523089]` | `[1.2134129645213867, 1.2134136998846528]` | `7.353632655671914e-07` | yes |
| 10 | 1/20 | `[1.3870088535516936, 1.3870161689880416]` | `[1.3870088535516931, 1.387009425126998]` | `5.715753042954663e-07` | yes |

The width range is `1.2226475645090283e-08` to `2.220013228610702e-05` bits/symbol (`MACHINE-VERIFIED` outward upper renders). Exact Arb balls and exact outward binary rationals are in `conservative_intervals.jsonl`; its line checksum ledger verifies all 20 rows.

## Rule-7 scope sentence

The correction sweep is exactly q=3 with `(n,d)` in `{6,7,8,9,10} x {1/2,1/5,1/10,1/20}`, using `G=S_3 x C_2` on symbol VALUES and reversal, exact total input-orbit mass snap `2^30`, exact total output-orbit mass snap `2^30`, full finite channel support for every dual candidate, direct full input alphabets of size `3^n`, and Arb 400-bit outward rounding. **NOT swept:** q=4 or q>=4; q=3 n>=11; other d; per-word tie-breaking as an admitted bound; non-G-invariant correction candidates; TNB rounding cause; Morozov-Duman, Pinto-Ribeiro, and every n->infinity/asymptotic capacity claim. Every replacement is only a finite-n theorem.

## Disproved premises and prospective/retrospective status

Disproved: (1) equal orbit-invariant float inputs automatically remain equal under per-word largest-remainder snapping; (2) all 20 frozen accepted `ba_word` vectors were invariant. Exact audit found 16 split and 4 invariant.

The hazard was discovered retrospectively during the parent q=4 startup gate. The separate correction box, input/output role separation, total-orbit-mass formulas, candidate order, direct-full checks, precision, endpoint rule, and resource stops were prospective and frozen before the correction runner existed.

## Next campaign and cost

Named next campaign: **q4-total-output-orbit-mass restart**. After owner integration, replay q=4 `n=5..10` at `d={1/2,1/5,1/10,1/20}`, then separately adjudicate q=3 n=11, using this corrected output-orbit-total-mass construction and conservative Arb rendering. **INFERENCE cost:** 6–10 single-thread CPU hours under the existing 96 GiB, 90-minute-per-row, and 10-hour-stage caps; q=4 n=10 has 11,685,028 sparse nonzero visits per BA iteration and dominates.
