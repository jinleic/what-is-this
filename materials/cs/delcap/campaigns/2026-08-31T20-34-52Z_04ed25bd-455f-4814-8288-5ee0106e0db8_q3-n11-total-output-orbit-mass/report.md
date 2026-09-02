# q=3, n=11 total-output-orbit-mass certification report

## Headline and verdict

**MACHINE-VERIFIED:** all four preregistered `q=3`, `n=11`, `d in {1/2,1/5,1/10,1/20}` rows are `CERTIFIED`. Every conservative interval lies strictly inside both published Tavakoli-Nguyen-Bose finite-`n` endpoints: the certified primal lower endpoint beats `LB+` and the certified dual upper endpoint beats `UB`, so every row's verdict is `CERT_LOWER_BEATS_LBplus+CERT_UPPER_BEATS_UB`. Every selected dual is the exact `ba_total_orbit_mass` candidate; exact orbit-uniform was evaluated over the same full input alphabet in every row and never selected.

This campaign closes the `q=3,n=11` cell that the completed `q=4` campaign explicitly excluded and named as its successor. It makes no `q>=4`, `q=3` with `n!=11`, other-`d`, or asymptotic-capacity claim.

Campaign: `campaigns/2026-08-31T20-34-52Z_04ed25bd-455f-4814-8288-5ee0106e0db8_q3-n11-total-output-orbit-mass/`.

Frozen production source SHA-256: `08ecbc493e1c5c674f04ea4a0d0c9c6efd9ac2254fe9d6513db90933b3cdf114`. Manifest SHA-256: `198cb863a104857dfd7b067e2035ebca9b08c6c0400e043c01508d95666a183a`.

## Certified table

Endpoints are conservative outward decimals backed by exact binary rationals and 400-bit Arb enclosures; widths are outward upper bounds in bits per symbol. The registered width target is `1/500`.

| n | d | certified interval | width upper | selected dual |
|---:|:---:|:---|---:|:---|
| 11 | 1/2 | `[0.4136415126778377, 0.4136475665528263]` | `6.053874988461827e-06` | `ba_total_orbit_mass` |
| 11 | 1/5 | `[0.9165946848818924, 0.9166017532024916]` | `7.068320598920295e-06` | `ba_total_orbit_mass` |
| 11 | 1/10 | `[1.2046904010391062, 1.2046925848595933]` | `2.1838204867062363e-06` | `ba_total_orbit_mass` |
| 11 | 1/20 | `[1.3817159992481873, 1.3817177042129642]` | `1.7049647766059455e-06` | `ba_total_orbit_mass` |

Widths run from `1.7049647766059455e-06` at `d=1/20` to `7.068320598920295e-06` at `d=1/5`. The tightest published margin is at `d=1/20`: the archived Arb lower-margin ball there has lower endpoint approximately `0.00129169637501813941901314966788`, and its upper-margin ball has lower endpoint approximately `0.12399667147213431625505939351440`. The exact archived Arb balls and exact outward binary rationals, not these display decimals, are the certificates. `report_only_binary64` fields in the rows are explicitly labeled `REPORT_ONLY_BINARY64_NOT_CERTIFICATE_ENDPOINTS`.

## Exact checks

`MACHINE-VERIFIED`, per row and in aggregate:

- the exact `q=3,n=11` census tuple is `(14884, 22450, 6148309, 3573542, 9721851)` in all four rows, with input/output orbit-size sums `177147 = 3^11` and `265720 = sum_{k=0}^{11} 3^k`, and independent Burnside counts `14884` and `22450`;
- exact input-orbit-total snap denominator `2^30`; exact output-orbit-total snap denominator `2^30` before the frozen full-support bump (post-bump total `1073764274` where the bump fired, and `265720` for exact orbit-uniform);
- both distributions sum to exactly `1`, and every expanded `S_3 x C_2` generator equality passes: `2,125,764` input-distribution checks and `9,565,920` output-distribution checks (induced output law plus both dual candidates);
- the orbit-weighted primal overlaps the direct full-input primal in every row (`full_vs_orbit_primal_overlap = PASS`), across `316,716,419` direct primal conditional entries and zero full-word/representative overlap failures;
- both fixed dual candidates agree between representative and direct full-alphabet maxima in every row: `1,417,176` direct full-input dual evaluations, `634,415,384` positive-`W`/positive-`D` terms, zero overlap failures, zero skipped positive-mass input words;
- the exactly zero-mass primal input words are accounted for rather than ignored: `828` skipped at `d=1/2` against `176,319` positive words, and zero skipped at the other three `d`;
- `ba_total_orbit_mass` needed the frozen full-support bump in three rows (`5448`, `1`, `4` zero orbit masses at `d=1/2,1/10,1/20`) and none at `d=1/5`; after the bump every orbit mass numerator is positive;
- nonnegative certified width, published-sandwich overlap, two-sided containment, and outward endpoint/width containment hold in all four rows.

Post-run immutable-resume validation returned `PASS` with exactly four rows. An independent owner check re-derived all four newline-inclusive canonical JSONL line hashes, confirmed the ledger binds the production source hash, re-verified every outward decimal against its archived exact binary rational, re-checked `width <= 1/500`, and confirmed `TABLE_q3_n11_conservative_intervals.csv` is an exact projection of the canonical rows (`MACHINE-VERIFIED`).

## Static controls and release gates

The pre-production static guard passed with zero BA iterations and zero rows: exact census anchor, independent Burnside counts, `6,480` exhaustive small-instance joint value/reversal equivariance cases, a breaking position plant (`A=1` versus `A=2`), a planted split output law rejected before KL with five failed generator equalities, an admitted exact orbit-uniform control over `1,092` expanded checks, an infinite zero-support dual, a rejected negative-width plant, and empty resume validating exactly zero rows. Details and its independent-review record are in `static_report.md`.

Production ran only under the frozen command: release key `DELCAP_Q3_N11_RELEASED=1`, niceness 10, all five thread pools pinned to one, `PYTHONDONTWRITEBYTECODE=1` enforced before any campaign import, and assertions enabled. The first `d=1/2` row served as the startup anchor and canonical row zero; no separate anchor row was written.

## Timing and resources

`COMPUTATIONAL-EVIDENCE`: one sequential low-priority process. Row stage `1151.072 s` wall / `1150.401 s` CPU; the exact `n=11` structure was built once in `7.455 s` and reused by all four rows; each row's 4,000-iteration float locator took about `85.3-86.0 s`. Maximum row-recorded RSS was `929,677,312` bytes against the 96 GiB cap. No row hit the 90-minute row wall or the six-hour stage wall, and no `NOT_REACHED_RESOURCE` row exists.

## Evidence and provenance

Canonical row evidence is `orbit_rows.jsonl`, bound line-by-line by `orbit_rows.line_checksums.json`. `summary.json` is the terminal runner summary, `TABLE_q3_n11_conservative_intervals.csv` the derived presentation table, and `run.log` / `production_stdout.log` the byte-identical run trace. `pre_statement.md`, `manifest.json`, `checksums_static.sha256`, `static_guard.json`, `static_guard_stdout.log`, `static_report.md`, and `checksums_release.sha256` remain byte-preserved; `manifest_final.json` and `checksums_final.sha256` close the campaign additively without rewriting that history. The row file and its hash ledger are each atomically replaced and `fsync`ed, but the pair is not crash-atomic: a mismatch hard-aborts resume and is never repaired or deleted automatically. No `__pycache__` was produced, because bytecode writing is refused before the first campaign import.

## Rule-7 scope sentence

The completed sweep is exactly `q=3`, `n=11`, `d in {1/2,1/5,1/10,1/20}`, using `G=S_3 x C_2` on symbol values and word reversal, exact total input-orbit mass snap `2^30`, exact total output-orbit mass snap `2^30`, fixed dual candidates `(ba_total_orbit_mass, uniform)`, direct full-alphabet KL comparison over all `3^11` inputs, full finite channel support, and Arb 400-bit outward rounding. **NOT swept:** any `q=4` or `q>=4` row; `q=3` with `n!=11`; other deletion probabilities; `ba_word`; non-`G`-invariant candidates; the falsified `S_n` position/type reduction; the cause of Tavakoli-Nguyen-Bose rounding; Morozov-Duman, Pinto-Ribeiro, and every `n -> infinity` or asymptotic-capacity claim. Each certified row is only a finite-`n` theorem. No published theorem or printed sandwich value is refuted: every certified interval lies inside and strictly improves both published endpoints.

## Next campaign

Named next campaign: **q=3,n=12 total-output-orbit-mass extension**. It requires its own prospective freeze before any BA iteration, because `n=12` is excluded here. Its exact structural targets are already fixed by the same closed forms: `44,530` input `S_3 x C_2` orbits over `3^12 = 531,441` words and `66,980` output orbits over `sum_{k=0}^{12} 3^k = 797,161` words. **INFERENCE cost:** roughly 3x this campaign's per-row cost, so about 1.0-1.5 single-thread CPU hours for its four rows under the same 96 GiB, 90-minute-per-row, and six-hour-stage caps.
