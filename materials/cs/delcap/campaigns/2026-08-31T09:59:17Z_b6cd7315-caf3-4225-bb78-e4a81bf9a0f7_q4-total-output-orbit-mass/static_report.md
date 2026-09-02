# q=4 total-output-orbit-mass restart — static freeze report

## Headline

The exact 24-row q=4 restart is preregistered, source-frozen, compile-clean, hash-validated, counterfactually guarded, and runnable. **No q=4 anchor or certificate row was computed:** `orbit_rows.jsonl` does not exist and resume validates exactly 0 rows (`COMPUTATIONAL-EVIDENCE`). The deferred anchor and run remain blocked until Main releases the heavy CPU slot.

Campaign: `campaigns/2026-08-31T09:59:17Z_b6cd7315-caf3-4225-bb78-e4a81bf9a0f7_q4-total-output-orbit-mass/`.

## Frozen construction

The only primal is input-orbit-total: snap total input-orbit masses to `2^30`, divide exactly by orbit size, induce D exactly, and compare orbit-weighted mutual information with direct evaluation over every positive-mass input word. The only dual candidates, in fixed order, are `(ba_total_orbit_mass, uniform)`: snap total output-orbit masses to `2^30`, apply the full-support bump in orbit coordinates, divide exactly by output-orbit size, and compare representative KL with direct KL over all `4^n` inputs. `ba_word` is forbidden. Every exact candidate must pass rational sum, positivity/support, `S_4 x C_2` invariance, full-versus-representative overlap, sandwich containment, nonnegative width, and 400-bit Arb outward-endpoint containment.

Frozen authoritative runner SHA-256: `469e31a7cb59712d0eea6e5d932681fce63ff8eecd78b093d2da37be89ae09d7` (`MACHINE-VERIFIED`). Frozen structural dependency SHA-256: `6d12cc77cefe4ff7cce1efa7770e69cb5c470ca6797f66446756bd4d2b141cdd`; frozen orbit dependency SHA-256: `d17c2b438840bd91e84b775e096d43c60c419ade8ff2ee1ff6e6095694eca9f9`; frozen sandwich dependency SHA-256: `9bee55a5e137ae065a83271586284fa9624cffbfbc7c8d416f7a62ad86a71db9` (all `MACHINE-VERIFIED`).

## Static guard evidence

`static_guard.json` is `PASS` (`MACHINE-VERIFIED`):

- exact q=4,n=5 census tuple `(input G-orbits, output G-orbits, representative entries, orbit-sum slots, visits/iteration)=(31,50,631,332,963)`; input orbit sizes sum to 1,024 and output orbit sizes sum to 1,365;
- 130,560 joint symbol-value/reversal equivariance checks passed;
- the position-permutation plant breaks as required with `A=1` versus `A=2`;
- a real q=4,n=5 output distribution plant has one mixed output orbit, five generator failures, and is rejected before KL; the exact orbit-total uniform control passes;
- zero support on a reachable output returns `+infinity`;
- the negative-width plant is rejected;
- the exact 24-row schema excludes q=3,n=11 and resume validates 0 landed rows.

Guard timing `0.382562 s` CPU is `COMPUTATIONAL-EVIDENCE`; it ran no BA iterations and no capacity row.

## Provenance correction

The row JSONL replacement and line-checksum-ledger replacement are each independently temporary-write + flush + `fsync` + atomic rename + directory `fsync`. The **pair is not crash-atomic**. A crash between replacements preserves a mismatch and causes a hard resume abort; the runner never accepts, truncates, or repairs it. Any failed temporary file is never unlinked: its exact path and existence state are reported to stderr and preserved as evidence. This correction is appended to `pre_statement.md` rather than silently replacing its earlier shorthand.

## Rule-7 scope sentence

The intended sweep is exactly q=4 with `(n,d)` in `{5,6,7,8,9,10} x {1/2,1/5,1/10,1/20}`, using `G=S_4 x C_2` on symbol VALUES and word reversal, exact total input-orbit mass snap `2^30`, exact total output-orbit mass snap `2^30`, fixed dual candidates `(ba_total_orbit_mass, uniform)`, direct full-alphabet KL comparison, full finite channel support, and Arb 400-bit outward rounding. **NOT swept:** q=3,n=11 or any q=3 row; q>=5; q=4,n>=11; other d; `ba_word`; non-G-invariant input or output candidates; the falsified S_n position/type reduction; the cause of Tavakoli-Nguyen-Bose rounding; Morozov-Duman, Pinto-Ribeiro, and every n->infinity/asymptotic capacity claim. Every future landed row is only a finite-n theorem.

## Deferred execution and cost

After Main releases the heavy slot, run the frozen `--anchors` mode first; it computes q=4,n=5,d=1/2 without writing a row. Only a passing anchor permits `--run`. Per-row stop is 90 minutes, whole row-stage stop is 10 hours, RSS stop is 96 GiB, BA is fixed at 4,000 iterations with a resource checkpoint every 25 iterations, and width target is `1/500` bits/symbol. **INFERENCE cost:** 6–10 single-thread CPU hours; q=4,n=10 dominates. Any resource miss is reported without changing the fixed box.
