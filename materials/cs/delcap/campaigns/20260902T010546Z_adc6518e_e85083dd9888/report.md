# q=3, n=12 continuation (box rows 2-3) — certification report

## Headline and verdict

**MACHINE-VERIFIED:** the pre-registered `q=3, n=12` four-row box is now
COMPLETE and every row is `CERTIFIED`, with BOTH published
Tavakoli-Nguyen-Bose finite-`n` endpoints strictly beaten in all four rows
(`CERT_LOWER_BEATS_LBplus+CERT_UPPER_BEATS_UB`).

This campaign certified box rows 2 and 3 (`d=1/10`, `d=1/20`). Box rows 0 and
1 (`d=1/2`, `d=1/5`) are inherited, byte-pinned evidence from the parent
campaign `20260901T121346Z_d1db3155_4354a957b8ba` (closed `CRASHED`; see the
crash provenance below), whose two rows were themselves complete, valid,
ledger-validated certificates.

The result is a **finite-`n` theorem for exactly the four registered cells.**
It makes no `q>=4`, no `q=3` with `n!=12`, no other-`d`, and **no
asymptotic-capacity claim whatsoever.**

## The completed four-row box

Endpoints are conservative outward decimals backed by exact binary rationals
and 400-bit Arb enclosures; widths are outward upper bounds in bits/symbol.
Registered width target `1/500`.

| box row | campaign | d | certified interval (outward) | width upper | selected dual |
|---:|:---|:---:|:---|---:|:---|
| 0 | parent `20260901T121346Z…` | 1/2 | `[0.40486950906723745, 0.4048828292730583]` | `1.3320205820702559e-05` | `ba_total_orbit_mass` |
| 1 | parent `20260901T121346Z…` | 1/5 | `[0.9060415197452321, 0.9060644026068896]` | `2.2882861657318862e-05` | `ba_total_orbit_mass` |
| 2 | this campaign | 1/10 | `[1.1968923784214915, 1.1968988416160002]` | `6.4631945083381036e-06` | `ba_total_orbit_mass` |
| 3 | this campaign | 1/20 | `[1.376928469567818, 1.3769335302038865]` | `5.060636068201896e-06` | `ba_total_orbit_mass` |

Per-row strict improvement over BOTH published endpoints (Arb lower bounds on
the margins, all `> 0`):

| box row | d | `LB+` margin >= | `UB` margin >= |
|---:|:---:|---:|---:|
| 0 | 1/2 | `9.774185e-02` | `3.875984e-01` |
| 1 | 1/5 | `2.095869e-02` | `3.619056e-01` |
| 2 | 1/10 | `5.290779e-03` | `2.295674e-01` |
| 3 | 1/20 | `1.314571e-03` | `1.287808e-01` |

The tightest published margin in the box is at `d=1/20` on the lower side
(`>= 1.314571e-03`). The exact archived Arb balls and exact outward binary
rationals - not these display decimals - are the certificates;
`report_only_binary64` fields are labeled
`REPORT_ONLY_BINARY64_NOT_CERTIFICATE_ENDPOINTS`.

## Pre-registered orbit counts: ASSERTED, not adjusted

`MACHINE-VERIFIED`, in the static guard and again inside every row:

* input `S_3 x C_2` orbits `44,530` over `3^12 = 531,441` words;
* output orbits `66,980` over `sum_{k=0}^{12} 3^k = 797,161` words;
* independent Burnside closed form reproduces both:
  `B(k) = (3^k + 3 + rev(k))/12`, series
  `B(0..12) = 1, 1, 2, 4, 10, 25, 70, 196, 574, 1681, 5002, 14884, 44530`,
  `B(12) = 44530`, `sum_{k=0}^{12} B(k) = 66980` (its `k<=11` prefix is
  `22450`, reproducing the frozen `n=11` anchor);
* constructive `rgs_orbit_data` agrees with Burnside per length;
* orbit-size sums `531441` and `797161` exactly;
* full census tuple `(44530, 66980, 30666848, 18311678, 48978526)`, identical
  in the parent static guard, this campaign's static guard, and all four rows
  (cross-process, cross-campaign determinism).

The `representative_entries`/`orbit_slots` components carried NO
pre-registered value by design (no pre-directory design probe was run for
`n=12`); they were measured once and then re-asserted byte-equal everywhere.

## Exact orbit-mass identities (internal consistency)

`MACHINE-VERIFIED` per row:

* input orbit-mass numerators sum to exactly `2^30 = 1,073,741,824`; the
  per-word expansion `p_x = mass(O)/(2^30 |O|)` sums to exactly `1`;
  `positive_input_words + zero_input_words = 531,441`;
* `ba_total_orbit_mass` numerators sum to exactly `2^30` when no bump fires
  (`d=1/5`) and to exactly `2^30 + 66,980 = 1,073,808,804` when the frozen
  full-support bump fires (`d=1/2, 1/10, 1/20`, with `30904`, `2`, `8` zero
  orbit masses before the bump); after the bump every numerator is positive;
* the exact orbit-uniform dual's numerators sum to exactly `797,161`, the
  total output word count, and its orbit sizes sum to the same;
* every distribution (induced output law and both dual candidates) has exact
  sum `1`.

## Exact checks, in aggregate over the four rows

* `1,553,383,596` direct primal conditional entries; `4,251,528` direct
  full-alphabet dual word evaluations; zero full-word/representative overlap
  failures anywhere;
* `6,377,292` expanded `S_3 x C_2` input-generator checks and `28,697,796`
  expanded output-generator checks (induced law plus both duals);
* zero-mass primal input words accounted for rather than ignored: `32,442`
  skipped at `d=1/2`, zero at the other three `d`;
* nonnegative certified width, published-sandwich overlap, two-sided
  containment, and outward endpoint/width containment in all four rows;
* every selected dual is the exact `ba_total_orbit_mass` candidate; exact
  orbit-uniform was evaluated over the same full input alphabet in every row
  and never selected.

## Controls, in-run, both directions

Run inside THIS campaign after init and before its first row
(`control_inrun.json`, `control_accept_n11.json`):

* **ACCEPT (known-true), direction MATCH:** the frozen `n=11` row
  `0:3:11:1/2` from
  `2026-08-31T20-34-52Z_04ed25bd-…_q3-n11-total-output-orbit-mass` was
  re-executed through THIS `n=12` code path and reproduced **BIT-EXACTLY on
  all 35 compared certificate fields**, including the census tuple
  `(14884, 22450, 6148309, 3573542, 9721851)`, all `14,884` input
  orbit-mass numerators, all `22,450` BA output numerators, its post-bump
  denominator `1,073,764,274`, every archived Arb ball
  (`primal_orbit_ball`, `lbplus`, `ub`, `delta_n`, both margins, both
  per-symbol balls), the full `conservative_interval` block
  (`[0.4136415126778377, 0.4136475665528263]`, width
  `6.053874988461827e-06`) with its exact binary rationals, the verdict
  `CERT_LOWER_BEATS_LBplus+CERT_UPPER_BEATS_UB`, and even the float BA
  locator values (`rate 4.5500566394747235`, `dual 4.5505592605796785`).
  Zero mismatches. The frozen row's own line hash was verified against the
  frozen per-line ledger before use.
* **REJECT plants, direction REJECT (all fired at their pre-registered
  assertion points):** `R1a` input orbit-mass total `2^30+1` rejected by the
  exact normalization assertion; `R1c` negative input numerator rejected by
  the nonnegativity assertion; `R1b` negative output numerator rejected by
  the dual positivity gate before any KL evaluation; `R2` dual-below-primal
  endpoints rejected by the negative-width assertion; `R3` planted split
  (non-`G`-invariant) word-level output law rejected before representative
  compression with `5` failing generator equalities, while the exact
  orbit-total-uniform control over the same structure was ADMITTED with
  `1,092` expanded checks.

The parent campaign ran the identical battery with identical outcomes
(`control_accept_n11_inherited_parent_08836f43.json`,
`static_guard_inherited_parent_20260901T123035Z.json`, preserved here as
inherited evidence under stamped names).

## Crash provenance of the parent, and the instrument defect it exposed

The parent campaign committed box rows 0 and 1, then died on box row 2 with
`BrokenPipeError(32)` raised inside `log()`'s `print` at
`2026-09-01T14:05:14Z`, archiving `FAILURE.json` with `completed_rows: 2`.
Root cause is the **launcher, not the instrument**: the launch piped the
runner's stdout through `tail`, deviating from the frozen command (which does
not pipe), and the read end was closed mid-run. The runner's own frozen
semantics then behaved exactly as designed - hard-abort, archive
`FAILURE.json`, forbid rerun inside that campaign - so the parent was frozen
and closed `CRASHED` with its two valid rows intact and byte-pinned
(`orbit_rows.jsonl` SHA-256
`3c5cf81a91512eeacdf94fc0131bf2564068ff53cd875734b6195e1514cf8bb8`).

**Disclosed instrument defect (inherited, not introduced):** `log()` writes to
stdout AND `run.log`, so a broken stdout is fatal mid-row even though all
evidence is already durable on disk. The defect is in the logging layer only;
no certificate arithmetic is implicated, and rows 0-1 predate any change.
This campaign excludes the failure mode by contract (no pipe on stdout;
supervised process with stdout to a file) rather than by patching frozen
mathematics. Two further inherited disclosures carried forward from the
pre-statement: the frozen largest-remainder `snap` docstring's own cycling
caveat, and the fact that `ru_maxrss`/wall figures are
`COMPUTATIONAL-EVIDENCE` while all CPU claims use in-process
`time.process_time()` (repository rule 17e).

Two owner-side verification mistakes were made and corrected during
close-out, and are recorded because they were mine, not the artifact's: an
initial re-check compared the archived outward width against the difference of
the two *display doubles* (legitimately ~2 ULP larger than the true Arb
endpoint difference the width bounds), and a second compared the post-bump BA
denominator against `2^30 + z_zeros` instead of `2^30 + n_out` (the frozen
bump adds `+1` to every orbit). Both were errors in the check, not in the
evidence; the corrected 180-assertion re-derivation passes.

## Runner identity

This campaign's runner
(`gate_c_q3_n12_continuation.py.asrun`, SHA-256
`5bafe94f1625a31de72284e69d8ca458d4032dc72576c1cbd56adfe7d6097a4f`) is the
parent runner
(`08836f43adbffcbf35da90feb140ca030e387605043802243dd3dac4e6cc6718`) with
exactly 11 mechanical host-side substitutions (campaign id, schema tag,
runner/table filenames, release key, prereg filename, and the two-row `d`
box with its `len(grid())` assertions). A line diff modulo the campaign id
confirms no arithmetic, snapping, candidate, comparison, or interval logic
differs. The three frozen dependencies are hash-pinned and byte-identical to
the `n=11`/`q=4` freeze
(`d17c2b43…`, `9bee55a5…`, `6d12cc77…`), as is the frozen `n=11` row evidence
used by the ACCEPT control.

## Timing and resources

`COMPUTATIONAL-EVIDENCE`, one sequential low-priority (`nice -n 15`,
single-thread-pinned) process per stage. This campaign: rows stage
`4299.939 s` wall / `4033.048 s` CPU for two rows; static guard `48.5 s` CPU;
control battery `6m55s`. Across the whole four-row box the row CPU total is
`7664.5 s = 2.13 CPU-hours`, inside the pre-registered 6-CPU-hour stage
budget, with maximum row-recorded RSS `4,398,972,928` bytes against the
96 GiB cap. No row hit the 90-minute row wall; no `NOT_REACHED_RESOURCE` row
exists in either campaign. The pre-registered budget did NOT bind.

## Evidence and provenance

Canonical rows: `orbit_rows.jsonl`, bound line-by-line by
`orbit_rows.line_checksums.json` (both re-derived independently at close-out:
newline-inclusive line hashes, canonical-JSON byte equality, source-hash
binding). `summary.json` is the terminal runner summary,
`TABLE_q3_n12_rows2_3_conservative_intervals.csv` the derived presentation
table, `run.log` the run trace, `static_guard.json` the pre-production guard,
`control_inrun.json` / `control_accept_n11.json` the in-run controls,
`instrument_manifest.json` the frozen pre-result instrument manifest,
`pre_statement_q3_n12_continuation.md` the byte-identical prereg copy, and
`PREREG_PROVENANCE.json` its source commit (`29eda4e`) and SHA-256
(`6e19324c78b084bc2fc497f8e3ace66f9010ae693212a09fa14f42c327a33de8`).
`sha256s.txt` pins the frozen run dir; `status.json` records the terminal
verdict.

## Rule-7 scope sentence

The completed sweep is exactly `q=3`, `n=12`, `d in {1/2,1/5,1/10,1/20}` (rows
2-3 here, rows 0-1 inherited from the parent), using `G=S_3 x C_2` on symbol
values and word reversal, exact total input-orbit mass snap `2^30`, exact
total output-orbit mass snap `2^30` with the frozen full-support bump, fixed
dual candidates `(ba_total_orbit_mass, uniform)`, direct full-alphabet KL
comparison over all `3^12 = 531,441` inputs, full finite channel support, and
Arb 400-bit outward rounding. **NOT swept:** any `q=4` or `q>=4` row; `q=3`
with `n!=12` (beyond the single pre-registered `n=11` ACCEPT control); other
deletion probabilities; `ba_word`; non-`G`-invariant candidates; the falsified
position/type reduction; the cause of Tavakoli-Nguyen-Bose rounding;
Morozov-Duman, Pinto-Ribeiro, and every `n -> infinity` or
asymptotic-capacity claim. Each certified row is only a finite-`n` theorem.
No published theorem is contradicted: all four certified intervals lie
strictly inside the published sandwich.

## Next gate

Named next campaign: **`q=3, n=13` total-output-orbit-mass extension**,
requiring its own prospective freeze. Its structural targets follow the same
closed forms: input orbits `B(13) = (3^13 + 3 + 2*3^7)/12 = 133,225` over
`3^13 = 1,594,323` words; output orbits
`sum_{k=0}^{13} B(k) = 200,205` over `sum_{k=0}^{13} 3^k = 2,391,484` words
(machine-recomputed at close-out; an earlier draft of this line carried
`132,860`/`199,840`, which was wrong and is corrected here).
**INFERENCE cost:** the measured `n=11 -> n=12` step was `1150 s -> 7664 s`
row CPU (6.7x) at 3x the word count and 5x the representative entries, so
`n=13` is estimated at `6-12` single-thread CPU-hours for four rows and
about `12-14 GiB` peak RSS - it needs a larger pre-registered budget than the
6-CPU-hour cap used here, and should be pre-registered as such rather than
resumed under changed caps.
