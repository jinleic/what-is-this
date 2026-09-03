# q=3, n=13 total-output-orbit-mass certification report

## Headline and verdict

**MACHINE-VERIFIED:** all four pre-registered `q=3, n=13`,
`d in {1/2,1/5,1/10,1/20}` rows are `CERTIFIED`, and **every row's
conservative interval strictly beats BOTH published Tavakoli-Nguyen-Bose
finite-`n` endpoints** (`CERT_LOWER_BEATS_LBplus+CERT_UPPER_BEATS_UB`).
Widths `1.52e-05` to `9.37e-05` bits/symbol, all inside the registered
`1/500` target. Zero `VALID_BUT_WIDTH_MISS` rows, zero resource rows: the
pre-registered budget did NOT bind.

This is a **finite-`n` theorem for exactly the four registered cells.** It
makes no `q>=4`, no `q=3` with `n!=13`, no other-`d`, and **no
asymptotic-capacity claim whatsoever.**

## Certified table

Endpoints are conservative outward decimals backed by exact binary rationals
and 400-bit Arb enclosures; widths are outward upper bounds in bits/symbol;
registered width target `1/500`.

| n | d | certified interval (outward) | width upper | selected dual |
|---:|:---:|:---|---:|:---|
| 13 | 1/2 | `[0.39713666568100003, 0.397161636764393]` | `2.4971083392840876e-05` | `ba_total_orbit_mass` |
| 13 | 1/5 | `[0.8966732662961678, 0.8967669676702715]` | `9.370137410351712e-05` | `ba_total_orbit_mass` |
| 13 | 1/10 | `[1.1898717660767768, 1.1898901620041096]` | `1.8395927332359186e-05` | `ba_total_orbit_mass` |
| 13 | 1/20 | `[1.3725688358567392, 1.3725840117520043]` | `1.5175895264739525e-05` | `ba_total_orbit_mass` |

Per-row strict improvement over BOTH published endpoints (Arb lower bounds on
the margins, all `> 0`):

| d | `LB+` margin `>=` | `UB` margin `>=` |
|:---:|---:|---:|
| 1/2 | `9.859301e-02` | `3.953196e-01` |
| 1/5 | `2.129289e-02` | `3.712030e-01` |
| 1/10 | `5.372100e-03` | `2.365761e-01` |
| 1/20 | `1.333958e-03` | `1.331304e-01` |

The tightest published margin in the box is the `d=1/20` lower side
(`>= 1.334e-03`). The exact archived Arb balls and exact outward binary
rationals - not these display decimals - are the certificates.
`report_only_binary64` fields are labeled
`REPORT_ONLY_BINARY64_NOT_CERTIFICATE_ENDPOINTS`. Every selected dual is the
exact `ba_total_orbit_mass` candidate; exact orbit-uniform was evaluated over
the full input alphabet in every row and never selected.

## Pre-registered orbit counts: ASSERTED, never adjusted

`MACHINE-VERIFIED` in the static guard and again in every row:

* input orbits `133,225 = B(13) = (3^13 + 3 + 2*3^7)/12` over
  `3^13 = 1,594,323` words;
* output orbits `200,205 = sum_{k=0}^{13} B(k)` over
  `sum_{k=0}^{13} 3^k = 2,391,484` words;
* the machine-checked Burnside series
  `B(0..13) = 1, 1, 2, 4, 10, 25, 70, 196, 574, 1681, 5002, 14884, 44530,
  133225`; its `k<=11` prefix `22,450` and `k<=12` prefix `66,980` reproduce
  the frozen `n=11` and `n=12` anchors exactly;
* constructive `rgs_orbit_data(3,k)` orbit counts equal `B(k)` for every
  `k = 0..13`;
* measured dense census tuple
  `(133225, 200205, 152963378, 93991184, 246954562)` — byte-identical across
  the static guard, the two production structure builds (control battery and
  row stage), and all four rows;
* per-rep closure identity: `sum_t |O_j| F_j(t) == |O_j| * 2^13` for every
  input orbit (frozen constructor assertion);
* the pre-registered projection (`~1.53e8` / `~9.4e7` / `~2.5e8`) matched the
  measured tuple to within 0.1% / 0.2%.

## Exact orbit-mass identities (per row)

* input numerators sum to exactly `2^30 = 1,073,741,824`; per-word expansion
  sums to exactly `1`; `positive + zero input words = 1,594,323`;
* `ba_total_orbit_mass` numerators sum to exactly `2^30` with no bump
  (`d=1/5`) and to exactly `2^30 + 200,205 = 1,073,942,029` where the frozen
  full-support bump fires (`131,380` / `4` / `18` zero orbit masses at
  `d=1/2, 1/10, 1/20`); post-bump every numerator is positive;
* orbit-uniform numerators sum to exactly `2,391,484`, the total output word
  count, equal to the orbit-size sum;
* induced output law and both dual candidates each sum to exactly `1`.

## Exact work, in aggregate over the four rows

* `7,446,791,969` direct primal conditional entries;
* `19,131,876` expanded `S_3 x C_2` input-generator checks and `86,093,424`
  expanded output-generator checks (induced law plus both duals);
* zero full-word/representative overlap failures anywhere;
* zero-mass primal input words accounted for rather than ignored: `289,980`
  skipped at `d=1/2`, zero at the other three `d`;
* nonnegative width, published-sandwich overlap, two-sided containment, and
  outward endpoint/width containment in all four rows.

## Controls, in-run, both directions

Archived in `control_inrun.json` and `control_accept_n12.json`, run after
init and before the first `n=13` row:

* **ACCEPT (known-true), direction MATCH:** the frozen `n=12` row
  `0:3:12:1/10` (from the `FROZEN-CERTIFIED` campaign
  `20260902T010546Z_adc6518e_e85083dd9888`, whose newline-inclusive line hash
  was verified against its frozen per-line ledger first) was re-executed with
  THIS runner on the same frozen dependency bytes and reproduced **BIT-EXACTLY
  on all 41 compared certificate fields**, zero mismatches: census tuple
  `(44530, 66980, 30666848, 18311678, 48978526)` with orbit-size sums
  `531,441`/`797,161`, all `44,530` input numerators, the BA dual's `66,980`
  output numerators and post-bump denominator `1,073,808,804` with bump zeros
  `2`, every archived Arb ball, the full `conservative_interval` block
  (`[1.1968923784214915, 1.1968988416160002]`, width
  `6.4631945083381036e-06`) with exact binary rationals, the verdict, and the
  float BA locator value `14.362708541093044`.
* **REJECT plants, direction REJECT, all fired at their named assertions:**
  `R1a` input total `2^30+1` (normalization assertion); `R1b` negative output
  numerator (dual positivity gate, before any KL); `R1c` negative input
  numerator (nonnegativity assertion); `R2` dual-below-primal endpoints
  (negative-width assertion); `R3` planted split output law on the `q=3,n=5`
  control structure (rejected before representative compression, failing
  generator equalities, while the exact orbit-total-uniform control over the
  same structure was ADMITTED).

## Static guard and memory-feasibility probe

`static_guard.json`, zero BA iterations, zero rows: dependency and evidence
hashes verified before import; Burnside series asserted; full `n=13`
structure built and its census pinned; the 32 GiB RSS cap was NOT tripped —
measured build peak `16,827,580,416` bytes (49.0% of cap), build `383.66 s`
wall — so the pre-registered budget-bound freeze path was not taken and the
row stage proceeded; negative-width plant rejected; frozen `n=12` ACCEPT
target validated; empty resume validated exactly zero rows.

## Instrument identity

Runner `gate_c_q3_n13_total_orbit.py.asrun` SHA-256
`569585824a13abf8d0e17c2d67522e5008985a7bc066cfb4f83c2be1c39bf0c4`, derived
from the `n=12` continuation runner
(`5bafe94f1625a31de72284e69d8ca458d4032dc72576c1cbd56adfe7d6097a4f`). A
per-function AST diff proves **22 certificate/ledger functions byte-identical**
(locator, kl_evaluation, certificate_row, conservative_interval, outward
endpoint helpers, input-distribution assertion, commit/validate_resume,
table, resource rows, plants, atomic writers, resource checks); the only
changed functions are `rows_stage` (row-count assertion `2 -> 4`) and
`assert_preregistered_census` (n=13 targets, n=12 control targets); the only
new functions are `load_frozen_n12_row`, `accept_control_n12`,
`pinned_dense_tuple`. The three frozen dependencies are hash-pinned and
byte-identical to the `n=11`/`n=12`/`q=4` freeze
(`d17c2b43…`, `9bee55a5…`, `6d12cc77…`).

## Timing and resources

`COMPUTATIONAL-EVIDENCE`, one sequential low-priority process (`nice -n 15`,
single-thread-pinned), stdout to a file per the pre-registered launcher rule —
no pipe anywhere; the `n=12` BrokenPipe failure mode never recurred. Stage:
`40,605.0 s` wall / `38,668.4 s` CPU = **10.74 CPU-hours of the registered
24**, wall 11.28 h of 26; per-row CPU `8,118-10,943 s` (2.25-3.04 h) against
the 5 h row wall; max row-recorded RSS `16,669,851,648` bytes against the
32 GiB cap (48.5%); the `n=13` structure was built twice (control battery and
row stage) at `419 s` each and re-asserted byte-equal. No resource stop; no
`NOT_REACHED_RESOURCE` row; the pre-registered budget did not bind. Including
controls and static guard, the campaign consumed `~11.9 CPU-hours`.

## Crash/failure status

None. Unlike the `n=12` parent, this campaign completed its stage in one
process with `exit 0`, `FAILURE.json` absent, terminal artifacts written by
the runner itself.

## Rule-7 scope sentence

The completed sweep is exactly `q=3`, `n=13`, `d in {1/2,1/5,1/10,1/20}`,
using `G=S_3 x C_2` on symbol values and word reversal, exact total
input-orbit mass snap `2^30`, exact total output-orbit mass snap `2^30` with
the frozen full-support bump, fixed dual candidates
`(ba_total_orbit_mass, uniform)`, direct full-alphabet KL comparison over all
`3^13 = 1,594,323` inputs, full finite channel support, and Arb 400-bit
outward rounding. **NOT swept:** any `q=4` or `q>=4` row; `q=3` with `n!=13`
(beyond the single pre-registered `n=12` ACCEPT control); other deletion
probabilities; `ba_word`; non-`G`-invariant candidates; the falsified
position/type reduction; the cause of Tavakoli-Nguyen-Bose rounding;
Morozov-Duman, Pinto-Ribeiro, and every `n -> infinity` or
asymptotic-capacity claim. Each certified row is only a finite-`n` theorem.
No published theorem is contradicted: all four certified intervals lie
strictly inside the published sandwich.

## Next gate: q=3, n=14 — evaluated honestly, opened only if it fits

Measured scaling of row CPU: `n=11` `287.6 s`/row -> `n=12` `1,916.1 s`/row
(`6.66x`) -> `n=13` `9,573.8 s`/row (`5.00x`). The `n=12 -> n=13` per-row
full-alphabet KL pass measured `~32-41 min` (`~640-820` words/s); `n=14`
doubles the word count per pass again (`3^14 = 4,782,969` inputs) and grows
the representative/orbit-slot operator by the same per-representative density
ratio as before (`~1.67x` per rep at `3x` reps):

* projected `n=14` census: input orbits `B(14) = (3^14 + 3 + 4*3^7)/12 =
  399,310` over `3^14 = 4,782,969` words; output orbits
  `sum_{k=0}^{14} B(k) = 599,515` over `sum_{k=0}^{14} 3^k = 7,174,453`
  words (machine-recomputed at close-out; an earlier draft carried
  `394,024` / `594,229`, which was wrong and is corrected here);
  `representative_entries ~ 5.1e8`, `orbit_slots ~ 3.1e8`;
* projected per-row CPU `~5.0e4 s (~13.9 h)`, four rows `~55.5 CPU-hours`;
  stage wall `~2.5 days` single-threaded;
* measured `n=13` peak RSS `16.67 GB` at `2.47e8` total sparse slots implies
  `n=14` peak `~35-40 GB` of live list/array structures — **above this
  workstation's pre-registered 32 GiB envelope**, with real OOM/thrash risk
  against sibling campaigns;
* the row wall cap alone would need `~20 h` and the CPU budget `~56 h`:
  more than twice this campaign's, in one process, on a shared machine.

**Decision (per the pre-registered rule "opened only if a fresh prereg budget
fits, otherwise recorded and stopped"):** `q=3, n=14` is **NOT opened** by
this agent. The projection above is recorded as the named next campaign with
its estimated cost. Opening it requires a fresh preregistration that at
minimum (a) raises the RSS envelope to `~64 GiB` with the memory-feasibility
probe retained and the resource-stop rule explicit, (b) budgets `~64 CPU-h`
and multi-day stage wall, and (c) justifies the single-thread wall
multiplication — or changes the instrument (e.g. a pre-registered chunked or
array-backed structure build) so the memory profile fits; any such
instrument change is a NEW freeze with its own controls, not a modification
of this one. The `n=13` box stands complete on its own.

Canonical evidence: `orbit_rows.jsonl` + `orbit_rows.line_checksums.json`
(row-keyed SHA-256 ledger), `control_accept_n12.json`, `control_inrun.json`,
`static_guard.json`, `instrument_manifest.json`,
`pre_statement_q3_n13.md` (byte-identical to committed
`cs/delcap/pre_statement_q3_n13.md`, sha256
`748391d512bd929da5de908322f38e47975c9022fec46ba1b326d86a6e3b90df`, commit
`7bea130`), `PREREG_PROVENANCE.json`, `summary.json`,
`TABLE_q3_n13_conservative_intervals.csv`, `run.log`, `sha256s.txt`,
`status.json`.
