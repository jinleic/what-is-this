# Pre-statement — q=3, n=13 total-output-orbit-mass extension

Written and committed BEFORE any `n=13` compute: no BA iteration, no `n=13`
structure build, no `n=13` census measurement, and no certificate row exists
when this file is written. Successor to the completed `q=3,n=12` four-row box
(campaigns `20260901T121346Z_d1db3155_4354a957b8ba`, verdict `CRASHED`, box
rows 0-1 certified and byte-pinned; and
`20260902T010546Z_adc6518e_e85083dd9888`, verdict `FROZEN-CERTIFIED`, box rows
2-3), named as the next gate in that campaign's report.

## Question and adjudication (one finite question)

Under the same corrected exact total-input-orbit-mass primal and
total-output-orbit-mass dual protocol, can each of the four registered
`(q,n)=(3,13)` deletion-channel rows be enclosed by a valid two-sided
conservative interval of outward width at most `1/500` bits/symbol, and does
each such interval strictly improve BOTH published Tavakoli-Nguyen-Bose
(arXiv:2607.19559) finite-`n` sandwich endpoints, as all four `n=12` rows did?

This is a finite-`n` target. **No asymptotic-capacity claim is registered,
sought, or permitted, at any `n`, in any direction.** A bounded instance is a
bounded instance: rows the budget does not reach are reported as unreached,
never as evidence.

## Frozen row box and order

Exactly four rows, in this order (the inherited deletion-probability grid):

1. `(q,n,d)=(3,13,1/2)`;
2. `(3,13,1/5)`;
3. `(3,13,1/10)`;
4. `(3,13,1/20)`.

No other row may be added to this campaign. Rows are attempted in this order
and committed one at a time, so a budget stop freezes a prefix of the box.

## Exact orbit counts to assert (pre-registered values)

Orbits are `S_3 x C_2`-orbits acting on symbol VALUES (transposition `(0 1)`
and 3-cycle `(0 1 2)`) and on whole-word reversal; position permutation is
forbidden (the falsified reduction). The constructive path is the frozen
restricted-growth-string canonicalization plus reversal merge
(`rgs_orbit_data`/`GOrbits`), with no `q^n` enumeration.

Independent closed form (Burnside), for `k >= 1`:

```
B(k) = (3^k + 3 + rev(k)) / 12,
rev(k) = 4 * 3^(k/2)     if k even,
rev(k) = 2 * 3^((k+1)/2) if k odd,
B(0) = 1.
```

Series, exact integer arithmetic, hand- and machine-checked at
pre-registration: `B(0..13) = 1, 1, 2, 4, 10, 25, 70, 196, 574, 1681, 5002,
14884, 44530, 133225`. Its `k<=11` prefix sum is `22,450` and its `k<=12`
prefix sum is `66,980`, reproducing the frozen `n=11` and `n=12` output-orbit
anchors exactly.

**Pre-registered assertion values for `n=13`:**

* input orbits `133,225 = B(13) = (3^13 + 3 + 2*3^7)/12` over
  `3^13 = 1,594,323` input words; assertions `len(reps) == 133225`,
  `len(reps) == burnside(13)`, `sum(orbit sizes) == 1594323`;
* output orbits `200,205 = sum_{k=0}^{13} B(k)` over
  `sum_{k=0}^{13} 3^k = 2,391,484` output words; assertions
  `n_out == 200205`, `n_out == sum(burnside(k) for k in 0..13)`,
  `sum(flat orbit sizes) == 2391484`;
* per-length agreement: the constructive `rgs_orbit_data(3,k)` orbit count
  equals `B(k)` for every `k = 0..13`;
* per-input-rep exact closure identity from the frozen constructor: for every
  input orbit `O_j`, `sum_t |O_j| * F_j(t) == |O_j| * 2^13` (each length-13
  word has exactly `2^13` subsequences of length `<= 13`, counted with
  multiplicity).

**A mismatch of any reproduced orbit count against these values is a
STOP-AND-REPORT hard abort (terminal verdict `CRASHED`, mismatch archived),
never a silent adjustment and never a post-hoc reinterpretation.**

## Census tuple convention

The census tuple is, in this fixed order,

```
(input G-orbits, output G-orbits, representative_entries, orbit_slots,
 total_sparse_slots = representative_entries + orbit_slots)
```

where `representative_entries` counts `(input representative, distinct output
subsequence)` incidences over all input orbit representatives, and
`orbit_slots` counts the nonzero aggregated orbit-to-orbit slots of the
marginal operator. The first two components are pre-registered above. The
dense components (`representative_entries`, `orbit_slots`,
`total_sparse_slots`) carry **no** pre-registered numeric value: no
pre-directory design probe is run for `n=13`. They are measured exactly once
by the static guard, pinned into the frozen `static_guard.json`, and then
re-asserted byte-equal by **every** production row — a cross-stage
determinism check of the structure construction. For reference, the frozen
`n=12` tuple is `(44530, 66980, 30666848, 18311678, 48978526)` and is asserted
during this campaign's `n=12` ACCEPT control.

## Frozen mathematical protocol (inherited verbatim)

* Symmetry group `G = S_3 x C_2` on symbol values and whole-word reversal.
* Float locator: 4,000 Blahut-Arimoto iterations in orbit coordinates; float
  values locate candidates only and are `COMPUTATIONAL-EVIDENCE`, never
  certificate endpoints.
* Primal: snap total input-orbit masses to denominator `2^30` (frozen exact
  largest-remainder snap); for input orbit `O`,
  `p_x = mass(O)/(2^30 * |O|)`; the full output law is induced exactly.
* Dual candidates in fixed order `(ba_total_orbit_mass, uniform)`; snap total
  BA output-orbit masses to `2^30`; apply the frozen full-support bump (+1 to
  every orbit) if any orbit mass is zero; for output orbit `T`,
  `D_y = mass(T)/(S * |T|)`. `ba_word` is forbidden.
* Candidate selection: minimum finite direct-full-alphabet Arb upper endpoint,
  fixed tie order above.
* Arithmetic: exact `fmpz`/`fmpq` integers and rationals for all distributions
  and channel weights; 400-bit Arb (`ctx.prec = 400`) outward-rounded interval
  arithmetic for every logarithm, mutual information, KL maximum, published
  bound, margin and endpoint. Toolchain pinned to the frozen instrument
  environment: python 3.14.3, python-flint 0.9.0, numpy 2.5.2, scipy 1.18.1
  (`cs/.venv`, single-thread pinned).
* **Conservative two-sided interval construction:** certified interval
  `= [primal_per_symbol.lower(), dual_per_symbol.upper()]` taken from the
  400-bit Arb balls; each presented endpoint is moved one binary64 ULP
  OUTWARD via `math.nextafter`, converted back to the EXACT binary rational,
  and machine-asserted to contain the Arb endpoint (outward lower `<=` Arb
  lower endpoint; outward upper `>=` Arb upper endpoint; outward width `>=`
  the exact Arb endpoint difference). Binary64 display values are labeled
  `REPORT_ONLY_BINARY64_NOT_CERTIFICATE_ENDPOINTS`. Note for auditors: the
  two outward endpoints are each moved a further ULP, so their difference may
  exceed the certified width by ~2 ULP; the width bounds the TRUE Arb
  endpoint difference, which is the quantity claimed.
* **Width target `1/500` bits/symbol.** A row is `CERTIFIED` only if its
  exact primal and admitted exact dual form a nonnegative-width interval whose
  conservative outward width is at most `1/500`. A mathematically valid wider
  interval is `VALID_BUT_WIDTH_MISS`, reported honestly as such. A strict
  lower improvement requires the certified primal lower endpoint to exceed the
  published `LB+` upper endpoint; a strict upper improvement requires the
  published `UB` lower endpoint to exceed the certified dual upper endpoint. A
  row with neither is `NO_STRICT_IMPROVEMENT` and is never tuned away.

Every admitted row must pass, before reporting a bound: exact input/output
normalization; positivity/full support where required; expanded
`S_3 x C_2` generator equality over all words; direct full-input primal
agreement with the orbit-weighted primal; direct full-input KL maximum
agreement with the representative maximum for both dual candidates;
nonnegative certified width; overlap with the published finite-`n` sandwich;
and conservative outward endpoint/width containment.

## Exact orbit-mass identities (asserted per row)

* input orbit-mass numerators sum to exactly `2^30 = 1,073,741,824`, and the
  per-word expansion sums to exactly `1`;
* `positive_input_words + zero_input_words == 1,594,323`;
* `ba_total_orbit_mass` numerators sum to exactly `2^30` when no bump fires
  and to exactly `2^30 + 200,205 = 1,073,942,029` when the frozen
  full-support bump fires (the bump adds `+1` to every one of the `200,205`
  output orbits); after any bump every numerator is positive;
* the exact orbit-uniform dual's numerators sum to exactly `2,391,484`, the
  total output word count, equal to its orbit-size sum;
* induced output law and both dual candidates each have exact sum `1`.

## In-run controls (both directions, after init, before the first n=13 row)

A control failure is a hard abort with the stop archived: no row, summary or
table may be written afterwards.

### ACCEPT (known-true, required direction MATCH)

The frozen `n=12` row `0:3:12:1/10` from the `FROZEN-CERTIFIED` campaign
`20260902T010546Z_adc6518e_e85083dd9888` is re-executed with THIS campaign's
runner on the same frozen dependency bytes and pinned toolchain, and must
reproduce it **BIT-EXACTLY on all certificate fields**: the census tuple
`(44530, 66980, 30666848, 18311678, 48978526)` with its `n=12` orbit-size
sums `531441` and `797161`; all `44,530` input orbit-mass
numerators and the denominator `1073741824`; the BA dual's `66,980` output
numerators, its `total_mass_denominator` and its
`zero_orbit_masses_before_full_support_bump`; every archived Arb ball
(`primal_orbit_ball`, `lbplus_arb_ball`, `ub_arb_ball`, `delta_n_arb_ball`,
both margin balls, both per-symbol balls); the entire `conservative_interval`
block including the exact binary rationals, i.e.
`[1.1968923784214915, 1.1968988416160002]` with width upper
`6.4631945083381036e-06`; `status = CERTIFIED`;
`verdict = CERT_LOWER_BEATS_LBplus+CERT_UPPER_BEATS_UB`;
`dual_from = ba_total_orbit_mass`; the published-containment flags; and the
float BA locator values. The frozen row's newline-inclusive line SHA-256 is
verified against its frozen per-line ledger before use. Excluded from the
comparison: wall/CPU/RSS/timing fields only. Archived as
`control_accept_n12.json`; committed to no row file.

### REJECT plants (each MUST be rejected, at the named assertion)

1. **R1a** input orbit-mass vector with total `2^30 + 1` — rejected by the
   exact normalization assertion before any distribution is built;
2. **R1b** output orbit-mass vector containing `-1` — rejected by the dual
   positivity gate before any KL evaluation;
3. **R1c** input orbit-mass vector containing `-1` — rejected by the
   nonnegativity assertion;
4. **R2** dual upper endpoint below primal lower endpoint — rejected by the
   `conservative_interval` negative-width assertion, producing no interval;
5. **R3** planted split (non-`G`-invariant) word-level output law on the
   `q=3,n=5` control structure — rejected BEFORE representative compression
   with failing generator equalities, while the exact orbit-total-uniform
   control over the same structure is ADMITTED.

Archived as `control_inrun.json`.

## Instrument identity and frozen dependencies

The runner is the `n=12` continuation runner (SHA-256
`5bafe94f1625a31de72284e69d8ca458d4032dc72576c1cbd56adfe7d6097a4f`) with
mechanical host-side changes only: campaign id and paths, schema tag,
runner/table filenames, release key, prereg filename, the `n=13` row box and
its count assertions, the `n=13` pre-registered orbit counts, the `n=12`
ACCEPT-control target, the census-anchor pinning path, and the resource caps
below. **No arithmetic, snapping, candidate, comparison, interval, resume or
ledger logic may change**, and a line diff modulo the campaign id is recorded
in the closing report to demonstrate exactly that.

Hash-pinned before import (unchanged from the `n=11`/`n=12`/`q=4` freezes):

* `gate_c_orbit_dependency.py.asrun`
  `d17c2b438840bd91e84b775e096d43c60c419ade8ff2ee1ff6e6095694eca9f9`;
* `gate_c_dhalf_dependency.py.asrun`
  `9bee55a5e137ae065a83271586284fa9624cffbfbc7c8d416f7a62ad86a71db9`;
* `gate_c_q4_structure_dependency.py.asrun`
  `6d12cc77cefe4ff7cce1efa7770e69cb5c470ca6797f66446756bd4d2b141cdd`.

Hash-pinned frozen ACCEPT-control evidence (read-only):
`20260902T010546Z_adc6518e_e85083dd9888/orbit_rows.jsonl` and its
`orbit_rows.line_checksums.json`.

## Launcher rule (the BrokenPipe lesson, pre-registered)

The `n=12` parent campaign died mid-row with `BrokenPipeError(32)` raised
inside the runner's `log()` `print` because the LAUNCHER piped the runner's
stdout through `tail`, deviating from the frozen command, and the read end was
closed. The instrument's mathematics was never implicated, and `log()` is left
FROZEN AND UNMODIFIED here. Instead the launch contract excludes the failure
mode:

* **the runner's stdout and stderr MUST be a regular file or a supervised
  process capture; piping them into any process whose reader can vanish
  (`| tail`, `| head`, `| grep`, a harness-managed pipe) is PROHIBITED;**
* stdin MUST be `/dev/null`;
* the process MUST outlive the launching shell (supervised process manager),
  because a stage is expected to run for hours;
* every stage is one sequential process; no parallel row workers.

Exact frozen production command, cwd `cs/` (stdout/stderr to a file, no pipe):

```sh
env PYTHONDONTWRITEBYTECODE=1 \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  DELCAP_Q3_N13_RELEASED=1 \
  nice -n 15 ./.venv/bin/python \
  delcap/campaigns/<THIS_ID>/gate_c_q3_n13_total_orbit.py.asrun \
  --run < /dev/null >> delcap/campaigns/<THIS_ID>/production_stdout.log 2>&1
```

The static-guard and `--controls` stages use the same rule with their own
capture files. Production additionally requires
`DELCAP_Q3_N13_RELEASED=1`, niceness at least 15, all five thread pools
pinned to `1`, `PYTHONDONTWRITEBYTECODE=1` enforced before any campaign
import, and Python assertions enabled (`-O` refused).

## Immutable resume with per-row byte pinning

* Canonical rows are newline-terminated canonical JSON (`sort_keys`, compact
  separators) in `orbit_rows.jsonl`, appended one row at a time in box order.
* Each row's newline-inclusive SHA-256 is recorded in
  `orbit_rows.line_checksums.json` under the key `ordinal:q:n:d`, together
  with the production source hash; the ledger is rewritten atomically with the
  row file (`O_EXCL` temp, `write`, `flush`, `fsync`, `os.replace`, directory
  `fsync`).
* Every resume re-reads and re-validates ALL committed rows: canonical-JSON
  byte equality, ordinal/box agreement, campaign id, schema, source hash, key
  bijection and per-line hash. Any mismatch, malformed prefix, wrong ordering,
  wrong source hash, or surviving `.*.tmp.*` staging file HARD-ABORTS and is
  never repaired or deleted automatically.
* The pair is deliberately not crash-atomic; recovery is exhaustive
  recomputation of the next missing ordinal only. Committed rows are never
  recomputed and never overwritten.
* Terminal `FAILURE.json`, `summary.json` or table evidence forbids rerun
  inside this campaign; a crashed stage is frozen and closed with its verdict,
  and any continuation is a separate pre-registered campaign (the pattern
  already used for `n=12`).

## Pre-registered hard budget and the resource-stop rule

Sized from measured rates, not guesses. Measured: `n=11` row CPU
`1150.401 s` for four rows, `n=12` row CPU `7664.515 s` for four rows
(`6.66x`), representative entries `6,148,309 -> 30,666,848` (`4.99x`, i.e.
`413.1 -> 688.7` per representative, ratio `1.667`), orbit slots
`3,573,542 -> 18,311,678` (`5.12x`), peak RSS `929,677,312 ->
4,398,972,928` bytes (`4.73x`), structure build `7.455 s -> 57.38 s`.
Projected for `n=13` (`INFERENCE`): `representative_entries ~ 1.53e8`,
`orbit_slots ~ 9.4e7`, `total ~ 2.5e8`; per-row CPU `~2.65 h`; four rows
`~10.6 CPU-hours`; peak RSS `~20.6 GiB`; structure build `~5 min`.

Registered caps, enforced in-process:

* **stage CPU budget 24 single-thread CPU-hours** (`86,400 s`, measured by
  in-process `time.process_time()` — the only CPU accounting this workstation
  does not distort, repository rule 17e);
* **per-row wall 5 hours** (`18,000 s`) and **stage wall 26 hours**
  (`93,600 s`), SIGALRM-enforced;
* **RSS 32 GiB** (`34,359,738,368` bytes) via in-process
  `resource.getrusage` (Darwin byte semantics asserted);
* niceness `>= 15`, one sequential process, all thread pools `1`.

**Resource-stop rule, explicit.** When any cap trips: the in-flight cell is
committed as a `NOT_REACHED_RESOURCE` row carrying the exact tripped cap and
the measured cost; every later cell is committed as
`NOT_REACHED_RESOURCE_UPSTREAM`; the runner summary status becomes
`FAILURE_TO_CERTIFY_RESOURCE`; the campaign is frozen and closed
`FROZEN-INCONCLUSIVE` with the exact frontier — which cells are certified,
which are unreached, and the measured cost of the unreached remainder. **The
campaign is never silently resumed under changed caps**; a larger budget is a
new pre-registered campaign.

**The static guard doubles as the memory-feasibility probe.** It builds the
full `n=13` structure with zero BA iterations and zero rows and records the
peak RSS. If that build trips the 32 GiB cap, the campaign freezes
`FROZEN-INCONCLUSIVE` (budget-bound, first-class negative) with the measured
peak and the exact remaining cost, and **no row is attempted** — the 10+
CPU-hour row stage is not burned against a known-infeasible memory profile.

## Verdict rule (pre-registered, exactly one terminal verdict)

* **`FROZEN-CERTIFIED`** iff all four rows are `CERTIFIED`, the pre-registered
  orbit counts hold, and every ACCEPT and REJECT control passes in its
  required direction. The report then states, per row, whether BOTH published
  endpoints are strictly beaten.
* **`FROZEN-INCONCLUSIVE`** iff the box is only partly certified — budget-bound
  (any `NOT_REACHED_RESOURCE*` row, or a guard-stage memory stop) or a
  mathematically valid `VALID_BUT_WIDTH_MISS` row — with the exact frontier
  and remaining cost stated. Certified rows inside a partial box remain
  certified finite-`n` theorems and are reported as such.
* **`FROZEN-NEGATIVE`** iff some row's certified interval strictly EXCLUDES
  the published sandwich, i.e. the published finite-`n` sandwich is falsified
  at that cell; report the direction and never generalize beyond that cell.
* **`CRASHED`** iff any assertion, control, orbit-count assertion or
  instrument invariant fails, or an unhandled exception aborts a stage; all
  partial artifacts are preserved, nothing deleted.
* `REJECTED` is reserved for a lifecycle refusal; `SUPERSEDED` if withdrawn by
  a later preregistration; `REHEARSAL` does not apply.

## Evidence labels

* Exact orbit-count, census, distribution, invariance,
  direct-full-alphabet equality, Arb inequality, containment, status and
  checksum statements: `MACHINE-VERIFIED` once their checks pass.
* Float BA locator values, wall/CPU/RSS figures and all cost projections:
  `COMPUTATIONAL-EVIDENCE`.
* Published Tavakoli-Nguyen-Bose formulas: `CITED-DEPENDENCY`, re-evaluated by
  the exact/Arb dependency.
* The Burnside arithmetic in this pre-statement: `MACHINE-VERIFIED` (exact
  integer arithmetic at write time).

## Inherited disclosures

* `log()` writes to stdout and `run.log`, so a broken stdout is fatal mid-row
  even though all evidence is already durable; kept frozen, excluded by the
  launcher rule above.
* The frozen largest-remainder `snap` carries its own documented cycling
  caveat plus a frozen 200-trial equality verification; inherited unmodified,
  any anomaly is a stop-and-report.
* `ru_maxrss` and wall figures are `COMPUTATIONAL-EVIDENCE`; CPU claims use
  in-process `time.process_time()` only.
* The direct full-input KL pass recomputes each word's exact subsequence
  dictionary; its cost scales as `~2^n` per word and dominates the budget. No
  shortcut is certified unless separately pre-registered; none is registered.

## Rule-7 scope sentence

The intended sweep is exactly `q=3, n=13` with `d in {1/2,1/5,1/10,1/20}`,
using `G=S_3 x C_2` on symbol values and reversal, exact total input-orbit
mass snap `2^30`, exact total output-orbit mass snap `2^30` with the frozen
full-support bump, fixed dual candidates `(ba_total_orbit_mass, uniform)`,
direct full-alphabet KL comparison over all `3^13 = 1,594,323` inputs, full
finite channel support, and Arb 400-bit outward rounding. **NOT swept:** any
`q=4` or `q>=4` row; `q=3` with `n!=13` (beyond the single pre-registered
`n=12` ACCEPT control); other deletion probabilities; `ba_word`;
non-`G`-invariant candidates; the falsified position/type reduction; the cause
of Tavakoli-Nguyen-Bose rounding; Morozov-Duman, Pinto-Ribeiro, and every
`n -> infinity` or asymptotic-capacity claim. Each certified row is only a
finite-`n` theorem, and an unreached row is only an unreached row.

Signed: DelcapN12 (assigned by Main). Lifecycle order: this prereg
path-scoped commit ->
`campaign.py init --gate delcap-q3-n13-total-output-orbit-mass-v1 --prereg
pre_statement_q3_n13.md` -> byte-identical prereg copy into the minted run dir
recording source commit + SHA-256 -> instrument manifest frozen at the final
runner bytes -> static guard (also the memory-feasibility probe) -> in-run
control battery -> main compute -> `campaign.py freeze` ->
`campaign.py close --verdict <one verdict>`.
