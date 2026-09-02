# Pre-statement — q=3, n=11 total-output-orbit-mass extension

Campaign ID:
`2026-08-31T20-34-52Z_04ed25bd-455f-4814-8288-5ee0106e0db8_q3-n11-total-output-orbit-mass`.

Status at this first-file freeze: **PRE-RESULT / ZERO BA ITERATIONS / ZERO CERTIFICATE ROWS**.

## Question and adjudication

The immediately preceding q=4 campaign explicitly excluded `q=3,n=11` and named it as the next separate extension. This campaign asks one finite question: under the corrected exact total-input-orbit-mass primal and total-output-orbit-mass dual protocol, can each of the four registered `q=3,n=11` deletion-channel rows be enclosed with a valid two-sided Arb certificate, and does that interval strictly improve either or both Tavakoli-Nguyen-Bose finite-n sandwich endpoints?

A row is `CERTIFIED` only if its exact primal and admitted exact dual form a nonnegative-width interval whose conservative outward width is at most `1/500` bits/symbol. A mathematically valid wider interval is `VALID_BUT_WIDTH_MISS`. A strict lower improvement requires the certified primal lower endpoint to exceed the published `LB+` upper endpoint. A strict upper improvement requires the published `UB` lower endpoint to exceed the certified dual upper endpoint. A row with neither strict inequality is reported honestly as `NO_STRICT_IMPROVEMENT`; it is never tuned away.

## Frozen row box and order

Exactly four rows, in this order:

1. `(q,n,d)=(3,11,1/2)`;
2. `(3,11,1/5)`;
3. `(3,11,1/10)`;
4. `(3,11,1/20)`.

No other row may be added to this campaign.

## Frozen mathematical protocol

- Symmetry group: `G=S_3 x C_2`, acting on symbol values and whole-word reversal. Position permutation is forbidden.
- Float locator: 4,000 Blahut-Arimoto iterations in orbit coordinates. Float values locate candidates only and are `COMPUTATIONAL-EVIDENCE`.
- Primal: snap total input-orbit masses to denominator `2^30`; for input orbit `O`, set `p_x=mass(O)/(2^30 |O|)`; induce the full output law exactly.
- Dual candidates, in fixed order: `(ba_total_orbit_mass, uniform)`. Snap total BA output-orbit masses to denominator `2^30`; apply the frozen full-support bump if any orbit mass is zero; for output orbit `T`, set `D_y=mass(T)/(S |T|)`. `ba_word` is forbidden.
- Candidate selection: minimum finite direct-full-alphabet Arb upper endpoint, with the fixed order above as the tie order.
- Arithmetic: exact integers/rationals for distributions and channel weights; 400-bit Arb for logarithms, mutual information, KL maxima, published bounds, margins, and interval endpoints.
- Presentation endpoints: move one binary64 ULP outward from the exact Arb endpoint, archive the resulting exact binary rational, and assert containment. Binary64 locator and display values are not certificate endpoints.
- Width target: `1/500` bits/symbol.

Every admitted row must check, before reporting a bound:

1. exact input and output normalization;
2. positivity/full support where required;
3. expanded `S_3 x C_2` generator equality;
4. direct full-input primal agreement with the orbit-weighted primal;
5. direct full-input KL maximum agreement with the representative maximum for both fixed dual candidates;
6. nonnegative certified width;
7. containment relative to the published finite-n sandwich;
8. conservative outward endpoint and width containment.

## Frozen dependencies

The runner must hash these immutable prior-campaign bytes before importing any of them:

- `../2026-08-31T09:59:17Z_b6cd7315-caf3-4225-bb78-e4a81bf9a0f7_q4-total-output-orbit-mass/gate_c_orbit_dependency.py.asrun`, SHA-256 `d17c2b438840bd91e84b775e096d43c60c419ade8ff2ee1ff6e6095694eca9f9`;
- `../2026-08-31T09:59:17Z_b6cd7315-caf3-4225-bb78-e4a81bf9a0f7_q4-total-output-orbit-mass/gate_c_dhalf_dependency.py.asrun`, SHA-256 `9bee55a5e137ae065a83271586284fa9624cffbfbc7c8d416f7a62ad86a71db9`;
- `../2026-08-31T09:59:17Z_b6cd7315-caf3-4225-bb78-e4a81bf9a0f7_q4-total-output-orbit-mass/gate_c_q4_structure_dependency.py.asrun`, SHA-256 `6d12cc77cefe4ff7cce1efa7770e69cb5c470ca6797f66446756bd4d2b141cdd`.

The new runner and manifest will be frozen additively after this pre-statement. The manifest must bind the runner hash, all dependency hashes, this exact row box, the dual tuple, precision, snaps, caps, and schema before any BA iteration.

## Census anchor and pre-design disclosure

Before this campaign directory existed, owner ran one low-priority, thread-pinned **structure-only** design probe at `q=3,n=11`. It ran no BA iteration, produced no input/dual candidate, and wrote no result artifact. It measured the anchor tuple

`(input G-orbits, output G-orbits, representative entries, orbit slots, total sparse slots)`

as

`(14884, 22450, 6148309, 3573542, 9721851)`,

with input/output orbit-size sums `177147=3^11` and `265720=sum_{k=0}^{11}3^k`. Owner independently checked the orbit counts by Burnside's lemma: the `S_3 x C_2` input-orbit count is `14884`, and summing the same closed form over output lengths `0..11` gives `22450`.

A second structure-only design probe evaluated the already-frozen published-sandwich formula at the four registered deletion probabilities solely to confirm that its domain supports `n=11`. It computed no primal or dual capacity candidate. These two pre-design probes are **RETROSPECTIVE DESIGN EVIDENCE**, not prospective controls and not campaign results. The row box, protocol, candidate order, resource caps, status rules, and all BA/result work are frozen here before their first execution.

## Static controls required before release

The frozen runner's `--static-guard` mode may build exact orbit structures but must run zero BA iterations and write zero row/summary/table artifacts. It must pass all of the following:

- manifest and all dependency hashes match before dependency import;
- the row box is exactly the four rows above and excludes q=4 and every `n!=11` cell;
- independent Burnside counts equal the measured `14884` input and `22450` output orbits;
- the full `q=3,n=11` structure matches the frozen census tuple and alphabet-size sums;
- exhaustive small-instance joint value/reversal equivariance passes, while a position-permutation plant breaks;
- a planted split output distribution is rejected before representative compression;
- an exact orbit-uniform output control is admitted;
- a zero-support dual returns infinity;
- a negative-width plant is rejected;
- empty resume validation reports exactly zero rows.

No failed control may be weakened, removed, or reinterpreted inside this campaign.

## Execution, resource, and crash contract

One sequential process only. Production requires:

- `DELCAP_Q3_N11_RELEASED=1`;
- process niceness at least 10;
- `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=NUMEXPR_NUM_THREADS=1`;
- `PYTHONDONTWRITEBYTECODE=1`;
- Python assertions enabled.

Caps: 96 GiB RSS, 90 minutes per row, six hours for the four-row stage. A resource stop becomes an explicit `NOT_REACHED_RESOURCE` row; every later cell becomes `NOT_REACHED_RESOURCE_UPSTREAM`. The campaign then finalizes as `FAILURE_TO_CERTIFY_RESOURCE`; it is not silently resumed under changed caps.

Canonical rows are newline-terminated canonical JSON. The row file and its per-line SHA-256 ledger are each written through a temporary file, flushed, `fsync`ed, and replaced atomically. The pair is not crash-atomic: any mismatch, malformed prefix, wrong ordering, wrong source hash, or preserved `.*.tmp.*` file hard-aborts resume and is never repaired or deleted automatically. Terminal `FAILURE.json`, `summary.json`, or table evidence forbids rerun. No existing artifact may be overwritten.

The exact production command, after static acceptance, is:

```sh
env PYTHONDONTWRITEBYTECODE=1 \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  DELCAP_Q3_N11_RELEASED=1 \
  nice -n 10 ./.venv/bin/python \
  delcap/campaigns/2026-08-31T20-34-52Z_04ed25bd-455f-4814-8288-5ee0106e0db8_q3-n11-total-output-orbit-mass/gate_c_q3_n11_total_orbit.py.asrun \
  --run
```

The first registered `d=1/2` row is the full startup anchor and, if accepted, is also canonical row zero; it is not recomputed in a separate anchor mode.

## Evidence labels

- Exact census, distribution, invariance, direct-full-alphabet equality, Arb inequality, containment, status, and checksum statements: `MACHINE-VERIFIED` after their checks pass.
- Float BA locator values, wall/CPU/RSS figures, and cost estimates: `COMPUTATIONAL-EVIDENCE`.
- Published Tavakoli-Nguyen-Bose formulas: `CITED-DEPENDENCY`, re-evaluated by the exact/Arb dependency.
- The two pre-directory probes above: `RETROSPECTIVE DESIGN EVIDENCE` only.

## Rule-7 scope sentence

The intended sweep is exactly `q=3,n=11` with `d in {1/2,1/5,1/10,1/20}`, using `G=S_3 x C_2` on symbol values and reversal, exact total input-orbit mass snap `2^30`, exact total output-orbit mass snap `2^30`, fixed dual candidates `(ba_total_orbit_mass, uniform)`, direct full-alphabet KL comparison, full finite channel support, and Arb 400-bit outward rounding. **NOT swept:** any q=4 row; q>=4; q=3 with n!=11; other deletion probabilities; `ba_word`; non-G-invariant candidates; the falsified position/type reduction; the cause of published rounding; Morozov-Duman, Pinto-Ribeiro, and every `n -> infinity` or asymptotic-capacity claim. Every accepted row will be only a finite-n theorem.

Signed: Main. No BA iteration or q=3,n=11 certificate row existed when this first file was written.
