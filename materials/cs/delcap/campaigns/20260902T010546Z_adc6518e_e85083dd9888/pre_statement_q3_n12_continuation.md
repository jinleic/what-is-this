# Pre-statement — q=3, n=12 rows 2-3 continuation (d=1/5 successor)

Written and committed BEFORE any continuation compute. Successor to campaign
`20260901T121346Z_d1db3155_4354a957b8ba` (gate `delcap-q3-n12-total-output-orbit-mass-v1`,
closed `CRASHED`), whose frozen bytes remain immutable and whose two
CERTIFIED rows (0: `d=1/2`, 1: `d=1/5`) are inherited evidence frozen at
`orbit_rows.jsonl` SHA-256 `3c5cf81a91512eeacdf94fc0131bf2564068ff53cd875734b6195e1514cf8bb8`.

## Crash provenance (exact)

The parent campaign's production runner committed ordinals 0 and 1, then died
on ordinal 2 with `BrokenPipeError(32)` inside `log()` at 2026-09-01T14:05:14Z,
writing `FAILURE.json` which, by the frozen parent rule, forbids rerun inside
that campaign. Root cause is the launcher, not the instrument: the parent's
launch piped the runner's stdout through `tail -10` (launcher deviation from
the frozen command, which does not pipe), and the OS closed the read end mid
run; the runner's own semantics then correctly hard-aborted the stage and
archived `FAILURE.json` with `completed_rows: 2`. The instrument arithmetic is
not implicated by this crash: rows 0 and 1 are valid, complete, exact
certificates, byte-validated by the parent's own ledger.

## Question (narrowed)

Complete the pre-registered four-row box: certify rows 2 (`d=1/10`) and 3
(`d=1/20`) of the SAME box `(q,n)=(3,12)`, `d in {1/2,1/5,1/10,1/20}`, under
the identical frozen protocol. Rows 0-1 are NOT recomputed; they are inherited.

## Frozen row box and order

Exactly two rows, in this order:

1. `(q,n,d)=(3,12,1/10)` (ordinal 0 of THIS campaign; box ordinal 2);
2. `(3,12,1/20)` (ordinal 1 of THIS campaign; box ordinal 3).

## Frozen mathematical protocol

Byte-identical to the parent preregistration
(`cs/delcap/pre_statement_q3_n12.md`, sha256
`d3b2b01493928e486f9c2bcce06aae20c93e3408cdca1e60010118a7e07e7081`): same
group `S_3 x C_2` (values + reversal), same exact snaps `2^30`, same dual
tuple `(ba_total_orbit_mass, uniform)` with `ba_word` forbidden, same
400-bit Arb outward interval construction, same width target `1/500`, same
verdict semantics per row, same STOP-AND-REPORT orbit-count assertions
(`44,530` input orbits over `531,441` words; `66,980` output orbits over
`797,161` words; census tuple
`(44530, 66980, 30666848, 18311678, 48978526)`), same exact orbit-mass
identity checks, same frozen dependencies. The runner is the SAME frozen
bytes (SHA-256 `08836f43adbffcbf35da90feb140ca030e387605043802243dd3dac4e6cc6718`)
with three mechanical host-side changes only:

1. `CAMPAIGN_ID`/`CAMPAIGN` point at the new run dir;
2. `SCHEMA` gate tag string becomes
   `delcap-q3-n12-continuation-d10-d20-v1`;
3. the row box enumerates only the two remaining `d` values, and ordinal
   labels carry the global box ordinal (2 and 3) alongside.

No arithmetic, snapping, candidate, comparison, or interval logic may change:
production re-executes `certificate_row` from the identical source, whose
hash remains `08836f43adbffcbf35da90feb140ca030e387605043802243dd3dac4e6cc6718`
(apart from the mechanical host constants above, re-frozen in the successor
instrument manifest).

## In-run controls (both directions, as in the parent)

* ACCEPT (known-true): the frozen parent row `0:3:12:1/2` (this continuation's
  own runner) is re-validated in-run by hash against the parent ledger — and,
  because the two rows here re-run the identical `certificate_row` code on the
  identical structure, a BIT-EXACT reproducer of row 0:3:11:1/2 (the parent's
  35-field ACCEPT control, already executed and archived with the parent's
  identical hash) is inherited by reference; its artifact
  `control_accept_n11.json` (parent-campaign bytes, hash-pinned) is copied
  unchanged into this campaign as inherited evidence.
* REJECT plants R1a/R1b/R1c/R2/R3 are re-fired in-run before row 2 with the
  same `q=3,n=5` control structure and the same required directions.
* The exact orbit-mass identity (input masses sum to exactly `2^30`; post-bump
  output numerators sum to `2^30 + z_bump`; per-word expansions sum to
  exactly `1`) is asserted inside every admitted row exactly as in the parent.

## Execution contract

Same resource budget as parent (96 GiB RSS, 90 min/row, 6 h stage wall,
6 CPU-hours stage budget, nice >= 15, single-thread pinning,
`PYTHONDONTWRITEBYTECODE=1`, assertions on, release key
`DELCAP_Q3_N12_CONT_RELEASED=1`). Exact production command (cwd `cs/`),
NO pipe on stdout — the crash root cause is deliberately excluded:

```sh
env PYTHONDONTWRITEBYTECODE=1 \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  DELCAP_Q3_N12_CONT_RELEASED=1 \
  nice -n 15 ./.venv/bin/python \
  delcap/campaigns/<THIS_ID>/gate_c_q3_n12_continuation.py.asrun \
  --run < /dev/null
```

stdout redirection to a file (not a pipe) is allowed. Immutable single-pass
row files with per-line SHA-256 ledger, atomic replacement, no-overwrite
terminal evidence, and hard-abort resume semantics are inherited exactly.

## Verdict rule

* `FROZEN-CERTIFIED` iff BOTH continuations rows are `CERTIFIED` (or valid
  width-miss rows are reported honestly alongside certified ones), all
  in-run controls pass in both directions, and the orbit-count assertions
  hold; the closing report then states the four-row box outcome — which rows
  are certified and whether BOTH published TNB endpoints are strictly beaten
  per row — combining parent rows 0-1 (inherited, byte-pinned) with rows 2-3.
* `FROZEN-INCONCLUSIVE` iff the stage cannot complete under budget: the exact
  frontier (which rows stand certified, remaining cost) is frozen as a bounded
  result.
* `CRASHED` on any assertion/control failure or unhandled exception
  (`FAILURE.json` written; no silent reinterpretation).
* No asymptotic-capacity claim is registered or permitted; a bounded outcome
  is a bounded outcome.

Signed: DelcapN12. Zero continuation BA iterations existed when this file was
written.
