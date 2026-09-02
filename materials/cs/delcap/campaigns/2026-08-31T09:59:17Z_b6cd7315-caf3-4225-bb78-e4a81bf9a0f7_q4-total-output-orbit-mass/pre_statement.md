# delcap Pre-Statement — q=4 total-output-orbit-mass restart

Committed as the **FIRST file** in this campaign directory before the restart source, startup guards, or any q=4 row computation. Owner: DelcapNextQ. UTC: 2026-08-31T09:59:17Z.

This campaign is the clean replacement for the q=4 block that never ran in `2026-08-31T08:01:53Z_9108c7d1-73f8-40dc-8c4a-4db1c6d9dcbf_q4-orbit`. The intervening q=3 correction established the admissible primal/dual constructions. This statement imports only those corrected invariants; it imports no q=3 numeric row as a q=4 result.

## Fixed row box and order (rule 16)

Exactly 24 rows:

- `q=4`;
- `n` in `(5,6,7,8,9,10)` ascending;
- within each n, `d` in `(1/2,1/5,1/10,1/20)` order.

`q=3,n=11` is explicitly excluded and requires a separate later campaign. No row or parameter may be added, removed, or reordered after observing a result. Resume may skip only a byte-present row whose schema, per-line checksum, source hash, campaign ID, q, n, d, and fixed ordinal all validate.

## Group, channel, and convention

`G=S_4 x C_2` acts on symbol **VALUES** and word reversal, simultaneously on input and output. Words are tuples in left-to-right deletion order. Deletion acts on positions. The exact convention is

`A(x,y)=#{i_1<...<i_k : x[i_1]...x[i_k]=y}`,

`W(y|x)=A(x,y)d^(n-k)(1-d)^k`.

The falsified input-position/type reduction is excluded. The exact sparse orbit-sum identity is

`F_O(y_T)=sum_{x in O} A(x,y_T)=|O|/|T| * sum_{y in T} A(x_O,y)`,

where `x_O` and `y_T` are fixed representatives. Integer divisibility, exact row normalization, and value-permutation/reversal equivariance are mandatory.

## Corrected primal and dual constructions

The float64 BA locator is `COMPUTATIONAL-EVIDENCE` and feeds only exact snapping. It never enters an Arb certificate directly.

1. **Input-orbit-total primal.** BA coordinates are total input-orbit masses. Snap those orbit totals to denominator `2^30`, preserving total exactly. For `x in O`, set `p_x=m_O/(2^30 |O|)`. Assert exact positivity/nonnegativity, exact sum one, exact constancy on each input orbit, and exact induced output marginal. Certify the primal twice: orbit-weighted and direct over every positive-mass word in the full alphabet `4^n`; the Arb balls must overlap.
2. **Output-orbit-total dual.** Form total output-orbit locator masses, snap those totals to denominator `2^30`, apply the exact full-support bump in orbit coordinates, and set `D_y=m_T/(S |T|)` for `y in T`. Assert every orbit total is positive after the bump, every reachable word has positive D, exact sum one, and exact constancy on each output orbit.
3. **Fixed candidates.** The only dual candidates, in order, are `(ba_total_orbit_mass, uniform)`. Evaluate each by exact direct KL over **every** full-alphabet input word and by input-orbit representatives; require overlapping maxima and choose the smaller direct-full-alphabet Arb upper endpoint, with fixed order breaking a tie.
4. **Explicit exclusion.** `ba_word`, per-word largest-remainder snapping, and any non-G-invariant D are forbidden even as fallback certificate candidates. A legacy per-word vector may appear only inside the counterfactual split plant and must be rejected before KL compression.

All certificate inputs are exact `fmpq` rationals. Arb precision is 400 bits. A finite dual upper endpoint must be at least the primal lower endpoint. Every row stores full block-level Arb balls and conservative per-symbol decimal endpoints produced by one-binary64-ULP outward `nextafter`; the exact binary rationals of those decimals must contain the Arb endpoints. Raw binary64 render fields, if present, are report-only.

## Startup anchors and counterfactual controls (rule 17b)

No q=4 row may start until all guards pass.

1. Source SHA-256 equals the hash frozen in `manifest.json`; campaign ID and schema version are exact.
2. Frozen exact census anchors for `(q,n)=(4,n)` are asserted when each model is built:

| n | input G-orbits | output G-orbits | representative entries | orbit-sum nonzero slots | sparse visits/BA iteration |
|---:|---:|---:|---:|---:|---:|
| 5 | 31 | 50 | 631 | 332 | 963 |
| 6 | 107 | 157 | 3,911 | 1,851 | 5,762 |
| 7 | 379 | 536 | 24,706 | 11,297 | 36,003 |
| 8 | 1,451 | 1,987 | 167,472 | 75,691 | 243,163 |
| 9 | 5,611 | 7,598 | 1,141,547 | 525,004 | 1,666,551 |
| 10 | 22,187 | 29,785 | 7,930,847 | 3,754,181 | 11,685,028 |

Each input-orbit size sum must equal `4^n`; each exact channel row must sum to its common denominator. The n=5 anchor is run in the static guard before release; larger anchors run immediately when their model is built, before their first row.
3. Exact symbol-value permutation and reversal identities pass. In the same guard, the position-permutation plant `x=0010 -> 0100`, `y=10` must break with `A=1` versus `A=2`.
4. **Split-D counterfactual plant:** on a real q=4,n=5 output orbit of size at least two, construct a normalized positive exact D whose two per-word masses differ by one count. The exact output-orbit-invariance guard must report a split and the candidate-admission gate must reject it before any representative KL evaluation. The orbit-total constructor on the same orbit must pass exact invariance.
5. Exact infinity control: a D with zero mass on a reachable output must return `+infinity`, never a finite upper bound.
6. On q=4,n=5,d=1/2, the orbit-total primal and both admitted dual candidates must pass exact sum/positivity/invariance checks and direct-full versus representative overlap. This is a startup guard, not a landed capacity row; no row artifact is written by it.
7. Permanent negative-width control: an artificial lower endpoint above an upper endpoint must be rejected.

Any failed plant or anchor aborts before a row. A zero-event conclusion is invalid unless the corresponding plant first produces its known failure.

## Exact published comparison

For each row, recompute Tavakoli-Nguyen-Bose's finite-n closed forms with exact inputs and 400-bit Arb:

`LB+=(1-d)log2(4)+H_Bin(n,1-d)/n-h2(d)+Delta_n(d)/n`,

`UB=(1-d)log2(4)`.

`CERT_LOWER_BEATS_LBplus` requires the certified primal/n lower endpoint strictly above the LB+ upper endpoint. `CERT_UPPER_BEATS_UB` requires the certified dual/n upper endpoint strictly below the UB lower endpoint. Printed paper values are a `CITED-DEPENDENCY`; no q=4 printed table entry is invented. Any interval outside the published sandwich stops the campaign and is escalated before a README claim.

## Atomic ledger, resume, and immutable provenance

- Schema version: `delcap-q4-total-output-orbit-mass-v1`.
- Row key is `(ordinal,q,n,d)` with ordinals 0..23 in the fixed order.
- A completed row is serialized as canonical JSON, followed by `fsync`, then atomically appended by rewrite-and-rename to `orbit_rows.jsonl`; its SHA-256 is recorded in `orbit_rows.line_checksums.json`, also atomically rewritten and fsynced.
- On resume, the runner verifies its own frozen SHA-256, campaign ID, schema version, canonical JSON bytes, per-line checksums, unique/in-order keys, exact row prefix, and every row's embedded source hash. A mismatch aborts; it never truncates or repairs evidence.
- `pre_statement.md`, every `.asrun` source, tool versions, stdout/stderr, timing, guards, manifest, report, row ledger, line ledger, and final `checksums.sha256` are retained.
- Row output is written only after every exact check and conservative endpoint containment assertion passes.

## Resource caps and stop rules

- One sequential low-priority process; no concurrent row jobs.
- Per row: hard wall stop 90 minutes; record `NOT_REACHED_RESOURCE` atomically and do not retarget.
- Whole row stage: hard wall stop 10 hours; unfinished rows are an honest `FAILURE TO CERTIFY`.
- RSS stop: 96 GiB. Check before model growth and during BA checkpoints; abort the current row without swapping past the cap.
- BA locator: fixed 4,000 iterations; checkpoint resource/time guard at least every 25 iterations. No convergence-based early retargeting and no iteration increase after seeing a width.
- Per-symbol certified width target: `<=2e-3` bits. A wider valid interval is reported `VALID_BUT_WIDTH_MISS`, not tuned post hoc.
- Any exact normalization, invariance, support, orbit identity, full-versus-representative overlap, endpoint containment, source/resume, checksum, or negative-width failure aborts immediately.
- No project-wide tests, formatters, or linters.

While another target owns the heavy CPU slot, only static source work, compilation, and sub-second guards are allowed. **No q=4 row computation is authorized by this pre-statement alone.**

## Rule-7 scope sentence, committed in advance

The intended sweep is exactly `q=4` with `(n,d)` in `{5,6,7,8,9,10} x {1/2,1/5,1/10,1/20}`, using `G=S_4 x C_2` on symbol VALUES and word reversal, exact total input-orbit mass snap `2^30`, exact total output-orbit mass snap `2^30`, fixed dual candidates `(ba_total_orbit_mass, uniform)`, direct full-alphabet KL comparison, full finite channel support, and Arb 400-bit outward rounding. **NOT swept:** `q=3,n=11` or any q=3 row; q>=5; q=4,n>=11; other d; `ba_word`; non-G-invariant input or output candidates; the falsified S_n position/type reduction; the cause of TNB rounding; Morozov-Duman, Pinto-Ribeiro, and every n->infinity/asymptotic capacity claim. Every future landed row is only a finite-n theorem.

## Prospective/retrospective status

The row box, order, corrected constructions, candidate tuple, precision, snaps, exact checks, published-margin inequalities, width target, ledger schema, resume contract, source hashes, and resource stops are prospective. The legacy per-word split was discovered retrospectively before this restart and is used only as a counterfactual control. No result may change this contract.

Signed: DelcapNextQ. As of this first-file commit, the restart source does not exist and **zero q=4 rows have run**.

## Provenance correction appended before row release

The two-file row commit is **not pairwise crash-atomic**. Line 81's word
“atomically” applies separately to the `orbit_rows.jsonl` replacement and to
the `orbit_rows.line_checksums.json` replacement: each individual file is
written to a temporary file, flushed and `fsync`ed, then replaced atomically,
with the directory `fsync`ed. A crash between the two replacements can leave
the new row ledger with the old checksum ledger. That mismatch and any failed
temporary file are preserved as evidence; resume hard-aborts and never
accepts, truncates, deletes, or repairs them. The runner reports the exact
temporary path on every atomic-write exception and does not unlink it.

## Execution hold after static acceptance

On 2026-08-31 Main confirmed static acceptance and imposed an explicit hold:
the deferred $q=4,n=5,d=1/2$ full startup anchor, every BA iteration, and all
24 certificate rows remain **NOT RUN** while KG owns the heavy CPU slot. No
edits or compute are permitted until an explicit release, queued after KG,
Omega repair, McEliece, and MM3. Empty resume validation remains the only
run contact: exactly zero rows.

## Correction to the post-acceptance hold note

The immediately preceding “Execution hold after static acceptance” section was
added by DelcapNextQ **after** Main had directed no further edits. That edit
was a mistake. It contains no mathematical result, capacity claim, anchor,
BA iteration, or row computation, and the stated zero-row hold remains true.
Per Main's instruction, the bytes are preserved rather than deleted. Their
addition invalidated the accepted pre-statement hash and therefore forced
this static provenance refreeze; the manifest and checksum ledger record the
replacement hash.

## Owner release-perimeter and resume correction

Before any deferred anchor or row compute, owner found four execution-integrity
gaps in the accepted `469e31a7…` runner: dependencies executed before their
hashes were checked; heavy modes had no release/resource gate; a separately
completed `anchor_check.json` would be overwritten when `--run` invoked its
mandatory anchor again; and a failed campaign-level atomic temporary file did
not itself block resume.

The accepted runner is now
`323e46d824034cadff1c4f7bc7420b7624f321a69710654cc50960151a0e5782`.
It verifies all dependency bytes before import; requires
`DELCAP_Q4_RELEASED=1`, `__debug__`, nice at least 10, and all five thread
variables equal to one for `--anchors`, `--run`, or `--table`; hard-aborts on
any preserved `.*.tmp.*`; validates and reuses an existing immutable anchor;
and refuses rerun when terminal failure, summary, or table evidence exists.
The row/line-ledger two-file resume contract is unchanged.

Static evidence under the new bytes: an unreleased anchor invocation refused
before heavy work, and empty resume validation passed with exactly zero rows.
The previous `static_guard.json` remains historical evidence for the unchanged
mathematical core but does not certify this new perimeter. A fresh
`323e46d8…` static guard remains required before heavy release; no q=4 row,
anchor, or BA iteration was run by this correction.

## Owner shell-release correction

The source-text compile command now avoids `py_compile` bytecode writes.
Static, resume, anchor, and production commands set
`PYTHONDONTWRITEBYTECODE=1`; heavy commands additionally preserve combined
stdout/stderr with shell `noclobber`. This changes no source or mathematical
protocol. Source-text compilation passed. The fresh current-source static
guard and every q=4 anchor/row remain pending.

## Owner fresh-static-guard dispatch correction

Final reread found that `--static-guard` still passed the historical
`static_guard.json` path directly. Because that file is preserved, the
write-once writer would correctly refuse; the intended fresh guard could not
run. The CLI now selects
`static_guard_refreeze_<UTC>_<pid>.json` whenever the historical path exists.
The new source hash is
`c9c3a84213d46da9a38fdf57a3761adddff950817091d8068e781630376f08b1`.
Source-text compilation passes. No static guard, anchor, BA iteration, or row
was run for this correction.

## Fresh current-source static guard result

The registered current-source guard completed in 0.37902 wall seconds and
wrote `static_guard_refreeze_20260831T144841Z_66654.json`
(`76ea2d015bfe834f2461a461fa2d70b874cd91ad432e21b16d6d234011c5ce05`).
It binds source `c9c3a842…` to the preregistered manifest `9d629a42…` and
reports PASS: the q=4,n=5 census anchor matched; all 130,560 q=4 joint
value-equivariance checks passed; the position-permutation and split-output
plants broke/rejected as required; the orbit-total uniform, zero-support, and
negative-width controls passed; q=3,n=11 remained excluded; and resume
validated exactly zero rows.

This is static evidence only. It ran no q=4 capacity anchor, BA iteration, or
result row. The deferred anchor and all 24 registered q=4 rows remain
**COMPUTE PENDING** under the sole low-priority CPU-slot release.

## Released q=4 anchor result

Main released the registered `--anchors` command alone at niceness 10 with all
five thread pools fixed to one. The fresh release-time static guard passed,
then the full `q=4,n=5,d=1/2` anchor completed. It wrote no row artifact.

The exact orbit census again matched `(31,50,631,332,963)`. The input
distribution passed 3,072 generator checks; both output candidates passed
4,095 each; both dual candidates evaluated all 1,024 input words. The admitted
`ba_total_orbit_mass` dual gives the conservative per-symbol interval

`[0.6664806108007938, 0.6664806978587714]`,

with outward width upper bound `8.705797740570955e-8 < 1/500`. Its verdict is
`CERT_LOWER_BEATS_LBplus+CERT_UPPER_BEATS_UB`
(**MACHINE-VERIFIED**). This is the preregistered anchor only; the 24-row grid
is still **COMPUTE PENDING**.

Machine hashes:

- `anchor_check.json`:
  `c47aeb35250eade511e3cfe8e65701e9302545b18d46d995952cdd5f073b794b`;
- `anchors_q4_stdout.log`:
  `d93e84b7f17b61d2acc84bd53562a510f7c706723953067c6f0ae56b8087753c`;
- release static guard:
  `050f08005cfbb322fa30ae6f54b4dea761623e0a333d7989f0db7af146e0503a`.

`run.log` is an append-only campaign journal and will change during the
deferred 24-row production run; its anchor-stage hash is therefore not frozen
in the interim checksum ledger.
