# PRE-STATEMENT — Two-rung corrected-core replay campaign

Written 2026-08-31T11:30:16Z before any rung-1/rung-2 evaluation.
Scope: two-rung replay with the ACCEPTED v11 interval_core (ce70f959…),
frozen in this campaign, replacing the shared v9/v10 interval core that was
live during the 2026-08-30T03:20:00 (rung 1) and 2026-08-30T05:20:00
(rung 2) Stage-(b) evaluations.

## Parent provenance

- Rung-1 source: `2026-08-30T03:20:00Z_c9d4e2a8/stage_b_rung1_certified.py`
  (alman25_2.371339 certified interval).
- Rung-2 source: `2026-08-30T05:20:00Z_e97c35ae_70c28fc61780/stage_b_rung2.py`
  (vxxz24_2.37155181 three-region single-p_comp verified interval).
- Interval core: canonical v11 `interval_core.py` (hash ce70f959…),
  partial from_two upper b.upper() fix, endpoint-based hulls, Arb-strict
  membership and __debug__ guard.
- The v9/v10 entropy campaign is NOT part of this scope; this campaign is
  pure Stage-(b) rung replay.

## Frozen inputs staged here (no live paths)

- `stage_b_rung1_certified.py` (live sys.path replaced with HERE)
- `stage_b_rung2.py` (live sys.path replaced with HERE)
- `alman25_float.py` (from rung-1 producer)
- `vxxz24_float.py` (from rung-2 producer)
- `interval_core.py` (from accepted v11)
- `W1.00_2.371339.mat` (rung-1 input)
- `K100_2.37155181.mat` (rung-2 input)
- `launcher.py` (frozen execution orchestrator)

All frozen files' sha256 hashes are recorded in `checksums.sha256` and
`launcher.py` snapshots the hashes alongside each run's output JSON.

## Launcher perimeter

- OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=
  NUMEXPR_NUM_THREADS=1 enforced and set in the subprocess environment.
- `__debug__` requirement per Main v11 D3; python -O refuses to run.
- Sequential only: rung1 → rung2, never parallel.
- Per-rung wall cap: 900 s via `Popen.communicate(timeout=900)`; the source
  code path also sets an internal SIGALRM at 900 s.
- Nice 10 is supplied by the launch command and asserted by the launcher.

## Rung-2 diagnostic-off (bounded invocation CONFIRMED)

Main audit round-2 confirms the frozen `stage_b_rung2.py` DOES expose
`without` as a separate argv-1 flag (lines 22-23, 270-273).  The
launcher invokes a THIRD run `stage_b_rung2.py without` with its own
900 s cap and its own reserved output path.  The with-Lemma result
participates in the rigorous ordering (rung1 < rung2_with and both <
2.371866); the diagnostic-off (without) run is recorded but does not
participate in the ordering.

Correction: “rigorous ordering” above is historical preregistration language
superseded by the computational-only scope in “Historical labels are SOURCE
SCHEMA ONLY” below. The launcher records a numerical ordering; it does not
promote that comparison to a formal omega claim.

## NO expected-equality gate to old endpoints

Main directive: "No expected-equality gate to old endpoints."  Therefore
this campaign does NOT pre-register `omega_cert_upper == 2.371339…` or any
specific numeric value as a pass/fail criterion.  The old certified values
are recorded HERE for comparison only (not asserted):
- rung1 frozen certified_result value: 2.371339 (OMEGA_PUB constant)
- rung2 frozen with_lemma1 value: 2.3715538358 (recorded in v9 archives)
The acceptance criterion per Main directive is strictly:
  - each old stage's own structural/containment guards pass under the
    corrected core (evaluated by the corrected-core run itself),
  - then compare rung1 < rung2 and both < 2.371866 AFTER both runs.

## Historical labels are SOURCE SCHEMA ONLY (Main round-3 correction)

Both overridden rung scripts contain DECISION-RELEVANT float
extractions (e.g. rung1's M_low/min float, R_low midpoint sign,
target/total_R float; rung2 analogous) and add `(cF_max+ceF_max)/M`
without a proved repair/Lipschitz theorem.  A corrected-core replay
faithfully measures the OLD-PROTOCOL outputs and their ordering but
CANNOT by itself certify an omega bound.

The scripts' historical `certified_result`/`omega_cert_upper`/
`omega_with_absorbed_defects` field names are preserved verbatim as
**SOURCE SCHEMA ONLY** — these are the names the old producer used, not
an assertion that the values are rigorous omega bounds under the
corrected core.  The corrected-core replay measures the old-protocol
outputs faithfully; it does not propagate their historical
certifications.

Formal two-rung / omega claims remain SUSPENDED.  Precise blockers:
1. outward end-to-end aggregate rewrite (all float extraction points
   replaced with Arb endpoint-pair semantics end to end), and
2. mathematical justification of the feasibility absorption
   `(cF_max+ceF_max)/M` repair (a Lipschitz/other bound formally
   proving the absorbed defect covers the omitted true defect).

Do not silently harden the rung scripts in this campaign (per Main
directive).  This campaign is a MEASUREMENT campaign only.

### This campaign's ACTUAL acceptance result (COMPUTATIONAL corrected-core diagnostic)

What the launcher ACCEPTS as a positive run (its "PASS") is strictly:
  - each child's own structural/containment guards complete without
    aborting (the child's internal asserts pass), and
  - stdout yields a JSON dict that the launcher parses.
Then the launcher records, preserves, and COMPARES the old-protocol
outputs (per Main round-2/3 directives) and reports the comparison.
The result is a COMPUTATIONAL corrected-core DIAGNOSTIC — the campaign
does NOT assert any formal omega claim.

### Rule 7 scope discipline update

No new scope was opened.  The rule-7 sentence above is already correct.

## Zero real point compute (static freeze only)

This campaign does NOT run either stage during the freeze (ZERO compute
pre-statement: no evaluate_point, no entropy solve, no objective
evaluation).  Both stages will be run only when Main issues the CPU
boundary release.  The launcher is the frozen entry point for that
release.

## Rule 7 scope sentence

The replay covers only the two historically registered Stage-(b) runs
(rung 1 = alman25 2.371339 certified interval over W1.00; rung 2 = vxxz24
2.37155181 over K100; 3-region single-p_comp, max_level=3, q=5) with the
v11 corrected interval_core — no other direction, step, box, region or
parameter block, multi-block/cross-block class, q/max-level regime,
unpublished omega parameters, other construction, or other rung.
Comparisons to 2.371866 and to each other's outputs are the only numeric
judgments performed; nothing here propagates via the entropy-repair
campaign.  No paper or record claim is made.

## Owner static repair after delegated v5 report

The delegated v5 report said the multiline-JSON parser, raw-stdout capture,
atomic-write hardening, and resolved-path plumbing were on disk. They were
not: owner reread `launcher.py` and found the rejected per-line parser and
the old write path unchanged. Main repaired the executed source directly.
The accepted launcher SHA-256 is
`988055d274ef983d480abd264a357ae00d574f820f3aacfda7b9c47dfaafbc12`.

Static behavioral check, without a real rung evaluation: loading the actual
launcher through `runpy` parsed an indented trailing JSON object, rejected
trailing junk, and compiled the source. Delegated `__pycache__` artifacts are
preserved but are neither frozen payloads nor evidence. Real-child acceptance
remains reserved for the CPU-bound release.

## Owner release-perimeter refreeze v7

Before release, the owner found four non-mathematical defects in static v6:
there was no Main release gate; arbitrary nested campaign-local output names
were accepted; existing temp paths were not part of the all-output preflight;
and the six internally named relative logs resolved against the authoritative
command's `cwd=cs`, outside the campaign, so that command would refuse before
rung 1.

The launcher now requires `OMEGA_TWO_RUNG_RELEASED=1`, locks all ten parsed,
raw, stderr, and comparison artifacts to exact direct campaign paths, resolves
internal logs under the launcher's `HERE`, rejects final or temp paths with
`lexists`, creates JSON temps exclusively, refuses a target that appears
during the write, and `fsync`s each child stderr file and its parent
directory. The shell stdout remains protected by `noclobber`.

Accepted launcher SHA-256:
`3da22d832bdbd762807c896f1a4b48e33d543215ced1a711f86b801f53fbc236`.
These changes do not touch either child, the corrected interval core, the
comparison semantics, or the **SUSPENDED** formal status. No child or interval
computation ran during this refreeze.

Static owner checks passed: source compilation; all nine frozen file/parameter
hashes; all ten exact reserved final/temp paths under the release environment;
and all four authoritative `cwd=cs` argv paths. No child, interval evaluation,
or campaign artifact was created.

The final release command also sets `PYTHONDONTWRITEBYTECODE=1`, preserving
the delegated historical `__pycache__` files without creating or replacing
runtime bytecode.

## Released-run abort and preregistered comparison recovery

The exact released command completed all three child invocations and wrote
their parsed/raw/stderr artifacts, then the launcher aborted before comparison
with `KeyError: 'omega_cert_upper'`. The frozen `_compare` expected a rung-1
field that the equally frozen rung-1 producer does not emit. No
`comparison_summary.json` exists. All child evidence and the abort log are
preserved; the launcher and existing outputs will not be edited, deleted, or
rerun.

The child outputs expose two distinct rung-1 numbers,
`omega_cert_upper_raw` and `omega_cert_with_absorbed_slack`. Choosing one
silently would change the failed preregistration. The bounded recovery
therefore records **both** comparisons against rung 2's
`omega_with_absorbed_defects` and against target `2.371866`, plus the
with/without-Lemma rung-2 delta. It makes no replacement claim for the missing
key and no formal omega claim.

Before that comparison runs, `recover_comparison.py` will be frozen with hard
SHA-256 expectations for the three parsed child artifacts and the abort log,
parse every JSON number as `Decimal` from its exact text, assert the missing
key and child schemas, require `OMEGA_TWO_RUNG_RECOVERY_RELEASED=1`, and write
only `comparison_summary_recovered.json` by exclusive temp, file/directory
`fsync`, and atomic rename. Formal status remains **SUSPENDED** for the two
pre-existing mathematical blockers. This is post-abort computational
diagnostic recovery only.

Frozen recovery verifier SHA-256:
`6a9121f5e73115727fb5a29ca0186e11d495ea9cd2956312cf0301a6922e5957`.
The released command is exactly
`OMEGA_TWO_RUNG_RECOVERY_RELEASED=1 PYTHONDONTWRITEBYTECODE=1 python3 recover_comparison.py`.

## Recovered comparison result

The write-once recovery completed. `comparison_summary_recovered.json`
SHA-256 is
`8160e846c4b64a167624dad8460e83bc29c4e69ac16115a7015da208a838adda`.
Rung 1 raw/absorbed values are `2.371340083602922` and
`2.3713411672715115`; rung 2 with Lemma is `2.3715538358544617`.
Both rung-1 comparisons are true, with gaps `0.0002137522515397` and
`0.0002126685829502`, respectively. All three values are below target
`2.371866`. The rung-2 without-Lemma diagnostic is
`2.3715518062057623`; with-minus-without is `0.0000020296486994`.

This recovered artifact is **COMPUTATIONAL-EVIDENCE** after the frozen
launcher's aggregate abort. It deliberately leaves the missing key
unresolved, replaces no frozen artifact, and leaves formal status
**SUSPENDED**.
