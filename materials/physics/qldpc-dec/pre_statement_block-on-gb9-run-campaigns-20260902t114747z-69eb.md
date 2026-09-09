# Pre-statement — re-launch GB9 (beam64_32res_640iters, num_results=32) on mini-0, hold IN-FLIGHT

Cycle 20260908T193320Z_ba7b32 · physics/qldpc-dec · recorded 2026-09-08, before
`campaign.py init` and before any decode of this run.

## Gate question

Inherited unchanged from `pre_statement.md` Revision GB9 (recorded
2026-09-02T10:56Z, pre-init, pre-sample): does arXiv:2512.07057 Table 1 row 4
`beam64_32res_640iters` = (beam_width 64, initial_iters 40, iters_per_round 30,
max_rounds 20, **num_results 32**) achieve the published **17x** LER
improvement over BP-OSD at p=1e-3 on the [[144,12,12]] gross code, as
reimplemented from published Algorithm 3?

This cycle's bounded objective is narrower than answering it: the previous
GB9 run (`campaigns/20260902T114747Z_69eba724_61e9fe30ba5d`) died silently
twice on studio.local and was closed **CRASHED** 2026-09-08T16:29:10Z with 4
of 20 beam64_32res shards decoded. The objective here is to **re-launch the
identical protocol as a fresh run on mini-0.local and hold it IN-FLIGHT with
a live claim and first-output validation**. The gate question is answered at
close, not this cycle.

## Design (inherited, frozen)

- Same runner `src/gateb_gb9_32res.py` (sha256 `3a2e2800…250ab` — verified
  byte-identical to the crashed run's snapshot), same pinned binary
  `src/beam_cpp/beam_nr_cpp` (sha256 `297db977a57b…`), same arms, same frozen
  GB7 stream (`sampling_seed(20260829, 1e-3, "Z", "bp30+osd")`, N=1e8), same
  frozen decision rule (`decide()`, Revision GB9).
- One deliberate change, **execution site only**: mini-0.local (Apple M4, 10
  cores, 16 GB, macOS 26.4.1) instead of studio.local, `--threads 10`. Host
  and thread count are manifest metadata, not gate quantities: every hard
  void is a sha256 equality check against frozen GB7 artifacts and is
  host-independent.
- **Declared artifact inheritance from the crashed run** (recorded
  pre-launch): `shots_bposd.bin` + `sample_inventory.json` (driver re-hashes
  the stream and asserts equality with the frozen GB7 sha); all 20 beam8
  shard predictions + metas (driver idempotency re-hashes every predictions
  file and checks binary sha, arm, source shard sha, and config line; the
  assembled beam8 predictions must then equal the frozen GB7
  `beam8_predictions.bin` byte-for-byte — the 1e8-shot K=1 regression, fully
  re-checked); beam64_32res shard predictions + metas for shards 00–03 (same
  per-shard validation). Shard `shots_XX.bin` files are **not** inherited:
  `write_shards` re-slices them deterministically from the verified stream on
  mini-0. Any hash drift aborts the run before any new decode — inheritance
  can only save compute, never change an outcome.
- No new sampling decisions; no decoder, circuit, DEM, source, or pre-statement
  changes.

## Acceptance criteria (falsifiable; this cycle = launch objective)

- **L1.** `campaign.py init` mints the run dir with this file's sha256;
  claim held on the workspace control plane.
- **L2.** Remote launch on mini-0 detached per the workspace runbook:
  `nohup caffeinate -s .venv/bin/python src/gateb_gb9_32res.py --run-dir
  campaigns/<run_id> --threads 10 --shards 20 > ~/gb9/gb9-run.log 2>&1 &`;
  exact command, remote PID, and log path recorded in ITERATE.md.
- **L3.** First-output validation, all in `gb9-run.log` / shard files:
  `shared stream: reused (hash verified)` and `stream identity with frozen
  GB7: PASS`; all 20 beam8 shard metas validated with zero shard re-decodes;
  `beam8: 223 failures; identity with frozen GB7: PASS` (223 = frozen GB7
  beam8 count); beam64_32res metas for shards 00–03 validated; shard 04
  decode starts with `beam_nr_cpp` visible to `pgrep` and the shard stderr
  showing exactly `config: beam_width=64 initial_iters=40 iters_per_round=30
  max_rounds=20 num_results=32`.
- **L4.** The process survives ssh disconnect (nohup) and holds
  PreventSystemSleep (caffeinate -s): verified by `pgrep` from a fresh ssh
  session after launch.

Failure of any L1–L4 is a failed launch: record the obstruction verbatim and
stop; that is a result, not a workaround opportunity.

## Verdict mapping (applies at close, per Revision GB9 — unchanged)

- `GAP_CONFIRMED` (width effect AND CI excludes 1 AND CI excludes 1/17) →
  close `FROZEN-NEGATIVE`.
- 1/17 not excluded AND improvement over the frozen GB7 beam64_640iters mask
  confirmed (McNemar p<0.05, b>c) → close `FROZEN-CERTIFIED`.
- Anything else at N=1e8 → close `FROZEN-INCONCLUSIVE`.
- Process death before the driver's `campaign.close()` → `CRASHED` (no
  verdict; this is the failure mode this re-launch engineers against with
  nohup + caffeinate + heartbeat monitoring).

**This cycle closes nothing.** The run stays open (IN-FLIGHT) with a live
claim; close happens when the decode completes — estimated ~3 weeks on
mini-0 (16 remaining shards, ~15 h/shard measured at 26 threads on the M3
Ultra; M4 has 10 cores).

## Scope boundary

- Launch + first-output validation only. No freeze, no close, no usage
  record this cycle (the run has not ended; accounting at close or by
  settle's estimator).
- The sealed crashed run dir is a read-only source of inherited artifacts;
  never written.
- No changes to decoder sources, binaries, circuits, DEMs, or the frozen
  GB7 campaign; X basis, p≠1e-3, and the authors' undisclosed heuristics
  remain out of scope exactly as Revision GB9 bounded them.
- All control-plane operations (`campaign.py`, state.json) run locally in
  the workspace; all decode compute on mini-0.local only (never
  mini-pro.local).
- `--threads 10`, the driver's own `os.nice(2)`; mini-0 is shared with the
  qlops short rebuild, which is capped/niced — this job has priority.

## Verification command (exact)

```bash
ssh mini-0.local 'pgrep -fl beam_nr_cpp; tail -n 40 ~/gb9/gb9-run.log'
```

plus claim keep-alive, run locally from `physics/qldpc-dec`:

```bash
python3 /Users/jinleic/jinleic-workspace/scripts/campaign.py heartbeat
```
