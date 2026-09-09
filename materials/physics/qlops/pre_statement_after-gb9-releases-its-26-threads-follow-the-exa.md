# Cycle preregistration — physics/qlops — 20260908T174625Z_a9a475

hmz FREEZE artifact and campaign prereg for this cycle's single bounded
campaign. Written 2026-09-08T17:5xZ, BEFORE `campaign.py init` and before any
sampling.

## Gate question

Has the ungrown/grown zero-level CCZ author schedule — the frozen 12-cell,
71,000,000-shot d=3/d=9 reproduction pinned in `pre_statement.md`
Revision 7 (sha256 `99b3985b261dc3de561917736ea25c1e552283a67b2353df50f137a954bb0668`)
— now been executed to completion against the immutable R7 source tar, and
with what preregistered outcome?

Precondition named by `state.json.next_action`: "After GB9 releases its 26
threads". Verified BEFORE this prereg was written:
`physics/qldpc-dec/campaigns/20260902T114747Z_69eba724_61e9fe30ba5d/status.json`
reads `{"verdict": "CRASHED", "ended_utc": "2026-09-08T16:29:10Z"}`; no
python/stim/caffeinate compute processes are present on this workstation,
mini-0.local, or mini-pro.local (checked ~17:47Z). GB9's threads are
released; the launch precondition is satisfied.

## Init-command deviation (recorded before init, not after)

Brief step 3's literal gate slug
`after-gb9-releases-its-26-threads-follow-the-exa` is mechanically
incompatible with the pinned adapter `src/zero_level_repo_repro.py`:

1. `validate_run_dir` recomputes the run-id `hash12` as
   `sha256(GATE \x1f agent \x1f prereg_sha256 \x1f stamp \x1f uuid8)[:12]`
   with the adapter constant `GATE = "zero-level-author-repro-r7"`.
   `scripts/campaign.py:mint_run_id` uses the identical binder with the gate
   string passed to `init`. Any init gate other than
   `zero-level-author-repro-r7` mints a run id the adapter deterministically
   REFUSES ("run_id hash12 does not bind gate/agent/prereg/stamp/uuid8").
2. The adapter captures the prereg bytes of `physics/qlops/pre_statement.md`
   (hardcoded `PREREG_PATH`) and refuses any manifest whose `prereg_sha256`
   differs ("must be re-init'd against the current pre_statement.md").
3. The adapter refuses any prereg whose latest `## Revision N` heading is
   not exactly 7.

All three incompatibilities are provable by inspection; burning the cycle's
one funded run to demonstrate them is prohibited by budget (runs_max 1).
The cycle therefore inits with the exact two-command handoff frozen in
`physics/qlops/README.md` ("Official-source Revision 7 — frozen readiness
smoke only") and `state.json.next_action`:

```sh
python3 ../../scripts/campaign.py init \
  --gate zero-level-author-repro-r7 \
  --prereg pre_statement.md
```

This file remains the hmz freeze/cycle-prereg artifact; the machine-bound
prereg is `pre_statement.md` Revision 7. Nothing in this deviation widens
scope: the launched computation is byte-for-byte the preregistered one.

## Falsifiable acceptance criteria (inherited, frozen in Revision 5–7)

- 12 cells: ungrown shots [5e6,5e6,5e6,5e6,5e6,1e7] at
  p=[1e-3,8e-4,6e-4,4e-4,2e-4,1e-4]; grown [1e6,5e6,5e6,5e6,1e7,1e7];
  71,000,000 author shots total; fixed preregistered seeds; 10,000-shot
  chunks; one primary seeded run per point.
- Full mode re-runs the eight-gate smoke battery FIRST; any smoke-gate or
  shipped-vs-rebuilt `flattened()` inequality ⇒ REJECTED refusal before a
  single author-count shot.
- 24 comparisons (12 LER + 12 acceptance) against the author rows embedded
  in the pinned sources, combined analytic binomial sigmas, familywise 1%
  over 24 dependent comparisons, Bonferroni-derived |z| ≤ 3.53 gate.
- Preregistered verdict mapping (adapter `verdict_from_zs`, unchanged):
  any |z| ≥ 5 or ≥ 3 values beyond 3.53 ⇒ FROZEN-NEGATIVE;
  otherwise undefined statistics or 1–2 mild failures ⇒ FROZEN-INCONCLUSIVE;
  all 24 defined |z| ≤ 3.53 ⇒ FROZEN-CERTIFIED.

## Verdict mapping for `campaign.py close --verdict`

- Adapter completes full mode, candidate FROZEN-CERTIFIED ⇒ close
  **FROZEN-CERTIFIED** (certifies the executable d3/d9 artifact ONLY).
- Adapter completes full mode, candidate FROZEN-NEGATIVE ⇒ close
  **FROZEN-NEGATIVE**.
- Adapter completes full mode, candidate FROZEN-INCONCLUSIVE ⇒ close
  **FROZEN-INCONCLUSIVE**.
- Adapter refuses (exit 1, REJECTED) before sampling ⇒ close **REJECTED**.
- Run interrupted before `results/analysis.json` exists (e.g. cycle
  wall-clock budget 14,400 s reached) ⇒ stop at the last completed
  per-point artifact, freeze the partial state, close
  **FROZEN-INCONCLUSIVE**, and report the exact completed fraction
  (cells and shots) — partial progress is a result, never a certification.
- Adapter/process crash (no refusal, no analysis) ⇒ close **CRASHED**.

## Scope boundary

- IN: the one pinned 71m-shot/12-cell full run; its smoke re-gate; per-point
  artifacts; the adapter's own analysis; hash-integrity re-verification.
- OUT: paper d=7 and physical c≈300 (remain NOT-REPRODUCED regardless of
  outcome; no d=7 reconstruction is attempted); Gate A/B conclusions
  (unchanged); any modification of the R5/R6/R7 frozen run dirs; any change
  to `pre_statement.md`, the adapter, or the source tar; physics/qldpc-dec
  GB9 crash forensics (noted as precondition evidence only).
- Compute placement: the preregistered handoff pins the local
  `physics/.venv` (Python 3.14.3, stim 1.16.0, pymatching 2.4.0,
  numpy 2.5.2 — verified identical to the R7 smoke environment) and the
  "26 threads" precondition only exists on this 28-thread workstation;
  the full run therefore executes LOCALLY with OMP/OPENBLAS/MKL/NUMEXPR/
  VECLIB thread caps = 1 and `/usr/bin/nice -n 10`, launched detached
  (`nohup caffeinate -s ... &`). Control plane strictly local.

## Exact verification command

From `/Users/jinleic/jinleic-workspace/physics/qlops`:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
/usr/bin/nice -n 10 ../.venv/bin/python src/zero_level_repo_repro.py \
  --source-tar campaigns/20260904T054537Z_818bca48_18f5ac3b849c/source/source_tar.gz \
  --run-dir campaigns/<printed-run-id> \
  --mode full
```

Audit re-verification (fresh, post-run): `sha256sum -c sha256s.txt` inside
the frozen run dir; source tar sha256 must remain
`d67dbe7482b391984da5e64aeff7668bdaee45c262e6fb2dfd41dfed3d7f3338`.
