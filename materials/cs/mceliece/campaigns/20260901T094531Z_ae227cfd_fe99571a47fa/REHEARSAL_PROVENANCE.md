# REHEARSAL provenance — 20260901T094531Z_ae227cfd_fe99571a47fa (closed REHEARSAL, non-evidential)

## Why this run is a REHEARSAL

Main-directed (IRC 2026-09-01 ~10:0xZ, after owner steering on shared-worktree
git ordering): pre-commit t=96 tool probes make the run non-evidential under
the user's strict ordering (prereg commit BEFORE any t=96 compute), even
though no threshold adapted to any probe observation.

## Full timeline (all times UTC, 2026-09-01)

- 09:45:31 — `campaign.py init` minted this run id; prereg file
  `cs/mceliece/pre_statement_t96.md` existed on disk (uncommitted);
  claim file `cs/.claims/mceliece.claim.json` acquired.
- 09:46 — inherited instrument copied byte-identically into `code/`
  (verified vs current src and vs the frozen 14-09Z /
  M12ALPHA_DIVONLY snapshots; see INSTRUMENT_INHERITANCE.md — that
  record reads "at init, before compute", accurate for the code
  copies).
- 09:46:57 — anchors stage started (bg_45): the three REGISTERED
  m<=11 anchor rows (11,2048,48,6211), (10,1024,40,5113), (6,64,3,1387)
  against the frozen 2026-08-30T14-09Z_57200ADD records. This stage
  contains NO t=96 content (the rows and their expected outputs were
  frozen 2026-08-30, i.e. before this session). Result delivered ~09:54:
  byte_equal 3/3 TRUE, extras_ok TRUE, alpha MEASURED true on all
  three. NOT counted as campaign evidence (rehearsal); re-run
  post-init in the successor run.
- 09:48–09:50 — t=96 TOOL-ONLY probes (cost measurement, no instance,
  no verification): one rng draw + Rabin test on it (~0.05 s CPU);
  the observation that an empty polynomial trivially lies in every
  subspace (dimension argument for build-cost estimation). No
  (12,3488,96,16384) instance was built and no prereg threshold was
  adapted. Output DISCARDED — used only to size the (unchanged,
  pre-registered) budget; the pre-statement's budget numbers were
  derived from the FROZEN t=64 rates, and were NOT edited after the
  probes.
- 09:52 — main-agent steering on git ordering (path-scoped commits).
- 09:54:10 — target claim heartbeated (agent McelieceM12T96).
- 09:56:49 — prereg COMMITTED as `9189ff7` (path-scoped `git commit
  --only -- cs/mceliece/pre_statement_t96.md`; file byte-identical to
  what init consumed, sha256
  9d1e28d56c9e9af777c0bae062b375607e2a7857be599c4e37071f4e800d3047);
  byte-identical copy frozen into this run dir as `pre_statement.md`.
- 10:04 — a (12,3488,96,16384) `Instance(...)` build was STARTED as a
  cost probe (~503 s wall, TERM'd by Main's restart directive BEFORE
  completion; no artifact was written; nothing from it is used).
- ~10:1x — Main directive: restart cleanly. The run was made
  non-evidential by the pre-commit probes regardless of their
  innocence in content. Build probe terminated gracefully (SIGTERM,
  confirmed dead); no deletion of any file.
- freeze: 14 files hashed into sha256s.txt (no __pycache__, no .pyc —
  PYTHONDONTWRITEBYTECODE=1 throughout).
- close: verdict REHEARSAL, claim released.

## Discarded outputs statement

Nothing produced inside this rehearsal run is evidential for the t=96
campaign or reused by it: not the anchor files (they will be REDONE
inside the successor run post-init), not the instrument-inheritance
record (RE-written there from scratch), not any timing observation.
The successor campaign starts from the committed prereg `9189ff7`,
re-inits via `campaign.py init --gate ... --prereg pre_statement_t96.md`,
and runs anchors + verdict + plants fresh, serially, after init.

## Stage scripts (provisional, preserved but NOT authoritative)

`verdict_m12_t96_stage.py`, `anchors_stage.py`,
`plants_m12_t96_stage.py` in this dir were written pre-commit as
mechanical re-instantiations of the M12ALPHA_DIVONLY stages. The
successor run re-creates its own copies (fresh edits if needed), so no
stale-script hazard carries over. They are frozen here purely as the
run-dir state at freeze time.

— McelieceM12T96, 2026-09-01.
