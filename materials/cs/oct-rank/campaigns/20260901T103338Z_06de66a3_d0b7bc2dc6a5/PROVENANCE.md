# PROVENANCE — N1 campaign (run 20260901T103338Z_06de66a3_d0b7bc2dc6a5)

- Minted by: `python3 scripts/campaign.py init --gate n1_exact_cp_completion
  --prereg n1_prereg.md` from `/Users/jinleic/jinleic-workspace/cs/oct-rank`,
  agent OctRankNext, 2026-09-01, after the path-scoped prereg commit.
- Prereg source commit: `d6c7e44d6901eb025835313fe8bf9c7b01ee3a81`
  (message: "oct-rank: preregister N1 exact-CP-completion campaign …").
- Prereg sha256 (workspace file, campaign.py-hashed, and this dir's
  pre_statement.md — all three identical):
  d569057eb8e26bf2f1a5eb9179de0676ada4b76dfd9f01473321de9e753f0f82
- Copy discipline (Main's order): campaign.py hashes the external prereg but
  does not copy it; `pre_statement.md` in this dir was copied byte-identically
  (`cmp` clean) from cs/oct-rank/n1_prereg.md AFTER init, BEFORE any
  control/main compute.
- Working tree note: cs/oct-rank/src/ is untracked repo state (pre-existing
  from earlier campaigns, zero modifications by this run — any needed N1
  code lives in THIS run dir and is frozen here).
- Pristine freeze: this dir must contain no __pycache__/*.pyc at freeze time;
  all scripts set PYTHONDONTWRITEBYTECODE=1 in-process before imports.
- Instrument heritage: rf_krawczyk.py (frozen,
  2026-09-01T04:30:00Z_routeF_kraw) is copied here as `n1_instrument_base.py`
  with an appended N1 driver (see n1_run.py) — the certification core
  (certify(), Sys, LADDER, controls semantics) is frozen code path-pinned in
  PROVENANCE via sha256 below; the N1-specific additions (seed classes,
  admission, exact Newton, merge seeds) are new code, frozen in n1_run.py.
