# Physics / open-problem attack repository

**Start here:** [`RESULTS.md`](RESULTS.md) is the authoritative current-results
index. This README is the repository map. Provisional claims and candidate
results live in a non-authoritative `/scratch/` subfolder until verified.

## Mission

Attack open problems in physics — prioritizing **quantum computing** (QEC, qLDPC
decoding, fault tolerance, quantum complexity) — with the same verification
discipline as the sibling `../math/` repository.

## Tracking contract (SSOT)

One fact, one owner. Four layers, no duplication:

| layer | owner | contents |
|---|---|---|
| authoritative results index | [`RESULTS.md`](RESULTS.md) | current headline result, status per problem, campaign ladder, reproduction, and adopted next targets |
| chronological ledger | [`PROGRESS.md`](PROGRESS.md) | dated, newest-first session entries; every claim carries its verification; retractions inline |
| current state per target | `<target>/README.md` | what is PROVED / CONDITIONAL / NUMERICAL right now, file inventory, how to run |
| run artifacts | `<target>/campaigns/<UTC-timestamp>_<uuid>_<code-hash12>/` | immutable frozen snapshots + committed results; inventoried in `<target>/campaigns/README.md` |
| generated compact state, per target | `<target>/state.json` | owner-maintained semantic fields (`problem`, `current_gate`, `headline`, `next_action`); mechanics refreshed by `scripts/campaign.py state refresh`; institute-managed fields (`status`, `lifecycle`, `budget`, `priority`) owned by `scripts/institute.py` — see `docs/knowledge-system/INSTITUTE.md` |

Generated, not narrative: when the one-line summary in `<target>/state.json` and an
owner document above disagree, the owner document wins — fix [`RESULTS.md`](RESULTS.md),
`<target>/README.md`, or [`PROGRESS.md`](PROGRESS.md) first, then re-run
`scripts/campaign.py state refresh` so the generated copy follows.

Naming rules:

* **One problem, one top-level folder**, short lowercase name (`qldpc-dec/`,
  `sdec/`, `msd/`, `qmaqcma/`). No dated top-level folders, no second folder
  for the same problem.
* Campaign directories are producer-generated
  (`cert3_<UTC>Z_<uuid>_<hash>`); never hand-created, never edited after
  launch. Disposable rehearsals go under `<target>/campaigns-smoke/`.
* Update the ledger and the owning README **in the same session** as the work;
  a result that is not in `PROGRESS.md` does not exist.
* Unverified/candidate material lives under `<target>/scratch/` (non-authoritative);
  promotion to a README requires the evidence labels below.

Original targets opened 2026-08-29 (**six, one owner session each**): `qldpc-dec/`,
`msd/`, `shadows/`, `na-compiler/`, `fss-bb/`, `qlops/`. Per-target state
lives in each `README.md`; RESULTS.md collects verdicts.

Shortlists: v1 in [`docs/QC_FRONTIER_2026-08-29.md`](docs/QC_FRONTIER_2026-08-29.md);
the **five-scout deep scan** in
[`docs/QC_SCOUT_SCAN_2026-08-29.md`](docs/QC_SCOUT_SCAN_2026-08-29.md) is
the selection basis for the six targets above.

Opened-target gate table (2026-08-29, five-scan ranking):

| rank | candidate folder / problem | first falsifiable gate |
|---:|---|---|
| 1 | `qldpc-dec/` — cross-paper decoder reproduction harness for BB codes (BP+OSD baseline replay, beam-search and GARI headline checks) | reproduce published [[144,12,12]] curve within 95% bootstrap CI; ≥2σ miss on a beam-search/GARI claim = refutation |
| 2 | `msd/` — magic-state supply reproduction + space-time scoreboard (zero-level CCZ, sqrt(T) catalysts, cultivation attribution) | re-run 22-qubit zero-level CCZ in stim; fit p_L=c·p², reproduce c≈300 within 2× |
| 3 | `shadows/` — contractive-shadow exact census (11,520 two-qubit Cliffords; 1.8^k base search) | exhaustive Clifford enumeration confirming/refuting the unpublished ≤4-of-9 lemma |
| 4 | `na-compiler/` — certified-optimal small-instance transport baselines (exact ILP vs MQT QMAP/Enola) | certified-optimum gap table on ≤10-qubit instances under one duration model |
| 5 | `fss-bb/` — BB finite-size scaling + ν anomaly | independent erasure-FSS re-run; ν drift to 4/3 vs ≈1.18 decision |
| 6 | `qlops/` — FTQC resource-estimate arithmetic audit (QLOPS sensitivity) | reproduce published RSA-2048 numbers from stated formulas; >2× shifts under decoder-latency sweep |

**Gate revisions after first-hand reads (authoritative text = each target's
`pre_statement.md`):** row 1's "published LER curve" does not exist in any of
the four decoder papers — `qldpc-dec/` Gate A was re-pinned to the published
*runtime* quantities (arXiv:2512.07057 Table II/III). Row 2's `c≈300` test is
only meaningful inside the quadratic window, so `msd/` added Gate-A-prime at
p ∈ {1e-4, 3e-4, 1e-3} (pre_statement rev1/rev2). Verdicts live in
[`RESULTS.md`](RESULTS.md), never in this table.

Superseded v1 rows: `sdec/` folded into `qldpc-dec/`. Parked (not
re-examined): `qmaqcma/`, `scars/`, `mbl2d/` from v1.

Opened 2026-09-06: [`stabrank-witness/`](stabrank-witness/README.md), an exact
magic-cat stabilizer-rank witness target. Its first gate calibrates the public
classical search before probing a five-term cat-eight construction. Current
state belongs to its README and RESULTS.md, not the historical ranking above.

Note: `../math/qec/` is already active on the **algebraic** qLDPC side
(non-CSS PBB programme, distance-reversal censuses). Any new `physics/qldpc-dec/`
target must therefore be scoped to **receiver-side** questions (decoders,
noise models, circuit-level thresholds) to stay disjoint from `math/qec/`.

## Discipline

Inherited from `../math/README.md` without modification:

1. Read primaries first-hand; quote verbatim. Secondary summaries are not evidence.
2. Every quantitative claim is machine-checked by a script in this repo.
3. Distinguish PROVED / CONDITIONAL / NUMERICAL at every step, always.
4. Verify by hand anything an optimiser reports before believing it.
5. Record retractions inline rather than silently editing claims away.
6. A theorem is not proved because its asserts pass — asserts check arithmetic at
   the points you chose.
7. **Never promote a sampled parameter check to a universal claim.** Where the
   quantity is monotone or affine in the parameter, solve for the boundary and
   assert the equivalence. Sampling a range is evidence, never a quantifier.

Physics-specific additions:

8. **Every simulation number labels its noise model.** A logical-error curve is
   meaningless without (a) circuit vs code noise, (b) idle/leakage/crosstalk
   assumptions, (c) a Monte-Carlo seed and confidence-interval recipe.
9. **Publication-figure reproduction ≠ re-derivation.** Tag a reproduced curve
   `[REPRODUCED]` (same inputs, same pipeline assumptions) versus `[DERIVED]`
   (independent implementation). They are different claims.
10. **Hardware-adjacent numbers always cite the platform and date** (qubit
    count, basis gates, error rates) — a threshold is a property of a model.
11. Complexity-theory claims follow `math/` rules: a reported result quoted
    without a primary-source read is `[REPORTED]`.

## Resource policy

- One low-priority process (`nice -n 10`) at a time.
- Pin OMP, OpenBLAS, MKL, vecLib, and NumExpr to one thread.
- Bounded searches with explicit stop conditions; no daemonized CPU loops.
- Monte-Carlo campaigns: fixed seed per run, seed recorded in the campaign
  inventory; repeated-seed policy defined before the first run.
- Run prior-art review before announcing novelty. Keep failed routes and
  corrections in the chronological ledger.

## Layout

| Path | Contents |
|---|---|
| `RESULTS.md` | Authoritative per-target results index (skeleton until first result). |
| `PROGRESS.md` | Chronological session ledger (newest first). |
| `<target>/` | One folder per problem, when open. |
| `<target>/README.md` | Per-target authoritative state. |
| `<target>/campaigns/` | Frozen run snapshots, producer-generated names. |
| `<target>/scratch/` | Non-authoritative candidate material. |
| `docs/` | Frontier scans, benchmarks, and background reading. |
| `gate_a_results.json`, `gate_b_results.json` | **NON-AUTHORITATIVE STRAY COPIES — do not cite.** See the note below. |

**Repo-root artifact duplication, unresolved 2026-08-30.** `gate_a_results.json` and
`gate_b_results.json` sit at the repository root, where the layout above defines no
artifact location. They are byte-identical to the frozen `qlops/` campaign artifacts —
md5 `5fef096d7b2f345ef39aeff311fb5f62` and `e9e0b38306cd933583ed8012453e5aba`
respectively, matching
`qlops/campaigns/20260829T223540Z_dcf693ab/artifacts/` — so they carry no information
that the frozen snapshot does not already own, and they constitute a second copy of a
fact whose sole owner is that campaign directory. The authoritative copies are the
campaign ones; **these two files are non-authoritative and must not be cited, read as
input, or updated.** Deletion is the correct resolution and is **pending explicit user
approval**, since this repository requires confirmation before removing any file; until
that approval arrives the duplication is recorded here rather than silently tolerated.

## Evidence discipline

1. Every physical/computational claim carries a source: arXiv ID, DOI, or an
   in-repo simulation with committed campaign artifacts.
2. Numbers copied from papers name the figure/table and the arXiv/DOI reference.
3. Not personally verified against a primary source → tag `[REPORTED]`.
4. Inferred rather than read → tag `[INFERENCE]`.
5. Conflicting values are recorded as conflicts, not silently resolved.
6. Numerical candidates remain `NUMERICAL`; only exact/verifiable evidence
   (Lean proof, DRAT certificate, exact rational arithmetic, machine-checked
   replay of a full campaign) can be promoted to `PROVED`.
7. **A no-go theorem is worth exactly its hypotheses.** Transcribe them; the
   surviving hypotheses are where a construction can live.
