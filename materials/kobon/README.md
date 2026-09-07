# Kobon triangle problem — exact decisions at small n

Canonical home of the Kobon campaign. Everything authoritative about this
target lives in this directory; `scratch/` holds only machine working data
(CNF instances, DRAT proofs, solver logs, venv, tools).

## Statement and convention

The Kobon triangle problem (Fujimura; on Wikipedia's list of unsolved problems
in mathematics) asks for the maximum number of non-overlapping triangles whose
sides lie on $n$ straight lines. This campaign decides the **broad
convention**, written $K_{\rm gen}(n)$: selected triangles have pairwise
disjoint interiors, other lines MAY cross a selected triangle, and parallel
lines and multiple intersection points are allowed.

**Campaign result:** $K_{\rm gen}(10) = 25$, **certified complete
2026-08-15** — machine-checkable certificates in both directions.

**Active front (opened 2026-08-15): decide $K_{\rm gen}(12)\in\{38,39\}$.**
UNSAT at target 39 is self-contained (any configuration with $\ge39$ disjoint
triangles contains a 39-subfamily, so UNSAT@39 ⇒ $K_{\rm gen}(12)\le38$ —
no reliance on the refuted published ceilings). Instances:
`scratch/kobon/n12/` — 15-cube cover (`simple`, `par_noc`, `par_conc`,
`c0`–`c11` from the constructively exhaustive D12 orbit enumeration),
dimensions: `simple` 65,594 vars / 1,937,254 clauses; `par_noc`
70,922 / 1,937,114; `par_conc` and `c0`–`c11` 60,194 /
1,925,750–1,925,815. `par_conc` is split 8-way on X-vars 545/927/1385
(byte-verified = parent + 3 units). A second 8-way split under
`par_conc_0` uses X-vars 295/775/1121: branches 1/3/5/7 are
DRAT-verified UNSAT.

A legacy compressed rerun produced `s VERIFIED` markers for four `par_noc`
branches, but those markers are **untrusted**: its wrapper passed the XZ
container itself to `drat-trim`, so input unit propagation could masquerade as
proof checking. The legacy wrapper now fails closed; these markers do not
change the bound.

The current replacement campaign is bound to schema-2 decision manifest
SHA-256 `48671f44f726dc4272438d0797beff936878c39b9b8d26ec110546b1136db7cf`.
Regeneration is byte-identical for the base instance and all 15 cubes; an exact
$D_{12}$ audit covers all 220 triples in 12 disjoint orbits. The proof-free
Kissat 4.0.4 discovery lanes `simple`, `par_noc`, `par_conc`, and `c0`–`c3`
are live; no cube has a verdict.

`n12_decide/certify_unsat.py` is the clean certificate path. It binds the
decision manifest, CNFs, tools, child environment, and its own source hash;
writes Kissat's binary DRAT directly inside XZ; streams decompression into
`drat-trim -i -w`; separates proof-consuming UNSAT from input-unit-propagation
UNSAT; preserves every attempt; and freshly rechecks terminal evidence on
resume. Solver production has an 8 GiB aggregate process-group RSS cap and a
64 GiB free-space reserve. The production status receipt
`n12_decide/certificate_status.final.json` records all 15 cases as
`UNATTEMPTED`, so this infrastructure is not a target-39 verdict. A real-tool
control in `n12_decide/certifier-real-smoke/` records
`PROOF_CONSUMING_BINARY_DRAT_XZ` with two nonzero core lemmas for a non-unit
UNSAT formula, `INPUT_UNIT_PROPAGATION_UNSAT` for the distinct input-UP
control, and `XZ_INTEGRITY_FAILED` for a one-byte-corrupted proof.

An independent all-degeneracy direct-gap lane now adds two proved accelerator
families: hereditary subarrangement bounds from the already established
$K(11)=32$ and $K(10)=25$, and local vertex-sector bounds (four faces per
simple line pair; two, or zero when nonadjacent, at a multipoint). The combined
$(12,39)$ instance has 374,381 variables / 830,030 clauses and is live under
both Kissat and CaDiCaL. Exact audit: 3,600 rational arrangements, 27,705
triangular faces and 93,000 line-pair checks, zero violations. The hash-bound
record `scratch/kobon/n12_gap_sector_deletion_t39.metadata.json` binds the
$K(10)=25$ input to the 2026-08-15 cube-cover certificate, not to the later
face monolith. This is a sound search accelerator, not a verdict; $K(12)$
remains open.

An additional **opt-in** shared-ray cut couples face adjacency to degeneracy:
two faces leaving their common vertex along the same ray of a shared line must
meet at the same remote endpoint, forcing a triple concurrency there. The
exact audit checked 1,591 such cases with zero violations and includes a
six-line four-sector sharpness control. Because solver performance was mixed,
the prior `--sector-bounds` encoding remains unchanged; the separate
`--shared-ray-bounds` lane has 374,381 variables / 841,910 clauses and is live
as `n12-gap-shared-ray-t39`. It has no verdict and does not change the open
status of \(K(12)\).

The stronger **endpoint-closure** lane now handles two faces that share only
one line but meet at a multipoint endpoint. Equal ray directions force their
other endpoints to coincide; a shared side forces the apexes to opposite
sides. A separate perturbation theorem says at most
\(\bar a_3^s(n)\) selected faces can avoid finite multipoints, so target 39 at
\(n=12\) forces at least two multipoint-incident faces. With four prime
\(K(4)=2\) clauses per four-line subset, the 11-/10-line deletion instance is
375,037 variables / 1,005,366 clauses (SHA-256
`6dca104aa120ee1f209aada8aca07a08bf3341f1335e6e99f4ae1f11bc2ba76b`);
adding the checked \(K(9)=21\) cuts gives 842,317 / 1,917,926 (SHA-256
`34406004c932ce77d33f76f8b7d19c7b0f5a2e3e241e94651846269bf750f2df`).
Both low-priority Kissat lanes are live without proof logging. Exact audits
found zero violations in 1,607 endpoint collisions, 84 new colliding same-ray
cases, 42 single-line shared segments, and 1,591 shared-pair segments. The
known \(n=8,9\) optima remain SAT; UNSAT controls show mixed performance.
These are sound discovery accelerators, not a verdict: \(K(12)\) remains open.

**Lower bound $K_{\rm gen}(12)\ge38$: PROVED** (exact rational) —
Kabanovitch's arrangement (Charade 1999; combinatorics via Savchuk's
LineOrder gallery table + SVG) reconstructed over $\mathbb{Q}$ with its two
triple points $\{0,3,5\},\{0,8,10\}$ imposed exactly;
`engine.verify_n12_lower_bound()` → `(True,'ok')` (38 triangles, all 703
pairs disjoint; no published exact coefficients existed — scout-confirmed).
Artifact: `scratch/kobon/n12/n12_lower_certificate.json`.

Why no located source proves this case (priority not exhaustively reviewed;
scouted 2026-08-13): the public tables' 25/38 entries at
$n=10,12$ trace to theorems about triangular *faces* of *simple* arrangements
(BBL arXiv:0706.0723, Blanc arXiv:0801.2845) or to the Clément–Bader draft
whose charging lemmas assume nondegenerate/pairwise-crossing configurations.
Exact rational counterexamples (in [`report.md`](report.md) §6) show those
charging steps are false once degeneracies and crossed triangles are admitted,
so no located source proves the broad-convention value at $n=10$. Scope
verdicts were independently re-derived by two literature scouts (2026-08-13).

## State (certified 2026-08-15; independently recertified 2026-08-23)

| claim | status | evidence |
|---|---|---|
| $K_{\rm gen}(10)\ge25$ | **PROVED** (exact rational) | integer certificate in `engine.py::N10_LOWER_BOUND_LINES/TRIANGLES`; `verify_n10_lower_bound()` → `(True,'ok')`; independently re-implemented by a reviewer agent (300/300 pairs separated) |
| CNF relaxation sound (UNSAT ⇒ upper bound) | **PROVED** (necessity of every clause) | hand proofs + 1,296 planted-degeneracy configs, 0 violations; 203,673 triangle pairs vs exact separating-axis oracle, 0 gaps; adversarial reviewer audit found no unsound clause |
| 11-cube cover exhaustive | **PROVED** | `dump_instances` docstring + CubeCoverAudit (D10 orbit transversal, Burnside $(120+5\cdot8)/20=8$; single WLOG σ-bit via $b\mapsto-b$) |
| $K(9)\le21$ (control) | **CERTIFIED** | `proof_n9_t22.drat` `s VERIFIED` (drat-trim, 121,436,534 resolution steps) |
| 9 of 11 cubes UNSAT at target 26 | **CERTIFIED** | all 9 shipped DRATs `s VERIFIED` by drat-trim (≈23.4 CPU-h; `scratch/kobon-audit/results_verified.tsv`) |
| `c0`, `par_conc` cubes UNSAT | **CERTIFIED** (verified 8-way case splits) | 16/16 subcube DRATs `s VERIFIED` (≈16.4 CPU-h, 44 GB; `scratch/kobon-audit/splits.tsv`); each subcube CNF byte-verified = parent + 3 unit clauses (tautological split, `AUDIT.md`) |
| Independent faces-only refutation at target 26 | **CERTIFIED** | `ladder_n10_t26_faces.drat` checked `s VERIFIED`: 131,314 variables / 792,462 clauses; 52,808 formula clauses and 40,395,481 lemmas in core; record `scratch/kobon/ladder_n10_t26_faces_certificate.json` |
| $K_{\rm gen}(10)=25$ | **THEOREM — fully certified** | every row above green; verdict sealed in [`AUDIT.md`](AUDIT.md) |

## Canonical release

The canonical release archive is **`math/kobon/release/kobon-2026-08.zip`**
(19,286,189 bytes, SHA-256
`51e3789d1bb7b571302e6ab990f5a31d0b21790a316484cb83fdd2a984fca70d`, pinned by
the sidecar `math/kobon/release/kobon-2026-08.zip.sha256`), together with its
unpacked mirror `math/kobon/release/kobon-2026-08/` (133-file
`MANIFEST.sha256`, all 133 verified). This is a ruling, not a convention:
`math/kobon/docs/INDEX.md` §"Canonical release snapshot" and
`math/kobon/docs/ARTIFACTS.md` both name it canonical, and
`math/PROGRESS.md` records its packaging (2026-08-21, 44/44 manifest
verified). All earlier packaging zips — `kobon-2026-08-v1` through
`kobon-2026-08-v8` and the three `kobon-2026-08-obsolete-structure*`
archives — are superseded packaging history, not alternate authorities, and
live in `math/kobon/release/superseded/` (see that directory's
`MANIFEST.md` for hashes, provenance, and per-file restore commands).
No zips are deleted; nothing here changes any result.

## Layout

| path | contents |
|---|---|
| [`README.md`](README.md) | this file — SSOT index and status |
| [`report.md`](report.md) | the theorem report: convention, encoding, machine-verified lemmas L1–L5, axioms A1–A7, rules R1–R4, accelerators, counterexamples to classical charging, results |
| [`engine.py`](engine.py) | canonical engine: CNF builder, exact rational oracles, LP straightening, embedded lower-bound certificates (`verify_n10_lower_bound`, `verify_n12_lower_bound`), `dump_instances` proof-family driver. sha256 now `0b0dfa2a…e990` (n=12 certificate added 2026-08-15; the audited n=10 bundle hash `1f524967…79fe` is preserved in `AUDIT.md`) |
| [`AUDIT.md`](AUDIT.md) | 2026-08-15 independent certification audit: truncation finding, byte-census method, reproducibility seal, in-flight verification fleet, closure criteria |
| `~/jinleic-workspace/scratch/kobon/` | working data: 12 proof-family CNFs (regenerate byte-identically from `engine.py`), 9 complete DRATs + 2 truncated, historical instances |
| `~/jinleic-workspace/scratch/kobon-audit/` | live audit working dir: tools (kissat 4.0.4, drat-trim), split subcube CNFs, replacement proofs (`resolve/`), fleet logs, `results.tsv`, `memwatch.log` |
| `~/Downloads/files-to-transfer-math-phys/workspace-files/scratch/kobon/` | frozen source bundle audited by AUDIT.md (byte-identical `engine.py`) |

## Verification (rerun anytime)

```sh
cd ~/jinleic-workspace
# lower bound, exact rational:
scratch/kobon-audit/venv/bin/python -c "import sys; sys.path.insert(0,'math/kobon'); import engine; print(engine.verify_n10_lower_bound())"
# regenerate the proof-family CNFs and diff against shipped:
scratch/kobon-audit/venv/bin/python -c "import sys; sys.path.insert(0,'math/kobon'); import engine; print(engine.dump_instances(10,26,'/tmp/kobon-regen'))"
# check any cube proof:
scratch/kobon-audit/tools/drat-trim/drat-trim scratch/kobon/kobon_n10_t26_proof_simple.cnf scratch/kobon/proof_n10_t26_simple.drat -w
```

## Process history (all retired 2026-08-15)

`kobon-drat-fleet` (9-proof recheck — complete, 9/9 `s VERIFIED`),
`ksplit-*` / `vsplit-*` (16 split solves + 16 verifications — complete,
16/16 `s VERIFIED`), `kobon-solve-c0` / `kobon-solve-parconc` (whole-cube
backup re-solves — stopped as superseded once their cubes were certified),
`kobon-memwatch` (RSS observability — retired).

## Discipline notes

- UNSAT results are believed only after `drat-trim` verification of a complete
  DRAT ending in the empty clause; solver exit codes are not evidence.
- SAT results would be believed only after LP straightening + exact rational
  re-verification (none occurred at target 26).
- Every WLOG (the single σ-bit, cube covers) carries an explicit group-action
  proof; a false x-flip justification and an unsound double-σ cube family were
  caught by adversarial audit and killed before use (report §5).
