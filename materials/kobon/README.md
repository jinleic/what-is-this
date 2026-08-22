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
DRAT-verified UNSAT. A compressed, completion-gated rerun now has
`par_noc` branches 1/3/5/7 independently `s VERIFIED`; the remaining
branches of `par_noc`, `simple`, `par_conc_0`, `c0`, and `c1` are active.
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

## State (2026-08-15, final)

| claim | status | evidence |
|---|---|---|
| $K_{\rm gen}(10)\ge25$ | **PROVED** (exact rational) | integer certificate in `engine.py::N10_LOWER_BOUND_LINES/TRIANGLES`; `verify_n10_lower_bound()` → `(True,'ok')`; independently re-implemented by a reviewer agent (300/300 pairs separated) |
| CNF relaxation sound (UNSAT ⇒ upper bound) | **PROVED** (necessity of every clause) | hand proofs + 1,296 planted-degeneracy configs, 0 violations; 203,673 triangle pairs vs exact separating-axis oracle, 0 gaps; adversarial reviewer audit found no unsound clause |
| 11-cube cover exhaustive | **PROVED** | `dump_instances` docstring + CubeCoverAudit (D10 orbit transversal, Burnside $(120+5\cdot8)/20=8$; single WLOG σ-bit via $b\mapsto-b$) |
| $K(9)\le21$ (control) | **CERTIFIED** | `proof_n9_t22.drat` `s VERIFIED` (drat-trim, 121,436,534 resolution steps) |
| 9 of 11 cubes UNSAT at target 26 | **CERTIFIED** | all 9 shipped DRATs `s VERIFIED` by drat-trim (≈23.4 CPU-h; `scratch/kobon-audit/results_verified.tsv`) |
| `c0`, `par_conc` cubes UNSAT | **CERTIFIED** (verified 8-way case splits) | 16/16 subcube DRATs `s VERIFIED` (≈16.4 CPU-h, 44 GB; `scratch/kobon-audit/splits.tsv`); each subcube CNF byte-verified = parent + 3 unit clauses (tautological split, `AUDIT.md`) |
| $K_{\rm gen}(10)=25$ | **THEOREM — fully certified** | every row above green; verdict sealed in [`AUDIT.md`](AUDIT.md) |

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
