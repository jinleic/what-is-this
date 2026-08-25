# Kobon campaign documentation index

This is the navigation root for the consolidated 2026-08 Kobon campaign.
Narrative documents live under `math/`; large working objects remain under
`scratch/` and are inventoried in `ARTIFACTS.md`.  Paths in backticks are
workspace-relative.

## Evidence classes

| Class | Meaning |
|---|---|
| **PROVED** | A mathematical argument is written with its hypotheses and scope; no solver verdict is being substituted for the proof. |
| **MACHINE-CHECKED** | A finite/exact check or independently checked certificate completed and its pass condition is identified. |
| **DISCOVERY-ONLY** | Search, heuristic, SAT model, solver-only UNSAT, novelty search, or diagnostic evidence that is not a proof. |
| **PENDING** | The relevant run/check/argument has no settled verdict, is conditional, queued, or in flight. |

A document may carry more than one class because it contains separately scoped
claims.  The least-settled claim does not invalidate a proved section, and a
proved analytic branch does not settle an unrelated live SAT run.

## Start here

| Document | One-line description | Evidence class |
|---|---|---|
| `math/kobon/docs/INDEX.md` | This navigation map, status vocabulary, and disk-drift warnings. | MACHINE-CHECKED inventory |
| `math/kobon/paper_kobon_2026-08.md` | Consolidated paper: capacity theorem, branch obstructions, verified 54 lower bound, machinery, literature, and limitations. | PROVED + MACHINE-CHECKED + PENDING |
| `math/kobon/paper/kobon_broad_capacity.{tex,pdf}` | Current 17-page research paper: convention collapse, capacity theorem, certified values, square-penalty refutation, sector identity, and checked coefficient-2/3 endpoint at n=7. | PROVED + MACHINE-CHECKED + PENDING frontiers |
| `math/kobon/docs/methods.md` | Exact description of the SAT encoding, rational geometry, cube/monolith workflow, DRAT policy, and fixed-arrangement packing prober. | PROVED + MACHINE-CHECKED |
| `math/kobon/docs/experiments.md` | Snapshot ledger of completed, interrupted, queued, and live experiments with actual DIMACS headers. | MACHINE-CHECKED + DISCOVERY-ONLY + PENDING |
| `math/kobon/docs/reproduce.md` | Exact workspace commands for the lower certificate, 15-solution sweep, obstruction checks, cube builds, max packing, and engine test. | MACHINE-CHECKED procedures |
| `math/kobon/docs/ARTIFACTS.md` | Size/hash/release-inclusion inventory for every computational artifact referenced by this hub. | MACHINE-CHECKED inventory |
| `math/kobon/verification/README.md` | Script-by-script purpose, exit-code contract, expected output, and relocation caveats. | MACHINE-CHECKED documentation |

## Canonical documents under `math/`

| Document | One-line description | Evidence class |
|---|---|---|
| `math/kobon/README.md` | Original canonical status page and rerun instructions for the certified `K_gen(10)=25` campaign; its state table is dated 2026-08-15. | MACHINE-CHECKED |
| `math/kobon/AUDIT.md` | Byte-level audit and completed recovery ledger for the n=10 DRAT cube family. | MACHINE-CHECKED |
| `math/kobon/report.md` | Full computational-proof report for `K_gen(10)=25`, including encoding, validation, capacity theorem, counterexamples, and artifact map. | PROVED + MACHINE-CHECKED |
| `math/kobon/paper_capacity.md` | Standalone proof of the crossing-refined capacity theorem, equality analysis, sharpness checks, literature comparison, and open-target consequences. | PROVED + MACHINE-CHECKED |
| `math/kobon/paper/kobon_broad_capacity.{tex,pdf}` | Current full paper and compiled PDF, including the 2026-08-23 square-penalty frontier. | PROVED + MACHINE-CHECKED + PENDING frontiers |
| `math/kobon/paper_kobon_2026-08.md` | Current short synthesis and scope ledger for the late-August campaign. | PROVED + MACHINE-CHECKED + PENDING |
| `math/PROGRESS.md` | Multi-campaign chronological ledger; Kobon entries begin at the 2026-08-15 canonicalization and continue through the 2026-08-21 release/census. | MACHINE-CHECKED + DISCOVERY-ONLY + PENDING |

## Verification material under `math/kobon/verification/`

| Document / script | One-line description | Evidence class |
|---|---|---|
| `README.md` | Exit and output contracts for all copied verifiers. | MACHINE-CHECKED documentation |
| `maiorana_verify.py` + `.log` | Fail-closed solution-1 checker and its successful seven-stage JSONL transcript. | MACHINE-CHECKED |
| `sweep_all_verify.py` + `.log` | All-15 hash/census sweep and its successful 15/15 transcript. | MACHINE-CHECKED |
| `appendix_independent_check.py` | Pure-stdlib exact-integer checker for both n=14 Q=3 obstruction legs. | MACHINE-CHECKED finite leg |
| `appendix_check_n18.py` | Pure-stdlib exact-integer checker for n=18 placements, orbits, and reflection counts. | MACHINE-CHECKED finite leg |
| `max_packing.py` | Exact overlap-graph/SAT optimum prober for one fixed rational arrangement. | MACHINE-CHECKED when its result says `optimality_proved=true` |

## Canonical release snapshot

The canonical archive is `math/kobon/release/kobon-2026-08.zip`; the unpacked
mirror is `math/kobon/release/kobon-2026-08/`.  Earlier `*-obsolete-*` and `-v1`
archives are superseded packaging history, not alternate authorities.

| Release document | One-line description | Evidence class |
|---|---|---|
| `release/kobon-2026-08/README.md` | Package contents, result-to-artifact map, omissions, core verification commands, and provenance. | MACHINE-CHECKED package map + PENDING live-object notes |
| `release/kobon-2026-08/MANIFEST.sha256` | Integrity manifest for every file in the canonical archive snapshot. | MACHINE-CHECKED |
| `release/kobon-2026-08/papers/paper_kobon_2026-08.md` | Packaged snapshot of the consolidated paper. | PROVED + MACHINE-CHECKED + PENDING |
| `release/kobon-2026-08/papers/paper_capacity.md` | Packaged snapshot of the full capacity paper. | PROVED + MACHINE-CHECKED |
| `release/kobon-2026-08/papers/kobon_broad_capacity.{tex,pdf}` | Packaged 17-page current paper. | PROVED + MACHINE-CHECKED + PENDING frontiers |
| `release/kobon-2026-08/papers/report.md` | Packaged snapshot of the n=10 report. | PROVED + MACHINE-CHECKED |
| `release/kobon-2026-08/certificates/c23_n7_*` | Hash-bound CNF, 34 MB DRAT, `s VERIFIED` transcript and manifest for the coefficient-2/3 endpoint at n=7. | MACHINE-CHECKED theorem |
| `release/kobon-2026-08/verification/square_penalty_counterexample.json` | Exact A(12,1) face/edge/sector census refuting the coefficient-1 square penalty. | MACHINE-CHECKED counterexample |
| `release/kobon-2026-08/proofs/Kgen14_ge54_certificate.md` | Provenance and three-way exact verification matrix for the 54 lower certificate. | MACHINE-CHECKED |
| `release/kobon-2026-08/proofs/escape_obstruction.md` | Packaged n=14 Q=3 analytic obstruction and hypothesis ledger. | PROVED + MACHINE-CHECKED finite leg |
| `release/kobon-2026-08/proofs/obstruction_n18.md` | Packaged n=18 Q=3 analytic obstruction and finite-check interface. | PROVED + MACHINE-CHECKED finite leg |
| `release/kobon-2026-08/proofs/obstruction_n20.md` | Hand analysis showing why n=20 target 117 is not saturated; residual checks remain. | PROVED conditional statements + PENDING branch |
| `release/kobon-2026-08/proofs/defect_lemma.md` | Exact cycle-level escape families demonstrating that the current n=20 machinery cannot close the defect branch. | PROVED method limitation + PENDING branch |
| `release/kobon-2026-08/docs/litrefs.md` | Source-audited literature quotes, corrected identifiers, and citation records. | MACHINE-CHECKED source audit |
| `release/kobon-2026-08/docs/orbit_cube_plan.md` | Planned n=14 concurrency-orbit queue; commands in it are explicitly not executed builds. | PENDING |
| `release/kobon-2026-08/docs/batch2_launch.md` | Fingerprinted n=18 CNFs and the exact q2 discovery-launch ledger. | MACHINE-CHECKED fingerprints + PENDING runs |
| `release/kobon-2026-08/docs/signature_algebra.md` | Hand-derived n=20 degeneracy lattice and cube-readiness plan with an explicit unrun verification queue. | PENDING |
| `release/kobon-2026-08/docs/PROGRESS.md` | Packaged snapshot of the workspace progress ledger. | MACHINE-CHECKED + DISCOVERY-ONLY + PENDING |

## Campaign-authored working documents under `scratch/kobon/`

These remain references rather than canonical publication copies.  Their status
words are preserved; newer canonical documents take precedence when they
conflict.

| Working document | One-line description | Evidence class |
|---|---|---|
| `scratch/kobon/litrefs.md` | Live literature pin file with verbatim quotations, URLs, attribution corrections, and current OEIS scope. | MACHINE-CHECKED source audit |
| `scratch/kobon/novelty_audit_capacity_theorem.md` | Source-by-source novelty search for the capacity inequality; explicitly not a correctness audit. | DISCOVERY-ONLY novelty evidence |
| `scratch/kobon/n14_face_status.md` | Separates broad straight-line, simple pseudoline, and face conventions at n=14 and audits Blanc/Savchuk scope. | PROVED cited results + PENDING broad case |
| `scratch/kobon/n11_chain.md` | Audits the n=11 target-33 narrowing, monolith scope, lower construction, and missing upper certificate. | MACHINE-CHECKED lower bound + PENDING upper bound |
| `scratch/kobon/n11/section_draft.md` | Contingent paper section for `K_gen(11)=32`; final DRAT checkbox remains open. | PENDING |
| `scratch/kobon/n11/equality_branches.md` | Shows exactly why the Theorem-M cycle mechanism is silent at n=11 and records exact local checks. | PROVED limitation + MACHINE-CHECKED microchecks |
| `scratch/kobon/n12/closeout_draft.md` | Arbitration draft for the proved n=12 lower bound and still-pending target-39 upper certificate. | MACHINE-CHECKED lower bound + PENDING upper bound |
| `scratch/kobon/n12/path_arbiter_rec.md` | CPU-free comparison of monolith/cube strategies and failure modes at n=12 target 39. | DISCOVERY-ONLY + PENDING |
| `scratch/kobon/n14/Kgen14_ge54_certificate.md` | Source copy of the pinned Maiorana provenance and independent-verification matrix. | MACHINE-CHECKED |
| `scratch/kobon/n14/escape_obstruction.md` | Full proof of the n=14 Q=3 no-concurrency obstruction, validation ledger, and reproduction map. | PROVED + MACHINE-CHECKED finite leg |
| `scratch/kobon/n14/concurrency_cubes.md` | Budget enumeration, cube coverage design, built-run ledger, and encoding API status. | MACHINE-CHECKED enumeration + DISCOVERY-ONLY + PENDING |
| `scratch/kobon/n14/orbit_cube_plan.md` | Next-orbit queue and exact budget checks; states that listed future builds were not performed. | PENDING |
| `scratch/kobon/n14/work/NOTES.md` | Savchuk-tooling probe results, generator bugs, exact stack controls, and strategic caveats. | DISCOVERY-ONLY + MACHINE-CHECKED controls |
| `scratch/kobon/n14/work/line-order/README.md` | Upstream straightening tool documentation used by the probe campaign. | DISCOVERY-ONLY tooling reference |
| `scratch/kobon/n14/work/line-order/CHANGELOG.md` | Upstream LineOrder release history. | DISCOVERY-ONLY tooling reference |
| `scratch/kobon/n14/work/line-order/TODO.md` | Upstream LineOrder limitations and open work. | PENDING tooling reference |
| `scratch/kobon/n14/work/kobon-cnf/README.md` | Upstream Savchuk CNF generator documentation. | DISCOVERY-ONLY tooling reference |
| `scratch/kobon/n14/work/kobon-cnf/CHANGELOG.md` | Upstream KobonCNF release history. | DISCOVERY-ONLY tooling reference |
| `scratch/kobon/n14/work/kobon-cnf/TODO.md` | Upstream KobonCNF open work. | PENDING tooling reference |
| `scratch/kobon/n18/obstruction_n18.md` | Full n=18 Q=3 analytic obstruction with an executed finite-check appendix. | PROVED + MACHINE-CHECKED finite leg |
| `scratch/kobon/n18/batch2_launch.md` | Actual n=18 CNF headers, sizes, hashes, process check, and q2 launch command. | MACHINE-CHECKED fingerprints + PENDING runs |
| `scratch/kobon/n20/signature_algebra.md` | Symbolic degeneracy-signature algebra and an intentionally unrun verification item. | PENDING |
| `scratch/kobon/n20/obstruction_n20.md` | Corrected hand derivation separating saturated consequences from the real three-unit defect slack. | PROVED conditional statements + PENDING branch |
| `scratch/kobon/n20/defect_lemma.md` | Exact affine cycle countermodels and the two residual statements needed for any n=20 closure. | PROVED method limitation + PENDING branch |

Third-party license/readme files inside
`scratch/kobon/n14/venv_tbl/site-packages/` are dependency documentation, not
campaign documents; they are intentionally excluded from this navigation map.

## Working artifact map

| Artifact family | Purpose | Evidence class | Detailed inventory |
|---|---|---|---|
| `scratch/kobon/n14/maiorana/` | Fifteen pinned rational 14-line witnesses, upstream hash manifest, and 15 proved max-packing outputs. | MACHINE-CHECKED | `ARTIFACTS.md`, Maiorana tables |
| `scratch/kobon/discoveries/` | Persisted rational witnesses for n=9,11,14,18,20 and exact max-packing outputs. | MACHINE-CHECKED where verifier fields pass; otherwise DISCOVERY-ONLY | `ARTIFACTS.md`, discovery tables |
| `scratch/kobon/n14/escape/` | Exact scripts and JSON diagnostics supporting Theorem M. | PROVED + MACHINE-CHECKED + DISCOVERY-ONLY support | `ARTIFACTS.md` |
| `scratch/kobon/n14/*.cnf`, root n=14 cubes/logs, `q2pos/`, `q3pos/` | Target-54/55/56 monoliths and structural/position cubes. | PENDING / DISCOVERY-ONLY until proof-checked | `experiments.md`; `ARTIFACTS.md` |
| `scratch/kobon/n18/*.cnf` and logs | Target-94 q0/q1/q2 discovery runs plus held q3 cases. | PENDING / DISCOVERY-ONLY | `experiments.md`; `ARTIFACTS.md` |
| `scratch/kobon/n20/*.cnf` | Target-117 q0--q4 cube batch. | PENDING | `experiments.md`; `ARTIFACTS.md` |
| `scratch/kobon/n11_monolith_t33.{cnf,drat}` | In-flight target-33 proof attempt; DRAT is not complete. | PENDING | `experiments.md`; `ARTIFACTS.md` |
| `scratch/kobon/max_packing.py` | Fixed-arrangement exact optimum prober. | MACHINE-CHECKED only with `optimality_proved=true` | `methods.md`; `reproduce.md` |
| `scratch/kobon/n14/concurrency_cases.py` | Target-controlled n=14 signature enumeration and cube/monolith builder. | MACHINE-CHECKED generator; generated run verdicts separate | `methods.md`; `reproduce.md` |
| `scratch/kobon/square_penalty_counterexample.py` + JSON | Exact rational A(12,1) construction, dual face predicates, edge census and sector-run excess replay. | MACHINE-CHECKED counterexample | `paper/kobon_broad_capacity.tex`; `reproduce.md` |
| `scratch/kobon/square_penalty_sat.py` + `c23_n7_*` | Rational-weight direct-gap violation generator and checked coefficient-2/3 n=7 certificate. | MACHINE-CHECKED finite theorem; general endpoint PENDING | `paper/kobon_broad_capacity.tex`; `reproduce.md` |
| `scratch/kobon/sector_bound_audit.py` + JSON | Exact-rational replay of simple/multipoint sectors, shared-pair rays, single-line endpoint closure, and shared-edge opposite-side parity over six degeneracy modes, including a four-sector sharpness control. | PROVED local lemmas + MACHINE-CHECKED implementation audit | `paper/kobon_broad_capacity.tex`; `reproduce.md` |
| `scratch/kobon/{reified_face_audit,chirotope_gp_audit}.py` + JSON | Fixed exact-witness acceptance for reified/local cuts; determinant-sign and zero-aware Grassmann--Plücker replay plus base-relaxation gap probe. | MACHINE-CHECKED implementation audits; GP propagation-only through \(n=7\) | `frontier_endpoint_research.json`; `reproduce.md` |
| `scratch/kobon/pappus_relaxation_probe.py` + JSON | Exhibits a guarded non-Pappus abstract-order model and checks that the projective implication clauses reject it. | MACHINE-CHECKED relaxation-gap probe; clauses not used in frontier | `paper_capacity.md`; `reproduce.md` |
| `scratch/kobon/n12_gap_sector_deletion_t39.{cnf,metadata.json}` | Direct-gap target-39 lane with hereditary one/two-line deletion cuts and multipoint sector cuts; metadata pins each smaller-\(n\) bound and the promotion gate. | PENDING; no solver verdict yet | `math/PROGRESS.md`; `ARTIFACTS.md` |
| `scratch/kobon/ladder_n10_t26_faces_{certificate.json,dratcheck.log}` | Hash-bound record and independent checker transcript for the faces-only target-26 monolith; the 14.05 GB DRAT is referenced by digest. | VERIFIED UNSAT; independent second certificate for \(K(10)=25\) | `paper/kobon_broad_capacity.tex`; `reproduce.md` |
| `scratch/kobon/n12/n12_face_cells_q1_t39.discovery.json` | Hash-bound Kissat result for the distilled FACE, unique-parallel-pair, no-concurrency target-39 branch; CNF regeneration is byte-identical. | DISCOVERY UNSAT only; exit 20 is normal, but no proof was logged and no theorem is promoted | `math/PROGRESS.md`; `ARTIFACTS.md` |
| `scratch/kobon/shared_ray_bound_experiment.json` + `n12_gap_shared_ray_deletion_t39.cnf` | Hash-bound proof, clauses, exact audit, A/B benchmarks and live opt-in target-39 solver lane for the shared-ray endpoint cut. | PROVED local lemma; VERIFIED implementation; live lane PENDING with no verdict | `math/PROGRESS.md`; `ARTIFACTS.md`; `reproduce.md` |
| `scratch/kobon/frontier_endpoint_research.json` + `n12_gap_endpoint_mi_k4_{deletion,sub9}_t39.cnf` | Hash-bound recent-source survey, endpoint/perturbation theorems, exact audits, mixed A/B controls, and two reproducible target-39 discovery inputs. | PROVED local cuts + MACHINE-CHECKED audits; both solver processes terminated without verdict or proof | `math/PROGRESS.md`; `ARTIFACTS.md`; `reproduce.md` |
| `scratch/kobon/{one_line_relocation.py,test_one_line_relocation.py,relocation_walk.py}` + relocation JSON censuses | Exact parameter-stratum enumeration of every one-line replacement around the 38-face record and three selected depth-two neighbors; 273,696 placements, maximum 38. | MACHINE-CHECKED exhaustive local basin; not a global upper bound | `frontier_relocation_research.json`; `math/PROGRESS.md`; `experiments.md` |
| `scratch/kobon/{cutpath_insertion.py,test_cutpath_insertion.py,cutpath_line0_n12.json}` | Exact tope-graph enumeration of every simple source-to-antipode insertion into the line-0 deletion base; 7,960 fixed-affine topologies, maximum 36. | MACHINE-CHECKED exhaustive simple fixed-base census; excludes vertex/parallel insertions and is not a global bound | `math/PROGRESS.md`; `experiments.md`; `NEXT_BREAKTHROUGHS.md` |
| `scratch/kobon/{covector_insertion.py,test_covector_insertion.py,covector_line0_n12.json}` | Generalized line-0 deletion search with open-edge, arbitrary disjoint base-vertex, and point-at-infinity parallel transitions; 94,609 topologies, maximum 38. | MACHINE-CHECKED exhaustive generalized simple-parent census; contains all 6,304 exact geometric strata, but is not global | `math/PROGRESS.md`; `experiments.md`; `NEXT_BREAKTHROUGHS.md` |
| `math/kobon/release/kobon-2026-08.zip` | Curated, hash-manifested release; large live solver objects deliberately referenced only. | MACHINE-CHECKED package integrity | `ARTIFACTS.md` |

## Disk/document discrepancies found during this inventory

1. `math/kobon/README.md` is a valid n=10 closeout page but its state heading is
   dated 2026-08-15; it is not the current whole-campaign status page.
2. The release README abbreviates the monolith location as
   `scratch/kobon/n14_monolith_t54/t55.cnf`.  On disk the actual paths are
   `scratch/kobon/n14_monolith_t54.cnf` and
   `scratch/kobon/n14/n14_monolith_t55.cnf` (with t56 beside t55).
3. The release README records the live n=11 DRAT as 78 GB.  At the experiment
   snapshot it had grown to **94,389,665,792 bytes** and was still in flight;
   the old number is a packaging-time snapshot, not a sealed size.
4. The same README calls `n14c-q*.cnf.drat` streams “70–76 GB, frozen.”  At the
   snapshot the three proof-producing processes were active; their streams were
   75,473,354,752 bytes (`q0tp1`), 81,948,311,552 bytes (`q1tp1d`), and
   464,519,168 bytes (`q1tp1s`).  None was complete or verified.
5. `scratch/kobon/n14/concurrency_cases.py` still says the multipoint extension
   is under independent audit and labels every UNSAT “discovery UNSAT.”  The
   capacity theorem is now audited/proved, but the conservative solver label
   remains correct for every live run lacking an independently checked DRAT.
6. No `.log` transcript exists beside either appendix obstruction checker.
   Their sources and documented expected terminal lines are present; only the
   Maiorana and sweep verifiers have copied transcript files.
7. `sweep_all_verify.py` and `max_packing.py` contain source-location-relative
   path calculations.  Their copies under `math/kobon/verification/` are exact
   archival copies, not silently rewritten relocatable variants; use the
   commands in `reproduce.md`.
8. `scratch/kobon/n14/q3pos_index.json` describes a numbered q3-position
   family, but the on-disk CNF set contains 33 files (`000`, `001`, and
   `003`--`033`): `q3pos-002.cnf` is absent.  The artifact inventory records
   the missing expected path rather than inventing a file or experiment.
