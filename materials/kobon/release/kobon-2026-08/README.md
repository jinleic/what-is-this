# Kobon campaign consolidated release — 2026-08-24

Self-contained snapshot of the current main results. Companion canonical zip:
`../kobon-2026-08.zip` (this directory archived). Two earlier drafts
(`*-obsolete-structure{,2}.zip`) are superseded packaging drafts kept only
for transparency; use the canonical zip.

## What is here

```
papers/           Consolidated paper draft + full capacity paper + campaign report
proofs/           Kgen(14)>=54 certificate; Theorem-M obstruction docs
                  (n=14 exact integer certificate; n=18 machine-checked legs;
                  n=20 conditional-ledger) + verification scripts
verification/     fail-closed exact verifiers and audit replays, including
                  Maiorana, endpoint/sector, reified-face, chirotope and Pappus
maiorana/         15 exact-rational witness JSONs + SHA256SUMS.txt manifest +
                  upstream LICENSE (CC BY 4.0 — Andrea Maiorana)
scripts/          campaign engine (exact rational core) + cube builder
docs/             literature pins, orbit/cube ledger, signature algebra,
                  campaign PROGRESS at time of packaging
discoveries/      exact certificates for known lower bounds
MANIFEST.sha256   integrity manifest of every file below
```

## Relation of artifacts to the results

| Result | Backing artifact |
|---|---|
| $K_{\rm gen}(14)\ge 54$ (verified) | `maiorana/` + `verification/maiorana_verify.{py,log}` + `verification/sweep_all_verify.{py,log}` + `proofs/Kgen14_ge54_certificate.md` |
| Theorem M ($n=14$, $Q=3$) | `proofs/escape_obstruction.md` + `proofs/appendix_independent_check.py` |
| $n=18$ $Q=3$ closure | `proofs/obstruction_n18.md` + `proofs/appendix_check_n18.py` |
| $n=20$ $Q=3$ conditional | `proofs/obstruction_n20.md` + `docs/signature_algebra.md` |
| Capacity theorem (main) | `papers/paper_capacity.md` §§1–3 (corrected identities) |
| Full machinery + ledger | `papers/paper_capacity.md` §§4–8, `papers/report.md` |
| Direct-gap face criterion | `scripts/gap_faces.py` + `verification/gap_face_verify.{py,json}` + `certificates/gap_n7_t12.*` |
| $n=12$, $Q=1$, no-concurrency branch UNSAT | `scripts/q1_gap_faces.py` + `certificates/n12_q1_gap_certificate.json` + CNF + checker transcript |
| Saturated-flower proposition / C--B audit | `scripts/flower_counterexamples.py` + `verification/flower_counterexamples.json` |
| Face square-penalty refutation | `scripts/square_penalty_counterexample.py` + `verification/square_penalty_counterexample.json` + `papers/kobon_broad_capacity.{tex,pdf}` |
| Coefficient-$2/3$ endpoint at $n=7$ | `scripts/square_penalty_sat.py` + `certificates/c23_n7_{certificate.json,violation.cnf,violation.drat,violation.dratcheck.log}` |
| Hereditary and vertex-sector face cuts | `scripts/{engine,gap_faces}.py` + `verification/sector_bound_audit.{py,json}` |
| Abstract-order stretchability gap (dual Pappus) | `verification/pappus_relaxation_probe.{py,json}`; probe only, clauses not used in frontier |
| Independent faces-only $K(10)=25$ certificate | `certificates/ladder_n10_t26_faces_{certificate.json,cnf,dratcheck.log}` + `scripts/faces_only.py`; the 14.05 GB proof is hash-bound externally |
| Shared-ray endpoint-concurrency lemma and opt-in lane | `scripts/{engine,gap_faces}.py` + `verification/sector_bound_audit.{py,json}` + `discoveries/shared_ray_bound_experiment.json`; live CNF referenced by hash |
| Endpoint closure, perturbation incidence, and live \(n=12\) lanes | `scripts/{engine,gap_faces}.py` + `verification/{sector_bound_audit,reified_face_audit,chirotope_gp_audit}.{py,json}` + `discoveries/frontier_endpoint_research.json`; live CNFs referenced by hash |

## Deliberately NOT included (referenced, not copied)

- Large or live SAT objects: the two hash-bound $n=12$ Q1 proofs:
  primary unsymmetrized `n12_q1_gap_t39.drat` (2,437,053,136 bytes,
  SHA-256 `11ce6f6f8c31621c5ecff355dfd23678f1e21e572ea5e43349397b2867ec5545`)
  and independent symmetry-broken `n12_q1_gap_sym_t39.drat`
  (759,328,450 bytes, SHA-256
  `c6f7eaed6355a9a4adbb3dc7aaa9b1ff1b9616ab8842e088733d107dcb1df088`).
  Both CNFs, the dual-certificate manifest and both `s VERIFIED` checker
  transcripts are included;
  the in-flight n=11 proof stream, frozen 70--82 GB n=14 cube proofs, and
  40 MB--1 GB frontier CNFs. The new endpoint target-39 CNFs are also
  referenced rather than copied: 21,685,650 and 39,495,970 bytes. Their
  generators and hash-bound discovery records are included.
- The complete independent $n=10$, target-26 face-monolith DRAT is referenced,
  not copied: `scratch/kobon/ladder_n10_t26_faces.drat`
  (14,050,806,440 bytes, SHA-256
  `bdbab47463806ee96163d03ec36a052e5818e2125202272182ee4cbe94629007`).
  Its byte-reproducible 22,125,971-byte CNF, hash-bound certificate JSON and
  `s VERIFIED` checker transcript are included.
- Live process state, `*.kissat.log` (rotate and enormous).

## Verify the core claim yourself

```
scratch/kobon-audit/venv/bin/python verification/maiorana_verify.py \
    maiorana/sol1_lines_rational.json
scratch/kobon-audit/venv/bin/python verification/sweep_all_verify.py
scratch/kobon-audit/venv/bin/python scripts/square_penalty_counterexample.py
scratch/kobon-audit/venv/bin/python verification/sector_bound_audit.py
scratch/kobon-audit/venv/bin/python verification/pappus_relaxation_probe.py
```

All five commands must exit 0. The verifier imports need the PySAT-capable
campaign environment and `scripts/engine.py` on the Python path; the packaged
sector replay also imports `scripts/square_penalty_counterexample.py`.

## Provenance note

`maiorana/` copies CC BY 4.0 data (LICENSE included) from
https://github.com/rufio72/kobon_triangles_k14 `e47c7cfd9661e54`.
All hashes in `maiorana/SHA256SUMS.txt` match upstream at capture time;
this workspace's verification is the independent confirmation cited in the
papers. The OEIS row recording 54 at $n=14$ is independently supported here
from below by the exact-rational Maiorana data. Its upper-bound provenance is
audited separately in `docs/litrefs.md`: simple-only BBL/Blanc bounds do not
apply, the Clément--Bader draft contains a false local proof step, and
Rasukaru's 2005 web argument is important but informal.

The regular simplicial family used in the square-penalty refutation is pinned
to Grünbaum, *A catalogue of simplicial arrangements in the real projective
plane*, Ars Math. Contemp. 2 (2009), 1–25,
doi:10.26493/1855-3974.88.e12. The unrelated arXiv identifiers 0904.1244 and
1011.1862 are explicitly rejected in `docs/litrefs.md`; Cuntz's correct
identifier is arXiv:1108.3000.

## v3 additions (2026-08-21/22)

```
papers/kobon_broad_capacity.{tex,pdf}   the consolidated research paper (18 pp)
docs/                                   documentation hub: INDEX, methods,
                                        experiments, reproduce, ARTIFACTS,
                                        NEXT_BREAKTHROUGHS
verification/max_packing.py             exact broad-optimum prober (overlap
                                        graph + SAT-certified optimality)
verification/broad_ladder.py            (n,T) ladder prober, all degeneracies free
certificates/ladder_n8_t16.dratcheck.log   s VERIFIED  -> K_gen(8) = 15
certificates/ladder_n9_t22.dratcheck.log   s VERIFIED  -> K_gen(9) = 21
certificates/*.maxpacking.json          exact optima of all 20 known witnesses
```

Large objects still referenced, not copied: `scratch/kobon/ladder_n8_t16.drat`
(50.5 MB), `ladder_n9_t22.drat` (328.5 MB), the monolith/cube CNFs
(240 MB-1 GB each) and the in-flight n=11 proof stream.

## v9 additions (2026-08-22)

```
scripts/gap_faces.py                    universal direct-gap face encoder
scripts/q1_gap_faces.py                 canonical Q=1/no-concurrency encoder
scripts/q1_deficit_cubes.py              exhaustive one-deficit cover builder
scripts/face_regions.py                  stronger FACE / distilled variants
scripts/flower_counterexamples.py        exact k=3,4 saturated flowers
verification/gap_face_verify.{py,json}   109,104-triple direct-gap oracle audit
verification/flower_counterexamples.json exact face/shared-side census
certificates/gap_n7_t12.{cnf,drat,dratcheck.log}
                                         compact direct-gap `s VERIFIED` control
certificates/n12_q1_gap_certificate.json  hash-bound dual branch certificate
certificates/n12_q1_gap_t39.{cnf,dratcheck.log}
                                         primary 2.44 GB proof VERIFIED
certificates/n12_q1_gap_sym_t39.{cnf,dratcheck.log}
                                         independent 759 MB proof VERIFIED
docs/novelty_audit_capacity_theorem.md    Rasukaru-2005 prior-art addendum
```

The raw 2.44 GB and 759 MB Q1 proofs are intentionally outside the zip; their
exact sizes and SHA-256 digests are bound by the included JSON manifest and
reproduced above. The
included CNF plus generator is sufficient to rerun Kissat and compare a fresh
proof through `drat-trim`.

## v10 additions (2026-08-23)

```
scripts/square_penalty_counterexample.py
                                         exact rational A(12,1) constructor,
                                         dual face oracles, edge/sector census
verification/square_penalty_counterexample.json
                                         F=30, Q=0, N3=15, N6=1;
                                         coefficient-1 square bound fails 90>89
scripts/square_penalty_sat.py             rational-weight direct-gap violation
                                         generator; selected-face count and
                                         every degeneracy free
certificates/c23_n7_violation.{cnf,drat,dratcheck.log}
                                         74,507 vars / 152,419 clauses;
                                         34,133,280-byte proof, s VERIFIED
certificates/c23_n7_certificate.json      hashes, scaled theorem and checker core
papers/kobon_broad_capacity.{tex,pdf}   18-page draft with edge-use and
                                         sector-run identities plus infinite
                                         A(2m,1), m>=6 counterexample theorem
papers/paper_capacity.md                 companion square-penalty post-mortem
docs/{litrefs,novelty_audit_capacity_theorem}.md
                                         corrected proof-path provenance
docs/PROGRESS.md                         newest-first verified frontier entry
```

The counterexample is geometry, not SAT evidence: every coordinate and
predicate uses exact rational arithmetic, `engine.verify_selection` returns
`ok`, and the independent edge census gives
`(M,E_1,E_2,E_0)=(51,12,39,0)`. The infinite family also proves that every
fixed coefficient greater than `2/3` in front of
`sum_p (k_p-2)^2` fails asymptotically.

## v11 additions (2026-08-23)

```
scripts/{engine,gap_faces}.py            hereditary deletion and local
                                         vertex-sector SAT accelerators
verification/sector_bound_audit.{py,json}
                                         3,600 exact arrangements;
                                         93,000 pair checks, zero violations
verification/pappus_relaxation_probe.{py,json}
                                         guarded non-Pappus abstract model
                                         and rejecting projective cut
```

The combined \(n=12,T=39\) discovery CNF is intentionally referenced rather
than copied. Its live solvers have no verdict; these additions do not change
the proved window for \(K(12)\).

## v12 additions (2026-08-24)

```
certificates/ladder_n10_t26_faces.cnf
                                         byte-reproducible 131,314-variable /
                                         792,462-clause face monolith
certificates/ladder_n10_t26_faces.dratcheck.log
                                         independent drat-trim transcript;
                                         3,464,683,123 resolution steps;
                                         s VERIFIED
certificates/ladder_n10_t26_faces_certificate.json
                                         SHA-256 pins, scope, regeneration,
                                         solver and checker statistics
discoveries/n12_gap_sector_deletion_t39.metadata.json
                                         explicit K(11)/K(10) provenance,
                                         live-discovery status and proof gate
```

The $n=12$ hereditary lane was already unconditional: its $K(10)=25$ input is
the 2026-08-15 checked eleven-cube cover. It never depended on the face
monolith. The newly checked monolith is a second independent certificate.
The $n=12$ lane itself still has no verdict and remains discovery-only.

## v13 additions (2026-08-24)

```
discoveries/n12_face_cells_q1_t39.discovery.json
                                         byte-identical distilled FACE Q=1 /
                                         no-concurrency instance provenance;
                                         Kissat exit 20, s UNSATISFIABLE
```

This is deliberately labelled `DISCOVERY_UNSAT_NO_PROOF`: the run emitted no
DRAT, so it is not promoted as a new certificate. The same Q=1/no-concurrency
geometric branch already has two checked direct-gap certificates. The
concurrency-inclusive \(n=12,T=39\) frontier remains open.

## v14 additions (2026-08-24)

```
scripts/{engine,gap_faces}.py            opt-in shared-ray endpoint clauses;
                                         sector-only default unchanged
verification/sector_bound_audit.{py,json}
                                         1,591 same-ray cases, zero violations;
                                         exact four-sector sharpness control
discoveries/shared_ray_bound_experiment.json
                                         proof, clause schema, review,
                                         mixed A/B benchmarks and live-lane pin
papers/kobon_broad_capacity.{tex,pdf}    shared-ray lemma and proof
```

The live `n12_gap_shared_ray_deletion_t39.cnf` is referenced rather than
copied: 374,381 variables / 841,910 clauses, SHA-256
`0aa9e81e1806d5ff7413669d7357454c9e08a496d2172d8b4a512b27115b7e7d`.
It has no verdict and no proof logging. Performance controls were mixed, so
the cut remains opt-in and is presented as solver diversity, not a speedup.

## v15 additions (2026-08-24)

```
scripts/{engine,gap_faces}.py            endpoint closure; simple-perturbation
                                         incidence; K(4) prime clauses;
                                         opt-in exact-face and chirotope probes
verification/sector_bound_audit.{py,json}
                                         1,607 endpoint collisions, 84 new
                                         same-ray cases, 42 single-line and
                                         1,591 shared-pair segments; 0 failures
verification/reified_face_audit.{py,json}
                                         exact K(10)/K(12) fixed-geometry PASS
verification/chirotope_gp_audit.{py,json}
                                         195,000 signs / 1,386,000 GP
                                         relations; propagation-only through n=7
discoveries/frontier_endpoint_research.json
                                         source survey, theorem proofs, clause
                                         deltas, A/B controls and live hashes
papers/kobon_broad_capacity.{tex,pdf}    endpoint and perturbation lemmas
```

The primary endpoint/multipoint/K(4) lane has 375,037 variables /
1,005,366 clauses, SHA-256
`6dca104aa120ee1f209aada8aca07a08bf3341f1335e6e99f4ae1f11bc2ba76b`.
The \(K(9)=21\) hereditary variant has 842,317 / 1,917,926, SHA-256
`34406004c932ce77d33f76f8b7d19c7b0f5a2e3e241e94651846269bf750f2df`.
Both are live discovery runs without proof logging. Exact face reification and
zero-aware Grassmann--Plücker were retained only as negative/propagation
experiments after mixed or adverse controls. \(K(12)\) remains open.
