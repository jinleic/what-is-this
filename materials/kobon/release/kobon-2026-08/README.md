# Kobon campaign consolidated release — 2026-08-21

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
verification/     fail-closed Maiorana certifier (+ exit-0 transcript logs);
                  replayable 15-solution sweep verifier (+ log)
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

## Deliberately NOT included (referenced, not copied)

- Byte-scale live solver objects: `scratch/kobon/n11_monolith_t33.drat` (78 GB,
  in-flight proof), `scratch/kobon/n14/n14c-q*.cnf.drat` (70-76 GB, frozen),
  `scratch/kobon/n14_monolith_t54/t55.cnf` (240 MB each), `scratch/kobon/n17..20` cube CNFs (~1 GB each).
  These regenerate from `scripts/concurrency_cases.py` /
  `scratch/kobon/nOO/cube_build.py` within minutes-hours; the DRATs are
  in-flight proof streams sealed only by live solver runs.
- Live process state, `*.kissat.log` (rotate and enormous).

## Verify the core claim yourself

```
scratch/kobon-audit/venv/bin/python verification/maiorana_verify.py \
    maiorana/sol1_lines_rational.json
scratch/kobon-audit/venv/bin/python verification/sweep_all_verify.py
```

Both must exit 0. The engine import inside `maiorana_verify.py` needs
`math/kobon/engine.py` (repository workspace, or `scripts/engine.py` here
placed on the Python path as package `engine`).

## Provenance note

`maiorana/` copies CC BY 4.0 data (LICENSE included) from
https://github.com/rufio72/kobon_triangles_k14 `e47c7cfd9661e54`.
All hashes in `maiorana/SHA256SUMS.txt` match upstream at capture time;
this workspace's verification is the independent confirmation cited in the
papers. The OEIS A006066 approved row `14 >= 53 54 [Bader]` is unchanged;
the proposed 54 row is unapproved OEIS history (see `docs/litrefs.md`).

## v3 additions (2026-08-21/22)

```
papers/kobon_broad_capacity.{tex,pdf}   the consolidated research paper (8 pp)
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
