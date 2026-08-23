# h10q deliverable bundle --- 2026-08-19

Two paper drafts plus the artifacts, generators and ledger they cite.

## Provenance labels

- **repo-regenerable**: a checked-in script in `scripts/` regenerates the artifact.
- **frozen authority (suite-asserted)**: byte-identical match is enforced by the test suites.
- **persisted evidence (no checked-in producer)**: reviewed data whose original
  session-side generator was not recovered into the repository.
- **generator/replay script**, **paper draft**, **project ledger**: as named.

## Contents

| file | provenance | bytes |
|---|---|---|
| `artifacts/l13h_all_closures.json` | frozen authority (suite-asserted) | 85082 |
| `artifacts/l15_density.json` | persisted evidence (no checked-in producer) | 53072 |
| `artifacts/l15_remainders.jsonl` | persisted evidence (no checked-in producer) | 234709 |
| `artifacts/l16_char.jsonl` | persisted evidence (no checked-in producer) | 131104 |
| `artifacts/l17_badroots.jsonl` | persisted evidence (no checked-in producer) | 23289601 |
| `artifacts/l17_badroots_closures.jsonl` | repo-regenerable | 25323967 |
| `artifacts/l17_classfactors.jsonl` | repo-regenerable | 79486 |
| `artifacts/l17_cofactor.jsonl` | persisted evidence (no checked-in producer) | 496513 |
| `artifacts/l17_cofactor_nonclosure.jsonl` | persisted evidence (no checked-in producer) | 201567 |
| `artifacts/l17_composite_w.jsonl` | repo-regenerable | 9636 |
| `artifacts/l17_horizon101.jsonl` | repo-regenerable | 143596 |
| `artifacts/l17_horizon103.jsonl` | repo-regenerable | 12003 |
| `artifacts/l17_horizon107.jsonl` | repo-regenerable | 496 |
| `artifacts/l17_horizon109.jsonl` | repo-regenerable | 133326 |
| `artifacts/l17_horizon109_crossk.jsonl` | repo-regenerable | 162073 |
| `artifacts/l17_horizon113.jsonl` | repo-regenerable | 170662 |
| `artifacts/l17_horizon127.jsonl` | repo-regenerable | 25131 |
| `artifacts/l17_horizon131.jsonl` | repo-regenerable | 200559 |
| `artifacts/l17_horizon137.jsonl` | repo-regenerable | 139772 |
| `artifacts/l17_matched.jsonl` | repo-regenerable | 112630 |
| `artifacts/l17_ratemodel_censored.jsonl` | persisted evidence (no checked-in producer) | 32723474 |
| `artifacts/l17_schinzel_audit.jsonl` | repo-regenerable | 385432 |
| `artifacts/l17_sieve.jsonl` | persisted evidence (no checked-in producer) | 1950792 |
| `artifacts/l17_stepii.jsonl` | repo-regenerable | 631803 |
| `artifacts/l18_divergence_model.jsonl` | repo-regenerable | 629784 |
| `artifacts/l18_fivewall.jsonl` | repo-regenerable | 123574 |
| `artifacts/l18_horizon_sweep_a.jsonl` | repo-regenerable | 522290 |
| `artifacts/l18_horizon_sweep_b.jsonl` | repo-regenerable | 392941 |
| `artifacts/l18_route131.jsonl` | repo-regenerable | 413188 |
| `artifacts/l18_schinzel_implies_h.jsonl` | repo-regenerable | 23440 |
| `artifacts/l19_cell179.jsonl` | repo-regenerable | 12013190 |
| `artifacts/l19_classexist.jsonl` | repo-regenerable | 405756 |
| `artifacts/l19_tauzero.jsonl` | repo-regenerable | 568166 |
| `artifacts/l20_admissible.jsonl` | repo-regenerable | 808659 |
| `artifacts/l21_irred_direct.jsonl` | repo-regenerable | 643542 |
| `artifacts/l21_irred_generic.jsonl` | repo-regenerable | 1818481 |
| `artifacts/l21_reducible_locus.jsonl` | repo-regenerable | 19466 |
| `artifacts/l22_elimination.jsonl` | repo-regenerable | 1053714 |
| `artifacts/l22_factor_tuple.jsonl` | repo-regenerable | 11463 |
| `artifacts/l22_fiber_geometry.jsonl` | repo-regenerable | 7158 |
| `artifacts/l22_infinity.jsonl` | repo-regenerable | 4870 |
| `artifacts/l22_reciprocal_cube.jsonl` | repo-regenerable | 17615 |
| `artifacts/l22_square_branch.jsonl` | repo-regenerable | 19003 |
| `artifacts/l6_witnesses.jsonl` | frozen authority (suite-asserted) | 3785 |
| `artifacts/l9_steered.jsonl` | frozen authority (suite-asserted) | 23025 |
| `artifacts/litscout_h10q.md` | persisted evidence (no checked-in producer) | 8438 |
| `ledger/CONDITIONAL.md` | project ledger | 19600 |
| `ledger/NOTES.md` | project ledger | 73220 |
| `ledger/README.md` | project ledger | 27786 |
| `ledger/RESULTS.md` | project ledger | 50612 |
| `ledger/THEOREMS.md` | project ledger | 140675 |
| `papers/companion-verification.tex` | paper draft | 35144 |
| `papers/main-conditional-forall6.tex` | paper draft | 9462 |
| `papers/sections/architecture.tex` | paper draft | 10004 |
| `papers/sections/classes.tex` | paper draft | 8324 |
| `papers/sections/existence.tex` | paper draft | 6483 |
| `papers/sections/frontier.tex` | paper draft | 9171 |
| `papers/sections/intro.tex` | paper draft | 9570 |
| `papers/sections/irreducibility.tex` | paper draft | 14175 |
| `papers/sections/record.tex` | paper draft | 11208 |
| `papers/sections/schinzel.tex` | paper draft | 8203 |
| `papers/sections/sieve.tex` | paper draft | 8499 |
| `papers/sections/walls.tex` | paper draft | 10010 |
| `scripts/h10q.py` | generator/replay script | 198192 |
| `scripts/l12_class.py` | generator/replay script | 16371 |
| `scripts/l13_filter.py` | generator/replay script | 40427 |
| `scripts/l13h_scan.py` | generator/replay script | 6566 |
| `scripts/l14_replay_all.py` | generator/replay script | 3595 |
| `scripts/l16_emergent.py` | generator/replay script | 5704 |
| `scripts/l17_badroots_closures.py` | generator/replay script | 7650 |
| `scripts/l17_classrisk.py` | generator/replay script | 17233 |
| `scripts/l17_composite_w.py` | generator/replay script | 20637 |
| `scripts/l17_horizon101.py` | generator/replay script | 13495 |
| `scripts/l17_horizon103.py` | generator/replay script | 13355 |
| `scripts/l17_horizon107.py` | generator/replay script | 7751 |
| `scripts/l17_horizon109.py` | generator/replay script | 17760 |
| `scripts/l17_horizon109_crossk.py` | generator/replay script | 6953 |
| `scripts/l17_horizon113.py` | generator/replay script | 20060 |
| `scripts/l17_horizon127.py` | generator/replay script | 15699 |
| `scripts/l17_horizon131.py` | generator/replay script | 14541 |
| `scripts/l17_horizon137.py` | generator/replay script | 21023 |
| `scripts/l17_matched.py` | generator/replay script | 10650 |
| `scripts/l17_ratemodel_padic.py` | generator/replay script | 24292 |
| `scripts/l17_schinzel_audit.py` | generator/replay script | 18870 |
| `scripts/l18_divergence_model.py` | generator/replay script | 42896 |
| `scripts/l18_fivewall.py` | generator/replay script | 25178 |
| `scripts/l18_horizon_sweep_a.py` | generator/replay script | 53037 |
| `scripts/l18_horizon_sweep_b.py` | generator/replay script | 30368 |
| `scripts/l18_route131.py` | generator/replay script | 39902 |
| `scripts/l18_schinzel_implies_h.py` | generator/replay script | 27822 |
| `scripts/l19_cell179.py` | generator/replay script | 41306 |
| `scripts/l19_classexist.py` | generator/replay script | 39638 |
| `scripts/l19_tauzero.py` | generator/replay script | 38085 |
| `scripts/l20_admissible.py` | generator/replay script | 40213 |
| `scripts/l21_irred_direct.py` | generator/replay script | 29699 |
| `scripts/l21_irred_generic.py` | generator/replay script | 30789 |
| `scripts/l21_reducible_locus.py` | generator/replay script | 8891 |
| `scripts/l22_elimination.py` | generator/replay script | 43490 |
| `scripts/l22_factor_tuple.py` | generator/replay script | 38235 |
| `scripts/l22_fiber_geometry.py` | generator/replay script | 47516 |
| `scripts/l22_infinity.py` | generator/replay script | 33552 |
| `scripts/l22_reciprocal_cube.py` | generator/replay script | 32683 |
| `scripts/l22_square_branch.py` | generator/replay script | 46150 |

## Verification

```sh
shasum -a 256 -c SHA256SUMS      # from inside this directory
```

## Reproduction (repo-regenerable chain, run from math/h10q)

```sh
python3 l17_ratemodel_padic.py            # p-adic rate model / bad-mass engine
python3 l17_badroots_closures.py          # closure-authority bad-root tables
python3 l17_matched.py                    # matched statistics artifact
python3 l17_classrisk.py                  # per-class window factors
python3 l17_schinzel_audit.py             # Schinzel local-condition audit
python3 l17_horizon101.py                 # off-grid closure replays (103/107/109 likewise)
python3 l17_horizon109_crossk.py          # ~17 min deterministic cross-k scan
python3 l19_tauzero.py                   # tau=0 wall criterion / class-side break
python3 l19_classexist.py                # uniform aligned-class construction
python3 l19_cell179.py                   # exact w=179 member closure
python3 l20_admissible.py                # uniform admissibility audit
nice -n 19 python3 l21_reducible_locus.py # explicit reducible locus / branch escape
nice -n 19 python3 l21_irred_generic.py # generic and per-fiber irreducibility
nice -n 19 python3 l21_irred_direct.py  # local no-go laws / reciprocal subfamily
nice -n 19 python3 l22_elimination.py # PROVED fixed-Z theorem; scans EVIDENCE
nice -n 19 python3 l22_reciprocal_cube.py # PROVED independent theorem; scans EVIDENCE
nice -n 19 python3 l22_infinity.py # PROVED local structure; route OPEN
nice -n 19 python3 l22_factor_tuple.py # PROVED criterion/no-go; compatibility OPEN
nice -n 19 python3 l22_square_branch.py # PROVED free-lambda theorem; six-count OPEN/excluded
nice -n 19 python3 l22_fiber_geometry.py # PROVED reductions; uniform lemma OPEN
python3 h10q.py                           # default suite (exit 0 required)
python3 h10q.py --extended                # extended suite (exit 0 required)
```
