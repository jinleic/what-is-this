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
| `artifacts/l19_tauzero.jsonl` | repo-regenerable | 568168 |
| `artifacts/l20_admissible.jsonl` | repo-regenerable | 808659 |
| `artifacts/l6_witnesses.jsonl` | frozen authority (suite-asserted) | 3785 |
| `artifacts/l9_steered.jsonl` | frozen authority (suite-asserted) | 23025 |
| `artifacts/litscout_h10q.md` | persisted evidence (no checked-in producer) | 8438 |
| `ledger/CONDITIONAL.md` | project ledger | 23633 |
| `ledger/NOTES.md` | project ledger | 66528 |
| `ledger/README.md` | project ledger | 21245 |
| `ledger/RESULTS.md` | project ledger | 48370 |
| `ledger/THEOREMS.md` | project ledger | 123532 |
| `papers/companion-verification.tex` | paper draft | 27586 |
| `papers/main-conditional-forall6.tex` | paper draft | 9431 |
| `papers/sections/architecture.tex` | paper draft | 9549 |
| `papers/sections/classes.tex` | paper draft | 8186 |
| `papers/sections/existence.tex` | paper draft | 6483 |
| `papers/sections/frontier.tex` | paper draft | 9153 |
| `papers/sections/intro.tex` | paper draft | 9286 |
| `papers/sections/record.tex` | paper draft | 10558 |
| `papers/sections/schinzel.tex` | paper draft | 7469 |
| `papers/sections/sieve.tex` | paper draft | 8461 |
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
python3 h10q.py                           # default suite (exit 0 required)
python3 h10q.py --extended                # extended suite (exit 0 required)
```
