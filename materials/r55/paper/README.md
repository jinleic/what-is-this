# Engström mixed identity at n=45 — paper and artifact bundle

**Result:** `MIXED_NO_CUT_IN_FROZEN_BASIS` — the recorded first exact identity
basis for hypothetical Ramsey(5,5,45) graphs (McKay–Radziszowski m=2,3,4 plus
Engström's non-separable degree-4 identity) is exhausted, with exact rational
primal witnesses for every instantiable route. **No claim is made about R(5,5).**

| paper | [`engstrom_mixed_n45.tex`](engstrom_mixed_n45.tex) → [`engstrom_mixed_n45.pdf`](engstrom_mixed_n45.pdf) (6 pages) |
|---|---|
| canonical artifact | `../data/engstrom_identity.json` (schema 3, 4.4 MB) |
| terminal record | `../notes/mixed_campaign_result_2026-08-21.md` |
| implementation plan | `../notes/engstrom_identity_mixed_plan_2026-08-18.md` |

## Headline numbers

| quantity | value |
|---|---|
| states in the relaxation | 3,215 |
| catalog graphs streamed | 8,500,211 (4348.2 s, budget 7200 s) |
| endpoint rank without / with F | 6 → 8 (**row is geometrically new**) |
| dual & primal optima with / without F | identical (**zero cut strength**) |
| exact route bounds = witnesses | 360 · 14310/349 · 45 |
| median F window width vs median \|F\| | 802,396 vs 664,870 (ratio 1.207) |
| ρ(h_x, F_x) on class (22,114) | 0.9464 (structure the box discards) |
| envelope classes: share of all window width | **98.6%** (77 of 130 classes) |
| envelope looseness vs catalog truth | **111×** (width 92,888 vs 836) |
| routes with envelope windows collapsed to points | **360 / 39.25 / 38.36** vs edges 315 / 1 / 1 |

The last row is the decisive finding: **perfect motif certification would still
leave a ≈39× gap**, so the bottleneck is not data quality but the absence of
vertex-to-vertex coupling. Three natural follow-on campaigns (catalog-cloud
hulls, structural caps on missing classes, exact envelope projection) are ruled
out by measurement; a level-2 pair lift is the one direction the evidence
supports. See the paper §5 and the terminal record §6.

## Reproduce

```bash
cd <repo>/math
./.venv/bin/python -m unittest discover -s r55/tests -p 'test_*mixed*.py' -v
./.venv/bin/python r55/src/mixed_deficiency_cone.py --output r55/data/engstrom_identity.json   # ~75 min
./.venv/bin/python r55/src/search_mixed_cuts.py     --analysis r55/data/engstrom_identity.json # seconds
./.venv/bin/python r55/src/verify_mixed_independent.py       r55/data/engstrom_identity.json   # sub-second
```

`verify_mixed_independent.py` imports **no** producer module (stdlib only): it
re-hashes every recorded input, re-derives the catalog counts, recomputes all
3,215 F intervals from the artifact's own windows, and re-derives every
certificate, witness, route status and the disposition.

## Verification status (explicit)

Verified: kernel exhaustive to n≤6 (33,867 graphs) plus 656-graph replay and the
n=49 type-combo replay; catalog envelope self-audit (0 violations); 3,215-state
projection equality against the frozen m4 artifact; 24 focused contract tests;
independent re-verification of the search and F-interval layers; a 12-mutation
corruption battery rejected item-by-item.

**Not** yet independently redone (Task 4 tier): a second, independently written
motif kernel re-sweeping the 8.5 M graphs; recomputation of the g and h
endpoints from the R1–R5 relation system; independent re-derivation of the
envelope-LP windows for the 77 catalog-missing stratum classes.

## Next campaign

Replace per-row interval boxes with the exact Minkowski-sum hull of the joint
(g, h, F) point clouds per stratum — strictly tighter, same trust base, no new
identity or catalog. Falsification criterion: if the lifted programme still
returns 360, 14310/349, 45, the interval-free geometry of this basis is also
closed. See the paper's §5 and the terminal record's §6.
