# Figure provenance: `reports/fig_rate_distance.pdf`, `reports/fig_trade_law.pdf`

Generated 2026-08-18 with matplotlib 3.11.1 (pre-existing in `.venv`; no package installs), PDF, 300 dpi, `bbox_inches="tight"`. Every numeric claim below was recomputed from the artifacts in this session; the plotting run asserts all counts before saving.

Certified-data inputs (nothing else feeds the ordinates except the explicitly grey "untrusted" markers):

- `results/partial_runs/exp039_nogo_module.json` (schema `exp039-nogo-module-v1`, encoding `sat-decide-v2`)
- `results/partial_runs/exp039/parent_*.json` — 134 parent certificates (same schema); all have `d_z_exact=true`, `T_is_exact=true`, and `k_ceiling = k_parent - T` (asserted)
- `results/partial_runs/exp037/row_*.json` — 308 per-row records, all with `verification.valid=true`
- `results/processed/exp040_saturation_probe.json`
- Catalogue rows only via `E27.load_catalogue()` (`experiments/exp027_delta_audit.py`, imported with the same `importlib` pattern as `experiments/exp039_nogo_module.py`); used for row positions, parent fingerprints, and — for 43 rows — the grey *untrusted* `d` claim, never for a certified ordinate.

Fingerprint lengths: summary rows in `exp039_nogo_module.json` carry 16-char prefixes; certificates carry the full 64-char fingerprint. Joined by prefix (map verified 1:1 over the 134 certs).

## `reports/fig_rate_distance.pdf` — marker-by-marker

x = `n` (log scale, i.e. `matplotlib xscale("log")`; warranted because n spans 36→360), y = certified distance ordinate (linear), marker area ∝ k. Counts in parentheses are asserted at plot time.

| marker | n_x marker count | x | y | artifact + exact keys |
|---|---|---|---|---|
| ● black filled — CSS parent, family closed (61) | 61 | `family_closed_parents[].n` | `family_closed_parents[].d_z_parent` | `exp039_nogo_module.json → family_closed_parents[]` keys `{fingerprint, n, d_z_parent, k_parent, T, num_witnesses}`; each joined to `exp039/parent_*.json` on fingerprint prefix and cross-checked (`family_closed=true`, matching `n`, `d_z_parent`) |
| ○ black open — CSS parent, not closed (73) | 73 | cert `n` | cert `d_z_parent` | `results/partial_runs/exp039/parent_*.json` keys `{n, d_z_parent, k_parent, family_closed:false}`; 73 = 134 certified − 61 closed |
| ▼ red filled — PBB capped a priori (249) | 249 | `rows_capped_a_priori[].n` | `rows_capped_a_priori[].d_z_parent` | `exp039_nogo_module.json → rows_capped_a_priori[]` keys `{catalogue_index, label, n, k_parent, k_pbb, T, k_ceiling, d_z_parent}`. Ordinate is the certified Theorem-H cap d ≤ d_Z(P); asserted ≤ the witness `upper_bound` for all 232 of the 249 that also have an EXP-037 record |
| △ red open — PBB theorem-permitted (30) | 30 | `rows_theorem_permits[].n` | `upper_bound` | `exp039_nogo_module.json → rows_theorem_permits[]` (`catalogue_index`, `n`, `k_pbb`), joined on `catalogue_index` to `results/partial_runs/exp037/row_*.json → upper_bound` with `verification.valid=true` (all 30 present) |
| ★ red star — certified reversal (7) | 7 | `falsification_gate.checks[].n` | `upper_bound` | `exp039_nogo_module.json → falsification_gate.checks[]` keys `{label, catalogue_index, n, k_parent, k_pbb, d_z_parent, d_pbb_certified_gt, T, k_ceiling, slack, hypothesis_verified, permitted_by_theorem_H}` (all 7: `hypothesis_verified=true`, `permitted_by_theorem_H=true`, `slack=0`); star ordinate = EXP-037 `upper_bound`, so the plotted star sits at the top of the certified interval `(d_pbb_certified_gt, upper_bound]`. Red link segments run from the parent point to each star; parent fingerprint recomputed per row via `E27.parent_matrices` + `E27.matrix_fingerprint` and matched to its cert (asserted equal `(n, k_parent, d_z_parent, T)`) |
| □ green open — Gross [[144,12,12]] (1) | 1 | 144 | 12 | Published Gross polynomials A = x³+y+y², B = y³+x+x² on ℓ=12, m=6 match the parent of catalogue rows 41 and 82. Its parent certificate is `results/partial_runs/exp039/parent_aebb649586c1f5ab*.json`, fingerprint `aebb649586c1f5ab68730884255244dba345b922148e469eacb7b0001865f52d`: `d_z_parent=12, T=12, k_parent=12, family_closed=true`; EXP-037 pool record `results/partial_runs/exp037/css_aebb649586c1f5ab.json`: `distance=12, exact=true, encoding_version="sat-decide-v2"` (min(d_X,d_Z)=12 exact, with `d_Z=12` from the certificate) |
| ▽ grey — parent uncertified, witness UB (46) | 46 | record `n` | `upper_bound` | `results/partial_runs/exp037/row_*.json` keys `{catalogue_index, n, k, upper_bound, verification.valid=true}` for rows absent from both EXP-039 classification lists (parent never certified; all such rows have n ∈ {180, 360}) |
| ✕ grey — parent uncertified, catalogue d **untrusted** (43) | 43 | catalogue `n` | catalogue `d` | `E27.load_catalogue()` row keys `n`, `d` — the only non-certified ordinate in the figure; explicitly grey and named untrusted (43 rows lack both an EXP-039 classification and an EXP-037 record) |

Coverage: 249 + 30 + 46 + 43 = 368 catalogue rows (asserted). Markers for the two CSS classes are parents, not rows (134 of the 202 distinct parents are certified; the 68 uncertified parents appear only implicitly through their 89 grey rows).

Star ordinates and the certified intervals: `12_6_0217` (144, UB 10, d>9), `phase2_58` (72, UB 6, d>4), `phase2_60` (72, UB 6, d>4), `phase2_71` (108, UB 6, d>5), `phase2_72` (108, UB 6, d>5), `9_6_0183` (108, UB 10, d>9), `phase2_88` (108, UB 6, d>5). Parent points sit at `d_z_parent` ∈ {4, 8}; every star is strictly above its parent.

Other layout notes: horizontal jitter is deterministic (`md5` of the row key, ±0.004–0.006 log₁₀ decades) and cosmetic; reversal labels are grouped for coincident stars; the provenance footnote is embedded inside the figure.

## `reports/fig_trade_law.pdf` — element-by-element

Left panel (scatter, 279 classified rows, diagonal y = x, cell counts overlaid):

- x = `T`, y = `k_parent − k_pbb` from `exp039_nogo_module.json → rows_capped_a_priori[] / rows_theorem_permits[]` — the quantity `dim Δ̄` via the identity `k_Q = k_P − dim Δ̄`, which EXP-040 asserts on the lattice (record key `k_drop` alongside `dim_Delta_bar`, and the code raises `AssertionError` if `k_Q ≠ k_P − dim_Delta_bar`).
- Asserted at plot time: all 249 capped rows have `k_parent − k_pbb < T` (strictly below diagonal); all 30 permitted rows have `k_parent − k_pbb = T` (on diagonal); k_ceiling = k_parent − T for every row.
- Red stars: the 7 `falsification_gate.checks[]` rows at (`T`, `k_parent − k_pbb`); EXP-040 confirms `certified_reversals[].dim_Delta_bar == T` for all 7 (`equality_dim_Delta_bar_eq_T=true`).

Right panel (aggregate bars, log scale; zero bars annotated at baseline):

- Black bars "catalogue rows (279)": less = 249, equal = 30, greater = 0 — computed from the same row keys as the scatter.
- Red hatched bars "probe perturbations (26,898)": sum over `results/processed/exp040_saturation_probe.json → small_lattice.lattices[].parent_records[].perturbations.delta_bar_relation_to_T` → `less_than_T = 23,634`, `equal_T = 3,264`, `greater_than_T` absent (= 0). Cross-checked: 23,634 + 3,264 = `small_lattice.aggregate.delta_positive_perturbations` = 26,898 (asserted).

In-figure caption (exact phrase): **"upper law: zero violations in 27,177 instances"** = 279 classified catalogue rows + 26,898 probe perturbations; decomposition 249 less / 30 equal / 0 greater and 23,634 less / 3,264 equal / 0 greater printed next to it.

## File outputs

- `reports/fig_rate_distance.pdf` (58,726 bytes)
- `reports/fig_trade_law.pdf` (39,544 bytes)

Nothing else was written into the repository; the plotting driver lived outside the repo tree and installed no packages.
