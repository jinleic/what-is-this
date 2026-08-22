# Manuscript claim bench — summary

Generated 2026-08-18 from `results/manuscript_claim_bench.json` (the machine-readable bench, 49 rows).
Paper audited: `/Users/jinleic/jinleic-workspace/math/qec/reports/paper_pbb_nogo.md` (571 lines; 11-page PDF).
Method: single process; every anchor read from disk this session; no SAT/Stim/pytest-suite runs (test names verified to exist by reading `tests/`, not executing).

**Status counts: 49 rows = 44 VERIFIED + 3 UNDERWRITTEN + 2 OVERWRITTEN + 0 UNSOURCED.**
Zero UNSOURCED for the two boxed results, Theorems G/H/I, and the six key claim families (W4 [5,12]; EXP-042/4,621; 61/133/134; 7/7 reversals; 368 rows; remaining-witness numbers).

## Two headline checks — reviewer re-run commands

1. **W4 certified `d_DEM_mech >= 5` on both circuits** (anchor `results/processed/exp040_circuit_distance.json`):
```bash
cd /Users/jinleic/jinleic-workspace/math/qec && .venv/bin/python -c "
import json; d=json.load(open('results/processed/exp040_circuit_distance.json'))
print(d['verdict'])
print('gross', d['circuits']['gross']['lb']['value'], d['circuits']['gross']['ub']['value'])
print('pbb  ', d['circuits']['pbb']['lb']['value'], d['circuits']['pbb']['ub']['value'])"
# expect: {'certified_d_dem_mech_both_ge_5': True, 'summary': 'certified d_DEM_mech >= 5 on both circuits'}
#         gross 5 12 / pbb 5 12
```
Verdict this session: **True** — paper §6(iii) says only "intervals … are tied" → **UNDERWRITTEN** (anchor supports `[5,12]` on both).

2. **EXP-042 `SMALL_CLASS_EMPTY`, 4,621 parents, `T < k_P/2` bucket empty** (anchor `results/processed/exp042_residual_probe.json`):
```bash
cd /Users/jinleic/jinleic-workspace/math/qec && .venv/bin/python -c "
import json; e=json.load(open('results/processed/exp042_residual_probe.json'))
print(e['verdict']); print(e['aggregate'])"
# expect: SMALL_CLASS_EMPTY / {... 'distinct_parent_fingerprints': 4621, 'bucket_T_lt_half_k': 0, ...}
```
Verdict this session: **SMALL_CLASS_EMPTY confirmed** — paper predates it (silent) → **UNDERWRITTEN**.

## Rows (compact)

| id | § | claim (short) | status | anchor |
|---|---|---|---|---|
| B01 | Abstr/S3 | Boxed 1: Theorem H no-go bound | VERIFIED | partial_runs/exp039_nogo_module.json + tests/test_pbb_nogo.py |
| B02 | Abstr/S3.2 | Boxed 2: Theorem I forced saturation | VERIFIED | tests/test_pbb_theorems.py |
| B03 | Abstr/S3.2 | dim Δ̄ ≤ k_P/2 via left kernel (k_P/2 on every BB parent) | VERIFIED | tests/test_pbb_theorems.py |
| B04 | Abstr/S1 | 368 rows / 202 parents; [[144,12,12]]; [[360,12,≤24]] | VERIFIED | processed/exp027_delta_audit.json |
| B05 | Abstr/S5 | 134/202 exact; 61 family-closed incl. Gross T=12 & 11 [[144,12,12]] parents | VERIFIED | exp039_nogo_module.json + parent_aebb649586c1f5ab*.json |
| B06 | Abstr/S5 | 249/368 capped a priori; 30 permitted | VERIFIED | exp039_nogo_module.json (counts) |
| B07 | Abstr/S5.2 | 7/7 reversals slack 0, hypotheses replayed | VERIFIED | exp039_nogo_module.json (falsification_gate) |
| B08 | Abstr/S7 | Probe 64 parents / 26,898 / 496 on law, 0 above | VERIFIED | processed/exp040_saturation_probe.json |
| B09 | S1 | Validity = symmetry of M (not vanishing) | VERIFIED | src/qec_research/symplectic/core.py |
| B10 | S2 | G(i) centralizer invariance | VERIFIED | tests/test_pbb_survival.py |
| B11 | S2 | G(ii) dimension identity + J.3 sector identities | VERIFIED | notes/theorem_j_xsector.md + tests/test_pbb_survival.py |
| B12 | S2 | G(iii) survival criterion | VERIFIED | tests/test_pbb_survival.py |
| B13 | S2 | Xcen 40→12 on phase2_58; shortcut false | VERIFIED | notes/theorem_j_xsector.md |
| B14 | S3 | Lemma 1 + Cor 2 (submodule, orbit absorption) | VERIFIED | tests/test_pbb_nogo.py |
| B15 | S3.1 | Exact T via UNSAT; dual GF(2) witness re-verification | VERIFIED | partial_runs/exp039/ (134 certs) + exp041_t_crosscheck.json |
| B16 | S3.1b | Symmetry break sound; 30 SAT/UNSAT agreements | VERIFIED | tests/test_translation_symmetry_break.py |
| B17 | S3.1b | 6.2× speedup (77.8→12.5 s) | VERIFIED | notes/failed_routes.md:776 (FR-022) |
| B18 | S3.2 | Lemma 2: k_P = 2 dim L on 202 parents | VERIFIED | tests/test_pbb_theorems.py |
| B19 | S3.2 | Lemma 3 ceiling on 368 rows | VERIFIED | tests/test_pbb_theorems.py |
| B20 | S3.2 | 133/134 scope; exception 9a7638586033 T=4, k_P=12, sandwich 4–6 | VERIFIED | exp039 certs + exp042 known_exception_context |
| B21 | S4 | Gross row + Corollary 3 family closure | VERIFIED | parent_aebb649586c1f5ab*.json (T=12=k_P) |
| B22 | S4 | 11 [[144,12,12]] parents T=12; 220.4 s / 3 witnesses | VERIFIED | parent_4c7eb964a6a3e38c*.json (220.43 s) |
| B23 | S4 | Genuinely non-CSS: 19-parity + one-hot CP-SAT INFEASIBLE | VERIFIED | certificates/exp034_target_full_lc_decision.json |
| B24 | S4 | Exact d=12 two-sided (≤5 MIM; 6: 72/72; w12 witness) | VERIFIED | certificates/pbb_12_6_0193_distance.json |
| B25 | S4 | Sector UNSAT seconds 502/674/422/576 | VERIFIED | partial_runs/pbb_12_6_0193_distance_sectors/sector_0{0..3}_pysat.json (501.98/673.75/421.97/575.90, sum 2173.60) |
| B26 | S5 | Totals 202/134/61/249/30 | VERIFIED | exp039_nogo_module.json (counts) |
| B27 | S5 | Per-length table incl. capped col + ≤144 subtotal 117/117/57/254/224 | VERIFIED | exp039 by_n + R5 recompute |
| B28 | S5 | 68 uncertified = n∈{180,360}; slowest 16,493 s | VERIFIED | exp039 by_n + certs (max 16492.95 s) |
| B29 | S5.1 | Cross-tab 279 rows; cells 109/67/1, 11/10, 25/13, 13/30 | VERIFIED | exp039 rows_* + tests/test_paper_claims.py |
| B30 | S5.2 | Saturation table + executable gate | VERIFIED | exp039 falsification_gate (7× slack 0) |
| B31 | S6 | 5-of-7 parent-profitable; 1.78 vs 2.00; ≤1.13× | VERIFIED | exp036_envelope_check.json + R7 recompute |
| B32 | S6 | 7/7 CSS-dominated; 162 candidates; factors 3.0–11.1× | VERIFIED | processed/exp036_envelope_check.json |
| B33 | S6 | X-sector: J.0/J.1; open claim 3.16M / 12 lattices / 368 rows, 0 violations | VERIFIED | notes/theorem_j_xsector.md + xsector artifacts |
| B34 | S6 | 60 pending / 35 candidate escapees | **OVERWRITTEN** | exp037_envelope_classification.json now {60, 30, 278} |
| B35 | S6iii | DEM intervals "tied"; 6–17× slower (10.3/16.6/6.3) | **UNDERWRITTEN** | exp040_circuit_distance.json + technical_report.md §5–6 |
| B36 | S6 | Heuristic bound ≥2.5× loose (40→16) | VERIFIED | partial_runs/exp036/row_0001.json |
| B37 | S7 | 279-row upper law, spread to −24 | VERIFIED | exp039 rows_* + R17 recompute |
| B38 | S7 | Probe 23,634 / 3,264 / 0; containment verified | VERIFIED | exp040_saturation_probe.json |
| B39 | S7 | 7 reversals at dim=T=k_P/2, M̄=Δ̄ | VERIFIED | exp040_saturation_probe.json (certified_reversals) |
| B40 | S7 | Conj C k-halving confinement, zero exceptions | VERIFIED | exp039 + tests |
| B41 | S7 | Conj D uniform +2 gap | VERIFIED | exp039 gate checks (R7) |
| B42 | S8 | 15/17/4 checks; 159 passing + 1 skipped | VERIFIED | tests/ (claims-matrix collect audit 160) |
| B43 | S8 | 5 reversals, 41 of 87; 155-row δ>0 audit | VERIFIED | exp036_delta_closure.json + exp027_delta_audit.json |
| B44 | S8 | CNF SHA-256 hash-binding + forgery regressions | VERIFIED | exp039 certs + forgery tests |
| B45 | S8 | Environment line (Python 3.13.14, numpy 2.5.2, …) | **OVERWRITTEN** | live venvs + exp040 env block: 3.13.9 / 2.4.6 |
| B46 | S8 | Reproduce commands exist | VERIFIED | experiments/exp039_nogo_module.py |
| B47 | S10 | Conclusion restates 61/202, 249/368, 7× slack 0 | VERIFIED | inherited anchors |
| P1 | post-paper | W4 certified d_DEM_mech ≥ 5 both circuits, ub 12 → [5,12] | **UNDERWRITTEN** | processed/exp040_circuit_distance.json (verdict true) |
| P2 | post-paper | EXP-042 SMALL_CLASS_EMPTY: 4,621 parents, T<k_P/2 bucket = 0 | **UNDERWRITTEN** | processed/exp042_residual_probe.json |

All artifact paths are under `/Users/jinleic/jinleic-workspace/math/qec/` (absolute paths in the JSON).

## Drift findings (all live-anchor quotes in the JSON)

1. **B34** — paper "35 candidate escapees"; live `exp037_envelope_classification.json` counts `{PENDING: 60, CANDIDATE_ESCAPEE: 30, DOMINATED: 278}` (was 35/273 at the 2026-08-17 claims-matrix read; sweep advanced). PENDING 60 still matches.
2. **B45** — paper env "Python 3.13.14 … numpy 2.5.2 in qec-codesign/.venv"; both live venvs and the newest artifact env blocks record **Python 3.13.9 / numpy 2.4.6** (Darwin 25.5.0, stim 1.16.0, CaDiCaL 1.9.5 all match).
3. **P1** — paper says intervals "tied"; anchor certifies the stronger `[5,12]`-on-both (quote the numbers in a revision).
4. **P2** — paper silent on EXP-042; Conjecture B′ residual class now known empty for all ℓm ≤ 24.

*(Correction 2026-08-18: an earlier draft flagged B25 as OVERWRITTEN — retracted. The per-sector wall times are persisted in `results/partial_runs/pbb_12_6_0193_distance_sectors/sector_*_pysat.json` and match the paper to rounding.)*

## Claim-echo audit (grep-level)

`5,12`, `[5,12]`, `4,621` (and `4{,}621`): **0 occurrences** — nothing to defend; both are post-paper facts anchored at P1/P2.
`61` (4×, lines 40/383/505/562), `134` (6×, lines 39/49/251/303/331/449), `368` (8×, lines 16/42/242/297/305/420/537/563), `249` (3×, lines 42/505/563): every occurrence defensible — `parents_family_closed=61`, `parents_T_exact=134` (133/134 scope recomputed), `catalogue_rows=368` (exp027 summary), `catalogue_rows_capped_a_priori=249` (exp039 counts).
`24` (3×): line 17 `[[360,12,≤24]]` catalogue-advertised bound (tests/test_conventions.py::test_bb_reference_parameters); line 318 n=180 capped count 24 (exp039 rows_capped_a_priori by-n); line 458 spread minimum −24 (R17 recompute over the 279 classified rows).

## Scope exclusion

**EXP-043 is out of scope** per parent steering (2026-08-18): `results/processed/exp043_sector_entanglement.json` does not exist — the streamed result was lost pre-broker-restart and no ledger row claims it; no anchor was invented. A grep for `ENTANGLEMENT|entangle` over `reports/paper_pbb_nogo.{md,tex}` returns no matches, so no paper phrase requires it.
