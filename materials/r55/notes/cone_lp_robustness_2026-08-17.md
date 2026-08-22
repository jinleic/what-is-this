# Cone LP-Robustness Audit (m=3 campaign, post-Task 3)

Date: 2026-08-17. Method: fully independent re-implementation (own `Fraction`
verifiers + own exact Bland-rule simplex with exact dual certification; no
producer imports). Artifact under audit: `r55/data/higher_identity_m3.json`
(SHA-256 `8295ab4d...`, disposition `M3_NO_CUT_IN_FROZEN_CONE`).

## 1. Independent verification: PASS, zero mismatches

- 3,215 states: sorted/unique, `deficiency == a+b`,
  `excess_balance == 2(a+b) - budget2(d)` with `budget2 = {20:12, 21:15, 22:16,
  23:15, 24:12}` recomputed from `structural_tables.json`; all g windows valid.
- All three route certificates hold at every state in exact arithmetic; minimum
  slack 0 on each (alpha provably minimal for the stored multipliers);
  `45*alpha` equals every stored bound.
- All three primal witnesses: weights sum 45, balance exactly 0,
  `g_lo <= 0 <= g_hi` in aggregate, values equal stored.
- Statuses, unavailable route-4 record, and disposition re-derived identically.
- All 58 input hashes/sizes re-verified.
- **LP optimality bonus**: the stored bounds (360, 14310/349, 45) are certified
  true LP optima of the frozen cone by an independent exact simplex; strong
  duality holds on all three routes (route 3 has an alternative optimal dual
  `(1, 1/16, 0, 0)` of equal value).

## 2. Slack structure (why the routes are robust)

| Route | Exact bound | Edge | Tight states | What a flip requires |
|---|---:|---|---|---|
| total_deficiency | 360 | <=315 | 729 (all d=22) | remove 2,093 states of demand>7 (=>270) |
| degree20_count | 14310/349 (~41.00) | <1 | 9 | remove all 561 d=20 states (=>0) |
| deficiency_ge8_count | 45 | <1 | 3,035 | remove all 3,035 def>=8 states |

Demand classes (route 1): `(a+b)+e/2 ∈ {8 (729 states), 15/2 (1364), 6 (1122)}`.

## 3. Cone sign structure

g windows: straddling 2,049; one-sided positive 1,055; one-sided negative 108;
width-0 exactly 8. Knowing the exact g sign would prune 1,163 states (36.2%)
but moves **only** route 2: 14310/349 → 585/19 (drop g<0) or 12195/302
(drop g>0). Routes 1 and 3 are unmoved even by exact g.

Structural certificates: d=20-only and d=24-only cones are LP-infeasible
(forced aggregate `g_hi <= -9630` resp. `g_lo >= +9630`); d=21-only or
d=23-only def>=8 bound drops to 675/16; d=22-only stays 45.

## 4. Ranked next-invariant effects (measured by exact LP re-solve)

1. m=4-moment-style added equality: route 2 41.00 → 36.97 (window halved) →
   35.89 (quarter) → 34.59 (exact midpoint); routes 1/3 zero reduction.
2. Sharper q envelopes (graphicality-realized): same channel as 1; route 1's
   obstruction sits on census-EXACT classes, so envelopes cannot move it at all.
3. Exact g-sign: route 2 → 30.79; routes 1/3 unchanged.
4. Parity/integrality: zero LP effect.
5. Degree-class exclusion (non-cone theorem, e.g. "no d=20 vertex exists in
   R(5,5,45)"): the **only** measured flip channel (route 2 → 0 < 1).

## 5. Consequences recorded for the next campaign

- The m=3 cone result is not a search artifact: it is the exact LP optimum of
  the frozen relaxation, certified twice by disjoint machinery.
- Any m=4 expectation must be calibrated: modeling an added global equality as
  g-window shrink moves route 2 only, and by <10 units. The genuinely new
  content of the m=4 row (per-state t, diamond, K3+K1, i4(dual) intervals) is
  not captured by that proxy and must be computed exactly before any claim.
- The decisive flip channels are degree-class exclusions; these cannot come
  from the cone itself (d=20-only cones are already infeasible) — they require
  catalog/gluing-level mathematics, exactly the AM46 bottleneck the m=4 plan's
  motivation section records.

Independent reproduction: the campaign final reviewer re-derived the §2
tight-state decomposition (729/1364/1122) and the budget2 table from the
artifact's own state rows on 2026-08-17, and hand-verified all three
certificates and witnesses; see `.superpowers/sdd/task-4-report.md`.
