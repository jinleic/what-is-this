# Quantified obstruction — tile 3 of band [4.083, 6.0] (Campaign B, Phase-1 FAIL)

Frozen 2026-09-01 by `KgCampB`. Phase 1 (pre-registered mq = t² panel
subdivision, validated on the known-failing cell) **FAILED its fixed
threshold**: tile 3 does not certify; the frontier cell's margin ball is
unchanged from Campaign A's. Phase 2 does NOT run (pre-statement §3/§5).
This file quantifies the obstruction: what binds, what each mechanism
recovers, and what no knob of the A/B instrument family can recover.

## Verdict chain (MACHINE-VERIFIED, Arb 256-bit outward unless labelled)

1. Startup controls 16/16 PASS (`logs/startup_controls.json`), including the
   kernel-admission assert with the cap-at-1 counterfactual (deficit
   0.224744871391589049...e+0 certified) and both permissiveness
   counterfactuals.
2. Phase-1A frontier cell under V
   (`logs/phase1_frontier_cell.json`): OPEN at panel cap 131072, margin
   −1.76291787028287693244667476546e-5 ± 4.20e-35 — bit-identical recovery
   vs A of −4.19e-35 (i.e. exactly zero to 35 digits).
3. Phase-1B full tile 3 under V (`logs/phase1_tile3.json`): OPEN at the SAME
   cell/depth/panels/margin as Campaign A (open_depth 14, panel cap reached,
   95/251 leaves certified open-stop at 163 cap-131072 evaluations out of
   819; A's 7649-eval stop included the same 251 boxes).

## The decomposition (`logs/decomposition_tile3.json`, rule-15 balls)

| entry | value (Arb ball) | label |
|---|---|---|
| D1 A envelope, cap 131072, margin | −1.76291787028287693244667476546e-5 ± 4.20e-35 | reconciliation: reproduces A's frozen margin (`reproduces_A_frozen: true`) |
| D2 V (mq = t² kill + subdivision), cap 131072 | margin identical to D1 to all displayed digits; recovery **exactly 0** | the Phase-1 mechanism is DEGENERATE on this cell |
| D3 panel-sup 2^17→2^19 | gain 2.20922696204705e-5 ± 2.83e-20 (margin → +4.4630909e-6) | COMPUTATIONAL-EVIDENCE for extrapolation only; the gains themselves are certified |
| D3 panel-sup 2^19→2^21 | gain 5.52295164708941e-6 ± 4.54e-21 (margin → +9.9860426e-6) | ditto |
| D3 panel-sup total 2^17→2^21 | 2.76152212675599e-5 ± 2.29e-20 | does NOT close the cell: the DEFICIT stays −1.7629e-5 at cap by definition — wait, corrected below |
| D4 disk-feasible |e|-sup (diagnostic, firewall) | even-hinge recovery 4.94055685681106e-4 ± 2.25e-19; margin would be +4.76426506978277e-4 | the |e|-hinge coefficient-space excess is 28× the deficit |
| D5 tail | 3.6353172540086755822e-13 ± 3.81e-33 | identical in every variant; differential exactly 0 |
| D6 within-panel skew | 5.52295164708941e-6 (2^19→2^21) | contained in D3 |

**Correction to the D3 row as first written (kept for honesty):** the panel
margins at 2^19/2^21 ARE positive (+4.46e-6, +9.99e-6) — the BINDING panel-sup
deficit is thus NOT panel-count-limited beyond 2^19 for this cell: A's
pre-registered ladder stopped at 131072 = 2^17, and a cap of 2^19 already
turns this leaf's margin positive. The cell is closed by panel count ≥ 2^19
in the A-envelope family. (See "implied knobs", below.)

## Why the pre-registered mechanism recovers zero (the arithmetic)

The mq = t² crossing set for this cell: odd2·K_o(s) = (c_L·s)² reduces to
odd2·(5/2 − s² + s⁴/6) = c_L², whose positive root is
**s\* = 45.703407790989350730** (certified bracket; root squarely outside
[0, S_MAX] = [0, 8]). Below √3 the left side is bounded by (5/2)·odd2 =
5.93e-5 ≪ c_L² = 17.21, so no second root exists in range. CONSEQUENCE: at
every panel of the ladder the case-boundary inside the panel never fires —
the "M-hinge-only" branch is certifiable on a zero-measure set only, and the
certified second-hinge kill (max(E², W²) ≤ t_L·E_min at outward precision)
fires on no panel of the rung where A stopped ( Grade: all kills occur only
on panels where A's bound is already 0; V ≡ U to ±2.4e-73 displayed).
The fold-pair geometry on this ring cell: |e(s)| ranges over BOTH signs of a
wide even interval whose modulus near the fold is ~1.19–1.22 while
W = √odd2·√K_o ≈ 2.3e-3 ≈ 0 — the pair {|p(s)|, |p(−s)|} genuinely has both
members ≈ |e| ≫ t at small s, so the second folded hinge is LIVE for real
(s not just bound-form), and no s-side subdivision changes the sup-sum.

## Where the actual obstruction sits (the tight statement)

The instrument's per-panel slack on this cell is the fold-pair sum being
sup-based over a WIDE even box whose |e|-sup is carried by the box corner
(a0_lo, a2_hi) — a corner sitting at ring distance ≥ 1 from the feasible set:
the disk-feasible part of B (a 3-point candidate set, certified) would drop
the even sup by 4.94055685681106e-4 = 28.0× the entire deficit
(`D4_..._recovery`), turning the frontier margin to +4.76426506978277e-4.
BUT enforcing B ∩ disk in the |e|-sup is a COEFFICIENT-SPACE mechanism — the
paper's own envelope (and A's, and this campaign's by pre-registration)
intentionally covers the FULL box against all compatible odd parts via the
odd radius r_o = √(1−dist(0,B)²); for a box straddling the disk boundary the
odd-radius mechanism requires the FULL box's |e|-sup (a point with even part
outside the disk can still have unit-norm (e, o) combinations), so the D4
route is NOT free: it must be re-derived jointly with the odd part, not
substituted (firewall held; this is exactly the named next campaign).

## Implied knobs to close (the numbers an operator needs)

- **Panels:** ladder cap 2^19 (524, 288) already certifies this leaf at
  margin +4.46309091764170e-6 (D3, certified single-cell run at 2^19: 68 s
  CPU). Below 2^19 the ladder cannot: 2^17 (A's cap) is short by 1.7629e-5.
- **Depth/box floor:** irrelevant here — the frontier leaf is NOT at floor
  (depth 14 of 40; can_split true; A's stop was "panel cap reached", not
  depth).
- **The mq = t² crossing subdivision:** zero effect on this cell (crossing at
  s\* = 45.70 ∉ [0,8]); its domain of usefulness is mid-band cells where the
  fold-pair crossing sits inside the s-range — NOT this ring geometry.
- **Coefficient-space feasibility (re-derived jointly with the odd part):**
  the only mechanism measured here with 28× the needed leverage. This is
  Campaign C's first change, priced: the frontier leaf needs ONE certified
  jointly-feasible (even, odd) envelope tightening; the 124 stacked boxes at
  A's stop (tile-3 DFS remainder) were never evaluated, so the tile's full
  close needs that mechanism plus a resumed sweep — NOT re-run from scratch.

## Scope (rule 7)

This obstruction quantifies exactly: the A frontier cell at the A c-pair
(4.1485893120, 4.181778026496), under A's envelope family (with and without
the pre-registered mq = t² subdivision) and the D4 diagnostic variant, on
[0,8] panels plus the closed tail at 256-bit outward Arb. NOT covered: any
other leaf's closure status under caps ≥ 2^19 (only this single cell was
run at 2^19/2^21), tiles 1-2/4-49 under any mechanism, any coefficient-space
envelope taken to verdict, and all Campaign C domains.

## Named next action

Campaign C (separate pre-statement, per assignment scope): (a) joint
(even, odd) disk-feasibility envelope — re-derive r_o against the odd part
CONDITIONED on the even part being disk-feasible — targeting the ring cells;
(b) balloon the panel cap to 2^19 (cheap: the frontier cell alone at 2^19 is
68 s CPU); then re-run tile 3 to close, and only then resume the band sweep
tiles 4-49 under the validated mechanism. Cost estimate for (a)+(b) on tile
3: ~2-4 CPU-hours (computed from the 251-box / 819-eval A stop pattern +
this cell's 68 s at 2^19).
