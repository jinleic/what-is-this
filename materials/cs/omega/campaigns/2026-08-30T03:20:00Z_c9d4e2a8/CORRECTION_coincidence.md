# CORRECTION (2026-08-30T04:00Z) — Lemma-1 gap was already in raw endpoint

Owner caught a suspicious near-equality: raw - published = 1.0836e-6 and
cert - raw = 1.0837e-6. Decisive test (coincidence_test.py, subprocess-isolated,
one boolean switch):

- WITHOUT Lemma-1 charged (assume dist_max exact-max-entropy):
    R_sum_low = 2.813020984577725
    omega_up  = 2.3713389007129155   → raw - published = -9.93e-8 (below pub)
- WITH Lemma-1 charged (what we shipped):
    R_sum_low = 2.8130185050939773
    omega_up  = 2.371340083602922    → raw - published = +1.0836e-6 (above pub)

DECISION: **CONSERVATIVELY DOUBLE-COUNTED.** The raw endpoint already encloses
the Lemma-1 gap; the earlier "eps_abs" added ~6.6e-11 on top of a quantity
already including a 1.08e-6 Lemma-1 effect — hence raw-cert ≈ raw-pub ≈ 1.08e-6.
Same quantity at two stages, not a coincidence.

## Owner-verified numerics at 200 bits

  cert (loose, double-counted) = 2.3713411672715115
  cert (corrected, tight)      = 2.3713400836689
  endpoint improvement          = 1.0836026e-6
  absorbed-slack TERM shrank by = 16,424.7x  (1.0836e-6 → 6.5978e-11)
  corrected cert − published   = 1.0836689e-6 (matches)
  VXXZ24 − corrected cert      = 2.1172633e-4 (margin, IMPROVED from 2.1064e-4)
  corrected cert < 2.37155181  : TRUE (improvement over VXXZ24 holds)
  corrected cert > 2.371177    : TRUE (still does NOT certify unreleased record)

## TIGHTER HONEST ENDPOINT

    omega ≤ 2.3713400836689  (raw interval endpoint 2.3713400836029
                              + 6.5978e-11 for c + ceq + Schönhage absorbed)

## Headline conclusions unchanged

- Still well below VXXZ24's 2.37155181 (margin 2.1173e-4, IMPROVED).
- Does NOT reach 2.371177 (the unreleased record; gate A stands).
- Alman25's rung genuinely improves on VXXZ24's, under rigorous arithmetic.
- Observed coincidence was not random — the same quantity at two stages
  (Lemma-1 in raw; absorbed defects on top).

## Attribution

Correction credited to owner jinleic (near-equality flag).
Decisive arithmetic: subprocess-isolated run of the same interval machinery
with a single boolean switch. Absorbed-slack term ratio 16,424.7x per owner.
