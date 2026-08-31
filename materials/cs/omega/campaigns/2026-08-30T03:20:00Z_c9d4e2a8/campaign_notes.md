# Stage (b) rung 1 — certified enclosure result

Run: 2026-08-30T03:20Z. Machine: M3 Ultra, python 3.14.3, python-flint 0.9.0.
Single process. Reproduces from command:
  python src/stage_b_rung1_certified.py
(with src/interval_core.py, src/alman25_float.py alongside, plus
scratch/alman_code/data/W1.00_2.371339.mat for input verification.)

## Certified statement (Addendum A3, option (a) — SLACK ABSORPTION)

At the exact released float64 parameter point (24 855 params, sha256
f2336913…b46f) of Alman et al. 2025 (SODA'25, arXiv:2404.16349, osf.io/mw5ak,
data/W1.00_2.371339.mat), the combination-loss program's certified enclosure
is:

    omega <= 2.371340083602922  (raw interval endpoint)
    omega <= 2.3713411672715115 (with feasibility slack absorbed, explicit)

with the slack decomposition:
  - Lemma-1 max-entropy residual of the SHIPPED Lagrange multipliers
    (dominant): eps_L1 = 2.062e-6 (regions 1-6 individually:
    1.33e-6, 1.05e-6, 7.39e-7, 6.28e-7, 2.06e-6, 1.64e-6);
  - stage-a inequality violation c_max = 1.137e-10 (absorbed / M = 5.4e-11);
  - stage-a equality violation ceq_max = 2.390e-11 (absorbed / M = 1.1e-11);
  - Schonhage-line half-ulp float noise -6.2e-15 (absorbed / M = 3.0e-15).

## Interpretation (per Main directive)

Even under the rigorous arithmetic treatment (Lemma-1 slack + outward-rounded
entropies + absorbed feasibility defects), the certified endpoint
2.3713412 sits WELL BELOW VXXZ24's 2.37155181 — a margin of 2.1e-4. So:

  **Alman25's rung omega <= 2.371339 genuinely improves on VXXZ24's
  omega <= 2.371552, and the improvement survives rigorous interval
  arithmetic with explicit feasibility slack.**

This is the first rigorous interval statement about any rung of the
combination-loss ladder. It does NOT certify the unpublished 2.371177 record
(that artifact remains unreleased — gate A finding, 2026-08-29).

## Evidence chain

- Stage (a) transcription: parameter count 24855 == 24855; float64 value
  reproduces the published 2.371339 to 9.95e-8 (display rounding);
  max c violation 1.137e-10, max ceq 2.390e-11 — inside their refine
  tolerance 1.1e-9. Sources read first-hand; verbatim quotes in campaign dir.
- Per-block containment (Addendum A3 R3): 18/18 num_block entropy blocks,
  6/6 hash penalty blocks, 5/5 Lagrange ceq blocks, 1/1 value line.
- Taint test (Addendum A3 R2, sweep form): every value-chain input responds;
  dist/dist_max respond through the entropy path into num_block (width
  3.09e-8 under 2^-30 taint), confirming the constraint-block interval path
  is live at every stage.
- Width-nonzero invariant: num_block widths 1.81e-20 (genuine, nonzero);
  Schonhage-line c_viol[-1] width 6.42e-15 (genuine, nonzero, equal to the
  float64 noise captured by outward rounding).
- Root cause fixed along the way: python-flint's arb(lo, hi) is (mid, RADIUS)
  — broadcast to Main and siblings (kg/, delcap/, rs-pe3d/).

## What is NOT claimed

- That omega <= 2.371177 (the unreleased record) is certified. Gate A
  finding stands: not independently checkable from arXiv:2608.16884 as
  published.
- That eps_abs is minimal. It is the certified slack achieved with the
  SHIPPED Lagrange multipliers producing Lemma-1 residuals of order 2e-6;
  tighter multipliers could shrink eps, but that would require optimizing
  them (not done here).
- That 2.3713412 is a *tighter* bound than 2.371339. It is
  LOOSER by eps_abs: the honest cost of absorbing the infeasibility of the
  released float64 point under exact arithmetic. This is the correct
  semantics per Addendum A3 — never hide the slack.

## Prior-rung comparison (rigorously certified)

  rung 2.371552 (VXXZ24)  <— RUNG 1 (this result, 2.3713412)  <— published
                            (interpolates)                       2.371339
  Improvement of rung 1 over prior: 2.11e-4, survives rigorous interval
  arithmetic with all feasibility defects folded in.
