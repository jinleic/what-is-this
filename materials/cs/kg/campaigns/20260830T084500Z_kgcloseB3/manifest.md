# manifest.md — campaign 20260830T084500Z_kgcloseB3 (agent KgGateBClose)

## Scope
Close the OPEN gate-B sub-bands of cs/kg via the paper's D.3 2D-Christoffel
route (arXiv:2608.11158, scratch/paper_full.txt lines 1961-1984, read
first-hand; inequality implemented = line 1984 pairing "Ubar(box, cL) <=
d(cU)_lo certifies J(c,p) <= d(c) for all c in [cL, cU]").

## Files
- `pre_statement.md` — pre-registered budget/box/outcomes; committed BEFORE
  the first computation (git d3ebb62); Addendum 1 pre-registered the
  adaptive-c-splitting after float calibration (git dad0eaf). Later, mid-run,
  TWO further pre-run-of-that-code amendments were committed to the file as
  ADDENDUM 2 (ordered-pair q-scheme) — flagged honestly: Addendum 2 was
  written when the 2x fold defect was diagnosed and BEFORE the b6/b7 runs
  that use it, but AFTER the b5 v2 run. Reader should treat band b5's
  q-scheme margins as ADDENDUM-2-era code (same pre-registration caveat).
- `code/core2.py` — certified Arb machinery (GL32 construction byte-audited
  port of the predecessor's core.py; psi/K_e/K_o per paper lines 1976-1980;
  tail_account = paper lines 1982-84 majorant in exact erfc/phi balls).
- `code/gate_B2.py` — final envelope + search (v2 scheme); See header
  docstring for the certified object: ordered-pair M/q pipeline with
  certified case split on (mq - t^2)_+.
- `code/gate_B2_v1_superseded.py` — v1, kept for provenance; defective
  (double-counted panels, no phi weight), NEVER used for any verdict.
- `logs/b5_progress.json`, `logs/b5_result.json` — frozen run output
  (band [3.5, 4.083], 8 geometric c-tiles).
- `sha256.txt` — checksums.

## NOT finished (budget cutoff at the 200-request soft bound)
Bands [1.0, 1.3], [1.45, 1.75], [1.75, 3.5], [6.0, 12.0] NOT run to verdict.
The anchor re-run and b5/b6 certify attempts consumed budget; the run is
frozen mid-campaign per the repo's checkpoint discipline. Verdicts below are
the honest per-band state.

## Post-run probe (bg_9, arrived at wrap-up) — corner-cell truth level

Float probe (COMPUTATIONAL-EVIDENCE, never a verdict): max J(3.75, p) over the
b5/b6 corner region a0 ∈ [0.833, 1.0] × a2 ∈ [−0.333, −0.167] = **0.14467**
(at (a0,a2) ≈ (0.9429, −0.3330), odd part ≈ 0.0114 nearly vanishing), while
d(3.75) = 0.15684. True spare ≈ **+0.0122** at the band midpoint. This pins
the stuck-cell story: certified envelope margin there was ≈ −2e-5, so the
envelope excess ≈ 0.012 — the fold-inflation + panel case-split granularity
account for essentially ALL of it, and no cell is anywhere near a violation
of the paper inequality (true slack is ~600× the paper's own typical margin).
