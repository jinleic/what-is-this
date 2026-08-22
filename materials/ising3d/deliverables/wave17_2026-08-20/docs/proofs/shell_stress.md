# Shell-route stress test: the uniform shell defect is not exhibited by any exact box *(wave 16, lead takeover)*

**Question.** Theorem PM's mildest upper-endpoint condition, PM-S (`proofs/paired_momentum.md`
sec. B2), certifies `K_c <= K*` from a **uniform per-mode relative defect**
`Ghat_L(k) <= (1 - eps)/(2K* lambda(k))` on all modes `lambda(k) >= Lambda`, with certified
thresholds `eps*(1/4, 2) = 0.02290...` and `eps*(63/250, 2) = 0.006130...`. Do the exact
finite boxes exhibit defects of the required size? If yes, the wave-17 dual-walk global
audit front was to be launched.

**Answer [COMPUTATION].** No. On every exact box on record (`2x2x2`, `2x2x3`, `2x2x4`,
`3x3x2`; full 2^N enumerations, exact Fraction intervals from
`results/bounds/paired_momentum.json`):

| box | min shell defect at `K*=1/4`, `Lambda=2` |
|---|---|
| 2x2x2 | -0.182744 |
| 2x2x3 | -0.141396 |
| 2x2x4 | -0.179753 |
| 3x3x2 | -0.020384 |

All are **negative** (`Ghat` *exceeds* the `1/(2K* lambda)` form on those shell modes), versus
the required `+0.0229`. At the `K* = 63/250 ~ 0.252` target the requirement is `+0.00613`;
the only positive shell minima on record are at the incumbent endpoint `I3/2`: `+0.00520`
on `3x3x2` (and negative on every smaller box), with the K-trend going *down* toward `-0.020`
at `K*=1/4`.

**Structural aggravation.** The required `eps*(1/4, 2) = 2.29%` exceeds the defect of the
comparison Gaussian field itself: `Ghat_sat(k) = 1/(I3 lambda(k))` gives a uniform defect
`1 - 2K*/I3 = 0.01080... < 0.0229`. Any field that saturates the certified two-point
comparison class therefore cannot satisfy PM-S at `Lambda=2, K*=1/4` — the route needs the
true critical `Ghat` to sit uniformly 2.2x *below* the saturant's defect, and the measured
defects point the other way.

**Verdict.** FALSIFIED-for-purpose: the dual-walk global audit front has no empirical or
structural support at `Lambda=2` and is NOT launched. This is a statement about the
premise's measurability, not about Theorem PM: PM-S with `Lambda = 3, 4` (
`eps*(1/4,3) = 4.40%`, `eps*(1/4,4) = 11.10%`-class) or genuinely new input remains an exact
open condition, and the massive-Watson bracket `mu*(1/4)^2 in (1/1024, 1/32)` is untouched.

**Provenance.** Lead takeover of the wave-16 ShellStress front after its agent stalled
(zero files on disk). All numbers recomputed from the certified intervals in
`results/bounds/paired_momentum.json` with exact `Fraction` arithmetic
(`results/shell_stress.json`, 250 mode rows); the stale agent probe
("R_2 = 0.194 at 2x2x4") measured a different, non-operative quantity (mass fraction, not
the uniform per-mode defect) and is not used. Independent replay:
`tests/test_shell_stress.py`.

## PM-S at Lambda = 3, 4 (wave 17, lead direct derivation)

Artifact: `results/shell_stress_lambda34.json` (exact Fraction derivation on certified
stored inputs only; no new graph data). Formula (certified): `eps*(K*,Lam) = delta*(K*) *
(2K*) / J_lo(Lam)` with `J_lo(3) = 0.1241778...`, `J_lo(4) = 0.0438377...` (same 96^3 grid),
`delta*(K*) = I3/(2K*) - 1` (stored exact).

| K* | eps*(.,3) | eps*(.,4) | factor vs band-average (3 / 4) |
|---|---|---|---|
| 63/250 | 1.18% | 3.34% | 4.06x / 11.50x |
| 1/4 | 4.40% | 12.46% | 4.03x / 11.41x |
| 49/200 | 12.45% | 35.27% | 3.95x / 11.18x |
| 6/25 | 20.50% | 58.08% | 3.87x / 10.95x |

The concentration factor `2K*/J_lo(Lam)` is (nearly) K*-independent: PM-S(Lambda) demands
the high-lambda shell carry a uniform per-mode defect 4.0x / 11.5x the whole-band average.
Measured on all four exact 3D boxes x four benchmark-free evaluation points, premise
decided per shell by the WEAKEST mode's defect-interval upper end (the uniform requirement):
all 48 (box, evaluation, Lambda in {2,3,4}) combinations are falsified — the weakest mode is
negative on 45 of 48, and the tightest combination (3x3x2, incumbent) ends at +0.005196
against `eps*(63/250,2) = 0.613%` and `eps*(63/250,3) = 1.177%`. While individual HIGH-lambda
modes do carry large positive defects (up to +0.384 on 3x3x2), the premise needs EVERY mode
above threshold, and each shell contains a negative-defect mode. The required excess is not
approached on any box, at any evaluation point, at any shell.

**Verdict [COMPUTATION]: FALSIFIED-for-purpose at Lambda = 2, 3, 4.** The uniform-shell
route family is closed. The box-measured defect concentrates in the INFRARED, the direction
opposite to what PM-S needs, and the audited LP-class floor `D(1/L)` forces any uniform
infrared floor to vanish as L grows — the surviving candidate shape is a **lambda-dependent
weighted premise** `DeltaHat >= <eps(lambda) w(lambda)>` (recorded as the PM-W design item;
requires a genuinely new two-point input, not a re-weighting). Theorem PM remains exact and
conditional; the upper endpoint 0.25273100... and its what-it-would-take statement are
unchanged. Ledger H417.
