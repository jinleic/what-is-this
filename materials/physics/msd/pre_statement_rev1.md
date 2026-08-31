# Gate-A-prime pre-statement: asymptotic-window check

**Revision 1; the original `pre_statement.md` and its frozen campaign remain
unchanged.** The first Gate-A grid used `p={5e-4,1e-3,3e-3,1e-2}` and produced
an accepted-shot logical-error range of roughly `0.18`--`0.93`, with fitted
slope `alpha=0.553` and `c=13.6` (`campaigns/20260829T233635Z_d9fe3c0d_1938285066db/`).
Those values are a saturation regime, not evidence about the asymptotic
quadratic coefficient. Its frozen machine verdict string remains
`REFUTE_FOR_THIS_RECONSTRUCTION`, but the scientific interpretation is
**UNTESTED**, not refuted or confirmed.

## Decisive finite experiment

Run exactly `10,000,000` shots at each of
`p ∈ {1e-4, 3e-4, 1e-3}` using the same committed Stim circuit, noise model,
detectors, and X-frame observables as `pre_statement.md`. Keep the run
single-core and low priority (`nice -n 10`), and stop after those three finite
batches. Freeze each circuit, the seed, counts, and fit in a new immutable
campaign directory; do not overwrite the original artifact.

The asymptotic prediction under test is `p_L ~= 300 p^2`, so the low-p point
has target `p_L(1e-4) ~= 3e-3`. Define its factor-of-two compatibility interval
as `[1.5e-3, 6e-3]`, equivalently `c_eff=p_L/p^2 ∈ [150,600]`. The
p=`3e-4` and p=`1e-3` rows are transition diagnostics; report their counts and
`c_eff`, but do not force a quadratic fit through a saturated point.

**Gate-A-prime reproduction-compatible:** the p=`1e-4` accepted-shot point
estimate lies in `[1.5e-3,6e-3]` (and its Wilson 95% interval is reported).
**Gate-A-prime surrogate refutation:** the p=`1e-4` point estimate lies outside
that interval after the required 10-million-shot batch. **Inconclusive:** the
accepted denominator is zero or a batch does not complete. These labels apply
to this committed Clifford/effective-output reconstruction only; they are not
universal claims about arXiv:2605.21867.

## Evidence handling

`src/run_campaign.py` is the only sampling harness and records
`NUMERICAL` rows. `src/estimate.py` reports accepted/failure counts and a
Wilson interval derived from those counts. The original Gate-A artifact and
statement remain immutable evidence for the saturation diagnosis. This
revision adds only the asymptotic-window test; it does not change the circuit,
noise injection, or decoder boundary.
