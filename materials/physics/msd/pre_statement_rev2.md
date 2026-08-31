# Gate-A-prime pre-statement: corrected asymptotic quantification

**Revision 2; append-only correction.** The original `pre_statement.md` and
`pre_statement_rev1.md` are frozen and unchanged. Revision 1 used the band
`p_L ∈ [0.0015, 0.006]` at `p=1e-4` because of a steering arithmetic error.
The implementation agent caught that error before any promotion leaf used the
band. The measured point `p_L=0.038356` lies outside both Revision 1's band
and the corrected band below, so the Gate-A-prime verdict is insensitive to
the mistake.

## Correct implied-pL table

The paper's reported Fig. 10 headline is `p_L ~= c p^2` with `[REPORTED]`
`c ~= 300` (arXiv:2605.21867). The implied values are:

| physical `p` | `300 p^2` target | factor-of-two target band |
| ---: | ---: | ---: |
| `1e-4` | `3e-6` | `[1.5e-6, 6e-6]` |
| `3e-4` | `2.7e-5` | `[1.35e-5, 5.4e-5]` |
| `1e-3` | `3e-4` | `[1.5e-4, 6e-4]` |

These are arithmetic consequences of the reported quadratic, not new paper
measurements. At `p=1e-4`, `10^7` shots and an approximately 0.9 acceptance
fraction would therefore provide order tens of accepted failures, not
thousands.

## Frozen experiment and exact criterion

Gate-A-prime ran exactly `10,000,000` shots at each
`p ∈ {1e-4, 3e-4, 1e-3}` using the unchanged Clifford/effective-output Stim
circuit, noise model, detector postselection, and X-frame observables from the
original statement. The frozen artifact is
`campaigns/20260829T234855Z_02ebe2db_40d024611d98/`.

**Gate-A-prime reproduction-compatible:** the `p=1e-4` accepted-shot point
estimate lies in the corrected band `[1.5e-6, 6e-6]`. **Gate-A-prime surrogate
refutation:** that point estimate lies outside the corrected band after the
required 10-million-shot batch. **Inconclusive:** the accepted denominator is
zero or a batch does not complete. These labels apply only to this committed
surrogate; they are not universal claims about the paper.

The measured rows are `[NUMERICAL]`:

* `p=1e-4`: 9,206,721 accepted, 353,133 accepted failures,
  `p_L=0.0383560`, `c_eff=3.8356e6`, Wilson 95% interval
  `[0.0382321, 0.0384803]`.
* `p=3e-4`: 7,808,694 accepted, 861,556 accepted failures,
  `p_L=0.110333`, `c_eff=1.2259e6`.
* `p=1e-3`: 4,430,706 accepted, 1,415,900 accepted failures,
  `p_L=0.319565`, `c_eff=3.1957e5`.

The falling `c_eff` ratio (`3.8e6 -> 1.2e6 -> 3.2e5`) and the log-log fit
exponent `alpha=0.9201` show floor-like/saturated behavior rather than the
quadratic asymptotic law at every reachable p in this surrogate. The measured
`p=1e-4` point is about `1.28e4` times the corrected `3e-6` target, and the
other two rows are about `4.09e3` and `1.07e3` times their corrected targets.
The campaign summary's verdict string remains
`SURROGATE_REFUTATION_AT_LOW_P`.

## Changelog

* **Revision 2:** corrected `300 p^2` arithmetic after the implementation
  agent caught Main's steering-arithmetic error; no frozen artifact or prior
  statement was overwritten. The verdict is unchanged because `0.038356`
  fails both the erroneous and corrected bands.
* **Revision 1:** asymptotic-window experiment statement, retained verbatim in
  `pre_statement_rev1.md` for provenance.
