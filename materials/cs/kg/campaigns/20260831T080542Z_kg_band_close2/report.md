# Report — Gate-B band close-out stops at the mandatory anchor

## Headline

**ANCHOR FAILURE; CAMPAIGN STOPPED; FAILURE TO CERTIFY.**  The immutable
inherited `gate_B2.py` does not reproduce its recorded single-cell margin
`+0.1930`.  On the pre-registered cell and panel count, the certified envelope
margin is instead strictly negative.  No band sweep was run, no band is passed
by this campaign, and this result is not a counterexample to paper inequality
(56).

Campaign: `campaigns/20260831T080542Z_kg_band_close2/`

Pre-statement first-file commit: `90ea1ce`

## Mandatory startup verdict

The exact pre-registered invocation used:

- orthonormal-probabilist Hermites;
- even cell `(a0,a2)=(0.8,-0.6)+[-0.0001,0.0001]^2`;
- c-pair `(c_L,c_U)=(1.30,1.45)`;
- `512` panels on `[0,8]`;
- the inherited exact `Phi(p_b)-Phi(p_a)` panel weights;
- the inherited closed `[8,infinity)` tail majorant;
- `256`-bit Arb arithmetic; and
- the corrected sign test `d(c_U).lower()-U(c_L,B).upper()>0`.

The pre-registered acceptance interval was `[+0.192,+0.194]`, corresponding to
the inherited rounded claim `+0.1930` **[REPORTED from this repo, not backed by
a frozen result artifact]**.  The actual result was

`d(1.45).lower() - U(1.30,B).upper()
 = [-0.149911438966541 +/- 2.64e-16]`

**[MACHINE-VERIFIED, Arb 256 bits; radius `2.64e-16`]**.

The assertion failed twice identically; the second invocation's complete
stderr is frozen in `logs/startup_failure.log`, and the exact imported sources
and driver are frozen in `asrun/`.  Per the pre-statement, this failure stopped
all subsequent computation.  The number is a certified lower bound on the
margin of an *upper envelope*.  Its negative sign proves only that this
instrument does not certify the cell.  It does not prove `J(1.30,p)>d(1.45)`,
and still less does it refute inequality (56).

## Counterfactual controls reached before the anchor

Execution reached the anchor assertion, so every preceding startup assertion
passed:

- the paper's reproducing-kernel coefficients
  `(sqrt(2/3),-1/sqrt(3))` have squared norm `1` and `p(0)^2=3/2`
  **[MACHINE-VERIFIED in Arb against exact known values]**;
- the exact `12 x 12` even-disk cover contains the pre-anchored count `132`
  intersecting root boxes **[MACHINE-VERIFIED after Rule-17b startup
  assertion]**;
- on `[0,1]`, `phi(1) < Phi(1)-Phi(0) < phi(0)`
  **[MACHINE-VERIFIED]**, so the planted right-endpoint "supremum" weight is
  rejected and the sound left-endpoint rectangle is accepted; and
- planted above-bound and overlapping-ball margins were rejected, while a
  planted separated positive margin was accepted **[MACHINE-VERIFIED]**.

These controls establish that the zero-band outcome was produced by an
instrument capable of nonzero decisions; the campaign did not merely execute
a permanently failing branch.

## Band table at the earned strength

| c-domain | verdict in this campaign | evidence |
|---|---|---|
| `[1.0,1.3]` | **NOT RUN — anchor stop** | no band artifact |
| `[1.30,1.45]` | **FAILURE TO CERTIFY at startup cell; full band NOT RUN** | margin above, MACHINE-VERIFIED |
| `[1.45,1.75]` | **NOT RUN — anchor stop** | no band artifact |
| `[1.75,3.5]` | **NOT RUN — anchor stop** | no band artifact |
| `[3.50,4.083]` | **NOT RUN — anchor stop** | no band artifact |
| `[4.083,6.0]` | **NOT RUN — anchor stop** | no band artifact |
| razor pair `(6.0,6.05)` | **NOT RUN — anchor stop** | no artifact |

No prior PASS was imported.  In particular, a single cell never establishes a
full band over all unit cubics.

## Premises disproved before or at the anchor

1. **The recorded `+0.1930` inherited anchor is reproducible from its frozen
   source.**  Disproved by the mandatory run above.  The only occurrences of
   `0.1930` in the target folder before this campaign are prose in the prior
   pre-statement and appended README; no frozen numerical anchor artifact
   exists **[source audit + MACHINE-VERIFIED rerun]**.
2. **Exact CDF panel masses are a new tightening to adopt.**  Disproved by
   source: inherited `cdf_weight` already computes `Phi(p_b)-Phi(p_a)` and
   `cell_envelope` already uses it **[source-read]**.
3. **`phi(p_b)` is a Gaussian-density supremum on a positive panel.**
   Disproved: Gaussian density decreases for `s>=0`, and the `[0,1]` planted
   control verifies the strict reverse ordering **[INFERENCE from derivative;
   MACHINE-VERIFIED control]**.  Using `phi(p_b)(p_b-p_a)` as an upper weight
   would be permissively unsound.
4. **Subdivision at `m q=t^2` is needed to unlock the second hinge.**
   Disproved at the semantic level before computation.  Lemma D.4 directly
   gives `f(x,y;t)<=f(E_P,W_P;t)` for the independently valid upper bounds
   `E_P,W_P`; no product lower bound is needed.  The inherited `A2-r2` is not
   a lower bound on `||e|^2-|o|^2|` because both terms are independent upper
   envelopes.  The inherited positive branch already equals the paper's
   direct monotone envelope, while its fallback only drops `-t` and is
   coarser **[INFERENCE, derived from the paper's stated formula and source]**.
5. **The prior README's `[3.50,4.083]` summary matches its frozen as-run JSON.**
   Disproved by source artifacts: `b5_result.json` has zero certified c-tiles
   and unresolved margins ranging from about `-0.0283` to `-0.0194`, not seven
   of eight and `-2e-5` to `-5e-5` **[source-read]**.

None of these is a disagreement with the paper.  They are failures of this
repo's inherited instrumentation/provenance.

## Rule-7 actual scope sentence

This stopped campaign swept exactly one even-coefficient cell,
`(a0,a2)=(0.8,-0.6)+[-0.0001,0.0001]^2`, with all compatible odd coefficients
covered by `sqrt(1-dist(0,B)^2)`, at the single c-pair `(1.30,1.45)`, using
`512` panels on `[0,8]`, the closed `[8,infinity)` tail, and `256`-bit Arb; it
also ran only the stated basis, `132`-root-cover, Gaussian-weight, and margin
counterfactual controls, and it did **not** search any complete c-band, any
other coefficient cell or radius, the razor pair `(6.0,6.05)`, c in
`[0.993405,1.0)` or above `1.45`, degrees above three, non-unit cubics, the
paper's C1/C3/splice certificates, Gate-C septic/Krivine families, or any
construction outside the paper's unit-cubic threshold family.

## Refutation gate

Not triggered.  No explicit unit cubic received a certified lower bound on
`J` above an upper bound on `d`.  Therefore nothing was escalated to Main as a
paper refutation.

## Named next campaign and cost

**Named next campaign: `Gate-B anchor provenance reconstruction and direct
Lemma-D.4 restart`.**  It should first reconstruct, from git history and exact
as-run sources if they exist, how the stale `+0.1930` was obtained; then start
a new pre-registered campaign whose anchor is a known value generated by an
independent analytic unit-cubic calculation, not an undocumented envelope
number, and whose panel bound is the direct `f(E_P,W_P;t)` formula with strict
endpoint comparisons.  Estimated cost: **less than 1 CPU-hour for provenance
and anchor reconstruction, then 2–6 CPU-hours for the six full bands if the
new anchor passes [INFERENCE from the frozen predecessor's 231-second partial
run and the fixed per-band budgets].**
