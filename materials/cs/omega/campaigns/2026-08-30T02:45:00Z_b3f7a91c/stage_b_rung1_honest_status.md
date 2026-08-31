# Stage (b) rung 1 — HONEST STATUS (2026-08-30T02:50Z)

## What is verified

Per-block containment (Addendum A3 R3) at the released point
W1.00_2.371339.mat (sha256 f23369136314cc51497d0c468c00d69a02b74d2b665073b9f934240ab428b46f):
- num_block entropies (all 18 region×dim blocks): interval encloses float64.
- hash penalties H(dist_max) - H(dist) (6 regions): encloses float64.
- Lagrange ceq residuals: 5 blocks checked, enclosing float64.
- Schonhage value line: interval encloses the float64 value.

Lemma-1 max-entropy slack of the SHIPPED Lagrange multipliers (upper endpoints,
region 1..6): 1.326e-6, 1.046e-6, 7.392e-7, 6.283e-7, 2.062e-6, 1.636e-6.
max = 2.062e-6 — this dominates eps_abs as Main predicted.

## What is NOT yet certified

The final value interval collapsed to a POINT at exactly the float64
7.783640596221247. A genuine outward-rounded enclosure of a transcendental
expression at a non-special rational point cannot have zero width. Therefore
the current interval machinery has a remaining leak (some sum/radius path
collapses to float64) and the omega_cert_upper = 2.37133890054342 it outputs
is NOT a certified statement — it is the float value with interval labels.

Per pre_statement Addendum A2/A3 and Main's standing directive: this is
TRANSCRIPTION-stage evidence, not a bound statement. Nothing is recorded as a
bound. Escalation to Main filed before any conclusion is recorded.

## Next engineering steps (well-defined, in priority order)
1. find the collapse: every G operation must propagate ball radius; the
   suspect is entropy_vec() returning a point (uses `float(...)` casts) and
   the term_frac/region_prop item() multiplications along the hashing lines;
2. recompute R_sum from entropies with the outward buffer at every endpoint;
3. assert R_sum_low > 0 AND < float estimate (impossibility checks);
4. report Omega_cert_upper + eps = max(lemma1_eps_max, absorbing the stage-a
   defects c=1.137e-10, ceq=2.390e-11, schonhage=-6.2e-15).
