# Pre-statement — m=12 alpha-guard campaign (bounded cell, t=64), mceliece target

**Committed 2026-09-01T03:28Z, BEFORE any m=12 computation.** Fixes the
instance, the instrument, the evidence fork, the anchors, the
counterfactual plants, the pass/fail criterion, and the measured abort
gate. Implements — and where infeasible, explicitly amends, with the
reason declared here — ADDENDUM 2 of `cs/mceliece/pre_statement.md`,
which is quoted verbatim below (primary first-hand; not paraphrased).

## 0. ADDENDUM 2 verbatim (from cs/mceliece/pre_statement.md as on disk 2026-09-01)

Clause (iii), verbatim:

> ## (iii) alpha stays MEASURED on the certified ladder; split by row
>
> alpha (the differential identity Pi F' + Pi' F = G^2 F^(2)) is the
> INDEPENDENT CHECK ON APON LEMMA 3 ITSELF — it is what would catch the
> lemma being wrong, and it is the evidence behind the headline phrase
> "confirmed by measurement".  Therefore:
>
>   * m <= 11 (the 13 certified instances): alpha remains a per-instance
>     MEASUREMENT, unchanged.  It is already done and affordable at these
>     sizes.  Its label does not change.
>   * m = 12, IF AND ONLY IF measuring alpha there proves unaffordable:
>     alpha may be DERIVED from the exact delta plus Apon Lemma 3, and that
>     row's alpha is then labelled CITED-DEPENDENCY.  This is a different
>     KIND of evidence, not a different confidence level, and the artifact
>     must state it in one sentence: at m <= 11 the mechanism is MEASURED;
>     at m = 12 it is INFERRED from a cited lemma given a measured delta.
>
> No global relabelling of alpha.  A reader must be able to see which rows
> test Apon and which assume him.

Clause (i), verbatim (the delta route this campaign inherits):

> ## (i) delta by the exact Lagrange route, with ABORT on degree failure
>
> F is constructed as f_j = SUM_i Y[j,i] * lam_i^{-1} * L_i, with
> L_i = Pi/(Z - a_i) / Pi'(a_i).  Off-diagonal vanishing L_i(a_l) = 0 (l != i)
> holds BY CONSTRUCTION: Pi/(Z-a_i) retains the factor (Z - a_l), and the
> division is exact (zero remainder, asserted per i).  Therefore
>
>     lam_l * f_j(a_l) = Y[j,l]   for ALL (j,l)
>
> follows EXACTLY from the n scalar checks L_i(a_i) = 1, at cost O(n D)
> instead of O(k n D).  This replaces both the earlier full k x n grid check
> (m = 6,7,8) and the row-wise pivot+8-sample check (m = 9 onward).
>
> LOAD-BEARING DEGREE CONDITION, checked per instance, ABORT (not warn) on
> failure: max_j deg(f_j) <= D.  This is what licenses the Lemma-2 degree
> argument (both sides of Apon (12) have degree <= n-1: LHS terms
> Pi/(Z-a_i) have degree n-1; RHS G^2 f_j has degree <= 2t + D = n-1).  The
> bound is essential, not decorative: f and f + Pi agree at every support
> point yet differ, so without deg <= D the agreement-at-n-points argument
> collapses.  A run whose degree check fails produces NO delta verdict.

Clause (ii) verbatim (anchor discipline inherited here):

> ## (ii) uniform re-run across ALL THIRTEEN instances, ONE code path
>
> The delta guard is re-run for every one of the 13 certified instances
> (m = 6,7,8,9,10,11) through a single implementation, so the ladder rests
> on one standard verified by one code path.  [...] If any of m = 6,7,8
> fails to reproduce under the new route, that is a FINDING: escalate to
> Main BEFORE writing anything.

Clause (iv) verbatim: "The m <= 11 labels, thresholds, conventions, flag
depth R = 4, pass/fail criteria and the enumerated claim of record are
unchanged by this addendum.  No universal ('m <= N complete') phrasing is
used for any N."

## 0a. Registered instance and population (fixed before compute — rule 16)

Exactly ONE m=12 instance is registered, the shipped `mceliece348864`
lattice row:

    (m, n, t, seed) = (12, 3488, 64, 16384)

- k = n - mt = 2720, D = n - 2t - 1 = 3359, 2t+3 = 131, |E| = 4096,
  candidate-irreducible population at t=64 is monic degree-64 polys over
  F_4096 (count (4096^64 - 4096^32)/64 for the exhausted set — NOT
  swept; G is one seeded irreducible, scalar-Rabin-confirmed).
- seed 16384 is fixed BY THIS SENTENCE, before compute. Build via the
  UNCHANGED `src/instance.py` seed stream. If any build guard fails,
  regenerate with seed 16385, 16386, ... ; every substitution is logged
  (pre_statement.md section-6 rule, unchanged). Guards are never adjusted.
- NOT registered, hence never computed in this campaign: t = 48, t = 96,
  any other t, any n != 3488, any other support ordering, any m > 12.
  They remain UNREACHED and are named in the rule-7 boundary (section 6).

## 0b. Flag-collection domain (fixed before compute)

The identity check is coefficient-wise mod Pi. Equivalent checks at the
support points are licensed by Apon's Lemma-2 degree argument (both sides
degree <= n-1; deg condition ABORT-gated). Domains, all PREFIXES of the
single seed-16384 support ordering, chosen now:

    FLAGCOL = {64, 128, 256, 512, 1024, 2048, 3488}

3488 = n, the full support. These are the ONLY domains of this campaign;
no adaptive extension, no re-run at other c.

## 1. The statutory fork, decided by MEASUREMENT, with the pre-declared amendment path

ADDENDUM 2(iii) forks: measure alpha at m=12 (label MEASURED) if and only
if affordable; else derive it from exact delta + Apon Lemma 3 (label
CITED-DEPENDENCY). This campaign executes the fork in THIS order:

  1. Attempt the measurement path on the registered instance with the
     unmodified instrument (`instance.py` guards, alpha block), without
     commitment — a measured slowness, not an assumption.
  2. Record the instrumented rate via in-process timing only
     (`time.process_time()` / `time.perf_counter()`; per Main's steering,
     `ps`/CPU% readings are NOT usable for rate claims on this box).
  3. If the FULL measurement would exceed the registered CPU budget
     (section 4), the fork resolves to the DERIVE clause, with the
     measured numbers as its justification. The resolution is then
     RECORDED, with the measured rate, as an explicit amendment to this
     pre-statement in the campaign dir BEFORE any derived verdict is
     computed (section 5 of the process rules: declared before compute,
     reason attached). Under no circumstances may the fork be resolved
     after seeing any alpha verdict number.

Budget set NOW, before calibration: 3 CPU-hours total (10,800 s
process_time) for the m=12 verdict block — construction + ramp calibration
+ alpha verdict + delta verdict. Exceeding this = unaffordable. The
calibration loop measures the per-full-check cost of the unmodified
instrument at m=12 scale first (one full check = one call of
`_guards`-route inTimingWrapper); its measured CPU seconds per full check
and the projected wall to complete 3488 of them are what resolve 3.
Also measured for context: at m=11 (2048,48,6211, k=1520, D=1951 = dim
3451 coeffs), the same instrument completed the full 3488-equivalent grid
(that's 2912-equivalent there) at 285.91s (frozen record), but because
m=12 has 2x the k and 1.7x the D and a larger n, the scaling is
super-linear in numpy gather cardinality (gather over `(k x nz)` blocks
per point) — NOT assumed linear in (n k D); it is measured.

## 1a. THE DERIVE INSTRUMENT, if reached: exact delta + Apon Lemma 3

Per clause (iii), the m=12 alpha is then: "INFERRED from a cited lemma
given a measured delta". Concretely, on the registered instance:

- Delta is measured EXACT per clause (i)'s chain (5 parts, section 3):
  lam_l * F(a_l) = Y[:,l] for ALL (j,l), exact by the Lagrange route.
- The cited lemma: Apon 2026/1810, Lemma 3: Pi F' + Pi' F = kappa * G^2
  F^(2) with kappa = 1 on the Gate-A construction (kappa = 1 fixed by the
  multipliers lam_i = G(a_i)^2 / Pi'(a_i), Apon eq. (1); square-free-ness
  of G holds because G is irreducible, so G, G' coprime — machine-checked:
  gcd(G, G') constant, recorded).
- The inference is F_2-EXACT (no sampling): the identity is a polynomial
  identity of degree <= n-1 on both sides (Lemma-2 degree argument);
  given it holds for the hidden F (Lemma 3), the exact delta chain
  guarantees our constructed F IS the hidden F of the construction, hence
  the identity holds for OUR F. Label: CITED-DEPENDENCY, statement
  recorded verbatim as the addendum requires: "at m <= 11 the mechanism
  is MEASURED; at m = 12 it is INFERRED from a cited lemma given a
  measured delta."
- MANDATORY SEPARATE MACHINE CHECKS at m=12 (cheap, run regardless):
  (a) gcd(G, G') is a nonzero constant (square-free-ness, the lemma's
  hypothesis); (b) deg(Pi) = 3488, Pi has n DISTINCT roots (support
  distinctness — machine-checked by the build's distinct-support assert);
  (c) max_j deg f_j <= D (the ABORT gate); (d) lam_i = G(a_i)^2/Pi'(a_i)
  recomputed and byte-equal to the build's lam vector.

## 2. Anchors before any m=12 number is reported (HARD GATE)

Re-run through THIS campaign's code path (frozen 14-09Z code, current src,
which differ ONLY by the inert `forced_G` ctor param — diff frozen
byte-for-byte in the campaign's evidence/) at least TWO of the m=6..11
instances, requiring identical outcomes vs the frozen 14-09Z records:

    (11, 2048, 48, 6211)  — the largest frozen row  (frozen 285.91 s)
    (10, 1024, 40, 5113)  — large frozen row
    (6,   64,  3, 1387)   — the smallest frozen row (frozen 0.03 s)

Required: every field of the new run's delta/alpha/beta/eps/gamma record
must equal the frozen record's corresponding field, EXCEPT `elapsed_s`
(timing is machine-dependent and NOT part of the claim) and any field
whose value is derived from `elapsed_s`. Any other field differing: NO
m=12 verdict is written; escalate to Main with the mismatch. The anchor
rows keep their m<=11 MEASURED alpha labels — no relabel, no downgrade,
no upgrade.

## 3. Counterfactual plants (HARD GATE, rule 14: the instrument must be able to ESTABLISH, not merely agree)

Before any m=12 number, plant corrupted instances and show rejection:

  - CF-1 (delta plant): rebuild the m=11 anchor with the SAME seed but a
    support shuffle that duplicates one element (well-defined: duplicate
    a_12 := a_11 in the list, reducing distinct count by 1) — the
    distinct-support assert inside Instance.build must fire (recorded as
    the expected REJECTION), OR, if that assert is not reached, the
    lagrange-unit check must fail with >= 1 failing i.
  - CF-2 (alpha plant): construct an m=11 instance whose F is replaced
    by F + Pi*Z^0 (adding the full-support vanishing polynomial; changes
    NOTHING at support points — a pure degree-violation of max_j deg f_j
    <= D). Required: the degree gate flags deg_ok = false and NO delta
    verdict is produced (ABORT semantics demonstrated); separately, this
    corrupted F must FAIL the alpha identity (the corrupted F has deg
    D+1, so the Lemma-2 argument does not apply; verify by direct check
    that the identity FAILS on at least one support point).
    NOTE: this plant demonstrates the instrument REJECTS a corrupted
    alpha; it cannot demonstrate the instrument REJECTS an alpha-failing
    instance that is not degree-degenerate — that direction is licensed
    by the identity being machine-checked at all n points (exhaustive,
    not sampled).

Both plants record their exact rejection certificate (which check, which
index) in the results.

## 4. Budget and measured abort gate

- Instrumentation: `time.process_time()` and `time.perf_counter()` inside
  the program; `nice -n 10`; BLAS/OMP threads = 1 (pre_statement.md
  section 4, unchanged). ps/CPU% NOT used for any claim (Main steering
  2026-09-01: ps-accounting on framework Python children is broken).
- Calibration run FIRST (before the commit of the verdict-critical
  budget-relevant amendment IF 1.3 resolves to derive): measure (a) build
  CPU+wall for (12,3488,64,16384); (b) ONE full alpha-check pass at
  m=12-scale matmul-dims (mock k x n grid, mock D) — measured CPU s per
  grid pass — then extrapolate to 3488 support points; (c) the m=11
  anchor (11,2048,48,6211) re-timed in-process as the rate sanity anchor.
- HARD ABORT GATE (registered now): if the extrapolated m=12
  full-alpha-measurement exceeds 3 CPU-hours (10,800 s process_time
  projection from measured rate), the FORK RESOLVES TO DERIVE with the
  measured numbers; the amendment is recorded in this campaign dir (a
  file `AMENDMENT_fork_resolution.md`, dated and checksummed) BEFORE the
  derived verdict is computed. If the extrapolation instead fits the
  budget, the path stays MEASUREMENT and no amendment is filed.
- TOTAL campaign wall budget: 4 wall-hours. If the total run exceeds 4
  wall-hours, freeze whatever sub-block has completed, checksum, and
  report a bounded verdict over what completed (rule 16: bounded
  sub-domain, never extrapolated).

## 5. Verdict shape and evidence labels (fixed now)

- delta verdict at (12,3488,64,16384): PASS iff the 5-part chain (exact
  division, n unit checks, 200 off-diagonal probes, 3 assembly checks,
  degree gate) returns delta_exact true. FAIL is a FINDING (escalate to
  Main before external writing). ABORT on degree failure = no verdict
  (recorded as such).
- alpha verdict at m=12: EITHER "MEASURED at full 3488-point grid
  (addendum path 1)" OR "CITED-DEPENDENCY: inferred from measured exact
  delta + Apon Lemma 3 (addendum path 2)" — decided by the section-1
  fork, resolution recorded (as an amendment if resolved) BEFORE compute
  of the derived verdict. Never a mixed or unlabeled kind. The one-
  sentence artifact sentence (addendum-mandated) is included verbatim in
  the results.
- Anchor rows keep m<=11 labels verbatim (MEASURED alpha, etc.).
- Boundaries for the claim of record (rule 7): see the scope sentence
  (section 6). In particular: "m=12 complete" is NOT claimed; the
  section-3.6 hole remains BOUNDED, never "closed"; no NIST-cell
  verification is claimed; no Asymptotic or cryptographic-scale claim is
  made.

## 6. Scope sentence (rule 7, exact)

Swept in this campaign: the single instance (12, 3488, 64) at seed 16384;
the flag-collection prefix chain FLAGCOL = {64, 128, 256, 512, 1024,
2048, 3488} of its single seeded support ordering; the 5-part exact
delta chain of ADDENDUM 2(i); the alpha instrument (unmodified 14-09Z
grid route, rate-measured; BCOMP only if the derive fork is reached and
only as the exact certificate supporting the CITED-DEPENDENCY
classification); the 3 registered anchor rows in byte-equality mode; 2
counterfactual plants. NOT swept: any t != 64 at m = 12 (t = 48, t = 96
remain UNREACHED), any n != 3488, any support ordering other than the
seed-16384 shuffle, any m > 12, any census row (N_fam / per-point rank /
held-count ladder) at m = 12, any sub-support construct, any
non-irreducible G, and any claim venue beyond one yielded verdict.
A single m=12 instance verdict is a BOUNDED FINITE STATEMENT about that
instance; it is not an m=12 completion, not a NIST-cell verification,
and does not touch the bounded (not closed) Apon section-3.6 hole.

## 7. Named next action (pre-declared)

If the m=12 verdict lands PASS: the named next action is a second m=12
row at t = 96 (the README ladder's probe-above row), under a fresh
pre-statement, budgeted from THIS campaign's measured m=12 rates. If the
verdict lands FAIL: escalate to Main before any external write. If the
fork resolves to derive: the named next action is to re-attempt
MEASURED alpha at m=12 with a vectorized-(over-points) rewrite of the
alpha grid in a successor campaign, not to relabel this one.

— MceliecelM12, 2026-09-01T03:28Z, committed BEFORE any m=12 compute.
