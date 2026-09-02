# Pre-statement — second m=12 row at t=96: alpha-first guard campaign (bounded cell 12/3488/96), mceliece target

**Committed 2026-09-01, BEFORE any t=96 computation.** Fixes the exact
`(m,n,t,seed)` cell, the instrument (INHERITED UNCHANGED from the certified
`M12ALPHA_DIVONLY` campaign — no second convention is created), the alpha
mode and its MEASUREMENT-GATE fork, the budget, the anchors, the plants,
and the verdict logic. Implements the NAMED NEXT ACTION of
`campaigns/2026-09-01T03-28-10Z_M12ALPHA_DIVONLY/` (t=96 row, fresh
pre-statement, budgeted from measured rates) and the README Gate-A ladder
row (12, 3488, 96 | k = 3488 − 1152 = 2336 > 0, 2t+3 = 195, D = 3488 − 2·96 − 1 = 3295).

## 0. Registered instance and population (fixed before compute — rule 16)

Exactly ONE m=12, t=96 instance is registered — the README ladder's
probe-above-ship row:

    (m, n, t, seed) = (12, 3488, 96, 16384)

- k = n − mt = 3488 − 1152 = 2336; D = n − 2t − 1 = 3488 − 193 = 3295;
  2t+3 = 195; |E| = 4096; support n = 3488 (as in every certified row).
- seed 16384 is fixed BY THIS SENTENCE, before compute. Build via the
  UNCHANGED `src/instance.py` seed stream (same `random.Random(seed)`
  stream as all 13 certified ladder rows and the t=64 row). If any build
  guard fails, regenerate with seed 16385, 16386, …; every substitution
  is logged here in an explicit dated amendment BEFORE the substitute
  verdict is computed (pre_statement.md section-6 rule, unchanged).
  Guards are never adjusted.
- NOT registered, hence never computed in this campaign: t = 48, any
  other t ≠ 96 at m = 12, any n ≠ 3488, any other support ordering, any
  m > 12, any census row (N_fam / held-count ladder) at any m. They
  remain UNREACHED and are named in the rule-7 boundary (section 6).

## 1. Instrument: INHERITED, byte-identical; no second convention

The instrument is the certified `M12ALPHA_DIVONLY` code path, copied
byte-for-byte into `code/` at init time (before any compute):

- `instance.py`, `fastfield.py`, `gfield.py`, `census.py` — identical to
  the current `src/` (which is itself byte-identical to both the
  `M12ALPHA_DIVONLY` snapshot and — except for instance.py's inert
  `forced_G` ctor parameter, an unused Gate-B tooling artifact — the
  frozen `2026-08-30T14-09Z_57200ADD` snapshot; byte-diffs recorded at
  init and re-verified at freeze).
- Guard chain semantics: the ADDENDUM-2(i) 5-part exact delta chain
  (Lagrange route with ABORT on degree failure), the alpha differential
  identity ΠF′ + Π′F = G²F⁽²⁾ MEASURED at support points by exact
  uint16 log-table evaluation (the identical operation sequence used at
  m = 6..11 and at the t=64 row), guards β/γ/ε (Wronskian pair, k
  identity, gcd/max-degree), gcd(G,G′) square-freeness, λ recomputation
  equality, and the F₂-linearity control on seeded selectors — ALL as
  fielded, unchanged. Any deviation discovered mid-run is a FINDING and
  escalates to Main before any verdict is written.

The delta/alpha verdict stage is a mechanical re-instantiation of the
t=64 verdict stage at the new registered cell (same checks, same order,
same abort semantics, new instance tuple). The stage script is itself
derived by minimal edit from
`campaigns/2026-09-01T03-28-10Z_M12ALPHA_DIVONLY/verdict_m12_stage.py`
and copied into `code/` before runs.

## 2. Alpha mode: MEASURED FIRST, fork only on a FAILED measurement gate

ADDENDUM-2(iii) (`cs/mceliece/pre_statement.md`, quoted verbatim in the
prior campaign pre-statement) forks alpha at m=12 between MEASURED
(identical instrument over the support grid) and DERIVED/CITED-DEPENDENCY
(exact delta + Apon Lemma 3) — measure iff measuring is affordable.

This campaign registers the measurement path as PRIMARY and declares,
before any compute:

1. ALPHA IS MEASURED FIRST. The identical m≤11 instrument's
   per-support-point alpha check runs on the FULL support grid
   (n = 3488 points) on the registered t=96 instance.
2. The MEASUREMENT GATE: measured t=64 calibration (t=64 row: 0.1504
   CPU s per support point, 700.53 CPU s for the full 3488-point grid,
   measured in-process; build 2175.43 CPU s) projects the t=96 alpha
   grid at ≈ same n with mildly larger D (3295 vs 3359 — the t=64 D is
   LARGER, so per-point vector widths at t=96 are ≤ t=64's). The
   registered budget below is set with ≥ 3x headroom over that
   projection. The per-point cost is re-measured in-process here
   (calibration checkpoint written BEFORE the full-grid commit point)
   and the full-grid CPU total is measured, not assumed.
3. DERIVE-branch entry condition (registered NOW, evaluated only
   against CPU totals MEASURED IN THIS PROCESS, never estimated): if
   the cumulative `time.process_time()` at the alpha-grid checkpoint
   has already exceeded 75% of the registered budget
   (0.75 × 32400 = 24,300 CPU s), the alpha grid is aborted and the
   DERIVE branch executes: alpha labelled CITED-DEPENDENCY, inferred
   from the exact 5-part delta chain + Apon Lemma 3 (κ = 1 on the
   Gate-A construction), with the mandatory separate machine checks
   (a)–(d) of the prior pre-statement section 1a (gcd(G,G′) constant;
   deg Π = n distinct roots; max_j deg f_j ≤ D ABORT gate; λ
   recomputed byte-equal), the one-sentence addendum artifact sentence
   recorded verbatim, and the resolution (measured numbers + abort
   point) committed as a dated checksummed `AMENDMENT_fork_resolution.md`
   BEFORE the delta chain is computed. Under no circumstances may the
   fork be resolved after seeing any alpha verdict number beyond the
   aborting checkpoint itself.
4. Expected resolution (projection, NOT a decision input): measured
   rates say MEASURED with ≥ 3x headroom; the DERIVE branch exists so
   the gate is honest, and it is entered ONLY on the measured numbers.

The one-sentence kind labeling, verbatim, either way: "at m ≤ 11 the
mechanism is MEASURED; at m = 12 [this row] it is MEASURED / it is
INFERRED from a cited lemma given a measured delta."

## 3. Anchors (HARD GATE, before any t=96 number is reported)

At least THREE anchors are re-run through THIS campaign's code snapshot
(`code/`, the inherited instrument) and compared field-by-field against
against the same frozen records used by the t=64 campaign
(`campaigns/2026-08-30T14-09Z_57200ADD/delta_uniform_13instances.json`):

    (11, 2048, 48, 6211)   — largest frozen row  (t=64 campaign anchor)
    (10, 1024, 40, 5113)   — large frozen row    (t=64 campaign anchor)
    (6,   64,  3, 1387)    — smallest frozen row (t=64 campaign anchor)

Required: every field of each new run's record byte-equal to the frozen
record's corresponding field, EXCEPT `elapsed_s` (machine-dependent
timing, not part of the claim) — the same exclusion the t=64 campaign
used — and any additive field this run records on top (disclosed per
anchor, not part of the frozen set). Any other field differing: NO
t=96 verdict is written; escalate to Main with the mismatch. PASS
requires all three anchors byte-equal (the comparison harness records
`byte_equal_vs_frozen=true` per anchor, field-by-field). Anchor rows
keep their m≤11 MEASURED alpha labels — no relabel, no downgrade, no
upgrade.

## 4. Counterfactual plants (HARD GATE, rule 14: the instrument must be able to ESTABLISH, not merely agree)

At least THREE adversarial plants, REJECTED by the same inherited
verifier (never by an outside ad-hoc checker), mirroring the t=64
campaign's three — re-instantiated at the m=11 anchor scale
((11,2048,48,6211), as the t=64 plants were) EXCEPT where the t=96 cell
is the natural base, each recording its exact rejection certificate.
The t=64 suite's amendment history is inherited verbatim: CF-2's
value-level invisibility prediction (f+Π passes all value checks; the
degree-gate ABORT is the rejection) and CF-3's scalar-multiply form
(row-XOR proven inert a priori in F₂ — never re-registered):

  - P-1 (duplicate support): rebuild the m=11 anchor instance with a
    support shuffle that duplicates one element (a_12 := a_11).
    Expected rejection: Π′(a)=0 at the duplicated position(s) — the
    build's distinct-support assert fires; failure mode recorded as
    the certificate.
  - P-2 (degree gate): replace F by F + Π on the m=11 anchor (changes
    nothing at support points; pure degree violation deg = D+1 > D).
    Expected: `deg_ok = false`, ABORT semantics — NO delta verdict is
    produced; separately, and recorded as the AMENDMENT-2-visible
    prediction, the value-level alpha check is INVISIBLE to this
    corruption (0 failures on the probed points).
  - P-3 (scalar multiply): row 0 of the m=11 anchor replaced by c·f_0
    with c = 3 ∈ F_4096 \ F₂ (3² = 5 ≠ 3, verified in-run). Expected:
    alpha value check FAILS on row 0 wherever the value relation is
    exercised; the c² ≠ c certificate (c = 3, c² = 5) recorded. The
    simultaneous delta-binary trip is recorded as part of the
    certificate if it fires (per AMENDMENT 3's note), but the
    alpha-level rejection is the plant's purpose.
  - P-4 (t=96-native plant, additional insurance at the registered
    cell): on the t=64-registered m=12 anchor… NOT APPLICABLE — no
    m=12 t=64 rerun is registered; instead the P-3 mechanism is
    cross-checked at the t=96 anchor base itself with a SMALL probe
    grid (first 128 support points only, bounded): corrupting the
    built instance's F row 0 by ×3 must produce ≥ 1 value-check
    failure on the probed points while the UNCORRUPTED instance is
    clean on the same points. This is an instrument-direction check
    at scale, not a second convention: same verifier, same ops.

If any expected rejection fails to materialize as registered: that is
a FINDING, escalated to Main BEFORE any verdict is written.

## 5. Budget (registered before compute)

- Instrumentation: `time.process_time()` / `time.perf_counter()` inside
  the program only; `nice -n 10`; BLAS/OMP threads = 1.
  `ps`/CPU%-based rates are NOT used for any claim (Main steering,
  inherited).
- Registered HARD budget: **9 CPU-hours = 32,400 CPU s `process_time`,
  total for the whole verdict block** (build + alpha grid + delta
  chain + controls). This is set at 3× the t=64 campaign's measured
  total (2893.37 CPU s) + margin for the larger D/k — the projection
  says t=96 costs ≈ t=64 (D is smaller; k, n unchanged), but the cap
  is deliberately generous while the abort checks above keep it
  honest. Exceeding the cap = the run freezes whatever sub-block has
  completed, checksums it, and reports a bounded verdict over what
  completed (rule 16: bounded sub-domain, never extrapolated).
- Wall budget: 6 wall-hours total. Same bounded-report rule on breach.
- The MEASUREMENT GATE (section 2.3) is judged on
  `time.process_time()` measured inside the process, at the alpha
  checkpoint.

## 6. Verdict shape, pass criterion, and evidence labels (fixed now)

The campaign terminates in exactly one of:

- **PASS (FROZEN-CERTIFIED, bounded to the cell):** on the registered
  instance (12, 3488, 96, 16384),
  (i) build guards all pass (G monic irreducible degree 96,
  scalar-Rabin-confirmed; support distinct; G nonzero on support);
  (ii) delta: the 5-part ADDENDUM-2(i) chain returns delta_exact =
  true — exact division Π/(Z−a_i) with zero remainder on all 3488,
  n unit checks L_i(a_i) = 1 (3488/3488), 200 off-diagonal probes with
  0 violations, 3 assembly spot-checks, degree gate max_j deg f_j ≤ D
  (= 3295) passing (ABORT, no verdict, if it fails);
  (iii) alpha: MEASURED, identity holds at ALL 3488 support points,
  0 failures (or, only if the section-2.3 gate fires,
  CITED-DEPENDENCY via the declared derive branch);
  (iv) controls: gcd(G,G′) constant, λ recomputed byte-equal,
  F₂-linearity on seeded selectors;
  (v) anchors: 3/3 byte-equal (section 3);
  (vi) plants: 3/3 REJECTED with recorded certificates (section 4);
  (vii) all within the registered budget (section 5).
  Runner exit code 0. Verdict written as `verdict.json` +
  `plants_result.json` + `anchor_*.json` inside the campaign dir,
  checksummed (freeze), inventory row appended (close).
- **FAIL:** any guard, chain part, anchor comparison, or plant
  expectation fails as registered. A FAIL is a FINDING: escalate to
  Main via hub BEFORE external writing; campaign freezes whatever
  completed and closes FROZEN-NEGATIVE (registered verdict logic) or
  FROZEN-INCONCLUSIVE (instrument/integrity failure only).
- **ABORT/BUDGET:** degree-gate failure = NO delta verdict (recorded
  as such); budget breach = freeze finished sub-blocks, close
  FROZEN-INCONCLUSIVE over the completed bounded portion, never
  extrapolate.

Evidence labels fixed now: exact finite-field arithmetic artifacts are
MACHINE-VERIFIED; wall/CPU times are COMPUTATIONAL-EVIDENCE only; any
derived alpha (only if the fork fires) is CITED-DEPENDENCY. No
universal over m or t is claimed; "m=12 complete" is NOT asserted
(one row per t, t ∈ {64} previously, {96} here); the Apon section-3.6
all-Δ-zero hole remains BOUNDED by the Gate-B rates, NOT closed, and
is not touched by this campaign; no NIST attack-cost figure is
verified; t=48 at m=12 remains UNREACHED.

## 7. Scope sentence (rule 7, exact)

Swept in this campaign: the single instance (12, 3488, 96) at seed
16384 — one seeded irreducible G (scalar-Rabin-confirmed), one seeded
support ordering; alpha over the FULL support grid (3488 points,
inherited identical instrument) + the 5-part exact delta chain +
build guards + 3 seeded F₂-linearity selectors + λ/gcd controls; the
3 registered anchor rows (11,2048,48,6211), (10,1024,40,5113),
(6,64,3,1387) in byte-equality mode vs the frozen 14-09Z records; the
3 registered counterfactual plants at the (11,2048,48,6211) base plus
the bounded 128-point t=96 probe-grid direction check (P-4). NOT
swept: t = 48 or any other t ≠ 96 at m = 12 (t = 64 was the PRIOR
campaign's single cell; it is not re-run here), any n ≠ 3488, any
support ordering other than the single seed-16384 shuffle, any m > 12,
any census row (N_fam / per-point rank / held-count ladder) at any
cell, any sub-support construct, any non-irreducible G, no G
population exhausted, any reuse/extension of the t=64 cell (its
records stay frozen in their own campaign dir). This campaign's
positive result is a BOUNDED FINITE statement about ONE new instance
at (m,n,t) = (12,3488,96): it is NOT "m=12 complete", NOT a
verification of any NIST-cell attack-cost claim, NOT a section-3.6
closure, and adds no universal quantifier over any family.

## 8. Named next action (pre-declared)

If this t=96 verdict lands PASS: the named next action is the
remaining m=12 row (t=48), or — at Main's discretion — the m=12
waterfall census (N_fam / held-count ladder), each under its own fresh
pre-statement budgeted from the now-two measured m=12 rows. If FAIL:
escalate to Main before any external write.

— McelieceM12T96, 2026-09-01, committed BEFORE any t=96 compute.
