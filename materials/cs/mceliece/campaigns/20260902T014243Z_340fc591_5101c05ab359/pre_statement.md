# Pre-statement — third m=12 row at t=48: alpha-first guard campaign (bounded cell 12/3488/48), mceliece target

**Committed 2026-09-02, BEFORE any t=48 computation.** Fixes the exact
`(m,n,t,seed)` cell, the instrument (INHERITED UNCHANGED from the
certified `M12ALPHA_DIVONLY` campaign — the identical code path that
produced the m≤11 ladder and the t=64/t=96 rows; no second convention),
the alpha mode and its MEASUREMENT-GATE fork budgeted from the two
MEASURED m=12 rates, the anchors, the plants (including a t=48-native
probe), the budget, and the verdict logic. Implements the named next
action after the certified t=96 campaign
(`campaigns/20260901T100025Z_bd55fdfd_8fbba9be5ab3/`, FROZEN-CERTIFIED):
the last unreached m=12 ladder cell (12, 3488, 48 | k = 3488 − 576 =
2912, D = 3488 − 97 = 3391, 2t+3 = 99).

## 0. Registered instance and population (fixed before compute — rule 16)

Exactly ONE m=12, t=48 instance is registered:

    (m, n, t, seed) = (12, 3488, 48, 16384)

- k = n − mt = 3488 − 576 = 2912; D = n − 2t − 1 = 3488 − 97 = 3391;
  2t+3 = 99; |E| = 4096; support n = 3488 (as in every certified row).
- seed 16384 is fixed BY THIS SENTENCE, before compute — the same seed
  as the t=64 and t=96 registered rows (the seed stream draws depend on
  t through G's degree, so the instances are distinct objects; the
  rule is uniform). Build via the UNCHANGED `src/instance.py` seed
  stream. If any build guard fails, regenerate with seed 16385,
  16386, …; every substitution is logged in an explicit dated
  amendment BEFORE the substitute verdict is computed. Guards are
  never adjusted.
- NOT registered, hence never computed in this campaign: any t ≠ 48 at
  m = 12 (t=64 and t=96 records stay frozen in their own campaign
  dirs — no re-run, no reuse), any n ≠ 3488, any other support
  ordering, any m > 12, any census row (N_fam / held-count ladder) at
  any cell, any sub-support construct, any non-irreducible G.

## 1. Instrument: INHERITED, byte-identical; no second convention

`instance.py`, `fastfield.py`, `gfield.py`, `census.py` copied
byte-for-byte from current `src/` into `code/` at init (before any
compute) and verified by `cmp -s` against BOTH frozen references used
by the two prior m=12 campaigns: the `M12ALPHA_DIVONLY` snapshot
(byte-expected identical) and the `2026-08-30T14-09Z_57200ADD` ladder
snapshot (expected identical except instance.py's inert `forced_G`
ctor parameter, lines 84–91 — the artifact disclosed by both prior
campaigns). Any other diff: STOP, escalate to Main before compute.
Guard-chain semantics and the verdict stage are mechanical
re-instantiations at the new cell, exactly as in the t=96 campaign.

## 2. Alpha mode: MEASURED FIRST, fork only on a FAILED measurement gate

Identical structure to the t=96 pre-statement (section 2 there), with
the abort threshold re-budgeted from the two MEASURED m=12 rates:

1. ALPHA IS MEASURED FIRST: the identical m≤11 instrument's
   per-support-point alpha check runs over the FULL support grid
   (n = 3488) on the registered t=48 instance.
2. Budget basis (MEASURED, in-process, both prior rows): t=64 cell —
   build 2175.43, alpha grid 700.53, total 2893.37 CPU-s;
   t=96 cell — build 3351.79, alpha grid 648.43, total 4017.55 CPU-s.
   The alpha per-point cost scales with the D-truncated vector width,
   so t=48 (D = 3391) sits BETWEEN the two: D(t=48) = 3391 vs
   3359 (t=64) vs 3295 (t=96). The t=48 build cost is expected between
   the two measured builds (G degree 48 is smaller; the GF(2)-linalg
   part of the build is unchanged; the interpolation XOR work is
   k×(D+1)-scaled, k = 2912 between 2720 and 2336). Projection (NOT a
   decision input): total ≤ ~4200 CPU-s.
3. Registered HARD budget: **9 CPU-hours = 32,400 CPU s
   `time.process_time()`, total for the verdict block.** Provision:
   a full 3x multiple of the expected cost, identical to the t=96
   campaign's cap, so the two cells stay comparable.
4. DERIVE-branch entry condition (registered NOW, judged only on
   measured cumulative CPU inside this process): at the alpha-grid
   checkpoint, if cumulative process_time() has already exceeded
   75% of the budget (0.75 × 32400 = 24,300 CPU s), abort the grid and
   execute the DERIVE branch: alpha labelled CITED-DEPENDENCY
   (inferred from the exact 5-part delta chain + Apon Lemma 3, κ = 1),
   with the mandatory separate machine checks (a) gcd(G,G′) constant,
   (b) deg Π = n with Π′ nonzero, (c) max_j deg f_j ≤ D ABORT gate,
   (d) λ recomputed byte-equal; the resolution recorded in a dated
   checksummed AMENDMENT file BEFORE the delta chain compute. Never
   resolved after seeing alpha verdict numbers beyond the aborting
   checkpoint. Expected resolution (projection only): MEASURED, as in
   both prior rows, with ≥ 7x headroom.

One-sentence kind labeling, verbatim, either way: "at m ≤ 11 the
mechanism is MEASURED; at m = 12 [this row] it is MEASURED / it is
INFERRED from a cited lemma given a measured delta."

## 3. Anchors (HARD GATE, before any t=48 number is reported)

The SAME three anchors as both prior m=12 campaigns, re-run through
THIS campaign's `code/` snapshot and compared field-by-field against
the frozen `2026-08-30T14-09Z_57200ADD` records:

    (11, 2048, 48, 6211)   — largest frozen row
    (10, 1024, 40, 5113)   — large frozen row
    (6,   64,  3, 1387)    — smallest frozen row

byte-equal for every field EXCEPT `elapsed_s` (timing exclusion) and
disclosed additive fields (`guards_build` plus max_deg_f/deg_ok,
cross-checked for internal consistency). All three must pass before
any t=48 verdict; any mismatch: escalate to Main, no verdict.
Anchor rows keep their m≤11 MEASURED alpha labels.

## 4. Counterfactual plants (HARD GATE, rule 14)

The t=96 suite re-instantiated verbatim (rounded to P-1..P-4 names),
all rejected by the same inherited verifier, at the m=11 anchor base
(11,2048,48,6211) except P-4:

  - P-1 duplicate support (a_12 := a_11): expect Π′(a)=0 at the
    duplicated positions; build distinct-support guard fires.
  - P-2 f_0 → f_0 + Π: degree gate ABORT (deg = D+1 > D), NO delta
    verdict; value-level invisibility prediction confirmed (0 failures
    on probed points).
  - P-3 row 0 → 3·f_0: alpha value check rejects (c = 3, c² = 5 ≠ 3,
    verified in-run); delta-binary side-trip recorded if it fires.
  - P-4 t=48-native bounded probe: the REGISTERED (12,3488,48,16384)
    instance must be alpha-clean on the first 128 support points and
    its row-0 ×3 corruption must fail ≥ 1 of them (same verifier ops,
    same corruption mechanism). Instrument-direction check at t=48
    scale; not a second convention.

Any expected rejection failing to materialize: FINDING, escalate to
Main before any verdict is written.

## 5. Budget and abort gates

- Instrumentation: `time.process_time()` / `time.perf_counter()` in
  process; `nice -n 10`; BLAS/OMP threads = 1; no `ps`-based claims
  (inherited Main steering).
- 9 CPU-h hard cap (section 2.3); on breach: freeze the completed
  sub-blocks, checksum, close FROZEN-INCONCLUSIVE over the bounded
  completed portion — never extrapolate.
- Wall budget: 24 h (the t=96 wall was ~9.5 h dominated by the build;
  the t=48 build is projected shorter, but the cap stays generous).

## 6. Verdict shape, pass criterion, evidence labels (fixed now)

- PASS (FROZEN-CERTIFIED, bounded to the cell): build guards pass;
  delta chain 5/5 exact (division 3488/3488 zero-remainder, unit
  checks 3488/3488, 0/200 offdiag, 3/3 assembly, degree gate
  max deg f ≤ D — expected max deg f = D = 3391 exactly as on both
  prior rows); alpha MEASURED full grid 0 failures (or, only if the
  section-2.4 gate fires, CITED-DEPENDENCY via the declared derive
  branch); controls pass (gcd(G,G′)=0, λ byte-equal, F2-linearity on
  4 seeded selectors, β pair exhibited); 3/3 byte-equal anchors; 4/4
  plants rejected; all within budget.
- FAIL: any registered check fails as expected-runs-false: FINDING,
  escalate to Main via hub BEFORE external writing; close
  FROZEN-NEGATIVE (registered-verdict failure) or
  FROZEN-INCONCLUSIVE (instrument/integrity failure only).
- ABORT: degree-gate failure = NO delta verdict (recorded as such);
  budget breach = bounded freeze per section 5.
- Evidence labels: exact arithmetic artifacts MACHINE-VERIFIED;
  wall/CPU times COMPUTATIONAL-EVIDENCE; a derived alpha (only if the
  fork fires) CITED-DEPENDENCY. No universal over m or t; with this
  row the m=12 ladder is THREE bounded instances, t ∈ {48, 64, 96} —
  still NOT "m=12 complete", still no closure of the Apon section-3.6
  all-Δ-zero hole (bounded by Gate-B rates, not closed), no NIST
  attack-cost verification.

## 7. Scope sentence (rule 7, exact)

Swept: the single instance (12, 3488, 48) at seed 16384 — one seeded
irreducible G (scalar-Rabin-confirmed), one seeded support ordering;
alpha over the FULL support grid (3488 points, inherited identical
instrument) + the 5-part exact delta chain + build guards + 4 seeded
F2-linearity selectors + λ/gcd controls; the 3 registered anchor rows
in byte-equality mode vs the frozen 14-09Z records; the 4 registered
plants (P-1/P-2/P-3 at the (11,2048,48,6211) base, P-4 bounded
128-point probe at the t=48 cell itself). NOT swept: t = 64 or t = 96
at m = 12 (frozen records in their own campaign dirs — this campaign
does not re-run or extend them), any other t, any n ≠ 3488, any other
support ordering, any m > 12, any census row, any sub-support
construct, any non-irreducible G, no G population exhausted. A
positive result is a BOUNDED FINITE statement about ONE new instance:
it completes the pre-registered m=12 ladder as three bounded
instances (t ∈ {48, 64, 96}, one seed each) and is NOT "m=12
complete", NOT a section-3.6 closure, NOT a NIST-cell verification.

## 8. Named next action (pre-declared)

On PASS: the m=12 ladder (t ∈ {48, 64, 96}) is fully measured under
this code path; the named next action is Main's choice between (a) the
m=12 waterfall census (N_fam / held-count ladder) at one or more of
the three cells, or (b) standing down with the three bounded
instances as the delivered record. On FAIL: escalate to Main before
any external write.

— McelieceM12T96, 2026-09-02, committed BEFORE any t=48 compute.
