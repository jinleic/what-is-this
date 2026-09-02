# PRE-STATEMENT — OmegaPcVersusPen: certified subdomain pinning of the
# p_comp-vs-penalty ORDER inside the dist[0] 1e-7 box (exact-kernel slice).
# Agent OmegaPcPen (worker), 2026-09-01. Committed to git BEFORE any
# parameter-dependent compute (Main strict order). No campaign init before
# this file only reads of frozen artifacts, tool inspection, and the
# source-wiring structural read happened.

## 0. Route selection, from the frozen ledger only

- Route "uniform branch-and-bound over D" is FROZEN INFEASIBLE-AT-COST
  (campaign `2026-09-01T03:15:00Z_OmegaGateCBnBGate_9d2e4a7c`: k* = 36 of
  45 coordinates at 90% of L1 width; ~2^223 subboxes; stage 2 never
  opened). NOT rerun here.
- Route "analytic branch-pinning by freezing eps" is FROZEN FALSIFIED
  (campaign `20260901T102048Z_b555441e_d3303018f163`, verdict
  FROZEN-NEGATIVE): both registered order inequalities
    G1 = cand_1.lo - cand_0.hi and G2 = cand_2.lo - cand_0.hi
  take outward intervals [-2.9602186479446094e-05, +2.6689872684820378e-05]
  resp. [-2.9602189567791676e-05, +2.668987624020422e-05] over the FULL
  un-subdivided D+ box — a width of 5.629e-05 = 355.9x the registered
  first-order signal 1.5816497000997742e-07 (187.2x in the negative lobe).
  eps_0 is certified BIT-CONSTANT over that box (2.1502086942121845e-06 at
  r=0 and r=1e-7) and fully dual-coupled (C2 boundary moves by exactly
  4*drift, coupling gap 4.41e-16; R_glob[0] responds -2*prop*4*drift, gap
  2.9e-16). eps_0 is therefore OFF the movers list, and the straddle is
  carried by the p_comp functional pc_0 (box motion ~1.538e-5) racing the
  penalty pen_0 (~1.142e-5). NOT rerun here.
- Registered successor (BnBGate manifest + RESULTS.md 2026-09-01): "a
  successor must pin a p_comp-vs-penalty ordering, or pin a different
  branch." THIS campaign pins the p_comp-vs-penalty ordering — the
  pc_0(R) >= pen_0(R) family — but on the EXACT-KERNEL SLICE of the box
  (the 21-dimensional feasible set itself), not on the full axis-aligned
  superset box where the falsified campaign measured. This is new
  technology in three concrete ways (no rerun of either dead route):
  (i)   the swept set is the true feasible subdomain — the frozen exact
        Rational kernel basis V (from CandidateWitness, sha-pinned) with
        certified feasible-t extreme points — not the hulled box;
  (ii)  the functional evaluated is the IRREDUCIBLE core of the ordering
        — a one-dimensional certified family mu(rho) := pc_0 re-term(rho)
        minus pen streak(rho) term(rho) evaluated on a certified box
        around the frozen ray delta* = r0·s (the LP midpoint candidate
        direction lifted to the kernel) with rho in [0, 1];
  (iii) the certified object is a pc>=pen ORDER on subdomains, which is
        the step the falsified campaign explicitly could not take.

## 1. Exact domain, subdomains, and the swept set (rule 16)

- Base point p*: the VXXZ24 released vector K100_2.37155181.mat (OSF
  7wgh2; params sha256 df75ae3acaa5b1388bd2f17cce266aec169d874c9fea1eb0
  722c3e25871da93d — first 64 hex prefix asserted IN-RUN before any use).
- Transcribed program: 3-region single-p_comp (vxxz24_float 31657c92…),
  max_level 3, q 5, K 1, IC v11 MID 300 OUTB 64. Every parameter block
  except the 45 dist[0] coordinates is held at p* (POINT intervals,
  width-0 asserted IN-RUN per leaf — the falsified campaign's P1
  invariance check is RE-USED VERBATIM as control C1a here).
- dist[0] box: the same certified box as the parent campaigns — 45
  exact-Arb dyadic endpoints IC.const(center_i) ± IC.const(radius),
  radius exactly binary64 1e-7 = 0x1.ad7f29abcaf48p-24. Call this box B.
- L_SUBDOMAIN (the certified subdomain — these are the precise registered
  subdomains of the pin): inside B, the exact affine kernel slice
    D = { delta in R^45 : A·delta = 0 exactly, |delta_i|_inf <= 1e-7 },
  with A the frozen 27x45 0/1 margin matrix (rank 24 over Q — re-proved
  in-run by exact Fraction elimination on the sha-pinned npy) and the
  frozen exact-Rational kernel basis V (21 cols, each an exact
  Fraction vector; sha256 7d69c1dd4c9355dfe6d6b3a762658ba10d93e4bcaea3b
  240ebbbd811de3eba60 over the JSON file asserted in-run).
- THE MAIN SUBDOMAIN: the certified feasible-t single-coordinate family
    D_mu := { p* + rho·s : rho in [0, 1] },
  where s := theta·v (theta = the 21 binary coefficients
  [1,-1,0,-1,1,0,1,0,0,1,-1,-1,-1,1,0,-1,0,0,-1,0,1], EXACTLY the
  CandidateWitnessray u) and
    r0 := 1.0607145909782878e-06
  is the certified feasible-t ray radius registered in §2. D_mu is
  D-reachable, 1-dimensional, and the certified mu-boxes in §2 are built
  such that for every real rho in [0,1] the point lies in D exactly.

## 2. The certified core bounds (NEW instrument, no float dependency)

Certified feasible-t ray radius: for the 21-coefficient binary ray
s = theta·v, the positivity wall is
  R := min( 1e-7 / max_i |s_i| , min over s_i < 0 of c_i / (-s_i) )
with c_i the exact Fraction values |p*_i| in dyadic form (all dist[0]
coordinates in the released vector are positive). Evaluate BOTH
branches in exact Fractions; take the min; ROUND DOWN to a 31-bit
floating dyadic, recording the exact rational bound R_exact and the
dyadic r0 <= R_exact with a certified gap. The float r0 1.0607145909e-06
is the EXPECTED value (from frozen protocol_static_corrected
positivity_ceiling arithmetic) but is RECOMPUTED here in-run by the
exact-Fraction instrument; the exact value is authoritative, not the
float transcription. If recomputation differs from the expectation by
more than 1e-12 relative, the campaign STOPS at E-FIRST (registered
below) — expectation fails, machinery records and exits.

The certified mu-box (e0, e1): with r0 fixed by the exact instrument,
    e0 := [IC.const(0), IC.const(0)]  (degenerate point at rho=0)
    e1 := IC.const(r0) ± 2^-100    (center r0, certified radius 2^-100,
        outward buffer 1e-30 enforced by the v11 buffered .low()/.up())
THE RAY BOX IS THE 1-DIMENSIONAL FAMILY WITH ENCLOSED ENDPOINT r0.
The certified intersection with the kernel constraint: A·(rho·s) = rho·A·s
= 0 for every real rho because A·v = 0 col-by-col over Z (asserted
in-run by exact integer multiply on the sha-pinned npy and the exact
basis, and because s is an integer combination of kernel columns,
A·s = 0 EXACTLY in Z^27).

## 3. The functional mu (pre-declared build; no existing instrument)

The certified generator set U(s, rho): read vxxz24_float source wiring
(no randomness): pen_0(R) = prop_0*(H(dm0) - H(d0(R)) + 2*eps_0) with
eps_0 POINT-constant per the frozen C1a-invariance fact (its dm0 and lam
argument blocks are p*-frozen; the in-run leaf-width-zero assertion IS
the certificate — C3 re-verifies live with the same threshold the frozen
campaign used: coupling gap <= 1e-12).

pc_0(R) is the p_comp[0] of GlobalStage — a functional of the moving
dist[0] box through the parts' complete_split chains and region
marginals: num = sum_i dist0[i] * (entropy of part cs[dimW]) for
zero-coord shapes, plus the weighted[k] normalized-entropy terms;
den = entropy(complete_split[0][0]) - entropy(J2M-marginal-0 of dist0).
Both generator families are read from the v11 dependency-chain core
(hash 31657cec92… asserted in-run) but stitched into a NEW certified
composite: the v11 core computes pc_0 only through the full workspace
pass; this campaign's mu-generator binds the workspace pass to a
one-parameter family and produces the certified DIFFERENCE
    mu(rho) := pc_0(p* + rho·s) - pen_0(p* + rho·s)
as a single outward interval function of rho.

Registered generator semantics (identical contract to the falsified
campaign's, meaning FULL compatibility with the frozen records): the
underlying aggregation identity is REUSED VERBATIM from the frozen
branchpin.py (sha d7aca589c5dafbf7942f2b4009849edc311071e0262606695c474
618f7ba3c56): R_glob[0] candidate structure cand_tt = nb_tt - (pc_0 or
pen_0), and the registered G1/G2-style ORDER QUESTION
    G1(rho) := cand_1.lo(rho) - cand_0.hi(rho)
    G2(rho) := cand_2.lo(rho) - cand_0.hi(rho)
are NOW computed as explicit functions of rho on the 1-DIMENSIONAL
kernel-reachable family D_mu above — NOT over the full D+ box. The
candidates are the R_glob[0] candidates from the frozen aggregation —
identical code path (stage_b_rung2 / branchpin.py aggregate()).

DELIVERABLE^: the certified statement
  "the pc-vs-pen ordering G1(R) > 0 AND G2(R) > 0 holds for ALL rho in
  [0, 1]" (subdomain pin PASS on D_mu), OR the exact failing
  sub-domain: the measured set of rho where the interval straddles or
  fails, with the factor versus the signal
  1.5816497000997742e-07 (numeric low and width both recorded).

Infeasibility-cohort guard (mandatory fairness to the frozen negative):
because mu's pc_0 term moves through part_prob_sum chains, a radius-0
point pass (rho=0) MUST REPRODUCE the frozen GATE-Z records and the
frozen candidates (frozen cand0 = [0.4968507617380737, 0.4968507617380738]
etc., within 1e-9 — the falsified campaign's frozen tolerance) BEFORE
any rho > 0 pass. Any drift kills the campaign at that control.

## 4. Controls (registered obstacle-rejection gates)

- GATE-Z (the GATE-Z equivalence reproduction of both frozen records,
  stage_b_rung2 with/without Lemma-1, R-branch and Om_raw tolerances
  1e-12 / 1e-6 as frozen): PASS gate, in-run, radius-0 point pass.
- C1 (known-true ACCEPT — local identity): the radius-0 candidate-gap
  identity cand1.lo - cand0.lo = -2*prop_0*eps_0 must reproduce to
  <= 1e-12 gap (frozen value -1.4334718842490268e-06 vs registered
  -1.4334724628081230e-06, gap 5.79e-13 — the falsified campaign's
  frozen C1, re-verified in-run as the accept control). PASS REQUIRED.
- C1a (point-invariance of the eps leaf): all non-dist0 leaves must
  carry width-0 point intervals in the rho>0 passes. Asserted by exact
  leaf introspection before P2. FAIL = STOP (wiring defect).
- C2 (structural REJECT — wrong-index dual lookup): region-1 dual lookup
  must return eps_1 = 1.143325024764798e-06 (frozen with-lemma1 record)
  and NOT eps_0. Both relations asserted. PASS REQUIRED.
- C3 (sign/branch REJECT — dual drift, +1e-3 broadcast on every glob-0
  lam entry, drift region 0): the eps_0 boundary must move by EXACTLY
  4*drift (coupling gap <= 1e-12; frozen registered semantics — the
  eps residual becomes strictly positive under drift) AND R_glob[0]
  must move by -2*prop*4*drift within gap <= 1e-12. This is a
  REJECT-side plant: any positive G1/G2 under this perturbed dual must
  be distinct from the base verdict; equal-to-base G1/G2 under the
  plant would mean the machinery is insensitive to the Lemma-1 coupling
  and the campaign reports REJECT (instrument rejected, no PASS
  possible).
- P-quantifier control (NEW, the know-true object for the order): at
  rho = r0 exactly (the certified endpoint), the certified G1/G2
  OUTWARD intervals are computed TWICE through two independent code
  paths: (i) the direct rho-box pass (§3), (ii) the subdomain pair-sum
  identity with the mu-generator difference (§3's stitcher). The two
  outward intervals must MATCH to a registered relative tolerance
  (min width/width < 1e-9, both must be nonempty intervals over the
  agreed quantifier domain) before the deliverable is recorded —
  a subdomain PASS that fails to reproduce through the stitcher is
  FALSIFIED (the order claim is not robust).
- Nothing else: no slope/gradient pass, no LP solver, no Monte Carlo,
  no subdivision of B, no re-run of B&B (the one-dimensional family is
  NOT a subdivision of the box — only two parent-sized passes exist in
  the main protocol (P1 at e0, P2 at e1) and each is a single interval
  pass with 45 moving leaves, exactly as heavy as the frozen
  falsification's single pass, measured ~3 s by its checkpoint).

## 5. Backend, precision, rounding, budget, kill logic

- Backend: python-flint (cs/.venv/bin/python), arb endpoint-pair Ival,
  MID 300, OUTB 64, EXACT-Arb dyadic box endpoints via IC.const.
  All shears of the mu-generator computed with buffered .low()/.up()
  (v11 semantics). sys.dont_write_bytecode = True BEFORE importing
  anything that lives inside any campaign/frozen directory; PYTHONDONT
  WRITEBYTECODE=1 in the environment; __debug__ asserted (the v11 core
  refuses -O).
- Every accepted claim uses outward-rounded intervals only; no float
  endpoint participates in any decision inequality. All averaging into
  float fields is display-only and separate from decisions.
- Budget: hard wall 30 minutes CPU under `nice -n 15`
  (OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1; several
  campaigns share the machine), single process. Checkpoint (write-first
  JSON under a DISTINCT per-attempt filename, freeze-hygiene rule)
  prior to each control after GATE-Z. The full P-phase costs ~2
  parent-sized passes (~6-10 s measured parent cost); the budget is
  30 min = 180x headroom.
- KILL logic (kill = verdict OPEN with the failing step named; the
  checkpoint and per-attempt stderr are preserved):
  - any of GATE-Z, C1, C1a, C2, C3 fails or the C1a leaf introspection
    detects nonzero-width non-dist0 leaves -> OPEN at that step
    (verdict OPEN, not FALSIFIED: the machinery is defective, the
    problem is not);
  - exact r0 recomputation disagrees with the expected float by >1e-12
    relative -> OPEN at E-FIRST;
  - the P-stage pair-path mismatch (P-quantifier control) fails ->
    FALSIFIED (the subdomain order claim is falsified for THIS
    instrument; recorded as a quantified negative with the factor
    versus signal, and eps_0 remains exonerated per the frozen record).
- TERMINAL verdicts (exactly one at close):
  - FROZEN-CERTIFIED (the registered warm outcome): all five controls
    pass AND both G1/G2 outward intervals are STRICTLY positive
    across the WHOLE certified rho family (both e0 and e1 passes; any
    overlap in the one-dimensional family is RULED OUT by the
    certification that mu(rho) is verified to be increasing, below):
    deliverable = "on the DIST-REACHABLE one-dimensional subdomain
    D_mu of the 21-dim slice, pc_0(p*+rho s) >= pen_0(p*+rho s)
    strictly, for all rho in [0, r0]; hence the p_comp candidate is
    the unique possible R_glob[0] minimizer on that subdomain" with
    the numeric margin recorded; this is certified LOCAL progress on
    the ordering question; the 21-dim kernel question and the global
    21-dimensional kernel statement remain UNCHANGED AND OPEN
    (D_mu is 1-dimensional; nothing here decides the remaining
    20-dim space); the 12.808x guardrail restated (below).
  - FROZEN-NEGATIVE: any of the six registered inequalities fails —
    the deliverable is the exact failing inequality, its outward
    interval, the factor versus the registered signal
    1.5816497000997742e-07 (numeric low/width quotient, recorded),
    and the named mover carrying the failure (the pc-vs-pen family
    likely carries it since eps is exonerated). The mover is
    REMOVED from the movers list as a first-class frozen negative
    (the route "pin the ordering on the kernel ray family at
    1.06e-6 radius" is measured-false at the recorded factor).
  - OPEN: budget/kill (named step).

## 6. Rule-7 sentence (scope of every statement here)

The swept set of THIS campaign is exactly the 1-dimensional polytope
segment D_mu = { p* + rho·s : rho in [0, r0] } with s the integer
combination of the frozen exact rational kernel basis (21 coefficients
in {-1,0,1}, A·s = 0 over Z) and r0 the certified feasible ray radius
(exact-Fraction instrument, target 1.0607145909782878e-06, recomputed
in-run), together with the radius-0 control point p* itself — all
inside the dist[0]-block 1e-7 box D+ around the region-0 glob dist of
the VXXZ24 K100_2.37155181 released vector, under the transcribed
3-region single-p_comp program at max_level 3, q = 5, with the 45-leaf
dist[0] box moved exactly as the parent campaigns moved it and every
other parameter block held at p* point values (width-0 asserted),
interval core v11 (ce70f959…), vxxz24_float (31657c92…). NOT swept: the
rest of the 21-dimensional kernel slice (the other basis directions,
their cones, any 2+-dimensional subdomain), any other radius, any other
region or parameter block, any other branch of R beyond R_glob[0]'s
three candidates, any subdivision, any gradient machinery, the
published 2.37155181 and the record 2.371177 (never searched;
parameters unpublished). Even total success caps at the 12.808x
guardrail: the certified pc-vs-pen ordering margin lives at the ~1e-5
scale while the gap between the certified rung-2 endpoint
2.3715538358544617 and the published 2.37155181 is 2.0258350805768544e-6
= 12.808x the signal 1.5816497000997742e-7 — so NO outcome of this
campaign produces an exponent improvement, a record claim, or any
statement about the VXXZ24 paper's validity; gate C (local optimality
over D) stays PARTIAL regardless of verdict.

## 7. Explicit not-swept list (mandatory)

- The 20-dimensional remainder of the kernel slice (any direction other
  than the single registered s), the opposite ray rho < 0, any rho
  outside [0, r0], any other subdomain of B outside D.
- Regions 1/2 glob blocks, all part-level blocks, all lam blocks beyond
  what C2/C3 read (they are query-only plants), the omega x
  single_mat_size Schoenage 2-coordinate class, region_prop.
- No B&B re-run, no LP solver run, no slope pass, no subdivision, no
  Monte Carlo, no float-precision shortcuts (exact-Arb endpoints only).
- The published 2.37155181 and the record 2.371177: untouched; nothing
  here bears on the validity of the VXXZ24 paper.

## 8. Adjudication, verbatim (fixed before compute)

- PASS requires ALL of: GATE-Z pass, C1 pass, C1a pass, C2 pass, C3
  pass, E-FIRST consistent (exact r0 within 1e-12 rel of expectation),
  P-quantifier control pass (dual-path agreement), and the SUBDOMAIN
  ORDER both strictly positive outward at BOTH endpoints e0, e1, with
  e0.e1 cross-verification. Everything else in the P-phase is FROZEN-
  NEGATIVE with the failing inequality, its interval, and the
  signal-factor recorded. OPEN only on budget/kill.
- The deliverable margin numbers: for the PASS branch, the numeric
  minimum of the four outward G-interval lows (G1@e0, G1@e1, G2@e0,
  G2@e1) and its factor 1/(signal) quoted as a dual:signal relation
  (display only). For the FAIL branch, the failing inequality's
  low/width and its width/signal factor.

## 9. Mandatory record-binding quotes (frozen ledger; preserved verbatim)

- BnBGate: "k* at >= 90% of the L1 width is 36 of 45 coordinates …
  ~2^223 subboxes to force width/signal < 1 … INFEASIBLE AT THIS COST —
  frozen as a first-class quantified negative".
- BranchPin manifest: "both order inequalities G1/G2 straddle zero over
  the full 1e-7 dist[0] box, interval [-2.9602e-05, +2.6690e-05] —
  355.9x the 1.5816e-07 signal in width"; "the straddle is carried by
  the pc_0 (p_comp) hull moving ~1.54e-5 against the pen hull moving
  ~1.14e-5 over the box"; "eps_0 is CERTIFIED CONSTANT over D
  (bit-identical to the frozen point value 2.1502086942121845e-06) and
  is therefore EXONERATED as the box-straddle driver"; "eps_0 is
  RADIUS-INVARIANT over D (bit-exact) — any future box analysis can
  treat the glob Lemma-1 residual as a frozen constant, removing it
  from the movers list; the movers that remain are pc_0, nb_* and
  (d_box) entropy terms".
- Gate C standing frame: certified box enclosure
  [2.3715431556360005, 2.3715624359361276] (rigorous, 122x looser than
  signal); midpoint LP candidate -1.577332901516852e-7 (COMPUTATIONAL-
  EVIDENCE ONLY, not rigorous); ANY-y enclosure width 1.168e-5
  straddles zero (MACHINE-VERIFIED); 12.808x guardrail; the
  21-dimensional kernel question UNCHANGED and OPEN; gate C PARTIAL.

## 10. Commit-bound preregistration self-checks (registered NOW)

The following data are pinned BY THIS FILE at commit time; the runner
re-derives each in-run and compares:
(a) sha256(interval_core.py) = ce70f95922f3d5098e892e1b7047570ebb526782
    6b51f2776c91b56d9934adb5;
(b) sha256(vxxz24_float.py) = 31657c92a1f841805e7355d5b5390f76e0d7ca2d0
    4b1a3f2ce9fc204d9d5d890;
(c) params mat sha256 prefix df75ae3a (full: df75ae3acaa5b1388bd2f17cce
    266aec169d874c9fea1eb0722c3e25871da93d);
(d) kernel basis JSON sha256 = 7d69c1dd4c9355dfe6d6b3a762658ba10d93e4bc
    aea3b240ebbbd811de3eba60;
(e) frozen FROZEN_R_DETAIL / FROZEN_RAW / FROZEN_RAW_OFF / FROZEN_C /
    FROZEN_EPS0 / FROZEN_EPS1 exactly as registered in the falsified
    campaign's branchpin.py (sha d7aca589…); reproduced in-run to 1e-12
    (branches) / 1e-6 (Om);
(f) expected candidate triple at the base point (CandW):
    cand0 [0.4968507617380737, 0.4968507617380738],
    cand1 [0.49684932826618944, 0.4968507617386526],
    cand2 [0.4968493282665137, 0.4968507617389769].
