# PRE-STATEMENT — gate C stage 4 gate: width-decomposition feasibility
# gate for branch-and-bound (agent OmegaPenaltyProbe, 2026-09-01T03:15Z)

Written and committed BEFORE any compute of this campaign. Parent: frozen
`campaigns/2026-08-31T02:40:00Z_OmegaGateCFinish_v1_4c8a2f1b/` (V1 closed
under the v11 core; LP/ANY-y FAILURE TO CERTIFY; width/signal
73.85977508485684). Steering: Main 2026-09-01 (feasibility-gate-first
order, thresholds fixed in advance, quantified obstruction if the gate
fails).

## 0. Name and scope

Campaign dir: `campaigns/2026-09-01T03:15:00Z_OmegaGateCBnBGate_9d2e4a7c`
— stamped with the ACTUAL UTC start (per the provenance rule adopted in
the parent). One diagnostic measurement of an EXISTING enclosure + a
pre-registered two-branch verdict. No new mathematics: every number
consumed is re-derived from the already MACHINE-VERIFIED slope bundle of
the parent run (`run_v1_v11_checkpoint.json`, sha
a3c8a08f78a46093f4d3c1aaffc6addc4eeabc6da9e5b6edf38f83e85cbc83fa) plus
one targeted recomputation of per-term components (below), labelled
MACHINE-VERIFIED where recomputed in exact interval arithmetic under the
same v11 core, COMPUTATIONAL-EVIDENCE for the derived subbox-count
extrapolations (they are arithmetic on measured quantities, not new
enclosures).

## 1. STAGE 1 — width decomposition (the one pass; cheap)

Fix the 45 coordinates K = {0..44} of the region-0 glob dist[0] block
(parent radius 1e-7, frozen centre, exact-Arb endpoints).

1.1 Per-coordinate width table: w_i = hi_i - lo_i from the frozen
`grad_lo`/`grad_hi` arrays; sort descending; cumulative-share curve;
report the smallest k with cumulative share >= 90%, and the shares at
k = 2, 4, 6, 10, 23.

1.2 Per-term / per-block attribution (targeted recomputation, one point
pass + one box pass reused from the parent state where possible): for
the four aggregation families — level-3 num_block per (r, tt), level-3
penalties (hash_penalty) and p_comp, level-2 nb2, glob num_block /
penalty / p_comp / mat_size — compute each term's contribution to the
L1 width sum 113.66562530936552 by re-running the same asrun aggregate
extraction used in the parent (byte-reuse of aggregate()) and recording
each block's grad-interval width L1 norm BEFORE the min_slope hulling,
plus the hull-added width per R branch (the "possible_minimizers"
hull cost) — the hull-overhead term is reported separately because it
is the part subdivision actually reduces.

1.3 Branch-straddle component: from the parent box branch data
(`slope_box.R_branch_hulls`), report per branch (l, r), lvl2, (G, r):
candidate_intervals, whether the candidates overlap (straddle of the
MINIMIZER IDENTITY), and the interval span of each candidate. Tie this
to the Lemma-1 straddle evidence frozen in GateCDiag/CandidateWitness
(R interval span ~4.2506006968e-6 across the seven components; Lemma-1
residual intervals straddle zero at the base point).

## 2. STAGE 1 VERDICT — thresholds fixed now (Main's rule, made exact)

- CONCENTRATED branch: some coordinate set of size k* <= 6 carries
  >= 90% of the L1 width sum (k* read off 1.1). Then proceed to STAGE 2
  with subdivision restricted to EXACTLY those k* coordinates; depth
  cap 8 (=> <= 2^8 = 256 subboxes per staged sweep, cost calibrated
  first); split rule: bisect the widest binding coordinate of the
  subbox at its midpoint (exact dyadic midpoint of the exact-Arb
  endpoints); a subbox is CLOSED when its certified enclosure width and
  its pinned branch structure produce a width/signal < 1 statement or
  the branch identity is decided.
- SPREAD branch (the obstruction branch): NO set of <= 6 coordinates
  reaches 90%. Then STOP — no branch-and-bound run. Freeze the
  QUANTIFIED OBSTRUCTION: the full sorted width distribution; the
  exact minimal k* at 90% (measured above; current pre-read says
  ~36 of 45 coordinates with a nearly uniform 0.66–4.64 spread and
  40/45 coordinates each carrying >= 1% — a healthy concentration
  number will be frozen as the obstruction's first line); the implied
  subbox count to force width/signal < 1 assuming width ~ linear in
  box diameter: ~2^(ceil(log2(73.86)) * k*) subboxes with the milder
  per-coordinate factor (log2 74 ≈ 6.2, so ~2^(6.2·k*) — at k*=36 this
  is ~2^223, and even the mildest honest reading, width only in the
  top half, is ~2^112); and the explicit sentence: branch-and-bound
  over D at this enclosure technology is INFEASIBLE AT THIS COST — a
  first-class frozen negative, not a failure of the loop.
- GATE Z (control): recompute the parent's width arrays from the SAME
  live modules in this campaign's process; require
  max|live − frozen| over the 45 lo/hi arrays <= 1e-12 absolute
  (float-extraction replay noise: the frozen arrays are binary64
  renders of Arb endpoints computed under a different process image;
  the 2026-09-01T02:40 attempt already showed ±1e-14-scale last-digit
  drift on re-run while every DECISION quantity — L1 sum to 1e-9,
  ordering, concentration — is unaffected). The distribution below
  consumes the FRESH live arrays bit-exactly (single source of truth
  for this campaign); the frozen arrays are the cross-run anchor for
  the tolerance above. A worst deviation > 1e-12, or an L1-sum
  disagreement > 1e-6, stops the campaign before any verdict
  (regression finding).

The measured facts of STAGE 1 are checked-in regardless of branch; the
verdict is computed from them mechanically. If the CONCENTRATED branch
fires, STAGE 2's additional pre-registration (calibration subbox cost,
depth cap already fixed, per-subbox branch pinning, three-way outcome
CERTIFY / FAILURE-TO-CERTIFY-at-cap / INFEASIBLE-AT-COST) is appended
to this campaign as pre_statement_addendum_stage2.md BEFORE its first
subbox compute; if the SPREAD branch fires, stage 2 never opens.

## 3. Scope discipline (verbatim per steering)

Any success is a statement about D = {delta : A·delta = 0 exactly,
|delta|_inf <= 1e-7} on the region-0 glob dist block of the VXXZ24
K100_2.37155181 released vector and NOTHING larger. The 12.808x
guardrail stands (gap to published 2.0258350805768544e-6 vs best
available signal 1.5816497000997742e-7): no outcome of this campaign
can bear on the published 2.37155181 or the record 2.371177, and
neither may be mentioned as reachable. The exact rational kernel basis
(rank 24, nullspace 21, integer-asserted in-run; NOT a float SVD
render) is carried unchanged from the parent if stage 2 opens. Reuse
interval_core v11 (ce70f959…) and the existing slope pass (c6b8cef2…);
build nothing that already works. No paper or record claim anywhere;
gate C stays PARTIAL; the six-probe negative is not a theorem unless
stage 2 CERTIFYs.

## 4. Budget and protocol

Stage 1: one nice -n 10 process, threads=1, MID=300, __debug__ on;
reuse of the parent checkpoint for 1.1 (no compute) plus ONE fresh
point+box pass family only if a term-level decomposition requires it
(worst case ~6 min). Freeze protocol: per-attempt stderr under DISTINCT
filenames (adopted 2026-09-01 rule); checksums at freeze; name = actual
UTC. Abort of any attempt writes its traceback to its own file.
