# PRE-STATEMENT — OmegaBranchSix: certified subdomain pins for the SIX
# remaining R branches on the exact-kernel ray family; agent OmegaPcPen,
# 2026-09-01. Committed to git BEFORE campaign init and BEFORE any
# parameter-dependent compute. Main steering: "one more attempt, then
# omega rests". No rerun of any dead route (B&B, whole-box R_glob[0] pin,
# kernel-ray R_glob[0] pin — all three FROZEN-NEGATIVE).

## 0. Route selection from the frozen ledger

The parent campaign 20260901T132102Z_f17bb4cb_117823c32d87 pinned ONE
branch (R_glob[0], span 3.2123e-5 = 203x signal at the whole box) and
FROZE the obstruction. Its own manifest names the successor: "certify
ordering on a DIFFERENT branch (the other six: six have whole-box spans
<= 2.06e-6 — two orders smaller than R_glob[0]'s)". THIS campaign does
exactly that — one fresh campaign, the last before omega rests.

## 1. The six branches, with their MEASURED whole-box spans (frozen
##    stage1_decomposition.json, campaign 2026-09-01T03:15:00Z…9d2e4a7c)

| branch | whole-box max candidate span | candidates' lows (box, frozen) | argmin |
|---|---|---|---|
| R_comp[0,3] | 3.825927815115371e-07 | 0.24530643888350231 / 0.24530642754901422 / 0.24530642751181178 | 2 |
| R_comp[1,3] | 3.3618547615188454e-07 | 0.24715958702781543 / 0.2471595763033872 / 0.24715957464433005 | 2 |
| R_comp[2,3] | 3.604135526968921e-07 | 0.2456556191252014 / 0.24565561011857182 / 0.24565561130106037 | 1 |
| R_comp[0,2] | 6.170856338139075e-07 | 0.5883306757881364 / 0.5883306699849327 / 0.5883306702158023 | 1 |
| R_glob[1] | 7.622166832144117e-07 | 0.4968521822710942 / 0.49685294448520734 / 0.4968521822710635 | 2 |
| R_glob[2] | 2.0511233057018963e-06 | 0.49684892842286427 / 0.49684892842208755 / 0.496850979549942 | 1 |

(versus R_glob[0]'s 3.2123059392596964e-05, the pinned-and-falsified
branch). Every other parameter block stays at p* points; the swept
subdomain is the SAME certified kernel ray family as the parent campaign:

- s = theta·v with theta the frozen 21-coefficient CandidateWitness
  direction; s in {-1,0,+1}^45; A·s = 0 over Z; sum(s) = 0; basis V,
  margin matrix A, params MAT hard-hash-verified in-run (E0);
- ray = {p* + rho·s : rho in [0, r0)}, r0 = 1e-7 (box-wall-bound, exact
  Fraction 944473296573929/9444732965739290427392; positivity wall
  1.6166317499672092e-05 not binding — parent PFIRST);

## 2. The certified subdomains and inequality tiers

Because the parent's P2 measured the ray moves every candidate hull by
~1.5e-7..1.4e-6 while the CENTER candidate gaps are:
  comp branches: ~1.14e-8 (R_comp[0,3] cand0-vs-cand2 at p*),
  R_glob[1]: 3.08e-14 (cand0 vs cand2 at p*),
  R_glob[2]: 7.77e-13 (cand0 vs cand1 at p*),
full argmin-identity collapse (Tier-1) is POSSIBLE ONLY where the center
gap exceeds the ray hull growth. The campaign therefore pre-registers
THREE tiers, each a certified statement on the ray family, adjudicated
per branch — the verdict language will name the tier actually certified:

- TIER-1 (full pin, the strongest): the argmin identity over the ray —
  for a branch with frozen center argmin t*, certify
  G_j := cand_t*.lo - cand_j.hi > 0 outward at BOTH rho=0 and rho=r0 for
  every j != t*, with the ray endpoints' intervals overlapping nothing
  (then t* is the UNIQUE possible minimizer on the ray family and the
  branch's hull collapses to the single-candidate width).
- TIER-2 (exclusion pin): certify a single weaker inequality per branch
  — cand_EXCLUDE.lo - cand_argmin.hi > 0 outward at both endpoints,
  which removes ONE candidate from the possible-minimizer set on the
  ray family (a genuine certified subdomain statement, weaker than
  Tier-1). Registered exclusion targets (from the frozen box lows, where
  the gapped candidate is far from the argmin):
    R_comp[0,3]: exclude cand_0 vs argmin cand_2 (center gap 1.14e-8);
    R_comp[1,3]: exclude cand_0 vs argmin cand_2 (1.27e-8);
    R_comp[2,3]: exclude cand_0 vs argmin cand_1 (9.01e-9);
    R_comp[0,2]: exclude cand_0 vs argmin cand_1 (5.80e-9);
    R_glob[1]:   exclude cand_1 (the pc candidate) vs cand_2/0 (3.08e-14
                 center gap — EXPECTED TO FAIL; recorded either way);
    R_glob[2]:   exclude cand_2 vs cand_1 (7.78e-13 — expected to fail;
                 recorded either way).
- TIER-3 (sign knowledge): if no gap survives positively even at rho=0,
  the branch records its exact failing inequality at both endpoints with
  the factor versus the signal — a certified no-pin for THIS ray family.

Per-branch verdict = the strongest tier whose inequalities hold outward
at BOTH endpoints e0, e1. The GLOBAL verdict is:
  FROZEN-CERTIFIED iff at least one branch certifies Tier-1 OR Tier-2
  (a genuinely pinned/excluded subdomain — the Main-steering outcome);
  FROZEN-NEGATIVE iff all six branches report Tier-3 (no pin anywhere);
  OPEN only on budget/kill.

## 3. Machinery (pre-declared build)

- The parent's pcpen.py aggregation (dependency-of-record, sha
  f07b70ae…f2c6 — this campaign's runner is a LINE-BY-LINE extension of
  it, hash-pinned in-run like the parent's E0) already computes the
  per-branch candidate intervals for comp branches via `cand =
  [nb[tt].lo - (pc.hi if tt==r else pen.hi)]` at l=3 and via nb2 at
  l=2, and via glob candidate triples for R_glob. Nothing new is built
  except: (a) the per-branch candidate-triple bookkeeping (the parent
  only recorded R_glob[0]'s triple), and (b) the PV-style membership
  proof and outward G_j interval assembly per tier.
- Backend identical: python-flint arb, v11 endpoint-pair Ival, MID 300,
  OUTB 64, exact-Arb dyadic endpoints; buffered .low()/.up() outward
  endpoints on every decision inequality; PYTHONDONTWRITEBYTECODE=1;
  nice -n 15; OMP/OPENBLAS/MKL threads = 1.
- In-run hash verification (E0, reused verbatim thresholds): v11 core
  ce70f959…, vxxz24_float 31657c92…, kernel basis 7d69c1dd…, margin
  matrix 423dbafd…, parent runner f07b70ae…f2c6 (the run-dir bytes),
  plus params prefix df75ae3a.

## 4. Controls (registered, each a stop-gate; both directions covered)

- GATE-Z (reproduction of BOTH frozen records at radius 0): PASS gate,
  tolerances 1e-12 per R branch, 1e-6 on Om_raw — the parent's frozen
  GATE-Z semantics reused verbatim.
- C1 (known-true ACCEPT, plant): the parent's radius-0 candidate-identity
  cand1.lo - cand0.lo = -2*prop_0*eps_0 at 1e-12 (frozen 5.79e-13) — the
  same known-true object, which doubles as the R_glob[0] sanity anchor.
- C3 (sign/branch perturbation REJECT): +1e-3 broadcast dual drift on
  region 1 (NOT region 0 — a fresh plant direction) must move eps_1 by
  exactly 4*drift (<=1e-12) and R_glob[1] must respond; a positive
  R_comp[0,3] Tier-2 gap that does NOT move under this drift would mean
  the machinery is insensitive — recorded as instrument-reject.
- C2 (wrong-index dual REJECT): region-2 lookup must return the frozen
  eps_2 live value, bit-distinct from eps_0 and eps_1 (fresh index — the
  parent's C2 pinned region 1).
- PV (ray-endpoint kernel membership, exact): A·delta=0, sum=0,
  |delta_i|<=r0, positivity — the parent's PV reused verbatim.
- A control FAIL stops the campaign with the failing control named
  verbatim (kill rule).

## 5. Budget, checkpoints, kill logic

- Hard CPU budget: 30 minutes wall, single process, nice -n 15,
  threads=1 (matches the parent's actual 9 s cost — 200x headroom).
- Checkpoint JSON written after EACH control and each per-branch
  P-pass, under per-attempt distinct log filenames (freeze hygiene).
- KILL = OPEN verdict naming the exact failed step (machinery defect —
  NOT a statement about the object); FROZEN-NEGATIVE only for genuine
  certified Tier-3-neutral outcomes.

## 6. Adjudication (fixed now; exactly one terminal verdict)

- FROZEN-CERTIFIED: >=1 branch pins Tier-1 or Tier-2 at BOTH ray
  endpoints with all controls passing. Deliverable: the certified branch
  list with the pin/exclusion inequality intervals and what they imply
  for the 21-dimensional kernel question — a Tier-1/Tier-2 pin means
  those branch hulls on the ray subdomain are collapsed/excludable,
  which SHRINKS the certified-gradient target set for any successor LP/
  ANY-y instrument. The 21-dimensional kernel question is UNCHANGED
  AND OPEN (this campaign decides nothing about it; it trims the
  hull-cost mover list inside it).
- FROZEN-NEGATIVE: every branch is Tier-3; deliverable: the exact
  failing inequality per branch with its factor versus the signal and
  the surviving movers list.
- OPEN: kill/budget, named step.

## 7. Rule-7 sentence (scope)

The measured set of THIS campaign: the candidate-order inequalities of
the SIX branches R_comp[0,3], R_comp[1,3], R_comp[2,3] (level-3 hashing,
pc/pen structure), R_comp[0,2] (level-2 symmetric), R_glob[1], R_glob[2]
(glob hashing) at two points of the certified kernel ray family
{p* + rho·s : rho in [0, 1e-7]} (rho = 0 and the exact dyadic endpoint
r0·s; PV-certified) on the region-0 glob dist[0] block of the VXXZ24
K100_2.37155181 released vector under the transcribed 3-region
single-p_comp program at max_level 3, q = 5, all other parameter blocks
at p* (width-0 asserted in-run). NOT swept: R_glob[0] itself (frozen
falsified), any other radius, any other subdomain (other kernel
directions, the box), any other region's dist blocks, any part-level or
lam blocks as swept sets (C3's region-1 drift and C2's region-2 lookup
are plants, not sweeps), any subdivision, Monte Carlo, LP solver, or
gradient pass. The published 2.37155181 and the record 2.371177 were
not touched (parameters unpublished; 12.808x guardrail survives any
outcome — no record claim anywhere). Gate C (local optimality over D)
stays PARTIAL and the 21-dimensional kernel question stays UNCHANGED
AND OPEN regardless of verdict.

## 8. Exact constant register (asserted in-run; display anchors)

- frozen R detail (with lemma1): the seven branches as in the parent
  campaign pcpen.py (FROZEN_R_DETAIL) — bit-anchored GATE-Z targets;
- frozen CandidateWitness base candidates for R_glob[0] (the C1 anchor);
- frozen eps_0 = 2.1502086942121845e-06, eps_1 = 1.143325024764798e-06;
- frozen branch spans: the six rows of §1 (display anchors, re-derived
  from the live box pass for the tier arithmetic — the certification is
  the live outward intervals, not the frozen floats).
