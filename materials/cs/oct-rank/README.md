# `oct-rank/` — real tensor rank of octonion multiplication

**Status: BENCHMARK (no campaign run yet).** A 1977-vintage window was narrowed
to seven values twelve days ago. Closing it further is a bounded,
certificate-producing job.

## The claim

Owner-read (**V-owner**, 2026-08-29) — arXiv:**2608.16649**, Hardik Jain,
*"Bounds on the real tensor rank of octonion multiplication"* (2026-08-17):

* tensor rank of algebra multiplication is $3$ for $\mathbb C$ and $8$ for
  $\mathbb H$ (both classical); for the octonions only a range was known:
  $\ge15$ (Fiduccia & Zalcstein, 1977) and $\le30$ (Cariow & Cariowa);
* **new:** $18\le \mathrm R_{\mathbb R}(T_{\mathbb O})\le25$;
* the lower bound *"peels the eight slices of $T_{\mathbb O}$ down to two and
  bounds the rank of the surviving pencil through the octonion norm"*, and
  generalizes to $\mathrm R_{\mathbb R}(T_A)\ge\frac52 n-2$ for every real
  normed division algebra $A$ of even dimension $n$ — sharp for $\mathbb C$
  ($n{=}2\Rightarrow3$) and $\mathbb H$ ($n{=}4\Rightarrow8$), giving $18$ at
  $n{=}8$;
* the upper bound is *"an explicit rank-$25$ decomposition certified by a
  Krawczyk argument, in exact rational arithmetic, to sit within $10^{-6}$ of
  an exact one"*;
* the same arguments give $\mathrm R_{\mathbb R}(\tau)=7$ exactly for a smaller
  three-slice quaternion tensor $\tau$;
* *"The Lean 4 kernel checks the lower bounds and the Krawczyk existence
  principle; the accompanying scripts check the certificate's finitely many
  exact-rational inequalities."*

Public code: `github.com/hxrdxkxvxd/octonion-rank` (**V-scout**, owner fetch
pending).

## Soundness note — a correction carried from the scan

The `MetaComplexity` scout proposed certifying non-$r$-decomposability by
"exhaustive enumeration of integer witness candidates". **That is unsound and is
not adopted.** Real tensor rank is a condition over $\mathbb R$ — semialgebraic,
not enumerable — so no search over integer or rational witnesses can establish
a *lower* bound. Consequently:

* the **upper** side is witness-producing and certifiable here (find a
  decomposition, prove existence by interval-Newton/Krawczyk);
* the **lower** side is algebra, and this target's contribution to it is
  *auditing* Jain's $\frac52n-2$ argument and its Lean artifacts, not searching.

A failed rank-$24$ search yields **no** claim whatsoever. This is stated in the
gates so it cannot be quietly upgraded later.

## Gates

1. **Gate A (re-verify the published certificate, hours).**
   Independently re-check the finitely many exact-rational inequalities of the
   rank-25 Krawczyk certificate in `python-flint` (exact rationals plus Arb
   balls), and independently re-derive the $10^{-6}$ Krawczyk radius claim.
   Separately, replay the Lean 4 artifacts for the $\frac52n-2$ lower bound and
   the Krawczyk existence principle, recording Lean's version and the exact
   theorem statements accepted.
   **Pass:** every inequality re-verified and the Lean kernel accepts;
   $18\le\mathrm R_{\mathbb R}(T_{\mathbb O})\le25$ becomes
   `MACHINE-VERIFIED` **as a re-verification**, `CITED-DEPENDENCY` for the
   mathematics. **Refutation:** any inequality that fails outward-rounded
   re-evaluation, or a Lean statement weaker than the paper's prose claim.
2. **Gate B (attack the upper bound, days).**
   Numerically search for a rank-$r$ real decomposition of $T_{\mathbb O}$ for
   $r=24$, then $23,\dots$: alternating least squares / Levenberg–Marquardt from
   many seeded starts on the $8\times8\times8$ structure constants, followed by
   **interval-Newton/Krawczyk certification** of an exact decomposition near the
   numerical one — the same method that established $25$.
   **Pass:** a Krawczyk-certified rank-$24$ (or lower) decomposition ⇒ an
   improved upper bound on a 1977-vintage question, shipped with its
   certificate and an independent replay script.
   **No-claim outcome:** search exhaustion produces *no* statement about
   $\mathrm R_{\mathbb R}(T_{\mathbb O})>24$. Record seeds, budget and the
   negative as `COMPUTATIONAL-EVIDENCE` only.
   Budget: fixed before launch — seed count, iteration cap, wall-clock cap.
3. **Gate C ($\tau$ and the generalization, days).**
   Re-verify $\mathrm R_{\mathbb R}(\tau)=7$ exactly (both directions), then
   test the $\frac52n-2$ bound's sharpness claim at $n=2,4$ by exact
   computation, and check whether the peeling argument admits a strictly better
   constant at $n=8$ specifically. Any improvement of $18$ is an algebra result
   requiring human audit before any claim.

## Pre-campaign requirements

- Owner first-hand read of the full arXiv:2608.16649: the peeling argument, the
  exact nondegeneracy hypotheses of the Krawczyk step, and the precise content
  of each Lean theorem.
- Fetch and inventory `github.com/hxrdxkxvxd/octonion-rank`: what the scripts
  check, what Lean checks, and what is assumed rather than proved.
- Owner read (or explicit `[REPORTED]` tag) for Fiduccia & Zalcstein 1977 and
  Cariow & Cariowa. Gates A–C must not depend on either.
- One-page pre-statement per gate: precision, the Krawczyk hypotheses to be
  verified, the search budget, and the pass/fail criteria above.

## Disjointness

No `../math/` or `../physics/` target concerns tensor rank; `../math/qec/` is
algebraic but on bivariate-bicycle qLDPC codes. In-repo neighbour `../mm3/`
shares exact-arithmetic tensor machinery — coordinate on utilities only, the
questions are distinct.

## Layout

- `src/`, `campaigns/`, `scratch/` — created at first use.

## Current state (agent OctRank, 2026-08-30)

Nothing above this section was modified by me — gate ladder, claim,
pre-campaign requirements, disjointness, layout are verbatim from the
owner's contract.

### Verdicts

| gate | verdict | evidence label |
|------|---------|----------------|
| A (rank-25 cert re-verify + Lean) | PASS — all exact-rational inequalities re-checked in python-flint (fmpq + outward arb); K, Banach, per-equation margins all positive; Lean 4 replay accepts lower bounds + Krawczyk principle; axiom audit recorded | MACHINE-VERIFIED (re-verification); math CITED-DEPENDENCY arXiv:2608.16649 |
| B (rank-24 attack) | FAILED-AT-24 — 256 seeded random starts (best rel 2.63e-03, norms O(30), all 300 LM rounds consumed) + 17 structured probes (compress-merge best 1.76e-02; alphabet snap rel ~3; continuation term-0 25-col rows re-converge to machine precision w/ norm explosion then reset; true-24 endpoint degrades 0.563→0.519) + 5-term appendix (16 rows, same clean negative). NO sub-1e-13 rank-24 candidate escaped escalation (one s=0.9 event escalated and correctly re-interpret as re-convergence) | COMPUTATIONAL-EVIDENCE about the declared search only; NO statement about R_R(T_O) ≥ 24 |
| C (tau + sharpness) | PASS — R_R(tau) <= 7 exact Krawczyk replay: K = 0.005351824, worst margin 9.946e-06 (rho = 1e-5), outward arb confirm; upper witnesses: T_C rank-3 Gauss identity, T_H rank-8 Hadamard+sparse identity, both exact-fmpq entrywise | MACHINE-VERIFIED (upper + identities); (5/2)n−2 sharpness at n=2,4 :: equality conditional on the CITED lower bounds (pencil Thm 2), Lean replay side |

The window 18 ≤ R_R(T_O) ≤ 25 remains OPEN at both ends (gate B negative
changes nothing on the window itself).

### Frozen campaigns

- `campaigns/2026-08-30T00:04:39ZZ_gateA_A6A4B0CB/` — gate A full record
  (VERDICT.md, code hashes, gate_a_verify_run.log, lean_replay/).
- `campaigns/2026-08-30T01:05:00ZZ_gateB_control/` — r=25 anchor/control
  (baseline for the ladder contrast) + first gate B scripts.
- `campaigns/2026-08-30T01:40:00ZZ_gateC_tau_sharpness/` — gate C record
  (tau upper, sharpness witnesses, flattening scaffolding + scope
  disclaimers, hashes).
- `campaigns/2026-08-30T01:55:00ZZ_gateB_ladder/` — gate B closure record
  (256-row random sweep, structured probes, continuation curves, appendix,
  hashes, full VERDICT.md).

### Personas / responsible code (all frozen under campaigns/)

- `src/gate_a_verify.py`, `src/gate_a_radius.py` — gate A exact machinery.
- `src/gate_b_search3.py` — random sweep engine (the declared budget).
- `src/gate_b_deflate.py`, `src/gate_b_continuation.py` — structured
  probes (owner-directed redirection).
- `src/gate_c_tau_upper.py`, `src/gate_c_sharpness_n2_n4.py`,
  `src/gate_c_lower_rank4.py` — gate C scripts (see the gate C campaign
  VERDICT.md for scope honesty notes on what is and is not a CP bound).

### How to run

    PY=/Users/jinleic/jinleic-workspace/cs/.venv/bin/python
    cd oct-rank/src
    $PY gate_a_verify.py          # gate A exact replay (flight check)
    $PY gate_a_radius.py          # independent bisection radius re-derivation
    $PY gate_b_ladder_report.py   # regenerate gate B summary table
    $PY gate_c_tau_upper.py       # tau rank-7 upper (exact Krawczyk replay)
    $PY gate_c_sharpness_n2_n4.py # sharpness witnesses (exact)
    $PY gate_c_lower_rank4.py     # flattening scaffolding (no CP bound)

Every script deterministic; environment: python 3.14.3, python-flint 0.9.0,
numpy 2.5.2, scipy 1.18.1, single-threaded BLAS.

### Honest gaps

- Gate B: no statement about R_R(T_O) >= 24; the ladder's negative
  applies to THIS search under THIS budget only.
- Gate C4 (`C3` in pre-statement): n=8 constant probe (could peeling stop
  at 5 or 6 slices, or pencil exceed 12 at n=8?) untouched — remains
  HUMAN-AUDIT-PENDING scratch territory only. Nothing filed.
- Fiduccia-Zalcstein 1977 and Cariow & Cariowa remain [REPORTED] (not
  read first-hand); no gate depended on them.


## Current state (agent OctRankGateC, 2026-08-30)

Append to the section above; nothing in it is edited. Artifact:
`campaigns/2026-08-30T12:15:25Z_gateC4_c598b05f-5e6d-4293-8871-37be3d5fff87_0481d816063f/`
(pre-registered plan in `pre_statement.md`, "Gate C addendum C4/C3",
committed before the probe ran).

### Verdicts

| Item | Statement | Evidence label |
|---|---|---|
| C1 tau lower | every polynomial identity the tau chain consumes re-derived exactly (sympy Z[...] expansion + fmpq probe grid, no floats); chain skeleton = Lean replay [CITED-DEPENDENCY] | MACHINE-VERIFIED (identities) |
| C1 tau | R_R(tau) = 7: upper Krawczyk replay re-run PASS (exact K < 1, margins 9.946e-06 > 0); lower peel(1)+pencil(6) via Lean [CITED-DEPENDENCY] + exact identities | MACHINE-VERIFIED (both replays) |
| C2 n=2 | R_R(T_C) = 3, now fully self-contained: exact 3-term witness (prior campaign) + OWN exact lower proof (det multiplicativity => x^2+y^2 = W(prod of linears) => (W f10 f21)^2 = -1, impossible; Sturm zero-root certificate for t^2+1, exact) | MACHINE-VERIFIED (both directions, no external lower dep) |
| C2 n=4 | R_R(T_H) = 8: exact 8-term witness + lower peel(2)+pencil(6) [CITED-DEPENDENCY: Lean] | MACHINE-VERIFIED (upper); lower CITED |
| C2 quote | paper, verbatim: "The first two are the exact ranks [7, 11], so on the classical cases the bound is sharp." (abstract: sharp for C and H) | REPORTED (read first-hand from arXiv:2608.16649v1 HTML) |
| C4 S1 | pencilRank(1, C) = 12 EXACTLY on the class the chain consumes: exact 12-term witness for (I_8, J_8) entrywise in fmpq (16/16 block checks) + Lean floor 12 [CITED-DEPENDENCY]; J_8^2 = -I entrywise | MACHINE-VERIFIED (witness); conjugacy R-similar J ~ J_8 stated (standard canonical form) |
| C4 S2 | stopping points: peel j (0..6) + pencil floor gives exactly 18 for every j (all seven values + symbolic identity) | MACHINE-VERIFIED |
| C4 S4 | key_bound inequalities EXACTLY TIGHT on the 12-term witness: UD=0, UA=I, UB=J entrywise; rank D = 4 = n/2; rank U = 8; dim ker U = 4 = r - n | MACHINE-VERIFIED |
| C4 verdict | **NEGATIVE**: no parametric choice inside the substitution-peeling + Thm-2-pencil machinery lifts 18. Rule 7: covers ONLY that framework (6 peels + pencil floor 12, the exact Oct18.lean machinery); NOT swept: non-peeling lower-bound routes for T_O, any rank bound on 3-slice L-families (L_u, L_v, L_w) beyond the peel-implied 13, non-similar normalizations | MACHINE-VERIFIED (all three probes) |

### Honest gaps (update, within this appended section only)

- Gate C4 (the "untouched" item above) is NOW FILED as a negative with
  rule-7 scope — see the C4 rows above and the frozen VERDICT.md.
- S3 remains OPEN and was NOT swept: is rank(L_u, L_v, L_w) >= 14 for
  independent u, v, w in R^8? Machinery gives 5 + 13 = 18; a >= 14 bound
  would improve 18 to 19. Anything here is HUMAN-AUDIT-PENDING scratch
  territory per pre-statement C3.
- Gate B honest gap and the Fiduccia-Zalcstein / Cariow status above are
  unchanged.

Run (from cs/oct-rank/src, single-threaded):

    PY=/Users/jinleic/jinleic-workspace/cs/.venv/bin/python
    $PY gate_c_tau_lower_exact.py       # tau lower: exact identities (PASS)
    $PY gate_c_sharpness_lower_exact.py # n=2 self-contained lower (PASS)
    $PY gate_c4_peel_n8.py              # C4 probe (NEGATIVE, exit 0)

Environment: python 3.14.3, flint 0.9.0, sympy 1.14.0; no floats in the
three new scripts; deterministic.


## Current state (agent OctRankS3, 2026-08-30)

Append to the sections above; nothing in them is edited. Artifact:
`campaigns/2026-08-30T19:27:14Z_d6ca8a53-0950-4a3e-a6ef-084d0a1f39d0_S3/`
(`pre_statement.md` committed BEFORE the first computation; scripts and
outputs byte-frozen with SHA256 in `code_hashes.sha256`; verdict with
rule-7 scope in `VERDICT.md`). This closes the S3 honest-gap item as an
exact NEGATIVE-with-obstruction. Per pre-statement C3 the campaign-level
candidate notes live in
`scratch/s3_three_slice_notes_HUMAN-AUDIT-PENDING.md` only.

### Verdicts

| Item | Statement | Evidence label |
|---|---|---|
| Anchors | n=2⇒3, n=4⇒8 reproduced exactly (upper witnesses entrywise fmpq: 3-term Gauss identity, frozen gate-C 8-term; lower-chain identities exact); tau = 7 CONTROL (the rule-14 pin: any S3 machinery giving ≥ 8 at (n=4, k=3) is refuted) | MACHINE-VERIFIED (all three) |
| S3 outcome | **NEGATIVE**: rank(L_u, L_v, L_w) ≥ 14 NOT established. Best certified bound remains the chain's 13 (1 + pencil floor 12, C4-verified tight); gap to 14 = 1. No claim about R_R(T_O); 18 ≤ R_R(T_O) ≤ 25 unchanged | MACHINE-VERIFIED (route ceilings); universal floor statement unchanged |
| Route A (Strassen commutator) | rank([L_ū L_v, L_ū L_w]) = 8 exactly on 173/173 swept triples (56 basis + 112 structured + 5 fixed); twin readings 8+8=16 or 8+4=12. NOT a certified floor: sweep ≠ universal quantifier; the naked twin overshoots the tau CONTROL (would say 8 for the n=4 3-slice tensor whose exact rank is 7); Strassen/Bläser/Lickteig 3-slice form not read first-hand [CITED-DEPENDENCY UNVERIFIED] | MACHINE-VERIFIED (per-triple ranks); universal 16 statement: COMPUTATIONAL-EVIDENCE only, NOT a bound |
| Route B (AFT pivots) | pivot residual independence verified exactly (×3 pivots ×56 triples); all bookkeepings = 13; double-peel 10; flattening 8 | MACHINE-VERIFIED (route ceiling 13) |
| Route C (residual-pencil floor) | residual pairs independent + irreducible quadratic C² − 2aC + bI = 0, a² < b verified exactly per triple; with C4's pencil floor = 12 exactly, the substitution route cannot exceed 1 + 12 = 13 | MACHINE-VERIFIED |
| Route D (flattenings) | flattening ranks 3 / 8 / 8 exactly on all 56 basis triples (capped at 8 ≪ 14: a miss, as pre-registered). First-page retraction: an earlier mode-(c) slicing bug reported 1 for e₀-triples; fixed script + rerun frozen | MACHINE-VERIFIED |
| Route E (exact infeasibility) | pre-registered BUDGET-FAIL (312-unknown Gröbner out of reach); only a necessary-only linear-feasibility diagnostic ran; NOT a bound | BUDGET-FAIL (as pre-registered) |
| Route F (post-hoc upper side) | EXACT: span{e₀..e₃} is a quaternion subalgebra; for x ∈ H′, T L_x T⁻¹ is block-diag(Lq_x, conj-equiv block) with T = diag(1,..,1,−1); the specific triple (1, i, j) satisfies rank ≤ 14 (7+7 block-concatenation of two tau copies) and ≥ 13 (chain) — triple sits in [13, 14]; worst case NOT settled (upper on one triple says nothing about the universal floor) | MACHINE-VERIFIED (block facts, exact); 14 upper is specific-triple only |
| Cert-digit landmine | tau_r7/rank25 certificate decimals are 19-significant-digit literals; float64 parse destroys ~4 digits (per-entry ~5e-19, post-contraction ~1.9e-16 measured exactly). Gate A/C replays unaffected only because their margins ≫ 1e-6. Anchor semantic fixed: certificate ⇒ exact witness within ρ, not entrywise-zero stored decimals | MACHINE-VERIFIED (exact diffs measured) |

### Rule-7 scope

Covered: the S3 question (3-slice L-family floor at n=8) over ALL
independent (u, v, w) ∈ ℝ⁸ via five pre-registered routes + one post-hoc
upper-side route; 173 triples swept for Route A; exact (fmpq) arithmetic
throughout, no floats in load-bearing steps. NOT swept: a universal
commutator-rank ≥ 8 statement (sweep only — rule 7/17a), first-hand
Strassen/Bläser/Lickteig form verification, non-L-family residuals,
other n, the (1, i, j) 13-or-14 decision, and anything touching
R_R(T_O) itself.

### Honest gaps (update, appended section only)

- S3: OPEN, now with a filed obstruction profile (Route A is the only
  live candidate; needs Strassen-form first-hand verification + a rank
  argument replacing the sweep). Best certified floor 13; gap = 1.
- The (1, i, j) triple: rank ∈ [13, 14] — deciding requires an exact
  14-term witness (13 would need rank-13 infeasibility): both open.

Run (from the campaign directory, single-threaded):

    PY=/Users/jinleic/jinleic-workspace/cs/.venv/bin/python
    cd oct-rank/campaigns/2026-08-30T19:27:14Z_d6ca8a53-0950-4a3e-a6ef-084d0a1f39d0_S3
    for s in s3_anchor s3_routeAB s3_routeC_substitution s3_routeB_AFT \
             s3_routeD_flattening s3_routeE_groebner s3_routeF_subalgebra \
             s3_routeF2_blocks s3_routeF3_global; do $PY $s.py; done
    # all 9 scripts exit 0; deterministic; python 3.14.3, flint 0.9.0

Two retractions from this campaign are recorded inline in VERDICT.md
(Route-A first-draft bookkeeping; Route-D first-version flattening bug)
with the superseded artifacts kept byte-frozen alongside the corrected
ones.
