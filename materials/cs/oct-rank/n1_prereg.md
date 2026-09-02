# N1 pre-statement — exact CP-completion / downward certification attack on rank((L_1, L_i, L_j))

Campaign: to be minted by `scripts/campaign.py init --gate n1_exact_cp_completion
--prereg cs/oct-rank/n1_prereg.md` AFTER this file is committed (path-scoped
commit, per Main's standing order of 2026-09-01). Owner agent: `OctRankNext`.
No N1-parameter-dependent computation has been run. Disclosed precommit
environment probes (Main-ordered disclosure; NO campaign value, all in /tmp,
nothing used by any N1 script): (a) msolve-0.8.0 source downloaded and its
configure inspected; (b) flint-3.6.0 built once in /tmp (make exit 0),
abandoned. The exact backends below use ONLY the existing venv
(python-flint 0.9.0 fmpq/arb, numpy 2.5.2, scipy 1.18.1, sympy 1.14.0).

## 0. Standing state this campaign acts on (all re-read first-hand this session)

Object: real tensor rank of the 3-slice octonion tensor
T = (L_1, L_i, L_j) on Cayley-Dickson octonions, basis (1,i,j,k,l,il,jl,kl),
product (a,b)(c,d) = (ac - conj(d)b, da + b conj(c)), L_x[c,b] = (x*e_b)_c.
Certification target (Route-F convention, machine-verified rank-equivalent by
two entrywise diagonal identities in the frozen rf_krawczyk.py anchors):
TF = blockdiag(tau, tau), tau = quaternion 3-slice tensor (I, L_i, L_j).

Frozen certified facts being consumed [MACHINE-VERIFIED, replayable]:
- FLOOR: rank(T) >= 13 (1 peel + pencil 12; s3/C4/routeAF frozen chain;
  Route A prov campaign closed the only live >13 route at 12 < 13).
- UPPER: rank(T) <= 14 (blockwise 7+7 via the frozen tau_r7 certificate).
- Hence rank(T) in {13, 14}: OPEN. Deciding = finding a rank-13 witness
  (=> =13) or a certified rank >= 14 argument (N2's job; not this campaign).
- Route F (2026-09-01T04:30:00Z) FAILURE-TO-CERTIFY: the two frozen
  numerical candidates, polished to rel 8.816619e-06 / 5.694968e-05 under the
  declared budget, are CERTIFIED-EXCLUDED: complete no-root-in-box statements
  at rho <= 1e-6 (exactly reproduced this session from frozen bytes; replay
  exit 0). They are NOT seeds for N1 — they sit in certified empty boxes.

## 1. Fixed conventions

- Rank attacked: r = 13. Field/ring: real numbers; every load-bearing number
  is an exact fmpq (dyadic parses of float64 candidates, rational Newton
  iterates); float64 is used only for (i) TRF polishing, (ii) the QRCP slice
  choice, (iii) diagnostics. Final adjudication numbers are exact
  (fmpq; outbound comparisons re-verified in outward-rounded arb, midpoint/
  radius form — the rf_krawczyk.py discipline, reused verbatim).
- Normalization: the CP torus gauge is handled by the frozen geometric-mean
  column balancing (`balance`) between polish rounds (gauge move; image
  preserved in exact arithmetic) — identical to Route F.
- Symmetry handling: none imposed. The exact block structure of TF
  (blockdiag(tau,tau), 96 exactly-zero cross-block entries among 192) is
  recorded as a descriptive anchor only; no ansatz forces it on the solver
  (structured seeds below are INITIALIZATIONS, not constraints).
- Search cap (fixed NOW): 40 seeds = 32 random + 7 declared 14->13 merge
  seeds + 1 frozen-excluded restart; 3 TRF rounds x 5000 nfev max per seed
  (identical solver settings to the frozen instrument: method='trf',
  tr_solver='exact', x_scale='jac', ftol=xtol=gtol=3e-15); hard CPU cap
  2 h (time.process_time, in-process accounting only, rule 17e); nice -n 10;
  single-threaded BLAS pinned in-process before numpy import.

## 2. Candidate generation (the NEW domain; all seeds fixed before commit)

Class R (32 random): numpy default_rng, base seed 20260901; seed k uses
PCG64(seed=20260901 + k), k = 0..31; each start = factors drawn
N(0,1)/sqrt(8) scaled by ||TF||_F^{1/3} (the frozen route_af_search scaling;
||TF||_F = 4.898979485566356).
Class M (7 merge seeds): from the frozen tau_r7 certificate factors
(certs/tower_cert/tau_r7/{A,B,C}.csv, parsed exactly as dyadics, block-
duplicated to the certified 14-term decomposition of TF) — mechanical merge
rule fixed now: for s in {0,1,2,3}: one 13-column layout where top-term s and
bottom-term s share column 12 of 13 (b = (B7[:,s]; B7[:,s]) stacked, c
likewise; a = A7[:,s]); the remaining 12 columns are the un-merged pairs.
For s in {0,1,2}: analogous merge of top-term (s+4) with bottom-term (s+4)
... concretely: seed M_j for j = 0..6 merges tau-column j with tau-column
(j+7) mod 14 of the block-duplicated 14-term layout. These seeds encode the
hypothesis that a 13-term solution mixes blocks near the certified 14-term
one. (The merge is an initialization only: the merged column's single
3-vector a generally fits NEITHER block term it replaces; TRF does the rest.)
Class F (1 seed): the frozen polished candidate (candidate_main_factors.csv
of 2026-09-01T04:30:00Z_routeF_kraw) — it sits in a certified-empty box; its
3 declared rounds are a recorded probe of whether a DEEPER root exists in the
same basin (the certified exclusion covers only rho in [1e-6, 1e-12] boxes;
larger neighbourhoods were never excluded).

## 3. Exact verification protocol (fixed NOW)

Stage E1 (per seed, after 3 TRF rounds): parse the polished factor vector to
exact dyadics (Fraction(float) semantics); compute the EXACT fmpq residual
g0 = CP13(x) - TF entrywise; record ||g0||_max exactly.
Stage E2 (admission): if a seed reaches ||g0||_max <= 1e-7 (exact fmpq
comparison), the seed is ADMITTED to exact iteration:
  E2a exact Newton (max 3 iterations): Y = J_S^{-1} exact (fmpq_mat.inv,
  slice rule preserved FREEZE-TO-FREEZE: qrcp with natural fallback as
  frozen), x <- x - Y g0 (exact; Krawczyk displacement without remainder —
  candidate-only, NOT a certification), re-parse, re-slice, re-invert each
  iteration. Stop early when g == 0 entrywise EXACTLY.
  E2b if any exact iterate x* has CP13(x*) == TF entrywise (fmpq, all 192):
  x* IS an exact rank-13 witness. Verdict: rank(T) = 13 EXACTLY
  [MACHINE-VERIFIED]; the witness bytes are frozen; independent replay script
  re-verifies entrywise from frozen bytes. Krawczyk is then NOT needed
  (witness is exact, not proximate) — recorded as such.
Stage E3 (proximate certification): if no exact-zero iterate is reached but
some admitted iterate satisfies the CONTainment window at any fixed rung of
the frozen 15-rung ladder (1e-2 .. 1e-12, same margins: |c0_c| + R2_c < rho
strict, exact fmpq, outward-arb re-verified), then containment PASS =
certified exact root on the slice = rank 13 [MACHINE-VERIFIED] (the
gate-A/Route-F margin algebra, reused verbatim).
Stage E4 (kill): if no seed is admitted (all ||g0||_max > 1e-7), or every
admitted seed's Newton diverges/bounces without exact zero AND the window
never opens at any rung for any admitted candidate: verdict FAILURE TO
CERTIFY. Full residual distribution, admission counts, and (for admitted
seeds) per-iteration ||g0||_max histories are frozen. NO rank claim in
either direction; {13,14} stays OPEN; N2 proceeds (unless the assignment's
logical-implication clause fires — it cannot here, see section 6).

## 4. Controls (all through the SAME instrument path, run BEFORE the main
sweep; behavior contract identical to the frozen instrument)

- C-TAU7 (accept known-valid decomposition): tau at r=7 with the frozen
  certificate factors through the identical square-slice path MUST certify
  containment (frozen result: at rho=1e-3). If it fails: instrument broken;
  campaign verdict INVALID INSTRUMENT; main sweep not run.
- C-TAU6 (reject globally-infeasible plant): tau at r=6 (true rank 7,
  machine-verified) MUST NOT contain at any rung; certified exclusion
  expected. Containment here = FALSE CERTIFICATE = quarantine; verdict
  INVALID INSTRUMENT, nothing else reported.
- P-POS (accept synthetic known-true): synthetic rank-13 tensor (exact dyadic
  factors, numpy default_rng(20260901), scale 0.7) + S-only 1e-8 plant
  perturbation (default_rng(20260902)), identical to the frozen R2/R3
  construction; MUST certify.
- P-WRNG (reject one-coordinate perturbation of the target — the assignment's
  named near-miss control): corrupted table T_cor = TF with entry [0,0,0]=2,
  probed with the best N1 candidate at the same slice; containment MUST be
  absent at every rung (a certified N1 root of TF has residual exactly 1 in
  coordinate [0,0,0] of T_cor). Containment firing = instrument broken =
  quarantine + INVALID INSTRUMENT.
- NOT re-run this campaign: P-POS parameters, ladder, and slice rule are
  byte-frozen from Route F (no re-registration drift); any deviation found in
  the copied instrument code (SHA256-pinned) aborts initialization.

## 5. Budget, environment, accounting

- Compute cap: 2 h CPU total (in-process time.process_time printed by every
  script); expected ~40-70 min (polish dominates). Exact steps: up to 4
  admitted seeds x 3 Newton iterations x one 192x192 fmpq inverse (~30 s
  measured, frozen datum) per iteration <= ~6 min; ladder only for admitted
  non-exact candidates (1 certification pass each, <= 1 min each).
- Environment: /Users/jinleic/jinleic-workspace/cs/.venv/bin/python
  (3.14.3, flint 0.9.0, numpy 2.5.2, scipy 1.18.1, sympy 1.14.0);
  OMP/OPENBLAS/MKL/VECLIB/NUMEXPR threads = 1 pinned in-process before numpy
  import; PYTHONDONTWRITEBYTECODE=1 (pristine-freeze compatible).
- Determinism: seeds, merge rule, round counts, admission thresholds,
  slice rule, ladder are all fixed above; every script exits 0 or a
  documented nonzero; no daemonized processes.

## 6. Relation to N2 (documented, per assignment item 2)

- If N1 produces a certified rank-13 witness (exact or containment): the
  question "is there a 13-decomposition" is answered YES, refuting the
  infeasibility hypothesis that N2 would attack; N2 becomes logically
  unnecessary for the {13,14} decision (rank is then =13) and will NOT be
  opened; the implication (witness => infeasibility question void) is this
  clause, recorded here before any compute.
- If N1 closes FAILURE-TO-CERTIFY: no implication either way (absence of a
  witness is NOT evidence about feasibility); N2 opens as a distinct
  campaign under its own committed pre-statement.

## 7. Adjudication vocabulary (fixed now)

- "rank(T) = 13 EXACTLY" [MACHINE-VERIFIED]: only via E2b exact-witness or
  E3 containment PASS (floor 13 inherited machine-verified; replay check
  part of every admission path).
- "FAILURE TO CERTIFY": everything else, with the section-3 kill data. The
  two frozen sentences (verbatim, in every failure wording): "Absence of a
  13-witness is NOT evidence for 14; absence of an impossibility argument is
  NOT evidence for 13."
- "INVALID INSTRUMENT": on any control break; campaign halts, main data not
  reported as evidence.
- "CRASHED": script failure left the run incomplete (status.json verdict
  FROZEN-INCONCLUSIVE with the failure history).
- Published window 18 <= R_R(T_O) <= 25 is untouched by every outcome of
  this campaign; nothing here bears on it; no published work is refuted.

## 8. Rule-7 prospective scope

Covered: the fixed conjugated Route-F tensor (blockdiag(tau,tau), proven
rank-equivalent to (L_1,L_i,L_j) by the two frozen entrywise diagonal
identities) at rank 13, from the 40 declared seeds through the declared
3-round TRF polish, the declared exact-Newton admission/iteration protocol,
one declared slice rule per system (frozen qrcp/natural-fallback), the frozen
15-rung ladder for admitted proximate candidates, and the four declared
controls; exact fmpq arithmetic on every load-bearing number; outward-arb
re-verification of every passing comparison; in-process CPU accounting only.
NOT covered: any other seed, round count, nfev cap, admission threshold, box
radius, or slice scheme; rank-14 statements in either direction; the full
8-slice T_O; border rank; complex decompositions; other octonion triples; the
commuting-extension Koiran model (untouched this campaign after the frozen
EXT evidence); Gröbner/elimination machinery (N2's domain, separate prereg).

## 9. Lifecycle and freeze plan

Sequence (Main's strict order): commit this file (path-scoped; no other
agent's paths) -> `campaign.py init --gate n1_exact_cp_completion --prereg
cs/oct-rank/n1_prereg.md` from the workspace root -> copy this file
byte-identically into the minted run dir (campaign.py hashes the file but
does not copy it) and record the source commit hash + sha256 in
run-dir PROVENANCE.md -> controls (section 4) -> main sweep -> freeze
(`campaign.py freeze`) with scripts, outputs, candidate factors, and verdict
notes in-run -> `campaign.py close` with the section-7 verdict mapping.
status.json verdict mapping: E2b/E3 success -> FROZEN-CERTIFIED; E4 with
clean kill data -> FROZEN-INCONCLUSIVE (question stays open; no negative
rank claim exists to freeze); control break -> FROZEN-INCONCLUSIVE with
INVALID INSTRUMENT wording; script crash -> FROZEN-INCONCLUSIVE (CRASHED).
