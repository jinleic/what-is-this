# Route F pre-statement — af-triple13-existence (agent OctRankRouteF)

Campaign: `campaigns/2026-09-01T04:30:00Z_routeF_kraw/` (directory named with
the ACTUAL UTC start of preregistration; the first committed compute starts
only after the git commit of this file). This file is committed BEFORE any
campaign computation. Rules 1-17e of cs/README.md and the oct-rank README
contract were read before writing. Rule 16: every domain, box, radius, cap
and selection rule below is FIXED NOW; no adjustment after seeing results;
every miss is reported.

## 0. Question and standing state (re-derived from frozen artifacts, not summaries)

Decide real tensor rank of the 3-slice tensor T_F = (L_1, L_i, L_j) on the
Cayley-Dickson octonions, basis (1, i, j, k, l, il, jl, kl), product
(a,b)(c,d) = (ac - conj(d) b, da + b conj(c)), L_x[c,b] = (x*e_b)_c
(convention byte-identical to gate A / S3 / Route AF).

Frozen state being acted on [all re-read first-hand this session]:
- rank(T_F) >= 13 [MACHINE-VERIFIED chain: 1 + pencil 12, S3 lower13 chain
  replay + C4-verified pencil floor].
- rank(T_F) <= 14 [MACHINE-VERIFIED modulo frozen tau facts: blockwise
  7 + 7, upper14_replay.out].
- So rank(T_F) in {13, 14}: OPEN. Route A is CLOSED-NEGATIVE at 12 < 13
  (Strassen 1983 primary, VERDICT.md, commit 076a4dd); this complements the
  S3 negative (best certified floor 13, gap to 14 exactly 1; substitution
  route proven unable to exceed 13 by exact per-triple checks).
- A numerical rank-13 search (5 CP + 3 commuting-extension seeds, 5000 evals)
  produced NO witness [COMPUTATIONAL-EVIDENCE about those seeds only]. Two
  frozen candidates are retrievable and are THE named seeds:
    S_CP  = campaigns/2026-08-31T08:02:18Z_routeAF/rank13_candidate.npz
            (rel = 1.3252060344282049e-3, seed 2, random+ALS)
    S_EXT = campaigns/2026-08-31T08:02:18Z_routeAF/rank13_extension_candidate.npz
            (rel = 9.99874493968644e-5, seed 2, commuting extension)
  Both files were SHA256-verified against that campaign's frozen
  checksums.sha256 ledger before any use: S_CP = 3a1be26e2ece...,
  S_EXT = 5d1bec8b63a6.... Conventions: CP x-vector layout is
  (a:3x13, b:8x13, c:8x13) row-major concatenated, T[p,b,c] = coefficient of
  e_c in e_p * e_b for p = 0,1,2; the extension z-vector is (S:13x13 row-major,
  d:13, e:13) with the same target.

## 1. Fixed convention and startup anchor (assignment item 5a)

At startup the octonion multiplication table is built by THREE independent
constructions and must agree AS EXACT INTEGER LISTS, entrywise:
  (i) the S3/RouteAF cd_mul Cayley-Dickson recursion (s3_routeAB.py);
  (ii) the gate-A index-table construction (t4 + h_mul, gate_a_verify.py);
  (iii) the upstream scratch/upstream_ref/verify/octonion_core.py
       octonion_tensor() table.
Additional anchors before any campaign compute: e0 unit entrywise (both L and
R actions); norm multiplicativity at 20 fixed pseudo-random integer points
(fmpq-exact, seeded); T_F has 24 nonzero integer entries, all in {-1,0,1};
slice flattening ranks of T_F are exactly (3, 8, 8); L_{conj(u)} L_u = N(u) I
entrywise on all 8 basis elements (route-A C2 anchor). PASS is a precondition
for everything that follows; any mismatch HALTS the campaign (escalate).

## 2. Instrument substitution (diagnosis accepted by Main 2026-09-01; authorized with 4 conditions, all registered here)

The handoff plan ("20-variable symmetric-real realification; symmetrized
20x20 Jacobian positive-definiteness as the hypothesis") CANNOT establish
rank 13 as sketched:

  (a) A Krawczyk box around a critical point of a LEAST-SQUARES functional
      certifies only that the gradient vanishes there (or that no critical
      point exists in the box). It NEVER certifies that the residual is
      exactly zero: a positive-definite strict local minimum with residual
      ~1e-6 passes every positive-definiteness test. Rank-13 existence needs
      a certified root of the CP equations themselves.
  (b) Interval arithmetic cannot prove a floating residual is exactly zero;
      the certification target must be a square system whose exact solutions
      ARE exact CP decompositions of the target tensor.

Substituted instrument (existence-carrying): Krawczyk/interval-Newton on the
SQUARE 192x192 system obtained from the rank-13 CP equations of T_F by
freezing a declared 55-dimensional coordinate slice. The margin algebra is
the SAME one gate A certified for the rank-25 certificate (exact fmpq
quantities + outward-rounded arb final comparisons; explicit midpoint/radius
control throughout — the python-flint arb(lo, hi) MIDPOINT/RADIUS landmine is
respected: every arb is built from an exact rational midpoint and an exact
rational radius, radii are NEVER formed by float endpoint subtraction).

### 2.1 The square system, counted explicitly (Main condition 1)

Equations: g(x) = CP_13(a,b,c) - T_F, one per tensor entry:
  g[p,b,c] = sum_{r=1..13} a[p,r] b[b,r] c[c,r] - T_F[p,b,c].
  Count: 3 * 8 * 8 = 192 equations.
Unknowns: x = (a: 3x13, b: 8x13, c: 8x13), (3 + 8 + 8) * 13 = 247 unknowns.
Gauge: CP parametrizations carry a 13-dimensional continuous column-scaling
gauge (a,b,c) ~ (a D^-1, D b, D c), D diagonal 13x13, and a discrete S13
permutation gauge; additionally H(+)H doubling structure is not a gauge.
Gauge is NOT subtracted algebraically; instead certification works on a
DECLARED AFFINE SLICE that freezes exactly 55 = 247 - 192 = (13 scaling +
42 excess) coordinates at fixed exact dyadic values. On the slice the system
is square: 192 equations in the remaining 192 unknowns (the "S-block"); the
frozen coordinates ("F-block") enter as exact constants. The arithmetic:
  247 (unknowns) - 55 (frozen) = 192 = 3 * 8 * 8 (equations). Balanced.
Any root of the square slice system, extended by the frozen F-coordinates,
is a full 247-coordinate exact CP decomposition of T_F. Hence a PASS implies
rank(T_F) <= 13; with the certified floor 13, rank(T_F) = 13 EXACTLY.

Slice rule (fixed NOW, mechanical, no post-hoc choice):
  Slice-1 (primary): column-pivoted QR (scipy.linalg.qr, mode='economic',
  pivoting=True) on the float64 Jacobian J at the certification seed; take
  the first 192 pivot columns in pivot order as S; F = the remaining 55
  columns in ascending natural order.
  Precondition for Slice-1 validity: float sigma_min(J restricted toCandidate
  S) >= 1e-9 diagnostic and rank(J) = 192 at tolerance 1e-9 at the seed.
  If the precondition fails: Slice-2 (fallback) S = columns 0..191 in natural
  order, F = columns 192..246. Same script applies the rule mechanically;
  both rules were fixed before any run; which rule fired is reported.
The float QR is used ONLY to choose a discrete column order. All
certification arithmetic on the chosen square block is exact (fmpq over
dyadic rationals). Pivots/singular values are recorded as
COMPUTATIONAL-EVIDENCE diagnostics; the load-bearing facts (Y exactly
invertible as an fmpq matrix; margin inequalities strict) are
MACHINE-VERIFIED exact.

Transversality diagnostic (exact): J_F = d g / d x_F at the seed, exactly on
the chosen slice; report exact fmpq rank (fmpq_mat.rref) and sigma_min at
1/10^k thresholds. If rank(J_F) < 55 the slice is structurally handicapped
(an off-slice root may exist arbitrarily close in F-coordinates); this is
REPORTED, and any subsequent failure is labelled slice-limited — it is not,
by itself, any statement about rank 13 feasibility off-slice.

### 2.2 Certification margins (gate-A algebra, restricted to the slice)

Exact data at the seed x0 (every float64 value parsed to its EXACT dyadic
rational via Fraction(str-repr semantics of float.hex / Fraction(float)):
  g0 = g(x0) exact (fmpq); J_S exact 192x192 (fmpq_mat);
  Y = J_S^{-1} exact (single fmpq_mat inverse; ~3 s measured for 192 —
  feasible);
  bounce c0 = Y g0 exact.
For a box radius rho (each S-coordinate independently |x_s - x0_s| <= rho):
intervals for x_S = x0_S + [-rho, rho]^192 (F fixed at x0_F exactly), and
  K_c = x0_c - c0_c + [-R2_c(rho), +R2_c(rho)]
with R2_c(rho) <= sum_{e} |Y_{c,e}| * Q_e (exact fmpq bound, coordinatewise
Cauchy contraction, exact abs entries) and the EXACT per-equation quadratic
form factor Q_e from trilinearity: for equation e = (p, b, c):
  q_e(d) = sum_{r=1..13} a0[p,r] * db[b,r] * dc[c,r],
  |q_e(d)| <= Q_e := sum_r |a0[p,r]| * rho^2    (since |db|, |dc| <= rho).
The bound is EXACT (no slack beyond the coordinatewise absolute-value
contraction; the linear-in-d pieces are absorbed exactly by the constant
Jacobian term, and there is NO first-derivative remainder because the CP
Jacobian is constant — the system is multilinear, all higher derivatives
vanish).

Containment test, per coordinate c:
  |c0_c| + R2_c(rho) < rho        (strict, exact fmpq comparison).
PASS at rung rho := all 192 comparisons strict AND the outward-arb
re-evaluation of the same comparison on both sides confirms strictness
(K upper bound < 1 form, cf. gate A [A3]). Then Krawczyk's fixed-point
theorem (Moore / Krawczyk; the exact-arithmetic form of the gate-A margin
algebra, cf. arXiv:2608.16649 Prop 6 semantics replayed in gate A) gives an
exact root of g on the slice inside Box(x0, rho): an exact rank-13 CP
decomposition — certificate: MACHINE-VERIFIED.

Rung ladder (FIXED, every rung reported, none skipped, none added, rungs
evaluated in the listed order with early-stop on first PASS):
  rho in {1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 3e-9, 1e-9, 3e-10,
          1e-10, 1e-11, 1e-12}.
The ladder extends Main's cited 1e-2..1e-8 window DOWNWARD ONLY (monotone
safe: smaller boxes are never wider) because the bounce-vs-quadratic-margin
crossover is not predictable a priori (|Y| ~ 1/sigma_min ~ 1e4..1e6 and the
crossover rung can sit below 1e-8). All Main's rungs 1e-2..1e-8 are present
verbatim.

### 2.3 Seed conditioning pipeline (FIXED; deterministic; no new seeds)

The frozen seeds carry residuals 1.3e-3 / 1.0e-4, far above the certificable
bounce magnitude (the bounce must beat rho *and* the quadratic margin at the
same rung). The declared conditioning pipeline (part of THIS
pre-registration, not a new search):
  P0: load S_CP / S_EXT. Parse to exact dyadic x0. Verify (float diagnostic)
  the stored rel values reproduce to 1e-12 (anchor on parse correctness).
  P1..P5: five rounds of trust-region-reflective least squares
  (scipy least_squares, method='trf', tr_solver='exact', x_scale='jac',
  analytic Jacobian, max_nfev=5000 per round, ftol=xtol=gtol=3e-15), each
  round followed by exact-rational-independent column balancing
  (geometric-mean equalization per rank term — the balance() convention of
  route_af_search.py), which is a gauge move: it preserves the CP image
  exactly in real arithmetic (float drift is inherent; the FINAL candidate is
  re-parsed to dyadic values before any exact step, and the parse note is
  part of the record).
  Candidate selection rule (fixed NOW): the round output with the smallest
  relative Frobenius residual; ties break toward the earlier round. S_CP and
  S_EXT pipelines are SEPARATE and BOTH are fully reported; the certification
  attempt uses the smaller-rel candidate of the two (tie: the CP-direct
  pipeline).
  Factor bound: every |factor entry| <= 1e4 required (route AF's factor-norm
  convention). If the candidate violates it the certification attempt STILL
  RUNS (no post-hoc rejection) and the violation is reported alongside the
  outcome.
  "Converged" trigger: rel <= 1e-11 (the same trigger the Route AF campaign
  froze). If no candidate reaches 1e-11, the certification attempt STILL
  RUNS on the selected candidate at every rung; the expected outcome is
  failure at every rung, reported as FAILURE TO CERTIFY with the full
  residual data — it makes NO rank claim in either direction.

NO other seeds, NO random restarts, NO further polish rounds, NO box
widening, NO re-slicing. If the ladder fails everywhere: FAILURE TO CERTIFY
with (i) the certification box (rho set that was tried), (ii) every rung's
bounce / margin data, (iii) the exact per-coordinate classification (which of
the 192 coordinates straddle: p-slice / b-row / c-col indices listed), (iv)
the exact J_F transversality diagnostic, (v) the full rel history per
pipeline. A near-miss is not convertible into any claim (Main condition 4).

## 3. Mandatory controls (all through the SAME instrument code)

C-tau-7 (Main condition 3, certify direction): the quaternion tau
(3x4x4, true rank exactly 7) through the identical square-slice path at
r = 7: 64 equations (4*4*3), 77 unknowns, frozen 13. Seed: the FROZEN
tau_r7 certificate factors parsed exactly as dyadics
(scratch/upstream_ref/certs/tower_cert/tau_r7/ A/B/C). Expected: containment
fires at some rung (existence of an exact rank-7 decomposition near the
stored certificate — independently re-certifying gate C's existence claim in
this instrument). Record the rung and the margins.

C-tau-6 (Main condition 3, fail direction + Main condition 2 near-miss):
tau at r = 6: 64 equations, 66 unknowns, frozen 2. Declared seed: the r=6
truncation of the tau_r7 witness (first 6 columns, exactly parsed), then the
SAME declared 5-round TRF polish (max_nfev 5000/round). The true rank of tau
is EXACTLY 7 [MACHINE-VERIFIED in-repo, both directions], so the r=6 system
has NO root anywhere; the instrument MUST NOT certify containment (PASS
absent at every rung). The strongest rejection recorded is the disjointness
branch: an outward interval evaluation of the Krawczyk map over the whole box
proving K(Box) and Box are DISJOINT at some rung — a certified no-root
statement over a NAMED box. If interval widths swamp the displacement and no
disjointness, the honest report is "containment correctly absent at every
rung; local exclusion not achieved at the rungs tried" plus the citation of
the global rank-7 fact; either way the r=6 plant verdict is REJECT-EXISTENCE
(no containment ever fires on a known-impossible system).

P-pos (positive control, route-F shape): synthetic rank-13 tensor
U = CP(a*, b*, c*) with factors drawn from numpy default_rng(seed=20260901)
normal(scale=0.7); plant seed = exact dyadic factors + 1e-6 *
default_rng(seed=20260902) normal perturbation; instrument MUST certify at
some rung (the exact factors are a root at known distance <= sqrt(247)*1e-6
in the infinity norm on each dyadic coordinate — the box radius 1e-2 covers
it with margin). Validates the positive direction end to end including the
exact inverse and margin comparison.

P-wrong (assignment item 5c, deliberately wrong candidate): the SAME main
square system on the CORRUPTED table T_cor = T_F with entry [0,0,0] set to 2
(integer, exactly representable), from the same main candidate seed at the
same slice. A truthful instrument MUST NOT certify containment (the root
certified for T_F does not solve the corrupted system — the residual at the
T_F root is exactly 1 in the [0,0,0] entry). Expected: containment absent at
every rung; disjointness expected to fire at neighboring rungs (the corrupt
entry shifts the residual by exactly 1.0 in one coordinate, c0 displaces
correspondingly). Report exactly which branch fires at which rungs. Behavior
contract: if containment fires on P-wrong, the instrument is BROKEN and every
main-campaign result is QUARANTINED — the campaign verdict becomes "INVALID
INSTRUMENT" and nothing else is reported.

P-wrong is a NEAR-MISS by construction (Main condition 2): the corrupted
table differs from the true table in exactly one entry (1 vs 2, relative
entry magnitude 1/1), one coordinate among 192 — deliberately a SMALL
residual difference, not an obviously-bad candidate; the T_F-certified root
is a least-squares-1e0-witness for the corrupted system in the corrupted
coordinate only.

Anchor replays (assignment item 5b, rule 14): the tau control through the
certification path MUST certify at r = 7 AND MUST fail (no containment) at
r = 6 (C-tau-7/C-tau-6 above: this is the known-case discrimination test). An
instrument that cannot distinguish the known case is not an instrument, and
nothing else in this campaign may be reported as evidence.

## 4. Upward side (assignment item 4), attempted at low cost, honestly bounded

No exact argument that 13 is IMPOSSIBLE for this slice class is known to this
campaign (none found in the Route AF / Route A / S3 audits, re-read
first-hand this session). What is attempted and reported:
  U1 (audit, cited not re-derived): routeF_theorem_audit.txt negatives — no
  covering additivity/multiplicativity theorem; Christandl-Jensen-Zuiddam
  Prop 22 inapplicable at a=3 over R; direct-sum additivity inapplicable;
  Strassen route closed at 12 < 13 by primary read.
  U2 (pre-committed necessary-condition probe, exact): the block 2x2
  structure of the (1,i,j) family (L_1 = I_8, L_i, L_j block-diagonal with
  quaternion blocks; s3_routeF2_blocks facts) forces on any rank-13
  decomposition the I-slice equations
  sum_r a[0,r] b[b,r] c[c,r] = delta_{bc}; the 13-term outer-product
  representation of I_8 has graded constraints (the 13 pairs
  (b_r, c_r) must cover the 8 diagonal positions and the pencil-12 witness
  structure C4 pinned). The probe computes the NECESSARY linear-algebraic
  conditions exactly (fmpq) and records the numbers; if the necessary
  conditions are satisfiable, they give NO impossibility — reported as such.
  U3: Koiran Thm 4 exact characterization (rank <= 13 iff a commuting
  diagonalizable extension exists) is noted as the exact upward natural
  system; NO Groebner/feasibility attack is in budget (pre-declared:
  195+ unknowns, out of scope for the low-cost clause).
ABSENCE OF AN IMPOSSIBILITY ARGUMENT IS NOT EVIDENCE FOR 13, and a
FAILURE-TO-CERTIFY on the downward instrument is NOT EVIDENCE FOR 14. Both
sentences are part of every verdict wording (Main condition 4).

## 5. Adjudication (fixed now; Main's vocabulary verbatim)

- A contraction PASS settles the question: any ladder rung PASSES containment
  on the main system with margins strict end to end (exacted as in 2.2 and
  outward-arb confirmed): verdict "rank(T_F) = 13 EXACTLY"
  [MACHINE-VERIFIED], because 13 is already the certified floor. A dedicated
  REPLAY script in this campaign re-derives Y, the bounce, the margins, and
  the containment inequalities from frozen bytes and exits 0.
- Anything else: "FAILURE TO CERTIFY" with the achieved contraction/margin
  data per rung, the rungs tried on the fixed 13-rung ladder, and which
  coordinate blocks straddle. No widening the box; no adding seeds; no
  converting a small residual into a claim.
- Any P-wrong / C-tau-6 behavioral break: instrument quarantined; campaign
  verdict "INVALID INSTRUMENT"; nothing else reported.
- The published window 18 <= R_R(T_O) <= 25 is untouched by every outcome in
  this campaign; 13 does not bear on it (the chain feeding 18 uses the
  3-family floor 13 only as already certified, and this campaign adds no
  lower-bound machinery).

## 6. Budget (fixed)

- Polish: 2 main pipelines x 5 rounds x 5000 nfev + 1 tau-6 pipeline (5
  rounds at 64x66) + utilities; TRF exact-solver; <= ~5 min CPU expected;
  hard cap 30 min CPU (time.process_time() accounting printed by every
  script; rule 17e: no ps-derived numbers anywhere).
- Exact steps: one 192x192 fmpq inversion per certified system (main + 4
  controls), margin evaluation at <= 13 rungs each; expected <= 2 min total.
- nice -n 10; single-threaded BLAS (OMP/OpenBLAS/MKL/vecLib/NumExpr = 1).
  Main confirmed 2026-09-01 that no slot serialization is needed (rule 17e
  context: the 2.4% readings were a ps accounting artifact).
- No daemonized processes; every script exits; no unbounded loops (polish
  rounds are capped by max_nfev; ladder length is fixed).

## 7. Outcome vocabulary (fixed; repeats Main condition 4 verbatim in force)

PASS = "rank 13 established, MACHINE-VERIFIED" (only via section 5 bullet 1).
Everything else = "FAILURE TO CERTIFY" + contraction/margin data + the rungs
tried + straddling blocks + the two absence sentences verbatim:
  "Absence of a 13-witness is NOT evidence for 14; absence of an
  impossibility argument is NOT evidence for 13."

## 8. Rule-7 prospective scope

Covered: the exact (L_1, L_i, L_j) tensor in the fixed Cayley-Dickson basis;
target rank 13; the two frozen seeds; the declared 5-round polish; one
declared slice rule per system; the fixed 13-rung ladder; the four declared
controls (C-tau-7, C-tau-6, P-pos, P-wrong); the audit-grade upward attempt
U1-U3 with its pre-committed no-impossibility conclusion unless an exact
proof emerges in-run.
NOT covered: any other seed (including any restart of any failure branch),
any other slice, any widened box, rank-14 certification in either direction,
border rank, complex decompositions, other octonion triples, any statement
making the (1,i,j) 13-vs-14 decision except a section-5 PASS or an explicit
FAILURE-TO-CERTIFY report, and anything touching R_R(T_O) itself
(18 <= R_R(T_O) <= 25 stays untouched; no published work is refuted).
