# Next actions

**Snapshot:** 2026-08-18, after wave 15 (ledger H380-H403 appended; 375 data rows, 368 unique IDs with 7 correction/mapping rows,
highest ID `H403` — 28 H-numbers from earlier agent ID blocks were never allocated, so row count and
highest ID do not coincide; all rows uniform 11 columns; `H383`/`H385` are marked DUPLICATEs of
`H380`/`H382` and `H384`/`H390` are corrections, per the canonical mapping in row `H386`; row
`H388`/`H389`/`H390` column misalignments — recorded by append-only row `H402`, rows not edited).
**LEDGER next free ID: `H417`.** (Wave-16 rows H404-H416 landed; four of them are lead-takeover rows after harness terminations of the agent executions.) Acceptance
suites: wave 10 `86/86`, wave 11 `90/90`, wave 12 `94/95` (the one failure is
`tests/test_ladder_all_l.py` hitting its own `NON-DECISIVE` 3600 s wall under heavy machine load; it
PASSES in 1565 s when run in isolation), wave 13 `102/103` (same load artefact; every wall-clock
budget gate has since been converted to `time.process_time()`), wave 14: fullsuite14 110 total,
109 passed, 1 failed (`tests/test_extensive_3d_r2.py`, dropped artifact-load line, lead
regression); standalone rerun after one-line restore: PASS, 555 s; full-suite re-run deferred to
wave-15 close-out. Repository:
`/Users/jinleic/jinleic-workspace/math/ising3d`. `K_c = 0.221654626` stays benchmark-only; the
ledger is append-only; `.venv` is shared and immutable for agents.

Waves 11-14 used one isolated `git worktree` per front, an independent read-only review, **and** a
lead re-verification inside the main checkout before integration. Keep all three steps: six of
eight wave-13 deliverables were corrected because of review, always for the same reason — a verifier
trusted a stored field, or a statement quantified beyond its evidence.

## Wave-17 candidate fronts (wave-16 landed: H404-H416)

The wave-14 fronts that finished are retired here (2x5 pair-product certificate, the
`SectorSaturationCyclicity` decision, `s<=5` in the `Z^3` radius-2 box, the Peierls-at-incumbent
limitation, the anchor-section KW component cleaning, and the `n = 36` sieve union); what replaces
them, in priority order:

1. **KW boundary lifts (Kac--Ward branch `11111`).** The anchor thin-torus section
   `(1,1,1,1,1,4,4)` is cleaned — exactly three nonsingular unit-chart `F_5` solutions
   (`det J` in `{2,4,1}`; the full `F_5^6` grid has 1553 solutions, the fourth nonsingular point
   chart-degenerate), each a genuine `Q_5` point of the 42-equation construction variety, and all
   three fail the independent `3x3x3` order-8/10/12 holdout mod 625; the orbit-union lemma
   collapses the diagonal torus action to one representative per orbit (`proofs/kw_components.md`).
   Next: lift the component census across the remaining `4^7-1` thin-torus sections (the full
   16384-section census) and either classify every `Q_5`/`Q(alpha)` component or prove each fails
   an independent finite holdout — the last step to the `Q`/`C` decision of the sole surviving
   branch.
2. **`L=8` ladder.** `SectorSaturationCyclicity` is refuted at `L=4..7` (`W_L = 0,2,10,66,364`;
   cyclic dims `14/14, 42/44, 142/152, 494/560, 1780/2144`), so the all-`L` route must now be
   argued with `W_L`, not cyclicity: either prove `W_L/K_L` stays below `1 - 2^L/K_L` for every
   `L` (rescuing `dim g_L >= 2^L` unconditionally) or find the first crossing. The next decisive
   datum is the `L=8` rung (`K_8 = 8384`): certify `W_8` and the `L=8` cyclic dimension on the
   three-probe architecture (tensor, duality pairing, state space;
   `proofs/sector_saturation_tensor.md`, `proofs/sector_saturation_pairing.md`).
3. **Char-0 quotient of the `2x3` algebra.** The individual-term DLA on `2x3` has dimension 1056
   with mod-`p` quotient witnesses `sp(14)` (dim 105), `sp(6)`, `gl(6)` certified at two primes,
   but no characteristic-zero simple quotient is known. Decide whether any modular witness lifts
   to an exact char-0 quotient — the same-basis exact-container architecture that settled the
   `2x3`/`2x4`/`3x3` two-sum algebras is the template.
4. **Octahedral `s<=6` (extensive charges).** `s<=5` in the `Z^3` radius-2 box is decided
   (quotient exactly `1`; 21,121,156 columns, two-prime CSR echelon meeting the analytic `I,h`
   bound — `proofs/extensive_3d_s5.md`). Next: the `s<=6` class and the full box, where the wall is
   combinatorial size; the order-48 point-group reduction (noted plausible in
   `proofs/extensive_3d_r2.md`, never done) is the natural symmetry cut.
5. **Paired-momentum gap (upper endpoint below `I_3/2`).** Both incumbent-route dead ends are now
   certified: the two-point infrared/GKS class is exactly optimal at `I_3/2` with an `O(1/L)`
   torus LP floor (`proofs/upper_infrared.md`, `proofs/mag_floor.md`), and the plain Peierls
   head+tail certificate provably cannot close (`K* = 0.2558326` sits `0.003102` above the
   incumbent, and its exponential tail converges only `1.42` above it — `proofs/upper_beyond.md`).
   What remains audited: paired-momentum/four-point Fourier constraints, random-current
   inequalities. Highest value, lowest probability.
6. **HT `v^24` feasibility.** The independent FLM box-class route reaches `v^22` and stops:
   `v^24` there requires the canonical `5x5x5` box (cross-section 25 > hard guard 22; projected
   three-array peak 242 GB — `notes/series_extension2.md`). Either generalise the cut-capped
   parity frontier that broke the same 242-GB wall on the primary route (`proofs/ht_v24.md`,
   352 MB) or produce a rigorous resource lower bound explaining the wall.

## Ranked after that

1. **Close the isotropic gap in Theorem F.** Theorem F is all-size but anisotropic; Theorem G covers
   only eight named isotropic graphs. The proved obstruction is precise: the isotropic claw is not
   certified above the ceiling, and equal-field localization is only linear in `m`. The repair that
   already worked once was replacing free sites by an open chain (distinct modes). Next: find a
   core that is certified above the ceiling at uniform fields *and* an isotropic decoupling whose
   remainder has distinct modes for every `n`, which would make the isotropic family theorem follow
   by the same tensor identity.
2. **Prove `r(T) = 417` and `r_all(T) = 445` exactly** (characteristic-zero upper bounds). That
   single computation would convert the withdrawn Theorem U2 into a theorem, pin the excess at
   `202*3^m + 4*2^m`, and discharge one of the two hypotheses of Theorem F2.
3. **`x^54` for the low-temperature series.** It is the exact order that decides the last 7 legacy
   truncation-unobservable rows. Blocker: the transfer engine's 22-spin cross-section guard. Either
   raise the guard with a memory-blocked transfer, or find a route that reaches `x^54` for the few
   box classes that matter.
4. **Beyond the `n = 36` sieve union.** The wave-14 two/three-block sieve closed the exact
   schedule family at `n = 35` (`proofs/saw_union4.md`); the remaining lever on the lower endpoint
   is an audited `c_37`/`c_38`, which would move the endpoint more than any further refinement at
   `n = 36`.
5. **`3x4` pair-product certificate.** The trace pipeline that cleared `2x5` (`dim 1024`,
   `C(1024,2) = 523776` pair slots) is the tool; `3x4` is `dim 4096` and needs its per-stage
   `process_time` budget arithmetic made explicit first (`proofs/pair_product_2x5.md`).
6. **Open `4x4` layer algebra.** Only `dim_Q g >= 1794`. The `3x3` route (sector projectors plus a
   same-basis exact upper container) is the template; the blocker is the `268591168`-dimensional
   container.
7. **Interlayer `c10`, and a verified closed form for `c8`.**

## Method notes worth reusing

- **Monic-in-parameter upgrade.** If a spectral invariant is a monic polynomial in `z` with
  coefficients in `Z[parameter]`, then specialization can only raise the derivative-gcd degree, so
  ONE evaluated point bounds the generic value, and one resultant confines the exceptions to a
  finite set. This turned three point certificates into coupling-generic theorems (waves 12-13).
  Check every new finite certificate for it before spending more compute.
- **Decoupling degenerations must keep modes distinct.** Equal parameters on decoupled pieces make
  their modes coincide and the pair-product count collapses to linear growth. Use pieces that are
  non-isomorphic or of opposite parity.
- **Direction discipline.** Modular ranks and gcd degrees are one-sided. A lower bound substituted
  into an upper bound proves nothing — that error cost Theorem U2. Certify a rational value only
  against an independently proved bound in the opposite direction.
- **Verifiers must cover every row of a universal claim** and must fail, not warn, when a case the
  note claims ends in a wall.

## Standing hygiene

- Run every new standalone verifier from the repository root with the repository interpreter before
  editing any summary.
- Distinguish an OBSERVED wall (measured time and RSS) from an unlaunched input-size preflight.
- The machine is shared with unrelated long-lived user processes; a suite failure that is a
  self-declared wall timeout must be rechecked in isolation before being called a regression.
- Scratch files `tmp_scout_simon.py`, `tmp_elim_simon.py`, `tmp_general_tm.py`, `tmp_trg_simon.py`,
  and `checkpoints/hypotheses_pre_overwrite_snapshot.csv` still await user deletion approval.
- Worktrees under `/private/tmp/ising3d-w1*-*` and their `omp/w1*` branches can be removed with
  `git worktree remove` once no comparison is needed; ask before deleting.

## Wave-17 candidates (from wave-16 outcomes, priority order)

1. **Alternation all-L proof**: closure/saturation uniform arguments from the certified set
   law (proofs/alternation_law.md gaps a/b/c); start with L=5 exact-Q closure to certify
   262656 over Q (currently modular probe only).
2. **x^56/x^58**: SAME box inventory as x^54, deeper truncation only — cheapest next
   coefficients; x^62 needs cross-section 30 (next guard lift).
3. **PM-S at Lambda=3,4** and new upper-endpoint input class (shell Lambda=2 structurally
   overshooting — proofs/shell_stress.md); massive-Watson bracket remains.
4. **2x6 spectral certificate**: calibrated (P,rho)-block route on record (blocks
   1056/992/1024/1024, kill-cap gcd<295415); rerun as front.
5. **5 surviving legacy LT rows**: test at x^56/x^58 when computed.
6. **s<=6 lower bound**: extract certified sum-of-block-ranks from the preserved
   /private/tmp/ising3d-w16-gelims6 pivot store, or close the question.
7. **3x3 layer DLA**: closure saturates 65535 = 2^16-1 (outside the 2xL alternation family);
   char-0 structure front.
8. **Kac-Ward**: other 4^7-1 thin-torus sections per H391 machinery.

## Operational state (2026-08-19, ~17:00 local)

- **Host congestion**: the workstation carries ~40 multi-day foreign jobs; new launches get
  ~2% CPU and large-RSS work (>= ~20 GB) thrashes. ALL heavy e15x producers and the full
  suite are ENVIRONMENT-BLOCKED; re-run when the box clears. Small-RSS work (e154, tests,
  proofs) proceeds slowly.
- x^60 producer (e155) profile-verified launchable (123 classes, 36 CRT boxes, peak
  49.12 GB) and got 36 min of CPU-sharestarved run; LAUNCH when congestion clears
  (est. 60-90 min at sensible share).
- x^56 independent audit (tests/test_lt_x56.py) needs a rerun under calmer load (~50-70 min).
- fullsuite acceptance run needed after wave-17 integration (previous run aborted at 6h51m).
- Gaussian2x6: e154 built from wave-16's terminated producer + 3 fixes (chi-product
  stabilizer test, false-symmetry assert, next_prime EVEN-return bug, cap design advisory).
  Controls/digest replay pass; heavy block certificates underway (small-footprint run OK).
