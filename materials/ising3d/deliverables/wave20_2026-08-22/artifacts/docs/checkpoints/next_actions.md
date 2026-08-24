# Next actions

**Snapshot:** 2026-08-22. Waves 18, 19, and 20 LANDED and were independently re-verified.
This file is the single source of truth for the ranked open queue, ledger/experiment allocation,
and environment state. Strategic reasoning lives in `notes/frontier_analysis.md`; do not duplicate
one into the other.

**Repository:** `/Users/jinleic/jinleic-workspace/math/ising3d` · **Interpreter:** `.venv/bin/python`
· `K_c = 0.221654626` stays benchmark-only · the ledger is append-only · `.venv` is shared and
immutable for agents.

## 1. Ledger and numbering state (read before appending anything)

- **LEDGER next free ID: `H584`.** Allocation to date: `H001`-`H505` historical through wave 19;
  `H506`-`H577` are the twelve wave-20 front blocks; `H578`-`H583` are the lead finite-memory-SAW
  block. All 78 rows are present, 11 fields wide, uniquely numbered, and verifier-reviewed.
- **Next free experiment number: `e228`.** `e191`-`e226` are the twelve wave-20 fronts;
  `e227` is the lead finite-memory SAW automaton.
- **Width anomalies (historical, do not edit):** `H388` has 12 fields; `H389`/`H390` have 6;
  `H393`-`H397` and `H399`-`H401` have 10. Recorded by append-only row `H402`. `H383`/`H385` are
  DUPLICATEs of `H380`/`H382`; `H384`/`H390` are corrections, per the canonical mapping in `H386`.
  `H360` predates the schema correction and is superseded by `H361`.
- **Field-convention drift, wave 17 only:** rows `H417`-`H427` are 11 fields wide but use
  `id, date, statement, topic, owner, status, artifact, artifact, detail, status-word, next-action`
  rather than the legacy header semantics. Wave-18 rows return to the legacy header.

## 2. Wave 20 landed — do not re-open these

- **Critical floor moved (headline).** The deterministic memory-12 simple-cubic finite-walk
  automaton has 41,424 states and an exact positive integer certificate
  `1000000 A w < 4747526 w` (minimum residual `18,990,104`), hence
  `mu<4.747526` and
  `K_c>0.2138368062108697304302094293879379891333`. The clean-room verifier independently
  rebuilds every automaton `k=4..12` and all 169,975 target transitions
  (`proofs/kc_lower_finite_memory.md`, `H578`-`H583`).
- **Connected-bipartite local-term classification completed.** Every `Delta>=3` graph is the
  nonradical `Q_c=1` branch, with dimension `2^(2n-2)-1` (odd `n`) or
  `2^(2n-2)-(-1)^r 2^(n-1)` (even `n`, colour size `r`). The full Seven v1 source is archived
  and hashed; every theorem hypothesis is proved locally. This removes the Hamiltonian-path
  restriction and, with wave 19, classifies every connected simple graph
  (`proofs/nonhamiltonian_bipartite.md`, `H572`-`H577`).
- **Every-coupling isotropic theorem on the non-bipartite locus.** The determinant-center-free
  invariant `Phi_N` reduces to a positive odd-Eulerian trace skew, excluding every full
  subset-product spectrum for every `0<t<1`; it vanishes identically on bipartite graphs, so
  rectangular/simple-cubic layers remain the exact blind locus
  (`proofs/isotropic_invariant.md`, `H524`-`H529`).
- **New exact structures.** `dim h_(K4,4)=137`, `dim h_(K4,5)=471`; the conjectured `W_L`
  defect is the unique possible `SL_L` top-wedge hyperplane, while every natural monomial map
  is too small; `c8` uniquely forces `W8^2/8!`, `c10` forces `W10^2/10!`
  (`H506`-`H517`, `H566`-`H571`).
- **Routes closed exactly.** No rank-one/`|+>` ladder-to-rectangle induction; Galois `S6`
  occurs in the free-fermion control; a finite 3D CTM has contraction-dependent spectra;
  local-term dimension cannot imply an all-auxiliary Lax no-go; the named four-point
  power-simplex relaxation projects to the old two-point class; susceptibility
  supermultiplicativity fails at `a_2=30<36`. The fixed Pluecker signature fails on both a
  cube and a fresh slab (`H518`-`H565`).
- **Verification.** All 13 new standalone verifiers passed in the lead checkout. The initial
  susceptibility verifier correctly failed on a false “strict” ratio statement; producer,
  proof, artifact, and verifier were corrected for the exact `5,5` plateau, then passed.

## 3. Wave 19 landed — do not re-open these

- **Non-bipartite trichotomy (now part of the wave-20 complete theorem).** Odd cycles give
  `2n(2n-1)`; every other connected non-bipartite graph gives `2^(2n-1)-2`. Wave 19 also found
  that Jordan--Wigner grade-compatible orders exist iff the bipartition classes differ by at
  most one, and that the Hamiltonian-grade proof could not cover `K1,3`. Wave 20 supersedes the
  resulting “hypothesis necessary” wording: the algebra itself is order-independent and the
  transvection proof classifies every branching bipartite graph, including `K1,3=72`.
  Every finite `a x b x c` box remains covered, with
  `2x6=3x4=4192256`, `3x5=268435455`, `4x4=1073709056`
  (`proofs/trichotomy_extend.md`, `proofs/nonhamiltonian_bipartite.md`).
- **The two-generator algebra does NOT inherit the trichotomy.** Branching is not sufficient:
  `K_{2,3}` gives `44 < 45` and `K_{3,3}` gives `63 < 66`, the latter explained by its
  automorphism-fixed container of dimension `64`. Fixed-path census: no branching case at `n=4`;
  `2/3` clear at `n=5` (`K_{2,3}` the exception); `13/14` at `n=6` (`K_{3,3}` the exception); all
  `63` clear at `n=7`, which is finite, so only `n0 >= 7` is claimed. Proved for the physical family:
  all `2xL` with `L>=3` (by citation), plus `3x2`, `3x3`, `3x4` (`proofs/twogen_allsize.md`,
  `H488`-`H493`).
- **All-`L` ladder no-go (from the retargeted `LadderW8` front).** For every `L>=3` the ladder
  transfer generators admit no simultaneous Majorana-quadratic representation:
  `dim <A_L,B_L> > 2L(4L-1)`. Contiguous seam table: exact dimensions at `L=3,4`; a word certificate
  `191` at `L=5`; `K_L-W_L` at `L=6..9`; sector-capped state certificates at `L=10..31`, each hitting
  rank `8L^2-2L+1` exactly; and the unconditional theorem `dim >= Q(L) > 8L^2-2L+1` for `L>=32`,
  where `Q(L)` counts nonempty partitions with all parts `>= 2` and total `<= L`. `L=2` is honestly
  inconclusive (`11 <= 28`, the integrable 4-cycle ring), in contrast with the local-term algebra
  where `L=2` IS excluded (`56 > 28`). A `[CONJECTURE]` closed form reproduces all seven certified
  `W_L` with zero fitted parameters: `dim U_L = (C(2L,L)+2^L)/2 - [4|L]`, `W_L = K_L - dim U_L`,
  predicting `W_10 = 38950` (`proofs/ladder_w8.md`, `H458`-`H463`).
- **Two proof routes exactly closed.** `ad_A` Vandermonde grade separation cannot isolate chords: every
  Ising edge monomial, path bond or chord and regardless of grade, has `ad_{iA}` spectrum
  `{0, ±4i}`. Group rigidity cannot deliver the isotropic spectral theorem: `exp(aH)exp(bX)` is
  pointwise conjugate into `Spin(2,C)` on an open set while generating `sl(2,C)`, so a varying
  conjugator absorbs the transverse derivative (`proofs/layer_group_spectral.md`, `H500`-`H505`).

## 4. Wave 18 landed — do not re-open these

- **All-`L` classification (headline).** For EVERY `L>=2` the open `2xL` local-term algebra is exactly
  `so_m(Q)+so_m(Q)` for even `L` and `sp_m(Q)+sp_m(Q)` for odd `L`, `m=2^(n-1)`, `n=2L`, with
  `dim = 2^(n-1)(2^(n-1)-(-1)^L)`. The saturation step left open since wave 16 is closed by an
  explicit production induction with invariant `I(L)`, a square-triple second move, a no-rank-loss
  transport lemma, and leaf-induction token connectivity. Full closure enumerated at `L=2..6`
  (`56, 1056, 16256, 262656, 4192256`, zero set-law violations, 19 process seconds and 437 MB at
  `L=6`); `L=7` is certified by construction plus an independent Clifford census as `67117056` and is
  explicitly NOT a full enumeration (`proofs/alll_saturation.md`, `H428`-`H433`).
- **Exact trichotomy (headline).** For any bipartite `Gamma` on `n` vertices with a Hamiltonian path,
  `N=2n`: the path gives grades `{2}` and `dim = n(2n-1)`; an even cycle gives `{2,N-2}` and
  `dim = 2n(2n-1)`; and `Delta(Gamma)>=3` gives every grade `= 2 mod 4` except `k=N`, with
  `dim = 2^(2n-2) - (-1)^(n/2) 2^(n-1)` for even `n` and `2^(2n-2) - 1` for odd `n`. In the branching
  branch the answer depends only on `n`. `C_6` (dim `132`) is the exact counterexample that fixes the
  hypotheses; `C_4 = 2x2` is in the cycle branch but accidentally has the full class. Non-bipartite
  `C_5`+chord gives `510`, the full even-grade algebra (`proofs/clifford_grade_classification.md`,
  `H480`-`H487`).
- **All-`L` quadratic no-go.** `dim g_L >= n(2n-1)+1` for every `L>=2`, so the local terms are never
  simultaneously Majorana-quadratic on the physical space; exact conservative mode floors
  `6, 18, 65, 216, 991, 3168, 15354` at `L=2..8`. Covers `L=2`, where the claw argument does not apply
  (`proofs/alll_quadratic_nogo.md`, `H434`-`H439`).
- **Isotropic gap: split verdict.** The core half SUCCEEDS — paw `119>65`, five-site tree `417>211`,
  `2x2`+pendant `475>211`, with monic upgrade and exceptional-degree bounds `913920`, `17677440`,
  `19641600`; controls `K_{1,3}` `65=65` and `K_{1,4}` `147<211` explain why the bare claw never
  sufficed. The decoupling half provably FAILS on the isotropic curve: crossing-factor minor
  `q(t)^(2c)-1`, and `r(K_{1,3} + m K_1) <= 136(2m+1) < 3^(m+4)-2^(m+4)`. Any isotropic all-size
  spectral proof must be non-localizing (`proofs/isotropic_allsize.md`, `H440`-`H445`).
- **Theorem U2 restored.** Exact char-0 `r(T)=417` and `r_all(T)=445` at four named physical points
  (derivative-gcd degrees `79` and `83`); the exact core upper bounds squeeze the modular lower bounds
  to certify `3893`, `4005`, `35597`, so hypothesis CZ is discharged and Theorem F2 is unconditional
  for `m=2,4`. `INJ(n-5)` for other chain lengths stays open
  (`proofs/rank_char0.md`, `H446`-`H451`).
- **Kac-Ward: no orbit reduction, plus a complete unit census.** The sound section group on the
  diagonal-normalised thin chart is trivial, so there are exactly `16384` singleton orbits and no
  quotient reduction exists. The census was then completed anyway by exact character evaluation over
  all `4^13 = 67108864` unit points: `2960` construction zeros in `2437` sections, each with an exact
  disposition (`2232` order-8 mod 5, `664` order-12 mod 5, `45` nonsingular mod 625, `19` singular
  lift obstructions). Classification `0 / 2437 / 13947` for
  `EMPTY_OVER_Q / HAS_LOCAL_POINT_BUT_FAILS_HOLDOUT / UNDECIDED`. Branch `11111` remains OPEN over
  `Q`: higher 5-adic digits, non-unit valuations, zero reductions, chart-degenerate loci and other
  primes are untouched (`proofs/kw_orbit_quotient.md`, `H452`-`H457`).
- **Upper endpoint: two classes exactly refuted.** At `v=6/25` the aggregate one-moment
  paired-momentum SDP has exact optimum `sup s = 1` and uniform order floor `eta* = 0`, so the
  `eta = 1/1024` challenge is infeasible with separating gap `1/1024`; the single-edge conditional
  domination random-current LP has exact optimum `p* = v^2 = 36/625` against the self-proved planar
  threshold `5/6`, short by `2909/3750`, and `p* < 1/16` uniformly. No improvement to the incumbent
  upper endpoint. Mode-resolved four-point/DLR, multi-edge or block current, and sourced-current
  classes are untouched (`proofs/upper_endpoint4.md`, `H464`-`H469`).

## 5. Ranked queue after wave 20

1. **The bipartite isotropic spectral theorem.** Wave 20 exactly solves the non-bipartite locus,
   but the trace-skew invariant vanishes on every rectangular grid. Need a non-reciprocity
   conjugacy invariant nonzero modulo `(y-1)x^2-(y+1)=0`; do not retry localization or group rigidity.
2. **Two-generator grids beyond ladders.** Prove the physical inequality for all `3xL` and then
   rectangles with `min(a,b)>=4`. Rank-one frozen rows and `|+>` compression are closed by an
   all-size depth-three anomaly; any induction needs a higher-rank invariant code or direct words.
3. **The exact `W_L` rank law.** The exceptional `4|L` middle sector now has its canonical
   invariant-theory model, but both rank inclusions remain open. Natural signed monomial maps are
   impossible. Target a direct hard-core/Krylov proof of sector two, or one exact `L=10`
   counterexample; do not fit another sequence.
4. **Memory-14 finite-walk automaton.** This is the cheapest likely endpoint gain. Preflight:
   roughly 620k states from the observed state-count ratio, not an observed wall. First replace
   tuple states and the unbounded closability cache by compact deterministic step codes.
5. **The upper endpoint.** `MR4-power-simplex-L1` joins the closed classes. A viable certificate
   must use phase-sensitive Fourier identities, genuinely spin-realizable four-point/DLR
   constraints, or multi-edge/sourced random currents.
6. **Kac-Ward `11111` beyond the unit census.** Higher `5`-adic digits, non-unit valuations,
   zero reductions, chart-degenerate loci and other primes remain. Nullstellensatz and orbit
   quotient routes are proved unavailable.
7. **`s<=6` extensive charges — structural idea required.** `s<=5` is exactly quotient one.
   Both symmetry-compression routes are closed; do not launch a larger version of either.
8. **HT `v^30`, LT `x^60`, LT `x^62`.** Resource frontiers, not mathematical breakthroughs.
   `x^60` is profile-verified at 49.12 GB; `x^62` needs a 30-spin cross-section; `v^30` needs a
   generalized cut-capped frontier. Launch only when the host clears.
9. **`3x4` pair-product certificate.** Dimension 4096; cost every stage with `process_time`
   before launch. This is a stronger finite certificate, not the bipartite all-size theorem.
10. **Evaluate the forced `W8/W10` interlayer data.** Universal closure by lower cumulants is
    impossible. Seek square-lattice Pfaffian evaluation or rigorous bounds for the actual infinite sums.
11. **Actual higher-dimensional/spectral Lax equations.** The local-algebra/Schur shortcut is
    refuted. Progress requires a concrete nonidentical/spectral/IRF/dynamical equation.
12. **Lee-Yang edge scaling.** Needs a no-earlier-root/interlacing/argument-principle theorem or
    correction-amplitude control, not larger fitted cubes.

## 6. Environment state (measured 2026-08-22, 14:16 local)

- Host: 28 cores, 96 GB; load average measured **~65** with nine users. About 19 GiB was free, but
  CPU contention is severe and sustained.
- Keep new work single-process, nice/low priority, under 2 GiB RSS, and budget with
  `time.process_time()`, never wall clock. Heavy series producers remain ENVIRONMENT-BLOCKED.
- A verifier's self-declared wall timeout is a resource artifact; recheck that test in isolation
  before calling a mathematical regression.

## 7. Method notes worth reusing

- **When an exact upper bound exists, construct a matching lower bound rather than refining the upper
  bound.** A matching pair is a theorem; a tighter upper bound is only a better estimate. This is how
  the all-`L` classification and the trichotomy were found, both in one wave, after the containment
  theorem had pinned the upper bound since wave 17.
- **Test every new dimension formula against every already-certified value before claiming novelty.**
  That check turned `3x3 = 65535` from an anomaly into the odd-`n` confirmation of the trichotomy.
- **`ad_A` grade separation is closed for the two-generator problem.** Although `ad_A` preserves
  Clifford grade, every Ising edge occupies exactly the same eigenvalues `{0,+4i,-4i}`; a
  Vandermonde projector recovers only three global bond components and cannot isolate chords.
- **Monic-in-parameter upgrade.** A spectral invariant monic in `z` over `Z[parameter]` can only gain
  derivative-gcd degree under specialization, so ONE evaluated point bounds the generic value and one
  resultant confines the exceptions to a finite set.
- **Direction discipline.** Modular ranks and gcd degrees are one-sided. Substituting a bound into the
  wrong side cost this program Theorem U2, which took four waves to restore.
- **Razor margins are one-sided.** `deg gcd <= 295414 < 295415` proves "at least one above the
  ceiling", never "exactly one".
- **Content digests, not byte digests**, for cross-artifact references.

## 8. Standing hygiene

- Run every new standalone verifier from the repository root with the repository interpreter before
  editing any summary.
- Distinguish an OBSERVED wall (measured process time and RSS) from an unlaunched input-size
  preflight; label the latter `PREFLIGHT`.
- Verifiers must fail, not warn, when a row of a universal claim is uncovered.
- Scratch files `tmp_scout_simon.py`, `tmp_elim_simon.py`, `tmp_general_tm.py`, `tmp_trg_simon.py`,
  and `checkpoints/hypotheses_pre_overwrite_snapshot.csv` still await user deletion approval.
- Worktrees under `/private/tmp/ising3d-w1*-*` and their `omp/w1*` branches can be removed with
  `git worktree remove` once no comparison is needed; ask before deleting.
