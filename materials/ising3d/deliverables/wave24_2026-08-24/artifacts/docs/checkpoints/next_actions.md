# Next actions

**Snapshot:** 2026-08-24. Waves 18 through 24 LANDED and were independently verified.
This file is the single source of truth for the ranked open queue, ledger/experiment allocation,
and environment state. Strategic reasoning lives in `notes/frontier_analysis.md`; do not duplicate
one into the other.

**Repository:** `/Users/jinleic/jinleic-workspace/math/ising3d` · **Interpreter:** `.venv/bin/python`
· `K_c = 0.221654626` stays benchmark-only · the ledger is append-only · `.venv` is shared and
immutable for agents.

## 1. Ledger and numbering state (read before appending anything)

- **LEDGER next free ID: `H658`.** Allocation to date: `H001`-`H505` historical through wave 19;
  `H506`-`H583` wave 20; `H584`-`H587` and `H590`-`H595` wave 21; `H596`-`H634`
  wave 22; `H635`-`H642` wave 23; `H643`-`H657` wave 24. `H588`-`H589` remain
  intentionally unallocated. `H653`-`H657` are append-only review corrections to
  `H644/H645/H649/H654/H655`; all wave-24 rows use the 11-field legacy schema.
- **Next unused experiment number: `e247`.** Wave 24 lands `e244` octahedral matchgate/sign
  structure, `e245` graph-uniform exceptional-scheme degree, and `e246` selected-row Callen
  Walsh compression. The pre-existing `e233_trace_resultant_positivity.py` remains an unlanded
  draft with no artifact/test/proof; do not claim it and do not reuse its number.
- **Width anomalies (historical, do not edit):** `H388` has 12 fields; `H389`/`H390` have 6;
  `H393`-`H397` and `H399`-`H401` have 10. Recorded by append-only row `H402`. `H383`/`H385` are
  DUPLICATEs of `H380`/`H382`; `H384`/`H390` are corrections, per the canonical mapping in `H386`.
  `H360` predates the schema correction and is superseded by `H361`.
- **Field-convention drift, wave 17 only:** rows `H417`-`H427` are 11 fields wide but use
  `id, date, statement, topic, owner, status, artifact, artifact, detail, status-word, next-action`
  rather than the legacy header semantics. Wave-18 rows return to the legacy header.

## 2. Waves 20-24 landed — do not re-open these

### Wave 24
Frozen pre-implementation contract and decision record:
`checkpoints/wave24_research_plan.md`.


- **Graph-uniform finite-scheme degree.** For every connected simple nonpath graph,
  `#E_G <= [2^n(2n+2m)]^(n+2)` and uniformly at fixed `n`,
  `#E_G <= [2^n(n^2+n)]^(n+2)`. Paths are excluded. The bound is algebraic control,
  not emptiness or a root list (`H643`-`H645`, corrections `H653`-`H654`, `e245`).
- **Exact even-star structure.** For `K!=0`, every nonempty even log-Walsh sector of an
  even star has strict alternating sign. At `0<v=tanh(K)<1`, the degree-six weight is a
  pure-even matchgate with ordinary CP rank two, all flattenings rank two, and Walsh-basis
  nonnegative CP rank `32`. The normalized checkerboard contraction is the cycle-space
  polynomial; no global 3D Pfaffian follows (`H646`-`H649`, correction `H655`, `e244`).
- **Selected-row Callen block classified.** At `0<v=tanh(K)<1`, the fixed contact-free
  selected-pivot/far-mark proper-even/odd leakage rank is
  `2^(d-1)-C(d,d/2)/2`; every closed combination has zero pair coefficient. Cubic inversion
  yields six exact higher-template relations in a `19x16` rank-13 block, but no pair closure
  (`H650`-`H652`, `e246`).
- **Verification.** Producers passed `13/13`, `9/9`, and `10/10`; all three clean-room
  verifiers passed after adversarial normalization, path-scope, provenance, and false-pass
  repairs. Separate theorem reviewers returned PASS. No endpoint moved and no full
  166-script suite was launched.

### Wave 23

- **Every-coupling full-spectrum no-go on the minimal open bipartite grid.** The three e238
  roots have unique exact coefficient sextics. Their Hermite signatures are `2,2,6`; the
  third has shifted sign pattern `+---+++` and at most two roots above `u=4`. Therefore the
  positive-definite open `2x3` transfer spectrum is not a full six-mode subset-product
  spectrum for any `0<t<1` (`H635`-`H639`, `e242`).
- **All-pivot termwise Callen closure is exponential.** On every connected six-regular simple
  graph, the nonzero one- and three-neighbour terms generated from one pair force exactly all
  even proper supports. For periodic cubic boxes, removing the empty normalization still leaves
  at least `ceil((2^(n-1)-[n even]-1)/n)` nontrivial translation orbits
  (`H640`-`H642`, `e243`).
- **Verification.** The e242 producer passed 10/10 in 814.373944 process CPU seconds. Its final
  current-source independent verifier passed the affine-Macaulay/cofactor stage in 713.290031
  process CPU seconds and the dyadic/Hermite/Descartes stage in 192.412972 process CPU seconds.
  The verifier-only bracket-alignment repair is explicit in artifact provenance; the producer
  was not rerun because the mathematical payload did not change. The e243 producer passed 8/8
  and its independently reconstructed closure verifier passed after adversarial repairs.
  No critical endpoint moved and no full 163-script suite was launched.

### Wave 22

- **Generic full-spectrum theorem on every connected simple nonpath graph.** Low-temperature
  valuations reduce to graph cuts; the Boolean cut criterion is exactly forest. High-temperature
  splittings reduce to unsigned token graphs; the additive-compound signing is balanced exactly
  for paths. Every nonpath graph is excluded near high temperature, cyclic graphs also near low
  temperature, and Chevalley makes the full-cube coupling set finite (`H600`-`H608`).
- **Finite open `2x3` continuum compressed to three algebraic trace candidates.** A stable
  rank-eight quotient and degree-971 norm polynomial leave exactly three physical `t` brackets.
  All are genuine trace branches and remained semialgebraically undecided at the wave-22 close;
  wave 23 above disposes them (`H628`-`H632`, `e238`, superseded in scope by `H635`-`H639`,
  `e242`).
- **New positive exact representations plus sharp limitations.** Checkerboard decimation is an
  exact octahedral 2-/4-/6-spin model with `c2>0,c4<0,c6>0`, not a visible pairwise closure
  (`H596`-`H599`). The complete Callen system reconstructs every finite free-boundary Gibbs
  vector but uses `2^n` correlators and therefore supplies no compression (`H614`-`H617`,
  correction `H623`).
- **Structural finite/all-size extensions.** Two-slice W-law injectivity is exact through `L=9`,
  with a `133x133` core showing raw singleton peeling cannot prove all `L` (`H610`-`H613`).
  Uniform phase-free graph-H1 twist averaging fails at girth on every cyclic simple graph
  (`H618`-`H622`). Open-box genus has the exact lower bound and planarity boundary recorded in
  `H624`-`H627`.
- **Verification.** All eight targeted fronts passed standalone lead verification. The degree-971
  verifier passed as separate elimination (`433.85` wall seconds) and roots (`110.00`) stages
  on the final source after an exact one-shot Sturm path exceeded 1,800 seconds; this is an
  observed resource wall, not a mathematical failure (`H633`).
  The SAW compact verifier additionally had both digest
  comparisons promoted into its pass predicate and reran PASS (`H609`). No critical endpoint
  moved and no full 161-script suite was launched.

### Wave 21

- **Finite bipartite spectral closures.** Exact determinant-centred log moments exclude full
  subset-product spectra on the open `2x2` and `2x3` grids throughout
  `I=[1/3-10^-9,1/3+10^-9]`. Independently, seven exact traces exclude the open `2x3` grid at
  `t=1/3`, with trace seven first decisive in that named resultant system. These do not prove
  every-coupling or all-size bipartite failure (`H584`-`H585`, `e228`, `e232`).
- **W-law mechanism scoped exactly.** The `Delta in {0,2}` slices observe every sector rank at
  `L=3..8`. The bonding/antibonding truncation retaining at most two `a` factors is closed by its
  first kernel at `(4,2)`, `Q a_0a_1a_2a_3`, rank `19<20`. The all-`L` rank law remains open
  (`H586`, `e230`).
- **Phase-sensitive upper class closed without endpoint movement.**
  `PS4`-convolution-block-Gram-L1 strictly strengthens `MR4`, but for every
  `0<K<=I_3/2` the exact Green sequence lifts on every even `L>=266`, has
  `p_0<=5168/(525L)`, and leaves the certified upper endpoint `I_3/2` unchanged (`H587`, `e231`).
- **Compact finite-memory SAW crosswalk landed without endpoint movement.** For every even
  `k=4,6,8,10,12`, the packed first-use automaton matches the ordered `e227` state and transition
  content; counts are `3/7`, `20/69`, `205/805`, `2722/11074`, and `41424/169975`.
  The producer passed `33/33` within the bounded cache/resource gate, and the independent verifier
  passed using distinct canonicalization. Memory 14 was explicitly not launched, so the certified
  interval is unchanged (`H590`-`H595`, `e229`).
- **Verification.** Producer artifacts passed `17/17`, `5/5`, `16/16`, `29/29`, and `33/33`;
  the lead observed PASS from all five targeted verifiers, including 21 checks and 0 failures for
  the W-law verifier. The compact SAW verifier proves only the `k<=12` crosswalk; it does not
  produce a `k=14` graph, certificate, or endpoint.

### Wave 20

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

## 5. Ranked queue after wave 24

1. **Empty the bounded bipartite exceptional schemes.** The new Heintz bound controls their
   cardinality but is far too crude to prove emptiness. Use shifted Stieltjes/localizing
   `H0/H1` minors of the mode polynomial, or a `2xL` replica-column recurrence. The smallest
   bounded pilot is open `2x4`; do not rerun open `2x3` or add same-point tests.
2. **Prove the W-law support-defect block theorem.** The 17-core is `16` complement-packing
   `(6,2)` rows plus one `(8,0)` row; each 133-core is the corresponding `128+5` boundary
   stratum. Seek an exact support-defect/empty-rung filtration, signed incidence complex, or
   uniform block determinant. Raw singleton peeling and naive `D^2=0` Koszul maps are closed.
3. **Enlarge the Callen block past the pair-row obstruction.** The fixed contact-free far-mark
   block has six higher-template relations but no pair-bearing dependency. A live route must
   add contacts, the full-neighbour row, radius, additional pivots, cross-orbit aggregation,
   or genuinely different nonlinear/auxiliary observables, with its own exact verifier.
4. **Exploit the local matchgate globally without hiding topology.** Every even star has the
   alternating log-sign law and the degree-six weight is locally Gaussian, but a nonplanar
   contraction is not thereby one Pfaffian. Locate the first crossing/spin-structure defect,
   build the sector-resolved surface bridge, or prove a controlled resummation. Any positive
   Walsh-basis auxiliary must beat the exact 32-atom barrier rather than undo decimation.
5. **Two-generator grids beyond ladders.** Prove the physical inequality for all `3xL` and then
   rectangles with `min(a,b)>=4`. Rank-one frozen rows and `|+>` compression are closed by the
   all-size depth-three anomaly; use a higher-rank invariant code or direct word family.
6. **Actual memory-14 finite-walk launch.** Compact `k<=12` crosswalks are exact, but no `A_14`,
   spectral radius, Collatz vector, or endpoint exists. A claim must construct all four, replay
   the componentwise integer inequality, and direction-round `atanh`; launch only after load clears.
7. **The upper endpoint.** `PS4` and earlier relaxation classes have zero uniform floor. A live
   certificate must use full cross-channel fourth moments, higher localizers, or multi-edge/sourced
   random currents that exclude the exact Green lift.
8. **Kac-Ward and surface routes beyond the closed classes.** Uniform phase-free H1 averaging
   fails at girth, fixed face/coordinate phases fail on the cube, and box genus grows. Remaining
   routes are fixed/nonuniform spin-structure sums, genuine surface 2-cell data, and branch
   `11111` outside the unit census. Never turn genus or `4^g` construction size into an arbitrary
   Pfaffian lower bound.
9. **`s<=6` extensive charges — structural idea required.** `s<=5` is exactly quotient one;
   both symmetry-compression routes are closed. Do not launch a larger version of either.
10. **Resource frontiers:** HT `v^30`, LT `x^60/x^62`, and the `3x4` pair-product certificate.
    These extend exact data but are not structural breakthroughs; cost every stage and wait for
    host clearance.
11. **Forced interlayer data, genuine spectral/IRF Lax equations, and Lee--Yang interlacing.**
    `W8/W10` need actual square-lattice evaluation/bounds; Lax progress needs a concrete equation;
    Lee--Yang needs a no-earlier-root or argument-principle theorem, not larger fits.

## 6. Environment state (measured 2026-08-24, 16:49 local)

- Host: 28 cores, 96 GB; load average measured **97.74 / 98.19 / 94.98** with nine users.
  CPU contention is severe. Keep heavy work single-process, nice/low priority, under 2 GiB,
  and budget by `time.process_time()`, never wall clock. Heavy series and memory-14 remain
  ENVIRONMENT-BLOCKED.
- Darwin `ru_maxrss` inherited a parent high-water mark during wave 24. New bounded producers
  use `mach_task_basic_info.resident_size_max`; this measures the task itself and keeps the
  two-GiB gate meaningful.
- Wave 24 leaves 166 auto-discovered test scripts. All three new targeted fronts and independent
  theorem reviews passed; no full 166-script suite was launched.

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
- **Bound the reduced incidence variety in ambient dimension, then audit the projection scope.**
  Homogenization gives `d^(n+2)` independently of the `2^n` equation count, but only after the
  inherited complex projection is finite. Path controls are cofinite and must stay excluded.
- **Matrix rank is not nonnegative tensor rank.** The octahedral tensor has rank two under every
  flattening yet nonnegative all-leg CP rank `32`; basis, sign, and atom-support constraints are
  load-bearing.
- **Diagonalize local XOR kernels before materializing supports.** The Callen neighbour cube turns
  one exponential-looking block into exact Walsh eigenvalues and exposes both its relations and
  its pair-row obstruction without generating the global support closure.
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
