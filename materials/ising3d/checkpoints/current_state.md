# Current state

**Snapshot:** 2026-08-24, wave 26 landed the first nonzero open-`2x4` trace-nine projection; the targeted producer and clean-room verifier pass.

Frozen Wave-26 plan and decision contract: `checkpoints/wave26_research_plan.md`.

- **Wave 26 (landed, one standalone producer/verifier front run by the lead):**
  1. **Trace nine now proves a finite open-`2x4` incidence.** The e248 `F4`-`F8` quotient has
     Hilbert vector `(1,5,12,19,22,19,12,5,1)`. Exact `96x96` multiplication-by-`F9`
     determinants are nonzero at `q=2` modulo `2147483647` and at `q=5/3` modulo
     `2147483629`; the synthetic mode polynomial with roots `4,...,11` gives rank 95 and
     determinant zero. The no-infinity leading system and multihomogeneous clearing therefore
     confine every full eight-mode incidence to a nonzero `q` polynomial of degree at most
     `182400` (`H669`, `e251`).
  2. **The projection is not yet empty.** Seven pointwise Lucas-product necessities hold strictly
     for the literal open-`2x4` targets after exact `q=1+x` coefficient tests, so none supplies a
     q-only shortcut. The primitive norm, all `q>1` roots, repeated-root-safe shifted
     `H0/H1` disposition, and open-`2x4` emptiness remain unresolved (`H670`).
  3. **Verification and allocation.** The producer passed ten hard gates in 86.942454 process CPU
     seconds and 174,309,376 bytes. The clean-room verifier used direct resultants, an independent
     reducer, and scalar elimination; it passed in 61.230270 seconds and 133,529,600 bytes. Three
     adversarial reviewers returned final PASS after repairing one displayed Lucas identity,
     denominator-scalar/count wording, q-residue guards, and finite-graph scope (`H671`).
     Experiments through `e251` and ledger through `H671` are allocated; `e252` and `H672` are
     next. The tree has 171 test scripts. No full suite, primitive-norm interpolation, endpoint
     movement, or thermodynamic calculation ran.
  4. **Bundle.** `deliverables/ising3d_wave26_final_bundle.zip` has 28 staged files plus manifest;
     all 29 members and every staged hash validate. Size 1,124,173 bytes, SHA-256
     `d07e9130b19c8cd5d4f1a95b71fb11d347e4c3b2edf73240ac50d9d48bb33dc1`.
     Snapshot docs predate this archive self-reference.


Frozen Wave-25 plan and decision contract: `checkpoints/wave25_research_plan.md`.

- **Wave 25 (landed, four standalone producer/verifier fronts run by the lead):**
  1. **First all-coupling boundary matchgate defect.** On the smallest nonplanar open cubic box
     `2x3x3`, eliminate either checkerboard color, retain the unique five-site opposite-color
     `3x3` face, and average the other four retained spins. The normalized five-leg Walsh tensor
     fails all five four-leg principal-Pfaffian equations in every one of `5!=120` ordinary
     unsigned leg orders for every `0<v=tanh(K)<1`. Four exact positive factor classes prove the
     sign; an alternating planar `C10` control passes exactly its ten dihedral orders. This is a fixed-boundary,
     fixed-basis defect, not a global Pfaffian-count or hidden-auxiliary no-go (`H658`-`H659`,
     `e247`).
  2. **Open `2x4` trace nine became exact and correctly dimensioned.** A `4^k`-state
     replica-column recurrence gives `H_k(q)=tr(B(q)^k)` for every open `2xL`; on `2x4` it
     constructs `H1` through `H9` with maximum 262,144 states and degree 162. The eight-mode
     rows have degrees `2,2,2,3,4,3`; `F4`-`F8` have a parameter-independent
     zero-dimensional leading ideal with 35 grevlex basis elements and 96 standard monomials.
     Trace eight is underdetermined once `q` is free; the `96x96` `F9` norm and physical branch
     disposition remain unresolved, so no `2x4` emptiness theorem is claimed (`H660`-`H662`,
     `e248`).
  3. **W-law cores now have an exact support-defect block mechanism.** The all-size elementary
     identity `|E|-|D|=L-2m` gives `|E|+|D|>=|L-2m|`. The literal e237 cores are exactly the
     boundary families `17=16+1` and `133=128+5`; particle-hole identifies the two `L=9`
     matrices. Each `133x133` core splits into four `29x29` shells plus the inherited signed
     `17x17` `L=8` core, explaining `3^5*2^288`. Smith torsion `6` in the first normalized shell
     exactly refutes the tempting eight-copy/dyadic recurrence. No `L>=10` rank or residual
     theorem follows (`H663`-`H666`, `e249`).
  4. **Contact-inclusive Callen aggregation still cannot close the pair.** On the finite
     `C3^3` periodic cubic graph, all 864 radius-one contact rows fold to four classes under
     translations and full cubic symmetry. After allowing normalization and the nearest-neighbour
     pair, a nonallowed `4x4` minor is
     `4(c1-c5)(c1+4c3+c5)` and is strictly positive for every `0<v<1`; hence no nonzero folded
     combination closes on the pair. Contacts, all pivots, and the five-neighbour/full-neighbour
     row are included; radius, nonlinear, auxiliary, and nonsymmetric zero-average relations are
     not (`H667`-`H668`, `e250`).
  5. **Verification and allocation.** Producers `e247/e248/e249/e250` passed `13/9/18/10`
     hard gates; all four clean-room verifiers passed. The largest bounded stage was the e248
     verifier at 21.271467 process CPU seconds and 1,245,003,776 bytes peak task RSS; its producer
     used 18.714901 seconds and 1,241,841,664 bytes. Adversarial review
     returned PASS on all four fronts after hardening ordinary-order scope, actual trace degrees,
     Lucas/mode linkage, semantic hashes, dependency locks, and Callen pivot provenance.
     Experiments through `e250` and ledger through `H668` are allocated; `e251` and `H669` are
     next. The tree has 170 test scripts. No full suite or memory-14 launch ran, and the critical
     interval is unchanged.
     Definitive Wave-25 bundle: 46 staged files plus manifest; ZIP 1,575,644 bytes,
     SHA-256 `d167e544c7693b10d9526f66987da5a28cb71d610d5f3f4819831bf525f5a29e`;
     every staged hash and all 47 members validate. Snapshot docs omit archive self-references.
     The two earlier Wave-25 ZIPs are superseded and retained pending user-approved deletion.


Frozen Wave-24 plan and decision contract: `checkpoints/wave24_research_plan.md`.

- **Wave 24 (landed, three standalone producer/verifier fronts run by the lead):**
  1. **Graph-uniform exceptional-scheme degree bound.** For every connected simple nonpath
     graph with `n` vertices and `m` edges, the physical full-spectrum exceptional set has at
     most `[2^n(2n+2m)]^(n+2)` points; uniformly at fixed `n`, at most
     `[2^n(n^2+n)]^(n+2)`. The polynomial representative has entry degree `2n+2m`.
     Paths are excluded, and the bound does not prove emptiness (`H643`-`H645`, `e245`).
  2. **All-even-star sign law and local matchgate structure.** For even star degree `d`, every
     nonempty log-Walsh sector of size `2r` has strict sign `(-1)^(r+1)` for `K!=0`.
     At `0<v=tanh(K)<1`, the degree-six weight signature is the pure-even
     sub-Pfaffian tensor of `v^2 J_6`.
     It has ordinary CP rank two and every flattening has rank two, but nonnegative all-leg
     CP rank exactly `32`; the positive spin-basis rank-two form merely restores the eliminated
     centre. The normalized checkerboard contraction is the positive cycle-space polynomial,
     not a global 3D Pfaffian theorem (`H646`-`H649`, `e244`).
  3. **Selected-row Callen Walsh theorem.** For even local degree `d=2m` and
     `0<v=tanh(K)<1`, the proper-even by odd
     leakage matrix has rank `2^(d-1)-C(d,m)/2` and left nullity `C(d,m)/2-1`; every
     leakage-canceling combination has zero pair-row coefficient. Cubic inversion leaves six
     exact higher-template relations in a `19x16` block of rank `13`, while deleting the pair
     row lowers the rank to `12`. This closes only the fixed contact-free far-mark block
     (`H650`-`H652`, `e246`).
  4. **Verification and allocation.** Producers `e244/e245/e246` passed `13/13`, `9/9`, and
     `10/10`; all three clean-room verifiers passed after review-driven normalization,
     path-scope, provenance, and false-pass repairs. Independent theorem reviewers returned
     PASS. Append-only corrections `H653`-`H657` supersede the pre-hardening scope,
     normalization, timing, resource-provenance, and external-source wording in
     `H644/H645/H649/H654/H655`. Experiments through `e246` and ledger through `H657`
     are allocated; `e247` and `H658` are next. The tree has 166 test scripts.
     `e233` remains an
     unlanded draft.
     The critical interval is unchanged, and no full 166-script suite was launched under the
     contended multi-user host.
     Wave-24 bundle: 43 staged files plus manifest; ZIP 814,269 bytes,
     SHA-256 `b716c320aa6e253dc26345497ff512ff94949fd2ae5cc5c3bff3666ce41f7608`;
     every staged hash and all 44 ZIP members validate. Bundle snapshot docs intentionally
     predate this self-referential archive record.

- **Wave 23 (landed, both standalone verifier fronts run by the lead):**
  1. **Open `2x3` every-coupling full-spectrum theorem.** The three e238 real roots have unique
     exact coefficient sextics from rank-seven adjugate certificates. Their Hermite signatures
     are `2,2,6`; the third has shifted signs `+---+++` and at most two roots above `u=4`.
     Therefore no physical `0<t<1` gives a full six-mode subset-product spectrum on this finite
     layer (`H635`-`H639`, `e242`).
  2. **All-pivot termwise Callen closure theorem.** On every connected six-regular simple graph,
     nonzero one- and three-neighbour terms generated from one pair force exactly all even proper
     supports. On periodic cubic boxes, removing the empty normalization still leaves at least
     `ceil((2^(n-1)-[n even]-1)/n)` nontrivial translation orbits. Selected-pivot schemes, row
     combinations, cross-orbit aggregation, and nonlinear/auxiliary variables remain open
     (`H640`-`H642`, `e243`).
  3. **Verification and allocation.** The e242 producer passed 10/10; its final current-source
     quotient/cofactor stage passed in 713.290031 process CPU seconds and its
     dyadic/Hermite/Descartes stage in 192.412972 process CPU seconds. The e243 producer passed
     8/8 and its independent verifier passed after adversarial repairs. Experiments `e242` and
     `e243` are landed; `e244` is next. Ledger `H635`-`H642` is allocated; `H643` is next.
     The tree has 163 test scripts. `e233` remains an unlanded draft. The critical interval is
     unchanged and no full 163-script suite was launched.
     Wave-23 bundle: 28 staged files plus manifest; ZIP 1,333,269 bytes,
     SHA-256 `fcbb60b2bf0797b7141193110e2d78d411221ab6942e90b1796e024a2a60e57a`;
     staged hashes and all 29 ZIP members validate. Bundle snapshot docs intentionally predate
     this self-referential archive record.


- **Wave 22 (landed, all eight accepted standalone verifier fronts run by the lead):**
  1. **All-size full-spectrum endpoint and generic theorem.** At `t->0`, ordered transfer
     valuations are the graph cut distribution; this is a Boolean subset sum iff a simple graph
     is a forest. At `t->1`, the first graph-sensitive bands are unsigned token graphs, while a
     full Gaussian requires additive compounds; their signing is balanced iff a connected graph
     is a path. Thus every finite connected simple nonpath graph is excluded near high
     temperature, every cyclic simple graph near low temperature, and the full-cube coupling set
     is finite by Chevalley constructibility. Exceptional couplings are not located or proved
     absent in general (`H600`-`H608`, `e235`-`e236`).
  2. **Open `2x3` continuum reduced to three exact candidates at wave 22.** A stable rank-eight
     quotient gives a primitive squarefree dense `E(q)` of degree `971`. Exactly three `q>1`
     roots map to the stored physical intervals. They were genuine seven-trace algebraic branches;
     wave 23 now disposes all three semialgebraically (`H628`-`H632`, `e238`, superseded in scope
     by `H635`-`H639`, `e242`).
  3. **Exact octahedral checkerboard decimation.** Eliminating a degree-six star gives
     `c2>0,c4<0,c6>0` for every real `K!=0`; arbitrary visible fields/pairs cannot reproduce
     the positive six-spin coefficient. Checkerboard elimination is exactly a 2-/4-/6-spin
     octahedral model, not a pairwise closure (`H596`-`H599`, `e234`).
  4. **Exact Callen reconstruction.** The six-neighbour local field has explicit
     `c1(v),c3(v),c5(v)`. On every finite loopless undirected free-boundary graph with uniform
     internal coupling and no field, the complete Callen system has formal rank `2^n-1` and
     reconstructs the unnormalised Gibbs moments up to scale. It remains exponential and gives no
     thermodynamic compression (`H614`-`H617`, correction `H623`, `e240`).
  5. **W-law finite exact sequence extended through `L=9`.** Every one of 49 sectors at
     `L=3..9` has exact-Q two-slice injectivity. The balanced `L=9` rows use 1,556 integral leaf
     pivots plus a `133x133` determinant `-2^288*3^5`; raw singleton peeling is thereby closed as
     an all-`L` proof method. The all-`L` lemma remains open (`H610`-`H613`, `e237`).
  6. **Uniform flat-holonomy averaging closed all sizes.** Averaging phase-free Hashimoto
     determinants over every graph-`H^1` twist kills each shortest cycle: for every connected
     simple cyclic graph the first mismatch with `P_G(v)^2` is `0` versus twice the girth-cycle
     count. All 32 open-cube classes and two finite phase repairs were decided exactly
     (`H618`-`H622`, `e239`).
  7. **Open-box topology.** For `min(a,b,c)>=2`,
     `gamma(P_a square P_b square P_c)>=max(0,ceil((abc-ab-bc-ca+4)/4))`; an open box is planar
     exactly when a side is one or it is `2x2xL` up to permutation. This is not an arbitrary
     Pfaffian-term lower bound (`H624`-`H627`, `e241`).
  8. **Verification and hygiene.** Producer checks are `14/14`, `9/9`, `7/7`, `13/13`,
     `9/9`, `14/14`, `13/13`, and `12/12` for `e234`-`e241` in experiment-number order;
     every independent verifier passed. The final trace verifier passed as affine-Macaulay
     elimination (`433.85` wall seconds) and hardened root-certificate (`110.00`) stages after a
     one-shot Sturm path exceeded its resource wall; `H633` records the final-source correction.
     The compact SAW verifier now asserts both content digests and reran PASS (`H609`).
     At the wave-22 close, ledger `H596`-`H634` was used, `H635` was next, experiments
     `e234`-`e241` were landed, `e242` was next unused, and 161 tests were present. The wave-23
     block above supersedes that allocation. `e233` remains an unlanded draft. The critical
     interval is unchanged.
     Wave-22 bundle: 47 staged files plus manifest; ZIP 1,248,489 bytes,
     SHA-256 `5767bfed50a02cdbeddebb2b686f6f55e5657d25b4bcfd9a623a39b38a13ac7c`;
     staged hashes and all 48 ZIP members validate.

- **Wave 21 (landed, all five targeted standalone verifiers run by the lead in this checkout):**
  1. **Finite bipartite log-moment interval.** On
     `I=[999999997/3000000000,1000000003/3000000000]=[1/3-10^-9,1/3+10^-9]`,
     the open `2x2` and `2x3` positive reciprocal spectra violate the necessary
     Stieltjes/Hankel inequality `p1*p3-p2^2>=0`; artifact `17/17`, lead verifier PASS
     (`H584`, `e228`).
  2. **Trace-seven resultant obstruction.** At `t=1/3`, the physical open `2x3` layer is
     not a full six-mode subset-product spectrum; traces four through six leave proper ideals
     and trace seven first gives `{1}` in the named elimination; artifact `5/5`, lead verifier
     PASS (`H585`, `e232`).
  3. **Finite W-law rank mechanism and truncation closure.** The `Delta in {0,2}` slices
     observe every exact sector rank for `L=3..8`; the at-most-two-`a` bonding/antibonding
     truncation first fails at `(L,m)=(4,2)` with kernel `Q a_0a_1a_2a_3` and rank `19<20`.
     Artifact `16/16`; the lead verifier reported 21 checks, 0 failures, PASS (`H586`, `e230`).
  4. **Phase-sensitive upper class closed at its exact scope.**
     `PS4`-convolution-block-Gram-L1 is strictly stronger than `MR4` at `L=4,K=6/25`, but
     for every `0<K<=I_3/2` the exact Green profile lifts on every even `L>=266` with
     `0<p_0<=5168/(525L)`.
     Its uniform floor is zero and the upper endpoint does not move; artifact `29/29`, lead
     verifier PASS (`H587`, `e231`).
  5. **Compact finite-memory SAW crosswalk.** For each even `k=4,6,8,10,12`, the packed
     first-use automaton matches the ordered `e227` state and transition content, with exact
     state/transition counts `3/7`, `20/69`, `205/805`, `2722/11074`, and `41424/169975`.
     The bounded producer passed `33/33`, and the distinct-canonicalization verifier passed.
     Memory 14 was explicitly not launched and the critical interval is unchanged
     (`H590`-`H595`, `e229`).
  6. Ledger `H584`-`H587` and `H590`-`H595`; `H588`-`H589` remain unallocated. Experiments
     `e228`-`e232` are used, 153 test scripts are present, the next free ledger ID is `H596`,
     and the next free experiment is `e233`.
- **Wave 20 (landed, all thirteen standalone verifiers run by the lead in this checkout):**
  1. **Critical floor moved.** Exact memory-12 finite-walk automaton: 41,424 states,
     `1000000 A w < 4747526 w` with minimum integer residual `18,990,104`, hence
     `mu<4.747526` and `K_c>0.2138368062108697304302094293879379891333` (`H578`-`H583`).
  2. **Connected-bipartite local-term classification completed.** Every `Delta>=3` graph is the
     nonradical `Q_c=1` branch, without a Hamiltonian hypothesis; full archived Seven source
     plus independent census (`H572`-`H577`).
  3. **Every-coupling isotropic theorem on the non-bipartite locus.** The odd-Eulerian trace
     invariant excludes every full subset-product spectrum for every `0<t<1`; it is identically
     blind on bipartite graphs (`H524`-`H529`).
  4. **New exact structures:** `K4,4=137`, `K4,5=471`; invariant-theory explanation of the
     conjectured `W_L` defect; all-`v` `c8` structure forces `W8`, and `c10` forces `W10`.
  5. **Exact route closures:** rank-one/`|+>` rectangle induction, Galois nonsolvability as an
     integrability discriminator, intrinsic finite 3D CTM spectrum, local-algebra-to-all-aux-Lax,
     `MR4-power-simplex-L1`, and susceptibility supermultiplicativity.
  6. Ledger `H506`-`H583`, experiments `e191`-`e227`; verifier review caught and corrected the
     unique coefficient-ratio plateau `5,5` before integration.
- **Wave 19 (landed):** local-term non-bipartite dichotomy; two-generator counterexamples
  `K2,3=44`, `K3,3=63`; all-`L>=3` ladder fixed-basis quadratic no-go; exact closure of the
  `ad_A` grade-separation and varying-conjugator rigidity routes (`H488`-`H505`).
- **Wave 18 (landed, all fronts independently re-verified by the lead in this checkout):**
  1. **All-`L` classification.** For EVERY `L>=2`, `g_{2xL}(Q) = so_m + so_m` (even `L`) or
     `sp_m + sp_m` (odd `L`), `m = 2^(n-1)`, `n = 2L`, `dim = 2^(n-1)(2^(n-1)-(-1)^L)`. The
     saturation step open since wave 16 is closed by an explicit production induction. Full closure
     enumerated at `L=2..6`, reaching `4192256` at `L=6` with zero set-law violations; `L=7` is
     certified as `67117056` by construction plus an independent Clifford census and is explicitly
     not a full enumeration (`H428`-`H433`).
  2. **Exact trichotomy.** For bipartite `Gamma` on `n` vertices with a Hamiltonian path: path gives
     `dim = n(2n-1)`; even cycle gives `2n(2n-1)`; `Delta >= 3` gives
     `2^(2n-2) - (-1)^(n/2) 2^(n-1)` (even `n`) or `2^(2n-2) - 1` (odd `n`), depending only on `n`.
     `C_6` at dim `132` is the exact counterexample fixing the hypotheses; `3x3 = 65535` is the
     odd-`n` instance, so it is not a separate regime (`H480`-`H487`).
  3. **All-`L` quadratic no-go**, mode floors `6, 18, 65, 216, 991, 3168, 15354` (`H434`-`H439`).
  4. **Isotropic gap split verdict**: three cores certified above ceiling, decoupling provably fails
     (`H440`-`H445`).
  5. **Theorem U2 restored**: exact char-0 `r(T)=417`, `r_all(T)=445`; F2 unconditional for `m=2,4`
     (`H446`-`H451`).
  6. **Kac-Ward**: no orbit reduction exists (16384 singleton orbits) and the complete `F_5`-unit
     census over `4^13` points kills all 2960 local points; branch `11111` still open over `Q`
     (`H452`-`H457`).
  7. **Upper endpoint**: two named certificate classes exactly refuted; incumbent unchanged
     (`H464`-`H469`).
- **Wave 17 (closed):** 2x6 spectral non-Gaussian certificate, the largest layer so certified
  (H425: pair-poly gcd <= 295414 < 295415 on both primes; H426: independent replay of the saved
  tail; H427: wording corrections). 2xL alternation: containment for ALL L + equality at L=2..5
  with form/type/exact closures (H419-H424; H422/H423 are review corrections). Theorem-PM shell
  family PM-S closed at Lambda = 2,3,4 (H417). LT series exact to x^56 = -979227570369/4, both
  producer and fresh-prime audit verified (H418). Acceptance gate: `fullsuite17b` 123/123 green
  (8h10m under host contention), including the standalone 38-check H425 verifier.
- **Wave 16 (landed earlier this session):** 2x4 char-0 DLA `so_128(Q) (+) so_128(Q)` exact
  (H411/H412); the 2xL alternation law bit-exact at L<=4 (H415/H416); exact LT x^54 =
  675410878105/9 (H413/H414); shell-route falsified-for-purpose (H404/H405); annihilator
  horizon L0~14 (H409/H410); s<=6 compression-class closed (H406); suite gate 115/115
  pre-wave-16.

**Repository:** `/Users/jinleic/jinleic-workspace/math/ising3d`  
**Interpreter:** `.venv/bin/python`  
**Research status:** active exact-research program; **the 3D Ising model is not solved**.

## Cold-resume orientation

The repository studies the simple-cubic nearest-neighbour Ising model using exact finite-volume arithmetic, the 2D Onsager/Kaufman solution as a control, exact finite-lattice series, rigorous critical bounds, and scoped algebraic/integrability obstructions.  Begin with:

1. `problem_specification.md` — lattice, boundary, and normalisation conventions.
2. `notes/CODE_API.md` — stable code APIs and Pauli representation.
3. `README.md` — result map and reproduction commands.
4. `checkpoints/verified_results.json` — machine-readable values and provenance.
5. `checkpoints/failed_routes.md` and `notes/hypotheses.csv` — negative results and corrected quantifiers.
6. `checkpoints/next_actions.md` — ranked open work.

The reusable code is under `src/ising/`; standalone computations are under `experiments/`; machine records are under `results/`; mathematical scopes are under `proofs/` and `notes/`.

## Strongest established facts

- **Finite arithmetic:** independent exact enumeration and transfer propagation agree on 17 lattice/boundary cases (`tests/test_tm_vs_enumeration.py`).  The natural layer transfer operator reproduces exact 3D partition functions on seven lattices (`tests/test_transfer_operator_identity.py`).
- **2D control:** Kaufman's finite-torus formula agrees with exact integer transfer matrices in eight cases to maximum relative error `3.63867276317e-60`; parity-resolved `log V` has maximum off-Majorana-bilinear coefficient `1.246e-15` (`results/onsager/onsager_2d.json`).
- **Exact series:** simple-cubic reduced free energy reaches `v^28` at HT (`proofs/ht_v28.md`) and `x^56` at LT (`proofs/lt_x56.md`, superseding the `x^34`/`x^52` stages in `proofs/lt_x34.md` and `proofs/lt_x54.md`).  The generalized cap-vector frontier reproduced `v^24`/`v^26` before giving `[v^28]phi=2102198327465307/28`, interaction `a_28=525549581866326/7`, with complete profile-wise CRT uniqueness and a post-computation Arisue--Fujiwara witness; the LT side reaches `x^56 = -979227570369/4` with per-box CRT certificates, an independent fresh-prime audit, and a post-computation Guttmann--Enting witness.  No D-finite recurrence at budget ≤11; complete second-order/linear-ODE/Mahler negative frontiers remain finite-budget results.
- **Rigorous critical interval (wave 20; unchanged in wave 21):** `K_c` belongs to
  `[0.2138368062108697304302094293879379891333, 0.2527310098586630030260020266135701299926]`.
  The lower endpoint is `atanh(1000000/4747526)` rounded strictly downward after the exact
  memory-12 Collatz certificate; the independent verifier rebuilds all 41,424 states and replays
  all 169,975 transitions. Wave 21's compact packed-state crosswalk independently reproduces the
  same ordered automata through `k=12`; it launched no `k=14` computation and therefore does not
  alter either endpoint. This supersedes the wave-14 SAW-union floor. The upper endpoint remains
  method-optimal only for the audited two-point infrared/GKS class, not globally. Wave 21's
  strictly stronger `PS4` phase-sensitive class also admits the exact Green sequence for every
  `0<K<=I_3/2` along all even `L>=266`, with zero uniform floor and no endpoint movement
  (`proofs/kc_lower_finite_memory.md`, `proofs/kc_lower_memory14.md`,
  `proofs/upper_phase_sensitive.md`).
- **Lee--Yang:** finite field zeros pass unit-circle checks on 13 lattices through `4x4x4`; the largest stored 64-site radial residual is `9.1044993726259e-92` (`results/lee_yang/lee_yang_analysis.json`).
- **Duality:** exact finite-volume 3D Ising--`Z_2` gauge duality passes 27 integer checks; the simple-cubic Ising model is not self-dual (`proofs/duality_3d.md`).
- **Algebraic theorem:** for `A=sum_v X_v`, `B=sum_(ij) Z_iZ_j`, DG1 holds on every graph.  The exact DG2 residual is the degree-linear plus neighbor-triple formula in `proofs/dolan_grady_defect.md`.  It vanishes iff every degree is `0` or `2`; its quartic part vanishes iff maximum degree is at most 2.
- **Spectral no-go (Theorem S, wave 5):** the `2x3` layer transfer operator at `tanh(K*/2)=1/3` is NOT spectrally Gaussian on 6 modes in any Majorana basis: exact rational Sylvester-inertia certificate; hypotheses (three lowest eigenvalues distinct) certified with relative gaps 53.0%/40.9%; forced value window `[0.0347241721739, 0.0347241728052]` contains no eigenvalue (inertia 8=8).  Controls: window-occupancy consistency only, never exact presence (`proofs/spectral_gaussianity.md`).
- **Wave-21 scoped spectral and ladder closures:** determinant-centred log moments exclude full
  subset products for the open `2x2` and `2x3` layers on the exact interval
  `[1/3-10^-9,1/3+10^-9]`; an independent trace-resultant obstruction excludes the open `2x3`
  layer at `t=1/3`, with trace seven first decisive in that named system. For ladders, the
  `Delta in {0,2}` slices recover every sector rank at `L=3..8`, while the at-most-two-`a`
  bonding truncation first fails at `(4,2)`. No every-coupling/all-size bipartite theorem or
  all-`L` `W_L` rank law follows (`proofs/bipartite_log_moments.md`,
  `proofs/trace_resultant.md`, `proofs/wlaw_rank_mechanism.md`).
- **Algebra structure (wave 5, since resolved):** `C4` exact-Q decomposition `A1^3+C^2`; `2x3` certified mod-p quotient witnesses `sp(14), sp(6), gl(6)` (saturation-certified; the char-0 Levi was later settled exactly by wave 6 — see Theorem OA above and `results/algebra_structure/char0_levi.json`); `2x4` dim 2952, centre 1, derived 2951 (computed, not literal; the complete char-0 structure was later settled by waves 7-8 — see Theorem OA-2x4 above and `results/algebra_structure/char0_complete_2x4.json`).  `structure.json` keeps e29-scoped fields with `characteristic_zero_resolution` pointers.
- **Interlayer expansion (through wave 20):** exact `c2,c4,c6,c8` through `v^12`; `c8` now has
  an all-`v` connected-correlation formula and uniquely forces `W8^2/8!`; `c10` has 16 profiles
  and forces `W10^2/10!`; generally `c_(2m,q)` contains the square of the connected `2m`-point
  datum divided by `(2m)!` (`proofs/interlayer_c10.md`).
- **DG+TD simultaneous no-go (wave 5):** both varieties empty for the minimal orbit basis; Groebner `{p^2,pq,q^2,pr,qr,r^2}`, all saturations `[1]` (`proofs/dg_simultaneous.md`).
- **Theorem OA (wave 6, headline):** the characteristic-zero `2x3` layer algebra has exact structure `263 = 1 + [C7+C3+C3+A8+A5]` (two independent exact-Q closures, Killing rank 262, C7-vs-B7 separated by an exact alternating form det 16384); combined with the fulltext Date–Roan classification (Levi of any finite-dim OA quotient = `sl2^n`), the layer algebra is **not an Onsager-algebra quotient in any disguise** — no generating pair satisfies Dolan–Grady (`proofs/char0_structure.md`, `proofs/onsager_quotient_nogo.md` sec. 9).
- **Theorem S′ (wave 6):** the two physical `P=prod X` sector spectra of the `2x3` layer are not the even/odd halves of any single 6-mode Gaussian, either assignment (`proofs/parity_sector_gaussianity.md`).  Theorem S itself extended to `t=1/5,2/5,1/2` on `2x3` and to the `2x4` layer at `t=1/3` via certified symmetry sectors (`results/spectral/couplings.json`).
- **Graded tetrahedron no-go (wave 6):** every Z2-grade-preserving auxiliary `R` is zero for ALL `q` outside `{0,±1}`: `det B(q) = -2 q^56 (q-1)^45 (q+1)^39` (`notes/tetra_graded.md`).
- **Kac–Ward (wave 6):** branch `00000` empty over Q (integer certificate `2F-3F=56`); `F_7` witness has no lift mod 49; `F_5`/`F_11` height bounds `~10^29` against rational points; family still UNRESOLVED (`notes/kw_hensel.md`).
- **Ladder obstruction (wave 6):** all single-Pauli words in `<A_L,B_L>` are row-swap fixed and pairwise commuting for every `L`; the single-string injection strategy class is dead; `D_{2L}>=2^L` still open beyond `L=6` (`proofs/ladder_growth_family.md`).
- **Peierls attempt (wave 6):** honest non-improvement `K_c <= 1.696898` with new exact contour counts `N(6..24)` and a proved method-limitation (`proofs/kc_upper_peierls.md`); incumbent upper bound untouched.
- **Series structure (wave 6):** complete algebraicity + first-order differential-algebraicity negative frontiers (HT all-negative; LT survivors all truncation-unobservable); interlayer `[v^8]c2=778` (two routes); `c3(v)=0` identically by layer-flip symmetry (`notes/series_structure.md`).
- **Complete local-term graph classification (waves 18-20):** connected bipartite graphs are
  exactly path, even-cycle, or `Delta>=3` quadratic-root branches; no Hamiltonian hypothesis
  remains. Connected non-bipartite graphs are odd cycles or the full noncentral even branch.
  Every open `2xL` algebra is `so_m+so_m` (even `L`) or `sp_m+sp_m` (odd `L`). The TWO-SUM
  algebra is different: exact anchors include `11,263,2952,8034` on `2x2,2x3,2x4,3x3` and
  `44,63,167,137,471` on `K2,3,K3,3,K3,4,K4,4,K4,5`; every ladder `L>=3` clears the physical
  quadratic ceiling, but no all-grid or universal graph threshold is proved
  (`proofs/nonhamiltonian_bipartite.md`, `proofs/kmn_twogen.md`, `proofs/ladder_w8.md`).
- **Theorem OA-2x4 (waves 7-8, headline):** complete characteristic-zero structure `2952 = Qz + so(F_000) + so(F_010) + so(F_100) + sl24 + sl28 + 2 sl4`, absolute type `D21+2B13+A23+A27+2A3`, radical = centre, structural Killing rank 2951; splitness of the orthogonal Q-forms deliberately not claimed; clean-room no-producer-import verifier re-derives all 13 controlling minors.  Corollary with Date--Roan: **the 2x4 layer algebra is not an Onsager-algebra quotient** (`proofs/char0_complete_2x4.md`, `proofs/levi_images_2x4.md`).
- **Kac--Ward (waves 7-8):** positive-real cone EMPTY by exact Positivstellensatz on the free `2x2x2` cube (`F = 48 + 8 sum M`); rational branch map **31/32 EMPTY_OVER_Q** via thin-box integer constants; branch `11111` unresolved with exact negative certificates (no degree-0 span certificate, no multiplier-degree-1 Nullstellensatz; small-prime points persist) (`proofs/kac_ward_positive_real.md`, `proofs/kac_ward_close.md`).
- **Ungraded tetrahedron (wave 8):** arbitrary `8x8` auxiliary `R` vanishes for all char-0 `q` outside `{0,+1,-1}`; complementary minors `2q^109(q-1)^92(q+1)^77(q^3+2q-1)` and `2q^110(q-1)^92(q+1)^77(q^2+1)` with Bezout clearance (`proofs/tetra_ungraded.md`).
- **Full commutant program (waves 7-8):** exact joint commutants `2x2..2x6 = 27,12,41,20,71` with split types (`M2+C^8`, `M2^2+C^12`, ..., `M7+M2^2+C^14` — the `2x6` value closed exactly in wave 9); `3x3 = 98`; support-resolved vacuum through size 7 on `2x4`/`3x4` (6 on `3x3`); pointwise-stabilizer criterion proved with 2682-orbit zero-mismatch converse (`proofs/full_commutant.md`, `proofs/commutant_family.md`, `proofs/charge_criterion.md`).
- **Spectral extensions (waves 7-8):** Theorem S now covers `3x3` (t=1/3, 1/2; twelve refined sectors) and `2x5` (t=1/3, 1/2; 1024-dim, strict congruence lemma), plus four exact continuum coupling components on `2x3`; the crossing question was settled in wave 9: the crossings are real attained algebraic values, and no absence interval can span one; three crossing-free absence intervals are proved, while crossing uniqueness/exhaustiveness remains unproved (see the wave-9 spectral bullet above); the proposed degree-four polynomial shortcut was exactly falsified in-house and converted to a regression gate (`results/spectral/`, `proofs/parametric_gaussianity.md`, `proofs/spectral_polynomial.md`).
- **Series (waves 7-8):** HT extended to `v^24 = 2135670379057/8` exactly (cut-capped parity frontier, 352 MB vs the 242-GB projection; five-prime CRT uniqueness; reserved as strict holdout); complete second-order/linear-ODE/Mahler negative frontiers, 138 rows, zero unexplained survivors (`proofs/ht_v24.md`, `notes/series_frontier2.md`).
- **Interlayer (waves 7-9):** exact `c4`, `c6` through `v^12`, each by two independent exact routes (`c6` second route: wave-9 open-slab even-subgraph FLM over 112 boxes with a proved finite-support lemma — `proofs/interlayer_closed.md`); exact closed form for `c4_q` in 2D `G`/connected-`U4` data, derived and verified coefficientwise through `v^12`; all-order theory: evenness theorem with sharpness witnesses, composition selection rule, `[v^2 w^(2m)] = 2` for all `m` (`proofs/interlayer_allorders.md`).
- **3x3 layer algebra (wave 9):** complete char-0 structure — `g = Qz + sl33 + sl48 + sl59 + sl3 + sl16 + sl30`, dim 8034, rad = centre, type `A32+A47+A58+A2+A15+A29` (`proofs/char0_3x3.md`).
- **HT series (wave 9):** extended to `v^26 = 115377512914251/26` exactly; external Arisue--Fujiwara witness v^2..v^26; square-free-log bug found/fixed with exact controls; `v^25 = 0` (`proofs/ht_v26.md`).
- **Spectral crossings (wave 9):** no absence interval can span a certified crossing — the forced value is attained at real algebraic t* (four exact c-jump brackets); absence IS proved on three connected closed rational intervals (crossing-free), widest 1/1048576 (`proofs/spectral_interval2.md`).
- **2x6 commutant (wave 9):** exact dim 71, type `M_7 + M_2^2 + C^14`; M_7 = End(common zero space) (`proofs/commutant_2x6.md`).
- **SAW/bounds (wave 9):** first-backtrack theorem `c_(m+n) <= c_m(c_n - a_(n-1))`; exact LP floors at L=4,6; wider Simon--Lieb family exactly rejects improvement; endpoints unchanged (`proofs/kc_interval2.md`).
- **Ladder (wave 9):** `D_16 >= 256` restored by patched 256-word family; defect proved not symmetry-removable; stationary two-letter class exactly obstructed at L=8 (`proofs/ladder_all_l.md`).
- **Cabled tetrahedron (wave 9):** generic Q(q) rank 256/nullity 0 for the 16x16 cabled system + three exact rational certificates; exceptional locus open (`proofs/tetra16.md`).
- **Kac--Ward 11111 (wave 9):** no degree-<=2 Nullstellensatz certificate; 16-dim Laurent torus reduction; Q(i) mu4 chart exhaustively empty; branch open (`proofs/kw_branch11111.md`).
- **Lee--Yang (wave 9):** exact 5x4x4/5x5x4 zeros (width 2^-260 isolation, unit-circle certified); fixed-edge correction class excluded; 5^3 wall quantified 50.7 GB (`proofs/lee_yang_zeros.md`).
- **K_c interval (historical wave-8 state, superseded at the lower endpoint):** wave 8 left `[0.2122119011661678, 0.2527310098586630]`; wave 10 strictly raised the lower endpoint as recorded in the current interval bullet above.  The upper endpoint remains method-optimal only for the audited two-point infrared class, not globally.
- **Wave-10 algebra:** the `3x3` algebra is now also proved not to be any finite-dimensional Onsager quotient (`proofs/oa_3x3.md`).  For open `4x4`, 1794 stored literal Pauli words prove only `dim_Q g>=1794`; the closure remained incomplete at depth 16, so exact dimension/type/centre/radical/Killing rank are open (`proofs/char0_4x4.md`).  The ladder certificate reaches `D_18>=512` at `L=9`; no all-`L` recursion is proved (`proofs/ladder_nonstat.md`).
- **Wave-10 Kac--Ward:** branch `11111` has a nonsingular `F_5` construction point (`det J=2 mod 5`), hence a `Q_5` component and a proper rational ideal.  Therefore no rational Nullstellensatz emptiness certificate exists at any multiplier degree.  Its local Hensel lift fails the independent `3x3x3` order-eight holdout by residue `250 mod 625`; the full signed/complex branch remains open (`proofs/kw_cubic.md`).
- **Wave-10 tetrahedron:** the canonical one-line-cabled `16384x256` system now has complete characteristic-zero rank-drop locus `{ -1,0,1 }`; the complementary-minor gcd is `q^584(q-1)^580(q+1)^546(q^2+1)^4`, with `q=+-i` cleared as full-rank spurious minor roots (`proofs/tetra16_locus.md`).
- **Wave-10 spectral crossings:** each of the four inherited fine brackets contains exactly one crossing by exact five-variable Krawczyk contraction and branch transport.  Global crossing exhaustiveness on `[1/4,7/20]` and one-variable minimal polynomials remain open (`proofs/crossing_alg.md`).
- **Wave-10 interlayer:** `c6_q` has a complete Bell(6)/XOR reduction to 13 anchored products of solved-2D `G`, connected `U4`, and connected `W6`, independently checked through `v^12`.  The infinite `W6` lattice sums are not evaluated at criticality (`proofs/c6_closed.md`).
- **Wave-10 Lee--Yang:** a semilinear `D8`-plus-spin fixed-circle evaluator reproduces exact `4x4x4` and `5x4x4` controls and reduces a `5x5` layer to 2,105,872 representatives with a 986,550,880-byte explicit live layout.  No isotropic `5^3` zero or root count is claimed (`proofs/lee_yang_53.md`).
- **Lee--Yang (wave 8):** 2D calibration now PASSES (envelope contains `-1/6`); 3D provably non-identifiable at current sizes (`sigma = 1/2` vs `-1/4` both fit exactly); no exponent reported (`proofs/lee_yang_prg.md`).
- **Scoped no-go results (residual items not covered by the wave-7/8 bullets above):** no tridiagonal pair exists for the tested 3D layer graphs (`proofs/tridiagonal_nogo.md`); the strict cubic-covariant two-weight Kac--Ward slice is inconsistent at `v^8` (the FULL 25-parameter family status is the wave-7/8 bullet above: 31/32 branches empty over Q, positive-real cone empty, branch `11111` open); natural local conserved Pauli charges are absent through tested range, now strengthened to the full-operator support-resolved theorem above; tetrahedron: the constant, six-space edge-rapidity, and ALL arbitrary-`8x8`-auxiliary branches fail exactly (the ungraded theorem above covers parity-breaking `R`; higher auxiliary dimension `>2` and nonidentical `L` factors remain undecided); no nontrivial Pauli flux centraliser exists for the individual term generator set.
- **Claim falsification:** Zhang 2007 and Degang Zhang 2021 fail independent critical/series/symmetry tests; both evaluable claimed free energies have logarithmic rather than external `alpha=0.110(1)` critical behavior (`notes/falsification.md`, `notes/cft_constraints.md`).
- **Wave-11 local commutant (headline, all-size):** for every `a,b` with `ab != 0`, the only compact finite-support operator on infinite `Z^3` commuting with `H_g = a sum_x X_x + b sum_(<xy>) Z_xZ_y` (finite-commutator-derivation sense) is a scalar; the simultaneous `A`/`B` centralizer is a corollary.  Proof: private outside neighbour `w=v+e_1` at a maximal-first-coordinate support site gives `Z_w`-slot `= b[Q,Z_v]`, then `Y_v`-slot `= 2ia Q_Z`, then support induction.  This supersedes the Pauli-string flux theorem in operator class, generator count, and volume.  Explicitly outside scope: extensive translation sums of local densities (the 1D integrable class), quasilocal/infinite-support charges, finite-torus `prod X`, `a=0` or `b=0` (`proofs/local_commutant.md`).
- **Wave-11 ordering-free spectral no-go:** at `t=1/3` the `2x3` layer spectrum is not any full six-mode Gaussian subset-product multiset, with no ordering, inertia window, positivity, or simplicity hypothesis.  A Gaussian forces `deg gcd(C_2,C_2') >= 2016-665 = 1351` on unordered eigenvalue-slot pair products; the monic integral model `R_int = D R`, `D = 3^15 5^4`, makes every modular gcd degree an upper bound, and the exact value is `385` at `p=1000003` and `p=2000003`, hence at least `1631 > 665` distinct pair products.  Controls: open chain `n=6` and a synthetic six-mode Gaussian give exactly `1351`, a four-mode synthetic exactly `55`.  Char-0 equality with `385` is not claimed (`proofs/pair_product_obstruction.md`).
- **Wave-11 series frontier:** the sole HT Euler-ODE `TRUNCATION_UNOBSERVABLE` survivor (order 2, degree 6, 21 columns) is refuted at its first unobservable order: exact residual `499752451221372610239349461368647432702945848808/11` at `v^24`, with `v^26`/`v^28` recorded confirmation-only and full column rank 21 over all 29 equations; eight frozen `e43` refutations plus one `NO_RELATION_AT_BUDGET` row replay as controls.  The HT side of the searched frontier now has no unexplained survivor; this is still finite-prefix and finite-ansatz, not non-D-finiteness (`proofs/series_frontier3.md`).
- **Wave-11 SAW union:** exact inclusion-exclusion of `H` (`75473476338501000000`), `L` (`644233352324156721`), and `C` (`333543290395947705`) with `H cap L = {}`, `|C cap L| = 50032713330104219`, `|C cap H| = 157640944278888` gives `a_35 >= 76401062626946721319` and `mu^36 <= 2941294455272074779839351`.  Only the first-backtrack/Fekete chaining is all-size; the union is a finite `a_35` certificate (`proofs/saw_union.md`).
- **Wave-12 spectral (coupling-generic):** the pair-product invariant is now generic in the coupling.  Because the pair polynomial is monic in `z` over `Z[t]` after the scalar `(2t)^3(1+t^2)^4`, specialization can only raise the derivative-gcd degree, so one evaluated point bounds the generic value, and `rho = Res_z(A/G, B/G) != 0` confines the exceptions to a finite set (crude bound `422472960`).  Hence for every `t>0` outside that finite set the `2x3` layer is not a full six-mode Gaussian spectrum.  The same invariant decides `2x4` at `t=1/3` (gcd degree `9329` vs floor `26335`, so at least `23311 > 6305` distinct pair products); `3x3` walled out at 90 s with no claim (`proofs/pair_product_scale.md`).
- **Wave-12 parity (Theorem GP):** for all `t` outside a finite set (bound `25534083`, including the scalar zeros `0, ±i`), the two physical `P=prod X` sector spectra of `2x3` are not the even/odd halves of any single six-mode Gaussian, in either assignment.  Proved parity class counts `A00=A11=(3^m+1)/2-2^m`, `A01=(3^m-1)/2` give the within-sector floor `195` at `m=6`; the `P=+1` sector has modular gcd degree `177` at `t=1/3` and at five further rational samples.  This is the ordering-free, coupling-generic replacement for the wave-6 inertia-based Theorem S' (`proofs/parity_pair_product.md`).
- **Wave-12 series grid:** 633 admissible cells over algebraic/first-order differential-algebraic/Euler-ODE/Mahler spaces give `322 NO_RELATION_AT_BUDGET` (stored exact nonzero minors), `159 CANDIDATE_REFUTED_BY_HOLDOUT` (first order plus exact residual), and `152 CANDIDATE_SURVIVES`, every one certified truncation-unobservable, so zero unexplained survivors.  The verifier replays all 633 cells, all 152 survivor certificates, and all 18 legacy LT rows from raw coefficients (`proofs/series_grid.md`).
- **Wave-12 extensive charges (the class Theorem LC excluded):** for translation-covariant `Q = sum_x tau_x(q)` with finite density `q`, the derivation vanishes iff the local commutator is a lattice divergence, so the nontrivial quotient dimension is `rank S_R - rank M_R - 1`.  The 1D chain control gives exactly `2R` for `R=0..9` (oriented energy current already at `R=1`), while `Z^2` at `R=1,2` and `Z^3` at `R=1` give exactly `1` — only the Hamiltonian density, with the larger two certified by two-prime ranks meeting the independent `I,h` upper bound `rank S_R - 2`.  A finite-box corollary extends each completed case to all but finitely many nonzero `b/a`.  `Z^2` at `R>=3` (`4^16` densities) and `Z^3` at `R>=2` (`4^27` densities) were never launched and are labelled input-size preflights, not observed walls; quasilocal and non-translation-covariant charges stay unclassified (`proofs/extensive_charges.md`).
- **Wave-13 all-size spectral theorem (Theorem F, headline):** for EVERY finite graph with a vertex of degree `>= 3` and EVERY `n >= 4`, at inhomogeneous parameters off an explicit hypersurface of degree `<= 4(2n+|E|)S(S-1)` with `S = C(2^n,2)`, the layer spectrum has at least `3^n - 2^n + 48*3^(n-4)` distinct unordered-slot pair products, hence is not a full `n`-mode Gaussian subset-product multiset; the excess is a constant fraction `16/27` of `3^n`.  Mechanism: Lemma 1 (Gaussian ceiling `3^n-2^n`), Lemma 5 (exact tensor identity `r = (3^m-2^m) r_all(core) + 2^m r(core)`), the claw seed `r_all >= 129 > 81`, and monic specialization (which transfers gcd upper bounds special -> generic only).  Corollary F1: the exceptional set has Lebesgue measure zero, so almost every random-field/random-bond branching layer is non-Gaussian at every size.  SCOPE: the physical isotropic curve is NOT covered at any `n`, and section 9 proves that is an obstruction of the method — the isotropic claw is not certified above the ceiling and equal-field localization is provably linear in `m` (`proofs/allsize_gaussian.md`).
- **Wave-13 isotropic genericity (Theorem G) and Theorem U:** on the isotropic curve the scaled operator is integral with entry degree `<= 2n+2|E|`, so the pair polynomial is monic over `Z[t]` and one named point plus one resultant give the conclusion for all but finitely many `t`; applied to `T u P_m` (`m = 0..4`, uniform-field no-go at `n = 5..9` with zero anisotropy) and to the `2x3`, `2x4`, `3x3` layers, with the `2x3` bound reproducing the independent `422472960`.  Withdrawn on review: exactness of the counts `3893`/`35597` and any certified multiplicative independence — substituting modular lower bounds into an upper bound proves nothing, so those survive only as the doubly conditional Theorem F2 (`proofs/allsize_gaussian.md`).
- **Wave-13 `3x3` layer (first degree-4 vertex):** modular pair-product gcd degree `59641` at `p = 1000003` and `2000003` against the Gaussian floor `111645`, so at least `71175 > 19171` distinct pair products; two independent implementations agree with matching pair-polynomial digests, and the `n = 9` open chain hits its floor `111645` exactly.  This closes the 90-second wall recorded in `proofs/pair_product_scale.md` §6 (`proofs/pair_product_3x3.md`).
- **Wave-13 LT series:** exact low-temperature reduced free energy extended from `x^32` to `x^52` (ten new coefficients, `x^34 = 666750` .. `x^52 = -46131904167/2`) with per-box CRT uniqueness certificates and a hashed post-computation Guttmann-Enting witness; `x^34` is the first order that can test the 18 legacy LT truncation-unobservable rows, of which **11 are now refuted** with exact residuals while 7 need `x^54` (blocked by the 22-spin cross-section guard) (`proofs/lt_x34.md`).
- **Wave-13 extensive charges beyond the wall:** support-size stratification is exact (a bond commutator grows the support by at most one site), so named classes are finite work: the `3x3x3` box class with `<= 4` non-identity sites per word has quotient exactly `1` (rank `822334 -> 822332` at two primes, meeting the analytic `I,h` bound), and `Z^2` radius 3 gives `1` up to support size `5`.  The FULL `Z^3` radius-2 box is NOT decided: `s <= 5` is an observed 5 GiB wall, and the order-48 point group still leaves `>= 3.75e14` invariant coefficients (`proofs/extensive_3d_r2.md`).
- **Wave-13 ladder reduction:** `A_L = 2L - 2N` and `B_L = D + F + Ddag` in the `X` eigenbasis with `ad_A` eigenvalues `-4, 0, +4` and all three at bracket depth 3; `dim g_L >= dim(g_L psi)` with a stabiliser-module argument moves the problem from `16^L` Pauli to `4^L` state coordinates, and `g_L psi` lies in an invariant space of dimension `2^(2L-3) + 3*2^(L-2)`.  CONDITIONAL: if that space is saturated then `dim g_L >= 2^(2L-3) >= 2^L` for all `L >= 3`, quadratically stronger than the conjecture.  Saturation is verified only at `L = 3`; state-space certificates give `dim g_L >= 2^L` for `L = 3..8`; `L = 9` reached rank `420/512` and is an explicit non-certificate (`proofs/ladder_alll_proof.md`).
- **Wave-13 interlayer `c8`:** exact through `v^12`, open-slab FLM agreeing independently through `v^12` and spin-DOS through `v^6`; the all-order identities emerge from an exhaustive evaluation of all 3,432 normalized gap multisets rather than being imposed (`proofs/interlayer_c8.md`).
- **Wave-13 local-charge literature scope:** the general absence of local conserved quantities in `d >= 2` for this model class is due to Chiba (arXiv:2412.18903, PRB 111 195130 (2025)); our finite-shape reconstruction (T1/T2/M, single-shape kernels, residual kernel dimensions 2/6/28/184/120) is tagged `[CONJECTURE]` plus `[COMPUTATION]`, not a rival proof.  An earlier draft of `proofs/allr_charges.md` overstated four items as theorems and was corrected (`proofs/allr_charges.md`).

- **Wave-14 2x5 pair-product (headline):** at `t=1/3` the open `2x5` layer transfer operator (dim 1024, 13 bonds) is not a full ten-mode Gaussian subset-product spectrum (Theorem P10): under the actual-LCM integral model `D = 3^26*5^7` the modular pair-product gcd degree is `124921` at both `p=1000003` and `p=2000003`, below the Gaussian floor `465751 = C(1024,2) - (3^10-2^10)`, hence at least `398855` distinct pair products against the Gaussian maximum `58025`; the `n=10` open-chain control hits its floor `465751` exactly and the `2x3`/`3x3` regressions (385, 59641) replay with matching digests; Theorem G10 upgrades to every positive real `t` outside a finite set of size `<= 51576065587200`.  Char-0 gcd equality with `124921` is not claimed (`proofs/pair_product_2x5.md`).
- **Wave-14 ladder (saturation refuted):** `SectorSaturationCyclicity` — the single missing step of the wave-13 conditional theorem — is refuted by certified falsity at `L=4,5,6,7`: `dim(g_L psi) = K_L - W_L` with certified 4B-invariant annihilator blocks `W_L = 0,2,10,66,364` and exact cyclic dims `14/14, 42/44, 142/152, 494/560, 1780/2144` at `L=3..7` (`dim K_L = (2^(2L-1)+3*2^L)/4` canonical by Burnside, computed two routes).  The certified bound `dim g_L >= K_L - W_L` at `L=3..7` subsumes the `2^L` certificates there; three independent probes agree (tensor `e131`, duality pairing `e132`, state-space `e130`); the one-rung cut is certified NOT `K`-preserving by an explicit rho-breaking orbit witness; nothing is claimed at `L >= 8`, `dim g_L` for `L >= 4` is unchanged, and the `W_L/K_L` asymptotics (0, 4.5%, 6.6%, 11.8%, 17.0%) are unresolved (`proofs/sector_saturation_tensor.md`, `proofs/sector_saturation_pairing.md`).  This supersedes the conditional statement in the wave-13 ladder bullet above.
- **Wave-14 extensive charges (`s<=5` decided):** the `Z^3` radius-2 support-size `s<=5` class `C(3,2,5)` is decided: `21,121,156` columns, `rank S = 14,757,412`, and column-oriented exact `F_p` elimination with compact CSR pivot storage gives `rank M = 14,757,410` at both `p=2147483647` and `p=2147483629` (1446.7 s / 1460.2 s process, peak RSS 2292.7 MiB), meeting the analytic `I,h` upper bound `rank S - 2`, so the quotient is exactly `1` — only the Hamiltonian density survives modulo shifts and identity.  The certified `s<=4` regression (`822332` at both primes) is reproduced first; the FULL radius-2 box and `s >= 6` remain untouched (`proofs/extensive_3d_s5.md`).
- **Wave-14 Peierls limitation (Theorem P):** the plain Peierls head+tail contour certificate cannot beat the incumbent `I_3/2`: with the exact cell-class census through `n<=10` (9,638,143 classes flood-filled, zero cavity-bearing; the smallest cavity needs exactly 11 cells, attained by an 11-cell tree witness) the rooted contour counts extend to `N(26)=20448` and `N(28)=33240`, and the certified partial sum at `x*=e^{-I_3}` already exceeds `1/2` (`V(x_-) = 0.5883195276...`, margin `0.08832`); `K* = 0.2558326 > I_3/2` is 3-way-bisection-certified (gap `0.003102`); the head-only crossing's feigned gain `0.01804` is annihilated by the proved tail mass `1.65288 > 1/2`; the Lagrange/Lagrange-tail converges only above `(1/2)log lambda = 1.6755`, a certified `1.42` above the incumbent, so no finite exact head can close the gap; re-certified grid endpoint `1.695292` (below wave-6's `1.696898`).  `K*` is a property of the certificate, not of `K_c`; the incumbent upper endpoint is unchanged (`proofs/upper_beyond.md`).
- **Wave-14 Kac--Ward `11111` (anchor components cleaned):** over the anchor thin-torus section `(1,1,1,1,1,4,4)` the unit chart has exactly three nonsingular `F_5` solutions (`det J` in `{2,4,1}`; the full `F_5^6` grid has `1553` solutions, the fourth nonsingular point being chart-degenerate), each lifting by Hensel to a genuine `Q_5` point of the 42-equation construction variety (all 42 construction residues `0 mod 625`), and all three lifts fail the independent `3x3x3` holdout at orders `8/10/12 mod 625` (residues `[250,55,266]`, `[456,350,132]`, `[430,480,286]`; the anchor `k=8` residue `250` reproduces wave-10).  Two independent routes (symbolic walk-monomial expansion and `108x108` matrix-power traces) agree at all five orders, and the orbit-union lemma collapses the diagonal torus action to one representative per orbit.  Branch `11111` stays `[UNRESOLVED]` over `Q`; the other `4^7-1` thin-torus sections are unanalysed (`proofs/kw_components.md`).
- **Wave-14 SAW sieve union:** `a_35 >= 1972465461070186257835` from the 404-member two/three-block height-schedule family closed by exact inclusion-exclusion — overlap correction `5030204334415261659120` against the naive disjoint-union volume `7002669795485447916955` (2,780 active sieve terms, 8,344 covered chains, 2,167 intersecting pairs; full member counts rest on the wave-13 820 certificates).  This set the wave-14 lower endpoint, superseded by the wave-20 memory-12 bound; the clean-room verifier independently re-derives `T(13)=142016661`, the `l<=11` band, the per-row sieve identities and the 40-dp floor (`proofs/saw_union4.md`).
- **Wave-14 ledger state (`H380..H390`), canonical mapping copied from ledger row H386:** "The five earlier appended rows H380-H385 are: H380 = 2x5 spectral no-go (canonical); H381 = SectorSaturationCyclicity refutation (canonical); H382 = Z^3 s<=5 class (canonical); H383 = DUPLICATE of H380 (2x5, mistaken re-append); H384 = CORRECTION to H381 quantifier; H385 = DUPLICATE of H382 (s<=5, mistaken re-append).  The H383/H385 rows are superseded by their canonical H380/H382 and by THIS row; they are NOT additional certification."  The remaining wave-14 rows are `H386` = the mapping correction itself, `H387` = the Peierls-at-incumbent limitation theorem, `H388` = the KW `11111` anchor-section component theorem (its column alignment is recorded by row `H402` — the artifact list sits in the `next_action` slot; append-only, to be marked by a future row), `H389` = the height-schedule sieve union, `H390` = CORRECTION scoping `H388` to the unit chart.  The ledger holds `362` data rows, all unique IDs, highest `H390`; **next free ID `H391`**.

Exact coefficient arrays, dimensions, residuals, methods, and scripts are recorded without rounding loss in `checkpoints/verified_results.json`.

## Essential corrections that must survive future summaries

1. **Dolan--Grady zero set:** “DG iff 2-regular” is false without a no-isolated-vertices assumption.  The exact condition is `deg(v) in {0,2}` for every vertex.  `C3` plus an isolated vertex is the decisive counterexample (`results/audit/audit_replication.json`).
2. **Claws and algebra size:** an induced claw is the exact local term-wise Majorana/quartic-defect obstruction, but is **not equivalent** to a large two-sum algebra.  `K_(1,3)`, `K_(1,4)`, and `K_4` have dimensions `24`, `33`, and `15` (`results/minimal_obstruction.json`).
3. **Flux theorem scope:** the proved centraliser is for Pauli strings commuting with every individual `X_i` and `Z_iZ_j`.  It is not the full commutant of the sums `A,B`; translation on `C3` is a counterexample (`proofs/algebraic_obstruction.md`).
4. **Tridiagonal certificates:** the no-go is exact, but a universal one-row witness is false.  Some geometries require the stored four-row certificate (`proofs/tridiagonal_nogo.md`).
5. **Kac--Ward scope:** the strict two-weight slice being ruled out never decided the full 30-weight family.  Current status: positive-real cone empty, 31/32 rational branches empty, branch `11111` open over `Q`/`C`.  Its nonsingular `Q_5` component proves the construction ideal is proper, so no Nullstellensatz certificate can establish branch emptiness at any degree; the recorded local component nevertheless fails an independent holdout (`proofs/kw_cubic.md`).
6. **Watson normalisation:** the correct gamma denominator is `32*pi^3`; the problem-stated `4*pi^3` expression is eight times too large (`proofs/kc_bounds.md`).
7a. **Cross-artifact digests must pin content, not bytes:** `test_levi_images_2x4.py` broke
   when `test_oa_quotient.py` regenerated `oa_quotient.json` with fresh `meta` timestamps.
   Cross-artifact references now digest the `data` section with sorted keys
   (`classification_artifact_content_sha256`); never pin whole-file bytes of a mutable artifact.
7. **`plus_boundary` length-one edge case:** an independent audit found and repaired a missing ghost bond when the transfer length was `c=1`; `tests/test_plus_boundary.py` now locks the corrected behavior.
8. **Ledger schema correction:** H360 was appended with a separate claim-tag field against the 11-column legacy header, shifting its semantics.  Append-only row H361 is the authoritative field-aligned statement; H360 remains historical and must not be deleted.
9. **Ledger field-convention drift, wave 17 only:** rows `H417`-`H427` are 11 fields wide but do not follow the legacy header semantics; they use `id, date, statement, topic, owner, status, artifact, artifact, detail, status-word, next-action`.  A parser that reads `mathematical_motivation` from those rows gets a topic word, not a sentence.  Width anomalies elsewhere (`H388` 12 fields, `H389`/`H390` 6, and `H393`-`H397`/`H399`-`H401` 10) are historical and recorded by append-only row `H402`.  Full inventory and the current ID allocation are in `checkpoints/next_actions.md` §1.

## Open mathematical state

No thermodynamic-limit exact 3D free energy, magnetisation, correlation function, spectrum, exact critical coupling, or exact exponent has been produced.  The following are specifically open rather than silently ruled out:

- an all-`L` proof or counterexample for **exponential** two-generator ladder/grid growth; the
  fixed-basis quadratic no-go is all-`L>=3`, and two-slice observability is exact at `L=3..8`,
  but the two-slice exact-sequence lemma and exact `W_L` rank law remain open;
- the every-coupling/all-size isotropic spectral theorem on **bipartite** branching grids; wave 21
  gives an explicit log-moment interval on open `2x2`/`2x3` and a trace-seven obstruction at
  `t=1/3` on open `2x3`, while the all-size non-bipartite trace invariant vanishes here;
- simultaneous nonlinear deformations of both `A` and `B` beyond the invariant F1--F4 supports;
- Kac--Ward branch `11111` over `Q`/`C` (signed/complex weights; sole surviving scalar branch);
- nonidentical/spectral/IRF/dynamical tetrahedron mechanisms and enlarged-Hilbert-space fluxes;
- memory-14 finite-walk continuation: the compact `e229`/`H590`-`H595` crosswalk is landed and
  verified through `k=12`, but memory 14 was explicitly not launched; an endpoint attempt must
  still construct `A_14`, verify its mass, replay an integer Collatz rotation componentwise, and
  apply directed `atanh` rounding; HT `v^30`, LT `x^60/x^62`, and `Z^3` radius-2 extensive
  charges at support size 6 also remain open;
- phase-sensitive input capable of improving the infrared upper endpoint; the audited `PS4`
  convolution-block class is strictly stronger than `MR4` but has zero uniform Green floor, so
  full cross-channel fourth moments, higher localizers, DLR equations, or sourced/multi-edge
  current identities are still required;
- scalar evaluation/bounds for the forced interlayer `W8/W10` sums and controlled Lee--Yang
  thermodynamic edge scaling.

The finite F3 deformation in `proofs/dg_deformation_nogo.md` cancels DG2 on the `3x3` torus but fails DG1; it is a positive cancellation result, not an integrable solution.

## Fresh full test run

Acceptance history, most recent first. Each run is `bash tests/run_all.sh` from the repository root
as a supervised process; the runner expands `tests/test_*.py` once at start, so a run covers the
scripts present when it was launched.

| wave | scripts | runner line | notes |
|---|---:|---|---|
| 21 | 153 present | no full-suite launch; **5/5 landed targeted standalone verifiers passed** | Producer artifacts passed `17/17`, `5/5`, `16/16`, `29/29`, and `33/33`; lead observed PASS from the bipartite log-moment, trace-resultant, W-law rank-mechanism (21 checks), phase-sensitive upper, and compact finite-memory SAW verifiers. The fifth verifier proves the packed `k<=12` crosswalk only; memory 14 was explicitly not launched and no endpoint result is claimed. |
| 20 | 148 present | no full-suite launch; **13/13 new standalone verifiers passed** | Includes the independent 41,424-state memory-12 rebuild (18 checks); one initial susceptibility verifier failure exposed a false “strictly decreasing” wording at the exact `5,5` plateau, then passed after producer/proof/artifact correction. |
| 17 | 123 | `FINAL: 123 total, 123 passed, 0 failed` | `fullsuite17b`, 8 h 10 m under heavy host contention (load ~40, hundreds of foreign jobs). Includes the standalone 38-check `H425` `2x6` verifier and the independent `x^56` fresh-prime audit. The earlier `fullsuite17` was snapshot-launched before `tests/test_gaussian_2x6.py` existed and was restarted so the gate would cover it. |
| 16 | — | no separate full-suite gate | wave-16 deliverables were each run through their own standalone verifier before their ledger rows were appended, and are covered by the wave-17 `fullsuite17b` expansion. |
| 15 | 115 | `FINAL: 115 total, 115 passed, 0 failed` | `fullsuite15b`, run pre-wave-16 on the integrated tree. |
| 14 | 110 | `FINAL: 110 total, 109 passed, 1 failed` | fullsuite14: 110 total, 109 passed, 1 failed (test_extensive_3d_r2.py, dropped artifact-load line, lead regression); standalone rerun after one-line restore: PASS, 555 s. Full suite re-run deferred to wave-15 close-out. |
| 13 | 103 | `FINAL: 103 total, 102 passed, 1 failed` | 14 h 25 m; the single failure is `tests/test_ladder_all_l.py` firing its wall-clock `TimeoutError: NON-DECISIVE: exceeded 3600s` under heavy unrelated machine load; it passes in isolated rerun (1565 s), and the gate has since been converted to `time.process_time()` so it can no longer fire on elapsed time |
| 12 | 95 | `FINAL: 95 total, 94 passed, 1 failed` | the same wall-clock gate in `tests/test_ladder_all_l.py`, same load; passes in isolation |
| 11 | 90 | `FINAL: 90 total, 90 passed, 0 failed` | 2 h 49 min |
| 10 | 86 | `FINAL: 86 total, 86 passed, 0 failed` | 3 h 56 min, run on a frozen tree while wave 11 worked in isolated worktrees |
| 9 | 75 | `FINAL: 75 total, 75 passed, 0 failed` | 81 min; the first wave-9 run scored 74/75 and exposed the artifact-coupling bug fixed in correction 7a above |

Independently of the suite, **every** wave-11 through wave-21 landed deliverable was run through its
own standalone verifier by the lead inside this checkout before its ledger rows were appended.
Wave 20 adds 13/13 and wave 21 adds 5/5; the per-test wall times for wave 13 are tabulated at the
end of `research_log.md`.
The audit-replication script still deliberately prints `FAIL` for two historical over-general
claims while exiting 0.

**Reading a suite failure.** Several verifiers carry their own wall budgets and declare
`NON-DECISIVE` when they expire. Under heavy machine load such a line is a resource artefact, not a
mathematical regression: recheck that single test in isolation before treating it as one.

## Reproduction hygiene

- Run from the repository root; scripts assume that location.
- Use `.venv/bin/python`, not the system interpreter.
- Experiments overwrite their corresponding JSON outputs and refresh timestamps.
- Do not use `0.221654626` to fit or select a derivation; it is comparison-only.
- Treat ranks over the two recorded primes as rigorous lower bounds over `Q`; two-prime agreement is a cross-check, not a universal rational-rank proof.
- Preserve result scopes exactly: finite, ansatz-bounded, numerical, external, or theorem.
