# HT `v^24` and LT `x^54` series walls: exact feasibility audit (READ-ONLY)

Wave-15 feasibility note, agent `HTv24Feasibility`.  READ-ONLY by contract: no series
computation was launched; the only executed code was a print-only exact-integer wrapper
(under `nice -n 19`, single-thread libraries).  Every integer below is either (i) quoted from
a cited artifact key or (ii) re-derived in that wrapper from the engine's own stored
metadata/formulas.  Claim tags follow the repository convention
([THEOREM]/[COMPUTATION]/[CONJECTURE]/[EXTERNAL]/[UNRESOLVED]).

Referenced engine sources are in the main checkout:
`src/ising/transfer_matrix/__init__.py` (dense int64 engine),
`src/ising/transfer_matrix/crt.py` (31-bit-prime CRT engine),
`src/ising/series/__init__.py` (Enting finite-lattice assembly).

## 0. Engine facts that every wall reduces to

**Storage model [COMPUTATION, from source].**  Both dense engines propagate one array of
shape `(1 << ns) x (deg+1)` int64, `ns` = cross-section sites:

* allocation: `__init__.py:176` (`vec = np.zeros((1 << ns, deg + 1), dtype=np.int64)`);
  same shape in the CRT engine, `crt.py:248`;
* `_shift_rows` (`__init__.py:92`) allocates a full second array `out` (`:95`);
* `_apply_interlayer` (`__init__.py:106`) / `_apply_interlayer_mod` (`crt.py:168`) do one
  full-buffer copy `out = v.copy()` per spin bit (`__init__.py:120`);
* so at most TWO big arrays are mechanically live at any instant; the per-layer scratch
  vectors (`_broken_in_layer`, `_weighted_downcount`, `__init__.py:76-89`) are `2^ns` int64
  each: `3 * 2^25 * 8 = 805,306,368` bytes at `ns=25`, subleading.
* The repository's own peak model is the conservative closed formula
  `bytes = 3 * 2^(a*b) * (retained+1) * 8` in `estimate_box_peak_bytes`
  (`crt.py:403-425`, docstring line 410: "Conservative peak for three simultaneous int64
  propagation arrays").  It calls `_shape_data` (`crt.py:412`) and therefore ITSELF refuses
  any cross-section above the guard.

**Guards [COMPUTATION, from source].**

* int64 engine: `if ns > 22: raise` (`__init__.py:149-150`) and
  `if ntot > 62: raise` (`__init__.py:151-153`).  The second guard is an exactness guard:
  coefficients can reach `2^N` (module docstring, `__init__.py:8-17`), so a 125-site box
  genuinely overflows int64; the int64 engine is structurally out for the boxes below.
* CRT engine: `_MAX_CROSS_SECTION = 22` (`crt.py:47`), enforced in `_shape_data`
  (`crt.py:151-155`).  Its rationale is OVERFLOW HEADROOM of the modular butterfly, not
  memory: docstring `crt.py:11-15` — one inter-layer butterfly multiplies a residue below a
  `2^31` prime by at most `2^ns`, so a complete `ns`-bit pass stays below `2^(31+22) = 2^53`
  before the single modular reduction (`crt.py:183`, after the whole pass).
  **[COMPUTATION]** The same bound at `ns=25` is `2^56 < 2^63 - 1`: int64 headroom factor
  `2^7 = 128`.  Raising the constant to 25 cannot change any output value; the CRT engine
  must then be used because `2^125 > 2^63 - 1`.
* The same pair (22, 62) is hardcoded in the HT assembly skip
  `unreachable_but_proven_zero` (`series/__init__.py:304-310`).

**Truncation asymmetry between the two series paths [COMPUTATION, from source].**

* LT: `_plus_broken_polynomial` (`series/__init__.py:197-216`) feeds `log_series`, which
  consumes only the low-degree prefix; the x^52 producer therefore routes N>62 boxes through
  `crt.box_broken_bond_poly(..., max_degree=order)` (`experiments/e124_lt_x34.py:180-188`).
  Truncation is coefficient-exact for this use.
* HT: `_free_broken_polynomial` (`series/__init__.py:191-193`) has NO degree parameter, and
  `broken_bond_to_even_subgraph` asserts `sum(c) == 2^N` (`series/__init__.py:130-132`):
  the conversion `P(v) = 2^{-N} sum_q c[q] (1-v)^q (1+v)^{n_b - q}` (docstring
  `series/__init__.py:105-116`) mixes ALL broken-bond degrees into every low `P` coefficient.
  Hence the HT dense path is intrinsically FULL-degree: 301 columns for `(5,5,5)`, never 25.

## 1. HT `v^24`

### 1.1 Why cross-section 25 exceeds the 22-spin guard

At order 24, `_ht_shapes(3, 24, 0)` uses `extent_budget = 24 // 2 = 12`
(`series/__init__.py:231-238`), so canonical `(5,5,5)` (coordinate-span sum `3*4 = 12`)
enters the inventory; its minimal degree is `2*12 = 24 <= 24`, so the proven-zero skip
(`series/__init__.py:304-310`) does NOT fire; `_ht_box_log((5,5,5), 24)`
(`series/__init__.py:218-223`) calls the dense engine, which raises
`cross-section (5, 5) has 25 sites: too large` (`__init__.py:149`).  A canonical sorted
shape puts the largest side last (`series/__init__.py` `_canonical_shape`), so the
cross-section is `5*5 = 25 > 22`.  **[COMPUTATION]** `v^22` was the last clean plain order:
extent budget 11, largest inventory box `(4,5,5)`, cross-section 20 (consistent with every
non-walled transfer record at `v^28` having `cross_section_sites <= 21`,
`results/series/ht_v28.json` `transfer_records_large_cross_section`).

### 1.2 The 242 GB figure, re-derived from the stored per-box sizes

Stored block `results/series/ht_v24.json` `data.old_engine`:
`cross_section_sites 25`, `spin_states 33554432` (`= 2^25`), `full_broken_bond_degree 300`,
`simultaneous_int64_arrays 3`, `projected_peak_bytes 242397216768`,
`crt_prime_count_for_full_spin_polynomial 5`.  Recomputation: the free `(5,5,5)` box has
`2*(5*4) = 40` in-layer bonds per layer (`layer_bonds`, `__init__.py:40-73`) and
`40*5 + 25*4 = 300` physical bonds (`_box_residues` formula, `crt.py:222`), hence 301
columns.  [COMPUTATION]

    3 * 2^25 * 301 * 8 = 242,397,216,768 bytes = 242.40 GB = 225.75 GiB

The figure is CONFIRMED, not refuted.  It is the engine's THREE-array conservative model
(section 0); the mechanical dense lower peak is two arrays:
`2 * 2^25 * 301 * 8 = 161,598,144,512` bytes.
Cross-confirmation that this closed formula is exactly what the series program quotes: the
same model reproduces ALL six stored `v^28` walled-class projections with zero mismatch
(`results/series/ht_v28.json` keys `walled`, `walled_cross_sections`,
`walled_projected_bytes`; `proofs/ht_v28.md`):

    (4,6,6) ->  140,525,961,216     (4,6,7) ->  165,490,458,624
    (5,5,5) ->  242,397,216,768     (5,5,6) ->  294,742,130,688
    (5,5,7) ->  347,087,044,608     (5,6,6) -> 11,467,562,680,320

### 1.3 Reduction (a): even-odd bipartite splitting of the layer

The 5x5 layer grid is bipartite with classes of 13 (even) and 12 (odd) sites
[COMPUTATION].  Three exact statements:

* **Persistent state: shrink factor exactly 1.**  The plane-cost shift multiplies the
  new-layer row `s'` by `x^{E(s')}`, where `E(s')` sums the 40 in-layer bonds and EVERY bond
  joins an even site to an odd site (checkerboard grid, `layer_bonds`
  `__init__.py:57-72`).  The operator factorization
  `T = D_plane * prod_even (I + x X_i) * prod_odd (I + x X_j)` is exact, but a dense
  application still materializes a full `2^25 x (deg+1)` intermediate: the odd-sublattice
  butterfly mixes all 12 odd bits jointly for every even setting, and the even pass then
  mixes all 13 even bits for every odd setting — there is no streaming schedule with
  sub-`2^25` live dense storage.  (Structural statement about THIS dense array engine, not
  a lower bound over all exact algorithms; the frontier engine of section 1.5 needs far less
  precisely because it does not track general states.)
* **Transient workspace: factor exactly 2 on the copy.**  A butterfly on bit `i` pairs rows
  differing only in bit `i`, so the `out = v.copy()` transient (`__init__.py:120`) can be a
  half-array slice.  Estimator-model peak: 3 -> 2.5 array units:
  `242,397,216,768 -> 201,997,680,640` bytes = 188.13 GiB (full degree), respectively
  `44,291,850,240 -> 36,909,875,200` bytes = 34.38 GiB (LT-style 55 columns).
  NOT below 8 GiB.
* **Adjacent exact factor 2: global spin-flip sector.**  For the FREE box, `s -> -s`
  commutes with every operator (`E(-s) = E(s)`, vertical-bond disagreement unchanged) and
  the all-ones start vector and the all-sum readout are flip-symmetric, so the computation
  may live in the `+1` sector: factor exactly 2:
  `121,198,608,384` bytes = 112.88 GiB full-degree.  Even a hypothetical truncated HT path
  would sit at `22,145,925,120` bytes = 20.63 GiB — still above 8 GiB (and truncation is
  unavailable on the HT path, section 0).  For the LT plus-boundary box even this factor is
  false: the ghost term sums `w_i * [spin i is DOWN]` (`_weighted_downcount`,
  `__init__.py:76-89`) and the frozen `+1` exterior breaks the symmetry.

**[COMPUTATION] Verdict (a): the layer state-space shrink factor is exactly 1 (bipartite
streaming) or 2 (flip sector, free boxes only); nothing reaches 8 GiB.**

### 1.4 Reduction (b): slab (row-cyclic) symmetry

* **For the canonical open FLM box this reduction is not a symmetry at all (shrink factor
  NOT APPLICABLE):**  making a section direction periodic adds wrap bonds (`layer_bonds`,
  `__init__.py:66-72`, the `elif pa` / `elif pb` branches): in-layer bonds go 40 -> 45 (one
  direction) or 40 -> 50 (both), and `W(5,5,5)` is defined by the OPEN box.  A cyclic
  section computes a different object; quoting a factor for it would be a category error.
* **Counterfactual orbit arithmetic anyway [COMPUTATION].**  Row rotations form `C_5`;
  Burnside gives orbit count `(2^25 + 4*2^5)/5 = 6,710,912` exactly; shrink factor
  `33,554,432 / 6,710,912 = 4.99998`.  Estimator-model peaks:
  `3 * 6,710,912 * 301 * 8 = 48,479,628,288` bytes = 45.15 GiB (full degree), and
  `3 * 6,710,912 * 55 * 8 = 8,858,403,840` bytes = 8.25003 GiB at LT-style truncation —
  MISSING the 8 GiB budget (8,589,934,592 B, `experiments/e124_lt_x34.py:51`
  `RSS_BUDGET_BYTES`) by 268,469,248 bytes even with the wrong object.
* **The e26 slab route (both directions periodic, torus section):**  `C_5 x C_5` orbits
  `= (2^25 + 24*2^5)/25 = 1,342,208` exactly; factor 24.99943; full-degree
  `3 * 1,342,208 * 301 * 8 = 9,696,110,592` bytes = 9.03 GiB — STILL over budget.  Its
  truncated 1.65 GiB variant is doubly disqualified: wrong object, AND the route's own error
  analysis caps its validity far below `x^54` — wrapping contamination of an L-periodic
  section starts at surface `4L` (`x^16` at L=4, hence `x^20` at L=5) and omitted
  z-extent-`h` clusters start at `x^{4(c+1)+2}`
  (`experiments/e26_lt_series_independent.py` docstring lines 20-25); reaching `x^54`
  error-free needs `c >= 12` while `4L = 20 << 54` kills the 5x5 section after `x^19`.

**[COMPUTATION] Verdict (b): does not drop any legal computation below 8 GiB; and for the
boxes the FLM actually needs it is not a symmetry of the weight.**

### 1.5 What actually broke the wall (recorded; orientation only)

The dense projections above are historical counterfactuals; the working route was an
algorithm change, the cut-capped parity frontier:

* `v^24` (`experiments/e68_ht_v24.py`): frontier peak 1,727,071 states; isolated production
  maximum RSS 352,305,152 bytes; whole-coefficient run peak RSS 7,305,396,224 bytes
  (`results/series/ht_v24.json` keys `frontier.peak_states`,
  `production_resource_measurement.maximum_resident_set_bytes`,
  `frontier.artifact_generation_peak_rss_bytes`; `proofs/ht_v24.md`).
* `v^26` (`experiments/e86_ht_v26.py`, one cut capped at four): worker peak states
  2,294,347 / 3,420,675 / 3,644,255 / 7,689,913, worker RSS <= 1.70 GB
  (`proofs/ht_v26.md` Resource telemetry).
* `v^28` (`experiments/e94_ht_v28.py`, general cap vectors): 27 profiles, max worker
  29,378,975 states / 7,280,050,176 bytes RSS; producer peak 22,577,053,696 bytes; elapsed
  10,584.737070 s (`proofs/ht_v28.md` Resource telemetry).

## 2. LT `x^54`

### 2.1 Which guard blocks it

Stored preflight (`results/series/lt_x34.json` `data.next_order_preflight`,
[UNRESOLVED]-tagged): `next_order 54`, `required_budget 15` (`= (54+6)//4`),
`required_canonical_shape [5,5,5]`, `cross_section_sites 25`,
`transfer_cross_section_guard 22`, `launched false`,
`status UNLAUNCHED_INPUT_SIZE_PREFLIGHT_WALL` (`experiments/e124_lt_x34.py:138-165`).
The blocking mechanism is the CRT engine's `_shape_data` raise (`crt.py:151-155`); the
order bound is the proved surface inequality `4*(a+b+c)-6 <= order`
(`proofs/lt_x34.md` section 1, proved in `notes/flm_derivation.md`), so `x^54` needs
`a+b+c <= 15`.  Inventory [COMPUTATION]: 455 ordered boxes `= C(15,3)` (vs `364 = C(14,3)`
at `x^52`, `proofs/lt_x34.md` section 1) and `83 + 19 = 102` canonical classes; the 19 new
side-sum-15 classes are tabulated in 2.2 — exactly ONE exceeds 22, namely `(5,5,5)` with
cross-section 25; all others are <= 20 and would run today.

### 2.2 What the per-box CRT uniqueness certificates need at that order

Engine rule: `primes_for_sites` (`crt.py:119-134`) starts at `sites // 30` descending
31-bit primes (`_descending_primes`, `crt.py:106-116`; deterministic 32-bit Miller-Rabin
`_is_prime_32`, `crt.py:73-102`) and increments until `product > 2^N`.  The recorded
artifact certificate is the strengthened TWO-SIDED condition `product > 2*2^N`
(`experiments/e124_lt_x34.py:187-188`; rows `uniqueness_certificate_passed`;
`proofs/lt_x34.md` section 2, [THEOREM] `0 <= [x^q] Xi_A <= 2^N`).  Prime-product bit
lengths are 31/62/93/124/155 for 1-5 primes [COMPUTATION].  The 19 side-sum-15 classes:

| class | N | cross | engine | primes | product bits | > 2^N | > 2*2^N |
|---|---|---|---|---|---|---|---|
| (1,1,13) | 13 | 1 | int64 | - | - | - | - |
| (1,2,12) | 24 | 2 | int64 | - | - | - | - |
| (1,3,11) | 33 | 3 | int64 | - | - | - | - |
| (1,4,10) | 40 | 4 | int64 | - | - | - | - |
| (1,5,9) | 45 | 5 | int64 | - | - | - | - |
| (1,6,8) | 48 | 6 | int64 | - | - | - | - |
| (1,7,7) | 49 | 7 | int64 | - | - | - | - |
| (2,2,11) | 44 | 4 | int64 | - | - | - | - |
| (2,3,10) | 60 | 6 | int64 | - | - | - | - |
| (2,4,9) | 72 | 8 | CRT | 3 | 93 | yes | yes |
| (2,5,8) | 80 | 10 | CRT | 3 | 93 | yes | yes |
| (2,6,7) | 84 | 12 | CRT | 3 | 93 | yes | yes |
| (3,3,9) | 81 | 9 | CRT | 3 | 93 | yes | yes |
| (3,4,8) | 96 | 12 | CRT | 4 | 124 | yes | yes |
| (3,5,7) | 105 | 15 | CRT | 4 | 124 | yes | yes |
| (3,6,6) | 108 | 18 | CRT | 4 | 124 | yes | yes |
| (4,4,7) | 112 | 16 | CRT | 4 | 124 | yes | yes |
| (4,5,6) | 120 | 20 | CRT | 4 | 124 | yes | yes |
| (5,5,5) | 125 | 25 | CRT | 5 | 155 | yes | yes |

So 9 classes use the int64 engine and 10 need CRT with 3,3,3,3,4,4,4,4,4,5 primes.
For `(5,5,5)`, `N = 125`: 4 primes give 124 bits `< 125`, so the rule lands on 5 primes —
the first five descending 31-bit primes
`2147483647, 2147483629, 2147483587, 2147483579, 2147483563`, product
`45671921168693645933699105804560590380377589537` (155 bits; same five primes recorded in
`results/series/ht_v24.json` `crt` block).  `155 > 126 = 1 + 125`: both the engine rule and
the strengthened two-sided certificate pass.  Time multiplies by the prime count; peak
memory does not (one prime at a time, `crt.py:286-301`).

### 2.3 Exact byte projection for `(5,5,5)`

Degrees [COMPUTATION]: physical bonds 300 (section 1.2); ghost bonds by the engine formula
(`crt.py:232-235`): `sum(ghost_inplane) = 4 corners*2 + 12 boundary-non-corner*1 + 9
interior*0 = 20`, hence `20*5 + 2*25 = 150`; full degree `300 + 150 = 450` (451 columns).
Truncated at the target order 54 (the e124 route, `e124:183`): retained 54, 55 columns.
Substituting `ns = 25` into the engine's own closed formula (`crt.py:425`):

    estimator model (3 arrays):  3 * 2^25 * 55 * 8 = 44,291,850,240 bytes = 41.25 GiB
    mechanical 2-array peak:     2 * 2^25 * 55 * 8 = 29,527,900,160 bytes = 27.50 GiB
    full-degree counterfactual:  3 * 2^25 * 451 * 8 = 363,193,171,968 bytes = 338.25 GiB
    (+ subleading scratch 805,306,368 bytes, section 0)

Anchors certifying the formula on STORED values at cs=20: `4x5x5` truncated
`3 * 2^20 * 53 * 8 = 1,333,788,672` (`results/series/lt_x34.json`
`maximum_estimated_peak_box.estimated_peak_bytes_at_target_truncation`, `e124:721-723`) and
full `3 * 2^20 * 366 * 8 = 9,210,691,584` (`e124:724`; `proofs/lt_x34.md` section 6).
The formula cannot be evaluated by `estimate_box_peak_bytes` itself at `ns=25` — the
preflight raises first (`crt.py:412` -> `:151-155`) — so 41.25 GiB is the closed formula
extended by one substitution, explicitly the "unlaunched input-size estimate" kind of
figure the producer records (`e124:134`).  41.25 GiB is 43% of the 96 GiB envelope and
below the 51.5 GB machine memory quoted in `proofs/ht_v24.md`.  Non-transfer overhead
(FLM inversion, cached Fractions) scales from the observed `x^52` process peak
1,809,235,968 bytes (`results/series/lt_x34.json` `resource_usage.process_peak_rss_bytes`)
— second-order here.
Runtime is deliberately NOT certified: the observed `x^52` all-box wall was 39.596147 s
(same key, with `4x5x5` dominant); `(5,5,5)` has 32x its states, 5/4 its primes, the same
layer count (5) — a guarded-lift run is plausibly tens of minutes as one governed process
[CONJECTURE — order-of-magnitude, not a measurement].

### 2.4 Why `x^54` is the order that matters

[COMPUTATION, quoted] Exactly seven truncation-unobservable ansatz survivors become
testable at `x^54` and no earlier: six algebraic rows (deg_f = 9..14, deg_x = 0; next
testable `u^27`) and one differential-algebraic row `(0,0,13)` (next testable `u^26`) in
`results/series/lt_x34.json` (re-audit rows with `required_coefficient_order_x: 54`) and
`proofs/lt_x34.md` sections 5-6.  The LT program's clean boundary closes exactly here.

## 3. Recommendation

**Cheapest exact next coefficient: LT `x^54`, via a guard lift, exact and in budget.**
The sole blocked input is `W(5,5,5)` (section 2.1).  Raising the single constant
`_MAX_CROSS_SECTION` 22 -> 25 (`crt.py:47`) after the value-safety audit of section 0
(`2^56 < 2^63 - 1`, headroom 128x) makes the EXISTING `e124` pipeline launchable with
projected peak <= 44,291,850,240 bytes = 41.25 GiB estimator model (29,527,900,160 bytes
mechanical two-array) — inside 96 GiB with 57% headroom.  The exactness stack is unchanged:
five primes, product 155 bits `> 2 * 2^125`, truncated prefix, sum-check explicitly
unavailable-but-unneeded exactly as recorded for the thirteen `x^52` CRT boxes
(`proofs/lt_x34.md` section 2).  A lift must still re-pass the engine's brute-force
cross-validation (`tests/test_tm_vs_enumeration.py`) on small boxes.

**HT `v^24` needs no new work, and the two classical reductions are certified dead ends
for the 8 GiB budget.**  The 242 GB projection is confirmed verbatim from stored metadata
(section 1.2); even-odd bipartite splitting shrinks the dense layer state by exactly 1
(and the flip sector by 2 for free boxes), landing at >= 20.6 GiB in the best case;
row-cyclic shrinking is not a symmetry of the open FLM box, and even tolerated as a
wrong-object counterfactual it lands at 8.25 / 9.03 / 45.15 GiB (section 1.4).  The
recorded working route was the algorithm change to the cut-capped parity frontier
(352 MB at `v^24`; section 1.5), already carried to `v^28`.

**No certified no-go of both walls is warranted, because an in-budget exact route exists.**
The genuinely open next HT wall is `v^30`: its new minimal-span box `(6,6,6)` has
cross-section 36, 15 cuts, dense projection `892,253,685,940,224` bytes [COMPUTATION], so
only the frontier family applies; per-state storage at the measured 204.0 / 300.1 / 247.8
bytes per state (`v^24`/`v^26`/`v^28` telemetry) is affordable only if peak-state growth
stays near the observed range, and no width-36 probe has been run [UNRESOLVED].  If a wave
wants one coefficient soon, schedule `x^54`; if it wants to invest an engine, the `v^30`
frontier at `(6,6,6)` is the next real risk point, not `v^24` and not `x^54`.

## 4. Scope and walls of this note

* Section 1.3's no-streaming statement is structural about THIS dense array engine family
  (full-polynomial layer transfer), not a lower-bound theorem over all exact algorithms.
* The `2^56 < 2^63` audit certifies only VALUE safety of a raised guard (no int64 overflow);
  it is not a re-run of the engine and not a claim about wall time.
* All "fits in budget" statements are input-size projections in the engine's own estimator
  model, anchored at cs=20 by two stored values; nothing at cs>=23 was launched.  This is
  the same epistemic status as the producer's own "unlaunched input-size estimate"
  (`e124:134`).
* No external value was used as an input, fit target, or selection criterion; `K_c` was not
  consulted anywhere.
