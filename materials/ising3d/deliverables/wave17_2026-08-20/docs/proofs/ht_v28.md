# Exact simple-cubic HT coefficient through v^28

## Scope

[COMPUTATION] This is a finite exact coefficient calculation, not a solution of
the three-dimensional Ising model.

## Exact order-28 box-class bound

[LEMMA] Enting's finite-lattice identity is
`L(A)=sum_(R<=A) product_i(A_i-R_i+1) W(R)`.  Every edge of an
`a*b*c` box crosses exactly one of its `a+b+c-3` coordinate cuts.  An even
subgraph crosses each cut evenly; an exact-bounding-box connected contribution
crosses each cut positively.  Hence it has at least
`2*(a+b+c-3)` edges.  At truncation order 28 this proves the exact class bound
`a+b+c <= 17`, not merely a heuristic cutoff.  The artifact enumerates all
147 canonical classes, of which
24 have sum 17 and first occur at v^28.

[COMPUTATION] The six and only six classes whose canonical open spin-transfer
cross-sections exceed 22 are `[[4, 6, 6], [4, 6, 7], [5, 5, 5], [5, 5, 6], [5, 5, 7], [5, 6, 6]]` with cross-sections
`[24, 24, 25, 25, 25, 30]`.  Their dense broken-bond projections are
`[140525961216, 165490458624, 242397216768, 294742130688, 347087044608, 11467562680320]` bytes, so their weights are instead computed
by the sparse parity frontier.

## Generalized cut-profile frontier

[LEMMA] For a requested positive half-crossing profile `alpha`, cap cut `i` at
`2*alpha_i` crossings and retain the complete coordinatewise quotient with
radix `alpha_i+1`.  The pending parity register has width `b*c` in a sweep of
an `a*b*c` box; it enforces even vertex degree exactly.  The formal-log
recurrence in this ordinary mixed-radix quotient includes all decompositions of
`alpha`, including a two-plus-two split of a four-crossing cut.  Counts above a
cap cannot enter the target coefficient, so the cap is exact.

[LEMMA] At degree 28, a sum-17 box has the all-ones profile; each sum-16 box
has exactly one exponent two; and the 5x5x5 box has either one exponent three
(one six-crossing cut) or two distinct exponents two (two four-crossing cuts).
The latter requires exponent radix four.  Box reflections and permutations of
equal sides give two one-six-crossing cube orbits (edge/middle, multiplicity
six each) and seven two-four-crossing orbits: same-axis opposite edge-edge
(multiplicity 3), same-axis middle-middle (3), same-axis edge-middle adjacent
(6), same-axis edge-middle nonadjacent (6), and different-axis edge-edge (12),
middle-middle (12), edge-middle (24).  The cut-coordinate distance distinction
is retained; no inequivalent profiles were merged.

## Primitive log profiles

[LEMMA] Let `P` be an integer-coefficient formal power series with constant
term one.  If a multi-index `alpha` is primitive
(`gcd(alpha_i)=1`), then `[x^alpha] log P` is an integer.  Indeed, in the
`k`-fold term of the logarithm expansion, cyclic rotation acts on the ordered
`k` factors.  A nontrivial stabilizer would repeat a shorter word and hence
divide every coordinate of `alpha`; primitivity makes every orbit have size
`k`, cancelling the displayed denominator `1/k` orbit by orbit.

[COMPUTATION] Every listed v^24, v^26, and v^28 cut profile has gcd one.  Thus
the reconstructed centered residues are integral formal-log coefficients before
their exact symmetry multiplicities are applied.

## CRT certificate

[LEMMA] Put `O_q=sum_(k=1..q) S(q,k)(k-1)!`.  For a profile `alpha`, label its
`q=sum_i alpha_i` half-crossing slots.  Ignoring vertex parity, a slot on a cut
with `E_i` crossing edges is bounded by `C(E_i,2*alpha_i)`.  Expanding the
formal logarithm over ordered set partitions gives the strict bound
`O_q product_i C(E_i,2*alpha_i)^alpha_i`.  The profile-wise values are stored
verbatim in the JSON artifact.  The product of the six independently primality
checked 31-bit primes is `98079699360994458463449574431304277015588525938982026813`, and it exceeds twice every
stored profile bound; centered reconstruction is therefore unique.

## Lower-order reproduction and independent controls

[COMPUTATION] Before any degree-28 extraction, the new cap-vector pipeline
recomputed the v^24/v^26 walled weights and fresh FLM assemblies.  It reproduced
`[v^24] phi = 2135670379057/8` and
`[v^26] phi = 115377512914251/26` exactly.  Independent
spin-transfer FLM controls for cap-2, cap-4, and the new cap-6/pair profile
families are `{'cube333_degree16': {'aggregate': '-295914', 'spin_transfer': '-295914', 'profile_orbit_count': 3}, 'cube444_degree18': {'aggregate': '8655072', 'spin_transfer': '8655072', 'profile_orbit_count': 1}, 'cube444_degree20': {'aggregate': '159872736', 'spin_transfer': '159872736', 'profile_orbit_count': 2}}`.

## Exact coefficient

[COMPUTATION] The injected walled weights are `{'[5, 5, 5]': {'24': '8163299968', '26': '310476131424', '28': '6185993200560'}, '[5, 5, 6]': {'26': '66848887584', '28': '3065159677392'}, '[4, 6, 6]': {'26': '37481891632', '28': '1815705969248'}, '[4, 6, 7]': {'28': '232312855136'}, '[5, 5, 7]': {'28': '412031318992'}, '[5, 6, 6]': {'28': '651743244920'}}`.
The exact finite-lattice assembly gives

`[v^28] phi = 2102198327465307/28`,

with interaction coefficient `a_28 = 525549581866326/7`.  The odd
coefficient v^27 vanishes and every coefficient through v^26 reproduces the
previous exact artifact.

[EXTERNAL] Arisue and Fujiwara, arXiv:hep-lat/0209002 Table I tabulates
`a_28 = 525549581866326/7`.  This value was used only after
computation as a falsification witness; it agrees exactly.  Recorded PDF sha256:
`bf79ec7918694768b57e6793ed4b0c8f581aaea6a705b4273ecab42641723afc`.

## Resource telemetry

[COMPUTATION] Frontier workers were isolated.  Their summary is
`{'profile_count': 27, 'max_peak_states': 29378975, 'max_worker_rss_bytes': 7280050176, 'sum_worker_wall_seconds': '7270.996971'}`.  The producer peak RSS was
`22577053696` bytes and elapsed wall time was
`10584.737070` seconds.  No disk checkpoint was used.

[COMPUTATION] The finalized worker samples RSS while both source and
destination frontier dictionaries coexist and stops with an `[UNRESOLVED]`
certificate at `36000000000` bytes rather than relying on an
allocator failure.  This post-production guard is not used as evidence for the
coefficient above.
