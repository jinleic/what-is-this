# Octahedral (order-48) compression of the extensive-charge classes: feasibility, execution, and the exact s<=6 wall

Artifacts: `experiments/e144_octahedral.py`, `tests/test_octahedral.py`,
`results/integrability/octahedral_s6.json`.

Reproduce:

```
PYTHONPATH=src .venv/bin/python experiments/e144_octahedral.py
PYTHONPATH=src .venv/bin/python tests/test_octahedral.py
```

**Status convention.** `[THEOREM]`/`[LEMMA]` are proved here (the cap lemma is
elementary linear algebra); `[COMPUTATION]` is the exact finite certificate;
`[CONJECTURE]` marks the identified-but-unimplemented follow-up route;
`[UNRESOLVED]` is what is deliberately not decided. All certified quantities
are exact integers or rational arithmetic on integers; no floats.

---

## 1. Setting

`H = sum_x X_x + sum_(<xy>) Z_x Z_y` on `Z^3`, `D = (1/2)[H, .]`, `pi` the
translation-orbit sum. The wave-12/13/14 framework gives, for the
support-size class `C(3,2,s)` (ordered-Pauli words supported in
`B_2 = {0,1,2}^3`, word support size `<= s`),

```
Q(C) = rank S_C - rank M_C - 1        (rank S = #anchored class words,
                                        rank M = rank of pi D|_C)
```

with the certified values `Q = 1` for `s <= 4` (wave 13) and `s <= 5`
(wave 14), by the sandwich `rank_Fp <= rank_Q <= rank S - 2` (the analytic
upper end holds because the Hamiltonian density `h = X_0 + sum_i Z_0 Z_(e_i)`
lies in the class and `pi D h = 0` exactly).

The group of this note is `G = O_h`, the 48 signed coordinate permutations
(the point group of the `Z^3` lattice), acting on the box about its centre;
it preserves `B_2`, the bond structure, and `H`, hence commutes with `D`,
commutes with translations, and preserves the support-size stratification
(Lemma A of `proofs/extensive_3d_r2.md`). `G` acts on word columns by site
relabeling (letters ride along) and on anchored rows by the induced action.

## 2. Feasibility: exact Burnside inventory

**[COMPUTATION — sizes]** `|C(3,2,6)| = sum_{k<=6} C(27,k) 3^k = 236,912,446`
columns (`k`-th layer `C(27,k) 3^k`; per-size values in the artifact), and
`rank S_6 = 1 + sum_{k<=6} A_k 3^k = 191,949,610` rows with
`A_k = C(27,k) - 3C(18,k) + 3C(12,k) - C(8,k)` (inclusion-exclusion and
brute force forced equal per `k <= 6` in the producer).

**[THEOREM B1 — column orbits]** The number of `G`-orbits on the word class
`C(3,2,s)` is `(1/48) sum_g fix_s(g)` where `fix_s(g)` counts letter
assignments constant on the site cycles of `g` with total cycle length `<= s`
(a DP over the cycle multiset: a taken cycle contributes its length and a
factor 3, one of `{X, Z, XZ}`). This yields exactly

| s | columns | column G-orbits | factor |
|---|---:|---:|---:|
| 4 | 1,503,766 | 36,035 | 41.73 |
| 5 | 21,121,156 | 464,957 | 45.43 |
| 6 | 236,912,446 | **5,048,368** | 46.93 |

*Validation.* The Burnside sums are integral at every `s <= 6`; an
independent brute-force orbit count by canonical forms agrees at `s <= 3`
(2,456 orbits); on the 8-site box `B_1 = {0,1}^3` (half-integer centre) the
same machinery is validated by FULL brute force through `s = 8` (1,996
orbits). The canonical-form representative enumeration (`iter_representative_columns`)
emits exactly the Burnside counts at `s <= 2,3,4` and both representative
echelons consumed exactly 464,957 resp. 5,048,368 columns.

**[THEOREM B2 — row orbits]** The induced `G`-action on anchored words (the
`rank S` rows) fixes a nonempty anchored word `w` iff `g.w = tau_t . w` for
a unique translation `t`, iff the support of `w` is invariant under the
affine map `A_t(y) = g y - t` (centred coordinates) with letters constant on
the finite `A_t`-cycles lying inside the box. Summing a per-`(g, t)` DP over
`t in {-2..2}^3` and dividing by 48:

| s | rows (rank S) | row G-orbits | factor |
|---|---:|---:|---:|
| 5 | 14,757,412 | 326,710 | 45.17 |
| 6 | 191,949,610 | **4,094,477** | 46.88 |

*Validation.* At `(g = e, t = 0)` the DP must and does return exactly
`1 + sum_k A_k 3^k` rows for every `s <= 6` (191,949,609 nonempty at
`s = 6`); brute force by re-anchoring canonical forms agrees at `s <= 3`
(1,024 orbits); `B_1` rows validated by brute force through `s = 4` (269).

## 3. The cap lemma and the decision

**[LEMMA CAP]** A column echelon over any field, fed only a set of `N`
columns, certifies rank at most `N`. Consequently an orbit-representative
echelon for `C(3,2,6)` certifies at most `N_orb = 5,048,368 < R* :=
rank S_6 - 2 = 191,949,608` pivots — i.e. at most 2.63% of the rank needed
for the quotient-1 certificate — no matter how the representatives are
chosen or how cheap the memory. *Proof:* the rank of a matrix is at most
its number of columns. ∎

**[COMPUTATION — representative echelons, two primes]** The representative
echelon was nevertheless RUN (the reduced working set peaks at 1.4 GiB
RSS, far under the 8 GiB gate), reproducing the wave-14 engine exactly:

* `s <= 5` calibration (where the exact full-class answer is known): the
  464,957 representative columns have rank **460,647 at both primes**
  (4,310 columns lost to inter-representative coupling, 98.97% of them
  independent) — i.e. the representative echelon recovers 460,647 /
  14,757,410 = **3.12%** of the certified target rank. The lead-stream
  digests of the two primes coincide.
* `s <= 6` headline: the 5,048,368 representative columns have rank
  **5,036,417 at both primes** (`p = 2147483647` and `p = 2147483629`;
  digests coincide). Since the submatrix has integer entries,
  `rank_Fp(submatrix) <= rank_Q(submatrix) <= rank_Q M_6`, so this is a
  genuine **certified lower bound** `rank_Q M_6 >= 5,036,417` and hence

```
[THEOREM — certified bracket]  1 <= Q(C(3,2,6)) <= 186,913,192
```

  (lower end: `h` is in the class, wave-14 argument, `pi D h = 0` re-verified
  exactly by the producer's defect traps; upper end:
  `rank S_6 - 5,036,417 - 1`).


**[COMPUTATION — regression]** The reimplemented pipeline reproduces the
frozen wave-13 `s <= 4` certificate exactly (full class, both primes:
1,503,766 columns, rank 822,332 = `rank S_4 - 2`), and extends the same
sandwich to `s <= 2` (3,241 columns, rank 560) and `s <= 3` (82,216 columns,
rank 29,747), i.e. `Q = 1` holds for every `s <= 5` in this pipeline.

**[COMPUTATION — the full-set wall]** The full `s <= 6` certificate needs
`R* = 191,949,608` stored pivots. Calibrated on the frozen wave-14 `s <= 5`
run (14,757,410 pivots, 313,579,672 stored nnz, peak RSS 2,404,089,856 B):

* lead table alone (2^28 slots, load 0.715, 12 B/slot): 3,221,225,472 B;
* CSR pivot store (`R* x 21.248964 avg nnz x 12 B`): 48,944,763,437 B;
* **total model 52,165,988,909 B = 48.6 GiB**, conservative whole-process
  model (`R* x 163 B/pivot`) 31,269,992,868 B = 29.1 GiB — both far above
  the 8 GiB = 8,589,934,592 B gate;
* projected time at the frozen `s <= 5` rate (14,599 columns/s):
  16,227 s per prime (9.02 h for both).

**[THEOREM — wall certificate]** Under the audited certificate architecture
(wave-12/13/14 column echelon with stored pivot rows), the `s <= 6` class
quotient cannot be decided within 8 GiB: representative compression is
capped at 2.63% of the required rank by Lemma CAP, and the full-set
certificate needs ~30-49 GiB by both exact memory models. The `s <= 6`
quotient and, a fortiori, the full radius-2 `Z^3` box remain
**[UNRESOLVED]**, now with an exact wall certificate (no observed resource
wall was needed: the wall is arithmetic).

## 4. What actually survives the quotient (the "loosely coupled rows")

* Columns: 5,048,368 of 236,912,446 survive as representatives (1/46.93).
* Rows: 4,094,477 of 191,949,610 row orbits (1/46.88) — but a plain echelon
  cannot exploit this: its pivots are individual anchored keys.
* **The rank target does not shrink at all.** The image of `M_6` is
  `G`-stable of dimension `R*` (or `R* + delta` for the unknown excess
  `delta >= 0`); `G`-symmetry never reduces the dimension that a rank
  certificate must certify — this is the precise sense in which the
  commutator-image rows are "loosely coupled": they span a `G`-stable space
  of near-full row dimension rather than a small orbit-structured one. The
  measured `s <= 5` stall (see artifact) and the `s <= 3` sub-cap
  representative ranks (2,292 of 2,456) confirm the coupling is real:
  representatives lose rank strictly below their count.
* **[CONJECTURE — identified follow-up]** A `G`-equivariant elimination
  storing one pivot row per ROW orbit (about
  `4,094,477 x 21.248964 x 12 B ~ 1.0 GiB` of CSR) plus a lead structure
  over lead orbits would plausibly fit the 8 GiB gate with the lead table
  (~3.0 GiB) intact. It requires on-the-fly exact reconstruction of stored
  pivot rows under every group element (anchor-key transforms with the
  reduction history) and a `G`-compatible pivot order. This is engineering
  with nontrivial correctness obligations and is deliberately NOT attempted
  here; it is the concrete route this front contributes to the program.

## 5. Scope, resources, and provenance

Feasibility arithmetic is exact integer work (validated by brute force on
two boxes); the echelons are the wave-14 engine (column-oriented, streamed,
compact CSR pivot store excluding leads, open-addressing lead map, static
descending (support size, packed blob) lead priority — the pivot count is
lead-order independent). The only clocks are `time.process_time()`; RSS via
`ru_maxrss`. All stages ran sequentially under `nice -n 19` with single-thread
BLAS environments on the wave-15 box. K_c is untouched by this front. No
all-size theorem is claimed: every statement here is about the finite class
`C(3,2,s)`, `s <= 6`, or an elementary all-size lemma (Lemma CAP) about
column echelons.
