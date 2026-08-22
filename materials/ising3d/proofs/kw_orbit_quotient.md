# Exact orbit audit and complete `F_5`-unit census for Kac--Ward branch `11111`

Artifacts: `experiments/e169_kw_section_group.py`,
`experiments/e170_kw_unit_census.py`, `experiments/e171_kw_orbit_quotient.py`,
`results/kac_ward/orbit_quotient.json`, and
`tests/test_kw_orbit_quotient.py`.

## 1. Scope and verdict

**[THEOREM]** On the gauge- and diagonal-normalized 13-coordinate thin chart,
the only exact symmetry that induces an action on the particular six-solve / seven-held
section fibration of `experiments/e156_kw_full_thin_census.py` is the identity.
Consequently the `4^7=16384` sections have exactly **16384 singleton orbits**:
the proposed orbit quotient is valid but gives reduction factor one.  This is a
negative orbit theorem, not a missing optimization.

**[COMPUTATION]** Exact character evaluation nevertheless completes the entire
`4^13=67108864`-point `F_5^times` grid.  The six construction primitives have
exactly **2960** common unit zeros, lying in **2437** of the 16384 sections.
Every one of those 2960 local points has an exact independent-holdout disposition.

**[UNRESOLVED]** Branch `11111` is not proved empty over `Q`.  The census fixes
the seven held coordinates *exactly* at their least positive representatives
in `{1,2,3,4}`; it does not cover other higher 5-adic digits in those held
coordinates, non-unit valuations, zero reductions, chart-degenerate points,
`Q_p` for `p != 5`, or arbitrary characteristic-zero points.

## 2. The exact acting groups

### 2.1 Direction-state gauge

**[LEMMA]** Let `U(i,j)` be the 30 allowed direction-transition weights.  The
raw direction-state gauge group is

\[
 (a_0,\ldots,a_5): U(i,j)\longmapsto a_i^{-1}U(i,j)a_j,
 \qquad (\mathbb G_m)^6/\mathbb G_m.
\]

**[LEMMA]** Every closed-walk monomial contains each entering and leaving direction factor
with equal multiplicity, so all gauge factors cancel.  Therefore every trace
equation in the 42-equation construction catalog and every holdout trace equation
is gauge invariant.

**[LEMMA]** The wave-4 tree
`((0,2),(0,3),(0,4),(0,5),(2,1))` is a connected tree on six direction states.
On its common nonzero chart, setting its five weights to one determines all ratios
`a_j/a_i` uniquely.  The only remaining common factor acts trivially.  Hence the
normalized chart has **no residual tree-gauge sign** or other nontrivial gauge action.

### 2.2 The diagonal torus

**[THEOREM]** Before diagonal normalization, the six primitive thin numerators
have support-difference rank 13 in the 16 free coordinates.  An exact integer
kernel basis is

```text
v1 = ( 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0,  1, 0, 0, 0)
v2 = ( 0, 1, 1, 1, 1,-1,-2, 0,-1, 0, 0, 0,  0,-1, 0, 0)
v3 = ( 1, 0,-1,-1,-1, 2, 2, 2, 2, 1, 1, 1,  0, 0,-1,-1)
```

**[COMPUTATION]** On `(u_px_px,u_py_py,u_pz_pz)` its character matrix is

```text
D = ((1,0,0), (0,1,0), (1,-1,1)),    det(D) = 1.
```

**[THEOREM]** Thus the map from the three torus parameters to the three diagonal rescalings
is a unimodular torus automorphism.  After all three diagonals are fixed to one,
its stabilizer is the identity over every unit group.  The orbit-union lemma in
`proofs/kw_components.md` remains correct on the 16-coordinate pre-normalized
torus, but it supplies no second quotient after passing to the 13-coordinate
normalized chart.

### 2.3 Cubic point group and the section fibration

**[LEMMA]** The full lattice point group consists of 48 signed axis permutations.
A signed axis permutation sends the construction equation labelled `(shape,k)`
to the equation with correspondingly permuted side lengths and the same `k`.
The exact 14-shape catalog is invariant only for axis permutation `id` or the
swap `x<->y`; arbitrary independent axis signs do not change side lengths.
Therefore the construction-catalog point group is

```text
P_cat = (C2)^3 semidirect <swap(x,y)>,    |P_cat| = 8*2 = 16.
```

**[LEMMA]** Relabelling directions bijects closed nonbacktracking walks and box edges, while
the unique tree-gauge restoration is a trace-preserving similarity.  Hence these
16 maps preserve emptiness and nonemptiness of the exact construction variety.
The other 32 cubic maps are not symmetries of this finite 14-shape catalog and
are not used.

**[LEMMA]** Write `pi` for projection to the seven e156 held coordinates.  A
catalog symmetry induces a map on sections exactly when every Laurent monomial
in `pi(g.x)` is independent of the six solve coordinates, equivalently when
`pi o g` factors through `pi`.  Exact Laurent-monomial reconstruction gives one
such element, the identity.  For each of the other 15 elements the table gives
a held output with a nonzero solve-coordinate exponent; this is an exact witness
that a whole fiber is not sent to a fiber.

| axis permutation | signs `(x,y,z)` | held output | solve dependence witness |
|---|---:|---|---|
| `012` | `++-` | `u_mz_px` | `u_pz_px` |
| `012` | `+-+` | `u_pz_mx` | `u_my_mx` |
| `012` | `+--` | `u_pz_mx` | `u_my_mx` |
| `012` | `-++` | `u_pz_mx` | `u_my_mx` |
| `012` | `-+-` | `u_pz_mx` | `u_my_mx` |
| `012` | `--+` | `u_pz_py` | `u_my_px,u_pz_px` |
| `012` | `---` | `u_pz_py` | `u_my_px` |
| `102` | `+++` | `u_pz_mx` | `u_py_px,u_py_pz` |
| `102` | `++-` | `u_pz_mx` | `u_py_px,u_py_mz` |
| `102` | `+-+` | `u_pz_mx` | `u_py_pz,u_my_px` |
| `102` | `+--` | `u_pz_mx` | `u_py_mz,u_my_px` |
| `102` | `-++` | `u_pz_mx` | `u_py_mz,u_my_px` |
| `102` | `-+-` | `u_pz_mx` | `u_py_pz,u_my_px` |
| `102` | `--+` | `u_pz_mx` | `u_py_px,u_py_mz` |
| `102` | `---` | `u_pz_mx` | `u_py_px,u_py_pz` |

**[COMPUTATION]** The stronger variety-restricted obstruction also holds.  For
each of the 15 nonidentity catalog elements, the artifact records two of the
2960 exact `F_5` construction zeros in one input section whose transformed
points are construction zeros in two distinct target sections.  Thus no
nonidentity element induces a section map even after restricting from the
ambient thin chart to the construction zero set.  The standalone verifier
recomputes all 15 paired counterexamples with its independent Laurent action.

**[THEOREM]** Combining the trivial diagonal stabilizer, the absence of residual
tree gauge, and the identity fibration stabilizer inside `P_cat`, the group
acting soundly on the section set is `G_sec={1}`.  Every orbit is a singleton,
so Burnside gives

```text
|{1,2,3,4}^7 / G_sec| = 4^7 = 16384.
```

**[THEOREM]** Testing one representative per orbit is sound, but it is exactly the unreduced
16384-section census.

## 3. Exact complete census without a nontrivial quotient

**[LEMMA]** Since `2` generates `F_5^times`, write each coordinate as
`x_j=2^{a_j}`, `a_j in Z/4Z`.  A Laurent monomial `x^m` evaluates as the
character `2^{m.a}`, so every exponent may be reduced modulo four.  The six
primitive equations are therefore six exact character sums on `(Z/4Z)^13`.
No floating point enters this reduction.

**[COMPUTATION]** `e170` evaluates chunks of `2^20` exponent vectors, first the
7-term primitive and then the 33-, 155-, 155-, 155-, and 277-term primitives.
The largest unsigned phase dot product is `13*3*3=117<256`; the largest signed
coefficient accumulation bound used by the code is `277*4*4=4432<2^31`.
The complete grid produced 2960 zeros.  The independent verifier uses generator
`3`, chunks of `2^19`, and a different equal-size primitive order, and obtains
the identical 2960-point set.

**[COMPUTATION]** The exact per-point kill table is

| disposition | number of local points | exact certificate |
|---|---:|---|
| `KILLED_H8_MOD5` | 2232 | nonzero order-8 matrix-power residue modulo 5 |
| `KILLED_H12_MOD5` | 664 | order 8 is zero, order-12 matrix-power residue is nonzero modulo 5 |
| `KILLED_H8_MOD625_NONSINGULAR` | 45 | unique six-variable Hensel lift; all 42 construction residues vanish and order 8 is nonzero modulo 625 |
| `KILLED_BY_SINGULAR_LIFT_OBSTRUCTION` | 19 | exhaustive affine correction trees modulo 25, 125, and at most 625 leave no construction-plus-order-8 lift |
| **total** | **2960** | every local point covered exactly once |

**[THEOREM]** For a direct nonzero residue modulo 5, every 5-adic lift of that
reduction fails the named holdout.  For each nonsingular survivor, multivariate
Hensel uniqueness in the six solve coordinates reduces the exact fixed-held
slice to its recorded lift, whose order-8 residue is nonzero modulo 625.  For
each singular survivor, every affine digit correction compatible with all 42
construction equations is enumerated until none also has zero order-8 residue.
Therefore, in each of the 2437 exact integer-held sections with local points,
no 5-adic unit point in that fixed section satisfies all recorded holdouts.

**[THEOREM]** In each of the other 13947 exact integer-held sections, exhaustive
evaluation of all `4^6=4096` solve-unit values finds no construction reduction.
Thus those sections have no 5-adic construction point with unit solve coordinates.
They are nevertheless labelled `UNDECIDED`, not `EMPTY_OVER_Q`, because this
local absence is not a characteristic-zero emptiness certificate.

## 4. Per-representative classification table

**[COMPUTATION]** Representatives are ordered lexicographically in values
`1,2,3,4`.  For a section `h=(h0,...,h6)`, put `d(1)=0`, `d(2)=1`, `d(3)=2`,
`d(4)=3` and

\[
i(h)=\sum_{j=0}^6 d(h_j)4^{6-j}\in\{0,\ldots,16383\}.
\]

**[COMPUTATION]** The following 2048-byte bitmap is the complete per-representative classification
table.  Decode the displayed hex to bytes; representative `i` has status
`HAS_LOCAL_POINT_BUT_FAILS_HOLDOUT` exactly when bit `i mod 8` of byte
`floor(i/8)` is one, and has status `UNDECIDED` when that bit is zero.  There
are no `EMPTY_OVER_Q` rows.  Its raw-byte SHA-256 is
`96efae5ef53fe4bc23e2a403d6ceca01b82b94a8f350dd1898d8bf5707b042ef`.  The same 16384 rows appear explicitly as JSON objects at
`data.classification_table` in the artifact, so every tuple, count, point-id list,
and certificate mechanism is machine-readable.

```text
0084640008480008005804080060000000800008044101c000801000404000004000208108800010160808090010000014410000000400443258000900502020
00b08a18008000000605013810d0062486a160001104000080354100104400080000000082120100380cc080020110a0281d020040120000008448c3c8010802
0a400006800000000205004408510002024920400000000100c0004000402080800600608802000008d0404001000008081480c0012000000204c06e04008020
02100c0010c020600000000000484002000000300070000000400c588000000104290000808210102428041840000000002d4000020808020001408104908000
8ee240000044008494210928000a80019481400008c400008820080000c8104000010000a4008040205c00132800000200008420a03000400008001208320400
ac300001020110083832140cc069240029140101020000022404180122c140000a4009800810000099428d420880080282420d02372010288880240362024000
0102406890008060257440c00042000064543904020400000e8400080040401004c1008400004800804042284431040084456008006100008415481200000020
00820140ac0020a412284020090d0000426901841000814031780480102000022d510000a0380008146800000050000001000000800025080100082201000840
04580180400260070441028000880001040d0090000880440008018048060240002005401408c000040020004000000000000000b41f800808000000040800c0
0008002092028006c14260a000020010000e006080011306004802008221001002b100140220d002400800281220000000040000000054144e09200040800400
0010000000412084000460010240041400040000401520540100000002550010043828240060882c00001204104800001300000400820020000074044000c024
2420340241a000002450000c00044f00a921000000088021881100042181400000028008152cc0000c000001c421804224080001000000042000800524180488
880080a01420083002089000c4a489408206800000202000020c80001020010028880080a9c40000000000c4206a010410080000205401000022400401c00059
0000600082204c000020e000060106041044201201480104301028024000020020310800c28021884514408042c000820620000069c007848330800082480244
4000008401c0004084020040818002cc40020100c0001024042209000080441060020c0017e50080202024008001080104500400800100010081000180410410
4810004010288080205800020040400830000400510801103040025210203009a8600100614200000a1801130c40400200000102025000040000050ac0420033
2008404a01800040204080e018888008040000000080040000804016018040020014040800000000800c100820020818081008040000000008081cac00000019
2000820580a08004816882a10020000000008000002000b00002820000200080104a000800800620000a01000148000023000002000010040849008000908000
0400024101090081044000000000000040002e600000080204588041000000480008800820144000000a020000040100008840804014830000ce000222070041
80001399008000810824e048008c00100100280c20000000080801882002400008848699001120001000042400010c232052042200018400c812042804000080
186052248000008408041000a08648800100042080080134000494201000021000110782620102a1000200000000004d88080000000089102800258000040010
001112002003358002411000200085240400c120200002010411000360103ec94000e24400080800814400080090084800000342000101000081830059080300
820090002804082018a2120080482000ca00000000000c04048a00822000041020208841006010a040e884270440800b01064000040011010082284104200800
439030a6001000001135108a08a000905105209a00008921440110680020c8020000610a0100e002000100194d4808008401081849000060c60c403a41411960
0040400b2046008400000048200008c001800405002d0800080000020088089000000400000414180000000c0004140404004010200e14080000010808010088
008400348810d0a904200000011482210000d804000002410002080c8004822080a80000006861040800602012dca00d00200005803621060000000100011111
08010000082000030000084000241241000c00000820060180000c1000803a480000004180144004004000000000089c00001100000008125284000800841814
0024a0d000000480400000044080284800100b0001410000000848002000000400126200208004818001a04080c2002900000f01054404282041000008200c21
000000202024808010010122804094000000400012061428400000c040229c249008808000201880400608089080000008000500010020000040040000402000
024021800c1320000001810600012a100000140c0004491404003c180400200520840002102a4080000020400300824005088010030494490000020800002a48
00001c81420221008210104000002244000000000203280400001433045100800a00600000928401002000100000844140102001080480402448208084e10400
00020010208210184092001a018841380241830110001428010011000004008000400028000025080000104044028349194100080425401200008218d4450100
```

**[COMPUTATION]** The exact status totals are

| status | representatives | exact membership rule |
|---|---:|---|
| `EMPTY_OVER_Q` | 0 | no certificate claimed |
| `HAS_LOCAL_POINT_BUT_FAILS_HOLDOUT` | 2437 | one-bits of the bitmap |
| `UNDECIDED` | 13947 | zero-bits of the bitmap |

**[COMPUTATION]** Thus the sections that remain are **exactly** the 13947 tuples whose decoded bit
is zero; their lexicographic indices are also listed explicitly at
`data.classification_summary.remaining_undecided_lex_indices`.

## 5. Scope of the theorem and what remains

**[THEOREM]** For every one of the 16384 slices with the seven held coordinates
fixed exactly at a tuple in `{1,2,3,4}^7`, no `Q_5` point with all six solve
coordinates 5-adic units satisfies both the 42 construction equations and the
recorded order-8/10/12 holdouts.

**[THEOREM]** The construction ideal itself is proper: the prior anchor Hensel
point is an exact `Q_5` construction point.  Consequently no rational
Nullstellensatz identity can certify that construction ideal empty at any
multiplier degree.  Future work cannot close the construction branch by merely
raising that certificate degree.

**[UNRESOLVED]** The global branch question remains open.  In addition to the
13947 rows conservatively labelled `UNDECIDED` over `Q`, this finite family of
least-representative held slices does not cover higher 5-adic held digits,
non-unit valuations, zero reductions, chart-degenerate loci, other primes, or
arbitrary characteristic-zero points.  A process-time cap or missing local
point is never promoted to an `EMPTY_OVER_Q` theorem.

## 6. Reproduction and independent verification

**[COMPUTATION]** Run from the repository root:

```bash
.venv/bin/python experiments/e169_kw_section_group.py
.venv/bin/python experiments/e170_kw_unit_census.py --cap-seconds 120
.venv/bin/python experiments/e171_kw_orbit_quotient.py --cap-seconds 120
.venv/bin/python tests/test_kw_orbit_quotient.py
```

**[COMPUTATION]** In the final recorded producer run, exact character enumeration took
`8.535136999999999` process seconds and the complete decision path
took `24.558994000000002` process seconds.  The e156 preflight had
projected `100002.627584` process seconds from a `6.103676`-second first section;
the observed projection-to-character-enumeration ratio is
`11716.581419138323`.  These are host-specific measured costs, not mathematical
claims.  The producer's canonical volatile-free body digest is
`088706b51dea683c801e3648ff010bb68f2da8997146a15169c1b1e8aa529142`.
The artifact also preserves the earlier complete PASS observation
`8.650371 / 24.405328` process seconds and ratio `11560.50157663758`, recorded
before adding the 15 variety-fibration counterexample rows.

**[COMPUTATION]** The standalone verifier imports none of e169/e170/e171.  It
rebuilds the 42 equations, recomputes the order-16 catalog group and the identity
section stabilizer, enumerates all 67108864 unit points independently, verifies
all 16384 classification rows, and replays four distinct kill mechanisms plus an
`UNDECIDED` no-point section.  It prints `PASS` only after universal coverage.
