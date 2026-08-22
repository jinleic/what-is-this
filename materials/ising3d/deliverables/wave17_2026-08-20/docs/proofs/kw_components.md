# Component cleaning of Kac--Ward branch `11111` at the anchor thin-torus section

Artifacts: `experiments/e139_kw_components.py`,
`results/kac_ward/components.json`, and `tests/test_kw_components.py`.
This note continues `proofs/kw_cubic.md` / `proofs/kw_branch11111.md` and
the full-family study `proofs/kac_ward_full_family.md`.

## 1. Scope, chart, and exact computational convention

**[COMPUTATION]** We work throughout in the gauge-tree chart `11111`: the five
`DIRECTION_GAUGE_TREE` weights are set to 1, leaving the 25 weights
`u_px_px, ..., u_mz_mz`.  The construction catalogue is the same 42 exact
integer equations as the branch study: orders `4,6,8` on the fourteen
shapes

```text
permutations of (2,2,1), permutations of (3,2,1),
(3,3,2), (2,2,3), (2,2,2), (3,2,2), (2,3,2).
```

The *thin torus* is the chart in which the three diagonal weights
`u_px_px = u_py_py = u_pz_pz = 1` and the nine reflected weights are the
rational functions

```text
u_mx_my = -1/u_my_px,      u_mx_py = -1/(u_my_mx*u_py_px),
u_mx_mx = 1/u_px_px,       u_mx_mz = -1/(u_mz_px*u_pz_mx),
u_mx_pz = -1/(u_mz_mx*u_pz_px),         u_mz_mz = 1/u_pz_pz,
u_my_mz = -1/(u_mz_py*u_py_pz*u_pz_my),  u_my_pz = -1/(u_mz_my*u_py_mz*u_pz_py),
u_my_my = 1/u_py_py.
```

Under the thin reduction the primary-shape equations (indices `0..27`) either
vanish identically or reduce to a 13-weight quotient; the six independent
numerator constraints at indices `28,29,32,35,38,41` define the
13-coordinate primitive system `S`.  All arithmetic below is exact integer
residue arithmetic over `Z/625Z` (or `Z/5Z`); no floating point and no
`2^62`-style truncation decides any statement.  `[COMPUTATION]` tags mark
the enumeration steps; `[THEOREM]` the statements they certify.

## 2. The anchor section and its solution census

Define the *anchor section* by holding the seven trailing thin coordinates
(`u_pz_mx, u_pz_py, u_pz_my, u_mz_px, u_mz_mx, u_mz_py, u_mz_my`)
at the integers `(1,1,1,1,1,4,4)` and solving the six primitives in the
six remaining coordinates (`u_py_px, u_py_pz, u_py_mz, u_my_px,
u_my_mx, u_pz_px`).

**[THEOREM]** (solution census, unit stratum).  Over `F_5`, restricted to the
**unit stratum** — all six solve coordinates in `F_5^x` (the `4^6 = 4096`-point
grid) — the primitive system at the anchor section has **exactly three
solutions with all six solve coordinates in `F_5^x`**:

```text
p0 = (1,2,3,3,1,1 ; 1,1,1,1,1,4,4)
p1 = (2,4,1,3,3,1 ; 1,1,1,1,1,4,4)
p2 = (3,2,3,1,1,1 ; 1,1,1,1,1,4,4)
```

(Semantics: first six entries = solve coordinates.)  This is a finite fact: the
unit grid was evaluated exactly, exactly three points satisfy all six
primitives mod 5, and the count is recomputed independently in the standalone
test.  **Scope, per adversarial review:** over the FULL `F_5^6` grid the six
primitive numerators have **1553** solutions, four of them nonsingular; the
fourth nonsingular point has a zero coordinate (`q = (1,3,0,3,2,3)`), its
Hensel lift keeps `u_py_mz = 0` exactly through `mod 5^6`, so it is
chart-degenerate (the chart substitution divides by `u_py_mz`) and is not a
variety point of this chart.  The 1549 remaining boundary solutions are
singular and outside the unit stratum; they are unanalysed here.  The census
claim is for the unit stratum only.  The `6 x 6` Jacobian of the slice at the
three unit points has determinants `[COMPUTATION]`


```text
det J(p0) = 2 (mod 5),   det J(p1) = 4 (mod 5),   det J(p2) = 1 (mod 5).
```

The first value `det J(p0) = 2 mod 5` reproduces the wave-10 record
(`results/kac_ward/branch11111_cubic.json`), recomputed here from the raw
derivative tables of the reduced slice.

## 3. Nonsingular reduction and the three unique `Q_5` components

**[THEOREM]** (Hensel, unit stratum, integer-held section).  Because the
reductions of `p0, p1, p2` are nonsingular, the multivariate Hensel lemma
applies: each of the three `F_5` unit points lifts (uniquely in the solve six)
to a point of the full 42-equation variety modulo `5^4 = 625`.  Any `Q_5` point
of the variety over the anchor section **whose seven held thin coordinates are
the integers `(1,1,1,1,1,4,4)` and whose twenty-five weights are 5-adic
units** reduces to one of the three unit points and therefore equals the
corresponding lift.  (The unit hypothesis is the wave-10 Laurent convention —
mod-625 residues are only meaningful on that stratum — and the held-as-integers
hypothesis matches the section definition; slices held at values merely
congruent mod 5 are different polynomials with their own lifts, not covered.)
The three lifts (solve six, mod 625) are `[COMPUTATION]`

```text
L0 = (226, 197, 98, 63, 161, 216)      (the wave-10 anchor lift)
L1 = (157, 459, 431, 228, 38, 496)
L2 = (393, 72, 218, 596, 446, 496)
```

At each lift the full 42 construction residues were evaluated exactly mod 625 and
are identically zero:

```text
construction_residues_mod_625 = [0]*42   for L0, L1, L2.
```

So `L0, L1, L2` are genuine `Q_5` points of the construction variety, and
together they *exhaust* it over the anchor section in the chart.

### 3.1 The exact `3x3x3` holdout

The holdout is the open `3x3x3` box `3x3x3` exact trace system at
orders `4,6,8,10,12`.  For an order `k` the equation is
\[
 E_k(x) \;=\; \mathrm{Tr}(T^k)(x) + 2k\,\log P_k ,
\]
where `T` is the 30-direction (108 directed-edge) `3x3x3` Kac--Ward
transfer matrix and `P` its even-subgraph polynomial (`log P_4 = 36`,
`log P_6 = 164`, `log P_8 = 663`, `log P_10 = 2280`,
`log P_12 = 972`).  The *residue* is `E_k(x) mod 625`.

**[THEOREM]** (holdout killing of all three components).  At each of the
three `Q_5` components the residues `(k=4,6,8,10,12)` are
`[COMPUTATION]`

```text
L0 :  (0, 0, 250,  55, 266)
L1 :  (0, 0, 456, 350, 132)
L2 :  (0, 0, 430, 480, 286).
```

Hence at **order 8, 10, and 12** the residue is nonzero mod 625 at
*every* one of the three components; in particular order 8, the threshold of
Theorem S in `proofs/kw_cubic.md`, fails at all three.  Combined with
Section 3, the precise quantifier is:

> In the gauge-and-diagonal chart over the anchor thin-torus section
> `(1,1,1,1,1,4,4)`, there is **no** `Q_5` point `x` of the
> 42-equation construction variety and **no** order `k in {8,10,12}` with
> `E_k(x) = 0` mod 625.  The three anchor-section components `L0, L1, L2`
> are the complete `Q_5`-point set of the section in the chart, and each
> fails the order-8/10/12 holdout.

The order-8 value at the anchor, `E_8(L0) ≡ 250 mod 625`
(`= 2*5^3`), reproduces exactly the wave-10 first holdout reproduction
(`holdout_residues_mod_625 = [0,0,250]` in `branch11111_cubic.json`);
this is a mandatory control of the present verification.

### 3.2 Cross-validation of the symbolic and matrix-power routes

**[THEOREM]** (route agreement).  The residues are produced by two independent
exact routes and agree:

1. *Symbolic route.* Orders `4,6,8` use the exact
   `full_weight_finite_system((3,3,3), (4,6,8))` equations evaluated at
   the 25-coordinate lift mod 625.  Orders `10,12` use the *walk expansion*
   of `Tr(T^k)`: an exact monomial multiset (a closed directed walk contributes
   the product of its weight symbols, gauge factors = 1; 23288 monomials at
   order 12), evaluated mod 625, minus `-2k log P_k`.  (The direct
   symbolic `translation_invariant_trace_expression` at orders `10,12` is
   intractable -- `[UNRESOLVED]` at scale -- so the walk expansion is the
   symbolic route there.)
2. *Matrix-power route.* `Tr(T^k) mod 625` via binary exponentiation of the
   108 x 108 transfer matrix over `Z/625Z`.

The two routes agree at every order `4,6,8,10,12` for all three
lifts; the walk and matrix routes additionally agree at `10,12`, so the
stored matrix values `55/350/480` and `266/132/286` are certified by an
independent computation.

## 4. The orbit-union lemma

**[THEOREM]** (diagonal torus action).  On the thin torus, written in the
sixteen free coordinates `free_variable_order` (thirteen non-diagonal plus the
three diagonal weights),

```text
u_px_px, u_py_px, u_py_py, u_py_pz, u_py_mz, u_my_px, u_my_mx,
u_pz_px, u_pz_mx, u_pz_py, u_pz_my, u_pz_pz, u_mz_px, u_mz_mx,
u_mz_py, u_mz_my,
```

the support-difference lattice of the six primitive numerators has rank 13, and
its `Z^16` kernel has basis `[COMPUTATION]`

```text
v1 = ( 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0,  1, 0, 0, 0)
v2 = ( 0, 1, 1, 1, 1,-1,-2, 0,-1, 0, 0, 0,  0,-1, 0, 0)
v3 = ( 1, 0,-1,-1,-1, 2, 2, 2, 2, 1, 1, 1,  0, 0,-1,-1).
```

The square restriction of `(v1,v2,v3)` to the three diagonal slots
`(u_px_px, u_py_py, u_pz_pz)` has determinant `±1`, so the three
`Z`-directions generate a genuine rank-3 lattice acting by
`x_j -> x_j * t1^{v1_j} * t2^{v2_j} * t3^{v3_j}` and the diagonal
coordinates may be normalized to 1 within an orbit
(`[THEOREM]` chart normalization).

For each of the six primitive numerators `S_i`, and for each of the three
kernel directions, the dot product `m . v` of every exponent vector `m` with
the kernel vector is constant.  Hence each `S_i` is *semi-invariant* of the
torus action, so its `F_5` zero set is a union of `F_5^times x F_5^times x F_5^times`
orbits.  Consequently **the orbit scan collapses to one representative per
orbit**: for the census, one representative per orbit decides the recovery pattern of
the section.

**[THEOREM]** (holdout semi-invariance).  The order-8 and order-10 holdout
numerators `H6, H8` on the thin torus are semi-invariant with exact
characters

```text
H4 = 0 identically on the thin torus,
H6 :  exponents dot (v1,v2,v3) = (4,-2,6),
H8 :  exponents dot (v1,v2,v3) = (5,-1,7).
```

Because `F_5^times` has no zero, a semi-invariant `P` with value `P(x)`
satisfies `P(t.x) = chi(t) P(x)` with `chi(t) != 0`, so "`P(x) = 0 mod 5`" is
an *orbit property*.  Therefore the holdout zero sets are orbit-unions, and
the mod-5 survivor search over the `4^3 = 64` translates of a solution
reduces to the zero pattern of a single representative.  The empirical
uniformity check confirms the statement directly: for each of the three anchor
solutions, the `(H6 == 0, H8 == 0)` pattern is constant across all 64
translates `[COMPUTATION]`.

## 5. The strongest honest conclusion (exact quantifier)

**[THEOREM]** (anchor-section component cleaning).  In the gauge-and-diagonal
chart, let `X` be a `Q_5` point of the 42-equation branch-`11111`
construction variety whose seven thin coordinates congruent to
`(1,1,1,1,1,4,4) mod 5`.  Then `X` is one of `L0, L1, L2` (Hensel
uniqueness, Section 3), and at each such `X` the `3x3x3` holdout
residues at orders `8, 10, 12` are nonzero mod 625 (Section 3.1).
Equivalently:

> **There is no `Q_5` point of branch `11111` over the anchor thin-torus
> section `(1,1,1,1,1,4,4)` that satisfies the order-eight/twelve
> `3x3x3` holdout.**

## 6. Cross-artifact controls and provenance

The standalone test `tests/test_kw_components.py` (run:
`PYTHONPATH=src .venv/bin/python tests/test_kw_components.py`) imports no
producer module.  It independently rebuilds the 42-equation catalogue, the six
primitives, the `F_5` census (count 3), the Jacobians (2, 4, 1), the
Hensel lifts, the order-8 full-weight residue (250 at the anchor) and the
walk-vs-matrix residues at all five orders, the kernel vectors, the unimodular
diagonal minor, `H4 = 0`, the characters, and the per-orbit uniformity; it
then asserts agreement with `results/kac_ward/components.json` and with the
wave-10 records `[0,0,250]`, `det = 2`, `[226,197,98,63,161,216]`
in `results/kac_ward/branch11111_cubic.json`.

Producer run time ~72 s, test run time ~156 s, both `time.process_time()`
gated, RSS < 6 GB; no floating point.

## 7. Limits and open questions `[UNRESOLVED]`

The component theorem is by construction *section-local*: the seven held
coordinates are pinned to the integers `(1,1,1,1,1,4,4)`.  The
following are NOT decided here:

* the `4^7 - 1` other thin-torus sections (the full 16 384-section census);
* whether the order-10/12 symbolic `translation_invariant_trace_expression`
  can be computed directly (currently intractable; the walk expansion stands in);
* the characteristic-zero component structure of the full 42-ideal over the whole
  thin torus, and hence the final `Q`/`Q_5` decision of branch `11111`.

The honest global status of the branch therefore remains `[UNRESOLVED]`; the
present theorem kills the three anchor-section components and supplies the orbit-union
collapse lemma needed by any future complete census.
