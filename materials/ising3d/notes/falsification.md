# Independent falsification of claimed exact 3-D Ising solutions

## Reproduction and provenance

Run, from the repository root,

```text
.venv/bin/python experiments/e06_falsify_claims.py
.venv/bin/python tests/test_falsification.py
```

The experiment uses 80-decimal-digit `mpmath` arithmetic for transcendental quantities and exact
Python integers or `Fraction` values for combinatorial series.  Its machine-readable output is
`results/falsification/falsification.json`; the script records its UTC timestamp and precision.
No coefficient below was copied into the computation from a published series.  The value
`0.221654626(5)` is used only as the requested external falsification target, with `(5)` interpreted
as a one-standard-deviation uncertainty of `5e-9`.

## Conclusions at a glance

| Claim tested | Independent result | Conclusion | Classification |
|---|---:|---|---|
| Z.-D. Zhang critical condition `sinh(2K)sinh(6K)=1` | `K_c = 0.240605912529801723748879456712...` | Disagrees with `0.221654626(5)` by about 3.79 million quoted standard deviations | Published objection reproduced numerically |
| Z.-D. Zhang magnetisation coefficient at `u^6` | `+14`, from exact plus-boundary transfer matrices | Supports the coefficient quoted by Perk, not Zhang's `-18` | Published series objection independently reproduced |
| Degang Zhang anisotropic critical condition | Three cyclic assignments give three different `beta_c` values | Falsified by cubic rotational symmetry alone | Published objection independently reproduced |
| Degang Zhang isotropic critical condition | `beta_c = 0.304688931718003...` | Disagrees with `0.221654626` | Published objection independently reproduced |
| Degang Zhang free-energy integral | First high-temperature failure is at `K^4`: candidate `0.75...`, exact `11/4` | Falsified by an exact finite-lattice series test | Published high-temperature-series objection reproduced; the bounded direct quadrature comparison is additional evidence computed here |

The words “additional evidence” do **not** assert literature priority.  The new work in this
repository is the independent computation and error-controlled comparison, not the observation
that the proposals had previously been criticised.

## 1. Z.-D. Zhang (2007, arXiv:0705.1045)

### 1.1 Claimed critical condition

The script solves

```text
sinh(2K) sinh(6K) = 1
```

without using the closed form as the initial answer.  At 80-digit working precision it obtains

```text
root of equation  = 0.2406059125298017237488794567121842115675921671928302598305090844200819
closed form       = 0.2406059125298017237488794567121842115675921671928302598305090844200819
x_c               = 0.6180339887498948482045868343656381177203091798057628621354486227052605
```

Here the closed form is

```text
K_c = -1/2 log((sqrt(5)-1)/2).
```

The independently solved root and the closed form differ by
`9.2748757996887e-82`; the equation residual is `-1.0542197943231e-81`.  Thus the algebraic
interpretation of the claimed condition is confirmed to far more than the requested 30 digits.
It is the condition itself, not its numerical solution, that fails the benchmark test:

```text
absolute discrepancy          = 0.0189512865298017237488794567121842...
relative discrepancy          = 0.0854991699104070300291384701900343...
relative discrepancy (percent)= 8.549916991040703... %
quoted standard deviations    = 3790257.3059603447497758913424...
```

### 1.2 Exact low-temperature magnetisation through `u^6`

The calculation adds a down-spin axis to the repository's exact broken-bond transfer matrix.  For
each plus-boundary box it constructs the exact integer table

```text
c_B[q,n] = number of configurations with q broken bonds and n down spins,
```

and asserts that `sum_n c_B[q,n]` is exactly
`ising.transfer_matrix.box_broken_bond_poly(B, plus_boundary=True)`.  Hence

```text
Z_B^+(x) = sum_(q,n) c_B[q,n] x^q,
rho_B(x) = (sum_(q,n) n c_B[q,n] x^q) / Z_B^+(x)
```

are formed using exact series division.  `rho_B` is the field derivative of `log Z_B^+`, so the
mixed backward difference in all three side lengths removes face, edge, and corner translation
multiplicities.  The code takes that difference at `(3,3,3)`, using all eight boxes whose side
lengths are independently 2 or 3.

Why this is exact through `x^12`: a single flipped spin has surface 6; a connected adjacent pair
has surface 10; any connected set of at least three simple-cubic sites has surface at least 14.
At total degree 12 the only additional connected Mayer contribution consists of two incompatible
single-spin excitations.  Every contribution through degree 12 therefore has bounding-box side at
most 2 and is included by the side-3 difference.  This gives

```text
rho_down(x) = x^6 + 6 x^10 - 7 x^12 + O(x^14),
m(x)        = 1 - 2 x^6 - 12 x^10 + 14 x^12 + O(x^14).
```

Since `u=x^2=e^{-4K}`, the independently computed exact result is

```text
m(u) = 1 - 2 u^3 - 12 u^5 + 14 u^6 + O(u^7).
```

The `u^6` coefficient is therefore `+14`.  It agrees with the coefficient attributed in the task
to Perk and contradicts Zhang's `-18`; the disagreement at this coefficient is 32.

## 2. Degang Zhang (arXiv:2110.11233)

### 2.1 Rotational-invariance test of the critical line

For each cyclic assignment of the same physical couplings `(0.5,1,1.5)`, the script solves the
claimed equation

```text
sinh(2 beta J) sinh(2 beta (J1+J2)) = 1.
```

The results are:

| `J` singled out by the formula | `J1` | `J2` | computed `beta_c` |
|---:|---:|---:|---:|
| 0.5 | 1.0 | 1.5 | 0.3503982204355165621737533493008828504... |
| 1.0 | 1.5 | 0.5 | 0.3046889317180031157684016855841993477... |
| 1.5 | 0.5 | 1.0 | 0.2937911956731810084108697749932641030... |

Their range is

```text
0.0566070247623355537628835743076187474...
```

A rotation of the simple-cubic lattice merely permutes the three coupling labels and cannot change
the physical critical inverse temperature.  The unequal results are therefore a complete
falsification of the claimed anisotropic critical line, independent of Monte Carlo data or any
series coefficient.  Perk's comment, arXiv:2202.03136, prints the same objection and rounded
values; the table above was nevertheless regenerated here directly from the equation.

In the isotropic case the same formula reduces to

```text
sinh(2H) sinh(4H) = 1,
beta_c = 0.3046889317180031157684016855841993477....
```

Relative to `0.221654626`, its absolute and relative discrepancies are
`0.0830343057180031157684...` and `0.374611201292966092972...` (about `37.4611%`).

### 2.2 Literal numerical evaluation of the claimed free energy

The code implements the stated `D_m(omega)`, `xi_m=arccosh(D_m)`, and prefactor literally, with
`H*=atanh(exp(-2H))`.  It divides `[0,pi]` into `4m` panels so the oscillatory terms are resolved.
For isotropic `H=H1=H2=K=0.001`, the finite-`m` trend is:

| `m` | claimed `phi_m(K)` |
|---:|---:|
| 2 | 0.6931486815606973096612965300368529984... |
| 4 | 0.6931486805606963101568915069077942665... |
| 8 | 0.6931486805606953101505581354654832659... |
| 16 | 0.6931486805606953101505581344654562660... |
| 32 | 0.6931486805606953101505581344654562660... |
| 64 | 0.6931486805606953101505581344654562660... |

The 32- and 64-mode values agree at the recorded 70 decimal digits.

#### Finite lattice and its error bound

A raw `L^3` torus is **not** an infinite-lattice reference through order `<2L`: it contains a
single winding even subgraph already at order `v^L`.  The experiment therefore evaluates all
eight periodic/antiperiodic seam twists and averages them.  This exactly projects onto zero
`Z_2` winding.  After projection, a nonlocal even graph must contain at least two winding
traversals, so its length is at least `2L`.  Thus the projected `L^3` coefficients are exact for
orders `n<2L`; for `L=3`, the largest certified order is 5.  This projection is essential to the
`<2L` statement.

The exact projected polynomial begins

```text
P_0(v) = 1 + 81 v^4 + 702 v^6 + ... .
```

The `v^6` term is allowed to have finite-size contamination, as expected.  Independent open-box
finite-lattice inclusion/exclusion gives the infinite-lattice bulk series

```text
(1/N) log P(v)
  = 3 v^4 + 22 v^6 + (375/2) v^8 + 1980 v^10 + 24044 v^12 + O(v^14).
```

The bounding-box proof used in the code is explicit: a connected Eulerian graph with bounding box
`(a,b,c)` needs at least `2[(a-1)+(b-1)+(c-1)]` edges, so every required box through order 12 has
at most 27 sites and is handled exactly.  After `v=tanh K`, the corresponding reduced-free-energy
coefficients are

```text
phi(K) = log(2)
       + (3/2) K^2
       + (11/4) K^4
       + (271/15) K^6
       + (123547/840) K^8
       + (7236721/4725) K^10
       + (5668303511/311850) K^12
       + O(K^14).
```

The remainder used for the numerical bound is deliberately conservative.  At infinite
temperature the bond variables have magnitude one and dependency-graph degree at most 10.  The
connected spanning-tree cumulant estimate used by the script is

```text
|[K^r] phi| <= 3 * 2^(r-1) * r^(r-2) * 11^(r-1) / r! < 64^r.
```

Consequently the absolute tail after order 12 is bounded by
`(64|K|)^13/(1-64|K|)` when `64|K|<1`.  Combining that thermodynamic tail with the exactly computed
difference between the projected `3^3` value and the order-12 bulk polynomial gives, at `K=0.001`,

```text
projected finite-lattice phi       = 0.6931486805626953314841158702642197262...
claimed m->infinity trend          = 0.6931486805606953101505581344654562660...
absolute discrepancy               = 2.0000213335577357987634602083941e-12
finite-size plus series-tail bound = 3.2689692353039426306921924092299e-16
discrepancy / bound                = 6118.2017620663767577...
```

Subtracting the finite-size/tail bound from the stabilised `m=64` discrepancy leaves
`1.9996944366342054e-12`.  This is not presented as a rigorous error bound for the separate
`m -> infinity` extrapolation: it is the result of the observed trend, whose `m=32` and `m=64`
values agree at the recorded precision.  The symmetry and exact Taylor-series falsifications do
not rely on that extrapolation.

## 3. Reusable general free-energy falsifier

`experiments/e06_falsify_claims.py` exports

```python
falsify_free_energy(phi_callable, name)
```

for any real candidate analytic from the right at `K=0`.  It evaluates the callable on two
80-digit grids, performs degree-8 polynomial interpolation, uses the coarse/fine difference as a
numerical stability diagnostic, and compares orders 0 through 5 with the exact zero-winding
`3x3x3` torus reference.  A raw torus would fail at order 3; the eight-twist projection makes the
reference exact through order 5, because the leading zero-winding wrap contribution has order 6.

Applied to the Degang Zhang integral (using `m=8`, whose mode-dependent corrections enter above the
tested orders), it obtains

```text
order       candidate coefficient                 exact coefficient
K^0         log(2)                                log(2)
K^1         -1.29e-25 (numerical zero)            0
K^2         1.50000000000000000000146...          3/2
K^3         -6.63e-18 (numerical zero)            0
K^4         0.750000000000016294925...            11/4 = 2.75
```

The first failure is therefore order `K^4`, with a coefficient discrepancy of 2.  This independently
reproduces the published high-temperature-series objection.  The direct `m=2,...,64` quadrature,
the exact twist-projected finite-lattice value, and the explicit error bound above are separate
numerical corroboration produced in this repository; they are not presented as a claim of new
literature priority.

## 4. Separation of reproduced objections and repository findings

### Published objections reproduced independently

1. The 2007 critical condition yields the wrong isotropic critical point.
2. The 2007 low-temperature magnetisation fails at `u^6`; exact finite-lattice arithmetic gives
   `+14`, not `-18`.
3. The 2021 anisotropic critical line changes under a rotation of the cubic lattice; the three
   values printed above reproduce the objection in Perk, arXiv:2202.03136.
4. The 2021 free energy fails the exact high-temperature series, here first at `K^4`.

### Additional computations performed here

1. A new exact transfer-matrix down-spin axis and finite-lattice difference derive the
   magnetisation coefficient rather than assuming it.
2. An eight-twist zero-winding projection states and enforces the correct finite-torus exactness
   bound: raw torus `<L`, projected torus `<2L`.
3. Literal oscillatory quadrature of the 2021 formula is compared with an exact projected finite
   lattice at `K=0.001`, with a conservative analytic tail bound over 6,000 times smaller than the
   observed discrepancy.
4. The reusable falsifier determines the first failing Taylor order of future candidate free
   energies without fitting any coefficient to the critical benchmark.
