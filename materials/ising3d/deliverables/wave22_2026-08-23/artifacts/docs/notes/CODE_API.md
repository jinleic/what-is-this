# Repository code API (stable contract for all contributors/agents)

Interpreter: `/Users/jinleic/jinleic-workspace/math/ising3d/.venv/bin/python`
Package is installed editable; `import ising...` works from anywhere.
Repo root: `/Users/jinleic/jinleic-workspace/math/ising3d`

## Conventions (see `problem_specification.md`)

* `K = beta J`; `v = tanh K`; `x = e^{-2K}`; `u = e^{-4K} = x^2`.
* Bit convention in transfer-matrix state integers: bit `i` = 1 means spin `+1`.
* A periodic direction of length 2 contributes TWO parallel bonds; length 1 contributes none.
* `Z = e^{K n_b} * sum_q c[q] x^q` where `q` counts UNSATISFIED bonds and `sum_q c[q] = 2^N`.
* `Z = 2^N (cosh K)^{n_b} P(v)` where `P(v)` is the even-subgraph generating polynomial.

## `ising.lattices`

```python
from ising.lattices import hyperrect, chain, square, cubic, Lattice
lat = cubic(3, 3, 3, periodic=True)      # or periodic=(True, True, False)
lat.n_sites, lat.n_bonds, lat.bonds, lat.bonds_by_dir, lat.describe()
```

## `ising.exact_enumeration`  (brute force, N <= 30)

```python
from ising.exact_enumeration import joint_dos, dos_bonds, even_subgraph_polynomial
g  = joint_dos(lat, with_field=True)   # shape (n_1+1,...,n_d+1, N+1); last axis = #up spins
gb = dos_bonds(lat)                    # index = number of SATISFIED bonds
P  = even_subgraph_polynomial(lat)     # list[int], P[k] = #even subgraphs with k edges
```

## `ising.transfer_matrix`  (independent implementation, exact int64)

```python
from ising.transfer_matrix import box_broken_bond_poly, torus_broken_bond_poly
c = box_broken_bond_poly((a, b, c), periodic=(False, False), plus_boundary=False)
c = torus_broken_bond_poly((a, b, c))       # cross-section a*b <= 12 for the trace
```
Both return `list[int]` indexed by the number of UNSATISFIED bonds (`x = e^{-2K}` powers).
`plus_boundary=True` adds ghost bonds to a frozen `+1` exterior so every site attains
coordination `2*dim` -- required by the low-temperature finite-lattice method.
Cross-validated against brute force in `tests/test_tm_vs_enumeration.py` (17 lattices, all pass).

## Verified reference data already in the repo

`tests/test_tm_vs_enumeration.py` passes. Known-good exact values you may regression-test against:

* 4x4 square torus, even-subgraph polynomial
  `[1,0,0,0,24,0,128,0,876,0,3584,0,13160,0,28032,0,39462,0,28032,0,13160,0,3584,0,876,0,128,0,24,0,0,0,1]`
* 2x2x2 cubic torus (doubled bonds), even-subgraph polynomial
  `[1,0,12,0,162,0,2012,0,11631,0,30744,0,41948,0,30744,0,11631,0,2012,0,162,0,12,0,1]`
* `P(1) = 2^{n_b - N + 1}` (cycle-space dimension) for any connected lattice -- cheap invariant.

## Benchmarks (NEVER fit to these; they are falsification targets)

* 2D square exact: `K_c = ln(1+sqrt2)/2 = 0.44068679350977147...`
* 3D sc: `K_c = 0.221654626(5)` (Ferrenberg-Xu-Landau 2018, arXiv:1806.03558)
* 3D sc exponents: `nu = 0.629971(4)`, `eta = 0.036298(2)` (conformal bootstrap, Kos et al. 2016)

## Rules for added code

1. Exact integer/rational arithmetic wherever the object is combinatorial. `mpmath` with an
   explicitly set and recorded precision elsewhere.
2. Every new module gets a `tests/test_<name>.py` that runs standalone with the venv python and
   prints `OK`/raises. No pytest requirement.
3. Every numerical claim written into `results/` carries a `provenance` field naming the script
   that produced it.
4. Never fit a constant to a benchmark and then present it as derived.
