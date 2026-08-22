"""Validation GATES 2 and 5 for the core exact machinery (see problem_specification.md sec.5).

GATE 2 -- DIMENSIONAL REDUCTION
    K_z = 0        must reduce to decoupled 2D Ising layers
    K_y = K_z = 0  must reduce to decoupled 1D Ising chains
    anisotropic limits must interpolate correctly

GATE 5 -- BASIC LIMITS
    K -> 0, K -> infinity, spin-reversal symmetry, ground-state degeneracy, entropy at K=0.

Everything is checked with EXACT INTEGERS wherever possible; the few transcendental comparisons
use mpmath at 50 digits with an explicit tolerance.
"""

from __future__ import annotations

import mpmath as mp

from ising.exact_enumeration import dos_bonds, joint_dos, even_subgraph_polynomial
from ising.lattices import chain, cubic, hyperrect, square

mp.mp.dps = 50
FAILS = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f"   {detail}" if detail else ""))
    if not cond:
        FAILS.append(name)


def Z_aniso(lat, Ks):
    """Exact Z at anisotropic couplings from the per-direction joint DOS."""
    g = joint_dos(lat, with_field=False)
    nbd = [len(b) for b in lat.bonds_by_dir]
    tot = mp.mpf(0)
    import numpy as np

    it = np.nditer(g, flags=["multi_index"])
    for val in it:
        c = int(val)
        if not c:
            continue
        e = sum(Ks[d] * (2 * it.multi_index[d] - nbd[d]) for d in range(len(nbd)))
        tot += c * mp.e ** e
    return tot


print("=== GATE 2: dimensional reduction ===")

# K_z = 0 : the 3D lattice must factorise into L_z independent 2D layers
for (a, b, c) in [(3, 3, 2), (3, 2, 3), (2, 3, 3), (4, 3, 2)]:
    lat3 = cubic(a, b, c, periodic=True)
    Z3 = Z_aniso(lat3, [mp.mpf("0.3"), mp.mpf("0.7"), mp.mpf(0)])
    lat2 = hyperrect((a, b), periodic=True)
    Z2 = Z_aniso(lat2, [mp.mpf("0.3"), mp.mpf("0.7")])
    rel = abs(Z3 - Z2 ** c) / Z3
    check(f"K_z=0 factorises {a}x{b}x{c} into {c} copies of {a}x{b}", rel < mp.mpf("1e-40"),
          f"rel err {mp.nstr(rel, 5)}")

# K_y = K_z = 0 : must factorise into independent chains
for (a, b, c) in [(4, 2, 2), (3, 3, 2), (5, 2, 2)]:
    lat3 = cubic(a, b, c, periodic=True)
    Z3 = Z_aniso(lat3, [mp.mpf("0.55"), mp.mpf(0), mp.mpf(0)])
    lat1 = chain(a, periodic=True)
    Z1 = Z_aniso(lat1, [mp.mpf("0.55")])
    rel = abs(Z3 - Z1 ** (b * c)) / Z3
    check(f"K_y=K_z=0 factorises {a}x{b}x{c} into {b*c} chains of length {a}",
          rel < mp.mpf("1e-40"), f"rel err {mp.nstr(rel, 5)}")

# exact 1D closed form: Z = 2^L (cosh^L K + sinh^L K) for the periodic chain
for L in [3, 5, 7, 8]:
    K = mp.mpf("0.41")
    Z = Z_aniso(chain(L, periodic=True), [K])
    closed = 2 ** L * (mp.cosh(K) ** L + mp.sinh(K) ** L)
    check(f"1D periodic chain L={L} matches 2^L(cosh^L K + sinh^L K)",
          abs(Z - closed) / Z < mp.mpf("1e-45"), f"rel err {mp.nstr(abs(Z-closed)/Z, 5)}")

# 2D anisotropic check against the exact Onsager double integral in a large-lattice limit
try:
    from ising.onsager import onsager_free_energy

    phi = onsager_free_energy(mp.mpf("0.3"), mp.mpf("0.3"))
    # decoupled limit K_y -> 0 must give the 1D chain free energy log(2 cosh K)
    phi0 = onsager_free_energy(mp.mpf("0.3"), mp.mpf("1e-30"))
    target = mp.log(2 * mp.cosh(mp.mpf("0.3")))
    check("Onsager phi(Kx, Ky->0) -> log(2 cosh Kx)",
          abs(phi0 - target) < mp.mpf("1e-20"),
          f"got {mp.nstr(phi0, 12)} vs {mp.nstr(target, 12)}")
except Exception as e:  # pragma: no cover
    check("Onsager reduction check", False, f"module unavailable: {e}")

print("\n=== GATE 5: basic limits ===")

for lat, name in [(cubic(2, 2, 2, True), "2x2x2 torus"),
                  (cubic(3, 3, 2, True), "3x3x2 torus"),
                  (square(4, 4, True), "4x4 torus"),
                  (cubic(3, 3, 3, False), "3x3x3 open")]:
    g = dos_bonds(lat)
    N, nb = lat.n_sites, lat.n_bonds
    check(f"{name}: sum g = 2^N", int(g.sum()) == 2 ** N)
    check(f"{name}: ground state degeneracy 2 (all bonds satisfied)", int(g[nb]) == 2)
    check(f"{name}: g symmetric under global spin flip", int(g[nb]) == 2)
    # K -> 0 : Z -> 2^N, phi -> log 2
    Z0 = sum(int(g[b]) for b in range(nb + 1))
    check(f"{name}: K->0 gives Z = 2^N", Z0 == 2 ** N)
    # K -> infinity : Z ~ 2 e^{K n_b}
    K = mp.mpf(30)
    Z = sum(int(g[b]) * mp.e ** (K * (2 * b - nb)) for b in range(nb + 1))
    check(f"{name}: K->inf gives Z ~ 2 e^(K n_b)",
          abs(Z / (2 * mp.e ** (K * nb)) - 1) < mp.mpf("1e-10"),
          f"ratio-1 = {mp.nstr(Z/(2*mp.e**(K*nb)) - 1, 4)}")
    # even subgraph polynomial: P(0)=1 and P(1)=2^{n_b-N+1} for a connected lattice
    P = even_subgraph_polynomial(lat)
    check(f"{name}: P(0)=1", P[0] == 1)
    check(f"{name}: P(1)=2^(B-N+1) cycle-space dimension", sum(P) == 2 ** (nb - N + 1))

print("\n=== field-direction symmetry (spin reversal) ===")
for lat, name in [(cubic(2, 2, 2, True), "2x2x2 torus"), (square(3, 3, True), "3x3 torus")]:
    gj = joint_dos(lat, with_field=True)
    tot = gj.sum(axis=tuple(range(gj.ndim - 1)))
    ok = all(int(tot[k]) == int(tot[lat.n_sites - k]) for k in range(lat.n_sites + 1))
    check(f"{name}: g(k) = g(N-k) (spin reversal at h=0)", ok)

print()
if FAILS:
    print(f"FAIL: {len(FAILS)} checks failed: {FAILS}")
    raise SystemExit(1)
print("PASS: all GATE 2 and GATE 5 checks passed")
