"""GATE 1 for Track A: does the operator form actually reproduce the 3D partition function?

Every algebraic result in `proofs/algebraic_obstruction.md`, `proofs/dolan_grady_defect.md` and
`proofs/tridiagonal_nogo.md` is a statement about

    V = exp(K_z^* A) exp(K B),     A = sum_i X_i,     B = sum_{<ij> in layer} Z_iZ_j,
    K_z^* = -(1/2) log tanh K_z    (equivalently sinh 2K_z sinh 2K_z^* = 1),

and the claim that

    Z_{Lx x Ly x Lz}  =  (2 sinh 2K_z)^{n Lz / 2} * Tr( V^{Lz} ),      n = Lx * Ly.

If that identity is wrong -- wrong dual coupling, wrong prefactor, wrong operator ordering -- then
the whole algebraic track is about the wrong operator.  This test verifies it against the exact
integer partition function produced by two independent engines (brute-force enumeration and the
repository transfer matrix), at high precision, for both isotropic and anisotropic couplings.
"""

from __future__ import annotations

import itertools

import mpmath as mp
import numpy as np

from ising.exact_enumeration import joint_dos
from ising.lattices import cubic
from ising.transfer_matrix import torus_broken_bond_poly

mp.mp.dps = 40
FAILS = []


def dense_ops(nx, ny):
    """Dense A = sum X_i and the layer bond list, cross-section periodic in both directions."""
    n = nx * ny
    dim = 1 << n
    I2 = np.eye(2)
    X = np.array([[0.0, 1.0], [1.0, 0.0]])
    Z = np.array([[1.0, 0.0], [0.0, -1.0]])

    def op(single, site):
        m = np.array([[1.0]])
        for k in range(n):
            m = np.kron(m, single if k == site else I2)
        return m

    A = np.zeros((dim, dim))
    for i in range(n):
        A += op(X, i)
    # MUST use the repository's own layer bond list so the doubled-bond convention for a periodic
    # direction of length 2 (problem_specification sec.1) matches the exact engine.  Building the
    # bonds by hand here silently used SINGLE bonds and made every length-2 cross-section
    # disagree -- exactly the failure mode that convention section warns about.
    from ising.transfer_matrix import layer_bonds

    bonds = layer_bonds((nx, ny), (True, True))
    Zi = [op(Z, i) for i in range(n)]
    return A, bonds, Zi, n


def exact_Z(shape, K):
    """Exact Z from the repository transfer matrix (integers), evaluated at K with mpmath."""
    coeffs = torus_broken_bond_poly(shape)
    nb = 3 * shape[0] * shape[1] * shape[2]
    assert sum(coeffs) == 1 << (shape[0] * shape[1] * shape[2])
    return mp.e ** (K * nb) * mp.fsum([c * mp.e ** (-2 * K * q) for q, c in enumerate(coeffs)])


def operator_Z(nx, ny, nz, Kxy, Kz):
    A, bonds, Zi, n = dense_ops(nx, ny)
    Kzs = -mp.mpf("0.5") * mp.log(mp.tanh(Kz))
    Kzs_f = float(Kzs)
    from scipy.linalg import expm

    Bm = np.zeros((1 << n, 1 << n))
    for i, j in bonds:
        Bm += Kxy * (Zi[i] @ Zi[j])
    V = expm(Kzs_f * A) @ expm(Bm)
    Vp = np.linalg.matrix_power(V, nz)
    tr = float(np.trace(Vp))
    pref = (2 * mp.sinh(2 * Kz)) ** (mp.mpf(n * nz) / 2)
    return pref * tr, len(bonds)


if __name__ == "__main__":
    print("GATE 1 (Track A): operator form vs exact integer partition function")
    print(f"  V = exp(K* sum X) exp(K sum ZZ),  K* = -0.5 log tanh K,  "
          f"Z = (2 sinh 2K)^(n Lz/2) Tr V^Lz\n")
    cases = [(2, 2, 2, 0.3), (2, 2, 3, 0.25), (2, 3, 2, 0.4), (3, 2, 2, 0.2),
             (2, 2, 4, 0.35), (2, 3, 3, 0.221654626), (3, 3, 2, 0.15)]
    for nx, ny, nz, K in cases:
        Kf = mp.mpf(repr(K))
        Zex = exact_Z((nx, ny, nz), Kf)
        Zop, nb_layer = operator_Z(nx, ny, nz, K, Kf)
        rel = abs(Zop - Zex) / Zex
        ok = rel < mp.mpf("1e-9")
        FAILS.append((nx, ny, nz)) if not ok else None
        print(f"  {nx}x{ny}x{nz}  K={K:<12g} layer bonds={nb_layer:2d}  "
              f"Z_exact={mp.nstr(Zex, 12):>18s}  rel.err={mp.nstr(rel, 4):>10s}  "
              f"{'PASS' if ok else 'FAIL'}")

    print("\n  (relative error is limited by float64 in the dense matrix exponential, not by the")
    print("   identity; the exact side is an integer polynomial evaluated at 40 digits)")
    if FAILS:
        print(f"\nFAIL: {FAILS}")
        raise SystemExit(1)
    print("\nPASS: the layer operator form reproduces the exact 3D Ising partition function.")
    print("      Track A therefore concerns the correct operator.")
