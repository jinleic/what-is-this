"""Regression test for the spectral (free-fermion-in-disguise) obstruction.

Two independent code paths are exercised, neither importing the other's machinery:

1. an INDEPENDENT dense reimplementation of the Gaussian criterion, written here from the
   definition rather than imported from `experiments/e37_*`, checked against the 2D Ising
   strip (must be Gaussian) and the smallest genuine 3D layer (must not be);
2. the EXACT rational certificate of `experiments/e38_*`, re-derived here on the smallest
   case so that the inertia argument is verified independently of the experiment script.

The point of the test is the CONTRAST. A test that only showed 3D failing would pass just as
well if the criterion were broken and rejected everything; the 2D control is what gives it
teeth.
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ising.transfer_matrix import layer_bonds  # noqa: E402

FAILS: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


# --------------------------------------------------------------------------------------
# 1. independent dense implementation of the Gaussian spectral criterion
# --------------------------------------------------------------------------------------
def spectrum(n: int, bonds, K: float) -> np.ndarray:
    """exp(K*A/2) exp(KB) exp(K*A/2), built here from scratch (no experiment imports)."""
    Kstar = -0.5 * np.log(np.tanh(K))
    c, s = np.cosh(0.5 * Kstar), np.sinh(0.5 * Kstar)
    ea = np.array([[1.0]])
    for _ in range(n):
        ea = np.kron(ea, np.array([[c, s], [s, c]]))
    dim = 1 << n
    diag = np.zeros(dim)
    for k in range(dim):
        sp = [1 - 2 * ((k >> (n - 1 - i)) & 1) for i in range(n)]
        diag[k] = sum(sp[i] * sp[j] for i, j in bonds)
    M = ea @ (np.exp(K * diag)[:, None] * ea)
    return np.linalg.eigvalsh(0.5 * (M + M.T))


def is_subset_sum(t: np.ndarray, tol: float):
    """Forced greedy: the smallest not-yet-generated value must itself be a generator."""
    t = np.sort(t)
    gens: list[float] = []
    while True:
        gen = np.zeros(1)
        for g in gens:
            gen = np.concatenate([gen, gen + g])
        gen = np.sort(gen)
        if len(gen) == len(t):
            return gens, float(np.max(np.abs(gen - t))) <= tol
        if len(gen) > len(t) or len(gens) > 20:
            return gens, False
        i = j = 0
        new = None
        while i < len(t):
            if j < len(gen) and abs(t[i] - gen[j]) <= tol:
                i += 1
                j += 1
                continue
            if j < len(gen) and gen[j] < t[i] - tol:
                return gens, False
            new = float(t[i])
            break
        if new is None:
            return gens, False
        gens.append(new)


def verdict(n: int, bonds, K: float):
    lam = spectrum(n, bonds, K)
    t = np.sort(np.log(lam) - np.log(lam).min())
    gens, ok = is_subset_sum(t, 1e-8 * max(1.0, float(t.max())))
    return len(gens), ok


print("1. independent dense reimplementation of the Gaussian criterion")
for n in (4, 6, 8):
    bonds = [(i, i + 1) for i in range(n - 1)]
    g, ok = verdict(n, bonds, 0.3)
    check(f"2D Ising strip n={n} IS spectrally Gaussian", ok and g == n,
          f"{g} generators for {n} modes")

for (nx, ny) in ((2, 3), (2, 4), (3, 3)):
    bonds = list(layer_bonds((nx, ny), (False, False)))
    n = nx * ny
    g, ok = verdict(n, bonds, 0.3)
    check(f"3D Ising layer {nx}x{ny} is NOT spectrally Gaussian", not ok,
          f"greedy stalled after {g} generators")

# the 2x2 layer is C_4: 2-regular, so Theorem DG predicts 2D-like behaviour
b22 = list(layer_bonds((2, 2), (False, False)))
check("2x2 layer is C_4 (2-regular), i.e. not a genuine 3D layer",
      sorted(len([1 for e in b22 if v in e]) for v in range(4)) == [2, 2, 2, 2],
      f"degree sequence from the shared bond builder: "
      f"{sorted(len([1 for e in b22 if v in e]) for v in range(4))}")


# --------------------------------------------------------------------------------------
# 2. exact rational certificate, re-derived here
# --------------------------------------------------------------------------------------
print()
print("2. exact rational certificate (Sylvester inertia, no floating point)")


def build_R(n: int, bonds, t: Fraction):
    q = (1 + t * t) / (2 * t)
    dim = 1 << n
    b = []
    for k in range(dim):
        sp = [1 - 2 * ((k >> (n - 1 - i)) & 1) for i in range(n)]
        b.append(sum(sp[i] * sp[j] for i, j in bonds))
    eps = b[0] % 2
    assert all((x - eps) % 2 == 0 for x in b)
    d = [q ** ((x - eps) // 2) for x in b]
    tp = [t ** h for h in range(n + 1)]
    P = [[tp[bin(k ^ l).count("1")] for l in range(dim)] for k in range(dim)]
    R = [[Fraction(0)] * dim for _ in range(dim)]
    for i in range(dim):
        for j in range(i, dim):
            R[i][j] = R[j][i] = sum(P[i][k] * d[k] * P[j][k] for k in range(dim))
    return R


def count_at(R, sigma: Fraction):
    n = len(R)
    a = [[R[i][j] - (sigma if i == j else 0) for j in range(n)] for i in range(n)]
    neg = 0
    for k in range(n):
        p = a[k][k]
        if p == 0:
            return None
        if p < 0:
            neg += 1
        inv = 1 / p
        row = a[k]
        for i in range(k + 1, n):
            f = a[i][k] * inv
            if f:
                ai = a[i]
                for j in range(k, n):
                    ai[j] -= f * row[j]
    return neg


def isolate(R, index: int, lo: Fraction, hi: Fraction, rel: Fraction):
    while not (lo > 0 and hi - lo <= rel * lo):
        span = hi - lo
        mid, c = None, None
        for k in range(30):
            trial = lo + span / 2 if k == 0 else lo + span / 2 + span * Fraction((-1) ** k,
                                                                                 3 ** (k + 2))
            if lo < trial < hi:
                c = count_at(R, trial)
                if c is not None:
                    mid = trial
                    break
        if mid is None:
            raise RuntimeError("no usable shift")
        if c <= index:
            lo = mid
        else:
            hi = mid
    return lo, hi


def forced_eigenvalue_present(n: int, bonds, t: Fraction, rel_exp: int = 6):
    """Is lambda_1*lambda_2/lambda_0 -- forced for any Gaussian -- actually in the spectrum?"""
    R = build_R(n, bonds, t)
    assert count_at(R, Fraction(0)) == 0, "R must be positive definite"
    tr = sum(R[i][i] for i in range(len(R)))
    rel = Fraction(1, 10 ** rel_exp)
    (l0, h0), (l1, h1), (l2, h2) = [isolate(R, i, Fraction(0), tr + 1, rel) for i in range(3)]
    plo, phi = l1 * l2 / h0, h1 * h2 / l0
    clo, chi = count_at(R, plo), count_at(R, phi)
    assert clo is not None and chi is not None
    return clo != chi, (clo, chi)


t = Fraction(1, 3)
present, counts = forced_eigenvalue_present(4, [(0, 1), (1, 2), (2, 3)], t)
check("CONTROL: 2D Ising strip n=4 forced-value window contains an eigenvalue", present,
      f"exact inertia counts {counts[0]} -> {counts[1]} across the window (consistency check: "
      f"differing counts locate SOME eigenvalue in the window, not exact presence)")

present, counts = forced_eigenvalue_present(6, list(layer_bonds((2, 3), (False, False))), t)
check("3D Ising layer 2x3 is MISSING the forced eigenvalue", not present,
      f"exact inertia counts identical ({counts[0]}) at both outward-rounded window ends, "
      f"so no eigenvalue lies inside")

print()
if FAILS:
    print(f"FAIL: {len(FAILS)} checks failed: {FAILS}")
    raise SystemExit(1)
print("PASS: the 2D control passes the Gaussian consistency checks and the 3D layer is certifiably not Gaussian.")
