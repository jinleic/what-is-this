"""Exact brute-force enumeration of finite Ising lattices with integer arithmetic.

Everything here is exact: the outputs are integer arrays counting spin configurations.
No floating point is used anywhere in the enumeration path.

Central object: the joint density of states

    g[b_1, ..., b_d, k]  =  # { sigma : (# satisfied bonds in direction i) = b_i,
                                        (# up spins) = k }

from which every zero-field or finite-field quantity follows exactly
(see ``problem_specification.md`` section 3).
"""

from __future__ import annotations

import numpy as np

from ..lattices import Lattice

__all__ = [
    "joint_dos",
    "dos_bonds",
    "even_subgraph_polynomial",
    "partition_polynomial_y",
    "z_at",
]

_CHUNK = 1 << 22


def joint_dos(lat: Lattice, with_field: bool = True, dtype=np.int64) -> np.ndarray:
    """Exact joint density of states.

    Returns an array of shape ``(n_1+1, ..., n_d+1, N+1)`` (last axis = number of up spins)
    if ``with_field`` else shape ``(n_1+1, ..., n_d+1)``.

    Raises on int64 overflow risk (N >= 63) and validates the total count.
    """
    N = lat.n_sites
    if N > 30:
        raise ValueError(f"brute force refuses N={N} > 30; use the transfer matrix")
    nbd = [len(b) for b in lat.bonds_by_dir]
    shape = tuple(n + 1 for n in nbd) + ((N + 1,) if with_field else ())
    total_bins = 1
    for s in shape:
        total_bins *= s
    counts = np.zeros(total_bins, dtype=np.int64)

    # strides for the flat key
    strides = []
    acc = 1
    for s in reversed(shape):
        strides.append(acc)
        acc *= s
    strides = list(reversed(strides))

    n_states = 1 << N
    for start in range(0, n_states, _CHUNK):
        stop = min(start + _CHUNK, n_states)
        s = np.arange(start, stop, dtype=np.uint32)
        key = np.zeros(stop - start, dtype=np.int64)
        for d, blist in enumerate(lat.bonds_by_dir):
            sat = np.zeros(stop - start, dtype=np.int32)
            for i, j in blist:
                # satisfied iff bits equal
                sat += (~((s >> np.uint32(i)) ^ (s >> np.uint32(j))) & np.uint32(1)).astype(np.int32)
            key += sat.astype(np.int64) * strides[d]
        if with_field:
            up = np.bitwise_count(s).astype(np.int64)
            key += up * strides[len(nbd)]
        counts += np.bincount(key, minlength=total_bins)

    g = counts.reshape(shape)
    assert int(g.sum()) == n_states, "density of states does not sum to 2^N"
    return g


def dos_bonds(lat: Lattice) -> np.ndarray:
    """Zero-field DOS in the TOTAL number of satisfied bonds; shape ``(n_bonds+1,)``."""
    g = joint_dos(lat, with_field=False)
    nb = lat.n_bonds
    out = np.zeros(nb + 1, dtype=np.int64)
    it = np.nditer(g, flags=["multi_index"])
    for val in it:
        out[sum(it.multi_index)] += int(val)
    assert int(out.sum()) == 1 << lat.n_sites
    return out


# ---------------------------------------------------------------------------
# exact polynomial conversions (python integers, no overflow, no floats)
# ---------------------------------------------------------------------------


def _poly_mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    out[i + j] += ai * bj
    return out


def _poly_pow(a: list[int], n: int) -> list[int]:
    result = [1]
    base = a
    while n:
        if n & 1:
            result = _poly_mul(result, base)
        base = _poly_mul(base, base)
        n >>= 1
    return result


def even_subgraph_polynomial(lat: Lattice) -> list[int]:
    """Exact ``P(v) = sum over even subgraphs E of v^{|E|}`` as an integer coefficient list.

    Derivation (``problem_specification.md`` sec. 3):
        Z = 2^N (cosh K)^{n_b} P(v),  and  Z = sum_b g(b) e^{K(2b-n_b)}
        =>  P(v) = 2^{-N} sum_b g(b) (1+v)^b (1-v)^{n_b-b}
    which is an identity between integer polynomials.
    """
    g = dos_bonds(lat)
    nb = lat.n_bonds
    N = lat.n_sites
    acc = [0] * (nb + 1)
    pos_pows = [[1]]
    neg_pows = [[1]]
    for _ in range(nb):
        pos_pows.append(_poly_mul(pos_pows[-1], [1, 1]))
        neg_pows.append(_poly_mul(neg_pows[-1], [1, -1]))
    for b in range(nb + 1):
        c = int(g[b])
        if not c:
            continue
        term = _poly_mul(pos_pows[b], neg_pows[nb - b])
        for i, t in enumerate(term):
            acc[i] += c * t
    denom = 1 << N
    out = []
    for a in acc:
        assert a % denom == 0, "non-integer even-subgraph coefficient: convention error"
        out.append(a // denom)
    assert out[0] == 1, "P(0) must be 1 (the empty subgraph)"
    # odd-length subgraphs on a bipartite lattice must vanish -- not asserted in general
    return out


def partition_polynomial_y(lat: Lattice) -> tuple[int, list[int]]:
    """Return ``(n_bonds, c)`` with ``Z = e^{-K n_b} * sum_b c[b] y^b``, ``y = e^{2K}``."""
    return lat.n_bonds, [int(x) for x in dos_bonds(lat)]


def z_at(lat: Lattice, K, mp=None):
    """Exact-arithmetic-friendly evaluation of Z at a given K (mpmath or float)."""
    import math

    g = dos_bonds(lat)
    nb = lat.n_bonds
    if mp is None:
        return sum(int(g[b]) * math.exp(K * (2 * b - nb)) for b in range(nb + 1))
    return mp.fsum([mp.mpf(int(g[b])) * mp.e ** (K * (2 * b - nb)) for b in range(nb + 1)])
