"""Layer-to-layer transfer matrix with EXACT integer polynomial arithmetic.

Independent of ``ising.exact_enumeration``: it never enumerates configurations of the full
lattice, only of a single layer, and it propagates an integer polynomial in the broken-bond
variable

    x  =  e^{-2K}        (each UNSATISFIED bond carries one power of x)

so that for a lattice with ``n_b`` bonds

    Z = e^{K n_b} * sum_q c[q] x^q,      c[q] in Z_{>=0},  sum_q c[q] = 2^N.

Every operation in the propagation is a *shift* or an *addition* of integer arrays; no two
unknown integers are ever multiplied.  Hence ``int64`` is exact provided no coefficient exceeds
2^63-1, and coefficients are bounded by 2^N, so ``N <= 62`` is safe.

Bit convention: bit ``i`` of the state integer is ``1`` for spin ``+1`` and ``0`` for spin ``-1``.

Optional ``plus_boundary`` adds "ghost" bonds so that every site reaches the full lattice
coordination number, the missing neighbours being frozen at ``+1``.  This is what makes
low-temperature cluster weights position-independent, as the finite-lattice method requires.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "box_broken_bond_poly",
    "torus_broken_bond_poly",
    "layer_bonds",
]


# ---------------------------------------------------------------------------
# layer combinatorics
# ---------------------------------------------------------------------------


def layer_bonds(cross: tuple[int, ...], periodic: tuple[bool, ...]):
    """Bonds inside one layer.  ``cross`` is the cross-section shape (1 or 2 entries)."""
    if len(cross) == 1:
        (a,) = cross
        (pa,) = periodic
        bonds = [(i, i + 1) for i in range(a - 1)]
        if pa and a >= 2:
            # a == 2 yields the second parallel bond (1,0): see problem_specification sec.1
            bonds.append((a - 1, 0))
        return bonds
    a, b = cross
    pa, pb = periodic
    bonds = []
    for x in range(a):
        for y in range(b):
            i = x * b + y
            if x + 1 < a:
                bonds.append((i, (x + 1) * b + y))
            elif pa and a >= 2:
                bonds.append((i, y))
            if y + 1 < b:
                bonds.append((i, x * b + y + 1))
            elif pb and b >= 2:
                bonds.append((i, x * b))
    return bonds


def _broken_in_layer(ns: int, bonds) -> np.ndarray:
    s = np.arange(1 << ns, dtype=np.uint32)
    out = np.zeros(1 << ns, dtype=np.int64)
    one = np.uint32(1)
    for i, j in bonds:
        out += (((s >> np.uint32(i)) ^ (s >> np.uint32(j))) & one).astype(np.int64)
    return out


def _weighted_downcount(ns: int, w) -> np.ndarray:
    """``sum_i w[i] * [spin i is DOWN]`` for every layer state."""
    s = np.arange(1 << ns, dtype=np.uint32)
    out = np.zeros(1 << ns, dtype=np.int64)
    one = np.uint32(1)
    for i, wi in enumerate(w):
        if wi:
            out += int(wi) * (1 - ((s >> np.uint32(i)) & one)).astype(np.int64)
    return out


# ---------------------------------------------------------------------------
# exact polynomial propagation
# ---------------------------------------------------------------------------


def _shift_rows(vec: np.ndarray, shifts: np.ndarray) -> np.ndarray:
    """Multiply row ``s`` by ``x^{shifts[s]}`` (degree axis is the LAST axis, truncating)."""
    deg = vec.shape[-1] - 1
    out = np.zeros_like(vec)
    for k in np.unique(shifts):
        k = int(k)
        sel = shifts == k
        if k == 0:
            out[sel] = vec[sel]
        elif k <= deg:
            out[sel, ..., k:] = vec[sel, ..., : deg + 1 - k]
    return out


def _apply_interlayer(vec: np.ndarray, ns: int) -> np.ndarray:
    """Apply ``prod_i ( I + x * X_i )`` where ``X_i`` flips layer spin ``i``.

    This is exactly the inter-layer Boltzmann operator ``prod_i [[1,x],[x,1]]``, in which a
    disagreeing vertical bond carries one factor ``x``.  Degree axis is last; ``vec`` has shape
    ``(2**ns, ..., deg+1)`` with the state index leading.
    """
    tail = vec.shape[1:]
    deg = vec.shape[-1] - 1
    v = vec
    for i in range(ns):
        hi = 1 << (ns - i - 1)
        lo = 1 << i
        v = v.reshape((hi, 2, lo) + tail)
        out = v.copy()
        out[..., 1:] += v[:, ::-1][..., :deg]
        v = out.reshape((1 << ns,) + tail)
    return v


def box_broken_bond_poly(
    shape,
    periodic=None,
    plus_boundary: bool = False,
) -> list[int]:
    """Exact ``sum_{sigma} x^{#unsatisfied bonds}`` for a box, open in the transfer direction.

    ``shape`` is ``(a, c)`` in 2D or ``(a, b, c)`` in 3D; the LAST entry is the transfer
    direction and must be open.  ``periodic`` gives per-direction BCs for the cross-section
    directions (default all open).
    """
    shape = tuple(int(s) for s in shape)
    dim = len(shape)
    if dim not in (2, 3):
        raise ValueError("dim must be 2 or 3")
    cross, c = shape[:-1], shape[-1]
    if periodic is None:
        per = (False,) * (dim - 1)
    else:
        per = tuple(bool(p) for p in periodic)[: dim - 1]
    ns = 1
    for s in cross:
        ns *= s
    if ns > 22:
        raise ValueError(f"cross-section {cross} has {ns} sites: too large")
    ntot = ns * c
    if ntot > 62:
        raise ValueError(f"{ntot} sites: int64 coefficients would overflow")

    lb = layer_bonds(cross, per)
    n_bonds = len(lb) * c + ns * (c - 1)

    inplane_coord = 2 * (dim - 1)
    deg_inplane = [0] * ns
    for i, j in lb:
        deg_inplane[i] += 1
        deg_inplane[j] += 1
    ghost_inplane = [inplane_coord - d for d in deg_inplane]
    assert all(g >= 0 for g in ghost_inplane)
    n_ghost = (sum(ghost_inplane) * c + 2 * ns) if plus_boundary else 0
    deg = n_bonds + n_ghost

    broken_layer = _broken_in_layer(ns, lb)
    if plus_boundary:
        g_in = _weighted_downcount(ns, ghost_inplane)
        g_z = _weighted_downcount(ns, [1] * ns)
    else:
        g_in = np.zeros(1 << ns, dtype=np.int64)
        g_z = np.zeros(1 << ns, dtype=np.int64)

    vec = np.zeros((1 << ns, deg + 1), dtype=np.int64)
    vec[:, 0] = 1
    # The first layer carries its outward z ghosts; so does the last.  When c == 1 the single
    # layer is BOTH the first and the last and must carry 2*g_z -- omitting the doubling
    # undercounts one ghost bond per site for every c == 1 box.
    first = broken_layer + g_in + g_z + (g_z if c == 1 else 0)
    vec = _shift_rows(vec, first)
    for k in range(1, c):
        vec = _apply_interlayer(vec, ns)
        shift = broken_layer + g_in + (g_z if k == c - 1 else 0)
        vec = _shift_rows(vec, shift)
    total = vec.sum(axis=0)
    assert total.min() >= 0, "int64 overflow"
    out = [int(t) for t in total]
    assert sum(out) == 1 << ntot, f"sum {sum(out)} != 2^{ntot}"
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def torus_broken_bond_poly(shape, block: int = 0) -> list[int]:
    """Fully periodic lattice: exact ``Tr T^{L_z}``.

    ``shape = (a, b, c)`` (3D) or ``(a, c)`` (2D).  Cross-section is periodic in all its
    directions.  Memory scales as ``4^{ns} * deg``; keep ``ns <= 12`` unless you have RAM.
    """
    shape = tuple(int(s) for s in shape)
    dim = len(shape)
    cross, c = shape[:-1], shape[-1]
    ns = 1
    for s in cross:
        ns *= s
    ntot = ns * c
    if ntot > 62:
        raise ValueError("int64 overflow risk")
    per = (True,) * (dim - 1)
    lb = layer_bonds(cross, per)
    n_bonds = len(lb) * c + ns * c
    deg = n_bonds
    broken_layer = _broken_in_layer(ns, lb)
    dimH = 1 << ns
    if block <= 0:
        block = dimH
    trace = np.zeros(deg + 1, dtype=np.int64)
    for c0 in range(0, dimH, block):
        c1 = min(c0 + block, dimH)
        vec = np.zeros((dimH, c1 - c0, deg + 1), dtype=np.int64)
        for j in range(c0, c1):
            vec[j, j - c0, 0] = 1
        for _ in range(c):
            vec = _shift_rows(vec, broken_layer)
            vec = _apply_interlayer(vec, ns)
        for j in range(c0, c1):
            trace += vec[j, j - c0]
    assert trace.min() >= 0
    out = [int(t) for t in trace]
    assert sum(out) == 1 << ntot, f"sum {sum(out)} != 2^{ntot}"
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out
