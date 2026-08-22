"""Lattice graph construction with explicit, documented boundary-condition conventions.

Conventions are fixed in ``problem_specification.md`` section 1:

* a periodic direction of length 1 contributes NO bonds (self-loops dropped);
* a periodic direction of length 2 contributes TWO parallel bonds per site pair.

These two rules are what make ``Tr T^{L}`` agree with brute-force enumeration for every ``L>=1``.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product


@dataclass(frozen=True)
class Lattice:
    """A finite hyper-rectangular lattice graph.

    Attributes
    ----------
    shape : tuple[int, ...]
        side lengths ``(L_1, ..., L_d)``
    periodic : tuple[bool, ...]
        per-direction boundary condition
    bonds_by_dir : tuple[tuple[tuple[int, int], ...], ...]
        for each direction, the tuple of ``(site_index, site_index)`` pairs.  Sites are indexed
        in row-major (C) order of the coordinate tuple.
    """

    shape: tuple[int, ...]
    periodic: tuple[bool, ...]
    bonds_by_dir: tuple[tuple[tuple[int, int], ...], ...]

    @property
    def dim(self) -> int:
        return len(self.shape)

    @property
    def n_sites(self) -> int:
        n = 1
        for L in self.shape:
            n *= L
        return n

    @property
    def bonds(self) -> tuple[tuple[int, int], ...]:
        out: list[tuple[int, int]] = []
        for bl in self.bonds_by_dir:
            out.extend(bl)
        return tuple(out)

    @property
    def n_bonds(self) -> int:
        return sum(len(b) for b in self.bonds_by_dir)

    def coord(self, index: int) -> tuple[int, ...]:
        out = []
        for L in reversed(self.shape):
            out.append(index % L)
            index //= L
        return tuple(reversed(out))

    def index(self, coord) -> int:
        idx = 0
        for c, L in zip(coord, self.shape):
            idx = idx * L + (c % L)
        return idx

    def degree_sequence(self) -> list[int]:
        deg = [0] * self.n_sites
        for i, j in self.bonds:
            deg[i] += 1
            deg[j] += 1
        return deg

    def describe(self) -> str:
        bc = "".join("P" if p else "O" for p in self.periodic)
        return f"{'x'.join(map(str, self.shape))}[{bc}] N={self.n_sites} B={self.n_bonds}"


def hyperrect(shape, periodic=True) -> Lattice:
    """Build a hyper-rectangular lattice.

    ``periodic`` may be a bool (applied to all directions) or a per-direction sequence.
    """
    shape = tuple(int(s) for s in shape)
    if any(s < 1 for s in shape):
        raise ValueError("side lengths must be >= 1")
    d = len(shape)
    if isinstance(periodic, bool):
        per = (periodic,) * d
    else:
        per = tuple(bool(p) for p in periodic)
        if len(per) != d:
            raise ValueError("periodic must match dimension")

    def index(coord):
        idx = 0
        for c, L in zip(coord, shape):
            idx = idx * L + c
        return idx

    bonds_by_dir: list[tuple[tuple[int, int], ...]] = []
    for direction in range(d):
        L = shape[direction]
        blist: list[tuple[int, int]] = []
        for coord in product(*(range(s) for s in shape)):
            c = coord[direction]
            if c + 1 < L:
                nxt = list(coord)
                nxt[direction] = c + 1
                blist.append((index(coord), index(tuple(nxt))))
            elif per[direction]:
                # wrap-around bond; dropped entirely when L == 1 (self loop)
                if L == 1:
                    continue
                nxt = list(coord)
                nxt[direction] = 0
                blist.append((index(coord), index(tuple(nxt))))
        bonds_by_dir.append(tuple(blist))
    return Lattice(shape=shape, periodic=per, bonds_by_dir=tuple(bonds_by_dir))


def chain(L: int, periodic: bool = True) -> Lattice:
    return hyperrect((L,), periodic)


def square(Lx: int, Ly: int, periodic=True) -> Lattice:
    return hyperrect((Lx, Ly), periodic)


def cubic(Lx: int, Ly: int, Lz: int, periodic=True) -> Lattice:
    return hyperrect((Lx, Ly, Lz), periodic)


def triangular(Lx: int, Ly: int, periodic: bool = True) -> Lattice:
    """Triangular lattice as a square lattice plus one diagonal per cell (control lattice)."""
    base = hyperrect((Lx, Ly), periodic)
    diag: list[tuple[int, int]] = []
    for x in range(Lx):
        for y in range(Ly):
            if periodic or (x + 1 < Lx and y + 1 < Ly):
                if (Lx == 1 and periodic) or (Ly == 1 and periodic):
                    continue
                diag.append((base.index((x, y)), base.index(((x + 1) % Lx, (y + 1) % Ly))))
    return Lattice(
        shape=base.shape,
        periodic=base.periodic,
        bonds_by_dir=base.bonds_by_dir + (tuple(diag),),
    )
