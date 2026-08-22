#!/usr/bin/env python3
"""All-radius structure of translation-covariant extensive charges of H = f*sum X + g*sum ZZ.

Exact finite evidence for `proofs/allr_charges.md`.  This is an independent
implementation: it does not import `e117_extensive_charges.py`, and the global
box quotients it recomputes are compared against the frozen
`results/integrability/extensive_charges.json` only as a consistency report.

Conventions (the repository's ordered-Pauli convention):

    word  w = (a, b)  ->  X^a Z^b ,   a, b finite subsets of Z^d
    ad_H  = (1/2)[H, .]  with  H = f * sum_x X_x + g * sum_{<xy>} Z_x Z_y

    ad_H (a,b) = f * sum_{p in b} (a XOR p, b)
               - g * sum_{p in a} sum_{q ~ p, q not in a} (a, b XOR p XOR q)

`pi` sums Pauli coefficients over translation orbits.  By the divergence
criterion of `proofs/extensive_charges.md`, the formal translation sum
Q(q) = sum_x tau_x(q) is conserved iff pi(ad_H q) = 0.  Everything below
computes with  Lbar = pi ad_H  on orbit space.

Blocks
------
A. Global box quotients  dim ker(Lbar|V_R) - 1, exactly, plus exact charge bases.
B. The direction-i leading symbol  pi Gamma^(i)  on a single shape sigma:
   (i)  every orbit row has at most two nonzero entries  [structure theorem],
   (ii) exact kernel by an O(rows) grounding/propagation solver,
   (iii) the Theorem M predicates T1 (type i) and T2 (thick a-free face),
   (iv) the combined all-direction leading symbol.
C. Single-shape charge kernels ker(Lbar|V_sigma), obtained exactly by first
   restricting to the (proved) combined leading-symbol kernel.
D. The mandatory d = 1 control: the same code must find nonzero kernels in d = 1
   and must certify that the T2 hypothesis is unsatisfiable in d = 1.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import platform
import resource
import sys
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "allr_charges.json"
E117_RESULT = ROOT / "results" / "integrability" / "extensive_charges.json"
SCRIPT = "experiments/e127_allr_charges.py"

RSS_WALL_BYTES = 7 * 1024**3

# (field, bond) coefficient pairs; both nonzero is the nondegenerate regime.
COEFFICIENTS = ((1, 1), (2, 3), (1, -1))
# Global boxes: (dimension, radius).  B_R has (R+1)^d sites.
GLOBAL_CASES = (
    (1, 0), (1, 1), (1, 2), (1, 3), (1, 4),
    (2, 0), (2, 1),
    (3, 0),
)
GLOBAL_KERNEL_CASES = ((1, 1), (1, 2), (1, 3), (2, 1))

# Leading-symbol shapes (4^volume enumeration + linear solving).
SHAPES_SYMBOL = (
    (2,), (3,), (4,), (5,), (6,), (7,), (8,), (9,),
    (2, 2), (3, 2), (2, 3), (4, 2), (2, 4), (3, 3),
    (2, 2, 2),
)
# Single-shape charge-kernel shapes.
SHAPES_FULL = (
    (2,), (3,), (4,), (5,), (6,), (7,),
    (2, 2), (3, 2), (2, 3), (4, 2), (2, 4), (3, 3),
    (2, 2, 2),
)
# Shapes where the graph solver is cross-validated against Gaussian elimination
# and where the leading-symbol identification is verified term by term.
SHAPES_CROSSCHECK = ((2,), (3,), (4,), (2, 2), (3, 2), (2, 3), (2, 2, 2))
# Shapes for the LOCAL pointwise row verification.  These are far beyond the
# 4^volume enumeration frontier, so Theorem M is tested where no census exists.
SHAPES_LOCAL = (
    (12,), (20,),
    (2, 2), (6, 4), (8, 5), (12, 7),
    (2, 2, 2), (3, 3, 3), (4, 3, 3), (5, 4, 3), (6, 5, 4),
)
LOCAL_SAMPLES = 200

GLOBAL_WALL_SECONDS = 900.0
SYMBOL_WALL_SECONDS = 900.0
FULL_WALL_SECONDS = 900.0


# --------------------------------------------------------------------------- #
# resource bookkeeping                                                        #
# --------------------------------------------------------------------------- #
def peak_rss_bytes() -> int:
    """Darwin reports ru_maxrss in bytes; Linux reports KiB."""
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


class ResourceWall(RuntimeError):
    def __init__(self, payload: dict[str, object]) -> None:
        super().__init__(str(payload.get("reason")))
        self.payload = payload


class WallClock:
    def __init__(self, phase: str, seconds: float) -> None:
        self.phase = phase
        self.seconds = seconds
        self.started = time.monotonic()

    def check(self, processed: int = 0) -> None:
        elapsed = time.monotonic() - self.started
        rss = peak_rss_bytes()
        if elapsed > self.seconds:
            raise ResourceWall(
                {
                    "claim_tag": "[UNRESOLVED]",
                    "phase": self.phase,
                    "reason": f"wall time {self.seconds:g} s reached",
                    "elapsed_seconds": round(elapsed, 6),
                    "processed": processed,
                    "peak_rss_mib": round(rss / 1024**2, 3),
                }
            )
        if rss > RSS_WALL_BYTES:
            raise ResourceWall(
                {
                    "claim_tag": "[UNRESOLVED]",
                    "phase": self.phase,
                    "reason": f"RSS wall {RSS_WALL_BYTES // 1024**3} GiB reached",
                    "elapsed_seconds": round(elapsed, 6),
                    "processed": processed,
                    "peak_rss_mib": round(rss / 1024**2, 3),
                }
            )

    def elapsed(self) -> float:
        return time.monotonic() - self.started


# --------------------------------------------------------------------------- #
# geometry                                                                    #
# --------------------------------------------------------------------------- #
def bits_of(mask: int):
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low


class Box:
    """Inner box of side lengths `shape`, offset by one inside a one-shell box.

    Bit `i` refers to `self.coords[i]`, a coordinate of the ENLARGED box; the
    inner sites are those with every coordinate in [1, shape[k]].
    """

    __slots__ = (
        "shape", "dimension", "ext_shape", "coords", "index",
        "inner_bits", "neighbors", "volume", "expand", "slice_masks",
        "growth",
    )

    def __init__(self, shape: tuple[int, ...]) -> None:
        self.shape = shape
        self.dimension = len(shape)
        self.ext_shape = tuple(s + 2 for s in shape)
        self.coords = tuple(itertools.product(*(range(s) for s in self.ext_shape)))
        self.index = {c: i for i, c in enumerate(self.coords)}
        self.inner_bits = tuple(
            self.index[c] for c in itertools.product(*(range(1, s + 1) for s in shape))
        )
        self.volume = len(self.inner_bits)
        neighbors: dict[int, tuple[int, ...]] = {}
        for bit in self.inner_bits:
            site = self.coords[bit]
            entries: list[int] = []
            for axis in range(self.dimension):
                for step in (-1, 1):
                    nb = list(site)
                    nb[axis] += step
                    entries.append(self.index[tuple(nb)])
            neighbors[bit] = tuple(entries)
        self.neighbors = neighbors
        # (a) subset -> inner-bit mask, for fast enumeration; built on demand
        #     because it has 2^volume entries.
        self.expand = None
        # (b) extreme-slice masks, for the O(1) spanning test
        slice_masks = []
        for axis in range(self.dimension):
            lo = sum(
                1 << bit for bit in self.inner_bits if self.coords[bit][axis] == 1
            )
            hi = sum(
                1 << bit
                for bit in self.inner_bits
                if self.coords[bit][axis] == shape[axis]
            )
            slice_masks.append((lo, hi))
        self.slice_masks = tuple(slice_masks)
        # (c) growth tables: for each (axis, sign) the target shape and the
        #     already-normalised target index of every enlarged-box bit.  A
        #     direction-`axis` growth term lands in shape sigma + e_axis with
        #     the origin shifted by one in that axis when the growth is at the
        #     low face; both cases are therefore pure re-indexing, with no
        #     bounding-box scan.
        growth: dict[tuple[int, int], tuple] = {}
        for axis in range(self.dimension):
            shape2 = list(shape)
            shape2[axis] += 1
            shape2 = tuple(shape2)
            for sign in (1, -1):
                offs = [1] * self.dimension
                if sign < 0:
                    offs[axis] = 0
                table = [-1] * len(self.coords)
                for bit, site in enumerate(self.coords):
                    rel = [site[k] - offs[k] for k in range(self.dimension)]
                    if any(r < 0 or r >= shape2[k] for k, r in enumerate(rel)):
                        continue
                    idx = 0
                    for k in range(self.dimension):
                        idx = idx * shape2[k] + rel[k]
                    table[bit] = idx
                growth[(axis, sign)] = (shape2, tuple(table))
        self.growth = growth

    def encode_growth(self, mask: int, table) -> int:
        out = 0
        for bit in bits_of(mask):
            out |= 1 << table[bit]
        return out

    def canon(self, a: int, b: int) -> tuple:
        support = a | b
        if not support:
            return ((), 0, 0)
        d = self.dimension
        mins = [10**9] * d
        maxs = [-(10**9)] * d
        for i in bits_of(support):
            site = self.coords[i]
            for k in range(d):
                c = site[k]
                if c < mins[k]:
                    mins[k] = c
                if c > maxs[k]:
                    maxs[k] = c
        shape2 = tuple(maxs[k] - mins[k] + 1 for k in range(d))
        a2 = 0
        b2 = 0
        for i in bits_of(a):
            site = self.coords[i]
            idx = 0
            for k in range(d):
                idx = idx * shape2[k] + site[k] - mins[k]
            a2 |= 1 << idx
        for i in bits_of(b):
            site = self.coords[i]
            idx = 0
            for k in range(d):
                idx = idx * shape2[k] + site[k] - mins[k]
            b2 |= 1 << idx
        return (shape2, a2, b2)

    def decode(self, a_code: int, b_code: int) -> tuple[int, int]:
        a = 0
        b = 0
        for idx, bit in enumerate(self.inner_bits):
            if (a_code >> idx) & 1:
                a |= 1 << bit
            if (b_code >> idx) & 1:
                b |= 1 << bit
        return a, b

    def ad_terms(self, a: int, b: int, field: int, bond: int) -> dict[tuple[int, int], int]:
        out: dict[tuple[int, int], int] = {}
        for p in bits_of(b):
            key = (a ^ (1 << p), b)
            v = out.get(key, 0) + field
            if v:
                out[key] = v
            else:
                out.pop(key, None)
        for p in bits_of(a):
            pbit = 1 << p
            for q in self.neighbors[p]:
                qbit = 1 << q
                if a & qbit:
                    continue
                key = (a, b ^ pbit ^ qbit)
                v = out.get(key, 0) - bond
                if v:
                    out[key] = v
                else:
                    out.pop(key, None)
        return out

    def gamma_keys(self, a: int, b: int, axis: int, bond: int) -> dict[tuple, int]:
        """Direction-`axis` box-GROWING bond terms, already in canonical form."""
        out: dict[tuple, int] = {}
        lo, hi = 1, self.shape[axis]
        cached: dict[int, tuple] = {}
        for p in bits_of(a):
            site = self.coords[p]
            pbit = 1 << p
            for edge, sign in ((hi, 1), (lo, -1)):
                if site[axis] != edge:
                    continue
                nb = list(site)
                nb[axis] += sign
                qbit = 1 << self.index[tuple(nb)]
                shape2, table = self.growth[(axis, sign)]
                enc_a = cached.get(sign)
                if enc_a is None:
                    enc_a = cached[sign] = (shape2, self.encode_growth(a, table))
                key = (
                    enc_a[0],
                    enc_a[1],
                    self.encode_growth(b ^ pbit ^ qbit, table),
                )
                v = out.get(key, 0) - bond
                if v:
                    out[key] = v
                else:
                    out.pop(key, None)
        return out

    def iter_words(self, spanning_only: bool):
        if self.expand is None:
            self.expand = tuple(
                sum(1 << self.inner_bits[k] for k in bits_of(sub))
                for sub in range(1 << self.volume)
            )
        expand = self.expand
        total = 1 << self.volume
        slices = self.slice_masks
        if not spanning_only:
            for a_sub in range(total):
                a = expand[a_sub]
                for b_sub in range(total):
                    yield a, expand[b_sub]
            return
        for a_sub in range(total):
            a = expand[a_sub]
            for b_sub in range(total):
                b = expand[b_sub]
                support = a | b
                for lo, hi in slices:
                    if not (support & lo) or not (support & hi):
                        break
                else:
                    yield a, b


# --------------------------------------------------------------------------- #
# exact solver for systems whose every row has at most two nonzero entries     #
# --------------------------------------------------------------------------- #
@dataclass
class TwoEntrySolution:
    max_entries_per_row: int
    forced_zero: set[int]
    free_variables: list[int]
    components: list[dict[int, Fraction]]      # consistent, ungrounded components
    grounded_components: int
    inconsistent_components: int

    @property
    def kernel_dimension(self) -> int:
        return len(self.free_variables) + len(self.components)

    @property
    def active_kernel_dimension(self) -> int:
        return len(self.components)

    def basis(self) -> list[dict[int, Fraction]]:
        out = [{i: Fraction(1)} for i in self.free_variables]
        out.extend(self.components)
        return out


def solve_two_entry(rows: list[dict[int, int]], n_variables: int) -> TwoEntrySolution:
    """Exact kernel of a linear system whose rows have <= 2 nonzero entries.

    Rows with one entry ground that variable; rows with two entries force a
    fixed ratio.  Every connected component of the ratio graph therefore
    contributes exactly one kernel dimension when it is consistent and
    ungrounded, and zero otherwise.  Untouched variables are free.
    """
    max_entries = 0
    adjacency: dict[int, list[tuple[int, Fraction]]] = {}
    grounded: set[int] = set()
    touched: set[int] = set()
    for row in rows:
        items = [(k, v) for k, v in row.items() if v]
        max_entries = max(max_entries, len(items))
        if not items:
            continue
        for k, _ in items:
            touched.add(k)
        if len(items) == 1:
            grounded.add(items[0][0])
        elif len(items) == 2:
            (i, ci), (j, cj) = items
            adjacency.setdefault(i, []).append((j, Fraction(-ci, cj)))
            adjacency.setdefault(j, []).append((i, Fraction(-cj, ci)))
        else:
            raise ValueError("row with more than two entries; solver not applicable")

    forced_zero: set[int] = set()
    components: list[dict[int, Fraction]] = []
    grounded_components = 0
    inconsistent_components = 0
    visited: set[int] = set()
    for start in sorted(touched):
        if start in visited:
            continue
        values: dict[int, Fraction] = {start: Fraction(1)}
        order = [start]
        visited.add(start)
        consistent = True
        head = 0
        while head < len(order):
            node = order[head]
            head += 1
            for other, ratio in adjacency.get(node, ()):
                target = values[node] * ratio
                if other in values:
                    if values[other] != target:
                        consistent = False
                else:
                    values[other] = target
                    visited.add(other)
                    order.append(other)
        has_ground = any(node in grounded for node in values)
        if consistent and not has_ground:
            components.append(values)
        else:
            if not consistent:
                inconsistent_components += 1
            if has_ground:
                grounded_components += 1
            forced_zero.update(values)
    free_variables = [i for i in range(n_variables) if i not in touched]
    return TwoEntrySolution(
        max_entries, forced_zero, free_variables, components,
        grounded_components, inconsistent_components,
    )


# --------------------------------------------------------------------------- #
# exact sparse elimination over Q                                             #
# --------------------------------------------------------------------------- #
def rank_and_kernel(columns: list[dict[object, Fraction]], want_kernel: bool, clock: WallClock):
    pivots: dict[object, tuple[dict[object, Fraction], dict[int, Fraction]]] = {}
    row_order: dict[object, int] = {}
    kernel: list[dict[int, Fraction]] = []
    for j, col in enumerate(columns):
        clock.check(j)
        vec: dict[object, Fraction] = {r: Fraction(v) for r, v in col.items() if v}
        combo: dict[int, Fraction] = {j: Fraction(1)} if want_kernel else {}
        while vec:
            for r in vec:
                if r not in row_order:
                    row_order[r] = len(row_order)
            row = min(vec, key=lambda r: row_order[r])
            if row not in pivots:
                pivots[row] = (vec, combo)
                vec = {}
                combo = {}
                break
            pvec, pcombo = pivots[row]
            factor = vec[row] / pvec[row]
            for r, val in pvec.items():
                nv = vec.get(r, Fraction(0)) - factor * val
                if nv:
                    vec[r] = nv
                else:
                    vec.pop(r, None)
            if want_kernel:
                for c, val in pcombo.items():
                    nv = combo.get(c, Fraction(0)) - factor * val
                    if nv:
                        combo[c] = nv
                    else:
                        combo.pop(c, None)
        else:
            if want_kernel and combo:
                kernel.append(combo)
    return len(pivots), kernel


# --------------------------------------------------------------------------- #
# word classification used by Theorem M                                       #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Classification:
    kind: str            # "i", "ii-lo", "ii-hi", "iii"; lo/hi = face met by a
    lo_size: int
    hi_size: int
    thick_z_face: bool   # type ii and the a-free extreme face has >= 2 support sites
    clause_b: bool       # type ii, a-free face is a single r in b\a, inward neighbour in a


def classify(box: Box, a: int, b: int, axis: int = 0) -> Classification:
    lo, hi = 1, box.shape[axis]
    support = a | b
    lo_sites = [i for i in bits_of(support) if box.coords[i][axis] == lo]
    hi_sites = [i for i in bits_of(support) if box.coords[i][axis] == hi]
    a_lo = any(a & (1 << i) for i in lo_sites)
    a_hi = any(a & (1 << i) for i in hi_sites)
    if a_lo and a_hi:
        kind = "i"
    elif a_lo:
        kind = "ii-lo"
    elif a_hi:
        kind = "ii-hi"
    else:
        kind = "iii"
    thick = False
    clause_b = False
    if kind in ("ii-lo", "ii-hi"):
        free = hi_sites if kind == "ii-lo" else lo_sites
        thick = len(free) >= 2
        if len(free) == 1:
            r = free[0]
            if (b >> r) & 1 and not (a >> r) & 1:
                site = list(box.coords[r])
                site[axis] += -1 if kind == "ii-lo" else 1
                clause_b = bool((a >> box.index[tuple(site)]) & 1)
    return Classification(kind, len(lo_sites), len(hi_sites), thick, clause_b)


# --------------------------------------------------------------------------- #
# LOCAL row verification with plain set arithmetic (no bitmasks, no            #
# enumeration): this is the machine form of the proof of Theorem M and it runs #
# at shapes far beyond the 4^volume enumeration frontier.                     #
# --------------------------------------------------------------------------- #
def s_shape(a: frozenset, b: frozenset) -> tuple[int, ...]:
    support = a | b
    d = len(next(iter(support)))
    return tuple(
        max(s[k] for s in support) - min(s[k] for s in support) + 1 for k in range(d)
    )


def s_normalize(a: frozenset, b: frozenset) -> tuple[frozenset, frozenset]:
    support = a | b
    d = len(next(iter(support)))
    mins = tuple(min(s[k] for s in support) for k in range(d))
    shift = lambda s: tuple(s[k] - mins[k] for k in range(d))
    return frozenset(map(shift, a)), frozenset(map(shift, b))


def s_step(site: tuple[int, ...], axis: int, sign: int) -> tuple[int, ...]:
    out = list(site)
    out[axis] += sign
    return tuple(out)


def s_growth(a: frozenset, b: frozenset, axis: int, shape: tuple[int, ...]):
    """Direction-`axis` box-growing bond images of a normalized shape-`shape` word."""
    out = []
    for p in a:
        for edge, sign in ((shape[axis] - 1, 1), (0, -1)):
            if p[axis] != edge:
                continue
            q = s_step(p, axis, sign)
            out.append((s_normalize(a, b ^ {p} ^ {q}), sign))
    return out


def s_preimages(a: frozenset, b: frozenset, axis: int, shape: tuple[int, ...]):
    """All shape-`shape` preimages of a normalized shape-(shape + e_axis) word.

    Both families are reconstructed by the same explicit rule: the extreme
    axis-face of the target must carry exactly ONE support site, it must be a
    pure Z there, and its inward neighbour must be an X site.  Removing it
    inverts the growth step.
    """
    tau = s_shape(a, b)
    support = a | b
    out = []
    for edge, inward in ((tau[axis] - 1, -1), (0, 1)):
        face = [s for s in support if s[axis] == edge]
        if len(face) != 1:
            continue
        s = face[0]
        if s not in b or s in a:
            continue
        nb = s_step(s, axis, inward)
        if nb not in a:
            continue
        cand = s_normalize(a, b ^ {s} ^ {nb})
        if s_shape(*cand) == shape:
            out.append((cand, -inward))
    return out


def s_classify(a: frozenset, b: frozenset, axis: int, shape: tuple[int, ...]) -> dict:
    support = a | b
    lo = [s for s in support if s[axis] == 0]
    hi = [s for s in support if s[axis] == shape[axis] - 1]
    a_lo = any(s in a for s in lo)
    a_hi = any(s in a for s in hi)
    kind = "i" if (a_lo and a_hi) else "ii-lo" if a_lo else "ii-hi" if a_hi else "iii"
    thick = False
    clause_b = False
    if kind in ("ii-lo", "ii-hi"):
        free = hi if kind == "ii-lo" else lo
        thick = len(free) >= 2
        if len(free) == 1:
            r = free[0]
            if r in b and r not in a:
                clause_b = s_step(r, axis, 1 if kind == "ii-hi" else -1) in a
    return {"kind": kind, "thick": thick, "clause_b": clause_b,
            "lo_size": len(lo), "hi_size": len(hi)}


CLASS_LABELS = ("i", "ii-thick", "ii-clause-b", "ii-thin-no-clause-b", "iii")


def s_label(cl: dict) -> str:
    if cl["kind"] == "i":
        return "i"
    if cl["kind"] == "iii":
        return "iii"
    if cl["thick"]:
        return "ii-thick"
    return "ii-clause-b" if cl["clause_b"] else "ii-thin-no-clause-b"


def s_build(shape, axis, want, rng, tries=400):
    """Construct a spanning shape-`shape` word whose direction-`axis` class is
    `want`.  The construction is heuristic; a candidate is ACCEPTED only after
    `s_shape` and `s_classify` confirm both the shape and the class, so a
    non-None result is always a genuine witness of that class."""
    d = len(shape)
    n = shape[axis]
    sites = [tuple(c) for c in itertools.product(*(range(s) for s in shape))]
    face_lo = [s for s in sites if s[axis] == 0]
    face_hi = [s for s in sites if s[axis] == n - 1]
    middle = [s for s in sites if 0 < s[axis] < n - 1]
    if want == "iii" and n < 3:
        return None
    if want == "ii-thick" and len(face_lo) < 2:
        return None
    for _ in range(tries):
        a: set = set()
        b: set = set()
        for s in middle:
            if rng.random() < 0.25:
                a.add(s)
            if rng.random() < 0.35:
                b.add(s)
        if want == "i":
            a.add(rng.choice(face_lo))
            a.add(rng.choice(face_hi))
            for s in face_lo + face_hi:
                if rng.random() < 0.3:
                    b.add(s)
        elif want == "iii":
            a -= set(face_lo) | set(face_hi)
            if not a:
                a.add(rng.choice(middle))
            b.add(rng.choice(face_lo))
            b.add(rng.choice(face_hi))
        else:
            # a meets the high face only; the low face carries the structure
            a.add(rng.choice(face_hi))
            a -= set(face_lo)
            for s in face_hi:
                if rng.random() < 0.3:
                    b.add(s)
            if want == "ii-thick":
                b.update(rng.sample(face_lo, 2))
            else:
                r = rng.choice(face_lo)
                b -= set(face_lo)
                b.add(r)
                inward = s_step(r, axis, 1)
                if want == "ii-clause-b":
                    a.add(inward)
                else:
                    a.discard(inward)
        # patch transverse spanning; for the thin-face classes the low face must
        # stay a single site, so it is excluded from the patch pool
        if want in ("ii-clause-b", "ii-thin-no-clause-b"):
            pool = middle + face_hi
        else:
            pool = sites
        for k in range(d):
            if k == axis:
                continue
            for value in (0, shape[k] - 1):
                if any(s[k] == value for s in a | b):
                    continue
                choices = [s for s in pool if s[k] == value]
                if not choices:
                    return None
                b.add(rng.choice(choices))
        fa, fb = frozenset(a), frozenset(b)
        if not (fa | fb) or s_shape(fa, fb) != shape:
            continue
        if s_label(s_classify(fa, fb, axis, shape)) == want:
            return fa, fb
    return None


def verify_word_rows(a: frozenset, b: frozenset, axis: int, shape) -> dict:
    """Assemble every orbit row of pi Gamma^(axis) that contains (a, b)."""
    key = (tuple(sorted(a)), tuple(sorted(b)))
    sizes: list[int] = []
    coefficients: set[int] = set()
    missing = False
    for (ta, tb), _sign in s_growth(a, b, axis, shape):
        row: dict[tuple, int] = {}
        for cand, _c in s_preimages(ta, tb, axis, shape):
            k = (tuple(sorted(cand[0])), tuple(sorted(cand[1])))
            row[k] = row.get(k, 0) - 1
        nonzero = {k: v for k, v in row.items() if v}
        sizes.append(len(nonzero))
        coefficients.update(nonzero.values())
        if key not in nonzero:
            missing = True
    return {
        "row_sizes": sizes,
        "coefficients": coefficients,
        "grounded": any(s == 1 for s in sizes),
        "missing_self": missing,
        "has_rows": bool(sizes),
    }


def local_row_case(shape: tuple[int, ...], samples: int, seed: int) -> dict:
    """Class-targeted pointwise verification of Theorem M at a single shape.

    For every direction and every one of the five classes, words of exactly that
    class are constructed, the full orbit rows containing them are assembled from
    the explicit preimage rule, and the predictions of the proof are checked:
    rows have at most two entries; type-i and non-clause-(b) type-ii words own a
    single-entry (grounding) row; clause-(b) words never do.
    """
    import random

    rng = random.Random(seed)
    d = len(shape)
    clock = WallClock(f"local_{'x'.join(map(str, shape))}", SYMBOL_WALL_SECONDS)
    volume = 1
    for s in shape:
        volume *= s
    max_row = 0
    coefficient_values: set[int] = set()
    counts = {label: 0 for label in CLASS_LABELS}
    unreachable = {label: 0 for label in CLASS_LABELS}
    violations: list[dict] = []
    for axis in range(d):
        for label in CLASS_LABELS:
            for _ in range(samples):
                clock.check(sum(counts.values()))
                built = s_build(shape, axis, label, rng)
                if built is None:
                    unreachable[label] += 1
                    continue
                a, b = built
                counts[label] += 1
                info = verify_word_rows(a, b, axis, shape)
                if info["row_sizes"]:
                    max_row = max(max_row, max(info["row_sizes"]))
                coefficient_values.update(info["coefficients"])
                record = {"axis": axis, "label": label,
                          "a": sorted(a), "b": sorted(b)}
                if info["missing_self"]:
                    violations.append({**record, "reason": "word missing from its own row"})
                if label == "iii":
                    if info["has_rows"]:
                        violations.append({**record, "reason": "type iii has a growth image"})
                    continue
                if not info["has_rows"]:
                    violations.append({**record, "reason": "active word without a growth image"})
                    continue
                if label in ("i", "ii-thick", "ii-thin-no-clause-b"):
                    if not info["grounded"]:
                        violations.append({**record, "reason": "expected grounding row absent"})
                elif label == "ii-clause-b":
                    if info["grounded"]:
                        violations.append({**record, "reason": "clause-(b) word grounded"})
    return {
        "claim_tag": "[COMPUTATION]",
        "shape": list(shape),
        "dimension": d,
        "volume": volume,
        "samples_per_class_per_direction": samples,
        "witnesses_per_class": counts,
        "unconstructible_per_class": unreachable,
        "max_row_size": max_row,
        "row_coefficient_values": sorted(coefficient_values),
        "violations": violations[:8],
        "violation_count": len(violations),
        "seconds": round(clock.elapsed(), 3),
    }

# --------------------------------------------------------------------------- #
# reference densities                                                         #
# --------------------------------------------------------------------------- #
def hamiltonian_density(dimension: int, field: int = 1, bond: int = 1) -> dict[tuple, int]:
    """h = field * X_0 + bond * sum_i Z_0 Z_{e_i}: the density whose translation
    sum is H itself, so pi ad_H h = 0 for every (field, bond)."""
    out: dict[tuple, int] = {}
    one = (1,) * dimension
    out[(one, 1, 0)] = field
    for axis in range(dimension):
        shape = list(one)
        shape[axis] = 2
        shape = tuple(shape)
        strides = []
        acc = 1
        for k in range(dimension - 1, -1, -1):
            strides.insert(0, acc)
            acc *= shape[k]
        out[(shape, 0, 1 | (1 << strides[axis]))] = bond
    return out


def chain_current_density(n: int) -> dict[tuple, int]:
    """The single-shape chain density of extent n found by the exact kernel:

        j_n = - X_0 X_1 ... X_{n-2} Z_0 Z_{n-1}  +  X_1 ... X_{n-1} Z_0 Z_{n-1}
            (Hermitian form:  Y_0 X_1..X_{n-2} Z_{n-1} - Z_0 X_1..X_{n-2} Y_{n-1} )
    """
    if n < 2:
        raise ValueError("n >= 2")
    b = 1 | (1 << (n - 1))
    a_lo = (1 << (n - 1)) - 1                 # sites 0..n-2
    a_hi = ((1 << n) - 1) ^ 1                 # sites 1..n-1
    return {((n,), a_lo, b): -1, ((n,), a_hi, b): 1}


def orbit_image(key: tuple, field: int, bond: int, cache: dict) -> dict[tuple, int]:
    shape, a_code, b_code = key
    if not shape:
        return {}
    box = cache.get(shape)
    if box is None:
        box = cache[shape] = Box(shape)
    a, b = box.decode(a_code, b_code)
    out: dict[tuple, int] = {}
    for (na, nb), coeff in box.ad_terms(a, b, field, bond).items():
        k = box.canon(na, nb)
        v = out.get(k, 0) + coeff
        if v:
            out[k] = v
        else:
            out.pop(k, None)
    return out


def apply_lbar(density: dict[tuple, int], field: int, bond: int, cache: dict) -> dict[tuple, int]:
    out: dict[tuple, int] = {}
    for key, coeff in density.items():
        for k, c in orbit_image(key, field, bond, cache).items():
            v = out.get(k, 0) + coeff * c
            if v:
                out[k] = v
            else:
                out.pop(k, None)
    return out


# --------------------------------------------------------------------------- #
# Block A: global box quotients and charge bases                              #
# --------------------------------------------------------------------------- #
def global_case(dimension: int, radius: int, field: int, bond: int, want_kernel: bool) -> dict:
    clock = WallClock(f"global_d{dimension}_R{radius}", GLOBAL_WALL_SECONDS)
    box = Box((radius + 1,) * dimension)
    orbit_keys: list[tuple] = []
    seen: set[tuple] = set()
    for a, b in box.iter_words(spanning_only=False):
        clock.check(len(seen))
        key = box.canon(a, b)
        if key not in seen:
            seen.add(key)
            orbit_keys.append(key)
    orbit_keys.sort(key=lambda k: (len(k[0]), k[0], k[1], k[2]))
    cache: dict[tuple, Box] = {}
    columns = []
    for key in orbit_keys:
        clock.check(len(columns))
        columns.append({k: Fraction(v) for k, v in orbit_image(key, field, bond, cache).items()})
    rank, kernel = rank_and_kernel(columns, want_kernel, clock)
    n = len(orbit_keys)
    payload: dict[str, object] = {
        "claim_tag": "[COMPUTATION]",
        "dimension": dimension,
        "radius": radius,
        "field": field,
        "bond": bond,
        "word_count": 1 << (2 * box.volume),
        "orbit_count": n,
        "rank_Lbar": rank,
        "kernel_dimension": n - rank,
        "quotient_dimension": n - rank - 1,
        "seconds": round(clock.elapsed(), 3),
        "peak_rss_mib": round(peak_rss_bytes() / 1024**2, 3),
    }
    if want_kernel:
        payload["kernel_report"] = kernel_report(orbit_keys, kernel)
    return payload


def kernel_report(orbit_keys, kernel) -> dict:
    """Theorem M predicates on the exact kernel of Lbar|V_R."""
    cache: dict[tuple, Box] = {}
    rows = []
    for vec in kernel:
        support = {orbit_keys[c]: v for c, v in vec.items() if v}
        shapes = sorted({k[0] for k in support if k[0]})
        if not shapes:
            rows.append({"identity_only": True})
            continue
        n1 = max(s[0] for s in shapes)
        entries = []
        for key in sorted(k for k in support if k[0] and k[0][0] == n1):
            shape = key[0]
            box = cache.get(shape) or cache.setdefault(shape, Box(shape))
            cl = classify(box, *box.decode(key[1], key[2]))
            entries.append(
                {
                    "shape": list(shape), "a": key[1], "b": key[2],
                    "coefficient": str(support[key]), "kind": cl.kind,
                    "lo_size": cl.lo_size, "hi_size": cl.hi_size,
                    "thick_z_face": cl.thick_z_face, "clause_b": cl.clause_b,
                }
            )
        rows.append(
            {
                "identity_only": False,
                "shapes": [list(s) for s in shapes],
                "max_sigma1": n1,
                "top_layer": entries,
                "violates_T1": any(e["kind"] == "i" for e in entries),
                "violates_T2": any(e["thick_z_face"] for e in entries),
                "top_all_clause_a_or_b": all(
                    e["kind"] == "iii" or e["clause_b"] for e in entries
                ),
            }
        )
    return {"vectors": rows}


# --------------------------------------------------------------------------- #
# Block B: leading symbol on a single shape                                   #
# --------------------------------------------------------------------------- #
def symbol_rows(box: Box, words, axis: int, bond: int, clock: WallClock):
    """Rows of pi Gamma^(axis) indexed by target orbit key.

    `Box.gamma_keys` already returns canonical (translation-normalised) keys,
    which is exactly the projection pi; see the growth tables in `Box`.
    """
    rows: dict[tuple, dict[int, int]] = {}
    gamma = box.gamma_keys
    for i, (a, b) in enumerate(words):
        if not (i & 0xFFFF):
            clock.check(i)
        for key, coeff in gamma(a, b, axis, bond).items():
            row = rows.setdefault(key, {})
            v = row.get(i, 0) + coeff
            if v:
                row[i] = v
            else:
                row.pop(i, None)
    return rows


def symbol_case(shape: tuple[int, ...], bond: int, crosscheck: bool) -> dict:
    clock = WallClock(f"symbol_{'x'.join(map(str, shape))}", SYMBOL_WALL_SECONDS)
    box = Box(shape)
    words = list(box.iter_words(spanning_only=True))
    d = len(shape)

    rows1 = symbol_rows(box, words, 0, bond, clock)
    target = (shape[0] + 1,) + tuple(shape[1:])
    target_ok = all(key[0] == target for key in rows1)
    sol1 = solve_two_entry(list(rows1.values()), len(words))

    classes = [classify(box, a, b) for a, b in words]
    forced = sol1.forced_zero
    t1 = [i for i, c in enumerate(classes) if c.kind == "i"]
    t2 = [i for i, c in enumerate(classes) if c.thick_z_face]
    surviving_active = [
        i for i in range(len(words))
        if i not in forced and i not in set(sol1.free_variables)
    ]
    kinds: dict[str, int] = {}
    violations = []
    for i in surviving_active:
        c = classes[i]
        kinds[c.kind] = kinds.get(c.kind, 0) + 1
        if c.kind != "iii" and not c.clause_b:
            violations.append({"a": words[i][0], "b": words[i][1], "kind": c.kind})

    # combined all-direction leading symbol (valid when sigma is maximal in every
    # direction, e.g. for single-shape densities)
    combined_rows: list[dict[int, int]] = []
    for axis in range(d):
        combined_rows.extend(symbol_rows(box, words, axis, bond, clock).values())
    sol_all = solve_two_entry(combined_rows, len(words))

    witness = None
    if sol1.components:
        comp = sol1.components[0]
        witness = [
            {
                "a_sites": [list(box.coords[j]) for j in bits_of(words[i][0])],
                "b_sites": [list(box.coords[j]) for j in bits_of(words[i][1])],
                "coefficient": str(v),
                "kind": classes[i].kind,
                "lo_size": classes[i].lo_size,
                "hi_size": classes[i].hi_size,
                "clause_b": classes[i].clause_b,
            }
            for i, v in sorted(comp.items())
        ]

    payload = {
        "claim_tag": "[COMPUTATION]",
        "shape": list(shape),
        "dimension": d,
        "bond": bond,
        "word_count": len(words),
        "target_shape": list(target),
        "all_images_in_target_shape": bool(target_ok),
        "row_count": len(rows1),
        "max_nonzero_entries_per_row": sol1.max_entries_per_row,
        "annihilated_word_count": len(sol1.free_variables),
        "active_word_count": len(words) - len(sol1.free_variables),
        "kernel_dimension": sol1.kernel_dimension,
        "active_kernel_dimension": sol1.active_kernel_dimension,
        "grounded_components": sol1.grounded_components,
        "inconsistent_components": sol1.inconsistent_components,
        "active_kernel_witness": witness,
        "T1_words": len(t1),
        "T1_all_forced_zero": all(i in forced for i in t1),
        "T2_words": len(t2),
        "T2_all_forced_zero": all(i in forced for i in t2),
        "surviving_active_count": len(surviving_active),
        "surviving_active_kinds": kinds,
        "clause_b_violations": violations,
        "combined_kernel_dimension": sol_all.kernel_dimension,
        "combined_active_kernel_dimension": sol_all.active_kernel_dimension,
        "combined_free_variable_count": len(sol_all.free_variables),
        "seconds": round(clock.elapsed(), 3),
        "peak_rss_mib": round(peak_rss_bytes() / 1024**2, 3),
    }
    if crosscheck and len(words) <= 4000:
        cols = []
        row_ids = {k: n for n, k in enumerate(rows1)}
        for i in range(len(words)):
            cols.append({})
        for key, row in rows1.items():
            for i, v in row.items():
                cols[i][row_ids[key]] = Fraction(v)
        rank, kern = rank_and_kernel(cols, True, clock)
        payload["gauss_rank"] = rank
        payload["gauss_kernel_dimension"] = len(words) - rank
        payload["graph_matches_gauss"] = (len(words) - rank) == sol1.kernel_dimension
        gauss_forced = {
            i for i in range(len(words))
            if all(v.get(i, Fraction(0)) == 0 for v in kern)
        }
        payload["graph_forced_matches_gauss"] = gauss_forced == forced
    if crosscheck:
        # Lemma B, verified term by term: the shape-(sigma + e_1) component of
        # pi ad_H restricted to V_sigma IS pi Gamma^(1) -- nothing else in the
        # commutator can reach the enlarged direction-1 extent.
        identification_ok = True
        for i, (a, b) in enumerate(words):
            direct: dict[tuple, int] = {}
            for (na, nb), coeff in box.ad_terms(a, b, 1, bond).items():
                k = box.canon(na, nb)
                if k[0] != target:
                    continue
                v = direct.get(k, 0) + coeff
                if v:
                    direct[k] = v
                else:
                    direct.pop(k, None)
            if direct != box.gamma_keys(a, b, 0, bond):
                identification_ok = False
                break
        payload["leading_symbol_identification_ok"] = identification_ok
    return payload


# --------------------------------------------------------------------------- #
# Block C: single-shape charge kernels                                        #
# --------------------------------------------------------------------------- #
def full_shape_case(shape: tuple[int, ...], field: int, bond: int) -> dict:
    clock = WallClock(f"full_{'x'.join(map(str, shape))}", FULL_WALL_SECONDS)
    box = Box(shape)
    words = list(box.iter_words(spanning_only=True))
    d = len(shape)

    combined: list[dict[int, int]] = []
    for axis in range(d):
        combined.extend(symbol_rows(box, words, axis, bond, clock).values())
    sol = solve_two_entry(combined, len(words))
    basis = sol.basis()

    # Lbar restricted to the (proved) combined leading-symbol kernel
    image_cache: dict[int, dict[tuple, int]] = {}
    columns: list[dict[object, Fraction]] = []
    for n, vec in enumerate(basis):
        clock.check(n)
        col: dict[object, Fraction] = {}
        for i, weight in vec.items():
            img = image_cache.get(i)
            if img is None:
                a, b = words[i]
                img = {}
                for (na, nb), coeff in box.ad_terms(a, b, field, bond).items():
                    k = box.canon(na, nb)
                    v = img.get(k, 0) + coeff
                    if v:
                        img[k] = v
                    else:
                        img.pop(k, None)
                image_cache[i] = img
            for k, coeff in img.items():
                nv = col.get(k, Fraction(0)) + weight * coeff
                if nv:
                    col[k] = nv
                else:
                    col.pop(k, None)
        columns.append(col)
    rank, kernel = rank_and_kernel(columns, True, clock)

    witness = None
    if kernel:
        vec = kernel[0]
        combined_word: dict[int, Fraction] = {}
        for j, w in vec.items():
            for i, v in basis[j].items():
                nv = combined_word.get(i, Fraction(0)) + w * v
                if nv:
                    combined_word[i] = nv
                else:
                    combined_word.pop(i, None)
        witness = [
            {
                "a_sites": [list(box.coords[j]) for j in bits_of(words[i][0])],
                "b_sites": [list(box.coords[j]) for j in bits_of(words[i][1])],
                "coefficient": str(v),
                "kind": classify(box, *words[i]).kind,
            }
            for i, v in sorted(combined_word.items())
        ]

    return {
        "claim_tag": "[COMPUTATION]",
        "shape": list(shape),
        "dimension": d,
        "field": field,
        "bond": bond,
        "word_count": len(words),
        "combined_symbol_kernel_dimension": len(basis),
        "rank_on_symbol_kernel": rank,
        "kernel_dimension": len(basis) - rank,
        "kernel_witness": witness,
        "seconds": round(clock.elapsed(), 3),
        "peak_rss_mib": round(peak_rss_bytes() / 1024**2, 3),
    }


# --------------------------------------------------------------------------- #
def check(name: str, passed: bool, detail: str) -> dict[str, object]:
    return {"name": name, "passed": bool(passed), "detail": detail}


def main() -> int:
    started = time.monotonic()
    checks: list[dict[str, object]] = []
    walls: list[dict[str, object]] = []
    cache: dict[tuple, Box] = {}

    # -- reference densities ------------------------------------------------ #
    reference: dict[str, object] = {}
    for dimension in (1, 2, 3):
        ok = [
            not apply_lbar(hamiltonian_density(dimension, f, g), f, g, cache)
            for f, g in COEFFICIENTS
        ]
        reference[f"hamiltonian_d{dimension}_conserved"] = ok
        checks.append(
            check(
                f"hamiltonian_density_conserved_d{dimension}",
                all(ok),
                "pi ad_H h = 0 exactly for every tested (field, bond)",
            )
        )
    chain_ok = {}
    for n in range(2, 9):
        j = chain_current_density(n)
        chain_ok[n] = [not apply_lbar(j, f, g, cache) for f, g in COEFFICIENTS]
    reference["chain_single_shape_charges_conserved"] = {
        str(k): v for k, v in chain_ok.items()
    }
    checks.append(
        check(
            "chain_single_shape_charge_family_conserved",
            all(all(v) for v in chain_ok.values()),
            "d=1: the single-shape family j_n (n = 2..8) is exactly conserved; "
            "these are the words a 1D-blind lemma would have to kill",
        )
    )
    j_kinds = []
    for n in (2, 3, 4, 5):
        box = Box((n,))
        for key in chain_current_density(n):
            cl = classify(box, *box.decode(key[1], key[2]))
            j_kinds.append(
                {"n": n, "kind": cl.kind, "clause_b": cl.clause_b,
                 "thick": cl.thick_z_face, "lo": cl.lo_size, "hi": cl.hi_size}
            )
    reference["chain_charge_word_types"] = j_kinds
    checks.append(
        check(
            "chain_charges_are_type_ii_clause_b_and_thin",
            all(
                k["kind"] in ("ii-lo", "ii-hi") and k["clause_b"] and not k["thick"]
                and k["lo"] == 1 and k["hi"] == 1
                for k in j_kinds
            ),
            "1D control: every word of an actual 1D charge is type ii, satisfies "
            "clause (b), and has THIN faces -- so 'type ii implies zero' is FALSE "
            "and the thick-face hypothesis of T2 is indispensable",
        )
    )

    # -- Block A ------------------------------------------------------------ #
    global_rows = []
    for dimension, radius in GLOBAL_CASES:
        want_kernel = (dimension, radius) in GLOBAL_KERNEL_CASES
        try:
            row = global_case(dimension, radius, 1, 1, want_kernel)
        except ResourceWall as wall:
            payload = dict(wall.payload)
            payload.update({"dimension": dimension, "radius": radius})
            walls.append(payload)
            continue
        global_rows.append(row)
        print(
            f"[global] d={dimension} R={radius} orbits={row['orbit_count']} "
            f"rank={row['rank_Lbar']} quotient={row['quotient_dimension']} "
            f"({row['seconds']:.1f}s)",
            flush=True,
        )

    expected = {
        (1, 0): 0, (1, 1): 2, (1, 2): 4, (1, 3): 6, (1, 4): 8, (1, 5): 10,
        (2, 0): 0, (2, 1): 1, (3, 0): 0,
    }
    mismatches = [
        (r["dimension"], r["radius"], r["quotient_dimension"])
        for r in global_rows
        if r["quotient_dimension"] != expected.get((r["dimension"], r["radius"]))
    ]
    checks.append(
        check(
            "global_quotients_are_2R_in_d1_and_1_for_d_ge_2",
            not mismatches,
            f"mismatches={mismatches}",
        )
    )

    cross = []
    if E117_RESULT.exists():
        stored = json.loads(E117_RESULT.read_text())
        for case in stored.get("data", {}).get("cases", []):
            dim, rad = case.get("dimension"), case.get("radius")
            quot = case.get("nontrivial_quotient_dimension")
            if quot is None:
                continue
            mine = next(
                (r for r in global_rows
                 if r["dimension"] == dim and r["radius"] == rad),
                None,
            )
            cross.append(
                {
                    "dimension": dim,
                    "radius": rad,
                    "arithmetic": case.get("arithmetic"),
                    "e117_quotient": quot,
                    "e117_translation_rank": case.get("translation_rank"),
                    "e117_commutator_rank": case.get("commutator_rank"),
                    "e127_quotient": None if mine is None else mine["quotient_dimension"],
                    "e127_orbit_count": None if mine is None else mine["orbit_count"],
                    "e127_rank": None if mine is None else mine["rank_Lbar"],
                    "agree": mine is None
                    or (
                        mine["quotient_dimension"] == quot
                        and mine["orbit_count"] == case.get("translation_rank")
                        and mine["rank_Lbar"] == case.get("commutator_rank")
                    ),
                }
            )
    checks.append(
        check(
            "cross_reference_e117_quotients",
            bool(cross) and all(r["agree"] for r in cross),
            f"{len(cross)} stored rows compared, "
            f"{sum(1 for r in cross if r['e127_quotient'] is not None)} recomputed "
            "here (orbit count, rank and quotient all matched)",
        )
    )

    # -- Block B ------------------------------------------------------------ #
    symbol_results = []
    for shape in SHAPES_SYMBOL:
        try:
            row = symbol_case(shape, 1, shape in SHAPES_CROSSCHECK)
        except ResourceWall as wall:
            payload = dict(wall.payload)
            payload["shape"] = list(shape)
            walls.append(payload)
            continue
        symbol_results.append(row)
        print(
            f"[symbol] {shape} words={row['word_count']} maxrow={row['max_nonzero_entries_per_row']} "
            f"act.ker={row['active_kernel_dimension']} comb.act.ker={row['combined_active_kernel_dimension']} "
            f"T1={row['T1_words']}/{row['T1_all_forced_zero']} T2={row['T2_words']}/{row['T2_all_forced_zero']} "
            f"surv={row['surviving_active_count']} ({row['seconds']:.1f}s)",
            flush=True,
        )

    checks.append(
        check(
            "leading_symbol_rows_have_at_most_two_entries",
            bool(symbol_results) and all(r["max_nonzero_entries_per_row"] <= 2 for r in symbol_results),
            "structure theorem: every orbit row of pi Gamma^(1) has <= 2 entries",
        )
    )
    checks.append(
        check(
            "leading_symbol_identification_lemma_B",
            any("leading_symbol_identification_ok" in r for r in symbol_results)
            and all(
                r["leading_symbol_identification_ok"]
                for r in symbol_results
                if "leading_symbol_identification_ok" in r
            ),
            "term-by-term: the shape-(sigma+e_1) component of pi ad_H on V_sigma "
            "equals pi Gamma^(1) -- no field or shrinking bond term reaches it",
        )
    )
    checks.append(
        check(
            "graph_solver_matches_gaussian_elimination",
            all(r.get("graph_matches_gauss", True) and r.get("graph_forced_matches_gauss", True)
                for r in symbol_results),
            "the O(rows) grounding solver reproduces exact Gaussian elimination "
            "(kernel dimension and forced-zero set) on the cross-checked shapes",
        )
    )
    checks.append(
        check(
            "T1_holds_on_every_shape",
            all(r["T1_all_forced_zero"] for r in symbol_results),
            "type-i words (X-support on both extreme 1-faces) are forced to zero",
        )
    )
    checks.append(
        check(
            "T2_holds_on_every_shape",
            all(r["T2_all_forced_zero"] for r in symbol_results),
            "type-ii words whose a-free extreme 1-face carries >= 2 support sites "
            "are forced to zero",
        )
    )
    checks.append(
        check(
            "surviving_active_words_satisfy_clause_b",
            all(not r["clause_b_violations"] for r in symbol_results),
            "every surviving active word satisfies clause (b) of Theorem M",
        )
    )
    d1 = [r for r in symbol_results if r["dimension"] == 1]
    dge2 = [r for r in symbol_results if r["dimension"] >= 2]
    checks.append(
        check(
            "d1_control_active_kernel_nonzero",
            bool(d1) and all(r["active_kernel_dimension"] > 0 for r in d1),
            "d=1 control: the leading symbol keeps a nonzero ACTIVE kernel at "
            "every chain shape, so no lemma may kill all active words",
        )
    )
    checks.append(
        check(
            "d1_control_T2_hypothesis_vacuous",
            bool(d1) and all(r["T2_words"] == 0 for r in d1),
            "d=1 control: no chain word has a thick face, so T2 is VACUOUS in d=1 "
            "-- this is exactly where d >= 2 is used",
        )
    )
    checks.append(
        check(
            "dge2_T2_hypothesis_has_content",
            bool(dge2) and all(r["T2_words"] > 0 for r in dge2),
            "d >= 2: thick-face words exist at every tested shape",
        )
    )
    checks.append(
        check(
            "leading_symbol_alone_does_not_close",
            bool(dge2) and any(r["combined_active_kernel_dimension"] > 0 for r in dge2),
            "obstruction: even the all-direction leading symbol retains a nonzero "
            "active kernel for d >= 2, so the subleading level is required",
        )
    )

    # -- Block B' : pointwise verification beyond the enumeration frontier --- #
    local_results = []
    for shape in SHAPES_LOCAL:
        try:
            row = local_row_case(shape, LOCAL_SAMPLES, seed=20260816)
        except ResourceWall as wall:
            payload = dict(wall.payload)
            payload["shape"] = list(shape)
            walls.append(payload)
            continue
        local_results.append(row)
        print(
            f"[local] {shape} vol={row['volume']} "
            f"maxrow={row['max_row_size']} coeffs={row['row_coefficient_values']} "
            f"classes={row['witnesses_per_class']} "
            f"violations={row['violation_count']} ({row['seconds']:.1f}s)",
            flush=True,
        )
    checks.append(
        check(
            "local_rows_verify_theorem_M_pointwise",
            bool(local_results) and all(r["violation_count"] == 0 for r in local_results)
            and all(r["max_row_size"] <= 2 for r in local_results),
            "independent set-arithmetic check at shapes beyond the 4^volume "
            "frontier: every row has <= 2 entries, every type-i and every "
            "non-clause-(b) type-ii word owns a single-entry (grounding) row, and "
            "no clause-(b) word is grounded",
        )
    )
    local_d1 = [r for r in local_results if r["dimension"] == 1]
    local_dge2 = [r for r in local_results if r["dimension"] >= 2]
    checks.append(
        check(
            "local_d1_control_no_thick_faces",
            bool(local_d1)
            and all(r["witnesses_per_class"]["ii-thick"] == 0 for r in local_d1),
            "d=1 control at large chain shapes: thick faces do not exist, so T2 "
            "never fires; clause-(b) words are present instead",
        )
    )
    checks.append(
        check(
            "local_d1_control_clause_b_words_present",
            bool(local_d1)
            and all(r["witnesses_per_class"]["ii-clause-b"] > 0 for r in local_d1),
            "d=1 control: surviving clause-(b) words exist at large chain shapes, "
            "matching the 2R tower",
        )
    )
    checks.append(
        check(
            "local_dge2_thick_faces_dominate",
            bool(local_dge2)
            and all(r["witnesses_per_class"]["ii-thick"] > 0 for r in local_dge2),
            "d >= 2 at large shapes: thick-face words are present, so T2 has content "
            "at every tested shape including three-dimensional ones",
        )
    )

    # -- Block C ------------------------------------------------------------ #
    full_results = []
    for shape in SHAPES_FULL:
        for field, bond in ((1, 1), (2, 3)):
            try:
                row = full_shape_case(shape, field, bond)
            except ResourceWall as wall:
                payload = dict(wall.payload)
                payload["shape"] = list(shape)
                walls.append(payload)
                continue
            full_results.append(row)
            print(
                f"[full] {shape} (f,g)=({field},{bond}) words={row['word_count']} "
                f"symker={row['combined_symbol_kernel_dimension']} "
                f"ker={row['kernel_dimension']} ({row['seconds']:.1f}s)",
                flush=True,
            )

    d1_full = [r for r in full_results if r["dimension"] == 1]
    dge2_full = [r for r in full_results if r["dimension"] >= 2]
    checks.append(
        check(
            "d1_single_shape_charge_kernel_is_one_dimensional",
            bool(d1_full) and all(r["kernel_dimension"] == 1 for r in d1_full),
            "d=1 control: every chain shape carries exactly ONE single-shape charge "
            "density (the free-fermion current family j_n)",
        )
    )
    checks.append(
        check(
            "dge2_no_single_shape_charge_density",
            bool(dge2_full) and all(r["kernel_dimension"] == 0 for r in dge2_full),
            "d >= 2: no nonzero single-shape charge density at any tested shape",
        )
    )

    payload = {
        "provenance": SCRIPT,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "conventions": {
            "word": "X^a Z^b with a, b finite subsets of Z^d",
            "operator": "Lbar = pi . (1/2)[H, .], H = field*sum X + bond*sum ZZ",
            "quotient": "orbit_count - rank(Lbar) - 1 (identity class removed)",
            "faces": "for direction i, the extreme i-faces of the support bounding box",
        },
        "coefficients_tested": [list(c) for c in COEFFICIENTS],
        "reference_densities": reference,
        "global_cases": global_rows,
        "e117_cross_reference": cross,
        "leading_symbol_cases": symbol_results,
        "local_row_cases": local_results,
        "single_shape_cases": full_results,
        "resource_walls": walls,
        "checks": checks,
        "all_checks_passed": all(c["passed"] for c in checks),
        "total_seconds": round(time.monotonic() - started, 3),
        "peak_rss_mib": round(peak_rss_bytes() / 1024**2, 3),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n")
    print(f"\nwrote {RESULT.relative_to(ROOT)}")
    for c in checks:
        print(f"  [{'PASS' if c['passed'] else 'FAIL'}] {c['name']}: {c['detail']}")
    print(f"all_checks_passed={payload['all_checks_passed']}")
    return 0 if payload["all_checks_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
