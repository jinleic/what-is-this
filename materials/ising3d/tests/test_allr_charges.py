#!/usr/bin/env python3
"""Standalone verifier for proofs/allr_charges.md and results/integrability/allr_charges.json.

This test is an INDEPENDENT reimplementation.  It never imports
`experiments/e127_allr_charges.py`; every decisive quantity is recomputed here
from plain set arithmetic on coordinate tuples, with its own exact linear
algebra (union-find ratio propagation AND, on the small shapes, unrelated
fraction-free Gaussian elimination).  It then compares against every stored row.

Run:  PYTHONPATH=src .venv/bin/python tests/test_allr_charges.py
"""
from __future__ import annotations

import itertools
import json
import random
import sys
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "allr_charges.json"

FAILURES: list[str] = []
CHECKS = 0


def expect(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        FAILURES.append(message)


def expect_eq(got, want, message: str) -> None:
    expect(got == want, f"{message}: got {got!r}, want {want!r}")


# --------------------------------------------------------------------------- #
# independent Pauli / orbit arithmetic (frozensets of integer coordinates)     #
# --------------------------------------------------------------------------- #
Site = tuple[int, ...]
Word = tuple[frozenset, frozenset]


def support(word: Word) -> frozenset:
    return word[0] | word[1]


def box_of(word: Word) -> tuple[tuple[int, ...], tuple[int, ...]]:
    sup = support(word)
    d = len(next(iter(sup)))
    mins = tuple(min(s[k] for s in sup) for k in range(d))
    maxs = tuple(max(s[k] for s in sup) for k in range(d))
    return mins, maxs


def shape_of(word: Word) -> tuple[int, ...]:
    mins, maxs = box_of(word)
    return tuple(maxs[k] - mins[k] + 1 for k in range(len(mins)))


def normalise(word: Word) -> Word:
    mins, _ = box_of(word)
    d = len(mins)

    def shift(s: Site) -> Site:
        return tuple(s[k] - mins[k] for k in range(d))

    return frozenset(map(shift, word[0])), frozenset(map(shift, word[1]))


def step(s: Site, axis: int, sign: int) -> Site:
    out = list(s)
    out[axis] += sign
    return tuple(out)


def neighbours(s: Site) -> list[Site]:
    out = []
    for axis in range(len(s)):
        for sign in (-1, 1):
            out.append(step(s, axis, sign))
    return out


def pi_ad(word: Word, field: int, bond: int) -> dict[Word, int]:
    """pi . (1/2)[H, X^a Z^b] with H = field*sum X + bond*sum ZZ.

    Field term:  p in b            ->  (a symdiff p, b)      coefficient +field
    Bond  term:  p in a, q ~ p, q not in a
                                   ->  (a, b symdiff p symdiff q) coefficient -bond
    Images are normalised, i.e. projected to translation orbits.
    """
    a, b = word
    out: dict[Word, int] = {}

    def add(w: Word, c: int) -> None:
        key = normalise(w) if support(w) else (frozenset(), frozenset())
        v = out.get(key, 0) + c
        if v:
            out[key] = v
        else:
            out.pop(key, None)

    for p in b:
        add((a ^ {p}, b), field)
    for p in a:
        for q in neighbours(p):
            if q in a:
                continue
            add((a, b ^ {p} ^ {q}), -bond)
    return out


def pi_ad_density(density: dict[Word, int], field: int, bond: int) -> dict[Word, int]:
    out: dict[Word, int] = {}
    for word, coeff in density.items():
        for image, c in pi_ad(word, field, bond).items():
            v = out.get(image, 0) + coeff * c
            if v:
                out[image] = v
            else:
                out.pop(image, None)
    return out


def growth_images(word: Word, axis: int, shape: tuple[int, ...], bond: int) -> dict[Word, int]:
    """Box-growing direction-`axis` bond terms of a normalised shape-`shape` word."""
    a, b = word
    out: dict[Word, int] = {}
    for p in a:
        for edge, sign in ((shape[axis] - 1, 1), (0, -1)):
            if p[axis] != edge:
                continue
            q = step(p, axis, sign)
            key = normalise((a, b ^ {p} ^ {q}))
            out[key] = out.get(key, 0) - bond
    return {k: v for k, v in out.items() if v}


def face_sites(word: Word, axis: int, shape: tuple[int, ...]):
    sup = support(word)
    lo = [s for s in sup if s[axis] == 0]
    hi = [s for s in sup if s[axis] == shape[axis] - 1]
    return lo, hi


def classify(word: Word, axis: int, shape: tuple[int, ...]) -> dict:
    a, b = word
    lo, hi = face_sites(word, axis, shape)
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
                clause_b = step(r, axis, 1 if kind == "ii-hi" else -1) in a
    return {"kind": kind, "thick": thick, "clause_b": clause_b,
            "lo_size": len(lo), "hi_size": len(hi)}


_SHAPE_CACHE: dict[tuple[int, ...], list] = {}


def enumerate_shape(shape: tuple[int, ...]) -> list[Word]:
    """All normalised words whose bounding box is exactly `shape`."""
    cached = _SHAPE_CACHE.get(shape)
    if cached is not None:
        return cached
    sites = [tuple(c) for c in itertools.product(*(range(s) for s in shape))]
    volume = len(sites)
    subsets = []
    for code in range(1 << volume):
        subsets.append(frozenset(sites[k] for k in range(volume) if (code >> k) & 1))
    lo_masks = []
    hi_masks = []
    for axis in range(len(shape)):
        lo_masks.append(frozenset(s for s in sites if s[axis] == 0))
        hi_masks.append(frozenset(s for s in sites if s[axis] == shape[axis] - 1))
    out: list[Word] = []
    for a in subsets:
        for b in subsets:
            sup = a | b
            if not sup:
                continue
            ok = True
            for axis in range(len(shape)):
                if not (sup & lo_masks[axis]) or not (sup & hi_masks[axis]):
                    ok = False
                    break
            if ok:
                out.append((a, b))
    _SHAPE_CACHE[shape] = out
    return out


# --------------------------------------------------------------------------- #
# independent exact solvers                                                   #
# --------------------------------------------------------------------------- #
class WeightedUnionFind:
    """Union-find with multiplicative weights, convention x_i = weight[i]*x_parent[i].

    After `value(i)` the node points straight at its root and `weight[i]` is the
    accumulated ratio x_i / x_root.
    """

    def __init__(self, n: int) -> None:
        self.parent = list(range(n))
        self.weight = [Fraction(1)] * n
        self.zero = [False] * n

    def value(self, i: int) -> tuple[int, Fraction]:
        path = []
        node = i
        while self.parent[node] != node:
            path.append(node)
            node = self.parent[node]
        root = node
        for child in reversed(path):
            parent = self.parent[child]
            self.weight[child] = self.weight[child] * self.weight[parent]
            if self.zero[parent]:
                self.zero[child] = True
            self.parent[child] = root
        return root, self.weight[i]

    def union(self, i: int, j: int, ratio: Fraction) -> None:
        """Impose x_j = ratio * x_i."""
        ri, wi = self.value(i)
        rj, wj = self.value(j)
        if ri == rj:
            if wj != ratio * wi:
                self.zero[ri] = True
            return
        self.parent[rj] = ri
        self.weight[rj] = ratio * wi / wj
        if self.zero[rj] or self.zero[ri]:
            self.zero[ri] = True
            self.zero[rj] = True

    def mark_zero(self, i: int) -> None:
        root, _ = self.value(i)
        self.zero[root] = True
        self.zero[i] = True


def solve_two_entry_rows(rows, n: int):
    """Kernel of a system whose rows have at most two entries, by union-find."""
    max_entries = 0
    touched = set()
    ground: list[int] = []
    pairs: list[tuple[int, int, Fraction]] = []
    for row in rows:
        items = [(k, v) for k, v in row.items() if v]
        max_entries = max(max_entries, len(items))
        for k, _ in items:
            touched.add(k)
        if len(items) == 1:
            ground.append(items[0][0])
        elif len(items) == 2:
            (i, ci), (j, cj) = items
            pairs.append((i, j, Fraction(-ci, cj)))
        elif len(items) > 2:
            raise AssertionError("row with more than two entries")
    uf = WeightedUnionFind(n)
    for i, j, ratio in pairs:
        uf.union(i, j, ratio)
    for i in ground:
        uf.mark_zero(i)
    roots: dict[int, list[int]] = {}
    for i in sorted(touched):
        root, _ = uf.value(i)
        roots.setdefault(root, []).append(i)
    forced: set[int] = set()
    components: list[dict[int, Fraction]] = []
    for root, members in roots.items():
        if uf.zero[root] or any(uf.zero[m] for m in members):
            forced.update(members)
            continue
        vec = {}
        for m in members:
            _, w = uf.value(m)
            vec[m] = w
        components.append(vec)
    free = [i for i in range(n) if i not in touched]
    return {
        "max_entries": max_entries,
        "forced_zero": forced,
        "free": free,
        "components": components,
        "kernel_dimension": len(free) + len(components),
        "active_kernel_dimension": len(components),
    }


def gauss_rank_kernel(columns, want_kernel: bool):
    """Exact rank and kernel over Q, by fraction elimination on sparse columns."""
    pivots: dict[object, tuple[dict, dict]] = {}
    order: dict[object, int] = {}
    kernel = []
    for j, col in enumerate(columns):
        vec = {r: Fraction(v) for r, v in col.items() if v}
        combo = {j: Fraction(1)} if want_kernel else {}
        while vec:
            for r in vec:
                if r not in order:
                    order[r] = len(order)
            row = min(vec, key=lambda r: order[r])
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
# reference densities, recomputed here                                        #
# --------------------------------------------------------------------------- #
def hamiltonian_density(d: int, field: int, bond: int) -> dict[Word, int]:
    origin = (0,) * d
    out: dict[Word, int] = {(frozenset({origin}), frozenset()): field}
    for axis in range(d):
        other = step(origin, axis, 1)
        out[(frozenset(), frozenset({origin, other}))] = bond
    return out


def chain_charge(n: int) -> dict[Word, int]:
    b = frozenset({(0,), (n - 1,)})
    a_lo = frozenset((k,) for k in range(n - 1))
    a_hi = frozenset((k,) for k in range(1, n))
    return {(a_lo, b): -1, (a_hi, b): 1}


# --------------------------------------------------------------------------- #
# blocks                                                                      #
# --------------------------------------------------------------------------- #
def check_reference(stored) -> None:
    couplings = [tuple(c) for c in stored["coefficients_tested"]]
    for d in (1, 2, 3):
        for field, bond in couplings:
            image = pi_ad_density(hamiltonian_density(d, field, bond), field, bond)
            expect(not image, f"h is not conserved in d={d} at (f,g)=({field},{bond})")
    for n in range(2, 9):
        for field, bond in couplings:
            image = pi_ad_density(chain_charge(n), field, bond)
            expect(not image, f"j_{n} is not conserved at (f,g)=({field},{bond})")
    # 1D DISCRIMINATOR: every word of a genuine 1D charge is type ii, satisfies
    # clause (b) and has THIN faces.  Any lemma of the form "type ii implies
    # zero", or T2 stated without the thick-face hypothesis, is therefore FALSE
    # in d = 1 and this assertion is what refutes it.
    for n in range(2, 9):
        for word in chain_charge(n):
            cl = classify(word, 0, (n,))
            expect(
                cl["kind"] in ("ii-lo", "ii-hi"),
                f"j_{n} word is not type ii (1D discriminator broken)",
            )
            expect(cl["clause_b"], f"j_{n} word fails clause (b)")
            expect(not cl["thick"], f"j_{n} word has a thick face; T2 would kill it")
            expect_eq((cl["lo_size"], cl["hi_size"]), (1, 1),
                      f"j_{n} face sizes")


def check_global(stored) -> None:
    for row in stored["global_cases"]:
        d = row["dimension"]
        R = row["radius"]
        box = [tuple(c) for c in itertools.product(range(R + 1), repeat=d)]
        volume = len(box)
        subsets = []
        for code in range(1 << volume):
            subsets.append(frozenset(box[k] for k in range(volume) if (code >> k) & 1))
        orbits: dict[Word, int] = {}
        keys: list[Word] = []
        for a in subsets:
            for b in subsets:
                key = normalise((a, b)) if (a | b) else (frozenset(), frozenset())
                if key not in orbits:
                    orbits[key] = len(keys)
                    keys.append(key)
        columns = []
        for key in keys:
            if not support(key):
                columns.append({})
                continue
            columns.append(
                {img: Fraction(c) for img, c in pi_ad(key, row["field"], row["bond"]).items()}
            )
        rank, kernel = gauss_rank_kernel(columns, True)
        expect_eq(len(keys), row["orbit_count"], f"orbit count d={d} R={R}")
        expect_eq(rank, row["rank_Lbar"], f"rank d={d} R={R}")
        expect_eq(len(keys) - rank - 1, row["quotient_dimension"],
                  f"quotient d={d} R={R}")
        expected = 2 * R if d == 1 else (0 if R == 0 else 1)
        expect_eq(row["quotient_dimension"], expected,
                  f"quotient value against the 2R tower, d={d} R={R}")
        if "kernel_report" in row:
            expect_eq(len(kernel), 1 + row["quotient_dimension"],
                      f"kernel basis size d={d} R={R}")
            for vec in kernel:
                sup = {keys[c]: v for c, v in vec.items() if v}
                shapes = {shape_of(k) for k in sup if support(k)}
                if not shapes:
                    continue
                n1 = max(s[0] for s in shapes)
                for key, coeff in sup.items():
                    if not support(key) or shape_of(key)[0] != n1:
                        continue
                    cl = classify(key, 0, shape_of(key))
                    expect(cl["kind"] != "i",
                           f"T1 violated by an exact charge, d={d} R={R}")
                    expect(not cl["thick"],
                           f"T2 violated by an exact charge, d={d} R={R}")
                    expect(cl["kind"] == "iii" or cl["clause_b"],
                           f"Theorem M clause violated by an exact charge, d={d} R={R}")


def check_symbol(stored) -> None:
    for row in stored["leading_symbol_cases"]:
        shape = tuple(row["shape"])
        bond = row["bond"]
        words = enumerate_shape(shape)
        expect_eq(len(words), row["word_count"], f"word count {shape}")
        index = {w: i for i, w in enumerate(words)}
        d = len(shape)

        by_axis = []
        for axis in range(d):
            rows: dict[Word, dict[int, int]] = {}
            for i, w in enumerate(words):
                for target, coeff in growth_images(w, axis, shape, bond).items():
                    rows.setdefault(target, {})[i] = (
                        rows.get(target, {}).get(i, 0) + coeff
                    )
            by_axis.append(rows)

        # every target of the direction-0 symbol has shape sigma + e_0
        target_shape = (shape[0] + 1,) + shape[1:]
        expect(all(shape_of(t) == target_shape for t in by_axis[0]),
               f"direction-0 growth leaves shape sigma+e_1 at {shape}")

        sol = solve_two_entry_rows(list(by_axis[0].values()), len(words))
        expect_eq(sol["max_entries"], row["max_nonzero_entries_per_row"],
                  f"max entries per row {shape}")
        expect(sol["max_entries"] <= 2, f"row with >2 entries at {shape}")
        expect_eq(sol["kernel_dimension"], row["kernel_dimension"],
                  f"symbol kernel dimension {shape}")
        expect_eq(sol["active_kernel_dimension"], row["active_kernel_dimension"],
                  f"symbol active kernel dimension {shape}")
        expect_eq(len(sol["free"]), row["annihilated_word_count"],
                  f"annihilated word count {shape}")

        classes = [classify(w, 0, shape) for w in words]
        t1 = [i for i, c in enumerate(classes) if c["kind"] == "i"]
        t2 = [i for i, c in enumerate(classes) if c["thick"]]
        expect_eq(len(t1), row["T1_words"], f"T1 word count {shape}")
        expect_eq(len(t2), row["T2_words"], f"T2 word count {shape}")
        expect(all(i in sol["forced_zero"] for i in t1),
               f"T1 fails at {shape}: a type-i word is not forced to zero")
        expect(all(i in sol["forced_zero"] for i in t2),
               f"T2 fails at {shape}: a thick-face type-ii word is not forced to zero")
        expect_eq(row["T1_all_forced_zero"], True, f"stored T1 flag {shape}")
        expect_eq(row["T2_all_forced_zero"], True, f"stored T2 flag {shape}")

        free = set(sol["free"])
        surviving = [
            i for i in range(len(words))
            if i not in sol["forced_zero"] and i not in free
        ]
        expect_eq(len(surviving), row["surviving_active_count"],
                  f"surviving active count {shape}")
        bad = [i for i in surviving
               if classes[i]["kind"] != "iii" and not classes[i]["clause_b"]]
        expect(not bad, f"a surviving active word violates clause (b) at {shape}")
        expect_eq(row["clause_b_violations"], [], f"stored clause-(b) violations {shape}")

        combined = []
        for axis in range(d):
            combined.extend(by_axis[axis].values())
        sol_all = solve_two_entry_rows(combined, len(words))
        expect_eq(sol_all["active_kernel_dimension"],
                  row["combined_active_kernel_dimension"],
                  f"combined active kernel dimension {shape}")
        expect_eq(sol_all["kernel_dimension"], row["combined_kernel_dimension"],
                  f"combined kernel dimension {shape}")

        # Lemma B, term by term, on the small shapes: nothing but direction-0
        # growth can reach shape sigma + e_1.
        if len(words) <= 4000:
            for w in words:
                direct = {
                    img: c
                    for img, c in pi_ad(w, 1, bond).items()
                    if support(img) and shape_of(img) == target_shape
                }
                expect_eq(direct, growth_images(w, 0, shape, bond),
                          f"Lemma B identification at {shape}")
            cols = [{} for _ in words]
            row_ids = {k: n for n, k in enumerate(by_axis[0])}
            for key, entries in by_axis[0].items():
                for i, v in entries.items():
                    cols[i][row_ids[key]] = Fraction(v)
            rank, kern = gauss_rank_kernel(cols, True)
            expect_eq(len(words) - rank, sol["kernel_dimension"],
                      f"union-find vs Gaussian elimination at {shape}")
            gauss_forced = {
                i for i in range(len(words))
                if all(v.get(i, Fraction(0)) == 0 for v in kern)
            }
            expect_eq(gauss_forced, sol["forced_zero"],
                      f"union-find vs Gaussian forced-zero set at {shape}")

        # 1D DISCRIMINATOR
        if d == 1:
            expect_eq(len(t2), 0,
                      f"d=1 must have NO thick-face words at {shape}; T2 is vacuous "
                      "on chains and any T2-free strengthening is false")
            expect(sol["active_kernel_dimension"] > 0,
                   f"d=1 active symbol kernel must be nonzero at {shape}")
            expect(len(surviving) > 0,
                   f"d=1 must retain surviving active words at {shape}")
        else:
            expect(len(t2) > 0, f"d>=2 must have thick-face words at {shape}")
        _ = index


def check_local(stored) -> None:
    rng = random.Random(20260817)
    for row in stored["local_row_cases"]:
        shape = tuple(row["shape"])
        d = len(shape)
        found = {label: 0 for label in
                 ("i", "ii-thick", "ii-clause-b", "ii-thin-no-clause-b", "iii")}
        max_row = 0
        for axis in range(d):
            for label in list(found):
                for _ in range(25):
                    built = build_word(shape, axis, label, rng)
                    if built is None:
                        continue
                    found[label] += 1
                    cl = classify(built, axis, shape)
                    sizes = []
                    for target in growth_images(built, axis, shape, 1):
                        pre = preimages(target, axis, shape)
                        sizes.append(len(pre))
                        expect(built in pre,
                               f"word absent from its own row at {shape}")
                    if sizes:
                        max_row = max(max_row, max(sizes))
                    grounded = any(s == 1 for s in sizes)
                    if label == "iii":
                        expect(not sizes, f"type-iii word has a growth image at {shape}")
                    elif label == "ii-clause-b":
                        expect(sizes, f"clause-(b) word has no growth image at {shape}")
                        expect(not grounded,
                               f"clause-(b) word grounded at {shape} (T2 overreach)")
                    else:
                        expect(sizes, f"active word has no growth image at {shape}")
                        expect(grounded,
                               f"{label} word not grounded at {shape}")
        expect(max_row <= 2, f"local row with more than two entries at {shape}")
        expect(row["max_row_size"] <= 2, f"stored max row size {shape}")
        expect_eq(row["violation_count"], 0, f"stored violation count {shape}")
        if d == 1:
            expect_eq(found["ii-thick"], 0,
                      f"d=1 thick-face words must not exist at {shape}")
            expect(found["ii-clause-b"] > 0,
                   f"d=1 clause-(b) words must exist at {shape}")
            expect_eq(row["witnesses_per_class"]["ii-thick"], 0,
                      f"stored d=1 thick count {shape}")
        else:
            expect(found["ii-thick"] > 0,
                   f"d>=2 thick-face words must exist at {shape}")
            expect(row["witnesses_per_class"]["ii-thick"] > 0,
                   f"stored d>=2 thick count {shape}")


def preimages(word: Word, axis: int, shape: tuple[int, ...]) -> dict[Word, int]:
    """All shape-`shape` preimages of a shape-(shape + e_axis) word."""
    a, b = word
    tau = shape_of(word)
    sup = a | b
    out: dict[Word, int] = {}
    for edge, inward in ((tau[axis] - 1, -1), (0, 1)):
        face = [s for s in sup if s[axis] == edge]
        if len(face) != 1:
            continue
        s = face[0]
        if s not in b or s in a:
            continue
        nb = step(s, axis, inward)
        if nb not in a:
            continue
        cand = normalise((a, b ^ {s} ^ {nb}))
        if shape_of(cand) == shape:
            out[cand] = out.get(cand, 0) - 1
    return {k: v for k, v in out.items() if v}


def build_word(shape, axis, want, rng, tries=300):
    d = len(shape)
    n = shape[axis]
    sites = [tuple(c) for c in itertools.product(*(range(s) for s in shape))]
    lo = [s for s in sites if s[axis] == 0]
    hi = [s for s in sites if s[axis] == n - 1]
    mid = [s for s in sites if 0 < s[axis] < n - 1]
    if want == "iii" and n < 3:
        return None
    if want == "ii-thick" and len(lo) < 2:
        return None
    for _ in range(tries):
        a: set = set()
        b: set = set()
        for s in mid:
            if rng.random() < 0.3:
                a.add(s)
            if rng.random() < 0.4:
                b.add(s)
        if want == "i":
            a.add(rng.choice(lo))
            a.add(rng.choice(hi))
        elif want == "iii":
            a -= set(lo) | set(hi)
            if not a:
                a.add(rng.choice(mid))
            b.add(rng.choice(lo))
            b.add(rng.choice(hi))
        else:
            a.add(rng.choice(hi))
            a -= set(lo)
            if want == "ii-thick":
                b.update(rng.sample(lo, 2))
            else:
                r = rng.choice(lo)
                b -= set(lo)
                b.add(r)
                inward = step(r, axis, 1)
                if want == "ii-clause-b":
                    a.add(inward)
                else:
                    a.discard(inward)
        pool = (mid + hi) if want in ("ii-clause-b", "ii-thin-no-clause-b") else sites
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
        word = (frozenset(a), frozenset(b))
        if not support(word) or shape_of(word) != shape:
            continue
        cl = classify(word, axis, shape)
        label = ("i" if cl["kind"] == "i" else "iii" if cl["kind"] == "iii"
                 else "ii-thick" if cl["thick"]
                 else "ii-clause-b" if cl["clause_b"] else "ii-thin-no-clause-b")
        if label == want:
            return word
    return None


def check_single_shape(stored) -> None:
    for row in stored["single_shape_cases"]:
        shape = tuple(row["shape"])
        field = row["field"]
        bond = row["bond"]
        d = len(shape)
        words = enumerate_shape(shape)
        expect_eq(len(words), row["word_count"], f"word count {shape}")
        combined = []
        for axis in range(d):
            rows: dict[Word, dict[int, int]] = {}
            for i, w in enumerate(words):
                for target, coeff in growth_images(w, axis, shape, bond).items():
                    entry = rows.setdefault(target, {})
                    entry[i] = entry.get(i, 0) + coeff
            combined.extend(rows.values())
        sol = solve_two_entry_rows(combined, len(words))
        basis = [{i: Fraction(1)} for i in sol["free"]] + sol["components"]
        expect_eq(len(basis), row["combined_symbol_kernel_dimension"],
                  f"combined symbol kernel dimension {shape}")
        cache: dict[int, dict[Word, int]] = {}
        columns = []
        for vec in basis:
            col: dict[Word, Fraction] = {}
            for i, weight in vec.items():
                img = cache.get(i)
                if img is None:
                    img = cache[i] = pi_ad(words[i], field, bond)
                for key, coeff in img.items():
                    nv = col.get(key, Fraction(0)) + weight * coeff
                    if nv:
                        col[key] = nv
                    else:
                        col.pop(key, None)
            columns.append(col)
        rank, kernel = gauss_rank_kernel(columns, True)
        expect_eq(len(basis) - rank, row["kernel_dimension"],
                  f"single-shape kernel dimension {shape} (f,g)=({field},{bond})")
        if d == 1:
            expect_eq(row["kernel_dimension"], 1,
                      f"d=1 single-shape kernel must be exactly 1 at {shape}: the "
                      "free-fermion current tower is real and no lemma may kill it")
            # the surviving vector must be proportional to j_n
            reference = chain_charge(shape[0])
            vec = kernel[0]
            recovered: dict[Word, Fraction] = {}
            for j, w in vec.items():
                for i, v in basis[j].items():
                    nv = recovered.get(words[i], Fraction(0)) + w * v
                    if nv:
                        recovered[words[i]] = nv
                    else:
                        recovered.pop(words[i], None)
            expect_eq(set(recovered), set(reference),
                      f"d=1 kernel support must be exactly j_{shape[0]}")
            if set(recovered) == set(reference):
                scales = {recovered[k] / reference[k] for k in reference}
                expect_eq(len(scales), 1,
                          f"d=1 kernel vector must be proportional to j_{shape[0]}")
        else:
            expect_eq(row["kernel_dimension"], 0,
                      f"d>=2 single-shape kernel must vanish at {shape}")
        # small shapes: also confirm without the leading-symbol reduction
        if len(words) <= 3600:
            full_cols = [
                {k: Fraction(v) for k, v in pi_ad(w, field, bond).items()}
                for w in words
            ]
            frank, _ = gauss_rank_kernel(full_cols, False)
            expect_eq(len(words) - frank, row["kernel_dimension"],
                      f"unreduced single-shape kernel dimension {shape}")


def main() -> int:
    started = time.monotonic()
    if not RESULT.exists():
        print(f"missing {RESULT}", file=sys.stderr)
        return 1
    stored = json.loads(RESULT.read_text())

    expect_eq(stored["provenance"], "experiments/e127_allr_charges.py", "provenance")
    expect(stored["all_checks_passed"] is True, "stored run reports a failing check")

    check_reference(stored)
    print(f"  reference densities OK ({time.monotonic() - started:.1f}s)", flush=True)
    check_global(stored)
    print(f"  global box quotients OK ({time.monotonic() - started:.1f}s)", flush=True)
    check_symbol(stored)
    print(f"  leading-symbol rows OK ({time.monotonic() - started:.1f}s)", flush=True)
    check_local(stored)
    print(f"  local pointwise rows OK ({time.monotonic() - started:.1f}s)", flush=True)
    check_single_shape(stored)
    print(f"  single-shape kernels OK ({time.monotonic() - started:.1f}s)", flush=True)

    print(f"\n{CHECKS} assertions in {time.monotonic() - started:.1f}s")
    if FAILURES:
        print(f"{len(FAILURES)} FAILURES:")
        for message in FAILURES[:40]:
            print(f"  - {message}")
        return 1
    print("OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
