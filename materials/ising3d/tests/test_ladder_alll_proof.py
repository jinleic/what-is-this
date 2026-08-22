#!/usr/bin/env python3
"""Clean-room verifier for experiments/e129_ladder_alll_proof.py.

This file does NOT import the producer.  It rebuilds, from the definitions:

  * the open 2xL ladder and the two generators A_L, B_L;
  * exact integer Pauli bracket arithmetic;
  * the hard-core-boson operators D, F, Ddag *defined directly by their action
    on X-basis configurations* (create a pair on an empty edge / hop / destroy a
    pair), and then checks the depth-3 identities
        D    = (1/32)[A,[A,B]] - (1/8)[A,B]
        Ddag = (1/32)[A,[A,B]] + (1/8)[A,B]
        F    = B - (1/16)[A,[A,B]]
    as operator identities on every basis vector;
  * the full Lie closure of <A,B> at L = 2,3 (dim 11 and 263) and the
    evaluation dimension dim(g_L . psi);
  * the closed form for dim K_L against brute-force orbit counting;
  * every row of the stored certificate table, recomputing the module rank.

Injectivity/hypothesis guard: the module bound is valid ONLY because every
generator S used in the closure satisfies S psi = 0.  The verifier recomputes
S psi for every generator recorded in the artifact recipe and FAILS if any of
them is nonzero; it also runs a negative control showing that the guard fires
on an element that does not annihilate psi.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "results" / "algebra_growth" / "ladder_alll_proof.json"
P = 2_147_483_647

FAILURES: list[str] = []
CHECKS = 0


def check(name: str, ok: bool, detail: str = "") -> None:
    global CHECKS
    CHECKS += 1
    if ok:
        print(f"  ok   {name} {detail}")
    else:
        print(f"  FAIL {name} {detail}")
        FAILURES.append(name)


# --------------------------------------------------------------------------
def edges(L: int) -> list[tuple[int, int]]:
    e = [(2 * i, 2 * i + 1) for i in range(L)]
    for i in range(L - 1):
        e += [(2 * i, 2 * i + 2), (2 * i + 1, 2 * i + 3)]
    return e


# ---------------------------- exact Pauli side ----------------------------
def pb(n: int, u: dict[int, Fraction], v: dict[int, Fraction]) -> dict[int, Fraction]:
    """[u, v] in the Q_(a|b) = X^a Z^b basis, exact rationals."""
    mask = (1 << n) - 1
    out: dict[int, Fraction] = {}
    for p1, c1 in u.items():
        a1, b1 = p1 & mask, p1 >> n
        for p2, c2 in v.items():
            a2, b2 = p2 & mask, p2 >> n
            s1 = bin(b1 & a2).count("1") & 1
            s2 = bin(b2 & a1).count("1") & 1
            if s1 == s2:
                continue
            key = p1 ^ p2
            val = out.get(key, Fraction(0)) + (2 if s1 == 0 else -2) * c1 * c2
            if val:
                out[key] = val
            else:
                out.pop(key, None)
    return out


def comb_pauli(*pairs) -> dict[int, Fraction]:
    out: dict[int, Fraction] = {}
    for coef, vec in pairs:
        for k, c in vec.items():
            val = out.get(k, Fraction(0)) + coef * c
            if val:
                out[k] = val
            else:
                out.pop(k, None)
    return out

def direct_hcb_pauli_pieces(L: int) -> tuple[dict[int, Fraction], dict[int, Fraction], dict[int, Fraction]]:
    """Independently expand D, F, Ddag in Q_(a|b)=X^a Z^b coordinates."""
    n = 2 * L

    def piece(selectors: tuple[tuple[int, int], ...]) -> dict[int, Fraction]:
        out: dict[int, Fraction] = {}
        for u, v in edges(L):
            zmask = ((1 << u) | (1 << v)) << n
            for su, sv in selectors:
                for cu, xu in ((Fraction(1, 2), 0), (Fraction(su, 2), 1 << u)):
                    for cv, xv in ((Fraction(1, 2), 0), (Fraction(sv, 2), 1 << v)):
                        key = zmask | xu | xv
                        out[key] = out.get(key, Fraction(0)) + cu * cv
        return {k: c for k, c in out.items() if c}

    # Q_(a|b) applies Z^b first, so the X projector selects the output
    # occupancy: P^-P^- for pair creation, P^+P^+ for annihilation.
    D = piece(((-1, -1),))
    F = piece(((+1, -1), (-1, +1)))
    Ddag = piece(((+1, +1),))
    return D, F, Ddag


def pauli_to_matrix_action(n: int, vec: dict[int, Fraction], x: int) -> dict[int, Fraction]:
    """Apply sum c Q_(a|b) to the X-basis vector |x>.

    In the X eigenbasis, X^a is diagonal with eigenvalue (-1)^{a.x} and Z^b
    flips exactly the sites of b.  Q_(a|b)|x> = X^a Z^b |x> = (-1)^{a.(x^b)}|x^b>.
    """
    mask = (1 << n) - 1
    out: dict[int, Fraction] = {}
    for code, c in vec.items():
        a, b = code & mask, code >> n
        y = x ^ b
        s = -1 if (bin(a & y).count("1") & 1) else 1
        val = out.get(y, Fraction(0)) + s * c
        if val:
            out[y] = val
        else:
            out.pop(y, None)
    return out


def direct_hcb_action(L: int, kind: str, x: int) -> dict[int, Fraction]:
    """D / F / Ddag defined straight from the configuration picture."""
    out: dict[int, Fraction] = {}
    for (u, v) in edges(L):
        bu, bv = (x >> u) & 1, (x >> v) & 1
        if kind == "D" and bu == 0 and bv == 0:
            y = x | (1 << u) | (1 << v)
        elif kind == "Ddag" and bu == 1 and bv == 1:
            y = x & ~((1 << u) | (1 << v))
        elif kind == "F" and bu != bv:
            y = x ^ (1 << u) ^ (1 << v)
        else:
            continue
        out[y] = out.get(y, Fraction(0)) + 1
    return {k: c for k, c in out.items() if c}


def test_structure_identities(Ls=(2, 3, 4)) -> None:
    print("[1] depth-3 identities for D, F, Ddag")
    for L in Ls:
        n = 2 * L
        A = {1 << v: Fraction(1) for v in range(n)}
        B = {((1 << u) | (1 << v)) << n: Fraction(1) for (u, v) in edges(L)}
        adAB = pb(n, A, B)
        adAAB = pb(n, A, adAB)
        Dform = comb_pauli((Fraction(1, 32), adAAB), (Fraction(-1, 8), adAB))
        Ddform = comb_pauli((Fraction(1, 32), adAAB), (Fraction(1, 8), adAB))
        Fform = comb_pauli((1, B), (Fraction(-1, 16), adAAB))
        Ddirect, Fdirect, Dddirect = direct_hcb_pauli_pieces(L)
        check(f"L={L} D formula in integer Pauli basis", Dform == Ddirect)
        check(f"L={L} Ddag formula in integer Pauli basis", Ddform == Dddirect)
        check(f"L={L} F formula in integer Pauli basis", Fform == Fdirect)
        check(f"L={L} B=D+F+Ddag in integer Pauli basis",
              comb_pauli((1, B), (-1, Ddirect), (-1, Fdirect), (-1, Dddirect)) == {})
        check(f"L={L} ad_A(D)=-4D", pb(n, A, Ddirect) == {k: -4 * c for k, c in Ddirect.items()})
        check(f"L={L} ad_A(F)=0", pb(n, A, Fdirect) == {})
        check(f"L={L} ad_A(Ddag)=+4Ddag",
              pb(n, A, Dddirect) == {k: 4 * c for k, c in Dddirect.items()})

        okD = okDd = okF = okA = True
        for x in range(1 << n):
            if pauli_to_matrix_action(n, Dform, x) != direct_hcb_action(L, "D", x):
                okD = False
            if pauli_to_matrix_action(n, Ddform, x) != direct_hcb_action(L, "Ddag", x):
                okDd = False
            if pauli_to_matrix_action(n, Fform, x) != direct_hcb_action(L, "F", x):
                okF = False
            got = pauli_to_matrix_action(n, A, x)
            want = {x: Fraction(n - 2 * bin(x).count("1"))} if n != 2 * bin(x).count("1") else {}
            if got != want:
                okA = False
        check(f"L={L} D formula agrees with direct configuration action", okD)
        check(f"L={L} Ddag formula agrees with direct configuration action", okDd)
        check(f"L={L} F formula agrees with direct configuration action", okF)
        check(f"L={L} A=2L-2N on every configuration", okA)


# ---------------------------- state-space side ----------------------------
class States:
    def __init__(self, L: int):
        self.L = L
        self.n = 2 * L
        self.N = 1 << self.n
        self.E = edges(L)
        idx = np.arange(self.N, dtype=np.int64)
        self.sel = []
        for (u, v) in self.E:
            bu, bv = (idx >> u) & 1, (idx >> v) & 1
            m = (1 << u) | (1 << v)
            self.sel.append((np.flatnonzero((bu == 1) & (bv == 1)),
                             np.flatnonzero((bu == 0) & (bv == 0)),
                             np.flatnonzero(bu != bv), m))
        self.pc = np.zeros(self.N, dtype=np.int16)
        for b in range(self.n):
            self.pc += ((idx >> b) & 1).astype(np.int16)

    def z(self):
        return np.zeros(self.N, dtype=np.int64)

    def D(self, v):
        o = self.z()
        for d1, _d0, _dm, m in self.sel:
            o[d1] += v[d1 ^ m]
        return o % P

    def Dd(self, v):
        o = self.z()
        for _d1, d0, _dm, m in self.sel:
            o[d0] += v[d0 ^ m]
        return o % P

    def F(self, v):
        o = self.z()
        for _d1, _d0, dm, m in self.sel:
            o[dm] += v[dm ^ m]
        return o % P

    def op(self, c):
        return {"D": self.D, "F": self.F, "d": self.Dd}[c]


def lie(S: States, word: str, v):
    if len(word) == 1:
        return S.op(word)(v)
    g = S.op(word[0])
    return (g(lie(S, word[1:], v)) - lie(S, word[1:], g(v))) % P


DN = {"D": 2, "F": 0, "d": -2}


def all_words(maxlen: int) -> list[str]:
    return ["".join(w) for l in range(1, maxlen + 1) for w in product("DFd", repeat=l)]


class Ech:
    def __init__(self, width: int, folds_signs=None):
        self.rows: dict[int, np.ndarray] = {}
        self.width = width
        self.sg = folds_signs

    def proj(self, v):
        if self.sg is None:
            return v % P
        return ((v * self.sg).reshape(-1, self.width).sum(axis=0)) % P

    def add(self, v) -> bool:
        r = self.proj(v)
        while True:
            nz = np.flatnonzero(r)
            if nz.size == 0:
                return False
            piv = int(nz[-1])
            if piv not in self.rows:
                self.rows[piv] = (r * pow(int(r[piv]), P - 2, P)) % P
                return True
            r = (r - int(r[piv]) * self.rows[piv]) % P

    def __len__(self):
        return len(self.rows)


# ------------------------- full closure at small L -------------------------
def closure_dim_and_eval(L: int) -> tuple[int, int]:
    n = 2 * L
    mask = (1 << n) - 1
    A = {1 << v: 1 for v in range(n)}
    B = {((1 << u) | (1 << v)) << n: 1 for (u, v) in edges(L)}

    def hb(u, v):
        out = {}
        for p1, c1 in u.items():
            a1, b1 = p1 & mask, p1 >> n
            for p2, c2 in v.items():
                a2, b2 = p2 & mask, p2 >> n
                s1 = bin(b1 & a2).count("1") & 1
                s2 = bin(b2 & a1).count("1") & 1
                if s1 == s2:
                    continue
                k = p1 ^ p2
                val = (out.get(k, 0) + (c1 * c2 if s1 == 0 else -c1 * c2)) % P
                if val:
                    out[k] = val
                else:
                    out.pop(k, None)
        return out

    rows: dict[int, dict[int, int]] = {}

    def add(vec):
        v = {k: c % P for k, c in vec.items() if c % P}
        while v:
            piv = max(v)
            if piv not in rows:
                inv = pow(v[piv], P - 2, P)
                rows[piv] = {k: (c * inv) % P for k, c in v.items()}
                return True
            f = v[piv]
            for k, c in rows[piv].items():
                nv = (v.get(k, 0) - f * c) % P
                if nv:
                    v[k] = nv
                else:
                    v.pop(k, None)
        return False

    basis = []
    for g in (A, B):
        if add(g):
            basis.append(g)
    i = 0
    while i < len(basis):
        cur = basis[i]
        for g in (A, B):
            w = hb(g, cur)
            if w and add(w):
                basis.append(w)
        i += 1
    dim = len(rows)

    erows: dict[int, dict[int, int]] = {}

    def eadd(vec):
        v = {k: c % P for k, c in vec.items() if c % P}
        while v:
            piv = max(v)
            if piv not in erows:
                inv = pow(v[piv], P - 2, P)
                erows[piv] = {k: (c * inv) % P for k, c in v.items()}
                return True
            f = v[piv]
            for k, c in erows[piv].items():
                nv = (v.get(k, 0) - f * c) % P
                if nv:
                    v[k] = nv
                else:
                    v.pop(k, None)
        return False

    for vec in basis:
        img: dict[int, int] = {}
        for code, c in vec.items():
            a, b = code & mask, code >> n
            # Q_(a|b)|0> = (-1)^{|a & b|} |b>
            s = -1 if (bin(a & b).count("1") & 1) else 1
            img[b] = (img.get(b, 0) + s * c) % P
        img = {k: c for k, c in img.items() if c}
        if img:
            eadd(img)
    return dim, len(erows)


# ------------------------------ ceiling K_L -------------------------------
def ceiling_closed_form(L: int) -> int:
    n = 2 * L
    total_even = 1 << (n - 1)
    fix_tau = 1 << L
    fix_rho = (1 << L) if L % 2 == 0 else (1 << (L - 1)) * 2
    fix_taurho = 1 << L
    return (total_even + fix_tau + fix_rho + fix_taurho) // 4


def direct_invariant_orbit_basis(L: int) -> tuple[int, dict[str, int]]:
    """Construct the disjoint orbit-sum basis of K_L directly."""
    n = 2 * L

    def tau(x):
        y = 0
        for i in range(L):
            y |= ((x >> (2 * i + 1)) & 1) << (2 * i)
            y |= ((x >> (2 * i)) & 1) << (2 * i + 1)
        return y

    def rho(x):
        y = 0
        for i in range(L):
            j = L - 1 - i
            y |= ((x >> (2 * i)) & 1) << (2 * j)
            y |= ((x >> (2 * i + 1)) & 1) << (2 * j + 1)
        return y

    seen: set[int] = set()
    sector_dims: dict[str, int] = {}
    for x in range(1 << n):
        weight = x.bit_count()
        if weight & 1 or x in seen:
            continue
        orbit = {x, tau(x), rho(x), tau(rho(x))}
        # Each orbit sum is invariant; distinct orbits have disjoint support
        # and therefore give a linearly independent direct basis of K_L.
        if any(y.bit_count() != weight for y in orbit):
            raise AssertionError("ladder symmetry changed particle number")
        seen |= orbit
        sector_dims[str(weight)] = sector_dims.get(str(weight), 0) + 1
    return sum(sector_dims.values()), sector_dims


def ceiling_brute(L: int) -> int:
    return direct_invariant_orbit_basis(L)[0]


# ---------------------- stabiliser guard + module rank ---------------------
def build_generators(S: States, word_len: int, frontier_cap: int = 64):
    """Rebuild the producer's F_p stabiliser recipe without importing it."""
    psi = S.z()
    psi[0] = 1
    words = all_words(word_len)
    images = {w: lie(S, w, psi) for w in words}

    negative_words = [w for w in words if sum(DN[c] for c in w) < 0]
    negative_words.sort(key=len)
    stab_words = list(dict.fromkeys(["d"] + [w for w in negative_words if len(w) <= 3][:12]))

    # Delta N = 0 acts by a scalar on the vacuum sector.
    zero_words = [w for w in words if sum(DN[c] for c in w) == 0]
    alphas = {}
    scalar_ok = True
    for w in zero_words:
        img = images[w]
        alphas[w] = int(img[0]) % P
        t = img.copy()
        t[0] = 0
        if t.any():
            scalar_ok = False
    ref = next((w for w in zero_words if alphas[w]), None)
    zero_combos = []
    for w in zero_words:
        if len(w) > 3:
            continue
        if alphas[w] == 0:
            zero_combos.append({w: 1})
        elif ref is not None and w != ref and len(ref) <= 3:
            f = (alphas[w] * pow(alphas[ref], P - 2, P)) % P
            zero_combos.append({w: 1, ref: (-f) % P})

    # Delta N = +2 kernel combinations, all formed in the exact field F_p.
    plus2 = [w for w in words if sum(DN[c] for c in w) == 2 and images[w].any()]
    rows: dict[int, np.ndarray] = {}
    cpv: dict[int, dict[str, int]] = {}
    kernel: list[dict[str, int]] = []
    for w in plus2:
        r = images[w] % P
        comb = {w: 1}
        placed = False
        while True:
            nz = np.flatnonzero(r)
            if nz.size == 0:
                break
            piv = int(nz[-1])
            if piv not in rows:
                inv = pow(int(r[piv]), P - 2, P)
                rows[piv] = (r * inv) % P
                cpv[piv] = {k: (c * inv) % P for k, c in comb.items()}
                placed = True
                break
            f = int(r[piv])
            r = (r - f * rows[piv]) % P
            for k, c in cpv[piv].items():
                nc = (comb.get(k, 0) - f * c) % P
                if nc:
                    comb[k] = nc
                else:
                    comb.pop(k, None)
        if not placed and comb:
            kernel.append(comb)
    kernel.sort(key=lambda c: (sum(2 ** len(w) for w in c), len(c)))
    raising_combos = kernel[:8]
    zero_combos = zero_combos[:16]

    def serialise(comb: dict[str, int]) -> list[dict[str, int]]:
        return [
            {"word": w, "coefficient_mod_p": int(c)}
            for w, c in sorted(comb.items())
        ]

    recipe = {
        "negative_stabiliser_words": stab_words,
        "raising_kernel_combos_selected": [serialise(c) for c in raising_combos],
        "zero_grade_reference_word": ref,
        "zero_grade_combos_selected": [serialise(c) for c in zero_combos],
        "frontier_cap": frontier_cap,
        "hopping_chain_length": max(4 * (1 << S.L), 64),
    }
    gens = [{"F": 1}] + [{w: 1} for w in stab_words] + raising_combos + zero_combos
    return gens, images, {
        "scalar_ok": scalar_ok,
        "n_kernel": len(kernel),
        "recipe": recipe,
    }


def apply_generator(S: States, gen: dict[str, int], v):
    acc = S.z()
    for w, c in gen.items():
        acc = (acc + c * lie(S, w, v)) % P
    return acc


def stabiliser_guard(S: States, gens) -> tuple[bool, int]:
    """Every generator must annihilate psi; this is the hypothesis of the
    module bound.  Returns (all_ok, number_checked)."""
    psi = S.z()
    psi[0] = 1
    ok = True
    for g in gens:
        if apply_generator(S, g, psi).any():
            ok = False
    return ok, len(gens)


def module_rank(S: States, L: int, gens, images, target: int,
                width_factor: int = 4, seed: int = 20260816,
                max_frontier: int = 24):
    rng = np.random.default_rng(seed)
    width = min(S.N, width_factor * (1 << L))
    w = 1
    while w * 2 <= width:
        w *= 2
    signs = rng.integers(0, 2, size=S.N, dtype=np.int64) * 2 - 1
    ech = Ech(w, signs if w < S.N else None)
    psi = S.z()
    psi[0] = 1
    ech.add(psi)
    queue = []
    for k, v in images.items():
        if v.any() and ech.add(v):
            queue.append(v)
    head = 0
    done = len(ech) >= target
    chain = max(4 * (1 << L), 64)
    while head < len(queue) and not done:
        vec = queue[head]
        head += 1
        cur = vec
        for _ in range(chain):
            cur = S.F(cur)
            if not cur.any():
                break
            if ech.add(cur):
                if len(queue) - head < max_frontier:
                    queue.append(cur)
                if len(ech) >= target:
                    done = True
                    break
        if done:
            break
        for g in gens:
            z = apply_generator(S, g, vec)
            if z.any() and ech.add(z):
                if len(queue) - head < max_frontier:
                    queue.append(z)
                if len(ech) >= target:
                    done = True
                    break
        if head > 4 * max_frontier:
            queue = queue[head:]
            head = 0
    return len(ech)

def independent_krylov_profile(L: int) -> dict:
    """Recompute the selected F-Krylov data in the direct invariant sectors."""
    S = States(L)
    psi = S.z()
    psi[0] = 1
    _, sector_dims = direct_invariant_orbit_basis(L)
    out = {}
    for m in range(1, L + 1):
        k = 2 * m
        word = "D" * m + "F"
        v = lie(S, word, psi)
        target = sector_dims[str(k)]
        full_dim = int(np.count_nonzero(S.pc == k))
        if not v.any():
            out[str(k)] = {
                "seed_word": word,
                "seed_nonzero": False,
                "F_krylov_rank_mod_p1": 0,
                "invariant_sector_dim_K": target,
                "full_particle_sector_dim": full_dim,
                "selected_seed_is_F_cyclic_in_K": False,
            }
            continue
        ech = Ech(S.N)
        cur = v.copy()
        rank = 0
        while ech.add(cur):
            rank += 1
            cur = S.F(cur)
            if not cur.any() or rank > 4000:
                break
        out[str(k)] = {
            "seed_word": word,
            "seed_nonzero": True,
            "F_krylov_rank_mod_p1": rank,
            "invariant_sector_dim_K": target,
            "full_particle_sector_dim": full_dim,
            "selected_seed_is_F_cyclic_in_K": rank == target,
        }
    return out


# ---------------------------------- main ----------------------------------
def main() -> int:
    t0 = time.time()
    raw_art = json.loads(ART.read_text())
    envelope_ok = (
        set(raw_art) == {"provenance", "data", "checks"}
        and isinstance(raw_art["provenance"], dict)
        and isinstance(raw_art["data"], dict)
        and isinstance(raw_art["checks"], list)
    )
    check("artifact uses provenance/data/checks envelope", envelope_ok)
    if not envelope_ok:
        print(f"\n{CHECKS} checks, {len(FAILURES)} failures, {time.time() - t0:.1f}s")
        print("FAIL: " + ", ".join(FAILURES))
        return 1
    provenance = raw_art["provenance"]
    art = raw_art["data"]
    check("artifact records global wall budget and sub-8GB observed RSS",
          provenance.get("global_wall_budget_seconds") == 7_200
          and 0 < provenance.get("peak_rss_bytes_process", 0) < 8 * (1 << 30))
    check("artifact checks envelope contains only passing producer checks",
          bool(raw_art["checks"]) and all(c.get("passed") is True for c in raw_art["checks"]))
    check("artifact records observed phase walls and RSS",
          bool(art.get("resource_measurements"))
          and all("wall_seconds_observed" in m and "rss_peak_bytes_process" in m
                  for m in art["resource_measurements"]))

    test_structure_identities()
    for row in art["T1_structure_identities"]:
        check(f"L={row['L']} producer recorded all exact structure checks",
              row["all_passed"] and all(row["checks"].values()))

    print("[2] direct orbit-sum construction of K_L")
    for L in range(2, 10):
        direct_dim, sector_dims = direct_invariant_orbit_basis(L)
        row = next(r for r in art["T4_ceiling"] if r["L"] == L)
        check(f"L={L} dim K_L closed form = direct invariant orbit basis",
              ceiling_closed_form(L) == direct_dim, f"({direct_dim})")
        check(f"L={L} producer stored direct dim K_L", row["bruteforce_dim_K"] == direct_dim)
        check(f"L={L} producer stored direct K_L sector dimensions",
              row.get("direct_sector_dimensions") == sector_dims)
    for L in range(3, 13):
        check(f"L={L} dim K_L >= 2^(2L-3) >= 2^L",
              ceiling_closed_form(L) >= (1 << (2 * L - 3)) >= (1 << L),
              f"({ceiling_closed_form(L)})")

    print("[3] full Lie closure and evaluation at psi (small L)")
    for L, expect_dim in ((2, 11), (3, 263)):
        dim, ev = closure_dim_and_eval(L)
        check(f"L={L} dim g = {expect_dim}", dim == expect_dim, f"got {dim}")
        check(f"L={L} evaluation bound dim(g.psi) <= dim g", ev <= dim, f"{ev} <= {dim}")
        check(f"L={L} ceiling dim(g.psi) <= dim K_L", ev <= ceiling_closed_form(L),
              f"{ev} <= {ceiling_closed_form(L)}")
        row = next(r for r in art["full_closure_crosscheck"] if r["L"] == L)
        check(f"L={L} producer p1 Lie/evaluation dimensions reproduced",
              row["dim_g_mod_p1"] == dim and row["dim_g_psi_mod_p1"] == ev,
              f"({row['dim_g_mod_p1']}, {row['dim_g_psi_mod_p1']}) vs ({dim}, {ev})")
        check(f"L={L} stored saturation flag agrees with direct evaluation",
              bool(row["saturates_ceiling"]) == (ev == ceiling_closed_form(L)))
    check("L=3 evaluation saturates the ceiling (14 = 14), cross-checked with dim g_3=263",
          closure_dim_and_eval(3) == (263, 14))

    print("[3b] exact selected-seed F-Krylov diagnostics")
    for L in (3, 4, 5):
        check(f"L={L} stored Krylov profile is independently reproduced",
              art["obstruction_krylov"][str(L)] == independent_krylov_profile(L))
    criterion = art["sector_saturation_cyclicity_criterion"]
    check("sector saturation/cyclicity hypothesis is explicitly named",
          criterion["hypothesis_name"] == "SectorSaturationCyclicity(L)"
          and "g_L.psi=K_L" in criterion["hypothesis"])

    print("[4] stabiliser hypothesis guard (this is what makes the bound valid)")
    for L in (3, 4, 5):
        S = States(L)
        gens, images, diag = build_generators(S, art["C5_certificate"][0]["word_len"])
        check(f"L={L} every Delta-N=0 Lie word acts on psi by a scalar", diag["scalar_ok"])
        ok, ngen = stabiliser_guard(S, gens)
        check(f"L={L} all {ngen} module generators annihilate psi", ok)
        # negative control: D does not annihilate psi, so the guard must reject it
        bad = apply_generator(S, {"D": 1}, S.z() + np.eye(1, S.N, 0, dtype=np.int64)[0])
        check(f"L={L} negative control: D psi != 0 (guard would reject D)", bool(bad.any()))
        # the module identity S(X psi) = [S,X] psi used throughout
        psi = S.z()
        psi[0] = 1
        Xpsi = lie(S, "DF", psi)
        lhs = S.F(Xpsi)
        rhs = lie(S, "FDF", psi)
        check(f"L={L} F(X psi) = [F,X] psi for X = [D,F]",
              bool(np.array_equal(lhs % P, rhs % P)))

    print("[5] certificate table: independently rebuild every recorded recipe")
    for row in art["C5_certificate"]:
        L = row["L"]
        if L > int(sys.argv[1]) if len(sys.argv) > 1 else False:
            continue
        check(f"L={L} certificate field is the recorded exact F_p",
              row["certificate_field"] == f"F_{P}")
        check(f"L={L} projection seed and widths are recorded",
              row["projection"]["seed"] == 20260816
              and row["projection"]["width"] == row["projection_width"]
              and row["projection"]["source_width"] == 1 << (2 * L))
        check(f"L={L} stabiliser word list has no duplicate operations",
              len(row["recipe"]["negative_stabiliser_words"])
              == len(set(row["recipe"]["negative_stabiliser_words"])))
        S = States(L)
        gens, images, diag = build_generators(
            S, row["word_len"], row["recipe"]["frontier_cap"]
        )
        check(f"L={L} all Delta-N=0 recipe words act by a scalar on psi", diag["scalar_ok"])
        check(f"L={L} deterministic stabiliser recipe is reproduced",
              diag["recipe"] == row["recipe"])
        check(f"L={L} all selected module generators annihilate psi",
              stabiliser_guard(S, gens)[0] and all(row["facts"].values()))
        rk = module_rank(
            S,
            L,
            gens,
            images,
            row["target_2L"],
            width_factor=max(1, row["projection"]["width"] // (1 << L)),
            seed=row["projection"]["seed"],
            max_frontier=row["recipe"]["frontier_cap"],
        )
        check(f"L={L} recorded modular module rank is reproduced",
              rk == row["module_rank_lower_bound"],
              f"{rk} vs {row['module_rank_lower_bound']}")
        if row["certified"]:
            check(f"L={L} certified module rank >= 2^L = {1 << L}",
                  rk >= row["target_2L"], f"got {rk}")
        else:
            check(f"L={L} non-certificate remains below 2^L = {1 << L}",
                  rk < row["target_2L"], f"got {rk}")
        check(f"L={L} artifact 'certified' flag agrees",
              bool(row["certified"]) == (rk >= row["target_2L"]))
        del S
    l9 = next(r for r in art["C5_certificate"] if r["L"] == 9)
    l9_note = art["l9_small_recipe_noncertificate"]
    check("L=9 non-certificate data records the observed 420/512 shortfall",
          l9_note["observed_mod_p1_lower_rank"] == l9["module_rank_lower_bound"] == 420
          and l9_note["target"] == l9["target_2L"] == 512
          and not l9["certified"])
    check("L=9 observed certificate wall respects its explicit budget",
          l9["wall_seconds_observed"] <= l9["wall_budget_seconds"])
    attempts = art["optional_improvements"]
    multi = art["multivector_evaluation_preflight"]
    check("three-state multi-vector preflight is explicitly non-closing",
          multi["states"] == [
              "all_plus_X_basis_vacuum",
              "all_minus_X_basis_full_configuration",
              "rung_staggered_X_basis: even rungs empty, odd rungs fully occupied",
          ]
          and multi["word_pool"]["candidate_operator_count"] == 120
          and multi["rank_lower_bound_mod_p1"] < multi["target"]
          and not attempts["multi_vector_evaluation"]["closed_L9_gap"])
    check("rung-lex triangular screen is explicitly non-closing",
          attempts["rung_lex_leader_triangularity"]["order"]
          == "(particle_number, rung_adapted_lexicographic)"
          and attempts["rung_lex_leader_triangularity"]["distinct_minimum_support_leaders"]
          < attempts["rung_lex_leader_triangularity"]["target"]
          and not attempts["rung_lex_leader_triangularity"]["closed_L9_gap"])
    check("longer multi-vector preflight is marked unlaunched rather than inferred",
          attempts["unlaunched_longer_word_multivector_preflight"]["launched"] is False)

    print(f"\n{CHECKS} checks, {len(FAILURES)} failures, {time.time() - t0:.1f}s")
    if FAILURES:
        print("FAIL: " + ", ".join(FAILURES))
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
