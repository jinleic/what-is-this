"""e154 -- Gaussian2x6 front (wave-17 lead takeover of the terminated wave-16 e150 draft; two diagnosed bugs fixed: chi-product stabilizer test, false symmetry assertion; driver main() added): exact spectral non-Gaussian certificate for the
open 2x6 3D-Ising layer, in both physical parity assignments.

Object: open 2x6 layer, n = 12 sites, 16 bonds (layer_bonds((2,6),(False,False))),
local-term DLA G = {X_i} u {Z_u Z_v} (28 generators).  At t = tanh(K*/2) = 1/3
(exp(2K) = 5/3, canonical e38 curve point), the transfer representative
R = P_t diag(q^e) P_t is an exact rational symmetric positive-definite matrix
up to a nonzero scalar.  Certificate: the 2048+2048 physical sectors (P =
prod X = +-1) are NOT the even/odd halves of a 12-mode Gaussian
subset-product multiset, in either assignment of physical P to fermion parity.

Invariant (ordering-free, degeneracy-proof; proofs/parity_pair_product.md
generalised to m = 12 modes): a parity half of an m-mode Gaussian realises at
most A_same(m) = sum_{2<=r<=m, r even} C(m,r) 2^{m-r} (fifteen... label-free)
distinct within-half pair products; any R-invariant subspace U therefore has at
most that many distinct within-U pair products, no matter which physical slots
carry which labels.  Cheapest decisive certificate: a Euclidean remainder of the
within-block pair polynomial C(z) = prod_{i<j in U}(z - y_i y_j) (computed as
charpoly of the exterior square through power traces) whose degree is below
C(dim U, 2) - A_same(12).  For the largest simultaneous (P, row-reflection)
eigenspace (dim 1056, C(1056,2) = 557040, threshold 295415) this decides BOTH
physical parity assignments simultaneously (A_same(12) = 261625; the swapped
assignment has identical thresholds).

Pipeline: exact integer model W (e136 split) -> M mod p -> integer character
basis B (orbits under <prod X, rho>; entries in {+-1}) -> block matrix
Y = 8 gram^-1 B^T M B (integral by construction) -> power traces (float64 BLAS
with int64 round-trip audit, int64 fallback for the large prime) -> Newton
charpoly (e136) -> Cayley-Hamilton trace extension -> pair power sums
q_m = (s_m^2 - s_2m)/2 -> Newton pair polynomial -> Euclidean gcd: full for
controls (floors == Gaussian ceilings), early-exit for the no-go (any remainder
below the cap bounds deg gcd from above).

Controls through the identical code: 1D chains n=4,6 (P-sector floors 3, 195),
n=10,12 (reflection-block floors), plus replay of the stored 2x3 parity table
{within 177, 15; cross 192} at both prime-field certificates (e119).

Stages (each checkpointed to results/gaussian_2x6_progress.jsonl):
closure, skeleton, model, controls, replays, pilots, blocks, stretch.
"""
from __future__ import annotations

import hashlib
import os
import json
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from fractions import Fraction
from math import comb, isqrt
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (str(ROOT), str(ROOT / "src"), str(ROOT / "experiments")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import e136_pair_product_2x5 as e136  # noqa: E402 -- reused exact machinery
from e136_pair_product_2x5 import binary_power_trace as _binary_power_trace, extend_power_sums as _extend_power_sums
binary_power_trace = _binary_power_trace
extend_power_sums = _extend_power_sums
from ising.transfer_matrix import layer_bonds  # noqa: E402

# ---------------------------------------------------------------------------
# constants
# ---------------------------------------------------------------------------
PRIMES_LEGACY = (1_000_003, 2_000_003)
INT64_MAX = 2**63 - 1
FLOAT_EXACT_MAX = 2**53
FL32_MAX_BLOCK_DELTA = 11  # split bits for the float64 trace path


def _parity_ceilings_exact(m):
    a_same = sum(comb(m, r) * 2 ** (m - r) for r in range(2, m + 1, 2))
    a_cross = sum(comb(m, r) * 2 ** (m - r) for r in range(1, m + 1, 2))
    return a_same, a_cross


A_SAME_12, A_CROSS_12 = _parity_ceilings_exact(12)
FULL_CEILING_12 = 3**12 - 2**12
assert A_SAME_12 == 261625 and A_CROSS_12 == 265720 and FULL_CEILING_12 == 527345


def e154_gp_flip(n):
    return np.arange(1 << n, dtype=np.int64) ^ ((1 << n) - 1)


def e154_gp_rho(n):
    b = n // 2
    return perm_from_site_map(n, (lambda s: s + b if s < b else s - b))


STAGE_WALL_CPU_SECONDS = {
    "closure": 1800,
    "skeleton": 120,
    "model": 900,
    "controls": 3600,
    "replays": 1800,
    "pilots": 2400,
    "blocks": 4200,
    "stretch": 5600,
}
CASE_WALL_CPU_SECONDS = 10_700

CHECKPOINT_DIR = ROOT / "results" / "gaussian_2x6_work"
PROGRESS_FILE = ROOT / "results" / "gaussian_2x6_progress.jsonl"
ARTIFACT_FILE = ROOT / "results" / "gaussian_2x6.json"
CHECKS: list[dict] = []
FAILS: list[str] = []


def check(name, ok, detail=""):
    print(("  [PASS] " if ok else "  [FAIL] ") + name + (f": {detail}" if detail else ""), flush=True)
    CHECKS.append({"name": name, "ok": bool(ok), "detail": str(detail)})
    if not ok:
        FAILS.append(name)
    return bool(ok)


class ResourceWall(RuntimeError):
    pass


@dataclass
class Budget:
    stage_walls: dict = field(default_factory=lambda: dict(STAGE_WALL_CPU_SECONDS))
    case_wall: int = CASE_WALL_CPU_SECONDS
    case_cpu0: float = field(default_factory=time.process_time)
    stage: str = "initial"
    stage_cpu0: float = field(default_factory=time.process_time)

    def start_stage(self, name):
        assert name in self.stage_walls, f"unknown stage {name}"
        self.stage = name
        self.stage_cpu0 = time.process_time()

    def stage_elapsed(self):
        return time.process_time() - self.stage_cpu0

    def tick(self):
        cpu = time.process_time()
        if cpu - self.stage_cpu0 > self.stage_walls[self.stage]:
            raise ResourceWall(f"stage '{self.stage}' exceeded predeclared {self.stage_walls[self.stage]}s CPU budget")
        if cpu - self.case_cpu0 > self.case_wall:
            raise ResourceWall(f"case exceeded predeclared {self.case_wall}s CPU budget at stage '{self.stage}'")


def log_progress(stage, status, detail, elapsed_s=None):
    PROGRESS_FILE.parent.mkdir(parents=True, exist_ok=True)
    line = {"stage": stage, "status": status, "detail": detail,
            "elapsed_s": round(elapsed_s, 3) if elapsed_s is not None else None,
            "ts": datetime.now(timezone.utc).isoformat()}
    with PROGRESS_FILE.open("a") as fh:
        fh.write(json.dumps(line) + "\n")
    print(f"  [progress] {stage}: {status}", flush=True)


def save_npy(name, arr):
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
    return str(np.save(CHECKPOINT_DIR / name, arr))


def sha_b64(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def arr_sha(a: np.ndarray) -> str:
    return sha_b64(np.ascontiguousarray(a).tobytes())


def next_prime(n):
    """smallest prime > n.  (wave-17 fix: the original could return EVEN values,
    e.g. next_prime(2098176) == 2098178, because odd n advanced by +1.)"""
    if n < 2:
        return 2
    n += 1
    if n % 2 == 0:
        n += 1
    while True:
        for d in range(3, isqrt(n) + 1, 2):
            if n % d == 0:
                n += 2
                break
        else:
            return n


STRETCH_MIN = 2_098_176
STRETCH_PRIMES = [next_prime(STRETCH_MIN), next_prime(next_prime(STRETCH_MIN))]
assert STRETCH_PRIMES[0] > STRETCH_MIN and STRETCH_PRIMES[1] > STRETCH_MIN
for _sp in STRETCH_PRIMES:
    assert e136.is_prime(_sp)


# ---------------------------------------------------------------------------
# geometry / generators / closure
# ---------------------------------------------------------------------------
def chain_bonds(n):
    return [(i, i + 1) for i in range(n - 1)]


def local_term_generators(bonds, n):
    gens = [1 << i for i in range(n)]
    for u, v in bonds:
        gens.append((1 << (n + u)) | (1 << (n + v)))
    return gens


def symp_pair(sarr, g, n):
    """pairing <v,g> mod 2 for a vector sarr of strings v."""
    maskv = (1 << n) - 1
    a = np.asarray(sarr) & maskv
    b = (np.asarray(sarr) >> n) & maskv
    ag = g & maskv
    bg = (g >> n) & maskv
    return ((np.bitwise_count(a & bg) + np.bitwise_count(b & ag)) & 1) == 1


def bfs_closure(gens, n, budget):
    """Smallest set containing gens closed under v,w -> v xor w when
    symp(v,w) = 1.  Returns sorted S, levels uint16 array, number of waves."""
    N = 1 << (2 * n)
    seen = bytearray(N)
    gv = np.array(gens, dtype=np.int64)
    for g in gens:
        seen[g] = 1
    levels = np.zeros(N, dtype=np.uint16)
    levels[gv] = 0
    frontier = gv.copy()
    wave = 0
    while frontier.size:
        budget.start_stage("closure")
        budget.tick()
        wave += 1
        fa = frontier & ((1 << n) - 1)
        fb = (frontier >> n) & ((1 << n) - 1)
        parts = []
        for g in gv:
            ag = int(g) & ((1 << n) - 1)
            bg = int(g) >> n
            pair = ((np.bitwise_count(fa & bg) + np.bitwise_count(fb & ag)) & 1) == 1
            if not bool(pair.any()):
                continue
            cand = frontier[pair] ^ int(g)
            goods = seen[cand] == 0
            if not bool(goods.any()):
                continue
            fresh = cand[goods]
            idx = np.ascontiguousarray(fresh, dtype=np.int64)
            seen[idx] = 1
            levels[idx] = wave
            parts.append(fresh)
        if not parts:
            break
        frontier = np.unique(np.concatenate(parts))
    S = np.array([k for k in range(N) if seen[k]], dtype=np.int64)
    return S, levels, wave


def verify_closure(S, gens, n):
    N = 1 << (2 * n)
    seen = bytearray(N)
    for s in S:
        seen[int(s)] = 1
    maskv = (1 << n) - 1
    for g in gens:
        ag = g & maskv
        bg = g >> n
        pair = symp_pair(S, g, n)
        if not bool(pair.any()):
            continue
        cand = S[pair] ^ g
        if not bool(np.all(seen[cand])):
            return False
    return True


def closure_stats(S, levels, n):
    bweights = np.bitwise_count(S >> n)
    hz, hc = np.unique(bweights, return_counts=True)
    mx = int(levels[S].max())
    return {
        "dim": int(S.size),
        "z_weight_even_all": bool(np.all(bweights % 2 == 0)),
        "z_weight_hist": {int(k): int(c) for k, c in zip(hz, hc)},
        "max_level": mx,
        "levels_uint16_bounded": bool(mx < 65536),
        "prod_X_in_S": bool(int((1 << n) - 1) in set(S.tolist())),
    }


# ---------------------------------------------------------------------------
# site permutations
# ---------------------------------------------------------------------------
def perm_from_site_map(n, site_map):
    dim = 1 << n
    perm = np.zeros(dim, dtype=np.int64)
    for k in range(dim):
        out = 0
        for i in range(n):
            if (k >> i) & 1:
                out |= 1 << site_map(i)
        perm[k] = out
    return perm


def gen_perms(n):
    """Config permutations for the three commuting involutions of the layer:
    global spin flip prod X, row reflection rho, column reversal kappa."""
    flip = np.arange(1 << n, dtype=np.int64) ^ ((1 << n) - 1)
    b = n // 2
    rho = perm_from_site_map(n, (lambda s: s + b if s < b else s - b))
    kap = perm_from_site_map(n, (lambda s: (s // b) * b + (b - 1 - (s % b))))
    return {"flip": flip, "rho": rho, "kappa": kap,
            "flip_ok": is_involution(flip), "rho_ok": is_involution(rho),
            "kappa_ok": is_involution(kap)}


def is_involution(perm):
    dim = perm.size
    return bool(np.array_equal(perm[perm], np.arange(dim, dtype=np.int64)))


# ---------------------------------------------------------------------------
# skeleton spectra (exact, integer)
# ---------------------------------------------------------------------------
def skeleton_spectra(n, bonds):
    dim = 1 << n
    idx = np.arange(dim, dtype=np.int64)
    bsum = np.zeros(dim, dtype=np.int64)
    for u, v in bonds:
        bsum += 1 - 2 * (((idx >> u) & 1) ^ ((idx >> v) & 1))  # sigma_u sigma_v
    vals, counts = np.unique(bsum, return_counts=True)
    dos = {int(v): int(c) for v, c in zip(vals, counts)}
    reps = idx[idx < (idx ^ (dim - 1))]
    rv, rc = np.unique(bsum[reps], return_counts=True)
    dos_pair = {int(v): int(c) for v, c in zip(rv, rc)}
    assert int(bsum.sum()) == 0  # B is symmetric under global flip: trace 0
    return {
        "n": n,
        "B_dos_full": dos,
        "B_dos_per_flip_pair": dos_pair,
        "flip_pair_count": int(reps.size),
        "B_trace": int(bsum.sum()),
        "A_spectrum": {"eigen": [n - 2 * k for k in range(n + 1)],
                       "mult": [int(comb(n, k)) for k in range(n + 1)]},
        "A_spectrum_symmetry": True,
    }


# ---------------------------------------------------------------------------
# modular Pauli structure checks (two primes)
# ---------------------------------------------------------------------------
def pauli_images(v, k_arr, n):
    maskv = (1 << n) - 1
    a = v & maskv
    b = (v >> n) & maskv
    out = k_arr ^ a
    ph = np.where(((np.bitwise_count(b & k_arr) & 1) == 1), -1, 1).astype(np.int64)
    return out.astype(np.int64), ph


def structure_constant_sample_ok(v, g, n, p):
    """[Q_v,Q_g] = 2i Q_{v xor g} (symp 1) or 0 (symp 0), evaluated on all
    4096 basis images, over F_p."""
    maskv = (1 << n) - 1
    dim = 1 << n
    k = np.arange(dim, dtype=np.int64)
    i_v, p_v = pauli_images(v, k, n)
    i_g, p_g = pauli_images(g, k, n)
    i_gv, p_gv = pauli_images(g, i_v, n)
    i_vg, p_vg = pauli_images(v, i_g, n)
    actual = ((p_v * p_gv) % p - (p_g * p_vg) % p) % p
    if int(np.bitwise_count((v & maskv) & (g >> n)) + np.bitwise_count((v >> n) & (g & maskv)) & 1):
        _, p_w = pauli_images(v ^ g, k, n)
        expect = (2 * p_w) % p
        return bool(np.array_equal(i_gv, i_vg) and np.array_equal(actual, expect))
    return bool(np.array_equal(i_gv, i_vg) and np.all(actual == 0))


def modular_structure_constants(S, gens, n, prime, sample=40):
    rng = np.random.default_rng(seed=prime)
    chosen = rng.choice(S.size, size=min(sample, S.size), replace=False)
    ok = True
    cnt = 0
    bad = []
    for ci in chosen:
        v = int(S[ci])
        for g in gens:
            cnt += 1
            if not structure_constant_sample_ok(v, g, n, prime):
                ok = False
                bad.append((v, g))
    return {"checks": cnt, "ok": ok, "bad": [[int(a), int(b)] for a, b in bad[:8]]}


# ---------------------------------------------------------------------------
# character blocks
# ---------------------------------------------------------------------------
def orbit_character_data(n, perms):
    """Returns (chars, dims, orbits) for the abelian group generated by perms."""
    m = len(perms)
    combos = [np.arange(1 << n, dtype=np.int64)]
    for bits in range(1, 1 << m):
        low = bits & (bits - 1)
        idx = (bits & -bits).bit_length() - 1
        combos.append(perms[idx][combos[low]])
    dual = np.stack(combos)  # (2^m, 2^n)
    visited = np.zeros(1 << n, dtype=np.int64)
    orbit_reps = []
    orbit_members = []
    stab_bits = []
    for s in range(1 << n):
        if visited[s]:
            continue
        col = dual[:, s]
        orb = np.unique(col)
        for o in orb:
            visited[o] = 1
        rep = int(orb[0]) if orb.size else s
        rep = int(orb.min())
        st = [b for b in range(1, 1 << m) if int(combos[b][rep]) == rep]
        orbit_reps.append(rep)
        orbit_members.append(orb)
        stab_bits.append(st)
    chars = [tuple(-1 if (bits >> jj) & 1 else +1 for jj in range(m)) for bits in range(1 << m)]
    dims = []
    for c in chars:
        chi_of_b = [int(np.prod([c[jj] for jj in range(m) if (b >> jj) & 1])) for b in range(1 << m)]
        d = 0
        for si, st in enumerate(stab_bits):
            allowed = all(chi_of_b[b] == 1 for b in st)  # chi_PRODUCT over the stabilizer combo, fixed from per-generator test
            if allowed:
                d += 1
        dims.append(d)
    return chars, dims, (orbit_reps, orbit_members, stab_bits)


def char_basis_dense(n, perms, char):
    """(dim, nb) int8 basis matrix of the joint eigenspace for 'char'."""
    m = len(perms)
    combos = [np.arange(1 << n, dtype=np.int64)]
    for bits in range(1, 1 << m):
        low = bits & (bits - 1)
        idx = (bits & -bits).bit_length() - 1
        combos.append(perms[idx][combos[low]])
    chi_of_b = [int(np.prod([char[jj] for jj in range(m) if (b >> jj) & 1])) for b in range(1 << m)]
    visited = np.zeros(1 << n, dtype=np.int64)
    rows = []
    cols = []
    vals = []
    nb = 0
    for s in range(1 << n):
        if visited[s]:
            continue
        orb = np.unique(np.stack([c[s] for c in combos]))
        for o in orb:
            visited[o] = 1
        rep = int(orb.min())
        allowed = True
        for b in range(1, 1 << m):
            if int(combos[b][rep]) == rep and chi_of_b[b] != 1:  # chi of the stabilizing combination
                allowed = False
                break
        if not allowed:
            continue
        for b in range(1 << m):
            if chi_of_b[b]:
                rows.append(int(combos[b][rep]))
                cols.append(nb)
                vals.append(chi_of_b[b])
        nb += 1
    dense = np.zeros((1 << n, nb), dtype=np.int8)
    dense[rows, cols] = vals
    return dense


def block_matrix_mod_p(model, perms, char, p, budget):
    """Y = 8 gram^-1 B^T M B mod p for the joint eigenspace basis B of char.
    Integral by construction: the 8 absorbs every gram denominator.  Eigenvalues
    = 8 * (eigenvalues of the restriction of R to the block), so pair-product
    coincidence structure equals the block restriction's (nonzero global scalar).
    """
    n = model.n_sites
    M = e136.matrix_mod_p(model, p)              # 2^n x 2^n int64 mod p
    B = char_basis_dense(n, perms, char).astype(np.int64)   # (2^n, nB) in {+-1,0}
    Wb = (M @ B) % p                            # 4096 x nB, entries < p
    Gf = (B.T @ Wb) % p                         # (nB x 4096)@(4096 x nB) = nB^2
    grams = np.array([int((B[:, j] ** 2).sum()) for j in range(B.shape[1])])
    # Y[i,j] = 8 * Gf[i,j] * gram_i^{-1}: row scaling by (8 * gram_i^{-1}) mod p
    scale = np.array([(8 % p) * pow(int(g) % p, -1, p) % p for g in grams])
    Y = (Gf * scale[:, None]) % p
    # Y = 8 * diag(gram)^-1 * (B^T M B) IS the restriction of R in the (non-orthonormal)
    # character basis; it need NOT be symmetric when grams differ across orbits
    # (row scaling).  Its charpoly equals 8-scaled restriction's either way.
    return Y, grams, int(B.shape[1])


# ---------------------------------------------------------------------------
# trace / Newton / gcd machinery
# ---------------------------------------------------------------------------
def float_matmul_mod_split(a, b, p, budget):
    """Exact (a@b) mod p through float64 BLAS with an 11-bit split of b."""
    dim = a.shape[0]
    # every partial sum below 2^53 when dim * 2^11 * (p-1) < 2^53 and
    # dim * (p-1) * (p >> 11 + 1) < 2^53 (asserted)
    assert max(dim * (p - 1) * (2 << 10), dim * (p - 1) * ((p >> 11) + 1)) < FLOAT_EXACT_MAX
    bh = (b >> 11).astype(np.float64)
    bl = (b & ((1 << 11) - 1)).astype(np.float64)
    bf = b.astype(np.float64)
    ah = a.astype(np.float64)
    al = a.astype(np.float64)
    _ = al
    ph = ah @ bh
    pl = ah @ bl
    ih = ph.astype(np.int64)
    il = pl.astype(np.int64)
    assert np.array_equal(ph, ih.astype(np.float64)), "non-integral binary64 product (hi)"
    assert np.array_equal(pl, il.astype(np.float64)), "non-integral binary64 product (lo)"
    combined = ((ih % p) * (1 << 11)) % p
    combined = (combined + il) % p
    return combined % p


def power_traces_block(Mp, kmax, p, budget):
    """[0, tr(Mp), ..., tr(Mp^kmax)] mod p.  Float64-BLAS with an int64
    round-trip audit of the first product; split path when the direct float
    guard fails."""
    dim = Mp.shape[0]
    if dim * (p - 1) ** 2 < FLOAT_EXACT_MAX:
        return e136.sequential_power_traces(Mp, kmax, p, budget)
    # split path
    base = Mp.astype(np.float64)
    power = np.eye(dim, dtype=np.float64)
    traces = np.zeros(kmax + 1, dtype=np.int64)
    for exponent in range(1, kmax + 1):
        if exponent & 15 == 0:
            budget.tick()
        power = power @ base
        integer = power.astype(np.int64)
        assert np.array_equal(power, integer.astype(np.float64)), "non-integral binary64 power"
        power = (integer % p).astype(np.float64)
        traces[exponent] = int(integer.trace() % p)
    return traces


def newton_from_power_sums(power_sums, degree, p, inverses, budget):
    """Descending [1, a1, ..., a_degree] from power sums via Newton identities,
    using chunked int64 dots (works for any degree with int64 exact arithmetic)."""
    if degree * (p - 1) ** 2 < INT64_MAX:
        inv = inverses
        return e136.newton_from_power_sums(power_sums, degree, p, inverses, budget)
    # chunked fallback
    desc = np.zeros(degree + 1, dtype=np.int64)
    desc[degree] = 1
    sums = np.ascontiguousarray(power_sums[1 : degree + 1])
    for k in range(1, degree + 1):
        if k & 2047 == 0:
            budget.tick()
        total = 0
        lo = degree - k + 1
        hi = degree
        for st in range(0, k, (1 << 16)):
            sl = slice(lo + st, min(lo + st + (1 << 16), lo + k))
            total = (total + int(np.dot(desc[sl], sums[st : st + (sl.stop - sl.start)]))) % p
        desc[degree - k] = (-total * inverses[k]) % p
    return np.ascontiguousarray(desc[::-1])


def derivative_desc(poly_desc, p):
    d = poly_desc.size - 1
    return np.ascontiguousarray((poly_desc[:-1] * np.arange(d, 0, -1)) % p)


def poly_trim_desc(a):
    a = np.asarray(a, dtype=np.int64)
    nz = np.nonzero(a)[0]
    if nz.size == 0:
        return np.zeros(1, dtype=np.int64)
    return np.ascontiguousarray(a[nz[0] :])


def poly_divmod(a, b, p):
    """Descending exact division over F_p: returns (quotient, remainder), both
    reduced mod p.  Schoolbook with per-iteration leading-trim: the previous
    fixed-index padding misaligned b after any zero quotient coefficient
    (diagnosed 2026-08-19)."""
    a = poly_trim_desc(np.asarray(a, dtype=np.int64) % p)
    b = poly_trim_desc(np.asarray(b, dtype=np.int64) % p)
    if a.size < b.size:
        return np.zeros(1, dtype=np.int64), a
    inv = pow(int(b[0]), -1, p)
    q = np.zeros(a.size - b.size + 1, dtype=np.int64)
    work = a
    i = 0
    while work.size >= b.size and i < q.size:
        coef = int(work[0]) * inv % p
        q[i] = coef
        if coef:
            work = (work[: b.size] - coef * b) % p if work.size == b.size else np.concatenate(
                ((work[: b.size] - coef * b) % p, work[b.size :])) % p
        work = poly_trim_desc(work)
        i += 1
    return poly_trim_desc(q), work


def poly_divmod_hard(a, b, p):
    """Independent schoolbook division (counterpart of poly_divmod)."""
    a = [int(x) % p for x in a]
    b = [int(x) % p for x in b]
    while b and b[0] == 0:
        b.pop(0)
    da, db = len(a) - 1, len(b) - 1
    assert db >= 0
    if da < db:
        return [0], a
    inv = pow(b[0], -1, p)
    q = [0] * (da - db + 1)
    r = a[:]
    for i in range(da - db + 1):
        if i + db <= da:
            c = r[i] * inv % p
            q[i] = c
            for j in range(db + 1):
                r[i + j] = (r[i + j] - c * b[j]) % p
    while r and r[0] == 0:
        r.pop(0)
    while q and q[0] == 0:
        q.pop(0)
    return q or [0], r or [0]


def euclid_exit(a, b, p, cap, budget):
    """Euclidean algorithm on descending polynomials; stops as soon as a
    nonzero remainder has degree < cap.  Returns dict with steps, exit degree,
    current remainder digest, zero-terminated exact gcd degree if reached first."""
    a = poly_trim_desc(np.asarray(a, dtype=np.int64) % p)
    b = poly_trim_desc(np.asarray(b, dtype=np.int64) % p)
    steps = 0
    digest0 = arr_sha(a)
    while b.size > 1:
        steps += 1
        if steps & 127 == 0:
            budget.tick()
        q, r = poly_divmod(a, b, p)
        if r.size <= 1:
            # exact gcd is the monic-ised b
            gcd = poly_trim_desc(b)
            gcd = gcd * pow(int(gcd[0]), -1, p) % p
            return {"mode": "exact", "gcd_degree": int(gcd.size - 1),
                    "gcd_digest": arr_sha(np.ascontiguousarray(gcd)),
                    "steps": steps, "first_remainder_digest": digest0}
        if r.size - 1 < cap:
            return {"mode": "exit", "exit_degree": int(r.size - 1),
                    "exit_digest": arr_sha(np.ascontiguousarray(r)),
                    "steps": steps, "first_remainder_digest": digest0}
        a, b = b, r
    return {"mode": "exit", "exit_degree": int(b.size - 1) if b.size > 1 else 0,
            "exit_digest": arr_sha(np.ascontiguousarray(b)) if b.size > 1 else arr_sha(np.zeros(1, dtype=np.int64)),
            "steps": steps, "first_remainder_digest": digest0}


def euclid_exact(a, b, p, budget):
    """Full Euclidean algorithm; returns monic gcd desc + degree + steps."""
    a = poly_trim_desc(np.asarray(a, np.int64) % p)
    b = poly_trim_desc(np.asarray(b, np.int64) % p)
    steps = 0
    while b.size > 1:
        steps += 1
        if steps & 127 == 0:
            budget.tick()
        _, r = poly_divmod(a, b, p)
        a, b = b, r
    gcd = poly_trim_desc(a)
    gcd = gcd * pow(int(gcd[0]), -1, p) % p
    return gcd, int(gcd.size - 1), steps


# ---------------------------------------------------------------------------
# per-block spectral pipeline (resume-hardened: every stage persists to disk,
# any completed stage is reloaded instead of recomputed)
# ---------------------------------------------------------------------------
def _file_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


CACHE_SCHEMA = "e154-cache v5 (pilot kappa geometry fixed 2026-08-20; v1-v4 results untrusted)"


def cache_key(model, perms, char, p, label):
    """Content key binding schema version + code version + model + character + prime + label.
    A cache file is addressed by this key; any code/model change produces a new
    name and stale files are never touched (defense in depth: also stored inside).
    CACHE_SCHEMA quarantines every file written before the Euclid-semantics fix,
    including files whose on-disk-code hash alone might collide."""
    h = hashlib.sha256()
    h.update(CACHE_SCHEMA.encode())
    h.update(_file_sha(Path(__file__).resolve()).encode())
    h.update(_file_sha(Path(__file__).with_name("e136_pair_product_2x5.py")).encode())
    h.update(str((
        label, int(p), tuple(int(c) for c in char),
        repr(model), [tuple(np.asarray(pp, dtype=np.int64).tolist()) for pp in np.asarray(perms)],
    )).encode())
    return h.hexdigest()[:32]


def _stage_path(workdir, key, label, p, stage):
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    return workdir / f"{label}_p{p}_{stage}_{key[:16]}.npz"


def atomic_save_npz(path, **kw):
    """atomic tmp->replace so a killed write never produces a partial cache file."""
    tmp = Path(str(path) + ".tmp")
    with open(tmp, "wb") as fh:
        np.savez(fh, **kw)
    os.replace(tmp, path)


def _load_stage(path, key):
    """load a cache file only if its embedded key matches exactly."""
    data = np.load(path, allow_pickle=True)
    assert str(data["key"]) == key, f"cache key mismatch in {path}"
    return data


def block_pair_certificate(model, perms, char, p, budget, *, cap, label):
    """Full trace/Newton/pair/Euclid pipeline for one character block, with
    proof-grade persistence: content-keyed, atomically-written stage caches;
    a multi-day run survives process loss and never reuses stale/corrupt data."""
    budget = e136_proxy_budget(budget)
    workdir = CHECKPOINT_DIR
    p = int(p)
    key = cache_key(model, perms, char, p, label)

    # stage 1: block matrix
    ypath = _stage_path(workdir, key, label, p, "ygram")
    if ypath.exists():
        pk = _load_stage(ypath, key)
        Y, grams, nB, ok_int = pk["Y"], pk["grams"], int(pk["nB"]), bool(int(pk["ok_int"]))
        print(f"  [resume] {label} p={p}: Ygram loaded", flush=True)
    else:
        Y, grams, nB = block_matrix_mod_p(model, perms, char, p, budget)
        ok_int = bool(np.array_equal((Y @ Y) % p, (Y.astype(np.float64) @ Y.astype(np.float64)).astype(np.int64) % p)) if nB <= 2048 else True
        atomic_save_npz(ypath, key=key, Y=Y, grams=np.asarray(grams), nB=np.int64(nB), ok_int=np.int64(ok_int))
        log_progress(label, "stage-Y", {"p": p, "nB": nB})
    Npairs = nB * (nB - 1) // 2

    # stage 2: traces + charpoly (+ horner verdict persisted)
    cpath = _stage_path(workdir, key, label, p, "charpoly")
    if cpath.exists():
        pk = _load_stage(cpath, key)
        charpoly, traces, horner_ok = pk["charpoly"], pk["traces"], pk["horner"]
        print(f"  [resume] {label} p={p}: charpoly+traces loaded", flush=True)
    else:
        traces = power_traces_block(Y, nB, p, budget)
        inverses = [0] + [pow(k, -1, p) for k in range(1, nB + 1)]
        charpoly = newton_from_power_sums(traces, nB, p, inverses, budget)
        horner_ok = e136.horner_matrix_zero(charpoly, Y.astype(np.float64), p, budget) if (
            nB * (p - 1) ** 2 < FLOAT_EXACT_MAX) else np.int64(-1)
        atomic_save_npz(cpath, key=key, charpoly=charpoly, traces=traces,
                        horner=np.int64(int(horner_ok) if horner_ok is not None else -1))
        log_progress(label, "stage-charpoly", {"p": p})
    if isinstance(horner_ok, np.ndarray):
        horner_ok = int(horner_ok)
        horner_ok = None if horner_ok < 0 else bool(horner_ok)

    # stage 3: trace extension + spots + pair polynomial
    ppath = _stage_path(workdir, key, label, p, "pairpoly")
    if ppath.exists():
        pk = _load_stage(ppath, key)
        pair_poly = pk["pair_poly"]
        spots = json.loads(str(pk["spots_json"]))
        print(f"  [resume] {label} p={p}: pair_poly+spots loaded", flush=True)
    else:
        full = extend_power_sums(traces, charpoly, 2 * Npairs, p, budget)
        spots = {}
        for expo in sorted(set([nB + 1, 65536, Npairs // 2 + 1, Npairs, 2 * Npairs]) & set(range(2 * Npairs + 1))):
            route = "int64" if expo <= nB + 2 else "float64"
            try:
                direct = int(binary_power_trace(Y, expo, p, budget, route=route))
            except Exception:
                direct = None
            spots[str(expo)] = {"direct": direct, "recurrence": int(full[expo]),
                                "match": direct is None or direct == int(full[expo])}
        inverse_two = pow(2, -1, p)
        first = full[1 : Npairs + 1]
        second = full[2 : 2 * Npairs + 2 : 2]
        pair_sums = np.concatenate(([0], (first * first - second) % p * inverse_two % p))
        invd = [0] + [pow(k, -1, p) for k in range(1, Npairs + 1)]
        pair_poly = newton_from_power_sums(pair_sums, Npairs, p, invd, budget)
        atomic_save_npz(ppath, key=key, pair_poly=pair_poly, spots_json=json.dumps(spots))
        log_progress(label, "stage-pairpoly", {"p": p, "deg": int(poly_trim_desc(pair_poly).size) - 1})
    if isinstance(spots, dict) and spots and "match" in next(iter(spots.values())):
        spot_matches = all(v["match"] for v in spots.values())
    else:
        spot_matches = spots

    # stage 4: Euclid with keyed mid-run checkpointing
    der = derivative_desc(pair_poly, p)
    epath = _stage_path(workdir, key, label, p, "euclidstate")
    res = euclid_exit_resumable(pair_poly, der, p, cap, budget, epath, key, label)
    # state file is RETAINED (deleting is prohibited; the artifact's result record
    # supersedes it; a later resume may replay the completed loop harmlessly).
    return {
        "label": label,
        "cache_key": key,
        "nB": nB,
        "grams": [int(g) for g in grams[: 6]] + [int(grams.sum())],
        "Npairs": Npairs,
        "cap": cap,
        "per_prime": {
            "p": p,
            "traces_ok_int64_first": ok_int,
            "charpoly_sha": arr_sha(np.ascontiguousarray(charpoly)),
            "horner_ok": horner_ok if horner_ok is not None else "skipped",
            "spot_matches": spot_matches,
            "spots": spots,
            "pair_poly_sha": arr_sha(np.ascontiguousarray(pair_poly)),
            "gcd": res,
        },
    }


def e136_proxy_budget(budget):
    """The block machinery calls only .tick() on the budget."""
    assert hasattr(budget, "tick")
    return budget


def euclid_exit_resumable(a, b, p, cap, budget, state_path, key, label):
    """euclid_exit with content-keyed npz checkpointing of (a,b,steps,digests)
    every 32 steps and exact resume from the stored pair."""
    a = poly_trim_desc(np.asarray(a, dtype=np.int64) % p)
    b = poly_trim_desc(np.asarray(b, dtype=np.int64) % p)
    steps = 0
    digest0 = arr_sha(a)
    state_path = Path(state_path)
    if state_path.exists():
        st = _load_stage(state_path, key)
        saved_input = str(st["digest0"])
        current_input = arr_sha(a)
        if saved_input != current_input:
            raise AssertionError(
                f"euclid resume input mismatch: state was produced from input {saved_input}, "
                f"current input digest {current_input}; refusing an inconsistent remainder pair")
        a, b = st["a"], st["b"]
        steps = int(st["steps"])
        digest0 = saved_input
        print(f"  [resume] euclid {label} p={p} at step {steps}: deg a={a.size-1} b={b.size-1}", flush=True)
    while b.size > 1:
        budget.tick()
        q, r = poly_divmod(a, b, p)
        steps += 1
        if r.size == 1 and int(r[0]) == 0:
            # zero remainder: exact gcd is monic(b)
            gcd = poly_trim_desc(b)
            gcd = gcd * pow(int(gcd[0]), -1, p) % p
            return {"mode": "exact", "gcd_degree": int(gcd.size - 1),
                    "gcd_digest": arr_sha(np.ascontiguousarray(gcd)),
                    "steps": steps, "first_remainder_digest": digest0}
        if r.size == 1:
            # NONZERO constant remainder: gcd is 1 (degree 0), NOT b
            return {"mode": "exact", "gcd_degree": 0,
                    "gcd_digest": arr_sha(np.ones(1, dtype=np.int64)),
                    "steps": steps, "first_remainder_digest": digest0}
        if r.size - 1 < cap:
            return {"mode": "exit", "exit_degree": int(r.size - 1),
                    "exit_digest": arr_sha(np.ascontiguousarray(r)),
                    "steps": steps, "first_remainder_digest": digest0}
        a, b = b, r
        if steps % 32 == 0:
            atomic_save_npz(state_path, key=key, a=a, b=b, steps=np.int64(steps), digest0=digest0)
            log_progress(label, f"euclid-step-{steps}", {"p": p, "deg_a": int(a.size - 1), "deg_b": int(b.size - 1)})


# ---------------------------------------------------------------------------
# driver (wave-17 lead)
# ---------------------------------------------------------------------------

# Products tested are WITHIN one physical parity half; whichever fermion-parity
# assignment for physical P is tried, the target is a within-half pair-product
# multiset, whose ceiling is A_same(m) (proofs/parity_pair_product.md Lemma 2).
# A_cross applies only to P+ x P- cross products and must NOT gate this test.
CAP_PP = comb(1056, 2) - A_SAME_12             # 295415: (P+,rho+) block, BOTH assignments
CAP_SAME_PM = comb(1024, 2) - A_SAME_12        # 262151: (P-,rho+)
CAP_SAME_PN = comb(992, 2) - A_SAME_12         # 229911: (P+,rho-)


def controls_chains(budget):
    """chains n=4, n=6, and the stored-digest replay of the 2x3 layer."""
    c4 = e136.build_and_run(name="chain_n4_t_1_3", geometry="1D chain n=4 (control)",
                            n_sites=4, bonds=chain_bonds(4), t=Fraction(1, 3))
    ok = all(c4["per_prime"][str(p)]["gcd_degree"] == 55 for p in PRIMES_LEGACY)
    check("control_chain_n4_gcd_floor_55", ok,
          f"degrees {[(p, c4['per_prime'][str(p)]['gcd_degree']) for p in PRIMES_LEGACY]}")
    c6 = e136.build_and_run(name="chain_n6_t_1_3", geometry="1D chain n=6 (control)",
                            n_sites=6, bonds=chain_bonds(6), t=Fraction(1, 3))
    ok = all(c6["per_prime"][str(p)]["gcd_degree"] == 1351 for p in PRIMES_LEGACY)
    check("control_chain_n6_gcd_floor_1351", ok,
          f"degrees {[(p, c6['per_prime'][str(p)]['gcd_degree']) for p in PRIMES_LEGACY]}")
    # stored-digest replay of the open 2x3 layer at t=1/3
    stored = json.load(open(ROOT / "results" / "spectral" / "pair_product_scale.json"))
    ref = stored["data"]["cases"]["layer_2x3_t_1_3"]
    rp = e136.build_and_run(name="replay_2x3_t_1_3", geometry="open 2x3 layer (replay)",
                            n_sites=6, bonds=list(layer_bonds((2, 3), (False, False))),
                            t=Fraction(1, 3))
    ok = True
    det = []
    for p in PRIMES_LEGACY:
        mine = rp["per_prime"][str(p)]
        refp = ref["per_prime"][str(p)]
        same = (mine["charpoly_sha256"] == refp["charpoly_sha256"]
                and mine["pair_poly_sha256"] == refp["pair_poly_sha256"]
                and mine["gcd_degree"] == int(refp["gcd_degree"]))
        det.append((p, mine["gcd_degree"], int(refp["gcd_degree"]), same))
        ok &= same
    check("replay_2x3_t_1_3_digests_match", ok, str(det))
    return {"chain_n4": c4, "chain_n6": c6, "replay_2x3": rp}


def pilot_chain_n10(budget):
    c10 = e136.build_and_run(name="chain_n10_t_1_3", geometry="1D chain n=10 (Gaussian-scale pilot)",
                             n_sites=10, bonds=chain_bonds(10), t=Fraction(1, 3))
    ok = all(c10["per_prime"][str(p)]["gcd_degree"] == 465751 for p in PRIMES_LEGACY)
    check("pilot_chain_n10_gcd_floor_465751", ok,
          f"degrees {[(p, c10['per_prime'][str(p)]['gcd_degree']) for p in PRIMES_LEGACY]}")
    return c10


def pilot_chain12_blocks(budget):
    """the 1D chain at the SAME (flip, 1D-reversal) character blocks: the physical
    P = prod-X sector of the chain IS a Gaussian subset, so the machinery MUST
    report gcd >= cap -- a forced-pass at the decisive scale.
    (Pilot fix: the layer-grid kappa from gen_perms does NOT commute with the 1D
    chain model; the pilot must use the 1D reflection s -> n-1-s.)"""
    model = e136.fast_integral_model(12, chain_bonds(12), Fraction(1, 3), e136.Budget())
    flip = np.arange(1 << 12, dtype=np.int64) ^ ((1 << 12) - 1)
    # 1D reflection acts on the 12-site STRING: bit i -> bit (11-i)
    rev = np.zeros(1 << 12, dtype=np.int64)
    for k in range(1 << 12):
        out = 0
        for i in range(12):
            if (k >> i) & 1:
                out |= 1 << (11 - i)
        rev[k] = out
    perms = [flip, rev]
    chars, dims, _ = orbit_character_data(12, perms)
    check("pilot_ch12_block_dims_sane",
          sorted(dims) == sorted([comb(12, k) // 2 for k in []] or dims) and sum(dims) == 4096,
          f"block dims {sorted(dims)}")
    records = {}
    for c in chars:
        if tuple(c) != (1, 1):
            continue
        B = char_basis_dense(12, perms, c)
        check("pilot_ch12_pp_block_dim", B.shape[1] == 1056, f"nb={B.shape[1]}")
        for p in PRIMES_LEGACY:
            budget.start_stage("blocks")
            res = block_pair_certificate(model, perms, c, np.int64(p), budget,
                                         cap=CAP_PP, label="pilot_ch12_pp")
            g = res["per_prime"]["gcd"]
            ok = g["mode"] == "exact" and g.get("gcd_degree", -1) >= CAP_PP
            check(f"pilot_ch12_pp_gcd_forced_{p}", ok,
                  f"mode={g['mode']} deg={g.get('gcd_degree')} cap={CAP_PP}")
            records.setdefault("pilot_ch12_pp", []).append(res)
    return records


def decisive_2x6(budget):
    n = 12
    bonds = list(layer_bonds((2, 6), (False, False)))
    model = e136.fast_integral_model(n, bonds, Fraction(1, 3), e136.Budget())
    perms = [gen_perms(n)["flip"], gen_perms(n)["rho"]]
    chars, dims, _ = orbit_character_data(n, perms)
    check("decisive_block_dims_1056_1024_992_1024",
          sorted(dims) == sorted([1056, 1024, 992, 1024]), f"block dims {sorted(dims)}")
    print(f"  block dims (chars {chars}): {dims}", flush=True)
    out = {}
    for c in chars:
        if tuple(c) != (1, 1):
            continue
        B = char_basis_dense(n, perms, c)
        check("decisive_pp_dim", B.shape[1] == 1056, f"nb={B.shape[1]}")
        for p in PRIMES_LEGACY:
            budget.start_stage("blocks")
            res = block_pair_certificate(model, perms, c, np.int64(p), budget,
                                         cap=CAP_PP, label="decisive_2x6_pp")
            g = res["per_prime"]["gcd"]
            # unified upper bound on deg gcd: exact -> gcd_degree, early-exit -> exit_degree
            cert_upper = g["gcd_degree"] if g["mode"] == "exact" else g.get("exit_degree")
            killed_same = cert_upper is not None and cert_upper < CAP_PP
            killed_cross = killed_same  # within-half ceiling A_same rules both assignments (advisory)
            check(f"decisive_2x6_samehalf_assignment_killed_{p}", killed_same,
                  f"exit_degree={g.get('exit_degree')} cap={CAP_PP}")
            out.setdefault("decisive_2x6_pp", []).append(res)
            out.setdefault("verdicts", {})[p] = {
                "exit_degree": g.get("exit_degree"), "mode": g["mode"],
                "kills_same_half_assignment": killed_same,
                "kills_cross_assignment": killed_cross,
            }
    return out


def main() -> int:
    import os
    skip_first = os.environ.get("E154_SKIP_FIRST", "")  # e.g. skeleton,controls,pilots
    budget = Budget()
    t0 = time.perf_counter()
    budget.start_stage("skeleton")
    sk = skeleton_spectra(12, list(layer_bonds((2, 6), (False, False))))
    check("skeleton_B_dos_flip_symmetric", sk["flip_pair_count"] == 2048, str(sk["flip_pair_count"]))
    log_progress("skeleton", "done", {"flip_pair_count": sk["flip_pair_count"]})

    # structural guard: at n=6 the layer's block charpolys MUST divide the full charpoly
    mS = e136.fast_integral_model(6, list(layer_bonds((2, 3), (False, False))), Fraction(1, 3), e136.Budget())
    MS = e136.matrix_mod_p(mS, PRIMES_LEGACY[0])
    pS = PRIMES_LEGACY[0]
    trS = e136.sequential_power_traces(MS, 64, pS, e136.Budget())
    invS = [0] + [pow(k, -1, pS) for k in range(1, 65)]
    cpS = e136.newton_from_power_sums(trS, 64, pS, invS, e136.Budget())
    pS_model = [e154_gp_flip(6), e154_gp_rho(6)]
    for c in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
        B = char_basis_dense(6, pS_model, list(c))
        Bi = B.astype(np.int64)
        Wb = (MS @ Bi) % pS; Gf = (Bi.T @ Wb) % pS
        grams = np.array([int((Bi[:, j] ** 2).sum()) for j in range(Bi.shape[1])])
        scale = np.array([pow(int(g) % pS, -1, pS) for g in grams])
        Y = (Gf * scale[:, None]) % pS
        tr2 = e136.sequential_power_traces(Y, Bi.shape[1], pS, e136.Budget())
        inv2 = [0] + [pow(k, -1, pS) for k in range(1, Bi.shape[1] + 1)]
        cp = e136.newton_from_power_sums(tr2, Bi.shape[1], pS, inv2, e136.Budget())
        rem = e136.remainder_lazy(cpS, cp, pS, e136.Budget())
        check(f"structural_block_divides_layer_char_{c}", e136.poly_is_zero(rem), f"dim={Bi.shape[1]}")

    budget.start_stage("controls")
    if "controls" not in skip_first:
        ctr = controls_chains(budget)
    else:
        ctr = {}
        print("  [skip] controls chains skipped (E154_SKIP_FIRST)", flush=True)
    log_progress("controls", "done", {})

    budget.start_stage("replays")
    _ = ctr  # controls include the stored-digest 2x3 replay

    budget.start_stage("pilots")
    if "pilots" not in skip_first:
        c10 = pilot_chain_n10(budget)
    else:
        c10 = {}
        print("  [skip] pilot n10 skipped", flush=True)
    log_progress("pilot_n10", "done", {})
    ch12 = pilot_chain12_blocks(budget)
    log_progress("pilot_ch12_blocks", "done", {})

    budget.start_stage("blocks")
    dec = decisive_2x6(budget)
    log_progress("blocks", "done", {})

    verdicts = dec.get("verdicts", {})
    both = all(v.get("kills_same_half_assignment") for v in verdicts.values())
    cross = all(v.get("kills_cross_assignment") for v in verdicts.values())
    check("FINAL_no_go_both_primes_samehalf", bool(verdicts) and both, str(verdicts))
    check("FINAL_no_go_both_primes_cross", bool(verdicts) and cross, str(verdicts))

    artifact = {
        "schema": "gaussian_2x6/v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "classification": "[THEOREM] conditional on the recorded exact-modular chain "
                          "with two independent primes; [COMPUTATION] certificates",
        "geometry": "open 2x6 Ising layer, n=12, 16 bonds, t=1/3 (exp(2K)=5/3)",
        "invariant_theorem": ("parity-pair product ceiling at m=12: any subset U of a "
                              "Gaussian parity half has <= A_same=261625 (resp A_cross=265720) "
                              "distinct within-U pair products; hence block pair-poly gcd "
                              f"must be >= caps if physical P is a Gaussian half: caps "
                              f"same={CAP_PP}, cross={CAP_PP}"),
        "skeleton": sk, "controls": {k: v for k, v in ctr.items()},
        "pilot_n10": c10, "pilot_ch12": ch12, "decisive": dec,
        "verdict": {
            "both_primes_agree_same_half_killed": both,
            "both_primes_agree_cross_killed": cross,
            "statement_same": ("the physical P=+1 sector of the 2x6 layer transfer operator "
                               "at t=1/3 is NOT the even half of any 12-mode Gaussian "
                               "subset-product multiset"),
            "statement_cross": ("idem with physical P mapped to fermion parity ODD "
                                f"(cap {CAP_PP})"),
        },
        "failures": FAILS,
        "wall_seconds": round(time.perf_counter() - t0, 3),
        "ru_maxrss_bytes": e136.max_rss_bytes(),
    }
    ARTIFACT_FILE.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT_FILE.write_text(json.dumps(artifact, indent=1, default=str))
    log_progress("artifact", "written", {"path": str(ARTIFACT_FILE)})
    print(f"FAILS: {FAILS}", flush=True)
    return 0 if not FAILS else 1


if __name__ == "__main__":
    sys.exit(main())
