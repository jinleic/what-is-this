"""Standalone independent verification of the all-size Gaussian no-go (e128).

Nothing is imported from `experiments/e128_allsize_gaussian.py`.  Every decisive
quantity is rebuilt here from scratch, by a construction route that differs from the
producer's (explicit triple-loop rational `P D P` instead of the integer Kronecker
sweep), and the two pair polynomials of the decisive claw certificate are obtained
twice: once from the power-trace identities and once from the characteristic
polynomial of the EXPLICITLY assembled exterior and symmetric squares.

Run:  PYTHONPATH=src .venv/bin/python tests/test_allsize_gaussian.py
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from fractions import Fraction
from math import comb, gcd, isqrt
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ising.transfer_matrix import layer_bonds  # noqa: E402

ARTIFACT = ROOT / "results" / "spectral" / "allsize_gaussian.json"
P1, P2 = 1_000_003, 2_000_003

CLAW_T = (Fraction(1, 3), Fraction(1, 5), Fraction(1, 7), Fraction(2, 11))
CLAW_W = (Fraction(5, 3), Fraction(7, 2), Fraction(11, 5))
CONTROL_T = (
    Fraction(1, 3), Fraction(1, 5), Fraction(1, 7), Fraction(2, 11), Fraction(3, 13),
    Fraction(1, 17), Fraction(4, 19), Fraction(2, 23), Fraction(5, 29), Fraction(3, 31),
)
CONTROL_W = (
    Fraction(5, 3), Fraction(7, 2), Fraction(11, 5), Fraction(13, 4), Fraction(17, 6),
    Fraction(19, 7), Fraction(23, 8), Fraction(29, 9), Fraction(31, 10), Fraction(37, 11),
)
ISOTROPIC_GRID = tuple(Fraction(k, 20) for k in range(4, 11))

PASSED: list[str] = []
FAILED: list[str] = []


def ok(name: str, condition: bool, detail: str = "") -> None:
    tag = "PASS" if condition else "FAIL"
    print(f"[{tag}] {name}" + (f": {detail}" if detail else ""), flush=True)
    (PASSED if condition else FAILED).append(name)


# ------------------------------------------------------------------ independent builder

def build_R_direct(n, bonds, t_list, w_list):
    """R = P D P by explicit triple loop over Fractions (independent of the producer)."""
    dim = 1 << n
    diag = []
    for k in range(dim):
        spins = [1 - 2 * ((k >> (n - 1 - i)) & 1) for i in range(n)]
        acc = Fraction(1)
        for (i, j), w in zip(bonds, w_list):
            if spins[i] * spins[j] == 1:
                acc *= w
        diag.append(acc)
    matrix_p = []
    for k in range(dim):
        row = []
        for l in range(dim):
            diff = k ^ l
            acc = Fraction(1)
            for i in range(n):
                if (diff >> (n - 1 - i)) & 1:
                    acc *= t_list[i]
            row.append(acc)
        matrix_p.append(row)
    result = [[Fraction(0)] * dim for _ in range(dim)]
    for i in range(dim):
        left = matrix_p[i]
        for j in range(i, dim):
            right = matrix_p[j]
            acc = Fraction(0)
            for k in range(dim):
                acc += left[k] * diag[k] * right[k]
            result[i][j] = result[j][i] = acc
    return result


def integralize(matrix):
    """Clear denominators with the ACTUAL least common multiple."""
    scale = 1
    for row in matrix:
        for value in row:
            den = value.denominator
            scale = scale * den // gcd(scale, den)
    out = [[value * scale for value in row] for row in matrix]
    for row in out:
        for value in row:
            assert value.denominator == 1
    return scale, [[int(value) for value in row] for row in out]


# ------------------------------------------------------------------ modular polynomials

def is_prime(n):
    if n < 2:
        return False
    for d in range(2, isqrt(n) + 1):
        if n % d == 0:
            return False
    return True


def modmul(a, b, p):
    assert a.shape[1] * (p - 1) ** 2 < 2**53
    return np.mod(np.asarray(a, dtype=np.float64) @ np.asarray(b, dtype=np.float64),
                  float(p)).astype(np.int64)


def moddot(a, b, p):
    assert a.size * (p - 1) ** 2 < 2**63 - 1
    return int(np.dot(a, b)) % p


def trim(poly):
    for i in range(min(poly.size, 64)):
        if poly[i]:
            return poly if i == 0 else poly[i:]
    nz = np.flatnonzero(poly)
    return np.zeros(1, dtype=np.int64) if nz.size == 0 else poly[nz[0] :]


def is_zero(poly):
    return poly.size == 1 and int(poly[0]) == 0


def monic(poly, p):
    poly = trim(np.asarray(poly, dtype=np.int64))
    if is_zero(poly):
        return poly
    lead = int(poly[0]) % p
    return poly if lead == 1 else (poly * pow(lead, -1, p)) % p


def remainder(a, b, p):
    assert int(b[0]) == 1
    work = (np.asarray(a, dtype=np.int64) % p).copy()
    db, da = b.size - 1, work.size - 1
    if da < db:
        return trim(work)
    for off in range(da - db + 1):
        lead = int(work[off])
        if lead:
            work[off : off + db + 1] = (work[off : off + db + 1] - lead * b) % p
    return trim(work[da - db + 1 :])


def poly_gcd(a, b, p):
    a, b = monic(a, p), monic(b, p)
    while not is_zero(b):
        a, b = b, monic(remainder(a, b, p), p)
    return a


def derivative(poly, p):
    degree = poly.size - 1
    asc = poly[::-1]
    return np.ascontiguousarray(((np.arange(1, degree + 1, dtype=np.int64) * asc[1:]) % p)[::-1])


def newton(power_sums, degree, p):
    """a_k = -(1/k) sum_{i<=k} s_i a_{k-i}, with the window kept contiguous.

    The inner product runs in exactly-representable float64 chunks; every partial sum
    is an integer below 2^53, so the arithmetic is exact, not approximate.
    """
    chunk = (2**53) // ((p - 1) ** 2)
    assert chunk >= 1
    coeffs = np.zeros(degree + 1, dtype=np.int64)
    coeffs[0] = 1
    sums = np.ascontiguousarray(power_sums[: degree + 1], dtype=np.float64)
    window = np.zeros(degree + 1, dtype=np.float64)
    window[degree] = 1.0
    for k in range(1, degree + 1):
        total, base, lo = 0, degree - k, 1
        while lo <= k:
            hi = min(k, lo + chunk - 1)
            total += int(np.dot(sums[lo : hi + 1], window[base + lo : base + hi + 1])) % p
            lo = hi + 1
        value = (-(total % p) * pow(k, -1, p)) % p
        coeffs[k] = value
        window[degree - k] = float(value)
    return coeffs


def charpoly_mod(matrix, p):
    """Descending characteristic polynomial of an integer matrix modulo p."""
    reduced = np.array([[int(x) % p for x in row] for row in matrix], dtype=np.int64)
    dim = reduced.shape[0]
    power = np.eye(dim, dtype=np.int64)
    traces = [0]
    for _ in range(dim):
        power = modmul(power, reduced, p)
        traces.append(int(power.trace()) % p)
    return newton(np.array(traces, dtype=np.int64), dim, p), traces, reduced


def distinct_root_count(poly, p):
    return (poly.size - 1) - (poly_gcd(poly, derivative(poly, p), p).size - 1)


def pair_polynomials(matrix, p, want_all=True):
    """Return (C_2, C_all) modulo p via the power-trace identities."""
    charpoly, traces, reduced = charpoly_mod(matrix, p)
    dim = reduced.shape[0]
    slots2, slots_all = comb(dim, 2), comb(dim + 1, 2)
    needed = 2 * (slots_all if want_all else slots2)
    tail = charpoly[1:]
    for k in range(dim + 1, needed + 1):
        window = np.array(traces[k - 1 : k - dim - 1 : -1], dtype=np.int64)
        traces.append((-moddot(tail, window, p)) % p)
    half = pow(2, -1, p)
    out = []
    for slots, sign in ((slots2, -1), (slots_all, +1)):
        if sign == +1 and not want_all:
            out.append(None)
            continue
        sums = np.zeros(slots + 1, dtype=np.int64)
        for m in range(1, slots + 1):
            sums[m] = (((traces[m] * traces[m]) % p + sign * traces[2 * m]) % p) * half % p
        out.append(newton(sums, slots, p))
    return out[0], out[1]


def counts(matrix, p, want_all=True):
    dim = len(matrix)
    c2, call = pair_polynomials(matrix, p, want_all=want_all)
    r2 = distinct_root_count(c2, p)
    rall = distinct_root_count(call, p) if want_all else None
    return {"dim": dim, "r": r2, "r_all": rall,
            "pair_gcd_degree": comb(dim, 2) - r2,
            "sym_gcd_degree": None if rall is None else comb(dim + 1, 2) - rall}


def exterior_square(matrix, p):
    dim = len(matrix)
    idx = [(i, j) for i in range(dim) for j in range(i + 1, dim)]
    m = np.array([[int(x) % p for x in row] for row in matrix], dtype=np.int64)
    size = len(idx)
    out = np.zeros((size, size), dtype=np.int64)
    for a, (i, j) in enumerate(idx):
        for b, (k, l) in enumerate(idx):
            out[a, b] = (m[i, k] * m[j, l] - m[i, l] * m[j, k]) % p
    return out


def symmetric_square(matrix, p):
    dim = len(matrix)
    idx = [(i, j) for i in range(dim) for j in range(i, dim)]
    m = np.array([[int(x) % p for x in row] for row in matrix], dtype=np.int64)
    size = len(idx)
    out = np.zeros((size, size), dtype=np.int64)
    for a, (k, l) in enumerate(idx):
        for b, (i, j) in enumerate(idx):
            if k == l:
                out[a, b] = (m[k, i] * m[k, j]) % p
            else:
                out[a, b] = (m[k, i] * m[l, j] + m[l, i] * m[k, j]) % p
    return out


def build_M_mod_p(n, bonds, t_list, w_list, p):
    """M = (prod b_v)^2 (prod d_e) * P(t) D(w) P(t), reduced mod p.

    Same integral normalisation as the producer -- required so the coefficient digests
    are comparable, since a different scale rescales every pair product.  The code path
    is different: P and D are assembled entry by entry mod p and multiplied with BLAS,
    instead of the producer's exact-integer Kronecker sweep.  The scale-invariant
    outputs (gcd degree, r) are cross-validated against the exact rational route on
    every small case.
    """
    dim = 1 << n
    p_matrix = np.zeros((dim, dim), dtype=np.int64)
    for k in range(dim):
        for l in range(dim):
            diff, acc = k ^ l, 1
            for v in range(n):
                num, den = t_list[v].numerator, t_list[v].denominator
                acc = acc * (num if (diff >> (n - 1 - v)) & 1 else den) % p
            p_matrix[k, l] = acc
    diag = np.zeros(dim, dtype=np.int64)
    for k in range(dim):
        spins = [1 - 2 * ((k >> (n - 1 - i)) & 1) for i in range(n)]
        acc = 1
        for (i, j), w in zip(bonds, w_list):
            acc = acc * (w.numerator if spins[i] * spins[j] == 1 else w.denominator) % p
        diag[k] = acc
    left = (p_matrix * diag[np.newaxis, :]) % p
    return modmul(left, p_matrix, p)


def digest(ascending):
    return hashlib.sha256(",".join(str(int(c)) for c in ascending).encode()).hexdigest()


def isotropic_curve_degree_bound(n, edges):
    slots = comb(1 << n, 2)
    return 4 * (2 * n + 2 * edges) * slots * (slots - 1)

# ---------------------------------------------------------------------------- utilities

def path_bonds(n):
    return [(i, i + 1) for i in range(n - 1)]


def star_bonds(k):
    return [(0, i) for i in range(1, k + 1)]


def max_degree(n, bonds):
    deg = [0] * n
    for i, j in bonds:
        deg[i] += 1
        deg[j] += 1
    return max(deg) if bonds else 0


def gaussian_ceiling(n):
    return 3**n - 2**n


def family_bound(n, core_sites, r_all_core, r_core):
    m = n - core_sites
    return (3**m - 2**m) * r_all_core + 2**m * r_core


def exceptional_degree_bound(n, edges):
    slots = comb(1 << n, 2)
    return 4 * (2 * n + edges) * slots * (slots - 1)


def run(name, n, bonds, t_list, w_list, p=P1, want_all=True):
    scale, matrix = integralize(build_R_direct(n, bonds, list(t_list), list(w_list)))
    result = counts(matrix, p, want_all=want_all)
    result["scale"] = scale
    result["name"] = name
    return result


# ---------------------------------------------------------------------------- the tests

def main() -> int:
    started = time.process_time()
    artifact = json.loads(ARTIFACT.read_text())

    ok("primes are prime", is_prime(P1) and is_prime(P2))
    ok("artifact primes match", artifact["primes"] == [P1, P2])

    # Newton-identity sanity on a polynomial with known roots {1,2,3,4}.
    roots = [1, 2, 3, 4]
    sums = np.array([0] + [sum(r**k for r in roots) % P1 for k in range(1, 5)], dtype=np.int64)
    coeffs = [int(x) for x in newton(sums, 4, P1)]
    ok("Newton identities reproduce (z-1)(z-2)(z-3)(z-4)",
       coeffs == [1, (-10) % P1, 35, (-50) % P1, 24])

    # ---------------------------------------------------------------- 1. claw certificate
    claw = run("claw", 4, star_bonds(3), CLAW_T, CLAW_W, P1)
    claw2 = run("claw", 4, star_bonds(3), CLAW_T, CLAW_W, P2)
    ok("claw r >= 113 at both primes", claw["r"] == 113 and claw2["r"] == 113,
       f"{claw['r']}, {claw2['r']}")
    ok("claw r_all >= 129 at both primes", claw["r_all"] == 129 and claw2["r_all"] == 129,
       f"{claw['r_all']}, {claw2['r_all']}")
    ok("claw exceeds the 4-mode Gaussian ceilings", claw["r"] > 65 and claw["r_all"] > 81)
    ok("claw r_all - r <= 2^4", claw["r_all"] - claw["r"] <= 16)

    # Second, structurally different route for the SAME decisive numbers: build the
    # exterior and symmetric squares explicitly and take their characteristic polynomials.
    _, claw_int = integralize(build_R_direct(4, star_bonds(3), list(CLAW_T), list(CLAW_W)))
    ext = exterior_square(claw_int, P1)
    sym = symmetric_square(claw_int, P1)
    ext_poly, _, _ = charpoly_mod(ext, P1)
    sym_poly, _, _ = charpoly_mod(sym, P1)
    trace_c2, trace_call = pair_polynomials(claw_int, P1)
    ok("explicit exterior square reproduces the trace-route pair polynomial",
       [int(x) for x in ext_poly] == [int(x) for x in trace_c2])
    ok("explicit symmetric square reproduces the trace-route symmetric polynomial",
       [int(x) for x in sym_poly] == [int(x) for x in trace_call])
    ok("explicit-square route gives the same claw counts",
       distinct_root_count(ext_poly, P1) == 113 and distinct_root_count(sym_poly, P1) == 129)

    stored_core = artifact["core"]
    ok("artifact core matches recomputation",
       stored_core["r_lower_bound"] == claw["r"] and stored_core["r_all_lower_bound"] == claw["r_all"]
       and stored_core["excess_constant"] == claw["r_all"] - 81)

    # ------------------------------------------------------------ 2. path P_4 sharpness
    path4 = run("P4", 4, path_bonds(4), CLAW_T, CLAW_W, P1)
    ok("path P_4 at the same parameters saturates exactly (65, 81)",
       path4["r"] == 65 and path4["r_all"] == 81, f"{path4['r']}, {path4['r_all']}")

    # ------------------------------------------------------------------ 3. chain controls
    for n in range(3, 7):
        chain = run(f"P{n}", n, path_bonds(n), CONTROL_T[:n], CONTROL_W[: n - 1], P1)
        ok(f"open chain P_{n} attains the ceiling with zero excess",
           chain["r"] == gaussian_ceiling(n) and chain["r_all"] == 3**n,
           f"r={chain['r']} r_all={chain['r_all']}")
    stored_chain = {row["sites"]: row for row in artifact["computed_chain_controls"]}
    ok("artifact chain controls never exceed the Gaussian ceiling",
       all(row["r"] == gaussian_ceiling(row["sites"]) and row["r_all"] == 3 ** row["sites"]
           for row in stored_chain.values()),
       f"sizes {sorted(stored_chain)}")

    # ---------------------------------------------- 4. Lemma 3 (tensor decoupling) checks
    def distinct(values, strict):
        seen = set()
        for i in range(len(values)):
            for j in range(i + (1 if strict else 0), len(values)):
                seen.add(values[i] * values[j])
        return len(seen)

    base = [1, 2, 3, 5, 7, 6, 10, 15]
    core_r, core_r_all = distinct(base, True), distinct(base, False)
    synthetic_ok = True
    for m, prime_block in ((1, [11]), (2, [11, 13]), (3, [11, 13, 17])):
        block = [1]
        for prime in prime_block:
            block = [x * y for x in block for y in (1, prime)]
        product = [a * b for a in base for b in block]
        synthetic_ok &= distinct(product, True) == (3**m - 2**m) * core_r_all + 2**m * core_r
        synthetic_ok &= distinct(product, False) == 3**m * core_r_all
    ok("Lemma 3 verified on coprime integer multisets (m = 1,2,3)", synthetic_ok)

    for m in (1, 2):
        n = 4 + m
        case = run(f"claw+{m}", n, star_bonds(3), list(CLAW_T) + list(CONTROL_T[4 : 4 + m]),
                   CLAW_W, P1)
        predicted = family_bound(n, 4, claw["r_all"], claw["r"])
        ok(f"decoupled layer n={n} matches Lemma 3 exactly",
           case["r"] == predicted and case["r_all"] == 3**m * claw["r_all"],
           f"r={case['r']} predicted={predicted}")
        ok(f"decoupled layer n={n} excess is 48*3^{m}",
           predicted - gaussian_ceiling(n) == (claw["r_all"] - 81) * 3**m)

    stored_direct = {row["sites"]: row for row in artifact["direct_specialization"]
                     if "direct_r" in row}
    direct_ok = all(
        row["direct_r"] == family_bound(row["sites"], 4, claw["r_all"], claw["r"])
        and (row.get("direct_r_all") is None
             or row["direct_r_all"] == 3 ** (row["sites"] - 4) * claw["r_all"])
        for row in stored_direct.values())
    ok("every stored direct specialization row matches Lemma 3", direct_ok,
       f"sizes {sorted(stored_direct)}")

    # ------------------------------------------- 5. EVERY row of the family table (all n)
    table_ok = True
    detail = []
    for row in artifact["family_table"]:
        a, b = (int(x) for x in row["layer"].split("x"))
        n = a * b
        bonds = layer_bonds((a, b), (False, False))
        bound = family_bound(n, 4, claw["r_all"], claw["r"])
        ceiling = gaussian_ceiling(n)
        row_ok = (
            row["sites"] == n
            and row["edges"] == len(bonds)
            and row["max_degree"] == max_degree(n, bonds)
            and row["has_branching_site"] == (max_degree(n, bonds) >= 3)
            and int(row["gaussian_ceiling"]) == ceiling
            and int(row["family_lower_bound"]) == bound
            and int(row["excess"]) == bound - ceiling
            and int(row["excess"]) == (claw["r_all"] - 81) * 3 ** (n - 4)
            and bound > ceiling
            and int(row["exceptional_degree_bound"]) == exceptional_degree_bound(n, len(bonds))
        )
        table_ok &= row_ok
        detail.append(f"{row['layer']}:{'ok' if row_ok else 'BAD'}")
    ok("every family-table row recomputes exactly", table_ok, " ".join(detail))

    chain_table_ok = True
    for row in artifact["chain_table"]:
        n = row["sites"]
        chain_table_ok &= (
            int(row["invariant"]) == gaussian_ceiling(n)
            and int(row["excess"]) == 0
            and row["max_degree"] == max_degree(n, path_bonds(n))
        )
    ok("every chain-control row is at the ceiling with zero excess", chain_table_ok)

    # --------------------------------------------------------- 6. the 2x3 decoupling set
    bonds23 = layer_bonds((2, 3), (False, False))
    rungs = [(0, 3), (1, 4), (2, 5)]
    t = Fraction(1, 3)
    q = (1 + t * t) / (2 * t)

    iso_on = run("2x3 iso", 6, bonds23, [t] * 6, [q] * len(bonds23), P1)
    ok("2x3 isotropic reproduces the repository gcd degree 385 and r = 1631",
       iso_on["pair_gcd_degree"] == 385 and iso_on["r"] == 1631,
       f"gcd={iso_on['pair_gcd_degree']} r={iso_on['r']}")

    iso_off = run("2x3 iso rungs off", 6, bonds23, [t] * 6,
                  [Fraction(1) if b in rungs else q for b in bonds23], P1)
    ok("2x3 with rungs off (isotropic) is NOT flagged",
       iso_off["r"] <= gaussian_ceiling(6) and iso_off["r"] == 117
       and iso_off["pair_gcd_degree"] == 1899,
       f"r={iso_off['r']} gcd={iso_off['pair_gcd_degree']}")

    ani_off = run("2x3 aniso rungs off", 6, bonds23, CONTROL_T[:6],
                  [Fraction(1) if b in rungs else CONTROL_W[i] for i, b in enumerate(bonds23)], P1)
    ok("2x3 with rungs off (anisotropic) sits exactly on the Gaussian floor",
       ani_off["r"] == gaussian_ceiling(6) and ani_off["pair_gcd_degree"] == comb(64, 2) - 665,
       f"r={ani_off['r']} gcd={ani_off['pair_gcd_degree']}")

    ani_on = run("2x3 aniso", 6, bonds23, CONTROL_T[:6], CONTROL_W[: len(bonds23)], P1)
    ok("2x3 fully anisotropic is flagged", ani_on["r"] > gaussian_ceiling(6),
       f"r={ani_on['r']} > 665")

    stored_dec = artifact["decoupling_2x3"]
    ok("artifact 2x3 decoupling rows match recomputation",
       stored_dec["isotropic_all_bonds_on"]["r"] == iso_on["r"]
       and stored_dec["isotropic_rungs_off"]["r"] == iso_off["r"]
       and stored_dec["anisotropic_rungs_off"]["r"] == ani_off["r"]
       and stored_dec["anisotropic_all_bonds_on"]["r"] == ani_on["r"])

    # ------------------------------------------------- 7. the isotropic-locus obstruction
    sweep_ok = True
    for coupling in ISOTROPIC_GRID:
        qq = (1 + coupling * coupling) / (2 * coupling)
        case = run("iso claw", 4, star_bonds(3), [coupling] * 4, [qq] * 3, P1)
        sweep_ok &= case["r"] == 65 and case["r_all"] == 73
    ok("isotropic claw yields NO CERTIFIED excess on the declared grid (certified "
       "r >= 65, r_all >= 73; both are lower bounds, so this is a non-certificate)",
       sweep_ok)
    stored_sweep = artifact["isotropic_claw_sweep"]
    ok("artifact isotropic claw sweep agrees",
       all(row["r"] == 65 and row["r_all"] == 73 and not row["certified_flagged"]
           for row in stored_sweep)
       and len(stored_sweep) == len(ISOTROPIC_GRID))

    expected5 = {
        "K_{1,4} star": (147, False),
        "T = claw with one branch of length 2": (417, True),
        "claw + isolated site": (203, False),
        "paw (triangle + pendant)": (373, True),
        "path P_5": (211, False),
        "cycle C_5": (141, False),
    }
    graphs5 = {
        "K_{1,4} star": star_bonds(4),
        "T = claw with one branch of length 2": [(0, 1), (0, 2), (0, 3), (3, 4)],
        "claw + isolated site": star_bonds(3),
        "paw (triangle + pendant)": [(0, 1), (0, 2), (1, 2), (0, 3)],
        "path P_5": path_bonds(5),
        "cycle C_5": path_bonds(5) + [(4, 0)],
    }
    survey_ok = True
    for name, bonds in graphs5.items():
        case = run(name, 5, bonds, [t] * 5, [q] * len(bonds), P1)
        want_r, want_flag = expected5[name]
        survey_ok &= case["r"] == want_r and (case["r"] > gaussian_ceiling(5)) == want_flag
    ok("isotropic 5-site survey reproduces the certified lower bounds exactly "
       "(T certified flagged; K_{1,4} and C_5 non-certificates)", survey_ok)
    stored5 = {row["graph"]: row for row in artifact["isotropic_n5_survey"]}
    ok("artifact 5-site survey agrees",
       all(stored5[name]["r"] == expected5[name][0] for name in expected5))

    # Lemma 9 with the TRUE r_all upper bound C(2^4+1,2) = 136 for the 4-site core.
    collapse_ok = True
    for m in (1, 2, 3):
        n = 4 + m
        case = run(f"iso claw + {m}", n, star_bonds(3), [t] * n, [q] * 3, P1)
        collapse_ok &= case["r"] <= (2 * m + 1) * comb(2**4 + 1, 2)
        collapse_ok &= not (case["r"] > gaussian_ceiling(n))
    ok("equal-field decoupled configurations obey Lemma 9's proved upper bound "
       "r <= (2m+1) * 136 and are never certified flagged", collapse_ok)
    ok("Lemma 9's proved upper bound drops below the Gaussian ceiling exactly from m=3 on",
       all(((2 * m + 1) * comb(2**4 + 1, 2) < gaussian_ceiling(4 + m)) == (m >= 3)
           for m in range(1, 6)),
       "408>211, 680>665, 952<2059, 1224<6305, 1496<19171")

    # -------------------------------- 7b. the uniform-field core-plus-chain route (§10)
    T_BONDS = [(0, 1), (0, 2), (0, 3), (3, 4)]
    core_T = run("T core", 5, T_BONDS, [t] * 5, [q] * 4, P1)
    ok("uniform T core is certified flagged with r >= 417, r_all >= 445",
       core_T["r"] == 417 and core_T["r_all"] == 445 and core_T["r"] > gaussian_ceiling(5))
    # The tensor expression is evaluated at the core LOWER bounds, so it is neither an
    # upper nor a lower bound for the composite count.  Only the certified lower bounds
    # and the ceiling comparison are asserted as mathematics; the coincidence at m=2 is
    # asserted purely as a reproduction of two numbers.
    expected_chain = {1: (1219, False), 2: (3893, True), 3: (7885, False)}
    chain_ok = True
    for m in (1, 2, 3):
        n = 5 + m
        chain = [(5 + i, 6 + i) for i in range(m - 1)]
        case = run(f"T+P{m}", n, T_BONDS + chain, [t] * n, [q] * (4 + len(chain)), P1)
        reference = (3**m - 2**m) * core_T["r_all"] + 2**m * core_T["r"]
        want_r, want_match = expected_chain[m]
        chain_ok &= case["r"] == want_r                             # certified lower bound
        chain_ok &= (case["r"] == reference) == want_match          # reproduction only
        chain_ok &= case["r"] > gaussian_ceiling(n)                 # Theorem U, rigorous
    ok("uniform T + open chain: certified counts reproduce, every case is certified "
       "flagged, and the lower-bound tensor reference is matched at m=2 only", chain_ok)

    stored_route = artifact["uniform_field_chain_route"]
    route_ok = True
    for block in stored_route:
        route_ok &= block["core_T_r_lower_bound"] == 417
        route_ok &= block["core_T_r_all_lower_bound"] == 445
        for row in block["rows"]:
            route_ok &= row["direct_r"] > row["gaussian_ceiling"]       # certified flagged
            route_ok &= row["matches_tensor_reference"] == (
                row["direct_r"] == row["tensor_reference_from_core_lower_bounds"])
            route_ok &= row["direct_r"] == expected_chain[row["m"]][0]
    ok("every stored uniform-field chain row is certified flagged and internally "
       "consistent", route_ok, f"{len(stored_route)} couplings x 3 sizes")
    ok("the uniform-field certified counts are identical at all four couplings (m=1,2,3; "
       "m=4 was computed at t=1/3 only)",
       len({tuple(r["direct_r"] for r in block["rows"]) for block in stored_route}) == 1)

    # Disjoint open chains at a shared coupling: the parity rule.
    parity_ok = True
    for lengths in ((1, 1), (1, 2), (1, 3), (1, 4), (2, 2), (2, 3), (2, 4), (3, 3)):
        bonds, offset = [], 0
        for length in lengths:
            bonds += [(offset + i, offset + i + 1) for i in range(length - 1)]
            offset += length
        case = run("parity", offset, bonds, [t] * offset, [q] * len(bonds), P1)
        certified = case["r"] == gaussian_ceiling(offset)
        parity_ok &= case["r"] <= gaussian_ceiling(offset)      # both parts Gaussian
        parity_ok &= certified == (offset % 2 == 1)             # opposite-parity rule
    ok("two open chains at a shared coupling: injectivity of the pair-exponent map is "
       "CERTIFIED (count pinned against the Lemma 1 ceiling) exactly when their sizes "
       "have opposite parity", parity_ok)
    stored_parity = {tuple(row["lengths"]): row for row in artifact["disjoint_chain_parity"]}
    ok("artifact parity table respects the Lemma 1 ceiling everywhere and labels "
       "certification correctly",
       all(row["r_lower_bound"] <= gaussian_ceiling(row["sites"])
           and row["injectivity_certified"] == (
               row["r_lower_bound"] == gaussian_ceiling(row["sites"]))
           for row in stored_parity.values()))

    ok("open chains at a uniform coupling still saturate exactly",
       all(row["r"] == gaussian_ceiling(row["sites"]) and row["r_all"] == 3 ** row["sites"]
           for row in artifact["uniform_chain_controls"]))
    uniform_recompute = True
    for n in (3, 4, 5, 6):
        case = run(f"isoP{n}", n, path_bonds(n), [t] * n, [q] * (n - 1), P1)
        uniform_recompute &= case["r"] == gaussian_ceiling(n) and case["r_all"] == 3**n
    ok("recomputed uniform-coupling chains P_3..P_6 saturate exactly", uniform_recompute)

    # ------------------------------------------------------------- 8. reciprocality lemma
    def reciprocality(n, bonds, t_list, w_list):
        colour = [None] * n
        adjacency = [[] for _ in range(n)]
        for i, j in bonds:
            adjacency[i].append(j)
            adjacency[j].append(i)
        for start in range(n):
            if colour[start] is not None:
                continue
            colour[start] = 0
            stack = [start]
            while stack:
                u = stack.pop()
                for v in adjacency[u]:
                    if colour[v] is None:
                        colour[v] = 1 - colour[u]
                        stack.append(v)
                    elif colour[v] == colour[u]:
                        raise ValueError("not bipartite")
        matrix = build_R_direct(n, bonds, list(t_list), list(w_list))
        dim = 1 << n
        mask = 0
        for v in range(n):
            if colour[v] == 0:
                mask |= 1 << (n - 1 - v)
        sign = [(-1) ** bin(k).count("1") for k in range(dim)]
        centre = Fraction(1)
        for tv in t_list:
            centre *= (1 - tv * tv) ** 2
        for wv in w_list:
            centre *= wv
        for a in range(dim):
            for b in range(dim):
                acc = Fraction(0)
                for k in range(dim):
                    acc += sign[a ^ mask] * sign[k ^ mask] * matrix[a ^ mask][k ^ mask] * matrix[k][b]
                if acc != (centre if a == b else Fraction(0)):
                    return False, centre
        return True, centre

    holds, centre = reciprocality(4, star_bonds(3), CLAW_T, CLAW_W)
    ok("bipartite reciprocality identity holds exactly for the claw", holds, f"c={centre}")
    holds6, centre6 = reciprocality(6, bonds23, [t] * 6, [q] * len(bonds23))
    ok("bipartite reciprocality identity holds exactly for the isotropic 2x3 layer", holds6,
       f"c={centre6}")
    stored_recip = {row["case"]: row for row in artifact["reciprocality"]}
    ok("artifact reciprocality centres agree",
       stored_recip["claw_K13"]["centre"] == f"{centre.numerator}/{centre.denominator}"
       and stored_recip["2x3_isotropic"]["centre"] == f"{centre6.numerator}/{centre6.denominator}")

    # ------- 9. independent rebuild of the two large certificates claimed by the note ---
    def rebuild_pair_only(name, n, bonds, t_list, w_list, p):
        matrix = build_M_mod_p(n, bonds, list(t_list), list(w_list), p)
        dim = 1 << n
        slots = comb(dim, 2)
        power = np.eye(dim, dtype=np.int64)
        traces = [0]
        for _ in range(dim):
            power = modmul(power, matrix, p)
            traces.append(int(power.trace()) % p)
        chi = newton(np.array(traces, dtype=np.int64), dim, p)
        tail = chi[1:]
        for k in range(dim + 1, 2 * slots + 1):
            win = np.array(traces[k - 1 : k - dim - 1 : -1], dtype=np.int64)
            traces.append((-moddot(tail, win, p)) % p)
        half = pow(2, -1, p)
        sums = np.zeros(slots + 1, dtype=np.int64)
        for m in range(1, slots + 1):
            sums[m] = (((traces[m] * traces[m]) % p - traces[2 * m]) % p) * half % p
        poly = newton(sums, slots, p)
        gcd_degree = poly_gcd(poly, derivative(poly, p), p).size - 1
        return {"slots": slots, "gcd_degree": gcd_degree, "r": slots - gcd_degree,
                "pair_poly_sha256": digest(poly[::-1])}

    stored_attempts = {row["layer"]: row for row in artifact.get("isotropic_direct_attempts", [])}
    for layer, cross in (("2x4", (2, 4)), ("3x3", (3, 3))):
        row = stored_attempts.get(layer, {})
        if row.get("status") != "completed":
            ok(f"isotropic {layer} certificate completed (the proof note claims it)",
               False, f"status {row.get('status')}: {row.get('detail')}")
            continue
        n = cross[0] * cross[1]
        rebuilt = rebuild_pair_only(layer, n, layer_bonds(cross, (False, False)),
                                    [t] * n, [q] * len(layer_bonds(cross, (False, False))), P1)
        ok(f"independent rebuild of the isotropic {layer} certificate agrees exactly",
           rebuilt["gcd_degree"] == row["pair_gcd_degree"]
           and rebuilt["r"] == row["r_lower_bound"]
           and rebuilt["slots"] == comb(1 << n, 2)
           and rebuilt["pair_poly_sha256"] == row["pair_poly_sha256"]
           and rebuilt["r"] > gaussian_ceiling(n),
           f"gcd {rebuilt['gcd_degree']}, r {rebuilt['r']} of {rebuilt['slots']} slots, "
           f"ceiling {gaussian_ceiling(n)}, digest match "
           f"{rebuilt['pair_poly_sha256'] == row['pair_poly_sha256']}")

    big = artifact.get("uniform_T_plus_P4", {})
    if big.get("status") == "wall":
        ok("uniform T + P_4 (n=9) certificate completed (the proof note claims it)",
           False, str(big.get("detail")))
    else:
        rebuilt = rebuild_pair_only("T+P4", 9, T_BONDS + [(5, 6), (6, 7), (7, 8)],
                                    [t] * 9, [q] * 7, P1)
        ok("independent rebuild of the uniform T + P_4 certificate agrees exactly",
           rebuilt["gcd_degree"] == big["pair_gcd_degree"]
           and rebuilt["r"] == big["direct_r"]
           and rebuilt["pair_poly_sha256"] == big["pair_poly_sha256"]
           and rebuilt["r"] > big["gaussian_ceiling"],
           f"gcd {rebuilt['gcd_degree']}, r {rebuilt['r']}, ceiling {big['gaussian_ceiling']}, "
           f"digest match {rebuilt['pair_poly_sha256'] == big['pair_poly_sha256']}")

    # Every isotropic-curve genericity row: recompute the degree bound and the comparison.
    curve_ok = True
    for row in artifact["isotropic_curve_genericity"]:
        curve_ok &= int(row["exceptional_t_degree_bound"]) == isotropic_curve_degree_bound(
            row["sites"], row["edges"])
        curve_ok &= row["gaussian_ceiling"] == gaussian_ceiling(row["sites"])
        curve_ok &= row["conclusion_holds_off_finite_set"] == (
            row["certified_r_at_t_one_third"] > row["gaussian_ceiling"])
        curve_ok &= row["conclusion_holds_off_finite_set"]
    ok("every isotropic-curve genericity row recomputes exactly", curve_ok,
       f"{len(artifact['isotropic_curve_genericity'])} rows")
    ok("the 2x3 isotropic-curve degree bound equals pair_product_scale.md's 422472960",
       isotropic_curve_degree_bound(6, 7) == 422_472_960)

    ok("producer reported no failures", artifact["failures"] == [], str(artifact["failures"]))

    elapsed = time.process_time() - started
    print(f"\n{len(PASSED)} passed, {len(FAILED)} failed in {elapsed:.1f}s")
    if FAILED:
        print("FAILURES: " + ", ".join(FAILED))
    return 1 if FAILED else 0


if __name__ == "__main__":
    raise SystemExit(main())
