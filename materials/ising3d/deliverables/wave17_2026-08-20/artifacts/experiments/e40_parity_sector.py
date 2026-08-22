"""Exact physical-parity sector certificate for the 2x3 Ising layer.

At t = tanh(K*/2) = 1/3 the symmetrized transfer operator is, up to a
positive scalar, an exact rational SPD matrix R.  This script decomposes R
under the physical spin flip P = product X_v and uses exact rational inertia
to test both possible assignments of the P sectors to fermionic even/odd
parity.  No floating-point number participates in a certificate.
"""

from __future__ import annotations

import json
import os
import platform
import signal
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ising.transfer_matrix import layer_bonds  # noqa: E402

def hard_timeout(_signum, _frame):
    raise TimeoutError("NON-DECISIVE: exact parity-sector computation exceeded 300 seconds")


def build_R(n: int, bonds: list[tuple[int, int]], t: Fraction):
    """Return the exact rational symmetrized transfer representative R."""
    q = (1 + t * t) / (2 * t)
    dim = 1 << n
    b = []
    for state in range(dim):
        spins = [1 - 2 * ((state >> (n - 1 - i)) & 1) for i in range(n)]
        b.append(sum(spins[i] * spins[j] for i, j in bonds))
    eps = b[0] % 2
    assert all((x - eps) % 2 == 0 for x in b)
    d = [q ** ((x - eps) // 2) for x in b]
    tp = [t**h for h in range(n + 1)]
    kernel = [[tp[bin(i ^ j).count("1")] for j in range(dim)] for i in range(dim)]
    matrix = [[Fraction(0)] * dim for _ in range(dim)]
    for i in range(dim):
        for j in range(i, dim):
            value = sum((kernel[i][k] * d[k] * kernel[j][k] for k in range(dim)), Fraction(0))
            matrix[i][j] = matrix[j][i] = value
    return matrix, q, eps


def physical_flip_map(n: int) -> list[int]:
    """Computational-basis permutation induced by P = product_v X_v."""
    mask = (1 << n) - 1
    return [mask ^ state for state in range(1 << n)]


def split_physical_parity(R: list[list[Fraction]]):
    """Return rational 32x32 matrices representing R on P=+/- sectors.

    For one representative a from every orbit {a, Pa}, use b_a^s=e_a+s e_Pa.
    The matrix below is (1/2) B_s^T R B_s.  Since B_s^T B_s=2I, this is the
    actual coordinate matrix of R in that unnormalised orthogonal basis.  It is
    also a positive scalar multiple of the congruence B_s^T R B_s, so inertia
    tests are unchanged.
    """
    dim = len(R)
    assert dim and dim & (dim - 1) == 0
    half = dim // 2
    mask = dim - 1
    sectors: dict[int, list[list[Fraction]]] = {}
    for sign in (1, -1):
        block = [[Fraction(0)] * half for _ in range(half)]
        for a in range(half):
            pa = mask ^ a
            for b in range(half):
                pb = mask ^ b
                block[a][b] = (R[a][b] + sign * R[a][pb] + sign * R[pa][b] + R[pa][pb]) / 2
        sectors[sign] = block
    return sectors


def count_at(R: list[list[Fraction]], sigma: Fraction):
    """Exact number of eigenvalues below sigma, or None on a zero pivot."""
    n = len(R)
    a = [[R[i][j] - (sigma if i == j else 0) for j in range(n)] for i in range(n)]
    negative = 0
    for k in range(n):
        pivot = a[k][k]
        if pivot == 0:
            return None
        if pivot < 0:
            negative += 1
        inverse = 1 / pivot
        row = a[k]
        for i in range(k + 1, n):
            factor = a[i][k] * inverse
            if factor:
                target = a[i]
                for j in range(k, n):
                    target[j] -= factor * row[j]
    return negative


def count_inside(R, lo: Fraction, hi: Fraction):
    span = hi - lo
    for k in range(40):
        trial = lo + span / 2 if k == 0 else lo + span / 2 + span * Fraction((-1) ** k, 3 ** (k + 2))
        if lo < trial < hi:
            count = count_at(R, trial)
            if count is not None:
                return trial, count
    raise RuntimeError("no nondegenerate rational shift inside bracket")


def count_outward(R, sigma: Fraction, direction: int):
    count = count_at(R, sigma)
    if count is not None:
        return sigma, count
    step = abs(sigma) if sigma else Fraction(1)
    for k in range(1, 40):
        trial = sigma + direction * step * Fraction(1, 10 ** (20 - k // 2))
        count = count_at(R, trial)
        if count is not None:
            return trial, count
    raise RuntimeError("no nondegenerate outward rational shift")


def isolate(R, index: int, lo: Fraction, hi: Fraction, rel: Fraction):
    """Enclose the zero-based index-th eigenvalue by exact inertia bisection."""
    while not (lo > 0 and hi - lo <= rel * lo):
        trial, count = count_inside(R, lo, hi)
        if count <= index:
            lo = trial
        else:
            hi = trial
    return lo, hi


def fraction_interval(pair: tuple[Fraction, Fraction]):
    return [str(pair[0]), str(pair[1])]


def certify_case(name: str, n: int, bonds: list[tuple[int, int]], t: Fraction, rel_exp: int = 8):
    started = time.time()
    R, q, eps = build_R(n, bonds, t)
    flip = physical_flip_map(n)
    dim = 1 << n
    signed_permutation = (
        sorted(flip) == list(range(dim))
        and all(flip[flip[i]] == i for i in range(dim))
        and all(flip[i] != i for i in range(dim))
    )
    commutes = all(R[flip[i]][flip[j]] == R[i][j] for i in range(dim) for j in range(dim))
    sectors = split_physical_parity(R)
    rel = Fraction(1, 10**rel_exp)

    enclosures: dict[int, list[tuple[Fraction, Fraction]]] = {}
    sector_rows = {}
    for sign in (1, -1):
        block = sectors[sign]
        trace_bound = sum((block[i][i] for i in range(len(block))), Fraction(0)) + 1
        assert count_at(block, Fraction(0)) == 0
        enclosures[sign] = [isolate(block, i, Fraction(0), trace_bound, rel) for i in range(3)]
        sector_rows[str(sign)] = {
            "physical_parity": "+1" if sign == 1 else "-1",
            "dimension": len(block),
            "positive_definite_count_below_zero": count_at(block, Fraction(0)),
            "lowest_three_exact_enclosures": [fraction_interval(x) for x in enclosures[sign]],
            "lowest_three_decimal_display": [[float(x[0]), float(x[1])] for x in enclosures[sign]],
            "lowest_three_disjoint": bool(
                enclosures[sign][0][1] < enclosures[sign][1][0]
                and enclosures[sign][1][1] < enclosures[sign][2][0]
            ),
        }

    assignments = []
    for even_sign, odd_sign in ((1, -1), (-1, 1)):
        even_min = enclosures[even_sign][0]
        odd_min, odd_second = enclosures[odd_sign][:2]
        relevant = (even_min, odd_min, odd_second)
        all_positive = all(interval[0] > 0 for interval in relevant)
        relevant_distinct = all(
            left[1] < right[0] or right[1] < left[0]
            for i, left in enumerate(relevant)
            for right in relevant[i + 1 :]
        )
        no_zero_modes_certified = even_min[1] < odd_min[0]
        smallest_factor_simple_certified = odd_min[1] < odd_second[0]
        gaussian_required_order = (
            no_zero_modes_certified and smallest_factor_simple_certified
        )
        assert all_positive and smallest_factor_simple_certified, (
            "cross-product computation requires positive intervals and distinct odd minima"
        )
        predicted_lo = odd_min[0] * odd_second[0] / even_min[1]
        predicted_hi = odd_min[1] * odd_second[1] / even_min[0]
        shift_lo, count_lo = count_outward(sectors[even_sign], predicted_lo, -1)
        shift_hi, count_hi = count_outward(sectors[even_sign], predicted_hi, 1)
        assignments.append(
            {
                "physical_to_fermionic_assignment": {
                    "+1": "even" if even_sign == 1 else "odd",
                    "-1": "odd" if even_sign == 1 else "even",
                },
                "fermionic_even_physical_sign": even_sign,
                "fermionic_odd_physical_sign": odd_sign,
                "all_three_relevant_values_pairwise_disjoint": relevant_distinct,
                "no_zero_modes_certified_by_even_min_below_odd_min": no_zero_modes_certified,
                "smallest_factor_simple_certified_by_disjoint_odd_minima": smallest_factor_simple_certified,
                "gaussian_required_minimum_order_certified": gaussian_required_order,
                "even_min_enclosure_exact": fraction_interval(even_min),
                "odd_min_enclosure_exact": fraction_interval(odd_min),
                "odd_second_enclosure_exact": fraction_interval(odd_second),
                "forced_even_value_interval_exact": [str(predicted_lo), str(predicted_hi)],
                "forced_even_value_interval_decimal_display": [float(predicted_lo), float(predicted_hi)],
                "forced_interval_relative_width_decimal_display": float(
                    (predicted_hi - predicted_lo) / predicted_lo
                ),
                "window_actually_tested_exact": [str(shift_lo), str(shift_hi)],
                "window_encloses_forced_interval": bool(
                    shift_lo <= predicted_lo and shift_hi >= predicted_hi
                ),
                "target_even_sector_count_below_lo": count_lo,
                "target_even_sector_count_below_hi": count_hi,
                "forced_value_certifiably_absent": bool(count_lo == count_hi),
                "interpretation_if_not_absent": (
                    "Differing counts show only that some eigenvalue lies in the positive-width "
                    "window; they do not prove exact equality with the forced value."
                ),
            }
        )

    return {
        "name": name,
        "n_sites": n,
        "n_bonds": len(bonds),
        "dimension": dim,
        "t_tanh_half_Kstar": str(t),
        "exp_2K": str(q),
        "bond_parity": eps,
        "physical_flip_is_fixed_point_free_signed_permutation": signed_permutation,
        "R_commutes_with_physical_flip_exactly": commutes,
        "sector_basis_entries": "{0,+1,-1}; b_a^s=e_a+s e_{P a}",
        "sector_matrix_convention": "(1/2) B_s^T R B_s, equal to the coordinate restriction because B_s^T B_s=2I",
        "relative_enclosure_width": f"1e-{rel_exp}",
        "sectors": sector_rows,
        "assignments": assignments,
        "both_assignments_absent": all(row["forced_value_certifiably_absent"] for row in assignments),
        "elapsed_seconds": round(time.time() - started, 2),
    }


def main() -> int:
    signal.signal(signal.SIGALRM, hard_timeout)
    signal.alarm(300)
    t = Fraction(1, 3)
    cases = [
        ("1D open chain n=6 (Gaussian control)", 6, [(i, i + 1) for i in range(5)]),
        ("2D open grid 2x3 (3D Ising layer)", 6, list(layer_bonds((2, 3), (False, False)))),
    ]
    rows = []
    for name, n, bonds in cases:
        print(f"certifying {name} ...", flush=True)
        row = certify_case(name, n, bonds, t)
        rows.append(row)
        for assignment in row["assignments"]:
            mapping = assignment["physical_to_fermionic_assignment"]
            print(
                f"  P+->{mapping['+1']}, P-->{mapping['-1']}: counts "
                f"{assignment['target_even_sector_count_below_lo']} -> "
                f"{assignment['target_even_sector_count_below_hi']}; absent="
                f"{assignment['forced_value_certifiably_absent']}"
            )

    control, layer = rows
    checks = [
        {
            "name": "physical_flip_sector_construction_exact",
            "passed": all(
                row["physical_flip_is_fixed_point_free_signed_permutation"]
                and row["R_commutes_with_physical_flip_exactly"]
                and all(row["sectors"][str(s)]["dimension"] == 32 for s in (1, -1))
                for row in rows
            ),
            "detail": "Computed P as bitwise complement, verified a fixed-point-free permutation and exact PR=RP, and constructed two rational 32x32 restrictions in each case.",
        },
        {
            "name": "all_three_relevant_lowest_values_certifiably_distinct",
            "passed": all(
                assignment["all_three_relevant_values_pairwise_disjoint"]
                for row in rows
                for assignment in row["assignments"]
            ),
            "detail": "For both assignment symmetries in both cases, exact rational inertia enclosures make the assigned even minimum and first two assigned odd values pairwise disjoint.",
        },
        {
            "name": "2x3_viable_assignment_factor_hypotheses_certified",
            "passed": bool(
                layer["assignments"][0]["no_zero_modes_certified_by_even_min_below_odd_min"]
                and layer["assignments"][0]["smallest_factor_simple_certified_by_disjoint_odd_minima"]
                and not layer["assignments"][1]["gaussian_required_minimum_order_certified"]
            ),
            "detail": "For P+=even, even_min<odd_min<odd_second certifies no zero modes and a simple smallest factor. The swapped assignment is independently rejected because its assigned odd minimum is below its even minimum.",
        },
        {
            "name": "chain_control_windows_contain_eigenvalues",
            "passed": all(
                not assignment["forced_value_certifiably_absent"]
                for assignment in control["assignments"]
            ),
            "detail": "For both assignment symmetries in the Gaussian n=6 chain control, exact inertia counts differ across the forced-value window. This is consistency only, not exact-presence proof.",
        },
        {
            "name": "2x3_forced_value_absent_in_both_physical_parity_assignments",
            "passed": layer["both_assignments_absent"],
            "detail": "The exact target-sector counts agree across both cross-product windows. For P+=even the lemma makes this value forced; the swapped assignment already fails minimum order, so its empty formal window is auxiliary.",
        },
    ]
    out = {
        "meta": {
            "provenance": "experiments/e40_parity_sector.py",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "python": platform.python_version(),
            "arithmetic": "Exact fractions for matrices, enclosures, interval arithmetic, and Sylvester-inertia counts. Decimal floats are display-only.",
            "scope": "One exact coupling t=1/3, one 2x3 open layer, physical P=product X only.",
        },
        "data": {"tanh_half_Kstar": str(t), "rows": rows},
        "checks": checks,
    }
    os.makedirs(ROOT / "results" / "spectral", exist_ok=True)
    with open(ROOT / "results" / "spectral" / "parity_sector.json", "w", encoding="utf-8") as handle:
        json.dump(out, handle, indent=2)
        handle.write("\n")

    for check in checks:
        print(f"[{'PASS' if check['passed'] else 'FAIL'}] {check['name']}: {check['detail']}")
    failures = [check["name"] for check in checks if not check["passed"]]
    if failures:
        print(f"FAIL: {failures}")
        return 1
    print("PASS: P+=even misses the forced value; the swapped assignment fails minimum order and also has an empty formal window.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
