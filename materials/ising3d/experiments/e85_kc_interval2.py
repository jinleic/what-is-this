#!/usr/bin/env python3
"""Wave-9 exact K_c interval research certificate.

The script has three deliberately separate routes:

1. A wider, exact Simon--Lieb transfer search using the audited full parity
   transfer at cross-sections of 18, 20, and 21 sites.
2. An exact finite-torus linear-program relaxation of the infrared/GKS
   constraint system.  Primal and dual Fraction certificates quantify the
   finite-volume magnetization floor for L=4 and L=6.
3. FKG-free finite SAW counting inequalities.  A direct injection produces a
   strict submultiplicative refinement; finite Kesten-ratio instances are then
   certified where Fekete plus the audited c_36 count suffices.

No decision depends on a floating-point critical-coupling benchmark.  Float64
is used only to locate Simon--Lieb roots before an exact integer residual is
computed at a rational candidate.
"""
from __future__ import annotations

from datetime import datetime, timezone
from decimal import Decimal, ROUND_FLOOR, localcontext
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
import time
from typing import Sequence

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ising.rigorous_bounds.simon_lieb import (  # noqa: E402
    atanh_rational_interval,
    enumerate_pq_polynomials,
    exact_criterion_residual,
    exact_pq_at_rational,
    solve_box_root_float,
)


SCRIPT = "experiments/e85_kc_interval2.py"
RESULT_PATH = ROOT / "results" / "bounds" / "kc_interval2.json"
UPPER_INFRARED_PATH = ROOT / "results" / "bounds" / "upper_infrared.json"
SAW_SOURCE = ROOT / "sources" / "fulltext" / "schram_barkema_bisseling2011.pdf"
DPS = 90
REPORT_PLACES = 40
C36 = 2941370856334701726560670
C36_SHA256 = "898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12"
INCUMBENT_K_LOWER = "0.2122119011661678393310862783954278184914"
CHALLENGE = (23, 110)

# This fixed certificate family is declared before its locators run.  It
# deliberately samples both reachable wider cross-sections (18 and 20 sites)
# and short/long prisms; every member is strictly beyond the prior <=16-site
# audit.
ROUTE1_SHAPES: tuple[tuple[int, int, int], ...] = (
    (3, 6, 4),
    (3, 6, 12),
    (4, 5, 4),
    (4, 5, 16),
)
ROUTE1_SMALL_DENOMINATOR = 10_000
RATIONAL_COSINES: dict[int, tuple[Fraction, ...]] = {
    4: (Fraction(1), Fraction(0), Fraction(-1), Fraction(0)),
    6: (
        Fraction(1),
        Fraction(1, 2),
        Fraction(-1, 2),
        Fraction(-1),
        Fraction(-1, 2),
        Fraction(1, 2),
    ),
}
LP_SIDES = (4, 6)
LP_POINTS = (Fraction(1, 5), Fraction(21, 100), Fraction(11, 50), Fraction(6, 25))
EXTERNAL_SAW_COUNTS = {
    30: 270569905525454674614,
    31: 1274191064726416905966,
    32: 5997359460809616886494,
    33: 28233744272563685150118,
    34: 132853629626823234210582,
    35: 625248129452557974777990,
    36: C36,
}


def _check(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    """Append a machine-readable check and fail at the first false assertion."""

    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    if not passed:
        raise AssertionError(f"{name}: {detail}")
    print(f"PASS {name}: {detail}")


def _fingerprint(value: int) -> dict[str, object]:
    magnitude = abs(value)
    encoded = magnitude.to_bytes(max(1, (magnitude.bit_length() + 7) // 8), "big")
    return {
        "sign": (value > 0) - (value < 0),
        "bit_length": magnitude.bit_length(),
        "sha256_magnitude_big_endian": hashlib.sha256(encoded).hexdigest(),
    }


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _parse_fraction(value: str | int) -> Fraction:
    return Fraction(value)


def _floor_decimal(value: str, places: int = REPORT_PLACES) -> str:
    with localcontext() as context:
        context.prec = max(DPS + 20, len(value) + 10)
        quantum = Decimal(1).scaleb(-places)
        return format(Decimal(value).quantize(quantum, rounding=ROUND_FLOOR), "f")


def _atanh_record(p: int, q: int) -> tuple[list[str], str]:
    interval = atanh_rational_interval(p, q, dps=DPS)
    return list(interval), _floor_decimal(interval[0])


def _box_boundary_multiplicities(shape: tuple[int, int, int]) -> tuple[int, ...]:
    """q_B(x)=6-deg_B(x), in z,y,x order with x the fast bit coordinate."""

    nx, ny, nz = shape
    values: list[int] = []
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                degree = (
                    int(x > 0)
                    + int(x + 1 < nx)
                    + int(y > 0)
                    + int(y + 1 < ny)
                    + int(z > 0)
                    + int(z + 1 < nz)
                )
                values.append(6 - degree)
    return tuple(values)


def _box_edge_count(shape: tuple[int, int, int]) -> int:
    nx, ny, nz = shape
    return (nx - 1) * ny * nz + nx * (ny - 1) * nz + nx * ny * (nz - 1)


# ---------------------------------------------------------------------------
# Route 1: full exact Simon--Lieb transfer beyond 16 cross-section sites
# ---------------------------------------------------------------------------


def _reduced_safe_point(
    shape: tuple[int, int, int], locator: float, denominator: int
) -> tuple[int, int, int]:
    """Turn a non-authoritative locator into a decided exact-safe rational."""

    raw_numerator = max(1, math.floor(locator * denominator))
    while raw_numerator:
        divisor = math.gcd(raw_numerator, denominator)
        numerator, reduced_denominator = raw_numerator // divisor, denominator // divisor
        residual = exact_criterion_residual(shape, numerator, reduced_denominator)
        if residual < 0:
            return numerator, reduced_denominator, residual
        raw_numerator -= 1
    raise AssertionError(f"no positive safe point proposed by locator for {shape}")


def _route1_controls(checks: list[dict[str, object]]) -> None:
    """Exercise both independent exact paths supplied by the stable API."""

    p_poly, q_poly = enumerate_pq_polynomials((2, 2, 2))
    p_value, q_value = exact_pq_at_rational((2, 2, 2), 1, 5)
    scaled_p = sum(value * 5 ** (len(p_poly) - 1 - degree) for degree, value in enumerate(p_poly))
    scaled_q = sum(value * 5 ** (len(q_poly) - 1 - degree) for degree, value in enumerate(q_poly))
    _check(
        checks,
        "route1_spin_enumeration_control",
        p_value == scaled_p and q_value == scaled_q,
        "complete 2x2x2 polynomial enumeration agrees with exact parity transfer at t=1/5",
    )
    p, q = 2, 9
    full_p, full_q = exact_pq_at_rational((2, 3, 2), p, q)
    residual = exact_criterion_residual((2, 3, 2), p, q)
    _check(
        checks,
        "route1_residual_stream_control",
        residual == p * full_q - q * full_p,
        "direct residual stream equals separately reconstructed pQ-qP on 2x3x2",
    )


def _route1_search(checks: list[dict[str, object]]) -> dict[str, object]:
    """Certify the predetermined full-transfer family with integer signs."""

    _route1_controls(checks)
    records: list[dict[str, object]] = []
    challenge_p, challenge_q = CHALLENGE
    for shape in ROUTE1_SHAPES:
        started = time.monotonic()
        locator_v, locator_k = solve_box_root_float(shape, iterations=27)
        denominator = ROUTE1_SMALL_DENOMINATOR
        exact_started = time.monotonic()
        safe_p, safe_q, safe_residual = _reduced_safe_point(shape, locator_v, denominator)
        challenge_residual = exact_criterion_residual(shape, challenge_p, challenge_q)
        exact_elapsed = time.monotonic() - exact_started
        safe_interval, safe_lower = _atanh_record(safe_p, safe_q)
        record = {
            "shape": list(shape),
            "sites": math.prod(shape),
            "internal_bonds": _box_edge_count(shape),
            "cross_section_sites": shape[0] * shape[1],
            "full_transfer_states": 1 << (shape[0] * shape[1]),
            "float_locator_only": {
                "v_root_locator": format(locator_v, ".12g"),
                "K_root_locator": format(locator_k, ".12g"),
                "status": "[COMPUTATION] binary64 locator only; exact signs decide every stored certificate",
            },
            "safe_exact_t": f"{safe_p}/{safe_q}",
            "safe_residual_tQ_minus_P": _fingerprint(safe_residual),
            "safe_K_interval": safe_interval,
            "safe_K_decimal_lower": safe_lower,
            "safe_point_beats_incumbent_exactly": safe_p**36 * C36 > safe_q**36,
            "challenge_exact_t": f"{challenge_p}/{challenge_q}",
            "challenge_residual_tQ_minus_P": _fingerprint(challenge_residual),
            "challenge_passes": challenge_residual < 0,
            "exact_residual_elapsed_seconds": exact_elapsed,
            "total_shape_elapsed_seconds": time.monotonic() - started,
        }
        _check(
            checks,
            "route1_safe_" + "x".join(map(str, shape)),
            safe_residual < 0,
            f"exact integer tQ-P<0 at t={safe_p}/{safe_q}",
        )
        records.append(record)

    best = max(records, key=lambda item: Decimal(str(item["safe_K_decimal_lower"])))
    global_improved = any(bool(record["safe_point_beats_incumbent_exactly"]) for record in records)
    _check(
        checks,
        "route1_fixed_family_has_wider_cross_sections",
        all(record["cross_section_sites"] > 16 for record in records),
        "every completed prism has cross-section strictly larger than the prior <=16 audit",
    )
    return {
        "status": (
            "[THEOREM from EXACT COMPUTATION] a wider Simon--Lieb certificate improves the lower endpoint"
            if global_improved
            else "[UNRESOLVED] the completed wider exact Simon--Lieb family remains below the c_36 endpoint"
        ),
        "selection": (
            "[COMPUTATION] fixed shapes (3x6)x{4,12} and (4x5)x{4,16}; "
            "no benchmark selects a shape or rational point"
        ),
        "state_scaling": "the audited full parity transfer uses exactly 2^(a*b) states for an a-by-b cross-section",
        "challenge": {
            "t": "23/110",
            "reason": "exactly above c_36^(-1/36), as checked by (23^36)*c_36 > 110^36",
            "all_completed_families_reject": not any(record["challenge_passes"] for record in records),
        },
        "certificates": records,
        "best_exact_safe_certificate": best,
        "global_interval_improved": global_improved,
    }


# ---------------------------------------------------------------------------
# Route 2: exact finite-volume infrared/GKS LP floor
# ---------------------------------------------------------------------------


def _signed_permutation_orbit(point: tuple[int, int, int], side: int) -> tuple[tuple[int, int, int], ...]:
    images = set()
    for permutation in itertools.permutations(range(3)):
        reordered = tuple(point[index] for index in permutation)
        for signs in itertools.product((-1, 1), repeat=3):
            images.add(tuple((sign * coordinate) % side for sign, coordinate in zip(signs, reordered, strict=True)))
    return tuple(sorted(images))


def _torus_orbits(side: int) -> tuple[tuple[tuple[int, int, int], ...], ...]:
    remaining = set(itertools.product(range(side), repeat=3))
    orbits: list[tuple[tuple[int, int, int], ...]] = []
    while remaining:
        point = min(remaining)
        orbit = _signed_permutation_orbit(point, side)
        remaining.difference_update(orbit)
        orbits.append(orbit)
    return tuple(orbits)


def _torus_green(side: int) -> dict[tuple[int, int, int], Fraction]:
    """Exact zero-mode-removed Green values for the rational-cosine tori."""

    cosines = RATIONAL_COSINES[side]
    sites = list(itertools.product(range(side), repeat=3))
    modes = [mode for mode in sites if mode != (0, 0, 0)]
    volume = side**3
    result: dict[tuple[int, int, int], Fraction] = {}
    for site in sites:
        total = Fraction(0)
        for mode in modes:
            lam = 3 - sum(cosines[coordinate] for coordinate in mode)
            phase = sum(mode[axis] * site[axis] for axis in range(3)) % side
            total += cosines[phase] / lam
        result[site] = total / volume
    return result


def _lp_problem(side: int, coupling: Fraction) -> dict[str, object]:
    """Build max S subject to Parseval, IR, GKS-I, |G|<=1, and energy input.

    Variables x_j are common values of Ghat(k) on signed-permutation momentum
    orbits.  Symmetrising any feasible physical Ghat over the cubic group
    preserves every displayed constraint and S, so this finite LP is a valid
    relaxation for the actual torus correlation function.
    """

    if side not in RATIONAL_COSINES or coupling <= 0:
        raise ValueError("side must have rational cosine data and K must be positive")
    cosines = RATIONAL_COSINES[side]
    volume = side**3
    all_orbits = _torus_orbits(side)
    mode_orbits = tuple(orbit for orbit in all_orbits if orbit != ((0, 0, 0),))
    site_orbits = tuple(orbit for orbit in all_orbits if orbit != ((0, 0, 0),))
    multiplicities = [len(orbit) for orbit in mode_orbits]
    dispersions = [
        3 - sum(cosines[coordinate] for coordinate in orbit[0]) for orbit in mode_orbits
    ]
    if any(value <= 0 for value in dispersions):
        raise AssertionError("nonzero torus mode must have positive dispersion")

    phase_sums: list[list[Fraction]] = []
    for site_orbit in site_orbits:
        site = site_orbit[0]
        row: list[Fraction] = []
        for mode_orbit in mode_orbits:
            row.append(
                sum(
                    cosines[sum(mode[axis] * site[axis] for axis in range(3)) % side]
                    for mode in mode_orbit
                )
            )
        phase_sums.append(row)

    rows: list[list[Fraction]] = []
    rhs: list[Fraction] = []
    labels: list[str] = []
    # Infrared ceilings x_j <= 1/(2K lambda_j).
    for index, dispersion in enumerate(dispersions):
        row = [Fraction(0) for _ in mode_orbits]
        row[index] = Fraction(1)
        rows.append(row)
        rhs.append(Fraction(1, 2) / (coupling * dispersion))
        labels.append(f"ceiling_mode_{index}")
    # G(z)>=0 and G(z)<=1 for one representative of every nonzero site orbit.
    for orbit, phase in zip(site_orbits, phase_sums, strict=True):
        site = orbit[0]
        gap = [Fraction(multiple) - coefficient for multiple, coefficient in zip(multiplicities, phase, strict=True)]
        rows.append(gap)
        rhs.append(Fraction(volume))
        labels.append("gks_nonnegative_" + "_".join(map(str, site)))
        rows.append([-value for value in gap])
        rhs.append(Fraction(0))
        labels.append("spin_upper_" + "_".join(map(str, site)))
    # M_L^2=1-S is nonnegative: sum multiplicity*x <= |T_L|.
    rows.append([Fraction(value) for value in multiplicities])
    rhs.append(Fraction(volume))
    labels.append("zero_mode_nonnegative")
    # Exact energy/infrared consequence G(e)>=1-(1-1/N)/(6K).
    neighbour_index = next(index for index, orbit in enumerate(site_orbits) if (1, 0, 0) in orbit)
    neighbour_phase = phase_sums[neighbour_index]
    gap_e = [Fraction(multiple) - coefficient for multiple, coefficient in zip(multiplicities, neighbour_phase, strict=True)]
    energy_rhs = Fraction(volume - 1, 1) / (6 * coupling)
    rows.append(gap_e)
    rhs.append(energy_rhs)
    labels.append("energy_nearest_neighbour")

    return {
        "side": side,
        "coupling": coupling,
        "volume": volume,
        "mode_orbits": mode_orbits,
        "site_orbits": site_orbits,
        "multiplicities": multiplicities,
        "dispersions": dispersions,
        "phase_sums": phase_sums,
        "a": rows,
        "b": rhs,
        "c": [Fraction(value) for value in multiplicities],
        "labels": labels,
        "green": _torus_green(side),
    }


def _solve_square(matrix: Sequence[Sequence[Fraction]], rhs: Sequence[Fraction]) -> list[Fraction]:
    """Exact Gauss--Jordan solve for a nonsingular square Fraction matrix."""

    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix) or len(rhs) != size:
        raise ValueError("expected a nonempty square system")
    augmented = [list(row) + [rhs[index]] for index, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((row for row in range(column, size) if augmented[row][column] != 0), None)
        if pivot is None:
            raise ValueError("singular matrix")
        if pivot != column:
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value / divisor for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [
                    value - factor * pivot_value
                    for value, pivot_value in zip(augmented[row], augmented[column], strict=True)
                ]
    return [row[-1] for row in augmented]


def _simplex_max(
    a: Sequence[Sequence[Fraction]], b: Sequence[Fraction], c: Sequence[Fraction]
) -> dict[str, object]:
    """Bland-rule exact primal simplex for max c*x, A*x<=b, x>=0.

    The returned dual is independently reconstructed from B^T y=c_B.  Exact
    primal feasibility, dual feasibility, and equal objectives are verified
    before the certificate is returned.
    """

    row_count = len(a)
    variable_count = len(c)
    if not row_count or any(len(row) != variable_count for row in a):
        raise ValueError("malformed constraint matrix")
    if len(b) != row_count or any(value < 0 for value in b):
        raise ValueError("initial slack basis requires nonnegative RHS")
    full_a = [
        list(row) + [Fraction(int(row_index == column)) for column in range(row_count)]
        for row_index, row in enumerate(a)
    ]
    tableau = [list(row) + [b[index]] for index, row in enumerate(full_a)]
    full_c = list(c) + [Fraction(0) for _ in range(row_count)]
    basis = [variable_count + row for row in range(row_count)]
    limit = 10_000
    for _ in range(limit):
        c_basis = [full_c[index] for index in basis]
        reduced = [
            full_c[column]
            - sum(c_basis[row] * tableau[row][column] for row in range(row_count))
            for column in range(variable_count + row_count)
        ]
        entering = next((column for column, value in enumerate(reduced) if value > 0), None)
        if entering is None:
            break
        ratios = [
            (tableau[row][-1] / tableau[row][entering], basis[row], row)
            for row in range(row_count)
            if tableau[row][entering] > 0
        ]
        if not ratios:
            raise ValueError("unbounded LP")
        _, _, leaving = min(ratios)
        divisor = tableau[leaving][entering]
        tableau[leaving] = [value / divisor for value in tableau[leaving]]
        for row in range(row_count):
            if row == leaving:
                continue
            factor = tableau[row][entering]
            if factor:
                tableau[row] = [
                    value - factor * pivot_value
                    for value, pivot_value in zip(tableau[row], tableau[leaving], strict=True)
                ]
        basis[leaving] = entering
    else:
        raise RuntimeError("simplex iteration cap reached")

    primal_full = [Fraction(0) for _ in range(variable_count + row_count)]
    for row, column in enumerate(basis):
        primal_full[column] = tableau[row][-1]
    primal = primal_full[:variable_count]
    objective = sum(value * coefficient for value, coefficient in zip(primal, c, strict=True))

    basis_matrix = [
        [full_a[row][basis[column]] for column in range(row_count)] for row in range(row_count)
    ]
    c_basis = [full_c[column] for column in basis]
    transpose = [
        [basis_matrix[column][row] for column in range(row_count)] for row in range(row_count)
    ]
    dual = _solve_square(transpose, c_basis)
    dual_objective = sum(value * bound for value, bound in zip(dual, b, strict=True))
    primal_feasible = all(value >= 0 for value in primal) and all(
        sum(coefficient * value for coefficient, value in zip(row, primal, strict=True)) <= bound
        for row, bound in zip(a, b, strict=True)
    )
    dual_feasible = all(value >= 0 for value in dual) and all(
        sum(dual[row] * full_a[row][column] for row in range(row_count)) >= full_c[column]
        for column in range(variable_count + row_count)
    )
    if not (primal_feasible and dual_feasible and objective == dual_objective):
        raise AssertionError("exact primal/dual LP certificate failed")
    return {
        "primal": primal,
        "dual": dual,
        "objective": objective,
        "basis": basis,
        "primal_feasible": primal_feasible,
        "dual_feasible": dual_feasible,
    }


def _serialize_lp_certificate(problem: dict[str, object], solved: dict[str, object]) -> dict[str, object]:
    volume = int(problem["volume"])
    objective = solved["objective"]
    assert isinstance(objective, Fraction)
    primal = solved["primal"]
    dual = solved["dual"]
    assert isinstance(primal, list) and isinstance(dual, list)
    return {
        "coupling": _fraction_text(problem["coupling"]),
        "objective_sum_nonzero_modes": _fraction_text(objective),
        "S_star": _fraction_text(objective / volume),
        "magnetization_square_floor": _fraction_text(Fraction(1) - objective / volume),
        "mode_orbits": [
            {"representative": list(orbit[0]), "multiplicity": multiplicity, "lambda": _fraction_text(dispersion)}
            for orbit, multiplicity, dispersion in zip(
                problem["mode_orbits"], problem["multiplicities"], problem["dispersions"], strict=True
            )
        ],
        "constraint_labels": list(problem["labels"]),
        "primal_x": [_fraction_text(value) for value in primal],
        "dual_y": [_fraction_text(value) for value in dual],
        "basis": list(solved["basis"]),
        "primal_dual_equal": True,
    }


def _route2_floors(checks: list[dict[str, object]]) -> dict[str, object]:
    upper_data = json.loads(UPPER_INFRARED_PATH.read_text(encoding="utf-8"))["data"]
    stored_green = {
        int(row["side"]): row for row in upper_data["finite_volume_audit"]["exact_torus_green_functions"]
    }
    side_records: list[dict[str, object]] = []
    for side in LP_SIDES:
        green = _torus_green(side)
        origin = green[(0, 0, 0)]
        minimum = min(green.values())
        minimum_sites = sorted(site for site, value in green.items() if value == minimum)
        stored = stored_green[side]
        _check(
            checks,
            f"route2_stored_green_L{side}",
            _fraction_text(origin) == str(stored["C_L_at_origin"])
            and _fraction_text(minimum) == str(stored["minimum_value"])
            and sum(green.values(), Fraction(0)) == 0,
            "direct exact mode sum reproduces stored C_L(0), min C_L, and zero spatial sum",
        )
        ceiling_feasible_threshold = (origin - minimum) / 2
        certificates: list[dict[str, object]] = []
        for coupling in LP_POINTS:
            problem = _lp_problem(side, coupling)
            solved = _simplex_max(problem["a"], problem["b"], problem["c"])
            certificate = _serialize_lp_certificate(problem, solved)
            certificates.append(certificate)
            objective = solved["objective"]
            assert isinstance(objective, Fraction)
            _check(
                checks,
                f"route2_primal_dual_L{side}_K{coupling.numerator}_{coupling.denominator}",
                bool(solved["primal_feasible"])
                and bool(solved["dual_feasible"])
                and objective == sum(
                    _parse_fraction(value) * _parse_fraction(bound)
                    for value, bound in zip(certificate["dual_y"], [_fraction_text(value) for value in problem["b"]], strict=True)
                ),
                "exact Fraction primal and dual are feasible with equal objectives",
            )
        # The all-ceiling profile is an exact witness that the finite LP reduces
        # to the naive C_L(0)/(2K) bound above this threshold.  Check it at the
        # largest fixed rational grid point without relying on floating point.
        coupling = LP_POINTS[-1]
        problem = _lp_problem(side, coupling)
        all_ceiling = [Fraction(1, 2) / (coupling * dispersion) for dispersion in problem["dispersions"]]
        all_ceiling_feasible = all(
            sum(coefficient * value for coefficient, value in zip(row, all_ceiling, strict=True)) <= bound
            for row, bound in zip(problem["a"], problem["b"], strict=True)
        )
        all_ceiling_objective = sum(
            Fraction(multiplicity) * value
            for multiplicity, value in zip(problem["multiplicities"], all_ceiling, strict=True)
        )
        _check(
            checks,
            f"route2_ceiling_witness_L{side}",
            coupling >= ceiling_feasible_threshold
            and all_ceiling_feasible
            and all_ceiling_objective / side**3 == origin / (2 * coupling),
            f"explicit all-ceiling profile is feasible and saturates the naive finite-volume sum at K={_fraction_text(coupling)}",
        )
        side_records.append(
            {
                "side": side,
                "volume": side**3,
                "C_L_at_origin": _fraction_text(origin),
                "minimum_C_L": _fraction_text(minimum),
                "minimum_sites": [list(site) for site in minimum_sites],
                "ceiling_profile_feasible_if_K_at_least": _fraction_text(ceiling_feasible_threshold),
                "finite_floor_certificates": certificates,
                "uniform_missing_lemma": (
                    "[UNRESOLVED] There must exist K<I3/2 and m>0 such that "
                    "M_L(K)^2>=m for every sufficiently large even L.  A finite-L "
                    "LP certificate proves only M_L(K)^2>=1-S_star(K,L) separately "
                    "for that L; it cannot replace this quantified uniform statement."
                ),
            }
        )
    _check(
        checks,
        "route2_threshold_increases_on_stored_even_tori",
        Fraction(side_records[0]["ceiling_profile_feasible_if_K_at_least"])
        < Fraction(side_records[1]["ceiling_profile_feasible_if_K_at_least"]),
        "exact C_L(0)-min C_L threshold is larger at L=6 than L=4",
    )
    return {
        "status": "[THEOREM from EXACT COMPUTATION] finite even-torus Parseval/infrared/GKS/energy constraint system gives the stated exact floors",
        "inequality_chain": (
            "M_L^2=1-|T_L|^-1 sum_{k!=0} Ghat(k) >= 1-S_star(K,L); "
            "S_star is the maximum of the exact finite LP with Ghat>=0, "
            "Ghat(k)<=1/(2K lambda(k)), 0<=G(z)<=1, M_L^2>=0, and "
            "G(e)>=1-(1-|T_L|^-1)/(6K)."
        ),
        "valid_even_tori": [4, 6],
        "excluded_stored_sides": {
            "2": "[FALSIFIED control in upper_infrared] the infrared ceiling is false under the degenerate L=2 test convention",
            "3": "[UNRESOLVED for this route] the audited Gaussian-domination statement is invoked only on even tori",
        },
        "records": side_records,
        "uniformity_conclusion": (
            "[UNRESOLVED] the finite exact floors do not move K_c without the displayed "
            "uniform-in-L lemma.  The exact all-ceiling witnesses above their finite-L "
            "thresholds identify why the finite constraint system reverts to the naive "
            "C_L(0)/(2K) barrier."
        ),
    }


# ---------------------------------------------------------------------------
# Route 3: auditable FKG-free finite SAW inequalities
# ---------------------------------------------------------------------------


def _saw_counts(max_steps: int, forbidden: tuple[int, int, int] | None = None) -> list[int]:
    """Exact rooted cubic SAW counts, optionally avoiding one fixed neighbour."""

    counts = [0 for _ in range(max_steps + 1)]
    origin = (0, 0, 0)
    steps = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    visited = {origin}

    def visit(point: tuple[int, int, int], depth: int) -> None:
        counts[depth] += 1
        if depth == max_steps:
            return
        for delta in steps:
            next_point = (point[0] + delta[0], point[1] + delta[1], point[2] + delta[2])
            if next_point not in visited and next_point != forbidden:
                visited.add(next_point)
                visit(next_point, depth + 1)
                visited.remove(next_point)

    visit(origin, 0)
    return counts


def _route3_saw(checks: list[dict[str, object]]) -> dict[str, object]:
    started = time.monotonic()
    small_counts = _saw_counts(9)
    avoid_counts = _saw_counts(8, forbidden=(1, 0, 0))
    expected_small = (1, 6, 30, 150, 726, 3534, 16926, 81390, 387966, 1853886)
    _check(
        checks,
        "route3_small_saw_enumeration",
        tuple(small_counts) == expected_small,
        "direct visited-set backtracking reproduces c_0 through c_9 exactly",
    )
    source_hash = hashlib.sha256(SAW_SOURCE.read_bytes()).hexdigest()
    _check(
        checks,
        "route3_external_saw_source_hash",
        source_hash == C36_SHA256,
        "cached Schram--Barkema--Bisseling PDF has the audited SHA-256",
    )
    counts = {index: value for index, value in enumerate(small_counts)} | EXTERNAL_SAW_COUNTS

    # Injection theorem: take a pair (omega,eta) of m- and n-step walks.  If
    # eta first backtracks from omega_m to omega_{m-1}, then its remaining
    # n-1-step tail must avoid omega_m in order that eta itself stay SAW.  There
    # are a_{n-1} such tails, and every resulting concatenation is invalid.
    # Thus at least c_m a_{n-1} of c_m c_n pairs fail to concatenate.
    small_injection_checks: list[dict[str, object]] = []
    for total in range(2, 10):
        for n in range(1, total):
            m = total - n
            lhs = counts[total]
            rhs = counts[m] * (counts[n] - avoid_counts[n - 1])
            passed = lhs <= rhs
            small_injection_checks.append(
                {"m": m, "n": n, "lhs_c_m_plus_n": lhs, "rhs": rhs, "passed": passed}
            )
            _check(
                checks,
                f"route3_injection_m{m}_n{n}",
                passed,
                "c_(m+n)<=c_m(c_n-a_(n-1)) by the explicit first-backtrack invalid-pair injection",
            )
    mixed_pairs = ((30, 6), (31, 5), (32, 4), (33, 3), (34, 2), (35, 1))
    mixed_injection_checks: list[dict[str, object]] = []
    for m, n in mixed_pairs:
        lhs = counts[m + n]
        rhs = counts[m] * (counts[n] - avoid_counts[n - 1])
        passed = lhs <= rhs
        mixed_injection_checks.append({"m": m, "n": n, "lhs_c_m_plus_n": lhs, "rhs": rhs, "passed": passed})
        _check(
            checks,
            f"route3_injection_m{m}_n{n}",
            passed,
            "external exact c_(m+n) obeys the FKG-free first-backtrack inequality",
        )

    # Fekete at N=36 gives mu^2 <= c_36^(1/18).  It proves individual finite
    # Kesten-ratio instances exactly when c_36*c_n^18 <= c_(n+2)^18.
    proven_ratio_rows: list[dict[str, object]] = []
    for n in range(8):
        lhs = C36 * counts[n] ** 18
        rhs = counts[n + 2] ** 18
        passes = lhs <= rhs
        row = {
            "n": n,
            "ratio": f"{counts[n + 2]}/{counts[n]}",
            "fekete_cross_product": _fingerprint(lhs - rhs),
            "proved_by_c36_fekete": passes,
        }
        proven_ratio_rows.append(row)
        if n <= 6:
            _check(
                checks,
                f"route3_kesten_instance_n{n}",
                passes,
                "exact c_36*c_n^18<=c_(n+2)^18 proves mu^2<=c_(n+2)/c_n",
            )
        else:
            _check(
                checks,
                "route3_first_fekete_failure_n7",
                not passes,
                "the same exact Fekete certificate cannot prove the n=7 ratio instance",
            )

    # The direct injection bound on mu is checked over every available m in
    # {30,...,36} and n in {1,...,8}.  Each comparison with c_36^(1/36) is
    # exact after raising both positive bases to the common power 36.
    injection_mu_rows: list[dict[str, object]] = []
    all_not_better = True
    for m in range(30, 37):
        for n in range(1, 9):
            base = counts[m] * (counts[n] - avoid_counts[n - 1])
            exponent = m + n
            not_better = base**36 >= C36**exponent
            all_not_better &= not_better
            injection_mu_rows.append(
                {
                    "m": m,
                    "n": n,
                    "base": base,
                    "exponent": exponent,
                    "base_power_36_minus_c36_power_exponent": _fingerprint(base**36 - C36**exponent),
                    "not_better_than_c36_root_exactly": not_better,
                }
            )
    _check(
        checks,
        "route3_injection_refinement_does_not_beat_c36",
        all_not_better,
        "every available exact first-backtrack mu certificate is at least c_36^(1/36)",
    )

    # n=30 is the earliest fully available high-n conditional ratio whose
    # value would improve the incumbent.  This is an exact comparison, not a
    # decimal estimate: (c_32/c_30)^18 < c_36.
    conditional_improvement = counts[32] ** 18 < counts[30] ** 18 * C36
    _check(
        checks,
        "route3_n30_ratio_would_improve_incumbent",
        conditional_improvement,
        "exactly c_32/c_30 < c_36^(1/18), so an n=30 finite-ratio proof would improve K_c",
    )
    return {
        "status": "[UNRESOLVED] no first-principles proof of the improving n=30 Kesten-ratio inequality was found",
        "source": {
            "path": "sources/fulltext/schram_barkema_bisseling2011.pdf",
            "sha256": source_hash,
            "external_counts": {str(key): value for key, value in EXTERNAL_SAW_COUNTS.items()},
        },
        "small_exact_counts": small_counts,
        "avoid_fixed_neighbour_counts_a_n": avoid_counts,
        "injection_theorem": {
            "status": "[THEOREM]",
            "statement": "For m,n>=1, c_(m+n)<=c_m*(c_n-a_(n-1)), where a_r counts r-step SAWs from 0 avoiding a specified nearest neighbour.",
            "proof": "Inject c_m*a_(n-1) pairs into non-concatenable (m,n) pairs by forcing eta's first edge to backtrack omega's final edge; its tail avoids eta's start so eta is SAW, while the concatenation revisits omega_(m-1).",
            "small_checks": small_injection_checks,
            "mixed_external_checks": mixed_injection_checks,
            "mu_refinement_comparisons": injection_mu_rows,
        },
        "fekete_certified_finite_ratio_instances": {
            "status": "[THEOREM from EXTERNAL EXACT COMPUTATION]",
            "logic": "mu<=c_36^(1/36); c_36*c_n^18<=c_(n+2)^18 implies mu^2<=c_(n+2)/c_n.",
            "rows_n_0_through_7": proven_ratio_rows,
            "proved_n": list(range(7)),
            "first_not_proved_by_this_certificate": 7,
        },
        "missing_improving_lemma": {
            "status": "[UNRESOLVED]",
            "statement": "Prove mu^2<=c_32/c_30 by a first-principles finite-n argument (or a theorem valid at n=30).",
            "exact_value_comparison": "c_32^18 < c_30^18*c_36, equivalently c_32/c_30 < c_36^(1/18)",
            "why_material": "It would yield v_c>=sqrt(c_30/c_32)>c_36^(-1/36), strictly improving the current lower endpoint.",
        },
        "elapsed_seconds": time.monotonic() - started,
    }


# ---------------------------------------------------------------------------
# Artifact assembly
# ---------------------------------------------------------------------------


def main() -> None:
    checks: list[dict[str, object]] = []
    started = time.monotonic()
    _check(
        checks,
        "challenge_strictly_above_incumbent_exactly",
        23**36 * C36 > 110**36,
        "(23/110)^36*c_36>1, so t=23/110 exceeds c_36^(-1/36) without decimals",
    )
    route1 = _route1_search(checks)
    route2 = _route2_floors(checks)
    route3 = _route3_saw(checks)
    best_route1 = route1["best_exact_safe_certificate"]
    assert isinstance(best_route1, dict)
    route1_improved = bool(route1["global_interval_improved"])
    final_lower = (
        str(best_route1["safe_K_decimal_lower"]) if route1_improved else INCUMBENT_K_LOWER
    )
    final_status = (
        "[THEOREM from EXACT COMPUTATION] widened Simon--Lieb certificate improves the lower endpoint"
        if route1_improved
        else "[UNRESOLVED] exact route-1 family stays below c_36; route-2 lacks a uniform-in-L floor; route-3 lacks the n=30 finite-ratio lemma."
    )
    classification = (
        "[THEOREM from EXACT COMPUTATION] route 1 supplies a strictly improved global lower endpoint"
        if route1_improved
        else "[UNRESOLVED] No route supplies a strictly improved global K_c endpoint; all three yield new exact finite certificates and precise missing lemmas."
    )
    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "arithmetic": (
                "Route 1 decisions are Python integer residual signs; route 2 LP arithmetic is Fraction; "
                "route 3 comparisons are Python integers. Float64 appears only in explicitly non-authoritative root locators; "
                "mpmath.iv with directed rounding formats atanh intervals."
            ),
            "mpmath_iv_dps": DPS,
            "elapsed_seconds": time.monotonic() - started,
        },
        "data": {
            "classification": classification,
            "model": "nearest-neighbour ferromagnetic Ising model on Z^3; K=beta*J and t=tanh(K)",
            "incumbent_interval": {
                "status": "[THEOREM from EXTERNAL EXACT COMPUTATION] unchanged",
                "exact_lower": "atanh(c_36^(-1/36)), c_36=2941370856334701726560670",
                "decimal_lower": INCUMBENT_K_LOWER,
                "decimal_upper": "0.2527310098586630030260020266135701299926",
            },
            "route1_wider_simon_lieb": route1,
            "route2_finite_volume_floor": route2,
            "route3_kesten_attempt": route3,
            "final_certified_interval": {
                "lower": final_lower,
                "upper": "0.2527310098586630030260020266135701299926",
                "improved": route1_improved,
                "reason": final_status,
            },
        },
        "checks": checks,
    }
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    RESULT_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"WROTE {RESULT_PATH.relative_to(ROOT)}")
    print(
        "PASS e85_kc_interval2: exact widened Simon--Lieb, finite-torus floors, and FKG-free SAW certificates; "
        f"final interval [{final_lower}, 0.2527310098586630030260020266135701299926]"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL e85_kc_interval2: {error}")
        raise
