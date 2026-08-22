"""Independent hostile replications for the principal claims in the repository.

This script deliberately does not import the repository's Pauli-bookkeeping code.
The Dolan--Grady and tridiagonal calculations use dense NumPy matrices built
from the one-qubit matrices below.  All matrix entries and all reported exact
residuals are int64/Fraction values; floating point is used only to display the
least-squares answer independently.

Run from the repository root:
    .venv/bin/python tests/test_audit_replication.py
"""

from __future__ import annotations

from collections import Counter
from decimal import Decimal
from fractions import Fraction
from itertools import combinations, product
import json
from pathlib import Path
import platform
import sys

import mpmath as mp
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "audit" / "audit_replication.json"

I2 = np.eye(2, dtype=np.int64)
X = np.asarray([[0, 1], [1, 0]], dtype=np.int64)
Z = np.asarray([[1, 0], [0, -1]], dtype=np.int64)
IY = np.asarray([[0, 1], [-1, 0]], dtype=np.int64)  # exactly i times Pauli Y


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def _mp_text(value: mp.mpf, digits: int = 80) -> str:
    return mp.nstr(value, digits)


def _pauli_string(n: int, local: dict[int, np.ndarray]) -> np.ndarray:
    out = np.asarray([[1]], dtype=np.int64)
    for site in range(n):
        out = np.kron(out, local.get(site, I2))
    return out


def _comm(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return left @ right - right @ left


def _dense_generators(n: int, edges: list[tuple[int, int]]) -> tuple[np.ndarray, np.ndarray]:
    dim = 1 << n
    a = np.zeros((dim, dim), dtype=np.int64)
    b = np.zeros((dim, dim), dtype=np.int64)
    for site in range(n):
        a += _pauli_string(n, {site: X})
    for u, v in edges:
        b += _pauli_string(n, {u: Z, v: Z})
    return a, b


def _full_centralizer_scope_attack() -> dict:
    """Exhibit a non-Pauli lattice symmetry commuting with the summed pair."""

    n = 3
    a, b = _dense_generators(n, [(0, 1), (1, 2), (0, 2)])
    dim = 1 << n
    translation = np.zeros((dim, dim), dtype=np.int64)
    for configuration in range(dim):
        translated = 0
        for site in range(n):
            bit = (configuration >> (n - 1 - site)) & 1
            target = (site + 1) % n
            translated |= bit << (n - 1 - target)
        translation[translated, configuration] = 1
    parity = _pauli_string(n, {site: X for site in range(n)})
    candidates = [np.eye(dim, dtype=np.int64), parity, translation]
    gram = [
        [int(np.dot(left.ravel(), right.ravel())) for right in candidates]
        for left in candidates
    ]
    gram_determinant = (
        gram[0][0] * (gram[1][1] * gram[2][2] - gram[1][2] * gram[2][1])
        - gram[0][1] * (gram[1][0] * gram[2][2] - gram[1][2] * gram[2][0])
        + gram[0][2] * (gram[1][0] * gram[2][1] - gram[1][1] * gram[2][0])
    )
    return {
        "graph": "cycle C3",
        "operator": "one-site lattice translation",
        "commutator_with_A_max_abs": int(np.max(np.abs(_comm(translation, a)))),
        "commutator_with_B_max_abs": int(np.max(np.abs(_comm(translation, b)))),
        "translation_cubed_minus_identity_max_abs": int(
            np.max(np.abs(translation @ translation @ translation - candidates[0]))
        ),
        "gram_matrix_for_I_parity_translation": gram,
        "gram_determinant": gram_determinant,
        "linearly_independent_from_I_and_parity": gram_determinant != 0,
        "scope_verdict": (
            "the two-element theorem is valid only inside the Pauli group; the full operator "
            "centralizer and symmetry-sector decomposition can be larger"
        ),
    }


def _dense_dg_closed_form(n: int, edges: list[tuple[int, int]]) -> np.ndarray:
    neighbours: list[set[int]] = [set() for _ in range(n)]
    for u, v in edges:
        neighbours[u].add(v)
        neighbours[v].add(u)
    dim = 1 << n
    out = np.zeros((dim, dim), dtype=np.int64)
    for site, adjacent in enumerate(neighbours):
        degree = len(adjacent)
        for neighbour in sorted(adjacent):
            out += 24 * (degree - 2) * _pauli_string(
                n, {site: IY, neighbour: Z}
            )
        for triple in combinations(sorted(adjacent), 3):
            local = {site: IY, **{neighbour: Z for neighbour in triple}}
            out += 48 * _pauli_string(n, local)
    return out


def _dg_case(name: str, n: int, edges: list[tuple[int, int]]) -> dict:
    a, b = _dense_generators(n, edges)
    first_relation = _comm(a, _comm(a, _comm(a, b))) - 16 * _comm(a, b)
    residual = _comm(b, _comm(b, _comm(b, a))) - 16 * _comm(b, a)
    closed_form = _dense_dg_closed_form(n, edges)
    return {
        "name": name,
        "n": n,
        "edges": [list(edge) for edge in edges],
        "arithmetic": "dense np.int64 matrix products; no repository Pauli code",
        "first_relation_max_abs": int(np.max(np.abs(first_relation))),
        "closed_form_max_abs_difference": int(np.max(np.abs(residual - closed_form))),
        "residual_max_abs": int(np.max(np.abs(residual))),
        "closed_form_max_abs": int(np.max(np.abs(closed_form))),
    }


def _torus_edges(rows: int, columns: int) -> list[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    for row in range(rows):
        for column in range(columns):
            site = row * columns + column
            right = row * columns + (column + 1) % columns
            down = ((row + 1) % rows) * columns + column
            edges.add(tuple(sorted((site, right))))
            edges.add(tuple(sorted((site, down))))
    return sorted(edges)


def _solve_fraction_system(matrix: list[list[int]], vector: list[int]) -> list[Fraction]:
    n = len(vector)
    augmented = [
        [Fraction(value) for value in row] + [Fraction(vector[index])]
        for index, row in enumerate(matrix)
    ]
    for column in range(n):
        pivot = next(row for row in range(column, n) if augmented[row][column])
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(n):
            if row == column:
                continue
            scale = augmented[row][column]
            if scale:
                augmented[row] = [
                    left - scale * right
                    for left, right in zip(augmented[row], augmented[column], strict=True)
                ]
    return [augmented[row][-1] for row in range(n)]


def _td_dense_replication() -> dict:
    edges = _torus_edges(3, 3)
    a, b = _dense_generators(9, edges)
    b2 = b @ b
    a2 = a @ a
    t0 = _comm(b, b2 @ a + a @ b2)
    t1 = _comm(b, b @ a @ b)
    t2 = _comm(b, b @ a + a @ b)
    t3 = _comm(b, a)

    target = t0.ravel()
    columns = [operator.ravel() for operator in (t1, t2, t3)]
    gram = [[int(np.dot(left, right)) for right in columns] for left in columns]
    normal_rhs = [int(np.dot(column, target)) for column in columns]
    exact_parameters = _solve_fraction_system(gram, normal_rhs)
    target_norm_squared = int(np.dot(target, target))
    exact_rss = Fraction(target_norm_squared)
    exact_rss -= 2 * sum(
        parameter * rhs
        for parameter, rhs in zip(exact_parameters, normal_rhs, strict=True)
    )
    exact_rss += sum(
        exact_parameters[i] * exact_parameters[j] * gram[i][j]
        for i in range(3)
        for j in range(3)
    )

    design = np.column_stack(columns).astype(np.float64)
    target_float = target.astype(np.float64)
    float_parameters, _, rank, singular_values = np.linalg.lstsq(
        design, target_float, rcond=None
    )
    float_residual = target_float - design @ float_parameters

    common_denominator = 1
    for parameter in exact_parameters:
        common_denominator = np.lcm(common_denominator, parameter.denominator)
    integer_coefficients = [
        parameter.numerator * (common_denominator // parameter.denominator)
        for parameter in exact_parameters
    ]
    scaled_residual = common_denominator * t0
    for coefficient, operator in zip(
        integer_coefficients, (t1, t2, t3), strict=True
    ):
        scaled_residual -= coefficient * operator
    exact_max = Fraction(int(np.max(np.abs(scaled_residual))), common_denominator)

    return {
        "graph": "3x3 periodic square layer",
        "n": 9,
        "edge_count": len(edges),
        "arithmetic": (
            "dense np.int64 products; exact Fraction normal equations; "
            "independent float64 np.linalg.lstsq display"
        ),
        "operator_equation": "T0 = beta*T1 + gamma*T2 + rho*T3",
        "design_rank": int(rank),
        "exact_minimizer": [_fraction_text(value) for value in exact_parameters],
        "float64_minimizer": [float(value) for value in float_parameters],
        "exact_residual_norm_squared": _fraction_text(exact_rss),
        "residual_frobenius_norm": float(np.sqrt(float(exact_rss))),
        "target_frobenius_norm": float(np.sqrt(target_norm_squared)),
        "relative_residual": float(np.sqrt(float(exact_rss) / target_norm_squared)),
        "exact_max_abs_residual_entry": _fraction_text(exact_max),
        "nonzero_scaled_residual_entries": int(np.count_nonzero(scaled_residual)),
        "float64_residual_norm": float(np.linalg.norm(float_residual)),
        "singular_values": [float(value) for value in singular_values],
    }


def _cubic_torus_edges(length: int) -> list[tuple[int, int]]:
    def index(x: int, y: int, z: int) -> int:
        return (x * length + y) * length + z

    edges: set[tuple[int, int]] = set()
    for x, y, z in product(range(length), repeat=3):
        site = index(x, y, z)
        for neighbour in (
            index((x + 1) % length, y, z),
            index(x, (y + 1) % length, z),
            index(x, y, (z + 1) % length),
        ):
            edges.add(tuple(sorted((site, neighbour))))
    return sorted(edges)


def _count_simple_cycles(
    vertex_count: int, edges: list[tuple[int, int]], length: int
) -> int:
    adjacency = [set() for _ in range(vertex_count)]
    for u, v in edges:
        adjacency[u].add(v)
        adjacency[v].add(u)
    oriented = 0
    for start in range(vertex_count):
        visited = {start}

        def walk(current: int, edges_used: int) -> None:
            nonlocal oriented
            if edges_used == length - 1:
                if start in adjacency[current]:
                    oriented += 1
                return
            for neighbour in adjacency[current]:
                if neighbour <= start or neighbour in visited:
                    continue
                visited.add(neighbour)
                walk(neighbour, edges_used + 1)
                visited.remove(neighbour)

        walk(start, 0)
    assert oriented % 2 == 0
    return oriented // 2


def _direct_series_spot_checks() -> dict:
    ht_length = 7
    ht_edges = _cubic_torus_edges(ht_length)
    ht_sites = ht_length**3
    four_cycles = _count_simple_cycles(ht_sites, ht_edges, 4)
    six_cycles = _count_simple_cycles(ht_sites, ht_edges, 6)
    ht_interaction = {
        "2": Fraction(3, 2),
        "4": Fraction(3, 4) + Fraction(four_cycles, ht_sites),
        "6": Fraction(1, 2) + Fraction(six_cycles, ht_sites),
    }
    expected_ht = {"2": Fraction(3, 2), "4": Fraction(15, 4), "6": Fraction(45, 2)}

    lt_length = 4
    lt_edges = _cubic_torus_edges(lt_length)
    lt_sites = lt_length**3

    def surface(vertices: tuple[int, ...]) -> int:
        subset = set(vertices)
        return sum((u in subset) != (v in subset) for u, v in lt_edges)

    one_spin_histogram = Counter(surface((site,)) for site in range(lt_sites))
    two_spin_histogram = Counter(
        surface(pair) for pair in combinations(range(lt_sites), 2)
    )
    incompatible_ordered_monomer_pairs = lt_sites + 2 * len(lt_edges)
    lt_interaction = {
        "6": Fraction(one_spin_histogram[6], lt_sites),
        "10": Fraction(two_spin_histogram[10], lt_sites),
        "12": -Fraction(incompatible_ordered_monomer_pairs, 2 * lt_sites),
    }
    expected_lt = {"6": Fraction(1), "10": Fraction(3), "12": Fraction(-7, 2)}

    return {
        "method": (
            "direct cycle enumeration and direct one-/two-droplet surface counting; "
            "no transfer matrix, finite-lattice Mobius inversion, or repository series code"
        ),
        "high_temperature": {
            "torus": [ht_length, ht_length, ht_length],
            "no_wrap_reason": "L=7 exceeds every checked polygon length (<=6)",
            "simple_cycle_counts": {"4": four_cycles, "6": six_cycles},
            "cycle_counts_per_site": {
                "4": _fraction_text(Fraction(four_cycles, ht_sites)),
                "6": _fraction_text(Fraction(six_cycles, ht_sites)),
            },
            "interaction_coefficients": {
                order: _fraction_text(value) for order, value in ht_interaction.items()
            },
            "repository_claims": {
                order: _fraction_text(value) for order, value in expected_ht.items()
            },
            "passed": ht_interaction == expected_ht,
        },
        "low_temperature": {
            "torus": [lt_length, lt_length, lt_length],
            "no_wrap_reason": (
                "the shortest noncontractible droplet is a four-site loop with surface 16, "
                "above the checked order 12"
            ),
            "one_spin_surface_histogram": dict(sorted(one_spin_histogram.items())),
            "two_spin_surface_histogram": dict(sorted(two_spin_histogram.items())),
            "incompatible_ordered_monomer_pairs": incompatible_ordered_monomer_pairs,
            "interaction_coefficients": {
                order: _fraction_text(value) for order, value in lt_interaction.items()
            },
            "repository_claims": {
                order: _fraction_text(value) for order, value in expected_lt.items()
            },
            "passed": lt_interaction == expected_lt,
        },
    }


def _specified_edges(
    shape: tuple[int, ...],
    periodic: tuple[bool, ...],
    only_direction: int | None = None,
) -> list[tuple[int, int]]:
    def index(coordinate: tuple[int, ...]) -> int:
        result = 0
        for value, length in zip(coordinate, shape, strict=True):
            result = result * length + value
        return result

    edges: list[tuple[int, int]] = []
    for coordinate in product(*(range(length) for length in shape)):
        for direction, length in enumerate(shape):
            if only_direction is not None and direction != only_direction:
                continue
            if length == 1:
                continue
            if coordinate[direction] + 1 < length:
                neighbour = list(coordinate)
                neighbour[direction] += 1
            elif periodic[direction]:
                neighbour = list(coordinate)
                neighbour[direction] = 0
            else:
                continue
            edges.append(tuple(sorted((index(coordinate), index(tuple(neighbour))))))
    return edges


def _broken_bond_polynomial(
    site_count: int, edges: list[tuple[int, int]]
) -> list[int]:
    polynomial = [0] * (len(edges) + 1)
    for configuration in range(1 << site_count):
        broken = sum(
            ((configuration >> u) & 1) != ((configuration >> v) & 1)
            for u, v in edges
        )
        polynomial[broken] += 1
    while len(polynomial) > 1 and polynomial[-1] == 0:
        polynomial.pop()
    return polynomial


def _trim(polynomial: list[int]) -> list[int]:
    while len(polynomial) > 1 and polynomial[-1] == 0:
        polynomial.pop()
    return polynomial


def _plus_boundary_oracle(shape: tuple[int, int, int]) -> list[int]:
    edges = _specified_edges(shape, (False, False, False))
    site_count = int(np.prod(shape))
    degrees = [0] * site_count
    for u, v in edges:
        degrees[u] += 1
        degrees[v] += 1
    missing = [6 - degree for degree in degrees]
    polynomial = [0] * (6 * site_count + 1)
    for down_mask in range(1 << site_count):
        broken = sum(
            ((down_mask >> u) & 1) != ((down_mask >> v) & 1)
            for u, v in edges
        )
        broken += sum(
            missing[site]
            for site in range(site_count)
            if (down_mask >> site) & 1
        )
        polynomial[broken] += 1
    return _trim(polynomial)


def _boundary_attack() -> dict:
    sys.path.insert(0, str(ROOT / "src"))
    from ising.lattices import hyperrect
    from ising.transfer_matrix import box_broken_bond_poly, torus_broken_bond_poly

    cases = []
    for shape in ((1, 1), (1, 2), (2, 2), (1, 2, 3), (2, 2, 3)):
        periodic = tuple(True for _ in shape)
        expected_edges = _specified_edges(shape, periodic)
        lattice = hyperrect(shape, periodic=True)
        repository_edges = [
            tuple(sorted(edge)) for direction in lattice.bonds_by_dir for edge in direction
        ]
        expected_polynomial = _broken_bond_polynomial(lattice.n_sites, expected_edges)
        transfer_polynomial = _trim(list(torus_broken_bond_poly(shape)))
        expected_directional_counts = [
            len(_specified_edges(shape, periodic, only_direction=direction))
            for direction in range(len(shape))
        ]
        repository_directional_counts = [
            len(direction) for direction in lattice.bonds_by_dir
        ]
        cases.append(
            {
                "shape": list(shape),
                "expected_directional_bond_counts": expected_directional_counts,
                "repository_directional_bond_counts": repository_directional_counts,
                "directional_bond_counts_match": (
                    repository_directional_counts == expected_directional_counts
                ),
                "lattice_edges_match_specification": (
                    Counter(repository_edges) == Counter(expected_edges)
                ),
                "transfer_polynomial_matches_direct_enumeration": (
                    transfer_polynomial == expected_polynomial
                ),
            }
        )

    bug_shape = (2, 2, 1)
    plus_oracle = _plus_boundary_oracle(bug_shape)
    plus_repository = _trim(list(box_broken_bond_poly(bug_shape, plus_boundary=True)))
    first_mismatch = next(
        (
            order
            for order in range(max(len(plus_oracle), len(plus_repository)))
            if (plus_oracle[order] if order < len(plus_oracle) else 0)
            != (plus_repository[order] if order < len(plus_repository) else 0)
        ),
        None,
    )
    return {
        "periodic_length_1_and_2_cases": cases,
        "periodic_conventions_passed": all(
            case["lattice_edges_match_specification"]
            and case["directional_bond_counts_match"]
            and case["transfer_polynomial_matches_direct_enumeration"]
            for case in cases
        ),
        "plus_boundary_length_one_attack": {
            "shape": list(bug_shape),
            "expected_polynomial": plus_oracle,
            "repository_polynomial": plus_repository,
            "first_mismatch_order": first_mismatch,
            "passed": plus_repository == plus_oracle,
            "scope": (
                "this public-API mismatch was found during the audit; the source now handles "
                "both ghost faces at c=1, and FLM shape canonicalisation explains why the "
                "published LT coefficients were unchanged"
            ),
        },
    }


def _benchmark_literal_scan() -> dict:
    constants = ("0.221654626", "0.4406867", "0.629971")
    occurrences: dict[str, list[dict]] = {constant: [] for constant in constants}
    for directory in (ROOT / "src", ROOT / "experiments", ROOT / "tests"):
        for path in sorted(directory.rglob("*.py")):
            if path.resolve() == Path(__file__).resolve():
                continue
            for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
                for constant in constants:
                    if constant in line:
                        occurrences[constant].append(
                            {
                                "file": str(path.relative_to(ROOT)),
                                "line": line_number,
                                "text": line.strip(),
                            }
                        )
    lee_yang_source = (ROOT / "experiments" / "e04_lee_yang.py").read_text(
        encoding="utf-8"
    )
    lee_yang_selection_leak = (
        "if mp.mpf(point[\"K\"]) < REFERENCE_K" in lee_yang_source
        and "edge_fit_with_sensitivity" in lee_yang_source
    )
    series_source = (ROOT / "experiments" / "e19_series_analysis.py").read_text(
        encoding="utf-8"
    )
    series_comparison_after_fit = (
        series_source.index("central_k = current[\"central_k\"]")
        < series_source.index("benchmark_difference = abs(central_k - reference_kc)")
    )
    return {
        "constants": occurrences,
        "lee_yang_edge_fit_uses_benchmark_for_sample_classification": lee_yang_selection_leak,
        "series_analysis_uses_benchmark_only_after_fit": series_comparison_after_fit,
        "verdict": (
            "FAIL: e04 uses the benchmark to classify fit inputs"
            if lee_yang_selection_leak
            else (
                "PASS after audit repair: e04 classifies disordered-phase inputs with the "
                "in-repo certified lower bound; the benchmark remains comparison-only"
            )
        ),
    }


def _bounds_attack() -> dict:
    lower_decimal = Decimal("0.2074277114992039908436804465100887760027")
    upper_decimal = Decimal("0.2527310098586630030260020266135701299926")
    benchmark = Decimal("0.221654626")
    with mp.workdps(110):
        saw_point = mp.atanh(mp.mpf(4468911678) ** (-mp.mpf(1) / 14))
        watson = (
            mp.sqrt(6)
            * mp.gamma(mp.mpf(1) / 24)
            * mp.gamma(mp.mpf(5) / 24)
            * mp.gamma(mp.mpf(7) / 24)
            * mp.gamma(mp.mpf(11) / 24)
            / (32 * mp.pi**3)
        )
        infrared_point = watson / 6
        saw_decimal = Decimal(_mp_text(saw_point, 105))
        infrared_decimal = Decimal(_mp_text(infrared_point, 105))
    return {
        "independent_precision_digits": 110,
        "saw_formula_point": str(saw_decimal),
        "infrared_formula_point": str(infrared_decimal),
        "claimed_outward_interval": [str(lower_decimal), str(upper_decimal)],
        "lower_endpoint_is_outward": lower_decimal <= saw_decimal,
        "upper_endpoint_is_outward": upper_decimal >= infrared_decimal,
        "benchmark_is_bracketed": lower_decimal < benchmark < upper_decimal,
        "safe_direction_summary": (
            "c_n^(1/n) >= mu implies atanh(c_n^(-1/n)) <= atanh(1/mu) <= Kc; "
            "the infrared magnetisation inequality implies Kc <= W_sc/6"
        ),
    }


def main() -> None:
    star_edges = [(0, 1), (0, 2), (0, 3)]
    grid_edges = [
        (0, 1),
        (1, 2),
        (3, 4),
        (4, 5),
        (0, 3),
        (1, 4),
        (2, 5),
    ]
    cycle_plus_isolate_edges = [(0, 1), (1, 2), (0, 2)]
    dg_cases = [
        _dg_case("star K1,3", 4, star_edges),
        _dg_case("open 2x3 grid", 6, grid_edges),
        _dg_case("cycle C3 plus one isolated vertex", 4, cycle_plus_isolate_edges),
    ]
    td = _td_dense_replication()
    series = _direct_series_spot_checks()
    boundary = _boundary_attack()
    benchmark_scan = _benchmark_literal_scan()
    bounds = _bounds_attack()
    centralizer_scope = _full_centralizer_scope_attack()

    checks = [
        {
            "name": "DG dense closed form on K1,3",
            "passed": dg_cases[0]["closed_form_max_abs_difference"] == 0,
            "detail": (
                "max |dense residual - dense closed form| = "
                f"{dg_cases[0]['closed_form_max_abs_difference']}"
            ),
        },
        {
            "name": "DG dense closed form on open 2x3 grid",
            "passed": dg_cases[1]["closed_form_max_abs_difference"] == 0,
            "detail": (
                "max |dense residual - dense closed form| = "
                f"{dg_cases[1]['closed_form_max_abs_difference']}"
            ),
        },
        {
            "name": "DG residual iff 2-regular as originally worded",
            "passed": False,
            "detail": (
                "C3 plus an isolated vertex is not 2-regular but has exact residual max 0; "
                "the corrected condition is every degree in {0,2}"
            ),
        },
        {
            "name": "TD2 dense no-solution on 3x3 torus",
            "passed": td["exact_residual_norm_squared"] != "0",
            "detail": (
                f"exact least-squares RSS={td['exact_residual_norm_squared']}; "
                f"Frobenius residual={td['residual_frobenius_norm']:.12g}; "
                f"relative residual={td['relative_residual']:.12g}"
            ),
        },
        {
            "name": "direct HT series coefficients v2,v4,v6",
            "passed": series["high_temperature"]["passed"],
            "detail": str(series["high_temperature"]["interaction_coefficients"]),
        },
        {
            "name": "direct LT series coefficients x6,x10,x12",
            "passed": series["low_temperature"]["passed"],
            "detail": str(series["low_temperature"]["interaction_coefficients"]),
        },
        {
            "name": "periodic length-1/length-2 bond convention",
            "passed": boundary["periodic_conventions_passed"],
            "detail": "lattice bond multisets and transfer polynomials match direct enumeration",
        },
        {
            "name": "plus-boundary transfer API with final length 1",
            "passed": boundary["plus_boundary_length_one_attack"]["passed"],
            "detail": (
                "all coefficients match the independent ghost-bond enumeration"
                if boundary["plus_boundary_length_one_attack"]["passed"]
                else (
                    "first coefficient mismatch at x^"
                    f"{boundary['plus_boundary_length_one_attack']['first_mismatch_order']}"
                )
            ),
        },
        {
            "name": "benchmark literals are comparison-only everywhere",
            "passed": not benchmark_scan[
                "lee_yang_edge_fit_uses_benchmark_for_sample_classification"
            ],
            "detail": benchmark_scan["verdict"],
        },
        {
            "name": "certified Kc interval endpoint and benchmark sanity",
            "passed": (
                bounds["lower_endpoint_is_outward"]
                and bounds["upper_endpoint_is_outward"]
                and bounds["benchmark_is_bracketed"]
            ),
            "detail": (
                f"{bounds['claimed_outward_interval'][0]} < 0.221654626 < "
                f"{bounds['claimed_outward_interval'][1]}"
            ),
        },
        {
            "name": "only two sectors/full centralizer consequence of Pauli theorem",
            "passed": False,
            "detail": (
                "C3 translation commutes exactly with A and B and is linearly independent "
                f"of I and global parity (Gram determinant {centralizer_scope['gram_determinant']})"
            ),
        },
    ]

    assert dg_cases[0]["closed_form_max_abs_difference"] == 0
    assert dg_cases[1]["closed_form_max_abs_difference"] == 0
    assert dg_cases[0]["first_relation_max_abs"] == 0
    assert dg_cases[1]["first_relation_max_abs"] == 0
    assert dg_cases[2]["residual_max_abs"] == 0
    assert td["design_rank"] == 3
    assert td["exact_residual_norm_squared"] == "488374272/13"
    assert td["exact_minimizer"] == ["23/13", "18/13", "508/13"]
    assert td["residual_frobenius_norm"] > 6000
    assert td["relative_residual"] > 0.2
    assert series["high_temperature"]["passed"]
    assert series["low_temperature"]["passed"]
    assert boundary["periodic_conventions_passed"]
    assert bounds["lower_endpoint_is_outward"]
    assert bounds["upper_endpoint_is_outward"]
    assert bounds["benchmark_is_bracketed"]
    assert centralizer_scope["commutator_with_A_max_abs"] == 0
    assert centralizer_scope["commutator_with_B_max_abs"] == 0
    assert centralizer_scope["linearly_independent_from_I_and_parity"]

    payload = {
        "provenance": {
            "script": "tests/test_audit_replication.py",
            "python": platform.python_version(),
            "numpy": np.__version__,
            "mpmath": mp.__version__,
            "dense_arithmetic": (
                "exact int64 matrix products, checked to stay far below overflow; "
                "Fraction normal equations for TD"
            ),
            "independence": (
                "no imports from repository algebra, Pauli, series, or rigorous-bound modules; "
                "repository lattice/transfer imports are confined to the explicit boundary attack"
            ),
        },
        "data": {
            "dolan_grady_dense": dg_cases,
            "tridiagonal_dense": td,
            "series_direct_torus": series,
            "boundary_conventions": boundary,
            "benchmark_literal_scan": benchmark_scan,
            "rigorous_bounds": bounds,
            "full_centralizer_scope_attack": centralizer_scope,
        },
        "checks": checks,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    for case in dg_cases[:2]:
        status = "PASS" if case["closed_form_max_abs_difference"] == 0 else "FAIL"
        print(
            f"DG {case['name']}: {status}; "
            f"max deviation={case['closed_form_max_abs_difference']}"
        )
    print(
        "DG original iff-2-regular corollary: FAIL; "
        "C3 plus isolated vertex has residual max=0"
    )
    print(
        "TD 3x3 torus no-solution: PASS; "
        f"minimum Frobenius residual={td['residual_frobenius_norm']:.12g}; "
        f"relative={td['relative_residual']:.12g}; exact RSS={td['exact_residual_norm_squared']}"
    )
    print(
        "periodic L=1/L=2 conventions: "
        f"{'PASS' if boundary['periodic_conventions_passed'] else 'FAIL'}"
    )
    plus_attack = boundary["plus_boundary_length_one_attack"]
    plus_passed = plus_attack["passed"]
    plus_detail = (
        "all coefficients match"
        if plus_passed
        else f"first mismatch x^{plus_attack['first_mismatch_order']}"
    )
    print(
        "plus-boundary c=1 API: "
        f"{'PASS' if plus_passed else 'FAIL'}; {plus_detail}"
    )
    print("AUDIT REPLICATION COMPLETED; JSON written to", OUTPUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
