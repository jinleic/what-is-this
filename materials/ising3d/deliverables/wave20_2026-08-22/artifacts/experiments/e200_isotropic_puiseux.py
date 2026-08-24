"""Exact high-temperature leading terms on the physical isotropic curve.

For a graph G, Q_G(y) counts spin configurations by aligned edges.  The
anti-reciprocal polynomial I_G(y)=Q_G(y)-y^m Q_G(1/y) is reconstructed in two
independent exact ways: spin enumeration and the Eulerian-subgraph expansion.
The computations here are small lemma tests; the all-size arguments are in the
proof note and are integrated by e202.
"""

from __future__ import annotations

from collections import deque
from fractions import Fraction
from math import comb

Polynomial = list[int]
Edge = tuple[int, int]

SAMPLE_GRAPHS: dict[str, tuple[int, tuple[Edge, ...]]] = {
    "path_4": (4, ((0, 1), (1, 2), (2, 3))),
    "square": (4, ((0, 1), (1, 2), (2, 3), (3, 0))),
    "triangle": (3, ((0, 1), (1, 2), (2, 0))),
    "paw": (4, ((0, 1), (1, 2), (2, 0), (0, 3))),
    "lollipop_2": (5, ((0, 1), (1, 2), (2, 0), (0, 3), (3, 4))),
    "cycle_5": (5, ((0, 1), (1, 2), (2, 3), (3, 4), (4, 0))),
    "complete_4": (4, tuple((left, right) for left in range(4) for right in range(left + 1, 4))),
}


def _check(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def trim(poly: Polynomial) -> Polynomial:
    result = list(poly)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result or [0]


def poly_add(left: Polynomial, right: Polynomial) -> Polynomial:
    size = max(len(left), len(right))
    result = [0] * size
    for index in range(size):
        result[index] = (left[index] if index < len(left) else 0) + (
            right[index] if index < len(right) else 0
        )
    return trim(result)


def poly_sub(left: Polynomial, right: Polynomial) -> Polynomial:
    return poly_add(left, [-coefficient for coefficient in right])


def poly_scale(poly: Polynomial, scalar: int) -> Polynomial:
    return trim([scalar * coefficient for coefficient in poly])


def poly_mul(left: Polynomial, right: Polynomial) -> Polynomial:
    if left == [0] or right == [0]:
        return [0]
    result = [0] * (len(left) + len(right) - 1)
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            result[left_degree + right_degree] += left_coefficient * right_coefficient
    return trim(result)


def poly_pow(poly: Polynomial, exponent: int) -> Polynomial:
    if exponent < 0:
        raise ValueError("polynomial exponent must be nonnegative")
    result = [1]
    base = trim(poly)
    power = exponent
    while power:
        if power & 1:
            result = poly_mul(result, base)
        power >>= 1
        if power:
            base = poly_mul(base, base)
    return result


def linear_power(constant: int, exponent: int) -> Polynomial:
    """Return coefficients of (y + constant)^exponent in ascending order."""
    return [comb(exponent, degree) * constant ** (exponent - degree) for degree in range(exponent + 1)]


def poly_eval(poly: Polynomial, value: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


def shift_from_one(poly: Polynomial) -> Polynomial:
    """Return coefficients in s of poly(1+s)."""
    result = [0] * len(poly)
    for old_degree, coefficient in enumerate(poly):
        for new_degree in range(old_degree + 1):
            result[new_degree] += coefficient * comb(old_degree, new_degree)
    return trim(result)


def fraction_record(value: Fraction) -> dict[str, int]:
    return {"numerator": value.numerator, "denominator": value.denominator}


def validate_graph(sites: int, edges: tuple[Edge, ...]) -> None:
    normalized: set[Edge] = set()
    for left, right in edges:
        if not (0 <= left < sites and 0 <= right < sites) or left == right:
            raise ValueError(f"invalid edge {(left, right)} on {sites} sites")
        edge = (min(left, right), max(left, right))
        if edge in normalized:
            raise ValueError(f"duplicate edge {edge}")
        normalized.add(edge)


def connected_components(sites: int, edges: tuple[Edge, ...]) -> int:
    adjacency = [[] for _ in range(sites)]
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    seen = [False] * sites
    components = 0
    for start in range(sites):
        if seen[start]:
            continue
        components += 1
        seen[start] = True
        queue = deque([start])
        while queue:
            vertex = queue.popleft()
            for neighbor in adjacency[vertex]:
                if not seen[neighbor]:
                    seen[neighbor] = True
                    queue.append(neighbor)
    return components


def is_bipartite(sites: int, edges: tuple[Edge, ...]) -> bool:
    adjacency = [[] for _ in range(sites)]
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    colors: list[int | None] = [None] * sites
    for start in range(sites):
        if colors[start] is not None:
            continue
        colors[start] = 0
        queue = deque([start])
        while queue:
            vertex = queue.popleft()
            for neighbor in adjacency[vertex]:
                if colors[neighbor] is None:
                    colors[neighbor] = 1 - int(colors[vertex])
                    queue.append(neighbor)
                elif colors[neighbor] == colors[vertex]:
                    return False
    return True


def alignment_counts(sites: int, edges: tuple[Edge, ...]) -> Polynomial:
    """Return Q_G(y)=sum_sigma y^(number of aligned edges)."""
    counts = [0] * (len(edges) + 1)
    for state in range(1 << sites):
        aligned = sum(((state >> left) & 1) == ((state >> right) & 1) for left, right in edges)
        counts[aligned] += 1
    return counts


def eulerian_counts(sites: int, edges: tuple[Edge, ...]) -> Polynomial:
    """Count even-degree edge subsets by cardinality."""
    counts = [0] * (len(edges) + 1)
    endpoint_masks = [(1 << left) ^ (1 << right) for left, right in edges]
    for edge_mask in range(1 << len(edges)):
        boundary = 0
        for index, endpoint_mask in enumerate(endpoint_masks):
            if (edge_mask >> index) & 1:
                boundary ^= endpoint_mask
        if boundary == 0:
            counts[edge_mask.bit_count()] += 1
    return counts


def reverse_alignment(counts: Polynomial) -> Polynomial:
    return list(reversed(counts))


def trace_skew_from_spins(counts: Polynomial) -> Polynomial:
    """Return I_G(y)=Q_G(y)-y^m Q_G(1/y)."""
    return trim([counts[index] - counts[-1 - index] for index in range(len(counts))])


def trace_skew_from_eulerian(
    sites: int, edge_count: int, even_subgraph_counts: Polynomial
) -> Polynomial:
    """Rebuild I_G from the exact high-temperature Eulerian expansion."""
    result = [Fraction(0)] * (edge_count + 1)
    scale = Fraction(2 ** (sites + 1), 2**edge_count)
    for cardinality, count in enumerate(even_subgraph_counts):
        if cardinality % 2 == 0 or count == 0:
            continue
        term = poly_mul(linear_power(-1, cardinality), linear_power(1, edge_count - cardinality))
        for degree, coefficient in enumerate(term):
            result[degree] += scale * count * coefficient
    if any(coefficient.denominator != 1 for coefficient in result):
        raise AssertionError("Eulerian reconstruction did not cancel to integer coefficients")
    return trim([coefficient.numerator for coefficient in result])


def physical_clear(poly: Polynomial, denominator_power: int) -> Polynomial:
    """Return (2t)^d poly((1+t^2)/(2t)), with d=denominator_power."""
    if denominator_power < len(trim(poly)) - 1:
        raise ValueError("denominator power is smaller than polynomial degree")
    result = [0]
    for y_degree, coefficient in enumerate(poly):
        if coefficient == 0:
            continue
        one_plus_t_squared = [0] * (2 * y_degree + 1)
        for choice in range(y_degree + 1):
            one_plus_t_squared[2 * choice] = comb(y_degree, choice)
        monomial_degree = denominator_power - y_degree
        two_t_power = [0] * (monomial_degree + 1)
        two_t_power[monomial_degree] = 2**monomial_degree
        term = poly_scale(poly_mul(one_plus_t_squared, two_t_power), coefficient)
        result = poly_add(result, term)
    return trim(result)


def local_leading_data(trace_skew: Polynomial) -> tuple[int | None, int]:
    shifted = shift_from_one(trace_skew)
    for degree, coefficient in enumerate(shifted):
        if coefficient:
            return degree, coefficient
    return None, 0


def graph_row(label: str, sites: int, edges: tuple[Edge, ...]) -> tuple[dict[str, object], list[bool]]:
    validate_graph(sites, edges)
    alignments = alignment_counts(sites, edges)
    eulerian = eulerian_counts(sites, edges)
    spin_skew = trace_skew_from_spins(alignments)
    eulerian_skew = trace_skew_from_eulerian(sites, len(edges), eulerian)
    odd_cardinalities = [degree for degree, count in enumerate(eulerian) if degree % 2 and count]
    odd_girth = min(odd_cardinalities) if odd_cardinalities else None
    shortest_count = eulerian[odd_girth] if odd_girth is not None else 0
    valuation, leading = local_leading_data(spin_skew)
    expected_leading = 0 if odd_girth is None else 2 ** (sites + 1 - odd_girth) * shortest_count
    physical = physical_clear(spin_skew, len(edges))
    t0 = Fraction(1, 3)
    y0 = (1 + t0 * t0) / (2 * t0)
    value_at_y0 = poly_eval(spin_skew, y0)
    physical_value_at_t0 = poly_eval(physical, t0)
    components = connected_components(sites, edges)
    bipartite = is_bipartite(sites, edges)
    row = {
        "tag": "[COMPUTATION]",
        "label": label,
        "sites": sites,
        "edges": [list(edge) for edge in edges],
        "edge_count": len(edges),
        "connected_components": components,
        "bipartite": bipartite,
        "alignment_coefficients_ascending": alignments,
        "reverse_alignment_coefficients_ascending": reverse_alignment(alignments),
        "trace_skew_coefficients_ascending": spin_skew,
        "eulerian_counts_by_cardinality": eulerian,
        "odd_girth": odd_girth,
        "shortest_odd_cycle_count": shortest_count,
        "y_minus_one_valuation": valuation,
        "leading_coefficient_at_y_one": leading,
        "expected_leading_coefficient": expected_leading,
        "value_at_y_5_over_3": fraction_record(value_at_y0),
        "physical_cleared_coefficients_ascending": physical,
    }
    audits = [
        sum(alignments) == 1 << sites,
        sum(eulerian) == 1 << (len(edges) - sites + components),
        spin_skew == eulerian_skew,
        (odd_girth is None) == bipartite,
        valuation == odd_girth,
        leading == expected_leading,
        physical_value_at_t0 == (2 * t0) ** len(edges) * value_at_y0,
        (value_at_y0 == 0) == bipartite,
    ]
    return row, audits


def run_puiseux() -> dict[str, object]:
    checks: list[dict[str, object]] = []
    rows: list[dict[str, object]] = []
    row_by_label: dict[str, dict[str, object]] = {}
    for label, (sites, edges) in SAMPLE_GRAPHS.items():
        row, audits = graph_row(label, sites, edges)
        rows.append(row)
        row_by_label[label] = row
        _check(
            checks,
            f"{label}: spin and Eulerian constructions agree",
            all(audits[:3]),
            f"Q coefficients sum to {sum(row['alignment_coefficients_ascending'])}; cycle-space count={sum(row['eulerian_counts_by_cardinality'])}",
        )
        _check(
            checks,
            f"{label}: leading high-temperature germ",
            all(audits[3:6]),
            f"bipartite={row['bipartite']}, valuation={row['y_minus_one_valuation']}, leading={row['leading_coefficient_at_y_one']}",
        )
        _check(
            checks,
            f"{label}: physical-curve clearing and exact witness point",
            all(audits[6:]),
            f"I_G(5/3)={row['value_at_y_5_over_3']['numerator']}/{row['value_at_y_5_over_3']['denominator']}",
        )

    triangle = row_by_label["triangle"]
    paw = row_by_label["paw"]
    lollipop = row_by_label["lollipop_2"]
    triangle_skew = list(triangle["trace_skew_coefficients_ascending"])
    paw_skew = list(paw["trace_skew_coefficients_ascending"])
    lollipop_skew = list(lollipop["trace_skew_coefficients_ascending"])
    leaf_once = poly_mul(triangle_skew, [1, 1])
    leaf_twice = poly_mul(leaf_once, [1, 1])
    _check(
        checks,
        "leaf extension multiplies the trace skew by y+1",
        leaf_once == paw_skew and leaf_twice == lollipop_skew,
        "triangle -> paw -> two-edge lollipop checked coefficientwise",
    )

    t = Fraction(1, 3)
    x = (1 + t) / (1 - t)
    y = (1 + t * t) / (2 * t)
    curve_value = (y - 1) * x * x - (y + 1)
    _check(
        checks,
        "exact witness lies on the physical curve",
        x == 2 and y == Fraction(5, 3) and curve_value == 0,
        f"(x,y)=({x},{y}), F_phys={curve_value}",
    )

    data = {
        "tag": "[COMPUTATION]",
        "route": "high-temperature Taylor/Puiseux leading term of a conjugacy-invariant trace skew",
        "definitions": {
            "tag": "[LEMMA]",
            "alignment_polynomial": "Q_G(y)=sum_sigma y^{a_G(sigma)}",
            "reverse_polynomial": "Q_G^*(y)=y^m Q_G(1/y)=sum_sigma y^{m-a_G(sigma)}",
            "trace_skew": "I_G(y)=Q_G(y)-Q_G^*(y)",
            "eulerian_identity": "I_G(y)=2^{n+1-m} sum_{F Eulerian, |F| odd}(y-1)^{|F|}(y+1)^{m-|F|}",
            "leading_term": "if odd girth is g and there are c_g shortest odd cycles, I_G(1+s)=2^{n+1-g} c_g s^g+O(s^{g+1})",
            "physical_substitution": "y=q(t)=(1+t^2)/(2t); stored J_G(t)=(2t)^m I_G(q(t))",
        },
        "graph_rows": rows,
        "finite_scope": (
            "the rows are exact lemma tests only; the all-size bipartite/non-bipartite proof uses the "
            "Eulerian expansion and does not extrapolate from this table"
        ),
    }
    return {"data": data, "checks": checks}


def main() -> int:
    result = run_puiseux()
    for item in result["checks"]:
        print(f"[{'PASS' if item['passed'] else 'FAIL'}] {item['name']}: {item['detail']}")
    if all(bool(item["passed"]) for item in result["checks"]):
        print("PASS e200 isotropic Puiseux leading terms")
        return 0
    print("FAIL e200 isotropic Puiseux leading terms")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
