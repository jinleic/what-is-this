"""Exact random-current single-edge domination LP for the upper endpoint.

For two independent sourceless currents, combine the currents edge by edge and
retain each colour parity plus whether the combined edge is present.  At
v=tanh(K), all five local weights are rational.  The smallest conditional
stochastic-domination LP has exact optimum p=v^2.  At v=6/25 it is exactly
infeasible against a self-contained planar Peierls threshold p=5/6.

This module is a producer component; e177 writes the combined JSON artifact.
"""

from __future__ import annotations

import itertools
import time
from collections import deque
from fractions import Fraction
from typing import Any

SCRIPT = "experiments/e176_upper_endpoint4_current.py"
V_TARGET = Fraction(6, 25)
PERCOLATION_THRESHOLD = Fraction(5, 6)
PROCESS_BUDGET_SECONDS = 30.0


def fstr(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def _record(checks: list[dict[str, Any]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def boundary(vertices: int, edges: list[tuple[int, int]], bits: tuple[int, ...]) -> frozenset[int]:
    parity = [0] * vertices
    for bit, (u, w) in zip(bits, edges):
        if bit:
            parity[u] ^= 1
            parity[w] ^= 1
    return frozenset(i for i, value in enumerate(parity) if value)


def parity_polynomial(
    vertices: int, edges: list[tuple[int, int]], sources: frozenset[int], v: Fraction
) -> Fraction:
    total = Fraction(0)
    for mask in range(1 << len(edges)):
        bits = tuple((mask >> i) & 1 for i in range(len(edges)))
        if boundary(vertices, edges, bits) == sources:
            total += v ** sum(bits)
    return total


def connected(vertices: int, edges: list[tuple[int, int]], present: tuple[bool, ...], source: int, target: int) -> bool:
    adjacency = [[] for _ in range(vertices)]
    for is_present, (u, w) in zip(present, edges):
        if is_present:
            adjacency[u].append(w)
            adjacency[w].append(u)
    seen = {source}
    queue = deque([source])
    while queue:
        u = queue.popleft()
        for w in adjacency[u]:
            if w not in seen:
                seen.add(w)
                queue.append(w)
    return target in seen


def local_states(v: Fraction) -> list[dict[str, Any]]:
    """Five exact states of one edge in a pair of coloured currents.

    Both currents are scaled by one cosh(K) factor on this edge.  State 00z is
    absent; 00p is positive with both colour parities even.  The other three
    states are necessarily present.
    """
    return [
        {"name": "00_zero", "parity": (0, 0), "present": False, "weight": 1 - v * v},
        {"name": "00_positive", "parity": (0, 0), "present": True, "weight": v * v},
        {"name": "11_positive", "parity": (1, 1), "present": True, "weight": v * v},
        {"name": "10_positive", "parity": (1, 0), "present": True, "weight": v},
        {"name": "01_positive", "parity": (0, 1), "present": True, "weight": v},
    ]


def enumerate_double_current_plaquette(v: Fraction) -> dict[str, Fraction | int]:
    """Exact five-state enumeration on the four-edge plaquette C4."""
    vertices = 4
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    states = local_states(v)
    started = time.process_time()
    total = Fraction(0)
    connection = Fraction(0)
    accepted = 0
    for indices in itertools.product(range(len(states)), repeat=len(edges)):
        if time.process_time() - started > PROCESS_BUDGET_SECONDS:
            raise RuntimeError("NON-DECISIVE: process-time budget expired in 5^4 plaquette enumeration")
        chosen = [states[index] for index in indices]
        bits1 = tuple(state["parity"][0] for state in chosen)
        bits2 = tuple(state["parity"][1] for state in chosen)
        if boundary(vertices, edges, bits1) or boundary(vertices, edges, bits2):
            continue
        weight = Fraction(1)
        for state in chosen:
            weight *= state["weight"]
        accepted += 1
        total += weight
        present = tuple(bool(state["present"]) for state in chosen)
        if connected(vertices, edges, present, 0, 2):
            connection += weight
    return {
        "accepted_state_assignments": accepted,
        "total_weight": total,
        "connection_0_2_weight": connection,
        "process_seconds_numerator_bound": int((time.process_time() - started) * 1_000_000) + 1,
    }


def conditional_open_table(v: Fraction) -> list[dict[str, str]]:
    states = local_states(v)
    rows = []
    for residual in ((0, 0), (0, 1), (1, 0), (1, 1)):
        eligible = [state for state in states if state["parity"] == residual]
        total = sum((state["weight"] for state in eligible), Fraction(0))
        open_weight = sum((state["weight"] for state in eligible if state["present"]), Fraction(0))
        rows.append(
            {
                "residual_parity": f"{residual[0]}{residual[1]}",
                "total_weight": fstr(total),
                "open_weight": fstr(open_weight),
                "conditional_open_probability": fstr(open_weight / total),
            }
        )
    return rows


def build_current_data() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Build the exact random-current LP obstruction and finite switching audit."""
    checks: list[dict[str, Any]] = []
    v = V_TARGET
    states = local_states(v)
    state_weight_sum = sum((state["weight"] for state in states), Fraction(0))
    expected_weight_sum = 1 + 2 * v + v * v
    _record(
        checks,
        "current_local_weight_partition",
        state_weight_sum == expected_weight_sum,
        f"sum={fstr(state_weight_sum)}=(1+v)^2",
    )

    conditional = conditional_open_table(v)
    conditional_values = {
        row["residual_parity"]: Fraction(row["conditional_open_probability"]) for row in conditional
    }
    expected_conditional = {"00": v * v, "01": Fraction(1), "10": Fraction(1), "11": Fraction(1)}
    _record(
        checks,
        "current_conditional_probabilities",
        conditional_values == expected_conditional,
        f"conditional={{{', '.join(f'{key}:{fstr(value)}' for key, value in conditional_values.items())}}}",
    )

    plaquette = enumerate_double_current_plaquette(v)
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    p_empty = parity_polynomial(4, edges, frozenset(), v)
    p_02 = parity_polynomial(4, edges, frozenset({0, 2}), v)
    total_expected = p_empty * p_empty
    connection_expected = p_02 * p_02
    _record(
        checks,
        "current_plaquette_total",
        plaquette["total_weight"] == total_expected,
        f"enumerated={fstr(plaquette['total_weight'])}, P_empty^2={fstr(total_expected)}",
    )
    _record(
        checks,
        "switching_identity_plaquette",
        plaquette["connection_0_2_weight"] == connection_expected,
        f"event={fstr(plaquette['connection_0_2_weight'])}, P_02^2={fstr(connection_expected)}",
    )
    switching_probability = connection_expected / total_expected
    correlation = p_02 / p_empty
    _record(
        checks,
        "switching_probability_is_correlation_squared",
        switching_probability == correlation * correlation,
        f"P(0<->2)={fstr(switching_probability)}=G(0,2)^2",
    )

    # Exact one-variable LP: maximize p subject to p>=0 and p no larger than
    # every conditional open probability.  The 00 row is the exact dual bottleneck.
    p_optimum = min(conditional_values.values())
    _record(
        checks,
        "current_lp_exact_optimum",
        p_optimum == v * v == Fraction(36, 625),
        f"p*={fstr(p_optimum)}, dual bottleneck=residual 00",
    )

    infeasibility_gap = PERCOLATION_THRESHOLD - p_optimum
    _record(
        checks,
        "current_percolation_challenge_infeasible",
        infeasibility_gap == Fraction(2909, 3750) > 0,
        f"5/6-p*={fstr(infeasibility_gap)}",
    )

    # Self-contained square-plane Peierls estimate at p=5/6.  With q=1/6,
    # at most 4*n*3^(n-1) length-n closed dual circuits surround a vertex.
    # Sum_{n>=4} n*r^n = r^4*(4-3r)/(1-r)^2, r=3q=1/2.
    q_closed = 1 - PERCOLATION_THRESHOLD
    r = 3 * q_closed
    tail_n_rn = r**4 * (4 - 3 * r) / (1 - r) ** 2
    circuit_union_bound = Fraction(4, 3) * tail_n_rn
    theta_lower = 1 - circuit_union_bound
    _record(
        checks,
        "planar_peierls_threshold",
        r == Fraction(1, 2)
        and circuit_union_bound == Fraction(5, 6)
        and theta_lower == Fraction(1, 6) > 0,
        f"dual-circuit union bound={fstr(circuit_union_bound)}, theta>={fstr(theta_lower)}",
    )

    data = {
        "class_name": "RC2-single-edge-conditional-domination",
        "scope": (
            "two independent sourceless random currents; only the one-edge conditional support law "
            "is used to dominate iid bonds. Multi-edge block events, source-dependent currents, and "
            "other switching inequalities are outside the class"
        ),
        "target_v": fstr(v),
        "local_state_derivation": {
            "normalization": "divide each current edge weight by cosh(K), v=tanh(K)",
            "states": [
                {
                    "name": state["name"],
                    "colour_parity": list(state["parity"]),
                    "present": state["present"],
                    "weight": fstr(state["weight"]),
                }
                for state in states
            ],
            "identities": [
                "00 absent: sech(K)^2=1-v^2",
                "00 positive even-even: 1-sech(K)^2=v^2",
                "11 positive odd-odd: tanh(K)^2=v^2",
                "10 and 01 positive: tanh(K)=v",
            ],
            "conditional_open_table": conditional,
        },
        "finite_switching_evaluation": {
            "graph": "four-cycle C4 embedded as a square plaquette of Z^3",
            "vertices": 4,
            "edges": [[u, w] for u, w in edges],
            "accepted_double_current_state_assignments": plaquette["accepted_state_assignments"],
            "P_empty": fstr(p_empty),
            "P_02": fstr(p_02),
            "total_double_current_weight": fstr(plaquette["total_weight"]),
            "connection_0_2_weight": fstr(plaquette["connection_0_2_weight"]),
            "connection_probability": fstr(switching_probability),
            "two_point_correlation": fstr(correlation),
            "identity": "P^{empty,empty}(0 connected 2)=<sigma_0 sigma_2>^2",
            "process_time_microseconds_strict_upper": plaquette["process_seconds_numerator_bound"],
        },
        "domination_lp": {
            "variable": "p",
            "objective": "maximize p",
            "constraints": ["p>=0", "p<=v^2 (residual colour parity 00)", "p<=1 (rows 01,10,11)"],
            "exact_optimum": fstr(p_optimum),
            "primal_witness": fstr(p_optimum),
            "dual_certificate": "the residual-00 row alone gives p<=v^2=36/625",
            "certificate_threshold": fstr(PERCOLATION_THRESHOLD),
            "exact_infeasibility_gap": fstr(infeasibility_gap),
        },
        "self_contained_percolation_implication": {
            "plane": "an embedded square-lattice plane in Z^3",
            "iid_open_threshold": fstr(PERCOLATION_THRESHOLD),
            "closed_probability": fstr(q_closed),
            "circuit_count_bound": "number of length-n closed dual circuits surrounding a vertex <=4*n*3^(n-1)",
            "dual_circuit_union_bound": fstr(circuit_union_bound),
            "percolation_probability_lower": fstr(theta_lower),
            "correlation_lower_if_feasible": fstr(theta_lower),
            "implication": (
                "p>=5/6 gives plane percolation by the displayed Peierls sum; stochastic domination "
                "and switching give a uniform Ising two-point floor, hence K_c<=atanh(6/25)"
            ),
        },
    }
    return data, checks


def main() -> int:
    _, checks = build_current_data()
    for check in checks:
        print(f"{'PASS' if check['passed'] else 'FAIL'} {check['name']} -- {check['detail']}")
    passed = all(check["passed"] for check in checks)
    print("PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
