"""Independent exact verifier for results/bounds/upper_endpoint4.json.

The verifier imports neither producer component.  It re-derives the endpoint
comparison with a different coarse atanh bound, solves both finite optimization
problems analytically over Fraction, reconstructs the C4 switching identity by
a dynamic program rather than the producer's Cartesian enumeration, and checks
every declared universal-claim coverage row.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from collections import defaultdict, deque
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "bounds" / "upper_endpoint4.json"
UPPER_INFRARED = ROOT / "results" / "bounds" / "upper_infrared.json"
V = Fraction(6, 25)
P0 = Fraction(5, 6)
EXPECTED_COVERAGE = {
    "TARGET_BELOW_INCUMBENT",
    "PM_PARSEVAL_BOUNDED",
    "PM_SDP_OPTIMUM",
    "PM_POSITIVE_FLOOR_INFEASIBLE",
    "RC_EDGE_WEIGHTS",
    "RC_SWITCHING_C4",
    "RC_DOMINATION_LP",
    "PERC_PLANAR_5_6",
}

FAILURES: list[str] = []
COVERED: set[str] = set()


def check(name: str, passed: bool, detail: str = "", covers: str | None = None) -> None:
    print(f"{'PASS' if passed else 'FAIL'} {name} -- {detail}")
    if not passed:
        FAILURES.append(name)
    if passed and covers is not None:
        COVERED.add(covers)


def fstr(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atanh_interval_test(v: Fraction, terms: int = 13) -> tuple[Fraction, Fraction]:
    """Independent shorter positive-series interval."""
    low = sum((v ** (2 * n + 1) / (2 * n + 1) for n in range(terms)), Fraction(0))
    tail = v ** (2 * terms + 1) / ((2 * terms + 1) * (1 - v * v))
    return low, low + tail


def parity_mask(vertices: int, edges: list[tuple[int, int]], subset_mask: int) -> int:
    out = 0
    for edge_index, (u, w) in enumerate(edges):
        if (subset_mask >> edge_index) & 1:
            out ^= 1 << u
            out ^= 1 << w
    return out


def parity_sum(edges: list[tuple[int, int]], source_mask: int, v: Fraction) -> Fraction:
    total = Fraction(0)
    for subset in range(1 << len(edges)):
        if parity_mask(4, edges, subset) == source_mask:
            total += v ** subset.bit_count()
    return total


def support_connects(edges: list[tuple[int, int]], support: int, source: int, target: int) -> bool:
    adjacency = [[] for _ in range(4)]
    for index, (u, w) in enumerate(edges):
        if (support >> index) & 1:
            adjacency[u].append(w)
            adjacency[w].append(u)
    queue = deque([source])
    seen = {source}
    while queue:
        u = queue.popleft()
        for w in adjacency[u]:
            if w not in seen:
                seen.add(w)
                queue.append(w)
    return target in seen


def current_dp(v: Fraction) -> tuple[Fraction, Fraction, int]:
    """Edge-by-edge DP keyed by both source masks and union support."""
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    states = [
        (0, 0, 0, 1 - v * v),
        (0, 0, 1, v * v),
        (1, 1, 1, v * v),
        (1, 0, 1, v),
        (0, 1, 1, v),
    ]
    dp: dict[tuple[int, int, int], Fraction] = {(0, 0, 0): Fraction(1)}
    for edge_index, (u, w) in enumerate(edges):
        toggle = (1 << u) | (1 << w)
        nxt: defaultdict[tuple[int, int, int], Fraction] = defaultdict(Fraction)
        for (boundary1, boundary2, support), weight0 in dp.items():
            for parity1, parity2, present, local_weight in states:
                key = (
                    boundary1 ^ (toggle if parity1 else 0),
                    boundary2 ^ (toggle if parity2 else 0),
                    support | ((1 << edge_index) if present else 0),
                )
                nxt[key] += weight0 * local_weight
        dp = dict(nxt)
    total = Fraction(0)
    event = Fraction(0)
    supports = 0
    for (boundary1, boundary2, support), weight in dp.items():
        if boundary1 == boundary2 == 0:
            supports += 1
            total += weight
            if support_connects(edges, support, 0, 2):
                event += weight
    return total, event, supports


def main() -> int:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    data = artifact["data"]
    paired = data["paired_momentum"]
    current = data["random_current"]

    check("artifact_envelope", set(artifact) == {"meta", "data"} and isinstance(data.get("checks"), list), "meta/data with data.checks")
    check(
        "stored_checks_all_pass",
        bool(data["checks"]) and all(set(row) == {"name", "passed", "detail"} and row["passed"] for row in data["checks"]),
        f"{len(data['checks'])} stored rows",
    )

    input_rows = artifact["meta"]["provenance"]["inputs"]
    hash_ok = True
    for relative, row in input_rows.items():
        path = ROOT / relative
        hash_ok &= path.is_file() and sha256(path) == row["sha256"] and path.stat().st_size == row["size_bytes"]
    check("provenance_hashes", hash_ok, f"{len(input_rows)} exact inputs")

    coverage_rows = data["coverage_rows"]
    check(
        "coverage_manifest_exact",
        len(coverage_rows) == len(set(coverage_rows)) and set(coverage_rows) == EXPECTED_COVERAGE,
        ",".join(coverage_rows),
    )

    upper = json.loads(UPPER_INFRARED.read_text(encoding="utf-8"))
    i3_lo, i3_hi = map(Fraction, upper["data"]["certified_constants"]["I3"])
    k_test_lo, k_test_hi = atanh_interval_test(V)
    # Different coarse proof: for n>=1, 1/(2n+1)<=1/3.
    coarse_k_hi = V + V**3 / (3 * (1 - V * V))
    stored_k_lo, stored_k_hi = map(Fraction, paired["target"]["K_interval"])
    check(
        "target_below_incumbent_independent",
        k_test_lo <= stored_k_hi and stored_k_lo <= k_test_hi and stored_k_hi <= coarse_k_hi and 2 * coarse_k_hi < i3_lo,
        f"2*coarse K_hi={fstr(2*coarse_k_hi)} < I3_lo",
        "TARGET_BELOW_INCUMBENT",
    )
    stored_gap_lo, stored_gap_hi = map(Fraction, paired["target"]["incumbent_gap_interval"])
    check(
        "target_gap_rederived",
        stored_gap_lo == i3_lo / 2 - stored_k_hi
        and stored_gap_hi == i3_hi / 2 - stored_k_lo
        and stored_gap_lo > 0,
        f"gap_lo={fstr(stored_gap_lo)}",
    )

    # Independent bounded-spin/Parseval lift of the paired SDP witness.
    sites = list(itertools.product(range(2), repeat=3))
    spins = {site: (-1) ** sum(site) for site in sites}
    amplitudes = {}
    for mode in sites:
        amplitudes[mode] = sum(spins[site] * ((-1) ** sum(a * b for a, b in zip(mode, site))) for site in sites)
    n = len(sites)
    nonzero_power = sum(value * value for mode, value in amplitudes.items() if mode != (0, 0, 0))
    r_value = Fraction(nonzero_power, n * n)
    check(
        "paired_parseval_witness",
        amplitudes[(0, 0, 0)] == 0 and sum(value * value for value in amplitudes.values()) == n * n and r_value == 1,
        f"R={fstr(r_value)}, Fourier power={nonzero_power}",
        "PM_PARSEVAL_BOUNDED",
    )

    s = Fraction(paired["primal_sdp"]["primal_witness"]["s"])
    q = Fraction(paired["primal_sdp"]["primal_witness"]["q"])
    determinant = q - s * s
    feasible = q >= 0 and s - q >= 0 and 1 - s >= 0 and determinant >= 0
    exact_optimum = feasible and s == 1 and Fraction(paired["primal_sdp"]["exact_optimum_s"]) == 1
    dual_point = paired["dual_feasibility_sdp"]["exact_optimal_point"]
    dual_eta = Fraction(dual_point["eta"])
    dual_a = Fraction(dual_point["a"])
    dual_b = Fraction(dual_point["b"])
    dual_c = Fraction(dual_point["c"])
    dual_z00 = Fraction(dual_point["Z"][0][0])
    dual_z01 = Fraction(dual_point["Z"][0][1])
    dual_z11 = Fraction(dual_point["Z"][1][1])
    dual_feasible = (
        dual_a >= 0
        and dual_b >= 0
        and dual_c >= 0
        and dual_z00 >= 0
        and dual_z11 >= 0
        and dual_z00 * dual_z11 - dual_z01 * dual_z01 >= 0
        and dual_c + dual_z00 == 1 - dual_eta
        and dual_b - dual_c + 2 * dual_z01 == -1
        and dual_a - dual_b + dual_z11 == 0
    )
    whole_pm = paired["whole_incumbent_interval_obstruction"]
    check(
        "paired_sdp_exact_optimum",
        exact_optimum
        and dual_feasible
        and dual_eta == 0
        and Fraction(paired["primal_sdp"]["exact_optimum_order_floor_eta"]) == 0
        and whole_pm["exact_optimum_s"] == "1/1"
        and whole_pm["exact_optimum_order_floor_eta"] == "0/1",
        f"primal s={s}, det={determinant}; rational dual eta={dual_eta}; uniform for K<=I3/2",
        "PM_SDP_OPTIMUM",
    )
    eta = Fraction(paired["certificate_feasibility"]["challenge_eta"])
    challenge_upper = Fraction(paired["certificate_feasibility"]["challenge_upper_s"])
    challenge_gap = Fraction(paired["certificate_feasibility"]["exact_infeasibility_gap"])
    check(
        "paired_positive_floor_infeasible",
        eta == Fraction(1, 1024) and challenge_upper == 1 - eta and s - challenge_upper == challenge_gap == eta,
        f"exact separating witness gap={fstr(challenge_gap)}",
        "PM_POSITIVE_FLOOR_INFEASIBLE",
    )

    # Re-derive the two-current local weights directly from hyperbolic parity classes.
    weights = {
        "00_zero": 1 - V * V,
        "00_positive": V * V,
        "11_positive": V * V,
        "10_positive": V,
        "01_positive": V,
    }
    conditional = {
        "00": weights["00_positive"] / (weights["00_zero"] + weights["00_positive"]),
        "01": Fraction(1),
        "10": Fraction(1),
        "11": Fraction(1),
    }
    stored_states = {row["name"]: Fraction(row["weight"]) for row in current["local_state_derivation"]["states"]}
    check(
        "current_edge_weights_independent",
        stored_states == weights and conditional == {"00": Fraction(36, 625), "01": Fraction(1), "10": Fraction(1), "11": Fraction(1)},
        f"worst conditional={fstr(conditional['00'])}",
        "RC_EDGE_WEIGHTS",
    )

    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    p_empty = parity_sum(edges, 0, V)
    p_02 = parity_sum(edges, (1 << 0) | (1 << 2), V)
    dp_total, dp_event, support_rows = current_dp(V)
    stored_switch = current["finite_switching_evaluation"]
    switch_ok = (
        p_empty == 1 + V**4
        and p_02 == 2 * V**2
        and dp_total == p_empty**2
        and dp_event == p_02**2
        and Fraction(stored_switch["total_double_current_weight"]) == dp_total
        and Fraction(stored_switch["connection_0_2_weight"]) == dp_event
        and Fraction(stored_switch["connection_probability"]) == (p_02 / p_empty) ** 2
    )
    check(
        "switching_identity_c4_independent_dp",
        switch_ok,
        f"P0={fstr(p_empty)}, P02={fstr(p_02)}, support rows={support_rows}",
        "RC_SWITCHING_C4",
    )

    p_opt = min(conditional.values())
    stored_lp = current["domination_lp"]
    lp_gap = P0 - p_opt
    uniform_current = current["whole_incumbent_interval_obstruction"]
    uniform_gap = P0 - Fraction(1, 16)
    check(
        "current_domination_lp_exact",
        p_opt == Fraction(36, 625)
        and Fraction(stored_lp["exact_optimum"]) == p_opt
        and Fraction(stored_lp["exact_infeasibility_gap"]) == lp_gap == Fraction(2909, 3750)
        and i3_hi / 2 < Fraction(49, 192)
        and uniform_gap == Fraction(37, 48)
        and uniform_current["atanh_one_quarter_lower"] == "49/192"
        and uniform_current["uniform_gap_lower_to_5_over_6"] == "37/48",
        f"target p*={fstr(p_opt)}, gap={fstr(lp_gap)}; all K<=I3/2 have gap>{fstr(uniform_gap)}",
        "RC_DOMINATION_LP",
    )

    q_closed = 1 - P0
    ratio = 3 * q_closed
    tail = ratio**4 * (4 - 3 * ratio) / (1 - ratio) ** 2
    circuit_bound = Fraction(4, 3) * tail
    theta = 1 - circuit_bound
    stored_perc = current["self_contained_percolation_implication"]
    check(
        "planar_peierls_sum_exact",
        ratio == Fraction(1, 2)
        and circuit_bound == Fraction(5, 6)
        and theta == Fraction(1, 6)
        and Fraction(stored_perc["dual_circuit_union_bound"]) == circuit_bound
        and Fraction(stored_perc["percolation_probability_lower"]) == theta,
        f"circuit<={fstr(circuit_bound)}, theta>={fstr(theta)}",
        "PERC_PLANAR_5_6",
    )

    check(
        "all_universal_rows_reverified",
        COVERED == EXPECTED_COVERAGE,
        f"covered={sorted(COVERED)}",
    )
    outcome = data["outcome"]
    check(
        "honest_negative_outcome",
        outcome["improved_endpoint"] is False
        and outcome["challenged_endpoint_strictly_below_incumbent"] is True
        and "does not exclude" in outcome["scope_limit"],
        outcome["scope_limit"],
    )

    if FAILURES:
        print(f"FAIL ({len(FAILURES)}): {', '.join(FAILURES)}")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
