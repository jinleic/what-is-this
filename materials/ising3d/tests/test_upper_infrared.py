"""Standalone independent audit of the infrared upper-bound certificate.

Every quantity is recomputed by a route that differs from the experiment:
the walk Green function is checked through its exact discrete Laplace equation
and through the multinomial return formula, the torus Green function through the
same equation with the zero-mode defect, the barrier identities by exact
rational sampling instead of symbolic algebra, and the L=2 Ising correlations
by even-subgraph high-temperature enumeration at the rational v=tanh(K).
"""

from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction
from math import factorial
import itertools
import json
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "bounds" / "upper_infrared.json"
INCUMBENT = ROOT / "results" / "bounds" / "kc_bounds.json"
IMPROVED = ROOT / "results" / "bounds" / "improved_bounds.json"
DPS = 100
WALK_STEPS = 30
BENCHMARK = Decimal("0.221654626")  # comparison only
RATIONAL_COSINES = {
    2: (Fraction(1), Fraction(-1)),
    3: (Fraction(1), Fraction(-1, 2), Fraction(-1, 2)),
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


def _ends(value: object) -> tuple[Decimal, Decimal]:
    text = str(value).strip()
    assert text.startswith("[") and text.endswith("]"), text
    low, high = text[1:-1].split(",", maxsplit=1)
    return Decimal(low.strip()), Decimal(high.strip())


def _recompute_constants() -> dict[str, tuple[Decimal, Decimal]]:
    previous = mp.iv.dps
    try:
        mp.iv.dps = DPS
        q = mp.iv
        i3 = (
            q.sqrt(6)
            * q.gamma(q.mpf(1) / 24)
            * q.gamma(q.mpf(5) / 24)
            * q.gamma(q.mpf(7) / 24)
            * q.gamma(q.mpf(11) / 24)
            / (96 * q.pi**3)
        )
        watson = 3 * i3
        return {
            "I3": _ends(i3),
            "Wsc": _ends(watson),
            "I3_over_2": _ends(i3 / 2),
            "nn": _ends(1 - 1 / watson),
            "atanh_one_fifth": _ends(
                q.log((1 + q.mpf(1) / 5) / (1 - q.mpf(1) / 5)) / 2
            ),
        }
    finally:
        mp.iv.dps = previous


def _walk_probabilities(steps: int) -> list[dict[tuple[int, int, int], Fraction]]:
    layers = [{(0, 0, 0): Fraction(1)}]
    for _ in range(steps):
        following: dict[tuple[int, int, int], Fraction] = {}
        for point, weight in layers[-1].items():
            for axis in range(3):
                for step in (-1, 1):
                    moved = list(point)
                    moved[axis] += step
                    key = (moved[0], moved[1], moved[2])
                    following[key] = following.get(key, Fraction(0)) + weight / 6
        layers.append(following)
    return layers


def _weighted_returns(
    couplings: tuple[Fraction, Fraction, Fraction], steps: int
) -> tuple[Fraction, list[Fraction]]:
    """Exact anisotropic return probabilities of the sub-model walk."""

    total = sum(couplings, Fraction(0))
    step_weights = [coupling / (2 * total) for coupling in couplings]
    states: dict[tuple[int, int, int], Fraction] = {(0, 0, 0): Fraction(1)}
    returns = [Fraction(1)]
    for _ in range(steps):
        following: dict[tuple[int, int, int], Fraction] = {}
        for point, weight in states.items():
            for axis in range(3):
                for direction in (-1, 1):
                    moved = list(point)
                    moved[axis] += direction
                    key = (moved[0], moved[1], moved[2])
                    following[key] = (
                        following.get(key, Fraction(0)) + weight * step_weights[axis]
                    )
        states = following
        returns.append(states.get((0, 0, 0), Fraction(0)))
    return total, returns


def _multinomial_return(pairs: int) -> Fraction:
    """P_{2n}(0->0) from the exact simple-cubic multinomial formula."""

    total = 0
    for i in range(pairs + 1):
        for j in range(pairs - i + 1):
            k = pairs - i - j
            total += factorial(2 * pairs) // (
                factorial(i) ** 2 * factorial(j) ** 2 * factorial(k) ** 2
            )
    return Fraction(total, 6 ** (2 * pairs))


def _torus_green_by_equation(side: int) -> dict[tuple[int, int, int], Fraction]:
    cosines = RATIONAL_COSINES[side]
    sites = list(itertools.product(range(side), repeat=3))
    modes = [mode for mode in sites if mode != (0, 0, 0)]
    volume = side**3
    values: dict[tuple[int, int, int], Fraction] = {}
    for site in sites:
        total = Fraction(0)
        for mode in modes:
            dispersion = 3 - sum(cosines[axis_index] for axis_index in mode)
            phase = sum(mode[axis] * site[axis] for axis in range(3)) % side
            total += cosines[phase] / dispersion
        values[site] = total / volume
    return values


def _cube_correlations_high_temperature(v: Fraction) -> dict[tuple[int, int, int], Fraction]:
    """Even-subgraph high-temperature correlations of the single-bond L=2 torus."""

    sites = list(itertools.product(range(2), repeat=3))
    index = {site: position for position, site in enumerate(sites)}
    edges = []
    for site in sites:
        for axis in range(3):
            if site[axis] == 0:
                other = list(site)
                other[axis] = 1
                edges.append((index[site], index[(other[0], other[1], other[2])]))
    even = Fraction(0)
    odd = {site: Fraction(0) for site in sites}
    for subset in range(1 << len(edges)):
        parity = 0
        size = 0
        for position, (left, right) in enumerate(edges):
            if subset >> position & 1:
                size += 1
                parity ^= (1 << left) | (1 << right)
        weight = v**size
        if parity == 0:
            even += weight
        else:
            bits = [bit for bit in range(len(sites)) if parity >> bit & 1]
            if len(bits) == 2 and bits[0] == 0:
                odd[sites[bits[1]]] += weight
    result = {site: odd[site] / even for site in sites}
    result[(0, 0, 0)] = Fraction(1)
    return result


def main() -> None:
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    incumbent = json.loads(INCUMBENT.read_text(encoding="utf-8"))
    improved = json.loads(IMPROVED.read_text(encoding="utf-8"))
    assert set(payload) == {"provenance", "data", "checks"}
    assert payload["provenance"]["script"] == "experiments/e72_upper_infrared.py"
    assert payload["provenance"]["precision"]["mpmath_iv_dps"] == DPS
    assert payload["provenance"]["precision"]["rounding"] == "outward"
    assert payload["checks"] and all(check["passed"] for check in payload["checks"])
    data = payload["data"]

    # ---- certified constants ------------------------------------------------
    constants = _recompute_constants()
    for key, recomputed in (
        ("I3", constants["I3"]),
        ("Wsc", constants["Wsc"]),
        ("I3_over_2", constants["I3_over_2"]),
        ("atanh_one_fifth", constants["atanh_one_fifth"]),
        ("witness_nearest_neighbour_1_minus_1_over_Wsc", constants["nn"]),
    ):
        stored = tuple(map(Decimal, data["certified_constants"][key]))
        assert stored[0] <= recomputed[0] <= recomputed[1] <= stored[1], key
        assert stored[1] - stored[0] < Decimal("1e-90"), key
    retained = Decimal(data["certified_conclusion"]["retained_upper_bound"])
    stored_upper = tuple(map(Decimal, data["certified_constants"]["I3_over_2"]))
    with localcontext() as context:
        context.prec = 140
        assert stored_upper[0] <= constants["I3"][1] / 2
        assert constants["I3"][0] / 2 <= stored_upper[1]
    assert retained >= stored_upper[1]
    assert retained - stored_upper[1] < Decimal("1e-40")
    old = incumbent["data"]["upper_bound_infrared"]
    assert retained == Decimal(old["certified_decimal_upper"])
    assert Decimal(old["kc_upper_interval"][0]) <= stored_upper[0]
    assert stored_upper[1] <= Decimal(old["kc_upper_interval"][1])
    print("certified constants and retained endpoint: PASS")

    # ---- honest no-improvement conclusion ----------------------------------
    assert data["certified_conclusion"]["strict_improvement"] is False
    assert data["certified_conclusion"]["new_upper_bound"] is None
    certified_lower = Decimal(
        improved["data"]["final_certified_interval"]["decimal_outward_rounded"][0]
    )
    assert certified_lower < BENCHMARK < retained
    assert Decimal(data["barrier_window"]["witness_survives_from"]) == certified_lower
    assert Decimal(data["barrier_window"]["witness_survives_to"]) == retained
    exclusion = tuple(
        map(Decimal, data["barrier_window"]["elementary_susceptibility_exclusion_below"])
    )
    assert exclusion[1] < certified_lower < retained
    print("no-improvement headline and barrier window ordering: PASS")

    # ---- walk Green function: exact Laplace equation with truncation defect --
    probabilities = _walk_probabilities(WALK_STEPS + 1)
    truncated = {}
    support = set()
    for layer in probabilities[: WALK_STEPS + 1]:
        support |= layer.keys()
    for site in support:
        truncated[site] = (
            sum(layer.get(site, Fraction(0)) for layer in probabilities[: WALK_STEPS + 1])
            / 3
        )

    def value(site: tuple[int, int, int]) -> Fraction:
        return truncated.get(site, Fraction(0))

    checked = 0
    for site in ((0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 0, 0), (2, 1, 0)):
        neighbours = Fraction(0)
        for axis in range(3):
            for step in (-1, 1):
                moved = list(site)
                moved[axis] += step
                neighbours += value((moved[0], moved[1], moved[2]))
        residual = 3 * value(site) - neighbours / 2
        expected = (Fraction(1) if site == (0, 0, 0) else Fraction(0)) - probabilities[
            WALK_STEPS + 1
        ].get(site, Fraction(0))
        assert residual == expected, site
        checked += 1
    assert checked == 6
    stored_bounds = data["barrier_witness"]["walk_lower_bounds_C_z"]
    for key, text in stored_bounds.items():
        site = tuple(json.loads(key))
        stored = Fraction(text)
        assert stored > 0, key
        assert stored == value(site), key
    print("witness real-space positivity and Laplace equation: PASS")

    # ---- the visit series is Wsc, checked by an independent formula ----------
    stored_sums = {
        int(steps): Fraction(text)
        for steps, text in data["barrier_witness"]["visit_partial_sums_of_Wsc"].items()
    }
    independent = sum(
        (_multinomial_return(pairs) for pairs in range(WALK_STEPS // 2 + 1)),
        Fraction(0),
    )
    assert stored_sums[WALK_STEPS] == independent
    assert stored_sums[10] < stored_sums[20] < stored_sums[WALK_STEPS]
    assert Decimal(independent.numerator) / Decimal(independent.denominator) < constants[
        "Wsc"
    ][0]
    nn_stored = tuple(map(Decimal, data["barrier_witness"]["nearest_neighbour_interval"]))
    with localcontext() as context:
        context.prec = 140
        # 1-1/Wsc is exactly the infrared saturation value 1-1/(6K) at K=I3/2.
        saturation = (
            Decimal(1) - Decimal(1) / (6 * (constants["I3"][1] / 2)),
            Decimal(1) - Decimal(1) / (6 * (constants["I3"][0] / 2)),
        )
        assert nn_stored[0] <= saturation[1] and saturation[0] <= nn_stored[1]
    assert Decimal("0.34") < nn_stored[0] <= nn_stored[1] < Decimal("0.341")
    print("witness visit series, Wsc identity and saturation value: PASS")

    # ---- barrier identities by exact rational sampling ----------------------
    i3_low = constants["I3"][0]
    i3_rational = Fraction(int(i3_low.scaleb(60)), 10**60)  # certified lower endpoint
    assert Fraction(1, 2) < i3_rational < Fraction(51, 100)
    # The excess coefficients 12*I3-4 and 15*I3-6 increase in I3 and
    # K+1-3*I3 decreases in I3, so verifying the signs at the certified lower
    # endpoint proves them for the true I3.
    assert 12 * i3_rational - 4 > 0
    assert 15 * i3_rational - 6 > 0
    assert i3_rational / 2 + 1 - 3 * i3_rational < 0
    for t in (
        Fraction(0),
        Fraction(1, 1000),
        Fraction(1, 7),
        Fraction(1),
        Fraction(9, 2),
        Fraction(101),
    ):
        infrared = (
            i3_rational
            + 2 * t * (3 * i3_rational - 1)
            + t**2 * (9 * i3_rational - 3)
        )
        threshold = infrared / (2 * (1 + Fraction(3, 2) * t**2))
        excess = threshold - i3_rational / 2
        claimed = (
            t
            * ((12 * i3_rational - 4) + (15 * i3_rational - 6) * t)
            / (2 * (2 + 3 * t**2))
        )
        assert excess == claimed
        assert (excess == 0) == (t == 0)
        assert excess >= 0
        for coupling in (
            Fraction(1, 6),
            Fraction(1, 5),
            i3_rational / 2,
            Fraction(1, 4),
        ):
            augmented = 2 * coupling * (
                1 + Fraction(3, 2) * t**2 + 6 * t * (1 - Fraction(1, 6) / coupling)
            ) - infrared
            factored = (2 * coupling - i3_rational) * (1 + 6 * t) + 3 * t**2 * (
                coupling + 1 - 3 * i3_rational
            )
            assert augmented == factored
            if coupling <= i3_rational / 2:
                assert augmented <= 0
    assert data["test_function_audit"]["minimizer"] == "t=0"
    assert data["test_function_audit"]["minimum"] == "I3/2"
    for candidate in data["test_function_audit"]["rejected_pointwise_candidates"].values():
        assert candidate["status"] == "not a bound"
        low, high = map(Decimal, candidate["number_if_block_variance_is_dropped"])
        assert high < stored_upper[0]
        assert low < high
    print("exact barrier identities and rejected candidates: PASS")

    # ---- anisotropic sub-model barrier -------------------------------------
    anisotropic = data["anisotropic_audit"]
    steps = anisotropic["termwise_verification_steps"]
    isotropic_total, isotropic_returns = _weighted_returns(
        (Fraction(1), Fraction(1), Fraction(1)), steps
    )
    assert anisotropic["rows"]
    for row in anisotropic["rows"]:
        couplings = tuple(Fraction(text) for text in row["couplings"])
        assert all(0 < coupling <= 1 for coupling in couplings)
        total, returns = _weighted_returns(couplings, steps)
        assert Fraction(row["coupling_sum"]) == total <= 3
        assert all(
            returns[n] / total >= isotropic_returns[n] / isotropic_total
            for n in range(steps + 1)
        )
        assert row["termwise_dominates_isotropic"] is True
        assert Fraction(row["partial_sum_T_alpha"]) == sum(returns, Fraction(0)) / total
        assert Fraction(row["partial_sum_T_alpha"]) > Fraction(
            row["partial_sum_T_isotropic"]
        )
    print("anisotropic sub-model barrier: PASS")

    # ---- general nonnegative-kernel barrier --------------------------------
    # q^T C q >= I3 * ||q||^2 for nonnegative kernels; the excess for the
    # two-site kernel delta_0 + delta_e is exactly C(e) > 0.
    two_site_excess = value((1, 0, 0))
    assert two_site_excess > 0
    assert Decimal(two_site_excess.numerator) / Decimal(
        two_site_excess.denominator
    ) < constants["I3"][0]
    print("nonnegative-kernel excess certificate: PASS")

    # ---- exact finite-volume torus data ------------------------------------
    stored_rows = {row["side"]: row for row in data["finite_volume_audit"]["exact_torus_green_functions"]}
    assert set(stored_rows) == {2, 3, 4, 6}
    previous_origin = None
    for side in (2, 3, 4, 6):
        green = _torus_green_by_equation(side)
        row = stored_rows[side]
        assert Fraction(row["C_L_at_origin"]) == green[(0, 0, 0)]
        assert Fraction(row["minimum_value"]) == min(green.values())
        assert Fraction(row["sum_over_sites"]) == sum(green.values()) == 0
        assert Fraction(row["minimum_value"]) < 0
        assert Fraction(row["C_L_at_origin"]) < Fraction(
            int(constants["I3"][0].scaleb(60)), 10**60
        )
        if previous_origin is not None:
            assert previous_origin < Fraction(row["C_L_at_origin"])
        previous_origin = Fraction(row["C_L_at_origin"])
        # torus Laplace equation with the exact zero-mode defect
        for site in itertools.product(range(side), repeat=3):
            neighbours = Fraction(0)
            for axis in range(3):
                for step in (-1, 1):
                    moved = list(site)
                    moved[axis] = (moved[axis] + step) % side
                    neighbours += green[(moved[0], moved[1], moved[2])]
            residual = 3 * green[site] - neighbours / 2
            expected = (Fraction(1) if site == (0, 0, 0) else Fraction(0)) - Fraction(
                1, side**3
            )
            assert residual == expected, (side, site)
        half = Fraction(row["finite_volume_half_value_not_a_bound"])
        assert half == green[(0, 0, 0)] / 2
        assert Decimal(half.numerator) / Decimal(half.denominator) < stored_upper[0]
    print("exact torus Green functions, defect equation and rejected halves: PASS")

    # ---- L=2 energy identity and ceiling violation -------------------------
    l2 = data["finite_volume_audit"]["l2_energy_identity"]
    correlations = _cube_correlations_high_temperature(Fraction(1, 4))
    nearest = (
        correlations[(1, 0, 0)] + correlations[(0, 1, 0)] + correlations[(0, 0, 1)]
    ) / 3
    assert Fraction(l2["G_at_origin"]) == 1
    assert Fraction(l2["nearest_neighbour_correlation"]) == nearest
    assert Fraction(l2["weighted_dispersion_sum"]) == Fraction(
        l2["three_times_one_minus_nearest"]
    ) == 3 * (1 - nearest)
    ceiling_test = data["finite_volume_audit"]["single_bond_l2_torus_ceiling_test"]
    assert ceiling_test["exponentiated_coupling_w_equals_exp_minus_2K"] == "3/5"
    cosines = RATIONAL_COSINES[2]
    sites = list(itertools.product(range(2), repeat=3))
    assert len(ceiling_test["modes"]) == 7
    for row in ceiling_test["modes"]:
        mode = tuple(row["mode"])
        dispersion = 3 - sum(cosines[axis_index] for axis_index in mode)
        assert dispersion == row["dispersion"] and dispersion > 0
        transform = sum(
            correlations[site]
            * cosines[sum(mode[axis] * site[axis] for axis in range(3)) % 2]
            for site in sites
        )
        assert Fraction(row["Ghat_exact"]) == transform
        ceiling_high = Decimal(row["ceiling_interval"][1])
        assert (
            Decimal(transform.numerator) / Decimal(transform.denominator) > ceiling_high
        )
        assert row["ceiling_violated"] is True
    print("L=2 energy identity and ceiling-violation audit: PASS")

    # ---- method optimality wording is scoped, not overclaimed ---------------
    theorem = data["method_optimality_theorem"]
    assert "zero-magnetization solution" in theorem["statement"]
    assert "not a no-go theorem" in theorem["scope"]
    assert data["certified_conclusion"]["benchmark_role"].startswith("comparison only")
    print("scope statements present: PASS")
    print("PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"FAIL: {exc}")
        raise
