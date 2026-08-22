"""Clean-room exact checks for the v^28 simple-cubic HT artifact.

This verifier deliberately does not import ``experiments/e94_ht_v28.py``.  It
rebuilds the finite-lattice class bound, hyperrectangle cut-profile orbit
multiplicities, CRT bounds, and lower-order spin-transfer controls from stable
``ising`` APIs and the raw JSON certificate.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb, factorial, gcd, prod
from pathlib import Path

from ising.series import broken_bond_to_even_subgraph, log_series
from ising.transfer_matrix.crt import hybrid_box_broken_bond_poly

ROOT = Path(__file__).resolve().parents[1]
ORDER = 28
EXTERNAL_A28 = Fraction(525_549_581_866_326, 7)


def _check(name: str, passed: bool, detail: str) -> None:
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}")
    if not passed:
        raise AssertionError(f"{name}: {detail}")


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _canonical(shape: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(sorted(int(side) for side in shape))


def _is_prime_32(value: int) -> bool:
    """Local deterministic Miller--Rabin verifier for the artifact's 31-bit primes."""

    if value < 2:
        return False
    for prime in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if value == prime:
            return True
        if value % prime == 0:
            return False
    odd_part = value - 1
    powers_of_two = 0
    while odd_part % 2 == 0:
        powers_of_two += 1
        odd_part //= 2
    for base in (2, 3, 5, 7, 11, 13, 17):
        witness = pow(base, odd_part, value)
        if witness in (1, value - 1):
            continue
        for _ in range(powers_of_two - 1):
            witness = witness * witness % value
            if witness == value - 1:
                break
        else:
            return False
    return True


def _spin_transfer_cube_weights(side: int, order: int) -> tuple[Fraction, ...]:
    """Clean-room FLM inversion for one cube using only public transfer/series APIs."""

    target = (int(side),) * 3
    shapes = tuple(
        sorted(
            product(range(1, side + 1), repeat=3),
            key=lambda shape: (sum(shape), shape),
        )
    )
    logs: dict[tuple[int, int, int], tuple[Fraction, ...]] = {}
    weights: dict[tuple[int, int, int], tuple[Fraction, ...]] = {}
    for shape in shapes:
        canonical = _canonical(shape)
        if canonical not in logs:
            polynomial = broken_bond_to_even_subgraph(
                hybrid_box_broken_bond_poly(canonical), canonical, order
            )
            logs[canonical] = log_series(polynomial, order)
        exact_weight = list(logs[canonical])
        for subshape, subweight in weights.items():
            if all(sub <= bound for sub, bound in zip(subshape, shape)):
                placements = prod(
                    bound - sub + 1 for sub, bound in zip(subshape, shape)
                )
                for degree in range(order + 1):
                    exact_weight[degree] -= placements * subweight[degree]
        weights[shape] = tuple(exact_weight)
    return weights[target]


def _cut_inventory(shape: tuple[int, int, int]) -> tuple[tuple[tuple[int, int], ...], tuple[int, ...]]:
    sides = tuple(int(side) for side in shape)
    cuts = tuple(
        (axis, position)
        for axis, side in enumerate(sides)
        for position in range(side - 1)
    )
    edge_counts = tuple(
        prod(sides[other] for other in range(3) if other != axis)
        for axis, _ in cuts
    )
    return cuts, edge_counts


def _profile_orbits(shape: tuple[int, int, int], extras: int) -> list[dict]:
    """Return exact reflection/permutation orbits of positive half profiles."""

    sides = tuple(int(side) for side in shape)
    cuts, _ = _cut_inventory(sides)
    lookup = {cut: index for index, cut in enumerate(cuts)}
    actions = []
    for axis_permutation in permutations(range(3)):
        if any(sides[axis] != sides[axis_permutation[axis]] for axis in range(3)):
            continue
        for reflected in product((False, True), repeat=3):
            actions.append(
                tuple(
                    lookup[
                        (
                            axis_permutation[axis],
                            sides[axis] - 2 - position if reflected[axis] else position,
                        )
                    ]
                    for axis, position in cuts
                )
            )

    profiles: set[tuple[int, ...]] = set()
    if extras == 0:
        profiles.add((1,) * len(cuts))
    elif extras == 1:
        for index in range(len(cuts)):
            profile = [1] * len(cuts)
            profile[index] = 2
            profiles.add(tuple(profile))
    elif extras == 2:
        for index in range(len(cuts)):
            profile = [1] * len(cuts)
            profile[index] = 3
            profiles.add(tuple(profile))
        for left, right in combinations(range(len(cuts)), 2):
            profile = [1] * len(cuts)
            profile[left] = profile[right] = 2
            profiles.add(tuple(profile))
    else:
        raise ValueError("test only needs order-28 excess profiles")

    output = []
    while profiles:
        representative = min(profiles)
        orbit = {
            tuple(profile[action[index]] for index in range(len(cuts)))
            for action in actions
            for profile in (representative,)
        }
        output.append(
            {
                "profile": list(representative),
                "multiplicity": len(orbit),
                "orbit": sorted(list(item) for item in orbit),
            }
        )
        profiles -= orbit
    return sorted(output, key=lambda row: (row["profile"], row["multiplicity"]))


def _ordered_partition_factor(slot_count: int) -> int:
    stirling = [[0] * (slot_count + 1) for _ in range(slot_count + 1)]
    stirling[0][0] = 1
    for row in range(1, slot_count + 1):
        for blocks in range(1, row + 1):
            stirling[row][blocks] = (
                stirling[row - 1][blocks - 1]
                + blocks * stirling[row - 1][blocks]
            )
    return sum(
        stirling[slot_count][blocks] * factorial(blocks - 1)
        for blocks in range(1, slot_count + 1)
    )


def _profile_bound(shape: tuple[int, int, int], profile: tuple[int, ...]) -> int:
    _, edge_counts = _cut_inventory(shape)
    slots = sum(profile)
    return _ordered_partition_factor(slots) * prod(
        comb(edge_count, 2 * exponent) ** exponent
        for edge_count, exponent in zip(edge_counts, profile)
    )


def _artifact_profile_rows(data: dict, shape: tuple[int, int, int], degree: int) -> list[dict]:
    wanted_shape = list(shape)
    return sorted(
        [
            row
            for row in data["profiles"]
            if row["shape"] == wanted_shape and row["degree"] == degree
        ],
        key=lambda row: (row["profile"], row["multiplicity"]),
    )


def _telemetry_ok(row: dict) -> bool:
    timeline = row["timeline"]
    return (
        bool(timeline)
        and row["peak_states"] == max(entry["states"] for entry in timeline)
        and row["state_count_sha256"]
        == hashlib.sha256(
            json.dumps(
                [entry["states"] for entry in timeline],
                sort_keys=True,
                separators=(",", ":"),
            ).encode("ascii")
        ).hexdigest()
        and int(row["max_rss_bytes"]) > 0
        and float(row["wall_seconds"]) > 0
    )


def main() -> None:
    payload = _load(ROOT / "results" / "series" / "ht_v28.json")
    data = payload["data"]
    series = data["series"]
    interaction = tuple(Fraction(value) for value in series["interaction_coefficients"])
    coefficient = tuple(Fraction(value) for value in series["coefficients"])

    extent = ORDER // 2
    canonical = sorted(
        {
            _canonical(shape)
            for shape in product(range(1, extent + 2), repeat=3)
            if sum(side - 1 for side in shape) <= extent
        }
    )
    minimal = [shape for shape in canonical if sum(shape) == 17]
    _check(
        "finite-lattice order-28 class bound",
        data["box_class_lemma"]["contributing_classes"] == "a+b+c <= 17"
        and data["box_class_lemma"]["contributing_canonical_classes"]
        == [list(shape) for shape in canonical]
        and data["box_class_lemma"]["minimal_span_classes"]
        == [list(shape) for shape in minimal]
        and len(canonical) == 147
        and len(minimal) == 24,
        "2*sum(side-1)<=28 gives exactly 147 canonical classes and 24 span-28 classes",
    )

    expected_walled = [(4, 6, 6), (4, 6, 7), (5, 5, 5), (5, 5, 6), (5, 5, 7), (5, 6, 6)]
    _check(
        "walled classes are complete",
        data["box_class_lemma"]["walled"] == [list(shape) for shape in expected_walled]
        and [shape for shape in canonical if shape[0] * shape[1] > 22]
        == expected_walled,
        "all and only canonical classes with open transfer section above 22 are injected",
    )

    # This is an independent public-API spin-transfer FLM inversion, not a
    # producer import or a monkey patch.  The controls exercise cap-2, ordinary
    # cap-4, and the new radix-4 cap-6/pair profiles respectively.
    spin333 = _spin_transfer_cube_weights(3, 16)
    spin444 = _spin_transfer_cube_weights(4, 20)
    controls = data["controls"]["independent_spin_transfer"]
    _check(
        "independent generalized-profile spin-transfer controls",
        Fraction(controls["cube333_degree16"]["aggregate"]) == spin333[16]
        and Fraction(controls["cube444_degree18"]["aggregate"]) == spin444[18]
        and Fraction(controls["cube444_degree20"]["aggregate"]) == spin444[20],
        "public spin-transfer FLM reproduces generic cap-6/pair, cap-2, and cap-4 aggregates",
    )

    expected_555 = _profile_orbits((5, 5, 5), 2)
    observed_555 = _artifact_profile_rows(data, (5, 5, 5), 28)
    _check(
        "cube degree-28 profile orbit certificate",
        [
            {"profile": row["profile"], "multiplicity": row["multiplicity"]}
            for row in observed_555
        ]
        == [
            {"profile": row["profile"], "multiplicity": row["multiplicity"]}
            for row in expected_555
        ]
        and len(observed_555) == 9
        and sum(row["multiplicity"] for row in observed_555) == 78,
        "two cap-6 and seven pair-cap-4 orbits cover all 12 + C(12,2) degree-28 profiles",
    )

    expected_556 = _profile_orbits((5, 5, 6), 1)
    expected_466 = _profile_orbits((4, 6, 6), 1)
    _check(
        "rectangular one-cap-four orbit certificates",
        [
            {"profile": row["profile"], "multiplicity": row["multiplicity"]}
            for row in _artifact_profile_rows(data, (5, 5, 6), 28)
        ]
        == [
            {"profile": row["profile"], "multiplicity": row["multiplicity"]}
            for row in expected_556
        ]
        and [
            {"profile": row["profile"], "multiplicity": row["multiplicity"]}
            for row in _artifact_profile_rows(data, (4, 6, 6), 28)
        ]
        == [
            {"profile": row["profile"], "multiplicity": row["multiplicity"]}
            for row in expected_466
        ],
        "all span-13 degree-28 profiles are covered exactly once up to box symmetry",
    )

    family_fields = {
        ((5, 5, 5), 24): "w555_24",
        ((5, 5, 5), 26): "w555_26",
        ((5, 5, 5), 28): "w555_28",
        ((5, 5, 6), 26): "w556_26",
        ((5, 5, 6), 28): "w556_28",
        ((4, 6, 6), 26): "w466_26",
        ((4, 6, 6), 28): "w466_28",
        ((4, 6, 7), 28): "w467_28",
        ((5, 5, 7), 28): "w557_28",
        ((5, 6, 6), 28): "w566_28",
    }
    observed_families = {
        (tuple(row["shape"]), row["degree"])
        for row in data["profiles"]
    }
    weights = data["weights"]
    expected_injections: dict[str, dict[str, str]] = {}
    for (shape, degree), field in family_fields.items():
        expected_injections.setdefault(str(list(shape)), {})[str(degree)] = weights[field]
    _check(
        "complete profile families reconstruct every injected weight",
        observed_families == set(family_fields)
        and all(
            [
                {"profile": row["profile"], "multiplicity": row["multiplicity"]}
                for row in _artifact_profile_rows(data, shape, degree)
            ]
            == [
                {"profile": row["profile"], "multiplicity": row["multiplicity"]}
                for row in _profile_orbits(
                    shape, degree // 2 - (sum(shape) - 3)
                )
            ]
            and sum(
                row["multiplicity"] * int(row["centered"])
                for row in _artifact_profile_rows(data, shape, degree)
            )
            == int(weights[field])
            for (shape, degree), field in family_fields.items()
        )
        and weights["injections"] == expected_injections,
        "all ten lower/order-28 profile families, orbit sums, and FLM injections agree",
    )

    crt = data["crt"]
    modulus = int(crt["modulus_product"])
    rows = data["profiles"]
    _check(
        "profile-wise CRT centered uniqueness",
        prod(crt["primes"]) == modulus
        and len(set(crt["primes"])) == len(crt["primes"])
        and all(_is_prime_32(prime) for prime in crt["primes"])
        and all(
            gcd(crt["primes"][left], crt["primes"][right]) == 1
            for left in range(len(crt["primes"]))
            for right in range(left)
        )
        and all(
            gcd(*tuple(row["profile"])) == 1
            and int(row["bound"])
            == _profile_bound(tuple(row["shape"]), tuple(row["profile"]))
            and modulus > 2 * int(row["bound"])
            and int(row["residue"]) == int(row["centered"]) % modulus
            and abs(int(row["centered"])) <= int(row["bound"])
            for row in rows
        ),
        f"M={modulus} exceeds twice every independently rebuilt profile bound",
    )

    v24 = _load(ROOT / "results" / "series" / "ht_v24.json")
    v26 = _load(ROOT / "results" / "series" / "ht_v26.json")
    v24_coefficients = tuple(Fraction(value) for value in v24["data"]["series"]["coefficients"])
    v26_coefficients = tuple(Fraction(value) for value in v26["data"]["series"]["coefficients"])
    _check(
        "new pipeline reproduces v24 and v26 exactly",
        tuple(Fraction(value) for value in data["reproductions"]["v24"]["coefficients"])
        == v24_coefficients
        and tuple(Fraction(value) for value in data["reproductions"]["v26"]["coefficients"])
        == v26_coefficients
        and coefficient[: len(v24_coefficients)] == v24_coefficients
        and coefficient[: len(v26_coefficients)] == v26_coefficients,
        "fresh generalized-register lower-order assemblies match both stored exact series",
    )

    _check(
        "v28 external witness and internal normalization",
        interaction[28] == EXTERNAL_A28
        and coefficient[28] == interaction[28] + Fraction(3, 28)
        and coefficient[27] == 0
        and series["v28"] == str(coefficient[28])
        and series["interaction_v28"] == str(interaction[28]),
        f"a_28={interaction[28]} and [v^28] phi={coefficient[28]}",
    )

    _check(
        "frontier telemetry and producer envelope",
        all(_telemetry_ok(row) for row in data["profiles"])
        and data["production_resource_measurement"]["maximum_resident_set_bytes"] > 0
        and float(data["production_resource_measurement"]["real_seconds"]) > 0
        and payload["provenance"]["script"] == "experiments/e94_ht_v28.py"
        and payload["provenance"]["interpreter"] == ".venv/bin/python"
        and not payload["provenance"]["benchmark_kc_used"]
        and all(check["passed"] for check in payload["checks"]),
        "raw profile telemetry, provenance, and all producer checks are coherent",
    )
    print("PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL test_ht_v28: {error}")
        raise
