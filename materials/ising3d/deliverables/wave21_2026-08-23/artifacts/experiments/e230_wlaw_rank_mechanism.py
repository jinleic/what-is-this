#!/usr/bin/env python3
"""Exact two-slice observability and a bonding-truncation countercertificate.

The producer does not launch a ladder closure.  It uses the already certified
exact annihilator bases W_L = U_L^perp for L=3..8 and tests a direct evaluation
map: retain only configurations whose absolute leg imbalance is at most two.
The modular calculations below are lower bounds until they meet the omitted
coordinate count, which is an independent rational ceiling.

A second, independent calculation changes each rung to the unnormalised basis
0, b=top+bottom, a=top-bottom, d=top+bottom occupied.  It proves the local
D/F/D+ action is block-tridiagonal in the number of a factors and gives the
first exact kernel of the tempting truncation to at most two a factors.
"""
from __future__ import annotations

import hashlib
import json
import resource
import sys
import time
from fractions import Fraction
from itertools import product
from math import comb
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e230_wlaw_rank_mechanism.py"
OUT_PATH = ROOT / "results" / "algebra_growth" / "wlaw_rank_mechanism.json"
SOURCE_PATH = ROOT / "results" / "ladder" / "l8_saturation.json"
L_VALUES = tuple(range(3, 9))
PRIMES = (999_983, 1_000_003)
RSS_CAP_BYTES = 2 * 1024**3
LABELS = ("0", "b", "a", "d")
PARTICLE_WEIGHT = {"0": 0, "b": 1, "a": 1, "d": 2}


def peak_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return value == divisor
        divisor += 1 if divisor == 2 else 2
    return True


def tau_config(config: int, L: int) -> int:
    image = 0
    for rung in range(L):
        image |= ((config >> (2 * rung)) & 1) << (2 * rung + 1)
        image |= ((config >> (2 * rung + 1)) & 1) << (2 * rung)
    return image


def rho_config(config: int, L: int) -> int:
    image = 0
    for rung in range(L):
        image |= ((config >> (2 * rung)) & 3) << (2 * (L - 1 - rung))
    return image


def leg_imbalance(config: int, L: int) -> int:
    top = sum((config >> (2 * rung)) & 1 for rung in range(L))
    return abs(2 * top - config.bit_count())


def configuration_orbits(L: int) -> list[dict[str, int]]:
    seen: set[int] = set()
    records: list[dict[str, int]] = []
    for config in range(1 << (2 * L)):
        if config.bit_count() % 2 or config in seen:
            continue
        reflected = rho_config(config, L)
        orbit = {
            config,
            tau_config(config, L),
            reflected,
            tau_config(reflected, L),
        }
        seen.update(orbit)
        representative = min(orbit)
        records.append(
            {
                "representative": representative,
                "orbit_size": len(orbit),
                "particle_sector": representative.bit_count(),
                "absolute_leg_imbalance": leg_imbalance(representative, L),
            }
        )
    return records


def sparse_modular_rank(
    rows: Iterable[dict[int, int]], prime: int
) -> tuple[int, tuple[int, ...]]:
    """Exact sparse row rank over F_prime; the greatest column is the pivot."""
    pivots: dict[int, dict[int, int]] = {}
    for source in rows:
        vector = {
            int(column): int(coefficient) % prime
            for column, coefficient in source.items()
            if int(coefficient) % prime
        }
        while vector:
            pivot = max(vector)
            coefficient = vector[pivot]
            old = pivots.get(pivot)
            if old is None:
                inverse = pow(coefficient, prime - 2, prime)
                vector = {
                    column: (entry * inverse) % prime
                    for column, entry in vector.items()
                    if (entry * inverse) % prime
                }
                pivots[pivot] = vector
                break
            for column, entry in old.items():
                updated = (vector.get(column, 0) - coefficient * entry) % prime
                if updated:
                    vector[column] = updated
                else:
                    vector.pop(column, None)
    return len(pivots), tuple(sorted(pivots))


def source_records(source: dict[str, object]) -> dict[int, dict[str, object]]:
    data = source["data"]
    records = {
        int(record["L"]): record for record in data["regression_records"]
    }
    records[7] = data["L7_record"]
    records[8] = data["L8_record"]
    return records


def source_exact(record: dict[str, object]) -> bool:
    validation = record["validation"]
    return (
        len(record["basis"]) == int(record["W_dim"])
        and sum(int(value) for value in record["W_sectors"].values())
        == int(record["W_dim"])
        and int(record["cyclic_Q_dim"]) + int(record["W_dim"])
        == int(record["K_dim"])
        and bool(validation["A_invariant_over_Q"])
        and bool(validation["basis_independent_over_Q"])
        and bool(validation["four_B_invariant_over_Q"])
        and bool(validation["sector_homogeneous"])
        and bool(validation["vacuum_orthogonal"])
        and not validation["outside_image_indices"]
    )


def pivot_digest(representatives: list[int]) -> str:
    text = ",".join(str(value) for value in representatives)
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def sector_observability_record(
    L: int,
    m: int,
    orbit_records: list[dict[str, int]],
    annihilator_record: dict[str, object],
) -> dict[str, object]:
    sector = 2 * m
    sector_orbits = [
        row for row in orbit_records if row["particle_sector"] == sector
    ]
    low_orbits = [
        row for row in sector_orbits if row["absolute_leg_imbalance"] <= 2
    ]
    high_orbits = [
        row for row in sector_orbits if row["absolute_leg_imbalance"] >= 4
    ]
    high_index = {
        row["representative"]: index for index, row in enumerate(high_orbits)
    }
    high_sizes = {
        row["representative"]: row["orbit_size"] for row in high_orbits
    }

    rows: list[dict[int, int]] = []
    sector_basis_count = 0
    for vector in annihilator_record["basis"]:
        vector_sector = (
            int(vector[0]["representative"]).bit_count() if vector else -1
        )
        if vector_sector != sector:
            continue
        sector_basis_count += 1
        restricted: dict[int, int] = {}
        for term in vector:
            representative = int(term["representative"])
            column = high_index.get(representative)
            if column is not None:
                # The configuration pairing is diagonal in the orbit-sum basis
                # with the nonzero weight |orbit|.
                restricted[column] = (
                    int(term["coefficient"]) * high_sizes[representative]
                )
        rows.append(restricted)

    W_by_sector = {
        int(key): int(value)
        for key, value in annihilator_record["W_sectors"].items()
    }
    W_dimension = W_by_sector.get(sector, 0)
    K_dimension = len(sector_orbits)
    U_dimension = K_dimension - W_dimension
    n = comb(L, m)
    comparison_target = n * (n + 1) // 2 - int(L % 4 == 0 and 2 * m == L)

    modular: dict[str, dict[str, object]] = {}
    for prime in PRIMES:
        rank, pivots = sparse_modular_rank(rows, prime)
        pivot_representatives = [
            high_orbits[index]["representative"] for index in pivots
        ]
        modular[str(prime)] = {
            "tag": "[COMPUTATION]",
            "rank_Fp_lower_bound_for_rank_Q": rank,
            "rational_rank_ceiling_omitted_coordinate_dimension": len(high_orbits),
            "lower_equals_ceiling": rank == len(high_orbits),
            "pivot_representatives_sha256": pivot_digest(pivot_representatives),
            "first_pivot_representative": (
                pivot_representatives[0] if pivot_representatives else None
            ),
            "last_pivot_representative": (
                pivot_representatives[-1] if pivot_representatives else None
            ),
        }

    exact_omitted_rank = all(
        int(record["rank_Fp_lower_bound_for_rank_Q"]) == len(high_orbits)
        for record in modular.values()
    )
    return {
        "tag": "[COMPUTATION] exact finite-row two-sided rank sandwich",
        "L": L,
        "m": m,
        "particle_sector": sector,
        "K_sector_dimension": K_dimension,
        "retained_abs_leg_imbalance_0_or_2_dimension": len(low_orbits),
        "omitted_abs_leg_imbalance_at_least_4_dimension": len(high_orbits),
        "annihilator_W_sector_dimension_exact_Q": W_dimension,
        "annihilator_basis_vectors_in_sector": sector_basis_count,
        "cyclic_U_sector_dimension_exact_Q": U_dimension,
        "prior_symmetric_square_value_for_comparison_only": comparison_target,
        "prior_value_matches_exact_U_rank": U_dimension == comparison_target,
        "weighted_annihilator_restriction": modular,
        "annihilator_restriction_rank_exact_Q": (
            len(high_orbits) if exact_omitted_rank else None
        ),
        "restriction_kernel_on_U_dimension_exact_Q": (
            0 if exact_omitted_rank else None
        ),
        "retained_restriction_rank_exact_Q": (
            U_dimension if exact_omitted_rank else None
        ),
        "retained_codomain_excess_over_U": len(low_orbits) - U_dimension,
        "all_checks_passed": (
            sector_basis_count == W_dimension
            and exact_omitted_rank
            and len(low_orbits) + len(high_orbits) == K_dimension
            and U_dimension == comparison_target
        ),
    }


# New rung basis.  Old local labels are 0=empty, 1=top, 2=bottom, 3=double.
NEW_TO_OLD: dict[str, tuple[tuple[int, Fraction], ...]] = {
    "0": ((0, Fraction(1)),),
    "b": ((1, Fraction(1)), (2, Fraction(1))),
    "a": ((1, Fraction(1)), (2, Fraction(-1))),
    "d": ((3, Fraction(1)),),
}
OLD_TO_NEW: dict[int, tuple[tuple[str, Fraction], ...]] = {
    0: (("0", Fraction(1)),),
    1: (("b", Fraction(1, 2)), ("a", Fraction(1, 2))),
    2: (("b", Fraction(1, 2)), ("a", Fraction(-1, 2))),
    3: (("d", Fraction(1)),),
}


def expand_new_pair(left: str, right: str) -> dict[tuple[int, int], Fraction]:
    output: dict[tuple[int, int], Fraction] = {}
    for old_left, left_coefficient in NEW_TO_OLD[left]:
        for old_right, right_coefficient in NEW_TO_OLD[right]:
            key = (old_left, old_right)
            output[key] = output.get(key, Fraction(0)) + (
                left_coefficient * right_coefficient
            )
    return output


def old_pair_to_new(
    vector: dict[tuple[int, int], Fraction]
) -> dict[tuple[str, str], Fraction]:
    output: dict[tuple[str, str], Fraction] = {}
    for (old_left, old_right), coefficient in vector.items():
        for new_left, left_coefficient in OLD_TO_NEW[old_left]:
            for new_right, right_coefficient in OLD_TO_NEW[old_right]:
                key = (new_left, new_right)
                output[key] = output.get(key, Fraction(0)) + (
                    coefficient * left_coefficient * right_coefficient
                )
    return {key: value for key, value in output.items() if value}


def scaled_horizontal_action(left: str, right: str) -> dict[tuple[str, str], int]:
    old_output: dict[tuple[int, int], Fraction] = {}
    for (old_left, old_right), coefficient in expand_new_pair(left, right).items():
        for mask in (1, 2):  # top-top and bottom-bottom edge toggles
            target = (old_left ^ mask, old_right ^ mask)
            old_output[target] = old_output.get(target, Fraction(0)) + coefficient
    transformed = old_pair_to_new(old_output)
    scaled = {key: 2 * value for key, value in transformed.items()}
    if any(value.denominator != 1 for value in scaled.values()):
        raise AssertionError((left, right, scaled))
    return {key: int(value) for key, value in scaled.items() if value}


def scaled_rung_action(label: str) -> dict[str, int]:
    old_output: dict[int, Fraction] = {}
    for old, coefficient in NEW_TO_OLD[label]:
        target = old ^ 3
        old_output[target] = old_output.get(target, Fraction(0)) + coefficient
    transformed: dict[str, Fraction] = {}
    for old, coefficient in old_output.items():
        for new, conversion in OLD_TO_NEW[old]:
            transformed[new] = transformed.get(new, Fraction(0)) + coefficient * conversion
    scaled = {key: 2 * value for key, value in transformed.items() if value}
    if any(value.denominator != 1 for value in scaled.values()):
        raise AssertionError((label, scaled))
    return {key: int(value) for key, value in scaled.items()}


def action_records() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    horizontal: list[dict[str, object]] = []
    for left, right in product(LABELS, repeat=2):
        action = scaled_horizontal_action(left, right)
        source_a = int(left == "a") + int(right == "a")
        source_particles = PARTICLE_WEIGHT[left] + PARTICLE_WEIGHT[right]
        outputs = []
        for (target_left, target_right), coefficient in sorted(action.items()):
            target_a = int(target_left == "a") + int(target_right == "a")
            target_particles = (
                PARTICLE_WEIGHT[target_left] + PARTICLE_WEIGHT[target_right]
            )
            outputs.append(
                {
                    "target": [target_left, target_right],
                    "coefficient": coefficient,
                    "delta_antibonding_count": target_a - source_a,
                    "delta_particle_number": target_particles - source_particles,
                }
            )
        horizontal.append({"source": [left, right], "outputs": outputs})

    rung: list[dict[str, object]] = []
    for label in LABELS:
        outputs = [
            {
                "target": target,
                "coefficient": coefficient,
                "delta_antibonding_count": int(target == "a") - int(label == "a"),
                "delta_particle_number": PARTICLE_WEIGHT[target]
                - PARTICLE_WEIGHT[label],
            }
            for target, coefficient in sorted(scaled_rung_action(label).items())
        ]
        rung.append({"source": label, "outputs": outputs})
    return horizontal, rung


def transform_vector(
    vector: dict[int, int], L: int, transform: str
) -> dict[int, int]:
    function = tau_config if transform == "tau" else rho_config
    output: dict[int, int] = {}
    for config, coefficient in vector.items():
        image = function(config, L)
        output[image] = output.get(image, 0) + coefficient
    return output


def antibonding_product(L: int) -> dict[int, int]:
    vector = {0: 1}
    for rung in range(L):
        next_vector: dict[int, int] = {}
        for config, coefficient in vector.items():
            next_vector[config | (1 << (2 * rung))] = coefficient
            next_vector[config | (1 << (2 * rung + 1))] = -coefficient
        vector = next_vector
    return vector


def pairing_with_annihilator(
    vector: dict[int, int], basis_vector: list[dict[str, int]], L: int
) -> int:
    value = 0
    for term in basis_vector:
        representative = int(term["representative"])
        reflected = rho_config(representative, L)
        orbit = {
            representative,
            tau_config(representative, L),
            reflected,
            tau_config(reflected, L),
        }
        value += int(term["coefficient"]) * sum(
            vector.get(config, 0) for config in orbit
        )
    return value


def transformed_orbit_census(L: int, particle_sector: int) -> dict[str, int]:
    seen: set[tuple[str, ...]] = set()
    total = 0
    omitted = 0
    for word in product(LABELS, repeat=L):
        if sum(PARTICLE_WEIGHT[label] for label in word) != particle_sector:
            continue
        if sum(label == "a" for label in word) % 2:
            continue
        reflected = tuple(reversed(word))
        canonical = min(word, reflected)
        if canonical in seen:
            continue
        seen.add(canonical)
        total += 1
        if sum(label == "a" for label in word) >= 4:
            omitted += 1
    return {"tau_plus_rho_orbit_dimension": total, "antibonding_at_least_4_dimension": omitted}


def make_check(
    checks: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def main() -> int:
    started = time.process_time()
    source_hash = file_sha256(SOURCE_PATH)
    source = json.loads(SOURCE_PATH.read_text())
    records = source_records(source)
    checks: list[dict[str, object]] = []

    source_ok = (
        all(bool(check["passed"]) for check in source["checks"])
        and set(records) == set(L_VALUES)
        and all(source_exact(records[L]) for L in L_VALUES)
    )
    make_check(
        checks,
        "C1_source_exact_annihilators_L3_8",
        source_ok,
        "the cited source has exact-Q invariant annihilator bases and matching cyclic sandwiches for every L=3..8",
    )
    make_check(
        checks,
        "C2_good_primes",
        all(is_prime(prime) and prime > 4 for prime in PRIMES),
        f"odd primes {PRIMES}; orbit weights 1,2,4 remain nonzero",
    )

    observability: list[dict[str, object]] = []
    orbit_counts: dict[int, int] = {}
    for L in L_VALUES:
        orbits = configuration_orbits(L)
        orbit_counts[L] = len(orbits)
        for m in range(L + 1):
            observability.append(
                sector_observability_record(L, m, orbits, records[L])
            )

    make_check(
        checks,
        "C3_complete_sector_coverage_L3_8",
        {(row["L"], row["m"]) for row in observability}
        == {(L, m) for L in L_VALUES for m in range(L + 1)},
        f"all {len(observability)} required (L,m) sectors are present",
    )
    make_check(
        checks,
        "C4_direct_K_orbit_census",
        all(orbit_counts[L] == int(records[L]["K_dim"]) for L in L_VALUES),
        f"direct orbit dimensions {orbit_counts}",
    )
    for prime in PRIMES:
        make_check(
            checks,
            f"C5_omitted_pairing_full_column_rank_mod_{prime}",
            all(
                row["weighted_annihilator_restriction"][str(prime)][
                    "lower_equals_ceiling"
                ]
                for row in observability
            ),
            "every modular lower bound meets the omitted-coordinate rational ceiling",
        )
    make_check(
        checks,
        "C6_two_slice_restriction_injective_exact_Q",
        all(
            row["all_checks_passed"]
            and row["restriction_kernel_on_U_dimension_exact_Q"] == 0
            and row["retained_restriction_rank_exact_Q"]
            == row["cyclic_U_sector_dimension_exact_Q"]
            for row in observability
        ),
        "duality gives ker(R_|Delta|<=2 on U)=0 in every exact sector L=3..8",
    )

    exceptional = {
        f"L{L}_m{m}": next(
            row for row in observability if row["L"] == L and row["m"] == m
        )
        for L, m in ((4, 2), (8, 4))
    }
    make_check(
        checks,
        "C7_exceptional_sector_4_2",
        exceptional["L4_m2"]["retained_restriction_rank_exact_Q"] == 20
        and exceptional["L4_m2"][
            "annihilator_restriction_rank_exact_Q"
        ]
        == exceptional["L4_m2"][
            "omitted_abs_leg_imbalance_at_least_4_dimension"
        ]
        == 1,
        "(4,2): retained rank 20 and omitted pairing rank 1/1",
    )
    make_check(
        checks,
        "C8_exceptional_sector_8_4",
        exceptional["L8_m4"]["retained_restriction_rank_exact_Q"] == 2484
        and exceptional["L8_m4"][
            "annihilator_restriction_rank_exact_Q"
        ]
        == exceptional["L8_m4"][
            "omitted_abs_leg_imbalance_at_least_4_dimension"
        ]
        == 433,
        "(8,4): retained rank 2484 and omitted pairing rank 433/433",
    )

    horizontal_actions, rung_actions = action_records()
    all_action_outputs = [
        output
        for record in horizontal_actions + rung_actions
        for output in record["outputs"]
    ]
    make_check(
        checks,
        "C9_bonding_action_complete_and_integral",
        len(horizontal_actions) == 16
        and len(rung_actions) == 4
        and all(isinstance(output["coefficient"], int) for output in all_action_outputs),
        "all 16 horizontal and 4 rung local source states were evaluated for 2B",
    )
    make_check(
        checks,
        "C10_bonding_action_block_tridiagonal",
        all(
            output["delta_antibonding_count"] in {-2, 0, 2}
            and output["delta_particle_number"] in {-2, 0, 2}
            for output in all_action_outputs
        )
        and next(
            record for record in horizontal_actions if record["source"] == ["b", "a"]
        )["outputs"]
        == []
        and next(
            record for record in horizontal_actions if record["source"] == ["a", "b"]
        )["outputs"]
        == [],
        "each D/F/D+ grade changes the antibonding count only by 0 or 2 in absolute value; ba and ab are killed",
    )

    transformed_early = {
        f"L{L}_k{sector}": transformed_orbit_census(L, sector)
        for L, sector in ((3, 0), (3, 2), (3, 4), (3, 6), (4, 0), (4, 2), (4, 4))
    }
    witness = antibonding_product(4)
    W4_sector4 = [
        vector
        for vector in records[4]["basis"]
        if int(vector[0]["representative"]).bit_count() == 4
    ]
    witness_pairings = [
        pairing_with_annihilator(witness, vector, 4) for vector in W4_sector4
    ]
    witness_invariant = (
        transform_vector(witness, 4, "tau") == witness
        and transform_vector(witness, 4, "rho") == witness
    )
    first_omitted = transformed_early["L4_k4"][
        "antibonding_at_least_4_dimension"
    ]
    earlier_omitted = [
        row["antibonding_at_least_4_dimension"]
        for key, row in transformed_early.items()
        if key != "L4_k4"
    ]
    make_check(
        checks,
        "C11_bonding_truncation_first_possible_row",
        first_omitted == 1 and all(value == 0 for value in earlier_omitted),
        "the #a<=2 truncation is the identity before (L,m)=(4,2), where its omitted K-line is aaaa",
    )
    make_check(
        checks,
        "C12_bonding_witness_exact_symmetry_and_sector",
        len(witness) == 16
        and witness_invariant
        and all(config.bit_count() == 4 for config in witness)
        and set(witness.values()) == {-1, 1},
        "a0 a1 a2 a3 expands to a nonzero tau/rho-invariant 16-term integer vector in sector four",
    )
    make_check(
        checks,
        "C13_bonding_witness_in_U_by_exact_duality",
        len(W4_sector4) == 2
        and witness_pairings == [0, 0]
        and source_exact(records[4]),
        "the witness pairs to zero with the complete exact W_4 basis, hence belongs to U_4^(4)=W_4^perp",
    )
    make_check(
        checks,
        "C14_bonding_truncation_exact_kernel",
        first_omitted == 1
        and witness_pairings == [0, 0]
        and exceptional["L4_m2"]["cyclic_U_sector_dimension_exact_Q"] == 20,
        "the omitted line is contained in U, so the #a<=2 restriction has exact kernel dimension 1 and exact rank 19",
    )

    rss = peak_rss_bytes()
    make_check(
        checks,
        "C15_resource_ceiling",
        rss < RSS_CAP_BYTES,
        f"producer peak RSS {rss} bytes is below {RSS_CAP_BYTES}",
    )
    names = [str(check["name"]) for check in checks]
    if len(names) != len(set(names)):
        raise AssertionError("check names are not unique")
    if not all(bool(check["passed"]) for check in checks):
        failed = [check["name"] for check in checks if not check["passed"]]
        raise AssertionError(f"failed checks: {failed}")

    data = {
        "headline": {
            "tag": "[LEMMA][COMPUTATION][THEOREM][UNRESOLVED]",
            "result": (
                "[COMPUTATION] For every sector L=3..8, restriction of U_L^(2m) "
                "to |n_top-n_bottom|<=2 is injective exactly over Q. "
                "[THEOREM] The bonding truncation #a<=2 first fails at (4,2), "
                "with kernel span{a0 a1 a2 a3}. [UNRESOLVED] No all-L rank claim."
            ),
        },
        "duality_criterion": {
            "tag": "[LEMMA]",
            "statement": (
                "Let K=E direct_sum H be the orthogonal split into retained and omitted "
                "configuration coordinates, U<=K, and W=U^perp. For restriction R_E:U->E, "
                "dim ker R_E = dim H - rank(W->H*, w |-> <w, .>|_H). Thus R_E is injective "
                "iff the weighted restriction of W has full omitted-column rank."
            ),
            "rank_identity": "rank(R_E|U)=dim(U)-dim(H)+rank(W->H*)",
            "modular_discipline": (
                "rank_Fp(W->H*) is used only as a lower bound for rank_Q; dim(H) is the "
                "independent rational ceiling, and exactness is asserted only when they match"
            ),
        },
        "two_slice_observability_exact_L3_8": observability,
        "exceptional_rows": exceptional,
        "bonding_antibonding": {
            "tag": "[LEMMA][THEOREM]",
            "basis": ["0=empty", "b=top+bottom", "a=top-bottom", "d=double"],
            "scaled_operator": "2 B_L; scaling by 2 does not change cyclic subspaces over Q",
            "rung_actions": rung_actions,
            "horizontal_actions": horizontal_actions,
            "triangular_lemma": (
                "[LEMMA] In the rungwise 0,b,a,d basis, 2D, 2F, and 2D+ have integer "
                "coefficients and change the number of a factors only by 0 or +/-2. "
                "Therefore the tau-plus block is even-a and block-tridiagonal in a-count."
            ),
            "first_truncation_failure": {
                "tag": "[THEOREM]",
                "map_class": (
                    "all linear evaluations Phi that factor as Phi=T after Q_<=2, where "
                    "Q_<=2 discards transformed coordinates with at least four a factors"
                ),
                "lexicographically_first_row": {"L": 4, "m": 2, "particle_sector": 4},
                "earlier_transformed_census": transformed_early,
                "omitted_K_dimension": 1,
                "omitted_basis_word": ["a", "a", "a", "a"],
                "physical_witness": [
                    {"configuration": config, "coefficient": coefficient}
                    for config, coefficient in sorted(witness.items())
                ],
                "pairings_with_complete_exact_W4_sector_basis": witness_pairings,
                "witness_in_U_exact_Q": True,
                "kernel_dimension_exact_Q": 1,
                "restriction_rank_exact_Q": 19,
                "closure": (
                    "[THEOREM] Every map in the named factor-through-Q_<=2 class kills the "
                    "nonzero vector a0 a1 a2 a3 in U_4^(4), so no such map is injective."
                ),
            },
        },
        "single_missing_all_L_lemma": {
            "tag": "[UNRESOLVED]",
            "statement": (
                "For every L,m, the weighted map W_L^(2m)->H_Lm* is onto and its kernel "
                "has dimension dim(E_Lm)-[C(L,m)(C(L,m)+1)/2-delta(L,m)]. This exact-sequence "
                "lemma would give both rank inclusions and two-slice injectivity. It is verified "
                "here only for L=3..8 and is not claimed for any larger L."
            ),
        },
        "coverage_contract": {
            "tag": "[COMPUTATION]",
            "required_L": list(L_VALUES),
            "required_sector_pairs": [[L, m] for L in L_VALUES for m in range(L + 1)],
            "required_exceptional_pairs": [[4, 2], [8, 4]],
            "required_horizontal_local_sources": 16,
            "required_rung_local_sources": 4,
            "all_L_computation_claimed": False,
        },
        "scope": [
            "[COMPUTATION] exact finite rows only: every sector L=3..8, with separate (4,2) and (8,4) records",
            "[LEMMA] the duality criterion and local bonding block-tridiagonality are general finite-dimensional/local identities",
            "[THEOREM] the factor-through-#a<=2 map class is closed by the exact (4,2) kernel witness",
            "[UNRESOLVED] the symmetric-square rank law and the two-slice exact sequence for L>=9",
            "[UNRESOLVED] no L=9 or L=10 closure was launched; K_c was not used",
        ],
    }
    envelope = {
        "meta": {
            "script": SCRIPT,
            "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "python": sys.version.split()[0],
            "artifact_schema": "meta/data/checks-v1",
            "exact_arithmetic": (
                "Python integers/Fraction; sparse ranks over two exact prime fields. "
                "Every modular lower bound promoted to rank_Q is matched by the omitted-column ceiling."
            ),
            "source_artifacts": {
                str(SOURCE_PATH.relative_to(ROOT)): source_hash,
            },
            "primes": list(PRIMES),
            "cpu_seconds": round(time.process_time() - started, 6),
            "peak_rss_bytes": rss,
            "rss_cap_bytes": RSS_CAP_BYTES,
        },
        "data": data,
        "checks": checks,
    }
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUT_PATH.with_suffix(".tmp")
    temporary.write_text(json.dumps(envelope, indent=2, sort_keys=True) + "\n")
    temporary.replace(OUT_PATH)
    print(
        f"[e230] wrote {OUT_PATH.relative_to(ROOT)}: "
        f"{len(checks)}/{len(checks)} checks, {len(observability)} exact sectors, "
        f"cpu={envelope['meta']['cpu_seconds']}s rss={rss}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
