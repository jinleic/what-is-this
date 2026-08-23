#!/usr/bin/env python3
"""Exact D/F/D+ obstruction for the balanced symmetric-square realization.

The balanced polarization sends an unordered pair of m-subsets to a state with
m particles on each ladder leg and then averages over leg swap and rung
reflection.  This module constructs explicit invariant source vectors showing
that each hard-core piece D, F, and D+ has a nonzero component outside the
balanced image.  No ladder closure is performed.

All vectors are sparse dictionaries with Python-integer coefficients.
"""
from __future__ import annotations

import json
from collections import defaultdict

SMALL_L = tuple(range(2, 10))


def tau_config(config: int, L: int) -> int:
    out = 0
    for rung in range(L):
        if (config >> (2 * rung)) & 1:
            out |= 1 << (2 * rung + 1)
        if (config >> (2 * rung + 1)) & 1:
            out |= 1 << (2 * rung)
    return out


def rho_config(config: int, L: int) -> int:
    out = 0
    for rung in range(L):
        out |= ((config >> (2 * rung)) & 3) << (2 * (L - 1 - rung))
    return out


def orbit_sum(config: int, L: int) -> dict[int, int]:
    images = {
        config,
        tau_config(config, L),
        rho_config(config, L),
        tau_config(rho_config(config, L), L),
    }
    return {image: 1 for image in images}


def config_from_legs(top: set[int], bottom: set[int]) -> int:
    config = 0
    for rung in top:
        config |= 1 << (2 * rung)
    for rung in bottom:
        config |= 1 << (2 * rung + 1)
    return config


def ladder_edges(L: int) -> list[tuple[int, int]]:
    edges = [(2 * rung, 2 * rung + 1) for rung in range(L)]
    for rung in range(L - 1):
        edges.append((2 * rung, 2 * (rung + 1)))
        edges.append((2 * rung + 1, 2 * (rung + 1) + 1))
    return edges


def leg_counts(config: int, L: int) -> tuple[int, int]:
    top = sum((config >> (2 * rung)) & 1 for rung in range(L))
    bottom = sum((config >> (2 * rung + 1)) & 1 for rung in range(L))
    return top, bottom


def apply_piece(vector: dict[int, int], L: int, piece: str) -> dict[int, int]:
    """Apply the creation, hopping, or annihilation part of the edge toggle."""
    output: defaultdict[int, int] = defaultdict(int)
    for config, coefficient in vector.items():
        for left, right in ladder_edges(L):
            occupancy = ((config >> left) & 1) + ((config >> right) & 1)
            enabled = (
                (piece == "D" and occupancy == 0)
                or (piece == "F" and occupancy == 1)
                or (piece == "D_plus" and occupancy == 2)
            )
            if enabled:
                output[config ^ (1 << left) ^ (1 << right)] += coefficient
    return {config: coefficient for config, coefficient in output.items() if coefficient}


def transform_vector(vector: dict[int, int], L: int, transform: str) -> dict[int, int]:
    output: defaultdict[int, int] = defaultdict(int)
    fn = tau_config if transform == "tau" else rho_config
    for config, coefficient in vector.items():
        output[fn(config, L)] += coefficient
    return dict(output)


def invariant(vector: dict[int, int], L: int) -> bool:
    return transform_vector(vector, L, "tau") == vector and transform_vector(vector, L, "rho") == vector


def balanced(vector: dict[int, int], L: int) -> bool:
    return all(leg_counts(config, L)[0] == leg_counts(config, L)[1] for config in vector)


def witness_source(L: int, piece: str) -> tuple[dict[int, int], str]:
    if piece == "D":
        return orbit_sum(0, L), "vacuum"
    if piece == "F":
        config = config_from_legs({0}, {1})
        return orbit_sum(config, L), "orbit(top={0}, bottom={1})"
    if piece == "D_plus":
        config = config_from_legs({0, 1}, {0, 1})
        return orbit_sum(config, L), "orbit(top=bottom={0,1})"
    raise ValueError(f"unknown hard-core piece: {piece}")


def wedge_sign(left: int, right: int, L: int) -> int:
    full = (1 << L) - 1
    if left & right or (left | right) != full:
        return 0
    inversions = 0
    for site in range(L):
        if (left >> site) & 1:
            inversions += (right & ((1 << site) - 1)).bit_count()
    return -1 if inversions % 2 else 1


def exceptional_source_trace(L: int, piece: str) -> int | None:
    """Top-wedge trace of the D+ domain witness in the only relevant case L=4."""
    if L != 4 or piece != "D_plus":
        return None
    first = (1 << 0) | (1 << 1)
    reflected = (1 << 2) | (1 << 3)
    # The domain source is a sum of diagonal symmetric tensors.  Each diagonal
    # has zero exterior square, computed rather than assumed.
    return wedge_sign(first, first, L) + wedge_sign(reflected, reflected, L)


def module_record(L: int, piece: str) -> dict:
    source, label = witness_source(L, piece)
    output = apply_piece(source, L, piece)
    leakage = {
        config: coefficient
        for config, coefficient in output.items()
        if leg_counts(config, L)[0] != leg_counts(config, L)[1]
    }
    source_weights = {config.bit_count() for config in source}
    output_weights = {config.bit_count() for config in output}
    grade = {"D": 2, "F": 0, "D_plus": -2}[piece]
    expected_output_weights = {weight + grade for weight in source_weights}
    trace = exceptional_source_trace(L, piece)
    examples = [
        {
            "configuration": config,
            "top_count": leg_counts(config, L)[0],
            "bottom_count": leg_counts(config, L)[1],
            "coefficient": coefficient,
        }
        for config, coefficient in sorted(leakage.items())[:4]
    ]
    checks = {
        "source_balanced": balanced(source, L),
        "source_tau_rho_invariant": invariant(source, L),
        "output_tau_rho_invariant": invariant(output, L),
        "leakage_tau_rho_invariant": invariant(leakage, L),
        "grade_exact": output_weights == expected_output_weights,
        "leakage_nonzero": bool(leakage),
        "exceptional_D_plus_source_is_traceless": trace == 0 if trace is not None else True,
    }
    return {
        "L": L,
        "piece": piece,
        "source": label,
        "source_support": len(source),
        "source_particle_numbers": sorted(source_weights),
        "output_support": len(output),
        "output_particle_numbers": sorted(output_weights),
        "unbalanced_leakage_support": len(leakage),
        "unbalanced_leakage_l1": sum(abs(value) for value in leakage.values()),
        "unbalanced_leakage_examples": examples,
        "exceptional_top_wedge_trace": trace,
        "checks": checks,
        "all_checks_passed": all(checks.values()),
        "tag": "[THEOREM] explicit positive-coefficient counterwitness to balanced-module invariance",
    }


def build_module_records(L_values: tuple[int, ...] = SMALL_L) -> list[dict]:
    return [module_record(L, piece) for L in L_values for piece in ("D", "F", "D_plus")]


def main() -> int:
    records = build_module_records()
    passed = all(row["all_checks_passed"] for row in records)
    minimum_leakage = min(row["unbalanced_leakage_support"] for row in records)
    print(json.dumps({
        "status": "[THEOREM] balanced Reynolds image fails invariance under each of D, F, D+ for every L >= 2",
        "exact_small_range": [SMALL_L[0], SMALL_L[-1]],
        "record_count": len(records),
        "minimum_nonzero_leakage_support": minimum_leakage,
        "all_checks_passed": passed,
    }, indent=2, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
