"""Standalone checks for the streamed symmetry-reduced Lee--Yang evaluator.

Run from the repository root:
    timeout 1800 .venv/bin/python tests/test_lee_yang_53.py
"""

from __future__ import annotations

import importlib.util
from fractions import Fraction
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from ising.lee_yang import _CRT_PRIMES, rational_field_polynomial_transfer

ROOT = Path(__file__).resolve().parents[1]
PRODUCER = ROOT / "experiments" / "e98_lee_yang_53.py"


def load_producer():
    if not PRODUCER.exists():
        return None
    spec = importlib.util.spec_from_file_location("e98_lee_yang_53", PRODUCER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def independent_q_value(coefficients: list[int], t: Fraction) -> Fraction:
    """Evaluate the circle polynomial directly from a full coefficient vector."""

    values = [int(value) for value in coefficients]
    degree = len(values) - 1
    assert values == values[::-1]
    if degree % 2:
        quotient = [values[0]]
        for coefficient in values[1:-1]:
            quotient.append(coefficient - quotient[-1])
        assert values[-1] == quotient[-1]
        values = quotient
        degree -= 1
    middle = degree // 2
    result = Fraction(values[middle])
    previous, current = Fraction(1), t
    for order in range(1, middle + 1):
        if order == 1:
            chebyshev = current
        else:
            chebyshev = 2 * t * current - previous
            previous, current = current, chebyshev
        result += 2 * values[middle + order] * chebyshev
    return result


def direct_spatial_spin_orbits_4x4() -> int:
    """Brute-force D8 plus complement orbits, independent of e98 maps."""

    side = 4
    state_count = 1 << (side * side)
    mask = state_count - 1

    def rotate90(value: int) -> int:
        out = 0
        for row in range(side):
            for col in range(side):
                source = row * side + col
                target = col * side + (side - 1 - row)
                out |= ((value >> source) & 1) << target
        return out

    def reflect_row(value: int) -> int:
        out = 0
        for row in range(side):
            for col in range(side):
                source = row * side + col
                target = (side - 1 - row) * side + col
                out |= ((value >> source) & 1) << target
        return out

    seen = bytearray(state_count)
    count = 0
    for state in range(state_count):
        if seen[state]:
            continue
        count += 1
        rotations = [state]
        for _ in range(3):
            rotations.append(rotate90(rotations[-1]))
        orbit = set(rotations)
        orbit.update(reflect_row(value) for value in rotations)
        orbit.update(mask ^ value for value in tuple(orbit))
        for value in orbit:
            seen[value] = 1
    assert all(seen)
    return count


def main() -> None:
    module = load_producer()
    assert module is not None, "e98 must provide the streamed reduced evaluator"

    symmetry = module.build_symmetry_data((4, 4))
    assert len(symmetry["representatives"]) == direct_spatial_spin_orbits_4x4()
    assert symmetry["formula"]["spatial_plus_spin_orbits"] == 4324
    assert int(symmetry["count0"].sum() + symmetry["count1"].sum()) == 1 << 16
    coefficients = rational_field_polynomial_transfer((4, 4, 4), 2, 3)
    states = np.arange(1 << 16, dtype=np.uint32)
    broken, down = module.state_statistics_for_states(states, (4, 4))
    identity_symmetry = {
        "cross": (4, 4),
        "ns": 16,
        "n_states": 1 << 16,
        "representatives": states,
        "state_to_orbit": states,
        "flip_indices": np.empty(0, dtype=np.uint32),
        "count0": np.ones(1 << 16, dtype=np.uint64),
        "count1": np.zeros(1 << 16, dtype=np.uint64),
        "broken": broken,
        "down": down,
    }
    raw_residue = module.modular_reduced_circle_residue(
        (4, 4, 4),
        2,
        3,
        Fraction(0),
        _CRT_PRIMES[0],
        identity_symmetry,
        strip_pairs=1 << 12,
    )
    assert raw_residue["residue"] == independent_q_value(coefficients, Fraction(0)).numerator % _CRT_PRIMES[0]

    for t in (Fraction(0), Fraction(1, 2)):
        expected = independent_q_value(coefficients, t)
        streamed = module.exact_reduced_circle_value(
            (4, 4, 4), 2, 3, t, strip_pairs=1 << 12
        )
        assert streamed["value"] == expected
        assert streamed["butterfly_strips"] == (
            streamed["crt_primes_used"] * 2 * 3 * module.butterfly_strip_count(16, 1 << 12)
        )
    # The target cube has odd N=125, so independently exercise the z+1
    # removal and the centered odd-degree circle identity on 3^3.
    odd_coefficients = rational_field_polynomial_transfer((3, 3, 3), 2, 3)
    for t in (Fraction(0), Fraction(1, 2)):
        odd_streamed = module.exact_reduced_circle_value(
            (3, 3, 3), 2, 3, t, strip_pairs=1 << 12
        )
        assert odd_streamed["value"] == independent_q_value(odd_coefficients, t)
    with TemporaryDirectory() as directory:
        disk_streamed = module.exact_reduced_circle_value(
            (3, 3, 3),
            2,
            3,
            Fraction(1, 2),
            strip_pairs=1 << 12,
            disk_dir=Path(directory),
        )
        assert disk_streamed["value"] == independent_q_value(odd_coefficients, Fraction(1, 2))
        assert not list(Path(directory).iterdir())

    forward = module.modular_reduced_circle_residue(
        (4, 4, 4),
        2,
        3,
        Fraction(1, 2),
        _CRT_PRIMES[0],
        symmetry,
        strip_pairs=1 << 12,
        strip_order="forward",
    )
    reverse = module.modular_reduced_circle_residue(
        (4, 4, 4),
        2,
        3,
        Fraction(1, 2),
        _CRT_PRIMES[0],
        symmetry,
        strip_pairs=1 << 12,
        strip_order="reverse",
    )
    assert forward["residue"] == reverse["residue"]
    assert forward["butterfly_strips"] == reverse["butterfly_strips"]
    print("PASS")


if __name__ == "__main__":
    main()
