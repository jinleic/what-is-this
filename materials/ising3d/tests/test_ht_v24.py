"""Independent exact checks for experiments/e68_ht_v24.py."""

from __future__ import annotations

import importlib.util
import json
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "e68_ht_v24", ROOT / "experiments" / "e68_ht_v24.py"
)
assert SPEC is not None and SPEC.loader is not None
E68 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(E68)


def _check(name: str, passed: bool, detail: str) -> None:
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}")
    if not passed:
        raise AssertionError(f"{name}: {detail}")


def _frontier_weight(side: int) -> int:
    result = E68.frontier_even_polynomial(side, modulus=E68.CRT_PRODUCT)
    subsets = E68.cut_polynomial_as_subsets(
        result["coefficients"], side, E68.CRT_PRODUCT
    )
    residue = E68.square_free_log_full(subsets, E68.CRT_PRODUCT)
    return residue if residue <= E68.CRT_PRODUCT // 2 else residue - E68.CRT_PRODUCT


def main() -> None:
    # This directly enumerates 3x3x3 edge frontiers.  The reference is computed
    # independently by the old spin-transfer/broken-bond finite-lattice route.
    frontier_weight = _frontier_weight(3)
    ordinary = E68.flm.high_temperature_free_energy(3, 12)
    _check(
        "independent lower-order enumeration",
        frontier_weight == 9_188
        and ordinary.box_weights[(3, 3, 3)][12] == frontier_weight,
        "edge-frontier and spin-transfer FLM both give the 3x3x3 degree-12 weight 9188",
    )

    payload = json.loads(
        (ROOT / "results" / "series" / "ht_v24.json").read_text(encoding="utf-8")
    )
    data = payload["data"]
    bound = int(data["minimal_cube"]["absolute_bound"])
    modulus = int(data["crt"]["modulus_product"])
    weight = int(data["minimal_cube"]["exact_weight"])
    _check(
        "CRT uniqueness certificate",
        E68.minimal_cube_weight_bound(5) == bound
        and E68.CRT_PRODUCT == modulus > 2 * bound
        and all(E68._is_prime_32(prime) for prime in E68.CRT_PRIMES)
        and len(set(E68.CRT_PRIMES)) == len(E68.CRT_PRIMES)
        and data["crt"]["deterministic_primality_checked"]
        and abs(weight) <= bound
        and int(data["crt"]["residue"]) == weight % modulus,
        f"M={modulus} exceeds 2B={2 * bound} and reconstructs weight {weight}",
    )

    canonical = json.loads(
        (
            ROOT
            / "results"
            / "series"
            / "extended2_sc_ht_free_energy.json"
        ).read_text(encoding="utf-8")
    )
    coefficient = tuple(
        Fraction(value) for value in data["series"]["coefficients"]
    )
    prefix = tuple(
        Fraction(value) for value in canonical["data"]["coefficients"]
    )
    _check(
        "canonical prefix and v24",
        coefficient[: len(prefix)] == prefix
        and coefficient[24] == Fraction(2_135_670_379_057, 8)
        and data["series"]["interaction_v24"] == "266958797382",
        "v^0..v^22 match exactly and v^24=2135670379057/8",
    )
    _check(
        "artifact envelope and checks",
        set(payload) == {"provenance", "data", "checks"}
        and payload["provenance"]["script"] == "experiments/e68_ht_v24.py"
        and payload["provenance"]["interpreter"] == ".venv/bin/python"
        and all(check["passed"] for check in payload["checks"]),
        "provenance/data/checks envelope is complete and every embedded check passes",
    )
    print("PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL test_ht_v24: {error}")
        raise
