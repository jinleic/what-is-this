"""Standalone exact checks for experiments/e86_ht_v26.py."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from fractions import Fraction
from math import gcd, prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "e86_ht_v26", ROOT / "experiments" / "e86_ht_v26.py"
)
assert SPEC is not None and SPEC.loader is not None
E86 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(E86)
sys.path.insert(0, str(ROOT / "experiments"))
from e18_series_extend import _hybrid_flm_engine


def _check(name: str, passed: bool, detail: str) -> None:
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}")
    if not passed:
        raise AssertionError(f"{name}: {detail}")


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _worker_telemetry_ok(row: dict) -> bool:
    timeline = row["timeline"]
    return (
        bool(timeline)
        and row["peak_states"] == max(entry["states"] for entry in timeline)
        and row["state_count_sha256"]
        == E86._sha256_json([entry["states"] for entry in timeline])
        and int(row["max_rss_bytes"]) > 0
        and float(row["wall_seconds"]) > 0
        and float(row["parent_wall_seconds"]) >= float(row["wall_seconds"])
    )


def main() -> None:
    # Independent lower-order route.  At degree 14 in 3x3x3 one cut has four
    # crossings and the other five have two.  The ordinary mixed logarithm
    # must retain a 2+2 split across two log factors.
    special333 = E86.extract_frontier_weights((3, 3, 3), special_cut=0)
    ordinary333 = E86.flm.high_temperature_free_energy(3, 14)
    w333 = ordinary333.box_weights[(3, 3, 3)][14]
    _check(
        "independent 3x3x3 four-crossing route",
        6 * special333["four_cut_centered"] == w333,
        f"6*c4={6 * special333['four_cut_centered']} equals FLM {w333}",
    )

    # A second lower-order route tests the distinct edge/middle cut orbits.
    cap444_edge = E86.extract_frontier_weights((4, 4, 4), special_cut=0)
    cap444_middle = E86.extract_frontier_weights((4, 4, 4), special_cut=1)
    with _hybrid_flm_engine():
        ordinary444 = E86.flm.high_temperature_free_energy(3, 20)
    w444 = ordinary444.box_weights[(4, 4, 4)][20]
    orbit444 = (
        6 * cap444_edge["four_cut_centered"]
        + 3 * cap444_middle["four_cut_centered"]
    )
    _check(
        "independent 4x4x4 two-orbit route",
        cap444_edge["minimal_centered"] == 8_655_072
        and orbit444 == w444,
        f"6*c_edge+3*c_middle={orbit444} equals FLM {w444}",
    )

    payload = _load(ROOT / "results" / "series" / "ht_v26.json")
    data = payload["data"]
    weights = data["weights"]
    crt = data["crt"]
    modulus = int(crt["modulus_product"])
    w655 = int(weights["w_655"])
    w664 = int(weights["w_664"])
    w555 = int(weights["w_555_26"])
    edge = int(weights["w_555_26_edge_coefficient"])
    middle = int(weights["w_555_26_middle_coefficient"])

    _check(
        "FLM class cancellation certificate",
        data["box_class_lemma"]["contributing_classes"] == "a+b+c <= 16"
        and data["box_class_lemma"]["minimal_span_class_count"] == 21
        and data["box_class_lemma"]["minimal_span_classes"]
        == [list(shape) for shape in E86.minimal_span_classes(E86.ORDER)]
        and data["box_class_lemma"]["contributing_canonical_classes"]
        == [
            list(shape)
            for shape in sorted(
                {
                    E86.flm._canonical_shape(shape)
                    for shape in E86.flm._ht_shapes(3, E86.ORDER, 0)
                }
            )
        ],
        "all and only canonical boxes with a+b+c <= 16 are represented",
    )

    expected_bounds = {
        "w_655": E86.frontier_log_bound((6, 5, 5)),
        "w_664": E86.frontier_log_bound((6, 6, 4)),
        "w_555_26_edge_coeff": E86.frontier_log_bound((5, 5, 5), 0),
        "w_555_26": 12 * E86.frontier_log_bound((5, 5, 5), 0),
    }
    certificate_rows = crt["per_weight_bounds"]
    _check(
        "CRT uniqueness certificates",
        tuple(crt["primes"]) == E86.CRT_PRIMES
        and prod(E86.CRT_PRIMES) == modulus
        and len(set(E86.CRT_PRIMES)) == len(E86.CRT_PRIMES)
        and all(E86._is_prime_32(prime) for prime in E86.CRT_PRIMES)
        and all(
            gcd(E86.CRT_PRIMES[left], E86.CRT_PRIMES[right]) == 1
            for left in range(len(E86.CRT_PRIMES))
            for right in range(left)
        )
        and all(
            int(certificate_rows[name]["bound"]) == bound
            and modulus > 2 * bound
            and abs(int(certificate_rows[name]["weight"])) <= bound
            for name, bound in expected_bounds.items()
        ),
        f"M={modulus} exceeds twice all independently rebuilt bounds",
    )
    _check(
        "CRT residue and cut-orbit reconstruction",
        w555 == 6 * (edge + middle)
        and int(weights["residues"]["w_655"]) == w655 % modulus
        and int(weights["residues"]["w_664"]) == w664 % modulus
        and int(weights["residues"]["w_555_26_edge"]) == edge % modulus
        and int(weights["residues"]["w_555_26_middle"]) == middle % modulus,
        "stored residues reconstruct every injected frontier weight",
    )

    v24 = _load(ROOT / "results" / "series" / "ht_v24.json")
    canonical = _load(ROOT / "results" / "series" / "extended2_sc_ht_free_energy.json")
    coefficient = tuple(Fraction(value) for value in data["series"]["coefficients"])
    v24_coefficient = tuple(
        Fraction(value) for value in v24["data"]["series"]["coefficients"]
    )
    prefix = tuple(Fraction(value) for value in canonical["data"]["coefficients"])
    _check(
        "prefix and v26 assembly",
        coefficient[: len(v24_coefficient)] == v24_coefficient
        and coefficient[: len(prefix)] == prefix
        and coefficient[25] == 0
        and coefficient[26].denominator == 26
        and coefficient[26]
        == Fraction(data["series"]["interaction_v26"]) + Fraction(3, 26),
        f"v^0..v^24 match and v^26={coefficient[26]}",
    )

    external = E86.EXTERNAL_ARISUE_FUJIWARA
    interaction = tuple(
        Fraction(value) for value in data["series"]["interaction_coefficients"]
    )
    _check(
        "external Arisue-Fujiwara Table I agreement",
        data["controls"]["external_arisue_fujiwara"]["agreement"]
        and all(
            interaction[degree] == Fraction(value)
            for degree, value in external["table_i_interaction"].items()
        )
        and data["controls"]["external_arisue_fujiwara"]["table_i_sha256"]
        == external["table_i_sha256"],
        "every tabulated interaction coefficient v^2..v^26 matches "
        "arXiv:hep-lat/0209002",
    )

    measured = data["frontier"]["measured"]
    _check(
        "frontier telemetry",
        all(
            row["probe_completed"]
            and not row["lower_bound_only"]
            and row["within_budget"]
            for row in data["frontier"]["projection"].values()
        )
        and _worker_telemetry_ok(measured["w655_worker"])
        and _worker_telemetry_ok(measured["w664_worker"])
        and _worker_telemetry_ok(measured["alt_sweep_worker"])
        and all(
            _worker_telemetry_ok(row)
            for row in measured["special_555_workers"]
        )
        and data["production_resource_measurement"]["maximum_resident_set_bytes"] > 0
        and float(data["production_resource_measurement"]["real_seconds"]) > 0,
        "state timelines, isolated RSS/walls, and producer RSS are coherent",
    )
    _check(
        "artifact envelope and embedded checks",
        set(payload) == {"provenance", "data", "checks"}
        and payload["provenance"]["script"] == "experiments/e86_ht_v26.py"
        and payload["provenance"]["interpreter"] == ".venv/bin/python"
        and not payload["provenance"]["benchmark_kc_used"]
        and data["series"]["v24_artifact_sha256"]
        == hashlib.sha256((ROOT / "results" / "series" / "ht_v24.json").read_bytes()).hexdigest()
        and all(check["passed"] for check in payload["checks"]),
        "provenance, artifact hashes, and all producer checks are present",
    )
    print("PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL test_ht_v26: {error}")
        raise
