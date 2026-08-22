"""Exact/lightweight sanity checks for the Onsager-quotient obstruction.

This script does not recompute the 2x3 Pauli closure.  It audits the two
stored finite-field quotient witnesses, verifies the cached primary-source
hashes, and checks small truncated-current algebras with exact integer
structure constants.
"""
from __future__ import annotations

import hashlib
import json
import platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "sources" / "manifest.yaml"
STRUCTURE = ROOT / "results" / "algebra_structure" / "structure.json"
RESULT = ROOT / "results" / "algebra_structure" / "oa_quotient.json"

EXPECTED_FULLTEXT = {
    "date_roan2000_onsager": {
        "path": "sources/fulltext/date_roan2000_onsager.pdf",
        "sha256": "c1334b58ed30c77b1592e628d8bf26c40d52ef63b9226ecad1ddec0e96771b06",
    },
    "roan1991_onsager": {
        "path": "sources/fulltext/roan1991_onsager.pdf",
        "sha256": "82988421997ae8fb388c79c4f46381def00ee15ecaaed04922bd7ab94c352313",
    },
}
EXPECTED_PRIMES = [2147483647, 2147483629]


def scalar(value: str) -> Any:
    value = value.strip()
    if value == "null":
        return None
    if value.startswith('"'):
        return json.loads(value)
    if re.fullmatch(r"-?[0-9]+", value):
        return int(value)
    return value


def manifest_entries(path: Path) -> dict[str, dict[str, Any]]:
    """Parse the scalar fields needed from the repository's simple YAML list."""
    entries: dict[str, dict[str, Any]] = {}
    current: dict[str, Any] | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        key_match = re.match(r"^\s*- key:\s*(.*?)\s*$", line)
        if key_match:
            key = str(scalar(key_match.group(1)))
            current = {"key": key}
            entries[key] = current
            continue
        field_match = re.match(r"^\s{4}([a-z][a-z0-9_]*):\s*(.*?)\s*$", line)
        if current is not None and field_match:
            current[field_match.group(1)] = scalar(field_match.group(2))
    return entries


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check_record(name: str, passed: bool, detail: str) -> dict[str, Any]:
    return {"name": name, "passed": bool(passed), "detail": detail}


def c7_witness(classification: dict[str, Any]) -> dict[str, Any] | None:
    matches = [
        item
        for item in classification["certified_quotient_images"]
        if item.get("field_type") == "C7"
    ]
    return matches[0] if len(matches) == 1 else None


# sl2 basis indices: e=0, f=1, h=2.  Values are (coefficient, basis index).
SL2_BRACKET: dict[tuple[int, int], tuple[int, int]] = {
    (0, 1): (1, 2),
    (1, 0): (-1, 2),
    (2, 0): (2, 0),
    (0, 2): (-2, 0),
    (2, 1): (-2, 1),
    (1, 2): (2, 1),
}


def truncated_current_lower_central(L: int) -> list[int]:
    """Dimensions of gamma_k(u C[u]/u^L tensor sl2), using exact brackets."""
    nilradical = {(degree, basis) for degree in range(1, L) for basis in range(3)}
    gamma = set(nilradical)
    dimensions: list[int] = []
    while gamma:
        dimensions.append(len(gamma))
        next_gamma: set[tuple[int, int]] = set()
        for degree_x, basis_x in nilradical:
            for degree_y, basis_y in gamma:
                degree = degree_x + degree_y
                bracket = SL2_BRACKET.get((basis_x, basis_y))
                if degree < L and bracket is not None and bracket[0] != 0:
                    next_gamma.add((degree, bracket[1]))
        gamma = next_gamma
    return dimensions


def main() -> int:
    entries = manifest_entries(MANIFEST)
    source_data: dict[str, Any] = {}
    source_ok = True
    for key, expected in EXPECTED_FULLTEXT.items():
        entry = entries.get(key, {})
        path = ROOT / expected["path"]
        actual_hash = sha256(path) if path.is_file() else None
        record_ok = (
            entry.get("obtained") == "full"
            and entry.get("local_path") == expected["path"]
            and entry.get("sha256") == expected["sha256"]
            and actual_hash == expected["sha256"]
        )
        source_ok &= record_ok
        source_data[key] = {
            "local_path": expected["path"],
            "manifest_sha256": entry.get("sha256"),
            "computed_sha256": actual_hash,
            "passed": record_ok,
        }

    abstract_status = {
        key: entries.get(key, {}).get("obtained")
        for key in ("davies1990_onsager", "davies1991_onsager")
    }
    abstract_ok = all(value == "abstract_only" for value in abstract_status.values())

    structure = json.loads(STRUCTURE.read_text(encoding="utf-8"))
    grid = structure["data"]["grid_2x3"]
    classifications = [grid["classification_prime_1"], grid["classification_prime_2"]]
    primes = [item["prime"] for item in classifications]
    witnesses = [c7_witness(item) for item in classifications]
    forms = [
        item["linkage_witnesses"]["invariant_form_witnesses"]["000"]
        for item in classifications
    ]
    witness_ok = all(
        witness == {
            "dimension": 105,
            "field_type": "C7",
            "scope": "F_p only; saturates sp(14) because the image preserves a nondegenerate alternating form and has dimension 105=dim sp(14)",
            "sector": "000",
        }
        for witness in witnesses
    ) and all(
        form == {"alternating": [1, 14], "kernel_dimension": 0, "symmetric": [0, None]}
        for form in forms
    )
    prime_ok = primes == EXPECTED_PRIMES and len(set(primes)) == 2
    caveat_ok = (
        all(item.get("levi_factors") is None for item in classifications)
        and all("no characteristic-zero lift is claimed" in item.get("field_scope", "") for item in classifications)
        and grid.get("levi_decomposition") == "NOT_COMPUTED_BY_THIS_EXPERIMENT"
        and grid.get("solvable_radical") == "NOT_COMPUTED_BY_THIS_EXPERIMENT"
        and grid.get("characteristic_zero_resolution", {}).get("resolved_elsewhere") is True
        and grid.get("characteristic_zero_resolution", {}).get("artifact") == "results/algebra_structure/char0_levi.json"
        and {"characteristic-zero Killing rank", "solvable radical", "Levi dimension", "Levi rank", "Levi type"}
        <= set(grid["field_provenance"]["not_computed"])
    )

    current_sanity = []
    current_ok = True
    for L in range(1, 7):
        actual = truncated_current_lower_central(L)
        expected = [3 * (L - k) for k in range(1, L)]
        passed = actual == expected
        current_ok &= passed
        current_sanity.append(
            {
                "claim_tag": "[COMPUTATION]",
                "truncation_order": L,
                "dimension": 3 * L,
                "evaluation_quotient_dimension": 3,
                "positive_u_ideal_dimension": 3 * (L - 1),
                "positive_u_ideal_lower_central_dimensions": actual,
                "expected_lower_central_dimensions": expected,
                "passed": passed,
            }
        )

    checks = [
        check_record(
            "fulltext_hashes_match_manifest",
            source_ok,
            "; ".join(f"{key}={record['computed_sha256']}" for key, record in source_data.items()),
        ),
        check_record(
            "abstract_only_sources_not_promoted",
            abstract_ok,
            ", ".join(f"{key}={value}" for key, value in abstract_status.items()),
        ),
        check_record(
            "two_distinct_recorded_primes",
            prime_ok,
            f"primes={primes}",
        ),
        check_record(
            "sp14_quotient_witness_at_both_primes",
            witness_ok,
            f"C7 dimensions={[witness.get('dimension') if witness else None for witness in witnesses]}; forms={forms}",
        ),
        check_record(
            "modular_data_not_promoted_to_characteristic_zero",
            caveat_ok,
            f"levi={grid.get('levi_decomposition')}; radical={grid.get('solvable_radical')}",
        ),
        check_record(
            "truncated_current_nilpotent_ideal_sanity",
            current_ok,
            f"orders=1..6; lower-central dimensions={[item['positive_u_ideal_lower_central_dimensions'] for item in current_sanity]}",
        ),
    ]

    payload = {
        "meta": {
            "status": "CONDITIONAL_THEOREM",
            "provenance": {
                "script": "experiments/e47_oa_quotient_sanity.py",
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "python": platform.python_version(),
                "arithmetic": "exact integers, exact JSON fields, and SHA-256 byte hashes; no floating point",
            },
            "scope": "External classification over C plus a conditional application to the 2x3 layer algebra; the characteristic-zero antecedent is not evaluated here.",
        },
        "data": {
            "source_audit": {
                "claim_tag": "[COMPUTATION]",
                "fulltext": source_data,
                "abstract_only": abstract_status,
            },
            "external_classification": {
                "claim_tag": "[EXTERNAL]",
                "base_field": "C",
                "closed_ideal_meaning": "I is closed iff OA/I has zero centre (Date-Roan definition, not a topology)",
                "off_fixed_point_factor": "(C[u]/u^L C[u]) tensor sl2",
                "fixed_point_factor": "OA/I_(t-1)^L or OA/I_(t+1)^L; solvable of dimension L+floor(L/2); only even L occurs in a closed ideal",
            },
            "levi_corollary": {
                "claim_tag": "[THEOREM]",
                "statement": "Every finite-dimensional semisimple quotient of OA, hence every Levi factor of an arbitrary finite-dimensional OA quotient, is a direct sum of sl2 copies.",
                "proof_dependency": "Date-Roan closed-ideal classification plus the radical argument in proofs/onsager_quotient_nogo.md",
            },
            "conditional_application": {
                "claim_tag": "[THEOREM, CONDITIONAL]",
                "antecedent": "The characteristic-zero semisimple quotient of the 2x3 layer algebra has a simple ideal not isomorphic to sl2.",
                "conclusion": "The 2x3 layer algebra is not isomorphic to any quotient of OA; no generating pair in it can satisfy the normalized Dolan-Grady presentation.",
                "antecedent_status_in_this_artifact": "UNRESOLVED; delegated characteristic-zero computation required",
            },
            "modular_consistency": {
                "claim_tag": "[COMPUTATION]",
                "primes": primes,
                "sector": "000",
                "field_type": "C7",
                "module_dimension": 14,
                "image_dimension": 105,
                "alternating_form_kernel_dimension": [form["kernel_dimension"] for form in forms],
                "scope": "Exact over the two recorded finite fields only; this does not prove the characteristic-zero antecedent.",
            },
            "truncated_current_sanity": current_sanity,
        },
        "checks": checks,
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    failed = [check["name"] for check in checks if not check["passed"]]
    if failed:
        print("FAIL: " + ", ".join(failed))
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
