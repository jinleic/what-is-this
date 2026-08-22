"""Standalone checks for the conditional Onsager-quotient obstruction."""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "e47_oa_quotient_sanity.py"
RESULT = ROOT / "results" / "algebra_structure" / "oa_quotient.json"
STRUCTURE = ROOT / "results" / "algebra_structure" / "structure.json"
MANIFEST = ROOT / "sources" / "manifest.yaml"
EXPECTED_SOURCES = {
    "date_roan2000_onsager": (
        "sources/fulltext/date_roan2000_onsager.pdf",
        "c1334b58ed30c77b1592e628d8bf26c40d52ef63b9226ecad1ddec0e96771b06",
    ),
    "roan1991_onsager": (
        "sources/fulltext/roan1991_onsager.pdf",
        "82988421997ae8fb388c79c4f46381def00ee15ecaaed04922bd7ab94c352313",
    ),
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            value.update(chunk)
    return value.hexdigest()


def scalar(text: str) -> Any:
    text = text.strip()
    if text == "null":
        return None
    if text.startswith('"'):
        return json.loads(text)
    return text


def manifest_scalars(text: str, key: str) -> dict[str, Any]:
    marker = f"  - key: {key}"
    start = text.find(marker)
    if start < 0:
        return {}
    next_start = text.find("\n  - key: ", start + len(marker))
    block = text[start : len(text) if next_start < 0 else next_start]
    fields: dict[str, Any] = {"key": key}
    for line in block.splitlines()[1:]:
        match = re.match(r"^\s{4}([a-z][a-z0-9_]*):\s*(.*?)\s*$", line)
        if match:
            fields[match.group(1)] = scalar(match.group(2))
    return fields


def lower_central_dimensions(L: int) -> list[int]:
    """Independently bracket the positive-u ideal of C[u]/u^L tensor sl2."""
    bracket_basis = {
        (0, 1): 2,
        (1, 0): 2,
        (2, 0): 0,
        (0, 2): 0,
        (2, 1): 1,
        (1, 2): 1,
    }
    ideal = {(degree, basis) for degree in range(1, L) for basis in range(3)}
    term = set(ideal)
    dimensions: list[int] = []
    while term:
        dimensions.append(len(term))
        term = {
            (degree_x + degree_y, bracket_basis[(basis_x, basis_y)])
            for degree_x, basis_x in ideal
            for degree_y, basis_y in term
            if degree_x + degree_y < L and (basis_x, basis_y) in bracket_basis
        }
    return dimensions


def main() -> int:
    failures: list[str] = []

    def check(name: str, condition: bool, detail: str) -> None:
        if condition:
            print(f"PASS: {name}")
        else:
            print(f"FAIL: {name}: {detail}")
            failures.append(name)

    run = subprocess.run(
        [sys.executable, str(SCRIPT)],
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )
    check(
        "experiment_runs",
        run.returncode == 0 and run.stdout.strip().splitlines()[-1:] == ["PASS"],
        f"returncode={run.returncode}, stdout={run.stdout!r}, stderr={run.stderr!r}",
    )
    if run.returncode != 0:
        return 1

    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    check(
        "result_envelope",
        set(artifact) == {"meta", "data", "checks"}
        and artifact["meta"].get("status") == "CONDITIONAL_THEOREM",
        f"keys={sorted(artifact)}, meta={artifact.get('meta')}",
    )
    checks_well_formed = bool(artifact["checks"]) and all(
        set(item) == {"name", "passed", "detail"}
        and isinstance(item["name"], str)
        and isinstance(item["passed"], bool)
        and isinstance(item["detail"], str)
        and item["passed"]
        for item in artifact["checks"]
    )
    check("result_checks_pass", checks_well_formed, repr(artifact["checks"]))

    manifest_text = MANIFEST.read_text(encoding="utf-8")
    source_ok = True
    source_details = []
    for key, (relative_path, expected_hash) in EXPECTED_SOURCES.items():
        fields = manifest_scalars(manifest_text, key)
        path = ROOT / relative_path
        actual_hash = digest(path) if path.is_file() else None
        pdf_magic = path.read_bytes()[:5] if path.is_file() else b""
        record_ok = (
            fields.get("obtained") == "full"
            and fields.get("local_path") == relative_path
            and fields.get("sha256") == expected_hash
            and actual_hash == expected_hash
            and pdf_magic == b"%PDF-"
        )
        source_ok &= record_ok
        source_details.append(f"{key}: manifest={fields.get('sha256')}, actual={actual_hash}")
    check("fulltext_manifest_hashes", source_ok, "; ".join(source_details))

    abstract_fields = [
        manifest_scalars(manifest_text, key)
        for key in ("davies1990_onsager", "davies1991_onsager")
    ]
    abstract_ok = all(
        fields.get("obtained") == "abstract_only"
        and fields.get("local_path") is None
        and fields.get("sha256") is None
        for fields in abstract_fields
    )
    check("abstract_access_is_honest", abstract_ok, repr(abstract_fields))

    structure = json.loads(STRUCTURE.read_text(encoding="utf-8"))
    grid = structure["data"]["grid_2x3"]
    classifications = [grid["classification_prime_1"], grid["classification_prime_2"]]
    primes = [classification["prime"] for classification in classifications]
    prime_ok = primes == [2147483647, 2147483629] and len(set(primes)) == 2
    check("two_prime_identity", prime_ok, f"primes={primes}")

    witness_ok = True
    witness_details = []
    for classification in classifications:
        c7 = [
            item
            for item in classification["certified_quotient_images"]
            if item.get("field_type") == "C7"
        ]
        form = classification["linkage_witnesses"]["invariant_form_witnesses"]["000"]
        local_ok = (
            len(c7) == 1
            and c7[0].get("dimension") == 105
            and c7[0].get("sector") == "000"
            and form == {"alternating": [1, 14], "kernel_dimension": 0, "symmetric": [0, None]}
        )
        witness_ok &= local_ok
        witness_details.append(f"p={classification['prime']}: C7={c7}, form={form}")
    check("sp14_witness_both_primes", witness_ok, "; ".join(witness_details))

    no_promotion = (
        all(classification.get("levi_factors") is None for classification in classifications)
        and all(
            "no characteristic-zero lift is claimed" in classification.get("field_scope", "")
            for classification in classifications
        )
        and grid.get("levi_decomposition") == "NOT_COMPUTED_BY_THIS_EXPERIMENT"
        and grid.get("solvable_radical") == "NOT_COMPUTED_BY_THIS_EXPERIMENT"
        and grid.get("characteristic_zero_resolution", {}).get("resolved_elsewhere") is True
        and "lower bounds over Q" in structure["data"].get("scope", "")
        and artifact["data"]["conditional_application"]["antecedent_status_in_this_artifact"].startswith("UNRESOLVED")
    )
    check("no_mod_p_to_char0_promotion", no_promotion, repr(grid.get("field_provenance")))

    derived_rows = []
    current_ok = True
    recorded = artifact["data"]["truncated_current_sanity"]
    for L in range(1, 7):
        actual = lower_central_dimensions(L)
        expected = [3 * (L - degree) for degree in range(1, L)]
        row = recorded[L - 1]
        local_ok = (
            actual == expected
            and row["truncation_order"] == L
            and row["positive_u_ideal_lower_central_dimensions"] == actual
            and row["evaluation_quotient_dimension"] == 3
        )
        current_ok &= local_ok
        derived_rows.append((L, actual))
    check("exact_truncated_current_sanity", current_ok, repr(derived_rows))

    if failures:
        print("FAIL")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
