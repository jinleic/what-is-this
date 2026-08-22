#!/usr/bin/env python3
"""Onsager-quotient corollary for the solved 3x3 layer algebra.

Composes two certified inputs into an unconditional no-go theorem:

1. [THEOREM, in-repo]  ``results/algebra_structure/char0_3x3.json`` — the exact
   characteristic-zero structure of the open-``3x3`` layer algebra
   (wave 9, clean-room verified by ``tests/test_algebra_3x3.py``):
   ``g = Qz + sl33 + sl48 + sl59 + sl3 + sl16 + sl30``, radical = centre = Qz,
   derived algebra semisimple of type ``A32+A47+A58+A2+A15+A29``.

2. [EXTERNAL, audited]  ``results/algebra_structure/oa_quotient.json`` — the
   Date--Roan closed-ideal classification (full-text sources with pinned
   SHA-256): every finite-dimensional semisimple quotient of the complex
   Onsager algebra is a direct sum of copies of ``sl2``; hence every Levi
   factor of an arbitrary finite-dimensional Onsager quotient is a direct sum
   of ``sl2`` copies.

Corollary proved here (see ``proofs/oa_3x3.md``): the complexified open-``3x3``
layer algebra is NOT a quotient of the Onsager algebra, and no generating pair
inside it satisfies the normalized Dolan--Grady relations.  The contradiction
is available through EVERY one of the six simple factors — none is ``sl2``
(smallest is ``sl3``).

This experiment verifies every antecedent from disk (recomputing source
hashes), checks the factor arithmetic exactly, and emits
``results/algebra_structure/oa_3x3.json``.
"""
from __future__ import annotations

import datetime
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STRUCTURE = ROOT / "results" / "algebra_structure" / "char0_3x3.json"
OA_AUDIT = ROOT / "results" / "algebra_structure" / "oa_quotient.json"
MANIFEST = ROOT / "sources" / "manifest.yaml"
RESULT = ROOT / "results" / "algebra_structure" / "oa_3x3.json"

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

EXPECTED_RANKS = [32, 47, 58, 2, 15, 29]
EXPECTED_DIMS = [1088, 2303, 3480, 8, 255, 899]
EXPECTED_TYPE = "A32 + A47 + A58 + A2 + A15 + A29"


def sha256_of(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            value.update(chunk)
    return value.hexdigest()


def main() -> int:
    checks: list[dict] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})
        print(f"{'PASS' if passed else 'FAIL'}: {name} ({detail})")

    # ---- antecedent 1: the 3x3 structure theorem, consumed with its own checks ----
    structure = json.loads(STRUCTURE.read_text(encoding="utf-8"))
    struct_checks_ok = bool(structure["checks"]) and all(
        c["passed"] for c in structure["checks"]
    )
    check(
        "structure_artifact_checks_all_passed",
        struct_checks_ok,
        f"{len(structure['checks'])} embedded checks in char0_3x3.json",
    )

    levi = structure["data"]["Levi_type"]
    radical = structure["data"]["solvable_radical"]
    derived = structure["data"]["derived_algebra"]
    check(
        "structure_claims_are_theorems",
        levi["claim_tag"] == "[THEOREM]"
        and radical["claim_tag"] == "[THEOREM]"
        and derived["claim_tag"] == "[THEOREM]"
        and derived["semisimple"] is True
        and radical["equals_center"] is True
        and radical["dimension_Q"] == 1,
        "Levi/radical/derived all [THEOREM]; derived semisimple; radical = centre, dim 1",
    )

    ranks = list(levi["factor_ranks"])
    dims = list(levi["factor_dimensions"])
    check(
        "levi_factor_data_matches_certified_values",
        ranks == EXPECTED_RANKS and dims == EXPECTED_DIMS and levi["type"] == EXPECTED_TYPE,
        f"ranks {ranks}, dims {dims}, type {levi['type']}",
    )

    # exact factor arithmetic: dim(sl_d) = d^2 - 1 with d = rank + 1, and d >= 3
    factor_rows = []
    all_non_sl2 = True
    arithmetic_ok = True
    for rank, dim in zip(ranks, dims):
        d = rank + 1
        arithmetic_ok &= dim == d * d - 1
        all_non_sl2 &= d >= 3
        factor_rows.append(
            {
                "absolute_type": f"A{rank}",
                "rational_form": f"sl({d},Q)",
                "dimension": dim,
                "is_sl2": d == 2,
            }
        )
    check(
        "factor_arithmetic_exact",
        arithmetic_ok and len(ranks) == 6,
        "all six satisfy dim = (rank+1)^2 - 1",
    )
    check(
        "every_factor_is_non_sl2",
        all_non_sl2,
        "smallest factor is sl(3,Q) (A2, dim 8); no factor has rank 1",
    )

    # ---- antecedent 2: the audited external classification ----
    oa = json.loads(OA_AUDIT.read_text(encoding="utf-8"))
    oa_checks_ok = bool(oa["checks"]) and all(c["passed"] for c in oa["checks"])
    levi_corollary = oa["data"]["levi_corollary"]
    check(
        "external_classification_artifact_intact",
        oa_checks_ok and levi_corollary["claim_tag"] == "[THEOREM]",
        "oa_quotient.json checks pass; levi_corollary tagged [THEOREM] over the audited sources",
    )
    stmt = levi_corollary["statement"]
    check(
        "levi_corollary_states_sl2_property",
        "direct sum of sl2 copies" in stmt and "Levi factor" in stmt,
        f"statement: {stmt[:120]}",
    )

    manifest_text = MANIFEST.read_text(encoding="utf-8")
    source_rows = {}
    sources_ok = True
    for key, (relpath, expected_hash) in EXPECTED_SOURCES.items():
        path = ROOT / relpath
        actual = sha256_of(path)
        ok = actual == expected_hash and key in manifest_text and expected_hash in manifest_text
        sources_ok &= ok
        source_rows[key] = {"path": relpath, "sha256": actual, "matches_pin": ok}
    check(
        "fulltext_sources_sha256_recomputed",
        sources_ok,
        "; ".join(f"{k}={v['sha256'][:16]}..." for k, v in source_rows.items()),
    )

    # ---- the corollary itself ----
    theorem = {
        "claim_tag": "[THEOREM]",
        "statement": (
            "The complexified open-3x3 layer algebra g_C is not isomorphic to any "
            "quotient of the Onsager algebra, and no pair of elements of g_C which "
            "generates g_C satisfies the normalized Dolan--Grady relations."
        ),
        "proof_chain": [
            "g_Q = Qz + s with s = sl33 + sl48 + sl59 + sl3 + sl16 + sl30 the derived "
            "algebra, semisimple, and radical(g) = centre(g) = Qz "
            "[THEOREM, char0_3x3.json; clean-room verified by tests/test_algebra_3x3.py]",
            "g surjects onto g/rad(g) = s; a semisimple algebra surjects onto each of "
            "its simple ideals (project along the complementary ideal); scalar "
            "extension to C is right exact",
            "hence g_C surjects onto sl(33,C) (equally: any of the six factors; the "
            "smallest, sl(3,C), already suffices) — a simple algebra of type A32, "
            "not sl2",
            "if g_C were an Onsager quotient, sl(33,C) would be a finite-dimensional "
            "semisimple quotient of the Onsager algebra, hence a direct sum of sl2 "
            "copies [EXTERNAL, Date--Roan, audited]; a simple such sum is sl2 itself; "
            "dim sl(33,C) = 1088 != 3 — contradiction",
            "a generating pair satisfying the normalized Dolan--Grady relations would "
            "induce a surjection from the universal Onsager presentation onto g_C, "
            "giving the same contradiction",
        ],
        "conditionality": (
            "unconditional given the in-repo structure theorem and the audited "
            "external Date--Roan classification; identical logical shape to the "
            "audited open-2x4 corollary (proofs/levi_images_2x4.md sec. 6)"
        ),
        "robustness": (
            "the contradiction is available through every one of the six simple "
            "factors independently; no single-factor certification is load-bearing"
        ),
    }
    check(
        "corollary_assembled",
        all_non_sl2 and struct_checks_ok and oa_checks_ok and sources_ok,
        "all antecedents verified from disk in this run",
    )

    payload = {
        "provenance": {
            "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "script": "experiments/e100_oa_3x3.py",
            "inputs": [
                "results/algebra_structure/char0_3x3.json",
                "results/algebra_structure/oa_quotient.json",
                "sources/manifest.yaml",
                "sources/fulltext/date_roan2000_onsager.pdf",
                "sources/fulltext/roan1991_onsager.pdf",
            ],
            "method": (
                "logical composition of the wave-9 3x3 structure theorem with the "
                "audited Date--Roan Levi-factor corollary; all antecedent artifacts "
                "and source hashes re-verified from disk at run time"
            ),
        },
        "data": {
            "onsager_quotient_no_go_3x3": theorem,
            "levi_factors": factor_rows,
            "antecedent_structure": {
                "dimension_Q": structure["data"]["dimension_Q"]["value"],
                "levi_type": levi["type"],
                "radical_equals_centre": True,
                "artifact": "results/algebra_structure/char0_3x3.json",
            },
            "external_sources": source_rows,
        },
        "checks": checks,
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(payload, indent=1, sort_keys=True), encoding="utf-8")
    print(f"WROTE {RESULT.relative_to(ROOT)}")

    if not all(c["passed"] for c in checks):
        return 1
    print("PASS e100_oa_3x3: the open-3x3 layer algebra is not an Onsager-algebra quotient")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
