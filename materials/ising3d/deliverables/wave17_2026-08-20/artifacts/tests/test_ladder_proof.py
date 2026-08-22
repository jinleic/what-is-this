"""Standalone checks for the independent two-leg-ladder bracket experiment."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from experiments.e28_ladder_proof import (  # noqa: E402
    PRIME_31,
    SECOND_PRIME_31,
    ModularEchelon,
    boundary_lex_key,
    chain_bonds,
    closure_profile,
    depth_certificate,
    ladder_bonds,
    recursive_candidate_words,
    rung_pauli_string,
    vector_for_word,
)

RESULT = ROOT / "results" / "algebra_growth" / "ladder_proof.json"


def main() -> int:
    try:
        expected = {2: 6, 3: 16, 4: 38, 5: 85, 6: 197}

        # Rebuild every base certificate from actual bracket words.  No stored row and
        # no e20 basis builder participates in these ranks.
        fresh = {
            length: depth_certificate(length, include_full_supports=False)
            for length in range(2, 7)
        }
        for length, dimension in expected.items():
            record = fresh[length]
            assert record["max_depth"] == 2 * length, record
            assert record["D_2L"] == {
                str(PRIME_31): dimension,
                str(SECOND_PRIME_31): dimension,
            }, record
            assert len(record["maximal_independent_words_mod_prime"]) == dimension, record
            assert all(
                row["depth"] <= 2 * length
                for row in record["maximal_independent_words_mod_prime"]
            ), record
            assert record["recursive_candidate"]["rank"] == {
                str(PRIME_31): 1 << length,
                str(SECOND_PRIME_31): 1 << length,
            }, record

        # Explicit finite refutation of the fixed boundary-lex pivot proposal.
        assert fresh[2]["boundary_lex_distinct_leaders"] == 5
        assert fresh[3]["boundary_lex_distinct_leaders"] == 9
        assert fresh[4]["boundary_lex_distinct_leaders"] == 15 < 16
        collision = fresh[4]["recursive_candidate"]["first_boundary_lex_collision"]
        assert collision == {
            "first_label": "0101",
            "second_label": "0110",
            "leader": 33928,
            "leader_pauli": "II II ZI XY",
            "first_word": "BBABBA",
            "second_word": "ABBBBA",
        }, collision

        # Independently expand the colliding pair and verify the claimed common leader.
        n, bonds = ladder_bonds(4)
        family = dict(recursive_candidate_words(4))
        first = vector_for_word(family["0101"], n, bonds)
        second = vector_for_word(family["0110"], n, bonds)
        first_leader = max(first, key=lambda pauli: boundary_lex_key(pauli, 4))
        second_leader = max(second, key=lambda pauli: boundary_lex_key(pauli, 4))
        assert first_leader == second_leader == 33928
        assert rung_pauli_string(first_leader, 4) == "II II ZI XY"
        pair_basis = ModularEchelon(PRIME_31)
        pair_basis.add(first)
        pair_basis.add(second)
        assert pair_basis.rank == 2  # a pivot collision, not vector proportionality

        # Required controls, recomputed with the same independent full-Pauli engine.
        for sites in range(2, 7):
            _, bonds = chain_bonds(sites)
            assert closure_profile(sites, bonds, PRIME_31)[-1] == sites * sites
        for sites in range(3, 7):
            _, bonds = chain_bonds(sites, periodic=True)
            assert closure_profile(sites, bonds, PRIME_31)[-1] == 3 * sites - 1
        n, bonds = ladder_bonds(2)
        assert closure_profile(n, bonds, PRIME_31)[-1] == 11

        with RESULT.open(encoding="utf-8") as handle:
            artifact = json.load(handle)
        assert artifact["provenance"]["script"] == "experiments/e28_ladder_proof.py"
        assert artifact["data"]["pivot_rule_verified_lengths"] == [2, 3]
        assert artifact["data"]["pivot_rule_counterexample_length"] == 4
        assert [
            row["D_2L"][str(PRIME_31)] for row in artifact["data"]["ladders"]
        ] == list(expected.values())
        assert all(artifact["checks"]) and all(
            check["passed"] for check in artifact["checks"]
        )
        print("PASS")
        return 0
    except Exception as exc:
        print(f"FAIL: {type(exc).__name__}: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
