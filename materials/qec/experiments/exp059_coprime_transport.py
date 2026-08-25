"""EXP-059: exact CRT transport between coprime odd BB lattices.

If gcd(ell,m)=1 and N=ell*m, the map pi -> (x,y) identifies
F_2[pi]/(pi^N-1) with F_2[x,y]/(x^ell-1,y^m-1).  Consequently every
coprime factorization of the same N has exactly the same BB search space up to
an explicit coordinate permutation.  This removes duplicated lattice screens.

The first payoff is N=105 (n=210): (15,7), (21,5), and (35,3) are one exact
search problem, not three independent frontiers.

Run:
  python experiments/exp059_coprime_transport.py run --N 105
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.codes.bicycle import poly_matrix  # noqa: E402
from qec_research.gf2.linalg import rank_np  # noqa: E402

SCHEMA = "exp059-coprime-transport-v1"
OUT = ROOT / "results" / "processed" / "exp059_coprime_transport.json"


def matrix_sha256(matrix: np.ndarray) -> str:
    array = np.ascontiguousarray(np.asarray(matrix, dtype=np.uint8) & 1)
    return hashlib.sha256(array.tobytes()).hexdigest()


def _check_coprime(ell: int, m: int) -> None:
    if ell <= 0 or m <= 0 or math.gcd(ell, m) != 1:
        raise ValueError(f"lattice ({ell},{m}) is not a positive coprime pair")


def crt_exponent(i: int, j: int, ell: int, m: int) -> int:
    """Unique e mod ell*m with e=i mod ell and e=j mod m."""
    _check_coprime(ell, m)
    N = ell * m
    i, j = int(i) % ell, int(j) % m
    return int(
        (i * m * pow(m, -1, ell) + j * ell * pow(ell, -1, m)) % N
    )


def exponents_to_terms(
    exponents: Iterable[int], ell: int, m: int
) -> list[tuple[int, int]]:
    _check_coprime(ell, m)
    N = ell * m
    return [(int(e) % ell, int(e) % m) for e in exponents if 0 <= int(e) < N]


def terms_to_exponents(
    terms: Iterable[tuple[int, int]], ell: int, m: int
) -> list[int]:
    return [crt_exponent(i, j, ell, m) for i, j in terms]


def transport_terms(
    terms: Iterable[tuple[int, int]],
    source_ell: int,
    source_m: int,
    target_ell: int,
    target_m: int,
) -> list[tuple[int, int]]:
    if source_ell * source_m != target_ell * target_m:
        raise ValueError("transport requires equal block size")
    exponents = terms_to_exponents(terms, source_ell, source_m)
    return exponents_to_terms(exponents, target_ell, target_m)


def coordinate_transport(
    source_ell: int,
    source_m: int,
    target_ell: int,
    target_m: int,
) -> np.ndarray:
    """mapping[src_index] = dst_index for the CRT-labelled coordinates."""
    if source_ell * source_m != target_ell * target_m:
        raise ValueError("transport requires equal block size")
    _check_coprime(source_ell, source_m)
    _check_coprime(target_ell, target_m)
    N = source_ell * source_m
    mapping = np.empty(N, dtype=np.int64)
    for i in range(source_ell):
        for j in range(source_m):
            e = crt_exponent(i, j, source_ell, source_m)
            target_i, target_j = e % target_ell, e % target_m
            mapping[i * source_m + j] = target_i * target_m + target_j
    if not np.array_equal(np.sort(mapping), np.arange(N)):
        raise RuntimeError("CRT coordinate transport is not a permutation")
    return mapping


def permutation_matrix(mapping: np.ndarray) -> np.ndarray:
    mapping = np.asarray(mapping, dtype=np.int64)
    N = len(mapping)
    matrix = np.zeros((N, N), dtype=np.uint8)
    matrix[mapping, np.arange(N)] = 1
    return matrix


def bb_matrices(
    ell: int,
    m: int,
    A_terms: Iterable[tuple[int, int]],
    B_terms: Iterable[tuple[int, int]],
) -> tuple[np.ndarray, np.ndarray]:
    A = poly_matrix(ell, m, list(A_terms))
    B = poly_matrix(ell, m, list(B_terms))
    return (
        np.hstack([A, B]).astype(np.uint8),
        np.hstack([B.T, A.T]).astype(np.uint8),
    )


def transport_certificate(
    source_ell: int,
    source_m: int,
    target_ell: int,
    target_m: int,
    A_terms: Iterable[tuple[int, int]],
    B_terms: Iterable[tuple[int, int]],
) -> dict[str, Any]:
    A_terms, B_terms = list(A_terms), list(B_terms)
    target_A = transport_terms(
        A_terms, source_ell, source_m, target_ell, target_m
    )
    target_B = transport_terms(
        B_terms, source_ell, source_m, target_ell, target_m
    )
    source_HX, source_HZ = bb_matrices(
        source_ell, source_m, A_terms, B_terms
    )
    target_HX, target_HZ = bb_matrices(
        target_ell, target_m, target_A, target_B
    )
    mapping = coordinate_transport(
        source_ell, source_m, target_ell, target_m
    )
    P = permutation_matrix(mapping)
    zero = np.zeros_like(P)
    Q = np.block([[P, zero], [zero, P]])
    expected_HX = (P @ source_HX @ Q.T) % 2
    expected_HZ = (P @ source_HZ @ Q.T) % 2
    hx_exact = bool(np.array_equal(target_HX, expected_HX))
    hz_exact = bool(np.array_equal(target_HZ, expected_HZ))
    source_k = int(source_HX.shape[1] - rank_np(source_HX) - rank_np(source_HZ))
    target_k = int(target_HX.shape[1] - rank_np(target_HX) - rank_np(target_HZ))
    valid = bool(
        hx_exact and hz_exact and source_k == target_k
        and np.array_equal(np.sort(mapping), np.arange(len(mapping)))
    )
    if not valid:
        raise RuntimeError("BB CRT transport failed exact matrix gates")
    return {
        "schema": SCHEMA,
        "source_lattice": [source_ell, source_m],
        "target_lattice": [target_ell, target_m],
        "source_n": int(source_HX.shape[1]),
        "target_n": int(target_HX.shape[1]),
        "source_A": [list(term) for term in A_terms],
        "source_B": [list(term) for term in B_terms],
        "target_A": [list(term) for term in target_A],
        "target_B": [list(term) for term in target_B],
        "mapping_sha256": matrix_sha256(mapping[None, :]),
        "HX_exact": hx_exact,
        "HZ_exact": hz_exact,
        "source_k": source_k,
        "target_k": target_k,
        "valid": valid,
    }


def coprime_factorizations(N: int) -> list[tuple[int, int]]:
    if N <= 0 or N % 2 == 0:
        raise ValueError("N must be positive and odd")
    out = []
    for m in range(1, math.isqrt(N) + 1):
        if N % m:
            continue
        ell = N // m
        if math.gcd(ell, m) == 1:
            out.append((ell, m))
    return sorted(out, key=lambda pair: pair[0])


def _all_monomials_transport(
    source: tuple[int, int], target: tuple[int, int]
) -> bool:
    N = source[0] * source[1]
    mapping = coordinate_transport(*source, *target)
    P = permutation_matrix(mapping)
    for exponent in range(N):
        source_term = exponents_to_terms([exponent], *source)
        target_term = exponents_to_terms([exponent], *target)
        source_matrix = poly_matrix(*source, source_term)
        target_matrix = poly_matrix(*target, target_term)
        if not np.array_equal(target_matrix, (P @ source_matrix @ P.T) % 2):
            return False
    return True


def coprime_family_certificate(N: int) -> dict[str, Any]:
    factorizations = coprime_factorizations(N)
    screen = [pair for pair in factorizations if pair[1] >= 3]
    pairwise = []
    for source, target in combinations(screen, 2):
        monomials_exact = _all_monomials_transport(source, target)
        mapping = coordinate_transport(*source, *target)
        inverse = coordinate_transport(*target, *source)
        inverse_exact = bool(np.array_equal(inverse[mapping], np.arange(N)))
        pairwise.append(
            {
                "source": list(source),
                "target": list(target),
                "all_monomials_exact": monomials_exact,
                "inverse_exact": inverse_exact,
                "mapping_sha256": matrix_sha256(mapping[None, :]),
                "valid": bool(monomials_exact and inverse_exact),
            }
        )
    all_valid = all(item["valid"] for item in pairwise)
    if not all_valid:
        raise RuntimeError("coprime family transport proof failed")
    return {
        "schema": SCHEMA,
        "N": N,
        "n": 2 * N,
        "factorizations": [list(pair) for pair in factorizations],
        "screen_lattices": [list(pair) for pair in screen],
        "representative": list(screen[0]) if screen else None,
        "pairwise": pairwise,
        "all_pairwise_transports_valid": all_valid,
        "theorem": (
            "For gcd(ell,m)=1, pi -> (x,y) is a ring isomorphism from "
            "F2[pi]/(pi^N-1); the stored coordinate permutation conjugates "
            "every monomial matrix, hence every BB HX/HZ pair and its distance."
        ),
    }


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _map_block_support(
    support: list[int] | None, mapping: np.ndarray, *, two_blocks: bool
) -> list[int] | None:
    if support is None:
        return None
    N = len(mapping)
    output = []
    for index in support:
        index = int(index)
        if two_blocks and index >= N:
            output.append(N + int(mapping[index - N]))
        else:
            output.append(int(mapping[index]))
    return sorted(output)


def _transport_record(
    record: dict[str, Any],
    source: tuple[int, int],
    target: tuple[int, int],
    mapping: np.ndarray,
) -> dict[str, Any]:
    transported = json.loads(json.dumps(record))
    transported["A"] = [
        list(term) for term in transport_terms(record["A"], *source, *target)
    ]
    transported["B"] = [
        list(term) for term in transport_terms(record["B"], *source, *target)
    ]
    transported["ell"], transported["m"] = target
    transported["n"] = 2 * target[0] * target[1]
    transported["ceiling_witness_support"] = _map_block_support(
        record.get("ceiling_witness_support"), mapping, two_blocks=False
    )
    transported["witness_support"] = _map_block_support(
        record.get("witness_support"), mapping, two_blocks=True
    )
    if isinstance(transported.get("cdcl"), dict):
        transported["cdcl"]["witness_support"] = _map_block_support(
            record["cdcl"].get("witness_support"), mapping, two_blocks=True
        )
        if transported["cdcl"].get("cnf_sha256") is not None:
            transported["cdcl"]["source_cnf_sha256"] = transported["cdcl"][
                "cnf_sha256"
            ]
            transported["cdcl"][
                "cnf_sha256_is_source_coordinate_order"
            ] = True
    transported["transported_from"] = {
        "lattice": list(source),
        "mapping_sha256": matrix_sha256(mapping[None, :]),
        "distance_preserved_by_exact_permutation": True,
    }
    return transported


def _logical_witness_valid(record: dict[str, Any]) -> bool:
    support = record.get("witness_support")
    if support is None:
        return False
    HX, HZ = bb_matrices(
        int(record["ell"]), int(record["m"]), record["A"], record["B"]
    )
    vector = np.zeros(int(record["n"]), dtype=np.uint8)
    vector[np.asarray(support, dtype=int)] = 1
    return bool(
        not np.any(HX @ vector % 2)
        and rank_np(np.vstack([HZ, vector])) == rank_np(HZ) + 1
        and int(vector.sum()) == int(record["witness_bound"])
    )

def _ceiling_witness_valid(record: dict[str, Any]) -> bool:
    support = record.get("ceiling_witness_support")
    if support is None:
        return False
    HX, HZ = bb_matrices(
        int(record["ell"]), int(record["m"]), record["A"], record["B"]
    )
    vector = np.zeros(int(record["n"]), dtype=np.uint8)
    vector[np.asarray(support, dtype=int)] = 1
    return bool(
        not np.any(HX @ vector % 2)
        and rank_np(np.vstack([HZ, vector])) == rank_np(HZ) + 1
        and int(vector.sum()) == int(record["ceiling"])
    )


def transport_screen_shard(
    source_path: Path,
    target_ell: int,
    target_m: int,
    output_path: Path,
) -> dict[str, Any]:
    """Transport one solved screen shard to an exactly equivalent lattice."""
    source_payload = json.loads(source_path.read_text(encoding="utf-8"))
    source = (int(source_payload["ell"]), int(source_payload["m"]))
    target = (int(target_ell), int(target_m))
    if source[0] * source[1] != target[0] * target[1]:
        raise ValueError("screen transport requires equal block size")
    mapping = coordinate_transport(*source, *target)
    records = [
        _transport_record(record, source, target, mapping)
        for record in source_payload["records"]
    ]
    witnesses = [
        record
        for record in records
        if record["verdict"]
        in {"dominated_by_witness", "dominated_by_cdcl_witness"}
    ]
    ceilings = [
        record for record in records
        if record["verdict"] == "dominated_by_ceiling"
    ]
    witness_valid = bool(
        all(_logical_witness_valid(record) for record in witnesses)
        and all(_ceiling_witness_valid(record) for record in ceilings)
    )
    if not witness_valid:
        raise RuntimeError("transported physical witness failed target checks")
    output = {
        **source_payload,
        "ell": target[0],
        "m": target[1],
        "n": 2 * target[0] * target[1],
        "records": records,
        "survivors": [
            record for record in records if record["verdict"] == "survivor"
        ],
        "no_reference": [
            record for record in records if record["verdict"] == "no_reference"
        ],
        "undecided": [
            record for record in records if record["verdict"] == "undecided"
        ],
        "transport": {
            "schema": SCHEMA,
            "source_lattice": list(source),
            "target_lattice": list(target),
            "source_shard_sha256": file_sha256(source_path),
            "theorem_certificate": str(OUT.relative_to(ROOT)),
            "theorem_certificate_sha256": file_sha256(OUT),
            "mapping_sha256": matrix_sha256(mapping[None, :]),
            "records_transported": len(records),
            "candidates_after_symmetry_preserved":
                source_payload["candidates_after_symmetry"],
            "orbit_total_preserved": source_payload["orbit_total"],
            "physical_witnesses_rechecked": len(witnesses) + len(ceilings),
            "all_physical_witnesses_valid": witness_valid,
            "valid": witness_valid,
        },
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    temporary.write_text(json.dumps(output, indent=1, sort_keys=True) + "\n")
    temporary.replace(output_path)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("run", "transport"))
    parser.add_argument("--N", type=int, default=105)
    parser.add_argument(
        "--source-shard",
        default="results/partial_runs/exp055_screen/15x7.json",
    )
    parser.add_argument("--targets", default="21x5,35x3")
    args = parser.parse_args()
    if args.command == "run":
        payload = coprime_family_certificate(args.N)
        OUT.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
        print(json.dumps({
            "N": payload["N"],
            "screen_lattices": payload["screen_lattices"],
            "representative": payload["representative"],
            "all_pairwise_transports_valid":
                payload["all_pairwise_transports_valid"],
        }, indent=1))
        return 0
    outputs = []
    for token in args.targets.split(","):
        ell, m = (int(value) for value in token.split("x"))
        path = ROOT / "results" / "partial_runs" / "exp055_screen" / f"{ell}x{m}.json"
        outputs.append(
            transport_screen_shard(
                ROOT / args.source_shard, ell, m, path
            )
        )
    print(json.dumps([
        {
            "lattice": [payload["ell"], payload["m"]],
            "records": len(payload["records"]),
            "transport_valid": payload["transport"]["valid"],
        }
        for payload in outputs
    ], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
