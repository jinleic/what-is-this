"""EXP-056: reciprocal-pole orbit certificates for odd-lattice BB distance.

The exact bounded question for a CSS BB code is reduced from one sector per
logical basis vector to one sector per translation-orbit generator of the
opposite logical quotient.  On an odd lattice the common pole ideal is
principal, so the two physical blocks require two sectors total.

For a Z candidate ``z`` and an X-logical orbit generator ``u``:

    exists nontrivial z of wt <= c
      iff exists block b and translation t with <z, t u_b> = 1
      iff exists block b and z' with <z', u_b> = 1 and wt(z') <= c.

The second equivalence translates ``z`` back.  A fixed-sector query must not use
the full translation origin anchor; it may only use the stabilizer subgroup of
its detector.  The implementation constructs that subgroup and independently
checks orbit coverage before accepting any UNSAT result.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.codes.bicycle import monomial_matrix, poly_matrix  # noqa: E402
from qec_research.distance.sat_decide import (  # noqa: E402
    DecisionInstance,
    cnf_sha256,
    decision_cnf_digest,
    decide_weight_bounded,
    group_weight,
    verify_witness_two_paths,
)
from qec_research.gf2.linalg import (  # noqa: E402
    nullspace_np,
    rank_bitset,
    rank_np,
    rref_np,
    rows_to_bitsets,
    solve_bitset,
)

SCHEMA = "exp056-odd-orbit-distance-v1"
SECTOR_SCHEMA = "exp056-odd-orbit-sector-v1"
STATE_DIR = ROOT / "results" / "partial_runs" / "exp056_odd_distance"
CERTIFICATE = ROOT / "results" / "certificates" / "exp056_wm_162_8_14_distance.json"
PARTIAL = ROOT / "results" / "partial_runs" / "exp056_wm_162_8_14_distance.json"
CLASS_STATE_DIR = STATE_DIR / "classes"
COUNTEREXAMPLE = (
    ROOT / "results" / "certificates" / "exp056_wm_162_8_14_counterexample.json"
)
SOLVER_NAME = "cadical195"

TARGET: dict[str, Any] = {
    "name": "wm-3x27-162-8-14",
    "source": "Wang--Mueller arXiv:2408.10001v4 Table 1",
    "source_distance_provenance": "BP-OSD distance_upperbound estimate",
    "ell": 3,
    "m": 27,
    "A": [(0, 0), (0, 10), (0, 14)],
    "B": [(0, 12), (1, 0), (2, 0)],
    "expected_k": 8,
    "expected_d": 14,
    "lower_cap": 13,
    # Deterministic reduced-pole representative; independently rechecked below.
    "witness_support": [
        0, 1, 4, 14, 84, 87, 94, 106, 111, 126, 133, 138, 153, 160,
    ],
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def matrix_sha256(matrix: np.ndarray) -> str:
    arr = np.ascontiguousarray(np.asarray(matrix, dtype=np.uint8) & 1)
    return sha256_bytes(arr.tobytes())


def canonical_json_sha256(value: Any) -> str:
    return sha256_bytes(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    )


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    tmp.replace(path)


def _independent_rows(matrix: np.ndarray) -> np.ndarray:
    rr, _ = rref_np(np.asarray(matrix, dtype=np.uint8) & 1)
    return rr[: rank_np(rr)]


def _intersection(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Basis of the intersection of two GF(2) row spaces."""
    left = _independent_rows(left)
    right = _independent_rows(right)
    if not len(left) or not len(right):
        return np.zeros((0, left.shape[1]), dtype=np.uint8)
    relations = nullspace_np(np.vstack([left, right]).T)
    rows = (relations[:, : len(left)] @ left) % 2
    return _independent_rows(rows)


def _span_nonzero(basis: np.ndarray) -> np.ndarray:
    dim = int(basis.shape[0])
    if dim == 0:
        return np.zeros((0, basis.shape[1]), dtype=np.uint8)
    if dim > 16:
        raise ValueError(f"pole dimension {dim} exceeds exact generator-search cap 16")
    coefficients = (
        (np.arange(1, 1 << dim, dtype=np.uint32)[:, None] >> np.arange(dim)) & 1
    ).astype(np.uint8)
    return (coefficients @ basis) % 2


def bar_permutation(ell: int, m: int) -> np.ndarray:
    idx = np.arange(ell * m).reshape(ell, m)
    return idx[(-np.arange(ell)) % ell][:, (-np.arange(m)) % m].reshape(-1)


def build_problem(record: dict[str, Any]) -> dict[str, Any]:
    ell, m = int(record["ell"]), int(record["m"])
    if ell % 2 == 0 or m % 2 == 0:
        raise ValueError("EXP-056's principal-pole proof requires an odd x odd lattice")
    A_terms = [tuple(map(int, term)) for term in record["A"]]
    B_terms = [tuple(map(int, term)) for term in record["B"]]
    A = poly_matrix(ell, m, A_terms)
    B = poly_matrix(ell, m, B_terms)
    HX = np.hstack([A, B]).astype(np.uint8)
    HZ = np.hstack([B.T, A.T]).astype(np.uint8)
    n = int(HX.shape[1])
    k = int(n - rank_np(HX) - rank_np(HZ))
    if np.any(HX @ HZ.T % 2):
        raise RuntimeError("rebuilt CSS checks do not commute")
    expected_k = record.get("expected_k")
    if expected_k is not None and k != int(expected_k):
        raise RuntimeError(f"rebuilt k={k} disagrees with expected k={expected_k}")
    return {
        "record": record,
        "ell": ell,
        "m": m,
        "block": ell * m,
        "n": n,
        "k": k,
        "A_terms": A_terms,
        "B_terms": B_terms,
        "A": A,
        "B": B,
        "HX": HX,
        "HZ": HZ,
    }


def kernel_weight_parity_proof(
    problem: dict[str, Any], requested_cap: int
) -> dict[str, Any]:
    """If every H_X column has odd degree, every vector in ker(H_X) is even."""
    column_degrees = np.asarray(problem["HX"], dtype=np.uint8).sum(axis=0)
    all_odd = bool(np.all(column_degrees % 2 == 1))
    effective = int(requested_cap) - (int(requested_cap) % 2) if all_odd else int(requested_cap)
    return {
        "schema": "exp056-kernel-weight-parity-v1",
        "valid": all_odd,
        "all_column_degrees_odd": all_odd,
        "column_degree_min": int(column_degrees.min()),
        "column_degree_max": int(column_degrees.max()),
        "requested_cap": int(requested_cap),
        "effective_even_cap": effective,
        "argument": (
            "Summing every H_X parity equation gives "
            "sum_j(deg_j mod 2)v_j = wt(v) mod 2 = 0."
        ),
    }


def _translate_block(vector: np.ndarray, ell: int, m: int, a: int, b: int) -> np.ndarray:
    return (np.asarray(vector, dtype=np.uint8) @ monomial_matrix(ell, m, a, b)) % 2


def _translate_two_blocks(
    vector: np.ndarray, ell: int, m: int, a: int, b: int
) -> np.ndarray:
    block = ell * m
    shift = monomial_matrix(ell, m, a, b)
    vector = np.asarray(vector, dtype=np.uint8)
    return np.concatenate([vector[:block] @ shift, vector[block:] @ shift]) % 2


def _embed(seed: np.ndarray, block_index: int) -> np.ndarray:
    zero = np.zeros_like(seed)
    return np.concatenate([seed, zero] if block_index == 0 else [zero, seed])


def _quotient_rank(rows: np.ndarray, modulo: np.ndarray) -> tuple[int, int]:
    ncols = int(modulo.shape[1])
    base_np = rank_np(modulo)
    numpy_rank = rank_np(np.vstack([modulo, rows])) - base_np
    base_bits = rank_bitset(rows_to_bitsets(modulo), ncols)
    bitset_rank = rank_bitset(rows_to_bitsets(np.vstack([modulo, rows])), ncols) - base_bits
    return int(numpy_rank), int(bitset_rank)


def _cyclic_generator(pole: np.ndarray, ell: int, m: int) -> tuple[np.ndarray, np.ndarray]:
    """Choose a deterministic vector whose translation orbit spans ``pole``."""
    dimension = int(pole.shape[0])
    candidates: list[tuple[tuple[Any, ...], np.ndarray, np.ndarray]] = []
    for seed in _span_nonzero(pole):
        orbit = np.asarray(
            [
                _translate_block(seed, ell, m, a, b)
                for a in range(ell)
                for b in range(m)
            ],
            dtype=np.uint8,
        )
        if rank_np(orbit) != dimension:
            continue
        unique = np.unique(orbit, axis=0)
        key = (len(unique), int(seed.sum()), seed.tobytes())
        candidates.append((key, seed.copy(), orbit))
    if not candidates:
        raise RuntimeError("pole ideal has no machine-verified cyclic generator")
    _, seed, orbit = min(candidates, key=lambda item: item[0])
    return seed, orbit


def _stabilizer_subgroup(seed: np.ndarray, ell: int, m: int) -> list[tuple[int, int]]:
    return [
        (a, b)
        for a in range(ell)
        for b in range(m)
        if np.array_equal(_translate_block(seed, ell, m, a, b), seed)
    ]


def _subgroup_orbit_representatives(
    ell: int, m: int, subgroup: list[tuple[int, int]]
) -> tuple[list[int], list[list[int]]]:
    unseen = set(range(ell * m))
    representatives: list[int] = []
    orbits: list[list[int]] = []
    while unseen:
        index = min(unseen)
        i, j = divmod(index, m)
        orbit = sorted(
            {((i + a) % ell) * m + ((j + b) % m) for a, b in subgroup}
        )
        representatives.append(orbit[0])
        orbits.append(orbit)
        unseen.difference_update(orbit)
    return representatives, orbits


def _translation_preserves_rowspace(
    rows: np.ndarray, ell: int, m: int, a: int, b: int
) -> bool:
    translated = np.asarray(
        [_translate_two_blocks(row, ell, m, a, b) for row in rows], dtype=np.uint8
    )
    rank = rank_np(rows)
    return rank_np(np.vstack([rows, translated])) == rank


def build_orbit_cover(problem: dict[str, Any]) -> dict[str, Any]:
    """Construct and machine-check the two-sector principal-pole cover."""
    ell, m, block, k = problem["ell"], problem["m"], problem["block"], problem["k"]
    A, B, HX, HZ = problem["A"], problem["B"], problem["HX"], problem["HZ"]
    left_a = _independent_rows(nullspace_np(A.T))
    left_b = _independent_rows(nullspace_np(B.T))
    I = _intersection(left_a, left_b)
    J = I[:, bar_permutation(ell, m)]
    pole_dim = int(I.shape[0])
    if 2 * pole_dim != k:
        raise RuntimeError(f"principal pole dimension 2*{pole_dim} != k={k}")

    x_pole = np.vstack([_embed(row, side) for side in (0, 1) for row in I])
    z_pole = np.vstack([_embed(row, side) for side in (0, 1) for row in J])
    x_kernel = not np.any(HZ @ x_pole.T % 2)
    z_kernel = not np.any(HX @ z_pole.T % 2)
    x_quotient_np, x_quotient_bits = _quotient_rank(x_pole, HX)
    z_quotient_np, z_quotient_bits = _quotient_rank(z_pole, HZ)

    seed, pole_orbit = _cyclic_generator(I, ell, m)
    subgroup = _stabilizer_subgroup(seed, ell, m)
    representatives, subgroup_orbits = _subgroup_orbit_representatives(ell, m, subgroup)
    subgroup_free = all(len(orbit) == len(subgroup) for orbit in subgroup_orbits)
    translations_preserve_code = all(
        _translation_preserves_rowspace(HX, ell, m, a, b)
        and _translation_preserves_rowspace(HZ, ell, m, a, b)
        for a in range(ell)
        for b in range(m)
    )

    sectors: list[dict[str, Any]] = []
    union: list[np.ndarray] = []
    for side in (0, 1):
        detector = _embed(seed, side)
        orbit = np.asarray(
            [_embed(row, side) for row in pole_orbit], dtype=np.uint8
        )
        orbit_rank_np, orbit_rank_bits = _quotient_rank(orbit, HX)
        anchor = representatives + [block + index for index in representatives]
        detector_fixed = all(
            np.array_equal(
                _translate_two_blocks(detector, ell, m, a, b), detector
            )
            for a, b in subgroup
        )
        sectors.append(
            {
                "sector": side,
                "detector": detector,
                "detector_sha256": matrix_sha256(detector[None, :]),
                "detector_weight": int(detector.sum()),
                "orbit": orbit,
                "orbit_unique_count": int(len(np.unique(orbit, axis=0))),
                "orbit_quotient_rank": orbit_rank_np,
                "orbit_quotient_rank_bitset": orbit_rank_bits,
                "stabilizer_size": len(subgroup),
                "stabilizer_shifts": [[a, b] for a, b in subgroup],
                "detector_fixed_by_stabilizer": detector_fixed,
                "anchor_coordinates": anchor,
                "anchor_orbit_count_per_block": len(representatives),
            }
        )
        union.extend(orbit)

    union_rows = np.asarray(union, dtype=np.uint8)
    union_np, union_bits = _quotient_rank(union_rows, HX)
    pairing = (union_rows.astype(np.int64) @ z_pole.T.astype(np.int64) % 2).astype(
        np.uint8
    )
    pairing_rank_np = rank_np(pairing)
    pairing_rank_bits = rank_bitset(rows_to_bitsets(pairing), pairing.shape[1])
    complete = bool(
        x_kernel
        and z_kernel
        and x_quotient_np == x_quotient_bits == k
        and z_quotient_np == z_quotient_bits == k
        and union_np == union_bits == k
        and pairing_rank_np == pairing_rank_bits == k
        and translations_preserve_code
        and subgroup_free
        and all(
            sector["orbit_quotient_rank"]
            == sector["orbit_quotient_rank_bitset"]
            == pole_dim
            and sector["detector_fixed_by_stabilizer"]
            for sector in sectors
        )
    )
    if not complete:
        raise RuntimeError("reciprocal-pole orbit cover failed its completeness gates")
    return {
        "schema": "exp056-css-pole-orbit-cover-v1",
        "complete": complete,
        "pole_dimension": pole_dim,
        "required_sectors": len(sectors),
        "I": I,
        "J": J,
        "x_pole": x_pole,
        "z_pole": z_pole,
        "x_pole_in_kernel": x_kernel,
        "z_pole_in_kernel": z_kernel,
        "x_pole_quotient_rank_numpy": x_quotient_np,
        "x_pole_quotient_rank_bitset": x_quotient_bits,
        "z_pole_quotient_rank_numpy": z_quotient_np,
        "z_pole_quotient_rank_bitset": z_quotient_bits,
        "translations_preserve_both_check_rowspaces": translations_preserve_code,
        "subgroup_action_free": subgroup_free,
        "union_quotient_rank_numpy": union_np,
        "union_quotient_rank_bitset": union_bits,
        "pairing_rank_numpy": pairing_rank_np,
        "pairing_rank_bitset": pairing_rank_bits,
        "sectors": sectors,
        "coverage_argument": (
            "The selected detector translation orbits span the full X-logical "
            "quotient. A nontrivial Z class pairs with some orbit member; "
            "inverse translation produces an equal-weight candidate pairing "
            "with its stored representative."
        ),
        "anchor_argument": (
            "A fixed detector sector may use only translations in the "
            "detector stabilizer. One representative from every subgroup "
            "orbit in each physical block gives a sound support clause."
        ),
    }


def _cover_metadata(cover: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in cover.items()
        if key not in {"I", "J", "x_pole", "z_pole", "sectors"}
    } | {
        "I_sha256": matrix_sha256(cover["I"]),
        "J_sha256": matrix_sha256(cover["J"]),
        "sectors": [
            {key: value for key, value in sector.items() if key not in {"detector", "orbit"}}
            | {"orbit_sha256": matrix_sha256(sector["orbit"])}
            for sector in cover["sectors"]
        ],
    }


def _permutation_map(
    problem: dict[str, Any], a: int, b: int, reflect_x: bool
) -> np.ndarray:
    """Old physical coordinate -> new coordinate for one target automorphism."""
    ell, m, block = problem["ell"], problem["m"], problem["block"]
    mapping = np.empty(problem["n"], dtype=int)
    for side in (0, 1):
        offset = side * block
        for i in range(ell):
            for j in range(m):
                new_i = ((-i if reflect_x else i) + a) % ell
                new_j = (j + b) % m
                mapping[offset + i * m + j] = offset + new_i * m + new_j
    return mapping


def _permute_vector(vector: np.ndarray, mapping: np.ndarray) -> np.ndarray:
    output = np.zeros_like(vector)
    output[mapping] = vector
    return output


def _permutation_preserves_rowspace(rows: np.ndarray, mapping: np.ndarray) -> bool:
    translated = np.asarray([_permute_vector(row, mapping) for row in rows], dtype=np.uint8)
    return rank_np(np.vstack([rows, translated])) == rank_np(rows)


def _coset_reducer(modulo: np.ndarray):
    reduced, pivots = rref_np(modulo)
    reduced = reduced[: rank_np(reduced)]
    pivots = pivots[: len(reduced)]

    def reduce(vector: np.ndarray) -> np.ndarray:
        output = np.asarray(vector, dtype=np.uint8).copy()
        for row, pivot in zip(reduced, pivots):
            if output[pivot]:
                output ^= row
        return output

    return reduce


def _mask_vector(mask: int, basis: np.ndarray) -> np.ndarray:
    rows = [basis[index] for index in range(len(basis)) if (mask >> index) & 1]
    if not rows:
        return np.zeros(basis.shape[1], dtype=np.uint8)
    return np.bitwise_xor.reduce(np.asarray(rows, dtype=np.uint8), axis=0)


def _act_mask(mask: int, basis_images: list[int]) -> int:
    output = 0
    for index, image in enumerate(basis_images):
        if (mask >> index) & 1:
            output ^= image
    return output


def build_class_orbit_cover(
    problem: dict[str, Any], pole_cover: dict[str, Any]
) -> dict[str, Any]:
    """Partition every nonzero Z-logical class under verified automorphisms.

    This is a deliberately stronger fallback than the two functional sectors:
    each SAT subproblem fixes all ``k`` logical syndrome bits.  The target has
    255 nonzero classes but only 20 orbits under translations and x-reflection.
    """
    if not pole_cover.get("complete"):
        raise RuntimeError("class cover requires a complete pole cover")
    k, n = problem["k"], problem["n"]
    z_pole, x_pole = pole_cover["z_pole"], pole_cover["x_pole"]
    reduce_coset = _coset_reducer(problem["HZ"])

    residual_to_mask: dict[bytes, int] = {}
    for mask in range(1 << k):
        residual = reduce_coset(_mask_vector(mask, z_pole))
        key = residual.tobytes()
        if key in residual_to_mask:
            raise RuntimeError("pole coordinates are not injective modulo Z stabilizers")
        residual_to_mask[key] = mask

    reflection = _permutation_map(problem, 0, 0, True)
    reflection_valid = bool(
        _permutation_preserves_rowspace(problem["HX"], reflection)
        and _permutation_preserves_rowspace(problem["HZ"], reflection)
    )
    group: list[dict[str, Any]] = []
    for reflect_x in ([False, True] if reflection_valid else [False]):
        for a in range(problem["ell"]):
            for b in range(problem["m"]):
                mapping = _permutation_map(problem, a, b, reflect_x)
                preserves = bool(
                    _permutation_preserves_rowspace(problem["HX"], mapping)
                    and _permutation_preserves_rowspace(problem["HZ"], mapping)
                )
                if not preserves:
                    raise RuntimeError("candidate class-cover automorphism does not preserve code")
                basis_images: list[int] = []
                for row in z_pole:
                    residual = reduce_coset(_permute_vector(row, mapping))
                    image = residual_to_mask.get(residual.tobytes())
                    if image is None:
                        raise RuntimeError("automorphism left the logical quotient")
                    basis_images.append(image)
                if rank_bitset(basis_images, k) != k:
                    raise RuntimeError("automorphism action on quotient is not invertible")
                group.append(
                    {
                        "a": a,
                        "b": b,
                        "reflect_x": reflect_x,
                        "mapping": mapping,
                        "basis_images": basis_images,
                    }
                )

    unseen = set(range(1, 1 << k))
    representatives: list[dict[str, Any]] = []
    covered: set[int] = set()
    while unseen:
        seed_mask = min(unseen)
        orbit = sorted({_act_mask(seed_mask, item["basis_images"]) for item in group})
        if not orbit or orbit[0] == 0:
            raise RuntimeError("nonzero class orbit reached zero")
        representative_mask = min(orbit)
        vector = _mask_vector(representative_mask, z_pole)
        syndrome = (
            x_pole.astype(np.int64) @ vector[:, None].astype(np.int64) % 2
        )[:, 0].astype(np.uint8)
        stabilizer = [
            item
            for item in group
            if _act_mask(representative_mask, item["basis_images"])
            == representative_mask
        ]
        coordinate_unseen = set(range(n))
        anchor: list[int] = []
        coordinate_orbits: list[list[int]] = []
        while coordinate_unseen:
            coordinate = min(coordinate_unseen)
            coordinate_orbit = sorted(
                {int(item["mapping"][coordinate]) for item in stabilizer}
            )
            anchor.append(coordinate_orbit[0])
            coordinate_orbits.append(coordinate_orbit)
            coordinate_unseen.difference_update(coordinate_orbit)
        representatives.append(
            {
                "class_index": len(representatives),
                "mask": representative_mask,
                "mask_hex": f"0x{representative_mask:0{(k + 3) // 4}x}",
                "orbit_masks": orbit,
                "orbit_size": len(orbit),
                "vector": vector,
                "vector_sha256": matrix_sha256(vector[None, :]),
                "syndrome": syndrome,
                "syndrome_bits": syndrome.astype(int).tolist(),
                "syndrome_nonzero": bool(syndrome.any()),
                "class_stabilizer_size": len(stabilizer),
                "orbit_stabilizer_exact": (
                    len(orbit) * len(stabilizer) == len(group)
                ),
                "stabilizer_mappings": [
                    item["mapping"].copy() for item in stabilizer
                ],
                "anchor_coordinates": anchor,
                "coordinate_orbit_count": len(coordinate_orbits),
                "coordinate_orbits_cover_all": (
                    sorted(index for orbit_row in coordinate_orbits for index in orbit_row)
                    == list(range(n))
                ),
            }
        )
        covered.update(orbit)
        unseen.difference_update(orbit)

    unique_syndromes = {
        tuple((x_pole.astype(np.int64) @ _mask_vector(mask, z_pole) % 2).tolist())
        for mask in range(1, 1 << k)
    }
    complete = bool(
        len(residual_to_mask) == 1 << k
        and len(group) in {problem["ell"] * problem["m"], 2 * problem["ell"] * problem["m"]}
        and covered == set(range(1, 1 << k))
        and len(unique_syndromes) == (1 << k) - 1
        and all(
            item["syndrome_nonzero"]
            and item["coordinate_orbits_cover_all"]
            and item["orbit_stabilizer_exact"]
            for item in representatives
        )
    )
    if not complete:
        raise RuntimeError("logical-class orbit cover failed completeness gates")
    return {
        "schema": "exp056-css-class-orbit-cover-v1",
        "complete": complete,
        "nonzero_classes": (1 << k) - 1,
        "automorphism_group_size": len(group),
        "x_reflection_admitted": reflection_valid,
        "class_orbits": len(representatives),
        "covered_class_masks": len(covered),
        "group": group,
        "representatives": representatives,
        "coverage_argument": (
            "Every nonzero pole-coordinate mask is present exactly once in the "
            "orbit partition. Each fixed-class query pins all k pairings; "
            "automorphisms transport its result across the stored orbit."
        ),
    }


def _class_cover_metadata(class_cover: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in class_cover.items()
        if key not in {"group", "representatives"}
    } | {
        "group_sha256": canonical_json_sha256(
            [
                {
                    "a": item["a"],
                    "b": item["b"],
                    "reflect_x": item["reflect_x"],
                    "mapping_sha256": canonical_json_sha256(item["mapping"].astype(int).tolist()),
                    "basis_images": item["basis_images"],
                }
                for item in class_cover["group"]
            ]
        ),
        "representatives": [
            {
                key: value
                for key, value in item.items()
                if key
                not in {"vector", "syndrome", "orbit_masks", "stabilizer_mappings"}
            }
            | {
                "orbit_masks_sha256": canonical_json_sha256(item["orbit_masks"]),
                "stabilizer_mappings_sha256": canonical_json_sha256(
                    [
                        mapping.astype(int).tolist()
                        for mapping in item["stabilizer_mappings"]
                    ]
                ),
            }
            for item in class_cover["representatives"]
        ],
    }


def _component_transform(problem: dict[str, Any], vector: np.ndarray) -> np.ndarray:
    """Binary x=1/even decomposition for a three-row lattice."""
    if problem["ell"] != 3:
        raise ValueError("component formulation currently requires ell=3")
    m = problem["m"]
    symbols = 2 * m
    constant = np.zeros(symbols, dtype=np.uint8)
    even_0 = np.zeros(symbols, dtype=np.uint8)
    even_1 = np.zeros(symbols, dtype=np.uint8)
    for side in (0, 1):
        rows = np.asarray(
            vector[side * 3 * m : (side + 1) * 3 * m], dtype=np.uint8
        ).reshape(3, m)
        target = slice(side * m, (side + 1) * m)
        constant[target] = rows[0] ^ rows[1] ^ rows[2]
        even_0[target] = rows[1] ^ rows[2]
        even_1[target] = rows[0] ^ rows[2]
    return np.concatenate([constant, even_0, even_1])


def _component_inverse(problem: dict[str, Any], transformed: np.ndarray) -> np.ndarray:
    if problem["ell"] != 3:
        raise ValueError("component formulation currently requires ell=3")
    m = problem["m"]
    symbols = 2 * m
    transformed = np.asarray(transformed, dtype=np.uint8)
    constant = transformed[:symbols]
    even_0 = transformed[symbols : 2 * symbols]
    even_1 = transformed[2 * symbols :]
    output = np.zeros(problem["n"], dtype=np.uint8)
    for side in (0, 1):
        source = slice(side * m, (side + 1) * m)
        rows = np.vstack(
            [
                constant[source] ^ even_0[source],
                constant[source] ^ even_1[source],
                constant[source] ^ even_0[source] ^ even_1[source],
            ]
        )
        output[side * 3 * m : (side + 1) * 3 * m] = rows.reshape(-1)
    return output


def build_component_decomposition(
    problem: dict[str, Any], pole_cover: dict[str, Any]
) -> dict[str, Any]:
    """Exact direct-sum decomposition of S_Z for ell=3.

    The x+x^-1 operator is zero on the constant x-sector and identity on the
    two-dimensional even sector.  The transform is local and invertible; the
    stabilizer rowspace splits into ranks 27 and 50 for the target.
    """
    if problem["ell"] != 3:
        raise ValueError("component decomposition currently requires ell=3")
    transformed_rows = np.asarray(
        [_component_transform(problem, row) for row in problem["HZ"]],
        dtype=np.uint8,
    )
    symbols = 2 * problem["m"]
    constant = _independent_rows(transformed_rows[:, :symbols])
    even = _independent_rows(transformed_rows[:, symbols:])
    zero_ce = np.zeros((len(constant), 2 * symbols), dtype=np.uint8)
    zero_ec = np.zeros((len(even), symbols), dtype=np.uint8)
    direct_sum = np.vstack(
        [np.hstack([constant, zero_ce]), np.hstack([zero_ec, even])]
    )
    transformed_rank = rank_np(transformed_rows)
    direct_rank = rank_np(direct_sum)
    rowspaces_equal = rank_np(np.vstack([transformed_rows, direct_sum])) == transformed_rank
    logical = np.asarray(
        [_component_transform(problem, row) for row in pole_cover["z_pole"]],
        dtype=np.uint8,
    )
    logical_constant_zero = not logical[:, :symbols].any()
    even_logical = logical[:, symbols:]
    even_quotient_rank = rank_np(np.vstack([even, even_logical])) - rank_np(even)
    constant_checks = _independent_rows(nullspace_np(constant))
    even_checks = _independent_rows(nullspace_np(even))
    round_trip = all(
        np.array_equal(
            _component_inverse(problem, _component_transform(problem, basis)),
            basis,
        )
        for basis in np.eye(problem["n"], dtype=np.uint8)
    )
    complete = bool(
        round_trip
        and transformed_rank == rank_np(problem["HZ"])
        and direct_rank == transformed_rank
        and rowspaces_equal
        and logical_constant_zero
        and even_quotient_rank == problem["k"]
        and constant_checks.shape[0] == symbols - len(constant)
        and even_checks.shape[0] == 2 * symbols - len(even)
    )
    if not complete:
        raise RuntimeError("three-row component decomposition failed")
    return {
        "schema": "exp056-three-row-component-decomposition-v1",
        "complete": complete,
        "symbols": symbols,
        "transformed_stabilizer": transformed_rows,
        "constant_stabilizer": constant,
        "even_stabilizer": even,
        "constant_checks": constant_checks,
        "even_checks": even_checks,
        "even_logical": even_logical,
        "constant_stabilizer_rank": rank_np(constant),
        "even_stabilizer_rank": rank_np(even),
        "transformed_stabilizer_rank": transformed_rank,
        "direct_sum_rank": direct_rank,
        "direct_sum_rowspace_equal": bool(rowspaces_equal),
        "transform_round_trip_basis": round_trip,
        "logical_constant_projection_zero": logical_constant_zero,
        "even_logical_quotient_rank": int(even_quotient_rank),
    }


def _component_metadata(components: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in components.items()
        if key
        not in {
            "transformed_stabilizer",
            "constant_stabilizer",
            "even_stabilizer",
            "constant_checks",
            "even_checks",
            "even_logical",
        }
    } | {
        "constant_checks_sha256": matrix_sha256(components["constant_checks"]),
        "even_checks_sha256": matrix_sha256(components["even_checks"]),
        "even_logical_sha256": matrix_sha256(components["even_logical"]),
    }


def bb_duality_proof(problem: dict[str, Any]) -> dict[str, Any]:
    """Exact coordinate isometry exchanging X and Z distances of every BB code."""
    block = problem["block"]
    q = np.eye(block, dtype=np.uint8)[bar_permutation(problem["ell"], problem["m"])]
    zero = np.zeros_like(q)
    permutation_matrix = np.block([[zero, q], [q, zero]])
    permutation = np.argmax(permutation_matrix, axis=1)
    identity = np.eye(problem["n"], dtype=np.uint8)
    maps_xz = np.array_equal((q @ problem["HX"] @ permutation_matrix) % 2, problem["HZ"])
    maps_zx = np.array_equal((q @ problem["HZ"] @ permutation_matrix) % 2, problem["HX"])
    bijection = np.array_equal(np.sort(permutation), np.arange(problem["n"]))
    involution = np.array_equal((permutation_matrix @ permutation_matrix) % 2, identity)
    valid = bool(maps_xz and maps_zx and bijection and involution)
    if not valid:
        raise RuntimeError("BB X/Z duality isometry failed")
    return {
        "schema": "exp056-bb-xz-duality-v1",
        "valid": valid,
        "permutation": permutation,
        "permutation_sha256": canonical_json_sha256(permutation.astype(int).tolist()),
        "permutation_is_bijection": bool(bijection),
        "permutation_is_involution": bool(involution),
        "maps_HX_to_HZ_exactly": bool(maps_xz),
        "maps_HZ_to_HX_exactly": bool(maps_zx),
        "conclusion": "d_X = d_Z; one certified side determines CSS distance",
    }


def default_witness(problem: dict[str, Any]) -> np.ndarray:
    vector = np.zeros(problem["n"], dtype=np.uint8)
    support = problem["record"].get("witness_support")
    if support is None:
        raise ValueError("problem record has no explicit upper-bound witness")
    vector[np.asarray(support, dtype=int)] = 1
    return vector


def _vector_bitset(vector: np.ndarray) -> int:
    return sum(int(bit) << index for index, bit in enumerate(vector))


def _css_word_valid_bitset(
    vector: np.ndarray, constraints: np.ndarray, stabilizers: np.ndarray
) -> bool:
    bits = _vector_bitset(vector)
    parity_clean = all(
        (row & bits).bit_count() % 2 == 0 for row in rows_to_bitsets(constraints)
    )
    outside = solve_bitset(rows_to_bitsets(stabilizers), len(vector), bits) is None
    return bool(parity_clean and outside)


def verify_upper_witness(problem: dict[str, Any], witness: np.ndarray) -> dict[str, Any]:
    witness = np.asarray(witness, dtype=np.uint8) & 1
    if witness.shape != (problem["n"],):
        raise ValueError("witness width mismatch")
    duality = bb_duality_proof(problem)
    x_witness = (witness @ np.eye(problem["n"], dtype=np.uint8)[duality["permutation"]]) % 2
    z_valid_np = bool(
        not np.any(problem["HX"] @ witness % 2)
        and rank_np(np.vstack([problem["HZ"], witness])) == rank_np(problem["HZ"]) + 1
    )
    x_valid_np = bool(
        not np.any(problem["HZ"] @ x_witness % 2)
        and rank_np(np.vstack([problem["HX"], x_witness])) == rank_np(problem["HX"]) + 1
    )
    z_valid_bits = _css_word_valid_bitset(witness, problem["HX"], problem["HZ"])
    x_valid_bits = _css_word_valid_bitset(x_witness, problem["HZ"], problem["HX"])
    result = {
        "weight": int(witness.sum()),
        "support": [int(index) for index in np.flatnonzero(witness)],
        "vector_sha256": matrix_sha256(witness[None, :]),
        "z_valid_numpy": z_valid_np,
        "z_valid_bitset": z_valid_bits,
        "x_support": [int(index) for index in np.flatnonzero(x_witness)],
        "x_vector_sha256": matrix_sha256(x_witness[None, :]),
        "x_valid_numpy": x_valid_np,
        "x_valid_bitset": x_valid_bits,
        "duality_preserves_weight": int(x_witness.sum()) == int(witness.sum()),
    }
    if not all(
        result[key]
        for key in (
            "z_valid_numpy",
            "z_valid_bitset",
            "x_valid_numpy",
            "x_valid_bitset",
            "duality_preserves_weight",
        )
    ):
        raise RuntimeError("upper-bound witness failed independent verification")
    return result


def _sector_instance(problem: dict[str, Any], sector: dict[str, Any]) -> DecisionInstance:
    return DecisionInstance(
        parity_rows=problem["HX"],
        pairing_rows=sector["detector"][None, :],
        groups=[[index] for index in range(problem["n"])],
        kind=f"exp056_css_z_orbit_sector_{sector['sector']}",
        meta={
            "ell": problem["ell"],
            "m": problem["m"],
            "sector": sector["sector"],
            "anchor_group": "detector_stabilizer",
        },
        symmetry_clause=list(sector["anchor_coordinates"]),
    )


def _verify_sector_witness(
    problem: dict[str, Any], instance: DecisionInstance, vector: np.ndarray
) -> dict[str, Any]:
    generic = verify_witness_two_paths(instance, vector)
    outside_np = rank_np(np.vstack([problem["HZ"], vector])) == rank_np(problem["HZ"]) + 1
    outside_bits = solve_bitset(
        rows_to_bitsets(problem["HZ"]), problem["n"], _vector_bitset(vector)
    ) is None
    valid = bool(generic["valid"] and outside_np and outside_bits)
    return {
        **generic,
        "outside_z_stabilizer_numpy": bool(outside_np),
        "outside_z_stabilizer_bitset": bool(outside_bits),
        "valid": valid,
    }


def decide_pole_monolithic(
    problem: dict[str, Any],
    cover: dict[str, Any],
    weight_cap: int,
    *,
    conflict_budget: int = 0,
    solver_name: str = SOLVER_NAME,
    xor_encoding: str = "sequential",
    card_encoding: str = "seqcounter",
) -> dict[str, Any]:
    """One anchored CNF using the full pole-pairing disjunction."""
    import pysat
    from pysat.solvers import Solver

    if not cover.get("complete"):
        raise RuntimeError("monolithic pole decision requires a complete cover")
    cnf, weight_literals = build_pole_monolithic_cnf(
        problem,
        cover,
        weight_cap,
        xor_encoding=xor_encoding,
        card_encoding=card_encoding,
    )
    digest = (
        canonical_json_sha256(
            {
                "cnf_sha256": cnf_sha256(cnf),
                "native_atmost": [weight_literals, int(weight_cap)],
            }
        )
        if card_encoding == "native"
        else cnf_sha256(cnf)
    )
    if card_encoding == "native" and solver_name != "minicard":
        raise ValueError("native cardinality requires solver_name='minicard'")
    started = time.perf_counter()
    with Solver(name=solver_name, bootstrap_with=cnf.clauses) as engine:
        if card_encoding == "native":
            engine.add_atmost(weight_literals, int(weight_cap))
        if conflict_budget > 0:
            engine.conf_budget(int(conflict_budget))
            answer = engine.solve_limited(expect_interrupt=False)
        else:
            answer = engine.solve()
        try:
            stats = dict(engine.accum_stats())
        except NotImplementedError:
            stats = {}
        model = engine.get_model() if answer else None
    wall = time.perf_counter() - started
    status = {True: "SAT", False: "UNSAT", None: "UNDECIDED_BUDGET"}[answer]
    record: dict[str, Any] = {
        "status": status,
        "weight_cap": int(weight_cap),
        "cnf_sha256": digest,
        "instance_sha256": digest,
        "xor_encoding": xor_encoding,
        "card_encoding": card_encoding,
        "solver": {
            "name": f"PySAT {solver_name}",
            "version": pysat.__version__,
            "backend": "pysat",
            "conflict_budget": int(conflict_budget),
            "clauses": len(cnf.clauses),
            "cnf_variables": cnf.nv,
            "wall_time_s": wall,
            "stats": {key: int(value) for key, value in stats.items()},
        },
    }
    if model is not None:
        assignment = {abs(literal): literal > 0 for literal in model}
        vector = np.asarray(
            [int(assignment.get(index + 1, False)) for index in range(problem["n"])],
            dtype=np.uint8,
        )
        instance = DecisionInstance(
            parity_rows=problem["HX"],
            pairing_rows=cover["x_pole"],
            groups=[[index] for index in range(problem["n"])],
            kind="exp056_css_z_pole_monolithic",
        )
        verification = _verify_sector_witness(problem, instance, vector)
        if not verification["valid"] or int(vector.sum()) > int(weight_cap):
            raise RuntimeError("monolithic pole witness failed CSS verification")
        record["vector"] = vector.astype(int).tolist()
        record["weight"] = int(vector.sum())
        record["verification"] = verification
    return record


def _xor_output(cnf, pool, literals: list[int], encoding: str) -> int | None:
    def add_gate(left: int, right: int) -> int:
        auxiliary = pool.id()
        cnf.extend(
            [
                [-left, -right, -auxiliary],
                [left, right, -auxiliary],
                [left, -right, auxiliary],
                [-left, right, auxiliary],
            ]
        )
        return auxiliary

    if not literals:
        return None
    if encoding not in {"sequential", "balanced"}:
        raise ValueError("monolithic XOR outputs support sequential or balanced encoding")
    work = list(literals)
    if encoding == "balanced":
        while len(work) > 1:
            next_level = [
                add_gate(work[index], work[index + 1])
                for index in range(0, len(work) - 1, 2)
            ]
            if len(work) % 2:
                next_level.append(work[-1])
            work = next_level
        return work[0]
    output = work[0]
    for literal in work[1:]:
        output = add_gate(output, literal)
    return output


def build_pole_monolithic_cnf(
    problem: dict[str, Any],
    cover: dict[str, Any],
    weight_cap: int,
    *,
    xor_encoding: str,
    card_encoding: str,
):
    from pysat.formula import CNF, IDPool

    n = problem["n"]
    cnf = CNF()
    pool = IDPool(start_from=n + 1)
    cnf.append([1, problem["block"] + 1])
    for row in problem["HX"]:
        _add_xor_constraint(
            cnf,
            pool,
            [int(index) + 1 for index in np.flatnonzero(row)],
            0,
            xor_encoding,
        )
    pairing_outputs = [
        _xor_output(
            cnf,
            pool,
            [int(index) + 1 for index in np.flatnonzero(row)],
            xor_encoding,
        )
        for row in cover["x_pole"]
    ]
    cnf.append([output for output in pairing_outputs if output is not None])
    weight_literals = list(range(1, n + 1))
    _add_cardinality_constraint(
        cnf, pool, weight_literals, weight_cap, card_encoding
    )
    return cnf, weight_literals


def _add_xor_constraint(
    cnf, pool, literals: list[int], rhs: int, encoding: str
) -> None:
    def add_gate(left: int, right: int) -> int:
        auxiliary = pool.id()
        cnf.extend(
            [
                [-left, -right, -auxiliary],
                [left, right, -auxiliary],
                [left, -right, auxiliary],
                [-left, right, auxiliary],
            ]
        )
        return auxiliary

    if encoding == "direct-small" and len(literals) <= 8:
        for assignment in range(1 << len(literals)):
            if assignment.bit_count() % 2 == rhs:
                continue
            cnf.append(
                [
                    -literal if (assignment >> index) & 1 else literal
                    for index, literal in enumerate(literals)
                ]
            )
        return
    if encoding not in {"sequential", "balanced", "direct-small"}:
        raise ValueError(f"unknown XOR encoding {encoding!r}")
    work = list(literals)
    if encoding == "balanced":
        while len(work) > 1:
            next_level = [
                add_gate(work[index], work[index + 1])
                for index in range(0, len(work) - 1, 2)
            ]
            if len(work) % 2:
                next_level.append(work[-1])
            work = next_level
        accumulator = work[0] if work else None
    else:
        accumulator: int | None = None
        for literal in work:
            accumulator = literal if accumulator is None else add_gate(accumulator, literal)
    if accumulator is None:
        if rhs:
            cnf.append([])
    else:
        cnf.append([accumulator] if rhs else [-accumulator])


def _add_cardinality_constraint(
    cnf, pool, literals: list[int], bound: int, encoding: str
) -> None:
    if encoding == "native":
        return
    from pysat.card import CardEnc, EncType

    choices = {
        "seqcounter": EncType.seqcounter,
        "totalizer": EncType.totalizer,
        "cardnetwrk": EncType.cardnetwrk,
        "mtotalizer": EncType.mtotalizer,
        "kmtotalizer": EncType.kmtotalizer,
    }
    if encoding not in choices:
        raise ValueError(f"unknown cardinality encoding {encoding!r}")
    cnf.extend(
        CardEnc.atmost(
            lits=literals,
            bound=int(bound),
            vpool=pool,
            encoding=choices[encoding],
        ).clauses
    )


def _add_class_symmetry_break(
    cnf,
    pool,
    physical_literals: list[int],
    representative: dict[str, Any],
    mode: str,
) -> None:
    if mode not in {"anchor", "lex", "anchor+lex"}:
        raise ValueError(f"unknown class symmetry break {mode!r}")
    if mode in {"anchor", "anchor+lex"}:
        cnf.append(
            [
                physical_literals[int(index)]
                for index in representative["anchor_coordinates"]
            ]
        )
    if mode not in {"lex", "anchor+lex"}:
        return
    identity = np.arange(len(physical_literals))
    for mapping in representative["stabilizer_mappings"]:
        if np.array_equal(mapping, identity):
            continue
        inverse = np.argsort(mapping)
        equal_prefix = pool.id()
        cnf.append([equal_prefix])
        for coordinate in range(len(physical_literals)):
            left = physical_literals[coordinate]
            right = physical_literals[int(inverse[coordinate])]
            if left == right:
                continue
            cnf.append([-equal_prefix, -left, right])
            next_equal = pool.id()
            cnf.extend(
                [
                    [-next_equal, equal_prefix],
                    [-next_equal, -left, right],
                    [-next_equal, left, -right],
                    [-equal_prefix, left, right, next_equal],
                    [-equal_prefix, -left, -right, next_equal],
                ]
            )
            equal_prefix = next_equal


def build_class_cnf(
    problem: dict[str, Any],
    pole_cover: dict[str, Any],
    representative: dict[str, Any],
    weight_cap: int,
    *,
    xor_encoding: str = "sequential",
    card_encoding: str = "seqcounter",
    symmetry_break: str = "anchor",
):
    """CNF for one fully fixed nonzero logical class."""
    from pysat.formula import CNF, IDPool

    n = problem["n"]
    cnf = CNF()
    pool = IDPool(start_from=n + 1)

    def add_xor(row: np.ndarray, rhs: int) -> None:
        _add_xor_constraint(
            cnf,
            pool,
            [int(index) + 1 for index in np.flatnonzero(row)],
            rhs,
            xor_encoding,
        )

    _add_class_symmetry_break(
        cnf, pool, list(range(1, n + 1)), representative, symmetry_break
    )
    for row in problem["HX"]:
        add_xor(row, 0)
    for row, rhs in zip(pole_cover["x_pole"], representative["syndrome"]):
        add_xor(row, int(rhs))
    _add_cardinality_constraint(
        cnf, pool, list(range(1, n + 1)), weight_cap, card_encoding
    )
    return cnf


def _sparse_independent_rows(rows: np.ndarray) -> np.ndarray:
    chosen: list[np.ndarray] = []
    current_rank = 0
    for row in np.asarray(rows, dtype=np.uint8):
        trial = np.asarray([*chosen, row], dtype=np.uint8)
        rank = rank_np(trial)
        if rank > current_rank:
            chosen.append(row.copy())
            current_rank = rank
    return np.asarray(chosen, dtype=np.uint8)


def build_generator_class_cnf(
    problem: dict[str, Any],
    representative: dict[str, Any],
    weight_cap: int,
    *,
    xor_encoding: str = "sequential",
    card_encoding: str = "seqcounter",
    symmetry_break: str = "anchor",
):
    """Affine coset encoder ``v = representative + lambda H_Z``.

    Independent original H_Z rows retain weight six, so every physical output
    is defined by a short XOR instead of imposing dense logical syndromes.
    """
    from pysat.formula import CNF, IDPool

    generator = _sparse_independent_rows(problem["HZ"])
    if len(generator) != rank_np(problem["HZ"]):
        raise RuntimeError("sparse stabilizer generator selection lost rank")
    offset = len(generator)
    n = problem["n"]
    cnf = CNF()
    pool = IDPool(start_from=offset + n + 1)
    for coordinate in range(n):
        inputs = [int(index) + 1 for index in np.flatnonzero(generator[:, coordinate])]
        output = offset + coordinate + 1
        _add_xor_constraint(
            cnf,
            pool,
            [*inputs, output],
            int(representative["vector"][coordinate]),
            xor_encoding,
        )
    physical_literals = list(range(offset + 1, offset + n + 1))
    _add_class_symmetry_break(
        cnf, pool, physical_literals, representative, symmetry_break
    )
    _add_cardinality_constraint(
        cnf, pool, physical_literals, weight_cap, card_encoding
    )
    return cnf, generator


def build_component_class_cnf(
    problem: dict[str, Any],
    components: dict[str, Any],
    representative: dict[str, Any],
    weight_cap: int,
    *,
    xor_encoding: str = "sequential",
    card_encoding: str = "seqcounter",
    symmetry_break: str = "anchor",
):
    """Pinned-class CNF in the exact constant/even x-sector coordinates."""
    from pysat.formula import CNF, IDPool

    n = problem["n"]
    symbols = components["symbols"]
    cnf = CNF()
    pool = IDPool(start_from=2 * n + 1)

    def add(literals: list[int], rhs: int) -> None:
        _add_xor_constraint(cnf, pool, literals, rhs, xor_encoding)

    transformed_rep = _component_transform(problem, representative["vector"])
    if transformed_rep[:symbols].any():
        raise RuntimeError("logical representative leaked into constant x-sector")
    even_rep = transformed_rep[symbols:]
    even_rhs = (
        components["even_checks"].astype(np.int64)
        @ even_rep[:, None].astype(np.int64)
        % 2
    )[:, 0].astype(np.uint8)

    for row in components["constant_checks"]:
        add([int(index) + 1 for index in np.flatnonzero(row)], 0)
    for row, rhs in zip(components["even_checks"], even_rhs):
        add(
            [symbols + int(index) + 1 for index in np.flatnonzero(row)],
            int(rhs),
        )

    m = problem["m"]
    physical_literals: list[int] = [0] * n
    for symbol in range(symbols):
        side, j = divmod(symbol, m)
        c_literal = symbol + 1
        e0_literal = symbols + symbol + 1
        e1_literal = 2 * symbols + symbol + 1
        coordinates = [
            side * 3 * m + j,
            side * 3 * m + m + j,
            side * 3 * m + 2 * m + j,
        ]
        p_literals = [n + coordinate + 1 for coordinate in coordinates]
        add([c_literal, e0_literal, p_literals[0]], 0)
        add([c_literal, e1_literal, p_literals[1]], 0)
        add([c_literal, e0_literal, e1_literal, p_literals[2]], 0)
        for coordinate, literal in zip(coordinates, p_literals):
            physical_literals[coordinate] = literal

    _add_class_symmetry_break(
        cnf, pool, physical_literals, representative, symmetry_break
    )
    _add_cardinality_constraint(
        cnf, pool, physical_literals, weight_cap, card_encoding
    )
    return cnf


def _verify_class_witness(
    problem: dict[str, Any],
    pole_cover: dict[str, Any],
    representative: dict[str, Any],
    vector: np.ndarray,
    weight_cap: int,
) -> dict[str, Any]:
    vector = np.asarray(vector, dtype=np.uint8) & 1
    syndrome = (
        pole_cover["x_pole"].astype(np.int64) @ vector[:, None].astype(np.int64) % 2
    )[:, 0].astype(np.uint8)
    kernel_np = not np.any(problem["HX"] @ vector % 2)
    syndrome_np = np.array_equal(syndrome, representative["syndrome"])
    outside_np = rank_np(np.vstack([problem["HZ"], vector])) == rank_np(problem["HZ"]) + 1
    bits = _vector_bitset(vector)
    kernel_bits = all(
        (row & bits).bit_count() % 2 == 0 for row in rows_to_bitsets(problem["HX"])
    )
    syndrome_bits = [
        (row & bits).bit_count() % 2 for row in rows_to_bitsets(pole_cover["x_pole"])
    ]
    syndrome_bitset = syndrome_bits == representative["syndrome_bits"]
    outside_bits = solve_bitset(
        rows_to_bitsets(problem["HZ"]), problem["n"], bits
    ) is None
    weight = int(vector.sum())
    valid = bool(
        kernel_np
        and kernel_bits
        and syndrome_np
        and syndrome_bitset
        and outside_np
        and outside_bits
        and weight <= int(weight_cap)
    )
    return {
        "kernel_numpy": bool(kernel_np),
        "kernel_bitset": bool(kernel_bits),
        "syndrome_numpy": bool(syndrome_np),
        "syndrome_bitset": bool(syndrome_bitset),
        "outside_z_stabilizer_numpy": bool(outside_np),
        "outside_z_stabilizer_bitset": bool(outside_bits),
        "weight": weight,
        "valid": valid,
    }


def decide_class_representative_cpsat(
    problem: dict[str, Any],
    pole_cover: dict[str, Any],
    representative: dict[str, Any],
    weight_cap: int,
    *,
    time_limit_s: float,
    workers: int = 1,
    formulation: str = "physical",
) -> dict[str, Any]:
    """Exact fixed-class feasibility model; no optimization objective."""
    from ortools.sat.python import cp_model

    model = cp_model.CpModel()

    def add_parity(selected: list[Any], rhs: int, tag: str) -> None:
        del tag
        if not selected:
            if rhs:
                model.add_bool_or([])
            return
        if rhs:
            model.add_bool_xor(selected)
        else:
            model.add_bool_xor([*selected, model.new_constant(1)])

    variables = [model.new_bool_var(f"v{index}") for index in range(problem["n"])]
    if formulation == "physical":
        for index, row in enumerate(problem["HX"]):
            add_parity(
                [variables[int(column)] for column in np.flatnonzero(row)],
                0,
                f"h{index}",
            )
        for index, (row, rhs) in enumerate(
            zip(pole_cover["x_pole"], representative["syndrome"])
        ):
            add_parity(
                [variables[int(column)] for column in np.flatnonzero(row)],
                int(rhs),
                f"l{index}",
            )
    elif formulation == "generator":
        generator = _sparse_independent_rows(problem["HZ"])
        coefficients = [
            model.new_bool_var(f"c{index}") for index in range(len(generator))
        ]
        for coordinate in range(problem["n"]):
            selected = [
                coefficients[int(index)]
                for index in np.flatnonzero(generator[:, coordinate])
            ]
            add_parity(
                [*selected, variables[coordinate]],
                int(representative["vector"][coordinate]),
                f"g{coordinate}",
            )
    elif formulation == "component":
        components = build_component_decomposition(problem, pole_cover)
        transformed = [
            model.new_bool_var(f"u{index}") for index in range(problem["n"])
        ]
        symbols = components["symbols"]
        transformed_rep = _component_transform(problem, representative["vector"])
        even_rep = transformed_rep[symbols:]
        even_rhs = (
            components["even_checks"].astype(np.int64)
            @ even_rep[:, None].astype(np.int64)
            % 2
        )[:, 0].astype(np.uint8)
        for index, row in enumerate(components["constant_checks"]):
            add_parity(
                [transformed[int(column)] for column in np.flatnonzero(row)],
                0,
                f"c{index}",
            )
        for index, (row, rhs) in enumerate(
            zip(components["even_checks"], even_rhs)
        ):
            add_parity(
                [
                    transformed[symbols + int(column)]
                    for column in np.flatnonzero(row)
                ],
                int(rhs),
                f"e{index}",
            )
        m = problem["m"]
        for symbol in range(symbols):
            side, j = divmod(symbol, m)
            c = transformed[symbol]
            e0 = transformed[symbols + symbol]
            e1 = transformed[2 * symbols + symbol]
            coordinates = [
                side * 3 * m + j,
                side * 3 * m + m + j,
                side * 3 * m + 2 * m + j,
            ]
            add_parity([c, e0, variables[coordinates[0]]], 0, f"p0_{symbol}")
            add_parity([c, e1, variables[coordinates[1]]], 0, f"p1_{symbol}")
            add_parity(
                [c, e0, e1, variables[coordinates[2]]], 0, f"p2_{symbol}"
            )
    else:
        raise ValueError("CP-SAT supports physical, generator, or component formulation")
    model.add(sum(variables) <= int(weight_cap))
    model.add_bool_or(
        [variables[int(index)] for index in representative["anchor_coordinates"]]
    )
    model_hash = sha256_bytes(str(model.proto).encode("utf-8"))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(time_limit_s)
    solver.parameters.num_workers = int(workers)
    started = time.perf_counter()
    status_code = solver.solve(model)
    wall = time.perf_counter() - started
    status_name = solver.status_name(status_code)
    if status_code in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        status = "SAT"
    elif status_code == cp_model.INFEASIBLE:
        status = "UNSAT"
    else:
        status = "UNDECIDED_BUDGET"
    output: dict[str, Any] = {
        "class_index": int(representative["class_index"]),
        "mask_hex": representative["mask_hex"],
        "orbit_size": int(representative["orbit_size"]),
        "status": status,
        "weight_cap": int(weight_cap),
        "instance_sha256": model_hash,
        "backend": "cpsat",
        "formulation": formulation,
        "solver": {
            "name": "OR-Tools CP-SAT",
            "backend": "cpsat",
            "status": status_name,
            "time_limit_s": float(time_limit_s),
            "workers": int(workers),
            "wall_time_s": wall,
            "conflicts": int(solver.num_conflicts),
            "branches": int(solver.num_branches),
        },
    }
    if status == "SAT":
        vector = np.asarray(
            [int(solver.value(variable)) for variable in variables], dtype=np.uint8
        )
        verification = _verify_class_witness(
            problem, pole_cover, representative, vector, weight_cap
        )
        if not verification["valid"]:
            raise RuntimeError("CP-SAT fixed-class model failed verification")
        output["vector"] = vector.astype(int).tolist()
        output["weight"] = verification["weight"]
        output["verification"] = verification
    return output


def decide_class_representative_cryptosat(
    problem: dict[str, Any],
    pole_cover: dict[str, Any],
    representative: dict[str, Any],
    weight_cap: int,
    *,
    formulation: str,
    card_encoding: str,
    symmetry_break: str,
    conflict_budget: int,
    time_limit_s: float,
) -> dict[str, Any]:
    """CryptoMiniSat with native XORs and a clausal cardinality encoding."""
    import pysat
    from pysat.formula import CNF, IDPool
    from pysat.solvers import Solver

    if card_encoding == "native":
        raise ValueError("CryptoMiniSat needs a clausal cardinality encoding")
    n = problem["n"]
    if formulation == "physical":
        output_offset = 0
        primary_variables = n
        xor_rows = [
            ([int(index) + 1 for index in np.flatnonzero(row)], 0)
            for row in problem["HX"]
        ]
        xor_rows.extend(
            (
                [int(index) + 1 for index in np.flatnonzero(row)],
                int(rhs),
            )
            for row, rhs in zip(pole_cover["x_pole"], representative["syndrome"])
        )
    elif formulation == "generator":
        generator = _sparse_independent_rows(problem["HZ"])
        output_offset = len(generator)
        primary_variables = output_offset + n
        xor_rows = []
        for coordinate in range(n):
            literals = [
                int(index) + 1
                for index in np.flatnonzero(generator[:, coordinate])
            ]
            literals.append(output_offset + coordinate + 1)
            xor_rows.append((literals, int(representative["vector"][coordinate])))
    else:
        raise ValueError("CryptoMiniSat supports physical or generator formulation")

    cnf = CNF()
    pool = IDPool(start_from=primary_variables + 1)
    physical_literals = list(range(output_offset + 1, output_offset + n + 1))
    _add_class_symmetry_break(
        cnf, pool, physical_literals, representative, symmetry_break
    )
    _add_cardinality_constraint(
        cnf, pool, physical_literals, weight_cap, card_encoding
    )
    fingerprint = canonical_json_sha256(
        {
            "cnf_sha256": cnf_sha256(cnf),
            "xors": xor_rows,
            "weight_cap": int(weight_cap),
        }
    )
    started = time.perf_counter()
    with Solver(name="cryptosat", bootstrap_with=cnf.clauses) as engine:
        for literals, rhs in xor_rows:
            engine.add_xor_clause(literals, value=bool(rhs))
        if conflict_budget > 0:
            engine.conf_budget(int(conflict_budget))
        if time_limit_s > 0:
            engine.time_budget(float(time_limit_s))
        answer = engine.solve_limited(expect_interrupt=False)
        model = engine.get_model() if answer else None
    wall = time.perf_counter() - started
    status = {True: "SAT", False: "UNSAT", None: "UNDECIDED_BUDGET"}[answer]
    output: dict[str, Any] = {
        "class_index": int(representative["class_index"]),
        "mask_hex": representative["mask_hex"],
        "orbit_size": int(representative["orbit_size"]),
        "status": status,
        "weight_cap": int(weight_cap),
        "instance_sha256": fingerprint,
        "backend": "cryptosat",
        "formulation": formulation,
        "symmetry_break": symmetry_break,
        "card_encoding": card_encoding,
        "solver": {
            "name": "PySAT CryptoMiniSat5",
            "version": pysat.__version__,
            "backend": "cryptosat",
            "conflict_budget": int(conflict_budget),
            "time_limit_s": float(time_limit_s),
            "clauses": len(cnf.clauses),
            "native_xors": len(xor_rows),
            "cnf_variables": cnf.nv,
            "wall_time_s": wall,
        },
    }
    if model is not None:
        assignment = {abs(literal): literal > 0 for literal in model}
        vector = np.asarray(
            [
                int(assignment.get(output_offset + index + 1, False))
                for index in range(n)
            ],
            dtype=np.uint8,
        )
        verification = _verify_class_witness(
            problem, pole_cover, representative, vector, weight_cap
        )
        if not verification["valid"]:
            raise RuntimeError("CryptoMiniSat fixed-class model failed verification")
        output["vector"] = vector.astype(int).tolist()
        output["weight"] = verification["weight"]
        output["verification"] = verification
    return output


def decide_class_representative_milp(
    problem: dict[str, Any],
    pole_cover: dict[str, Any],
    representative: dict[str, Any],
    weight_cap: int,
    *,
    time_limit_s: float,
) -> dict[str, Any]:
    """HiGHS exact-MIP feasibility control for the sparse coset generator."""
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import lil_matrix

    generator = _sparse_independent_rows(problem["HZ"])
    rank = len(generator)
    n = problem["n"]
    total = rank + n + n
    matrix = lil_matrix((n + 2, total), dtype=float)
    lower = np.empty(n + 2, dtype=float)
    upper = np.empty(n + 2, dtype=float)
    for coordinate in range(n):
        for index in np.flatnonzero(generator[:, coordinate]):
            matrix[coordinate, int(index)] = 1.0
        matrix[coordinate, rank + coordinate] = 1.0
        matrix[coordinate, rank + n + coordinate] = -2.0
        rhs = float(representative["vector"][coordinate])
        lower[coordinate] = rhs
        upper[coordinate] = rhs
    matrix[n, rank : rank + n] = 1.0
    lower[n], upper[n] = -np.inf, float(weight_cap)
    for coordinate in representative["anchor_coordinates"]:
        matrix[n + 1, rank + int(coordinate)] = 1.0
    lower[n + 1], upper[n + 1] = 1.0, np.inf
    matrix = matrix.tocsr()

    variable_lower = np.zeros(total, dtype=float)
    variable_upper = np.ones(total, dtype=float)
    variable_upper[rank + n :] = 2.0
    fingerprint = sha256_bytes(
        b"".join(
            [
                matrix.data.tobytes(),
                matrix.indices.tobytes(),
                matrix.indptr.tobytes(),
                lower.tobytes(),
                upper.tobytes(),
                variable_upper.tobytes(),
            ]
        )
    )
    started = time.perf_counter()
    result = milp(
        c=np.zeros(total, dtype=float),
        integrality=np.ones(total, dtype=np.uint8),
        bounds=Bounds(variable_lower, variable_upper),
        constraints=LinearConstraint(matrix, lower, upper),
        options={
            "time_limit": float(time_limit_s),
            "mip_rel_gap": 0.0,
            "presolve": True,
        },
    )
    wall = time.perf_counter() - started
    if result.status == 2:
        status = "UNSAT"
    elif result.status == 0 and result.x is not None:
        status = "SAT"
    else:
        status = "UNDECIDED_BUDGET"
    output: dict[str, Any] = {
        "class_index": int(representative["class_index"]),
        "mask_hex": representative["mask_hex"],
        "orbit_size": int(representative["orbit_size"]),
        "status": status,
        "weight_cap": int(weight_cap),
        "instance_sha256": fingerprint,
        "backend": "milp",
        "formulation": "generator",
        "solver": {
            "name": "SciPy milp / HiGHS",
            "backend": "milp",
            "status": int(result.status),
            "message": str(result.message),
            "time_limit_s": float(time_limit_s),
            "wall_time_s": wall,
            "nodes": int(getattr(result, "mip_node_count", 0) or 0),
            "gap": (
                None
                if getattr(result, "mip_gap", None) is None
                else float(result.mip_gap)
            ),
        },
    }
    if status == "SAT":
        vector = np.rint(result.x[rank : rank + n]).astype(np.uint8)
        verification = _verify_class_witness(
            problem, pole_cover, representative, vector, weight_cap
        )
        if not verification["valid"]:
            raise RuntimeError("HiGHS fixed-class model failed verification")
        output["vector"] = vector.astype(int).tolist()
        output["weight"] = verification["weight"]
        output["verification"] = verification
    return output


def decide_class_representative(
    problem: dict[str, Any],
    pole_cover: dict[str, Any],
    representative: dict[str, Any],
    weight_cap: int,
    *,
    conflict_budget: int = 0,
    solver_name: str = SOLVER_NAME,
    xor_encoding: str = "sequential",
    card_encoding: str = "seqcounter",
    formulation: str = "physical",
    components: dict[str, Any] | None = None,
    backend: str = "pysat",
    time_limit_s: float = 60.0,
    backend_workers: int = 1,
    symmetry_break: str = "anchor",
) -> dict[str, Any]:
    """Decide one fixed-class bounded-weight query with CDCL."""
    import pysat
    from pysat.solvers import Solver
    if backend in {"cpsat", "milp"} and symmetry_break != "anchor":
        raise ValueError("lex symmetry breaking is implemented only for PySAT")
    if backend == "cpsat":
        if formulation not in {"physical", "generator", "component"}:
            raise ValueError(
                "CP-SAT control supports physical, generator, or component formulation"
            )
        return decide_class_representative_cpsat(
            problem,
            pole_cover,
            representative,
            weight_cap,
            time_limit_s=time_limit_s,
            formulation=formulation,
            workers=backend_workers,
        )
    if backend == "milp":
        if formulation != "generator":
            raise ValueError("MILP control currently supports generator formulation only")
        return decide_class_representative_milp(
            problem,
            pole_cover,
            representative,
            weight_cap,
            time_limit_s=time_limit_s,
        )
    if backend == "cryptosat":
        return decide_class_representative_cryptosat(
            problem,
            pole_cover,
            representative,
            weight_cap,
            formulation=formulation,
            card_encoding=card_encoding,
            symmetry_break=symmetry_break,
            conflict_budget=conflict_budget,
            time_limit_s=time_limit_s,
        )
    if backend != "pysat":
        raise ValueError(f"unknown fixed-class backend {backend!r}")

    generator = None
    if formulation == "physical":
        cnf = build_class_cnf(
            problem, pole_cover, representative, weight_cap,
            xor_encoding=xor_encoding, card_encoding=card_encoding,
            symmetry_break=symmetry_break,
        )
    elif formulation == "component":
        components = components or build_component_decomposition(problem, pole_cover)
        cnf = build_component_class_cnf(
            problem, components, representative, weight_cap,
            xor_encoding=xor_encoding, card_encoding=card_encoding,
            symmetry_break=symmetry_break,
        )
    elif formulation == "generator":
        cnf, generator = build_generator_class_cnf(
            problem, representative, weight_cap,
            xor_encoding=xor_encoding, card_encoding=card_encoding,
            symmetry_break=symmetry_break,
        )
    else:
        raise ValueError(f"unknown fixed-class formulation {formulation!r}")
    if formulation == "component":
        output_offset = problem["n"]
    elif formulation == "generator":
        assert generator is not None
        output_offset = len(generator)
    else:
        output_offset = 0
    native_weight_literals = list(
        range(output_offset + 1, output_offset + problem["n"] + 1)
    )
    digest = (
        canonical_json_sha256(
            {
                "cnf_sha256": cnf_sha256(cnf),
                "native_atmost": [native_weight_literals, int(weight_cap)],
            }
        )
        if card_encoding == "native"
        else cnf_sha256(cnf)
    )
    if card_encoding == "native" and solver_name != "minicard":
        raise ValueError("native cardinality requires solver_name='minicard'")
    started = time.perf_counter()
    with Solver(name=solver_name, bootstrap_with=cnf.clauses) as engine:
        if card_encoding == "native":
            engine.add_atmost(native_weight_literals, int(weight_cap))
        if conflict_budget > 0:
            engine.conf_budget(int(conflict_budget))
            answer = engine.solve_limited(expect_interrupt=False)
        else:
            answer = engine.solve()
        try:
            stats = dict(engine.accum_stats())
        except NotImplementedError:
            stats = {}
        model = engine.get_model() if answer else None
    wall = time.perf_counter() - started
    status = {True: "SAT", False: "UNSAT", None: "UNDECIDED_BUDGET"}[answer]
    output: dict[str, Any] = {
        "class_index": int(representative["class_index"]),
        "mask_hex": representative["mask_hex"],
        "orbit_size": int(representative["orbit_size"]),
        "status": status,
        "weight_cap": int(weight_cap),
        "cnf_sha256": digest,
        "instance_sha256": digest,
        "backend": "pysat",
        "xor_encoding": xor_encoding,
        "card_encoding": card_encoding,
        "formulation": formulation,
        "symmetry_break": symmetry_break,
        "solver": {
            "name": f"PySAT {solver_name}",
            "version": pysat.__version__,
            "backend": "pysat",
            "conflict_budget": int(conflict_budget),
            "clauses": len(cnf.clauses),
            "cnf_variables": cnf.nv,
            "wall_time_s": wall,
            "stats": {key: int(value) for key, value in stats.items()},
        },
    }
    if model is not None:
        assignment = {abs(literal): literal > 0 for literal in model}
        offset = output_offset
        vector = np.asarray(
            [
                int(assignment.get(offset + index + 1, False))
                for index in range(problem["n"])
            ],
            dtype=np.uint8,
        )
        verification = _verify_class_witness(
            problem, pole_cover, representative, vector, weight_cap
        )
        if not verification["valid"]:
            raise RuntimeError("fixed-class SAT model failed independent verification")
        output["vector"] = vector.astype(int).tolist()
        output["weight"] = verification["weight"]
        output["verification"] = verification
    return output


def decide_class_orbit_bundle(
    problem: dict[str, Any],
    pole_cover: dict[str, Any],
    class_cover: dict[str, Any],
    weight_cap: int,
    *,
    conflict_budget: int = 0,
    solver_name: str = SOLVER_NAME,
    representative_indexes: list[int] | None = None,
    xor_encoding: str = "sequential",
    card_encoding: str = "seqcounter",
    formulation: str = "physical",
    backend: str = "pysat",
    time_limit_s: float = 60.0,
    backend_workers: int = 1,
    symmetry_break: str = "anchor",
) -> dict[str, Any]:
    """Decide the class-orbit partition; strict subsets never certify UNSAT."""
    if not pole_cover.get("complete") or not class_cover.get("complete"):
        raise RuntimeError("refusing a decision from an incomplete class cover")
    total = int(class_cover["class_orbits"])
    selected = (
        list(range(total))
        if representative_indexes is None
        else sorted(set(map(int, representative_indexes)))
    )
    if any(index < 0 or index >= total for index in selected):
        raise ValueError("class representative index outside cover")
    records: list[dict[str, Any]] = []
    witness: dict[str, Any] | None = None
    undecided = False
    started = time.perf_counter()
    components = (
        build_component_decomposition(problem, pole_cover)
        if formulation == "component"
        else None
    )
    for index in selected:
        record = decide_class_representative(
            problem,
            pole_cover,
            class_cover["representatives"][index],
            weight_cap,
            conflict_budget=conflict_budget,
            solver_name=solver_name,
            xor_encoding=xor_encoding,
            card_encoding=card_encoding,
            formulation=formulation,
            components=components,
            backend=backend,
            time_limit_s=time_limit_s,
            backend_workers=backend_workers,
            symmetry_break=symmetry_break,
        )
        records.append(record)
        if record["status"] == "SAT":
            witness = record
            break
        if record["status"] == "UNDECIDED_BUDGET":
            undecided = True
    full_coverage = selected == list(range(total))
    if witness is not None:
        status = "SAT"
    elif undecided:
        status = "UNDECIDED_BUDGET"
    elif full_coverage:
        status = "UNSAT"
    else:
        status = "UNSAT_SUBSET"
    output: dict[str, Any] = {
        "status": status,
        "weight_cap": int(weight_cap),
        "required_class_orbits": total,
        "selected_class_orbits": selected,
        "class_orbits_solved": len(records),
        "records": records,
        "cnf_bundle_sha256": canonical_json_sha256(
            {
                "schema": "exp056-class-orbit-bundle-v1",
                "weight_cap": int(weight_cap),
                "xor_encoding": xor_encoding,
                "card_encoding": card_encoding,
                "formulation": formulation,
                "backend": backend,
                "time_limit_s": float(time_limit_s),
                "symmetry_break": symmetry_break,
                "component_decomposition": (
                    canonical_json_sha256(_component_metadata(components))
                    if components is not None
                    else None
                ),
                "class_cover": canonical_json_sha256(_class_cover_metadata(class_cover)),
                "selected": selected,
                "hashes": [
                    record.get("instance_sha256", record.get("cnf_sha256"))
                    for record in records
                ],
            }
        ),
        "wall_time_s": time.perf_counter() - started,
    }
    if witness is not None:
        output.update(
            {
                "class_index": witness["class_index"],
                "vector": witness["vector"],
                "weight": witness["weight"],
                "verification": witness["verification"],
            }
        )
    return output


def decide_orbit_bundle(
    problem: dict[str, Any],
    cover: dict[str, Any],
    weight_cap: int,
    *,
    conflict_budget: int = 0,
    solver_name: str = SOLVER_NAME,
) -> dict[str, Any]:
    """Decide one bounded Z-distance question using all orbit generators."""
    if not cover.get("complete"):
        raise RuntimeError("refusing a decision from an incomplete orbit cover")
    records: list[dict[str, Any]] = []
    witness_record: dict[str, Any] | None = None
    undecided = False
    started = time.perf_counter()
    for sector in cover["sectors"]:
        instance = _sector_instance(problem, sector)
        record = decide_weight_bounded(
            instance,
            int(weight_cap),
            solver_name=solver_name,
            conflict_budget=int(conflict_budget),
        )
        records.append(record)
        if record["status"] == "SAT":
            vector = np.asarray(record["vector"], dtype=np.uint8)
            verification = _verify_sector_witness(problem, instance, vector)
            if not verification["valid"]:
                raise RuntimeError("orbit-sector SAT witness failed CSS verification")
            witness_record = {
                "vector": record["vector"],
                "weight": group_weight(instance, vector),
                "verification": verification,
                "sector": sector["sector"],
            }
            break
        if record["status"] == "UNDECIDED_BUDGET":
            undecided = True
    if witness_record is not None:
        status = "SAT"
    elif undecided:
        status = "UNDECIDED_BUDGET"
    else:
        status = "UNSAT"
    hashes = [record["cnf_sha256"] for record in records]
    output: dict[str, Any] = {
        "status": status,
        "weight_cap": int(weight_cap),
        "required_sectors": int(cover["required_sectors"]),
        "sectors_solved": len(records),
        "records": records,
        "cnf_bundle_sha256": canonical_json_sha256(
            {
                "schema": SECTOR_SCHEMA,
                "weight_cap": int(weight_cap),
                "cover": canonical_json_sha256(_cover_metadata(cover)),
                "hashes": hashes,
            }
        ),
        "wall_time_s": time.perf_counter() - started,
    }
    if witness_record is not None:
        output.update(witness_record)
    return output


def _problem_identity(problem: dict[str, Any], cover: dict[str, Any], cap: int) -> dict[str, Any]:
    return {
        "schema": SECTOR_SCHEMA,
        "name": problem["record"]["name"],
        "ell": problem["ell"],
        "m": problem["m"],
        "n": problem["n"],
        "k": problem["k"],
        "A": [list(term) for term in problem["A_terms"]],
        "B": [list(term) for term in problem["B_terms"]],
        "HX_sha256": matrix_sha256(problem["HX"]),
        "HZ_sha256": matrix_sha256(problem["HZ"]),
        "cover_sha256": canonical_json_sha256(_cover_metadata(cover)),
        "weight_cap": int(cap),
    }


def _sector_path(sector: int) -> Path:
    return STATE_DIR / f"sector_{sector:02d}.json"


def run_sector(
    problem: dict[str, Any],
    cover: dict[str, Any],
    sector_index: int,
    cap: int,
    *,
    conflict_budget: int,
    solver_name: str,
    force: bool = False,
) -> dict[str, Any]:
    sector = cover["sectors"][sector_index]
    identity = {**_problem_identity(problem, cover, cap), "sector": sector_index}
    path = _sector_path(sector_index)
    if path.exists() and not force:
        stored = json.loads(path.read_text(encoding="utf-8"))
        if stored.get("identity") == identity and stored.get("decision", {}).get("status") != "UNDECIDED_BUDGET":
            return stored
    instance = _sector_instance(problem, sector)
    decision = decide_weight_bounded(
        instance,
        cap,
        conflict_budget=conflict_budget,
        solver_name=solver_name,
    )
    if decision["status"] == "SAT":
        vector = np.asarray(decision["vector"], dtype=np.uint8)
        decision["css_verification"] = _verify_sector_witness(problem, instance, vector)
        if not decision["css_verification"]["valid"]:
            raise RuntimeError("persisted sector witness failed verification")
    payload = {
        "identity": identity,
        "utc": utc_now(),
        "sector": {
            key: value for key, value in sector.items() if key not in {"detector", "orbit"}
        }
        | {"orbit_sha256": matrix_sha256(sector["orbit"])},
        "decision": decision,
    }
    atomic_write_json(path, payload)
    return payload


def _validate_sector_record(
    problem: dict[str, Any], cover: dict[str, Any], cap: int, sector_index: int
) -> dict[str, Any] | None:
    path = _sector_path(sector_index)
    if not path.exists():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    expected = {**_problem_identity(problem, cover, cap), "sector": sector_index}
    if payload.get("identity") != expected:
        return None
    instance = _sector_instance(problem, cover["sectors"][sector_index])
    decision = payload.get("decision", {})
    if decision.get("cnf_sha256") != decision_cnf_digest(instance, cap):
        return None
    if decision.get("status") == "SAT":
        vector = np.asarray(decision.get("vector", []), dtype=np.uint8)
        if vector.shape != (problem["n"],):
            return None
        if not _verify_sector_witness(problem, instance, vector)["valid"]:
            return None
    return payload


def _class_protocol(
    problem: dict[str, Any],
    cover: dict[str, Any],
    class_cover: dict[str, Any],
) -> dict[str, Any]:
    parity = kernel_weight_parity_proof(problem, int(TARGET["lower_cap"]))
    if not parity["valid"]:
        raise RuntimeError("class protocol requires the even-kernel parity proof")
    return {
        "schema": "exp056-class-certificate-protocol-v1",
        "encoding_version": "exp056-sparse-generator-direct-small-v1",
        "requested_exclusion_cap": int(TARGET["lower_cap"]),
        "effective_solver_cap": int(parity["effective_even_cap"]),
        "parity_proof": parity,
        "solver": "kissat404",
        "backend": "pysat",
        "formulation": "generator",
        "xor_encoding": "direct-small",
        "card_encoding": "seqcounter",
        "symmetry_break": "anchor",
        "conflict_budget": 0,
        "cover_sha256": canonical_json_sha256(_cover_metadata(cover)),
        "class_cover_sha256": canonical_json_sha256(
            _class_cover_metadata(class_cover)
        ),
    }


def _class_instance_digest(
    problem: dict[str, Any],
    representative: dict[str, Any],
    protocol: dict[str, Any],
) -> str:
    cnf, _ = build_generator_class_cnf(
        problem,
        representative,
        int(protocol["effective_solver_cap"]),
        xor_encoding=protocol["xor_encoding"],
        card_encoding=protocol["card_encoding"],
        symmetry_break=protocol["symmetry_break"],
    )
    return cnf_sha256(cnf)


def _class_identity(
    problem: dict[str, Any],
    cover: dict[str, Any],
    class_cover: dict[str, Any],
    class_index: int,
) -> dict[str, Any]:
    protocol = _class_protocol(problem, cover, class_cover)
    representative = class_cover["representatives"][class_index]
    return {
        "schema": "exp056-class-record-v1",
        "name": problem["record"]["name"],
        "ell": problem["ell"],
        "m": problem["m"],
        "n": problem["n"],
        "k": problem["k"],
        "A": [list(term) for term in problem["A_terms"]],
        "B": [list(term) for term in problem["B_terms"]],
        "HX_sha256": matrix_sha256(problem["HX"]),
        "HZ_sha256": matrix_sha256(problem["HZ"]),
        "protocol_sha256": canonical_json_sha256(protocol),
        "class_index": int(class_index),
        "class_mask_hex": representative["mask_hex"],
        "class_vector_sha256": representative["vector_sha256"],
        "instance_sha256": _class_instance_digest(
            problem, representative, protocol
        ),
    }


def _class_path(class_index: int) -> Path:
    return CLASS_STATE_DIR / f"class_{class_index:02d}.json"


def _validate_class_record(
    problem: dict[str, Any],
    cover: dict[str, Any],
    class_cover: dict[str, Any],
    class_index: int,
) -> dict[str, Any] | None:
    path = _class_path(class_index)
    if not path.exists():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    identity = _class_identity(problem, cover, class_cover, class_index)
    if payload.get("identity") != identity:
        return None
    decision = payload.get("decision", {})
    if decision.get("instance_sha256") != identity["instance_sha256"]:
        return None
    if decision.get("status") not in {"SAT", "UNSAT", "UNDECIDED_BUDGET"}:
        return None
    if decision.get("status") == "SAT":
        vector = np.asarray(decision.get("vector", []), dtype=np.uint8)
        if vector.shape != (problem["n"],):
            return None
        representative = class_cover["representatives"][class_index]
        if not _verify_class_witness(
            problem,
            cover,
            representative,
            vector,
            int(_class_protocol(problem, cover, class_cover)["effective_solver_cap"]),
        )["valid"]:
            return None
    return payload


def _class_replay_valid(record: dict[str, Any] | None) -> bool:
    if record is None or record.get("decision", {}).get("status") != "UNSAT":
        return False
    replay = record.get("replay", {})
    identity = record.get("identity", {})
    return bool(
        replay.get("status") == "UNSAT"
        and replay.get("instance_sha256") == identity.get("instance_sha256")
        and replay.get("encoding_version")
        == "exp056-sparse-generator-direct-small-v1"
    )


def run_class_record(
    problem: dict[str, Any],
    cover: dict[str, Any],
    class_cover: dict[str, Any],
    class_index: int,
    *,
    force: bool = False,
    replay: bool = False,
) -> dict[str, Any]:
    protocol = _class_protocol(problem, cover, class_cover)
    identity = _class_identity(problem, cover, class_cover, class_index)
    path = _class_path(class_index)
    stored = _validate_class_record(
        problem, cover, class_cover, class_index
    )
    if replay:
        if stored is None or stored["decision"]["status"] != "UNSAT":
            raise RuntimeError(
                f"class {class_index}: replay requires a validated UNSAT record"
            )
    elif stored is not None and not force:
        return stored

    representative = class_cover["representatives"][class_index]
    decision = decide_class_representative(
        problem,
        cover,
        representative,
        int(protocol["effective_solver_cap"]),
        conflict_budget=0,
        solver_name=protocol["solver"],
        xor_encoding=protocol["xor_encoding"],
        card_encoding=protocol["card_encoding"],
        formulation=protocol["formulation"],
        backend=protocol["backend"],
        symmetry_break=protocol["symmetry_break"],
    )
    if decision.get("instance_sha256") != identity["instance_sha256"]:
        raise RuntimeError(f"class {class_index}: rebuilt CNF digest drifted")
    if replay:
        if decision["status"] != "UNSAT":
            raise RuntimeError(
                f"class {class_index}: replay returned {decision['status']}"
            )
        assert stored is not None
        stored["replay"] = {
            "status": "UNSAT",
            "instance_sha256": decision["instance_sha256"],
            "encoding_version": protocol["encoding_version"],
            "solver": decision["solver"],
            "utc": utc_now(),
        }
        atomic_write_json(path, stored)
        return stored

    payload = {
        "identity": identity,
        "protocol": protocol,
        "representative": {
            key: value
            for key, value in representative.items()
            if key
            not in {
                "vector",
                "syndrome",
                "orbit_masks",
                "stabilizer_mappings",
            }
        },
        "decision": decision,
        "utc": utc_now(),
    }
    atomic_write_json(path, payload)
    if decision["status"] == "SAT":
        atomic_write_json(
            COUNTEREXAMPLE,
            {
                "schema": "exp056-distance-counterexample-v1",
                "utc": utc_now(),
                "target": {
                    key: value
                    for key, value in TARGET.items()
                    if key != "witness_support"
                },
                "identity": identity,
                "decision": decision,
            },
        )
    return payload


def assemble(*, require_exact: bool = False) -> dict[str, Any]:
    problem = build_problem(TARGET)
    cover = build_orbit_cover(problem)
    class_cover = build_class_orbit_cover(problem, cover)
    class_protocol = _class_protocol(problem, cover, class_cover)
    duality = bb_duality_proof(problem)
    upper = verify_upper_witness(problem, default_witness(problem))

    sector_cap = int(TARGET["lower_cap"])
    sector_records = [
        _validate_sector_record(problem, cover, sector_cap, sector)
        for sector in range(cover["required_sectors"])
    ]
    sector_missing = [
        index for index, record in enumerate(sector_records) if record is None
    ]
    sector_statuses = [
        record["decision"]["status"]
        for record in sector_records
        if record is not None
    ]
    sector_all_unsat = bool(
        not sector_missing
        and sector_statuses == ["UNSAT"] * cover["required_sectors"]
    )
    sector_replayed = bool(
        sector_all_unsat
        and all(record.get("replay", {}).get("status") == "UNSAT"
                for record in sector_records if record is not None)
    )

    class_records = [
        _validate_class_record(problem, cover, class_cover, class_index)
        for class_index in range(class_cover["class_orbits"])
    ]
    class_missing = [
        index for index, record in enumerate(class_records) if record is None
    ]
    class_statuses = [
        record["decision"]["status"]
        for record in class_records
        if record is not None
    ]
    class_all_unsat = bool(
        not class_missing
        and class_statuses == ["UNSAT"] * class_cover["class_orbits"]
    )
    class_replayed = bool(
        class_all_unsat and all(_class_replay_valid(record) for record in class_records)
    )
    counterexample = any(
        status == "SAT" for status in [*sector_statuses, *class_statuses]
    )
    lower_exact = bool(sector_replayed or class_replayed)
    expected_d = int(TARGET["expected_d"])
    excluded_through = int(class_protocol["requested_exclusion_cap"])
    exact = bool(
        not counterexample
        and lower_exact
        and upper["weight"] == expected_d
        and excluded_through + 1 == expected_d
        and duality["valid"]
        and cover["complete"]
        and class_cover["complete"]
    )
    if exact:
        classification = "CERTIFIED_EXACT"
    elif counterexample:
        classification = "COUNTEREXAMPLE"
    elif class_all_unsat or sector_all_unsat:
        classification = "REPLAY_PENDING"
    else:
        classification = "PARTIAL"

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "utc": utc_now(),
        "target": {
            key: value for key, value in TARGET.items() if key != "witness_support"
        },
        "identity": _problem_identity(problem, cover, sector_cap),
        "source_distance_was_estimate": True,
        "orbit_cover": _cover_metadata(cover),
        "class_cover": _class_cover_metadata(class_cover),
        "duality": {
            key: value for key, value in duality.items() if key != "permutation"
        },
        "upper_bound": upper,
        "lower_bound": {
            "excluded_through": excluded_through,
            "class_route": {
                "protocol": class_protocol,
                "required_classes": class_cover["class_orbits"],
                "missing_classes": class_missing,
                "statuses": class_statuses,
                "all_classes_unsat": class_all_unsat,
                "all_classes_replayed": class_replayed,
                "records": class_records,
            },
            "pole_sector_route": {
                "weight_cap": sector_cap,
                "required_sectors": cover["required_sectors"],
                "missing_sectors": sector_missing,
                "statuses": sector_statuses,
                "all_orbit_sectors_unsat": sector_all_unsat,
                "all_orbit_sectors_replayed": sector_replayed,
                "records": sector_records,
            },
        },
        "verdict": {
            "exact": exact,
            "d": expected_d if exact else None,
            "d_X": expected_d if exact else None,
            "d_Z": expected_d if exact else None,
            "classification": classification,
            "lower_route": (
                "class_orbits" if class_replayed
                else ("pole_sectors" if sector_replayed else None)
            ),
        },
    }
    atomic_write_json(CERTIFICATE if exact else PARTIAL, payload)
    if require_exact and not exact:
        raise RuntimeError(
            "EXP-056 remains noncanonical: "
            f"classification={classification}, class_missing={class_missing}"
        )
    return payload


def _parse_sectors(text: str, count: int) -> list[int]:
    if text == "all":
        return list(range(count))
    sectors = sorted({int(token) for token in text.split(",") if token})
    if any(sector < 0 or sector >= count for sector in sectors):
        raise ValueError(f"sector selection {sectors} outside 0..{count - 1}")
    return sectors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    run_parser = sub.add_parser("run")
    run_parser.add_argument("--sectors", default="all")
    run_parser.add_argument("--conflict-budget", type=int, default=0)
    run_parser.add_argument("--solver", default=SOLVER_NAME)
    run_parser.add_argument("--force", action="store_true")
    class_run_parser = sub.add_parser("run-classes")
    class_run_parser.add_argument("--indexes", default="all")
    class_run_parser.add_argument("--force", action="store_true")
    class_run_parser.add_argument("--replay", action="store_true")
    sub.add_parser("assemble")
    probe_parser = sub.add_parser("probe")
    probe_parser.add_argument("--cap", type=int, default=int(TARGET["lower_cap"]))
    probe_parser.add_argument("--conflict-budget", type=int, default=1_000_000)
    probe_parser.add_argument("--solver", default=SOLVER_NAME)
    monolithic_parser = sub.add_parser("probe-monolithic")
    monolithic_parser.add_argument("--cap", type=int, default=int(TARGET["lower_cap"]))
    monolithic_parser.add_argument("--conflict-budget", type=int, default=1_000_000)
    monolithic_parser.add_argument("--solver", default=SOLVER_NAME)
    monolithic_parser.add_argument(
        "--xor-encoding", choices=("sequential", "balanced"), default="sequential"
    )
    monolithic_parser.add_argument(
        "--card-encoding",
        choices=("seqcounter", "totalizer", "cardnetwrk", "mtotalizer", "kmtotalizer", "native"),
        default="seqcounter",
    )
    class_probe_parser = sub.add_parser("probe-classes")
    class_probe_parser.add_argument("--cap", type=int, default=int(TARGET["lower_cap"]))
    class_probe_parser.add_argument("--conflict-budget", type=int, default=100_000)
    class_probe_parser.add_argument("--solver", default=SOLVER_NAME)
    class_probe_parser.add_argument("--indexes", default="all")
    class_probe_parser.add_argument(
        "--xor-encoding",
        choices=("sequential", "balanced", "direct-small"),
        default="sequential",
    )
    class_probe_parser.add_argument(
        "--card-encoding",
        choices=("seqcounter", "totalizer", "cardnetwrk", "mtotalizer", "kmtotalizer", "native"),
        default="seqcounter",
    )
    class_probe_parser.add_argument(
        "--formulation",
        choices=("physical", "component", "generator"),
        default="physical",
    )
    class_probe_parser.add_argument(
        "--symmetry-break",
        choices=("anchor", "lex", "anchor+lex"),
        default="anchor",
    )
    class_probe_parser.add_argument(
        "--backend",
        choices=("pysat", "cpsat", "milp", "cryptosat"),
        default="pysat",
    )
    class_probe_parser.add_argument("--time-limit", type=float, default=60.0)
    class_probe_parser.add_argument("--backend-workers", type=int, default=1)
    args = parser.parse_args()

    if args.command == "assemble":
        payload = assemble(require_exact=False)
        print(json.dumps(payload["verdict"], indent=2))
        return 0 if payload["verdict"]["exact"] else 2

    problem = build_problem(TARGET)
    cover = build_orbit_cover(problem)
    if args.command == "run-classes":
        class_cover = build_class_orbit_cover(problem, cover)
        indexes = _parse_sectors(args.indexes, class_cover["class_orbits"])
        statuses: list[str] = []
        for class_index in indexes:
            record = run_class_record(
                problem,
                cover,
                class_cover,
                class_index,
                force=args.force,
                replay=args.replay,
            )
            status = record["decision"]["status"]
            statuses.append(status)
            print(
                json.dumps(
                    {
                        "class_index": class_index,
                        "status": status,
                        "replay_verified": _class_replay_valid(record),
                        "wall_time_s": record["decision"]["solver"]["wall_time_s"],
                    }
                ),
                flush=True,
            )
        payload = assemble(require_exact=False)
        print(json.dumps(payload["verdict"], indent=2))
        return 0 if all(status == "UNSAT" for status in statuses) else 2
    if args.command == "probe-monolithic":
        result = decide_pole_monolithic(
            problem,
            cover,
            args.cap,
            conflict_budget=args.conflict_budget,
            solver_name=args.solver,
            xor_encoding=args.xor_encoding,
            card_encoding=args.card_encoding,
        )
        print(
            json.dumps(
                {
                    "status": result["status"],
                    "weight_cap": result["weight_cap"],
                    "wall_time_s": result["solver"]["wall_time_s"],
                    "conflicts": result["solver"]["stats"].get("conflicts"),
                    "cnf_sha256": result["cnf_sha256"],
                },
                indent=2,
            )
        )
        return 0 if result["status"] != "UNDECIDED_BUDGET" else 2
    if args.command == "probe":
        result = decide_orbit_bundle(
            problem,
            cover,
            args.cap,
            conflict_budget=args.conflict_budget,
            solver_name=args.solver,
        )
        print(json.dumps({key: value for key, value in result.items() if key != "records"}, indent=2))
        print(json.dumps([
            {
                "sector": index,
                "status": record["status"],
                "wall_time_s": record["solver"]["wall_time_s"],
                "conflicts": record["solver"]["stats"].get("conflicts"),
            }
            for index, record in enumerate(result["records"])
        ], indent=2))
        return 0 if result["status"] != "UNDECIDED_BUDGET" else 2
    if args.command == "probe-classes":
        class_cover = build_class_orbit_cover(problem, cover)
        indexes = (
            None
            if args.indexes == "all"
            else [int(token) for token in args.indexes.split(",") if token]
        )
        result = decide_class_orbit_bundle(
            problem,
            cover,
            class_cover,
            args.cap,
            conflict_budget=args.conflict_budget,
            solver_name=args.solver,
            representative_indexes=indexes,
            xor_encoding=args.xor_encoding,
            card_encoding=args.card_encoding,
            formulation=args.formulation,
            backend=args.backend,
            time_limit_s=args.time_limit,
            backend_workers=args.backend_workers,
            symmetry_break=args.symmetry_break,
        )
        print(
            json.dumps(
                {key: value for key, value in result.items() if key != "records"},
                indent=2,
            )
        )
        print(
            json.dumps(
                [
                    {
                        "class_index": record["class_index"],
                        "orbit_size": record["orbit_size"],
                        "status": record["status"],
                        "wall_time_s": record["solver"]["wall_time_s"],
                        "conflicts": record["solver"].get(
                            "stats", {}
                        ).get("conflicts", record["solver"].get("conflicts")),
                    }
                    for record in result["records"]
                ],
                indent=2,
            )
        )
        return 0 if result["status"] != "UNDECIDED_BUDGET" else 2

    cap = int(TARGET["lower_cap"])
    for sector in _parse_sectors(args.sectors, cover["required_sectors"]):
        record = run_sector(
            problem,
            cover,
            sector,
            cap,
            conflict_budget=args.conflict_budget,
            solver_name=args.solver,
            force=args.force,
        )
        decision = record["decision"]
        print(
            json.dumps(
                {
                    "sector": sector,
                    "status": decision["status"],
                    "wall_time_s": decision["solver"]["wall_time_s"],
                    "conflicts": decision["solver"]["stats"].get("conflicts"),
                    "cnf_sha256": decision["cnf_sha256"],
                }
            ),
            flush=True,
        )
    payload = assemble(require_exact=False)
    print(json.dumps(payload["verdict"], indent=2))
    return 0 if payload["verdict"]["exact"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
