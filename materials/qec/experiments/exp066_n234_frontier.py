"""EXP-066: fail-closed odd-lattice frontier screen at n=234.

The bounded target is exactly the two nonempty weight-3 frontiers (13,9) and
(39,3).  The existing EXP-055 screen is run first.  Residual records receive a
deterministic, deeper information-set reduction of all 255 pole logical
classes; every successful reduction is accepted only after rebuilding the CSS
matrices and physically checking the witness. Residuals are grouped by the
complete automorphism action instead of multiplying an exact solver per class.

For (39,3), the script also audits the full group automorphism set

    Z_39 x Z_3 ~= Z_13 x F_3^2,
    Aut ~= F_13^* x GL(2,3),

rather than only the diagonal unit subgroup used by the census. This exact
576-map action partitions the initial hard set into exact-decision bundles.
Every transported verdict requires a matrix-verified qubit permutation and a
fresh physical witness check against the target's rebuilt matrices.

Run:
  uv run python experiments/exp066_n234_frontier.py run --force
  uv run python experiments/exp066_n234_frontier.py validate
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "experiments" / filename
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


E55 = _load("exp055_for_exp066", "exp055_odd_lattice_sweep.py")
E56 = _load("exp056_for_exp066", "exp056_odd_distance.py")

from qec_research.gf2.linalg import rank_np  # noqa: E402

SCHEMA = "exp066-n234-frontier-v1"
LATTICES = ((13, 9), (39, 3))
K_MIN = 8
K_MAX = 24
TIME_LIMIT_S = 120.0
DEEP_TRIES = 2_000
DEEP_WORKERS = 8
SUMMARY = ROOT / "results" / "processed" / "exp066_n234_frontier.json"


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _crt_13_3(residue_13: int, residue_3: int) -> int:
    """CRT inverse for Z_39 ~= Z_13 x Z_3."""
    residue_13 %= 13
    return (residue_13 + 13 * ((residue_3 - residue_13) % 3)) % 39


def _gl2_3() -> tuple[tuple[int, int, int, int], ...]:
    matrices = tuple(
        (a, b, c, d)
        for a in range(3)
        for b in range(3)
        for c in range(3)
        for d in range(3)
        if (a * d - b * c) % 3
    )
    if len(matrices) != 48:
        raise RuntimeError("GL(2,3) enumeration is incomplete")
    return matrices


def _map_coordinate(
    coordinate: tuple[int, int],
    unit_13: int,
    matrix: tuple[int, int, int, int],
) -> tuple[int, int]:
    x, y = coordinate
    a, b, c, d = matrix
    x13 = unit_13 * (x % 13) % 13
    x3 = (a * (x % 3) + b * y) % 3
    mapped_y = (c * (x % 3) + d * y) % 3
    return _crt_13_3(x13, x3), mapped_y


@lru_cache(maxsize=1)
def z39_z3_automorphisms() -> tuple[dict[str, Any], ...]:
    """Enumerate and machine-check all 576 automorphisms of Z_39 x Z_3."""
    automorphisms: list[dict[str, Any]] = []
    unique: set[bytes] = set()
    for unit_13 in range(1, 13):
        for matrix in _gl2_3():
            mapping = np.asarray(
                [
                    3 * mapped_x + mapped_y
                    for x in range(39)
                    for y in range(3)
                    for mapped_x, mapped_y in [
                        _map_coordinate((x, y), unit_13, matrix)
                    ]
                ],
                dtype=np.int16,
            )
            bijective = bool(
                sorted(mapping.astype(int).tolist()) == list(range(117))
            )
            generator_images = (
                _map_coordinate((1, 0), unit_13, matrix),
                _map_coordinate((0, 1), unit_13, matrix),
            )
            homomorphism = True
            for x in range(39):
                for y in range(3):
                    image = _map_coordinate((x, y), unit_13, matrix)
                    for generator, generator_image in zip(
                        ((1, 0), (0, 1)), generator_images
                    ):
                        shifted = ((x + generator[0]) % 39, (y + generator[1]) % 3)
                        expected = (
                            (image[0] + generator_image[0]) % 39,
                            (image[1] + generator_image[1]) % 3,
                        )
                        if (
                            _map_coordinate(shifted, unit_13, matrix)
                            != expected
                        ):
                            homomorphism = False
                            break
                    if not homomorphism:
                        break
                if not homomorphism:
                    break
            key = mapping.tobytes()
            if key in unique:
                raise RuntimeError("duplicate Z_39 x Z_3 automorphism")
            unique.add(key)
            automorphisms.append(
                {
                    "unit_13": unit_13,
                    "matrix_mod_3": list(matrix),
                    "coordinate_map": mapping,
                    "bijective": bijective,
                    "homomorphism": homomorphism,
                }
            )
    if len(automorphisms) != 576:
        raise RuntimeError("full Z_39 x Z_3 automorphism count is not 576")
    if not all(item["bijective"] and item["homomorphism"] for item in automorphisms):
        raise RuntimeError("invalid map in Z_39 x Z_3 automorphism audit")
    return tuple(automorphisms)


def automorphism_audit() -> dict[str, Any]:
    automorphisms = z39_z3_automorphisms()
    return {
        "schema": "exp066-z39-z3-automorphism-audit-v1",
        "decomposition": "Z_39 x Z_3 ~= Z_13 x F_3^2",
        "expected_automorphisms": 12 * 48,
        "automorphisms": len(automorphisms),
        "unique_coordinate_permutations": len(
            {item["coordinate_map"].tobytes() for item in automorphisms}
        ),
        "all_bijective": all(item["bijective"] for item in automorphisms),
        "all_homomorphisms": all(
            item["homomorphism"] for item in automorphisms
        ),
        "coordinate_maps_sha256": canonical_sha256(
            [item["coordinate_map"].astype(int).tolist() for item in automorphisms]
        ),
        "argument": (
            "Sylow factors are characteristic: Z_39 x Z_3 is Z_13 x F_3^2. "
            "Every automorphism is therefore one nonzero F_13 scalar and one "
            "GL(2,3) matrix; the 12*48 explicit bijective homomorphisms exhaust it."
        ),
    }


def _support_pair(record: dict[str, Any], ids: dict[tuple, int]) -> tuple[int, int]:
    left = ids[E55._canon_translate(record["A"], 39, 3)]
    right = ids[E55._canon_translate(record["B"], 39, 3)]
    return (left, right) if left <= right else (right, left)


def _support_action_maps(
    ids: dict[tuple, int], representatives: list[tuple]
) -> tuple[tuple[int, ...], ...]:
    maps: list[tuple[int, ...]] = []
    for item in z39_z3_automorphisms():
        unit = int(item["unit_13"])
        matrix = tuple(map(int, item["matrix_mod_3"]))
        maps.append(
            tuple(
                ids[
                    E55._canon_translate(
                        [_map_coordinate(term, unit, matrix) for term in support],
                        39,
                        3,
                    )
                ]
                for support in representatives
            )
        )
    return tuple(maps)


def automorphism_bundles(
    residual_records: list[dict[str, Any]],
) -> dict[str, Any]:
    """Partition residual comparison classes under the full exact action."""
    ids, representatives, _ = E55._orbit_tables(39, 3, 3)
    action_maps = _support_action_maps(ids, representatives)
    pair_to_index: dict[tuple[int, int], int] = {}
    for index, record in enumerate(residual_records):
        pair = _support_pair(record, ids)
        if pair in pair_to_index:
            raise RuntimeError("duplicate residual support pair")
        pair_to_index[pair] = index
    residual_pairs = set(pair_to_index)

    orbit_cache: dict[tuple[int, int], frozenset[tuple[int, int]]] = {}
    for pair in residual_pairs:
        images = set()
        for action in action_maps:
            left, right = action[pair[0]], action[pair[1]]
            image = (left, right) if left <= right else (right, left)
            if image in residual_pairs:
                images.add(image)
        if pair not in images:
            raise RuntimeError("automorphism orbit omitted its seed")
        orbit_cache[pair] = frozenset(images)
    for pair, orbit in orbit_cache.items():
        if any(orbit_cache[image] != orbit for image in orbit):
            raise RuntimeError("residual automorphism relation is not a partition")

    unseen = set(residual_pairs)
    bundles = []
    while unseen:
        seed = min(unseen)
        orbit = orbit_cache[seed]
        indexes = sorted(pair_to_index[pair] for pair in orbit)
        bundles.append(
            {
                "bundle_index": len(bundles),
                "representative_record": indexes[0],
                "record_indexes": indexes,
                "size": len(indexes),
            }
        )
        unseen.difference_update(orbit)
    covered = sorted(index for bundle in bundles for index in bundle["record_indexes"])
    if covered != list(range(len(residual_records))):
        raise RuntimeError("automorphism bundles do not cover every residual record")
    return {
        "audit": automorphism_audit(),
        "residual_classes": len(residual_records),
        "bundles": len(bundles),
        "bundle_sizes": sorted(bundle["size"] for bundle in bundles),
        "partition_sha256": canonical_sha256(bundles),
        "records": bundles,
        "use": (
            "Exact code-equivalence bundles for witness transport and residual "
            "solver scheduling; every transported support is physically checked "
            "against the target matrices."
        ),
    }


def _problem_matrices(record: dict[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    return E55.E53.bb_from_terms(
        int(record["ell"]), int(record["m"]), record["A"], record["B"]
    )


def _translated_support(
    support: set[tuple[int, int]], delta: tuple[int, int]
) -> set[tuple[int, int]]:
    return {
        ((x + delta[0]) % 39, (y + delta[1]) % 3)
        for x, y in support
    }


def _translation_delta(
    source: set[tuple[int, int]], target: set[tuple[int, int]]
) -> tuple[int, int] | None:
    for dx in range(39):
        for dy in range(3):
            delta = (dx, dy)
            if _translated_support(source, delta) == target:
                return delta
    return None


def _permuted_rows(matrix: np.ndarray, mapping: np.ndarray) -> np.ndarray:
    output = np.zeros_like(matrix)
    output[:, mapping] = matrix
    return output


def _transport_mapping(
    source: dict[str, Any], target: dict[str, Any]
) -> tuple[np.ndarray, dict[str, Any]]:
    """Find and independently verify a source-to-target qubit permutation."""
    if (
        (int(source["ell"]), int(source["m"]))
        != (int(target["ell"]), int(target["m"]))
        or (int(source["ell"]), int(source["m"])) != (39, 3)
    ):
        raise RuntimeError("EXP-066 transport is specific to the (39,3) lattice")
    source_a = {(int(x), int(y)) for x, y in source["A"]}
    source_b = {(int(x), int(y)) for x, y in source["B"]}
    target_a = {(int(x), int(y)) for x, y in target["A"]}
    target_b = {(int(x), int(y)) for x, y in target["B"]}
    source_hx, source_hz = _problem_matrices(source)
    target_hx, target_hz = _problem_matrices(target)
    target_hx_key = E55.canon(target_hx)
    target_hz_key = E55.canon(target_hz)

    for automorphism_index, item in enumerate(z39_z3_automorphisms()):
        coordinate_map = item["coordinate_map"].astype(int)
        mapped_a = {
            divmod(int(coordinate_map[3 * x + y]), 3) for x, y in source_a
        }
        mapped_b = {
            divmod(int(coordinate_map[3 * x + y]), 3) for x, y in source_b
        }
        for swap_blocks in (False, True):
            delta_a = _translation_delta(
                mapped_b if swap_blocks else mapped_a, target_a
            )
            delta_b = _translation_delta(
                mapped_a if swap_blocks else mapped_b, target_b
            )
            if delta_a is None or delta_b is None:
                continue
            mapping = np.empty(234, dtype=int)
            for old_coordinate in range(117):
                mapped_x, mapped_y = divmod(
                    int(coordinate_map[old_coordinate]), 3
                )
                target_a_coordinate = (
                    3 * ((mapped_x + delta_a[0]) % 39)
                    + (mapped_y + delta_a[1]) % 3
                )
                target_b_coordinate = (
                    3 * ((mapped_x + delta_b[0]) % 39)
                    + (mapped_y + delta_b[1]) % 3
                )
                if swap_blocks:
                    mapping[old_coordinate] = 117 + target_b_coordinate
                    mapping[117 + old_coordinate] = target_a_coordinate
                else:
                    mapping[old_coordinate] = target_a_coordinate
                    mapping[117 + old_coordinate] = 117 + target_b_coordinate
            if sorted(mapping.tolist()) != list(range(234)):
                continue
            if (
                E55.canon(_permuted_rows(source_hx, mapping)) != target_hx_key
                or E55.canon(_permuted_rows(source_hz, mapping)) != target_hz_key
            ):
                continue
            return mapping, {
                "schema": "exp066-code-transport-v1",
                "automorphism_index": automorphism_index,
                "unit_13": int(item["unit_13"]),
                "matrix_mod_3": item["matrix_mod_3"],
                "swap_blocks": swap_blocks,
                "delta_A": list(delta_a),
                "delta_B": list(delta_b),
                "qubit_permutation_sha256": canonical_sha256(mapping.tolist()),
                "HX_rowspace_preserved": True,
                "HZ_rowspace_preserved": True,
            }
    raise RuntimeError("no matrix-verified automorphism transport exists")


def transport_witness(
    source: dict[str, Any],
    target: dict[str, Any],
    support: list[int],
) -> dict[str, Any]:
    mapping, metadata = _transport_mapping(source, target)
    vector = np.zeros(int(source["n"]), dtype=np.uint8)
    vector[np.asarray(support, dtype=int)] = 1
    transported = np.zeros_like(vector)
    transported[mapping] = vector
    target_support = [int(index) for index in np.flatnonzero(transported)]
    verification = verify_witness(target, target_support)
    if verification["weight"] != len(support):
        raise RuntimeError("automorphism transport did not preserve Hamming weight")
    return {
        "support": target_support,
        "verification": verification,
        "transport": metadata,
    }


def verify_witness(record: dict[str, Any], support: list[int]) -> dict[str, Any]:
    hx, hz = _problem_matrices(record)
    n = int(record["n"])
    if len(set(map(int, support))) != len(support):
        raise RuntimeError("witness support contains duplicates")
    if any(int(index) < 0 or int(index) >= n for index in support):
        raise RuntimeError("witness support lies outside the physical code")
    vector = np.zeros(n, dtype=np.uint8)
    vector[np.asarray(support, dtype=int)] = 1
    commutes = not np.any(hx.astype(np.int64) @ vector.astype(np.int64) % 2)
    outside = rank_np(np.vstack([hz, vector])) == rank_np(hz) + 1
    if not commutes or not outside:
        raise RuntimeError("candidate-local physical witness failed verification")
    return {
        "weight": int(vector.sum()),
        "commutes_with_HX": True,
        "outside_Z_stabilizer": True,
        "vector_sha256": E56.matrix_sha256(vector[None, :]),
    }


def _deep_seed(record: dict[str, Any]) -> int:
    digest = canonical_sha256(
        {
            "schema": SCHEMA,
            "ell": record["ell"],
            "m": record["m"],
            "A": record["A"],
            "B": record["B"],
            "threshold": record["threshold"],
        }
    )
    return int(digest[:16], 16)


def deepen_record(
    record: dict[str, Any], *, tries: int = DEEP_TRIES
) -> dict[str, Any]:
    ell, m = int(record["ell"]), int(record["m"])
    _, hz = _problem_matrices(record)
    pole = E55.intersect(
        E55.ann_basis(record["A"], ell, m),
        E55.ann_basis(record["B"], ell, m),
    )
    physical_pole = pole[:, E55.bar_permutation(ell, m)]
    seed = _deep_seed(record)
    started = time.perf_counter()
    result = E55.reduced_witness_bound(
        physical_pole,
        hz,
        ell,
        m,
        tries=int(tries),
        keep=E55.WITNESS_KEEP,
        seed=seed,
        target=int(record["threshold"]),
    )
    support = result.get("witness_support")
    verification = verify_witness(record, support) if support is not None else None
    if verification is not None and verification["weight"] != int(result["bound"]):
        raise RuntimeError("deep-reduction weight does not match its support")
    return {
        "schema": "exp066-deep-information-set-reduction-v1",
        "tries": int(tries),
        "seed": seed,
        "classes_tested": int(result.get("classes_tested", 0)),
        "bound": result.get("bound"),
        "witness_support": support,
        "verification": verification,
        "wall_time_s": time.perf_counter() - started,
    }


def _finalize_residual(
    record: dict[str, Any], deep: dict[str, Any]
) -> dict[str, Any]:
    output = dict(record)
    output["initial_verdict"] = record["verdict"]
    output["deep_reduction"] = deep
    bound = deep.get("bound")
    if bound is not None:
        output["witness_bound"] = int(bound)
        output["witness_support"] = deep["witness_support"]
    if bound is not None and int(bound) <= int(record["threshold"]):
        output.update(
            {
                "verdict": "dominated_by_deep_reduction",
                "solver_calls": 0,
                "wall_s": round(float(deep["wall_time_s"]), 3),
            }
        )
    else:
        output.update({"verdict": "undecided", "solver_calls": 0})
    return output


def _finalize_transported(
    source: dict[str, Any],
    target: dict[str, Any],
    deep: dict[str, Any],
    bundle_index: int,
    source_residual_index: int,
) -> dict[str, Any]:
    bound = deep.get("bound")
    support = deep.get("witness_support")
    if bound is None or support is None or int(bound) > int(target["threshold"]):
        raise RuntimeError("transport source does not dominate the target threshold")
    transported = transport_witness(source, target, support)
    output = dict(target)
    output.update(
        {
            "initial_verdict": target["verdict"],
            "verdict": "dominated_by_automorphism_transport",
            "solver_calls": 0,
            "witness_bound": transported["verification"]["weight"],
            "witness_support": transported["support"],
            "automorphism_transport": {
                "bundle_index": bundle_index,
                "source_residual_index": source_residual_index,
                "source": {
                    "A": source["A"],
                    "B": source["B"],
                    "witness_support": support,
                    "deep_reduction": deep,
                },
                **transported["transport"],
            },
            "wall_s": round(float(deep["wall_time_s"]), 3),
        }
    )
    return output


def _screen_protocol() -> dict[str, Any]:
    return E55._screen_protocol(
        E55._file_sha256(E55.OUT),
        E55._reference_fingerprint(),
        K_MIN,
        K_MAX,
        TIME_LIMIT_S,
    )


def _lattice_payload(
    ell: int,
    m: int,
    records: list[dict[str, Any]],
    protocol: dict[str, Any],
    wall_s: float,
) -> dict[str, Any]:
    verdicts = Counter(record["verdict"] for record in records)
    return {
        "schema": "exp055-screen-v3",
        "utc": E56.utc_now(),
        "protocol": protocol,
        "exp066": {
            "schema": SCHEMA,
            "deep_tries": DEEP_TRIES,
            "bundle_witness_transport": True,
        },
        "ell": ell,
        "m": m,
        "n": 2 * ell * m,
        "k_range": [K_MIN, K_MAX],
        "time_limit_s": TIME_LIMIT_S,
        "candidates_after_symmetry": len(records),
        "orbit_total": sum(int(record["orbit"]) for record in records),
        "verdicts": dict(sorted(verdicts.items())),
        "survivors": [record for record in records if record["verdict"] == "survivor"],
        "no_reference": [
            record for record in records if record["verdict"] == "no_reference"
        ],
        "undecided": [record for record in records if record["verdict"] == "undecided"],
        "solver_calls": sum(int(record.get("solver_calls", 0)) for record in records),
        "records": records,
        "wall_s": round(wall_s, 1),
    }


def screen_lattice(
    ell: int,
    m: int,
    *,
    deep_tries: int = DEEP_TRIES,
    deep_workers: int = DEEP_WORKERS,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    started = time.perf_counter()
    candidates = E55.enumerate_candidates(ell, m, K_MIN, K_MAX)
    initial = [
        E55._screen_one((ell, m, candidate, TIME_LIMIT_S, False))
        for candidate in candidates
    ]
    residual_indexes = [
        index
        for index, record in enumerate(initial)
        if record["verdict"] == "solver_required"
    ]
    residual = [initial[index] for index in residual_indexes]
    records = list(initial)
    if residual:
        if (ell, m) != (39, 3):
            raise RuntimeError("unexpected EXP-066 residual outside (39,3)")
        bundle_records = automorphism_bundles(residual)["records"]
        representative_indexes = [
            int(bundle["representative_record"]) for bundle in bundle_records
        ]
        with ThreadPoolExecutor(
            max_workers=min(deep_workers, len(representative_indexes))
        ) as pool:
            representative_deep = list(
                pool.map(
                    lambda index: deepen_record(residual[index], tries=deep_tries),
                    representative_indexes,
                )
            )
        deep_by_residual = dict(zip(representative_indexes, representative_deep))

        open_bundles = [
            bundle
            for bundle, deep in zip(bundle_records, representative_deep)
            if deep.get("bound") is None
            or int(deep["bound"]) > int(
                residual[int(bundle["representative_record"])]["threshold"]
            )
        ]
        alternate_indexes = [
            int(index)
            for bundle in open_bundles
            for index in bundle["record_indexes"]
            if int(index) != int(bundle["representative_record"])
        ]
        if alternate_indexes:
            with ThreadPoolExecutor(
                max_workers=min(deep_workers, len(alternate_indexes))
            ) as pool:
                alternate_deep = list(
                    pool.map(
                        lambda index: deepen_record(
                            residual[index], tries=deep_tries
                        ),
                        alternate_indexes,
                    )
                )
            deep_by_residual.update(zip(alternate_indexes, alternate_deep))

        for bundle in bundle_records:
            source_index = next(
                (
                    int(index)
                    for index in bundle["record_indexes"]
                    if deep_by_residual.get(int(index), {}).get("bound")
                    is not None
                    and int(deep_by_residual[int(index)]["bound"])
                    <= int(residual[int(index)]["threshold"])
                ),
                None,
            )
            for residual_index in map(int, bundle["record_indexes"]):
                record_index = residual_indexes[residual_index]
                if source_index is None:
                    deep = deep_by_residual.get(residual_index)
                    if deep is None:
                        deep = deepen_record(
                            residual[residual_index], tries=deep_tries
                        )
                    records[record_index] = _finalize_residual(
                        residual[residual_index], deep
                    )
                else:
                    records[record_index] = _finalize_transported(
                        residual[source_index],
                        residual[residual_index],
                        deep_by_residual[source_index],
                        int(bundle["bundle_index"]),
                        source_index,
                    )
    payload = _lattice_payload(
        ell, m, records, _screen_protocol(), time.perf_counter() - started
    )
    E56.atomic_write_json(E55.SCREEN_DIR / f"{ell}x{m}.json", payload)
    return payload, residual


def _global_lattices_through_n234() -> list[tuple[int, int]]:
    census = json.loads(E55.OUT.read_text(encoding="utf-8"))
    lattices = [
        (int(record["ell"]), int(record["m"]))
        for record in census["lattices"]
        if record["frontier_weight3"] and int(record["n"]) <= 234
    ]
    if len(lattices) != 22 or lattices[-2:] != list(LATTICES):
        raise RuntimeError("unexpected nonempty-frontier lattice set through n=234")
    return lattices


def assemble_global_screen() -> dict[str, Any]:
    lattices = _global_lattices_through_n234()
    args = argparse.Namespace(
        k_min=K_MIN,
        k_max=K_MAX,
        time_limit=TIME_LIMIT_S,
        lattices=",".join(f"{ell}x{m}" for ell, m in lattices),
        min_n=18,
    )
    E55.assemble_screen(args)
    return json.loads(E55.SCREEN_OUT.read_text(encoding="utf-8"))


def _records_identity(records: list[dict[str, Any]]) -> str:
    return canonical_sha256(
        [
            {
                "A": record["A"],
                "B": record["B"],
                "k_parent": record["k_parent"],
                "orbit": record["orbit"],
                "threshold": record["threshold"],
                "threshold_source": record["threshold_source"],
            }
            for record in records
        ]
    )


def _validate_lattice_aggregates(lattice: dict[str, Any]) -> None:
    E55._validate_screen_shard_aggregates(lattice)
    E55._validate_screen_shard_records(lattice)


def _validate_global_screen(
    global_screen: dict[str, Any],
    target_shards: dict[tuple[int, int], dict[str, Any]] | None = None,
) -> dict[str, Any]:
    protocol = _screen_protocol()
    expected_lattices = _global_lattices_through_n234()
    lattices = global_screen.get("lattices", [])
    actual_lattices = [
        (int(lattice["ell"]), int(lattice["m"])) for lattice in lattices
    ]
    if (
        global_screen.get("schema") != "exp055-odd-lattice-screen-v2"
        or global_screen.get("protocol") != protocol
        or actual_lattices != expected_lattices
    ):
        raise RuntimeError("global screen identity is stale")
    if target_shards is None:
        target_shards = {
            pair: json.loads(
                (E55.SCREEN_DIR / f"{pair[0]}x{pair[1]}.json").read_text(
                    encoding="utf-8"
                )
            )
            for pair in LATTICES
        }
    embedded_by_lattice = {
        (int(lattice["ell"]), int(lattice["m"])): lattice
        for lattice in lattices
    }
    if any(
        embedded_by_lattice.get(pair) != shard
        for pair, shard in target_shards.items()
    ):
        raise RuntimeError(
            "global screen does not embed the validated EXP-066 shard"
        )

    verdicts: Counter[str] = Counter()
    survivors: list[dict[str, Any]] = []
    no_reference: list[dict[str, Any]] = []
    undecided: list[dict[str, Any]] = []
    for lattice in lattices:
        if (
            lattice.get("schema") != "exp055-screen-v3"
            or lattice.get("protocol") != protocol
            or int(lattice.get("n", -1))
            != 2 * int(lattice["ell"]) * int(lattice["m"])
        ):
            raise RuntimeError("embedded screen shard identity is stale")
        _validate_lattice_aggregates(lattice)
        verdicts.update(lattice["verdicts"])
        survivors.extend(lattice["survivors"])
        no_reference.extend(lattice["no_reference"])
        undecided.extend(lattice["undecided"])

    with_reference = sum(
        value for key, value in verdicts.items() if key != "no_reference"
    )
    dominated = sum(
        value for key, value in verdicts.items() if key.startswith("dominated")
    )
    expected_verdict = {
        "complete": True,
        "candidates_after_symmetry": sum(
            int(lattice["candidates_after_symmetry"]) for lattice in lattices
        ),
        "orbits_represented": sum(
            int(lattice["orbit_total"]) for lattice in lattices
        ),
        "solver_calls": sum(int(lattice["solver_calls"]) for lattice in lattices),
        "verdicts": dict(verdicts),
        "with_reference": with_reference,
        "dominated": dominated,
        "survivors": len(survivors),
        "no_reference": len(no_reference),
        "undecided": len(undecided),
        "all_referenced_decided": (
            dominated + len(survivors) == with_reference
        ),
        "all_referenced_dominated": dominated == with_reference,
    }
    expected_scope = {
        "lattices_expected": len(expected_lattices),
        "lattices_completed": len(lattices),
        "missing": [],
        "weight_A": 3,
        "weight_B": 3,
        "k_range": [K_MIN, K_MAX],
        "n_max": 234,
        "reference_rule": (
            "max independently exact-certified d with n_ref <= n_candidate "
            "and k_ref >= k_candidate; estimates excluded"
        ),
    }
    if global_screen.get("scope") != expected_scope:
        raise RuntimeError("global screen scope is inconsistent")
    if global_screen.get("verdict") != expected_verdict:
        raise RuntimeError("global screen verdict is inconsistent")
    if (
        global_screen.get("survivors") != survivors
        or global_screen.get("no_reference") != no_reference
        or global_screen.get("undecided") != undecided
    ):
        raise RuntimeError("global screen verdict lists are inconsistent")
    return expected_verdict


def validate_summary(summary: dict[str, Any]) -> None:
    if summary.get("schema") != SCHEMA:
        raise RuntimeError("unexpected EXP-066 summary schema")
    persisted_verdict = summary.get("verdict", {})
    if (
        persisted_verdict.get("exact_promotions") != []
        or persisted_verdict.get("exact_candidates_requiring_promotion") != []
    ):
        raise RuntimeError("EXP-066 has no replayable exact promotion")
    if summary.get("protocol") != _screen_protocol():
        raise RuntimeError("EXP-066 protocol is stale")
    if summary.get("automorphism_bundles", {}).get("audit") != automorphism_audit():
        raise RuntimeError("EXP-066 automorphism audit is stale")
    target_shards: dict[tuple[int, int], dict[str, Any]] = {}

    all_records: list[dict[str, Any]] = []
    initial_residual_39x3: list[dict[str, Any]] = []
    expected_shards: dict[str, dict[str, Any]] = {}
    for ell, m in LATTICES:
        path = E55.SCREEN_DIR / f"{ell}x{m}.json"
        shard = json.loads(path.read_text(encoding="utf-8"))
        target_shards[(ell, m)] = shard
        expected_hash = summary["shards"][f"{ell}x{m}"]["sha256"]
        if E55._file_sha256(path) != expected_hash:
            raise RuntimeError("EXP-066 shard hash mismatch")
        if (
            shard.get("schema") != "exp055-screen-v3"
            or shard.get("protocol") != summary["protocol"]
            or (shard.get("ell"), shard.get("m")) != (ell, m)
        ):
            raise RuntimeError("EXP-066 shard identity mismatch")
        if shard.get("exp066") != {
            "schema": SCHEMA,
            "deep_tries": DEEP_TRIES,
            "bundle_witness_transport": True,
        }:
            raise RuntimeError("EXP-066 shard method metadata is stale")
        _validate_lattice_aggregates(shard)
        expected = E55.enumerate_candidates(ell, m, K_MIN, K_MAX)
        if len(expected) != len(shard["records"]):
            raise RuntimeError("EXP-066 shard omitted comparison classes")
        if _records_identity(shard["records"]) != summary["shards"][f"{ell}x{m}"]["records_sha256"]:
            raise RuntimeError("EXP-066 record identity mismatch")
        expected_shards[f"{ell}x{m}"] = {
            "path": str(path.relative_to(ROOT)),
            "sha256": E55._file_sha256(path),
            "records_sha256": _records_identity(shard["records"]),
        }
        for expected_record, record in zip(expected, shard["records"]):
            expected_identity = {
                "A": [list(term) for term in expected_record["A"]],
                "B": [list(term) for term in expected_record["B"]],
                "k_parent": expected_record["k_parent"],
                "orbit": expected_record["orbit"],
            }
            if any(
                record[key] != expected_identity[key]
                for key in expected_identity
            ):
                raise RuntimeError("EXP-066 record order does not match enumeration")
            threshold, source = E55.domination_threshold(2 * ell * m, record["k_parent"])
            if (record["threshold"], record["threshold_source"]) != (threshold, source):
                raise RuntimeError("EXP-066 record threshold is stale")
            verdict = record["verdict"]
            if verdict.startswith("dominated"):
                if threshold <= 0:
                    raise RuntimeError("domination verdict has no exact reference")
                support = record.get("witness_support") or record.get(
                    "ceiling_witness_support"
                )
                if support is None:
                    raise RuntimeError("domination record lacks a physical witness")
                verification = verify_witness(record, support)
                if verification["weight"] > int(record["threshold"]):
                    raise RuntimeError("domination witness exceeds exact threshold")
                transport = record.get("automorphism_transport")
                if transport is not None:
                    source_record = {
                        "ell": ell,
                        "m": m,
                        "n": 2 * ell * m,
                        "A": transport["source"]["A"],
                        "B": transport["source"]["B"],
                    }
                    source_support = transport["source"]["witness_support"]
                    source_verification = verify_witness(
                        source_record, source_support
                    )
                    if source_verification["weight"] > int(record["threshold"]):
                        raise RuntimeError(
                            "transport source exceeds the exact threshold"
                        )
                    rebuilt_transport = transport_witness(
                        source_record, record, source_support
                    )
                    persisted_mapping = {
                        key: value
                        for key, value in transport.items()
                        if key
                        not in {
                            "bundle_index",
                            "source_residual_index",
                            "source",
                        }
                    }
                    if (
                        rebuilt_transport["support"] != support
                        or rebuilt_transport["transport"] != persisted_mapping
                    ):
                        raise RuntimeError(
                            "automorphism witness transport is stale"
                        )
            elif verdict == "no_reference":
                if threshold != 0:
                    raise RuntimeError(
                        "no_reference verdict has an admissible exact reference"
                    )
            elif verdict == "undecided":
                if threshold <= 0:
                    raise RuntimeError("undecided verdict has no comparison target")
            elif verdict == "survivor":
                raise RuntimeError(
                    "EXP-066 survivor lacks a replayable lower-bound certificate"
                )
            else:
                raise RuntimeError("unknown EXP-066 verdict")
            if ell == 39 and record.get("initial_verdict") == "solver_required":
                initial_residual_39x3.append(record)
        all_records.extend(shard["records"])

    rebuilt_bundles = automorphism_bundles(initial_residual_39x3)
    if summary["automorphism_bundles"] != rebuilt_bundles:
        raise RuntimeError("EXP-066 automorphism bundles are stale")
    if summary.get("shards") != expected_shards:
        raise RuntimeError("EXP-066 shard summary is inconsistent")
    verdicts = Counter(record["verdict"] for record in all_records)
    survivors = verdicts.get("survivor", 0)
    undecided = verdicts.get("undecided", 0)
    dominated = sum(
        value for key, value in verdicts.items() if key.startswith("dominated")
    )
    local_with_reference = (
        len(all_records) - verdicts.get("no_reference", 0)
    )
    initial_counts = Counter(
        record.get("initial_verdict", record["verdict"])
        for record in all_records
    )
    expected_scope = {
        "lattices": [list(pair) for pair in LATTICES],
        "n": 234,
        "k_range": [K_MIN, K_MAX],
        "classes": len(all_records),
        "represented_pairs": sum(
            int(record["orbit"]) for record in all_records
        ),
    }
    expected_initial = {
        "dominated": sum(
            value
            for key, value in initial_counts.items()
            if key.startswith("dominated")
        ),
        "solver_required": initial_counts.get("solver_required", 0),
    }
    expected_local_verdict = {
        "complete": dominated + survivors == local_with_reference,
        "dominated": dominated,
        "survivors": survivors,
        "undecided": undecided,
        "all_referenced_dominated": dominated == local_with_reference,
        "exact_promotions": [],
        "exact_candidates_requiring_promotion": [],
    }
    if summary.get("scope") != expected_scope:
        raise RuntimeError("EXP-066 scope summary is inconsistent")
    if summary.get("initial") != expected_initial:
        raise RuntimeError("EXP-066 initial summary is inconsistent")
    if summary.get("routes") != dict(sorted(verdicts.items())):
        raise RuntimeError("EXP-066 route summary is inconsistent")
    if summary.get("verdict") != expected_local_verdict:
        raise RuntimeError("EXP-066 verdict summary is inconsistent")
    global_screen = json.loads(E55.SCREEN_OUT.read_text(encoding="utf-8"))
    if E55._file_sha256(E55.SCREEN_OUT) != summary.get("global_screen_sha256"):
        raise RuntimeError("EXP-066 global screen hash mismatch")
    expected_global_verdict = _validate_global_screen(
        global_screen, target_shards
    )
    if summary.get("global_verdict") != expected_global_verdict:
        raise RuntimeError("EXP-066 global verdict summary is inconsistent")


def run(*, force: bool = False) -> int:
    if SUMMARY.exists() and not force:
        summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
        validate_summary(summary)
        print(json.dumps(summary["verdict"], indent=1))
        return 0 if summary["verdict"]["all_referenced_dominated"] else 1

    protocol = _screen_protocol()
    shards: dict[str, Any] = {}
    initial_residuals: dict[str, list[dict[str, Any]]] = {}
    for ell, m in LATTICES:
        payload, residual = screen_lattice(ell, m)
        key = f"{ell}x{m}"
        shards[key] = payload
        initial_residuals[key] = residual
        print(
            json.dumps(
                {
                    "lattice": key,
                    "classes": payload["candidates_after_symmetry"],
                    "represented_pairs": payload["orbit_total"],
                    "verdicts": payload["verdicts"],
                }
            ),
            flush=True,
        )

    bundles = automorphism_bundles(initial_residuals["39x3"])
    global_screen = assemble_global_screen()
    records = [record for payload in shards.values() for record in payload["records"]]
    initial_counts = Counter(
        record.get("initial_verdict", record["verdict"]) for record in records
    )
    initial = {
        "dominated": sum(
            value for key, value in initial_counts.items() if key.startswith("dominated")
        ),
        "solver_required": initial_counts.get("solver_required", 0),
    }
    verdicts = Counter(record["verdict"] for record in records)
    dominated = sum(value for key, value in verdicts.items() if key.startswith("dominated"))
    survivors = verdicts.get("survivor", 0)
    undecided = verdicts.get("undecided", 0)
    local_with_reference = len(records) - verdicts.get("no_reference", 0)
    exact_candidates = [
        {
            "ell": record["ell"],
            "m": record["m"],
            "A": record["A"],
            "B": record["B"],
            "k": record["k_parent"],
            "upper_bound": record.get("witness_bound"),
        }
        for record in records
        if record["verdict"] == "survivor"
    ]
    summary = {
        "schema": SCHEMA,
        "utc": E56.utc_now(),
        "protocol": protocol,
        "scope": {
            "lattices": [list(pair) for pair in LATTICES],
            "n": 234,
            "k_range": [K_MIN, K_MAX],
            "classes": len(records),
            "represented_pairs": sum(int(record["orbit"]) for record in records),
        },
        "initial": initial,
        "automorphism_bundles": bundles,
        "routes": dict(sorted(verdicts.items())),
        "shards": {
            key: {
                "path": str((E55.SCREEN_DIR / f"{key}.json").relative_to(ROOT)),
                "sha256": E55._file_sha256(E55.SCREEN_DIR / f"{key}.json"),
                "records_sha256": _records_identity(payload["records"]),
            }
            for key, payload in shards.items()
        },
        "verdict": {
            "complete": dominated + survivors == local_with_reference,
            "dominated": dominated,
            "survivors": survivors,
            "undecided": undecided,
            "all_referenced_dominated": dominated == local_with_reference,
            "exact_promotions": [],
            "exact_candidates_requiring_promotion": exact_candidates,
        },
        "global_screen_sha256": E55._file_sha256(E55.SCREEN_OUT),
        "global_verdict": global_screen["verdict"],
    }
    E56.atomic_write_json(SUMMARY, summary)
    validate_summary(summary)
    print(json.dumps(summary["verdict"], indent=1))
    return 0 if summary["verdict"]["all_referenced_dominated"] else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--force", action="store_true")
    subparsers.add_parser("validate")
    subparsers.add_parser("automorphisms")
    args = parser.parse_args()
    if args.command == "run":
        return run(force=args.force)
    if args.command == "validate":
        summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
        validate_summary(summary)
        print(json.dumps(summary["verdict"], indent=1))
        return 0
    print(json.dumps(automorphism_audit(), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
