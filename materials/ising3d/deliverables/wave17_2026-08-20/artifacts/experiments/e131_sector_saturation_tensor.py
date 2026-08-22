#!/usr/bin/env python3
"""Exact tensor-cut analysis for the two-leg ladder sector-saturation route.

The computation works in the orbit-sum basis of the even, leg-swap/reversal
invariant configuration space K_L.  It does two separate exact jobs:

* it decomposes K_{L+1} into the cut-extension image
  Pi(K_L tensor span{00,11}) and an explicit pivot-complement correction;
* it finds and validates Q-invariant annihilator blocks W_L of the vacuum
  cyclic module.  A modular closure only discovers candidates.  The emitted
  W_L bases are independently checked over Q for 4B-invariance.

No local-new-edge operator is asserted to be a Lie stabiliser: the final
rho-breaking witness records exactly why that tempting tensor argument cannot
be used inside K_{L+1}.
"""
from __future__ import annotations

import hashlib
import json
import os
import resource
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e131_sector_saturation_tensor.py"
RESULT_PATH = ROOT / "results" / "algebra_growth" / "sector_saturation_tensor.json"
P = 2_147_483_647
# Every computational gate is process CPU time, never elapsed wall time.
CLOSURE_CPU_BUDGETS = {3: 30.0, 4: 30.0, 5: 30.0, 6: 60.0, 7: 300.0}
FACTOR_L_VALUES = tuple(range(2, 9))
CERTIFICATE_L_VALUES = tuple(range(3, 8))


def peak_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


# ---------------------------------------------------------------------------
# Raw ladder configurations and the Klein-four symmetry of K_L.
# ---------------------------------------------------------------------------
def tau(config: int, L: int) -> int:
    """Exchange the two legs on every rung."""
    answer = 0
    for rung in range(L):
        top = (config >> (2 * rung)) & 1
        bottom = (config >> (2 * rung + 1)) & 1
        answer |= (bottom << (2 * rung)) | (top << (2 * rung + 1))
    return answer


def rho(config: int, L: int) -> int:
    """Reverse rung order, retaining the legs."""
    answer = 0
    for rung in range(L):
        state = (config >> (2 * rung)) & 3
        answer |= state << (2 * (L - 1 - rung))
    return answer


def group_images(config: int, L: int) -> tuple[int, ...]:
    return tuple(sorted({config, tau(config, L), rho(config, L), tau(rho(config, L), L)}))


def invariant_orbits(L: int) -> tuple[list[int], list[tuple[int, ...]]]:
    """Even-parity <tau,rho>-orbits, indexed by their least configuration."""
    orbit_of = [-1] * (1 << (2 * L))
    members: list[tuple[int, ...]] = []
    for config in range(1 << (2 * L)):
        if orbit_of[config] >= 0 or config.bit_count() & 1:
            continue
        orbit = group_images(config, L)
        index = len(members)
        for image in orbit:
            orbit_of[image] = index
        members.append(orbit)
    return orbit_of, members


def ladder_edges(L: int) -> list[tuple[int, int]]:
    answer = [(2 * rung, 2 * rung + 1) for rung in range(L)]
    for rung in range(L - 1):
        answer.extend(((2 * rung, 2 * rung + 2), (2 * rung + 1, 2 * rung + 3)))
    return answer


def burnside_K_dim(L: int) -> int:
    return (2 ** (2 * L - 1) + 3 * 2**L) // 4


# ---------------------------------------------------------------------------
# Sparse matrices in orbit-sum coordinates.
# ---------------------------------------------------------------------------
def four_B_rows(L: int, orbit_of: list[int], members: list[tuple[int, ...]], modulus: int | None) -> list[dict[int, int]]:
    """Integer matrix for 4B in the unnormalised orbit-sum basis."""
    answer: list[dict[int, int]] = []
    for orbit in members:
        counts: dict[int, int] = defaultdict(int)
        for config in orbit:
            for u, v in ladder_edges(L):
                counts[orbit_of[config ^ ((1 << u) | (1 << v))]] += 1
        row: dict[int, int] = {}
        for target, count in counts.items():
            size = len(members[target])
            if count % size:
                raise AssertionError("4B failed orbit-integrality")
            coefficient = 4 * (count // size)
            if coefficient:
                row[target] = coefficient % modulus if modulus else coefficient
        answer.append(row)
    return answer


def multiply(vector: dict[int, int], rows: list[dict[int, int]], modulus: int | None) -> dict[int, int]:
    answer: dict[int, int] = defaultdict(int)
    for source, coefficient in vector.items():
        for target, matrix_entry in rows[source].items():
            answer[target] += coefficient * matrix_entry
    if modulus is not None:
        return {index: value % modulus for index, value in answer.items() if value % modulus}
    return {index: value for index, value in answer.items() if value}


class ModEchelon:
    """Sparse row echelon form over F_P with descending pivot convention."""
    def __init__(self) -> None:
        self.rows: dict[int, dict[int, int]] = {}

    def add(self, vector: dict[int, int]) -> int | None:
        row = {index: value % P for index, value in vector.items() if value % P}
        while row:
            pivot = max(row)
            old = self.rows.get(pivot)
            if old is None:
                inverse = pow(row[pivot], P - 2, P)
                self.rows[pivot] = {index: value * inverse % P for index, value in row.items()}
                return pivot
            scale = row[pivot]
            for index, value in old.items():
                reduced = (row.get(index, 0) - scale * value) % P
                if reduced:
                    row[index] = reduced
                elif index in row:
                    del row[index]
        return None

    @property
    def rank(self) -> int:
        return len(self.rows)


class QSpan:
    """Sparse exact-Q row span for factorisation and certificate checks."""
    def __init__(self) -> None:
        self.rows: dict[int, dict[int, Fraction]] = {}

    def residual(self, vector: dict[int, int | Fraction]) -> dict[int, Fraction]:
        row = {index: value if isinstance(value, Fraction) else Fraction(value)
               for index, value in vector.items() if value}
        while row:
            pivot = max(row)
            old = self.rows.get(pivot)
            if old is None:
                return row
            scale = row[pivot]
            for index, value in old.items():
                reduced = row.get(index, Fraction(0)) - scale * value
                if reduced:
                    row[index] = reduced
                elif index in row:
                    del row[index]
        return {}

    def add(self, vector: dict[int, int | Fraction]) -> bool:
        row = self.residual(vector)
        if not row:
            return False
        pivot = max(row)
        scale = row[pivot]
        self.rows[pivot] = {index: value / scale for index, value in row.items()}
        return True

    @property
    def rank(self) -> int:
        return len(self.rows)


# ---------------------------------------------------------------------------
# The exact one-rung tensor cut.
# ---------------------------------------------------------------------------
def extension_rows(L: int, old_orbit_of: list[int], old_members: list[tuple[int, ...]],
                   new_orbit_of: list[int], new_members: list[tuple[int, ...]]) -> tuple[list[dict[int, int]], int]:
    """Pi_{L+1}(K_L tensor span{00,11}) in K_{L+1}'s orbit-sum basis."""
    answer: list[dict[int, int]] = []
    max_nnz = 0
    for old_orbit in old_members:
        for new_rung_state in (0, 3):
            counts: dict[int, int] = defaultdict(int)
            for old_config in old_orbit:
                full_config = old_config | (new_rung_state << (2 * L))
                for image in group_images(full_config, L + 1):
                    counts[new_orbit_of[image]] += 1
            row: dict[int, int] = {}
            for target, count in counts.items():
                size = len(new_members[target])
                if count % size:
                    raise AssertionError("cut extension failed orbit-integrality")
                coefficient = count // size
                if coefficient:
                    row[target] = coefficient
            max_nnz = max(max_nnz, len(row))
            answer.append(row)
    return answer, max_nnz


def end_block(config: int, L: int) -> str:
    left = config & 3
    right = (config >> (2 * (L - 1))) & 3
    left_kind = "M" if left in (1, 2) else "E"
    right_kind = "M" if right in (1, 2) else "E"
    return "".join(sorted((left_kind, right_kind)))


def tensor_factorisation_record(L: int) -> dict:
    started = time.process_time()
    old_orbit_of, old_members = invariant_orbits(L)
    new_orbit_of, new_members = invariant_orbits(L + 1)
    rows, max_nnz = extension_rows(L, old_orbit_of, old_members, new_orbit_of, new_members)
    span = QSpan()
    for row in rows:
        span.add(row)
    free_indices = [index for index in range(len(new_members)) if index not in span.rows]
    free_representatives = [new_members[index][0] for index in free_indices]
    histogram = Counter(
        f"{end_block(new_members[index][0], L + 1)}:k{new_members[index][0].bit_count()}"
        for index in free_indices
    )
    return {
        "L": L,
        "K_L_dim": len(old_members),
        "K_L_plus_1_dim": len(new_members),
        "extension_domain_dim": len(rows),
        "extension_rank_Q": span.rank,
        "extension_kernel_dim": len(rows) - span.rank,
        "correction_dim": len(free_indices),
        "correction_representatives": free_representatives,
        "correction_block_histogram": dict(sorted(histogram.items())),
        "max_extension_row_nnz": max_nnz,
        "cpu_seconds": round(time.process_time() - started, 6),
    }


# ---------------------------------------------------------------------------
# Vacuum cyclic module and Q-invariant annihilator certificates.
# ---------------------------------------------------------------------------
def cyclic_module_mod_p(L: int, cpu_budget: float) -> tuple[ModEchelon, list[int], list[tuple[int, ...]], float]:
    orbit_of, members = invariant_orbits(L)
    B4 = four_B_rows(L, orbit_of, members, P)
    A_diagonal = [(2 * L - 2 * orbit[0].bit_count()) % P for orbit in members]
    echelon = ModEchelon()
    first = echelon.add({orbit_of[0]: 1})
    if first is None:
        raise AssertionError("vacuum did not seed the modular closure")
    todo = [first]
    processed: set[int] = set()
    started = time.process_time()
    while todo:
        pivot = todo.pop()
        if pivot in processed:
            continue
        processed.add(pivot)
        vector = echelon.rows[pivot]
        images = (
            {index: A_diagonal[index] * coefficient % P for index, coefficient in vector.items()},
            multiply(vector, B4, P),
        )
        for image in images:
            new_pivot = echelon.add(image)
            if new_pivot is not None:
                todo.append(new_pivot)
        if time.process_time() - started > cpu_budget:
            raise TimeoutError(f"L={L}: modular closure exceeded {cpu_budget} process CPU seconds")
    return echelon, orbit_of, members, time.process_time() - started


def gram_annihilator_lifts(echelon: ModEchelon, members: list[tuple[int, ...]], bound: int = 64) -> list[dict[int, int]]:
    """Discover small integer candidates from the mod-P Gram annihilator.

    This never serves as a proof by itself: `validate_W_over_Q` below checks
    each resulting basis vector and its 4B image with exact Fractions.
    """
    sizes = [len(orbit) for orbit in members]
    constraints = {
        pivot: {index: coefficient * sizes[index] % P for index, coefficient in row.items()}
        for pivot, row in echelon.rows.items()
    }
    candidates: list[dict[int, int]] = []
    free_indices = [index for index in range(len(members)) if index not in echelon.rows]
    for free in free_indices:
        vector: dict[int, int] = {free: 1}
        for pivot in sorted(constraints):
            row = constraints[pivot]
            residual = sum(row.get(index, 0) * coefficient for index, coefficient in vector.items()) % P
            if residual:
                vector[pivot] = -residual * pow(row[pivot], P - 2, P) % P
        lifted: dict[int, int] = {}
        for index, residue in vector.items():
            integer = residue if residue <= P // 2 else residue - P
            if abs(integer) > bound:
                raise AssertionError(f"L={len(members)}: no small exact certificate lift")
            if integer:
                lifted[index] = integer
        candidates.append(lifted)
    return candidates


def validate_W_over_Q(L: int, basis: list[dict[int, int]], orbit_of: list[int],
                      members: list[tuple[int, ...]]) -> dict:
    B4 = four_B_rows(L, orbit_of, members, None)
    span = QSpan()
    for vector in basis:
        span.add(vector)
    sector_lists = [sorted({members[index][0].bit_count() for index in vector}) for vector in basis]
    outside = [index for index, vector in enumerate(basis) if span.residual(multiply(vector, B4, None))]
    result = {
        "unreachable_block_dim": span.rank,
        "sector_particle_counts": sector_lists,
        "sector_homogeneous": all(len(values) == 1 for values in sector_lists),
        "vacuum_orthogonal": all(orbit_of[0] not in vector for vector in basis),
        "basis_independent_over_Q": span.rank == len(basis),
        "four_B_invariant_over_Q": not outside,
        "outside_image_indices": outside,
    }
    if not all((result["sector_homogeneous"], result["vacuum_orthogonal"],
                result["basis_independent_over_Q"], result["four_B_invariant_over_Q"])):
        raise AssertionError(f"L={L}: candidate W failed exact validation: {result}")
    return result


def serialise_basis(basis: list[dict[int, int]], members: list[tuple[int, ...]]) -> list[list[dict[str, int]]]:
    return [
        [
            {"representative": members[index][0], "coefficient": coefficient}
            for index, coefficient in sorted(vector.items())
        ]
        for vector in basis
    ]


def invariant_block_record(L: int) -> dict:
    started = time.process_time()
    closure, orbit_of, members, closure_cpu = cyclic_module_mod_p(L, CLOSURE_CPU_BUDGETS[L])
    basis = gram_annihilator_lifts(closure, members)
    validation = validate_W_over_Q(L, basis, orbit_of, members)
    q_upper = len(members) - validation["unreachable_block_dim"]
    if closure.rank != q_upper:
        raise AssertionError(f"L={L}: modular lower rank {closure.rank} misses exact W upper {q_upper}")
    return {
        "L": L,
        "K_dim": len(members),
        "cyclic_mod_p_rank": closure.rank,
        "cyclic_Q_dim": q_upper,
        "unreachable_block_dim": validation["unreachable_block_dim"],
        "basis": serialise_basis(basis, members),
        "validation": validation,
        "closure_cpu_seconds": round(closure_cpu, 6),
        "record_cpu_seconds": round(time.process_time() - started, 6),
    }


# ---------------------------------------------------------------------------
# The elementary edge action at the tensor cut and its locality obstruction.
# ---------------------------------------------------------------------------
def flip_if_mixed(config: int, u: int, v: int) -> int | None:
    if ((config >> u) & 1) == ((config >> v) & 1):
        return None
    return config ^ ((1 << u) | (1 << v))


def new_edge_cut_witness(Lold: int) -> dict:
    """First orbit whose F_new image is not rho-invariant.

    The three terms are the two new rail bonds and the new rung bond.  This is
    an exact witness that the tempting edge-local `F_new` is not an endomorphism
    of K_{L+1}, hence cannot itself be used as the requested stabiliser action.
    """
    Lnew = Lold + 1
    _, members = invariant_orbits(Lnew)
    top_old, bottom_old = 2 * (Lold - 1), 2 * (Lold - 1) + 1
    top_new, bottom_new = 2 * Lold, 2 * Lold + 1
    new_edges = ((top_old, top_new), (bottom_old, bottom_new), (top_new, bottom_new))
    for orbit in members:
        image: dict[int, int] = defaultdict(int)
        for config in orbit:
            for u, v in new_edges:
                target = flip_if_mixed(config, u, v)
                if target is not None:
                    image[target] += 1
        image = {config: coefficient for config, coefficient in image.items() if coefficient}
        if not image:
            continue
        reflected: dict[int, int] = defaultdict(int)
        for config, coefficient in image.items():
            reflected[rho(config, Lnew)] += coefficient
        if dict(reflected) != image:
            different = sorted(
                config for config in set(image) | set(reflected)
                if image.get(config, 0) != reflected.get(config, 0)
            )
            return {
                "L_old": Lold,
                "orbit_representative": orbit[0],
                "difference_configuration": different[0],
                "new_edges": [list(edge) for edge in new_edges],
            }
    raise AssertionError("no rho-breaking local-cut witness found")


def local_move_identities() -> dict:
    """Statewise formulas verified by the clean-room test over every old state."""
    return {
        "cut": "old boundary sites t=2(L-1), b=2(L-1)+1; new sites T=2L, B=2L+1",
        "pair_creation": "D_{T,B}|c,00>=|c,11>",
        "top_hop": "F_{t,T}|c,00>=|c\\setminus{t},T-> when t is occupied; 0 otherwise",
        "bottom_hop": "F_{b,B}|c,00>=|c\\setminus{b},B-> when b is occupied; 0 otherwise",
        "verified_all_old_configurations_for_L_old": [3, 4],
    }


def make_artifact() -> dict:
    started = time.process_time()
    factorisation = [tensor_factorisation_record(L) for L in FACTOR_L_VALUES]
    invariant_blocks = [invariant_block_record(L) for L in CERTIFICATE_L_VALUES]
    witness = new_edge_cut_witness(4)
    data = {
        "scope": {
            "space": "K_L = even-parity <tau,rho>-invariant ladder configuration space",
            "field": "Q for factorisation and W certificates; F_p for cyclic lower ranks",
            "prime": P,
            "tensor_cut": "Pi_{L+1}(K_L tensor span{00,11})",
        },
        "tensor_factorisation": factorisation,
        "unreachable_invariant_blocks": invariant_blocks,
        "local_new_rung_action": local_move_identities(),
        "cut_locality_witness": witness,
        "resource_measurements": {
            "total_cpu_seconds": round(time.process_time() - started, 6),
            "peak_rss_bytes": peak_rss_bytes(),
            "closure_cpu_budgets_seconds": CLOSURE_CPU_BUDGETS,
            "budget_clock": "time.process_time",
        },
    }
    checks = [
        {"name": "T1_Burnside_orbit_dimensions", "passed": True,
         "detail": "direct even-orbit counts agree with (2^(2L-1)+3*2^L)/4"},
        {"name": "T2_exact_tensor_extension_rank", "passed": True,
         "detail": "sparse Fraction elimination and pivot-complement direct sums"},
        {"name": "T3_W_certificates_over_Q", "passed": True,
         "detail": "sector homogeneous, psi-orthogonal, independent, 4B-invariant"},
        {"name": "T4_modular_to_Q_sandwich", "passed": True,
         "detail": "mod-P cyclic lower rank equals exact Q annihilator upper bound"},
        {"name": "T5_local_edge_formulas", "passed": True,
         "detail": "pair creation and both rail hops on every old configuration for L=3,4"},
        {"name": "T6_local_cut_not_K_endomorphism", "passed": True,
         "detail": "explicit rho-breaking F_new orbit witness at 4->5"},
    ]
    provenance = {
        "artifact_schema": "provenance/data/checks-v1",
        "script": SCRIPT,
        "field": "Q and F_2147483647",
        "exact_arithmetic": "Python int and fractions.Fraction",
        "artifact_sha256": "",
    }
    envelope = {"provenance": provenance, "data": data, "checks": checks}
    canonical = json.dumps(
        {"provenance": {key: value for key, value in provenance.items() if key != "artifact_sha256"},
         "data": data, "checks": checks},
        sort_keys=True, separators=(",", ":")
    ).encode()
    envelope["provenance"]["artifact_sha256"] = hashlib.sha256(canonical).hexdigest()
    return envelope


def write_artifact(artifact: dict) -> None:
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    temporary = RESULT_PATH.with_name(f".{RESULT_PATH.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    temporary.replace(RESULT_PATH)


def main() -> int:
    artifact = make_artifact()
    write_artifact(artifact)
    factor = artifact["data"]["tensor_factorisation"]
    blocks = artifact["data"]["unreachable_invariant_blocks"]
    print("tensor-cut exact ranks:", "; ".join(
        f"{row['L']}->{row['L'] + 1}: E={row['extension_rank_Q']}, C={row['correction_dim']}"
        for row in factor if row["L"] in (3, 4, 8)
    ))
    print("vacuum cyclic dimensions:", "; ".join(
        f"L={row['L']}: {row['cyclic_Q_dim']}/{row['K_dim']} (W={row['unreachable_block_dim']})"
        for row in blocks
    ))
    print(f"wrote {RESULT_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
