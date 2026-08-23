#!/usr/bin/env python3
"""Exact structural model behind the conjectural ladder symmetric-square law.

This module does not compute a ladder closure.  It constructs the parameter-free
candidate

    Sym^2(exterior^m Q^L),

identifies the unique possible SL_L-invariant symmetric top-wedge trace, and
computes the kernel forced by rung-reflection averaging.  The latter is a
route-closure certificate: the obvious balanced-leg/Reynolds map cannot be the
conjectural isomorphism.

All arithmetic is Python integer arithmetic.  Explicit subset and unordered-
pair enumeration is used as a check of every closed formula for 2 <= L <= 9.
"""
from __future__ import annotations

import json
from itertools import combinations
from math import comb

SMALL_L = tuple(range(2, 10))


def subset_masks(L: int, m: int) -> list[int]:
    """Return the m-subsets of [0,L), in a deterministic numeric order."""
    return sorted(sum(1 << r for r in choice) for choice in combinations(range(L), m))


def reverse_mask(mask: int, L: int) -> int:
    out = 0
    for r in range(L):
        if (mask >> r) & 1:
            out |= 1 << (L - 1 - r)
    return out


def fixed_subset_count(L: int, m: int) -> int:
    """Number of m-subsets fixed by reversal, from its 1/2-cycle structure."""
    half = L // 2
    if L % 2 == 0:
        return comb(half, m // 2) if m % 2 == 0 else 0
    return comb(half, m // 2)


def candidate_defect(L: int, m: int) -> int:
    return int(L % 4 == 0 and 2 * m == L)


def wedge_sign(left: int, right: int, L: int) -> int:
    """Coefficient of the ordered volume form in e_left wedge e_right."""
    full = (1 << L) - 1
    if left & right or (left | right) != full:
        return 0
    inversions = 0
    for a in range(L):
        if (left >> a) & 1:
            inversions += (right & ((1 << a) - 1)).bit_count()
    return -1 if inversions % 2 else 1


def top_wedge_record(L: int, m: int, subsets: list[int]) -> dict:
    """Exact finite record for the canonical top-wedge bilinear functional."""
    top_degree = 2 * m == L
    admissible = 0
    nonzero_symmetric = 0
    if top_degree:
        full = (1 << L) - 1
        index = set(subsets)
        for left in subsets:
            right = full ^ left
            if right not in index or left > right:
                continue
            admissible += 1
            coefficient = wedge_sign(left, right, L) + wedge_sign(right, left, L)
            if coefficient:
                nonzero_symmetric += 1
    exchange_sign = (-1 if m % 2 else 1) if top_degree else None
    symmetric_functional_dim = int(top_degree and m % 2 == 0)
    return {
        "top_degree": top_degree,
        "exchange_sign": exchange_sign,
        "torus_weight_admissible_unordered_pairs": admissible,
        "nonzero_symmetrized_wedge_coefficients": nonzero_symmetric,
        "symmetric_invariant_functional_dim": symmetric_functional_dim,
        "defect_matches_functional_dim": candidate_defect(L, m) == symmetric_functional_dim,
        "tag": "[LEMMA] exact exterior top-wedge calculation; uniqueness is proved in proofs/wlaw.md",
    }


def reflection_pair_counts(L: int, m: int, subsets: list[int]) -> dict:
    """Enumerate reversal orbits on unordered pairs of m-subsets."""
    pairs = [(subsets[i], subsets[j]) for i in range(len(subsets)) for j in range(i, len(subsets))]
    seen: set[tuple[int, int]] = set()
    orbit_count = 0
    fixed_pair_count = 0
    for pair in pairs:
        if pair in seen:
            continue
        reflected = tuple(sorted((reverse_mask(pair[0], L), reverse_mask(pair[1], L))))
        seen.add(pair)
        seen.add(reflected)
        orbit_count += 1
        if reflected == pair:
            fixed_pair_count += 1
    return {
        "unordered_pair_count": len(pairs),
        "fixed_unordered_pair_count": fixed_pair_count,
        "reversal_orbit_count": orbit_count,
    }


def balanced_orbit_count(L: int, m: int, subsets: list[int]) -> int:
    """Count <tau,rho>-orbits of balanced (m,m) leg configurations directly.

    A configuration is encoded as (top_mask, bottom_mask); tau swaps the legs
    and rho reverses both masks.  This is the configuration-side count, kept
    independent of the unordered-pair enumeration above.
    """
    canonicals: set[tuple[int, int]] = set()
    for top in subsets:
        for bottom in subsets:
            images = (
                (top, bottom),
                (bottom, top),
                (reverse_mask(top, L), reverse_mask(bottom, L)),
                (reverse_mask(bottom, L), reverse_mask(top, L)),
            )
            canonicals.add(min(images))
    return len(canonicals)


def structure_record(L: int, m: int) -> dict:
    subsets = subset_masks(L, m)
    n = len(subsets)
    fixed_explicit = sum(reverse_mask(mask, L) == mask for mask in subsets)
    fixed_formula = fixed_subset_count(L, m)
    pair_counts = reflection_pair_counts(L, m, subsets)

    symmetric_square_dim = n * (n + 1) // 2
    fixed_pair_formula = (fixed_formula * fixed_formula + n) // 2
    reflection_plus_dim = (n * n + 2 * n + fixed_formula * fixed_formula) // 4
    reflection_minus_dim = (n * n - fixed_formula * fixed_formula) // 4
    defect = candidate_defect(L, m)
    candidate_dim = symmetric_square_dim - defect
    balanced_orbits = balanced_orbit_count(L, m, subsets)
    top = top_wedge_record(L, m, subsets)

    formula_matches_explicit = (
        n == comb(L, m)
        and fixed_explicit == fixed_formula
        and pair_counts["unordered_pair_count"] == symmetric_square_dim
        and pair_counts["fixed_unordered_pair_count"] == fixed_pair_formula
        and pair_counts["reversal_orbit_count"] == reflection_plus_dim
        and balanced_orbits == reflection_plus_dim
        and reflection_plus_dim + reflection_minus_dim == symmetric_square_dim
        and top["defect_matches_functional_dim"]
    )
    return {
        "L": L,
        "m": m,
        "particle_sector": 2 * m,
        "subset_dimension": n,
        "reflection_fixed_subsets_explicit": fixed_explicit,
        "reflection_fixed_subsets_formula": fixed_formula,
        "symmetric_square_dimension": symmetric_square_dim,
        "candidate_defect": defect,
        "candidate_dimension": candidate_dim,
        "reflection_fixed_unordered_pairs": pair_counts["fixed_unordered_pair_count"],
        "reflection_plus_dimension": reflection_plus_dim,
        "reflection_minus_dimension": reflection_minus_dim,
        "balanced_configuration_orbit_count": balanced_orbits,
        "balanced_polarization_rank": reflection_plus_dim,
        "balanced_polarization_kernel": reflection_minus_dim,
        "natural_equivariant_injection_possible": reflection_minus_dim == 0,
        "top_wedge": top,
        "formula_matches_explicit": formula_matches_explicit,
        "tag": "[THEOREM] dimensions of the canonical symmetric-square candidate and the balanced polarization",
    }


def build_structure_records(L_values: tuple[int, ...] = SMALL_L) -> list[dict]:
    return [structure_record(L, m) for L in L_values for m in range(L + 1)]


def main() -> int:
    records = build_structure_records()
    all_formulas = all(row["formula_matches_explicit"] for row in records)
    all_interior_obstructed = all(
        row["balanced_polarization_kernel"] > 0
        for row in records
        if 0 < row["m"] < row["L"]
    )
    delta_rows = [
        [row["L"], row["m"]]
        for row in records
        if row["top_wedge"]["symmetric_invariant_functional_dim"] == 1
    ]
    print(json.dumps({
        "status": "[THEOREM] exact small checks support the proved formulas; [UNRESOLVED] no physical U_L isomorphism",
        "all_formula_checks_passed": all_formulas,
        "all_interior_polarization_kernels_nonzero": all_interior_obstructed,
        "delta_rows_L2_9": delta_rows,
        "record_count": len(records),
    }, indent=2, sort_keys=True))
    return 0 if all_formulas and all_interior_obstructed else 1


if __name__ == "__main__":
    raise SystemExit(main())
