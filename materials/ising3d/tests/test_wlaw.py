#!/usr/bin/env python3
"""Clean-room verifier for the W-law structural route-closure certificate.

Independently re-derives every stored number in results/algebra_growth/wlaw.json:
subset/pair combinatorics from a Pascal table and generating functions, wedge
signs by explicit sort parity, hard-core leakage in a block (non-interleaved)
site encoding, sector dimensions by direct 4^L orbit sweeps, and the seven
held-out controls straight from their source artifacts.

This file must not import experiments/e191*, e192*, or e193*.

Run:
    .venv/bin/python tests/test_wlaw.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "wlaw.json"
REGRESSION = ROOT / "results" / "algebra_growth" / "ladder_w8_regression.json"
W9 = ROOT / "results" / "ladder" / "w9_saturation.json"

SMALL_L = tuple(range(2, 10))
EXPECTED_CHECK_NAMES = {
    "C1_structure_formulas_exact_L2_9",
    "C2_top_wedge_delta_location",
    "C3_reflection_kernel_and_strict_ceiling",
    "C4_earliest_exact_countercertificate",
    "C5_module_leakage_D",
    "C5_module_leakage_F",
    "C5_module_leakage_D_plus",
    "C6_exceptional_witness_is_traceless",
    "C7_old_controls_integrity",
    "C8_old_controls_parameter_free_match",
    "C9_resource_ceiling",
}

CHECKS = 0
FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    global CHECKS
    CHECKS += 1
    tag = "ok" if condition else "FAIL"
    print(f"[test_wlaw] {tag:4} {name} {detail}")
    if not condition:
        FAILURES.append(name)


# ---------------------------------------------------------------------------
# Independent combinatorics (no math.comb: Pascal table and generating funcs).
# ---------------------------------------------------------------------------
def pascal(limit: int) -> list[list[int]]:
    table = [[1]]
    for row_index in range(1, limit + 1):
        previous = table[-1]
        row = [1]
        for k in range(1, row_index):
            row.append(previous[k - 1] + previous[k])
        row.append(1)
        table.append(row)
    return table


PASCAL = pascal(64)


def binom(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    return PASCAL[n][k]


def fixed_subsets_genfunc(L: int, m: int) -> int:
    """[x^m] (1+x)^(L mod 2) (1+x^2)^(L//2), by polynomial convolution."""
    poly = [1]
    for _ in range(L % 2):
        nxt = [0] * (len(poly) + 1)
        for degree, coefficient in enumerate(poly):
            nxt[degree] += coefficient
            nxt[degree + 1] += coefficient
        poly = nxt
    for _ in range(L // 2):
        nxt = [0] * (len(poly) + 2)
        for degree, coefficient in enumerate(poly):
            nxt[degree] += coefficient
            nxt[degree + 2] += coefficient
        poly = nxt
    return poly[m] if m < len(poly) else 0


def reverse_subset(subset: tuple[int, ...], L: int) -> tuple[int, ...]:
    return tuple(sorted(L - 1 - r for r in subset))


def wedge_sign_sorted(left: tuple[int, ...], right: tuple[int, ...], L: int) -> int:
    """Sign of the permutation sorting list(left)+list(right); 0 if degenerate."""
    if set(left) & set(right) or len(left) + len(right) != L:
        return 0
    sequence = list(left) + list(right)
    parity = 0
    for i in range(len(sequence)):
        for j in range(i + 1, len(sequence)):
            if sequence[i] > sequence[j]:
                parity ^= 1
    return -1 if parity else 1


def delta_defect(L: int, m: int) -> int:
    return int(L % 4 == 0 and 2 * m == L)


def T_sym(L: int, m: int) -> int:
    n = binom(L, m)
    return n * (n + 1) // 2


# ---------------------------------------------------------------------------
# Independent ladder configurations: block encoding, top r -> r, bottom -> L+r.
# ---------------------------------------------------------------------------
def block_edges(L: int) -> list[tuple[int, int]]:
    edges = [(r, L + r) for r in range(L)]
    edges += [(r, r + 1) for r in range(L - 1)]
    edges += [(L + r, L + r + 1) for r in range(L - 1)]
    return edges


def block_tau(config: frozenset[int], L: int) -> frozenset[int]:
    return frozenset(s + L if s < L else s - L for s in config)


def block_rho(config: frozenset[int], L: int) -> frozenset[int]:
    return frozenset(L - 1 - s if s < L else L + (L - 1 - (s - L)) for s in config)


def block_orbit(config: frozenset[int], L: int) -> set[frozenset[int]]:
    reflected = block_rho(config, L)
    return {config, block_tau(config, L), reflected, block_tau(reflected, L)}


def block_leg_counts(config: frozenset[int], L: int) -> tuple[int, int]:
    top = sum(1 for s in config if s < L)
    return top, len(config) - top


def block_apply(vector: dict[frozenset[int], int], L: int, piece: str) -> dict[frozenset[int], int]:
    output: dict[frozenset[int], int] = {}
    for config, coefficient in vector.items():
        for u, v in block_edges(L):
            occupancy = (u in config) + (v in config)
            enabled = (
                (piece == "D" and occupancy == 0)
                or (piece == "F" and occupancy == 1)
                or (piece == "D_plus" and occupancy == 2)
            )
            if enabled:
                image = frozenset(config ^ {u, v})
                output[image] = output.get(image, 0) + coefficient
    return {config: coefficient for config, coefficient in output.items() if coefficient}


def block_transform(vector: dict[frozenset[int], int], L: int, which: str) -> dict[frozenset[int], int]:
    fn = block_tau if which == "tau" else block_rho
    output: dict[frozenset[int], int] = {}
    for config, coefficient in vector.items():
        image = fn(config, L)
        output[image] = output.get(image, 0) + coefficient
    return output


def to_interleaved(config: frozenset[int], L: int) -> int:
    encoded = 0
    for s in config:
        encoded |= 1 << (2 * s) if s < L else 1 << (2 * (s - L) + 1)
    return encoded


def config_orbit_census(L: int) -> tuple[dict[int, int], dict[int, int]]:
    """Direct sweep over all 4^L configurations.

    Returns (even-sector orbit counts by particle number,
             balanced (m,m) orbit counts by m).
    """
    def reverse_mask(mask: int) -> int:
        out = 0
        for r in range(L):
            if (mask >> r) & 1:
                out |= 1 << (L - 1 - r)
        return out

    sector_seen: set[tuple[int, int]] = set()
    balanced_seen: set[tuple[int, int]] = set()
    sector_counts: dict[int, int] = {}
    balanced_counts: dict[int, int] = {}
    for top in range(1 << L):
        top_bits = top.bit_count()
        top_reversed = reverse_mask(top)
        for bottom in range(1 << L):
            bottom_bits = bottom.bit_count()
            total = top_bits + bottom_bits
            if total % 2:
                continue
            bottom_reversed = reverse_mask(bottom)
            canonical = min(
                (top, bottom),
                (bottom, top),
                (top_reversed, bottom_reversed),
                (bottom_reversed, top_reversed),
            )
            if canonical not in sector_seen:
                sector_seen.add(canonical)
                sector_counts[total] = sector_counts.get(total, 0) + 1
            if top_bits == bottom_bits and canonical not in balanced_seen:
                balanced_seen.add(canonical)
                balanced_counts[top_bits] = balanced_counts.get(top_bits, 0) + 1
    return sector_counts, balanced_counts


# ---------------------------------------------------------------------------
# V1: envelope schema, meta, provenance hashes, stored check table.
# ---------------------------------------------------------------------------
def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def stage_V1(envelope: dict) -> None:
    check("V1_top_level_schema", set(envelope) == {"meta", "data", "checks"},
          "envelope has exactly meta/data/checks")
    meta = envelope["meta"]
    check("V1_meta_script", meta.get("script") == "experiments/e193_wlaw_certificate.py")
    check("V1_meta_rss_budget", int(meta.get("peak_rss_bytes", 1 << 62)) < 2 * 1024**3,
          f"peak_rss_bytes={meta.get('peak_rss_bytes')}")
    sources = meta.get("source_artifacts", {})
    expected_sources = {
        "results/algebra_growth/ladder_w8_regression.json": REGRESSION,
        "results/ladder/w9_saturation.json": W9,
    }
    hashes_ok = set(sources) == set(expected_sources) and all(
        sources[name] == sha256_of(path) for name, path in expected_sources.items()
    )
    check("V1_source_artifact_hashes", hashes_ok,
          "stored sha256 of both source artifacts matches the files on disk")
    stored = envelope["checks"]
    names = [row["name"] for row in stored]
    check("V1_check_names", set(names) == EXPECTED_CHECK_NAMES and len(names) == len(set(names)),
          f"{len(names)} stored checks")
    check("V1_all_stored_checks_passed", all(row["passed"] is True for row in stored))
    data = envelope["data"]
    check("V1_record_counts",
          len(data["structure_records_exact_L2_9"]) == 52
          and len(data["module_records_exact_L2_9"]) == 24
          and len(data["old_exact_controls_L3_9"]) == 7)
    check("V1_law_left_unresolved",
          "[UNRESOLVED]" in data["headline"]["law_status"]
          and any("[UNRESOLVED]" in row for row in data["scope"]),
          "artifact does not promote the rank law beyond finite evidence")


# ---------------------------------------------------------------------------
# V2: structure records, all fields re-derived independently.
# ---------------------------------------------------------------------------
def recompute_structure(L: int, m: int, balanced_census: dict[int, int]) -> dict:
    subsets = list(combinations(range(L), m))
    n = len(subsets)
    fixed_explicit = sum(reverse_subset(s, L) == s for s in subsets)
    fixed_formula = fixed_subsets_genfunc(L, m)
    sym_dim = n * (n + 1) // 2

    fixed_pairs = 0
    for i in range(n):
        for j in range(i, n):
            pair = tuple(sorted((subsets[i], subsets[j])))
            reflected = tuple(sorted((reverse_subset(subsets[i], L), reverse_subset(subsets[j], L))))
            if pair == reflected:
                fixed_pairs += 1
    plus_dim = (sym_dim + fixed_pairs) // 2  # Burnside over {1, rho}
    minus_dim = sym_dim - plus_dim

    top_degree = 2 * m == L
    admissible = 0
    nonzero_symmetrized = 0
    if top_degree:
        universe = set(range(L))
        for i in range(n):
            complement = tuple(sorted(universe - set(subsets[i])))
            if len(complement) == m and subsets[i] < complement:
                admissible += 1
                symmetrized = (
                    wedge_sign_sorted(subsets[i], complement, L)
                    + wedge_sign_sorted(complement, subsets[i], L)
                )
                if symmetrized:
                    nonzero_symmetrized += 1
    functional_dim = int(top_degree and m % 2 == 0)
    defect = delta_defect(L, m)

    return {
        "L": L,
        "m": m,
        "particle_sector": 2 * m,
        "subset_dimension": n,
        "reflection_fixed_subsets_explicit": fixed_explicit,
        "reflection_fixed_subsets_formula": fixed_formula,
        "symmetric_square_dimension": sym_dim,
        "candidate_defect": defect,
        "candidate_dimension": sym_dim - defect,
        "reflection_fixed_unordered_pairs": fixed_pairs,
        "reflection_plus_dimension": plus_dim,
        "reflection_minus_dimension": minus_dim,
        "balanced_configuration_orbit_count": balanced_census[m],
        "balanced_polarization_rank": plus_dim,
        "balanced_polarization_kernel": minus_dim,
        "natural_equivariant_injection_possible": minus_dim == 0,
        "top_wedge_top_degree": top_degree,
        "top_wedge_exchange_sign": ((-1) ** m if top_degree else None),
        "top_wedge_admissible": admissible,
        "top_wedge_nonzero_symmetrized": nonzero_symmetrized,
        "top_wedge_functional_dim": functional_dim,
    }


def stage_V2(envelope: dict, balanced_by_L: dict[int, dict[int, int]]) -> None:
    rows = envelope["data"]["structure_records_exact_L2_9"]
    all_fields_ok = True
    closed_forms_ok = True
    first_bad = ""
    for stored in rows:
        L, m = stored["L"], stored["m"]
        mine = recompute_structure(L, m, balanced_by_L[L])
        n, f = mine["subset_dimension"], mine["reflection_fixed_subsets_formula"]
        closed_forms_ok = closed_forms_ok and (
            n == binom(L, m)
            and mine["reflection_fixed_unordered_pairs"] == (f * f + n) // 2
            and mine["reflection_plus_dimension"] == (n * n + 2 * n + f * f) // 4
            and mine["reflection_minus_dimension"] == (n * n - f * f) // 4
            and mine["balanced_configuration_orbit_count"] == mine["reflection_plus_dimension"]
        )
        top = stored["top_wedge"]
        row_ok = (
            all(stored[key] == mine[key] for key in (
                "particle_sector", "subset_dimension",
                "reflection_fixed_subsets_explicit", "reflection_fixed_subsets_formula",
                "symmetric_square_dimension", "candidate_defect", "candidate_dimension",
                "reflection_fixed_unordered_pairs", "reflection_plus_dimension",
                "reflection_minus_dimension", "balanced_configuration_orbit_count",
                "balanced_polarization_rank", "balanced_polarization_kernel",
                "natural_equivariant_injection_possible",
            ))
            and stored["formula_matches_explicit"] is True
            and top["top_degree"] == mine["top_wedge_top_degree"]
            and top["exchange_sign"] == mine["top_wedge_exchange_sign"]
            and top["torus_weight_admissible_unordered_pairs"] == mine["top_wedge_admissible"]
            and top["nonzero_symmetrized_wedge_coefficients"] == mine["top_wedge_nonzero_symmetrized"]
            and top["symmetric_invariant_functional_dim"] == mine["top_wedge_functional_dim"]
            and top["defect_matches_functional_dim"]
            is (mine["candidate_defect"] == mine["top_wedge_functional_dim"])
        )
        if not row_ok and not first_bad:
            first_bad = f"(L={L}, m={m})"
        all_fields_ok = all_fields_ok and row_ok
    check("V2_structure_rows_recomputed", all_fields_ok and len(rows) == 52,
          first_bad or "all 52 rows agree field-by-field")
    check("V2_closed_forms", closed_forms_ok,
          "binomial, fixed-pair, plus/minus, and balanced-orbit closed forms all agree")
    delta_rows = [(row["L"], row["m"]) for row in rows
                  if row["top_wedge"]["symmetric_invariant_functional_dim"] == 1]
    check("V2_delta_rows", delta_rows == [(4, 2), (8, 4)], f"delta rows {delta_rows}")
    interior = [row for row in rows if 0 < row["m"] < row["L"]]
    check("V2_interior_kernel_and_ceiling",
          all(row["balanced_polarization_kernel"] > 0
              and row["balanced_polarization_rank"] < row["candidate_dimension"]
              for row in interior),
          f"{len(interior)} interior rows")


# ---------------------------------------------------------------------------
# V3: module leakage records re-derived in the block encoding.
# ---------------------------------------------------------------------------
def recompute_module(L: int, piece: str) -> dict:
    if piece == "D":
        seed = frozenset()
    elif piece == "F":
        seed = frozenset({0, L + 1})
    elif piece == "D_plus":
        seed = frozenset({0, 1, L, L + 1})
    else:
        raise ValueError(piece)
    source = {config: 1 for config in block_orbit(seed, L)}
    output = block_apply(source, L, piece)
    leakage = {
        config: coefficient
        for config, coefficient in output.items()
        if block_leg_counts(config, L)[0] != block_leg_counts(config, L)[1]
    }
    grade = {"D": 2, "F": 0, "D_plus": -2}[piece]
    source_weights = {len(config) for config in source}
    output_weights = {len(config) for config in output}
    examples = sorted(
        (to_interleaved(config, L),) + block_leg_counts(config, L) + (coefficient,)
        for config, coefficient in leakage.items()
    )[:4]
    checks = {
        "source_balanced": all(
            block_leg_counts(config, L)[0] == block_leg_counts(config, L)[1] for config in source
        ),
        "source_tau_rho_invariant": (
            block_transform(source, L, "tau") == source
            and block_transform(source, L, "rho") == source
        ),
        "output_tau_rho_invariant": (
            block_transform(output, L, "tau") == output
            and block_transform(output, L, "rho") == output
        ),
        "leakage_tau_rho_invariant": (
            block_transform(leakage, L, "tau") == leakage
            and block_transform(leakage, L, "rho") == leakage
        ),
        "grade_exact": output_weights == {weight + grade for weight in source_weights},
        "leakage_nonzero": bool(leakage),
    }
    trace = None
    if L == 4 and piece == "D_plus":
        trace = (
            wedge_sign_sorted((0, 1), (0, 1), 4)
            + wedge_sign_sorted((2, 3), (2, 3), 4)
        )
    return {
        "source_support": len(source),
        "source_particle_numbers": sorted(source_weights),
        "output_support": len(output),
        "output_particle_numbers": sorted(output_weights),
        "unbalanced_leakage_support": len(leakage),
        "unbalanced_leakage_l1": sum(abs(coefficient) for coefficient in leakage.values()),
        "examples": examples,
        "checks": checks,
        "trace": trace,
    }


def stage_V3(envelope: dict) -> None:
    rows = envelope["data"]["module_records_exact_L2_9"]
    all_ok = True
    first_bad = ""
    seen = set()
    for stored in rows:
        L, piece = stored["L"], stored["piece"]
        seen.add((L, piece))
        mine = recompute_module(L, piece)
        stored_examples = [
            (row["configuration"], row["top_count"], row["bottom_count"], row["coefficient"])
            for row in stored["unbalanced_leakage_examples"]
        ]
        expected_trace_check = mine["trace"] == 0 if mine["trace"] is not None else True
        row_ok = (
            stored["source_support"] == mine["source_support"]
            and stored["source_particle_numbers"] == mine["source_particle_numbers"]
            and stored["output_support"] == mine["output_support"]
            and stored["output_particle_numbers"] == mine["output_particle_numbers"]
            and stored["unbalanced_leakage_support"] == mine["unbalanced_leakage_support"]
            and stored["unbalanced_leakage_l1"] == mine["unbalanced_leakage_l1"]
            and stored_examples == mine["examples"]
            and stored["exceptional_top_wedge_trace"] == mine["trace"]
            and all(stored["checks"][key] is value for key, value in mine["checks"].items())
            and stored["checks"]["exceptional_D_plus_source_is_traceless"] is expected_trace_check
            and stored["all_checks_passed"] is True
            and mine["unbalanced_leakage_support"] > 0
        )
        if not row_ok and not first_bad:
            first_bad = f"(L={L}, {piece})"
        all_ok = all_ok and row_ok
    check("V3_module_rows_recomputed",
          all_ok and seen == {(L, p) for L in SMALL_L for p in ("D", "F", "D_plus")},
          first_bad or "all 24 rows agree, leakage nonzero everywhere")
    exceptional = next(row for row in rows if row["L"] == 4 and row["piece"] == "D_plus")
    check("V3_exceptional_traceless_leak",
          exceptional["exceptional_top_wedge_trace"] == 0
          and exceptional["unbalanced_leakage_support"] > 0,
          "the 4|L middle-sector witness is traceless and still leaks")


# ---------------------------------------------------------------------------
# V4: the seven held-out controls, re-derived from the source artifacts.
# ---------------------------------------------------------------------------
def control_rows_from_sources() -> tuple[list[dict], bool]:
    regression = json.loads(REGRESSION.read_text())
    w9 = json.loads(W9.read_text())
    rows = []
    for record in regression["data"]["regression_closures"]:
        rows.append({
            "L": int(record["L"]),
            "cyclic": int(record["cyclic_dim"]),
            "K": int(record["K_dim"]),
            "W": int(record["W_dim"]),
            "low_ranks": {int(k): int(v) for k, v in record["low_ranks"].items()},
            "deficits": {int(k): int(v) for k, v in record["W_sectors"].items()},
        })
    data = w9["data"]
    rows.append({
        "L": 9,
        "cyclic": int(data["W9"]["cyclic_Q_dim"]),
        "K": int(data["K9_dimension_record"]["K_dim_formula"]),
        "W": int(data["W9"]["value"]),
        "low_ranks": {int(k): int(v) for k, v in data["closures"]["q1"]["low_ranks"].items()},
        "deficits": {int(k): int(v) for k, v in data["W9"]["W_sectors"].items()},
    })
    source_ok = (
        all(item["passed"] for item in regression["data"]["checks"])
        and all(item["passed"] for item in w9["checks"])
        and data["closures"]["q1"]["low_ranks"] == data["closures"]["p1"]["low_ranks"]
    )
    return rows, source_ok


def stage_V4(envelope: dict, sector_census_by_L: dict[int, dict[int, int]]) -> None:
    stored_rows = {row["L"]: row for row in envelope["data"]["old_exact_controls_L3_9"]}
    source_rows, source_ok = control_rows_from_sources()
    check("V4_source_artifact_internal_checks", source_ok,
          "both source artifacts pass their own stored checks; w9 two-prime ranks agree")

    all_ok = True
    match_count = 0
    first_bad = ""
    for source in source_rows:
        L = source["L"]
        stored = stored_rows[L]
        ranks = {}
        for sector, value in source["low_ranks"].items():
            ranks[sector] = value
            ranks[2 * L - sector] = value
        deficits = {sector: source["deficits"].get(sector, 0) for sector in range(0, 2 * L + 1, 2)}
        K_sectors = {sector: ranks[sector] + deficits[sector] for sector in ranks}
        predicted = {2 * m: T_sym(L, m) - delta_defect(L, m) for m in range(L + 1)}
        census = sector_census_by_L[L]

        integrity = (
            sum(ranks.values()) == source["cyclic"]
            and sum(deficits.values()) == source["W"]
            and sum(K_sectors.values()) == source["K"]
            and source["cyclic"] + source["W"] == source["K"]
            and source["K"] == 2 ** (2 * L - 3) + 3 * 2 ** (L - 2)
            and all(census.get(sector, 0) == dim for sector, dim in K_sectors.items())
        )
        row_ok = (
            stored["cyclic_dimension"] == source["cyclic"]
            and stored["annihilator_dimension"] == source["W"]
            and stored["K_dimension"] == source["K"]
            and {int(k): v for k, v in stored["certified_rank_by_sector"].items()} == ranks
            and {int(k): v for k, v in stored["certified_deficit_by_sector"].items()} == deficits
            and {int(k): v for k, v in stored["certified_K_dimension_by_sector"].items()} == K_sectors
            and {int(k): v for k, v in stored["parameter_free_candidate_rank_by_sector"].items()}
            == predicted
            and stored["integrity_passed"] is True
            and stored["all_sector_matches"] is True
            and integrity
        )
        sector_matches = {sector: ranks[sector] == predicted[sector] for sector in ranks}
        row_ok = row_ok and (
            {int(k): v for k, v in stored["sector_matches"].items()} == sector_matches
            and all(sector_matches.values())
        )
        match_count += len(sector_matches)
        if not row_ok and not first_bad:
            first_bad = f"L={L}"
        all_ok = all_ok and row_ok
    check("V4_controls_recomputed", all_ok, first_bad or "all 7 control rows agree with sources")
    check("V4_sector_match_count", match_count == 49
          and envelope["data"]["control_protocol"]["control_sector_values"] == 49,
          f"{match_count} certified sector values matched, parameter-free")
    protocol = envelope["data"]["control_protocol"]
    check("V4_zero_fit_protocol",
          protocol["fit_parameters"] == 0 and protocol["fit_totals"] == 0
          and protocol["fit_sector_values"] == 0 and protocol["control_totals"] == 7)


# ---------------------------------------------------------------------------
# V5: earliest countercertificate and stored check table recomputation.
# ---------------------------------------------------------------------------
def stage_V5(envelope: dict) -> None:
    earliest = envelope["data"]["earliest_countercertificate"]
    regression = json.loads(REGRESSION.read_text())
    l3 = next(r for r in regression["data"]["regression_closures"] if r["L"] == 3)
    certified_rank_2 = int(l3["low_ranks"]["2"])
    n, f = binom(3, 1), fixed_subsets_genfunc(3, 1)
    ceiling = (n * n + 2 * n + f * f) // 4
    kernel = (n * n - f * f) // 4
    candidate = T_sym(3, 1) - delta_defect(3, 1)
    check("V5_earliest_arithmetic",
          n == 3 and f == 1 and ceiling == 4 and kernel == 2 and candidate == 6
          and certified_rank_2 == 6 and ceiling < certified_rank_2)
    check("V5_earliest_record",
          earliest["L"] == 3 and earliest["m"] == 1 and earliest["particle_sector"] == 2
          and earliest["candidate_dimension"] == candidate
          and earliest["natural_reflection_equivariant_rank_ceiling"] == ceiling
          and earliest["forced_kernel_dimension"] == kernel
          and earliest["certified_U_sector_rank"] == certified_rank_2)

    # Interior strict ceiling, formula-level, well beyond the enumerated window.
    strict = True
    for L in range(2, 31):
        for m in range(1, L):
            n_lm = binom(L, m)
            f_lm = fixed_subsets_genfunc(L, m)
            plus = (n_lm * n_lm + 2 * n_lm + f_lm * f_lm) // 4
            if plus >= T_sym(L, m) - delta_defect(L, m):
                strict = False
    check("V5_interior_ceiling_strict_L2_30", strict,
          "(n^2+2n+f^2)/4 < T(L,m)-delta at every interior (L,m), L=2..30")


def stage_V6(envelope: dict) -> None:
    data = envelope["data"]
    stored = {row["name"]: row["passed"] for row in envelope["checks"]}
    structure_rows = data["structure_records_exact_L2_9"]
    module_rows = data["module_records_exact_L2_9"]
    controls = data["old_exact_controls_L3_9"]
    interior = [row for row in structure_rows if 0 < row["m"] < row["L"]]
    recomputed = {
        "C1_structure_formulas_exact_L2_9": all(r["formula_matches_explicit"] for r in structure_rows),
        "C2_top_wedge_delta_location": [
            (r["L"], r["m"]) for r in structure_rows
            if r["top_wedge"]["symmetric_invariant_functional_dim"] == 1
        ] == [(4, 2), (8, 4)],
        "C3_reflection_kernel_and_strict_ceiling": all(
            r["balanced_polarization_kernel"] > 0
            and r["balanced_polarization_rank"] < r["candidate_dimension"]
            for r in interior
        ),
        "C4_earliest_exact_countercertificate": (
            data["earliest_countercertificate"]["candidate_dimension"] == 6
            and data["earliest_countercertificate"]["certified_U_sector_rank"] == 6
            and data["earliest_countercertificate"]["natural_reflection_equivariant_rank_ceiling"] == 4
            and data["earliest_countercertificate"]["forced_kernel_dimension"] == 2
        ),
        "C6_exceptional_witness_is_traceless": next(
            r for r in module_rows if r["L"] == 4 and r["piece"] == "D_plus"
        )["exceptional_top_wedge_trace"] == 0,
        "C7_old_controls_integrity": len(controls) == 7 and all(r["integrity_passed"] for r in controls),
        "C8_old_controls_parameter_free_match": all(r["all_sector_matches"] for r in controls),
        "C9_resource_ceiling": int(envelope["meta"]["peak_rss_bytes"]) < 2 * 1024**3,
    }
    for piece in ("D", "F", "D_plus"):
        recomputed[f"C5_module_leakage_{piece}"] = all(
            r["all_checks_passed"] and r["unbalanced_leakage_support"] > 0
            for r in module_rows if r["piece"] == piece
        )
    agree = all(stored[name] is bool(value) or stored[name] == bool(value)
                for name, value in recomputed.items())
    check("V6_stored_check_table_recomputed", agree and len(stored) == len(recomputed),
          "every stored check value re-derives from the stored data")


def main() -> int:
    envelope = json.loads(ARTIFACT.read_text())
    stage_V1(envelope)

    sector_census_by_L: dict[int, dict[int, int]] = {}
    balanced_by_L: dict[int, dict[int, int]] = {}
    for L in SMALL_L:
        sector_census_by_L[L], balanced_by_L[L] = config_orbit_census(L)

    stage_V2(envelope, balanced_by_L)
    stage_V3(envelope)
    stage_V4(envelope, sector_census_by_L)
    stage_V5(envelope)
    stage_V6(envelope)

    print(f"[test_wlaw] {CHECKS} checks, {len(FAILURES)} failures")
    if FAILURES:
        print("[test_wlaw] FAILED: " + ", ".join(FAILURES))
        return 1
    print("[test_wlaw] PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
