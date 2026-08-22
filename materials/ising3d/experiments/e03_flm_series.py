"""Generate and validate exact finite-lattice Ising free-energy series.

Run from the repository root:
    .venv/bin/python experiments/e03_flm_series.py
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from fractions import Fraction
from math import comb
from pathlib import Path

from ising.lattices import cubic
from ising.series import (
    FLMSeries,
    high_temperature_free_energy,
    low_temperature_free_energy,
    onsager_square_ht_series,
)
from ising.transfer_matrix import torus_broken_bond_poly

SCRIPT = "experiments/e03_flm_series.py"
RESULT_DIR = Path(__file__).resolve().parents[1] / "results" / "series"


def _fraction_strings(values):
    return [str(value) for value in values]


def _nonzero_coefficients(values):
    return {str(degree): str(value) for degree, value in enumerate(values) if value}


def _canonical_boxes(boxes):
    return sorted({tuple(sorted(shape)) for shape in boxes}, key=lambda shape: (sum(shape), shape))


def _count_four_cycles(lattice):
    adjacency = [set() for _ in range(lattice.n_sites)]
    for left, right in lattice.bonds:
        adjacency[left].add(right)
        adjacency[right].add(left)
    opposite_pair_count = sum(
        comb(len(adjacency[left] & adjacency[right]), 2)
        for left in range(lattice.n_sites)
        for right in range(left + 1, lattice.n_sites)
    )
    assert opposite_pair_count % 2 == 0
    return opposite_pair_count // 2


def _record_check(checks, name, passed, detail):
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    if not passed:
        raise AssertionError(f"{name}: {detail}")
    print(f"PASS {name}: {detail}")


def _all_exact(series: FLMSeries):
    return all(isinstance(value, Fraction) for value in series.coefficients) and all(
        isinstance(value, Fraction)
        for weight in series.box_weights.values()
        for value in weight
    )


def _series_payload(series: FLMSeries):
    canonical = _canonical_boxes(series.boxes)
    return {
        "dimension": series.dimension,
        "variable": series.variable,
        "achieved_order": series.order,
        "coefficients": _fraction_strings(series.coefficients),
        "nonzero_coefficients": _nonzero_coefficients(series.coefficients),
        "interaction_coefficients": _fraction_strings(series.interaction_coefficients),
        "nonzero_interaction_coefficients": _nonzero_coefficients(
            series.interaction_coefficients
        ),
        "ordered_box_count": len(series.boxes),
        "box_list": [list(shape) for shape in series.boxes],
        "canonical_box_count": len(canonical),
        "canonical_box_list": [list(shape) for shape in canonical],
        "symmetry_note": (
            "All ordered boxes enter the Moebius inversion; transfer polynomials are cached "
            "by sorted side lengths because the isotropic box polynomial is permutation invariant."
        ),
        "order_bound": {
            "statement": series.order_bound,
            "budget": series.bound_budget,
            "bound_slack": series.bound_slack,
        },
    }


def main():
    timestamp = datetime.now(timezone.utc).isoformat()
    provenance = {
        "script": SCRIPT,
        "timestamp": timestamp,
        "precision": "exact Python int and fractions.Fraction; no floating-point arithmetic",
    }

    square_flm = high_temperature_free_energy(2, 24)
    square_onsager = onsager_square_ht_series(24)
    sc_ht = high_temperature_free_energy(3, 16)
    sc_ht_prefix = high_temperature_free_energy(3, 14)
    sc_ht_enlarged = high_temperature_free_energy(3, 16, bound_slack=1)
    sc_lt = low_temperature_free_energy(3, 20)
    sc_lt_enlarged = low_temperature_free_energy(3, 20, bound_slack=1)

    common_checks = []
    nonzero_square_orders = [
        degree for degree in range(1, 25) if square_onsager[degree]
    ]
    _record_check(
        common_checks,
        "square_onsager_exact_v24",
        square_flm.coefficients == square_onsager
        and nonzero_square_orders == list(range(2, 25, 2)),
        "FLM equals the independently integrated Onsager series at 12 nonzero orders v^2..v^24",
    )

    ht_checks = list(common_checks)
    ht_extra = set(sc_ht_enlarged.boxes) - set(sc_ht.boxes)
    ht_unreachable = [
        shape
        for shape in ht_extra
        if shape[0] * shape[1] * shape[2] > 62
    ]
    _record_check(
        ht_checks,
        "sc_ht_enlarged_box_bound_v16",
        bool(ht_extra)
        and sc_ht.coefficients == sc_ht_enlarged.coefficients
        and all(
            not any(sc_ht_enlarged.box_weights[shape]) for shape in ht_extra
        )
        and ht_unreachable == [(4, 4, 4)],
        f"adding {len(ht_extra)} ordered boxes of coordinate-span budget 9 changes no "
        "coefficient through v^16; 54 are transfer-evaluated and 4x4x4 is zero by "
        "the proven order-18 lower bound",
    )
    _record_check(
        ht_checks,
        "sc_ht_order_extension_v16",
        sc_ht.coefficients[:15] == sc_ht_prefix.coefficients,
        "the v^16 computation reproduces every coefficient through v^14",
    )
    square_lattice = cubic(5, 5, 5, periodic=True)
    four_cycles = _count_four_cycles(square_lattice)
    _record_check(
        ht_checks,
        "sc_ht_v4_direct_cycle_count",
        four_cycles == 3 * square_lattice.n_sites
        and sc_ht.interaction_coefficients[4]
        == Fraction(four_cycles, square_lattice.n_sites),
        f"directly enumerated {four_cycles} four-cycles on the 5^3 torus: 3 per site, matching [v^4] log P = 3",
    )
    _record_check(
        ht_checks,
        "sc_ht_exact_rationals",
        _all_exact(sc_ht),
        "all coefficients and all finite-lattice weights through v^16 are Fraction values",
    )

    torus_reason = ""
    try:
        torus_broken_bond_poly((4, 4, 4), block=1)
    except ValueError as error:
        torus_reason = str(error)
    _record_check(
        ht_checks,
        "torus_4x4x4_reachability_recorded",
        "int64 overflow" in torus_reason,
        "no numerical torus comparison: the exact engine rejects N=64 (" + torus_reason + ")",
    )

    lt_checks = list(common_checks)
    lt_extra = set(sc_lt_enlarged.boxes) - set(sc_lt.boxes)
    _record_check(
        lt_checks,
        "sc_lt_enlarged_box_bound_x20",
        bool(lt_extra)
        and sc_lt.coefficients == sc_lt_enlarged.coefficients
        and all(not any(sc_lt_enlarged.box_weights[shape]) for shape in lt_extra),
        f"adding {len(lt_extra)} ordered boxes of side-sum budget 7 changes no coefficient through x^20",
    )
    _record_check(
        lt_checks,
        "sc_lt_single_spin_and_prefix",
        sc_lt.coefficients[6] == 1
        and sc_lt.coefficients[10] == 3
        and sc_lt.coefficients[12] == Fraction(-7, 2),
        "derived prefix is x^6 + 3*x^10 - (7/2)*x^12; x^6 is one six-bond spin flip per site",
    )
    _record_check(
        lt_checks,
        "sc_lt_exact_rationals",
        _all_exact(sc_lt),
        "all coefficients and all finite-lattice weights through x^20 are Fraction values",
    )

    ht_data = _series_payload(sc_ht)
    ht_data.update(
        {
            "lattice": "simple cubic",
            "normalization": "phi = log(2) + sum_n coefficients[n] * v^n",
            "interaction_normalization": (
                "phi = log(2) + 3*log(cosh K) + "
                "sum_n interaction_coefficients[n] * v^n"
            ),
            "order_bound_proof": (
                "A connected even graph has an Euler circuit.  In each coordinate the circuit "
                "must travel from the minimum to the maximum and back, using at least twice the "
                "coordinate span; summing gives 2*sum_i(side_i-1)."
            ),
            "highest_order_limitation": (
                "v^18 requires the 4x4x4 free box (64 sites).  The public exact transfer engine "
                "rejects boxes above 62 sites because its coefficients use int64."
            ),
            "finite_torus_comparison": {
                "performed": False,
                "reason": "4x4x4 exact torus is unreachable: " + torus_reason,
            },
        }
    )
    lt_data = _series_payload(sc_lt)
    lt_data.update(
        {
            "lattice": "simple cubic",
            "normalization": "phi = 3*K + sum_n coefficients[n] * x^n",
            "order_bound_proof": (
                "For a connected droplet spanning a,b,c, each occupied coordinate column gives "
                "two exposed faces.  Its xy, xz, yz projections have at least a+b-1, a+c-1, "
                "b+c-1 cells, so the surface has at least 4*(a+b+c)-6 broken bonds.  A connected "
                "Mayer cluster has total surface no smaller than the exterior surface of its union."
            ),
        }
    )

    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    outputs = {
        RESULT_DIR / "sc_ht_free_energy.json": {
            "provenance": provenance,
            "data": ht_data,
            "checks": ht_checks,
        },
        RESULT_DIR / "sc_lt_free_energy.json": {
            "provenance": provenance,
            "data": lt_data,
            "checks": lt_checks,
        },
    }
    for path, payload in outputs.items():
        path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(f"WROTE {path.relative_to(RESULT_DIR.parents[1])}")

    print(
        "PASS e03_flm_series: square HT v^24; simple-cubic HT v^16; "
        "simple-cubic LT x^20"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL e03_flm_series: {error}")
        raise
