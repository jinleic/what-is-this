"""Acceptance checks for CRT transfer propagation and extended Ising series.

Run from the repository root:
    .venv/bin/python tests/test_series_extend.py
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import comb
from pathlib import Path

import mpmath as mp

from ising.series.analysis import (
    d_finiteness_scan,
    dlog_pade_scan,
    exponential_series,
    first_order_differential_approximant,
)
from ising.transfer_matrix import (
    box_broken_bond_poly as int64_box_broken_bond_poly,
)
from ising.transfer_matrix import (
    torus_broken_bond_poly as int64_torus_broken_bond_poly,
)
from ising.transfer_matrix.crt import (
    box_broken_bond_poly as crt_box_broken_bond_poly,
)
from ising.transfer_matrix.crt import (
    primes_for_sites,
    torus_broken_bond_poly as crt_torus_broken_bond_poly,
)

ROOT = Path(__file__).resolve().parents[1]
SERIES_DIR = ROOT / "results" / "series"


def _load(name: str) -> dict:
    return json.loads((SERIES_DIR / name).read_text(encoding="utf-8"))


def _fractions(payload: dict) -> tuple[Fraction, ...]:
    return tuple(Fraction(value) for value in payload["data"]["coefficients"])


def main() -> None:
    mp.mp.dps = 80
    comparisons = [
        ("box", (2, 2), {}),
        ("box", (3, 3), {}),
        ("box", (2, 2, 2), {}),
        ("box", (2, 3, 2), {}),
        ("box", (3, 3, 2), {}),
        ("box", (2, 2, 3), {"plus_boundary": True}),
        ("box", (2, 3, 1), {"plus_boundary": True}),
        ("box", (3, 2, 2), {"periodic": (True, False)}),
        ("box", (2, 31), {}),
        ("torus", (4, 4), {}),
        ("torus", (2, 2, 2), {}),
        ("torus", (2, 3, 2), {"block": 2}),
    ]
    for kind, shape, kwargs in comparisons:
        if kind == "box":
            expected = int64_box_broken_bond_poly(shape, **kwargs)
            actual = crt_box_broken_bond_poly(shape, **kwargs)
        else:
            expected = int64_torus_broken_bond_poly(shape, **kwargs)
            actual = crt_torus_broken_bond_poly(shape, **kwargs)
        assert actual == expected, (kind, shape, kwargs)
    assert len(comparisons) >= 8
    assert primes_for_sites(64).product > 1 << 64

    old_ht = _load("sc_ht_free_energy.json")
    old_lt = _load("sc_lt_free_energy.json")
    new_ht = _load("extended_sc_ht_free_energy.json")
    new_lt = _load("extended_sc_lt_free_energy.json")
    assert new_ht["data"]["achieved_order"] >= 20
    assert new_lt["data"]["achieved_order"] >= 28
    assert _fractions(new_ht)[: len(_fractions(old_ht))] == _fractions(old_ht)
    assert _fractions(new_lt)[: len(_fractions(old_lt))] == _fractions(old_lt)

    required_checks = {
        "crt_matches_int64_at_least_8_lattices",
        "sc_ht_enlarged_box_self_consistency",
        "sc_lt_enlarged_box_self_consistency",
        "sc_ht_existing_prefix_exact",
        "sc_lt_existing_prefix_exact",
        "direct_even_subgraph_crosscheck",
    }
    checks = new_ht["checks"] + new_lt["checks"]
    passed = {check["name"] for check in checks if check["passed"]}
    assert required_checks <= passed, required_checks - passed

    direct = new_ht["data"]["independent_direct_crosscheck"]
    assert direct["connected_even_subgraph_counts"] == {
        "4": 3,
        "6": 22,
        "8": 237,
        "10": 3000,
    }
    assert direct["ordered_overlapping_pair_counts"] == {
        "8": 99,
        "10": 2040,
    }

    control = d_finiteness_scan(exponential_series(18), holdout=3)
    assert control.found, "the ODE guesser did not recognize exp(z)"
    assert any(
        fit.order == 1 and fit.degree == 1 and fit.passed for fit in control.fits
    ), control

    algebraic = tuple(
        Fraction(comb(degree + 2, 2) * 2**degree)
        for degree in range(10)
    )
    dlog_control = dlog_pade_scan(algebraic)
    assert any(
        abs(fit.singularity - mp.mpf("0.5")) < mp.mpf("1e-30")
        and abs(fit.exponent - 3) < mp.mpf("1e-30")
        for fit in dlog_control
    ), "Dlog-Pade did not recover (1-2z)^-3"

    branch_series = [Fraction(1)]
    branch_exponent = Fraction(3, 2)
    for degree in range(1, 10):
        branch_series.append(
            branch_series[-1]
            * Fraction(degree - 1 - branch_exponent, degree)
        )
    differential_control = first_order_differential_approximant(
        branch_series, 1, 0, 0
    )
    assert differential_control is not None
    assert abs(differential_control.singularity - 1) < mp.mpf("1e-30")
    assert (
        abs(differential_control.singular_exponent - mp.mpf("1.5"))
        < mp.mpf("1e-30")
    )

    print(
        "PASS test_series_extend: "
        f"CRT={len(comparisons)} lattices; HT=v^{new_ht['data']['achieved_order']}; "
        f"LT=x^{new_lt['data']['achieved_order']}; holonomic control ODE found"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL test_series_extend: {error}")
        raise
