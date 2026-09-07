# -*- coding: utf-8 -*-
"""Compact regression test for the qlops reproduction layer (pytest-style,
also runnable as a script).  Defends the observable contracts of gate A:

  1. Table-1 component sums recover the printed SEC cycle times exactly.
  2. Eq.(1)/(2) reproduce the published Table-5 neutral-atom numbers to <=
     0.05% (paper prints 4-6 sig figs).
  3. RSA-2048 Sec. 3.5 headline numbers reproduce exactly under exact-
     rational ceiling; ceil on floats must NOT be used (documents the bug).
  4. Every Litinski 15-to-1 Table-6 unit-qubit count matches the closed
     formula exactly, and Total = unit * integer.
  5. Gate-B latency bound: for m<=2, |Q(m tr)/Q(tr) - 1| < 2 strictly
     (structural property of Eq. 2's ceil).

Run: python3 src/test_reproduction.py  (or pytest src/test_reproduction.py)
"""

import math
import sys
from decimal import Decimal as D
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import paper_data as P  # noqa: E402
import formulas as F    # noqa: E402


def _q(k, d, t_sec, t_r):
    r = math.ceil(D(t_r) / D(t_sec))
    return int(k) / ((r + int(d)) * float(D(t_sec)))


def test_sec_component_sums():
    cur = (P.SC_HW["current"]["t_prep"] + P.SC_HW["current"]["t_ro"]
           + 4 * P.SC_HW["current"]["t_1q"] + 4 * P.SC_HW["current"]["t_2q"])
    fut = (P.SC_HW["future"]["t_prep"] + P.SC_HW["future"]["t_ro"]
           + 4 * P.SC_HW["future"]["t_1q"] + 4 * P.SC_HW["future"]["t_2q"])
    assert abs(cur - 0.86) < 1e-12, cur
    assert abs(fut - 0.40) < 1e-12, fut
    assert abs(P.T_SEC_SC["current"] - cur * 1e-6) < 1e-18
    assert abs(P.T_SEC_SC["future"] - fut * 1e-6) < 1e-18


def test_table5_neutral_atom_rows():
    for code, t5 in P.T5.items():
        n, k, d = t5["code"]
        t4 = P.GB_TABLE4[code]
        q = _q(k, d, t4["t_sec"], t4["tr"])
        assert abs(q / t5["na_qlops"] - 1) <= 5e-4, (code, q, t5["na_qlops"])
        p0 = float(F.p0_from_pL(D(t4["pL"]), k, d))
        assert abs(p0 / t5["p0"] - 1) <= 5e-4, (code, p0, t5["p0"])
        assert t5["na_qubits"] == 2 * n  # N = 2n convention


def test_rsa_exact_vs_float_ceil():
    k = P.RSA_SC["n_cold"] + P.RSA_SC["n_active"]
    d, t_sec, t_r = P.RSA_SC["d"], P.RSA_SC["t_sec"], P.RSA_SC["t_r"]
    # exact-rational ceiling:
    r_exact = math.ceil(Fraction(1, 10**5) / Fraction(1, 10**6))
    assert r_exact == 10
    q_exact = k / ((r_exact + d) * 1e-6)
    assert abs(q_exact / P.RSA_SC["qlops"] - 1) <= 1.3e-5
    # float ceiling reproduces the WRONG round count (11) -- the pitfall:
    assert math.ceil(t_r / t_sec) == 11
    q_float = k / ((11 + d) * 1e-6)
    assert abs(q_float / P.RSA_SC["qlops"] - 1) > 0.02  # 2.8% off


def test_rsa_2048_block():
    q_sc = 1411 / ((10 + 25) * 1e-6)
    # exact 40314285.71... ; printed to 5 sig figs -> |dev| <= ~1.24e-5
    assert abs(q_sc / P.RSA_SC["qlops"] - 1) <= 1.3e-5, q_sc
    # and the printed value IS the correct 5-sig-fig rounding:
    assert float(f"{q_sc:.5g}") == P.RSA_SC["qlops"]  # 5-sig-fig round
    data = 1280 * 430 + 131 * (2 * 25**2 - 1)
    assert data == 714019 == P.RSA_SC["data_qubits"]
    dens = q_sc / data
    assert abs(dens - 56.4611) / 56.4611 < 1e-4
    q_na = 6128 / 900e-6
    assert abs(q_na - 6.8089e6) / 6.8089e6 < 1e-5
    na_data = 6128 * (2 * 27**2 - 1)           # 8,932,496
    assert abs(q_na / na_data - 0.762593) / 0.762593 < 1e-4
    # Eq.(3)/(4)
    r3 = (q_sc * 4.96) / (q_na * 5.6)
    assert abs(r3 - 5.244) / 5.244 < 5e-2     # 2.8%; paper rounds inputs
    r4 = (q_sc * 4.96 / 6.5e9) / (q_na * 5.6 / 3.0e9)


def test_table6_litinski_units():
    import re
    for row in P.T6:
        m = re.match(r"15to1_(\d+),(\d+),(\d+)", row["proto"])
        dX, dZ, dm = map(int, m.groups())
        unit = 2 * (dX + 4 * dZ) * 3 * dX + 4 * dm
        assert unit == row["unit"], row
        assert row["total"] % unit == 0, row
        quo = row["total"] // unit
        assert quo >= 1


def test_latency_axis_bounded_by_2():
    # analytic: (ceil(m r)+d)/(ceil(r)+d) <= (2r+1+d)/(r+d) < 2 for m<=2
    for code, t4 in P.GB_TABLE4.items():
        n, k, d = (int(x) for x in code.split("-")[0].split(","))
        r0 = math.ceil(D(t4["tr"]) / D(t4["t_sec"]))
        q0 = k / ((r0 + d) * float(D(t4["t_sec"])))
        for m in (0.5, 2.0):
            rm = math.ceil(D(str(m)) * D(t4["tr"]) / D(t4["t_sec"]))
            qm = k / ((rm + d) * float(D(t4["t_sec"])))
            assert abs(qm / q0 - 1) < 2.0 - 1e-9, (code, m)
    # RSA-SC sanity at m=2 (10->20 us): 1411/(50 us) vs 1411/(45?) ... exact:
    q1 = 1411 / ((10 + 25) * 1e-6)
    q2 = 1411 / ((20 + 25) * 1e-6)
    assert abs((q2 / q1) - (35 / 45)) < 1e-12  # ceil clean at both ends


def test_zero_level_sensitivity_uses_common_units_and_full_envelope():
    """Revision 4: accepted-output comparison uses one physical-qubit×SEC basis."""
    import zero_level_provenance as Z
    s = Z.sensitivity()

    legacy = s["legacy_frozen_gate_b_row"]
    assert abs(legacy["advantage_factor"][0] - 193.3) < 0.1
    assert abs(legacy["advantage_factor"][1] - 29946.5) < 0.1
    assert legacy["usable_for_common_basis_magnitude"] is False

    baseline = s["litinski_accepted_baseline_rows"]
    assert baseline[0]["cycles_including_postselection"] == P.T6[0]["cycles"]
    assert baseline[0]["spacetime_per_accepted_T"] == (
        P.T6[0]["unit"] * P.T6[0]["cycles"])

    variants = s["zero_level_variants"]
    assert set(variants) == {
        "distillation_only_output_patches_omitted",
        "with_output_patches_d3",
        "with_output_patches_d7",
    }
    assert variants["with_output_patches_d3"]["qubits"] == 73
    assert variants["with_output_patches_d7"]["qubits"] == 313
    assert variants["with_output_patches_d7"]["spacetime_per_attempt"] == 939

    acceptance_keys = {
        "p=1e-3_optimistic", "p=1e-3_pessimistic", "p=1e-4"}
    expected_combinations = {
        f"{variant}|accept_{acceptance}"
        for variant in variants
        for acceptance in acceptance_keys
    }
    assert set(s["corrected_combinations"]) == expected_combinations

    worst_key = "with_output_patches_d7|accept_p=1e-3_pessimistic"
    worst = s["worst_case"]
    assert worst["combination"] == worst_key
    assert math.isclose(worst["advantage_factor"][0], 33.76693322683706)
    assert math.isclose(worst["advantage_factor"][1], 5086.233185367411)
    assert worst["advantage_factor"][0] == min(
        combo["advantage_factor"][0]
        for combo in s["corrected_combinations"].values())
    assert s["verdict_survives_corrections"] == (
        worst["advantage_factor"][0] > P.GATE_B["comparability_falsifier"])

    assert Z.rotated_patch_qubits(3) == 17
    assert Z.rotated_patch_qubits(7) == 97


if __name__ == "__main__":
    fns = [v for kk, v in sorted(globals().items())
           if kk.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"{len(fns)} tests passed")
