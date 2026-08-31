"""Audit the fixed AI-claim mapping for Gate C.

The campaign does not pretend that the pilot upper bound is an exact global
rho. Instead it records a rigorous universal interval: for every nonzero
M in V, delta(M) <= 3N, hence rho_inst >= 1/(3N). For eta <= 1/32,
eta**(e**3) < 1/90 = 1/(3*30), using exact rational lower bounds on e.
Thus the fixed Gate-C claim is consistent with this pilot at the level of a
proved interval, while the exact rho value remains OPEN.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import math
import platform
import sys
import uuid
from fractions import Fraction
from pathlib import Path

import mpmath

from .rs import Inst

ROOT = Path(__file__).resolve().parents[1]


def _code_hash() -> str:
    h = hashlib.sha256()
    for p in sorted((ROOT / "src").glob("*.py")):
        h.update(p.name.encode())
        h.update(p.read_bytes())
    return h.hexdigest()


def _exp_lower_bound() -> tuple[Fraction, Fraction, int]:
    # The positive Taylor partial sum is an exact lower bound on e.
    e_lower = sum((Fraction(1, math.factorial(k)) for k in range(9)), Fraction(0))
    assert e_lower > Fraction(2718, 1000)
    e3_lower = e_lower**3
    assert e3_lower > 20
    return e_lower, e3_lower, 20


def run() -> Path:
    q = 31
    s = (2, 3, 5)
    N = math.prod(s)
    e_lower, e3_lower, exponent_lower = _exp_lower_bound()
    lower_rho = Fraction(1, 3 * N)
    scale = 10**60
    rows = []
    for eta_den in (32, 64, 128):
        eta = Fraction(1, eta_den)
        # mpmath is display-only; the exact comparison uses the rational
        # exponent lower bound above.
        mpmath.mp.dps = 100
        display_value = mpmath.power(mpmath.mpf(eta.numerator) / eta.denominator, mpmath.e**3)
        below = Fraction(int(display_value * scale), scale)
        # Since 0<eta<1 and e^3>20, eta^(e^3)<eta^20.
        exact_upper = eta**exponent_lower
        assert exact_upper < lower_rho
        assert below < lower_rho
        rows.append(
            {
                "eta": f"{eta.numerator}/{eta.denominator}",
                "rho_claim_decimal_100dps": mpmath.nstr(display_value, 100),
                "rho_claim_floor_10^60": f"{below.numerator}/{below.denominator}",
                "strict_exact_upper_from_e3_gt_20": f"{exact_upper.numerator}/{exact_upper.denominator}",
                "universal_rho_lower_bound": f"{lower_rho.numerator}/{lower_rho.denominator}",
                "claim_le_universal_lower_bound": True,
                "verdict": "consistent-with-form on proved interval",
                "evidence_label": "COMPUTATIONAL-EVIDENCE",
            }
        )

    now = dt.datetime.now(dt.timezone.utc).replace(microsecond=0)
    name = f"{now.strftime('%Y-%m-%dT%H-%M-%SZ')}_{uuid.uuid4().hex[:8]}_gateC"
    out = ROOT / "campaigns" / name
    out.mkdir(parents=True, exist_ok=False)
    manifest = {
        "campaign": name,
        "utc": now.isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "q": q,
        "s": list(s),
        "N": N,
        "mapping": {
            "epsilon": "eta",
            "k": 3,
            "exp_k": "e^3",
            "source_status": "AI candidate proof explicitly unverified",
        },
        "exact_bound": {
            "e_lower_partial_sum_terms": 9,
            "e_lower": f"{e_lower.numerator}/{e_lower.denominator}",
            "e3_lower": f"{e3_lower.numerator}/{e3_lower.denominator}",
            "e3_gt": 20,
            "rho_inst_lower": f"{lower_rho.numerator}/{lower_rho.denominator}",
            "reason": "delta(M) <= N+N+N=3N for every M in V",
        },
        "rows": rows,
        "global_rho_status": "OPEN; pilot has only a finite-candidate upper bound",
        "evidence_label": "COMPUTATIONAL-EVIDENCE",
        "code_sha256": _code_hash(),
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    (out / "gate_c.json").write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n")
    (out / "codehash.txt").write_text(manifest["code_sha256"] + "\n")
    (out / "run.log").write_text(
        "exact e^3 lower bound > 20; exact rho lower bound = 1/90\n"
        + "\n".join(
            f"eta={r['eta']} claim_floor={r['rho_claim_floor_10^60']} verdict={r['verdict']}"
            for r in rows
        )
        + "\n"
    )
    return out


if __name__ == "__main__":
    print(run())
