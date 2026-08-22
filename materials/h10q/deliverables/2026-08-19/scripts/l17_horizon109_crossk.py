#!/usr/bin/env python3
"""Deterministic cross-check probe for the fresh H109 aligned classes.

The probe space is the 25 aligned L11 classes discovered by
l17_horizon109.py, with k=0..50 in the fixed tuple order below.  It records
only engine results (no elapsed times or host-dependent values).  The
FactorBudget/PrimalityBound labels are refusals, never closure evidence.
"""
from __future__ import annotations

import json
import os
import signal
import sys
import time
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import h10q  # noqa: E402
import l13_filter  # noqa: E402

CELL = (109, (-1, 1))
Z = F(-109)
MAX_K = 50
GUARD_SECONDS = 20
SEED = 0
RUN_DATE = "2026-08-20"

# (a, eps, f, q1, N), in the exact order used by the original 1,275-row
# probe.  This list is intentionally self-contained: replay does not consume
# l17_horizon109.jsonl or any other generated artifact.
H109_ALIGNED_CLASSES = (
    ("147", 1, 1, 31, "290428320"),
    ("147", 1, 1, 71, "871284960"),
    ("147", -1, 1, 89, "2032998240"),
    ("147", -1, 1, 193, "2613854880"),
    ("147", 1, 1, 193, "1452141600"),
    ("289", -1, 1, 73, "11449761120"),
    ("289", 1, 1, 73, "34349283360"),
    ("289", -1, 1, 83, "583937817120"),
    ("289", 1, 1, 83, "11449761120"),
    ("289", -1, 1, 113, "34349283360"),
    ("289", 1, 1, 113, "11449761120"),
    ("289", -1, 1, 157, "11449761120"),
    ("289", 1, 1, 157, "34349283360"),
    ("289", -1, 1, 173, "34349283360"),
    ("289", 1, 1, 173, "11449761120"),
    ("289", -1, 1, 193, "11449761120"),
    ("289", 1, 1, 193, "34349283360"),
    ("71", -1, 1, 43, "962112480"),
    ("71", 1, 1, 43, "6734787360"),
    ("71", -1, 1, 97, "6734787360"),
    ("71", 1, 1, 97, "962112480"),
    ("71", -1, 1, 113, "962112480"),
    ("71", 1, 1, 113, "6734787360"),
    ("71", -1, 1, 193, "962112480"),
    ("71", 1, 1, 193, "962112480"),
)

META = {
    "type": "meta",
    "schema": "l17-horizon109-crossk-v1",
    "version": 1,
    "cell": [109, [-1, 1]],
    "family": "L11-fresh",
    "probe_space": {
        "class_count": len(H109_ALIGNED_CLASSES),
        "k_min": 0,
        "k_max": MAX_K,
        "ordering": "embedded tuple order; a=147 then 289 then 71, q1 ascending, eps ascending",
        "guard_seconds": GUARD_SECONDS,
    },
    "seed": SEED,
    "date": RUN_DATE,
    "engine": "math/h10q/h10q.py proven-primality only",
    "nondeterministic_fields": [],
    "determinism_note": "No elapsed, host, or random fields are emitted; all rows are engine verdicts and refusal labels.",
    "row_count_without_meta": len(H109_ALIGNED_CLASSES) * (MAX_K + 1),
}


class ProbeTimeout(Exception):
    pass


def _alarm_handler(signum, frame):
    raise ProbeTimeout()


def guarded(fn):
    old_handler = signal.getsignal(signal.SIGALRM)
    signal.signal(signal.SIGALRM, _alarm_handler)
    signal.alarm(GUARD_SECONDS)
    try:
        return fn(), None
    except Exception as exc:
        return None, "refused:" + type(exc).__name__
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)


def validate_class_list():
    assert len(H109_ALIGNED_CLASSES) == 25
    assert CELL not in h10q._L11_CLASSES
    assert h10q.factorint(abs(Z.numerator)) == {109: 1}
    assert h10q._l11_escape_primes(Z) == [109]
    for a_text, eps, f, q1, N_text in H109_ALIGNED_CLASSES:
        a = F(a_text)
        tau = (1 + 2 * a * a) / (1 + 4 * a * a)
        cert = h10q._l10_class_cert(a, Z, tau, eps, q1)
        assert cert is not None and cert.get("ok")
        assert int(cert["N"]) == int(N_text)
        assert f == 1
        time.sleep(0.1)


def probe():
    h10q._selftest()
    validate_class_list()
    rows = []
    found = []
    for a_text, eps, f, q1, N_text in H109_ALIGNED_CLASSES:
        a = F(a_text)
        tau = (1 + 2 * a * a) / (1 + 4 * a * a)
        N = int(N_text)
        for k in range(MAX_K + 1):
            time.sleep(0.1)
            Q = q1 + k * N
            rec = {
                "a": str(a),
                "eps": eps,
                "f": f,
                "q1": q1,
                "N": str(N),
                "k": k,
                "Q": str(Q),
            }
            prime, err = guarded(lambda: h10q._is_prime(Q))
            rec["prime"] = err if err else bool(prime)
            if err or not prime:
                rows.append(rec)
                continue
            b = F(eps * f * Q)
            tied, err = guarded(lambda: h10q._l7_tied_status(a, b, Z, tau))
            rec["tied"] = err if err else str(tied)
            c0, err = guarded(lambda: h10q._sun_h(a, b, Z**3))
            if err:
                rec["ramified_empty"] = err
            else:
                A = 1 + 4 * a * a
                delta = 1 - A * tau * tau
                alpha = -delta * A
                M0 = 16 - delta * c0 * c0 - 32 * A * b * ((a - 1) / 2) ** 2
                x0, d0 = alpha * M0, alpha * 2 * b
                ramified, err = guarded(lambda: h10q.ramified(x0, d0))
                rec["ramified_empty"] = err if err else ramified == []
                rec["M0_num_digits"] = len(str(abs(M0.numerator)))
                rec["M0_den_digits"] = len(str(M0.denominator))
                if rec["tied"] == "True" and rec["ramified_empty"] is True:
                    smooth, err = guarded(lambda: l13_filter.smooth_emergent(a, Z, tau, b))
                    rec["smooth"] = err if err else smooth[0]
                    if err is None and smooth[0] == "cofactor-big":
                        em, (rn, rd), _ = smooth[1]
                        if not em:
                            ladder, err = guarded(
                                lambda: l13_filter.cofactor_decide(
                                    a, Z, tau, b, detail=smooth[1]
                                )
                            )
                            rec["ladder"] = err if err else ladder[0]
                            if err is None and ladder[0] == "zero":
                                found.append(rec)
                    elif err is None:
                        rec["ladder"] = "not-run"
            rows.append(rec)
            if found:
                break
        if found:
            break
    return rows, found


def main():
    rows, found = probe()
    assert len(rows) == META["row_count_without_meta"] or found
    output = (
        Path(os.environ["H109_CROSSK_OUT"]).expanduser().resolve()
        if "H109_CROSSK_OUT" in os.environ
        else ROOT / "data" / "l17_horizon109_crossk.jsonl"
    )
    with output.open("w", encoding="utf-8", newline="\n") as out:
        out.write(json.dumps(META) + "\n")
        for row in rows:
            out.write(json.dumps(row) + "\n")
    print("WROTE", output)
    print("ROWS", len(rows), "META_PLUS_ROWS", len(rows) + 1)
    print("CROSSCHECKED_CLOSURES", len(found))


if __name__ == "__main__":
    main()
