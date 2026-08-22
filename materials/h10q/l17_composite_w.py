#!/usr/bin/env python3
"""Structural audit and exact arithmetic closure at the pseudo-cell (21, -1).

A composite integer is not a place of Q.  This script first records that
structural rejection, then asks the separate computational question that does
make sense: can the existing prime-place machinery close z=-21 after choosing
an actual prime factor f=3 of z?  Every primality and factorization decision is
made by h10q.py; a refusal is recorded and is never treated as evidence.
"""
from __future__ import annotations

import json
import os
import sys
import time
from fractions import Fraction as F
from pathlib import Path


def find_root() -> Path:
    candidates = []
    if os.environ.get("H10Q_ROOT"):
        candidates.append(Path(os.environ["H10Q_ROOT"]))
    here = Path(__file__).resolve().parent
    candidates.extend((here, Path.cwd(), Path.cwd() / "math" / "h10q"))
    for candidate in candidates:
        if (candidate / "h10q.py").is_file():
            return candidate.resolve()
    raise RuntimeError("set H10Q_ROOT to the math/h10q directory")


ROOT = find_root()
sys.path.insert(0, str(ROOT))

import h10q  # noqa: E402
import l12_class  # noqa: E402
import l13_filter  # noqa: E402


CELL = (21, (-1, 1))
W, UT = CELL
Z = F(W) * F(*UT)
A = F(1)
TAU = F(3, 5)
F_ESC = 3
EPS = 1
Q1 = 11
OUT = ROOT / "data" / "l17_composite_w.jsonl"
REPORT = Path(os.environ.get("COMPOSITE_W_REPORT", "/tmp/l17_composite_w.report"))

SOURCE_MAP = {
    "architecture_place": [
        "THEOREMS.md:193-215 (w is an odd place and m_w is its maximal ideal)",
        "THEOREMS.md:246-254 (target residue field, v_w, Hensel lift, ramification at S union {w})",
        "THEOREMS.md:504-525 (W1 is explicitly w-adic at an odd place)",
    ],
    "model_field": [
        "h10q.py:854-875 (_GFq accepts only an odd prime power)",
        "h10q.py:962-969 (_field uses F_q: prime q or a prime-power extension)",
        "h10q.py:2540-2551 (_l9_grid draws w from primerange)",
    ],
    "character_and_valuation": [
        "h10q.py:333-355 (v_p/unit/Legendre implementation)",
        "h10q.py:359-376 (Hilbert symbol documents v as a prime)",
        "h10q.py:847-850 (_sun_target_tau uses w mod 4, valid for the prime target)",
        "THEOREMS.md:477-501 (W0 works over a finite field and then specializes by prime w mod 4)",
    ],
    "l11_root": [
        "THEOREMS.md:1241-1258 (escape is through a numerator prime p=1 mod 4 and a root modulo p)",
        "h10q.py:3128-3130 (_l11_escape_primes factors numerator(z) and returns prime factors 1 mod 4)",
    ],
    "l10_class": [
        "h10q.py:2617-2649 (certificate depends on a,z,tau,eps,q1 and freezes prime support)",
        "h10q.py:2650-2688 (q1, not w, must be prime; Hilbert symbols are evaluated prime by prime)",
    ],
    "l12_escape": [
        "THEOREMS.md:1325-1332 (escape requires a shared prime f dividing z)",
        "h10q.py:3390-3415 (conservative symbols factor all rational data into prime places)",
        "h10q.py:3418-3420 (escape rows use f dividing z)",
        "l12_class.py:70-115 (computational certificate requires both q and f prime and v_f(z)>0)",
    ],
    "l13e_composite_A": [
        "THEOREMS.md:1415-1429 (composite squarefree A wall is factor-by-factor at p dividing A)",
        "h10q.py:4262-4308 (replay factors A and multiplies symbols over its prime support)",
    ],
}

CLOSURE_KEYS = {
    "family", "cell", "a", "eps", "f", "q1", "N", "k_zero", "Q",
    "rung", "tied", "ramified_empty", "source",
}


def emit(handle, record) -> None:
    handle.write(json.dumps(record) + "\n")
    handle.flush()


def refusal(exc: BaseException) -> str:
    return "refused:" + type(exc).__name__


def temporary_l12_cert(f: int, eps: int, q: int):
    """Call the table-backed certificate with no persistent table mutation."""
    sentinel = object()
    previous = l12_class._L12_ESCAPE.get(CELL, sentinel)
    l12_class._L12_ESCAPE[CELL] = (f, eps, q)
    try:
        return l12_class.l12_class_cert(W, *UT)
    finally:
        if previous is sentinel:
            l12_class._L12_ESCAPE.pop(CELL, None)
        else:
            l12_class._L12_ESCAPE[CELL] = previous


def model_field_status(w: int) -> tuple[str, str | None]:
    try:
        h10q._field(w)
    except Exception as exc:
        return "rejected", f"{type(exc).__name__}:{exc}"
    return "accepted", None


def candidate_records() -> list[dict]:
    """All odd squarefree composite w=1 mod 4 below 100."""
    records = []
    for w in range(3, 100, 2):
        factors = h10q.factorint(w)
        if w % 4 != 1 or len(factors) < 2 or any(e != 1 for e in factors.values()):
            continue
        time.sleep(0.1)
        status, detail = model_field_status(w)
        roots = [r for r in range(1, w) if (1 + 4 * r * r) % w == 0]
        z = F(-w)
        component_taus = {str(p): str(h10q._sun_target_tau(A, p)) for p in factors}
        composite_tau = str(h10q._sun_target_tau(A, w))
        records.append({
            "type": "candidate",
            "label": "PROVED",
            "w": w,
            "factorization": {str(p): e for p, e in sorted(factors.items())},
            "squarefree": True,
            "w_mod_4": w % 4,
            "prime_place": h10q._is_prime(w),
            "architecture_admissible_as_one_place": False,
            "model_field_status": status,
            "model_field_detail": detail,
            "l11_roots_mod_w": roots,
            "l11_escape_primes_from_z": h10q._l11_escape_primes(z),
            "composite_mod4_tau_shortcut": composite_tau,
            "component_prime_tau_shortcuts": component_taus,
            "composite_shortcut_matches_components": all(
                tau == composite_tau for tau in component_taus.values()
            ),
        })
    assert [record["w"] for record in records] == [21, 33, 57, 65, 69, 77, 85, 93]
    return records


def replay_closure(cert: dict, cofactor_info: dict) -> tuple[dict, dict]:
    """Recompute every decisive kernel, class, and L13 fact before emission."""
    k = 0
    Q = Q1 + k * int(cert["N"])
    b = F(EPS * F_ESC * Q)
    assert Q == Q1 == 11
    assert h10q._is_prime(Q)
    assert h10q._is_prime(F_ESC)
    assert h10q._is_prime(W) is False
    assert h10q.factorint(W) == {3: 1, 7: 1}
    assert h10q.vp(Z, F_ESC) == 1

    a0 = A
    tau = TAU
    A0 = 1 + 4 * a0 * a0
    delta = 1 - A0 * tau * tau
    alpha = -delta * A0
    Z3 = Z ** 3
    D = 1 - Z3 - a0 * a0 * Z3 * Z3
    s = (a0 - 1) / 2
    c0 = h10q._sun_h(a0, b, Z3)
    assert c0 is not None
    M0 = 16 - delta * c0 * c0 - 32 * A0 * b * s * s
    x0, d0 = alpha * M0, alpha * 2 * b
    P = h10q._l10_P(a0, Z3, D, A0, delta, s)
    Pb = sum(coefficient * b ** i for i, coefficient in enumerate(P))
    assert M0 == Pb / (b ** 4 * D * D * A0 * A0)

    class_replay = temporary_l12_cert(F_ESC, EPS, Q1)
    assert class_replay == cert and cert["ok"]
    assert all(value == 1 for value in cert["syms"].values())
    conservative_syms = h10q._l12_syms(a0, Z, tau, b)
    aligned_syms = l13_filter.aligned(a0, Z, tau, b)
    assert conservative_syms is not None and all(value == 1 for value in conservative_syms.values())
    assert aligned_syms is not None and all(value == 1 for value in aligned_syms.values())

    smooth_status, smooth_detail = l13_filter.smooth_emergent(a0, Z, tau, b)
    assert smooth_status == "cofactor-big"
    emergent, (rn, rd), smooth = smooth_detail
    assert emergent == [] and not smooth and rd == 1
    verdict, replay_info = l13_filter.cofactor_decide(
        a0, Z, tau, b, detail=smooth_detail, deadline=time.time() + 600
    )
    assert verdict == "zero"
    assert replay_info == cofactor_info
    assert replay_info["rung"] == "factorint"
    assert replay_info["places"]
    assert all(place[1] == 1 and place[2] == "proved" for place in replay_info["places"])
    for place in replay_info["places"]:
        assert h10q._is_prime(int(place[0]))

    tied = h10q._l7_tied_status(a0, b, Z, tau)
    ramified = h10q.ramified(x0, d0)
    assert tied is True and ramified == []

    Pi = h10q._l13_make_P(Z, "sq")
    irreducibility_prime = None
    for p in h10q.primerange(7, 1000):
        if p == 5 or Pi[-1] % p == 0:
            continue
        time.sleep(0.1)
        if h10q._l13_irred8(Pi, p):
            irreducibility_prime = p
            break
    assert len(Pi) - 1 == 8 and irreducibility_prime == 11

    closure = {
        "family": "ESC",
        "cell": [W, list(UT)],
        "a": str(a0),
        "eps": EPS,
        "f": F_ESC,
        "q1": Q1,
        "N": str(cert["N"]),
        "k_zero": k,
        "Q": str(Q),
        "rung": replay_info["rung"],
        "tied": str(tied),
        "ramified_empty": not ramified,
        "source": "l17_composite_w.py",
    }
    assert len(closure) == 13 and set(closure) == CLOSURE_KEYS

    replay = {
        "type": "kernel-replay",
        "label": "PROVED",
        "interpretation": "arithmetic certificate for z=-21 after selecting the actual prime place f=3; not a composite-place certificate",
        "cell": [W, list(UT)],
        "z": str(Z),
        "w_factorization": {"3": 1, "7": 1},
        "selected_target_prime": F_ESC,
        "a": str(a0),
        "b": str(b),
        "tau": str(tau),
        "A": str(A0),
        "delta": str(delta),
        "alpha": str(alpha),
        "D": str(D),
        "c0": str(c0),
        "M0": str(M0),
        "x0": str(x0),
        "d0": str(d0),
        "kernel_identity": True,
        "class_cert_replayed": True,
        "conservative_syms": conservative_syms,
        "aligned_syms": aligned_syms,
        "smooth_status": smooth_status,
        "small_emergent": emergent,
        "cofactor_rn": str(rn),
        "cofactor_rd": str(rd),
        "cofactor_verdict": verdict,
        "cofactor_info": replay_info,
        "tied": tied,
        "ramified_empty": not ramified,
        "P_degree": len(Pi) - 1,
        "P_irreducibility_prime": irreducibility_prime,
        "P_irreducible_over_Q": True,
        "P_coefficients": Pi,
    }
    return closure, replay


def main() -> None:
    h10q._selftest()
    assert Z == F(-21)
    assert h10q.factorint(W) == {3: 1, 7: 1}
    assert h10q._is_prime(W) is False
    assert CELL not in h10q._L11_CLASSES

    candidates = candidate_records()
    candidate21 = next(record for record in candidates if record["w"] == W)
    assert candidate21["model_field_status"] == "rejected"
    assert candidate21["l11_roots_mod_w"] == []
    assert candidate21["l11_escape_primes_from_z"] == []
    assert candidate21["composite_shortcut_matches_components"] is False

    # Literal template substitution f=w fails first: f is required to be prime.
    time.sleep(0.1)
    composite_f_cert = temporary_l12_cert(W, EPS, Q1)
    assert composite_f_cert is None

    # Prime-place decomposition is the meaningful repair: choose f=3 | z.
    time.sleep(0.1)
    cert = temporary_l12_cert(F_ESC, EPS, Q1)
    assert cert is not None and cert["ok"]
    assert cert["N"] == 272160
    assert cert["S"] == [2, 3, 5, 7]
    assert all(value == 1 for value in cert["syms"].values())

    Q = Q1
    b = F(EPS * F_ESC * Q)
    time.sleep(0.1)
    smooth_status, smooth_detail = l13_filter.smooth_emergent(A, Z, TAU, b)
    assert smooth_status == "cofactor-big"
    time.sleep(0.1)
    verdict, cofactor_info = l13_filter.cofactor_decide(
        A, Z, TAU, b, detail=smooth_detail, deadline=time.time() + 600
    )
    assert verdict == "zero" and cofactor_info["rung"] == "factorint"

    closure, replay = replay_closure(cert, cofactor_info)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8") as out:
        emit(out, {
            "type": "meta",
            "schema": "l17-composite-w-v1",
            "label": "PROVED",
            "engine": "math/h10q/h10q.py proven-primality only",
            "selftest": True,
            "cell": [W, list(UT)],
            "z": str(Z),
            "verdict": "no composite place; arithmetic closure succeeds after prime-factor decomposition",
            "refusals_are_not_evidence": True,
            "closure_row_fields": 13,
        })
        emit(out, {"type": "source-map", "label": "PROVED", "sources": SOURCE_MAP})
        for record in candidates:
            emit(out, record)
        emit(out, {
            "type": "attempt-step",
            "label": "PROVED",
            "step": "model-field",
            "cell": [W, list(UT)],
            "status": "rejected:not-an-odd-prime-power",
            "detail": candidate21["model_field_detail"],
            "first_structural_break": True,
        })
        emit(out, {
            "type": "attempt-step",
            "label": "PROVED",
            "step": "L11-root",
            "cell": [W, list(UT)],
            "status": "no-root-mod-composite-w",
            "roots": candidate21["l11_roots_mod_w"],
            "escape_primes_1_mod_4": candidate21["l11_escape_primes_from_z"],
        })
        emit(out, {
            "type": "attempt-step",
            "label": "PROVED",
            "step": "L12-literal-f-equals-w",
            "cell": [W, list(UT)],
            "f": W,
            "eps": EPS,
            "q1": Q1,
            "status": "rejected:composite-f",
            "cert": None,
        })
        emit(out, {
            "type": "attempt-step",
            "label": "PROVED",
            "step": "L12-prime-factor-escape",
            "cell": [W, list(UT)],
            "f": F_ESC,
            "eps": EPS,
            "q1": Q1,
            "status": "aligned",
            "N": str(cert["N"]),
            "S": cert["S"],
            "ks": cert["ks"],
            "syms": cert["syms"],
            "Q0": str(cert["Q0"]),
            "excluded": cert["excluded"],
        })
        emit(out, {
            "type": "member",
            "label": "PROVED",
            "family": "ESC-prime-factor",
            "cell": [W, list(UT)],
            "a": str(A),
            "eps": EPS,
            "f": F_ESC,
            "q1": Q1,
            "N": str(cert["N"]),
            "k": 0,
            "Q": str(Q),
            "Q_prime_proven": True,
            "b": str(b),
            "smooth_status": smooth_status,
            "verdict": verdict,
            "rung": cofactor_info["rung"],
            "cofactor_info": cofactor_info,
        })
        emit(out, replay)
        emit(out, closure)
        emit(out, {
            "type": "summary",
            "label": "PROVED",
            "status": "ARITHMETIC-CLOSED_STRUCTURAL-NOT-A-COMPOSITE-PLACE",
            "cell": [W, list(UT)],
            "selected_prime_place": F_ESC,
            "Q": str(Q),
            "closure_schema_fields": len(closure),
        })

    report_lines = [
        "COMPOSITE-w STRUCTURAL + COMPUTATIONAL INVESTIGATION",
        "",
        "STRUCTURAL VERDICT [PROVED from source]: A squarefree composite integer w is not an odd place of Q, so the architecture does not admit it as one target place. The first break is the residue-field/model step: _field(w) accepts a prime or prime power, while the grid itself enumerates w with primerange. For composite z, the correct operation is to factor its numerator and choose an actual prime place p dividing z.",
        "",
        "ARITHMETIC VERDICT [PROVED by exact replay]: The downstream rational closure machinery does handle composite support. At the bookkeeping pseudo-cell [21,[-1,1]], z=-21, the literal template choice f=w=21 is rejected because L12 requires f prime. Choosing the actual prime factor f=3 gives an aligned class at a=1, tau=3/5, eps=1, q1=11, N=272160. Its k=0 member Q=11, b=33 has smooth_emergent=cofactor-big and cofactor_decide=zero at rung factorint; the kernel identity, all conservative and aligned symbols, tied status True, empty ramification, and degree-eight irreducibility certificate mod 11 are replayed in-script. The emitted closure object has exactly 13 fields.",
        "",
        "INTERPRETATION [PROVED]: This is a certificate for the rational input z=-21 after selecting the genuine prime place 3. It is not evidence for a new composite place. Since z is also divisible by 7, it lies at multiple ordinary prime places; the architecture's union over places already handles that situation without defining a composite valuation.",
        "",
        "CANDIDATE CHECK [PROVED]: The odd squarefree composite w=1 mod 4 below 100 are 21, 33, 57, 65, 69, 77, 85, 93 (replayed in data/l17_composite_w.jsonl). None is admissible as one place or as one finite residue field. The suggested values 21, 33, 57, 69, 77, 93 have no root of 1+4r^2 modulo w and no numerator prime factor 1 mod 4. Values 65 and 85 do have roots modulo w because all their prime factors are 1 mod 4, but they still are products of places, not fields/places themselves.",
        "",
        "PRIMALITY ASSUMPTIONS AND NON-ASSUMPTIONS [PROVED from source]:",
        "1. Place/model: THEOREMS.md:193-215 and 246-254 formulate m_w, kappa_w, v_w, Hensel lifting, and ramification for an odd place w. THEOREMS.md:504-525 repeats this in W1. h10q.py:854-875 and 962-969 implement only finite fields of prime-power order; h10q.py:2540-2551 constructs grid labels from primerange. h10q.py:333-376 implements v_p, Legendre, and Hilbert at prime p; passing composite w to those routines would not define the required local objects.",
        "2. W0 branch shortcut: THEOREMS.md:477-501 proves the selector over a finite field and only then writes the rational specialization by prime w mod 4; h10q.py:847-850 mirrors it. For w=21 the composite congruence shortcut selects the w=1 mod 4 branch, while both prime components 3 and 7 select the w=3 mod 4 branch. Thus w mod 4 is not a substitute for componentwise residue fields.",
        "3. L11 root: THEOREMS.md:1241-1258 requires a numerator prime p=1 mod 4 and solves modulo that p. h10q.py:3128-3130 factors numerator(z) and returns those prime factors. The horizon templates may solve modulo W only because their W is prime. For squarefree composite W, a root modulo W exists exactly when it exists at every prime factor; that stronger CRT condition is unnecessary for L11, which needs one suitable prime factor of z.",
        "4. _l10_class_cert: h10q.py:2617-2688 has no w argument and makes no primality assumption on z's numerator. It factors the supports of alpha and delta, freezes actual primes, and requires q1 prime at lines 2660-2663 because L9a and Legendre(A,q1) need it. Therefore its symbol bookkeeping already accepts composite support in z; the prime assumption belongs to q1 and each Hilbert place, not to an integer label w.",
        "5. L12 escape: THEOREMS.md:1325-1332 and h10q.py:3418-3420 require a shared prime f dividing z. h10q.py:3390-3415 factors all fixed rational data into prime places. l12_class.py:70-115 enforces q and f prime plus v_f(z)>0. Hence f=w is valid only in the prime-horizon templates; for composite w one must choose an actual prime factor f of z.",
        "6. L13e: THEOREMS.md:1415-1429 and h10q.py:4262-4308 concern composite squarefree A=1+4a^2, not composite w. The wall is proved by factoring A and multiplying local symbols factor by factor. It already removes primality of A from the wall, but it neither creates a composite residue field nor licenses a composite place label.",
        "",
        "WHAT WOULD BE NEEDED [STRUCTURAL CONSEQUENCE]: Keep target keys prime-place based: factor a composite multiplier W, choose p|W, and carry the full rational z separately. Run W0/W1, L11 roots, valuations, characters, and Hensel lifting separately in each F_p/Q_p; for an L12 escape choose a prime f|numerator(z). If a product-modulus API were nevertheless desired, it would need an explicit product of residue fields with componentwise characters and local symbols—never the current Legendre/Hilbert routines at composite w. _l10_class_cert and the L13 cofactor ladder need no conceptual change once supplied rational z and genuine prime places.",
        "",
        "SOURCE MAP:",
    ]
    for topic, sources in SOURCE_MAP.items():
        report_lines.append(f"- {topic}: " + "; ".join(sources))
    report_lines.extend([
        "",
        f"ARTIFACT: {OUT}",
        "REPLAY: python3 math/h10q/l17_composite_w.py (from workspace, or set H10Q_ROOT)",
        "ENGINE POLICY: h10q.py proven-primality only; refusals are never evidence.",
    ])
    REPORT.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(json.dumps(closure, sort_keys=True))
    print("WROTE", OUT)
    print("WROTE", REPORT)


if __name__ == "__main__":
    main()
