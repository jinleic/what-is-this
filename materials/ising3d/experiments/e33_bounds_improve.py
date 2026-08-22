"""Audit and improve the certified simple-cubic Ising critical interval.

The improved lower endpoint uses the published exact 36-step simple-cubic
self-avoiding-walk count and submultiplicativity.  All transcendental endpoints
are independently enclosed with ``mpmath.iv`` directed rounding.  The finite-volume
Simon--Lieb and Peierls computations below are diagnostics, never promoted to
thermodynamic theorems.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from pathlib import Path

import mpmath as mp
import numpy as np


SCRIPT = "experiments/e33_bounds_improve.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "bounds" / "improved_bounds.json"
INCUMBENT_PATH = ROOT / "results" / "bounds" / "kc_bounds.json"
SIMON_LIEB_PATH = ROOT / "results" / "bounds" / "simon_lieb_bounds.json"
SAW_SOURCE_PATH = ROOT / "sources" / "fulltext" / "schram_barkema_bisseling2011.pdf"
SAW_SOURCE_SHA256 = "898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12"
PEIERLS_SOURCE_PATH = ROOT / "sources" / "fulltext" / "bonati2014_peierls.pdf"
PEIERLS_SOURCE_SHA256 = "2bd9daa39aad312c9a1141b404a9de19020c85e97d0b779644ddfe71567acd36"
DPS = 80
DECIMAL_PLACES = 40
BENCHMARK = Decimal("0.221654626")  # comparison only
C36 = 2941370856334701726560670
C36_STEPS = 36
MACDONALD_REPORTED_MU_UPPER = (47114, 10000)


def _interval_endpoints(value: object) -> tuple[str, str]:
    text = str(value).strip()
    if not (text.startswith("[") and text.endswith("]")):
        raise ValueError(f"unexpected interval representation: {text}")
    lower, upper = text[1:-1].split(",", maxsplit=1)
    return lower.strip(), upper.strip()


def _floor_decimal(value: str, places: int = DECIMAL_PLACES) -> str:
    with localcontext() as context:
        context.prec = max(DPS + 20, len(value) + 10)
        quantum = Decimal(1).scaleb(-places)
        return format(Decimal(value).quantize(quantum, rounding=ROUND_FLOOR), "f")


def _ceil_decimal(value: str, places: int = DECIMAL_PLACES) -> str:
    with localcontext() as context:
        context.prec = max(DPS + 20, len(value) + 10)
        quantum = Decimal(1).scaleb(-places)
        return format(Decimal(value).quantize(quantum, rounding=ROUND_CEILING), "f")


def _atanh_inverse_rational_interval(
    numerator: int, denominator: int, *, dps: int = DPS
) -> tuple[str, str]:
    """Enclose atanh(denominator/numerator) by directed interval arithmetic."""

    previous = mp.iv.dps
    try:
        mp.iv.dps = dps
        mu = mp.iv.mpf(numerator) / denominator
        v = 1 / mu
        return _interval_endpoints(mp.iv.log((1 + v) / (1 - v)) / 2)
    finally:
        mp.iv.dps = previous


def _count_root_interval(
    count: int, steps: int, *, dps: int = DPS
) -> tuple[str, str]:
    """Enclose ``count**(1/steps)`` by directed interval arithmetic."""

    previous = mp.iv.dps
    try:
        mp.iv.dps = dps
        return _interval_endpoints(
            mp.iv.exp(mp.iv.log(mp.iv.mpf(count)) / steps)
        )
    finally:
        mp.iv.dps = previous


def _atanh_inverse_count_root_interval(
    count: int, steps: int, *, digits: int = 90, dps: int = DPS
) -> tuple[str, str, str, str]:
    """Certify ``atanh(count**(-1/steps))`` from an exact rational bracket."""

    scale = 10**digits
    target = scale**steps
    lo, hi = 0, scale
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if count * mid**steps < target:
            lo = mid
        else:
            hi = mid
    if not (count * lo**steps < target <= count * hi**steps):
        raise AssertionError("failed to bracket inverse count root exactly")

    previous = mp.iv.dps
    try:
        mp.iv.dps = dps
        v_lo = mp.iv.mpf(lo) / scale
        v_hi = mp.iv.mpf(hi) / scale
        k_lo = mp.iv.log((1 + v_lo) / (1 - v_lo)) / 2
        k_hi = mp.iv.log((1 + v_hi) / (1 - v_hi)) / 2
        lower = _interval_endpoints(k_lo)[0]
        upper = _interval_endpoints(k_hi)[1]
    finally:
        mp.iv.dps = previous
    return f"{lo}/{scale}", f"{hi}/{scale}", lower, upper


def _peierls_root_interval(*, dps: int = DPS) -> tuple[str, str, str, str]:
    """Bracket the Bonati 3D Peierls equality using interval-safe rational bisection.

    The theorem bounds the minus-spin density by
    ``F(x)=3*x^3*(9-11*x+4*x^2)/(16*(1-x)^3)``, with
    ``x=9*exp(-4K)``.  We bracket the unique root F(x)=1/2 using exact
    rational arithmetic for the sign, then transform the rational bracket with
    directed interval logarithms.  This candidate is intentionally weak.
    """

    scale = 10**60

    def sign_at(p: int) -> int:
        # sign of F(p/scale)-1/2; all denominators are positive on (0,1).
        q = scale
        return 3 * p**3 * (9 * q * q - 11 * p * q + 4 * p * p) - 8 * (
            q - p
        ) ** 3 * q * q

    lo = 0
    hi = scale
    for _ in range(220):
        mid = (lo + hi) // 2
        if sign_at(mid) < 0:
            lo = mid
        else:
            hi = mid
        if hi - lo <= 1:
            break
    if not (sign_at(lo) < 0 <= sign_at(hi)):
        raise AssertionError("failed to bracket Peierls root exactly")

    previous = mp.iv.dps
    try:
        mp.iv.dps = dps
        # K(x)=-log(x/9)/4 decreases with x.  Evaluate both rational endpoints.
        k_from_hi = -mp.iv.log((mp.iv.mpf(hi) / scale) / 9) / 4
        k_from_lo = -mp.iv.log((mp.iv.mpf(lo) / scale) / 9) / 4
        lower = _interval_endpoints(k_from_hi)[0]
        upper = _interval_endpoints(k_from_lo)[1]
    finally:
        mp.iv.dps = previous
    return f"{lo}/{scale}", f"{hi}/{scale}", lower, upper


def _simon_lieb_orbit_diagnostic(
    *,
    shape: tuple[int, int, int] = (4, 4, 32),
    origin: tuple[int, int, int] = (1, 1, 15),
    v: float = 0.20234669789378,
) -> dict[str, object]:
    """Recompute the boundary split with an independent spin-basis transfer."""


    nx, ny, nz = shape
    layer_sites = nx * ny
    state_count = 1 << layer_sites
    indices = np.arange(state_count, dtype=np.uint32)
    coupling = math.atanh(v)
    layer_energy = np.zeros(state_count, dtype=np.int16)
    for y in range(ny):
        for x in range(nx):
            site = x + nx * y
            spin = 2 * ((indices >> site) & 1).astype(np.int16) - 1
            if x + 1 < nx:
                neighbor = 2 * ((indices >> (site + 1)) & 1).astype(np.int16) - 1
                layer_energy += spin * neighbor
            if y + 1 < ny:
                neighbor = 2 * ((indices >> (site + nx)) & 1).astype(np.int16) - 1
                layer_energy += spin * neighbor

    half_weight = np.exp(coupling * layer_energy / 2)
    same = math.exp(coupling)
    different = math.exp(-coupling)

    def transfer(vector: np.ndarray) -> tuple[np.ndarray, float]:
        work = half_weight * vector
        for bit in range(layer_sites):
            stride = 1 << bit
            view = work.reshape(-1, 2 * stride)
            low = view[:, :stride].copy()
            high = view[:, stride:].copy()
            view[:, :stride] = same * low + different * high
            view[:, stride:] = different * low + same * high
        result = half_weight * work
        scale = float(np.max(np.abs(result)))
        if scale == 0:
            raise AssertionError("zero transfer vector")
        return result / scale, math.log(scale)

    left = []
    left_logs = []
    current = half_weight / half_weight.max()
    log_scale = math.log(float(half_weight.max()))
    left.append(current)
    left_logs.append(log_scale)
    for _ in range(1, nz):
        current, increment = transfer(current)
        log_scale += increment
        left.append(current)
        left_logs.append(log_scale)

    right = [None] * nz
    right_logs = [0.0] * nz
    current = half_weight / half_weight.max()
    log_scale = math.log(float(half_weight.max()))
    right[-1] = current
    right_logs[-1] = log_scale
    for z in range(nz - 2, -1, -1):
        current, increment = transfer(current)
        log_scale += increment
        right[z] = current
        right_logs[z] = log_scale

    spins = np.empty((layer_sites, state_count), dtype=np.int8)
    for bit in range(layer_sites):
        spins[bit] = 2 * ((indices >> bit) & 1).astype(np.int8) - 1
    ox, oy, oz = origin
    origin_bit = ox + nx * oy
    correlations = np.empty((nz, layer_sites), dtype=np.float64)

    marked = left[oz] * spins[origin_bit]
    marked_log = left_logs[oz]
    for z in range(oz, nz):
        if z > oz:
            marked, increment = transfer(marked)
            marked_log += increment
        correlations[z] = (
            math.exp(marked_log - left_logs[z])
            * (spins @ (marked * right[z]))
            / np.dot(left[z], right[z])
        )

    marked = right[oz] * spins[origin_bit]
    marked_log = right_logs[oz]
    for z in range(oz - 1, -1, -1):
        marked, increment = transfer(marked)
        marked_log += increment
        correlations[z] = (
            math.exp(marked_log - right_logs[z])
            * (spins @ (left[z] * marked))
            / np.dot(left[z], right[z])
        )

    orbit_data: dict[tuple[int, int, int], dict[str, object]] = {}
    total = 0.0
    lateral = 0.0
    caps = 0.0
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                correlation = float(correlations[z, x + nx * y])
                lateral_q = int(x == 0) + int(x + 1 == nx) + int(y == 0) + int(
                    y + 1 == ny
                )
                cap_q = int(z == 0) + int(z + 1 == nz)
                contribution = v * (lateral_q + cap_q) * correlation
                total += contribution
                lateral += v * lateral_q * correlation
                caps += v * cap_q * correlation
                if contribution:
                    key = (min(x, y), max(x, y), z)
                    if key not in orbit_data:
                        orbit_data[key] = {
                            "representative_xyz": list(key),
                            "multiplicity": 0,
                            "q": lateral_q + cap_q,
                            "correlation": correlation,
                            "contribution": 0.0,
                        }
                    orbit_data[key]["multiplicity"] += 1  # type: ignore[operator]
                    orbit_data[key]["contribution"] += contribution  # type: ignore[operator]

    ordered = sorted(
        orbit_data.values(),
        key=lambda item: float(item["contribution"]),
        reverse=True,
    )
    largest = ordered[0]
    return {
        "method": "independent binary64 spin-basis transfer; recomputed on every experiment run; descriptive only",
        "v": format(v, ".14f"),
        "kappa": repr(total),
        "lateral_faces_contribution": repr(lateral),
        "end_caps_contribution": repr(caps),
        "lateral_fraction": repr(lateral / total),
        "largest_rooted_orbit": {
            "representative_xyz": largest["representative_xyz"],
            "multiplicity": largest["multiplicity"],
            "q": largest["q"],
            "correlation": repr(float(largest["correlation"])),
            "kappa_contribution": repr(float(largest["contribution"])),
        },
        "ten_largest_orbits_fraction": repr(
            sum(float(item["contribution"]) for item in ordered[:10]) / total
        ),
        "interpretation": "near-origin lateral boundary orbits dominate; caps are negligible, hence lengthening 4x4 prisms has saturated while the 4x4 cross-section remains the bottleneck",
    }


def _simon_lieb_diagnostic(payload: dict[str, object]) -> dict[str, object]:
    boxes = payload["data"]["boxes"]  # type: ignore[index]
    selected = {tuple(box["shape"]): box for box in boxes}  # type: ignore[index]
    roots = [
        Decimal(selected[(4, 4, length)]["numerical_root"]["v"])  # type: ignore[index]
        for length in (8, 16, 32)
    ]
    # A three-point exponential fit is descriptive metadata only.  It is not a bound.
    with mp.workdps(DPS):
        r8, r16, r32 = (mp.mpf(str(value)) for value in roots)
        ratio = (r32 - r16) / (r16 - r8)
        rho = mp.root(ratio, 8)
        amplitude = (r16 - r8) / (rho**16 - rho**8)
        limit = r8 - amplitude * rho**8
        limit_k = mp.atanh(limit)
    orbit_diagnostic = _simon_lieb_orbit_diagnostic()

    return {
        "strongest_box": {
            "shape": [4, 4, 32],
            "origin": [1, 1, 15],
            "certified_K_lower": selected[(4, 4, 32)]["certified_endpoint"][  # type: ignore[index]
                "reported_decimal_lower"
            ],
        },
        "exact_boundary_decomposition": {
            "identity": "kappa=v*sum_x q_B(x)<sigma_o sigma_x>_B; q_B is the exact number of crossing edges at x",
            "symmetry_grouping": "the rooted box stabilizer only permutes equal positive summands, so grouping cannot alter kappa or its root",
            "separating_surface_test": "for the theorem as implemented every boundary edge must be included; deleting or down-weighting a positive crossing term is not licensed by Simon--Lieb",
            "provable_tail_test": "positive end-cap or axial-tail terms may be upper-bounded, but the exact transfer already includes them; replacement by an upper bound cannot lower kappa",
        },
        "boundary_orbit_diagnostic_at_certified_v": orbit_diagnostic,
        "descriptive_length_limit_fit": {
            "model": "v_L = v_inf + a*rho^L fitted to L=8,16,32; numerical diagnostic only",
            "v_inf": mp.nstr(limit, 70),
            "rho": mp.nstr(rho, 70),
            "K_inf": mp.nstr(limit_k, 70),
            "rigor": "not a certificate and not used in any endpoint",
        },
    }


def main() -> None:
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})
        if not passed:
            raise AssertionError(f"{name}: {detail}")
        print(f"PASS {name}: {detail}")

    incumbent = json.loads(INCUMBENT_PATH.read_text(encoding="utf-8"))
    simon_lieb = json.loads(SIMON_LIEB_PATH.read_text(encoding="utf-8"))
    old_lower, old_upper = incumbent["data"]["final_certified_interval"][
        "decimal_outward_rounded"
    ]

    source_hash = hashlib.sha256(SAW_SOURCE_PATH.read_bytes()).hexdigest()
    v_root_lo, v_root_hi, saw_lo, saw_hi = _atanh_inverse_count_root_interval(
        C36, C36_STEPS
    )
    macdonald_lo, macdonald_hi = _atanh_inverse_rational_interval(
        *MACDONALD_REPORTED_MU_UPPER
    )
    mu_lo, mu_hi = _count_root_interval(C36, C36_STEPS)
    new_lower = _floor_decimal(saw_lo)
    p_root_lo, p_root_hi, peierls_lo, peierls_hi = _peierls_root_interval()
    peierls_upper = _ceil_decimal(peierls_hi)
    new_upper = str(old_upper)
    old_width = Decimal(old_upper) - Decimal(old_lower)
    new_width = Decimal(new_upper) - Decimal(new_lower)
    lower_improvement = Decimal(new_lower) - Decimal(old_lower)
    simon_diagnostic = _simon_lieb_diagnostic(simon_lieb)
    orbit = simon_diagnostic["boundary_orbit_diagnostic_at_certified_v"]
    orbit_total = Decimal(str(orbit["kappa"]))  # type: ignore[index]
    orbit_lateral = Decimal(str(orbit["lateral_faces_contribution"]))  # type: ignore[index]
    orbit_caps = Decimal(str(orbit["end_caps_contribution"]))  # type: ignore[index]

    check(
        "exact_saw_source_sha256",
        source_hash == SAW_SOURCE_SHA256,
        f"cached full text sha256={source_hash}",
    )
    peierls_source_hash = hashlib.sha256(PEIERLS_SOURCE_PATH.read_bytes()).hexdigest()
    check(
        "peierls_source_sha256",
        peierls_source_hash == PEIERLS_SOURCE_SHA256,
        f"cached full text sha256={peierls_source_hash}",
    )
    check(
        "published_c36_improves_lower_endpoint",
        Decimal(new_lower) > Decimal(old_lower),
        f"{new_lower} > {old_lower}",
    )
    check(
        "directed_rounding_order",
        Decimal(new_lower) <= Decimal(saw_lo) <= Decimal(saw_hi),
        "reported lower endpoint is rounded downward from the 80-dps interval",
    )
    check(
        "peierls_candidate_is_weaker",
        Decimal(peierls_upper) > Decimal(old_upper),
        f"Peierls gives only Kc <= {peierls_upper}",
    )
    check(
        "simon_lieb_boundary_split",
        abs(orbit_total - orbit_lateral - orbit_caps) < Decimal("2e-15")
        and orbit_caps / orbit_total < Decimal("0.000004"),
        "recomputed split has >99.9996% lateral mass and <0.0004% cap mass",
    )
    check(
        "benchmark_inside_final_interval",
        Decimal(new_lower) < BENCHMARK < Decimal(new_upper),
        "comparison-only benchmark lies strictly inside the certified interval",
    )

    data = {
        "model": "nearest-neighbor ferromagnetic Ising model on the simple cubic lattice",
        "coupling_convention": "K=beta*J; v=tanh(K)",
        "incumbent_interval": {
            "source": str(INCUMBENT_PATH.relative_to(ROOT)),
            "decimal_outward_rounded": [old_lower, old_upper],
        },
        "candidates": [
            {
                "id": "published_exact_saw_c36",
                "side": "lower",
                "classification": "[EXTERNAL COMPUTATION] plus [THEOREM] plus [COMPUTATION] directed rounding",
                "theorem": "block submultiplicativity c_(36q+r)<=c_36^q*c_r directly makes the SAW susceptibility majorant geometric for v<c_36^(-1/36)",
                "primary_source": {
                    "key": "schram2011_exact_saw",
                    "citation": "R. D. Schram, G. T. Barkema, and R. H. Bisseling, Exact enumeration of self-avoiding walks (2011)",
                    "doi": "10.1088/1742-5468/2011/06/P06019",
                    "url": "https://webspace.science.uu.nl/~bisse101/Articles/schram11.pdf",
                    "exact_location": "Section IV, first paragraph, and Table I: Z_36=2941370856334701726560670",
                    "local_path": str(SAW_SOURCE_PATH.relative_to(ROOT)),
                    "sha256": source_hash,
                },
                "certificate": {
                    "published_exact_count": {"steps": C36_STEPS, "c_n": C36},
                    "mu_upper_exact": "2941370856334701726560670^(1/36)",
                    "mu_upper_interval_directed_rounding": [mu_lo, mu_hi],
                    "logical_chain": [
                        "[EXTERNAL COMPUTATION] c_0,...,c_36, including c_36=2941370856334701726560670, by the exact length-doubling inclusion-exclusion algorithm",
                        "[THEOREM] for n=36q+r, repeated submultiplicativity gives c_n<=c_36^q*c_r",
                        "[THEOREM] chi(v)<=sum_n c_n*v^n<=(sum_(r=0)^35 c_r*v^r)/(1-c_36*v^36)",
                        "therefore chi is finite and v<v_c whenever v<c_36^(-1/36)",
                        "[THEOREM] atanh is strictly increasing on (0,1)",
                        "therefore K_c>=atanh(c_36^(-1/36))",
                        "[REMARK] Fekete equivalently gives mu<=c_36^(1/36), but the direct certificate does not need the limit defining mu",
                    ],
                    "v_exact_rational_bracket": [v_root_lo, v_root_hi],
                    "K_interval_directed_rounding": [saw_lo, saw_hi],
                    "certified_decimal_lower": new_lower,
                    "mpmath_iv_dps": DPS,
                },
                "improves_incumbent": True,
            },
            {
                "id": "macdonald_mu_4_7114_source_audit",
                "side": "lower",
                "classification": "[EXTERNAL] abstract claim; rejected as repository certificate",
                "theorem": "publisher abstract calls mu<4.7114 a rigorous upper bound, but the proof was not accessible for audit",
                "primary_source": {
                    "key": "macdonald2000_saw",
                    "citation": "D. MacDonald et al., Self-avoiding walks on the simple cubic lattice, J. Phys. A 33 (2000) 5973-5983",
                    "doi": "10.1088/0305-4470/33/34/303",
                    "exact_location": "abstract, final sentence only",
                    "verbatim": "A new, improved rigorous upper bound for the connective constant µ<4.7114 is obtained.",
                },
                "certificate": {
                    "reported_mu_upper": "47114/10000",
                    "implied_K_interval_if_the_proof_is_verified": [macdonald_lo, macdonald_hi],
                    "gap": "full text and the theorem/data proving 4.7114 were not obtained; an abstract assertion alone is insufficient for a headline certificate",
                    "decision": "not used in any certified endpoint",
                },
                "improves_incumbent": False,
                "reason": "source audit incomplete, so the stronger numerical candidate is rejected rather than promoted",
            },
            {
                "id": "finite_count_saw_c14",
                "side": "lower",
                "classification": "[THEOREM] plus [COMPUTATION] imported incumbent",
                "theorem": "SAW domination and submultiplicativity mu<=c_14^(1/14)",
                "certificate": incumbent["data"]["lower_bound_saw"],
                "improves_incumbent": False,
                "reason": "superseded by the published exact c_36 count",
            },
            {
                "id": "simon_lieb_4x4x32",
                "side": "lower",
                "classification": "[THEOREM] finite-box certificate; [COMPUTATION] diagnostics",
                "theorem": "modified Simon--Lieb finite-set criterion kappa_B(K)<1 implies K<K_c",
                "certificate": simon_lieb["data"][
                    "best_certified_simon_lieb_lower_bound"
                ],
                "diagnosis": simon_diagnostic,
                "improves_incumbent": False,
                "reason": "the 4x4 transverse boundary, not the end caps, controls kappa; exact symmetry grouping and positive tail domination cannot change the criterion",
            },
            {
                "id": "griffiths_simon_finite_region",
                "side": "lower",
                "classification": "[LEMMA] route-equivalence diagnosis",
                "theorem": "the tested free finite-region Griffiths/Simon separator inequality reduces to the same nonnegative boundary mass kappa_B",
                "certificate": {
                    "scope": "free rectangular regions with the modified Simon--Lieb crossing kernel",
                    "result": "no independent sharper criterion was identified; the strongest imported exact region is 4x4x32",
                },
                "improves_incumbent": False,
                "reason": "renaming/grouping the same positive crossing sum is not a mathematical sharpening",
            },
            {
                "id": "infrared_watson",
                "side": "upper",
                "classification": "[THEOREM] plus imported directed-rounding certificate",
                "theorem": "Frohlich--Simon--Spencer reflection positivity / infrared bound and Parseval: K_c<=I_3/2=W_sc/6",
                "certificate": incumbent["data"]["upper_bound_infrared"],
                "sharpness_control": "the nonzero-momentum Gaussian-domination constant and the exact spin sum rule; the Watson integral evaluation is already exact and is not the source of slack",
                "improves_incumbent": False,
                "reason": "this is the incumbent upper endpoint",
            },
            {
                "id": "peierls_surface_count",
                "side": "upper",
                "classification": "[THEOREM] plus [COMPUTATION] exact-rational root bracket",
                "theorem": "Bonati higher-dimensional Peierls contour estimate: <N_->/N <= 3*x^3*(9-11*x+4*x^2)/(16*(1-x)^3), x=9*exp(-4K)",
                "primary_source": {
                    "key": "bonati2014_peierls",
                    "citation": "C. Bonati, The Peierls argument for higher dimensional Ising models (2014)",
                    "arxiv": "1401.7894",
                    "local_path": "sources/fulltext/bonati2014_peierls.pdf",
                    "sha256": "2bd9daa39aad312c9a1141b404a9de19020c85e97d0b779644ddfe71567acd36",
                    "exact_location": "3D equations (17)-(23)",
                },
                "certificate": {
                    "x_root_exact_rational_bracket": [p_root_lo, p_root_hi],
                    "K_interval_directed_rounding": [peierls_lo, peierls_hi],
                    "certified_decimal_upper": peierls_upper,
                    "criterion": "the displayed minus-spin-density upper bound is <1/2 for K above the root, proving spontaneous magnetization and K_c no larger than the root",
                    "mpmath_iv_dps": DPS,
                    "exact_surface_prefix_from_repo": {
                        "source": "ising.transfer_matrix.box_broken_bond_poly((3,3,4), plus_boundary=True)",
                        "coefficients_area_0_to_28_nonzero": {
                            "0": 1,
                            "6": 36,
                            "10": 75,
                            "12": 555,
                            "14": 250,
                            "16": 2102,
                            "18": 5697,
                            "20": 8412,
                            "22": 27806,
                            "24": 58192,
                            "26": 122144,
                            "28": 276078,
                        },
                        "scope": "finite plus-boundary broken-bond configurations, not infinite-volume connected contours",
                    },
                },
                "improves_incumbent": False,
                "reason": "K_c<=0.7488... is far weaker; finite exact surface coefficients cannot replace the theorem's infinite tail without a proved contour-tail domination",
            },
        ],
        "final_certified_interval": {
            "exact": "[atanh(2941370856334701726560670^(-1/36)), W_sc/6]",
            "decimal_outward_rounded": [new_lower, new_upper],
            "lower_improved": True,
            "upper_improved": False,
            "incumbent_width": str(old_width),
            "final_width": str(new_width),
            "lower_endpoint_improvement": str(lower_improvement),
            "benchmark": str(BENCHMARK),
            "benchmark_role": "comparison-only sanity check; never used to select, fit, or round an endpoint",
        },
    }
    result = {
        "provenance": {
            "script": SCRIPT,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "precision": {
                "mpmath_iv_dps": DPS,
                "reported_decimal_places": DECIMAL_PLACES,
                "combinatorial": "exact Python integers/rationals",
                "simon_lieb_orbit_diagnostic": "numpy binary64; descriptive only and excluded from endpoints",
            },
        },
        "data": data,
        "checks": checks,
    }
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    RESULT_PATH.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        "PASS final certified interval:",
        f"[{new_lower}, {new_upper}]",
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"FAIL: {exc}")
        raise
