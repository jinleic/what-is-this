"""Exact-vs-certified Gate-D cross-check for BB [[36,4,4]] code capacity.

This is small-code, same-family reference evidence only.  It is not a
[[144,12,12]] headline result and must not be interpreted as one.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

from gatec_exact.exact_ml import ExactMLReference, dem_matrices_raw
from gatec_exact.run_gatec import build_instance, codecap_circuit_zmem
from qldpc_dec.seeds import derive_seed

from .ais import AISEngine
from .certificate import paired_bootstrap
from .gf2 import solve_gf2
from .model import GateDModel

SCOPE = "[[36,4,4]] code-capacity"
DEFAULT_BASE_SEED = 20260830


def _bits_to_key(bits: np.ndarray) -> int:
    """Little-endian logical-bit vector to the ExactMLReference class key."""
    return sum(int(bit) << i for i, bit in enumerate(np.asarray(bits).ravel()))


def _fraction_payload(value: Fraction) -> dict[str, str]:
    """Lossless JSON representation; strings avoid JSON integer precision loss."""
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
    }


def _log_fraction(value: Fraction) -> float:
    if value <= 0:
        raise ValueError("class masses must be positive")
    return math.log(value.numerator) - math.log(value.denominator)


def _summary(values: list[float]) -> dict[str, int | float | None]:
    if not values:
        return {"count": 0, "min": None, "median": None, "mean": None, "max": None}
    arr = np.asarray(values, dtype=np.float64)
    return {
        "count": len(values),
        "min": float(arr.min()),
        "median": float(np.median(arr)),
        "mean": float(arr.mean()),
        "max": float(arr.max()),
    }


def _rate(count: int, total: int) -> float | None:
    return count / total if total else None


def _lighten_representatives(
    reps: np.ndarray,
    model: GateDModel,
    generator_order: np.ndarray,
    max_sweeps: int,
) -> None:
    """Greedily lighten arbitrary representatives using only the trivial kernel."""
    for rep in reps:
        for _ in range(max_sweeps):
            improved = False
            for generator in generator_order:
                start = model.row_ptr[generator]
                stop = model.row_ptr[generator + 1]
                mechs = model.row_mechs[start:stop]
                if mechs.size == 0:
                    continue
                delta = float(model.lam[mechs].sum()) - 2.0 * float(
                    np.dot(rep[mechs], model.lam[mechs])
                )
                if delta < -1e-9:
                    rep[mechs] ^= 1
                    improved = True
            if not improved:
                break


def _record_strict_failure(
    failures: list[dict[str, object]],
    shot: int,
    category: str,
    message: str,
    strict: bool,
) -> None:
    failure = {"shot": shot, "category": category, "message": message}
    failures.append(failure)
    if strict:
        raise AssertionError(f"shot {shot}: {message}")


def run_crosscheck(
    p: float = 0.03,
    shots: int = 200,
    T: int = 64,
    K: int = 64,
    B: int = 2500,
    *,
    q0: float = 0.02,
    delta: float = 0.05,
    base_seed: int = DEFAULT_BASE_SEED,
    lighten: bool = True,
    lighten_sweeps: int = 8,
    strict: bool = False,
) -> dict[str, object]:
    """Run the deterministic 16-class exact-vs-AIS comparison."""
    if not 0.0 < p < 0.75:
        raise ValueError("p must satisfy 0 < p < 0.75")
    if shots <= 0 or T <= 0 or K <= 0 or B <= 0:
        raise ValueError("shots, T, K, and B must all be positive")
    if not 0.0 < q0 < 1.0:
        raise ValueError("q0 must satisfy 0 < q0 < 1")
    if not 0.0 < delta < 1.0:
        raise ValueError("delta must satisfy 0 < delta < 1")
    if lighten_sweeps < 0:
        raise ValueError("lighten_sweeps must be non-negative")

    total_start = time.perf_counter()

    code_start = time.perf_counter()
    inst = build_instance()
    if (inst["n"], inst["k"], inst["dx"], inst["dz"]) != (36, 4, 4, 4):
        raise AssertionError("exact cross-check is restricted to BB [[36,4,4]]")
    q = 2.0 * p / 3.0
    circuit = codecap_circuit_zmem(q, inst)
    dem = circuit.detector_error_model()
    H, A, pvec = dem_matrices_raw(dem)
    if not np.allclose(pvec, q):
        raise AssertionError("raw DEM mechanism probabilities do not match q=2p/3")
    code_s = time.perf_counter() - code_start

    model_start = time.perf_counter()
    lam = np.log1p(-pvec) - np.log(pvec)
    exact = ExactMLReference(H, A)
    model = GateDModel(H, A, lam)
    model_s = time.perf_counter() - model_start

    num_classes = 1 << inst["k"]
    expected_class_keys = list(range(num_classes))
    selector_bits = (
        (np.arange(num_classes, dtype=np.uint64)[:, None]
         >> np.arange(inst["k"], dtype=np.uint64)[None, :])
        & 1
    ).astype(np.uint8)
    logical_offsets = ((selector_bits @ model.log_shifts) & 1).astype(np.uint8)
    generator_order = np.argsort(model.K_weights, kind="stable")
    H64 = H.astype(np.uint64)
    A64 = A.astype(np.uint64)

    sampler_seed = derive_seed(
        base_seed, "gateD", "exact_crosscheck", "sampler", p
    )
    sample_start = time.perf_counter()
    syndromes, observables = circuit.compile_detector_sampler(seed=sampler_seed).sample(
        shots, separate_observables=True
    )
    syndromes = np.asarray(syndromes, dtype=np.uint8)
    observables = np.asarray(observables, dtype=np.uint8)
    sampling_s = time.perf_counter() - sample_start

    exact_hist = {str(key): 0 for key in expected_class_keys}
    ais_hist = {str(key): 0 for key in expected_class_keys}
    representative_index_hist = {str(index): 0 for index in expected_class_keys}
    records: list[dict[str, object]] = []
    strict_failures: list[dict[str, object]] = []
    exact_margins: list[float] = []
    ais_margins: list[float] = []
    certificate_lower_bounds: list[float] = []
    certified_ais_margins: list[float] = []

    agreement_count = 0
    canonical_tiebreak_agreement_count = 0
    certified_count = 0
    certified_and_exact_count = 0
    certified_wrong_count = 0
    valid_count = 0
    exact_s = 0.0
    representatives_s = 0.0
    ais_s = 0.0
    bootstrap_s = 0.0

    loop_start = time.perf_counter()
    for shot in range(shots):
        shot_start = time.perf_counter()
        syndrome = syndromes[shot]
        true_class_key = _bits_to_key(observables[shot])
        ais_seed = derive_seed(
            base_seed, "gateD", "exact_crosscheck", "ais", p, shot
        )
        bootstrap_seed = derive_seed(
            base_seed, "gateD", "exact_crosscheck", "bootstrap", p, shot
        )
        record: dict[str, object] = {
            "shot": shot,
            "status": "pending",
            "syndrome": syndrome.astype(int).tolist(),
            "true_class_key": true_class_key,
            "seeds": {"ais": ais_seed, "bootstrap": bootstrap_seed},
        }

        exact_start = time.perf_counter()
        histogram = exact.class_histogram(syndrome)
        masses = exact.class_masses(histogram, p)
        exact_class_key, exact_best_mass = ExactMLReference.ml_decision(masses)
        runner_up_key = max(
            (key for key in expected_class_keys if key != exact_class_key),
            key=lambda key: masses[key],
        )
        exact_margin_ratio = exact_best_mass / masses[runner_up_key]
        exact_margin_nats = _log_fraction(exact_margin_ratio)
        maximizer_keys = [
            key for key, mass in enumerate(masses) if mass == exact_best_mass
        ]
        shot_exact_s = time.perf_counter() - exact_start
        exact_s += shot_exact_s
        exact_margins.append(exact_margin_nats)
        exact_hist[str(exact_class_key)] += 1
        record["exact"] = {
            "decision_class_key": exact_class_key,
            "maximizer_class_keys": maximizer_keys,
            "runner_up_class_key": runner_up_key,
            "margin_nats": exact_margin_nats,
            "margin_ratio": _fraction_payload(exact_margin_ratio),
            "unnormalized_class_masses_by_key": [
                _fraction_payload(mass) for mass in masses
            ],
            "matches_sampled_logical": exact_class_key == true_class_key,
        }

        reps_start = time.perf_counter()
        e0 = solve_gf2(H, syndrome)
        if e0 is None:
            message = "sampled syndrome has no syndrome-consistent representative"
            record.update({
                "status": "invalid_syndrome",
                "syndrome_coverage": {
                    "all_representatives_consistent": False,
                    "inconsistent_representative_indices": expected_class_keys,
                },
                "class_coverage": {
                    "complete": False,
                    "actual_class_key_by_representative_index": [],
                    "missing_class_keys": expected_class_keys,
                    "duplicate_class_keys": {},
                },
                "strict_failure": message,
            })
            _record_strict_failure(
                strict_failures, shot, "syndrome_coverage", message, strict
            )
            shot_reps_s = time.perf_counter() - reps_start
            representatives_s += shot_reps_s
            record["runtime_s"] = {
                "exact": shot_exact_s,
                "representatives": shot_reps_s,
                "ais": 0.0,
                "bootstrap": 0.0,
                "total": time.perf_counter() - shot_start,
            }
            records.append(record)
            continue

        unlightened_reps = e0[None, :] ^ logical_offsets
        prelight_logical_bits = (
            (A64 @ unlightened_reps.T.astype(np.uint64)) & 1
        ).T.astype(np.uint8)
        prelight_class_keys = [
            _bits_to_key(bits) for bits in prelight_logical_bits
        ]
        e0_class_key = _bits_to_key((A64 @ e0.astype(np.uint64)) & 1)
        expected_mapping = [
            e0_class_key ^ selector_key for selector_key in expected_class_keys
        ]

        reps = unlightened_reps.copy()
        if lighten and lighten_sweeps:
            _lighten_representatives(
                reps, model, generator_order, max_sweeps=lighten_sweeps
            )

        rep_syndromes = ((H64 @ reps.T.astype(np.uint64)) & 1).T.astype(np.uint8)
        inconsistent_indices = np.flatnonzero(
            np.any(rep_syndromes != syndrome[None, :], axis=1)
        ).astype(int).tolist()
        syndrome_coverage_ok = not inconsistent_indices

        logical_bits = ((A64 @ reps.T.astype(np.uint64)) & 1).T.astype(np.uint8)
        actual_class_keys = [_bits_to_key(bits) for bits in logical_bits]
        class_counts = np.bincount(actual_class_keys, minlength=num_classes)
        missing_class_keys = np.flatnonzero(class_counts == 0).astype(int).tolist()
        duplicate_class_keys = {
            str(key): int(count)
            for key, count in enumerate(class_counts)
            if count > 1
        }
        class_preserved = actual_class_keys == prelight_class_keys
        selector_mapping_ok = prelight_class_keys == expected_mapping
        class_coverage_ok = (
            not missing_class_keys
            and not duplicate_class_keys
            and class_preserved
            and selector_mapping_ok
        )
        shot_reps_s = time.perf_counter() - reps_start
        representatives_s += shot_reps_s

        record["representatives"] = {
            "lightened": lighten,
            "lighten_sweeps": lighten_sweeps if lighten else 0,
            "e0_weight": int(e0.sum()),
            "e0_class_key": e0_class_key,
            "selector_key_by_representative_index": expected_class_keys,
            "expected_class_key_by_representative_index": expected_mapping,
            "prelight_class_key_by_representative_index": prelight_class_keys,
            "actual_class_key_by_representative_index": actual_class_keys,
            "weights_by_representative_index": reps.sum(axis=1).astype(int).tolist(),
        }
        record["syndrome_coverage"] = {
            "all_representatives_consistent": syndrome_coverage_ok,
            "inconsistent_representative_indices": inconsistent_indices,
        }
        record["class_coverage"] = {
            "complete": class_coverage_ok,
            "all_16_keys_present_once": not missing_class_keys
            and not duplicate_class_keys,
            "lightening_preserved_class": class_preserved,
            "logical_shift_selector_mapping_ok": selector_mapping_ok,
            "actual_class_key_by_representative_index": actual_class_keys,
            "missing_class_keys": missing_class_keys,
            "duplicate_class_keys": duplicate_class_keys,
        }

        shot_failures: list[str] = []
        if not syndrome_coverage_ok:
            message = (
                "representatives broke syndrome coverage at indices "
                f"{inconsistent_indices}"
            )
            shot_failures.append(message)
            _record_strict_failure(
                strict_failures, shot, "syndrome_coverage", message, strict
            )
        if not class_coverage_ok:
            message = (
                "representatives broke logical-class coverage or mapping: "
                f"missing={missing_class_keys}, duplicates={duplicate_class_keys}, "
                f"class_preserved={class_preserved}, "
                f"selector_mapping_ok={selector_mapping_ok}"
            )
            shot_failures.append(message)
            _record_strict_failure(
                strict_failures, shot, "class_coverage", message, strict
            )
        if shot_failures:
            record["status"] = "invalid_representatives"
            record["strict_failures"] = shot_failures
            record["runtime_s"] = {
                "exact": shot_exact_s,
                "representatives": shot_reps_s,
                "ais": 0.0,
                "bootstrap": 0.0,
                "total": time.perf_counter() - shot_start,
            }
            records.append(record)
            continue

        ais_start = time.perf_counter()
        engine = AISEngine(
            model.row_mechs,
            model.row_ptr,
            lam,
            T=T,
            K=K,
            q0=q0,
            seed=ais_seed,
        )
        ais_output = engine.run(reps)
        shot_ais_s = time.perf_counter() - ais_start
        ais_s += shot_ais_s

        bootstrap_start = time.perf_counter()
        certificate = paired_bootstrap(
            ais_output["logW"],
            B=B,
            delta=delta,
            rng=np.random.default_rng(bootstrap_seed),
        )
        shot_bootstrap_s = time.perf_counter() - bootstrap_start
        bootstrap_s += shot_bootstrap_s

        best_rep_index = int(certificate["best"])
        ais_class_key = actual_class_keys[best_rep_index]
        canonical_tiebreak_agreement = ais_class_key == exact_class_key
        agreement = ais_class_key in maximizer_keys
        certified = bool(certificate["certified"])
        certified_and_exact = certified and agreement
        certified_wrong = certified and not agreement
        exact_regret_ratio = exact_best_mass / masses[ais_class_key]
        exact_regret_nats = _log_fraction(exact_regret_ratio)

        logz_by_key: list[float | None] = [None] * num_classes
        for rep_index, class_key in enumerate(actual_class_keys):
            logz_by_key[class_key] = float(certificate["logZ"][rep_index])

        valid_count += 1
        canonical_tiebreak_agreement_count += int(canonical_tiebreak_agreement)
        agreement_count += int(agreement)
        certified_count += int(certified)
        certified_and_exact_count += int(certified_and_exact)
        certified_wrong_count += int(certified_wrong)
        ais_hist[str(ais_class_key)] += 1
        representative_index_hist[str(best_rep_index)] += 1
        ais_margins.append(float(certificate["margin_nats"]))
        certificate_lower_bounds.append(float(certificate["worst_lo_q"]))
        if certified:
            certified_ais_margins.append(float(certificate["margin_nats"]))

        record.update({
            "status": "ok",
            "ais": {
                "decision_representative_index": best_rep_index,
                "decision_class_key": ais_class_key,
                "logZ_by_class_key": logz_by_key,
                "margin_nats": float(certificate["margin_nats"]),
                "exact_regret_nats": exact_regret_nats,
                "exact_regret_ratio": _fraction_payload(exact_regret_ratio),
                "matches_sampled_logical": ais_class_key == true_class_key,
            },
            "certificate": {
                "certified": certified,
                "worst_lower_bound_nats": float(certificate["worst_lo_q"]),
                "bonferroni_level": float(certificate["bonf_level"]),
                "bootstrap_failure_probability": float(certificate["boot_p_fail"]),
            },
            "agreement": agreement,
            "canonical_tiebreak_agreement": canonical_tiebreak_agreement,
            "certified_and_exact": certified_and_exact,
            "certified_wrong": certified_wrong,
            "runtime_s": {
                "exact": shot_exact_s,
                "representatives": shot_reps_s,
                "ais": shot_ais_s,
                "bootstrap": shot_bootstrap_s,
                "total": time.perf_counter() - shot_start,
            },
        })
        records.append(record)

        if certified_wrong:
            _record_strict_failure(
                strict_failures,
                shot,
                "certified_wrong",
                (
                    f"certified AIS class {ais_class_key} is outside the "
                    f"exact ML maximizer set {maximizer_keys}"
                ),
                strict,
            )

    loop_s = time.perf_counter() - loop_start
    total_s = time.perf_counter() - total_start
    invalid_count = shots - valid_count

    return {
        "scope": SCOPE,
        "scope_note": (
            "Same-family BB [[36,4,4]] small-code reference only; "
            "not a [[144,12,12]] headline verdict."
        ),
        "code": {
            "family": "BB (l=6, m=3)",
            "n": inst["n"],
            "k": inst["k"],
            "dx": inst["dx"],
            "dz": inst["dz"],
            "sector": "Z-memory",
            "num_detectors": int(H.shape[0]),
            "num_mechanisms": int(H.shape[1]),
            "num_logical_classes": num_classes,
        },
        "noise": {
            "model": "code-capacity depolarizing, X-error marginal",
            "p_circuit_equiv": p,
            "q_x_marginal": q,
            "dem_mechanism_probability_min": float(pvec.min()),
            "dem_mechanism_probability_max": float(pvec.max()),
        },
        "sampling": {
            "shots_requested": shots,
            "shots_with_valid_class_coverage": valid_count,
            "shots_invalid": invalid_count,
            "base_seed": base_seed,
            "sampler_seed": sampler_seed,
            "seed_derivation": "qldpc_dec.seeds.derive_seed",
        },
        "exact_config": {
            "arithmetic": "fractions.Fraction class masses and argmax",
            "rank_H": exact.rank,
            "affine_space_log2": exact.dim,
            "solutions_per_syndrome": 1 << exact.dim,
        },
        "ais_config": {
            "T": T,
            "K": K,
            "q0": q0,
            "beta_schedule": "power2",
            "B": B,
            "delta": delta,
            "bonferroni_competitors": num_classes - 1,
            "num_trivial_kernel_generators": int(model.G),
            "mean_generator_weight": float(model.K_weights.mean()),
            "lighten_representatives": lighten,
            "lighten_sweeps": lighten_sweeps if lighten else 0,
        },
        "strict_mode": strict,
        "strict_failures": strict_failures,
        "exact_decision_histogram": exact_hist,
        "ais_decision_histogram": ais_hist,
        "ais_representative_index_histogram": representative_index_hist,
        "agreement_count": agreement_count,
        "agreement_rate": _rate(agreement_count, valid_count),
        "canonical_tiebreak_agreement_count": canonical_tiebreak_agreement_count,
        "canonical_tiebreak_agreement_rate": _rate(
            canonical_tiebreak_agreement_count, valid_count
        ),
        "certified_count": certified_count,
        "certified_rate": _rate(certified_count, valid_count),
        "certified_and_exact_count": certified_and_exact_count,
        "certified_and_exact_rate": _rate(certified_and_exact_count, valid_count),
        "certified_wrong_count": certified_wrong_count,
        "certified_wrong_rate_among_certified": _rate(
            certified_wrong_count, certified_count
        ),
        "margins": {
            "exact_ml_nats": _summary(exact_margins),
            "ais_point_estimate_nats": _summary(ais_margins),
            "certificate_worst_lower_bound_nats": _summary(
                certificate_lower_bounds
            ),
            "certified_ais_point_estimate_nats": _summary(
                certified_ais_margins
            ),
        },
        "runtime": {
            "code_and_dem_build_s": code_s,
            "exact_and_gateD_model_build_s": model_s,
            "gateD_basis_build_s": float(model.basis_seconds),
            "sampling_s": sampling_s,
            "shot_loop_s": loop_s,
            "exact_s": exact_s,
            "representatives_s": representatives_s,
            "ais_s": ais_s,
            "bootstrap_s": bootstrap_s,
            "total_s": total_s,
            "total_s_per_requested_shot": total_s / shots,
        },
        "per_shot": records,
    }


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Exact-vs-certified Gate-D cross-check for the BB [[36,4,4]] "
            "code-capacity reference"
        )
    )
    parser.add_argument("--p", type=float, default=0.03)
    parser.add_argument("--shots", type=int, default=200)
    parser.add_argument("--T", type=int, default=64)
    parser.add_argument("--K", type=int, default=64)
    parser.add_argument("--B", type=int, default=2500)
    parser.add_argument("--q0", type=float, default=0.02)
    parser.add_argument("--delta", type=float, default=0.05)
    parser.add_argument("--base-seed", type=int, default=DEFAULT_BASE_SEED)
    parser.add_argument("--lighten-sweeps", type=int, default=8)
    parser.add_argument(
        "--no-lighten",
        dest="lighten",
        action="store_false",
        help="disable deterministic trivial-kernel representative lightening",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help=(
            "raise on broken syndrome/class coverage or any certified "
            "decision that disagrees with exact ML"
        ),
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.set_defaults(lighten=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    report = run_crosscheck(
        p=args.p,
        shots=args.shots,
        T=args.T,
        K=args.K,
        B=args.B,
        q0=args.q0,
        delta=args.delta,
        base_seed=args.base_seed,
        lighten=args.lighten,
        lighten_sweeps=args.lighten_sweeps,
        strict=args.strict,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(args.out)


if __name__ == "__main__":
    main()
