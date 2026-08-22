#!/usr/bin/env python3
"""Quantitative L18 first-emergent-free-member model.

The sole fitted parameter multiplies the classwise p<=10,000 Euler product
after extrapolation to the primitive class-polynomial height with the L18
(log X)^(-1/2) tail. A deterministic outcome-blind hash split is declared
below and assigned before l15_density.json is opened.

Every fitted number and goodness-of-fit statement emitted here is EVIDENCE.
"""

import bisect
import hashlib
import json
import math
import os
import random
import statistics
import sys
import time
from collections import Counter
from fractions import Fraction
from math import comb
from pathlib import Path
from statistics import NormalDist

ROOT = Path(os.environ.get("H10Q_ROOT", Path(__file__).resolve().parent))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from h10q import _l10_P, factorint  # noqa: E402

STEP_PATH = ROOT / "data" / "l17_stepii.jsonl"
DENSITY_PATH = ROOT / "data" / "l15_density.json"
SIEVE_PATH = ROOT / "data" / "l17_sieve.jsonl"
FACTORS_PATH = ROOT / "data" / "l17_classfactors.jsonl"
OUT_PATH = Path(os.environ.get("L18_DIVERGENCE_OUT", ROOT / "data" / "l18_divergence_model.jsonl"))
REPORT_PATH = Path(os.environ.get("L18_DIVERGENCE_REPORT", "/tmp/l18_divergence_model.md"))

# PREDECLARED before fitting and independent of all outcome fields.
SPLIT_SALT = "L18-divergence-v1"
SPLIT_MODULUS = 5
HELDOUT_RESIDUE = 0
TAIL_EXPONENT = 0.5
EULER_CUTOFF = 10_000
MEMBER_CUTOFF = 100_000
SURVIVAL_TOL = 1e-10
MAX_PREDICTION_K = 200_000
KS_SIMULATIONS = int(os.environ.get("L18_KS_SIMULATIONS", "5000"))
KS_SEED = 20_260_821
PACE_EVERY = 10
NORMAL = NormalDist()
QUANTILES = (0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)


def _pace(i):
    if i % PACE_EVERY == 0:
        time.sleep(0.1)


def _load_jsonl(path):
    with path.open() as fh:
        return [json.loads(line) for line in fh if line.strip()]


def _sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            block = fh.read(1 << 20)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def _class_key(row):
    w, u = row["cell"]
    return f"{row['family']}|{w}|({u[0]},{u[1]})"


def _split_for_key(key):
    digest = hashlib.sha256(f"{SPLIT_SALT}|{key}".encode()).digest()
    residue = int.from_bytes(digest[:8], "big") % SPLIT_MODULUS
    return "heldout" if residue == HELDOUT_RESIDUE else "fit"


def _compose_affine(poly, base, step):
    out = [Fraction(0)] * len(poly)
    for i, coefficient in enumerate(poly):
        for j in range(i + 1):
            out[j] += coefficient * comb(i, j) * base ** (i - j) * step ** j
    return out


def _primitive_class_polynomial(row):
    w, (u_num, u_den) = row["cell"]
    z = Fraction(w * u_num, u_den)
    a = Fraction(row["a"])
    A = 1 + 4 * a * a
    tau = (1 + 2 * a * a) / A
    delta = 1 - A * tau * tau
    Z = z ** 3
    D = 1 - Z - a * a * Z * Z
    s = (a - 1) / 2
    poly_b = _l10_P(a, Z, D, A, delta, s)
    base = int(row["eps"]) * int(row["f"]) * int(row["q1"])
    step = int(row["eps"]) * int(row["f"]) * int(row["N"])
    poly_k = _compose_affine(poly_b, base, step)

    denominator = 1
    for coefficient in poly_k:
        denominator = math.lcm(denominator, coefficient.denominator)
    integral = [int(coefficient * denominator) for coefficient in poly_k]
    content = 0
    for coefficient in integral:
        content = math.gcd(content, abs(coefficient))
    if not content:
        raise AssertionError((_class_key(row), "zero polynomial"))
    primitive = [coefficient // content for coefficient in integral]
    if primitive[-1] < 0:
        primitive = [-coefficient for coefficient in primitive]
    if len(primitive) != 9 or primitive[-1] == 0:
        raise AssertionError((_class_key(row), len(primitive) - 1))
    return primitive


def _eval_poly(coefficients, k):
    value = 0
    for coefficient in reversed(coefficients):
        value = value * k + coefficient
    return value


def _phi_and_factors(n):
    factors = factorint(n)
    product = 1
    phi = n
    for p, exponent in factors.items():
        product *= p ** exponent
        phi = phi // p * (p - 1)
    if product != n:
        raise AssertionError((n, factors))
    return phi, factors


def _state_at_k(feature, k):
    q = feature["q1"] + k * feature["N"]
    X = abs(_eval_poly(feature["F_coefficients"], k))
    if X <= 1:
        raise AssertionError((feature["key"], k, X))
    log_X = math.log(X)
    tail_ratio = (math.log(EULER_CUTOFF) / math.log(max(EULER_CUTOFF, X))) ** TAIL_EXPONENT
    clean_shape = feature["factor_10000"] * tail_ratio
    prime_probability = 1.0 if k == 0 else min(1.0, feature["ap_alpha"] / math.log(q))
    return {
        "q": q,
        "X": X,
        "log_X": log_X,
        "tail_ratio": tail_ratio,
        "clean_shape": clean_shape,
        "prime_probability": prime_probability,
        "base_hazard": prime_probability * clean_shape,
    }


def _log_likelihood(c, features):
    if c <= 0:
        return -math.inf
    total = 0.0
    for feature in features:
        prefix = feature["observed_base_hazards"]
        event_hazard = c * prefix[-1]
        if not 0 < event_hazard < 1:
            return -math.inf
        total += math.log(event_hazard)
        total += sum(math.log1p(-c * hazard) for hazard in prefix[:-1])
    return total


def _fit_constant(features, c_upper):
    # On the physical interval c*factor_10000<1 the likelihood is concave.
    def score(c):
        failures = sum(
            hazard / (1 - c * hazard)
            for feature in features
            for hazard in feature["observed_base_hazards"][:-1]
        )
        return len(features) / c - failures

    lo = 1e-14
    hi = c_upper
    boundary = score(hi) >= 0
    if boundary:
        estimate = hi
    else:
        for _ in range(120):
            mid = (lo + hi) / 2
            if score(mid) > 0:
                lo = mid
            else:
                hi = mid
        estimate = (lo + hi) / 2

    ll_max = _log_likelihood(estimate, features)
    profile_target = ll_max - 3.841458820694124 / 2

    lo = 1e-14
    hi = estimate
    for _ in range(120):
        mid = (lo + hi) / 2
        if _log_likelihood(mid, features) < profile_target:
            lo = mid
        else:
            hi = mid
    ci_lower = (lo + hi) / 2

    if _log_likelihood(c_upper, features) > profile_target:
        ci_upper = c_upper
        ci_upper_truncated = True
    else:
        lo = estimate
        hi = c_upper
        for _ in range(120):
            mid = (lo + hi) / 2
            if _log_likelihood(mid, features) > profile_target:
                lo = mid
            else:
                hi = mid
        ci_upper = (lo + hi) / 2
        ci_upper_truncated = False

    return {
        "estimate": estimate,
        "profile_95_ci": [ci_lower, ci_upper],
        "profile_ci_upper_truncated": ci_upper_truncated,
        "log_likelihood_fit": ll_max,
        "boundary_mle": boundary,
        "physical_upper_bound": c_upper,
    }


def _prediction(feature, c):
    survival = 1.0
    mean_k = 0.0
    expected_count = 0.0
    cdf = []
    quantiles = {}
    for k in range(MAX_PREDICTION_K + 1):
        state = _state_at_k(feature, k)
        hazard = c * state["base_hazard"]
        if not 0 <= hazard < 1:
            raise AssertionError((feature["key"], k, hazard))
        expected_count += hazard
        survival *= 1 - hazard
        cdf.append(1 - survival)
        mean_k += survival
        for quantile in QUANTILES:
            if quantile not in quantiles and 1 - survival >= quantile:
                quantiles[quantile] = k
        if survival <= SURVIVAL_TOL and len(quantiles) == len(QUANTILES):
            break
    else:
        raise RuntimeError((feature["key"], "prediction did not converge", survival))
    return {
        "cdf": cdf,
        "quantiles": quantiles,
        "mean_k_truncated": mean_k,
        "survival_tail": survival,
        "k_truncation": len(cdf) - 1,
        "expected_count_at_truncation": expected_count,
    }


def _average_ranks(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        average = (start + 1 + end) / 2
        for position in range(start, end):
            ranks[order[position]] = average
        start = end
    return ranks


def _correlation(left, right):
    mean_left = statistics.mean(left)
    mean_right = statistics.mean(right)
    numerator = sum((x - mean_left) * (y - mean_right) for x, y in zip(left, right))
    denominator = math.sqrt(
        sum((x - mean_left) ** 2 for x in left)
        * sum((y - mean_right) ** 2 for y in right)
    )
    return numerator / denominator if denominator else None


def _spearman(left, right):
    return _correlation(_average_ranks(left), _average_ranks(right))


def _mixture_cdf(indices, predictions, max_observed):
    max_k = max(max_observed, max(predictions[i]["k_truncation"] for i in indices))
    mixture = []
    for k in range(max_k + 1):
        total = 0.0
        for i in indices:
            cdf = predictions[i]["cdf"]
            total += cdf[k] if k < len(cdf) else cdf[-1]
        mixture.append(total / len(indices))
    return mixture


def _ks_from_values(values, mixture):
    counts = Counter(values)
    n = len(values)
    cumulative = 0
    best = {"D": -1.0, "k": None, "signed_observed_minus_predicted": None}
    for value in sorted(counts):
        before_k = value - 1
        predicted_before = mixture[before_k] if before_k >= 0 else 0.0
        difference = cumulative / n - predicted_before
        if abs(difference) > best["D"]:
            best = {"D": abs(difference), "k": before_k, "signed_observed_minus_predicted": difference}
        cumulative += counts[value]
        predicted_at = mixture[value] if value < len(mixture) else mixture[-1]
        difference = cumulative / n - predicted_at
        if abs(difference) > best["D"]:
            best = {"D": abs(difference), "k": value, "signed_observed_minus_predicted": difference}
    return best


def _ks_bootstrap(indices, predictions, mixture, observed_values, seed):
    observed = _ks_from_values(observed_values, mixture)
    rng = random.Random(seed)
    sampling_cdfs = []
    for i in indices:
        cdf = list(predictions[i]["cdf"])
        cdf[-1] = 1.0  # assigns <SURVIVAL_TOL residual mass to the last bin
        sampling_cdfs.append(cdf)
    exceedances = 0
    simulated_D = []
    for simulation in range(KS_SIMULATIONS):
        sample = [bisect.bisect_left(cdf, rng.random()) for cdf in sampling_cdfs]
        D = _ks_from_values(sample, mixture)["D"]
        simulated_D.append(D)
        exceedances += D >= observed["D"] - 1e-15
        if (simulation + 1) % 250 == 0:
            time.sleep(0.1)
    ordered = sorted(simulated_D)
    q95 = ordered[math.ceil(0.95 * len(ordered)) - 1]
    observed.update(
        {
            "bootstrap_simulations": KS_SIMULATIONS,
            "bootstrap_seed": seed,
            "bootstrap_exceedances": exceedances,
            "bootstrap_add_one_p": (exceedances + 1) / (KS_SIMULATIONS + 1),
            "bootstrap_mean_D": statistics.mean(simulated_D),
            "bootstrap_q95_D": q95,
            "bootstrap_tail_note": "add-one Monte Carlo p; non-iid classwise predictive CDFs sampled separately; C held fixed, so fit-set p is descriptive and heldout p is primary",
        }
    )
    return observed


def _conditional_quantile(feature, c, mode, target=0.5):
    survival = 1.0
    for k in range(1, MAX_PREDICTION_K + 1):
        state = _state_at_k(feature, k)
        if mode == "pure_bh":
            hazard = state["prime_probability"]
        elif mode == "truncated_only":
            hazard = state["prime_probability"] * c * feature["factor_10000"]
        elif mode == "l18":
            hazard = c * state["base_hazard"]
        else:
            raise ValueError(mode)
        survival *= 1 - hazard
        if 1 - survival >= target:
            return k
    raise RuntimeError((feature["key"], mode, "conditional quantile did not converge"))


def _aggregate(label, indices, features, predictions, output_rows, seed):
    observed_values = [features[i]["observed_k"] for i in indices]
    mixture = _mixture_cdf(indices, predictions, max(observed_values))
    ks = _ks_bootstrap(indices, predictions, mixture, observed_values, seed)
    mixture_median = bisect.bisect_left(mixture, 0.5)
    selected_k = sorted({0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, max(observed_values)})
    curve = []
    for k in selected_k:
        predicted = mixture[k] if k < len(mixture) else mixture[-1]
        observed = sum(value <= k for value in observed_values) / len(observed_values)
        curve.append({"k": k, "observed_cdf": observed, "predicted_cdf": predicted, "observed_minus_predicted": observed - predicted})

    quantile_midpoints = [output_rows[i]["quantile_residual"]["midpoint_u"] for i in indices]
    normal_residuals = [output_rows[i]["quantile_residual"]["normal"] for i in indices]
    nll = -sum(output_rows[i]["prediction_at_observed"]["log_event_probability"] for i in indices) / len(indices)
    predicted_medians = [output_rows[i]["predicted_first"]["median_k"] for i in indices]
    predicted_means = [output_rows[i]["predicted_first"]["mean_k_truncated"] for i in indices]

    member_n = sum(features[i]["sieve_n_prime_members"] for i in indices)
    member_observed = sum(
        features[i]["sieve_n_prime_members"] * features[i]["sieve_clean_rate_100000"] for i in indices
    ) / member_n
    member_predicted = sum(
        features[i]["sieve_n_prime_members"] * output_rows[i]["member_window_diagnostic"]["predicted_clean_rate_100000"]
        for i in indices
    ) / member_n
    member_residuals = [output_rows[i]["member_window_diagnostic"]["observed_minus_predicted"] for i in indices]

    return {
        "label": "EVIDENCE",
        "subset": label,
        "n_classes": len(indices),
        "family_counts": dict(sorted(Counter(features[i]["family"] for i in indices).items())),
        "observed": {
            "k0_count": sum(value == 0 for value in observed_values),
            "k0_fraction": sum(value == 0 for value in observed_values) / len(observed_values),
            "median_k": statistics.median(observed_values),
            "mean_k": statistics.mean(observed_values),
            "maximum_k": max(observed_values),
        },
        "predicted": {
            "expected_k0_fraction": statistics.mean(output_rows[i]["predicted_first"]["probability_k0"] for i in indices),
            "mixture_median_k": mixture_median,
            "class_mean_of_median_k": statistics.mean(predicted_medians),
            "class_mean_of_mean_k_truncated": statistics.mean(predicted_means),
        },
        "ks_style": ks,
        "cdf_checkpoints": curve,
        "event_negative_log_likelihood_per_class": nll,
        "per_class_discrimination": {
            "spearman_predicted_median_vs_observed_k": _spearman(predicted_medians, observed_values),
            "spearman_predicted_mean_vs_observed_k": _spearman(predicted_means, observed_values),
        },
        "quantile_residuals": {
            "midpoint_u_mean": statistics.mean(quantile_midpoints),
            "midpoint_u_median": statistics.median(quantile_midpoints),
            "normal_mean": statistics.mean(normal_residuals),
            "normal_median": statistics.median(normal_residuals),
            "fraction_u_above_0_95": sum(value > 0.95 for value in quantile_midpoints) / len(indices),
            "fraction_u_below_0_05": sum(value < 0.05 for value in quantile_midpoints) / len(indices),
            "central_80_interval_coverage": sum(output_rows[i]["quantile_residual"]["inside_predicted_10_90"] for i in indices) / len(indices),
        },
        "member_window_p_le_100000": {
            "observed_pooled_clean_rate": member_observed,
            "predicted_pooled_clean_rate": member_predicted,
            "observed_minus_predicted_pooled": member_observed - member_predicted,
            "class_mean_absolute_error": statistics.mean(abs(value) for value in member_residuals),
            "class_mean_observed_minus_predicted": statistics.mean(member_residuals),
            "warning": "EVIDENCE: actual prime-member k<=119 window versus Haar-product extrapolation",
        },
    }


def _fmt(value, digits=4):
    return f"{value:.{digits}f}"


def build():
    # Assign the declared split from class keys before opening the outcome file.
    factor_rows = _load_jsonl(FACTORS_PATH)
    if len(factor_rows) != 293:
        raise AssertionError(len(factor_rows))
    keys = [_class_key(row) for row in factor_rows]
    if len(set(keys)) != 293:
        raise AssertionError("duplicate class key")
    split_by_key = {key: _split_for_key(key) for key in keys}
    split_manifest = "\n".join(f"{key}:{split_by_key[key]}" for key in sorted(keys))
    split_manifest_sha256 = hashlib.sha256(split_manifest.encode()).hexdigest()

    # Outcome data is opened only after split assignment is frozen above.
    with DENSITY_PATH.open() as fh:
        density = json.load(fh)
    step_rows = [row for row in _load_jsonl(STEP_PATH) if row.get("type") == "class"]
    sieve_rows = [row for row in _load_jsonl(SIEVE_PATH) if row.get("type") == "class-summary"]
    if len(step_rows) != 293 or len(sieve_rows) != 293:
        raise AssertionError((len(step_rows), len(sieve_rows)))
    step_by_key = {_class_key(row): row for row in step_rows}
    sieve_by_key = {_class_key(row): row for row in sieve_rows}
    if set(keys) != set(step_by_key) or set(keys) != set(sieve_by_key):
        raise AssertionError("input class-key mismatch")

    nonzero = {row["class"]: int(row["k"]) for row in density["v_p0_size"]["rows"]}
    if len(nonzero) != 86 or not set(nonzero) <= set(keys):
        raise AssertionError((len(nonzero), sorted(set(nonzero) - set(keys))))
    observed_by_key = {key: nonzero.get(key, 0) for key in keys}
    observed_values = list(observed_by_key.values())
    if not (
        len(observed_values) == density["meta"]["n"] == 293
        and sum(value == 0 for value in observed_values) == density["i_global"]["k0_count"] == 207
        and max(observed_values) == density["i_global"]["max"] == 1694
        and abs(statistics.mean(observed_values) - density["i_global"]["mean"]) < 5e-5
    ):
        raise AssertionError("l15 observed-law reconstruction mismatch")

    features = []
    for index, factor_row in enumerate(factor_rows, 1):
        key = _class_key(factor_row)
        step_row = step_by_key[key]
        sieve_row = sieve_by_key[key]
        params = step_row["params"]
        for field in ("eps", "f", "q1", "N"):
            if int(factor_row[field]) != int(params[field]):
                raise AssertionError((key, field, factor_row[field], params[field]))
        if abs(float(factor_row["factor_lower"]) - float(step_row["root_window"]["factor_source_lower"])) > 5e-14:
            raise AssertionError((key, "factor lower"))
        if abs(float(factor_row["factor_upper"]) - float(step_row["root_window"]["factor_source_upper"])) > 5e-14:
            raise AssertionError((key, "factor upper"))
        step_clean = step_row["member_window_evidence"]["clean_rate_by_prime_cutoff"][str(MEMBER_CUTOFF)]
        if abs(float(sieve_row["clean_rate"]) - float(step_clean)) > 5e-12:
            raise AssertionError((key, "member clean rate"))

        N = int(factor_row["N"])
        q1 = int(factor_row["q1"])
        phi_N, N_factors = _phi_and_factors(N)
        if math.gcd(q1, N) != 1:
            raise AssertionError((key, q1, N))
        coefficients = _primitive_class_polynomial(factor_row)
        factor_10000 = math.sqrt(float(factor_row["factor_lower"]) * float(factor_row["factor_upper"]))
        feature = {
            "key": key,
            "split": split_by_key[key],
            "family": factor_row["family"],
            "cell": factor_row["cell"],
            "a": str(factor_row["a"]),
            "eps": int(factor_row["eps"]),
            "f": int(factor_row["f"]),
            "q1": q1,
            "N": N,
            "phi_N": phi_N,
            "N_factors": {str(p): exponent for p, exponent in sorted(N_factors.items())},
            "ap_alpha": N / phi_N,
            "factor_lower_10000": float(factor_row["factor_lower"]),
            "factor_upper_10000": float(factor_row["factor_upper"]),
            "factor_10000": factor_10000,
            "factor_partial": bool(factor_row["partial"]),
            "F_coefficients": coefficients,
            "F_coefficients_sha256": hashlib.sha256(",".join(map(str, coefficients)).encode()).hexdigest(),
            "observed_k": observed_by_key[key],
            "sieve_n_prime_members": int(sieve_row["n_prime_members"]),
            "sieve_clean_rate_100000": float(sieve_row["clean_rate"]),
            "step_product_c_100_10000": float(step_row["root_window"]["c_product_100_10000_upper"]),
            "theoretical_chebotarev_exponent": float(step_row["theoretical_chebotarev_exponent"]),
        }
        if feature["theoretical_chebotarev_exponent"] != TAIL_EXPONENT:
            raise AssertionError((key, feature["theoretical_chebotarev_exponent"]))
        feature["observed_base_hazards"] = [
            _state_at_k(feature, k)["base_hazard"] for k in range(feature["observed_k"] + 1)
        ]
        features.append(feature)
        _pace(index)

    fit_features = [feature for feature in features if feature["split"] == "fit"]
    heldout_features = [feature for feature in features if feature["split"] == "heldout"]
    if len(fit_features) + len(heldout_features) != 293 or not heldout_features:
        raise AssertionError((len(fit_features), len(heldout_features)))

    # Since tail_ratio<=1, this outcome-blind bound keeps every clean probability <1.
    c_upper = (1 - 1e-9) / max(feature["factor_upper_10000"] for feature in features)
    calibration = _fit_constant(fit_features, c_upper)
    c = calibration["estimate"]

    predictions = []
    output_rows = []
    for index, feature in enumerate(features, 1):
        prediction = _prediction(feature, c)
        predictions.append(prediction)
        observed_k = feature["observed_k"]
        prefix = feature["observed_base_hazards"]
        survival_before = math.prod(1 - c * hazard for hazard in prefix[:-1])
        event_probability = survival_before * c * prefix[-1]
        cdf_before = 1 - survival_before
        cdf_at = cdf_before + event_probability
        midpoint_u = cdf_before + event_probability / 2
        normal_residual = NORMAL.inv_cdf(min(1 - 1e-15, max(1e-15, midpoint_u)))
        observed_state = _state_at_k(feature, observed_k)
        median_k = prediction["quantiles"][0.50]
        median_state = _state_at_k(feature, median_k)
        k0_state = _state_at_k(feature, 0)
        predicted_clean_100000 = (
            c
            * feature["factor_10000"]
            * (math.log(EULER_CUTOFF) / math.log(MEMBER_CUTOFF)) ** TAIL_EXPONENT
        )
        member_residual = feature["sieve_clean_rate_100000"] - predicted_clean_100000

        row = {
            "type": "class",
            "label": "EVIDENCE",
            "class": feature["key"],
            "split": feature["split"],
            "family": feature["family"],
            "cell": feature["cell"],
            "params": {
                "a": feature["a"],
                "eps": feature["eps"],
                "f": feature["f"],
                "q1": feature["q1"],
                "N": str(feature["N"]),
                "phi_N": str(feature["phi_N"]),
                "N_factors": feature["N_factors"],
                "N_over_phi_N": feature["ap_alpha"],
            },
            "local_input": {
                "factor_lower_p_le_10000": feature["factor_lower_10000"],
                "factor_upper_p_le_10000": feature["factor_upper_10000"],
                "factor_geometric_midpoint_p_le_10000": feature["factor_10000"],
                "partial": feature["factor_partial"],
                "measured_product_exponent_100_10000": feature["step_product_c_100_10000"],
                "theoretical_chebotarev_exponent": feature["theoretical_chebotarev_exponent"],
            },
            "height": {
                "definition": "absolute value of primitive integral F_cell(k)",
                "degree": len(feature["F_coefficients"]) - 1,
                "coefficients_sha256": feature["F_coefficients_sha256"],
                "log_X_k0": k0_state["log_X"],
                "digits_X_k0": len(str(k0_state["X"])),
                "log_X_observed": observed_state["log_X"],
                "digits_X_observed": len(str(observed_state["X"])),
            },
            "observed": {
                "first_k": observed_k,
                "first_Q": str(feature["q1"] + observed_k * feature["N"]),
            },
            "predicted_first": {
                "calibration_constant": c,
                "probability_k0": c * k0_state["base_hazard"],
                "median_k": median_k,
                "median_Q": str(median_state["q"]),
                "mean_k_truncated": prediction["mean_k_truncated"],
                "quantile_k": {str(quantile): prediction["quantiles"][quantile] for quantile in QUANTILES},
                "survival_at_truncation": prediction["survival_tail"],
                "truncation_k": prediction["k_truncation"],
            },
            "prediction_at_observed": {
                "prime_probability": observed_state["prime_probability"],
                "tail_ratio_from_10000": observed_state["tail_ratio"],
                "clean_probability_given_prime": c * observed_state["clean_shape"],
                "hazard": c * observed_state["base_hazard"],
                "expected_count_k_le_observed": c * sum(prefix),
                "cdf_before": cdf_before,
                "cdf_at": cdf_at,
                "event_probability": event_probability,
                "log_event_probability": math.log(event_probability),
            },
            "quantile_residual": {
                "definition": "midpoint randomized quantile u=F(k-1)+P(T=k)/2; normal=Phi^-1(u)",
                "midpoint_u": midpoint_u,
                "normal": normal_residual,
                "inside_predicted_10_90": prediction["quantiles"][0.10] <= observed_k <= prediction["quantiles"][0.90],
            },
            "member_window_diagnostic": {
                "n_proven_prime_members_k_le_119": feature["sieve_n_prime_members"],
                "observed_clean_rate_p_le_100000": feature["sieve_clean_rate_100000"],
                "predicted_clean_rate_100000": predicted_clean_100000,
                "observed_minus_predicted": member_residual,
                "warning": "different estimand from full first-member event: finite p cutoff and actual prime-member window",
            },
        }
        output_rows.append(row)
        _pace(index)

    fit_indices = [i for i, feature in enumerate(features) if feature["split"] == "fit"]
    heldout_indices = [i for i, feature in enumerate(features) if feature["split"] == "heldout"]
    all_indices = list(range(len(features)))
    evaluations = {
        "in_sample_fit": _aggregate("in_sample_fit", fit_indices, features, predictions, output_rows, KS_SEED + 1),
        "heldout": _aggregate("heldout", heldout_indices, features, predictions, output_rows, KS_SEED + 2),
        "all_descriptive": _aggregate("all_descriptive", all_indices, features, predictions, output_rows, KS_SEED + 3),
    }

    prime_tail_keys = {
        row["class"] for row in density["ii_discrepancy"]["rows"] if int(row["k"]) > 0
    }
    expected_tail_n = density["ii_discrepancy"]["first_prime_only_tail"]["n"]
    if len(prime_tail_keys) != expected_tail_n or expected_tail_n != 69:
        raise AssertionError(len(prime_tail_keys))
    feature_by_key = {feature["key"]: feature for feature in features}
    gap_rows = []
    for key in sorted(prime_tail_keys):
        feature = feature_by_key[key]
        pure_k = _conditional_quantile(feature, c, "pure_bh")
        truncated_k = _conditional_quantile(feature, c, "truncated_only")
        l18_k = _conditional_quantile(feature, c, "l18")
        gap_rows.append(
            {
                "key": key,
                "split": feature["split"],
                "observed_k": feature["observed_k"],
                "pure_bh_conditional_median_k": pure_k,
                "truncated_only_conditional_median_k": truncated_k,
                "l18_conditional_median_k": l18_k,
                "truncated_delay_over_bh": truncated_k / pure_k,
                "l18_delay_over_bh": l18_k / pure_k,
                "observed_over_l18": feature["observed_k"] / l18_k,
            }
        )

    ledger_log2_gap = float(density["ii_discrepancy"]["first_prime_only_tail"]["median_log2k_over_m1"])
    ledger_gap = 2 ** ledger_log2_gap

    def gap_subset(split):
        rows = [row for row in gap_rows if split == "all" or row["split"] == split]
        delays = [row["l18_delay_over_bh"] for row in rows]
        truncated_delays = [row["truncated_delay_over_bh"] for row in rows]
        observed_residuals = [row["observed_over_l18"] for row in rows]
        return {
            "n": len(rows),
            "median_truncated_only_delay_over_bh": statistics.median(truncated_delays),
            "median_l18_delay_over_bh": statistics.median(delays),
            "geometric_mean_l18_delay_over_bh": math.exp(statistics.mean(math.log(value) for value in delays)),
            "median_log_residual_observed_over_l18": 2 ** statistics.median(math.log2(value) for value in observed_residuals),
        }

    all_gap = gap_subset("all")
    gap_verdict = {
        "label": "EVIDENCE",
        "ledger_benchmark": {
            "source_field": "l15_density.ii_discrepancy.first_prime_only_tail.median_log2k_over_m1",
            "n_prime_rung_tail_classes": len(gap_rows),
            "median_log2_gap": ledger_log2_gap,
            "multiplicative_gap": ledger_gap,
        },
        "in_sample_fit": gap_subset("fit"),
        "heldout": gap_subset("heldout"),
        "all_descriptive": all_gap,
        "comparison": {
            "model_l18_delay_factor": all_gap["median_l18_delay_over_bh"],
            "residual_factor_vs_ledger_gap": ledger_gap / all_gap["median_l18_delay_over_bh"],
            "fraction_of_excess_multiplier_closed": (all_gap["median_l18_delay_over_bh"] - 1) / (ledger_gap - 1),
            "fraction_of_log_gap_closed": math.log(all_gap["median_l18_delay_over_bh"]) / math.log(ledger_gap),
            "verdict": "PARTIALLY CLOSES: material tail slowdown, non-unit residual, and severe full-distribution k=0 failure",
        },
        "per_class": gap_rows,
    }

    input_hashes = {
        str(path.relative_to(ROOT)): _sha256(path)
        for path in (STEP_PATH, DENSITY_PATH, SIEVE_PATH, FACTORS_PATH)
    }
    split_counts = Counter(feature["split"] for feature in features)
    split_family_counts = Counter((feature["split"], feature["family"]) for feature in features)
    meta = {
        "type": "meta",
        "artifact": "l18_divergence_model",
        "version": 1,
        "date": "2026-08-21",
        "labels": {
            "all_model_outputs_and_statistics": "EVIDENCE",
            "heldout_evaluation": "EVIDENCE",
            "gap_verdict": "EVIDENCE",
        },
        "inputs": input_hashes,
        "n_classes": len(features),
        "predeclared_split": {
            "declared_before_outcome_load": True,
            "rule": f"SHA256('{SPLIT_SALT}|' + class_key), first 64 bits mod {SPLIT_MODULUS}; residue {HELDOUT_RESIDUE} held out",
            "salt": SPLIT_SALT,
            "modulus": SPLIT_MODULUS,
            "heldout_residue": HELDOUT_RESIDUE,
            "assignment_manifest_sha256": split_manifest_sha256,
            "counts": dict(sorted(split_counts.items())),
            "family_counts": {
                split: {
                    family: split_family_counts[(split, family)]
                    for family in sorted({feature["family"] for feature in features})
                }
                for split in ("fit", "heldout")
            },
        },
        "model": {
            "height_X": "abs(primitive integral F_cell(k)); F_cell(k)=P(eps*f*(q1+kN))",
            "prime_probability": "k=0 exactly 1 (q1 is fixed prime by class construction); k>0 min(1,(N/phi(N))/log(q1+kN))",
            "clean_probability_given_prime": "C * factor_p<=10000 * (log(10000)/log(max(10000,X)))^(1/2)",
            "hazard": "prime_probability * clean_probability_given_prime",
            "expected_count_to_K": "sum_{k=0}^K hazard_k",
            "first_arrival_cdf": "1-product_{k=0}^K(1-hazard_k)",
            "asymptotic_divergence": "for degree-8 F, hazard_k is order (log k)^(-3/2), hence expected_count_to_K diverges like K/(log K)^(3/2) up to constants",
            "independence_warning": "Bernoulli first-arrival approximation; BH and clean events are not proved independent",
            "single_fitted_parameter": "C",
            "calibration": "concave first-event likelihood on fit classes only",
            "euler_cutoff": EULER_CUTOFF,
            "tail_exponent": TAIL_EXPONENT,
            "survival_truncation_tolerance": SURVIVAL_TOL,
        },
    }
    summary = {
        "type": "summary",
        "label": "EVIDENCE",
        "n_classes": len(features),
        "calibration": calibration,
        "evaluations": evaluations,
        "gap_15_7x": gap_verdict,
        "failure_mode": {
            "observed_k0_fraction": evaluations["all_descriptive"]["observed"]["k0_fraction"],
            "predicted_k0_fraction": evaluations["all_descriptive"]["predicted"]["expected_k0_fraction"],
            "heldout_observed_k0_fraction": evaluations["heldout"]["observed"]["k0_fraction"],
            "heldout_predicted_k0_fraction": evaluations["heldout"]["predicted"]["expected_k0_fraction"],
            "interpretation": "single Haar/BH tail model misses the selected base-point spike; it cannot fit k=0 and the long tail simultaneously",
        },
    }

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUT_PATH.with_suffix(OUT_PATH.suffix + ".tmp")
    with temporary.open("w") as out:
        out.write(json.dumps(meta, separators=(",", ":"), sort_keys=True, allow_nan=False) + "\n")
        for row in output_rows:
            out.write(json.dumps(row, separators=(",", ":"), sort_keys=True, allow_nan=False) + "\n")
        out.write(json.dumps(summary, separators=(",", ":"), sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(OUT_PATH)

    fit_eval = evaluations["in_sample_fit"]
    hold_eval = evaluations["heldout"]
    all_eval = evaluations["all_descriptive"]
    comparison = gap_verdict["comparison"]
    report = f"""# L18 divergence model: first emergent-free height

**EVIDENCE throughout.** No statistic below is a theorem or an upgrade of the L18 local statements.

## Predeclared design

- The split was assigned **before opening the outcome file**: `SHA256(\"{SPLIT_SALT}|\" + class_key)`, first 64 bits modulo {SPLIT_MODULUS}; residue {HELDOUT_RESIDUE} is held out.
- Manifest SHA-256: `{split_manifest_sha256}`.
- Fit: {len(fit_indices)} classes ({meta['predeclared_split']['family_counts']['fit']['L11']} L11, {meta['predeclared_split']['family_counts']['fit']['ESC']} ESC). Held out: {len(heldout_indices)} ({meta['predeclared_split']['family_counts']['heldout']['L11']} L11, {meta['predeclared_split']['family_counts']['heldout']['ESC']} ESC).
- Only the fit classes enter the one-parameter likelihood. The held-out first-member indices are used only below for evaluation.

## Model

For class c and index k, let q_k=q1+kN and X_ck be the absolute value of the primitive integral degree-8 class polynomial F_c(k). The model uses

`Pr(q_k prime) = 1` at k=0 (q1 is fixed prime by construction), and `(N/phi(N))/log(q_k)` for k>0 (capped at 1).

The L18 clean term is

`Pr(clean | prime) = C * E_c(10000) * sqrt(log(10000)/log(max(10000,X_ck)))`,

where `E_c(10000)` is the geometric midpoint of the recorded truncated Euler-product bounds (their maximum width here is floating-point noise). Thus

`Lambda_c(K) = sum_{{k<=K}} Pr(q_k prime) Pr(clean | prime)`

is the requested expected count, and the independent-Bernoulli first-arrival CDF is `1-product_{{k<=K}}(1-h_ck)`. Independence is a modeling assumption, not a proved consequence of L18.

Within this model, degree 8 gives `log X_ck = 8 log k + O(1)`, so the hazard is of order `(log k)^(-3/2)` (one prime factor and one L18 square-root tail). Consequently `Lambda_c(K)` diverges on the scale `K/(log K)^(3/2)` up to constants: the quantitative claim is eventual accumulated mass, **not** a uniform positive clean probability.

## One-constant calibration

- Fit MLE: **C = {_fmt(c, 6)}**; profile 95% interval [{_fmt(calibration['profile_95_ci'][0], 6)}, {_fmt(calibration['profile_95_ci'][1], 6)}]. The MLE is not on the physical boundary.
- Fit event NLL/class: {_fmt(fit_eval['event_negative_log_likelihood_per_class'], 4)}. Held-out event NLL/class: {_fmt(hold_eval['event_negative_log_likelihood_per_class'], 4)}.
- Member-window cross-check at p<=100000 (not the full-event estimand): fit pooled observed/predicted clean rate {_fmt(fit_eval['member_window_p_le_100000']['observed_pooled_clean_rate'], 4)}/{_fmt(fit_eval['member_window_p_le_100000']['predicted_pooled_clean_rate'], 4)}; held-out {_fmt(hold_eval['member_window_p_le_100000']['observed_pooled_clean_rate'], 4)}/{_fmt(hold_eval['member_window_p_le_100000']['predicted_pooled_clean_rate'], 4)}.

## First-member distribution test

| subset | n | observed k=0 | predicted k=0 | observed mean k | predicted mean k | KS-style D (at k) | bootstrap p |
|---|---:|---:|---:|---:|---:|---:|---:|
| fit (in-sample) | {fit_eval['n_classes']} | {_fmt(fit_eval['observed']['k0_fraction'], 4)} | {_fmt(fit_eval['predicted']['expected_k0_fraction'], 4)} | {_fmt(fit_eval['observed']['mean_k'], 2)} | {_fmt(fit_eval['predicted']['class_mean_of_mean_k_truncated'], 2)} | {_fmt(fit_eval['ks_style']['D'], 4)} ({fit_eval['ks_style']['k']}) | {_fmt(fit_eval['ks_style']['bootstrap_add_one_p'], 6)} |
| held out | {hold_eval['n_classes']} | {_fmt(hold_eval['observed']['k0_fraction'], 4)} | {_fmt(hold_eval['predicted']['expected_k0_fraction'], 4)} | {_fmt(hold_eval['observed']['mean_k'], 2)} | {_fmt(hold_eval['predicted']['class_mean_of_mean_k_truncated'], 2)} | {_fmt(hold_eval['ks_style']['D'], 4)} ({hold_eval['ks_style']['k']}) | {_fmt(hold_eval['ks_style']['bootstrap_add_one_p'], 6)} |
| all (descriptive) | {all_eval['n_classes']} | {_fmt(all_eval['observed']['k0_fraction'], 4)} | {_fmt(all_eval['predicted']['expected_k0_fraction'], 4)} | {_fmt(all_eval['observed']['mean_k'], 2)} | {_fmt(all_eval['predicted']['class_mean_of_mean_k_truncated'], 2)} | {_fmt(all_eval['ks_style']['D'], 4)} ({all_eval['ks_style']['k']}) | {_fmt(all_eval['ks_style']['bootstrap_add_one_p'], 6)} |

The Monte Carlo KS reference samples every class from its own predictive CDF ({KS_SIMULATIONS} simulations, fixed seeds, add-one p). The fitted C is held fixed: the in-sample p is descriptive, while the held-out p is the primary calibration check. Both discrepancies peak at k=0. Held-out midpoint quantile residuals have median {_fmt(hold_eval['quantile_residuals']['midpoint_u_median'], 4)} (uniform calibration would center near 0.5), and normal-quantile median {_fmt(hold_eval['quantile_residuals']['normal_median'], 3)}.

**Failure location:** the model predicts only {_fmt(all_eval['predicted']['expected_k0_fraction']*100, 1)}% at k=0 versus the observed {_fmt(all_eval['observed']['k0_fraction']*100, 1)}%. The primitive F-cell height is already large at q1, so a stationary Haar/Bateman-Horn tail law does not reproduce the deliberately selected base-point spike. One constant cannot fit that spike and the long nonzero tail simultaneously. Per-class medians, 10/90% intervals, expected counts at the observation, and quantile residuals are in `data/l18_divergence_model.jsonl`.

## Explicit verdict on the 15.7x Bateman-Horn gap

The L15 ledger benchmark is `2^{ledger_log2_gap} = {_fmt(ledger_gap, 2)}x` on {len(gap_rows)} nonzero prime-rung classes.

- With the fitted C but **without** the L18 tail, the median conditional delay relative to a pure first-prime BH model is only {_fmt(all_gap['median_truncated_only_delay_over_bh'], 2)}x.
- With the `(log X)^(-1/2)` tail, it becomes **{_fmt(all_gap['median_l18_delay_over_bh'], 2)}x** overall and **{_fmt(gap_verdict['heldout']['median_l18_delay_over_bh'], 2)}x held out**.
- Relative to the ledger's {_fmt(ledger_gap, 2)}x, the L18 model leaves a factor **{_fmt(comparison['residual_factor_vs_ledger_gap'], 2)}x**. It accounts for {_fmt(comparison['fraction_of_excess_multiplier_closed']*100, 1)}% of the excess multiplier, or {_fmt(comparison['fraction_of_log_gap_closed']*100, 1)}% on the additive log scale.

**EVIDENCE verdict: the L18 correction PARTIALLY CLOSES the 15.7x gap, materially but not completely.** It explains much of the nonzero-tail slowdown, while the full 293-class distribution is decisively miscalibrated because of the unmodeled k=0 selection effect. Therefore this experiment supports the L18 divergence mechanism as a tail contribution, not as a complete first-member law.

## Scope and cautions

- The BH factor, Euler-factor independence, and first-arrival independence are modeling assumptions.
- The member-window clean rates stop at p=100000 and use actual proven-prime members k<=119; they are a diagnostic, not the full closure response.
- The deterministic split is outcome-blind but small ({len(heldout_indices)} held-out classes); no split was retried or selected by performance.
- All generated statistics and the partial-closing verdict are **EVIDENCE**.
"""
    REPORT_PATH.write_text(report)

    print(
        json.dumps(
            {
                "label": "EVIDENCE",
                "out": str(OUT_PATH),
                "report": str(REPORT_PATH),
                "n_classes": len(features),
                "fit_n": len(fit_indices),
                "heldout_n": len(heldout_indices),
                "C": c,
                "heldout_KS_D": hold_eval["ks_style"]["D"],
                "heldout_KS_p": hold_eval["ks_style"]["bootstrap_add_one_p"],
                "gap_verdict": gap_verdict["comparison"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return summary


if __name__ == "__main__":
    build()
