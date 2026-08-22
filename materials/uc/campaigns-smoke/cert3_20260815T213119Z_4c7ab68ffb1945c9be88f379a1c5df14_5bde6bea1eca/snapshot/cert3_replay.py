"""Independent mathematical replay for compact ``cert3`` proof traces."""

import argparse
import json
import math
import os
import sys
from collections import Counter
from fractions import Fraction

# Replay must execute exactly the hashed sources beside this file: refuse to
# materialize bytecode caches that a later audit could not distinguish.
sys.dont_write_bytecode = True

from flint import ctx
from mpmath import mp


HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# Import in the worker's order.  In particular, cert3 performs the final native
# precision setup after importing cert2 and diag_exhaust.
import cert2
import cert3
import arbcore
import bound_kkt
import diag_exhaust
import entropy


_NATIVE_PREC = 160
_WORK_PREC = 80
_DEFAULT_FACE_BOX_BUDGET = 100000
_REPLAY_TALLY_NAMES = (
    "processed",
    "infeasible",
    "corner",
    "ratio",
    "center",
    "center_mixed",
    "center_mixed_swap",
    "center_w",
    "face",
    "residual",
    "split",
)


def _verify_local_modules():
    expected = os.path.realpath(HERE)
    for module in (arbcore, bound_kkt, cert2, cert3, diag_exhaust, entropy):
        location = getattr(module, "__file__", None)
        if location is None:
            raise ImportError("local module %s has no source path" % module.__name__)
        actual = os.path.dirname(os.path.realpath(location))
        if actual != expected:
            raise ImportError(
                "local module %s resolved outside replay directory: %s" %
                (module.__name__, location))


_verify_local_modules()


def _runtime_error(event_number, event, message):
    return RuntimeError(
        "trace event %d (code %d): %s" % (event_number, event, message))


def _normalise_root_box(root_box):
    try:
        coordinates = tuple(root_box)
    except TypeError as exc:
        raise RuntimeError("root_box must contain five coordinate intervals") from exc
    if len(coordinates) != 5:
        raise RuntimeError("root_box must contain exactly five coordinate intervals")

    result = []
    for coordinate, pair in enumerate(coordinates):
        try:
            endpoints = tuple(pair)
        except TypeError as exc:
            raise RuntimeError(
                "root_box coordinate %d is not an interval" % coordinate) from exc
        if len(endpoints) != 2:
            raise RuntimeError(
                "root_box coordinate %d must have two endpoints" % coordinate)
        try:
            lo, hi = float(endpoints[0]), float(endpoints[1])
        except (TypeError, ValueError, OverflowError) as exc:
            raise RuntimeError(
                "root_box coordinate %d has non-numeric endpoints" % coordinate) from exc
        if not (math.isfinite(lo) and math.isfinite(hi)):
            raise RuntimeError(
                "root_box coordinate %d has non-finite endpoints" % coordinate)
        if not (0.0 <= lo <= hi <= 1.0):
            raise RuntimeError(
                "root_box coordinate %d is outside [0,1] or reversed" % coordinate)
        result.append((lo, hi))
    return tuple(result)


def _replay_parameters(t_text, parameters):
    if parameters is None:
        parameters = {}
    if not isinstance(parameters, dict):
        raise RuntimeError("parameters must be a dictionary when supplied")

    if "t_decimal" in parameters and parameters["t_decimal"] != t_text:
        raise RuntimeError("parameters t_decimal does not match replay target")

    work_prec = parameters.get("work_prec_bits", _WORK_PREC)
    if isinstance(work_prec, bool) or work_prec != _WORK_PREC:
        raise RuntimeError("replay requires work_prec_bits=80")

    face_box_budget = parameters.get(
        "face_box_budget", _DEFAULT_FACE_BOX_BUDGET)
    if (isinstance(face_box_budget, bool)
            or not isinstance(face_box_budget, int)
            or face_box_budget <= 0):
        raise RuntimeError("face_box_budget must be a positive integer")

    try:
        face_min_width = float(parameters.get(
            "face_min_width", cert3.DEFAULT_FACE_MIN_WIDTH))
    except (TypeError, ValueError, OverflowError) as exc:
        raise RuntimeError("face_min_width must be a positive finite number") from exc
    if not math.isfinite(face_min_width) or face_min_width <= 0.0:
        raise RuntimeError("face_min_width must be a positive finite number")

    return parameters, face_min_width, face_box_budget


def _worker_setup(t_text, parameters):
    # Match cert3_par.worker_main: form t and both global covers at the modules'
    # native 160-bit precision, then switch to the campaign's 80-bit precision.
    ctx.prec = _NATIVE_PREC
    mp.dps = 60
    try:
        t = mp.mpf(t_text)
    except (TypeError, ValueError) as exc:
        raise RuntimeError("t_text is not a valid decimal target") from exc
    gmax = arbcore.get_rh_gmax()
    rho_gmax = cert2.get_rho_gmax()
    ctx.prec = _WORK_PREC

    lambdas, _ = cert3.lambda_family(t)
    expected_lambdas = parameters.get("lambda_family")
    actual_lambdas = [repr(float(value)) for value in lambdas]
    if expected_lambdas is not None and expected_lambdas != actual_lambdas:
        raise RuntimeError(
            "lambda_family does not match values recomputed at 80-bit precision")

    # This is certify's setup after the worker has supplied its lambda family.
    t_arb = cert3._arb(t_text)
    t_upper = t_arb.upper()
    bound_kkt.RH_GMAX = gmax
    _, lamhat = cert3.lambda_family(t)
    center_lams = (0.0, float(lamhat))
    return t_arb, t_upper, gmax, rho_gmax, tuple(lambdas), center_lams


def _require_clear(value, event_number, event, rule):
    if value is None or not (value >= cert3.ZERO):
        raise _runtime_error(
            event_number, event, "%s did not certify a nonnegative bound" % rule)


def _validate_expected_tallies(counts, expected_tallies):
    if expected_tallies is None:
        return
    if not isinstance(expected_tallies, dict):
        raise RuntimeError("expected_tallies must be a dictionary when supplied")
    for name in _REPLAY_TALLY_NAMES:
        expected = expected_tallies.get(name, 0)
        if isinstance(expected, bool) or not isinstance(expected, int) or expected < 0:
            raise RuntimeError(
                "expected tally %s must be a nonnegative integer" % name)
        actual = counts[name]
        if actual != expected:
            raise RuntimeError(
                "trace tally mismatch for %s: expected %d, replayed %d" %
                (name, expected, actual))


def replay_trace(t_text, root_box, trace_path, parameters=None,
                 expected_tallies=None):
    """Replay and mathematically verify a ``cert3-trace-v1`` DFS trace.

    Only the discharge rule named by each terminal event is recomputed.  A split
    is intrinsically sound and therefore need not demonstrate that optional
    clearing rules failed first.
    """
    t_text = str(t_text)
    parameters, face_min_width, face_box_budget = _replay_parameters(
        t_text, parameters)
    root = _normalise_root_box(root_box)
    (t_arb, t_upper, gmax, rho_gmax,
     lambdas, center_lams) = _worker_setup(t_text, parameters)

    stack = [root]
    counts = Counter()
    try:
        trace_name = os.fspath(trace_path)
    except TypeError as exc:
        raise RuntimeError("trace_path must be a filesystem path") from exc

    with open(trace_name, "rb") as trace:
        while True:
            encoded = trace.read(1)
            if not encoded:
                break
            event_number = counts["processed"] + 1
            event = encoded[0]
            if event not in cert3.TRACE_TALLY_NAMES:
                raise _runtime_error(event_number, event, "unknown trace event")
            if not stack:
                raise _runtime_error(
                    event_number, event, "extra event after DFS stack was exhausted")

            raw = stack.pop()
            counts["processed"] += 1
            box = diag_exhaust.mean_contract(raw, t_upper)

            if event == cert3.TRACE_MEAN_INFEASIBLE:
                if box is not None:
                    raise _runtime_error(
                        event_number, event,
                        "mean contractor did not prove infeasibility")
                counts["infeasible"] += 1
                continue
            if box is None:
                raise _runtime_error(
                    event_number, event,
                    "mean contractor proved infeasibility but event claims another outcome")

            if event == cert3.TRACE_CORNER_INFEASIBLE:
                if diag_exhaust.phi_corner(box, t_arb) is not None:
                    raise _runtime_error(
                        event_number, event,
                        "phi_corner did not prove infeasibility")
                counts["infeasible"] += 1
                continue

            if event == cert3.TRACE_CORNER:
                _require_clear(
                    diag_exhaust.phi_corner(box, t_arb),
                    event_number, event, "corner rule")
            elif event == cert3.TRACE_RATIO:
                w_window = cert3.hull(
                    cert3._arb(repr(box[4][0])),
                    cert3._arb(repr(box[4][1])))
                if not cert2.ratio_rule(
                        box[:4], t_arb, W=w_window, rho_gmax=rho_gmax):
                    raise _runtime_error(
                        event_number, event, "ratio rule did not clear the box")
            elif event == cert3.TRACE_CENTER:
                _require_clear(
                    cert3.centered5_best(
                        box, t_arb, center_lams, gmax),
                    event_number, event, "centered5 rule")
            elif event == cert3.TRACE_CENTER_MIXED:
                _require_clear(
                    cert3.centered5_mixed(
                        box, t_arb, center_lams, gmax),
                    event_number, event, "centered5_mixed rule")
            elif event == cert3.TRACE_CENTER_MIXED_SWAP:
                _require_clear(
                    cert3.centered5_mixed_swap(
                        box, t_arb, center_lams, gmax),
                    event_number, event, "centered5_mixed_swap rule")
            elif event == cert3.TRACE_CENTER_W:
                _require_clear(
                    cert3.centered_w_best(
                        box, t_arb, center_lams, gmax),
                    event_number, event, "centered_w rule")
            elif event == cert3.TRACE_FACE:
                pinned_lams = []
                if box[3][0] >= cert3.Q2PIN:
                    pinned_lams = [
                        lam for lam in lambdas
                        if cert3.q2_pin(box, t_arb, lam, gmax)
                    ]
                if not pinned_lams:
                    raise _runtime_error(
                        event_number, event,
                        "q2 pin did not certify any lambda")
                face_box = (box[0], box[1], box[2], box[4])
                completed, _ = cert3.face_bb(
                    face_box, t_arb, tuple(pinned_lams), gmax,
                    min_width=face_min_width,
                    box_budget=face_box_budget,
                    deadline=None,
                )
                if not completed:
                    raise _runtime_error(
                        event_number, event,
                        "face_bb did not return COMPLETE within the bound budget")
            elif event == cert3.TRACE_RESIDUAL:
                raise _runtime_error(
                    event_number, event,
                    "residual leaves are invalid in a COMPLETE proof trace")
            elif cert3.TRACE_SPLIT_BASE <= event < cert3.TRACE_SPLIT_BASE + 5:
                coordinate = event - cert3.TRACE_SPLIT_BASE
                lo, hi = box[coordinate]
                left, right = cert3._split(box, coordinate)
                left_mid = left[coordinate][1]
                right_mid = right[coordinate][0]
                if (left_mid != right_mid
                        or not (lo < left_mid < hi)):
                    raise _runtime_error(
                        event_number, event,
                        "split midpoint is not strict and representable")
                stack.extend((left, right))
                counts["split"] += 1
                continue
            else:
                # TRACE_TALLY_NAMES and this dispatcher must stay in lockstep.
                raise _runtime_error(event_number, event, "unsupported trace event")

            counts[cert3.TRACE_TALLY_NAMES[event]] += 1

    if stack:
        raise RuntimeError(
            "trace ended after %d events with %d DFS boxes still pending" %
            (counts["processed"], len(stack)))

    replayed = {name: int(counts[name]) for name in _REPLAY_TALLY_NAMES}
    _validate_expected_tallies(replayed, expected_tallies)
    return replayed


def _cli_root(launch, result):
    index = result.get("slice")
    if isinstance(index, bool) or not isinstance(index, int):
        raise RuntimeError("result slice must be an integer")
    runs = launch.get("runs")
    if not isinstance(runs, list) or not 0 <= index < len(runs):
        raise RuntimeError("result slice is absent from launch")
    run = runs[index]
    if run.get("slice") != index or run.get("run_id") != result.get("run_id"):
        raise RuntimeError("result run does not match launch slice")
    try:
        w_lo = Fraction(run["w_lo"])
        w_hi = Fraction(run["w_hi"])
    except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
        raise RuntimeError("launch slice has invalid exact w endpoints") from exc
    return ((Fraction(0), Fraction(1)),) * 4 + ((w_lo, w_hi),)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Replay one cert3 binary proof trace")
    parser.add_argument("launch_json")
    parser.add_argument("result_json")
    parser.add_argument("trace_bin")
    args = parser.parse_args(argv)

    with open(args.launch_json) as handle:
        launch = json.load(handle)
    with open(args.result_json) as handle:
        result = json.load(handle)
    if launch.get("launch_sha256") != result.get("launch_sha256"):
        raise RuntimeError("launch and result hashes do not match")
    if launch.get("t_decimal") != result.get("t_decimal"):
        raise RuntimeError("launch and result targets do not match")

    replayed = replay_trace(
        launch["t_decimal"],
        _cli_root(launch, result),
        args.trace_bin,
        parameters=launch.get("parameters"),
        expected_tallies=result.get("tallies"),
    )
    print(json.dumps(replayed, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    main()
