"""Regression test for the EXP-031 verdict polarity bug (FR-013).

The mismatch branch once returned "POSITIVE" while its reason described a
counterexample, making the machine-readable verdict field ambiguous between
"conjecture holds" and "conjecture falsified".  A counterexample must yield
the negative verdict.
"""
import importlib.util
import sys
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "exp031_depth_criterion",
    Path(__file__).resolve().parents[1] / "experiments" / "exp031_depth_criterion.py",
)
_MOD = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _MOD  # dataclasses require the module to be registered
_SPEC.loader.exec_module(_MOD)


def test_mismatch_is_negative_never_positive():
    verdict, reason = _MOD.decide_verdict(n_mismatches=1, n_decided=208, n_undecided=3)
    assert verdict == "NEGATIVE"
    assert "FALSIFIED" in reason
    # The historical bug: POSITIVE on the counterexample path.
    assert verdict != "POSITIVE"


def test_zero_mismatch_large_sample_is_positive():
    verdict, reason = _MOD.decide_verdict(n_mismatches=0, n_decided=208, n_undecided=3)
    assert verdict == "POSITIVE"
    assert "theorem candidate" in reason


def test_small_sample_is_inconclusive():
    verdict, _ = _MOD.decide_verdict(n_mismatches=0, n_decided=42, n_undecided=169)
    assert verdict == "INCONCLUSIVE"


def test_mismatch_dominates_small_sample():
    # Even with too few decided instances, a real counterexample decides the verdict.
    verdict, _ = _MOD.decide_verdict(n_mismatches=2, n_decided=42, n_undecided=169)
    assert verdict == "NEGATIVE"


def test_certified_infeasible_at_ignores_unknown():
    trace = [{"depth": 9, "status": "INFEASIBLE"}, {"depth": 10, "status": "UNKNOWN"}]
    assert _MOD.certified_infeasible_at(trace, 9)
    assert not _MOD.certified_infeasible_at(trace, 10)   # timeout is not a certificate
    assert not _MOD.certified_infeasible_at(trace, 11)


def test_refuted_prediction_is_mismatch_even_when_undecided():
    """FR-014: INFEASIBLE at the predicted depth refutes the prediction."""
    trace = [{"depth": 9, "status": "INFEASIBLE"},
             {"depth": 10, "status": "INFEASIBLE"},
             {"depth": 11, "status": "INFEASIBLE"}]
    c = _MOD.classify_prediction(decided=False, certified_depth=None,
                                 predicted_depth=10, w_max=9, trace=trace)
    assert c["matches_predicted_depth"] is False          # not None
    assert c["prediction_refuted_in_class"] is True
    assert c["two_value_law_holds"] is False              # 9 and 10 both infeasible
    assert c["two_value_refuted_in_class"] is True


def test_unknown_at_predicted_depth_stays_undecided():
    trace = [{"depth": 9, "status": "INFEASIBLE"}, {"depth": 10, "status": "UNKNOWN"}]
    c = _MOD.classify_prediction(decided=False, certified_depth=None,
                                 predicted_depth=10, w_max=9, trace=trace)
    assert c["matches_predicted_depth"] is None
    assert c["two_value_law_holds"] is None


def test_decided_paths_unchanged():
    trace = [{"depth": 9, "status": "INFEASIBLE"}, {"depth": 10, "status": "OPTIMAL"}]
    c = _MOD.classify_prediction(decided=True, certified_depth=10,
                                 predicted_depth=10, w_max=9, trace=trace)
    assert c["matches_predicted_depth"] is True
    assert c["two_value_law_holds"] is True
    assert c["prediction_refuted_in_class"] is False
