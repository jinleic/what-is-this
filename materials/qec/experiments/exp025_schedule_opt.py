"""EXP-025 v2: fixed-depth translation-invariant schedule optimisation.

V1 was retired as development-only before any PBB result.  V2 is prospectively
frozen here.  The default invocation prepares and audits the schedule population
only; production is executable-gated until the imported supervised baseline
runner passes the EXP-025 resume smoke.  Partial work is never canonicalized.
"""
from __future__ import annotations

# --------------------------- frozen protocol v2 ---------------------------
PROTOCOL_VERSION = "exp025-preregistered-2026-08-13-v2"
SCHEDULE_CLASS = "translation-invariant"
BASE_SEED = 20_260_813
P = 0.002
ROUNDS = 12
STAGE1_K_UNIFORM = 16
STAGE1_SHOTS = 1_000
STAGE2_TOP_K = 3
STAGE2_SHOTS = 10_000
MC_WORKERS_MAX = 6
DEADLINE_S = 30.0
CHECKPOINT_EVERY = 25
ALPHA = 0.05
GROSS_DEPTH = 7
PBB_DEPTH = 8
EXPECTED_GROSS_SCHEDULES = 8_496
EXPECTED_PBB_SCHEDULES = 9_968
PBB_CODE_ID = "12_6_0193"
EXP016_SCHEDULES_PER_CODE = 5
ENUMERATION_LIMIT = 200_000
ENUMERATION_TIME_LIMIT_S = 600.0
NOISELESS_VALIDATION_SHOTS = 4_000

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
from typing import Any

for _thread_var in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_var] = "1"

import numpy as np  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.artifacts import canonical_route  # noqa: E402
from qec_research.circuits.bb_syndrome import DEPTH7_SCHEDULE, _directional_neighbors  # noqa: E402
from qec_research.circuits.bicycle_schedule import (  # noqa: E402
    bb_supports_and_orbits, pbb_supports_and_orbits,
)
from qec_research.circuits.scheduling import (  # noqa: E402
    enumerate_orbit_schedules, verify_schedule,
)
from qec_research.codes.bicycle import (  # noqa: E402
    BRAVYI_BB, PBBSpec, build_bb,
)
from qec_research.decoders.bposd_dem import clopper_pearson  # noqa: E402

CATALOGUE_PATH = ROOT / "third_party/qcode-discovery/results/campaign7_publication_merged.jsonl"
EXP016_PATH = ROOT / "results/raw/exp016_schedule_controlled.json"
V1_AMENDMENT_PATH = ROOT / "results/partial_runs/exp025_v1_feasibility_amendment.json"
PARTIAL_CATALOG_PATH = ROOT / "results/partial_runs/exp025_schedule_catalog-v2.json"
PARTIAL_STATE_PATH = ROOT / "results/partial_runs/exp025_schedule_opt-v2-partial.json"
HARNESS_GATE_PATH = ROOT / "results/partial_runs/exp025_supervised_harness_gate-v2.json"
QUARANTINE_PATH = ROOT / "results/quarantine/exp025_schedule_opt-v2-failing.json"
RAW_RESULT_PATH = ROOT / "results/raw/exp025_schedule_opt.json"
BEST_SCHEDULES_PATH = ROOT / "results/raw/exp025_best_schedules.json"
PROCESSED_RESULT_PATH = ROOT / "results/processed/exp025_schedule_opt.json"
REPORT_PATH = ROOT / "notes/agent_reports/exp025_schedule_opt.md"
EXP026_PATH = ROOT / "experiments/exp026_decoder_codesign.py"
CODE_KEYS = ("gross", "pbb")


def protocol() -> dict[str, Any]:
    return {
        "version": 2,
        "protocol_version": PROTOCOL_VERSION,
        "revision": {
            "reason": "v1 infeasible before any PBB production row",
            "v1_data_reused": False,
            "v1_feasibility_amendment": str(V1_AMENDMENT_PATH.relative_to(ROOT)),
        },
        "schedule_class": SCHEDULE_CLASS,
        "schedule_populations": {
            "gross": {
                "count": EXPECTED_GROSS_SCHEDULES,
                "exhaustive_scope": "translation-invariant only",
                "depth": GROSS_DEPTH,
                "depth_claim": "exact(OPTIMAL) within translation-invariant class",
                "unrestricted_depth_note": (
                    "unrestricted depth-6 exclusion is external ASC evidence, not proved "
                    "by this repository enumeration"
                ),
            },
            "pbb": {
                "count": EXPECTED_PBB_SCHEDULES,
                "exhaustive_scope": "translation-invariant only",
                "depth": PBB_DEPTH,
                "depth_claim": (
                    "class-free exact: lower bound max-check-weight=8 plus valid "
                    "depth-8 witness; candidate population remains translation-invariant"
                ),
            },
        },
        "sampling": {
            "method": "fresh uniform without replacement from each complete TI population",
            "K_per_code": STAGE1_K_UNIFORM,
            "rng": "numpy.random.default_rng",
            "seed": BASE_SEED,
            "extras": (
                "five persisted EXP-016 schedules per code and published Gross default; "
                "exact duplicate hashes run once"
            ),
        },
        "stage1": {
            "shots_per_schedule": STAGE1_SHOTS,
            "selection": (
                "three schedules per code by ascending operational-failure count then "
                "schedule hash; published Gross default is not excluded"
            ),
        },
        "stage2": {
            "shots_per_selected_schedule": STAGE2_SHOTS,
            "selected_per_code": STAGE2_TOP_K,
            "gross_default_added_if_absent": True,
            "fresh_heldout_seeds": True,
            "early_stopping": False,
        },
        "circuit": {
            "p": P, "rounds": ROUNDS, "basis": "Z",
            "noise": (
                "DEPOLARIZE1(p) at every data/ancilla idle slot; DEPOLARIZE2(p) "
                "after each controlled-P gate; prep/measurement flip p"
            ),
        },
        "decoder": {
            "candidate": "EXP-026 v2 CANDIDATES[0] baseline_bposd",
            "bp": "min-sum max_iter=30 scale=0.625",
            "osd": "osd0 order=0",
            "undecomposed_dem": True,
            "deadline_s_per_shot": DEADLINE_S,
            "checkpoint_every_shots_max": CHECKPOINT_EVERY,
            "operational_failure": (
                "logical mismatch OR invalid syndrome correction OR decoder timeout"
            ),
            "persistent_supervised_worker": True,
        },
        "statistics": {
            "intervals": "two-sided equal-tailed Clopper-Pearson 95%",
            "cross_circuit_pairing": False,
            "test": "independent one-sided Fisher exact",
            "multiplicity": "Holm across the three screened candidate comparisons",
            "gross_default_question": "separate heldout independent exact/Holm family",
            "winner_question": "best heldout PBB versus best heldout Gross, independent Fisher",
            "alpha": ALPHA,
        },
        "verdict": {
            "BREAKTHROUGH_CANDIDATE": (
                "best heldout PBB significantly lower operational-failure rate than "
                "best heldout Gross with no resource relaxation"
            ),
            "NEGATIVE": "best heldout PBB significantly higher than best heldout Gross",
            "INCONCLUSIVE": "neither directional independent exact test is significant",
        },
        "workers_max": MC_WORKERS_MAX,
        "thread_limits": {name: 1 for name in (
            "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
        )},
    }


def protocol_hash() -> str:
    return hashlib.sha256(
        json.dumps(protocol(), sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def slot_hash(slot: dict[tuple[int, int], int]) -> str:
    payload = ",".join(
        f"{check}:{qubit}:{slot_index}"
        for (check, qubit), slot_index in sorted(slot.items())
    ).encode()
    return hashlib.sha256(payload).hexdigest()[:16]


def slot_to_json(slot: dict[tuple[int, int], int]) -> list[list[int]]:
    return [[int(c), int(q), int(t)] for (c, q), t in sorted(slot.items())]


def slot_from_json(rows: list[list[int]]) -> dict[tuple[int, int], int]:
    return {(int(c), int(q)): int(t) for c, q, t in rows}


def seed_for(stage: str, code_key: str, schedule_hash_value: str) -> int:
    payload = f"{BASE_SEED}|{stage}|{code_key}|{schedule_hash_value}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little") % (2**63 - 2) + 1


def load_pbb_spec() -> tuple[PBBSpec, dict[str, Any]]:
    rows = [json.loads(line) for line in CATALOGUE_PATH.read_text().splitlines()]
    matches = [row for row in rows if row.get("code_id") == PBB_CODE_ID]
    if len(matches) != 1:
        raise AssertionError(f"expected one {PBB_CODE_ID} row, found {len(matches)}")
    row = matches[0]
    if (row["n"], row["k"], row["d"]) != (144, 12, 12):
        raise AssertionError("PBB target parameters changed")
    if not row.get("d_is_exact") or not row.get("milp_exact_deep"):
        raise AssertionError("PBB distance-12 exactness metadata missing")
    spec = PBBSpec(
        row["ell"], row["m"],
        [tuple(term) for term in row["A_terms"]],
        [tuple(term) for term in row["B_terms"]],
        [tuple(term) for term in row["C_terms"]],
        [tuple(term) for term in row["D_terms"]], name=PBB_CODE_ID,
    )
    return spec, row


def code_context(code_key: str) -> dict[str, Any]:
    if code_key == "gross":
        code, supports, orbits = bb_supports_and_orbits(BRAVYI_BB["[[144,12,12]]"])
        return {
            "key": code_key, "code": code, "supports": supports, "orbits": orbits,
            "depth": GROSS_DEPTH, "expected_count": EXPECTED_GROSS_SCHEDULES,
            "label": "CSS Gross [[144,12,12]]", "catalogue": None,
        }
    if code_key == "pbb":
        spec, row = load_pbb_spec()
        code, supports, orbits = pbb_supports_and_orbits(spec)
        return {
            "key": code_key, "code": code, "supports": supports, "orbits": orbits,
            "depth": PBB_DEPTH, "expected_count": EXPECTED_PBB_SCHEDULES,
            "label": f"non-CSS PBB [[144,12,12]] {PBB_CODE_ID}", "catalogue": row,
        }
    raise ValueError(code_key)


def published_gross_slot(context: dict[str, Any]) -> dict[tuple[int, int], int]:
    hx, hz = build_bb(BRAVYI_BB["[[144,12,12]]"])
    x_neighbors, z_neighbors = _directional_neighbors(hx, hz)
    dim = hx.shape[0]
    orbit_by_edge = context["orbits"]
    x_orbits = [orbit_by_edge[(0, int(x_neighbors[d, 0]))] for d in range(6)]
    z_orbits = [orbit_by_edge[(dim, int(z_neighbors[d, 0]))] for d in range(6)]
    orbit_slot: dict[int, int] = {}
    for time_slot, (x_direction, z_direction) in enumerate(DEPTH7_SCHEDULE):
        if x_direction != "idle":
            orbit_slot[x_orbits[int(x_direction)]] = time_slot
        if z_direction != "idle":
            orbit_slot[z_orbits[int(z_direction)]] = time_slot
    slot = {edge: orbit_slot[orbit] for edge, orbit in orbit_by_edge.items()}
    if not verify_schedule(context["supports"], slot)["valid"]:
        raise AssertionError("published Gross schedule failed independent verification")
    return slot


def exp016_runs(code_key: str) -> list[dict[str, Any]]:
    prior = json.loads(EXP016_PATH.read_text())
    needle = "CSS-BB" if code_key == "gross" else PBB_CODE_ID
    entries = [entry for entry in prior["codes"] if needle in entry["label"]]
    if len(entries) != 1 or len(entries[0]["runs"]) != EXP016_SCHEDULES_PER_CODE:
        raise AssertionError(f"EXP-016 {code_key} persistence contract changed")
    return entries[0]["runs"]


def frozen_sample_indices() -> dict[str, list[int]]:
    """Draw both fresh samples in one fixed, documented RNG stream."""
    rng = np.random.default_rng(BASE_SEED)
    return {
        "gross": [int(i) for i in rng.choice(
            EXPECTED_GROSS_SCHEDULES, STAGE1_K_UNIFORM, replace=False
        )],
        "pbb": [int(i) for i in rng.choice(
            EXPECTED_PBB_SCHEDULES, STAGE1_K_UNIFORM, replace=False
        )],
    }


def build_schedule_catalog() -> dict[str, Any]:
    indices = frozen_sample_indices()
    result: dict[str, Any] = {
        "experiment": "EXP-025", "protocol": protocol(),
        "protocol_sha256": protocol_hash(), "schedule_class": SCHEDULE_CLASS,
        "uniform_indices": indices, "codes": {}, "verdict": "PREPARED_PARTIAL",
    }
    for code_key in CODE_KEYS:
        context = code_context(code_key)
        started = time.perf_counter()
        schedules, complete = enumerate_orbit_schedules(
            context["supports"], context["code"].n, context["depth"], context["orbits"],
            limit=ENUMERATION_LIMIT, time_limit_s=ENUMERATION_TIME_LIMIT_S,
        )
        if not complete or len(schedules) != context["expected_count"]:
            raise RuntimeError(
                f"{code_key} TI enumeration incomplete: complete={complete}, "
                f"count={len(schedules)}"
            )
        by_hash = {slot_hash(slot): (index, slot) for index, slot in enumerate(schedules)}
        if len(by_hash) != len(schedules):
            raise AssertionError("schedule hash collision in complete TI population")
        chosen: dict[str, dict[str, Any]] = {}

        def add(slot: dict[tuple[int, int], int], origin: str, index: int) -> None:
            if not verify_schedule(context["supports"], slot)["valid"]:
                raise AssertionError(f"invalid {code_key} candidate")
            schedule_hash_value = slot_hash(slot)
            item = chosen.setdefault(schedule_hash_value, {
                "schedule_hash": schedule_hash_value, "enumeration_index": int(index),
                "depth": context["depth"], "schedule_class": SCHEDULE_CLASS,
                "origins": [], "slot_map": slot_to_json(slot),
            })
            if item["slot_map"] != slot_to_json(slot):
                raise AssertionError("truncated hash collision")
            if origin not in item["origins"]:
                item["origins"].append(origin)

        for index in indices[code_key]:
            add(schedules[index], "fresh_uniform_K16_seed_20260813", index)
        for prior in exp016_runs(code_key):
            slot = slot_from_json(prior["slot_map"])
            if slot_hash(slot) != prior["schedule_hash"] or prior["schedule_hash"] not in by_hash:
                raise AssertionError("EXP-016 candidate not recovered in complete TI population")
            current_index, current_slot = by_hash[prior["schedule_hash"]]
            if slot != current_slot:
                raise AssertionError("EXP-016 slot differs from regenerated TI schedule")
            add(slot, "exp016_persisted_sample", current_index)
        default_hash = None
        if code_key == "gross":
            default = published_gross_slot(context)
            default_hash = slot_hash(default)
            index, regenerated = by_hash[default_hash]
            if default != regenerated:
                raise AssertionError("published Gross schedule regeneration mismatch")
            add(default, "published_gross_default", index)
        result["codes"][code_key] = {
            "label": context["label"], "depth": context["depth"],
            "schedule_class": SCHEDULE_CLASS,
            "population_count": len(schedules),
            "population_exhaustive_scope": "translation-invariant only",
            "enumeration_status": "OPTIMAL",
            "depth_claim": protocol()["schedule_populations"][code_key]["depth_claim"],
            "default_schedule_hash": default_hash,
            "candidate_count": len(chosen), "candidates": list(chosen.values()),
            "enumeration_wall_s": time.perf_counter() - started,
        }
    result["wall_s"] = sum(code["enumeration_wall_s"] for code in result["codes"].values())
    return result


def validate_catalog_scope(catalog: dict[str, Any]) -> None:
    if catalog.get("protocol_sha256") != protocol_hash():
        raise RuntimeError("catalog protocol mismatch")
    if catalog.get("schedule_class") != SCHEDULE_CLASS:
        raise RuntimeError("schedule class missing or changed")
    for code_key, expected in (
        ("gross", EXPECTED_GROSS_SCHEDULES), ("pbb", EXPECTED_PBB_SCHEDULES),
    ):
        code = catalog["codes"][code_key]
        if code.get("schedule_class") != SCHEDULE_CLASS:
            raise RuntimeError(f"{code_key} schedule class mismatch")
        if code.get("population_count") != expected:
            raise RuntimeError(f"{code_key} TI population count mismatch")
        if code.get("population_exhaustive_scope") != "translation-invariant only":
            raise RuntimeError(f"{code_key} exhaustiveness scope widened")
    gross_claim = catalog["codes"]["gross"]["depth_claim"]
    if "within translation-invariant class" not in gross_claim:
        raise RuntimeError("Gross depth claim improperly widened")
    pbb_claim = catalog["codes"]["pbb"]["depth_claim"]
    if "class-free exact" not in pbb_claim:
        raise RuntimeError("PBB class-free depth claim missing")


def select_stage2(stage1: dict[str, dict[str, Any]], catalog: dict[str, Any]) -> dict[str, list[str]]:
    """Apply the frozen failure-count/hash tie-break; no v1 data is accepted."""
    selection: dict[str, list[str]] = {}
    for code_key in CODE_KEYS:
        expected = {item["schedule_hash"] for item in catalog["codes"][code_key]["candidates"]}
        if set(stage1.get(code_key, {})) != expected:
            raise RuntimeError(f"incomplete/asymmetric stage1 coverage for {code_key}")
        rows = list(stage1[code_key].values())
        for row in rows:
            if row.get("protocol_sha256") != protocol_hash() or row.get("shots") != STAGE1_SHOTS:
                raise RuntimeError("stage1 row is not complete v2 data")
        rows.sort(key=lambda row: (row["operational_failures"], row["schedule_hash"]))
        selection[code_key] = [row["schedule_hash"] for row in rows[:STAGE2_TOP_K]]
    return selection


def fisher_greater(
    first_failures: int, first_shots: int, second_failures: int, second_shots: int,
) -> float:
    """P(first rate > second rate), using independent samples only."""
    return float(fisher_exact(
        [[first_failures, first_shots - first_failures],
         [second_failures, second_shots - second_failures]],
        alternative="greater",
    ).pvalue)


def independent_exact_test(
    first: dict[str, Any], second: dict[str, Any], *, alternative: str,
) -> dict[str, Any]:
    """Test two distinct circuits; paired/McNemar inputs are expressly rejected."""
    if first.get("circuit_sha256") == second.get("circuit_sha256"):
        raise ValueError("independent cross-circuit test requires distinct circuit hashes")
    if first.get("paired") or second.get("paired"):
        raise ValueError("McNemar/paired statistics are forbidden across distinct circuits")
    if alternative == "first_greater":
        pvalue = fisher_greater(
            first["operational_failures"], first["shots"],
            second["operational_failures"], second["shots"],
        )
    elif alternative == "first_less":
        pvalue = fisher_greater(
            second["operational_failures"], second["shots"],
            first["operational_failures"], first["shots"],
        )
    else:
        raise ValueError(alternative)
    return {
        "test": "one-sided Fisher exact on independent circuits",
        "alternative": alternative, "p_value": pvalue,
    }


def holm(pvalues: dict[str, float], alpha: float = ALPHA) -> dict[str, dict[str, Any]]:
    ordered = sorted(pvalues.items(), key=lambda item: (item[1], item[0]))
    result: dict[str, dict[str, Any]] = {}
    still_rejecting = True
    count = len(ordered)
    for index, (name, pvalue) in enumerate(ordered):
        threshold = alpha / (count - index)
        rejected = still_rejecting and pvalue <= threshold
        if not rejected:
            still_rejecting = False
        result[name] = {
            "p_value": pvalue, "holm_threshold": threshold,
            "holm_adjusted_p_value": min(1.0, pvalue * (count - index)),
            "rejected": rejected,
        }
    return result


def summarize_records(records: list[dict[str, Any]]) -> dict[str, Any]:
    shots = len(records)
    failures = sum(bool(row["operational_failure"]) for row in records)
    lo, hi = clopper_pearson(failures, shots) if shots else (0.0, 1.0)
    return {
        "shots": shots, "operational_failures": failures,
        "operational_failure_rate": failures / shots if shots else None,
        "clopper_pearson_95_ci": [lo, hi],
        "logical_mismatches": sum(bool(row["logical_mismatch"]) for row in records),
        "invalid_corrections": sum(bool(row["invalid_correction"]) for row in records),
        "timeouts": sum(bool(row["timeout"]) for row in records),
    }


def assert_complete_coverage(
    state: dict[str, Any], catalog: dict[str, Any], selection: dict[str, list[str]],
) -> None:
    """Canonicalization gate: every declared row must be complete and symmetric."""
    if state.get("protocol_sha256") != protocol_hash():
        raise RuntimeError("state protocol mismatch")
    for code_key in CODE_KEYS:
        expected_stage1 = {
            item["schedule_hash"] for item in catalog["codes"][code_key]["candidates"]
        }
        actual_stage1 = set(state.get("stage1", {}).get(code_key, {}))
        if actual_stage1 != expected_stage1:
            raise RuntimeError(f"incomplete/asymmetric stage1 coverage for {code_key}")
        expected_stage2 = set(selection[code_key])
        default_hash = catalog["codes"][code_key]["default_schedule_hash"]
        if code_key == "gross" and default_hash is not None:
            expected_stage2.add(default_hash)
        actual_stage2 = set(state.get("stage2", {}).get(code_key, {}))
        if actual_stage2 != expected_stage2:
            raise RuntimeError(f"incomplete/asymmetric stage2 coverage for {code_key}")
        for stage, expected_shots in (("stage1", STAGE1_SHOTS), ("stage2", STAGE2_SHOTS)):
            for row in state[stage][code_key].values():
                if row.get("shots") != expected_shots:
                    raise RuntimeError(f"partial {stage} row")
                if row.get("protocol_sha256") != protocol_hash():
                    raise RuntimeError(f"foreign {stage} row")
                if row.get("fatal_errors", 0):
                    raise RuntimeError(f"failing {stage} row")


def determine_verdict(best_pbb: dict[str, Any], best_gross: dict[str, Any]) -> dict[str, Any]:
    pbb_less = independent_exact_test(best_pbb, best_gross, alternative="first_less")
    pbb_greater = independent_exact_test(best_pbb, best_gross, alternative="first_greater")
    pbb_rate = best_pbb["operational_failures"] / best_pbb["shots"]
    gross_rate = best_gross["operational_failures"] / best_gross["shots"]
    if pbb_rate < gross_rate and pbb_less["p_value"] < ALPHA:
        verdict = "BREAKTHROUGH_CANDIDATE"
    elif pbb_rate > gross_rate and pbb_greater["p_value"] < ALPHA:
        verdict = "NEGATIVE"
    else:
        verdict = "INCONCLUSIVE"
    return {
        "verdict": verdict, "pbb_less_test": pbb_less,
        "pbb_greater_test": pbb_greater,
    }


def import_supervised_harness():
    spec = importlib.util.spec_from_file_location("exp026_for_exp025", EXP026_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load EXP-026 supervised harness")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    required = ("SupervisedDecoder", "evaluate_candidate", "CANDIDATES", "_pure_rows")
    if not all(hasattr(module, name) for name in required):
        raise RuntimeError("EXP-026 supervised harness interface incomplete")
    baseline = module.CANDIDATES[0]
    if baseline.get("id") != "baseline_bposd" or module.CHECKPOINT_EVERY > CHECKPOINT_EVERY:
        raise RuntimeError("EXP-026 baseline/checkpoint contract mismatch")
    return module


def harness_gate_status() -> dict[str, Any]:
    return json.loads(HARNESS_GATE_PATH.read_text()) if HARNESS_GATE_PATH.exists() else {
        "experiment": "EXP-025", "protocol_sha256": protocol_hash(),
        "passed": False,
        "reason": (
            "EXP-025 resume smoke has not yet proved exact-prefix resume, <=25-shot "
            "checkpoints, timeout restart, and invalid-correction accounting"
        ),
        "verdict": "PRODUCTION_BLOCKED",
    }


def require_production_gate() -> None:
    gate = harness_gate_status()
    if gate.get("protocol_sha256") != protocol_hash() or not gate.get("passed"):
        raise RuntimeError(
            "EXP-025 v2 production is blocked until --harness-smoke writes a passing "
            "protocol-bound gate"
        )


def prepare() -> dict[str, Any]:
    catalog = build_schedule_catalog()
    validate_catalog_scope(catalog)
    write_json(PARTIAL_CATALOG_PATH, catalog)
    if not PARTIAL_STATE_PATH.exists():
        write_json(PARTIAL_STATE_PATH, {
            "experiment": "EXP-025", "protocol": protocol(),
            "protocol_sha256": protocol_hash(), "schedule_class": SCHEDULE_CLASS,
            "stage1": {"gross": {}, "pbb": {}},
            "selection": None, "stage2": {"gross": {}, "pbb": {}},
            "coverage_complete": False, "verdict": "PREPARED_PARTIAL",
        })
    return catalog


def run_harness_smoke() -> dict[str, Any]:
    """Focused executable gate only; never creates a scientific result."""
    module = import_supervised_harness()
    # Pure contract checks precede any expensive real-circuit smoke.  The EXP-026
    # focused suite proves timeout restart; EXP-025 additionally requires exact
    # prefix resume with a protocol-bound checkpoint before production.
    records = [
        {"shot_id": index, "operational_failure": False, "logical_mismatch": False,
         "invalid_correction": False, "timeout": False, "fatal": None,
         "latency_s": 0.0, "worker_restarted": False}
        for index in range(CHECKPOINT_EVERY)
    ]
    module.validate_record_prefix(records, CHECKPOINT_EVERY + 1)
    if module.validate_record_prefix(records, CHECKPOINT_EVERY + 1) != CHECKPOINT_EVERY:
        raise AssertionError("EXP-026 exact-prefix resume contract failed")
    gate = {
        "experiment": "EXP-025", "protocol_sha256": protocol_hash(),
        "interface": "experiments.exp026_decoder_codesign",
        "baseline_candidate": module.CANDIDATES[0],
        "checkpoint_every_shots": module.CHECKPOINT_EVERY,
        "deadline_s": module.DEADLINE_S,
        "exact_prefix_resume_contract": True,
        "timeout_restart_verified_by": "tests/test_exp026_decoder_codesign.py",
        "real_circuit_smoke": False,
        "passed": False,
        "reason": (
            "contract import and exact-prefix validation passed, but a real-circuit "
            "interrupted/resumed baseline decode has not been executed by this scaffold"
        ),
        "verdict": "PRODUCTION_BLOCKED",
    }
    write_json(HARNESS_GATE_PATH, gate)
    return gate


def canonicalize(state: dict[str, Any], catalog: dict[str, Any]) -> Path:
    """Route complete clean output canonically; reject every partial/failing state."""
    validate_catalog_scope(catalog)
    selection = state.get("selection")
    if selection is None:
        raise RuntimeError("stage2 selection missing")
    assert_complete_coverage(state, catalog, selection)
    route = canonical_route(clean=True, full_coverage=True)
    if route != "canonical":
        raise AssertionError(route)
    write_json(PROCESSED_RESULT_PATH, state)
    return PROCESSED_RESULT_PATH


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phase", choices=("prepare", "harness-smoke", "production"), default="prepare"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.phase == "prepare":
        catalog = prepare()
        print(
            f"prepared EXP-025 v2 protocol={protocol_hash()} "
            f"Gross={catalog['codes']['gross']['candidate_count']} "
            f"PBB={catalog['codes']['pbb']['candidate_count']} candidates; "
            f"wrote {PARTIAL_CATALOG_PATH.relative_to(ROOT)}"
        )
        return
    if args.phase == "harness-smoke":
        gate = run_harness_smoke()
        print(json.dumps(gate, indent=2))
        return
    require_production_gate()
    raise RuntimeError(
        "production gate passed but the production driver is intentionally not invoked "
        "by the no-production scaffold command"
    )


if __name__ == "__main__":
    main()
