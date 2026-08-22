"""Envelope check: every certified reversal is CSS-dominated at equal n.

Scope is *all* certified reversals, from both sources:

  * EXP-036's tie-branch closures over EXP-027's 87 undecided rows, re-derived
    from validated evidence and required to carry replay stamps hash-bound to
    freshly rebuilt CNFs; and
  * EXP-027's own two certified reversals (rows 331/339), which lie outside the
    undecided set and so are invisible to the EXP-036 pass.  Their PBB upper
    bound is not taken from the immutable artifact but re-derived from EXP-037's
    independently verified witness, and the two must agree or the run fails.

A dominating CSS code must be certified *exact* by the same standard in either
source: a verified minimum-weight witness at ``U`` plus UNSAT at ``U - 1`` on
both symplectic sides.  EXP-036 supplies these as tie-branch DOMINATION parents;
EXP-037 supplies them as ``css_*.json`` pool records, whose stored witness is
re-verified through two independent GF(2) paths here before use.  Domination is
witnessed by ``k_css >= k_pbb`` and ``d_css_exact >= U_b(reversal)``.

Fail-closed: if any reversal lacks an explicit dominator, the artifact routes to
quarantine and the claim MUST NOT be made.  This keeps the "CSS envelope stands"
sentence machine-checked instead of prose.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_SPEC = importlib.util.spec_from_file_location(
    "exp036_delta_closure", ROOT / "experiments" / "exp036_delta_closure.py"
)
assert _SPEC and _SPEC.loader
E36 = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(E36)

_SPEC37 = importlib.util.spec_from_file_location(
    "exp037_envelope_classification",
    ROOT / "experiments" / "exp037_envelope_classification.py",
)
assert _SPEC37 and _SPEC37.loader
E37 = importlib.util.module_from_spec(_SPEC37)
_SPEC37.loader.exec_module(E37)

OUTPUT = ROOT / "results" / "processed" / "exp036_envelope_check.json"
QUARANTINE = ROOT / "results" / "quarantine" / "exp036_envelope_check.json"


def exp037_exact_css() -> list[dict]:
    """Certified-exact CSS BB codes from EXP-037's pool.

    Same exactness standard as EXP-036's tie branch (verified witness at ``U``
    plus UNSAT at ``U - 1`` on both sides).  The stored witness is re-verified
    through two independent GF(2) paths here; a record that fails is dropped
    rather than trusted, and one that claims exactness without a distance is a
    hard error.
    """

    import numpy as np

    rows = E36.E27.load_catalogue()
    out: list[dict] = []
    for record in E37.load_css_pool():
        if not record.get("exact"):
            continue
        if record.get("distance") is None:
            raise RuntimeError(
                f"css record {record.get('css_fingerprint')} claims exact "
                "without a distance"
            )
        index = int(record["source_catalogue_index"])
        _, instances = E37.css_instances(rows[index])
        instance = instances[record["witness_side"]]
        vector = np.asarray(record["witness_vector"], dtype=np.uint8)
        checked = E37.verify_witness_two_paths(instance, vector)
        if not checked["valid"] or int(vector.sum()) != int(
            record["certified_upper_bound"]
        ):
            continue
        out.append(
            {
                "catalogue_index": index,
                "row": f"row_{index:04d}",
                "n": int(record["n"]),
                "k_css": int(record["k"]),
                "d_css_exact": int(record["distance"]),
                "parent_fingerprint": record["css_fingerprint"],
                "source": "exp037_pool",
            }
        )
    return out


def exp027_reversals(artifact: dict, rows: list[dict]) -> list[dict]:
    """EXP-027's own certified reversals, outside EXP-036's undecided set.

    The immutable artifact's ``pbb_d_upper_bound`` is a cross-check only: the
    bound actually used is EXP-037's independently verified witness weight, and
    the two must agree.
    """

    out: list[dict] = []
    for item in artifact["per_item_results"]:
        if item.get("verdict") != "CERTIFIED_REVERSAL":
            continue
        index = int(item["catalogue_index"])
        path = E37.row_record_path(index)
        if not path.exists():
            raise RuntimeError(
                f"row {index}: EXP-027 reversal has no EXP-037 witness record; "
                "run exp037 for this row before claiming domination"
            )
        record = json.loads(path.read_text(encoding="utf-8"))
        fresh = int(record["upper_bound"])
        if fresh != int(item["pbb_d_upper_bound"]):
            raise RuntimeError(
                f"row {index}: EXP-037 witness weight {fresh} disagrees with "
                f"EXP-027 upper bound {item['pbb_d_upper_bound']}"
            )
        out.append(
            {
                "catalogue_index": index,
                "row": f"row_{index:04d}",
                "label": item["label"],
                "n": int(item["n"]),
                "k_pbb": int(rows[index]["k"]),
                "d_pbb_upper": fresh,
                "source": "exp027_artifact",
            }
        )
    return out


def build_payload() -> tuple[dict, bool]:
    rows = E36.E27.load_catalogue()
    artifact = json.loads(E36.EXP027_ARTIFACT.read_text(encoding="utf-8"))
    by_index = {
        item["catalogue_index"]: item for item in artifact["per_item_results"]
    }
    exact_css: list[dict] = exp037_exact_css()
    reversals: list[dict] = exp027_reversals(artifact, rows)
    for index, row in E36.undecided_rows():
        state = E36.verified_row_state(index, row)
        if state is None:
            continue
        if state.get("verdict") not in {
            "DOMINATION_PROVED",
            "CERTIFIED_REVERSAL",
        }:
            continue
        # Nothing stored is trusted: re-derive bounds and verdict from
        # re-verified witnesses and call statuses, and require replay stamps
        # hash-bound to freshly rebuilt CNFs.
        instances = E36.build_instances(row)
        bounds = E36.validated_bounds(instances, state)
        verdict, _ = E36.derive_verdict(bounds)
        if verdict != state["verdict"]:
            raise RuntimeError(
                f"row {index}: stored verdict {state['verdict']!r} does not "
                f"re-derive from validated evidence (got {verdict!r})"
            )
        if not E36.replay_stamps_valid(instances, verdict, bounds, state):
            # Unstamped rows contribute nothing; a reversal without stamps
            # will fail the certified-set equality check below.
            continue
        if verdict == "DOMINATION_PROVED":
            U_p = bounds["U_p"]
            parent_lower_gt = bounds["parent_lower_gt"]
            # Exactness from validated evidence only: witness at U_p plus
            # both-sides UNSAT at >= U_p - 1 pins d_css = U_p.
            if U_p is None or parent_lower_gt is None or parent_lower_gt < U_p - 1:
                continue
            exact_css.append(
                {
                    "catalogue_index": index,
                    "row": f"row_{index:04d}",
                    "n": state["n"],
                    "k_css": by_index[index]["parent_k"],
                    "d_css_exact": U_p,
                    "parent_fingerprint": state["parent_fingerprint"],
                    "source": "exp036_tie_branch",
                }
            )
        else:
            reversals.append(
                {
                    "catalogue_index": index,
                    "row": f"row_{index:04d}",
                    "label": state["label"],
                    "n": state["n"],
                    "k_pbb": int(rows[index]["k"]),
                    "d_pbb_upper": bounds["U_b"],
                    "source": "exp036_tie_branch",
                }
            )
    # Non-vacuity: the reversal set must equal the union of every certified
    # reversal set on record.  Subtracting a source would let the claim shrink
    # to whatever happens to be dominated.
    aggregate_path = (
        E36.PROCESSED if E36.PROCESSED.exists() else E36.PARTIAL
    )
    aggregate = json.loads(aggregate_path.read_text(encoding="utf-8"))
    certified_set = {
        item["catalogue_index"]
        for item in aggregate["items"]
        if item["verdict"] == "CERTIFIED_REVERSAL"
    } | {
        int(item["catalogue_index"])
        for item in artifact["per_item_results"]
        if item.get("verdict") == "CERTIFIED_REVERSAL"
    }
    found_set = {reversal["catalogue_index"] for reversal in reversals}
    if not certified_set or found_set != certified_set:
        raise RuntimeError(
            "envelope check is vacuous or inconsistent: certified reversal "
            f"set {sorted(certified_set)} vs re-derived stamped set "
            f"{sorted(found_set)}"
        )
    checks = []
    all_dominated = True
    for reversal in reversals:
        dominators = [
            candidate
            for candidate in exact_css
            if candidate["n"] == reversal["n"]
            and candidate["k_css"] >= reversal["k_pbb"]
            and candidate["d_css_exact"] >= reversal["d_pbb_upper"]
        ]
        if not dominators:
            all_dominated = False
        # Canonical choice: the strongest dominator by k*d^2, so the reported
        # comparison is well defined rather than dependent on iteration order.
        best = max(
            dominators,
            key=lambda c: (c["k_css"] * c["d_css_exact"] ** 2, c["k_css"], c["d_css_exact"]),
            default=None,
        )
        checks.append(
            {
                **reversal,
                "dominated": bool(dominators),
                "dominator": best,
                "num_dominators": len(dominators),
                "kd2_over_n_pbb": reversal["k_pbb"] * reversal["d_pbb_upper"] ** 2 / reversal["n"],
                "kd2_over_n_css": (
                    best["k_css"] * best["d_css_exact"] ** 2 / reversal["n"] if best else None
                ),
            }
        )
    payload = {
        "experiment": "EXP-036-envelope",
        "schema": "exp036-envelope-check-v2",
        "protocol": __doc__,
        "generated_utc": E36.utc_now(),
        "num_reversals_checked": len(reversals),
        "num_exact_css_candidates": len(exact_css),
        "all_reversals_css_dominated_at_equal_n": all_dominated,
        "checks": sorted(checks, key=lambda c: c["catalogue_index"]),
        "exact_css_candidates": sorted(
            exact_css, key=lambda c: c["catalogue_index"]
        ),
    }
    return payload, all_dominated


def main() -> int:
    payload, all_dominated = build_payload()
    destination = OUTPUT if all_dominated else QUARANTINE
    E36.atomic_write_json(destination, payload)
    print(
        json.dumps(
            {
                "all_reversals_css_dominated_at_equal_n": all_dominated,
                "reversals": payload["num_reversals_checked"],
                "candidates": payload["num_exact_css_candidates"],
                "wrote": str(destination.relative_to(ROOT)),
            },
            indent=2,
        )
    )
    return 0 if all_dominated else 2


if __name__ == "__main__":
    raise SystemExit(main())
