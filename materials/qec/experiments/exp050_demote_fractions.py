"""EXP-050: exact demote-fraction census over quotient classes.

EXP-047 decided the per-parent demote/immune verdict with early-stop on the
first realized class.  The flagship's fraction (4095/4095 of its nonzero
quotient classes demote) suggested collapse directions are generic, but the
exact fraction was unknown for every other parent.  This experiment enumerates
ALL nonzero quotient classes for every distinct catalogue parent with
``2^k_P - 1 <= CLASS_CAP`` and counts demoting classes exactly, using the
identical demotion test as EXP-047 (``yvec not in rowspace(Lpre @ S_z @ Rk)``).

Certification built in:
* flagship ``12_6_0193`` must reproduce 4095/4095;
* EXP-047's stored witness class per demoting parent must lie in the
  demote-set (mismatch => hard failure);
* immune parents (J-E) must yield fraction exactly 0.

Sharding: ``run --stride S --offset O``; results concatenate to
``results/processed/exp050_demote_fractions.json`` via ``assemble``.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

_SPEC = importlib.util.spec_from_file_location(
    "exp047_exact_demotion_decision",
    ROOT / "experiments" / "exp047_exact_demotion_decision.py",
)
assert _SPEC and _SPEC.loader
E47 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E47
_SPEC.loader.exec_module(E47)
E27 = E47.E27  # noqa: N816
E39 = E47.E39  # noqa: N816
from qec_research.gf2.linalg import rank_np  # noqa: E402

SCHEMA = "exp050-demote-fractions-v1"
STATE_DIR = ROOT / "results" / "partial_runs" / "exp050"
OUT = ROOT / "results" / "processed" / "exp050_demote_fractions.json"
MAX_K = 12  # stated census scope: only parents with k_P <= 12 (<= 4095 classes)
# EXP-050-EXT widens to k_P <= 20 via --max-k (1,048,575 classes/parent at 20).
FLAGSHIP = "12_6_0193"
FLAGSHIP_FRACTION = (4095, 4095)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(f"{path.suffix}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


def fraction_parent(fp: str, entry: dict, row: dict, witness_class: int | None):
    """Enumerate every nonzero quotient class exactly for one parent."""

    t0 = time.time()
    ell, m = int(entry["ell"]), int(entry["m"])
    dim = ell * m
    _, HX, HZ = E27.parent_matrices(row)
    Kb = E47.nullspace_np(HX)
    k_P = int(entry["k_parent"])
    Q, R, kq = E47.quotient_and_projector(HZ, Kb)
    if kq != k_P:
        raise RuntimeError(f"{entry['members'][0]['label']}: quotient dim {kq} != k_P {k_P}")
    Rk = R[:, :k_P]
    Lpre = E47.nullspace_np(HX.T)

    demote = 0
    total = (1 << k_P) - 1
    witness_in_set = None
    for yint in range(1, total + 1):
        yvec = np.array([(yint >> i) & 1 for i in range(k_P)], np.uint8)
        z = (yvec @ Q) % 2
        Sz = E47.shift_matrix_of(z, ell, m)
        proj = (Lpre @ Sz % 2) @ Rk % 2
        if rank_np(np.vstack([proj, yvec[None, :]])) != rank_np(proj):
            demote += 1
            if witness_class is not None and yint == witness_class:
                witness_in_set = True
    if witness_class is not None:
        witness_in_set = bool(witness_in_set)
    return {
        "schema": SCHEMA,
        "fingerprint": fp,
        "label": entry["members"][0]["label"],
        "ell": ell,
        "m": m,
        "n": 2 * dim,
        "k_parent": k_P,
        "classes_total": total,
        "classes_demote": demote,
        "fraction": demote / total if total else None,
        "exp047_witness_in_demote_set": witness_in_set,
        "wall_s": round(time.time() - t0, 3),
        "utc": utc_now(),
    }


def run(args: argparse.Namespace) -> int:
    rows = E27.load_catalogue()
    parents = E39.distinct_parents(rows)

    items = [
        (fp, entry, rows[entry["members"][0]["catalogue_index"]])
        for fp, entry in sorted(
            parents.items(), key=lambda kv: (1 << kv[1]["k_parent"], kv[0])
        )
        if entry["k_parent"] <= args.max_k
    ]
    items = items[args.offset :: args.stride]
    print(f"parents selected: {len(items)}", flush=True)
    from collections import Counter

    per_k = Counter(entry["k_parent"] for _, entry, _ in items)
    print(f"class budget by k_P: {dict(sorted(per_k.items()))}", flush=True)
    print(
        f"total classes to test: {sum((1 << entry['k_parent']) - 1 for _, entry, _ in items)}",
        flush=True,
    )
    for fp, entry, row in items:
        path = STATE_DIR / f"fraction_{fp}.json"
        if path.exists() and not args.force:
            print(f"{entry['members'][0]['label']}: cached", flush=True)
            continue
        payload = fraction_parent(fp, entry, row, None)
        label = payload["label"]
        if label == FLAGSHIP:
            ok = (payload["classes_demote"], payload["classes_total"]) == FLAGSHIP_FRACTION
            if not ok:
                raise RuntimeError(
                    f"flagship fraction mismatch: {payload['classes_demote']}"
                    f"/{payload['classes_total']} != {FLAGSHIP_FRACTION}"
                )
        atomic_write_json(path, payload)
        print(
            f"{label}: {payload['classes_demote']}/{payload['classes_total']} "
            f"demote ({payload['fraction']:.6f}) k_P={payload['k_parent']} in {payload['wall_s']}s",
            flush=True,
        )
    return 0


def assemble(args: argparse.Namespace) -> int:
    records = [
        json.loads(p.read_text(encoding="utf-8"))
        for p in sorted(STATE_DIR.glob("fraction_*.json"))
        if json.loads(p.read_text(encoding="utf-8")).get("schema") == SCHEMA
    ]
    demoting = [r for r in records if r["classes_demote"] > 0]
    immune47 = {
        lbl
        for _, lbl in json.loads(
            (ROOT / "results/processed/exp047_exact_demotion_decision.json").read_text(
                encoding="utf-8"
            )
        )["immunity_list"]
    }
    immune_mismatch = [
        r["label"]
        for r in records
        if (r["label"] in immune47) != (r["classes_demote"] == 0)
    ]
    flagship = next((r for r in records if r["label"] == FLAGSHIP), None)
    payload = {
        "schema": SCHEMA,
        "experiment": "exp050_demote_fractions",
        "utc": utc_now(),
        "parents": len(records),
        "parents_demoting": len(demoting),
        "parents_immune": len(records) - len(demoting),
        "fraction_min_over_demoting": min(
            (r["fraction"] for r in demoting), default=None
        ),
        "fractions_distinct": sorted({round(r["fraction"], 6) for r in records}),
        "immune_mismatch_vs_exp047": immune_mismatch,
        "flagship_check": (
            (flagship["classes_demote"], flagship["classes_total"]) if flagship else None
        ),
        "records": records,
    }
    atomic_write_json(OUT, payload)
    print(json.dumps({k: v for k, v in payload.items() if k != "records"}, indent=2))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("run")
    p.add_argument("--stride", type=int, default=1)
    p.add_argument("--offset", type=int, default=0)
    p.add_argument("--max-k", type=int, default=MAX_K)
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=run)
    p = sub.add_parser("assemble")
    p.set_defaults(func=assemble)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
