"""EXP-049: dedicated certification of the three X-U demotion-immune parents.

The EXP-039 fleet sweeps catalogue parents serially at nice 15; the three
X-U parents (``phase2_109``, ``15_6_0256``, ``30_6_0289``) are the last open
X-side cells of the catalogue partition.  This driver runs each one directly,
reusing EXP-039's ``certify_parent`` so results, fingerprints, and the atomic
on-disk certificates are byte-identical to what the fleet would produce
(idempotent: if the fleet lands first, the stored record is reused).

Usage (internal): one parent label per process.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

_SPEC = importlib.util.spec_from_file_location(
    "exp039_nogo_module", ROOT / "experiments" / "exp039_nogo_module.py"
)
assert _SPEC and _SPEC.loader
E39 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E39
_SPEC.loader.exec_module(E39)

TRIO_LABELS = {"15_6_0256", "phase2_109", "30_6_0289"}


def run(args: argparse.Namespace) -> int:
    rows = E39.E27.load_catalogue()
    parents = E39.distinct_parents(rows)
    bounds = E39.load_pool_lower_bounds()
    found = False
    for fp, entry in parents.items():
        labels = {m.get("label") for m in entry["members"]}
        if args.label not in labels:
            continue
        found = True
        payload = E39.certify_parent(entry, bounds.get(fp, 2))
        print(
            f"{args.label}: parent {fp[:16]} n={payload['n']} "
            f"k_P={payload['k_parent']} d_Z={payload['d_z_parent']} "
            f"T={payload['T']} exact={payload['T_is_exact']} "
            f"family_closed={payload['family_closed']}",
            flush=True,
        )
    if not found:
        raise SystemExit(f"label not found among parents: {args.label}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True, choices=sorted(TRIO_LABELS))
    return run(parser.parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
