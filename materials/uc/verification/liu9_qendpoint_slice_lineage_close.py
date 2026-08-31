#!/usr/bin/env python3
"""Lineage-resolved closure of one qendpoint-lift slice residual (generic K).

Thin parameterization of `liu9_qendpoint_slice1_lineage_close.py`: identical
verifier — same `resolve_lineage`, same certified bound stack, same Arb
precision 480, same generation cap 120, same tube rho=1/1728, same digest
check of the frozen input — with `--source`/`--output` free so the queue can
close slices 02..16.  Slice 1's committed verifier stays untouched and
remains the reproduction entry for its committed report.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import liu9_qendpoint_slice1_lineage_close as base


RESULTS = HERE / "results"
INPUT_DEFAULT = RESULTS / "liu9-qendpoint-slice-01-of-16.json"
OUTPUT_DEFAULT = RESULTS / "liu9-qendpoint-slice01-lineage-close.json"


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=INPUT_DEFAULT)
    parser.add_argument("--output", type=Path, default=OUTPUT_DEFAULT)
    args = parser.parse_args(argv)
    report = base.run(args.source, args.output)
    closing = report["closing_run"]
    print("LIU H2 SLICE RESIDUAL LINEAGE CLOSURE (generic)")
    print("verdicts %s" % closing["verdicts"])
    print("max generations %d, elapsed %.0fs"
          % (closing["max_generations_needed"], closing["elapsed_seconds"]))
    print("claim_status %s" % report["claim_status"])
    print("report_sha256 %s" % report["report_sha256"])
    return 0 if closing["unresolved_lineages"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
