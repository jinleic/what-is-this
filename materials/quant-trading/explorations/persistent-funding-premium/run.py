#!/usr/bin/env python3
"""Directory entry point for the preregistered persistent-funding-premium test.

The canonical implementation is the campaign-level runner:
    quant-trading/explorations/run_persistent_funding_premium.py
It verifies the frozen contract hash, verifies every upstream archive
checksum, then resolves the single holdout verdict (and the preregistered
replication only on a pass). This shim exists so the direction directory
carries a run.py per the campaign output contract without duplicating the
frozen logic.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

CANONICAL = (
    Path(__file__).resolve().parents[1] / "run_persistent_funding_premium.py"
)


def load_canonical():
    spec = importlib.util.spec_from_file_location(
        "run_persistent_funding_premium", CANONICAL
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load canonical runner: {CANONICAL}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


if __name__ == "__main__":
    raise SystemExit(load_canonical().main(sys.argv[1:]))
