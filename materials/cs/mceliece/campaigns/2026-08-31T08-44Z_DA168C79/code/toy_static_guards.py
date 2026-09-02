"""Static guard runner — campaign 2026-08-31T08-44Z_DA168C79.

Runs ONLY trivial structural checks (label-set cardinality, multiindex
counts, one toy public record build, index bijectivity).  It constructs
NO census matrices for the 24-cell population and runs NO registered
control population.  CPU slot remains Delcap-owned.
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import toy_census as tc


def main():
    assert os.environ.get(tc.RELEASE_KEY, "0") != "1", \
        "release flag must NOT be set for the static guard run"
    guards = tc.static_guards()
    out = {
        "module": "toy_census.py",
        "static_only": True,
        "release_flag_check": "absent (required)",
        "guards": guards,
        "p2_p3_scan_implementation": "skeleton present; phase-2 system "
                                     "builder gated behind released runner",
    }
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
