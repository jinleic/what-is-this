#!/usr/bin/env python3
"""Fail-closed tests for the exhaustive entropy-bridge control.

A passing control is only evidence if a broken construction would fail it.
These tests run the real check on every family with at most three coordinates
and then break the two load-bearing ingredients — the ``s*`` clipping law and
the prefix conditional probability — and require the control to reject.
"""

from __future__ import annotations

import importlib.util
import unittest
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE_PATH = HERE / "entropy_bridge_exhaustive.py"


def load():
    spec = importlib.util.spec_from_file_location("uc_entropy_bridge_test", MODULE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load %s" % MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    from flint import ctx

    ctx.prec = module.PRECISION_BITS
    return module


def first_rejection(module, max_coordinates):
    """Return the number of families checked before the first rejection."""
    checked = 0
    for coordinate_count in range(1, max_coordinates + 1):
        for family in module.enumerate_families(coordinate_count):
            checked += 1
            try:
                module.check_family(family, coordinate_count)
            except module.CheckError:
                return checked
    return None


class EntropyBridgeExhaustiveTests(unittest.TestCase):
    def test_every_small_family_passes(self):
        module = load()
        families = 0
        union_closed = 0
        for coordinate_count in range(1, 4):
            for family in module.enumerate_families(coordinate_count):
                summary = module.check_family(family, coordinate_count)
                families += 1
                union_closed += int(summary["union_closed"])
                self.assertIn(
                    summary["iid_class"], ("positive", "nonnegative", "equality")
                )
                self.assertIn(
                    summary["coupled_class"], ("positive", "nonnegative", "equality")
                )
        self.assertEqual(families, 273)
        self.assertEqual(union_closed, 137)

    def test_dropping_the_max_branch_of_sstar_is_rejected(self):
        module = load()
        module.sstar = lambda p, r: min(p + r, module.HALF)
        self.assertIsNotNone(first_rejection(module, 3))

    def test_constant_prefix_probability_is_rejected(self):
        module = load()
        module.prefix_parameter = lambda family, prefix, coordinate: Fraction(1, 2)
        self.assertIsNotNone(first_rejection(module, 3))

    def test_double_counting_the_union_bit_is_rejected(self):
        module = load()
        original = module.sstar
        module.sstar = lambda p, r: min(original(p, r) + Fraction(1, 100), Fraction(1))
        self.assertIsNotNone(first_rejection(module, 3))


if __name__ == "__main__":
    unittest.main()
