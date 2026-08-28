#!/usr/bin/env python3
"""Mutation tests for the exhaustive-cover certificate.

A control that cannot fail proves nothing.  Each test below rejects one
specific way the cover argument could be wrong, and the four the task names --
an omitted endpoint, an unenclosed origin segment, an inflated radius, and a
kappa above the certified ceiling -- are covered by
`test_dropping_an_endpoint_cell_is_caught`,
`test_a_cover_missing_the_centre_is_rejected`,
`test_an_inflated_radius_is_rejected`, and
`test_kappa_above_the_ceiling_is_rejected` respectively.

Run directly; this module is not importable as `uc.verification.*`:

    ./.venv/bin/python -I -B uc/verification/test_liu9_chart_cover.py
"""

from __future__ import annotations

import os
import sys
import unittest
from fractions import Fraction as F

for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

HERE = os.path.dirname(os.path.abspath(__file__))
UC = os.path.dirname(HERE)
for _path in (UC, HERE):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from flint import arb  # noqa: E402

from liu9_binding import (  # noqa: E402
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_boundary_layer import _arbf  # noqa: E402
from liu9_chart_cover import (  # noqa: E402
    CERTIFIED_KAPPA,
    Axis,
    Cover,
    Partition,
    build_axis,
    build_cover,
    build_q_partition,
    cell_hessians,
    partition_of_axis,
    sweep_cover,
    verify_axis_is_exhaustive,
    verify_cover_is_exhaustive,
    verify_partition_is_exhaustive,
    certify_centre_identity,
    collapsed_centre_gap,
    collapsed_centre_jet,
    certify_gd_antisymmetry,
)
from liu9_chart_centered import (  # noqa: E402
    pencil_is_psd,
    smooth_chart_ceiling,
)

MP = solve_equation_parameters(90)
ARB = certify_equation_parameters(MP)


class CoverGeometryTest(unittest.TestCase):
    """The exact-rational tiling checks."""

    def test_the_honest_cover_is_exhaustive(self):
        cover = build_cover(F(1, 32), F(1, 512), F(1, 4))
        report = verify_cover_is_exhaustive(cover)
        self.assertEqual(report["s"]["count"], 16)
        self.assertEqual(report["d"]["count"], 16)
        self.assertEqual(report["q"]["count"], 256)
        self.assertEqual(cover.cell_count, 16 * 16 * 256)
        self.assertTrue(report["contains_chart_centre"])

    def test_dropping_an_endpoint_cell_is_caught(self):
        """The outermost cell of the cover cannot be silently omitted.

        This is the 'omitted radial endpoint' mutation.  Shortening the count
        by one leaves the union ending at 15/512 instead of 1/32, so the box is
        not covered and the sweep would be certifying a strictly smaller
        region than the one it reports.
        """
        honest = build_axis(-F(1, 32), F(1, 32), F(1, 512))
        verify_axis_is_exhaustive(honest)
        truncated = Axis(honest.lower, honest.upper, honest.count - 1,
                         honest.half)
        with self.assertRaises(AssertionError) as caught:
            verify_axis_is_exhaustive(truncated)
        self.assertIn("endpoint cell is missing", str(caught.exception))

    def test_dropping_an_interior_cell_breaks_abutment(self):
        """A gap strictly inside the cover must not pass.

        Only one interior cell is displaced, so both endpoints still land
        exactly on the box boundary and the endpoint check stays silent.  The
        abutment identity is then the only thing standing between a gapped
        cover and a false certificate, which is what this asserts.
        """
        class Gapped(Axis):
            def centre(self, index: int) -> F:
                shift = F(1, 4096) if index == 8 else F(0)
                return super().centre(index) + shift

        gapped = Gapped(-F(1, 32), F(1, 32), 16, F(1, 512))
        first_lower, _ = gapped.bounds(0)
        _, last_upper = gapped.bounds(gapped.count - 1)
        self.assertEqual(first_lower, -F(1, 32))
        self.assertEqual(last_upper, F(1, 32))
        with self.assertRaises(AssertionError) as caught:
            verify_axis_is_exhaustive(gapped)
        self.assertIn("do not abut", str(caught.exception))

    def test_an_inflated_radius_is_rejected(self):
        """Claiming a larger box than the cells tile must fail.

        The cells are sized for 1/32; asserting they cover 1/16 leaves the
        union short at both ends.
        """
        cells = build_axis(-F(1, 32), F(1, 32), F(1, 512))
        inflated = Axis(-F(1, 16), F(1, 16), cells.count, cells.half)
        with self.assertRaises(AssertionError) as caught:
            verify_axis_is_exhaustive(inflated)
        self.assertIn("endpoint cell is missing", str(caught.exception))

    def test_a_non_dividing_half_width_is_rejected(self):
        """A cell size that does not tile the span exactly cannot be used."""
        with self.assertRaises(ValueError):
            build_axis(-F(1, 32), F(1, 32), F(1, 300))

    def test_a_cover_missing_the_centre_is_rejected(self):
        """The 'unenclosed origin segment' mutation.

        An annulus-shaped region would pass every pointwise PSD test and prove
        nothing, because the mean-value step needs the segment from the chart
        centre.  Moving the covered band off the centre must be refused.
        """
        offset = build_cover(F(1, 32), F(1, 512), F(1, 4))
        shifted = Cover(
            offset.radius, offset.delta, offset.q_floor,
            build_axis(F(1, 64), F(3, 64), F(1, 512)),
            offset.d_axis,
            offset.q_axis,
        )
        with self.assertRaises(AssertionError) as caught:
            verify_cover_is_exhaustive(shifted)
        self.assertIn("excludes the chart centre", str(caught.exception))


class QPartitionTest(unittest.TestCase):
    """The octave tail on the q axis, and why a uniform grid is not enough."""

    def test_the_octave_tail_is_exhaustive_and_ascending(self):
        partition = build_q_partition(F(1, 4096), F(1, 4), 64, F(1, 512))
        report = verify_partition_is_exhaustive(partition)
        self.assertEqual(report["count"], 1536)
        self.assertEqual(partition.lower, F(1, 4096))
        self.assertEqual(partition.upper, F(4095, 4096))
        self.assertEqual(partition.widest_relative(), F(1, 64))
        lowers = [cell[0] for cell in partition.cells]
        self.assertEqual(lowers, sorted(lowers),
                         "the tail must ascend; a flat reversal also "
                         "reverses the pieces inside each octave")

    def test_a_uniform_q_grid_is_detectably_insufficient(self):
        """The finding that a uniform q axis hides.

        Both `gap.hdd` and `dist.hdd` carry an exact factor q(1-q), so the
        pencil's dd entry vanishes linearly at the endpoints while a
        fixed-width cell's enclosure error does not.  A uniform grid therefore
        stops certifying well before q reaches the 1/4096 that ingredient
        (iii) needs, and this test pins that failure rather than letting a
        future edit quietly reintroduce it.
        """
        kappa = _arbf(CERTIFIED_KAPPA)
        uniform = cell_hessians(
            ARB, F(0), F(0), F(1, 512), F(1, 128), F(1, 128) + F(1, 512),
        )
        self.assertFalse(pencil_is_psd(*uniform, kappa),
                         "a uniform-width cell at q=1/128 must fail")
        geometric = cell_hessians(
            ARB, F(0), F(0), F(1, 512),
            F(1, 128), F(1, 128) + F(1, 128) / 64,
        )
        self.assertTrue(pencil_is_psd(*geometric, kappa),
                        "the octave cell at the same q must succeed")

    def test_the_dd_entry_carries_an_exact_q_factor(self):
        """m22/(q(1-q)) must be the same constant across four decades.

        If this ever stops holding, the octave grid is no longer the right
        device and the whole small-q argument has to be revisited.
        """
        from liu9_chart_centered import (
            Jet3, distance_squared_jet, gap_jet,
        )
        kappa = _arbf(CERTIFIED_KAPPA)
        ratios = []
        for q_value in (F(1, 4), F(1, 64), F(1, 1024), F(1, 4096)):
            q = _arbf(q_value)
            centre = Jet3.variable(ARB.x, "s")
            zero = Jet3.variable(arb(0), "d")
            gap = gap_jet(ARB, q, centre, zero, arb(0))
            distance = distance_squared_jet(ARB, q, centre, zero, arb(0))
            m22 = gap.hdd - kappa * distance.hdd
            ratios.append(m22 / (_arbf(q_value) * (1 - _arbf(q_value))))
        for ratio in ratios[1:]:
            self.assertTrue(
                ratio.overlaps(ratios[0]),
                "the q(1-q) factor is not exact; the octave grid rests on it",
            )
        self.assertGreater(float(ratios[0].lower()), 0.6)

    def test_a_non_power_of_two_tail_is_rejected(self):
        """The octave tail must land exactly on q_min."""
        with self.assertRaises(ValueError):
            build_q_partition(F(1, 3000), F(1, 4), 64, F(1, 512))

    def test_the_tail_must_stay_inside_the_open_interval(self):
        with self.assertRaises(ValueError):
            build_q_partition(F(1, 2), F(1, 4), 64, F(1, 512))


class CoverPencilTest(unittest.TestCase):
    """The interval arithmetic behind one cell."""

    def test_kappa_above_the_ceiling_is_rejected(self):
        """The sweep must refuse a kappa the smooth chart cannot support."""
        ceiling = smooth_chart_ceiling(ARB)
        too_large = F(1, 2)
        self.assertGreater(float(too_large), float(ceiling.upper()))
        cover = build_cover(F(1, 256), F(1, 512), F(1, 4))
        with self.assertRaises(AssertionError) as caught:
            sweep_cover(ARB, cover, too_large)
        self.assertIn("exceeds the smooth-chart ceiling",
                      str(caught.exception))

    def test_the_certified_kappa_holds_on_a_sampled_cell(self):
        """A cell at the far corner must actually pass at the reported kappa."""
        gap, distance = cell_hessians(
            ARB, F(1, 32) - F(1, 512), F(1, 32) - F(1, 512), F(1, 512),
            F(1, 4), F(1, 4) + F(1, 512),
        )
        self.assertTrue(pencil_is_psd(gap, distance, _arbf(CERTIFIED_KAPPA)))

    def test_a_coarser_cell_is_detectably_weaker(self):
        """delta=1/256 must fail where delta=1/512 succeeds.

        This is the measurement that fixed the cell size.  If a later edit
        widened the cells, this test says the enclosure genuinely stops
        certifying rather than merely getting looser.
        """
        kappa = _arbf(CERTIFIED_KAPPA)
        fine = cell_hessians(
            ARB, F(1, 32) - F(1, 512), F(0), F(1, 512),
            F(1, 2), F(1, 2) + F(1, 512),
        )
        coarse = cell_hessians(
            ARB, F(1, 32) - F(1, 256), F(0), F(1, 256),
            F(1, 2), F(1, 2) + F(1, 256),
        )
        self.assertTrue(pencil_is_psd(*fine, kappa))
        self.assertFalse(pencil_is_psd(*coarse, kappa))

    def test_the_expansion_point_is_the_cell_centre_not_the_chart_centre(self):
        """The whole gain comes from re-centering; verify it actually happens.

        If `cell_hessians` still expanded at the chart centre, a far cell and a
        near cell with the same q would return identical Hessians.  They must
        not.
        """
        near = cell_hessians(
            ARB, F(0), F(0), F(1, 512), F(1, 2), F(1, 2) + F(1, 512),
        )
        far = cell_hessians(
            ARB, F(1, 32) - F(1, 512), F(0), F(1, 512),
            F(1, 2), F(1, 2) + F(1, 512),
        )
        self.assertFalse(
            near[0].ss.overlaps(far[0].ss)
            and near[0].dd.overlaps(far[0].dd)
            and float(near[0].ss.mid()) == float(far[0].ss.mid()),
            "the two cells returned the same enclosure; the expansion point "
            "is not tracking the cell",
        )

    def test_a_widened_q_cell_loses_the_certificate(self):
        """q resolution, not geometry, is what limits the cell size."""
        kappa = _arbf(CERTIFIED_KAPPA)
        tight = cell_hessians(
            ARB, F(0), F(0), F(1, 512), F(1, 2), F(1, 2) + F(1, 512),
        )
        wide = cell_hessians(
            ARB, F(0), F(0), F(1, 512), F(1, 2), F(1, 2) + F(1, 32),
        )
        self.assertTrue(pencil_is_psd(*tight, kappa))
        self.assertFalse(pencil_is_psd(*wide, kappa))
        tight_width = float(tight[0].ss.upper()) - float(tight[0].ss.lower())
        wide_width = float(wide[0].ss.upper()) - float(wide[0].ss.lower())
        self.assertGreater(wide_width, 4 * tight_width)

    def test_a_small_cover_certifies_end_to_end(self):
        """The inherited radius must reproduce through the new code path."""
        cover = build_cover(F(1, 256), F(1, 512), F(1, 4))
        verify_cover_is_exhaustive(cover)
        state = sweep_cover(ARB, cover, CERTIFIED_KAPPA)
        self.assertEqual(state.processed, cover.cell_count)
        self.assertEqual(state.failures, [])
        self.assertIsNotNone(state.worst_margin)
        self.assertGreater(state.worst_margin, 0.0)


class CentreIdentityTest(unittest.TestCase):
    """Step B's hypothesis, and the honesty of the deficit it reports."""

    def test_the_collapsed_form_matches_gap_jet_at_point_q(self):
        """The q-free transcription must be the same function.

        `collapsed_centre_gap` hand-collapses the chart at its centre.  If that
        derivation were wrong the module would report a confident, q-free, and
        false constant, so it is checked against the general `gap_jet` at q
        spanning four decades.
        """
        collapsed = collapsed_centre_gap(ARB)
        from liu9_chart_centered import Jet3, gap_jet
        for q_value in (F(1, 4096), F(1, 64), F(1, 2), F(4095, 4096)):
            gap = gap_jet(
                ARB, _arbf(q_value), Jet3.variable(ARB.x, "s"),
                Jet3.variable(arb(0), "d"), arb(0),
            )
            self.assertTrue(
                gap.v.overlaps(collapsed),
                f"the collapsed centre form disagrees with gap_jet at q="
                f"{q_value}",
            )

    def test_the_collapsed_form_carries_no_q(self):
        """Two calls cannot differ, because no q is passed in at all."""
        first = collapsed_centre_gap(ARB)
        second = collapsed_centre_gap(ARB)
        self.assertEqual(first.str(30), second.str(30))
        self.assertLess(float(abs(first).upper()), 1e-60)

    def test_the_reported_deficit_is_not_silently_zero(self):
        """The module must admit the residual rather than round it away.

        Claiming an exact zero here would be the single easiest way to
        overstate the result, so the deficit is required to be positive and
        of the size the binding solve actually leaves.
        """
        report = certify_centre_identity(ARB, F(1, 32))
        self.assertGreater(report["step_b_deficit"], 0.0)
        self.assertLess(report["step_b_deficit"], 1e-60)
        self.assertGreater(report["worst_centre_gs"], 0.0)

    def test_an_interval_q_centre_evaluation_is_useless(self):
        """Why the collapsed form exists at all.

        If a later edit replaced `collapsed_centre_gap` with a plain
        interval-q call, this is the width it would silently accept.
        """
        from liu9_chart_centered import Jet3, gap_jet
        q = _arbf(F(1, 4)).union(_arbf(F(1, 4) + F(1, 512)))
        gap = gap_jet(
            ARB, q, Jet3.variable(ARB.x, "s"), Jet3.variable(arb(0), "d"),
            arb(0),
        )
        width = float(gap.v.upper()) - float(gap.v.lower())
        self.assertGreater(
            width, 1e-4,
            "interval-q evaluation at the centre is expected to be wide; if "
            "it became tight, the collapsed form is no longer needed",
        )
        self.assertLess(float(abs(collapsed_centre_gap(ARB)).upper()), 1e-60)


class AllQCentreBoundTest(unittest.TestCase):
    """The Step B deficit must hold for every q, not just the probe set."""

    def test_the_collapsed_jet_reproduces_gs_at_point_q(self):
        """The q-free s-gradient must be the same number gap_jet gives.

        `gs` dominates the deficit, so if this transcription were wrong the
        reported all-q bound would be confidently wrong rather than merely
        loose.
        """
        from liu9_chart_centered import Jet3, gap_jet
        jet = collapsed_centre_jet(ARB)
        for q_value in (F(1, 4096), F(1, 64), F(1, 2), F(4095, 4096)):
            gap = gap_jet(
                ARB, _arbf(q_value), Jet3.variable(ARB.x, "s"),
                Jet3.variable(arb(0), "d"), arb(0),
            )
            self.assertTrue(
                gap.gs.overlaps(jet.gs),
                f"collapsed gs disagrees with gap_jet at q={q_value}",
            )
            self.assertTrue(gap.v.overlaps(jet.v))

    def test_the_gd_bracket_is_the_zero_polynomial(self):
        report = certify_gd_antisymmetry()
        self.assertTrue(report["identically_zero"])
        for power in ("q^0", "q^1", "q^2"):
            self.assertEqual(report["coefficients"][power], "0")

    def test_a_mis_signed_antisymmetry_bracket_is_detected(self):
        """The cancellation must come from the weights, not be assumed.

        If the component weights were paired the wrong way round the bracket
        would be -2q(1-q), not zero, and the d-gradient would NOT vanish.  The
        zero-polynomial test must be able to see that difference.
        """
        def wrong(q):
            return (1 - q) * (-q) - q * (1 - q)

        at0, at1, at2 = (wrong(F(k)) for k in (0, 1, 2))
        quadratic = (at2 - 2 * at1 + at0) / 2
        linear = at1 - at0 - quadratic
        self.assertFalse(at0 == 0 and linear == 0 and quadratic == 0,
                         "the mis-signed bracket must not look like zero")
        self.assertEqual(quadratic, F(2))
        self.assertEqual(linear, F(-2))

    def test_the_deficit_no_longer_depends_on_the_probe_set(self):
        """Shrinking the probe set must not change the reported bound.

        The probes are a transcription guard.  If the deficit still moved with
        them, it would still be a sampled quantity wearing an all-q label.
        """
        wide = certify_centre_identity(ARB, F(1, 32))
        narrow = certify_centre_identity(ARB, F(1, 32), probe_q=(F(1, 2),))
        self.assertEqual(wide["step_b_deficit"], narrow["step_b_deficit"])
        self.assertTrue(wide["deficit_is_all_q"])
        self.assertEqual(wide["q_free_gs_bound"], narrow["q_free_gs_bound"])

    def test_the_deficit_scales_with_the_reach(self):
        """It is |value| + |gs|*R, so a larger box must report a larger bound."""
        small = certify_centre_identity(ARB, F(1, 256))["step_b_deficit"]
        large = certify_centre_identity(ARB, F(1, 32))["step_b_deficit"]
        self.assertGreater(large, small)
        self.assertLess(large, 1e-60)


if __name__ == "__main__":
    unittest.main(verbosity=2)
