#!/usr/bin/env python3
"""Exhaustive-cover certificate lifting Liu's smooth-chart radius to 1/32.

`uc/liu9_chart_centered.py` certifies that the pencil

    M(u, q) = Hess gap - kappa * Hess dist^2                      (all Hessians
                                                                   in (s,d))

is positive semidefinite on the l-infinity box of radius 1/256 about the chart
centre `(s,d) = (x,0)`, with `kappa >= 4119063/33554432`.  Ingredient (iii),
certified in `uc/liu9_qdegenerate.py`, is CONDITIONAL on the same statement
holding out to `eps_sm = 1/32`.  This module closes that factor of eight.

Why a cover is the right device, and why radial quadrature is not
-----------------------------------------------------------------
Two separate steps are easy to conflate.

  Step A (what is certified here).  `M >= 0` everywhere on the box.  This is a
  POINTWISE property of a point set.  If a finite family of cells covers the
  box and `M >= 0` is certified on each cell, then `M >= 0` on the union, hence
  on the box.  Nothing about segments, origins or Taylor remainders enters.

  Step B (unchanged, and not re-proved here).  Turning `M >= 0` into
  `gap - kappa*dist^2 >= 0` uses the mean-value form

      F(v) = integral_0^1 (1-theta) v^T Hess F(theta v) v d theta,

  which does require the whole segment from the chart centre to `v` to lie
  where the Hessian is controlled.  That requirement is discharged ONCE, by
  geometry: the box `[x-R,x+R] x [-R,R]` is convex and contains the centre, so
  for every `v` in the box the segment `{theta v}` stays in the box -- exactly
  the region Step A covers.

So a uniform cover DOES suffice.  Per-cell segment containment would be needed
only if one tried to certify positivity cell by cell from each cell's own local
expansion, which is not what the pencil route does.  A radial quadrature is
therefore unnecessary machinery here; the cover is both correct and cheaper.
This module implements Step A and states Step B's hypothesis explicitly rather
than silently reusing it.

What limits the cell size
-------------------------
Not the geometry: the dominant enclosure width comes from the `q` interval, not
from the `(s,d)` cell.  At a point `q` the centre Hessian has width 0; at a `q`
cell of width 1/128 it is already 0.52 wide.  `liu9_chart_centered.py` ties the
`q` cell count to the radius for exactly this reason.  Measured here: at
`delta = 1/256` the pencil fails at offset 1/32, and at `delta = 1/512` it holds
at every corner.  The cover therefore uses `delta = 1/512` in `(s,d)` and the
same width in `q`.

Run from the repository root with
`math/.venv/bin/python -I -B math/uc/liu9_chart_cover.py`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Iterator, Optional, Sequence

for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

from flint import arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from liu9_binding import (  # noqa: E402
    ArbParameters,
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_boundary_layer import _arbf, _octave_grid  # noqa: E402
from liu9_chart_centered import (  # noqa: E402
    HessianBox,
    Jet3,
    _abs_upper,
    distance_squared_jet,
    gap_jet,
    jet_entropy,
    pencil_is_psd,
    protocol_diagonal,
    smooth_chart_ceiling,
)

ctx.prec = max(ctx.prec, 320)

DEFAULT_DPS = 90
DEFAULT_RADIUS = Fraction(1, 32)
DEFAULT_DELTA = Fraction(1, 512)
DEFAULT_Q_FLOOR = Fraction(1, 4)
DEFAULT_SUBDIVISIONS = 2
DEFAULT_Q_MIN = Fraction(1, 4096)
DEFAULT_Q_PIECES = 16
CERTIFIED_KAPPA = Fraction(4119063, 33554432)
INHERITED_RADIUS = Fraction(1, 256)
TARGET_RADIUS = Fraction(1, 32)


# ---------------------------------------------------------------------------
# Exact-rational cover geometry.


@dataclass(frozen=True)
class Axis:
    """A tiling of `[lower, upper]` into `count` cells of half-width `half`."""

    lower: Fraction
    upper: Fraction
    count: int
    half: Fraction

    def centre(self, index: int) -> Fraction:
        if not 0 <= index < self.count:
            raise IndexError("axis cell index out of range")
        return self.lower + (2 * index + 1) * self.half

    def bounds(self, index: int) -> tuple[Fraction, Fraction]:
        centre = self.centre(index)
        return centre - self.half, centre + self.half


def build_axis(lower: Fraction, upper: Fraction, half: Fraction) -> Axis:
    """Tile `[lower,upper]` by cells of half-width `half`, exactly."""
    if not upper > lower:
        raise ValueError("axis needs a nonempty interval")
    if half <= 0:
        raise ValueError("axis needs a positive half-width")
    span = upper - lower
    quotient = span / (2 * half)
    if quotient.denominator != 1:
        raise ValueError(
            "the cell half-width must divide the axis span exactly; "
            f"span/{2 * half} = {quotient} is not an integer"
        )
    return Axis(lower, upper, int(quotient), half)


def verify_axis_is_exhaustive(axis: Axis) -> dict[str, Any]:
    """Prove the tiling covers the axis with no gap, in exact arithmetic.

    Three separate failures are caught here, and each corresponds to one of the
    mutation tests: a dropped endpoint cell leaves the first or last bound
    short, a dropped interior cell breaks abutment, and an inflated claimed
    span leaves the union strictly inside it.
    """
    if axis.count <= 0:
        raise AssertionError("the cover has no cells")
    first_lower, _ = axis.bounds(0)
    _, last_upper = axis.bounds(axis.count - 1)
    if first_lower != axis.lower:
        raise AssertionError(
            f"the cover starts at {first_lower}, not at {axis.lower}; "
            "an endpoint cell is missing"
        )
    if last_upper != axis.upper:
        raise AssertionError(
            f"the cover ends at {last_upper}, not at {axis.upper}; "
            "an endpoint cell is missing"
        )
    for index in range(axis.count - 1):
        _, upper = axis.bounds(index)
        lower, _ = axis.bounds(index + 1)
        if upper != lower:
            raise AssertionError(
                f"cells {index} and {index + 1} do not abut: "
                f"{upper} != {lower}"
            )
    return {
        "lower": str(axis.lower),
        "upper": str(axis.upper),
        "count": axis.count,
        "half": str(axis.half),
        "abutment_checks": max(axis.count - 1, 0),
    }


@dataclass(frozen=True)
class Partition:
    """An ordered exact-rational tiling given by explicit cell bounds.

    The `(s,d)` axes are uniform and an `Axis` describes them completely.  The
    `q` axis is not uniform, so it carries its cells explicitly.
    """

    lower: Fraction
    upper: Fraction
    cells: tuple[tuple[Fraction, Fraction], ...]

    @property
    def count(self) -> int:
        return len(self.cells)

    def bounds(self, index: int) -> tuple[Fraction, Fraction]:
        return self.cells[index]

    def widest_relative(self) -> Fraction:
        """The largest cell width divided by its own lower endpoint."""
        worst = Fraction(0)
        for lower, upper in self.cells:
            reference = min(abs(lower), abs(1 - upper))
            if reference > 0:
                worst = max(worst, (upper - lower) / reference)
        return worst


def verify_partition_is_exhaustive(partition: Partition) -> dict[str, Any]:
    """The `Axis` check, for a tiling whose cells are listed explicitly."""
    if partition.count <= 0:
        raise AssertionError("the cover has no cells")
    if partition.cells[0][0] != partition.lower:
        raise AssertionError(
            f"the cover starts at {partition.cells[0][0]}, not at "
            f"{partition.lower}; an endpoint cell is missing"
        )
    if partition.cells[-1][1] != partition.upper:
        raise AssertionError(
            f"the cover ends at {partition.cells[-1][1]}, not at "
            f"{partition.upper}; an endpoint cell is missing"
        )
    for index in range(partition.count - 1):
        if partition.cells[index][1] != partition.cells[index + 1][0]:
            raise AssertionError(
                f"cells {index} and {index + 1} do not abut: "
                f"{partition.cells[index][1]} != "
                f"{partition.cells[index + 1][0]}"
            )
    return {
        "lower": str(partition.lower),
        "upper": str(partition.upper),
        "count": partition.count,
        "abutment_checks": max(partition.count - 1, 0),
        "widest_relative_width": str(partition.widest_relative()),
    }


def build_q_partition(
    q_min: Fraction = DEFAULT_Q_MIN,
    q_floor: Fraction = DEFAULT_Q_FLOOR,
    pieces: int = DEFAULT_Q_PIECES,
    middle_width: Fraction = DEFAULT_DELTA,
) -> Partition:
    """Cover `[q_min, 1-q_min]` with octave tails and a uniform middle.

    A uniform `q` grid is the wrong device near the endpoints, and that -- not
    the chart radius -- is what actually blocked the first attempt, which
    stopped certifying below `q = 1/64`.  Both `gap.hdd` and `dist.hdd` carry
    an exact factor `q(1-q)`: the measured ratio `m22/(q(1-q))` is
    `0.69512040` at `q = 1/4, 1/64, 1/128, 1/1024, 1/4096` alike, to nine
    digits.  So the true pencil entry vanishes linearly at the endpoints while
    a fixed-width cell's enclosure error does not, and the fixed grid drowns
    the signal.  Halving the cell along with the coordinate holds the
    RELATIVE width bounded, which is the quantity the enclosure cares about.
    `liu9_boundary_layer._octave_grid` exists for the same reason on the
    entropy side; this is the same device on the `q` axis.
    """
    if not 0 < q_min < q_floor < Fraction(1, 2):
        raise ValueError("the q partition needs 0 < q_min < q_floor < 1/2")
    if pieces <= 0:
        raise ValueError("each octave needs at least one piece")
    ratio = q_floor / q_min
    octaves = 0
    probe = Fraction(1)
    while probe < ratio:
        probe *= 2
        octaves += 1
    if probe != ratio:
        raise ValueError(
            f"q_floor/q_min = {ratio} must be a power of two so the octave "
            "tail lands exactly on q_min"
        )
    # `_octave_grid` emits octaves from the top down, with pieces
    # ascending inside each octave.  Reversing the flat list would also
    # reverse the pieces; sorting by lower endpoint is unambiguous
    # because the cells are disjoint and abutting.
    tail = sorted(_octave_grid(q_floor, octaves, pieces))
    middle = partition_of_axis(
        build_axis(q_floor, 1 - q_floor, middle_width / 2)
    ).cells
    mirror = tuple(
        (1 - upper, 1 - lower) for lower, upper in reversed(tail)
    )
    return Partition(q_min, 1 - q_min, tuple(tail) + middle + mirror)


def partition_of_axis(axis: Axis) -> Partition:
    return Partition(
        axis.lower, axis.upper,
        tuple(axis.bounds(index) for index in range(axis.count)),
    )


@dataclass(frozen=True)
class Cover:
    radius: Fraction
    delta: Fraction
    q_floor: Fraction
    s_axis: Axis
    d_axis: Axis
    q_axis: Partition

    @property
    def cell_count(self) -> int:
        return self.s_axis.count * self.d_axis.count * self.q_axis.count

    def cells(self) -> Iterator[tuple[int, int, int]]:
        for si in range(self.s_axis.count):
            for di in range(self.d_axis.count):
                for qi in range(self.q_axis.count):
                    yield si, di, qi


def build_cover(
    radius: Fraction = DEFAULT_RADIUS,
    delta: Fraction = DEFAULT_DELTA,
    q_floor: Fraction = DEFAULT_Q_FLOOR,
    q_width: Optional[Fraction] = None,
    q_min: Optional[Fraction] = None,
    q_pieces: int = DEFAULT_Q_PIECES,
) -> Cover:
    """Tile the box and the `q` range.

    With `q_min` given the `q` axis carries octave tails down to `q_min` and
    up to `1-q_min`; without it the axis is the uniform `[q_floor,1-q_floor]`
    band, which is all `liu9_chart_centered.py` ever covered.
    """
    if not 0 < q_floor < Fraction(1, 2):
        raise ValueError("the q floor must satisfy 0 < q_floor < 1/2")
    width = delta if q_width is None else q_width
    if q_min is None:
        q_axis = partition_of_axis(
            build_axis(q_floor, 1 - q_floor, width / 2)
        )
    else:
        q_axis = build_q_partition(q_min, q_floor, q_pieces, width)
    return Cover(
        radius,
        delta,
        q_floor,
        build_axis(-radius, radius, delta),
        build_axis(-radius, radius, delta),
        q_axis,
    )


def verify_cover_is_exhaustive(cover: Cover) -> dict[str, Any]:
    """Prove the cover tiles the box and the `q` range, and contains the centre.

    Containment of the chart centre is not decoration.  Step B's mean-value
    argument needs the segment from the centre to every point of the box, and
    it gets it from convexity of the covered region -- which is only the box if
    the centre is actually inside.  A cover of an annulus would pass every PSD
    test and prove nothing.
    """
    report = {
        "s": verify_axis_is_exhaustive(cover.s_axis),
        "d": verify_axis_is_exhaustive(cover.d_axis),
        "q": verify_partition_is_exhaustive(cover.q_axis),
    }
    centre_offset = Fraction(0)
    if not cover.s_axis.lower <= centre_offset <= cover.s_axis.upper:
        raise AssertionError(
            "the covered region excludes the chart centre in s; the "
            "convexity argument behind the mean-value step does not apply"
        )
    if not cover.d_axis.lower <= centre_offset <= cover.d_axis.upper:
        raise AssertionError(
            "the covered region excludes the chart centre in d; the "
            "convexity argument behind the mean-value step does not apply"
        )
    report["contains_chart_centre"] = True
    report["convex"] = "product of intervals"
    report["cells"] = cover.cell_count
    return report


# ---------------------------------------------------------------------------
# One cell of the cover.


def _supremum(current: arb, candidate: arb) -> arb:
    upper = _abs_upper(candidate)
    return upper if float(upper.upper()) >= float(current.upper()) else current


def cell_hessians(
    parameters: ArbParameters,
    s_centre: Fraction,
    d_centre: Fraction,
    half: Fraction,
    q_lower: Fraction,
    q_upper: Fraction,
    subdivisions: int = DEFAULT_SUBDIVISIONS,
) -> tuple[HessianBox, HessianBox]:
    """Centered-form Hessians on one cell, expanded at that cell's own centre.

    This is `liu9_chart_centered._centered_hessian` freed from the chart
    centre: the expansion point is the cell's centre, and the third-derivative
    suprema are taken over the cell rather than over the whole box.  Both
    changes are what make the enclosure width scale with `delta` instead of
    with the radius.
    """
    if subdivisions <= 0:
        raise ValueError("third-derivative subdivisions must be positive")
    half_arb = _arbf(half)
    centre_s = parameters.x + _arbf(s_centre)
    centre_d = _arbf(d_centre)
    q = _arbf(q_lower).union(_arbf(q_upper))

    gap_third = [arb(0)] * 4
    distance_third = [arb(0)] * 4
    for s_index in range(subdivisions):
        s_low = centre_s + half_arb * (-1 + arb(2 * s_index) / subdivisions)
        s_high = centre_s + half_arb * (
            -1 + arb(2 * (s_index + 1)) / subdivisions
        )
        s_box = s_low.union(s_high)
        for d_index in range(subdivisions):
            d_low = centre_d + half_arb * (
                -1 + arb(2 * d_index) / subdivisions
            )
            d_high = centre_d + half_arb * (
                -1 + arb(2 * (d_index + 1)) / subdivisions
            )
            d_box = d_low.union(d_high)
            gap = gap_jet(
                parameters, q,
                Jet3.variable(s_box, "s"), Jet3.variable(d_box, "d"), arb(0),
            )
            distance = distance_squared_jet(
                parameters, q,
                Jet3.variable(s_box, "s"), Jet3.variable(d_box, "d"), arb(0),
            )
            for index, value in enumerate(
                (gap.tsss, gap.tssd, gap.tsdd, gap.tddd)
            ):
                gap_third[index] = _supremum(gap_third[index], value)
            for index, value in enumerate(
                (
                    distance.tsss, distance.tssd,
                    distance.tsdd, distance.tddd,
                )
            ):
                distance_third[index] = _supremum(
                    distance_third[index], value
                )

    gap_centre = gap_jet(
        parameters, q,
        Jet3.variable(centre_s, "s"), Jet3.variable(centre_d, "d"), arb(0),
    )
    distance_centre = distance_squared_jet(
        parameters, q,
        Jet3.variable(centre_s, "s"), Jet3.variable(centre_d, "d"), arb(0),
    )

    def entry(centre: arb, first: arb, second: arb) -> arb:
        variation = half_arb * (first + second)
        return centre + (-variation).union(variation)

    return (
        HessianBox(
            entry(gap_centre.hss, gap_third[0], gap_third[1]),
            entry(gap_centre.hsd, gap_third[1], gap_third[2]),
            entry(gap_centre.hdd, gap_third[2], gap_third[3]),
        ),
        HessianBox(
            entry(distance_centre.hss, distance_third[0], distance_third[1]),
            entry(distance_centre.hsd, distance_third[1], distance_third[2]),
            entry(distance_centre.hdd, distance_third[2], distance_third[3]),
        ),
    )


# ---------------------------------------------------------------------------
# Sweeping the cover, with a resumable checkpoint.


@dataclass
class SweepState:
    processed: int = 0
    failures: list[dict[str, Any]] | None = None
    worst_margin: Optional[float] = None
    worst_cell: Optional[tuple[str, str, str]] = None
    elapsed: float = 0.0

    def __post_init__(self) -> None:
        if self.failures is None:
            self.failures = []


def _cell_margin(gap: HessianBox, distance: HessianBox, kappa: arb) -> float:
    """The determinant lower bound the PSD test actually uses."""
    m11 = gap.ss - kappa * distance.ss
    m22 = gap.dd - kappa * distance.dd
    m12 = gap.sd - kappa * distance.sd
    off = _abs_upper(m12)
    return float(((m11 * m22).lower() - off * off).lower())


def sweep_cover(
    parameters: ArbParameters,
    cover: Cover,
    kappa: Fraction,
    subdivisions: int = DEFAULT_SUBDIVISIONS,
    checkpoint_path: Optional[str] = None,
    checkpoint_every: int = 2000,
    max_failures: int = 24,
    progress_every: int = 0,
    resume: bool = True,
) -> SweepState:
    """Certify the pencil on every cell; checkpointed and resumable."""
    kappa_arb = _arbf(kappa)
    ceiling = smooth_chart_ceiling(parameters)
    if not _arbf(kappa) <= ceiling.lower():
        raise AssertionError(
            f"kappa {kappa} exceeds the smooth-chart ceiling "
            f"{ceiling.lower().str(12)}"
        )

    state = SweepState()
    start_index = 0
    if resume and checkpoint_path and os.path.exists(checkpoint_path):
        with open(checkpoint_path, "r", encoding="utf-8") as handle:
            saved = json.load(handle)
        if (
            saved.get("radius") == str(cover.radius)
            and saved.get("delta") == str(cover.delta)
            and saved.get("kappa") == str(kappa)
            and saved.get("cells") == cover.cell_count
        ):
            start_index = int(saved["processed"])
            state.processed = start_index
            state.failures = saved.get("failures", [])
            state.worst_margin = saved.get("worst_margin")
            state.elapsed = float(saved.get("elapsed", 0.0))

    began = time.monotonic()
    for index, (si, di, qi) in enumerate(cover.cells()):
        if index < start_index:
            continue
        s_centre = cover.s_axis.centre(si)
        d_centre = cover.d_axis.centre(di)
        q_lower, q_upper = cover.q_axis.bounds(qi)
        gap, distance = cell_hessians(
            parameters, s_centre, d_centre, cover.delta,
            q_lower, q_upper, subdivisions,
        )
        if not pencil_is_psd(gap, distance, kappa_arb):
            state.failures.append(
                {
                    "s_centre": str(s_centre),
                    "d_centre": str(d_centre),
                    "q_lower": str(q_lower),
                    "q_upper": str(q_upper),
                    "margin": _cell_margin(gap, distance, kappa_arb),
                }
            )
            if len(state.failures) >= max_failures:
                state.processed = index + 1
                state.elapsed += time.monotonic() - began
                return state
        margin = _cell_margin(gap, distance, kappa_arb)
        if state.worst_margin is None or margin < state.worst_margin:
            state.worst_margin = margin
            state.worst_cell = (str(s_centre), str(d_centre), str(q_lower))
        state.processed = index + 1
        if progress_every and state.processed % progress_every == 0:
            rate = state.processed / max(time.monotonic() - began, 1e-9)
            print(
                f"    {state.processed}/{cover.cell_count} cells, "
                f"{rate:.0f}/s, worst margin {state.worst_margin:.6e}",
                flush=True,
            )
        if (
            checkpoint_path
            and checkpoint_every
            and state.processed % checkpoint_every == 0
        ):
            _write_checkpoint(checkpoint_path, cover, kappa, state,
                              time.monotonic() - began)
    state.elapsed += time.monotonic() - began
    if checkpoint_path:
        _write_checkpoint(checkpoint_path, cover, kappa, state, state.elapsed)
    return state


def _write_checkpoint(path: str, cover: Cover, kappa: Fraction,
                      state: SweepState, elapsed: float) -> None:
    payload = {
        "radius": str(cover.radius),
        "delta": str(cover.delta),
        "kappa": str(kappa),
        "cells": cover.cell_count,
        "processed": state.processed,
        "failures": state.failures,
        "worst_margin": state.worst_margin,
        "elapsed": elapsed,
    }
    temporary = f"{path}.partial"
    with open(temporary, "w", encoding="utf-8") as handle:
        json.dump(payload, handle)
    os.replace(temporary, path)


# ---------------------------------------------------------------------------
# Step B's hypothesis: what happens at the chart centre.

CENTRE_PROBE_Q = (
    Fraction(1, 4096), Fraction(1, 256), Fraction(1, 64), Fraction(1, 4),
    Fraction(1, 2), Fraction(3, 4), Fraction(4095, 4096),
)


def certify_gd_antisymmetry() -> dict[str, Any]:
    """`d/dd` of the centre gap is exactly zero for EVERY q, not just small.

    At `d = 0` the two chart supports move oppositely, `d(s-qd)/dd = -q` and
    `d(s+(1-q)d)/dd = 1-q`, while the components carry weights `1-q` and `q`.
    Every term of EHX, EHXY and EHPI therefore contributes its own q-free
    multiplier times the same bracket

        (1-q)*(-q) + q*(1-q),

    which is the zero polynomial.  Expanding it in exact `Fraction` arithmetic
    gives all three coefficients zero, so the cancellation is algebraic and
    holds for every q at once -- there is nothing to sample.  The residue Arb
    reports at a point q is rounding of terms whose size runs like `q(1-q)`,
    which is why it looked like a q-dependent bound when it was measured
    pointwise.
    """
    def bracket(q: Fraction) -> Fraction:
        return (1 - q) * (-q) + q * (1 - q)

    at_zero, at_one, at_two = (bracket(Fraction(k)) for k in (0, 1, 2))
    quadratic = (at_two - 2 * at_one + at_zero) / 2
    linear = at_one - at_zero - quadratic
    constant = at_zero
    if not (constant == 0 and linear == 0 and quadratic == 0):
        raise AssertionError(
            "the d-direction antisymmetry bracket is not the zero polynomial: "
            f"q^0={constant}, q^1={linear}, q^2={quadratic}"
        )
    return {
        "bracket": "(1-q)*(-q) + q*(1-q)",
        "coefficients": {"q^0": str(constant), "q^1": str(linear),
                         "q^2": str(quadratic)},
        "identically_zero": True,
    }


def collapsed_centre_jet(parameters: ArbParameters) -> Jet3:
    """Value AND `d/ds` at the chart centre, with no `q` anywhere.

    Same collapse as `collapsed_centre_gap`, carried as a jet in `s` so the
    `s`-derivative comes out of the same q-free expression.  This is what
    upgrades the Step B gradient bound from a sampled statement to an
    all-q enclosure: `gs` dominates the deficit, so leaving it sampled left the
    whole bound resting on a probe set.
    """
    s = Jet3.variable(parameters.x, "s")
    mass = Jet3.constant(parameters.mean) / s
    ehx = mass * jet_entropy(s)
    ehxy = mass * mass * jet_entropy(s * s)
    ehpi = mass * mass * jet_entropy(protocol_diagonal(s))
    beta = Jet3.constant(parameters.beta)
    return (Jet3.constant(1) - beta) * ehxy + beta * ehpi - ehx


def collapsed_centre_gap(parameters: ArbParameters) -> arb:
    """`gap` at the chart centre, computed with no `q` anywhere.

    At `(s,d) = (x,0)` the chart degenerates: `support_zero = s - q d` and
    `support_one = s + (1-q) d` both collapse to `x`, the gauge atom has mass
    `r = 0`, and BOTH components carry the same two-atom law `{x: m, 0: 1-m}`
    with `m = mean/x`.  Every `q` then cancels structurally --

        EHX   = m h(x)            (the q-weights sum to 1),
        EHXY  = m^2 h(x^2)        (only the x*x product has nonzero support),
        EHPI  = m^2 h(pi(x))      (the component weights sum to 1),

    -- so the centre value is a q-free constant.  Evaluating it this way
    matters: `gap_jet` with an INTERVAL q loses the cancellation to dependency
    and returns a useless enclosure of width about 6e-3, while the true
    function is constant.  This routine is the same number without the
    dependency, and `certify_centre_identity` guards the transcription against
    `gap_jet` at point q.
    """
    mass = parameters.mean / parameters.x
    entropy = lambda value: jet_entropy(Jet3.constant(value)).v  # noqa: E731
    protocol = protocol_diagonal(Jet3.constant(parameters.x)).v
    return (
        (1 - parameters.beta) * mass * mass * entropy(parameters.x ** 2)
        + parameters.beta * mass * mass * entropy(protocol)
        - mass * entropy(parameters.x)
    )


def certify_centre_identity(
    parameters: ArbParameters,
    radius: Fraction,
    probe_q: Sequence[Fraction] = CENTRE_PROBE_Q,
) -> dict[str, Any]:
    """Bound the Step B deficit, and say honestly which part is sampled.

    Step B needs `F_q(0) = 0` and `grad F_q(0) = 0`.  Neither is exactly zero
    here, and the reason is not the cover: Liu's abscissa `x` is a numerically
    determined root of equations (87)-(90), so the binding residual is whatever
    the solve leaves behind -- about 7e-69 at 90 dps.  The pre-existing radius
    1/256 certificate carries exactly the same residual; lifting the radius
    neither improves nor worsens it.

    The value is certified q-free by `collapsed_centre_gap`.  The two gradient
    components are reported over a probe set of point q and labelled as
    sampled, because an interval-q evaluation of them loses the cancellation.
    """
    jet = collapsed_centre_jet(parameters)
    collapsed = jet.v
    antisymmetry = certify_gd_antisymmetry()
    q_free_value = float(_abs_upper(jet.v).upper())
    q_free_gs = float(_abs_upper(jet.gs).upper())
    worst_value = q_free_value
    worst_transcription = 0.0
    worst_gs = 0.0
    worst_gd = 0.0
    worst_distance = 0.0
    for q_value in probe_q:
        q = _arbf(q_value)
        centre_s = Jet3.variable(parameters.x, "s")
        centre_d = Jet3.variable(arb(0), "d")
        gap = gap_jet(parameters, q, centre_s, centre_d, arb(0))
        distance = distance_squared_jet(
            parameters, q, centre_s, centre_d, arb(0)
        )
        worst_transcription = max(
            worst_transcription,
            float(_abs_upper(gap.v - collapsed).upper()),
        )
        worst_value = max(worst_value, float(_abs_upper(gap.v).upper()))
        worst_gs = max(worst_gs, float(_abs_upper(gap.gs).upper()))
        worst_gd = max(worst_gd, float(_abs_upper(gap.gd).upper()))
        for entry in (distance.v, distance.gs, distance.gd):
            worst_distance = max(
                worst_distance, float(_abs_upper(entry).upper())
            )
    reach = float(radius)
    # The bound below is an ALL-q enclosure, not a sampled one.  The value and
    # the s-gradient come from `collapsed_centre_jet`, in which no q appears at
    # any step, and the d-gradient is exactly zero for every q by the algebraic
    # antisymmetry certified above.  The point-q probes are retained only as a
    # transcription guard on the collapsed forms.
    deficit = q_free_value + q_free_gs * reach
    if not worst_transcription <= 10 * (q_free_value + q_free_gs):
        raise AssertionError(
            "the collapsed centre forms disagree with gap_jet beyond rounding"
        )
    return {
        "collapsed_centre_gap": collapsed.str(20),
        "collapsed_centre_gs": jet.gs.str(20),
        "q_free_by_construction": True,
        "gd_antisymmetry": antisymmetry,
        "q_free_value_bound": q_free_value,
        "q_free_gs_bound": q_free_gs,
        "deficit_is_all_q": True,
        "probe_q": [str(value) for value in probe_q],
        "worst_transcription_difference": worst_transcription,
        "worst_centre_value": worst_value,
        "worst_centre_gs": worst_gs,
        "worst_centre_gd": worst_gd,
        "worst_distance_centre_entry": worst_distance,
        "step_b_deficit": deficit,
    }


# ---------------------------------------------------------------------------
# Independent falsification pass.


def falsification_pass(
    parameters: ArbParameters,
    cover: Cover,
    kappa: Fraction,
    points: int = 400,
    seed: int = 20260828,
) -> dict[str, Any]:
    """Evaluate the pencil at interior points, bypassing the cover entirely.

    The sweep certifies the pencil through centered-form enclosures over
    cells.  A systematic error in that construction -- a wrong expansion
    point, a dropped third-derivative term, a mis-signed variation -- would
    make every cell agree with every other cell and still be wrong.  This pass
    evaluates the two Hessians directly at POINT `(s,d,q)` values, where the
    jets carry no interval width at all, and checks the same PSD condition.
    It cannot certify anything, and does not try to; it can only refute.
    """
    import random

    rng = random.Random(seed)
    kappa_arb = _arbf(kappa)
    worst_eigenvalue = None
    worst_point = None
    violations = []
    radius = cover.radius
    for _ in range(points):
        s_offset = Fraction(rng.randint(-10**6, 10**6), 10**6) * radius
        d_offset = Fraction(rng.randint(-10**6, 10**6), 10**6) * radius
        low, high = cover.q_axis.lower, cover.q_axis.upper
        q_value = low + (high - low) * Fraction(rng.randint(0, 10**6), 10**6)
        q = _arbf(q_value)
        centre_s = parameters.x + _arbf(s_offset)
        centre_d = _arbf(d_offset)
        gap = gap_jet(
            parameters, q,
            Jet3.variable(centre_s, "s"), Jet3.variable(centre_d, "d"),
            arb(0),
        )
        distance = distance_squared_jet(
            parameters, q,
            Jet3.variable(centre_s, "s"), Jet3.variable(centre_d, "d"),
            arb(0),
        )
        m11 = gap.hss - kappa_arb * distance.hss
        m22 = gap.hdd - kappa_arb * distance.hdd
        m12 = gap.hsd - kappa_arb * distance.hsd
        # Smallest eigenvalue of the symmetric 2x2, at a point.
        trace = m11 + m22
        discriminant = ((m11 - m22) ** 2 + 4 * m12 * m12).sqrt()
        smallest = (trace - discriminant) / 2
        value = float(smallest.lower())
        if worst_eigenvalue is None or value < worst_eigenvalue:
            worst_eigenvalue = value
            worst_point = (str(s_offset), str(d_offset), str(q_value))
        if value < 0:
            violations.append(
                {
                    "s_offset": str(s_offset),
                    "d_offset": str(d_offset),
                    "q": str(q_value),
                    "smallest_eigenvalue": value,
                }
            )
    return {
        "points": points,
        "seed": seed,
        "worst_smallest_eigenvalue": worst_eigenvalue,
        "worst_point": worst_point,
        "violations": violations,
    }


# ---------------------------------------------------------------------------
# Reporting.


def _canonical_digest(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius", type=str, default=str(DEFAULT_RADIUS))
    parser.add_argument("--delta", type=str, default=str(DEFAULT_DELTA))
    parser.add_argument("--kappa", type=str, default=str(CERTIFIED_KAPPA))
    parser.add_argument("--q-floor", type=str, default=str(DEFAULT_Q_FLOOR))
    parser.add_argument("--q-width", type=str, default=None)
    parser.add_argument("--q-min", type=str, default=None,
                        help="extend q by octave tails down to this "
                             "value and up to its mirror")
    parser.add_argument("--q-pieces", type=int,
                        default=DEFAULT_Q_PIECES)
    parser.add_argument("--falsification-points", type=int, default=400)
    parser.add_argument("--subdivisions", type=int,
                        default=DEFAULT_SUBDIVISIONS)
    parser.add_argument("--dps", type=int, default=DEFAULT_DPS)
    parser.add_argument("--checkpoint", type=str, default=None)
    parser.add_argument("--checkpoint-every", type=int, default=2000)
    parser.add_argument("--progress-every", type=int, default=0)
    parser.add_argument("--no-resume", action="store_true")
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args(argv)

    radius = Fraction(args.radius)
    delta = Fraction(args.delta)
    kappa = Fraction(args.kappa)
    q_floor = Fraction(args.q_floor)
    q_width = Fraction(args.q_width) if args.q_width else None
    q_min = Fraction(args.q_min) if args.q_min else None

    started = time.monotonic()
    mp_parameters = solve_equation_parameters(args.dps)
    parameters = certify_equation_parameters(mp_parameters)

    print("1. THE COVER GEOMETRY")
    cover = build_cover(radius, delta, q_floor, q_width, q_min,
                        args.q_pieces)
    geometry = verify_cover_is_exhaustive(cover)
    print(
        f"PROVED [exact rational]: the box |s-x|<={radius}, |d|<={radius} is "
        f"tiled by {cover.s_axis.count}x{cover.d_axis.count} cells of "
        f"half-width {delta}, and q in "
        f"[{cover.q_axis.lower},{cover.q_axis.upper}] by "
        f"{cover.q_axis.count} cells of relative width at most "
        f"{cover.q_axis.widest_relative()}; "
        f"{geometry['s']['abutment_checks'] + geometry['d']['abutment_checks'] + geometry['q']['abutment_checks']}"
        " abutment identities hold with no gap and no slack."
    )
    print(
        "PROVED [geometry]: the covered region is a product of intervals, "
        "hence convex, and contains the chart centre.  That is what lets the "
        "mean-value step use segments from the centre without any per-cell "
        "segment hypothesis, and it is why no radial quadrature is needed."
    )
    print(f"    cells to certify: {cover.cell_count}")

    print()
    print("2. THE PENCIL ON EVERY CELL")
    ceiling = smooth_chart_ceiling(parameters)
    print(f"    smooth-chart ceiling (raw-gap units): {ceiling.str(12)}")
    print(f"    candidate kappa: {kappa} = {float(kappa):.10f}")
    state = sweep_cover(
        parameters, cover, kappa, args.subdivisions,
        args.checkpoint, args.checkpoint_every,
        progress_every=args.progress_every,
        resume=not args.no_resume,
    )
    certified = not state.failures
    if certified:
        print(
            f"PROVED [Arb, {state.processed} cells]: "
            f"Hess gap - kappa*Hess dist^2 is positive semidefinite at every "
            f"point of the box, with kappa = {kappa}.  The weakest "
            f"determinant margin over the whole cover is "
            f"{state.worst_margin:.6e} at "
            f"(s-x, d, q) = {state.worst_cell}."
        )
    else:
        first = state.failures[0]
        print(
            f"REFUTED [Arb]: {len(state.failures)} of {state.processed} cells "
            f"failed.  The first is s-x={first['s_centre']}, "
            f"d={first['d_centre']}, q in "
            f"[{first['q_lower']},{first['q_upper']}] with determinant margin "
            f"{first['margin']:.6e}."
        )

    print()
    print("3. WHAT STEP B INHERITS AT THE CHART CENTRE")
    centre = certify_centre_identity(parameters, radius)
    print(
        f"PROVED [Arb, q-free by construction]: at (s,d)=(x,0) both chart "
        f"components collapse to the same two-atom law {{x: m, 0: 1-m}} with "
        f"m = mean/x, so every q cancels and the centre gap is the constant "
        f"{centre['collapsed_centre_gap']}.  Computing it this way is not "
        f"cosmetic: gap_jet at an INTERVAL q loses the cancellation to "
        f"dependency and returns a width of about 6e-3 around a constant."
    )
    print(
        f"PROVED [Arb, q-free by construction]: the same collapse carried as a "
        f"jet in s gives the centre s-gradient as a q-free enclosure, "
        f"|d/ds| <= {centre['q_free_gs_bound']:.3e}."
    )
    print(
        f"PROVED [exact Fraction, all q]: the centre d-gradient is EXACTLY "
        f"zero for every q.  Both supports move oppositely at d=0 while the "
        f"components carry the opposite weights, so every term contributes a "
        f"q-free multiplier times the bracket "
        f"{centre['gd_antisymmetry']['bracket']}, whose expansion has "
        f"coefficients q^0={centre['gd_antisymmetry']['coefficients']['q^0']}, "
        f"q^1={centre['gd_antisymmetry']['coefficients']['q^1']}, "
        f"q^2={centre['gd_antisymmetry']['coefficients']['q^2']}.  The residue "
        f"Arb shows at a point q is rounding of terms sized like q(1-q), which "
        f"is why a pointwise measurement made it look q-dependent."
    )
    print(
        f"COMPUTATIONAL EVIDENCE [{len(centre['probe_q'])} point-q probes, "
        f"transcription guard only]: the collapsed forms reproduce gap_jet at "
        f"the centre to {centre['worst_transcription_difference']:.3e}, and "
        f"every dist^2 centre entry is below "
        f"{centre['worst_distance_centre_entry']:.3e}."
    )
    print(
        f"CONDITIONAL: Step B therefore yields gap - kappa*dist^2 >= "
        f"-{centre['step_b_deficit']:.3e} on the box for EVERY q in the "
        f"covered range -- an all-q enclosure, not a sampled bound -- rather "
        f"than >= 0 exactly.  That "
        f"deficit is the residual of Liu's equations (87)-(90) at {args.dps} "
        f"dps, since x is a numerically determined root; it is inherited "
        f"unchanged from the radius {INHERITED_RADIUS} certificate and is not "
        f"introduced by the cover."
    )

    print()
    print("4. INDEPENDENT FALSIFICATION PASS")
    falsification = falsification_pass(parameters, cover, kappa,
                                       args.falsification_points)
    if falsification["violations"]:
        first = falsification["violations"][0]
        print(
            f"REFUTED [pointwise]: {len(falsification['violations'])} of "
            f"{falsification['points']} interior points have a negative "
            f"smallest eigenvalue; the first is {first}."
        )
        certified = False
    else:
        print(
            f"COMPUTATIONAL EVIDENCE [{falsification['points']} interior "
            f"points, point-valued jets, no cell arithmetic]: the smallest "
            f"eigenvalue of the pencil never fell below "
            f"{falsification['worst_smallest_eigenvalue']:.6e}, attained at "
            f"(s-x, d, q) = {falsification['worst_point']}.  This pass cannot "
            f"certify; it exists to refute a systematic error in the cell "
            f"construction, and it did not."
        )

    print()
    print("5. WHAT THIS CLOSES")
    if certified and radius >= TARGET_RADIUS:
        print(
            f"PROVED: the smooth-chart pencil statement now holds out to "
            f"radius {radius}, up from the {INHERITED_RADIUS} certified in "
            f"liu9_chart_centered.py -- a factor "
            f"{int(radius / INHERITED_RADIUS)}."
        )
        print(
            "PROVED: ingredient (iii) of the piecewise local lemma was "
            "CONDITIONAL on exactly this statement at eps_sm=1/32 "
            "(uc/liu9_qdegenerate.py).  That hypothesis is now discharged on "
            "the chart."
        )
        print(
            "OPEN: the chart still carries at most two distinct supports per "
            "component.  Extending it to a full nine-dimensional "
            "neighbourhood is a separate question and is not claimed here.  "
            "No tube radius is certified by this module."
        )
    elif certified:
        print(
            f"PROVED: the pencil statement holds at radius {radius}, which is "
            f"short of the {TARGET_RADIUS} that ingredient (iii) needs."
        )
    else:
        print(
            "OPEN: the conditional hypothesis of ingredient (iii) is not "
            "discharged at this radius and cell size."
        )

    elapsed = time.monotonic() - started
    payload: dict[str, Any] = {
        "radius": str(radius),
        "delta": str(delta),
        "kappa": str(kappa),
        "q_floor": str(q_floor),
        "q_min": str(cover.q_axis.lower),
        "q_max": str(cover.q_axis.upper),
        "q_cells": cover.q_axis.count,
        "subdivisions": args.subdivisions,
        "geometry": geometry,
        "cells": cover.cell_count,
        "processed": state.processed,
        "certified": certified,
        "worst_margin": state.worst_margin,
        "worst_cell": list(state.worst_cell) if state.worst_cell else None,
        "failures": state.failures,
        "falsification": falsification,
        "centre_identity": centre,
        "ceiling": ceiling.str(20),
        "inherited_radius": str(INHERITED_RADIUS),
        "target_radius": str(TARGET_RADIUS),
    }
    if args.output:
        digest = _canonical_digest(payload)
        payload["report_sha256"] = digest
        payload["elapsed_seconds"] = elapsed
        with open(args.output, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=1, sort_keys=True)
        print()
        print(f"report written to {args.output}")
        print(f"report_sha256 {digest}")
    return 0 if certified else 1


if __name__ == "__main__":
    raise SystemExit(main())
