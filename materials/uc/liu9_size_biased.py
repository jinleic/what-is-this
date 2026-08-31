#!/usr/bin/env python3
"""Two-variable size-biased kernel certificate for Liu's q-one raw gap."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from fractions import Fraction
from pathlib import Path
from typing import Optional, Sequence

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

from flint import arb, ctx

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from liu9_binding import certify_equation_parameters, solve_equation_parameters
from liu9_objective import _ArbOps, arb_fraction

from liu9_pilot import GradientUnavailable, Jet, _JetOps
from liu9_chart_centered import Jet3, compose_unary, jet_entropy
ctx.prec = max(ctx.prec, 320)
_OPS = _ArbOps(monotone_corners=True)
OUTPUT_DEFAULT = (
    HERE / "verification/results/liu9-size-biased-qone-kernel.json")


def _canonical_digest(payload: dict) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode()).hexdigest()




def _mu_point(value: arb) -> arb:
    if value == 0:
        return arb(1)
    if value == 1:
        return arb(0)
    return -((1 - value) * (1 - value).log()) / value


def _mu_range(value: arb) -> arb:
    return _mu_point(arb(value.upper())).union(_mu_point(arb(value.lower())))


def direct_kernel_range(x: arb, y: arb, mean: arb, beta: arb) -> arb:
    hx, hy = _OPS.entropy(x), _OPS.entropy(y)
    kernel = (
        (1 - beta) * _OPS.entropy(_OPS.xy_argument(x, y))
        + beta * _OPS.entropy(_OPS.pi_argument(x, y))
    )
    return mean * kernel - (y * hx + x * hy) / 2


def centered_direct_kernel_range(
    x: arb, y: arb, mean: arb, beta: arb,
) -> arb:
    ops = _JetOps(2)
    x_jet = Jet(x, (arb(1), arb(0)))
    y_jet = Jet(y, (arb(0), arb(1)))
    kernel = (
        (1 - beta) * ops.entropy(ops.xy_argument(x_jet, y_jet))
        + beta * ops.entropy(ops.pi_argument(x_jet, y_jet))
    )
    value = mean * kernel - (
        y_jet * ops.entropy(x_jet) + x_jet * ops.entropy(y_jet)
    ) / 2
    x_mid = (x.lower() + x.upper()) / 2
    y_mid = (y.lower() + y.upper()) / 2
    center = direct_kernel_range(arb(x_mid), arb(y_mid), mean, beta)
    dx = arb(x.lower() - x_mid).union(arb(x.upper() - x_mid))
    dy = arb(y.lower() - y_mid).union(arb(y.upper() - y_mid))
    return center + value.gradient[0] * dx + value.gradient[1] * dy


def kernel_range(x: arb, y: arb, mean: arb, beta: arb) -> tuple[arb, arb]:
    hx, hy = _OPS.entropy(x), _OPS.entropy(y)
    product = _OPS.xy_argument(x, y)
    protocol = _OPS.pi_argument(x, y)
    t_kernel = _mu_range(x) + _mu_range(y) - _mu_range(product)
    protocol_gain = _OPS.entropy(protocol) - _OPS.entropy(product)
    entropy_term = y * hx + x * hy
    remainder = -x * y * t_kernel + beta * protocol_gain
    derivative_mean = entropy_term + remainder
    kernel = mean * derivative_mean - entropy_term / 2
    return kernel, derivative_mean
def _mu_jet(value: Jet) -> Jet:
    if not (value.value.lower() > 0 and value.value.upper() < 1):
        raise GradientUnavailable("mu derivative needs an interior interval")
    derivative = (
        value.value + (1 - value.value).log()
    ) / value.value**2
    return Jet(
        _mu_range(value.value),
        tuple(derivative * entry for entry in value.gradient),
    )


def _kernel_jet(x: Jet, y: Jet, mean: arb, beta: arb) -> Jet:
    ops = _JetOps(2)
    hx, hy = ops.entropy(x), ops.entropy(y)
    product = ops.xy_argument(x, y)
    protocol = ops.pi_argument(x, y)
    t_kernel = _mu_jet(x) + _mu_jet(y) - _mu_jet(product)
    protocol_gain = ops.entropy(protocol) - ops.entropy(product)
    entropy_term = y * hx + x * hy
    derivative_mean = entropy_term - x * y * t_kernel + beta * protocol_gain
    return mean * derivative_mean - entropy_term / 2


def centered_kernel_range(x: arb, y: arb, mean: arb, beta: arb) -> arb:
    x_jet = Jet(x, (arb(1), arb(0)))
    y_jet = Jet(y, (arb(0), arb(1)))
    gradient = _kernel_jet(x_jet, y_jet, mean, beta).gradient
    x_mid = (x.lower() + x.upper()) / 2
    y_mid = (y.lower() + y.upper()) / 2
    center, _ = kernel_range(arb(x_mid), arb(y_mid), mean, beta)
    dx = arb(x.lower() - x_mid).union(arb(x.upper() - x_mid))
    dy = arb(y.lower() - y_mid).union(arb(y.upper() - y_mid))
    return center + gradient[0] * dx + gradient[1] * dy




def _mu_jet3(argument: Jet3) -> Jet3:
    x = argument.v
    if not (x.lower() > 0 and x.upper() < 1):
        raise ArithmeticError("mu Jet3 needs an interior interval")
    logarithm = (1 - x).log()
    value = _mu_range(x)
    first = (x + logarithm) / x**2
    second = (
        -x**2 - 2 * x * logarithm + 2 * x + 2 * logarithm
    ) / (x**3 * (x - 1))
    third = (
        2 * x**3 + 6 * x**2 * logarithm - 9 * x**2
        - 12 * x * logarithm + 6 * x + 6 * logarithm
    ) / (x**4 * (1 - x)**2)
    return compose_unary(argument, value, first, second, third)


def kernel_jet3(x_value: arb, y_value: arb, mean: arb, beta: arb) -> Jet3:
    x = Jet3.variable(x_value, "s")
    y = Jet3.variable(y_value, "d")
    product = x * y
    protocol = product * (1 + (1 - x) * (1 - y))
    t_kernel = _mu_jet3(x) + _mu_jet3(y) - _mu_jet3(product)
    protocol_gain = jet_entropy(protocol) - jet_entropy(product)
    entropy_term = y * jet_entropy(x) + x * jet_entropy(y)
    derivative_mean = entropy_term - x * y * t_kernel + beta * protocol_gain
    return mean * derivative_mean - entropy_term / 2


def direct_kernel_jet3(
    x_value: arb, y_value: arb, mean: arb, beta: arb,
) -> Jet3:
    x = Jet3.variable(x_value, "s")
    y = Jet3.variable(y_value, "d")
    product = x * y
    protocol = product * (1 + (1 - x) * (1 - y))
    kernel = (
        (1 - beta) * jet_entropy(product) + beta * jet_entropy(protocol))
    return mean * kernel - (
        y * jet_entropy(x) + x * jet_entropy(y)) / 2


def certify_direct_optimizer_box(subdivisions: int = 32) -> dict:
    parameters = certify_equation_parameters(solve_equation_parameters(100))
    lo, hi = Fraction(17, 25), Fraction(7, 10)
    radius = max(
        parameters.x.upper() - arb_fraction(lo),
        arb_fraction(hi) - parameters.x.lower(),
        key=float,
    )
    names = ("tsss", "tssd", "tsdd", "tddd")
    suprema = {name: arb(0) for name in names}
    for i in range(subdivisions):
        x = (
            arb_fraction(lo + (hi - lo) * i / subdivisions)
            .union(arb_fraction(lo + (hi - lo) * (i + 1) / subdivisions))
        )
        for j in range(subdivisions):
            y = (
                arb_fraction(lo + (hi - lo) * j / subdivisions)
                .union(arb_fraction(
                    lo + (hi - lo) * (j + 1) / subdivisions))
            )
            jet = direct_kernel_jet3(
                x, y, parameters.mean, parameters.beta)
            for name in names:
                value = getattr(jet, name)
                upper = max(abs(value.lower()), abs(value.upper()))
                if upper > suprema[name]:
                    suprema[name] = arb(upper)
    center = direct_kernel_jet3(
        parameters.x, parameters.x, parameters.mean, parameters.beta)
    vxx = radius * (suprema["tsss"] + suprema["tssd"])
    vxy = radius * (suprema["tssd"] + suprema["tsdd"])
    vyy = radius * (suprema["tsdd"] + suprema["tddd"])
    hxx = center.hss + (-vxx).union(vxx)
    hxy = center.hsd + (-vxy).union(vxy)
    hyy = center.hdd + (-vyy).union(vyy)
    offdiag = max(abs(hxy.lower()), abs(hxy.upper()))
    determinant = hxx.lower() * hyy.lower() - offdiag**2
    if not hxx > 0 or not determinant > 0:
        raise AssertionError("direct optimizer Hessian is not positive definite")
    return {
        "determinant_lower": float(determinant.lower()),
        "hxx": str(hxx), "hxy": str(hxy), "hyy": str(hyy),
        "center_value": str(center.v),
        "center_gradient": [str(center.gs), str(center.gd)],
    }


def local_hessian(lo: Fraction, hi: Fraction) -> tuple[arb, arb, arb]:
    parameters = certify_equation_parameters(solve_equation_parameters(100))
    box = arb_fraction(lo).union(arb_fraction(hi))
    jet = kernel_jet3(box, box, parameters.mean, parameters.beta)
    return jet.hss, jet.hsd, jet.hdd


def _cell(index: int, cells: int) -> arb:
    lo = arb_fraction(Fraction(index, cells))
    hi = arb_fraction(Fraction(index + 1, cells))
    return lo.union(hi)


def scan_direct(cells: int) -> dict:
    parameters = certify_equation_parameters(solve_equation_parameters(100))
    unresolved = []
    for i in range(cells):
        x = _cell(i, cells)
        for j in range(i, cells):
            y = _cell(j, cells)
            value = direct_kernel_range(
                x, y, parameters.mean, parameters.beta)
            try:
                centered = centered_direct_kernel_range(
                    x, y, parameters.mean, parameters.beta)
                if centered.lower() > value.lower():
                    value = centered
            except GradientUnavailable:
                pass
            if not value >= 0:
                unresolved.append((
                    i, j, float(value.lower()), float(value.upper())))
    return {
        "cells_per_axis": cells,
        "boxes": cells * (cells + 1) // 2,
        "unresolved_count": len(unresolved),
        "unresolved_pairs": [(item[0], item[1]) for item in unresolved],
        "unresolved": unresolved[:100],
    }


def scan(cells: int) -> dict:
    parameters = certify_equation_parameters(solve_equation_parameters(100))
    unresolved = []
    derivative_unresolved = []
    minimum_lower = None
    for i in range(cells):
        x = _cell(i, cells)
        for j in range(i, cells):
            y = _cell(j, cells)
            value, derivative = kernel_range(
                x, y, parameters.mean, parameters.beta)
            if not value >= 0:
                try:
                    centered = centered_kernel_range(
                        x, y, parameters.mean, parameters.beta)
                    if centered.lower() > value.lower():
                        value = centered
                except GradientUnavailable:
                    pass
            lower = value.lower()
            if minimum_lower is None or lower < minimum_lower:
                minimum_lower = lower
            if not value >= 0:
                unresolved.append(
                    (i, j, float(lower), float(value.upper())))
            if not derivative >= 0:
                derivative_unresolved.append(
                    (i, j, float(derivative.lower()), float(derivative.upper())))
    return {
        "cells_per_axis": cells,
        "boxes": cells * (cells + 1) // 2,
        "minimum_lower": float(minimum_lower),
        "unresolved_count": len(unresolved),
        "derivative_unresolved_count": len(derivative_unresolved),
        "unresolved_by_first_index": {
            str(index): sum(item[0] == index for item in unresolved)
            for index in sorted({item[0] for item in unresolved})
        },
        "unresolved_index_bounds": (
            [min(item[0] for item in unresolved),
             max(item[0] for item in unresolved),
             min(item[1] for item in unresolved),
             max(item[1] for item in unresolved)]
            if unresolved else None
        ),
        "unresolved_pairs": [(item[0], item[1]) for item in unresolved],
        "unresolved": unresolved[:100],
        "derivative_unresolved": derivative_unresolved[:100],
    }


def scan_direct_near_zero(y_cells: int, first_power: int,
                          last_power: int) -> dict:
    parameters = certify_equation_parameters(solve_equation_parameters(100))
    unresolved = []
    for power in range(first_power, last_power + 1):
        x = arb_fraction(Fraction(1, 2 ** (power + 1))).union(
            arb_fraction(Fraction(1, 2**power)))
        for j in range(y_cells):
            y = _cell(j, y_cells)
            value = direct_kernel_range(
                x, y, parameters.mean, parameters.beta)
            try:
                centered = centered_direct_kernel_range(
                    x, y, parameters.mean, parameters.beta)
                if centered.lower() > value.lower():
                    value = centered
            except GradientUnavailable:
                pass
            if not value >= 0:
                unresolved.append((
                    power, j, float(value.lower()), float(value.upper())))
    return {
        "checked": (last_power - first_power + 1) * y_cells,
        "unresolved_count": len(unresolved),
        "unresolved": unresolved,
    }


def scan_direct_zero_corner(first_power: int, last_power: int) -> dict:
    parameters = certify_equation_parameters(solve_equation_parameters(100))
    unresolved = []
    for x_power in range(first_power, last_power + 1):
        x = arb_fraction(Fraction(1, 2 ** (x_power + 1))).union(
            arb_fraction(Fraction(1, 2**x_power)))
        for y_power in range(first_power, x_power + 1):
            y = arb_fraction(Fraction(1, 2 ** (y_power + 1))).union(
                arb_fraction(Fraction(1, 2**y_power)))
            value = direct_kernel_range(
                x, y, parameters.mean, parameters.beta)
            try:
                centered = centered_direct_kernel_range(
                    x, y, parameters.mean, parameters.beta)
                if centered.lower() > value.lower():
                    value = centered
            except GradientUnavailable:
                pass
            if not value >= 0:
                unresolved.append((
                    x_power, y_power,
                    float(value.lower()), float(value.upper())))
    return {
        "unresolved_count": len(unresolved),
        "unresolved": unresolved,
    }


def scan_direct_near_one(x_cells: int, first_power: int,
                         last_power: int) -> dict:
    parameters = certify_equation_parameters(solve_equation_parameters(100))
    unresolved = []
    for power in range(first_power, last_power + 1):
        y = (
            arb(1) - arb_fraction(Fraction(1, 2**power))
        ).union(
            arb(1) - arb_fraction(Fraction(1, 2 ** (power + 1)))
        )
        for i in range(x_cells):
            x = _cell(i, x_cells)
            value = direct_kernel_range(
                x, y, parameters.mean, parameters.beta)
            try:
                centered = centered_direct_kernel_range(
                    x, y, parameters.mean, parameters.beta)
                if centered.lower() > value.lower():
                    value = centered
            except GradientUnavailable:
                pass
            if not value >= 0:
                unresolved.append(
                    (power, i, float(value.lower()), float(value.upper())))
    return {"unresolved_count": len(unresolved), "unresolved": unresolved}


def scan_direct_one_corner(first_power: int, last_power: int) -> dict:
    parameters = certify_equation_parameters(solve_equation_parameters(100))
    unresolved = []
    for x_power in range(first_power, last_power + 1):
        x = (
            arb(1) - arb_fraction(Fraction(1, 2**x_power))
        ).union(
            arb(1) - arb_fraction(Fraction(1, 2 ** (x_power + 1)))
        )
        for y_power in range(x_power, last_power + 1):
            y = (
                arb(1) - arb_fraction(Fraction(1, 2**y_power))
            ).union(
                arb(1) - arb_fraction(Fraction(1, 2 ** (y_power + 1)))
            )
            value = direct_kernel_range(
                x, y, parameters.mean, parameters.beta)
            try:
                centered = centered_direct_kernel_range(
                    x, y, parameters.mean, parameters.beta)
                if centered.lower() > value.lower():
                    value = centered
            except GradientUnavailable:
                pass
            if not value >= 0:
                unresolved.append((
                    x_power, y_power,
                    float(value.lower()), float(value.upper())))
    return {"unresolved_count": len(unresolved), "unresolved": unresolved}


def scan_direct_cross_corner(first_power: int, last_power: int) -> dict:
    parameters = certify_equation_parameters(solve_equation_parameters(100))
    unresolved = []
    for x_power in range(first_power, last_power + 1):
        x = arb_fraction(Fraction(1, 2 ** (x_power + 1))).union(
            arb_fraction(Fraction(1, 2**x_power)))
        for y_power in range(first_power, last_power + 1):
            y = (
                arb(1) - arb_fraction(Fraction(1, 2**y_power))
            ).union(
                arb(1) - arb_fraction(Fraction(1, 2 ** (y_power + 1)))
            )
            value = direct_kernel_range(
                x, y, parameters.mean, parameters.beta)
            try:
                centered = centered_direct_kernel_range(
                    x, y, parameters.mean, parameters.beta)
                if centered.lower() > value.lower():
                    value = centered
            except GradientUnavailable:
                pass
            if not value >= 0:
                unresolved.append((
                    x_power, y_power,
                    float(value.lower()), float(value.upper())))
    return {"unresolved_count": len(unresolved), "unresolved": unresolved}


def scan_near_one(x_cells: int, first_power: int,
                  last_power: int) -> dict:
    parameters = certify_equation_parameters(solve_equation_parameters(100))
    unresolved = []
    checked = 0
    for power in range(first_power, last_power + 1):
        y = (
            arb(1) - arb_fraction(Fraction(1, 2**power))
        ).union(
            arb(1) - arb_fraction(Fraction(1, 2 ** (power + 1)))
        )
        for i in range(x_cells):
            x = _cell(i, x_cells)
            value, _ = kernel_range(x, y, parameters.mean, parameters.beta)
            try:
                centered = centered_kernel_range(
                    x, y, parameters.mean, parameters.beta)
                if centered.lower() > value.lower():
                    value = centered
            except GradientUnavailable:
                pass
            checked += 1
            if not value >= 0:
                unresolved.append(
                    (power, i, float(value.lower()), float(value.upper())))
    return {
        "checked": checked,
        "unresolved_count": len(unresolved),
        "unresolved": unresolved,
    }


def scan_near_one_corner(first_power: int, last_power: int) -> dict:
    parameters = certify_equation_parameters(solve_equation_parameters(100))
    unresolved = []
    for x_power in range(first_power, last_power + 1):
        x = (
            arb(1) - arb_fraction(Fraction(1, 2**x_power))
        ).union(
            arb(1) - arb_fraction(Fraction(1, 2 ** (x_power + 1)))
        )
        for y_power in range(x_power, last_power + 1):
            y = (
                arb(1) - arb_fraction(Fraction(1, 2**y_power))
            ).union(
                arb(1) - arb_fraction(Fraction(1, 2 ** (y_power + 1)))
            )
            value, _ = kernel_range(x, y, parameters.mean, parameters.beta)
            try:
                centered = centered_kernel_range(
                    x, y, parameters.mean, parameters.beta)
                if centered.lower() > value.lower():
                    value = centered
            except GradientUnavailable:
                pass
            if not value >= 0:
                unresolved.append((
                    x_power, y_power,
                    float(value.lower()), float(value.upper())))
    return {
        "checked": (
            (last_power - first_power + 1)
            * (last_power - first_power + 2) // 2
        ),
        "unresolved_count": len(unresolved),
        "unresolved": unresolved,
    }


def certify_optimizer_box(subdivisions: int = 32) -> dict:
    parameters = certify_equation_parameters(solve_equation_parameters(100))
    lo, hi = Fraction(17, 25), Fraction(7, 10)
    radius = max(
        parameters.x.upper() - arb_fraction(lo),
        arb_fraction(hi) - parameters.x.lower(),
        key=float,
    )
    third_names = ("tsss", "tssd", "tsdd", "tddd")
    suprema = {name: arb(0) for name in third_names}
    for i in range(subdivisions):
        x = (
            arb_fraction(lo + (hi - lo) * i / subdivisions)
            .union(arb_fraction(lo + (hi - lo) * (i + 1) / subdivisions))
        )
        for j in range(subdivisions):
            y = (
                arb_fraction(lo + (hi - lo) * j / subdivisions)
                .union(arb_fraction(
                    lo + (hi - lo) * (j + 1) / subdivisions))
            )
            jet = kernel_jet3(x, y, parameters.mean, parameters.beta)
            for name in third_names:
                value = getattr(jet, name)
                upper = max(abs(value.lower()), abs(value.upper()))
                if upper > suprema[name]:
                    suprema[name] = arb(upper)
    center = kernel_jet3(
        parameters.x, parameters.x, parameters.mean, parameters.beta)
    variation_xx = radius * (suprema["tsss"] + suprema["tssd"])
    variation_xy = radius * (suprema["tssd"] + suprema["tsdd"])
    variation_yy = radius * (suprema["tsdd"] + suprema["tddd"])
    hxx = center.hss + (-variation_xx).union(variation_xx)
    hxy = center.hsd + (-variation_xy).union(variation_xy)
    hyy = center.hdd + (-variation_yy).union(variation_yy)
    offdiag = max(abs(hxy.lower()), abs(hxy.upper()))
    determinant = hxx.lower() * hyy.lower() - offdiag**2
    if not hxx > 0 or not determinant > 0:
        raise AssertionError("optimizer-box Hessian is not positive definite")
    return {
        "box": [str(lo), str(hi)],
        "subdivisions": subdivisions,
        "center_value": str(center.v),
        "center_gradient": [str(center.gs), str(center.gd)],
        "hxx": str(hxx),
        "hxy": str(hxy),
        "hyy": str(hyy),
        "determinant_lower": float(determinant.lower()),
        "third_suprema": {
            name: float(value.upper()) for name, value in suprema.items()
        },
    }


def certify_qone_kernel() -> dict:
    parameters = certify_equation_parameters(solve_equation_parameters(100))
    cells = 1024
    bulk = scan(cells)
    optimizer_lo, optimizer_hi = Fraction(17, 25), Fraction(7, 10)
    bulk_discharged = 0
    for i, j in bulk["unresolved_pairs"]:
        small = i < cells // 32
        mirror = j == cells - 1
        local = (
            Fraction(i, cells) >= optimizer_lo
            and Fraction(i + 1, cells) <= optimizer_hi
            and Fraction(j, cells) >= optimizer_lo
            and Fraction(j + 1, cells) <= optimizer_hi
        )
        if not (small or mirror or local):
            raise AssertionError(f"uncovered bulk cell {(i, j)}")
        bulk_discharged += 1
    if bulk_discharged != bulk["unresolved_count"]:
        raise AssertionError("bulk unresolved cells were not all discharged")

    x0 = arb_fraction(Fraction(1, 32))
    small_margin = (
        (parameters.mean - arb(1) / 2)
        * _OPS.entropy(x0) / x0 - arb(1) / 2
    )
    if not small_margin > 0:
        raise AssertionError("small-support analytic margin is not positive")

    optimizer = certify_optimizer_box()
    mirror = scan_near_one(cells, 10, 300)
    mirror_discharged = 0
    for _, index, _, _ in mirror["unresolved"]:
        if index not in (0, 1, cells - 1):
            raise AssertionError(f"uncovered mirror strip x-cell {index}")
        mirror_discharged += 1
    if mirror_discharged != mirror["unresolved_count"]:
        raise AssertionError("mirror unresolved cells were not all discharged")
    corner = scan_near_one_corner(10, 300)
    if corner["unresolved_count"]:
        raise AssertionError("near-one geometric corner is unresolved")

    delta = arb_fraction(Fraction(1, 2**300))
    logarithm = arb(300) * arb(2).log()
    coefficient = (parameters.mean - arb(1) / 2) * (1 - delta)
    tail_margin = (
        2 * coefficient * logarithm - 3 * parameters.mean
        - parameters.beta * parameters.mean * delta * (2 * logarithm + 1)
    )
    if not tail_margin > 0:
        raise AssertionError("near-one analytic tail margin is not positive")
    edge_entropy = min(
        _OPS.entropy(arb_fraction(Fraction(1, 32))).lower(),
        _OPS.entropy(arb_fraction(Fraction(1, 1024))).lower(),
        key=float,
    )
    strip_tail_margin = (
        (parameters.mean - arb(1) / 2) * edge_entropy
        - 20 * (_OPS.entropy(delta) + delta)
    )
    if not strip_tail_margin > 0:
        raise AssertionError("near-one strip tail margin is not positive")
    corner_scale = arb_fraction(Fraction(1, 1024))
    wedge_coefficient = (
        (parameters.mean - arb(1) / 2) * (1 - corner_scale))
    wedge_tail_margin = (
        wedge_coefficient * logarithm - 3 * parameters.mean
        - parameters.beta * parameters.mean * delta * (2 * logarithm + 1)
    )
    if not wedge_tail_margin > 0:
        raise AssertionError("asymmetric near-one tail margin is not positive")
    corner_logarithm = arb(10) * arb(2).log()
    wedge_t_monotonicity = (
        wedge_coefficient * corner_logarithm - parameters.mean
        - parameters.beta * parameters.mean * corner_scale
    )
    wedge_r_monotonicity = (
        wedge_coefficient * (corner_logarithm - 1) - parameters.mean
        - 2 * parameters.beta * parameters.mean * delta * logarithm
    )
    if not wedge_t_monotonicity > 0 or not wedge_r_monotonicity > 0:
        raise AssertionError("asymmetric wedge monotonicity is not positive")

    return {
        "claim_status": "PROVED",
        "identity": (
            "F(P)=M E_{QxQ} D_M(x,y)/(xy), Q(dx)=xP(dx)/M"
        ),
        "bulk_cells_per_axis": cells,
        "bulk_boxes": bulk["boxes"],
        "bulk_machine_unresolved_raw": bulk["unresolved_count"],
        "bulk_analytic_discharged": bulk_discharged,
        "bulk_unresolved_after_composition": 0,
        "small_boundary_cells": cells // 32,
        "small_boundary_margin": str(small_margin),
        "optimizer": optimizer,
        "mirror_strip": {
            "first_power": 10,
            "last_power": 300,
            "checked": mirror["checked"],
            "machine_unresolved": mirror["unresolved_count"],
            "analytic_discharged": mirror_discharged,
            "unresolved_after_composition": 0,
        },
        "mirror_corner": {
            "checked": corner["checked"],
            "machine_unresolved": corner["unresolved_count"],
        },
        "tail_delta": "2^-300",
        "tail_margin": str(tail_margin),
        "wedge_tail_margin": str(wedge_tail_margin),
        "wedge_t_monotonicity": str(wedge_t_monotonicity),
        "wedge_r_monotonicity": str(wedge_r_monotonicity),
        "strip_tail_margin": str(strip_tail_margin),
        "analytic_lemmas": {
            "small_support": (
                "For 0<=x<=1/32 and arbitrary 0<=y<=1, "
                "y[1+(1-x)(1-y)]<=1 gives pi(x,y)<=x<1/2 and "
                "protocol_gain>=0; "
                "mu decreasing gives T<=mu(y); concavity makes h(x)/x "
                "decreasing. Hence D_m>=xy[(m-1/2)h(1/32)/(1/32)-1/2], "
                "the certified small_boundary_margin."
            ),
            "optimizer": (
                "The defining equality and stationarity equations (87)-(90) "
                "give D_m(x*,x*)=grad D_m(x*,x*)=0 exactly. The certified "
                "positive-definite Hessian on [17/25,7/10]^2 makes x* its "
                "minimum throughout the local box."
            ),
            "near_one_strip_tail": (
                "At y=1, D_m(x,1)=(m-1/2)h(x). On "
                "x in [1/32,1-1/1024], entropy continuity and T<=mu(y) "
                "bound the total change for 1-y<=delta by "
                "20[h(delta)+delta], below strip_tail_margin. The omitted "
                "small-x and near-(1,1) corners are covered separately."
            ),
            "asymmetric_near_one_tail": (
                "If t=min(1-x,1-y)<=delta and "
                "s=max(1-x,1-y)<=1/1024, T<=t[3+log(s/t)] and "
                "protocol_gain>=-h(st). The certified wedge_t_monotonicity "
                "shows the lower bound is minimized at t=delta. Writing "
                "s=delta*r, wedge_r_monotonicity includes both -m/r and the "
                "full -beta*m*delta*(2L-log r) derivative term and proves "
                "increase for r>=1; the r=1 value is wedge_tail_margin>0. "
                "This covers the asymmetric infinite wedge beyond the finite "
                "geometric shells."
            ),
            "near_one_tail": (
                "For s=1-x<=t=1-y<=delta, the exact positive series for "
                "T gives T<=s[3+log(t/s)]; entropy continuity gives "
                "protocol_gain>=-h(st), with h(z)<=z(log(1/z)+1). "
                "After h(s)>=s log(1/s), the lower bracket is increasing in "
                "u=log(t/s) and its u=0 value is tail_margin>0."
            ),
            "mean_monotonicity": (
                "Writing D_M=M E-A/2 with A=yh(x)+xh(y)>=0, D_m>=0 implies "
                "E=(D_m+A/2)/m>=0; therefore D_M>=D_m for every M>=m."
            ),
        },
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=int, default=256)
    parser.add_argument("--certify", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT_DEFAULT)
    args = parser.parse_args(argv)
    if not args.certify:
        print(scan(args.cells))
        return 0
    report = certify_qone_kernel()
    report["report_sha256"] = _canonical_digest(report)
    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    print("LIU H2 SIZE-BIASED Q-ONE KERNEL")
    print("PROVED bulk=%d tail_margin=%s"
          % (report["bulk_boxes"], report["tail_margin"]))
    print("report_sha256 %s" % report["report_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
