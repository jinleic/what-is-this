#!/usr/bin/env python3
"""Pre-registered startup anchors and counterfactual controls."""
from __future__ import annotations

import json
import os

from flint import arb, fmpq

import core2
import gate_B2_inherited as inherited

PREC = 256


def strict_margin_lower(d_ball: arb, upper_ball: arb) -> arb:
    """Certified lower endpoint of d-upper; positive iff separated."""
    with core2.Prec(PREC):
        return d_ball.lower() - upper_ball.upper()


def count_root_boxes(n_side: int) -> int:
    """Exact rational grid; count closed boxes meeting the closed unit disk."""
    count = 0
    with core2.Prec(PREC):
        w = arb(fmpq(2, n_side))
        for i in range(n_side):
            for j in range(n_side):
                cx = arb(-1) + arb(fmpq(2 * i + 1, n_side))
                cy = arb(-1) + arb(fmpq(2 * j + 1, n_side))
                box = inherited.box_from_center(cx, cy, w)
                if inherited.box_intersects_disk(*box):
                    count += 1
    return count


def main() -> None:
    with core2.Prec(PREC):
        # Basis convention and the paper's reproducing-kernel anchor.
        a0 = (arb(2) / 3).sqrt()
        a2 = -(arb(1) / 3).sqrt()
        norm_residual = a0 * a0 + a2 * a2 - 1
        e0 = a0 + a2 * core2.psi(2, arb(0))
        p0_residual = e0 * e0 - arb(3) / 2
        assert norm_residual.contains(0), norm_residual
        assert p0_residual.contains(0), p0_residual

        # Rule-17b count anchor, independently fixed in pre_statement.md.
        root_count = count_root_boxes(12)
        assert root_count == 132, root_count

        # Counterfactual panel weights on [0,1].  phi(1) is a lower, not
        # upper, rectangle because phi decreases on the positive half-line.
        exact_weight = inherited.cdf_weight(arb(0), arb(1))
        right_rect = core2.phi_density(arb(1))
        left_rect = core2.phi_density(arb(0))
        assert right_rect.upper() < exact_weight.lower(), (right_rect, exact_weight)
        assert exact_weight.upper() < left_rect.lower(), (exact_weight, left_rect)

        # Counterfactual margin signs, including an overlapping ball.
        d1 = arb(1)
        rejected_above = strict_margin_lower(d1, arb('1.0001'))
        rejected_overlap = strict_margin_lower(d1, arb(1, arb('0.001')))
        accepted_below = strict_margin_lower(d1, arb('0.9999'))
        assert rejected_above < 0
        assert rejected_overlap < 0
        assert accepted_below > 0

        # Historical single-cell anchor, using certified endpoints rather
        # than the inherited driver's permissive margin extraction.
        eps = arb('0.0001')
        box = (arb('0.8') - eps, arb('0.8') + eps,
               arb('-0.6') - eps, arb('-0.6') + eps)
        old_upper = inherited.cell_envelope(arb('1.30'), box, 512)
        d_anchor = core2.d_of_c(arb('1.45'))
        old_margin_lower = strict_margin_lower(d_anchor, old_upper)
        assert old_margin_lower > arb('0.192'), old_margin_lower
        assert old_margin_lower < arb('0.194'), old_margin_lower

        result = {
            "precision_bits": PREC,
            "basis": {
                "a0": a0.str(40),
                "a2": a2.str(40),
                "norm_residual": norm_residual.str(12),
                "p0_at_zero_squared_residual": p0_residual.str(12),
            },
            "root_box_count_n12": root_count,
            "weight_counterfactual_0_1": {
                "right_phi_b_lower_rectangle": right_rect.str(40),
                "exact_cdf_mass": exact_weight.str(40),
                "left_phi_a_upper_rectangle": left_rect.str(40),
            },
            "margin_counterfactuals": {
                "above_rejected_lower_margin": rejected_above.str(30),
                "overlap_rejected_lower_margin": rejected_overlap.str(30),
                "below_accepted_lower_margin": accepted_below.str(30),
            },
            "inherited_anchor": {
                "c_pair": ["1.30", "1.45"],
                "center": ["0.8", "-0.6"],
                "radius_each_axis": "0.0001",
                "n_panels": 512,
                "upper": old_upper.str(40),
                "upper_radius": old_upper.rad().str(20),
                "d": d_anchor.str(40),
                "d_radius": d_anchor.rad().str(20),
                "certified_margin_lower": old_margin_lower.str(40),
            },
            "status": "PASS",
        }

    out = os.path.join(os.path.dirname(__file__), "..", "logs", "startup_controls.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
        f.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
