"""Scaled QMAP zoned architecture with exact QMAP coordinates and timing.

Coordinates follow QMAP ``Architecture::exactSLMLocation``:
  x = site_separation.first * col + location.x
  y = site_separation.second * row + location.y
Site separations and operation durations equal QMAP's sqare architecture;
zone extents scale down to the instance's layer capacity.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass

from .config import (
    ENTANGLEMENT_SITE_SEPARATION,
    QMAP_SQUARE_ARCH,
    STORAGE_SITE_SEPARATION,
)


@dataclass(frozen=True)
class Site:
    """A trap site: zone ('s'|'e'), slm index, row, column."""

    zone: str
    slm: int
    r: int
    c: int

    def xy(self, arch: "Architecture") -> tuple[float, float]:
        return arch.site_xy(self)


class Architecture:
    """One storage SLM and one paired entanglement zone (two SLMs).

    The entanglement zone sits above storage with the same x/y offsets as
    QMAP's ``square_architecture.json`` (offset [1, 309]); apart from zone
    sizes every parameter equals the shipped machine.
    """

    def __init__(
        self,
        n_storage_rows: int,
        n_storage_cols: int,
        n_ent_rows: int,
        n_ent_cols: int,
    ):
        assert n_storage_rows >= 2 and n_storage_cols >= 2
        assert n_ent_rows >= 1 and n_ent_cols >= 1
        self.n_storage_rows = n_storage_rows
        self.n_storage_cols = n_storage_cols
        self.n_ent_rows = n_ent_rows
        self.n_ent_cols = n_ent_cols
        self.ent_offset = (1.0, 309.0)

    # ---------------- geometry ----------------
    def site_xy(self, site: Site) -> tuple[float, float]:
        if site.zone == "s":
            return (
                site.c * STORAGE_SITE_SEPARATION[0],
                site.r * STORAGE_SITE_SEPARATION[1],
            )
        # The two entanglement SLMs interleave columns so that each pair of
        # sites row r/col c is (12,10) apart, mirroring QMAP's locations
        # x=1 and x=3 with column pitch 12.
        slm_xoff = 1.0 + 2.0 * site.slm
        return (
            self.ent_offset[0] + slm_xoff + site.c * ENTANGLEMENT_SITE_SEPARATION[0],
            self.ent_offset[1] + site.r * ENTANGLEMENT_SITE_SEPARATION[1],
        )

    def dist(self, a: Site, b: Site) -> float:
        ax, ay = self.site_xy(a)
        bx, by = self.site_xy(b)
        return math.hypot(ax - bx, ay - by)

    def n_sites(self, zone: str) -> int:
        if zone == "s":
            return self.n_storage_rows * self.n_storage_cols
        return 2 * self.n_ent_rows * self.n_ent_cols

    def storage_sites(self) -> list[Site]:
        return [
            Site("s", 0, r, c)
            for r in range(self.n_storage_rows)
            for c in range(self.n_storage_cols)
        ]

    def entanglement_sites(self) -> list[Site]:
        return [
            Site("e", slm, r, c)
            for slm in (0, 1)
            for r in range(self.n_ent_rows)
            for c in range(self.n_ent_cols)
        ]

    def initial_placement(self, n_qubits: int) -> list[Site]:
        """QMAP VertexMatchingPlacer.makeInitialPlacement: row-major fill of
        the first storage SLM starting at row 0, column 0."""
        assert n_qubits <= self.n_sites("s")
        return [
            Site("s", 0, i // self.n_storage_cols, i % self.n_storage_cols)
            for i in range(n_qubits)
        ]

    # ---------------- serialization ----------------
    def to_quimap_json(self) -> str:
        """mqt.qmap-loadable JSON with scaled dimensions; QMAP durations."""
        arch = dict(QMAP_SQUARE_ARCH)
        arch["storage_zones"] = [
            {
                "zone_id": 0,
                "slms": [
                    {
                        "id": 0,
                        "site_separation": [4, 4],
                        "r": self.n_storage_rows,
                        "c": self.n_storage_cols,
                        "location": [0, 0],
                    }
                ],
                "offset": [0, 0],
                "dimension": [
                    self.n_storage_cols * 4,
                    self.n_storage_rows * 4,
                ],
            }
        ]
        arch["entanglement_zones"] = [
            {
                "zone_id": 0,
                "slms": [
                    {
                        "id": 1,
                        "site_separation": [12, 10],
                        "r": self.n_ent_rows,
                        "c": self.n_ent_cols,
                        "location": [1, 309],
                    },
                    {
                        "id": 2,
                        "site_separation": [12, 10],
                        "r": self.n_ent_rows,
                        "c": self.n_ent_cols,
                        "location": [3, 309],
                    },
                ],
                "offset": [1, 309],
                "dimension": [
                    3 + self.n_ent_cols * 12,
                    self.n_ent_rows * 10,
                ],
            }
        ]
        arch["rydberg_range"] = [[[1, 309], [3 + self.n_ent_cols * 12, 309 + self.n_ent_rows * 10]]]
        return json.dumps(arch)

    def manifest(self) -> dict:
        return {
            "storage_rows": self.n_storage_rows,
            "storage_cols": self.n_storage_cols,
            "ent_rows": self.n_ent_rows,
            "ent_cols": self.n_ent_cols,
            "site_separation_storage": STORAGE_SITE_SEPARATION,
            "site_separation_ent": ENTANGLEMENT_SITE_SEPARATION,
            "rydberg_range": [[1, 309], [3 + self.n_ent_cols * 12, 309 + self.n_ent_rows * 10]],
            "ent_offset": self.ent_offset,
        }


def architecture_for(n_qubits: int, max_layer: int) -> Architecture:
    """Scaled architecture: storage >= 2n sites, entanglement capacity >=
    widest ASAP layer (2 sites per CZ: one SLM per operand)."""
    cols = max(n_qubits, 2)
    rows = max(2, math.ceil((2 * n_qubits) / cols))
    ent_cols = max(math.ceil(max_layer / 2), 2)
    ent_rows = 1
    return Architecture(rows, cols, ent_rows, ent_cols)
