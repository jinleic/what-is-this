"""Seed ledger: deterministic per-run seeds, recorded commitments.

Policy (physics/README.md resource rules): fixed seed per run, recorded in
the campaign inventory; repeated-seed policy defined before first run.
We derive per-(point, decoder) seeds from a base seed with a stable hash
so the whole campaign reproduces from (base_seed, campaign config).
"""
from __future__ import annotations

import hashlib
import json


def derive_seed(base_seed: int, *labels: object) -> int:
    """Stable 63-bit seed from a base seed and label strings."""
    payload = json.dumps([base_seed, *map(str, labels)]).encode()
    digest = hashlib.sha256(payload).digest()
    return int.from_bytes(digest[:8], "big") % (2**63)


def sampling_seed(base_seed: int, p: float, basis: str, arm: str) -> int:
    return derive_seed(base_seed, "sampler", p, basis, arm)


def decode_seed(base_seed: int, p: float, basis: str, arm: str, shot: int) -> int:
    return derive_seed(base_seed, "decoder", p, basis, arm, shot)
