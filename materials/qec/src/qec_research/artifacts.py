"""Artifact-write routing.

Single source of truth for where a verification run's output lands.  The
canonical artifact may only be produced by a run that both PASSES every check
and covers the FULL catalogue; anything else lands in a clearly non-canonical
location and never overwrites canonical evidence.
"""
from __future__ import annotations


def canonical_route(clean: bool, full_coverage: bool) -> str:
    """'canonical' | 'partial_runs' (clean but incomplete) | 'quarantine' (failing)."""
    if clean and full_coverage:
        return "canonical"
    return "partial_runs" if clean else "quarantine"
