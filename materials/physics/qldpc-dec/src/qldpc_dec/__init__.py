"""cross-paper decoder reproduction harness for bivariate bicycle codes.

Modules:
    code         BB code construction + stim circuit builder (documented fallback)
    circuits     circuit/DEM loading conventions for the harness
    harness      BBCodeHarness: shared decode interface (contract with fss-bb)
    bp_osd       BP+OSD baseline (ldpc package), gate A
    beam_search  Beam-search decoder reimplementation (arXiv:2512.07057), gate B
    nms          Normalized-min-sum decoder + 24-ensemble, gate C
    bootstrap    Binomial bootstrap CIs
    runner       Campaign runner: immutable snapshots under campaigns/<UTC>_<uuid>_<hash>/
    seeds        Seed ledger helpers

Run conventions (fixed in pre_statement.md, 2026-08-29):
    - circuits: IonQ stim circuits from github.com/ionq-publications/beamsearchdecoder
    - LER definition: per-basis (errors/shots)/12, total = X + Z
"""
from .harness import BBCodeHarness

__all__ = ["BBCodeHarness"]
