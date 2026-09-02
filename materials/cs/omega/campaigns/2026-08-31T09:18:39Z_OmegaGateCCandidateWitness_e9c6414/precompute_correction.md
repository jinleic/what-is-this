# Pre-compute protocol correction (2026-08-31)

No base or moved objective evaluation has run.

The initially frozen `gate_c_candidate_witness.py.asrun` is **SUPERSEDED —
STATIC GUARDS ONLY**. It saved the first strict-negative rung but continued
through later rungs, contrary to the pre-statement decision rule (“freeze the
first such rung and all attempted earlier rungs”). No result was produced by
that source.

The corrected as-run source is
`gate_c_candidate_witness_corrected.py.asrun`. After saving a strict-negative
rung it breaks immediately, and `decision.evaluated_divisors` is the actual
attempted prefix rather than the full registered ladder. Direction, ladder,
objective arithmetic, evidence gates, command, and resource cap are unchanged.
