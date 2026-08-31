"""Depth-three [[8,3,2]] encoder transcribed from arXiv:2605.21867 Fig. 4(a).

The rendered figure contains ten CNOT symbols despite the surrounding text
calling the construction an "eight CNOT" encoder.  The circuit below follows
the symbols (and their control/target dots), with four |0> inputs and four
|+> inputs.  Its ideal output is checked against the eight commuting state
stabilizers in ``smoke.py``.
"""

from __future__ import annotations

try:  # Support both ``python src/file.py`` and package imports.
    from .layout import ENCODER_LAYERS
except ImportError:  # pragma: no cover - exercised by the CLI invocation.
    from layout import ENCODER_LAYERS

ZERO_INPUTS = (0, 2, 5, 7)
PLUS_INPUTS = (1, 3, 4, 6)


def apply_encoder(builder, offset: int = 0, *, inverse: bool = False) -> None:
    """Append the figure's encoder (or its inverse) to a noise builder.

    All encoder gates are CNOTs, hence inverse execution is the same gate list
    in reverse depth order.  ``builder`` must provide ``cx_layer``.
    """

    layers = ENCODER_LAYERS[::-1] if inverse else ENCODER_LAYERS
    for layer in layers:
        builder.cx_layer([(offset + control, offset + target) for control, target in layer])


def input_basis(builder, offset: int = 0) -> None:
    """Prepare the four |0> and four |+> inputs used by Fig. 4(a)."""

    builder.reset_many([offset + q for q in ZERO_INPUTS], basis="Z")
    builder.reset_many([offset + q for q in PLUS_INPUTS], basis="X")
