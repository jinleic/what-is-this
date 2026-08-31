"""task0_witness.py — PRIORITY ZERO: witness test on Sun's 13-gate left circuit.

Expands Sun's printed 13-addition left SLP (arXiv:2604.27645 §4) and tests,
over Z, exactly which 9->23 factor maps it computes — coefficient by
coefficient against: (1) Sun's OWN Appendix-6 blocks (U/V/W), (2) the
55-paper's U, V, W blocks (my tensor_data.py). The raw blocks are re-cut
from the fetched HTML artifact itself (no hand transcription anywhere).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from tensor_data import U_BLOCK_PRINTED, V_BLOCK_PRINTED, W_BLOCK_PRINTED
from tensor_data_sun56 import SUN_LEFT_SLP, SUN_RIGHT_SLP, SUN_PRODUCTS


def expand(slp, base):
    env = {}
    for i in range(9):
        env[f"{base}{i}"] = {f"{base}{i}": 1}
    for name, terms in slp:
        acc: dict[str, int] = {}
        for sign, atom in terms:
            for k, v in env[atom].items():
                acc[k] = acc.get(k, 0) + sign * v
        env[name] = acc
    return env


def left_factor_vectors():
    """The 23 exact Z^9 vectors Sun's left circuit computes (with the printed
    signs of the product list)."""
    L = expand(SUN_LEFT_SLP, "A")
    vecs = []
    for (ls, la), _ in SUN_PRODUCTS:
        t = L[la] if la.startswith("u") else {la: 1}
        if ls == -1:
            t = {k: -v for k, v in t.items()}
        vecs.append(tuple(t.get(f"A{i}", 0) for i in range(9)))
    return vecs


def block_cols(block):
    """block: 9 rows (entries) x 23 cols (products) -> 23 column vectors."""
    return [tuple(block[i][r] for i in range(9)) for r in range(23)]


def canon(v):
    neg = tuple(-x for x in v)
    return v if v <= neg else neg


def multiset_pairs(A, B):
    """Test whether A is B under a product permutation, with per-column sign
    freedom reported separately (exact and sign-exact variants)."""
    from collections import Counter
    exact = Counter(A) == Counter(B)
    signfree = Counter(canon(a) for a in A) == Counter(canon(b) for b in B)
    return bool(exact), bool(signfree)


def main():
    sun_blocks = json.load(open(HERE.parent / 'scratch' / 'sun56_blocks.json'))
    SU, SV, SW = sun_blocks["U"], sun_blocks["V"], sun_blocks["W"]

    out = {}
    sun_vecs = left_factor_vectors()
    out["sun_left_n_vectors"] = len(sun_vecs)

    for name, blk in (("sunU", SU), ("sunV", SV), ("sunW", SW),
                      ("paperU", U_BLOCK_PRINTED), ("paperV", V_BLOCK_PRINTED),
                      ("paperW", W_BLOCK_PRINTED)):
        cols = block_cols(blk)
        ex, sf = multiset_pairs(sun_vecs, cols)
        out[f"sun_left_equals_{name}_cols"] = {"exact_order_agnostic": ex,
                                               "up_to_sign_per_column": sf}

    # adjacency-free statement: print Sun's vectors vs paper V columns
    def fmt(v):
        return "".join(str(x) for x in v)
    out["sun_left_vectors"] = [fmt(v) for v in sun_vecs]
    out["paperV_columns"] = [fmt(v) for v in block_cols(V_BLOCK_PRINTED)]
    out["paperU_columns"] = [fmt(v) for v in block_cols(U_BLOCK_PRINTED)]
    out["paperW_columns"] = [fmt(v) for v in block_cols(W_BLOCK_PRINTED)]

    print(json.dumps(out, indent=2))
    (HERE.parent / 'scratch' / 'witness_result.json').write_text(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    main()
