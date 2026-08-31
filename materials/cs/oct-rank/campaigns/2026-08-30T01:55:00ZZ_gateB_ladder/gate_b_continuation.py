"""gate_b_continuation.py — BORDER-RANK PROBE (Main-directed reframing).

Path: start from the exact-verified rank-25 certificate x0, scale term t's
A-column by s in {0.5,0.2,0.1,0.05,0.02,0.01,0}, LM-refine AT EACH s over
ALL 25 columns (25-column configuration — NOT rank 24), then at s=0 remove
the column entirely and run the true 24-column LM.

Recorded JOINTLY at each s: rel residual, per-factor column norms, max norm,
norm of the scaled term (pre-LM), bookkeeping labels ("25col_scaled" vs
"24col_removed").  No rank claim is made from 25-column rows.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gate_b_search3 import TF, N1, N2, N3, residual, factor_stats  # noqa
from gate_b_deflate import load_cert_x0, split, lm_refine  # noqa: E402

R = 25
TERM = int(os.environ.get("TERM", "0"))
OUT = os.environ.get("OUTFILE", "gate_b_continuation_term0.jsonl")


def main():
    x0 = load_cert_x0()
    rows = []
    for s in [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.0]:
        aa, bb, cc = split(x0)
        aa[:, TERM] *= s
        term_norm = float(np.linalg.norm(bb[:, TERM]) * np.linalg.norm(cc[:, TERM]) * np.linalg.norm(aa[:, TERM]))
        xs = np.concatenate([aa.ravel(), bb.ravel(), cc.ravel()])
        rel_start = float(np.linalg.norm(residual(xs, R)) / TF)
        x, rel, nfev, status, mx, per, ratio, prod = lm_refine(xs, R, 500)
        row = {"phase": "25col_scaled", "s": s, "term": TERM,
               "term_a_col_norm": float(np.linalg.norm(aa[:, TERM])),
               "term_norm_preLM": term_norm, "rel_start": rel_start,
               "rel_after_LM": rel, "nfev": nfev, "status": status,
               "max_factor_norm": mx, "factor_norms": per,
               "max_column_product": prod, "ratio": ratio}
        rows.append(row)
        print(json.dumps(row), flush=True)
        if rel > 1e-9 and s < 0.05:
            break
    # true 24-column endpoint at s=0
    aa, bb, cc = split(x0)
    aa = np.delete(aa, TERM, axis=1)
    bb = np.delete(bb, TERM, axis=1)
    cc = np.delete(cc, TERM, axis=1)
    x24 = np.concatenate([aa.ravel(), bb.ravel(), cc.ravel()])
    rel_start = float(np.linalg.norm(residual(x24, 24)) / TF)
    x, rel, nfev, status, mx, per, ratio, prod = lm_refine(x24, 24, 600)
    row = {"phase": "24col_removed", "s": None, "term": TERM,
           "rel_start": rel_start, "rel_after_LM": rel, "nfev": nfev,
           "status": status, "max_factor_norm": mx, "factor_norms": per,
           "max_column_product": prod, "ratio": ratio}
    rows.append(row)
    print(json.dumps(row), flush=True)
    with open(OUT, "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
    print("CONTINUATION_DONE " + OUT)


if __name__ == "__main__":
    main()
