"""Rank-13 starts from the certified two-block shared-factor tau upper side.

Seven starts delete one matched tau term (block-swap symmetry removes the
other seven).  Seven starts replace a matched pair by one cross-block term,
retaining both diagonal contributions and introducing only off-diagonal error.
Each start is refined through pivoted square CP systems.  Floating results are
discovery evidence only.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from route_af_refine import refine
from route_af_search import T, TF, RANK, balance, factor_stats, residual_for

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CAMPAIGN = ROOT / "campaigns" / "2026-08-31T08:02:18Z_routeAF"
CERT = ROOT / "scratch" / "upstream_ref" / "certs" / "tower_cert" / "tau_r7"
TRIGGER = 1e-11


def load(name: str) -> np.ndarray:
    with (CERT / name).open(newline="") as handle:
        return np.asarray([[float(v) for v in row] for row in csv.reader(handle)], dtype=float)


def rank14() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    at, bt, ct = load("A.csv"), load("B.csv"), load("C.csv")
    assert at.shape == (3, 7) and bt.shape == ct.shape == (4, 7)
    a = np.concatenate([at, at], axis=1)
    b = np.zeros((8, 14)); c = np.zeros((8, 14))
    b[:4, :7] = bt; c[:4, :7] = ct
    b[4:, 7:] = bt; c[4:, 7:] = ct
    rel = float(np.linalg.norm(np.einsum("pr,br,cr->pbc", a, b, c) - T) / TF)
    assert rel < 1e-12, rel
    return a, b, c


def join(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    return balance(np.concatenate([a.ravel(), b.ravel(), c.ravel()]))


def drop_start(term: int) -> np.ndarray:
    a, b, c = rank14()
    keep = [s for s in range(14) if s != term]
    return join(a[:, keep], b[:, keep], c[:, keep])


def merge_start(term: int) -> np.ndarray:
    a, b, c = rank14()
    mate = term + 7
    keep = [s for s in range(14) if s not in (term, mate)]
    am = a[:, term]
    b1, b2 = b[:4, term], b[4:, mate]
    c1, c2 = c[:4, term], c[4:, mate]
    cross12 = np.linalg.norm(b1) ** 2 * np.linalg.norm(c2) ** 2
    cross21 = np.linalg.norm(b2) ** 2 * np.linalg.norm(c1) ** 2
    scale = (cross21 / cross12) ** (1.0 / 8.0)
    bm = np.concatenate([scale * b1, b2 / scale])
    cm = np.concatenate([c1 / scale, scale * c2])
    aa = np.column_stack([a[:, keep], am])
    bb = np.column_stack([b[:, keep], bm])
    cc = np.column_stack([c[:, keep], cm])
    assert aa.shape == (3, RANK)
    return join(aa, bb, cc)


def main() -> None:
    a14, b14, c14 = rank14()
    anchor_rel = float(np.linalg.norm(np.einsum("pr,br,cr->pbc", a14, b14, c14) - T) / TF)
    print("AF_STRUCTURED_ANCHOR " + json.dumps({
        "tau_terms_per_block": 7, "rank14_stored_decimal_rel": anchor_rel,
        "semantics": "stored decimals approximate; tau Krawczyk cert supplies exact nearby factors",
    }))
    starts = [("drop", t, drop_start(t)) for t in range(7)]
    starts += [("merge", t, merge_start(t)) for t in range(7)]
    best = None
    for kind, term, x0 in starts:
        rel0 = float(np.linalg.norm(residual_for(T, x0, RANK)) / TF)
        print("AF_STRUCTURED_START " + json.dumps({
            "kind": kind, "term": term, "rel": rel0, **factor_stats(x0)
        }), flush=True)
        x, selected, trace = refine(T, x0, rounds=3, max_nfev=1500)
        rel = float(np.linalg.norm(residual_for(T, x, RANK)) / TF)
        row = {"kind": kind, "term": term, "rel_before": rel0, "rel": rel,
               "selected_count": int(len(selected)), **factor_stats(x)}
        print("AF_STRUCTURED_RESULT " + json.dumps(row), flush=True)
        if best is None or rel < best[0]:
            best = (rel, x, selected, row, trace)
        if rel < TRIGGER:
            break
    assert best is not None
    rel, x, selected, row, trace = best
    out = CAMPAIGN / "rank13_structured_candidate.npz"
    np.savez(out, x=x, selected=selected, target=T, rel=rel,
             row=json.dumps(row), trace=json.dumps(trace),
             convention="T[p,b,c]=sum_s A[p,s]B[b,s]C[c,s]; globally conjugated Route-F target is the shared-first-factor block duplication tau boxtimes (1 x I2)")
    print("AF_STRUCTURED_BEST " + json.dumps({"path": str(out), **row}))
    if rel >= TRIGGER:
        raise SystemExit(2)


if __name__ == "__main__":
    main()

