"""Campaign W runner — weight-exactly-4 census at PN4 (421,(4,5,7),t=(1,1,1)).

Instrument: vectorized (numpy, exact residues mod q) left-kernel test of the
V-basis' off-S columns per support; zero float anywhere. Semantics match the
frozen gate-B/P2 instrument exactly:
  - nonzero-support count: supports S with dim(V∩F^S) > 0;
  - candidate class: for each nonzero support, map the kernel basis to
    V-vectors, rref in F_q^N (canonical for the subspace), normalize each
    row (first nonzero coord = 1), dedupe across supports.
  Closure identity (hand-verified at PN6): 42 wt-2 candidates -> supports
  42 + 42*82 = 3486. Anchors assert BOTH quantity types.

Modes:
  --anchors-only    validation ladder (13/13 gate-B candidate counts, P2
                    support+candidate anchors, cross-check vs unvectorized
                    frozen path, PN4 w<=3 re-anchor + 35-line w=4 sanity).
  --kill-demo       census the first KILL_AT supports, checkpoint, hard-exit.
  --run             census with checkpoint/resume (immutable-ledger).
  --validate-resume print validated ledger state as JSON.

Ledger: census_ledger.jsonl canonical-JSON lines, each sha256-stamped
(delcap validate_resume pattern). A checksum failure or index-order gap
hard-aborts; nothing is ever truncated or repaired.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import os
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RS = HERE.parents[1]
sys.path.insert(0, str(RS))
from src.campaign import _normalize  # noqa: E402  (frozen normalization)
from src.delta import DeltaEngine  # noqa: E402
from src.field import rref  # noqa: E402
from src.rs import Inst  # noqa: E402

Q = 421
S_PN4 = (4, 5, 7)
T_PN4 = (1, 1, 1)
N_PN4 = 140
W3_TOTAL = 457_450          # C(140,1)+C(140,2)+C(140,3)  (pre-registered)
W4_TOTAL = 15_329_615       # C(140,4)                     (pre-registered)
CKPT_EVERY = 400_000
LEDGER = HERE / "census_ledger.jsonl"


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def append_rec(kind: str, data: dict) -> None:
    body = {"kind": kind, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                               time.gmtime()), **data}
    canon = json.dumps(body, sort_keys=True, separators=(",", ":"))
    body["_sha256"] = sha(canon.encode())
    line = json.dumps(body, sort_keys=True, separators=(",", ":")) + "\n"
    existing = LEDGER.read_bytes() if LEDGER.exists() else b""
    tmp = LEDGER.with_name(LEDGER.name + ".tmp")
    with open(tmp, "wb") as f:
        f.write(existing + line.encode())
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, LEDGER)
    dfd = os.open(HERE, os.O_RDONLY)
    try:
        os.fsync(dfd)
    finally:
        os.close(dfd)


def load_valid_ledger() -> list[dict]:
    if not LEDGER.exists():
        return []
    out = []
    for ln in LEDGER.read_text().splitlines():
        if not ln.strip():
            continue
        d = json.loads(ln)
        s = d.pop("_sha256")
        got = sha(json.dumps(d, sort_keys=True, separators=(",", ":")).encode())
        if got != s:
            raise SystemExit(f"HARD ABORT: ledger checksum-invalid; keep "
                             f"{LEDGER} and inspect; never truncate/repair")
        out.append(d)
    last = -1
    for d in out:
        ni = d.get("next_index")
        if ni is None:
            continue
        if ni < last:
            raise SystemExit(f"HARD ABORT: ledger index-order gap at {ni}")
        last = ni
    return out


# ---------------- instrument ----------------

def build_instruments(q: int, s: tuple, t: tuple):
    I = Inst(q, s, t)
    eng = DeltaEngine(I)
    rows = []
    for i in range(3):
        rows.extend(eng.lift_rref(i))
    vb, _ = rref(rows, q)
    return I, eng, np.array(vb, dtype=np.int64)


def rank_mod_np(mat_cols: np.ndarray, q: int, max_rank: int) -> int:
    """Exact rank mod q; mat_cols: (r, c) with r <= max_rank (cols = off-S)."""
    m = mat_cols.T.copy() % q          # constraint rows
    rr, piv, rows, cols = 0, 0, *m.shape
    while rr < rows and piv < cols:
        nz = np.nonzero(m[rr:, piv])[0]
        if nz.size == 0:
            piv += 1
            continue
        i = rr + nz[0]
        if i != rr:
            m[[rr, i]] = m[[i, rr]]
        inv = pow(int(m[rr, piv]), -1, q)
        m[rr] = (m[rr] * inv) % q
        colv = m[:, piv].copy()
        colv[rr] = 0
        m = (m - np.outer(colv, m[rr])) % q
        rr += 1
        piv += 1
    return rr


def kernel_V_basis(sub: np.ndarray, B: np.ndarray, q: int):
    """Reduced (rref, F_q^N) basis of V ∩ F^S for off-columns `sub` = B[:,off]:
    kernel basis in coeff space -> V-vectors -> rref -> normalized rows."""
    rows, r = sub.shape
    m = sub.T.copy() % q              # (r, N-|S|) coefficient constraints? no:
    # x @ sub = 0 with x in F_q^r: constraints rows = sub.T  (cols r)
    m = sub.T.copy() % q              # (N-|S|, r)
    rr, piv_cols = 0, []
    n_rows, n_cols = m.shape
    for cc in range(n_cols):
        nz = np.nonzero(m[rr:, cc])[0]
        if nz.size == 0:
            continue
        i = rr + nz[0]
        if i != rr:
            m[[rr, i]] = m[[i, rr]]
        inv = pow(int(m[rr, cc]), -1, q)
        m[rr] = (m[rr] * inv) % q
        colv = m[:, cc].copy()
        colv[rr] = 0
        m = (m - np.outer(colv, m[rr])) % q
        piv_cols.append(cc)
        rr += 1
        if rr == n_rows:
            break
    free = [c for c in range(n_cols) if c not in piv_cols]
    vecs = []
    Bn = np.asarray(B, dtype=np.int64) % q     # (r, N): basis rows
    for f in free:
        x = np.zeros(n_cols, dtype=np.int64)   # coeff vector over the r rows
        x[f] = 1
        for i, pc in enumerate(piv_cols):
            x[pc] = (-m[i, f]) % q
        v = (x @ Bn) % q                        # (N,)
        vecs.append(v)
    if not vecs:
        return []
    rr2, _ = rref([list(map(int, v)) for v in vecs], q)
    return [_normalize(list(map(int, v)), q) for v in rr2]


def census_w3(I: Inst, B: np.ndarray, w: int, collect_cands: bool = False):
    """Full w<=w window: nonzero-support count (by weight), total, and
    optionally the distinct normalized candidate count."""
    n = I.N
    r = B.shape[0]
    by_w = {}
    cands = set()
    seen = 0
    for k in range(1, w + 1):
        cnt_k = 0
        for comb in itertools.combinations(range(n), k):
            mask = np.zeros(n, dtype=bool)
            mask[list(comb)] = True
            off = np.nonzero(~mask)[0]
            sub = B[:, off]
            rank = rank_mod_np(sub, I.q, r)
            if rank < r:
                cnt_k += 1
                if collect_cands:
                    for v in kernel_V_basis(sub, B, I.q):
                        cands.add(v)
            seen += 1
        by_w[k] = cnt_k
    return {"by_w": by_w, "total_nonzero": sum(by_w.values()),
            "supports_seen": seen, "candidates": len(cands),
            "cand_set": cands if collect_cands else None}


# ---------------- validation ladder ----------------

def anchors() -> dict:
    t0 = time.perf_counter()
    gateb = json.loads((RS / "campaigns" /
                        "2026-08-30T11-43-59Z_50595CC9_fffe84b0_gateB" /
                        "payload.json").read_text())
    out = []
    for row in gateb["rows"]:
        q, s, t = row["q"], tuple(row["s"]), tuple(row["t"])
        I, eng, B = build_instruments(q, s, t)
        res = census_w3(I, B, 3, collect_cands=True)
        want_c = row["candidate_count_total"]
        assert res["candidates"] == want_c, \
            f"ANCHOR FAIL {row['id']}: candidates {res['candidates']} != {want_c}"
        out.append({"id": row["id"], "candidates": res["candidates"],
                    "by_w": res["by_w"]})
        print(f"  anchor {row['id']}: candidates {res['candidates']} = {want_c}"
              f"  by_w {res['by_w']}", flush=True)
    # P2 anchors (support counts + candidate counts)
    p2 = json.loads((RS / "scratch" / "p2_checkpoint.json").read_text())
    p2rows = {r["id"]: r for r in p2["rows"] if isinstance(r, dict)}
    for pid, q, s, want_sup, want_c in [
            ("PN6", 43, (2, 6, 7), 3486, 42),
            ("PP2", 43, (2, 3, 7), 917, 77)]:
        I, eng, B = build_instruments(q, s, (1, 1, 1))
        res = census_w3(I, B, 3, collect_cands=True)
        assert res["total_nonzero"] == want_sup, \
            f"ANCHOR FAIL {pid}: supports {res['total_nonzero']} != {want_sup}"
        assert res["candidates"] == want_c, \
            f"ANCHOR FAIL {pid}: candidates {res['candidates']} != {want_c}"
        out.append({"id": pid, "supports": want_sup, "candidates": want_c,
                    "by_w": res["by_w"]})
        print(f"  anchor {pid}: supports {want_sup}, candidates {want_c}, "
              f"by_w {res['by_w']}", flush=True)
    n = len(out)
    append_rec("anchors", {"n": n, "detail": out,
                           "wall_s": round(time.perf_counter() - t0, 1)})
    return {"n": n, "detail": out}


def crosscheck_unvectorized() -> dict:
    """Vectorized vs frozen unvectorized path, exact agreement on BOTH
    support count and candidate count (the assignment's milestone ii)."""
    out = []
    for (q, s, t) in [(13, (2, 2, 4), (1, 1, 1)), (5, (2, 2, 2), (1, 1, 1)),
                      (31, (2, 3, 5), (1, 1, 1))]:
        I, eng, B = build_instruments(q, s, t)
        vec = census_w3(I, B, 3, collect_cands=True)
        # unvectorized: the frozen instrument, verbatim
        uni_cands = set()
        uni_by_w = {}
        for k in (1, 2, 3):
            cnt = 0
            for comb in itertools.combinations(range(I.N), k):
                basis = eng.V_inter_support(set(comb))
                if basis:
                    cnt += 1
                    for v in basis:
                        uni_cands.add(_normalize(list(v), I.q))
            uni_by_w[k] = cnt
        assert sum(uni_by_w.values()) == vec["total_nonzero"], \
            f"CROSSCHECK supports FAIL {(q, s, t)}"
        assert uni_cands == vec["cand_set"], \
            f"CROSSCHECK candidates FAIL {(q, s, t)}"
        out.append({"inst": [q, list(s), list(t)],
                    "supports": vec["total_nonzero"],
                    "candidates": len(uni_cands)})
        print(f"  crosscheck {(q, s, t)}: supports {vec['total_nonzero']}, "
              f"candidates {len(uni_cands)} — EXACT agreement", flush=True)
    append_rec("crosscheck", {"detail": out})
    return {"detail": out}


def w3_reanchor_and_w4_sanity() -> dict:
    """At PN4 itself: (a) complete w<=3 census empty (457,450 supports, 0
    nonzero — re-derives M's C-PN4 under the W kernel); (b) all 35
    direction-0 lines nonempty at w=4 (theorem T part 2)."""
    I, eng, B = build_instruments(Q, S_PN4, T_PN4)
    assert B.shape[0] == 68 and I.N == N_PN4
    res = census_w3(I, B, 3)
    assert res["supports_seen"] == W3_TOTAL, \
        f"A3 FAIL: {res['supports_seen']} != {W3_TOTAL}"
    assert res["total_nonzero"] == 0, f"w3 re-anchor FAIL: {res['total_nonzero']}"
    lines_ok = 0
    for l in range(I.num_lines(0)):
        idxs = I.line_indices(0, l)
        off = np.array([c for c in range(N_PN4) if c not in set(idxs)])
        rank = rank_mod_np(B[:, off], Q, 68)
        assert rank < 68, f"w4 sanity FAIL: line {l} empty at w=4"
        lines_ok += 1
    append_rec("w3_w4_sanity", {"w3_supports": W3_TOTAL, "w3_nonzero": 0,
                                "w4_lines_checked": lines_ok})
    return {"w3_supports": W3_TOTAL, "w3_nonzero": 0, "w4_lines": lines_ok}


# ---------------- the census ----------------

def census_w4(mode: str) -> None:
    I, eng, B = build_instruments(Q, S_PN4, T_PN4)
    r = B.shape[0]
    assert r == 68
    ledger = load_valid_ledger()
    done_idx, hits, final_seen = 0, [], None
    for d in ledger:
        if d["kind"] == "final":
            final_seen = d
        if "next_index" in d and d.get("kind") in ("progress", "final"):
            if d["next_index"] > done_idx and d["kind"] == "progress":
                done_idx = d["next_index"]
        if d["kind"] == "hits":
            hits = d["new_hits"]
    if final_seen:
        print("census already FINAL;", final_seen)
        return
    if done_idx >= W4_TOTAL:
        return finalize(hits)
    t0 = time.perf_counter()
    idx = done_idx
    since = 0
    print(f"census resume/start at index {done_idx}; {len(hits)} hits so far",
          flush=True)
    Bn = B % Q
    for comb in itertools.islice(itertools.combinations(range(N_PN4), 4),
                                 done_idx, None):
        mask = np.zeros(N_PN4, dtype=bool)
        mask[list(comb)] = True
        off = np.nonzero(~mask)[0]
        sub = Bn[:, off]
        rank = rank_mod_np(sub, Q, r)
        if rank < r:
            dim = r - rank
            basis = kernel_V_basis(sub, Bn, Q)
            hits.append({"S": list(comb), "dim": dim,
                         "basis": [list(v) for v in basis]})
        idx += 1
        since += 1
        if idx % 100_000 == 0:
            rate = (idx - done_idx) / (time.perf_counter() - t0)
            eta = (W4_TOTAL - idx) / max(rate, 1.0) / 3600
            print(f"idx {idx}/{W4_TOTAL} hits {len(hits)} rate {rate:.1f}/s "
                  f"ETA {eta:.2f}h", flush=True)
        if since >= CKPT_EVERY or idx == W4_TOTAL:
            append_rec("progress", {"next_index": idx,
                                    "hits_so_far": len(hits),
                                    "last_support": list(comb)})
            append_rec("hits", {"up_to_index": idx, "new_hits": hits})
            since = 0
    append_rec("final", {"next_index": idx, "hits_total": len(hits)})
    finalize(hits)


def finalize(hits: list) -> None:
    (HERE / "census_hits.json").write_text(json.dumps({"hits": hits}, indent=1))
    print(f"census_hits.json written: {len(hits)} nonzero w=4 supports",
          flush=True)


def validate_resume() -> dict:
    ledger = load_valid_ledger()
    idxs = [d["next_index"] for d in ledger if "next_index" in d]
    print(json.dumps({"valid_lines": len(ledger),
                      "last_index": max(idxs, default=None)}, sort_keys=True))
    return {"valid_lines": len(ledger), "last_index": max(idxs, default=None)}


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "--run"
    if mode == "--anchors-only":
        a = anchors()
        c = crosscheck_unvectorized()
        w = w3_reanchor_and_w4_sanity()
        print(json.dumps({"anchors": a["n"], "crosscheck": c, "sanity": w},
                         indent=1, default=str))
    elif mode == "--kill-demo":
        KILL_AT = 410_000
        I, eng, B = build_instruments(Q, S_PN4, T_PN4)
        r = B.shape[0]
        hits, idx = [], 0
        t0 = time.perf_counter()
        for comb in itertools.combinations(range(N_PN4), 4):
            mask = np.zeros(N_PN4, dtype=bool)
            mask[list(comb)] = True
            off = np.nonzero(~mask)[0]
            if rank_mod_np(B[:, off], Q, r) < r:
                hits.append({"S": list(comb), "idx": idx,
                             "dim": r - rank_mod_np(B[:, off], Q, r)})
            idx += 1
            if idx == KILL_AT:
                rate = idx / (time.perf_counter() - t0)
                append_rec("progress", {"next_index": idx,
                                        "hits_so_far": len(hits),
                                        "last_support": list(comb),
                                        "note": "kill-demo checkpoint"})
                append_rec("hits", {"up_to_index": idx, "new_hits": hits})
                print(f"KILL-DEMO: checkpointed idx {idx}, hits {len(hits)}, "
                      f"rate {rate:.1f}/s -> hard exit(137) now", flush=True)
                os._exit(137)
        raise SystemExit("kill-demo passed its kill point without killing?!")
    elif mode == "--validate-resume":
        validate_resume()
    elif mode == "--run":
        census_w4("--run")
    else:
        raise SystemExit(f"unknown mode {mode}")
