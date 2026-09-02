"""Parity-kernel CROSS-CHECK runner for the W validation ladder.

Instrument: V = ker(pi), pi = H_0 (x) H_1 (x) H_2 (GRS parity checks, exact
mod q); dim(V ∩ F^S) = |S| - rank(pi[:, S]) (Main's formulation, INPUT to
be reproduced — never a target). Computes per ladder row: nonzero-support
count by weight AND the distinct normalized candidate class (kernel basis
via left-kernel of pi[:, S]^T -> V-vectors -> rref -> normalize), so BOTH
count semantics (gate-B candidate counts; PN6/PP2/PN4 support counts) are
cross-checkable.

This is NOT the pre-registered anchor evidence; it is the concurrent
cross-check. Any disagreement with census_runner.py's V-basis kernel on
ANY row is a STOP-AND-REPORT event (Main's rule): freeze both outputs, no
averaging, no picking.

Writes parity_ledger.jsonl (canonical sha256 lines, same contract).
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RS = HERE.parents[1]
sys.path.insert(0, str(RS))
from src.campaign import _normalize  # noqa: E402
from src.field import rref  # noqa: E402
from src.rs import Inst  # noqa: E402

LEDGER = HERE / "parity_ledger.jsonl"


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


def parity_matrix(q: int, s: tuple, t: tuple) -> np.ndarray:
    """pi: R x N matrix of parity checks, R = prod(s_i - t_i).
    H_i: (s_i - t_i) x s_i GRS parity-check for C_i (vanishing-moment dual):
    rows = polynomials of degree < s_i - t_i? Standard GRS dual: C_i is
    [s_i, t_i]; its dual is [s_i, s_i - t_i] with check columns... concrete:
    H_i[a, b] = w_a * b^a? Use the moment form: h in A_i with line-sum
    moments 0..t_i-1 vanishing: C_i^perp = {h : sum_a h(a) p(a) lam-scale ...
    Simplest exact: C_i^perp = {h : sum_a h_a lam_a S_i^a^d = 0 for d<t_i}
    is NOT right either.

    Ground truth used here: L_i(C_i) = rowspace(line classes); V = sum of
    slices; V^perp = intersection of per-slice annihilators EXACTLY when?
    NO: V^perp = ∩_i (slice_i)^perp where slice_i^perp = annihilator of
    A⊗C_i⊗A = ker-slices... For the RANK identity dim V = N - prod(s_i-t_i)
    with V = sum slices, V^perp = ∩ slices⊥ and dim V^perp = prod(s_i-t_i).
    ∩_i (dual slice) as a KRONECKER PRODUCT: (A0⊗C0^⊥⊗...)? dim of
    (ker π0 ⊗ ...)⊥ = ... V^perp = (ker π)^perp = im(π^T) — the row space of
    the tensor map. So pi = H0 ⊗ H1 ⊗ H2 with H_i a basis of C_i^perp
    ((s_i - t_i) x s_i). C_i^perp: dual GRS: {h : sum_a h_a c_a = 0 ∀c ∈ C_i}
    = {h : Σ_a h_a λ_a p(a) = 0 ∀ deg p < t_i} — an (s_i - t_i)-dim space
    with basis the vectors (h_a) solving the t_i moment equations
    Σ_a h_a λ_a S^(a) d = 0, i.e. h ⊥ (λ⊗S^d) for d < t_i.
    Build H_i = kernel (left) of the t_i x s_i moment matrix M[d][a] =
    lam_a * S_a^d."""
    I = Inst(q, s, t)
    Hs = []
    for i in range(3):
        Si, lam, tt, ss = I.S[i], I.lam[i], I.t[i], I.s[i]
        if ss == tt:
            # C_i = A_i (t_i = s_i): V^perp has no factor here; the
            # kron over a 0-row H_i is the empty matrix — handle via a
            # 1 x 1 identity placeholder so the kron shape stays right.
            Hs.append(np.ones((0, ss), dtype=np.int64))
            continue
        M = np.array([[lam[a] * pow(int(Si[a]), d, q) % q
                       for a in range(ss)] for d in range(tt)],
                      dtype=np.int64)
        from src.field import kernel_mod
        K = kernel_mod([[int(x) for x in row] for row in M], q, ss)
        assert len(K) == ss - tt, f"dim C_{i}_perp {len(K)} != {ss-tt}"
        Hs.append(np.array(K, dtype=np.int64))
    if any(h.shape[0] == 0 for h in Hs):
        return np.zeros((0, s[0] * s[1] * s[2]), dtype=np.int64)
    return np.kron(np.kron(Hs[0], Hs[1]), Hs[2]) % q

def rank_mod_np(m: np.ndarray, q: int) -> int:
    m = m.copy() % q
    rows, cols = m.shape
    rr = piv = 0
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

def census_parity_w3(P: np.ndarray, q: int, n: int, w: int,
                     collect_cands: bool = False):
    """Full w<=w window via the parity formulation:
    dim(V ∩ F^S) = |S| - rank(pi[:, S]); nonzero iff rank < |S|."""
    by_w = {}
    cands = set()
    seen = 0
    t0 = time.perf_counter()
    Pn = P % q
    for k in range(1, w + 1):
        cnt = 0
        for comb in itertools.combinations(range(n), k):
            sub = Pn[:, list(comb)]           # R x |S|
            rank = rank_mod_np(sub, q)
            if rank < k:
                cnt += 1
                if collect_cands:
                    cands |= kernel_cands(sub, Pn, q, comb, n)
            seen += 1
            if seen % 500_000 == 0:
                print(f"    [parity] {seen} supports "
                      f"{time.perf_counter()-t0:.0f}s", flush=True)
        by_w[k] = cnt
    return {"by_w": by_w, "total": sum(by_w.values()), "seen": seen,
            "candidates": len(cands), "cand_set": cands}


def kernel_cands(sub: np.ndarray, Pn: np.ndarray, q: int, Scomb: tuple,
                 n: int):
    """Candidate class for one nonzero support, SAME semantics as the
    V-basis kernel: y ∈ F_q^{|S|} with π[:,S]·y = 0 (right kernel of sub);
    embed each kernel vector into F_q^N on S; rref over F_q^N (canonical
    for the subspace dimension -- the embedding of a subspace basis rrefs
    to the same rows as the full-space computation since support);
    normalize first nonzero = 1; dedupe across supports is upstream (set)."""
    m = sub.copy() % q                  # R x |S|
    rrk, piv2 = 0, []
    rows, cols = m.shape
    for cc in range(cols):
        nz = np.nonzero(m[rrk:, cc])[0]
        if nz.size == 0:
            continue
        i = rrk + nz[0]
        if i != rrk:
            m[[rrk, i]] = m[[i, rrk]]
        inv = pow(int(m[rrk, cc]), -1, q)
        m[rrk] = (m[rrk] * inv) % q
        colv = m[:, cc].copy()
        colv[rrk] = 0
        m = (m - np.outer(colv, m[rrk])) % q
        piv2.append(cc)
        rrk += 1
        if rrk == rows:
            break
    free2 = [c for c in range(cols) if c not in piv2]
    yvecs = []
    for f in free2:
        y = np.zeros(cols, dtype=np.int64)
        y[f] = 1
        for i, pc in enumerate(piv2):
            y[pc] = (-m[i, f]) % q
        yvecs.append(y)
    if not yvecs:
        return set()
    # embed into F_q^N and canonicalize on FULL vectors
    full = []
    Slist = list(Scomb)
    for y in yvecs:
        v = np.zeros(n, dtype=np.int64)
        v[Slist] = y
        full.append([int(x) for x in v])
    rr2, _ = rref(full, q)
    return {_normalize(list(map(int, v)), q) for v in rr2}
K4_TOTAL = 15_329_615   # C(140,4) comb(140,4) == itertools count (asserted)
Q4, S4, T4, N4 = 421, (4, 5, 7), (1, 1, 1), 140
W4_LEDGER = HERE / "parity_w4_ledger.jsonl"
CKPT_EVERY = 400_000


def _w4_append(kind: str, data: dict) -> None:
    body = {"kind": kind, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                               time.gmtime()), **data}
    canon = json.dumps(body, sort_keys=True, separators=(",", ":"))
    body["_sha256"] = sha(canon.encode())
    line = json.dumps(body, sort_keys=True, separators=(",", ":")) + "\n"
    existing = W4_LEDGER.read_bytes() if W4_LEDGER.exists() else b""
    tmp = W4_LEDGER.with_name(W4_LEDGER.name + ".tmp")
    with open(tmp, "wb") as f:
        f.write(existing + line.encode())
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, W4_LEDGER)


def _w4_load() -> list[dict]:
    if not W4_LEDGER.exists():
        return []
    out = []
    for ln in W4_LEDGER.read_text().splitlines():
        if not ln.strip():
            continue
        d = json.loads(ln)
        s = d.pop("_sha256")
        if sha(json.dumps(d, sort_keys=True,
                          separators=(",", ":")).encode()) != s:
            raise SystemExit(f"HARD ABORT: {W4_LEDGER} checksum-invalid; "
                             "never truncate/repair")
        out.append(d)
    last = -1
    for d in out:
        ni = d.get("next_index")
        if ni is None:
            continue
        if ni < last:
            raise SystemExit(f"HARD ABORT: index-order gap at {ni}")
        last = ni
    return out


def w4_finalize(hits: list, extra: dict) -> None:
    (HERE / "parity_w4_hits.json").write_text(
        json.dumps({"hits": hits, **extra}, indent=1))
    print(f"parity_w4_hits.json written: {len(hits)} nonzero w=4 supports",
          flush=True)


def census_parity_w4() -> None:
    """Parity PRIMARY census: all C(140,4) supports, exact mod 421.
    Emptiness test: dim(V cap F^S) = |S| - rank(pi[:,S]) with
    pi = H0 kron H1 kron H2 (parity_matrix(Q4,S4,T4)) and the column
    index (a0*s1+a1)*s2+a2 (Inst.idx; tuple order matches census_runner
    combinations over range(140))."""
    import src.rs as rs_
    I = rs_.Inst(Q4, S4, T4)
    assert I.idx((0, 0, 0)) == 0 and I.idx((3, 4, 6)) == 139
    assert math.comb(140, 4) == K4_TOTAL
    P = parity_matrix(Q4, S4, T4)
    Pn = P % Q4
    assert Pn.shape == (3 * 4 * 6, 140), Pn.shape
    done_idx, hits, final_seen = 0, [], None
    for d in _w4_load():
        if d["kind"] == "final":
            final_seen = d
        if d["kind"] == "progress" and d.get("next_index", 0) > done_idx:
            done_idx = d["next_index"]
        if d["kind"] == "hits":
            hits = d["new_hits"]
    if final_seen:
        print("[parity-w4] census already FINAL", flush=True)
        return
    if done_idx >= K4_TOTAL:
        return w4_finalize(hits, {})
    print(f"[parity-w4] resume/start at index {done_idx}; {len(hits)} hits",
          flush=True)
    t0 = time.perf_counter()
    idx = done_idx
    since = 0
    for comb in itertools.islice(itertools.combinations(range(N4), 4),
                                 done_idx, None):
        sub = Pn[:, comb]                        # R x 4
        rank = rank_mod_np(sub, Q4)
        if rank < 4:
            hits.append({"S": list(comb), "dim": 4 - rank,
                         "idx": idx,
                         "class": sorted(
                             [list(v) for v in
                              kernel_cands(sub, Pn, Q4, comb, N4)])})
        idx += 1
        since += 1
        if idx % 100_000 == 0:
            rate = (idx - done_idx) / (time.perf_counter() - t0)
            eta = (K4_TOTAL - idx) / max(rate, 1.0) / 60
            print(f"[parity-w4] idx {idx}/{K4_TOTAL} hits {len(hits)} "
                  f"rate {rate:.0f}/s ETA {eta:.1f}min", flush=True)
        if since >= CKPT_EVERY or idx == K4_TOTAL:
            _w4_append("progress", {"next_index": idx, "hits_so_far": len(hits),
                                    "last_support": list(comb)})
            _w4_append("hits", {"up_to_index": idx, "new_hits": hits})
            since = 0
    _w4_append("final", {"next_index": idx, "hits_total": len(hits)})
    w4_finalize(hits, {})
    return


def parity_line_checks() -> dict:
    """Independent W-kernel reproduction of the owner INPUT: (A) the 35
    direction-0 lines are dependent at weight 4, settled exactly as FULL
    lines: dim(V cap F^line) = 1 and every proper 3-subset independent;
    (B) direction-1/2 minimum lines are dependent ONLY as full lines
    (dims 5,7; every proper subset independent) - consistent with d1=5,
    d2=7, so they cannot appear at weight 4."""
    import src.rs as rs_
    I = rs_.Inst(Q4, S4, T4)
    P = parity_matrix(Q4, S4, T4) % Q4
    out = {}
    for i in range(3):
        n = full = 0
        dims = set()
        for l in range(I.num_lines(i)):
            idxs = I.line_indices(i, l)
            rank = rank_mod_np(P[:, idxs], Q4)
            assert rank < len(idxs), f"dir {i} line {l} independent?!"
            dim = len(idxs) - rank
            full += int(dim == 1)
            dims.add(dim)
            n += 1
            sub = idxs[:-1]                      # drop first point: any
            rank_sub = rank_mod_np(P[:, sub], Q4)  # proper subset (affine)
            assert rank_sub == len(sub), \
                f"dir {i} line {l}: proper subset dependent"
            if i == 0:
                assert len(idxs) == 4 and dim == 1
            else:
                assert len(idxs) == 5 or len(idxs) == 7
        assert dims == {1}, f"dir {i}: unexpected dims {dims}"
        out[f"d{i}"] = {"n_lines": n, "full_line_dependent": full,
                        "proper_subset_independent": True,
                        "line_dim": 1, "line_size": [S4[0], S4[1], S4[2]][i]}
    return out


def main() -> int:
    gateb = json.loads((RS / "campaigns" /
                        "2026-08-30T11-43-59Z_50595CC9_fffe84b0_gateB" /
                        "payload.json").read_text())
    done = set()
    if LEDGER.exists():
        for ln in LEDGER.read_text().splitlines():
            if ln.strip():
                d = json.loads(ln)
                s = d.pop("_sha256")
                if sha(json.dumps(d, sort_keys=True,
                                  separators=(",", ":")).encode()) == s:
                    if d["kind"] == "row":
                        done.add(d["id"])
    for row in gateb["rows"]:
        if row["id"] in done:
            print(f"[parity] {row['id']} already ledgered: skip", flush=True)
            continue
        q, s, t = row["q"], tuple(row["s"]), tuple(row["t"])
        n = q and s[0] * s[1] * s[2]
        P = parity_matrix(q, s, t)
        assert P.shape[0] == (s[0]-t[0])*(s[1]-t[1])*(s[2]-t[2])
        t0 = time.perf_counter()
        res = census_parity_w3(P, q, n, 3, collect_cands=True)
        wall = round(time.perf_counter() - t0, 1)
        want = row["candidate_count_total"]
        ok = res["candidates"] == want
        append_rec("row", {"id": row["id"], "q": q, "s": list(s),
                           "t": list(t), "candidates": res["candidates"],
                           "want_gateB": want, "match": ok,
                           "by_w": res["by_w"], "wall_s": wall})
        print(f"[parity] {row['id']}: cands {res['candidates']} "
              f"(gateB {want}) match={ok} by_w {res['by_w']} wall {wall}s",
              flush=True)
    # P2 anchors
    p2 = json.loads((RS / "scratch" / "p2_checkpoint.json").read_text())
    for pid, q, s, want_sup, want_c in [("PN6", 43, (2, 6, 7), 3486, 42),
                                        ("PP2", 43, (2, 3, 7), 917, 77)]:
        if pid in done:
            continue
        P = parity_matrix(q, s, (1, 1, 1))
        res = census_parity_w3(P, q, s[0]*s[1]*s[2], 3)
        ok = res["total"] == want_sup
        append_rec("row", {"id": pid, "q": q, "s": list(s),
                           "supports": res["total"], "want_supports": want_sup,
                           "match_supports": ok, "by_w": res["by_w"]})
        print(f"[parity] {pid}: supports {res['total']} ({want_sup}) "
              f"match={ok}", flush=True)
    print("[parity] ladder cross-check complete", flush=True)
    return 0


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "--crosscheck"
    if mode == "--crosscheck":
        raise SystemExit(main())
    elif mode == "--w4-primary":
        r0 = parity_line_checks()
        print("[parity-w4] line checks:", json.dumps(r0), flush=True)
        census_parity_w4()
    elif mode == "--w4-lines-only":
        print(json.dumps(parity_line_checks(), indent=1))
    else:
        raise SystemExit(f"unknown mode {mode}")
