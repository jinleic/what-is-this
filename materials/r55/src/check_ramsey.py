#!/usr/bin/env python3
"""Independent validator for Ramsey graph catalogs in graph6 format.

A Ramsey(s,t,n)-graph is a graph on n vertices with no clique of size s and
no independent set of size t.

This file is part of the trusted base of the R(5,5) campaign: it is written
directly from the graph6 format specification (nauty formats.txt) and imports
no graph library. Keep it small enough to audit by hand. `--selftest` checks
the parser by encode/decode roundtrip and the clique search against a
brute-force subset enumeration on random small graphs.
"""

import argparse
import itertools
import json
import random
import sys
from multiprocessing import Pool


# popcount: int.bit_count needs Python >= 3.10; 16-bit table works everywhere
try:
    popcount = int.bit_count
except AttributeError:
    _PC16 = bytes(bin(i).count("1") for i in range(1 << 16))

    def popcount(x):
        c = 0
        while x:
            c += _PC16[x & 0xFFFF]
            x >>= 16
        return c


# ---------------------------------------------------------------- graph6 ---

def parse_graph6_line(line):
    """Decode one graph6 line -> (n, adj) with adj a list of int bitmasks."""
    s = line.strip()
    if s.startswith(">>graph6<<"):
        s = s[10:]
    data = [ord(c) - 63 for c in s]
    if not data or any(b < 0 or b > 63 for b in data):
        raise ValueError(f"invalid graph6 bytes in {s!r}")
    if data[0] < 63:
        n, idx = data[0], 1
    elif len(data) > 1 and data[1] < 63:
        n, idx = (data[1] << 12) | (data[2] << 6) | data[3], 4
    else:
        n = 0
        for k in range(2, 8):
            n = (n << 6) | data[k]
        idx = 8
    adj = [0] * n
    need = n * (n - 1) // 2
    bits = 0
    i, j = 0, 1  # upper triangle in column order: (0,1),(0,2),(1,2),(0,3),...
    for b in data[idx:]:
        for shift in (5, 4, 3, 2, 1, 0):
            if bits == need:
                break
            if (b >> shift) & 1:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
            bits += 1
            i += 1
            if i == j:
                i, j = 0, j + 1
    if bits < need:
        raise ValueError(f"graph6 line too short for n={n}: {s!r}")
    return n, adj


def encode_graph6(n, adj):
    """Inverse of parse_graph6_line, for roundtrip self-testing (n <= 62)."""
    assert n <= 62
    out = [n + 63]
    acc, nb = 0, 0
    for j in range(1, n):
        for i in range(j):
            acc = (acc << 1) | ((adj[i] >> j) & 1)
            nb += 1
            if nb == 6:
                out.append(acc + 63)
                acc, nb = 0, 0
    if nb:
        out.append((acc << (6 - nb)) + 63)
    return "".join(chr(c) for c in out)


# ----------------------------------------------------------- clique test ---

def has_clique(adj, cand, k):
    """True iff the graph has a k-clique inside the candidate bitmask."""
    if k == 0:
        return True
    while cand:
        if popcount(cand) < k:
            return False
        v = (cand & -cand).bit_length() - 1
        cand &= cand - 1
        if has_clique(adj, cand & adj[v], k - 1):
            return True
    return False


def check_ramsey_graph(n, adj, s, t):
    """Return (no_Ks, no_It) for the graph."""
    full = (1 << n) - 1
    comp = [full & ~adj[v] & ~(1 << v) for v in range(n)]
    return (not has_clique(adj, full, s), not has_clique(comp, full, t))


# ------------------------------------------------------------ file check ---

_WORK = {}


def _init_worker(s, t):
    _WORK["st"] = (s, t)


def _check_line(arg):
    lineno, line = arg
    s, t = _WORK["st"]
    n, adj = parse_graph6_line(line)
    no_ks, no_it = check_ramsey_graph(n, adj, s, t)
    edges = sum(popcount(a) for a in adj) // 2
    degs = [popcount(a) for a in adj]
    return (lineno, n, edges, min(degs), max(degs), no_ks, no_it)


def check_file(path, s, t, workers):
    with open(path) as f:
        lines = [(i + 1, ln) for i, ln in enumerate(f) if ln.strip()]
    agg = {
        "file": path, "s": s, "t": t, "graphs": 0, "violations": [],
        "n_values": {}, "edge_min": None, "edge_max": None,
        "deg_min": None, "deg_max": None,
    }
    with Pool(workers, initializer=_init_worker, initargs=(s, t)) as pool:
        for lineno, n, e, dmin, dmax, no_ks, no_it in pool.imap_unordered(
                _check_line, lines, chunksize=256):
            agg["graphs"] += 1
            agg["n_values"][n] = agg["n_values"].get(n, 0) + 1
            agg["edge_min"] = e if agg["edge_min"] is None else min(agg["edge_min"], e)
            agg["edge_max"] = e if agg["edge_max"] is None else max(agg["edge_max"], e)
            agg["deg_min"] = dmin if agg["deg_min"] is None else min(agg["deg_min"], dmin)
            agg["deg_max"] = dmax if agg["deg_max"] is None else max(agg["deg_max"], dmax)
            if not (no_ks and no_it):
                agg["violations"].append(
                    {"line": lineno, "has_K_s": not no_ks, "has_I_t": not no_it})
    agg["all_ramsey"] = not agg["violations"]
    return agg


# -------------------------------------------------------------- selftest ---

def _brute_max_clique(n, adj):
    for k in range(n, 0, -1):
        for sub in itertools.combinations(range(n), k):
            if all(adj[u] >> v & 1 for u, v in itertools.combinations(sub, 2)):
                return k
    return 0


def selftest(trials=400, seed=20260813):
    rng = random.Random(seed)
    for trial in range(trials):
        n = rng.randint(1, 10)
        adj = [0] * n
        for u in range(n):
            for v in range(u + 1, n):
                if rng.random() < rng.choice((0.2, 0.5, 0.8)):
                    adj[u] |= 1 << v
                    adj[v] |= 1 << u
        # parser roundtrip
        n2, adj2 = parse_graph6_line(encode_graph6(n, adj))
        assert (n2, adj2) == (n, adj), f"roundtrip failed on trial {trial}"
        # clique search vs brute force, on graph and complement
        full = (1 << n) - 1
        comp = [full & ~adj[v] & ~(1 << v) for v in range(n)]
        for A in (adj, comp):
            w = _brute_max_clique(n, A)
            assert has_clique(A, full, w) or w == 0
            assert not has_clique(A, full, w + 1)
    # C5 is a Ramsey(3,3,5)-graph
    c5 = [0] * 5
    for u in range(5):
        v = (u + 1) % 5
        c5[u] |= 1 << v
        c5[v] |= 1 << u
    assert check_ramsey_graph(5, c5, 3, 3) == (True, True)
    # Paley(17): quadratic residues mod 17 -> the unique Ramsey(4,4,17)-graph
    qr = {pow(x, 2, 17) for x in range(1, 17)}
    p17 = [0] * 17
    for u in range(17):
        for v in range(u + 1, 17):
            if (u - v) % 17 in qr or (v - u) % 17 in qr:
                p17[u] |= 1 << v
                p17[v] |= 1 << u
    assert check_ramsey_graph(17, p17, 4, 4) == (True, True)
    assert check_ramsey_graph(17, p17, 4, 5) == (True, True)   # no I4 implies no I5
    assert check_ramsey_graph(17, p17, 3, 3) == (False, False)  # has K3 and I3
    print(f"selftest OK ({trials} random graphs + C5 + Paley17)")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--file")
    ap.add_argument("-s", type=int, help="forbidden clique size")
    ap.add_argument("-t", type=int, help="forbidden independent set size")
    ap.add_argument("--expect", type=int, help="expected graph count")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    if not (args.file and args.s and args.t):
        ap.error("--file, -s, -t required unless --selftest")
    agg = check_file(args.file, args.s, args.t, args.workers)
    if args.expect is not None:
        agg["expected"] = args.expect
        agg["count_matches"] = agg["graphs"] == args.expect
    print(json.dumps(agg, indent=1, sort_keys=True))
    ok = agg["all_ramsey"] and agg.get("count_matches", True)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
