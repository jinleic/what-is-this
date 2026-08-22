#!/usr/bin/env python3
"""Min-degree case-split census of R(4,5,21,e=107) (the AM 'D3' stratum).

Decomposition (proofs; delta = min degree, Delta = max degree of a census graph G):
  * delta >= 7: deleting a min-degree vertex leaves an R(4,5,20) graph with
    107-delta <= E(4,5,20) = 100 edges (order-20 census, gate-1 validated).
  * delta <= 10: min <= mean = 214/21 < 10.2.
  * CASE 1 (delta in {7,8,9}): G-u lies in the published complete classes
    r45extreme/r4520.{100,99,98}.g6, so G is a one-vertex extension of one of
    those 1+822+304,848 graphs by a new vertex u with |N(u)| = 107 - e'.
    Extension constraints (exact iff, given G-u in R(4,5,20)): N(u) contains
    no triangle of G-u; N(u) hits every I4 of G-u.
  * CASE 2 (delta >= 10): the audited max-degree gluing spec applies with the
    within-case prune GLUE_DMIN=10. Delta = d in {11,12,13} (d >= ceil(2e/n) = 11,
    d <= R(3,5)-1 = 13). Within the case, summing per-vertex cone lower bounds:
      d=11: e(B)+3 <= e(A) <= e(B)+6   (else no valid cone exists)
      d=12: e(B)+13 <= e(A) <= e(B)+15
      d=13: forces e(B) <= 3 < 5 = min e over R(4,4,7)  ->  case EMPTY.
    Cases are disjoint (min degree <= 9 vs >= 10); their union is the census.

Engines: glue_census2_c (frozen, audited) with GLUE_DMIN=10, node budget 1e6;
pairs exceeding the budget go to sat_pair.solve_pair (same dmin). Every output
graph is re-verified here (K4-free, I5-free, order, edge count, case min-degree
invariant) with an independent checker before entering the raw file.

Output: outD/case_split_21_107/ (isolated; nothing outside is touched).
Final step canonicalizes with labelg and only then diffs against the held-out
published file.
"""
import os
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

SRC = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(SRC, "..")
DATA = os.path.join(BASE, "data")
OUT = os.path.join(BASE, "outD", "case_split_21_107")
G2 = os.path.join(SRC, "glue_census2_c")
LABELG = "/opt/homebrew/bin/labelg"
WORKERS = 8
NODE_LIMIT = "1000000"

sys.path.insert(0, SRC)
from check_ramsey import parse_graph6_line, check_ramsey_graph  # noqa: E402
import sat_pair  # noqa: E402

from pysat.card import CardEnc, EncType  # noqa: E402
from pysat.formula import IDPool  # noqa: E402
from pysat.solvers import Cadical195  # noqa: E402


# ---------- shared small-graph helpers (independent of the engines) ----------

def parse_g6_bytes(line):
    n = line[0] - 63
    nbits = n * (n - 1) // 2
    nbytes = (nbits + 5) // 6
    bits = 0
    for c in line[1:1 + nbytes]:
        bits = (bits << 6) | (c - 63)
    bits >>= (6 * nbytes - nbits)
    adj = [0] * n
    pos = nbits - 1
    for j in range(1, n):
        for i in range(j):
            if (bits >> pos) & 1:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
            pos -= 1
    return n, adj


def encode_g6(n, adj):
    bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append((adj[i] >> j) & 1)
    while len(bits) % 6:
        bits.append(0)
    out = [n + 63]
    for i in range(0, len(bits), 6):
        v = 0
        for b in bits[i:i + 6]:
            v = (v << 1) | b
        out.append(v + 63)
    return bytes(out)


def bits_iter(x):
    while x:
        b = x & -x
        yield b.bit_length() - 1
        x ^= b


def complement(n, adj):
    full = (1 << n) - 1
    return [full & ~a & ~(1 << v) for v, a in enumerate(adj)]


def has_k4(n, adj):
    for u in range(n):
        for du in bits_iter(adj[u] >> (u + 1)):
            v = u + 1 + du
            B = adj[u] & adj[v] & ~((1 << (v + 1)) - 1)
            for w in bits_iter(B):
                if B & adj[w] & ~((1 << (w + 1)) - 1):
                    return True
    return False


def has_k5(n, adj):
    for u in range(n):
        for du in bits_iter(adj[u] >> (u + 1)):
            v = u + 1 + du
            B = adj[u] & adj[v] & ~((1 << (v + 1)) - 1)
            for w in bits_iter(B):
                C = B & adj[w] & ~((1 << (w + 1)) - 1)
                for x in bits_iter(C):
                    if adj[x] & C & ~((1 << (x + 1)) - 1):
                        return True
    return False


def k4_list(n, adj):
    out = []
    for u in range(n):
        for du in bits_iter(adj[u] >> (u + 1)):
            v = u + 1 + du
            B = adj[u] & adj[v] & ~((1 << (v + 1)) - 1)
            for w in bits_iter(B):
                C = B & adj[w] & ~((1 << (w + 1)) - 1)
                for x in bits_iter(C):
                    out.append((u, v, w, x))
    return out


def verify21(adj, case):
    n = 21
    e = sum(bin(a).count("1") for a in adj) // 2
    degs = [bin(a).count("1") for a in adj]
    assert e == 107, e
    assert not has_k4(n, adj) and not has_k5(n, complement(n, adj))
    if case == 1:
        assert min(degs) <= 9
    else:
        assert min(degs) >= 10
    return True


# ------------------------------- case 1 -------------------------------------

def case1_chunk(args):
    """Extension AllSAT over one chunk of base graphs.
    Returns (g6 strings, solver accum stats)."""
    lines, delta = args
    out, stats = [], {"conflicts": 0, "decisions": 0}
    for raw in lines:
        n0, g = parse_g6_bytes(raw)
        tris = [(u, v, w) for u in range(20) for v in range(u + 1, 20)
                for w in range(v + 1, 20)
                if (g[u] >> v) & 1 and (g[u] >> w) & 1 and (g[v] >> w) & 1]
        i4s = k4_list(20, complement(20, g))
        pool = IDPool()
        S = [pool.id(("s", i)) for i in range(20)]
        cnf = CardEnc.equals(S, delta, vpool=pool,
                             encoding=EncType.seqcounter).clauses
        for (u, v, w) in tris:
            cnf.append([-S[u], -S[v], -S[w]])
        for q in i4s:
            cnf.append([S[q[0]], S[q[1]], S[q[2]], S[q[3]]])
        with Cadical195(bootstrap_with=cnf) as s:
            while s.solve():
                m = set(x for x in s.get_model() if x > 0)
                sel = [i for i in range(20) if S[i] in m]
                adj = [a << 1 for a in g]          # shift: new vertex = 0
                adj = [0] + [a for a in adj]
                for i, a in enumerate(g):
                    adj[1 + i] = a << 1
                for i in sel:
                    adj[0] |= 1 << (1 + i)
                    adj[1 + i] |= 1
                verify21(adj, 1)
                out.append(encode_g6(21, adj).decode())
                s.add_clause([-S[i] if S[i] in m else S[i] for i in range(20)])
            acc = s.accum_stats()
            for key in stats:
                stats[key] += acc.get(key, 0)
    return out, stats


def run_case1():
    import tarfile
    tf = tarfile.open(os.path.join(DATA, "r45extreme.tar.gz"))
    sols = []
    c1_stats = {"conflicts": 0, "decisions": 0}
    t0 = time.time()
    for eprime, delta in ((100, 7), (99, 8), (98, 9)):
        raws = [l for l in
                tf.extractfile(f"r45extreme/r4520.{eprime}.g6").read().splitlines()
                if l.strip()]
        chunks = [(raws[i:i + 2000], delta) for i in range(0, len(raws), 2000)]
        with ProcessPoolExecutor(max_workers=WORKERS) as ex:
            for part, st in ex.map(case1_chunk, chunks):
                sols.extend(part)
                for key in c1_stats:
                    c1_stats[key] += st.get(key, 0)
        print(f"case1 delta={delta}: {len(raws)} base graphs done, "
              f"cumulative raw solutions={len(sols)}, "
              f"sat_conflicts={c1_stats['conflicts']}, "
              f"sat_decisions={c1_stats['decisions']}, {time.time()-t0:.0f}s",
              flush=True)
    return sols


# ------------------------------- case 2 -------------------------------------

WINDOWS = {11: (3, 6), 12: (13, 15)}   # e(A) - e(B) in [lo, hi]; d=13 empty


def sat_tail_chunk(args):
    """Resolve a chunk of deferred pairs exactly.
    Returns (g6 strings, solver accum stats)."""
    d, pairs = args
    out, stats = [], {}
    for (H_adj, K_adj, key) in pairs:
        out.extend(sat_pair.solve_pair(21, 107, d, H_adj, K_adj, key, 10,
                                       stats=stats))
    return out, stats



def bucket_catalogs(d):
    q = 20 - d
    hpath = os.path.join(DATA, f"r35_{d}.g6")
    kpath = os.path.join(DATA, f"r44_{q}.g6")
    hs, ks = {}, {}
    hidx, kidx = {}, {}
    for path, table, index in ((hpath, hs, hidx), (kpath, ks, kidx)):
        with open(path, "rb") as f:
            for i, line in enumerate(l for l in f.read().splitlines() if l.strip()):
                n0, adj = parse_g6_bytes(line)
                e = sum(bin(a).count("1") for a in adj) // 2
                table.setdefault(e, []).append(line)
                index.setdefault(e, []).append((i, adj))
    return hs, ks, hidx, kidx


def run_case2(d):
    lo, hi = WINDOWS[d]
    q = 20 - d
    hs, ks, hidx, kidx = bucket_catalogs(d)
    tmp = os.path.join(OUT, f"tmp_d{d}_r2")
    os.makedirs(tmp, exist_ok=True)
    env = dict(os.environ, GLUE_DMIN="10", GLUE_NODE_LIMIT=NODE_LIMIT)
    jobs = []
    for ea, hlines in sorted(hs.items()):
        for eb, klines in sorted(ks.items()):
            if not (lo <= ea - eb <= hi):
                continue
            hf = os.path.join(tmp, f"H_{ea}.g6")
            kf = os.path.join(tmp, f"K_{eb}.g6")
            if not os.path.exists(hf):
                open(hf, "wb").write(b"\n".join(hlines) + b"\n")
            if not os.path.exists(kf):
                open(kf, "wb").write(b"\n".join(klines) + b"\n")
            nk = len(klines)
            step = max(1, nk // (WORKERS * 3))
            for ci, klo in enumerate(range(0, nk, step)):
                jobs.append((ea, eb, hf, kf, len(hlines), klo,
                             min(klo + step, nk), ci))
    print(f"case2 d={d}: {len(jobs)} bucket-chunks", flush=True)

    def launch(job):
        ea, eb, hf, kf, nh, klo, khi, ci = job
        og = os.path.join(tmp, f"g_{ea}_{eb}_{ci}.g6")
        oc = os.path.join(tmp, f"g_{ea}_{eb}_{ci}.csv")
        p = subprocess.Popen([G2, "21", "107", str(d), hf, kf, "0", str(nh),
                              og, oc, str(klo), str(khi)], env=env)
        return p, og, oc, ea, eb

    sols, deferred = [], []
    dfs_nodes = [0]
    active, ji = [], 0
    ndone = 0
    t0 = time.time()
    while ji < len(jobs) or active:
        while ji < len(jobs) and len(active) < WORKERS:
            active.append(launch(jobs[ji])); ji += 1
        nxt = []
        for p, og, oc, ea, eb in active:
            if p.poll() is None:
                nxt.append((p, og, oc, ea, eb)); continue
            assert p.returncode == 0, (og, p.returncode)
            for line in open(og, "rb").read().splitlines():
                if line.strip():
                    n0, adj = parse_g6_bytes(line)
                    verify21(adj, 2)
                    sols.append(line.decode())
            for row in open(oc):
                if not row.strip():
                    continue
                if row.startswith("#nodes,"):
                    dfs_nodes[0] += int(row.split(",")[1])
                    continue
                if row.startswith("#"):
                    continue
                dd, h_local, k_local, ns = (int(x) for x in row.split(",")[:4])
                if ns == -1:
                    deferred.append((ea, eb, h_local, k_local))
            ndone += 1
            if ndone % 40 == 0:
                print(f"  d={d}: {ndone}/{len(jobs)} chunks, raw={len(sols)}, "
                      f"deferred={len(deferred)}, {time.time()-t0:.0f}s", flush=True)
        active = nxt
        time.sleep(0.05)
    print(f"case2 d={d}: DFS pass done: raw={len(sols)}, deferred={len(deferred)}, "
          f"dfs_nodes={dfs_nodes[0]}, {time.time()-t0:.0f}s", flush=True)

    # SAT tail: resolve deferred pairs exactly (parallel; solve_pair returns
    # graph6 STRINGS - parse before verifying)
    tasks = []
    for (ea, eb, hl, kl) in deferred:
        tasks.append((hidx[ea][hl][1], kidx[eb][kl][1], (d, ea, hl, eb, kl)))
    chunk_sz = max(1, len(tasks) // (WORKERS * 6) or 1)
    work = [(d, tasks[i:i + chunk_sz]) for i in range(0, len(tasks), chunk_sz)]
    ndone2 = 0
    sat_stats = {}
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for part, st in ex.map(sat_tail_chunk, work):
            for key, v in st.items():
                sat_stats[key] = sat_stats.get(key, 0) + v
            for g6s in part:
                n0, adj = parse_g6_bytes(g6s.encode())
                verify21(adj, 2)
                sols.append(g6s)
            ndone2 += 1
            if ndone2 % 20 == 0 or ndone2 == len(work):
                print(f"  d={d}: SAT {ndone2}/{len(work)} chunks, raw={len(sols)}, "
                      f"{time.time()-t0:.0f}s", flush=True)
    print(f"case2 d={d}: complete: raw={len(sols)}, "
          f"dfs_nodes={dfs_nodes[0]}, sat_conflicts={sat_stats.get('conflicts', 0)}, "
          f"sat_decisions={sat_stats.get('decisions', 0)}, {time.time()-t0:.0f}s",
          flush=True)
    return sols


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    raw = []
    raw += run_case1()
    raw += run_case2(12)
    raw += run_case2(11)
    raw_path = os.path.join(OUT, "r45_21_107.raw.g6")
    with open(raw_path, "w") as f:
        for line in raw:
            f.write(line + "\n")
    print(f"TOTAL raw (verified) solutions: {len(raw)}  [{time.time()-t0:.0f}s]",
          flush=True)
    # canonicalize, dedup, and only now unseal the published reference
    if not final_diff(raw_path):
        sys.exit(1)  # fail closed: automation must not see rc 0 on a bad census


def final_diff(raw_path, expected=31, outdir=None):
    """Canonical diff against the published census. True iff sets match AND
    the published count is the expected one (guards empty/truncated refs).
    `outdir` defaults to OUT; tests MUST pass an isolated directory so the
    official artifacts are never rewritten by exercising the verdict."""
    outdir = outdir or OUT
    os.makedirs(outdir, exist_ok=True)
    canon = os.path.join(outdir, "r45_21_107.canon")
    subprocess.run([LABELG, "-q", raw_path, canon + ".tmp"], check=True)
    lines = sorted({l for l in open(canon + ".tmp").read().splitlines() if l.strip()})
    with open(canon, "w") as f:
        f.write("\n".join(lines) + "\n")
    os.remove(canon + ".tmp")
    import tarfile
    tf = tarfile.open(os.path.join(DATA, "r45extreme.tar.gz"))
    pub_raw = os.path.join(outdir, "pub.g6")
    open(pub_raw, "wb").write(tf.extractfile("r45extreme/r4521.107.g6").read())
    pub_canon = os.path.join(outdir, "pub.canon")
    subprocess.run([LABELG, "-q", pub_raw, pub_canon + ".tmp"], check=True)
    plines = sorted({l for l in open(pub_canon + ".tmp").read().splitlines() if l.strip()})
    with open(pub_canon, "w") as f:
        f.write("\n".join(plines) + "\n")
    os.remove(pub_canon + ".tmp")
    print(f"iso classes: mine={len(lines)} published={len(plines)} "
          f"EQUAL={lines == plines}", flush=True)
    ok = lines == plines and len(plines) == expected  # AM Table 1: 31
    print("D3 RESULT:", "MATCH" if ok else "MISMATCH", flush=True)
    return ok


if __name__ == "__main__":
    main()
