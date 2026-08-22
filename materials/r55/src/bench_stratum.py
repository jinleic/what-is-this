#!/usr/bin/env python3
"""Frozen gate-3 benchmark harness: reproduce one edge-extremal stratum
R(4,5,n,e) end-to-end and diff canonically against the published census.

Method (audited spec, same trust base as the D3 gate-2 run):
  * Every census graph G has a max-degree vertex v with deg(v) = d; N(v) is
    triangle-free and I5-free (else K4/I5 with v)  => N(v) in R(3,5,d), d<=13;
    the non-neighborhood is K4-free and I4-free    => in R(4,4,q), q=n-1-d<=17.
  * Min-degree bound (theorem): deleting a min-degree vertex leaves an
    R(4,5,n-1) graph with e - delta <= E(4,5,n-1), so
    delta >= dmin_true = e - E(4,5,n-1).  E derived from the validated archive.
  * For each feasible d: glue every (H,K) pair from the validated catalogs
    with the fixed C engine (glue_census2_c, GLUE_DMIN=dmin_true); pairs that
    exceed the node budget are deferred and resolved exactly by sat_pair
    (AllSAT + trusted checker).  q=0 slices are refuted by Mantel's bound.
  * Every emitted graph is re-verified here (edge count, K4-free, I5-free,
    degree window) before entering the raw file.
  * Only then is the published member unsealed, both sides canonicalized with
    labelg, and the sets compared.  Any mismatch => exit 1 (fail closed).

Usage: bench_stratum.py <n> <e> [--workers 8] [--node-limit 1000000]
                                [--outdir OUTD/bench_strata]
Prints one BENCH-JSON line with the frozen metrics.
"""
import argparse
import json
import os
import subprocess
import sys
import tarfile
import time
from concurrent.futures import ProcessPoolExecutor

SRC = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SRC)
BASE = os.path.dirname(SRC)
DATA = os.path.join(BASE, "data")
G2 = os.path.join(SRC, "glue_census2_c")
LABELG = "/opt/homebrew/bin/labelg"

from check_ramsey import parse_graph6_line, check_ramsey_graph  # noqa: E402
import sat_pair  # noqa: E402

# E(4,5,m): derived from the validated r45extreme archive (gate 1), not memory.
E45 = {10: 33, 11: 40, 12: 48, 13: 53, 14: 60, 15: 66, 16: 72, 17: 79,
       18: 85, 19: 92, 20: 100, 21: 107, 22: 114, 23: 122}


def parse_g6_bytes(line):
    return parse_graph6_line(line.decode() if isinstance(line, bytes) else line)


def load_catalog(path):
    return [l for l in open(path, "rb").read().splitlines() if l.strip()]


def verify(n, e, dmin, d, adj):
    degs = [bin(a).count("1") for a in adj]
    assert sum(degs) // 2 == e, "edge count"
    assert min(degs) >= dmin, "min degree"
    assert max(degs) == d, "max degree"
    assert check_ramsey_graph(n, adj, 4, 5) == (True, True), "Ramsey property"


def sat_chunk(args):
    """Resolve deferred pairs exactly. Yields (g6_string, None) per solution
    and one trailing (None, stats) marker with the solver's accum stats."""
    n, e, d, dmin, pairs = args
    out, stats = [], {}
    for H_adj, K_adj, key in pairs:
        for g6s in sat_pair.solve_pair(n, e, d, H_adj, K_adj, key, dmin,
                                       stats=stats):
            out.append((g6s, None))
    out.append((None, stats))
    return out


def run_stratum(n, e, workers, node_limit, outdir):
    dmin = e - E45[n - 1]
    dlo = max(-(-2 * e // n), dmin)
    dhi = min(13, n - 1)
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    env = dict(os.environ, GLUE_DMIN=str(dmin), GLUE_NODE_LIMIT=str(node_limit))
    raw, metrics = [], {}
    for d in range(dlo, dhi + 1):
        q = n - 1 - d
        if q > 17:
            continue
        if q < 1:
            # root adjacent to everything: eH = e-d on d triangle-free vertices
            assert e - d > d * d // 4, "q=0 slice not refuted by Mantel"
            metrics[f"d{d}"] = {"pairs": 0, "note": "q=0 refuted by Mantel"}
            continue
        hp = os.path.join(DATA, f"r35_{d}.g6")
        kp = os.path.join(DATA, f"r44_{q}.g6")
        assert os.path.exists(hp) and os.path.exists(kp), (hp, kp)
        hl, kl = load_catalog(hp), load_catalog(kp)
        og = os.path.join(outdir, f"b{n}_{e}_d{d}.g6")
        oc = os.path.join(outdir, f"b{n}_{e}_d{d}.csv")
        # slice the K file across workers
        step = max(1, len(kl) // (workers * 3) or 1)
        jobs = [(klo, min(klo + step, len(kl)))
                for klo in range(0, len(kl), step)]
        procs, done, deferred, td = [], 0, [], time.time()
        ji = 0
        outfiles = []
        while ji < len(jobs) or procs:
            while ji < len(jobs) and len(procs) < workers:
                klo, khi = jobs[ji]
                o1, o2 = f"{og}.{ji}", f"{oc}.{ji}"
                p = subprocess.Popen(
                    [G2, str(n), str(e), str(d), hp, kp, "0", str(len(hl)),
                     o1, o2, str(klo), str(khi)], env=env)
                procs.append((p, o1, o2))
                outfiles.append((o1, o2))
                ji += 1
            nxt = []
            for p, o1, o2 in procs:
                if p.poll() is None:
                    nxt.append((p, o1, o2))
                    continue
                assert p.returncode == 0, (o1, p.returncode)
                done += 1
            procs = nxt
            time.sleep(0.02)
        nsol_dfs, dfs_nodes = 0, 0
        for o1, o2 in outfiles:
            for line in open(o1, "rb").read().splitlines():
                if line.strip():
                    _, adj = parse_g6_bytes(line)
                    verify(n, e, dmin, d, adj)
                    raw.append(line.decode())
                    nsol_dfs += 1
            for row in open(o2):
                if not row.strip():
                    continue
                if row.startswith("#nodes,"):
                    dfs_nodes += int(row.split(",")[1])
                    continue
                _, h_i, k_i, ns = (int(x) for x in row.split(",")[:4])
                if ns == -1:
                    deferred.append((h_i, k_i))
        # SAT tail: exact resolution of deferred pairs
        nsol_sat = 0
        sat_stats = {}
        if deferred:
            tasks = []
            for h_i, k_i in deferred:
                _, H_adj = parse_g6_bytes(hl[h_i])
                _, K_adj = parse_g6_bytes(kl[k_i])
                tasks.append((H_adj, K_adj, (d, h_i, k_i)))
            csz = max(1, len(tasks) // (workers * 4) or 1)
            work = [(n, e, d, dmin, tasks[i:i + csz])
                    for i in range(0, len(tasks), csz)]
            with ProcessPoolExecutor(max_workers=workers) as ex:
                for part in ex.map(sat_chunk, work):
                    for g6s, st in part:
                        if g6s is None:
                            for key, v in st.items():
                                sat_stats[key] = sat_stats.get(key, 0) + v
                            continue
                        _, adj = parse_g6_bytes(g6s.encode())
                        verify(n, e, dmin, d, adj)
                        raw.append(g6s)
                        nsol_sat += 1
        metrics[f"d{d}"] = {
            "pairs": len(hl) * len(kl), "deferred": len(deferred),
            "dfs_raw": nsol_dfs, "sat_raw": nsol_sat,
            "dfs_nodes": dfs_nodes,
            "sat_conflicts": sat_stats.get("conflicts", 0),
            "sat_decisions": sat_stats.get("decisions", 0),
            "wall_s": round(time.time() - td, 1)}
        print(f"  ({n},{e}) d={d}: pairs={len(hl)*len(kl)} "
              f"deferred={len(deferred)} raw={nsol_dfs}+{nsol_sat} "
              f"[{time.time()-td:.1f}s]", flush=True)
    return raw, metrics, dmin, time.time() - t0


def canonical_diff(n, e, raw, outdir):
    raw_path = os.path.join(outdir, f"b{n}_{e}.raw.g6")
    with open(raw_path, "w") as f:
        f.write("\n".join(raw) + ("\n" if raw else ""))
    canon = raw_path + ".canon"
    subprocess.run([LABELG, "-q", raw_path, canon], check=True)
    mine = sorted({l for l in open(canon).read().splitlines() if l.strip()})
    tf = tarfile.open(os.path.join(DATA, "r45extreme.tar.gz"))
    pub_raw = os.path.join(outdir, f"pub_{n}_{e}.g6")
    open(pub_raw, "wb").write(
        tf.extractfile(f"r45extreme/r45{n}.{e}.g6").read())
    pub_canon = pub_raw + ".canon"
    subprocess.run([LABELG, "-q", pub_raw, pub_canon], check=True)
    pub = sorted({l for l in open(pub_canon).read().splitlines() if l.strip()})
    return mine, pub


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("e", type=int)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--node-limit", type=int, default=1000000)
    ap.add_argument("--outdir",
                    default=os.path.join(BASE, "outD", "bench_strata"))
    a = ap.parse_args()
    raw, metrics, dmin, wall = run_stratum(
        a.n, a.e, a.workers, a.node_limit, a.outdir)
    mine, pub = canonical_diff(a.n, a.e, raw, a.outdir)
    ok = bool(mine) and mine == pub
    rec = {"n": a.n, "e": a.e, "dmin_true": dmin, "workers": a.workers,
           "node_limit": a.node_limit, "raw": len(raw), "classes": len(mine),
           "published": len(pub), "match": ok, "wall_s": round(wall, 1),
           "per_d": metrics}
    print("BENCH-JSON " + json.dumps(rec), flush=True)
    print(f"BENCH RESULT ({a.n},{a.e}):", "MATCH" if ok else "MISMATCH",
          flush=True)
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
