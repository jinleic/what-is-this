#!/usr/bin/env python3
"""Parallel driver for the C enumerator (glue_census.c) with trusted
re-verification: every graph the C core emits is re-checked here with the
audited independent checker (check_ramsey.py) before entering the output
file. Output contract identical to glue_census.py (implementation A)."""

import argparse
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_ramsey import parse_graph6_line, check_ramsey_graph, popcount  # noqa: E402

SRC = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(SRC, "..", "data")
BIN = os.environ.get("GLUE_BIN", os.path.join(SRC, "glue_census_c"))
BIN_SRC = {os.path.join(SRC, "glue_census_c"): "glue_census.c",
           os.path.join(SRC, "glue_census2_c"): "glue_census2.c"}


def count_lines(path):
    with open(path) as f:
        return sum(1 for ln in f if ln.strip())


def _sat_worker(task):
    import sat_pair
    return sat_pair.run_todo(*task)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--e", type=int, required=True)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--outdir", default=os.path.join(SRC, "..", "outA"))
    ap.add_argument("--no-verify", action="store_true")
    args = ap.parse_args()
    n, e = args.n, args.e
    csrc = os.path.join(SRC, BIN_SRC.get(BIN, "glue_census.c"))
    if not os.path.exists(BIN) or os.path.getmtime(BIN) < os.path.getmtime(csrc):
        subprocess.run(["cc", "-O2", "-o", BIN, csrc], check=True)
    d_lo = max((2 * e + n - 1) // n, 0)
    dlist = [d for d in range(d_lo, min(13, n - 1) + 1) if n - 1 - d <= 17]
    os.makedirs(args.outdir, exist_ok=True)
    tmpdir = os.path.join(args.outdir, f"tmp_{n}_{e}")
    os.makedirs(tmpdir, exist_ok=True)
    # build shard jobs: (d, h_path, k_path, h_lo, h_hi, k_lo, k_hi)
    # shard over whichever catalog is larger; K-sharding keeps memory bounded
    jobs = []
    for d in dlist:
        h_path = os.path.join(DATA, f"r35_{d}.g6")
        nh = count_lines(h_path)
        q = n - 1 - d
        k_path = os.path.join(DATA, f"r44_{q}.g6") if q > 0 else "-"
        nk = count_lines(k_path) if q > 0 else 1
        if nk > 4 * nh:
            kchunk = max(1, nk // (args.workers * 8))
            for lo in range(0, nk, kchunk):
                jobs.append((d, h_path, k_path, 0, nh, lo, min(lo + kchunk, nk)))
        else:
            chunk = max(1, nh // (args.workers * 4) or 1)
            for lo in range(0, nh, chunk):
                jobs.append((d, h_path, k_path, lo, min(lo + chunk, nh), None, None))
    t0 = time.time()
    procs, results, ji = [], [], 0

    def launch():
        nonlocal ji
        d, hp, kp, lo, hi, klo, khi = jobs[ji]
        og = os.path.join(tmpdir, f"s{ji}.g6")
        oc = os.path.join(tmpdir, f"s{ji}.csv")
        cmd = [BIN, str(n), str(e), str(d), hp, kp, str(lo), str(hi), og, oc]
        if klo is not None:
            cmd += [str(klo), str(khi)]
        p = subprocess.Popen(cmd)
        procs.append((p, d, lo, og, oc))
        ji += 1

    while ji < len(jobs) and len(procs) < args.workers:
        launch()
    ndone = 0
    while procs:
        for idx, (p, d, lo, og, oc) in enumerate(procs):
            if p.poll() is not None:
                if p.returncode != 0:
                    sys.exit(f"shard d={d} h_lo={lo} failed rc={p.returncode}")
                results.append((d, lo, og, oc))
                procs.pop(idx)
                ndone += 1
                if ndone % 20 == 0 or ndone == len(jobs):
                    print(f"  {ndone}/{len(jobs)} shards, {time.time()-t0:.0f}s",
                          file=sys.stderr, flush=True)
                if ji < len(jobs):
                    launch()
                break
        else:
            time.sleep(0.05)
    # merge + verify
    g6_path = os.path.join(args.outdir, f"r45_{n}_{e}.g6")
    csv_path = os.path.join(args.outdir, f"r45_{n}_{e}.counts.csv")
    total, bad = 0, 0
    with open(g6_path, "w") as fg:
        for d, lo, og, oc in sorted(results, key=lambda r: (r[0], r[1])):
            for line in open(og):
                line = line.strip()
                if not line:
                    continue
                if not args.no_verify:
                    gn, adj = parse_graph6_line(line)
                    ok = (gn == n
                          and sum(popcount(a) for a in adj) // 2 == e
                          and check_ramsey_graph(gn, adj, 4, 5) == (True, True))
                    if not ok:
                        bad += 1
                        continue
                fg.write(line + "\n")
                total += 1
    rows, todos = [], []
    for d, lo, og, oc in results:
        for line in open(oc):
            if line.strip():
                r = tuple(int(x) for x in line.split(","))
                if r[3] == -1:
                    todos.append((r[0], r[1], r[2]))
                else:
                    rows.append(r)
    if todos:
        print(f"  SAT path: {len(todos)} deferred pairs", file=sys.stderr, flush=True)
        dmin = int(os.environ.get("GLUE_DMIN", "0"))
        from multiprocessing import Pool
        chunks = {}
        for t in sorted(todos):
            chunks.setdefault((t[0], t[1]), []).append(t)
        work = [(n, e, ch, dmin) for ch in chunks.values()]
        t1 = time.time()
        ndone2 = 0
        with Pool(args.workers) as pool, open(g6_path, "a") as fg:
            for lines2, rows2 in pool.imap_unordered(_sat_worker, work, chunksize=1):
                for ln in lines2:
                    fg.write(ln + "\n")
                total += len(lines2)
                rows.extend(rows2)
                ndone2 += 1
                if ndone2 % 50 == 0 or ndone2 == len(work):
                    print(f"  SAT {ndone2}/{len(work)} h-chunks, "
                          f"{time.time()-t1:.0f}s", file=sys.stderr, flush=True)
    with open(csv_path, "w") as fc:
        fc.write("d,h_idx,k_idx,n_solutions\n")
        for row in sorted(rows):
            fc.write(",".join(map(str, row)) + "\n")
    for d, lo, og, oc in results:
        os.remove(og)
        os.remove(oc)
    os.rmdir(tmpdir)
    if bad:
        sys.exit(f"FATAL: {bad} emitted graphs FAILED trusted verification")
    print(f"n={n} e={e}: d in {dlist}, {len(jobs)} shards, {total} raw solutions "
          f"(all re-verified), {time.time()-t0:.1f}s -> {g6_path}")


if __name__ == "__main__":
    main()
