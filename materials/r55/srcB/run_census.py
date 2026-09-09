#!/usr/bin/env python3
"""Driver for implementation B of the R(4,5,n,e) census (frozen spec).

For each ladder stratum: chunk (d, K-range) tasks across cores, merge raw
.g6 + counts, run the independent full verifier, dedup with nauty shortg -k,
and compare canonical sets against the published catalog via labelg.
"""
import csv
import math
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRCB = os.path.join(BASE, "srcB")
OUTB = os.path.join(BASE, "outB")
DATA = os.path.join(BASE, "data")
GLUE = os.path.join(SRCB, "glue2")
VERIFY = os.path.join(SRCB, "verify")
SHORTG = "/opt/homebrew/bin/shortg"
LABELG = "/opt/homebrew/bin/labelg"
NCPU = os.cpu_count() or 8

LADDER = [
    (12, 48, 1),
    (13, 53, 2),
    (13, 52, 10),
    (14, 60, 1),
    (15, 66, 1),
    (16, 72, 5),
    (16, 71, 138),
    (17, 79, 1),
    (17, 78, 86),
    (21, 107, 31),
]


def catalog_size(q):
    if q == 0:
        return 1
    with open(os.path.join(DATA, f"r44_{q}.g6")) as f:
        return sum(1 for line in f if line.strip())


def run_chunk(args):
    n, e, d, k0, k1, g6p, csvp = args
    r = subprocess.run(
        [GLUE, str(n), str(e), str(d), str(k0), str(k1), g6p, csvp, DATA],
        capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"glue failed d={d} k[{k0},{k1}): {r.stderr}")
    return g6p, csvp


def run_stratum(n, e, expected):
    t0 = time.time()
    dmin = max((2 * e + n - 1) // n, n - 1 - 17, 1)
    dmax = min(13, n - 1)
    tmpdir = os.path.join(OUTB, f"tmp_{n}_{e}")
    os.makedirs(tmpdir, exist_ok=True)
    tasks = []
    for d in range(dmin, dmax + 1):
        q = n - 1 - d
        if q < 0 or q > 17:
            continue
        nk = catalog_size(q)
        nchunks = min(nk, NCPU * 8) if nk > 1 else 1
        step = math.ceil(nk / nchunks)
        for ci, k0 in enumerate(range(0, nk, step)):
            k1 = min(k0 + step, nk)
            g6p = os.path.join(tmpdir, f"d{d}_c{ci}.g6")
            csvp = os.path.join(tmpdir, f"d{d}_c{ci}.csv")
            tasks.append((n, e, d, k0, k1, g6p, csvp))

    with ProcessPoolExecutor(max_workers=NCPU) as ex:
        results = list(ex.map(run_chunk, tasks))

    raw_path = os.path.join(OUTB, f"r45_{n}_{e}.g6")
    counts_path = os.path.join(OUTB, f"r45_{n}_{e}.counts.csv")
    rows = []
    raw = 0
    with open(raw_path, "w") as out:
        for g6p, csvp in results:
            with open(g6p) as f:
                for line in f:
                    if line.strip():
                        out.write(line)
                        raw += 1
            with open(csvp) as f:
                for row in csv.reader(f):
                    if row:
                        rows.append((int(row[0]), int(row[1]), int(row[2]), int(row[3])))
    rows.sort()
    with open(counts_path, "w") as f:
        f.write("d,h_idx,k_idx,n_solutions\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]},{r[2]},{r[3]}\n")
    assert sum(r[3] for r in rows) == raw

    # independent full verification of every output graph
    v = subprocess.run([VERIFY, str(n), str(e), raw_path],
                       capture_output=True, text=True)
    verify_line = v.stdout.strip()
    if v.returncode != 0:
        raise RuntimeError(f"VERIFIER REJECTED SOLUTIONS ({n},{e}): "
                           f"{verify_line}\n{v.stderr[:2000]}")

    # dedup count via shortg -u (spec), canonical sets via labelg on both sides
    sg = subprocess.run([SHORTG, "-u", raw_path], capture_output=True, text=True)
    shortg_iso = None
    for line in sg.stderr.splitlines():
        if "graphs produced" in line:
            shortg_iso = int(line.split()[1])
    canon_path = os.path.join(tmpdir, "canon.g6")
    subprocess.run([LABELG, "-q", raw_path, canon_path], check=True)
    with open(canon_path) as f:
        mine = sorted({line.strip() for line in f if line.strip()})
    pub_src = os.path.join(DATA, "r45extreme", f"r45{n}.{e}.g6")
    pub_canon = os.path.join(tmpdir, "pub_canon.g6")
    subprocess.run([LABELG, "-q", pub_src, pub_canon], check=True)
    with open(pub_canon) as f:
        pub = sorted({line.strip() for line in f if line.strip()})
    iso = len(mine)
    if shortg_iso is not None and shortg_iso != iso:
        raise RuntimeError(f"shortg count {shortg_iso} != labelg set size {iso}")
    match = mine == pub
    secs = time.time() - t0
    return dict(n=n, e=e, raw=raw, iso=iso, expected=expected,
                match=match, secs=secs, verify=verify_line)


def main():
    os.makedirs(OUTB, exist_ok=True)
    only = None
    if len(sys.argv) == 3:
        only = (int(sys.argv[1]), int(sys.argv[2]))
    for n, e, expected in LADDER:
        if only and (n, e) != only:
            continue
        print(f"=== ({n},{e}) expected {expected} ===", flush=True)
        r = run_stratum(n, e, expected)
        print(f"  raw={r['raw']} iso={r['iso']} expected={expected} "
              f"match={r['match']} {r['secs']:.1f}s [{r['verify']}]", flush=True)
        with open(os.path.join(OUTB, "ladder_report.csv"), "a") as f:
            f.write(f"{r['n']},{r['e']},{r['raw']},{r['iso']},{r['expected']},"
                    f"{r['match']},{r['secs']:.2f},{r['verify']}\n")
        if not r["match"]:
            print("MISMATCH — stopping per spec discipline", flush=True)
            break


if __name__ == "__main__":
    main()
