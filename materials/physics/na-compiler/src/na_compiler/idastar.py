"""Gate-B Revision 1: complete BFS rows on 4x4 beyond the 30k state cap,
via (a) bidirectional BFS for single-target optimal distances and (b) a
compile-on-demand bitstate C enumerator for full per-layer frontier data.

Revision 1 (pre_statement.md, appended BEFORE any run) pinned:
- h(s) = ceil( max_a d_manhattan(site(a), goal_a) / D_max ), where D_max is
  the max Manhattan grid-step displacement any single atom can achieve in one
  legal batch on the fixed grid ((S-1)+(S-1) for S x S);
- admissibility (R1.3): one batch displaces each atom once, so the max
  remaining Manhattan distance drops by at most D_max per batch;
  #batches >= h(s). Consistency: any batch changes h by at most 1.
- Complete row (R1.4) = explicit witness of k legal batches (verified against
  the FROZEN legality code, i.e. `na_compiler.conflict.compatible` plus the
  §4 free-target rule via `exhaustive.legal_batches` membership) AND proof
  that no k-1-batch sequence reaches the target (exhaustion, pruned only by
  the consistent heuristic).

This module reuses frozen legality semantics; it introduces NO new legality
relation. The fast generator (`legal_batches_fast`) and the C enumerator
reproduce `exhaustive.legal_batches` exactly; both were equality-tested
against the frozen implementation on thousands of random states (n=3..6,
grids 3x3/4x4) before the campaign runs, and the C enumerator's per-layer
frontier sizes were required to match the Python BFS on (4,4x4) and (5,4x4)
before any larger run was archived.
"""

from __future__ import annotations

import subprocess
from math import ceil
from pathlib import Path

from .conflict import Move, compatible
from .exhaustive import apply_batch, legal_batches

# --------------------------------------------------------------------------- #
# Revision 1 pinned heuristic
# --------------------------------------------------------------------------- #


def dmax_grid(grid_side: int) -> int:
    """Max Manhattan grid-step travel of one atom in one legal batch on SxS."""
    return 2 * (grid_side - 1)


def manhattan_heuristic(
    state: tuple[int, ...], goal: tuple[int, ...], grid_side: int
) -> int:
    """h(state) = ceil(max_a Manhattan(site(a), goal_a) / D_max).

    Admissible + consistent for the batch metric (R1.2/R1.3)."""
    dmax = dmax_grid(grid_side)
    worst = 0
    for a, s in enumerate(state):
        x0, y0 = divmod(s, grid_side)
        x1, y1 = divmod(goal[a], grid_side)
        d = abs(x0 - x1) + abs(y0 - y1)
        if d > worst:
            worst = d
    return ceil(worst / dmax) if worst else 0


# --------------------------------------------------------------------------- #
# Fast exact replica of the frozen batch generator (equality-tested)
# --------------------------------------------------------------------------- #
def legal_batches_fast(state: tuple[int, ...], grid_side: int):
    """Bit-for-bit replica of `exhaustive.legal_batches` output (same set of
    batches, same `(atom, target)` pair encoding)."""
    ns = grid_side * grid_side
    occ = set(state)
    xs = [s // grid_side for s in state]
    ys = [s % grid_side for s in state]
    cands: list[tuple[int, int, int, int]] = []
    for a in range(len(state)):
        for t in range(ns):
            if t not in occ:
                cands.append((a, t, t // grid_side, t % grid_side))
    nc = len(cands)
    res: list[tuple[tuple[int, int], ...]] = []
    cur: list[tuple[int, int]] = []

    def rec(start_idx: int) -> None:
        res.append(tuple(cur))
        for idx in range(start_idx, nc):
            a, t, xt, yt = cands[idx]
            ok = True
            for a2, t2 in cur:
                if a == a2 or t == t2:
                    ok = False
                    break
                xa2 = xs[a2]
                ya2 = ys[a2]
                x2 = t2 // grid_side
                y2 = t2 % grid_side
                xac = xs[a]
                yac = ys[a]
                if ((xac == xa2) != (xt == x2)) or ((yac == ya2) != (yt == y2)):
                    ok = False
                    break
                if ((xac < xa2) != (xt < x2)) or ((yac < ya2) != (yt < y2)):
                    ok = False
                    break
            if ok:
                cur.append((a, t))
                rec(idx + 1)
                cur.pop()

    rec(0)
    return res


# --------------------------------------------------------------------------- #
# Witness verification against FROZEN legality machinery
# --------------------------------------------------------------------------- #
def verify_witness(
    start: tuple[int, ...],
    target: tuple[int, ...],
    path: list[tuple[tuple[int, int], ...]],
    grid_side: int,
) -> None:
    """Raise unless every batch is legal under frozen §4 semantics and the
    path maps start -> target in len(path) batches."""
    state = tuple(start)
    for k, batch in enumerate(path):
        if batch not in legal_batches(state, grid_side):
            raise ValueError(f"batch {k} not legal at {state}")
        state = apply_batch(state, batch)
    if state != tuple(target):
        raise ValueError(f"witness ends at {state}, not target {tuple(target)}")


# --------------------------------------------------------------------------- #
# IDA* with the pinned heuristic — certifies optimal distance to ONE target
# --------------------------------------------------------------------------- #
def ida_star_distance(
    n_atoms: int,
    grid_side: int,
    target: tuple[int, ...],
    max_bound: int = 12,
    node_cap: int | None = None,
) -> dict:
    """Optimal #batches identity -> `target` via IDA* with h = manhattan.

    Returns {distance, nodes, witnessed, complete}: `complete=True` iff the
    search exhausted all bounds below the returned distance (optimality
    proof) and a witness path exists (verified against frozen legality).
    """
    from .exhaustive import legal_batches as lb

    start = tuple(range(n_atoms))
    target = tuple(target)
    if start == target:
        return {"distance": 0, "nodes": 1, "witnessed": True, "complete": True,
                "witness": []}
    nodes = 0
    h0 = manhattan_heuristic(start, target, grid_side)

    def dfs(state, g, bound, path, path_states):
        nonlocal nodes
        nodes += 1
        if node_cap is not None and nodes > node_cap:
            raise RuntimeError("node cap exceeded")
        f = g + manhattan_heuristic(state, target, grid_side)
        if f > bound:
            return f
        if state == target:
            return -1  # found
        best = 1 << 30
        for b in lb(state, grid_side):
            if not b:
                continue
            s2 = list(state)
            for a, t in b:
                s2[a] = t
            s2 = tuple(s2)
            if s2 in path_states:
                continue
            path.append(b)
            path_states.add(s2)
            r = dfs(s2, g + 1, bound, path, path_states)
            if r == -1:
                return -1
            path_states.remove(s2)
            path.pop()
            if r < best:
                best = r
        return best

    witness = None
    bound = h0
    complete = False
    while bound <= max_bound:
        path: list = []
        path_states = {start}
        try:
            r = dfs(start, 0, bound, path, path_states)
        except RuntimeError:
            return {"distance": None, "nodes": nodes, "witnessed": False,
                    "complete": False, "witness": None}
        if r == -1:
            witness = list(path)
            # exhaustion of all shallower bounds already done by construction
            # (previous iterations returned min f >= their bound without find)
            complete = True
            break
        if r >= (1 << 30):
            # exhausted bound without finding AND no frontier: space exhausted
            complete = True
            break
        bound = r
    else:
        complete = False
    witnessed = witness is not None
    if witnessed:
        verify_witness(start, target, witness, grid_side)
    return {"distance": bound if witnessed else None, "nodes": nodes,
            "witnessed": witnessed, "complete": complete and witnessed,
            "witness": witness}


# --------------------------------------------------------------------------- #
# Bidirectional BFS single-target distance
# --------------------------------------------------------------------------- #
def bidir_distance(
    n_atoms: int,
    grid_side: int,
    target: tuple[int, ...],
    max_depth: int = 16,
) -> tuple[int | None, dict]:
    """Min #batches identity -> target, meet-in-the-middle.

    The move graph is undirected: the four compatibility conditions are
    iff-statements symmetric under reversing every move, and the free-target
    condition mirrors to a free-source condition (both property-tested
    against the frozen machinery on thousands of edges before use). Hence
    reverse expansion from `target` uses the SAME frozen `legal_batches`.
    Stop criterion: both sides expanded to depths lf+lb >= best candidate.
    """
    ident = tuple(range(n_atoms))
    target = tuple(target)
    if ident == target:
        return 0, {"lf": 0, "lb": 0, "visits": 1}
    df: dict[tuple, int] = {ident: 0}
    db: dict[tuple, int] = {target: 0}
    ff, fb = [ident], [target]
    lf = lbv = 0
    visits = 0

    def expand(frontier, dist, dcur):
        nonlocal visits
        new, seen = [], set()
        for st in frontier:
            visits += 1
            for b in legal_batches_fast(st, grid_side):
                s2 = list(st)
                for a, t in b:
                    s2[a] = t
                s2 = tuple(s2)
                if s2 not in dist and s2 not in seen:
                    seen.add(s2)
                    dist[s2] = dcur + 1
                    new.append(s2)
        return new

    while ff and fb:
        if len(ff) <= len(fb):
            ff = expand(ff, df, lf)
            lf += 1
        else:
            fb = expand(fb, db, lbv)
            lbv += 1
        meet = set(df) & set(db)
        if meet:
            cand = min(df[m] + db[m] for m in meet)
            if lf + lbv >= cand:
                stats = {"lf": lf, "lb": lbv, "visits": visits,
                         "frontier_f": len(ff), "frontier_b": len(fb)}
                return cand, stats
        if lf + lbv > max_depth:
            return None, {"lf": lf, "lb": lbv, "visits": visits}
    return None, {"lf": lf, "lb": lbv, "visits": visits}


# --------------------------------------------------------------------------- #
# Compile-on-demand bitstate C enumerator for complete frontier data
# --------------------------------------------------------------------------- #
C_BFS_SOURCE = r"""
/* Gate-B Revision 1 bitstate enumerator.
   Semantics = na_compiler.exhaustive.legal_batches (verified: equality with
   the frozen Python generator on random states, then per-layer frontier
   match vs Python BFS at (4,4x4) and (5,4x4) before larger runs).
   Output: JSON {"frontiers":[[...],[...],...],"visited":N,"complete":1,
   "layer_wall_s":[...], "diameter":D} on stdout; progress on stderr. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

static int G, N, NS;
static uint64_t cap, hmask;
static uint64_t *keys;
static uint8_t *dists;
static uint8_t *used;

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static inline uint64_t hh(uint64_t k) {
    k ^= k >> 33; k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33; return k;
}

static inline int ht_insert(uint64_t k, uint8_t d) {
    uint64_t i = hh(k) & hmask;
    while (used[i]) {
        if (keys[i] == k) return 0;
        i = (i + 1) & hmask;
    }
    used[i] = 1; keys[i] = k; dists[i] = d;
    return 1;
}

static inline uint8_t ht_get(uint64_t k) {
    uint64_t i = hh(k) & hmask;
    while (used[i]) {
        if (keys[i] == k) return dists[i];
        i = (i + 1) & hmask;
    }
    return 0xff;
}

static uint64_t *front;
static uint64_t fn, fcap;
static void push(uint64_t k) {
    if (fn == fcap) {
        fcap = fcap ? fcap * 2 : (1u << 16);
        front = realloc(front, fcap * sizeof(uint64_t));
        if (!front) { fprintf(stderr, "OOM frontier\n"); exit(3); }
    }
    front[fn++] = k;
}

static int xs[8], ys[8];
static uint8_t cur_dist;

/* one-parent expansion: iterative exact replica of the frozen recursion
   extend(cur, start_idx) over candidate moves (atom-major by target). */
static void gen(uint64_t pk) {
    uint64_t occ_set = 0;
    for (int a = 0; a < N; a++) {
        int s = (int)((pk >> (8 * a)) & 0xff);
        xs[a] = s / G;
        ys[a] = s % G;
        occ_set |= 1ULL << s;
    }
    /* candidate list */
    int ca[160], ct[160], cxt[160], cyt[160];
    int nc = 0;
    for (int a = 0; a < N; a++) {
        for (int t = 0; t < NS; t++) {
            if (!((occ_set >> t) & 1)) {
                ca[nc] = a; ct[nc] = t;
                cxt[nc] = t / G; cyt[nc] = t % G;
                nc++;
            }
        }
    }
    uint8_t cura[16], curt[16];
    int stack[17];
    int frame = 0;
    stack[0] = 0;
    for (;;) {
        if (stack[frame] >= nc) {
            if (frame == 0) break;
            frame--;
            stack[frame]++;
            continue;
        }
        int ci = stack[frame];
        int a = ca[ci], t = ct[ci], xt = cxt[ci], yt = cyt[ci];
        int ok = 1;
        for (int k = 0; k < frame; k++) {
            if (cura[k] == a || curt[k] == t) { ok = 0; break; }
            int xa2 = xs[cura[k]], ya2 = ys[cura[k]];
            int x2 = curt[k] / G, y2 = curt[k] % G;
            int xac = xs[a], yac = ys[a];
            if (((xac == xa2) != (xt == x2)) || ((yac == ya2) != (yt == y2))) { ok = 0; break; }
            if (((xac <  xa2) != (xt <  x2)) || ((yac <  ya2) != (yt <  y2))) { ok = 0; break; }
        }
        if (ok) {
            cura[frame] = (uint8_t)a;
            curt[frame] = (uint8_t)t;
            uint64_t child = pk;
            for (int i = 0; i <= frame; i++) {
                child &= ~(0xffULL << (8 * cura[i]));
                child |= (uint64_t)curt[i] << (8 * cura[i]);
            }
            if (child != pk && ht_insert(child, cur_dist + 1)) push(child);
            frame++;
            stack[frame] = ci + 1;
        } else {
            stack[frame]++;
        }
    }
}

static int key_lsb_lex(uint64_t x, uint64_t y) {
    /* compare states as (atom0 site, atom1 site, ...): atom0 = low byte */
    for (int a = 0; a < N; a++) {
        int bx = (int)((x >> (8 * a)) & 0xff);
        int by = (int)((y >> (8 * a)) & 0xff);
        if (bx != by) return bx < by ? -1 : 1;
    }
    return 0;
}
int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: cbfs N G capacity_log2\n");
        return 2;
    }
    N = atoi(argv[1]);
    G = atoi(argv[2]);
    NS = G * G;
    cap = 1ULL << atoi(argv[3]);
    hmask = cap - 1;
    keys = calloc(cap, sizeof(uint64_t));
    dists = calloc(cap, sizeof(uint8_t));
    used = calloc(cap, sizeof(uint8_t));
    if (!keys || !dists || !used) { fprintf(stderr, "OOM table\n"); return 3; }

    uint64_t start = 0;
    for (int a = 0; a < N; a++) start |= (uint64_t)a << (8 * a);
    ht_insert(start, 0);
    push(start);

    double t0 = now_s();
    uint64_t visited = 1;
    uint8_t d = 0;
    int complete = 1;
    printf("{\"frontiers\":[");
    int first_layer = 1;
    double lw[64];
    int ln = 0;
    while (1) {
        uint64_t cn = fn;
        if (!first_layer) printf(",");
        first_layer = 0;
        printf("%llu", (unsigned long long)cn);
        if (cn == 0) break;               /* trailing empty layer = BFS stop */
        uint64_t *curf = front;
        front = NULL; fn = 0; fcap = 0;
        for (uint64_t i = 0; i < cn; i++) {
            cur_dist = ht_get(curf[i]);
            gen(curf[i]);
        }
        free(curf);
        visited += fn;
        double t1 = now_s();
        lw[ln++] = t1 - t0;
        fprintf(stderr, "{\"layer\":%d,\"frontier\":%llu,\"cum_s\":%.3f}\n",
                d + 1, (unsigned long long)fn, t1 - t0);
        d++;
    }
    /* hardest states: scan table for keys at dist == diameter, keep the K
       smallest under (atom0 site, atom1 site, ...) lexicographic order. */
    {
        const int K = 10;
        uint64_t best[K]; int nbest = 0;
        int diam = (int)d - 1;
        for (uint64_t i = 0; i < cap; i++) {
            if (!used[i] || dists[i] != diam) continue;
            if (nbest == K && key_lsb_lex(keys[i], best[K-1]) >= 0) continue;
            int pos = nbest < K ? nbest : K - 1;
            while (pos > 0 && key_lsb_lex(keys[i], best[pos - 1]) < 0) {
                best[pos] = best[pos - 1];
                pos--;
            }
            best[pos] = keys[i];
            if (nbest < K) nbest++;
        }
        printf("],\"hardest\":[");
        for (int i = 0; i < nbest; i++) {
            printf(i ? "," : "");
            printf("[");
            for (int a = 0; a < N; a++)
                printf(a ? ",%d" : "%d", (int)((best[i] >> (8 * a)) & 0xff));
            printf("]");
        }
        long long ndiam = 0;
        for (uint64_t i = 0; i < cap; i++)
            if (used[i] && dists[i] == diam) ndiam++;
        printf("],\"n_diameter_states\":%lld", ndiam);
    }
    printf(",\"visited\":%llu,\"complete\":1,\"layer_wall_s\":[",
           (unsigned long long)visited);
    for (int i = 0; i < ln; i++)
        printf(i ? ",%.3f" : "%.3f", lw[i]);
    printf("],\"diameter\":%d}\n", (int)d - 1);  /* d = #expansions incl. empty layer */
    return 0;
}
"""


def c_bfs_complete(
    n_atoms: int,
    grid_side: int,
    workdir: str | Path,
    capacity_log2: int,
    timeout_s: int = 3600,
) -> dict:
    """Compile (once) the bitstate enumerator and run the COMPLETE BFS from
    identity. Returns parsed JSON result + compile/run provenance."""
    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    cpath = wd / "cbfs.c"
    binpath = wd / "cbfs"
    cpath.write_text(C_BFS_SOURCE)
    if not binpath.exists():
        cc = subprocess.run(
            ["cc", "-O2", "-std=c11", "-o", str(binpath), str(cpath)],
            capture_output=True, text=True,
        )
        if cc.returncode != 0:
            raise RuntimeError(f"cc failed: {cc.stderr[:2000]}")
    import json
    import time

    t0 = time.perf_counter()
    run = subprocess.run(
        ["nice", "-n", "10", str(binpath), str(n_atoms), str(grid_side),
         str(capacity_log2)],
        capture_output=True, text=True, timeout=timeout_s,
    )
    wall = time.perf_counter() - t0
    if run.returncode != 0:
        raise RuntimeError(f"cbfs failed rc={run.returncode}: {run.stderr[-2000:]}")
    out = json.loads(run.stdout)
    out["wall_s"] = wall
    out["progress_stderr"] = run.stderr
    return out


def cost_astar_complete(
    n_atoms: int,
    grid_side: int,
    target: tuple[int, ...],
    workdir: str | Path,
    capacity_log2: int,
    timeout_s: int = 3600,
) -> dict:
    """Certify global minimum additive transport cost to one target with A*.

    The compiled search uses the frozen batch graph and the pinned edge cost
    and consistent heuristic from the Revision 1 cost-search erratum.
    """
    if len(target) != n_atoms:
        raise ValueError("target length must equal n_atoms")
    if len(set(target)) != n_atoms:
        raise ValueError("target sites must be distinct")
    if any(site < 0 or site >= grid_side * grid_side for site in target):
        raise ValueError("target site outside grid")

    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    source = Path(__file__).with_name("cost_astar.c")
    binpath = wd / "cost_astar"
    if not binpath.exists():
        cc = subprocess.run(
            ["cc", "-O2", "-std=c11", str(source), "-o", str(binpath), "-lm"],
            capture_output=True,
            text=True,
        )
        if cc.returncode != 0:
            raise RuntimeError(f"cc failed: {cc.stderr[:2000]}")

    import json
    import time

    t0 = time.perf_counter()
    run = subprocess.run(
        [
            "nice",
            "-n",
            "10",
            str(binpath),
            str(n_atoms),
            str(grid_side),
            str(capacity_log2),
            *(str(site) for site in target),
        ],
        capture_output=True,
        text=True,
        timeout=timeout_s,
    )
    wall = time.perf_counter() - t0
    if run.returncode != 0:
        raise RuntimeError(
            f"cost A* failed rc={run.returncode}: {run.stderr[-2000:]}"
        )
    out = json.loads(run.stdout)
    out["driver_wall_s"] = wall
    out["progress_stderr"] = run.stderr
    return out


__all__ = [
    "dmax_grid",
    "manhattan_heuristic",
    "legal_batches_fast",
    "verify_witness",
    "ida_star_distance",
    "bidir_distance",
    "c_bfs_complete",
    "cost_astar_complete",
    "C_BFS_SOURCE",
]
