/* Gate-B Revision 1 cost search.
 *
 * State/batch semantics: exact port of na_compiler.exhaustive.legal_batches.
 * Cost: for each legal batch B,
 *   move_batch_duration_us(max Euclidean moved distance, 4um site spacing)
 *   + 2*15us*|B| (one load + one store per moved atom).
 * A* heuristic (pre_statement Revision 1 cost-search erratum):
 *   30us * #misplaced + max_a t(straight-line remaining distance).
 * It is admissible and consistent (see pre_statement.md); therefore the first
 * finalized target is the global optimum. Single-threaded; campaign runs use
 * nice -n 10.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int G, N, NS;
static uint8_t goal_site[8];

/* hash table */
static uint64_t cap, hmask;
static uint64_t *keys, *parents;
static double *gcost;
static uint8_t *steps, *used, *closed;
static uint32_t *hpos;
static uint64_t discovered = 0;

/* indexed binary min heap of hash-table slots */
static uint32_t *heap;
static uint32_t hn = 0, hcap = 0;
#define NO_POS UINT32_MAX

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static inline uint64_t hh(uint64_t k) {
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;
    return k;
}

static inline uint32_t slot_for(uint64_t k, int *is_new) {
    uint64_t i = hh(k) & hmask;
    while (used[i]) {
        if (keys[i] == k) {
            *is_new = 0;
            return (uint32_t)i;
        }
        i = (i + 1) & hmask;
    }
    used[i] = 1;
    keys[i] = k;
    gcost[i] = INFINITY;
    hpos[i] = NO_POS;
    discovered++;
    *is_new = 1;
    return (uint32_t)i;
}

static inline double move_duration(double d_um) {
    const double TMAX = 200.0;
    const double DMAX = 110.0;
    const double JERK = 32.0 * DMAX / (TMAX * TMAX * TMAX);
    if (d_um <= DMAX)
        return 2.0 * cbrt(4.0 * d_um / JERK);
    return TMAX + (d_um - DMAX) / 1.1;
}

static inline double heuristic(uint64_t key) {
    int misplaced = 0;
    int max_d2 = 0;
    for (int a = 0; a < N; a++) {
        int s = (int)((key >> (8 * a)) & 0xff);
        int t = goal_site[a];
        if (s == t)
            continue;
        misplaced++;
        int dx = s / G - t / G;
        int dy = s % G - t % G;
        int d2 = dx * dx + dy * dy;
        if (d2 > max_d2)
            max_d2 = d2;
    }
    return 30.0 * misplaced + move_duration(4.0 * sqrt((double)max_d2));
}

static inline double fscore(uint32_t slot) {
    return gcost[slot] + heuristic(keys[slot]);
}

static inline int heap_less(uint32_t a, uint32_t b) {
    double fa = fscore(a), fb = fscore(b);
    if (fa != fb)
        return fa < fb;
    if (gcost[a] != gcost[b])
        return gcost[a] < gcost[b];
    return keys[a] < keys[b];
}

static void heap_swap(uint32_t i, uint32_t j) {
    uint32_t a = heap[i], b = heap[j];
    heap[i] = b;
    heap[j] = a;
    hpos[a] = j;
    hpos[b] = i;
}

static void heap_sift_up(uint32_t i) {
    while (i) {
        uint32_t p = (i - 1) >> 1;
        if (!heap_less(heap[i], heap[p]))
            break;
        heap_swap(i, p);
        i = p;
    }
}

static void heap_sift_down(uint32_t i) {
    for (;;) {
        uint32_t l = i * 2 + 1;
        if (l >= hn)
            break;
        uint32_t r = l + 1;
        uint32_t m = l;
        if (r < hn && heap_less(heap[r], heap[l]))
            m = r;
        if (!heap_less(heap[m], heap[i]))
            break;
        heap_swap(i, m);
        i = m;
    }
}

static void heap_push_or_decrease(uint32_t slot) {
    if (hpos[slot] == NO_POS) {
        if (hn == hcap) {
            hcap = hcap ? hcap * 2 : (1u << 18);
            heap = realloc(heap, (size_t)hcap * sizeof(uint32_t));
            if (!heap) {
                fprintf(stderr, "OOM heap\n");
                exit(3);
            }
        }
        uint32_t i = hn++;
        heap[i] = slot;
        hpos[slot] = i;
        heap_sift_up(i);
    } else {
        heap_sift_up(hpos[slot]);
    }
}

static uint32_t heap_pop(void) {
    uint32_t out = heap[0];
    hpos[out] = NO_POS;
    hn--;
    if (hn) {
        heap[0] = heap[hn];
        hpos[heap[0]] = 0;
        heap_sift_down(0);
    }
    return out;
}

/* exact legal-batch generator globals for one expanded parent */
static int xs[8], ys[8];
static int ca[160], ct[160], cxt[160], cyt[160], nc;
static uint8_t cura[16], curt[16];
static uint64_t parent_key;
static uint32_t parent_slot;
static uint64_t generated = 0;

static void emit_batch(int depth) {
    generated++;
    uint64_t child = parent_key;
    int max_d2 = 0;
    for (int i = 0; i < depth; i++) {
        int a = cura[i], t = curt[i];
        int s = (int)((parent_key >> (8 * a)) & 0xff);
        int dx = s / G - t / G;
        int dy = s % G - t % G;
        int d2 = dx * dx + dy * dy;
        if (d2 > max_d2)
            max_d2 = d2;
        child &= ~(0xffULL << (8 * a));
        child |= (uint64_t)t << (8 * a);
    }
    double edge = move_duration(4.0 * sqrt((double)max_d2)) + 30.0 * depth;
    double ng = gcost[parent_slot] + edge;
    int is_new = 0;
    uint32_t cs = slot_for(child, &is_new);
    if (closed[cs]) {
        if (ng + 1e-10 < gcost[cs]) {
            /* A lower path to a closed node would falsify heuristic consistency. */
            fprintf(stderr, "ERROR: inconsistent heuristic reopen %.12f < %.12f\n",
                    ng, gcost[cs]);
            exit(4);
        }
        return;
    }
    if (ng + 1e-10 < gcost[cs]) {
        gcost[cs] = ng;
        steps[cs] = (uint8_t)(steps[parent_slot] + 1);
        parents[cs] = parent_key;
        heap_push_or_decrease(cs);
    }
}

static void rec_batches(int start_idx, int depth) {
    for (int idx = start_idx; idx < nc; idx++) {
        int a = ca[idx], t = ct[idx], xt = cxt[idx], yt = cyt[idx];
        int ok = 1;
        for (int k = 0; k < depth; k++) {
            if (cura[k] == a || curt[k] == t) {
                ok = 0;
                break;
            }
            int xa2 = xs[cura[k]], ya2 = ys[cura[k]];
            int x2 = curt[k] / G, y2 = curt[k] % G;
            int xac = xs[a], yac = ys[a];
            if (((xac == xa2) != (xt == x2)) ||
                ((yac == ya2) != (yt == y2)) ||
                ((xac < xa2) != (xt < x2)) ||
                ((yac < ya2) != (yt < y2))) {
                ok = 0;
                break;
            }
        }
        if (!ok)
            continue;
        cura[depth] = (uint8_t)a;
        curt[depth] = (uint8_t)t;
        emit_batch(depth + 1);
        rec_batches(idx + 1, depth + 1);
    }
}

static void expand_slot(uint32_t slot) {
    parent_slot = slot;
    parent_key = keys[slot];
    uint64_t occ = 0;
    for (int a = 0; a < N; a++) {
        int s = (int)((parent_key >> (8 * a)) & 0xff);
        xs[a] = s / G;
        ys[a] = s % G;
        occ |= 1ULL << s;
    }
    nc = 0;
    for (int a = 0; a < N; a++) {
        for (int t = 0; t < NS; t++) {
            if ((occ >> t) & 1)
                continue;
            ca[nc] = a;
            ct[nc] = t;
            cxt[nc] = t / G;
            cyt[nc] = t % G;
            nc++;
        }
    }
    rec_batches(0, 0);
}

static void print_state(uint64_t k) {
    printf("[");
    for (int a = 0; a < N; a++)
        printf(a ? ",%d" : "%d", (int)((k >> (8 * a)) & 0xff));
    printf("]");
}

int main(int argc, char **argv) {
    if (argc != 4 + atoi(argv[1])) {
        fprintf(stderr, "usage: cost_astar N G capacity_log2 goal0 ... goalN-1\n");
        return 2;
    }
    N = atoi(argv[1]);
    G = atoi(argv[2]);
    NS = G * G;
    int clog = atoi(argv[3]);
    cap = 1ULL << clog;
    hmask = cap - 1;
    for (int a = 0; a < N; a++)
        goal_site[a] = (uint8_t)atoi(argv[4 + a]);

    keys = calloc(cap, sizeof(uint64_t));
    parents = calloc(cap, sizeof(uint64_t));
    gcost = malloc(cap * sizeof(double));
    steps = calloc(cap, sizeof(uint8_t));
    used = calloc(cap, sizeof(uint8_t));
    closed = calloc(cap, sizeof(uint8_t));
    hpos = malloc(cap * sizeof(uint32_t));
    if (!keys || !parents || !gcost || !steps || !used || !closed || !hpos) {
        fprintf(stderr, "OOM tables\n");
        return 3;
    }
    memset(hpos, 0xff, cap * sizeof(uint32_t));

    uint64_t start = 0, goal = 0;
    for (int a = 0; a < N; a++) {
        start |= (uint64_t)a << (8 * a);
        goal |= (uint64_t)goal_site[a] << (8 * a);
    }
    int is_new = 0;
    uint32_t ss = slot_for(start, &is_new);
    gcost[ss] = 0.0;
    parents[ss] = start;
    heap_push_or_decrease(ss);

    double t0 = now_s();
    uint64_t expanded = 0;
    uint32_t goal_slot = NO_POS;
    while (hn) {
        uint32_t s = heap_pop();
        if (closed[s])
            continue;
        closed[s] = 1;
        expanded++;
        if (keys[s] == goal) {
            goal_slot = s;
            break;
        }
        expand_slot(s);
        if ((expanded % 100000) == 0) {
            fprintf(stderr,
                    "{\"expanded\":%llu,\"discovered\":%llu,\"heap\":%u,\"cum_s\":%.3f}\n",
                    (unsigned long long)expanded,
                    (unsigned long long)discovered, hn, now_s() - t0);
        }
    }
    if (goal_slot == NO_POS) {
        fprintf(stderr, "goal unreachable\n");
        return 5;
    }

    /* Reconstruct path of states. Positive edges imply no cycle in optimum. */
    uint64_t path[256];
    int plen = 0;
    uint64_t k = goal;
    while (1) {
        if (plen >= 256) {
            fprintf(stderr, "path too long\n");
            return 6;
        }
        path[plen++] = k;
        if (k == start)
            break;
        int dummy = 0;
        uint32_t sl = slot_for(k, &dummy);
        k = parents[sl];
    }

    double wall = now_s() - t0;
    printf("{\"target\":");
    print_state(goal);
    printf(",\"opt_cost_us\":%.12f,\"opt_cost_batches\":%u,"
           "\"expanded\":%llu,\"generated\":%llu,"
           "\"discovered\":%llu,\"wall_s\":%.6f,\"witness_states\":[",
           gcost[goal_slot], (unsigned)steps[goal_slot],
           (unsigned long long)expanded,
           (unsigned long long)generated,
           (unsigned long long)discovered, wall);
    for (int i = plen - 1; i >= 0; i--) {
        if (i != plen - 1)
            printf(",");
        print_state(path[i]);
    }
    printf("]}\n");
    return 0;
}
