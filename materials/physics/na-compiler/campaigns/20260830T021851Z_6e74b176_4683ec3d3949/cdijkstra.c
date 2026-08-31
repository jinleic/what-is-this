/* cdijkstra — global µs-optimal transport cost over the full reachable
 * NA-transport state space, exact-replica semantics.
 *
 * Provenance: batch legality + child generation are a line-for-line replica
 * of the frozen gate-B enumerator
 * campaigns/20260830T123706Z_4ed93f70_9bc5ba162924/cbfs.c gen() (its
 * equality with na_compiler.exhaustive.legal_batches was verified in that
 * campaign). Edge weight per frozen config.py:
 *   w(batch) = move_batch_duration_us(dmax) + 2*15.0*|batch|
 *   dmax = max euclidean distance in the batch (grid pitch 4.0 µm)
 *
 * Costs are ordered via fixed point (quantum 2^-20 µs ≈ 9.5e-7 µs — orders
 * below the 1e-6 comparison band). The FINAL check compares the extracted
 * witness replayed in exact Python arithmetic (verify_cost_rows.py style)
 * against the frozen claim at 1e-9 — the fixed point is only for ordering.
 *
 * Output on stdout (JSON): visited count, per-target {cost_us, reached},
 * and the witness (list of batches, each a list of [atom,site] pairs).
 * Progress on stderr. Single core, memory ~cap*(8+8+8+1) bytes.
 *
 * Usage: cdijkstra N G caplog2 t1 t2 ...  (each target: comma-separated
 * sites for atom0,atom1,..,atomN-1)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define TRANSFER_US 15.0
#define D_MAX_UM 110.0
#define T_D_MAX_US 200.0
static const double JERK = 32.0 * 110.0 / (200.0 * 200.0 * 200.0);
static const double VMAX = 110.0 / 200.0 * 2.0;
#define FPQ 1048576.0  /* 2^-20 µs */
#define INF_COST UINT64_MAX
#define NOPARENT UINT64_MAX

static int G, N, NS;
static uint64_t cap, hmask;
static uint64_t *keys;
static uint64_t *cost_fp;
static uint64_t *parent;
static uint8_t  *used;

static int ntargets;
static uint64_t target_keys[8];
static int target_hit[8];

static double now_s(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static inline double move_dur_us(double d) {
    if (d <= D_MAX_UM) return 2.0 * cbrt(4.0 * d / JERK);
    return T_D_MAX_US + (d - D_MAX_UM) / VMAX;
}

static inline uint64_t hh(uint64_t k) {
    k ^= k >> 33; k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL; k ^= k >> 33; return k;
}

static inline uint64_t slot_for(uint64_t k, int *fresh) {
    uint64_t i = hh(k) & hmask;
    while (used[i]) {
        if (keys[i] == k) { *fresh = 0; return i; }
        i = (i + 1) & hmask;
    }
    used[i] = 1; keys[i] = k; cost_fp[i] = INF_COST; parent[i] = NOPARENT;
    *fresh = 1; return i;
}
static inline uint64_t find_slot(uint64_t k) {
    uint64_t i = hh(k) & hmask;
    while (used[i]) { if (keys[i] == k) return i; i = (i + 1) & hmask; }
    return UINT64_MAX;
}

typedef struct { uint64_t c, k; } HE;
static HE *heap; static uint64_t hn, hcap;
static void hpush(uint64_t c, uint64_t k) {
    if (hn == hcap) {
        hcap = hcap ? hcap * 2 : (1u << 23);
        heap = realloc(heap, hcap * sizeof(HE));
        if (!heap) { fprintf(stderr, "OOM heap\n"); exit(3); }
    }
    uint64_t i = hn++;
    while (i > 0) { uint64_t p = (i - 1) >> 1; if (heap[p].c <= c) break; heap[i] = heap[p]; i = p; }
    heap[i].c = c; heap[i].k = k;
}
static HE hpop(void) {
    HE top = heap[0]; hn--;
    if (!hn) return top;
    HE last = heap[hn]; uint64_t i = 0;
    for (;;) {
        uint64_t l = 2*i+1, r = l+1, m = hn;
        uint64_t best = last.c;
        if (l < hn && heap[l].c < best) { m = l; best = heap[l].c; }
        if (r < hn && heap[r].c < best) { m = r; best = heap[r].c; }
        if (m == hn) break;
        heap[i] = heap[m]; i = m;
    }
    heap[i] = last; return top;
}

static int xs[8], ys[8];
static inline double site_dist_us_frame(int s1, int s2) {
    int dx = s1 / G - s2 / G, dy = s1 % G - s2 % G;
    return 4.0 * sqrt((double)(dx*dx + dy*dy));
}

static int ca[160], ct[160], nc;
static uint8_t cura[16], curt[16];
static int stack[17];

static void expand(uint64_t pk, uint64_t pcost) {
    uint64_t occ = 0;
    for (int a = 0; a < N; a++) {
        int s = (int)((pk >> (8 * a)) & 0xff);
        xs[a] = s / G; ys[a] = s % G; occ |= 1ULL << s;
    }
    nc = 0;
    for (int a = 0; a < N; a++)
        for (int t = 0; t < NS; t++)
            if (!((occ >> t) & 1)) { ca[nc] = a; ct[nc] = t; nc++; }

    int frame = 0; stack[0] = 0;
    for (;;) {
        if (stack[frame] >= nc) {
            if (frame == 0) break;
            frame--; stack[frame]++; continue;
        }
        int ci = stack[frame];
        int a = ca[ci], t = ct[ci];
        int ok = 1;
        for (int k = 0; k < frame; k++) {
            if (cura[k] == (uint8_t)a || curt[k] == (uint8_t)t) { ok = 0; break; }
            int a2 = (int)cura[k], t2 = (int)curt[k];
            int x1 = xs[a], y1 = ys[a], x2 = xs[a2], y2 = ys[a2];
            int u1 = t / G, v1 = t % G, u2 = t2 / G, v2 = t2 % G;
            if (((x1 == x2) != (u1 == u2)) || ((y1 == y2) != (v1 == v2))) { ok = 0; break; }
            if (((x1 <  x2) != (u1 <  u2)) || ((y1 <  y2) != (v1 <  v2))) { ok = 0; break; }
        }
        if (!ok) { stack[frame]++; continue; }
        cura[frame] = (uint8_t)a; curt[frame] = (uint8_t)t;
        uint64_t child = pk;
        for (int i = 0; i <= frame; i++)
            child = (child & ~(0xffULL << (8 * (int)cura[i])))
                  | ((uint64_t)curt[i] << (8 * (int)cura[i]));
        if (child != pk) {
            double dmax = 0.0;
            for (int i = 0; i <= frame; i++) {
                int aa = (int)cura[i], tt = (int)curt[i];
                double d = site_dist_us_frame(xs[aa] * G + ys[aa], tt);
                if (d > dmax) dmax = d;
            }
            double w = move_dur_us(dmax) + 2.0 * TRANSFER_US * (double)(frame + 1);
            uint64_t nc_fp = pcost + (uint64_t)llround(w * FPQ);
            int fresh;
            uint64_t slot = slot_for(child, &fresh);
            if (fresh || nc_fp < cost_fp[slot]) {
                cost_fp[slot] = nc_fp; parent[slot] = pk; hpush(nc_fp, child);
            }
        }
        frame++; stack[frame] = ci + 1;
    }
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: cdijkstra N G caplog2 [t1 t2 ...] (target = comma sites, atoms in order)\n");
        return 2;
    }
    N = atoi(argv[1]); G = atoi(argv[2]);
    NS = G * G;
    cap = 1ULL << atoi(argv[3]); hmask = cap - 1;
    keys = calloc(cap, sizeof(uint64_t));
    cost_fp = calloc(cap, sizeof(uint64_t));
    parent = calloc(cap, sizeof(uint64_t));
    used = calloc(cap, 1);
    if (!keys || !cost_fp || !parent || !used) { fprintf(stderr, "OOM tables\n"); return 3; }
    ntargets = argc - 4;
    if (ntargets > 8) { fprintf(stderr, "max 8 targets\n"); return 2; }
    for (int t = 0; t < ntargets; t++) {
        uint64_t k = 0;
        char *s = argv[4 + t], *e;
        for (int a = 0; a < N; a++) {
            long v = strtol(s, &e, 10);
            if (e == s || v < 0 || v >= NS) { fprintf(stderr, "bad target token '%s'\n", s); return 2; }
            k |= (uint64_t)v << (8 * a);
            s = (*e == ',') ? e + 1 : e;
        }
        target_keys[t] = k; target_hit[t] = 0;
    }

    uint64_t start = 0;
    for (int a = 0; a < N; a++) start |= (uint64_t)a << (8 * a);
    int fresh;
    uint64_t sslot = slot_for(start, &fresh);
    cost_fp[sslot] = 0;
    hpush(0, start);

    double t0 = now_s();
    uint64_t settled = 0, stale = 0;
    int all_done = (ntargets == 0);
    while (hn) {
        HE e = hpop();
        uint64_t slot = find_slot(e.k);
        if (cost_fp[slot] != e.c) { stale++; continue; }
        settled++;
        for (int t = 0; t < ntargets; t++)
            if (!target_hit[t] && e.k == target_keys[t]) target_hit[t] = 1;
        all_done = 1;
        for (int t = 0; t < ntargets; t++) if (!target_hit[t]) { all_done = 0; break; }
        if (all_done) break;
        expand(e.k, e.c);
        if ((settled & 0xFFFFF) == 0)
            fprintf(stderr, "{\"settled\":%llu,\"heap\":%llu,\"wall_s\":%.1f}\n",
                    (unsigned long long)settled, (unsigned long long)hn, now_s() - t0);
    }
    double wall = now_s() - t0;

    printf("{\"n\":%d,\"g\":%d,\"settled\":%llu,\"stale\":%llu,\"heap_left\":%llu,"
           "\"wall_s\":%.2f,\"complete\":%d,\"targets\":[",
           N, G, (unsigned long long)settled, (unsigned long long)stale,
           (unsigned long long)hn, wall, all_done ? 1 : 0);
    for (int t = 0; t < ntargets; t++) {
        uint64_t tk = target_keys[t];
        uint64_t slot = find_slot(tk);
        if (t) printf(",");
        printf("{\"sites\":[");
        for (int a = 0; a < N; a++) printf(a ? ",%d" : "%d", (int)((tk >> (8*a)) & 0xff));
        int reached = (slot != UINT64_MAX) && (cost_fp[slot] != INF_COST) && target_hit[t];
        printf("],\"cost_us\":%.9f,\"reached\":%d", 
               reached ? (double)cost_fp[slot] / FPQ : -1.0, reached ? 1 : 0);
        if (reached) {
            printf(",\"batches\":[");
            uint64_t walk[80]; int wl = 0;
            uint64_t cur = tk;
            while (cur != start && wl < 80) { walk[wl++] = cur; cur = parent[find_slot(cur)]; }
            int firstb = 1;
            for (int i = wl - 1; i >= 0; i--) {
                uint64_t childk = walk[i];
                uint64_t pk = parent[find_slot(childk)];
                if (!firstb) printf(",");
                firstb = 0;
                printf("[");
                int fb = 1;
                for (int a = 0; a < N; a++) {
                    int cs = (int)((childk >> (8*a)) & 0xff);
                    int ps = (int)((pk >> (8*a)) & 0xff);
                    if (cs != ps) { if (!fb) printf(","); fb = 0; printf("[%d,%d]", a, cs); }
                }
                printf("]");
            }
            printf("]");
        }
        printf("}");
    }
    printf("]}\n");
    return 0;
}
