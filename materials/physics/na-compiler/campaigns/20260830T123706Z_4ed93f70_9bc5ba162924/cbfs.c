
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
