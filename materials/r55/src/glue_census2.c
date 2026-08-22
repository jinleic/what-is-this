/* v2 of the census enumerator (see glue_census.c for the contract, which is
 * identical): adds forward checking via precomputed binary-constraint
 * compatibility bitsets over the candidate universe, plus fail-first dynamic
 * variable ordering. Ternary constraints (K-triangles, K-independent-triples)
 * are checked when their last member is assigned.
 *
 * Also applies a minimum-degree bound dmin (10th arg via env GLUE_DMIN or
 * computed by the harness): deleting any vertex of a Ramsey(4,5,n,e)-graph
 * leaves a Ramsey(4,5,n-1)-graph, so e - E(4,5,n-1) <= deg(v) for every v.
 * This is a proved necessary condition (pure pruning; results unchanged).
 *
 * Build: cc -O2 -o glue_census2 glue_census2.c
 * Usage: glue_census2 n e d r35.g6 r44.g6 h_lo h_hi out.g6 out.csv [k_lo k_hi]
 *        env GLUE_DMIN=<int> enables the min-degree prune (default 0 = off)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXN 32
#define MAXH 13
#define MAXQ 17
#define HTAB (1 << MAXH)
#define MAXU HTAB              /* candidate universe cap */
#define MAXW (MAXU / 64)       /* bitset words */

typedef unsigned long long u64;
typedef unsigned int uint;
typedef unsigned short u16;
typedef unsigned char u8;

/* ------------------------------------------------------------- graph6 -- */

static int g6_parse(const char *s, uint adj[MAXN]) {
    int n = s[0] - 63, idx = 1;
    if (n < 0 || n > 62) { fprintf(stderr, "bad g6 n\n"); exit(1); }
    memset(adj, 0, sizeof(uint) * MAXN);
    int i = 0, j = 1, need = n * (n - 1) / 2, bits = 0;
    while (bits < need) {
        int b = s[idx++] - 63;
        if (b < 0 || b > 63) { fprintf(stderr, "bad g6 byte\n"); exit(1); }
        for (int sh = 5; sh >= 0 && bits < need; sh--, bits++) {
            if ((b >> sh) & 1) { adj[i] |= 1u << j; adj[j] |= 1u << i; }
            if (++i == j) { i = 0; j++; }
        }
    }
    return n;
}

/* per-pair solution buffer: flushed only if the pair completes un-aborted */
#define PAIRBUF_CAP (4 << 20)
static char pairbuf[PAIRBUF_CAP];
static size_t pairbuf_len;

static void g6_buffer(int n, const uint adj[]) {
    char buf[80];
    int len = 0;
    buf[len++] = (char)(n + 63);
    int acc = 0, nb = 0;
    for (int j = 1; j < n; j++)
        for (int i = 0; i < j; i++) {
            acc = (acc << 1) | ((adj[i] >> j) & 1);
            if (++nb == 6) { buf[len++] = (char)(acc + 63); acc = 0; nb = 0; }
        }
    if (nb) buf[len++] = (char)((acc << (6 - nb)) + 63);
    buf[len++] = '\n';
    if (pairbuf_len + len > PAIRBUF_CAP) return;  /* overflow -> pair aborts */
    memcpy(pairbuf + pairbuf_len, buf, len);
    pairbuf_len += len;
}

/* --------------------------------------------------------------- HSide -- */

static int d_h, h_edges, h_minp, h_caph[MAXH];
static uint h_adj[MAXH];
static u8 mis[HTAB], indep_ok[HTAB], hit2_ok[HTAB], hit3_ok[HTAB], hit4_ok[HTAB];

static int NU, NW;             /* universe size, bitset words */
static u16 umask[MAXU];
static u8 upc[MAXU];
static u64 *Redge, *Ripair;    /* NU rows x NW words each */
static u64 cumpc[MAXH + 1][MAXW]; /* candidates with pc <= c */
static u64 avoid_h[MAXH][MAXW];   /* universe masks NOT containing h */

static void build_hside(const uint adj[], int d) {
    d_h = d;
    memcpy(h_adj, adj, sizeof(uint) * MAXH);
    int full = (1 << d) - 1;
    uint closed[MAXH];
    h_edges = 0;
    for (int v = 0; v < d; v++) {
        closed[v] = adj[v] | (1u << v);
        h_edges += __builtin_popcount(adj[v]);
        h_caph[v] = d - 1 - __builtin_popcount(adj[v]);
    }
    h_edges /= 2;
    static u8 hasedge[HTAB];
    mis[0] = 0; hasedge[0] = 0;
    for (int m = 1; m <= full; m++) {
        int v = __builtin_ctz(m), rest = m & (m - 1);
        int a = mis[rest], b = 1 + mis[m & ~closed[v]];
        mis[m] = (u8)(b > a ? b : a);
        hasedge[m] = (u8)(hasedge[rest] || (adj[v] & rest));
    }
    NU = 0;
    for (int m = 0; m <= full; m++) {
        int c = mis[full & ~m];
        indep_ok[m] = !hasedge[m];
        hit4_ok[m] = (u8)(c <= 3);
        hit3_ok[m] = (u8)(c <= 2);
        hit2_ok[m] = (u8)(c <= 1);
        if (hit4_ok[m]) { umask[NU] = (u16)m; upc[NU] = (u8)__builtin_popcount(m); NU++; }
    }
    NW = (NU + 63) / 64;
    h_minp = NU ? 99 : -1;
    for (int i = 0; i < NU; i++) if (upc[i] < h_minp) h_minp = upc[i];
    memset(cumpc, 0, sizeof cumpc);
    for (int c = 0; c <= d; c++)
        for (int i = 0; i < NU; i++)
            if (upc[i] <= c) cumpc[c][i >> 6] |= 1ull << (i & 63);
    memset(avoid_h, 0, sizeof avoid_h);
    for (int h = 0; h < d; h++)
        for (int i = 0; i < NU; i++)
            if (!((umask[i] >> h) & 1)) avoid_h[h][i >> 6] |= 1ull << (i & 63);
    /* relations */
    for (int i = 0; i < NU; i++) {
        u64 *re = Redge + (size_t)i * NW, *ri = Ripair + (size_t)i * NW;
        memset(re, 0, NW * 8);
        memset(ri, 0, NW * 8);
        int mi = umask[i];
        for (int j = 0; j < NU; j++) {
            int mj = umask[j];
            if (indep_ok[mi & mj]) re[j >> 6] |= 1ull << (j & 63);
            if (hit3_ok[mi | mj]) ri[j >> 6] |= 1ull << (j & 63);
        }
    }
}

/* --------------------------------------------------------------- KSide -- */

typedef struct { u8 a, b, c, type; } Tern;   /* type: 0 triangle, 1 itriple */

static int q_g, n_g, e_g, glue_d, dmin_g;
static int lov[MAXQ], need_h[MAXH];
static uint k_adj[MAXQ];
static int k_deg[MAXQ], k_caps[MAXQ], k_edges, k_skip;
static u8 btype[MAXQ][MAXQ];                 /* 0 none, 1 edge, 2 ipair */
static Tern tern[1000];
static int ntern, tid_of[MAXQ][160], ntid[MAXQ];

static void build_kside(const uint adj[], int q) {
    q_g = q;
    memcpy(k_adj, adj, sizeof(uint) * MAXQ);
    k_edges = 0; k_skip = 0; ntern = 0;
    for (int v = 0; v < q; v++) {
        k_deg[v] = __builtin_popcount(adj[v]);
        k_edges += k_deg[v];
        k_caps[v] = glue_d - k_deg[v];
        if (k_caps[v] < 0) k_skip = 1;
        ntid[v] = 0;
    }
    k_edges /= 2;
    for (int a = 0; a < q; a++)
        for (int b = a + 1; b < q; b++) {
            int e_ab = (adj[a] >> b) & 1;
            btype[a][b] = btype[b][a] = (u8)(e_ab ? 1 : 2);
            for (int c = b + 1; c < q; c++) {
                int e_ac = (adj[a] >> c) & 1, e_bc = (adj[b] >> c) & 1;
                int tri = e_ab && e_ac && e_bc, ind = !e_ab && !e_ac && !e_bc;
                if (!tri && !ind) continue;
                tern[ntern] = (Tern){(u8)a, (u8)b, (u8)c, (u8)(tri ? 0 : 1)};
                tid_of[a][ntid[a]++] = ntern;
                tid_of[b][ntid[b]++] = ntern;
                tid_of[c][ntid[c]++] = ntern;
                ntern++;
            }
        }
}

/* ----------------------------------------------------------------- DFS -- */

static u64 alive[MAXQ + 1][MAXQ][MAXW];
static int assigned[MAXQ], Sv[MAXQ], used[MAXH];
static long nsol_pair, nodes_pair, node_limit;
static int aborted_pair;
static FILE *out_g6;

static void emit_solution(void) {
    uint adj[MAXN];
    memset(adj, 0, sizeof adj);
    int d = d_h, q = q_g;
    for (int h = 0; h < d; h++) {
        adj[0] |= 1u << (1 + h);
        adj[1 + h] |= 1u;
        uint a = h_adj[h];
        while (a) { int b = __builtin_ctz(a); adj[1 + h] |= 1u << (1 + b); a &= a - 1; }
    }
    for (int k = 0; k < q; k++) {
        uint a = k_adj[k];
        while (a) { int b = __builtin_ctz(a); adj[1 + d + k] |= 1u << (1 + d + b); a &= a - 1; }
        uint m = (uint)Sv[k];
        while (m) {
            int h = __builtin_ctz(m);
            adj[1 + d + k] |= 1u << (1 + h);
            adj[1 + h] |= 1u << (1 + d + k);
            m &= m - 1;
        }
    }
    g6_buffer(n_g, adj);
    nsol_pair++;
    if (pairbuf_len + 80 > PAIRBUF_CAP) aborted_pair = 1;
}

static void dfs(int level, int remaining, uint fullmask) {
    int q = q_g;
    if (aborted_pair) return;
    if (++nodes_pair > node_limit) { aborted_pair = 1; return; }
    if (level == q) {
        /* terminal: the H-side min-degree contract must hold exactly;
         * checking only at interior nodes (below) lets terminal assignments
         * emit used[h] < need_h[h].  Fixed 2026-08-15; see PROGRESS.md. */
        if (dmin_g)
            for (int h = 0; h < d_h; h++)
                if (used[h] < need_h[h]) return;
        emit_solution();
        return;
    }
    /* min-degree deficit prune: cone edges still owed to H vertices */
    if (dmin_g) {
        int deficit = 0;
        for (int h = 0; h < d_h; h++) {
            int nd = need_h[h] - used[h];
            if (nd > 0) deficit += nd;
        }
        if (deficit > remaining) return;
    }
    /* fail-first: pick unassigned var with fewest alive candidates */
    int best = -1, bestcnt = 1 << 30;
    for (int v = 0; v < q; v++) {
        if (assigned[v]) continue;
        int cnt = 0;
        const u64 *av = alive[level][v];
        for (int w = 0; w < NW; w++) cnt += __builtin_popcountll(av[w]);
        if (cnt == 0) return;
        if (cnt < bestcnt) { bestcnt = cnt; best = v; }
    }
    int j = best;
    /* budget bounds for the remaining vars other than j */
    int rest_min = 0, rest_max = 0;
    for (int v = 0; v < q; v++)
        if (!assigned[v] && v != j) { rest_min += lov[v]; rest_max += k_caps[v]; }
    int lo = remaining - rest_max, hi = remaining - rest_min;
    if (lo < lov[j]) lo = lov[j];
    if (hi > k_caps[j]) hi = k_caps[j];
    if (lo > hi) return;
    assigned[j] = 1;
    const u64 *avj = alive[level][j];
    for (int w = 0; w < NW; w++) {
        u64 bits = avj[w];
        while (bits) {
            int ci = (w << 6) + __builtin_ctzll(bits);
            bits &= bits - 1;
            int p = upc[ci];
            if (p < lo || p > hi) continue;
            uint m = umask[ci];
            if (m & fullmask) continue;
            /* ternary constraints whose other two members are assigned */
            int ok = 1;
            for (int t = 0; t < ntid[j] && ok; t++) {
                Tern *T = &tern[tid_of[j][t]];
                int a, b;
                if (T->a == j) { a = T->b; b = T->c; }
                else if (T->b == j) { a = T->a; b = T->c; }
                else { a = T->a; b = T->b; }
                if (!assigned[a] || !assigned[b]) continue;
                if (T->type == 0) { if (Sv[a] & Sv[b] & (int)m) ok = 0; }
                else { if (!hit2_ok[Sv[a] | Sv[b] | (int)m]) ok = 0; }
            }
            if (!ok) continue;
            /* degree caps on H side */
            uint newfull = fullmask, x = m;
            int bad = 0;
            while (x) {
                int h = __builtin_ctz(x);
                if (++used[h] > h_caph[h]) bad = 1;
                else if (used[h] == h_caph[h]) newfull |= 1u << h;
                x &= x - 1;
            }
            if (!bad) {
                /* forward-check: intersect future alive sets; also kill
                 * candidates touching newly saturated H vertices */
                uint newbits = newfull & ~fullmask;
                int empty = 0;
                for (int v = 0; v < q && !empty; v++) {
                    if (assigned[v]) continue;
                    const u64 *src = alive[level][v];
                    u64 *dst = alive[level + 1][v];
                    const u64 *rel = (btype[j][v] == 1 ? Redge : Ripair)
                                     + (size_t)ci * NW;
                    u64 any = 0;
                    for (int w2 = 0; w2 < NW; w2++) {
                        dst[w2] = src[w2] & rel[w2];
                        any |= dst[w2];
                    }
                    uint nb = newbits;
                    while (nb && any) {
                        int hh = __builtin_ctz(nb);
                        nb &= nb - 1;
                        any = 0;
                        for (int w2 = 0; w2 < NW; w2++) {
                            dst[w2] &= avoid_h[hh][w2];
                            any |= dst[w2];
                        }
                    }
                    if (!any) empty = 1;
                }
                if (!empty) {
                    Sv[j] = (int)m;
                    dfs(level + 1, remaining - p, newfull);
                }
            }
            x = m;
            while (x) { used[__builtin_ctz(x)]--; x &= x - 1; }
        }
    }
    assigned[j] = 0;
}

/* ---------------------------------------------------------------- main -- */

int main(int argc, char **argv) {
    if (argc != 10 && argc != 12) {
        fprintf(stderr, "usage: %s n e d r35.g6 r44.g6 h_lo h_hi out.g6 out.csv [k_lo k_hi]\n", argv[0]);
        return 1;
    }
    n_g = atoi(argv[1]); e_g = atoi(argv[2]); glue_d = atoi(argv[3]);
    int h_lo = atoi(argv[6]), h_hi = atoi(argv[7]);
    long k_lo = 0, k_hi = -1;
    if (argc == 12) { k_lo = atol(argv[10]); k_hi = atol(argv[11]); }
    const char *dm = getenv("GLUE_DMIN");
    dmin_g = dm ? atoi(dm) : 0;
    const char *nl = getenv("GLUE_NODE_LIMIT");
    node_limit = nl ? atol(nl) : (1L << 62);
    int q = n_g - 1 - glue_d;
    Redge = malloc((size_t)MAXU * MAXW * 8);
    Ripair = malloc((size_t)MAXU * MAXW * 8);
    if (!Redge || !Ripair) { fprintf(stderr, "alloc\n"); return 1; }
    static char hline[100000][32];
    int nh = 0;
    FILE *f = fopen(argv[4], "r");
    if (!f) { perror(argv[4]); return 1; }
    while (fgets(hline[nh], 32, f)) if (hline[nh][0] > ' ') nh++;
    fclose(f);
    /* K catalog slice: store raw adjacency, KSide built per pair (cheap) */
    static uint (*kadjs)[MAXQ];
    long nk = 0, k_base = 0;
    if (q == 0) {
        kadjs = malloc(sizeof(uint) * MAXQ);
        memset(kadjs[0], 0, sizeof(uint) * MAXQ);
        nk = 1;
    } else {
        long cap = (k_hi >= 0) ? (k_hi - k_lo) : 1600000;
        kadjs = malloc(sizeof(uint) * MAXQ * (size_t)cap);
        char line[64];
        f = fopen(argv[5], "r");
        if (!f) { perror(argv[5]); return 1; }
        uint adj[MAXN];
        long idx = 0;
        k_base = (k_hi >= 0) ? k_lo : 0;
        while (fgets(line, 64, f)) {
            if (line[0] <= ' ') continue;
            if (k_hi < 0 || (idx >= k_lo && idx < k_hi)) {
                int qq = g6_parse(line, adj);
                if (qq != q) { fprintf(stderr, "K order mismatch\n"); return 1; }
                memcpy(kadjs[nk++], adj, sizeof(uint) * MAXQ);
            }
            idx++;
            if (k_hi >= 0 && idx >= k_hi) break;
        }
        fclose(f);
    }
    out_g6 = fopen(argv[8], "w");
    FILE *out_csv = fopen(argv[9], "w");
    if (!out_g6 || !out_csv) { perror("out"); return 1; }
    int verbose = getenv("GLUE_VERBOSE") != NULL;
    /* GLUE_CSV_ALL: emit a CSV row for EVERY pair that enters the DFS,
     * including completed zero-solution pairs (the dominant cost population;
     * required for fitting handoff policies).  Default stays sparse
     * (deferred + solution-bearing rows only) so frozen pipelines are
     * unchanged.  Pairs killed by the pre-DFS feasibility checks above never
     * enter the DFS and are never emitted in either mode. */
    int csv_all = getenv("GLUE_CSV_ALL") != NULL;
    long cum_nodes = 0;
    clock_t tstart = clock();
    uint adj[MAXN];
    for (int h = h_lo; h < h_hi && h < nh; h++) {
        int d = g6_parse(hline[h], adj);
        if (d != glue_d) { fprintf(stderr, "H order mismatch\n"); return 1; }
        build_hside(adj, d);
        if (h_minp < 0 || h_minp > 90) continue;
        int h_ok = 1;
        for (int hh = 0; hh < d; hh++) {
            int nd = dmin_g - 1 - __builtin_popcount(adj[hh]);
            need_h[hh] = nd > 0 ? nd : 0;
            if (need_h[hh] > h_caph[hh]) h_ok = 0;
        }
        if (!h_ok) continue;
        for (long k = 0; k < nk; k++) {
            build_kside(kadjs[k], q);
            if (k_skip) continue;
            int c_target = e_g - glue_d - h_edges - k_edges;
            if (c_target < 0 || c_target > glue_d * q) continue;
            int okcaps = 1, smin = 0, smax = 0;
            for (int i = 0; i < q; i++) {
                lov[i] = h_minp;
                if (dmin_g && dmin_g - k_deg[i] > lov[i]) lov[i] = dmin_g - k_deg[i];
                if (k_caps[i] < lov[i]) { okcaps = 0; break; }
                smin += lov[i];
                smax += k_caps[i];
            }
            if (!okcaps || c_target < smin || c_target > smax) continue;
            if (q == 0 && c_target != 0) continue;
            /* init alive sets: candidates with lov_v <= pc <= cap_v */
            for (int v = 0; v < q; v++) {
                int cv = k_caps[v] > d_h ? d_h : k_caps[v];
                memcpy(alive[0][v], cumpc[cv], NW * 8);
                if (lov[v] > 0)
                    for (int w = 0; w < NW; w++)
                        alive[0][v][w] &= ~cumpc[lov[v] - 1][w];
                assigned[v] = 0;
            }
            memset(used, 0, sizeof used);
            uint init_full = 0;
            for (int hh = 0; hh < d_h; hh++)
                if (h_caph[hh] == 0) init_full |= 1u << hh;
            nsol_pair = 0;
            nodes_pair = 0;
            aborted_pair = 0;
            pairbuf_len = 0;
            dfs(0, c_target, init_full);
            cum_nodes += nodes_pair;
            if (aborted_pair) {
                /* pair exceeded the node budget: hand to the SAT path */
                fprintf(out_csv, "%d,%d,%ld,-1,%ld\n", glue_d, h, k_base + k, nodes_pair);
            } else {
                if (nsol_pair)
                    fwrite(pairbuf, 1, pairbuf_len, out_g6);
                if (nsol_pair || csv_all)
                    fprintf(out_csv, "%d,%d,%ld,%ld,%ld\n", glue_d, h, k_base + k,
                            nsol_pair, nodes_pair);
            }
            if (verbose && (k & 1023) == 1023)
                fprintf(stderr, "  h=%d k=%ld cum_nodes=%ld t=%.0fs\n", h, k_base + k,
                        cum_nodes, (double)(clock() - tstart) / CLOCKS_PER_SEC);
        }
    }
    fprintf(out_csv, "#nodes,%ld\n", cum_nodes);
    fclose(out_g6);
    fclose(out_csv);
    return 0;
}
