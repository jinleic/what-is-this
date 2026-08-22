/* C port of implementation A's enumerator core (census_gluing_spec.md).
 *
 * Enumerates Ramsey(4,5,n,e)-graph candidates by gluing R(3,5,d) x R(4,4,q),
 * q = n-1-d, over a max-degree center vertex, for ONE fixed d and a range of
 * H-catalog indices (parallelism = many processes over h-ranges).
 *
 * Output: raw graph6 lines (vertex 0 = v, 1..d = H catalog order, d+1.. = K
 * catalog order) and counts rows "d,h_idx,k_idx,count". No final Ramsey
 * verification here: the Python harness re-checks every emitted graph with
 * the audited independent checker before accepting it.
 *
 * Build: cc -O2 -o glue_census glue_census.c
 * Usage: glue_census n e d r35.g6 r44.g6 h_lo h_hi out.g6 out.csv [k_lo k_hi]
 *        (q = 0 supported by passing "-" as r44.g6; the optional K-range
 *         shards over the R(4,4) catalog instead — only that slice is kept
 *         in memory, counts.csv keeps global k indices)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 32
#define MAXH 13
#define MAXQ 17
#define HTAB (1 << MAXH)

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

static void g6_write(FILE *f, int n, const uint adj[]) {
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
    fwrite(buf, 1, len, f);
}

/* --------------------------------------------------------------- HSide -- */

static int d_h;                    /* |V(H)| = gluing degree d */
static uint h_adj[MAXH];
static int h_edges, h_minp, h_caph[MAXH];
static u8 mis[HTAB], indep_ok[HTAB], hit2_ok[HTAB], hit3_ok[HTAB], hit4_ok[HTAB];
static u16 by_pc_mask[HTAB];       /* masks grouped by popcount */
static int by_pc_off[MAXH + 2];    /* offsets into by_pc_mask per pc value */

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
    mis[0] = 0;
    static u8 hasedge[HTAB];
    hasedge[0] = 0;
    for (int m = 1; m <= full; m++) {
        int v = __builtin_ctz(m), rest = m & (m - 1);
        int a = mis[rest], b = 1 + mis[m & ~closed[v]];
        mis[m] = (u8)(b > a ? b : a);
        hasedge[m] = (u8)(hasedge[rest] || (adj[v] & rest));
    }
    int cnt_pc[MAXH + 1];
    memset(cnt_pc, 0, sizeof cnt_pc);
    for (int m = 0; m <= full; m++) {
        int c = mis[full & ~m];
        indep_ok[m] = !hasedge[m];
        hit4_ok[m] = (u8)(c <= 3);
        hit3_ok[m] = (u8)(c <= 2);
        hit2_ok[m] = (u8)(c <= 1);
        if (hit4_ok[m]) cnt_pc[__builtin_popcount(m)]++;
    }
    by_pc_off[0] = 0;
    for (int p = 0; p <= d; p++) by_pc_off[p + 1] = by_pc_off[p] + cnt_pc[p];
    int fill[MAXH + 1];
    for (int p = 0; p <= d; p++) fill[p] = by_pc_off[p];
    for (int m = 0; m <= full; m++)
        if (hit4_ok[m]) by_pc_mask[fill[__builtin_popcount(m)]++] = (u16)m;
    h_minp = -1;
    for (int p = 0; p <= d; p++)
        if (cnt_pc[p]) { h_minp = p; break; }
}

/* --------------------------------------------------------------- KSide -- */

typedef struct { u8 type, i1, i2; } Check;   /* types: 0 (2,2)  1 (3,2)  2 (1,3)  3 (2,3) */

typedef struct {
    int q, edges;
    uint adj[MAXQ];
    int deg[MAXQ], caps[MAXQ], perm[MAXQ], skip;
    Check *checks[MAXQ];
    int nchecks[MAXQ];
    Check store[1400];
} KSide;

static int glue_d;

static void build_kside(KSide *K, const uint adj[], int q) {
    K->q = q;
    memcpy(K->adj, adj, sizeof(uint) * MAXQ);
    K->edges = 0;
    K->skip = 0;
    for (int v = 0; v < q; v++) {
        K->deg[v] = __builtin_popcount(adj[v]);
        K->edges += K->deg[v];
        K->caps[v] = glue_d - K->deg[v];
        if (K->caps[v] < 0) K->skip = 1;
    }
    K->edges /= 2;
    /* perm: caps ascending (stable) */
    for (int i = 0; i < q; i++) K->perm[i] = i;
    for (int i = 1; i < q; i++)
        for (int j = i; j > 0 && K->caps[K->perm[j]] < K->caps[K->perm[j - 1]]; j--) {
            int t = K->perm[j]; K->perm[j] = K->perm[j - 1]; K->perm[j - 1] = t;
        }
    int pos[MAXQ];
    for (int i = 0; i < q; i++) pos[K->perm[i]] = i;
    int cnt[MAXQ];
    memset(cnt, 0, sizeof cnt);
    int total = 0;
    /* first pass: count, second: fill */
    for (int pass = 0; pass < 2; pass++) {
        if (pass == 1) {
            int off = 0;
            for (int i = 0; i < q; i++) {
                K->checks[i] = K->store + off;
                off += cnt[i];
                K->nchecks[i] = 0;
            }
            total = off;
            if (total > 1400) { fprintf(stderr, "check overflow\n"); exit(1); }
        }
        for (int a = 0; a < q; a++)
            for (int b = a + 1; b < q; b++) {
                int e_ab = (K->adj[a] >> b) & 1;
                int i = pos[a] < pos[b] ? pos[a] : pos[b];
                int j = pos[a] < pos[b] ? pos[b] : pos[a];
                if (pass == 0) cnt[j]++;
                else {
                    Check *c = &K->checks[j][K->nchecks[j]++];
                    c->type = e_ab ? 0 : 1;
                    c->i1 = (u8)i; c->i2 = 0;
                }
                for (int cc = b + 1; cc < q; cc++) {
                    int e_ac = (K->adj[a] >> cc) & 1, e_bc = (K->adj[b] >> cc) & 1;
                    int tri = e_ab && e_ac && e_bc;
                    int ind = !e_ab && !e_ac && !e_bc;
                    if (!tri && !ind) continue;
                    int p1 = pos[a], p2 = pos[b], p3 = pos[cc], t;
                    if (p1 > p2) { t = p1; p1 = p2; p2 = t; }
                    if (p2 > p3) { t = p2; p2 = p3; p3 = t; }
                    if (p1 > p2) { t = p1; p1 = p2; p2 = t; }
                    if (pass == 0) cnt[p3]++;
                    else {
                        Check *c = &K->checks[p3][K->nchecks[p3]++];
                        c->type = tri ? 2 : 3;
                        c->i1 = (u8)p1; c->i2 = (u8)p2;
                    }
                }
            }
        if (pass == 0) { /* prefix layout happens at pass 1 start */ }
    }
    (void)total;
}

/* ----------------------------------------------------------------- DFS -- */

static int q_g, n_g, e_g, c_target;
static KSide *K_g;
static int sufmin[MAXQ + 1], sufmax[MAXQ + 1];
static u16 S[MAXQ];
static int used[MAXH];
static long nsol_pair;
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
    u16 Sk[MAXQ];
    for (int i = 0; i < q; i++) Sk[K_g->perm[i]] = S[i];
    for (int k = 0; k < q; k++) {
        uint a = K_g->adj[k];
        while (a) { int b = __builtin_ctz(a); adj[1 + d + k] |= 1u << (1 + d + b); a &= a - 1; }
        uint m = Sk[k];
        while (m) {
            int h = __builtin_ctz(m);
            adj[1 + d + k] |= 1u << (1 + h);
            adj[1 + h] |= 1u << (1 + d + k);
            m &= m - 1;
        }
    }
    g6_write(out_g6, n_g, adj);
    nsol_pair++;
}

static void dfs(int i, int remaining, uint fullmask) {
    if (i == q_g) { emit_solution(); return; }
    int k = K_g->perm[i];
    int hi = K_g->caps[k], t = remaining - sufmin[i + 1];
    if (t < hi) hi = t;
    int lo = h_minp, t2 = remaining - sufmax[i + 1];
    if (t2 > lo) lo = t2;
    Check *chk = K_g->checks[i];
    int nchk = K_g->nchecks[i];
    for (int p = hi; p >= lo; p--) {
        const u16 *mm = by_pc_mask + by_pc_off[p];
        int cnt = by_pc_off[p + 1] - by_pc_off[p];
        for (int ci = 0; ci < cnt; ci++) {
            uint m = mm[ci];
            if (m & fullmask) continue;
            int ok = 1;
            for (int c = 0; c < nchk; c++) {
                switch (chk[c].type) {
                case 0: if (!indep_ok[S[chk[c].i1] & m]) ok = 0; break;
                case 1: if (!hit3_ok[S[chk[c].i1] | m]) ok = 0; break;
                case 2: if (S[chk[c].i1] & S[chk[c].i2] & m) ok = 0; break;
                default: if (!hit2_ok[S[chk[c].i1] | S[chk[c].i2] | m]) ok = 0;
                }
                if (!ok) break;
            }
            if (!ok) continue;
            S[i] = (u16)m;
            uint newfull = fullmask, x = m;
            int bad = 0;
            while (x) {
                int h = __builtin_ctz(x);
                if (++used[h] > h_caph[h]) bad = 1;
                else if (used[h] == h_caph[h]) newfull |= 1u << h;
                x &= x - 1;
            }
            if (!bad) dfs(i + 1, remaining - p, newfull);
            x = m;
            while (x) { used[__builtin_ctz(x)]--; x &= x - 1; }
        }
    }
}

/* ---------------------------------------------------------------- main -- */

int main(int argc, char **argv) {
    if (argc != 10 && argc != 12) {
        fprintf(stderr, "usage: %s n e d r35.g6 r44.g6 h_lo h_hi out.g6 out.csv [k_lo k_hi]\n", argv[0]);
        return 1;
    }
    n_g = atoi(argv[1]); e_g = atoi(argv[2]); glue_d = atoi(argv[3]);
    int h_lo = atoi(argv[6]), h_hi = atoi(argv[7]);
    long k_lo = 0, k_hi = -1;  /* -1: no K sharding */
    if (argc == 12) { k_lo = atol(argv[10]); k_hi = atol(argv[11]); }
    int q = n_g - 1 - glue_d;
    /* load H catalog lines */
    static char hline[100000][32];
    int nh = 0;
    FILE *f = fopen(argv[4], "r");
    if (!f) { perror(argv[4]); return 1; }
    while (fgets(hline[nh], 32, f)) if (hline[nh][0] > ' ') nh++;
    fclose(f);
    /* load + precompute K catalog (slice [k_lo,k_hi) if sharded) */
    static KSide *ks;
    long nk = 0, k_base = 0;
    if (q == 0) {
        ks = malloc(sizeof(KSide));
        uint none[MAXQ] = {0};
        build_kside(&ks[0], none, 0);
        nk = 1;
    } else {
        long cap = (k_hi >= 0) ? (k_hi - k_lo) : 1600000;
        ks = malloc(sizeof(KSide) * (size_t)cap);
        if (!ks) { fprintf(stderr, "K alloc failed\n"); return 1; }
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
                build_kside(&ks[nk++], adj, q);
            }
            idx++;
            if (k_hi >= 0 && idx >= k_hi) break;
        }
        fclose(f);
    }
    out_g6 = fopen(argv[8], "w");
    FILE *out_csv = fopen(argv[9], "w");
    if (!out_g6 || !out_csv) { perror("out"); return 1; }
    q_g = q;
    uint adj[MAXN];
    for (int h = h_lo; h < h_hi && h < nh; h++) {
        int d = g6_parse(hline[h], adj);
        if (d != glue_d) { fprintf(stderr, "H order mismatch\n"); return 1; }
        build_hside(adj, d);
        if (h_minp < 0) continue;
        for (long k = 0; k < nk; k++) {
            KSide *K = &ks[k];
            if (K->skip) continue;
            c_target = e_g - glue_d - h_edges - K->edges;
            if (c_target < 0 || c_target > glue_d * q) continue;
            int okcaps = 1;
            for (int i = 0; i < q; i++)
                if (K->caps[i] < h_minp) { okcaps = 0; break; }
            if (!okcaps) continue;
            sufmin[q] = sufmax[q] = 0;
            for (int i = q - 1; i >= 0; i--) {
                sufmin[i] = sufmin[i + 1] + h_minp;
                sufmax[i] = sufmax[i + 1] + K->caps[K->perm[i]];
            }
            if (c_target < sufmin[0] || c_target > sufmax[0]) continue;
            K_g = K;
            memset(used, 0, sizeof used);
            uint init_full = 0;
            for (int hh = 0; hh < d_h; hh++)
                if (h_caph[hh] == 0) init_full |= 1u << hh;
            nsol_pair = 0;
            dfs(0, c_target, init_full);
            if (nsol_pair)
                fprintf(out_csv, "%d,%d,%ld,%ld\n", glue_d, h, k_base + k, nsol_pair);
        }
    }
    fclose(out_g6);
    fclose(out_csv);
    return 0;
}
