/* glue.c — implementation B of the R(4,5,n,e) census by vertex gluing.
 *
 * Per the frozen spec notes/census_gluing_spec.md:
 *   vertex 0 = v, 1..d = H vertices (catalog order), d+1..n-1 = K vertices
 *   (catalog order).  For each (H,K) catalog pair, DFS over cone sets
 *   S_k subseteq V(H), k in catalog order, under the complete constraint set.
 *
 * Usage: glue n e d kstart kend out.g6 out.csv datadir
 *   Enumerates stratum (n,e) at glue degree d, for K catalog lines
 *   [kstart,kend) (q = n-1-d; if q == 0 the synthesized empty K graph is
 *   "line 0", so pass kstart=0 kend=1).
 * csv rows: d,h_idx,k_idx,n_solutions  (no header; driver merges/sorts).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define MAXD 13
#define MAXQ 17
#define MAXVN 32
#define MAXPAIR 140   /* > C(16,2)=120 pairs (i,j) below a largest index k<=16 */

/* ---------------- graph6 ---------------- */

static int g6_read(const char *line, uint32_t adj[]) {
    int n = line[0] - 63;
    if (n < 1 || n > 62) return -1;
    for (int i = 0; i < n; i++) adj[i] = 0;
    const char *p = line + 1;
    int bitpos = 0;
    for (int j = 1; j < n; j++)
        for (int i = 0; i < j; i++) {
            int c = p[bitpos / 6] - 63;
            if (c < 0 || c > 63) return -1;
            if ((c >> (5 - bitpos % 6)) & 1) {
                adj[i] |= 1u << j;
                adj[j] |= 1u << i;
            }
            bitpos++;
        }
    return n;
}

static void g6_write(FILE *f, const uint32_t adj[], int n) {
    char buf[192];
    int len = 0;
    buf[len++] = (char)(63 + n);
    int acc = 0, nb = 0;
    for (int j = 1; j < n; j++)
        for (int i = 0; i < j; i++) {
            acc = (acc << 1) | (int)((adj[i] >> j) & 1u);
            if (++nb == 6) { buf[len++] = (char)(63 + acc); acc = 0; nb = 0; }
        }
    if (nb) buf[len++] = (char)(63 + (acc << (6 - nb)));
    buf[len++] = '\n';
    fwrite(buf, 1, (size_t)len, f);
}

/* ---------------- catalog loading ---------------- */

typedef struct {
    char **lines;
    int n;
} Catalog;

static Catalog load_catalog(const char *path) {
    Catalog c = {NULL, 0};
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(2); }
    size_t cap = 1024;
    c.lines = malloc(cap * sizeof(char *));
    char buf[512];
    while (fgets(buf, sizeof buf, f)) {
        size_t L = strlen(buf);
        while (L && (buf[L-1] == '\n' || buf[L-1] == '\r')) buf[--L] = 0;
        if (!L) continue;
        if ((size_t)c.n == cap) { cap *= 2; c.lines = realloc(c.lines, cap * sizeof(char *)); }
        c.lines[c.n++] = strdup(buf);
    }
    fclose(f);
    return c;
}

/* ---------------- H-side data ---------------- */

typedef struct {
    int d, eH;
    uint16_t adjH[MAXD];
    int capH[MAXD];       /* d - 1 - deg_H(h): remaining cone capacity of h */
    int HcapTotal;
    uint16_t full;
    uint16_t forb0;       /* h with capH == 0, forbidden from the start */
    uint8_t *mis;         /* mis[m] = max independent subset of H inside m */
    uint8_t *indep;       /* indep[m] = 1 iff m is independent in H */
    uint16_t *cands;      /* masks hitting every I4 of H, sorted by popcount */
    uint8_t *cpc;         /* popcounts of cands */
    int ncands;
    int cum[MAXD + 2];    /* cum[s] = first index with popcount >= s */
    int minCand;
} Hdat;

static void buildH(Hdat *H, const char *g6, int d) {
    uint32_t adj[MAXVN];
    int nn = g6_read(g6, adj);
    if (nn != d) { fprintf(stderr, "H catalog n=%d expected %d\n", nn, d); exit(2); }
    H->d = d;
    H->full = (uint16_t)((1u << d) - 1);
    int esum = 0;
    H->HcapTotal = 0;
    H->forb0 = 0;
    for (int i = 0; i < d; i++) {
        H->adjH[i] = (uint16_t)adj[i];
        int deg = __builtin_popcount(adj[i]);
        esum += deg;
        H->capH[i] = d - 1 - deg;
        H->HcapTotal += H->capH[i];
        if (H->capH[i] == 0) H->forb0 |= (uint16_t)(1u << i);
    }
    H->eH = esum / 2;
    int M = 1 << d;
    H->mis = malloc((size_t)M);
    H->indep = malloc((size_t)M);
    H->mis[0] = 0;
    H->indep[0] = 1;
    for (int m = 1; m < M; m++) {
        int i = __builtin_ctz((unsigned)m);
        int rest = m & (m - 1);
        int a = H->mis[rest];
        int b = 1 + H->mis[rest & ~(int)H->adjH[i]];
        H->mis[m] = (uint8_t)(a > b ? a : b);
        H->indep[m] = (uint8_t)(H->indep[rest] && !(H->adjH[i] & rest));
    }
    /* candidate masks: hit every independent 4-set of H
       <=> complement contains no I4 <=> mis[full^m] <= 3 */
    int cnt[MAXD + 2];
    memset(cnt, 0, sizeof cnt);
    int nc = 0;
    for (int m = 0; m < M; m++)
        if (H->mis[H->full ^ m] <= 3) { cnt[__builtin_popcount((unsigned)m)]++; nc++; }
    H->ncands = nc;
    H->cands = malloc((size_t)nc * sizeof(uint16_t));
    H->cpc = malloc((size_t)nc);
    int pos[MAXD + 2];
    int run = 0;
    for (int s = 0; s <= d + 1; s++) {
        H->cum[s] = run;
        pos[s] = run;
        if (s <= d) run += cnt[s];
    }
    H->cum[d + 1] = nc;
    for (int m = 0; m < M; m++)
        if (H->mis[H->full ^ m] <= 3) {
            int s = __builtin_popcount((unsigned)m);
            H->cands[pos[s]] = (uint16_t)m;
            H->cpc[pos[s]] = (uint8_t)s;
            pos[s]++;
        }
    H->minCand = d + 1;
    for (int s = 0; s <= d; s++)
        if (cnt[s]) { H->minCand = s; break; }
}

/* ---------------- K-side data ---------------- */

typedef struct {
    int q, eK;
    uint32_t adjK[MAXQ];
    int capK[MAXQ];          /* d - deg_K(k), indexed by actual K vertex */
    int ord[MAXQ];           /* DFS position p -> actual K vertex */
    int capP[MAXQ];          /* d - deg_K(ord[p]) */
    int capSuf[MAXQ + 1];    /* suffix sums of capP over positions */
    int feasible;            /* all capK >= 0 */
    /* constraints grouped by their largest DFS position p, checked when the
       vertex at position p is assigned; stored partners are ACTUAL vertices */
    int nE[MAXQ];  uint8_t Ej[MAXQ][MAXQ];                       /* K-edges */
    int nT[MAXQ];  uint8_t Ti[MAXQ][MAXPAIR], Tj[MAXQ][MAXPAIR]; /* K-triangles */
    int nP[MAXQ];  uint8_t Pj[MAXQ][MAXQ];                       /* K-indep pairs */
    int nI[MAXQ];  uint8_t Ii[MAXQ][MAXPAIR], Ij[MAXQ][MAXPAIR]; /* K-indep triples */
} Kdat;

static void buildK(Kdat *K, const char *g6, int q, int d) {
    uint32_t adj[MAXVN];
    int nn = g6_read(g6, adj);
    if (nn != q) { fprintf(stderr, "K catalog n=%d expected %d\n", nn, q); exit(2); }
    K->q = q;
    int esum = 0;
    K->feasible = 1;
    int deg[MAXQ];
    for (int i = 0; i < q; i++) {
        K->adjK[i] = adj[i];
        deg[i] = __builtin_popcount(adj[i]);
        esum += deg[i];
        K->capK[i] = d - deg[i];
        if (K->capK[i] < 0) K->feasible = 0;
    }
    K->eK = esum / 2;
    /* DFS order: most-constrained first = degree descending, tie by index */
    for (int p = 0; p < q; p++) K->ord[p] = p;
    for (int a = 1; a < q; a++) {           /* insertion sort, stable */
        int v = K->ord[a], b = a;
        while (b > 0 && deg[K->ord[b - 1]] < deg[v]) {
            K->ord[b] = K->ord[b - 1];
            b--;
        }
        K->ord[b] = v;
    }
    int pos[MAXQ];
    for (int p = 0; p < q; p++) {
        pos[K->ord[p]] = p;
        K->capP[p] = K->capK[K->ord[p]];
    }
    K->capSuf[q] = 0;
    for (int p = q - 1; p >= 0; p--) K->capSuf[p] = K->capSuf[p + 1] + K->capP[p];
    for (int p = 0; p < q; p++) K->nE[p] = K->nT[p] = K->nP[p] = K->nI[p] = 0;
    for (int j = 0; j < q; j++)
        for (int k = j + 1; k < q; k++) {
            int p = pos[j] > pos[k] ? pos[j] : pos[k];
            int other = pos[j] > pos[k] ? k : j;
            if ((adj[j] >> k) & 1) K->Ej[p][K->nE[p]++] = (uint8_t)other;
            else                   K->Pj[p][K->nP[p]++] = (uint8_t)other;
        }
    for (int i = 0; i < q; i++)
        for (int j = i + 1; j < q; j++)
            for (int k = j + 1; k < q; k++) {
                int aij = (adj[i] >> j) & 1, aik = (adj[i] >> k) & 1,
                    ajk = (adj[j] >> k) & 1;
                int all = aij && aik && ajk, none = !aij && !aik && !ajk;
                if (!all && !none) continue;
                int pi = pos[i], pj = pos[j], pk = pos[k];
                int p = pi, o1 = j, o2 = k;
                if (pj > p) { p = pj; o1 = i; o2 = k; }
                if (pk > p) { p = pk; o1 = i; o2 = j; }
                if (all) {
                    K->Ti[p][K->nT[p]] = (uint8_t)o1;
                    K->Tj[p][K->nT[p]] = (uint8_t)o2;
                    K->nT[p]++;
                } else {
                    K->Ii[p][K->nI[p]] = (uint8_t)o1;
                    K->Ij[p][K->nI[p]] = (uint8_t)o2;
                    K->nI[p]++;
                }
            }
}

/* ---------------- DFS ---------------- */

static int gN, gD, gQ, gCone;
static const Hdat *gH;
static const Kdat *gK;
static uint16_t gS[MAXQ];
static int gUsedH[MAXD];
static uint16_t gForb;
static int gSum, gHcapRem;
static long gPairSol;
static FILE *gOut;
static uint32_t gBase[MAXVN];   /* v + H + K adjacency, cone edges absent */

static void emit_solution(void) {
    uint32_t adj[MAXVN];
    memcpy(adj, gBase, sizeof(uint32_t) * (size_t)gN);
    for (int k = 0; k < gQ; k++) {
        uint16_t s = gS[k];
        int kv = 1 + gD + k;
        while (s) {
            int h = __builtin_ctz(s);
            s &= (uint16_t)(s - 1);
            adj[1 + h] |= 1u << kv;
            adj[kv] |= 1u << (1 + h);
        }
    }
    g6_write(gOut, adj, gN);
    gPairSol++;
}

static void dfs(int p) {
    if (p == gQ) { emit_solution(); return; }
    int k = gK->ord[p];   /* actual K vertex assigned at this position */
    int need = gCone - gSum;
    int lo = need - gK->capSuf[p + 1];
    if (lo < 0) lo = 0;
    int hi = need - gH->minCand * (gQ - 1 - p);
    if (hi > gK->capP[p]) hi = gK->capP[p];
    if (hi > gD) hi = gD;
    if (lo > hi) return;
    const uint16_t *cands = gH->cands;
    const uint8_t *cpc = gH->cpc;
    int i0 = gH->cum[lo], i1 = gH->cum[hi + 1];
    for (int idx = i0; idx < i1; idx++) {
        uint16_t m = cands[idx];
        if (m & gForb) continue;
        int pc = cpc[idx];
        if (need - pc > gHcapRem - pc) continue;   /* H-side capacity */
        /* K4 (2,2): S_j cap S_k independent in H, over K-edges */
        for (int t = 0; t < gK->nE[p]; t++)
            if (!gH->indep[gS[gK->Ej[p][t]] & m]) goto reject;
        /* K4 (1,3): empty triple intersection over K-triangles */
        for (int t = 0; t < gK->nT[p]; t++)
            if (gS[gK->Ti[p][t]] & gS[gK->Tj[p][t]] & m) goto reject;
        /* I5 (3,2): union hits every I3 of H, over K-independent pairs */
        for (int t = 0; t < gK->nP[p]; t++)
            if (gH->mis[gH->full ^ (gS[gK->Pj[p][t]] | m)] > 2) goto reject;
        /* I5 (2,3): union hits every I2 of H, over K-independent triples */
        for (int t = 0; t < gK->nI[p]; t++)
            if (gH->mis[gH->full ^ (gS[gK->Ii[p][t]] | gS[gK->Ij[p][t]] | m)] > 1)
                goto reject;
        /* accept: update H usage / forbidden and recurse */
        {
            gS[k] = m;
            gSum += pc;
            gHcapRem -= pc;
            uint16_t added = 0, s = m;
            while (s) {
                int h = __builtin_ctz(s);
                s &= (uint16_t)(s - 1);
                if (++gUsedH[h] == gH->capH[h]) added |= (uint16_t)(1u << h);
            }
            gForb |= added;
            dfs(p + 1);
            gForb &= (uint16_t)~added;
            s = m;
            while (s) {
                int h = __builtin_ctz(s);
                s &= (uint16_t)(s - 1);
                gUsedH[h]--;
            }
            gSum -= pc;
            gHcapRem += pc;
        }
    reject:;
    }
}

/* ---------------- main ---------------- */

int main(int argc, char **argv) {
    if (argc != 9) {
        fprintf(stderr, "usage: glue n e d kstart kend out.g6 out.csv datadir\n");
        return 2;
    }
    int n = atoi(argv[1]), e = atoi(argv[2]), d = atoi(argv[3]);
    int kstart = atoi(argv[4]), kend = atoi(argv[5]);
    const char *g6path = argv[6], *csvpath = argv[7], *datadir = argv[8];
    int q = n - 1 - d;
    if (d < 1 || d > MAXD || q < 0 || q > MAXQ) {
        fprintf(stderr, "bad d=%d (q=%d)\n", d, q);
        return 2;
    }
    char path[1024];
    snprintf(path, sizeof path, "%s/r35_%d.g6", datadir, d);
    Catalog hc = load_catalog(path);
    Hdat *HS = malloc((size_t)hc.n * sizeof(Hdat));
    for (int i = 0; i < hc.n; i++) buildH(&HS[i], hc.lines[i], d);

    Catalog kc = {NULL, 0};
    if (q > 0) {
        snprintf(path, sizeof path, "%s/r44_%d.g6", datadir, q);
        kc = load_catalog(path);
    } else {
        kc.n = 1;   /* synthesized empty graph on 0 vertices */
    }
    if (kstart < 0) kstart = 0;
    if (kend > kc.n) kend = kc.n;

    FILE *out = fopen(g6path, "w");
    FILE *csv = fopen(csvpath, "w");
    if (!out || !csv) { fprintf(stderr, "cannot open outputs\n"); return 2; }
    setvbuf(out, NULL, _IOFBF, 1 << 20);
    gOut = out;
    gN = n; gD = d; gQ = q;
    long total = 0;

    Kdat K;
    for (int kidx = kstart; kidx < kend; kidx++) {
        if (q > 0) {
            buildK(&K, kc.lines[kidx], q, d);
            if (!K.feasible) continue;
        } else {
            memset(&K, 0, sizeof K);
            K.q = 0;
            K.feasible = 1;
            K.capSuf[0] = 0;
        }
        gK = &K;
        for (int hidx = 0; hidx < hc.n; hidx++) {
            const Hdat *H = &HS[hidx];
            int cone = e - d - H->eH - K.eK;
            if (cone < 0 || cone > d * q) continue;
            if (cone > H->HcapTotal || cone > K.capSuf[0]) continue;
            if (H->minCand > d) continue;           /* no valid S_k at all */
            if (q > 0) {
                int bad = 0;
                for (int k = 0; k < q; k++)
                    if (K.capK[k] < H->minCand) { bad = 1; break; }
                if (bad) continue;
                if (cone < H->minCand * q) continue;
            }
            gH = H;
            /* base adjacency: v(0) - H(1..d) complete; H edges; K edges */
            memset(gBase, 0, sizeof gBase);
            for (int h = 0; h < d; h++) {
                gBase[0] |= 1u << (1 + h);
                gBase[1 + h] |= 1u;
                uint16_t s = H->adjH[h];
                while (s) {
                    int h2 = __builtin_ctz(s);
                    s &= (uint16_t)(s - 1);
                    gBase[1 + h] |= 1u << (1 + h2);
                }
            }
            for (int k = 0; k < q; k++) {
                uint32_t s = K.adjK[k];
                while (s) {
                    int k2 = __builtin_ctz(s);
                    s &= s - 1;
                    gBase[1 + d + k] |= 1u << (1 + d + k2);
                }
            }
            if (getenv("GLUE_DEBUG"))
                fprintf(stderr, "dfs pair h=%d k=%d cone=%d eH=%d eK=%d minCand=%d Hcap=%d Kcap=%d\n",
                        hidx, kidx, cone, H->eH, K.eK, H->minCand,
                        H->HcapTotal, K.capSuf[0]);
            gCone = cone;
            gPairSol = 0;
            gSum = 0;
            gHcapRem = H->HcapTotal;
            gForb = H->forb0;
            memset(gUsedH, 0, sizeof gUsedH);
            if (q == 0) {
                if (cone == 0) emit_solution();
            } else {
                dfs(0);
            }
            if (gPairSol > 0) {
                fprintf(csv, "%d,%d,%d,%ld\n", d, hidx, kidx, gPairSol);
                total += gPairSol;
            }
        }
    }
    fclose(out);
    fclose(csv);
    fprintf(stderr, "glue n=%d e=%d d=%d k[%d,%d): %ld solutions\n",
            n, e, d, kstart, kend, total);
    return 0;
}
