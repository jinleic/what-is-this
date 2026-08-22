/* verify.c — independent full-graph checker for implementation B.
 *
 * Reads a graph6 file and verifies every line is a graph on n vertices with
 * exactly e edges, no K4, and no independent 5-set.  Deliberately written
 * separately from glue.c's constraint machinery (safety net per spec).
 *
 * Usage: verify n e file.g6
 * Prints "checked <N> graphs, failures <F>"; exit 1 if F > 0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int parse_g6(const char *s, size_t len, uint32_t adj[], int *n_out) {
    if (len < 1) return -1;
    int n = s[0] - 63;
    if (n < 1 || n > 30) return -1;
    int nbits = n * (n - 1) / 2;
    size_t nchars = (size_t)((nbits + 5) / 6);
    if (len != 1 + nchars) return -1;
    for (int i = 0; i < n; i++) adj[i] = 0;
    /* edge (i,j), i<j, is bit number t = j*(j-1)/2 + i (column-major upper
       triangle), stored MSB-first in 6-bit chars */
    for (int j = 1; j < n; j++) {
        int rowbase = j * (j - 1) / 2;
        for (int i = 0; i < j; i++) {
            int t = rowbase + i;
            int c = s[1 + t / 6] - 63;
            if (c < 0 || c > 63) return -1;
            if (c & (1 << (5 - t % 6))) {
                adj[i] |= 1u << j;
                adj[j] |= 1u << i;
            }
        }
    }
    *n_out = n;
    return 0;
}

static int has_k4(const uint32_t adj[], int n) {
    for (int a = 0; a < n; a++) {
        uint32_t na = adj[a] & ~((1u << (a + 1)) - 1);   /* neighbors > a */
        uint32_t bb = na;
        while (bb) {
            int b = __builtin_ctz(bb);
            bb &= bb - 1;
            uint32_t m = na & adj[b] & ~((1u << (b + 1)) - 1);
            uint32_t cc = m;
            while (cc) {
                int c = __builtin_ctz(cc);
                cc &= cc - 1;
                if (m & adj[c] & ~((1u << (c + 1)) - 1)) return 1;
            }
        }
    }
    return 0;
}

static int has_i5(const uint32_t adj[], int n) {
    uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
    uint32_t cadj[32];
    for (int i = 0; i < n; i++) cadj[i] = full & ~adj[i] & ~(1u << i);
    for (int a = 0; a < n; a++) {
        uint32_t m1 = cadj[a] & ~((1u << (a + 1)) - 1);
        uint32_t bb = m1;
        while (bb) {
            int b = __builtin_ctz(bb);
            bb &= bb - 1;
            uint32_t m2 = m1 & cadj[b] & ~((1u << (b + 1)) - 1);
            uint32_t cc = m2;
            while (cc) {
                int c = __builtin_ctz(cc);
                cc &= cc - 1;
                uint32_t m3 = m2 & cadj[c] & ~((1u << (c + 1)) - 1);
                uint32_t dd = m3;
                while (dd) {
                    int dv = __builtin_ctz(dd);
                    dd &= dd - 1;
                    if (m3 & cadj[dv] & ~((1u << (dv + 1)) - 1)) return 1;
                }
            }
        }
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "usage: verify n e file.g6\n");
        return 2;
    }
    int want_n = atoi(argv[1]), want_e = atoi(argv[2]);
    FILE *f = fopen(argv[3], "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", argv[3]); return 2; }
    char buf[512];
    long checked = 0, fails = 0;
    while (fgets(buf, sizeof buf, f)) {
        size_t L = strlen(buf);
        while (L && (buf[L-1] == '\n' || buf[L-1] == '\r')) buf[--L] = 0;
        if (!L) continue;
        checked++;
        uint32_t adj[32];
        int n;
        if (parse_g6(buf, L, adj, &n) != 0) {
            fprintf(stderr, "line %ld: unparsable\n", checked);
            fails++;
            continue;
        }
        int bad = 0;
        if (n != want_n) { fprintf(stderr, "line %ld: n=%d\n", checked, n); bad = 1; }
        long degsum = 0;
        for (int i = 0; i < n; i++) degsum += __builtin_popcount(adj[i]);
        if (degsum / 2 != want_e) {
            fprintf(stderr, "line %ld: e=%ld\n", checked, degsum / 2);
            bad = 1;
        }
        if (has_k4(adj, n)) { fprintf(stderr, "line %ld: K4\n", checked); bad = 1; }
        if (has_i5(adj, n)) { fprintf(stderr, "line %ld: I5\n", checked); bad = 1; }
        if (bad) fails++;
    }
    fclose(f);
    printf("checked %ld graphs, failures %ld\n", checked, fails);
    return fails ? 1 : 0;
}
