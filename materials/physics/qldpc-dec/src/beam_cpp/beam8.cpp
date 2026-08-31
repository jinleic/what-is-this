// beam8.cpp — C++ port of physics/qldpc-dec/src/qldpc_dec/beam_search.py
// (BeamSearchDecoder, "beam8" arm of Gate B, arXiv:2512.07057 Algorithm 3).
//
// Equivalence contract with the pinned Python arm (see README.md):
//   * Input: stim textual DEM containing ONLY straight-line "error(p) D.. L.."
//     lines (the canonical str(dem) of every committed BB_144 DEM in this repo
//     has exactly 8784 such lines and no REPEAT/shift instructions — verified).
//     Raw convention: column j = j-th error instruction; per-column D/L ids
//     ascending (matches scipy sorted_indices used by dem_matrices.py).
//   * lam_j = log((1-p)/max(p,1e-12)). np.log may differ from std::log by <=1
//     ulp, so the harness exports the exact float64 lam array the Python
//     decoder uses (--lam); the C++ log is only used when no export is given
//     (and --check-lam reports any bit differences).
//   * Decode pipeline is a bit-exact replication of BeamSearchDecoder.decode
//     with the pinned parameters (beam 8 / 30 / 20 / 10 / num_results=1),
//     including: numpy float64 pairwise-summation order (np_pairwise_sum),
//     stable-sort beam pruning with creation order as tie-break, first-min
//     scans, (path, val)-order expansion, and the first converged+valid
//     early return. No RNG (the Python decoder has none; seeds only feed the
//     Stim sampler, which runs on the harness side).
//
// Build: clang++ -O3 -std=c++17 -Wall -Wextra beam8.cpp -o beam8_cpp
// Run:   ./beam8_cpp --dem dem.txt [--lam lam.npy] --shots shots.bin
//                   --out preds.bin [--limit N] [--repeats K] [...]
//        (exact commands in README.md / benchmark.sh)

#include <algorithm>
#include <atomic>
#include <array>
#include <cassert>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

namespace {

constexpr uint32_t kOutMagic = 0x31503842u;  // "B8P1" little-endian

[[noreturn]] void fail(const std::string& msg) { throw std::runtime_error(msg); }

// ---------------------------------------------------------------------------
// numpy float64 pairwise-sum replication (numpy/_core/src/umath/loops_utils.h,
// pairwise_sum_double; PW_BLOCKSIZE == 128). np.add.reduce / np.sum / np.mean
// on contiguous float64 all funnel through this exact accumulation order, and
// the Python decoder's mu-sums, score sums and msg means all hit it (lengths
// 16..35 and ~8784 here, so both the unrolled-8 and the split branches run).
// ---------------------------------------------------------------------------
double np_sum(const double* a, int64_t n) {
    if (n < 8) {
        double res = 0.0;
        for (int64_t i = 0; i < n; ++i) res += a[i];
        return res;
    }
    if (n <= 128) {
        double r[8];
        for (int j = 0; j < 8; ++j) r[j] = a[j];
        int64_t i = 8;
        const int64_t lim = n - (n % 8);
        for (; i < lim; i += 8) {
            for (int j = 0; j < 8; ++j) r[j] += a[i + j];
        }
        double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; ++i) res += a[i];
        return res;
    }
    int64_t n2 = n / 2;
    n2 -= n2 % 8;
    return np_sum(a, n2) + np_sum(a + n2, n - n2);
}

// ---------------------------------------------------------------------------
// DEM parsing (straight-line error lines only; contract above)
// ---------------------------------------------------------------------------

struct Dem {
    int32_t ndet = 0, nobs = 0, nerr = 0;
    std::vector<std::vector<int32_t>> col_dets, col_obs;  // ascending per col
    std::vector<double> p;
};

Dem load_dem_text(const std::string& path) {
    std::ifstream in(path);
    if (!in) fail("cannot open DEM file: " + path);
    Dem d;
    std::string line;
    int lineno = 0;
    while (std::getline(in, line)) {
        ++lineno;
        const size_t b = line.find_first_not_of(" \t\r\n");
        if (b == std::string::npos) continue;
        if (line.compare(b, 5, "error") != 0)
            fail("dem line " + std::to_string(lineno) +
                 ": unsupported instruction (port contract: straight-line "
                 "'error(...)' lines only): " + line.substr(b, 40));
        const size_t po = line.find('(', b);
        if (po == std::string::npos) fail("dem line " + std::to_string(lineno) + ": missing '('");
        const char* pbegin = line.c_str() + po + 1;
        char* pend = nullptr;
        const double p = std::strtod(pbegin, &pend);
        if (pend == pbegin) fail("dem line " + std::to_string(lineno) + ": bad probability");
        d.p.push_back(p);
        std::vector<int32_t> dets, obs;
        size_t i = (size_t)(pend - line.c_str());
        while (i < line.size() && line[i] != ')') ++i;   // skip past probability ')'
        if (i >= line.size()) fail("dem line " + std::to_string(lineno) + ": missing ')'");
        ++i;  // consume ')'
        while (i < line.size()) {
            if (std::isspace((unsigned char)line[i])) { ++i; continue; }
            const char c = line[i];
            if (c != 'D' && c != 'L')
                fail("dem line " + std::to_string(lineno) + ": bad target '" + line.substr(i, 12) + "'");
            const int32_t id = (int32_t)std::strtol(line.c_str() + i + 1, nullptr, 10);
            if (c == 'D') { dets.push_back(id); if (id + 1 > d.ndet) d.ndet = id + 1; }
            else          { obs.push_back(id);  if (id + 1 > d.nobs) d.nobs = id + 1; }
            while (i < line.size() && !std::isspace((unsigned char)line[i])) ++i;
        }
        std::sort(dets.begin(), dets.end());
        std::sort(obs.begin(), obs.end());
        if (std::adjacent_find(dets.begin(), dets.end()) != dets.end() ||
            std::adjacent_find(obs.begin(), obs.end()) != obs.end())
            fail("dem line " + std::to_string(lineno) + ": duplicate target (would change scipy semantics)");
        d.col_dets.push_back(std::move(dets));
        d.col_obs.push_back(std::move(obs));
    }
    d.nerr = (int32_t)d.p.size();
    if (d.nerr == 0) fail("dem has no error instructions");
    return d;
}

// ---------------------------------------------------------------------------
// minimal .npy reader (float64 '<f8', 1-D) for the harness-exported lam
// ---------------------------------------------------------------------------
std::vector<double> load_lam_npy(const std::string& path, int32_t expect_n) {
    std::ifstream in(path, std::ios::binary);
    if (!in) fail("cannot open lam npy: " + path);
    unsigned char hdr[10];
    in.read((char*)hdr, 10);
    if (in.gcount() != 10 || std::memcmp(hdr, "\x93NUMPY", 6) != 0)
        fail("lam npy: bad magic");
    const int major = hdr[6];
    uint32_t hlen32 = 0;
    if (major == 1) {
        std::memcpy(&hlen32, hdr + 8, 2);  // u16 LE at bytes 8..10
    } else {
        fail("lam npy: unsupported npy version (want v1)");
    }
    std::string dict((size_t)hlen32, '\0');
    in.read(dict.data(), (std::streamsize)hlen32);
    if ((int)in.gcount() != (int)hlen32) fail("lam npy: short header");
    if (dict.find("'<f8'") == std::string::npos && dict.find("\"<f8\"") == std::string::npos &&
        dict.find("'<d'") == std::string::npos)
        fail("lam npy: expected '<f8' dtype, header: " + dict);
    // parse shape: 1-D float array of length expect_n
    const size_t sp = dict.find("'shape'");
    if (sp == std::string::npos) fail("lam npy: no shape");
    const size_t op = dict.find('(', sp);
    const size_t cp = dict.find(')', sp);
    if (op == std::string::npos || cp == std::string::npos || cp < op)
        fail("lam npy: bad shape");
    const std::string shape = dict.substr(op + 1, cp - op - 1);
    size_t fs = shape.find_first_not_of(' ');
    size_t fe = shape.find_last_not_of(' ');
    const std::string first_dim = (fs == std::string::npos) ? "" : shape.substr(fs, fe - fs + 1);
    const size_t comma = first_dim.find(',');
    if (comma == std::string::npos || first_dim.find_first_not_of("0123456789 ") == std::string::npos)
        fail("lam npy: expected 1-D shape like (N,), got '" + shape + "'");
    const int64_t nfile = std::stoll(first_dim.substr(0, comma));
    if (nfile != expect_n)
        fail("lam npy: shape " + std::to_string(nfile) + " != expected " + std::to_string(expect_n));
    const int64_t n = expect_n;
    std::vector<double> v((size_t)n);
    in.read((char*)v.data(), (std::streamsize)(n * 8));
    if ((int64_t)in.gcount() != n * 8) fail("lam npy: truncated payload");
    return v;
}

// ---------------------------------------------------------------------------
// CSR pattern matrix built from per-column row lists (counting sort keeps
// ascending column order per row == scipy csc->csr sorted_indices)
// ---------------------------------------------------------------------------
struct CSR {
    std::vector<int32_t> indptr, indices;
    int32_t rows = 0, cols = 0;
};

CSR rows_from_cols(const std::vector<std::vector<int32_t>>& cols, int32_t nrows, int32_t ncols) {
    CSR m;
    m.rows = nrows;
    m.cols = ncols;
    m.indptr.assign((size_t)nrows + 1, 0);
    for (const auto& c : cols)
        for (int32_t r : c) ++m.indptr[(size_t)r + 1];
    for (int32_t r = 0; r < nrows; ++r) m.indptr[(size_t)r + 1] += m.indptr[(size_t)r];
    m.indices.resize(m.indptr.back());
    std::vector<int32_t> cur(m.indptr.begin(), m.indptr.end() - 1);
    for (int32_t j = 0; j < ncols; ++j)
        for (int32_t r : cols[(size_t)j]) m.indices[(size_t)cur[(size_t)r]++] = j;
    return m;
}

// row parity of a 0/1 vector e against a target (active-row variant below)
bool csr_rows_match(const CSR& H, const std::vector<uint8_t>& e, const std::vector<uint8_t>& target) {
    for (int32_t i = 0; i < H.rows; ++i) {
        int par = 0;
        for (int32_t k = H.indptr[(size_t)i]; k < H.indptr[(size_t)i + 1]; ++k)
            par ^= e[(size_t)H.indices[(size_t)k]];
        if ((uint8_t)par != target[(size_t)i]) return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// little-endian binary IO
// ---------------------------------------------------------------------------
uint64_t get_u64(const uint8_t* p) { uint64_t v; std::memcpy(&v, p, 8); return v; }
void put_u64(uint8_t* p, uint64_t v) { std::memcpy(p, &v, 8); }
uint32_t get_u32(const uint8_t* p) { uint32_t v; std::memcpy(&v, p, 4); return v; }
void put_u32(uint8_t* p, uint32_t v) { std::memcpy(p, &v, 4); }

// ---------------------------------------------------------------------------
// masked min-sum BP + beam search (the port itself)
// ---------------------------------------------------------------------------

struct ShotStat {
    bool seed_conv = false;
    bool returned_early = false;     // first converged+valid child -> return
    bool converged_invalid = false;
    bool fallback = false;           // best-effort branch taken
    int32_t bp_runs = 0;
    int64_t bp_iters = 0;
    int32_t rounds = 0;
};

struct Engine {
    int32_t ndet = 0, nerr = 0, nobs = 0;
    std::vector<double> lam;                 // nerr
    CSR H;                                   // ndet x nerr pattern
    std::vector<std::vector<int32_t>> H_col; // CSC of H: dets per col (ascending)
    CSR A;                                   // nobs x nerr pattern
    // pinned parameters (CLI-overridable; see README)
    int32_t beam_width = 8, initial_iters = 30, iters_per_round = 20;
    int32_t max_rounds = 10;
};

// One beam path. Mirrors the path dict in Python decode().
struct Path {
    std::vector<double> msgs;             // nerr (warm-start E_j->i(0))
    std::vector<std::pair<int32_t, uint8_t>> masked;  // (col, val), insertion order
    int32_t next = -1;
    double score = 0.0;
};

struct Bp {
    // immutable per engine
    const Engine* eng = nullptr;

    // per-construction (masked set dependent)
    std::vector<uint8_t> active;   // nerr
    std::vector<int32_t> act;      // list of active cols (ascending)
    std::vector<uint8_t> s;        // ndet, effective syndrome
    std::vector<int32_t> pair_det, pair_err;
    int32_t npairs = 0;
    std::vector<std::vector<int32_t>> det_pairs, err_pairs;  // compression

    // run scratch
    std::vector<double> nu, mu, sum_llr;
    std::vector<uint8_t> e_hat;
    std::vector<double> min1, min2, sg;
    std::vector<int32_t> argmin, deg;
    void build(const Engine& eng_in, const std::vector<std::pair<int32_t, uint8_t>>& masked,
               const std::vector<uint8_t>& syndrome);

    // returns converged; fills e_hat, sum_llr, nu; outputs iters_run
    bool run(const double* edge_msgs, int32_t max_iters, int32_t& iters_run);
};

void Bp::build(const Engine& e, const std::vector<std::pair<int32_t, uint8_t>>& masked,
                  const std::vector<uint8_t>& syndrome) {
    eng = &e;
    const int32_t ndet = e.ndet, nerr = e.nerr;
    active.assign((size_t)nerr, 1);
    for (const auto& [j, v] : masked) active[(size_t)j] = 0;
    act.clear();
    for (int32_t j = 0; j < nerr; ++j)
        if (active[(size_t)j]) act.push_back(j);

    // effective syndrome: masked-1 nodes flip their incident detector bits
    s = syndrome;
    for (const auto& [j, v] : masked)
        if (v)
            for (int32_t d : e.H_col[(size_t)j]) s[(size_t)d] ^= 1;

    // pair compression in CSR row order (== err_pair_idx_of(i) order)
    pair_det.clear();
    pair_err.clear();
    det_pairs.assign((size_t)ndet, {});
    for (int32_t i = 0; i < ndet; ++i)
        for (int32_t k = e.H.indptr[(size_t)i]; k < e.H.indptr[(size_t)i + 1]; ++k) {
            const int32_t j = e.H.indices[(size_t)k];
            if (!active[(size_t)j]) continue;
            det_pairs[(size_t)i].push_back((int32_t)pair_det.size());
            pair_det.push_back(i);
            pair_err.push_back(j);
        }
    npairs = (int32_t)pair_det.size();
    err_pairs.assign((size_t)nerr, {});
    for (int32_t idx = 0; idx < npairs; ++idx)
        err_pairs[(size_t)pair_err[(size_t)idx]].push_back(idx);

    nu.assign((size_t)npairs, 0.0);
    mu.assign((size_t)npairs, 0.0);
    sum_llr.assign((size_t)nerr, 0.0);
    e_hat.assign((size_t)nerr, 0);
    min1.assign((size_t)ndet, 0.0);
    min2.assign((size_t)ndet, 0.0);
    sg.assign((size_t)ndet, 0.0);
    argmin.assign((size_t)ndet, 0);
    deg.assign((size_t)ndet, 0);
}

// np.sum(mu[ks]) / nu[ks].mean() replication: numpy fancy indexing COPIES the
// selected elements into a fresh contiguous array, then pairwise-sums THAT.
// NOTE: an error's edge pairs are scattered in pair order (grouped by
// detector), so the sum must run over the gathered copy, never a raw slice.
double np_sum_pairs(const std::vector<double>& mu, const std::vector<int32_t>& ks) {
    static thread_local std::vector<double> g;
    g.resize(ks.size());
    for (size_t q = 0; q < ks.size(); ++q) g[q] = mu[(size_t)ks[q]];
    return np_sum(g.data(), (int64_t)ks.size());
}

static inline double lam_of(const Engine& e, int32_t j) { return e.lam[(size_t)j]; }

// convergence: parity over ACTIVE columns equals the effective syndrome s,
// full-length vector compare (== np.array_equal(reduced, self.s))
bool parity_active_converged(const Engine& e, const Bp& bp) {
    // reduced = (H_active @ e_hat[act]) % 2 over ALL ndet rows
    // (H_active columns are only active ones, so iterate columns)
    static thread_local std::vector<uint8_t> red;
    red.assign((size_t)e.ndet, 0);
    for (int32_t idx = 0; idx < bp.npairs; ++idx) {
        if (bp.e_hat[(size_t)bp.pair_err[(size_t)idx]])
            red[(size_t)bp.pair_det[(size_t)idx]] ^= 1;
    }
    for (int32_t i = 0; i < e.ndet; ++i)
        if (red[(size_t)i] != bp.s[(size_t)i]) return false;
    return true;
}
// One masked-BP run. Bit-exact vs _MaskedBP.run.
bool Bp::run(const double* edge_msgs, int32_t max_iters, int32_t& iters_run) {
    const Engine& e = *eng;
    const int32_t ndet = e.ndet, nerr = e.nerr;

    if (edge_msgs)
        for (int32_t idx = 0; idx < npairs; ++idx) nu[(size_t)idx] = edge_msgs[pair_err[(size_t)idx]];
    else
        for (int32_t idx = 0; idx < npairs; ++idx) nu[(size_t)idx] = 0.0;

    std::fill(sum_llr.begin(), sum_llr.end(), 0.0);
    std::fill(e_hat.begin(), e_hat.end(), (uint8_t)0);
    bool converged = false;
    iters_run = 0;

    for (int32_t t = 1; t <= max_iters; ++t) {
        iters_run = t;
        // ---- check -> error (A1) ----
        for (int32_t i = 0; i < ndet; ++i) {
            const auto& ks = det_pairs[(size_t)i];
            const int32_t dk = (int32_t)ks.size();
            deg[(size_t)i] = dk;
            double mn1 = std::numeric_limits<double>::infinity();
            double mn2 = std::numeric_limits<double>::infinity();
            int32_t am = -1;
            double sp = 1.0;
            for (int32_t q = 0; q < dk; ++q) {
                const double v = nu[(size_t)ks[(size_t)q]];
                const double av = std::fabs(v);
                if (q == 0 || av < mn1) {             // first-lowest, strict <  : stable argmin
                    if (q != 0) mn2 = mn1;
                    else if (dk > 1) mn2 = std::numeric_limits<double>::infinity();
                    mn1 = av;
                    am = ks[(size_t)q];
                } else if (q == 1 || av < mn2) {
                    mn2 = av;
                }
                const double sgn = v >= 0.0 ? 1.0 : -1.0;
                sp *= sgn;
            }
            if (dk == 0) {
                min1[(size_t)i] = std::numeric_limits<double>::infinity();
                min2[(size_t)i] = std::numeric_limits<double>::infinity();
                argmin[(size_t)i] = -1;
                sg[(size_t)i] = 1.0;
                continue;
            }
            min1[(size_t)i] = mn1;
            min2[(size_t)i] = mn2;
            argmin[(size_t)i] = am;
            sg[(size_t)i] = sp;
        }
        // mu update (A1 application), in pair order (== numpy idx loop order)
        for (int32_t idx = 0; idx < npairs; ++idx) {
            const int32_t i = pair_det[(size_t)idx];
            if (deg[(size_t)i] == 0) { mu[(size_t)idx] = 0.0; continue; }
            const double m = (idx == argmin[(size_t)i]) ? min2[(size_t)i] : min1[(size_t)i];
            const double own = nu[(size_t)idx];
            const double own_sign = own >= 0.0 ? 1.0 : -1.0;
            const double kappa = sg[(size_t)i] / own_sign;
            mu[(size_t)idx] = kappa * (s[(size_t)i] == 0 ? 1.0 : -1.0) * m;
        }
        // ---- error -> check (A2) ----
        for (int32_t j : act) {
            const auto& ks = err_pairs[(size_t)j];
            const double tot = np_sum_pairs(mu, ks);
            const int64_t n = (int64_t)ks.size();
            for (int64_t q = 0; q < n; ++q)
                nu[(size_t)ks[(size_t)q]] = lam_of(e, j) + tot - mu[(size_t)ks[(size_t)q]];
        }
        // ---- posteriors (A3) ----
        for (int32_t j : act) {
            const auto& ks = err_pairs[(size_t)j];
            const double post = lam_of(e, j) + np_sum_pairs(mu, ks);
            sum_llr[(size_t)j] += post;
            e_hat[(size_t)j] = post > 0.0 ? 0 : 1;
        }
        // ---- convergence on the active submatrix ----
        if (parity_active_converged(e, *this)) { converged = true; break; }
    }
    return converged;
}

// ---------------------------------------------------------------------------
// BeamSearchDecoder.decode port (bit-exact; same return points)
// ---------------------------------------------------------------------------

ShotStat decode_shot(Engine& e, const std::vector<uint8_t>& syndrome,
                     std::vector<uint8_t>& obs_out) {
    ShotStat st{};
    const int32_t nerr = e.nerr, ndet = e.ndet;
    obs_out.assign((size_t)e.nobs, 0);

    std::vector<uint8_t> syn64((size_t)ndet);
    for (int32_t i = 0; i < ndet; ++i) syn64[(size_t)i] = (uint8_t)(syndrome[(size_t)i] ? 1 : 0);

    // results: (weight, e_full) list; num_results=1 -> return first valid
    double best_wt = 0.0;
    std::vector<uint8_t> best_e;
    bool have_result = false;

    // ---------- seed: standard BP (no masked nodes) ----------
    Bp seed;
    seed.build(e, {}, syn64);
    int32_t seed_iters = 0;
    ++st.bp_runs;
    const bool conv = seed.run(nullptr, e.initial_iters, seed_iters);
    st.bp_iters += seed_iters;
    if (conv) {
        st.seed_conv = true;
        // return e_hat directly (masked empty -> e_full == e_hat)
        // observable prediction: (A @ e_hat) % 2
        for (int32_t o = 0; o < e.nobs; ++o) {
            int par = 0;
            for (int32_t k = e.A.indptr[(size_t)o]; k < e.A.indptr[(size_t)o + 1]; ++k)
                par ^= seed.e_hat[(size_t)e.A.indices[(size_t)k]];
            obs_out[(size_t)o] = (uint8_t)par;
        }
        return st;
    }

    // seed warm-start messages: per active err j, mean of nu over its pairs
    std::vector<double> seed_msgs((size_t)nerr, 0.0);
    for (int32_t j : seed.act) {
        const auto& ks = seed.err_pairs[(size_t)j];
        double acc = 0.0;
        // np.mean == np.add.reduce / n  (pairwise order then one division)
        acc = np_sum_pairs(seed.nu, ks);
        // pairwise-sum replication: np.mean computes sum via pairwise then /n
        seed_msgs[(size_t)j] = acc / (double)ks.size();
    }
    // nxt = int(np.argmin(np.abs(sum_llr)))  -- FIRST minimum
    int32_t nxt = 0;
    {
        double best = std::fabs(seed.sum_llr[0]);
        for (int32_t j = 1; j < nerr; ++j) {
            const double v = std::fabs(seed.sum_llr[(size_t)j]);
            if (v < best) { best = v; nxt = j; }
        }
    }
    Path root;
    root.msgs = std::move(seed_msgs);
    root.next = nxt;
    root.score = 0.0;
    std::vector<Path> paths;
    paths.push_back(std::move(root));

    const int32_t max_rounds = e.max_rounds;
    const int32_t beam_width = e.beam_width;

    for (int32_t rnd = 0; rnd < max_rounds; ++rnd) {
        st.rounds = rnd + 1;
        std::vector<Path> next_set;
        for (const Path& path : paths) {
            for (uint8_t val = 0; val <= 1; ++val) {
                std::vector<std::pair<int32_t, uint8_t>> masked = path.masked;
                masked.emplace_back(path.next, val);
                Bp eng2;
                eng2.build(e, masked, syn64);
                int32_t it2 = 0;
                ++st.bp_runs;
                const bool conv2 = eng2.run(path.msgs.data(), e.iters_per_round, it2);
                st.bp_iters += it2;
                if (conv2) {
                    // e_full = e_h with masked entries overwritten
                    std::vector<uint8_t> e_full((size_t)nerr);
                    for (int32_t j = 0; j < nerr; ++j) e_full[(size_t)j] = eng2.e_hat[(size_t)j];
                    for (const auto& [j, v] : masked) e_full[(size_t)j] = v;
                    // full parity check vs ORIGINAL syndrome
                    bool ok = true;
                    for (int32_t i = 0; i < ndet && ok; ++i) {
                        int par = 0;
                        for (int32_t k = e.H.indptr[(size_t)i]; k < e.H.indptr[(size_t)i + 1]; ++k)
                            par ^= e_full[(size_t)e.H.indices[(size_t)k]];
                        if ((uint8_t)par != syndrome[(size_t)i]) ok = false;
                    }
                    if (ok) {
                        // wt = float(np.dot(e_full, self.weight)); results.append
                        double wt = 0.0;
                        for (int32_t j = 0; j < nerr; ++j)
                            if (e_full[(size_t)j]) wt += std::fabs(e.lam[(size_t)j]);
                        if (!have_result || wt < best_wt) {
                            best_wt = wt;
                            best_e = e_full;
                        }
                        have_result = true;
                        // num_results == 1 -> return immediately
                        st.returned_early = true;
                        // obs_out from best_e
                        for (int32_t o = 0; o < e.nobs; ++o) {
                            int par = 0;
                            for (int32_t k = e.A.indptr[(size_t)o]; k < e.A.indptr[(size_t)o + 1]; ++k)
                                par ^= best_e[(size_t)e.A.indices[(size_t)k]];
                            obs_out[(size_t)o] = (uint8_t)par;
                        }
                        return st;
                    }
                    st.converged_invalid = true;
                }
                // build child path (scoring / pruning) even when converged-invalid
                std::vector<int32_t> unmasked;
                unmasked.reserve((size_t)nerr - masked.size());
                {
                    std::vector<uint8_t> fl((size_t)nerr, 0);
                    for (const auto& [j, v] : masked) (void)v, fl[(size_t)j] = 1;
                    for (int32_t j = 0; j < nerr; ++j)
                        if (!fl[(size_t)j]) unmasked.push_back(j);
                }
                if (unmasked.empty()) continue;
                // score = float(np.sum(np.abs(s_lr[unmasked])) / max(it2, 1))
                //       = pairwise_sum(|s_lr[unmasked]|) / max(it2,1) -- GAP!
                // np.sum operates on the GATHERED array (fancy indexing copies),
                // so the pairwise tree runs over the compacted sequence.
                double acc = 0.0;
                {
                    // materialize |s_lr[unmasked]| in unmasked order, then pairwise
                    static thread_local std::vector<double> g;
                    g.resize(unmasked.size());
                    for (size_t q = 0; q < unmasked.size(); ++q)
                        g[q] = std::fabs(eng2.sum_llr[(size_t)unmasked[q]]);
                    acc = np_sum(g.data(), (int64_t)g.size());
                }
                const double score = acc / (double)std::max(it2, 1);
                // nextj = first unmasked j minimizing |s_lr[j]| (argmin-first)
                int32_t nextj = -1;
                {
                    double best = std::numeric_limits<double>::infinity();
                    for (int32_t j : unmasked) {
                        const double v = std::fabs(eng2.sum_llr[(size_t)j]);
                        if (v < best) { best = v; nextj = j; }
                    }
                }
                Path child;
                child.msgs.assign((size_t)nerr, 0.0);
                for (int32_t j : eng2.act) {
                    const auto& ks = eng2.err_pairs[(size_t)j];
                    child.msgs[(size_t)j] = np_sum_pairs(eng2.nu, ks) / (double)ks.size();
                }
                child.masked = std::move(masked);
                child.next = nextj;
                child.score = score;
                next_set.push_back(std::move(child));
            }
        }
        if (next_set.empty()) break;
        // next_set.sort(key=score, reverse=True) then top beam_width.
        // Python sort is stable -> equal keys keep creation order; replicate
        // with a stable sort by descending score.
        std::stable_sort(next_set.begin(), next_set.end(),
                         [](const Path& a, const Path& b) { return a.score > b.score; });
        if ((int32_t)next_set.size() > beam_width) next_set.resize((size_t)beam_width);
        paths = std::move(next_set);
        // if all(p["next"] < 0 for p in paths): break
        bool all_neg = true;
        for (const Path& p : paths)
            if (p.next >= 0) { all_neg = false; break; }
        if (all_neg) break;
    }
    if (have_result) {
        st.fallback = false;
        // min(results) by weight -- with num_results=1 we already returned,
        // but keep the branch for num_results>1 future-proofing
        for (int32_t o = 0; o < e.nobs; ++o) {
            int par = 0;
            for (int32_t k = e.A.indptr[(size_t)o]; k < e.A.indptr[(size_t)o + 1]; ++k)
                par ^= best_e[(size_t)e.A.indices[(size_t)k]];
            obs_out[(size_t)o] = (uint8_t)par;
        }
        return st;
    }
    // no valid solution: best-effort (paths[0].masked)
    st.fallback = true;
    {
        std::vector<uint8_t> e_full((size_t)nerr, 0);
        if (!paths.empty())
            for (const auto& [j, v] : paths[0].masked) (void)v, e_full[(size_t)j] = v;
        // note: masked vals are only informative on '1'; v==0 sets 0 anyway
        for (int32_t o = 0; o < e.nobs; ++o) {
            int par = 0;
            for (int32_t k = e.A.indptr[(size_t)o]; k < e.A.indptr[(size_t)o + 1]; ++k)
                par ^= e_full[(size_t)e.A.indices[(size_t)k]];
            obs_out[(size_t)o] = (uint8_t)par;
        }
    }
    return st;
}

}  // namespace

// ---------------------------------------------------------------------------
// CLI + shot IO + batch loop
// ---------------------------------------------------------------------------

struct Args {
    std::string dem, lam, shots, out;
    int64_t limit = -1;        // decode first N shots (-1: all)
    int32_t repeats = 1;       // repeat whole batch (timing); decode count scaled
    int32_t threads = 1;       // harness-level parallelism: split shoot range
    int32_t beam_width = 8, initial_iters = 30, iters_per_round = 20;
    int32_t max_rounds = 10;   // num_results pinned to 1 (companion arm)
    bool check_lam = false;
    bool verbose = false;
};

Args parse_args(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        const std::string s = argv[i];
        auto need_val = [&](const char* name) -> std::string {
            if (i + 1 >= argc) fail(std::string("missing value for ") + name);
            return argv[++i];
        };
        auto eqval = [&](const char* name) -> std::string {
            const size_t eq = s.find('=');
            if (eq == std::string::npos) fail(std::string("expected ") + name + "=VALUE");
            return s.substr(eq + 1);
        };
        if (s == "--dem") a.dem = need_val("--dem");
        else if (s.rfind("--dem=", 0) == 0) a.dem = eqval("--dem");
        else if (s == "--lam") a.lam = need_val("--lam");
        else if (s.rfind("--lam=", 0) == 0) a.lam = eqval("--lam");
        else if (s == "--shots") a.shots = need_val("--shots");
        else if (s.rfind("--shots=", 0) == 0) a.shots = eqval("--shots");
        else if (s == "--out") a.out = need_val("--out");
        else if (s.rfind("--out=", 0) == 0) a.out = eqval("--out");
        else if (s == "--limit") a.limit = std::stoll(need_val("--limit"));
        else if (s.rfind("--limit=", 0) == 0) a.limit = std::stoll(eqval("--limit"));
        else if (s == "--repeats") a.repeats = std::stoi(need_val("--repeats"));
        else if (s.rfind("--repeats=", 0) == 0) a.repeats = std::stoi(eqval("--repeats"));
        else if (s == "--threads") a.threads = std::stoi(need_val("--threads"));
        else if (s.rfind("--threads=", 0) == 0) a.threads = std::stoi(eqval("--threads"));
        else if (s == "--check-lam") a.check_lam = true;
        else if (s == "--verbose") a.verbose = true;
        else if (s.rfind("--beam-width=", 0) == 0) a.beam_width = std::stoi(eqval("--beam-width"));
        else if (s.rfind("--initial-iters=", 0) == 0) a.initial_iters = std::stoi(eqval("--initial-iters"));
        else if (s.rfind("--iters-per-round=", 0) == 0) a.iters_per_round = std::stoi(eqval("--iters-per-round"));
        else if (s.rfind("--max-rounds=", 0) == 0) a.max_rounds = std::stoi(eqval("--max-rounds"));
        else fail("unknown argument: " + s);
    }
    if (a.dem.empty() || a.shots.empty() || a.out.empty())
        fail("usage: beam8_cpp --dem DEM --shots SHOTS --out OUT [--lam lam.npy] "
             "[--limit N] [--repeats K] [--threads T] [--check-lam] [--verbose]");
    return a;
}

int main(int argc, char** argv) {
    try {
        const Args args = parse_args(argc, argv);

        // ---- DEM ----
        Dem dem = load_dem_text(args.dem);
        Engine eng;
        eng.ndet = dem.ndet;
        eng.nerr = dem.nerr;
        eng.nobs = dem.nobs;

        // ---- lam ----
        if (!args.lam.empty()) {
            eng.lam = load_lam_npy(args.lam, dem.nerr);

        } else {
            eng.lam.resize((size_t)dem.nerr);
            for (int32_t j = 0; j < dem.nerr; ++j) {
                const double pj = dem.p[(size_t)j] < 1e-12 ? 1e-12 : dem.p[(size_t)j];
                eng.lam[(size_t)j] = std::log((1.0 - dem.p[(size_t)j]) / pj);
            }
        }

        // ---- matrices ----
        eng.H = rows_from_cols(dem.col_dets, dem.ndet, dem.nerr);
        eng.A = rows_from_cols(dem.col_obs, dem.nobs, dem.nerr);
        eng.H_col = dem.col_dets;

        // ---- shots ----
        std::ifstream sf(args.shots, std::ios::binary);
        if (!sf) fail("cannot open shots file: " + args.shots);
        uint64_t hdr[3];
        sf.read((char*)hdr, 24);
        if (sf.gcount() != 24) fail("shots file too short");
        const uint64_t ndet_f = get_u64((const uint8_t*)hdr + 0);
        const uint64_t nobs_f = get_u64((const uint8_t*)hdr + 8);
        const uint64_t nshot_f = get_u64((const uint8_t*)hdr + 16);
        if (ndet_f != (uint64_t)eng.ndet)
            fail("shots ndet " + std::to_string(ndet_f) + " != DEM ndet " + std::to_string(eng.ndet));
        const int64_t nshot = (int64_t)nshot_f;
        const int64_t limit = args.limit >= 0 ? std::min<int64_t>(args.limit, nshot) : nshot;
        const size_t det_bytes = (size_t)nshot * (size_t)((ndet_f + 63) / 64 * 8);  // u64-LE rows
        const size_t obs_bytes = (size_t)nshot * (size_t)((nobs_f + 63) / 64 * 8);
        std::vector<uint8_t> det_bits(det_bytes), obs_bits(obs_bytes);
        sf.read((char*)det_bits.data(), (std::streamsize)det_bytes);
        if ((size_t)sf.gcount() != det_bytes) fail("shots file truncated (det)");
        sf.read((char*)obs_bits.data(), (std::streamsize)obs_bytes);
        if ((size_t)sf.gcount() != obs_bytes) fail("shots file truncated (obs)");

        if (args.verbose)
            std::fprintf(stderr,
                         "dem: dets=%d errs=%d obs=%d | shots=%lld limit=%lld repeats=%d threads=%d\n",
                         eng.ndet, eng.nerr, eng.nobs, (long long)nshot, (long long)limit,
                         args.repeats, args.threads);

        std::vector<ShotStat> stats((size_t)limit);
        std::mutex stats_mu;
        const auto t0 = std::chrono::steady_clock::now();
        std::vector<uint8_t> preds((size_t)limit * (size_t)((nobs_f + 63) / 64 * 8), 0);

        for (int32_t rep = 0; rep < args.repeats; ++rep) {
            std::vector<std::thread> pool;
            std::atomic<int64_t> next_shot{0};
            auto worker = [&](int) {
                Engine local = eng;  // copy-on-write semantics: each thread owns one
                for (;;) {
                    const int64_t k = next_shot.fetch_add(1);
                    if (k >= limit) break;
                    std::vector<uint8_t> syn((size_t)eng.ndet);
                    const uint8_t* row = det_bits.data() + (size_t)k * ((ndet_f + 63) / 64 * 8);
                    for (uint64_t d = 0; d < ndet_f; ++d) {
                        // stim u64-LE packing: detector d = bit (d%64) of word d/64;
                        // word w occupies bytes 8w..8w+7 little-endian.
                        const uint64_t byte_idx = (d >> 6) * 8 + ((d & 63) >> 3);
                        syn[(size_t)d] = (row[byte_idx] >> (d & 7)) & 1u;
                    }
                    std::vector<uint8_t> ob((size_t)nobs_f, 0);  // unpacked per-obs
                    ShotStat st = decode_shot(local, syn, ob);
                    // decode_shot fills obs_out UNPACKED: one byte per observable.
                    for (uint64_t o = 0; o < nobs_f; ++o)
                        if (o < (uint64_t)ob.size() && ob[(size_t)o]) {
                            const uint64_t byte_idx = (o >> 6) * 8 + ((o & 63) >> 3);
                            preds[(size_t)k * ((nobs_f + 63) / 64 * 8) + byte_idx] |=
                                (uint8_t)(1u << (o & 7));
                        }
                    stats[(size_t)k] = st;
                }
            };
            const int tc_threads = std::max(1, args.threads);
            for (int t = 0; t < tc_threads; ++t) pool.emplace_back(worker, t);
            for (auto& th : pool) th.join();
        }
        const auto t1 = std::chrono::steady_clock::now();
        const double wall_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        // ---- output ----
        std::ofstream out(args.out, std::ios::binary);
        if (!out) fail("cannot write output: " + args.out);
        const uint32_t obs_words = (uint32_t)((nobs_f + 63) / 64);
        const uint32_t nobs32 = (uint32_t)nobs_f;
        const uint32_t nlim32 = (uint32_t)limit;
        std::vector<uint8_t> outbuf;
        outbuf.reserve(12 + (size_t)limit * 8 * obs_words);
        uint8_t hb[12];
        put_u32(hb + 0, kOutMagic);
        put_u32(hb + 4, nobs32);
        put_u32(hb + 8, nlim32);
        outbuf.insert(outbuf.end(), hb, hb + 12);
        // pack predictions: u64 little-endian words, word w holds bits 64w..64w+63,
        // bit b of the word = observable 64w+b (matches equiv.py reader)
        const size_t row_bytes = (size_t)((nobs_f + 63) / 64 * 8);  // u64-word padded row
        std::vector<uint8_t> prow((size_t)obs_words * 8, 0);
        for (int64_t k = 0; k < limit; ++k) {
            std::fill(prow.begin(), prow.end(), 0);
            const uint8_t* src = preds.data() + (size_t)k * row_bytes;
            for (uint64_t o = 0; o < nobs_f; ++o)
                if ((src[o >> 3] >> (o & 7)) & 1u) prow[(size_t)(o >> 3)] |= (uint8_t)(1u << (o & 7));
            outbuf.insert(outbuf.end(), prow.begin(), prow.end());
        }
        out.write((const char*)outbuf.data(), (std::streamsize)outbuf.size());
        out.close();

        // ---- stats to stderr ----
        ShotStat agg{};
        for (const auto& st : stats) {
            agg.seed_conv |= st.seed_conv;
            agg.returned_early |= st.returned_early;
            agg.converged_invalid |= st.converged_invalid;
            agg.fallback |= st.fallback;
            agg.bp_runs += st.bp_runs;
            agg.bp_iters += st.bp_iters;
            agg.rounds += st.rounds;
        }
        const double ms_per_shot = wall_ms / (double)args.repeats / (double)limit;
        std::fprintf(stderr,
                     "decoded %lld shots x %d reps in %.1f ms => %.3f ms/shot\n",
                     (long long)limit, args.repeats, wall_ms, ms_per_shot);
        std::fprintf(stderr,
                     "stats: seed_conv=%d early_return=%d conv_invalid=%d fallback=%d bp_runs=%lld bp_iters=%lld rounds=%d\n",
                     (int)agg.seed_conv, (int)agg.returned_early, (int)agg.converged_invalid,
                     (int)agg.fallback, (long long)agg.bp_runs, (long long)agg.bp_iters, agg.rounds);
        return 0;
    } catch (const std::exception& ex) {
        std::fprintf(stderr, "beam8_cpp: error: %s\n", ex.what());
        return 1;
    }
}
