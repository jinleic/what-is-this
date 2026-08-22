#!/usr/bin/env python3
"""W_9: the ladder pairing-kernel (annihilator) dimension at L = 9, exactly.

Wave-16 extension of the sector-saturation certificate chain
(`proofs/sector_saturation_pairing.md`, `experiments/e142_sector_saturation_l8.py`)
from L = 8 to L = 9:

    dim K_9 = (2^17 + 3*2^9)/4 = 33152   (Burnside, three ways)
    dim U(A,B) psi = K_9 - W_9           (cyclic module of the vacuum)
    dim g_9 >= dim U(A,B) psi            (evaluation bound, Theorem 2 of
                                          proofs/ladder_alll_proof.md)

Two structural improvements over e142, both proved in
`proofs/ladder_l9.md` and verified here against every certified L = 3..8
value:

(1) *Sector-pure closure.*  A_L = 2L - 2N is a scalar on each
    particle sector, so on K_L every sector projector P_k is a polynomial
    in A_L and the cyclic module equals the closure of psi under
    {P_k . B}.  Closure rows stay sector-pure; A-images never need to be
    generated.

(2) *Particle-hole mirror.*  P = prod_v Z_v (bit complement) commutes with
    B_L, maps sector k to 2L - k, and P psi = |all-ones>.  Since
    D = (1/32)[A,[A,B]] - (1/8)[A,B] lies in the enveloping algebra and
    D^L psi reaches the fully occupied sector with a positive perfect-
    matching coefficient, U contains |all-ones>, hence P U = U and
    U_{2L-k} = P(U_k), W_{2L-k} = P(W_k).  The closure therefore runs on
    sectors k <= L only, mirror-injecting high-sector image components via
    P; kernels for k > L are the P-images of the low-sector lifts.

Certificate chain (identical sandwich to e142):
  (S1) modular rank r_q of the mirror closure at TWO primes
       P1 = 2147483647 (>= 2^30) and Q1 = 999983 (lift-safe);
  (S2) exact-Q pairing-kernel basis: per-sector triangular lift of the
       Q1 echelon (e142.kernel_lifts_dense, same code path as the L = 8
       certificate) + mirrored copies, validated over Q by
       e142.validate_W_over_Q (vacuum-orthogonal, sector-homogeneous,
       independent, A- and 4B-invariant);
  (S3) r_q + dim W = K_9 at both primes  =>  equality over Q.

Per-field claim tags: "EXACT-Q" for everything certified by the Q
validation (W_9, per-sector deficits, cyclic dim); "MODULAR-TWO-PRIME"
for the ranks; the growth model is [COMPUTATION]+[CONJECTURE].

Checkpoint discipline: the closure appends one JSON line per progress
event to results/ladder/w9_progress.jsonl and rewrites a full resume
checkpoint (echelon rows + todo queue) at results/ladder/_w9_work/
checkpoint_L<9>_<q>.npz after every block, so a kill loses at most one block.
Sectors finalise only at the closure fixpoint (images flow k -> k+-2),
so the per-sector {sector, dim, deficit, elapsed} lines are emitted at
completion from the final state.

Usage (main venv python):
  PYTHONPATH=src .venv/bin/python experiments/e151_ladder_w9.py regress
  PYTHONPATH=src .venv/bin/python experiments/e151_ladder_w9.py closure 999983 --save-rows
  PYTHONPATH=src .venv/bin/python experiments/e151_ladder_w9.py closure 2147483647
  PYTHONPATH=src .venv/bin/python experiments/e151_ladder_w9.py assemble
"""
from __future__ import annotations

import hashlib
import json
import os
import resource
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))
import e142_sector_saturation_l8 as e142  # noqa: E402  (frozen wave-15 producer)

SCRIPT = "experiments/e151_ladder_w9.py"
RESULT_PATH = ROOT / "results" / "ladder" / "w9_saturation.json"
PROGRESS_PATH = ROOT / "results" / "ladder" / "w9_progress.jsonl"
WORK = ROOT / "results" / "ladder" / "_w9_work"
L9 = 9
P1 = 2_147_483_647   # >= 2^30 (rank certificate)
Q1 = 999_983         # lift-safe prime (kernel_lifts_dense int64 contract)
CLOSURE_CPU_BUDGET = float(os.environ.get("E151_CLOSURE_BUDGET", "6500"))
PROGRESS_BLOCK = int(os.environ.get("E151_PROGRESS_BLOCK", "256"))
CHECKPOINT_BLOCK = int(os.environ.get("E151_CHECKPOINT_BLOCK", "512"))
CHECKPOINT_SECONDS = float(os.environ.get("E151_CHECKPOINT_SECONDS", "300"))
EXPECTED = {  # wave-14/e142 certified (cyclic, K, W, W sectors)
    3: (14, 14, 0, {}),
    4: (42, 44, 2, {"4": 2}),
    5: (142, 152, 10, {"4": 5, "6": 5}),
    6: (494, 560, 66, {"4": 15, "6": 36, "8": 15}),
    7: (1780, 2144, 364, {"4": 35, "6": 147, "8": 147, "10": 35}),
    8: (6562, 8384, 1822, {"10": 448, "12": 70, "4": 70, "6": 448, "8": 786}),
}


def peak_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


def log(message: str) -> None:
    print(f"[e151 {time.strftime('%H:%M:%S')}] {message}", flush=True)


def progress(event: dict) -> None:
    PROGRESS_PATH.parent.mkdir(parents=True, exist_ok=True)
    with PROGRESS_PATH.open("a") as handle:
        handle.write(json.dumps(event, sort_keys=True) + "\n")


# ---------------------------------------------------------------------------
# Mirror sector closure (pilots pilot_w9.py / pilot_mirror.py, verified L=3..8).
# ---------------------------------------------------------------------------
class MirrorClosure:
    """Sector-pure closure of psi under {P_k B} on sectors k <= L with
    particle-hole mirror-injection of high-sector image components.

    Semantics identical to e142.FastClosure / cyclic_module_mod_p (the
    pivot of every stored row is its maximum sector-local coordinate, rows
    are pivot-normalised, each stored row is processed exactly once) with
    two provable restrictions: rows are sector-pure (A acts as a sector
    scalar, so its images are always in span) and only low sectors are
    stored (high components are injected as their P-images).
    """

    def __init__(self, L: int, q: int) -> None:
        self.L, self.q = L, q
        self.orbit_of, self.members = e142.invariant_orbits(L)
        n = len(self.members)
        self.n = n
        self.sector_of = np.array([m[0].bit_count() for m in self.members], dtype=np.int64)
        self.sector_dims = Counter(self.sector_of.tolist())
        self.sectors = sorted(self.sector_dims)
        self.loc2glob: dict[int, np.ndarray] = {}
        self.glob2loc: dict[int, np.ndarray] = {}
        for s in self.sectors:
            coords = np.nonzero(self.sector_of == s)[0].astype(np.int64)
            self.loc2glob[s] = coords
            back = np.full(n, -1, dtype=np.int64)
            back[coords] = np.arange(coords.size, dtype=np.int64)
            self.glob2loc[s] = back
        B4 = e142.four_B_rows(L, self.orbit_of, self.members, q)
        indptr = np.zeros(n + 1, dtype=np.int64)
        idxs: list[int] = []
        data: list[int] = []
        for src, row in enumerate(B4):
            for tgt, c in sorted(row.items()):
                idxs.append(tgt)
                data.append(c)
            indptr[src + 1] = len(idxs)
        self.B_idx = np.array(idxs, dtype=np.int64)
        self.B_dat = np.array(data, dtype=np.int64)
        self.B_ptr = indptr
        allbits = (1 << (2 * L)) - 1
        self.P_idx = np.array(
            [self.orbit_of[self.members[i][0] ^ allbits] for i in range(n)], dtype=np.int64
        )
        if not np.all(self.sector_of[self.P_idx] == 2 * L - self.sector_of):
            raise AssertionError("particle-hole map does not pair sectors k <-> 2L-k")
        self.low = [s for s in self.sectors if s <= L]
        self.rows: dict[int, dict[int, tuple[np.ndarray, np.ndarray]]] = {s: {} for s in self.low}
        self.adds = 0
        self.reduction_steps = 0
        self.processed: set[tuple[int, int]] = set()
        self.todo: list[tuple[int, int]] = []

    # -- echelon ------------------------------------------------------------
    def add(self, s: int, idx: np.ndarray, vals: np.ndarray) -> int | None:
        q = self.q
        w = np.zeros(self.loc2glob[s].size, dtype=np.int64)
        w[idx] = vals % q
        rows = self.rows[s]
        self.adds += 1
        while True:
            nz = np.nonzero(w)[0]
            if nz.size == 0:
                return None
            pivot = int(nz[-1])
            stored = rows.get(pivot)
            if stored is None:
                inverse = pow(int(w[pivot]), -1, q)
                rows[pivot] = (nz.astype(np.int32).copy(), (w[nz] * inverse % q).copy())
                return pivot
            sidx, svals = stored
            self.reduction_steps += 1
            w[sidx] = (w[sidx] - int(w[pivot]) * svals) % q

    # -- checkpointing -----------------------------------------------------
    def checkpoint_path(self) -> Path:
        return WORK / f"checkpoint_L{self.L}_{self.q}.npz"

    def save_checkpoint(self, elapsed_cpu: float, wall: float) -> Path:
        path = self.checkpoint_path()
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {"q": self.q, "L": self.L, "adds": self.adds,
                   "reduction_steps": self.reduction_steps,
                   "elapsed_cpu": elapsed_cpu, "wall": wall,
                   "todo_sectors": np.array([t[0] for t in self.todo], dtype=np.int64),
                   "todo_pivots": np.array([t[1] for t in self.todo], dtype=np.int64),
                   "processed_sectors": np.array([p[0] for p in self.processed], dtype=np.int64),
                   "processed_pivots": np.array([p[1] for p in self.processed], dtype=np.int64)}
        for s in self.low:
            rows = self.rows[s]
            pivots = sorted(rows)
            payload[f"pivots_{s}"] = np.array(pivots, dtype=np.int64)
            payload[f"idx_{s}"] = np.concatenate(
                [rows[p][0].astype(np.int64) for p in pivots]) if pivots else np.zeros(0, np.int64)
            payload[f"ptr_{s}"] = np.concatenate(
                [[0], np.cumsum([rows[p][0].size for p in pivots])]).astype(np.int64) \
                if pivots else np.zeros(1, np.int64)
            payload[f"vals_{s}"] = np.concatenate(
                [rows[p][1] for p in pivots]) if pivots else np.zeros(0, np.int64)
        temporary = path.with_name(path.stem + ".tmp.npz")
        np.savez_compressed(temporary, **payload)
        temporary.replace(path)
        return path

    def load_checkpoint(self) -> bool:
        path = self.checkpoint_path()
        if not path.exists():
            return False
        with np.load(path) as data:
            self.adds = int(data["adds"])
            self.reduction_steps = int(data["reduction_steps"])
            self.todo = list(zip(data["todo_sectors"].tolist(), data["todo_pivots"].tolist()))
            self.processed = set(zip(data["processed_sectors"].tolist(),
                                     data["processed_pivots"].tolist()))
            for s in self.low:
                pivots = data[f"pivots_{s}"].tolist()
                ptr = data[f"ptr_{s}"]
                idxall = data[f"idx_{s}"]
                valsall = data[f"vals_{s}"]
                rows: dict[int, tuple[np.ndarray, np.ndarray]] = {}
                for column, p in enumerate(pivots):
                    lo, hi = int(ptr[column]), int(ptr[column + 1])
                    rows[int(p)] = (idxall[lo:hi].astype(np.int32), valsall[lo:hi].copy())
                self.rows[s] = rows
        log(f"resumed q={self.q} from {path.name}: "
            f"{len(self.processed)} processed, todo {len(self.todo)}")
        return True

    # -- closure -----------------------------------------------------------
    def run(self, cpu_budget: float) -> dict:
        q = self.q
        started = time.process_time()
        wall_started = time.monotonic()
        resumed = self.load_checkpoint()
        if not resumed:
            first = self.add(0, np.array([0], dtype=np.int64), np.array([1], dtype=np.int64))
            if first != 0:
                raise AssertionError("vacuum did not seed the mirror closure")
            self.todo = [(0, 0)]
        next_progress = max(len(self.processed) + PROGRESS_BLOCK,
                            PROGRESS_BLOCK * (1 if not resumed else 2))
        next_checkpoint = len(self.processed) + CHECKPOINT_BLOCK
        last_checkpoint_cpu = time.process_time() - started
        while self.todo:
            s, pivot = self.todo.pop()
            if (s, pivot) in self.processed:
                continue
            self.processed.add((s, pivot))
            idx_l, vals = self.rows[s][pivot]
            glob = self.loc2glob[s][idx_l.astype(np.int64)]
            image = np.zeros(self.n, dtype=np.int64)
            for j in range(glob.size):
                src = int(glob[j])
                lo, hi = self.B_ptr[src], self.B_ptr[src + 1]
                image[self.B_idx[lo:hi]] += int(vals[j]) * self.B_dat[lo:hi]
                image[self.B_idx[lo:hi]] %= q
            targets = np.nonzero(image)[0]
            if targets.size:
                tgt_sectors = self.sector_of[targets]
                for t in np.unique(tgt_sectors):
                    t = int(t)
                    mask = tgt_sectors == t
                    tg = targets[mask]
                    vals_part = image[tg].copy()
                    if t > self.L:  # mirror-inject into sector 2L - t
                        t = 2 * self.L - t
                        tg = self.P_idx[tg]
                    tl = self.glob2loc[t][tg]
                    new_pivot = self.add(t, tl, vals_part)
                    if new_pivot is not None:
                        self.todo.append((t, new_pivot))
            cpu_now = time.process_time() - started
            if len(self.processed) >= next_progress:
                progress({"q": q, "event": "progress", "processed": len(self.processed),
                          "low_rank": sum(len(self.rows[s]) for s in self.low),
                          "sector_ranks": {str(s): len(self.rows[s]) for s in self.low},
                          "adds": self.adds, "steps": self.reduction_steps,
                          "cpu_s": round(cpu_now, 1),
                          "wall_s": round(time.monotonic() - wall_started, 1)})
                log(f"closure q={q}: processed {len(self.processed)}, "
                    f"low-rank {sum(len(self.rows[s]) for s in self.low)}, cpu {cpu_now:.0f}s")
                next_progress = len(self.processed) * 2
            if (len(self.processed) >= next_checkpoint
                    or cpu_now - last_checkpoint_cpu > CHECKPOINT_SECONDS):
                path = self.save_checkpoint(cpu_now, time.monotonic() - wall_started)
                progress({"q": q, "event": "checkpoint", "processed": len(self.processed),
                          "path": str(path.relative_to(ROOT)), "cpu_s": round(cpu_now, 1)})
                next_checkpoint = len(self.processed) + CHECKPOINT_BLOCK
                last_checkpoint_cpu = cpu_now
            if cpu_now > cpu_budget:
                path = self.save_checkpoint(cpu_now, time.monotonic() - wall_started)
                state = {"q": q, "processed": len(self.processed),
                         "low_rank": sum(len(self.rows[s]) for s in self.low),
                         "sector_ranks": {str(s): len(self.rows[s]) for s in self.low},
                         "adds": self.adds, "steps": self.reduction_steps,
                         "cpu_s": round(cpu_now, 1), "checkpoint": str(path.relative_to(ROOT))}
                progress({"q": q, "event": "budget_stop", **state})
                raise TimeoutError(f"closure q={q} exceeded {cpu_budget}s cpu: {state}")
        low_ranks = {s: len(self.rows[s]) for s in self.low}
        rank_total = sum(low_ranks[s] * (1 if s == 2 * self.L - s else 2) for s in self.low)
        stats = {
            "prime": q,
            "low_ranks": {str(s): low_ranks[s] for s in self.low},
            "rank_total": rank_total,
            "processed_rows": len(self.processed),
            "adds": self.adds,
            "reduction_steps": self.reduction_steps,
            "max_row_nnz": max((idx.size for s in self.low for idx, _ in self.rows[s].values()),
                               default=0),
            "cpu_seconds": round(time.process_time() - started, 1),
            "wall_seconds": round(time.monotonic() - wall_started, 1),
            "claim_tag": "MODULAR-TWO-PRIME-LOWER-BOUND",
        }
        progress({"q": q, "event": "closure_done", **stats})
        for s in self.low:
            progress({"q": q, "event": "sector_final", "sector": s,
                      "dim": self.sector_dims[s], "rank": low_ranks[s],
                      "deficit": self.sector_dims[s] - low_ranks[s],
                      "mirrored_sector": 2 * self.L - s,
                      "mirrored_deficit": self.sector_dims[2 * self.L - s] - low_ranks[s],
                      "cpu_s": stats["cpu_seconds"]})
        log(f"closure q={q} DONE: rank {rank_total} of K_{self.L}="
            f"{sum(self.sector_dims.values())} in {stats['cpu_seconds']}s cpu")
        return stats

    def save_rows(self) -> Path:
        path = WORK / f"rows_{self.q}.npz"
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {"q": self.q, "L": self.L, "n": self.n}
        for s in self.low:
            rows = self.rows[s]
            pivots = sorted(rows)
            payload[f"pivots_{s}"] = np.array(pivots, dtype=np.int64)
            payload[f"idx_{s}"] = np.concatenate(
                [rows[p][0].astype(np.int64) for p in pivots]) if pivots else np.zeros(0, np.int64)
            payload[f"ptr_{s}"] = np.concatenate(
                [[0], np.cumsum([rows[p][0].size for p in pivots])]).astype(np.int64) \
                if pivots else np.zeros(1, np.int64)
            payload[f"vals_{s}"] = np.concatenate(
                [rows[p][1] for p in pivots]) if pivots else np.zeros(0, np.int64)
        temporary = path.with_name(path.stem + ".tmp.npz")
        np.savez_compressed(temporary, **payload)
        temporary.replace(path)
        return path

    def load_rows(self) -> None:
        path = WORK / f"rows_{self.q}.npz"
        with np.load(path) as data:
            for s in self.low:
                pivots = data[f"pivots_{s}"].tolist()
                ptr = data[f"ptr_{s}"]
                idxall = data[f"idx_{s}"]
                valsall = data[f"vals_{s}"]
                rows: dict[int, tuple[np.ndarray, np.ndarray]] = {}
                for column, p in enumerate(pivots):
                    lo, hi = int(ptr[column]), int(ptr[column + 1])
                    rows[int(p)] = (idxall[lo:hi].astype(np.int32), valsall[lo:hi].copy())
                self.rows[s] = rows
        log(f"loaded stored echelon rows for q={self.q} from {path.name}")

    # -- kernel lift -------------------------------------------------------
    def lift_low_sectors(self) -> tuple[list[dict[int, int]], dict]:
        sizes = [len(o) for o in self.members]
        candidates: list[dict[int, int]] = []
        detail: dict[str, dict] = {}
        for s in self.low:
            rows_s = self.rows[s]
            loc2glob = self.loc2glob[s]
            dim_s = int(loc2glob.size)
            rows_local = {
                int(p): {int(i): int(v) for i, v in zip(idx, vals)}
                for p, (idx, vals) in rows_s.items()
            }
            sizes_local = [int(sizes[g]) for g in loc2glob]
            cands = e142.kernel_lifts_dense(rows_local, dim_s, sizes_local, self.q, bound=1 << 20)
            n_free = dim_s - len(rows_s)
            if len(cands) != n_free:
                raise AssertionError(f"sector {s}: {len(cands)} lifts != {n_free} free coords")
            for cand in cands:
                candidates.append({int(loc2glob[i]): int(c) for i, c in cand.items()})
            detail[str(s)] = {"dim": dim_s, "rank": len(rows_s), "free": n_free}
        return candidates, detail

    def mirror_candidates(self, low_candidates: list[dict[int, int]]) -> list[dict[int, int]]:
        out: list[dict[int, int]] = []
        for cand in low_candidates:
            g0 = next(iter(cand))
            s = int(self.sector_of[g0])
            if s == 2 * self.L - s:
                continue
            out.append({int(self.P_idx[i]): c for i, c in cand.items()})
        return out


# ---------------------------------------------------------------------------
# Stages.
# ---------------------------------------------------------------------------
def stage_regress() -> dict:
    """Mirror closure regression against the certified L = 3..8 values."""
    started = time.process_time()
    records = []
    for L in (3, 4, 5, 6, 7, 8):
        engine = MirrorClosure(L, Q1)
        stats = engine.run(1800.0)
        low_cands, detail = engine.lift_low_sectors()
        full = low_cands + engine.mirror_candidates(low_cands)
        validation = e142.validate_W_over_Q(L, full, engine.orbit_of, engine.members)
        sectors = Counter(v[0] for v in validation["sector_particle_counts"])
        W_sectors = {str(k): sectors[k] for k in sorted(sectors)}
        cyclic, K, W, Wsec = EXPECTED[L]
        record = {
            "L": L,
            "K_dim": K,
            "rank_total": stats["rank_total"],
            "cyclic_Q_dim": K - validation["unreachable_block_dim"],
            "W_dim": validation["unreachable_block_dim"],
            "W_sectors": W_sectors,
            "validation": {k: v for k, v in validation.items() if k != "sector_particle_counts"},
            "max_abs_coefficient": max((abs(c) for v in full for c in v.values()), default=0),
            "closure_cpu_seconds": stats["cpu_seconds"],
            "claim_tag": "EXACT-Q",
        }
        if (record["rank_total"], K, record["W_dim"], W_sectors) != (cyclic, K, W, Wsec):
            raise AssertionError(f"regression mismatch at L={L}: {record}")
        records.append(record)
        log(f"regress L={L}: cyclic {cyclic}/{K}, W={W}, sectors {W_sectors} OK "
            f"({stats['cpu_seconds']}s)")
    out = {"records": records, "cpu_seconds": round(time.process_time() - started, 1)}
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "regression.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    return out


def stage_closure(prime: int, save_rows: bool) -> dict:
    started = time.process_time()
    engine = MirrorClosure(L9, prime)
    stats = engine.run(CLOSURE_CPU_BUDGET)
    if save_rows:
        path = engine.save_rows()
        stats["rows_file"] = str(path.relative_to(ROOT))
    stats["build_cpu_seconds"] = round(time.process_time() - started - stats["cpu_seconds"], 1)
    stats["peak_rss_bytes"] = peak_rss_bytes()
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / f"closure_{prime}.json").write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    return stats


def growth_model(W9: int, K9: int) -> dict:
    """Fit W_L ~ C r^L on the seven exact points, predict the crossing."""
    import mpmath as mp
    mp.mp.dps = 50
    Ls = list(range(3, 10))
    Ws = [0, 2, 10, 66, 364, 1822, W9]
    Ks = [e142.burnside_K_dim(L) for L in Ls]
    fit_L = [4, 5, 6, 7, 8, 9]
    fit_W = Ws[1:]
    x = [mp.mpf(L) for L in fit_L]
    y = [mp.log(mp.mpf(W)) for W in fit_W]
    n = mp.mpf(len(x))
    sx, sy = sum(x), sum(y)
    sxx = sum(t * t for t in x)
    sxy = sum(a * b for a, b in zip(x, y))
    slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    intercept = (sy - slope * sx) / n
    r_hat = mp.exp(slope)
    C_hat = mp.exp(intercept)
    residuals = [float(mp.log(mp.mpf(W)) - (intercept + slope * mp.mpf(L)))
                 for L, W in zip(fit_L, fit_W)]
    max_rel_residual = max(abs(float(1 - mp.exp(res))) for res in residuals)
    all_ratios = [mp.mpf(Ws[i + 1]) / Ws[i] for i in range(1, len(Ws) - 1)]
    prior_ratios = all_ratios[:4]             # 5.0, 6.6, 5.5151..., 5.0055... (through L=8)
    new_ratio = all_ratios[-1]                # W9 / 1822 (the object under test)
    ratio_min, ratio_max = min(prior_ratios), max(prior_ratios)
    prior_last = prior_ratios[-1]
    prediction = 1822 * prior_last
    deviation_from_prior_trend = abs(mp.mpf(W9) - prediction) / prediction

    def model_W(L, r, C):
        return C * mp.power(r, L)

    def solve_crossing(r):
        """First L where the model W_L >= K_L, by sign-change + bisection."""
        def f(L):
            K = (mp.power(2, 2 * L - 1) + 3 * mp.power(2, L)) / 4
            return model_W(L, r, C_hat) / K - 1
        lo, hi = mp.mpf(9), mp.mpf(9.5)
        if f(lo) >= 0:
            return lo
        while f(hi) < 0:
            lo, hi = hi, hi + (hi - mp.mpf(8))
            if hi > 100:
                return None
        return mp.findroot(f, (lo, hi))

    L_star = solve_crossing(r_hat)

    def solve_crossing_anchored(r):
        """Crossing L for the anchored model W_L = W9 * r^(L-9)."""
        def f(L):
            K = (mp.power(2, 2 * L - 1) + 3 * mp.power(2, L)) / 4
            return mp.mpf(W9) * mp.power(r, L - 9) / K - 1
        lo, hi = mp.mpf(9), mp.mpf(9.5)
        if f(lo) >= 0:
            return lo
        while f(hi) < 0:
            lo, hi = hi, hi + (hi - mp.mpf(8))
            if hi > 100:
                return None
        return mp.findroot(f, (lo, hi))

    L_star_min = solve_crossing_anchored(ratio_min)
    L_star_max = solve_crossing_anchored(ratio_max)
    L_star_new = solve_crossing_anchored(new_ratio)

    def vacuity_2powL(model):
        for L in range(10, 40):
            if model(L) > mp.mpf(e142.burnside_K_dim(L)) - mp.power(2, L):
                return L
        return None

    lsq_model = lambda L: model_W(L, r_hat, C_hat)
    anchored_model = lambda L: mp.mpf(W9) * mp.power(new_ratio, L - 9)

    table = []
    for L in range(10, 21):
        mw = model_W(mp.mpf(L), r_hat, C_hat)
        K = e142.burnside_K_dim(L)
        table.append({
            "L": L, "K_dim": K, "model_W": float(mw),
            "model_W_over_K": float(mw / K),
            "model_cyclic_K_minus_W": float(mp.mpf(K) - mw),
            "2powL": 2 ** L,
            "certificate_2powL_vacuous": bool(mw > mp.mpf(K) - mp.power(2, L)),
        })
    return {
        "points_L": Ls,
        "W_values": Ws,
        "K_values": Ks,
        "W_over_K": [round(W / K, 6) for W, K in zip(Ws, Ks)],
        "consecutive_ratios": [float(r) for r in all_ratios],
        "prior_ratios_through_L8": [float(r) for r in prior_ratios],
        "new_ratio_W9_over_W8": float(new_ratio),
        "fit": {
            "model": "W_L = C * r^L, least squares on log W over L=4..9",
            "mpmath_dps": 50,
            "r_hat": float(r_hat),
            "C_hat": float(C_hat),
            "log_residuals": residuals,
            "max_relative_residual": max_rel_residual,
        },
        "ratio_band_through_L8": {"min": float(ratio_min), "max": float(ratio_max),
                                  "prior_last": float(prior_last)},
        "prediction_from_prior_last_ratio": float(prediction),
        "W9_deviation_from_prior_trend": float(deviation_from_prior_trend),
        "W9_deviation_from_prior_trend_exceeds_20pct": bool(deviation_from_prior_trend > mp.mpf("0.2")),
        "crossing_W_ge_K": {
            "continuous_L_star_lsq_r_hat": float(L_star),
            "first_integer_L_lsq_r_hat": int(mp.ceil(L_star)),
            "continuous_L_star_anchored_ratio_min": float(L_star_min) if L_star_min is not None else None,
            "first_integer_L_anchored_ratio_min": int(mp.ceil(L_star_min)) if L_star_min is not None else None,
            "continuous_L_star_anchored_ratio_max": float(L_star_max),
            "first_integer_L_anchored_ratio_max": int(mp.ceil(L_star_max)),
            "continuous_L_star_anchored_new_ratio": float(L_star_new),
            "first_integer_L_anchored_new_ratio": int(mp.ceil(L_star_new)),
        },
        "certificate_2powL_vacuity_L_lsq_r_hat": vacuity_2powL(lsq_model),
        "certificate_2powL_vacuity_L_anchored_new_ratio": vacuity_2powL(anchored_model),
        "extrapolation_table_L10_L20": table,
        "claim_tag": "[COMPUTATION]+[CONJECTURE] (seven exact points do not constitute a law; "
                     "the fitted geometric model has no proof)",
    }


def stage_assemble() -> dict:
    started = time.process_time()
    # K_9 three ways.
    K9_record = e142.K_dimension_record(L9)
    K9 = K9_record["K_dim_formula"]
    if K9 != 33152:
        raise AssertionError(f"K_9 != 33152: {K9_record}")
    # Closures.
    closure_q1 = json.loads((WORK / f"closure_{Q1}.json").read_text())
    closure_p1 = json.loads((WORK / f"closure_{P1}.json").read_text())
    if closure_q1["rank_total"] != closure_p1["rank_total"]:
        raise AssertionError("prime rank disagreement")
    # Exact-Q kernel from the stored Q1 echelon.
    engine = MirrorClosure(L9, Q1)
    engine.load_rows()
    t_lift = time.process_time()
    low_cands, detail = engine.lift_low_sectors()
    mirror_cands = engine.mirror_candidates(low_cands)
    full = low_cands + mirror_cands
    lift_cpu = time.process_time() - t_lift
    validation = e142.validate_W_over_Q(L9, full, engine.orbit_of, engine.members)
    W9 = validation["unreachable_block_dim"]
    cyclic = K9 - W9
    if cyclic != closure_q1["rank_total"]:
        raise AssertionError(f"sandwich fails: cyclic {cyclic} != modular rank "
                             f"{closure_q1['rank_total']}")
    sectors = Counter(v[0] for v in validation["sector_particle_counts"])
    W_sectors = {str(k): sectors[k] for k in sorted(sectors)}
    if sum(W_sectors.values()) != W9:
        raise AssertionError("sector split does not sum to W_9")
    if not all(W_sectors.get(str(2 * L9 - int(k))) == v for k, v in W_sectors.items()):
        raise AssertionError("W_9 sector split not palindromic under k -> 18-k")
    basis = e142.serialise_basis(full, engine.members)
    canonical = json.dumps(basis, sort_keys=True, separators=(",", ":"))
    basis_sha256 = hashlib.sha256(canonical.encode()).hexdigest()
    regression = json.loads((WORK / "regression.json").read_text())
    model = growth_model(W9, K9)
    low_deficits = {}
    for s in engine.low:
        rank_s = len(engine.rows[s])
        low_deficits[str(s)] = {
            "dim_K_sector": int(engine.sector_dims[s]),
            "rank": rank_s,
            "deficit": int(engine.sector_dims[s]) - rank_s,
            "mirrored_sector": 2 * L9 - s,
            "mirrored_dim": int(engine.sector_dims[2 * L9 - s]),
            "mirrored_deficit": int(engine.sector_dims[2 * L9 - s]) - rank_s,
        }
    checks = [
        {"name": "C_K9_three_ways",
         "passed": K9_record["K_dim_formula"] == K9_record["K_dim_direct_orbits"] == 33152
                   and sum(K9_record["K_dim_sector_burnside"].values()) == 33152
                   and sum(K9_record["K_dim_sector_direct"].values()) == 33152,
         "detail": "dim K_9 = (2^17+3*2^9)/4 = 33152 by Burnside formula, direct orbits, per-sector Burnside"},
        {"name": "C_L9_two_prime_rank_agreement",
         "passed": closure_q1["rank_total"] == closure_p1["rank_total"] == cyclic,
         "detail": f"mirror closure ranks {closure_q1['rank_total']} (q=999983) and "
                   f"{closure_p1['rank_total']} (p=2147483647) equal the exact cyclic dim"},
        {"name": "C_L9_sandwich",
         "passed": all(validation[k] for k in ("vacuum_orthogonal", "sector_homogeneous",
                                               "basis_independent_over_Q", "four_B_invariant_over_Q",
                                               "A_invariant_over_Q"))
                   and cyclic + W9 == K9,
         "detail": f"modular lower rank + exact-Q validated kernel => dim U(A,B)psi = "
                   f"{cyclic} = K_9 - W_9 exactly over Q"},
        {"name": "C_L9_sector_palindrome",
         "passed": all(W_sectors.get(str(2 * L9 - int(k))) == v for k, v in W_sectors.items()),
         "detail": f"W_9 sector split {W_sectors} palindromic under k -> 18-k (particle-hole mirror theorem)"},
        {"name": "C_L9_evaluation_certificate_form",
         "passed": cyclic >= 2 ** L9,
         "detail": f"dim g_9 >= dim U(A,B)psi = {cyclic} >= 2^9; certificate form "
                   f"dim g_9 >= K_9 - W_9 = {K9} - {W9} preserved"},
        {"name": "C_regression_L3_L8",
         "passed": all((r["rank_total"], r["W_dim"], r["cyclic_Q_dim"]) ==
                       (EXPECTED[r["L"]][0], EXPECTED[r["L"]][2], EXPECTED[r["L"]][0])
                       for r in regression["records"]),
         "detail": "mirror engine reproduces the certified cyclic dims and W dims 14/14, 42/44, "
                   "142/152, 494/560, 1780/2144, 6562/8384 with W = 0,2,10,66,364,1822"},
        {"name": "C_W_pattern_sequence",
         "passed": [0, 2, 10, 66, 364, 1822, W9] == [r["W_dim"] for r in regression["records"]] + [W9],
         "detail": f"W sequence L=3..9 = 0, 2, 10, 66, 364, 1822, {W9}"},
        {"name": "C_growth_horizon_present",
         "passed": model["crossing_W_ge_K"]["first_integer_L_lsq_r_hat"] is not None
                   and model["W9_deviation_from_prior_trend"] is not None,
         "detail": "growth model fitted on seven exact points with honest [CONJECTURE] tag"},
    ]
    data = {
        "scope": {
            "space": "K_9 = even-parity <tau,rho>-invariant ladder configuration space, "
                     "orbit-sum basis, 33152 orbits of 2^17 even configurations",
            "certificate": "sector-pure mirror closure rank at two primes + exact-Q validated "
                           "pairing kernel (e142.validate_W_over_Q) => equality",
            "engines": "experiments/e151_ladder_w9.py MirrorClosure (sector-pure rows, "
                       "particle-hole mirror on sectors k <= 9); kernel lift via "
                       "e142.kernel_lifts_dense on the q=999983 echelon",
            "claim_tags_legend": {
                "EXACT-Q": "certified over Q by the validated kernel basis",
                "MODULAR-TWO-PRIME-LOWER-BOUND": "rank over F_q at the listed primes",
                "[COMPUTATION]+[CONJECTURE]": "exact data with an unproved fitted model",
            },
        },
        "K9_dimension_record": K9_record,
        "regression_records": regression["records"],
        "closures": {"q1": closure_q1, "p1": closure_p1},
        "W9": {
            "value": W9,
            "claim_tag": "EXACT-Q",
            "cyclic_Q_dim": cyclic,
            "cyclic_claim_tag": "EXACT-Q",
            "certificate_form": f"dim g_9 >= K_9 - W_9 = {K9} - {W9} = {cyclic}",
            "W_sectors": W_sectors,
            "low_sector_deficits": low_deficits,
            "max_abs_coefficient": max((abs(c) for v in full for c in v.values()), default=0),
            "lift_cpu_seconds": round(lift_cpu, 1),
            "validation": {k: v for k, v in validation.items() if k != "sector_particle_counts"},
            "basis_sha256": basis_sha256,
            "basis_size": len(basis),
        },
        "growth_model": model,
        "resource_measurements": {
            "assemble_cpu_seconds": round(time.process_time() - started, 1),
            "closure_cpu_seconds_q1": closure_q1["cpu_seconds"],
            "closure_cpu_seconds_p1": closure_p1["cpu_seconds"],
            "peak_rss_bytes": peak_rss_bytes(),
            "budget_clock": "time.process_time",
        },
    }
    provenance = {
        "artifact_schema": "provenance/data/checks-v1",
        "script": SCRIPT,
        "field": "Q and F_2147483647 x F_999983",
        "exact_arithmetic": "Python int, numpy int64 (mod q, per-source reduced accumulation), "
                            "fractions.Fraction",
        "artifact_sha256": "",
    }
    envelope = {"provenance": provenance, "data": data, "checks": checks}
    canonical = json.dumps(
        {"provenance": {k: v for k, v in provenance.items() if k != "artifact_sha256"},
         "data": data, "checks": checks},
        sort_keys=True, separators=(",", ":"),
    ).encode()
    envelope["provenance"]["artifact_sha256"] = hashlib.sha256(canonical).hexdigest()
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    temporary = RESULT_PATH.with_name(f".{RESULT_PATH.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(envelope, indent=2, sort_keys=True) + "\n")
    temporary.replace(RESULT_PATH)
    # The basis itself is large; store it next to the artifact.
    (RESULT_PATH.parent / "w9_basis.json").write_text(
        json.dumps({"basis": basis, "sha256": basis_sha256,
                    "claim_tag": "EXACT-Q"}, indent=0, sort_keys=True) + "\n")
    log(f"assembled: W_9 = {W9}, cyclic {cyclic}/{K9}, sectors {W_sectors}")
    log(f"crossing r_hat: L* = "
        f"{model['crossing_W_ge_K']['continuous_L_star_lsq_r_hat']:.2f} (continuous), "
        f"first integer L = {model['crossing_W_ge_K']['first_integer_L_lsq_r_hat']}")
    return envelope


def main() -> int:
    args = sys.argv[1:]
    if not args or args[0] == "all":
        stage_regress()
        stage_closure(Q1, save_rows=True)
        stage_closure(P1, save_rows=False)
        stage_assemble()
        return 0
    stage = args[0]
    if stage == "regress":
        stage_regress()
    elif stage == "closure":
        prime = int(args[1])
        stage_closure(prime, save_rows="--save-rows" in args)
    elif stage == "assemble":
        stage_assemble()
    else:
        raise SystemExit(f"unknown stage {stage}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
