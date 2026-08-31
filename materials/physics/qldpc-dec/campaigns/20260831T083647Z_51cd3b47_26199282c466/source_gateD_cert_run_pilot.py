"""Gate-D certified-decoding pilot runner (arXiv:2608.25545, spacetime route).

Modes
  --mode mini   exact-ML sanity check: rotated surface code d=5, code
                capacity, 260 syndromes at p=0.10; AIS-CRN decisions +
                paired-bootstrap certificates vs exact coset enumeration
                (2^12 states/class; the paper's own Sec 8.1 ground truth).
  --mode gross  50-shot pilot on the frozen BB_144 Z-memory DEM
                (p=0.003 rescale circuit): beam_search seed -> 13 candidate
                classes (1 + 12 single-logical shifts, spacetime convention)
                -> lightened representatives -> CRN-AIS (T=64, K=64, q0=0.02)
                -> paired bootstrap (delta=0.05, Bonferroni /12, B=2500).

Outputs (machine-readable, scratch only, incremental flush):
  scratch/gateD_cert/mini_result.json
  scratch/gateD_cert/gross_result.json + gross_shots.jsonl
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import stim  # noqa: E402

from qldpc_dec.dem_matrices import dem_to_matrices  # noqa: E402
from qldpc_dec.seeds import derive_seed  # noqa: E402
from qldpc_dec.beam_search import BeamSearchDecoder  # noqa: E402

from .gf2 import gf2_rank, solve_gf2  # noqa: E402
from .model import GateDModel  # noqa: E402
from .ais import AISEngine  # noqa: E402
from .certificate import paired_bootstrap  # noqa: E402
from .mini_exact import rotated_surface_code, sample_shots, _independent_rows  # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRATCH = os.path.join(BASE, "scratch", "gateD_cert")
GROSS_CIRCUIT = os.path.join(
    BASE, "circuits", "BB_144_144_12_memory_Z_p0.003_sr12_derived_p1e-3_rescale.stim"
)
BASE_SEED = 20260830


def _jdump(path: str, obj: dict) -> None:
    with open(path, "w") as f:
        json.dump(obj, f, indent=1)


# ---------------------------------------------------------------- mini ----
def run_mini(n_shots: int = 260, T: int = 32, Kch: int = 128, p: float = 0.10,
             d: int = 5) -> dict:
    t_all = time.time()
    H, L = rotated_surface_code(d)
    n = d * d
    rank = gf2_rank(H)
    Hr = _independent_rows(H)
    assert gf2_rank(Hr) == rank
    lam_scalar = float(np.log((1 - p) / p))
    lam_vec = np.full(n, lam_scalar)

    shift = solve_gf2(np.vstack([H, L]).astype(np.uint8),
                      np.concatenate([np.zeros(H.shape[0], dtype=np.uint8),
                                      np.ones(1, dtype=np.uint8)]).astype(np.uint8))
    assert shift is not None, "logical shift solve failed"

    # exact machinery: 2^rank coset states
    Us = ((np.arange(1 << rank)[:, None] >> np.arange(rank)[None, :]) & 1).astype(np.uint8)
    coset_par = (Us @ Hr) % 2                        # (2^rank, n)

    nz_g, nz_m = np.nonzero(Hr)
    o = np.argsort(nz_g, kind="stable")
    row_ptr = np.concatenate([[0], np.cumsum(np.bincount(nz_g, minlength=rank))]).astype(np.int64)
    row_mechs = nz_m[o].astype(np.int64)

    syns, errs, obs = sample_shots(H, L, p, n_shots,
                                   seed=derive_seed(BASE_SEED, "gateD", "mini", "d5", "sampler"))
    agree = certs = cert_agree = n_used = 0
    exact_ml_success = 0
    ais_success = 0
    margin_exact = []
    t0 = time.time()
    for si in range(syns.shape[0]):
        s = syns[si]
        e0 = solve_gf2(H.astype(np.uint8), s.astype(np.uint8))
        assert e0 is not None
        reps = np.vstack([e0, e0.astype(np.uint8) ^ shift])       # (2, n)
        # exact ML (true ground truth)
        Ez = [((reps[c][None, :] ^ coset_par) @ lam_vec) for c in range(2)]
        logZ_exact = np.array([np.logaddexp.reduce(-Ez[c]) for c in range(2)])
        ml_exact = int(np.argmax(logZ_exact))
        exact_ml_success += int(ml_exact == int(obs[si, 0]))
        # AIS + certificate
        eng = AISEngine(row_mechs, row_ptr, lam_vec, T=T, K=Kch, q0=0.02,
                        seed=derive_seed(BASE_SEED, "gateD", "mini", "ais", si))
        out = eng.run(reps)
        cb = paired_bootstrap(out["logW"], B=2000, delta=0.05,
                              rng=np.random.default_rng(derive_seed(BASE_SEED, "gateD", "mini", "boot", si)))
        n_used += 1
        certs += int(cb["certified"])
        agree += int(cb["best"] == ml_exact)
        cert_agree += int(cb["certified"] and cb["best"] == ml_exact)
        ais_success += int(cb["best"] == int(obs[si, 0]))
        margin_exact.append(float(logZ_exact.max() - np.sort(logZ_exact)[-2]))
    dt = time.time() - t0
    res = {
        "mode": "mini_exact_vs_ais",
        "code": f"rotated-surface-d{d}-code-capacity",
        "p": p, "n_syndromes": n_used, "n_classes": 2, "coset_dim": int(rank),
        "ais": {"T": T, "K": Kch, "q0": 0.02, "B": 2000, "delta": 0.05},
        "agreement_ais_vs_exactML": agree,
        "agreement_rate": agree / max(n_used, 1),
        "certified": certs,
        "certified_rate": certs / max(n_used, 1),
        "certified_and_correct": cert_agree,
        "exact_ml_block_success": exact_ml_success / max(n_used, 1),
        "ais_block_success": ais_success / max(n_used, 1),
        "median_exact_margin_nats": float(np.median(margin_exact)),
        "wall_s_total": dt, "wall_s_per_syndrome": dt / max(n_used, 1),
        "wall_s_total_all": time.time() - t_all,
        "seed_base": BASE_SEED,
    }
    _jdump(os.path.join(SCRATCH, "mini_result.json"), res)
    print(json.dumps(res, indent=1))
    return res


# --------------------------------------------------------------- gross ----
def run_gross(
    n_shots: int = 50,
    T: int = 64,
    Kch: int = 64,
    q0: float = 0.02,
    B: int = 2500,
    delta: float = 0.05,
    out_dir: str | None = None,
) -> dict:
    t_all = time.time()
    output_dir = out_dir or SCRATCH
    os.makedirs(output_dir, exist_ok=True)
    circ = stim.Circuit.from_file(GROSS_CIRCUIT)
    dem = circ.detector_error_model(decompose_errors=True,
                                    ignore_decomposition_failures=True)
    Hs, As, lam = dem_to_matrices(dem, merge=False)
    H = Hs.toarray().astype(np.uint8)
    nz_rows = np.nonzero(np.diff(As.tocsr().indptr))[0]
    L = As.toarray()[nz_rows].astype(np.uint8)

    zb = np.load(os.path.join(SCRATCH, "gross_sparse_basis.npz"))
    K = zb["K"]
    shifts = np.load(os.path.join(SCRATCH, "gross_l_shifts.npy"))
    model = GateDModel(H, L, lam, sparse_basis=K, log_shifts=shifts)
    print(f"gross model: G={model.G} gens, n={model.n} mechs, "
          f"mean gen weight {model.K_weights.mean():.2f}", flush=True)

    sampler = dem.compile_sampler(seed=derive_seed(BASE_SEED, "gateD", "gross", "sampler"))
    dets, obs, _ = sampler.sample(n_shots)
    obs = np.asarray(obs, dtype=np.uint8)

    beam = BeamSearchDecoder(dem, beam_width=8, initial_iters=30,
                             iters_per_round=20, max_rounds=10)
    jsonl_path = os.path.join(output_dir, "gross_shots.jsonl")
    rows = []
    for si in range(n_shots):
        t0 = time.time()
        syn = dets[si].astype(np.uint8)
        rec = {"shot": si,
               "seed_ais": derive_seed(BASE_SEED, "gateD", "ais", si),
               "seed_boot": derive_seed(BASE_SEED, "gateD", "boot", si)}
        try:
            e0 = beam.decode(syn.astype(np.int64)).astype(np.uint8)
            rec["e0_weight"] = int(e0.sum())
            consistent = np.array_equal((H @ e0.astype(np.uint64)) % 2, syn)
            rec["beam_consistent"] = bool(consistent)
            if not consistent:
                rec["error"] = "beam seed not syndrome-consistent; certificate skipped"
                rows.append(rec)
                with open(jsonl_path, "a") as f:
                    f.write(json.dumps(rec) + "\n")
                print(f"shot {si}: INCONSISTENT beam seed, skipped", flush=True)
                continue
            reps = model.representatives(e0, lighten=True, max_sweeps=6)
            Hr_chk = (H @ reps.T.astype(np.uint64)) % 2          # (dets, 13)
            if not np.array_equal(Hr_chk.T.astype(np.uint8),
                                  np.broadcast_to(syn, Hr_chk.T.shape)):
                rec["error"] = "representative coset broke the syndrome"
                rows.append(rec)
                with open(jsonl_path, "a") as f:
                    f.write(json.dumps(rec) + "\n")
                print(f"shot {si}: reps inconsistent, skipped", flush=True)
                continue
            eng = AISEngine(model.row_mechs, model.row_ptr, lam, T=T, K=Kch,
                            q0=q0, seed=rec["seed_ais"])
            out = eng.run(reps)
            cb = paired_bootstrap(out["logW"], B=B, delta=delta,
                                  rng=np.random.default_rng(rec["seed_boot"]))
            # observable prediction of the decided class = L @ rep_best
            pred = (L @ reps[cb["best"]].astype(np.uint64)) % 2
            pred12 = pred[nz_rows] if len(pred) > L.shape[0] else pred
            true12 = obs[si][nz_rows] if obs.shape[1] > L.shape[0] else obs[si]
            rec.update({
                "decision_class": int(cb["best"]),          # 0 = beam class
                "certified": bool(cb["certified"]),
                "margin_nats": float(cb["margin_nats"]),
                "boot_p_fail": float(cb["boot_p_fail"]),
                "obs_correct": bool(np.array_equal(pred12.astype(np.uint8), true12.astype(np.uint8))),
                "wall_s": time.time() - t0,
            })
        except Exception as exc:  # keep the batch alive; record and continue
            rec["error"] = f"{type(exc).__name__}: {exc}"
            rec["wall_s"] = time.time() - t0
        rows.append(rec)
        with open(jsonl_path, "a") as f:
            f.write(json.dumps(rec) + "\n")
        print(f"shot {si}: class={rec.get('decision_class')} "
              f"cert={rec.get('certified')} margin={rec.get('margin_nats', float('nan')):.2f} "
              f"wall={rec.get('wall_s', 0):.1f}s", flush=True)

    ok = [r for r in rows if "decision_class" in r]
    certs = sum(int(r["certified"]) for r in ok)
    cert_ok = sum(int(r["certified"] and r.get("obs_correct")) for r in ok)
    res = {
        "mode": "gross_pilot",
        "circuit": os.path.basename(GROSS_CIRCUIT),
        "num_generators": int(model.G),
        "num_mechanisms": int(model.n),
        "num_detectors": int(model.num_dets),
        "mean_generator_weight": float(model.K_weights.mean()),
        "n_classes": int(reps.shape[0]),
        "ais": {"T": T, "K": Kch, "q0": q0, "B": B, "delta": delta,
                "bonf_competitors": int(reps.shape[0] - 1)},
        "n_shots": n_shots,
        "n_decoded": len(ok),
        "n_beam_inconsistent": sum(1 for r in rows if not r.get("beam_consistent", False)),
        "certified_fraction": certs / max(len(ok), 1),
        "certified": certs,
        "mean_margin_nats": float(np.mean([r["margin_nats"] for r in ok])) if ok else None,
        "median_wall_s": float(np.median([r["wall_s"] for r in ok])) if ok else None,
        "certified_and_obs_correct": cert_ok,
        "certified_obs_fail_rate": 1 - (cert_ok / certs) if certs else None,
        "decision_hist": {str(c): sum(1 for r in ok if r["decision_class"] == c)
                          for c in range(reps.shape[0])},
        "basis_seconds": model.basis_seconds,
        "wall_s_total": time.time() - t_all,
        "seed_base": BASE_SEED,
        "versions": {"numpy": np.__version__, "stim": stim.__version__},
    }
    _jdump(os.path.join(output_dir, "gross_result.json"), res)
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["mini", "gross"], required=True)
    ap.add_argument("--shots", type=int, default=None)
    ap.add_argument("--T", type=int, default=None)
    ap.add_argument("--K", type=int, default=None)
    ap.add_argument("--out-dir", default=None)
    a = ap.parse_args()
    if a.mode == "mini":
        run_mini(n_shots=a.shots or 260, T=a.T or 32, Kch=a.K or 128)
    else:
        run_gross(
            n_shots=a.shots or 50,
            T=a.T or 64,
            Kch=a.K or 64,
            out_dir=a.out_dir,
        )
