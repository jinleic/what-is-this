"""Gate C runner: exact ML vs beam8 / BP+OSD / NMS-ensemble24.

Scope: [[36,4,4]] BB code (l=6,m=3, Bravyi monomial family), X- or Z-memory
sector, code-capacity depolarizing noise with opposite-Pauli marginal
q = 2p/3. This is a small-scale exact-reference companion for qldpc-dec; it
is not the [[144,12,12]] headline configuration and cannot confirm or refute
the headline gates. See README.md in this directory.

For each sampled shot (syndrome s, true observable o):
  * exact degenerate-ML: full affine-coset enumeration (ExactMLReference),
    giving exact class masses Z_c; ML class c* = argmax_c Z_c (exact
    rational arithmetic, no truncation);
  * each decoder arm predicts the observable from s (harness decoders on the
    SAME stim DEM);
  * mismatch := predicted class != c*; regret := ln(Z_{c*} / Z_chosen) >= 0.
Decoder logical failure (vs the true o) is also recorded, so the report can
separate "decoder disagrees with ML" from "decoder (and ML) wrong".

Usage:
  python -m gatec_exact.run_gatec --p 1e-2 --shots 2000
  python -m gatec_exact.run_gatec --p 3e-2 --shots 2000
Outputs a JSON summary under qldpc-dec/scratch/gatec/ and prints a report.
"""
from __future__ import annotations

import argparse
import json
import hashlib
import platform
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import stim

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gatec_exact.exact_ml import (
    ExactMLReference,
    bb_distances,
    bb_logical_ops,
    bb_parity_checks,
    codecap_circuit,
    dem_matrices_raw,
    validate_exact_vs_bruteforce,
)
from qldpc_dec.seeds import derive_seed

SCRATCH = Path(__file__).resolve().parent.parent.parent / "scratch" / "gatec"


def build_instance():
    """Code and matrices for the BB [[36,4,4]] exact reference."""
    Hx, Hz = bb_parity_checks()
    Lz, Lx = bb_logical_ops(Hx, Hz)
    dx, dz = bb_distances(Hx, Hz)
    n = Hx.shape[1]
    k = Lz.shape[0]
    if not (n == 36 and k == 4 and (dx, dz) == (4, 4)):
        raise AssertionError(f"unexpected code parameters [[{n},{k},{dx},{dz}]]")
    return {"Hx": Hx, "Hz": Hz, "Lz": Lz, "Lx": Lx, "n": n, "k": k, "dx": dx, "dz": dz}


def codecap_circuit_zmem(q: float, inst) -> stim.Circuit:
    """Z-memory circuit: X errors, Z checks, and Z logicals."""
    return codecap_circuit(q, inst["Hz"], inst["Lz"], basis="Z")


def codecap_circuit_xmem(q: float, inst) -> stim.Circuit:
    """X-memory circuit: Z errors, X checks, and X logicals."""
    return codecap_circuit(q, inst["Hx"], inst["Lx"], basis="X")


def decoder_arms(dem: stim.DetectorErrorModel):
    """The three harness decoder arms bound to this DEM (compiled once)."""
    from qldpc_dec.beam_search import BeamSearchBatchDecoder
    from qldpc_dec.bp_osd import BpOsdBatchDecoder
    from qldpc_dec.nms import NmsEnsembleDecoder

    arms = {
        "beam8": BeamSearchBatchDecoder(
            dem, beam_width=8, initial_iters=30, iters_per_round=20,
            max_rounds=10, num_results=1,
        ),
        "bp30+osd": BpOsdBatchDecoder(dem, max_iter=30, bp_method="ms",
                                      schedule="parallel", osd_method="osd_cs",
                                      osd_order=10),
        "nms-ens24": NmsEnsembleDecoder(dem, ensemble_size=24, alpha=0.96875,
                                        max_iters=400, base_seed=0),
    }
    return arms


def obs_to_key(obs_row: np.ndarray, k: int) -> int:
    return int(np.dot(np.asarray(obs_row).ravel(), 1 << np.arange(k)))


def run_point(
    p: float,
    shots: int,
    base_seed: int,
    csv_path: Path | None = None,
    basis: str = "Z",
) -> dict:
    basis = basis.upper()
    if basis not in {"X", "Z"}:
        raise ValueError(f"basis must be X or Z, got {basis!r}")
    inst = build_instance()
    q = 2.0 * p / 3.0
    circ = (
        codecap_circuit_zmem(q, inst)
        if basis == "Z"
        else codecap_circuit_xmem(q, inst)
    )
    dem = circ.detector_error_model()
    H, A, pvec = dem_matrices_raw(dem)
    if not np.allclose(pvec, q):
        raise AssertionError("DEM priors do not match q")

    ref = ExactMLReference(H, A)
    validation_seed = (
        derive_seed(base_seed, "gatec", "validate")
        if basis == "Z"
        else derive_seed(base_seed, "gatec", "validate", basis)
    )
    t0 = time.perf_counter()
    val = validate_exact_vs_bruteforce(
        ref, seed=validation_seed, num_syndromes=2, wmax=5
    )
    val_s = time.perf_counter() - t0

    seed = (
        derive_seed(base_seed, "gatec", "sampler", p)
        if basis == "Z"
        else derive_seed(base_seed, "gatec", "sampler", basis, p)
    )
    det, obs = circ.compile_detector_sampler(seed=seed).sample(
        shots, separate_observables=True
    )

    arms = decoder_arms(dem)
    from qldpc_dec.dem_matrices import dem_to_matrices

    _bp_H, _bp_A, bp_priors = dem_to_matrices(dem, merge=True)
    configured_bp_priors = np.asarray(arms["bp30+osd"].decoder.channel_probs)
    np.testing.assert_array_equal(configured_bp_priors, bp_priors)
    preds = {}
    times = {}
    for name, dec in arms.items():
        t0 = time.perf_counter()
        preds[name] = dec.decode_batch(det)
        times[name] = time.perf_counter() - t0

    t0 = time.perf_counter()
    ml_decisions = np.zeros(shots, dtype=np.int64)
    ml_repents = {name: np.zeros(shots) for name in preds}
    ml_wrong = 0
    for i in range(shots):
        hist = ref.class_histogram(det[i])
        masses = ref.class_masses(hist, p)
        cstar, _ = ExactMLReference.ml_decision(masses)
        ml_decisions[i] = cstar
        true_c = obs_to_key(obs[i], inst["k"])
        if cstar != true_c:
            ml_wrong += 1
        for name, pr in preds.items():
            chosen = obs_to_key(pr[i], inst["k"])
            ml_repents[name][i] = 0.0 if chosen == cstar else ExactMLReference.regret_nats(masses, chosen)
    ml_ms = 1e3 * (time.perf_counter() - t0) / shots

    out = {
        "circuit_sha256": hashlib.sha256(str(circ).encode()).hexdigest(),
        "dem_sha256": hashlib.sha256(str(dem).encode()).hexdigest(),
        "dem_dimensions": {
            "detectors": dem.num_detectors,
            "errors": dem.num_errors,
            "observables": dem.num_observables,
        },
        "bposd_channel": {
            "model": "merged DEM heterogeneous error_channel",
            "columns": int(len(bp_priors)),
            "configured_exactly": True,
            "sha256_float64": hashlib.sha256(
                np.asarray(bp_priors, dtype="<f8").tobytes()
            ).hexdigest(),
            "reference": {
                "package": "stimbposd==0.1.0",
                "commit": "7921f5eb1b358ff616f9822280c9961e83df06cb",
            },
        },
        "code": {"family": "BB (l=6, m=3)", "n": inst["n"], "k": inst["k"],
                 "dx": inst["dx"], "dz": inst["dz"],
                 "monomials": {"A": "x^3 + y + y^2", "B": "y + x + x^3"}},
        "noise": {
            "model": "code-capacity depolarizing",
            "memory_basis": basis,
            "physical_pauli_error": "X" if basis == "Z" else "Z",
            "p_circuit_equiv": p,
            "q_error_marginal": q,
        },
        "shots": shots,
        "sampler_seed": seed,
        "exact_method": {"affine_space_log2": ref.dim,
                         "solutions_per_shot": 1 << ref.dim,
                         "rank": ref.rank,
                         "arithmetic": "exact Fractions for class argmax",
                         "bruteforce_crosscheck": val,
                         "bruteforce_crosscheck_s": val_s},
        "ml_logical_failures": ml_wrong,
        "ml_ler": ml_wrong / shots,
        "exact_ml_ms_per_shot": ml_ms,
        "decoders": {},
    }
    for name in sorted(preds):
        pr = preds[name]
        fails = int((pr != obs).any(axis=1).sum())
        pred_keys = np.array([obs_to_key(pr[i], inst["k"]) for i in range(shots)])
        mism = int((pred_keys != ml_decisions).sum())
        ml_right = (ml_decisions == np.array([obs_to_key(obs[i], inst["k"]) for i in range(shots)]))
        mism_mlright = int(((pred_keys != ml_decisions) & ml_right).sum())
        mism_mlwrong = int(((pred_keys != ml_decisions) & ~ml_right).sum())
        out["decoders"][name] = {
            "logical_failures": fails,
            "ler": fails / shots,
            "ml_mismatches": mism,
            "ml_mismatch_rate": mism / shots,
            "ml_mismatch_when_ml_right": mism_mlright,
            "ml_mismatch_when_ml_wrong": mism_mlwrong,
            "mean_regret_nats": float(np.mean(ml_repents[name])),
            "max_regret_nats": float(np.max(ml_repents[name])),
            "ms_per_shot": 1e3 * times[name] / shots,
        }
    if csv_path is not None:
        import csv as _csv
        with open(csv_path, "w", newline="") as f:
            w = _csv.writer(f)
            w.writerow(["shot", "true_obs", "ml_class"]
                       + [f"{name}_pred" for name in sorted(preds)]
                       + [f"{name}_regret" for name in sorted(preds)])
            for i in range(shots):
                row = [i, obs_to_key(obs[i], inst["k"]), int(ml_decisions[i])]
                for name in sorted(preds):
                    row.append(obs_to_key(preds[name][i], inst["k"]))
                for name in sorted(preds):
                    row.append(f"{ml_repents[name][i]:.6f}")
                w.writerow(row)
        out["per_shot_csv"] = str(csv_path)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=float, default=1e-2)
    ap.add_argument("--shots", type=int, default=2000)
    ap.add_argument("--base-seed", type=int, default=20260830)
    ap.add_argument("--basis", choices=("X", "Z"), default="Z")
    ap.add_argument("--csv", action="store_true", help="also write per-shot CSV to scratch")
    args = ap.parse_args()

    basis_tag = f"{args.basis.lower()}mem"
    csv_path = None
    if args.csv:
        csv_path = (
            SCRATCH
            / f"gatec_{basis_tag}_p{args.p:g}_{args.shots}shots_pershot.csv"
        )
    report = run_point(
        args.p,
        args.shots,
        args.base_seed,
        csv_path=csv_path,
        basis=args.basis,
    )
    report["host"] = {"machine": platform.machine(), "node": platform.node(),
                      "python": sys.version.split()[0]}
    SCRATCH.mkdir(parents=True, exist_ok=True)
    out_path = SCRATCH / f"gatec_{basis_tag}_p{args.p:g}_{args.shots}shots.json"
    out_path.write_text(json.dumps(report, indent=2))

    print(f"=== Gate C exact-ML reference — BB[[36,4,4]] "
          f"{args.basis}-memory, p={args.p:g} "
          f"(q={report['noise']['q_error_marginal']:.4g}), "
          f"{args.shots} shots ===")
    print(f"code: n={report['code']['n']} k={report['code']['k']} "
          f"d=({report['code']['dx']},{report['code']['dz']})  "
          f"exact space: 2^{report['exact_method']['affine_space_log2']} "
          f"solutions/shot, {report['exact_ml_ms_per_shot']:.1f} ms/shot")
    print(f"exact-ML logical failures: {report['ml_logical_failures']}/{args.shots} "
          f"(block failure rate {report['ml_ler']:.4g})")
    print(f"{'decoder':<12} {'failure':>9} {'ML-mismatch':>12} "
          f"{'mean regret [nats]':>20} {'max regret':>11} {'ms/shot':>9}")
    for name, decoder in report["decoders"].items():
        print(f"{name:<12} {decoder['ler']:>9.4f} "
              f"{decoder['ml_mismatch_rate']:>12.4f} "
              f"{decoder['mean_regret_nats']:>20.4f} "
              f"{decoder['max_regret_nats']:>11.4f} "
              f"{decoder['ms_per_shot']:>9.3f}")
    print(f"\nreport written: {out_path}")


if __name__ == "__main__":
    main()
