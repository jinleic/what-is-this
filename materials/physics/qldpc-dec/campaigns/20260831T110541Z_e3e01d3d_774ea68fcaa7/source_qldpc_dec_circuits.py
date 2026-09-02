"""Circuit/DEM loading conventions for the harness.

Primary convention (all gates): the published stim circuits from the IonQ
beam-search repo (arXiv:2512.07057), committed under ../circuits/ with
sha256. See pre_statement.md. DEMs are derived consistently:
    dem = circuit.detector_error_model(decompose_errors=True,
                                       ignore_decomposition_failures=True)
Same call as the IBM Relay-BP repo's dem.rs (arXiv:2506.01779). Without
ignore_decomposition_failures, the Y-type errors emit 3-detector
components that fail graphlike split in stim for this circuit family.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import stim

CIRCUIT_DIR = Path(__file__).resolve().parent.parent.parent / "circuits"
DEM_DIR = Path(__file__).resolve().parent.parent.parent / "dem"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_circuit(p: float, basis: str = "Z") -> stim.Circuit:
    """Load the committed IonQ circuit for BB[[144,12,12]], 12 rounds.

    basis in {"X", "Z"}; primary convention: committed IonQ circuits
    (github.com/ionq-publications/beamsearchdecoder). p=3e-4 has no upstream
    circuit (verified 404); basis-specific p1e-3-rescale artifacts are used.
    """
    tag = {"X": "X", "Z": "Z"}[basis.upper()]
    candidates = [
        CIRCUIT_DIR / f"BB_144_144_12_memory_{tag}_p{p:g}_sr12_ionq.stim",
        CIRCUIT_DIR / f"BB_144_144_12_memory_{tag}_p{p}_sr12_ionq.stim",
    ]
    if p == 3e-4:
        # no upstream circuit at p=3e-4 (IonQ repo, verified 404 on 2026-08-29);
        # basis-specific artifact rescales every noise arg (0.001)->(0.0003)
        candidates.append(
            CIRCUIT_DIR
            / f"BB_144_144_12_memory_{tag}_p0.0003_sr12_derived_p1e-3_rescale.stim"
        )
    for path in candidates:
        if path.exists():
            return stim.Circuit.from_file(path)
    raise FileNotFoundError(f"no committed circuit for p={p}, basis={basis}: {candidates}")


def circuit_sha(p: float, basis: str = "Z") -> str:
    tag = basis.upper()
    candidates = [
        CIRCUIT_DIR / f"BB_144_144_12_memory_{tag}_p{p:g}_sr12_ionq.stim",
        CIRCUIT_DIR / f"BB_144_144_12_memory_{tag}_p{p}_sr12_ionq.stim",
    ]
    if p == 3e-4:
        candidates.append(
            CIRCUIT_DIR
            / f"BB_144_144_12_memory_{tag}_p0.0003_sr12_derived_p1e-3_rescale.stim"
        )
    for path in candidates:
        if path.exists():
            return sha256_file(path)
    raise FileNotFoundError(f"no committed circuit for p={p}, basis={basis}")


def load_dem(p: float, basis: str = "Z", decompose: bool = True) -> stim.DetectorErrorModel:
    """Deterministic DEM from the committed circuit.

    decompose=True -> decomposition attempted, per-component failures kept
    whole (ignore_decomposition_failures=True); this is the convention
    used by the IBM Relay-BP loader (relay crates/relay_bp/src/dem.rs,
    arXiv:2506.01779).
    """
    circuit = load_circuit(p, basis)
    if decompose:
        return circuit.detector_error_model(
            decompose_errors=True, ignore_decomposition_failures=True
        )
    return circuit.detector_error_model()


def dem_sha(p: float, basis: str = "Z") -> str:
    """Stable hash of the derived DEM (reproducibility tag)."""
    dem = load_dem(p, basis)
    return hashlib.sha256(str(dem).encode()).hexdigest()


if __name__ == "__main__":
    for p in (1e-3, 5e-4):
        for basis in ("X", "Z"):
            dem = load_dem(p, basis)
            print(
                f"p={p:g} basis={basis}: dets={dem.num_detectors} "
                f"errors={dem.num_errors} obs={dem.num_observables} sha={dem_sha(p, basis)[:12]}"
            )
