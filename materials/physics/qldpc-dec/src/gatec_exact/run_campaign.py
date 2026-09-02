"""Freeze corrected Gate-C exact-reference campaigns for one memory basis."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import shutil
from pathlib import Path

import numpy as np
import stim

from gatec_exact.run_gatec import run_point
from qldpc_dec.runner import CAMPAIGNS_ROOT, Campaign

TARGET = Path(__file__).resolve().parents[2]
SRC = TARGET / "src"
POINTS = (1e-2, 3e-2)
SOURCE_PATHS = {
    "gatec_exact_run_campaign.py": Path(__file__).resolve(),
    "gatec_exact_run_gatec.py": Path(__file__).resolve().with_name("run_gatec.py"),
    "gatec_exact_exact_ml.py": Path(__file__).resolve().with_name("exact_ml.py"),
    "qldpc_dec_bp_osd.py": SRC / "qldpc_dec" / "bp_osd.py",
    "qldpc_dec_dem_matrices.py": SRC / "qldpc_dec" / "dem_matrices.py",
    "qldpc_dec_beam_search.py": SRC / "qldpc_dec" / "beam_search.py",
    "qldpc_dec_nms.py": SRC / "qldpc_dec" / "nms.py",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def run_campaign(
    *,
    basis: str,
    shots: int,
    base_seed: int,
    points: tuple[float, ...] = POINTS,
    smoke: bool = False,
) -> Path:
    basis = basis.upper()
    if basis not in {"X", "Z"}:
        raise ValueError(f"basis must be X or Z, got {basis!r}")
    if shots <= 0 or not points:
        raise ValueError("shots and points must be nonzero")

    prior_campaign = {
        "Z": "20260830T184327Z_53a194ad_e663d5477de1",
        "X": "20260831T085112Z_df905ed6_9b30cfe313ed",
    }[basis]
    config = {
        "run": "gateC-exact-GB3-corrected-smoke" if smoke else "gateC-exact-GB3-corrected",
        "basis": basis,
        "points": list(points),
        "shots_per_point": shots,
        "base_seed": base_seed,
        "scope": "BB [[36,4,4]] code-capacity exact reference",
        "bposd_prior_model": "merged DEM heterogeneous error_channel",
        "correction": "pre_statement.md Revision GB3",
        "supersedes_bposd_column_in": prior_campaign,
    }
    root = TARGET / "campaigns-smoke" if smoke else CAMPAIGNS_ROOT()
    campaign = Campaign(config, root=root)

    source_hashes: dict[str, str] = {}
    for label, source in SOURCE_PATHS.items():
        snapshot = campaign.dir / f"source_{label}"
        shutil.copy2(source, snapshot)
        source_hashes[snapshot.name] = sha256_file(snapshot)
    campaign.write_manifest(
        {
            "source_sha256": source_hashes,
            "versions": {
                "numpy": np.__version__,
                "stim": stim.__version__,
                "ldpc": importlib.metadata.version("ldpc"),
            },
            "bposd_reference": {
                "package": "stimbposd==0.1.0",
                "commit": "7921f5eb1b358ff616f9822280c9961e83df06cb",
                "constructor_argument": "error_channel=list(priors)",
            },
        }
    )

    reports: dict[str, dict[str, object]] = {}
    for p in points:
        stem = f"gatec_{basis.lower()}mem_p{p:g}_{shots}shots"
        report = run_point(
            p,
            shots,
            base_seed,
            basis=basis,
            csv_path=campaign.dir / f"{stem}_pershot.csv",
        )
        report_path = campaign.dir / f"{stem}.json"
        report_path.write_text(json.dumps(report, indent=2))
        campaign.append_result(report)
        reports[f"{p:g}"] = {
            "exact_ml": {
                "failures": int(report["ml_logical_failures"]),
                "rate": float(report["ml_ler"]),
            },
            "decoders": {
                name: {
                    "failures": int(values["logical_failures"]),
                    "rate": float(values["ler"]),
                    "ml_mismatches": int(values["ml_mismatches"]),
                }
                for name, values in report["decoders"].items()
            },
            "circuit_sha256": report["circuit_sha256"],
            "dem_sha256": report["dem_sha256"],
            "sampler_seed": int(report["sampler_seed"]),
        }

    summary = {
        "gate": "C-exact-small-BB-GB3-corrected",
        "basis": basis,
        "shots_per_point": shots,
        "points": reports,
        "bposd_prior_model": config["bposd_prior_model"],
        "supersedes_bposd_column_in": prior_campaign,
        "evidence_scope": (
            "exact/beam/NMS regression plus corrected BP comparator; "
            "BB [[36,4,4]] code-capacity only"
        ),
    }
    campaign.close(summary)
    return campaign.dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--basis", choices=("X", "Z"), required=True)
    parser.add_argument("--shots", type=int, default=2000)
    parser.add_argument("--base-seed", type=int, default=20260830)
    parser.add_argument("--points", type=float, nargs="+", default=list(POINTS))
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    shots = 20 if args.smoke and args.shots == 2000 else args.shots
    path = run_campaign(
        basis=args.basis,
        shots=shots,
        base_seed=args.base_seed,
        points=tuple(args.points),
        smoke=args.smoke,
    )
    print(path)


if __name__ == "__main__":
    main()
