#!/usr/bin/env python3
"""Recover both unambiguous comparisons after the frozen launcher key abort."""
from __future__ import annotations

from decimal import Decimal
import hashlib
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
RELEASE_ENV = "OMEGA_TWO_RUNG_RECOVERY_RELEASED"
OUTPUT = HERE / "comparison_summary_recovered.json"
EXPECTED = {
    "rung1_replay.json": "62e764f21554bc7355a45ab51bfb59690d3e9f59f2c1fc36383f861c9b9b88b5",
    "rung2_with_replay.json": "c5687a48d6a6e30038452997d1fda0b7e2a9c2259994242dce35d999a0fd6967",
    "rung2_without_replay.json": "3bd20c1bf2009e6a20e9b850bf3a50b00f7f305a386f07163cb17b6b58e079bf",
    "run_v11_replay_stdout.log": "af88e102fbc65d3d9ff1910550fa79904503417e965c1f42a983187caeda7330",
}
TARGET = Decimal("2.371866")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_decimal_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle, parse_float=Decimal)


def decimal_text(value: Decimal) -> str:
    return str(value)


def atomic_write_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    if os.path.lexists(path) or os.path.lexists(temporary):
        raise RuntimeError(f"write-once output occupied: {path} or {temporary}")
    payload = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()
    descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(descriptor, "wb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    if os.path.lexists(path):
        raise RuntimeError(f"target appeared during recovery write: {path}")
    os.replace(temporary, path)
    directory = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def main() -> None:
    if not __debug__:
        raise RuntimeError("python -O is refused")
    if os.environ.get(RELEASE_ENV) != "1":
        raise RuntimeError(f"{RELEASE_ENV}=1 required")
    actual = {name: sha256(HERE / name) for name in EXPECTED}
    if actual != EXPECTED:
        raise RuntimeError(f"frozen replay evidence mismatch: {actual}")

    abort_log = (HERE / "run_v11_replay_stdout.log").read_text(encoding="utf-8")
    if "KeyError: 'omega_cert_upper'" not in abort_log:
        raise RuntimeError("expected frozen-launcher abort is absent")

    rung1 = load_decimal_json(HERE / "rung1_replay.json")
    rung2_with = load_decimal_json(HERE / "rung2_with_replay.json")
    rung2_without = load_decimal_json(HERE / "rung2_without_replay.json")

    assert rung1["rung"] == "alman25_2.371339"
    assert rung1["params_sha256"] == "f23369136314cc51497d0c468c00d69a02b74d2b665073b9f934240ab428b46f"
    assert "omega_cert_upper" not in rung1
    assert isinstance(rung1["omega_cert_upper_raw"], Decimal)
    assert isinstance(rung1["omega_cert_with_absorbed_slack"], Decimal)

    for child in (rung2_with, rung2_without):
        assert child["rung"] == "vxxz24_2.37155181"
        assert child["params_sha256"] == "df75ae3acaa5b1388bd2f17cce266aec169d874c9fea1eb0722c3e25871da93d"
        assert isinstance(child["omega_with_absorbed_defects"], Decimal)
    assert rung2_with["include_lemma1"] is True
    assert rung2_without["include_lemma1"] is False

    r1_raw = rung1["omega_cert_upper_raw"]
    r1_absorbed = rung1["omega_cert_with_absorbed_slack"]
    r2_with = rung2_with["omega_with_absorbed_defects"]
    r2_without = rung2_without["omega_with_absorbed_defects"]

    result = {
        "campaign": HERE.name,
        "status": "RECOVERED_COMPUTATIONAL_DIAGNOSTIC_AFTER_LAUNCHER_KEY_ABORT",
        "evidence_label": "COMPUTATIONAL-EVIDENCE",
        "formal_status": "SUSPENDED",
        "formal_blockers": [
            "outward end-to-end aggregate rewrite removing decision-relevant float extraction",
            "proof that the feasibility-absorption term covers the omitted true defect",
        ],
        "launcher_abort": {
            "missing_frozen_key": "omega_cert_upper",
            "comparison_summary_json_written": False,
            "all_three_children_completed_and_parsed": True,
            "abort_log_sha256": EXPECTED["run_v11_replay_stdout.log"],
        },
        "input_sha256": actual,
        "rung1_available_values": {
            "omega_cert_upper_raw": decimal_text(r1_raw),
            "omega_cert_with_absorbed_slack": decimal_text(r1_absorbed),
        },
        "rung2_with_lemma": {
            "omega_with_absorbed_defects": decimal_text(r2_with),
        },
        "rung2_without_lemma_diagnostic": {
            "omega_with_absorbed_defects": decimal_text(r2_without),
            "with_minus_without": decimal_text(r2_with - r2_without),
        },
        "comparisons": {
            "rung1_raw_less_than_rung2_with": r1_raw < r2_with,
            "rung1_absorbed_less_than_rung2_with": r1_absorbed < r2_with,
            "rung1_raw_gap_to_rung2_with": decimal_text(r2_with - r1_raw),
            "rung1_absorbed_gap_to_rung2_with": decimal_text(r2_with - r1_absorbed),
            "rung1_raw_below_target": r1_raw < TARGET,
            "rung1_absorbed_below_target": r1_absorbed < TARGET,
            "rung2_with_below_target": r2_with < TARGET,
            "target": decimal_text(TARGET),
        },
        "interpretation": (
            "Both available rung-1 quantities are numerically below the rung-2 "
            "with-Lemma quantity and all are below the comparison target. The "
            "failed frozen key remains unresolved; this recovery chooses neither "
            "rung-1 quantity as its replacement and establishes no formal omega bound."
        ),
    }
    atomic_write_json(OUTPUT, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
