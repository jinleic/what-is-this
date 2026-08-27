#!/usr/bin/env python3
"""Independent 256-bit Arb check for the Cambie bridge point-mass margin.

This file does not import any UC campaign or certifier module. It checks only
the self-contained value

    f(t) = (1-alpha) h(2t-t^2) + alpha - h(t)

at the exact rational Campaign I values. The monotonicity argument on [psi,t],
the elementary u <= psi case, and the sequential entropy-to-UC construction
are separate human-audited steps in ``uc/AUDIT.md`` and the revised paper.
"""

from __future__ import annotations

import argparse
import datetime as _datetime
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import stat
import sys
import uuid

from flint import arb, ctx, fmpq
import flint


SCHEMA = "uc-cambie-bridge-strictness-report-v1"
PRECISION_BITS = 256
EXPECTED_FLINT_VERSION = "0.9.0"
EXPECTED_FLINT_MODULE_SHA256 = (
    "2e5f8f1768d14eccd7961353c635195bd557f297edff8b8de09e2d211f03ec2d"
)
T_EXACT = Fraction(955165028125263, 2500000000000000)
ALPHA_EXACT = Fraction(356069, 10000000)
LOWER_TEXT = (
    "0.001292783842665984988348621200827023977631257987517743548011480508341589953"
)
UPPER_TEXT = (
    "0.001292783842665984988348621200827023977631257987517743548011480508341589955"
)
HERE = Path(__file__).resolve().parent
SCHEMA_PATH = HERE / "report.schema.json"


class CheckError(RuntimeError):
    pass


def canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def object_hash(value: dict[str, object], omitted: str) -> str:
    reduced = {key: item for key, item in value.items() if key != omitted}
    return hashlib.sha256(canonical_bytes(reduced)).hexdigest()


def file_hash(path: Path) -> dict[str, object]:
    info = path.lstat()
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode):
        raise CheckError(f"not a non-symlink regular file: {path}")
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
            size += len(chunk)
    after = path.lstat()
    if (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns) != (
        after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns
    ):
        raise CheckError(f"file changed while hashing: {path}")
    return {"sha256": digest.hexdigest(), "size": size}


def exact_ball(value: Fraction) -> arb:
    return arb(fmpq(value.numerator, value.denominator))


def entropy(value: arb) -> arb:
    one = arb(1)
    if not (value > 0 and value < 1):
        raise CheckError("entropy argument was not certified inside (0,1)")
    return -(value * value.log() + (one - value) * (one - value).log()) / arb(2).log()


def run(output: Path | None) -> int:
    if not (
        sys.flags.isolated == 1
        and sys.flags.no_user_site == 1
        and sys.flags.ignore_environment == 1
        and sys.flags.dont_write_bytecode == 1
    ):
        raise CheckError("invoke with an isolated no-bytecode interpreter: python -I -B")
    version = getattr(flint, "__version__", "unknown")
    if version != EXPECTED_FLINT_VERSION:
        raise CheckError(
            f"python-flint version {version!r} differs from pinned {EXPECTED_FLINT_VERSION!r}"
        )
    module_path = Path(flint.__file__).resolve()
    module_fingerprint = file_hash(module_path)
    if module_fingerprint.get("sha256") != EXPECTED_FLINT_MODULE_SHA256:
        raise CheckError("python-flint module bytes differ from the audited environment")
    ctx.prec = PRECISION_BITS
    t = exact_ball(T_EXACT)
    alpha = exact_ball(ALPHA_EXACT)
    one = arb(1)
    q = 2 * t - t * t
    if not (t > 0 and t < 1 and alpha > 0 and alpha < 1 and q > 0 and q < 1):
        raise CheckError("exact inputs or 2t-t^2 were not certified in (0,1)")
    value = (one - alpha) * entropy(q) + alpha - entropy(t)
    lower_exact = Fraction(LOWER_TEXT)
    upper_exact = Fraction(UPPER_TEXT)
    lower = exact_ball(lower_exact)
    upper = exact_ball(upper_exact)
    if not (lower > 0 and value > lower and value < upper):
        raise CheckError(
            "Arb did not certify the pinned positive bracket for f(t): "
            + value.str(90, more=True)
        )
    psi = (arb(3) - arb(5).sqrt()) / 2
    if not t > psi:
        raise CheckError("Arb did not certify t > psi")

    # module_path and its fingerprint were pinned before arithmetic.
    report: dict[str, object] = {
        "$schema": str(SCHEMA_PATH),
        "arb_interval": value.str(90, more=True),
        "checks": {
            "f_t_gt_lower_gt_zero": True,
            "f_t_lt_upper": True,
            "t_gt_psi": True,
        },
        "claim_status": "MACHINE-VERIFIED",
        "exact_inputs": {
            "alpha": f"{ALPHA_EXACT.numerator}/{ALPHA_EXACT.denominator}",
            "t": f"{T_EXACT.numerator}/{T_EXACT.denominator}",
        },
        "finished_utc": _datetime.datetime.now(_datetime.timezone.utc).isoformat(),
        "formula": "(1-alpha)*h(2*t-t^2)+alpha-h(t)",
        "human_audited_dependencies": [
            "f'(u)<0 on [psi,t] for the point-mass function",
            "the elementary u<=psi point-mass case",
        ],
        "human_audited_dependency": (
            "Sequential entropy-to-union-closed proof in uc/AUDIT.md and uc/paper/main.tex"
        ),
        "limitations": [
            "This computation checks one exact point value only; it does not prove the derivative sign on an interval.",
            "This computation is independent of the frozen Campaign I arithmetic modules but still depends on python-flint/Arb.",
        ],
        "outcome": "PASS",
        "positive_bracket": {
            "lower_exclusive": LOWER_TEXT,
            "upper_exclusive": UPPER_TEXT,
        },
        "precision_bits": PRECISION_BITS,
        "report_type": "cambie_bridge_strictness",
        "runtime": {
            "flint_module": str(module_path),
            "flint_module_fingerprint": module_fingerprint,
            "python": sys.version,
            "python_executable": str(Path(sys.executable).resolve()),
            "python_flint": version,
        },
        "schema": SCHEMA,
        "verifier": file_hash(Path(__file__).resolve()),
    }
    report["report_sha256"] = object_hash(report, "report_sha256")
    if output is None:
        print(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=True))
    else:
        output = output.expanduser().resolve()
        output.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
        temporary = output.with_name(f".{output.name}.{uuid.uuid4().hex}.tmp")
        payload = json.dumps(report, indent=2, sort_keys=True, ensure_ascii=True).encode("ascii") + b"\n"
        descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        try:
            view = memoryview(payload)
            while view:
                written = os.write(descriptor, view)
                if written <= 0:
                    raise CheckError(f"short report write: {output}")
                view = view[written:]
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
        os.replace(temporary, output)
        print(canonical_bytes({
            "claim_status": "MACHINE-VERIFIED",
            "outcome": "PASS",
            "report": str(output),
            "report_sha256": report["report_sha256"],
        }).decode("ascii"))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    try:
        return run(parser.parse_args().output)
    except (CheckError, OSError) as exc:
        print(canonical_bytes({
            "claim_status": "FAILED", "error": str(exc), "outcome": "FAIL_CLOSED"
        }).decode("ascii"), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
