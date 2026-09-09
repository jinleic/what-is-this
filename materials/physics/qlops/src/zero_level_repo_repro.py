# -*- coding: utf-8 -*-
"""zero_level_repo_repro.py — control-plane adapter for gate
`zero-level-author-repro-r8` (physics/qlops): reproduce the c ≈ 300 zero-level
CCZ Monte-Carlo study of Itogawa et al. (arXiv:2605.21867) from the authors'
own pinned repository, without ever editing or trusting a single float of it.

Primary pin (see pre_statement.md "Revisions 5–6" for the full table):
  repo    FujitsuResearch/Zero-level_CCZ_Distillation
  commit  1b59e223590492e224bd8623a4e0bcba59029e01
  tree    9a5d89401fd40554635fb0e9d0b4da818670c2bc
  tar     codeload tar.gz, sha256 d67dbe7482b391984da5e64aeff7668bdaee45c262e6fb2dfd41dfed3d7f3338,
          size 1072583 bytes (advisory; mismatch ⇒ refusal before any write)

Scope findings preserved as preregistered facts (F-Z1/F-Z2), never silently
equated: F-Z1 — the PAPER simulates the grown variant at d=7 with 1e7–1e8
trials/point; F-Z2 — the REPOSITORY executes distance_expand=9 and embeds
published counts taken at 1e6–1e7 shots/point.  This adapter reproduces the
executable d3/d9 artifact only.  The paper's d7 claim stays NOT-REPRODUCED
regardless of the outcome here; no d7 reconstruction is attempted anywhere.

Hard rules enforced by this file:
  * The run dir must already exist at the target's exact campaigns parent and
    carry a campaign.py-minted, still-live RUNNING manifest for this gate whose
    prereg_sha256 matches the captured prereg bytes.  Frozen/closed runs are
    refused; this adapter never creates, freezes, or closes the run dir.
  * Only pinned members (matched by repo-relative path AND git blob sha1)
    are extracted from one validated in-memory capture of the source tar; the
    same bytes are retained as provenance and the source path is never reopened.
  * Postselection masks are recovered ONLY by executing the pinned driver's
    circuit-construction prefix against the pinned builder — never inferred
    from shipped .stim text.
  * The R6 ungrown p=0 smoke arm samples only the circuit rebuilt by that
    pinned driver/builder path.  It neither reads the historical
    832_text_0.stim nor makes a shipped-equality claim; builder-only recovery
    is refused for every other arm.
  * Sampling always uses circuit.compile_detector_sampler(seed=seed) with a
    fixed recorded chunk schedule (schedule is part of run identity);
    sample() itself is never given a seed.  One primary seeded run per point;
    there is no unseeded/seeded duplication and no mid-point resume.
  * detector_error_model(decompose_errors=True) always; the PyMatching graph
    must satisfy matching.num_fault_ids == circuit.num_observables (== 3).
    Accepted ⇔ no builder-postselected detector fires; error ⇔ any of the 3
    predicted observables differs from the actual ones.
  * decode_batch decodes ONLY the accepted rows and is bitwise cross-checked
    against scalar decode on the exact accepted smoke rows.  A decode_batch
    ValueError, or shipped != rebuilt (flattened) on any shipped author-sampled
    circuit, is a source/API semantic failure: REJECTED refusal before any
    scientific artifact is written.  Sole exception (Revision 8): the four
    F-Z4 split grown cells in REBUILT_ORACLE_CELLS are sampled from the
    pinned-driver rebuild itself (rebuilt oracle); their shipped artifacts
    are opened for provenance only and their factual inequality is recorded,
    never refused.

Determinism disclaimer (restated in every artifact): seeded streams are
reproducible only for this stim build, this machine, and this exact chunk
schedule.  No cross-version, cross-machine, or cross-schedule stream
identity is claimed or testable.

No top-level side effects: importing this module performs no I/O.
"""

from __future__ import annotations

import hashlib
import io
import json
import math
import os
import platform
import re
import sys
import tarfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from statistics import NormalDist

import numpy as np

# --------------------------------------------------------------- identity pins

GATE = "zero-level-author-repro-r8"
ADAPTER_VERSION = "r8"
TARGET = "physics/qlops"
ADAPTER_DIR = Path(__file__).resolve().parent
TARGET_DIR = ADAPTER_DIR.parent
PREREG_PATH = TARGET_DIR / "pre_statement.md"
TEST_PATH = ADAPTER_DIR / "test_zero_level_repo_repro.py"
PREREG_REVISION_MARKER = "## Revision 8"

# Revision 8 (F-Z4 amendment): the four grown cells whose shipped circuits
# differ from the pinned-driver rebuild by the frozen 23-hunk topology delta
# (run 20260908T174922Z_dc3b41ff_d04f0ae58ca0) are recovered with
# REBUILT-ORACLE semantics: the pinned-driver rebuild IS the reference and
# sampled circuit; the shipped artifact is demoted to provenance (blob pin,
# factual flattened inequality) and never sampled.  Every other cell keeps
# shipped-oracle semantics; inequality there remains a REJECTED refusal.
REBUILT_ORACLE_CELLS: frozenset[tuple[str, str]] = frozenset({
    ("grown", "0.0008"), ("grown", "0.0006"),
    ("grown", "0.0004"), ("grown", "0.0001"),
})

REPO = "FujitsuResearch/Zero-level_CCZ_Distillation"
COMMIT = "1b59e223590492e224bd8623a4e0bcba59029e01"
TREE_SHA = "9a5d89401fd40554635fb0e9d0b4da818670c2bc"
TAR_SHA256 = "d67dbe7482b391984da5e64aeff7668bdaee45c262e6fb2dfd41dfed3d7f3338"
TAR_SIZE = 1_072_583
REPO_ROOT_NAME = f"Zero-level_CCZ_Distillation-{COMMIT}"

REQUIRED_VERSIONS = {"stim": "1.16.0", "pymatching": "2.4.0", "numpy": "2.5.2"}

# Pinned tar members: role -> (repo-relative path, git blob sha1).
# Ticket truncated prefixes (builders f6c25a55/7d06f315, surfaces 1dc06f18/
# 3af84832, std_calc d2df580b/0e0aaf3e, plot 35d58268, LICENSE b9d9364d) are
# asserted as prefixes in the test file.
MEMBER_PINS: dict[str, tuple[str, str]] = {
    "driver_ungrown": ("stim/ungrown/4_surface3_d3.py",
                       "6205e13302cbe0285f3e704b9a66da1df3ee6173"),
    "driver_grown": ("stim/grown/3_surface3_d3_expand.py",
                     "592d1a8422ae2aad47b37d66e694c4eeb2b0f019"),
    "builder_ungrown": ("stim/ungrown/stim_builder.py",
                        "f6c25a555c99ff69ec251fef42e33ba4566621d0"),
    "builder_grown": ("stim/grown/stim_builder.py",
                      "7d06f3154130fab26850028c056ea2734a24be77"),
    "surface_ungrown": ("stim/ungrown/surface_func.py",
                        "1dc06f181feba621f8451750dda9cd34897c9dcb"),
    "surface_grown": ("stim/grown/surface_func.py",
                      "3af84832018afc4959c646f077326c88769332d3"),
    "std_calc_ungrown": ("stim/ungrown/std_calc.py",
                         "d2df580b8e51ac41b383ed0e8ad4b0be1978a91e"),
    "std_calc_grown": ("stim/grown/std_calc.py",
                       "0e0aaf3e54a82eb5a311eac51b975debe2409279"),
    "plot": ("plot_LER_suc/CCZ_plot.py",
             "35d582688622481aa2fce3f76e22c1bae3b8b0d2"),
    "license": ("LICENSE.txt",
                "b9d9364d999630ddd675f259576e6ed66711f3e5"),
    # Shipped circuit text.  The 12 named p points are author-sampled oracles;
    # ungrown p=0 remains pinned historical source but is not an R6 oracle.
    "shipped_ungrown_0": ("stim/ungrown/832_text_0.stim",
                          "5ffcb9fa182bdc3cfac648272caff1ca4b9ad1dd"),
    "shipped_ungrown_0.0001": ("stim/ungrown/832_text_0.0001.stim",
                               "3b8113e2c31d262555dda61476908b524db67ef0"),
    "shipped_ungrown_0.0002": ("stim/ungrown/832_text_0.0002.stim",
                               "7dbad92ef8efe0d83fc9e443003067a54a342dcd"),
    "shipped_ungrown_0.0004": ("stim/ungrown/832_text_0.0004.stim",
                               "fc600e2bf4b578c3adad6d876f03daad1728b5d3"),
    "shipped_ungrown_0.0006": ("stim/ungrown/832_text_0.0006.stim",
                               "48dea0bd09dc9dc149ce841fb0f17a30d0639679"),
    "shipped_ungrown_0.0008": ("stim/ungrown/832_text_0.0008.stim",
                               "bd3b0682760a417c11d66231dd55bce92d91321b"),
    "shipped_ungrown_0.001": ("stim/ungrown/832_text_0.001.stim",
                              "782579515559fba45dea0b362167f81d5e26bbd5"),
    "shipped_grown_0.0001": ("stim/grown/832_text_0.0001.stim",
                             "cecca6083c5575da45071615fcef6a5001356804"),
    "shipped_grown_0.0002": ("stim/grown/832_text_0.0002.stim",
                             "a12007a39c2622ba0ab745a14544ce4cef30a394"),
    "shipped_grown_0.0004": ("stim/grown/832_text_0.0004.stim",
                             "3a28a475a9f8893e55638df97ccebba96c81125c"),
    "shipped_grown_0.0006": ("stim/grown/832_text_0.0006.stim",
                             "ea17f178c42b46c117d6826ddd87c14b861fdba7"),
    "shipped_grown_0.0008": ("stim/grown/832_text_0.0008.stim",
                             "27fa0d4283ee05c991098e0d4d3a9efc3854d274"),
    "shipped_grown_0.001": ("stim/grown/832_text_0.001.stim",
                            "3f3f7e9f4b58f3be7c769120948c970f217f8f13"),
}

VARIANTS = ("ungrown", "grown")
P_LABELS = ("0.001", "0.0008", "0.0006", "0.0004", "0.0002", "0.0001")

# Author rows embedded verbatim in the pinned sources (driver docstrings,
# cross-checked byte-for-byte against plot_LER_suc/CCZ_plot.py):
# (p, shots, logical error rate, acceptance).
AUTHOR_ROWS: dict[str, tuple[tuple[str, int, float, float], ...]] = {
    "ungrown": (
        ("0.001", 5_000_000, 0.0002703683888719046, 0.418688),
        ("0.0008", 5_000_000, 0.00019523066744384006, 0.4978726),
        ("0.0006", 5_000_000, 0.00010763579824934929, 0.5927396),
        ("0.0004", 5_000_000, 4.65009456137416e-05, 0.705362),
        ("0.0002", 5_000_000, 1.1667975146736694e-05, 0.8399058),
        ("0.0001", 10_000_000, 3.1638657726672983e-06, 0.9166002),
    ),
    "grown": (
        ("0.001", 1_000_000, 0.00034852566383247535, 0.275446),
        ("0.0008", 5_000_000, 0.0002219381364258208, 0.3568562),
        ("0.0006", 5_000_000, 0.00011960131162771752, 0.4615334),
        ("0.0004", 5_000_000, 5.422499759000011e-05, 0.5975104),
        ("0.0002", 10_000_000, 1.3586052894515038e-05, 0.7728514),
        ("0.0001", 10_000_000, 3.866977783985162e-06, 0.8792396),
    ),
}
AUTHOR_TOTAL_SHOTS = sum(r[1] for rows in AUTHOR_ROWS.values() for r in rows)
# Executable/authoritative schedule: ungrown 5×5e6 + 1e7 (p=1e-4), grown
# 1e6 + 3×5e6 + 2×1e7 — 71,000,000 author shots total.  SOURCE CONFLICT
# (F-Z3): the ungrown driver's embedded comment says 5e6 for p=1e-4, but the
# printed LER is only consistent with 1e7 shots (29 errors / 9,166,002
# accepted); the row is treated as 1e7-shot and the conflict is named.
TICKET_BUDGET_SHOTS = 71_000_000
assert AUTHOR_TOTAL_SHOTS == TICKET_BUDGET_SHOTS

# Authors' unweighted-through-origin c fits (== least squares of LER = c*p^2,
# equivalent to their scipy curve_fit on a*x**2 with no weights); recomputed
# from the embedded rows and asserted by the test file.
AUTHOR_FITS = {"ungrown": 282.1598546872, "grown": 346.5564346243}
# F-Z3 audit record: per-point source conflicts named, never hidden.
AUTHOR_ROW_SOURCE_CONFLICTS: dict[tuple[str, str], dict] = {
    ("ungrown", "0.0001"): {
        "comment_says_shots": 5_000_000,
        "arithmetic_implies_shots": 10_000_000,
        "implied_accepted": 9_166_002,
        "implied_errors": 29,
        "nominal_reading_errors": 14.5,
        "resolution": ("treat the row as 1e7 shots for the executable "
                       "schedule and both author sigmas; the comment/row "
                       "conflict is recorded here and in pre_statement.md "
                       "Revision 5"),
    },
}


# Author analytic sigmas (pinned std_calc.py, both variants):
#   LER:        sqrt(ler*(1-ler) / (shots*acceptance))
#   acceptance: sqrt(acc*(1-acc) / shots)


def ler_sigma(shots: int, ler: float, acceptance: float) -> float:
    if shots <= 0 or acceptance <= 0:
        return math.inf
    return math.sqrt(ler * (1.0 - ler) / (shots * acceptance))


def acc_sigma(shots: int, acceptance: float) -> float:
    if shots <= 0:
        return math.inf
    return math.sqrt(acceptance * (1.0 - acceptance) / shots)


# Fixed seeds (preregistered; part of run identity).  Smoke repeats reuse the
# same seed — that identity is the point of the repeat gate.
FULL_SEEDS: dict[str, dict[str, int]] = {
    "ungrown": {"0.001": 1101, "0.0008": 1102, "0.0006": 1103,
                "0.0004": 1104, "0.0002": 1105, "0.0001": 1106},
    "grown": {"0.001": 2201, "0.0008": 2202, "0.0006": 2203,
              "0.0004": 2204, "0.0002": 2205, "0.0001": 2206},
}
SMOKE_SEEDS = {"ungrown_noiseless": 3301, "grown_0.001": 3302}

CHUNK_SMOKE = 5_000
CHUNK_FULL = 10_000
SMOKE_NOISELESS_SHOTS = 10_000
SMOKE_GROWN_SHOTS = 20_000

# Familywise 1% over all 24 comparisons (12 LER + 12 acceptance).  The 24
# statistics share runs, so the dependence-valid Bonferroni bound is the
# preregistered derivation: exact z = Φ⁻¹(1 - 0.01/48) = 3.529296088834735,
# rounded upward to 3.53.  (Exact Šidák assumes independence:
# z = 3.528022653482286 — it rounds up to the same gate; 3.53 is used.)
SIDAK_FAMILY = 24
SIDAK_ALPHA = 0.01
DECISIVE_Z = 5.0


def _two_sided_z(alpha_per_comparison: float) -> float:
    return NormalDist().inv_cdf(1.0 - alpha_per_comparison / 2.0)


def bonferroni_z_exact(family: int = SIDAK_FAMILY,
                       alpha: float = SIDAK_ALPHA) -> float:
    return _two_sided_z(alpha / family)


def sidak_z_exact(family: int = SIDAK_FAMILY,
                  alpha: float = SIDAK_ALPHA) -> float:
    per = 1.0 - (1.0 - alpha) ** (1.0 / family)
    return _two_sided_z(per)


Z_THRESHOLD = 3.53
assert math.ceil(bonferroni_z_exact() * 100.0) / 100.0 == Z_THRESHOLD
assert math.ceil(sidak_z_exact() * 100.0) / 100.0 == Z_THRESHOLD

STREAM_NOTE = (
    "Seeded detector streams are deterministic only for this stim build, this "
    "machine, and this exact chunk schedule. No stream identity is claimed or "
    "implied across stim versions, machines, or schedules.")

# Driver prefix handling: execute the pinned construction code only.
BUILD_LINE = "circuit = circuit_builder.build()"
ERROR_RATE_RE = re.compile(r"^error_rate = \[[^\]\n]*\]\s*$", re.MULTILINE)
# The prefix must never touch the filesystem or sample.
PREFIX_FORBIDDEN = ("open(", "output_svg", "output_circuit_text",
                    "compile_detector_sampler", ".sample(")


class Refusal(Exception):
    """Pre-sampling refusal (REJECTED path); CLI maps to exit code 1."""


# ------------------------------------------------------------------ hash utils

def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def git_blob_sha1(data: bytes) -> str:
    """Git object name of a blob with exactly this content."""
    header = f"blob {len(data)}".encode("ascii") + b"\x00"
    return hashlib.sha1(header + data).hexdigest()


# -------------------------------------------------------------- atomic writers

def _fsync_dir(path: Path) -> None:
    try:
        fd = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _tmp_name(dest: Path) -> Path:
    return dest.parent / f".{dest.name}.{os.getpid()}.{os.urandom(4).hex()}.tmp"


def atomic_write_bytes(dest: Path, data: bytes) -> None:
    tmp = _tmp_name(dest)
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(fd, "wb") as fh:
        fh.write(data)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, dest)
    _fsync_dir(dest.parent)


def atomic_write_json(dest: Path, payload: object) -> None:
    atomic_write_bytes(dest,
                       (json.dumps(payload, indent=2, sort_keys=True)
                        + "\n").encode("utf-8"))


# ------------------------------------------------------------ environment gate

def validate_env() -> dict:
    import stim
    import pymatching
    versions = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "stim": stim.__version__,
        "pymatching": pymatching.__version__,
        "numpy": np.__version__,
    }
    for pkg, want in REQUIRED_VERSIONS.items():
        got = versions[pkg]
        if got != want:
            raise Refusal(
                f"environment pin violated: {pkg}=={want} required, found {got}; "
                "seeded-stream identity is version-bound, refusing")
    return versions
def refuse_tar_inside_run_dir(tar_path: Path, run_dir: Path) -> None:
    if Path(tar_path).resolve().is_relative_to(Path(run_dir).resolve()):
        raise Refusal(
            "source tar lives inside the run dir; refusing (provenance must "
            "not be nested where run artifacts are written)")




# ------------------------------------------------------------- source tar gate

def validate_source_tar(path: Path) -> dict:
    """Advisory hash + size gate. Runs before anything is written anywhere."""
    if not path.is_file():
        raise Refusal(f"source tar absent: {path}")
    size = path.stat().st_size
    if size != TAR_SIZE:
        raise Refusal(f"source tar size {size} != pinned {TAR_SIZE}: {path}")
    digest = sha256_file(path)
    if digest != TAR_SHA256:
        raise Refusal(
            f"source tar sha256 {digest} != pinned {TAR_SHA256}: {path}; "
            "refusing to substitute an unpinned artifact")
    return {"path": str(path), "size": size, "sha256": digest}


def _check_member_name(name: str) -> None:
    pure = PurePosixPath(name)
    if pure.is_absolute() or ".." in pure.parts or "\\" in name:
        raise Refusal(f"tar member name rejected (traversal): {name!r}")


def _repo_rel(member_name: str) -> str | None:
    parts = PurePosixPath(member_name).parts
    if len(parts) < 2 or parts[0] != REPO_ROOT_NAME:
        return None
    return "/".join(parts[1:])


def extract_pinned_members(tar_source, dest_root: Path,
                           pins: dict[str, tuple[str, str]]) -> dict:
    """Scan the whole tar, validate every pin, then (and only then) write.

    tar_source may be a Path (its pristine hash is recorded and re-checked
    after extraction) or a binary file object (tests).  Writes are idempotent:
    an existing destination with identical bytes is kept, with different
    bytes it is a refusal.  Nothing outside dest_root is ever written.
    """
    path_pins = {rel: role for role, (rel, _blob) in pins.items()}
    blob_pins = {blob: role for role, (_rel, blob) in pins.items()}
    pristine = None
    if isinstance(tar_source, (str, Path)):
        tar_path = Path(tar_source)
        pristine = {"sha256": sha256_file(tar_path),
                    "size": tar_path.stat().st_size}
        source: object = tar_path
    else:
        source = tar_source

    found: dict[str, tuple[str, bytes]] = {}
    seen_rels: set[str] = set()
    ignored = 0
    tar_context = (
        tarfile.open(source, mode="r:*")
        if isinstance(source, (str, Path))
        else tarfile.open(fileobj=source, mode="r:*")
    )
    with tar_context as tf:
        while True:
            member = tf.next()
            if member is None:
                break
            _check_member_name(member.name)
            if member.isdir():
                continue
            if not member.isfile():
                raise Refusal(
                    f"tar member is not a regular file: {member.name!r}")
            rel = _repo_rel(member.name)
            if rel is None or rel not in path_pins:
                ignored += 1
                continue
            if rel in seen_rels:
                raise Refusal(f"duplicate pinned member in tar: {rel}")
            seen_rels.add(rel)
            data = tf.extractfile(member).read()
            blob = git_blob_sha1(data)
            role = path_pins[rel]
            expected_blob = pins[role][1]
            if blob != expected_blob:
                raise Refusal(
                    f"pinned member {rel}: git blob sha1 {blob} != pinned "
                    f"{expected_blob}")
            owner = blob_pins.get(blob, role)
            if owner != role:
                raise Refusal(
                    f"content of {rel} matches pin of {owner}: mislabeled tar")
            found[role] = (rel, data)

    missing = sorted(set(pins) - set(found))
    if missing:
        raise Refusal(f"pinned members absent from tar: {missing}")

    dest_root.mkdir(parents=True, exist_ok=True)
    written, kept = [], []
    for role in sorted(found):
        rel, data = found[role]
        dest = dest_root / REPO_ROOT_NAME / rel
        if dest.exists():
            if sha256_file(dest) == sha256_bytes(data):
                kept.append({"role": role, "path": rel,
                             "sha256": sha256_bytes(data),
                             "bytes": len(data)})
                continue
            raise Refusal(
                f"refusing to overwrite differing run artifact: {dest}")
        dest.parent.mkdir(parents=True, exist_ok=True)
        atomic_write_bytes(dest, data)
        written.append({"role": role, "path": rel,
                        "sha256": sha256_bytes(data), "bytes": len(data)})

    recheck = None
    if pristine is not None:
        tar_path = Path(tar_source)
        recheck = {"sha256": sha256_file(tar_path),
                   "size": tar_path.stat().st_size}
        if recheck != pristine:
            raise Refusal("source tar changed during extraction")

    return {"pristine": pristine, "recheck": recheck,
            "written": written, "kept_idempotent": kept,
            "ignored_members": ignored,
            "pins": {role: {"path": pins[role][0], "blob_sha1": pins[role][1],
                            "sha256": sha256_bytes(found[role][1]),
                            "bytes": len(found[role][1])}
                     for role in sorted(found)}}


# --------------------------------------------------------- campaign run gate

MANIFEST_KEYS = {"run_id", "target", "gate", "agent", "created_utc",
                 "prereg_sha256", "status"}
RUN_ID_RE = re.compile(
    r"^(?P<stamp>\d{8}T\d{6}Z)_(?P<uuid8>[0-9a-f]{8})_"
    r"(?P<hash12>[0-9a-f]{12})$")
TS_Z = "%Y-%m-%dT%H:%M:%SZ"


def assert_run_live(run_dir: Path, captured_manifest: bytes) -> None:
    """Refuse if campaign.py froze/closed the run or changed its manifest."""
    manifest_path = Path(run_dir) / "manifest.json"
    try:
        live_manifest = manifest_path.read_bytes()
    except OSError as exc:
        raise Refusal(f"campaign manifest unreadable at recheck: {exc}") from exc
    if live_manifest != captured_manifest:
        raise Refusal("campaign manifest bytes changed since validation; "
                      "refusing to continue")
    for marker in ("status.json", "sha256s.txt"):
        if (Path(run_dir) / marker).exists():
            raise Refusal(
                f"campaign is no longer live: {marker} exists; refusing to "
                "accept or write run evidence")


def validate_run_dir(run_dir: Path, prereg_source, *,
                     manifest_bytes: bytes | None = None) -> dict:
    """Validate an exact, live campaign.py run; never create or close it.

    prereg_source is the preregistration CAPTURED ONCE at prepare start
    (bytes; a Path is accepted for tests and is read once here).  prepare_run
    also supplies its one captured manifest byte string for later live checks.
    """
    run_dir = Path(run_dir).resolve()
    if not run_dir.is_dir():
        raise Refusal(f"run dir does not exist (never created by adapter): "
                      f"{run_dir}")
    expected_parent = (TARGET_DIR / "campaigns").resolve()
    if run_dir.parent.resolve() != expected_parent:
        raise Refusal(
            f"run dir parent {run_dir.parent.resolve()} != exact campaigns "
            f"parent {expected_parent}")
    manifest_path = run_dir / "manifest.json"
    if not manifest_path.is_file():
        raise Refusal(f"missing campaign manifest: {manifest_path}")
    if manifest_bytes is None:
        try:
            manifest_bytes = manifest_path.read_bytes()
        except OSError as exc:
            raise Refusal(f"manifest unreadable {manifest_path}: {exc}") from exc
    elif not isinstance(manifest_bytes, bytes):
        raise Refusal("captured campaign manifest must be bytes")
    try:
        manifest = json.loads(manifest_bytes.decode("utf-8"))
    except (UnicodeError, ValueError) as exc:
        raise Refusal(f"unparsable manifest {manifest_path}: {exc}") from exc
    if not isinstance(manifest, dict) or set(manifest) != MANIFEST_KEYS:
        keys = sorted(manifest) if isinstance(manifest, dict) else \
            type(manifest).__name__
        raise Refusal(
            f"manifest keys {keys} != canonical {sorted(MANIFEST_KEYS)}")
    if manifest["gate"] != GATE:
        raise Refusal(f"manifest gate {manifest['gate']!r} != {GATE!r}")
    if manifest["target"] != TARGET:
        raise Refusal(f"manifest target {manifest['target']!r} != {TARGET!r}")
    if manifest["status"] != "RUNNING":
        raise Refusal(
            f"manifest status {manifest['status']!r} != 'RUNNING' (a live "
            "campaign.py init is required; frozen/closed runs are refused)")
    run_id = manifest["run_id"]
    match = RUN_ID_RE.fullmatch(run_id) if isinstance(run_id, str) else None
    if match is None:
        raise Refusal(f"manifest run_id malformed: {run_id!r}")
    if run_dir.name != run_id:
        raise Refusal(
            f"run dir basename {run_dir.name!r} != manifest run_id {run_id!r}")
    agent = manifest["agent"]
    if not isinstance(agent, str) or not agent:
        raise Refusal(f"manifest agent malformed: {agent!r}")
    created_utc = manifest["created_utc"]
    try:
        created = datetime.strptime(created_utc, TS_Z)
    except (TypeError, ValueError) as exc:
        raise Refusal(
            f"manifest created_utc malformed: {created_utc!r}") from exc
    if created.strftime(TS_Z) != created_utc:
        raise Refusal(f"manifest created_utc is not canonical: {created_utc!r}")
    stamp = match.group("stamp")
    if created.strftime("%Y%m%dT%H%M%SZ") != stamp:
        raise Refusal(
            f"manifest created_utc {created_utc!r} does not match run_id "
            f"stamp {stamp!r}")
    if isinstance(prereg_source, (bytes, bytearray)):
        prereg_bytes = bytes(prereg_source)
    else:
        prereg_path = Path(prereg_source)
        if not prereg_path.is_file():
            raise Refusal(f"prereg absent: {prereg_path}")
        prereg_bytes = prereg_path.read_bytes()
    captured_sha = sha256_bytes(prereg_bytes)
    if manifest["prereg_sha256"] != captured_sha:
        raise Refusal(
            "manifest prereg_sha256 does not match captured prereg bytes "
            f"({manifest['prereg_sha256']} != {captured_sha}); the campaign "
            "must be re-init'd against the current pre_statement.md")
    binder = "\x1f".join(
        (GATE, agent, captured_sha, stamp, match.group("uuid8")))
    expected_hash12 = sha256_bytes(binder.encode("utf-8"))[:12]
    if match.group("hash12") != expected_hash12:
        raise Refusal(
            f"manifest run_id hash12 {match.group('hash12')!r} does not bind "
            "gate/agent/prereg/stamp/uuid8 as minted by campaign.py")
    try:
        prereg_text = prereg_bytes.decode("utf-8")
    except UnicodeError as exc:
        raise Refusal(f"captured prereg is not UTF-8: {exc}") from exc
    if PREREG_REVISION_MARKER not in prereg_text:
        raise Refusal(
            f"prereg lacks {PREREG_REVISION_MARKER}; refusing to run against "
            "an out-of-date preregistration")
    revisions = [int(n) for n in re.findall(
        r"^## Revision ([1-9][0-9]*)\b", prereg_text, re.MULTILINE)]
    if 8 not in revisions:
        raise Refusal("prereg has no '## Revision 8' heading")
    if revisions != sorted(set(revisions)):
        raise Refusal(
            f"Revision headings not strictly increasing: {revisions}")
    if revisions[-1] != 8:
        raise Refusal(f"latest prereg revision is {revisions[-1]}; expected "
                      "exactly 8 (no later revision may exist)")
    assert_run_live(run_dir, manifest_bytes)
    return manifest


def assert_prereg_unchanged(prereg_path: Path, captured: bytes) -> None:
    """TOCTOU guard: the live prereg must still equal the validated bytes."""
    try:
        live = Path(prereg_path).read_bytes()
    except OSError as exc:
        raise Refusal(f"prereg unreadable at recheck: {exc}") from exc
    if live != captured:
        raise Refusal("prereg bytes changed since validation; refusing to "
                      "continue scientific work (TOCTOU guard)")


# --------------------------------------------- pinned driver execution (masks)

def split_driver_prefix(driver_text: str) -> str:
    """Everything up to and including the single build() line, verbatim."""
    lines = driver_text.split("\n")
    hits = [i for i, ln in enumerate(lines) if ln.strip() == BUILD_LINE]
    if len(hits) != 1:
        raise Refusal(
            f"expected exactly one {BUILD_LINE!r} line, found {len(hits)}; "
            "driver layout deviates from pin")
    prefix = "\n".join(lines[:hits[0] + 1]) + "\n"
    for token in PREFIX_FORBIDDEN:
        if token in prefix:
            raise Refusal(
                f"driver prefix contains forbidden token {token!r}; refusing "
                "to execute file-writing or sampling code")
    return prefix


def parameterize_error_rate(prefix_text: str, p: float) -> str:
    """Rewrite ONLY the driver's own per-run `error_rate = [...]` assignment.

    The authors operated the repo by editing this one line per point; the
    construction code around it is executed byte-identically.
    """
    hits = ERROR_RATE_RE.findall(prefix_text)
    if len(hits) != 1:
        raise Refusal(
            f"expected exactly one error_rate assignment, found {len(hits)}")
    return ERROR_RATE_RE.sub(f"error_rate = [{p!r}]", prefix_text, count=1)


def _exec_driver_prefix(prefix_text: str, variant_dir: Path, p: float):
    """Exec the parameterized prefix with pinned-module import hygiene."""
    text = parameterize_error_rate(prefix_text, p)
    dir_str = str(variant_dir)
    baseline = frozenset(sys.modules)
    previous_dont_write_bytecode = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    sys.path.insert(0, dir_str)
    try:
        namespace: dict = {"__name__": "zero_level_pinned_driver_prefix"}
        exec(compile(text, "<pinned-driver-prefix>", "exec"), namespace)  # noqa: S102 - pinned, hashed source
        builder = namespace.get("circuit_builder")
        built = namespace.get("circuit")
        if builder is None or built is None:
            raise Refusal("driver prefix did not define the expected objects")
        for mod_name in ("stim_builder", "surface_func"):
            mod = sys.modules.get(mod_name)
            mod_file = getattr(mod, "__file__", "") or ""
            if not mod_file or not Path(mod_file).resolve().is_relative_to(
                    variant_dir.resolve()):
                raise Refusal(
                    f"imported {mod_name} from {mod_file or '<builtin>'}, not "
                    f"the pinned extract under {variant_dir}")
    finally:
        sys.dont_write_bytecode = previous_dont_write_bytecode
        sys.path.remove(dir_str)
        for mod_name in ("stim_builder", "surface_func"):
            if mod_name not in baseline:
                sys.modules.pop(mod_name, None)
    return builder, built, sha256_bytes(text.encode("utf-8"))

def shipped_filename(p_label: str) -> str:
    return f"832_text_{p_label}.stim"


class RecoveredPoint:
    """Circuit, builder mask, matcher, rebuilt circuit, and provenance.

    ``flattened_equal`` is ``None`` only for the tightly restricted
    builder-only ungrown p=0 smoke arm, where no shipped oracle is read.
    """

    def __init__(self, variant: str, p_label: str, circuit, mask: list[int],
                 matcher, rebuilt, flattened_equal: bool | None,
                 provenance: dict):
        self.variant = variant
        self.p_label = p_label
        self.circuit = circuit          # the circuit actually sampled
        self.mask = mask                # builder postselection detector ids
        self.matcher = matcher          # None for the noiseless p=0 arm
        self.rebuilt = rebuilt
        self.flattened_equal = flattened_equal
        self.provenance = provenance


def recover_point(source_root: Path, variant: str, p_label: str, *,
                  builder_only: bool = False,
                  rebuilt_oracle: bool = False) -> RecoveredPoint:
    """Recover a pinned circuit and builder mask with explicit provenance.

    Default recovery reads the shipped circuit and refuses flattened
    inequality here, before DEM construction or sampling.  Builder-only
    recovery is permitted solely for the ungrown p=0 smoke arm; it samples
    the rebuilt circuit without opening or making any claim about
    ``832_text_0.stim``.  Rebuilt-oracle recovery (Revision 8) is permitted
    solely for the F-Z4 split cells in ``REBUILT_ORACLE_CELLS``: it opens
    the shipped artifact for provenance only (blob pin and factual
    flattened inequality are recorded), never refuses on that inequality,
    and samples the pinned-driver rebuild, which IS the reference circuit
    for those cells.
    """
    if builder_only and rebuilt_oracle:
        raise Refusal("builder_only and rebuilt_oracle are mutually "
                      "exclusive recovery modes")
    if builder_only and (variant, p_label) != ("ungrown", "0"):
        raise Refusal(
            "builder-only recovery is restricted to the ungrown p=0 smoke "
            "arm; shipped equality remains mandatory for author-sampled arms")
    if rebuilt_oracle and (variant, p_label) not in REBUILT_ORACLE_CELLS:
        raise Refusal(
            "rebuilt-oracle recovery is restricted to the Revision-8 F-Z4 "
            f"split cells {sorted(REBUILT_ORACLE_CELLS)}; "
            f"{variant}@{p_label} keeps shipped-oracle semantics")

    import stim
    import pymatching

    p = float(p_label)
    vdir = source_root / REPO_ROOT_NAME / "stim" / variant
    driver_role = f"driver_{variant}"
    driver_rel, _blob = MEMBER_PINS[driver_role]
    driver_text = (vdir / Path(driver_rel).name).read_text(encoding="utf-8")
    prefix = split_driver_prefix(driver_text)
    builder, built, prefix_exec_sha = _exec_driver_prefix(prefix, vdir, p)

    if not isinstance(built, list) or len(built) != 1:
        raise Refusal(
            f"{variant}: build() returned {type(built).__name__} len="
            f"{len(built) if isinstance(built, list) else '-'}, expected a "
            "single-element list (one error rate)")
    rebuilt = built[0]

    postselect = builder.postselct_numbers()
    mask = [i for i, flag in enumerate(postselect) if flag]

    circuit = rebuilt
    flattened_equal: bool | None = None
    sampled_source = "rebuilt_from_pinned_driver_builder"
    shipped_path = None
    if not builder_only:
        shipped_path = vdir / shipped_filename(p_label)
        shipped = stim.Circuit(shipped_path.read_text(encoding="utf-8"))
        flattened_equal = rebuilt.flattened() == shipped.flattened()
        if not flattened_equal and not rebuilt_oracle:
            raise Refusal(
                f"{variant}@{p_label}: shipped != rebuilt (flattened()); "
                "source semantic failure — refusing before sampling")
        if not rebuilt_oracle:
            circuit = shipped
            sampled_source = "shipped_author_circuit"

    if circuit.num_observables != 3 or rebuilt.num_observables != 3:
        raise Refusal(
            f"{variant}@{p_label}: expected 3 observables, got sampled="
            f"{circuit.num_observables}, rebuilt={rebuilt.num_observables}")
    if not (circuit.num_detectors == rebuilt.num_detectors
            == len(postselect)):
        raise Refusal(
            f"{variant}@{p_label}: detector count mismatch sampled="
            f"{circuit.num_detectors} rebuilt={rebuilt.num_detectors} "
            f"mask={len(postselect)}")

    dem = circuit.detector_error_model(decompose_errors=True)
    dem_sha = sha256_bytes(str(dem).encode("utf-8"))
    matcher = None
    n_error_instr = sum(1 for instr in dem if instr.type == "error")
    if p == 0.0:
        if n_error_instr != 0:
            raise Refusal(
                f"{variant}@{p_label}: noiseless circuit carries "
                f"{n_error_instr} DEM error mechanisms")
    else:
        matcher = pymatching.Matching.from_detector_error_model(dem)
        if matcher.num_fault_ids != circuit.num_observables:
            raise Refusal(
                f"{variant}@{p_label}: matching.num_fault_ids "
                f"({matcher.num_fault_ids}) != num_observables "
                f"({circuit.num_observables})")

    recovery_mode = ("builder_only_p0" if builder_only
                     else "rebuilt_oracle_r8" if rebuilt_oracle
                     else "shipped_equality_required")
    provenance = {
        "variant": variant, "p": p, "p_label": p_label,
        "recovery_mode": recovery_mode,
        "sampled_circuit_source": sampled_source,
        "dem_source": sampled_source,
        "shipped_oracle_used": not builder_only and not rebuilt_oracle,
        "driver_role": driver_role, "driver_blob_sha1": _blob_of(driver_role),
        "builder_role": f"builder_{variant}",
        "builder_blob_sha1": _blob_of(f"builder_{variant}"),
        "surface_role": f"surface_{variant}",
        "surface_blob_sha1": _blob_of(f"surface_{variant}"),
        "prefix_exec_sha256": prefix_exec_sha,
        "mask": mask,
        "mask_sha256": sha256_bytes(
            json.dumps(mask, separators=(",", ":")).encode("ascii")),
        "num_detectors": circuit.num_detectors,
        "num_observables": circuit.num_observables,
        "dem_sha256": dem_sha,
        "dem_error_instructions": n_error_instr,
        "shipped_rebuilt_flattened_equal": flattened_equal,
        "mask_source": "executed pinned driver prefix + pinned stim_builder"
                       ".postselct_numbers(); never inferred from .stim",
    }
    if not builder_only:
        provenance.update({
            "shipped_path": str(shipped_path),
            "shipped_blob_sha1": _blob_of(
                f"shipped_{variant}_{p_label}"),
        })
    return RecoveredPoint(variant, p_label, circuit, mask, matcher, rebuilt,
                          flattened_equal, provenance)


def _blob_of(role: str) -> str:
    return MEMBER_PINS[role][1]


# ------------------------------------------------------------------- sampling

def chunk_schedule(total_shots: int, chunk: int) -> dict:
    if total_shots <= 0 or chunk <= 0:
        raise Refusal(f"bad schedule: total={total_shots} chunk={chunk}")
    num_chunks = math.ceil(total_shots / chunk)
    final = total_shots - chunk * (num_chunks - 1)
    return {"total": total_shots, "chunk": chunk,
            "num_chunks": num_chunks, "final_chunk": final,
            "rule": "every chunk == 'chunk' except the final remainder"}


def _chunk_digest(syndrome, actual) -> str:
    h = hashlib.sha256()
    h.update(np.ascontiguousarray(syndrome).tobytes())
    h.update(np.ascontiguousarray(actual).tobytes())
    return h.hexdigest()


def sample_counts(circuit, mask: list[int], matcher, shots: int, chunk: int,
                  seed: int, crosscheck_accepted_rows: bool = False,
                  _sampler_factory=None) -> dict:
    """One primary seeded run over a fixed chunk schedule.

    The seed is bound at compile_detector_sampler; sample() is never seeded.
    decode_batch is used for all predictions; when requested, it is bitwise
    cross-checked against scalar decode on exactly the accepted rows.
    """
    if _sampler_factory is None:  # test seam; production path below
        _sampler_factory = lambda circuit_, seed_: (
            circuit_.compile_detector_sampler(seed=seed_))
    sampler = _sampler_factory(circuit, int(seed))
    sched = chunk_schedule(shots, chunk)
    accepted = errors = discarded = 0
    chunk_records = []
    crosscheck = None
    if crosscheck_accepted_rows:
        crosscheck = {"rows": 0, "mismatches": 0}
    noiseless = matcher is None
    noiseless_detector_bits_fired = 0
    noiseless_actual_flip_rows = 0

    for i in range(sched["num_chunks"]):
        n = sched["chunk"] if i < sched["num_chunks"] - 1 else sched["final_chunk"]
        syndrome, actual = sampler.sample(n, separate_observables=True)
        record = {"n": int(n),
                  "syndrome_obs_sha256": _chunk_digest(syndrome, actual)}
        if noiseless:
            fired = int(np.count_nonzero(syndrome))
            noiseless_detector_bits_fired += fired
            actual_flip_rows = int(np.count_nonzero(np.any(actual, axis=1)))
            noiseless_actual_flip_rows += actual_flip_rows
            errors += actual_flip_rows
            record.update({"accepted": n, "errors": actual_flip_rows,
                           "detector_bits_fired": fired,
                           "actual_observable_flip_rows": actual_flip_rows})
            if fired != 0:
                raise Refusal(
                    "noiseless arm produced nonzero detectors "
                    f"({fired} bits in chunk {i}); circuit contradicts p=0")
            if actual_flip_rows != 0:
                raise Refusal(
                    "noiseless arm produced nonzero actual-observable flips "
                    f"({actual_flip_rows} rows in chunk {i}); circuit "
                    "contradicts p=0")
            accepted += n
        else:
            # Upstream semantics: rejected shots are skipped BEFORE decode.
            # Only accepted syndromes reach the matcher; rejected rows can
            # neither cause API failure nor spend decoder work.
            fired = np.any(syndrome[:, mask], axis=1) if mask else \
                np.zeros(n, dtype=bool)
            acc_rows = ~fired
            n_acc = int(np.count_nonzero(acc_rows))
            n_err = 0
            if n_acc:
                acc_syndrome = np.ascontiguousarray(syndrome[acc_rows])
                acc_actual = actual[acc_rows]
                try:
                    prediction = matcher.decode_batch(acc_syndrome)
                except ValueError as exc:
                    raise Refusal(
                        f"decode_batch ValueError (chunk {i}): {exc}; "
                        "decoder API semantic failure — refusing before "
                        "any scientific artifact") from exc
                n_err = int(np.count_nonzero(
                    np.any(prediction != acc_actual, axis=1)))
                if crosscheck is not None:
                    for j in range(n_acc):
                        crosscheck["rows"] += 1
                        if not np.array_equal(matcher.decode(acc_syndrome[j]),
                                              prediction[j]):
                            crosscheck["mismatches"] += 1
            accepted += n_acc
            errors += n_err
            discarded += int(np.count_nonzero(fired))
            record.update({"accepted": n_acc, "errors": n_err})
        chunk_records.append(record)

    acceptance = accepted / shots
    ler = (errors / accepted) if accepted else None
    result = {"status": "ok", "reason": None, "seed": int(seed),
              "shots": shots, "schedule": sched, "accepted": accepted,
              "errors": errors, "discarded": discarded,
              "acceptance": acceptance, "ler": ler,
              "chunk_records": chunk_records, "crosscheck": crosscheck}
    if noiseless:
        result["detector_bits_fired"] = noiseless_detector_bits_fired
        result["actual_observable_flip_rows"] = noiseless_actual_flip_rows
    return result


# ------------------------------------------------------------------- analysis

def zscore(repro: float, author: float, sigma_combined: float) -> float | None:
    if not math.isfinite(sigma_combined) or sigma_combined <= 0:
        return None
    return (repro - author) / sigma_combined


def verdict_from_zs(zs: list[float | None]) -> str:
    """Preregistered mapping with decisive-mismatch precedence.

    Over the DEFINED statistics: any |z| >= 5 or >= 3 values beyond 3.53
    decides FROZEN-NEGATIVE first (undefined statistics must not mask a
    decisive mismatch); otherwise undefined statistics or 1-2 mild failures
    give FROZEN-INCONCLUSIVE; CERTIFIED only when every defined |z| <= 3.53.
    """
    present = [z for z in zs if z is not None]
    failures = [z for z in present if abs(z) > Z_THRESHOLD]
    if any(abs(z) >= DECISIVE_Z for z in present) or len(failures) >= 3:
        return "FROZEN-NEGATIVE"
    if len(present) != len(zs) or failures:
        return "FROZEN-INCONCLUSIVE"
    return "FROZEN-CERTIFIED"


def fit_through_origin(p_values, lers) -> float:
    """c in LER = c*p^2, unweighted through origin: sum(p^2*LER)/sum(p^4)."""
    denom = math.fsum(p ** 4 for p in p_values)
    numer = math.fsum(p * p * ler for p, ler in zip(p_values, lers))
    if denom <= 0:
        raise Refusal("degenerate through-origin fit (all p == 0)")
    return numer / denom


def fit_sigma_propagated(p_values, lers, sigmas) -> float:
    """Linear propagation: Var(c) = sum((p^2 / sum(p^4))^2 * sigma^2)."""
    denom = math.fsum(p ** 4 for p in p_values)
    return math.sqrt(math.fsum(
        (p * p / denom) ** 2 * s * s for p, s in zip(p_values, sigmas)))


def _author_row(variant: str, p_label: str):
    for row in AUTHOR_ROWS[variant]:
        if row[0] == p_label:
            return row
    raise Refusal(f"no author row for {variant}@{p_label}")


def _validate_point_rows(rows: list[dict], *,
                         expected_identity: dict | None = None,
                         require_complete: bool = True) -> None:
    """Validate the canonical production row, including recovery and ledger.

    The same validator is used immediately for a resumed/new point and for the
    final exact-12 analysis.  A copied scalar summary is never resumable.
    """
    expected_cells = {(v, label) for v in VARIANTS for label in P_LABELS}
    row_fields = {
        "identity", "variant", "p", "p_label", "status", "oracle",
        "flattened_equal", "recovery", "ended_utc", "reason", "seed",
        "shots", "schedule", "accepted", "errors", "discarded",
        "acceptance", "ler", "chunk_records", "crosscheck",
    }
    identity_fields = {
        "gate", "adapter_version", "run_id", "prereg_sha256", "repo",
        "commit", "tree_sha", "tar_sha256", "tar_size",
        "source_members_sha256", "env", "stream_note", "mode", "variant",
        "p_label", "shots", "seed", "chunk", "schedule",
    }
    point_identity_fields = {
        "mode", "variant", "p_label", "shots", "seed", "chunk", "schedule",
    }
    recovery_fields = {
        "variant", "p", "p_label", "recovery_mode",
        "sampled_circuit_source", "dem_source", "shipped_oracle_used",
        "driver_role", "driver_blob_sha1", "builder_role",
        "builder_blob_sha1", "surface_role", "surface_blob_sha1",
        "prefix_exec_sha256", "mask", "mask_sha256", "num_detectors",
        "num_observables", "shipped_path", "shipped_blob_sha1",
        "dem_sha256", "dem_error_instructions",
        "shipped_rebuilt_flattened_equal", "mask_source",
    }
    seen = set()
    shared_base_identity = None
    for row in rows:
        if not isinstance(row, dict):
            raise Refusal("point row is not a JSON object")
        key = (row.get("variant"), row.get("p_label"))
        if key in seen:
            raise Refusal(f"duplicate point row: {key}")
        if key not in expected_cells:
            raise Refusal(f"unknown point row: {key}")
        seen.add(key)
        if set(row) != row_fields:
            raise Refusal(
                f"point row fields for {key} are not the exact production "
                f"schema: got {sorted(row)}, expected {sorted(row_fields)}")

        author = _author_row(*key)
        shots = author[1]
        seed = FULL_SEEDS[key[0]][key[1]]
        sched = chunk_schedule(shots, CHUNK_FULL)
        if shots % CHUNK_FULL or sched["final_chunk"] != CHUNK_FULL:
            raise Refusal(f"full schedule for {key} is not all 10,000-shot "
                          "chunks")

        ident = row["identity"]
        if not isinstance(ident, dict) or set(ident) != identity_fields:
            raise Refusal(f"identity fields mismatch for {key}")
        if expected_identity is not None and ident != expected_identity:
            raise Refusal(f"identity mismatch for {key}")
        if (ident["gate"] != GATE
                or ident["adapter_version"] != ADAPTER_VERSION
                or ident["repo"] != REPO or ident["commit"] != COMMIT
                or ident["tree_sha"] != TREE_SHA
                or ident["tar_sha256"] != TAR_SHA256
                or ident["tar_size"] != TAR_SIZE
                or ident["stream_note"] != STREAM_NOTE
                or ident["mode"] != "full" or ident["variant"] != key[0]
                or ident["p_label"] != key[1] or ident["shots"] != shots
                or ident["seed"] != seed or ident["chunk"] != CHUNK_FULL
                or ident["schedule"] != sched):
            raise Refusal(f"identity pins mismatch for {key}")
        if (not isinstance(ident["run_id"], str)
                or RUN_ID_RE.fullmatch(ident["run_id"]) is None):
            raise Refusal(f"identity run_id malformed for {key}")
        for name in ("prereg_sha256", "source_members_sha256"):
            value = ident[name]
            if not isinstance(value, str) or re.fullmatch(
                    r"[0-9a-f]{64}", value) is None:
                raise Refusal(f"identity {name} malformed for {key}")
        env = ident["env"]
        if (not isinstance(env, dict)
                or set(env) != {"python", "platform", "machine", "stim",
                                "pymatching", "numpy"}
                or any(not isinstance(value, str) or not value
                       for value in env.values())
                or any(env[pkg] != version
                       for pkg, version in REQUIRED_VERSIONS.items())):
            raise Refusal(f"identity environment pins mismatch for {key}")
        base_identity = {
            name: value for name, value in ident.items()
            if name not in point_identity_fields
        }
        if shared_base_identity is None:
            shared_base_identity = base_identity
        elif base_identity != shared_base_identity:
            raise Refusal("point rows mix different run identities")

        expected_oracle = ("rebuilt" if key in REBUILT_ORACLE_CELLS
                           else "shipped")
        if row["oracle"] != expected_oracle:
            raise Refusal(f"row/identity disagreement for {key}: oracle "
                          f"{row['oracle']!r} != {expected_oracle!r}")
        if (row["shots"] != shots or row["seed"] != seed
                or row["schedule"] != sched or row["status"] != "ok"
                or row["p"] != float(key[1])
                or row["reason"] is not None
                or row["crosscheck"] is not None):
            raise Refusal(f"row/identity disagreement for {key}")
        if expected_oracle == "shipped":
            if row["flattened_equal"] is not True:
                raise Refusal(f"row/identity disagreement for {key}: "
                              "shipped-oracle row must record "
                              "flattened_equal True")
        elif not isinstance(row["flattened_equal"], bool):
            raise Refusal(f"row/identity disagreement for {key}: "
                          "rebuilt-oracle row must record the factual "
                          "flattened_equal bool")
        try:
            ended = datetime.strptime(row["ended_utc"], TS_Z)
        except (TypeError, ValueError) as exc:
            raise Refusal(f"ended_utc malformed for {key}") from exc
        if ended.strftime(TS_Z) != row["ended_utc"]:
            raise Refusal(f"ended_utc is not canonical for {key}")

        recovery = row["recovery"]
        if not isinstance(recovery, dict) or set(recovery) != recovery_fields:
            raise Refusal(f"recovery fields mismatch for {key}")
        roles = {
            "driver": f"driver_{key[0]}",
            "builder": f"builder_{key[0]}",
            "surface": f"surface_{key[0]}",
            "shipped": f"shipped_{key[0]}_{key[1]}",
        }
        if expected_oracle == "shipped":
            expected_mode = "shipped_equality_required"
            expected_source = "shipped_author_circuit"
            expected_oracle_used = True
            expected_fe: bool | None = True
        else:
            expected_mode = "rebuilt_oracle_r8"
            expected_source = "rebuilt_from_pinned_driver_builder"
            expected_oracle_used = False
            expected_fe = row["flattened_equal"]
        if (recovery["variant"] != key[0]
                or recovery["p"] != float(key[1])
                or recovery["p_label"] != key[1]
                or recovery["recovery_mode"] != expected_mode
                or recovery["sampled_circuit_source"] != expected_source
                or recovery["dem_source"] != expected_source
                or recovery["shipped_oracle_used"] is not expected_oracle_used
                or recovery["driver_role"] != roles["driver"]
                or recovery["driver_blob_sha1"] != _blob_of(roles["driver"])
                or recovery["builder_role"] != roles["builder"]
                or recovery["builder_blob_sha1"] != _blob_of(roles["builder"])
                or recovery["surface_role"] != roles["surface"]
                or recovery["surface_blob_sha1"] != _blob_of(roles["surface"])
                or recovery["shipped_blob_sha1"] != _blob_of(
                    roles["shipped"])
                or recovery["num_observables"] != 3
                or recovery["shipped_rebuilt_flattened_equal"] != expected_fe
                or recovery["mask_source"] !=
                "executed pinned driver prefix + pinned stim_builder"
                ".postselct_numbers(); never inferred from .stim"):
            raise Refusal(f"recovery pins mismatch for {key}")
        expected_shipped = (
            TARGET_DIR / "campaigns" / ident["run_id"] / "source"
            / REPO_ROOT_NAME / "stim" / key[0] / shipped_filename(key[1])
        ).resolve()
        if recovery["shipped_path"] != str(expected_shipped):
            raise Refusal(f"recovery shipped_path mismatch for {key}")
        mask = recovery["mask"]
        num_detectors = recovery["num_detectors"]
        if (not isinstance(mask, list)
                or any(type(index) is not int for index in mask)
                or mask != sorted(set(mask))
                or type(num_detectors) is not int or num_detectors <= 0
                or any(index < 0 or index >= num_detectors for index in mask)):
            raise Refusal(f"recovery detector mask malformed for {key}")
        expected_mask_sha = sha256_bytes(
            json.dumps(mask, separators=(",", ":")).encode("ascii"))
        if recovery["mask_sha256"] != expected_mask_sha:
            raise Refusal(f"recovery mask hash mismatch for {key}")
        for name in ("prefix_exec_sha256", "mask_sha256", "dem_sha256"):
            value = recovery[name]
            if not isinstance(value, str) or re.fullmatch(
                    r"[0-9a-f]{64}", value) is None:
                raise Refusal(f"recovery {name} malformed for {key}")
        if (type(recovery["dem_error_instructions"]) is not int
                or recovery["dem_error_instructions"] <= 0):
            raise Refusal(f"recovery DEM evidence malformed for {key}")

        accepted, errors = row["accepted"], row["errors"]
        if (type(accepted) is not int or type(errors) is not int
                or not (0 <= errors <= accepted <= shots)):
            raise Refusal(f"invariant violated for {key}: "
                          "0 <= errors <= accepted <= shots")
        if (row["discarded"] != shots - accepted
                or type(row["discarded"]) is not int):
            raise Refusal(f"invariant violated for {key}: discarded")
        if (type(row["acceptance"]) is not float
                or row["acceptance"] != accepted / shots):
            raise Refusal(f"invariant violated for {key}: acceptance")
        if accepted == 0:
            if row["ler"] is not None:
                raise Refusal(f"invariant violated for {key}: ler defined "
                              "at zero accepted shots")
        elif (type(row["ler"]) is not float
              or row["ler"] != errors / accepted):
            raise Refusal(f"invariant violated for {key}: ler")

        ledger = row["chunk_records"]
        if not isinstance(ledger, list) or len(ledger) != sched["num_chunks"]:
            raise Refusal(f"chunk ledger length mismatch for {key}")
        ledger_n = ledger_accepted = ledger_errors = 0
        for index, record in enumerate(ledger):
            if (not isinstance(record, dict)
                    or set(record) != {"n", "syndrome_obs_sha256",
                                       "accepted", "errors"}):
                raise Refusal(f"chunk ledger schema mismatch for {key} "
                              f"chunk {index}")
            expected_n = (sched["chunk"]
                          if index < sched["num_chunks"] - 1
                          else sched["final_chunk"])
            n = record["n"]
            chunk_accepted = record["accepted"]
            chunk_errors = record["errors"]
            if (type(n) is not int or n != expected_n
                    or type(chunk_accepted) is not int
                    or type(chunk_errors) is not int
                    or not (0 <= chunk_errors <= chunk_accepted <= n)):
                raise Refusal(f"chunk ledger count invariant violated for "
                              f"{key} chunk {index}")
            digest = record["syndrome_obs_sha256"]
            if not isinstance(digest, str) or re.fullmatch(
                    r"[0-9a-f]{64}", digest) is None:
                raise Refusal(f"chunk ledger digest malformed for {key} "
                              f"chunk {index}")
            ledger_n += n
            ledger_accepted += chunk_accepted
            ledger_errors += chunk_errors
        if (ledger_n != shots or ledger_accepted != accepted
                or ledger_errors != errors):
            raise Refusal(f"chunk ledger sums mismatch for {key}")

    if require_complete and (len(rows) != 12 or seen != expected_cells):
        missing = sorted(expected_cells - seen)
        raise Refusal(f"expected exactly 12 unique point rows covering all "
                      f"(variant, p) cells; got {len(rows)}, missing "
                      f"{missing}")


def analyze_points(point_rows: list[dict]) -> dict:
    """Full-mode final analysis: 24 z-scores, fits, preregistered verdict."""
    _validate_point_rows(point_rows)
    comparisons = []
    zs: list[float | None] = []
    for variant in VARIANTS:
        p_vals, lers, sigmas = [], [], []
        for row in sorted((r for r in point_rows
                           if r["variant"] == variant),
                          key=lambda r: -r["p"]):
            author = _author_row(variant, row["p_label"])
            author_ler, author_acc = author[2], author[3]
            author_shots = author[1]
            # Revision 8: a rebuilt-oracle row is comparable on the strength
            # of its own clean recovery+sampling; its factual
            # flattened_equal=False records F-Z4 and is not a defect.
            clean = (row["status"] == "ok"
                     and (row["oracle"] == "rebuilt"
                          or row["flattened_equal"]))
            if not clean:
                comparisons.append({"variant": variant, "p": row["p"],
                                    "comparable": False,
                                    "reason": row.get("reason")
                                    or "point not clean"})
                zs.append(None)
                zs.append(None)
                continue
            sa_ler = ler_sigma(author_shots, author_ler, author_acc)
            sa_acc = acc_sigma(author_shots, author_acc)
            sr_acc = acc_sigma(row["shots"], row["acceptance"])
            # acceptance z is always defined; only the LER z may be undefined
            z_acc = zscore(row["acceptance"], author_acc,
                           math.hypot(sa_acc, sr_acc))
            if row["accepted"] and row["ler"] is not None:
                sr_ler = ler_sigma(row["shots"], row["ler"],
                                   row["acceptance"])
                z_ler = zscore(row["ler"], author_ler,
                               math.hypot(sa_ler, sr_ler))
            else:
                sr_ler = None
                z_ler = None  # LER undefined at zero accepted shots
            zs.extend([z_ler, z_acc])
            entry = {
                "variant": variant, "p": row["p"], "comparable": True,
                "shots": row["shots"], "accepted": row["accepted"],
                "errors": row["errors"], "ler": row["ler"],
                "acceptance": row["acceptance"],
                "author": {"shots": author_shots, "ler": author_ler,
                           "acceptance": author_acc},
                "z_ler": z_ler, "z_acceptance": z_acc,
            }
            conflict = AUTHOR_ROW_SOURCE_CONFLICTS.get(
                (variant, row["p_label"]))
            if conflict:
                entry["author_row_source_conflict"] = conflict
            comparisons.append(entry)
            p_vals.append(row["p"])
            lers.append(row["ler"])
            sigmas.append(sr_ler)
    fits = {}
    for variant in VARIANTS:
        rows = [c for c in comparisons
                if c["variant"] == variant and c["comparable"]]
        # Reproduced fit over DEFINED-LER points only (a zero-accepted point
        # has no defined LER); labelled incomplete unless all six points
        # contributed.  Never crash before the verdict.
        fit_rows = [c for c in rows
                    if c["accepted"] and c["ler"] is not None]
        pv = [c["p"] for c in fit_rows]
        lv = [c["ler"] for c in fit_rows]
        sv = [ler_sigma(c["shots"], c["ler"], c["acceptance"])
              for c in fit_rows]
        fits[variant] = {
            "author_c": AUTHOR_FITS[variant],
            "reproduced_c": fit_through_origin(pv, lv) if pv else None,
            "reproduced_c_sigma": (fit_sigma_propagated(pv, lv, sv)
                                   if pv else None),
            "n_points_comparable": len(rows),
            "n_points_used": len(pv),
            "fit_complete": len(fit_rows) == 6,
            "fit_points_missing": 6 - len(pv),
        }
    verdict = verdict_from_zs(zs)
    return {
        "comparisons": comparisons,
        "n_comparisons": len(zs),
        "n_undefined": sum(1 for z in zs if z is None),
        "n_failures_beyond_threshold": sum(
            1 for z in zs if z is not None and abs(z) > Z_THRESHOLD),
        "n_decisive": sum(1 for z in zs
                          if z is not None and abs(z) >= DECISIVE_Z),
        "z_threshold": Z_THRESHOLD,
        "threshold_rule": ("familywise 1% Bonferroni over 24 comparisons, "
                           "rounded upward: exact z 3.529296088834735 -> 3.53 "
                           "(dependence-valid; exact Šidák 3.528022653482286 "
                           "rounds up to the same gate)"),
        "bonferroni_z_exact": bonferroni_z_exact(),
        "sidak_z_exact": sidak_z_exact(),
        "fits": fits,
        "fit_note": ("fit/slope is descriptive; authors print no fit "
                     "uncertainty; reproduced sigma is linear propagation "
                     "of accepted-conditioned binomial LER sigmas"),
        "candidate_verdict": verdict,
        "scope_note": ("verdict binds the executable d3/d9 artifact only; "
                       "the paper's grown d=7 claim remains NOT-REPRODUCED"),
    }



# ------------------------------------------------------------- run scaffolding

def source_members_digest(snapshot: dict) -> str:
    canon = json.dumps(snapshot["pins"], sort_keys=True,
                        separators=(",", ":")).encode("utf-8")
    return sha256_bytes(canon)


def prepare_run(run_dir: Path, tar_path: Path, prereg_path: Path) -> dict:
    """Capture control/source bytes once, validate, then write provenance."""
    run_dir = Path(run_dir).resolve()
    tar_path = Path(tar_path).resolve()
    prereg_path = Path(prereg_path).resolve()
    if not prereg_path.is_file():
        raise Refusal(f"prereg absent: {prereg_path}")
    prereg_bytes = prereg_path.read_bytes()
    manifest_path = Path(run_dir) / "manifest.json"
    try:
        manifest_bytes = manifest_path.read_bytes()
    except OSError as exc:
        raise Refusal(f"campaign manifest unreadable: {exc}") from exc
    manifest = validate_run_dir(
        run_dir, prereg_bytes, manifest_bytes=manifest_bytes)
    env = validate_env()
    refuse_tar_inside_run_dir(tar_path, run_dir)
    try:
        tar_bytes = Path(tar_path).read_bytes()
    except OSError as exc:
        raise Refusal(f"source tar unreadable: {tar_path}: {exc}") from exc
    tar_size = len(tar_bytes)
    if tar_size != TAR_SIZE:
        raise Refusal(
            f"source tar size {tar_size} != pinned {TAR_SIZE}: {tar_path}")
    tar_sha = sha256_bytes(tar_bytes)
    if tar_sha != TAR_SHA256:
        raise Refusal(
            f"source tar sha256 {tar_sha} != pinned {TAR_SHA256}: "
            f"{tar_path}; refusing to substitute an unpinned artifact")
    pristine = {"path": str(tar_path), "size": tar_size, "sha256": tar_sha}

    source_root = run_dir / "source"
    extraction = extract_pinned_members(
        io.BytesIO(tar_bytes), source_root, MEMBER_PINS)
    # The archived object is the exact byte string validated and extracted.
    tar_copy = source_root / "source_tar.gz"
    if tar_copy.exists():
        if tar_copy.read_bytes() != tar_bytes:
            raise Refusal(f"existing tar copy differs: {tar_copy}")
    else:
        atomic_write_bytes(tar_copy, tar_bytes)

    snapshots = snapshot_code_and_prereg(run_dir, prereg_bytes)
    assert_prereg_unchanged(prereg_path, prereg_bytes)
    assert_run_live(run_dir, manifest_bytes)
    identity = {
        "gate": GATE,
        "adapter_version": ADAPTER_VERSION,
        "run_id": manifest["run_id"],
        "prereg_sha256": manifest["prereg_sha256"],
        "repo": REPO, "commit": COMMIT, "tree_sha": TREE_SHA,
        "tar_sha256": TAR_SHA256, "tar_size": TAR_SIZE,
        "source_members_sha256": source_members_digest(extraction),
        "env": env,
        "stream_note": STREAM_NOTE,
    }
    return {"manifest": manifest, "manifest_bytes": manifest_bytes,
            "env": env, "pristine": pristine, "extraction": extraction,
            "identity": identity, "snapshots": snapshots,
            "source_root": source_root, "prereg_path": prereg_path,
            "prereg_bytes": prereg_bytes, "source_tar_bytes": tar_bytes}


def snapshot_code_and_prereg(run_dir: Path, prereg_bytes: bytes) -> dict:
    """Snapshot the code files and the CAPTURED prereg bytes (never re-read
    the live prereg here: snapshots must match what was validated)."""
    if not TEST_PATH.is_file():
        raise Refusal(f"snapshot source absent: {TEST_PATH}")
    snap_dir = run_dir / "provenance" / "snapshots"
    items = []
    for label, data in (("zero_level_repo_repro.py",
                         Path(__file__).read_bytes()),
                        ("test_zero_level_repo_repro.py",
                         TEST_PATH.read_bytes()),
                        ("pre_statement.md", prereg_bytes)):
        dest = snap_dir / label
        if dest.exists():
            if sha256_file(dest) != sha256_bytes(data):
                raise Refusal(f"existing snapshot differs: {dest}")
        else:
            dest.parent.mkdir(parents=True, exist_ok=True)
            atomic_write_bytes(dest, data)
        items.append({"name": label, "sha256": sha256_bytes(data),
                      "bytes": len(data)})
    return {"items": items}


def write_summary(run_dir: Path, payload: dict) -> Path:
    verdict = (payload.get("candidate_verdict")
               or payload.get("smoke_verdict") or "N/A")
    lines = [f"# {GATE} — {payload.get('mode')} summary", "",
             f"- run: {payload.get('run_id')}",
             f"- ended_utc: {payload.get('ended_utc')}",
             f"- candidate verdict: **{verdict}**",
             "",
             "Machine-readable artifacts: `results/` JSON files; this file is",
             "descriptive. Freezing/closing remains the campaign harness's act.",
             ""]
    out = run_dir / "results" / "summary.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_bytes(out, "\n".join(lines).encode("utf-8"))
    return out


def write_inventory(run_dir: Path) -> Path:
    """Audit inventory of every file the run dir now contains (except itself)."""
    entries = []
    root = run_dir.resolve()
    for path in sorted(root.rglob("*")):
        if path.is_file():
            rel = path.relative_to(root).as_posix()
            if rel == "results/inventory.json":
                continue
            entries.append({"path": rel, "sha256": sha256_file(path),
                            "bytes": path.stat().st_size})
    out = run_dir / "results" / "inventory.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_json(out, {"run_dir": str(root), "files": entries})
    return out


# ------------------------------------------------------------------ mode: full

def run_full(run_dir: Path, ctx: dict) -> dict:
    base_identity = ctx["identity"]
    manifest_bytes = ctx["manifest_bytes"]
    # Full mode executes the smoke battery FIRST: a certified result must
    # never bypass the author-sampled circuit-equality, DEM/API,
    # scalar-vs-batch, and deterministic-replay gates.  Any failure is a
    # REJECTED refusal before a single author-count point is sampled.
    assert_run_live(run_dir, manifest_bytes)
    smoke_gates, smoke_evidence = run_smoke_gates(ctx["source_root"])
    failed_gates = [g["gate"] for g in smoke_gates if not g["pass"]]
    if failed_gates:
        raise Refusal(
            f"smoke gates failed inside full mode: {failed_gates}; refusing "
            "to sample author-count points")
    assert_run_live(run_dir, manifest_bytes)
    points_dir = run_dir / "results" / "points"
    points_dir.mkdir(parents=True, exist_ok=True)

    point_specs = []
    for variant in VARIANTS:
        for p_label in P_LABELS:
            shots = _author_row(variant, p_label)[1]
            seed = FULL_SEEDS[variant][p_label]
            point_identity = dict(base_identity, mode="full", variant=variant,
                                  p_label=p_label, shots=shots, seed=seed,
                                  chunk=CHUNK_FULL,
                                  schedule=chunk_schedule(shots, CHUNK_FULL))
            point_specs.append({
                "variant": variant, "p_label": p_label, "shots": shots,
                "seed": seed, "identity": point_identity,
                "path": points_dir / f"{variant}_{p_label}.json",
            })

    # Validate every existing point before any new author-count sampling.  A
    # corrupt later row must not allow earlier missing rows to spend shots.
    resumed = {}
    for spec in point_specs:
        path = spec["path"]
        if not path.exists():
            continue
        assert_run_live(run_dir, manifest_bytes)
        try:
            row = json.loads(path.read_bytes().decode("utf-8"))
        except (OSError, UnicodeError, ValueError) as exc:
            raise Refusal(
                f"{path.name}: existing point JSON unreadable: {exc}") from exc
        _validate_point_rows(
            [row], expected_identity=spec["identity"],
            require_complete=False)
        assert_run_live(run_dir, manifest_bytes)
        resumed[(spec["variant"], spec["p_label"])] = row

    rows = []
    for spec in point_specs:
        assert_run_live(run_dir, manifest_bytes)
        variant = spec["variant"]
        p_label = spec["p_label"]
        key = (variant, p_label)
        if key in resumed:
            rows.append(resumed[key])
            continue

        rebuilt_oracle = key in REBUILT_ORACLE_CELLS
        recovered = recover_point(ctx["source_root"], variant, p_label,
                                  rebuilt_oracle=rebuilt_oracle)
        if not rebuilt_oracle and not recovered.flattened_equal:
            raise Refusal(
                f"{variant}@{p_label}: shipped != rebuilt (flattened()); "
                "source semantic failure — refusing before sampling")
        counts = sample_counts(
            recovered.circuit, recovered.mask, recovered.matcher,
            spec["shots"], CHUNK_FULL, spec["seed"])
        # Canonical row schema: status/flattened_equal live at the top;
        # recovery provenance and the complete sampling ledger are mandatory
        # for both fresh writes and resumed acceptance.  ``oracle`` names the
        # Revision-8 semantics: "shipped" (shipped artifact sampled, equality
        # mandatory) or "rebuilt" (pinned-driver rebuild sampled; shipped
        # artifact provenance-only, factual inequality recorded).
        row = {"identity": spec["identity"],
               "variant": variant, "p": float(p_label),
               "p_label": p_label, "status": "ok",
               "oracle": "rebuilt" if rebuilt_oracle else "shipped",
               "flattened_equal": bool(recovered.flattened_equal),
               "recovery": recovered.provenance, "ended_utc": _now()}
        row.update(counts)
        _validate_point_rows(
            [row], expected_identity=spec["identity"],
            require_complete=False)
        assert_run_live(run_dir, manifest_bytes)
        atomic_write_json(spec["path"], row)  # resumable per-point artifact
        rows.append(row)

    assert_prereg_unchanged(Path(ctx["prereg_path"]), ctx["prereg_bytes"])
    analysis = analyze_points(rows)
    analysis["identity"] = base_identity
    analysis["mode"] = "full"
    analysis["smoke_gates"] = smoke_gates
    analysis["smoke_evidence"] = dict(smoke_evidence, gates=smoke_gates)
    analysis["run_id"] = base_identity["run_id"]
    analysis["ticket_budget_shots"] = TICKET_BUDGET_SHOTS
    analysis["author_total_shots"] = AUTHOR_TOTAL_SHOTS
    analysis["ended_utc"] = _now()
    assert_run_live(run_dir, manifest_bytes)
    atomic_write_json(run_dir / "results" / "analysis.json", analysis)
    assert_run_live(run_dir, manifest_bytes)
    write_summary(run_dir, analysis)
    assert_run_live(run_dir, manifest_bytes)
    write_inventory(run_dir)
    return analysis


# ----------------------------------------------------------------- mode: smoke

def run_smoke_gates(source_root: Path) -> tuple[list[dict], dict]:
    """The advertised smoke battery, shared by smoke and full modes.

    The ungrown p=0 arm is builder-only: it samples the rebuilt pinned
    driver/builder circuit and makes no claim about the historical
    832_text_0.stim.  The grown p=1e-3 arm remains a shipped-vs-rebuilt
    flattened-equality oracle.  Both arms exercise DEM/API semantics,
    deterministic seeded sampling on the fixed chunk schedule, and same-seed
    replay; the grown arm also checks scalar-vs-batch bitwise identity on
    accepted rows.  Any failure raises Refusal (REJECTED) before scientific
    artifacts.
    """
    gates: list[dict] = []

    # Arm A — ungrown, noiseless (p=0), 10k shots, rebuilt circuit only.
    rec0 = recover_point(source_root, "ungrown", "0", builder_only=True)
    p0_provenance_ok = (
        rec0.circuit is rec0.rebuilt
        and rec0.flattened_equal is None
        and rec0.provenance.get("recovery_mode") == "builder_only_p0"
        and rec0.provenance.get("sampled_circuit_source")
        == "rebuilt_from_pinned_driver_builder"
        and rec0.provenance.get("dem_source")
        == "rebuilt_from_pinned_driver_builder"
        and rec0.provenance.get("shipped_oracle_used") is False
        and rec0.provenance.get("shipped_rebuilt_flattened_equal") is None
        and "shipped_path" not in rec0.provenance
        and "shipped_blob_sha1" not in rec0.provenance
    )
    gates.append(_gate(
        rec0.provenance, "ungrown_noiseless_builder_only_rebuilt",
        p0_provenance_ok))
    gates.append(_gate(
        rec0.provenance, "ungrown_noiseless_rebuilt_dem_zero_errors",
        rec0.matcher is None
        and rec0.provenance.get("dem_error_instructions") == 0))
    counts0 = sample_counts(rec0.circuit, rec0.mask, rec0.matcher,
                            SMOKE_NOISELESS_SHOTS, CHUNK_SMOKE,
                            SMOKE_SEEDS["ungrown_noiseless"])
    repeat0 = sample_counts(rec0.circuit, rec0.mask, rec0.matcher,
                            SMOKE_NOISELESS_SHOTS, CHUNK_SMOKE,
                            SMOKE_SEEDS["ungrown_noiseless"])
    gates.append(_gate(
        counts0, "ungrown_noiseless_rebuilt_zero_detector_observable_bits",
        counts0["accepted"] == SMOKE_NOISELESS_SHOTS
        and counts0["errors"] == 0
        and counts0.get("detector_bits_fired") == 0
        and counts0.get("actual_observable_flip_rows") == 0))
    gates.append(_gate(
        {"a": counts0["chunk_records"], "b": repeat0["chunk_records"]},
        "ungrown_noiseless_rebuilt_same_seed_repeat_bitwise",
        counts0["chunk_records"] == repeat0["chunk_records"]))

    # Arm B — grown, p=1e-3, 20k noisy shots + scalar/batch bitwise check.
    rec1 = recover_point(source_root, "grown", "0.001")
    gates.append(_gate(
        rec1.provenance, "grown_0.001_shipped_rebuilt_flattened_equal",
        rec1.flattened_equal))
    gates.append(_gate(rec1.provenance, "grown_0.001_dem_fault_ids",
                       rec1.matcher is not None))
    counts1 = sample_counts(rec1.circuit, rec1.mask, rec1.matcher,
                            SMOKE_GROWN_SHOTS, CHUNK_SMOKE,
                            SMOKE_SEEDS["grown_0.001"],
                            crosscheck_accepted_rows=True)
    repeat1 = sample_counts(rec1.circuit, rec1.mask, rec1.matcher,
                            SMOKE_GROWN_SHOTS, CHUNK_SMOKE,
                            SMOKE_SEEDS["grown_0.001"])
    xcheck = counts1.get("crosscheck") or {}
    gates.append(_gate(counts1, "grown_0.001_scalar_batch_bitwise",
                       xcheck.get("rows", 0) > 0
                       and xcheck.get("mismatches", 1) == 0))
    gates.append(_gate(
        {"a": counts1["chunk_records"], "b": repeat1["chunk_records"]},
        "grown_0.001_same_seed_repeat_bitwise",
        counts1.get("chunk_records") == repeat1.get("chunk_records")))

    failed = [g["gate"] for g in gates if not g["pass"]]
    if failed:
        raise Refusal(f"smoke semantic gates failed: {failed}; REJECTED "
                      "before scientific artifacts")
    return gates, {
        "builder_only_p0_recovery": rec0.provenance,
        "grown_author_sampled_recovery": rec1.provenance,
        "counts0": counts0, "repeat0": repeat0,
        "counts1": counts1, "repeat1": repeat1,
    }


def run_smoke(run_dir: Path, ctx: dict) -> dict:
    manifest_bytes = ctx["manifest_bytes"]
    assert_run_live(run_dir, manifest_bytes)
    gates, evidence = run_smoke_gates(ctx["source_root"])
    assert_prereg_unchanged(Path(ctx["prereg_path"]), ctx["prereg_bytes"])
    report = {"mode": "smoke", "run_id": ctx["identity"]["run_id"],
              "identity": ctx["identity"],
              "smoke_seed_note": ("smoke repeats re-use the same seed; the "
                                  "repeat gate verifies fixed-schedule "
                                  "bitwise count identity; any failing "
                                  "gate is a refusal, never a FAIL report"),
              "gates": gates,
              "smoke_verdict": "PASS",
              "counts": evidence,
              "ended_utc": _now()}
    assert_run_live(run_dir, manifest_bytes)
    (run_dir / "results").mkdir(parents=True, exist_ok=True)
    atomic_write_json(run_dir / "results" / "smoke.json", report)
    assert_run_live(run_dir, manifest_bytes)
    write_summary(run_dir, report)
    assert_run_live(run_dir, manifest_bytes)
    write_inventory(run_dir)
    return report


def _gate(provenance: dict, name: str, passed) -> dict:
    return {"gate": name, "pass": bool(passed),
            "detail": {k: v for k, v in provenance.items()
                       if k in ("variant", "p", "p_label", "recovery_mode",
                                "sampled_circuit_source", "dem_source",
                                "shipped_oracle_used", "mask_sha256",
                                "num_detectors", "num_observables",
                                "shipped_rebuilt_flattened_equal",
                                "dem_sha256", "dem_error_instructions",
                                "status", "reason", "detector_bits_fired",
                                "actual_observable_flip_rows", "a", "b")}}


# ------------------------------------------------------------------------ CLI

def _now() -> str:
    return datetime.now(timezone.utc).strftime(TS_Z)


def build_parser():
    import argparse
    parser = argparse.ArgumentParser(
        prog="zero_level_repo_repro.py",
        description="Gate %s: author-repository reproduction of the "
                    "zero-level CCZ simulations (arXiv:2605.21867). "
                    "Run inside a live campaign.py run dir; never creates, "
                    "freezes, or closes it." % GATE)
    parser.add_argument("--source-tar", required=True, type=Path,
                        help="codeload tar.gz of the pinned commit "
                             "(sha256 %s…)" % TAR_SHA256[:12])
    parser.add_argument("--run-dir", required=True, type=Path,
                        help="existing campaign.py run dir (status RUNNING)")
    parser.add_argument("--mode", required=True, choices=("smoke", "full"),
                        help="smoke = gates only; full = 12 author-count "
                             "points + analysis")
    return parser


def main(argv: list[str] | None = None) -> int:
    ns = build_parser().parse_args(argv)
    run_dir = ns.run_dir.expanduser().resolve()
    tar_path = ns.source_tar.expanduser().resolve()
    try:
        ctx = prepare_run(run_dir, tar_path, PREREG_PATH)
        report = run_smoke(run_dir, ctx) if ns.mode == "smoke" \
            else run_full(run_dir, ctx)
    except Refusal as exc:
        print(f"REJECTED: {exc}", file=sys.stderr)
        return 1
    verdict = report.get("candidate_verdict") or report.get("smoke_verdict")
    print(f"mode      : {ns.mode}")
    print(f"run dir   : {run_dir}")
    print(f"outcome   : {verdict}")
    print(f"artifacts : {run_dir / 'results'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
