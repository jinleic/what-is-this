# -*- coding: utf-8 -*-
"""Deterministic, sampling-free tests for zero_level_repo_repro.py.

These tests never build the upstream circuits, never sample, never touch the
network, and never write outside pytest/tmp directories.  They pin:
  * git blob sha1 hashing and the full member pin table (incl. the ticket's
    truncated blob prefixes and the two full driver SHAs);
  * the embedded author-row arithmetic (integer error-count reconstruction,
    the 71m-shot executable schedule, and the named F-Z3 comment/row
    conflict on ungrown p=1e-4) against the authors' own std formulas;
  * through-origin c fits (denominator sum(p^4)) and the propagated fit
    uncertainty (sum((p^2/sum(p^4))^2 sigma^2), sqrt at the end);
  * the exact Bonferroni / Šidák two-sided criticals and the 3.53 gate;
  * the verdict mapping INCLUDING precedence of a decisive mismatch over
    undefined statistics, and acceptance scoring at zero accepted shots;
  * vectorized postselection/error counting with fake arrays and a fake
    matcher, including accepted-only decode and noiseless actual observables;
  * local pinned imports suppressing bytecode and restoring interpreter state;
  * builder-only p=0 recovery never opens a shipped oracle and records no
    equality claim, while author-sampled recovery refuses inequality itself;
  * tar traversal/pin/no-overwrite handling and immutable prepare-time capture;
  * campaign.py parent, basename, mint binder, timestamp, and live markers;
  * strict full-row recovery/chunk ledgers, scalar-only refusal before new
    point sampling, and complete-row resume after the smoke-first battery;
  * CLI return codes (0 ok / 1 REJECTED refusal / 2 argparse).
"""

import io
import json
import math
import sys
import tarfile
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))
import zero_level_repo_repro as Z  # noqa: E402

# ----------------------------------------------------------------- constants

EMPTY_BLOB = "e69de29bb2d1d6434b8b29ae775ad8c2e48c5391"
HELLO_NL_BLOB = "ce013625030ba8dba906f756967f9e9ca394464a"
HELLO_BLOB = "b6fc4c620b67d95f953a5c1c1230aaab5db5a1b0"

TICKET_BLOB_PREFIXES = {
    "driver_ungrown": "6205e13302cbe0285f3e704b9a66da1df3ee6173",
    "driver_grown": "592d1a8422ae2aad47b37d66e694c4eeb2b0f019",
    "builder_ungrown": "f6c25a55",
    "builder_grown": "7d06f315",
    "surface_ungrown": "1dc06f18",
    "surface_grown": "3af84832",
    "std_calc_ungrown": "d2df580b",
    "std_calc_grown": "0e0aaf3e",
    "plot": "35d58268",
    "license": "b9d9364d",
}

# Published std lists from plot_LER_suc/CCZ_plot.py (author std_calc formula,
# evaluated by the authors at the NOMINAL comment budgets).
PLOT_LER_STD = {
    "ungrown": [1.1362894640086647e-05, 8.854978406193676e-06,
                6.026123395572438e-06, 3.6310293778435086e-06,
                1.6668438679594455e-06, 8.308705429272706e-07],
    "grown": [3.556505232335755e-05, 1.1151573302085898e-05,
              7.1987224382979435e-06, 4.2602026578644015e-06,
              1.3258539944309626e-06, 6.6317993647621e-07],
}
PLOT_ACC_STD = {
    "ungrown": [0.00022063016958521336, 0.00022360477372777174,
                0.00021972681520098544, 0.00020387567238687405,
                0.00016399027234952686, 0.0001236481082426739],
    "grown": [0.0004467387391798477, 0.00021424745157016923,
              0.00022294408298245547, 0.000219313347469706,
              0.0001324960805148741, 0.00010304238244132363],
}
# Corrected ungrown p=1e-4 author sigmas under the authoritative 1e7 reading.
UNGROWN_1E4_SIGMA_LER_10M = 5.875141951920214e-07
UNGROWN_1E4_SIGMA_ACC_10M = 8.743241581928296e-05

BONFERRONI_Z_EXACT = 3.529296088834735
SIDAK_Z_EXACT = 3.528022653482286


# ------------------------------------------------------------ git blob pins

def test_git_blob_sha1_known_vectors():
    assert Z.git_blob_sha1(b"") == EMPTY_BLOB
    assert Z.git_blob_sha1(b"hello\n") == HELLO_NL_BLOB
    assert Z.git_blob_sha1(b"hello") == HELLO_BLOB


def test_member_pin_table_matches_ticket_and_is_wellformed():
    assert set(Z.MEMBER_PINS) == {
        "driver_ungrown", "driver_grown", "builder_ungrown", "builder_grown",
        "surface_ungrown", "surface_grown", "std_calc_ungrown",
        "std_calc_grown", "plot", "license",
        "shipped_ungrown_0", "shipped_ungrown_0.0001", "shipped_ungrown_0.0002",
        "shipped_ungrown_0.0004", "shipped_ungrown_0.0006",
        "shipped_ungrown_0.0008", "shipped_ungrown_0.001",
        "shipped_grown_0.0001", "shipped_grown_0.0002", "shipped_grown_0.0004",
        "shipped_grown_0.0006", "shipped_grown_0.0008", "shipped_grown_0.001",
    }
    for role, (rel, blob) in Z.MEMBER_PINS.items():
        assert len(blob) == 40 and all(c in "0123456789abcdef" for c in blob)
        if role in TICKET_BLOB_PREFIXES:
            assert blob.startswith(TICKET_BLOB_PREFIXES[role]), role
    # the two drivers are pinned at full length by the ticket
    assert Z.MEMBER_PINS["driver_ungrown"][1] == TICKET_BLOB_PREFIXES["driver_ungrown"]
    assert Z.MEMBER_PINS["driver_grown"][1] == TICKET_BLOB_PREFIXES["driver_grown"]
    assert Z.COMMIT == "1b59e223590492e224bd8623a4e0bcba59029e01"
    assert Z.TREE_SHA == "9a5d89401fd40554635fb0e9d0b4da818670c2bc"
    assert Z.TAR_SHA256.startswith("d67dbe7482b391984da5e64aeff7668")
    assert Z.TAR_SIZE == 1_072_583
    assert Z.GATE == "zero-level-author-repro-r8"
    assert Z.ADAPTER_VERSION == "r8"
    assert Z.PREREG_REVISION_MARKER == "## Revision 8"
    assert Z.REBUILT_ORACLE_CELLS == frozenset({
        ("grown", "0.0008"), ("grown", "0.0006"),
        ("grown", "0.0004"), ("grown", "0.0001"),
    })


# ------------------------------------------------- author rows and arithmetic

def test_author_row_schedule_is_71m_and_conflict_named():
    shots = {v: [r[1] for r in Z.AUTHOR_ROWS[v]] for v in Z.VARIANTS}
    assert shots["ungrown"] == [5_000_000, 5_000_000, 5_000_000,
                                5_000_000, 5_000_000, 10_000_000]
    assert shots["grown"] == [1_000_000, 5_000_000, 5_000_000,
                              5_000_000, 10_000_000, 10_000_000]
    assert Z.AUTHOR_TOTAL_SHOTS == 71_000_000 == Z.TICKET_BUDGET_SHOTS
    conflict = Z.AUTHOR_ROW_SOURCE_CONFLICTS[("ungrown", "0.0001")]
    assert conflict["comment_says_shots"] == 5_000_000
    assert conflict["arithmetic_implies_shots"] == 10_000_000
    assert conflict["implied_accepted"] == 9_166_002
    assert conflict["implied_errors"] == 29
    assert conflict["nominal_reading_errors"] == 14.5  # impossible fraction


def test_author_rows_reconcile_to_integer_error_counts():
    for variant in Z.VARIANTS:
        for label, shots, ler, acc in Z.AUTHOR_ROWS[variant]:
            accepted = shots * acc
            errors = ler * accepted
            if (variant, label) == ("ungrown", "0.0001"):
                assert accepted == 9_166_002
                assert abs(errors - 29.0) < 1e-6  # authoritative 1e7 reading
            else:
                assert abs(errors - round(errors)) < 1e-6, (variant, label)


def test_author_std_formulas_reproduce_published_lists_and_10m_correction():
    # The authors' plotted stds come from their std_calc.py formulas at each
    # row's own budget (grown [1e6,5e6,5e6,5e6,1e7,1e7]); for ungrown p=1e-4
    # the published std used the stale NOMINAL 5e6 comment budget while the
    # published LER is 1e7-based — F-Z3; the executable row is 1e7-shot below.
    for variant in Z.VARIANTS:
        for i, (label, shots, ler, acc) in enumerate(Z.AUTHOR_ROWS[variant]):
            published_shots = 5_000_000 if (variant, label) == (
                "ungrown", "0.0001") else shots
            assert Z.ler_sigma(published_shots, ler, acc) == \
                PLOT_LER_STD[variant][i]
            assert Z.acc_sigma(published_shots, acc) == \
                PLOT_ACC_STD[variant][i]
    row = Z._author_row("ungrown", "0.0001")
    assert row[1] == 10_000_000
    assert Z.ler_sigma(row[1], row[2], row[3]) == UNGROWN_1E4_SIGMA_LER_10M
    assert Z.acc_sigma(row[1], row[3]) == UNGROWN_1E4_SIGMA_ACC_10M


# ------------------------------------------------------------------ fits

def test_through_origin_fit_uses_sum_p4_and_matches_author_constants():
    p = [float(r[0]) for r in Z.AUTHOR_ROWS["ungrown"]]
    ler_u = [r[2] for r in Z.AUTHOR_ROWS["ungrown"]]
    ler_g = [r[2] for r in Z.AUTHOR_ROWS["grown"]]
    assert abs(Z.fit_through_origin(p, ler_u) / 282.1598546872 - 1) < 1e-12
    assert abs(Z.fit_through_origin(p, ler_g) / 346.5564346243 - 1) < 1e-12
    # model identification: perfect quadratic data recover the constant
    assert abs(Z.fit_through_origin([1e-3, 1e-4], [3e-6, 3e-8]) - 3.0) < 1e-9
    # denominator is sum(p^4), NOT sum(p^2): the synthetic case below is
    # exact for c=1 and would explode under a sum(p^2) denominator
    assert abs(Z.fit_through_origin([2.0, 1.0], [4.0, 1.0]) - 1.0) < 1e-12


def test_fit_sigma_propagated_is_sqrt_of_sum_weight_sq_sigma_sq():
    p = [1e-3, 1e-4]
    sig = [2.0, 5.0]
    denom = p[0] ** 4 + p[1] ** 4
    w = [p[0] ** 2 / denom, p[1] ** 2 / denom]
    expected = math.sqrt((w[0] * 2.0) ** 2 + (w[1] * 5.0) ** 2)
    assert abs(Z.fit_sigma_propagated(p, [1e-6, 1e-8], sig) - expected) < 1e-15
    # sigma^2 (not sigma) enters: doubling every sigma doubles the result
    doubled = Z.fit_sigma_propagated(p, [1e-6, 1e-8], [4.0, 10.0])
    assert abs(doubled - 2.0 * expected) < 1e-15


def _base_identity(run_id=None, prereg_sha256="a" * 64):
    if run_id is None:
        stamp = "20260903T000000Z"
        uuid8 = "deadbeef"
        binder = "\x1f".join(
            (Z.GATE, "tester", prereg_sha256, stamp, uuid8))
        run_id = f"{stamp}_{uuid8}_{Z.sha256_bytes(binder.encode())[:12]}"
    return {
        "gate": Z.GATE, "adapter_version": Z.ADAPTER_VERSION,
        "run_id": run_id, "prereg_sha256": prereg_sha256,
        "repo": Z.REPO, "commit": Z.COMMIT, "tree_sha": Z.TREE_SHA,
        "tar_sha256": Z.TAR_SHA256, "tar_size": Z.TAR_SIZE,
        "source_members_sha256": "b" * 64,
        "env": {"python": "3.test", "platform": "test", "machine": "test",
                **Z.REQUIRED_VERSIONS},
        "stream_note": Z.STREAM_NOTE,
    }


def _oracle_of(variant, p_label):
    return ("rebuilt" if (variant, p_label) in Z.REBUILT_ORACLE_CELLS
            else "shipped")


def _recovery(variant, p_label, identity, source_root=None):
    if source_root is None:
        source_root = (Z.TARGET_DIR / "campaigns" / identity["run_id"]
                       / "source")
    mask = [0, 2]
    rebuilt = _oracle_of(variant, p_label) == "rebuilt"
    return {
        "variant": variant, "p": float(p_label), "p_label": p_label,
        "recovery_mode": ("rebuilt_oracle_r8" if rebuilt
                          else "shipped_equality_required"),
        "sampled_circuit_source": ("rebuilt_from_pinned_driver_builder"
                                   if rebuilt else "shipped_author_circuit"),
        "dem_source": ("rebuilt_from_pinned_driver_builder"
                       if rebuilt else "shipped_author_circuit"),
        "shipped_oracle_used": not rebuilt,
        "driver_role": f"driver_{variant}",
        "driver_blob_sha1": Z.MEMBER_PINS[f"driver_{variant}"][1],
        "builder_role": f"builder_{variant}",
        "builder_blob_sha1": Z.MEMBER_PINS[f"builder_{variant}"][1],
        "surface_role": f"surface_{variant}",
        "surface_blob_sha1": Z.MEMBER_PINS[f"surface_{variant}"][1],
        "prefix_exec_sha256": Z.sha256_bytes(
            f"prefix:{variant}:{p_label}".encode()),
        "mask": mask,
        "mask_sha256": Z.sha256_bytes(
            json.dumps(mask, separators=(",", ":")).encode("ascii")),
        "num_detectors": 4, "num_observables": 3,
        "shipped_path": str(
            source_root / Z.REPO_ROOT_NAME / "stim" / variant
            / Z.shipped_filename(p_label)),
        "shipped_blob_sha1":
            Z.MEMBER_PINS[f"shipped_{variant}_{p_label}"][1],
        "dem_sha256": Z.sha256_bytes(
            f"dem:{variant}:{p_label}".encode()),
        "dem_error_instructions": 1,
        "shipped_rebuilt_flattened_equal": not rebuilt,
        "mask_source": ("executed pinned driver prefix + pinned stim_builder"
                        ".postselct_numbers(); never inferred from .stim"),
    }


def _chunk_records(variant, p_label, shots, accepted, errors):
    sched = Z.chunk_schedule(shots, Z.CHUNK_FULL)
    remaining_accepted = accepted
    remaining_errors = errors
    records = []
    for index in range(sched["num_chunks"]):
        n = sched["chunk"] if index < sched["num_chunks"] - 1 \
            else sched["final_chunk"]
        chunk_accepted = min(n, remaining_accepted)
        chunk_errors = min(chunk_accepted, remaining_errors)
        remaining_accepted -= chunk_accepted
        remaining_errors -= chunk_errors
        records.append({
            "n": n,
            "syndrome_obs_sha256": Z.sha256_bytes(
                f"{variant}:{p_label}:{index}".encode()),
            "accepted": chunk_accepted, "errors": chunk_errors,
        })
    assert remaining_accepted == remaining_errors == 0
    return records


def _point_row(variant, i, *, identity_base=None, source_root=None,
               accepted=None, errors=None):
    label, shots, ler, acc = Z.AUTHOR_ROWS[variant][i]
    seed = Z.FULL_SEEDS[variant][label]
    if accepted is None:
        accepted = int(shots * acc)
    if errors is None:
        errors = int(round(ler * accepted))
    identity_base = dict(identity_base or _base_identity())
    identity = dict(
        identity_base, mode="full", variant=variant, p_label=label,
        shots=shots, seed=seed, chunk=Z.CHUNK_FULL,
        schedule=Z.chunk_schedule(shots, Z.CHUNK_FULL))
    return {
        "identity": identity,
        "variant": variant, "p": float(label), "p_label": label,
        "status": "ok", "oracle": _oracle_of(variant, label),
        "flattened_equal": _oracle_of(variant, label) == "shipped",
        "recovery": _recovery(variant, label, identity_base, source_root),
        "ended_utc": "2026-09-03T00:00:00Z",
        "reason": None, "seed": seed, "shots": shots,
        "schedule": Z.chunk_schedule(shots, Z.CHUNK_FULL),
        "accepted": accepted, "errors": errors,
        "discarded": shots - accepted,
        "acceptance": accepted / shots,
        "ler": errors / accepted if accepted else None,
        "chunk_records": _chunk_records(
            variant, label, shots, accepted, errors),
        "crosscheck": None,
    }


def test_full_analysis_schema_uses_top_level_fields_and_labels_fits():
    rows = [_point_row("ungrown", i) for i in range(6)]
    rows += [_point_row("grown", i) for i in range(6)]
    analysis = Z.analyze_points(rows)  # must not raise on the canonical rows
    assert analysis["n_comparisons"] == 24
    assert all(c["comparable"] for c in analysis["comparisons"])
    for variant in Z.VARIANTS:
        fit = analysis["fits"][variant]
        assert fit["fit_complete"] is True
        assert fit["n_points_used"] == 6
        assert math.isfinite(fit["reproduced_c"])
        assert math.isfinite(fit["reproduced_c_sigma"])


def test_zero_accepted_rows_score_acceptance_label_fit_incomplete():
    rows = [_point_row("ungrown", 0, accepted=0, errors=0)]
    rows += [_point_row("ungrown", i) for i in range(1, 6)]
    rows += [_point_row("grown", i) for i in range(6)]
    analysis = Z.analyze_points(rows)  # never crashes before the verdict
    comp = analysis["comparisons"][0]
    author = Z._author_row("ungrown", "0.001")
    expected_sigma = math.hypot(
        Z.acc_sigma(author[1], author[3]),
        Z.acc_sigma(rows[0]["shots"], 0.0))
    assert abs(comp["z_acceptance"] - (0.0 - author[3]) / expected_sigma) < 1e-9
    assert comp["z_ler"] is None
    fit = analysis["fits"]["ungrown"]
    assert fit["fit_complete"] is False
    assert fit["n_points_used"] == 5
    assert fit["fit_points_missing"] == 1
    assert math.isfinite(fit["reproduced_c"])
    # decisive acceptance mismatch wins over the undefined LER statistic
    assert analysis["candidate_verdict"] == "FROZEN-NEGATIVE"
    assert analysis["n_undefined"] >= 1


# ---------------------------------------------------- criticals and verdicts

def test_exact_bonferroni_and_sidak_criticals_round_to_3_53():
    assert abs(Z.bonferroni_z_exact() - BONFERRONI_Z_EXACT) < 1e-12
    assert abs(Z.sidak_z_exact() - SIDAK_Z_EXACT) < 1e-12
    assert math.ceil(Z.bonferroni_z_exact() * 100) / 100 == 3.53
    assert math.ceil(Z.sidak_z_exact() * 100) / 100 == 3.53
    # the preregistered derivation is the dependence-valid Bonferroni bound
    assert Z.Z_THRESHOLD == 3.53


def test_verdict_boundaries_and_decisive_precedence():
    assert Z.verdict_from_zs([0.0] * 24) == "FROZEN-CERTIFIED"
    assert Z.verdict_from_zs([3.53] * 24) == "FROZEN-CERTIFIED"  # |z| <= 3.53
    assert Z.verdict_from_zs([3.5300001] + [0.0] * 23) == "FROZEN-INCONCLUSIVE"
    assert Z.verdict_from_zs([4.999] + [0.0] * 23) == "FROZEN-INCONCLUSIVE"
    assert Z.verdict_from_zs([5.0] + [0.0] * 23) == "FROZEN-NEGATIVE"  # >= 5
    assert Z.verdict_from_zs([3.6, 3.6] + [0.0] * 22) == "FROZEN-INCONCLUSIVE"
    assert Z.verdict_from_zs([3.6, 3.6, 3.6] + [0.0] * 21) == "FROZEN-NEGATIVE"
    # precedence: a decisive value is NOT downgraded by undefined statistics
    assert Z.verdict_from_zs([None] * 24) == "FROZEN-INCONCLUSIVE"
    assert Z.verdict_from_zs([None] * 23 + [5.0]) == "FROZEN-NEGATIVE"
    assert Z.verdict_from_zs([None] * 22 + [3.6, 3.6]) == "FROZEN-INCONCLUSIVE"
    assert Z.verdict_from_zs([None] * 21 + [3.6, 3.6, 3.6]) == "FROZEN-NEGATIVE"


# -------------------------------------------- pinned import bytecode hygiene

def test_exec_driver_prefix_suppresses_bytecode_and_restores_flag(tmp_path):
    variant_dir = tmp_path / "variant"
    variant_dir.mkdir()
    (variant_dir / "surface_func.py").write_text(
        "VALUE = 7\n", encoding="utf-8")
    (variant_dir / "stim_builder.py").write_text(
        "import surface_func\n"
        "class Builder:\n"
        "    def build(self):\n"
        "        return [surface_func.VALUE]\n",
        encoding="utf-8")
    prefix = (
        "from stim_builder import Builder\n"
        "error_rate = [0.25]\n"
        "circuit_builder = Builder()\n"
        "circuit = circuit_builder.build()\n"
    )
    saved_modules = {
        name: sys.modules.pop(name, None)
        for name in ("stim_builder", "surface_func")
    }
    original_flag = sys.dont_write_bytecode
    sys.dont_write_bytecode = False
    try:
        _builder, built, _digest = Z._exec_driver_prefix(
            prefix, variant_dir, 0.001)
        assert built == [7]
        assert sys.dont_write_bytecode is False
        assert not any(path.name == "__pycache__"
                       for path in variant_dir.rglob("*"))
        assert not any(variant_dir.rglob("*.pyc"))
        sys.dont_write_bytecode = True
        with pytest.raises(Z.Refusal, match="expected objects"):
            Z._exec_driver_prefix(
                prefix.replace("circuit = circuit_builder.build()\n",
                               "circuit = None\n"),
                variant_dir, 0.001)
        assert sys.dont_write_bytecode is True
        assert not any(path.name == "__pycache__"
                       for path in variant_dir.rglob("*"))
        assert not any(variant_dir.rglob("*.pyc"))
    finally:
        sys.dont_write_bytecode = original_flag
        for name, module in saved_modules.items():
            sys.modules.pop(name, None)
            if module is not None:
                sys.modules[name] = module


class _RecoveryInstruction:
    def __init__(self, instruction_type):
        self.type = instruction_type


class _RecoveryDem:
    def __init__(self, error_instructions):
        self.instructions = [
            _RecoveryInstruction("error") for _ in range(error_instructions)
        ]

    def __iter__(self):
        return iter(self.instructions)

    def __str__(self):
        return "fake detector error model"


class _RecoveryCircuit:
    def __init__(self, flattened_value, *, error_instructions=0,
                 dem_allowed=True):
        self.flattened_value = flattened_value
        self.num_detectors = 4
        self.num_observables = 3
        self.error_instructions = error_instructions
        self.dem_allowed = dem_allowed
        self.dem_calls = 0

    def flattened(self):
        return self.flattened_value

    def detector_error_model(self, *, decompose_errors):
        assert decompose_errors is True
        self.dem_calls += 1
        assert self.dem_allowed, "inequality must refuse before DEM creation"
        return _RecoveryDem(self.error_instructions)


class _RecoveryBuilder:
    def postselct_numbers(self):
        return [True, False, True, False]


def _write_recovery_driver(tmp_path, variant):
    source_root = tmp_path / "source"
    vdir = source_root / Z.REPO_ROOT_NAME / "stim" / variant
    vdir.mkdir(parents=True)
    driver_name = Path(Z.MEMBER_PINS[f"driver_{variant}"][0]).name
    (vdir / driver_name).write_text(
        "error_rate = [1e-4]\n"
        "circuit = circuit_builder.build()\n",
        encoding="utf-8")
    return source_root, vdir


def test_builder_only_p0_uses_rebuilt_circuit_without_shipped_oracle(
        tmp_path, monkeypatch):
    source_root, vdir = _write_recovery_driver(tmp_path, "ungrown")
    rebuilt = _RecoveryCircuit("rebuilt")
    builder = _RecoveryBuilder()
    monkeypatch.setattr(
        Z, "_exec_driver_prefix",
        lambda prefix, variant_dir, p: (builder, [rebuilt], "a" * 64))

    class NoShippedStim:
        @staticmethod
        def Circuit(_text):
            pytest.fail("builder-only p=0 must not construct a shipped circuit")

    class FakeMatching:
        @staticmethod
        def from_detector_error_model(_dem):
            pytest.fail("noiseless builder-only p=0 must not build a matcher")

    class FakePyMatching:
        Matching = FakeMatching

    monkeypatch.setitem(sys.modules, "stim", NoShippedStim)
    monkeypatch.setitem(sys.modules, "pymatching", FakePyMatching)
    shipped_path = vdir / Z.shipped_filename("0")
    assert not shipped_path.exists()

    recovered = Z.recover_point(
        source_root, "ungrown", "0", builder_only=True)

    assert recovered.circuit is rebuilt
    assert recovered.rebuilt is rebuilt
    assert recovered.mask == [0, 2]
    assert recovered.matcher is None
    assert recovered.flattened_equal is None
    assert rebuilt.dem_calls == 1
    provenance = recovered.provenance
    assert provenance["recovery_mode"] == "builder_only_p0"
    assert provenance["sampled_circuit_source"] == \
        "rebuilt_from_pinned_driver_builder"
    assert provenance["dem_source"] == "rebuilt_from_pinned_driver_builder"
    assert provenance["mask"] == [0, 2]
    assert provenance["mask_source"] == (
        "executed pinned driver prefix + pinned stim_builder"
        ".postselct_numbers(); never inferred from .stim")
    assert len(provenance["mask_sha256"]) == 64
    assert len(provenance["dem_sha256"]) == 64
    assert provenance["shipped_oracle_used"] is False
    assert provenance["shipped_rebuilt_flattened_equal"] is None
    assert provenance["dem_error_instructions"] == 0
    assert "shipped_path" not in provenance
    assert "shipped_blob_sha1" not in provenance


def test_author_sampled_recovery_cannot_bypass_or_survive_inequality(
        tmp_path, monkeypatch):
    with pytest.raises(Z.Refusal, match="restricted to the ungrown p=0"):
        Z.recover_point(
            tmp_path, "grown", "0.001", builder_only=True)

    source_root, vdir = _write_recovery_driver(tmp_path, "grown")
    (vdir / Z.shipped_filename("0.001")).write_text(
        "shipped author circuit", encoding="utf-8")
    rebuilt = _RecoveryCircuit("rebuilt", dem_allowed=False)
    shipped = _RecoveryCircuit("shipped", dem_allowed=False)
    builder = _RecoveryBuilder()
    monkeypatch.setattr(
        Z, "_exec_driver_prefix",
        lambda prefix, variant_dir, p: (builder, [rebuilt], "b" * 64))
    opened = []

    class FakeStim:
        @staticmethod
        def Circuit(text):
            opened.append(text)
            return shipped

    class FakeMatching:
        @staticmethod
        def from_detector_error_model(_dem):
            pytest.fail("inequality must refuse before matcher construction")

    class FakePyMatching:
        Matching = FakeMatching

    monkeypatch.setitem(sys.modules, "stim", FakeStim)
    monkeypatch.setitem(sys.modules, "pymatching", FakePyMatching)
    with pytest.raises(Z.Refusal, match=r"shipped != rebuilt .*refusing before sampling"):
        Z.recover_point(source_root, "grown", "0.001")
    assert opened == ["shipped author circuit"]
    assert rebuilt.dem_calls == shipped.dem_calls == 0


def test_rebuilt_oracle_r8_samples_rebuild_and_records_inequality(
        tmp_path, monkeypatch):
    # rebuilt-oracle is refused outside the F-Z4 split cells, and the two
    # non-shipped recovery modes are mutually exclusive.
    with pytest.raises(Z.Refusal, match="restricted to the Revision-8"):
        Z.recover_point(tmp_path, "grown", "0.001", rebuilt_oracle=True)
    with pytest.raises(Z.Refusal, match="restricted to the Revision-8"):
        Z.recover_point(tmp_path, "ungrown", "0.0008", rebuilt_oracle=True)
    with pytest.raises(Z.Refusal, match="mutually exclusive"):
        Z.recover_point(tmp_path, "grown", "0.0008", builder_only=True,
                        rebuilt_oracle=True)

    source_root, vdir = _write_recovery_driver(tmp_path, "grown")
    (vdir / Z.shipped_filename("0.0008")).write_text(
        "shipped split-generation circuit", encoding="utf-8")
    rebuilt = _RecoveryCircuit("rebuilt", error_instructions=1)
    shipped = _RecoveryCircuit("shipped-DIFFERENT", dem_allowed=False)
    builder = _RecoveryBuilder()
    monkeypatch.setattr(
        Z, "_exec_driver_prefix",
        lambda prefix, variant_dir, p: (builder, [rebuilt], "c" * 64))
    opened = []

    class FakeStim:
        @staticmethod
        def Circuit(text):
            opened.append(text)
            return shipped

    class FakeMatcher:
        num_fault_ids = 3

    class FakeMatching:
        seen_dems = []

        @staticmethod
        def from_detector_error_model(dem):
            FakeMatching.seen_dems.append(dem)
            return FakeMatcher()

    class FakePyMatching:
        Matching = FakeMatching

    monkeypatch.setitem(sys.modules, "stim", FakeStim)
    monkeypatch.setitem(sys.modules, "pymatching", FakePyMatching)

    recovered = Z.recover_point(source_root, "grown", "0.0008",
                                rebuilt_oracle=True)

    # the shipped artifact is opened for provenance but never sampled, never
    # DEM'd, and its factual inequality is recorded — never refused.
    assert opened == ["shipped split-generation circuit"]
    assert recovered.circuit is rebuilt
    assert recovered.rebuilt is rebuilt
    assert recovered.flattened_equal is False
    assert rebuilt.dem_calls == 1
    assert shipped.dem_calls == 0
    assert len(FakeMatching.seen_dems) == 1
    provenance = recovered.provenance
    assert provenance["recovery_mode"] == "rebuilt_oracle_r8"
    assert provenance["sampled_circuit_source"] == \
        "rebuilt_from_pinned_driver_builder"
    assert provenance["dem_source"] == "rebuilt_from_pinned_driver_builder"
    assert provenance["shipped_oracle_used"] is False
    assert provenance["shipped_rebuilt_flattened_equal"] is False
    assert provenance["shipped_path"].endswith(
        Z.shipped_filename("0.0008"))
    assert provenance["shipped_blob_sha1"] == \
        Z.MEMBER_PINS["shipped_grown_0.0008"][1]
    assert provenance["dem_error_instructions"] == 1


def test_rebuilt_oracle_rows_validate_and_stay_comparable():
    rows = [_point_row("ungrown", i) for i in range(6)]
    rows += [_point_row("grown", i) for i in range(6)]
    analysis = Z.analyze_points(rows)
    # all twelve cells comparable; the four rebuilt-oracle grown cells are
    # comparable on their own clean recovery despite flattened_equal=False.
    assert all(c["comparable"] for c in analysis["comparisons"])
    split = [c for c in analysis["comparisons"]
             if (c["variant"], str(c["p"])) in
             {("grown", p) for _, p in Z.REBUILT_ORACLE_CELLS}]
    assert len(split) == 4
    for row in rows:
        key = (row["variant"], row["p_label"])
        if key in Z.REBUILT_ORACLE_CELLS:
            assert row["oracle"] == "rebuilt"
            assert row["flattened_equal"] is False
        else:
            assert row["oracle"] == "shipped"
            assert row["flattened_equal"] is True
    # a rebuilt-oracle cell mislabelled shipped must fail validation
    bad = _point_row("grown", 1)  # grown@0.0008 is a split cell
    bad["oracle"] = "shipped"
    bad["flattened_equal"] = True
    with pytest.raises(Z.Refusal, match="oracle"):
        Z._validate_point_rows([bad])


def test_smoke_wires_builder_only_p0_and_shipped_grown_oracle(monkeypatch):
    calls = []
    rebuilt0 = object()
    rec0 = Z.RecoveredPoint(
        "ungrown", "0", rebuilt0, [0, 2], None, rebuilt0, None,
        {
            "variant": "ungrown", "p": 0.0, "p_label": "0",
            "recovery_mode": "builder_only_p0",
            "sampled_circuit_source": "rebuilt_from_pinned_driver_builder",
            "dem_source": "rebuilt_from_pinned_driver_builder",
            "shipped_oracle_used": False,
            "shipped_rebuilt_flattened_equal": None,
            "dem_error_instructions": 0,
        })
    rec1 = Z.RecoveredPoint(
        "grown", "0.001", object(), [0, 2], object(), object(), True,
        {
            "variant": "grown", "p": 0.001, "p_label": "0.001",
            "recovery_mode": "shipped_equality_required",
            "sampled_circuit_source": "shipped_author_circuit",
            "dem_source": "shipped_author_circuit",
            "shipped_oracle_used": True,
            "shipped_rebuilt_flattened_equal": True,
            "dem_error_instructions": 1,
        })

    def fake_recover(source_root, variant, p_label, *, builder_only=False):
        calls.append((variant, p_label, builder_only))
        return rec0 if builder_only else rec1

    def fake_counts(circuit, mask, matcher, shots, chunk, seed,
                    crosscheck_accepted_rows=False, _sampler_factory=None):
        result = {
            "accepted": shots, "errors": 0,
            "chunk_records": [{"seed": seed, "shots": shots}],
            "crosscheck": ({"rows": shots, "mismatches": 0}
                           if crosscheck_accepted_rows else None),
        }
        if matcher is None:
            result.update({
                "detector_bits_fired": 0,
                "actual_observable_flip_rows": 0,
            })
        return result

    monkeypatch.setattr(Z, "recover_point", fake_recover)
    monkeypatch.setattr(Z, "sample_counts", fake_counts)
    gates, evidence = Z.run_smoke_gates(Path("/unused"))

    assert calls == [
        ("ungrown", "0", True),
        ("grown", "0.001", False),
    ]
    assert [gate["gate"] for gate in gates] == [
        "ungrown_noiseless_builder_only_rebuilt",
        "ungrown_noiseless_rebuilt_dem_zero_errors",
        "ungrown_noiseless_rebuilt_zero_detector_observable_bits",
        "ungrown_noiseless_rebuilt_same_seed_repeat_bitwise",
        "grown_0.001_shipped_rebuilt_flattened_equal",
        "grown_0.001_dem_fault_ids",
        "grown_0.001_scalar_batch_bitwise",
        "grown_0.001_same_seed_repeat_bitwise",
    ]
    assert all(gate["pass"] for gate in gates)
    assert evidence["builder_only_p0_recovery"][
        "shipped_rebuilt_flattened_equal"] is None
    assert evidence["grown_author_sampled_recovery"][
        "shipped_rebuilt_flattened_equal"] is True


# ------------------------------------------- vectorized counting with fakes

class FakeSampler:
    """Duck sampler: replays preloaded chunks, records requested sizes."""

    def __init__(self, syndrome_chunks, actual_chunks):
        self._syn = [np.array(s, dtype=bool) for s in syndrome_chunks]
        self._act = [np.array(a, dtype=bool) for a in actual_chunks]
        self.requested = []

    def sample(self, n, separate_observables=True):
        self.requested.append(n)
        syn, act = self._syn.pop(0), self._act.pop(0)
        assert len(syn) == n
        return syn, act


class FakeMatcher:
    """Vectorized zero predictor; records exactly the rows it is given."""

    def __init__(self, fail_batch=False, flip=False):
        self.seen_rows = []
        self._fail = fail_batch
        self._flip = flip

    def decode_batch(self, syndrome):
        self.seen_rows.append(syndrome.shape[0])
        if self._fail:
            raise ValueError("boom")
        out = np.zeros((syndrome.shape[0], 3), dtype=np.uint8)
        if self._flip:
            out[:, 0] = 1
        return out

    def decode(self, row):
        out = np.zeros(3, dtype=np.uint8)
        if self._flip:
            out[0] = 1
        return out


def _synthetic_chunks(rows_per_chunk, detectors=6):
    rng = np.random.default_rng(20260903)
    syn = [rng.random((n, detectors)) < 0.5 for n in rows_per_chunk]
    act = [rng.random((n, 3)) < 0.02 for n in rows_per_chunk]
    return syn, act


def test_sample_counts_matches_bruteforce_and_decodes_accepted_only():
    mask = [1, 4]
    syn, act = _synthetic_chunks([7, 5])
    matcher = FakeMatcher()
    sampler = FakeSampler(syn, act)
    counts = Z.sample_counts(object(), mask, matcher, 12, 7, seed=99,
                             crosscheck_accepted_rows=True,
                             _sampler_factory=lambda *_: sampler)
    syndrome = np.concatenate([np.array(s) for s in syn])
    actual = np.concatenate([np.array(a) for a in act])
    acc = ~syndrome[:, mask].any(axis=1)
    exp_accepted = int(acc.sum())
    exp_errors = int((acc & actual.any(axis=1)).sum())
    assert counts["accepted"] == exp_accepted
    assert counts["errors"] == exp_errors
    assert counts["discarded"] == 12 - exp_accepted
    assert counts["acceptance"] == exp_accepted / 12
    assert counts["ler"] == exp_errors / exp_accepted
    # the decoder saw ONLY the accepted rows, per chunk
    per_chunk_acc = [int((~np.array(s)[:, mask].any(axis=1)).sum())
                     for s in syn]
    assert matcher.seen_rows == per_chunk_acc
    assert counts["crosscheck"] == {"rows": exp_accepted, "mismatches": 0}
    # fixed chunk schedule was honoured (identity component)
    assert sampler.requested == [7, 5]
    assert counts["schedule"]["num_chunks"] == 2
    assert counts["schedule"]["final_chunk"] == 5


def test_sample_counts_decode_batch_valueerror_is_a_refusal():
    syn, act = _synthetic_chunks([6])
    with pytest.raises(Z.Refusal, match="decode_batch ValueError"):
        Z.sample_counts(object(), [0], FakeMatcher(fail_batch=True), 6, 6,
                        seed=1,
                        _sampler_factory=lambda *_: FakeSampler(syn, act))


def test_sample_counts_noiseless_arm_requires_zero_detectors_and_observables():
    zeros_syn = [np.zeros((4, 5), dtype=bool)]
    zeros_act = [np.zeros((4, 3), dtype=bool)]
    counts = Z.sample_counts(object(), [0, 2], None, 4, 4, seed=5,
                             _sampler_factory=lambda *_:
                             FakeSampler(zeros_syn, zeros_act))
    assert counts["accepted"] == 4 and counts["errors"] == 0
    assert counts["actual_observable_flip_rows"] == 0
    assert counts["detector_bits_fired"] == 0
    assert counts["chunk_records"][0]["detector_bits_fired"] == 0
    assert counts["chunk_records"][0]["actual_observable_flip_rows"] == 0
    dirty = [np.zeros((4, 5), dtype=bool)]
    dirty[0][2, 3] = True
    with pytest.raises(Z.Refusal, match="noiseless arm"):
        Z.sample_counts(object(), [0, 2], None, 4, 4, seed=5,
                        _sampler_factory=lambda *_:
                        FakeSampler(dirty, zeros_act))
    flipped_actual = [np.zeros((4, 3), dtype=bool)]
    flipped_actual[0][1, 2] = True
    with pytest.raises(Z.Refusal, match="actual-observable flips"):
        Z.sample_counts(object(), [0, 2], None, 4, 4, seed=5,
                        _sampler_factory=lambda *_:
                        FakeSampler(zeros_syn, flipped_actual))


def test_chunk_schedule_rule_and_identity():
    assert Z.chunk_schedule(10_000_000, 10_000) == {
        "total": 10_000_000, "chunk": 10_000, "num_chunks": 1000,
        "final_chunk": 10_000,
        "rule": "every chunk == 'chunk' except the final remainder"}
    assert Z.chunk_schedule(20_000, 5_000)["num_chunks"] == 4
    odd = Z.chunk_schedule(7, 3)
    assert odd["num_chunks"] == 3 and odd["final_chunk"] == 1
    assert Z.chunk_schedule(71_000_000, 10_000) != Z.chunk_schedule(
        71_000_000, 5_000)


# ------------------------------------------------------------- tar handling

ROOT = f"Zero-level_CCZ_Distillation-{Z.COMMIT}/"


def _fake_pins():
    contents = {
        "driver": b"error_rate = [1e-4]\nprint('x')\n",
        "shipped": b"QUBIT_COORDS(0, 0) 0\n",
    }
    pins = {
        "driver_ungrown": ("stim/ungrown/4_surface3_d3.py",
                           Z.git_blob_sha1(contents["driver"])),
        "shipped_ungrown_0": ("stim/ungrown/832_text_0.stim",
                              Z.git_blob_sha1(contents["shipped"])),
    }
    return pins, contents


def _tar_bytes(members):
    bio = io.BytesIO()
    with tarfile.open(fileobj=bio, mode="w") as tf:
        for name, data, kind in members:
            info = tarfile.TarInfo(name)
            if kind == "dir":
                info.type = tarfile.DIRTYPE
                tf.addfile(info)
            elif kind == "symlink":
                info.type = tarfile.SYMTYPE
                info.linkname = "/etc/passwd"
                tf.addfile(info)
            else:
                info.size = len(data)
                tf.addfile(info, io.BytesIO(data))
    return bio.getvalue()


def _good_members(pins, contents):
    return [
        (ROOT + "stim/ungrown/", None, "dir"),
        (ROOT + "stim/ungrown/4_surface3_d3.py", contents["driver"], "file"),
        (ROOT + "stim/ungrown/832_text_0.stim", contents["shipped"], "file"),
        (ROOT + "plot_LER_suc/CCZ_plot.py", b"unpinned repo member\n", "file"),
        ("README.md", b"outside the repo root\n", "file"),
    ]


def test_extraction_writes_only_pinned_members_and_is_idempotent(tmp_path):
    pins, contents = _fake_pins()
    blob = _tar_bytes(_good_members(pins, contents))
    dest = tmp_path / "source"
    report = Z.extract_pinned_members(io.BytesIO(blob), dest, pins)
    rels = sorted(w["path"] for w in report["written"])
    assert rels == ["stim/ungrown/4_surface3_d3.py",
                    "stim/ungrown/832_text_0.stim"]
    assert (dest / Z.REPO_ROOT_NAME / "stim/ungrown/832_text_0.stim"
            ).read_bytes() == contents["shipped"]
    assert not (dest / Z.REPO_ROOT_NAME / "plot_LER_suc").exists()
    report2 = Z.extract_pinned_members(io.BytesIO(blob), dest, pins)
    assert report2["written"] == []
    assert len(report2["kept_idempotent"]) == 2


def test_tar_traversal_and_bad_members_are_rejected(tmp_path):
    pins, contents = _fake_pins()
    pinned = (ROOT + "stim/ungrown/4_surface3_d3.py", contents["driver"],
              "file")
    bad_names = [
        (ROOT + "../evil.py", b"x", "file"),
        ("/abs/evil.py", b"x", "file"),
        (ROOT + "stim/ungrown/link", b"", "symlink"),
    ]
    for member in bad_names:
        blob = _tar_bytes([pinned, member])
        with pytest.raises(Z.Refusal):
            Z.extract_pinned_members(io.BytesIO(blob), tmp_path / "d1", pins)
    dup = _tar_bytes([(ROOT + "stim/ungrown/832_text_0.stim",
                       contents["shipped"], "file")] * 2)
    with pytest.raises(Z.Refusal, match="duplicate"):
        Z.extract_pinned_members(io.BytesIO(dup), tmp_path / "d2", pins)
    wrong = _tar_bytes([(ROOT + "stim/ungrown/832_text_0.stim",
                         b"tampered\n", "file")])
    with pytest.raises(Z.Refusal, match="blob sha1"):
        Z.extract_pinned_members(io.BytesIO(wrong), tmp_path / "d3", pins)
    for d in ("d1", "d2", "d3"):
        assert not (tmp_path / d).exists() or \
            not any((tmp_path / d).rglob("*"))


def test_extraction_refuses_missing_pins_and_differing_existing(tmp_path):
    pins, contents = _fake_pins()
    partial = _tar_bytes([(ROOT + "stim/ungrown/4_surface3_d3.py",
                           contents["driver"], "file")])
    with pytest.raises(Z.Refusal, match="absent from tar"):
        Z.extract_pinned_members(io.BytesIO(partial), tmp_path / "m", pins)
    blob = _tar_bytes(_good_members(pins, contents))
    dest = tmp_path / "n"
    Z.extract_pinned_members(io.BytesIO(blob), dest, pins)
    victim = dest / Z.REPO_ROOT_NAME / "stim/ungrown/832_text_0.stim"
    victim.write_bytes(b"locally mutated\n")
    with pytest.raises(Z.Refusal, match="refusing to overwrite"):
        Z.extract_pinned_members(io.BytesIO(blob), dest, pins)
    assert victim.read_bytes() == b"locally mutated\n"


def test_source_tar_gates_and_run_dir_nesting(tmp_path):
    fake = tmp_path / "tar.gz"
    fake.write_bytes(b"not the advisory tar")
    with pytest.raises(Z.Refusal, match="size"):
        Z.validate_source_tar(fake)
    same_size = tmp_path / "same_size.bin"
    same_size.write_bytes(b"\0" * Z.TAR_SIZE)
    with pytest.raises(Z.Refusal, match="sha256"):
        Z.validate_source_tar(same_size)
    # a source tar inside the run dir would risk provenance self-overwrite
    run_root = tmp_path / "run"
    nested = run_root / "source" / "source.tar.gz"
    nested.parent.mkdir(parents=True)
    nested.write_bytes(b"x")
    with pytest.raises(Z.Refusal, match="inside the run dir"):
        Z.refuse_tar_inside_run_dir(nested, run_root)
    outside = tmp_path / "outside.tar.gz"
    outside.write_bytes(b"x")
    Z.refuse_tar_inside_run_dir(outside, run_root)  # no refusal outside


def test_prepare_uses_one_validated_tar_capture_for_extract_and_archive(
        tmp_path, monkeypatch):
    pins, contents = _fake_pins()
    captured = _tar_bytes(_good_members(pins, contents))
    replacement = b"path replaced after immutable capture"
    tar_path = tmp_path / "source.tar"
    tar_path.write_bytes(captured)
    prereg = _prereg(tmp_path)
    target_dir = tmp_path / "physics" / "qlops"
    monkeypatch.setattr(Z, "TARGET_DIR", target_dir)
    manifest = _canonical_manifest(prereg, uuid8="facefeed")
    run_dir = _write_campaign_run(target_dir, manifest)
    monkeypatch.setattr(Z, "TAR_SIZE", len(captured))
    monkeypatch.setattr(Z, "TAR_SHA256", Z.sha256_bytes(captured))
    monkeypatch.setattr(Z, "MEMBER_PINS", pins)
    monkeypatch.setattr(
        Z, "validate_env",
        lambda: {"python": "3.test", "platform": "test", "machine": "test",
                 **Z.REQUIRED_VERSIONS})
    source_reads = []
    real_read_bytes = Path.read_bytes

    def count_source_reads(path):
        if path.resolve() == tar_path.resolve():
            source_reads.append(path)
        return real_read_bytes(path)

    monkeypatch.setattr(Path, "read_bytes", count_source_reads)
    real_extract = Z.extract_pinned_members

    def replace_path_then_extract(tar_source, dest_root, active_pins):
        tar_path.write_bytes(replacement)
        return real_extract(tar_source, dest_root, active_pins)

    monkeypatch.setattr(Z, "extract_pinned_members",
                        replace_path_then_extract)
    ctx = Z.prepare_run(run_dir, tar_path, prereg)
    extracted = (run_dir / "source" / Z.REPO_ROOT_NAME
                 / "stim/ungrown/4_surface3_d3.py")
    archived = run_dir / "source" / "source_tar.gz"
    assert len(source_reads) == 1
    assert real_read_bytes(tar_path) == replacement
    assert extracted.read_bytes() == contents["driver"]
    assert archived.read_bytes() == captured
    assert ctx["source_tar_bytes"] == captured
    assert ctx["pristine"]["sha256"] == Z.sha256_bytes(captured)


# ------------------------------------------------------ campaign manifest

def _prereg(tmp_path):
    path = tmp_path / "pre_statement.md"
    path.write_text(
        "# pre_statement.md\n\n"
        "## Revision 5 — zero-level original\nold prereg body\n\n"
        "## Revision 6 — zero-level oracle replacement\nreplacement body\n\n"
        "## Revision 7 — zero-level writer fix\nprereg body\n\n"
        "## Revision 8 — rebuilt-oracle amendment for the F-Z4 split cells\n"
        "amendment body\n",
        encoding="utf-8")
    return path


def _canonical_manifest(prereg_path, *, stamp="20260903T000000Z",
                        uuid8="deadbeef", agent="tester",
                        created_utc="2026-09-03T00:00:00Z"):
    prereg_sha = Z.sha256_file(prereg_path)
    binder = "\x1f".join((Z.GATE, agent, prereg_sha, stamp, uuid8))
    run_id = f"{stamp}_{uuid8}_{Z.sha256_bytes(binder.encode())[:12]}"
    return {
        "run_id": run_id,
        "target": "physics/qlops",
        "gate": Z.GATE,
        "agent": agent,
        "created_utc": created_utc,
        "prereg_sha256": prereg_sha,
        "status": "RUNNING",
    }


def _write_campaign_run(target_dir, manifest):
    run_dir = target_dir / "campaigns" / manifest["run_id"]
    run_dir.mkdir(parents=True, exist_ok=False)
    (run_dir / "manifest.json").write_text(
        json.dumps(manifest), encoding="utf-8")
    return run_dir


def test_canonical_manifest_and_live_campaign_gates(tmp_path, monkeypatch):
    prereg = _prereg(tmp_path)
    target_dir = tmp_path / "physics" / "qlops"
    monkeypatch.setattr(Z, "TARGET_DIR", target_dir)
    canonical = _canonical_manifest(prereg)
    run_dir = _write_campaign_run(target_dir, canonical)
    manifest_path = run_dir / "manifest.json"
    captured = manifest_path.read_bytes()
    manifest = Z.validate_run_dir(run_dir, prereg)
    assert manifest["gate"] == Z.GATE
    Z.assert_run_live(run_dir, captured)

    def refused(mutation):
        bad = dict(canonical)
        for key in mutation.get("_drop", ()):
            bad.pop(key, None)
        bad.update({key: value for key, value in mutation.items()
                    if key != "_drop"})
        manifest_path.write_text(json.dumps(bad), encoding="utf-8")
        with pytest.raises(Z.Refusal):
            Z.validate_run_dir(run_dir, prereg)
        manifest_path.write_bytes(captured)

    refused({"gate": "other-gate"})
    refused({"status": "FROZEN-CERTIFIED"})
    refused({"status": "CLOSED"})
    refused({"run_id": "20260903T000000Z_deadbeef"})
    refused({"prereg_sha256": "0" * 64})
    refused({"target": "physics/msd"})
    refused({"created_utc": "2026-09-03 00:00:00"})
    refused({"extra_key": True})
    refused({"_drop": ["gate"]})
    refused({"_drop": ["status"]})

    wrong_parent = tmp_path / "elsewhere" / canonical["run_id"]
    wrong_parent.mkdir(parents=True)
    (wrong_parent / "manifest.json").write_bytes(captured)
    with pytest.raises(Z.Refusal, match="exact campaigns parent"):
        Z.validate_run_dir(wrong_parent, prereg)

    wrong_basename = target_dir / "campaigns" / "wrong-basename"
    wrong_basename.mkdir()
    (wrong_basename / "manifest.json").write_bytes(captured)
    with pytest.raises(Z.Refusal, match="basename"):
        Z.validate_run_dir(wrong_basename, prereg)

    bad_binder = dict(canonical)
    bad_binder["run_id"] = canonical["run_id"][:-1] + (
        "0" if canonical["run_id"][-1] != "0" else "1")
    bad_binder_dir = _write_campaign_run(target_dir, bad_binder)
    with pytest.raises(Z.Refusal, match="does not bind"):
        Z.validate_run_dir(bad_binder_dir, prereg)

    refused({"created_utc": "2026-09-03T00:00:01Z"})

    (run_dir / "sha256s.txt").write_text("frozen\n", encoding="utf-8")
    with pytest.raises(Z.Refusal, match="sha256s.txt"):
        Z.validate_run_dir(run_dir, prereg)
    (run_dir / "sha256s.txt").unlink()
    (run_dir / "status.json").write_text("{}\n", encoding="utf-8")
    with pytest.raises(Z.Refusal, match="status.json"):
        Z.validate_run_dir(run_dir, prereg)
    (run_dir / "status.json").unlink()

    manifest_path.write_bytes(captured + b"\n")
    with pytest.raises(Z.Refusal, match="manifest bytes changed"):
        Z.assert_run_live(run_dir, captured)
    manifest_path.write_bytes(captured)

    plain = tmp_path / "plain.md"
    plain.write_text("no revision marker\n", encoding="utf-8")
    plain_run = _write_campaign_run(target_dir, _canonical_manifest(plain))
    with pytest.raises(Z.Refusal, match="Revision 8"):
        Z.validate_run_dir(plain_run, plain)
    missing = target_dir / "campaigns" / "does_not_exist"
    with pytest.raises(Z.Refusal, match="never created"):
        Z.validate_run_dir(missing, prereg)
    assert not missing.exists()


# ------------------------------------------- full/smoke wiring and CLI codes

class _StubRecovered:
    circuit = None
    mask = []
    matcher = None

    def __init__(self, variant, p_label, identity, source_root):
        # factual F-Z4 inequality on the Revision-8 rebuilt-oracle cells
        self.flattened_equal = (variant, p_label) not in Z.REBUILT_ORACLE_CELLS
        self.provenance = _recovery(
            variant, p_label, identity, source_root)


def _full_context(run_dir, prereg, manifest):
    prereg_bytes = prereg.read_bytes()
    identity = _base_identity(
        manifest["run_id"], manifest["prereg_sha256"])
    return {
        "identity": identity,
        "source_root": run_dir / "source",
        "prereg_path": prereg,
        "prereg_bytes": prereg_bytes,
        "manifest_bytes": (run_dir / "manifest.json").read_bytes(),
    }


def test_smoke_orchestration_creates_results_before_atomic_writes(
        tmp_path, monkeypatch):
    prereg = _prereg(tmp_path)
    target_dir = tmp_path / "physics" / "qlops"
    monkeypatch.setattr(Z, "TARGET_DIR", target_dir)
    manifest = _canonical_manifest(prereg)
    run_dir = _write_campaign_run(target_dir, manifest)
    ctx = _full_context(run_dir, prereg, manifest)
    monkeypatch.setattr(
        Z, "run_smoke_gates",
        lambda _: ([{"gate": "stub", "pass": True, "detail": {}}],
                   {"counts1": {"crosscheck": {"rows": 1, "mismatches": 0}}}))

    report = Z.run_smoke(run_dir, ctx)

    assert report["smoke_verdict"] == "PASS"
    assert (run_dir / "results" / "smoke.json").is_file()
    assert (run_dir / "results" / "summary.md").is_file()
    assert (run_dir / "results" / "inventory.json").is_file()


def test_full_mode_runs_smoke_then_writes_and_resumes_complete_rows(
        tmp_path, monkeypatch):
    order = []
    prereg = _prereg(tmp_path)
    target_dir = tmp_path / "physics" / "qlops"
    monkeypatch.setattr(Z, "TARGET_DIR", target_dir)
    manifest = _canonical_manifest(prereg)
    run_dir = _write_campaign_run(target_dir, manifest)
    ctx = _full_context(run_dir, prereg, manifest)

    def fake_smoke(source_root):
        order.append("smoke")
        return [{"gate": "stub", "pass": True, "detail": {}}], {
            "counts1": {"crosscheck": {"rows": 5, "mismatches": 0}}}

    monkeypatch.setattr(Z, "run_smoke_gates", fake_smoke)

    def fake_recover(source_root, variant, p_label, *, rebuilt_oracle=False):
        order.append(f"recover:{variant}:{p_label}")
        assert rebuilt_oracle == (
            (variant, p_label) in Z.REBUILT_ORACLE_CELLS)
        return _StubRecovered(
            variant, p_label, ctx["identity"], source_root)

    # Unique seed -> exact (variant, p_label) author row; shots alone is
    # ambiguous across the 5m/10m cells.
    seed_to_point = {
        seed: (variant, p_label, Z._author_row(variant, p_label))
        for variant, table in Z.FULL_SEEDS.items()
        for p_label, seed in table.items()
    }

    def fake_sample_counts(circuit, mask, matcher, shots, chunk, seed,
                           crosscheck_accepted_rows=False,
                           _sampler_factory=None):
        order.append(f"sample:{seed}")
        variant, p_label, (_label, _shots, ler, acc) = \
            seed_to_point[int(seed)]
        accepted = int(shots * acc)
        errors = int(round(ler * accepted))
        return {
            "status": "ok", "reason": None, "seed": seed, "shots": shots,
            "schedule": Z.chunk_schedule(shots, chunk),
            "accepted": accepted, "errors": errors,
            "discarded": shots - accepted,
            "acceptance": accepted / shots,
            "ler": errors / accepted,
            "chunk_records": _chunk_records(
                variant, p_label, shots, accepted, errors),
            "crosscheck": None,
        }

    monkeypatch.setattr(Z, "recover_point", fake_recover)
    monkeypatch.setattr(Z, "sample_counts", fake_sample_counts)
    analysis = Z.run_full(run_dir, ctx)
    assert order[0] == "smoke"
    assert order[1] == "recover:ungrown:0.001"
    assert sum(1 for item in order if item.startswith("sample:")) == 12
    assert sum(1 for item in order if item.startswith("recover:")) == 12
    assert analysis["candidate_verdict"] == "FROZEN-CERTIFIED"
    assert analysis["smoke_gates"] == [
        {"gate": "stub", "pass": True, "detail": {}}]
    assert analysis["smoke_evidence"]["counts1"]["crosscheck"] == \
        {"rows": 5, "mismatches": 0}
    assert (run_dir / "results" / "analysis.json").is_file()

    # A second invocation accepts every complete row without recovery/sample.
    order.clear()
    resumed = Z.run_full(run_dir, ctx)
    assert order == ["smoke"]
    assert resumed["candidate_verdict"] == "FROZEN-CERTIFIED"

    # A failing smoke battery stops full mode before directories or sampling.
    order.clear()

    def failing_smoke(source_root):
        order.append("smoke")
        return [{"gate": "stub", "pass": False, "detail": {}}], {}

    monkeypatch.setattr(Z, "run_smoke_gates", failing_smoke)
    manifest2 = _canonical_manifest(prereg, uuid8="cafebabe")
    run_dir2 = _write_campaign_run(target_dir, manifest2)
    ctx2 = _full_context(run_dir2, prereg, manifest2)
    with pytest.raises(Z.Refusal, match="smoke gates failed inside full"):
        Z.run_full(run_dir2, ctx2)
    assert order == ["smoke"]
    assert not (run_dir2 / "results" / "points").exists()
    assert not (run_dir2 / "results" / "analysis.json").exists()


def test_full_resume_refuses_scalar_only_row_before_sampling(
        tmp_path, monkeypatch):
    prereg = _prereg(tmp_path)
    target_dir = tmp_path / "physics" / "qlops"
    monkeypatch.setattr(Z, "TARGET_DIR", target_dir)
    manifest = _canonical_manifest(prereg, uuid8="0123abcd")
    run_dir = _write_campaign_run(target_dir, manifest)
    ctx = _full_context(run_dir, prereg, manifest)
    full_row = _point_row(
        "grown", 5, identity_base=ctx["identity"],
        source_root=ctx["source_root"])
    scalar_only = dict(full_row, recovery={}, chunk_records=[])
    points = run_dir / "results" / "points"
    points.mkdir(parents=True)
    (points / "grown_0.0001.json").write_text(
        json.dumps(scalar_only), encoding="utf-8")
    sample_calls = []
    monkeypatch.setattr(
        Z, "run_smoke_gates",
        lambda _source: ([{"gate": "stub", "pass": True, "detail": {}}], {}))
    monkeypatch.setattr(
        Z, "recover_point",
        lambda *args: pytest.fail("resume validation must precede recovery"))

    def unexpected_sample(*args, **kwargs):
        sample_calls.append(1)
        pytest.fail("resume validation must precede sampling")

    monkeypatch.setattr(Z, "sample_counts", unexpected_sample)
    with pytest.raises(Z.Refusal, match="recovery fields"):
        Z.run_full(run_dir, ctx)
    assert sample_calls == []


def test_cli_return_codes(tmp_path, monkeypatch):
    # argparse misuse -> 2
    with pytest.raises(SystemExit) as exc:
        Z.main(["--source-tar", str(tmp_path / "t.gz")])
    assert exc.value.code == 2
    # refusal (absent run dir) -> 1, REJECTED on stderr
    rc = Z.main(["--source-tar", str(tmp_path / "t.gz"),
                 "--run-dir", str(tmp_path / "nope"), "--mode", "smoke"])
    assert rc == 1
    assert not (tmp_path / "nope").exists()          # nothing was created
    # success path -> 0
    monkeypatch.setattr(Z, "prepare_run",
                        lambda *a, **k: {"identity": {}, "source_root": None})
    monkeypatch.setattr(Z, "run_smoke", lambda *a, **k: {"smoke_verdict":
                                                         "PASS"})
    rc = Z.main(["--source-tar", str(tmp_path / "t.gz"),
                 "--run-dir", str(tmp_path / "run"), "--mode", "smoke"])
    assert rc == 0
    # refusal inside the mode runner -> 1
    def refusing_smoke(*a, **k):
        raise Z.Refusal("smoke semantic gates failed: ['x']")
    monkeypatch.setattr(Z, "run_smoke", refusing_smoke)
    rc = Z.main(["--source-tar", str(tmp_path / "t.gz"),
                 "--run-dir", str(tmp_path / "run"), "--mode", "smoke"])
    assert rc == 1


def test_validate_run_dir_prereg_byte_gates(tmp_path, monkeypatch):
    target_dir = tmp_path / "physics" / "qlops"
    monkeypatch.setattr(Z, "TARGET_DIR", target_dir)
    prereg = _prereg(tmp_path)
    good = prereg.read_bytes()
    run_dir = _write_campaign_run(target_dir, _canonical_manifest(prereg))
    with pytest.raises(Z.Refusal, match="captured prereg bytes"):
        Z.validate_run_dir(run_dir, good + b"tail")

    rev9 = tmp_path / "rev9.md"
    rev9.write_bytes(good + b"\n## Revision 9 (\n")
    rev9_run = _write_campaign_run(
        target_dir, _canonical_manifest(rev9, uuid8="00000009"))
    with pytest.raises(Z.Refusal, match="expected exactly 8"):
        Z.validate_run_dir(rev9_run, rev9.read_bytes())

    bad = tmp_path / "bad_order.md"
    bad.write_bytes(good + b"\n## Revision 4 (\n")
    bad_run = _write_campaign_run(
        target_dir, _canonical_manifest(bad, uuid8="00000004"))
    with pytest.raises(Z.Refusal, match="strictly increasing"):
        Z.validate_run_dir(bad_run, bad.read_bytes())


def test_snapshot_uses_captured_prereg_bytes(tmp_path):
    prereg = _prereg(tmp_path)
    captured = prereg.read_bytes()
    run_dir = tmp_path / "snap"
    run_dir.mkdir()
    Z.snapshot_code_and_prereg(run_dir, captured)
    snap = run_dir / "provenance" / "snapshots" / "pre_statement.md"
    assert snap.read_bytes() == captured
    prereg.write_bytes(captured + b"mutated after capture")
    Z.snapshot_code_and_prereg(run_dir, captured)  # idempotent vs captured
    assert snap.read_bytes() == captured


def test_analyze_points_row_set_gates():
    rows = [_point_row("ungrown", i) for i in range(6)]
    rows += [_point_row("grown", i) for i in range(6)]
    assert Z.analyze_points(rows)["n_comparisons"] == 24
    with pytest.raises(Z.Refusal, match="12 unique"):
        Z.analyze_points([])
    with pytest.raises(Z.Refusal, match="12 unique"):
        Z.analyze_points(rows[:-1])
    with pytest.raises(Z.Refusal, match="duplicate"):
        Z.analyze_points(rows + [dict(rows[0])])
    corrupt = [dict(r) for r in rows]
    corrupt[3]["errors"] = corrupt[3]["accepted"] + 1
    with pytest.raises(Z.Refusal, match="invariant"):
        Z.analyze_points(corrupt)
    badid = [dict(r) for r in rows]
    badid[2] = dict(badid[2], identity=dict(badid[2]["identity"], seed=999))
    with pytest.raises(Z.Refusal, match="identity"):
        Z.analyze_points(badid)
    bad_recovery_hash = list(rows)
    bad_recovery_hash[0] = json.loads(json.dumps(rows[0]))
    bad_recovery_hash[0]["recovery"]["prefix_exec_sha256"] = "f" * 63
    with pytest.raises(Z.Refusal, match="prefix_exec_sha256 malformed"):
        Z.analyze_points(bad_recovery_hash)
    bad_recovery_pin = list(rows)
    bad_recovery_pin[0] = json.loads(json.dumps(rows[0]))
    bad_recovery_pin[0]["recovery"]["driver_blob_sha1"] = "0" * 40
    with pytest.raises(Z.Refusal, match="recovery pins"):
        Z.analyze_points(bad_recovery_pin)
    bad_recovery = list(rows)
    bad_recovery[0] = json.loads(json.dumps(rows[0]))
    bad_recovery[0]["recovery"].pop("dem_sha256")
    with pytest.raises(Z.Refusal, match="recovery fields"):
        Z.analyze_points(bad_recovery)
    short_ledger = list(rows)
    short_ledger[0] = json.loads(json.dumps(rows[0]))
    short_ledger[0]["chunk_records"].pop()
    with pytest.raises(Z.Refusal, match="ledger length"):
        Z.analyze_points(short_ledger)
    bad_digest = list(rows)
    bad_digest[0] = json.loads(json.dumps(rows[0]))
    bad_digest[0]["chunk_records"][0]["syndrome_obs_sha256"] = "not-a-hash"
    with pytest.raises(Z.Refusal, match="digest malformed"):
        Z.analyze_points(bad_digest)
    bad_sum = list(rows)
    bad_sum[0] = json.loads(json.dumps(rows[0]))
    bad_sum[0]["chunk_records"][0]["accepted"] -= 1
    with pytest.raises(Z.Refusal, match="ledger sums"):
        Z.analyze_points(bad_sum)
    bad_reason = list(rows)
    bad_reason[0] = dict(rows[0], reason="scalar-only excuse")
    with pytest.raises(Z.Refusal, match="row/identity"):
        Z.analyze_points(bad_reason)


def test_smoke_summary_candidate_line_is_pass(tmp_path):
    payload = {"mode": "smoke", "run_id": "r", "ended_utc": "e",
               "smoke_verdict": "PASS"}
    out = Z.write_summary(tmp_path, payload)
    assert "PASS" in out.read_text()


if __name__ == "__main__":
    import inspect
    import tempfile
    fns = [(k, v) for k, v in sorted(globals().items())
           if k.startswith("test_")]
    failures = 0
    for name, fn in fns:
        params = inspect.signature(fn).parameters
        kwargs = {}
        mp = None
        if "tmp_path" in params:
            td = tempfile.mkdtemp()
            kwargs["tmp_path"] = Path(td)
        if "monkeypatch" in params:
            class _MP:
                def __init__(self):
                    self._undo = []
                def setattr(self, target, name, value):
                    old = getattr(target, name)
                    setattr(target, name, value)
                    self._undo.append(("attr", target, name, old, True))
                def setitem(self, target, key, value):
                    present = key in target
                    old = target.get(key)
                    target[key] = value
                    self._undo.append(("item", target, key, old, present))
                def undo(self):
                    for kind, target, key, old, present in reversed(self._undo):
                        if kind == "attr":
                            setattr(target, key, old)
                        elif present:
                            target[key] = old
                        else:
                            target.pop(key, None)
            mp = _MP()
            kwargs["monkeypatch"] = mp
        try:
            fn(**kwargs)
        except Exception as exc:  # noqa: BLE001 - report and continue
            failures += 1
            print(f"FAIL {name}: {exc}")
        finally:
            if mp is not None:
                mp.undo()
    print(f"{len(fns) - failures}/{len(fns)} tests passed")
    sys.exit(1 if failures else 0)
