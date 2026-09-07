# -*- coding: utf-8 -*-
"""Exact-rational confirmation runner for the global real-versus-unitary
noisy-shadow formulas of arXiv:2608.18935v1 (Hingane & Koh, "Real Classical
Shadows with Noise", posted 2026-08-19).

Scope (preregistered as pre_statement.md Revision 3, appended before any
campaign run): the paper proves analytically (Cor. 3.11, Eq. (44)) that the
variance ratio Var_U/Var_O = (R - x)/(1 - x) >= R >= rho_S(d) > 1 on the whole
admissible beta domain 1 < beta <= d, so the channel-resolved table below is
REGRESSION-GRADE EXACT CONFIRMATION of a known no-flip result -- not a
discovery search.  All 132 unique (noise model, d, p, anchor) cells are
recomputed from two independent formula paths and compared:

  direct path     Eq. (31) invariants t, r, m -> Eq. (33) constants
                  -> second moments E[o^2] = C(a t + 2 r) (Eq. (32), i.e.
                  Eqs. (28)/(29) with the shared m^2 added back)
                  -> second-moment ratio R = E_U/E_O;
                  variances = E[o^2] - m^2 (exactly Eqs. (28)/(29));
                  variance ratio via Eq. (44) with x = m^2/E_O.
  criterion path  Eq. (36) kappa = a_O t/(2 r)
                  -> Eq. (37) ratio (rho_L kappa + rho_S)/(kappa + 1).

Every verdict-path quantity is a fractions.Fraction.  No RNG, no floats, no
matrix allocation anywhere.  Python stdlib only.

Grid (frozen by Revision 3):
  depolarizing       d in {4, 8, 16}; p in {1, 3/4, 1/2, 1/4, 1/10, 1/100,
                     p*(d)} deduplicated per d (the d=4 crossing point
                     p*(4)=3/4 coincides with a ladder value; provenance is
                     carried by the row's "p_roles", not by a duplicate row);
                     3 anchors -> 60 rows.
  amplitude damping  d in {4, 8, 16}; p in {1, 9/10, 4/5, 3/4, 1/2, 1/4,
                     1/10, 1/100}; beta = (1+p)^n with d = 2^n;
                     3 anchors -> 72 rows.
  total              132 unique cells.

Anchors (Eq. (31) invariants (t, r, m), exact for every d):
  maximally_mixed_pauli_z  rho = I/d, O_0 = Z_1:            (d, 1, 0)
  ghz_projector            rho = |GHZ><GHZ|, O_0 = P - I/d:
                                                      ((d-1)/d, ((d-1)/d)^2,
                                                        (d-1)/d)
  ghz_pauli_x              rho = |GHZ><GHZ|, O_0 = X^xn:    (d, 1, 1)

Terminal verdict mapping (frozen by Revision 3, clarified by Revisions 4-7):
  identities invalid                       -> FROZEN-INCONCLUSIVE
  valid and any ratio <= 1                 -> FROZEN-NEGATIVE
                                               (scientific: FLIP-FOUND)
  all 132 valid and both ratios > 1        -> FROZEN-CERTIFIED
                                               (scientific:
                                                NO-FLIP-CONFIRMED-IN-SCOPE)

CLI:
  python3 src/shadows_noisy_variance.py --run-dir DIR
      Runs the 132-cell table into an EXISTING, INITIALIZED, RUNNING campaign
      directory DIR created by scripts/campaign.py.  The runner NEVER creates
      the run directory or manifest.  Before inventory generation it
      snapshots the executed module and the preregistration byte-for-byte into
      DIR (shadows_noisy_variance.snapshot.py, pre_statement.snapshot.md), so
      a campaign freeze binds code and prereg, and writes deterministic
      results.json, summary.json and inventory.json.  The source/prereg
      SHA256 fields in results.json/summary.json are taken from exactly the
      bytes snapshotted, hence match the snapshot copies and the current
      files.  Refusals exit 2 and write nothing.
  python3 src/shadows_noisy_variance.py --check
      Computes the whole table and all regression identities, prints the
      verdict, writes nothing.

Exit codes: 0 run/check completed (verdict recorded in output), 1 completed
but regression identities failed (instrument error -- artifacts are still
written so the failure is auditable), 2 refusal (bad arguments or manifest
state).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from fractions import Fraction
from pathlib import Path

F = Fraction

# --------------------------------------------------------------- frozen grid

DIMENSIONS = (4, 8, 16)

ANCHORS = ("maximally_mixed_pauli_z", "ghz_projector", "ghz_pauli_x")

NOISE_MODELS = ("depolarizing", "amplitude_damping")

DEPOLARIZING_P_LADDER = (
    F(1), F(3, 4), F(1, 2), F(1, 4), F(1, 10), F(1, 100),
)
AMPLITUDE_DAMPING_P_LADDER = (
    F(1), F(9, 10), F(4, 5), F(3, 4), F(1, 2), F(1, 4), F(1, 10), F(1, 100),
)

# SHA256 of the canonical sorted
# "model|d|p|anchor|comma-separated-roles\\n" rows frozen in Revision 5.
# This independent digest makes a same-cardinality ladder substitution fail.
EXPECTED_GRID_SHA256 = "7d573827348828f932bd008836c718c72ec09af12b0257cf24dc5703893e732a"

# ----------------------------------------------------------------- metadata

PAPER = {
    "arxiv_id": "2608.18935",
    "version": "v1",
    "title": "Real Classical Shadows with Noise",
    "authors": "Atharva Hingane, Dax Enshan Koh",
    "posted": "2026-08-19",
    "read_first_hand": "2026-09-04 (arXiv HTML v1; equations transcribed verbatim)",
    "equations": {
        "28": "Var_O(a^) exact orthogonal single-shot variance (real basis, global)",
        "29": "Var_U(a^) exact unitary single-shot variance",
        "33": "second-moment constants C_O, a_O, C_U, a_U",
        "36": "kappa = a_O t / (2 r), the single dimensionless criterion parameter",
        "37": "E_U[o^2]/E_O[o^2] = (rho_L kappa + rho_S)/(kappa + 1); rho_S, rho_L",
        "39": "beta*(d) = d(d^2+7d+8)/(2(d^2+3d+4)), threshold for rho_L > 2",
        "44": "Var_U/Var_O = (R - x)/(1 - x) >= R with x = m^2/E_O[o^2]",
    },
    "supporting": {
        "20": "beta conventions: depolarizing beta-1 = p(d-1); amplitude damping (1+p)^n - 1",
        "31": "invariants t = tr(O_0^2), r = tr(rho O_0^2), m = tr(O_0 rho)",
        "32": "second moments E[o^2] = C (a t + 2 r)",
        "41-43": "rho_S = C_U/C_O, rho_L = C_U a_U/(C_O a_O), rho_L - 2 numerator",
    },
}

ROOT = Path(__file__).resolve().parent.parent
MODULE_PATH = Path(__file__).resolve()
PREREG_PATH = ROOT / "pre_statement.md"
PREREG_REVISION = 7

# Byte-for-byte source loaded for this script invocation. The run snapshots
# these bytes and refuses if the on-disk module diverges before payloads are
# assembled.
MODULE_BYTES_AT_IMPORT = MODULE_PATH.read_bytes()


def _assert_module_stable() -> None:
    if MODULE_PATH.read_bytes() != MODULE_BYTES_AT_IMPORT:
        raise RunnerRefused(
            "REFUSED: shadows_noisy_variance.py changed on disk since this "
            "runner was loaded; rerun on the final source")

MODULE_SNAPSHOT_NAME = "shadows_noisy_variance.snapshot.py"
PREREG_SNAPSHOT_NAME = "pre_statement.snapshot.md"

MANIFEST_NAME = "manifest.json"
EXPECTED_TARGET = "physics/shadows"
EXPECTED_GATE = "noisy-variance-r7"
MANIFEST_STATUS_RUNNING = "RUNNING"

RESULTS_SCHEMA = "shadows-noisy-results/1"
SUMMARY_SCHEMA = "shadows-noisy-summary/1"
INVENTORY_SCHEMA = "shadows-noisy-inventory/1"

MODEL_SECTION = {
    "depolarizing": "arXiv:2608.18935v1 section 3.4.1",
    "amplitude_damping": "arXiv:2608.18935v1 section 3.4.2",
}


class RunnerRefused(Exception):
    """Raised when the CLI must refuse to run (missing/invalid manifest...)."""


# ------------------------------------------------------------------ helpers

def _check_dimension(d) -> int:
    if not isinstance(d, int) or isinstance(d, bool):
        raise ValueError(f"dimension must be an int, got {d!r}")
    if d < 2 or (d & (d - 1)) != 0:
        raise ValueError(f"dimension must be a power of two >= 2, got {d}")
    return d


def _check_probability(p) -> Fraction:
    p = F(p)
    if not (F(0) < p <= F(1)):
        raise ValueError(f"channel parameter p must lie in (0, 1], got {p}")
    return p


def _check_beta(d: int, beta) -> Fraction:
    beta = F(beta)
    if not (F(1) < beta <= d):
        raise ValueError(
            f"beta must lie in the invertible scope (1, d] with d={d}, got {beta}")
    return beta


def _n_qubits(d: int) -> int:
    n = 0
    while (1 << n) < d:
        n += 1
    if (1 << n) != d:
        raise ValueError(f"dimension {d} is not a power of two")
    return n


def frac_text(x: Fraction) -> str:
    x = F(x)
    return f"{x.numerator}/{x.denominator}"


def frac_obj(x: Fraction) -> dict:
    x = F(x)
    return {"numerator": x.numerator, "denominator": x.denominator,
            "text": frac_text(x)}


def to_jsonable(obj):
    """Deterministic JSON projection: every Fraction becomes
    {numerator, denominator, text} with the canonical reduced 'n/d' text."""
    if isinstance(obj, Fraction):
        return frac_obj(obj)
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, bool) or obj is None or isinstance(obj, (int, str)):
        return obj
    raise TypeError(f"unserializable value {obj!r} of type {type(obj)!r}")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(Path(path).read_bytes())


# --------------------------------------------------- paper formulas (exact)

def rho_s(d) -> Fraction:
    """Eq. (37)/(41): rho_S(d) = (d+1)(d+4)/(d+2)^2 = C_U/C_O.

    The bounded-kappa floor of the second-moment ratio; strictly > 1 for
    every d >= 2."""
    d = _check_dimension(d)
    return F((d + 1) * (d + 4), (d + 2) ** 2)


def rho_l(d, beta) -> Fraction:
    """Eq. (37)/(42): rho_L(d, beta) =
    2(d+1)(d+4)(d^2+d-2beta) / ((d+2)^2 (d^2+3d-4beta)).

    The kappa -> infinity ceiling of the second-moment ratio; strictly
    increasing in beta on (1, d], equal to 2 exactly at beta*(d) (Eq. (39))."""
    d = _check_dimension(d)
    beta = F(beta)
    return F(2 * (d + 1) * (d + 4) * (d * d + d - 2 * beta),
             (d + 2) ** 2 * (d * d + 3 * d - 4 * beta))


def beta_star(d) -> Fraction:
    """Eq. (39): beta*(d) = d(d^2+7d+8)/(2(d^2+3d+4)).

    The unique beta with rho_L(d, beta) = 2 (the Eq. (43) numerator vanishes
    there); always strictly inside (1, d)."""
    d = _check_dimension(d)
    return F(d * (d * d + 7 * d + 8), 2 * (d * d + 3 * d + 4))


def depolarizing_beta(d, p) -> Fraction:
    """Section 3.4.1: beta = 1 + p(d-1) for E(Pi_b) = p Pi_b + (1-p) I/d."""
    d = _check_dimension(d)
    p = _check_probability(p)
    return _check_beta(d, F(1) + p * (d - 1))


def depolarizing_p_star(d) -> Fraction:
    """Depolarizing strength whose beta hits Eq. (39): p* = (beta*(d)-1)/(d-1).

    The crossing point where rho_L = 2; p*(4) = 3/4 and p*(8) = 15/23 are
    verbatim in the paper's Section 3.3 numerics."""
    d = _check_dimension(d)
    return _check_probability(F(beta_star(d) - 1, d - 1))


def amplitude_damping_beta(d, p) -> Fraction:
    """Section 3.4.2: one-qubit amplitude damping with Kraus
    K_0 = diag(1, sqrt(p)), K_1 = sqrt(1-p)|0><1|, composed over
    n = log2(d) qubits, gives beta = (1+p)^n."""
    d = _check_dimension(d)
    p = _check_probability(p)
    return _check_beta(d, F(1 + p) ** _n_qubits(d))


def anchor_invariants(d, name: str) -> tuple:
    """Eq. (31) invariants (t, r, m) = (tr(O_0^2), tr(rho O_0^2), tr(O_0 rho))
    for the three frozen anchor (state, observable) pairs.

    Derivations (exact, no matrices):
      maximally_mixed_pauli_z: rho = I/d, O_0 = Z_1.  O_0^2 = I, so
        t = tr(I) = d, r = tr((I/d) I) = 1, m = tr(Z_1 I/d) = tr(Z_1)/d = 0.
      ghz_projector: rho = P = |GHZ_n><GHZ_n|, O_0 = P - I/d (traceless
        part).  With c = 1 - 1/d: O_0|GHZ> = (1 - 1/d)|GHZ> gives m = c;
        O_0^2 = (1 - 2/d) P + I/d^2 gives t = (1 - 2/d) + 1/d = c and
        r = <GHZ|O_0^2|GHZ> = (1 - 1/d)^2 = c^2.
      ghz_pauli_x: rho = |GHZ_n><GHZ_n|, O_0 = X^xn.  O_0^2 = I gives
        t = d, r = 1; |GHZ_n> is a +1 eigenstate of X^xn gives m = 1.
    """
    d = _check_dimension(d)
    c = F(d - 1, d)
    if name == "maximally_mixed_pauli_z":
        return (F(d), F(1), F(0))
    if name == "ghz_projector":
        return (c, c * c, c)
    if name == "ghz_pauli_x":
        return (F(d), F(1), F(1))
    raise ValueError(f"unknown anchor {name!r}; expected one of {ANCHORS}")


def _constants(d: int, beta):
    """Eq. (33): (C_O, a_O, C_U, a_U); requires 1 < beta <= d (validated)."""
    beta = _check_beta(d, beta)
    c_o = F((d - 1) * (d + 2), (d + 4) * (beta - 1))
    a_o = F(d * (d + 3) - 4 * beta, 2 * d * (beta - 1))
    c_u = F(d * d - 1, (d + 2) * (beta - 1))
    a_u = F(d + d * d - 2 * beta, d * (beta - 1))
    return c_o, a_o, c_u, a_u


def second_moment_orthogonal(d, beta, t, r) -> Fraction:
    """E_O[o^2] = C_O (a_O t + 2 r): Eq. (32) LHS, i.e. Eq. (28) plus the
    shared m^2 term (the estimator is unbiased in both ensembles)."""
    d = _check_dimension(d)
    c_o, a_o, _, _ = _constants(d, beta)
    return c_o * (a_o * F(t) + 2 * F(r))


def second_moment_unitary(d, beta, t, r) -> Fraction:
    """E_U[o^2] = C_U (a_U t + 2 r): Eq. (32) LHS, i.e. Eq. (29) plus m^2."""
    d = _check_dimension(d)
    _, _, c_u, a_u = _constants(d, beta)
    return c_u * (a_u * F(t) + 2 * F(r))


def kappa(d, beta, t, r) -> Fraction:
    """Eq. (36): kappa = a_O t/(2 r) = (d(d+3)-4beta) t/(4 d (beta-1) r)."""
    d = _check_dimension(d)
    beta = _check_beta(d, beta)
    t = F(t)
    r = F(r)
    if r <= 0:
        raise ValueError(f"Eq. (36) requires r > 0, got r = {r}")
    return F((d * (d + 3) - 4 * beta) * t, 4 * d * (beta - 1) * r)


def criterion_ratio(d, beta, kappa_value) -> Fraction:
    """Eq. (37): E_U[o^2]/E_O[o^2] = (rho_L kappa + rho_S)/(kappa + 1)."""
    d = _check_dimension(d)
    kappa_value = F(kappa_value)
    return (rho_l(d, beta) * kappa_value + rho_s(d)) / (kappa_value + 1)


def variance_ratio_from(ratio, x) -> Fraction:
    """Eq. (44): Var_U/Var_O = (R - x)/(1 - x), x = m^2/E_O[o^2] in [0, 1)."""
    ratio = F(ratio)
    x = F(x)
    if not (F(0) <= x < F(1)):
        raise ValueError(f"Eq. (44) requires x in [0, 1), got x = {x}")
    return (ratio - x) / (F(1) - x)


# --------------------------------------------------------------- grid build

def _p_values(model: str, d: int):
    """Frozen p ladder for one (model, d), deduplicated preserving the frozen
    ladder order, with the Eq. (39) crossing point p*(d) appended for the
    depolarizing model when it is not already a ladder value."""
    if model == "depolarizing":
        values = list(DEPOLARIZING_P_LADDER)
        p_star = depolarizing_p_star(d)
        if p_star not in values:
            values.append(p_star)
        return values, p_star
    if model == "amplitude_damping":
        return list(AMPLITUDE_DAMPING_P_LADDER), None
    raise ValueError(f"unknown noise model {model!r}; expected one of {NOISE_MODELS}")


def _beta_for(model: str, d: int, p) -> Fraction:
    if model == "depolarizing":
        return depolarizing_beta(d, p)
    if model == "amplitude_damping":
        return amplitude_damping_beta(d, p)
    raise ValueError(f"unknown noise model {model!r}")


def build_rows() -> list:
    """The frozen 132-cell grid; every value an exact Fraction.

    Row order (deterministic): noise models in NOISE_MODELS order, then d
    ascending, then p in frozen-ladder order (crossing point appended when
    new), then anchors in ANCHORS order.
    """
    rows = []
    for model in NOISE_MODELS:
        for d in DIMENSIONS:
            p_values, p_star = _p_values(model, d)
            for p in p_values:
                roles = ["grid"]
                if p_star is not None and p == p_star:
                    roles.append("crossing_point")
                beta = _beta_for(model, d, p)
                for anchor in ANCHORS:
                    t, r, m = anchor_invariants(d, anchor)
                    e_o = second_moment_orthogonal(d, beta, t, r)
                    e_u = second_moment_unitary(d, beta, t, r)
                    var_o = e_o - m * m          # Eq. (28)
                    var_u = e_u - m * m          # Eq. (29)
                    ratio = e_u / e_o            # direct second-moment ratio
                    variance_ratio = var_u / var_o
                    k = kappa(d, beta, t, r)     # Eq. (36)
                    crit = criterion_ratio(d, beta, k)  # Eq. (37)
                    x = m * m / e_o
                    eq44_ratio = variance_ratio_from(ratio, x)
                    in_scope = F(1) < beta <= d
                    criterion_identity = crit == ratio
                    variance_identity = eq44_ratio == variance_ratio
                    dominance = variance_ratio >= ratio
                    positive = var_o > 0 and var_u > 0
                    rows.append({
                        "noise_model": model,
                        "model_source": MODEL_SECTION[model],
                        "d": d,
                        "n_qubits": _n_qubits(d),
                        "p": p,
                        "p_roles": roles,
                        "anchor": anchor,
                        "beta": beta,
                        "rho_s": rho_s(d),
                        "rho_l": rho_l(d, beta),
                        "t": t,
                        "r": r,
                        "m": m,
                        "second_moment_orthogonal": e_o,
                        "second_moment_unitary": e_u,
                        "variance_orthogonal": var_o,
                        "variance_unitary": var_u,
                        "kappa": k,
                        "x_mean_fraction": x,
                        "second_moment_ratio": ratio,
                        "criterion_ratio": crit,
                        "variance_ratio": variance_ratio,
                        "eq44_variance_ratio": eq44_ratio,
                        "criterion_identity_holds": criterion_identity,
                        "eq44_identity_holds": variance_identity,
                        "variance_ratio_geq_second_moment_ratio": dominance,
                        "beta_in_scope": in_scope,
                        "valid": bool(criterion_identity and variance_identity
                                      and dominance and positive and in_scope
                                      and x < 1),
                    })
    return rows


def grid_signature(rows: list) -> str:
    """Digest the exact frozen cell keys and p-role assignments."""
    ordered = sorted(
        rows,
        key=lambda row: (
            row["noise_model"], row["d"], row["p"], row["anchor"]
        ),
    )
    text = "".join(
        f"{row['noise_model']}|{row['d']}|{frac_text(row['p'])}|"
        f"{row['anchor']}|{','.join(row['p_roles'])}\n"
        for row in ordered
    )
    return sha256_bytes(text.encode("utf-8"))


def grid_summary(rows: list) -> dict:
    """Shape facts over the frozen grid (used by the regression gate)."""
    dep = sum(1 for row in rows if row["noise_model"] == "depolarizing")
    ad = sum(1 for row in rows if row["noise_model"] == "amplitude_damping")
    keys = {(row["noise_model"], row["d"], row["p"], row["anchor"])
            for row in rows}
    crossing = {
        (row["d"], row["p"], row["anchor"])
        for row in rows
        if row["noise_model"] == "depolarizing"
        and "crossing_point" in row["p_roles"]
    }
    expected_crossing = {
        (d, depolarizing_p_star(d), anchor)
        for d in DIMENSIONS
        for anchor in ANCHORS
    }
    signature = grid_signature(rows)
    return {
        "row_count": len(rows),
        "expected_row_count": 132,
        "model_counts": {"depolarizing": dep, "amplitude_damping": ad},
        "expected_model_counts": {"depolarizing": 60, "amplitude_damping": 72},
        "unique_cells": len(keys),
        "grid_signature_sha256": signature,
        "expected_grid_signature_sha256": EXPECTED_GRID_SHA256,
        "count_ok": len(rows) == 132 == len(keys),
        "split_ok": dep == 60 and ad == 72,
        "exact_grid_ok": signature == EXPECTED_GRID_SHA256,
        "crossing_roles_ok": crossing == expected_crossing,
    }


# --------------------------------------------------- regression identities

def regression_identities(rows: list) -> list:
    """Paper-verbatim and test-pinned exact identities.

    These guard the transcription: if any fails, the instrument is broken and
    the preregistered verdict is FROZEN-INCONCLUSIVE regardless of ratios."""
    checks = []

    def add(name, expected, actual):
        def show(v):
            if isinstance(v, Fraction):
                return frac_text(v)
            if isinstance(v, tuple):
                return [show(x) for x in v]
            return v
        checks.append({"name": name, "expected": show(expected),
                       "actual": show(actual), "ok": expected == actual})

    # Eq. (37)/(41) rho_S pinned values (also test_noisy_variance.py).
    for d, expected in zip(DIMENSIONS, (F(10, 9), F(27, 25), F(85, 81))):
        add(f"rho_s(d={d})", expected, rho_s(d))
    # Eq. (39) beta* pinned values.
    for d, expected in zip(DIMENSIONS, (F(13, 4), F(128, 23), F(752, 77))):
        add(f"beta_star(d={d})", expected, beta_star(d))
    # Depolarizing crossing points; p*(4) = 3/4 and p*(8) = 15/23 are verbatim
    # in the paper's Section 3.3 numerics.
    for d, expected in zip(DIMENSIONS, (F(3, 4), F(15, 23), F(45, 77))):
        add(f"depolarizing_p_star(d={d})", expected, depolarizing_p_star(d))
    # Paper's own numeric cross-check: "At d=4 and p=0.5, rho_L = 1.852".
    add("rho_l(d=4, beta=5/2)", F(50, 27),
        rho_l(4, depolarizing_beta(4, F(1, 2))))
    add("depolarizing_beta(d=4, p=1/2)", F(5, 2), depolarizing_beta(4, F(1, 2)))
    # Eq. (43): rho_L - 2 vanishes exactly at beta*(d).
    for d in DIMENSIONS:
        add(f"rho_l(d={d}, beta=beta_star(d))", F(2), rho_l(d, beta_star(d)))
    # Anchor invariants at d=4 (test-pinned tuples; Eq. (31) definitions).
    add("anchor_invariants(4, maximally_mixed_pauli_z)",
        (F(4), F(1), F(0)), anchor_invariants(4, "maximally_mixed_pauli_z"))
    add("anchor_invariants(4, ghz_projector)",
        (F(3, 4), F(9, 16), F(3, 4)), anchor_invariants(4, "ghz_projector"))
    add("anchor_invariants(4, ghz_pauli_x)",
        (F(4), F(1), F(1)), anchor_invariants(4, "ghz_pauli_x"))
    # Amplitude-damping boundary: beta = d exactly at p = 1 (section 3.4.2).
    for d in DIMENSIONS:
        add(f"amplitude_damping_beta(d={d}, p=1)", F(d),
            amplitude_damping_beta(d, F(1)))
    # Grid shape.
    grid = grid_summary(rows)
    add("grid_row_count", grid["expected_row_count"], grid["row_count"])
    add("grid_unique_cells", grid["expected_row_count"], grid["unique_cells"])
    for model in NOISE_MODELS:
        add(f"grid_count[{model}]", grid["expected_model_counts"][model],
            grid["model_counts"][model])
    add("exact_frozen_grid_signature", EXPECTED_GRID_SHA256,
        grid["grid_signature_sha256"])
    add("exact_crossing_roles", True, grid["crossing_roles_ok"])
    # Per-row paper identities: Eq. (37) equals the direct ratio everywhere
    # (Prop. 3.10) and Cor. 3.11 dominance holds everywhere.
    add("eq37_equals_direct_all_rows", True,
        all(row["criterion_identity_holds"] for row in rows))
    add("eq44_equals_direct_variance_ratio_all_rows", True,
        all(row["eq44_identity_holds"] for row in rows))
    add("variance_ratio_geq_second_moment_ratio_all_rows", True,
        all(row["variance_ratio_geq_second_moment_ratio"] for row in rows))
    add("beta_in_scope_all_rows", True, all(row["beta_in_scope"] for row in rows))
    add("variances_positive_all_rows", True,
        all(row["variance_orthogonal"] > 0 and row["variance_unitary"] > 0
            for row in rows))
    return checks


# ------------------------------------------------------------------ verdict

def classify_verdict(identities_valid: bool, grid_valid: bool,
                     rows_valid: bool, any_leq_one: bool) -> tuple[str, str]:
    """Map instrument state and observations to control-plane verdicts."""
    if not (identities_valid and grid_valid and rows_valid):
        return "INCONCLUSIVE-IDENTITIES-INVALID", "FROZEN-INCONCLUSIVE"
    if any_leq_one:
        return "FLIP-FOUND", "FROZEN-NEGATIVE"
    return "NO-FLIP-CONFIRMED-IN-SCOPE", "FROZEN-CERTIFIED"


def evaluate(rows: list, identities: list) -> dict:
    """Frozen verdict mapping from preregistration Revisions 3-7."""
    identities_valid = all(check["ok"] for check in identities)
    grid = grid_summary(rows)
    grid_valid = (grid["count_ok"] and grid["split_ok"]
                  and grid["exact_grid_ok"] and grid["crossing_roles_ok"])
    rows_valid = all(row["valid"] for row in rows)
    second = [row["second_moment_ratio"] for row in rows]
    variances = [row["variance_ratio"] for row in rows]
    any_leq_one = any(v <= 1 for v in second + variances)

    def cell_of(row):
        return {"noise_model": row["noise_model"], "d": row["d"],
                "p": frac_text(row["p"]), "anchor": row["anchor"]}

    def extremum(values):
        best = min(values)
        worst = max(values)
        return {
            "min": frac_obj(best),
            "min_cell": cell_of(rows[values.index(best)]),
            "max": frac_obj(worst),
            "max_cell": cell_of(rows[values.index(worst)]),
        }

    scientific, terminal = classify_verdict(
        identities_valid, grid_valid, rows_valid, any_leq_one)

    flips = [
        {"noise_model": row["noise_model"], "d": row["d"],
         "p": frac_text(row["p"]), "anchor": row["anchor"],
         "second_moment_ratio": frac_text(row["second_moment_ratio"]),
         "variance_ratio": frac_text(row["variance_ratio"])}
        for row in rows
        if row["second_moment_ratio"] <= 1 or row["variance_ratio"] <= 1
    ]
    return {
        "regression_identities_valid": identities_valid,
        "all_rows_valid": rows_valid,
        "row_count": grid["row_count"],
        "model_counts": grid["model_counts"],
        "any_ratio_leq_one": any_leq_one,
        "flip_cells": flips,
        "second_moment_ratio_extremes": extremum(second),
        "variance_ratio_extremes": extremum(variances),
        "scientific_verdict": scientific,
        "terminal_verdict_by_prereg_mapping": terminal,
    }


def verdict_line(verdict: dict) -> str:
    return (f"verdict={verdict['terminal_verdict_by_prereg_mapping']} "
            f"scientific={verdict['scientific_verdict']} "
            f"rows={verdict['row_count']} "
            f"min_second_moment_ratio="
            f"{verdict['second_moment_ratio_extremes']['min']['text']} "
            f"min_variance_ratio="
            f"{verdict['variance_ratio_extremes']['min']['text']}")


# ------------------------------------------------------------- manifest/CLI

def load_manifest(run_dir: Path) -> dict:
    """Require an existing, initialized, RUNNING campaign manifest.

    Never creates anything; every failure raises RunnerRefused."""
    run_dir = Path(run_dir)
    if not run_dir.is_dir():
        raise RunnerRefused(
            f"REFUSED: run directory does not exist: {run_dir} "
            "(the runner never creates run directories)")
    manifest_path = run_dir / MANIFEST_NAME
    if not manifest_path.is_file():
        raise RunnerRefused(
            f"REFUSED: {manifest_path} missing -- campaign not initialized")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RunnerRefused(
            f"REFUSED: {manifest_path} is not readable JSON: {exc}") from exc
    if not isinstance(manifest, dict):
        raise RunnerRefused(f"REFUSED: {manifest_path} is not a JSON object")
    required = ("run_id", "target", "gate", "agent", "created_utc",
                "prereg_sha256", "status")
    missing = [key for key in required if not manifest.get(key)]
    if missing:
        raise RunnerRefused(
            f"REFUSED: {manifest_path} is not a campaign.py manifest; "
            f"missing fields: {missing}")
    if manifest["run_id"] != run_dir.name:
        raise RunnerRefused(
            f"REFUSED: manifest run_id {manifest['run_id']!r} does not match "
            f"directory {run_dir.name!r}")
    if manifest["target"] != EXPECTED_TARGET:
        raise RunnerRefused(
            f"REFUSED: manifest target {manifest['target']!r}, "
            f"expected {EXPECTED_TARGET!r}")
    if manifest["gate"] != EXPECTED_GATE:
        raise RunnerRefused(
            f"REFUSED: manifest gate {manifest['gate']!r}, "
            f"expected {EXPECTED_GATE!r}")
    status = manifest["status"]
    if status != MANIFEST_STATUS_RUNNING:
        raise RunnerRefused(
            f"REFUSED: manifest status is {status!r}, "
            f"expected {MANIFEST_STATUS_RUNNING!r}")
    current_prereg_sha = sha256_file(PREREG_PATH)
    if manifest["prereg_sha256"] != current_prereg_sha:
        raise RunnerRefused(
            "REFUSED: current preregistration SHA256 does not match the "
            "campaign.py manifest")
    return manifest


def _dumps(payload: dict) -> bytes:
    text = json.dumps(to_jsonable(payload), indent=1, sort_keys=True,
                      ensure_ascii=True)
    return (text + "\n").encode("utf-8")


def prereg_revision_from_bytes(data: bytes) -> int:
    """Return the latest strictly increasing Revision heading, else -1."""
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return -1
    revisions = [
        int(value)
        for value in re.findall(r"^## Revision ([1-9][0-9]*) \(", text, re.MULTILINE)
    ]
    if not revisions or revisions != sorted(set(revisions)):
        return -1
    return revisions[-1]


def run_payloads(run_dir: Path):
    """Bind run inputs, compute the table, and assemble artifact payloads."""
    run_dir = Path(run_dir)
    manifest = load_manifest(run_dir)
    _assert_module_stable()
    module_bytes = MODULE_BYTES_AT_IMPORT
    prereg_bytes = PREREG_PATH.read_bytes()
    module_sha = sha256_bytes(module_bytes)
    prereg_sha = sha256_bytes(prereg_bytes)
    if manifest["prereg_sha256"] != prereg_sha:
        raise RunnerRefused(
            "REFUSED: captured preregistration does not match the campaign "
            "manifest SHA256")
    prereg_revision = prereg_revision_from_bytes(prereg_bytes)
    if prereg_revision != PREREG_REVISION:
        raise RunnerRefused(
            f"REFUSED: latest preregistration revision is {prereg_revision}; "
            f"expected {PREREG_REVISION}")

    rows = build_rows()
    identities = regression_identities(rows)
    verdict = evaluate(rows, identities)
    identities_valid = verdict["regression_identities_valid"]
    _assert_module_stable()
    if PREREG_PATH.read_bytes() != prereg_bytes:
        raise RunnerRefused(
            "REFUSED: preregistration changed during evaluation")
    prereg = {
        "path": str(PREREG_PATH.relative_to(ROOT)),
        "revision": prereg_revision,
        "sha256": prereg_sha,
    }
    snapshots = {
        "module": {"path": str(MODULE_PATH.relative_to(ROOT)),
                   "snapshot": MODULE_SNAPSHOT_NAME,
                   "sha256": module_sha, "bytes": len(module_bytes)},
        "prereg": {"path": prereg["path"],
                   "snapshot": PREREG_SNAPSHOT_NAME,
                   "sha256": prereg_sha, "bytes": len(prereg_bytes)},
    }

    results = {
        "schema": RESULTS_SCHEMA,
        "campaign": manifest["run_id"],
        "campaign_manifest": manifest,
        "paper": PAPER,
        "preregistration": prereg,
        "module_sha256": module_sha,
        "snapshots": snapshots,
        "arithmetic": {
            "verdict_path": "fractions.Fraction only",
            "rng": "none",
            "floats": "none",
            "matrix_allocation": "none",
        },
        "grid": grid_summary(rows),
        "rows": rows,
        "scientific_verdict": verdict["scientific_verdict"],
        "terminal_verdict_by_prereg_mapping":
            verdict["terminal_verdict_by_prereg_mapping"],
    }
    results_bytes = _dumps(results)

    summary = {
        "schema": SUMMARY_SCHEMA,
        "campaign": manifest["run_id"],
        "campaign_manifest_status": manifest["status"],
        "paper": {"arxiv_id": PAPER["arxiv_id"], "version": PAPER["version"],
                  "equations": PAPER["equations"],
                  "supporting": PAPER["supporting"]},
        "preregistration": prereg,
        "module_sha256": module_sha,
        "snapshots": snapshots,
        "grid": grid_summary(rows),
        "regression_identities": identities,
        "verdict": verdict,
        "results_sha256": sha256_bytes(results_bytes),
        "results_bytes": len(results_bytes),
        "scientific_verdict": verdict["scientific_verdict"],
        "terminal_verdict_by_prereg_mapping":
            verdict["terminal_verdict_by_prereg_mapping"],
    }
    summary_bytes = _dumps(summary)

    inventory = {
        "schema": INVENTORY_SCHEMA,
        "campaign": manifest["run_id"],
        "manifest": {"path": MANIFEST_NAME, "run_id": manifest["run_id"],
                     "target": manifest["target"], "gate": manifest["gate"],
                     "status": manifest["status"],
                     "sha256": sha256_file(run_dir / MANIFEST_NAME)},
        "module": {"path": snapshots["module"]["path"],
                   "snapshot": MODULE_SNAPSHOT_NAME, "sha256": module_sha},
        "preregistration": {"path": prereg["path"], "revision": prereg["revision"],
                            "snapshot": PREREG_SNAPSHOT_NAME,
                            "sha256": prereg_sha},
        "artifacts": [
            {"name": MANIFEST_NAME,
             "sha256": sha256_file(run_dir / MANIFEST_NAME),
             "bytes": (run_dir / MANIFEST_NAME).stat().st_size},
            {"name": MODULE_SNAPSHOT_NAME, "origin": snapshots["module"]["path"],
             "sha256": module_sha, "bytes": len(module_bytes)},
            {"name": PREREG_SNAPSHOT_NAME, "origin": snapshots["prereg"]["path"],
             "sha256": prereg_sha, "bytes": len(prereg_bytes)},
            {"name": "results.json", "sha256": sha256_bytes(results_bytes),
             "bytes": len(results_bytes)},
            {"name": "summary.json", "sha256": sha256_bytes(summary_bytes),
             "bytes": len(summary_bytes)},
        ],
        "row_count": verdict["row_count"],
        "model_counts": verdict["model_counts"],
        "scientific_verdict": verdict["scientific_verdict"],
        "terminal_verdict_by_prereg_mapping":
            verdict["terminal_verdict_by_prereg_mapping"],
        "summary_sha256": sha256_bytes(summary_bytes),
    }
    inventory_bytes = _dumps(inventory)
    return manifest, identities_valid, verdict, [
        (MODULE_SNAPSHOT_NAME, module_bytes),
        (PREREG_SNAPSHOT_NAME, prereg_bytes),
        ("results.json", results_bytes),
        ("summary.json", summary_bytes),
        ("inventory.json", inventory_bytes),
    ]


def write_artifacts(run_dir: Path) -> int:
    manifest, identities_valid, verdict, artifacts = run_payloads(run_dir)
    run_dir = Path(run_dir)
    _assert_module_stable()
    for name, _ in artifacts:
        target = run_dir / name
        if target.exists():
            raise RunnerRefused(
                f"REFUSED: {target} already exists -- this campaign directory "
                "already carries artifacts; refusing to overwrite")
    for name, payload in artifacts:
        (run_dir / name).write_bytes(payload)
    print(f"RUN {run_dir}: wrote " + ", ".join(name for name, _ in artifacts))
    print(verdict_line(verdict))
    return 0 if identities_valid else 1


def check() -> int:
    rows = build_rows()
    identities = regression_identities(rows)
    verdict = evaluate(rows, identities)
    for identity in identities:
        if not identity["ok"]:
            print(f"IDENTITY FAIL: {identity['name']}: expected "
                  f"{identity['expected']} got {identity['actual']}")
    print(f"CHECK: identities={len(identities)} "
          f"failed={sum(1 for c in identities if not c['ok'])}")
    print(verdict_line(verdict))
    return 0 if verdict["regression_identities_valid"] else 1


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Exact 132-cell real-vs-unitary noisy-shadow variance "
                    "confirmation runner (arXiv:2608.18935v1).")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--run-dir", type=Path, default=None,
                       help="existing initialized RUNNING campaign directory "
                            "to write the module/prereg snapshots plus "
                            "results.json/summary.json/inventory.json into "
                            "(never created by this runner)")
    group.add_argument("--check", action="store_true",
                       help="compute the grid and identities; write nothing")
    args = parser.parse_args(argv)
    try:
        if args.check:
            return check()
        return write_artifacts(args.run_dir)
    except RunnerRefused as exc:
        print(str(exc), file=sys.stderr)
        return 2
    except (ValueError, ZeroDivisionError) as exc:
        print(f"REFUSED: invalid input or arithmetic domain: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
