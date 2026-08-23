#!/usr/bin/env python3
"""Integrate the all-auxiliary Lax method-limitation certificate."""

from __future__ import annotations

import hashlib
import json
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from e212_allaux_intertwiner import run_intertwiner_audit
from e213_allaux_controls import run_control_audit

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "integrability" / "allaux_lax.json"
CPU_BUDGET_SECONDS = 45.0
RSS_CAP_BYTES = 2_000_000_000


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def budget_tick(started: float, stage: str) -> None:
    used = time.process_time() - started
    if used > CPU_BUDGET_SECONDS:
        raise RuntimeError(f"process-time budget exceeded at {stage}: {used}")
    rss = max_rss_bytes()
    if rss >= RSS_CAP_BYTES:
        raise MemoryError(f"RSS cap exceeded at {stage}: {rss}")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scope_guard() -> dict[str, object]:
    paths = {
        "tetra_ungraded": ROOT / "proofs" / "tetra_ungraded.md",
        "tetra16_locus": ROOT / "proofs" / "tetra16_locus.md",
        "local_term_trichotomy": ROOT
        / "proofs"
        / "clifford_grade_classification.md",
    }
    texts = {name: path.read_text() for name, path in paths.items()}
    lower = {name: text.lower() for name, text in texts.items()}
    checks = {
        "ungraded_theorem_explicitly_excludes_higher_auxiliary_dimension": (
            "does not classify higher auxiliary dimensions" in lower["tetra_ungraded"]
        ),
        "cabled_theorem_explicitly_excludes_unrestricted_C4_leg": (
            "does not establish an rlll solution for an unrestricted"
            in lower["tetra16_locus"]
            and "fixed-order one-line-cabled" in lower["tetra16_locus"]
        ),
        "trichotomy_is_a_physical_local_term_lie_algebra_statement": (
            "g_gamma = lie_r" in lower["local_term_trichotomy"]
            and "act on the full `2^n`-dimensional spin hilbert space"
            in lower["local_term_trichotomy"]
        ),
    }
    return {
        "tag": "[THEOREM]",
        "source_sha256": {name: sha256_file(path) for name, path in paths.items()},
        "computed_scope_checks": checks,
        "all_scope_checks_passed": all(checks.values()),
        "preserved_8x8_statement": (
            "[THEOREM] The existing arbitrary-8x8 identical-L Ising RLLL kernel is zero "
            "for characteristic-zero q outside {0,+1,-1}; this certificate neither "
            "reproves nor enlarges that theorem."
        ),
        "preserved_16x16_statement": (
            "[THEOREM] The existing 16x16 result is only the named fixed-order one-line "
            "cabling; it is not an unrestricted C4-leg theorem."
        ),
    }


def formal_problem() -> dict[str, object]:
    return {
        "tag": "[LEMMA]",
        "field_assumption": "characteristic-zero field k",
        "spaces": (
            "A=A1 tensor A2 tensor A3 is the auxiliary product and "
            "Q=Q4 tensor Q5 tensor Q6 is the physical product"
        ),
        "ordered_products": {
            "F": "L_145 L_246 L_356",
            "G": "L_356 L_246 L_145",
        },
        "rlll_equation": "(R tensor I_Q) F = G (R tensor I_Q)",
        "coefficient_map": (
            "Phi_L: End_k(A) -> End_k(A tensor Q), "
            "R |-> (R tensor I_Q)F-G(R tensor I_Q)"
        ),
        "solution_space": "ker(Phi_L)",
        "one_operator_intertwiner_identity": (
            "ker(Phi_L)=Hom_{k[t]}((A tensor Q,F),(A tensor Q,G)) "
            "intersection (End_k(A) tensor I_Q)"
        ),
        "missing_implication": (
            "[UNRESOLVED] No identity deriving T rho(x)=sigma(x) T for every physical "
            "local-term-algebra element x follows from the single displayed RLLL equation."
        ),
        "dimension_change_requirement": (
            "[LEMMA] Changing dim(Ai) changes the domains of the local L factors and "
            "requires a newly specified L-family; the fixed binary-leg 8x8 coefficient "
            "matrix has no canonical arbitrary-d continuation."
        ),
    }


def schur_theorem() -> dict[str, object]:
    return {
        "tag": "[THEOREM]",
        "named_assumptions": [
            "B is a unital associative k-algebra",
            "rho:B->End(P) and sigma:B->End(Q) are finite-dimensional simple modules",
            "for the one-dimensional conclusion End_B(P)=k (absolute Schur condition)",
            "the candidate T satisfies T rho(b)=sigma(b) T for every b in B",
        ],
        "zero_criterion": (
            "If P and Q are nonisomorphic simple B-modules, Hom_B(P,Q)=0."
        ),
        "equivalent_case": (
            "If P and Q are isomorphic and End_B(P)=k, then Hom_B(P,Q)=kS for "
            "any fixed module isomorphism S; every nonzero solution is invertible."
        ),
        "tensor_multiplicity": (
            "Hom_B(U tensor P,V tensor Q)=Hom_k(U,V) tensor Hom_B(P,Q) when B "
            "acts trivially on U,V."
        ),
        "commutant_corollary": (
            "If End_B(P)=k, then the commutant of I_U tensor rho(B) is "
            "End_k(U) tensor I_P and has dimension dim(U)^2."
        ),
        "double_commutant_warning": (
            "[LEMMA] Over an algebraically closed field Burnside/double-commutant "
            "arguments may identify rho(B)=End(P), but this yields scalar physical "
            "commutants or auxiliary-square multiplicities, not absence of R."
        ),
        "correct_application_gate": (
            "[THEOREM] Schur gives an RLLL no-go only after the RLLL component equations "
            "are proved to impose an all-elements intertwiner between named nonisomorphic "
            "simple representations. Algebra dimension or irreducibility alone is insufficient."
        ),
    }


def counterfamily_theorem() -> dict[str, object]:
    return {
        "tag": "[THEOREM]",
        "scope": "every integer auxiliary/quantum dimension d>=2 and every chain length n>=1 over Q",
        "strict_permutation_rlll": {
            "tag": "[THEOREM]",
            "spaces": "six copies of V=Q^d in the tetrahedral 123,145,246,356 placement",
            "local_operator": "L=P_(first,second) tensor I_third",
            "auxiliary_operator": "R=P_12 tensor I_3",
            "identity": "R_123 L_145 L_246 L_356=L_356 L_246 L_145 R_123",
            "proof_identity": (
                "P_12 P_14 P_24=P_24 P_14 P_12 and P_35 commutes with all "
                "three factors; this is an equality in S6 and hence in every tensor "
                "permutation representation."
            ),
            "invertibility": "R is a nonscalar permutation with determinant plus or minus 1",
            "local_entry_envelope": (
                "Relative to the first leg, the entries of L are E_sr on the active "
                "second leg tensor I on the spectator third leg, so they contain full End(V)."
            ),
        },
        "rational_yang_rll": {
            "tag": "[THEOREM]",
            "spaces": "V=Q^d; auxiliary copies Va,Vb and quantum copy Vq",
            "operators": {
                "L_aq(u)": "u I + P_aq",
                "L_bq(v)": "v I + P_bq",
                "R_ab(u-v)": "(u-v) I + P_ab",
            },
            "identity": (
                "R_ab(u-v)L_aq(u)L_bq(v)=L_bq(v)L_aq(u)R_ab(u-v) in Z[u,v][S3]"
            ),
            "invertible_sample": (
                "At u=2,v=-1, det R=4^{d(d+1)/2} 2^{d(d-1)/2}, so R is "
                "nonscalar and invertible for every d>=2."
            ),
        },
        "local_entry_formula": "L_rs(u)=u delta_rs I+E_sr",
        "entry_envelope": (
            "The off-diagonal entries and their products contain every E_rs; on n sites "
            "their associative envelope is End(V^{tensor n}), dimension d^{2n}, irreducible."
        ),
        "entry_derived_full_lie_control": (
            "Take onsite sl_d matrix units/diagonal differences from L entries and one "
            "nearest-neighbour H_j H_{j+1}, H=E_11-E_22, on every edge. Their Lie "
            "algebra is sl_{d^n}, dimension d^{2n}-1, by induction using simplicity of "
            "sl_m and sl_m tensor I + I tensor sl_d + sl_m tensor sl_d=sl_{md}."
        ),
        "method_limitation": (
            "An exponentially large irreducible local algebra, even one inside the local "
            "L-entry product envelope, does not by itself forbid a local R/L or strict "
            "RLLL intertwiner at any auxiliary dimension. A compatibility theorem tying "
            "the particular physical Ising generators to all RLLL component equations "
            "is indispensable."
        ),
        "physical_model_caveat": (
            "[LEMMA] The full-control set certifies the logical insufficiency of algebra "
            "size; it does not assert that adding arbitrary controls to the XXX Hamiltonian "
            "preserves its commuting transfer family."
        ),
    }


def known_model_audit(controls: dict[str, object]) -> dict[str, object]:
    return {
        "tag": "[EXTERNAL]",
        "identification": (
            "The rational permutation solution is the standard GL_d/XXX Yang R-matrix."
        ),
        "exact_content_is_internal": (
            "[THEOREM] The RLL identity, invertible sample, entry-envelope result, and "
            "all finite controls in this artifact are exact internal derivations and do "
            "not rely on the external name."
        ),
        "decomposition_audit": (
            "[COMPUTATION] For the integrable spin-1/2 XXX density, splitting the Pauli "
            "part into individual XX,YY,ZZ bond summands gives the stored exact closures "
            "at n=3,4,5,6. These are finite controls only, not an all-n theorem."
        ),
        "finite_rows": controls["xxx_pauli_summand_closures"],
        "interpretation": (
            "[LEMMA] The size of a local-term algebra depends on which summands are treated "
            "as independent generators, while the RLL identity is an identity of L and R."
        ),
    }


def make_artifact(
    intertwiner: dict[str, object],
    controls: dict[str, object],
    checks: list[dict[str, object]],
    started: float,
) -> dict[str, object]:
    sources = (
        "experiments/e212_allaux_intertwiner.py",
        "experiments/e213_allaux_controls.py",
        "experiments/e214_allaux_certificate.py",
    )
    return {
        "meta": {
            "experiment": "e214_allaux_certificate",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "command": ".venv/bin/python experiments/e214_allaux_certificate.py",
            "working_directory": str(ROOT),
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "arithmetic": (
                "exact Z[u,v][S3] group algebra, Python integers/Fraction, exact Q "
                "row reduction, and GF(2) symplectic Pauli closure"
            ),
            "budgets": {
                "process_time_seconds": CPU_BUDGET_SECONDS,
                "rss_cap_bytes": RSS_CAP_BYTES,
            },
            "observed_process_time_seconds": str(time.process_time() - started),
            "observed_peak_rss_bytes": max_rss_bytes(),
            "source_sha256": {
                source: sha256_file(ROOT / source) for source in sources
            },
        },
        "data": {
            "verdict": {
                "tag": "[THEOREM]",
                "headline": (
                    "REFUTED METHOD: exponential dimension and irreducibility of a physical "
                    "local-term algebra do not imply absence of local Lax/R-type intertwiners "
                    "at any auxiliary dimension."
                ),
                "positive_theorem": (
                    "Schur yields zero only for an all-elements intertwiner between named "
                    "nonisomorphic simple modules; tensor-factor irreducibility instead leaves "
                    "an auxiliary-square commutant."
                ),
                "ising_conclusion": (
                    "[UNRESOLVED] Arbitrary higher-dimensional Ising L factors remain "
                    "unclassified. The exact 8x8 and fixed cabled 16x16 theorems retain "
                    "exactly their published scope."
                ),
            },
            "formal_rlll_intertwiner": formal_problem(),
            "schur_and_tensor_commutant_theorem": schur_theorem(),
            "all_dimension_counterfamily": counterfamily_theorem(),
            "known_integrable_model_audit": known_model_audit(controls),
            "exact_intertwiner_controls": intertwiner,
            "exact_algebra_controls": controls,
            "existing_theorem_scope_guard": scope_guard(),
        },
        "checks": checks,
    }


def run_all() -> tuple[dict[str, object], list[dict[str, object]]]:
    started = time.process_time()
    intertwiner, intertwiner_checks = run_intertwiner_audit()
    budget_tick(started, "intertwiner audit")
    controls, control_checks = run_control_audit()
    budget_tick(started, "algebra controls")
    guard = scope_guard()

    integration_checks = [
        {
            "name": "all-d group-algebra identity and every finite RLL control agree",
            "passed": bool(
                intertwiner["universal_yang_rll"]["identity_passed"]  # type: ignore[index]
            )
            and all(
                bool(row["passed"])
                for row in intertwiner["finite_sparse_controls"]  # type: ignore[index]
            ),
            "detail": "universal Z[u,v][S3] residual=0; exact d=2,3,4 residuals=0",
        },
        {
            "name": "strict all-d permutation RLLL identity and finite controls agree",
            "passed": bool(
                intertwiner["universal_permutation_rlll"]["identity_passed"]  # type: ignore[index]
            )
            and all(
                bool(row["passed"])
                for row in intertwiner["finite_permutation_rlll_controls"]  # type: ignore[index]
            ),
            "detail": "exact S6 identity; nonscalar invertible R at d=2,3,4",
        },
        {
            "name": "large irreducible entry envelopes coexist with the certified RLL family",
            "passed": all(
                bool(row["passed"])
                for row in controls["local_l_entry_envelopes"]  # type: ignore[index]
            )
            and all(
                bool(row["passed"])
                for row in controls["tensor_product_entry_envelopes"]  # type: ignore[index]
            ),
            "detail": "full End bases at d=2,3,4 locally and (d,n)=(2,3),(3,2)",
        },
        {
            "name": "Schur hypotheses, not algebra size, control exact intertwiner nullity",
            "passed": bool(controls["schur_controls"]["all_passed"])  # type: ignore[index]
            and all(
                bool(row["passed"])
                for row in controls["tensor_factor_commutants"]  # type: ignore[index]
            ),
            "detail": "nonisomorphic nullity 0; equivalent nullity 1; tensor nullity a^2",
        },
        {
            "name": "existing 8x8 and cabled 16x16 theorem boundaries are machine-guarded",
            "passed": bool(guard["all_scope_checks_passed"]),
            "detail": "; ".join(
                name
                for name, passed in guard["computed_scope_checks"].items()  # type: ignore[union-attr]
                if passed
            ),
        },
    ]
    checks = [*intertwiner_checks, *control_checks, *integration_checks]
    if not all(bool(check["passed"]) for check in checks):
        raise AssertionError(checks)
    artifact = make_artifact(intertwiner, controls, checks, started)
    budget_tick(started, "artifact assembly")
    return artifact, checks


def main() -> int:
    artifact, checks = run_all()
    ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(
        f"PASS e214 ({len(checks)} checks, artifact={ARTIFACT}, "
        f"peak_rss={max_rss_bytes()})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
