# -*- coding: utf-8 -*-
"""Formula translation module -- pure arithmetic translation of every closed
formula in arXiv:2507.12024v2 plus the Litinski-2019b pieces its Table 6 uses.

Float policy
------------
Distances, qubit counts, toffoli counts and the integer factors of Sec. 3.5
are Python ints or Fraction -- exact, no rounding.  Physical times and
probabilities enter as the DECIMAL literals printed by the source (they are
inherently measured quantities); arithmetic on them uses the platform double
(IEEE-754 binary64), which carries ~15-16 significant digits -- far more than
the 4-6 digits the paper prints, so faithful reproduction of printed values
dominates float error by >= 6 orders of magnitude.  `ceil` and integer
rounding are applied exactly as the paper's formulas state, on exact ints
wherever the input allows (Fraction), else on the float, and any observed
boundary sensitivity is reported rather than hidden.  No tolerance or epsilon
is applied inside a formula; tolerances belong to the verification layer.
"""

from fractions import Fraction
import math

# --- Eq. (1) (S2.E1) -------------------------------------------------------
def p0_from_pL(pL, k, d):
    """p0 = 1-(1-pL)^(1/(k*d)); paper's Eq.(1). Uses exact Fraction power
    unless pL is float, in which case float math (see float policy)."""
    if isinstance(pL, Fraction):
        return 1 - (1 - pL) ** Fraction(1, int(k) * int(d))
    return 1.0 - (1.0 - float(pL)) ** (1.0 / (int(k) * int(d)))


# --- Eq. (2) (S2.E2) -------------------------------------------------------
def qlops_eq2(k, d, t_sec, t_r):
    """Q = k / ((ceil(tr/t_sec)+d) * t_sec), per Eq.(2).
    k and d ints -> ceil on exact Fraction ratio when possible."""
    if isinstance(t_r, Fraction) and isinstance(t_sec, Fraction):
        rounds = -((-(t_r / t_sec)).__floor__())  # ceil via Fraction
    else:
        rounds = math.ceil(float(t_r) / float(t_sec))
    return int(k) / ((rounds + int(d)) * float(t_sec))


def qlops_ignoring_latency(k, t_sec):
    """Sec. 3.5 neutral-atom variant: Q = k * 1/t_SEC (reaction ignored)."""
    return int(k) / float(t_sec)


def qlops_density(qlops, n_physical_qubits):
    return qlops / n_physical_qubits


# --- Sec. 3.5 arithmetic ---------------------------------------------------
def rsa_sc_logical_qubits(n_cold, n_active):
    return int(n_cold) + int(n_active)


def rsa_sc_dataloader_qubits(n_cold, qubits_per_cold_qubit=430):
    """1280 x 430 = 550,400 (the paper's per-cold-qubit footprint, S3.SS5.p2)."""
    return int(n_cold) * int(qubits_per_cold_qubit)


def rsa_sc_active_qubits(n_active, d):
    """131 * (2d^2 - 1) for the active qubits, d=25 -> 162,499 (S3.SS5.p2)."""
    d = int(d)
    return int(n_active) * (2 * d * d - 1)


def surface_code_patch_physical_qubits(d, with_ancilla=True):
    """d^2 data + (d^2 - 1) ancilla = 2d^2-1; with_ancilla=False -> d^2.
    Used by both the active-block formula and the Zhou et al. patches."""
    d = int(d)
    return d * d + (d * d - 1) if with_ancilla else d * d


def rsa_ratio_eq3(q_sc, t_sc_days, q_na, t_na_days):
    """Q_SC*t_SC / (Q_atom*t_atom) with t in days -- days cancel; Eq.(3)."""
    return (q_sc * t_sc_days) / (q_na * t_na_days)


def rsa_ratio_eq4(q_sc, t_sc_days, n_sc, q_na, t_na_days, n_na):
    """(Q_SC*t_SC/n_SC_T) / (Q_atom*t_atom/n_NA_T), Eq.(4)."""
    return ((q_sc * t_sc_days) / n_sc) / ((q_na * t_na_days) / n_na)


def lower_bound_time_days(n_toffoli, qlops, non_clifford_per_toffoli=10,
                          seconds_per_day=86400.0):
    """10*nT/Q in days (S3.SS5.p6: 'a Toffoli requires ~10 Clifford gates')."""
    return (non_clifford_per_toffoli * n_toffoli / qlops) / seconds_per_day


# --- Litinski 15-to-1 (basis of Table 6) ----------------------------------
def litinski_15to1_unit_qubits(dX, dZ, dm):
    """Litinski 2019b S3.p4 (Fig. 11): 2*(dX+4dZ)*3dX + 4dm physical qubits."""
    dX, dZ, dm = int(dX), int(dZ), int(dm)
    return 2 * (dX + 4 * dZ) * 3 * dX + 4 * dm


def litinski_15to1_plain_cycles(dm):
    """Litinski 2019b S3.p4: time cost 6*dm code cycles, pre post-selection."""
    return 6 * int(dm)


def litinski_pfail_from_cycles(cycles, dm):
    """Invert cycles = 6*dm/(1-p_fail) for the paper's printed cycles."""
    dm = int(dm)
    return 1.0 - (6.0 * dm) / float(cycles)


def distillation_unit_count(code_cycle_s_computation, code_cycle_s_distillation,
                            cycles_per_output, n_logical):
    """The paper's Table-6 'Total' convention, S3.T6 caption: 'the number of
    qubits in each distillation unit multiplied by the number of distillation
    units needed... takes into account the syndrome extraction cycle
    difference and the post-selection rate, and is rounded up to an integer.'
    One logical op needs 1 magic state per logical qubit per logical cycle;
    per SEC of the computation cycle clock, a unit produces
    (t_computation/t_distillation) / cycles_per_output states."""
    rate = (float(code_cycle_s_computation) / float(code_cycle_s_distillation)) \
        / float(cycles_per_output)
    return math.ceil(int(n_logical) / rate)


# --- surface-code logical-error model (Litinski 2019b Eq. 7) --------------
def litinski_pL(pphys, d):
    """Litinski 2019b Eq.(7): pL = 0.1*(100*p_phys)^((d+1)/2).
    NOT used for the reproduction (the paper replaces this with its own
    Fig. 1 fit); provided for gate-B cross-checks only."""
    return 0.1 * ((100.0 * float(pphys)) ** ((int(d) + 1) / 2.0))
