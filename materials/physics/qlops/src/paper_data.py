# -*- coding: utf-8 -*-
"""Transcription registry: evry formula & constant of arXiv:2507.12024 (v2).

Every constant in this file was read first-hand from the arXiv HTML v2 full
text (https://arxiv.org/html/2507.12024v2, fetched 2026-08-29) of

    "Benchmarking fault-tolerant quantum computing hardware via QLOPS",
    Linghang Kong, Fang Zhang, Jianxin Chen,
    ACM TQC 7(2), Article 14 (DOI 10.1145/3797968); arXiv:2507.12024v2.

Anchor strings like "S3.T5" are the LaTeX ids of the arXiv HTML rendering, so
each row can be re-verified by grepping the source page.  Entries split into
  * published tables  (inputs to the reproduction), and
  * derived constants (glue arithmetic over published rows, marked DERIVED).
"""

DERIVED = "[DERIVED]"   # computed from this paper's own published numbers
REPORTED = "[REPORTED]"  # printed by the paper, not reproducible from prose

WEBSITE = dict(
    arxiv="2507.12024", version="v2",
    v2_dateline="arXiv:2507.12024v2 [quant-ph] 22 Apr 2026",
    journal="ACM TQC 7(2):14, DOI 10.1145/3797968",
    authors="Linghang Kong, Fang Zhang, Jianxin Chen",
    fetched="2026-08-29, https://arxiv.org/html/2507.12024v2",
)

# ---------------------------------------------------------------------------
# Eq.(1) (S2.E1):  p0 = 1 - (1-pL)^(1/(k*d))
# Eq.(2) (S2.E2):  Q   = k / ((ceil(tr/t_SEC) + d) * t_SEC)
# QLOPS density (S2.p10): QLOPS / N, N = physical qubits of the code patch
#   (data + ancilla).  For GB codes N = 2n (S3.T5 "Physical qubit number"
#   column equals 2n for every row; verified programmatically in gate A).
# ---------------------------------------------------------------------------

# ----- Table 1 (S3.T1): superconducting platform (Willow-based; Future =
# "expected parameters in ten years") --------------------------------------
SC_HW = {
    "current": dict(coherence_us=80,    t_1q=0.025, t_2q=0.04,
                    p_1q=5e-4, p_2q=2e-3, t_ro=0.5,  p_ro=7e-3,
                    t_prep=0.1, p_prep=5e-3),
    "future":  dict(coherence_us=1000,  t_1q=0.02,  t_2q=0.03,
                    p_1q=1e-4, p_2q=5e-4, t_ro=0.1,  p_ro=2e-3,
                    t_prep=0.1, p_prep=1e-3),
}

# S3.SS1.p2: t_SEC = prep + measure + 4*1Q + 4*CZ; paper states the sums
# 0.86 us (current) and 0.40 us (future) and names the components.
T_SEC_SC = {"current": 0.86e-6, "future": 0.40e-6}
# check (component sums, in us): 0.1+0.5+0.1+4*0.025+4*0.04 = 0.86
T_SEC_SC_COMPONENTS_us = {"current": 0.86, "future": 0.40}

# ----- Table 2 (S3.T2): PyMatching t_r (s), Xeon Gold 6248R (dual-socket,
# 3.00 GHz base / 4.00 GHz max).  (a) current-hardware runs; (b) future. ----
T_R_SURFACE_CURRENT = {
    11: 1.2443e-5, 13: 2.0780e-5, 15: 3.2834e-5, 17: 9.7008e-5,  # note (i)
    19: 4.8795e-5, 21: 6.9903e-5, 23: 1.3085e-4, 25: 1.7153e-4,
    27: 2.2188e-4,
}
T_R_SURFACE_FUTURE = {
    5: 2.0955e-7, 7: 5.5807e-7, 9: 1.1682e-6,
    11: 2.1191e-6, 13: 3.4904e-6, 15: 5.4921e-6,
}
# GATE-A FINDINGS on this table (see campaigns/20260829T223540Z_dcf693ab/):
#   (F1) d=17 entry (9.7008e-5) is out-of-trend AND inconsistent with the
#        published Table-5 QLOPS: Table 5's current-SC rows for d=17/19
#        reproduce EXACTLY (dev <= 4e-8) only when t_r is taken from the
#        NEXT d row (d=19 -> 4.8795e-5, d=21 -> 6.9903e-5).  The printed
#        entries appear shifted by one table slot from d>=17 onward; the
#        sequence 3.2834, 4.8795, 6.9903, 9.7008, 13.085, 17.153, 22.188
#        (e-5, d=15..27) is smooth if re-anchored at d=15.
#   (F2) Table 5's [[288,12,18]] current-SC row implies ceil(tr/t_SEC)=259
#        while tr/t_SEC = 2.2188e-4/0.86e-6 = 258.0 exactly -> the paper
#        (or its float pipeline of the same trap documented in
#        test_rsa_exact_vs_float_ceil) took the ceiling one unit high;
#        impact 0.35% on Q, within the 5% gate-A tolerance.
# note (i): the d=17 row (9.7008e-5) is out-of-trend vs its neighbors;
# recorded verbatim as printed.  Full analysis: GATE-A FINDINGS above.
T_R_SURFACE_NOTES = {
    17: "out-of-trend in source table; entry reproduced as printed",
}

# ----- Table 3 (S3.T3): neutral-atom platform -------------------------------
NA_HW = {
    "current": dict(coherence_s=1.5, t_1q=0.5, t_2q=0.2, p_1q=1e-3,
                    p_2q=5e-3, p_unintended=1e-3, t_ro=500, p_ro=2e-3,
                    p_prep=7e-3, p_move=1e-3),
    "future":  dict(coherence_s=20,  t_1q=0.5, t_2q=0.2, p_1q=1e-4,
                    p_2q=1e-3, p_unintended=2e-4, t_ro=50,  p_ro=2e-4,
                    p_prep=2e-4, p_move=1e-4),
}
# Movement model (S3.SS2.p1, citing Xu et al. 2024):
#   t_move = sqrt(6*DX/a_p)  with lattice spacing 5 um, a_p = 0.02 um/us^2.
ATOM_MOVE = dict(lattice_um=5, a_p=0.02, formula="sqrt(6*DX/a_p)")

# ----- Table 4 (S3.T4): GB-code simulation results (future-NA platform) ----
# keys are "<n>,<k>,<d>-<decode>"; decode Z = Z-syndromes only,
# ALL = X+Z syndromes (BP-LSD, min-sum, lsd_cs order 10).
GB_TABLE4 = {
    "72,12,6-Z":   dict(pL=8.372e-4,   tr=6.33e-4, t_sec=2.677e-3),
    "90,8,10-Z":   dict(pL=1.614e-4,   tr=2.070e-3, t_sec=2.911e-3),
    "108,8,10-Z":  dict(pL=7.40e-5,    tr=2.194e-3, t_sec=2.799e-3),
    "144,12,12-Z": dict(pL=7.10e-5,    tr=4.644e-3, t_sec=2.901e-3),
    "288,12,18-Z": dict(pL=1.198e-6,   tr=2.8281e-2, t_sec=3.550e-3),
    "72,12,6-ALL": dict(pL=3.286e-4,   tr=6.6876e-2, t_sec=2.677e-3),
}

# ----- Table 5 (S3.T5): published comparison values ------------------------
# For each GB code: row 1 = neutral-atom (future) numbers; rows 2-3 =
# superconducting device matched at the same p0 (R2 current hw, R3 future).
# d_sc = the surface-code distance the paper chose (odd integers).
T5 = {
    "72,12,6-Z":   dict(code=(72, 12, 6), p0=1.1633e-5, d_sc=6,
        na_qubits=144,   na_qlops=640.35,     na_dens=4.4468,
        sc_cur_qubits=4044,  sc_cur_qlops=367197.06, sc_cur_dens=90.8005,
        sc_fut_qubits=588,   sc_fut_qlops=5000000.00, sc_fut_dens=8503.4014),
    "90,8,10-Z":   dict(code=(90, 8, 10), p0=2.0177e-6, d_sc=10,
        na_qubits=180,   na_qlops=249.80,     na_dens=1.3878,
        sc_cur_qubits=4616,  sc_cur_qlops=125707.10, sc_cur_dens=27.2329,
        sc_fut_qubits=776,   sc_fut_qlops=2222222.22, sc_fut_dens=2863.6884),
    "108,8,10-Z":  dict(code=(108, 8, 10), p0=9.2503e-7, d_sc=10,
        na_qubits=216,   na_qlops=259.85,     na_dens=1.2030,
        sc_cur_qubits=5768,  sc_cur_qlops=92102.23,  sc_cur_dens=15.9678,
        sc_fut_qubits=776,   sc_fut_qlops=2222222.22, sc_fut_dens=2863.6884),
    "144,12,12-Z": dict(code=(144, 12, 12), p0=4.9307e-7, d_sc=12,
        na_qubits=288,   na_qlops=295.45,     na_dens=1.0259,
        sc_cur_qubits=8652,  sc_cur_qlops=138153.35, sc_cur_dens=15.9678,
        sc_fut_qubits=1164,  sc_fut_qlops=3333333.33, sc_fut_dens=2863.6884),
    "288,12,18-Z": dict(code=(288, 12, 18), p0=5.5451e-9, d_sc=18,
        na_qubits=576,   na_qlops=130.01,     na_dens=0.2257,
        sc_cur_qubits=17484, sc_cur_qlops=48788.42,  sc_cur_dens=2.7905,
        sc_fut_qubits=2892,  sc_fut_qlops=1764705.88, sc_fut_dens=610.2026),
    "72,12,6-ALL": dict(code=(72, 12, 6), p0=4.5646e-6, d_sc=15,
        na_qubits=144,   na_qlops=144.59,     na_dens=1.0041,
        sc_cur_qubits=5388,  sc_cur_qlops=258397.93, sc_cur_dens=47.9580,
        sc_fut_qubits=1164,  sc_fut_qlops=3333333.33, sc_fut_dens=2863.6884),
}
# 5000000.00 and 2222222.22 are exact 5-significant-digit renders of
# 5,000,000 = 2*10^6*... (S3.T5 rows for d=5 / d=7 / d=11 at future SC) --
# see formulas.py docstring for the exact reconstruction that reproduces
# 5000000.00 = 10^7/2... ; treat these two columns as rounded.

# ----- Section 3.5 (S3.SS5): RSA-2048 case study ---------------------------
# Superconducting system (Gidney, arXiv:2505.15917, Table 5 data):
RSA_SC = dict(
    n_cold=1280, n_active=131, d=25, t_sec=1e-6, t_r=1e-5,
    qlops=4.0314e7, data_qubits=714019, qlops_per_data=56.4611,
    t_days=4.96, n_toffoli=6.5e9, factories=6, msd_rate_relation=
    "6 factories produce ~1 Toffoli state per logical cycle (S3.SS6.p3)",
)
# Neutral-atom system (Zhou et al., ISCA'25): 19e6 physical qubits;
# patches of d^2 data + (d^2-1) ancilla; 192 MSD factories of 3x12 patches;
# 6128 computational logical qubits; correlated decoding -> O(1) SEC per
# transversal op; reaction time IGNORED per S3.SS5.p3 ("we simply ignore the
# reaction time... assuming 1 syndrome extraction cycle following each
# transversal gate according their Section IV.2").
RSA_NA = dict(
    tops=19e6, d=27, factories=192, factory_patches="3x12",
    n_logical=6128, t_sec=900e-6, t_r_used=0.0,
    qlops=6.8089e6, t_days=5.6, n_toffoli=3.0e9,
    qlops_per_data=0.7626,
)
RSA_PUBLISHED = dict(
    eq3_ratio=5.244,    # Q_SC*t_SC / (Q_atom*t_atom)          (S3.E3)
    eq4_ratio=2.4204,   # same normalized by Toffoli counts     (S3.E4)
    non_clifford_per_toffoli=10,  # "a Toffoli via gate teleportation
    # requires ~10 Clifford gates (Gidney & Fowler 2019)"; "10 would still
    # be a good estimation" (S3.SS5.p6)
    underest_sc=270,    # t / (10 nT / Q) for superconducting   (S3.SS5.p6)
    underest_na=110,    # same for neutral atoms                (S3.SS5.p6)
)

# ----- Table 6 (S3.T6): Litinski 15-to-1 magic-state cost -----------------
# Protocol "(15-to-1)a,b,c" == Litinski (dX,dZ,dm) (Litinski 2019b, arXiv
# 1905.06903v3 Figs. 5/11); unit-space formula = 2*(dX+4dZ)*3dX + 4dm
# physical qubits (Litinski S3.p4, "taking physical measurement ancillas
# into account"); plain time cost = 6*dm cycles, published cycles =
# 6*dm/(1-p_fail) with p_fail the protocol's post-selection failure.
# unit + total columns reproduced below; p_dist and cycles are REPORTED.
# NOTE the source table lists TWO rows with protocol (15-to-1)9,3,3 for
# [[72,12,6]] ALL (different p_dist); suffix "-2" disambiguates.
T6 = [
    dict(code="72,12,6-Z",   proto="15to1_9,3,3",     p0=1.1633e-5,
         p_dist=2.3317e-6, unit=1146,  cycles=18.6423, total=2292),
    dict(code="72,12,6-Z",   proto="15to1_23,9,9",    p0=1.1633e-5,
         p_dist=9.0607e-6, unit=8178,  cycles=58.0795, total=155382),
    dict(code="72,12,6-Z",   proto="15to1_7,3,3",     p0=1.1633e-5,
         p_dist=1.1269e-5, unit=810,   cycles=18.6403, total=30780),
    dict(code="90,8,10-Z",   proto="15to1_11,3,3",    p0=2.0177e-6,
         p_dist=9.5818e-7, unit=1530,  cycles=18.7980, total=1530),
    dict(code="90,8,10-Z",   proto="15to1_29,9,11",   p0=2.0177e-6,
         p_dist=1.9588e-6, unit=11354, cycles=70.1846, total=90832),
    dict(code="90,8,10-Z",   proto="15to1_9,3,3",     p0=2.0177e-6,
         p_dist=1.3258e-6, unit=1146,  cycles=18.8261, total=19482),
    dict(code="108,8,10-Z",  proto="15to1_11,5,3",    p0=9.2503e-7,
         p_dist=6.8826e-7, unit=2058,  cycles=18.6779, total=2058),
    dict(code="108,8,10-Z",  proto="15to1_29,11,11",  p0=9.2503e-7,
         p_dist=5.5174e-7, unit=12746, cycles=68.1534, total=76476),
    dict(code="108,8,10-Z",  proto="15to1_11,3,5",    p0=9.2503e-7,
         p_dist=9.8766e-8, unit=1538,  cycles=31.1579, total=43064),
    dict(code="144,12,12-Z", proto="15to1_11,3,5",    p0=4.9307e-7,
         p_dist=2.5036e-7, unit=1538,  cycles=30.9210, total=3076),
    dict(code="144,12,12-Z", proto="15to1_31,11,11",  p0=4.9307e-7,
         p_dist=4.7231e-7, unit=13994, cycles=68.2899, total=125946),
    dict(code="144,12,12-Z", proto="15to1_11,3,5",    p0=4.9307e-7,
         p_dist=9.8766e-8, unit=1538,  cycles=31.1579, total=64596),
    dict(code="288,12,18-Z", proto="15to1_15,5,5",    p0=5.5451e-9,
         p_dist=2.1524e-9, unit=3170,  cycles=30.1425, total=3170),
    dict(code="288,12,18-Z", proto="15to1_39,17,15",  p0=5.5451e-9,
         p_dist=4.4065e-9, unit=25098, cycles=90.6157, total=100392),
    dict(code="288,12,18-Z", proto="15to1_13,5,5",    p0=5.5451e-9,
         p_dist=1.3677e-9, unit=2594,  cycles=30.0962, total=57068),
    dict(code="72,12,6-ALL", proto="15to1_9,3,3",     p0=4.5646e-6,
         p_dist=2.3317e-6, unit=1146,  cycles=18.6423, total=1146),
    dict(code="72,12,6-ALL", proto="15to1_25,11,9",   p0=4.5646e-6,
         p_dist=4.3460e-6, unit=10386, cycles=57.4921, total=135018),
    dict(code="72,12,6-ALL", proto="15to1_9,3,3",     p0=4.5646e-6,
         p_dist=1.3258e-6, unit=1146,  cycles=18.8261, total=29796),
]

# ----- Gate-B inputs -------------------------------------------------------
# 7-T baseline = the Litinski Table-6 rows above.
# Zero-level CCZ replacement (arXiv:2605.21867, Itogawa, Hirano, Akahoshi,
# Fujii; v1 21 May 2026, the only version).  Constants were transcribed from
# the ABSTRACT on 2026-08-29; the FULL TEXT was read first-hand 2026-09-03
# (arxiv.org/html/2605.21867v1) and every printed value below was located in
# the body, so the transcription half of the [REPORTED] tag is DISCHARGED:
#   c = 300      : Sec. IV, "Least-squares fitting (solid lines) shows that
#                  both cases satisfy p_L ~= 300 p^2", and the Fig. 10
#                  caption; also abstract, Sec. I, Sec. V.  It is a fit, not
#                  a table value (the paper prints no numbered tables).
#   phys_qubits  : 22, abstract + Sec. I + Sec. III + Sec. V.  SCOPE: 22 is
#                  the DISTILLATION CIRCUIT only; the 3 output surface-code
#                  patches are additional and are NOT inside the 22.
#   depth        : 24, abstract + Sec. I + Sec. III + Sec. V ("three rounds
#                  of surface-code syndrome extraction").
#   log_qubits   : 3, ditto (output patches, teleported from the [[8,3,2]]
#                  block by AIT lattice surgery, Sec. III.2).
# The REPRODUCTION half of the tag stays open and is now known to be
# unreachable through ../msd/: msd's reconstruction has a distance-1 output
# frame -- per-observable failure is linear in p (log-log slope 0.992,
# coefficient ~167 across p=1e-5..3e-4) -- so it measures its own floor, not
# c*p^2, and is a SCOPE-LIMITED-NON-TEST.  Reproducing 300 would need a
# distance-carrying (expanded d3/d7) rebuild; the paper ships no Stim file,
# no repository and no ancillary files (Appendix A is figure-only).
ZERO_LEVEL_CCZ = dict(
    arxiv="2605.21867", c=300.0, p_L_formula="300 p^2", p_order=2,
    phys_qubits=22, depth=24, syndrome_rounds=3, log_qubits=3,
    # Historical Gate-B arithmetic, retained only so gate_b_sweep.py can
    # reproduce the frozen artifact. This mixes physical-qubit*circuit-layer
    # units with Litinski physical-qubit*syndrome-cycle units and is not a
    # common-basis resource metric. See zero_level_provenance.py Revision 4.
    legacy_gate_b_spacetime_mixed_units=22 * 24,
    legacy_gate_b_spacetime_tag=DERIVED,
    claim="5-10x space-time reduction vs previous methods",
    # Conditions the abstract-level transcription omitted (Sec. IV):
    fit_p_range=(1e-4, 1e-3),          # six sampled p values
    fit_trials_per_p="1e7-1e8",
    fit_coefficient_uncertainty=None,  # no uncertainty/CI for fitted c
    fit_residuals=None,                # no fit residuals printed
    pointwise_error_bars=(
        "Figures 10 and 11 draw pointwise error bars; their definition and "
        "confidence level are not stated"),
    decoder=None,                      # no decoder named anywhere
    noise_model=("single-qubit depolarizing on gates/idles, flip on "
                 "init/measurement, two-qubit depolarizing; all probability p"),
    clifford_approximation=("T/T+ replaced by identity in Stim; paper states "
                            "this makes p_L 'slightly optimistic', at most a "
                            "modest correction to the prefactor"),
    # Postselected protocol: 300 p^2 applies to ACCEPTED runs only (Fig. 11).
    acceptance_p1e3=(0.30, 0.40),
    acceptance_p1e4=0.90,
    reduction_claim_evaluated_at_p=1e-3,
    tag=REPORTED,
)

# Sensitivity axes fixed by ../README + pre_statement.md:
GATE_B = dict(latency_multipliers=[0.5, 0.75, 1.0, 1.5, 2.0],
              comparability_falsifier=2.0)
