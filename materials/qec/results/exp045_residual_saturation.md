# EXP-045 residual saturation: B' on parent 9a7638586033
- Parent phase2_98 k_P=12 T=4 (exact), d_Z(P)=6 exact.
- dim V = 112 (kernel of skew(AC^T+BD^T) over the monomial span; matches EXP-044).
- cap rule: dim_V=112 > 22: enumerate every GF(2) subset sum of the canonical nullspace_np basis with support weight 0..w_max=6; chosen so # = sum C(112, k) over k in 0..6 = 2,533,006,645 instances; calibrated per-instance rate keeps the run far below the 4 h single-thread wall cap; hard guard at 14400 s sets truncated=true with partial counts.
- instances_total = 2,533,006,645; instances_increase = 0.
- dim_bar_histogram_on_increase = {}.
- ALL instances: bar_Delta histogram {'0': 1266490205, '2': 1266516440}.
- Exact full-space certificate: dim A = 2 < T = 4, M_bar not in A; absorption (hence increase) impossible over every v in V.
- Verdict: DECIDED_ON_PARENT_WITHIN_BUDGET; B' holds vacuously on this parent (no INCREASE instance exists).
- Wall 3874.826 s; enumeration 3873.63 s; integrity verified=True; cross-validation 400 dense draws, 0 pivot mismatches.
