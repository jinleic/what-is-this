# Gate B control record — rank-25 anchor seed

- UTC freeze: 2026-08-30T01:05Z. Campaign dir frozen after write; no post-hoc edits.
- Script: `gate_b_search3.py` (sha256 6b16640a8af43145... full hash in `code_hashes.sha256`),
  the version running as `oct-rank-ladder` (pid 96931, started 17:42:08 local 2026-08-29).
- Purpose: rank-25 is the *positive control* (published rank is 25). It validates the
  ALS+LM pipeline before the ladder spends budget on ranks 24..20.
- Control run config: BASE_SEED=20260829, STARTS=1, ALS_MAX=2000, LM_MAX=300,
  CERT_TOL=1e-13, WCAP_H=8, PCG64 stream `[20260829, 25]`.
- Control result (from `gate_b_rank25_control_results.json`):
  - `best_rel = 1.848031286166752e-14` (machine precision, ||T||_F = 4.0).
  - factor column norms `.best_factor_norms = [2.3036, 9.019, 9.463]` — **bounded, O(1)..O(10)**.
  - `ratio = 5.12e14` (norm/residual, large only because residual ~ 1e-14).
  - trajectory: ALS rel 0.019 after 2000 sweeps (stall stop), then LM rounds 287 to rel 1.85e-14, status=3 (xtol).
  - norms grew modestly during LM (9.26 → 9.46) but stayed bounded — consistent with an
    exact rank-25 point existing nearby (as published and Krawczyk-certified).
- Evidence label: `COMPUTATIONAL-EVIDENCE` (numerical anchor; the mathematical fact
  R_R(T_O) <= 25 is `MACHINE-VERIFIED` separately via Gate A's Krawczyk replay).
- This control baseline (`rel ~ 1e-14 with norms O(1)`) is the contrast class the ladder
  ranks 24..20 are reported against: joint (residual, max factor norm) reporting only.
- Environment: python 3.14.3, numpy 2.5.2, scipy 1.18.1, singlethreaded BLAS, nice -n 10.
