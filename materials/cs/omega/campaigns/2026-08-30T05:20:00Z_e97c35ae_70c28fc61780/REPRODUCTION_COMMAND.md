# Reproduction

From this directory (paths resolved relative to cs/omega/):

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 10 \
      /Users/jinleic/jinleic-workspace/cs/.venv/bin/python \
      campaigns/2026-08-30T05:20:00Z_e97c35ae_70c28fc61780/stage_b_rung2.py

append `without` for the no-Lemma-1 arm.

Inputs: ../../../../scratch/rmmcode/data/K100_2.37155181.mat (inside repo,
sha256 df75ae3acaa5b1388bd2f17cce266aec169d874c9fea1eb0722c3e25871da93d).
Scripts in this dir are byte-copies of src/vxxz24_float.py,
src/stage_b_rung2.py, src/interval_core.py at run time
(sha of concatenation 70c28fc61780, embedded in the dir name).

Recorded outputs in this dir:
- stage_b_rung2_with_lemma1.json    (primary, pre-committed semantics)
- stage_b_rung2_without_lemma1.json (double-count diagnostic arm)

Both re-run byte-identical on 2026-08-30T05:3xZ (omega_cert_upper_raw
2.3715538358350807 / 2.3715518061863814, R_sum and M identical).

Environment: macOS arm64, python 3.14.3, python-flint 0.9.0 (Arb),
numpy/scipy for .mat load + float64 reference only. MID=300, OUTB=64
(interval_core). No randomness anywhere; single deterministic process.
