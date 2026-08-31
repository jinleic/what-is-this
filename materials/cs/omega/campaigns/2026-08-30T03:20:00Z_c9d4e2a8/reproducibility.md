# Reproducibility

Certified result verified byte-identical on re-run (2026-08-30T03:30Z, same
machine, same process isolation):
  omega_cert_upper_raw  = 2.371340083602922   (bit-identical)
  omega_cert_with_slack = 2.3713411672715115  (bit-identical)
  eps_abs               = 1.083668589495582e-06
  params_sha256         = f23369136314cc51497d0c468c00d69a02b74d2b665073b9f934240ab428b46f

Deterministic run; no randomness, no threading; everything flows from
the released parameters' float64 bytes through dyadic conversion to
python-flint ball arithmetic with fixed precision (MID=300, OUTB=64).
