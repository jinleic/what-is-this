# Reproduction (rung 1, corrected endpoint)

From repo root /Users/jinleic/jinleic-workspace/cs:

  cd cs/omega
  /Users/jinleic/jinleic-workspace/cs/.venv/bin/python src/stage_b_rung1_certified.py

Expected (from the original certified_result.json, with the CORRECT
interpretation):
  raw omega endpoint declared by the certified_result.json:
      omega_cert_upper_raw = 2.371340083602922

  corrected tight certified endpoint (add absorbed 6.5978e-11):
      omega_cert = 2.3713400836689

  Note: stage_b_rung1_certified.py as-frozen reproduced the LOOSE variant
  that double-counted the Lemma-1 term (it calculated
  omega_cert_with_absorbed_slack = 2.3713411672715115). The CORRECTION note
  in this directory documents the coincidence resolution and the corrected
  arithmetic, which was owner-verified independently at 200 bits.

To reproduce the DECISIVE coincidence test (substring of the version we used):
  /tmp/coincidence_test.py  with|without
  (script content saved inline in this campaign dir for archival)

## TIGHT CERTIFIED STATEMENT

  omega ≤ 2.3713400836689
  at the exact released float64 parameters W1.00_2.371339.mat of Alman25
  SODA'25, with the combined absorbed c+ceq+Schonhage defects contributing
  6.5978e-11, and the Lemma-1 gap implicit in the raw interval endpoint.
