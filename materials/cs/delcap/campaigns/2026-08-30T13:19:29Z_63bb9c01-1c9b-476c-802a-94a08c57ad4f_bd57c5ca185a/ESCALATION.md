# ESCALATION — TNB Table I C_(q,n) column vs certified enclosures

FROZEN, escalated to Main, NOT written into any README or index.

Their C_(q,n) column is their Blahut-Arimoto value for the finite-block
capacity, printed to 3 dp. Our certified LOWER endpoint is the exact mutual
information of an explicit rational input distribution, evaluated in
outward-rounded Arb at 400 bits, so this comparison needs only the primal
certificate and does not depend on the dual at all.

Rows where cert_lo exceeds the printed C by MORE than one 3-dp ulp (1.0e-3),
i.e. beyond what any printing convention can explain: 9 of 15.
Rows where the excess is within one ulp (explainable as truncation printing): 5.

  q=2 n=3 d=1/5: printed C=0.676 certified [0.680473561, 0.680473562] cert_lo-printed=+0.004474
  q=2 n=5 d=1/10: printed C=0.786 certified [0.787495955, 0.787495959] cert_lo-printed=+0.001496
  q=2 n=5 d=1/5: printed C=0.613 certified [0.620298060, 0.620298068] cert_lo-printed=+0.007298
  q=2 n=10 d=1/10: printed C=0.728 certified [0.730086455, 0.730086624] cert_lo-printed=+0.002086
  q=2 n=10 d=1/5: printed C=0.531 certified [0.542096537, 0.542096780] cert_lo-printed=+0.011097
  q=3 n=3 d=1/10: printed C=1.329 certified [1.330307006, 1.330307014] cert_lo-printed=+0.001307
  q=3 n=3 d=1/5: printed C=1.101 certified [1.105231693, 1.105231699] cert_lo-printed=+0.004232
  q=3 n=5 d=1/10: printed C=1.279 certified [1.280432698, 1.280432727] cert_lo-printed=+0.001433
  q=3 n=5 d=1/5: printed C=1.020 certified [1.026499717, 1.026499762] cert_lo-printed=+0.006500

NO THEOREM OF THEIRS IS CONTRADICTED: their stated ordering is
LB1 <= LB2 <= LB+ <= C_{q,n} <= UB, and our certified interval lies strictly
inside their own [LB+, UB] in all 15 compared rows, which CONFIRMS the
sandwich. What the numbers say is that the printed C column under-states the
true finite-block capacity. The excess grows monotonically in d (1/20: up to
1.0e-3; 1/10: 6.8e-4 to 2.1e-3; 1/5: 4.2e-3 to 1.1e-2) and in n at fixed d
(q=2, d=1/5: n=3 4.5e-3, n=5 7.3e-3, n=10 1.1e-2), which is the signature of
an under-converged Blahut-Arimoto run. That reading is INFERENCE; the
certified intervals are MACHINE-VERIFIED.
