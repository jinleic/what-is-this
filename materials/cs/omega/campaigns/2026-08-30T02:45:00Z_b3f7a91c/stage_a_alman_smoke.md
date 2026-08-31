/** Gate B, rung 1 (Alman25 omega <= 2.371339): stage-(a) smoke log.

Command: .venv/bin/python src/alman25_float.py <abs>/scratch/alman_code/data/W1.00_2.371339.mat

stdout:
registered params: 24855  file params: 24855
loaded omega in file: 2.3713389005
max c violation      : 1.137411e-10
max |ceq| violation  : 2.389974e-11
value (Schonhage LHS): 7.7836405962
target 2^(L-1) log(q+2) = 7.7836405962
GetFeasibility maxViol (c, |ceq| combined): 1.137411e-10

Detail run (same program, json digest at /tmp/stage_a_alman.json):
{
 "rung": "alman25_2.371339",
 "param_file": "W1.00_2.371339.mat",
 "params_sha256": "f23369136314cc51497d0c468c00d69a02b74d2b665073b9f934240ab428b46f",
 "n_params": 24855,
 "registered": 24855,
 "count_match": true,
 "omega_in_file": 2.3713389005434182,
 "published": 2.371339,
 "abs_delta_omega_published": 9.94565816370141e-08,
 "max_c_violation": 1.1374110912366195e-10,
 "max_ceq_violation": 2.389974329553013e-11,
 "combined_max_violation": 1.1374110912366195e-10,
 "verify_tol_1p1e6_ok": true,
 "refine_tol_1p1e9_ok": true
}

Stage-(a) verdict: REPRODUCED within 1.1e-9 (their own refine tolerance).
Interpretation per pre_statement Addendum A2: transcription faithful; stage (b)
(interval) failure, if any, is evidence about the bound, not the transcription.
*/
