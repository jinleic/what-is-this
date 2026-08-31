# Diagnostic protocol correction (2026-08-31, before successful run)

The first attempted invocation of `scratch/gate_c_diag_penalty.py` completed the
122.5 s slope pass but crashed before producing a diagnostic result:
`gm.p_comp[r].v` is a scalar `Slope`, not a one-element list. No V1 value or
penalty-zero verdict was produced.

The pre-statement's phrase “Part-level … (the only place lam_sum enters the
aggregation; the glob blocks carry no lam_sum)” is **RETRACTED**. Glob penalties
also depend on `gm.lam_sum` through their Lemma-1 epsilon term. The successful
one-line diagnostic therefore zeros every hash-penalty RHS in the replay's
R aggregation (all level-3 Part penalties and all three glob penalties), while
leaving p_comp untouched. This is the literal Main Amendment-2 test (“zero the
penalty contribution in the replay”). It remains DIAGNOSTIC ONLY and cannot
enter any certified path.

## Post-diagnostic arithmetic correction

The pre-statement's approximate inferred R deficit `1.8866342744` is
**RETRACTED** (decimal transcription error). From the successfully reproduced
point values, the exact float-level identity is

    R_failing = 0.9303495043602614
    R_correct = 2.8170035674609757
    deficit   = 1.8866540631007143
    deficit / M_low = 0.9008714859433105

with `M_low = 2.0942543887102634`. This is the recorded `+0.9008715` V1
residual to float display.
