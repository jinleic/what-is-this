# N6 precompute amendment — frozen N3 T5 split semantics

Date: 2026-09-04. Run: `20260904T052006Z_e5f2adce_0e7798f95678`.
This amendment is bound before any N6 plant, continuation, or target computation.
It corrects one factor-name/transformation error in `n6_prereg.md` while preserving
the preregistered term, split vector, ordering, and no-search rule.

## Error and exact correction

`n6_prereg.md` (SHA-256
`f6c7648a1a3a39ae84ded246b5dc45bbc4f7891890cc3821909d70e46b455564`),
lines 46-49, incorrectly says that N3 T5 splits the conjugated **B-factor** and
maps the split vector by `S0^-1`. The frozen N3 source instead splits the
conjugated **C/Q-factor**:

```
(a, P, Q) -> (a, P, Q + W), (-a, P, W),   W = e_3
```

The two children sum identically to the original tensor term. N3 then applies
its asymmetric pullback

```
b_orig = S0^-1 P,    c_orig = S0^T Q.
```

Consequently, in original-frame coordinates the added C-vector is
`S0^T e_3`, i.e. row 3 (zero-based) of the serialized `S0`. N6 must use that
C-vector, leave the B-vector unchanged, and preserve the frozen child ordering.

## Frozen authority and fixed branch

The authority is frozen N3 run
`20260904T034332Z_f5a61843_0333458ed434`:

* `n3_instrument.py`, SHA-256
  `d529cce2938af06645529f624c58a877a9fc3c45f6bc1efc8affafa16cca94ef`,
  lines 600-626 and 659-692, specifies the asymmetric pullback and C/Q split;
* `witness_tf_c12_qi.json`, SHA-256
  `251e00577b43b400aa97d306bb923a98eb79ca7850ac67652d45288c2e570037`,
  supplies the twelve terms and `S0_rows`;
* `VERDICT.md`, SHA-256
  `78987c16dd8b2435e9a9b58b5b0bdfd3154faf5ba4d245afa7052f30c4a3fdb6`,
  lines 33-37, records that the registered first-valid branch selected term
  `r=0` immediately.

There is no new split choice or fallback in N6: term zero is used exactly.
Before tracking, N6 must separately require 192/192 exact substitution for the
original twelve-term point, its coefficientwise conjugate, the corrected
thirteen-term T5 point, and its coefficientwise conjugate. Each of those four
points must reject `TF[0,0,0] += 1`, with its exact mismatch-index set recorded.
Any failure is `INVALID-INSTRUMENT`, and TF continuation must not run.

This amendment supersedes only the erroneous B-factor/S0-inverse wording in
N6 preregistration. All other preregistered gates and stopping rules remain
unchanged. Frozen N3 bytes remain untouched.
