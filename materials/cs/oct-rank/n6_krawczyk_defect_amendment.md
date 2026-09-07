# N6 precompute amendment — Krawczyk remainder defect and mandatory re-audit

Run `20260904T052006Z_e5f2adce_0e7798f95678`; gate
`n6-qi-real13-descent`; preregistration commit `f8187e1`, SHA-256
`f6c7648a1a3a39ae84ded246b5dc45bbc4f7891890cc3821909d70e46b455564`.
Bound before any N6 scientific computation. This amendment supersedes the
preregistration/provenance description of the Route-F core as “validated.” N6
is blocked until the corrected-core re-audit below passes and receives review.

## 1. Defect

For one CP summand `f(a,b,c)=abc`, with center `(a0,b0,c0)` and displacement
`(da,db,dc)`, direct expansion gives

```
f(x0+d)-f(x0)-Df(x0)d
  = da*db*c0 + da*dc*b0 + db*dc*a0 + da*db*dc.
```

Consequently the exact componentwise absolute remainder bound is

```
Pa*Pb*|c0| + Pa*Pc*|b0| + Pb*Pc*|a0| + Pa*Pb*Pc.
```

Frozen source SHA-256
`781e51ec1a72cfb8acf7e7c7b1e44a96cac011ec6f46dc61e88f5b7e01b87035`
instead computed the first two terms as
`Pa*(Pb*|b0| + Pc*|c0|)`: the `b0` and `c0` weights are swapped. This can
under-bound the Taylor remainder when one of the `b`/`c` coordinates is held
free. No containment or exclusion produced by those source bytes is accepted
without corrected replay.

The original bytes are preserved in the live N6 run as
`rf_krawczyk_core_buggy_781e51.py` with the same SHA-256. The live corrected
copy changes only that coefficient placement and has SHA-256
`7daff00369e2ad3b368428e7a349b15b8ad37b397cacfa845ed56e03d2873e1b`.
Frozen runs remain byte-untouched.

## 2. Independent formula control

Before any corrected certificate is accepted,
`remainder_formula_control.py` (SHA-256
`613b2c4ff8c63fd8513fd5244c6dfc4487b5b1f29e5585a807691687b17ff7ba`)
must pass. It does not import the Krawczyk core. It checks direct polynomial
subtraction against an independently written expansion, exhausts every
solved/free mask and rational box corner/zero choice, and checks the `|Y|W`
triangle propagation.

Its planted asymmetric case is `a0=2,b0=3,c0=11, da=5,db=7,dc=0`.
The exact remainder and corrected bound are both `385`; the historical formula
returns `105`, so the historical formula is required to be rejected as an
actual under-bound. A reciprocal `da,dc` case and a case exercising all
quadratic/cubic terms are also mandatory.

## 3. Frozen-use inventory and corrected replay

`krawczyk_defect_reaudit.py` (SHA-256
`31d115af836155664f5b38565d6fadb568d18ad85243da2e586c6013ba54a437`)
must inventory every Python file whose complete bytes hash to the defective
source and re-adjudicate every historical call to `certify`:

1. frozen Route-F `MAIN`, `TAU7`, `TAU6`, `PPOS`, and `PWRNG`;
2. Route-F replay's separately parsed decimal `MAIN` center and `TAU7` call;
3. frozen N1 controls `TAU7`, `TAU6`, `PPOS`, and `PWRNG`; and
4. whether N1 admitted any main candidate to an exact/Krawczyk adjudication.

The audit reconstructs Route-F's historical binary-float dyadic `MAIN` center,
the replay's distinct exact-decimal center, exact tau factors, both deterministic
five-round Route-F and three-round N1 tau-6 centers, deterministic P-POS, and
the common P-WRNG center from frozen inputs. It runs with
`PYTHONDONTWRITEBYTECODE=1`, nice priority at least 10, CPU limit 1800 seconds,
and RSS limit 4 GiB. It writes append-only live N6 audit artifacts; it never
imports or edits a frozen run.

Acceptance requires both positive controls to retain strict containment with
all-coordinate outward-Arb confirmation; every registered negative/local probe
to retain no containment plus at least one exact box-exclusion coordinate; all
historical classifications to be unchanged; and the frozen N1 result to show
zero admitted main candidates. Any failure blocks N6 and requires an explicit
claim correction. Exact radii or other evidence details that move are recorded
even when the classification is unchanged.

## 4. Re-reviewed theorem obligations

For a selected square chart, all free coordinates are fixed exactly at the
center. `J_S` must be inverted exactly and `Y*J_S=I` checked. The centered
linear term therefore cancels. The corrected componentwise `W` bounds every
quadratic/cubic equation remainder for the solved-coordinate infinity box;
exact multiplication by `|Y|` bounds its image. Strict
`|Yg(x0)|+|Y|W < rho` in every solved coordinate maps the closed box into its
interior and gives box-local existence by the fixed-point argument. Conversely,
`|Yg(x0)| > rho+|Y|W` in any coordinate excludes a zero from that fixed box.
The QRCP rank/pivot thresholds and `rank(J_F)` are diagnostic/proposal data;
the exact inverse identity is authoritative. Outward Arb rechecks every strict
containment margin. Exclusion is box-local only and supports no global
infeasibility or rank lower-bound inference.

No N6 continuation, plant, target solve, or rank claim may run until the source
review and corrected frozen-use re-audit both pass. If cleared, N6 must use and
hash-check only the corrected core, retain the separate full-path positive and
corrupted-box controls from `n6_precompute_amendment.md`, and independently
replay any target containment.

## 5. Audited-snapshot to final N6-core cutover

The corrected source SHA-256 `7daff003...` above is preserved byte-for-byte as
`rf_krawczyk_core_corrected_audited_7daff0.py` in the N6 run. Before any N6
plant, continuation, or target computation, the final live core was extended
only to complete every rung in N6's explicit 15-rung ladder, serialize exact
per-coordinate margins, and make any outward-Arb disagreement fatal to the
complete classification. Final live-core SHA-256:
`3560f6c38a93f059dc1c0ecea15a39fe68508d208f2352cb872a33563f763fc0`.
The corrected remainder expression and audited mathematical logic are
unchanged. Final driver/replay SHA-256 values are respectively
`0b4be99ebba9d7d5323a2fe7c202b4aecc66c72dec8aad3f86519babe6a2b609`
and
`c7dbe09c54f227b984af3bed8edefa2d04c41012b4b6edbf4eb476c24e12d919`.
