# GENERAL-R H2: THE EXACT ONE-STEP IDENTITY, AND WHY THE INSERTION
# INDUCTION DOES NOT CLOSE

**Status of this document.**  An earlier version of it asserted that an
induction over atom insertion closes Liu's Hypothesis 2 on the whole
finite-atom paired class.  That assertion was unfounded: it rested on a
"4-slot S_asym" object that was invented, on a certifier that computed
nothing, and on a cover that scanned cell centres in floats.  Those
artifacts are deleted (see Addendum 11 of
`uc/LIU9_BLOCK_COPOSITIVE_2026-08-29.md`).  This rewrite keeps only what is
proved, states the refutation of the induction explicitly, and names the
remaining wall.

## 1. Setup

A paired law with `r` atoms is
`P0 = sum_{i<=r} a_i delta_{x_i}`, `P1 = sum_{i<=r} a_i delta_{y_i}` with
shared masses `a_i > 0` summing to 1, mixed as
`P_mix = (1-q) P0 + q P1`.  Write `Kpi(s,t) = h(pi(s,t))` with
`pi(s,t) = s t (1 + (1-s)(1-t))`, and

    c_ch^{(r)} = beta sum_{i,j<=r} a_i a_j [Kpi(x_i,x_j) + Kpi(y_i,y_j)
                                            - 2 Kpi(x_i,y_j)].

**PROVED (MACHINE-VERIFIED, symbolic, zero residual).**  The raw gap
splits exactly as

    gap = F(P_mix) + q(1-q) c_ch^{(r)},     F = discharge / M_mix,

for `r = 1` (Family A: `uc/liu9_scalar_margin_region1.py`,
`...region1b.py`), `r = 2` (`uc/liu9_scalar_margin_region2.py`) and
`r = 3` (`uc/liu9_paired_class_c.py`, sympy multilinear expansion in the
canonical free log basis, zero residual).

**PROVED.**  `F(P_mix) >= 0` for every mixture mean `>= m` -- the
universal q = 1 theorem, 524,800-box centered Arb cover, with
`dC_M/dM >= 0`.

## 2. The one-step insertion identity (PROVED)

Insert one atom pair `(u, v)` at shared relative mass `t`, rescaling the
existing masses by `(1-t)`.  Then, exactly,

    c_ch^{(r+1)} = (1-t)^2 c_ch^{(r)} + t^2 D(u,v) + 2 t (1-t) X(u,v)

    D(u,v) = beta [Kpi(u,u) + Kpi(v,v) - 2 Kpi(u,v)]
    X(u,v) = beta sum_{i<=r} a_i [Kpi(x_i,u) + Kpi(y_i,v)
                                  - Kpi(x_i,v) - Kpi(y_i,u)]

so the channel increment is **exactly quadratic in `t`**.  Verified as a
zero residual over a free basis on the `h(pi(.,.))` values at
`r = 1, 2, 3`, with three mutations that must and do fail
(`uc/liu9_channel_insertion.py`,
`uc/verification/results/liu9-channel-insertion.json`).

Combined with section 1 this gives the exact one-step gap identity

    gap^{(r+1)} - gap^{(r)} = [F^{(r+1)} - F^{(r)}]
        + q(1-q) [ (t^2 - 2t) c_ch^{(r)} + t^2 D(u,v) + 2t(1-t) X(u,v) ].

Both brackets are exact.  Neither is signed.

## 3. The induction does not close (REFUTED)

Three separate facts, each certified, block the induction that the earlier
draft claimed:

1. **The channel increment is sign-indefinite.**  `D(u,v)/beta` is
   enclosed strictly below zero by Arb at exact rational points:
   `[-0.633399360587514216356264 +/- 2.41e-26]` at `(u,v) = (7/20, 39/40)`,
   `-0.0769642406944118250747776` at `(9/10, 1/10)`, and exactly `0` at
   `(1/2, 1/2)`.  This is the kernel-level restatement of the recorded
   fact that `k_pi` is not PSD.
2. **S_asym is not the increment.**  The certified pencil
   `S_asym = (1-q)S(y1) + qS(y2) + q(1-q)Chat(y1,y2)` of
   `uc/liu9_second_order.py` shares none of its six protocol coefficients
   with `D`, and additionally carries product entropies `Hprod(.,.)`,
   the single entropies `h(y1)`, `h(y2)` and the first-order constant
   `A_const`.  So the certified `S_asym >= 0` theorem cannot be cited to
   sign the channel increment.
3. **The smooth-chart margin does not dominate the channel rate.**  For
   the single-pair split `P0 = delta(x*) + delta(x*+d)`,
   `P1 = delta(x*) + delta(x*-d)` with shared masses `(1/2, 1/2)`, the
   channel rate converges to `E_K/d^2 = -1.5926157`, giving
   `|rate|/kappa = 12.97` against the smooth-chart constant
   `kappa = 4119063/33554432` (`uc/probe_single_pair_ek.py`, Addendum 10).

Consequently there is no per-step nonnegativity to compose, and the
"weak-limit composition" step of the earlier draft is moot: there is
nothing certified to pass to the limit.

## 4. What survives, and what the wall is

Surviving, in the ledger's own labels:

* the exact block diagonalization `gap = F + q(1-q) c_ch` (PROVED, r <= 3);
* the general channel formula for the finite paired class (PROVED);
* `F(P_mix) >= 0` for mean `>= m` (PROVED);
* `S_asym >= 0` on the full cube (PROVED in `liu9_second_order.py`;
  independently replicated over 1,388,611 cell pairs in
  `uc/liu9_sasym_cover_check.py`);
* Family A on `(0, 1]` and the r = 2 curve families (PROVED);
* r = 3: 7644 exact pinch configurations with Arb-certified positive gap
  and `2^-12` neighbourhood boxes (MACHINE-VERIFIED finite);
* the full-half-space tube at `rho = 1/1701` (PROVED).

The wall is unchanged from what
`uc/verification/results/liu9-block-kernel.json` already recorded as
`OPEN_REDUCTION`: **prove that the F-margin dominates the negative part of
the channel**, i.e. copositivity of the explicit `2x2` size-biased block
kernel on positive paired measures with shared masses.  Section 2 sharpens
what has to be dominated -- an exactly quadratic-in-`t` increment whose
`t^2` coefficient can reach `-0.633 beta` -- but sharpening is not
closing.

## 5. Honest label table

| statement | label | evidence |
|---|---|---|
| `gap = F + q(1-q) c_ch`, r <= 3 | PROVED | `liu9-block-copositive-diagonalization.json`, `liu9-paired-class-c.json` |
| general channel formula | PROVED | `liu9-paired-class-c.json` |
| `F(P_mix) >= 0`, mean `>= m` | PROVED | universal q = 1 theorem |
| `S_asym >= 0` on the cube | PROVED | `liu9_second_order.py`; replicated in `liu9-sasym-cover-check.json` |
| one-step channel identity, quadratic in `t` | PROVED | `liu9-channel-insertion.json` |
| channel increment `>= 0` | **REFUTED** | `liu9-channel-insertion.json` (Arb, exact points) |
| `S_asym` = channel increment | **REFUTED** | `liu9-channel-insertion.json` (coefficient comparison) |
| insertion induction closes H2 | **REFUTED** | sections 2-3 of this document |
| r = 3 enumerated pinch grid | MACHINE-VERIFIED finite | `liu9-scalar-margin-r3.json` |
| general-r H2 | **OPEN** | blocked on block copositivity |

## 6. Certificates cited

* `uc/verification/results/liu9-block-copositive-diagonalization.json`
* `uc/verification/results/liu9-scalar-margin-region1.json`
* `uc/verification/results/liu9-scalar-margin-region1b.json`
* `uc/verification/results/liu9-scalar-margin-region2.json`
* `uc/verification/results/liu9-scalar-margin-region2-cover.json`
* `uc/verification/results/liu9-scalar-margin-region2-hulls.json`
* `uc/verification/results/liu9-scalar-margin-region2-curve2.json`
* `uc/verification/results/liu9-paired-class-c.json`
* `uc/verification/results/liu9-channel-insertion.json`
* `uc/verification/results/liu9-sasym-cover-check.json`
* `uc/verification/results/liu9-scalar-margin-r3.json`
* `uc/verification/results/liu9-block-kernel.json` (the open reduction)
