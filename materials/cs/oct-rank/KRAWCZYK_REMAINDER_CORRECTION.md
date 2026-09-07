# Append-only correction — Route-F/N1 Krawczyk remainder

Date: 2026-09-04. Corrective audit run:
`campaigns/20260904T052006Z_e5f2adce_0e7798f95678`. Frozen campaigns are not
modified.

## Defect and correction

Frozen source SHA-256
`781e51ec1a72cfb8acf7e7c7b1e44a96cac011ec6f46dc61e88f5b7e01b87035`
used

```
Pa*(Pb*|b0| + Pc*|c0|) + Pb*Pc*|a0| + Pa*Pb*Pc
```

for the centered remainder of `abc`. Direct expansion requires

```
Pa*Pb*|c0| + Pa*Pc*|b0| + Pb*Pc*|a0| + Pa*Pb*Pc.
```

The old expression can under-bound when a `b` or `c` coordinate is held free.
It is therefore unsound, not merely loose. The original bytes remain preserved;
the corrected live core has SHA-256
`7daff00369e2ad3b368428e7a349b15b8ad37b397cacfa845ed56e03d2873e1b`.
An independent exact control evaluates an asymmetric pair with exact/correct
remainder `385` versus historical bound `105`, exhausts all eight solved/free
masks and rational box corners, and checks `|Y|W` propagation. Its SHA-256 is
`613b2c4ff8c63fd8513fd5244c6dfc4487b5b1f29e5585a807691687b17ff7ba`.
Independent static review found the corrected residual, Jacobian, chart,
remainder, propagation, containment, exclusion, and Arb logic sound, with no
additional theorem-critical defect (`krawczyk_bound_review.json`).

The exact defective-source inventory is:

1. frozen Route-F `rf_krawczyk.py`;
2. frozen N1 `n1_instrument_base.py`; and
3. the deliberately preserved live copy `rf_krawczyk_core_buggy_781e51.py`.

## Fresh corrected classifications

Corrected re-audit result SHA-256:
`f48d264de479f711cec75ad94a62d92a0b968229302cce386e5cf129f1fa1b9d`.
Verdict: `PASS-CLAIMS-UNCHANGED`; Python 3.14.3, NumPy 2.5.2, nice 10,
`PYTHONDONTWRITEBYTECODE=1`, `RLIMIT_CPU=1800`, RSS cap 4 GiB, 187.83424 CPU s.

| use | fresh corrected classification |
|---|---|
| Route-F MAIN, historical binary-dyadic center | `NO-CONTAINMENT-ANY-RUNG`; first box-local exclusion at `rho=1/10000`; all 192 coordinates exclude at `rho<=1/1000000` |
| Route-F replay MAIN, its distinct exact-decimal center | `NO-CONTAINMENT-ANY-RUNG`; first box-local exclusion at `rho=1/10000` |
| Route-F/N1 TAU7 | strict `CONTAINMENT` at `rho=1/1000`, all-coordinate outward-Arb PASS |
| Route-F TAU6, deterministic five-round reconstruction | `NO-CONTAINMENT-ANY-RUNG`; first box-local exclusion at `rho=1/100`; all 48 coordinates exclude at `rho<=1/1000` |
| N1 TAU6, deterministic three-round reconstruction | `NO-CONTAINMENT-ANY-RUNG`; first/full 48-coordinate box-local exclusion at `rho=1/100` |
| Route-F/N1 PPOS, deterministic reconstruction | strict `CONTAINMENT` at `rho=1/10000`, all-coordinate outward-Arb PASS |
| Route-F/N1 PWRNG | `NO-CONTAINMENT-ANY-RUNG`; first box-local exclusion at `rho=1/1000`, all 192 coordinates exclude at `rho<=1/10000` |

The decisive classifications and first containment/exclusion rungs are unchanged.
Some formula-sensitive intermediate counts are corrected and supersede the
frozen counts: MAIN exclusions at `rho=1/10000`, `44 -> 36`; TAU7 contained
coordinates at `rho=1/10`, `1 -> 0`, and `rho=1/100`, `45 -> 36`; Route-F
TAU6 exclusions at `rho=1/100`, `40 -> 39`; PPOS contained coordinates at
`rho=1/1000`, `150 -> 143`. PWRNG rung counts are unchanged.

Route-F's serialized MAIN center is exactly reconstructed under its historical
binary-float-to-dyadic semantics; its decimal replay center and TAU7 center are
also replayed from serialized exact text. TAU6 and PPOS did not serialize their
centers. Those controls are therefore explicitly labeled deterministic
reconstructions, not byte-exact historical-box replays. They use the frozen
source/inputs/seeds and matching Python/NumPy/SciPy environment; their free
columns, chart rule, exact bounce maximum, ranks, and complete classifications
match the frozen record, and Route-F TAU6's final residual matches exactly.

All exclusions above are only for the named fixed-free boxes. They are not
global infeasibility certificates and give no rank lower bound. In particular,
TAU6 global rank infeasibility, where used, comes from the separate rank-seven
lower chain, not from these box exclusions.

## Independent claims that do not consume the defective bound

* N1 recorded `n_admitted=0` and an empty admitted-outcomes list, so no N1 main
  candidate ever reached exact Newton or `certify`; its 40-seed residual/no-
  admission observation is independent. Fresh corrected controls re-establish
  the gate around that observation.
* Route-F's exact block-equivalence identities, direct rank-14 concatenation,
  and frozen lower-13 chain do not consume source hash `781e51...`; hence the
  specific triple remains open in `{13,14}`. The corrected MAIN replay
  independently re-establishes the narrower failure-to-certify result only for
  the selected CP center under its historical and replay rational parsings.
* N3/N4/N5 exact substitutions, flattenings, and algebraic identities do not
  import this core. N5's `rank_C(T_H)=7` closeout is unaffected.
* The published octonion window `18 <= R_R(T_O) <= 25` is unaffected.

The historical Route-F/N1 source bytes and outputs remain preserved as a defect
record. For affected containment/exclusion evidence, this corrective audit—not
the historical `behavior_ok` flag or frozen detailed rung counts—is the current
authority.

## Collateral frozen-record corrections

The automated string `PASS-CLAIMS-UNCHANGED` is narrow: it compares the
formula-sensitive classifications that were actually executed. Independent
record review found four additional scope/provenance defects:

1. Route-F preregistered the 13-rung ladder beginning at `rho=1e-2`, but the
   shared Route-F/N1 core executed two extra unregistered rungs, `rho=1` and
   `rho=0.1`. The corrected audit inventoried all 15 actual rows. Both extras
   were non-decisive: every first containment/exclusion rung reported above is
   inside the registered 13-rung suffix.
2. N1 preregistered PWRNG as using the “best N1 candidate,” although it ran
   before the sweep. The code actually loaded Route-F's selected-CP candidate
   CSV, and that file was absent from N1 `seed_input_hashes`. PWRNG is
   reclassified as an **unregistered-center, box-local diagnostic**, not a
   valid N1 preregistered control. The corrected result establishes the actual
   diagnostic box only; it does not repair the registration defect.
3. N1's frozen verdict names R18 as the best seed; authoritative
   `n1_results.json` identifies R22. The recorded numeric minimum is unchanged.
4. Route-F polished CP and EXT but called `certify` only on selected CP.
   Replay's exact-decimal center is a different rational parsing of that same
   selected CP CSV, not an EXT center. Historical README/root-results/N2-prereg
   language claiming certified exclusions around “both” polished candidates
   is retracted. EXT received no Krawczyk adjudication and was not serialized.

These are append-only corrections; the frozen Route-F/N1 files remain intact.
N6 binds its complete 15-rung list explicitly in
`n6_record_scope_amendment.md`, uses its own preregistered same-center plant
controls, and does not inherit the N1 PWRNG or selected-CP/EXT ambiguity.

## N6 final-core cutover (precompute clarification)

The SHA-256 `7daff003...` named above is now the immutable **audited corrected
snapshot**, preserved as
`campaigns/20260904T052006Z_e5f2adce_0e7798f95678/rf_krawczyk_core_corrected_audited_7daff0.py`.
Before any N6 plant, continuation, or target computation, the live N6 core was
extended only to traverse all 15 preregistered rungs after first containment,
serialize exact per-coordinate margins, and make any outward-Arb disagreement
fatal to the complete classification. That final live
`rf_krawczyk_core.py` has SHA-256
`3560f6c38a93f059dc1c0ecea15a39fe68508d208f2352cb872a33563f763fc0`.
The corrected remainder formula and all audited chart/residual/Jacobian logic
are unchanged from the `7daff003...` snapshot.

N6 executable source hashes bound before launch are:

* `n6_qi_real13_descent.py`:
  `0b4be99ebba9d7d5323a2fe7c202b4aecc66c72dec8aad3f86519babe6a2b609`;
* `n6_replay.py`:
  `c7dbe09c54f227b984af3bed8edefa2d04c41012b4b6edbf4eb476c24e12d919`;
* `remainder_formula_control.py`:
  `613b2c4ff8c63fd8513fd5244c6dfc4487b5b1f29e5585a807691687b17ff7ba`.

The T5 split preregistration wording is separately corrected and bound in
`n6_t5_split_amendment.md`, SHA-256
`88378ef0dfb88921d97fa22a5d3a71e036ecda38cca4c8906ba01b433b08f361`.
