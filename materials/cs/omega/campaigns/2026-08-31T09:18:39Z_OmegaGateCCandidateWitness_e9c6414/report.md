# OmegaGateCCandidateWitness final report

## Headline

**All five exact dyadic points on the predeclared LP ray were evaluated with
the full branch-safe Arb objective; every full difference interval overlaps
zero. Verdict: OPEN. No improving witness and no local-optimality/no-improvement
claim.**

Pre-statement commit: `d739c45`. The initially frozen static source
`8ad6062` is preserved as superseded/static-only. A pre-compute correction
aligned the stop rule, and corrected source commit `1ec82a2` is the as-run
source. No objective evaluation occurred before that correction.

## Exact witness controls

With columns of V as kernel coordinates and rows of A as the 27 frozen margin
constraints, exact integer multiplication used

`u=[1,-1,0,-1,1,0,1,0,0,1,-1,-1,-1,1,0,-1,0,0,-1,0,1]`

and obtained the frozen 45-entry `d=V*u` in `protocol_static_corrected.json`.
All 27 entries of `A*d` and `sum(d)` equal zero exactly
(**MACHINE-VERIFIED**). Every direction entry lies in `{-1,0,1}`. All exact
rational moved coordinates for the fixed ladder `R/{16,8,4,2,1}` are frozen;
positivity, unchanged normalization/margins, `|eps*d|_inf<=R`, and unchanged
non-target parameters passed using exact integer/rational arithmetic
(**MACHINE-VERIFIED**).

The base run encloses both registered anchors. Its decomposition is:

- `R in [2.8170035674609752, 2.817007818061672]`
  (**MACHINE-VERIFIED**);
- `M in [2.0942543887102625, 2.0942543887102634]`
  (**MACHINE-VERIFIED**);
- raw Omega in `[2.371551806186382, 2.371553835835082]`
  (**MACHINE-VERIFIED**);
- equality residual in
  `[3.076684957718638e-6, 3.076684957718639e-6]`
  (**MACHINE-VERIFIED**);
- full absorbed defect in
  `[1.4691097749046072e-6, 1.4691097749046076e-6]`
  (**MACHINE-VERIFIED**);
- full objective in `[2.3715532752961566, 2.371555304944857]`
  (**MACHINE-VERIFIED**).

## All five common-expression differences

The consumed difference uses the shared-target cross product registered in
the pre-statement, followed by the Arb defect difference. Every interval below
is **MACHINE-VERIFIED** and every strict-negative test is false.

| step | raw(new-base) | full(new-base) | upper < 0? |
|---|---|---|---|
| `R/16` | `[-2.025168026027451e-6, 2.033229938243452e-6]` | `[-2.023173821748832e-6, 2.035224142522071e-6]` | no |
| `R/8` | `[-2.0206855735124868e-6, 2.037191647907985e-6]` | `[-2.016695072801931e-6, 2.041182148618541e-6]` | no |
| `R/4` | `[-2.011720657664221e-6, 2.045848749605194e-6]` | `[-2.0037375640758724e-6, 2.0538318431935425e-6]` | no |
| `R/2` | `[-1.9937907824447623e-6, 2.063181787671388e-6]` | `[-1.977822503045151e-6, 2.0791500670709997e-6]` | no |
| `R` | `[-1.9579308580745946e-6, 2.0979232018237934e-6]` | `[-1.9259922068297443e-6, 2.1298618530686437e-6]` | no |

## Why the midpoint-negative signal did not survive

The midpoint LP value at R was `-1.577332901516852e-7`
(**COMPUTATIONAL-EVIDENCE**). It selected a direction but could not resolve
point-level branch semantics. At the exact base point, Lemma-1 residual
intervals straddle zero. Consequently multiple R candidates remain possible:
for example, `R_comp[0,3]` is
`[0.24530661807075185,0.2453066195454394]` and `R_glob[0]` is
`[0.49684932826618944,0.4968507617380738]`
(**MACHINE-VERIFIED**). Summed over all seven components, R spans about
`4.2506006968e-6`, which becomes about `2.03e-6` in raw Omega. The direct
equality residual is also about `3.076685e-6`.

At R the full difference width is `4.055854059898388e-6`
(**MACHINE-VERIFIED**), `25.713367520566205` times the absolute midpoint
candidate (**INFERENCE**). Shared-target subtraction removes the common T
expression, but it cannot cancel genuine interval uncertainty about which
penalty/minimum branch is active or the full defect term. Therefore the
negative float midpoint is consistent with every result but establishes no
sign.

## Verdict and exact scope

**OPEN:** none of the five registered rungs has strict
`upper(full_new_minus_base)<0`. Overlap is not no-improvement; this run neither
exhibits a witness nor proves local optimality.

This campaign swept exactly the five points `p_star+eps*d` for
`eps=R/{16,8,4,2,1}` on the one exact integer direction d, changing only the
45 region-0 glob `dist[0]` coordinates of the frozen VXXZ24
`K100_2.37155181` vector under the transcribed 3-region single-p_comp program
at `max_level=3`, `q=5`. It did not search another direction, negative or
intermediate steps, a box, another parameter block/region, multi-block or
cross-block classes, another q/max-level regime, unpublished record
parameters, or any other construction. No paper or record claim is made.

Executed with every thread pool pinned to one and `nice -n 10`; script elapsed
`7.93340802192688` s and wall time was `8.38` s
(**COMPUTATIONAL-EVIDENCE**). Full branch/candidate intervals are frozen in
`results.json`.
