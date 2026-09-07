# N5 pre-audit amendment — independent verifier budget boundaries

Run: `20260904T044841Z_861bb779_61abb80d3b13`. Date: 2026-09-04.
Preregistration SHA-256:
`187e2c8d244df75a20bd83e3119a8c622ef3f5f2940e54d535acfc1f715e85ee`.

The authoritative `n5_instrument.py` completed its single registered pass and
emitted witness SHA-256
`6e3bd5ae1f00c97d089975720e5907c0e81474e913eed38bca04b4c9ebd45b2f`.
A subsequent static review found one non-arithmetic protocol defect in the
separate independent verifier: it installed equal soft/hard 7200-second
`RLIMIT_CPU` values but did not call a one-second-margin budget check at its
phase boundaries, so a cap hit could terminate without its own orderly
budget-abort JSON. No Hamilton-table, witness, C/Cinv, transport, substitution,
or corruption-control defect was found.

The first independent audit is retained as non-adjudicative history:

* pre-fix source `independent_resubstitute_pre_budget_fix.py` SHA-256
  `a5b2313f594f34790d09c4752e81969c9e071996a46b70824485976017863775`;
* pre-fix output `n5_independent_audit_pre_budget_fix.json` SHA-256
  `aeaa124f435c503b7c312794bed65c1e7b135ca64df3a25c3def8973606801ad`.

Before the adjudicative independent rerun, the verifier was changed only to
add `budget_check(stage_name)` with the registered one-second safety margin,
an orderly `n5_independent_budget_abort.json`, and calls at runtime, witness
load/parse, independent target/substitution, and final serialization
boundaries. Fixed verifier SHA-256:
`dbce8c0378d6aa473fe479746086525f28d3ff20c89eea734082a64a153429ed`.
No scientific knob, coefficient, target entry, witness byte, evidence label,
or verdict rule changed.

Adjudication rule: rerun only the fixed independent verifier with
`PYTHONDONTWRITEBYTECODE=1 nice -n 10`; certification requires its fresh
64/64 exact result and singleton corrupted-target rejection. A failure or
budget abort forces FROZEN-INCONCLUSIVE. The preserved pre-fix pass carries no
weight.
