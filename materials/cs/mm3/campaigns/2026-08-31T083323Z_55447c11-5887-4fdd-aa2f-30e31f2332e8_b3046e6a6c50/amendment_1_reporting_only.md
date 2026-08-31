# Amendment 1 — reporting semantics only (before compute)

**Committed before any target computation in this campaign. Domain, action, order, caps, arithmetic, and escalation policy are unchanged.**

The producer's unused `render_report()` headline is retracted before use. It says that `paper55` would be the unique 55 *lower-bound* minimizer if both new decompositions have minimum at least 56. That wording is false because the completed campaign already has two lower-bound-55 data triples: the all-monomial triples of `paper55` and `sun56`. Only `paper55` carries a verified 55-addition upper-bound circuit; `sun56` remains certified lower bound 55 / published upper bound 56.

Therefore this campaign will run only the byte-hashed `--run` computation. It will not invoke that source's `--freeze` reporting branch. Freeze will byte-copy the exact as-run source and dependencies, write a report that distinguishes (a) certified lower bounds, (b) known upper-bound witnesses, and (c) unresolved gaps, then generate `manifest.json` and `checksums.sha256`. No numerical result may be altered by this reporting correction.

Correct conditional wording:

* if each new minimum is at least 56: **no new lower-bound-55 data triple enters from `mws59` or `stapleton60`; `paper55` remains the only verified 55-addition circuit in the expanded swept set, while `sun56` remains LB 55 / UB 56**;
* if a new minimum is 55: **a new lower-bound-55 data triple exists, but it is not a 55-addition circuit unless an explicit upper-bound witness is synthesized and verified**;
* any value at most 54 still triggers the pre-registered halt and Main escalation before it is recorded.
