# Aborted campaign report

**Verdict: ABORTED BEFORE ANCHORS; no target claim and no new count.**

The first startup-anchor loop loaded all five decompositions and then stopped before a Brent, floor, survivor, or population result because `META["perminov58"]["split"]` is intentionally `None`; Perminov's 58 is stored in the pinned JSON's `complexity.reduced` field rather than a left/right/output split. The instrument incorrectly attempted to iterate the missing split and raised `TypeError`.

This is an implementation/provenance-path defect, not evidence about any decomposition. The byte-identical failed source and verbatim run log are frozen here. Per rule 5 the campaign is retained rather than deleted or silently resumed. The fix must read and assert `cr58_cn122_ZT_reduced.json["complexity"]["reduced"] == 58`, receive a new code hash, and launch under a new pre-statement.
