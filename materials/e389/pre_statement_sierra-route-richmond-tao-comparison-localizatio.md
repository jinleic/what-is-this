# Preregistration — Sierra route: Richmond/Tao comparison, or accept the conditional frontier

Target: `math/e389` (Erdős Problem #389: consecutive-product divisibility).
Gate: `sierra-route-richmond-tao-comparison-localizatio`.
Cycle: `20260908T150309Z_2b7759`. Date frozen: 2026-09-08 (UTC).

This file is frozen before any campaign compute. It is both the hmz freeze
artifact and the campaign preregistration.

## 1. Gate question

state.json `next_action` asks: pursue the **Sierra route** — the effective
per-class smooth-count bound route named by the Wright indicator-corollary
repair discussion (`LOCALIZATION.md` §3.1 and §10 item 3), labelled
"Richmond/Tao-style comparison theorem route" in
`math/PROGRESS.md` (E389-SHIU-EFFECTIVE-MARGIN, next-action note) — **or accept
the conditional frontier** as the honest record for now?

The deliverable is a **route decision memo**, not a new theorem. It must:

1. compare the remaining live routes on record — at minimum
   (i) the effective per-class bound route (§3.1/§10 item 3),
   (ii) the Weyl-sum route of Input III (§4, §10 item 4, (46)), and
   (iii) the Input VII paired moving-center estimate (§8, §10 item 1, (94)),
   and for completeness
   (iv) the Input-I citation-complete frontier with its two open gaps
   (§2, §10 item 5);
2. assess **what the Sierra-route Richmond/Tao comparison would buy**, including
   an explicit attempt to identify the named theorem(s) in the literature,
   with the search trail preserved;
3. conclude exactly one of:
   (a) a specific next mathematical action, with justification; or
   (b) that the conditional frontier should be accepted for now.

## 2. Falsifiable acceptance criteria

The memo (`campaigns/<run_id>/DECISION_MEMO.md`) is accepted iff all hold:

- **A1 (route coverage).** Each route (i)–(iv) is named, and for each the memo
  states the target statement, the status of its instrument on record, and the
  exact obstruction, with `LOCALIZATION.md` section citations.
- **A2 (metered requirement).** For route (i) the memo states the certified
  requirement: an effective per-class bound of shape $\Psi\ll c\,y\,\rho(2)^{e}$
  with $e>e_{\rm req}=0.880084778655420\ldots$ at $M=\sqrt{2N}$; the certified
  verdicts on the two audited mechanisms (§3.1 Wright REFUTED-for-purpose,
  break-even $u^\ast\approx2.2125>2$; §3.2 Shiu $R\ge2(M/\phi)\ge2$ at every
  depth, coverage stop $M=2\sqrt N$); and the per-class sieve-level form of the
  gate, $\ln z' > M/\phi$ against $\log(y/M)/20$ available.
- **A3 (identification trail).** The memo records the attempt to identify the
  "Richmond/Tao comparison" theorem, with the concrete trail (arXiv API author
  and full-text queries, Semantic Scholar, web search, grep of the audited
  Wright source `shiu_for_journal.tex`, erdosproblems.com/389, OEIS A375071,
  Tao's recent math.NT bibliography) and its outcome. A3 does not prescribe
  the outcome, only that the trail exists and the outcome is stated.
- **A4 (decision).** Exactly one of (a)/(b) is concluded; the adopted next
  bounded objective is written into the memo and, at close, into
  `state.json.next_action`. `current_gate` is updated only if the gate actually
  moved.
- **A5 (scope).** No new theorem proving; no sieve estimates beyond quoting the
  certified numbers already on record; no edits outside `math/e389/`; the
  domain ledgers receive only the proposed entries in
  `campaigns/<run_id>/LEDGER_ENTRIES.md` (driver appends serially at settle).

## 3. Verdict mapping

- **FROZEN-CERTIFIED**: a decision (a) or (b) is reached and recorded, with
  A1–A5 verified fresh by the verification command in §5.
- **FROZEN-INCONCLUSIVE**: the routes were compared but the record does not
  support concluding either (a) or (b) within the bound (the decision itself
  is undecidable from the record, not merely unpalatable).
- **FROZEN-NEGATIVE**: the record contradicts a premise of the gate question
  (e.g. the surviving route §10 item 3 names does not exist as recorded, or a
  cited certified verdict does not appear in `LOCALIZATION.md`).

## 4. Scope boundary

Analysis/decision work only. Explicitly excluded: proving or disproving (45),
(46), (94); closing the two Input-I gaps; any new sieve computation; editing
`LOCALIZATION.md`, `THEOREMS.md`, `README.md` or any file outside
`math/e389/`; editing `math/RESULTS.md`/`math/PROGRESS.md` directly (LEDGER_ENTRIES.md
only); touching any other target, cycle, or closed/frozen run dir; hand-editing
run ids, claim files, `manifest.json`, `status.json`.

## 5. Exact verification command

After the memo is written, from the workspace root:

```sh
python3 math/e389/campaigns/<run_id>/verify_decision.py
```

which must exit 0. It checks, fresh from disk:

1. `DECISION_MEMO.md` exists and contains every required anchor: the four
   route names; the certified constants `0.880084778655420`, `2.2125`,
   `17/32`, `165/328`; the decision keyword; and the adopted next action;
2. the cited record anchors exist in `math/e389/LOCALIZATION.md`
   (the §3.1 REFUTED-for-purpose verdict, §3.2 verdict, §10 items 1–5);
3. the identification-trail artifact exists and records the searched channels;
4. `state.json` is valid JSON and its `next_action` equals the adopted action
   recorded in the memo (checked at audit time, after step 9 of the contract).
