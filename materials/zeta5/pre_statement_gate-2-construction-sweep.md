# Preregistration — gate-2-construction-sweep (math/zeta5)

Frozen before compute, 2026-09-08, institute cycle 20260908T150309Z_2b7759.
This file is both the hmz FREEZE artifact and the campaign prereg. One run is
funded (state.json budget runs_used 0/2).

## Gate question

Does any construction in a small, explicitly bounded family of Apéry-style
candidate linear forms for $\zeta(5)$ survive the README gate-2 retention
criteria with a certified exponential error bound strictly better than the
gate-1 baseline (Zudilin recursion: $\mu_2 = 0.33753726443403620704$,
$\log|\mu_2| = -1.0860793616267493417$, honesty deficit $3.9139206383732506583$
against the $2D_n^5$ inclusions (6))?

## Frozen family bounds (the entire swept ladder; nothing outside it)

Exactly the ladder enumerated by `src/gate2_sweep.py::build_ladder`, 44 probes:

- **F1 — bounded deformations of the very-well-poised series (7) of
  arXiv:math/0206178.** Slope pair $(s_1,s_2)\in\{-1,0,1\}^2$ on the two
  numerator Pochhammer products, denominator exponent $w\in\{4,5,6\}$.
  27 probes. Terms exactly rational (sympy `Rational` end-to-end); partial
  sums over $k=1..40$; alternating weight $(-1)^k(1/2)^k$.
- **F2 — Vasilyev-type binomial family.** $t\in\{1,2,3\}$,
  $s\in\{0,1,2,3,4\}$: $\sum_{k=1}^{30}\binom{n}{k}^t\binom{n+t}{k}^s/k^5$
  (exactly rational). 15 probes.
- **F3 — Ball–Rivoal-type odd-weight family.** $r\in\{1,2\}$:
  $\sum_{k=1}^{30}\big(\binom{2k}{k}\binom{n}{k}/\binom{n+k}{k}\big)^r/(2k+1)^5$
  (exactly rational). 2 probes.

Every probe: 80 exact rational terms $a_0..a_{79}$; recurrence guessing with
`src/recsearch.py::guess_recurrence` at orders 2 and 3, polynomial-coefficient
degree $\le 14$, `start=1`, `extra=8` — FLINT exact nullspace over $\mathbb Q$
with full-range exact re-verification of every candidate (no PSLQ anywhere in
the guess path). Rate screening: Poincaré–Perron constant-coefficient surrogate
(min $|root|$ at $n_{\rm ref}=40$), plus the multi-anchor stability audit at
$n_{\rm ref}\in\{10,30,40,79\}$ for every promoted probe.

## Falsifiable acceptance criteria (retention, per README gate 2)

A probe is **retained** only if ALL THREE hold:

- (a) exact symbolic recurrence identity (guessed AND exactly verified on the
  full 80-term range; symbolic certification is the gate-3 layer — a
  full-range-exact guess satisfies (a) at the sweep level and is labelled
  CONDITIONAL per trust discipline);
- (b) integrality / denominator control with linear-form content in
  $\{1,\zeta(3),\zeta(5)\}$ — a degenerate sequence whose numerator product is
  constant in $k$ (shrinking rational, no form content) FAILS (b) trivially;
- (c) certified exponential error bounds strictly better than the baseline:
  surrogate min $|root| < \mu_2 = 0.33753726443403620704$ **at every audit
  anchor** $n_{\rm ref}\in\{10,30,40,79\}$; a rate that crosses 1 or $\mu_2$
  at any anchor is an evaluation artifact and the promotion is withdrawn.

**Retention count = 0** means zero survivors better than baseline — the
expected-plausible outcome per the README ("the campaign fails here; record
and stop"). No survivor will be manufactured: rate promotions that fail the
audit are recorded as withdrawn, not retained.

## Verdict mapping (binding)

- **FROZEN-CERTIFIED** — at least one probe satisfies (a)+(b)+(c) with all
  anchors below baseline.
- **FROZEN-NEGATIVE** — all 44 probes are evaluated through the guess + audit
  pipeline and retention count is 0.
- **FROZEN-INCONCLUSIVE** — the tooling fails, the ladder cannot be evaluated
  as frozen, the record diverges semantically from `data/GATE2_SWEEP.json`
  (different probe count, different verdict, or a spot-checked probe's stored
  candidate not reproducing exactly), or scope must be exceeded to answer.

## Scope boundary

- In scope: the 44-probe ladder above; re-execution of the existing tooling
  (`zudilin_rec.py` baseline re-certification, `recsearch.py --selftest`
  guesser validation, `gate2_sweep.py` sweep + audit); exact spot-check
  reproduction of recorded probes; comparison against the pre-existing record
  `data/GATE2_SWEEP.json` (sha256
  `d7d71318fa62d2c3276eff473409b0ad8f68c968bf2da6673c46e53fb8476cb7`,
  verified intact at freeze time).
- Out of scope: any new family, parameter, order/degree widening, or a fresh
  ladder; creative-telescoping symbolic certification (gate-3 layer, spent
  only on a survivor); edits to `src/` or `data/` other than the sweep tool's
  own deterministic rewrite of `data/GATE2_SWEEP.json`; any attack on gates
  other than gate 2.

## Exact verification command

```sh
cd /Users/jinleic/jinleic-workspace/math && \
  ./.venv/bin/python -B zeta5/src/gate2_sweep.py; \
  echo "exit=$?"; \
  shasum -a 256 zeta5/data/GATE2_SWEEP.json
```

Pass: exit 0; stdout `VERDICT: ZERO_RETAINED ...` with
`survivors_stability_audited == []`; record sha256 reproduces
`d7d71318fa62d2c3276eff473409b0ad8f68c968bf2da6673c46e53fb8476cb7`
byte-identically (a byte-level divergence with identical semantic content —
same 44 probes, same verdicts, same retention outcome — is documented in
AUDIT.md and does not by itself change the verdict).

Supporting commands (must also exit 0):

```sh
cd /Users/jinleic/jinleic-workspace/math && ./.venv/bin/python -B zeta5/src/zudilin_rec.py
cd /Users/jinleic/jinleic-workspace/math && ./.venv/bin/python -B zeta5/src/recsearch.py --selftest
```
