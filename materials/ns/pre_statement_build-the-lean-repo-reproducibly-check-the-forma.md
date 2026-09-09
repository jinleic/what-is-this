# Preregistration — math/ns — gate `build-the-lean-repo-reproducibly-check-the-forma`

**Cycle:** 20260908T193320Z_ba7b32 · **Target:** math/ns · **Owner label:** NsRunner
**Frozen at:** 2026-09-08T19:37Z (before any evidence is gathered)

## Gate question

OpenAI claims (2026-09-08, post + PDFs + public Lean repo
`openai/NavierStokesAndEuler`) that Clay Navier–Stokes alternatives **(C) and
(D) are resolved**, with a public Lean formalization. Independent verification
question, decomposed:

1. **Statement match.** Does the formalized main theorem match the CMI
   statement (C) (Fefferman) — solution class satisfying (1)(2)(3), energy
   conditions (6)(7), domain $\mathbb{R}^3$ (not a torus/box), Schwartz-decay
   initial data (4), force regularity class (5) — and the claimed (D), or is
   something weaker/different proved?
2. **Proof integrity.** Is the Lean development free of load-bearing
   `sorry`/`admit`, and does `#print axioms` on the main theorems show only
   the standard Lean/Mathlib axioms (no nonstandard axioms, no load-bearing
   `native_decide`)?
3. **Reproducible build.** Does the pinned repo build from scratch on
   mini-pro.local with the repo-pinned elan toolchain (`lake exe cache get`
   then `lake build`)?

## Falsifiable acceptance criteria

- **A1 (statement):** the Lean statement is quoted verbatim into the run dir
  and each CMI clause — (1)(2)(3) solution class, (4) Schwartz data, (5) force
  class, (6)(7) energy/global conditions, $\mathbb{R}^3$ vs torus, and the
  (D) analogue — is marked MATCH / WEAKENED / DIFFERENT, with the precise
  divergence named for every non-MATCH.
- **A2 (integrity):** `grep -R "sorry\|admit\|axiom"` over the Lean source is
  reported in full; `#print axioms` output for every main theorem is captured;
  any `native_decide` site is listed and judged load-bearing or not.
- **A3 (build):** the commit hash is pinned in the run dir; the toolchain
  (lean-toolchain file) is recorded; build logs are copied into the run dir;
  exit status of `lake build` is recorded.

## Verdict mapping (terminal verdicts of campaign.py)

- **FROZEN-CERTIFIED** — build passes (A3 green) **AND** statement matches CMI
  (C) verbatim on every clause (A1 all MATCH) **AND** no load-bearing sorry
  and no nonstandard/auditable-gap axiom (A2 green).
- **FROZEN-NEGATIVE** — the claim fails on the merits: statement is
  weakened/different on any clause (solution class, force class, energy
  inequality direction, domain, data class), **OR** a load-bearing
  sorry/admit/nonstandard axiom is found, **OR** the repo/paper does not exist
  as claimed. The exact divergence is named in the verdict line.
- **FROZEN-INCONCLUSIVE** — build/infrastructure prevents a conclusion
  (toolchain unavailable, build cannot finish inside the window, network
  failures, RAM exhaustion) with no merits verdict reachable.

## Scope boundary

- Verification only. No attempt to re-derive the mathematics, no refereeing of
  the paper's PDE content beyond what the Lean statement and build settle.
- The workspace obstruction map (README: Lemma 1 Liouville step, DSS
  reductions, Lemma 4 Leray covariance) is confronted in UNDERSTAND.md as
  context (route through/around/contradict), but that confrontation is **not**
  part of the verdict mapping — the verdict rests on A1–A3.
- One run only this cycle. Wall ceiling 4 h total with sub-gates:
  statement match ≤1 h (cheapest-first, before any long build);
  sorry/axiom audit ≤30 min; build ≤2.5 h (detached on mini-pro; if it cannot
  finish, STOP and close FROZEN-INCONCLUSIVE with the log);
  UNDERSTAND.md ≤30 min.
- Host: mini-pro.local only for remote compute. Control plane strictly local.

## Exact verification commands

```sh
# on mini-pro.local, repo cloned at pinned commit:
cd ~/ns-verify/NavierStokesAndEuler
git rev-parse HEAD                          # pinned into run dir
cat lean-toolchain                          # recorded into run dir
grep -RInE 'sorry|admit|axiom' --include='*.lean' .   # A2
lake exe cache get && lake build            # A3 (detached, nohup caffeinate -s)
# A2 (axioms): lake env lean on a scratch file importing the main theorems
#   with #print axioms <MainTheorem>
```

Locally: `curl -sL` the post + both PDFs; `pdftotext` the main paper; extract
the CMI-relevant statement and compare against the Lean statement verbatim.
