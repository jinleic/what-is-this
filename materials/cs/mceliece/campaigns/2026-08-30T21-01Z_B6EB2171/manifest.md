# Campaign manifest — `2026-08-30T21-01Z_B6EB2171`

**Gate B: failure rate of Apon's `Δ_{p,q} ≠ 0` genericity condition**
(Vedenev/Apon McEliece dispute; ePrint 2026/1810 Lemma 9 + §3.6 hole.
Apon's Lemma 9 injectivity step is conditional on SOME pair p<q having
`Δ_{p,q} = f_p f_q' + f_q f_p' ≠ 0`; §3.6 states the all-zero case is NOT
covered. This campaign MEASURED the failure rate of that condition over
exact binary-Goppa instances, adversarially and by sampling.)

Agent: `MceliececelGateB`. Owner-notified before launch, status-checked
mid-run, and cleared to freeze with the three-layer split (Main, 2026-08-30
~21:00Z). Pre-registration committed BEFORE any computation: see
`pre_statement.md` (git `02166a6`, addenda `03ecafd`, `b8e8273`; base
17:35Z, addendum 1 17:52Z, addendum 2 after self-test, all pre-computation
for their respective stages).

## HEADLINE (three separate results — never merge them)

1. **EXHAUSTIVE THEOREM at (m,t) = (6,2).** Every monic irreducible
   degree-2 polynomial over F_64 was enumerated (exactly (q²−q)/2 = 2016
   polynomials; count verified independently by Main) and built at
   complete support n = 64: **0/2016 degenerate**. Population exhausted —
   not a sample, no confidence interval attached. Evidence label:
   **MACHINE-VERIFIED** (exact F_64 arithmetic; completeness = the
   enumeration plus the exact per-instance verdict instrument).
2. **SAMPLED layer: 0/250.** Ten (m,t) cells × 25 registered seeds,
   m ∈ {6..11}, smallest admissible t per m: **0 degeneracy events, 0
   guard failures, 0 build errors (after in-window repair)**. One-sided
   exact 95% Clopper–Pearson upper bound **0.01191** (Arb ball, radius
   1.2e-16; tail(upper) = 0.05 verified). Evidence label:
   **COMPUTATIONAL-EVIDENCE** over the named sampling distribution — this
   is a rate over the sample, NOT over the space.
3. **ADVERSARIAL special-G families: 0/127.** Z³+c (42 irreducible of 64),
   Z³+Z+c (21 of 64), and the whole depressed-quadratic family Z²+Z+c at
   m = 7 (64 of 128, Tr(c)=1): all NONDEGENERATES. CP95 upper bound
   **0.02331** (same certified-arithmetic route). Evidence label:
   **COMPUTATIONAL-EVIDENCE** over these named families.
4. **No forcing configuration found.** No instance of the Gate-A
   construction anywhere in the searched space has all Wronskians zero.
   Apon's §3.6 hole therefore remains **OPEN** — we did not close it; we
   measured how hard it is to hit. Neither disputant has even this much.
5. **V1 invariance verified:** support-ordering and Y-basis changes leave
   the verdict unchanged (machine-checked), confirming ADDENDUM A1.1's
   algebraic reduction: at full support the degeneracy verdict is a pure
   function of (m, t, G).

## LAYERS (exactly as registered — no merging)

| layer | what it is | N | degenerate | bound/label |
|---|---|---|---|---|
| E1 (theorem) | exhaustive all irreducible G at (6,64,2), full support | 2016 | 0 | none — MACHINE-VERIFIED population exhaustion |
| A (sampled) | seeds 1000..1024/cell; cells (6,2),(6,3),(7,2),(7,3),(8,2),(8,3),(9,2),(9,3),(10,3),(11,3) | 250 | 0 | CP95 ≤ 0.01191; rate over the sampling distribution |
| E2 (adversarial) | structured-G families listed above | 127 | 0 | CP95 ≤ 0.02331; rate over the named families |
| V1 | invariance verification (3 perms + 2 basis ops) | 5/5 agree | — | MACHINE-VERIFIED invariance on these probes |

Also settled during registration: **t = 1 at full support is
construct-impossible** (the unique root of a degree-1 G lies in E = full
support), and **k = 1 is unreachable for m ≤ 11** (2^m ≢ 1 mod m).

## Instrument (rule-14 discharge)

Verdicts come from a NEW exact engine (`code/engine.py`): a per-point
rank-scan of the k×2 matrix [f(aᵢ) | f′(aᵢ)] reconstructed exactly from
(Y, λ) via the interpolation identity λᵢ·fⱼ(aᵢ) = Y[j,i], with the
COMPLETENESS LEMMA that at full support the single x = Z²-substituted
polynomial Δ̃ = f_e·g_o + f_o·g_e has deg_x ≤ D−1 < n = 2^m, so
"all n pointwise minors vanish" ⟺ "all pair Wronskians are the zero
polynomial" (Frobenius x ↦ x² injective). DEGENERATE verdicts are never
produced by the scan alone in the campaign (zero events); the
counterfactual plants in the self-test prove the engine DOES return
DEGENERATE on planted all-Wronskian-zero families (route all_squares and
pointscan_complete), i.e. it is capable of establishing both directions.

Three exact instruments cross-checked on **3000 property cases + 5
counterfactual plants, 0 disagreements** (`code/selftest.py`): engine
rank-scan, brute-Wronskian polynomials (gcd/arithmetic over E[Z] via the
frozen `gfield.py`), and known-answer cases. The self-test **caught a real
iff-bug** in the first cascade draft (gcd + coprime-squares is only
SUFFICIENT for degeneracy; the general coprime-degenerate condition is
(f_e,f_o) ∥ (g_e,g_o) over E(x)) — fixed BEFORE any campaign number
existed; ADDENDUM 2 records it. A v2 synthetic-harness artifact (Y
inconsistent with its own F) was also found and repaired the same way;
no number from any self-test enters the rate.

## Inherited quantities — disclosure (owner's mid-run directive)

* The frozen 13-instance Gate-A ladder was **rebuilt byte-exactly from
  seeds** before use (`reproduction_check.json`): G match 13/13, α-guard
  agrees, β-pair brute-reconfirmed, cascade/engine verdicts agree.
* The Gate-B verdict path re-derives `λ_i F(a_i) = Y[j,i]` (frozen
  ADDENDUM 2 Lagrange legality: odd-degree coefficient of fⱼ is
  Y[j,i]·G(aᵢ)², exactly the value the engine reads); **the Ψ = C·J
  factorization is NOT on this path** (it is the waterfall-census
  factorization) — stated because it was flagged as an inherited quantity.
* Support-ordering and Y-basis invariance were PROVEN (A1.1) AND
  machine-checked (V1, T5): degeneracy is a function of (m, t, G) at
  full support.

## Incident log (inline, per rule 5)

1. `(m,t) = (6,2)` was first registered as a random-seed cell; promoted to
   EXHAUSTIVE enumeration in ADDENDUM 1 BEFORE that cell ran. The
   registered 25 seeds for it became redundant and the exhaustive run
   supersedes them; per-cell records for (6,2) seeds live in
   `layerA_state.json` and are NOT part of the headline seeds-only
   denominator (m6_t2 shows 25 = exhaustive superset).
2. A mid-run edit to `engine.py` (restoring the dropped `V = ...` line,
   which my own earlier edit had clobbered) raced the running Layer-A
   process: 9 instances at m11 recorded BUILD_ERROR from the stale module.
   They were REBUILT after the fix and all verdict NONDEGENERATE (log
   `rebuild_of: BUILD_ERROR_batch`); state dedup: seed 1018 appeared twice
   (rebuild + fresh restart), one copy removed, verdicts agreed.
3. One edit accidentally landed on the FROZEN campaign's `instance.py`
   (hash-matched resolve); **reverted immediately, sha256 re-verified
   byte-identical** (5248caf7…). No frozen artifact changed.
4. Main's cross-check values (0.01305 @ N=228 pre-dedup, 0.02331 @ 127)
   reproduce exactly; the final N after dedup is 250 → 0.01191.

## Files

* `pre_statement.md` — pre-registration + both addenda (committed pre-run).
* `code/engine.py` — exact verdict instrument; `code/selftest.py` — the
  3000-case + 5-plant battery; `code/run_layerA.py` / `code/run_layerB.py`
  — runners; `code/cp_bounds.py` — Certified-arithmetic CP bounds (Arb
  balls + mpmath cross-check); `code/gfield.py`, `code/fastfield.py`,
  `code/census.py`, `code/instance.py` — byte-copies of the frozen
  57200ADD sources (sha256 in `code_hashes_imported.sha256`,
  scratch copy) with the SINGLE gate-B edit (`forced_G`, default-path
  byte-equivalence machine-checked both directions on instance 1387).
* `layerA_state.json` / `layerA_log.jsonl` — per-instance records
  (checkpointed after EVERY instance; survivors of the mid-run timeline
  all present). `layerB_state/log` — exhaustive + adversarial records.
  `reproduction_check.json` — frozen-ladder reproduction records.
  `cp_bounds_arb.json` — Arb-ball intervals + mpmath cross-checks.
  `final_tally.json` — the numbers in this manifest, derived from state.
* `sha256s.txt` — checksums of everything in this directory.

## Run environment

macOS (Darwin 25.5.0) Apple M3 Ultra; CPython 3.14.3; python-flint 0.9.0
(Arb); numpy 2.5.2; `nice -n 10`, threads pinned to 1. All verdicts EXACT
finite-field arithmetic (F_{2^m} tables + E[Z] poly ops); NO floats
anywhere in any verdict; the only floating numbers in this campaign are
the displayed CP bounds, each backed by a printed ball enclosure.

— MceliececelGateB, 2026-08-30T21:0xZ.
