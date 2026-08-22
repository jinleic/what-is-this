# Next actions

**Snapshot:** 2026-08-21. Wave 17 CLOSED and audited. Wave 18 LANDED (eight fronts, §3). Wave 19 IN
FLIGHT (three fronts, §2). This file is the single source of truth for the ranked open queue, the
ledger ID allocation, and the environment state. Strategic reasoning lives in
`notes/frontier_analysis.md`; do not duplicate one into the other.

**Repository:** `/Users/jinleic/jinleic-workspace/math/ising3d` · **Interpreter:** `.venv/bin/python`
· `K_c = 0.221654626` stays benchmark-only · the ledger is append-only · `.venv` is shared and
immutable for agents.

## 1. Ledger and numbering state (read before appending anything)

- **LEDGER next free ID: `H506`.** Allocation to date: `H001`-`H427` historical; `H428`-`H469`
  wave-18 fronts; `H470`-`H479` lead block (unused, reserved); `H480`-`H487` wave-18
  `CliffordGradeClass`; `H488`-`H505` reserved by the three wave-19 fronts in §2 and possibly not yet
  written.
- **Next free experiment number: `e191`.** `e157`-`e181` are wave 18; `e182`-`e190` are reserved by
  wave 19.
- **Width anomalies (historical, do not edit):** `H388` has 12 fields; `H389`/`H390` have 6;
  `H393`-`H397` and `H399`-`H401` have 10. Recorded by append-only row `H402`. `H383`/`H385` are
  DUPLICATEs of `H380`/`H382`; `H384`/`H390` are corrections, per the canonical mapping in `H386`.
  `H360` predates the schema correction and is superseded by `H361`.
- **Field-convention drift, wave 17 only:** rows `H417`-`H427` are 11 fields wide but use
  `id, date, statement, topic, owner, status, artifact, artifact, detail, status-word, next-action`
  rather than the legacy header semantics. Wave-18 rows return to the legacy header.

## 2. Wave 19 in flight

| front | question | ledger | experiments | owned note |
|---|---|---|---|---|
| `TwoGenAllSize` | all-size lower bound on the TWO-GENERATOR algebra `<A,B>`, via `ad_A` grade separation | `H488`-`H493` | `e182`-`e184` | `proofs/twogen_allsize.md` |
| `TrichotomyExtend` | drop the Hamiltonian-path hypothesis; settle the non-bipartite branch; harvest `4x4`, `3x4`, `3x5`, 3D-box corollaries | `H494`-`H499` | `e185`-`e187` | `proofs/trichotomy_extend.md` |
| `LayerGroupSpectral` | isotropic all-size SPECTRAL no-go by the group route (the non-localizing route wave 18 proved is required) | `H500`-`H505` | `e188`-`e190` | `proofs/layer_group_spectral.md` |

`LadderW8` (`H458`-`H463`, `e172`-`e174`, `proofs/ladder_w8.md`) may still be running from wave 18.

## 3. Wave 18 landed — do not re-open these

- **All-`L` classification (headline).** For EVERY `L>=2` the open `2xL` local-term algebra is exactly
  `so_m(Q)+so_m(Q)` for even `L` and `sp_m(Q)+sp_m(Q)` for odd `L`, `m=2^(n-1)`, `n=2L`, with
  `dim = 2^(n-1)(2^(n-1)-(-1)^L)`. The saturation step left open since wave 16 is closed by an
  explicit production induction with invariant `I(L)`, a square-triple second move, a no-rank-loss
  transport lemma, and leaf-induction token connectivity. Full closure enumerated at `L=2..6`
  (`56, 1056, 16256, 262656, 4192256`, zero set-law violations, 19 process seconds and 437 MB at
  `L=6`); `L=7` is certified by construction plus an independent Clifford census as `67117056` and is
  explicitly NOT a full enumeration (`proofs/alll_saturation.md`, `H428`-`H433`).
- **Exact trichotomy (headline).** For any bipartite `Gamma` on `n` vertices with a Hamiltonian path,
  `N=2n`: the path gives grades `{2}` and `dim = n(2n-1)`; an even cycle gives `{2,N-2}` and
  `dim = 2n(2n-1)`; and `Delta(Gamma)>=3` gives every grade `= 2 mod 4` except `k=N`, with
  `dim = 2^(2n-2) - (-1)^(n/2) 2^(n-1)` for even `n` and `2^(2n-2) - 1` for odd `n`. In the branching
  branch the answer depends only on `n`. `C_6` (dim `132`) is the exact counterexample that fixes the
  hypotheses; `C_4 = 2x2` is in the cycle branch but accidentally has the full class. Non-bipartite
  `C_5`+chord gives `510`, the full even-grade algebra (`proofs/clifford_grade_classification.md`,
  `H480`-`H487`).
- **All-`L` quadratic no-go.** `dim g_L >= n(2n-1)+1` for every `L>=2`, so the local terms are never
  simultaneously Majorana-quadratic on the physical space; exact conservative mode floors
  `6, 18, 65, 216, 991, 3168, 15354` at `L=2..8`. Covers `L=2`, where the claw argument does not apply
  (`proofs/alll_quadratic_nogo.md`, `H434`-`H439`).
- **Isotropic gap: split verdict.** The core half SUCCEEDS — paw `119>65`, five-site tree `417>211`,
  `2x2`+pendant `475>211`, with monic upgrade and exceptional-degree bounds `913920`, `17677440`,
  `19641600`; controls `K_{1,3}` `65=65` and `K_{1,4}` `147<211` explain why the bare claw never
  sufficed. The decoupling half provably FAILS on the isotropic curve: crossing-factor minor
  `q(t)^(2c)-1`, and `r(K_{1,3} + m K_1) <= 136(2m+1) < 3^(m+4)-2^(m+4)`. Any isotropic all-size
  spectral proof must be non-localizing (`proofs/isotropic_allsize.md`, `H440`-`H445`).
- **Theorem U2 restored.** Exact char-0 `r(T)=417` and `r_all(T)=445` at four named physical points
  (derivative-gcd degrees `79` and `83`); the exact core upper bounds squeeze the modular lower bounds
  to certify `3893`, `4005`, `35597`, so hypothesis CZ is discharged and Theorem F2 is unconditional
  for `m=2,4`. `INJ(n-5)` for other chain lengths stays open
  (`proofs/rank_char0.md`, `H446`-`H451`).
- **Kac-Ward: no orbit reduction, plus a complete unit census.** The sound section group on the
  diagonal-normalised thin chart is trivial, so there are exactly `16384` singleton orbits and no
  quotient reduction exists. The census was then completed anyway by exact character evaluation over
  all `4^13 = 67108864` unit points: `2960` construction zeros in `2437` sections, each with an exact
  disposition (`2232` order-8 mod 5, `664` order-12 mod 5, `45` nonsingular mod 625, `19` singular
  lift obstructions). Classification `0 / 2437 / 13947` for
  `EMPTY_OVER_Q / HAS_LOCAL_POINT_BUT_FAILS_HOLDOUT / UNDECIDED`. Branch `11111` remains OPEN over
  `Q`: higher 5-adic digits, non-unit valuations, zero reductions, chart-degenerate loci and other
  primes are untouched (`proofs/kw_orbit_quotient.md`, `H452`-`H457`).
- **Upper endpoint: two classes exactly refuted.** At `v=6/25` the aggregate one-moment
  paired-momentum SDP has exact optimum `sup s = 1` and uniform order floor `eta* = 0`, so the
  `eta = 1/1024` challenge is infeasible with separating gap `1/1024`; the single-edge conditional
  domination random-current LP has exact optimum `p* = v^2 = 36/625` against the self-proved planar
  threshold `5/6`, short by `2909/3750`, and `p* < 1/16` uniformly. No improvement to the incumbent
  upper endpoint. Mode-resolved four-point/DLR, multi-edge or block current, and sourced-current
  classes are untouched (`proofs/upper_endpoint4.md`, `H464`-`H469`).

## 4. Ranked queue after wave 19

1. **The `W_L` law.** `W_8 = 1822` (wave 15) and `W_9 = 8586` (wave 16, `H409`) are already certified;
   the sequence is `0, 2, 10, 66, 364, 1822, 8586` for `L=3..9` with
   `K_L = (2^(2L-1)+3*2^L)/4`. Another rung is not decisive; a PROVED closed form, recurrence, or
   one-sided bound is. Note the direction discipline: an UPPER bound on `W_L` strengthens
   `dim >= K_L - W_L`, a LOWER bound only strengthens the cyclicity refutation. Seven terms do not
   determine a law.
2. **`s<=6` extensive charges — needs a NEW idea.** `s<=5` in the `Z^3` radius-2 box is decided
   (quotient exactly `1`). BOTH compression routes are closed negatively: octahedral/orbit by a cap
   lemma (`H397`, 46.93x compression reaching only 2.63% of required pivots) and the second route by
   rank non-additivity over orbit blocks (`H406`). Do not re-attempt symmetry compression.
3. **HT `v^30`, LT `x^60`, LT `x^62`.** Resource-blocked, not mathematically blocked. `x^60` is
   profile-verified launchable (123 classes, 36 CRT boxes, peak 49.12 GB); `x^62` needs the canonical
   `(5,6,6)` box, a 30-spin cross-section, i.e. the next transfer-guard lift; `v^30` needs the
   cut-capped parity frontier generalised, since the dense FLM route stalls at `v^24`. Launch only
   when the host clears.
4. **`3x4` pair-product spectral certificate.** The pipeline that cleared `2x5` and `2x6` is the tool;
   `3x4` is dim 4096 and needs explicit per-stage `process_time` budget arithmetic first.
5. **Kac-Ward branch `11111` beyond the unit census.** Remaining: higher 5-adic digits, non-unit
   valuations, zero reductions, chart-degenerate loci, `Q_p` for `p != 5`. Emptiness can never be
   certified by Nullstellensatz at any degree, and no orbit quotient exists, so the decision must be
   structural or by exhaustive holdout in a new chart.
6. **Interlayer `c10`, and a verified closed form for `c8`.**
7. **Lee-Yang edge scaling.** Provably non-identifiable at accessible sizes; needs a better estimator
   or a rigorous no-earlier-root argument, not larger cubes.

## 5. Environment state (measured 2026-08-21, 21:03 local)

- Host: 28 cores, 96 GB. Load average **~54**, ten interactive users, multi-day foreign jobs; free
  memory about 1.4 GB. This is worse than earlier in the session.
- Therefore: keep new work under ~2 GB peak RSS and budget with `time.process_time()`, never wall
  clock. All heavy series producers stay ENVIRONMENT-BLOCKED.
- A verifier failure that is a self-declared wall timeout is a resource artefact. Recheck in isolation
  before calling it a regression.

## 6. Method notes worth reusing

- **When an exact upper bound exists, construct a matching lower bound rather than refining the upper
  bound.** A matching pair is a theorem; a tighter upper bound is only a better estimate. This is how
  the all-`L` classification and the trichotomy were found, both in one wave, after the containment
  theorem had pinned the upper bound since wave 17.
- **Test every new dimension formula against every already-certified value before claiming novelty.**
  That check turned `3x3 = 65535` from an anomaly into the odd-`n` confirmation of the trichotomy.
- **`ad_A` preserves Clifford grade** when `A` is a sum of bilinears, so grade decomposition plus a
  Vandermonde argument on `ad_A^k(B)` separates `B` into eigencomponents inside the generated algebra.
  This is the tool wave 19 is using on the two-generator problem.
- **Monic-in-parameter upgrade.** A spectral invariant monic in `z` over `Z[parameter]` can only gain
  derivative-gcd degree under specialization, so ONE evaluated point bounds the generic value and one
  resultant confines the exceptions to a finite set.
- **Direction discipline.** Modular ranks and gcd degrees are one-sided. Substituting a bound into the
  wrong side cost this program Theorem U2, which took four waves to restore.
- **Razor margins are one-sided.** `deg gcd <= 295414 < 295415` proves "at least one above the
  ceiling", never "exactly one".
- **Content digests, not byte digests**, for cross-artifact references.

## 7. Standing hygiene

- Run every new standalone verifier from the repository root with the repository interpreter before
  editing any summary.
- Distinguish an OBSERVED wall (measured process time and RSS) from an unlaunched input-size
  preflight; label the latter `PREFLIGHT`.
- Verifiers must fail, not warn, when a row of a universal claim is uncovered.
- Scratch files `tmp_scout_simon.py`, `tmp_elim_simon.py`, `tmp_general_tm.py`, `tmp_trg_simon.py`,
  and `checkpoints/hypotheses_pre_overwrite_snapshot.csv` still await user deletion approval.
- Worktrees under `/private/tmp/ising3d-w1*-*` and their `omp/w1*` branches can be removed with
  `git worktree remove` once no comparison is needed; ask before deleting.
