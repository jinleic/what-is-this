# Magic-cat stabilizer-rank witness

## Current state (Main, 2026-09-06)

**C for unchanged local simulated annealing; mathematical rank remains open.**
The [registered assay](campaigns/20260906T025240Z_725b49e8_31de21aa6d69/)
closed **FROZEN-INCONCLUSIVE**. The public search recovered exact
`cat_4`/rank-2, `cat_6`/rank-3 and `cat_8`/rank-6 witnesses. Seven rank-5
seeds completed without a witness; seed 7 hit the frozen 60-second cell wall
limit. No limits, schedules or seeds were changed.

Assay plus supervision: **73.105411 CPU s / 72.049315 wall s**; largest child
RSS **100,204,544 bytes**. The independent checker reconstructs genuine
stabilizers and checks the full target exactly over Q(i). After freeze, the
six-term cat-eight certificate replayed all 256 coordinates successfully.
Corrupted phase, target and support controls were rejected.

Post-assay exact diagnostics on the already returned tuples agree with all
seven numerical residuals. The best returned five-state span has squared
distance **205/2016** from the normalized target. This is a distance to that
particular span, **not** a global floor or a rank lower bound. No new record,
generalization or AI-value claim follows. No LLM search is queued.

Replay the frozen known six-term witness from the workspace root, without
running the search or writing into its frozen directory:

```sh
physics/stabrank-witness/.venv/bin/python -I -B \
  physics/stabrank-witness/campaigns/20260906T025240Z_725b49e8_31de21aa6d69/exact_verify.py \
  physics/stabrank-witness/campaigns/20260906T025240Z_725b49e8_31de21aa6d69/cat8-r6-seed0.json \
  physics/stabrank-witness/scratch/cat8-replay.json
```

Dependencies: Python 3.12.12, SymPy 1.14.0 for exact replay; search dependencies
and all effective SA settings are pinned in the campaign registration.
The next scientifically justified step needs a material structural ansatz,
not more seeds or inference from variation among local minima.

## Claim and frontier at admission

**A for a bounded headroom assay; headroom is not yet demonstrated.** Find an
exact expression of the eight-qubit Qassim–Pashayan–Gosset magic cat state as a
linear combination of at most five stabilizer states. The strongest located
published bound is six terms; the exact optimum is unknown.

Define `T = (|0> + exp(i*pi/4)|1>)/sqrt(2)`, `T_perp = Z T`, and
`cat_m = (T^tensor(m) + T_perp^tensor(m))/sqrt(2)`. This is the edge-type orbit,
not the distinct face-type Bravyi–Kitaev T state. In computational coordinates,
the unnormalized target is zero at odd Hamming weight and `i^(weight/2)` at
even weight. Removing the parity ancilla by CNOTs gives the compressed target
on `m-1` qubits: amplitude `i^ceil(weight(y)/2)`.

A five-term witness would imply `chi(T^tensor(n)) = O(5^(n/6))` by chaining cat
states with stabilizer `cat_2` contractions: `l` copies yield `cat_(6*l+2)` with
at most `5^l` terms. Together with `chi(T^tensor(n)) <= 2 chi(cat_n)` and
projection monotonicity, the exponent is `log2(5)/6 = 0.386988...`, below the
located `log2(3)/4 = 0.396240...`. This direct cat-contraction argument does not
rely on the regular-exponential-growth assumption in the general-code Theorem 4.

Primary sources:
- [Qassim–Pashayan–Gosset, Quantum 5, 606 (2021)](https://arxiv.org/abs/2106.07740),
  magic-cat contractions and Appendix A (`chi(cat_8) <= 6`).
- [Labib–Russo (2026)](https://arxiv.org/abs/2605.28586), current related orbit
  results and exact-verification workflow.
- [Public stabrank leaderboard](https://unitaryfoundation.github.io/stabrank/)
  and [search engine](https://github.com/unitaryfoundation/stabrank).

## Ten admission gates

1. **External value:** finite witness improves an asymptotic simulation bound,
   not an internal score.
2. **Precise statement:** exact rank upper bound for the pinned magic-cat state;
   amplitudes and Clifford compression are explicit above.
3. **Strongest comparator:** public six-term cat construction and contemporary
   non-AI simulated annealing, not a naive random-state baseline.
4. **Accessible artifacts:** public formula and wheel `stabrank==0.1.1`;
   isolated local environment, package files pinned before compute.
5. **Verifier:** independent affine-support/quadratic-phase reconstruction and
   exact Gaussian-rational span membership; a float residual is never proof.
6. **Headroom:** unknown. Reproduce known-answer cells before probing rank five;
   a failed search is not a rank lower bound.
7. **Search leverage:** finite exact witness in dimension 128; nonlocal sharing
   of stabilizer states is the representation, not open-ended prose proving.
8. **Novelty:** only an exact witness below the located bound can become a
   discovery, after refreshed prior art. Reproductions are labeled separately.
9. **Bounded cost and AI role:** one chain/thread at a time, fixed seeds and
   iteration limits; no LLM search until a separate matched-budget experiment
   is justified. This assay is non-AI.
10. **Failure value:** distinguish an uncalibrated instrument from a calibrated
    search miss; preserve exact candidate tuples, residuals and resource costs.

## Frozen assay design

The machine-readable registration is `scratch/headroom-preregistration.json`;
its hash is bound by `../../scripts/campaign.py init` before any target search.
The candidate ladder is compressed `cat_4` at rank 2, `cat_6` at rank 3,
`cat_8` at rank 6, then the one discovery cell `cat_8` at rank 5. The first
three are calibration/reproduction controls. A failed calibration stops the
ladder; it does not imply absent mathematical headroom.

Use the public Pauli-expansion simulated annealer with every exposed argument
pinned, including `num_chains=1`; the upstream multi-chain example would
violate this target's one-thread policy. Sequential seeds are fixed in the
registration. No schedule tuning, retries or cell extensions after results.
An independently exact-verified witness is success. All other results remain
bounded search observations. Repeated nonzero numerical residuals are neither
exact algebraic floors nor lower bounds.

The verifier accepts only affine support `x0 + span(W)` over F2 and phase
`i^(sum l_j*y_j + 2*sum q_jk*y_j*y_k)`, reconstructs the vectors, and solves for
coefficients over Q(i). It checks the full parity-expanded QPG target, not the
library's target implementation. Canonical certificates are extracted by the
untrusted runner; the verifier does not import the search engine. A deliberately
corrupted target coefficient must be rejected.

No holdout/generalization or AI-value claim arises from this fixed-state assay.
Any later policy-evolution study needs a separate public instance panel,
untouched evaluation set, strong fixed SA and non-LLM mutation controls,
matched candidate/CPU budgets, and separately reported model calls and cost.

## Publication and stopping boundary

An exact five-term cat-eight witness would justify a research note plus an
independent formal/mathematical review. This session authorizes neither external
submission nor a large search. If the instrument fails, first isolate the exact
calibration blocker. If a calibrated bounded search finds nothing, stop the
registered lane; new structures require a new protocol. Do not claim rank six
is optimal, and do not infer LLM headroom from variation among local minima.
