# Millennium / open-problem attack repository

Single source of truth. Everything authoritative lives here; nothing in
`scratch/` is canonical.

## Tracking contract (SSOT)

One fact, one owner. Three layers, no duplication:

| layer | owner | contents |
|---|---|---|
| chronological ledger | [`PROGRESS.md`](PROGRESS.md) | dated, newest-first session entries; every claim carries its verification; retractions inline |
| current state per target | `<target>/README.md` | what is PROVED / CONDITIONAL / NUMERICAL right now, file inventory, how to run |
| run artifacts | `<target>/campaigns/<UTC-timestamp>_<uuid>_<code-hash12>/` | immutable frozen snapshots + committed results; inventoried in `<target>/campaigns/README.md` |

Naming rules:

* **One problem, one top-level folder**, short lowercase name (`uc/`, `r55/`,
  `zeta5/`, `ns/`, `ccf/`). No dated top-level folders, no second folder for
  the same problem.
* Campaign directories are producer-generated
  (`cert3_<UTC>Z_<uuid>_<hash>`); never hand-created, never edited after
  launch. Disposable rehearsals go under `<target>/campaigns-smoke/`.
* Update the ledger and the owning README **in the same session** as the work;
  a result that is not in `PROGRESS.md` does not exist.

> **Naming drift, 2026-08-15:** `ramsey-r55/` appeared beside `r55/` (fresh
> catalog downloads; `results/`/`scripts/` empty; writer session unknown).
> $R(5,5)$ work belongs in `r55/`. Merge `ramsey-r55/data` into `r55/data`
> once its writer is done; do not start new work there.


## Targets

| directory | problem | state |
|---|---|---|
| [`uc/`](uc/README.md) | **Union-closed sets conjecture** (Frankl 1979) — *active* | refuted the reduction behind the published record and replaced it: **Theorem B‴** (support $\le2$ pair-orbits, 5 parameters, complete writeup `uc/thmB3_proof.py`), Margin Lemma tight at the obstruction, exact bridge to Cambie Q2, $rh\le2x$ proved. $\psi=\frac{3-\sqrt5}{2}$ still the largest constant with a complete proof. Now attacking $t=\psi+10^{-4}$ by certified 5-D branch-and-bound: **replayable proof traces** (one byte per DFS node; independent checker re-proves every step — `uc/cert3_replay.py`); immutable 8-slice campaign running since 2026-08-15, inventory in [`uc/campaigns/README.md`](uc/campaigns/README.md) |
| [`r55/`](r55/README.md) | **Ramsey number $R(5,5)\le45$** — *active* (opened 2026-08-13) | campaign scaffolded; gate 1 (independent validation of McKay's Ramsey catalogs: full (3,5), (4,4) sets, the 352,366 Ramsey(4,5,24)-graphs, 328 known Ramsey(5,5,42)-graphs) in progress; recon on Angeltveit–McKay $\le46$ proof structure and SAT-certification prior art |
| [`zeta5/`](zeta5/README.md) | **Irrationality of $\zeta(5)$** — *active* (opened 2026-08-15) | gate 1 passed: Zudilin's third-order recursion (arXiv:math/0206178) reproduced and certified exactly (convergents, integrality (6) to $n{=}300$, signs, roots, rates); the honest baseline **deficit is $5+\log\lvert\mu_2\rvert=3.914$** against the $D_n^5$ threshold; recurrence-search tool validated by blind re-derivation of the recursion from raw $q_n$ data. Gate 2: the triage's four-week falsifiable construction sweep |
| [`kobon/`](kobon/README.md) | **Kobon triangle problem, $n=10$** (Fujimura) — *active* (opened 2026-08-13) | deciding $K_{\rm gen}(10)=25$ under the broad convention (parallels, multipoints, crossed triangles allowed — no located source proves this case; the classical charging lemmas fail over $\mathbb{Q}$, exact counterexamples in `kobon/report.md` §6). Lower bound 25: exact integer certificate, verified. Upper bound: 11-cube necessary-clause CNF cover; 9 cube DRATs re-verifying, 2 truncated proofs being regenerated (whole-cube re-solves + tautological 8-way splits) |
| [`qec/`](qec/README.md) | **Automated co-design of quantum LDPC codes, syndrome-extraction circuits, and decoders under circuit-level noise** — *active* | non-CSS PBB programme: parent-domination theorems (δ=0 proved 213/368; two certified δ>0 distance reversals); basis-independent depth-≥8 separation (Thm C6) vs Gross depth 7; end-to-end circuit-level negative for the PBB [[144,12,12]] target; target `12_6_0193` proven genuinely non-CSS under row ops + qubit permutations + arbitrary independent local Cliffords (exact GF(2) XOR contradiction certificate); independent target distance certified 6≤d≤12, four CP-SAT orbit sectors running; 21-entry failed-route ledger |
| [`h10q/`](h10q/README.md) | **Hilbert's Tenth Problem over $\mathbb{Q}$** — *active* (opened 2026-08-12, canonicalized 2026-08-15) | frontier mapped and fully executable (definability + elliptic pillars); **L6 witness-tie**: universal definition of $\mathbb{Z}$ in $\mathbb{Q}$ with **6 quantifiers, conditional** on one explicitly scoped global norm-assembly lemma — W0–W2 proved and machine-checked, 61/61 guarded canonical certificates for $w<300$; unconditional records untouched ($\forall_{10}$ refereed Daans, $\forall_7$ unrefereed Sun) |
| [`ns/`](ns/README.md) | Clay Navier–Stokes, Fefferman (C) on $\mathbb{R}^3$ — *paused* | obstruction map; blocked on two publisher-inaccessible primaries (NRS, ESS) |
| [`ccf/`](ccf/README.md) | Córdoba–Córdoba–Fontelos model — *paused* | technique validation only; prize-irrelevant |

`uc/` was chosen over `ns/` on one criterion: **every primary is on arXiv.** The
Navier–Stokes attempt stalled because NRS and ESS could not be retrieved, and an
obstruction map that rests on unread primaries cannot be pushed further
responsibly.

## Discipline

1. Read primaries first-hand; quote verbatim. Secondary summaries are not evidence.
2. Every quantitative claim is machine-checked by a script in this repo.
3. Distinguish PROVED / CONDITIONAL / NUMERICAL at every step, always.
4. Verify by hand anything an optimiser reports before believing it.
5. Record retractions inline rather than silently editing claims away.
6. A theorem is not proved because its asserts pass — asserts check arithmetic at
   the points you chose. State the measure class and the cone explicitly; when
   writing "concave"/"convex", write the inequality out and check its direction
   against the definition; carry a corollary's hypotheses into every restatement.
7. **Never promote a sampled parameter check to a universal claim.** Where the
   quantity is monotone or affine in the parameter, solve for the boundary and
   assert the equivalence. Sampling a range is evidence, never a quantifier.

---

# Navier–Stokes target (paused) — detail

## Canonical target

**Clay Millennium Problem: Navier–Stokes existence and smoothness — Fefferman
alternative (C)**, i.e. finite-time breakdown for the 3D incompressible Navier–Stokes
equations with constant viscosity $\nu>0$ on $\mathbb{R}^3$ with Schwartz-decay initial
data, with an admissible force permitted. **(D)** on the torus $\mathbb{R}^3/\mathbb{Z}^3$
is **out of scope**: every argument here uses the Euclidean dilation $x\mapsto\lambda x$,
which is not a self-map of the fixed torus, and Liouville differs there (periodic
harmonic ⟹ *constant*, not zero). (C) alone suffices for the prize.

CMI rules §5(b): *"In the case of the P versus NP Problem and the Navier-Stokes Problem,
a resolution in either direction will be evaluated by the standard evaluation
procedure"* — a breakdown construction is fully prize-eligible.

Active work and all progress gates: **[`ns/README.md`](ns/README.md)**.

## Where it stands

The current object is an **obstruction map**: establish exactly which shapes a
finite-time breakdown is permitted to have, since the holes are the only places a
construction can live.

- **Lemma 0** (immediate): demanding invariance under $u\mapsto\lambda u(\lambda^2t,\lambda x)$
  forces $\alpha=\beta=1/2$. Scaling-invariant self-similarity *is* Leray's ansatz.
- **Lemma 1** (proved here, checked in [`ns/scaling.py`](ns/scaling.py)): for the
  first-kind family, **assuming** $\alpha+\beta=1$, $\alpha>0$ and the standard pressure
  scaling, $\beta\ne1/2$ forces $\Delta U\equiv0$, hence $U\equiv0$ by Liouville.
  $\alpha+\beta=1$ is a **hypothesis, not a consequence** — powers of $\tau$ are linearly
  independent, so terms may cancel in groups; the full exponent classification is not
  attempted.
- **Lemma 2** (Nečas–Růžička–Šverák 1996; Tsai 1998): $\beta=1/2$ is Leray's ansatz and
  is excluded. Tsai Thm 2 needs only distributional NS plus a *local* energy bound.
  **Exact backward self-similar blowup is closed.**
- **Lemma 3** (Chae 2007): asymptotically self-similar blowup is excluded under $L^p$
  convergence hypotheses — the paper explicitly leaves weaker pointwise convergence open.
- **DSS reduction** (derived and checked in [`ns/dss.py`](ns/dss.py)): $\lambda$-DSS is
  exactly $S$-periodicity in $\tau=-\log(-t)$ with $S=2\log\lambda$; the steady case is
  Leray. **Then closed for Clay data:** $L^3$ is scale-invariant, so exact DSS makes
  $\lVert u(\cdot,t)\rVert_{L^3}$ $\tau$-periodic; finite at one Schwartz slice gives
  $u\in L^\infty_tL^3_x$, and Escauriaza–Seregin–Šverák removes the singularity **for
  every $\lambda$** (Chae–Wolf Remark 1.2). The large-$\lambda$ hole exists only
  *outside* $L^3$, which Fefferman's condition (4) forbids.
- **Lemma 4** (proved here, checked in [`ns/forcing.py`](ns/forcing.py)): (C) lets the
  claimant choose a force $f$, but (5) with $\alpha=m=K=0$ gives $|f|\le C$. For an
  *exactly* self-similar $u$ a bounded $f$ is subcritical — but that alone does **not**
  give $f\equiv0$, since the pressure need not share the scaling ($u\equiv0$ with
  $f=\nabla\phi$, $p=\phi$ is a counterexample). Both branches instead go through the
  **Leray projection**: every term of $N(u,p)$ carries $\lambda^3$, so the two forms of the
  equation differ by a pure gradient, and $\mathbb{P}$ — degree-0 symbol, hence commuting
  with dilations — gives $(\mathbb{P}f)(x,t)=\lambda^3(\mathbb{P}f)(\lambda x,\lambda^2t)$.
  Bounded plus covariant forces $\mathbb{P}f\equiv0$, so $f$ is a pure gradient and the
  equation is **effectively unforced** — never $f\equiv0$. **This closes the *exact*
  profile routes with forcing.** It says
  nothing about *asymptotic* or localized profile routes with forcing, where $u$ is only
  approximately covariant — those stay **open**.
- **Surviving hole:** **localized / asymptotically DSS with unbounded $L^3$ norm, i.e.
  $\limsup_{t\uparrow T}\lVert u(t)\rVert_{L^3}=\infty$.** Chae–Wolf Thm 1.5's convergence lives on shrinking
  balls only, so it yields no global $L^3$ bound and ESS cannot be invoked. The sharp
  separator is $\sup_t\lVert u(t)\rVert_{L^3}<\infty$. Their $\lambda\approx1$ restriction
  is **essential to the method** — the proof degenerates DSS into exact self-similarity by
  compactness as $\lambda_j\to1$. Beyond this: Type II (no construction) and
  non-self-similar forced blowup (neither built nor excluded).

**Verification.** Fefferman, Tsai, Chae, Chae–Wolf and KNŠŠ were read as primaries.
**Nečas–Růžička–Šverák and Escauriaza–Seregin–Šverák could not be retrieved** (publisher
blocks) and are tagged `[UNVERIFIED]`; the DSS closure rests on ESS, so that tag
propagates. Discharging it is the next task.

Direct consequence: every profile-based program, including anything descended from
Chen–Hou or the DeepMind unstable-singularity work, is **structurally inapplicable** to
the Clay statement — not merely unbridged.

## Layout

| Path | Contents |
|---|---|
| `PROGRESS.md` | Dated status ledger. Read first to resume. |
| `ns/` | **Canonical target.** Obstruction map, lemmas, Clay-tied experiments. |
| `ns/README.md` | The statement, Lemmas 1–3, the surviving hole, next deliverable. |
| `ns/scaling.py` | Symbolic verification of Lemma 1. |
| `docs/CMI_RULES.md` | Verbatim prize rules; corrections issued. |
| `docs/GAP_MAP.md` | Model-equation ladder; technique vs conceptual gaps. |
| `docs/SURVEY.md` | All seven problems: proximity scores and path verdicts. |
| `docs/SOURCES.md` | Bibliography; what each reference actually proves. |
| `ccf/` | **Paused.** Technique validation only — see below. |
| `.venv/` | numpy, scipy, mpmath, sympy, python-flint. |

## Scope discipline

`ccf/` is a 1D inviscid model equation with **no viscosity bridge and no domain
bridge**. It is not progress on the Clay problem and is never reported as such. It is
retained because it produced reusable machinery — an exactly-invertible Hilbert-transform
basis closed under $H$, $d/dy$ and $\int_0^y$ — and one hard lesson about validating an
operator on the decay class you actually need rather than a convenient one. It is
**paused** pending a Clay-tied use.

No experiment enters `ns/` unless tied to a named lemma or hole in the obstruction map.

## Running

```sh
cd math
./.venv/bin/python ns/scaling.py                  # Lemma 1
./.venv/bin/python zeta5/src/zudilin_rec.py       # zeta(5) gate-1 baseline
./.venv/bin/python zeta5/src/recsearch.py --selftest
python3 h10q/h10q.py                              # H10/Q suite (~5 s; --extended ~22 s)
./.venv/bin/python ccf/tests/test_hilbert.py      # paused subproject
./.venv/bin/python ccf/tests/test_tail_basis.py
```

## Evidence discipline

1. Every mathematical claim carries a source: arXiv ID, DOI, or an in-repo derivation.
2. Numbers copied from papers name the figure or table they came from.
3. Not personally verified against a primary source → tag `[REPORTED]`.
4. Inferred rather than read → tag `[INFERENCE]`.
5. Conflicting values are recorded as conflicts, not silently resolved.
6. A model-equation result is never described as progress on the Clay problem.
7. **A no-go theorem is worth exactly its hypotheses.** Transcribe them; the surviving
   hypotheses are where a construction can live.
