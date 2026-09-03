# Millennium / open-problem attack repository

**Start here:** [`RESULTS.md`](RESULTS.md) is the authoritative results index **for the
union-closed programme** — `uc/` and `LIU_H1/` — which is what its own title states. It is
**not** repo-wide: per-target current state is owned by each `<target>/README.md`, and the
[Targets](#targets) table below is the repo-wide map. This README is that map; everything
authoritative lives here, and nothing in `scratch/` is canonical.

Sibling repositories `../physics/` and `../cs/` carry one `RESULTS.md` row per target. This
repository does not, because its results index grew as the union-closed project's own
document. The scope is stated rather than silently implied; do not cite `RESULTS.md` for
any target outside `uc/`/`LIU_H1/`.

## Tracking contract (SSOT)

One fact, one owner. Four layers, no duplication:

| layer | owner | contents |
|---|---|---|
| results index, union-closed programme only | [`RESULTS.md`](RESULTS.md) | headline result, theorem/lemma status, campaign ladder, reproduction, and adopted next targets — for `uc/` and `LIU_H1/` |
| chronological ledger | [`PROGRESS.md`](PROGRESS.md) | dated, newest-first session entries; every claim carries its verification; retractions inline |
| current state per target | `<target>/README.md` | what is PROVED / CONDITIONAL / NUMERICAL right now, file inventory, how to run |
| run artifacts | `<target>/campaigns/<UTC-timestamp>_<uuid>_<code-hash12>/` | immutable frozen snapshots + committed results; inventoried in `<target>/campaigns/README.md` |
| generated compact state, per target | `<target>/state.json` | owner-maintained semantic fields (`problem`, `current_gate`, `headline`, `next_action`); mechanics refreshed by `scripts/campaign.py state refresh` |

Generated, not narrative: when the one-line summary in `<target>/state.json` and an
owner document above disagree, the owner document wins — fix [`RESULTS.md`](RESULTS.md),
`<target>/README.md`, or [`PROGRESS.md`](PROGRESS.md) first, then re-run
`scripts/campaign.py state refresh` so the generated copy follows.

Naming rules:

* **One problem, one top-level folder**, short lowercase name (`uc/`, `r55/`,
  `zeta5/`, `ns/`, `ccf/`). No dated top-level folders, no second folder for
  the same problem.
* Campaign directories are producer-generated
  (`cert3_<UTC>Z_<uuid>_<hash>`); never hand-created, never edited after
  launch. Disposable rehearsals go under `<target>/campaigns-smoke/`.
* Update the ledger and the owning README **in the same session** as the work;
  a result that is not in `PROGRESS.md` does not exist.

> **Naming drift, 2026-08-15 — resolved:** `ramsey-r55/` appeared beside
> `r55/` (writer: a parallel agent session). Merged 2026-08-15: its data files
> were byte-identical duplicates of `r55/data/` (SHA-256 verified); its unique
> deliverables now live at `r55/src/structural_constraints.py`,
> `r55/src/n49_gluing_sat.py`, `r55/notes/structural_constraints_2026-08-15.md`,
> `r55/data/structural_tables.json`. `ramsey-r55/` was deleted with user
> confirmation after the merge (2026-08-15). Incident during the merge: this
> README was briefly overwritten and restored verbatim from the 15:03:06 APFS
> snapshot (evidence copy:
> `r55/notes/incident_2026-08-15_readme_clobber_recovered_copy.md`; details in
> `PROGRESS.md`).


## Targets

| directory | problem | state |
|---|---|---|
| [`e389/`](e389/README.md) | **Erdős Problem #389: consecutive-product divisibility** — *active* (opened 2026-08-25) | **PROVED:** exact valuation/zone formulas; exact band/local small-prime controls; compensation-good window separation; guaranteed prefix-repair cone; no fixed CRT class can suffice (Dirichlet); shift-spike/first-term identity; **half-window translation evicts every large blocker with $\lceil m/2\rceil$ the sharp minimal offset**; **the powersmooth necessity bound $p^{v_p(w)}<\sqrt{2w}$**; **the $\log 2$ short-cofactor density, giving witness upper density $\le1-\log 2$**; **the reduction of every witness shift to one adjoined integer**; **the small-prime tier is asymptotically free, so Erdős #389 reduces exactly to one run conjecture and no band or CRT class can obstruct**; **congruence-forced powersmooth mass costs a modulus above the searched range (Mertens), so a size-forcing class names the witness rather than predicting it**. **Every arithmetic progression of length $\approx2\sqrt{2\,\mathrm{top}}$ contains a non-witness — Bertrand and one inversion, no prime localisation — giving effective class eviction at $k<4M^2$ (beating Linnik's $M^{5.18}$), block eviction at $M\le\sqrt N/2$ (beating GRH, which reaches only $\sqrt N/\log^{2+}N$), and the price $r\ge\frac12\delta\sqrt N$ for any congruence programme covering density $\delta$; sharp to a factor 2 by an $m$-smooth obstruction at $M^2\ge N$**. **$R(1)$ — the compensation-good integers have positive lower density $\ge0.0293$ (Mertens + Vinogradov), proving the first rung of Conjecture R by bounding the Bonferroni overlap from below on the *prime* side**; **the reduction is finitary — one scale with $\#(G_m\cap[1,X])>4\pi(m)mX^{1-1/(24\log m)}$ yields a witness**. **TYPE-I/TYPE-II PAIR THEOREM:** Pascadi--Yang gives positive density with $P^+(n)<P^+(n+1)<x^{3/8+\varepsilon}$, so the short-cofactor pair is unconditionally closed; the $5/8\to3/8$ transfer is Theorem 19.1. **DOUBLE-MASK FRONTIER:** exact Fourier expansion shows $15/16$ is one-sided only and the bilinear lift is $5/4$. Theorem 23.1 carries the mask through Poisson: dual mass stays controlled, but the moving center expands off-diagonal support from $x^{1/2}$ to $x^{7/8}$. **CERTIFIED FINITE:** covering certificates evict every class to $\theta=0.5702$ at $m=27$, $N=10^{12}$, within $0.05\%$ of the measure ceiling; an exhaustive census of $[10^7,2\cdot10^7)$ evicts every class of every tested modulus to $\theta=0.990$ for $m=27,51$; all 27 OEIS witnesses; repaired schema-2 atlas and 101,646-case diff; any $m=27$ witness has $k>1{,}000{,}050{,}000$; no $m=27$ witness in the contiguous $10^9$ candidates from the published $m=26$ shift target; all 12 former survivor eviction systems now classified exactly at 13 digits. **OPEN:** Conjecture R($L$) for $L\ge2$. For $R(2)$ Input VII is now exactly the paired moving-center Kloosterman estimate (94) and its one-mask cross variant; the post-Cauchy support is $x^{7/8}$ against spectral level $x^{1/2}$. After it, only factors of $\ell$, ultra-small factors of $n$, and small repeated primes remain; selected-prime repetition and all repeated factors above $x^\zeta$ are already negligible. Above $\sqrt N$ the gap is now one effective constant: Balog--Pomerance and Shiu already bound survivors for every individual modulus with a $3.26$ margin, and no averaged theorem can be used (see `e389/LOCALIZATION.md`). |
| [`LIU_H1/`](LIU_H1/AUDIT.md) | **Liu's Hypothesis 1** (union-closed programme) — *active* | State is owned by [`RESULTS.md`](RESULTS.md), which **declares** `LIU_H1/` in scope alongside `uc/` (see README top) and carries 14 path-anchored references to it; the folder additionally holds `AUDIT.md`, `paper/`, `literature/manifest.sha256`, and `verification/logs/`. Listed here because the Targets table is the repo-wide map and omitting a declared in-scope programme left it invisible to that map. |
| [`uc/`](uc/README.md) | **Union-closed sets conjecture** (Frankl 1979) — *active* | **CANDIDATE THEOREM:** Campaign I's replayed 488,465,854-node Arb certificate plus a human-audited support reduction/entropy bridge gives explicit candidate frequency \(0.3820660112501052\ge\psi+10^{-4}\); not end-to-end machine verified. **LIU H1 CANDIDATE:** complete human-audited Gram proof with independent exact algebra checks. **LIU H2 PROVED AND THE CONSTANT UNCONDITIONAL (2026-09-02):** H2 holds for all Borel probability measures and every conditionally i.i.d. coupling at the exact root-defined constants (two certificates, both independently audited), so Liu's Proposition 3 gives the unconditional union-closed constant \(c'=1-m^*=0.38270908791873502993\ldots\) (certified enclosure \(\pm4.95\times10^{-46}\); Liu's printed \(0.382709087918741\) is wrong in its last two digits) — strongest bound in the ledger, pending external refereeing; see `RESULTS.md` Fourth result. **GATE B PROVED/CERTIFIED:** two Cartesian-power constructions prove \(c_{\rm cl}^\star=\infty\), with exact rational, Arb, and clamp-free dyadic certificates. The numerator obeys \((3/640)n-3/80\le\Lambda(n)\le4n/5\). **NEW LOCAL RECORD:** a normalized eight-coordinate, fifteen-row core at defect \(14/45<2/5\) has an exact negative fixed order; the sharp Plackett--Luce clone law yields a finite symbolic cap/Reimer family with the same defect and exact \(A_+\in[-0.000084288238,-0.000084288237]\). The clone repeats columns, so the separating frontier remains \(1144/1875\); \(c_{\rm loc}\) as defect tends to zero remains open. Start at [`RESULTS.md`](RESULTS.md). |
| [`r55/`](r55/README.md) | **Ramsey number $R(5,5)\le45$** — *active* | Gates 1–3 and eleven ≤46 strata remain exactly reproduced. The new fixed-five involution lane exhausts 48,865,656 exact `W` row domains: 705 `D5` Ramsey-support orbits reduce to 155 `W`-square orbits, 108 signed-square orbits, then zero Ramsey-good orbits after 73,336 switching-gauge completions, all with both `K5` and `I5`. **Certification is incomplete:** pinned VeriPB/CakePB accepts a sound deletion-free DFS certificate for source 0 only (1/550 W-negative supports); 549 W-negative and 155 deeper cases remain uncertified, so this is candidate-level evidence, not a theorem-grade branch closure or an $R(5,5)$ bound. |
| [`zeta5/`](zeta5/README.md) | **Irrationality of $\zeta(5)$** — *active* (opened 2026-08-15) | gate 1 passed: Zudilin's third-order recursion (arXiv:math/0206178) reproduced and certified exactly (convergents, integrality (6) to $n{=}300$, signs, roots, rates); the honest baseline **deficit is $5+\log\lvert\mu_2\rvert=3.914$** against the $D_n^5$ threshold; recurrence-search tool validated by blind re-derivation of the recursion from raw $q_n$ data. Gate 2: the triage's four-week falsifiable construction sweep |
| [`kobon/`](kobon/README.md) | **Kobon triangle problem, $n=10$** (Fujimura) — **CLOSED 2026-08-15** | **Theorem: $K_{\rm gen}(10)=25$** under the broad convention (parallels, multipoints, crossed triangles allowed — no located source proves this case; the classical charging lemmas fail over $\mathbb{Q}$, exact counterexamples in `kobon/report.md` §6). Lower bound: exact integer certificate, re-verified. Upper bound: 11-cube necessary-clause CNF cover, every cube with an independently `drat-trim`-verified UNSAT proof (9 shipped DRATs + 2 cubes re-proven via byte-verified tautological 8-way splits, 16/16 subproofs `s VERIFIED`) |
| [`qec/`](qec/README.md) | **Automated co-design of quantum LDPC codes, syndrome-extraction circuits, and decoders under circuit-level noise** — *active* | Non-CSS PBB programme, two proved sectors. **Z-side:** exact trade law $d_Q>d_Z(P)\Rightarrow k_Q=k_P-T(P)$ (Theorems G/H/I), 63/202 parents family-closed including Gross and every catalogue $[[144,12,12]]$ parent; all catalogue rows through $n\le180$ CSS-envelope-dominated; seven certified distance reversals, all envelope-dominated. **X-side:** collapse is exactly classified by $I=\operatorname{Ann}_R(A,B)$ (Theorem J-G, 202/202), odd lattices cannot collapse, no weight-$\le3$ pair is mixed (653,022,021-pair census), and weight 4 is sharp. **Odd fixed-point screen:** 65-lattice algebraic census through $n=360$; exact comparison closure through $n=234$ with 4,862 classes / 150,581 pairs and all 4,658 referenced classes dominated. EXP-067 independently certifies Liang et al.'s two published $[[234,8,18]]$ rows by replayed rooted connected-cluster enumeration. Target `12_6_0193` is genuinely non-CSS with exact $d=12$; circuit evidence remains negative and the fault-tolerant mixed-stabilizer construction remains open. |
| [`h10q/`](h10q/README.md) | **Hilbert's Tenth Problem over $\mathbb{Q}$** — *active* (opened 2026-08-12, canonicalized 2026-08-15) | **conditional record only:** classical Schinzel H implies a $6$-universal definition of $\mathbb Z$ (L19–L22); H10/$\mathbb Q$ remains open.  Unconditionally, L23 isolates the fixed-family dispersion barrier, L24 kills both unscaled diagonals at $2$, L26–L30 close the reciprocal local lift and lower the coupled cover to degree $4$, and L31 proves one exact globally good reciprocal member at \(w=13\).  L32 gives automatic-\(\Phi\) local bases; L33 proves a positive unweighted two-large pair main; L34 shows that the unrestricted Chebotarev root-pair envelope at \(\vartheta=0.49\) is \(0.974960\ldots<1\), so a conditioned prime-weighted root-pair upper bound by it would imply AP1 without needing exact odd valuations.  That estimate, a globally controlled quartic root, AP1, and Schinzel-H removal remain unproved; refereed $\forall_{10}$ and unrefereed $\forall_7$ records are untouched. |
| [`ising3d/`](ising3d/README.md) | **Exact 3D Ising model program** — *active* (relocated 2026-08-14) | **not solved**; 26 verified waves. Complete connected-graph local-term classification; all-`L` ladder fixed-basis no-go; every-coupling non-bipartite spectral theorem; finite open-`2x3` full-spectrum exceptional set empty. Wave 25 established the first nonplanar fixed-Walsh boundary matchgate defect, exact open-`2x4` traces through nine, W-law Smith shells, and the radius-one cubic contact no-go. Wave 26 proves the open-`2x4` trace-nine incidence lies on one nonzero `q` polynomial of degree at most `182400`, using two exact 96-dimensional modular quotient determinants plus a singular eight-mode control; the primitive norm and emptiness remain open. Certified `K_c` interval: `[0.2138368062108697304302094293879379891333, 0.2527310098586630030260020266135701299926]`. |
| [`ns/`](ns/README.md) | Clay Navier–Stokes, Fefferman (C) on $\mathbb{R}^3$ — *paused* | obstruction map; blocked on two publisher-inaccessible primaries (NRS, ESS) |
| [`ccf/`](ccf/README.md) | Córdoba–Córdoba–Fontelos model — *paused* | technique validation only; prize-irrelevant |

`uc/` was chosen over `ns/` on one criterion: **every primary is on arXiv.** The
Navier–Stokes attempt stalled because NRS and ESS could not be retrieved, and an
obstruction map that rests on unread primaries cannot be pushed further
responsibly.

### Non-target directories

Two top-level directories are **not** problem targets and carry no Targets row. Both were
previously absent from every table in this repository, which contradicted the claim above
that everything authoritative is mapped here; recorded 2026-08-30.

| directory | what it is | why it has no Targets row |
|---|---|---|
| `dist/` | distribution bundles for `r55/` — five dated `.zip` evidence archives plus [`dist/MANIFEST.md`](dist/MANIFEST.md), whose head records the `MIXED_PAIR_LIFT_CLOSED` disposition | packaged evidence, not a problem; `r55/README.md` remains the owner of that programme's state |
| `prime-gaps/` | a **vendored external checkout**, not authored here: `PrimeGapsLib` from `https://github.com/morluto/PrimeGapsLib.git` (Axiom Math), with its own `.git/` history, `lean-toolchain`, and a `lakefile.toml` requiring `PrimeNumberTheoremAnd` at a pinned rev | third-party code under an independent history; nothing in it is a claim of this repository, and its `STATUS.md` is that project's, not ours |

Do not cite either directory as a result of this repository, and do not add campaign
artifacts under them.

## Target-selection surveys

* [`docs/AI_FRONTIER_TARGETS_2026-08-25.md`](docs/AI_FRONTIER_TARGETS_2026-08-25.md)
  audits the July-August 2026 AI-mathematics results, ranks seven checker-backed
  frontier targets, selects the union-closed closure-defect programme, and
  records its first certified falsification.
* [`docs/SURVEY.md`](docs/SURVEY.md) scores the six unsolved Millennium
  problems and records why none currently offers a better checkable foothold.

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
python3 h10q/h10q.py                              # H10/Q suite (58.56 s; --extended 179.74 s)
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
