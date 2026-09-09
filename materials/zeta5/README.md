# ζ(5) irrationality campaign

**Canonical target.** Prove $\zeta(5)\notin\mathbb{Q}$.

Fallback outcomes that still publish (per the 2026-08-13 triage): a new
*certified* Apéry-style construction with materially better arithmetic (faster
decay or a proved denominator gain), a certified integrality theorem, or a
sharp obstruction ruling out a family of constructions. A recurrence alone, or
numerical convergence alone, publishes nothing.

## Status quo (source-verified 2026-08-15)

- **Open.** The corridor's landmarks:
  - Apéry 1979: $\zeta(3)\notin\mathbb{Q}$ — second-order recursion; decay
    $-4\log(1+\sqrt2)=-3.5255$ beats the $D_n^3$ threshold $-3$ ($D_n=\mathrm{lcm}(1..n)$).
  - Ball–Rivoal 2001: infinitely many $\zeta(2k{+}1)$ are irrational.
  - Zudilin 2001: at least one of $\zeta(5),\zeta(7),\zeta(9),\zeta(11)$ is
    irrational (DOI 10.1070/RM2001v056n04ABEH000427). No individual odd value
    past 3 is settled.
  - Zudilin 2002, arXiv:math/0206178 (primary, read first-hand): a third-order
    Apéry-like recursion for $\zeta(5)$, presented explicitly as a fast-computation
    algorithm — its inclusions "do not allow oneself to prove the irrationality".
    Reproduced and certified here (gate 1).
- Recent corridor movement `[REPORTED — abstracts read, papers not yet studied]`:
  the 2-adic analogue $\zeta_2(5)$ **is** irrational — Calegari–Dimitrov–Tang
  (arXiv:2408.15403 machinery) and an Apéry-style reproof with
  $\mu(\zeta_2(5))\le20.35$ by Lai–Sprang–Zudilin (arXiv:2505.05005, rev. 2026-05).
  Apéry-style machinery is alive in this corridor; the archimedean case stays open.
- arXiv:2407.07121 claimed $\zeta(5)\notin\mathbb{Q}$; **withdrawn by its author
  2025-05-19**. No accepted claim exists (checked 2026-08-15).
- Provenance: target chosen over $\zeta(7)$ / Catalan / $\gamma$ / the CDT constant
  by the 2026-08-13 triage (session `019ffb23-…`, report `AperyTarget.md` in the
  migration bundle): CDT already solved, $\zeta(7)$ strictly weaker starting
  point, $\gamma$ has no comparable recurrence.

## The race (PROVED, elementary)

If $\zeta(5)=A/B$ and $\ell_n=q_n\zeta(5)-p_n$ with $q_n\in\mathbb{Z}$,
$2D_n^5p_n\in\mathbb{Z}$, $\ell_n\neq0$, then $B\cdot2D_n^5\ell_n$ is a nonzero
integer, so $2D_n^5|\ell_n|\ge1/B$ for every $n$. Hence
$2D_n^5|\ell_n|\to0$ (with $\ell_n\ne0$) proves irrationality. Since
$\tfrac1n\log D_n\to1$ (PNT), an attack in this form needs

$$\lim_n\tfrac1n\log|\ell_n| < -5.$$

Zudilin's recursion achieves $\log|\mu_2|=-1.08607936$.
**Deficit: $5+\log|\mu_2|=3.9139$** (with the sharper, empirically verified
inclusions (6); $5.9139$ against the proved (14)-type inclusions). Equivalently:
need $|\mu_2|<e^{-5}=0.006738$, have $0.337537$. Apéry's $\zeta(3)$ margin, for
scale: needed $<-3$, achieved $-3.5255$.

So any Apéry-style attack must (a) find constructions with radically faster
decay, (b) prove denominator gains (Krattenthaler–Rivoal-type denominator
conjectures) that shrink the exponent 5, or (c) eliminate the $\zeta(3)$
contaminant from joint forms — the same $q_n$ here carries
$\tilde\ell_n=q_n\zeta(3)-\tilde p_n$, a simultaneous-approximation structure.
Gate 2 sweeps (a) systematically; (b)/(c) get attention only on a surviving
candidate.

## Go/no-go gates (failing a gate stops the campaign)

1. **DONE 2026-08-15 — baseline certified** ([`src/zudilin_rec.py`](src/zudilin_rec.py),
   record [`data/BASELINE.json`](data/BASELINE.json)).
   Recursion (1) of math/0206178 transcribed (PDF exponents reconstructed; the
   transcription is pinned by degree-9 homogeneity, the printed characteristic
   polynomial, the printed $q_2$, and the paper's convergent table) and verified:
   printed convergents $n\le7$ exact; printed error bounds at $n=3..7,10,20,50$;
   integrality (6) ($q_n\in\mathbb{Z}$, $2D_n^5p_n\in\mathbb{Z}$,
   $2D_n^3\tilde p_n\in\mathbb{Z}$) machine-verified to $n=300$; sign patterns (3);
   root decimals (5) and the $\mu_i=\lambda_j\lambda_k$ pairing; rates (4)–(5) at
   $n=300$ [NUMERICAL]. Scoping fact: the $n{=}1$ instance of (1) holds for
   $\{q_n\}$ but **not** $\{p_n\}$ — the recursion is imposed for $n\ge2$ only.
   Search tool validated blind ([`src/recsearch.py`](src/recsearch.py)): from raw
   $q_0..q_{63}$ it re-derives (1) as the *unique* order-3/degree-9 recurrence
   (exact primitive-vector match) and finds nothing at order 2 (deg ≤ 12) or
   order 3, deg ≤ 8.
2. **CLOSED 2026-09-08 — FROZEN-NEGATIVE, zero survivors** (run
   [`20260908T150622Z_6a444c77_79aa3a46d868`](campaigns/20260908T150622Z_6a444c77_79aa3a46d868/),
   prereg [`pre_statement_gate-2-construction-sweep.md`](pre_statement_gate-2-construction-sweep.md)).
   The frozen 44-probe ladder (F1 very-well-poised deformations of series (7),
   slopes {−1,0,1}² × w ∈ {4,5,6}; F2 Vasilyev-type t = 1..3, s = 0..4; F3
   Ball–Rivoal-type r = 1,2; 80 exact terms each; orders 2–3, degree ≤ 14)
   retains nothing: 33/44 no recurrence in window, 7/44 slower than baseline,
   and all 4 fast-looking promotions withdrawn on audit (3 degenerate
   content-free, 1 rate artifact). Baseline stands:
   $\mu_2 = 0.33753726443403620704$, deficit $3.9139206383732506583$. The
   campaign records and stops here per the rule below; resumption requires a
   NEW family ladder (fresh triage decision). The gate as frozen:
   enumerate low-order polynomial recurrences from exact data of candidate
   constructions (deformations of the very-well-poised series (7) of
   math/0206178, Vasilyev-type integrals, Ball–Rivoal-type series; guess with
   `recsearch`, certify by creative telescoping); a candidate is retained only
   if all three hold (exact symbolic recurrence identity; integrality /
   denominator control; certified exponential error bounds). The honored rule:
   **zero survivors better than the baseline deficit ⇒ the campaign fails
   here; record and stop** — which is exactly what the run executed.
3. **Arithmetic sharpening — only on a gate-2 survivor.** Denominator-gain
   proofs and contaminant elimination toward $\lim\tfrac1n\log|\ell_n|<-5$.

## Trust discipline (inherits [`math/README.md`](../README.md); additions)

- A numerically guessed recurrence is exact on its verified range only;
  validity for all $n$ is **CONDITIONAL** until certified symbolically
  (creative telescoping / algebraic proof).
- PSLQ/LLL output is a **candidate**, never evidence.
- Rate claims from finite $n$ are **NUMERICAL**; only Poincaré–Perron applied
  to a certified recurrence upgrades them (and Poincaré alone does not give
  the constant or exclude degenerate solutions — track which solution you hold).
- Integrality observed to $n=N$ (like (6) here) is labelled with its $N$;
  the paper itself only proves the weaker (14)-type inclusions.

## Layout

| path | contents |
|---|---|
| [`src/zudilin_rec.py`](src/zudilin_rec.py) | gate-1 baseline: exact recursion (1), table/integrality/sign/root/rate checks, deficit; writes `data/BASELINE.json`; run = selftest |
| [`src/recsearch.py`](src/recsearch.py) | recurrence guesser (FLINT exact nullspace, full-range verification) + PSLQ candidate detector (`--selftest` inside) |
| [`data/`](data/) | machine-written verification records |

## Running

```sh
cd math
./.venv/bin/python zeta5/src/zudilin_rec.py         # gate-1 baseline, N=300, <1 s
./.venv/bin/python zeta5/src/recsearch.py --selftest
```
