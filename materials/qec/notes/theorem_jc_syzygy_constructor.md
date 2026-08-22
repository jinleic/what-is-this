# Theorem J-C — the constructive syzygy demotion

**Main, 2026-08-19.** Companion to `notes/theorem_j_xsector.md` (J.0–J.6) and
`notes/j5_attempt.md` (J.5 refutation and its census). The refutation of J.5 is hereby
upgraded from *witnesses + sampled census* to an **explicit deterministic constructor
with a machine-complete catalogue verdict** (EXP-046).

All numbers are read from `results/processed/exp046_generic_syzygy.json`
(schema `exp046-generic-syzygy-v1`), produced by `experiments/exp046_generic_syzygy.py`
(GF(2)/numpy only, deterministic, no RNG, no sampling; 202 parents, ~228 s wall).
Cross-consistency with the earlier sampled census
`results/partial_runs/xsector/exp044_n180_hunt.json`: **zero mismatches** (the 49
census-positive parents are exactly inside the EXP-046 collapse set; no inversions).

---

## 1. The constructor (proved part)

Fix a BB parent $P=(A,B)$ over $R=\mathbb F_2[x,y]/(x^\ell-1,y^m-1)$, and write
$r:=\mathrm{wt}(A)+\mathrm{wt}(B)$ for the check-row weight (for every row of
$H_X=[\,A\;B\,]$, by circulant translation).

**Theorem J-C(i) — a kernel family of maximally degenerate perturbations.**
For every $c\in\ker(A)\subseteq\mathbb F_2^{\ell m}$ (the annihilator of the polynomial
$a$ in the convolution convention of the check matrix), the perturbation
$$C=\operatorname{circ}(c),\qquad D=0$$
satisfies $M=AC^{\mathsf T}+BD^{\mathsf T}=0$: it is **valid** and *every* row is a
single-row demotion candidate (L3 of `j5_attempt.md`). Moreover
$$\dim_{GF(2)}\ker A\ \ge\ \ell m-\operatorname{rank}[\,A\;B\,]\ =\ k_P/2 ,$$
so the family is nonempty on every parent with $k_P>0$.

*Proof.* $M=AC^{\mathsf T}$ is the circulant whose row $0$ is the vector $A(c)^{\mathsf
T}$; row-specifically $M[p,q]=(a\cdot c^{\mathsf T})[p-q]$ in the group algebra, and
$A(c)=0$ annihilates every row. Validity requires $M$ symmetric; $M=0$ is. Any $x$ in
rowspace $H_X$ meets the demotion-candidate space $W=\operatorname{rowspace}H_X\cap\ker[C\;D]$
for free since $M=0\Rightarrow W=\operatorname{rowspace}H_X$ entirely (L4). A is singular
whenever $k_P>0$ because $\operatorname{rank}A\le\operatorname{rank}[A\;B]=\ell m-k_P/2$.
$\square$

**Theorem J-C(ii) — the demotion predicate is a GF(2) rank test.**
Row $0$ of $H_X$ (weight $r$) is a *demoted parent X-stabilizer* — a genuine X-logical
of $Q$ of weight $r$ — **iff** its partner
$$z_0=(c,0)\ \notin\ S_Z+\Delta_C,\qquad \Delta_C=\{\lambda[\,C\;0\,]:\lambda[\,A\;B\,]=0\},$$
and $z_0\in\ker H_X$ holds automatically ($Az_0=0$). Consequently
$d_X(Q)\le r$, with $k_Q=k_P-\dim\bar\Delta$ and $\dim\bar\Delta=0$ in the
k-preserving case. The same construction with a general syzygy pair $(c,d)\in\ker[A\;B]$,
$C=\operatorname{circ}(c),D=\operatorname{circ}(d)$, has $M=0$ and
$z_0=(c,d)\in\ker H_X$ both automatic.

*Proof.* J.2/J.4(ii) demotion pairing with $\lambda=e_0$; $M=0$ conditions are (i).
$\square$

**Corollary (open-window collapse).** If $r<d_X(P)$ (the *open window*, L2) and the
rank test in (ii) passes for some $c\in\ker A$, then
$$d_X(Q)\ \le\ r\ <\ d_X(P),$$
i.e. a strictly distance-reducing PBB sibling exists at the price of the k-drop
$\dim\bar\Delta$ (usually zero).

---

## 2. The machine-complete catalogue verdict (EXP-046)

Constructor applied to **all 202 distinct BB parents** of the catalogue, deterministic
three-stage fallback: D=0 basis vectors of $\ker A$ → syzygy-basis vectors of
$\ker[A\;B]$ → pairwise syzygy-basis sums. Every accepted witness passes the dual
verification battery (path A: $z_0\notin S_Z+\Delta$ rank test; path B:
$x_0\notin S_X(Q)$ nullspace probe) and consistency invariants
($\mathrm{wt}(z_0)\ge$ bound, per J.4(ii)).

| scope | parents | demotion verified | collapse ($d_X(Q)\le6$ strictly below bound) |
|---|---|---|---|
| all catalogue parents | 202 | **192** | 99/100 bounded open-window |
| bounded parents (cert or pool) | 140 | — | **99**, of which **98 k-preserving** ($\bar\Delta=0$) |
| bounded, window closed ($r\ge d$) | 34 | 34 (at equality) | impossible |

- **All 192 demotions come from single D=0 basis vectors** — no search depth at all.
- Demotion cannot be claimed on **10 parents** (`phase2_75/76/77/83/84/87`,
  `9_6_0175`, `phase2_109`, `15_6_0256`, `30_6_0289`). 7 of them have certified
  $d=6=r$: window closed, demotion at equality is not a collapse channel anyway
  (and single-row demotion itself is blocked there by full basis+pairs search).
  3 are unbounded (window unknown). **1 is the genuine boundary instance:**
  `9_6_0175` ($\ell{=}9,m{=}6,k_P{=}8$, certified $d_Z(P)=8>6=r$):
  - its row-overlap bound is $\le 2$, so **any 2-row sum has weight $\ge 8=d$** —
    *multi-row* demotion below $d$ is arithmetically impossible; the single-row
    channel is the only possible collapse mechanism;
  - the entire syzygy kernel (dim $58$) was probed through basis vectors + all
    $\binom{58}{2}=1653$ pairs plus the full ker($A$) D=0 family: **no demotion**.
  - So `9_6_0175` is single-row-demote-immune to machine-tested depth; an
    exhaustive $2^{58}$ certificate remains open (Row 2's gap is now *exactly* the
    quotient question $\mu(\mathrm{Syz})\cap\mathrm{im}\Delta=?$).
- Collapse bounds among the 99: $d\in\{8:34,\ 10:26,\ 12:31,\ 14:7,\ 16:1\}$,
  $n\in\{72..360\}$.

**Theorem J-C(iii) (catalogue-complete, machine-verified).** *Every* bounded
catalogue parent with an open window admits an explicit, deterministically
constructible, k-preserving syzygy collapse — **with the single exception of
`9_6_0175`,** whose only possible channel (single-row) is immune to the tested
depth of the syzygy kernel.

---

## 3. Consequences

1. **The X-side hazard is now a computable invariant, not a search artifact.** The
   decision "does this parent admit a syzygy collapse" costs one GF(2) nullspace +
   rank test per candidate — it belongs in every qLDPC candidate audit pipeline
   (cf. Theorem F two-sided certification).
2. The collapse is structurally *generic*: 98/99 collapses at full $k$, mechanism a
   pure-$R$ annihilator phenomenon. The earlier census' "49/49 bounded parents"
   verdict extends to **99/100 bounded open-window parents**.
3. **The exception isolates the boundary of the theory.** `9_6_0175` has minimal
   nullity ($\dim\ker A=4=k_P/2$, half the flagship's $12=k_P$); whether the
   kernel-small parents are exactly the immune ones is the next structural
   question (all 10 failures have $\dim\ker A\in\{4,12\}$, and the nullity
   histogram over 202 parents sits at $\{2,4,6,8,10,12,20,24,60\}$ — the small-kernel
   hypothesis is falsifiable on the nullity=2 parents, all of which demote).
4. Machine-readable: every one of the 192 witnesses ship full term lists
   (C[, D]) in the artifact; independent reproduction is one GF(2) rank test away.

---
*Provenance: experiment `experiments/exp046_generic_syzygy.py`; artifact
`results/processed/exp046_generic_syzygy.json`; row-overlap certificates computed inline
2026-08-19 (max pairwise support intersection of $H_X$ rows per parent).*
