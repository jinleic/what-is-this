# Theorem J-E — the exact syzygy-collapse characterization

**Main, 2026-08-19.** Upgrade of Theorem J-C (`notes/theorem_jc_syzygy_constructor.md`)
and completion of the X sector of the PBB perturbation program
(`notes/theorem_j_xsector.md`, `notes/j5_attempt.md`). The demotion decision is not only
polynomial-time — it has a **closed algebraic characterization**, proved for one
direction, proved in general via Nakayama, and machine-confirmed exhaustively on all
202 catalogue parents (EXP-047).

Artifacts: `results/processed/exp046_generic_syzygy.json`,
`results/processed/exp047_exact_demotion_decision.json`;
scripts `experiments/exp046_generic_syzygy.py`, `experiments/exp047_exact_demotion_decision.py`;
invariant check inline (32 s, 202/202 zero mismatches).

---

## 1. Statement

Let $R=\mathbb F_2[x,y]/(x^\ell-1,y^m-1)$ and $P=(A,B)$ a BB parent. Put

$$M\ :=\ \ker H_X\,/\,S_Z\qquad\text{(GF(2)-dimension }k_P\text{, an }R\text{-module)},$$

$$L_{\rm pre}\ :=\ \{\lambda\in R:\ \lambda[A\;B]=0\}\ \subset\ R\qquad\text{(an }R\text{-submodule)}.$$

**Theorem J-E.** The following are equivalent:
1. No single-row syzygy demotion exists for $P$: for every pair $(C,D)$ with
   $M=AC^{\mathsf T}+BD^{\mathsf T}=0$, no check row becomes a logical of $Q$.
2. Every class $y\in M$ satisfies $y\in L_{\rm pre}\!\cdot y$.
3. $L_{\rm pre}\!\cdot M=M$.
4. $1\in L_{\rm pre}+\operatorname{Ann}_R(M)$.

*Proof.* Single-row demotion for the perturbation $z=(c,d)\in\ker H_X$ (i.e. $M=0$,
J-C(i)) is J.4(ii)'s condition $z\notin S_Z+\Delta_z$ with
$\Delta_z=L_{\rm pre}\cdot z$. By $R$-equivariance of
$\mu:\ker H_X\to M$ (all four spaces are shift-invariant), $\mu(\Delta_z)=L_{\rm
pre}\!\cdot y$ for $y=\mu(z)$, and the pose is constant on $S_Z$-cosets
($\Delta(S_Z)\subseteq S_Z$), demote-classes are exactly
$\{y\in M:y\notin L_{\rm pre}\!\cdot y\}$: (2$\iff$1). (2$\Rightarrow$3) trivial;
(3$\Rightarrow$4) is Nakayama's lemma for the finite module $M$ over the commutative
ring $R$ (there is $x\in L_{\rm pre}$ with $(1-x)\in\operatorname{Ann}_R(M)$);
(4$\Rightarrow$2): $y=1\cdot y=(\lambda+a)\cdot y=\lambda\cdot y\in L_{\rm pre}\!\cdot y$.
$\square$

**Corollary (decidability and cost).** The single-row syzygy demotion problem for a BB
parent is decided by *one* GF(2) rank computation
($\operatorname{rank}(\mu(L_{\rm pre}\cdot Q))\ {\stackrel{?}{=}}\ k_P$, $Q$ a quotient
basis) — or, refinement-wise, the demote-classes are exactly
$\{y:y\notin L_{\rm pre}\!\cdot y\}$, computable class-by-class with one rank test each.
Every demotion realizes an X-logical of $Q$ of weight $\le\mathrm{wt}(A)+\mathrm{wt}(B)$.

**Corollary J-E′ ($\lambda$-partners of a fixed seed direction).** J-E(1) is stated
for demoted check rows ($\lambda=e_i$); the same invariant also covers the
$\lambda$-combined partners of a *fixed* seed direction. Fix an $M=0$ perturbation
$[C D]$ with seed partner class $y$, and suppose immunity holds, i.e. $1=l+a$ with
$l\in L_{\rm pre}$, $a\in\operatorname{Ann}_R(M)$. Then for every $\lambda\in R$,
$$\lambda y\;=\;(\lambda l)\,y\;\in\;L_{\rm pre}\!\cdot y\;=\;\mu(\Delta_z)\,,$$
so the J.4 partner $z=\lambda[C D]$ satisfies $z\in S_Z+\Delta_z$: no
$\lambda$-combined partner escapes the absorption of the seed class. Equivalently,
for the fixed seed, both the $\lambda=e_i$ cases *and* all coefficient combinations
are decided by the seed's own demote-class membership. Two cautions: the demote-set
$\{y:y\notin L_{\rm pre}\!\cdot y\}$ is **not** $R$-invariant as a set ($\lambda=0$
maps every class to the never-demoting class $0$), and EXP-047 classes parameterize
seed partners $z_y$; this corollary makes no claim about demoted vectors beyond the
fixed-seed partner family. **Boundary:** $M\neq0$ multi-row $\lambda[A B]$
light-partner demotions are *not* covered and are certified separately:
L3 equality-overlap + EXP-046 row-overlap bounds + exact light enumerations (\S3).

---

## 2. Machine verdicts (EXP-047, exhaustive)

Exhaustive class enumeration ($2^{k_P}$ quotient classes; 256 to 16.7M per parent; the
only >4M parents early-stop on the first realized class — exp046 tells us they have one):

| scope | parents | verdict |
|---|---|---|
| all catalogue parents | 202 | **192 demotion-realized exact, 10 demotion-immune exact** |
| flagship `12_6_0193` | — | **all 4095/4095 nonzero classes demote** (100 %) |
| consistency vs EXP-046 | 202 | zero upgrades, zero regressions |
| invariant `1∈Lpre+Ann_R(M)` ⟺ immunity | 202 | **zero mismatches** |

The proven-immune parents:
`phase2_75, phase2_76, phase2_77, phase2_83, phase2_84, phase2_87, phase2_109`
($\ell\in\{9,15\}$), `9_6_0175`, `15_6_0256`, `30_6_0289`.

---

## 3. The complete X-side classification of the catalogue

Combine J-E demotion-immunity with L3 (single-row demotion exists only when $M=0$) and
the row-overlap bounds of EXP-046 (max pairwise $H_X$-row support intersection per
parent; certified in the same session):

- **Class X-C (192 parents): a syzygy demotion exists, explicitly constructed.**
  With an open window ($\mathrm{wt}(A)+\mathrm{wt}(B)<d_X(P)$) this is a strict
  collapse $d_X(Q)\le\mathrm{wt}(A)+\mathrm{wt}(B)$ — on bounded parents: **99 of 100**,
  all but `9_6_0175`; 98/99 of the collapses are k-preserving. Strict-collapse
  certificates exist for **99** of the 192 (see the severity paragraph below);
  demotion-realization itself is exact for all 192.
> **RESOLVED (2026-08-20, EXP-051).** The prior basis+pairs scan is superseded by
> a one-shot CP-SAT channel query per parent: membership $x\in\operatorname{rowspace}(H_X)$
> encoded against the dual kernel basis, weight $\le d_P-1$, every original
> translated row forbidden by an exclusion clause. `INFEASIBLE` is the
> completeness certificate; any feasible witness is a multi-row light channel
> (each witness independently re-verified: physical weight, rowspace rank
> identity, exclusion from all original rows). **All eight X-M parents below are
> now CERTIFIED** ($9_6\_0175$ at cap 7, the seven phase2-style parents at cap 5),
> and the two X-U parents are decided in the other direction: both exhibit a
> verified weight-8 multi-row channel, so their X-monotonicity is not decidable
> by the light-exclusion route (syzygy immunity stays exact).

- **Class X-M (7 parents), CERTIFIED X-monotone (EXP-051)**: syzygy-demote-immune, and either
  $d=6=\mathrm{wt}(A)+\mathrm{wt}(B)$ with demotic weights from multi-row
  $\ge 6$ (`phase2_75/76/77/83/84/87`: overlap $\le3$, so $\lambda$-sums weigh
  $\ge6=d$: no strict drop ever), or
  **`9_6_0175` ($d=8$, overlap $\le2$): multi-row light stabilizers below 8 do not
  exist (certified by the light-enumeration: 54 lights = single rows), so
  $d_X(Q)\ge8=d_X(P)$ for EVERY valid perturbation — the unique parent with both a
  certified open window and provable X-monotonicity.**
- ~~**Class X-U (3 parents: `phase2_109`, `15_6_0256`, `30_6_0289`)**: no syzygy
  demotion at all (exact); no multi-row demotion below weight 8 (overlap $\le2$);
  verdicts for $d_X(Q)$ vs $d_X(P)$ await the parent's certified distance ($d_P>8$
  would leave a possible channel at weight exactly 8).~~ **RESOLVED 2026-08-20
  (EXP-049 dedicated certs):** `phase2_109` has $d_P=6=\mathrm{wt}(A)+\mathrm{wt}(B)$
  → window-closed, provably X-monotone (joins Class X-M; total X-M = 8);
  `15_6_0256` has $d_P=10$, $T=8=k_P$ and `30_6_0289` has $d_P=10$, $T=16=k_P$:
  both are Theorem-H family-closed (no PBB over them exceeds $d_Z(P)$ with
  $k>0$) — that rules out $d_Q>d_Z(P)$ but does **not** certify X-monotonicity;
  both remain in **Class X-U**. **2026-08-20 (EXP-051): the open sub-question is
  decided** — BOTH parents exhibit an independently verified weight-8 multi-row
  light $X$-stabilizer (rank-identity + row-exclusion checks persisted), so the
  light-exclusion route to X-monotonicity does not apply; their syzygy immunity
  stays exact and they remain X-U *by nature* (monotonicity of $d_X$ genuinely
  undecided). **Final X-side partition: 192 demotion-realized / 8 CERTIFIED
  X-monotone / 2 X-U with verified weight-8 channels.**
- **2026-08-20 — demote-set totality (EXP-050):** exact census of all nonzero
  quotient classes on the 188 parents with $k_P\le 12$: demote-fraction is EXACTLY
  $1$ on all 179 demoting parents and exactly $0$ on the 9 immune ones (zero mixed;
  zero mismatch vs EXP-047). The "measure-1" reading is superseded: see
  `notes/theorem_je2_demote_trichotomy.md` (Theorem J-E′′, proposed).

**Constructed collapse severity (upper bounds, EXP-046 witnesses).** The constructor
fixes its candidate at the first perturbation witness, so `witness.x0_weight`
($=\mathrm{wt}(\text{row }0)$, here 6 for all 192 demote parents) is a **constructed
upper bound** $d_X(Q)\le 6$, not a minimized demoted weight and not $d_X(Q)$ itself.
Under certified parent distances ($d_{\rm bound}$, 133 of 192): **99 parents have a
certified strict collapse** $d_X(Q)\le 6<d_X(P)$ (gap distribution
$2{:}34,\ 4{:}26,\ 6{:}31,\ 8{:}7,\ 10{:}1$; the largest constructed drop halves a
$d_P=16$ parent to a $d_X(Q)\le6$ sibling at full $k$); **29 sit at bound-parity**
(6 = d_bound: channel exists, strictness undecided); **5 have vacuous UBs** (d_bound
$\le6$). The remaining 59 await $d_P$ certificates (the exp039 sweep). A true
minimum-depth census is a separate minimization problem — deliberately not claimed.

The boundary instance `9_6_0175` is thus settled *in the safe direction*: it is the
**first catalogue parent whose entire PBB family is provably X-monotone** — the exact
dual of the flagship's collapse ubiquity (4095/4095).

---

## 4. Consequences

1. The catalogue's PBB X-side theory is **closed** modulo three parents' missing
   distance bounds. Together with the Z side (L1 universal monotonicity; Theorems
   G/H/I caps; EXP-045 residual), the perturbation program has named, decidable
   verdicts everywhere except one compute-fed gap (T/d certifications at $n\ge180$).
2. Theorem J-E(4) is a *one-rank-computation* audit for any proposed qLDPC parent:
   collapse-hazard is decidable in polynomial time, no search at all.
3. The immune family splits demote/immune at an algebraic boundary
   ($\dim\ker A$, $L_{\rm pre}$ structure, annihilator) — the "small-kernel
   hypothesis" is refuted (34 nullity-2 parents demote); immunity is an ideal-action
   property, not a dimension property.
4. For search practitioners: when demote-classes are 100 % (flagship), sampling the
   syzygy family gives a collapsing sibling with probability 1·(nonzero class); the
   $\sim2^{-34}$ relative density measured earlier was about the *set* Syz ⊂ V, not
   about inner-class abundance — inside Syz, collapse is generic.
