# Theorem J — the X-sector of the PBB no-go

**XSectorAnalyst, 2026-08-18.** Companion to `reports/paper_pbb_nogo.md` (Theorems G/H/I).

Theorems G/H/I of the paper are statements about the **Z-sector**. HostileReview
observed that the X-side is not untouched: the pure-X centralizer shrinks under the
perturbation (`phase2_58`: 40 → 12 — verified exactly below). The published proof of
Theorem G(ii) justifies the dimension identity with the sentence *"the X-sector of the
quotient is unchanged by the same argument applied to the (unmodified) second block"*.
That sentence is **wrong as an argument** — the second block is unmodified, but the
first block's Z-part $[C\;D]$ also binds pure-X vectors — although the identity it
supports is correct. This note derives the exact X-side analogue (**Theorem J**),
machine-checks it against Theorem G on four catalogue rows, and settles what can be
settled about the pure-X distance. Everything numeric below was computed this session;
provenance paths are given inline.

---

## 1. Conventions (as in the paper, §1)

$R=\mathbb{F}_2[x,y]/(x^\ell-1,y^m-1)$, $n=N=2\ell m$ qubits, and with $A,B,C,D\in R$
acting as $\ell m\times\ell m$ matrices,

$$H_Q=\begin{pmatrix} A & B &\big|& C & D\\ 0&0&\big|& B^{\mathsf T}& A^{\mathsf T}\end{pmatrix},\qquad
H_X=[\,A\;B\,],\quad H_Z=[\,B^{\mathsf T}\;A^{\mathsf T}\,].$$

Rows of $H_Q$: $r_\lambda=(\lambda A,\lambda B\,|\,\lambda C,\lambda D)$ and
$s_\mu=(0,0\,|\,\mu B^{\mathsf T},\mu A^{\mathsf T})$, $\lambda,\mu\in\mathbb{F}_2^{\ell m}$.
Validity of the perturbation is symmetry of $M=AC^{\mathsf T}+BD^{\mathsf T}$
(Lemma 0 / FR-006 correction; `proofs/pbb_structure.md`).

For a code with check rowspace $S\subseteq\mathbb{F}_2^{2N}$ and centralizer
$C=S^\perp$, write $V_X=\mathbb{F}_2^{N}\times\{0\}$ and define

- **pure-X centralizer** $\mathrm{Xcen}=\{x:\ (x|0)\in C\}$,
- **pure-X stabilizer space** $S_X=\{x:\ (x|0)\in S\}$,
- **pure-X logicals** $\bar X=\mathrm{Xcen}/S_X$ (a subspace of the logical space $L=C/S$),

and dually $\mathrm{Zcen},S_Z,\bar Z$. The paper's $S_Z(P)=\operatorname{rowspace}H_Z$ and
$\Delta=\{\lambda[C\;D]:\lambda[A;B]=0\}$, $\bar\Delta=(\Delta+S_Z)/S_Z$.

---

## 2. Theorem J (X-sector analogue)

**(J.0) Reversal symmetry of the parent.** $d_X(P)=d_Z(P)$ for every BB parent.

*Proof.* Let $\tau_0$ be coefficient reversal ($\tau_0 f \tau_0 = f^{\mathsf T}$ for the
shift matrices) and $\sigma$ the block swap. $\rho=\sigma\circ\tau_0$ maps
$\mathrm{Xcen}(P)=\ker_{\rm col} H_Z$ bijectively onto $\ker_{\rm col}[\,A\;B\,]=\mathrm{Zcen}(P)$
and $\operatorname{rowspace}H_X$ onto $\operatorname{rowspace}H_Z$, preserving weight
($B^{\mathsf T}x_1+A^{\mathsf T}x_2=0 \Leftrightarrow [\,B\;A\,](\rho x)^{\mathsf T}=0$).
Hence it maps pure-X logicals of $P$ to pure-Z logicals of $P$ weight-preservingly. $\square$

**(J.1) Pure-X centralizer.** For every valid perturbation,
$$\mathrm{Xcen}(Q)\;=\;\{\,x:\ H_Zx^{\mathsf T}=0\ \wedge\ [C\;D]x^{\mathsf T}=0\,\}
\;=\;\mathrm{Xcen}(P)\cap\ker[C\;D]\;=\;\ker_{\rm col}\begin{pmatrix}H_Z\\ [C\;D]\end{pmatrix}.$$

*Proof.* $(x|0)$ commutes with every $s_\mu$ iff $x\cdot(\mu B^{\mathsf T},\mu A^{\mathsf T})^{\mathsf T}
=(x_1B+x_2A)\mu^{\mathsf T}=0$ for all $\mu$, i.e. $H_Zx^{\mathsf T}=0$ — the parent condition.
It commutes with every $r_\lambda$ iff $x\cdot(\lambda C,\lambda D)^{\mathsf T}
=(Cx_1^{\mathsf T}+Dx_2^{\mathsf T})^{\mathsf T}\!\lambda^{\mathsf T}=0$ for all $\lambda$, i.e.
$[C\;D]x^{\mathsf T}=0$ — **new** constraint from the perturbation. $\square$

Consequences, with $R_{CD}=\operatorname{rowspace}[C\;D]$:

$$\dim\mathrm{Xcen}(Q)=\dim\mathrm{Xcen}(P)-\rho_X,\qquad
\rho_X:=\operatorname{rank}\begin{pmatrix}H_Z\\ [C\;D]\end{pmatrix}-\operatorname{rank}H_Z
=\dim\frac{R_{CD}+S_Z}{S_Z}.$$

Since $\Delta=L\cdot[C\;D]\subseteq R_{CD}$ ($L$ the left kernel of $[\,A\;B\,]$):

$$\boxed{\ \rho_X\ \ge\ \dim\bar\Delta\ }\qquad\text{[proved]}$$

— the X-centralizer loses **at least** as much dimension as the Z-logical sector. The loss
can be huge with zero Z-side cost: on `12_6_0193`, $\rho_X=46$ while $\bar\Delta=0$.

**(J.2) Pure-X stabilizer.** $S_X(Q)=\{\,\lambda[\,A\;B\,]:\ \lambda[C\;D]\in S_Z+\Delta\,\}$.

*Proof.* An element of $\operatorname{rowspace}H_Q$ is $(\lambda[\,A\;B\,]\,|\,\lambda[C\;D]+\mu H_Z)$;
it is pure X iff $\lambda[C\;D]\in\operatorname{rowspace}H_Z=S_Z$. Two preimages of the same
$x=\lambda[\,A\;B\,]$ differ by $\nu\in L$ and their $[C\;D]$-parts differ by $\nu[C\;D]\in\Delta$,
so membership is well defined modulo $\Delta$. ($S_X(Q)\subseteq\mathrm{Xcen}(Q)$ holds
automatically — group elements commute — and is asserted in code on every row tested.) $\square$

Note $S_X(Q)\subseteq \operatorname{rowspace}H_X\cap\ker[C\;D]$ always; the gap is the
**demotion space** of J.4.

**(J.3) Sector dimension identities.** *For every stabilizer code on $N$ qubits
(CSS or not): $\ \dim\bar X=\dim\bar Z=k$.* In particular, for $Q$:

$$\dim\mathrm{Xcen}(Q)=k_Q+\dim S_X(Q),\qquad \dim\mathrm{Zcen}(Q)=k_Q+\dim S_Z(Q),$$
$$\dim\mathrm{Xcen}(P)-\dim\mathrm{Xcen}(Q)=\rho_X,\qquad
\dim S_X(P)-\dim S_X(Q)=\rho_X-\dim\bar\Delta .$$

*Proof.* $V_X$ is Lagrangian ($V_X^{\perp}=V_X$). For any $S\subseteq C=S^\perp$:
$(C+V_X)^{\perp}=C^{\perp}\cap V_X^{\perp}=S\cap V_X=S_X$, so
$\dim(C\cap V_X)=\dim C+N-\dim(C+V_X)=(N+k)+N-(2N-\dim S_X)=k+\dim S_X$, and
$\dim\bar X=\dim(C\cap V_X)-\dim S_X=k$. Dually for $Z$. The two display lines combine
this with J.1/J.2 and Theorem G(ii). $\square$

*Machine checks.* 299/299 random non-CSS stabilizer codes (`span`-built check matrices):
$\dim\bar X=k$ exactly. All 4 target rows and **all 368 catalogue rows**: the four J.3
identities hold on every row (`results/partial_runs/xsector/crosstable.json`,
`catalogue_scan.json`; the scan also independently re-verifies Theorem G(ii)'s
$k_Q=k_P-\dim\bar\Delta$ on 368/368).

**Repair of Theorem G(ii).** The published proof's "X-sector unchanged" step is false as
stated — $\mathrm{Xcen}$ itself shrinks by $\rho_X\ge\bar\Delta$ (`phase2_58`: $40\to12$).
The **conclusion** of G(ii) survives verbatim: by J.3 the pure-X quotient keeps dimension
$k_Q$ because the pure-X **stabilizer** shrinks simultaneously by exactly $\rho_X-\bar\Delta$,
and the Z-side identity $\dim\bar Z(Q)=\dim\ker[\,A\;B\,]-\dim(S_Z+\Delta)=k_P-\bar\Delta=k_Q$
*is* Theorem G(ii) restated. The correct replacement sentence for the paper:
*"both the pure-X centralizer and the pure-X stabilizer shrink (by $\rho_X$ and
$\rho_X-\bar\Delta$ respectively); the pure-X logical dimension therefore stays $k_Q$, and
the full drop $\bar\Delta$ is accounted for on the Z-side."*

**On the identity "$\dim\bar X(Q)+\dim\bar Z(Q)=k_Q$".** As stated it is **false in
general**: each sector has dimension exactly $k_Q$, so the sum is $2k_Q$ (already the CSS
parent gives $\dim\bar X(P)+\dim\bar Z(P)=2k_P$). The correct, cross-checked statement is
$\dim\bar X(Q)=\dim\bar Z(Q)=k_Q=k_P-\bar\Delta$, verified on all four rows below.

**(J.4) Decomposition of the pure-X logicals; survival and demotion.**
$\mathrm{Xcen}(Q)\subseteq\mathrm{Xcen}(P)$ splits the X-logicals of $Q$ into two disjoint
classes:

- **Survivors** — parent X-logicals with $[\,C\;D\,]x^{\mathsf T}=0$. Each has weight $\ge d_X(P)$.
- **Demoted** — $x\in\operatorname{rowspace}H_X\cap\ker[C\;D]$ with $x\notin S_X(Q)$:
  parent X-*stabilizers* that stop being stabilizers of $Q$.

Hence $\ d_X(Q)=\min(\text{survivor min},\ w_{\rm dem})$ with $w_{\rm dem}$ the minimum
demoted weight ($\infty$ if none), and:

**(i) X-survival criterion [proved].** *If some minimum-weight X-logical $x^*$ of $P$
satisfies $[C\;D]x^{*\mathsf T}=0$, then $d_Q\le d_X(Q)\le d_X(P)=d_Z(P)$ (J.0).* This is
the exact mirror of Theorem G(iii), logically independent of it, and caps $d_Q$ by a
pure-X mechanism Theorem H never sees.

**(ii) Demotion pairing [proved].** *$x=\lambda[\,A\;B\,]$ is demoted $\iff$
$z:=\lambda[\,C\;D\,]\in\ker[\,A\;B\,]\smallsetminus(S_Z+\Delta)$, i.e. $z$ is a
Z-logical of $P$ **and** of $Q$; moreover $[(x|0)]=[(0|z)]$ in $L(Q)$ (their sum is the
row $r_\lambda$), so the demoted class has both a pure-X and a pure-Z representative.*
In particular $\mathrm{wt}(z)\ge d_Z(P)$.

*Proof.* $[\,A\;B\,]z^{\mathsf T}=[\,A\;B\,][C\;D\,]^{\mathsf T}\lambda^{\mathsf T}
=M\lambda^{\mathsf T}=[C\;D\,]x^{\mathsf T}=0$; nonmembership is J.2. $\square$

**(iii) Single-row barrier [proved].** For $\lambda=e_i$ (a single first-block row):
$\mathrm{wt}(x)=\mathrm{wt}(A)+\mathrm{wt}(B)$ and $\mathrm{wt}(z)=\mathrm{wt}(C)+\mathrm{wt}(D)=\delta$,
so a single-row demotion below $d_X(P)$ requires
$$\mathrm{wt}(A)+\mathrm{wt}(B)\ <\ d_X(P)\ \le\ d_Z(P)\ \le\ \delta .$$
Catalogue-style light perturbations ($\delta\le6$, rows of weight $\ge6$) can never demote
a single row at all.

**(J.6) Rank cost of killing the X-side [proved].** Let $M_X(P)$ be the span of the
translation orbits of all minimum-weight X-logicals of $P$ and $T_X(P)=\dim(M_X+S_X)/S_X$.
If $d_Q>d_Z(P)=d_X(P)$ then every minimum-weight X-logical is killed or demoted-heavy;
in particular the survivor set is empty of minima, so $M_X\cap\mathrm{Xcen}(Q)\subseteq S_X(P)$,
hence
$$ T_X(P)\ \le\ \rho_X .$$
*Proof.* $\dim(M_X\cap\mathrm{Xcen}(Q))\ge\dim M_X-\rho_X$ (codimension of
$\mathrm{Xcen}(Q)$ in $\mathrm{Xcen}(P)$ is $\rho_X$); if the intersection lies in
$S_X(P)$ then $\dim M_X-\rho_X\le\dim(M_X\cap S_X(P))$, i.e. $T_X(P)\le\rho_X$. $\square$

---

## 3. Verdicts

**Direction A — "is $d_X(Q)\ge d_X(P)$ always?"**

By J.4, $d_X(Q)<d_X(P)$ **iff** some demoted parent X-stabilizer has weight $<d_X(P)$.
The unconditional statement is **open**: no proof is known, and no counterexample exists
in any computation performed:

| search | lattices ($\ell m$) | perturbations | valid | demotion instances | violations |
|---|---|---|---|---|---|
| lattice hunt 1 (≤3 monomials, menu parents) | 6–12 | ~155k | 3,228 | 0 | 0 |
| lattice hunt 2 (δ ≤ 5–6, menu parents) | 6–12 | 16,087,835 | 1,860,331 | 920,678 | **0** |
| lattice hunt 3 (**complete** $(C,D)$ space) | 4, 6, 8 (complete); 9 (δ ≤ 7) | 19,682,092 | 4,603,076 | 2,925,746 | **0** |
| EXP-044 catalogue parents (valid space $V$ sampled; §6) | 36–144 (**all** $n\le144$ parents) | 5,483 checked | 5,483 | 2,334 + 1,675 no-demotion | **0** |
| catalogue rows as published (all 368) | 72–360 | 368 | 368 | **255** rows with demotion | **0** |

(For $\ell m\in\{4,6,8\}$ the perturbation space is **complete**: every nonempty
$(C,D)$ pair of monomial sets up to simultaneous translation, any weight; parents: all
$A,B$ of weight ≤3 — at $\ell m\le6$ that is every parent — plus menu shapes at $\ell m=8$.
Catalogue comparison: $d_X(Q)$ exact on 276 of 368 rows and $d_Z(P)=d_X(P)$ certified on
279; **221 rows carry both quantities, with 0 violations** of $d_X(Q)\ge d_X(P)$ (17 at
equality). All 7 certified reversals additionally obey $d_X(Q)>d_Z(P)$
($12,12,36$ exact on `phase2_58/60/88`; $24>8$ exact on `9_6_0183`, `12_6_0217`;
automatic on `phase2_71/72` from $d_Q>d_Z(P)$).

*Correction note (this session).* The validity filter in the first hunt run used
numpy-uint8 matrix products, i.e. counts mod 256 XOR-ed instead of GF(2) products; the
effect was to **falsely reject** some truly valid $(C,D)$ (e.g. count pairs $(0,2)$,
$(1,3)$ in $AC^{\mathsf T}+BD^{\mathsf T}$). All hunts above were re-run with exact
`gf2.linalg.matmul` after the fix (valid counts roughly doubled); no accepted instance
was ever invalid, and no demotion or distance computation was affected — only the
sampled subset was distorted, which the re-runs remove. Artifacts:
`results/partial_runs/xsector/hunt{,2,3}*.json`, `catalogue_scan.json`,
`results/processed/exp044_xsector_collapse_hunt.json`.)

Proved regimes, both GF(2)-checkable per row:

- **[proved]** $d_X(Q)\ge d_X(P)$ whenever $d_Q>d_Z(P)$ (then $d_X(Q)\ge d_Q>d_Z(P)=d_X(P)$) — vacuous on reversals;
- **[proved]** $d_X(Q)\ge d_X(P)$ whenever there is no demotion below $d_X(P)$ (in particular when the demotion space is $0$; holds on `phase2_58/60/88`);
- **[proved]** a violation forces the weight flip $\mathrm{wt}(\lambda[C\;D\,])\ge d_Z(P)>\mathrm{wt}(\lambda[\,A\;B\,])$ (J.4(ii)); single-row violations force $\delta\ge d_Z(P)$ (J.4(iii)).

**Overall verdict: [no counterexample found — open in general; holds in every computed instance, including exhaustively at $\ell m\le8$].**

**Direction B — "is $d_X(Q)<d_X(P)$ possible?"**

**[open — not exhibited].** The event is exactly a light demotion (J.4); every candidate
mechanism is pinched by the proved necessary conditions above. On all computed instances
the empirical law is sharp: $\min_{\text{instances}}(w_{\rm dem}-d_X(P))=0$ over both
hunts and the catalogue (17 catalogue rows sit at $d_X(Q)=d_X(P)$ exactly) — equality
occurs, violation never.

---

## 4. Numerical cross-table (computed this session)

Dims via GF(2) rank/nullspace on rebuilt $A,B,C,D$ (`results/partial_runs/xsector/crosstable.json`);
$T$, $d_Z(P)$ from EXP-039 parent certificates (`results/partial_runs/exp039/parent_*.json`, keys
`T`, `T_is_exact`, `d_z_parent`, `d_z_exact`); $d_X(P)=d_Z(P)$ by J.0; $d_X(Q)$ exact by
enumerating $\mathrm{Xcen}(Q)/S_X(Q)$ ($2^{\dim\mathrm{Xcen}(Q)}$ cosets).

| row | $\ell,m$ | $n$ | $k_P$ | $k_Q$ | $\bar\Delta$ | $T$ | $\mathrm{Xcen}(P)$ | $\mathrm{Xcen}(Q)$ | $\rho_X$ | $S_X(P)$ | $S_X(Q)$ | demotion | $d_X(P)$ | $d_X(Q)$ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `phase2_58` | 6,6 | 72 | 8 | 4 | 4 | 4 | **40** | **12** | 28 | 32 | 8 | 0 | 4 (cert.) | 12 (exact) |
| `phase2_60` | 6,6 | 72 | 8 | 4 | 4 | 4 | **40** | **12** | 28 | 32 | 8 | 0 | 4 (cert.) | 12 (exact) |
| `phase2_88` | 9,6 | 108 | 4 | 2 | 2 | 2 | 56 | 2 | 54 | 52 | 0 | 0 | 4 (cert.) | 36 (exact) |
| `12_6_0193` | 12,6 | 144 | 12 | 12 | **0** | 12 | 78 | 32 | 46 | 66 | 20 | 6 | 12 (cert.) | ≥12 (UNVERIFIED above; $\dim\mathrm{Xcen}(Q)=32>22$) |

Consistency checks, all passing on all four rows:

1. **HostileReview's 40→12 confirmed**: $\dim\mathrm{Xcen}(P)=40\to\dim\mathrm{Xcen}(Q)=12$ on `phase2_58`/`60` — a shrink of $\rho_X=28\ge\bar\Delta=4$.
2. $k_Q=k_P-\bar\Delta$ (Theorem G(ii)) — 4/4.
3. $\dim\bar X(Q)=\dim\bar Z(Q)=k_Q$ (J.3) — 4/4 (e.g. `phase2_58`: $12-8=4=4=k_Q$; `12_6_0193`: $32-20=12$ with $\bar\Delta=0$).
4. $S_X$ shrink $=\rho_X-\bar\Delta$ — 4/4 (`phase2_58`: $32-8=24=28-4$; `12_6_0193`: $66-20=46=46-0$).
5. Forced saturation (Theorem I) cross-check: all three reversals sit on the law
   ($\bar\Delta=T$: 4, 4, 2 with $T=k_P/2$); `12_6_0193`'s parent is family-closed ($T=12=k_P$)
   and its perturbation pays nothing ($\bar\Delta=0$, $k_Q=k_P$).
6. The X-side shrink does **not** track the Z-side loss: $\rho_X$ spans $28\to54$ while
   $\bar\Delta\in\{0,2,4\}$; `12_6_0193` loses 46 X-centralizer dimensions with **zero**
   logical qubits lost. Catalogue-wide, $\rho_X$ reaches 176 (`30_6_0273`,
   $\bar\Delta=2$) and the demotion space dimension reaches 20.
7. Reversal rows must and do satisfy $d_X(Q)>d_Z(P)$: 12, 12, 36 > 4 on the three rows above,
   and $24 > 8$ exact on `9_6_0183` and `12_6_0217` (catalogue scan) — the perturbations
   killed **every** light parent X-logical (and demoted nothing light), consistent with J.4(i).

`12_6_0193` is the sharp display case: demotion space has dimension 6 (so demoted
X-logicals of $Q$ **exist**), yet $d_Q=12$ certified (paper §4) forces every one of them
to weight $\ge12=d_X(P)$ — demotion without violation.

---

## 5. Implications for the paper's $d_Q>d_Z(P)$ framing

The no-go theorems G/H/I, read literally, never consult the X-sector: the survival
criterion is pure-Z, the absorption mechanism is pure-Z, and the closure corollary
($T=k_P\Rightarrow$ family closed) follows from the Z-side plus the dimension identity
alone. **A full-distance no-go in the paper's current formulation therefore does not
require the X-side — with one exception: the proof of the dimension identity G(ii)
itself.** That proof currently contains an incorrect step ("X-sector unchanged");
Theorem J.3 is the drop-in repair, with no new hypotheses, and it shows the repair is
pure bookkeeping: centralizer and stabilizer shrink by $\rho_X$ and $\rho_X-\bar\Delta$,
the pure-X quotient keeps dimension $k_Q$, and the whole loss $\bar\Delta$ is accounted
for on the Z-side. Recommendation: patch G(ii)'s proof with the J.3 sentence and cite
this note; nothing else in the paper changes.

The X-side becomes *load-bearing* the moment one wants more than pricing — i.e. a
**per-row domination engine** complementary to Theorem H, or a full-distance no-go on
parents with $T(P)<k_P$ (where H permits increases, and the 7 certified reversals prove
increases happen). There the exact missing lemma is the **parent-level X-survival
statement**:

> **(Missing Lemma X)** For a BB parent $P$, exhibit a parent-only quantity
> $T_X(P)$ (analogue of $T(P)$: $\dim(M_X+S_X)/S_X$ over the translation module $M_X$ of
> minimum-weight X-logicals) and prove: *for every valid $[C\;D]$ over $P$, either some
> minimum-weight Z-logical survives (Theorem H's premise fails, capping $d_Q\le d_Z(P)$)
> or some minimum-weight X-logical $x^*$ has $[C\;D]x^{*\mathsf T}=0$ (capping
> $d_Q\le d_X(P)=d_Z(P)$ by J.4(i)).*

Two caveats make this the honest frontier rather than a formality. First, the
unconditional form ("some X-logical always survives") is **false**: the reversals
`phase2_58`/`60`/`88` kill every minimum-weight X-logical ($d_X(Q)=12,12,36\gg4$) while
absorbing every minimum-weight Z-logical — the disjunction must be proved, not either
disjunct. Second, the module argument that powers Theorem H does not transplant: Z-side
absorption is a *subspace* condition on the logicals themselves, while X-side survival
is *kernel avoidance* ($x^*\notin\ker[C\;D]$), and sums of minimum-weight X-logicals are
not minimum-weight; the only parent-level fragment that survives is the rank bound
$T_X(P)\le\rho_X$ (J.6), which is valid but weak. Settling Lemma X — or refuting it with
a computed lattice — is the precise open problem standing between the paper's
price-theorems and a genuine full-distance no-go over general BB parents.

---

## 6. EXP-044: parent-level collapse hunt over the catalogue (follow-up)

**Question (J.5).** Does *any* valid perturbation of *any* catalogue parent exhibit
$d_X(Q)<d_X(P)$? Artifact: `results/processed/exp044_xsector_collapse_hunt.json`
(schema `exp044-xsector-collapse-hunt-v1`); script
`results/partial_runs/xsector/exp044_collapse_hunt.py`.

**Method.** Valid perturbations form the linear space
$V=\ker\nu$, $\nu(C,D)=\mathrm{skew}(AC^{\mathsf T}+BD^{\mathsf T})$; measured over all
202 parents, $\dim V\in[29,282]$ ($0.78\,\ell m$ at the large end), so **full
enumeration of $V$ is infeasible at every parent** ($2^{29}$ already at $n=36$) — this
is the honest reason the budget is sampled, not enumerated. Per parent the run checks:
the catalogue's own $(C,D)$ (≤6 per parent), uniform random elements of $V$, sums of
≤3 basis vectors of $V$, and — for $\ell m\le72$ parents only — rejection-sampled light
$(C,D)$ of weight ≤8. Each instance gets the exact J.4 demotion check: full
$2^{\dim\ker M}$ span enumeration when $\dim\ker M\le18$, else a 4096-sample probe
(recorded as *capped*, never as safe). The decrease threshold is the **certified**
$d_Z(P)$ from EXP-039 (all 116 covered parents carry an exact certificate), sound
because $d_X(P)=d_Z(P)$ (J.0).

**Budget and outcome (26-minute single-thread cap).** 116 of 202 parents covered —
exactly **all** parents with $n\le144$ ($n=36$: 3, $n=72$: 8, $n=108$: 68, $n=144$: 37);
the 86 parents with $n\ge180$ were not reached (per-instance exact GF(2) chains cost
0.1–0.5 s at those sizes; the intended ≥10k instances would have needed roughly double
the cap — stated here rather than inflated). 5,483 valid instances checked:
1,675 with no demotion at all, 2,334 with exact $w_{\rm dem}\ge d_Z(P)$,
1,474 capped-unresolved (large $\ker M$; probe found nothing), **0 decrease witnesses**.
Across the 109 parents with an exact $w_{\rm dem}$, $\min(w_{\rm dem}-d_Z(P))=0$ —
equality attained, violation never.

**Verdict [computed]: `X_COLLAPSE_NEVER_OBSERVED_WITHIN_COMPUTED_BUDGET`.**
No X-side distance-decrease witness exists among (i) all 368 published catalogue rows,
(ii) 6,466,635 valid perturbations on small lattices including the complete spaces at
$\ell m\le8$, and (iii) 5,483 sampled valid perturbations across every $n\le144$
catalogue parent. The unconditional statement remains open (§3); the parent-scale
evidence base now covers the entire $n\le144$ half of the catalogue.

---

### Provenance

- Algebra + cross-table: `results/partial_runs/xsector/verify_algebra.py` → `crosstable.json` (this session).
- Catalogue scan: `results/partial_runs/xsector/catalogue_scan.py` → `catalogue_scan.json` (this session).
- Lattice hunts (validity filter corrected to GF(2) matmul, re-run): `hunt_demotion.py`, `hunt2.py`, `hunt3.py` → `hunt.json`, `hunt2.json`, `hunt3.json` (this session).
- EXP-044 parent-level collapse hunt: `results/partial_runs/xsector/exp044_collapse_hunt.py` → `results/processed/exp044_xsector_collapse_hunt.json` (this session).
- Parent certificates: `results/partial_runs/exp039/parent_*.json` and pool `results/partial_runs/exp037/` (read only).
- Random-stabilizer identity check and demotion-pairing spot checks: inline, this session.
