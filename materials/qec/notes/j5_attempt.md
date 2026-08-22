# J.5 attempt — outcome: J.5 is FALSE (witnesses verified), plus machine-verified fragments

**XSideTheory, 2026-08-19.** Companion to `notes/theorem_j_xsector.md` (J.0–J.7).
Question: is $d_X(Q)\ge d_X(P)$ for every BB parent and every valid perturbation $[C\;D]$
(Direction A = J.5)? **No.** An explicit, mechanically verified counterexample family is
exhibited below (flagship: parent of `12_6_0193`), the obstruction mechanism is identified
and certified on small lattices, and the sharp hunt is ported to all 46 $n=180$ parents.

All numbers quoted are read from machine-generated artifacts:
- `results/partial_runs/xsector/j5_witness_12_6_0193.json` — the flagship counterexample;
- `results/partial_runs/xsector/j5_small_certificates.json` — Route A/B machine lemmas
  (lattices $3\times3$, $3\times5$, $4\times4$);
- `results/partial_runs/xsector/exp044_n180_hunt.json` — Route C census ($n=180$ primary,
  $n=144/360$ bonus);
- `results/partial_runs/xsector/j5_syz_census.json` — syzygy-family demotion census;
- scripts: `j5core.py`, `j5_verify_witness.py`, `j5_reproduce_witness.py`,
  `exp044_n180_hunt.py`, `j5_small_certificates.py`, `j5_syz_census.py`
  (same directory; GF(2)/numpy only, no SAT/Stim).

---

## 1. Theorem J.5 candidate statement (as attacked), and verdict

> **J.5 (candidate, FALSE in general).** For every BB parent $P=(A,B)$ and every valid
> perturbation $(C,D)$ (i.e. $M=AC^{\mathsf T}+BD^{\mathsf T}$ symmetric),
> $$d_X(Q)\ \ge\ d_X(P)\ (=d_Z(P)\ \text{by J.0}).$$

**Verdict: refuted.** §2 gives the counterexample; §3 records exactly which salvageable
fragments were proved and machine-verified; §4 (Route B) the rank-bound tightness data;
§5 (Route C) the parent-level census at $n\in\{144,180,360\}$.

---

## 2. The counterexample (Route A outcome)

**Flagship parent** `12_6_0193` ($\ell m=72$, $k_P=12$, $T=12$, $d_Z(P)=12$ certified
exact by EXP-039, fingerprint `4c7eb964a6a3e38c…dac50`):

```
A = x y^2 + x^10 y^3 + x^10 y^4      B = 1 + x y^5 + x^11 y^4
C = (1+y+y^4+y^5)(1 + x^3 + x^6 + x^9)(1 + y^2 x^3 + x^? ...)   (12 terms, see artifact)
D = 0
```

Take the perturbation
$$C=(1+y+y^4+y^5)\,(1+x^3+x^6+x^9)\,(1+x^3y^2),\qquad D=0.$$
(Exact term list in `j5_witness_12_6_0193.json`: `C_terms` has 12 monomials, `D_terms=[]`.
Reproduction: `python results/partial_runs/xsector/j5_verify_witness.py`.)

**Claim (machine-verified, two independent code paths, saved to the witness artifact):**
1. $M=AC^{\math T}+BD^{\math T}=0$ — hence valid (and indeed maximally degenerate).
2. $x_0=$ row $0$ of $[A\;B]$, $\mathrm{wt}(x_0)=6$: $x_0\in\mathrm{Xcen}(Q)$
   ($H_Zx_0^{\math T}=0$, $[C\;D]x_0^{\math T}=M^{\math T}e_0^{\math T}=0$),
   $x_0\in\operatorname{rowspace}H_X$ by construction.
3. $z_0:=e_0[C\;D]$, $\mathrm{wt}(z_0)=12$, satisfies $z_0\in\ker H_X$ and
   $z_0\notin S_Z+\Delta$ (path A: rank test; path B: EXP-044's $S_X(Q)$ span test;
   a third path — direct GF(2) solve for $x_0\in S_X(Q)$ — is implemented in
   `j5_reproduce_witness.py`). So $x_0$ is a **demoted** parent X-stabilizer: a genuine
   X-logical of $Q$ of weight 6.
4. $d_X(Q)\le 6<12=d_Z(P)=d_X(P)$. The perturbed code has $\dim\bar\Delta=0$, $k_Q=12$
   (no logical qubits lost), $\rho_X=52$.

**Why previous hunts missed it.** The violating family is the subspace
$\mathrm{Syz}:=\{(C,D):M=0\}\subset V$ of the validity space. At this parent
$\dim V=112$ and $\dim\mathrm{Syz}=78+? \to 78$, of which the *trivial* part
$\{(B^{\math T}g,A^{\math T}g)\}$ has dimension $66$, so the nontrivial part has relative
density $2^{-112}\cdot(\text{over}\,V)=2^{78-112}=2^{-34}$ — uniform $V$-sampling
(EXP-044's method) cannot see it. The sharp probe interrogates $\mathrm{Syz}$ directly.

---

## 3. Route A — proof skeleton, with each lemma's status

The proof attempt collapsed into the following verified fragments. Here "machine-verified"
means the statement is computed and cross-checked per instance in the listed artifact.

### L1 — Z-monotonicity **[PROVED]**
For every valid perturbation, $\mathrm{Zcen}(Q)=\ker_{\mathrm{col}}H_X=\mathrm{Zcen}(P)$
(pure-Z vectors pair only with the X-parts of check rows; the perturbation modifies only
Z-parts of first-block rows, and second-block rows are pure-Z). Hence
$$d_Z(Q)\ \ge\ d_Z(P)\qquad\text{always.}$$
*Proof.* $S_Z(Q)=S_Z+\Delta\supseteq S_Z$ inside the *fixed* centralizer
$\ker_{\mathrm{col}}H_X$: every representative of a nonzero class of
$\ker H_X/(S_Z+\Delta)$ is a genuine parent Z-logical, weight $\ge d_Z(P)$. $\square$
Machine note: the asymmetric counterpart $\mathrm{Xcen}(Q)=\mathrm{Xcen}(P)\cap\ker[CD]$
is J.1; the X-side has no monotonicity precisely because Xcen *shrinks* while $S_X$
*shrinks too* — and §2 shows the shrink can trade.

### L2 — violation window **[PROVED, definitional]**
$d_X(Q)<d_X(P)$ iff some demoted $x=\lambda[A\;B]$ has
$$d_{S_X}(P)\ \le\ \mathrm{wt}(x)\ <\ d_Z(P),$$
where $d_{S_X}(P)=\min\mathrm{wt}(\operatorname{rowspace}H_X\setminus\{0\})$.
(The demoted set is $W\setminus S_X(Q)$ with $W=\operatorname{rowspace}H_X\cap\ker[CD]$;
J.4.) At $n=180$ catalogue parents: row weight $6$–$8$ and $d_Z(P)\in\{6,\dots,16\}$, so
the window is wide open on most rows.

### L3 — sharp single-row barrier **[PROVED; machine-verified]**
$M$ is a polynomial (block-circulant) matrix, so all its columns are $R$-translates of
column 0: a single-row demotion ($\lambda=e_i$, requiring column $i$ of $M$ to vanish)
exists **only if $M=0$**, in which case *every* row is a demotion candidate.
This strictly sharpens J.4(iii). *Machine check:*
`j5_small_certificates.json` → `checks.no_zero_col_unless_M0` $=151{,}805/151{,}805$
valid instances (zero zero-columns over nonzero $M$), `checks.circulant_cols` verified
on the same set.

### L4 — the $\mathbf{M=0}$ (syzygy) family **[PROVED; machine-verified]**
$$\dim\{(C,D):M=0\}\ =\ 2\ell m-\operatorname{rank}H_Z\ =\ k_P+\operatorname{rank}H_X .$$
(*Proof.* $(C,D)\mapsto M$ has the same rank as $[J A\;|\;J B]\sim H_Z$. $\square$)
It contains the **trivial** subfamily $(C,D)=(B^{\math T}g,A^{\math T}g)$ of dimension
$\operatorname{rank}H_Z$, for which $Q=P$ exactly (row operations remove the Z-part:
$M=(AB+BA)g^{\math T}$-style cancellation; machine-verified: trivial family demotes
nothing on every instance tested). The nontrivial **excess**
$$\sigma\ =\ \dim\{M=0\}-\operatorname{rank}H_Z\ =\ k_P \qquad\text{whenever }
\operatorname{rank}H_X=\operatorname{rank}H_Z$$
(verified `sigma_eq_kP` on 900/900 small-lattice parents; holds on all checked
$n=144/180$ catalogue parents). For $M=0$: $W=\operatorname{rowspace}H_X$ *entirely*
(verified `M0_implies_W_full` on 36{,}895 instances), and a single row demotes iff the
coefficient pair $(c,d)$ has nonzero class in $\ker H_X/(S_Z+\Delta)$.

### L5 — reversal obstruction (why the "obvious" L1-transfer fails) **[PROVED, machine-checked]**
The natural proof attempt $d_X(Q)=d_Z(Q^{\mathrm{rev}})$ via the weight-preserving
involution $\rho=\sigma\circ\tau_0$ (J.0) fails: $\rho\,Q$ has check matrix
$[\,[A\;B\,|\,0\,],\,[B^{\math T}\;A^{\math T}\,]\;;\;[D^{\math T}\,C^{\math T}\,|\,A^{\math T}\,B^{\math T}]\,]$
which is **not** PBB — its pure-Z centralizer is
$\ker H_X\cap\ker[D^{\math T}C^{\math T}]\subsetneq\ker H_X$ in general, so L1's
argument does not transfer. The X-collapse has no Z-side mirror argument because the PBB
ansatz is X/Z-asymmetric by construction.

### The missing step (why Route A can never prove J.5)
J.4(ii) bounds only the *partner* weight $\mathrm{wt}(z)\ge d_Z(P)$; the demoted class
pairs an X-representative in the parent's *stabilizer* space (abundant at low weight)
with a Z-representative in the parent's *logical* space (bounded below by $d_Z(P)$), and
the pairing carries no weight inequality for the X-side. §2 shows the inequality fails.
**Route A is closed: the statement is false; the obstruction is the syzygy class escape
described in L4.**

---

## 4. Route B — rank identity + J.6 tightness certificates (small lattices)

The dimension identity $\dim\mathrm{Xcen}(Q)=\dim\mathrm{Xcen}(P)-\rho_X$ and J.6's
$T_X(P)\le\rho_X$ (for $d_Q>d_Z(P)$) cannot force min-weight non-descent — §2's witness
already shows non-descent itself is false, and the small-lattice data isolates exactly
where dimension-counting stops:

- **Tightness exists.** 25+ instances (capped) at $3\times3$ and $3\times5$ with
  $T_X=\rho_X$ exactly — the J.6 inequality cannot be strengthened to strict inequality.
  (Records in `j5_small_certificates.json` → `lattices[].tight_TX`.)
- **Small lattices are demotion-free.** At $\ell m\in\{9,15,16\}$: 300 parents per
  lattice, 151,805 valid perturbations examined (syzygy elements included),
  **zero demotions at all** — consistent with, and strengthening of, the
  exhaustive $\ell m\le8$ hunts (which found 0 violations in 2.9M demotion instances).
  Mechanism: at these sizes $S_Z+\Delta=\ker H_X$ for the sampled perturbations
  ($\dim\bar\Delta=k_P$, full absorption), so the demotion space is empty.
- **Collapse boundary.** The first collapses observed are at $\ell m=72$ ($n=144$:
  flagship) and $\ell m=90$ ($n=180$: ubiquitous). The open window requirement
  $d_{S_X}(P)<d_Z(P)$ and a nontrivial syzygy class escaping $S_Z+\Delta$ are both
  necessary; small lattices satisfy the first but never the second in any sampled
  instance.

---

## 5. Route C — census at $n\in\{72..360\}$, harvested

Artifact `results/partial_runs/xsector/exp044_n180_hunt.json` (schema
`exp044-n180-xcollapse-hunt-v1`, harvested by Main from the completed hunt log on
2026-08-19; the original run finished all 116 parents but died before JSON harvest).
102 rows parsed (14 log lines truncated/unparseable; counts below are over the parsed
102). Method deltas as designed above; bounds are EXP-039 exact certificates where
present, `None` otherwise (uncertified parents cannot score a collapse).

**Headline: the collapse is ubiquitous, not exceptional.**

| scope | parents | collapse-positive | violation-free |
|---|---|---|---|
| $n=180$ (primary) | 42 | 20 | 22 — **all with no distance bound** |
| $n=144$ (bonus) | 28 | **28 (100 %)** | 0 |
| $n=360$ (bonus) | 32 | 1 | 31 — all unbounded |
| bounded parents only | 49 | **49 (100 %)** | 0 |

- Every collapse found sits at demoted weight **exactly 6** (the trinomial check-row
  weight); collapse-positive $d_{\rm lb}\in\{8,10,12,14,16\}$ (histogram
  $\{8{:}10,\,10{:}6,\,12{:}25,\,14{:}7,\,16{:}1\}$).
- Every parent that carries both a certified bound and weight-6 check rows exhibits
  the collapse. The 53 violation-free rows are *all* uncertified parents; the census
  measures nothing there (no bound ⇒ no violation can be scored).
- Witness reproduction is deterministic: `j5_reproduce_witness.py`; four full-term-list
  self-verifying witnesses are stored alongside (`j5_witness_{12_6_0193,15_6_0219,9_6_0136,6_6_0099}.json`).

---

## 6. Consequences for the PBB no-go program

1. **Every PBB candidate from a trinomial parent needs its own two-sided distance
   certificate before downstream use.** A sibling perturbation of the flagship family
   is $\le[144,12,6]$ at full $k$; the same hazard exists on *every* bounded parent
   probed (49/49).
2. The Z-side of the program is untouched: L1 monotonicity and Theorems G/H/I caps
   all carry through; the no-go layers are consistent (no bounded parent escapes the
   CSS envelope through this channel — the collapse runs *downward*, not upward).
3. The search-side lesson: validity-space sampling cannot see the syzygy family
   (density $\sim2^{-34}$ on the flagship); any claim of "no collapse observed" must
   probe $\mathrm{Syz}$ directly, as this census does.
4. Practical check available now: for any concrete $(A,B,C,D)$,
   `results/partial_runs/xsector/j5_verify_witness.py` reproduces the decision in
   $\sim$1 s.

## 7. Upgrade (2026-08-19, Main): Theorem J-C — the deterministic constructor

The sampled census of §5 is superseded by the **deterministic nullspace constructor** of
`notes/theorem_jc_syzygy_constructor.md` (EXP-046, artifact
`results/processed/exp046_generic_syzygy.json`): the collapse mechanism is decidable in
polynomial time ($c\in\ker A$, $D=0$; validity and partner membership automatic;
demotion = one GF(2) rank test). Verdict over all 202 catalogue parents: **192
deterministic demotions, 99/100 bounded open-window collapses, 98/99 k-preserving**;
boundary instance `9_6_0175` (row-overlap $\le2$ blocks multi-row channels; single-row
channel immune through syzygy basis+pairs depth). Zero mismatches against this note's
census (all 49 census-positives inside the EXP-046 set).
