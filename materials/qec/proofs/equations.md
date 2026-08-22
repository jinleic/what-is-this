# Equations of the PBB structure theory

Every equation below is numbered, stated in the exact form used by
`proofs/pbb_structure.md`, and **machine-verified** by
`experiments/exp021_verify_equations.py`.  The verifier reports PASS/FAIL per
equation over the seven Bravyi BB instances, the 368-code published PBB
catalogue, and randomly generated codes, and it includes **negative controls**
that must fail — a check that cannot fail proves nothing.

Epistemic tags: **[proved]** analytic proof given in `pbb_structure.md`;
**[exact]** verified by exact GF(2) arithmetic or a solver returning proven
`OPTIMAL`/`INFEASIBLE`.

---

## Part I — Conventions

**(E1) Symplectic form.**  For $v=(x\mid z),\;v'=(x'\mid z')\in\mathbb F_2^{2n}$,

$$\langle v,v'\rangle_s \;=\; x\cdot z'^{\mathsf T} \;+\; z\cdot x'^{\mathsf T} \pmod 2 .$$

**(E2) Stabiliser validity.**  A check matrix $H=(H_x\mid H_z)$ is valid iff

$$H_x H_z^{\mathsf T} \;+\; H_z H_x^{\mathsf T} \;=\; 0 \pmod 2 .$$

**(E3) Code dimension and distance.**

$$k \;=\; n-\operatorname{rank}_{\mathbb F_2}(H),\qquad
d \;=\; \min\{\, \mathrm{wt}_s(v) \;:\; v\in S^{\perp_s}\setminus S \,\},$$

where $\mathrm{wt}_s(v)=\#\{j: (x_j,z_j)\neq(0,0)\}$.

**(E4) Ring and regular representation.**

$$R=\mathbb F_2[x,y]/(x^{\ell}-1,\,y^{m}-1),\qquad \dim=\ell m,\qquad n=2\dim,$$
$$M(x^ay^b)_{g,\,g+(a,b)}=1,\qquad g\in G=\mathbb Z_\ell\times\mathbb Z_m .$$

**(E5) Ring facts.**

$$\textbf{(R1)}\quad M(f)M(h)=M(h)M(f),
\qquad\qquad
\textbf{(R2)}\quad M(f)^{\mathsf T}=M(f^{*}),\;\; f^{*}(x,y)=f(x^{-1},y^{-1}).$$

---

## Part II — Bivariate-bicycle (CSS) codes

**(E6) Construction.**

$$H_X=[\,A\;\;B\,],\qquad H_Z=[\,B^{\mathsf T}\;\;A^{\mathsf T}\,],\qquad n=2\ell m .$$

**(E7) CSS commutation.** **[proved]**  By **(R1)**,

$$H_X H_Z^{\mathsf T} \;=\; AB+BA \;=\; 0 .$$

**(E8) Dimension.** **[exact]**

$$k\big(\mathrm{BB}(A,B)\big) \;=\; n-\operatorname{rank}H_X-\operatorname{rank}H_Z
\;=\; 2\big(\ell m-\operatorname{rank}[\,A\;\;B\,]\big).$$

**(E9) Self-duality (Proposition 4).** **[proved]**  Let $\sigma$ permute qubits by

$$\sigma:\;(\mathrm L,g)\mapsto(\mathrm R,-g),\qquad (\mathrm R,g)\mapsto(\mathrm L,-g).$$

Then $\sigma\big(\operatorname{rowspace}H_X\big)=\operatorname{rowspace}H_Z$, hence

$$d_X\big(\mathrm{BB}(A,B)\big)\;=\;d_Z\big(\mathrm{BB}(A,B)\big).$$

---

## Part III — Perturbed bivariate-bicycle codes

**(E10) Construction.**  With $P=[\,C\;\;D\,]$,

$$H \;=\; \begin{pmatrix} H_X & P\\ 0 & H_Z\end{pmatrix},
\qquad
H_x=\begin{pmatrix}A&B\\0&0\end{pmatrix},\quad
H_z=\begin{pmatrix}C&D\\ B^{\mathsf T}&A^{\mathsf T}\end{pmatrix}.$$

$C=D=0$ recovers the parent $\mathrm{BB}(A,B)$.

**(E11) Lemma 0 — validity.** **[proved]**

$$H_xH_z^{\mathsf T}=\begin{pmatrix}M&0\\0&0\end{pmatrix},
\qquad M:=AC^{\mathsf T}+BD^{\mathsf T},$$

so the rows of $H$ pairwise commute **iff**

$$\boxed{\;M=M^{\mathsf T}\;}\qquad\text{equivalently}\qquad M+M^{\mathsf T}=0 .$$

In polynomial form $ac^{*}+bd^{*}$ is invariant under $f\mapsto f^{*}$ — it is
**not** required to vanish.

**(E12) Proposition 1 — pure-$Z$ centraliser.** **[proved]**

$$\langle (0\mid w),\,(h_x\mid h_z)\rangle_s = w\cdot h_x
\;\;\Longrightarrow\;\;
\{w : (0\mid w)\in S^{\perp_s}\} \;=\; \ker H_X ,$$

independent of $C,D$ — identical to the parent's.

**(E13) Proposition 2 — pure-$Z$ stabiliser subgroup.** **[proved]**  A general
row-space element is $(u,w)H=(uH_X,\;uP+wH_Z)$, whose $x$-part vanishes iff
$uH_X=0$.  Hence

$$S_Z \;=\; \operatorname{rowspace}(H_Z)\;+\;\{\,uP \;:\; u\in \mathrm{LN}(H_X)\,\},
\qquad \mathrm{LN}(H_X)=\{u: uH_X=0\}.$$

**(E14) Proposition 3 — dimension bookkeeping.** **[proved]**  With
$r_X=\operatorname{rank}H_X$, $r_Z=\operatorname{rank}H_Z$,

$$\delta \;:=\; \dim S_Z - r_Z \;\in\; [\,0,\;\ell m - r_X\,],$$
$$\operatorname{rank}H \;=\; r_X+r_Z+\delta,
\qquad
\boxed{\;k(Q) \;=\; k\big(\mathrm{BB}(A,B)\big)-\delta\;}$$

**The perturbation can only remove logical qubits.**

**(E15) Theorem 1 — distance ceiling at $\delta=0$.** **[proved]**  If $\delta=0$
then $S_Z=\operatorname{rowspace}(H_Z)$, the nontrivial pure-$Z$ logicals of $Q$
are *exactly* the parent's, and

$$d(Q)\;\le\;d_Z\big(\mathrm{BB}(A,B)\big).$$

**(E16) Corollary 1 — parent domination at $\delta=0$.** **[proved]**  Combining
(E14), (E15) and (E9), the parent $P=\mathrm{BB}(A,B)$ satisfies

$$n(P)=n(Q),\qquad k(P)=k(Q),\qquad d(P)=d_Z(P)\;\ge\;d(Q),$$

so $P$ **weakly dominates $Q$ in $[[n,k,d]]$** while having lower maximum check
weight and a purely CSS syndrome circuit.

**(E17) Theorem 2 — CSS shadow, arbitrary $\delta$.** **[proved]**  Define the CSS
code $Q'$ with $X$-checks $H_X$ and $Z$-checks $S_Z$.  Then $H_XS_Z^{\mathsf T}=0$ and

$$n(Q')=n(Q),\qquad k(Q')=k(Q),\qquad d(Q)\;\le\;d_Z(Q').$$

**Scope.**  This is a *pure-$Z$ ceiling*.  It does **not** bound $d_X(Q')$, so it
does **not** assert $d(Q)\le d(Q')$.  Counterexamples exist: `phase2_58`,
`phase2_60` ($[[72,4,6]]$, $\delta=4$) have $d_Z(Q')=6=d(Q)$ but $d_X(Q')=4$.

---

## Part IV — Circuit-level consequences

**(E18) Lemma C1 — measurement-correctness criterion.** **[proved]**  In a
one-ancilla-per-check circuit with $J(a,b)=\{j: P^a_j,P^b_j \text{ anticommute}\}$,
the circuit measures the intended generators **iff** for every pair $(a,b)$

$$\#\{\,j\in J(a,b)\;:\; a \text{ acts before } b\,\}\;\equiv\;0 \pmod 2 .$$

$|J(a,b)|$ is even because $g_a,g_b$ commute, so the condition is symmetric in $a,b$.

**(E19) Proposition C2 — one-ancilla depth lower bound.** **[proved]**

$$T\;\ge\;\max\Big(\max_a|\operatorname{supp}(g_a)|,\;\max_j \deg j\Big).$$

**(E20) Theorem C3 — depth 7 forced for weight-6 BB, TI schedule class.**
**[exact, translation-invariant class]**  For
$[[72,12,6]]$, $[[108,8,10]]$, $[[144,12,12]]$, $[[288,12,18]]$: CP-SAT returns
`INFEASIBLE` at $T=6$ under (E18) **within the translation-invariant (orbit)
schedule model** and a valid schedule at $T=7$ (achievement is class-free).
Within that class this *derives* the depth-7 cycle of arXiv:2308.07915 rather
than assuming it; the unrestricted depth-6 exclusion is external
(ASC, arXiv:2603.21499).

**(E21) Theorem C4 — depth separation (basis-independent via E21′).**
**[proved + exact]**  Every one of the 14 catalogue PBB $[[144,12,12]]$ codes has
$\max_a|\operatorname{supp}(g_a)|\ge 8$ **in its published generating set**, so by
(E19)

$$T_{\mathrm{PBB}}\;\ge\;8 \quad\text{vs}\quad T_{\mathrm{Gross}}=7 ,$$

for any one-ancilla schedule **measuring those generators**, translation-invariant
or not.

**(E21′) Basis-independence — PROVED (EXP-023).** **[proved + exact computation]**
For every one of the 14 members, with $V_w=\operatorname{span}\{v\in S:
\operatorname{wt}_s(v)\le w\}$:

$$\operatorname{rank}(V_6)\;=\;\operatorname{rank}(V_7)\;=\;66\;<\;132=n-k ,$$

by a two-sided sandwich: the theorem gate — *no element of $S$ with nonzero
$X$-part has weight $\le 7$* — is CP-SAT INFEASIBLE for all 14 (certified-
complete $X$-codeword enumeration + exact coset tests; independent second
encoding agrees; positive SAT controls pass), so $V_7\subseteq S_Z$ whose rank
is exactly 66; the published weight-6 $Z$-rows give the matching lower bound.
Hence **every generating set contains a generator of weight $\ge8$** and (E21)
holds for any one-ancilla schedule measuring **any** generating set.  Gross
control: $\operatorname{rank}(V_5)=0$ and the weight-6 rows span rank 132, so
its optimal max generator weight is exactly 6.  Machine verification lives in
EXP-023's artifacts and `tests/test_exp023_light_gensets.py` (not exp021).
The $\{uP\}$ caution below is resolved for these codes: the kernel terms do
not enlarge $S_Z$ beyond rank 66 ($\delta=0$, consistent with Theorem 1).

**(E22) Scope limit of (E19) — the incidence-counting route fails.** **[exact]**
(E19) is derived from the premise *"each check-qubit incidence $(g,j)$ needs its
own two-qubit gate touching $j$."*  That **premise** is false once
ancilla-ancilla gates are allowed.  The circuit

$$\texttt{R }a_1a_2;\;\; \mathrm{CX}\,q_0{\to}a_1;\;\; \mathrm{CX}\,a_1{\to}a_2;\;\;
\mathrm{CX}\,q_1{\to}a_1;\;\; \mathrm{CX}\,q_2{\to}a_2;\;\; \texttt{M }a_1a_2$$

measures $Z_0Z_1$ and $Z_0Z_2$ using **one** gate at $q_0$ although $\deg(q_0)=2$.

**What this is not.**  The circuit uses $T=3$ two-qubit layers while
$\max_j\deg(j)=2$, so it *satisfies* $T\ge\max_j\deg(j)$.  It is **not** a
counterexample to the inequality — only to the derivation.  Whether
$T\ge\max_j\deg(j)$ holds for arbitrary multi-ancilla schemes is **open**.

**A second gap in the same route.**  Incidence counts vary across the family:
4 of the 14 $[[144,12,12]]$ members (including `12_6_0193`) have
$\sum_g|\operatorname{supp}(g)|=1008$ over $144$ qubits — average degree
**exactly 7** — so counting yields only $T\ge7$ for them and cannot prove the
family-wide claim.  (E21)'s $T\ge8$ rests instead on
$\max_a|\operatorname{supp}(g_a)|\ge8$ (min 8, max 10 over all 14) via (E19),
valid for one-ancilla schedules measuring the published generators.

---

## Verification map

| eq | statement | how verified | instances |
|---|---|---|---|
| E1–E3 | conventions | definition check + negative control | random |
| E5 | (R1), (R2) | exact GF(2) | random polynomials |
| E7 | CSS commutation | exact GF(2) | 7 Bravyi + random |
| E8 | BB dimension | formula vs direct rank | 7 Bravyi + random |
| E9 | $d_X=d_Z$ | $\sigma$ row-space map + certified distances | 7 Bravyi |
| E11 | Lemma 0 both directions | valid **and** invalid $(C,D)$ | catalogue + random |
| E12 | pure-$Z$ centraliser | $\ker H_X$ equality, $C,D$-independence | catalogue |
| E13 | $S_Z$ decomposition | mutual span containment | catalogue |
| E14 | $k=k_{\rm BB}-\delta$, $\delta\ge0$ | formula vs direct rank | 368 catalogue |
| E15/E17 | distance ceilings | certified sector minima | small codes |
| E19 | depth bound | vs built circuits | survey codes |
| E21 | max check weight $\ge 8$ **for the published generators** | exact count | 14 PBB $[[144,12,12]]$ |
| E21′ | basis-independence | **OPEN** — procedure in exp022 | — |
| E22 | forwarding refutes the *premise* (not the inequality); counting gives only $\ge7$ | Stim simulation + arithmetic | explicit circuit |
