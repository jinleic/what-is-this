# The conditional record: $\forall_6$ for $\mathbb Z$ in $\mathbb Q$, modulo one hypothesis per cell

Draft, 2026-08-18. Status: **draft theorem document** — the verified content lives in
`h10q.py` (frozen tables, replayed every run), `THEOREMS.md` (exact statements L6–L13),
`RESULTS.md` (scope table), and `data/`. Every number below is sourced to a file and table.
Naming, formulas, and conventions follow `THEOREMS.md` exactly ($K=\mathbb Q$, $S=\{2\}$,
$\pi=2$, $u=1$; $A=1+4a^2$, $\delta_\tau=1-A\tau^2$, $\alpha=-\delta_\tau A$,
$Z=z^3$, $D_z=1-Z-a^2Z^2$).

---

## 1. The statement

**Theorem C (conditional).** Assume **H** (§2). Then:

1. $\mathbb Q\setminus\mathbb Z$ is Diophantine in $\mathbb Q$ with **6 unknowns**, realized
   by the L6 formula
   $$F(z):=\exists b,s\;\bigl[(1+2s,b)\in\Phi_1^{\{2\}}\ \wedge\
   \Theta_*\bigl(1+2s,\,b,\,h(1+2s,b,z^3),\,s\bigr)\bigr],$$
   where $\Theta_*(a,b,c,s)$ is the two-branch tied block
   $\exists y\,r\ \bigl[c^2-Ay^2-16Br^2=16-16ABs^2\bigr]\vee
   \bigl[c^2-Ay^2-16ABr^2=16A-16A^2Bs^2\bigr]$, $B=2b$.
2. $\mathbb Z$ is definable in $\mathbb Q$ by a formula with **6 universal quantifiers**
   ($\mathbb Z = \mathbb Q\setminus(F\vee\mathfrak m_2)$).
3. $\operatorname{efd}_{\mathbb Q}(\mathbb Q\setminus\mathbb Z)\le 5$.
4. (Transferred) The $\forall_9\exists_6$-theory of $\mathbb Q$ is undecidable — the same
   code-transfer as the recorded $\forall_9\exists_7$ row (`THEOREMS.md`, consequences
   table; Daans Cor 6.2 + Sun 2021). Daans' Question 5.7 would narrow to $2\le m\le 6$;
   the lower bound $m\ge2$ (DDF §8, Daans Thm 2.6) is untouched.

**Count (audited in `THEOREMS.md` L6).** Base unknowns $(b,s)$: 2. The $\Phi$-membership of
$(1+2s,b)$ is $\exists_3$ (Daans Lemma 5.2 pulled back along the polynomial map
$(b,s)\mapsto(1+2s,b)$, which adds no witness). The tied block is $\exists_2$
($y,r$; $s$ is a base coordinate). Daans–Dittmann–Fehm, *Existential rank and essential
dimension of diophantine sets*, arXiv **2102.06941v5**, **Theorem 1.4** — over a field
finitely generated over a perfect subfield, the intersection of $\exists_{m_1}$ and
$\exists_{m_2}$ subsets of the same space is $\exists_{m_1+m_2-1}$ — gives $3+2-1=4$
witnesses for the conjunction on the common base, hence $2+4=6$ after projecting away
$(b,s)$. If either set had existential rank 0 (the only unmet hypothesis of Thm 1.4), the
trivial conjunction bound is $\le 3$, so the count is never worse (frozen audit note,
L6 "Count [AUDITED]"). The branch disjunction multiplies two polynomials in the same
$(y,r)$ and adds nothing (L11a(3)); denominator clearing and the $\delta_\tau\ne0$
condition add nothing because $A\equiv5\bmod8$ is not a square in $\mathbb Q_2$ (L6,
"No hidden denominator or inequation variable"). Adjoining the one omitted place
$\mathfrak m_2$ as a finite $\exists_3$ disjunct keeps the count at 6 (L6, Status).
In efd form the count is $\operatorname{rk}^\exists=\operatorname{efd}+1\le6$
(DDF 4.16 + 5.10 + Cor 4.12, recorded as `THEOREMS.md` A3), i.e. the tie rigidifies one
$(a,b)$-direction and drops the generic fibre bound from 6 to 5.

*Citation note.* The drafting instruction for this document carried the arXiv id
**1612.03992** for Theorem 1.4; that id resolves to an unrelated astrophysics preprint
(checked 2026-08-18). The citation recorded in `NOTES.md` (References) and used here is
**arXiv 2102.06941v5**, Thm 1.4 and Cor 5.11.

**Soundness (no hypothesis consumed).** Any solution of $F(z)$ is a $\Psi_0$- or
$\Psi_{2a/A}$-solution of Sun's block for $(a,b)=(1+2s,b)$, and soundness of the block is
**$\tau$-uniform**: the Sun identities (2.2)–(2.4) hold for arbitrary $\tau$
($\operatorname{Nrd}(\gamma_\tau)=1$ identically) and Prop 2.1 ($c$ integral at every
finite ramified place of $(A,B)$) never mentions $\tau$ (`THEOREMS.md` L11a(2); machine
checks `(A)`, `(B)` in `_verify_sun`). The branch menu used by the class side is L11b's
table of the four square classes of $\langle-1,A\rangle$
($\tau\in\{0,\,2a/A,\,1,\,\tau^\dagger=(1+2a^2)/A\}$); the square branch $\tau^\dagger$
with $\alpha=(2a^2)^2$ is the one carrying all 293 verified classes below. Hence
$F(z)\Rightarrow z\in\bigcup_{w\ \mathrm{odd}}\mathfrak m_w$, unconditionally.

**Completeness = H.** The converse — the L6 assembly lemma — is precisely H of §2: for
each target cell, a verified aligned class provides a candidate at every place, and an
emergent-free member of it makes the attached tied conic globally soluble (W2,
Hasse–Minkowski), producing a witness of $F(z)$. This document changes nothing about
soundness; it packages the remaining gap as a clean analytic statement.

---

## 2. The hypothesis H

**Setting of one cell.** A *cell* is a pair $(w,z)$ with $w$ an odd prime and
$z\in\mathbb Q$, $v_w(z)\ge1$. On the square branch $\tau=\tau^\dagger=(1+2a^2)/A$
($\delta_\tau=-4a^4/A$, $\alpha=(2a^2)^2$ a perfect square), a *verified aligned class*
for the cell is data $(a,\varepsilon,f,q,N)$ with:
- $a$ admissible ($a=1+2s$ a 2-adic odd unit; $A\equiv5\bmod8$);
- $f=1$ (*L11 family*: available iff $z$ has a numerator prime $p\equiv1\bmod4$;
  `THEOREMS.md` L11h) or $f\mid\operatorname{num}(z)$ prime (*L12 escape family*;
  `THEOREMS.md` L12b/L13a);
- $b=\varepsilon f Q$, $Q$ prime, $Q\equiv q\pmod N$, with controlled set
  $S=\{2,3,5,7\}\cup\operatorname{supp}(\alpha)\cup\operatorname{supp}(\delta_\tau)
  \cup\{f\}$
  (at $a=1$: $S=\{2,3,5,7,f\}$) and $N$ built from the exponent lemma
  (`h10q.py::_l10_exponent`) so that every frozen symbol $(x_0,d_0)_p$, $p\in S$, the wild
  symbol $(x_0,d_0)_Q=(A\mid Q)$, and the real symbol are all $+1$ on the class;
- a finite *excluded set* of members ($Q\mid D_z\cdot\operatorname{num}(z)\cdot a\cdot
  \delta_\tau$) at which the class claim is void.

Here $x_0=\alpha M_0$, $d_0=2\alpha b$, $M_0=16-\delta_\tau c_0^2-32Ab\,s_{\rm tie}^2$,
$c_0=h(a,b,Z)$, and the value polynomial is
$$P(b)=16D_z^2A^2b^4-\delta_\tau a^4Z^4N_g(b)^2-32A^3s^2D_z^2b^5,\qquad
N_g=16a^4b^2-A(b-1)^4,$$
with the kernel identity $x_0=\alpha\,P(b)/(D_z^2A^2b^4)$ asserted at every replay point
(`l12_class.py`, module docstring and assertion; the kernel's `_l10_P` differs by the
inert square factor $A^2$ — `THEOREMS.md` L13, Normalization note).

**Hypothesis H (one sentence, per cell).**

> **H.** For every odd prime $w$ and every rational $z$ with $v_w(z)\ge1$ there exists a
> verified aligned class for the cell $(w,z)$ — prime-$b$ data $(a,\varepsilon,1,q,N)$
> (L11 family) or escape data $(a,\varepsilon,f,q,N)$ with $f\mid\operatorname{num}(z)$
> (L12 family) — and a prime member $Q\equiv q\pmod N$ of that class, outside its finite
> excluded set, at which the *bad* emergent set is empty: every odd prime
> $l\notin S\cup\{Q\}$ with $v_l\bigl(P(\varepsilon fQ)\bigr)$ **odd** has
> Hilbert symbol $(x_0,d_0)_l=+1$.

> **Terminology (corrected 2026-08-19).** An earlier wording of H required every
> outside valuation to be *even*. That is strictly stronger than the kernel/L13/L14
> contract and is contradicted by the 197 accepted prime-rung closures, in which a
> single outside prime has odd valuation and reciprocity-forced symbol $+1$
> (`l13_filter.py` `smooth_emergent`/`cofactor_decide`; `THEOREMS.md` L13f, L14, L17).
> The operational contract above — *no place with symbol $-1$* — is the one the
> machine verifies, the one the 293/293 instance verification establishes, and the
> one implied by Schinzel prime values (L19).

Two uniform equivalences make the final clause exactly the analytic content:

- **The $v_l$ criterion.** At any odd $l\notin S\cup\operatorname{supp}(b)$,
  $v_l(x_0)$ is odd **iff** $v_l(P(b))$ is odd — the denominators $D_z^2A^2b^4$ are
  squares and $v_l(\alpha)=0$ since $\operatorname{supp}(\alpha)\subseteq S$
  (`THEOREMS.md` L10 identity and L11f; `THEOREMS.md` L13a defines emergent places by
  exactly this dichotomy).
- **The parity law.** The set of places where the quaternion algebra $(x_0,d_0)$
  ramifies has even cardinality (`h10q.py::ramified`; `THEOREMS.md` L9 product-over-$T$,
  L7c). With frozen, wild, and real symbols all $+1$ on the class, the number of emergent
  $-1$s is therefore **automatically even**, so "zero emergent places" is equivalent to
  "all Hilbert symbols $+1$", and the associated tied conic is globally soluble
  (W2 + Hasse–Minkowski). Measured exactness: the bad count is even at every fully
  factored member in both recorded sweeps — 39/39 (`data/l12b_class_sample.json`) and
  44/44, histogram $\{0{:}9,\ 2{:}27,\ 4{:}6,\ 6{:}2\}$ (`data/l12_density.json`).

Per family:

- **F1 — L11 family** ($z$ has a numerator prime $\equiv1\bmod4$): the class itself is
  proved to exist (`THEOREMS.md` L11h, constructive sufficiency); H asks only for one
  emergent-free prime member.
- **F2 — L12 escape family** (all numerator primes $\equiv3\bmod4$ or $2$): coprimality
  of $b$ to the cell data is provably fatal (`THEOREMS.md` L12a, generalized wall), so the
  existence of the escape class is part of H's per-cell assertion; H then asks for one
  emergent-free prime member.
- **F3 — directly witnessed cells**: cells with an independently certified soluble witness
  carry no hypothesis at all.

H's finitary content (class data and certificates) is machine-checkable per cell; its only
non-finitary content is the emergent-free member.

---

## 3. Evidence on the 353-cell grid

Grid: the 353 cells $(w,(u_1,u_2))$ with $w$ an odd prime $<100$ and $u$ in the 15-value
pool `_L9_U_POOL` (`h10q.py::_l9_grid`). Sources: kernel tables, `data/` artifacts,
`THEOREMS.md` scope statements, `RESULTS.md` scope-table rows.

| Object | Count | Source (file + table) | Status |
|---|---|---|---|
| Aligned classes, L11 family | **190/353** | `h10q.py::_L11_CLASSES`; `data/l11_classes.jsonl` (190 rows); replayed in `_verify_L11` | **PROVED** (L11h constructive sufficiency) + **VERIFIED-instance** (every frozen row re-verified per run); necessity scope: now **proved for all admissible $a$ via L13e** (the wall extends to composite squarefree $A$ — `data/l13_compositewall.{md,json}`) |
| Depth audit of L11 classes | 162 members / 30 cells | `data/audit_l11h_members.jsonl` | **VERIFIED-instance** (deciders agree everywhere; 5 zero-bad members found there, parity law exact; refusals logged, never evidence — `NOTES.md`, L13 log item 2) |
| Walled cells | **163** | complement of `_L11_CLASSES` in `_l9_grid`; wall `THEOREMS.md` L12a + L13e | **PROVED** (all admissible $a$: prime, composite squarefree kernel, or rational) |
| — of which escape classes | **103** | `h10q.py::_L12_ESCAPE`; cert `l12_class.py::l12_class_cert` | **PROVED** (frozen side, L13a): minimal controlled set $S_{\min}=\{2,3,5,7,f\}$, $\gcd(q,N)=1$ |
| — sampled members of the escapes | **883 primes, 0 violations** | `data/l12b_class_sample.json` (103 rows) | **VERIFIED-instance**: median modulus $N=24\,202\,080$, 0 failures, 0 Cauchy-sign failures, 0 excluded-member hits; 174 distinct excluded primes named per row; 12 rows $\times$ 2 members replayed in-suite each run (40 rows extended) |
| — of which direct witnesses only | **60** | `h10q.py::_L13_ESCAPE2`; `data/l12_escape2.jsonl` (60 rows: branch sq; $a\in\{3,5,33\}$, 40/10/10; `fully_soluble`/`verified` true on all 60) | **PROVED per row**: 55 Hasse–Minkowski + 5 steered certificates (`THEOREMS.md` L13b) |
| H instance-verified on the grid | **293/293 classes** | `h10q.py::_L13_H_CLASSES` + `_verify_L14` (default 60-seed / extended all); `data/l13h_all_closures.json`, `data/l13h_replay.jsonl`; engine upgrade: BLS-relaxed Pocklington (`test_bls.py`) | **PROVED per row**: every aligned class carries a certified emergent-free member (197 prime-rung + 96 factorint; lead replay 293/293, ramified empty 271, 22 refusals logged). **H on the grid: CLOSED**; the uniform all-$w$ form remains the open analytic input |
| H remainder + counting laws on the grid | 293 exact row audits | `data/l15_remainders.jsonl` (0 refusals, 0 alarms: $R$ squarefree in all 293; $R{=}1$ never; 197 $R{=}p$ exactly, 96 factorint $e\in\{2,3,4\}$; parity law exact on all 15 odd-$|E|$ rows); `data/l15_density.json` ($k_{\text{zero}}{=}0$ in 70.65%; tail median 28; $p_0\approx0.040$; resistant subclass $a{=}1$, $w\equiv3(4)$) | **PROVED per row** (exact recomputation under discipline); the statistical claims are **EVIDENCE** for the uniform all-$w$ shape of H |
| Residual-cell witnesses | **7** | `h10q.py::_L13_RESIDUAL` (previous 6 + the $(89,(-2,1))$ closure: $a{=}3$, $b=89/367$, $\tau=19/37$, frozen 2026-08-18) | **PROVED per row** (tied-status True, replayed every run) |
| Directly witnessed cells, total | **65** | 60 + 7 above minus the two prior algebraic closures already in the 347 (one escape-class cell, three L11 cells, 60 witness-only walled cells, + cell $(89,-2)$) | **PROVED per row** |
| Global steered coverage | **347/353** | `h10q.py::_L9_STEERED` (345) + recorded wider-$s$-pool hits $(29,-2)$, $(53,(1,3))$ (`THEOREMS.md` L11, Consequence) | **VERIFIED-instance** (345-row table replayed byte-identically every run; the 2 recorded hits lead-replayed) |
| Total witness coverage | **353/353 — COMPLETE** | union of the above: $347 + (31,-2) + 4$ residual cells + $(89,-2)$ via L13f | **PROVED per row** (every cell carries its own replayable certificate; the grid is closed since 2026-08-18) |
| Certified emergent-free members | **24** earlier; **293 class-complete L14** | `data/l12_density.json` (9 of 44 factored members, sweep 1) + `data/l13_zerobad_more.jsonl` (15 of 57, sweep 2); even-only histograms — superseded on the grid by the L14 row below | **EVIDENCE** (positive empirical density of step (ii)); 15 distinct cells, 16 frozen rows `h10q.py::_L13_ZERO_BAD` replayed each run — **PROVED per row** |
| Emergent $-1$s are real | 467/959 spot-probes | `data/l12b_class_sample.json` | **VERIFIED-instance** (step (ii) is genuinely outside the class claim) |
| No structural shortcut | 706/706 rows irreducible deg 8 | `THEOREMS.md` L13c; `data/l12_squarehunt.json`; `data/l12_factorP.json` | **PROVED per row** (sympy `factor_list`; in-suite Frobenius re-certification) + **EVIDENCE** (137.5M square-class tests, 0 hits; Galois $D_8\wr C_2$ consistent, 28,240 certified samples) |
| Cell $(89,(-2,1))$ — **CLOSED 2026-08-18** | 325 certified soluble candidates | `data/l13_filter_run2.json` (236 in the structured box: 133 proved-prime-rung + 103 factorint-rung) + dense-box block (89: 46 + 43); witness $(a,b)=(3,\,89/367)$, $\tau=19/37$, frozen `_L13_RESIDUAL` row 7 | **PROVED** (L13f cofactor ladder: zero-bad ⟺ stripped remainder is a square OR a single proved prime — parity law forces its symbol $+1$; lead-replayed all five rungs; INCONSISTENT = 0 across 5,767 decided rows; 1,896 + 1,637 + 937 refusals recorded with reasons, never evidence) |

The closure does **not** touch hypothesis H below: it certifies ONE witness for the
grid's last cell, whereas H requires an emergent-free member of *each aligned class*
uniformly (`THEOREMS.md` L13f).

---

## 4. The blocking gap

H restricts the **square class of values of an irreducible degree-8 polynomial at primes
in an arithmetic progression**: for each cell, $P\in\mathbb Q[t]$ is irreducible of
degree 8 over $\mathbb Q$ (`THEOREMS.md` L13c: all 706 (cell, branch) rows certified; the
content's squarefree part is confined to $\{\pm1,\pm5\}$, so the constant factor never
introduces an emergent place), and H asks for one prime $Q\equiv q\pmod N$ (outside a
finite excluded set) such that every odd multiplicity $v_l(P(\varepsilon fQ))$ with
$l\notin S\cup\{Q\}$ is even — equivalently, the squarefree kernel of $P(\varepsilon fQ)$
is supported on the fixed finite set $S\cup\{Q\}$. One member per class suffices; no
density, no infinite subsequence, no squarefreeness of the whole value is required.

**Upgrade (L19, 2026-08-19): the analytic input is now a NAMED classical conjecture.**
For a fixed verified aligned class, write $F(t)=P(\varepsilon f(q_1+Nt))=c\,G(t)$ with
$c>0$ supported on $S$ and $G\in\mathbb Z[t]$ primitive with positive leading
coefficient. Then **Schinzel's Hypothesis H for the pair $\{q_1+Nt,\;G(t)\}$ implies the
per-class clause of H — indeed infinitely many qualifying members** (PROVED implication;
the conclusion remains CONDITIONAL on Schinzel). Mechanism: at a simultaneous prime
value, stripping $S$ leaves the single odd-valuation place $R=G(t)$; the moving prime
satisfies $v_Q(P(b))=0$; every other outside place has even valuation and symbol $+1$;
the class certificate pins the frozen, moving and infinite symbols to $+1$; and Hilbert
reciprocity then forces $(x_0,d_0)_R=+1$, i.e. ladder verdict *zero*. The Schinzel
hypotheses themselves are verified exactly on all 293 canonical pairs (primitive
irreducible degree-8 $G$ with replayed Frobenius certificates, positive leading
coefficients, pair-product fixed divisor 1 by an exact degree-9 finite-difference
certificate, $G$ an $S$-unit for every integer $t$, frozen signs uniform in $t$):
`data/l18_schinzel_implies_h.jsonl`, generator `l18_schinzel_implies_h.py`,
cross-checked on 24 prime-rung rows. **Not supplied by Schinzel:** existence of an
aligned class for every cell (that is the separate finite/structural side, characterized
by L11h and verified 293/293 on the grid).

Honesty points:

- The polynomial cannot be split or collapsed: $P$ is irreducible degree 8 on every grid
  row, and no $P=cR^2$ or $P=cR^2S$ shortcut exists anywhere on the grid
  (`THEOREMS.md` L13c), so the gap is exactly as phrased and no finer.
- The condition at primes is strictly stronger than the classical "squarefree values at
  integers" questions: here the argument itself is prime and must sit in a fixed residue
  class modulo a scale-dependent modulus $N$ (median $2.4\times10^7$ over the 103 sampled
  escape classes; `data/l12b_class_sample.json`).

**Literature status (LitScout verified report, 2026-08-18; every item in this list was
checked at the source by the recon agent).** The required input sits outside the
unconditionally solved range at every angle:

- *Integer argument, unconditional:* squarefree (and square-class) values of
  single-variable polynomials are understood only in degree $\le3$ — Erdős,
  J. London Math. Soc. 28 (1953) 416–425; Hooley, "On the power-free values of
  polynomials", Mathematika 14 (1967) 21–26 (asymptotic); for binary forms, Greaves,
  "Power-free values of binary forms", Q. J. Math. Oxford (2) 43 (1992) 45–65. The
  univariate **quartic** case is already open; our degree is 8.
- *Conditional ceiling:* Granville, "ABC allows us to count squarefrees", IMRN 1998 (19)
  991–1009 (under $\mathsf{abc}$, any degree, integers); Poonen, Duke Math. J. 118 (2003)
  353–373 (multivariable analogue).
- *Values at primes:* degree 1 — Baker–Pollack, "Clusters of primes with square-free
  translates" (Rev. Mat. Iberoam., submitted; pollack-math.net/ClustersPrimesSquarefree.pdf);
  degree 3 — Helfgott, "Square-free values of $f(p)$, $f$ cubic", arXiv:1112.3820
  (squarefree values at primes), and Reuss, "Power-Free Values of Polynomials",
  arXiv:1307.2802 (asymptotics for $(d-1)$-free values $f(p)$, degree $d\ge3$).
  **Nothing for degree 8; nothing prescribing a square class.**
- *Square class at a prime argument in an arithmetic progression* (H's shape):
  degree $\le2$ only — Krumm, "Squarefree parts of polynomial values", JTNB 28 (2016)
  699–724, arXiv:1407.4890, Thm 1.3 with Props 3.4–3.5 (unconditional); degree 3 is
  conditional on the **Parity Conjecture for elliptic curves over $\mathbb Q$** (Krumm
  Prop 3.8 — rank parity, not the sieve "parity problem").
- *Conic-bundle route* (the L8 "Direction" attack on the same assembly): HSW,
  Compositio Math. 150 (2014) 2095–2111 (arXiv:1304.3333), Thm 3.1/Cor 3.4 — unconditional
  only when degenerate fibres have abelian constant field over $\mathbb Q$; the general
  case is conditional on Schinzel's Hypothesis H (Colliot-Thélène–Skorobogatov–
  Swinnerton-Dyer, J. reine angew. Math. 495 (1998) 1–28; the CT–Sansuc Schinzel paper
  is Acta Arith. 41 (1982) 33–53). The value divisor here is irreducible of degree 8 with non-abelian
  Galois group (best evidence: $D_8\wr C_2$; `THEOREMS.md` L13c) — outside the unconditional
  case.
- *Negative reservation:* do not cite Carella arXiv:2310.16952 (math.GM, unrefereed,
  claims contradicting the quartic open-problem status); it is not part of the record.

Within the verified literature, the closest published anchors to H's analytic input are
the degree-3 results of Helfgott and Reuss at primes; the degree-8 square-class input H
needs is unconditional neither at integers nor at primes, and the theorem of §1 stays
conditional as stated.

---

## 5. Provenance and discipline

- Authorities: `h10q.py` frozen tables (`_L9_STEERED`, `_L11_CLASSES`, `_L12_ESCAPE`,
  `_L13_ESCAPE2`, `_L13_RESIDUAL`, `_L13_ZERO_BAD`), replayed on every suite run;
  `THEOREMS.md` L6–L13; `RESULTS.md` scope table; `data/` artifacts cited per row above.
- Engine discipline: proven primality only — deterministic MR below the A014233
  13-base bound; above it, Pocklington (trial division to $10^6$, early abort at
  $F\ge n^{1/3}$) with the **Brillhart–Lehmer–Selfridge relaxation** (two-factor
  discriminant test, CP Thm 4.1.5), upgraded 2026-08-18 (`test_bls.py` unit check);
  `FactorBudget`/`PrimalityBound` refusals are logged and are never evidence.
- Nothing here rests on Sun's unrefereed §§3–8 chain: soundness is inherited from the
  audited bridge/block chain, and completeness is isolated as H (`THEOREMS.md` L6, Status:
  "independent of whether Sun's own §§3–8 chain survives refereeing").
- Assembly remains **OPEN** in the unconditional sense; the $\forall_6$ count and
  $\operatorname{efd}\le5$ are conditional on H, and H is certified at 5 fixed members
  (`_L13_ZERO_BAD`) with 9 measured instances and no structural shortcut.
