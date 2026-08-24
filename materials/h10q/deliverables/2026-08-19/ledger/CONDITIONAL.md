# The conditional record under classical Schinzel H: $\forall_6$ for $\mathbb Z$ in $\mathbb Q$

Status: **CONDITIONAL on classical Schinzel's Hypothesis H alone**.
Classical Schinzel H is unproved.  The implication is proved in
`THEOREMS.md` L19–L22; the six-unknown architecture and count are
`THEOREMS.md` L6.  The intermediate per-cell statement H is retained in
§2, but it is now a proved consequence of classical Schinzel H, not an
additional conjectural premise.

Conventions follow `THEOREMS.md`: $K=\mathbb Q$, $S=\{2\}$,
$\pi=2$, $u=1$,
$$
A=1+4a^2,\quad \delta_\tau=1-A\tau^2,\quad
\alpha=-\delta_\tau A,\quad Z=z^3,\quad D_z=1-Z-a^2Z^2.
$$
(Source: `THEOREMS.md`, L6 and L22.)

---

## 1. The statement

**Theorem C (CONDITIONAL).** Assume classical Schinzel H.  Then:

1. $\mathbb Q\setminus\mathbb Z$ is Diophantine over $\mathbb Q$ with
   **six unknowns**, realized by
   $$
   F(z):=\exists b,s\;\bigl[(1+2s,b)\in\Phi_1^{\{2\}}\ \wedge\
   \Theta_*\bigl(1+2s,b,h(1+2s,b,z^3),s\bigr)\bigr].
   $$
   (Source for the formula and count: `THEOREMS.md`, L6;
   chain closure: `THEOREMS.md`, L22d.)
2. $\mathbb Z$ has a definition over $\mathbb Q$ with **six universal
   quantifiers**, since
   $\mathbb Z=\mathbb Q\setminus(F\vee\mathfrak m_2)$
   (`THEOREMS.md`, L6 Status and L22d).
3. $\operatorname{efd}_{\mathbb Q}
   (\mathbb Q\setminus\mathbb Z)\le5$
   (`THEOREMS.md`, A3 and L6 Count).
4. The same code transfer gives undecidability of the
   $\forall_9\exists_6$-theory of $\mathbb Q$, and Daans' range becomes
   $2\le m\le6$ under the same assumption
   (`THEOREMS.md`, consequences table; Daans Corollary 6.2 and Sun
   2021 are recorded there).

The finite branch menu is unchanged:
$$
\tau\in
\left\{0,\ \frac{2a}{A},\
\tau^\dagger=\frac{1+2a^2}{A}\right\}.
$$
For each branch, $\Theta_*$ uses the same two tied variables $(y,r)$ and
the cleared equation
$$
\exists y,r\quad
c^2-Ay^2-16\delta_\tau^{-1}Br^2
=16\delta_\tau^{-1}
\left(1-\frac{\delta_\tau}{16}\Delta_\tau(a,b,s)\right),
\qquad B=2b.
$$
The fixed branch used in the closure is
$\tau^\dagger=(1+2a^2)/A$ with
$\delta_{\tau^\dagger}=-4a^4/A$.  A finite disjunction multiplies the
cleared branch equations in the **same** $(y,r)$, so it costs no new
witness (`THEOREMS.md`, L11a and L22).

**Exact count.** The base variables are $(b,s)$.  The pulled-back
$\Phi$-membership has existential rank at most $3$, and the tied block
has existential rank at most $2$.  Daans–Dittmann–Fehm,
arXiv:2102.06941v5, Theorem 1.4, gives
$3+2-1=4$ witnesses for their conjunction on the common base, hence
$2+4=6$ after projection.  If either rank is zero, the elementary bound
is no larger.  Adding the finite $\mathfrak m_2$ disjunct, clearing
denominators, and adding the fixed finite branch menu do not increase
the count (`THEOREMS.md`, L6 “Count [AUDITED]” and “No hidden
denominator”; L11a).  Equivalently,
$\operatorname{rk}^{\exists}=\operatorname{efd}+1\le6$
(`THEOREMS.md`, A3).

*Citation correction.* The earlier drafting id 1612.03992 is an
unrelated astrophysics preprint; the rank theorem used here is
arXiv:2102.06941v5, Theorem 1.4 and Corollary 5.11
(`NOTES.md`, References; citation audit recorded in `NOTES.md`, L19
session).

**Soundness is unconditional.** Every branch in the displayed menu is a
specialization of Sun's $\tau$-uniform block: identities (2.2)–(2.4)
hold for arbitrary $\tau$, and Proposition 2.1 does not mention
$\tau$.  L11a proves the cleared finite disjunction preserves that
soundness.  Hence
$$
F(z)\Longrightarrow
z\in\bigcup_{w\ {\rm odd}}\mathfrak m_w
$$
without Schinzel H (`THEOREMS.md`, L6 and L11a; replay authority
`h10q.py::_verify_sun`).

**Completeness under classical Schinzel H.** Fix a target cell $(w,z)$.
L22 proves vertical irreducibility on the fixed
$\tau^\dagger$, $\delta=-4a^4/A$ branch for the fixed nonzero value
$Z=z^3$; quantitative Hilbert irreducibility in L20's odd
character-admissible $a$-progression selects a concrete irreducible
specialization (`/tmp/l22_elimination.md`; `THEOREMS.md`, L20–L22).
Together, L20 supplies the aligned $f=w$ class and conditions (b)–(d),
while L22 completes condition (a).  L19 applies classical Schinzel H to the
resulting linear/octic pair and supplies an emergent-free member.  That
is precisely intermediate H of §2; W2 and Hasse–Minkowski make the tied
conic globally soluble, producing a witness of $F(z)$
(`THEOREMS.md`, L6, L19–L22).  Thus
$$
\boxed{\text{classical Schinzel H}\Longrightarrow
       \text{intermediate H}\Longrightarrow\text{Theorem C}.}
$$

---

## 2. Intermediate per-cell H, and its derivation from Schinzel H

**One cell.** A cell is $(w,z)$ with $w$ an odd prime,
$z\in\mathbb Q$, and $v_w(z)\ge1$.  Work only on
$$
\tau=\tau^\dagger=\frac{1+2a^2}{A},\qquad
\delta_{\tau^\dagger}=-\frac{4a^4}{A},\qquad
\alpha=(2a^2)^2.
$$
(Source: `THEOREMS.md`, L20 and L22.)  A verified aligned class consists
of $(a,\varepsilon,f,q_1,N)$ with admissible odd
$a=1+2s$, $b=\varepsilon f(q_1+Nt)$, and controlled set
$$
S=\{2,3,5,7\}\cup\operatorname{supp}(\alpha)
\cup\operatorname{supp}(\delta_{\tau^\dagger})\cup\{f\},
$$
such that all frozen, moving and real Hilbert symbols are $+1$ and a
finite excluded set is named.  L20's uniform construction takes
$f=w$ (`THEOREMS.md`, L20; replay
`data/l19_classexist.jsonl`).

The value polynomial on this fixed branch is
$$
P(b)=16D_z^2A^2b^4-\delta_{\tau^\dagger}a^4Z^4N_g(b)^2
     -32A^3s^2D_z^2b^5,\qquad
N_g=16a^4b^2-A(b-1)^4.
$$
The kernel identity is
$x_0=\alpha P(b)/(D_z^2A^2b^4)$ and
$d_0=2\alpha b$ (`THEOREMS.md`, L13 normalization note and L22;
`l12_class.py`, asserted kernel identity).

**Intermediate H (per cell).**

> **H.** For every cell $(w,z)$ there is one verified aligned class on
> $\tau^\dagger=(1+2a^2)/A$,
> $\delta_{\tau^\dagger}=-4a^4/A$, and a prime member
> $Q=q_1+Nt$ outside its finite excluded set such that every odd prime
> $l\notin S\cup\{Q\}$ with
> $v_l(P(\varepsilon fQ))$ odd has
> $(x_0,d_0)_l=+1$.
> (Source: `THEOREMS.md`, L19–L22.)

“Emergent-free” means **no place with Hilbert symbol $-1$**; it does not
mean that every outside valuation is even.  The stronger historical
wording is contradicted by the $197$ prime-rung closures among the
$293$ recorded classes, each with one reciprocity-forced outside
symbol $+1$ (`data/l13h_all_closures.json`;
`data/l15_remainders.jsonl`; `THEOREMS.md`, L13f–L15).

The final clause is exact.  For
$l\notin S\cup\operatorname{supp}(b)$,
$v_l(x_0)$ is odd exactly when $v_l(P(b))$ is odd, because the displayed
denominator is a square and $\alpha$ is an $l$-unit
(`THEOREMS.md`, L10–L13).  Hilbert reciprocity makes the number of
remaining $-1$ symbols even; when none remain, W2 and
Hasse–Minkowski give global solubility (`THEOREMS.md`, L7c, L9 and W2).

**PROVED implication: classical Schinzel H $\Rightarrow$ H.**

1. L20 supplies a nonempty odd character-admissible progression of
   $a$-values for the cell, with $f=w$ (`THEOREMS.md`, L20;
   `data/l19_classexist.jsonl`).
2. L22 proves $P(a,z^3,b)$ irreducible in $\mathbb Q(a)[b]$ for the
   fixed nonzero $z^3$.  Quantitative Hilbert irreducibility selects an
   $a$ in that progression whose specialized octic is irreducible
   (`/tmp/l22_elimination.md`; `THEOREMS.md`, L21d and L22).
3. For this $a$, the normalized
   $F(t)=P(\varepsilon w(q_1+Nt))=cG(t)$ has positive lead; the square
   class of $c$ is $S$-supported; $\gcd(q_1,N)=1$; the pair has no fixed
   prime divisor; and the $S$-adic unit and frozen-sign conditions hold.
   Thus Schinzel conditions (a), (b), (c), and (d) are all **PROVED**
   (`/tmp/l20_admissible.md`; closure of (a):
   `/tmp/l22_elimination.md`; `THEOREMS.md`, L20 and L22).
4. Classical Schinzel H for
   $\{q_1+Nt,G(t)\}$ gives infinitely many simultaneous positive prime
   values.  After stripping $S$, there is one outside odd-valuation
   place $R=G(t)$; $v_Q(P(b))=0$; every other outside symbol is $+1$;
   reciprocity forces $(x_0,d_0)_R=+1$
   (`THEOREMS.md`, L19; replay
   `data/l18_schinzel_implies_h.jsonl`).

Consequently H is **CONDITIONAL**, because classical Schinzel H is
unproved, but it is no longer an independent hypothesis.  No finite
evidence is used in this implication.

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
| H instance-verified on the grid | **293/293 classes** | `h10q.py::_L13_H_CLASSES` + `_verify_L14`; `data/l13h_all_closures.json`, `data/l13h_replay.jsonl` | **PROVED per row**: every recorded aligned class has a certified emergent-free member.  This is finite corroboration; the all-cell implication now comes from classical Schinzel H through L19–L22, not from extrapolating the grid (`THEOREMS.md`, L22d). |
| H remainder + counting laws on the grid | 293 exact row audits | `data/l15_remainders.jsonl` (0 refusals, 0 alarms: $R$ squarefree in all 293; $R{=}1$ never; 197 $R{=}p$ exactly, 96 factorint $e\in\{2,3,4\}$; parity law exact on all 15 odd-$|E|$ rows); `data/l15_density.json` ($k_{\text{zero}}{=}0$ in 70.65%; tail median 28; $p_0\approx0.040$; resistant subclass $a{=}1$, $w\equiv3(4)$) | **PROVED per row** (exact recomputation under discipline); the statistical claims are **EVIDENCE** for the uniform all-$w$ shape of H |
| Residual-cell witnesses | **7** | `h10q.py::_L13_RESIDUAL` (previous 6 + the $(89,(-2,1))$ closure: $a{=}3$, $b=89/367$, $\tau=19/37$, frozen 2026-08-18) | **PROVED per row** (tied-status True, replayed every run) |
| Directly witnessed cells, total | **65** | 60 + 7 above minus the two prior algebraic closures already in the 347 (one escape-class cell, three L11 cells, 60 witness-only walled cells, + cell $(89,-2)$) | **PROVED per row** |
| Global steered coverage | **347/353** | `h10q.py::_L9_STEERED` (345) + recorded wider-$s$-pool hits $(29,-2)$, $(53,(1,3))$ (`THEOREMS.md` L11, Consequence) | **VERIFIED-instance** (345-row table replayed byte-identically every run; the 2 recorded hits lead-replayed) |
| Total witness coverage | **353/353 — COMPLETE** | union of the above: $347 + (31,-2) + 4$ residual cells + $(89,-2)$ via L13f | **PROVED per row** (every cell carries its own replayable certificate; the grid is closed since 2026-08-18) |
| Certified emergent-free members | **24** earlier; **293 class-complete L14** | `data/l12_density.json` (9 of 44 factored members, sweep 1) + `data/l13_zerobad_more.jsonl` (15 of 57, sweep 2); even-only histograms — superseded on the grid by the L14 row below | **EVIDENCE** (positive empirical density of step (ii)); 15 distinct cells, 16 frozen rows `h10q.py::_L13_ZERO_BAD` replayed each run — **PROVED per row** |
| Emergent $-1$s are real | 467/959 spot-probes | `data/l12b_class_sample.json` | **VERIFIED-instance** (step (ii) is genuinely outside the class claim) |
| No structural shortcut | 706/706 rows irreducible deg 8 | `THEOREMS.md` L13c; `data/l12_squarehunt.json`; `data/l12_factorP.json` | **PROVED per row** (sympy `factor_list`; in-suite Frobenius re-certification) + **EVIDENCE** (137.5M square-class tests, 0 hits; Galois $D_8\wr C_2$ consistent, 28,240 certified samples) |
| Cell $(89,(-2,1))$ — **CLOSED 2026-08-18** | 325 certified soluble candidates | `data/l13_filter_run2.json` (236 in the structured box: 133 proved-prime-rung + 103 factorint-rung) + dense-box block (89: 46 + 43); witness $(a,b)=(3,\,89/367)$, $\tau=19/37$, frozen `_L13_RESIDUAL` row 7 | **PROVED** (L13f cofactor ladder: zero-bad ⟺ stripped remainder is a square OR a single proved prime — parity law forces its symbol $+1$; lead-replayed all five rungs; INCONSISTENT = 0 across 5,767 decided rows; 1,896 + 1,637 + 937 refusals recorded with reasons, never evidence) |

The grid is corroboration, not the source of the uniform implication.
L20 supplies a class for every cell, L22 supplies fixed-fibre
irreducibility for every nonzero rational $Z$, and L19 supplies a member
from classical Schinzel H (`THEOREMS.md`, L19–L22).  No count in this
table is promoted beyond its cited finite scope.

---

## 4. The closed non-Schinzel gap and the remaining conditional input

### 4.1 What L22 closed

The former algebraic obstacle was vertical irreducibility of the
degree-$8$ class polynomial on the **fixed**
$$
\tau^\dagger=(1+2a^2)/A,\qquad
\delta_{\tau^\dagger}=-4a^4/A
$$
branch.  L22a proves it for every fixed nonzero rational $Z$ by
$A$-adic parity, the $a=0$ and $a=\infty$ endpoint allocations, and
exact elimination of degrees $2$ and $4$
(`/tmp/l22_elimination.md`; `data/l22_elimination.jsonl`).  The
exception $Z=0$ is explicitly reducible and outside the cell family;
$Z=1$ is handled by an exact $a=1$, mod-$11$ certificate
(`/tmp/l22_elimination.md`, §5).

L22b is independent.  At
$a=1,A=5,s=0,\tau=3/5=\tau^\dagger,\delta=-4/5$, the reciprocal trace
has a uniform mod-$8$ norm obstruction, and the lift reduces to the
rank-$0$ curve $y^2=x(x+5)(x+80)$
(`/tmp/l22_reciprocal_cube.md`;
`data/l22_reciprocal_cube.jsonl`).  This $a=1$ value is only a
specialization witness; quantitative HIT chooses the actual admissible
$a$ in L20's progression (`THEOREMS.md`, L21d and L22b).

The other L22 routes retain their exact labels:

- higher-order infinity analysis proves four quadratic local clusters
  and stops at possible factor degrees $2,4,6$; L22a eliminates them
  globally (`/tmp/l22_infinity.md`; `data/l22_infinity.jsonl`);
- the factor-tuple criterion is PROVED and its prime-value conclusion is
  CONDITIONAL, but its no-go shows parity controls only the product of
  outside signs (`/tmp/l22_factor_tuple.md`;
  `data/l22_factor_tuple.jsonl`);
- the free-$\lambda$ class theorem is PROVED algebraically but **NOT
  APPLICABLE** to this count, because the selected branch runs through
  the infinite L11c family and would require an additional witness
  (`/tmp/l22_square_branch.md`; `agent://BranchAudit`);
- the Noether-pencil descent is PROVED, while one explicit polynomial
  nonvanishing lemma remains OPEN in that auxiliary route.  It is no
  longer a premise of Theorem C (`/tmp/l22_fiber_geometry.md`;
  `data/l22_fiber_geometry.jsonl`).

Thus the non-Schinzel chain is closed:
$$
\text{fixed cell}\xrightarrow{\mathrm{L22+HIT}}
\text{irreducible admissible }a
\xrightarrow{\mathrm{L20}}\text{aligned class}
\xrightarrow{\mathrm{L19+Schinzel}}\text{member}.
$$
(Source: `THEOREMS.md`, L19–L22.)  Conditions (a)–(d) are all PROVED;
the final arrow is CONDITIONAL only because classical Schinzel H is
unproved.

### 4.2 Literature provenance, not a remaining irreducibility premise

The earlier literature search remains useful provenance for why the
Schinzel input has not been replaced unconditionally.  It is **not** a
claim that a degree-$8$ irreducibility certificate is still missing.
All bibliographic details below were source-checked in
`data/litscout_h10q.md`.

- For integer arguments, unconditional squarefree-value theorems cover
  degrees at most $3$: Erdős, *J. London Math. Soc.* 28 (1953),
  416–425; Hooley, “On the power-free values of polynomials,”
  *Mathematika* 14 (1967), 21–26; and, for binary forms, Greaves,
  “Power-free values of binary forms,” *Q. J. Math. Oxford* (2) 43
  (1992), 45–65 (`data/litscout_h10q.md`).
- Under $\mathsf{abc}$, Granville, “ABC allows us to count
  squarefrees,” *IMRN* 1998 (19), 991–1009, handles arbitrary degree at
  integer arguments; Poonen, *Duke Math. J.* 118 (2003), 353–373,
  gives the multivariable analogue (`data/litscout_h10q.md`).
- At prime arguments, the recorded anchors are Baker–Pollack,
  “Clusters of primes with square-free translates,”
  pollack-math.net/ClustersPrimesSquarefree.pdf, in degree $1$;
  Helfgott, “Square-free values of $f(p)$, $f$ cubic,”
  arXiv:1112.3820, in degree $3$; and Reuss, “Power-Free Values of
  Polynomials,” arXiv:1307.2802, for $(d-1)$-free values
  (`data/litscout_h10q.md`).
- For prescribed square class in a progression, Krumm, “Squarefree
  parts of polynomial values,” *JTNB* 28 (2016), 699–724,
  arXiv:1407.4890, Theorem 1.3 and Propositions 3.4–3.5, proves the
  degree-at-most-$2$ case.  Proposition 3.8 treats degree $3$
  conditionally on the elliptic-curve Parity Conjecture, not the sieve
  parity problem (`data/litscout_h10q.md`).
- For the conic-bundle provenance, Harpaz–Skorobogatov–Wittenberg,
  *Compositio Math.* 150 (2014), 2095–2111, arXiv:1304.3333,
  Theorem 3.1 and Corollary 3.4, is unconditional under its abelian
  constant-field hypothesis.  The general route uses Schinzel H in
  Colliot-Thélène–Skorobogatov–Swinnerton-Dyer, “Rational points and
  zero-cycles on fibred varieties: Schinzel's hypothesis and
  Salberger's device,” *J. reine angew. Math.* 495 (1998), 1–28, and
  the Colliot-Thélène–Sansuc Schinzel paper in *Acta Arith.* 41 (1982),
  33–53 (`data/litscout_h10q.md`).
- Carella arXiv:2310.16952 remains excluded from the record for the
  reasons documented in the source-checked literature report
  (`data/litscout_h10q.md`).

### 4.3 L23–L27: unconditional frontier theorems, but no member theorem

L23 now proves a genuine unconditional theorem for the selected
linear/octic sequence.  After absorbing the finitely many exceptional
primes into the progression, the bad-root sieve has exact dimension
$\kappa=1/2$ and Bombieri--Vinogradov level
$$D=X^{1/2}/(\log X)^B.$$
The semilinear lower sieve with $z=X^{49/100}$ gives
$$
\#\{X<t\le2X:\ Q(t)\ {\rm prime},\
G(t)\ {\rm has\ no\ bad\ root\ prime}<z\}
\gg\frac{X}{(\log X)^{3/2}}.
$$
This is unconditional **small-prime cleanliness**, not intermediate H:
$|G(t)|=X^{8+o(1)}$ can retain up to $16$ factors above $z$, and
reciprocity permits $0,2,\ldots,16$ bad odd-valuation factors.

The first missing analytic input is therefore a parity-sensitive
two-large-bad-divisor estimate in the sector
$$p_1,p_2\ge z,\qquad p_1p_2>D,$$
with both primes of bad sign and odd valuation.  Equivalently one needs
$R_{\rm bad}\le1$ for sign-decorated divisors of the octic at prime
arguments.  Even Elliott--Halberstam leaves up to eight large factors;
the strongest cited general almost-prime result gives only $P_{12}$
after a fixed-AP adaptation (`THEOREMS.md`, L23a–b;
`data/l23_half_sieve.jsonl`).

The algebraic alternatives tested in L23 do not bypass that sector.
The untwisted Capell square-in-$K$ collapse is impossible on the
selected refined protocol; the associated fixed-cell conic bundle has
non-split rank $10$ and
$\operatorname{Br}(X)/\operatorname{Br}(\mathbb Q)=\mathbb Z/2$, outside
the checked low-rank fibration theorems.  Exact norm matching and
polynomial-in-$b$ sections miss $\Phi$; the elementary
$P=\pm2b\,\square$ and $P=A\,\mathrm{Norm}$ routes hit dyadic walls;
the first admissible factor ansatz reduces to genus-$4$ curves; and
moving $a$ has smallest **certified** squareclass degree $18$ in the
audit, with no universal minimum claim (`THEOREMS.md`, L23c–f).

L24 closes the tempting self-coupled diagonal alternative even more
sharply.  Both formal tied-rank-$\le1$, finite-flat degree-$8$ covers
are uniformly empty on $\Phi$ over $\mathbb Q_2$.  Orientation I has
normalized coefficient valuations $(4,0,0)$; orientation II reduces
to the complete contradiction
$$0\equiv8+2Ab(s^2-1)\pmod{16}.$$
Thus the diagonal $\le5$ count is vacuous and supplies no member.  The
$4{,}026{,}282$ zero-hit scan is corroboration only; the mod-$16$
exhaustion is the proof (`THEOREMS.md`, L24;
`data/l24_diagonal_{geometry,arithmetic,local,search}.jsonl`).

L25–L27 supply genuinely non-diagonal local and algebraic advances.
L25 proves that two fixed scalings cover complementary dyadic parities,
but its target formula is not aligned.  L26's even-pullback reciprocal
tie gives
$$P_{\rm rec}(b)=b^4T(b+b^{-1}),\qquad\deg T=4,$$
automatic dyadic splitting, and the trace-character identity
$(2b|p)=(2(u+2)|p)$; the quartic squareclass and global member remain
open.  L27's fixed shear $(2X+28\rho,sX+\rho)$ covers every dyadic
parity and every standard aligned odd target, including $3,5$, but its
global degree-$8$ cover has no proved rational point
(`THEOREMS.md`, L25–L27;
`data/l25_scaled_coupling.jsonl`;
`data/l26_reciprocal_tie.jsonl`;
`data/l27_triangular_shear.jsonl`).

Consequently classical Schinzel H has **not** become proved, weaker, or
removable in Theorem C.  The exact next targets are the L26 quartic
trace squareclass/member problem, a target-specific rational point on
the L27 global cover, or the two-large-bad-divisor dispersion estimate.
Theorem C remains **CONDITIONAL on classical Schinzel H alone**.

---

## 5. Provenance and discipline

- Authorities: `h10q.py` frozen tables; `THEOREMS.md` L6–L27;
  `RESULTS.md`; and the cited `data/` artifacts.  The L22 proof reports
  remain `/tmp/l22_elimination.md` and
  `/tmp/l22_reciprocal_cube.md`; all L23–L27 producer/data pairs are
  inventoried in `README.md`.
- Engine discipline remains proven-primality only.  The exact
  Miller–Rabin, Pocklington and Brillhart–Lehmer–Selfridge scopes are
  recorded in `README.md`, Engine note, and `test_bls.py`; refusals are
  logged and never treated as evidence.
- The HalfSieveAudit, FibrationAudit and DiagonalAudit verdicts were
  SOUND at confidences $0.90$, $0.88$ and $0.90$ respectively.
  ReciprocalTieAudit's three scope findings were applied; the corrected
  TriangularShearAudit returned no high-confidence finding.  Their
  corrections are incorporated in `THEOREMS.md`, L23–L27 and retain all
  finite rows at their declared scopes.
- The $4{,}026{,}282$ diagonal attempts and every other bounded no-hit
  scan are EVIDENCE only.  Uniform diagonal emptiness is the dyadic
  theorem, not extrapolation from the scan.
- Nothing here depends on Sun's unrefereed §§3–8 chain.  Soundness is
  inherited from the audited block; completeness of the conditional
  implication uses L19–L22 and classical Schinzel H.  L23–L27 sharpen
  the frontier but are not replacements for that conjectural input.
- **Final status:** the six-universal-quantifier record and
  $\operatorname{efd}\le5$ are established conditionally on classical
  Schinzel H alone (`THEOREMS.md`, L22d and L24e–L27).  Classical
  H is unproved, so H10/$\mathbb Q$ and the unconditional six-quantifier
  statement remain open.
