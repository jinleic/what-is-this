# Precise statements at the frontier (session output, 2026-08-12)

Status labels: [PUBLISHED] = located in print with citation; [FOLKLORE?] = we prove it here,
no published statement located (priority search: Mazur Astérisque 228, Cornelissen–Zahidi
math/0006140, Koenigsmann 1309.0441, Daans 2301.02107, DDF21); [AUDIT] = our analysis of a
published proof.

---

## L1 [PUBLISHED — Cornelissen–Zahidi]. No one-witness definition of Z in Q

**Statement.** There is no $f \in \mathbb{Q}[t, y]$ with
$\mathbb{Z} = \{t \in \mathbb{Q} : \exists y \in \mathbb{Q},\, f(t,y) = 0\}$.
More generally, $\mathbb{Z}$ is not the coordinate projection of $X(\mathbb{Q})$ for any curve
$X/\mathbb{Q}$.

**Source.** Cornelissen–Zahidi, arXiv math/0006140: Example 2.2(a) proves *Mazur's Conjecture 2.1
is true for curves* (genus $\ge 2$: Faltings; genus 0: $\mathbb{P}^1(\mathbb{Q})$ dense in
$\mathbb{P}^1(\mathbb{R})$; genus 1: closed infinite subgroups of $E(\mathbb{R}) \cong S^1$ or
$S^1 \times \mathbb{Z}/2$ have finite index), and Remark 2.4 identifies existentially definable
sets with projections of rational points. Projection from a smooth projective model is a closed
map on compact real loci, so the closure of the projected rational points has finitely many
connected components; $\overline{\mathbb{Z}} = \mathbb{Z}$ has infinitely many.

**Consequence for the stratification.** The smallest open case of "is $\mathbb{Z}$ existentially
definable in $\mathbb{Q}$?" is **two witness variables / witness varieties of dimension 2** —
exactly where Mazur's conjecture (real topology of surfaces) and strong Bombieri–Lang
(Koenigsmann 1011.3424, Cor 23) are the only known obstructions, both conjectural.

**Numerical illustration** (in `h10q.py`): on $E: y^2 = x^3 - 2$, the elliptic-log positions
$\{n\alpha\}$, $\alpha = 0.2682929110$, of $nP$ fill all 24 bins of $E(\mathbb{R})$ — closure has
1 component; contrast $\overline{\mathbb{Z}} = \mathbb{Z}$.

---

## L2′ [REPAIRED 2026-08-12; FOLKLORE?]. Cofinite-S structure theorem (model-relative)

**Setting.** $E/\mathbb{Q}$ an elliptic curve given by a **fixed** Weierstrass equation whose
coefficients $a_1,\dots,a_6$ are integral at every prime $q \in Q_0$, where $Q_0$ is a finite set
of primes, $S$ = its complement, $R = \mathbb{Z}[S^{-1}]$ (so $x \in R \iff v_q(x) \ge 0\ \forall
q \in Q_0$). **No minimality is assumed.** $P \in E(\mathbb{Q})$ of infinite order. All statements
are relative to this model; see the model-dependence remark. Let
$$D = \{n \ge 1 : x(nP) \in R\}.$$

**Lemma (integral-model formal group).** For each $q \in Q_0$,
$$G_q := \{T \in E(\mathbb{Q}_q) : v_q(x(T)) < 0\} \cup \{O\}$$
is an **open subgroup of finite index** in $E(\mathbb{Q}_q)$.

*Proof.* If $v_q(x(T)) < 0$ then comparing valuations in the Weierstrass equation (all $a_i$
integral) forces $v_q(x) = -2k$, $v_q(y) = -3k$, $k \ge 1$; so $z = -x/y$ has $v_q(z) = k \ge 1$
and $G_q$ is the image of $\hat{E}(q\mathbb{Z}_q)$ under the standard parametrization. The formal
group law $F(z_1, z_2)$ has coefficients in $\mathbb{Z}[a_1,\dots,a_6]$ (Silverman IV.2.3 — the
computation uses only integrality of the $a_i$; minimality enters nowhere), so
$\hat{E}(q\mathbb{Z}_q)$ is a group and $G_q$ is a subgroup. It is open: at affine points
$v_q(x) < 0$ is an open condition, and the filtration sets $\{v_q(x) \le -2k\} \cup \{O\}$ are a
neighborhood basis of $O$. Since $E(\mathbb{Q}_q)$ is compact (closed in
$\mathbb{P}^2(\mathbb{Q}_q)$), an open subgroup has finite index. $\square$

**Theorem (L2′).** Let $m_q \ge 1$ be the (finite) order of $P + G_q$ in $E(\mathbb{Q}_q)/G_q$,
and $M = \mathrm{lcm}_{q \in Q_0}\, m_q$. Then:

(a) *Exact apparition.* $\{n \in \mathbb{Z} \setminus \{0\} : v_q(x(nP)) < 0\} =
    m_q\mathbb{Z} \setminus \{0\}$, for each $q \in Q_0$.

(b) *Structure.* $D = \{n \ge 1 : m_q \nmid n\ \forall q \in Q_0\}$ — a union of congruence
    classes mod $M$, computable from the finite data $\{m_q\}$; membership in $D$ is decidable,
    and $D$ is (eventually) periodic.

(c) *Emptiness criterion.* $D = \emptyset \iff$ some $m_q = 1 \iff x(P) \notin R \iff 1 \notin D$.
    If $D \ne \emptyset$ then $D \supseteq \{n \ge 1 : n \equiv \pm 1 \pmod{M}\}$, an infinite
    union of arithmetic progressions of positive density.

(d) *Equidistribution.* If $D \ne \emptyset$, then $\{nP : n \in D\}$ contains
    $\{(1+jM)P : j \ge 0\}$, which is dense in the coset $P + E(\mathbb{R})^0$ (Weyl, elliptic
    log irrational since $P$ is non-torsion); hence $\{x(nP) : n \in D\}$ is dense in the
    $x$-range of a full connected component of $E(\mathbb{R})$ — never discrete.

*Proof.* (a) If $m_q \mid n$, $n \ne 0$, then $nP \in G_q \setminus \{O\}$ ($P$ non-torsion), so
$v_q(x(nP)) < 0$. Conversely $v_q(x(nP)) < 0$ means $nP \in G_q$, and
$\{n : nP \in G_q\}$ is the kernel of $\mathbb{Z} \to E(\mathbb{Q}_q)/G_q$, $= m_q\mathbb{Z}$.
$m_q$ is finite because the quotient is finite (Lemma). (b) is (a) intersected over $Q_0$.
(c) $m_q = 1 \iff P \in G_q \iff v_q(x(P)) < 0$; and if all $m_q \ge 2$ then
$n \equiv \pm 1 \pmod{M}$ misses every $m_q\mathbb{Z}$. (d) as in the previous version:
$\overline{\langle MP \rangle} \supseteq E(\mathbb{R})^0$, Weyl equidistribution of
$\{n\alpha\}$. $\blacksquare$

**Remark (model-dependence — repaired).** The earlier version of this section asserted that a
coordinate change $u$ supported on $Q_0$ "does not change $R$-integrality of $x$". **That claim is
false and is retracted** (advisory, 2026-08-12): under $x = u^2 x' + r$, any $v_q(u) \ne 0$ shifts
$v_q(x)$ by $2v_q(u)$ at deep points, and $q \in Q_0$ is precisely not inverted in $R$. $D$ is
genuinely model-relative: rescaling by $v_q(u) = -1$ replaces $G_q$ by a deeper filtration
subgroup and multiplies $m_q$ by a power of $q$. Concrete verified instance (h10q.py,
`_verify_L2prime`): $y^2 = x^3 - 2$, $P = (3,5)$ vs. the 5-rescaled model
$Y^2 = X^3 - 2\cdot 5^6$, $P' = (75, 625)$ — both have 5-integral coefficients, and with
$Q_0 = \{5\}$ the sets $D$ differ (the rescaled model's $m_5$ picks up a factor 5). The *theorem*
— openness, finite index, exact APs, decidable structure, equidistribution — holds for **every**
$Q_0$-integral model; only the specific moduli $m_q$ (hence $D$) depend on the choice. Nothing is
gained by choosing the global minimal model beyond canonicity.

**Remarks (unchanged in substance).**
- *Why this is the exact cliff.* Poonen's construction (math/0306277) needs
  $\{n : nP \in E(\mathbb{Z}[S^{-1}])\}$ to be (essentially) the primes — density 0, spread out,
  archimedean-discretizable. L2′(b) says that with only finitely many primes withheld the
  integral index set is a boolean combination of congruences — positive density, equidistributed,
  **decidable**. No Poonen-style model can exist over cofinite $S$: an infinite excluded prime
  set is necessary for every encoding of this shape. The known undecidability results (density-1
  $S$, complementary partitions) are not "one prime away" from $\mathbb{Q}$ — the mechanism
  degenerates qualitatively.
- *Scope.* L2′ constrains the raw $x$-integrality predicate on $\langle P \rangle$. It does not
  rule out discrete Diophantine sets over $\mathbb{Z}[S^{-1}]$, $S$ cofinite, via *other*
  predicates — that stronger statement is a Mazur-conjecture variant for these rings and remains
  open. Under Mazur/CZ for $\mathbb{Q}$ no such set exists over $\mathbb{Q}$ itself.
- *Priority.* No published statement located in this exact form (searched: Mazur 1995, CZ 2000,
  Poonen 2003, Koenigsmann 2014, Daans 2024, DDF21). The ingredients are standard; treat as
  expository formalization pending a deeper search.
- *Verified computationally* (h10q.py, `_verify_L2prime`): exact-AP claim (a) for both models at
  every prime $q < 50$ appearing; the congruence formula (b) against brute force; the emptiness
  criterion (c); and the model-dependence exhibit above.

---

## A1 [AUDIT]. Anatomy of the 10-quantifier record (Daans, JLMS 2024, arXiv 2301.02107)

**The formula** (Thm 5.6 via Lemma 5.5, Eq (5.46)): for suitable $\pi, u$ ($S = \mathrm{Odd}(\pi)$,
$|S|$ odd, dyadics in $S$, $u$ a unit at $S$ with $X^2 - X - u^2$ irreducible at all $v \in S$),
$$x \in \bigcup_{v \notin S} \mathfrak{m}_v \iff \exists a, b: (a,b) \in \Phi_u^S \wedge
\frac{a^2 x^2 g(a,b)}{1 - x - a^2 x^2} \in \bigcap_{v \in \Delta[a^2, b\pi)_K} \mathcal{O}_v,$$
then $\mathcal{O}_S$ = complement, universalized by negation (Lemma 5.1 + §2.4).

**Exact budget.** Naively **12** existential variables:
| role | count | source |
|---|---|---|
| quaternion parameters $a, b$ | 2 | Eq 5.46 |
| $\Phi_u^S$ membership: $(b^2{+}1)/b,\ (a{-}u)/\pi \in \bigcap_{v\in S}\mathcal{O}_v$ | 3 | Lemma 5.2 (unit trick Prop 4.6 packs the pair) |
| splitting witness $y$ (Thm 4.9: $z \in \bigcap_\Delta \mathcal{O}_v \iff \exists y,\ y, z{-}y \in S(Q)$) | 1 | Eq 4.26 |
| trace-norm certificate for $y \in S(Q)$ — Prop 3.6(iii): $x_1^2 + x_1(c-2x_1) - a(c-2x_1)^2 - b(x_2^2 + x_2x_3 - a x_3^2) = d$ | 3 | Prop 3.6 |
| same for $z - y \in S(Q)$ | 3 | Prop 3.6 |

**Two fusions** via DDF21 Thm 1.4 (over fields finitely generated over a perfect subfield,
$D_1 \cap D_2$ existential with $m_1 + m_2 - 1$): once inside Prop 4.10
($1 + (3 + 3 - 1) = 6$, uniform in $(a,b)$), once at top level ($2 + (3 + 6 - 1) = 10$).
Iterating pairwise fusion on the three quadratic-form blocks gives $12 - 2 = 10$ under ANY
association order — the count 10 is the fixed point of this proof shape.

**Slack channels (each would give ≤ 9):**
1. **Stronger fusion.** A $k$-fold simultaneous fusion better than iterated pairwise
   ($m_1{+}m_2{+}m_3-3$ instead of $-2$) for quadratic-form-type blocks. DDF21 Thm 1.4's proof is
   the place to look; whether the fused pair can share MORE than one coordinate is not addressed
   in Daans.
2. **Cheaper primitive.** A 2-witness certificate for $s \in S(Q)$ (currently 3 via the ternary
   norm form; a conic parameterization argument would need to handle the anisotropic locus), or a
   pair-membership $(z_1, z_2) \in (\bigcap_\Delta \mathcal{O}_v)^2$ at cost $\le 5$ uniform.
3. **Restructured $\Phi$.** Since $S \subseteq \Delta[a^2, b\pi)_K$ for $(a,b) \in \Phi_u^S$, the
   three ring-membership conditions live in nested intersections; a single combined membership
   statement for a triple of rational functions could collapse $\Phi$ into the main block.

**Lower bound state.** $m \ge 2$: Thm 2.6 ($K \setminus R$ is never $\exists_1$ for infinite
proper subrings — thin-set/Hilbert irreducibility). DDF21 §8 (quoted by Daans §2):
*"we do not have any example of an ∃-definable subset of a global field which we can show is not
$\exists_2$-definable."* So **Question 5.7 (Daans): $2 \le m \le 10$**, and BOTH ends are walls:
$m \ge 3$ needs lower-bound technology beyond thin sets; $m \le 9$ needs one of the slack channels.

**The dimension-2 wall, unified.** Three independent lower-bound technologies all stop at 2
witnesses: (i) real topology of curves (Faltings + Lie groups) — kills $\exists_1$ definitions of
$\mathbb{Z}$ [CZ]; (ii) thin sets / Hilbert irreducibility — kills $\exists_1$ complements
[DDF21/Daans 2.6]; (iii) nothing at all is known against $\exists_2$ in either direction, and the
conjectural obstructions (Mazur for surfaces, strong BL) are exactly statements about
2-dimensional witness varieties.

---

## Consequences table († = Sun 2607.28606, unrefereed; chain audited 2026-08-12, no substantive gap found — A2)

| Object | Status |
|---|---|
| $\mathbb{Z}$ in $\mathbb{Q}$: $\forall_m$ | $2 \le m \le 7$† ($\le 10$ refereed: Daans, JLMS 2024) |
| $\mathbb{Z}$ in $\mathbb{Q}$: $\exists_k$ | $k = 1$ impossible [CZ]; $k \ge 2$ fully open; any $k$ ⟹ H10/$\mathbb{Q}$ undecidable |
| $\mathbb{Q}\setminus\mathbb{Z}$: $\exists_k$ | $k = 1$ impossible [Daans 2.6]; $k = 7$† ($k = 10$ refereed) |
| $\forall_9\exists_m$-theory of $\mathbb{Q}$ | undecidable for $m = 7$† ($m = 10$ refereed: Daans Cor 6.2 + Sun 2021) |
| r.e. subsets of $\mathbb{Q}$ | $\exists_{10}\forall_{10}$-definable [Daans Cor 6.2]; $\forall$-part $\to 7$† |
| H10/$\mathbb{Q}$ (full existential theory) | open in both directions |

---

## A2 [AUDIT, revised 2026-08-12]. The claimed record 10 → 7 (Sun, arXiv 2607.28606)

**Claim** (posted 2026-07-30, unrefereed): $\bigcup_{v \notin S_0}\mathfrak{m}_v$ is Diophantine
with **7** unknowns over any global field; hence $\mathbb{Z}$ is $\forall_7$-definable in
$\mathbb{Q}$ and (with Sun 2021) the $\forall_9\exists_7$-theory of $\mathbb{Q}$ is undecidable.
Acknowledgment in the paper: *"This approach was initially motivated by AI's analysis of Daans's
paper."*

**Architecture** (as posted). Daans' bridge (his Lemma 5.5 = Sun Eq. 9.2) is kept verbatim:
$z \in \bigcup_{w \notin S}\mathfrak{m}_w \iff \exists (a,b) \in \Phi_u^S\,
[h(a,b,z) \in \bigcap_{v \in \Delta(Q_{a,b})}\mathcal{O}_v]$, with
$\Phi_u^S = \{(a,b) : b \in \bigcap_{v\in S}\mathcal{O}_v^\times,\ a - u \in \bigcap_{v\in S}
\mathfrak{m}_v\}$ ($\exists^3$, Daans Lemma 5.2), $Q_{a,b} = [a^2, b\pi)_K$, and
$h(a,b,z) = a^2z^2g(a,b)/(1 - z - a^2z^2)$. Daans' $\exists_6$ intersection block is replaced by
$\Theta(a,b,c) = \bigvee_{\tau \in \Lambda}\Psi_\tau(a,b,c)$ where, in char $\ne 2$
($A = 1+4a^2$, $B = b\pi$, $\delta_\tau = 1 - A\tau^2$):
$$\Psi_\tau(a,b,c):\ \exists y,r,s\ [\delta_\tau \ne 0 \wedge
\delta_\tau(c^2 - Ay^2) - 16B(r^2 - As^2) = 16],$$
$\Lambda$ a **finite set of frozen constants** produced by weak approximation (his §§3–4). A
solution yields $\alpha = z + wj$ with $\mathrm{Nrd}(\alpha) = 1$ and $\mathrm{Trd}(\alpha) +
\mathrm{Trd}(\alpha\gamma_\tau) = c$ (identities 2.2–2.4), so $c$ is integral at every ramified
place (Lemma 2.1) — soundness. $A = 1+4a^2 > 0$ splits the real place. All $\Psi_\tau$ share the
same 3 witnesses, so $\Theta$ is $\exists^3$; count $2 + (3 + 3 - 1) = 7$ with the $-1$ from DDF
Thm 1.4 (non-constructive — the 7-unknown polynomial is an existence statement, like Daans' 10).
His §8 = an affine-point lemma (8.1: a nondegenerate quadric minus a hyperplane keeps a point
over an infinite field) + classical **Hasse–Minkowski for the resulting quaternary form**
(O'Meara 66:1; char 2 via Pollak/Wu). Completeness is *family-level*: for $x \in \mathfrak{m}_w$,
$w \notin S \cup E_{\mathrm{exc}}$, §§5–7 must select $(a,b) \in \Phi_u^S$ (character-sum
selection at $w$, Weil bound; global norm approximation) making some frozen $\Psi_\tau$ locally
solvable everywhere; $E_{\mathrm{exc}}$ finite, patched by $\exists^3$-definable
$\mathfrak{m}_v$'s.

**Computationally corroborated** (h10q.py `_verify_sun`; exact-arithmetic Hasse–Minkowski
engine over a **proven-primality** factoring layer — deterministic MR (A014233, bases 2–41,
exclusive bound $\psi_{13}$) + Pocklington $n{-}1$ certificates, no probabilistic acceptance;
$K = \mathbb{Q}$, $S = \{2\}$, $\pi = 2$, $u = 1$; all scopes exact, all bounded):
- (A) identities (2.2)–(2.4), 150 random exact instances;
- (B) **soundness** of the block: $\Psi_\tau$ solvable $\Rightarrow$ $c$ integral at every finite
  ramified place — 0 violations across all sampling (≈3000 instances this session);
- (C) **fixed-$\tau$ incompleteness**: $(a,b,c) = (2, -7/3, -9)$, $\Delta = \{3,7\}$, $c$
  $\Delta$-integral, yet his explicit rational specializations (9.5) and (9.6) are BOTH
  unsolvable — so the $\tau_0/\tau_1$ remark cannot replace the $\Lambda$-family: the §§3–4
  freezing is load-bearing, and any reading of $\Theta$ as (9.5)∨(9.6) breaks the theorem;
- (D) Lemma 2.2's residue-field combinatorics, exhaustive over $\mathbb{F}_q$, all $q < 50$;
- (E) bridge existence: explicit frozen witnesses $(a, b, \tau)$ with $\Delta = \{2, w\}$ for
  **every odd prime $w < 300$** (61 targets; $w < 100$ re-verified on every run, full table via
  `--extended`, 8.1 s). **Eight targets ($w = 137, 139, 173, 181, 241, 251, 269, 281$) required
  a frozen $\tau \notin \{0, \tau_1\}$** within our $(a,b)$ grids (e.g. $\tau = 1/4$ five
  times) — the $\Lambda$-family license is exercised by real targets, reinforcing (C).

These items corroborate the chain at its mechanically accessible points; **the local–global
completeness proof itself (§§3–8) is hand-traced, not machine-verified** — see below.

**Completeness chain — traced (revision 2 of this audit).** The posted text has **stale
cross-references** at the decisive step: Prop 9.1's proof says *"Section 11 then produces a
global point"* (the paper ends at §9) and §7/§9 cite *"Sections 7 and 8"* for the target-place
local point and the exceptional set. These are **editorial defects** (draft renumbering), not
missing mathematics: the chain is present and traceable in the posted sections —
- **§3** (Prop 3.1, parameterized Hensel): finite freezing data at each $v \in S$ — smooth local
  points for the fixed-parameter equation on neighborhoods $U_{v,i}$;
- **§4**: weak approximation assembles the finite global family $\Lambda$; $\Theta$ is $\exists^3$;
- **§5** (target place, char $\ne$ 2): Lemma 5.1 (character sum, Weil bound, needs
  $|\kappa_w| > 25$) selects $\bar a_w$ with $A_w$ a nonsquare and $-A_w(1 - A_w\bar\tau^2)$ a
  **square**; weak approximation gives global $a$; with $v_w(b) = 1$, Daans' valuation
  computation yields $v_w(c) = 6m - 2 \ge 4$ (5.12), the equation reduces mod $\mathfrak{m}_w$
  to $-\bar A\bar\delta_\tau\bar y^2 = 16$ — solvable by the square-ness, unit derivative,
  **Hensel lifts: the local point at $w$** (the intended target of "Sections 7 and 8");
- **§7**: global $(a,b)$ via weak + norm approximation (Lemma 7.1), keeping $(a,b) \in \Phi_u^S$,
  $v_w(b) = 1$, ramification confined to $S \cup \{w\}$; split places get points from the
  standard matrix construction;
- **§8**: local points everywhere ⟹ Hasse–Minkowski isotropy of the quaternary form + Lemma 8.1
  (affine-point extraction, $\ell \ne 0$) ⟹ a global $\Psi_\tau$-solution — the intended
  "Section 11".

**Verdict (corrected).** No substantive missing implication found by this audit; the defects are
editorial. Computational corroboration at every mechanically accessible point (items A–E above).
Referee-level items that remain are ordinary verification of written content: §3's dyadic
freezing case, §7's simultaneous norm-approximation constraints, and a v2 fixing the
cross-references. **Status: 7 = established modulo refereeing (unrefereed, editorially
defective); refereed anchor: 10 (Daans, JLMS 2024).** Question 5.7 range: $2 \le m \le 7$,
with $\le 10$ the refereed bound.

**Remaining slack if 7 stands:** (i) collapse $\Phi_u^S$ into the trace-sum block (A1 channel 3,
worth $\le 6$); (ii) an $\exists_2$ certificate for either block — the naive route is now
**certified dead**: freezing $s = 0$ in $\Psi_\tau$ stays sound but provably loses completeness
(24 of 224 $\Delta$-integral instances on the witness grid, first at $(a,b,c,\tau) =
(1,-3,2,0)$; h10q.py `_ternary_probe`) — a 2-witness certificate cannot be a specialization of
$\Psi$; (iii) better-than-pairwise fusion (constructive account announced in Becher–Daans).
Lower bound unchanged: $m \ge 2$; nothing is known against $\exists_2$ (DDF §8; their best
positive example: $\mathrm{rk}^\exists(K^{(2)} + K^{(2)}) = 2$ for Hilbertian real $K$, Cor 8.11).

---

## L4 [PROVED HERE, 2026-08-12]. Sun's Lemma 5.1: the constant 25 is really 3

**Sun's Lemma 5.1** (arXiv 2607.28606). Let $k$ be a finite field, $\mathrm{ch}(k) \ne 2$,
$\chi$ the quadratic character ($\chi(0) = 0$). For any $\tau \in k \setminus \{0, \pm 1\}$
there exists $a \in k$ with
$$\chi(1 + 4a^2) = -1 \quad\text{and}\quad \chi\bigl(-(1+4a^2)\,(1 - \tau^2(1+4a^2))\bigr) = 1,$$
**provided $|k| > 25$** (character-sum/Weil estimate). Accordingly his $E_{\mathrm{exc}}$
includes *all* places with $|\kappa_w| \le 25$.

**L4 (exhaustion; the exhaustive check over a fixed finite field IS a proof for that field).**
For every odd prime power $q \le 25$ and **every** admissible $\tau \in \mathbb{F}_q \setminus
\{0, \pm 1\}$, the required $a$ exists — with one structural exception:
- $q = 3$: **vacuous** — $\mathbb{F}_3 = \{0, \pm 1\}$ admits no $\tau$ at all;
- $q \in \{5, 7, 9, 11, 13\}$: holds, min. solution count 2;
- $q \in \{17, 19, 23, 25\}$: holds, min. solution count 4;
- sanity $25 < q \le 49$ ($q = 27, 29, 31, 37, 41, 43, 47, 49$): holds for all $\tau$,
  consistent with the Weil range.

(Verified by h10q.py `_lemma51_exhaust`: prime fields directly, $\mathbb{F}_9 =
\mathbb{F}_3[x]/(x^2{+}1)$, $\mathbb{F}_{25} = \mathbb{F}_5[x]/(x^2{+}2)$, $\mathbb{F}_{27} =
\mathbb{F}_3[x]/(x^3{+}2x{+}1)$, $\mathbb{F}_{49} = \mathbb{F}_7[x]/(x^2{+}1)$; 0.02 s.)

**Consequence.** Superseded by **L5** (next section): the table above is the finite base
data; L5's exact identity $4N(\tau) = \#E_{\tau^{-2}}(\mathbb{F}_q)$ proves the selection
for *every* odd prime power at once and voids Sun's threshold entirely. The min counts
observed here (2 for $q \le 13$, 4 for $17 \le q \le 25$) are exactly
$\min_\tau \#E_\lambda(\mathbb{F}_q)/4$ with minima $8$ and $16$.

---

## L5 [PROVED HERE, 2026-08-12]. Sun's Lemma 5.1 = Hasse for the Legendre family; the threshold "|k| > 25" is void

**Theorem (L5).** Let $q$ be an odd prime power, $\chi$ the quadratic character of
$\mathbb{F}_q$ ($\chi(0) = 0$), $\tau \in \mathbb{F}_q \setminus \{0, \pm 1\}$, and put
$\lambda = \tau^{-2}$, $E_\lambda : y^2 = x(x-1)(x-\lambda)$ (Legendre curve; elliptic since
$\lambda \ne 0, 1$). Then
$$N(\tau) := \#\bigl\{a \in \mathbb{F}_q : \chi(1+4a^2) = -1,\
\chi\bigl(-(1+4a^2)(1-\tau^2(1+4a^2))\bigr) = 1\bigr\} \;=\; \frac{\#E_\lambda(\mathbb{F}_q)}{4}.$$
By Hasse, $\#E_\lambda(\mathbb{F}_q) \ge (\sqrt q - 1)^2 > 0$, so $N(\tau) \ge (\sqrt q - 1)^2/4
> 0$ for **every** odd prime power $q$: the hypothesis $|k| > 25$ in Sun's Lemma 5.1 (arXiv
2607.28606) can be deleted outright. The only residual obstruction in his §5 selection is
$\tau$-admissibility itself: $\mathbb{F}_3 = \{0, \pm 1\}$ carries no $\tau$ (L4); char 2 is
his §6 ($|\kappa| \ge 4$).

**Proof.** Write $u_a = 1 + 4a^2$, $w_a = -u_a(1 - \tau^2 u_a)$, and
$z = \#\{a : u_a = 0\} = 1 + \chi(-1) \in \{0, 2\}$ (the equation $a^2 = -1/4$).

*(i) Indicator.* For each $a$, $\tfrac14(1-\chi(u_a))(1+\chi(w_a))$ is $1$ exactly on the set
counted by $N(\tau)$, except: if $u_a = 0$ then $w_a = 0$ and the product is $\tfrac14$; and
if $w_a = 0 \ne u_a$ then $u_a = \tau^{-2}$ is a nonzero square, so $1 - \chi(u_a) = 0$ and
the product is already $0$. Hence
$$4N = q - T_1 + T_2 - T_3 - z, \qquad T_1 = \sum_a \chi(u_a),\quad T_2 = \sum_a \chi(w_a),
\quad T_3 = \sum_a \chi(u_a w_a).$$

*(ii) Hyperbola lemma.* For $c \ne 0$: $\sum_b \chi(b^2 + c) = -1$. Indeed
$\#\{(b,y) : y^2 - b^2 = c\} = q - 1$ (put $t = y - b \in \mathbb{F}_q^\times$, then
$y + b = c/t$: a bijection with $\mathbb{F}_q^\times$), while the same count equals
$\sum_b (1 + \chi(b^2 + c)) = q + \sum_b \chi(b^2+c)$.

*(iii) $T_1 = -1$:* $\chi(4) = 1$ and $b = 2a$ is a bijection, so
$T_1 = \sum_b \chi(b^2 + 1) = -1$ by (ii).

*(iv) Fibre count.* $\#\{a : 1 + 4a^2 = u\} = 1 + \chi(u-1)$ for every $u$: for $u = 1$ only
$a = 0$, matching $1 + \chi(0) = 1$; otherwise $a^2 = (u-1)/4$ has $1 + \chi(u-1)$ roots.

*(v) $T_2 = -1 + \sum_u \chi(u(u-1)(u-\lambda))$.* By (iv),
$T_2 = \sum_u (1 + \chi(u-1))\,\chi(\tau^2 u^2 - u)$. First piece: completing the square,
$\tau^2 u^2 - u = b^2 - 1/(4\tau^2)$ under the affine bijection $b = \tau u - 1/(2\tau)$, so
it equals $-1$ by (ii) with $c = -1/(4\tau^2) \ne 0$. Second piece:
$\chi((u-1)(\tau^2u^2-u)) = \chi(\tau^2)\,\chi\bigl(u(u-1)(u-\tau^{-2})\bigr) =
\chi(u(u-1)(u-\lambda))$.

*(vi) $T_3 = -1 - z$.* Terms with $u_a = 0$ vanish ($\chi(0) = 0$). For $u_a \ne 0$,
$\chi(u_a w_a) = \chi(-u_a^2(1 - \tau^2 u_a)) = \chi(\tau^2 u_a - 1)$. Summing over all $a$
instead adds $z$ terms $\chi(\tau^2 \cdot 0 - 1) = \chi(-1)$, and $z\chi(-1) = z$; so
$T_3 = \sum_a \chi(4\tau^2 a^2 + \tau^2 - 1) - z = -1 - z$, by (ii) with the bijection
$b = 2\tau a$ and $c = \tau^2 - 1 \ne 0$ — this is where $\tau \ne \pm 1$ enters.

*(vii) Assembly.* $4N = q - (-1) + (-1 + \sum_u \chi(u(u-1)(u-\lambda))) - (-1 - z) - z
= q + 1 + \sum_u \chi(u(u-1)(u-\lambda))$. The affine count of $E_\lambda$ is
$q + \sum_u \chi(u(u-1)(u-\lambda))$; adding the point at infinity,
$\#E_\lambda(\mathbb{F}_q) = 4N(\tau)$. $\blacksquare$

**Machine verification** (h10q.py `_lemma51_identity`): the identity, Hasse, $4 \mid \#E$,
and each step evaluation (iii)-(vi), exhaustively over all odd prime powers $q \le 250$:
6494 admissible $(q, \tau)$ pairs, 0 failures (in-suite $q \le 49$; `--extended` $q \le 250$;
non-prime $q$ via $\mathbb{F}_{p^k}$ models with machine-found irreducible moduli).

**Corollaries.**
1. **L4's table explained**: $\min_\tau \#E_\lambda(\mathbb{F}_q) = 8$ ($q \in \{5..13\}$) and
   $16$ ($q \in \{17..25\}$) yield min counts $2$ and $4$; $4 \mid \#E_\lambda$ (full rational
   2-torsion, $(0,0), (1,0), (\lambda,0)$) is why the count is divisible.
2. **The Lemma-5.1 component of $E_{\mathrm{exc}}$ collapses over every global field**: the
   blanket small-residue-field exclusion $\{w : |\kappa_w| \le 25\}$ that Lemma 5.1's Weil
   threshold forced reduces to $\{w : |\kappa_w| = 3\}$ ($\tau$-vacuity; char 2 needs
   $|\kappa_w| \ge 4$ per his §6). **Scope**: his §5 puts *more* into $E_{\mathrm{exc}}$ —
   the support of $\pi$, poles of the $\Lambda$-parameters, and zeros of their pairwise
   differences; those exclusions are data-dependent and are NOT touched by L5. What L5
   proves is that no place is excluded *for selection-failure reasons* beyond
   $|\kappa_w| = 3$.
3. **Density**: $N(\tau) = q/4 + O(\sqrt q)$ with the error the Frobenius trace of
   $E_\lambda$ — the fluctuation across the family is Sato-Tate; a valid $a$ is not merely
   existent but has density $\to 1/4$.

**Remark.** Sun's route to 7 unknowns replaces Poonen's elliptic-curve machinery by
quaternion algebras; L5 exhibits the Hasse bound for the Legendre family at the heart of the
selection lemma. The elliptic curve is not avoided — it is relocated.

---

## A3 [REFRAMING]. The quantifier race IS a fibre-dimension problem (DDF, arXiv 2102.06941)

**Theorem (DDF 4.16 + 5.10 + Cor 4.12).** For $K$ finitely generated over a perfect field (so
$K = \mathbb{Q}$) and any Diophantine $D \subseteq K^n$ that is not quantifier-free definable:
$$\mathrm{rk}^\exists_K(D) = \mathrm{efd}_K(D) + 1,$$
where $\mathrm{efd}$ is the essential fibre dimension. **(Cor 6.17)** $\mathrm{efd}_K(D) \le d$
iff $D = f(V(K))$ for a $K$-morphism $f: V \to \mathbb{A}^n_K$ **all of whose fibres have
dimension $\le d$**.

**Consequences for Question 5.7.** The minimal number of unknowns for $\mathbb{Q}\setminus
\mathbb{Z}$ equals $\mathrm{efd}_\mathbb{Q}(\mathbb{Q}\setminus\mathbb{Z}) + 1$. So:
- Daans' 10 $\iff$ a parametrization with 9-dimensional fibres; Sun's claimed 7 $\iff$
  $\mathrm{efd} \le 6$; a record of 6 $\iff$ **covering $\mathbb{Q}\setminus\mathbb{Z}$ by the
  image of a morphism with 5-dimensional fibres**.
- The natural fibration behind Sun's shape has fibre dimension exactly 6 over a generic point
  ($(a,b)$: 2, $\Phi$-witnesses: 2 after fusion with the block's quadric surface: 2 more) — 7 is
  the fixed point of this geometry, not of the bookkeeping. Beating it means removing a fibre
  dimension, e.g. rigidifying $(a,b)$ as a function of $z$ (channel 3 in geometric form).
- Lower bounds: $m \ge 2 \iff \mathrm{efd} \ge 1$, which is exactly DDF §8's thin-set argument.
  Proving $m \ge 3$ means proving no fibre-dimension-1 cover exists — no known technology.
**The record race and the lower-bound wall are the same question: pin
$\mathrm{efd}_\mathbb{Q}(\mathbb{Q}\setminus\mathbb{Z}) \in [1, 6]$.**

---

## L6 [NEW ARCHITECTURE, 2026-08-13]. The witness-tie $a = 1 + 2s$: a 6-unknown shape for $\bigcup_{w \notin S}\mathfrak{m}_w$ over $\mathbb{Q}$

**Setting.** $K=\mathbb Q$, $S=\{2\}$, $\pi=2$, $u=1$. Apply Sun's bridge
(= Daans Lemma 5.5) to $z^3$ (a polynomial pullback, so the witness count is
unchanged), with $A=1+4a^2$, $B=2b$, and
$c=h(a,b,z^3)=a^2z^6g(a,b)/(1-z^3-a^2z^6)$,
$\delta_\tau = 1 - A\tau^2$, and the block
$$\Psi_\tau(a,b,c):\ \exists y,r,s\ [\delta_\tau(c^2 - Ay^2) - 16B(r^2 - As^2) = 16].$$

**The tie and the two branches.** Demand that the block witness $s$ *be* the
$\Phi$-parameter: $a:=1+2s$. Define the 2-witness Diophantine block
$$
\begin{aligned}
\Theta_*(a,b,c,s):\quad \exists y,r\quad&
\bigl[c^2-Ay^2-16Br^2=16-16ABs^2\bigr]\\
&{}\vee
\bigl[c^2-Ay^2-16ABr^2=16A-16A^2Bs^2\bigr].
\end{aligned}
$$
The first line is tied $\Psi_0$; the second is tied
$\Psi_{2a/A}$ multiplied by $A$. The candidate definition, over base tuple
$(b,s,z)$, is
$$F(z):=\exists b,s\ \bigl[(1+2s,b)\in\Phi_1^{\{2\}}\ \wedge\
\Theta_*(1+2s,b,h(1+2s,b,z^3),s)\bigr].$$
Here $F(z)\Rightarrow z\in\bigcup_{w\notin S}\mathfrak m_w$ is proved below;
the reverse implication is exactly the open L6-assembly lemma. All displayed
rational expressions are cleared as in the next paragraph.

**No hidden denominator or inequation variable.** On $\Phi_1^{\{2\}}$ we have
$a\equiv1\pmod2$, hence $a$ is a $2$-adic unit and
$A=1+4a^2\equiv5\pmod8$, so $A$ is not a square in $\mathbb Q_2$ and in
particular is nonzero. If $1-z^3-a^2z^6=0$ for $z\in\mathbb Q$, then
$u=z^3$ solves $1-u-a^2u^2=0$, whose discriminant is $A$; this contradicts
$A\notin\mathbb Q^2$. Thus the denominator of $c=h(a,b,z^3)$ clears without
an inverse witness. The $\tau=0$ branch has
$\delta=1$. In the $\tau=2a/A$ branch, $\delta=1/A$ and multiplying the tied
equation by the already nonzero $A$ gives the second polynomial equation in
$\Theta_*$. No $\delta\ne0$ condition or inverse witness is added.


**Count [AUDITED].** Base unknowns $b, s$: 2. $\Phi$-membership of $(1+2s, b)$: at most
$\exists_3$ (Daans Lemma 5.2 pulled back along the polynomial map
$(b,s) \mapsto (1+2s,b)$, which adds no witness). The tied block is $\exists_2$
($y,r$; $s$ is a base coordinate). DDF Theorem 1.4 applies to the two arbitrary
Diophantine subsets of the common $(b,s,z)$-space and gives $3+2-1=4$ witnesses
for their intersection, hence $2+4=\mathbf{6}$ after projecting away $(b,s)$.
If either set has existential rank zero (the only unmet hypothesis of DDF 1.4),
the trivial conjunction bound is at most 3, hence no worse. The two branches in
$\Theta_*$ preserve $\exists_2$ (multiply their two defining polynomials);
denominator clearing and the finite exceptional-place patch add
no unknowns. An independent adversarial audit agreed with this count;
the DDF v5 Theorem 1.4 statement was then rechecked directly. In efd form (A3),
the tie rigidifies one $(a,b)$-direction and drops the generic fibre bound from 6 to 5.


**Soundness — inherited, zero new obligations.** A tied solution $(b,s,y,r)$ is
a $\Psi_0$- or $\Psi_{2a/A}$-solution $(y,r,s)$ for $(a,b)=(1+2s,b)$; and
$(1+2s,b)\in\Phi_1^{\{2\}}$ is
part of the formula, whose $\exists_3$-certificate already witnesses $s \in \mathcal{O}_2$
(the $a$-side congruence $a \equiv 1 \bmod \mathfrak{m}_2$ *is* $v_2(s) \ge 0$) and
$b \in \mathcal{O}_2^\times$. So Sun's soundness chain (identities 2.2–2.4, Lemma 2.1) applies
verbatim: $c$ integral at every finite ramified place of $(A, B)$, and the bridge gives
$z \in \bigcup\mathfrak{m}_w$.

**W0 [PROVED HERE — a canonical two-branch selector, including $\mathbb{F}_3$].**
Let $k$ be any finite field of odd order and let $\chi$ be its quadratic character.
For every $d\in k^\times$, counting the solutions of
$u^2-v^2=d$ in two ways gives
$q-1=\sum_v(1+\chi(v^2+d))$, hence
$\sum_v\chi(v^2+d)=-1$. Taking $d=1$ and using the bijection $v=2a$ yields
$$\sum_{a\in k}\chi(1+4a^2)=-1.$$
Consequently some $a$ makes $A=1+4a^2$ a nonsquare (and hence nonzero). For any
such $a$, set
$$
\tau_*(a)=
\begin{cases}
0,&\chi(-1)=-1,\\[2pt]
2a/A,&\chi(-1)=1.
\end{cases}
$$
If $\tau_*=0$, then $\chi(-\delta_{\tau_*}A)=\chi(-A)=1$. In the other branch,
$\delta_{\tau_*}=1-A(2a/A)^2=1/A$, so
$-\delta_{\tau_*}A=-1$ is a square. Thus $\delta_{\tau_*}\ne0$ and Sun's target
condition $\chi(-\delta_{\tau_*}A)=1$ always holds. Since $s\mapsto a=1+2s$ is
a bijection over $k$, the witness tie costs no residue classes. Over $\mathbb{Q}$,
the branch is simply $\tau_*=0$ for $w\equiv3\pmod4$ and
$\tau_*=2a/A$ for $w\equiv1\pmod4$. This removes the earlier appeal to L5 and
also removes its $\mathbb{F}_3$ vacancy. The identity is proved above; all
61 odd primes $w<300$ are additionally checked exactly by `h10q.py::_verify_L6`.


**W1 [PROVED HERE — the $w$-adic selection is tie-free].** Let $w$ be an odd
place with $v_w(A)=v_w(\delta_\tau)=0$, $v_w(b)=1$, $v_w(z)\ge1$, and
$v_w(s)\ge0$. The bridge formula gives either $c=0$ or
$v_w(c)=2v_w(a)+6v_w(z)-2\ge6v_w(z)-2\ge4$. Then the tied conic
$$-\delta_\tau A\,y^2-16B\,r^2=16-\delta_\tau c^2-16ABs^2$$
is solvable over $\mathbb Q_w$ **iff**
$\chi_w(-\bar\delta_\tau\bar A)=1$ — the same condition as Sun's Lemma 5.1,
with $s$ absent. *Proof.* The RHS is a $w$-unit $\equiv16$
(the subtracted terms have valuations at least $8$ and $1$, respectively).
A binary form $\alpha y^2+\beta r^2$ with
$v_w(\alpha)=0$, $v_w(\beta)=1$
represents a $w$-unit $\gamma$ iff $\chi_w(\gamma/\alpha) = 1$ (unit part must be represented
mod $w$ by $\alpha y^2$, then Hensel; the $\beta$-term cannot reach unit classes with
$\chi_w(\cdot/\alpha) = -1$ at any valuation of $r$). Here $\gamma/\alpha \equiv
16/(-\delta_\tau A)$ and $\chi_w(16) = 1$. $\square$ Machine check:
`h10q.py::_verify_L6` runs 2000 deterministic exact local instances with 0 mismatches
(an independent kernel screen ran 4638, also 0 mismatches). **Consequence: W0
selects a valid $\bar a$ and the canonical branch $\tau_*$ over every odd residue
field, including $\mathbb F_3$. The tie costs nothing at the target:
$\bar s=(\bar a-1)/2$ is determined mod $w$, and the imposed $2$-adic and
$w$-adic congruences on $s$ are simultaneously satisfiable by weak
approximation. This does **not** settle W2's global norm constraint.**

**W2 [PROVED HERE — exact global obstruction].** Put
$E_\tau=\mathbb Q(\sqrt{-\delta_\tau AB})$. The tied conic is equivalent to
the norm equation
$$
N_{E_\tau/\mathbb Q}
\bigl(\delta_\tau A\,y+4r\sqrt{-\delta_\tau AB}\bigr)
=-\delta_\tau A\,M_\tau,\qquad
M_\tau=16-\delta_\tau c^2-16ABs^2.
$$
Indeed the left side is
$\delta_\tau^2A^2y^2+16\delta_\tau ABr^2
=-\delta_\tau A(-\delta_\tau Ay^2-16Br^2)$.
Thus global conic solvability is **exactly** the disjunction
$$M_\tau=0\quad\text{or}\quad
-\delta_\tau A M_\tau\in N_{E_\tau/\mathbb Q}(E_\tau^\times).$$
The zero case is genuine and not absorbable: $M_\tau=0$ is solved by
$(y,r)=(0,0)$, yet $0\notin N_{E_\tau/\mathbb Q}(E_\tau^\times)$. For a
quadratic extension, Hasse's norm theorem makes the second disjunct
equivalent to membership in every local norm group.
The two formula branches simplify to
$$
\begin{array}{c|c|c}
\tau& E_\tau&M_\tau\\ \hline
0&\mathbb Q(\sqrt{-AB})&16-c^2-16ABs^2\\
2a/A&\mathbb Q(\sqrt{-B})&16-c^2/A-16ABs^2.
\end{array}
$$
At odd places $v\nmid2w$ with $v(\delta_\tau A)=v(B)=0$ but
$v(M_\tau)>0$ odd and $\chi_v(-\delta_\tau AB)=-1$, the conic can die locally
even though Sun's quaternary form survives (the freed $s$-direction was exactly
the escape). The needed L6-assembly lemma is: *for the target residue class
from W0 and every $z$ with $v_w(z)\ge1$, one can choose $(s,b)$ with
$(1+2s,b)\in\Phi_1^{\{2\}}$, $s\equiv(\bar a_w-1)/2\bmod\mathfrak m_w$,
$v_w(b)=1$, and one of the two branches so that $M_\tau=0$ or
$-\delta_\tau A M_\tau\in N_{E_\tau/\mathbb Q}(E_\tau^\times)$.* This is a Sun-§7-style
simultaneous weak/norm-approximation problem with one extra moving part
($s$ enters $M_\tau$ quadratically). It is not proved here; the certified-dead
$s=0$ freeze (A1 channel 2) shows the risk is real. The norm identity is checked
on 200 deterministic exact rational instances by `h10q.py::_verify_L6`.

**E6 [MACHINE EVIDENCE, NOT A PROOF OF ASSEMBLY].** Under the proven-primality engine,
canonical tied certificates $(s,b,\tau_*(a))$ exist for **all 61 odd primes
$w<300$**, with $A$ a nonsquare $w$-unit, $\Delta=\{2,w\}$ exactly, $c$
$\Delta$-integral, and the tied conic globally solvable. The frozen table in
`h10q.py::_L6_WITNESSES` is rechecked end-to-end by `_verify_L6` (24/24 below
100 by default; 61/61 with `--extended`): 32 use $\tau=0$, 29 use
$\tau=2a/A$. The guarded search in `l6_search.py --canonical` enforces
$v_w(A)=0$ and $\chi_w(A)=-1$ **before** its $b$ loop, uses $|s|\le3$ with
denominators $\le13$ and $b=mw$ with odd $|m|\le15$, and tests only the W0
branch. An earlier unguarded run was discarded: at $w=5$ it accepted $A=5$,
which is not a $w$-unit and therefore did not test W0/W1. After the guard was
added, all 61 targets were searched afresh and independently reverified; no
uncertified factorization was accepted.


**Status.** Architecture/count: **audited — 6 is correct conditional bookkeeping**.
Target-place theory and exact norm obstruction: **proved** (W0–W2).
Global tied-conic assembly: **open**, precisely scoped above; E6 is bounded
evidence, not its proof. If L6-assembly holds,
the main union outside $S$ is $\exists_6$; adjoining the finitely many omitted
$\mathfrak m_v$'s (each $\exists_3$) by a finite disjunction stays $\exists_6$, so
$\mathbb{Q}\setminus\mathbb{Z}$ is $\exists_6$, $\mathbb{Z}$ is $\forall_6$, and
$\mathrm{efd}_\mathbb{Q}(\mathbb{Q}\setminus\mathbb{Z})\le5$. This is independent
of whether Sun's own §§3–8 chain survives refereeing: the tie reuses the bridge but
replaces his block-completeness step.


---

## L7 [STRUCTURAL COROLLARIES + BOUNDED PROBE, 2026-08-15]. What the open assembly problem is, exactly

Nothing here changes L6's status: the assembly lemma remains **open** and the
count 6 remains conditional. L7 pins down the *shape* of the remaining
problem. Machine checks: `h10q.py::_verify_L7`.

**L7a [tie-cost dictionary — corollary of W2 + standard local theory].**
Let $v$ be an odd place with $v(\delta_\tau A)=v(B)=0$ and $M_\tau\ne0$. Then:
(i) if $v(M_\tau)$ is even, the tied conic
$-\delta_\tau Ay^2-16Br^2=M_\tau$ is soluble over $\mathbb Q_v$;
(ii) if $v(M_\tau)$ is odd, it is soluble iff $\chi_v(-\delta_\tau AB)=1$;
(iii) if **additionally** $v(A)=0$ — so that all three coefficients of
$\langle-\delta_\tau A,-16B,16AB\rangle$ are $v$-units — the
*untied* quaternary block
$-\delta_\tau Ay^2-16Br^2+16ABs^2=16-\delta_\tau c^2$ is soluble over
$\mathbb Q_v$ regardless of the right-hand side.
*Proof.* (i)–(ii): a binary form of unit discriminant over $\mathbb Q_v$
represents exactly the classes prescribed by $\chi_v$ of its discriminant;
with even valuation the target is a unit times a square and Hensel applies;
with odd valuation solubility is equivalent to isotropy of the reduced binary
form, i.e. $\chi_v(-\delta_\tau AB)=1$. (iii) Chevalley–Warning gives a
nontrivial zero of the unit ternary over the residue field; the reduction is
nonsingular (diagonal, unit coefficients, $v$ odd), so Hensel lifts it and
the ternary is isotropic over $\mathbb Q_v$; a regular isotropic form splits
off a hyperbolic plane and is therefore universal. $\square$

**L7a$'$ [coefficient-bad places].** The unit hypothesis in (iii) is
**necessary**, and coefficient-bad places
($v(A)\ne0$ or $v(B)\ne0$) carry *real* obstructions beyond the wild set,
all machine-witnessed in-suite:
- $a{=}1,b{=}11$ ($A{=}5$, $\tau{=}2/5$, $\delta_\tau A{=}1$): the ternary
  $\langle-1,-352,1760\rangle$ is anisotropic over $\mathbb Q_5$ and the
  untied block fails at $v{=}5\mid A$ (checked: RHS for $z{=}1$ not
  represented).
- guarded $w{=}5$ shape $s{=}3$ ($a{=}7$, $A{=}197$), $b{=}15$, $z{=}15$,
  $\tau{=}2a/A$: at $v{=}3\mid B$ the **tied** conic is insoluble although
  $v_3(M_\tau)=0$ — an obstruction at a non-wild place.
- the bad set is **not** "$v\mid AB$": valuations cancel in the product.
  $a{=}1,b{=}1/5$ gives $v_5(A)=1$, $v_5(B)=-1$, $v_5(AB)=0$, yet the
  ternary $\langle-1,-32/5,32\rangle$ is anisotropic over $\mathbb Q_5$ and
  the untied RHS fails there.
**Consequence: the tie's cost, correctly supported.** At odd places with
$v(A)=v(B)=0$, Sun's free $s$-direction absorbed every odd-valuation prime
of the target and the tie $s=(a-1)/2$ re-exposes exactly the wild places of
$M_\tau$ (i)–(ii). The full candidate obstruction set of a tied conic is
$$\{2,\infty\}\ \cup\ \{v:\ v(M_\tau)\ \text{odd}\}\ \cup\ \{v:\ v(A)\ne0
\ \text{or}\ v(B)\ne0\},$$
and the third constituent is nonempty in practice (the $v{=}3$ example).
Dictionary + untied-universality verified on machine-selected unit places
every run (counts printed by the suite).

**L7b [no rational-function witness — corollary, conditional on L6 soundness].**
Fix any $(s,b)$ with $(1+2s,b)\in\Phi_1^{\{2\}}$ and either branch $\tau$.
There are **no** rational functions $Y(z),R(z)\in\mathbb Q(z)$ with
$-\delta_\tau A\,Y(z)^2-16B\,R(z)^2=M_\tau(z)$ identically.
*Proof.* Such a pair has finitely many poles, so it specializes to a rational
solution at all but finitely many $z$, in particular at some non-target
$z=\pm2^k$ (infinitely many available). A solution there, together with the
$\Phi$-certificate of $(s,b)$, would put $z\in\bigcup_{w\ \mathrm{odd}}\mathfrak m_w$
by L6 soundness (Sun identities 2.2–2.4 + Lemma 2.1, applied verbatim) —
false for $z=\pm2^k$. $\square$
Machine corroboration (bounded probe): the 8 canonical witnesses with
$w\le23$ $\times$ $z=\pm2^k$, $|k|\le3$, both signs and branches (224 cells)
— **0 solutions found**: 205 certified failures (Hasse–Minkowski), 19
budget refusals left unknown; an in-suite subset re-runs every time.
**Consequence (scope-limited).** L7b kills exactly the *uniform
rational-function* witness: no fixed $(s,b,\tau)$ admits a section
$Y,R\in\mathbb Q(z)$. It does **not** by itself rule out a fixed pair
covering every relevant rational fiber pointwise — conic bundles can have
points in all rational fibers without a section, and the two branches may
alternate with $z$. Separately, the rescue table *certifies*
(Hasse–Minkowski) that each probed canonical pair fails **both** branches at
some $z\in\mathfrak m_w$ (e.g. $w{=}5$, $z{=}-5$; 7 cells re-certified
in-suite every run, covering all five probed pairs $w=3,5,7,11,13$), so
those particular pairs do not cover $\mathfrak m_w$; whether some *other*
fixed pair could remains open.

**L7c [reciprocity parity — corollary of W2 + Hilbert reciprocity].**
For $M_\tau\ne0$, the tied conic is soluble over $\mathbb Q_v$ iff
$t_v:=\bigl(-16\delta_\tau AB,\ -\delta_\tau A M_\tau\bigr)_v=1$, and
$\prod_v t_v=1$. Hence the set of obstructed places of a failing tied conic
has **even cardinality $\ge2$**: a witness never fails at one place alone,
and any repair must flip places in pairs.
*Proof.* $\gamma\ne0$ is represented by $\langle\alpha,\beta\rangle$ over
$\mathbb Q_v$ iff $(-\alpha\beta,\alpha\gamma)_v=1$; take
$\alpha=-\delta_\tau A$, $\beta=-16B$, $\gamma=M_\tau$. The product formula
is Hilbert reciprocity. $\square$
Verified on 8 disproved cells in-suite every run; 42 more in the recorded
2026-08-15 session sweep (not persisted in-suite).

**Bounded completeness probe [evidence only, NOT a proof].** Sweep of
2026-08-15, proven-primality engine throughout:
- Coverage of the *fixed* canonical witness $(s_w,b_w)$ over
  $z=wu$, $u$ in a 22-value menu (integers, negatives, non-integral):
  $w{=}3$: 6/14, $w{=}5$: 6/16, $w{=}7$: 2/19, $w{=}11$: 2/22, $w{=}13$: 6/22
  (budget-refused cells excluded from both counts) — fixed witnesses cover a
  *minority*, consistent with L7b.
- Witness-switching over the free $\Phi$ space ($v_2(s)\ge0$, $v_2(b)=0$,
  no $w$-conditions; $|s|$ menu of 15, odd $b$ menu of 16, both branches):
  **71/71** previously uncovered cells acquired a certified witness —
  every one of the 93 tested target cells is completeness-positive, against
  0/224 non-target cells. First-hit witnesses scatter across $s$, $b$, and
  both branches (26/71 use $s=0$); no single small family suffices.
  Ten frozen spot-checks (`_L7_RESCUES`) re-verify in-suite.

**What would close the lemma.** For each $(w,z)$, $z\in\mathfrak m_w$:
exhibit $(s,b)\in\Phi_1^{\{2\}}$ and a branch with empty obstruction set —
by L7a/W1 the candidate bad set is $\{2,\infty\}\cup\{\text{wild places of }
M_\tau\}\cup\{v: v(A)\ne0\ \text{or}\ v(B)\ne0\}$ (the target $w$ is soluble
$z$-uniformly for guarded witnesses by W1; coefficient-bad places are genuine
candidates — see the $v{=}3$ example), by L7c repairs come in pairs, and
L7b bars any single rational-function section (the probe shows the probed
fixed pairs miss cells; pointwise coverage by some other fixed pair remains
open). The probe suggests the free $b$-direction always suffices; turning
that into a proof (Dirichlet/weak-approximation steering of
$E_\tau=\mathbb Q(\sqrt{-\delta_\tau AB})$ against the moving wild set of
$M_\tau$ *and* the moving coefficient-bad set $\{v: v(A)\ne0\ \text{or}\
v(B)\ne0\}$) is precisely the
open problem.
L8 below rules out the exact square-class matching mechanisms for this
steering and sharpens the candidate bad set to odd-valuation places.


---

## L8 [THE ALIGNMENT WALL, 2026-08-15]. Exact-matching mechanisms are empty; the rest is evidence

Nothing here changes L6's status: the assembly lemma remains **open** and the
count 6 remains conditional. L8a-c are *unconditional* (proofs below);
L8d and the closing direction are **evidence and conjecture, clearly
scoped**. Machine checks: `h10q.py::_verify_L8`.

**Setting (exact value identity).** For admissible $(s,b)$ ($v_2(s)\ge0$,
$v_2(b)=0$), write $a=1+2s$, $A=1+4a^2$, $B=2b$,
$N_g:=16a^4b^2-A(b-1)^4=Ab^2\,g(a,b)$, $D_z:=1-z^3-a^2z^6$,
$$X:=a^2z^6N_g,\qquad Y:=4Ab^2D_z,\qquad P:=1-ABs^2 .$$
Then with $u_0:=c^2-16P$ and $u_1:=c^2-16AP$ one has the polynomial
identities
$$(Ab^2D_z)^2\,u_0=X^2-P\,Y^2,\qquad (Ab^2D_z)^2\,u_1=X^2-AP\,Y^2,$$
and W2 becomes: branch $\tau{=}0$ is soluble iff $M_0=0$ or
$A\,u_0\in N_{E_0/\mathbb Q}(E_0^\times)$ with $E_0=\mathbb Q(\sqrt{-2Ab})$;
branch $\tau{=}2a/A$ is soluble iff $M_1=0$ or
$A\,u_1\in N_{E_1/\mathbb Q}(E_1^\times)$ with $E_1=\mathbb Q(\sqrt{-2b})$.
So each branch value $u_\tau$ is a norm from its *source* field
$\mathbb Q(\sqrt{D_\tau})$, $D_0=P$, $D_1=AP$, while solubility asks for
$A\,u_\tau\in N_{E_\tau/\mathbb Q}(E_\tau^\times)$ against the *target*
field $E_\tau$ — and in **both** branches
$$D_\tau\cdot\operatorname{disc}(E_\tau)\equiv -2AbP =: \varkappa
\pmod{\mathbb Q^{\times2}}.$$

**L8a [absorption — sharpens L7a$'$].** Let $v$ be odd with $M_\tau\ne0$
(the case $M_\tau=0$ being trivial) and $v(\delta_\tau A)$, $v(B)$,
$v(M_\tau)$ **all even**. Then the tied conic is
soluble over $\mathbb Q_v$. *Proof.* Scale $y,r$ and divide by $v^{v(M)}$:
the form becomes $\langle u_1,u_2\rangle$ with unit coefficients
representing a unit. A nondegenerate binary form over $\mathbb F_v$
represents every nonzero class (isotropic: obvious; anisotropic: it is the
norm form of $\mathbb F_{v^2}$, and the norm is surjective); Hensel lifts.
$\square$ Consequence: the obstruction support is contained in
$\{2,\infty\}\cup\{v\ \text{odd}:\ v(A),\ v(B),\ \text{or}\ v(M_\tau)\
\text{odd}\}$ — strictly sharper than L7a$'$'s
$\{v:v(A)\ne0\ \text{or}\ v(B)\ne0\}\cup\mathrm{wild}(M_\tau)$.

**L8b [2-adic parity wall].** For every admissible $(s,b)$:
$v_2(ABs^2)\ge1$, hence $v_2(P)=0$ and
$$v_2(\varkappa)=v_2(-2AbP)=1 .$$
In particular $\varkappa\notin\mathbb Q^{\times2}$; both
$\operatorname{disc}$ classes $-2b$, $-2Ab$ have $v_2=1$, so $E_\tau$ is
always a **field** (the same 2-adic fact that protects soundness), and
$E_\tau\ne\mathbb Q(\sqrt{D_\tau})$ for every admissible witness: when
$D_\tau\notin\mathbb Q^{\times2}$ the biquadratic
$\mathbb Q(\sqrt{D_\tau},\sqrt{\operatorname{disc}E_\tau})$
always contains the 2-adically ramified third field
$\mathbb Q(\sqrt\varkappa)$ (when $D_\tau\in\mathbb Q^{\times2}$,
$\mathbb Q(\sqrt{D_\tau})=\mathbb Q$ and
$\varkappa\equiv\operatorname{disc}E_\tau$, so the claim is immediate).
$\square$

**L8c [exact-matching gauges are empty — corollary of L8b].** The two
square-class identities that would trivialize the wild value-primes
*identically* are (i)
$\operatorname{disc}(E_\tau)\equiv D_\tau$ (source $=$ target field, i.e.
$\varkappa\equiv\square$; then $u_\tau$ is a global norm from $E_\tau$ and
solubility reduces to the single residual condition $A\in N(E_\tau)$) and
(ii) $\operatorname{disc}(E_\tau)\equiv\square$
($E_\tau$ split, every class a norm). Neither has an admissible point: (i)
contradicts $v_2(\varkappa)=1$ (L8b), (ii) contradicts
$v_2(\operatorname{disc}E_\tau)=1$. Independently of L8b, for $s\ne0$ (i)
is rationally parametrized by $b=-\theta^2/(2As^2(1-\theta^2))$ or
$b=1/(2As^2(1-\theta^2))$ (all $s\ne0$ solutions, via the Pell conic
$\sigma^2-(2s\rho)^2=1$), and a $v_2(\theta)$ case check shows $v_2(b)\ne0$
on both roots for *every* $\theta$; for $s=0$ the gauge reads
$-2Ab\in\mathbb Q^{\times2}$, which forces $v_2(b)$ odd — verifier L8c
sweeps both roots.
$\square$ For a general fixed class $j$ with
$\operatorname{disc}(E_\tau)\equiv j\,D_\tau$ along a family, wild-prime
killing needs $\chi_v(j)=1$ at every odd-multiplicity value-prime $v$; a
single value-prime inert in $\mathbb Q(\sqrt j)$ defeats it. Whether such
primes must occur for every $j\notin\{\square,D_\tau\square\}$ and every
admissible family is **not proved here** (the value-primes of a thin family
need not equidistribute), so L8c rules out exactly the two identities above
and no more.

**L8d [alignment evidence — scoped, not a proof].** The frozen rescues do
not exploit any support-controlled mechanism: all 10 have wild
odd-multiplicity primes between $1.6\times10^5$ and $4.5\times10^{21}$
whose symbols are all $+1$ (verifier L8d) — the ten frozen witnesses
succeed by *alignment*. In the recorded 2026-08-15 session sweep
(grid in NOTES; not persisted in-suite), 13 of 161 admissible pairs
succeeded ($\approx8\%$) and every observed obstruction set had even size
(L7c live). Consistent with,
but not proving, the expectation that no uniform square-class mechanism
exists beyond the two ruled out in L8c.

**Direction (conjectural).** A natural route to conditional assembly: a
Schinzel-type Hasse principle (Colliot-Thélène–Sansuc style) for the slice
conic bundle $X_{z,s^*}\to\mathbb P^1_b$ with fibers
$\langle A,\,2b,\,-u_0(b)\rangle$ ($\tau=0$) and
$\langle A,\,2Ab,\,-u_1(b)\rangle$ ($\tau=2a/A$) — the form
$\langle1,\,\delta_\tau AB,\,-Au_\tau\rangle$ scaled by $A$ and reduced
mod squares. Its vertical Brauer
data begins with the $b{=}0$ residue $A$ and the residue of $-2b$ in the
octic field $\kappa(\pi)$ of the degree-8 irreducible value divisor —
computed **non-square at sampled parameters** (finite-field splitting
witnesses at $(s,Z)\in\{(1,3),(1,5),(2,3),(-1,7),(1,\tfrac13)\}$, primes
$101\le p\le137$; NOTES), so the bundle is genuinely ramified there and any
CTS-style argument must engage that residue. Executing this (local
solubility of the slice at all places, vertical Brauer analysis, Schinzel
input) is the concrete next attack; nothing about it is claimed proved.

## L9 [RECIPROCITY STEERING, 2026-08-16]. One free wild prime is never an obstruction

**L9 [steering lemma — PROVED].** Let $x,d\in\mathbb Q^\times$ and
$T=\{2,\infty\}\cup\{\text{odd }q:\ v_q(x)\ne0\text{ or }v_q(d)\ne0\}$.
Then $\prod_{v\in T}(x,d)_v=1$, so if $(x,d)_v=1$ for every $v\in T$
except possibly a single place $q_0$, then $(x,d)_{q_0}=1$ as well.
*Proof.* Hilbert reciprocity gives $\prod_v(x,d)_v=1$ over all places; at
an odd place $v\notin T$ both arguments are $v$-adic units, so
$(x,d)_v=1$; the product collapses to $T$. $\square$

**Consequence for W2 (exact).** For an admissible tied pair the conic is
soluble iff $M_\tau=0$ or $(x,d)_v=1$ at every place, with
$x=\alpha M_\tau$, $d=\alpha B$, $\alpha=-\delta_\tau A$ (W2). By L9 it
suffices to verify the symbols on $T\setminus\{q_0\}$ for any chosen
$q_0\in T$ with $v_{q_0}(d)=0$: the huge wild prime of the norm value
*rides free*. In particular, restricting the witness search to candidates
whose norm value factors as
$(\text{trial-smooth part})\times(\text{single proven prime power})$
makes $T$ fully computable at trial-division cost plus a bounded number
of proven-primality checks (at most one per surviving cofactor of $x$ and
of $d$; 251 of the 345 frozen rows need two) —
no hard factoring, hence no budget refusals — and eliminates the
"alignment luck" at the wild prime that L8d documented: with one wild
prime and aligned controlled places, its symbol is *forced* $+1$. This
retro-explains part of L8d: any rescue whose norm value has exactly one
wild odd-multiplicity prime coprime to $d$ *had* to show symbol $+1$
there; only rescues with $\ge2$ wild primes carry alignment content
beyond reciprocity at those primes.

**L9-E [steered coverage — bounded evidence, engine `l9_steer.py`,
authority `_L9_STEERED`, verifier L9-T].** 2026-08-16 runs. Grid: the
fixed 15-value pool $u\in\pm\{1,2,3,5\}\cup
\{\pm\tfrac13,\tfrac23,-\tfrac15,\tfrac73,\tfrac17,-\tfrac27\}$ filtered
per target by $v_w(u)\ge0$ (equivalently $v_w(wu)\ge1$, i.e. $z$ really
lies in $\mathfrak m_w$) — 11 values retained at $w=3$, 14 at $w=5$, 13
at $w=7$, 15 at each of the other 21 targets, so **353 cells** in all
(note $u=\pm w$ is retained and has $v_w(u)=1$, so the pool is not a pool
of $w$-units):
**345/353 cells (97.7%) received an admissible witness whose tied conic
is PROVED soluble** — 344 by the steered mechanism (exactly one wild
odd-multiplicity prime, symbol concluded from $T\setminus\{q_0\}$ by L9,
then cross-replayed; largest forced prime has 45 digits), 1 fully smooth.
Every certificate resolves $d$ completely and $x$ into proven prime
powers times even-exponent blobs coprime to $\operatorname{supp}d$
(such blob primes have $v_q(x)$ even, $v_q(d)=0$, hence symbol $+1$
without identification).  51 hits are single-path (the independent
`_l7_tied_status` replay hit its factoring budget); the other 294 carry
both proofs.  The full table is re-verified in-suite on every run
(`_verify_L9`, ~4 s) and `data/l9_steered.jsonl` must match its
serialization byte-for-byte.  The 8 uncovered cells
$(29,-2),(31,-2),(41,\tfrac17),(53,\tfrac13),(61,2),(61,-\tfrac27),
(67,2),(89,-2)$ are search-exhaustions, not certified failures: each had
$1{,}600$–$4{,}100$ *certified insoluble* candidates and no aligned one
within deep budgets ($\le120{,}000$ candidates, $420$ s, widened pools).
Their bad-place profiles are spread over $\{2,3,5,7,\dots\}$ with no
single dominating wall (NOTES) — consistent with alignment depth, not a
local obstruction of the slice.  Coverage compares to the 2026-08-15
fixed-pair sweep (13/161 pairs): searching *witnesses per cell* under the
steering filter is the difference.  Bounded evidence ONLY; the assembly
lemma stays OPEN.

**L9a [prime-b symbol — PROVED, corrected hypotheses].** Fix a cell
$(w,z)$, an admissible $s$ (so $a=1+2s$ is a $2$-adic odd unit and
$A=1+4a^2\equiv5\bmod8$, whence $A\notin\mathbb Q^{\times2}$), a branch
$\tau$, and $b=\varepsilon q_1$ with $\varepsilon=\pm1$ and $q_1$ an odd
prime satisfying the unit conditions **factor by factor**:
$$v_{q_1}(\delta_\tau)=v_{q_1}(A)=v_{q_1}(a)=v_{q_1}(z)=v_{q_1}(D_z)=0 .$$
Then $v_{q_1}(d)=1$ and, clearing denominators in $x=\alpha M_\tau$, the
numerator at $b\equiv0\bmod q_1$ is
$-\alpha\delta_\tau a^4z^{12}N_g(0)^2$ with $N_g(0)=-A$ a $q_1$-unit while
the cleared denominator $A^2b^4D_z^2$ contributes $v_{q_1}=4$: so
$v_{q_1}(x)=-4$ is even and the unit part of $x$ at $q_1$ is
$-\alpha\delta_\tau a^4z^{12}/D_z^2\equiv A\bmod\mathbb Q_{q_1}^{\times2}$
(using $-\alpha\delta_\tau=\delta_\tau^2A$). Hence
$(x,d)_{q_1}=\left(\tfrac{A}{q_1}\right)$ (the Legendre symbol of
$\mathrm{num}(A)\,\mathrm{den}(A)$), which by quadratic reciprocity is
constant on classes of $q_1$ modulo $4\,\mathrm{num}(A)\,\mathrm{den}(A)$,
and $+1$-classes exist because $A$ is not a square. $\square$

*The factor-by-factor form is necessary, not cosmetic.* On the branch
$\tau=2a/A$ one has $\delta_\tau=1/A$, so $2\delta_\tau A=2$ and a
"$q_1\nmid2\delta_\tau A$" hypothesis does **not** force $q_1\nmid A$ —
an earlier draft of this statement was false for that reason. Recorded
counterexample to the earlier form (now an in-suite boundary regression,
`_verify_L9` (b)): $w=3$, $z=6$, $s=0$, $a=1$, $A=5$, $\tau=2/5$,
$\delta=1/5$, $\varepsilon=+1$, $q_1=b=5$ — every condition of the old
phrasing holds, yet $N_g(0)=-A$ is *not* a $q_1$-unit,
$v_5(x)=-5$ (odd, not $-4$) and $(x,d)_5=-1$.

**Reduction direction (partly conjectural, exactly scoped; modulus recipe
corrected 2026-08-16).** Restrict the witness family to
$b=\varepsilon q_1$, $q_1$ prime in a fixed class $C$ modulo $N$, with
base point $b_0=\varepsilon q_1^{(0)}$, $\varepsilon$ fixed, and $q_1$
large enough that $\operatorname{sign}x$ — hence the $\infty$-symbol — is
constant. **The frozen sets are computed in L10-0 below**, and they are
much smaller than the set an earlier draft of this paragraph demanded.
Write $S_0=\{2\}\cup\operatorname{supp}(\alpha)$ for the *fixed* finite
places, $\alpha=-\delta_\tau A$; the places whose symbols must all be
$+1$ are $\Sigma=\{\infty\}\cup S_0\cup\{q_1\}$. Require $N$ divisible by
$8$, by $4\,\mathrm{num}(A)\,\mathrm{den}(A)$, and by $p^{k_p}$ for
$p\in S_0$ only, with $k_p$ the explicit Taylor exponent of L10-0, so
that $b\mapsto\bigl(v_p(x),[x],v_p(d),[d]\bigr)$ is constant on
$b\equiv b_0\bmod p^{k_p}$; require $\gcd(q_1^{(0)},N)=1$, without which
the progression contains no second prime at all. The moving place $q_1$
is deliberately **not** frozen into $N$ (that would contradict the gcd
condition): its symbol is $(A\mid q_1)$ by L9a and is fixed by the
residue class modulo $4\,\mathrm{num}(A)\,\mathrm{den}(A)$.
*Two superseded drafts, both recorded as in-suite regressions.* The first
demanded only $p^{v_p(x(b_0))+1}$ over data primes — refuted below. The
second demanded $\{p\le H\}$ together with the whole fixed-data support,
on the worry that data primes can dwarf any trial bound $H$ (cell
$(73,5)$ has $D_z=-59\cdot40077920921911$). That worry is **void**:
L10-0 shows such primes are not controlled places at all, only candidate
wild ones. Note $S_0$ can still contain large primes — through
$\operatorname{supp}(\alpha)$, e.g. $\alpha=-A$ on the $\tau=0$ branch.

*Why the explicit recipe had to be replaced.* An earlier draft demanded
only $p^{v_p(x(b_0))+1}$ for primes $p$ of $2\delta_\tau A$ and of the
cell data. That is false: a controlled symbol can flip inside the
advertised class at a prime that divides neither the data nor $x(b_0)$.
Recorded counterexample (now an in-suite regression, `_verify_L9` (b2)):
cell $(w,u)=(73,5)$ with $s=0$, $b_0=-29$, $\tau=2/5$ — the old recipe
gives least modulus $N=689120=2^5\cdot5\cdot59\cdot73$; the prime
$q=244637629=29+355N$ lies in the same class with the same sign and the
same $(A|q)=+1$, yet at $p=59$ one has $v_{59}(x(b_0))=0$ with symbol
$+1$ while $v_{59}(x(-q))=1$ with symbol $-1$.

On the corrected family the support of $d=\alpha B$ is the fixed
controlled set plus exactly $\{q_1\}$ — no uncontrolled primes enter $d$
— every controlled symbol is frozen by the class, and
$(x,d)_{q_1}=(A|q_1)=+1$ by L9a and the choice of $C$. The only
remaining places are the wild odd-multiplicity primes of $x$'s cofactor.
Assembly on the cell therefore reduces to: (i) existence of one aligned
triple $(s,\varepsilon,C)$ — a *finite* condition, **settled by L10
below** on the listed cells: reachable on the canonical branch
$\tau=2a/A$, and provably unreachable on $\tau=0$ for prime $A$;
(ii) a Schinzel-type input: infinitely
many primes $q_1\in C$ for which the cofactor of $x(\varepsilon q_1)$
beyond its (class-frozen) $S_0$-part is a *single* proven prime — a
two-condition Schinzel-H/Bunyakovsky hypothesis on the pair
$\{t,\ F_{\text{cell}}(t)\}$ — so that L9 forces its symbol. **Status
after L10:** (i) is proved on the 164 listed $w\equiv1\bmod4$ cells and
proved *impossible* on the $\tau=0$ branch for prime $A$ and admissible
$a$; (ii) is not
claimed proved anywhere. Each engine hit is nonetheless an
*unconditional* pointwise certificate for its own cell.

## L10 [CLASS-SIDE DICHOTOMY, 2026-08-16]. The prime-$b$ family aligns on one branch and is walled on the other

Throughout $x=\alpha P(b)/(D_z^2A^2b^4)$ and $d=2\alpha b$, with
$\alpha=-\delta_\tau A$, $Z=z^3$, $N_g=16a^4b^2-A(b-1)^4$ and
$$P(b)=16D_z^2A^2b^4-\delta_\tau a^4Z^4N_g(b)^2-32A^3s^2D_z^2b^5 .$$

**L10-0 [the frozen set is small — PROVED].** Two sets, deliberately
distinct:
$$S_0=\{2\}\cup\operatorname{supp}(\alpha)\ \ \text{(fixed; frozen into
$N$ by Taylor exponents)},\qquad
\Sigma=\{\infty\}\cup S_0\cup\{q_1\}\ \ \text{(all symbols must be $+1$)} .$$
The moving place $q_1$ is **never** put into $N$ — that would contradict
$\gcd(q_1,N)=1$ and empty the progression; it is controlled by L9a
($(x,d)_{q_1}=(A\mid q_1)$) together with the residue class of $q_1$
modulo $4\,\mathrm{num}(A)\,\mathrm{den}(A)$, and $\infty$ is controlled
by $q_1>Q_0$. Now $\operatorname{supp}(d)=S_0\cup\{q_1\}$, so at every
odd $p\notin\operatorname{supp}(d)$ we have $v_p(d)=0$ and
$(x,d)_p=(d\mid p)^{v_p(x)}=+1$ whenever $v_p(x)$ is even. Hence only
$\Sigma$ matters — **not** $\{p\le H\}$, **not** the full data support.
This supersedes the earlier over-large recipe: the huge fixed primes that
worried it (e.g. $40077920921911\mid D_z$ at cell $(73,5)$) are not
controlled places at all, only candidate wild ones.
*Explicit exponents.* Writing $P(b_0+h)=\sum_j c_jh^j$ (exact Taylor,
$c_j=\sum_{i\ge j}P_i\binom{i}{j}b_0^{\,i-j}$), any
$k_p>v_p(c_0)+e_p-\min_{i\ge1}v_p(c_i)$ with $e_2=2$, $e_p=0$ ($p$ odd),
$k_2\ge3$, $k_p\ge1$, forces $v_p(P(b))=v_p(c_0)$ and
$P(b)/P(b_0)\in(\mathbb Q_p^\times)^2$ on $b\equiv b_0\bmod p^{k_p}$;
since $D_z^2A^2b^4$ is a square this freezes $v_p(x)$ and $[x]_p$, while
$k_2\ge3,k_p\ge1$ freeze $[q_1]_p$ and hence $[d]_p$, for $p\in S_0$.

**L10a [canonical branch: $\alpha=-1$ — identity PROVED, frozen-set
corollary RETRACTED].** On $\tau=2a/A$,
$\delta_\tau=1-A\cdot4a^2/A^2=(A-4a^2)/A=1/A$, so $\alpha=-\delta_\tau A=-1$ **exactly**. Hence $\operatorname{supp}(\alpha)=\varnothing$.
**The inference drawn here originally — that therefore $S_0=\{2\}$,
$\Sigma=\{2,\infty,q_1\}$, and "no place dividing $A$ among them" — is
FALSE and is retracted** (see L11e/L11f). $\alpha$ is an $A$-unit but
$\delta_\tau=1/A$ is not: $v_A(\delta_\tau)=-1$ forces odd $v_A(x)$, so
$p\mid A$ stays a *controlled* place. The correct controlled set is
$S=\{2,3,5,7\}\cup\operatorname{supp}(\alpha)\cup
\operatorname{supp}(\delta_\tau)$ (L11f), which on this branch contains
$\operatorname{supp}(A)$. What survives of L10a is only the identity
$\delta_\tau=1/A$, $\alpha=-1$.

**L10b [$\tau=0$ wall — PROVED for prime $A$, admissible $a$].** On
$\tau=0$ we have $\delta_\tau=1$, $\alpha=-A$. Let $a$ be **odd** (the
standing admissibility hypothesis: $a=1+2s$ is an odd $2$-adic unit —
without it the wall is FALSE, e.g. $a=2$ gives $A=17\equiv1\bmod8$ and
$(-2\varepsilon\mid A)=+1$; frozen as an in-suite regression) with
$A=1+4a^2$ **prime**, $A\nmid 2zD_z$, and $b=\varepsilon q_1$ with $q_1\ne A$ prime
satisfying L9a's factor-by-factor units. Then $N_g\equiv16a^4b^2
\not\equiv0\bmod A$, so $v_A(c)=-1$, $v_A(M)=-2$, $v_A(x)=-1$,
$v_A(d)=1$, and $u_x=-A^2M\equiv(Ac)^2$ is a square mod $A$. Since
$A\equiv5\bmod8$, $\epsilon(A)=(A-1)/2=2a^2$ is even, so the symbol
formula gives $(x,d)_A=(-2\varepsilon\mid A)(q_1\mid A)$; and
$A\equiv1\bmod4$ gives $(A\mid q_1)=(q_1\mid A)$. Therefore
$$(x,d)_A\cdot(A\mid q_1)=(-2\varepsilon\mid A)=-1$$
for **both** signs, because $(-1\mid A)=+1$ ($A\equiv1\bmod4$) and
$(2\mid A)=-1$ ($A\equiv5\bmod8$). Two frozen symbols are permanently
anti-correlated: **no aligned class exists on $\tau=0$ for such a
witness.** Machine-checked on 80 in-suite instances (250 extended).
*Scope, exactly.* The derivation needs $a$ odd and $A$ prime, integral.
For rational or non-squarefree $A$ the identity is re-verified in Jacobi
form on 200 in-suite instances per run (600 extended) —
**evidence, not proof** (the cross
terms $(A/p^{k}\mid p)$ and rational $A$ are outside the argument). And
"the branch admits no aligned class for *any* admissible $s$" is **not**
claimed: $s$ ranges over an infinite set, the scan is finite. What is
claimed is: no aligned class for prime $A$ with $a$ odd (proved), and
none found in the scanned pool otherwise (11 $s$-values $\times$ 2 signs
$\times$ 300 primes per cell, 189 cells).

**L10c [SUPERSEDED — 130 of the 164 rows were invalid].** The count below
is *not* correct: the certificate froze $S_0=\{2\}\cup
\operatorname{supp}(\alpha)$, which on this branch ($\delta=1/A$) omits the
places $p\mid A$, where $v_p(x)$ is odd and the symbol is systematically
$-1$. Only 34 rows keep valid symbols, and those carry stale moduli under
L11f's complete controlled set; the corrected and complete statement
is L11f/L11h, whose table replaces this one. The original text is kept below
for provenance.
For each
of the **164 canonical grid cells with $w\equiv1\bmod4$** (the $u$-pool is
single-sourced as `_L9_U_POOL`; the verifier asserts the table's key set
*equals* that cell set), `_L10_CLASSES` records
$(s,\varepsilon,q_1,N)$ with: every frozen symbol $+1$; $q_1>Q_0$, the
Cauchy bound of $P$, so the asymptotic $\infty$-sign is the actual one;
L9a's units holding factor by factor; and $\gcd(q_1,N)=1$. That last
condition is what makes the class non-degenerate — if $q_1\mid N$ the
progression contains no further prime and Dirichlet says nothing. With
it, **Dirichlet's theorem** supplies infinitely many primes in the class,
all but the finitely many dividing the fixed cell/branch data carrying
the same frozen symbols — dropping finitely many primes leaves the
infinitude intact, and the verifier *proves* the reason for each
exclusion instead of skipping it. Moduli are small: $N\in\{40,640\}$.
**Step (i) of the
L9 reduction is proved on these 164 cells**; assembly there rests on the
Schinzel condition alone.

**Consequence (exact) — SUPERSEDED, see L11g/L11h.** W0 forces the branch
from the target: $\tau=0$ iff $w\equiv3\bmod4$. The text below concluded that
the prime-$b$ reduction is established on the $w\equiv1\bmod4$ half (164 of
353 cells); **that count is wrong** — 130 of those rows were invalid, and the
correct answer is L11h's 190 cells, characterized by a numerator prime
$\equiv1\bmod4$ rather than by $w\bmod4$. On the other 189 it
is **dead for prime $A$** (L10b) and **unobserved** otherwise — the scan
covered 11 $s$-values $\times$ 2 signs $\times$ 300 primes per cell and
found no aligned class — but *no universal impossibility is claimed*,
since $s$ ranges over an infinite set and the non-prime-$A$ identity is
evidence only. On present knowledge those cells need a different witness
certificates from L9-E — any wall here is in the class-based reduction,
not in solubility.
**Corrected and generalised by L11 (2026-08-16):** the conclusion above is
right but the reason is not branch-specific. L11g shows the same pairing
walls **every** branch, so no branch choice rescues these cells; and L11h
characterizes the reachable ones exactly (numerator prime $\equiv1\bmod4$).
L10b stays true verbatim as the $\tau=0$ instance of L11g. **L10c does
not:** its certificate omitted the controlled places $p\mid A$, and
**130 of its 164 rows are invalid** — see L11e/L11f.

## L11 [BRANCH COMPLETION, 2026-08-16]. Branches are free — but the wall is branch-independent, and the reachable cells are exactly characterized

L10 left the $w\equiv3\bmod4$ half of the grid needing "a different witness
family", and this section began by asserting the deficient object was the
**branch set**. Branches are indeed free (L11a–L11c). But the premise was
wrong: L11g shows the wall is branch-independent, and L11h characterizes
exactly which cells the prime-$b$ family can serve. The intermediate claim
L11e was **false and is retracted below**, along with 130 rows of L10c.

**L11a [branch completion is free — PROVED].** Let $\tau=p/q$ be a rational
function of $(a,b,s)$ over $\mathbb Q$ **whose denominator $q$ is proved
nonvanishing on $\Phi_1^{\{2\}}$** (all branches used here have
$q\in\{1\}\cup\{2dA:d\in\mathbb Q^\times\}$, and $A=1+4a^2\ge1$, so this
holds by inspection). Then the branch $\Psi_\tau$, tied at
$s=(a-1)/2$, may be adjoined to $\Theta_*$ at **zero cost** in all three
currencies:
1. *Well-defined, no inverse witness.* On $\Phi_1^{\{2\}}$, $a$ is an odd
   $2$-adic unit, so $A=1+4a^2\equiv5\bmod8$ is not a square in
   $\mathbb Q_2$; hence $\delta_\tau=1-A\tau^2\ne0$ identically. Clearing
   the denominator of $\delta_\tau$ multiplies the equation by a manifestly
   nonzero polynomial, exactly as the existing $\tau=2a/A$ branch does.
2. *Sound.* Soundness is **$\tau$-uniform**: the Sun identities (2.2)–(2.4)
   hold for arbitrary $\tau$ — $\operatorname{Nrd}(\gamma_\tau)
   =((1+A\tau^2)^2-4A\tau^2)/\delta_\tau^2=1$ identically — and Prop 2.1
   ($c$ integral at every finite ramified place of $(A,B)$) never mentions
   $\tau$. The bridge from there to $z\in\bigcup\mathfrak m_w$ is a statement
   about $c=h(a,b,z^3)$ and $\Phi$-membership alone.
3. *No witnesses.* $P_1=0\vee\dots\vee P_k=0$ is $P_1\cdots P_k=0$: one
   polynomial in the same two witnesses $(y,r)$. The count argument of L6
   already uses this for $k=2$; it is insensitive to $k$. **Count 6 and the
   efd bound 5 are unchanged.**

**L11b [the four square classes — PROVED].** $\alpha=-\delta_\tau A$ takes
each class of $\langle-1,A\rangle\subset\mathbb Q^\times/\square$ at an
explicit rational $\tau$:

| $\tau$ | $\delta_\tau$ | $\alpha=-\delta_\tau A$ | class | $E_\tau=\mathbb Q(\sqrt{\alpha B})$ |
|---|---|---|---|---|
| $0$ | $1$ | $-A$ | $-A$ | $\mathbb Q(\sqrt{-2Ab})$ |
| $2a/A$ | $1/A$ | $-1$ | $-1$ | $\mathbb Q(\sqrt{-2b})$ |
| $1$ | $-4a^2$ | $4a^2A$ | $A$ | $\mathbb Q(\sqrt{2Ab})$ |
| $\tau^\dagger=\dfrac{1+2a^2}{1+4a^2}$ | $-\dfrac{4a^4}{A}$ | $4a^4=(2a^2)^2$ | $\mathbf 1$ | $\mathbb Q(\sqrt{2b})$ |

*Proof of the last row.* $\tau^\dagger=(A+1)/(2A)$, so
$\delta=1-A\frac{(A+1)^2}{4A^2}=\frac{4A-(A+1)^2}{4A}=-\frac{(A-1)^2}{4A}$,
and $A-1=4a^2$ gives $\delta=-4a^4/A$, $\alpha=4a^4=(2a^2)^2$. $\square$

**L11c [the square-branch family — PROVED].** $\alpha$ is a square iff
$\alpha=\lambda^2$ with $A+\lambda^2=\mu^2$ (then $\tau=\mu/A$), i.e. iff
$(\mu-\lambda)(\mu+\lambda)=A$. Every rational factorisation is admissible,
so the square branches form a one-parameter family indexed by a free
rational $d$:
$$\tau_d=\frac{A+d^2}{2dA},\qquad \lambda_d=\frac{A-d^2}{2d},\qquad
\alpha=\lambda_d^2,\qquad \tau_1=\tau^\dagger .$$
A *finite* subfamily may be adjoined; an infinite one may not (it would
need $d$ as a witness).

**L11d [the target-place *character* is unconditional on $\tau^\dagger$ —
PROVED as a character identity; NOT a local-solubility claim].** On
$\tau^\dagger$, $\alpha=(2a^2)^2$, so for every odd $w\nmid 2a$
$$\chi_w(\alpha)=\chi_w\bigl((2a^2)^2\bigr)=+1 ,$$
with no condition on $w\bmod4$ and none on $\chi_w(A)$. Machine-checked:
$122$ $(w,a)$ pairs in-suite, of which $32$ are pairs where the canonical
$-A$ branch has $\chi_w=-1$; the search-side sweep saw $1914/1914$ against
$44$–$51\%$ for the other three branches.

*Scope — this is not yet W1.* W1's local point at $w$ (Sun Lemma 5.1)
needs $-A\delta_\tau$ to be a square **unit** at $w$, i.e. it carries the
valuation hypotheses $v_w(A)=v_w(\delta_\tau)=0$ on top of the character
condition. Those fail on exactly the branch used by the table: at $a=1$,
$A=5$, $\delta=-4/5$, the cell $w=5$ has $v_5(A)=1$ and $v_5(\delta)=-1$.
So **no claim is made here that $\tau^\dagger$ satisfies W1**, and W0 is
*not* asserted to be redundant; what is proved is that the character
obstruction — the only part of W1 that produced the $w\bmod4$ split — is
identically absent on the square branch. Supplying a replacement local
lemma covering $v_w(A)\ne0$ is open work, and it is *not* needed for L11e,
which is a statement about frozen symbols at $\{2,q_1,\infty\}$ computed
from actual values.

**L11e [RETRACTED 2026-08-16 — the claim was FALSE].** This slot claimed an
aligned class at **all 353 cells** on $\tau^\dagger$, and that step (i) of
the L9 reduction was therefore complete. It is false, and it was caught by
the adversarial audit of this very layer before anything downstream was
built on it. The certificate checked the symbols at its *advertised*
controlled set $\{2\}\cup\operatorname{supp}(\alpha)$ plus $q_1$ and
$\infty$, but L10-0 frees an off-set place only where $v_p(x)$ is **even**,
and that hypothesis was never verified. It fails systematically at $p\mid
A$. Recorded regression, exact values: cell $(w,u)=(3,1)$, $s=0$, $a=1$,
$A=5$, $\tau^\dagger=3/5$, $\delta=-4/5$, $\alpha=4$, $\varepsilon=+1$,
$q_1=41$, $N=160$ gives $c_0=9311592816/6345775$,
$d_0=328$, $v_5(x_0)=-5$ and
$$(x_0,d_0)_5=-1\quad\text{while}\quad(x_0,d_0)_2=(x_0,d_0)_{41}
=(x_0,d_0)_\infty=+1 .$$
**293 of the 353 rows were bad this way. The same defect invalidates 130 of
the 164 rows of L10c**, which shares the certificate; the other 34 keep valid
symbols but carry stale moduli once the complete controlled set of L11f is
used, so no row of that table can be cited as printed.

**L11f [the controlled set, complete — PROVED].** For a prime-$b$ witness
$b=\varepsilon q$ on branch $\tau$, take
$$S=\{2,3,5,7\}\cup\operatorname{supp}(\alpha)\cup
\operatorname{supp}(\delta_\tau).$$
Every place at which $v_p(x)$ is forced odd by fixed data lies in $S$.
*Proof.* $x=\alpha P(b)/(D_z^2A^2b^4)$ and the denominator is a square, so
$[x]=[\alpha P(b)]$. Let $p$ be odd, $p\notin S$. From $\alpha=-\delta_\tau
A$ we get $v_p(A)=v_p(\alpha)-v_p(\delta_\tau)=0$, and $v_p(a)\ge0$ (else
$v_p(A)=2v_p(a)<0$). The three terms of
$P=16D_z^2A^2b^4-\delta_\tau a^4Z^4N_g^2-32A^3s^2D_z^2b^5$ have $p$-contents
$2v_p(D_z)$, $4v_p(a)+4v_p(Z)$ and $2v_p(s)+2v_p(D_z)$ — all **even** — and
the middle one is attained uniquely at its constant coefficient
$-\delta_\tau a^4Z^4A^2$, so ties cannot cancel that witness. Hence the
content of $P$ is an even power of $p$; divide it out. What remains is a
nonzero polynomial of degree $\le8$ modulo $p$, and a nonzero degree-$8$
form over $\mathbb F_p$ has at most $8$ roots, so for $p\ge11$ ($p-1>8$) it
cannot vanish at every unit residue. Therefore no prime outside $S$ carries
a permanent odd valuation, and the only remaining candidates $p\in\{3,5,7\}$
are frozen by fiat. Finally $p\nmid N$, since
$\operatorname{supp}(A)\subseteq\operatorname{supp}(\alpha)\cup
\operatorname{supp}(\delta_\tau)$. $\square$
*(Derived by the 2026-08-16 adversarial audit; $p=3$ genuinely occurs — for
$z=5$, $q\equiv41\bmod160$, $P(q)$ is divisible by 3 at every prime $q$
because $P$ vanishes on $\mathbb F_3^\times$, while $3\nmid\operatorname{cont}(P)$.)*

**L11g [the wall is branch-INDEPENDENT — PROVED].** Let $a$ be admissible
with $A=1+4a^2$ prime, let $\tau$ be **any** branch of the square-class
group $\langle-1,A\rangle$, and let $b=\varepsilon q$ with $q$ an odd prime
satisfying L9a's units. If $v_A(z)\le0$ then
$$(x,d)_A\cdot(x,d)_q=-1 ,$$
so no aligned class exists — on *any* branch.
*Proof (square branch, $a=1$, $A=5$; the others are identical).*
$[x]=[5V]$ and $[d]=[2q]$ with $V=(Z^2N_g)^2+5(10q^2D_z)^2$. If $v_5(Z)=0$
then $V\equiv(Z^2N_g)^2\not\equiv0\bmod5$ is a square unit, and if
$v_5(Z)<0$ the first summand has the strictly lower, even valuation; either
way $V\in(\mathbb Q_5^\times)^2$, so $[x]_5=[5]$ and $v_5(x)$ is odd. Since
$v_5(d)=0$, $(x,d)_5=(d\mid5)=(2q\mid5)=(2\mid5)(q\mid5)=-(q\mid5)$. L9a
gives $(x,d)_q=(5\mid q)$, and $5\equiv1\bmod4$ makes that $(q\mid5)$.
Multiplying, $(x,d)_5(x,d)_q=-(q\mid5)^2=-1$. $\square$
This **subsumes L10b**, which is the $\tau=0$ instance: there the same
pairing appears as $\prod_{p\mid A}(x,d)_p\cdot(A\mid q_1)=(-2\varepsilon\mid
A)=-1$. The common cause is admissibility: $A\equiv5\bmod8$ forces
$(2\mid A)=-1$. Machine-checked: $195/195$ instances on each of the four
branches $\tau\in\{0,\,2a/A,\,1,\,\tau^\dagger\}$, with $v_5(x)$ odd in all
$195$; and $(x,d)_5=+1$ at $60/60$ of the cells with $v_5(z)>0$.

**So branch completion does not help the class side at all.** L11's opening
premise — that the deficient object was the branch set — is refuted by
L11g. The deficient object is the **prime-$b$ family** relative to the
*witness* $a$: the wall is escaped only by making $A$ share a prime with
$z$.

**L11h [complete characterization of the reachable cells — PROVED].** For
the prime-$b$ family, on every branch, an aligned class exists at the cell
$(w,u)$, $z=wu$, **iff $z$ has a numerator prime $p\equiv1\bmod4$.**
*Necessity*, **for witnesses with $A$ prime, is L11g**: escaping the pairing
requires $v_p(z)>0$ for some $p\mid A$, and every odd prime dividing
$A=1+4a^2$ — i.e. dividing $n^2+4m^2$ for $a=m/n$ — is $\equiv1\bmod4$.
*Composite or rational $A$ necessity is no longer open evidence*: L13e
proves the same wall for every admissible $a$ (the middle term of $P$ alone
carries no $A$ and strictly dominates it at odd $p\mid A$, so $[x]_p=[p]$;
the pairing telescopes identically) — the $878{,}400$-certificate search
(spanning $244$ admissible $a$, $228$ with composite or rational $A$, **0
aligned classes found**) merely confirms what the derivation now covers.
*Sufficiency* is constructive: given $p\equiv1\bmod4$, solve $n^2\equiv-4m^2
\bmod p$ and choose odd representatives (possible since $p$ is odd), giving
an admissible $a=m/n$ with $v_p(A)>0$ and $A\equiv5\bmod8$; verified for
the first 14 such $p$, and the resulting class search succeeds at
**190/190** of the cells the criterion predicts. `_L11_CLASSES` freezes
those 190 rows, each carrying its escape prime $p$, its constructed $a$, its
branch, and its modulus. The other **163 cells are unreachable by this family on
every branch** — proved for prime $A$, evidence-only for composite or
rational $A$ (see the scope split above).

**Consequence, stated exactly.** Step (i) of the L9 reduction holds at
$190$ of $353$ cells and is **impossible** at the remaining $163$ for the
prime-$b$ family. This replaces both L10c's $164$ (of which $130$ rows were
invalid) and L11e's retracted $353$. Assembly is **OPEN**, now for two
independent reasons: the Schinzel condition (ii) at the 190, and the
absence of any admissible family at the 163. The count stays $\mathbf 6$
and stays conditional. Global steered coverage — a different and weaker
statement — improved to $347/353$ this session: a wider $s$-pool closed
$(53,\tfrac13)$, and the square-branch family closed $(29,-2)$ at
$s=-8$, $d=1/3$, $b=-29/5$, which no $(s,b)$ reached on the two original
branches. That last hit is the one place where branch completion paid.

## L12 [THE GENERALIZED WALL, 2026-08-16]. Every $b$ coprime to the cell data fails, on every branch; the one door is non-coprimality

L11g walled the *prime*-$b$ family and L11h characterized what it can reach.
The natural next move — a cofactor $m$ with $(m\mid A)=-1$, chosen to flip the
sign — **does not work**, and the reason generalizes into a much stronger
no-go.

**L12a [the generalized wall — PROVED].** Let $a$ be admissible with
$A=1+4a^2$ prime, let $\tau$ be any branch, let $v_A(z)\le0$, and let $b$ be
any admissible rational whose odd-valuation places are coprime to the cell
data $\{2,z,D_z,a,\delta_\tau,A\}$. Write $b_{\mathrm{sf}}$ for the squarefree
part of its odd part. Then
$$\prod_{p\in\{A\}\cup\operatorname{oddsupp}(b)}(x,d)_p=-1 .$$
*Proof.* Fix an odd $p$ with $v_p(b)$ odd, $p$ coprime to the data. For
numerator primes ($v_p(b)=k>0$): $P(b)\equiv P(0)
= -\delta_\tau a^4Z^4A^2\not\equiv0$, so $v_p(x)=-4k$ is even while
$v_p(d)=k$ is odd, and the unit part of $x$ is
$$u_x\equiv\frac{\alpha\cdot(-\delta_\tau a^4Z^4A^2)}{D_z^2A^2}
=\frac{\delta_\tau^2Aa^4Z^4}{D_z^2}\equiv A\pmod{\square},$$
whence $(x,d)_p=(A\mid p)=(p\mid A)$ by reciprocity ($A\equiv1\bmod4$).
For denominator primes ($v_p(b)=-k<0$) the $p$-integral congruence
$P(b)\equiv P(0)$ no longer applies; the dual dominant-term computation gives
$v_p(x)=+4k$ (the top-degree term $\delta_\tau a^4Z^4\operatorname{Ng}(b)^2$
dominates, with $\operatorname{Ng}(b)=b^4\cdot(\text{unit}\equiv1\bmod p)$),
and the same $u_x\equiv A$ holds, so $(x,d)_p=(A\mid p)$ again. Multiplying
over those places gives $(b_{\mathrm{sf}}\mid A)$.
At $A$ itself: $[x]_A=[A]$ holds on every branch ($v_A(x)$ odd with unit a
residue), so $(x,d)_A=(u_d\mid A)$ by the parity formula with
$(A,A)_A=1$ ($A\equiv1\bmod4$), where $u_d=u_\alpha\cdot2\varepsilon
b_{\mathrm{sf}}$ up to squares and $(u_\alpha\mid A)=+1$ on every named
branch ($u_\alpha\equiv-1$ resp. $(1+2a^2)^2\bmod A$). Hence
$(x,d)_A=(2\varepsilon b_{\mathrm{sf}}\mid A)=-(b_{\mathrm{sf}}\mid A)$,
since $(2\mid A)=-1$ and $(-1\mid A)=+1$. (On the two branches whose
$\delta_\tau$ carries $A$ in the denominator, $v_A(d)=0$ and the intermediate
form $(x,d)_A=(d\mid A)$ holds literally; on $\tau\in\{0,1\}$,
$v_A(d)=1$ and $(d\mid A)$ is undefined — the branch-uniform derivation above
is the correct one.) The product telescopes to $-1$. $\square$
Machine-checked on 231 instances (default) covering **prime, semiprime,
three-prime, squarefull and rational $b$** and all four branches, with
$(x,d)_p=(A\mid p)$ asserted at every odd place of $b$.

**Consequence — the $m$-escape is impossible.** Taking $b=\varepsilon mq$ with
$(m\mid A)=-1$ does flip $(x,d)_A\cdot(x,d)_q$ to $+1$, but the same
computation at each $p\mid m$ contributes $\prod_{p\mid m}(x,d)_p=(m\mid A)$,
and $-(m\mid A)\cdot(m\mid A)=-1$. **The $-1$ relocates; it never cancels.**
Squarefull parts contribute $+1$ on both sides. No choice of odd cofactor,
no number of prime factors, and no branch changes this. **L12a subsumes
L11g** (prime $b$ is the case $b_{\mathrm{sf}}=q$) and closes the entire
coprime-$b$ direction.

**L12b [the one door — CONSTRUCTIVE, 103 of 163].** The only hypothesis left
to attack is coprimality: if $b$ shares a prime $f$ with the cell data then
$P(b)\not\equiv P(0)$ at $f$, $u_x\equiv A$ fails, and the telescoping breaks.
It does break. Taking $b=\varepsilon fq$ with $f\mid z$, **103 of the 163
cells that L11h walls acquire every controlled symbol $+1$** on the
conservative controlled set — $\{2,3,5,7\}$ together with the places of
$\alpha,\delta_\tau,A,a,z,D_z,b$, i.e. every fixed datum frozen. Frozen in
`_L12_ESCAPE` and replayed each run.

**Scope, exactly — superseded by L13a.** As frozen on 2026-08-16 this was
base-point alignment only, with class constancy asserted via the exponent
lemma but unsampled (conservative moduli $\sim10^{13}$). The minimal
controlled set $S_{\min}=\{2,3,5,7,f\}$ derived the next day shrank the
moduli to median $2.4\times10^7$, and the class alignment is now **verified**
(L13a). Assembly is **OPEN**, the count stays $\mathbf 6$ and stays
conditional.

## L13 [AFTERMATH, 2026-08-17/18]. The class side is verified; the shortcut side is dead; witness coverage 353/353 — COMPLETE (L13f, 2026-08-18)

**L13a [verified class alignment of the 103 escapes — PROVED, frozen side].**
For every row of `_L12_ESCAPE` (cell $(w,z)$ walled by L11h, $b=\varepsilon fq$
on the square branch at $a=1$), the minimal controlled set is
$S_{\min}=\{2,3,5,7,f\}$: at every prime member $Q\equiv q\pmod N$ of the
exponent-lemma modulus ($N$ built from `_l10_exponent` on $S_{\min}$;
median $2.4\times10^7$), every frozen symbol, the wild symbol
$\operatorname{hilbert}(x,d,Q)=(A\mid Q)$, and the real place are $+1$ —
verified on **883 prime members** (zero violations, zero Cauchy-sign
failures, zero excluded-member hits; `l12_class.py`,
`data/l12b_class_sample.json`). The class claim **excludes** the finitely
many members $Q$ dividing $D_z\cdot\operatorname{num}(z)\cdot a\cdot\delta$
(174 primes, named per row) and does **not** cover the emergent places
($v_\ell(P(b))$ odd, $\ell\notin S_{\min}\cup\{Q\}$): 959 probed emergent
entries (827 small-prime spot entries + 132 exact entries from the 39 fully
factored members), 927 distinct places after de-duplication, of which 452
carry symbol $-1$ — and the count is even at every fully factored member —
so step (ii) is genuinely outside the class claim, and the parity
bookkeeping is exact.

**L13b [witness coverage 353/353 — PROVED per row, COMPLETE since L13f].** Sixty walled cells
without an aligned class carry fully soluble witnesses $b=\varepsilon fq$
sharing a prime with the cell data ($f\mid z$ in 57 rows; in 3 rows —
$(59,1)$, $(79,1)$, $(79,\tfrac13)$ — it is $q=w$ itself that carries the
$z$-factor), on the square branch at $a\in\{3,5,33\}$ — 55 by Hasse-Minkowski
(`_l7_tied_status` True), 5 by full steered certificates; all 60 replayed
independently by the lead (`_L13_ESCAPE2`, `data/l12_escape2.jsonl`). Four
residual cells close the same way, all tied-status True on the square branch
(`_L13_RESIDUAL`): $(67,(2,1))$ at $a=1$, $b=4757=67\cdot71$;
$(41,\tfrac17)$ at $a=25$, $b=4985$; $(61,-\tfrac27)$ at $a=25$,
$b=55571$; $(61,(2,1))$ at $a=25$, $b=-71431$. Combined with the prior
347-cell steered coverage, grid witness coverage was 352/353 with the
single remaining cell $(89,-2)$ an L11h class cell: its aligned class
existed (step (i) done), but no decidable member had been found — its
class modulus is $1.65\times10^7$ with every early member exceeding the
primality engine's refusal threshold. **Closed on 2026-08-18 (L13f):** the
cell now carries a frozen residual witness $(a,b)=(3,\,89/367)$ on the
square branch ($\tau=19/37$), tied-status True, one of 325 certified
soluble candidates found by the cofactor ladder across two boxes; see
L13f. Grid witness coverage is **353/353 — COMPLETE**.

**L13c [no shortcut exists — PROVED per row; EVIDENCE beyond].** The degree-8
polynomial $P$ is **irreducible over $\mathbb Q$ on all 706 (cell, branch)
rows** (sympy `factor_list` certified; 24 rows per default run and 96 per
extended run re-certified in-suite by Frobenius irreducibility mod
$p<1000$, `_l13_irred8`; an irreducibility prime exists for a positive-density
set of $p$ (8-cycle Frobenius, density $\approx1/8$)), and its content's squarefree part is confined to
$\{\pm1,\pm5\}$ — the constant factor never introduces an emergent place.
Consequently **no $P=cR^2$ constant-square or $P=cR^2S$ weak-form shortcut
exists anywhere on the grid**: the emergent step (ii) cannot be collapsed
structurally. The search evidence goes further: 137.5 million exact
square-class tests over four branches and $a\in\{1,3,5,9\}$ found **zero**
admissible $b$ with $\operatorname{sqf}(P(b))$ frozen-supported, and 40
designed families (progressions, odd-unit scalings, Pell) are all clean —
consistent with the genus-3 heuristics for $y^2=f\cdot P(b)$
(`data/l12_squarehunt.json`). Frobenius cycle-type sampling (28,240
certified unramified primes) matches $D_8\wr C_2$ exactly (all 10 classes,
no others; $(7,1)$ and $(3,\dots)$ never occur) — **evidence, not proof**,
that $\operatorname{Gal}(P)\cong D_8\wr C_2$; $S_8$ is excluded with
overwhelming Chebotarev evidence.

**L13d [emergent-free members exist — certified instances].** Aligned
classes have certified members with $\operatorname{ramified}(x,d)=\varnothing$
and tied status True: **24 certified members across 15 distinct cells** — 16
frozen rows in `_L13_ZERO_BAD` replayed every run (incl. L11 cells
$(5,\tfrac17)$, $(7,5)$, $(11,-5)$, $(17,-5)$, $(13,-\tfrac27)$,
$(37,-\tfrac27)$ and L12 escape cells $(3,\pm1)$, $(3,\pm2)$,
$(3,\tfrac13)$...). Sampling density: 15 zero-bad among 57 fully factored
members in the dedicated sweep (`data/l13_zerobad_more.jsonl`) plus 9 of 44
in the probe (`data/l12_density.json`) — the step-(ii) condition has
positive empirical density (~25%), nowhere zero.

**L13e [the wall extends to composite squarefree $A$ — PROVED, audited
derivation].** The L12a pairing never used $A$ prime: at every odd $p\mid A$
the middle term $-\delta_\tau a^4Z^4\operatorname{Ng}^2$ strictly dominates
$P$ ($v_p(\operatorname{Ng})=0$ always), so $[x]_p=[p]$ — $v_p(x)$ odd with
unit part in the residue class of $A$'s unit — and every prime dividing
$1+4a^2$ is $\equiv1\bmod4$ (reciprocity factor-by-factor, Jacobi cross
terms cancelling). Hence $(x,d)_p=(2b\mid p)$ at every $p\mid A$ of odd
valuation, $(x,d)_\ell=(A\mid\ell)$ at every odd place of $b$, and the
product telescopes to $(2\mid\ker A)^{1+v_2(b)}=-1$ for every admissible
$v_2(b)=0$ witness — composite or rational $A$ included (even-exponent
primes drop out; $v_2(b)=1$ flips the sign: the exact admissibility
boundary). **Consequence: L11h's necessity (reachable $\Rightarrow z$ has a
numerator prime $\equiv1\bmod4$) is a theorem for ALL admissible $a$**,
and the 878,400-certificate search merely confirms it. `h10q.py` L13e
replay: 40 instances per run across $a\in\{3,9,15,25,33,45,51,63,65,77\}$
and all four branches (120 extended); full derivation and 316-instance sweep
in `data/l13_compositewall.{md,json}` (CompositeWall agent; lead-verified).

**L13f [the cofactor ladder; the last cell falls — witness mechanism PROVED
per row, 2026-08-18].** Emergent-freeness of a member $b$ never needs
$P(b)$ fully factored: strip the primes of $S\cup\operatorname{supp}(b)$
from $P(b)$ (numerator and denominator separately; strip-subsets decided
exactly by `_l10_supp`), then a member is **zero-bad whenever the remainder
is (i) a perfect square (no emergent places at all) or (ii) a single proved
prime $R$ with $v_R$ odd** — because with frozen/wild/$\infty$ symbols all
$+1$ the parity law (even bad count from reciprocity) forces
$(x,d)_R=+1$ on the sole surviving emergent place. Both rungs are
factorization-free at engine scale: rung (i) is `isqrt`-arithmetic; rung
(ii) is kernel-proved primality (deterministic MR / Pocklington) plus one
Hilbert-symbol verification. Method validated on ground truth: all 16
frozen `_L13_ZERO_BAD` rows return ladder-zero, including the
$(5,(1,7))$-$k{=}65$ row whose 48-digit remainder is a proved prime; the
live parity-law cross-check fired **inconsistent = 0** across 5,767
decided rows. Applied to cell $(89,-2)$: over 309,392 structured
candidates ($a\in\{1,3,5,9,17,25,33\}$, three branches,
$b=\pm f^e q^{\pm1}$, $q\le 8000$; `l13_filter.py stage A–C`,
`data/l13_filter_run2.json`) **236 certified soluble candidates** (133
prime-rung + 103 exact-factorisation rung), plus 89 more in a dense
$|b|<2\times10^4$ four-branch box (DenseScan dense box, stage D; 46 + 43). Of
the first-stage subfamily obstructions: 6,034 aligned candidates share an
emergent $-1$ pair (verdict *bad*, exact), 1,896 cofactors exceed 72
digits (recorded as refusals with reason, never evidence), 516 carry an
unproved Jacobi remainder (evidence tier only). The **frozen** witness
$(a,b,z,\tau)=(3,\,89/367,-178,19/37)$ was independently replayed by the
lead through the full five-rung chain — `smooth_emergent` cofactor-big →
29-digit remainder **proved prime** (Pocklington) → $v_R=1$ →
$\operatorname{hilbert}(x_0,d_0,R)=+1$ → `_l7_tied_status` **True** — and
is suite-asserted in `_L13_RESIDUAL` (7 rows) on every run. **Consequence:
every one of the 353 grid cells carries a replayable certified soluble
assembly witness (E2 ResidualGridWitness COMPLETE); the per-candidate hit
prior of the earlier square-subfamily scan (~$10^{-11}$, 0/290,928) was
not evidence, and the ladder's yield (236/4,481 among stored aligned
cofactor-big rows ≈ 5.3%) measures the prime-cofactor mechanism, not
luck.** The remaining open obligation is untouched: hypothesis **H** (an
emergent-free member of *each aligned class* — the analytic input to the
conditional record) stays OPEN; the ladder closes *this grid*, not H,
since H requires the class-member statement for all cells uniformly rather
than one witness per cell.
*Wording note (2026-08-19).* H is existential in both the class and the
member: it needs **one selected aligned class per cell** carrying an
emergent-free member, not a member of *every* aligned class. Read "each
aligned class" in this and neighbouring rows as "the selected class of each
cell"; the grid ledger verifies all 293, which is stronger than H needs.

*Wording note (2026-08-19).* H is existential in both the class and the
member: it needs **one selected aligned class per cell** carrying an
emergent-free member, not a member of *every* aligned class. Read
"each aligned class" in this and neighbouring rows as "the selected
class of each cell"; the grid ledger verifies all 293, which is
stronger than H needs.


**Normalization note.** The kernel's `_l10_P` evaluates to
$A^2\cdot b^4D_z^2\cdot M_0$ — its internal $\operatorname{Ng}$ carries one
factor $A$ versus the prose formula; $A^2$ is a square, so every Hilbert
symbol, valuation parity, and factor-degree statement is invariant. The
identity is asserted at every L13 replay point.

**L14 [hypothesis H, instance-verified on the entire grid — 293/293 classes,
PROVED per row, 2026-08-18].** Every one of the 293 aligned classes (190
L11 + 103 L12 escapes) carries a certified **emergent-free member**, frozen
in `h10q.py::_L13_H_CLASSES` and replayed on every suite run (60 seeded in
default, all 293 in extended: ladder verdict *zero* on each; `ramified(x0,d0)`
empty on 242, 51 budget refusals on the auxiliary cross-check — logged,
never evidence). Method chain of the day: (i) the L13f **cofactor ladder**
turns emergent-freeness into a factorization-free decision whenever the
stripped remainder of $P(b)$ is a square or a *single proved prime* (parity
law forces its symbol $+1$); (ii) the kernel prover gained the
**Brillhart–Lehmer–Selfridge relaxation** ($F\ge n^{1/3}$, exact
two-factor discriminant test — `test_bls.py` unit-checked; trial to $10^6$
with early abort), which converted hundreds of evidence-tier
*zero-jacobi* remainders into PROVED prime rungs; (iii) three scan waves:
frozen classes at $k\le200$ (97 classes), upgraded prover rerun (+20),
then **alternate aligned classes** — H needs only one class per cell;
Hensel-root $a$-lifts (L11) and $f$/$\varepsilon$ variants (ESC) closed
164 more, including every $w\equiv3\bmod4$ cell via the $a=1$ path, the
last three L11 holdouts via widened lotteries at $q\le2347$, and the
hardest ESC cell (67,(−1,1)) at $f=67$, $\varepsilon=-1$, $q=811$, $k=0$.
Member-level yield across waves: 197 prime-rung + 96 factorint closures
(293/293 = 100%); zero alarms (CLASS-CONTRADICTION/INCONSISTENT never
fired; Jacobi-tier rows never claimed). **Status:** step (ii) — one
emergent-free member of *each aligned class* — is now **verified
instance-wise on the whole 353-cell grid**; combined with the complete
witness coverage (L13f), the grid is hypothesis-free. What remains of H
is precisely the general statement for all cells $(w,z)$ (all odd $w$),
which is the analytic input the literature currently reaches only to
degree 2 unconditionally (Krumm, LitScout-pinned): the grid evidence
says the emergent-free member exists *in every tested class* — the
uniform theorem over all $w$ is exactly the open frontier.

**L15 [remainder law and counting law for H on the grid — PROVED per row,
statistics exact, 2026-08-18].** Independent exact audit of all 293 closure
members (`data/l15_remainders.jsonl`; every record re-verified end-to-end
through the L13f ladder; zero refusals, zero alarms):

- **Exact remainder law.** Strip $P(b)$ at the frozen set $S$, the
  $\{2,3,5,7\}$-shell and all primes $\le10^6$ (the `smooth_emergent`
  convention): the remainder $R$ is **squarefree in every one of the 293
  records**, and $R=1$ (pure-square rung) **never occurs**. Exactly two
  shapes: **prime** (197/293: $R$ a single prime $p$, never $p\cdot C^2$;
  digit count median 48, max 119) and **factorint** (96/293: $R$ a
  squarefree product of $e\in\{2,3,4\}$ emergent primes — counts 76/15/5;
  15 odd-$|E|$ rows all with every Hilbert symbol $+1$, the parity law
  holding to the letter). Deciding tension: the strip leaves *odd
  valuations only at unfrozen emergent primes*; frozen valuations of
  $P(b)$ itself are odd somewhere in 275/293 rows (absorbed by design).
- **Counting law** (`data/l15_density.json`, deterministic 2.8 s
  playback). First-member index $k_{\text{zero}}$ over the 293 classes:
  $k=0$ for **207/293 = 70.65%** (L11 65.3%, ESC 80.6%); median 0, mean
  25.1; tail ($k>0$, 86 classes) median 28, mean 85.6, max 1694. The
  per-prime-member emergent-free rate (MLE) is $p_0\approx0.040$ (L11
  0.035, ESC 0.072), i.e. closure waits $\approx15.7\times$ longer than
  the pure Bateman–Horn first-prime prediction for the same progression
  (geomean admissible density $\prod_{p\le1000}(1-1/p)=0.369$ over the
  prime-rung classes). The resistant subclass is $a=1$, $w\equiv3\bmod4$
  (26 classes: $k=0$ only 30.8%, both outliers 1694/1511 are here); the
  $a>1$ tower closes at $k=0$ in 80.6%.

**L16 [emergent-symbol structure — measured, per-row checks PROVED,
layering statement EVIDENCE, 2026-08-18].** Lead probe over 7 sampled
classes × $k\le600$ prime members (`data/l16_char.jsonl`, full
sign-decorated emergent lists, not just the ladder's −1 subset):

- The product character $\chi(k)=\prod_{p\text{ emergent}}(x_0,d_0)_p$ is
  identically $+1$ on aligned members (parity law; the frozen symbols are
  $+1$ by class construction). The discriminating content is per-prime.
- **Small layer is siev-form**: for every small emergent prime $p$
  recurring $\ge4$ times ($\le10^6$), in **all 7 sampled classes**, the
  symbol $(x_0,d_0)_p = (2\alpha b/p)$ is a function of the residue
  $k\bmod p$ alone (flag `residue_determined_mod_p` true on every row;
  analytically predictable since $b=b(k)$ is affine in $k$ — measured
  here with zero counterexamples). Hence the small-emergent obstruction
  at $p$ is a union of root classes of $k\bmod p$ with the bad sign
  lifted by **valuation parity**: the bad event is $\{v_p(P(b(k))) \text{
  odd}\}$ ∩ {bad sign}, and for a simple root its exact $p$-adic mass is
  $1/(p{+}1)$ (odd rungs $v_p=1,3,\dots$ of the Hensel tower). Singular
  ($P'\equiv0$) roots must be lifted recursively until each branch's
  odd/even mass is certified (with a residual bound) or excluded
  (product labeled partial). A sieve with $p$-adic local factors — not
  mere $\bmod p$ exclusion.
- **Two-layer obstruction split**: the fraction of prime members with
  *clean* small layer (no $-1$ emergent $\le10^6$) ranges 35–58% per
  sampled class; the full closure rate is $\approx4\%$ (L15) — so the
  big-cofactor layer ($p>10^6$, emergent primes forced singly or in
  small products by parity) carries its own symbol obstruction; the
  two-layer model rate ~ $0.45\times0.09$ is consistent with L15's
  measured $p_0\approx0.040$.
- Parity bookkeeping checked row-by-row: odd small-$-1$ count implies
  an odd cofactor-$-1$ count (parity law; every row consistent under the
  audit of L15's 293/293).

**L17 [sieve structure at full-grid scale + two-layer rate band —
per-row checks PROVED, model validated; rate claims band-only, 2026-08-19].**

- **Residue-determination, whole grid (FullGridSieve).** For all 293
  classes (prime members $k\le119$, 8,760 members): every recurring
  small emergent prime's symbol is $k\bmod p$-determined — **293/293,
  zero counterexamples** (`data/l17_sieve.jsonl`).
- **Mechanistic $p$-adic small-layer model, first fit (RateModel).**
  Exact odd-valuation bad masses $m_p$ in $\mathbb Z_p$ (simple root
  $1/(p{+}1)$; 3,746 step roots recursively lifted with certified
  residual bounds; 56 classes partial with 67 excluded (class,$p$)
  pairs — labeled): the product $\prod_p(1-m_p)$ tracks member-level
  small-layer clean rates. **The hit-conditioned matched evaluation is
  exact; the class-factor offset is a window artifact, characterized.**
  On the identical 8,760-member sample with closure-authority bad-root
  residues reconciled 5,375/5,375 (`data/l17_badroots_closures.jsonl`,
  v2 with exact resolved-mass prefixes on the 67 excluded singular
  pairs): marginal odd-hit expectation **4481.3 vs observed 4478
  (ratio 0.9993)**; hit-conditioned union expectation **3477.9 vs
  observed 3475**; marginal-preserving permutation null
  (2,000 sims, seed 20260819): mean 3475.3, observed $3475$,
  $z=-0.09$, empirical tails $0.61/0.53$ (observed sits at the null
  mean; two-sided tail clips to 1.0) — no dependence/overlap signal
  (`data/l17_matched.jsonl`, regenerated by `l17_matched.py`). The
  class-weighted Haar-product
  comparison still reads $\approx+0.0195$/3.6$\sigma$: that estimand compares
  $\mathbb Z_p$-uniform masses against a $k\in[0,119]$ window — the
  window-vs-Haar gap, characterized as such, not a model failure. The
  product model is the validated predictor of member-level small-layer
  obstruction (`data/l17_ratemodel_censored.jsonl`,
  `data/l17_sieve.jsonl`, `data/l17_badroots_closures.jsonl`).
- **Cofactor layer: proved statements are conditional/reciprocity;
  the rest is closure-conditioned description (CofactorSign).** Among
  the 293 CLOSURE rows (selected for all-$+1$ symbols): all 197
  prime-rung cofactor primes have symbol $+1$ auto-forced by
  reciprocity (proved, given the frozen/small $+1$ row structure), and
  the parity product is $+1$ on 293/293 rows (proved). The individual
  $+1$ signs of the 217 factorint cofactor primes and the residue
  tables (mod 4 / mod 8 / mod $d_0$-support) are **closure-conditioned
  descriptive data only — no distributional claim can be drawn from
  them** (`data/l17_cofactor.jsonl`); a nonclosure cofactor sample is
  required for structure claims beyond parity/reciprocity.
- **Two-layer rate: band-only.** On the predeclared wave-1 cohort with
  right-censoring and unknown trials counted both ways: per-family
  calibrated cofactor factor $c\in[0.017,0.80]$ — the refusal/
  zero-jacobi unknown mass dominates and the band does not resolve;
  predicted wave 1+2 closures $[114.3, 293.0]$ vs 117 observed;
  aggregated rate band $[0.010,0.472]$ brackets L15's measured
  $p_0=0.0400$. **Ledger discipline:** only the band is recorded; the
  L15 15.7× Bateman–Horn gap localizes mainly in the unresolved
  unknown-trial mass.

- **Horizon extension, off-grid (L17-horizon, 2026-08-19).** The
  verified horizon extends beyond the grid box ($w\le97$): four off-grid
  cells closed and independently replayed by the lead —
  $[101,[-1,1]]$ (L11: $a=5$, $\tau=51/101$, $\varepsilon=+1$, $q_1=13$,
  $N=339360$, member $k=16$, $Q=5429773$ prime, factorint rung);
  $[103,[-1,1]]$ (ESC: $f=103$, $\varepsilon=-1$, $q_1=41$, member
  $k=0$, $Q=41$, prime rung); $[107,[-1,1]]$ (ESC: $f=107$,
  $\varepsilon=+1$, $q_1=31$, $k=0$, $Q=31$, prime rung);
  $[109,[-1,1]]$ (L11: $a=71$, $\tau=10083/20165$, $\varepsilon=+1$,
  $q_1=193$, $k=0$, $Q=193$, prime rung; auxiliary tied/ramified
  budget-refused — exactly the accepted L14 closure form; a guarded
  cross-check of 1,275 class/member probes through $k=50$ found no
  tied=True+ramified-empty pair: `data/l17_horizon109_crossk.jsonl`).
  Off-grid machinery: fresh-L11 roots $1+4r^2\equiv0\ (w)$ +
  `_l10_class_cert`, or temporary `_L12_ESCAPE` registration +
  `l12_class_cert`, then the L13 ladder
  (`data/l17_horizon{101,103,107,109}.jsonl`, replay scripts
  `l17_horizon*.py` persisted in `math/h10q`; cross-k scan
  `data/l17_horizon109_crossk.jsonl`).

**Status.** Assembly's finitary form on the grid is **VERIFIED** (353/353
witnesses + 293/293 class members); the conditional $\forall_6$ record
rests on H's general (all-$w$ form, one emergent-free member per class),
which remains **OPEN** — the emergent obstruction now has a
matched-validated mechanistic $p$-adic sieve model (small layer);
the cofactor layer is proved at the reciprocity level with its
nonclosure sign distribution measurement-limited (11/400 resolved);
the two-sided rate band holds with the refused/undecidable mass as the
dominant uncertainty; no structural shortcut (L13c).



## L18 [STEP-(ii) DENSITY LAW + HORIZON WAVE 2, 2026-08-19]. The local masses are bounded, the clean product decays like $(\log X)^{-1/2}$, and the off-grid frontier splits into closures and one protocol wall

- **Local sieve bound (PROVED).** For every aligned class, with
  $F(k)=P(\varepsilon f(q_1+kN))$ the degree-8 class polynomial and
  $u(k)=2\alpha\varepsilon f(q_1+kN)$ the signature unit, the local
  emergent mass $m_p$ (Haar mass in $\mathbb Z_p$ of $\{k: v_p(F(k))$
  odd and $(u(k)|p)=-1\}$) satisfies
  $m_3=m_5=m_7=0$ (frozen alignment) and $m_p\le 8/(p+1)<1$ for all
  $p\ge11$, singular roots included with multiplicity. Independently
  re-verified by the lead on every stored pair: 102,439 local pairs
  with $p\le10^4$ across the 293 classes, **zero violations**
  (`data/l17_stepii.jsonl`, generator delivered by StepIIDensity).
- **No cutoff-uniform positive mass (PROVED for the 293 current
  classes).** Each class polynomial carries a replayed degree-8
  irreducibility certificate and a replayed witness that the signature
  unit is a nonsquare at a simple bad root, so the Chebotarev exponent
  of the bad-signature condition is $c=1/2$ in every class; hence
  $\prod_{p\le X}(1-m_p)\asymp(\log X)^{-1/2}$ and the small-layer clean
  mass tends to $0$. The hoped-for uniform positive lower bound does
  **not** exist. Measured exponents (EVIDENCE): root-harmonic 0.4778,
  product 0.4801, same-member sieve ($10^4\to10^5$) 0.4590.
- **Consequence for H (framing, not a theorem).** With
  $\#\{$prime members up to $X\}\sim X/(\log X)^2$ and clean rate
  $\asymp(\log X)^{-1/2}$, the expected number of emergent-free members
  diverges at rate $X/(\log X)^{5/2}$: H is a divergence statement, and
  the shortcut through a uniform constant is closed.
- **Horizon wave 2 (PROVED per row, lead-replayed).** Off-grid closures
  now cover $w\in\{101,103,107,109,113,127,137\}$ at $u=-1$: L11 route
  for $w=101,109,137$; escape route $f=w$ for $w=103,107,113,127$.
  New rows: $[113,[-1,1]]$ ESC $f{=}113$, $q_1{=}41$, $k{=}0$ (auxiliary
  tied/ramified budget-refused, accepted L14 form); $[127,[-1,1]]$ ESC
  $f{=}127$, $q_1{=}149$, $k{=}0$, tied True, ramified empty;
  $[137,[-1,1]]$ L11 $a{=}87$, $\tau=15139/30277$, $q_1{=}149$, $k{=}0$,
  tied True, ramified empty
  (`data/l17_horizon{113,127,137}.jsonl`, scripts `l17_horizon*.py`).
- **A protocol wall at $w=131$ (PROVED within the protocol).** For
  $w=131$ the L11 route is empty since $(-1|131)=-1$ kills
  $1+4r^2\equiv0$, and on the canonical escape route
  $\mathrm{Hilb}_5=-(5|q_1)$, so the frozen place 5 and the wild place
  $q_1$ cannot both be $+1$: over 336 proven-prime class attempts with
  $q_1\le1000$, 326 eligible candidates were obstructed, 0 aligned, 0
  refusals. The cell is **open**, the obstruction is a property of the
  protocol (`data/l17_horizon131.jsonl`). Open question: characterize
  the primes $w$ with this 5-symbol collision.
- **Composite $w$ is out of scope (structural verdict).** $w$ indexes a
  place; the construction breaks at the residue-field constructor for
  composite modulus. The downstream arithmetic runs after decomposing
  into genuine primes — illustrated by the pseudo-cell $[21,[-1,1]]$
  closing with $f=3$, $q_1=11$, $k=0$, factorint rung, tied True,
  ramified empty (`data/l17_composite_w.jsonl`). Scope clarification,
  not a horizon extension.
## L19 [SCHINZEL UPGRADE + THE CANONICAL 5-WALL AND ITS BREAK, 2026-08-19]. The analytic input is a named classical conjecture; the off-grid protocol wall is protocol-specific

- **Schinzel H implies the per-class clause of H (PROVED implication;
  conclusion CONDITIONAL).** Fix a verified aligned class and write
  $F(t)=P(\varepsilon f(q_1+Nt))=c\,G(t)$ with $c>0$ having the **square class** of $c$ supported on $S$ (literal $S$-support is FALSE for
  cells whose $z$ has denominator primes outside $S$ — counterexample
  $(w,z)=(3,3/11)$, $a=1$, $q_1=19$: $c=9072/15692141883605$ with outside part
  $11^{-12}$ — but the square class is what the parity/reciprocity step needs,
  since even outside valuations contribute symbol $+1$; corrected 2026-08-19)
  prime has $v_Q(P(b))=0$ since $P(b)\equiv-\delta a^4Z^4A^2 \bmod Q$);
  all other outside places have even valuation and symbol $+1$; the
  class certificate pins the frozen, moving and infinite symbols to
  $+1$ for every $t$; Hilbert reciprocity forces $(x_0,d_0)_R=+1$.
  Exact verification of the Schinzel hypotheses on all 293 canonical
  pairs (primitive irreducible degree-8 $G$ with replayed Frobenius
  certificates; positive leading coefficient; pair-product fixed
  divisor 1 by exact degree-9 finite differences; $G$ an $S$-unit for
  every integer $t$; 283 rows by the Taylor exponent, 10 alternate ESC
  rows by an exact normalized-polynomial congruence) plus a 24-row
  prime-rung cross-check (`data/l18_schinzel_implies_h.jsonl`,
  `l18_schinzel_implies_h.py`; lead-replayed).
  **Consequence:** the member-existence side of the conditional record
  is Schinzel/Bunyakovsky, not a bespoke hypothesis. Class *existence*
  per cell remains the separate structural side (L11h; verified
  293/293 on the grid).
- **Terminology correction (load-bearing).** Operational
  emergent-freeness is *no place with symbol $-1$*, not *no odd
  valuation outside $S$*: the kernel appends a place only when the
  valuation is odd AND the symbol is $-1$
  (`l13_filter.py` `smooth_emergent`/`cofactor_decide`), and the 197
  prime-rung closures have exactly one odd-valuation outside prime with
  reciprocity-forced $+1$. `CONDITIONAL.md` has been corrected; the
  stronger wording would not follow from Schinzel prime values.
- **The canonical 5-wall (PROVED for the canonical protocol).** With
  $e=v_5(D)$ one has $v_5(x_0)=-3-2e$ odd and $v_5(d_0)=0$, so
  $\mathrm{Hilb}_5=(d_0|5)=-(w|5)(5|q_1)$, while the wild place gives
  $\mathrm{Hilb}_{q_1}=(5|q_1)$. Both can be $+1$ iff $(w|5)=-1$.
  Combined with "L11 route available iff $w\equiv1\bmod4$", the
  canonical protocol is obstructed exactly for
  $w\equiv3\bmod4$ and $(w|5)=+1$, i.e. $w\equiv11,19\bmod 20$.
  Confirmed on new data: canonically obstructed $w=131,139,151,179$;
  canonically closed
  $w=101,103,107,109,113,127,137,149,157,163,167,173$
  (`data/l18_horizon_sweep_{a,b}.jsonl`, `data/l17_horizon131.jsonl`).
- **The 5-wall is protocol-specific, not intrinsic (PROVED).** The
  obstructed cell $[131,[-1,1]]$ closes off the canonical protocol with
  $a=7$ (so $A=197$), $\tau=0$, $f=131$, $\varepsilon=+1$, $q_1=Q=41$,
  $N=45046449133189418880$, $k=0$, prime rung, tied True, ramified
  empty.

  *Why no wall excludes this row (stated exactly).* **L10b** ($\tau=0$
  wall, L10b above, lines 1016-1021) assumes verbatim "$b=\varepsilon
  q_1$ with $q_1\ne A$ prime satisfying L9a's factor-by-factor units".
  Every *other* L10b hypothesis does hold here — $a=7$ is odd, $A=197$
  is prime, $A\nmid2zD_z$ — and the single failing clause is the shape
  of $b$: $b=131\cdot41$ is not $\pm$ a prime. **L12** (the generalized
  wall, L12 header line 1275) assumes "$b$ coprime to the cell data";
  $f=w=131$ divides $\operatorname{num}(z)$, so that hypothesis fails
  too. This is precisely the L12b non-coprimality door, and the
  $b$-shape $b=\varepsilon fQ$ with $f\mid\operatorname{num}(z)$ is the
  orthodox escape shape already accepted in 103 grid rows: **the
  novelty of this row is the $\tau=0$ alignment alone**, not a new
  $b$-family.

  *Certificate authority.* The row verdict uses the kernel and nothing
  else: `h10q._l7_tied_status`, `h10q.ramified`, and the identity
  $M_0=P(b)/(b^4D_z^2A^2)$ asserted against `h10q._sun_h` and
  `h10q._l10_P` (`l18_route131.py:620-631`) — the same authority as
  every canonical row. The *class* certificate is
  `fixed_factor_class_cert`, a local generalization of
  `l12_class.l12_class_cert` with two restrictions removed ($a$ need
  not be $1$; $f$ need not be prime, every prime of $f$ being frozen),
  same exponent-lemma proof on wider inputs
  (`l18_route131.py:118-132`, replayed at `close_target`); it is
  **not** the byte-identical kernel class function, and is labelled as
  such here.

  The same
  $w$-column also closes at $[131,[5,1]]$ ($a=1$, $q_1=29$, factorint
  rung). Bounded searches: composite-$f$ 378 attempts 0 aligned
  (EVIDENCE); enlarged-$S$ proved universally ineffective (168 probes);
  $w=179$ probe at $a=7,\tau=0$ found 4 aligned classes but no closure
  within budget (EVIDENCE) (`data/l18_route131.jsonl`,
  `l18_route131.py`; lead-replayed).
- **The canonical 5-wall is ALWAYS breakable on the class side (PROVED).**
  The exact break criterion, derived factor by factor: for $p^e\Vert A$,
  $(x,d)_p=1$ when $e$ is even and $(-2\varepsilon f q_1\mid p)$ when $e$
  is odd, whence
  $$\prod_{p\mid A}(x,d)_p\cdot(A\mid q_1)=(-2\varepsilon\mid A)(f\mid A).$$
  For admissible odd $a$ one has $A\equiv5\bmod8$ and
  $(-2\varepsilon\mid A)=-1$ for both signs, so **L10b's anti-correlation
  breaks exactly when $(f\mid A)=-1$** — the extra Legendre factor
  $(f\mid A)$ is precisely what the prime-$b$ family lacks. With $f=w$,
  $z=-w$ one also gets $v_w(x)=0$, $v_w(d)=1$, $(x,d)_w=(-A\mid w)$, and
  the same condition $(w\mid A)=-1$ makes that symbol $+1$.
  *Existence, uniformly in $w$:* for prime $w\equiv3\bmod4$ the character
  sum $\sum_{r\bmod w}(1+4r^2\mid w)=-1$ has no zero terms, so exactly
  $(w+1)/2$ residues give $(A\mid w)=-1$; CRT with $a\equiv0$ mod
  $\operatorname{rad}(w^3+2)$ gives an odd lift preserving that character
  and forcing $\gcd(A,D)=1$ (since $4D\equiv(Z-2)^2\bmod p$ for $p\mid A$,
  $Z=-w^3$); a second CRT plus Dirichlet supplies infinitely many primes
  $q_1$ aligning every frozen place, the $2$-adic criterion
  $\varepsilon wq_1\equiv a\bmod4$, and the real place. Hence **for every
  prime $w\equiv11,19\bmod20$ an aligned $\tau=0$, $f=w$ class exists.**
  Audit: 16,744 candidates / 14,142 eligible / 908 fully aligned bases,
  **0 symbol, 0 derivation and 0 status mismatches**
  (`data/l19_tauzero.jsonl`, `l19_tauzero.py`; lead-replayed, 14.7 s).
  *Scope, exactly.* This is the class side. Member-level closure still
  needs the L19 Schinzel input, so unconditional closure of an
  obstructed cell remains **OPEN**. And the **fixed** $a=7$ subfamily is
  not universal: it is breakable iff $(w\mid197)=-1$ — bounded check
  breaks $131,139,151,179,199,211,271,359$ and leaves
  $191,239,251,311,331,379$ obstructed *for that $a$*. The three $k=0$
  closures are EVIDENCE for the family, not a universal member theorem.
- **Off-grid horizon (PROVED per row, all lead-replayed).** Closed:
  $w\in\{101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,
  179\}$ at $u=-1$, plus $[131,[5,1]]$ — **every off-grid cell tried is
  now closed; none is known to be intrinsically obstructed.** The last
  one, $w=179$, closes at $a=7$, $\varepsilon=-1$, $f=179$, $q_1=Q=251$,
  $N=1459660664008632960$, $k=0$, prime rung, exactly as the fixed-$a$
  criterion $(179\mid197)=-1$ predicts; protocol 1,872 class attempts,
  178 aligned classes, all searched over $k\le60$, **decided fraction
  10,282/10,858 $=94.7\%$**, largest cofactor attempted 222 digits, 576
  refusals never counted as evidence (`data/l19_cell179.jsonl`,
  `l19_cell179.py`; lead-replayed, 168 s). Routes: L11 for
  $101,109,137,149\ (a{=}127),157\ (a{=}457)$; canonical ESC for
  $103,107,113,127,163,167,173$; **non-canonical $a{=}7,\tau{=}0,f{=}w$
  for every canonical wall met so far** — $131$
  ($N=45046449133189418880$), $139$ ($N=1236465506993834880$), $151$
  ($N=16835857850229022080$), each at $q_1=41$, $k=0$, prime rung, tied
  True, ramified empty. $w=179$ was the last open one and is now **closed**
  too ($a=7$, $\varepsilon=-1$, $q_1=Q=251$; horizon bullet below,
  `data/l19_cell179.jsonl`)
  (`data/l18_horizon_sweep_{a,b}.jsonl`, `data/l18_route131.jsonl`).
  **Pattern:** the canonical two-route wall at $w\equiv11,19\bmod20$ has
  been broken in every instance tested, so no off-grid cell is yet known
  to be intrinsically obstructed.
- **Divergence model (EVIDENCE).** With degree-8 height
  $\log X_k=8\log k+O(1)$, the Bateman-Horn $\times$ L18 hazard is of
  order $(\log k)^{-3/2}$ and the expected count grows like
  $K/(\log K)^{3/2}$: accumulated divergence, explicitly not uniform
  positive clean mass. Calibrated on a predeclared SHA split
  (239 fit / 54 held out): $C=0.9011$ (95% CI $[0.794,1.017]$); the
  L18 tail explains a $8.5\times$ conditional delay versus pure
  Bateman-Horn against the ledger's $15.67\times$, closing 51% of the
  excess multiplier and 78% of the log gap, residual $1.84\times$; the
  stationary model misses the selected $k=0$ spike (observed 70.65% vs
  predicted 16.67%, KS $D=0.497$) (`data/l18_divergence_model.jsonl`,
  `l18_divergence_model.py`).

## L20 [CLASS EXISTENCE IS A THEOREM, 2026-08-19]. Clause (ii) leaves the hypothesis

**Theorem (uniform class existence; PROVED).** Let $(w,z)$ be any cell of
the family, $w$ an odd prime with $v_w(z)\ge1$. Take the fixed factor
$f=w$ (always available, since $v_w(z)\ge1$ forces $w\mid\operatorname{num}(z)$).
Choose an odd $a$ with $A=1+4a^2$ and $(A\mid w)=-1$; such $a$ exists
because
$$\sum_{r\bmod w}\Bigl(\frac{1+4r^2}{w}\Bigr)=-1,$$
a standard quadratic character sum with no vanishing terms, so exactly
$(w+1)/2$ residues give character $-1$. Put
$S=\{2,3,5,7\}\cup\operatorname{supp}(\alpha)\cup\operatorname{supp}(\delta)\cup\{w\}$
and $M=4A\prod_{p\in S}p$. Then the system
$$q_1\in(\mathbb Z/M\mathbb Z)^\times,\qquad
  (2wq_1\mid p)=+1\ \ \text{for every odd }p\in S\setminus\{w\},$$
has exactly $\varphi(M)/2^{\#\{p\in S\ \text{odd},\,p\ne w\}}$ classes —
in particular it is **nonempty** — and Dirichlet supplies infinitely many
primes $q_1$ in it. The exponent-lemma certificate then yields an
**aligned class** for the cell.

*Why $f=w$ never fails.* With $Z=z^3$, $v_w(Z)\ge3$ gives
$D=1-Z-a^2Z^2\equiv1\bmod w$, and $(A\mid w)=-1$ forces $w\nmid aA$, so
$\delta$ is a $w$-unit: all fixed-factor guards hold, with no
denominator, gcd or branch exception (including $w=5$, where $a=3$,
$A=37$ works).

*Local symbols, exactly.* On the square branch with $\varepsilon=+1$: $x$
is a square at $w$; at $2$ one has $v_2(x)=6$, $v_2(d)=3$,
$u_d\equiv b$, $u_x\equiv1-2As^2b \pmod 8$, and the $2$-adic Hilbert
exponent vanishes for every odd $b$ class (so the symbol is $+1$ whether
or not $s=(a-1)/2$ is even); $d$ is a square at every other frozen odd
$p$; and since every $p\mid A$ is $1$ mod $4$,
$(A\mid q_1)=(2w\mid A)=(2\mid A)(w\mid A)=(-1)(-1)=+1$.

**Verification.** 103/103 canonical escape rows: residue systems proved
nonempty, 1,113,000 residues directly enumerated across 7 exact moduli,
and 103/103 generalized certificates reproduce the recorded
$N,S,k_s,Q_0$, excluded set and `ok` of `data/l12b_class_sample.json`.
**353/353 grid cells** carry a clean generalized certificate with $f=w$.
Composite $A=325=5^2\cdot13$ replayed; 16/16 odd $(a,b)$ residue pairs
checked against `h10q.hilbert`; zero refusals
(`data/l19_classexist.jsonl`, `l19_classexist.py`; lead-replayed, 73 s).

**Independent corroboration of the 5-wall.** With $a$ *fixed* to $1$ the
construction collides at exactly $w=11,19,31,59,71,79$ — precisely the
grid primes $w\equiv11,19\bmod20$, the canonical 5-wall set derived
independently in L19 by the Hilbert-symbol route. Every collision is
repaired by an alternate $a$; the bounded pool $a\in\{1,3,5,7\}$ leaves
zero empty cells.

**Scope, exactly.** This removes clause (ii) — *existence of one verified
aligned class per cell* — from the hypothesis, for every cell, not just
the 103 escapes. It does **not** prove: irreducibility or fixed-divisor
$1$ for the newly constructed pair, a Schinzel prime value, an
emergent-free member, the L6 assembly lemma, or H10 over $\Q$. Member
existence remains **CONDITIONAL** on Schinzel H via L19.
