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
Here $F(z)\Rightarrow z\in\bigcup_{w\notin S}\mathfrak m_w$ is proved
unconditionally.  At the L6 stage the reverse implication was the open
assembly lemma; L22d now proves it **CONDITIONAL on classical Schinzel H
alone** (`THEOREMS.md`, L19–L22).  All displayed rational expressions
are cleared as in the next paragraph.

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


**Current status (updated by L22).** Architecture and exact count:
**PROVED**; target-place theory and norm obstruction: **PROVED**
(W0–W2).  The global assembly implication is now **PROVED from classical
Schinzel H** by L19–L22.  Hence, conditionally on that unproved classical
conjecture, the main union outside $S$ is $\exists_6$; adjoining the
finitely many omitted $\exists_3$ maximal ideals by finite disjunction
keeps $\exists_6$, so $\mathbb Z$ is $\forall_6$ and
$\operatorname{efd}_{\mathbb Q}(\mathbb Q\setminus\mathbb Z)\le5$
(`CONDITIONAL.md`, §1; `THEOREMS.md`, L22d).  E6 remains bounded
EVIDENCE and is not used in that implication.  This is independent of
whether Sun's own §§3–8 chain survives refereeing.


---

## L7 [STRUCTURAL COROLLARIES + BOUNDED PROBE, 2026-08-15]. What the then-open assembly problem looked like

**Historical checkpoint, SUPERSEDED as frontier by L19–L22.** L7 did
not close L6 at that stage; its structural corollaries and machine checks
remain valid (`h10q.py::_verify_L7`).

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

**Historical checkpoint, SUPERSEDED as frontier by L19–L22.** L8 did
not close L6 at that stage.  L8a–c remain unconditional theorems and
L8d remains scoped EVIDENCE (`h10q.py::_verify_L8`).

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

**Historical consequence, SUPERSEDED by L20 and L22.** At L11, step (i)
held for $190$ of the $353$ recorded cells and the prime-$b$ family was
impossible at the remaining $163$ (`data/l11_classes.jsonl`;
`THEOREMS.md`, L11h).  L20 later supplies the non-coprime $f=w$ class
for every cell, and L22 supplies vertical irreducibility
(`/tmp/l22_elimination.md`).  Thus this family wall is not a current
assembly gap.
Separate historical steering evidence also closed the recorded
$(53,\tfrac13)$ and $(29,-2)$ cells by wider-parameter searches; those
are finite per-row certificates, not part of the uniform proof
(`data/l9_steer_run.jsonl`; `THEOREMS.md`, L11 session record).

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

**Historical scope, superseded first by L13a and then by L20–L22.** On
the recorded L12 date this was base-point alignment only.  L13a verified
the minimal controlled set and class constancy; L20 later constructed an
aligned $f=w$ class for every cell, and L22 closed irreducibility
(`data/l12b_class_sample.json`; `data/l19_classexist.jsonl`;
`/tmp/l22_elimination.md`).

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
fired; Jacobi-tier rows were never claimed).  **Historical L14 status,
SUPERSEDED as frontier by L19–L22:** step (ii) was verified
instance-wise on the whole $353$-cell grid, while the all-cell statement
was then open (`data/l13h_all_closures.json`).  L20 now supplies one
class for every cell, L22 supplies vertical irreducibility, and L19
derives a member from classical Schinzel H (`THEOREMS.md`, L19–L22).
The grid remains finite corroboration only.

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

**Historical L17 status, SUPERSEDED as frontier by L19–L22.** Assembly's
finitary form on the grid was verified at $353/353$ witnesses and
$293/293$ class members (`data/l13h_all_closures.json`); the all-cell H
statement was then open.  The mechanistic sieve and rate statements
below retain their PROVED/EVIDENCE scopes, but the current implication
uses L20 class existence, L22 irreducibility, and L19 plus classical
Schinzel H (`THEOREMS.md`, L19–L22).



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

## L20 [CLASS EXISTENCE IS A THEOREM, 2026-08-19]. Clause (ii) is removed; L22 later completes condition (a)

**Theorem (uniform class existence; PROVED).** Let $(w,z)$ be any cell of
the family, $w$ an odd prime with $v_w(z)\ge1$. Take the fixed factor
$f=w$ (always available, since $v_w(z)\ge1$ forces $w\mid\operatorname{num}(z)$).
Choose an odd $a$ with $A=1+4a^2$ and $(A\mid w)=-1$; such $a$ exists
because
$$\sum_{r\bmod w}\Bigl(\frac{1+4r^2}{w}\Bigr)=-1,$$
a standard quadratic character sum.  It has
$1+\bigl(\frac{-1}{w}\bigr)$ zero residues, and exactly
$$\frac{w-\bigl(\frac{-1}{w}\bigr)}2$$
residues give character $-1$: $(w+1)/2$ when $w\equiv3\bmod4$ and
$(w-1)/2$ when $w\equiv1\bmod4$.  In particular the set is nonempty. Put
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
the 103 recorded escapes (`data/l19_classexist.jsonl`). L20 itself does
not prove irreducibility or a Schinzel prime value; L22 below supplies
fixed-fibre irreducibility, and member existence remains **CONDITIONAL**
on classical Schinzel H through L19.

## L21 [THE REDUCIBLE LOCUS AND WHY THE CONSTRUCTION AVOIDS IT, 2026-08-19]

At the **historical L21 checkpoint, SUPERSEDED by L22**, the L20
admissibility audit had left one uniform gap: irreducibility of the
degree-$8$ class polynomial $P$.  Blanket irreducibility is false on the
explicit locus below, while the construction's fixed branch avoids it.

**L21a [reducible locus; PROVED].** Suppose $s=0$ (equivalently $a=1$) and
$\delta_\tau=\sigma^2$ is a rational square. Then the $b^5$ term
$-32A^3s^2D^2b^5$ vanishes and $P$ becomes a difference of two squares:
$$P(b)=\bigl(4DAb^2\bigr)^2-\bigl(\sigma a^2Z^2N_g(b)\bigr)^2
      =\bigl(4DAb^2-\sigma a^2Z^2N_g\bigr)\bigl(4DAb^2+\sigma a^2Z^2N_g\bigr),$$
a genuine $4\times4$ factorization over $\mathbb Q$. The case $\tau=0$
($\delta=1$, $A=5$) gives $P=(20Db^2-Z^2N_g)(20Db^2+Z^2N_g)$. Verified
exactly on 56/56 applicable instances, with 8/8 control rows ($s\ne0$)
correctly inapplicable (`data/l21_reducible_locus.jsonl`, `l21_reducible_locus.py`).
*Discovery credit:* the $\tau=0$, $a=1$ instance was found by the IrredGeneric
agent's reducibility hunt; the general $\sigma$ form and the proof are the
lead's, independently replayed.

**L21b [the constructed branch escapes, uniformly; PROVED].** On
$\tau^\dagger=(1+2a^2)/A$,
$$\delta_{\tau^\dagger}=1-A\Bigl(\tfrac{1+2a^2}{A}\Bigr)^2
  =\frac{A-(1+2a^2)^2}{A}=\frac{1+4a^2-1-4a^2-4a^4}{A}=-\frac{4a^4}{A}<0$$
for every nonzero $a$. A negative rational is never a rational square, so the
L21a degeneration **cannot occur on the constructed branch, for any $a$ and
any cell.** Moreover $\delta_{\tau^\dagger}\cdot(-A)=(2a^2)^2$, so the square
class of $\delta_{\tau^\dagger}$ is that of $-A$: the relevant quadratic field
is the *imaginary* $\mathbb Q(\sqrt{-A})$. Checked 9/9; zero constructed-branch
instances fall in the L21a locus across the full $(a,z)$ sample.

**L21c [near-palindromic structure; PROVED on sample].** $P$ is palindromic up
to a single monomial:
$$P(b)+32A^3s^2D^2b^5\ \text{satisfies}\ P_i=P_{8-i}\ (0\le i\le8),$$
and $P_0=P_8$ unconditionally. Verified on 288/288 rows spanning
$a\in\{1,3,5,7,9,11,25,-3,-7\}$, eight $z$ values and four branches. Hence for
$s=0$ the polynomial is exactly palindromic, so $P(b)=b^4Q(b+1/b)$ with
$\deg Q=4$ and the Galois group embeds in $C_2\wr S_4$ — consistent with the
recorded $D_8\wr C_2$ evidence (L13c) — and the whole degree-8 irreducibility
question reduces there to a degree-4 question plus one quadratic condition.

**Historical L21 frontier (SUPERSEDED by L22).** L21 reduced the remaining
algebraic question to irreducibility on the fixed branch
$\tau^\dagger=(1+2a^2)/A$, $\delta_{\tau^\dagger}=-4a^4/A$, after excluding
the displayed square-$\delta$ degeneration. L22 below proves that fixed-fibre
statement for every nonzero rational $Z$ (`/tmp/l22_elimination.md`).

**L21d [generic irreducibility on the constructed branch; PROVED].** Let
$H=(A/4)P\in\mathbb Q[a,Z][b]$ on $\tau^\dagger$. Its leading and constant
coefficients are both $L=a^8A^2Z^4$. If $H/L$ factored over $\mathbb Q(a,Z)$,
the monic factors would lie in the integrally closed localization
$\mathbb Q[a,Z,1/L]$, and specializing $(a,Z)=(1,27)$ would preserve their
degrees; but the primitive specialized octic is irreducible mod $17$
($x^{17^8}\equiv x$ and $\gcd(x^{17^4}-x,q)=1$). Hence **$P$ is irreducible of
degree 8 over $\mathbb Q(a,Z)$.**

*From generic to per-cell, correctly.* A two-variable thin set may contain a
whole vertical line $Z=Z_{\text{cell}}$, so generic irreducibility alone does
not settle a fixed cell. The repair is one-variable: for each of the 353
target $z$, an exact modular certificate at $a=1$ proves $P(a,z)$ irreducible
in $\mathbb Q(a)[b]$; then quantitative Hilbert irreducibility (Cohen–Serre,
$O_z(\sqrt B\log B)$ bad integers up to $B$) against the admissible
progression's $B/M+O(1)$ members shows the progression contains a good $a$.
So **for every one of the 353 cells, an admissible $a$ with $P$ irreducible
provably exists** — and in practice the *first* admissible $a$ always works:
48/48 rows (32 grid, 16 off-grid) had first index 1.
Hunt: 1024/1024 constructed-branch rows proved irreducible mod a prime, **0
reducible**; the only reducible rows found anywhere were the 128 off-branch
$a=1$ square-$\delta$ cases of L21a (`data/l21_irred_generic.jsonl`,
`l21_irred_generic.py`; lead-replayed, 101 s).

**L21e [structural no-go laws and the $a=1$ quartic; PROVED].** Three negative
results that close off the cheap routes, plus one positive:
- **Mod $w$ is useless.** Uniformly, $P\equiv16A^2b^4(1-2As^2b)\pmod w$, so the
  factorization type is $1^4$ or $1^4\!\cdot\!1$ and $P$ is *never* an
  irreducible octic mod $w$. (This refutes the lead's suggested
  $(A\mid w)=-1\Rightarrow$ mod-$w$ irreducibility route.)
- **The $w$-Newton polygon** has several slopes, never a single slope of
  denominator 8; and for odd $p\mid A$ with $v_p(Z)=0$ it has exactly two
  length-4 slopes $\mp v_p(A)/2$, ruling out the $p\mid A$ single-slope route.
- **Reciprocity is exact:** $P-b^8P(1/b)=32A^3s^2D^2(b^3-b^5)$, so $P$ is
  reciprocal **iff** $a=1$, with no nontrivial rational reciprocal twist
  (this sharpens L21c to an iff).
- **Positive:** for $a=1$ and $(5\mid w)=-1$ the degree-4 trace polynomial is
  *uniformly* irreducible by residue obstructions over $\mathbb Q(\sqrt{-5})$,
  and the octic is irreducible whenever
  $\Xi(Z)=(125D^2+64Z^4)(125D^2+1024Z^4)$ is a nonsquare — which holds on all
  189 canonical $a=1$ rows.
353 pairs, 0 mod-$w$ anomalies, 0 Newton anomalies, 0 trace-norm anomalies,
353 constructed-branch modular certificates, 0 reducible
(`data/l21_irred_direct.jsonl`, `l21_irred_direct.py`; lead-replayed, 4 s).

**Historical L21 endpoint (SUPERSEDED by L22).** At L21 the vertical
statement was proved only for the 353 recorded target values, using exact
certificates and quantitative Hilbert irreducibility
(`data/l21_irred_generic.jsonl`). L22 now proves it for every fixed
$Z\in\mathbb Q^\times$ on
$\tau^\dagger=(1+2a^2)/A$, $\delta_{\tau^\dagger}=-4a^4/A$; no new-cell
certificate remains in the current chain (`/tmp/l22_elimination.md`,
`data/l22_elimination.jsonl`).

## L22 [FIXED-FIBER IRREDUCIBILITY ON THE CANONICAL BRANCH — PROVED]. The last non-Schinzel obligation closes

Throughout this block the branch is **fixed**:
$$
\tau=\tau^\dagger=\frac{1+2a^2}{A},\qquad
\delta_{\tau^\dagger}=-\frac{4a^4}{A},\qquad A=1+4a^2.
$$
(Source: `/tmp/l22_elimination.md`, §§1–5; exact replay
`data/l22_elimination.jsonl`.)  Put
$$
D=1-Z-a^2Z^2,\qquad N_g=16a^4b^2-A(b-1)^4,
$$
and, on this same $\tau^\dagger$ branch,
$$
H=\frac A4P
 =a^8Z^4N_g^2+4A^3D^2b^4-2A^4(a-1)^2D^2b^5.
$$
(Source: `/tmp/l22_elimination.md`, lines 5–20; kernel identity independently
checked in `data/l22_elimination.jsonl`.)

### L22a. Uniform fixed-$Z$ theorem and exact elimination

**Theorem (PROVED).** For every **fixed** nonzero rational $Z$,
$H(a,Z,b)$, and therefore $P(a,Z,b)=(4/A)H(a,Z,b)$, is irreducible in
$\mathbb Q(a)[b]$.  In particular this applies to every cell, because
$Z=z^3\ne0$ when $v_w(z)\ge1$ (`/tmp/l22_elimination.md`, §§5–6;
`data/l22_elimination.jsonl`, record `fixed-Z-theorem`).

The word *fixed* is indispensable: $Z^4$ is then a unit of
$\mathbb Q[a]$, which is what makes the endpoint divisor enumeration
below valid.  Also indispensable are $Z\ne0$, the exact equality of the
two endpoint coefficients, the irreducibility of $A=1+4a^2$ over
$\mathbb Q$, $\gcd(a,A)=1$, primitivity and Gauss's lemma, and the fact
that the residue fields at $a=0$ and $a=\infty$ are $\mathbb Q$
(`agent://ElimAudit`; proof source `/tmp/l22_elimination.md`, §§1–5).

**Exact local data.** The constant and leading $b$-coefficients are both
$$L=a^8A^2Z^4.$$
For $Z\notin\{0,1\}$, $H$ is primitive in $\mathbb Q[a][b]$: a common
divisor must divide $L$, while the $b^4$ coefficient at $a=0$ excludes
$a$ and reduction modulo $A$ excludes $A$ (`/tmp/l22_elimination.md`,
§1).  Direct sparse expansion gives
$$
\deg_a(H_i)_{i=0}^8=(12,12,14,14,16,13,14,12,12),
$$
$$
v_a(H_i)_{i=0}^8=(8,8,8,8,0,0,8,8,8),\qquad
v_A(H_i)_{i=0}^8=(2,2,1,1,0,1,1,2,2).
$$
(Source for both vectors and the cancellation leaving
$1024a^{13}Z^4b^5$: `/tmp/l22_elimination.md`, §1;
`data/l22_elimination.jsonl`, records `local-data` and
`global-degree-filter`.)

At the prime $A$, the Newton polygon has two length-$4$ sides of slopes
$-1/2$ and $+1/2$.  A factor inherits a sub-multiset of these slopes,
and integral endpoint valuations force its degree to be even.  Thus a
reducible octic has a factor of degree $2$ or $4$
(`/tmp/l22_elimination.md`, §1; `data/l22_elimination.jsonl`,
`global-degree-filter`).

At $a=0$, for $Z\notin\{0,1\}$, the sides are
$$4@(-2),\qquad1@0,\qquad3@(8/3),$$
and the four-root residual after $b=a^2x$ is
$$R_0(x)=Z^4+4(1-Z)^2x^4.$$
It has no rational linear factor because both summands are nonnegative
on $\mathbb Q$, and the denominator-$3$ side is an indivisible local
block.  Hence a quadratic global factor must take two small roots; a
quartic takes either all four small roots or the three large roots and
the unit root.  In a $4+4$ split, after swapping factors, one factor is
therefore the all-small factor (`/tmp/l22_elimination.md`, §1;
`data/l22_elimination.jsonl`, `global-degree-filter`).

At infinity the two exact rescalings are
$$
a^{-12}H(a,x/a)\longrightarrow16Z^4(4x^2-1)^2,\qquad
a^{-20}H(a,ax)\longrightarrow16Z^4x^4(x^2-4)^2.
$$
(Source: `/tmp/l22_elimination.md`, §1;
`data/l22_elimination.jsonl`, `local-data`.)

**Endpoint parameterization.** Let $F$ be a primitive degree-$m$ factor,
where $m\in\{2,4\}$, containing only small roots at $a=0$.  Unique
factorization of the common endpoint $L$ gives
$$
\operatorname{lc}_b(F)=cA^\beta,\qquad
F(0)=d\,a^{2m}A^\delta,
$$
with $c,d\in\mathbb Q^\times$ and
$0\le\beta,\delta\le2$.  If $\ell$ roots of $F$ are large at infinity,
endpoint comparison gives
$$2m+2\delta-2\beta=2\ell-m.$$
The exact enumeration leaves three quadratic endpoint types and, for
$m=4$, only $(\beta,\delta,\ell)=(2,0,4)$
(`/tmp/l22_elimination.md`, §2; `data/l22_elimination.jsonl`,
`endpoint-enumeration`).

**Degree $2$ is impossible.** Infinity forces the constant/leading ratio
of the quadratic residual to be $\rho=\pm16$.  A monic factorization of
the zero-place residual has the form
$$
(x^2+ex+\rho)(x^2-ex+s)=x^4+K,\qquad
K=\frac{Z^4}{4(1-Z)^2}>0.
$$
Coefficient comparison gives $e(s-\rho)=0$ and
$s-e^2+\rho=0$.  If $e=0$, then $K=-\rho^2<0$; otherwise
$s=\rho$ and $e^2=2\rho=\pm32$, impossible over $\mathbb Q$
(`/tmp/l22_elimination.md`, §3; `data/l22_elimination.jsonl`,
`degree-2-elimination`).

**Degree $4$ is impossible.** The sole endpoint type is
$$\operatorname{lc}_b(F)=cA^2,\qquad F(0)=da^8,$$
and the infinity residual forces $d=256c$.  Matching the all-small
residual at $a=0$ then forces
$$Z^4=1024(1-Z)^2.$$
Factoring the difference of squares gives
$$Z^2+32Z-32=0\quad\text{or}\quad Z^2-32Z+32=0,$$
whose discriminants are respectively
$1152=576\cdot2$ and $896=64\cdot14$, neither a rational square.
Therefore no rational $Z$ permits a quartic factor
(`/tmp/l22_elimination.md`, §4; `data/l22_elimination.jsonl`,
`degree-4-elimination`).

**Exceptional values and specializations.**

- $Z=0$ is genuinely reducible:
  $$H=2A^3b^4\bigl(2-A(a-1)^2b\bigr).$$
  It is excluded because a cell has $z\ne0$
  (`/tmp/l22_elimination.md`, §5;
  `data/l22_elimination.jsonl`, exceptional-fibre record `Z=0`).
- At $Z=1$ the $a=0$ polygon degenerates.  At $a=1$ the primitive
  coefficient vector is
  $$(25,-200,540,-760,1546,-760,540,-200,25),$$
  and its reduction modulo $11$ satisfies the exact degree-$8$ Rabin
  criterion.  Localization at the nonzero endpoint $L$ makes any
  putative vertical factorization specialize with positive degrees,
  contradicting that certificate (`/tmp/l22_elimination.md`, §5;
  `data/l22_elimination.jsonl`, exceptional-fibre record `Z=1`).
- A specialization with $D=0$ makes $H=a^8Z^4N_g^2$, but for fixed
  $Z\ne0$ there are at most two such $a$-values.  Vertical
  irreducibility is a statement over $\mathbb Q(a)$, and the admissible
  progression omits these finitely many values
  (`/tmp/l22_elimination.md`, §5;
  `data/l22_elimination.jsonl`, `fixed-Z-theorem`).

This proves the theorem symbolically; the finite-field tables in
`data/l22_elimination.jsonl` are **EVIDENCE** only and are not used to
infer any rational fibre.

### L22b. Independent reciprocal all-$Z$ theorem on the same branch

**Independent theorem (PROVED).** Fix
$$
a=1,\quad A=5,\quad s=0,\quad
\tau=3/5=\tau^\dagger,\quad\delta_{\tau^\dagger}=-4/5.
$$
On this explicitly named branch, the normalized polynomial
$P(1,Z,b)$ is irreducible over $\mathbb Q$ for every
$Z\in\mathbb Q^\times$ (`/tmp/l22_reciprocal_cube.md`, Status and Final
conclusion; `data/l22_reciprocal_cube.jsonl`, summary).  The
$\tau^\dagger,\delta=-4/5$ qualification is load-bearing: the
square-$\delta$ branches of L21a are reducible.

Put $D=1-Z-Z^2$ and $u=b+b^{-1}$.  For rational $Z$, $D$ cannot
vanish because $Z^2+Z-1$ has nonsquare discriminant $5$
(`/tmp/l22_reciprocal_cube.md`; `agent://ReciprocalAudit`).  Exact
expansion on this same branch gives
$$
P(b)=b^4T(u),\qquad
T(u)=400D^2+\frac45Z^4\bigl(16-5(u-2)^2\bigr)^2.
$$
Over $E=\mathbb Q(t)$, $t^2=-5$,
$$
\frac54T=Q_-(u)Q_+(u),\qquad
Q_\pm=Z^2\bigl(16-5(u-2)^2\bigr)\pm10Dt.
$$
(Source: `/tmp/l22_reciprocal_cube.md`, §1;
`data/l22_reciprocal_cube.jsonl`, symbolic-identity record.)
Trace reducibility is equivalent to
$80Z^2+50Dt$ being a square in $E$, and therefore forces
$$125D^2+64Z^4\in\mathbb Q^{\times2}.$$
Writing $Z=p/q$ in lowest terms and
$d=q^2-pq-p^2$, its cleared numerator is
$$125d^2+64p^4\equiv5\pmod8,$$
which is never a square.  Thus $T$ is irreducible for every rational
$Z\ne0$ (`/tmp/l22_reciprocal_cube.md`, §1;
`data/l22_reciprocal_cube.jsonl`, trace-theorem record).

For a root $u$ of $T$, the reciprocal lift is irreducible if
$u^2-4$ is nonsquare in $\mathbb Q(u)$.  Its norm is a rational square
times
$$
\Xi(Z)=\bigl(125D^2+64Z^4\bigr)
       \bigl(125D^2+1024Z^4\bigr).
$$
If $\Xi$ were square, setting $h=5D/(8Z^2)$ would give a rational point
with positive $x$ on
$$E:\ y^2=x(x+5)(x+80),\qquad x=25h^2.$$
A complete $2$-isogeny descent gives rank $0$ and
$$E(\mathbb Q)=\{O,(0,0),(-5,0),(-80,0)\},$$
which has no positive $x$; hence $\Xi$ is nonsquare and the lift is
irreducible (`/tmp/l22_reciprocal_cube.md`, §§2–3;
`data/l22_reciprocal_cube.jsonl`, reciprocal-lift,
two-isogeny-descent and elliptic-torsion records).  The proof uses norm
nonsquareness only in the sufficient direction.

The value $a=1$ is a specialization witness, not the final cell witness:
it need not satisfy $(5\mid w)=-1$.  If
$P(a,Z,b)$ on
$\tau^\dagger=(1+2a^2)/A$, $\delta_{\tau^\dagger}=-4a^4/A$
factored in $\mathbb Q(a)[b]$, localization at its nonzero endpoints
would specialize that factorization at $a=1$, contradicting the theorem
above.  Quantitative Hilbert irreducibility then selects the actual
irreducible $a$ inside L20's odd character-admissible progression; the
bad count is $O_Z(\sqrt B\log B)$ while the progression contributes
$B/M+O(1)$ members (`THEOREMS.md`, L21d; fixed-branch specialization
source `/tmp/l22_reciprocal_cube.md`).  This is an independent proof of
the vertical input used in the chain.

### L22c. Side lemmas and their exact scope

1. **Infinity clusters (PROVED local structure; not the closure).**
   Over $\mathbb Q(Z)((1/a))$, $Z\ne0$, higher-order Newton analysis
   resolves the repeated residuals into four irreducible quadratic
   clusters.  A proper global factor is therefore locally permitted only
   in degrees $2$, $4$, or $6$; the $a=0$ place usually narrows this to a
   $4+4$ split but has explicit square exceptional loci
   (`/tmp/l22_infinity.md`, L22-I1–I5;
   `data/l22_infinity.jsonl`).  This method alone stopped there; L22a
   globally eliminates the surviving degrees.
2. **Factor tuples (PROVED criterion and no-go; CONDITIONAL prime
   values).** If a primitive reducible $G$ is written as a product of
   irreducibles, classical Schinzel H need be applied only to its
   odd-exponent factors together with the linear moving-prime
   polynomial.  An all-plus residue in the explicit resultant-character
   period is necessary and sufficient for progression refinement.
   However, the admissible tuple
   $(8t+5,8t+7,8t+23)$ has individual outside signs $(-1,-1)$ at every
   simultaneous-prime value although their product is $+1$
   (`/tmp/l22_factor_tuple.md`, §§1–4;
   `data/l22_factor_tuple.jsonl`, summary).  Reciprocity controls only
   the product: the route weakens irreducibility when compatibility is
   separately known, but otherwise relocates rather than removes the
   obligation.  It is not used in the six-count closure.
3. **Free-$\lambda$ square branch (PROVED off-chain theorem; NOT
   APPLICABLE to the six-count).** Every rational solution of
   $X^2-r^2=A$ is
   $$X=(\lambda+A/\lambda)/2,\quad
     r=(\lambda-A/\lambda)/2,\quad
     \tau=X/A,\quad\delta=-r^2/A,$$
   and the associated pencil is irreducible over
   $\mathbb Q(\lambda)[b]`; Hilbert irreducibility plus the enlarged
   controlled set constructs an aligned class for each cell
   (`/tmp/l22_square_branch.md`, §§1–5;
   `data/l22_square_branch.jsonl`).  But the selected $\lambda$ varies
   through the infinite L11c family after the cell is known.  L11c
   permits only a predeclared finite branch menu at the current witness
   count; making $\lambda$ a witness would raise the count by one
   (`/tmp/l22_square_branch.md`, Six-count warning;
   `agent://BranchAudit`).  **This theorem is not, and must never be
   cited as, the six-count closure.**
4. **Geometric Noether reduction (PROVED reduction; OPEN auxiliary
   nonvanishing lemma).** The norm-pencil model gives a descent from
   absolute irreducibility of $S_r$,
   $r=(1-Z)/Z^2$, to vertical irreducibility of $H$.  Its generic member
   is geometrically integral, and fixed-fibre reducibility is detected
   by the gcd $N(r)$ of the $483$-minors of an explicit
   $903\times483$ integer matrix (`/tmp/l22_fiber_geometry.md`, §§2–5;
   `data/l22_fiber_geometry.jsonl`).  The statement
   $N((1-z^3)/z^6)\ne0$ for every cell remains **OPEN** in that auxiliary
   route (`/tmp/l22_fiber_geometry.md`, §7).  L22a proves the required
   vertical theorem directly, so this open polynomial nonvanishing is
   not a premise of the record.

### L22d. Chain consequence

Fix a cell $(w,z)$ and put $Z=z^3\ne0$.  L22a proves
$P(a,Z,b)$ irreducible in $\mathbb Q(a)[b]$ on the fixed branch
$\tau^\dagger=(1+2a^2)/A$,
$\delta_{\tau^\dagger}=-4a^4/A$.  Quantitative Hilbert
irreducibility inside L20's nonempty odd character-admissible
$a$-progression therefore supplies a concrete admissible $a$ with an
irreducible specialization (`/tmp/l22_elimination.md`;
`THEOREMS.md`, L20 and L21d).

For that $a$, L20 takes $f=w$ and constructs an aligned class.  In the
Schinzel admissibility package, condition (a) — positive leading
coefficient and irreducibility — is now **PROVED** by L20 plus L22;
conditions (b), (c), and (d) were already **PROVED** uniformly by L20.
Thus all four conditions are **PROVED** (`/tmp/l20_admissible.md`;
`/tmp/l22_elimination.md`; `data/l22_elimination.jsonl`).

Classical Schinzel H for the resulting pair
$\{q_1+Nt,G(t)\}$ then supplies infinitely many emergent-free members by
L19, so the intermediate per-cell hypothesis H holds.  The L6 assembly
and the unchanged finite branch menu give Theorem C with exactly six
universal quantifiers (`CONDITIONAL.md`, §§1–2; `THEOREMS.md`, L6,
L19–L20).  Therefore
$$
\boxed{\text{classical Schinzel H}\ \Longrightarrow\
       \text{intermediate H}\ \Longrightarrow\ \text{Theorem C}.}
$$
This is a **PROVED implication with CONDITIONAL conclusion**.  Classical
Schinzel H is unproved; neither it nor Hilbert's Tenth Problem over
$\mathbb Q$ is claimed solved.

## L23 [UNCONDITIONAL FRONTIER SHARPENED, 2026-08-22]. A genuine half-sieve, and the exact algebraic and large-divisor barriers

Throughout L23, work on the fixed canonical branch
$$
\tau^\dagger=\frac{1+2a^2}{A}=\frac{A+1}{2A},\qquad
\delta_{\tau^\dagger}=-\frac{4a^4}{A},\qquad
\alpha=-\delta_{\tau^\dagger}A=4a^4=(2a^2)^2,
$$
where
$$
A=1+4a^2,\quad s=\frac{a-1}{2},\quad
D=1-Z-a^2Z^2,\quad
N_g=16a^4b^2-A(b-1)^4,
$$
and
$$
\begin{aligned}
P(b)&=16D^2A^2b^4+\frac{4a^8}{A}Z^4N_g^2
      -32A^3s^2D^2b^5,\\
H&=\frac A4P
  =a^8Z^4N_g^2+4A^3D^2b^4-8A^4s^2D^2b^5.
\end{aligned}
$$
Put
$$
X=a^4Z^2N_g,\qquad Y=2ADb^2,\qquad
L=1-2As^2b.
$$
The exact identity used throughout is
$$
\boxed{H=X^2+ALY^2.}
$$
On this branch the tied quaternion has square class
$$
\boxed{(x_0,d_0)=(P(b),2b)}
$$
when $P(b)\ne0$; when $P(b)=0$, equivalently $M=0$, the original tied
equation is solved by $(y,r)=(0,0)$ and no quaternion with a zero slot is
invoked (`l23_norm_section.py`; `data/l23_norm_section.jsonl`).

### L23a. The half-dimensional bad-root sieve

**Theorem (PROVED; unconditional small-prime cleanliness).** Fix one
L20/L22 selected aligned class and write
$$
Q(t)=q_1+Nt,\qquad b(t)=\varepsilon fQ(t),
$$
with $\gcd(q_1,N)=1$.  Let $G(t)$ be the primitive irreducible
degree-$8$ part of $P(\varepsilon fQ(t))$, assume $QG$ is
fixed-divisor-free, and put
$$u(t)=2\alpha\varepsilon fQ(t).$$
If $\theta$ is a root of $G$, assume $u(\theta)$ is nonsquare in
$K=\mathbb Q(\theta)$.  Put the fixed-data, discriminant and resultant
primes, and every $p\le9$, into a finite set $S$.  Fixed-divisor-freeness
supplies a residue at each $p\in S$ on which $QG$ is nonzero; absorb all
of those finitely many choices into the progression modulus.

For $p\notin S$, define
$$
r_-(p)=\#\{r\bmod p:\ G(r)=0,\ 
                 \bigl(\tfrac{u(r)}p\bigr)=-1\}.
$$
A simple root has Haar odd-valuation mass $1/(p+1)$, so
$$
m_p=\frac{r_-(p)}{p+1},\qquad
g_{\rm odd}(p)=\frac{p\,r_-(p)}{p^2-1},\qquad
g_{\rm root}(p)=\frac{r_-(p)}{p-1}.
$$
The mod-$p$ root condition deliberately oversieves values with even
positive valuation.  Its excess
$$g_{\rm root}(p)-g_{\rm odd}(p)=\frac{r_-(p)}{p^2-1}$$
is summable.  Since $K(\sqrt{u(\theta)})/K$ is a nontrivial quadratic
extension, Chebotarev gives
$$
\sum_{p<z}g_{\rm root}(p)\log p=\frac12\log z+O(1),\qquad
\prod_{p<z}(1-g_{\rm root}(p))
   \sim\frac{C}{(\log z)^{1/2}}.
$$
Thus the sieve dimension is exactly
$$\boxed{\kappa=\frac12.}$$

For every squarefree sifting modulus $d$, the bad-root conditions are a
union of at most $8^{\omega(d)}$ reduced classes for $Q(t)$ modulo $Nd$.
The resulting residue count is divisor-bounded.  Divisor-weighted
Bombieri--Vinogradov therefore supplies level
$$
D_{\rm BV}=\frac{X^{1/2}}{(\log X)^B}.
$$
The dimension-$1/2$ beta sieve is the semilinear sieve, whose lower
sifting limit is $\beta(1/2)=1$.  Taking
$$
z=X^{49/100},\qquad
\frac{\log D_{\rm BV}}{\log z}\longrightarrow\frac{50}{49}>1
$$
gives
$$
\boxed{\#\{X<t\le2X:\ Q(t)\ {\rm prime},\
G(t)\ {\rm has\ no\ bad\ root\ prime}\ p<z\}
\gg\frac{X}{(\log X)^{3/2}}.}
$$
Any fixed $z=X^{1/2-\epsilon}$ works.  The implied constant is
ineffective through Chebotarev/Siegel.

**Strict scope.** This is an all-$p<z$ theorem, not an all-prime member
theorem.  The mod-$p$ oversieve keeps every sifting modulus squarefree;
it does not distinguish $v_p(G(t))=1$ from $v_p(G(t))\ge2$.  Encoding
exact odd valuation would introduce $p^2$ and reduce the underlying
product level to $X^{1/4-o(1)}$.  No finite local replay is used as the
analytic proof (`l23_half_sieve.py`; `data/l23_half_sieve.jsonl`;
`agent://HalfSieveAudit`).

### L23b. The first missing estimate is beyond Bombieri--Vinogradov

**Barrier theorem (PROVED).** For $t\asymp X$,
$$|G(t)|=X^{8+o(1)}.$$
After sifting to $z=X^{49/100}$, a value can retain as many as
$$\left\lfloor\frac8{49/100}\right\rfloor=16$$
prime factors at least $z$.  Hilbert reciprocity says only that the
number of bad odd-valuation primes is one of
$$0,2,\ldots,16;$$
it does not separate $0$ from $2$.  Parity alone would force zero only
for $z>|G(t)|^{1/2}=X^{4+o(1)}$, which at sifting limit $1$ would require
$D>X^{4+o(1)}$.  Bombieri--Vinogradov gives only
$X^{1/2-o(1)}$, and even Elliott--Halberstam at $X^{1-o(1)}$ leaves up
to eight large factors.

Every failing small-prime-clean member contains bad primes
$p_1,p_2\ge z$ with $p_1p_2\ge X^{0.98}>D_{\rm BV}$.  The first
unavailable analytic term is
$$
\sum_{\substack{p_1,p_2\ge z\\p_1p_2>D_{\rm BV}}}
\#\left\{t\asymp X:
\begin{array}{l}
Q(t)\ {\rm prime},\quad p_1p_2\mid G(t),\\
p_1,p_2\ {\rm have\ bad\ signs\ and\ odd\ valuations}
\end{array}\right\}.
$$
This parity-sensitive two-large-bad-divisor sector lies beyond both BV
and the beta-sieve fundamental lemma.  L31d sharpens the requested
estimate: the unsigned root-class expansion has a positive main term,
so a power-saving bound for the whole sector is stronger than the
natural expectation.  What is actually needed is its evaluation inside
a fixed-family Buchstab/Hilbert-detector asymptotic with positive
zero-bad main term.

Equivalently, one needs an $R_{\rm bad}\le1$ theorem for
sign-decorated divisors of an octic at prime arguments.  The threshold is
exact: $R=2$ permits two distinct bad odd-valuation primes whose symbols
multiply to $+1$ but obstruct at both places.  Kao's Theorem 1/Table 1
gives only $P_{12}$ for an irreducible octic at prime arguments
(improving Irving's $P_{14}$), and even its use here requires a fixed-AP
adaptation.  Therefore L23a does **not** remove classical Schinzel H.

### L23c. Capell's square-in-the-octic-field route is unavailable on the selected protocol

Let $P$ be irreducible at the chosen rational $(a,Z)$, let
$K=\mathbb Q(\theta)$ with $P(\theta)=0$, and normalize $P$ to be monic.
Capell's criterion gives the exact equivalence
$$
\boxed{2\theta\in K^{\times2}
\iff P(u^2/2)\ \text{is reducible over }\mathbb Q,}
$$
while
$$
N_{K/\mathbb Q}(2\theta)
=2^8\frac{P(0)}{\operatorname{lc}(P)}=256
$$
is necessary only and is not a square certificate.

Let $w$ be the target prime, $e=v_w(z)\ge1$, and $Z=z^3$.  The left
Newton edge joins $(0,12e)$ to $(4,0)$.  If $e$ is odd, its separable
quartic residual gives an unramified place with
$v_w(\theta)=3e$ odd, so $2\theta$ is nonsquare for **every**
admissible $a$.  If $w\nmid s$, reduction modulo $w$ gives
$$
P(b)\equiv16A^2b^4(1-2As^2b)\pmod w
$$
and the horizontal side has the simple root
$$
\beta=\frac1{2As^2},\qquad
2\beta=\frac1{As^2},\qquad
\left(\frac{2\beta}{w}\right)=\left(\frac Aw\right)=-1.
$$
Thus the untwisted Capell collapse is impossible for every $e\ge1$ on
that refined stratum.

The L20 residue choice can be refined freely.  The exact count is
$$
\#\left\{a\bmod w:
\left(\frac{1+4a^2}{w}\right)=-1\right\}
=\frac{w-\left(\frac{-1}{w}\right)}2
=
\begin{cases}
(w+1)/2,&w\equiv3\pmod4,\\
(w-1)/2,&w\equiv1\pmod4.
\end{cases}
$$
It is at least $2$ for every odd prime, so one may delete
$a\equiv1\pmod w$.  CRT, the L20 symbol calculation, and the L21d
$O_Z(\sqrt B\log B)$ thin-set bound remain valid in the refined
progression.  Hence the selected L20/L22 protocol has $2\theta$
nonsquare and the conic bundle below has rank $10$ on every cell.

**Strict scopes.** For a pre-existing class with even $v_w(z)$ and
$w\mid s$, the uniform all-$a$ local statement remains **OPEN**; it is
not inferred from a scan.  The recorded $353$ L20 grid rows are
nevertheless settled exactly: $348$ by the odd left edge, $2$ by the
horizontal residue, and the $3$ remaining $a=1$ rows by the reciprocal
trace-norm lemma below, with zero refusals.  No universal $2$-adic
Capell no-go is claimed: the first residual can be inseparable and the
$b^5$ term changes the hull.  On the reciprocal slice $a=1$, for every
rational $Z\ne0$,
$$
P(b)=b^4T(b+b^{-1}),\qquad
T(v)=400D^2+\frac45Z^4(16-5(v-2)^2)^2.
$$
If $2\theta$ were square, one of the two necessary norm squareclasses
would be
$$125D^2+1024Z^4\quad\text{or}\quad125D^2+64Z^4.$$
For reduced $Z=p/q$ and $d=q^2-pq-p^2$, $d$ is odd and both cleared
quantities
$$125d^2+1024p^4,\qquad125d^2+64p^4$$
are $5\bmod8$, hence nonsquares
(`l23_half_sieve.py`, `l23_fibration.py`;
`data/l23_{half_sieve,fibration}.jsonl`).

### L23d. The fixed-cell conic bundle has non-split rank ten

**Premise.** Fix a rational $(a,Z)$ for which $P_{a,Z}$ is irreducible
of degree $8$.  L22 proves irreducibility in $\mathbb Q(a)[b]$; the
actual rational $a$ used here is supplied by the L21d quantitative-HIT
selection in the L20 progression.  Under that premise the member
problem is the rational-point problem on
$$
\boxed{X_{a,Z}:\ U^2-2bV^2=P_{a,Z}(b)W^2.}
$$

The determinant squareclass is $2bP(b)$.  The natural regular conic
bundle therefore has ten geometric degenerate fibres, grouped into
closed degrees
$$1+8+1,$$
with component/splitting squareclasses
$$
A,\qquad2\theta\in K^\times/K^{\times2},\qquad A.
$$
Indeed
$$
P_0=P_8=4Aa^8Z^4=A(2a^4Z^2)^2,\qquad
N_{K/\mathbb Q}(2\theta)=256.
$$
The Faddeev residue calculation has the sole relation
$e_0=e_\infty$.  Consequently
$$
\boxed{
\begin{array}{c|c|c}
&\text{non-split rank}&\operatorname{Br}(X)/\operatorname{Br}(\mathbb Q)\\
\hline
2\theta\notin K^{\times2}&10&\mathbb Z/2,\ \text{generated by }(A,b)\\
2\theta\in K^{\times2}&2&0
\end{array}}
$$
and $K_X^2=8-10=-2$.  Here rank is the degree of the non-split closed
locus, not the number of geometric singular fibres.

By L23c the selected protocol is in the first row.  The exact
few-fibre theorems do not decide it: Harpaz--Wei--Wittenberg Theorem
1.3 first fails at rank $\le2$ (and its rank-$3$ clauses); their
Theorem 1.4 retains the homogeneous Schinzel hypothesis at the
degree-$8$ point; Browning--Schindler first fails at rank $\le3$;
Shute first fails at residue-field degree at most $3$; and the
Harpaz--Skorobogatov--Wittenberg route does not cover the non-split
degree-$8$ degenerate fibre over a non-rational closed point.  Classical
Chatelet theorems also do not apply: the quadratic algebra
$\mathbb Q(b)(\sqrt{2b})$ moves with $b$ and there are ten bad fibres.
No checked unconditional fibration theorem closes this rank-$10$
surface (`l23_fibration.py`; `data/l23_fibration.jsonl`;
`agent://FibrationAudit`).

### L23e. Exact norm-section and squareclass no-go theorems

**Exact matching (PROVED empty on $\Phi$).** Put $d=2b$.  In the complete
ansatz where $U,V$ are coefficient-linear in $X,Y$, exact source/target
matching is
$$d=-ALq^2,$$
and the residual scalar condition is
$$
(A,d)=(A,-AL)=(A,L)=1.
$$
Writing $L=R^2-AS^2$, the resulting identity is
$$
U=\frac{2(SX-LY)}R,\qquad
V=\frac{2(X+ASY)}{AqR},\qquad
\boxed{P=U^2-dV^2}.
$$
The two rational charts
$$
b_-=-\frac{\vartheta^2}{2As^2(1-\vartheta^2)},\quad
q=\frac{\vartheta}{As},\qquad
b_+=\frac1{2As^2(1-\vartheta^2)},\quad
q=\frac1{As\vartheta}
$$
are exactly the L8c Pell gauge.  Their complete $2$-adic valuation
analysis, including $s=0$ and the endpoints, gives $v_2(b)\ne0$ at
every point.  Thus exact matching, even with its residual scalar
condition restored, is $\Phi$-empty.

**No polynomial or generic $b$-line section (PROVED).** If
$$P(b)=U(b)^2-2bV(b)^2,\qquad U,V\in\mathbb Q(a,Z)[b],$$
then at $b=0$
$$
U(0)^2=P(0)=4a^8AZ^4=(2a^4Z^2)^2A,
$$
impossible because $A=1+4a^2$ has odd $A$-valuation.  More generally,
the endpoint unit of the normalized generic conic at both $b=0$ and
$b=\infty$ has squareclass $1/A$, so the conic has no
$\mathbb Q(a,Z)(b)$-point.  This rules out Laurent, polynomial and
rational sections over the generic $b$-line, but not specialized
rational points or nonlinear sections outside the classified ansatz.
The target-split section $2b=t^2$ exists algebraically and has
$v_2(b)=2v_2(t)-1$, so it is also $\Phi$-empty.

**Dyadic squareclass eliminations (PROVED).** On $\Phi$, $a$ and $b$
are $2$-adic units.  If $m=v_2(Z)$, exact dominance in
$H=X^2+ALY^2$ gives
$$
\boxed{v_2(P)=4+4\min(m,0)\equiv0\pmod4.}
$$
Hence neither
$$P=2b\,y^2\qquad\text{nor}\qquad P=-2b\,y^2$$
has a $\Phi$-admissible rational point.  If instead
$P=A\,N_{\mathbb Q(\sqrt{2b})/\mathbb Q}(\xi)$, then
$$
(P,2b)_2=(A,2b)_2=-1,
$$
so this norm condition forces nonsplitting rather than splitting.

For fixed $(a,b)$,
$$
P(Z)=C(1-Z-a^2Z^2)^2+EZ^4,\quad
C=16A^2b^4L,\quad E=\frac{4a^8}{A}N_g^2,
$$
and
$$
\operatorname{Disc}_Z(P)=16C^3E^2(A^2C+16E).
$$
On $\Phi$ every displayed factor is nonzero.  Thus the fixed-$P$ square
locus $y^2=P(Z)$ is a smooth genus-$1$ curve, never the hoped-for
genus-$0$ degeneration.

**First $\Phi$-admissible factor ansatz (PROVED reduction; OPEN points).**
The condition
$$
L=-A\rho^2,\qquad
b=\frac{1+A\rho^2}{2As^2}
$$
is $\Phi$-admissible exactly when $\rho$ and $s$ are $2$-adic units
(apart from the ordinary guards), and it gives the genuine factorization
$$
H=(X-A\rho Y)(X+A\rho Y).
$$
It does not make the factors norms or force all Hilbert signs to be
$+1$.  With $T=1+A\rho^2$,
$U_0=2T-A(a-1)^2$, and
$$
V(\rho)=64a^8A(a-1)^4T^2
 \pm2A^4(a-1)^4\rho T^2-a^4U_0^4,
$$
the zero-factor route is birational to
$$
\omega^2=\pm2A\rho V(\rho).
$$
The branch polynomial $\rho V(\rho)$ is generically squarefree of
degree $9$; its smooth completion has genus $4$.  This proves that the
ansatz is not a rational parametrization.  Isolated rational points for
a fixed cell remain **OPEN**, and any such point must still satisfy the
aligned-class congruences separately
(`l23_norm_section.py`, `l23_rational_section.py`,
`l23_squareclass.py`; corresponding `data/l23_*.jsonl`).

**Fixed twists (PROVED obstruction, not a closure).** For fixed
$j\in\mathbb Q^\times$ and monic normalization $p$,
$$
j\,2\theta\in K^{\times2}
\iff p(u^2/(2j))\ \text{factors}.
$$
At a simple value-prime this changes the required sign to
$(j\mid\ell)$; it relocates the character instead of making it $+1$.
If one simultaneously imposes a polynomial norm identity
$P=c(E^2-jO^2)$, evaluation at $\theta$ puts $\sqrt j$ in $K$, so the
twisted criterion implies the forbidden untwisted one.  A single
rational $j$ also cannot have odd valuation at infinitely many target
primes.  No fixed twist removes the member problem.

### L23f. Moving $a$ does not lower the certified squareclass degree in the tested families

When $a=a(Q)$ is nonconstant, $P=4H/A$ and the same polynomial
squareclass is represented by
$$G_{\rm move}=AH.$$
This is a complexity test, not a reuse of L20's fixed aligned class,
whose modulus freezes $a$.  For $b=\kappa Q$ and $\deg a=r\ge1$,
exact leading-term comparison gives
$$
\begin{array}{c|c}
\text{case}&\deg(AH)\\
\hline
r\ge2&18r+4\\
r=1,\ 4c^2\ne\kappa^2&22\\
r=1,\ 4c^2=\kappa^2&21 .
\end{array}
$$
Those are raw degrees only.

The missing $a(0)=0$ branch corrects the earlier apparent minimum:
$$Q^4\Vert H,\qquad A\,H/Q^4$$
represents the same squareclass.  On the exact diagonal
$$
Q(t)=7+12t,\qquad a(Q)=Q,\qquad b(Q)=Q,\qquad Z=27,
$$
the reduced representative has degree $18$ and is squarefree with
irreducible factor degrees $2$ and $16$ (the degree-$16$ factor is
irreducible modulo $43$).  On the separate
$a(0)\ne0$ diagonal $a=(Q+3)/2$, the raw degree is $21$ with factor
degrees $2$ and $19$.  These are certified examples and obstructions
for the tested polynomial diagonals, **not** a universal minimum over
all rational substitutions; in particular there is no minimum-$21$
claim (`l23_multivar.py`; `data/l23_multivar.jsonl`).

### L23g. Exact unconditional frontier

L23 proves that small-prime cleanliness is abundant and closes the
untwisted Capell, fixed-twist, direct norm-section, elementary
squareclass and tested moving-parameter escape routes at their stated
scopes.  It proves **no globally good member theorem**.  Classical
Schinzel H therefore remains the sole conjectural input in L22d.  The
first analytic target is the parity-sensitive two-large-bad-divisor
estimate in L23b; the first algebraic target is a genuinely
non-diagonal family outside the dyadic and norm walls.

## L24 [DIAGONAL ABSORPTION CLOSED, 2026-08-22]. Both exact degree-eight images are empty on $\Phi$ over $\mathbb Q_2$

### L24a. The self-coupled diagonal identities and the nominal five-count

Put
$$
a=1+2s,\quad A=1+4a^2,\quad B=2b,\quad Z=z^3,\quad
D=1-Z-a^2Z^2,
$$
$$
N_g=16a^4b^2-A(b-1)^4,\qquad
c=\frac{a^2Z^2N_g}{Ab^2D}.
$$
On the square branch write
$$
X=\frac{\lambda+A/\lambda}{2},\qquad
\rho=\frac{\lambda-A/\lambda}{2},\qquad
X^2-\rho^2=A,
$$
so $\tau=X/A$ and $\delta=-\rho^2/A$.  The exact tied equation is
$$
\delta(c^2-Ay^2)-16B(r^2-As^2)=16.
$$

There are two self-couplings.

**Orientation I** ($X=y,\rho=r$):
$$
\begin{aligned}
C_1&=y^2-r^2-A=0,\\
C_2&=-r^2(c^2-Ay^2)-16AB(r^2-As^2)-16A=0.
\end{aligned}
$$
With $u=r^2$ this eliminates to
$$
\boxed{Q_I(u)=
Au^2+(A^2-c^2-16AB)u+16A^2Bs^2-16A=0,}
$$
and $y^2=A+u$.

**Orientation II** ($X=r,\rho=y$):
$$
\begin{aligned}
C'_1&=r^2-y^2-A=0,\\
C'_2&=-y^2(c^2-Ay^2)-16AB(r^2-As^2)-16A=0.
\end{aligned}
$$
With $u=y^2$ this eliminates to
$$
\boxed{Q_{II}(u)=
Au^2-(c^2+16AB)u+16A^2B(s^2-1)-16A=0,}
$$
and $r^2=A+u$.

These are polynomial finite covers after clearing only a square.  Put
$$
E=Ab^2D,\qquad K=a^2Z^2N_g,\qquad T=AE^2.
$$
Then $c=K/E$ and the cleared equations are
$$
\begin{aligned}
F_I={}&Tu^2+\bigl((A-32b)T-K^2\bigr)u
       +16(2Abs^2-1)T=0,\\
F_{II}={}&Tu^2-(K^2+32bT)u
       +16(2Ab(s^2-1)-1)T=0.
\end{aligned}
$$
The only occurrence of $c$ after clearing is
$c^2E^2=K^2$; no inverse witness or inequation variable is introduced.
The guards are genuine: on $\Phi$, $A,b,D,N_g,Z$ and
$\lambda=X+\rho$ are nonzero.

For either orientation write the normalized $u$-quadratic as
$$u^2+(L_0-c^2/A)u+M_0=0,$$
where
$$
(L_0,M_0)=
\begin{cases}
(A-16B,\ 16(ABs^2-1)),&I,\\
(-16B,\ -16AB(1-s^2)-16),&II.
\end{cases}
$$
Substituting
$u=(\lambda^2-A)^2/(4\lambda^2)$ gives the even
$A$-reciprocal octic
$$
\mathcal G(\lambda)=
(\lambda^2-A)^4+
4(L_0-c^2/A)\lambda^2(\lambda^2-A)^2+
16M_0\lambda^4.
$$
After localizing at $A$, the coordinate algebra is free of rank
$4\cdot2=8$: both complete-intersection covers are finite flat of
degree $8$, including over the branch divisor.

Before local admissibility is imposed, DDF therefore bounds the tied
image rank by $1$, and the formal count would be
$2+(3+1-1)\le5$.  This algebraic rank calculation is PROVED.  L24b
shows that its image on $\Phi$ is empty, so the nominal count is
vacuous.

### L24b. Final dyadic theorem: both orientations are uniformly empty

**Theorem (PROVED).** The rational image of each diagonal cover is
empty on $\Phi$ over $\mathbb Q_2$.

On $\Phi$,
$$
s\in\mathbb Z_2,\qquad b\in\mathbb Z_2^\times,\qquad
a\in1+2\mathbb Z_2,\qquad A\equiv5\pmod{32}.
$$
Since $b-1$ is even, $v_2(N_g)\ge4$.  Moreover
$$
v_2(D)=
\begin{cases}
0,&v_2(Z)\ge0,\\
2v_2(Z),&v_2(Z)<0,
\end{cases}
$$
so
$$\boxed{v_2(c)\ge4}$$
for every guarded rational $Z$.  The guard is automatic on this locus:
$D=0$ would make $A$ the rational square discriminant of
$a^2Z^2+Z-1$.

For orientation I,
$$
\frac{Q_I(u)}A
=u^2+\left(A-32b-\frac{c^2}A\right)u
16(2Abs^2-1).
$$
Its coefficient valuations, from constant through quadratic, are
exactly $(4,0,0)$.  If $t=v_2(u)$, the term valuations
$(2t,t,4)$ force $t\in\{0,4\}$.  A square unit has
$u\equiv1\pmod8$, so $A+u\equiv6\pmod8$ has odd valuation; if
$t=4$, then $A+u\equiv5\pmod8$.  Neither is a square.  Thus
orientation I has no $\mathbb Q_2$-point over any $\Phi$ base.

For orientation II, write
$$
\frac{Q_{II}(u)}A=u^2-32hu+16e,
$$
where
$$
h=b+\frac{c^2}{32A}\in\mathbb Z_2^\times,\qquad
e=2Ab(s^2-1)-1\in\mathbb Z_2^\times.
$$
The coefficient valuations are $(4,5,0)$, so the unique Newton slope
forces $v_2(u)=2$.  If the required square existed, write
$u=4t^2$ with $t$ odd.  Division by $16$ and reduction modulo $16$
would give
$$
0\equiv t^4-8ht^2+e
\equiv8+2Ab(s^2-1)\pmod{16}.
$$
This would require $v_2(s^2-1)=2$, impossible: that valuation is $0$
for even $s$ and at least $3$ (or infinite) for odd $s$.  This single
congruence covers both parities of $s$.  Hence orientation II also has
no $\mathbb Q_2$-point over any $\Phi$ base.

Therefore
$$
\boxed{\operatorname{im}(I)(\mathbb Q)\cap\Phi
=\operatorname{im}(II)(\mathbb Q)\cap\Phi=\varnothing.}
$$
The **diagonal** five-unknown route is PROVED FALSE, not open, and its
$\le5$ count is vacuous.  This conclusion does not apply to the fixed
canonical $\Theta^\ast$ architecture, does not change
$\operatorname{efd}_{\mathbb Q}\le5$ under classical Schinzel H, and
does not alter the conditional six-count.

### L24c. Odd target places are not what kills the diagonal

Let $t=v_w(Z)>0$ and $k=v_w(b)$ at an odd target prime.  On the
standard stratum $0<k<t$ the reductions are
$$
u(A+u)=16\quad(I),\qquad u^2=16\quad(II).
$$
Orientation II has a nonsingular total-space $w$-adic point for every
odd $w$ on the standard $k=1$ stratum.  For $w\ge5$ one can choose a
fibre-regular local point while keeping that chosen local base fixed;
this local theorem does not assert compatibility with an independently
preselected W0/W2 aligned residue class.  At $w=3$ the fibre Jacobian
drops rank on the $u=4$ arm, but
$$\frac{\partial C_1}{\partial s}=-16a,\qquad
\frac{\partial C_2}{\partial y}=4Ay^3$$
restores full total-space rank.  Separately, the critical stratum
$k=t$ is empty at $w=5$, and the pole stratum $k>t$ is empty in both
orientations.  These are exact local classifications.  They do not
weaken L24b: the obstruction closing the rational route is genuinely
dyadic (`l24_diagonal_local.py`; `data/l24_diagonal_local.jsonl`).

### L24d. Geometry and the bounded corroborating search

The guarded generic cover is integral of degree $8$.  Its Galois closure
has group
$$V_4\wr C_2,\qquad |V_4\wr C_2|=32,$$
with three irreducible affine branch components.  It has no
$\mathbb Q(s,b,Z)$-point and hence no base-rational section.  On the
actual-cell slice $(s,Z)=(2,27)$, exact ramification
$$16\cdot4+1\cdot2+8\cdot2+1\cdot2=84$$
gives genus $35$.  The independent-$c$ reduction is genus $3$, and
the actual $z$-pullback on the recorded $(s,b)=(1,3)$ slice has genus
$53$.  Absolute, non-base-preserving rationality or unirationality of
the two-dimensional total surfaces remains **OPEN**; it is not inferred
from the high-genus fibres
(`l24_diagonal_geometry.py`, `l24_diagonal_arithmetic.py`;
corresponding JSONL artifacts).

Finally, the exact bounded search covered
$$
370\ \text{cells},\quad2{,}013{,}141\ \text{guarded base pairs},\quad
4{,}026{,}282\ \text{orientation attempts}
$$
and found zero hits, with $2{,}296$ fixed-base local certificates
(`l24_diagonal_search.py`; `data/l24_diagonal_search.jsonl`).  This
zero-hit scan is **CORROBORATION ONLY**.  Every-cell emptiness is proved
by the mod-$16$ theorem in L24b, never by finite exhaustion.

### L24e. Chain consequence

The L22 implication is unchanged:
$$
\boxed{\text{classical Schinzel H}
\Longrightarrow\text{intermediate H}
\Longrightarrow\text{Theorem C}.}
$$
L23 and L24 sharpen the unconditional frontier and close tempting
escape routes; they do not remove Schinzel H, improve the conditional
record, or solve H10/$\mathbb Q$.  The exact next targets are the
two-large-bad-divisor estimate of L23b or a genuinely non-diagonal
algebraic family escaping the dyadic and norm walls.

## L25 [SCALED SELF-COUPLINGS, 2026-08-23]. The dyadic coupling wall is artifactual, not universal

L24 closed the two unscaled diagonal orientations.  It thereby proved
the **unit coupling** route false, not every non-diagonal coupling.
The next exact step is to scale the tied norm witnesses by fixed rational
numbers depending on no new free variable.  Work on the L11c square
branch
$$
X=\frac{\lambda+A/\lambda}{2},\qquad
\rho=\frac{\lambda-A/\lambda}{2},\qquad X^2-\rho^2=A,
$$
so that $\delta=-\rho^2/A$, and use the L6 variables
$s=(a-1)/2$, $B=2b$, and $c=h(a,b,Z)$, $Z=z^3$.  On $\Phi$ we have
$$
s\in\mathbb Z_2,\qquad b\in\mathbb Z_2^\times,\qquad
a=1+2s,\qquad A=1+4a^2\equiv5\pmod{32},\qquad v_2(c)\ge4.
$$

### L25a. Exactly two fixed scalings are needed

For fixed rational $\kappa_1,\kappa_2$, the type-I coupling
$(y,r)=(\kappa_1X,\kappa_2\rho)$ and the type-II coupling
$(y,r)=(\kappa_1\rho,\kappa_2X)$ introduce no new free variable: the
witness pair is a fixed linear transform of the two branch coordinates.
Soundness is inherited from L11c and L6: substituting a rational
branch and rational witness pair into the tied equation satisfies the
same polynomial soundness proof as the L24 specializations.  With
$u=\rho^2$, exact elimination gives
$$
\begin{aligned}
\text{I}:&\quad
A\kappa_1^2u^2+
 \bigl(A^2\kappa_1^2-c^2-16AB\kappa_2^2\bigr)u
 +16A^2Bs^2-16A=0,\\
\text{II}:&\quad
A\kappa_1^2u^2-
 \bigl(c^2+16AB\kappa_2^2\bigr)u
 +16A^2B(s^2-\kappa_2^2)-16A=0.
\end{aligned}
$$
These identities were replayed on $1200$ exact instances
(`l25_scaled_coupling.py::identity_replay`).

**Theorem (PROVED).** On $\Phi$ over $\mathbb Q_2$:

1. The unscaled type-I coupling $(X,\rho)$ has no admissible root:
   L24's orientation-I argument applies, with coefficient valuations
   $(4,0,0)$ and impossible simultaneous squareness of $u$ and $A+u$.
2. The scaled coupling
   $$(y,r)=(2X,\rho)$$
   is $2$-adically admissible exactly for **even** $s$, uniformly in
   $b$, in the unit classes of the fixed scalars, and in every $c$ with
   $v_2(c)\ge4$.
3. The complementary scaled coupling
   $$(y,r)=(2X,\rho/2)$$
   is $2$-adically admissible exactly for **odd** $s$, with the same
   uniformity.
4. Consequently the disjunction of these two fixed formulas is
   $2$-adically nonempty on **every** $\Phi$ stratum.  L24b's dyadic
   wall is therefore not a property of the coupled architecture but of
   the two unscaled formulas.

**Proof.** For $(\kappa_1,\kappa_2)=(2\mu,\nu)$ with $\mu$ a $2$-adic
unit and $\nu\in\{1,1/2\}$ times a $2$-adic unit, the type-I equation
is
$$
4A\mu^2u^2+
\bigl(4A^2\mu^2-16AB\nu^2-c^2\bigr)u+32A^2bs^2-16A=0.
$$
Divide by $4A$ and put $u=4t^2$.  If a $2$-adic unit $t$ exists, then
$u$ and
$$A+u=A+4t^2\equiv A+4\equiv1\pmod8$$
are both $2$-adic squares.  The equation for $t$ is
$$
g(t)=4\mu^2t^4+Lt^2+M_0=0,
$$
where
$$
L=A\mu^2-4B\nu^2-\frac{c^2}{4A},\qquad M_0=2Abs^2-1.
$$
Both $L$ and $M_0$ are $2$-adic units.  For odd $t$,
$$g'(t)=16\mu^2t^3+2Lt$$
has valuation exactly $1$.  By the multivariate/Hensel criterion, a
root modulo $8$ therefore lifts to a $2$-adic unit root.  Since
$\mu^2\equiv1\bmod8$, $t^2\equiv1\bmod8$, and $v_2(c^2/(4A))\ge6$,
$$
g(1)\equiv4+L+2Abs^2-1\equiv3+L+2Abs^2\pmod8.
$$
For $\nu=1$, $L\equiv A\bmod8$, hence
$$g(1)\equiv3+A+2Abs^2\equiv2Abs^2\pmod8,$$
so the lift exists exactly when $s$ is even.  For $\nu=1/2$,
$L\equiv A-2b\bmod8$, hence
$$g(1)\equiv3+A-2b+2Abs^2\equiv2b(As^2-1)\pmod8,$$
and this is $0\bmod8$ exactly when $s$ is odd (then $s^2\equiv1\bmod8$
and $A\equiv5\bmod8$).  The normalized constant term before the
$u=4t^2$ substitution is $8Abs^2-4$, of valuation $2$; the only other
possible root valuation is $0$, and a unit root has $u\equiv1\bmod8$
and $A+u\equiv6\bmod8$, never a square.  Thus the quartic condition is
not only sufficient but necessary.  The complete residue assertion is
machine-replayed on all $s\bmod8$ classes and the declared scalar,
$b$, and $c$ residue/valuation classes
(`l25_scaled_coupling.py::complete_dyadic_classification`;
`data/l25_scaled_coupling.jsonl`).  $\square$

### L25b. The remaining target and real conditions

At an odd target prime $w$, the unramified-cell reduction has
$c\equiv0$ by $v_w(c)\ge1$.  For the two canonical scalings, choose the
freely chosen unit-$b$ local stratum
$$
s\equiv0,\qquad A\equiv5,\qquad X\equiv3,\qquad \rho\equiv2,\qquad
b\equiv\kappa^{-2}\pmod w,\quad \kappa\in\{1,1/2\}.
$$
Then both cover equations hold modulo $w$, and the Jacobian determinant
of the $(\rho,X)$ block is the nonzero integer
$$
-16AX\rho\,(X^2+\rho^2-4B\kappa^2)=-2400.
$$
Reduction modulo any prime not dividing $2400$, hence every odd $w\ge7$,
gives a smooth point and Hensel's lemma supplies a $w$-adic point.  For
$w\in\{3,5\}$ the complete bounded residue scans in
`data/l25_scaled_coupling.jsonl` find points but no smooth point, so
the streamlined residue proof does not lift there; higher-level
valuative analysis is open.  This is **not** the L24 aligned
low-positive target stratum and is not compatibility with a preselected
L20 residue class
(`l25_scaled_coupling.py::target_residue_certificate`;
`universal_target_formula_certificate`).

Over $\mathbb R$, two guarded bridge samples already give a positive
real root, one per parity:
$$
(a,b,Z)=(1,3,-10)\quad\text{and}\quad(3,4,-2/5).
$$
In both samples the eliminated quadratic has positive leading
coefficient, negative linear coefficient and nonnegative discriminant,
hence a positive real root $u$; then $u$ and $A+u$ are positive reals.
These are named certificates, not an all-parameter real theorem
(`l25_scaled_coupling.py::real_viability_certificate`).

### L25c. What the scaling does **not** prove

The local escape is not a global point theorem.  A rational coupled
branch requires the quadratic in $u$ to have a rational root and the
two values $u,A+u$ simultaneously rational squares.  The bounded search
in `data/l25_scaled_coupling.jsonl` found no such point in its declared
range and is EVIDENCE ONLY.  No rational-parametrization theorem,
uniform base-congruence class, global solubility of the residual genus
curve, improved conditional count, or unconditional member is claimed.
L24 remains true verbatim, and classical Schinzel H remains the sole
conjectural input in L22d.

**Exact frontier after L25.** The one-line obstruction that killed the
unscaled diagonal is transformed into an explicit linear congruence on
$s$ by fixed scalar changes.  Closing a five-shape route would now need
a separate argument governing the parity of a freely-free base parameter
$s$ **and** global rationality of the paired square conditions on the
scaled square-branch cover.  Neither follows from Schinzel H alone as
currently phrased; conversely, no displayed solvable equation contains an
unremoved $2$-adic coupling obstruction.

## L26 [RECIPROCAL EVEN-PULLBACK TIE, 2026-08-23]. The octic acquires a quartic trace model

This layer changes the tie rather than the prime-value argument.  Apply the
Daans--Sun bridge to
$$
Z=z^2,\qquad
D=1-Z-a^2Z^2,
$$
and retain
$$
A=1+4a^2,\qquad B=2b,\qquad
N_g=16a^4b^2-A(b-1)^4,\qquad
c=\frac{a^2Z^2N_g}{Ab^2D}.
$$
In the canonical bridge domain $b\ne0,1$, replace the third
$\Psi_\tau$ witness by the rational tie
$$
\boxed{\eta=\frac{Z(b+1)}{Db}.}
$$

### L26a. Guards, soundness and formal count

The pullback is exact:
$$
v_w(z)>0\iff v_w(z^2)>0
$$
at every finite place, so the bridge target set is unchanged.  On
$\Phi_1^{\{2\}}$, $a$ and $b$ are $2$-adic units and
$A\equiv5\bmod8$.  If $D=0$, then rational $Z$ would solve
$a^2Z^2+Z-1=0$ and force its discriminant $A$ to be a rational square,
contrary to its $2$-adic squareclass.  Thus the new tie adds no
denominator guard beyond the canonical bridge guards $b\ne0,1$.

Every tied solution is literally a $\Psi_\tau(a,b,c)$ solution with
third coordinate $\eta$, so Sun's soundness and the Daans bridge apply
unchanged.  The free base remains $(a,b)$, $\Phi$ has existential rank at
most $3$, and the tied block retains only $(y,r)$.  DDF fusion therefore
has the same **formal** six-variable bound
$$
2+(3+2-1)=6.
$$
This is a count for the candidate formula, not a completeness theorem.

### L26b. Exact reciprocal identity

Put
$$
M=16-\delta_\tau c^2-32Ab\eta^2,\qquad
P_{\rm rec}=D^2A^2b^4M.
$$
Direct substitution gives
$$
\boxed{
P_{\rm rec}
=16A^2D^2b^4-\delta_\tau a^4Z^4N_g^2
 -32A^3Z^2b^3(b+1)^2.}
$$
Every term is reciprocal of weight $8$.  With
$$
u=b+b^{-1},\qquad R(u)=16a^4-A(u-2)^2,
$$
one obtains
$$
\boxed{P_{\rm rec}(b)=b^4T_\tau(u)}
$$
for the quartic trace polynomial
$$
\boxed{
T_\tau(u)=16A^2D^2-\delta_\tau a^4Z^4R(u)^2
          -32A^3Z^2(u+2).}
$$
These are polynomial identities, not factorization or member claims
(`l26_reciprocal_tie.py::identity_replay`;
`data/l26_reciprocal_tie.jsonl`).

### L26c. The dyadic and target places

On $\Phi$, write $r=v_2(Z)$.  Exact dominance gives
$$
v_2(D)=
\begin{cases}
0,&r\ge0,\\
2r,&r<0,
\end{cases}
\qquad v_2(Z/D)\ge0.
$$
Because $b+1$ is even,
$$v_2(\eta)\ge1.$$
Also $v_2(N_g)\ge4$ and hence $v_2(c)\ge4$.  On the square branch
$$
\tau^\dagger=\frac{1+2a^2}{A},\qquad
\delta=-\frac{4a^4}{A},\qquad
\alpha=-\delta A=4a^4,
$$
we have $v_2(\delta)=2$ and
$$
\frac M{16}\in1+8\mathbb Z_2.
$$
Thus $M$ and $P_{\rm rec}=D^2A^2b^4M$ are squares in
$\mathbb Q_2$.  The tied symbol is automatically $+1$ at $2$.

At a standard guarded odd target $w$, assume
$$
v_w(z)\ge1,\quad v_w(b)=1,\quad
v_w(a)=v_w(A)=v_w(\delta)=0,\quad
\left(\frac Aw\right)=-1.
$$
Then $v_w(Z)=2v_w(z)\ge2$, $D$ and $N_g$ are units, and
$$
v_w(\eta)=v_w(Z)-1\ge1,\qquad
v_w(c)=2v_w(Z)-2\ge2,\qquad M\equiv16\pmod w.
$$
Since $\alpha=4a^4$ is a square unit, the target binary conic has a
smooth residue point and Hensel lifts.  The exact row
$$
(w,a,b,Z,D,\eta,c)=
\left(3,1,3,9,-89,-\frac{12}{89},-\frac{576}{445}\right)
$$
replays the valuation pattern.  No coefficient-bad target such as
$w\mid A$ is covered by this statement.

### L26d. The bad character descends to the trace quartic

Away from $b\equiv-1\pmod p$,
$$
u+2=\frac{(b+1)^2}{b}
$$
gives the exact character identity
$$
\boxed{\left(\frac{2b}{p}\right)
=\left(\frac{2(u+2)}p\right).}
$$
Thus, after putting the finitely many primes dividing the fixed value
$T_\tau(-2)$ into the controlled set, the emergent sign is a character
of the quartic trace root.  For the canonical branch
$$
\tau^\dagger=\frac{1+2a^2}{A},\qquad
\delta=-\frac{4a^4}{A},
$$
L28 proves uniformly that $T_{\tau^\dagger}$ is irreducible, so
$K=\mathbb Q(\theta)$ has degree $4$, and that
$2(\theta+2)$ is nonsquare.  Hence
$$K\bigl(\sqrt{2(\theta+2)}\bigr)/K$$
is always the nontrivial quadratic bad-sign extension on the canonical
L26 domain.  Equivalently, Capell gives the uniformly irreducible cover
$$
T_{\tau^\dagger}(v^2/2-2).
$$

The guarded specialization $(a,z,Z)=(1,3,9)$ has reciprocal octic and
trace quartic reductions irreducible modulo $37$.  That row remains
one-fibre EVIDENCE only; L28 proves the trace field and L29 classifies
the distinct reciprocal lift uniformly by its dyadic stratum.

### L26e. Strict frontier

L26 replaces the generic octic character field by a quartic trace field
and removes the dyadic symbol on this tie.  L28 closes trace
irreducibility and bad-character nonsquareness; L29 closes the distinct
reciprocal-lift squareclass and adds a compressed pullback.  L31
subsequently proves one globally good fixed-field member at \(w=13\).
The route still does **not** produce a uniform/per-target member
theorem, reduce the squarefree degree in the rational argument \(b\),
remove Schinzel H, or change the six-count.  The remaining obstruction
is the uniform member/fixed-family detector problem.

## L27 [ONE-PIECE TRIANGULAR SHEAR, 2026-08-23]. Dyadic parity and every aligned odd target are locally solved

Return to the standard pullback $Z=z^3$ and tie
$s=(a-1)/2$.  On the square branch write
$$
X^2-\rho^2=A,\qquad \delta=-\frac{\rho^2}{A}.
$$
Use the fixed base-dependent shear
$$
\boxed{
\binom yr=
\begin{pmatrix}2&28\\s&1\end{pmatrix}
\binom X\rho,
\qquad
y=2X+28\rho,\quad r=sX+\rho.}
$$
Its determinant is
$$2-28s=2(1-14s),$$
which is nonzero on $\Phi$ because $s\in\mathbb Z_2$.  It uses no new
witness or inverse.

### L27a. Exact cover and finite-flat rank

Put $u=\rho^2$ and $V=X\rho$.  After multiplying the tied equation by
$A$, exact elimination gives
$$
\boxed{\begin{aligned}
C_2={}&788Au^2+
\bigl(4A^2-c^2-16AB(s^2+1)\bigr)u\\
&+\bigl(112Au-32ABs\bigr)V-16A=0.
\end{aligned}}
$$
With $\lambda=X+\rho\ne0$,
$$
X=\frac{\lambda+A/\lambda}{2},\qquad
\rho=\frac{\lambda-A/\lambda}{2}.
$$
Writing $L_0=4A^2-c^2-16AB(s^2+1)$, the exact degree-$8$ eliminant is
$$
\boxed{\begin{aligned}
\mathcal G_{28}(\lambda)={}&
197A(\lambda^2-A)^4
+L_0\lambda^2(\lambda^2-A)^2\\
&+28A(\lambda^2-A)^2(\lambda^4-A^2)
-32ABs\lambda^2(\lambda^4-A^2)
-64A\lambda^4,
\end{aligned}}
$$
and $\mathcal G_{28}=4\lambda^4C_2$.  Its leading and constant
coefficients are $225A$ and $169A^5$, respectively.  Hence, after the
already-guarded bridge denominators are cleared, this is an integral
finite-flat degree-$8$ cover.  The associated rank-$1$ contribution
would give the same nominal five-count as L24--L25 **only if** a
completeness theorem supplied a rational root.

### L27b. A single formula covers all dyadic parities

On $\Phi$, put $\rho=2t$ and choose the analytic odd square-root branch
$$X(t)^2=A+4t^2.$$
For $t\in1+4\mathbb Z_2$, the right side is $9\bmod32$.  Dividing $C_2$
by $16A$ gives
$$
G(t)=788t^4+56t^3X(t)+Lt^2-4BsX(t)t-1,
$$
where
$$L=A-4B(s^2+1)-\frac{c^2}{4A}$$
is odd.  Since $B=2b$ with $b$ odd and $v_2(c)\ge4$,
$$
G(1)\equiv4+5-1\equiv0\pmod8.
$$
Moreover $X'(t)=4t/X(t)$ and every term of $G'(t)$ except $2Lt$ has
valuation at least $2$, while
$$v_2(2Lt)=1.$$
Thus
$$v_2(G(1))\ge3>2v_2(G'(1))=2,$$
and strong Hensel supplies a root.  No parity split in $s$ remains.
The $512$ residue rows in the artifact exhaust the finite congruence
lemma; the valuation argument is the uniform proof.

### L27c. Every standard aligned odd target has a smooth residue

Let $w$ be odd and impose the L20 target stratum
$$
v_w(z)\ge1,\qquad v_w(B)=1,\qquad
A\in\mathbb Z_w^\times,\qquad
\left(\frac Aw\right)=-1.
$$
Then $D,N_g$ are units and
$$v_w(c)=6v_w(z)-2\ge4.$$
Modulo $w$, the two cover equations become
$$
X^2-\rho^2=A,\qquad
\rho^2(2X+28\rho)^2=16.
$$
Choose a sign $\epsilon\in\{\pm1\}$ and put
$$
\rho(2X+28\rho)=4\epsilon,\qquad x=\rho^2.
$$
Then
$$
X=\frac{2\epsilon}{\rho}-14\rho,
$$
and
$$
\begin{aligned}
N_{0,\epsilon}(x)&=195x^2-56\epsilon x+4=Ax,\\
N_{1,\epsilon}(x)&=N_{0,\epsilon}(x)-x=(A-1)x.
\end{aligned}
$$
It therefore suffices to find
$$
\chi(x)=\chi(N_{1,\epsilon}(x))=+1,\qquad
\chi(N_{0,\epsilon}(x))=-1.
$$
These conditions reconstruct $\rho$, $a$, $A=1+4a^2$ and $X$.

For $\epsilon=+1$ and $w\ge197$, consider
$$
I(x)=\frac{(1+\chi(x))(1+\chi(N_1(x)))(1-\chi(N_0(x)))}8.
$$
Expanding $I$ leaves two quadratic, two cubic, one quartic and one
quintic squarefree character sums.  The quadratic sums have absolute
value $1$; the Weil bounds for the remaining sums are respectively
$2\sqrt w,2\sqrt w,3\sqrt w,4\sqrt w$.  Hence
$$
\sum_x I(x)\ge\frac{w-11\sqrt w-2}{8}.
$$
The roots of $N_0$ contribute spurious total weight at most $1$.
The $(X,\rho)$-Jacobian is, up to the unit $A$,
$$
J=8\rho(2X+28\rho)
\bigl(X^2+28X\rho+\rho^2\bigr).
$$
It can vanish at at most the two values satisfying $195x^2=4$.
Consequently the number of smooth desired residues is at least
$$
\boxed{\frac{w-11\sqrt w-26}{8}>0\qquad(w\ge197).}
$$
For the $43$ odd primes $w<197$, exact exhaustion over both signs
finds a nonzero-Jacobian row at every prime
(`l27_triangular_shear.py::target_character_theorem`;
`data/l27_triangular_shear.jsonl`).  These are points on the **standard
ramified L20 stratum**, not the free unit-$b$ stratum of L25.

### L27d. Real place and strict global scope

For the guarded bridge sample
$$
(a,b,Z)=(1,-1,-10)
$$
the exact residual has opposite signs at
$$
\lambda=-\frac83,\qquad\lambda=-\frac52:
$$
$$
\frac{121384526885}{1167998976}>0,\qquad
-\frac{21088635}{506944}<0.
$$
The intermediate value theorem gives one real coupled point.

L27 therefore removes, in one formula, the L25 parity split, the
$w=3,5$ residue gaps, and the aligned-target mismatch.  It does **not**
give a rational root of $\mathcal G_{28}$, control every remaining
finite place simultaneously, prove a member theorem, produce a
five-count, or remove Schinzel H.  The global bridge-specialized
degree-$8$ cover is now the sole obstruction in this sheared route.

### L27e. The natural linear-\(B\) section does not exist

The eliminant is even in \(\lambda\).  Put \(m=\lambda^2\) and separate
the terms involving \(B=2b\):
$$
\mathcal G_{28}(m)
=G_0(A,m)-c^2m(m-A)^2
-16ABm(m-A)
\bigl((s+1)^2m-(s-1)^2A\bigr),
$$
where \(G_0\) is independent of \(B\) and \(c\).  Thus the tempting
strategy of first cancelling all linear \(B\)-dependence has only three
possibilities.  The factor \(m=0\) is \(\lambda=0\), outside the conic
parameterization.  At \(m=A\), the full eliminant is
\(-64A^3\ne0\).  For \(s\ne\pm1\), the last factor would give
$$
\lambda^2=A\left(\frac{s-1}{s+1}\right)^2,
$$
forcing \(A\) to be a rational, hence \(2\)-adic, square; this
contradicts \(A\equiv5\bmod8\).  At \(s=1\) it reduces to \(m=0\), and
at \(s=-1\) it equals \(-4A\).

Therefore no rational point on the canonical cover can arise by
annihilating its linear \(B\)-coefficient
(`l27_triangular_shear.py::b_coefficient_no_go`;
`data/l27_triangular_shear.jsonl`).  This excludes exactly that section
ansatz; it is **not** a no-point theorem for the full degree-\(8\)
cover.


### L27f. The four canonical hyperbola sections are dyadically impossible

The rational hyperbola \(X^2-\rho^2=A\) has the four canonical
parameter values
\[
\lambda\in\{1,-1,A,-A\}.
\]
For \(\lambda=1,A\), respectively,
\[
X=1+2a^2,\qquad \rho=\mp2a^2;
\]
changing the sign of \(\lambda\) negates both \(X\) and \(\rho\), so
\(u=\rho^2\) and \(V=X\rho\) are unchanged.  On \(\Phi\), \(a,b\) are
\(2\)-adic units and \(c\in16\mathbb Z_2\).  Substitution in the exact
formula of L27a gives, for both signs of \(\rho\),
\[
\boxed{C_2\equiv128\pmod{256}.}
\]
In particular none of \(\lambda=\pm1,\pm A\) is a root of
\(\mathcal G_{28}\).  The congruence is an identity on odd residue
classes: the producer exhausts \(2\cdot32^2=2048\) pairs
\((\operatorname{sgn}\rho,a,b)\bmod64\), while the omitted
\(-c^2\rho^2\) term has valuation at least \(10\)
(`l27_triangular_shear.py::canonical_hyperbola_sections_no_go`;
`data/l27_triangular_shear.jsonl`).

### L27g. The coordinate section \(y=0\) is impossible over
\(\mathbb Q_3\)

If \(y=2X+28\rho=0\), then \(X=-14\rho\), and the hyperbola equation
forces
\[
195\rho^2=A=1+4a^2,\qquad
(2a)^2-195\rho^2=-1.
\]
This last conic has no \(\mathbb Q_3\)-point.  Indeed the valuations of
its two left terms are \(2v_3(2a)\) and \(1+2v_3(\rho)\), of opposite
parity, so they cannot cancel.  To equal a unit, the first must have
valuation \(0\); reduction modulo \(3\) would then make a square equal
to \(-1\), impossible.  This excludes only the coordinate-hyperplane
section \(y=0\), not the full degree-\(8\) cover.


## L28 [UNIFORM TRACE-FIELD CLOSURE, 2026-08-24]. The L26 quartic and its bad-sign cover are irreducible at \(2\)

Retain the canonical L26 hypotheses
$$
v_2(a)=0,\qquad Z=z^2\ne0,\qquad
A=1+4a^2,\qquad D=1-Z-a^2Z^2,
$$
and put \(t=v_2(Z)\in2\mathbb Z\).  Thus
\(A\equiv5\bmod8\), \(D\ne0\), and
$$
v_2(D)=
\begin{cases}
0,&t\ge0,\\
2t,&t<0.
\end{cases}
$$
Let \(T\) be the canonical trace quartic from L26b.

### L28a. The trace quartic is uniformly irreducible

Set \(w=u-2\).  Exact expansion gives
$$
T(w+2)=d_4w^4+d_2w^2+d_1w+d_0,
$$
where
$$
\begin{aligned}
d_4&=4a^8AZ^4,\\
d_2&=-128a^{12}Z^4,\\
d_1&=-32A^3Z^2,\\
d_0&=16A^2D^2+\frac{1024a^{16}Z^4}{A}-128A^3Z^2.
\end{aligned}
$$
The nonzero coefficient valuations are
$$
\begin{array}{c|rrrr}
 &d_0&d_1&d_2&d_4\\ \hline
t\ge0&4&5+2t&7+4t&2+4t\\
t<0&4+4t&5+2t&7+4t&2+4t .
\end{array}
$$
Consequently the Newton polygon has one segment, of slope
\(t-\tfrac12\) for \(t\ge0\) and \(-\tfrac12\) for \(t<0\).
Every root therefore has half-integral valuation, so \(T\) has no
linear factor over \(\mathbb Q_2\).

It remains to exclude a product of two quadratics.  Divide by \(d_4\)
and write
$$
f(w)=w^4+pw^2+qw+r,
$$
where
$$
\begin{aligned}
p&=-\frac{32a^4}{A},\\
q&=-\frac{8A^2}{a^8Z^2},\\
r&=\frac{4AD^2}{a^8Z^4}
   +\frac{256a^8}{A^2}
   -\frac{32A^2}{a^8Z^2}.
\end{aligned}
$$
If
$$
f=(w^2+\alpha w+\beta)(w^2-\alpha w+\gamma),
$$
then \(q\ne0\) forces \(\alpha\ne0\), and
\(Y=\alpha^2\) is a nonzero square root of the Ferrari resolvent
$$
\boxed{
H(Y)=Y^3+2pY^2+(p^2-4r)Y-q^2.}
$$
Put \(h=p^2-4r\).  Direct valuation gives
$$
\begin{array}{c|rrrr}
 &v_2(p)&v_2(q)&v_2(r)&v_2(h)\\ \hline
t\ge0&5&3-2t&2-4t&4-4t\\
t<0&5&3-2t&2&4 .
\end{array}
$$
After the displayed power of \(2\) is removed,
$$
\frac{h}{2^{\,4-4t}}\equiv3\pmod8\quad(t\ge0),
\qquad
\frac h{16}\equiv3\pmod8\quad(t<0),
$$
while
$$
\frac{q^2}{2^{\,6-4t}}\equiv1\pmod8.
$$

The Newton polygon of \(H\) leaves only the valuations in the second
column below.  If \(Y\) were a square, its normalized unit would be
\(1\bmod8\).  Substitution gives the listed nonzero normalized residue:
$$
\begin{array}{c|c|c|c}
\text{case}&v_2(Y)&\text{power dividing }H(Y)&
H(Y)/2^{\text{power}}\pmod8\\ \hline
t=0&2&6&3\\
t>0&2&6-4t&2\\
t>0&2-2t&6-6t&4\\
t<0&2&6&4\\
t<0&2-4t&6-4t&2 .
\end{array}
$$
Here \(t>0\) means \(t\ge2\), and \(t<0\) means \(t\le-2\), because
\(Z=z^2\); this is exactly what makes the omitted terms vanish modulo
\(8\) in the two residue-\(4\) rows.  Thus \(H\) has no square root in
\(\mathbb Q_2\), and \(f\) has no quadratic factor.  Together with the
linear-factor exclusion,
$$
\boxed{T\text{ is irreducible over }\mathbb Q_2,\text{ hence over }\mathbb Q}
$$
for every canonical L26 parameter
(`l28_trace_field.py::local_irreducibility_replay`;
`data/l28_trace_field.jsonl`).

### L28b. The bad-sign squareclass is uniformly nontrivial

Let \(\theta\) be a root of \(T\) and
\(K=\mathbb Q(\theta)\).  Since the leading coefficient is
\(4a^8AZ^4\),
$$
\begin{aligned}
N_{K/\mathbb Q}\bigl(2(\theta+2)\bigr)
&=\frac{16T(-2)}{4a^8AZ^4}\\
&=\frac{64E}{a^8Z^4A^2},
\end{aligned}
$$
where
$$
\boxed{
E=A^3D^2+64a^8Z^4(a^4-A)^2.}
$$
For a \(2\)-adic unit \(a\),
$$v_2(a^4-A)=2,$$
and therefore
$$
\frac{E}{A^3D^2}
=1+\frac{64a^8Z^4(a^4-A)^2}{A^3D^2}
\in1+2^{10}\mathbb Z_2.
$$
The correction is a \(2\)-adic square.  Every factor outside \(E\) in
the norm formula is also a square, so
$$
\boxed{
N_{K/\mathbb Q}\bigl(2(\theta+2)\bigr)
\in A\,\mathbb Q_2^{\times2}.}
$$
Because \(A\equiv5\bmod8\), this norm is nonsquare.  Hence
\(2(\theta+2)\) is nonsquare in \(K\), and
$$
[K(\sqrt{2(\theta+2)}):\mathbb Q]=8.
$$
Equivalently, Capell's criterion now gives the unconditional-on-the-
parameters conclusion
$$
\boxed{T(v^2/2-2)\text{ is irreducible over }\mathbb Q_2
\text{ and over }\mathbb Q.}
$$
The producer replays \(2{,}600\) exact unit-congruence rows, \(175\)
exact norm rows, and an independent mod-\(41\) degree-\(8\) sample
certificate (`l28_trace_field.py`;
`data/l28_trace_field.jsonl`).

TraceFieldAudit independently rederived every coefficient, Newton
polygon, Ferrari-resolvent case, norm squareclass and Capell implication;
verdict SOUND at confidence \(0.99\), with no finding
(`agent://TraceFieldAudit`).

### L28c. Consequence and strict frontier

The conditional qualifications in L26d are now removed: the trace
algebra is always a quartic field and its bad-sign quadratic extension
is always nontrivial on the canonical even-pullback domain.  This
closes the first algebraic target named in L26e.

At the L28 stage this did **not** decide the distinct reciprocal-lift
squareclass \(\theta^2-4\).  L29 now classifies that squareclass
completely and supplies a dyadically compressed pullback.  Neither
result produces a globally good \(b\), solves the two-large-bad-divisor
parity problem, removes Schinzel H, or improves the formal six-count.
The L26 route is now blocked by the global member/parity step rather than by
trace-field or bad-character nontriviality.

## L29 [RECIPROCAL LIFT CLASSIFIED; DYADIC COMPRESSION, 2026-08-24].
The missing local squareclass has four exact strata

Retain the square-\(Z\) L26 hypotheses and notation of L28:
\[
v_2(a)=0,\qquad Z=\zeta^2\ne0,\qquad t=v_2(Z)\in2\mathbb Z,
\]
and let \(\theta\) be a root of the irreducible trace quartic \(T\).
Put \(K=\mathbb Q_2(\theta)\).

### L29a. Complete classification of \(\theta^2-4\)

Write \(w=\theta-2\).  L28a gives
\[
v_2(w)=
\begin{cases}
\frac12-t,&t\ge0,\\[2pt]
\frac12,&t<0.
\end{cases}
\]
For \(t\ge2\),
\[
\theta^2-4=w(w+4)=w^2\left(1+\frac4w\right),
\qquad
v_2(4/w)=t+\frac32>2.
\]
Strong Hensel at \(1\) makes \(1+4/w\) a square in \(K\).  Hence
\(\theta^2-4\) is a square for every \(t\ge2\).

The remaining strata require the reciprocal lift itself.  Set
\(b=1+x\).  Since
\[
b+b^{-1}=2+\frac{x^2}{1+x},
\]
the monic translate of \(P_{\rm rec}(b)=b^4T(b+b^{-1})\) is
\[
\boxed{
S(x)=x^8+p\,x^4(1+x)^2+q\,x^2(1+x)^3+r(1+x)^4,}
\]
with \(p,q,r\) from L28a.  For \(t\le0\), its Newton polygon is the
single segment of slope \(-1/4\).  Thus every proper factor has degree
divisible by \(4\); a factorization must be \(4+4\), with both monic
factors \(2\)-Eisenstein.

Direct reduction of the displayed formula gives the following
ascending coefficient vectors modulo \(64\):
\[
\begin{array}{c|rrrrrrrrr}
&x^0&x^1&x^2&x^3&x^4&x^5&x^6&x^7&x^8\\ \hline
t\le-4&20&16&56&16&52&0&32&0&1\\
t=-2&52&16&56&16&20&0&32&0&1\\
t=0&52&16&48&56&60&56&32&0&1 .
\end{array}
\]
For completeness, the only input is
\[
v_2(p)=5,\qquad
q\equiv\begin{cases}56\pmod{64},&t=0,\\0\pmod{64},&t\le-2,\end{cases}
\]
and
\[
r\equiv\begin{cases}
20\pmod{64},&t\le-4,\\
52\pmod{64},&t=-2,0.
\end{cases}
\]
The last congruences follow immediately from the formula for \(r\):
for \(t<0\), divide \(D\) by \(-a^2Z^2\); the resulting unit is
\(1\bmod16\) for \(t\le-4\) and \(5\bmod16\) for \(t=-2\).  For \(t=0\),
odd squares modulo \(16\) give the same residue \(52\).

At \(t=0\), no product of two monic Eisenstein quartics has the third
vector modulo \(16\).  At \(t\le-4\), no such product has the first
vector modulo \(64\).  These are finite congruence lemmas, not a
computer-factorization oracle: choose the four lower coefficients of
the first monic quartic; the four highest product coefficients uniquely
determine the second.  Exhaustion checks respectively \(2048\) and
\(524{,}288\) possible first factors, with zero survivors.
Every exact \(\mathbb Q_2\)-factor is integral and Eisenstein by the
first Newton polygon, so its coefficients would reduce to one of these
pairs; no approximate-factor or precision-bound inference is used.
Therefore \(\theta^2-4\) is nonsquare for \(t=0\) and for \(t\le-4\).

The exceptional stratum \(t=-2\) goes the other way, by an ordinary
resultant Newton polygon.  Put
\[
\phi=x^4+6x^3+6,\qquad
\mathcal R(Y)=\operatorname{Res}_x(S(x),Y-\phi(x)).
\]
The resultant is monic of degree \(8\), and its roots are exactly
\(\phi(\alpha)\) as \(\alpha\) ranges over the roots of \(S\).  At
\(t=-2\), write
\[
p=32P,\qquad q=128Q,\qquad r=52+64R,
\qquad P,Q,R\in\mathbb Z_2.
\]
Computing \(\mathcal R\) as the characteristic polynomial of
multiplication by \(\phi(x)\) in
\(\mathbb Q[p,q,r][x]/(S)\), then making this substitution, gives
\[
\begin{array}{c|rrrrrrrrr}
j&0&1&2&3&4&5&6&7&8\\ \hline
v_2([Y^j]\mathcal R)
&21&19&\ge17&15&10&\ge11&6&5&0 .
\end{array}
\]
The three vertices are exact and independent of \(P,Q,R\):
\[
(0,21),\qquad(4,10),\qquad(8,0).
\]
All other points lie strictly above the two joining segments.  Hence
the Newton polygon of \(\mathcal R\) has slopes
\[
-\frac{11}{4}\quad\text{and}\quad-\frac52,
\]
each of horizontal length \(4\).

If \(S\) were irreducible over \(\mathbb Q_2\) and \(\alpha\) were one
root, every \(\mathbb Q_2\)-embedding of
\(\mathbb Q_2(\alpha)\) into \(\overline{\mathbb Q}_2\) would preserve
the unique extended valuation.  Since \(\phi\) has rational
coefficients, all eight conjugates \(\phi(\sigma\alpha)\) would have
the same valuation, forcing \(\mathcal R\) to have a one-segment Newton
polygon.  The two displayed slopes contradict this.  Thus \(S\) is
reducible; its first Newton polygon permits only \(4+4\), so it is a
product of two quartics.

Since \(T\) is irreducible of degree \(4\), this factorization is
equivalent to the quadratic
\[
b^2-\theta b+1
\]
splitting over \(K\), hence to \(\theta^2-4\in K^{\times2}\).

Combining the four cases,
\[
\boxed{
\theta^2-4\in K^{\times2}
\iff t=-2\ \text{or}\ t\ge2;}
\qquad
\boxed{
\theta^2-4\notin K^{\times2}
\iff t=0\ \text{or}\ t\le-4.}
\]
The producer replays \(5120\) unit rows, a division-free symbolic
resultant characteristic polynomial plus \(1024\) unit-congruence rows
at \(t=-2\), and the exact finite factor searches
(`l29_reciprocal_frontier.py`;
`data/l29_reciprocal_frontier.jsonl`).

### L29b. A target-preserving pullback forces the split stratum

For \(z\ne0\), replace the square pullback by
\[
\boxed{Z=\frac{8z^2}{1+z^2}.}
\]
The denominator is nonzero over \(\mathbb Q\).  For every odd prime
\(p\), a three-case valuation check gives
\[
\boxed{v_p(Z)>0\iff v_p(z)>0.}
\]
Indeed the denominator is a unit when \(v_p(z)>0\), has the same
valuation as \(z^2\) when \(v_p(z)<0\), and can only turn a unit input
into a pole when \(v_p(z)=0\).  At \(2\),
\[
v_2(Z)=
\begin{cases}
3+2v_2(z),&v_2(z)>0,\\
2,&v_2(z)=0,\\
3,&v_2(z)<0,
\end{cases}
\]
because an odd \(2\)-adic unit has square \(1\bmod8\).  Hence
\(v_2(Z)\ge2\) on every nonzero rational fibre.

Although this \(Z\) need not be a rational square, the \(t>0\)
Ferrari-resolvent calculation of L28 uses only \(t\ge2\), not its
parity.  Its two possible square-root cases still have normalized
residues \(2\) and \(4\bmod8\).  The norm calculation of L28b is
unchanged.  Consequently this compressed pullback has, uniformly,
\[
\begin{aligned}
&T\text{ irreducible over }\mathbb Q_2,\\
&2(\theta+2)\notin\mathbb Q_2(\theta)^{\times2},\\
&\theta^2-4\in\mathbb Q_2(\theta)^{\times2}.
\end{aligned}
\]
Thus the bad-sign cover remains degree \(8\), while the distinct
reciprocal octic is always \(4+4\) over \(\mathbb Q_2\).  The odd target
set, denominator guards, soundness, and formal witness count are
unchanged.  At \(z=0\), one has \(c=\eta=0\), and the canonical
binary equation is solved by \(y=2/a^2,\ r=0\); only the quartic
classification assumes \(Z\ne0\).

### L29c. Freezing the moving quadratic field

For a target prime \(w\), the slice
\[
\boxed{b=w\rho^2}
\]
is compatible with \(v_w(b)=1\) and with \(b\) a \(2\)-adic unit when
\(\rho\) is a unit at both places.  More importantly,
\[
(M,2b)_v=(M,2w)_v
\]
at every place.  The global member problem becomes the fixed-field norm
family
\[
\boxed{U^2-2wV^2=P_{\rm rec}(w\rho^2)W^2.}
\]
This removes the moving quadratic extension but not parity.  It has no
generic rational section: at \(\rho=0\),
\[
P_{\rm rec}(0)=A(2a^4Z^2)^2,\qquad
(A,2w)_w=\left(\frac Aw\right)=-1
\]
on the guarded target stratum, so the proper special conic has no
\(\mathbb Q_w\)-point.

The slice is nevertheless arithmetically live.  Eight named target
fibres \(w=3,5,7,11,13,17,19,23\) have exact globally soluble members.
On the fixed row \((w,a,Z)=(3,5,9)\), the \(33\) positive odd
\(w\)-adic units \(\rho<100\) give \(24\) fully factored decisions:
\(11\) with no bad place, \(12\) with two, and \(1\) with four; \(9\)
factor/primality budget refusals are excluded.  This is EVIDENCE, not a
density or member theorem.
It shows both the promise and the surviving parity wall.

### L29d. Strict frontier

L29 closes the local reciprocal-lift question and supplies two new
exact architectures.  It does **not** produce a globally good member
uniformly, prove the inert-prime or two-large-divisor estimate, find a
rational point on the L27 cover, remove Schinzel H, improve the
six-count, or solve H10/\(\mathbb Q\).

Recent source checks do not silently fill that gap.  Loughran--Matthiesen
(`arXiv:1904.12845`) first misses the non-split fibre over a non-rational
closed point; Shute (`arXiv:2209.08949`) first misses the generic
high-degree fixed-field slice; Diao's random-binary-form theorems
(`arXiv:2506.18065`) are almost-all coefficient results, not statements
for this fixed low-dimensional family; and Wang's newest averaged
dynamical Chowla theorem (`arXiv:2608.16108`, 2026-08-17) does not
control the fixed sign-decorated octic sequence.


## L30 [A LOWER-DEGREE QUARTIC SECTION; BRANCH AND TRACE RIGIDITY, 2026-08-24].
Fixing one norm witness halves the L27 cover degree

Return to the standard L24--L27 pullback and square branch
\[
a=1+2s,\qquad A=1+4a^2,\qquad B=2b,\qquad Z=z^3,
\]
\[
X^2-\rho^2=A,\qquad \delta=-\frac{\rho^2}{A}.
\]
The bridge value is still
\[
c=\frac{a^2Z^2N_g}{Ab^2D},\qquad
D=1-Z-a^2Z^2,\qquad
N_g=16a^4b^2-A(b-1)^4.
\]

### L30a. The constant-two section is quartic, not octic

Impose the nonlinear specialization
\[
\boxed{y=2,\qquad r=sX+\rho.}
\]
It uses no new witness or inverse.  The two exact equations are
\[
X^2-\rho^2=A
\]
and
\[
\boxed{
C_2=-\rho^2(c^2-4A)
-16AB\bigl((sX+\rho)^2-As^2\bigr)-16A=0.}
\]
Put
\[
\lambda=X+\rho,\qquad m=\lambda^2,
\]
so
\[
X=\frac{\lambda+A/\lambda}{2},\qquad
\rho=\frac{\lambda-A/\lambda}{2}.
\]
Define
\[
Q_s(m)=(s+1)^2m^2-2A(s^2+1)m+(s-1)^2A^2.
\]
Direct elimination gives
\[
\boxed{
H_2(m)=-(c^2-4A)(m-A)^2-16ABQ_s(m)-64Am
=4mC_2.}
\]
Thus the cover is quadratic in \(m\) and even quartic in \(\lambda\),
instead of the L27 degree-\(8\) cover.  Its leading coefficient is
\[
[m^2]H_2=4A-c^2-16AB(s+1)^2.
\]
On \(\Phi\), \(v_2(c)\ge4\), \(v_2(B)=1\), and \(A\) is a unit, so this
coefficient has valuation exactly \(2\).  The cover is genuinely
quartic at every \(\Phi\)-base.

The \(m\)-discriminant has the exact lower-degree form
\[
\boxed{\operatorname{disc}_m(H_2)=256A^2J,}
\]
\[
\boxed{
J=4A+16-c^2-16AB(s^2+1)+16A^2B^2s^2.}
\]
These identities are replayed on \(400\) exact rational instances
(`l30_quartic_frontier.py`;
`data/l30_quartic_frontier.jsonl`).

### L30b. Uniform dyadic point

On \(\Phi\), write
\[
h(m)=\frac{H_2(m)}A
=\left(4-\frac{c^2}{A}\right)(m-A)^2
-32bQ_s(m)-64m.
\]
At \(m=1\),
\[
Q_s(1)
=-16a^2\left(s(1+2a^2)-a^2(s^2+1)\right).
\]
Since \(a\) is odd, \(a^4\equiv1\pmod{16}\), while
\(c^2/A\in2^8\mathbb Z_2\).  Hence
\[
v_2(h(1))\ge9.
\]
Put
\[
g(t)=\frac{h(1+8t)}{2^8}.
\]
Then \(g(0)\in2\mathbb Z_2\), whereas
\[
g'(0)=
\left(4-\frac{c^2}{A}\right)\frac{1-A}{16}
-bQ_s'(1)-2
\in\mathbb Z_2^\times.
\]
Indeed \(Q_s'(1)\) is even and the first term is
\(-a^2\) modulo \(2\).  Ordinary Hensel therefore gives
\(t\in\mathbb Z_2\) with \(H_2(1+8t)=0\).
Every unit \(1+8t\) is a square in \(\mathbb Q_2\), so
\(m=\lambda^2\) and the constant-two cover has a
\(\mathbb Q_2\)-point on every \(\Phi\)-stratum.  The producer exhausts
\(2048\) residue rows; the displayed valuation argument is the uniform
proof.

### L30c. Every aligned odd target has a selected fibre point

At a standard target \(w\),
\[
v_w(c)\ge1,\qquad v_w(B)=1,\qquad
\left(\frac Aw\right)=-1.
\]
Modulo \(w\), the two equations reduce to
\[
X^2-\rho^2=A,\qquad 4A(\rho^2-4)=0.
\]
For \(w\ge5\), choose \(\rho=2\).  It is enough to find \(a\) with
\[
\chi(1+4a^2)=-1,\qquad
\chi(5+4a^2)=+1.
\]
For \(w\ne5\), expand
\[
\frac{(1-\chi(1+4a^2))(1+\chi(5+4a^2))}{4}.
\]
The two quadratic character sums are both \(-1\).  The remaining
squarefree quartic sum has absolute value at most \(3\sqrt w\), and the
two zero loci contribute at most \(2\) spurious units.  Therefore
\[
\boxed{N_w\ge\frac{w-3\sqrt w}{4}-2>0\qquad(w\ge23).}
\]
Exact exhaustion supplies rows for
\(w=5,7,11,13,17,19\).  At each such row the fibre Jacobian in
\((X,\rho)\) has determinant
\[
16AX\rho\ne0.
\]

For \(w=3\), the first residue point is singular in the fibre
coordinates, but the exact guarded aligned fibre
\[
(a,b,z,Z)=(5,3,3,27)
\]
is soluble.  At \(m=1\), its genuine quadratic satisfies
\[
v_3(H_2(1))=3,\qquad v_3(H_2'(1))=1.
\]
Strong Hensel gives a root \(m\in1+9\mathbb Z_3\), hence a
\(\mathbb Q_3\)-square, and therefore a rational \(\lambda\) over
\(\mathbb Q_3\).  The underlying residue geometry is also smooth on the
total space: at
\[
(a,A,X,\rho)=(2,2,0,1)\pmod3
\]
the Jacobian in the columns \((a,\rho)\) has determinant nonzero.

The selected \(a\) matters.  The nearby guarded aligned control
\[
(a,b,z,Z)=(1,3,3,27)
\]
is \(\mathbb Q_3\)-empty: here
\[
v_3(c)=4,\qquad v_3([m^2]H_2)=0,\qquad v_3(J)=1,
\]
so \(\operatorname{disc}_m(H_2)=256A^2J\) is nonsquare.  Thus L30
supplies a fibre point for every odd target prime, including \(3\), but
does not claim uniform solubility for every preselected \(a\).

### L30d. Real place and exact remaining global problem

The guarded sample
\[
(a,b,Z)=\left(1,-1,\frac18\right)
\]
has
\[
c=-\frac{64}{275},\qquad
4A-c^2-16AB=\frac{13608404}{75625}>0.
\]
Since \(s=0\), the equation is
\[
\rho^2(4A-c^2-16AB)=16A,
\]
and hence
\[
\rho^2=\frac{1512500}{3402101}>0,\qquad
X^2=\frac{18523005}{3402101}>0.
\]
This proves a nonempty real open stratum.

The remaining global condition is exact: find a rational
\(\lambda\ne0\) with
\[
\boxed{H_2(\lambda^2)=0}
\]
while meeting the other controlled places.  A bounded exact probe over
\(a\in\{1,3,5,7\}\) and rational numerator/denominator grid
\((8,5)\) checked \(13{,}224\) guarded rows and found no rational
\(\lambda\).  This is EVIDENCE only.  The degree drop \(8\to4\) is
proved; a rational-point theorem, a complete five-count, and any
unconditional member conclusion remain open.

### L30e. The full square-branch pencil cannot hide a generic section

Let \(d\in\mathbb Q(a,b,Z)^\times\) be any rational tie for the L11c
square-branch parameter and put
\[
\ell_d=\frac{A-d^2}{2d},\qquad
\mu_d=\frac{A+d^2}{2d}.
\]
Then
\[
\mu_d^2-\ell_d^2=A,\qquad
\delta_d=-\frac{\ell_d^2}{A},\qquad
\alpha_d=-A\delta_d=\ell_d^2.
\]
The tied equation is
\[
\ell_d^2y^2-16Br^2=M_d,\qquad
M_d=16+\frac{\ell_d^2}{A}c^2-16ABs^2.
\]
Consequently every square branch has the **same** target norm field
\[
\mathbb Q(\sqrt B)=\mathbb Q(\sqrt{2b});
\]
varying \(d\) does not create a quadratic twist.

Clear the bridge square by
\[
P_d=A^2b^4D^2M_d.
\]
One obtains
\[
P_d=16A^2D^2b^4
+\frac{\ell_d^2a^4}{A}Z^4N_g^2
-32A^3s^2D^2b^5
\]
and the source-pencil identity
\[
\boxed{
AP_d=(\ell_da^2Z^2N_g)^2
+A(1-2Abs^2)(4ADb^2)^2.}
\]

There is no generic rational norm section for **any**
\(d\in\mathbb Q(a,Z)(b)^\times\).  Indeed let
\(k=\operatorname{ord}_b(\ell_d)\).  If \(d\) has a zero or pole at
\(b=0\), direct dominance gives \(k<0\); if \(d\) is a unit, then
\(\ell_d(0)=0\) would force \(d(0)^2=A\), impossible over
\(\mathbb Q(a,Z)\).  Thus \(k\le0\).  The middle term of \(P_d\) has
the unique least \(b\)-valuation \(2k\), and its leading unit has
squareclass \(A\).  By contrast, in
\(\mathbb Q(a,Z)((b))(\sqrt{2b})\), the two terms of a norm
\[
Y^2-2bR^2
\]
have valuations of opposite parity.  A norm of even valuation is led
by \(Y^2\) and therefore has square leading unit.  Contradiction.
This rules out rational, Laurent, pole, and nonlinear generic
\(d\)-sections; it does not rule out specialized rational members.
The symmetries
\[
d\sim-d\sim A/d\sim-A/d
\]
also yield the same \(\delta_d,M_d,P_d\), so these are not independent
branch ansatzes.

The sole elementary nonunit dyadic class not eliminated by the coarse
dominance screen is also not a uniform escape.  Taking
\[
d=2a,\qquad \ell_d=\frac1{4a},
\]
the guarded \(\Phi\)-base
\[
(a,b,z,Z)=(1,-15,-9,-729)
\]
has \(v_2(c)\ge4\) but
\[
\boxed{(M_d,2b)_2=-1.}
\]
Thus \(d=2a\) is not uniformly dyadically admissible.  This is one
exact counterexample to a uniform tie, not a no-member theorem for
specialized \(d=2a\) fibres.

### L30f. All \(B\)-independent linear-shear cancellations are rigid

For a general linear shear
\[
\binom yr=
\begin{pmatrix}p&q\\r&t\end{pmatrix}\binom X\rho
\]
whose entries are independent of \(B=2b\), the coefficient of \(B\)
in the even degree-\(8\) eliminant factors as
\[
\boxed{-16ABmH_L(m),}
\]
\[
H_L(m)=(r+t)^2m^2
-2A(t^2-r^2+2s^2)m
+A^2(r-t)^2.
\]
Its discriminant is
\[
\boxed{16A^2s^2(s^2+t^2-r^2).}
\]
For \(s\ne0\), a rational nonzero root would have
\[
\frac mA=
\left(\frac{\kappa\pm s}{r+t}\right)^2,\qquad
\kappa^2=s^2+t^2-r^2.
\]
If \(r+t=0\), the same conclusion is \(m=A(r/s)^2\).
But \(m=\lambda^2\), while \(A\equiv5\bmod8\) is not a
\(\mathbb Q_2\)-square.  Hence no such rational root exists.  This
strictly extends L27e from its one fixed shear to every
\(B\)-independent linear shear.  The edge \(s=0\) remains a residual
conic.  If shear entries depend on \(b\), they no longer have a literal
linear-\(B\) coefficient; the b-adic norm-section theorem of L30e,
not this factorization, is the applicable generic obstruction.

### L30g. Trace-base descent is exact, but the reciprocal lift is separate

Return here to the canonical L26 **square pullback** \(Z=z^2\)
(the L29 compressed pullback is also covered on its \(v_2(Z)\ge2\)
stratum).  On that canonical branch,
\[
P_{\rm rec}=D^2A^2b^4M=b^4T(u),\qquad u=b+b^{-1},
\]
so in fact
\[
\boxed{M=\frac{T(u)}{D^2A^2}.}
\]
Also
\[
\frac{u+2}{b}=\left(\frac{b+1}{b}\right)^2.
\]
Thus the conic descends on \(b\ne0,-1\) to
\[
U^2-2(u+2)V^2=T(u)W^2.
\]
After the base change
\[
u+2=w\rho^2
\]
this becomes the fixed-field degree-\(8\) bundle
\[
\boxed{
U^2-2wV^2=T(w\rho^2-2)W^2.}
\]
Its leading coefficient is
\[
A(2a^4Z^2w^2)^2.
\]
If \(\theta\) is a root of \(T\), a root of the composed polynomial
squares to \((\theta+2)/w\).  L28 gives
\[
N_{\mathbb Q_2(\theta)/\mathbb Q_2}((\theta+2)/w)
\in A\mathbb Q_2^{\times2},
\]
which is nonsquare.  Capell therefore proves
\[
T(w\rho^2-2)
\]
irreducible over \(\mathbb Q_2\), hence over \(\mathbb Q\), for every
nonzero rational \(w\).

This auxiliary bundle does **not** remove the rational reciprocal lift.
A rational \(b\) additionally requires
\[
u^2-4\in\mathbb Q^2,
\]
equivalently
\[
\boxed{\lambda^2=w(w\rho^2-4).}
\]
Its parametrization is
\[
\rho=q+\frac1{wq},\qquad
\lambda=wq-\frac1q,\qquad
b=wq^2.
\]
Thus imposing the lift recovers exactly the L29 fixed-squareclass slice
and its parity wall.  Treating a rational point on the trace-base bundle
as a reciprocal member without this second conic would be circular.

The auxiliary bundle has eight geometric simple poles with constant
residue, which makes the 2025 degree-\(8\) conic-bundle work of
Casarotti--Gammelgaard--Massarenti (`arXiv:2511.17213`) relevant.
Their direct Mestre criterion nevertheless fails on this exact diagonal
subfamily: its two natural orientations have leading squareclasses
\(-2wA\) and \(A\), respectively.  Frei--Sofos
(`arXiv:2604.07047`) supplies a new analytic Hilbert-symbol detector,
but its main estimate is an \(L^2\) average over the full coefficient
box, not a pointwise theorem for this fixed form on a prime progression.
The first analytic input remains L23b's explicit two-large-bad-divisor
estimate.  The 2026 finite-extension result
(`arXiv:2607.25287`) does not descend a point to \(\mathbb Q\).

### L30h. Strict frontier

L30 replaces the L27 degree-\(8\) equation by a degree-\(4\) cover that
is uniform at \(2\), fibre-regular at every aligned \(w\ge5\), and has
an exact strong-Hensel fibre at \(w=3\).  It proves two broad
generic-section rigidity theorems and identifies exactly why the
tempting trace-base/Mestre route is not yet a reciprocal member theorem.
It does **not** prove a rational root of \(H_2(\lambda^2)\), a globally
good member, the fixed-family dispersion estimate, removal of
classical Schinzel H, a five- or six-count improvement, or
H10/\(\mathbb Q\).

## L31 [ONE EXACT RECIPROCAL MEMBER; THE CUBE AND DISPERSION WALLS, 2026-08-24].
The first fixed-field member is promoted from evidence, while the other three requested closures remain open

L31 separates four logically different questions.  It proves one exact
globally good member on the L29 fixed-field slice.  It also finds an exact
rational point on the constant-two quartic immediately before the mandatory
cube pullback, sharpens the cube obstruction to two genus-\(2\) curves, and
identifies the exact analytic statement that would replace classical
Schinzel H.  The last three are reductions, not closures.

### L31a. Four canonical quartic sections are uniformly absent

Retain L30's notation.  For \(m=1\), put
\[
T=s(1+2a^2)-a^2(s^2+1).
\]
If \(a\) is odd, then \(T\) is odd.  Direct substitution gives
\[
\frac{H_2(1)}A
=64(a^4-1)+512ba^2T-\frac{16a^4c^2}{A}.
\]
On \(\Phi\), \(a^4-1\in16\mathbb Z_2\), \(b\) is a unit, and
\(c\in16\mathbb Z_2\).  Hence
\[
\frac{H_2(1)}4\equiv128\pmod{256}.
\]
At \(m=A^2\), the identity
\[
\frac{Q_s(A^2)}{A^2}\equiv16\pmod{32}
\]
gives the same congruence after division by \(4A^2\).  Therefore
\[
\boxed{\lambda\in\{\pm1,\pm A\}\Longrightarrow C_2\equiv128\pmod{256}}
\]
on every \(\Phi\)-base.  None is a root.  The producer exhausts
\(2\cdot32^2\cdot8=16{,}384\) residue rows; the displayed congruences
are the uniform proof.

This extends L27f to the different constant-two equation.  It excludes
four sections only, not the full quartic.

### L31b. A unique constant-\(c\) bridge section gives an exact off-cube point

Set \(a=1\), \(A=5\), and write
\[
\upsilon=\frac{(b-1)^2}{b}.
\]
For a finite nonzero constant \(c\), rationality of the bridge quadratic
in \(Z\) is equivalent to
\[
q^2=5+\frac{4(16-5\upsilon^2)}{5c}
=5+\frac{64}{5c}-\frac4c\upsilon^2,
\qquad Z=\frac2{q+1}.
\]
As a polynomial in \(\upsilon\), the right side is a square in
\(\mathbb Q[\upsilon]\) exactly when
\[
\boxed{c=-\frac{64}{25}},
\qquad
q=\pm\frac54\upsilon.
\]
Indeed the square of an affine polynomial has a linear term unless one
coefficient vanishes; the nonconstant coefficient here cannot vanish,
so the constant term must be zero.  This proves uniqueness inside the
constant-\(c\) ansatz.

For this value of \(c\), the constant-two equation is solved by every
\(\lambda\ne0\) with \(\lambda^2\ne5\) after setting
\[
\boxed{
b=\frac{2101}{25000}
-\frac{2\lambda^2}{(\lambda^2-5)^2}.}
\]
The two bridge maps are
\[
Z_+=\frac{8b}{5b^2-6b+5},
\qquad
Z_-=\frac{8b}{-5b^2+14b-5}.
\]
At \(\lambda=3\), choose the positive bridge root.  Then
\[
\boxed{
\left(a,b,Z,\lambda,c\right)
=\left(
1,-\frac{3253}{3125},
\frac{2033125}{6101423},
3,-\frac{64}{25}
\right)}
\]
and
\[
X=\frac73,\qquad \rho=\frac23,\qquad
H_2(9)=C_1=C_2=0.
\]
All bridge guards hold and \(a,b\) are \(2\)-adic units.  This is a
genuine rational point on the bridge-specialized quartic with free
rational \(Z\).

It is not an L30 point.  Exactly,
\[
Z=\frac{5^4\cdot3253}{1009\cdot6047},
\]
so \(Z\notin\mathbb Q^{\times3}\).  Moreover
\[
v_{3253}(b)=v_{3253}(Z)=1,\qquad v_{3253}(c)=0.
\]
Thus the prime shared by \(b\) and \(Z\) is not on L30's target stratum,
where \(Z=z^3\) and \(v_w(c)=6v_w(z)-2\ge4\).

Imposing \(Z=z^3\) on the two bridge maps requires, respectively,
\[
\boxed{d^2=-64z^6+96z^3+64}
\]
or
\[
\boxed{d^2=96z^6-224z^3+64.}
\]
They factor as
\[
-32(2z^3+1)(z^3-2),\qquad
32(3z^3-1)(z^3-2),
\]
so both degree-\(6\) polynomials are squarefree and the curves have
genus \(2\).  An external exact calculation with Magma V2.29-9 ran
`RationalPointsGenus2` on both curves.  In each case it returned
\[
\{(0:-8:1),(0:8:1)\}
\]
with completeness flag `true` and Jacobian rank bounds \([0,1]\).
Magma's intrinsic help defines that flag to mean that the returned set
contains all rational points; it documents two-cover descent,
rank-zero elliptic subcovers, finite-index Mordell--Weil computation,
Chabauty, and the Mordell--Weil sieve as the available proof methods.
Consequently both curves have only \(z=0\), which violates the bridge
guard.  The unique constant-\(c\) bridge-square section is therefore
closed.

This is a scoped external-CAS proof, recorded with exact request,
version, output, and help semantics in
`data/l31_magma_genus2.json`; it has no checked-in standalone Magma
certificate.  It is not a no-point theorem for the full L30 quartic.

### L31c. An exact globally good reciprocal member on the fixed field

Take
\[
\boxed{w=13,\quad a=3,\quad z=13,\quad Z=169,\quad q_0=1,}
\]
and put
\[
b=wq_0^2=13,\qquad
\rho=q_0+\frac1{wq_0}=\frac{14}{13},\qquad
\lambda=wq_0-\frac1{q_0}=12.
\]
Then
\[
u=b+b^{-1}=\frac{170}{13},\qquad
\lambda^2=w(w\rho^2-4)=144.
\]
Thus the mandatory reciprocal lift of L30g is satisfied explicitly; no
rational trace-base point is being assumed.

The canonical L26 data are
\[
\begin{aligned}
A&=37,\\
D&=-257217=-3\cdot83\cdot1033,\\
N_g&=-548208=-2^4\cdot3^6\cdot47,\\
c&=\frac{277941456}{3172343},\\
\eta&=-\frac{182}{257217},\\
\delta&=-\frac{324}{37}.
\end{aligned}
\]
The tied norm value is
\[
\boxed{
M=\frac{225318830729979494224}{3351232116513117}
=\frac{2^4Q}{3^2\,37^3\,83^2\,1033^2},}
\]
where
\[
Q=14082426920623718389
\]
is prime.  A recursive Pocklington certificate uses
\[
\begin{aligned}
Q-1&=2^2\cdot3\cdot1173535576718643199,\\
1173535576718643199-1
&=2\cdot3^2\cdot10601\cdot6150025556911,\\
6150025556911-1
&=2\cdot3^3\cdot5\cdot17\cdot29\cdot46202581,\\
46202581-1&=2^2\cdot3^2\cdot5\cdot283\cdot907.
\end{aligned}
\]
The first three Pocklington bases are \(2,3,3\); every displayed
certified factor exceeds the square root of its parent.  The producer
checks each Fermat congruence and gcd condition and independently
replays the repository's proven-primality engine.

Every local symbol is \(+1\):
\[
\begin{array}{c|c|l}
v&(M,26)_v&\text{reason}\\ \hline
2&+1&M/16\equiv1\pmod8\\
3,83,1033&+1&v(M)\text{ is even}\\
13&+1&M\equiv3=4^2\pmod{13}\\
37&+1&(26\mid37)=+1\\
Q&+1&(26\mid Q)=(2\mid Q)(13\mid Q)=(-1)(-1)\\
\infty&+1&M>0 .
\end{array}
\]
There are no other possible ramified places.  Hence
\[
\operatorname{Ram}(M,26)=\varnothing.
\]
Hasse--Minkowski gives a rational point on
\[
U^2-26V^2=MW^2.
\]
Here \(W\ne0\), since \(26\) is not a rational square.  Therefore this
is an exact globally good reciprocal member on the fixed-field slice.
Equivalently,
\[
T(170/13)=M(DA)^2=\frac{2^4Q}{37}.
\]

**Strict scope.**  This proves one member at \(w=13\).  It does not
prove a member for every target.  For \(q_0=r/s\), the general slice is
the sign-decorated degree-\(16\) binary form
\[
G_w(r,s)=s^{16}P_{\rm rec}(wr^2/s^2),
\]
with every \((G_w(r,s),2w)_v\) required to be \(+1\).  Primes belonging
only to \(q_0\) enter \(M\) with even valuation away from fixed support
and therefore cannot cancel emergent bad places.  If an odd prime
\(p\mid A\) occurs to odd multiplicity under the standard unit guards,
the fixed local requirement is \((2w\mid p)=+1\).  Thus the general
problem remains a degree-\(16\) member theorem, not a consequence of
this row.

### L31d. The exact two-large-divisor dispersion boundary

For a selected L20/L22 linear--octic pair, and every good prime \(p\),
put
\[
\mathcal C_p^-=
\left\{r\bmod p:G(r)=0,\quad
\left(\frac{u(r)}p\right)=-1\right\}.
\]
For distinct \(p,q\), the simultaneous bad-root classes are exactly
\[
\mathcal C_{pq}^-=
\operatorname{CRT}(\mathcal C_p^-\times\mathcal C_q^-).
\]
The raw von Mangoldt pair sum is therefore
\[
\mathcal S_{2,\Lambda}^{\rm root}
=
\sum_{\substack{p<q,\ p,q\ge z\\D_{\rm BV}<pq\le X^{8+o(1)}}}
\ \sum_{c\in\mathcal C_{pq}^-}
\ \sum_{\substack{X<t\le2X\\t\equiv c\ (pq)}}\Lambda(Q(t)).
\]
Writing the inner progression count as its expected term plus
\(E_{pq,c}\), after dyadic truncation in this finite range, exhibits
the exact problem: every modulus in the first pair range satisfies
\[
pq\ge z^2=X^{0.98}>D_{\rm BV}
=\frac{X^{1/2}}{(\log X)^B}.
\]
The mod-\(p\) condition is an oversieve.  For a simple root, exact
\(v_p(G(t))=1\) occupies \(p-1\) classes modulo \(p^2\).  Pairing two
such conditions uses modulus \(p^2q^2\), reducing the plain BV product
range to \(pq\le X^{1/4-o(1)}\).

To distinguish counting scales, let \(\mathcal A(X)\) be L23a's
unweighted set of \(Q(t)\)-prime members with no bad root prime below
\(z\), and define
\[
\mathcal B_2^{\rm odd}(X)=
\#\left\{t\in\mathcal A(X):
\begin{array}{l}
\text{there exist distinct }p,q\ge z,\\
v_p(G(t)),v_q(G(t))\text{ odd, and both signs bad}
\end{array}\right\}.
\]
This is a prime count.  Its von Mangoldt-weighted analogue has scale
\(\log X\) times larger; neither is the raw mod-\(p\) oversieve
\(\mathcal S_{2,\Lambda}^{\rm root}\) above.

There is an additional load-bearing warning.  If \(R_p\) is the number
of roots and \(C_p\) their signed character sum, then
\[
\#(\mathcal C_p^-\times\mathcal C_q^-)
=\frac14(R_pR_q-R_pC_q-C_pR_q+C_pC_q).
\]
Cancellation of signed character sums leaves the positive
\(R_pR_q/4\) term.  L23a's Chebotarev law and partial summation give
\[
\sum_{X^{0.49}\le p\le X^8}g_{\rm root}(p)
=\frac12\log\!\left(\frac8{0.49}\right)+o(1).
\]
Thus the local mass is a positive constant, not \(o(1)\), and the raw
pair main is of order \(X\) in von Mangoldt weight
(\(X/\log X\) in prime count).  It cannot be discarded by character
cancellation or by an AP-error estimate alone.  The natural missing
theorem is a
fixed-family Buchstab or analytic-Hilbert-detector asymptotic with
positive zero-bad main term
\[
\boxed{\mathcal N_{\rm good}(X)
\sim C_{a,Z}\frac{X}{(\log X)^{3/2}},\qquad C_{a,Z}>0.}
\]
The two-large sector is the first term of that decomposition beyond the
available distribution level.  The stronger prime-count bound
\[
\mathcal B_2^{\rm odd}(X)
=o\!\left(\frac{X}{(\log X)^{3/2}}\right)
\]
would suffice but is neither proved nor the natural expected
asymptotic.  Its von Mangoldt-weighted analogue is
\(o(X/\sqrt{\log X})\).

Bombieri--Vinogradov, the beta sieve, elementary polynomial congruence
counts, and the cited almost-prime theorems do not supply this
fixed-family asymptotic.  Frei--Sofos gives an \(L^2\) coefficient-box
theorem, not the required pointwise result for this octic on a fixed
prime progression.  Hence the requested dispersion estimate remains
open.

### L31e. AP1 is the exact parity replacement for Schinzel H

Fix one selected aligned class and a \(Q(t)\)-prime member.  Let
\[
\mathcal B(t)=
\left\{\ell\notin S\cup\{Q(t)\}:
v_\ell(P(b(t)))\text{ is odd and }(x_0,d_0)_\ell=-1
\right\}
\]
and put \(R_{\rm bad}(t)=\#\mathcal B(t)\).  Alignment makes every
symbol at \(S\cup\{Q(t),\infty\}\) equal to \(+1\); unsupported even
valuations also give \(+1\).  Hilbert reciprocity therefore gives
\[
\boxed{R_{\rm bad}(t)\equiv0\pmod2.}
\]
Define:
\[
\boxed{\mathrm{AP1}:\
\text{for every cell, one selected aligned class has a \(Q\)-prime
member with }R_{\rm bad}\le1.}
\]
Inside the L19--L22 protocol,
\[
\boxed{\mathrm{AP1}\iff\text{intermediate per-cell H}.}
\]
The forward implication is parity: \(R_{\rm bad}\le1\) and even forces
\(R_{\rm bad}=0\).  The reverse implication is immediate.  This is the
weakest threshold available from parity; \(P_2\) is insufficient.
Indeed, at \(t=0\) the admissible tuple
\[
(8t+5,8t+7,8t+23)
\]
has individual signs
\[
\left(\frac{10}{7}\right)
=\left(\frac{10}{23}\right)=-1
\]
but product \(+1\).

L23a supplies many members with no bad root prime below \(X^{0.49}\),
but L23b permits
\[
R_{\rm bad}\in\{0,2,\ldots,16\}.
\]
Therefore L23 does not imply AP1.  If a uniform AP1 theorem were proved,
the chain would become
\[
\mathrm{AP1}\Longrightarrow H\Longrightarrow\text{Theorem C}
\]
without classical Schinzel H.  AP1 is not proved here, so the actual
chain remains
\[
\boxed{\text{classical Schinzel H}\Longrightarrow
H\Longrightarrow\text{Theorem C}.}
\]

### L31f. Strict frontier

The four requested outcomes now have exact statuses:

1. a rational root of \(H_2(\lambda^2)=0\) is proved immediately before
   the cube pullback; the unique constant-\(c\) bridge-square section is
   externally proved to have only \(z=0\) after imposing the cube, and
   no root satisfying \(Z=z^3\) and every controlled place is known;
2. one globally good fixed-field reciprocal member is proved at
   \(w=13\), but no uniform-in-\(w\) member theorem is known;
3. the two-large-divisor problem is reduced to the displayed
   fixed-family Buchstab/detector asymptotic beyond BV, which remains
   unproved;
4. AP1 is proved equivalent to the intermediate member hypothesis, but
   AP1 itself remains unproved, so classical Schinzel H is not removed.

No unconditional quantifier record or conclusion about
H10/\(\mathbb Q\) changes
(`l31_frontier_push.py`; `data/l31_frontier_push.jsonl`;
`data/l31_magma_genus2.json`).
