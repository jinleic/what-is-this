# Hilbert's Tenth Problem over Q — problem map (2026-08)

**Problem.** Is there an algorithm that, given $f \in \mathbb{Q}[x_1,\dots,x_n]$, decides whether
$f = 0$ has a solution in $\mathbb{Q}^n$? Equivalently: is the existential (positive) theory of the
field $\mathbb{Q}$ decidable? Hilbert's original Problem 10 (over $\mathbb{Z}$) was answered
negatively by Matiyasevich–Robinson–Davis–Putnam (1970). The rational version is the live
descendant and is **open in both directions**.

## Status ladder

| Ring | Status | Source |
|---|---|---|
| $\mathbb{Z}$ | undecidable (MRDP) | 1970 |
| $\mathbb{Z}[S^{-1}]$, some recursive $S$ of prime density 1 | undecidable | Poonen, arXiv math/0306277, Thm 1.3 |
| complementary partitions $S_1 \sqcup S_2$ of all primes | undecidable on both sides | Eisenträger–Everest–Shlapentokh, arXiv 1012.4878 |
| $\mathcal{O}_K$, all number fields $K$ (⇒ all infinite f.g. rings) | undecidable, unconditional | Koymans–Pagano, arXiv 2412.01768 (Dec 2024); independently Alpöge–Bhargava–Ho–Shnidman, arXiv 2501.18774 (Jan 2025) |
| $\mathbb{Q}$ | **open** | Koymans–Pagano expository, arXiv 2602.04468, Question 6.2 |

Key asymmetry: undecidability holds on rings "arbitrarily close" to $\mathbb{Q}$ (invert a
density-1 set of primes), yet the methods collapse exactly at $\mathbb{Q}$.

## Pillar 1 — definability of Z in Q (the direct route to undecidability)

If $\mathbb{Z}$ were existentially (Diophantine) definable in $\mathbb{Q}$, MRDP would transfer and
H10/$\mathbb{Q}$ would be undecidable. What is actually known:

- **Poonen, arXiv math/0703907.** $\mathbb{Z}$ is $\forall\exists$-definable in $\mathbb{Q}$ with
  2 universal + 7 existential quantifiers. Machinery: for $a,b \in \mathbb{Q}^\times$, the
  quaternion algebra $H_{a,b}$, its ramification set $\Delta_{a,b}$, and
  $S_{a,b} = \{2x_1 : \exists x_2,x_3,x_4,\ x_1^2 - a x_2^2 - b x_3^2 + a b x_4^2 = 1\}$
  (reduced traces of norm-1 elements).
  - Lemma 2.1: $p \notin \Delta \Rightarrow S(\mathbb{Q}_p) = \mathbb{Q}_p$; $p \in \Delta
    \Rightarrow \mathrm{red}_p^{-1}(U_p) \subseteq S(\mathbb{Q}_p) \subseteq \mathbb{Z}_p$, where
    $U_q = \{s \in \mathbb{F}_q : x^2 - sx + 1 \text{ irreducible}\}$.
  - Lemma 2.3: $U_q \ne \emptyset$; $q > 11 \Rightarrow U_q + U_q = \mathbb{F}_q$ (hence the
    constant $N = 2\cdot3\cdot5\cdot7\cdot11 = 2310$).
  - Lemma 2.4: for $a$ or $b > 0$: $T_{a,b} := S_{a,b} + S_{a,b} + \{0,\dots,2309\}
    = \bigcap_{p \in \Delta_{a,b}} \mathbb{Z}_{(p)}$.
  - Lemma 2.6: $\bigcap_{a,b > 0} T_{a,b} = \mathbb{Z}$ ($p=2$ via $(7,7)$; odd $p$ via
    $(p, b)$ with $b$ a positive nonsquare mod $p$).
  - Thm 4.1: explicit formula with parameterization $a' = a^2+b^2+1$, $b' = a^2+a+1+b^2$
    (Lemma 4.2), quantifiers $\forall a,b\, \exists x_1..x_4,y_2,y_3,y_4$:
    $Q_1^2 + \prod_{n=0}^{2309} Q_n^2 = 0$ with
    $Q_1 = x_1^2 - a'x_2^2 - b'x_3^2 + a'b'x_4^2 - 1$,
    $Q_n = (n - t - 2x_1)^2 - 4a'y_2^2 - 4b'y_3^2 + 4a'b'y_4^2 - 4$.
- **Koenigsmann, arXiv 1011.3424 (Annals 2016).**
  - Thm 1: $\mathbb{Z}$ is *purely universally* definable in $\mathbb{Q}$ ($\Pi_1$).
  - Cor 2: $\mathbb{Q} \setminus \mathbb{Z}$ is Diophantine in $\mathbb{Q}$.
  - Cor 21: a $\forall\exists$ definition with just **one** universal quantifier.
  - Cor 3: the $\forall\exists$-theory of $\mathbb{Q}$ is undecidable.
  - Cor 23: under a strong Bombieri–Lang variant, **no** existential definition of $\mathbb{Z}$
    in $\mathbb{Q}$ exists.
  - Architecture: Poonen's semilocal gadget → uniform Diophantine $\mathbb{Z}_{(p)}$ via quadratic
    reciprocity → existential definition of Jacobson radicals of semilocal subrings → Hilbert
    reciprocity supplies a missing ramified place → complementation flips $\Sigma$ to $\Pi$.
- **Quantifier record (universal definitions of $\mathbb{Z}$ in $\mathbb{Q}$):**
  Koenigsmann preprint 418 → Daans 2021: 37 → Sun–Zhang 2021: 32 → Daans 2023 (JLMS 2024): 10
  (refereed) → **Sun, arXiv 2607.28606 (2026-07-30): 7** — unrefereed; completeness chain
  audited 2026-08-12: traceable in the posted §§3–8, no substantive gap found; the "Section 11"
  / "Sections 7 and 8" citations at the decisive step are stale draft numbering (editorial),
  intended §§5–6 and §8 (audit: THEOREMS.md A2; machine-checked parts: h10q.py `_verify_sun`).
  Daans' count assembly: $2+3+6-1=10$ (Thm 5.6; fusion via DDF21 Thm 1.4, non-constructive).
  Sun's claimed assembly: $2+(3+3-1)=7$ — Daans' $\exists_6$ intersection block replaced by
  $\Theta = \bigvee_{\tau\in\Lambda}\Psi_\tau$, a finite frozen-$\tau$ disjunction of
  3-witness trace-sum equations, soundness via norm-one quaternion traces (verified), 
  completeness via §5 target-place Hensel + §§3–4 frozen smooth points + §7 parameter
  selection + §8 Hasse–Minkowski/Lemma 8.1 (traced at statement level; refereeing pending).
    Session contributions to this record (2026-08-12, THEOREMS.md L4/L5/A2): Sun's Lemma 5.1
    proved as an exact identity — solution count $= \#E_\lambda(\mathbb{F}_q)/4$, $\lambda =
    \tau^{-2}$ Legendre — so by Hasse his threshold $|k| > 25$ is VOID for every odd prime
    power. The *Lemma-5.1 component* of $E_{\mathrm{exc}}$ shrinks from $\{|\kappa_w| \le 25\}$
    to $\{|\kappa_w| = 3\}$ ($\tau$-vacuity); his data-dependent exclusions ($\pi$-support,
    $\Lambda$-poles, parameter-difference zeros) remain. Explicit witnesses for every odd
    prime target $w < 300$ under a proven-primality engine (8 targets need
    $\tau \notin \{0, \tau_1\}$); the $s{=}0$ ternary specialization certified sound-but-incomplete.
  - **L6 witness-tie architecture (2026-08-13; THEOREMS.md L6):** identify Sun's
    block witness with Daans' congruence parameter, $a=1+2s$. Then the base has
    $(b,s)$ and the tied block only $(y,r)$; DDF Thm 1.4 gives
    $2+(3+2-1)=6$. The count survived an independent adversarial audit and a
    direct recheck of DDF v5 Thm 1.4. Soundness is inherited.
    - **W0–W1, proved:** the tied target conic is soluble iff
      $\chi_w(-\delta_\tau A)=1$.  The character identity
      $\sum_a\chi(1+4a^2)=-1$ selects a nonsquare $A$, and the two
      selector branches $\tau=0$ and $\tau=2a/A$ cover every odd residue
      field, including $\mathbb F_3$ (`THEOREMS.md`, L6).
      Theorem C's fixed finite menu also contains
      $\tau^\dagger=(1+2a^2)/A$, $\delta=-4a^4/A$; L11a proves that
      adjoining it costs no witness (`THEOREMS.md`, L11a and L22).
    - **W2, proved:** its exact global criterion is
      $M_\tau=0$ or
      $-\delta_\tau A M_\tau\in
      N_{\mathbb Q(\sqrt{-\delta_\tau AB})/\mathbb Q}(E_\tau^\times)$,
      with the zero case represented by $(y,r)=(0,0)$
      (`THEOREMS.md`, L6/W2).
    - **Historical L6 endpoint — SUPERSEDED by L19–L22:** choosing
      global $(s,b)$ was originally OPEN, and the $61$ guarded examples
      below were EVIDENCE only (`data/l6_witnesses.jsonl`).  L20 now
      supplies an aligned class for every cell, L22 supplies vertical
      irreducibility, and L19 supplies a member under classical Schinzel
      H (`THEOREMS.md`, L19–L22).
    Thus **six is now an established CONDITIONAL count under classical
    Schinzel H alone**, not an unconditional record
    (`CONDITIONAL.md`, §§1–2).
    - **L7 structure (2026-08-15; THEOREMS.md L7, corollaries + bounded probe):**
      at odd places with $v(A)=v(B)=0$ the tie's exact cost is the wild
      (odd-valuation) places of $M_\tau$, where the untied quaternary is
      universal (L7a); the unit hypothesis is necessary and coefficient-bad
      places ($v(A)\ne0$ or $v(B)\ne0$) carry real obstructions
      (machine-witnessed thrice, incl. a cancellation case with
      $v_5(AB)=0$), so the obstruction support is
      $\{2,\infty\}\cup\mathrm{wild}(M_\tau)\cup\{v:v(A)\ne0\text{ or
      }v(B)\ne0\}$; no fixed $(s,b,\tau)$ has a rational-function witness
      $Y,R\in\mathbb{Q}(z)$ (L7b, conditional on L6 soundness; non-target
      probe of 8 canonical witnesses $w\le23$, 224 cells total: 0 solutions
      found, 205 certified failures, 19 refusals), and probed canonical pairs
      certifiably fail at some $z\in\mathfrak m_w$, so they do not cover
      pointwise (some other fixed pair remains open); obstruction sets
      are always even, so repairs flip places in pairs (L7c). Free-$\Phi$
      switching probe: 93/93 target cells and 71/71 recorded rescues
      were positive (`THEOREMS.md`, L7).  This is a historical bounded
      layer, superseded as the global frontier by L19–L22.

  Minimal count posed as Daans Question 5.7: the unconditional
  unrefereed range is $2\le m\le7$ and the refereed upper bound is $10$;
  under classical Schinzel H, L22 narrows it to $2\le m\le6$
  (`THEOREMS.md`, consequences table and L22d).  The lower bound is
  Daans Theorem 2.6.  Equivalently, the conditional essential-dimension
  range is $[1,5]$ (`THEOREMS.md`, A3; `CONDITIONAL.md`, §1).

- **Lower bounds stop at one witness.** DDF21 §8: no known example of an $\exists$-definable
  subset of a global field provably not $\exists_2$ (best positive case: sums of two squares
  need exactly 2, DDF Cor 8.11). $\exists$-definability of $\mathbb{Z}$ in $\mathbb{Q}$ is not
  excluded unconditionally. See THEOREMS.md (A1).

## Pillar 2 — the elliptic route (how O_K fell, and why Q resists)

- **Poonen's criterion** (arXiv 2602.04468, Thm 5.2; orig. Poonen 2002, refined by Shlapentokh to
  rank $> 0$): if $R$ is a maximal order with $K = \mathrm{Frac}\,R$ and $E/\mathbb{Q}$ is an
  elliptic curve with $\mathrm{rank}\,E(\mathbb{Q}) = \mathrm{rank}\,E(K) = 1$, then $\mathbb{Z}$
  is $R$-Diophantine, hence H10($R$) is undecidable.
- **Rank stability, now unconditional:**
  - Koymans–Pagano (arXiv 2412.01768): curves with full rational 2-torsion; 2-Selmer control in
    quadratic extensions via 2-descent + additive combinatorics (simultaneous prime values of four
    linear forms); no rank growth in the needed quadratic steps.
  - ABHS (arXiv 2501.18774, Thm 1.1): for *every* quadratic extension $K/F$ of number fields there
    is an abelian variety $A/F$ with $\mathrm{rank}\,A(F) = \mathrm{rank}\,A(K) > 0$; via
    $(1-\zeta)$-Selmer groups of superelliptic Jacobians $y^2 = x^\ell + n$ and unit-equation
    combinatorics.
- **The integrality mechanism** (Poonen, arXiv math/0306277): fix rank-1 $E/\mathbb{Q}$,
  generator $P$; write $x(nP) = A_n / B_n^2$. Facts (all verified numerically in `h10q.py`):
  - $\{n : r \mid B_n\} = n_r\mathbb{Z}$ (rank of apparition);
  - $v_r(B_{k n_r}) = v_r(B_{n_r}) + v_r(k)$ for odd $r$ (formal group);
  - $\log B_n = (c - o(1))\,n^2$, $c$ related to the canonical height $\hat h(P)$;
  - elliptic Zsigmondy: fresh (primitive) prime divisors at every index.
  Poonen inverts a density-1 prime set $S$ engineered so that $\{n : nP \text{ is }
  \mathbb{Z}[S^{-1}]\text{-integral}\}$ is essentially the primes; Vinogradov equidistribution of
  prime multiples in $E(\mathbb{R})$ then yields a Diophantine *model* of $(\mathbb{Z},+,\cdot)$,
  discrete in the archimedean topology.

**Why it stops at $\mathbb{Q}$:** the predicate "$nP$ is $S$-integral" degenerates as
$S \to$ all primes — over $\mathbb{Q}$ every point is trivially integral and there is no valuation
left to carve anything out. Any variant that only omits finitely many primes yields integral index
sets of positive density, whose points equidistribute in $E(\mathbb{R})$ — hence *not* discrete, so
the model construction fails structurally, matching:

- **Mazur's conjecture** (1992): for any variety $X/\mathbb{Q}$, the closure of $X(\mathbb{Q})$ in
  $X(\mathbb{R})$ has finitely many connected components. Consequences:
  - no infinite discrete Diophantine subset of $\mathbb{Q}$;
  - **Cornelissen–Zahidi** (arXiv math/0006140): Mazur ⇒ there is **no Diophantine model of
    $\mathbb{Z}$ over $\mathbb{Q}$ at all**.
- **Koenigsmann Cor 23:** strong Bombieri–Lang ⇒ no existential definition of $\mathbb{Z}$.

**Crucial nuance:** "no Diophantine model of $\mathbb{Z}$" does **not** imply H10/$\mathbb{Q}$ is
decidable — undecidability could enter without interpreting $\mathbb{Z}$. Conversely no one has any
positive algorithm even for 2 variables: deciding rational points on plane curves already needs an
effective Mordell/Faltings or BSD-type input (genus 1). Both directions are genuinely open.

## The fork (the "clear path" and its walls)

1. **Undecidable side:** produce a Diophantine model of $\mathbb{Z}$ in $\mathbb{Q}$ — must defeat
   Mazur + strong BL (real-topology and sparsity obstructions). All known encodings (discrete sets
   via integrality on rank-1 curves) die exactly at $\mathbb{Q}$.
2. **Structural side:** prove Mazur-type statements (kills route 1 but does not decide H10/Q).
3. **Quantitative side (tractable):** shrink the definability gap. Established records:
   $\forall^{10}$ (Daans, refereed), claimed/audited $\forall^7$ (Sun, unrefereed);
   L6 gives an audited $\forall^6$ *architecture* with one open assembly lemma.
   Stratify existential definability by witness count (see RESULTS.md).

## Session 2026-08-15 — L8: the alignment wall (derivation log)

Attack sequence and falsified hypotheses, kept for provenance:

1. **Value identities.** With $X=a^2z^6N_g$, $Y=4Ab^2D_z$, $P=1-ABs^2$:
   $(Ab^2D_z)^2u_0=X^2-PY^2$ and $(Ab^2D_z)^2u_1=X^2-APY^2$, where
   $u_0=-M_0$, $u_1=-AM_1$ exactly. W2 becomes $A u_\tau\in N(E_\tau)$ (or
   $M_\tau=0$), $E_0=\mathbb Q(\sqrt{-2Ab})$, $E_1=\mathbb Q(\sqrt{-2b})$.
   Both branches: $D_\tau\operatorname{disc}E_\tau\equiv-2AbP$.
2. **Falsified: fixed-gauge slices.** $b=\mp w$ with the parity-matched
   branch (the pattern of all 10 frozen rescues, which use $|m|=1$) fails
   for most $(u,s)$: obstruction profiles show huge inert odd-multiplicity
   primes of the value (up to $10^{21}$) — wild, not frame, places.
3. **Falsified: Pell steering with $t\ne\square$.** Gauge
   $-2AbP=t\rho^2$ is rationally parametrized
   ($b=-t\theta^2/(2As^2(1-t\theta^2))$ etc.); success rate rises
   ($\approx$4–8× baseline) but wild leaks persist: $v\nmid t$ does not give
   $\chi_v(t)=1$.
4. **Falsified: full collapse.** $t=\square$ plus
   $1-\theta^2=A^{\varepsilon}\kappa\,\square$ makes $E_\tau=\mathbb
   Q(\sqrt\kappa)$, $D_\tau\equiv\kappa$, the value literally a norm, and
   the branch condition collapses to $A\in N(\mathbb Q(\sqrt\kappa))$ —
   then $v_2(\kappa)$ odd (forced by $v_2(b)=0$) with $A\equiv5\ (8)$ gives
   $(A,\kappa)_2=-1$. Sharper: the collapse family never even reaches
   $v_2(b)=0$ (case check on $v_2(\theta)$), and the wall is elementary:
   admissibility forces $v_2(P)=0$, so $v_2(-2AbP)=1$ — THEOREMS L8b/L8c.
5. **Computed residue data** (sympy in `math/.venv`, logged here; not
   in-suite): over $\mathbb Q(b)$, $\operatorname{NUM}_\tau(b)$ is an
   irreducible octic for generic $(s,Z)$, the class $(Au_1,-2b)$ has residue
   $A$ at $b=0$, and $-2b$ is **not** a square in the octic residue field:
   finite-field splitting witnesses at $(s,Z)\in
   \{(1,3),(1,5),(2,3),(-1,7),(1,\tfrac13)\}$, $p\in\{101..137\}$ — mixed
   QR patterns (e.g. $s{=}1,Z{=}3,p{=}107$: factor degrees $(1,2,5)$ with
   QR $(F,F,T)$). So the slice bundle is genuinely ramified at the value
   divisor; any CTS/Schinzel attack must engage it.
6. **Scope discipline.** A general fixed-$j$ class-matching family is NOT
   ruled out (value-primes of a thin family need not equidistribute); L8c
   claims exactly the two identities $j\equiv\square$, $j\equiv D_\tau$ and
   no more (advisory upheld during the session).

7. **Witness-density sweep provenance (not persisted in-suite).** Grid:
   $w\in\{5,13\}$ with $b=-w$, branch $\tau=2a/A$, and $w\in\{7,11\}$ with
   $b=+w$, branch $\tau=0$; cells $u\in\{1,-2,3,\tfrac13,-\tfrac15,
   \tfrac23\}$; witnesses $s\in\{0,1,-1,2,-3,\tfrac13,\tfrac25\}$;
   budget-refused pairs excluded from both counts. Outcome: 13 of 161
   admissible $(u,s)$ pairs soluble on the parity-matched branch
   ($\approx8\%$); every failing pair's obstruction set had even size
   (sizes 2, 4, 6 observed — L7c live). This is the figure cited by
   THEOREMS L8d / RESULTS / README.

## Session 2026-08-16 — L9: reciprocity steering (derivation log)

1. **The lemma.** $\prod_{v\in T}(x,d)_v=1$ for
   $T=\{2,\infty\}\cup\operatorname{supp}x\cup\operatorname{supp}d$
   (Hilbert reciprocity + unit-unit symbols at odd $v\notin T$). Search
   consequence: candidates whose norm value $x=\alpha M_\tau$ resolves as
   (trial-smooth)$\times$(one proven prime power) need no hard factoring
   and no luck at the wild prime — its symbol is forced when the
   controlled places align. This retro-explains the single-wild L8d
   rescues; only multi-wild rescues carry alignment content beyond
   reciprocity.
2. **Blob extension (load-bearing at large $w$).** $x$'s denominator
   carries $D_z^2$-type squares whose cofactor roots are random ~$10^{20}$
   composites; requiring their primality killed every candidate at
   $w\approx89$. Fix: an unresolved composite cofactor $C^e$ with $e$
   EVEN and $\gcd(C,\operatorname{supp}d)=1$ is admissible — each of its
   primes has $v_q(x)$ even, $v_q(d)=0$, symbol $+1$ structurally.
   Certificates remain complete (exact reconstruction assert).
3. **Forcing vs checking (advisory upheld).** The verifier certifies
   'steered' rows from $T\setminus\{q_0\}$ and only *replays* the wild
   symbol as a consistency assert — the mechanism label matches what the
   code proves. An earlier draft checked the full support and mislabeled
   it as forcing.
4. **Prime-$b$ reduction (advisory upheld, then repaired).** A naive
   $b\equiv b_0\bmod N$ family lets uncontrolled primes into
   $d=\alpha B$. Repair: $b=\varepsilon q_1$ prime; then $d$'s new
   support is exactly $\{q_1\}$ and clearing denominators at
   $b\equiv0\bmod q_1$ gives numerator $\to-\alpha\delta a^4z^{12}A^2$,
   denominator valuation 4, so $v_{q_1}(x)=-4$ and unit part $\equiv A$:
   $(x,d)_{q_1}=(A|q_1)$, class-controlled by reciprocity;
   $A\equiv5\bmod8$ (admissible $a$ odd 2-adic unit) so $A$ is never a
   square and $+1$-classes exist. Machine-checked as L9a — but only
   under unit hypotheses imposed FACTOR BY FACTOR: the audit's P1
   showed $2\delta_\tau A=2$ on the branch $\tau=2a/A$, so a product-form
   guard admits $q_1\mid A$; at $w=3$, $z=6$, $q_1=b=5=A$ the conclusion
   fails outright ($v_5(x)=-5$ odd, symbol $-1$). Second audit P1: the
   advertised class modulus $p^{v_p(x(b_0))+1}$ over data primes does not
   freeze controlled symbols — cell $(73,5)$, $s=0$, $b_0=-29$,
   $\tau=2/5$, least $N=689120$, prime $q=29+355N$ in the same class and
   sign with $(A|q)=+1$, yet $v_{59}(x)$ moves $0\to1$ and the symbol
   flips $+1\to-1$. Corrected recipe: $p^{k_p}$ for every finite
   $p\in S^\dagger$ with $k_p$ from local constancy and no closed formula.
   A second draft took $S^\dagger=\{p\le H\}\cup\operatorname{supp}(\text{fixed
   data},d)$, $H=10^4$ (bound renamed from $B$, which is $2b$), on the
   worry that data primes dwarf $H$ — this cell has
   $D_z=-59\cdot40077920921911$. **L10-0 voided that worry**: those primes
   have $v_p(d)=0$, so they are wild candidates, not controlled places. The
   modulus freezes only $S_0=\{2\}\cup\operatorname{supp}(\alpha)$; the
   symbol set is $\Sigma=\{\infty\}\cup S_0\cup\{q_1\}$, and $q_1$ is held
   by L9a plus its residue class, never by $N$ ($\gcd(q_1,N)=1$).
   Both counterexamples and the
   $D_z$ factorization stay as in-suite regressions. Class EXISTENCE was
   unverified when this item was written; **L10 settled it** — proved on
   the 164 $w\equiv1\bmod4$ cells, proved impossible on $\tau=0$ for prime
   $A$. The engine's own hits remain pointwise (arbitrary admissible
   rational $b$) and are not class constructions.
5. **Run provenance.** Grid 353 cells: odd primes $w<100$ against a fixed
   15-value $u$-pool filtered per target by $v_w(u)\ge0$ (equivalently
   $v_w(wu)\ge1$) — 11 retained at $w=3$, 14 at $w=5$, 13 at $w=7$, 15 at
   the other 21 targets; $u=\pm w$ is retained, so this is not a
   $w$-unit pool. First pass: 283 covered (7305 s). Deep rescue of
   the 70 failures (6 parallel workers, widened pools, 420 s/cell):
   +62, total **345/353**; mechanisms 344 steered / 1 smooth; 51
   single-path (independent `_l7_tied_status` replay budget-refused), 294
   double-proof; largest forced wild prime 45 digits. Artifacts
   `data/l9_steer_run.jsonl`, `data/l9_rescue_{0..5}.jsonl`.
6. **Residual cells (8, named).** $(29,-2),(31,-2),(41,\tfrac17),
   (53,\tfrac13),(61,2),(61,-\tfrac27),(67,2),(89,-2)$: search-exhausted
   with 1.6k–4.1k certified-insoluble candidates each and no aligned one.
   Bad-place profiles over 150 sampled failures per cell are spread —
   e.g. $(29,-2)$: 2:88, 5:60, 3:55, 7:26, 37:25, 19:16, $\infty$:13 —
   no single dominating local wall; consistent with alignment depth.
   Next tool: CRT-targeted $b$-classes forcing $\{2,3,5,7\}$ jointly.

## Session 2026-08-16 — L10: class-side dichotomy (derivation log)
1. **The frozen set collapses.** $d=2\alpha b$, so
   $\operatorname{supp}(d)=\{2\}\cup\operatorname{supp}(\alpha)\cup\{q_1\}$.
   At any other odd $p$, $v_p(d)=0$ and $(x,d)_p=(d|p)^{v_p(x)}$, which is
   $+1$ for even $v_p(x)$. So the third advisory's worry — fixed primes
   above the trial bound, e.g. $40077920921911\mid D_z$ — dissolves: those
   are wild candidates, never controlled places. The over-large recipe
   $\{p\le H\}\cup$ data-support is superseded, not patched.
2. **Value model.** $x=\alpha P(b)/(D_z^2A^2b^4)$, $d=2\alpha b$ with
   $P(b)=16D_z^2A^2b^4-\delta_\tau a^4Z^4N_g^2-32A^3s^2D_z^2b^5$ (degree 8
   through $N_g^2$ — the octic value divisor of L8d). $D_z^2A^2b^4$ is a
   square, so $[x]_p=[\alpha P(b)]_p$.
3. **Exponents (fourth advisory).** The first `_l10_taylor` did synthetic
   division wrong and returned $P$'s constant coefficient instead of
   $P(b_0)$ — i.e. it expanded about $0$. Every $k_p$, hence every modulus,
   was void. Replaced by the exact binomial shift
   $c_j=\sum_{i\ge j}P_i\binom{i}{j}b_0^{i-j}$ with an in-suite invariant
   ($\sum_jc_jh^j=P(b_0+h)$ at four offsets). Rebuilt moduli came out
   *smaller*: $\{40,640\}$ where the void ones gave $\{40,10240\}$.
4. **Fifth advisory: L9a units.** On $\alpha=-1$ the guard
   $q_1\notin\operatorname{supp}(2\alpha)$ is vacuous, so
   $(x,d)_{q_1}=(A|q_1)$ could be asserted with $q_1\mid z$ or
   $q_1\mid D_z$. `_l10_class_cert` now rejects $v_{q_1}\ne0$ on each of
   $\delta_\tau,A,a,z,D_z$ separately — the same factor-by-factor
   discipline L9a needed.
5. **The dichotomy.** $\tau=2a/A\Rightarrow\delta_\tau=1/A\Rightarrow
   \alpha=-1$: no place divides $\alpha$, three symbols to align, and 164 classes were found — exactly the canonical
   $w\equiv1\bmod4$ cell set. $\tau=0\Rightarrow\alpha=-A$: at $p=A$ prime,
   $v_A(x)=-1$, $v_A(d)=1$, $u_x\equiv(Ac)^2$ square,
   $\epsilon(A)=2a^2$ even (needs $a$ ODD — an audit P1: $a=2$ gives
   $A=17\equiv1\bmod8$ and the identity yields $+1$), so
   $(x,d)_A=(-2\varepsilon|A)(q_1|A)$ and
   $(A|q_1)=(q_1|A)$; the product is $(-2\varepsilon|A)=-1$ for both signs
   since $A\equiv5\bmod8$. Anti-correlation, hence no aligned class.
6. **What is NOT proved (sixth advisory).** Rational and non-squarefree
   $A$: identity re-verified in-suite each run (200/600 instances), Jacobi
   cross-terms not handled. And
   no statement is made for all admissible $s$ — the scan is finite
   (11 $s$-values $\times$ 2 signs $\times$ 300 primes per cell).
7. **Consequence (scoped).** W0 ties $\tau$ to $w\bmod4$, so the reduction
   is established on the 164 $w\equiv1\bmod4$ cells. On the 189 others it
   is dead for prime $A$ and merely unobserved otherwise — the scan is
   finite and $s$ is not, so no universal claim. They keep
   their L9-E pointwise certificates; only the class-based *reduction* is
   affected. Next family to try: $b=\varepsilon q_1q_2$ (two primes), or
   a $b$ carrying a fixed square factor to move $\operatorname{supp}(\alpha)$
   off the anti-correlated place.

## References
- Poonen, math/0703907; math/0306277. Koenigsmann, 1011.3424. Daans, 2301.02107.
- Daans–Dittmann–Fehm, 2102.06941v5 (Thm 1.4, Cor 5.11; witness fusion).
- Sun, 2607.28606. Koymans–Pagano, 2412.01768; expository 2602.04468. ABHS, 2501.18774.
- Cornelissen–Zahidi, math/0006140. Mazur, Astérisque 228 (1995).
- Eisenträger–Everest–Shlapentokh, 1012.4878. "As easy as Q", 1601.07158.
- Sun–Zhang 2021 (32 quantifiers); Z[i] with 20 unknowns, 2510.18794.

## Session 2026-08-16 - L11: branch completion, a false claim, and the wall that is really there

1. **The hypothesis.** L10 ended with "those cells need a different witness
   family". I guessed the deficient object was the *branch set*: `_verify_sun`
   checks Sun's identities (2.2)-(2.4) at a **random** tau and Prop 2.1 never
   mentions tau, so soundness is tau-uniform and a new branch carries no
   soundness debt; and a disjunction of conics in the same $(y,r)$ is their
   product, so branches are free in the count too. Both parts are true (L11a).
2. **The branch.** $\alpha=-\delta_\tau A$ with $\delta_\tau=1-A\tau^2$. Forcing
   $v_A(\alpha)$ even needs $A\mid\operatorname{den}(\tau)$; with
   $\tau=m/(An')$, $\alpha\equiv m^2-An'^2$, the norm form of
   $\mathbb Q(\sqrt A)$. Taking $m=(A+1)/2$, $n'=1$ gives $(2a^2)^2$ - a perfect
   square. Generally $\alpha=\lambda^2$ iff $A+\lambda^2=\mu^2$, so
   $\tau_d=(A+d^2)/(2dA)$, $\lambda_d=(A-d^2)/(2d)$ for free rational $d\ne0$.
3. **The false step.** With supp$(\alpha)$ trivial and $\alpha>0$ the certificate
   reported frozen set $\{2\}$ and a free real place, and the search returned
   353/353 cells in about a second, uniformly at $s=0$. I checked the symbols
   independently over *all places of $d$* and $\infty$ - 0 failures - and froze
   it. **That check was the bug.** L10-0 frees an off-set place only where
   $v_p(x)$ is EVEN, and I never verified that hypothesis; the places of $x$
   were never examined at all.
4. **What the audit found.** $\delta$ carries $A$ in its denominator on every
   $\tau\ne0$ branch ($1/A$ canonical, $-4a^4/A$ square), so $v_p(x)$ is odd at
   $p\mid A$. Exact regression: cell $(3,1)$, $a=1$, $A=5$, $\tau^\dagger=3/5$,
   $q_1=41$: $c_0=9311592816/6345775$, $d_0=328$, $v_5(x_0)=-5$,
   $(x_0,d_0)_5=-1$ while $2$, $41$, $\infty$ are all $+1$. **293 of 353 L11
   rows were bad; the same defect invalidates 130 of L10c's 164 rows.**
5. **Second omission, also from the audit.** Prime-restricted fixed divisors:
   for $z=5$, $q\equiv41\bmod160$, $3\mid P(q)$ at *every* prime $q$ because $P$
   vanishes on $\mathbb F_3^\times$, though $3\nmid\operatorname{cont}(P)$. This
   is bounded: $\deg P=8$, so such a $p$ needs $p-1\le8$. Hence the complete
   controlled set is $\{2,3,5,7\}\cup\operatorname{supp}(\alpha)\cup
   \operatorname{supp}(\delta)$ (L11f), with the audit's content argument -
   term contents $2v_p(D_z)$, $4v_p(a)+4v_p(Z)$ (attained uniquely at the
   constant coefficient), $2v_p(s)+2v_p(D_z)$, all even - proving sufficiency.
6. **The real theorem (L11g).** The defect is not a bug, it is a wall, and it is
   **branch-independent**: $(x,d)_A\cdot(x,d)_q=-1$ whenever $v_A(z)\le0$, on all
   four branches (608 instances in-suite). Proof for $a=1$: $[x]=[5V]$ with
   $V=(Z^2N_g)^2+5(10q^2D_z)^2$ a square in $\mathbb Q_5$, so $[x]_5=[5]$ and
   $(x,d)_5=(2q\mid5)=-(q\mid5)$; L9a gives $(x,d)_q=(5\mid q)=(q\mid5)$; the
   product is $-1$. The common cause with L10b is admissibility:
   $A\equiv5\bmod8$ forces $(2\mid A)=-1$. **Branch completion buys nothing on
   the class side.**
7. **The characterization (L11h).** Escaping needs $v_p(z)>0$ for some $p\mid A$,
   and every odd prime dividing $A=1+4a^2=(n^2+4m^2)/n^2$ is $\equiv1\bmod4$;
   conversely for any such $p$, solving $n^2\equiv-4m^2\bmod p$ with odd
   representatives gives an admissible $a$ with $v_p(A)>0$. So: **reachable iff
   $z$ has a numerator prime $\equiv1\bmod4$ - 190 of 353 cells**, all 190
   certified by targeted construction. Necessity is proved only for $A$ prime
   (L11g's hypothesis); for composite or rational $A$ it is evidence - 878,400
   certificates over the first 10 walled cells, 244 admissible $a$ of which 228
   give composite/rational $A$, all four branches, 0 hits. Same scope split
   L10b already carries. This replaces L10c's 164 and L11e's 353.
8. **Where branch completion did pay.** The assembly (not class) side: the
   square-branch family closed the residual cell $(29,-2)$ at $s=-8$, $d=1/3$,
   $b=-29/5$, verified three ways (steered certificate, independent
   Hasse-Minkowski, symbols over the full support of $x$ and $d$) - no $(s,b)$
   reached it on the two original branches. With $(53,\tfrac13)$ closed by a
   wider $s$-pool, steered coverage went $345\to347$ of 353.
9. **A third hypothesis miss, caught by my own suite.** The extended run failed
   at cell $(59,-3)$ on $\tau=0$: the wall check had looped $q\in\{41,59\}$
   without enforcing L9a's units, and $q=59=w$ divides $z$. The theorem states
   the guard; the verifier omitted it. Fixed, and the instance count rose to
   2136.
10. **Method note.** Four advisories and one adversarial audit were upheld this
   session; two of them killed claims I had already written and frozen. The
   independent-verification habit that caught nothing here (recomputing symbols
   over supp$(d)$) was itself the blind spot: it re-derived the same wrong set.
   Checking the *hypotheses* of the lemma being invoked, not just recomputing
   its conclusion, is the lesson.
11. **Next.** For the 190: step (ii), the Schinzel condition, is all that
    remains. For the 163: the prime-$b$ family is dead, so a different $b$-shape
    is required - $b$ with two prime factors, or $b$ carrying a fixed square
    factor, changes the $q$-side symbol and may break the pairing that L11g
    exploits.

## Session 2026-08-16 - L12: the wall generalizes, and the door is named

1. **The obvious repair.** L11g's pairing is $(x,d)_A(x,d)_q=-1$ with
   $(x,d)_A=(d\mid A)=(2\varepsilon q\mid A)$. Put a fixed odd $m$ into $b$:
   $b=\varepsilon mq$. Then $(x,d)_A=-(m\mid A)(q\mid A)$, L9a still gives
   $(x,d)_q=(A\mid q)=(q\mid A)$, so the product is $-(m\mid A)$ - choose a
   non-residue $m$ and the wall flips to $+1$. And $b$ is an existentially
   quantified witness, so "prime $b$" was never part of the definition: the
   change costs nothing.
2. **Why it fails.** An advisory caught it before I ran the search. At $p\mid m$
   the roles of $m$ and $q$ swap and the *same* $u_x$ computation applies:
   $v_p(b)=1$, $P(b)\equiv P(0)$, $v_p(x)=-4$ even, $v_p(d)$ odd, and
   $u_x\equiv\delta^2Aa^4Z^4/D_z^2\equiv A$. So $(x,d)_p=(A\mid p)=(p\mid A)$,
   and $\prod_{p\mid m}(x,d)_p=(m\mid A)$. Total:
   $-(m\mid A)\cdot(m\mid A)=-1$. **The $-1$ relocates rather than cancels.**
3. **The generalization (L12a).** Nothing in that argument used the shape of
   $b$. For every admissible $b$ whose odd places are coprime to the cell data,
   $\prod_{\{A\}\cup\operatorname{oddsupp}(b)}(x,d)_p=-1$. Verified on 231
   in-suite instances spanning prime, semiprime, three-prime, squarefull and
   rational $b$, on all four branches, with $(x,d)_p=(A\mid p)$ asserted at
   every odd place. **The entire coprime-$b$ direction is closed** - L11g is
   the case $b_{\mathrm{sf}}=q$.
4. **The door (L12b).** Exactly one hypothesis remains: coprimality. With
   $b=\varepsilon fq$, $f\mid z$, the step $P(b)\equiv P(0)$ fails at $f$ and
   the telescoping breaks. It does: **103 of the 163 L11h-walled cells** get
   every controlled symbol $+1$ on the conservative set (every fixed datum
   frozen). Two of twelve spot-checked witnesses were additionally globally
   soluble; at the rest the surviving $-1$ sits on wild places, i.e. step (ii).
5. **What is NOT claimed.** Base-point alignment only. Class constancy follows
   from the exponent lemma but the conservative moduli are $\sim10^{13}$, so
   sampling the progression is out of reach and I do not claim it verified. I
   ran the sampling anyway and it returned 0 usable primes within budget -
   recorded here so the gap is visible rather than papered over. Shrinking the
   controlled set needs L11f's content argument redone at places dividing the
   data; that is open.
6. **Next.** Redo the content argument at $p\mid z$ to shrink the modulus, then
   sample the class; and push the escape search past 103/163 (wider $a$, $f$
   drawn from $D_z$ as well as $z$, larger $q$).


## Session 2026-08-17 - L13: eight agents, class side verified, 352/353 witnesses

Parallel session: eight agents (class audit, minimal controlled set, escape
widening, P-factorization, square-hunt, density probe, literature, adversarial
audit) plus four lead-run witness hunts. Every agent deliverable was
independently replayed by the lead before freezing; nothing rests on an
agent's word alone.

1. **Normalization fact (three agents independently).** The kernel's `_l10_P`
   evaluates to $A^2 b^4 D_z^2 M_0$ - its internal $\mathrm{Ng}$ carries one
   factor $A$ versus the prose formula. $A^2$ is a square: symbols, valuation
   parities, factor degrees all invariant. Asserted at every L13 replay point.
2. **L11h classes hold at depth (ClassAudit).** 162 members audited across 30
   sampled rows: frozen symbols $+1$ at every member (largest $Q=1.04\times
   10^{11}$), wild and $\infty$ aligned, parity law exact; 5 zero-bad members
   found (one deep, cell (5,1/7) at $k=65$); deciders agree everywhere;
   refusals logged, never evidence.
3. **Minimal controlled set (EscapeSet).** For the 103 escapes at $a=1$:
   $S_{\min}=\{2,3,5,7,f\}$ (at $a=1$: $\alpha=4$ square, $\delta=-4/5$).
   Moduli collapse from the conservative $10^{13}$ to median $2.4\times10^7$.
   `l12_class.py` freezes `l12_class_cert`; 883 prime members verified, zero
   violations; excluded-member set (174 primes) named per row. Emergent
   symbols genuinely hit $-1$ (467/959 spot-probes) - step (ii) is real and
   outside the class claim; count always even (39 fully factored members).
4. **The 60 walled-no-escape cells (EscapeWiden; verified by the lead).** All
   60 carry fully soluble witnesses $b=\varepsilon fq$ ($a\in\{3,5,33\}$):
   55 Hasse-Minkowski, 5 steered certificates. Includes residual cell
   $(31,-2)$ ($a=3$, $b=1457$).
5. **Residual-cell hunts (lead).** $(67,(2,1))$: $a=1$, $b=4757=67\cdot71$.
   Then at the L11h $a$-values: $(41,\tfrac17)$ $a=25$ $b=4985$;
   $(61,-\tfrac27)$ $a=25$ $b=55571$; $(61,(2,1))$ $a=25$ $b=-71431$
   (parallel six-process hunt). All tied-status True, independently replayed.
   **Grid witness coverage 352/353; the last cell is $(89,-2)$** (class cell;
   its members exceed the primality engine; off-class scans through $a=41$,
   $q\le1200$ have not hit).
6. **No shortcut (FactorP + SquareHunt).** $P$ is irreducible degree 8 over
   $\mathbb Q$ on all 706 (cell, branch) rows (sympy; in-suite Frobenius
   certificates mod $p<1000$ - the first extended run FAILED at (31,-3/1)
   'can' because an irreducibility prime need not exist below 100: 8-cycle
   density $\approx1/8$; range widened, comment frozen). Content sqf
   $\in\{\pm1,\pm5\}$. 137.5M square-class tests, 40 designed families: zero
   hits. Galois: $D_8\wr C_2$ evidence only (28,240 certified unramified
   samples, exactly its 10 cycle types, $S_8$ excluded with overwhelming
   Chebotarev evidence) - recorded as EVIDENCE, not proved.
7. **Density (DensityProbe).** 2262 members across 40 classes: alignment
   perfect everywhere; 9 zero-bad among 44 fully factored members
   (`ramified(x,d)` empty, tied True 9/9). Step (ii) has positive empirical
   density - this is the exact place the analytic input is needed.
8. **Audit (L12Audit).** L12a/L12b confirmed sound (confidence 0.95); two P3
   proof-text imprecisions fixed in THEOREMS (branch-specific $(d|A)$ step
   replaced by the branch-uniform $(u_d|A)$ derivation; $v_p(x)=-4$ corrected
   to $-4v_p(b)$ with the denominator-prime duality).
9. **What would close $(89,-2)$.** Wider $q$ at $a\in\{17,25\}$, the $\tau=1$
   branch, or a class member below the primality threshold (its modulus
   $N=1.65\times10^7$, base $q=59$: the $k=1$ member already factors
   $P$ at $\sim10^{60}$).

## Session 2026-08-18 - L13e, the frozen witnesses, and the last cell

1. **CompositeWall upgrade.** The L12a pairing never used the primality of A:
   at odd p | A the middle term -delta_tau a^4 Z^4 Ng^2 strictly dominates P
   (v_p(Ng) = 0 always), so [x]_p = [p]; every prime dividing 1+4a^2 is 1
   mod 4, reciprocity factor-by-factor, Jacobi cross terms cancel. The wall
   prod over p | A (odd valuation) and odd places of b = -1 now holds for
   EVERY admissible a; v_2(b) = 1 flips the sign - the exact admissibility
   boundary. L11h's necessity (reachable iff z has a 1-mod-4 numerator
   prime) is now a THEOREM for all admissible a; the 878,400 certificates
   merely confirm. Derivation + 316-instance sweep in
   data/l13_compositewall.{md,json}; kernel replay block (40/120 instances)
   in _verify_L13(e).
2. **Witness freezing.** The two previously unfrozen closures (29,(-2,1))
   (a=-15, b=-29/5, tau_d with d=1/3) and (53,(1/3)) (a=-15, b=9, square
   branch) are now frozen in _L13_RESIDUAL - both tied-status True,
   (53,(1/3)) also carries a steered certificate with wild prime
   253408554167897333140242269675351. The L13Audit finding "rests on an
   un-verifiable-in-suite input" is now closed: the whole 352/353 coverage
   is replayable in-suite.
3. **Zero-bad growth (ZeroBadExtend).** 11 new frozen zero-bad members (11
   new cells) on top of the 5: 24 certified members across 15 distinct
   cells; 15/57 fully factored members in its sweep - step-(ii) density is
   ~25% wherever factorization is feasible, and it never vanishes
   simultaneously at the class members that DO factor.
4. **L13Audit findings:** all six applied - banner/docstring 352/353, the
   f|z witness descriptive clause (3 of 60 rows carry q=w), emergent text
   arithmetic (959/927/452), (53,(1/3)) now frozen, and L13b/RESULTS
   arithmetic (64 directly witnessed cells = 60 escape2 + 4 residual).
5. **The last cell (89,-2).** Closer89 attack scan, all paced at 0.1s with
   single subprocesses: phase 1 (its own class members) refused; phase 2
   (tau=1, a in {1,3,5,9,17,25,33}, f in supp(num z) u supp(num D_z),
   e in {1,2}, eps, q < 3000) refusals dominated by M numerators ~65
   digits above _MR_LIMIT - of 5346 tested: 655 decided (all False), 4691
   refused. Remainder tracked per phase in data/l13_cell89.json. The cell's
   L11h class itself (a=17, N=16521960) has members too large for the
   primality engine from k=2.
6. **Discipline.** CPU-capped 50% at all times: one subprocess at a time,
   0.1s sleeps between heavy calls. Extended + default suites both green
   (44.2s / 5m).

## Session 2026-08-18 - L13f: the square/smooth zero-bad filter, and an honest scan of the last cell

1. **l13_filter.py development.** The emergent-free test at a member b
   does NOT need P(b) factored: strip the primes of the frozen set
   {2,3,5,7} u supp(alpha) u supp(delta) u supp(b) from P(b), numerator
   and denominator separately - the member is zero-bad iff both stripped
   remainders are perfect squares (`emergent_free`). Pure integer
   arithmetic (isqrt): NO factorization, NO engine refusals. The smooth
   variant `smooth_emergent` (y = 1e6) decides exactly when the stripped
   remainder is y-smooth (statuses `zero` / `align-fail` / `bad`) and
   returns `cofactor-big` otherwise - that subfamily is deferred to a
   bounded stage-2 pass (`cofactor_decide`), never guessed at.
2. **Validation A1.** All 16 frozen `_L13_ZERO_BAD` rows replay True
   under the perfect-square test.
3. **Validation A2.** Against the ClassAudit ground truth
   (`data/audit_l11h_members.jsonl`): agreement on every row where the
   test is decidable - `emergent_free` agrees with `bad_count == 0`
   everywhere; `smooth_emergent` agrees wherever the remainder is
   y-smooth. The deep zero-bad row (5,(1,7)) at k=65 has a y-non-smooth
   stripped remainder, so the smooth decider skips it BY DESIGN. Its
   48-digit remainder is PROVABLY prime (kernel _is_prime, Pocklington,
   0.3 s), so the parity law forces its single emergent symbol to +1 -
   the member is zero-bad by the prime-cofactor mechanism, which the
   stage-2 `cofactor_decide` ladder formalizes as a regression row.
4. **Phase B QS1 outcome (honest scope).** The scan decided 290,928
   structured candidates for cell (89,(-2,1)) - a in {1,3,5,9,17,25,33},
   branches tau in {sq,can,one}, b = +-f^e*q and +-f^e/q with f in the
   cell's per-a supp pools, e in {1,2}, q prime <= 3000 - with **0
   zero-bad hits**. The per-candidate hit prior is ~1e-11, so 0 hits is
   NOT evidence about the cell and this is not a "cell scanned" claim;
   only the exact counters are recorded. The scan says nothing about the
   non-smooth cofactor subfamily (stage-2 `cofactor_decide` pending).
   Cell (89,(-2,1)) remains the unique open cell; grid coverage 352/353
   unchanged. QS2 (a=17/sq, q in 3001..8000) crashed on a benign
   NameError in the hi-q loop (stale helper call); fix pending, no
   result claimed.

## Session 2026-08-18 (continued) - L13f stage 2: the cofactor ladder closes the last cell; GRID COMPLETE 353/353

1. **The stage-2 ladder.** `cofactor_decide` in l13_filter.py: a member is
   zero-bad whenever its stripped remainder is (i) a perfect square or
   (ii) a single PROVED prime with odd valuation - the parity law then
   forces the sole emergent symbol to +1. No full factorization. Kernel
   primality only (MR / Pocklington); composites <= 72 digits go to the
   budgeted kernel factorint rung; bigger cofactors are refusals with
   reason codes, never evidence. Ground truth: 16/16 frozen zero-bad rows
   ladder-zero (incl. the (5,(1,7)) k=65 prime-remainder row);
   INCONSISTENT = 0 across 5,767 decided rows.
2. **Box results (AttackA, data/l13_filter_run2.json).** 309,392
   structured candidates (QS2 fixed + QS1 re-run): 10,515 aligned,
   6,034 smooth-bad (exact emergent -1 pair), 4,481 cofactor-big stored
   -> 236 CERTIFIED SOLUBLE (133 proved-prime rung + 103 factorint rung),
   516 unproved-Jacobi (evidence tier), 1,896 cofactors >72 digits and
   1,637 factorint refusals (1,408 FactorBudget / 229 PrimalityBound).
   Dense box (DenseScan, stage D, 1,286 aligned |b|<20001 rows): 89 more
   certified soluble (46 + 43). Every hit cross-checked tied=True by
   _l7_tied_status.
3. **Dense-box structural finding.** In the dense box every aligned row
   has a=17 (the L11h class a-value): 2,294 sq + 747 can; all other
   (a, branch) cells are 100% align-fail - the alignment wall at this
   cell is total outside the class's a-value (the tau=0 branch and
   two-prime b-shapes, covered here for the first time, added zero
   symbol-pass rows).
4. **Lead replay and freeze.** The frozen witness (89,(-2,1)): a=3,
   b=89/367, tau=19/37, independently replayed rung by rung - smooth
   cofactor-big -> 29-digit remainder PROVED prime (Pocklington) ->
   v_R=1 -> hilbert(x0,d0,R)=+1 -> full emergent recount {43, 90947,
   204917, 471137} all +1 -> _l7_tied_status True. Frozen as row 7 of
   _L13_RESIDUAL; both suites exit 0 (default ~40s, extended ~99s).
   **Grid witness coverage 353/353 - COMPLETE.**
5. **Scope honesty.** The ladder closes THIS GRID (one witness per cell),
   not hypothesis H: H requires an emergent-free member of each ALIGNED
   CLASS uniformly, which remains the open analytic input to the
   conditional record - and the ~10^-11/candidate square-subfamily scan
   (item 4 of the previous section) was never evidence. Refusal counters
   above are recorded with reasons and are not claims.

## Session 2026-08-18 (continued) - the prime-prover upgrade and the per-class H scan

1. **Hypothesis H scan, wave 1 (`l13h_scan.py`).** Per-class zero-bad
   member certification across all 293 aligned classes (190 L11 + 103 ESC):
   members Q = q1 + kN decided by the L13f ladder, stop-on-first-zero.
   Result: **97/293 classes closed** (L11 64/190, ESC 33/103), every
   closure a ladder verdict 'zero' with tied-status cross-check True.
   LadderQA regression: PERFECT agreement (29/29 prior zero rows 'zero',
   16/16 prior bad rows bad, 0 alarms). Dominant wave-1 blocker:
   316+ zero-jacobi EVIDENCE rows - members whose stripped remainder is
   almost surely prime but whose primality the kernel's sqrt(n)-Pocklington
   could not prove at 60-118 digits.
2. **Kernel upgrade: BLS-relaxed Pocklington (h10q.py).** `_pocklington`
   now implements the Brillhart-Lehmer-Selfridge F >= n^{1/3} relaxation
   (Crandall-Pomerance Thm 4.1.5): with the same per-q gcd conditions a
   composite would have exactly two prime factors uF+1, vF+1, detected by
   the exact discriminant test on n-1 = F(c2 F + c1), D = c1^2 - 4 c2;
   also trial division widened to primes <= 1e6 with an early abort once
   F >= n^{1/3} (exact integer cube root). Unit check `test_bls.py`:
   3 pinned-regime proofs (old code raised PrimalityBound) + composite
   cousins; both suites exit 0 (default ~51 s, extended ~118 s).
   Engine discipline unchanged: every proof is exact (deterministic MR /
   Pocklington/BLS); jacobi-tier rows remain evidence, never claims.
3. **Wave 2 (running).** Re-scan of the 196 open classes with the
   upgraded prover; conversion accounting (jacobi -> proved) per class
   will be recorded in the data files and the ledger.

## Session 2026-08-18 (continued) - L14: hypothesis H verified on the entire grid (293/293)

1. **Wave 2 results.** BLS-upgraded prover: +20 newly closed (117/293
   total), every new closure a jacobi-> PROVED conversion; zero alarms.
2. **Wave 3: alternate classes.** H needs ONE class per cell. Driver
   `l13h_alt.py`: L11 Hensel-root a-lifts (w | 1+4a^2), ESC f/eps
   variants; 164 closures: AltL11W 53/56, AltL11E 48/55, AltESC 63/65.
   All closures self-validated (per-member alignment recompute;
   close_cell ramified/tied cross-checks; refusals logged, none claimed).
3. **Final 12.** AltHardL11 10/10 (w mod4=3 cells via the a=1 path;
   last L11 holdouts (73,(2,1)) a=23 q337 k83; (89,(-1,1)) a=17 q1097
   k=0; (89,(7,3)) a=17 q2347 k=0 factorint); AltHardESC 2/2
   ((59,(7,3)) f=7 eps=-1 q709 k51; (67,(-1,1)) f=67 eps=-1 q811 k=0);
   DeepK12 banked 7 more before parking. Incidents resolved in-band
   (double-writer on AltHardL11's file, clean single-writer rewrite).
4. **Lead integration.** Merged 293 unique cell closures across 16
   wave files into data/l13h_all_closures.json; lead-wide replay
   (`l14_replay_all.py`, brent-budget raised for the ramified recheck):
   **293/293 ladder-zero verified; ramified empty 271, refused 22**
   (cross-check limitation, logged). Output data/l13h_replay.jsonl.
5. **Freeze.** `h10q.py::_L13_H_CLASSES` (293 rows) + `_verify_L14`:
   default 60-seed exit 0 (~90 s total), extended all-293 exit 0
   (~230 s, in-suite replay 140 s). RESULTS/THEOREMS/CONDITIONAL updated;
   PROGRESS entry added. Waves' QA: LadderQA 45/45 exact agreement
   (29 zero + 16 bad), replays by every wave agent, zero
   CLASS-CONTRADICTION/INCONSISTENT across all runs.
6. **Data files.** data/l13h_{l11_a,l11_b,esc}{,_2}.jsonl (waves 1-2),
   data/l13h_alt_{l11_a,l11_b,esc,l11_hard,esc_hard}.jsonl + deepk_l11
   (wave 3), data/l13h_qa.{jsonl,md}, data/l13h_all_closures.json,
   data/l13h_replay.jsonl, test_bls.py (BLS unit check).
7. **Frontier after L14.** H on the grid is CLOSED instance-wise; the
   uniform all-w statement is the remaining analytic obstruction
   (literature: prescribed square class unconditional to degree 2
   (Krumm), degree 3 conditional on elliptic Parity, degree 8 open).
   Member-density histograms across waves are EVIDENCE for it.

## Session 2026-08-18 (continued) - L15: the remainder and counting laws; the H conjecture goes quantitative

1. **Exact remainder audit (PrimeRungAudit agent).** Independent
   recomputation of all 293 closure members (`data/l15_remainders.jsonl`,
   294 lines): every record re-verified end-to-end through the ladder;
   exact factorization of the stripped remainder under discipline; **0
   refusals, 0 alarms**. Law: $R$ squarefree everywhere; $R=1$ never
   occurs; 197 prime rungs all $R=p$ exactly (median 48-digit $p$, max
   119); 96 factorint rungs with $e\in\{2,3,4\}$ emergent primes
   (76/15/5); all 15 odd-$|E|$ rows land exactly on the parity law
   (every symbol $+1$).
2. **Counting law (DensityModel agent).** Deterministic 2.8-s playback
   (`data/l15_density.json`): $k_{\text{zero}}=0$ in 207/293 (70.65%;
   L11 65.3% / ESC 80.6%); median 0; tail (86 classes) median 28, max
   1694; emergent-free rate per prime member $p_0\approx0.040$ (L11
   0.035, ESC 0.072) ~15.7× slower than pure Bateman–Horn on the same
   AP; resistant subclass $(a=1, w\equiv3\bmod4)$: $k{=}0$ only 30.8%;
   both outliers (1694, 1511) are in it (L11, $N=10080$, cells (79,(5,1)),
   (83,(-5,1))). Small-member phenomenon: observed $P(k{=}0)$ boosts
   12× over the $1/\log N$ baseline.
3. **Quants for the frontier.** The uniform all-$w$ form of H now has
   testable shapes: (a) "$p_0>0$ uniformly away from the $(a{=}1,
   w\equiv3\bmod4)$ ridge, with geometric tail"; (b) the emergent-
   product symbol $\chi(k)$ as a function of $k$ — L16 test runs
   (data/l16_char.jsonl): if $\chi$ is a periodic character of $k$,
   H reduces to a Schinzel-primes statement in a co-sequence.
4. Ledger rows: THEOREMS L15, RESULTS L15, CONDITIONAL evidence table
   updated; README inventory extended (l15_/l16_ artifacts).

## Session 2026-08-18 (continued) - L16: emergent-symbol structure; the sieve-shaped small layer

5. **L15Char failure + design correction.** The first L16 agent died in
   setup; re-reading `smooth_emergent` before the rerun exposed a design
   flaw in the brief: the ladder's `em` list holds only the −1 emergent
   primes, and by the parity law (frozen symbols +1 by construction) the
   product $\chi(k)$ over ALL emergent primes is identically $+1$ on
   aligned members — the original $\chi$-probe would have measured a
   trivial constant. Lead ran the corrected probe itself
   (`l16_emergent.py`, 8.4 s): full sign-decorated emergent lists
   per prime-member $k\le600$ on 7 sampled classes →
   `data/l16_char.jsonl`.
6. **Results.** (a) every recurring small emergent prime's symbol is
   residue-determined ($k\bmod p$) in all 7 classes, zero
   counterexamples — analytically expectable ($b$ affine in $k$) and now
   measured; (b) clean-small-layer rates 35–58% vs closure $\sim4\%$ ⇒
   the big-cofactor layer carries its own $\sim20\times$ obstruction
   share; (c) parity bookkeeping confirmed on every row (odd small-neg
   ⇔ odd cofactor-neg); (d) probe cross-check reproduces the closure
   shape at every sampled frozen $k_0$ (small layer clean, cofactor the
   documented big prime/factorint).
7. **Interpretation for the frontier.** After the strip, the small
   obstruction is a sieve condition on $k$ (union of forbidden residue
   classes per emergent prime); the remaining analytic content of H
   concentrates in the cofactor layer, which parity forces to be sparse
   (one prime, or a small product of primes, on 33% / 67% of closures).

## Session 2026-08-19 - L17: mechanistic p-adic sieve (matched fit), cofactor proved = reciprocity-level, band-only rates; advisory-driven corrections absorbed

1. **Three agents (FullGridSieve, RateModel, CofactorSign), artifacts
   `data/l17_{sieve,ratemodel_censored,cofactor}.jsonl`.** Residue-
   determination extended to **293/293 classes** (8,760 prime members,
   $k\le119$, zero counterexamples); exact odd-valuation bad-mass
   tables ($1/(p{+}1)$ per simple root; step roots recursively lifted
   with certified residual bounds; 56 partial classes, 67 excluded
   pairs labeled); cofactor proved statements are reciprocity-level
   (197/197 auto-$+1$ singletons given frozen/small $+1$; parity
   product $+1$ on 293/293).
2. **Advisory corrections absorbed** (all correct; ledger updated):
   (a) recurring-only model overprediction is mechanical omission, not
   cross-prime correlation; (b) matched-support test conditioned on
   the identical member sample; (c) splits by final rung = outcome
   leakage — RateModel rebased to the predeclared wave-1 cohort with
   right-censoring; (d) unknown (refusal/zero-jacobi) trials counted
   both ways, **two-sided band only** (c ∈ [0.017, 0.80] per family);
   (e) raw $\prod f_p/p$ factors → exact $p$-adic masses
   ($1/(p{+}1)$ per simple root) with recursive lifting for singular
   roots; (f) malformed script hard-reset with py_compile + smoke
   gates; (g) "no bias / indistinguishable from noise" was 3.5-SE
   wrong at the mean — the $+0.019$ offset is GENUINE; the Q-
   conditioning explanation I proposed had the sign *backwards* (the
   Q-root residue is safe: $2\alpha b\equiv0$ there ⇒ `bad_sign=False`);
   resolution: the offset is the **Haar($\mathbb Z_p$) vs $k$-window
   $[0,119]$ gap** of the class factors (characterized on the matched
   sample, not a dependence); (h) the older pairwise co-occurrence
   dependence tests are formally **underpowered** — and the powered
   marginal-preserving permutation test now gives $z=-0.09$,
   empirical tails $0.61/0.53$ (no dependence); (i) cofactor "complete picture"
   retracted: closure rows are selected for all-$+1$ symbols, so
   individual factorint signs + residue tables are closure-conditioned
   description only — a nonclosure cofactor sample is required for any
   structure claim beyond the proved singleton/product reciprocity.
3. **Ledger positions**: sieve structure at grid scale PROVED;
   matched model validation EVIDENCE (marginal ratio 0.9993; union
   residual −2.9/8760; permutation z=−0.09, P=1.0); cofactor proved =
   reciprocity-level; rate predictions band-only. THEOREMS L17,
   RESULTS L17; matched statistics persisted in
   `data/l17_matched.jsonl`.
4. **Methodology saga (matched test).** Successive instrument errors
   and theoretical candidates were each caught by review, then resolved
   by the matched characterization: (a) recurring-prime-only support
   (mechanical omission); (b) wrong-sign prime-Q conditioning theory
   (the Q-root residue is safe, not bad: $2\alpha b\equiv0$);
   (c) cohort-vs-closure PROGRESSION mixing in the lead's matcher —
   sieve rows are rebuilt from the closure-authority classes while the
   first bad-roots table indexed the wave-1 cohort (e.g. ESC
   (3,(-3,1)): cohort q1=41/eps=+1 vs closure q1=59/eps=−1). Final
   matched verdict: roots reconciled 5,375/5,375; marginal ratio
   1.0000; per-class hit-conditioned residuals mean −0.00016/SD 0.012;
   Final verdict (artifact v2, resolved-mass prefixes restored):
   marginal ratio 0.9993; hit-conditioned union residual −2.9/8760;
   permutation null z=−0.09, empirical tails 0.61/0.53 — the
   original-estimand ≈+0.0195/3.6σ
   is the Haar($\mathbb Z_p$) vs $k$-window $[0,119]$ gap of the class
   factors (characterized, not a model failure).
5. **CofactorNonclosure: power-gated INCONCLUSIVE.** 400 nonclosure
   rows (60 classes) sampled from the 2,983 deferred-refusal pool: 11
   exact factorizations / 389 refusals (FactorBudget 309,
   PrimalityBound 77, digit-cap 3). Parity structurally consistent
   (0 odd-$n_{neg}$ violations over the 11); even-parity histogram
   (0:2, 2:8, 4:1) descriptive, null $p=0.784$ descriptive — the
   predeclared power gate ($\ge150$ resolved) FAILS, so the statement
   is INCONCLUSIVE: cofactor-sign distributions beyond the proved
   reciprocity-level statements are measurement-limited by the prover
   (`data/l17_cofactor_nonclosure.jsonl`, global-summary carries the
   gate + histogram). Labels the question open, not answered.

## Session 2026-08-19 (continued) - L17-horizon: the off-grid machinery; Schinzel audit; class-risk factors; reproducibility hardening

1. **Off-grid machinery exists and works.** Horizon agents discovered the
   fresh-closure recipe: L11 roots $1+4r^2\equiv0\ (w)$ give $a$-pools and
   canonical $\tau$; `_l10_class_cert` certifies alignment; an off-grid
   L12 escape requires temporarily registering
   `_L12_ESCAPE[(w,(-1,1))]=(w,\varepsilon,q_1)$ before
   `l12_class.l12_class_cert`. Four cells closed and independently
   replayed by the lead: `[101,[-1,1]]` (L11 $a{=}5$, $k{=}16$, factorint
   rung), `[103,[-1,1]]` (ESC $f{=}103$, $k{=}0$), `[107,[-1,1]]` (ESC
   $f{=}107$, $k{=}0$), `[109,[-1,1]]` (L11 $a{=}71$, $k{=}0$; auxiliary
   tied/ramified budget-refused — the accepted L14 closure form). H109's
   optional 1,275-probe cross-$k$ scan ($k\le50$, 25 aligned classes)
   found no tied=True+ramified-empty pair: auxiliary, not a defect
   (`data/l17_horizon109_crossk.jsonl`).
2. **Schinzel audit** (`data/l17_schinzel_audit.jsonl`, generator
   `l17_schinzel_audit.py`): all 293 canonical closure pairs
   $\{q_1(t),F_{\text{cell}}(t)\}$ pass local conditions — $F$ irreducible
   (certified), value-gcd 1, pair-product gcd 1; counts done with
   $q_1(t)\not\equiv0\ (p)$ (advisory: NOT $t\not\equiv0$; for $p\mid N$
   the class is a nonzero constant so the count is $p$). 6
   `_L13_RESIDUAL` rational-$b$ witnesses segregated as pointwise
   certificates. Open side remains degree-8 prime values (Schinzel/H).
3. **Class factors** (`data/l17_classfactors.jsonl`, generator
   `l17_classrisk.py`): per-class $p\le10^4$ window factor intervals from
   the v2 resolved-mass roots; closure-effort cross-axis: Spearman
   $\rho=-0.178$, medians $0.5742$ (cls$=0$, $n=82$) vs $0.5662$
   (cls$\ge1$, $n=108$) — direction consistent, descriptive only.
4. **Resolved-mass patch (the decisive correction).** RateModel's v2:
   `_singular_bad_mass` now persists `resolved_mass` per excluded root
   (67 pairs; max unresolved residual $6.9\times10^{-15}$), so the
   p≤10000 window interval $[\,\prod(1-k-r),\;\prod(1-k)\,]$ is valid
   full-grid. Matched stats: marginal 0.9993 (4481.3/4478), union
   residual $-2.9$/8760, permutation $z=-0.09$ tails 0.61/0.53, estimand
   $\approx+0.0195$/3.63σ. Review process caught and fixed: artificial
   $t$-counting in the Schinzel test; discarded resolved masses; the
   vacuous integer-support p-value; stale strings (five-file sweep).
5. **Reproducibility hardened (scoped)**. Repo-regenerable chain:
   `l17_ratemodel_padic.py` → `l17_badroots_closures.py` → `l17_matched.py`;
   `l17_classrisk.py`; `l17_schinzel_audit.py`;
   `l17_horizon{101,103,107,109}.py`; `l17_horizon109_crossk.py`
   (byte-identical regeneration verified); probes `l16_emergent.py`,
   `l14_replay_all.py`. Explicitly OUT OF SCOPE (persisted evidence, no
   checked-in producer): `l17_sieve.jsonl`, cohort `l17_badroots.jsonl`,
   `l17_ratemodel_censored.jsonl`, `l17_cofactor.jsonl`,
   `l17_cofactor_nonclosure.jsonl` — the older session-side scripts were
   not recovered into the repo (one is corrupted in /tmp); README carries
   the reproducibility-scope section saying exactly this.

## Session 2026-08-19 (continued) - L18: the step-(ii) density law, horizon wave 2, and the deliverable bundle

1. **The uniform-bound hope is dead, and that is progress.** StepIIDensity
   proved $m_3=m_5=m_7=0$ and $m_p\le8/(p+1)$ for $p\ge11$ (singular roots
   with multiplicity), then closed the door: every one of the 293 class
   polynomials has a replayed irreducibility certificate and a replayed
   nonsquare witness at a simple bad root, so the Chebotarev exponent is
   $c=1/2$ and $\prod_{p\le X}(1-m_p)\asymp(\log X)^{-1/2}$. The lead
   re-verified the local bound independently on all 102,439 stored pairs
   with $p\le10^4$: zero violations. Measured exponents 0.4778 /
   0.4801 / 0.4590 agree with the theoretical $1/2$.
2. **Reframing H.** Prime members up to $X$ grow like $X/(\log X)^2$ and
   each is clean with probability $\asymp(\log X)^{-1/2}$, so the
   expected count of emergent-free members diverges at rate
   $X/(\log X)^{5/2}$. H is a divergence statement; the shortcut through
   a uniform positive constant is provably unavailable.
3. **Horizon wave 2**: $w=113$ (ESC $f{=}113$, $q_1{=}41$, $k{=}0$;
   auxiliary budget refusals, accepted L14 form), $w=127$ (ESC
   $f{=}127$, $q_1{=}149$, $k{=}0$, tied True/ramified empty), $w=137$
   (L11 $a{=}87$, $\tau=15139/30277$, $q_1{=}149$, $k{=}0$, tied
   True/ramified empty) - all lead-replayed. Off-grid set is now
   $\{101,103,107,109,113,127,137\}$.
4. **A new wall shape at $w=131$** (PROVED within the protocol): L11 is
   empty because $(-1|131)=-1$; the canonical escape route forces
   $\mathrm{Hilb}_5=-(5|q_1)$, so the frozen place 5 and the wild place
   $q_1$ cannot both be $+1$. 336 attempts, 326 eligible obstructed, 0
   aligned, 0 refusals, member ladder never entered. Sharp open
   question: which primes $w$ have this 5-symbol collision?
5. **Composite $w$**: out of scope by construction ($w$ indexes a place;
   the break is at the residue-field constructor). Downstream arithmetic
   generalizes after prime decomposition - pseudo-cell $[21,[-1,1]]$
   closes at the honest prime $f=3$. Recorded as scope clarification.
6. **Deliverables**: two paper drafts staged under
   `deliverables/2026-08-19/papers/` - the main conditional-$\forall_6$
   paper (17 pp, sections authored per-layer and lead-integrated) and
   the companion verification/reproducibility paper (6 pp). Bundle built
   by `deliverables/make_bundle.py` with MANIFEST + SHA256SUMS; the ZIP
   is written outside the staged tree.

## Session 2026-08-19 (continued) - L19: the Schinzel upgrade, the canonical 5-wall, and its break

1. **The conditional record's analytic input is now a named classical
   conjecture.** SchinzelImpliesH proved: for a fixed verified aligned
   class, writing $F(t)=P(\varepsilon f(q_1+Nt))=c\,G(t)$ with $c$
   $S$-supported and $G$ primitive of positive lead, Schinzel's
   Hypothesis H for $\{q_1+Nt,G(t)\}$ gives infinitely many members with
   ladder verdict zero. The mechanism is clean: the $S$-strip leaves the
   single odd place $R=G(t)$, the moving prime has $v_Q(P(b))=0$ because
   $P(b)\equiv-\delta a^4Z^4A^2\bmod Q$, all other outside places are
   even with symbol $+1$, the class certificate pins frozen/moving/real
   symbols to $+1$ for every $t$, and Hilbert reciprocity forces the
   last symbol $+1$. Hypotheses verified exactly on 293/293 canonical
   pairs; 24-row prime-rung cross-check; lead-replayed.
2. **A terminology bug was found and fixed.** CONDITIONAL.md defined the
   emergent set as "every outside valuation even", which is strictly
   stronger than the kernel contract (a place is appended only when the
   valuation is odd AND the symbol is $-1$) and is contradicted by the
   197 accepted prime-rung closures. Corrected in CONDITIONAL.md with an
   explicit note; the operational contract is what the machine verifies
   and what Schinzel implies. (The stronger wording would NOT follow
   from prime values - they deliberately create one odd place.)
3. **The canonical 5-wall (PROVED).** With $e=v_5(D)$: $v_5(x_0)=-3-2e$
   odd, $v_5(d_0)=0$, so $\mathrm{Hilb}_5=-(w|5)(q_1|5)$ while the wild
   place gives $(5|q_1)$; both $+1$ iff $(w|5)=-1$. With L11 available
   iff $w\equiv1\bmod4$, the canonical protocol is obstructed exactly at
   $w\equiv11,19\bmod20$. Audited on all 53 primes in $[101,400]$:
   26 L11 / 13 ESC / 14 obstructed, zero mismatches.
4. **And the wall is not intrinsic.** Route131 closed $[131,[-1,1]]$
   with $a=7$ ($A=197$), $\tau=0$, $f=131$, $q_1=41$, $k=0$ - the L12b
   non-coprimality door, consistent with the L10 $\tau=0$ coprime wall.
   HorizonSweepA then closed both remaining walls ($139$, $151$) with
   the same shape. Off-grid closed set:
   $\{101,103,107,109,113,127,131,137,139,149,151,157,163,167,173\}$,
   plus $[131,[5,1]]$; only $w=179$ still open (bounded probe).
5. **Divergence model (EVIDENCE).** Predeclared 239/54 split,
   $C=0.9011$; the L18 tail explains an $8.5\times$ delay against the
   ledger's $15.67\times$ Bateman-Horn gap - 51% of the excess
   multiplier, 78% of the log gap - with the residual $1.84\times$ and a
   clear failure at the selected $k=0$ spike (70.65% observed vs 16.67%
   predicted). Honest partial success.

### Citation audit and wall-scope repair (2026-08-19, advisory-driven)

- All 13 bibliography entries rebuilt with complete author/title/journal/
  volume/year/pages plus DOI and arXiv ids; parenthetical ledger prose
  removed; inline bracket references converted to `\cite`. Zero undefined
  citations or references in either paper.
- **Real error caught:** the ten-quantifier record (arXiv:2301.02107)
  carried the title of a *different* Daans paper. Correct:
  *Universally defining Z in Q with 10 quantifiers*, JLMS (2) **109**
  (2024), e12864. The other paper, *Universally defining finitely
  generated subrings of global fields*, is Doc. Math. **26** (2021),
  1851-1869 (arXiv:1812.04372). Also verified: Crelle 495 (1998) 1-28 is
  CT-Skorobogatov-Swinnerton-Dyer, *Rational points and zero-cycles on
  fibred varieties: Schinzel's hypothesis and Salberger's device*;
  Cornelissen-Zahidi is Contemp. Math. 270 (2000), 253-260.
- Schinzel-Sierpinski (Acta Arith. 4 (1958)) and Bateman-Horn (Math.
  Comp. 16 (1962)) added and cited, since the record now names them.
- **Wall scope stated exactly** (THEOREMS L19, paper Prop.): for the
  tau=0/a=7 closure, L10b's failing hypothesis is verbatim
  "b = eps*q1 with q1 != A prime" - b = 131*41 is not +/- a prime -
  while a odd, A=197 prime and A doesn't divide 2zD_z all hold; L12's
  failing hypothesis is coprimality of b to the cell data. Row verdict:
  kernel authority (`_l7_tied_status`, `ramified`, identity against
  `_l10_P`/`_sun_h`). Class certificate: labelled local generalization
  of `l12_class.l12_class_cert` (a != 1 and composite f allowed, every
  prime of f frozen), not the byte-identical kernel function. The
  novelty is the tau=0 alignment alone; the b-shape is the orthodox
  escape shape used in 103 grid rows.

## Session 2026-08-19 (continued) - L20/L21: class existence proved, irreducibility narrowed

The arc began by asking whether L19's per-class Schinzel implication could
actually be fed a class for every cell. L20 removed that structural
uncertainty. For any cell, take the fixed factor $f=w$, choose odd $a$ with
$A=1+4a^2$ and $(A\mid w)=-1$ using
$\sum_{r\bmod w}(1+4r^2\mid w)=-1$, solve the explicit residue system for
$q_1$ modulo $M=4A\prod_{p\in S}p$, and invoke Dirichlet. The system is
nonempty uniformly, so one verified aligned class exists for **EVERY**
cell and clause (ii) is gone. The replay made the finite shadow exact:
103/103 canonical rows reproduced after 1,113,000 residue visits, and
353/353 grid certificates were clean. A second derivation supplied a useful
adversarial check: if $a$ is frozen at $1$, collisions occur exactly at
$w=11,19,31,59,71,79$. That is precisely the 5-wall set obtained earlier
from the independent Hilbert-symbol route. Two derivations found the same
obstruction, and varying $a$ repaired every collision.

At the **historical L20 checkpoint, SUPERSEDED by L22**, the attempt to
promote the selected class to a fully admissible Schinzel pair returned
the split verdict: (b), (c), (d), and positive lead PROVED;
**irreducibility then OPEN — a status now superseded by L22**
(`data/l20_admissible.jsonl`). One premise also
broke during normalization.
Literal $S$-support of the content $c$ is **FALSE**: for
$(w,z)=(3,3/11)$, $a=1$, and $q_1=19$, its outside part is $11^{-12}$.
What the parity argument actually uses, and what the corrected proof gives,
is that the *square class* of $c$ is supported on $S$; the outside
valuation is even and therefore harmless. This changed the statement
rather than suppressing the counterexample.

At that same **historical, superseded checkpoint**, the adversarial chain
audit correctly returned NO: Schinzel H alone did not yet imply the
six-quantifier definition because the vertical input had been overclaimed.
L22 later proves that input (`agent://ElimAudit`;
`THEOREMS.md`, L22). The audit also caught six record defects.
The displayed $\Theta^\ast$ contained only the $\tau=0$ and $2a/A$ branches
although the construction used
$\tau^\dagger=(1+2a^2)/A$; the repair is to add that finite branch at zero
variable cost by L11a. A bounded set of 706 grid irreducibility checks had
been described as “each cell”, alongside the obsolete demand that every
outside valuation be even; the corrected contract keeps the finite scope
and allows odd valuations whose Hilbert symbol is $+1$. “A member of each
aligned class” was stronger than needed and was replaced by one selected
class and member per cell. The L13f square-or-single-prime rung had been
written as an if-and-only-if, but the 96 L14 factorint closures and 103
cell-box factorization-rung candidates show that it is only a sufficient
factorization-free shortcut. The ledger simultaneously called $w=179$
closed and open; the exact replay closes it, so every tried off-grid cell
is closed. Finally, “$G(t)$ is an $S$-unit” had the terminology backwards:
the proved condition is that $G(t)$ is a $p$-adic unit for every $p\in S$,
that is, coprime to $S$. Those six corrections fixed the scope of the
chain without strengthening any hypothesis by prose.

L21 then attacked the remaining word “irreducible” directly, and the first
result was a counterexample to the blanket claim. **L21a (PROVED)** says
that when $s=0$ (so $a=1$) and
$\delta_\tau=\sigma^2$ is a rational square,
$$
P(b)=\bigl(4DAb^2-\sigma a^2Z^2N_g\bigr)
     \bigl(4DAb^2+\sigma a^2Z^2N_g\bigr),
$$
a genuine $4\times4$ factorization. The identity held on 56/56 applicable
rows and all 8/8 controls were correctly inapplicable: blanket
irreducibility is **FALSE**. **L21b (PROVED)** explains why this does not
kill the construction. On
$\tau^\dagger=(1+2a^2)/A$ one has
$\delta_{\tau^\dagger}=-4a^4/A<0$, never a rational square, for every
$a$ and cell; 9/9 exact checks agreed. Its square class is that of $-A$,
so the quadratic field governing this branch is the imaginary
$\mathbb Q(\sqrt{-A})$.

The coefficient symmetry turned out to be real but almost, rather than
literally, reciprocal. **L21c (PROVED)** gives that
$P+32A^3s^2D^2b^5$ is palindromic and that $P_0=P_8$ always, with
288/288 exact checks; on the reciprocal slice the Galois group embeds in
$C_2\wr S_4$. L21e sharpened the boundary to an identity,
$$
P(b)-b^8P(1/b)=32A^3s^2D^2(b^3-b^5),
$$
so $P$ itself is reciprocal if and only if $a=1$.

The successful irreducibility proof was therefore global first and
fiberwise second. **L21d (PROVED)** localizes at the common endpoint
coefficient $L=a^8A^2Z^4$: a factorization over
$\mathbb Q(a,Z)$ would survive the specialization $(a,Z)=(1,27)$, while
the resulting octic has an exact irreducibility certificate modulo $17$.
Thus $P$ is irreducible of degree 8 over $\mathbb Q(a,Z)$. The proof then
respected the vertical-line objection instead of pretending that a
two-variable thin set cannot contain a whole fixed-$z$ line. Exact
one-variable certificates prove $P(a,z)$ irreducible in
$\mathbb Q(a)[b]$ for 353/353 target $z$; Cohen–Serre's
$O_z(\sqrt B\log B)$ bound for bad integers loses to the admissible
progression's $B/M+O(1)$ members, so a good $a$ exists for every certified
cell. The first admissible $a$ worked on 48/48 rows. In the wider hunt,
1024/1024 constructed-branch rows were irreducible and none was reducible;
the only reducible cases anywhere were the 128 off-branch,
$a=1$, square-$\delta$ rows already explained by L21a.

The cheap local proofs were also closed honestly. **L21e (PROVED)** gives
$P\equiv16A^2b^4(1-2As^2b)\pmod w$, so reduction modulo $w$ never
certifies octic irreducibility, and the $w$-Newton polygon never has one
denominator-8 slope. At an odd $p\mid A$ with $v_p(Z)=0$ there are exactly
two length-4 slopes $\mp v_p(A)/2$. The positive residue of this analysis
is the $a=1$ trace criterion: when $(5\mid w)=-1$, the degree-4 trace
polynomial is uniformly irreducible over $\mathbb Q(\sqrt{-5})$, and the
octic is irreducible whenever
$$
\Xi(Z)=(125D^2+64Z^4)(125D^2+1024Z^4)
$$
is nonsquare, as it was on 189/189 rows.

**SUPERSEDED by L22.** This paragraph records the exact L21 endpoint:
for a cell outside the $353$ certified values, the **historical,
superseded** chain still required a fibre certificate in addition to
Schinzel H
(`data/l21_irred_generic.jsonl`).  L22 now proves the fixed-fibre theorem
for every nonzero rational $Z$ on
$\tau^\dagger=(1+2a^2)/A$, $\delta=-4a^4/A$; the current chain assumes
classical Schinzel H alone (`/tmp/l22_elimination.md`;
`THEOREMS.md`, L22d).
 
## Session 2026-08-22 — L22: two fixed-branch closures and the Schinzel-only conditional record

1. **The local analyses found the shape of the problem without closing
   it.**  Higher-order Newton analysis at infinity split the repeated
   residuals into four irreducible quadratic clusters and left only
   factor degrees $2,4,6$ locally possible
   (`/tmp/l22_infinity.md`; `data/l22_infinity.jsonl`).  The geometric
   route constructed an exact Noether-matrix reduction but left
   $N((1-z^3)/z^6)\ne0$ OPEN
   (`/tmp/l22_fiber_geometry.md`;
   `data/l22_fiber_geometry.jsonl`).  The factor-tuple route proved a
   useful criterion and a no-go: reciprocity determines the product of
   outside signs, not their individual values
   (`/tmp/l22_factor_tuple.md`;
   `data/l22_factor_tuple.jsonl`).

2. **An apparent free-parameter closure was rejected.**  The
   $\lambda$-pencil proof and generalized aligned-class theorem are
   algebraically PROVED, but HIT chooses $\lambda$ after the cell and
   ranges through the infinite L11c family.  That branch parameter
   would be an additional witness, so the theorem is **NOT APPLICABLE**
   to the six-count (`/tmp/l22_square_branch.md`;
   `data/l22_square_branch.jsonl`).  BranchAudit identified the
   overclaim with confidence $0.91$; the report and summary were patched
   to preserve the theorem while rejecting it as the closure
   (`agent://BranchAudit`).  SquareBranch is therefore recorded only as
   an off-chain side theorem.

3. **Primary closure: fixed-$Z$ endpoint elimination.**  On the fixed
   canonical branch
   $\tau^\dagger=(1+2a^2)/A$, $\delta=-4a^4/A$, set
   $H=(A/4)P$.  The exact $A$-adic slopes make every global factor degree
   even; the $a=0$ and $a=\infty$ residual allocations eliminate degree
   $2$; and degree $4$ forces
   $Z^4=1024(1-Z)^2$, whose two quadratic branches have nonsquare
   discriminants $1152$ and $896$
   (`/tmp/l22_elimination.md`, §§1–4;
   `data/l22_elimination.jsonl`).  The excluded fibre $Z=0$ is explicitly
   reducible, while $Z=1$ has a primitive $a=1$ specialization
   irreducible modulo $11$
   (`/tmp/l22_elimination.md`, §5).  Thus for every fixed
   $Z\in\mathbb Q^\times$, $P(a,Z,b)$ is irreducible in
   $\mathbb Q(a)[b]$.

4. **Independent closure: reciprocal trace and elliptic descent.**  On
   the explicitly named specialization
   $a=1,A=5,s=0,\tau=3/5=\tau^\dagger,\delta=-4/5$,
   $P=b^4T(b+b^{-1})$.  Trace reducibility forces
   $125D^2+64Z^4$ square, but its reduced numerator is $5\bmod8$.
   For the reciprocal lift, a square
   $\Xi=(125D^2+64Z^4)(125D^2+1024Z^4)$ maps to
   $E:y^2=x(x+5)(x+80)$ with positive $x$; complete
   $2$-isogeny descent gives rank $0$ and only the four
   $2$-torsion points, none with positive $x$
   (`/tmp/l22_reciprocal_cube.md`, §§1–3;
   `data/l22_reciprocal_cube.jsonl`).  Hence the same branch
   specialization is irreducible for every rational $Z\ne0$.  The
   $a=1$ point is only a specialization witness, not a claim that
   $(5\mid w)=-1$ for every cell.

5. **Audits.**  ElimAudit independently rebuilt the endpoints, all three
   local polygons, the endpoint enumeration, both degree eliminations,
   and the exceptional fibres; verdict SOUND, confidence $0.90$, with no
   P0/P1 finding (`agent://ElimAudit`).  ReciprocalAudit independently
   re-derived the trace identity, mod-$8$ obstruction, norm direction,
   complete isogeny descent and torsion bound; verdict SOUND, confidence
   $0.92$ (`agent://ReciprocalAudit`).  Its one scope finding was
   load-bearing: every statement must name
   $\tau=3/5=\tau^\dagger,\delta=-4/5$, because sibling square-$\delta$
   branches are reducible.  The producer was patched, and the current
   report and JSON summary carry the branch explicitly
   (`/tmp/l22_reciprocal_cube.md`;
   `data/l22_reciprocal_cube.jsonl`).

6. **Chain consequence.**  For a fixed cell, L22 plus quantitative HIT
   selects an irreducible $a$ in L20's odd character-admissible
   progression; L20 supplies the aligned $f=w$ class and conditions
   (b)–(d), while L22 completes condition (a).  L19 plus classical
   Schinzel H supplies a member.  Therefore
   $$\text{classical Schinzel H}\Rightarrow H\Rightarrow\text{Theorem C}$$
   with the same finite branch menu and the same six universal
   quantifiers (`THEOREMS.md`, L19–L22;
   `CONDITIONAL.md`, §§1–2).  This is a PROVED implication with a
   CONDITIONAL conclusion: Schinzel H remains unproved, and
   H10/$\mathbb Q$ remains open.

## Session 2026-08-22 — L23/L24: the half-sieve frontier and the dyadic end of the diagonal

1. **The function-field Capell idea.**  On the canonical branch
   $\tau^\dagger=(1+2a^2)/A$ one has
   $\alpha=4a^4=(2a^2)^2$, so the tied quaternion is exactly
   $(P(b),2b)$ and the member surface is
   $$U^2-2bV^2=P(b)W^2.$$
   If $\theta$ is the root of the irreducible octic, making
   $2\theta$ a square in $K=\mathbb Q(\theta)$ would split the
   degree-$8$ discriminant point.  Capell makes the hope precise:
   $$2\theta\in K^2\iff P(u^2/2)\ \text{is reducible},\qquad
     N_{K/\mathbb Q}(2\theta)=256.$$
   The square norm is only necessary.  At the target prime, odd
   $v_w(z)$ gives an odd left-edge valuation for every $a$; after the
   free refinement $a\not\equiv1\pmod w$, the horizontal root
   $\beta=(2As^2)^{-1}$ has
   $2\beta=1/(As^2)$ in the nonsquare class $1/A$.  Thus the selected
   protocol cannot realize the Capell collapse.  A pre-existing class
   with even $v_w(z)$ and $w\mid s$ remains OPEN uniformly; all $353$
   recorded rows are settled only because their exact local or
   reciprocal certificates cover them.

2. **The false starts were corrected rather than promoted.**  The
   identity
   $$H=(A/4)P=X^2+ALY^2,\qquad L=1-2As^2b$$
   first suggested exact norm matching.  The complete coefficient-linear
   classification requires both $2b=-ALq^2$ and $(A,L)=1$, and its
   full $2$-adic analysis shows every point misses $\Phi$.  A
   polynomial identity $P=U^2-2bV^2$ already fails at
   $P(0)=(2a^4Z^2)^2A$.  On $\Phi$,
   $v_2(P)=4+4\min(v_2(Z),0)$, so
   $P=\pm2b\,\square$ is empty; $P=A\,\mathrm{Norm}$ forces
   $(P,2b)_2=-1$; and the fixed-$P$ square locus is genus $1$, not
   genus $0$.  The first admissible factor ansatz
   $$L=-A\rho^2,\qquad H=(X-A\rho Y)(X+A\rho Y)$$
   is real, but it leaves scoped bad Hilbert symbols and genus-$4$
   residual curves.  Fixed twists merely relocate the character.
   Letting $a$ move also corrected an overstatement: on the
   $a(0)=0$ branch, $Q^4\Vert H$ and the certified squareclass
   representative has degree $18$ with factor degrees $2,16$.
   There is no universal minimum-$21$ claim.

3. **The half-sieve theorem and its barrier.**  For the selected
   linear/octic pair, the bad-root density has exact Chebotarev
   dimension $\kappa=1/2$.  Bombieri--Vinogradov gives
   $D=X^{1/2}/(\log X)^B$, and the semilinear lower sieve at
   $z=X^{49/100}$ proves
   $$\#\{t\asymp X:Q(t)\ {\rm prime},\
     \text{no bad root prime}<z\}\gg X/(\log X)^{3/2}.$$
   This is unconditional small-prime cleanliness.  It is not a member
   theorem: $|G(t)|=X^{8+o(1)}$ can retain $16$ large factors, and
   reciprocity permits $0,2,\ldots,16$ bad ones.  The first missing
   estimate is the parity-sensitive sector with
   $p_1,p_2\ge z$, $p_1p_2>D$, both signs bad and both valuations odd.
   Even Elliott--Halberstam leaves eight factors; Kao reaches only
   total $P_{12}$ after an AP adaptation, while reciprocity needs
   $R_{\rm bad}\le1$.  Schinzel H was not removed.

4. **The fixed surface explains the algebraic wall.**  The conic bundle
   has ten geometric degenerate fibres in closed degrees $1+8+1$ with
   splitting classes $(A,2\theta,A)$.  Since
   $N(2\theta)=256$, Faddeev leaves the sole relation
   $e_0=e_\infty$.  On the selected nonsquare scope,
   $$\operatorname{rank}_{\rm nonsplit}=10,\qquad
     \operatorname{Br}(X)/\operatorname{Br}(\mathbb Q)
       =\mathbb Z/2\langle(A,b)\rangle.$$
   The hypothetical $2\theta$-square case would have rank $2$ and
   trivial Brauer quotient, but L23 rules it out on the selected
   protocol.  Harpaz--Wei--Wittenberg first fails at rank $\le2$;
   Browning--Schindler fails at rank $\le3$; Shute first fails at
   residue degree $8$; and the non-split degree-$8$ closed fibre lies
   outside the applicable HSW hypothesis.  No checked unconditional
   fibration theorem supplies the member.

5. **The diagonal looked like a five-count before arithmetic was
   imposed.**  Identifying the square-branch hyperbola coordinates with
   the tied variables gives two exact systems.  With $u=r^2$ in
   orientation I and $u=y^2$ in orientation II, their eliminants are
   $$
   \begin{aligned}
   Q_I&=Au^2+(A^2-c^2-16AB)u+16A^2Bs^2-16A,\\
   Q_{II}&=Au^2-(c^2+16AB)u+16A^2B(s^2-1)-16A.
   \end{aligned}
   $$
   In each case $u$ and $A+u$ must both be squares.  Clearing by the
   existing square $(Ab^2D)^2$ yields a finite-flat degree-$8$ cover,
   so its tied image has rank at most $1$ and the formal intersection
   count is $2+(3+1-1)\le5$.  This was the L24 discovery, not yet a
   definition: the image still had to meet $\Phi$.

6. **The dyadic refutation is complete in both orientations.**  On
   $\Phi$, $A\equiv5\pmod{32}$ and $v_2(c)\ge4$.  Orientation I has
   normalized coefficient valuations $(4,0,0)$, forcing
   $v_2(u)=0$ or $4$; neither permits both $u$ and $A+u$ to be squares.
   Orientation II is
   $$u^2-32hu+16e=0,\qquad h,e\in\mathbb Z_2^\times.$$
   It forces $u=4t^2$ with $t$ odd, and division by $16$ gives the
   exhaustive congruence
   $$0\equiv8+2Ab(s^2-1)\pmod{16},$$
   which would require the impossible $v_2(s^2-1)=2$.  Thus both
   diagonal images are uniformly empty on $\Phi$ and the nominal
   five-count is vacuous.  Odd target places are locally soluble in
   orientation II's standard $k=1$ stratum, so the obstruction is
   genuinely at $2$.  The generic cover has Galois closure
   $V_4\wr C_2$ of order $32$, genus-$35$ actual-cell slices, a
   genus-$53$ actual-$z$ slice, and no base-rational section.  The
   $370$-cell, $4{,}026{,}282$-orientation zero-hit search is
   corroboration only; the mod-$16$ exhaustion is the proof.

7. **Audits and the character-count repair.**  HalfSieveAudit,
   FibrationAudit and DiagonalAudit returned SOUND at confidences
   $0.90$, $0.88$ and $0.90$.  Their applied fixes included:
   divisor-bounded (not bounded-multiplicity) BV residue counts; the
   literal identity $2\beta=1/(As^2)$, not $1/A$; no universal
   first-polygon $2$-adic Capell claim; the degree-$18$ correction;
   the non-split non-rational fibre as the HSW failure; one mod-$16$
   proof covering both parities of $s$; and the rule that finite no-hit
   scans prove nothing uniformly.  The L20 character count is now
   recorded in its only correct general form:
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
   The old blanket $(w+1)/2$ count was false outside the
   $w\equiv3\pmod4$ scope.  The final status is unchanged:
   classical Schinzel H remains the sole conjectural input to the
   conditional six-count, and H10/$\mathbb Q$ remains open.

## Session 2026-08-23 — L25–L27: scaling, a reciprocal tie, and a one-piece aligned shear

1. **L25 separated an accidental dyadic wall from the coupled
   architecture.**  The fixed type-I scalings $(2X,\rho)$ and
   $(2X,\rho/2)$ are admissible exactly on even and odd $s$,
   respectively.  The target formula with determinant $-2400$ is on a
   freely chosen unit-$b$ stratum, not an aligned L20 target.  Audit
   corrected its Jacobian label to the $(\rho,X)$ block.  The even real
   bridge sample was also corrected from the canonically rejected
   $b=1$ value to the guarded sample $(a,b,Z)=(1,3,-10)$.

2. **L26 is a new reciprocal six-count tie.**  With $Z=z^2$ and
   $\eta=Z(b+1)/(Db)$,
   $$
   P_{\rm rec}(b)=b^4T(b+b^{-1}),\qquad \deg T=4.
   $$
   On $\tau^\dagger$, $P_{\rm rec}$ is a $\mathbb Q_2$-square for every
   guarded $\Phi$ base, standard W1 targets lift, and
   $$
   \left(\frac{2b}{p}\right)
   =\left(\frac{2(u+2)}p\right).
   $$
   The character therefore descends to a quartic trace algebra.  At this
   checkpoint the field and nontrivial quadratic-extension language was
   conditional on $T$ being irreducible and $2(\theta+2)$ being
   nonsquare; L28 (2026-08-24) subsequently proves both uniformly.  The
   one $(a,z)=(1,3)$ modular certificate remains EVIDENCE only for the
   distinct reciprocal lift.

3. **L27 is a genuinely non-diagonal one-piece local solution.**  The
   fixed shear
   $$
   (y,r)=(2X+28\rho,\ sX+\rho)
   $$
   has determinant $2(1-14s)$ and one strong-Hensel argument covers
   every dyadic $s$-parity.  On the standard aligned target stratum,
   putting $x=\rho^2$ reduces existence to the character pattern
   $$
   \chi(x)=\chi(195x^2-57x+4)=+1,\qquad
   \chi(195x^2-56x+4)=-1.
   $$
   Weil bounds give at least
   $(w-11\sqrt w-26)/8>0$ smooth residues for $w\ge197$; both signs were
   exhausted for all $43$ odd primes below $197$.  The guarded real row
   uses $(a,b,Z)=(1,-1,-10)$.  The residual global cover is an exact
   octic in $\lambda=X+\rho$ with leading $225A$ and constant $169A^5$.

4. **Creative alternatives were tested, not promoted.**  The even
   pullback $Z=z^4$ with the $\alpha\sim A$ branch removes both displayed
   L23c target Capell tests, but that branch is itself locally
   incompatible with a standard ramified target; the chain does not
   move.  A parity-aware divisor-spin expansion is exact, but on a
   bad--bad pair reciprocity makes the product character $+1$, reinforcing
   rather than cancelling the pair sector.  A rationally parameterized
   $\Phi$ family reduces another route to an explicit low-degree
   multinorm/weak-approximation problem, whose Brauer obstruction and
   global points remain OPEN.

5. **Audit and scope.**  ReciprocalTieAudit caught three load-bearing
   qualifications: W1 unit guards, the canonical $b\ne1$ bridge domain,
   and conditional trace-field language.  TriangularShearAudit
   independently rederived the octic, strong-Hensel proof, character
   count, zero-weight correction, threshold, L20 compatibility and real
   bracket, with no remaining high-confidence finding.  Both new
   producers replay deterministically.  No global member, five-count,
   unconditional quantifier record, or H10/$\mathbb Q$ solution is
   claimed; classical Schinzel H remains the sole conjectural input to
   the established conditional six-count.

## Session 2026-08-24 — L28: the reciprocal trace field closes at \(2\)

1. **The conditional trace field became uniform.**  On the canonical
   even pullback, put \(t=v_2(Z)\in2\mathbb Z\) and shift \(w=u-2\).
   The trace quartic is
   $$
   4a^8AZ^4w^4-128a^{12}Z^4w^2-32A^3Z^2w+d_0.
   $$
   Its Newton polygon has one half-integral-slope segment, so it has no
   linear factor over \(\mathbb Q_2\).  A quadratic factor would force a
   square root \(Y=\alpha^2\) of the Ferrari resolvent
   $$
   Y^3+2pY^2+(p^2-4r)Y-q^2.
   $$
   The resolvent Newton polygon leaves five valuation cases.  Their
   normalized residues are \(3,2,4,4,2\bmod8\), never zero.  Therefore
   the quartic is irreducible over \(\mathbb Q_2\), hence over
   \(\mathbb Q\), for every canonical L26 parameter.

2. **The bad-sign squareclass also became uniform.**  For a trace root
   \(\theta\),
   $$
   N(2(\theta+2))
   =\frac{64\{A^3D^2+64a^8Z^4(a^4-A)^2\}}
          {a^8Z^4A^2}.
   $$
   Since \(v_2(a^4-A)=2\), the braced expression divided by \(A^3D^2\)
   lies in \(1+2^{10}\mathbb Z_2\).  The norm therefore has squareclass
   \(A\equiv5\bmod8\), so it is nonsquare.  Capell now proves
   \(T(v^2/2-2)\) irreducible uniformly.  The producer replays \(2{,}600\)
   exact unit rows, \(175\) norm rows, and one independent mod-\(41\)
   sample (`l28_trace_field.py`; `data/l28_trace_field.jsonl`).

3. **The first global follow-up was pruned exactly, not by a scan.**
   Writing L27's even eliminant in \(m=\lambda^2\), its coefficient
   linear in \(B=2b\) factors as
   $$
   -16ABm(m-A)\bigl((s+1)^2m-(s-1)^2A\bigr).
   $$
   The factors \(m=0\) and \(m=A\) are invalid or nonzero on the full
   eliminant.  The last factor would make \(A\) a square when
   \(m=\lambda^2\); the \(s=\pm1\) endpoints do not escape.  Thus the
   natural linear-\(B\) cancellation section is impossible.  This is
   only an ansatz no-go, not a theorem that the L27 cover has no rational
   point (`l27_triangular_shear.py::b_coefficient_no_go`).

4. **The independent multinorm idea remains pre-theorem.**  The checked-in
   record contains only the verbal reduction in the preceding session,
   not an explicit torsor equation, birational map and inverse, splitting
   fields, or Brauer evaluation.  Without those data there is no
   reproducible claim to promote.  The exact prerequisite is an explicit
   low-degree model followed by either a rational parameterization or a
   verified Brauer--Manin/weak-approximation calculation.  It remains a
   creative route, not evidence for a member.

5. **Independent audit found no defect.**  TraceFieldAudit rederived the
   shifted coefficients, both Newton polygons, the necessity of the
   square Ferrari-resolvent root for every quadratic factorization, all
   five residue exclusions, the norm squareclass, and the Capell
   inference.  Verdict: SOUND, confidence \(0.99\), no findings
   (`agent://TraceFieldAudit`).

6. **Strict consequence at the L28 checkpoint (squareclass clause
   superseded by L29 below).**  L28 removed trace irreducibility and
   \(2(\theta+2)\) from the blocker list but did not then decide
   \(\theta^2-4\).  It did not produce a globally good \(b\), solve the
   parity sector, improve the count, remove Schinzel H, or solve
   H10/\(\mathbb Q\); those latter statements remain current.

## Session 2026-08-24 — L29: reciprocal lift classified and the dyadic input compressed

1. **The distinct reciprocal lift now has an exact four-stratum law.**
   For a square \(Z\), \(t=v_2(Z)\) is even and
   \[
   \theta^2-4\in\mathbb Q_2(\theta)^{\times2}
   \iff t=-2\ \text{or}\ t\ge2.
   \]
   The positive stratum is strong Hensel applied to
   \(w^2(1+4/w)\), \(w=\theta-2\).  At \(t=0\) and \(t\le-4\), the
   shifted reciprocal octic has one slope \(-1/4\); every possible
   factorization would be two Eisenstein quartics, but exhaustive
   necessary congruences modulo \(16\) and \(64\) have respectively
   \(2048\) and \(524{,}288\) candidates and zero survivors.

2. **The exceptional \(t=-2\) split is proved by an ordinary resultant.**
   For
   \[
   S(x)=P_{\rm rec}(1+x)/\operatorname{lc}(P_{\rm rec}),\qquad
   \phi=x^4+6x^3+6,
   \]
   the polynomial
   \(\mathcal R(Y)=\operatorname{Res}_x(S,Y-\phi)\) has exact Newton
   vertices
   \[
   (0,21),\quad(4,10),\quad(8,0).
   \]
   Its two slopes are \(-11/4\) and \(-5/2\).  Were \(S\) irreducible,
   every conjugate \(\phi(\sigma\alpha)\) would have the same valuation;
   the resultant would have one slope.  Hence \(S\) is reducible, and
   its first polygon forces \(4+4\).  This is equivalent to
   \(b^2-\theta b+1\) splitting in the quartic trace field.

3. **Audit changed the proof, not the claim.**  A first independent
   pass correctly objected that the original higher-Newton-polygon
   presentation had not stated its key-polynomial/residual-type
   hypotheses.  Rather than defend that presentation, the proof and
   producer were replaced by the elementary resultant argument above.
   A second independent pass validated the embedding-invariance
   contradiction and both slopes (confidence \(0.82\)).  Its modular
   caveat is met: the producer enumerates every possible first monic
   Eisenstein quartic residue; an exact \(\mathbb Q_2\)-factorization
   would reduce to one of them.

4. **A new rational pullback removes every adverse dyadic stratum.**
   \[
   Z=\frac{8z^2}{1+z^2}
   \]
   has no rational denominator zero and satisfies
   \(v_p(Z)>0\iff v_p(z)>0\) at every odd prime.  At \(2\), its valuation
   is \(3+2v_2(z)\), \(2\), or \(3\) according as \(v_2(z)>0\), \(=0\),
   or \(<0\).  Thus \(v_2(Z)\ge2\) on every nonzero rational fibre.
   L28's trace and bad-sign proofs extend to these \(t\ge2\) values,
   while L29 makes the distinct reciprocal lift locally \(4+4\).
   Soundness, odd targets and the formal witness count are unchanged.

5. **The moving norm field can be frozen, but parity survives.**  On
   \(b=w\rho^2\),
   \[
   (M,2b)_v=(M,2w)_v,\qquad
   U^2-2wV^2=P_{\rm rec}(w\rho^2)W^2.
   \]
   The slice has no generic section because its \(\rho=0\) fibre has
   squareclass \(A\) and \((A,2w)_w=(A|w)=-1\).  It is nevertheless
   active: eight named target fibres are globally soluble.  On
   \((w,a,Z)=(3,5,9)\), \(24\) of the \(33\) positive odd
   \(w\)-adic units \(\rho<100\) were exactly decided: bad-place counts
   \(0,2,4\) occur \(11,12,1\) times; \(9\) refusals are excluded.  This
   is evidence
   for a fixed-field attack and simultaneous evidence that obstruction
   pairs remain.

6. **Two more natural L27 sections are now closed.**  Uniformly on
   \(\Phi\), \(\lambda=\pm1,\pm A\) gives
   \(C_2\equiv128\bmod256\).  The coordinate choice \(y=0\) forces
   \((2a)^2-195\rho^2=-1\), which has no \(\mathbb Q_3\)-point.  Together
   with the linear-\(B\) cancellation theorem, these are three scoped
   no-gos, not a no-point theorem for the full cover.

7. **The newest source scan changes routing, not status.**  Diao's
   random-binary-form Chowla/norm theorem (arXiv:2506.18065) is an
   almost-all coefficient result; Wang's 2026-08-17 dynamical Chowla
   theorem on average (arXiv:2608.16108) is also averaged.
   Loughran--Matthiesen (1904.12845) first misses the non-rational
   non-split fibre, and Shute (2209.08949) first misses the generic
   high-degree fixed-field slice.  None supplies the missing pointwise
   member theorem.

8. **Strict consequence.**  The local trace/lift questions are closed.
   The first reciprocal target is now a global member on the frozen
   field slice or the original two-large-divisor estimate.  No
   unconditional count improvement, removal of Schinzel H, or solution
   of H10/\(\mathbb Q\) follows
   (`l29_reciprocal_frontier.py`;
   `data/l29_reciprocal_frontier.jsonl`; `THEOREMS.md`, L29).

## Session 2026-08-24 — L30: the coupled equation drops from degree eight to four

1. **The unconventional constant-witness section works locally.**  On
   the L24--L27 square branch, set
   \[
   y=2,\qquad r=sX+\rho,\qquad X^2-\rho^2=A.
   \]
   With \(m=(X+\rho)^2\), the tied equation becomes the exact quadratic
   \[
   H_2(m)=-(c^2-4A)(m-A)^2-16ABQ_s(m)-64Am=0,
   \]
   where
   \[
   Q_s(m)=(s+1)^2m^2-2A(s^2+1)m+(s-1)^2A^2.
   \]
   Hence the \(\lambda=X+\rho\) cover is even quartic.  Its leading
   coefficient has valuation \(2\) on every \(\Phi\)-base, so the degree
   drop from L27's octic is genuine there.

2. **The dyadic theorem is one ordinary Hensel step.**  After dividing
   by \(A\), put \(m=1+8t\) and normalize by \(2^8\).
   The value at \(t=0\) is even because \(Q_s(1)\) is divisible by
   \(16\), while the derivative is odd because \(Q_s'(1)\) is even and
   \(A-1=4a^2\).  Thus \(m\in1+8\mathbb Z_2\) is a square and the cover
   has a \(\mathbb Q_2\)-point uniformly.  The artifact replays \(2048\)
   finite residue rows.

3. **Every aligned odd target is also solved.**  Modulo the target,
   \[
   \rho^2=4,\qquad X^2=A+4.
   \]
   For \(w\ge5\), the required residue conditions are
   \[
   \chi(1+4a^2)=-1,\qquad \chi(5+4a^2)=+1.
   \]
   A squarefree quartic Weil bound gives
   \[
   N_w\ge (w-3\sqrt w)/4-2>0\qquad(w\ge23),
   \]
   and exact rows close \(5,7,11,13,17,19\).  At \(w=3\), the exact
   aligned fibre \((a,b,z)=(5,3,3)\) has
   \(v_3(H_2(1))=3>2v_3(H_2'(1))=2\), so strong Hensel gives
   \(m\in1+9\mathbb Z_3\), a square.  The control \((1,3,3)\) is empty
   because \(v_3(\operatorname{disc}_m H_2)=1\).  Thus target existence
   holds at \(3\), but is not uniform in a preselected \(a\).

4. **The real place is nonempty, but the rational point is not found.**
   The guarded row \((a,b,Z)=(1,-1,1/8)\) gives positive rational
   \(\rho^2\) and \(X^2\).  The exact bounded probe checked
   \(13{,}224\) guarded rows and found no rational \(\lambda\); this is
   EVIDENCE only.  The new global target is precisely
   \(H_2(\lambda^2)=0\) with all remaining controlled places.

5. **Varying the whole L11c square-branch pencil cannot manufacture a
   generic section.**  For
   \[
   \ell_d=(A-d^2)/(2d),
   \]
   every branch has target norm field \(\mathbb Q(\sqrt{2b})\).
   At \(b=0\), the cleared value has even valuation and leading
   squareclass \(A\) for every rational \(d(b)\).  A norm
   \(Y^2-2bR^2\) of even valuation has square leading unit.  This
   contradiction covers rational, Laurent, pole and nonlinear generic
   sections, while leaving specialized members open.
   The marginal choice \(d=2a\), \(\ell_d=1/(4a)\), also fails
   uniform dyadic admissibility: the exact guarded row
   \((a,b,z)=(1,-15,-9)\) has \((M_d,2b)_2=-1\).  This deprioritizes
   that uniform tie without excluding specialized fibres.

6. **The fixed L27 cancellation no-go generalizes to every
   \(B\)-independent linear shear.**  Its \(B\)-coefficient is
   \[
   -16ABmH_L(m),
   \]
   with
   \[
   H_L=(r+t)^2m^2-2A(t^2-r^2+2s^2)m+A^2(r-t)^2.
   \]
   The discriminant is
   \(16A^2s^2(s^2+t^2-r^2)\); for \(s\ne0\), every rational root has
   \(m/A\) square, impossible because \(m\) is square and
   \(A\equiv5\bmod8\).  The \(s=0\) residual edge remains.

7. **Trace descent is exact but does not bypass the reciprocal lift.**
   Return here to L26's square pullback \(Z=z^2\) (or L29's compressed
   \(v_2(Z)\ge2\) stratum).  Then
   \[
   M=T(u)/(DA)^2,\qquad [2b]=[2(u+2)],\qquad u=b+b^{-1}.
   \]
   The base change \(u+2=w\rho^2\) gives the fixed-field degree-\(8\)
   auxiliary
   \[
   U^2-2wV^2=T(w\rho^2-2)W^2.
   \]
   L28's norm squareclass proves the degree-\(8\) polynomial
   irreducible over \(\mathbb Q_2\).  But rational \(b\) separately
   requires
   \[
   \lambda^2=w(w\rho^2-4).
   \]
   Parameterizing this conic recovers exactly \(b=wq^2\), so omitting it
   would merely hide L29's parity wall.

8. **Latest research gives techniques, not a theorem for this fixed
   family.**  Frei--Sofos (arXiv:2604.07047) introduce an analytic
   Hilbert-symbol detector and level lowering, but their main result is
   \(L^2\) over the full coefficient box, not pointwise on the fixed
   H10 form and prime progression.  The degree-\(8\) conic-bundle work
   arXiv:2511.17213 makes the trace auxiliary geometrically relevant,
   but its direct Mestre criterion fails in both diagonal orientations
   and the reciprocal-lift conic remains.  arXiv:2607.25287 kills
   Brauer--Manin obstructions only after suitable finite extensions.
   The first analytic estimate is still L23b's two-large-bad-divisor
   sector.

9. **Strict consequence and audit.**  L30 is a genuine algebraic/local
   advance: the best coupled equation is now degree \(4\), not \(8\).
   Independent review rederived the algebra, dyadic Hensel, \(w\ge5\)
   character count, branch/shear rigidity, trace descent, lift
   parameterization and \(d=2a\) counterexample.  Its initial \(w=3\)
   objection was retracted after the exact \((5,3,3)\) fibre was checked;
   its surviving square-pullback scope correction is incorporated in
   L30g.  L30 is not a rational-point theorem, member theorem, count
   improvement, removal of Schinzel H, or solution of H10/\(\mathbb Q\)
   (`l30_quartic_frontier.py`;
   `data/l30_quartic_frontier.jsonl`; `THEOREMS.md`, L30).

## Session 2026-08-24 — L31: one reciprocal member is proved; the other walls become exact

1. **The fixed-field slice now has one fully proved global member.**
   Take
   \[
   (w,a,z,Z,q_0,b)=(13,3,13,169,1,13).
   \]
   The separate reciprocal lift is explicit:
   \[
   \rho=q_0+\frac1{wq_0}=\frac{14}{13},\qquad
   \lambda=wq_0-\frac1{q_0}=12.
   \]
   Thus \(b+b^{-1}=170/13\) and
   \(\lambda^2=13(13\rho^2-4)=144\); no trace-base point is being
   smuggled through the lift.

2. **The norm certificate is exact rather than scan evidence.**  The
   bridge data are
   \[
   A=37,\quad D=-3\cdot83\cdot1033,\quad
   N_g=-2^4 3^6 47,
   \]
   \[
   c=\frac{277941456}{3172343},\qquad
   \eta=-\frac{182}{257217},\qquad
   \delta=-\frac{324}{37}.
   \]
   The norm value is
   \[
   M=\frac{2^4Q}{3^2 37^3 83^2 1033^2},\qquad
   Q=14082426920623718389.
   \]
   A recursive Pocklington chain proves \(Q\) prime.  The only possible
   symbol places are \(2,3,13,37,83,1033,Q,\infty\), and every
   \((M,26)_v\) is \(+1\).  Hasse--Minkowski therefore gives the
   rational norm point.  This promotes the \(w=13\) L29 evidence row;
   the other seven rows remain at their recorded finite scope.

3. **The general reciprocal problem gets harder, not easier.**  For
   \(q_0=r/s\), clearing the mandatory slice gives
   \[
   G_w(r,s)=s^{16}P_{\rm rec}(wr^2/s^2).
   \]
   It is a sign-decorated degree-\(16\) binary-form norm problem.
   Primes belonging only to \(q_0\) enter \(M\) with even valuation away
   from fixed support, so they cannot cancel emergent bad places.
   Odd-multiplicity primes of \(A\) impose the fixed condition
   \((2w\mid p)=+1\).  One exact row is not a per-target theorem.

4. **The constant-two quartic has a rational point immediately before
   the cube pullback.**  On \(a=1\), the bridge discriminant for constant
   \(c\) is
   \[
   q^2=5+\frac{64}{5c}-\frac4c
   \left(\frac{(b-1)^2}{b}\right)^2.
   \]
   It is a polynomial square in the reciprocal coordinate exactly at
   \(c=-64/25\).  Then
   \[
   b=\frac{2101}{25000}
   -\frac{2\lambda^2}{(\lambda^2-5)^2}.
   \]
   At \(\lambda=3\),
   \[
   (b,Z)=\left(-\frac{3253}{3125},
   \frac{2033125}{6101423}\right)
   \]
   satisfies every bridge guard, \(\Phi\) at \(2\), and
   \(H_2(9)=0\).

5. **The near miss identifies the actual wall.**  Exactly,
   \[
   Z=\frac{5^4\cdot3253}{1009\cdot6047},
   \]
   so it is not a rational cube.  At the shared prime \(3253\),
   \(v(b)=v(Z)=1\) but \(v(c)=0\), whereas L30's target stratum has
   \(Z=z^3\) and \(v_w(c)\ge4\).  Imposing the cube on the two rational
   bridge maps gives
   \[
   d^2=-64z^6+96z^3+64
   \quad\text{or}\quad
   d^2=96z^6-224z^3+64.
   \]
   Both are squarefree genus-\(2\) curves.  Magma V2.29-9
   `RationalPointsGenus2` returned exactly
   \(\{(0:-8:1),(0:8:1)\}\) for each, with completeness flag `true`
   and rank bounds \([0,1]\).  Its intrinsic help defines `true` to
   mean all rational points are returned.  Hence this unique
   constant-\(c\) section has only \(z=0\) and violates the guard.
   The exact request/output are recorded as an external-CAS audit; no
   checked-in standalone Magma certificate is claimed.  Other
   constant-two sections remain open.

6. **The four canonical hyperbola parameters are now excluded for
   L30 itself.**  Directly,
   \[
   C_2\equiv128\pmod{256}
   \]
   at \(\lambda=\pm1,\pm A\) on every \(\Phi\)-base.  This is a
   \(16{,}384\)-row residue replay plus a uniform congruence proof.  It
   is a scoped section no-go, not a global no-point theorem.

7. **The two-large-divisor target needs a positive-main asymptotic.**
   For bad root sets \(\mathcal C_p^-\), the pair classes are exactly
   \(\operatorname{CRT}(\mathcal C_p^-\times\mathcal C_q^-)\), and the
   first pair modulus is
   \[
   pq\ge X^{0.98}>D_{\rm BV}.
   \]
   Exact odd valuation requires \(p^2q^2\), leaving product level
   \(X^{1/4-o(1)}\).  More importantly, character expansion leaves the
   positive \(R_pR_q/4\) term, and
   \(\sum_{X^{0.49}\le p\le X^8}1/p\) has constant size.  Therefore
   signed cancellation alone cannot make the unsigned pair sector
   negligible.  The natural target is a fixed-family
   Buchstab/analytic-Hilbert-detector asymptotic with positive zero-bad
   constant, not an AP-error-only bound.

8. **AP1 is the exact Schinzel replacement.**  In one selected aligned
   class, let \(R_{\rm bad}\) count odd-valuation outside places with
   symbol \(-1\).  Every controlled place is \(+1\), so Hilbert
   reciprocity makes \(R_{\rm bad}\) even.  Thus
   \[
   \mathrm{AP1}:\ R_{\rm bad}\le1
   \quad\Longleftrightarrow\quad
   R_{\rm bad}=0
   \quad\Longleftrightarrow\quad H
   \]
   inside the L19--L22 protocol.  L23a allows
   \(R_{\rm bad}=0,2,\ldots,16\), so it does not prove AP1.

9. **Strict consequence.**  The singular request for a fixed-field
   reciprocal member is met at \(w=13\).  A uniform member theorem, a
   cube-compatible controlled \(H_2\)-root, the fixed-family detector
   asymptotic, and AP1 remain OPEN.  Classical Schinzel H remains the
   sole conjectural input; no count or H10/\(\mathbb Q\) status changes
   (`l31_frontier_push.py`; `data/l31_frontier_push.jsonl`;
   `THEOREMS.md`, L31).

## Session 2026-08-24 — L32: automatic-\(\Phi\) target maps

1. **The dyadic base can be parameterized without a membership
   search.**  For every \(t,u\in\mathbb Q\),
   \[
   \sigma(t)=\frac{t}{1+2t^2},\qquad
   a(t)=1+2\sigma(t),\qquad
   b_\kappa(u)=\frac{u^2+\kappa}{5u^2+\kappa}
   \]
   with \(\kappa\in\{1,2,-2\}\) satisfy
   \[
   v_2(\sigma(t))\ge0,\qquad v_2(a(t))=v_2(b_\kappa(u))=0.
   \]
   The proof is a complete three-case valuation split, not the finite
   replay.

2. **The \(a\)-map contains every regular constant-two target.**  With
   \(d=2t^2+1\), \(n=2t^2+2t+1\), the relevant squareclasses are
   \[
   F_5=d^2+4n^2,\qquad F_9=5d^2+4n^2.
   \]
   Their degrees are \(4,4\), their resultant is \(2^{28}\), and their
   odd discriminant primes are \(17\) and \(5,89\).  The indicator
   expansion has Weil error
   \((3+3+7)\sqrt w\), giving
   \[
   N_w\ge\frac{w-13\sqrt w}{4}-4>0\qquad(w\ge211).
   \]
   Exact exhaustion closes every \(5\le w<211\); all primes below
   \(5000\) are a cross-check.

3. **A three-map \(b\)-menu hits every odd target.**  Since
   \[
   (-1|w)(-2|w)(2|w)=1,
   \]
   at least one numerator \(u^2+\kappa\) has a simple root modulo
   \(w\), while its denominator is \(-4\kappa\).  A lift gives
   \(v_w(b_\kappa)=1\).  Combined with item 2 and L30, this supplies one
   explicit rational-map base locally soluble at both \(2\) and every
   \(w\ge5\).

4. **The proportional-trace simplification is now classified.**  The
   ansatz
   \[
   (b-1)^2/b=4ra^2
   \]
   has exactly the branches \(b=t^2/r,r/t^2\) and
   \[
   c=\frac{16a^6Z^2(1-Ar^2)}{AD}.
   \]
   If \(v_w(b)=1\) and \(a\) is a \(w\)-unit, both branches force
   \(v_w(r)=-1\).  Hence every finite fixed \(r\)-menu misses all but
   finitely many targets.

5. **The two linear simplifications are not uniform.**  Although
   \(b=1\pm2a\) gives \(N_g=\pm64a^5\), a target divisor forces
   \(A\equiv2\pmod w\).  Regularity then holds only for
   \(w\equiv5,19\bmod24\).

6. **The global wall moved, but did not close.**  Base selection at
   \(2\) and \(w\) is no longer the issue for L30.  The numerators and
   denominators of the rational maps create uncontrolled emergent
   primes, and no rational \(H_2\)-root satisfying all controlled
   places follows.  The detector asymptotic and AP1 remain open
   (`l32_local_parameter_frontier.py`;
   `data/l32_local_parameter_frontier.jsonl`; `THEOREMS.md`, L32).

## Session 2026-08-24 — L33: the unweighted pair main is a theorem

1. **A positive pair main can be proved before prime weighting.**  For
   \(0<\alpha<\beta<1/2\), let \(\mathcal U_2\) count simultaneous
   bad-root classes for primes \(p,q\in[X^\alpha,X^\beta]\), with
   unweighted \(t\in(X,2X]\).  L23 Chebotarev gives
   \[
   \sum\frac{r_-(p)}p
   =\frac12\log(\beta/\alpha)+o(1).
   \]
   CRT gives \(r_-(p)r_-(q)\) classes and \(X/(pq)+O(1)\) points per
   class.  Since \(2\beta<1\), the total class-counting error is
   \(O(X^{2\beta})=o(X)\).  Therefore
   \[
   \boxed{\mathcal U_2(X;\alpha,\beta)
   \sim\frac18\log^2(\beta/\alpha)X.}
   \]

2. **The first beyond-BV range already contains this main.**  Taking
   \(\alpha=0.49\) and \(0.49<\beta<1/2\) puts every pair modulus above
   \(X^{0.98}\), while retaining a positive order-\(X\) asymptotic.
   Thus the unsigned root-pair main is no longer merely a formal
   expectation.

3. **The remaining estimate is exactly the hard transfer.**  L33 is
   unweighted and uses the mod-\(p\) root oversieve.  Replacing the
   class count by \(\Lambda(Q(t))\), imposing exact odd valuations via
   \(p^2q^2\), and conditioning on L23 small-prime cleanliness are all
   outside the theorem.  Those are precisely the pointwise
   fixed-family detector problem; no negligible unsigned-pair estimate
   is compatible with the proved unweighted main
   (`l33_unweighted_pair_main.py`;
   `data/l33_unweighted_pair_main.jsonl`; `THEOREMS.md`, L33).

## Session 2026-08-24 — L34: a \(2.5\%\) pair-moment margin would prove AP1

1. **Negligibility is stronger than needed.**  Because
   \(R_{\rm bad}\) is even,
   \[
   \mathbf1_{\{R_{\rm bad}=0\}}
   \ge1-\binom{R_{\rm bad}}2.
   \]
   Let \(\mathcal P_2\) use exact bad primes and let
   \(\mathcal P_2^{\rm root}\) use all bad-sign root primes, including
   even valuations.  Since \(R_{\rm bad}\le R_{\rm root}\),
   \[
   \mathcal N_{\rm good}
   \ge\#\mathcal A_z-\mathcal P_2
   \ge\#\mathcal A_z-\mathcal P_2^{\rm root}.
   \]
   Any strict relative root-pair bound below \(1\) proves AP1; exact
   odd valuations are unnecessary for this upper-bound route.

2. **The unrestricted octic envelope is below \(1\) in the BV
   window.**  The large-bad-prime local mass is
   \[
   \mu(\vartheta)=\tfrac12\log(8/\vartheta).
   \]
   Ignoring the additional divisor-product constraint, the unordered
   pair envelope \(\mu(\vartheta)^2/2\) is below \(1\) exactly when
   \[
   \vartheta>8e^{-2\sqrt2}=0.4728459724\ldots.
   \]
   At L23's \(\vartheta=0.49\),
   \[
   \mu^2/2=0.9749604961\ldots.
   \]

3. **The exact sufficient analytic theorem now has a constant.**  A
   pointwise conditioned estimate
   \[
   \mathcal P_2^{\rm root}(X)
   \le(0.9749604961\ldots+o(1))\#\mathcal A_z(X)
   \]
   would imply AP1, intermediate H, and the six-variable theorem
   without classical Schinzel H.

4. **The estimate is not proved.**  The \(50/49\) semilinear-sieve
   ratio is too close to its lower limit to give this sharp constant,
   and pair moduli start at \(X^{0.98}\).  The root oversieve avoids
   \(p^2q^2\), but its prime-weighted moment after small-prime
   conditioning is unavailable.  L33's positive unweighted main also
   rules out an unsigned cancellation-to-zero strategy.  The
   unconditional untied
   seven-variable theorem remains separate
   (`l34_ap1_pair_threshold.py`;
   `data/l34_ap1_pair_threshold.jsonl`; `THEOREMS.md`, L34).
