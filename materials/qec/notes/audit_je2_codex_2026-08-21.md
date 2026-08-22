# Independent proof audit — Theorem J-E′′ fixed-point structure (2026-08-21)

- **Scope**: Theorem J-E's per-class premise as *used* (not re-derived) and
  the stable-image lemma ($S = e_I M = I^\infty M$), plus the trichotomy,
  the demote-fraction formula, the chain-computability claim, and the
  ideal-splitting premise of `notes/theorem_je2_demote_trichotomy.md`.
- **Auditor**: codex CLI `0.147.0`, `codex exec --sandbox read-only`
  (independent model; 27,141 tokens; 128 s). `fable` delegation per the
  standing preference was attempted first and **unavailable**: OAuth session
  expired and the Keychain refresh requires an interactive unlock; codex
  served as the independent model instead. Prompt was self-contained, no
  repository access, no credentials shared.
- **Verdict (verbatim)**: "FLAWED (AUDIT-6 omits the required 'nonzero
  classes' sampling assumption; AUDIT-7 falsely says rank equality is
  insufficient despite containment being automatic; AUDIT-2 must use the
  local element $i_\beta$, not assert that a global $1-i$ is a unit)."
- **Findings (verbatim)**:

> - **AUDIT-1 — Correct.** If $y=iy$, projection gives $y_\beta=i_\beta y_\beta$. Conversely, choose $i_\beta\in I_\beta$ for each factor and set $i=\sum_\beta i_\beta\in I$; orthogonality of the idempotents gives $iy=y$. Thus a single global $i$ exists.
> - **AUDIT-2 — Conclusion correct; proof needs local notation.** In a local factor, $y_\beta=i_\beta y_\beta$ for some $i_\beta\in I_\beta\subseteq\mathfrak m_\beta$. Hence $1_\beta-i_\beta$ is a unit and $y_\beta=0$. No single coefficient shared across factors is needed. If the displayed $i$ means a global element of $I$, claiming that global $1-i$ is a unit is unjustified; only $1_\beta-i_\beta$ is known to be a unit.
> - **AUDIT-3 — Correct.** The product-ideal argument is valid. More directly, $I_\beta=e_\beta I=R_\beta=e_\beta R$ implies $e_\beta\in e_\beta I\subseteq I$.
> - **AUDIT-4 — Correct.** Proper factors contribute only zero by AUDIT-2. On full factors, $e_\beta\in I$, so every element is fixed. Therefore $S=e_I M$.
> - **AUDIT-5 — Correct.** Module and ideal powers decompose factorwise. In a commutative local Artinian ring, $\mathfrak m_\beta$ is nilpotent; hence every proper $I_\beta\subseteq\mathfrak m_\beta$ is nilpotent. Taking the maximum exponent over finitely many factors leaves precisely $e_I M$. Finite generation of $M_\beta$ is not needed for this particular nilpotence argument.
> - **AUDIT-6 — Correct only with an unstated sampling convention.** The trichotomy is valid when $M\neq0$. The displayed fraction is the fraction among nonzero classes $M\setminus\{0\}$. Among all elements of $M$, the fraction would instead be $1-\frac{2^{\dim S}}{2^{k_P}}$. Minimal fix: explicitly say "demote fraction among nonzero quotient classes."
> - **AUDIT-7 — Computation correct, but the rank claim is overstated.** Bilinearity shows that a basis of $I$, together with a basis of $F_r$, spans $IF_r$. Also $F_{r+1}\subseteq F_r$ holds automatically because $I^{r+1}\subseteq I^r$. Consequently, for an exact implementation, rank equality alone already implies $F_{r+1}=F_r$. Explicit containment checking is a useful implementation sanity check, but not mathematically necessary.
> - **AUDIT-8 — Correct.** For $J\triangleleft\prod_\beta R_\beta$, let $J_\beta$ be its projection. Given $a_\beta\in J_\beta$, choose $j^{(\beta)}\in J$ with that component. Then $e_\beta j^{(\beta)}\in J$ has only component $a_\beta$. Summing over the finitely many factors proves $J=\prod_\beta J_\beta$.
> - Additional checks: characteristic 2 matters for the binary-vector-space dimensions and counts $2^d$, not for AUDIT-1–5 or 8; the presentation by $x^\ell-1, y^m-1$ matters only for producing the finite commutative ring and the shift-module setup. Quotient modules decompose normally: $M=\bigoplus_\beta e_\beta M$; being a quotient of $\ker H_X$ causes no problem. There is no separate commuting $I$-action: $I$ acts by restriction of the $R$-module action; associativity and the central orthogonal idempotents suffice.

## Applied fixes (all three, before "proved" was retained)

1. **AUDIT-2**: local step rephrased with the *local* element: "$y_\beta =
   i\,y_\beta$ for $i\in I_\beta$ gives $(1-i)y_\beta=0$ **in $R_\beta$**,
   where $1-i$ is a *local* unit" (was: global-$(1-i)$-unit wording).
2. **AUDIT-6**: the fraction formula is now explicitly **among nonzero
   quotient classes**, with the all-elements variant recorded.
3. **AUDIT-7**: stop rule reworded in note, code comments, and ledger:
   containment $F_{r+1}\subseteq F_r$ is *automatic* (ideal powers nest);
   rank equality certifies stabilization; the containment assertion remains
   in EXP-052 as an implementation tripwire, not as a mathematical
   requirement.

No theorem-statement change beyond the $k_P>0$ hypothesis (added per a
review flag: at $M=0$ cases (i) and (ii) coincide vacuously).

## Audit prompt (verbatim)

```text
You are auditing a pure-mathematics proof in quantum error correction (commutative algebra over GF(2)). Be adversarial: find any gap, unstated assumption, or false step. If a step is right, say why briefly. End with a verdict line: VERDICT: SOUND / VERDICT: FLAWED (list flaws).

SETUP (all claims to audit are labeled AUDIT-n).

Let R = GF(2)[x,y]/(x^l - 1, y^m - 1) (char 2). R is a finite commutative ring, hence Artinian: R = prod_{beta} R_beta, a finite product of local Artinian rings (idempotents e_beta). WARNING accepted as true: the coarser blocks GF(2)[x]/(p^a) tensor GF(2)[y]/(q^b) need NOT be local (e.g. GF(4) tensor GF(4) = GF(4) x GF(4)); the proof must work at LOCAL factor granularity beta and must never rely on the coarse blocks being local.

Let H_X = [A B] (binary circulant matrix built from polynomials A,B in R; rows are the lm translates of the single row-vector (A,B)); notice ker H_X is shift-stable because right multiplication by any polynomial commutes with the circulant action. Let S_Z = rowspace(H_Z) where H_Z = [B^T A^T] rows are all lm translates of (b-bar, a-bar) (bar = reciprocal polynomial x->x^{-1}, y->y^{-1}); hence S_Z is shift-stable as an R-submodule, so M := ker H_X / S_Z is a finitely generated R-module, dim_GF(2) M = k_P. Let I = L_pre := {lambda in R : lambda A = 0 and lambda B = 0} (an ideal of R; equivalently the left-null space of the stacked matrix [A;B] read as polynomial annihilation).

EXTERNAL PREMISE (do not re-derive; audit only its USE): for y in M (a quotient class), the "single-row syzygy demotion" of class y is decided by: demote(y) iff y is NOT in I*y (the set of products i*y for i in I). This is Theorem J-E/Demote-test as implemented (rank-based membership test over GF(2) in the quotient coordinates).

Define the demote fixed set S = { y in M : y in I y }.

AUDIT-1 (componentwise membership): y in I y iff for every beta, y_beta in I_beta y_beta, where I_beta = e_beta I is the component ideal and y_beta = e_beta y.

AUDIT-2 (local-factor step): in a local Artinian ring R_beta with maximal ideal m_beta, if I_beta is a PROPER ideal, then y_beta in I_beta y_beta implies y_beta = 0. (Claimed proof: y_beta = i y_beta for that i in I with (1-i) y_beta = 0; (1-i) is a unit because i in I_beta subset m_beta; hence y_beta = 0. Question: does "y in Iy" give a SINGLE i working for all beta simultaneously, and does the argument need it?)

AUDIT-3 (idempotent membership): if I_beta = R_beta then e_beta in I. (Claimed proof: ideals of a finite direct product split componentwise as sets, i.e. every ideal of prod R_beta is a product of ideals; hence I = prod I_beta; therefore when I_beta = R_beta, the element (0,...,1,...,0) = e_beta is in I.)

AUDIT-4 (fixed set is idempotent image): S = e_I M where e_I = sum over full factors (I_beta = R_beta) of e_beta.

AUDIT-5 (fixed set is stable image): I^infinity M (stabilization of the descending chain F_{r+1} = I F_r, F_0 = M) equals e_I M. (Claimed proof: I^r M = prod I_beta^r M_beta; proper I_beta subset m_beta is nilpotent-acting on the finitely generated module M_beta, so eventually 0; full I_beta^r M_beta = M_beta for all r.)

AUDIT-6 (trichotomy and fraction): ASSUME k_P > 0 (M nonzero; at M = 0 cases (i) and (ii) coincide vacuously): then exactly one of (i) S = M, (ii) S = 0, (iii) 0 != S != M; in case (iii) demote fraction = 1 - (2^{dim S} - 1)/(2^{k_P} - 1), where dim = GF(2)-dim.

AUDIT-7 (chain computation validity): dim S can be computed inductively in GF(2) coordinates: take any GF(2)-basis B of I (as a GF(2)-linear space) and GF(2)-coordinates for M; then GF(2)-rowspan{ i * y : i in B, y in F_r } = F_{r+1} as a GF(2)-subspace; termination must be recognized by subspace containment F_{r+1} <= F_r PLUS rank equality (rank equality alone is insufficient).

AUDIT-8 (ideal splitting premise): the proof uses "every ideal of a finite direct product of rings is a product of ideals". Confirm or refute (elementary).

Also check: (a) does anything depend on characteristic 2 or the specific shape x^l-1, y^m-1? (b) does M being a quotient module (rather than submodule of a free module) invalidate any product-decomposition argument? (c) any hidden use of commutativity of the I-action versus the R-action on M?

Be terse, adversarial, exact. If you find a genuine gap, state the minimal fix if one exists.
```
