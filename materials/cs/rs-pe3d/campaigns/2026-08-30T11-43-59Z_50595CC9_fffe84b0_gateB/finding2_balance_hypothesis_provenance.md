# Gate B, finding 2 — the balance hypothesis is present in EVERY formal statement and absent only from the informal restatement

Agent `RsPe3dGateB`, 2026-08-30. Source-verified this session by first-hand
read of both PDFs. This section is independent of the numerical sweep and
stands on its own.

Artifacts and checksums (my own downloads, this session):

| file | sha256 | provenance |
|------|--------|------------|
| TR26-150 PDF | `fb51d7dc84860f115d88ec47e92799e9c949ac90cbe714e8df7d3904c930af87` | `https://eccc.weizmann.ac.il/report/2026/150/download` |
| AI candidate proof PDF | `e2d9b4b61f9ed8589f90e5ed6bf146eaf642a85969c859589ad4ed405778330c` | `https://mitalibafna.github.io/higher-dimension-product-expansion.pdf`, linked from TR26-150's abstract as "(link)" — recovered from the PDF's `/URI` annotation |

Text line numbers below refer to my extraction `pdftotext -layout` of each PDF
(reproducible: same tool, same input hash).

## The four statements side by side

**(1) Conjecture 4.2 — formal, higher-dimensional. Balance PRESENT.**
TR26-150 text lines 740–748, verbatim:

> Conjecture 4.2 (Product Expansion of Punctured Reed–Solomon Codes). For every η ∈ (0, 1), **balance
> parameter β ⩾ 1**, and integer r ⩾ 2, there is ρ = ρ(r, η, β) ∈ (0, 1] such that the following holds. Let q be
> prime, and let S1 , . . . , Sr ⩽ F×_q be multiplicative subgroups of pairwise coprime orders si = |Si | satisfying
>
>                                          β^{−1} ⩽ s_i/s_j ⩽ β    (i, j ∈ [r]).

**(2) Theorem 4.11 — proved, author-verified, 2-dimensional. Balance PRESENT
as $K$, and quantified in the constant.**
TR26-150 text lines 1280–1285, verbatim:

> Theorem 4.11 (Two-dimensional product expansion of subgroup Reed–Solomon codes). Fix ε ∈ (0, 1) and
> **K ⩾ 1**. There is a constant **ρ = Ω(ε^6 /K^3)** ∈ (0, 1] with the following property. Let q be prime, and let
> S1 , S2 ⩽ F×_q be multiplicative subgroups of coprime orders n_i = |S_i | satisfying **K^{−1} ⩽ n1/n2 ⩽ K**.

**(3) Lemma 5.4 — the paper's own parameter selection for its PCP/code
construction. Balance PRESENT, instantiated at $\beta = 4$.**
TR26-150 text line 1632, verbatim fragment (item 2, "Domain sizes"):

> 1/4 ⩽ s_i/s_j ⩽ 4      (i, j ∈ [k]).

and line 861: "By Lemma 5.4, the evaluation domains have pairwise coprime,
**4-balanced** orders. Hence Conjecture 4.2 gives product expansion at least
ρ(i, 1/128, 4) ⩾ ρ(k)." Corollary 8.3 likewise fixes ρ = ρ(k, 1/128, 4).
Section 6 states the dependence outright (line 1604): "our product expansion
conjecture, **whose hypotheses require the side lengths of Ω to be pairwise
coprime and balanced**."

**(4) The AI candidate proof — unverified, higher-dimensional. Balance
PRESENT, and it is the dominant cost in the constant.**
AI PDF text lines 4–8 and 51–56, verbatim (double-column extraction artifacts
retained):

> Theorem. Fix d ≥ 2, 0 < ε < 1, and **K ≥ 1** [...] Let S_i ≤ F×_q have pairwise coprime orders n_i = |S_i |
> satisfying **K^{−1} ≤ n_i /n_j ≤ K for all i, j ∈ [d]**.
> [...] There is a universal constant C_0 ≥ 1 for which one may take
>
>       ρ = (1/(2dK^d)) · ( ε / (C_0 d K)^{10d} )^{3d·6^{d−2}}.

Its Proposition 5 (AI PDF line 252) makes the same hypothesis internal:
"assume that they are **K-comparable**, meaning K^{−1} ≤ n_i/n_j ≤ K for all
i, j ∈ [e]", and the proof's slicing estimate uses it repeatedly
("K-comparability absorbs the change", line ~264).

## The one statement that omits it

**Conjecture 1.4 — informal.** TR26-150 text lines 224–226, verbatim:

> Conjecture 1.4 (Informal). Let q be a prime, k ⩾ 2 and C = (C_1 , . . . , C_k ) be a k-tuple of Reed–Solomon
> codes evaluated on multiplicative subgroups of F⋆_q of pairwise coprime orders, with rate(C_i ) < 1 − ε. Then
> C is ρ-product-expanding with **ρ being only a function of ε and k**.

and the abstract's summary sentence: "We conjecture the analogous statement in
higher dimensions."

## What this establishes

- **The gap is between the informal restatement and every formal statement,
  not between a proved theorem and a bolder conjecture.** Conjecture 4.2 (the
  formal object the paper actually uses in Corollaries 4.3/4.4/4.5 and 8.3),
  Theorem 4.11, Lemma 5.4's parameter selection, and even the unverified AI
  candidate proof all carry a balance/comparability hypothesis. Only
  Conjecture 1.4's informal phrasing "ρ being only a function of ε and k"
  drops it.
- **Consequently the target README's framing needs one correction of emphasis**
  (recorded here, applied in the README current-state section): the README says
  the balance hypothesis is one "which the informal higher-dimensional
  conjecture drops" — that is exactly right, and the sharper statement is that
  the paper's own formal conjecture does NOT drop it, so no correction to
  Conjecture 4.2 is needed on this axis. The correction, if any is wanted, is
  to Conjecture 1.4's informal wording: it should read "ρ being only a function
  of ε, k and the balance parameter β", matching Conjecture 4.2. That is a
  wording correction to an informal restatement, not a mathematical claim
  against the paper.
- **The AI proof's $K$-dependence is explicit and severe.**
  **RETRACTED FIGURE (rule 5, retained inline with cause), 2026-08-30:** I first
  wrote here that at $d = 3$ the AI display gives
  $\rho = \frac{1}{6K^3}\left(\frac{\varepsilon}{(6C_0K)^{30}}\right)^{162}$ and
  hence $\rho \sim K^{-4863}$. **That is wrong**, and so is the independent
  "confirmation" of $K^{-4863}$ I received from `Main` (which read the printed
  exponent as $3^d 6^{d-2} = 162$). Cause: the printed ASCII exponent
  "$3d\,6^{d-2}$" is genuinely ambiguous between $3d$ and $3^d$, and I resolved
  it by guessing rather than by using the proof's own identity.

  **Corrected value, pinned by the AI proof's own recurrence** (AI PDF text
  lines 245, 489–493, 933): the proof defines $\gamma_1 = 1$,
  $\gamma_e = 1/(3 \cdot 6^{e-2})$ for $e \ge 2$, states
  $\gamma_e = \gamma_{e-1}/6$ for $e \ge 3$, and then writes "Since
  $d/\gamma_d = 3d\,6^{d-2}$". So the exponent **equals $d/\gamma_d$**, and that
  disambiguates it with no guessing:
  $\gamma_2 = 1/3 \Rightarrow d/\gamma_d = 6$ at $d=2$, which matches
  $3d\,6^{d-2} = 6$ and **contradicts** $3^d 6^{d-2} = 9$;
  $\gamma_3 = 1/18 \Rightarrow d/\gamma_d = 54$ at $d=3$, matching
  $3d\,6^{d-2} = 54$, not $3^d 6^{d-2} = 162$. The linear reading is therefore
  the correct one, at every $d$.

  Hence at $d = 3$:
  $$ \rho = \frac{1}{6K^3}\left(\frac{\varepsilon}{(3C_0K)^{30}}\right)^{54}, \qquad \rho \sim K^{-(3 + 30 \cdot 54)} = K^{-1623} $$
  up to constants — $K^{-3}$ from the prefactor, $K^{-1620}$ from the bracket.
  Checked exactly in `Fraction` arithmetic (no floats).
  Label: [DERIVED] by substituting $d = 3$ into the displayed formula and using
  the proof's own $\gamma$ recurrence to fix the exponent; the underlying proof
  is UNVERIFIED by the authors, so this is a statement about the claim, never
  about the truth. The qualitative point is unchanged and if anything
  strengthened by being correct: the AI candidate's $\rho$ carries an enormous
  negative power of the balance parameter, so that claim is very far from
  balance-free.
- Evidence labels: all quotes above are [REPRODUCED] from PDFs whose hashes are
  recorded; the $d=3$ substitution is [DERIVED]; nothing here is
  MACHINE-VERIFIED (it is a reading of text, not a computation).

## Escalation status

This section makes **no** claim that any TR26-150 statement is wrong. It
records which statements carry the hypothesis and which does not. Per the
campaign order, the only wording that could bear on a published claim — that
Conjecture 1.4's informal phrasing omits β — was reported to `Main` before
being written into the target README.
