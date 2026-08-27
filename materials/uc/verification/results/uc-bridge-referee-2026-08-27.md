# Independent UC bridge/support-reduction referee report — 2026-08-27

Invocation: read-only `scout` subagent (Read/Glob/Grep only), instructed to
adversarially re-derive rather than trust status labels, and to produce a
counterexample or a named missing lemma for every alleged defect.  Scope was
restricted to the two human links that remained pending ordinary external
review:

- **Link 1** — the sequential entropy-to-union-closed proof of
  Corollary `cor:uc` (`paper/main.tex`, statement plus proof).
- **Link 2** — the fold/unfold equivalence, optimization collapse, and
  Theorem B''' support reduction (`paper/main.tex`, Definition `def:pair`
  through Corollary `cor:5param`).

## Verdict on the pre-repair text

- Blocking findings: **none**.
- Proof defects: **none**.  No counterexample and no reversed implication was
  found in either link.
- Eight MINOR exposition defects, all repaired in the current sources.

## Independently re-derived transitions

Link 1: the three-case `s*` bounds `max(p,r) <= s <= min(p+r,1)`; nonnegativity
and Bernoulli marginals of the four coupling masses; `Pr(A_i or C_i)=s`;
prefix-only dependence of the conditionals, hence uniformity of `A` and `C` by
induction; `(P_i,R_i)` in `Pi(mu_i,mu_i)` so its cost dominates `C(mu_i)`;
means equal to element frequencies by the tower property; the direction of both
conditioning inequalities via `sigma((A u B)_{<i}) subset sigma(A_{<i},B_{<i})`;
coordinatewise legitimacy of Theorem `thm:main`; and the first-coordinate
strictness split at `psi`.  The report also cross-checked
`verification/cambie_bridge_strictness.py` against the manuscript: same
formula, same exact rationals, same positive bracket, honest exclusion of the
derivative step.

Link 2: both directions of fold/unfold with cost and marginal preservation;
weak-* closedness of `Pi(mu,mu)` and attainment of the minimum; nonempty
compact fibers; both inequalities of the `inf F = inf Phi` collapse; Steps 1-3
compactness/continuity; the concavity **sign** in Step 4; the direction of the
Bauer application in Step 5; the one-moment extreme-point support bound in
Step 6; and the `cor:5param` equivalence in both directions.

## Findings and repairs

| ID | Location | Finding | Repair |
|---|---|---|---|
| E1 | `cor:uc` strictness | "constant `u`" asserted without reason | The proof now says coordinate 1 has an empty prefix, so `P_1,Q_1,R_1` equal the deterministic `u=|F(1)|/|F|`, positive by the choice of element 1. |
| E2 | `cor:uc` strictness | `g'(u)<0` omitted its enabling inequality | Now derived: `2u-u^2=1-(1-u)^2` increases, `psi^2-3psi+1=0` gives `(1-psi)^2=psi`, so `2u-u^2>=1-psi>1/2` and `h'(2u-u^2)<0`, while `u<=t_cert<1/2` gives `h'(u)>0`. |
| E3 | `cor:uc` coupling | the "three-case identity" was cited but never stated | The three exhaustive cases for `s*` are now written out before the bound is used. |
| E4 | Theorem B''' Step 2 | pure tensors were called a subalgebra | Corrected to the linear span of products, which is the unital separating subalgebra Stone--Weierstrass needs. |
| E5 | Theorem B''' Step 2 | continuity of `Phi_alpha` was routed through `eq:Fid` | Continuity is now taken directly from the definition: `L`, `int c dnu`, and `E` are each weak-* continuous. |
| E6 | Theorem B''' Step 4 | "PSD gives convexity" was asserted | The cross-term derivation is now explicit, using `<W~,eta x eta> = ||int v deta||^2 >= 0` at `eta=mu_0-mu_1`. |
| E7 | Margin section | sink degeneracy conflated `Phi_alpha` on `P(Delta)` with `F` on `P([0,1])` | The two statements are now separate, with the pair-orbit version stated for `nu` supported in `{0,1}^2`. |
| E8 | `cor:uc` statement | the empty family was not excluded | The statement now reads "finite nonempty union-closed family". |

## Machine control added in the same pass

`verification/entropy_bridge_exhaustive.py` rebuilds the Link-1 construction
from the corollary statement alone and exhausts every nonempty family on
`n <= 4` coordinates: 65,808 families, 1,631,880 coupled prefix states.  Exact
rational checks cover the coupling masses, both marginals, the `s*` OR
identity, prefix-by-prefix uniformity, `law(P_i)=law(R_i)` with mean equal to
the element frequency, and the chain rule.  256-bit Arb certifies both
data-processing inequalities, each enclosure being provably nonnegative or a
structural tie inside `+/-2^-200`.  All 5,096 enumerated union-closed families
satisfy the `t_cert` conclusion, the extreme case being exactly `1/2`.
Report: `results/entropy-bridge-exhaustive.json`,
canonical SHA-256 `6681f8faf13d9344f8e7bf3a5d7785ed289f1e82c5e9badd1e7967c03e3ebfef`.

Mutation testing (`verification/test_entropy_bridge_exhaustive.py`, 4/4 pass)
shows the control is not vacuous: dropping the `max` branch of `s*`, forcing a
constant prefix probability, or inflating the union bit each makes it reject.
Replacing the `1/2` clip by another admissible cap does **not** make it reject,
because such a coupling is still valid; the identification of the cost function
with Cambie's `s*` is the separate three-case proof in `bridge_uc.py`.

Beyond the exhaustive range the same verifier checks seeded deterministic
samples on `5 <= n <= 8`: 7,650 uniform families contributing 1,682,137
further coupled prefix states, and 5,300 random union-closed families with up
to 146 members.  All pass, and the smallest maximal frequency observed among
them is again exactly `1/2`.  Total runtime 110 s on one core.

## Residual assumptions

Model review is not journal peer review.  This audit inherits Theorem
`thm:main` itself (Margin Lemma plus the Arb certificate), the Theorem A'
identities `Q=2BL-E` and PSD `W`, the Arb primitives, and the bibliographic
accuracy of the Bauer/Winkler/Pinelis/van Neerven citations.  The
`cambie_bridge_strictness.py` value was cross-checked for formula and
constants, not re-executed inside the read-only audit; it is separately
machine-verified in `results/cambie-bridge-strictness.json`.
