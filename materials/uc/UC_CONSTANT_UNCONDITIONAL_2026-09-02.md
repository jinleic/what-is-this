# The union-closed constant c' = 1 - m is unconditional - 2026-09-02

Claim audited here: **every union-closed family `F` of subsets of a finite set, other than
`{}` and `{{}}`, has an element contained in at least `c' |F|` of its members, where**

    c' = 1 - m = 0.3827090879187350299303129020989726116264...   (certified enclosure +/- 4.95e-46)

is the root-defined constant of Liu, arXiv:2306.08824v1, equations (87)-(90), i.e. `m = p* x*` with
`x*` the unique root in `(0,1)` of `x^4 - 2x^3 + 3x^2 - 1` (Liu's (87)), `p* = h(x*)/h(x*^2)` (Liu's
(88)), and `beta*` from (89)-(90). Liu's Theorem 13 obtains this constant **conditionally** on two
hypotheses (Section V-A positive semidefiniteness, Section V-B minimiser structure). Both are now
bypassed: the inequality his Proposition 3 needs is proved outright for every conditionally i.i.d.
coupling, with the uniform constant `M/m`.

Labels: PROVED / MACHINE-VERIFIED / HUMAN-AUDITED / CITED-DEPENDENCY / COMPUTATIONAL-EVIDENCE / OPEN.

## 1. The two machine-verified inputs

| artifact | claim | status | file sha256 |
|---|---|---|---|
| `uc/verification/results/liu9-h2-twovar.json` | `A >= 0` (and `R >= 0`) on all of `[0,1]^2`, with `beta`, `m` as Arb balls (radius `~2e-51`) around the exact root-defined constants; `x*` isolated to `6e-70` | CERTIFIED_PROVED | `a028ff17ae8492211a2ef9de5721de7f3a0afd8bd45bf7d8d9b77b7c90d2e6eb` |
| `uc/verification/results/liu9-h2-boundary.json` | `Phi(s,t) = h(pi(s,t))^2 - h(pi(s,s)) h(pi(t,t)) >= 0` on all of `[0,1]^2` (parameter-free), hence `B >= 0` for `beta > 0` | PROVED (independently audited, `uc/H2_PHI_AUDIT_2026-09-02.md`) | `e4d6d96df3ae435800cef5739fb1ae45b98febf539ab418df116ceb03bd481c9` |
| `uc/verification/results/liu9-h2-general-lift.json` | two-component general-measure theorem at the exact `(beta*, m*)`; kernel regularity L0; `B >= 0 <= Phi >= 0` (L4); the exact-constant discipline (balls from the 70-digit bracket of `x*`, cross-checked against `twovar`'s definitions) | PROVED (independently audited; first version's parameter instantiation rejected and repaired, Section 5) | `eb4874d39904593bddc3f2ae638d1bad6756a19d4ebf6e06a60427b6bb7ac1ce` (module `b265e726d110098a09c00de331c4433454973ba741966b3018b77afcc225b5c3`, internal `3ff1e797cd5dbcc4273d9f002b9e540b58829d237b2d344f266e0bea46c02bc8`) |
| `uc/verification/results/liu9-h2-mixture-theorem.json` | **the mixture theorem** below, at the exact `(beta*, m*)` | PROVED | `329f7e2d71af8cd78d1a921c72b9ae4d05b71113134eb3341d69932f19ea24b3` (module `3098a1ca30a0f16582df02e1fb7dfd8fa27970c031167125e696e58ac43dfbab`, internal `b108b781224dff601ffcf3f562703d3b620e7da6d89cf6d000b6ccbd06a4773b`) |

Natural log throughout; `h(u) = -u log u - (1-u) log(1-u)`, `pi(x,y) = xy(1+(1-x)(1-y))`, `K = h o pi`.

**Mixture theorem (PROVED; `liu9-h2-mixture-theorem.json`).** For every Borel probability law `mu`
on `[0,1]` and every conditionally i.i.d. coupling `Gamma = int nu_u (x) nu_u dP_U(u)` with barycentre
`mu = int nu_u dP_U`, writing `M = <mu,x>`,

    (1-beta) <mu (x) mu, h(xy)> + beta int <nu_u (x) nu_u, K> dP_U(u)  >=  (M/m) <mu, h>.

Proof (three identities in the free algebra of pairings, sympy residual 0 at `K = 2,3,4`
components, Arb residual `7e-119` at 60 random mixtures with up to 6 components; the general case is
the same computation with `int dP_U` by Fubini, every kernel being continuous on the compact square):

    (I)   beta sum_k w_k K_kk - beta <mu(x)mu,K> = (beta/2) sum_kl w_k w_l <(nu_k - nu_l)^(x)2, K>
    (II)  <mu(x)mu, (1-beta)h(xy) + beta K> - <mu,h> = <mu(x)mu, R> + ((M-m)/m) <mu,h>
    (III) <mu(x)mu, R> + channel = sum_kl w_k w_l <nu_k(x)nu_l, A> + sum_k w_k <nu_k(x)nu_k, B> + Var_w(<nu_k, phi>)

with `R = D_m/m = P2 + Q2`, `phi = sqrt(Q2(x,x))`, `A = P2 + phi(x)phi`, `B = Q2 - phi(x)phi`. The
right side of (III) is `>= 0` because `A >= 0`, `B >= 0` pointwise and a variance is nonnegative.
`R >= 0` is not needed. **Sharp along a whole family:** for `mu_w = (1-w) delta_0 + w delta_{x*}`,
`0 <= w <= 1` (one component), `numerator = w^2 h(x*^2)` (the cross terms vanish, `K(x*,x*) = h(x*^2)`),
`<mu_w,h> = w h(x*)`, `M = w x*`, and `m* = x* h(x*)/h(x*^2)` give `numerator = (M/m*) <mu_w,h>` EXACTLY
for every `w` - equality at every mean, not only at Liu's minimiser `w = m*/x*` (Arb at the exact
balls: the residual enclosures contain 0 at `w = 3/10, 1/2, 9/10, 1`). So the constant `C = M/m*` of the
theorem is attained, and `c' = 1 - m*` is the best constant this two-protocol inequality at `beta*` can
yield; the laws on `{0,1}` give `numerator = <mu,h> = 0`. **Perturbing the constants cannot
be assumed harmless**: with the binding's 100-digit rationals declared exact, the bound is certified
FALSE at Liu's minimiser (by `-5.95e-102`), because the interior zero of `R` at `(x*,x*)` is exact and
nondegenerate (lowering `m` or moving `beta` breaks `R >= 0` to first order; raising `m` alone does not). Every theorem-validity sign claim in both artifacts is therefore evaluated at
Arb balls containing the exact constants; the rationals are test points for parameter-free identities
and for the deliberate truncation counterexample only.

## 2. From a union-closed family to the constant (Liu's Proposition 3, re-derived)

Everything in this section is elementary and HUMAN-AUDITED (re-derived line by line from Liu
Sections I-II; nothing is taken from the paper without re-proof).

**Setting.** `F` union-closed, `|F| >= 2` (if `F = {A}` with `A` nonempty, every element of `A` has
frequency 1). `X^n` uniform on `F` as 0/1 vectors, so `H(X^n) = log|F| > 0`. Fix `c < c'` and
suppose `P(X_i = 1) < c` for every `i`.

**Protocols (Liu Definition 1).** A protocol assigns to each `(s,t) in [0,1]^2` a law `Pi_{s,t}` on
`{0,1}^2` with marginal means `s` and `t`. Given the law of `X^n` and a protocol, generate
`(X~^n, Y~^n)` sequentially: at step `i`, with `s = P(X_i = 1 | X^{i-1} = X~^{i-1})` and
`t = P(X_i = 1 | X^{i-1} = Y~^{i-1})`, draw `(X~_i, Y~_i) ~ Pi_{s,t}` with fresh randomness. Then
`P(X~_i = 1 | X~^{i-1}, Y~^{i-1}) = s`, so by induction `X~^n ~ X^n` and likewise `Y~^n ~ X^n`
(only pasts of positive probability are ever visited).

**Entropy accounting.** With `(X~ v Y~)` the coordinatewise OR and `S_i := P(X_i = 1 | X^{i-1} = X~^{i-1})`,
`T_i := P(X_i = 1 | X^{i-1} = Y~^{i-1})`:

    H(X~^n v Y~^n) = sum_i H((X~ v Y~)_i | (X~ v Y~)^{i-1})                       (chain rule)
                   >= sum_i H((X~ v Y~)_i | X~^{i-1}, Y~^{i-1})                    (conditioning reduces entropy)
                    = sum_i E[ h( Pi_{S_i,T_i}(0,0) ) ]                             (binary OR, h(u) = h(1-u))
    H(X^n)         = sum_i H(X_i | X^{i-1}) = sum_i E[ h(S_i) ].

`X~^n v Y~^n` is a union of two members of `F`, hence takes values in `F`, hence
`H(X~^n v Y~^n) <= log|F| = H(X^n)` (the uniform law maximises entropy on a support).

**The two protocols.**
* `Pi^(1)` (Gilmer): `Pi_{s,t} = Ber(s) (x) Ber(t)`, `Pi_{s,t}(0,0) = (1-s)(1-t)`. Then `X~^n` and
  `Y~^n` are independent, so `(S_i, T_i) ~ mu_i (x) mu_i` where `mu_i` is the law of
  `P(X_i = 1 | X^{i-1})`.
* `Pi^(2)` (Liu Example 5 with `f(x) = x(1-x)`): common randomness `U` uniform on `[0,1]` and local
  randomness; `X = 0` with probability `(1-s) + s(1-s) sigma(U)`, `sigma = +1` on `U > 1/2`,
  `-1` on `U <= 1/2`, and `Y` likewise from `t` with the SAME `U`. These are probabilities
  (`(1-s)(1+s) = 1 - s^2 <= 1`, `(1-s)^2 >= 0`), the mean of `X` is `s` (`E sigma = 0`), and
  `Pi_{s,t}(0,0) = E[((1-s) + s(1-s)sigma)((1-t) + t(1-t)sigma)] = (1-s)(1-t) + s(1-s)t(1-t)
  = (1-s)(1-t)(1 + st)`, which in `x = 1-s`, `y = 1-t` is `pi(x,y)`. Each `X~_j` is a function of
  `(X~^{j-1}, U_j, local X-randomness)` and each `Y~_j` of `(Y~^{j-1}, U_j, local Y-randomness)` by the
  same rule, so given `U^{i-1}` the pasts `X~^{i-1}`, `Y~^{i-1}` are i.i.d.; hence `(S_i, T_i)` is
  conditionally i.i.d. given `U^{i-1}`: its law is `Gamma_i = int nu_u (x) nu_u dP_{U^{i-1}}(u)` with
  barycentre `mu_i` (Liu, paragraph after Example 5). The fresh `U_i` makes
  `P((X~ v Y~)_i = 0 | X~^{i-1}, Y~^{i-1}) = Pi_{S_i,T_i}(0,0)` exactly.

**Convex combination.** With `beta = beta*` in `(0,1)`:

    (1-beta) H(X~^(1)n v Y~^(1)n) + beta H(X~^(2)n v Y~^(2)n)
        >= sum_i [ (1-beta) <mu_i (x) mu_i, h((1-s)(1-t))> + beta E_{Gamma_i} h((1-s)(1-t)(1+st)) ]
         = sum_i numerator_i                                   (in x = 1-s variables: h(xy), h(pi))
        >= sum_i (M_i/m) <mu_i, h>                             (mixture theorem)
         > ((1-c)/m) sum_i E[h(S_i)] = ((1-c)/m) H(X^n),

because `M_i = 1 - P(X_i = 1) > 1 - c` (and `h(x) = h(1-x)`, so `<mu_i,h>` in either variable is
`E[h(S_i)]`; coordinates with `E[h(S_i)] = 0` contribute `0 >= 0` and the strict inequality survives
because `H(X^n) > 0` forces some coordinate with `E[h(S_i)] > 0`). The left side is at most
`(1-beta) H(X^n) + beta H(X^n) = H(X^n)`. So `H(X^n) > ((1-c)/m) H(X^n)` with `H(X^n) > 0`, i.e.
`(1-c)/m < 1`, i.e. `c > 1 - m = c'`, contradicting `c < c'`. Hence for every `c < c'` some `i` has
`P(X_i = 1) >= c`; the frequencies are finitely many fixed numbers, so `max_i P(X_i = 1) >= c'`. QED.

**Why the strictness is free.** Liu certifies `c'` through "the nine-parameter optimum is `>= 1`",
which gives `H(X~ v Y~) >= H(X^n)` and needs a separate argument for `C > 1` (his Proposition 2/3
phrasing "for any `c < c*` there exists `C > 1`"). The quantitative form `numerator >= (M/m) ehx`
supplies `C = (1-c)/m = (1-c)/(1-c') > 1` uniformly for every `c < c'`, with no mean constraint and no
limiting argument.

**What is no longer needed.** Liu's Section V-A hypothesis (positive semidefiniteness of the
protocol kernel on a codimension-2 subspace; Lemma 11 at `l = 1`) enters only Theorem 12, the
Krein-Milman reduction of the conditionally i.i.d. couplings to three shared-mass atoms. The mixture
theorem handles every coupling directly, so the reduction - and with it the repo's own H1 theorem
(`LIU_H1/`) - is not load-bearing for `c'`. Liu's Section V-B hypothesis (the minimiser structure) is
replaced by the theorem `gap >= 0` for `M >= m` (equality at his minimiser). Theorem 9 (two-mixture
sufficiency) and the weak-closure `C_3(mu)` are likewise not used: the actual induced coupling is a
mixture, and the inequality is proved for it.

**Base of logarithms.** Liu's `h` is in bits; the inequality `numerator >= (M/m) ehx` is homogeneous
in the base, and the entropy accounting above holds in any base.

## 3. The value of c'

`m` is the binding's `mean = p* x*`; `1 - c - m = 0` exactly as stored rationals, and
`certify_equation_parameters` gives

    c' in [0.382709087918735029930312902098972611626381433 +/- 4.95e-46].

So `c' > 0.38270908791873` (MACHINE-VERIFIED). **Liu's printed (93), `c' ~ 0.382709087918741`, is wrong
in its last two digits** (`...8741` vs the true `...873503`); the ledger previously repeated that
decimal and must cite the certified value instead. `c'` exceeds `c* = 0.3823455` (Liu's (11), the
Sawin-Yu-Cambie constant) and this repository's certified `psi + 3e-4 = 0.3822660112501052`
(Campaign K). No universal priority claim is made: this ledger's literature search is bounded.

## 4. Scope, labels, and what an external referee should check

* PROVED with MACHINE-VERIFIED algebra: the mixture theorem (Section 1) and the three inequalities
  `A >= 0`, `B >= 0`, `Phi >= 0` it rests on.
* HUMAN-AUDITED: Section 2 (the protocol construction, the entropy accounting, the conditionally
  i.i.d. structure of `(S_i, T_i)`, the union-closed entropy bound). Each step is standard
  (Gilmer 2022, Sawin 2022, Liu 2023 Proposition 3) and was re-derived here; an independent
  adversarial read is recorded below.
* CITED-DEPENDENCY: none load-bearing. Liu's Theorem 9, Lemma 11, Theorem 12, and Section V-A are
  not used.
* Constants: `beta*` and `m*` are the root-defined solutions of Liu's (87)-(90). The certificates
  `A >= 0`, `R >= 0` hold at those exact constants - the tight point `(x*,x*)` is handled by exact
  algebra in `liu9-h2-twovar.json` (`R(x*,x*) = 0` from `h(pi(x*,x*)) = h(x*^2)` and `m* = x* h(x*)/h(x*^2)`;
  `grad R = 0` from the definition of `beta*`) - and NOT for perturbed values in the balls. The balls
  are how the exact constants are represented in Arb, never a licence to perturb them.
* The result is for the constant only; it says nothing about the conjecture's `1/2`.

## 5. Independent review

Two blank-context adversarial reviewers (code-reviewer agents, read-only), by different routes.

**`LiftSkeptic` on `liu9_h2_general_lift.py` - first verdict UNSOUND, final verdict SOUND.** It
confirmed the functional, the mean constraint, every step of the composition L1-L4, the absence of
any density step, the product-coupling identity, and the mutations, and then found the one real
defect of the first delivery: the module instantiated the theorem at the binding's 100-digit
rationals and declared them exact, while `R >= 0` is a theorem at the root-defined constants only.
Its 900-bit Arb attack (14,233 configurations) produced the certified negative
`[-5.9548013684e-102 +/- 4.85e-222]` at Liu's minimiser with those rationals, and the module's own
W1 had accepted `|gap| <= 1e-90` regardless of sign. Repairs: the constants are now defined
mathematically and propagated as Arb balls from a re-certified 70-digit bracket of `x*`, with
cross-checks against the binding and against `twovar`'s definitions; every theorem-validity sign
claim is evaluated at the balls; the rationals are serialised as exact fractions and used only for
parameter-free identities; W1a is PROVED by algebra with the ball enclosure as consistency check;
W1b records the counterexample as a certified negative and a seventh mutation; `liu9_binding.py` and
`liu9_objective.py` are pinned. A second-round overclaim ("false for every perturbation") was also
caught - raising `m` keeps `R, A, B >= 0` by L2 - and corrected. Final: SOUND, hashes recomputed by
the reviewer.

**`ConstantSkeptic` on `liu9_h2_mixture_theorem.py` and this note - SOUND, no BLOCKER/MAJOR/MINOR.**
It re-derived (I)-(III) for general `K` by hand, ran the exact sympy check at `K = 1, 5, 8` as well,
wrote out the continuum (Markov-kernel) case with the measurability facts needed and found no gap,
confirmed `numerator_ehx` against `liu9_objective._formula`, verified every step (a)-(h) of the
Proposition 3 chain including the fresh-`U_i` point and the `|F| <= 1` cases, re-derived Example 5
and the mean direction, recomputed the constants independently at 170 digits (`c' = 1 - m* =
0.382709087918735029930312902098972611626381432505757587...`, Liu's printed (93) exceeding it by
`5.97e-15`), and ran its own attack: 27,006 random mixtures with up to 8 components and 6 atoms per
component, atoms at 0, 1 and `10^-k`, weights up to `9e18 : 1`, **no mean filter**, all recomputed at
170 digits - no certified negative; smallest nonzero-entropy margin certified positive at
`3.57e-53`. Four nits (two X1 mutants sharing one error polynomial; a provenance sentence inherited
from the lift; "polynomial" vs "rational-function" identities; the `setdefault` thread environment)
are applied in the final versions.

Both final verdicts are for the hashes in the table of Section 1.
