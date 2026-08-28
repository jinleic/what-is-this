# Gate B — closure-defect scalar route

## Current status

**Start at [`CURRENT_STATUS.md`](CURRENT_STATUS.md)** for the verified current
facts and the next frontier alone. This file is the full evidence map and
reproduction index.

**RESOLVED in the repository's stated scope (2026-08-25):**

\[
c_{\rm cl}^\star=+\infty.
\]

The certified normalized \(n=7\) family \(\mathcal B\) has
\(A_+(\mathcal B)<-7/250\) and join-success probability \(17/81\). The
Cartesian powers \(\mathcal F_k=\mathcal B^{\boxtimes k}\) remain normalized,
cap-\(2/5\), and Reimer-admissible. The objective tensorizes while closure
success multiplies:

\[
A_+(\mathcal F_k)=kA_+(\mathcal B),\qquad
\varepsilon_\vee(\mathcal F_k)=1-(17/81)^k.
\]

Therefore

\[
-A_+(\mathcal F_k)/\varepsilon_\vee(\mathcal F_k)>7k/250\to\infty.
\]

The conclusion now has two arithmetically independent finite-base
certificates. The sharper 256-bit Arb computation proves
\(A_+(\mathcal B)<-7/250\). A pure-Python 160-bit dyadic checker instead
relaxes every nondegenerate Bellman action interval to \([0,1]\), avoids every
clamp decision, and still proves
\(A_+(\mathcal B)<-1/40\). Either negative bound is sufficient for the
Cartesian-power divergence.

The conclusion now has a second, combinatorially independent base.  The
25-row \(3+3\) weight-cell family
\[
\mathcal D=\{(0,0),(0,1),(1,0),(1,2),(2,1)\}
\]
has exact closure success \(181/625\).  A third-party-free dyadic checker
evaluates all 720 coordinate orders without a symmetry quotient and proves
\(A_+(\mathcal D)<-1/80\).  Hence
\[
\frac{-A_+(\mathcal D^{\boxtimes k})}
{\varepsilon_\vee(\mathcal D^{\boxtimes k})}
>\frac{k/80}{1-(181/625)^k}>\frac{k}{80}\to\infty.
\]
The same \(n=6\) base also has a 256-bit Arb certificate
\(A_+(\mathcal D)<-17/1250\), produced by the generic standalone Bellman
evaluator.  The two bases are therefore each certified by both arithmetic
routes; the all-720-order dyadic route remains the one that avoids symmetry and
clamp classification simultaneously.  Thus an error specific to the 45-row
family or its order-orbit reduction cannot reverse the Gate B conclusion.

The negativity of both bases now also has an **exact rational** certificate
that uses no interval library and no relaxation of the action set. It evaluates
the true one-sided Bellman recursion over the exact feasible interval
\([s^*,U]\) with Python fractions and returns a two-sided enclosure:
\[
\begin{aligned}
A_+(\mathcal D)&\in[-0.01367210773217777356185991566,\,
 -0.01367210773217777356185991299],\\
A_+(\mathcal B)&\in[-0.02864918679446831781602698193,\,
 -0.02864918679446831781602697645],
\end{aligned}
\]
each of width below \(6\times10^{-27}\), computed over **all** \(720\) and
**all** \(5040\) coordinate orders with no automorphism quotient. Both Arb
values lie strictly inside these enclosures. The three certified scalar
primitives are a bit-by-bit binary logarithm, repeated ceiling square roots of
two, and a concavity case split for the constrained entropy maximum.

Because the enclosure is two-sided, an inverted clamp comparison, a mis-stated
feasible interval, or a sign error breaks the sandwich instead of silently
moving one bound. The upper endpoints alone imply the theorem.

The divergence is also now quantified. With
\(c_{\rm cl}^\star(n)\) the supremum over admissible families on at most \(n\)
coordinates, the powers give \(c_{\rm cl}^\star(n)>\frac1{250}n-\frac7{250}\),
while \(Q\ge0\), \(C_+\ge0\) give \(-A_+\le\log_2m\le n\) and hence
\(c_{\rm cl}^\star(n)\le n/\varepsilon_0\) on families with
\(\varepsilon_\vee\ge\varepsilon_0\). So under a positive defect floor the
growth is exactly of order \(n\), and the only escape is
\(\varepsilon_\vee\to0\).

No fixed scalar closure-defect coefficient can repair \(A_+\) on all
cap/Reimer families. This is not a union-closed counterexample and does not
resolve Frankl's conjecture.

## The local regime (2026-08-27)

The defects in the two published constructions tend to one, and that is
structural, not incidental: the closure defect *multiplies through success*
under products, so every Cartesian power has
\(\varepsilon_\vee\ge\varepsilon_\vee(\text{base})\) and the mechanism that
proves \(c_{\rm cl}^\star=+\infty\) cannot exhibit a single family with
\(\varepsilon_\vee<444/625\). The surviving question is therefore

\[
c_{\rm loc}=\lim_{\varepsilon\to0^+}
\sup\Big\{\tfrac{-A_+(\mathcal F)}{\varepsilon_\vee(\mathcal F)}:
\mathcal F\text{ cap/Reimer},\ A_+(\mathcal F)<0,\
0<\varepsilon_\vee(\mathcal F)\le\varepsilon\Big\}.
\]

**A third certified base lowers the frontier.** A search witness with 45 rows
on seven coordinates, all degrees at the cap 18 and incidence 126, has exact
closure defect \(1336/2025=0.6598\ldots\) and a certified negative enclosure

\[
A_+\in[-0.00085518533372314171100789579,\,
 -0.00085518533372314171100789029],
\]

of width \(5.5\times10^{-27}\), so \(A_+<-1/1200\). Its automorphism group is
**trivial**, so all 5040 orders carry distinct enclosures and no symmetry
quotient is available even in principle; it is also not a weight-cell family.
Its powers give \(\varepsilon_\vee=1-(689/2025)^k\) and ratio \(>k/1200\), a
third infinite construction. Certificate:
[`certificates/gate_b_lowdefect_rational_v1.json`](certificates/gate_b_lowdefect_rational_v1.json).
Hence, writing
\(e^\star=\min\{\varepsilon_\vee(\mathcal F):\mathcal F\text{ admissible},\
A_+(\mathcal F)<0\}\),

\[
e^\star\ \le\ \frac{1336}{2025}\ =\ 0.6597\ldots
\qquad\text{(certified, down from }64/81=0.7901\ldots\text{)}.
\]

**The measured trend runs against a local counterexample.** Minimizing
\(A_+\) under a hard defect cap at \(n=7\), \(m=45\) gives a monotone curve:

| defect cap | least \(A_+\) found | negative |
|---|---|---|
| 0.79 | \(-0.01628\) | yes |
| 0.70 | \(-0.01490\) | yes |
| 0.66 | \(-0.00086\) | yes |
| 0.64 | \(+0.01171\) | no |
| 0.52 | \(+0.04399\) | no |
| 0.40 | \(+0.06940\) | no |
| 0.35 | \(+0.07951\) | no |

As the defect falls, \(-A_+\) falls faster, so the ratio *collapses* rather
than diverging: it is \(0.0363\) at the published base and \(0.0013\) at the new
one. Independently, the 38 minimum-defect witnesses over \(n\in\{6,7\}\),
\(m\ge15\) all have \(A_+>0\), from \(+0.0588\) to \(+0.5607\), and the searched
minimum defect of an admissible family never fell below \(2/9\) at any
\(n\le8\) although the combinatorial floor \(2/m^2\) is three orders of
magnitude smaller.

**Why this is not a proof, precisely.** A family with \(\varepsilon_\vee=0\) is
union-closed with all degrees \(\le\lfloor2m/5\rfloor<m/2\), i.e. a Frankl
counterexample. So proving a universal positive defect floor at cap \(2/5\) --
which is what would empty the local regime -- implies every union-closed family
has an element in more than \(2m/5\) of its sets, improving the published
frontier \(\psi=0.381966\ldots\) by about \(0.018\). That is the blocker, and it
is recorded in [PROOF.md](PROOF.md) with Lemmas 4-6 and Proposition 7 rather
than presented as evidence of impossibility. The local question stays **OPEN**;
what changed is that it is now measured, bracketed, and certified at its upper
end.

The ratio uses the extended-nonnegative convention at zero defect: a family
with \(A_+<0\) and \(\varepsilon_\vee=0\) has value \(+\infty\). The
constructed families all have \(0<\varepsilon_\vee<1\), so the proof does not
use this convention.

## Evidence map

| Artifact | Status | Role |
|---|---|---|
| [DEFINITIONS.md](DEFINITIONS.md) | exact statement | Freezes \(A_+\), \(\varepsilon_\vee\), cap, Reimer, normalization, and ordered-pair conventions. |
| [PROOF.md](PROOF.md) | proved | Complete Cartesian-product construction and tensorization proof. |
| [candidates/n7_block_extremizer.json](candidates/n7_block_extremizer.json) | exact data | Machine-readable certified base family. |
| [verify_gate_b.py](verify_gate_b.py) | standalone checker | No imports from other `math/uc` modules; exact filters, 256-bit Arb base certificate, direct square combinatorics, and actual product Bellman checks for consecutive and alternating orders. |
| [verify_gate_b_dyadic.py](verify_gate_b_dyadic.py) | independent standalone checker | Standard-library-only 160-bit dyadic intervals with rational transcendental remainders. It uses a clamp-free Bellman relaxation and proves \(A_+(\mathcal B)<-1/40\), independently sufficient for unboundedness. |
| [certificates/gate_b_unbounded_dyadic_v1.json](certificates/gate_b_unbounded_dyadic_v1.json) | independently certified | Exact dyadic endpoints for \(Q\), the relaxed Bellman upper bound, entropy, and \(A_+\), plus the exact family, defect, symmetry partition, and all-power formulas. |
| [verify_gate_b_n6_dyadic.py](verify_gate_b_n6_dyadic.py) | independently certified second base | Reconstructs the separate 25-row \(n=6\) family, evaluates the clamp-free relaxation on all 720 orders without a symmetry quotient, and proves \(A_+<-1/80\). It shares only the audited dyadic interval primitive with the \(n=7\) checker. |
| [certificates/gate_b_unbounded_n6_dyadic_v1.json](certificates/gate_b_unbounded_n6_dyadic_v1.json) | independently certified | Exact \(n=6\) rows, filters, defect \(444/625\), all-order count, interval bounds, and the second all-power consequence. |
| [verify_gate_b_n6_arb.py](verify_gate_b_n6_arb.py) | independently certified second base | Reconstructs the 25-row family separately and evaluates its 10 exact order orbits with the generic 256-bit Arb Bellman checker, proving \(A_+<-17/1250\). |
| [certificates/gate_b_unbounded_n6_arb_v1.json](certificates/gate_b_unbounded_n6_arb_v1.json) | certified | Exact \(n=6\) family, 72-element automorphism group, Arb intervals, rational bound, and all-power consequence. |
| [verify_gate_b_rational.py](verify_gate_b_rational.py) | exactly certified, both bases | Standard-library-only exact rational evaluator. No interval class, no third-party arithmetic, and no action-set relaxation: it maximizes over the true feasible interval \([s^*,U]\) and returns a two-sided enclosure of \(A_+\). Resumable per coordinate order with a schema-checked append-only checkpoint. |
| [certificates/gate_b_unbounded_rational_v1.json](certificates/gate_b_unbounded_rational_v1.json) | exactly certified | All 720 \(n=6\) and all 5040 \(n=7\) orders, exact rational enclosures of \(Q\), \(\mathbb E_\pi C_{+,\pi}\), \(\log_2m\) and \(A_+\) for both bases, the frozen rational targets each enclosure passes, and the linear-growth block. |
| [experiments/rational_certificate_checkpoint.jsonl](experiments/rational_certificate_checkpoint.jsonl) | append-only log | 5,760 exact per-order enclosures. A complete resume reports `new_evaluations: 0`; foreign-schema and truncated lines are counted and ignored, never trusted. |
| [certificates/gate_b_lowdefect_rational_v1.json](certificates/gate_b_lowdefect_rational_v1.json) | exactly certified third base | The low-defect search witness: 45 rows on seven coordinates, all degrees at the cap, closure defect \(1336/2025\), trivial automorphism group, all 5040 orders distinct, \(A_+<-1/1200\) with a \(5.5\times10^{-27}\) enclosure. Certifies \(e^\star\le1336/2025\). |
| [search_local_defect.py](search_local_defect.py) | exact for \(n\le4\), searched above | Minimizes the closure defect alone, with no entropy evaluation. `--mode exact` enumerates every admissible family for \(n\le4\); `--mode search` anneals with degree-preserving, free-replacement and closure-repair proposals; `--mode frontier` reports the best re-verified witness per size across every record schema. |
| [experiments/local_defect_checkpoint.jsonl](experiments/local_defect_checkpoint.jsonl) | append-only log | One record per completed restart. Schema strictness governs work skipping only: a witness family from an older algorithm is still valid evidence once re-verified, and the schema of origin is reported. |
| [candidates/local_defect_frontier_n7.json](candidates/local_defect_frontier_n7.json) | searched upper bounds | Least closure defect found per family size at \(n=7\), each witness re-checked for cap, integer incidence, distinct rows and activity. The minimum never approaches the \(2/m^2\) pair floor. |
| [hunt_local_ratio.py](hunt_local_ratio.py) | searched | Three instruments for the local regime: `hunt` minimizes \(A_++\lambda\varepsilon_\vee\) from random starts, `descend` minimizes the defect from a certified negative family under a hard negativity constraint, and `sweep` minimizes \(A_+\) under a hard defect cap, seeded from the best family already known under that cap. |
| [candidates/local_sweep_n7.json](candidates/local_sweep_n7.json) | searched frontier curve | Least \(A_+\) against the defect cap at \(m=45\), every point recomputed over all 5040 orders. Monotone; crosses zero between caps \(0.64\) and \(0.66\). |
| [candidates/local_descend_n7.json](candidates/local_descend_n7.json) | searched | 64 negative families found by walking down from \(\mathcal B\); lowest defect \(1336/2025\), later certified exactly. |
| [audit_tensorization_exact.py](audit_tensorization_exact.py) | exact falsification search | Repeats the product-lemma attack in exact rational arithmetic. For every global order of every product it demands that the certified enclosures of the two sides intersect; non-intersection would be a counterexample, not rounding. Resumable at case granularity. |
| [candidates/tensorization_exact_audit.json](candidates/tensorization_exact_audit.json) | exact adversarial | Zero intersection failures, zero exact defect-identity failures, and a worst fixed-order Bellman gap of exactly zero — replacing the float64 \(1.8\times10^{-15}\) evidence for Lemma 3. |
| [experiments/tensorization_exact_checkpoint.jsonl](experiments/tensorization_exact_checkpoint.jsonl) | append-only log | One record per product case; a complete resume reports `new_evaluations: 0`. |
| [audit_n6_square_direct.py](audit_n6_square_direct.py) | resumable falsification search | Materializes the second base's 625-row square as plain 12-bit masks and compares generic Bellman/iid evaluations against the induced block sums for five hostile global orders. |
| [candidates/n6_square_direct_audit.json](candidates/n6_square_direct_audit.json) | direct numerical | Exact square constraints and defect; zero Bellman product gap and \(1.8\times10^{-15}\) worst iid gap across five orders. |
| [experiments/n6_square_direct_checkpoint.jsonl](experiments/n6_square_direct_checkpoint.jsonl) | append-only log | Five deterministic order records; a complete resume reports zero new evaluations and leaves the artifact byte-identical. |
| [certificates/gate_b_unbounded_arb_v4.json](certificates/gate_b_unbounded_arb_v4.json) | certified | Authoritative interval output. It verifies the full cell-union row definition, directly certifies the base ratio, separates actual square-order checks from identity-derived values, and records exact power ratios. The v1--v3 JSON files remain intermediate checkpoints. |
| [search_n8_block_symmetric.py](search_n8_block_symmetric.py) | complete finite search | Checkpointed exhaustive \(S_k\times S_{8-k}\) class on \(n=8\); \(k\in\{1,2,3,4\}\) exhausts every block-symmetric class. Exact class/defects; **float64 objective, discovery only**. |
| [candidates/n8_k1_census.json](candidates/n8_k1_census.json) | complete numerical | 65,535 raw masks, 136 admissible canonical families, two negative; maximum ratio \(0.03433549\ldots\), below the certified \(n=7\) base. |
| [experiments/n8_k1_checkpoint.jsonl](experiments/n8_k1_checkpoint.jsonl) | append-only log | All 136 evaluated records; resumable without overwriting prior records. |
| [candidates/n8_k2_census.json](candidates/n8_k2_census.json) | complete numerical | The \(S_2\times S_6\) class: 21 cells, 2,097,151 raw masks, 2,272 canonical admissible families, 21 with negative float64 objective. Its ordering is *not* reliable — see the exact re-ranking below. |
| [experiments/n8_k2_checkpoint.jsonl](experiments/n8_k2_checkpoint.jsonl) | append-only log | One record per canonical family; resumable. |
| [certify_class_negatives.py](certify_class_negatives.py) | exact re-ranking | Re-evaluates census candidates in exact rational arithmetic with the same primitives as the rational verifier, and records a two-sided enclosure per family. Reports separation rather than requiring it, since `DEFINITIONS.md` makes normalization a reporting condition. Resumable per family. |
| [candidates/n8_k2_exact_frontier.json](candidates/n8_k2_exact_frontier.json) | exact | All 36 candidates with float \(A_+<0.01\) evaluated exactly: 21 float-negative, **20 certified negative**, one float negative exactly *positive* (\(2553/3200\)), and **no** float-non-negative family exactly negative. Largest float/exact gap \(2.03\times10^{-3}\). Certifies \(e^\star\le139/245\) and ratio \(\ge0.0580391\). |
| [experiments/n8_k2_exact_checkpoint.jsonl](experiments/n8_k2_exact_checkpoint.jsonl) | append-only log | One exact record per family; a complete resume reports `new_evaluations: 0`. |
| [certificates/gate_b_n8_rational_v1.json](certificates/gate_b_n8_rational_v1.json) | exactly certified, two bases | `n8lo`: 75 rows on eight coordinates, defect \(1144/1875\), \(A_+<-1/500\) — moves the local frontier below the \(n=7\) value. `n8hi`: defect \(468/625\), \(A_+<-693/25000\), certified ratio \(\ge77/2080=0.0370192\ldots\) — the first here to beat the published \(0.0362591\ldots\). Both use the **exact symmetry quotient** (\(|\mathrm{Aut}|=1440\), 28 orbits of size 1440); the quotient's two hypotheses are tested, not assumed. |
| [bound_local_regime.py](bound_local_regime.py) | proved, standalone | Independent stdlib-only checker for the size ceiling \(\log_2m\le4n/5\) (integer test \(m^5\le2^{4n}\)), the union-growth and downset-capacity defect floors, and the maximality of the certified bases. Imports nothing from the search or objective modules. |
| [candidates/local_regime_bounds.json](candidates/local_regime_bounds.json) | proved + exhaustive audit | Verdict `ALL_BOUNDS_HOLD`: ceiling holds for every admissible size at every \(3\le n\le16\); zero violations over all 366 admissible families with \(n\le4\); the capacity floor is attained with slack exactly \(0\); zero admissible additions for all three \(n\le7\) bases. |
| [audit_n5_complete.py](audit_n5_complete.py) | complete finite audit | Resumable independent float64 evaluation of all 5,172 exact nontrivial \(n=5\) cap/Reimer coordinate orbits representing 463,343 labeled families. |
| [candidates/n5_complete_audit.json](candidates/n5_complete_audit.json) | complete numerical | Zero negative or near-zero values; minimum \(A_+=0.006474313148642441\ldots\). |
| [experiments/n5_complete_checkpoint.jsonl](experiments/n5_complete_checkpoint.jsonl) | append-only log | All 5,172 evaluated nontrivial \(n=5\) records; exact filters/defects and float64 objective values. |
| [audit_tensorization.py](audit_tensorization.py) | resumable falsification search | 292 ordered product cases (205 with unequal block dimensions), every global order of each product, attacking fixed-order \(Q\) and \(C_+\) additivity and the exact defect product. |
| [candidates/n8_k3_census.json](candidates/n8_k3_census.json) | complete numerical | The \(S_3\times S_5\) class: 24 cells, 16,777,215 raw masks, 13,470 canonical admissible families, 54 with negative float64 objective. Ordering not reliable; see the exact re-ranking. |
| [experiments/n8_k3_checkpoint.jsonl](experiments/n8_k3_checkpoint.jsonl) | append-only log | One record per canonical family; resumable. |
| [candidates/n8_k3_exact_frontier.json](candidates/n8_k3_exact_frontier.json) | exact | All 110 candidates with float \(A_+<0.01\) evaluated exactly: 54 float-negative, **52 certified negative**, three float negatives exactly non-negative, and **one family the census called non-negative is exactly negative** (defect \(94/125\)) — the float evaluator fails in both directions. Largest gap \(3.24\times10^{-3}\). |
| [experiments/n8_k3_exact_checkpoint.jsonl](experiments/n8_k3_exact_checkpoint.jsonl) | append-only log | One exact record per family; a complete resume reports `new_evaluations: 0`. |
| [certificates/gate_b_n8_max_rational_v1.json](certificates/gate_b_n8_max_rational_v1.json) | exactly certified | `n8max`: 80 rows at the **maximum** admissible size for \(n=8\), every degree at the cap 32, incidence \(256=8\cdot32\) — the largest any \(n=8\) family can carry. Defect \(2343/3200\), \(A_+<-19/625\), certified ratio \(\ge2432/58575=0.0415194\), the best among separating families. |
| [certificates/gate_b_n8_clone_rational_v1.json](certificates/gate_b_n8_clone_rational_v1.json) | exactly certified, two bases | The cloned-coordinate witnesses, admissible but not separating (coordinates 0 and 1 coincide). `n8clone_lo`: defect \(139/245\), the lowest with \(A_+<0\) anywhere. `n8clone_hi`: defect \(317/490\), \(A_+<-3/80\), certified ratio \(\ge147/2536\), and **asymptotic slope \(3/640=0.0046875\) against the published \(1/250\)** — the certified growth constant, improved by \(75/64\). |
| [audit_clone_power.py](audit_clone_power.py) | proved + instantiated | Admissibility of every power of a cloned base, two ways: the exact integer argument (cap exactly met for all \(k\); Reimer reduces to the single base-level inequality \(m^m\le2^{2I}\), i.e. \(I\ge R_m\)), and direct instantiation of the 4,900-row square with no product-aware shortcut. |
| [candidates/clone_power_audit.json](candidates/clone_power_audit.json) | proved | Verdict `EVERY_POWER_ADMISSIBLE`: degrees equal the cap at every power to 40, Reimer strict throughout, and the instantiated square's defect is \(210171/240100=1-(173/490)^2\) exactly. |
| [candidates/n8_k4_census.json](candidates/n8_k4_census.json) | complete numerical | The \(S_4\times S_4\) class: 25 cells, 33,554,431 raw masks, 11,553 canonical admissible families, 73 screen-negative. Closes block-symmetric coverage at \(n=8\). |
| [experiments/n8_k4_checkpoint.jsonl](experiments/n8_k4_checkpoint.jsonl) | append-only log | One record per canonical family; resumable. |
| [candidates/n8_k4_exact_frontier.json](candidates/n8_k4_exact_frontier.json) | exact | All 125 candidates with screen value \(<0.01\) evaluated exactly: 73 screen-negative, **73 certified negative**, zero screen errors in either direction. |
| [experiments/n8_k4_exact_checkpoint.jsonl](experiments/n8_k4_exact_checkpoint.jsonl) | append-only log | One exact record per family; a complete resume reports zero new evaluations. |
| [certificates/gate_b_n8_k4_rational_v1.json](certificates/gate_b_n8_k4_rational_v1.json) | exactly certified, two bases | `n8tiny`: fifteen rows, defect \(4/9=0.4444\), \(A_+<-1/2100\) — **the lowest defect with \(A_+<0\) anywhere**, and it needed only a dominant set, not a large family. `n8best`: 75 rows, separating, defect \(284/375\), \(A_+<-177/5000\), slope \(177/40000=0.004425\) against the published \(1/250\) — the growth constant improved **without** leaving the separating class, and ratio \(531/11360\), the best among separating families. |
| [candidates/n8_k1_exact_frontier.json](candidates/n8_k1_exact_frontier.json) | exact | The \(k=1\) class re-ranked exactly: its published ratio \(0.034335490331\) is confirmed to every digit, so the screen is exact for that class. |
| [audit_clone_saturation.py](audit_clone_saturation.py) | exact | The cloning transformation measured exactly. Duplicating a coordinate is a bijection commuting with union, so it preserves \(m\), \(\varepsilon_\vee\), every degree and hence the cap, while raising incidence against a fixed \(R_m\): **admissibility is preserved unconditionally and the defect cannot move.** What moves is \(A_+\), downward. |
| [candidates/clone_saturation_audit.json](candidates/clone_saturation_audit.json) | exact | Verdict `CLONING_IS_DEFECT_FREE_AND_SATURATES`. Collapsing `n8tiny`'s four identical coordinates leaves an admissible **five**-coordinate family of the same size 15 and same defect \(4/9\) whose objective is *positive*, \(+0.010695694\); three clones drive it to \(-0.000478465\). So its clones buy negativity, not feasibility — unlike `n8clone_lo`/`n8clone_hi`, which collapse to inadmissible seven-coordinate families. The gain saturates: deltas \(-7.22\times10^{-3},-2.71\times10^{-3},-1.25\times10^{-3},-6.59\times10^{-4}\), ratios \(0.375,0.462,0.527\). |
| [candidates/tensorization_audit.json](candidates/tensorization_audit.json) | adversarial numerical | Zero counterexamples; worst fixed-order gaps \(1.8\times10^{-15}\) (Bellman) and \(8.9\times10^{-16}\) (iid), i.e. float64 rounding. |
| [audit_square_direct.py](audit_square_direct.py) | resumable direct check | Materializes the 2,025-row square as plain 14-bit masks and runs the generic evaluators, with no factored-fiber or product-aware shortcut. |
| [candidates/square_direct_audit.json](candidates/square_direct_audit.json) | direct numerical | Five declared global orders; exact square combinatorics and zero product gap in every order. |
| [audit_power_admissibility.py](audit_power_admissibility.py) | exact constructed-family audit | Builds \(\mathcal F_k\) for \(k\le3\) and recomputes size, counts, cap, incidence, Reimer, normalization, and the missing-join count by brute force over all ordered pairs. |
| [candidates/power_admissibility_audit.json](candidates/power_admissibility_audit.json) | exact | At \(k=3\): 91,125 rows, all counts 36,450 equal to the cap bound, incidence 765,450 against threshold 750,668, and 8,227,000,000 of 8,303,765,625 ordered joins missing, i.e. exactly \(1-(17/81)^3=526528/531441\). |
| [test_gate_b.py](test_gate_b.py) | regression suite | 42 tests: exact constraints, 5,400 fixed-order comparisons over all 225 products of nonempty \(n=2\) families, both bases under Arb, dyadic and exact rational arithmetic, all-720-order \(n=6\) certification, Arb-versus-rational bracketing of every certified primitive, 65-point sampled domination of the constrained maximum, per-orbit enclosure agreement, sharpness against the dyadic relaxation, state-canonicalization invariance, orbit-versus-all-order certificate equality, foreign-schema and truncated-checkpoint rejection, exact finite-to-general checks for both power sequences, the exact-integer growth corollary, symbolic Bellman algebra, all direct-square and tensorization audits, symmetry/order coverage, and artifact contracts. |
| [EXPERIMENTS.md](EXPERIMENTS.md) | chronological ledger | Commands, parameters, failures, scope labels, and resource use. |
| [paper/main.tex](paper/main.tex) | paper draft | Paper-quality statement, both infinite constructions, six finite-base certificate boundaries, the exact rational proposition, the \(\Theta(n)\) growth section, computational controls, and limitations. |

Authoritative SHA-256 values after final review:

```text
405b0b74129c6192be791788b22ffbfbcb921e407e1d0002d01fb8193d131c72  verify_gate_b.py
7d8a1484cf76a6fb2cf274c3874b52896662e85f005844601d962a40fb7f32a9  verify_gate_b_dyadic.py
c43cc21d6f0da26eca071319fb73c978c532d0917ca1f7cd7e7dfb2ae68878c7  verify_gate_b_n6_dyadic.py
ed6814505eb89f072e4efe8cf50a828aa3b514e024ec17c8933d5dae4a725b60  verify_gate_b_n6_arb.py
0f024845a7054a46afdacf46405955ad5e73829092b0aa2d014ca98d664f5419  verify_gate_b_rational.py
dbe8064eed8d68dc3b49a7823242b599188dfcc82d98e892cf19b332fcc971ea  audit_tensorization.py
1880fdc7752815f7d47054f0b862dd5be6ae3c13df003290c3e5c64dfa702e98  audit_tensorization_exact.py
a2a9477a71153242726311df9792006a1f614cee6b178a60b9c886fd3425740d  audit_square_direct.py
3836ddcabe26e839ec95b1334019e7ce3ea44bb95a05f35b854801164b1701a9  audit_n6_square_direct.py
6f963f4eaf05408a10bbc180854df30aa1052a5771eae3bc1f03000fcd3ef016  audit_power_admissibility.py
a79651c87edd36332fe1d0c4ff3b59b8efaad96f5194121d99a3dbd7e525b3b7  audit_n5_complete.py
fc07a078ee1098f7f3942044ac21765f14477b67af89ef6bd1ed826c4a1db81f  search_n8_block_symmetric.py
6434f78c2188489e2229b3089fde27ab97ae2e0ac7dcfa467684a1889a030fc1  test_gate_b.py
9deabde710eff1e9d49b7c198edecc8d8af4e3d606a9e69feb179eb52df4df36  candidates/n7_block_extremizer.json
14f6875418849011ccc350d5a0d1da7a57eae8bfe695182559184de57d281588  candidates/n5_complete_audit.json
d65d5f5e2735295ae5988cd1077814dc353d8609f7d951d45d150672e25b112d  candidates/n8_k1_census.json
0f48dc396827fad4043f148a552f2b2d1a4bceb7092d6ce4fd3a4b9cef2801da  candidates/tensorization_audit.json
456520716b0c5cf30d252ff96f16a92d03019a56b98585fa37c7051bbbdc97b1  candidates/tensorization_exact_audit.json
92659ae392a18b1146a1520dc63c80d28d70469aad16f30499dd2cfdb0cfd7cc  candidates/square_direct_audit.json
87dd06084361e14832c348b6a97fd433c8c58f16438d2895d432ba8aa9adbfb2  candidates/n6_square_direct_audit.json
287ee465e2376f12441641bcc5795463d0b867dd27821135531f43cb8ecdca0a  candidates/power_admissibility_audit.json
9c147dba19b8d3f00055c43917927cb5f70b07daa87d94528fad521564d09882  certificates/gate_b_unbounded_arb_v4.json
d4019a9bca35e0945ab75c05c9a72dcf52311f022429d179f1b193acadef937b  certificates/gate_b_unbounded_dyadic_v1.json
1c39d372bde1086da7f97099c02f99783d1576b9c10d72e033b2d8f38354613a  certificates/gate_b_unbounded_n6_dyadic_v1.json
3b4650fc00e17817cd34824643533320bd6718ea4cb341d404ff2587030f94da  certificates/gate_b_unbounded_n6_arb_v1.json
0ddb6ad4c61007a9059d1df205ae1f972d6b7ca551f56eccbc24da685d002816  certificates/gate_b_unbounded_rational_v1.json
90fa33d4849c568bd686941fc951e11bdc6aebc64aa94202cf4bc6d14f6631c6  experiments/n5_complete_checkpoint.jsonl
012051ad8feb2aaac3f410318957f14e9c9836c9fad7111c6fdc91940d52ed47  experiments/n8_k1_checkpoint.jsonl
2c46e8dae7b8268bd28e10357374617a21e434a3da44da1e8bac7136058cdd9f  experiments/tensorization_audit_checkpoint.jsonl
45644badf44ef8507130ece9c14208b0d2fa638970b394fc986cf9042f93fcc7  experiments/tensorization_exact_checkpoint.jsonl
e2101226694ed579f70bb772cefefa8589603683b25c8a92cc2967449c7784b3  experiments/rational_certificate_checkpoint.jsonl
0599d2b038a82fa7298c4ca0d3f5d8b6d6761451271591e4dac449037b07ac39  experiments/square_direct_checkpoint.jsonl
fcc0c916413271dc269fa9dbd1f5a50819fc972eeed651dce7d58fbde31ef775  experiments/n6_square_direct_checkpoint.jsonl
2695bb03ae5d5a5a222b71201330b6c37f0d2306bdb0db49e35248062d522cb1  experiments/power_admissibility_checkpoint.jsonl
1a8596c93a5fdea34f1dd4e1a3335046a9c3f32438383f7ef623bb9948a241f9  PROOF.md
2d295690e1159e60859401ab96433e810d778b409e91afe4edb01dd8ca4f08d8  DEFINITIONS.md
54c35eb40f2b12bdfb4ca745f6db14aa5a4f338bb9d5d845eee35fef70373518  EXPERIMENTS.md
472448dc0eeb28e077fc6afb52c0d775f04e9cfa0b98d4b8db1ef64d098e7f14  paper/main.tex
4b57f1ad03f9976995c1d5022225681ba8aacc25ef3a651b4202c0de3004eefd  paper/main.pdf
e207d5e9e452e02d7205c58691a47faca8cb23fd28da56275bfc912faf6e8c20  search_local_defect.py
8f9b22b439761a3eab7fb84822bc9276c2f01d0de10097a1f0ef2938effc4db9  hunt_local_ratio.py
d33e2b1ecc6d83a89ed5f12b591c4f9c6f552454c6146434fc67bc76c6914175  certificates/gate_b_lowdefect_rational_v1.json
a74724856c5e0ac4de4ac4f726cecd654176374bceffe5580c3ccf2b074d3c8c  experiments/local_defect_checkpoint.jsonl
e44a0bb9be2e834db4abaf0dc56fa01d95d1083cfb5b5b0b1753e8d6780f772d  experiments/local_ratio_checkpoint.jsonl
8e20e4f3a883e68ed53b4e5b037ee6885d0f6bd18f99f0dc803e75294cc0498b  candidates/local_defect_frontier_n7.json
8d8a3c9057c0565f675ae5701534f81e57bf98ef103f852d60476273fbfc5c4a  candidates/local_sweep_n7.json
6c73a58a563736a87c4b40316aee4e31d3c600833666a6d9286ced44df31fdb0  candidates/local_descend_n7.json
da4b6ca9327ca2035530fa44198e173d9d29d09bb8180b12d9ff1d6494972fd8  bound_local_regime.py
ce58cfdf129bfe7219a67d2c943c3551d786c2d5a96832ccc6bf3bf9914183d1  certify_class_negatives.py
b6ad2baa8bd74ca1fc8d5459a7d7c91948d0a89baff1090d053576827672c9ca  audit_clone_power.py
cc63961f987ef54195a1fe7d7ba3586bb4b2eca04718ae0e54d0294b4a8e2827  CURRENT_STATUS.md
d5f6a658329c238370ee7995f486d24909977cb4aa1e19fd785d2cf14a5c2f2b  certificates/gate_b_n8_rational_v1.json
b57909e5ef81f3b6933fa53951905b3f4606ed8b143255e01202ffec3a3f457e  certificates/gate_b_n8_clone_rational_v1.json
79d8208483fe1711336933b2f9122bc3dc3556743c6704f33befebf4d9aeebab  certificates/gate_b_n8_max_rational_v1.json
1135b741eadafd8efd6185b0b3c253c2af345b87f48b45dafe7bcbd9bbe60d8c  certificates/gate_b_n8_k4_rational_v1.json
f8e30d7ad468f9a8531886319e5db6a57ac9bfb47c0409d0662651fd317a0a5e  candidates/local_regime_bounds.json
6e87bdb0a8de4c05f60a990651f15d85a0f55181e9dfdb936666b706eca157ba  candidates/clone_power_audit.json
daa3fbe2fec959e0f159073ab4ddfe41b42ec0f3d0c4bc05f52b46dc24482740  candidates/n8_k2_census.json
5cc25f5187baa9925cb2716f5098e352dacad4d0ac59e95bc175caefe0c588fd  candidates/n8_k3_census.json
d72f074e28e99ca341cf65bc4a065822e65f95f8c98c4ae2cc659da25a58e661  candidates/n8_k4_census.json
d02d386978ebe8dd0f61e13afc188e60ceaddd9ee441a75ee36453093cfc37ce  candidates/n8_k1_exact_frontier.json
d6d4cccf5504f035e7b753654ab140786141bae438d8c3db4cb3bd684d450d16  candidates/n8_k2_exact_frontier.json
b7b9e46def0924f00401b272ad4fa2f33246ee8a4e32e9c40ee08354a97455d1  candidates/n8_k3_exact_frontier.json
98161d4c94df67d0fffa6bd12ffff3a9be92c39029755d7a46502f7bbded5743  candidates/n8_k4_exact_frontier.json
5ddc0b92c3c6360ced3d6515825090c377ea44f36b378a774a25c83e14ceb3db  experiments/n8_k2_checkpoint.jsonl
895858365dc5713f8743758ff1be36c9babb5d9a0980c8404a21e45ad5af9aa4  experiments/n8_k3_checkpoint.jsonl
e4a0656319778405e8d7caa6f4251a786961b36bb1532289cb4b7ab5d0347d3f  experiments/n8_k4_checkpoint.jsonl
59b7f7d60ef7c1b8b0995f1a2c2f6793429a22e1d7d1dc85e2e3b2f08f4468bf  experiments/n8_k1_exact_checkpoint.jsonl
55a3c155ef97a75e2e0a0b3a35948fa38a15f34428cab7ea776ba8aa394d357e  experiments/n8_k2_exact_checkpoint.jsonl
f93f157b834ee6841ff64f4f703c5f42d41b0a7bbc36db03f7c711efa9595b60  experiments/n8_k3_exact_checkpoint.jsonl
35ee76ed0dcfc713567762ca2ed8e24bad6d9d582dfa711b05333de8d00e5d8a  experiments/n8_k4_exact_checkpoint.jsonl
94939dd84e291fbcc434d0d9184998eb94da567305b45956f4924c12c6cfd05b  audit_clone_saturation.py
8c7621a039c976aee624a5e91064f0209885fd32f5994cb9d735bf172268ab4f  candidates/clone_saturation_audit.json
```

## Reproduce in the existing environment

Run from `math/`. Every expensive command is single-threaded and low priority.

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0

nice -n 10 ./.venv/bin/python -B uc/gate_b/test_gate_b.py
nice -n 10 ./.venv/bin/python -B uc/gate_b/verify_gate_b.py
nice -n 10 /Library/Developer/CommandLineTools/usr/bin/python3 -B \
  uc/gate_b/verify_gate_b_dyadic.py
nice -n 10 /Library/Developer/CommandLineTools/usr/bin/python3 -B \
  uc/gate_b/verify_gate_b_n6_dyadic.py
nice -n 10 ./.venv/bin/python -B uc/gate_b/verify_gate_b_n6_arb.py
nice -n 10 /Library/Developer/CommandLineTools/usr/bin/python3 -B \
  uc/gate_b/verify_gate_b_rational.py --bases n6,n7 --orders all
nice -n 10 /Library/Developer/CommandLineTools/usr/bin/python3 -B \
  uc/gate_b/verify_gate_b_rational.py --bases n7lo --orders all
nice -n 19 /Library/Developer/CommandLineTools/usr/bin/python3 -B \
  uc/gate_b/verify_gate_b_rational.py --bases n8clone_lo,n8clone_hi,n8max \
  --orders orbits
nice -n 19 /Library/Developer/CommandLineTools/usr/bin/python3 -B \
  uc/gate_b/audit_clone_power.py --base n8clone_hi
nice -n 19 /Library/Developer/CommandLineTools/usr/bin/python3 -B \
  uc/gate_b/audit_clone_saturation.py --max-clones 4
nice -n 19 ./.venv/bin/python -B uc/gate_b/search_n8_block_symmetric.py --k 3 \
  --checkpoint uc/gate_b/experiments/n8_k3_checkpoint.jsonl \
  --result uc/gate_b/candidates/n8_k3_census.json
nice -n 19 /Library/Developer/CommandLineTools/usr/bin/python3 -B \
  uc/gate_b/certify_class_negatives.py \
  --census uc/gate_b/experiments/n8_k3_checkpoint.jsonl --threshold 0.01 \
  --checkpoint uc/gate_b/experiments/n8_k3_exact_checkpoint.jsonl \
  --report uc/gate_b/candidates/n8_k3_exact_frontier.json
nice -n 19 ./.venv/bin/python -B uc/gate_b/search_n8_block_symmetric.py --k 4 \
  --checkpoint uc/gate_b/experiments/n8_k4_checkpoint.jsonl \
  --result uc/gate_b/candidates/n8_k4_census.json
nice -n 19 /Library/Developer/CommandLineTools/usr/bin/python3 -B \
  uc/gate_b/certify_class_negatives.py \
  --census uc/gate_b/experiments/n8_k4_checkpoint.jsonl --threshold 0.01 \
  --checkpoint uc/gate_b/experiments/n8_k4_exact_checkpoint.jsonl \
  --report uc/gate_b/candidates/n8_k4_exact_frontier.json
nice -n 19 /Library/Developer/CommandLineTools/usr/bin/python3 -B \
  uc/gate_b/verify_gate_b_rational.py --bases n8tiny,n8best --orders orbits
nice -n 10 /Library/Developer/CommandLineTools/usr/bin/python3 -B \
  uc/gate_b/audit_tensorization_exact.py --quiet
nice -n 19 /Library/Developer/CommandLineTools/usr/bin/python3 -B \
  uc/gate_b/bound_local_regime.py --bases n6,n7,n7lo --max-dimension 4
nice -n 19 /Library/Developer/CommandLineTools/usr/bin/python3 -B \
  uc/gate_b/verify_gate_b_rational.py --bases n8lo,n8hi --orders orbits
nice -n 19 ./.venv/bin/python -B uc/gate_b/search_n8_block_symmetric.py --k 2 \
  --checkpoint uc/gate_b/experiments/n8_k2_checkpoint.jsonl \
  --result uc/gate_b/candidates/n8_k2_census.json
nice -n 19 /Library/Developer/CommandLineTools/usr/bin/python3 -B \
  uc/gate_b/certify_class_negatives.py \
  --census uc/gate_b/experiments/n8_k2_checkpoint.jsonl --threshold 0.01
nice -n 10 ./.venv/bin/python -B uc/gate_b/search_local_defect.py \
  --mode exact --dimension 4
nice -n 19 ./.venv/bin/python -B uc/gate_b/search_local_defect.py \
  --mode frontier --dimension 7
nice -n 10 ./.venv/bin/python -B uc/gate_b/audit_n6_square_direct.py
nice -n 10 ./.venv/bin/python -B uc/gate_b/search_n8_block_symmetric.py --k 1
nice -n 10 c++ -O3 -std=c++17 \
  uc/shapley_n5_global_coupling_enumerate.cpp \
  -o /tmp/shapley_n5_global_coupling_enum
nice -n 10 ./.venv/bin/python -B uc/gate_b/audit_n5_complete.py \
  /tmp/shapley_n5_global_coupling_enum
nice -n 10 ./.venv/bin/python -B uc/gate_b/audit_tensorization.py
nice -n 10 ./.venv/bin/python -B uc/gate_b/audit_square_direct.py
nice -n 10 ./.venv/bin/python -B uc/gate_b/audit_power_admissibility.py
```

The \(n=8\), \(n=5\), tensorization, square, and exact-rational runs resume from
append-only checkpoints and must report zero new evaluations when the checked-in
checkpoints are complete. The rational verifier also recomputes a deterministic
sample of stored records on every resume and aborts on any mismatch.

The local-regime searches are resumable the same way, one append-only record per
completed restart, and every emitted family is re-checked for cap, exact integer
incidence, distinct rows and activity before it is allowed to set a bound. Their
coverage label is **SEEDED HEURISTIC** except `--mode exact`, which enumerates
every admissible family and is complete for \(n\le4\). `--mode frontier` reports
the best re-verified witness per size across every record schema, because schema
strictness governs work skipping rather than the validity of a witness.

Reproduce the original Gate A computations independently:

```sh
nice -n 10 ./.venv/bin/python -B uc/shapley_n7_falsifier_cert.py
nice -n 10 ./.venv/bin/python -B uc/shapley_n7_block_symmetric.py
nice -n 10 ./.venv/bin/python -B uc/shapley_adaptive_coupling.py
nice -n 10 c++ -O3 -std=c++20 uc/shapley_n6_small_enumerate.cpp \
  -o /tmp/shapley_n6_small_enum
nice -n 10 ./.venv/bin/python -B uc/shapley_n6_small_complete.py \
  /tmp/shapley_n6_small_enum
nice -n 10 ./.venv/bin/python -B uc/shapley_n6_block_symmetric.py
```

## Clean-environment verifier

The clean-room run verified on this workstation uses the system Python 3.9
interpreter, for which `python-flint==0.6.0` has a published arm64 wheel.

```sh
cd /path/to/jinleic-workspace/math
/Library/Developer/CommandLineTools/usr/bin/python3 -m venv \
  /tmp/uc-gate-b-py39
/tmp/uc-gate-b-py39/bin/python -m pip install 'python-flint==0.6.0'
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  nice -n 10 /tmp/uc-gate-b-py39/bin/python -B \
  uc/gate_b/verify_gate_b.py
```

The expected terminal verdict is `PROVED_GATE_B_UNBOUNDED`. The verifier prints
the certified base, direct-square combinatorics, actual square Bellman checks,
and structured all-power admissibility witnesses in JSON. Python 3.14 is not a
clean-install path for the frozen `python-flint==0.6.0` release; its source
distribution is incompatible with current Cython. Use the verified Python 3.9
wheel path above rather than silently changing the arithmetic dependency.

The independent certificate needs no third-party package at all. It runs on the
pristine system interpreter with an empty environment:

```sh
OMP_NUM_THREADS=1 nice -n 10 \
  /Library/Developer/CommandLineTools/usr/bin/python3 -B \
  uc/gate_b/verify_gate_b_dyadic.py
```

That command reproduced `certificates/gate_b_unbounded_dyadic_v1.json`
byte-identically in 4.0 s, so the unboundedness conclusion has one replay path
with no numerical dependency whatsoever.

The same pristine interpreter also reproduces the second-base certificate.
That checker evaluates all 720 six-coordinate orders rather than trusting an
automorphism quotient and proves \(A_+(\mathcal D)<-1/80\).

The exact rational certificate has the same property and needs no interval
arithmetic at all. In a stripped environment:

```sh
env -i OMP_NUM_THREADS=1 HOME="$HOME" PATH=/usr/bin:/bin nice -n 10 \
  /Library/Developer/CommandLineTools/usr/bin/python3 -B \
  uc/gate_b/verify_gate_b_rational.py --bases n6,n7 --orders all --fresh \
  --write-certificate /tmp/rational_replay.json
```

That reproduces `certificates/gate_b_unbounded_rational_v1.json` byte-identically
in 250 s. Use `--fresh` for the byte comparison: the certificate honestly
records whether each order was newly evaluated or replayed, so a resumed run
differs from the canonical fresh run in exactly those bookkeeping fields
(`new_evaluations`, `resumed_records`, `rechecked_records`) and in nothing
mathematical. Dropping `--fresh` on the committed checkpoint is the resumability
check: it must print `new_evaluations: 0` and recompute its sampled records
without a mismatch, which takes about one second.

## Evidence labels

- **PROVED:** product identities, admissibility of every power, and divergence.
- **CERTIFIED:** exact base constraints and the 256-bit Arb upper bound
  \(A_+(\mathcal B)<-7/250\).
- **INDEPENDENTLY CERTIFIED:** the standard-library dyadic relaxation proves
  \(A_+(\mathcal B)<-1/40\) without Arb or Bellman clamp classification. Its
  exact product consequence gives ratio \(>k/40\to\infty\).
- **CERTIFIED SECOND BASE:** the generic 256-bit Arb recurrence independently
  proves \(A_+(\mathcal D)<-17/1250\), giving ratio
  \(>(17k/1250)/(1-(181/625)^k)\).
- **INDEPENDENTLY CERTIFIED SECOND BASE:** all 720 orders of the distinct
  25-row \(n=6\) family give the standard-library dyadic relaxation
  \(A_+(\mathcal D)<-1/80\), and exact product arithmetic gives ratio
  \(>k/80\to\infty\).
- **EXACTLY CERTIFIED, BOTH BASES:** exact rational two-sided enclosures over
  all 720 and all 5040 coordinate orders against the true feasible action set,
  with no interval library and no relaxation. They give
  \(A_+(\mathcal D)<-17/1250\) and \(A_+(\mathcal B)<-7/250\), i.e. the sharp
  form of both base lemmas from one standard-library artifact, and they contain
  the Arb values.
- **PROVED (growth rate):** \(c_{\rm cl}^\star(n)>n/250-7/250\) from the powers,
  and \(-A_+\le\log_2m\le n\) from \(Q\ge0\), \(C_+\ge0\), so the ratio is
  \(\Theta(n)\) on families with defect bounded below.
- **COMPLETE EXACT / COMPLETE NUMERICAL:** the full nontrivial \(n=5\)
  cap/Reimer census,
  the declared \(n=6\) controls, and the declared \(n=8\)
  \(S_1\times S_7\) class. Enumeration, filters, defects, and stated orbit
  coverage are exact; objective comparisons are float64.
- **OPEN:** the newly separated local-stability problem with
  \(\varepsilon_\vee\to0\), and any non-scalar UC-specific charge capable of
  advancing the Frankl programme.
