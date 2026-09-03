# Research status: UC certificate and Liu Hypothesis 1

**Status date:** 2026-08-29  
**Scope:** the union-closed explicit-constant candidate, the Gate B
closure-defect programme, the Liu H2 scalar-margin reduction, and e389
Input-II constants.

## Evidence vocabulary

- **MACHINE-VERIFIED:** the named executable check completed on identified bytes.
- **HUMAN-AUDITED:** the mathematical argument was independently reconstructed with its hypotheses and endpoints.
- **COMPUTATIONAL-EVIDENCE:** numerical, sampled, truncated, or finite-dimensional evidence that is not a continuum proof.
- **CITED-DEPENDENCY:** an identified external theorem is used without formalization here.
- **OPEN:** a claim or external review remains unresolved.
- **FAILED:** a check or wording claim was rejected; it contributes no positive proof evidence.

## UC candidate bound

### Human-audited result

For the exact rationals

\[
\alpha_0=\frac{356069}{10^7},\qquad
 t_{\rm cert}=\frac{955165028125263}{2500000000000000}
 =0.3820660112501052,
\]

the audited chain is:

1. symmetric self-couplings fold exactly to pair-orbit laws;
2. the exact functional attains its minimum on at most two pair-orbits;
3. the Margin Lemma gives the explicit relaxed functional \(\Phi_{\rm rel}\), with \(\Phi_{\rm ex}-\Phi_{\rm rel}=(1-\alpha)(S^2-E)\ge0\);
4. the entropy-zero case \(L=0\) gives equality zero and requires no quotient;
5. the finite certificate proves \(\Phi_{\rm rel}\ge0\) on the feasible five-parameter domain;
6. a self-contained sequential Bernoulli-coupling and entropy-chain proof converts the functional inequality to the union-closed frequency conclusion.

No fatal circularity, inequality-direction error, uncovered feasible region, endpoint omission, or float/rational unsoundness was found. The resulting theorem remains a **candidate** pending ordinary external mathematical review; the analytic chain is not end-to-end machine checked.

### Machine verification

- Frozen Campaign I: 488,465,854 trace events, eight COMPLETE slices, zero residual/stack/budget fields.
- Independent structural audit: all trace bytes, opcodes, DFS topology, termination, and tallies passed; report [`uc/verification/results/structural-traces.json`](uc/verification/results/structural-traces.json), canonical SHA-256 `7270e5418be72b09e607457438bc9f8109154c31b491abca9546f4ebee3cf5a3`.
- Independent point-mass strictness: 256-bit Arb proves the required endpoint margin positive; report [`uc/verification/results/cambie-bridge-strictness.json`](uc/verification/results/cambie-bridge-strictness.json).
- External campaign identity is pinned by [`uc/verification/campaign-lock.json`](uc/verification/campaign-lock.json), raw SHA-256 `d254a23a4df7cb3a4c6aae2453883324fba2381c1890c4f011d61c7971b84b48`.
- Fresh direct arithmetic replay is complete: all eight isolated processes exited 0, every expected stdout hash/tally matched, residual was zero, and the total was 488,465,854. Pre/post structural input maps are identical. Summary: [`uc/verification/results/direct-isolated-replay/summary.json`](uc/verification/results/direct-isolated-replay/summary.json), canonical SHA-256 `6ef126ced337b035e5c5f22a130b1be3697666e60f16586b75540ef6fc734b16`.
- The hardened fresh-copy clean-room replay also passed all eight slices. Its aggregator accepted 488,465,854 total nodes and zero residuals: [`uc/verification/results/clean-room-arithmetic/composite.json`](uc/verification/results/clean-room-arithmetic/composite.json), canonical SHA-256 `57ca3f52c054bd9bccfefc349f440674c74fe7ef671d06c103ef3f52dc33d53d`. An earlier four-minute batch was stopped with exit 143 after verifier hardening; no result from it is accepted.
- Secure source-independent arithmetic replay is complete.  Eight byte-zero
  workers imported no frozen arithmetic, independently rederived centered
  gradients by interval AD, rediscovered every face proof, and reproduced all
  `488465854` events, `244232923` splits, `244232931` leaves, and zero
  residuals.  Trusted report-lock raw SHA-256:
  `af86480901c2c497739916bb410f1e117e4eaf3eefb82575ce494b7d8c04e730`;
  accepted composite:
  [`composite.json`](uc/verification/results/independent-arithmetic/secure-full/composite.json),
  canonical SHA-256
  `4acbd3b935bda7e51ed387e42e0598debf22f4a85b7a24ad976c3eaca4f23243`.
- An independent read-only analytic referee found no blocker in the human UC
  chain.  Two major exposition defects and two minor precision defects were
  repaired; preserved report:
  [`uc-analytic-referee-2026-08-26.md`](uc/verification/results/uc-analytic-referee-2026-08-26.md).
- The exhaustive entropy-bridge control is complete.  A verifier written from
  the corollary statement alone — importing no campaign, certificate, or
  bridge source — rebuilds the sequential coupling and exhausts **every**
  nonempty family on \(n\le4\) coordinates: 65,808 families, 1,631,880 coupled
  prefix states.  Exact rational arithmetic confirms the coupling masses, both
  Bernoulli marginals, the \(s^\star\) OR identity, prefix-by-prefix
  uniformity of \(A\) and \(C\), \(\mathcal L(P_i)=\mathcal L(R_i)\) with mean
  the element frequency, and the chain rule; 256-bit Arb certifies both
  data-processing inequalities, every enclosure being provably nonnegative or
  a structural tie inside \(\pm2^{-200}\).  All 5,096 enumerated union-closed
  families satisfy the \(t_{\rm cert}\) conclusion, the extreme case being
  exactly \(1/2\).  Report:
  [`entropy-bridge-exhaustive.json`](uc/verification/results/entropy-bridge-exhaustive.json),
  canonical SHA-256
  `6681f8faf13d9344f8e7bf3a5d7785ed289f1e82c5e9badd1e7967c03e3ebfef`.
  Mutation tests (`uc/verification/test_entropy_bridge_exhaustive.py`, 4/4)
  show the control rejects a broken \(s^\star\) clip, a constant prefix
  probability, and an inflated union bit.
- A second independent read-only referee re-derived the two remaining human
  links — the sequential entropy bridge and the fold/unfold plus Theorem
  B\('''\) support reduction — and returned **no proof defect and no
  blocker**.  Eight minor exposition defects were found and all eight were
  repaired in `uc/paper/main.tex`: the empty-prefix justification for the
  constant first-coordinate probability, the enabling inequality
  \((1-\psi)^2=\psi\Rightarrow 2u-u^2\ge1-\psi>1/2\) behind \(g'<0\), the
  previously unstated three-case \(s^\star\) split, the Stone--Weierstrass
  span-versus-tensor slip, a cleaner continuity route for \(\Phi_\alpha\), the
  explicit PSD cross-term derivation of concavity, the \(\Phi_\alpha\)/\(F\)
  sink-law distinction, and the exclusion of the empty family.  Preserved
  report:
  [`uc-bridge-referee-2026-08-27.md`](uc/verification/results/uc-bridge-referee-2026-08-27.md).

### UC unresolved obligations

- **OPEN (external blocker):** ordinary external specialist/journal review of
  the human support reduction, Margin Lemma, endpoint estimate, and
  entropy-chain proof.  Two independent read-only model referees have now
  found no blocker and no proof defect; neither is peer review, and this
  cannot be closed inside the repository.
- **RESOLVED as a finite control:** the entropy-to-union-closed bridge is
  exhaustively machine-checked on every family with \(n\le4\) coordinates.
  The universal statement remains the manuscript induction.
- **RESOLVED subject to explicit trust:** an independently implemented
  interval verifier completed all eight traces from byte zero.  It shares the
  authenticated trace partition and Arb primitives, but no frozen Python
  arithmetic or derivative formulas.
- **PARTIALLY RESOLVED:** a 343-file content-addressed runtime seal plus live
  eight-process image attestation binds the observed interpreter/native bytes.
  Hardware-backed attestation, hostile same-inode-race exclusion, and a
  malicious OS/kernel/tool stack remain outside the claim.
- **OPEN:** formal proof-assistant verification, optimality of \(t_{\rm cert}\), any bound at \(c^*\), and Frankl's \(1/2\).
- **COLLECTION IN PROGRESS (2026-08-29):** Campaign K at
  \(t=0.3822660112501052\ge\psi+3\times10^{-4}\) now has **all 16 of 16
  slices committed COMPLETE** with residual = stack = budget_time = 0
  (805,542,368 processed boxes; the three slices lost to the 2026-08-22
  broker disconnection — 10, 12, 13 — were relaunched 2026-08-27 under the
  *same* launch manifest and finished 2026-08-28; the frozen commit path
  cannot replace an existing artifact). The frozen collector
  (`6a2ae5ea065ffbc1…`) is replaying all sixteen traces as of 2026-08-29
  (first observed line: `REPLAY slice=0 PASS` in 1212 s). A
  `COMPOSITE CERTIFICATE` there would raise the certified constant to
  \(0.3822660112501052\); every earlier partial-collection record, including
  the pre-completion rejection naming exactly the three missing
  result/trace pairs, is preserved verbatim in
  [`data/uc-certificate-psi-plus-3e-4-2026-08-29/`](/Users/jinleic/jinleic-workspace/data/uc-certificate-psi-plus-3e-4-2026-08-29/collect_rejection_precompletions.txt).
  Nothing weaker than that collector's exit-0 certificate line may be
  reported as a certificate.

### UC supporting artifacts

- Full proof audit: [`uc/AUDIT.md`](uc/AUDIT.md)
- Reproduction contract: [`uc/REPRODUCIBILITY.md`](uc/REPRODUCIBILITY.md)
- Verification tooling/results: [`uc/verification/`](uc/verification/)
- Revised manuscript: [`uc/paper/main.tex`](uc/paper/main.tex), [`uc/paper/main.pdf`](uc/paper/main.pdf)
- Originality report and sources: [`uc/LITERATURE_ORIGINALITY.md`](uc/LITERATURE_ORIGINALITY.md), [`uc/literature/`](uc/literature/)
- Exact proof/evidence chain: [`uc/PROOF.md`](uc/PROOF.md)

## Gate B closure-defect obstruction

Separate from the candidate bound above. With the repository's frozen
definitions, \(c_{\rm cl}^\star=\sup\{-A_+/\varepsilon_\vee\}=+\infty\).
There are now two combinatorially distinct Cartesian-power constructions:
the normalized 45-row \(n=7\) family \(\mathcal B\), with closure success
\(17/81\), and the normalized 25-row \(n=6\) family \(\mathcal D\), with
closure success \(181/625\).  Both stay cap/Reimer-admissible while \(A_+\) is
additive and closure success is multiplicative.

- **PROVED:** the product lemmas, admissibility of every power, and divergence
  for both sequences.
- **BOTH BASES CERTIFIED THREE TIMES:** 256-bit Arb gives
  \(A_+(\mathcal B)<-7/250\) and \(A_+(\mathcal D)<-17/1250\).
  Standard-library dyadic interval checkers with clamp-free relaxations give
  \(A_+(\mathcal B)<-1/40\) and \(A_+(\mathcal D)<-1/80\).  The latter
  evaluates all 720 orders without a symmetry quotient.  Both dependency-free
  certificates reproduce byte-identically on the pristine system interpreter.
  Thus independently implemented interval verification is satisfied for both
  Gate B and, under its stated trace/Arb/runtime/custody trust boundary, the
  Campaign I certificate.
- **CERTIFIED EXACTLY (2026-08-27):** an exact rational evaluator with no
  third-party package, no interval class, no clamp classification and no
  action-set relaxation pins the *true* objective over **all** 720 and **all**
  5040 coordinate orders:
  \(A_+(\mathcal D)\in[-0.013672107732177773561859915661,
  -0.013672107732177773561859912998]\) and
  \(A_+(\mathcal B)\in[-0.028649186794468317816026981933,
  -0.028649186794468317816026976455]\), widths below \(6\times10^{-27}\).
  The upper endpoints re-derive the sharp \(-17/1250\) and \(-7/250\) bounds,
  and both Arb values lie strictly inside.  Its all-order value lists partition
  into 10 classes of size 72 and 21 of size 240, recovering the automorphism
  data instead of assuming it.  Certified primitives: bit-by-bit \(\log_2\),
  ceiling square-root towers for \(2^x\), and a concavity case split.
- **MACHINE-CHECKED SUPPORT:** symbolic verification of the shared Bellman
  algebra, non-factored direct-square audits for both bases over five global
  orders each, a 292-case float product falsification search with 205
  unequal-dimension instances, and a 35-case exact rational repeat over all
  15,010 global orders whose worst fixed-order Bellman discrepancy is exactly
  zero.
- **PROVED (growth rate, closed 2026-08-27):** \(c_{\rm cl}^\star(n)>n/250-7/250\)
  from the powers. In the other direction incidence is at most
  \(n\lfloor2m/5\rfloor\) while Reimer demands at least \(m\log_2m/2\), so
  **every** admissible family obeys the size ceiling \(\log_2m\le4n/5\) — the
  exact integer test \(m^5\le2^{4n}\), verified for every admissible size at
  every \(3\le n\le16\) and close to sharp (\(8.7977\) against \(8.8\) at
  \(n=11\)). Hence \(-A_+\le4n/5\) unconditionally and the numerator's growth
  is \(\Theta(n)\), bounded on both sides; only the denominator is open. The
  lower constant improved on 2026-08-27 by exactly \(75/64\), from slope
  \(1/250\) to \(3/640=0.0046875\), giving
  \((3/640)n-3/80\le\Lambda(n)\le4n/5\). The improvement comes from a
  cloned-coordinate base: **no lemma in the product argument uses separation**,
  and Reimer for every power reduces to the single base-level integer inequality
  \(m^m\le2^{2I}\), which is exactly \(I\ge R_m\). Powers audited by direct
  instantiation, verdict `EVERY_POWER_ADMISSIBLE`.
- **CERTIFIED (2026-08-27), complete \(n=8\) block-symmetric frontier:** all four
  \(S_k\times S_{8-k}\) classes enumerated; 52,494,332 raw masks, 27,431
  canonical admissible families, 150 screen-negative and 147 exactly negative.
  Current records: defect \(14/45\) in the displayed class and \(1144/1875\)
  under separation; ratio \(147/2536\) displayed and \(531/11360\) separating;
  slope \(3/640\) displayed and \(177/40000\) separating. All ten registered
  bases plus the finite symbolic \(14/45\) clone family have exact rational
  certificates.
- **PROVED (2026-08-27), defect floors:** the cap forces
  \(\mathbb E|X\vee Y|\ge\frac85\bar s\) while a successful union has size at
  most \(M\), giving \(\varepsilon_\vee\ge(\frac85\bar s-M)/(n-M)\); downset
  counting gives \(\varepsilon_\vee\ge1-m^{-2}\sum_AN(A)^2\), attained with
  equality. Hence the ratio is \(O(n)\) without a dominant set, so the local
  regime requires a set of size \(\ge\frac45\log_2m\). Also: 25 and 45 are the
  largest admissible sizes at \(n=6,7\), all three \(n\le7\) bases attain them
  with every degree at the cap, and **no set can be added to any of them**.
- **PROVED (2026-08-27), local lemmas:** failures come in pairs, so
  \(\varepsilon_\vee\ge2/m^2\) and \(-A_+/\varepsilon_\vee\le m^2\log_2m/2\); the
  defect multiplies through success, so **no Cartesian power of either base has
  defect below its base** and the theorem's mechanism cannot enter the local
  regime at all.
- **PROVED/CERTIFIED (2026-08-28), sharp clone closure:** an appended coordinate
  preserves defect iff its one-rows are join-consistent on every successful
  join. Clone multiplicities induce the exact Plackett--Luce first-appearance
  law; clone infimum/supremum equal fixed-order minimum/maximum, and a finite
  clone multiset is negative iff one fixed order is negative. This retracts the
  earlier universal-decrease and geometric-convergence claims.
- **EXACT WITNESS / HEURISTIC SEARCH COVERAGE (2026-08-28), defect
  \(14/45<2/5\):** the corrected minimum-fixed-order search found the normalized
  core
  \[
  (0,1,2,4,5,8,10,43,64,190,192,193,245,254,255)\subseteq2^{[8]}
  \]
  with degrees \((6,6,6,6,4,5,6,6)\), incidence \(45\ge R_{15}=30\), full set,
  and defect \(14/45\). Its average is positive,
  \(A_+\in[0.230029313333,0.230029313334]\), but order
  \((6,1,2,0,3,4,5,7)\) is exactly negative:
  \([-0.003540262562,-0.003540262561]\).
- **CERTIFIED EXACTLY (2026-08-28), negative below \(2/5\):** ratio-41
  multiplicities
  \[
  (2825761,4750104241,115856201,68921,1681,41,194754273881,1)
  \]
  give a finite symbolic family with dimension \(199623130728\), \(m=15\),
  incidence \(1197738780965\), defect \(14/45\), and
  \(A_+\in[-0.000084288238,-0.000084288237]\). All 40320 original fixed-order
  enclosures are reweighted exactly; float and sampled values only selected the
  candidate/order. Certificate verdict
  `PROVED_SUB_TWO_FIFTHS_NEGATIVE_EXACT_RATIONAL`.
- **Convention:** the core is active and separating, but the cloned witness
  repeats columns. Thus \(14/45\) is the displayed cap/Reimer record; the
  separating negative frontier remains \(1144/1875\).
- **The earlier \(14/45\) extension closure remains exact.** The first
  defect-minimizing core has all 5040 orders positive, and its unique new
  join-consistent extension has all 40320 orders positive. The obstruction was
  family-specific; changing rows exposed the negative order.
- **REIMER ENDPOINT (2026-08-28):** Chase--Lovett [REPORTED] refute cap +
  normalization + dominant-set floors and every pair-defect-continuous
  approximate Reimer inequality with vanishing normalized error:
  \(h(\psi)/2-\psi>0.0977433528\) while \(\varepsilon_\vee\to0\). At exact zero
  defect, Reimer's cited average-set-size theorem and the dominant set are
  automatic. A floor reaching zero is still the \(2/5\) frequency theorem.
- **OPEN:** the new witness has fixed defect \(14/45\), not a sequence tending
  to zero; \(c_{\rm loc}\), the separating frontier below \(1144/1875\), and
  proof-assistant formalization remain open.
- **METHODOLOGICAL:** six of ten registered bases satisfy
  \(\max_\pi Q_\pi<\log_2m\), an exact every-order iid-chain barrier. Cambie's
  two-atom upper obstruction is exactly enclosed at
  \(0.3823455333667027\ldots\), but its matching lower direction remains
  **OPEN / COMPUTATIONAL-EVIDENCE** in the broader UC audit.
- Artifacts: [`uc/gate_b/CURRENT_STATUS.md`](uc/gate_b/CURRENT_STATUS.md),
  [`uc/gate_b/README.md`](uc/gate_b/README.md),
  [`uc/gate_b/PROOF.md`](uc/gate_b/PROOF.md),
  [`uc/gate_b/EXPERIMENTS.md`](uc/gate_b/EXPERIMENTS.md),
  [`uc/gate_b/bound_local_regime.py`](uc/gate_b/bound_local_regime.py).
  New local artifacts:
  [`uc/gate_b/barrier_sawin_cambie.py`](uc/gate_b/barrier_sawin_cambie.py),
  [`uc/gate_b/audit_clone_limit.py`](uc/gate_b/audit_clone_limit.py),
  [`uc/gate_b/search_subtwofifths_q.py`](uc/gate_b/search_subtwofifths_q.py),
  [`uc/gate_b/certificates/gate_b_subtwofifths_clone_rational_v1.json`](uc/gate_b/certificates/gate_b_subtwofifths_clone_rational_v1.json),
  and the three
  [`subtwofifths_min_a`](uc/gate_b/candidates/subtwofifths_min_a_n9.json)
  search reports indexed in the Gate B README.

## Liu Hypothesis 1

### Human-audited theorem

Define

\[
R(s,t)=(1-s)(1-t)[st-(1+st)\ln(1+st)]
-[z+(1-z)\ln(1-z)],
\quad z=(1-s)(1-t)(1+st).
\]

For every finite real signed Borel measure \(\nu\) on \([0,1]\),

\[
\iint R(s,t)\,d\nu(s)d\nu(t)\le0.
\]

On Liu's tangent subspace annihilating \(1,s,s(1-s)\), the entropy and \(R\) quadratic forms agree. The unprojected theorem is stronger than that restricted equality and implies Liu's Section V-A Hypothesis 1. The proof uses an exact Taylor remainder, two rank-one power-series kernels, an explicit Lorentz factorization of \(D_r\), the Hilbert-ball reciprocal Gram series, quantitative corner extensions, nonnegative integration, and finite-signed-measure approximation.

The proof was independently re-derived and reviewed. Its theorem envelope matches the submission manuscript. It is ordinary mathematics, not a formal proof-assistant artifact.

- Two fresh read-only falsification passes—an ephemeral independent Codex
  referee and a separate internal re-derivation—found no blocker or major
  defect.  Four minor presentation findings were repaired: the
  Horn--Johnson locator, an explicit fixed-domain PSD pullback for the axes,
  joint \(r,s,t\) continuity/measurability, and the inventory of infinite Gram
  series.  Preserved reports:
  [`codex_referee_2026-08-26.md`](LIU_H1/verification/logs/codex_referee_2026-08-26.md)
  and
  [`internal_referee_2026-08-26.md`](LIU_H1/verification/logs/internal_referee_2026-08-26.md).

### Liu machine verification and reproduction

- Independent CAS-free checker: 16 exact identities pass using only integer/Fraction arithmetic.
- Canonical isolated run: exact SymPy identities, finite Arb compression upper bounds, and numerical diagnostics completed successfully.
- All four canonical Liu scripts were reproduced in an isolated `uv` environment with CPython 3.14.3 and exact pinned dependencies.
- Exact stdout logs, commands, versions, hashes, and observed runtimes are in [`LIU_H1/REPRODUCIBILITY.md`](LIU_H1/REPRODUCIBILITY.md) and [`LIU_H1/verification/logs/`](LIU_H1/verification/logs/).
- The 12-page LaTeX manuscript builds without undefined references or citations.

### Liu failed/contradictory checks resolved

- “Equivalent to unprojected \(R\preceq0\)” was rejected; the correct direction is restricted equality plus a stronger sufficient theorem.
- Roundoff-scale grid eigenvalues are not exact zero modes.
- The numerical ratio `0.9999999906` is near-tight evidence, not proof of exact operator tightness.
- Finite Arb compressions do not prove the continuum theorem.
- Liu's “codimension two” is correct inside the probability simplex; homogeneous signed tangent directions impose three constraints.
- The printed decimal `0.382709087918741` is conditional and is a rounding of the equation-defined value.

### LIU-H2-PROVED-AND-THE-CONSTANT-UNCONDITIONAL (2026-09-02) - EVERYTHING BELOW IN THIS SECTION IS HISTORY

**PROVED [machine-verified algebra + human-audited chain; two internal adversarial reviews SOUND; external refereeing open].**
Liu Hypothesis 2 holds for the paired class, for all Borel probability
measures, and for every conditionally i.i.d. coupling, at the exact
root-defined constants \((\beta^*,m^*)\): \(\mathrm{numerator}\ge(M/m^*)\,E_\mu[h]\)
with no mean constraint. Consequently Liu's Proposition 3 applies with
\(C=(1-c)/m^*>1\) for every \(c<c'\), both of Liu's hypotheses are bypassed, and
**the union-closed constant \(c'=1-m^*=0.38270908791873502993\ldots\) is
unconditional** (certified enclosure \(\pm4.95\times10^{-46}\); Liu's printed
\(0.382709087918741\) is wrong in its last two digits). Authority: `RESULTS.md`
Fourth result and open-problems row 509;
[`uc/UC_CONSTANT_UNCONDITIONAL_2026-09-02.md`](uc/UC_CONSTANT_UNCONDITIONAL_2026-09-02.md);
certificates `liu9-h2-general-lift.json` (sha256 `eb4874d39904593bddc3f2ae638d1bad6756a19d4ebf6e06a60427b6bb7ac1ce`) and
`liu9-h2-mixture-theorem.json` (sha256 `329f7e2d71af8cd78d1a921c72b9ae4d05b71113134eb3341d69932f19ea24b3`) composing
`liu9-h2-twovar.json` and `liu9-h2-boundary.json`. The block-copositivity,
scalar-margin, tube, and cover framings recorded below were the road; none is
the live obligation, and none should be cited as such.

### LIU-H2-SIZE-BIASED-QONE-THEOREM (2026-08-29) - THE ENTIRE Q=1 FACE IS PROVED; FULL H2 REDUCED TO ONE EXACT BLOCK-COPOSITIVITY INEQUALITY

**The near-maximal local theorem is unchanged.** The q-degenerate seam
certifies \(q_*=1/2254\) and
\[
\boxed{\rho=1/1701},
\]
within \(2.68\times10^{-7}\) of the current seam proof ceiling. Reports:
[`liu9-qdegenerate-maximal.json`](uc/verification/results/liu9-qdegenerate-maximal.json),
`sha256 353131d9124516841fe32da2a8f09829d401f9304890a604d91ff544ca693107`,
and the composed full-half-space tube
`sha256 4651ad293d5205a9ba416a3474fbc0b4ee6e553c762d4790f2ce7ccc83612742`.

**PROVED [the exact endpoint blocker].** Mean feasibility contracts
\[
 a_1\in[13/32,7/16],\ a_2\in[0,1/32],\
 b_1\in[1/32,1/16],\ b_3\in[13/32,7/16],\ b_5\in[15/16,1]
\]
to \(b_1>0.0579468\) and \(b_5>0.9968846\). Combining every
\(s\log s\) singular term before interval evaluation and using a nonpositive
mean multiplier proves the **raw gap**
\[
\mathrm{gap}\ge0.01232508028745.
\]
The uncontracted mutation remains negative, so the exact mean cut is
load-bearing. Report
[`liu9-endpoint-support.json`](uc/verification/results/liu9-endpoint-support.json),
`sha256 6510dc6c8b55b0132ee6993cee7f0141a11046013c1de5e0f06b51b1b71df843`.
An implementation-independent reconstruction proves the weaker positive lower
\(0.0120289787046\) and the exact inward-q identity to
\(1.57\times10^{-102}\); report
`sha256 7b97e2d4e3936c248ef70f76431fbcab082186cf781a810341c9e59317221837`.

**PROVED [universal q=1 theorem, stronger than the 230-box quotient].**
For any probability law \(P\) on \([0,1]\) with mean \(M\), put
\(Q(dx)=xP(dx)/M\). Exact entropy algebra gives
\[
F(P)=M\,\mathbb E_{Q\otimes Q} C_M(X,Y),
\]
where
\[
C_M(x,y)=\frac{2M-1}{2}\left(\frac{h(x)}x+\frac{h(y)}y\right)
-M\{\mu(x)+\mu(y)-\mu(xy)\}
+\beta M\frac{h(\pi(x,y))-h(xy)}{xy}.
\]
The new certificate proves \(C_m\ge0\) pointwise on the whole square:
524,800 symmetric bulk boxes, an exact small-support inequality, a
positive-definite Hessian box around \((x^*,x^*)\), geometric endpoint shells,
and finite symmetric/asymmetric tail inequalities. Moreover
\(\partial_M C_M=[(1-\beta)h(xy)+\beta h(\pi(x,y))]/(xy)\ge0\).
Consequently \(F(P)\ge0\) for **every** law of mean \(M\ge m=px\), not merely
three-atom laws. Report
[`liu9-size-biased-qone-kernel.json`](uc/verification/results/liu9-size-biased-qone-kernel.json),
`sha256 eb645792526578d115da8356f7d9378812164273900a60607691a8f32661f19a`;
an independent read-only referee found no blocker after the corrected Hessian
radius and complete infinite-tail partition.

**PROVED [finite inward-q algebra].** With \(r=1-q\), the full raw gap is
exactly \(G_0+rG_1+r^2G_2\); entropy arguments do not depend on \(q\), so there
is no asymptotic remainder. Tests reject a missing linear term, decorrelated
\(r,r^2\), an inward-rounded \(1-q\), a positive feasibility multiplier, and
unsafe zero-support/zero-mass reductions. The current endpoint suite has
19/19 passing mutations; the q-one kernel suite has 4/4.

**OPEN [one exact analytical inequality; no residual-count claim].** For
\(\alpha_i=xP_i\), the full two-component gap has the exact block form
\[
G=\langle\alpha_0^{\otimes2},B_{00}\rangle
+\langle\alpha_1^{\otimes2},B_{11}\rangle
+2\langle\alpha_0\otimes\alpha_1,B_{01}\rangle,
\]
\[
\begin{aligned}
B_{00}&=\frac{(1-q)^2}{M}C_M+\beta q(1-q)k_\pi,\\
B_{11}&=\frac{q^2}{M}C_M+\beta q(1-q)k_\pi,\\
B_{01}&=\frac{q(1-q)}{M}C_M-\beta q(1-q)k_\pi.
\end{aligned}
\]
Thus unconditional H2 is now equivalent to copositivity of this \(2\times2\)
block kernel on the positive paired measures with shared masses. Neither
\(k_\pi\) PSD nor rank-one square completion is valid. The tempting
\(m\)-frozen sufficient inequality is explicitly false for the
**mean-feasible** family \(P_0=\delta_x,P_1=\delta_1,q>0.7237\), but its true
raw gap remains positive; this is a counterexample to that proof shortcut, not
to H2. Extensive exact-feasibility searches still found no negative raw gap.
Exact reduction and the positive-raw-gap shortcut discriminator:
[`liu9-block-kernel.json`](uc/verification/results/liu9-block-kernel.json),
`sha256 58de58843487ffbeb93162698c7ff54f3e27cbb833ecf6b46bc6db171a77cfe2`.
The previous endpoint blocker and the missing finite \(1-q\) remainder are
therefore resolved, and the block-copositivity framing is **superseded
(2026-08-29)**: the displayed 2x2 block kernel was refuted as an identity of
the raw gap, and the exact structure is the scalar-margin diagonalization
\(\mathrm{gap}=F(P_{\rm mix})+\beta q(1-q)\,E_{h\circ\pi}(P_0-P_1)\) with the
kernel \(h\circ\pi\) (rank-2 PSD as \(xy+x(1-x)y(1-y)\), but no \(k_\pi\)
bilinear enters the gap). Under this reduction: Family A (\(P_0=\delta_x\),
\(P_1=\delta_1\)) is **PROVED nonnegative on all of \((0,1]\)** for every
mean-feasible \(q\) (Arb certificates, worst margin \(+1.05\times10^{-3}\),
small-x strip \(\ge 0.02\,x\log(1/x)\)), and the two-atom class has a
**certified curve-II pinch cover** (interval hulls, 12,712 cells, 0 negative,
worst \(+3.4\times10^{-3}\)).  **Correction 2026-08-30:** curve-I is
point-level quadratic extraction, NOT a cover, so the two-atom class is only
partially certified; and the S_asym-identification, 4-slot \(S_4\)
"certification", and the \(r=3\) "MACHINE-VERIFIED" pinch band are
**RETRACTED** (static-dict verdict; invented pencil scanned at 400 random
floats; a 1,568-cell "Arb cover" that is a cell-centre scan with a hardcoded
\(\pm10^{-10}\) envelope behind dead interval code with boundary cells
uncovered; and a single fixed support geometry with an unused cell variable).
See PROGRESS.md `RETRACTION-H2-SASYM-CHAIN`.  Consequently the general-r
composition chain is unsupported and the insertion calculus is OPEN.  Genuine
artifacts: `uc/liu9_block_copositive.py`, `uc/liu9_paired_class_c.py`,
`uc/liu9_scalar_margin_region1{,b}.py`,
`uc/liu9_scalar_margin_region2_curve2.py`.

**RETRACTED 2026-09-01 [box cover].** The `H2-BOX-COVER` artifact
`uc/liu9_h2_box_cover.py` -> `liu9-h2-box-cover.json`, which claimed the
11,025-cell grid became a cover of a positive-measure set (all 9,747
non-vacuous cells clearing at half-width \(2^{-14}\), covered volume
\(3.4986\times10^{-24}\)), is **withdrawn in full**. Its per-cell
certificate is not a lower bound: `radius_certificate` and
`SOUNDNESS_ARGUMENT` step (vi) assert
\(\min_P q_{A(s,q)}\ge q_{A(s,q)}(a^*)-R\) while justifying it with "the
fiber quadratic's min over \(P\) is \(\le\) its value at the \(P\)-point",
an inverted inequality; \(R\) is a single `evaluate_arb` bracket at the
fixed exact-rational masses \(a^*\), so it bounds the gap's oscillation at
one mass point and never the movement of the fiber minimizer in the mass
directions. The module's own sound alternative,
\(c_{0,\rm lo}-\sum_j\mathrm{width}(\mathrm{band}_j)\), is **vacuous** at
every affordable width (\(21651.78\) at \(2^{-10}\), \(1353.19\) at
\(2^{-14}\), \(1.3215\) at \(2^{-24}\), against a \(1.6\times10^{-3}\)
margin) because the \(1/400\)-spaced Cramér extraction amplifies every
interval width by \({\sim}400^{2}\); the sound form would need
\(h\approx2^{-34}\). The certified set for the paired \(r=3\) scalar margin
therefore returns to **measure zero**: the exact-fiber sweep
`liu9-h2-qp-sweep.json` (9,747 non-vacuous cells, worst
\(+1.6305\times10^{-3}\), exact rational supports and \(q\)) is unaffected
and remains the only valid grid authority. The module now refuses to run
without `--retracted-replay`. See PROGRESS.md `RETRACTION-H2-BOX-COVER`.

**COMPUTATIONAL-EVIDENCE [structure that survives].** Two facts measured
while retracting, both float64 Bernstein on exact rational triangulations
and therefore DISCOVERY, not certificates. (i) **Support monotonicity is
not uniform, but its existential form is.** An earlier version of this
paragraph claimed all six partials weakly nonnegative on all 12,383 exact
mass-polygon triangles; that was a bug of mine — I had used
\(\partial_x\pi=y(2-2x-y+xy)\) instead of the correct
\(\partial_x\pi=y[1+(1-2x)(1-y)]=y(2-y-2x+2xy)\) — caught by an
independent Arb port in a subagent and adjudicated against
`liu9_boundary_layer.closed_form_partial` and a 25-digit finite difference
of `gap_mp`, which agree on a negative value
\(-0.117029399303\) at cell `g4|q8`, geometry
\((1/8,1/4,3/8\,|\,1/8,1/4,7/8)\), \(q=1\), masses \((0,0,1)\). The
corrected census over the same 12,383 triangles is: \(k=0,3\) — 10,878
weakly nonnegative, 0 undecided; \(k=1,4\) — 9,640 nonnegative, 1,238
undecided; \(k=2,5\) — 7,758 nonnegative, 3,120 undecided; 1,505
degenerate-zero per direction (the \(q\in\{0,1\}\) columns). **No
direction is ever certified weakly nonpositive.** All six simultaneously
definite on \(7{,}058/12{,}383=56.997\%\). The \(q\)-partial census is
independent of the bug (computed from gap values) and stands at
\(86.58\%\). Crucially, the **existential** statement — there exists \(k\) with
\(w_k>0\) and Bernstein lower bound of
\(\partial\,\mathrm{gap}/\partial b_k\ge0\) — is now **MACHINE-VERIFIED in
Arb** (\(\texttt{ctx.prec}\ge320\), exact point supports and \(q\), no
boxes) on \(12{,}383/12{,}383\) triangles with \(\texttt{fail}=0\)
(`liu9-h2-support-envelope.json`, digest `60c5ce99...`, independently
re-run by me: exit 0, 54.1 s). That Arb census also settles the universal
question negatively: \(\texttt{nonpos}=0\) in every direction of both the
partial and the weight-cleared census (nothing is ever certified
decreasing), yet no direction is sign-definite across the whole grid, so
universal support monotonicity is FALSE; the worst nondegenerate margin is
exactly \(0\) (cell `g2|q8`, \(k=3\)). (ii) The mean
constraint is **active** at the minimum: of the 9,747 non-vacuous cells,
9,670 (\(99.21\%\)) attain their certified minimum on the mean face and
all 9,747 attain it at a polygon vertex; a nine-variable projected
continuum search returns minimum gap \(+1.3536588647\times10^{-3}\) at mean
slack \(2.15\times10^{-6}\). Together these give an **activation lemma** in
its correct existential, value-level form: the feasible set is compact and
the gap continuous, so a minimizer exists; from any minimizer with
\(M>m\), lower a support with \(w_k>0\) and
\(\partial\,\mathrm{gap}/\partial b_k\ge0\) (one exists at every point by
(i)); masses and \(q\) are untouched, so feasibility holds, the gap is
nonincreasing, and \(M\) decreases — driving all supports to 0 would give
\(M=0<m\), so the path meets \(\{M=m\}\). Hence
\(\min\{M\ge m\}=\min\{M=m\}\); no strictness and no uniform monotonicity
are needed. On that surface \(M\) is affine in \(q\), so for
\(\Delta=M_1-M_0\ne0\) the mean equation pins \(q^\*=(m-M_0)/\Delta\) and
clearing \(\Delta^2\) gives
\(\Delta^2G(q^\*)=u^2F_0+v^2F_1+uvT\) with \(u=M_1-m\), \(v=m-M_0\ge0\) —
an exact elimination of the \(q\) direction (identity independently
re-verified: \(\max|\Delta^2G(q^\*)-\Phi|=5.72\times10^{-17}\) over 400 random
admissible draws). The correct mean-aware bound is
\(T\ge-(u/v)F_0-(v/u)F_1\) at the single ratio fixed by the endpoint means;
the all-\(q\) Cauchy bound \(T\ge-2\sqrt{F_0F_1}\) is **false** on this
domain, refuted exactly by \(P_0=\delta_{3/4}\), \(P_1=\delta_0\)
(\(F_0=0.1218420986\), \(F_1=0\), \(T=-0.4949046680\), negative only outside
the feasible window \(q\le1-4m/3=0.176945\)). In size-biased form
\(\Phi=(u^2I_{00}+2uvI_{01}+v^2I_{11})/m+uvc\) with
\(c=\beta(B_{00}+B_{11}-2B_{01})\), so the sole unresolved absorption is the
sign-indefinite \(uvc\); at fixed supports
\(\Phi=a^{T}(u^2L_x+uvR+v^2L_y)a\) is a homogeneous ternary **quartic** in
the masses, the target for a rational degree-4 Handelman/Bernstein
copositivity certificate — not joint SOS, since \(h(x_ix_j)\) and
\(h(\pi(x_i,x_j))\) remain transcendental in the supports. On the
\(\Delta=0\) stratum the mean is \(q\)-free and the four exact cases reduce,
via the universal size-biased theorem, to the single new obligation
\(4F_0F_1-T^2\ge0\).

**REFUTED 2026-09-01 [radial activation].** The cleaner deformation
\(b(t)=(1-t)b\) needs only the scalar Euler condition
\(\sum_jb_j\,\partial\,\mathrm{gap}/\partial b_j\ge0\). Exact per-triangle
minimization of that Euler sum — it is an exact mass quadratic, so three
vertices, three edge criticals and the interior critical give the true
minimum, fit residual \(5.55\times10^{-16}\) — yields a strictly negative
minimum on \(3{,}925/12{,}383\) triangles, worst \(-0.10240072439008191\) at
cell `g4|q8`, geometry \((1/8,1/4,3/8\,\vert\,1/8,1/4,7/8)\), \(q=1\). So
radial shrinking can strictly increase the gap. Only a greedy
one-coordinate path survives, and its honest form is a differential
inclusion \(db/ds=-e_k/w_k\) (giving \(dM/ds=-1\), so \(s\) runs over the
finite interval \([0,M-m]\), and \(d\,\mathrm{gap}/ds\le0\) a.e.) under the
active existential hypothesis: at every feasible point with \(M>m\) there is
\(k\) with \(w_k>0\), \(b_k>0\), \(\partial\,\mathrm{gap}/\partial b_k\ge0\).
The reduction \(\min\{M\ge m\}=\min\{M=m\}\) is **CONDITIONAL** on that
hypothesis plus existence/continuation of the flow; the grid census neither
filters \(k\) by the active conditions nor covers off-grid points, and the
path leaves the \(k/8\) grid immediately. Not PROVED. Recorded
as refuted: Lagrange dualization of the mean constraint (nonconvex duality
gap, Lagrangian minimum \(-0.0210\) at mass vertex \((0,0)\) with
\(\lambda=0.9979\) where the primal cell minimum is
\(+1.63\times10^{-3}\)); and brute-force boxing, since covering
\([0,1]^7\) at half-width \(2^{-w}\) from the \(1/8\) grid costs
\((2^{w}/2^{4})^{7}\) boxes per macro cell — \({\sim}5.6\times10^{14}\) at
\(w=11\), and \({\sim}4.4\times10^{12}\) even on the 6-dimensional surface
\(\{M=m\}\).

### LIU-H2-SECOND-ORDER-LOCAL-TUBE (2026-08-29) - THE LOCAL NINE-VARIABLE GAP IS CLOSED; THE GLOBAL COMPLEMENT IS NOW THE EXACT BLOCKER

**PROVED [exact expansion].** For the simultaneous mean-preserving insertion
\((p-\varepsilon\bar y/x,\varepsilon,q,x,y_1,0,x,y_2,0)\), the pencil is
**exactly**
\[\varepsilon\bar L+\varepsilon^2P_2.\]
Epsilon changes masses only, so every entropy argument is independent of it:
there is no \(\varepsilon^2\log(1/\varepsilon)\), higher term, or remainder.
Automatic polynomial differentiation through the objective SSOT agrees with
the independently derived hand form; worst Arb residual
\(2.813e-68\). The independent verifier
reconstructed 36 configurations at five epsilon scales with maximum residual
\(3.285\times10^{-111}\).

**PROVED [second order at the double root].** \(S(y)=2L(y)+Q(y)\ge0\) on
\([0,1]\), over \(1,666\) exactly abutting cells.
The zero layer has \(Q/y^2\ge0.091066\);
near \(x\), \(S''\ge1.406691\).
\(S(0)=S(x)=S'(x)=0\) are exact structural/optimizer identities; Arb centre
width is dependency error and is not used as a deficit.

**PROVED [simultaneous third-atom interaction].** The distance has the
load-bearing quadratic term
\(\varepsilon^2q(1-q)(y_1-y_2)^2\); dropping it is a rejected mutation.
The interaction is genuinely negative on part of the domain, so first-order
decoupling cannot be reused. A divided-difference Arb cover plus exact
minimization of the quadratic in \(q\) proves \(S_{\rm asym}\ge0\) on
\([0,1]^2\times[0,1]\): \(1,388,611\)
cell pairs, of which \(544,770\) need
the interior \(q\)-vertex. Hence the exact pencil is nonnegative for every
\(0\le\varepsilon\le1/2\).

**REFUTED [sharp endpoint, pencil only].** Positivity does *not* hold on the
whole feasible split segment. At \(y_1=y_2=1\),
\(\varepsilon=px\), the raw gap is structurally zero because numerator and
\(E H(X)\) both vanish on the \(\{0,1\}\) law, while the pencil is at most
\(-0.007245231\). This is not a
counterexample to Hypothesis 2; it is the exact geometric ceiling on this
pencil route.

**PROVED [strict mean-feasible normal -- closes the half-space domain].**
Writing the mean excess as \(\delta\ge0\) moves mass from support \(0\) to
\(x\) in both components. Exact bilinearity gives
\[
 P(\varepsilon,\delta)=P(\varepsilon,0)
 +\delta[h(x)+2\varepsilon\bar R]
 +\delta^2[K(x,x)-\kappa x^2].
\]
\(\bar R\) is a convex combination of
\(R(y)=K(x,y)-(y/x)K(x,x)\). Arb proves the worst linear normal coefficient
is \(\ge0.235013\) on all
1,666 cells and the quadratic coefficient is
\(\ge0.633528\). The
\(q\)-weighted component-mean cross terms cancel exactly. Thus the chart,
insertion, and layer pieces cover **all** \(\mathrm{mean}\ge px\) tube
points, not merely the active boundary.

**PROVED [piecewise local theorem].** Atom relabelling designates the smaller
fragment of every two-way split as the inserted atom, so
\(\varepsilon\le\tfrac12\) without loss modulo the already proved
zero-mass, coincident-support, permutation, and inactive-law gauges. The
second-order result therefore composes with the zero-support, mirror,
interior-chart and q-degenerate certificates. The certified local tube is
\[\boxed{\rho=1/1701},\]
with the q-degenerate seam as bottleneck. All five component digests replay;
piecewise report `sha256 4651ad293d5205a9ba416a3474fbc0b4ee6e553c762d4790f2ce7ccc83612742`.

**OPEN [exact global blocker -- the prompt's claimed implication was missing
this step].** A positive local tube does **not** prove Liu's *global* minimizer
hypothesis. The deterministic complement run at the conservative smaller-tube radius \(\rho=1/1728\)
stopped on its 100,001-box budget with \(49,767\)
residual boxes: \(49,718\) wholly outside
the tube, \(49\) straddling,
\(48,996\) pinned at \(q=1\), and zero
objective-cleared boxes. Saved gap lower bounds range from
\(-4.158883\) to
\(-0.077149\). An exact wholly-outside retained cell dominated by the \(q=1\) frontier is
\[
 a_1\in[0,1/2],\quad a_2\in[1/2,1],\quad q\in[1/2,1],
\]
\[
 b_0,b_2,b_4,b_1,b_3\in[1/4,1/2],\quad b_5\in[1/2,1].
\]
Its certified distance-squared range begins at \(1.3524936813\times10^{-4}\),
strictly outside \(\rho^2=1/4096^2\), while the retained direct enclosure is
only gap \(\ge-2.062290692\), objective \(\ge0.008247823\). That exact cell,
not a sampled point, names where the current global proof fails.

This is a deterministic frontier,
not a proof failure or a counterexample: report
[`liu9-complement-residual-rho1-4096-100k.json`](uc/verification/results/liu9-complement-residual-rho1-4096-100k.json),
`sha256 1c5f1b29368c22ec3cd215fab058e8a10436113a2a111629694ec394a71613ec`.

**COMPUTATIONAL EVIDENCE [not proof].** Four constrained 9-D differential-
evolution runs, \(177,871\) evaluations total, found no
negative raw gap. They converged to a positive constant-support stationary
family and did not reliably locate every boundary zero, so they cannot promote
Hypothesis 2. Search report `sha256 847e09422e55765f0e8cc1b777e454bbada80263ef794e0dd3762fdc86cb864a`.

**NOVEL [absence after stated search].** Liu's paper and all ten checked citing
works contain no second-order transverse expansion, simultaneous-insertion
cross term, local mass theorem, or sharp raw-gap-zero/negative-pencil endpoint.
Primary source and bounded absence audit:
[`liu9-second-order-novelty.json`](uc/verification/results/liu9-second-order-novelty.json),
`sha256 f4785245764c57ce2663c1531a1e9d2a92fe7eb897584307917491e33659b871`.

Artifacts: [`uc/liu9_second_order.py`](uc/liu9_second_order.py), 21 mutation
tests [`test_liu9_second_order.py`](uc/verification/test_liu9_second_order.py),
6 integration tests
[`test_liu9_piecewise_tube.py`](uc/verification/test_liu9_piecewise_tube.py),
independent verifier
[`independent_second_order_check.py`](uc/verification/independent_second_order_check.py),
report [`liu9-second-order.json`](uc/verification/results/liu9-second-order.json)
`sha256 0b1b5bfc11313af7c97224e06b14532fd3da647141bc6b1be7045cde536d2958`.

### LIU-H2-TRANSVERSE-POSITIVITY (2026-08-29) - THE FEASIBLE HALF-SPACE QUESTION IS ANSWERED FOR THE ONE-NEW-SUPPORT FAMILY

**The decisive open question was whether a feasible-half-space nine-variable
theorem holds.** For the transverse family the ambient refutation actually
exploited, it does, and the reason is structural.

**PROVED [structural mechanism].** Liu's eq. (87), \(x^{*2}+x^{*2}(1+\bar
x^{*2})=1\), IS the protocol identity \(\mathrm{prot}(x,x)=1-x^2\) written in
protocol form -- verified against the primary source. Expanding it gives
\(x^4-2x^3+3x^2-1=0\), a polynomial that **never appears in the paper**. What
the paper does not record is the consequence: since \(h\) is symmetric about
\(1/2\), (87) forces \(h(\mathrm{prot}(x,x))=h(x^2)\), hence \(K(x,x)=h(x^2)\)
**for every** \(\beta\), hence \(pK(x,x)=h(x)\) and \(A=h(x)\).

**PROVED.** The mean-preserving transverse coefficient
\(L(y)=2pK(x,y)-h(y)-(y/x)A-\kappa y^2(y-x)^2\) has a **double root** at
\(y=x\). \(L(x)=0\) is structural -- certified as a rational-function identity
in exact `Fraction` with the entropy values as free symbols, so it uses no
property of any constant. \(L'(x)=0\) needs (87) *and* the \(\beta\)-system
(89)-(90); perturbing \(\beta\) by \(10^{-6}\) moves the slope to
\(2.12\times10^{-7}\), so the cancellation is Liu's optimizer and not a
transcription accident. \(L''(x)=1.2289613682\), reproduced independently.

**PROVED, universal enclosure.** \(L(y)\ge-4.562\times10^{-136}\) for **every**
\(y\in(0,1)\), replacing the previous 65-point scan. The deficit is the
*square* of the residual of Liu's own equations. \(L\) is **not** globally
convex (\(L''<0\) on roughly \([0.01,0.43]\)), so the proof splits \([0,1]\) into
five exactly abutting pieces: a \(y\log(1/y)\) boundary layer where
\(\mu(u)=\frac{1-u}{u}\log\frac1{1-u}\in[1-u,1]\) by termwise series
comparison (no \(0/0\) is ever evaluated) plus one interval tail cell containing
\(y=0\); \(11{,}777\) bulk cells, weakest margin \(7.32\times10^{-3}\); and a
\(4{,}096\)-cell double-root window with \(L''\ge0.4975\).

**PROVED, and free.** The **asymmetric** family decouples exactly:
\(L_{\text{asym}}(y_1,y_2,q)=(1-q)L(y_1)+qL(y_2)\), because the objective's
mean is the \(q\)-weighted mixture (so \(y_1\ne y_2\) is feasible with
correction \(\bar y\)) and every term is linear or bilinear in the measure. So
the bound extends to all of \((0,1)^2\times[0,1]\) with no extra work.

**NOVEL [absence after stated search].** The entropy-symmetry collapse, the
\(\beta\)-independence of the kernel diagonal, and any transverse or
first-variation analysis of Hypothesis 2 appear in neither Liu's paper nor any
of its ten citing works; the root is uncatalogued in OEIS. Hypothesis 2 itself
remains open in the literature, which is visibly confused about it -- three
different constants (0.38234 / 0.38237 / 0.38271) are attributed to Liu by
different papers.

**OPEN.** This is a **first-order** statement about a family with at most one
new support point per component. Simultaneous third-atom insertions and the \(\varepsilon^2\) terms were
**OPEN in this entry and are superseded by `LIU-H2-SECOND-ORDER-LOCAL-TUBE`
above**. The local tube is now proved; global Hypothesis 2 remains open only
because its complement has not been certified.

Module [`uc/liu9_transverse.py`](uc/liu9_transverse.py), 35 mutation tests
[`test_liu9_transverse.py`](uc/verification/test_liu9_transverse.py),
independent derivation
[`independent_transverse_check.py`](uc/verification/independent_transverse_check.py),
report [`liu9-transverse.json`](uc/verification/results/liu9-transverse.json)
`sha256 153fb77bbe2ae4eff027c9d3d80c9e1d905c4480db288be386649c21671534f8`.

### Liu originality and open scope

A theorem-level search through 2026-08-25 found no earlier proof of Liu's scale-one kernel hypothesis or the stronger unprojected residual theorem. The closest predecessors are AHS's unperturbed product-entropy kernel and Liu's sufficiently-small-scale perturbation lemma. Classical ingredients such as Schur closure and the Hilbert-ball kernel are not claimed as new. The problem-specific Taylor regrouping and Lorentz factorization are plausibly new, but universal priority remains **OPEN**.

An arXiv `all:"union-closed"` feed refresh at 2026-08-26T15:10:15Z still
reported 103 records and no item posted after the full-audit cutoff; no new
Liu-H1 collision appeared.  This one-index refresh does not broaden the
universal-priority claim.

**Superseded 2026-09-02:** Liu's Section V-B hypothesis is PROVED and the constant \(c'=1-m^*=0.38270908791873502993\ldots\) is unconditional (see LIU-H2-PROVED-AND-THE-CONSTANT-UNCONDITIONAL above and `RESULTS.md`, Fourth result); the printed `0.382709087918741` is wrong in its last two digits. The paragraph that stood here read: "Liu's Section V-B global-minimizer hypothesis remains OPEN. Therefore `0.382709087918741` remains conditional and is not an unconditional union-closed bound from this work." — historical.

The 2026-08-27 measurement sharpens where that gap sits. At tube radius
\(\rho=0.1\) the complement branch-and-bound leaves 71, 76, and 79 residual
boxes after 26,247, 866,551, and 3,000,001 processed boxes: a 114-fold budget
increase moves the frontier by eight boxes, so subdivision does not close the
complement and the earlier \(2.7\times10^4\)-box extrapolation is refuted as a
completion estimate. Of the 79 survivors, 70 touch the entropy-singular
support corner \(b_j=0\) and 57 of those also pin \(q\) at a degenerate
endpoint \(\{0,1\}\). An independent reclassification on 2026-08-28
([`uc/liu9_survivors.py`](uc/liu9_survivors.py)) corrects an earlier count in
this file: **61** boxes, not nine, are coordinate-identical across the two
localized runs; nine is the number of *mirror-only* survivors. The obstruction
is exactly the two strata the conjectured piecewise local lemma addresses, not
a diffuse volume. Reports:
[`liu9-complement-residual-rho0.1.json`](uc/verification/results/liu9-complement-residual-rho0.1.json),
[`liu9-complement-residual-rho0.1-3M.json`](uc/verification/results/liu9-complement-residual-rho0.1-3M.json),
[`liu9-survivor-classification.json`](uc/verification/results/liu9-survivor-classification.json).

**MACHINE VERIFIED: the residual frontier is exactly the union of the two
boundary strata.** Every one of the 79 survivors is covered either by the
zero-face hypothesis \(2(1-\beta)\mu-1\ge0.11105875229\ldots>0\) (70 boxes)
or is mirror-only and entirely outside the tube (9 boxes). The nine are
`large-044`, `-046` through `-052`, and `-054`.

**PROVED: ingredient (ii), the zero-support layer.**
([`uc/liu9_boundary_layer.py`](uc/liu9_boundary_layer.py).) Differentiating
Liu's own transcription gives the exact identity
\(\partial\,\mathrm{gap}/\partial b_j=w_j G_j(V)\), and for every point of the
raw box with \(0<b_j\le y_0\),
\[G_j(V)\ \ge\ \bigl[2(1-\beta)\mu(V)-1\bigr]\log(1/b_j)-2(1-\beta)\log\tfrac1{1-y_0},\]
using only the mean constraint. With \(y_0=1/32\),
\(\Lambda=2(1-\beta)(\text{mean}-y_0)-1\in[0.0548120372\ldots]\) and
\(m_0=\Lambda(\log(1/y_0)+1)-K\in[0.1876317632\ldots]\), every mean-feasible
point satisfies
\[\mathrm{gap}(V)\ \ge\ \mathrm{gap}(V|_{\text{layer}:=0})+m_0\sum_{j\in S}w_jb_j,\]
so the layer is *reduced to the face* \(b_j=0\), where \(|h'''|\le1022.94\)
is finite for every surviving support in \([y_0,1-y_0]\).

**PROVED: ingredient (iv), the mirror stratum \(b_j\to1\).**
([`uc/liu9_mirror_layer.py`](uc/liu9_mirror_layer.py).) The face is \(C^3\)
only on \([y_0,1-y_0]\): \(h'''\) blows up at 1 as it does at 0, and a
small-mass atom at support 1 costs only \(g(1)=0.0956\) per unit weight, so it
is well inside the tube.  That value is motivation only: every bound divides by
the *infimum* of \(g\) over the whole stratum, \(g(1-t_0)\), which at
\(t_0=1/32\) is \(24.16\%\) smaller than \(g(1)\) — the unsquared cost
\(b(b-x)\) drops \(12.92\%\) and \(g\) is its square — and gives critical
radius \(0.17951\) rather than the \(0.20613\) the endpoint would suggest, a
\(14.83\%\) overclaim.  `certify_face_infimum` certifies that
the infimum sits at the inner edge, and a mutation test rejects the
substitution. The singular coefficient is
\(1-2(1-\beta)W_1-2\beta A_1\), matched to the closed form to \(10^{-40}\)
at six configurations including the discriminating \(W_1=\tfrac12,A_1=1\)
case where it equals \(-\beta\). The mean does not sign it, but the *tube*
does: \(W_1\le\rho^2/g(1-t_0)\) with \(g(b)=[b(b-x)]^2\), so
\(M=1-2\beta-2(1-\beta)\rho^2/g(1-t_0)>0\) below an explicit critical radius.
A weak-duality bound \(\Theta\) controls the non-mirror cross sum, giving
\(G_j\ge M\log(1/t)-K_1\) and the reduction
\(\mathrm{gap}(V|_{\text{mirror}:=1})-\mathrm{gap}(V)\ge m_1\sum_jw_jt_j\).
The reduction raises supports *to* 1, so the mean rises and feasibility is
preserved for free. Certified at \((\rho,t_0)=(1/10,1/64)\) with
\(\kappa\le0.393693\) and **REFUTED** at \((1/10,1/32)\), so the threshold
is sharp. Report:
[`liu9-mirror-layer.json`](uc/verification/results/liu9-mirror-layer.json).

**PROVED, and it corrects both of the above: the two \(\kappa\) columns were
being compared against the wrong constant.** \(\mathrm{gap}=EHX\cdot(\Phi-1)\)
is an exact identity wherever \(EHX>0\), and
\(H_*=EHX(P_*)=p\,h(x)=0.5526667300\ldots\). The ceiling
\(0.7004675091\ldots\) printed by [`uc/liu9_tube.py`](uc/liu9_tube.py) is in
\((\Phi-1)/\mathrm{dist}^2\) units, while every layer \(\kappa\) derived from
\(\partial\,\mathrm{gap}/\partial b_j\) is in raw-gap units; the raw-gap
ceiling is \(0.3871250877\ldots\). Both modules now convert explicitly. The
conclusions survive with a larger margin.

**PROVED, and it removes ingredient (iii) as a second-order obstruction: the
\(q(1-q)\) factors cancel exactly.** The objective's quadratic form is
\(A|s|^2+Cq(1-q)|d|^2\) and the squared distance is
\(|s|^2+q(1-q)|d|^2\), so the correct comparison is positive semidefiniteness
of the pencil \(\mathrm{Hess\,gap}-\kappa\,\mathrm{Hess\,dist}^2\), in which
\(q(1-q)\) cancels for every \(0<q<1\); at \(q\in\{0,1\}\) the \(d\)
direction is a joint gauge. Comparing \(\lambda_{\min}\) against
\(\lambda_{\max}\) instead reports a spurious collapse. Audit:
[`uc/LIU_H2_INGREDIENT_I.md`](uc/LIU_H2_INGREDIENT_I.md).

**PROVED on the chart, with the obstruction to enlarging it quantified:
ingredient (i).** ([`uc/liu9_smooth_chart.py`](uc/liu9_smooth_chart.py).) On
the exact active-mean chart both \(\mathrm{gap}\) and \(\mathrm{dist}^2\)
vanish to second order at the centre, so Taylor's mean-value form gives
\(f(v)=\tfrac12v^{\mathsf T}\mathrm{Hess}f(\xi)v\) exactly and **no third
derivative is ever required** -- which is what defeated the one-piece cubic
bound. The jet Hessian reproduces \(H_*\mathrm{diag}(A,Cq(1-q))\) and
\(\mathrm{diag}(2px^2,2(p^2+px^2)q(1-q))\) to \(1.9\times10^{-67}\), an
independent confirmation of the published curvatures. Certified:
\(\mathrm{gap}\ge\kappa\,\mathrm{dist}^2\) on
\(|s-x|\le1/2048,\ |d|\le1/2048\), for every \(q\in[1/4,3/4]\), with
\(\kappa\ge234627/1048576=0.22376\ldots\), which is \(57.8\%\) of the
q-free ceiling. That is exactly the explicit radius, one-sided uniform bound,
and explicit raw-gap \(\kappa\) the audit named as missing. **OPEN, with the
blocker localized to arithmetic rather than geometry:** naive interval
enclosure of the Hessian inflates linearly in the box radius with constant
\(\approx811\) (width \(14.36\) at radius \(1/64\) falling to \(0.396\)
at \(1/2048\)), so reaching the \(|s-x|\approx0.15\) that a
\(\rho=1/10\) tube needs would take \(\approx9.7\times10^7\) cells. The
missing tool is Taylor-model or centered-form arithmetic, or an analytic
third-derivative bound on the chart. Report:
[`liu9-smooth-chart.json`](uc/verification/results/liu9-smooth-chart.json).

**PROVED on the pure-d endpoint chart: ingredient (iii).**
([`uc/liu9_qdegenerate.py`](uc/liu9_qdegenerate.py).) The inner core
\(|d|\le1/32\) has \(q\)-uniform \(\kappa=1/3\); the endpoint annulus
\(1/32\le|d|\le1/4\) with \(\min(q,1-q)\le1/4096\) has \(\kappa=1/20\),
from \(H_*D(Q_d)\ge\tfrac14\delta(Q_d)^2\) over 384 exact dyadic cells,
normalized \(D\ge0.00035275034\). The **seam is PROVED** under
\(\rho^2\le p^2\varepsilon_{\rm sm}^2q_*(1-q_*)\), giving the historical \(\rho=1/4096\), now superseded by
\(\rho=1/1701\), so no pure-\(d\) tube point falls between the two \(q\) regimes. The
transcription is symmetric only under the *simultaneous* swap
\((q,P_0,P_1)\mapsto(1-q,P_1,P_0)\); \(q\mapsto1-q\) alone is **REFUTED**,
changing the gap by \(2.155\times10^{-2}\).

**The centered form cut the ingredient (i) obstruction by an order of
magnitude.** ([`uc/liu9_chart_centered.py`](uc/liu9_chart_centered.py).) A
third-order interval jet, with \(h'''(u)=(1-2u)/(u^2(1-u)^2)\) — whose sign
was caught by the module's own 20-digit Richardson gate rather than assumed —
bounds each Hessian entry by its exact centre value plus the third-derivative
enclosure times the box radius. Legitimate here precisely because the two
boundary strata are separately closed, so every chart support stays inside
\([0.659,0.722]\) and \(h'''\) is bounded. The inflation constant falls from
\(811.29\) to \(72.20\), the certified radius rises \(8\times\) to
\(1/256\) with \(\kappa\ge4119063/33554432=0.1227576\), and the covering
estimate for \(|s-x|\approx0.15\) falls from \(9.7\times10^7\) to
\(1.9\times10^5\) cells.

**The factor of eight is CLOSED, and both premises behind the plan for it were
wrong.** [`uc/liu9_chart_cover.py`](uc/liu9_chart_cover.py) certifies the
pencil \(\mathrm{Hess\,gap}-\kappa\,\mathrm{Hess\,dist}^2\succeq0\) on
\(|s-x|,|d|\le1/32\) for every \(q\in[1/4096,4095/4096]\) at the same
\(\kappa\), over \(393{,}216\) exhaustively abutting cells, weakest
determinant margin \(2.48\times10^{-5}\) at
\((s-x,d,q)=(15/512,-15/512,1/4096)\). Report
[`liu9-chart-cover.json`](uc/verification/results/liu9-chart-cover.json).

First, **no radial quadrature is needed**: a uniform cover suffices. PSD is a
*pointwise* property, so certifying it cell by cell certifies it on the union;
and the box is convex and contains the chart centre, so the mean-value step
draws its segments from geometry once, not per cell. Per-cell segment
containment would matter only if positivity were certified from each cell's own
local expansion, which is not what the pencil route does.

Second, **the radius was never the binding gap — the \(q\) range was.**
Ingredient (iii) names \(q_{\rm int}\le q_*=1/4096\), not \(1/4\), and both
the old certificate and the first cover attempt covered only \(q\in[1/4,3/4]\).
A uniform \(q\) grid stops certifying below \(q=1/64\): \(\mathrm{gap.hdd}\)
and \(\mathrm{dist.hdd}\) both carry an exact factor \(q(1-q)\) — the measured
ratio \(m_{22}/(q(1-q))\) is \(0.69512040\) at
\(q=1/4,1/64,1/128,1/1024,1/4096\) alike, to nine digits — so the true entry
vanishes linearly at the endpoints while a fixed-width cell's enclosure error
does not. Octave cells of relative width \(1/64\) restore it, the same device
[`liu9_boundary_layer._octave_grid`](uc/liu9_boundary_layer.py) already uses on
the entropy side.

Step B inherits an explicit deficit rather than an exact zero, and the bound is
an **all-\(q\) enclosure, not a sampled one**. The chart centre collapses to
\(\{x\!:\!m,\,0\!:\!1-m\}\) in both components, so every \(q\) cancels
structurally; carrying that same collapse as a jet in \(s\) yields the centre
value and \(s\)-gradient as \(q\)-free enclosures,
\(|c_0|\le6.74\times10^{-69}\) and \(|\partial_s|\le2.016\times10^{-68}\).
The \(d\)-gradient is **exactly** zero for every \(q\): at \(d=0\) the two
supports move oppositely while the components carry the opposite weights, so
each term contributes a \(q\)-free multiplier times \((1-q)(-q)+q(1-q)\), whose
expansion has all three coefficients zero in exact `Fraction` arithmetic. The
pointwise Arb residue that scales like \(q(1-q)\) is rounding of the individual
terms, which is why measuring it pointwise made it look \(q\)-dependent.

Hence \(\mathrm{gap}-\kappa\,\mathrm{dist}^2\ge-7.361\times10^{-69}\) on the
box for every \(q\) in the covered range. Liu's abscissa \(x\) is a numerically
determined root of (87)–(90), and that residual is the whole deficit; it is
inherited unchanged from the \(1/256\) certificate, not introduced by the
cover.

**REFUTED: the ambient nine-variable extension.**
[`uc/liu9_ninevar.py`](uc/liu9_ninevar.py) exhibits
\((p-\varepsilon,\varepsilon,q,x,y,0,x,y,0)\) at
\((y,\varepsilon,q)=(1/32,1/1024,1/2)\) with raw gap \(-5.5338\times10^{-4}\)
and pencil \(-5.5348\times10^{-4}\); exact `Fraction` arithmetic puts the
linear coefficient below \(-1/2\) and the finite pencil below \(-1/2048\).
Leading order \(\varepsilon\), with the \(\varepsilon\log(1/\varepsilon)\)
coefficient exactly \(0\). Reproduced independently against `gap_mp`.

**CONDITIONAL, and this is the live gap.** That counterexample is
mean-INFEASIBLE (mean-target \(-6.4408\times10^{-4}\)), so it does not touch
Liu's Hypothesis 2, which lives in \(\{\mathrm{mean}\ge px\}\). The exact
mean-preserving control at the same \(y\) has linear coefficient above
\(1/64\), and across a 65-point \(y\) scan the mean-preserving family is (**superseded 2026-08-29** by the universal enclosure in `uc/liu9_transverse.py`)
positive at **every** point (0 negative) while the ambient split is negative at
46. Its weakest coefficient \(7.865\times10^{-8}\) sits at
\(y=707/1024\approx x\), where the inserted atom merges with the existing one
and the perturbation degenerates. That first-order statement is **superseded** by the exact quadratic theorem
above. A local tube radius now follows; the remaining gap is the global
complement, not the chart-to-nine-variable extension.

### Liu submission package

- Proof audit: [`LIU_H1/AUDIT.md`](LIU_H1/AUDIT.md)
- Manuscript: [`LIU_H1/paper/main.tex`](LIU_H1/paper/main.tex), [`LIU_H1/paper/main.pdf`](LIU_H1/paper/main.pdf)
- Reproduction: [`LIU_H1/REPRODUCIBILITY.md`](LIU_H1/REPRODUCIBILITY.md)
- Independent checker/logs: [`LIU_H1/verification/`](LIU_H1/verification/)
- Literature report and saved sources: [`LIU_H1/LITERATURE_ORIGINALITY.md`](LIU_H1/LITERATURE_ORIGINALITY.md), [`LIU_H1/literature/`](LIU_H1/literature/)

## Submission readiness

| Project | Mathematical readiness | Computational readiness | Literature readiness | Current verdict |
|---|---|---|---|---|
| UC explicit constant | Human audit plus two independent read-only referees, the second finding no proof defect in the entropy bridge or the support reduction; eight minor exposition defects repaired; external specialist/journal review open | Structural, two same-source full replays, and secure source-independent byte-zero replay pass; the entropy bridge is additionally exhausted for every family with \(n\le4\) and sampled to \(n=8\); 343-file/runtime-image evidence pinned | Bounded audit plus 2026-08-26 arXiv refresh, no collision found | Independently reproduced subject to listed trace/Arb/runtime/custody assumptions; strong candidate pending external mathematical review |
| Union-closed constant \(c'=1-m^*\) (unconditional, 2026-09-02) | Machine-verified algebra (two certificates, constants as balls) plus a human-audited information-theoretic chain re-derived by two adversarial readers; external referee open | Both certificates byte-identical across runs and reviewer reproductions; 27,006-configuration independent attack found no negative | Bounded search; no priority claim | PROVED in-repo; strongest bound in this ledger pending external review |
| Liu Hypothesis 1 | Complete human-audited proof; two fresh read-only referee passes found no blocker or major defect; external journal/expert review remains open | Clean isolated reproduction and independent exact checker pass | Bounded audit plus 2026-08-26 arXiv refresh, no collision found | Submission-ready candidate manuscript with explicit non-formal/non-priority qualifications |

## Concrete next actions

1. Ask an external union-closed/entropy specialist to referee the UC support reduction and entropy bridge.  This is the only remaining blocker on the UC theorem and it cannot be discharged inside the repository; two independent model referees have now found no proof defect.
2. Ask a kernel/functional-analysis specialist to referee the Liu Taylor--Lorentz--Gram proof.
3. Submit the Liu-H1 manuscript with the candidate/priority qualifications retained. Liu Hypothesis 2 is now PROVED and the constant \(c'\) unconditional (2026-09-02); H1 is no longer load-bearing for \(c'\), so the H1 manuscript should say so rather than claim the constant.
4. Campaign K is fully committed (16/16 slices COMPLETE with residual = stack
   = `budget_time` = 0, 805,542,368 boxes).  Run the frozen collector to term;
   an exit-0 `COMPOSITE CERTIFICATE` there would raise the certified constant
   to \(0.3822660112501052\).  Nothing weaker may be reported as a certificate.
5. Liu H2 is PROVED (2026-09-02) and the constant \(c'\) is unconditional.
   Next: external refereeing of `uc/UC_CONSTANT_UNCONDITIONAL_2026-09-02.md`
   and the two certificates; then, as DISCOVERY only, whether a three-protocol
   combination (Liu Lemma 8) or another protocol family raises \(c'\).
6. Treat formal proof-assistant verification, hardware-backed execution
   attestation, and removal of the remaining trace/Arb/OS trust assumptions as
   strengthening work; the independent interval implementation itself is
   complete.
