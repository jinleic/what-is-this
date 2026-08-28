# Gate B experiment ledger

All production commands below ran from `math/` on 2026-08-25 with
`nice -n 10` and

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0
```

No two heavy searches ran concurrently.

## E1 — repository and history audit

- Inspected the tracking documents, closure-defect code, prior candidates,
  certificates, paper scope, and campaign layout.
- `git status --short -- math` reported `?? math/`; `git log --all -- math`
  returned no commits. The entire math tree remains untracked, so no historical
  Gate B commit exists to replay.
- Gate B had a plan but no implementation, candidate store, checkpoint, test,
  standalone verifier, or paper section.
- The authoritative displayed supremum had no small-defect restriction.

Status: **EXACT REPOSITORY FACTS**.

## E2 — independent reproduction of the certified n=7 base

Command:

```sh
nice -n 10 ./.venv/bin/python -B uc/shapley_n7_falsifier_cert.py
```

Wall time: 307.75 s.

Reproduced 45 rows, seven counts of 18, incidence 126, 1,600 missing ordered
joins, defect `64/81`, 240 automorphisms, 21 order representatives, 5,152
Bellman states, and 982 clamp-ambiguous states. At every ambiguous state the
unconstrained concave maximum is retained, so the computed ball remains a
rigorous upper certificate.

```text
A_plus <= -0.0286491867944683178160...
A_plus + epsilon_join/50 <= -0.0128467176586658486802...
```

All 5,040 independent float policy replays passed.

Status: **PROVED/CERTIFIED FINITE**.

## E3 — independent reproduction of the n=7 declared-class census

Command:

```sh
nice -n 10 ./.venv/bin/python -B uc/shapley_n7_block_symmetric.py
```

Wall time: 241.60 s.

Reproduced 72/680/2,238 accepted raw families in the three declared classes,
2,980 full-S7 orbits, 16 negative `A_plus` values, four still negative after
the `1/50` correction, and maximum numerical ratio
`0.036259127036737...` at the certified base.

Status: **COMPLETE EXACT** for the declared classes, constraints,
canonicalization, order coverage, and defects; **COMPLETE NUMERICAL** for the
float64 objective.

## E4 — definition and threshold audit

The old n=7 producer computed
`ceil(m*log2(m)/2 - 1e-12)`. An independent integer implementation defined the
threshold as the least `r` with `2^(2r) >= m^m`. The values agreed for every
`3 <= m <= 128`; therefore no n=7 family changed. The producer was moved to
the integer implementation to make its exact evidence label literal.

Normalization was checked against the code: active, pairwise-distinct columns.
It is not part of the displayed Gate B supremum, but the certified base is
normalized. `ALPHA` is the repository-global constant
`356069/10000000`; it has no dimension dependence.

Status: **EXACT**, with a non-result-changing implementation correction.

## E5 — boundedness attempt and tensorization obstruction

The proposed direction was an inequality

```text
-A_plus(F) <= c * epsilon_join(F).
```

Before launching a near-UC search, the behavior under Cartesian product was
checked. Fixed-order conditional laws split by block. A backward Bellman
induction showed that cross-block prefix information cancels from the action
slope, giving exact additivity of both `Q` and `C_plus`, hence of `A_plus`.
Closure success multiplies.

This refutes the boundedness direction: any admissible negative base has
Cartesian powers with numerator linear in the power and denominator at most
one. The certified n=7 family supplies such a normalized base.

Status: **PROVED STRUCTURAL OBSTRUCTION**. A larger finite coefficient search
could no longer answer Gate B.

## E6 — complete n=8 S1 x S7 search

The first launch failed immediately because a script in `uc/gate_b/` did not
put its parent `uc/` directory on `sys.path`. No checkpoint was written. The
import path was fixed and the same bounded command was rerun:

```sh
nice -n 10 ./.venv/bin/python -B \
  uc/gate_b/search_n8_block_symmetric.py --k 1
```

Wall time after the fix: 18.72 s.

Results:

- all `2^16-1 = 65,535` raw masks exhausted;
- 136 cap/Reimer families and 136 canonical S8 orbits;
- 130 normalized families and two negative `A_plus` values;
- eight exact order representatives under the guaranteed S7 action;
- maximum numerical ratio `0.034335490331443...`, at a normalized 70-row
  family with all eight counts 28 and exact defect `283/350`.

The ratio is below the certified n=7 base. Complete append-only records are in
`experiments/n8_k1_checkpoint.jsonl`; the summary is
`candidates/n8_k1_census.json`.
The final resume smoke completed in 0.60 s with `new_evaluations: 0` and
re-emitted byte-identical summary content.

Status: **COMPLETE EXACT / COMPLETE NUMERICAL** for this declared class only.
It is not a complete arbitrary-family n=8 optimum.

The k=2 class was not promoted to a heavy run after the product theorem settled
the universal question. The script remains bounded and resumable if a separate
finite-class catalogue is later desired.

## E7 — standalone theorem certificates

Initial command:

```sh
nice -n 10 ./.venv/bin/python -B uc/gate_b/verify_gate_b.py \
  --write-certificate uc/gate_b/certificates/gate_b_unbounded_arb.json
```

Wall time: 0.66 s. It independently rebuilt the base Arb proof and directly
enumerated the 2,025-row square's combinatorics: 3,920,000 missing pairs and
defect `6272/6561 = 1-(17/81)^2`.

The first adversarial review correctly noted that the initial JSON stored
tensorization only as a formula string. The verifier was strengthened rather
than treating that checkpoint as final. The v2 command preserved the first
artifact and added two actual 2,025-row product Bellman calculations:

```sh
nice -n 10 ./.venv/bin/python -B uc/gate_b/verify_gate_b.py \
  --write-certificate uc/gate_b/certificates/gate_b_unbounded_arb_v2.json
```

Wall time: 35.80 s. The consecutive and alternating 14-coordinate orders agree
with twice the base-order optimum inside Arb residual balls of radius below
`1.4e-70`. The alternating calculation visits 2,985,831 states and has 405,641
ambiguous clamps; retaining the unconstrained maximum guarantees an upper
bound even when the exact boundary case is not resolved.

A second independent proof review returned `VERDICT: ACCEPT` and found no
mathematical blocker. It identified two evidence-labeling issues: v2 placed
identity-derived square `Q`, `C_plus`, and `A_plus` values beside actual square
checks, and prose described the ambiguous-state result as a full enclosure
rather than the rigorous upper certificate actually needed. Both were
corrected. The final command preserved v1/v2 and wrote v3:

```sh
nice -n 10 ./.venv/bin/python -B uc/gate_b/verify_gate_b.py \
  --write-certificate uc/gate_b/certificates/gate_b_unbounded_arb_v3.json
```

Wall time: 36.68 s. V3 separates `actual_square_fixed_orders` from
`derived_from_proved_identities`, removes a tautological residual assertion,
and labels the Bellman result as an upper certificate. It also fixes
`alpha = 356069/10000000` as dimension-independent and machine-checks cap,
incidence, the Reimer witness, and defect arithmetic for powers `k=1,...,8`.
The all-k theorem remains the exact Bellman induction in `PROOF.md`, not an
extrapolation from sampled orders.

Status: **PROVED/CERTIFIED; INDEPENDENT REVIEW ACCEPTED**. V1 and v2 are
preserved intermediate checkpoints; v3 is authoritative.

## E8 — adversarial regression suite

Command:

```sh
nice -n 10 ./.venv/bin/python -B uc/gate_b/test_gate_b.py
```

Initial result: 9 tests passed in 2.08 s wall. After the v3 evidence-label
review and Python 3.9 compatibility change, the final suite passed all 10 tests
in 2.73 s wall.

The tests cover exact Reimer thresholds and constraints, direct square counts,
ordered-pair defect products, Arb `Q` and `C_plus` additivity on a hand example
and four deterministic random products, independent evaluator agreement,
project/candidate $\alpha$ equality, v3 certificate structure, exact partition
of all 40,320 n=8 orders, block canonicalization against a brute S8 orbit, and
the completed n=8 result contract.

Status: **FINAL PASS — 10/10**.

## E9 — literature control

Chase and Lovett's *Approximate union closed conjecture*
(<https://arxiv.org/abs/2211.11689>) uses the same ordered-pair closure-success
notion and proves a frequency statement when almost all joins remain inside.
It does not define the repository-specific `A_plus` Bellman functional or
supply the product repair considered here. Its local small-defect theorem
reinforces the distinction between Gate B's unrestricted supremum and the
separate `epsilon_join -> 0` stability question.

Status: **CITED CONTEXT**, not proof evidence for the Gate B theorem.

## E10 — clean-environment replay

The first clean attempt used the workstation's default Python 3.14. Installing
the frozen `python-flint==0.6.0` fell back to its source distribution and failed
because that release's Cython code still imports removed `PyInt_*` APIs. A
second attempt used the system Python 3.9 wheel, but exposed one verifier
compatibility bug: `int.bit_count()` is unavailable on Python 3.9.

The verifier gained a local `popcount` fallback and the clean run was repeated
without changing any arithmetic dependency:

```sh
/Library/Developer/CommandLineTools/usr/bin/python3 -m venv \
  /tmp/uc-gate-b-py39-20260825
/tmp/uc-gate-b-py39-20260825/bin/python -m pip install \
  'python-flint==0.6.0'
nice -n 10 /tmp/uc-gate-b-py39-20260825/bin/python -B \
  uc/gate_b/verify_gate_b.py
```

The arm64 wheel installed successfully. The initial clean replay took 60.15 s;
after the v3 review corrections, the final clean replay again printed
`PROVED_GATE_B_UNBOUNDED` with the same mathematical intervals and exact counts
in 59.91 s.

Status: **CLEAN-ENVIRONMENT PASS** on Python 3.9 with the frozen dependency.

## E11 — adversarial theorem audit and strengthened v4 certificate

The Gate B statement, original one-sided Bellman source, n=7 producer,
standalone checker, product proof, candidate, certificate, tests, and paper
were re-read from source. The sign convention is unchanged:

```text
A_plus = (1-alpha) Q + alpha C_plus - log2(|F|),
epsilon_join = missing ordered pairs with replacement / |F|^2.
```

An independent exact reconstruction from the seven weight cells reproduced the
stored 45 rows, seven coordinate counts of 18, incidence 126, 1,600 missing
ordered joins, defect `64/81`, success `17/81`, 240 automorphisms, and 21
equal-size order orbits. The original n=7 Arb certificate was rerun before any
new search and reproduced the same substantive output bit-for-bit in 300.34 s:

```text
A_plus upper = -0.0286491867944683178160269789255044875...
required coefficient lower = 0.0362591270367489647359091452025916169...
all 5,040 float policy replays: PASS
```

Independent decimal/rational arithmetic reproduced the n=7 certified ratio,
the n=8 numerical ratio `0.034335490331443001`, and the exact power bounds
`(7k/250)/(1-(17/81)^k)`.

The adversarial audit found and fixed three rigor defects that did not change
the theorem.

1. The displayed supremum had no convention for a negative zero-defect family.
   `DEFINITIONS.md` and the paper now use the extended-nonnegative convention;
   every constructed power has strictly positive defect, so the proof is
   unaffected.
2. The checker verified the observed block-cell set but did not prove that the
   stored rows were the full union of those cells. It now reconstructs the
   family from the cell predicate and requires exact row-list equality.
3. A ratio derived from the one-sided `A_plus` upper certificate must itself be
   labeled only as a lower certificate. V4 exports
   `repair_ratio_lower_certificate` and no unsupported upper endpoint.

The Bellman product argument was expanded algebraically: the four transition
weights sum exactly to one, a common other-block continuation adds the same
constant to every action, and it cancels from
`-V00+V10+V01-V11`. A targeted counterexample test then exhausted all 225
ordered products of nonempty two-coordinate families and every one of their
24 global orders: 5,400 fixed-order `Q` and `C_plus` comparisons, plus exact
defect products. No counterexample was found.

Smaller-dimension controls were rerun with exact scope labels.

- Complete nontrivial n=4 cap-2/5 scan: 358 labeled families, no negative
  `A_plus`; 22.67 s wall. Lower dimensions occur as zero-column embeddings.
- Complete nontrivial n=5 producer: 463,343 labeled families / 5,172 coordinate
  orbits. The new append-only evaluator completed 5,172 records in 58.57 s
  after a two-record resume pilot; minimum
  `A_plus = 0.006474313148642441...`, zero negative and zero near-zero values.
  A subsequent resume reported `new_evaluations: 0`.
- Complete n=6 sizes 3 through 7: 262,585 labeled families / 1,081 orbits,
  no negative values; 23.35 s wall.
- Complete declared n=6 two-block classes: 450 full-S6 orbits, two numerical
  negatives; 9.76 s wall. This is not a complete arbitrary-family n=6 census.

The strengthened standalone verifier wrote schema-2
`gate_b_unbounded_arb_v4.json`, certifying the full cell-union definition and
the ratio lower bound `>181/5000`. A fresh Python 3.9 environment with only
`python-flint==0.6.0` reproduced the pre-final v4 bytes in 63.50 s. A bounded
read-only `fable` review attempt emitted only an initialization event and
exited with status 1 before returning a verdict; no external-review claim is
made from that attempt.

Evidence status: **THEOREM HUMAN-PROVED; BASE AND RATIO MACHINE-CERTIFIED;
FINITE CENSUS OBJECTIVES NUMERICAL**. The finite controls are falsification
evidence only. The all-\(k\) conclusion uses the exact Bellman induction, not
finite-\(n\) optimality.

## E12 — post-fix adversarial replay

After every E11 correction, the complete Gate B regression suite passed 12/12
tests in 1.925 s test time / 2.81 s wall. This includes the 5,400 fixed-order
tensorization attack, ordered-pair-with-replacement check, full-cell-union
contract, v4 ratio-lower-certificate contract, direct square, and n=8 scope
checks.

Both append-only finite audits resumed without work:

```text
n=5 completed_evaluations=5172 new_evaluations=0
n=8 completed_evaluations=136  new_evaluations=0
```

The final v4 checker was then replayed in the fresh Python 3.9 environment. It
printed `PROVED_GATE_B_UNBOUNDED` in 57.19 s, and `cmp` proved its generated JSON
byte-identical to the checked-in
`certificates/gate_b_unbounded_arb_v4.json`. `latexmk` rebuilt the five-page
paper with no error. The SHA-256 manifest in `README.md` was recomputed only
after the final verifier, candidate, v4 certificate, tests, proof, definitions,
n=5 artifacts, TeX source, and PDF were fixed.

Status: **ADVERSARIAL REPLAY PASS**. No hard blocker remains. The theorem's
human-proof boundary remains explicit: the all-power Bellman tensorization is
ordinary mathematics backed by exact algebra and finite falsification tests,
not a formal-proof-assistant certificate.

## E13 — second adversarial round: dependency-independent certificate

The E11/E12 round left one substantive rigor gap: every rigorous number came
from a single arithmetic stack (`python-flint` Arb) and a single Bellman
implementation whose clamp classification decides which of three formulas is
used at each state. A sign error, an inverted clamp comparison, or a
`python-flint` bug would not have been caught by any existing check. This round
removed that single point of failure.

`verify_gate_b_dyadic.py` is a second standalone checker with no third-party
import at all. It uses fixed-point dyadic intervals with denominator
\(2^{160}\) and explicit outward rounding, gets \(\ln\) from a 64-term rational
atanh series with an exact tail bound, and gets \(\exp\) from range reduction
plus a 32-term Taylor tail bound and repeated squaring. It also avoids all
clamp decisions by *relaxing* every nondegenerate action interval to
\([0,1]\), which is a superset of the feasible one-sided interval and therefore
yields an upper bound on \(C_{+,\pi}\).

```sh
nice -n 10 ./.venv/bin/python -B uc/gate_b/verify_gate_b_dyadic.py \
  --write-certificate uc/gate_b/certificates/gate_b_unbounded_dyadic_v1.json
```

Wall time: 4.02 s. Result:

```text
Q               in [5.44735497738032198816069988441895782571019992263,
                    5.44735497738032198816069988441895782571019992583]
E_pi W_pi       <  5.93388505570408299883193502354321107163753306829
A_plus          <  -0.02717429110348639668622783049078199477185032339
                <  -1/40
```

Because the base bound is all the theorem needs, the powers give
`ratio > k/40 -> infinity` with no Arb input. Cross-checks: the dyadic \(Q\)
interval strictly contains the Arb \(Q\) ball, the dyadic \(\ln\)/\(\exp\)
enclosures strictly contain the Arb ones on sampled points, and the relaxed
per-order value strictly exceeds the Arb one-sided certificate for all 21
order orbits.

Dependency-free clean replay on the pristine system interpreter, no venv and no
packages:

```sh
nice -n 10 /Library/Developer/CommandLineTools/usr/bin/python3 -B \
  uc/gate_b/verify_gate_b_dyadic.py
```

Output was byte-identical to the checked-in certificate. A fresh
`python-flint==0.6.0` Python 3.9 environment also reproduced
`gate_b_unbounded_arb_v4.json` byte-identically in the same round (76.65 s for
both replays including environment creation).

Status: **INDEPENDENTLY CERTIFIED**, with the previous Arb certificate retained
as the sharper finite record.

## E14 — targeted falsification searches for the product lemmas

Two new resumable searches attack the load-bearing identities where the earlier
equal-dimension test could not.

```sh
nice -n 10 ./.venv/bin/python -B uc/gate_b/audit_tensorization.py
nice -n 10 ./.venv/bin/python -B uc/gate_b/audit_square_direct.py
```

`audit_tensorization.py` evaluated 292 ordered product cases in 66.41 s, 205 of
them with unequal block dimensions, including dead columns, duplicate columns,
singletons, full power sets, middle layers, and 120 seeded random pairs. For
every case it compared fixed-order \(Q_\pi\) and \(C_{+,\pi}\) on the product
against the induced block sums for *every* global order, and checked the exact
ordered-defect product. Worst gaps: \(1.78\times10^{-15}\) (Bellman),
\(8.88\times10^{-16}\) (fixed-order iid), \(1.15\times10^{-14}\) (order-averaged
iid); zero exact-defect failures; zero counterexamples. A resume reported
`new_evaluations: 0`.

`audit_square_direct.py` attacks the remaining representational objection: the
Arb square check uses a factored-fiber representation that already encodes the
product structure. This script materializes the 2,025-row square as plain
14-bit masks and runs the generic prefix-fiber evaluators on it. Five declared
global orders, 868.52 s total:

```text
consecutive        11.743725681706515 vs 11.743725681706515   3,047,424 states
blockswap          11.743725681706515 vs 11.743725681706515   3,047,424 states
reversed           11.826243246195125 vs 11.826243246195125   3,071,736 states
alternating        11.743725681706515 vs 11.743725681706515   4,946,031 states
seeded_interleaved 11.833153208527014 vs 11.833153208527015   3,325,789 states
```

The three distinct values show the test tracks genuine order dependence rather
than a single constant. Exact square combinatorics were recomputed here too:
2,025 rows, all counts 810, cap bound 810 met with equality, incidence 11,340
against Reimer threshold 11,122, 3,920,000 missing ordered joins, defect
`6272/6561`, active and separating.

A third audit closed the last "formula versus family" gap. Until now the
admissibility and defect statements for \(k\ge3\) were only integer arithmetic
derived from the proved formulas; no power above the square had ever been
instantiated.

```sh
nice -n 10 ./.venv/bin/python -B uc/gate_b/audit_power_admissibility.py
```

Wall time: 6.15 s, of which 6.04 s is the cube. It builds \(\mathcal F_k\) for
\(k\le3\) and recomputes everything from the constructed rows, counting missing
ordered joins by brute force through a dense membership table:

```text
k=1  45 rows      1600/2025 missing            64/81
k=2  2025 rows    3920000/4100625 missing      6272/6561
k=3  91125 rows   8227000000/8303765625        526528/531441 = 1-(17/81)^3
```

At \(k=3\) all 21 coordinate counts equal 36,450, which is exactly the cap
bound \(\lfloor2\cdot91125/5\rfloor\); incidence is 765,450 against Reimer
threshold 750,668; and the family is active and separating. Nothing here uses
the product lemmas, so the cube independently confirms Lemma 1 and the
admissibility claim one power beyond the square.

Status: **NO COUNTEREXAMPLE FOUND** under a deliberately hostile search
formulation; the identities remain proved in `PROOF.md`.

## E15 — symbolic verification of the load-bearing algebra

Both certificates and the product proof reuse three algebraic facts. They are
now machine-checked in sympy inside the regression suite rather than trusted:

1. the four transition weights sum identically to one;
2. the step objective equals `const + s*D + h(s)`, so a common additive shift of
   all four continuation values leaves \(D\) unchanged and adds exactly that
   shift to the value;
3. the unrestricted maximizer of `s*D + h(s)` is `1/(1+2^-D)`, the maximum is
   `log2(1+2^D)`, and the second derivative `1/(s(s-1)ln2)` is negative on
   \((0,1)\), so the objective is strictly concave.

Fact 3 is exactly the formula used by both the Arb "unconstrained" branch and
the dyadic relaxation, so a sign or normalization slip there would now fail a
test.

The suite also now machine-checks the finite-to-general step itself, in exact
integers: \(2\cdot45^k\equiv0\pmod5\) with
\(\lfloor2\cdot45^k/5\rfloor=18\cdot45^{k-1}\) for \(k\le64\), so the cap is met
with equality and never by rounding; the Reimer inequality
\(2^{2I}>m^m\) with \(m=45^k\) and \(I=126k45^{k-1}\) reduces exactly to the
single \(k\)-free witness \(45^{45}<2^{252}\), with the two exponent identities
\(2I=252\,k45^{k-1}\) and \(km=45\,k45^{k-1}\) checked for \(k\le64\) and the
direct \(m^m\) comparison checked for \(k\le2\); and both certified ratio
bounds \((k\delta)/(1-(17/81)^k)\) are strictly increasing and exceed
\(k\delta\) for \(\delta\in\{7/250,1/40\}\) over \(k\le200\), with
\(k=100{,}000\) exceeding 1,000.

A first version of that test asserted \(2^{2I}>m^m\) directly for \(k\le32\).
With \(m=45^{32}\approx10^{52}\) the right side has more than \(10^{53}\) bits,
so the process ran 1,800 s without finishing and was cancelled. The recorded
lesson is that the general Reimer claim must be reduced to its \(k\)-free
witness before being machine-checked; the direct comparison is only tractable
for the smallest powers.

The complete suite is 24 tests, 12.5 s. Independent re-derivation of the
claimed extremal ratios with 80-digit decimal arithmetic reproduced
`-A_+/eps` `= 0.0362591270367489647...` at the certified \(n=7\) base and
`0.03433549033144272...` at the \(n=8\) class maximum, matching the stored
float64 values to \(3\times10^{-16}\).

Smaller-\(n\) controls were rerun in this round: complete nontrivial \(n=4\)
cap-2/5 scan (358 families, zero negative, 24.30 s), complete nontrivial
\(n=5\) resume (`new_evaluations: 0`, minimum
\(A_+=0.006474313148642441\ldots\), 1.23 s), complete \(n=6\) sizes 3--7
(262,585 labeled families, 1,081 orbits, zero negative, 24.68 s), and the
declared \(n=6\) two-block classes (450 orbits, two negatives, 10.38 s).

Status: **ADVERSARIAL ROUND PASS.** The conclusion now survives removal of the
Arb dependency, removal of the clamp classification, removal of the factored
square representation, a hostile unequal-dimension product search, and direct
instantiation of the cube.

## E16 — claim-strength audit of every quoted decimal bound

Every decimal constant quoted in prose was re-checked against the certificate
endpoint it summarizes, in the correct direction. One overstatement was found
and fixed: the progress ledger asserted

```text
A_plus(B) < -0.0271742911034863966862278304908
```

while the certified dyadic upper endpoint is
`-0.0271742911034863966862278304907819947718503233903471735`. The quoted
truncation is *more* negative than the endpoint, so it claimed strictly more
than the certificate proves, by about \(1.8\times10^{-32}\). It now reads
`-0.0271742911034863966862278304907`, which the endpoint implies. The paper's
corresponding display previously ended in `\ldots`, which is ambiguous about
whether the certified digits continue; the ellipsis was removed so the
displayed rational is exactly what is implied.

The remaining quoted bounds were verified valid in the direction used:
`Q<5.4473549773803219881606998844190` and
`E_pi W_pi<5.9338850557040829988319350235433` in the ledger,
the same two truncated one digit shorter in `PROOF.md` and the paper, and the
ledger's `A_plus < -0.02717429110348639668622783049078199477185032339`, which
is less negative than the certified endpoint and therefore implied.

A second defect surfaced while re-verifying the SHA-256 manifest: the three new
audit summaries embedded run-dependent fields (`new_evaluations`, and for the
tensorization search `runtime_seconds`), so a resume rewrote a checked-in
artifact and broke its hash even with zero new work. The persisted summaries now
contain only run-independent content; the counters are printed to stdout. Two
consecutive resumes of all three audits now leave all three artifacts
byte-identical, matching the convention already used by the \(n=8\) census.

Status: **TWO DEFECTS FOUND AND FIXED.** No claim now exceeds its certificate,
and every checked-in audit artifact is reproducible byte-for-byte on resume.

## E17 — second infinite construction from the six-coordinate falsifier

The highest-value remaining robustness gap was candidate dependence.  The two
existing arithmetic certificates both reconstructed the same 45-row
seven-coordinate family and used its 21-order symmetry quotient.  A
candidate-specific reconstruction or orbit error could therefore invalidate
both finite inputs at once even though their interval arithmetic was disjoint.

The older \(n=6\) two-block census contained a separate 25-row negative family,
the union of the \(3+3\) cells
\[
(0,0),(0,1),(1,0),(1,2),(2,1).
\]
An exploratory run of the clamp-free dyadic recurrence evaluated all 720
coordinate orders and found enough margin to target the exact rational bound
\(A_+<-1/80\).  This was promoted to the committed checker
`verify_gate_b_n6_dyadic.py`, which reconstructs the family from the cells,
uses no order quotient, and verifies:

```text
size=25; counts=(10,10,10,10,10,10); incidence=60; Reimer threshold=59
missing ordered joins=444/625; closure success=181/625
720 orders; 584,784 total Bellman states; 10 root intervals each repeated 72 times
Q < 4.6163838598571102600389234608954804813132320743559048317
W < 5.0343539516297632216898156407608692922319055595635882173
A_plus < -0.0125897106568747589615082453240509575649305961627765509
       < -1/80
```

The final system-Python run was one low-priority process and took 96.20 s.
The corresponding two-test contract replay under the project environment took
42.70 s.  The existing, separately implemented 256-bit Arb falsifier was then
replayed in 13.04 s; it reconstructed the same rows and exact defect and gave
\(A_+<-0.0136721077321777735\ldots\).  That route was promoted to the compact
committed checker `verify_gate_b_n6_arb.py`, which reuses only the generic
standalone Arb evaluator, writes a deterministic certificate in 0.21 s, and
proves the rational bound \(A_+<-17/1250\).
The isolated Python 3.9 environment reproduced the new Arb JSON
byte-identically in 0.41 s; the same environment reproduced the original
\(n=7\) Arb v4 certificate byte-identically in 59.77 s.

For \(\mathcal G_k=\mathcal D^{\boxtimes k}\), exact product arithmetic gives
\[
\varepsilon_\vee(\mathcal G_k)=1-(181/625)^k,\qquad
A_+(\mathcal G_k)<-k/80.
\]
Every coordinate count is \(10\cdot25^{k-1}=(2/5)25^k\).  The average row size
is \(12k/5\), and the one \(k\)-free witness \(25^5<2^{24}\) proves Reimer
strictly for every power.  Thus this second family alone yields
\[
\frac{-A_+(\mathcal G_k)}{\varepsilon_\vee(\mathcal G_k)}
>\frac{k/80}{1-(181/625)^k}>\frac{k}{80}\to\infty.
\]

A direct-square counterexample search then materialized
\(\mathcal D^{\boxtimes2}\) as 625 plain 12-bit row masks and passed those masks
to the generic float64 prefix-fiber evaluators.  It directly counted 357,864
missing ordered joins out of 390,625, exactly
\(1-(181/625)^2\), and verified cap equality, Reimer, activity, and separation.
Five hostile global orders (consecutive, block-swapped, reversed, alternating,
and seeded interleaved) had zero Bellman product gap and worst iid gap
\(1.78\times10^{-15}\).  The 13.20 s search is append-only and resumable; its
immediate replay reported `new_evaluations: 0` and left both artifacts
byte-identical.

The first isolated-system replay of the new dyadic checker exposed a packaging
defect: Python's `-I` mode did not place the script directory on `sys.path`, so
the sibling arithmetic module could not be imported.  The checker now inserts
its resolved sibling directory explicitly.  With an otherwise empty
environment, system Python 3.9 then reproduced both dyadic JSON certificates
byte-identically (8.9 s for \(n=7\), 96.68 s for all 720 \(n=6\) orders).
The expanded suite passes 29/29 tests in 54.930 s test time / 55.84 s wall.

Status: **SECOND INFINITE CONSTRUCTION CERTIFIED.**  The Gate B conclusion no
longer depends on the \(n=7\) candidate or any order-orbit quotient.  Both
finite bases now survive the same two arithmetically disjoint stacks: an exact
clamp-classified 256-bit Arb recurrence and a clamp-free standard-library
dyadic relaxation.

## E18 — exact rational two-sided certificate for both bases (2026-08-27)

Motivation: after E17 the load-bearing sign of \(A_+\) still rested on two
*numerical* stacks. The Arb route needs `python-flint` and a three-way clamp
classification; the dyadic route needs a 160-bit interval class, hand-derived
atanh/exp tail bounds, and — more seriously — it relaxes every feasible action
interval to \([0,1]\), so it certifies a relaxation of \(C_+\) rather than
\(C_+\). Neither route can distinguish "the true optimum is negative" from
"some upper bound of it is negative".

Command:

```sh
nice -n 10 /Library/Developer/CommandLineTools/usr/bin/python3 -B \
  uc/gate_b/verify_gate_b_rational.py --bases n6,n7 --orders all \
  --progress 500 \
  --write-certificate uc/gate_b/certificates/gate_b_unbounded_rational_v1.json
```

Wall time: 250.35 s in one low-priority process on the pristine system
interpreter with no third-party package. Peak observed share of total machine
CPU: 2.4 %.

The evaluator keeps the exact feasible interval \([s^*(p,r),U(p,r)]\) and needs
only three certified primitives, all integer arithmetic:

1. bit-by-bit binary logarithm from \(\log_2y=(b+\log_2(y^2/2^b))/2\), with
   monotone rounding of the state and the tail in \([0,2^{-N}]\);
2. `math.isqrt` towers \(r_{j+1}=\lceil\sqrt{r_j}\rceil\ge2^{2^{-(j+1)}}\) for
   \(2^x\) from above;
3. a concavity case split for \(\max_{s\in[\ell,u]}[h(s)+cs]\): exact at an
   endpoint when the derivative has a sign there, otherwise the Fenchel value
   \(\log_2(1+2^c)\).

Results, over **all** coordinate orders with no automorphism quotient:

| base | orders | Bellman states | distinct enclosures | multiplicity |
|---|---:|---:|---:|---:|
| \(\mathcal D\), \(n=6\) | 720 | 584,784 | 10 | 72 each |
| \(\mathcal B\), \(n=7\) | 5,040 | 13,381,680 | 21 | 240 each |

The multiplicities were *not* assumed: they are the observed partition of the
all-order value list, and they independently reproduce the 72- and 240-element
automorphism groups and the 10 and 21 orbit counts.

```text
A_+(D) in [-0.0136721077321777735618599156610223516934,
           -0.0136721077321777735618599129978279356672]
A_+(B) in [-0.0286491867944683178160269819326732064410,
           -0.0286491867944683178160269764548230805862]
```

Both enclosures have width below \(6\times10^{-27}\); both upper endpoints beat
the sharp frozen targets \(-17/1250\) and \(-7/250\); and both Arb values from
E2 and E17 lie strictly inside, which retroactively shows the clamp-ambiguous
relaxation in those runs cost less than \(10^{-26}\).

Resumability: 5,760 append-only records with an explicit `schema` field. A
complete rerun reports `new_evaluations: 0`, recomputes a deterministic sample
and aborts on mismatch, and counts foreign-schema or truncated lines instead of
trusting them. An intermediate run of an earlier record schema was discarded
for exactly that reason rather than reused.

Growth rate: the same enclosures give the exact linear statement
\(c_{\rm cl}^\star(n)>\frac1{250}n-\frac7{250}\) from \(\mathcal B\) and
\(\frac{17}{7500}n-\frac{17}{1250}\) from \(\mathcal D\), while \(Q\ge0\),
\(C_+\ge0\) and \(m\le2^n\) give \(-A_+\le\log_2m\le n\). Hence on families
with \(\varepsilon_\vee\ge\varepsilon_0\) the ratio is \(\Theta(n)\).

Exact product-lemma repeat: `audit_tensorization_exact.py` ran 35 ordered
products over all 15,010 global orders in 62.18 s (20 new cases resumed onto 15
earlier ones, demonstrating case-level resume). Zero intersection failures,
zero exact defect-identity failures, worst \(Q\) discrepancy
\(1.4\times10^{-27}\), worst fixed-order Bellman discrepancy exactly zero. An
immediate rerun reported `new_evaluations: 0` in 0.09 s.

Clean-room replay: in a stripped environment
(`env -i OMP_NUM_THREADS=1 HOME=... PATH=/usr/bin:/bin`) the pristine system
interpreter reproduced `certificates/gate_b_unbounded_rational_v1.json`
byte-identically with `--fresh` in 251.47 s. Without `--fresh` the same
environment replayed the committed checkpoint in 0.92 s, reporting
`new_evaluations: 0`, `resumed_records` 720 and 5040, and two recomputed
records per base with no mismatch. Only the bookkeeping fields differ between
those two modes.

Suite: 42 tests pass in 60.2 s under the project environment, including
Arb-versus-rational bracketing of every primitive, 65-point sampled domination
of the constrained maximum, per-orbit enclosure agreement on both bases,
sharpness against the dyadic relaxation, state-canonicalization invariance,
orbit-versus-all-order certificate equality, foreign-schema and
truncated-checkpoint rejection, an exact one-case audit replay, and the
exact-integer growth corollary. The paper draft rebuilds to eight pages with no
undefined reference.

Status: **EXACTLY CERTIFIED, BOTH BASES.** The finite input to the Gate B
theorem is now an exact rational fact about the true objective, established
without any interval library, any clamp classification, any action-set
relaxation, and any symmetry quotient.

## E19 — the local regime: is the ratio bounded when the defect is small? (2026-08-27)

Motivation: E18 closed the arithmetic of the finite input, so the highest-value
remaining gap is no longer a certificate but the one escape the theorem leaves
open. Every family in the Gate B construction has defect at least `444/625`
(Lemma 5, proved this round), so the published proof says nothing about

    c_loc = lim_{e -> 0+} sup { -A_+/eps_vee : admissible, A_+ < 0, 0 < eps_vee <= e }.

Four elementary facts were proved first and are now in
[PROOF.md](PROOF.md): failures come in pairs so `eps_vee >= 2/m^2` (this also
sharpens a `1/m^2` remark that was already in the document); the defect
multiplies under products so the power mechanism cannot enter the local regime
at all; hence `-A_+/eps_vee <= m^2 log2(m)/2` and reaching defect `e` forces
`m >= sqrt(2/e)`; and a family with zero defect would be a Frankl
counterexample, so proving a positive universal defect floor at cap `2/5` is
strictly harder than beating the published union-closed frontier by `~0.018`.
That last point is the blocker, and it is why this round measures the regime
instead of trying to prove it empty.

### Instruments

`search_local_defect.py` minimizes the closure defect alone, which costs
`O(m^2)` set lookups and no entropy evaluation. `--mode exact` enumerates every
admissible family for `n <= 4`. `hunt_local_ratio.py` has three modes:
`hunt` minimizes `A_+ + lam*eps_vee` from random starts, `descend` minimizes the
defect from a certified negative family under a hard negativity constraint, and
`sweep` minimizes `A_+` under a hard defect cap.

Commands (one worker, `nice -n 19`, total machine CPU below 50 % throughout):

```sh
OMP_NUM_THREADS=1 nice -n 10 ./.venv/bin/python -B \
  uc/gate_b/search_local_defect.py --mode exact --dimension 4
OMP_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -B \
  uc/gate_b/search_local_defect.py --mode search --dimension 7 \
  --sizes 25,33,40,45 --restarts 6 --steps 120000
OMP_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -B \
  uc/gate_b/search_local_defect.py --mode frontier --dimension 7
OMP_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -B \
  uc/gate_b/hunt_local_ratio.py --mode descend --dimension 7 \
  --restarts 3 --steps 1600 --screen 48
OMP_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -B \
  uc/gate_b/hunt_local_ratio.py --mode sweep --dimension 7 --sizes 45 \
  --restarts 2 --steps 1000 --screen 24 \
  --defect-caps 0.79,0.74,0.70,0.66,0.62,0.58,0.52,0.46,0.40,0.35
```

### Result 1 — the defect does not want to be small

Exact enumeration of *every* admissible family on `n <= 4` gives a minimum
defect of `2/9` at `m = 3` and `6/25` at `m = 5`; no admissible family on four
coordinates has defect below `6/25 = 0.24`. Searched upper bounds on
`eps_*(n, m)` for `n = 5..8` and every feasible size never went below `2/9`.
In the size range that Lemma 6 actually forces, the minima are

| n | m | best defect found | pair floor `2/m^2` | ratio to floor |
|---|---|---|---|---|
| 6 | 25 | `184/625 = 0.2944` | `0.0032` | 92x |
| 7 | 45 | `692/2025 = 0.3417` | `0.00099` | 346x |
| 8 | 63 | `1468/3969 = 0.3699` | `0.00050` | 734x |
| 8 | 80 | `163/400 = 0.4075` | `0.00031` | 1304x |

The combinatorial floor decays like `1/m^2` while the achieved minimum does not
decay at all, and the gap widens monotonically. Note the `n = 7`, `m = 45` row:
at the published base's own parameters there are admissible families with defect
`0.3417`, less than half the base's `64/81 = 0.7901`. Low-defect admissible
families are plentiful; the question is their sign.

### Result 2 — low defect forces a large *positive* objective

Evaluating `A_+` over all `n!` orders on the 38 minimum-defect witnesses for
`n = 6, 7` and `m >= 15` gave `A_+ > 0` in every case, from `+0.0588` to
`+0.5607`. The contrast at the certified bases' own sizes is the finding:

| family | defect | `A_+` |
|---|---|---|
| `n=6, m=25` minimum-defect witness | `0.2944` | `+0.05877257` |
| `n=6, m=25` certified base `D` | `0.7104` | `-0.01367211` |
| `n=7, m=45` minimum-defect witness | `0.3417` | `+0.08470357` |
| `n=7, m=45` certified base `B` | `0.7901` | `-0.02864919` |

Driving the defect down drives the objective up, at both dimensions.

### Result 3 — the negativity region terminates well away from zero

A first `hunt` sweep over `lam` in `{0.04, 0.06, 0.10}`, 48 restarts at
`n = 6, 7`, found **zero** negative families. That was a search failure, not a
landscape fact: it never rediscovered either certified base, and its best score
at `n = 6` (`+0.0311` at `lam = 0.04`) was much worse than `D`'s own `+0.0147`.
Diagnosis: a 12-order screen is too noisy relative to a signal of size `0.03`,
and random starts are infeasible. Recorded as a failed route rather than as
evidence.

The `descend` mode fixes both by starting *at* a certified negative family and
using the exact all-order objective at `n = 6`. Results:

| n | restarts | negatives visited | lowest defect with `A_+ < 0` | `A_+` there | ratio |
|---|---|---|---|---|---|
| 6 | 5 | 1 | `444/625 = 0.7104` (= base `D`, no improvement) | `-0.01367211` | `0.0192456` |
| 7 | 3 | 64 | `1336/2025 = 0.6598` | `-0.00085519` | `0.0012962` |

So the lowest defect at which negativity survives moved down from `0.7901` to
`0.6598` at `n = 7`, and at that point `A_+` is only `-0.00086`: the ratio
collapsed by a factor of 28, from `0.0363` to `0.0013`. This is the opposite of
the behaviour a local counterexample needs. Combined with Result 2, the
negativity region at reachable dimensions is an interval of large defect, and
`-A_+` vanishes as its lower endpoint is approached.

### Result 4 — the frontier curve, after fixing the seeding

The first `sweep` run seeded every cap from the minimum-defect witness and could
not climb back into the negativity region: at cap `0.79` it reported `+0.037`
although a negative family was feasible there. Recorded as a defect of the
instrument. The fix was to seed each cap from the *best family already known*
under that cap, pooled from both checkpoints and re-verified, which also makes
the curve monotone by construction. Reseeded records use the `sweep2:` key
namespace so the loose points are never silently reused. The curve at `n = 7`,
`m = 45`, all points recomputed over all 5040 orders:

| defect cap | least `A_+` found | defect attained | negative |
|---|---|---|---|
| 0.79 | `-0.01628321` | `314/405` | yes |
| 0.74 | `-0.01489899` | `56/81` | yes |
| 0.70 | `-0.01489899` | `56/81` | yes |
| 0.66 | `-0.00085519` | `1336/2025` | yes |
| 0.64 | `+0.01171204` | `1222/2025` | no |
| 0.58 | `+0.02129440` | `1154/2025` | no |
| 0.52 | `+0.04399390` | `206/405` | no |
| 0.46 | `+0.05599499` | `184/405` | no |
| 0.40 | `+0.06939557` | `796/2025` | no |
| 0.35 | `+0.07951401` | `236/675` | no |

### Result 5 — the low-defect witness is now certified exactly

The witness at defect `1336/2025` was promoted from a float64 find to an exact
fact. `verify_gate_b_rational.py` gained an `ExplicitBase` type, because the two
published bases are reconstructed from weight cells and a search witness has no
such description, while every exact check that matters is independent of how the
rows arose. Command and result:

```sh
OMP_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -B \
  uc/gate_b/verify_gate_b_rational.py --bases n7lo --orders all \
  --write-certificate uc/gate_b/certificates/gate_b_lowdefect_rational_v1.json
```

```text
verdict PROVED_GATE_B_UNBOUNDED_EXACT_RATIONAL
defect      1336/2025          counts [18]*7   incidence 126 >= 124
A_+ lower  -0.0008551853337231417110078957932
A_+ upper  -0.0008551853337231417110078902901
width       0.0000000000000000000000000055031
targets     a_plus_upper_lt = -1/1250, -1/1200
orders      5040 of 5040, automorphism_count 1, 12,716,352 Bellman states
```

Run time 159 s; a rerun resumes with `new_evaluations: 0`, 5040 resumed records
and 2 sampled rechecks. The trivial automorphism group is worth noting: unlike
both published bases this family admits no symmetry quotient even in principle,
so all 5040 enclosures are distinct and were computed independently.

This makes `e* <= 1336/2025 = 0.6597...` a certified statement, improving the
previous `64/81 = 0.7901...`, and supplies a third infinite construction with
`eps_vee = 1-(689/2025)^k` and ratio above `k/1200`. Its certified ratio lower
bound is `27/21376 = 0.001263...` and the enclosure pins the value to
`0.0012962...`, twenty-eight times below the published base's `0.0362591`, which
is the point: the family is far *less* negative and far *closer* to union-closed.

Suite: 55 tests pass in 61.9 s, including thirteen new ones for the local regime
(Lemma 4 parity on 400 random families, exact product multiplicativity of the
defect, cross-implementation agreement of the integer Reimer threshold, the
published defects and calibration ratios of both bases, the exact `n <= 4`
minima, the witness-versus-schema contract, seeder coverage of every feasible
size, the hard defect cap of the sweep, and two pinning the new certificate's
defect, trivial automorphism group and negative enclosure).

Status: **the local question stays OPEN, and it is now bracketed.** Certified:
negativity survives to defect `1336/2025`. Searched and not found: any negative
family below defect `0.64` at `m = 45`, or any negative family at all among the
minimum-defect witnesses. Blocked: by Proposition 7, proving the regime empty
implies a `2/5` Frankl bound, about `0.018` beyond the published frontier. The
measurements point away from a local counterexample without excluding one.

## E20 — the size ceiling, and the eighth coordinate (2026-08-27)

Motivation: E19 bracketed the local regime but its instruments all worked at
`n = 7`, `m = 45`. Two questions were left. Is there anything to prove about the
repair ratio that does not go through Frankl? And is `n = 7` actually where the
frontier lives?

### Result 1 — the numerator's growth is settled, both sides

Command:

```sh
OMP_NUM_THREADS=1 nice -n 19 python3 -B uc/gate_b/bound_local_regime.py \
    --bases n6,n7,n7lo --max-dimension 4
```

Incidence is at most `n*floor(2m/5)` and Reimer demands at least
`m log2 m / 2`, so every admissible family satisfies `log2 m <= 4n/5`, tested
exactly as the integer inequality `m^5 <= 2^(4n)`. With `Q, C_+ >= 0` this gives
the first unconditional `n`-explicit ceiling `-A_+ <= 4n/5`, and the certified
powers supply the matching linear lower bound: the numerator is `Theta(n)`.
Verified for every admissible size at every `3 <= n <= 16`. The ceiling is close
to sharp -- the largest admissible size at `n = 11` is 445, i.e. `log2 m =
8.7977` against `8.8`. The repository already carried this inequality per base
as the `reimer_witness` string without recognising it as a general ceiling.

### Result 2 — two defect floors, one of them attained

`eps_vee >= ((8/5) sbar - M) / (n - M)`, because the cap forces
`E|X or Y| >= (8/5) sbar` while a successful union has size at most `M`; and
`eps_vee >= 1 - m^-2 sum_A N(A)^2` from downset counting. Exhaustive audit over
all 366 admissible families with `n <= 4`: zero violations, and the capacity
floor is attained with slack exactly `0`. Consequence (Corollary 12): the ratio
is `O(n)` on every family *without* a dominant set, so approaching the local
regime **requires** a set of size at least `(8/5) sbar >= (4/5) log2 m`. The
`n = 6` base is below that threshold and its floor `7/25` is active; both `n = 7`
bases are above it, and every `n = 8` witness below contains `[8]` itself.

### Result 3 — the certified bases admit no addition whatsoever

25 and 45 are the *largest* admissible sizes at `n = 6` and `n = 7`, all three
certified bases attain them with every degree exactly at the cap, and the
verifier finds **zero** admissible additions over all `2^n` candidate sets. So
the one amplification route that lowers the defect and raises `log2 m` at the
same time -- adding the missing unions -- is empty here, not merely unpromising.
A first greedy attempt confirmed this before the proof: every base reported
"no admissible addition (cap/Reimer block)" at step 0.

### Result 4 — the combinatorial floor is not the obstruction

At `n = 7`, `m = 45` the recorded defect frontier over *all* admissible families
is `692/2025 = 0.3417`, while the best certified negative sits at
`1336/2025 = 0.6598`: a gap of `0.318`. Measuring `Q` along the way explains it.
At fixed `m = 45`, `log2 45 = 5.4919` is fixed and `Q` rises as the defect falls
-- `5.4474` at `0.79`, `5.4758` at `0.66`, `5.5608` at `0.34` -- crossing the
entropy and turning `A_+` positive. Negativity is a race against `log2 m` and is
won only at the maximum admissible size, which is why every base sits there.

### Result 5 — the eighth coordinate, and a stronger base

That diagnosis says the frontier should move when the size ceiling moves, so the
untried `S_2 x S_6` class at `n = 8` was enumerated:

```sh
OMP_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -B \
    uc/gate_b/search_n8_block_symmetric.py --k 2 \
    --checkpoint uc/gate_b/experiments/n8_k2_checkpoint.jsonl \
    --result uc/gate_b/candidates/n8_k2_census.json
```

21 cells, 2,097,151 raw masks, 2,272 canonical admissible families, 419 s. Then
every candidate with float `A_+ < 0.01` was re-evaluated in exact rational
arithmetic:

```sh
OMP_NUM_THREADS=1 nice -n 19 python3 -B uc/gate_b/certify_class_negatives.py \
    --census uc/gate_b/experiments/n8_k2_checkpoint.jsonl --threshold 0.01
```

36 families evaluated exactly, 20 certified negative. Certified improvements:

| quantity | old | new | family |
|---|---|---|---|
| lowest defect, separating | `1336/2025 = 0.659753` | `1144/1875 = 0.610133` | `n8lo`, `m = 75` |
| lowest defect, any admissible | `1336/2025` | `139/245 = 0.567347` | `m = 70` |
| highest repair ratio, separating | `0.0362591` | `0.0370512` | `n8hi`, `m = 75` |
| highest repair ratio, any admissible | `0.0362591` | `0.0580391` | `m = 70` |

Sizes 70 and 75 are **infeasible at `n = 7`** (`7*28 = 196 < R_70 = 215`); the
eighth coordinate supplies the missing incidence, and in the two 70-row families
it is an exact duplicate of the first. Coordinate cloning changes neither defect
nor size, so it is exactly the mechanism Lemma 9 predicts: raising `n` is the
only way to raise the size ceiling. Normalization is a reporting condition in
`DEFINITIONS.md`, not a condition in the displayed supremum, so both conventions
are reported.

### Result 6 — the float census produced a false negative

The 80-row family at defect `2553/3200` reads `-0.00026687` in float64 and
`+0.000454098` exactly. Float/exact gaps reach `2e-3` elsewhere in the class,
far beyond float64 noise, so the census ordering is not an ordering of the exact
values. Both frozen rational targets for the new bases were first written from
float values and both were rejected by the verifier before being corrected
against the exact enclosures -- the intended behaviour of the frozen-target
gate, observed working. The mirror check was then run: all 15 families with
float `A_+` in `[0, 0.01)` are exactly non-negative, so no negative is hiding in
the positive tail.

Suite: 66 tests pass, including eleven new ones for the size ceiling, both
defect floors, the maximality of the certified bases, the dominant-set
equivalence, cross-implementation agreement of the standalone checker, the two
`n = 8` bases, the infeasibility of size 75 at seven coordinates, and the false
negative.

Status: **the local frontier moved from `0.6598` to `0.6101` (separating) and
`0.5673` (any admissible), and the finite ratio witness beat the published base
for the first time.** The numerator's growth law is closed at `Theta(n)`. The
local question itself stays OPEN and Frankl-equivalent by Proposition 7.

## E21 — the growth constant improves, and separation turns out to be free (2026-08-27)

Motivation: E20 left a specific question. The two non-separating 70-row
witnesses had the largest `-A_+` per coordinate seen anywhere,
`0.037548/8 = 0.004694` against the `n=7` base's `0.028649/7 = 0.004093`, so
they *would* improve the certified asymptotic growth constant — if their
Cartesian powers are admissible and the product lemmas apply to a base that is
not separating. That had to be settled by reading the proof, not by assuming it.

### Result 1 — no lemma uses separation

Read against `PROOF.md`: Lemma 1 uses only that the two block join events are
independent under two independent uniform product rows. Lemmas 2 and 3 use only
that a uniform product row has independent blocks, that the conditional
marginals, feasible interval and four transition probabilities at a coordinate
depend on that block's prefixes alone, and that the four transition
probabilities sum to one. The additivity corollary uses only
`|H| = |F| |G|`. Separation appears exactly once, in the last bullet of the
admissibility section, where it is *concluded* for the powers of a normalized
base — never assumed. Recorded as Proposition 16.

### Result 2 — one base-level Reimer check certifies every power

For a base with size `m`, incidence `I` and every degree equal to `c`, the power
`F^k` has size `m^k`, degree `c m^(k-1)` and incidence `k I m^(k-1)`.

* Cap: `c <= floor(2m/5) <= 2m/5`, so the degree is an integer at most `2m^k/5`,
  hence at most `floor(2m^k/5)`. For `m = 70` it is *exactly* the cap, because
  `70/5 = 14` is an integer.
* Reimer: the requirement is `k I m^(k-1) >= ceil(k m^k log2(m)/2)`, and an
  integer dominating a real number dominates its ceiling, so it suffices that
  `2I >= m log2 m`, i.e. `m^m <= 2^(2I)` — which is exactly the base condition
  `I >= R_m`, since `R_m` is the least `r` with `2^(2r) >= m^m`.

A first attempt at this got the algebra wrong (`m^(2I) < 2^(4I)`, which is false
for `m > 4`) and the audit correctly returned `FAILED`; the corrected witness
`70^70 <= 2^448` holds with `log2(70^70) = 429.05` against `448`.

```sh
OMP_NUM_THREADS=1 nice -n 19 python3 -B uc/gate_b/audit_clone_power.py \
    --base n8clone_hi
```

Verdict `EVERY_POWER_ADMISSIBLE`: degrees equal the cap at every power up to 40,
Reimer strict throughout, and the instantiated square — 4,900 rows rebuilt from
the product and recounted with no product-aware shortcut — has defect
`210171/240100 = 1-(173/490)^2` exactly.

### Result 3 — the certified growth constant improves by 75/64

Registering the two cloned witnesses as `ClonedCoordinateBase` (which reports
separation instead of asserting it, and carries the duplicate column pair in its
facts) and certifying:

```sh
OMP_NUM_THREADS=1 nice -n 19 python3 -B uc/gate_b/verify_gate_b_rational.py \
    --bases n8clone_lo,n8clone_hi --orders orbits \
    --write-certificate uc/gate_b/certificates/gate_b_n8_clone_rational_v1.json
```

`n8clone_hi`: defect `317/490`, `A_+ <= -3/80`, certified ratio `>= 147/2536 =
0.0579653`, asymptotic slope `3/640 = 0.0046875` against the published
`1/250 = 0.004` — an improvement by exactly `75/64 = 1.1719`. With Corollary 10
the numerator is now pinned between `(3/640)n - 3/80` and `4n/5`.

The mechanism is not the one to guess. Cloning a coordinate changes neither the
size, the defect, nor the join structure; it only adds incidence. What it buys is
*feasibility*: size 70 is impossible at seven coordinates, where `7*28 = 196`
falls short of `R_70 = 215`. The eighth coordinate is spent entirely on clearing
Reimer, and the reward is a base whose `-A_+` per coordinate beats anything
available at seven. The cost is separation, which the definitions do not charge
for.

### Result 4 — the S_3 x S_5 class, and the best separating ratio

```sh
OMP_NUM_THREADS=1 nice -n 19 ./.venv/bin/python -B \
    uc/gate_b/search_n8_block_symmetric.py --k 3 \
    --checkpoint uc/gate_b/experiments/n8_k3_checkpoint.jsonl \
    --result uc/gate_b/candidates/n8_k3_census.json
```

24 cells, 16,777,215 raw masks, 13,470 canonical admissible families, 54
float-negative, 2,743 s. Exact re-ranking of all 110 candidates with float
`A_+ < 0.01` took 1,194 s and certified 52. Its best, `n8max`, sits at the
**maximum** admissible size `m = 80` with every degree at the cap 32 and
incidence `256 = 8*32`, the largest any `n=8` family can carry: defect
`2343/3200`, `A_+ <= -19/625`, certified ratio `>= 2432/58575 = 0.0415194`, the
best among separating families, improving `n8hi`'s `0.0370192`. Its slope
`19/5000 = 0.0038` does *not* improve the growth constant. `k=3` moved neither
the defect frontier nor the slope.

### Result 5 — the census screen fails in both directions

`k=3` settles a question `k=2` left open. Three of its 54 float negatives are
exactly non-negative, and — new — **one family the census called non-negative is
exactly negative**: defect `94/125`, float `+0.00003418`, exact
`-0.00005157`. So the census both raises false alarms and *misses* witnesses.
Largest census/exact gap in the class is `3.24e-3`, above `k=2`'s `2.03e-3`.
E22 below identifies the cause, which is **not** float64 error.

The practical consequence is already in force: the re-ranking threshold is
`+0.01`, not `0`, precisely so that missed witnesses are recovered. Across both
re-ranked classes the float census reported 75 negatives; the exact route
certifies 72, of which one was never reported as negative at all.

Suite: 77 tests pass, including nine new ones for the cloned bases, the
separation flag (an impostor declaring separation must fail loudly), the
equivalence of `I >= R_m` with `m^m <= 2^(2I)` on all eight registered bases, the
power audit, the slope improvement, `n8max` at maximum incidence, and the
two-directional float failure.

Status: **the growth constant is improved and the numerator is bracketed
`(3/640)n - 3/80 <= Lambda(n) <= 4n/5`.** The local question stays OPEN and
Frankl-equivalent.


## E22 — the census discrepancy is an order-set bug, not float error (2026-08-27)

Motivation: E20 and E21 both recorded "float/exact gaps up to `3.24e-3`" and
attributed them to float64. That attribution was never tested, and a `3e-3`
error in a quantity assembled from values of size `6` is about eleven orders of
magnitude worse than float64 noise. So it was tested.

### Result 1 — the float evaluator is not at fault

Evaluating `Q` and `C_+` for `n8lo` and `n8hi` over the *same* 28 automorphism
orbit representatives:

| base | Q gap (float vs exact) | C_+ gap | orders disagreeing beyond 1e-9 |
|---|---|---|---|
| `n8lo` | `-1.24e-14` | `0.0` | 0 of 28 |
| `n8hi` | `-2.93e-14` | `0.0` | 0 of 28 |

So `shapley_n6_shared_bellman.one_sided_costs` and `shapley_join_loss.shapley_iid`
agree with the exact rational route to `1e-14`. Nothing needs fixing there.

### Result 2 — the census averages over the wrong orders

For `n8lo`, comparing the three order sets with the *same* float evaluator:

```text
C_+ over census block reps : 6.683920061576
C_+ over Aut orbit reps    : 6.707147997083
C_+ over ALL 40320 orders  : 6.707147997083
```

The orbit route is exact to the last digit. The census is off by `-2.32e-2`,
which times `alpha = 0.0356069` is `-8.27e-4` in `A_+` — exactly the observed
`n8lo` gap (`-0.002859451` census against `-0.002032376` true).

### Result 3 — why, precisely

`block_order_representatives(k)` returns one order per `S_k x S_(8-k)` pattern.
That is a valid representative system only when the family's automorphism group
**is** that block group. Measured for `n8lo`:

* `|Aut| = 1440`, which is exactly `|S_2 x S_6|` — so a size check would pass;
* but only **120** of those 1440 elements preserve the block partition
  `{0,1} | {2..7}`, so `Aut` is a *different* group of the same order;
* consequently the 28 block-pattern orders hit only **12 of the 28** true
  automorphism orbits, double-counting six of them with multiplicities
  `[(1,6),(2,2),(3,1),(4,1),(5,1),(6,1)]` and missing sixteen entirely.

At `k=1` the block group *is* the automorphism group, and the `k=1` census value
is exactly right: re-ranked exactly, its negative reads `-0.027762696467` and
ratio `0.034335490331`, matching the published figures to every digit. The bug
is class-dependent, which is why it went unnoticed.

### Result 4 — a rigorous screening bound

Any average over a subset of orders lies within
`alpha * (max_pi C_+ - min_pi C_+)` of the true average. Measured over all
40320 orders per base:

| base | min C_+ | max C_+ | alpha * spread |
|---|---|---|---|
| `n8lo` | 6.674547 | 6.751629 | `0.002745` |
| `n8hi` | 6.677357 | 6.702007 | `0.000878` |
| `n8clone_hi` | 6.504946 | 6.648994 | `0.005129` |
| `n8max` | 6.684842 | 6.822257 | `0.004893` |

Worst bound `5.13e-3` against the re-ranking threshold `+0.01`: a `1.95x`
margin. That is why the threshold is `+0.01` and not `0`, and it is now the
documented reason rather than a guess.

Decision: **the census stays as a screen and its docstring now says so.** Its
evaluation logic is not changed, because fixing it properly needs the true
automorphism group per family — 40320 permutation tests each, times 13,470
families — while the two-stage architecture (cheap screen, exact re-rank) already
produces the correct answer and is what every label already comes from. What was
wrong was the *description*, in three documents, and that is corrected.

Suite: 79 tests pass, including two new ones that pin the diagnosis: the census
order set covers 12 of 28 orbits while the honest system covers 28 of 28, and
the screening bound is below `1/100` for every registered base.

## E23 — cloning is a defect-free amplifier, and it has a ceiling (2026-08-27)

Motivation: a review nit pointed out that a registry comment cited `7*28=196 <
R_70` while talking about size 75 (the right numbers there are `7*30=210 <
R_75=234`; the conclusion was unaffected). Checking it exposed something the
narrative had wrong. The comment's *story* — "the extra coordinate is what
admits the size" — cannot be what `n8tiny` is doing, because

```text
m=15 n=7: cap=6  n*cap=42  R_15=30  feasible=True
m=15 n=5: cap=6  n*cap=30  R_15=30  feasible=True
```

Size 15 is feasible down to five coordinates. So the clones in `n8tiny` are not
buying feasibility. They are buying something else.

### Result 1 — cloning is exactly free in size and defect

Duplicating a coordinate maps each row to itself with one bit repeated, a
bijection commuting with union. So `m` is unchanged, `eps_vee` is unchanged,
every degree is unchanged and the clone's degree equals its twin's (cap still
met), and the incidence rises while `R_m` depends only on `m`.
**Admissibility is preserved unconditionally and the defect cannot move.**

### Result 2 — but it is not objective-free: it flips the sign

Collapsing `n8tiny`'s duplicate columns leaves an admissible five-coordinate
family, same size 15, same defect `4/9`, incidence `30` exactly meeting both
`R_15 = 30` and the ceiling `5*6 = 30`. Its exact objective is **positive**:

```text
A_+(collapsed, n=5) = +0.010695694     A_+(n8tiny, n=8) = -0.000478465
```

So the clones, not the ground set, make the record low-defect witness a witness.
The earlier framing had this backwards. For `n8clone_lo` and `n8clone_hi` the
original story does hold — they collapse to *inadmissible* seven-coordinate
families, incidence 196 against `R_70 = 215`, with the cap satisfied — so the two
uses of a clone are genuinely different and both occur among the registered
bases.

### Result 3 — the amplification saturates geometrically

```sh
OMP_NUM_THREADS=1 nice -n 19 python3 -B uc/gate_b/audit_clone_saturation.py \
    --max-clones 4
```

| clones | n | `A_+ <=` | delta | defect | admissible |
|---|---|---|---|---|---|
| 0 | 5 | `+0.010695694` | | `4/9` | yes |
| 1 | 6 | `+0.003478204` | `-0.007217490` | `4/9` | yes |
| 2 | 7 | `+0.000772362` | `-0.002705842` | `4/9` | yes |
| 3 | 8 | `-0.000478465` | `-0.001250827` | `4/9` | yes |
| 4 | 9 | `-0.001137229` | `-0.000658764` | `4/9` | yes |

Delta ratios `0.375, 0.462, 0.527` (and `0.577` at the fifth clone, measured
separately at `n=10` with 210 orbit representatives). The improvement is
geometric and the total budget is bounded, roughly `0.013` for this family.
Verdict `CLONING_IS_DEFECT_FREE_AND_SATURATES`.

### Interpretation

This answers the original plan's amplification question in its sharpest
available form. *Can the numerator be amplified faster than the closure defect?*
Yes — at **zero** defect cost, since cloning cannot move the defect at all. But
the amplification is finite, so cloning flips a family whose objective is
positive but small and cannot rescue one that is far positive: the lowest-defect
admissible family known at `n=7, m=45` has `A_+ = +0.0847`, an order of magnitude
outside the budget. That is the honest reason the local frontier sits at `4/9`
and not at the combinatorial floor.

The actionable consequence for `n = 9, 10` is concrete: rank low-defect
admissible families by `A_+`, and clone only those already within about `0.01`
of zero. Cloning cheap positives is the highest-yield move available, and
cloning expensive ones is wasted work.

Suite: 87 tests pass, including four new ones pinning that cloning preserves
size, defect and admissibility at every rung; that the collapsed core is
admissible yet positive; that the other two bases' cores are inadmissible for
Reimer but not for the cap; and that every delta is negative and strictly
smaller than the one before.



## Failed or superseded routes

1. **Increase `1/50` to another finite coefficient.** Superseded: Cartesian
   powers defeat every fixed coefficient.
2. **Treat a finite n=7/n=8 maximum as Gate B.** Rejected by scope; finite
   searches remain discovery evidence only.
3. **Assume cross-factor adaptive policies might make `C_plus`
   superadditive.** Falsified by the pointwise backward induction; the other
   factor contributes the same constant to every child and cancels from the
   action.
4. **Interpret "near-UC" as an unstated admissibility condition.** Rejected.
   The authoritative supremum has no such restriction. The stronger local
   question is recorded separately rather than silently substituted.
5. **Float-based exact Reimer threshold.** Numerically harmless on the audited
   range, but replaced by integer arithmetic.
6. **Treat certificate formula strings or identity-derived values as independent
   square evidence.** Rejected by adversarial review; v3 keeps the actual
   product-order calculations separate and labels derived values explicitly.
