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
\(A_+<-0.0136721077321777735\ldots\), consistent with the deliberately looser
dyadic bound.

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

Status: **SECOND INFINITE CONSTRUCTION CERTIFIED.**  The Gate B conclusion no
longer depends on the \(n=7\) candidate or any order-orbit quotient.  The new
checker shares the audited dyadic transcendental primitive with the earlier
checker; this is a combinatorially independent route, not a third arithmetic
implementation.

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
