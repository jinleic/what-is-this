# Gate B — closure-defect scalar route

## Current status

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
Thus an error specific to the 45-row family or its order-orbit reduction cannot
reverse the Gate B conclusion.

No fixed scalar closure-defect coefficient can repair \(A_+\) on all
cap/Reimer families. This is not a union-closed counterexample and does not
resolve Frankl's conjecture.

The defects in this construction tend to one. A stronger question restricted
to \(\varepsilon_\vee\to0\) remains open, but that restriction was not present
in the authoritative Gate B supremum. See [DEFINITIONS.md](DEFINITIONS.md).

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
| [certificates/gate_b_unbounded_arb_v4.json](certificates/gate_b_unbounded_arb_v4.json) | certified | Authoritative interval output. It verifies the full cell-union row definition, directly certifies the base ratio, separates actual square-order checks from identity-derived values, and records exact power ratios. The v1--v3 JSON files remain intermediate checkpoints. |
| [search_n8_block_symmetric.py](search_n8_block_symmetric.py) | complete finite search | Checkpointed exhaustive \(S_1\times S_7\) class on \(n=8\). Exact class/defects; float64 objective. |
| [candidates/n8_k1_census.json](candidates/n8_k1_census.json) | complete numerical | 65,535 raw masks, 136 admissible canonical families, two negative; maximum ratio \(0.03433549\ldots\), below the certified \(n=7\) base. |
| [experiments/n8_k1_checkpoint.jsonl](experiments/n8_k1_checkpoint.jsonl) | append-only log | All 136 evaluated records; resumable without overwriting prior records. |
| [audit_n5_complete.py](audit_n5_complete.py) | complete finite audit | Resumable independent float64 evaluation of all 5,172 exact nontrivial \(n=5\) cap/Reimer coordinate orbits representing 463,343 labeled families. |
| [candidates/n5_complete_audit.json](candidates/n5_complete_audit.json) | complete numerical | Zero negative or near-zero values; minimum \(A_+=0.006474313148642441\ldots\). |
| [experiments/n5_complete_checkpoint.jsonl](experiments/n5_complete_checkpoint.jsonl) | append-only log | All 5,172 evaluated nontrivial \(n=5\) records; exact filters/defects and float64 objective values. |
| [audit_tensorization.py](audit_tensorization.py) | resumable falsification search | 292 ordered product cases (205 with unequal block dimensions), every global order of each product, attacking fixed-order \(Q\) and \(C_+\) additivity and the exact defect product. |
| [candidates/tensorization_audit.json](candidates/tensorization_audit.json) | adversarial numerical | Zero counterexamples; worst fixed-order gaps \(1.8\times10^{-15}\) (Bellman) and \(8.9\times10^{-16}\) (iid), i.e. float64 rounding. |
| [audit_square_direct.py](audit_square_direct.py) | resumable direct check | Materializes the 2,025-row square as plain 14-bit masks and runs the generic evaluators, with no factored-fiber or product-aware shortcut. |
| [candidates/square_direct_audit.json](candidates/square_direct_audit.json) | direct numerical | Five declared global orders; exact square combinatorics and zero product gap in every order. |
| [audit_power_admissibility.py](audit_power_admissibility.py) | exact constructed-family audit | Builds \(\mathcal F_k\) for \(k\le3\) and recomputes size, counts, cap, incidence, Reimer, normalization, and the missing-join count by brute force over all ordered pairs. |
| [candidates/power_admissibility_audit.json](candidates/power_admissibility_audit.json) | exact | At \(k=3\): 91,125 rows, all counts 36,450 equal to the cap bound, incidence 765,450 against threshold 750,668, and 8,227,000,000 of 8,303,765,625 ordered joins missing, i.e. exactly \(1-(17/81)^3=526528/531441\). |
| [test_gate_b.py](test_gate_b.py) | regression suite | 24 tests: exact constraints, 5,400 fixed-order comparisons over all 225 products of nonempty \(n=2\) families, independent evaluator agreement, dyadic-vs-Arb enclosure dominance, order-orbit value invariance, symbolic Bellman algebra, exact integer cap/Reimer/ratio checks for the general power claim, symmetry/order coverage, and artifact contracts. |
| [EXPERIMENTS.md](EXPERIMENTS.md) | chronological ledger | Commands, parameters, failures, scope labels, and resource use. |
| [paper/main.tex](paper/main.tex) | paper draft | Paper-quality statement, proof, both computational lemmas, controls, and limitations. |

Authoritative SHA-256 values after final review:

```text
405b0b74129c6192be791788b22ffbfbcb921e407e1d0002d01fb8193d131c72  verify_gate_b.py
6784c4d90abf8292a901185504c6727baff8ba8e894526659d9855263ddd0d46  verify_gate_b_dyadic.py
dbe8064eed8d68dc3b49a7823242b599188dfcc82d98e892cf19b332fcc971ea  audit_tensorization.py
a2a9477a71153242726311df9792006a1f614cee6b178a60b9c886fd3425740d  audit_square_direct.py
6f963f4eaf05408a10bbc180854df30aa1052a5771eae3bc1f03000fcd3ef016  audit_power_admissibility.py
a79651c87edd36332fe1d0c4ff3b59b8efaad96f5194121d99a3dbd7e525b3b7  audit_n5_complete.py
365a2ccb1482e7a68f1b4dab10fc813e0167fa83ddf87812ad10af915934612d  search_n8_block_symmetric.py
b007967afa6da21ec1c7506bb3d383b0656f2fb37c7825b6bf4500d633b3ab8d  test_gate_b.py
9deabde710eff1e9d49b7c198edecc8d8af4e3d606a9e69feb179eb52df4df36  candidates/n7_block_extremizer.json
14f6875418849011ccc350d5a0d1da7a57eae8bfe695182559184de57d281588  candidates/n5_complete_audit.json
d65d5f5e2735295ae5988cd1077814dc353d8609f7d951d45d150672e25b112d  candidates/n8_k1_census.json
0f48dc396827fad4043f148a552f2b2d1a4bceb7092d6ce4fd3a4b9cef2801da  candidates/tensorization_audit.json
92659ae392a18b1146a1520dc63c80d28d70469aad16f30499dd2cfdb0cfd7cc  candidates/square_direct_audit.json
287ee465e2376f12441641bcc5795463d0b867dd27821135531f43cb8ecdca0a  candidates/power_admissibility_audit.json
9c147dba19b8d3f00055c43917927cb5f70b07daa87d94528fad521564d09882  certificates/gate_b_unbounded_arb_v4.json
d4019a9bca35e0945ab75c05c9a72dcf52311f022429d179f1b193acadef937b  certificates/gate_b_unbounded_dyadic_v1.json
90fa33d4849c568bd686941fc951e11bdc6aebc64aa94202cf4bc6d14f6631c6  experiments/n5_complete_checkpoint.jsonl
012051ad8feb2aaac3f410318957f14e9c9836c9fad7111c6fdc91940d52ed47  experiments/n8_k1_checkpoint.jsonl
2c46e8dae7b8268bd28e10357374617a21e434a3da44da1e8bac7136058cdd9f  experiments/tensorization_audit_checkpoint.jsonl
0599d2b038a82fa7298c4ca0d3f5d8b6d6761451271591e4dac449037b07ac39  experiments/square_direct_checkpoint.jsonl
2695bb03ae5d5a5a222b71201330b6c37f0d2306bdb0db49e35248062d522cb1  experiments/power_admissibility_checkpoint.jsonl
c5caf2c7a5d6a3e04057641d23150262183e230b5ee89910c42aaf27371313a7  PROOF.md
2d295690e1159e60859401ab96433e810d778b409e91afe4edb01dd8ca4f08d8  DEFINITIONS.md
cb41ff8284189575fd17cc6d9015cabf4260405ab5acb9e3db9bb30a4110cb1a  EXPERIMENTS.md
62f550287852662a236a1d65e5d869d63781f9f55483a3da80bbc154d2904f6a  paper/main.tex
31da3136e01f37b4fb2f6ccc3e913e483f12acf6d3abf1208425cd781cfd1298  paper/main.pdf
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

The \(n=8\), \(n=5\), tensorization, and square audits resume from append-only
checkpoints and must report zero new evaluations when the checked-in
checkpoints are complete.

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

## Evidence labels

- **PROVED:** product identities, admissibility of every power, and divergence.
- **CERTIFIED:** exact base constraints and the 256-bit Arb upper bound
  \(A_+(\mathcal B)<-7/250\).
- **INDEPENDENTLY CERTIFIED:** the standard-library dyadic relaxation proves
  \(A_+(\mathcal B)<-1/40\) without Arb or Bellman clamp classification. Its
  exact product consequence gives ratio \(>k/40\to\infty\).
- **INDEPENDENTLY CERTIFIED SECOND BASE:** all 720 orders of the distinct
  25-row \(n=6\) family give the standard-library dyadic relaxation
  \(A_+(\mathcal D)<-1/80\), and exact product arithmetic gives ratio
  \(>k/80\to\infty\).
- **COMPLETE EXACT / COMPLETE NUMERICAL:** the full nontrivial \(n=5\)
  cap/Reimer census,
  the declared \(n=6\) controls, and the declared \(n=8\)
  \(S_1\times S_7\) class. Enumeration, filters, defects, and stated orbit
  coverage are exact; objective comparisons are float64.
- **OPEN:** the newly separated local-stability problem with
  \(\varepsilon_\vee\to0\), and any non-scalar UC-specific charge capable of
  advancing the Frankl programme.
