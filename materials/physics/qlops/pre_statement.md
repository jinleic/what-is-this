# pre_statement.md — `qlops/` gates A and B

Written before the first campaign run. Source paper read first-hand (arXiv HTML v2,
timestamp 22 Apr 2026 v2 line on the PDF page; DOI 10.1145/3797968, ACM TQC 7(2):14).
Authors: Kong, Zhang, Chen (Zhongguancun Lab / Tsinghua).

## Gate A — reproduce

**Quantity set** (each reproduced by pure arithmetic from the paper's own stated
inputs; full list with references in README "Transcription"):

1. Eq. (1): `p0 = 1-(1-pL)^(1/(k*d))` for the six GB codes of Table 4 -> Table 5
   `p0` column (6 numbers).
2. Eq. (2): `Q = k/((ceil(tr/tSEC)+d)*tSEC)` for the six GB rows of Table 5
   (neutral-atom QLOPS, 6 numbers) and the QLOPS-density column (12 numbers).
3. Table 5 superconducting rows: QLOPS and QLOPS density at matched p0, computed
   with Table 1 SEC times (0.86 us current / 0.40 us future, sums verified from
   components) and Table 2 reaction times (12 + 6 = 18 numbers). Requires the
   surface-code distance column as an input (not reproducible from prose; see
   "Known non-reproducible components").
4. Sec. 3.5 RSA-2048 numbers: Q_SC = 4.0314e7 (1411 logical, tSEC 1 us, tr 10 us,
   d=25), data-qubit count 1280*430 + 131*(2*25^2-1) = 714,019, ratio 56.4611;
   Q_atom = 6.8089e6 (6128 logical, tSEC 900 us), 0.7626 density; Eq. (3) ratio
   5.244; Eq. (4) ratio 2.4204; underestimation ratios ~110 / ~270.
5. Table 6 structural arithmetic: unit-qubit formula
   `2*(dX+4*dZ)*3*dX+4*dm` (13 distinct protocols, 18 rows), `Total = unit * n_units`
   integer factorization (18 rows), and inversion `cycles = 6*dm/(1-p_fail)` ->
   implied `p_fail` per row.

**Tolerance.** The paper prints no tolerance. Definition used: a number is
*PROVED-reproduced* if |computed/published - 1| <= 5% (relative). Justification:
every quantity above is arithmetic on constants the paper prints to 4-6 significant
digits; Monte-Carlo sampling error (pL measured from ~1e5+ shots) does not enter
because pL itself is taken as the published input. Expected deviations are
rounding-boundary artifacts of the `ceil` in Eq. (2) (observed in hand-checks:
<= 0.35%). A deviation > 5% therefore indicates a transcription/typo discrepancy in
the paper or in my reading, and must be recorded as a finding, per-number:

- reproduced within 5%: **PROVED-reproduced** (with exact delta in the artifact)
- not reproducible from stated inputs (missing constant / unstated convention):
  **NOT-REPRODUCED (blocked-on-missing-input)** -- recorded precisely, not guessed
- reproduced only after a correction the paper's own data implies: **PROVED-reproduced [with correction]**, correction documented.

**Gate A passes** if every reproducible number is PROVED-reproduced (the
non-reproducible set is reported separately); **fails** if any reproducible number
deviates > 5% and no documented correction resolves it.

**Known non-reproducible-from-prose components** (identified before coding):
- Fig. 1 exponential-fit constants: not printed in prose; only recoverable as
  feasibility ranges from Table 5's published distances.
- Table 5 "Distance" column (surface-code distance matching each p0): requires the
  paper's surface-code simulations; treated as input, consistency-checked.
- Table 6 "Distillation error" and "Number of cycles": produced by the paper's
  5-qubit density-matrix simulation under platform noise inputs that are not
  printed; cycles are checked structurally (6*dm/(1-p_fail) inversion) and p_out
  column is tagged [REPORTED].
- Table 6 "Total number of qubits" unit-count convention ("syndrome extraction
  cycle difference and post-selection rate"): not fully specified; checked as
  integer factorization.

## Gate B — sensitivity

**Axis 1 (decoder latency).** Multiply every reaction time `tr` (Tables 2a, 2b, 4)
by m in {0.5, 0.75, 1.0, 1.5, 2.0}; recompute every Eq. (2) QLOPS and density with
all other inputs fixed. Falsifier: any output shifts by more than 2x vs baseline
=> cross-paper comparability falsified on this axis.
Analytic note (to be confirmed numerically): for m<=2,
Q(m*tr)/Q(tr) = (ceil(m*r)+d)/(ceil(r)+d) < (2r+1+d)/(r+d) <= 2 (+epsilon), so a
> 2x shift is structurally impossible for Eq. (2) at fixed d; expected verdict:
robust, with ceil-boundary steps quantified.

**Axis 2 (magic-state protocol swap).** Baseline: Litinski 15-to-1 one-level
constants as used by the paper (Table 6: unit qubits 810..25,098, cycles 18.6..90.6,
per-protocol p_out). Replacement: zero-level CCZ constants of arXiv:2605.21867
(abstract-level, read 2026-08-29): c ~ 300 (p_L = 300 p^2), 22 physical qubits,
circuit depth 24, 3 logical qubits, space-time 22*24 = 528 qubit*cycles per CCZ.
These are tagged **[REPORTED]** until `../msd/` gate A reproduces c ~ 300 within 2x.
Test: under the paper's own rule ("infidelity of the output magic states smaller or
equal to the logical error rate" p_out <= p0), does the swapped protocol satisfy
p_out <= p0 for each of the six GB targets at the platform's physical error rates?
Resource axis: per-magic-state space-time factor 15-to-1 vs zero-level; falsifier:
output (resource or feasibility) shifts > 2x => comparability falsified on this
axis. Note the direction is expected to be opposite on the two axes (resources
shrink ~30x, accuracy limit worsens), which itself is the interesting result.

**Labels.** PROVED / CONDITIONAL / NUMERICAL / [REPRODUCED] / [DERIVED] /
[REPORTED] / [INFERENCE] as per repo contract. All runs low-priority, single core,
bounded; computation here is seconds of arithmetic.
