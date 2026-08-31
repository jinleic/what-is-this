# Gate-A pre-statement: zero-level CCZ

**Target.** Reproduce or refute the zero-level CCZ-distillation headline
constant from arXiv:2605.21867 (Itogawa et al., May 2026) with a bounded,
single-core Stim experiment.  The paper's Fig. 10 reports a least-squares
headline `p_L ~= 300 p^2` for the zero-level construction.  This statement is
committed before the first campaign run.

## Exact quantity and pass/fail rule

For each `p` in `{5e-4, 1e-3, 3e-3, 1e-2}`, run at least `100000` independent
shots of the committed circuit.  A shot is **accepted** iff every explicit
verification/syndrome detector is zero.  An accepted shot is a logical failure
iff at least one of the three output-X frame observables is one.  Estimate
`p_L = logical_failures / accepted` and fit positive points to

```text
log(p_L) = log(c) + alpha*log(p).
```

**Gate-A PASS:** `150 <= c <= 600` (the requested factor-of-two window around
300), with the fitted `alpha` reported but not used as an additional pass
condition. **Gate-A REFUTE for this reconstruction:** a completed grid with
`c < 150` or `c > 600`.  If zero failures leave too few positive points, the
result is `BENCHMARK`/inconclusive rather than a fabricated coefficient.
Neither outcome is a universal refutation of the paper: the decoder and
surface-code boundary are not specified in the paper text, so this gate tests
the executable model below.

## Committed executable construction

* **Data block:** physical wires `0..7` are cube vertices labelled by
  `i = 4x+2y+z`.  The [[8,3,2]] stabilizers and logicals are the supports in
  `src/code_832.py`.  Prepare `|0>` on `{0,2,5,7}` and `|+>` on `{1,3,4,6}`.
* **Encoder:** use the three rendered Fig. 4(a) layers from
  `src/encoder.py`: `(4->0, 6->7, 1->5, 3->2)`, then
  `(0->2, 7->5, 3->1, 4->6)`, then `(5->4, 1->0)`.  This is ten visible
  CNOT symbols.  The prose calls it an eight-CNOT construction; that
  text/figure mismatch is retained as `[REPORTED]`, not silently resolved.
  The ideal output is checked against the eight commuting `+1` stabilizers
  (five code stabilizers plus three logical-X generators).
* **Fault-tolerant preparation:** wires `8..15` are initialized to
  `|+>^8`; apply data-to-ancilla transversal CNOT, reverse the encoder on the
  ancilla block, measure its four Fig. 4(a) `|+>` positions `{1,3,4,6}` in X,
  discard any `-1`, reset those four wires to `|+>`, re-encode, and apply the
  second data-to-ancilla transversal CNOT.  This follows Sec. III.1's
  double-check description.
* **CCZ layer:** use identity in place of
  `T0 T1^dag T2^dag T3 T4^dag T5 T6 T7^dag`, as Sec. IV does for its Clifford
  Stim estimate.  Therefore this is not a direct non-Clifford CCZ simulation.
* **Syndrome:** wires `16..19` are syndrome ancillas.  The committed circuit
  records a four-check scaffold (one `X^8`, then the three weight-four Z
  checks) and repeats it three times.  It includes a four-wire GHZ preparation
  before each explicit parity extraction; because Appendix A has no
  machine-readable gate list, the reset-to-parity boundary is `[INFERENCE]`.
  Every ideal check is +1 and is a detector.
* **AIT/teleportation boundary:** wires `20,21` are reusable parity meters.
  For logical-Z pairs `(0,4)`, `(0,2)`, `(0,1)`, measure the three-body parity
  `Z_source_a Z_source_b Z_output_k` with one meter, where output
  representatives are reused wires `8,9,10`.  These three surgery outcomes
  are random and are not postselected.  Measure all source data wires in X;
  the three output observables are `output_k-X xor source-LX_k-X`, the ideal
  teleportation frame relation.  The unexpanded d=3/d=7 surface patches and
  this wire reuse are `[INFERENCE]` model boundaries, not claims about the
  paper's physical patch decoder.
* **Noise:** after each one-qubit unitary and each idle wire, apply
  `DEPOLARIZE1(p)`; after each CNOT apply `DEPOLARIZE2(p)`; apply a basis-flip
  with probability `p` at reset and measurement.  This is the Sec. IV model
  as transcribed; the implementation makes idles explicit on all inactive
  core wires at each tick.

## Evidence and falsification

`src/smoke.py` produces a committed campaign snapshot containing the exact
Stim text, metadata, source hash, ideal determinism check, postselection probe,
and a 1000-shot `p=1e-3` pilot.  `src/run_campaign.py` produces immutable grid
snapshots.  Quantitative rows are labelled `NUMERICAL`; ideal static checks are
labelled `REPRODUCED`.  The paper gives no named decoder in the body or
Appendix A; this model therefore reports detector postselection plus the
explicit X-frame parity, and does not label it as the paper's decoder.

## Later gates (not part of Gate-A PASS)

* **Gate B sketch:** arXiv:2606.27358 describes a sqrt(T) (`Z^(1/8)`)
  catalyst, a nine-qubit brickwork Clifford circuit with eight controlled-CNOT
  interactions, and nine d=3 rotated patches grown to d=7; its abstract reports
  leakage near `1e-6` at `p=1e-3` after about seven expected attempts.  A later
  gate must first pin the full patch schedule and leakage decoder.
* **Gate C sketch:** the parent README tracks arXiv:2512.13908 separately.  No
  Gate-C implementation is included in this Gate-A campaign.
