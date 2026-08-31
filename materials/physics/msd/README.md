# `msd/` — zero-level CCZ reproduction target

**Status: SCOPE-LIMITED-NON-TEST for the paper's `c≈300·p²` claim. The floor
is now explained at mechanism level and the output-patch unblock route is
CLOSED by construction (2026-08-30). Verdicts live in `../RESULTS.md`.**

This target owns a Clifford-only Stim reconstruction of the zero-level CCZ
distillation proposal, with explicit boundaries where the papers do not
provide a machine-readable surface-code circuit or decoder.

### Why the reconstruction cannot test `c·p²` (settled 2026-08-30)

Each observable is built at `src/gate_a.py:237-239` as
`output_x_records[k]` XOR the **raw** `source_x_records[q]` over `LX[k]`,
i.e. a bare parity of the distance-2 [[8,3,2]] source block. Enumerating
the unpatched DEM yields exactly **7 weight-1 undetectable logical
mechanisms** — the 7 nonzero observable subsets, one per source wire 0..6.
Their summed probability reproduces the measured floor: 0.038781 vs
0.038356 measured at p=1e-4 (1.1%), 0.114423 vs 0.110333 at 3e-4,
0.360266 vs 0.319565 at 1e-3. So the "logical error rate" of this
surrogate **is** bare single-fault mass.

Fault locations of those mechanisms split **522 source-DATA / 775
output-rep / 167 other**, so armoring only the output side can shift the
coefficient by ≈2× at most and can never change the exponent. The
attempted expanded-patch build confirmed the prediction: d=3 slope
**0.9902** (c 957) and 2.7× worse than base — it removed weight-1 readout
channels but added anchorless-frame ones; d=7 is worse still (slope
0.474). Code kept as a negative control in `src/expanded.py` +
`src/run_expanded.py` (default behavior unchanged; circuit hashes match
the frozen campaigns).

The only remaining unblock is fault-tolerant readout of the **source**
frame — a full rebuild whose decoder and scheduling would be [INFERENCE]
against a paper that names no decoder. Even a future slope→2 result would
bind only this reconstruction.

## Gate A — refute or confirm

Target: arXiv:2605.21867, *Zero-level CCZ Distillation* (Itogawa et al.,
May 2026). Its Fig. 10 reports `[REPORTED]` `p_L ~= 300 p^2`; the paper's
Sec. IV describes single-qubit gates/idles followed by `DEPOLARIZE1(p)`,
two-qubit gates followed by `DEPOLARIZE2(p)`, and reset/measurement flips with
probability `p`. The body does not name a decoder. It reports `10^7`–`10^8`
trials at six values in the `10^-4`–`10^-3` range, while the rendered Fig. 10
axis extends farther; this text/figure grid discrepancy is retained as
`[REPORTED]`.

The committed pre-statement requires at least `100000` shots at each
`p ∈ {5e-4, 1e-3, 3e-3, 1e-2}`. A shot is accepted iff all explicit
verification/syndrome detectors are zero. An accepted shot fails iff any of
the three output-X frame observables is one. Fit positive points to
`log(p_L) = log(c) + alpha log(p)`. **PASS** is `150 <= c <= 600`; outside
that interval is **REFUTE_FOR_THIS_RECONSTRUCTION**. Too few positive points
is inconclusive, not a zero-error claim. See `pre_statement.md` for the
full one-page quantity, construction, and falsification scope.
The original grid is frozen in
`campaigns/20260829T233635Z_d9fe3c0d_1938285066db/`. Its accepted-shot
logical-error rates are `[NUMERICAL]` saturation values rather than an
asymptotic window: the fitted slope is `alpha=0.553` and the fitted
`c=13.6`. Consequently its frozen machine verdict string remains
`REFUTE_FOR_THIS_RECONSTRUCTION`, but the scientific interpretation is
**UNTESTED** for the paper's `c ~= 300` claim, not refuted and not confirmed.
Gate-A-prime is frozen in
`campaigns/20260829T234855Z_02ebe2db_40d024611d98/`. Revision 2 corrects the
literal `300 p^2` implication at `p=1e-4` to target `p_L=3e-6`, with
factor-of-two band `[1.5e-6,6e-6]`. Revision 1's `[0.0015,0.006]` band
encoded Main's steering arithmetic error; the implementation agent caught it
before any promotion leaf used it. The measured `p_L=0.038356` lies outside
both bands, so the verdict is insensitive to the error. Its corrected
low-p verdict remains `SURROGATE_REFUTATION_AT_LOW_P`; falling
`c_eff=3.84e6 -> 1.23e6 -> 3.20e5` and fitted `alpha=0.9201` indicate
floor-like rather than quadratic behavior at every reachable p
`[NUMERICAL]`. See `pre_statement_rev2.md`.

The circuit is 22 wires: data `0..7`, an eight-wire verification block
`8..15`, syndrome wires `16..19`, and reusable surgery meters `20..21`.
Three verification wires (`8..10`) are reused as effective output
representatives after the second double-check CNOT. The source code block,
the three rendered encoder layers, the double-check, identity T/T-dagger
surrogate, explicit syndrome detectors, and logical-ZZ parity measurements
are in `src/`. The rendered Fig. 4(a) contains ten CNOT symbols although its
nearby prose says eight; the figure-transcribed ten-CNOT sequence is used and
the mismatch is `[REPORTED]`.

The unexpanded d=3/d=7 surface-code patches, three-wire output reuse, and
GHZ-to-parity scheduling boundary are `[INFERENCE]` model choices. The
logical-X frame is explicit rather than a claim about an unnamed paper
decoder. Therefore a future Gate-A verdict applies only to this committed
reconstruction, not universally to arXiv:2605.21867.

## Gate B and Gate C (later)

Gate B tracks arXiv:2606.27358, *Cultivating logical catalysts for
fault-tolerant dyadic phase rotations* (Xu and Wang, June 2026). The abstract
describes a `sqrt(T)=Z^(1/8)` catalyst, a nine-qubit brickwork Clifford
circuit with eight controlled-CNOT interactions, and nine d=3 rotated
patches grown to d=7. It reports `[REPORTED]` leakage near `1e-6` at
`p=1e-3` after about seven expected attempts. Its full patch schedule and
leakage decoder are not implemented here.

Gate C is deferred to arXiv:2512.13908 and will ask whether a pure-Pauli
Clifford model can reproduce that paper's reported cultivation numbers.
No Gate-C claim is made by this target.

## Evidence and inventory

All quantitative run rows are labelled `NUMERICAL`; paper headlines are
labelled `[REPORTED]`; ideal static checks are labelled `REPRODUCED`.
The final smoke snapshot is
`campaigns/20260829T234326Z_8fb8e8ec_4e77386cfa1b/` and contains summary,
the exact ideal Stim circuit, and the `p=1e-3` noisy Stim circuit. Its
256-shot encoder and full-circuit ideal checks have zero violations
`[REPRODUCED]`; the 1,000-shot `p=1e-3` pilot has 454 accepted shots and 132
accepted logical failures `[NUMERICAL]`. The corrected Gate-A-prime snapshot
is `campaigns/20260829T234855Z_02ebe2db_40d024611d98/`; its p=`1e-4` row is
the low-p surrogate-refutation datum described above `[NUMERICAL]`.

* `pre_statement.md` — committed quantity, gates, and falsification rule.
* `pre_statement_rev1.md` — Gate-A-prime asymptotic-window protocol.
* `pre_statement_rev2.md` — corrected Gate-A-prime arithmetic and verdict.
* `src/code_832.py` — code stabilizer/logical/CCZ constants.
* `src/layout.py` — 22-wire role map and figure-transcribed encoder layers.
* `src/encoder.py` — exact input basis and encoder/inverse helpers.
* `src/gate_a.py` — Stim circuit builder and circuit-level noise injection.
* `src/estimate.py` — detector postselection estimator and log-log fit.
* `src/smoke.py` — static, postselection, and 1,000-shot pilot harness.
* `src/run_campaign.py` — finite p-grid runner and immutable artifact writer.
* `campaigns/` — frozen run snapshots; `scratch/` — non-authoritative search.

## How to run

From `physics/` (the prepared `.venv` is required):

```bash
nice -n 10 .venv/bin/python msd/src/smoke.py --pilot-shots 1000
nice -n 10 .venv/bin/python msd/src/run_campaign.py --shots 100000
nice -n 10 .venv/bin/python msd/src/run_campaign.py --label gate_a_prime \
  --shots 10000000 --p-values 1e-4,3e-4,1e-3
```

The default campaign stops after exactly four p points and 100,000 shots per
point. The Gate-A-prime command stops after exactly three p points and
10,000,000 shots per point. Use `--artifact-dir` only with a new empty
directory; snapshots are immutable and contain `summary.json` plus one Stim
circuit per noise point.

## Layout

* `campaigns/` — immutable run snapshots named
  `<UTC-timestamp>_<uuid>_<circuit-hash>`.
* `scratch/` — exploratory files, not evidence.
