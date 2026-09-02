# Run record — gate H-MINLINE-DGE3-PROOF — verdict FROZEN-CERTIFIED

Campaign: `20260901T111136Z_99745efb_610bb61dc5a8`, agent `RsPe3dTDge`.
Prereg committed pre-compute at git commit `f3d1b03636bc97fb68640f3b01edef7126ddebb6`
(sha256 `558ba6fe9bb61874ee340fc0111e0e2753e067b0ae31e903ae4b0053432db84c`,
byte-identical copy in `pre_statement.md`, source commit recorded in
`provenance.txt`). All controls were executed IN-RUN after `init`, per
prereg section 0 and Main's steering; the verdict rests only on those
in-run results (`controls_results.json`).

## Pre-commit scratch pass (NON-EVIDENTIAL, PRESERVED, disclosed per prereg s0)

Before the prereg commit the author ran a scratch linear-algebra pass
(P2 circuit-kernel all-nonzero on worked instances, P2 constant-vector
conclusion construction, P3 direct + converse on a directed GF(7)
coincidence, P4 pair formula on four hand-picked instances, an M6
sharpness instance, and a reconciliation pre-check of the 15 RS census
rows). Per Main's steering instruction that disclosure was made IN THE
PREREG (section 0) BEFORE the commit; outputs are non-evidential and the
scratch code is preserved in this run dir as
`scratch_precheck.py` (reconstructed transcript of the pre-commit cells;
byte-exactness of a REPL transcript is not claimed — the in-run
`controls_runner.py` computed everything again under the committed byte
constraints, and THOSE runs are the frozen evidence). Every numeric target
the scratch pass checked was SUPPLIED in the assignment (batch context /
Run-2 record), so no threshold could adapt to unseen data; the prereg
states this explicitly (section 0 item 2). One consequence of the pre-pass
is institutional rather than numeric: the 15th RS row family was found to
involve EMPTY parity checks (the zero-map d_i = 1 clause), which the
committed prereg therefore pinned as "count = N" — that reading is also
the assignment's own P5 wording.

## What was discharged

1. **`theorem_t_dge.md` — P1-P5 ALL DISCHARGED with no gap.** P3
   (no-parallel product lemma, both directions, arbitrary-field functional
   extraction via Steinitz basis extension); P1 (spark(H) = d by induction
   with regrouping B = H_2 (x) ... (x) H_n, B's columns verified nonzero,
   spark by IH, dual-isolation below d with the strict inequality r <= k <
   d, plus the hyperplane transport for the upper bound); P2 (Case 1
   r < d_A => single group => fiber along B with d_B = d circuit; Case 2
   r = d: subcase 2a impossible via functional separation; subcase 2b =>
   circuit 1-dim kernel everywhere-nonzero alpha, NO zero-coefficient case
   split because 0 lies in span(alpha), all c_k b_{v_k}/alpha_k equal one
   nonzero vector w, so all b_{v_j} parallel; P3' then forces all v_j
   equal whenever B has no distinct parallel pair); P4 (regime lemma
   completing (B) and (C), pair count (prod A_i - N)/2 with
   A_i = sum m^2, fiber subcount with the combinatorial identity
   A_i = s_i + 2 sum binom(m,2), fiber-regime carrier count with unique
   witness axis); P5 (generalized-Vandermonde spark s_i - t_i + 1, MDS
   circuit counts binom(s_i,d), class structures per d_i, the LINE
   corollary stated ONLY for t_i = 1, Run-1 reconciliation 24 = 16 + 8,
   the t_i > 1 line-failure example (4,4,4),t=(2,1,1) at 64 circuit
   carriers and zero line carriers, the zero-map d_i = 1 clause, and the
   Lambda-invariance note at Id).
2. **M1 6/6 EXACT** — [(3,5),(3,5)]->100, [(3,5),(4,5)]->50, [(3,4)]^3->192,
   [(3,4),(4,4),(4,4)]->64, [(4,6)]^2->180, [(4,6),(5,6)]->90; every one
   with off = 0, below-d emptiness exhaustively asserted, spark(product)
   asserted = d, and the P4 circuit-count formula asserted EXACTLY.
3. **M2 7/7 EXACT against the P4 pair formula with MEASURED parallel-class
   multiplicities** (prereg-pinned A_i = sum m^2 semantics): all-distinct
   [2,2] classes N=25: 300 pairs (10 fibers... measured 100 fibers + 200
   off); Main's pinned anchor [2,2](4,4): 120 = 48 fibers + 72 off
   (Run-2's closing-record reading reproduced EXACTLY); [2,2]+clean
   N=100: 1200 = 400 + 800; single-degenerate rows 50/0, 96/0, all
   fibers; MULTI-class spark-2 rows 12/0 (off = 0, single-degenerate
   regime) and 66 = 30 + 36 (multi x one, both degenerate).
4. **M3 21/21 randomized non-MDS GF(7) configs, ZERO violations** (55
   draws; draws with an infinite-spark (square full-rank) factor are
   resampled and unrecorded — disclosed in-run): spark(product) = min
   spark everywhere; 12 fiber-regime rows with off = 0; 9 diagonal-regime
   rows each satisfying the P4 pair formula with measured multiplicities
   AND exhibiting >= 1 off-fiber diagonal (M6 sharpness inside M3).
5. **M4 plants both REJECTED, fail-loud, pristine ACCEPT:** non-tensor
   column plant (flat-0 column := sum of product columns at flats 13 and
   21) — the planted triple {0,13,21} is an off-fiber dependent d-set
   (witness recorded), a regime-1 violation the machinery refuses;
   duplicated factor column plant — measured spark 2 vs registered 3,
   config mismatch fires; the pristine anchor passes the identical census
   path.
6. **M5 15/15 EXACT** — all 13 (q,s,t) reconstructions plus the two
   trivial zero-map rows (count N), and the Run-1 anchor
   (13,(2,2,4),t=(1,1,1)): 24 = 16 fibers + 8 diagonals with the 8
   off-fiber pairs enumerated and the flat-(0,12) witness pair
   REPRODUCED under the pinned sorted-axis convention.
7. **M6 boundary sharpness 6/6** — every 2-degenerate configuration has
   >= 1 off-fiber diagonal (witnesses recorded: e.g. [0,6], [0,5], [0,24]);
   every single-degenerate configuration has EXACTLY 0 diagonals. The
   hypothesis of (B) is sharp in both directions.

## Inherited-defect note (carried from Run 2, disclosed in prereg s0)

Run 2's committed prereg stated the regime-2 count as
$\sum_{\emptyset \ne D} 2^{|D|-1} \prod_{i \notin D} s_i$, which
implicitly assumes single-class (fully proportional) degenerate columns
and UNDERCOUNTS multi-class spark-2 factors (e.g. it gives 10 vs the true
120 at [2,2](4,4), and 20 vs 300 on the all-distinct 5-column row). This
run pinned Run 2's own closing-record semantics — unordered distinct
dependent pairs = (prod A_i - N)/2 with A_i = sum_k m_{i,k}^2 — which
sibling `RsPe3dMinline` confirmed by IRC before the prereg commit and Main
pinned in the batch context. Run 2's prereg BYTES are frozen and were not
modified; the correction lives here and in this run's prereg only.

## Instrument defects during THIS run (all pre-freeze, all disclosed)

- D-1 (caught by first in-run execution): `is_fiber` and the
  fiber-census circuit test initially built the r x k transpose instead
  of the k column vectors, flipping dependence detection for r < k
  submatrices (M1 halted with the formula mismatch 100 vs 0). Root cause:
  `columns_dependent` expects a LIST of coordinate vectors; call sites
  passed a row-matrix. Fixed by constructing true column vectors at both
  call sites; all six families re-run from scratch afterward, PASS.
- D-2 (caught by execution): M3's first battery version called
  `min(sparks)` on a draw containing a square full-rank factor whose
  spark is infinite (None). Resolution: such draws are RESAMPLED and not
  recorded (disclosed in-run and above); the 21 recorded configs all have
  finite measured sparks.
- D-3 (caught by execution): the M4 non-tensor plant first replaced a
  column by the sum of a NEIGHBOR pair, which produces NO below-d or
  off-fiber signature (the plant was invisible — an instrument-design
  defect, not a theorem defect). Re-designed per the prereg's intent to
  plant the sum of two OTHER product columns, making the planted triple
  {0,13,21} an exhaustive, witnessable off-fiber carrier; verified the
  triple dependence explicitly before re-running. All families re-run
  from scratch afterward, PASS.
- D-4 (caught by execution): M5 originally measured sparks on EMPTY
  parity checks (s_i = t_i zero-map rows) and crashed on the empty
  column list. Resolution: the zero-map clause is detected before
  measurement and counted as N per the prereg. Also the M2/M6 ndeg book
  keeping initially read a stale dict key, misclassifyng row 'd';
  fixed to measured sparks. All families re-run from scratch, PASS.
- No floating point anywhere; PYTHONDONTWRITEBYTECODE=1 was set for
  every execution (no __pycache__ in the run dir at freeze).

## Verdict rationale

The prereg's fixed rule: FROZEN-CERTIFIED iff P1-P5 all discharged with
no gap AND all control families pass exactly in-run. Both arms hold:

- `theorem_t_dge.md` discharges every obligation; the six audit points
  named in the batch context (existence of dual-isolation functionals
  from linear independence over an arbitrary field; the r = d step's
  1-dimensional everywhere-nonzero circuit kernel; the no-case-split role
  of 0 in span(alpha); separation of functionals over an arbitrary field;
  the regrouping of B = H_2 (x) ... (x) H_n with nonzero columns and the
  IH spark; the d = 2 single-degenerate case via P3) are each written out
  explicitly in the artifact.
- M1: 6/6 exact; M2: 7/7 exact (measured multiplicities); M3: 21/21, zero
  violations; M4: plants rejected + pristine accepted; M5: 15/15 exact;
  M6: sharp 6/6. `controls_results.json` embeds verdict_fields with
  all_families_pass = true.

Therefore **Theorem T-DGE is promoted: FROZEN-CERTIFIED**, with the LINE
corollary of P5 stated only for t_i = 1 exactly as pre-registered. The
Run-1 FROZEN-NEGATIVE verdict (line form falsified on
(13,(2,2,4),t=(1,1,1))) stands, and is EXPLAINED by this theorem's
regime-(C): two tied degenerate axes produce 8 genuine off-line
diagonals — the falsified claim was the line phrasing in regime 2, not
the circuit theorem.

## Scope

Per the prereg's scope sentence: exact censuses at the registered
configs/rows only; no weight->d claims; Lambda tested at Id; non-Vandermonde
Gospels and Conjecture 4.2's rho content untouched.
