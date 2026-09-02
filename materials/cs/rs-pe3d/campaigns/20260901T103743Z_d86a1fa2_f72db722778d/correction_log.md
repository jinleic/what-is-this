# Correction log — Run 1 (gate H-MINLINE) — campaigns/20260901T103743Z_d86a1fa2_f72db722778d

All corrections are logged BEFORE freeze; the frozen artifacts carry the
corrected state. Frozen prereg bytes were NOT modified (prereg = the committed
`cs/rs-pe3d/prereg/H_MINLINE_2026-09-01.md`, sha256
`04294bc3e4757c61975e4c495c1fcf542839856575b3ca11bf74e128c62682cc`, commit
`66fda27`); every correction below is recorded with its reason and evidence.

## C-1 (pre-compute) — order-4 subgroup list in the pinned prereg text

The prereg's "PINNED axes" paragraph describes the order-4 subgroup via the
doubled generator 2 and lists `2^(2k) mod 13 for k=0..3 = {1,4,3,12}` and a
sorted list `[1,3,4,12]`. Both element lists are WRONG as written: 1, 4, 3, 12
are not closed under multiplication mod 13 (4·3=12, 12·3=36=10∉set) and
`[1,3,4,12]` is not a group. The unique order-4 subgroup of the cyclic
F_13^* (2 is a primitive root, 2^6 = 64 = 12 has order 2, so the squares of
F_13^* have order 6 — the order-4 subgroup is <2^3> = <8>) is

    {1, 5, 8, 12}   (5^2 = 25 = 12; 5^3 = 60 = 8; 5^4 = 40 = 1)

sorted ascending `[1, 5, 8, 12]`, which is what every instrument uses
(`subgroup(13,4,...)` in both the run instrument and the frozen `src/field.py`
agree: `[1, 5, 8, 12]`).

Why the instance is unaffected: a cyclic group has exactly ONE subgroup of
each order, so "the order-4 subgroup" is unambiguous; no census number, no
classification, and no witness depends on the mis-listed elements. The
flat-(0,12) mapping is convention-driven (row-major over sorted axis lists)
and unchanged.

Evidence: `python -c "from src.field import subgroup, primitive_root; print(subgroup(13,4,primitive_root(13)))"`
→ `[1, 5, 8, 12]`; the run instrument asserts
`E0.axes == [[1,12],[1,12],[1,5,8,12]]` at startup.

## C-1b (in-run ROOT-CAUSE, route-A semantics) — kron-Vandermonde rows are
NOT V's annihilator; true-annihilator route substituted; all counts reproduce

The prereg's route A described V as ker(H_0 ⊗ H_1 ⊗ H_2) with H_i the GRS
parity check of C_i, and the instrument built H exactly that way. Post-run
audit (triggered by the witness-membership anomaly behind C-2) found:
rowspace(kron) = C0^⊥ ⊗ C1^⊥ ⊗ C2^⊥ (dim (s0−1)(s1−1)(s2−1) = 3 at E0) is
NOT V^⊥: pairing a lifted full-line indicator (an honest V word) against the
kron constancy row gives 2 ≠ 0 mod 13. The TRUE annihilator at t = (1,1,1)
is V^⊥ = ∩_i L_i^⊥, also dim 3, but with DIFFERENT rows — per-axis-2-slice
difference patterns [[12,1],[1,12]] on each (a0,a1) plate. ker(kron) = {v :
axis-2 slice sums constant} is a DIFFERENT 13-dim space from V. The two
spaces nevertheless have dim(X ∩ F^S) EQUALING on every support of weight
≤ 4 of all three instances (exhaustive 120 + 276 + 7140, plus 800 random E0
spot checks at weights 3–4) — a dimensional coincidence recorded here and
NOT relied upon.

RESOLUTION: route A was REBUILT as the space-exact dual (annihilator rows =
independent nullspace of the frozen-src basis route G, per instance), all
three censuses RE-RUN with per-support equality asserted against the basis
route: identical results E0 24/16/8, C1 12/0, C2 24/0, same 8 off-line pairs,
same witness word in the kernel. Machine artifact of the superseding run:
`census_trueann_run1.json`; census counts are frozen from THAT run.
`census_runner.py` kron path retained for the record, output superseded.
Classification: instrument defect (route-A semantics), caught by cross-check
before any frozen byte; no pre-registered numeric target changed.

## C-2 (in-run, caught by cross-check) — witnesses.json hand-transcription defects

Two successive hand-transcription errors were made while writing
`witnesses.json` by hand (wrong summand supports/values, incl. a spurious
entry at flat 4; and later flat-12 value 11 instead of 12). Both were caught
by the independent second-convention cross-check and by
`record_crosscheck.py` (which asserts every recorded word lies in the kernel
of an independently built check matrix and re-censuses all three instances).
FINAL STATE: `witnesses.json` was REGENERATED PROGRAMMATICALLY from machine
artifacts (no hand-transcribed numbers remain), and the FULLY machine-generated
`census_trueann_run1.json` (from `census_trueann_runner.py`) is the primary
witness record for freeze. There the summands are

    sum0 (over flats): [0,0,0,0, 12,0,0,0, 0,0,0,0, 12,0,0,0], support {4,12}
      = (axis-0-constant) (x) h (x) delta_1, h = (0,-1) over axis-1 list
    sum1: [1,0,0,0, 1,0,0,0, 0,...], support {0,4}
      = g (x) (axis-1-constant) (x) delta_1, g = (1,0) over axis-0 list
    sum = sum0 + sum1, support {0,12}, values (1, 12)

with flat 4 = (1,12,1) cancelling EXACTLY (12+1 = 13 = 0); all three words
pair to zero against the TRUE annihilator (C-1b), and membership was
independently re-confirmed by rank(G + [word]) == rank(G).

## C-3 (cosmetic) — prereg said "anti-diagonal", reality is the plate DIAGONAL

The prereg's pre-registered expectation worded the witness plane matrix as
"[[1,0],[0,-1]] ... ANTI-diagonal pair". Under the pinned convention the
witness lies on the x2 = 1 plate with matrix (rows = axis-0 index over
[1,12], cols = axis-1 index over [1,12], values × δ_1):

        x1 = 1   x1 = 12
    x0=1:    1        0
    x0=12:   0       12 (= -1)

i.e. the DIAGONAL of that plate, support {(1,1,1), (12,12,1)} = flats {0,12}.
The "anti-diagonal" wording was the transposed-reading slip; the substantive
pre-registered expectation — reproduction of the prior census's flat-(0,12)
pair — holds exactly. No numeric target affected. Superseded detail: the
primary record is now `census_trueann_run1.json` (C-1b).
## C-5 (in-run, caught by record cross-check) — witnesses.json summand
transcription defect (RD-2) and missing JSON closer (RD-1)

The machine witness in `census_hminline_run1.json` (`witness.summands`,
unchanged since the single census run) is:

    sum0 = 1 ⊗ h ⊗ δ_1  (h = (0,-1) on axis-1 list [1,12]), support {4, 12}
    sum1 = g ⊗ 1 ⊗ δ_1  (g = (1, 0) on axis-0 list [1,12]), support {0,  4}
    sum  = sum0 + sum1, support {0, 12}, values (flat 0: 1, flat 12: 12)

with flat 4 = (1,12,1) cancelling exactly (12 + 1 = 13 ≡ 0). The first
hand-transcription in `witnesses.json` wrote sum0 support {12} only and sum1
support {0,4} with 1 at flat 12 — both WRONG; the corrected record (written
from the reconstructed machine words) carries:

    sum0 : [0,0,0,0,12,0,0,0, 0,0,0,0,12,0,0,0]
    sum1 : [1,0,0,0, 1,0,0,0, 0,0,0,0, 0,0,0,0]
    sum  : [1,0,0,0, 0,0,0,0, 0,0,0,0,12,0,0,0]

and passes the independent `record_crosscheck.py` (all three words in ker(H)
under the kron check — see C-1b for why that kernel COINCIDES with V here at
weight ≤ 4 — pairwise sum exact, axis-constancy of both summands verified,
full 24/16/8 + 12/0 + 24/0 recount reproduced). The record ALSO initially
lacked one closing brace (RD-1) — caught when the cross-check failed to
parse; both defects are hand-transcription defects ONLY: the machine
artifacts were correct throughout. SUPERSEDED: the primary machine witness
record for freeze is `census_trueann_run1.json`; see C-6.

## C-6 — FINAL supersession: census_trueann_run1.json is the frozen claim

Per C-1b, all three censuses were re-run with the space-exact dual routes
(`census_trueann_runner.py`, asserting per-support equality of the
true-annihilator and basis routes): results identical to every prior run —
E0 24/16/8, C1 12/0, C2 24/0; same 8 off-line pairs; same witness words,
now certified against the TRUE annihilator. `witnesses.json` additionally
regenerated programmatically from machine artifacts. FROZEN CLAIM RESTS ON:
census_trueann_run1.json + witnesses.json + record_crosscheck.py + controls
recorded in census_trueann_run1.json:controls_C3.

## C-4 — instrument defect record (edit-tool body-row mishaps during authoring)

During authoring of `prereg/H_MINLINE_2026-09-01.md` and `census_runner.py`,
several edit-tool operations produced orphaned/duplicated lines (e.g. a stray
non-ASCII character line in census_runner.py; duplicated lines in the prereg).
All were caught by immediate re-reads, non-ASCII grep, and the final full-file
read-throughs BEFORE the prereg commit and BEFORE the first instrument run;
no wrong number, control, or adjudication clause ever existed in a committed
or executed byte. Listed here for completeness per defect-disclosure policy.
The executed instrument differs from its first written form only in:
(a) removal of the stray character line, (b) the witness proportionality
check using support-position indexing (the first form indexed flat ids into a
length-2 extracted basis — a genuine indexing bug caught by the first run's
traceback and fixed before any data were recorded), (c) removal of orphaned
lines from a partial edit. Same policy for `census_trueann_runner.py` (its
import block needed one repair after an edit-tool slip, caught by an
import traceback before any run). The FINAL run is defect-free end to end.
