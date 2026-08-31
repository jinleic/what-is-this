# Pre-statement: exact census gates for contractive-unitary classical shadows

Frozen: 2026-08-29, before any campaign run. Owner: shadows.

Basis. Written after a first-hand read (2026-08-29) of the published open-access
version: Wu, Wang, Yao, Zhai, You, Zhang, *Contractive unitary and classical
shadow tomography*, npj Quantum Information 12, 86 (2026),
DOI 10.1038/s41534-026-01227-w (= arXiv:2412.01850, v1 posted 2024-11-28),
including the full Methods section and the tables/2 landing page (Table 2
content read verbatim: rows [XX YX ZX; XY YY ZY; XZ YZ ZZ]). The follow-up
paper arXiv:2608.18935 (Hingane & Koh, "Real Classical Shadows with Noise",
2026-08-19) was read at abstract level for the follow-up campaign only; no
gate below depends on it.

## Published-side conventions (pinned, first-hand, 2026-08-29)

- Conjugation direction: Methods states the two-qubit anchor as
  U = exp(i pi/4 ZeZ1 Z2) acting as O -> U-dagger O U, transforming
  {X1Z2, Y1Z2, Z1X2, Z1Y2} -> {Y1, -X1, I1Y2, -I1X2} and leaving the five
  remaining size-2 operators' sizes unchanged. Hand re-derivation (class
  level, rule: P anticommutes with G=Z1Z2 => U-dP U has class iP*G, else P):
  X1Z2 -> +Y1, Y1Z2 -> -X1, Z1X2 -> -Y2 (class [Y2] agrees with paper's +Y2;
  Z-scheme sign conventions differ, phase-irrelevant), Z1Y2 -> -X2;
  XX, XY, YX, YY, ZZ fixed. The main-text sign list differs from Methods on
  two signs; signs never enter a size census, so the class-level map is
  authoritative here.
- Symplectic anchor (derived, to be machine-asserted): in (x1,z1,x2,z2)
  GF(2) label coordinates, the class map of U has generator images
  S(X1)=x1+z1+z2 = 1101, S(Z1)=z1 = 0100, S(X2)=x2+z1+z2 = 0111,
  S(Z2)=z2 = 0001 (columns). Symplectic (S^T J S = J) by hand; code MUST
  re-assert: S symplectic AND all four contracted classes + five unchanged
  classes match the Methods list exactly.
- Phase invariance: a size census depends only on Pauli classes (up to the
  4-phases); Clifford unitaries mod global phase =
  (Pauli-class translation p in GF(2)^4) x (symplectic S), conjugation on
  classes is P_c -> S P_c for all 16 translations p. So the exhaustive size
  census lives on Sp(4,2) (720 elements), and the mod-phase Clifford count
  is 720 x 16 = 11,520 = 2^(n^2+2n) prod_j (4^j - 1) at n=2 [DERIVED,
  Appleson-Gross/Hostens-Dehaene standard structure; the paper itself never
  states 11,520].
- Eq. (3) (size after U_ct): m = N_XY if N_XY odd, else m = k, where N_XY
  counts X/Y letters of the size-k input string. Identities inserted as I
  letters become Z letters when N_XY is odd (Methods, identity-insertion
  discussion).
- Eq. (4): w(O)_ct = (1/2)[3^-k + (-1)^k 9^-k] + (1/2)[(5/9)^k - 9^-k];
  dominant asymptote 1/(2 x 1.8^k); shadow norm ||O||^2_ct = 1/w.
- Eq. (5): pi(m)_ct = 3^-k C(k,m) 2^m (1-(-1)^m)/2
           + 3^-k (sum over l of C(k,2l) 2^(2l)) delta_{m,k}.
- Eq. (7) (q inserted identities, size k~ = k-q):
  w = (1/2)[3^-k~ + (-1)^k~ 9^-k~] + (1/2)[(5/9)^k~ - 9^-k~] / 3^(k-k~).
- Fig. 1b (k=50, read verbatim): U_ct ensemble = broad peak near m/k = 2/3
  (odd N_XY) + delta peak at m/k = 1 (even N_XY), total weights of the two
  contributions EQUAL; maximally scrambled operator peaks near m/k = 3/4.
- Lemma (Methods, "Two-qubit contractive unitary", read verbatim): at most 4
  of the 9 size-2 Paulis can be contracted to size 1 by ANY Clifford;
  achievability by exp(i pi/4 ZZ); Table 2 = 3x3 arrangement whose rows and
  columns are pairwise-anticommuting triples; 5 contractions would force 3
  pairwise-anticommuting size-1s on one site, whose preimages fill a
  row/column, and the remaining contracted preimages must commute with that
  triple while avoiding its row/column - impossible in 3x3. [INFERENCE, to be
  noted in the artifact: the same phrasing does not mention the extra
  pairwise-anticommuting diagonal triple {XX, YY, ZZ} of Table 2; that case
  dies by the same contradiction and the exhaustive census below settles all
  cases at once regardless.]

## Gate A -- exhaustive two-qubit contraction census (paper lemma)

- Quantity: enumerate EXACTLY, in GF(2) integer/bit arithmetic (no floating
  point anywhere in the census): (i) all 65536 binary 4x4 matrices, keep
  symplectic ones; expect exactly 720 (hard assert); (ii) all 11,520
  mod-phase Clifford labels (S, p), p in GF(2)^4 (hard assert 11,520);
  (iii) for every one of the 720 class actions S, the multiset
  n(S) = #{of the 9 size-2 classes c : size(S c) <= 1}.
  Deliverable: full histogram of n over all 11,520 labels (each histogram
  bin count must be 16 x bin over S), plus the set of achieving actions
  (expect: includes the U_ct anchor pattern).
- Pass (CONFIRMED): max_S n(S) == 4, |Sp(4,2)| == 720, total labels
  == 11520, and the anchor S validates the paper's exact 4-contracted /
  5-unchanged pattern.
- Fail / falsification: any S with n(S) >= 5 (refutes the lemma's upper
  bound), or no S with n(S) == 4 anchor-equivalent (refutes achievability),
  or any count assertion fails (census itself invalid - reported as error,
  not as refutation).
- Evidence label: [PROVED] if passes - the enumeration is exhaustive and
  exact, so this is a universal statement, not a sampled one.

## Gate B -- exact re-derivation of the Pauli weight formula

All quantities exact via Python fractions.Fraction unless tagged otherwise.
Independent means = computing from the ensemble definition
(uniform token strings in {X,Y,Z}^k + the composed contractive class map
S_ct^(k) built factor-by-factor from the anchor S), NOT from the paper's
formulas.

1. (k <= 10, exhaustive) For every one of the 3^k token strings: applied
   S_ct^(k) per-string size equals the Eq. (3) rule (N_XY odd -> N_XY, even
   -> k). Evidence: [PROVED for k <= 10 by exhaustion].
2. (k <= 10) The exact per-string size histogram (integer counts / 3^k,
   Fractions) equals Eq. (5) pi(m)_ct for every m. [PROVED for k <= 10].
3. (k <= 16) Eq. (5)-summed weight  sum_m pi(m)/3^m  (own Fractions
   evaluation) equals Eq. (4) exactly for every k in 1..16.
   [DERIVED, exact rational identity verified to k=16].
4. (k <= 16) Independent identity-insertion check of Eq. (7): for k <= 7,
   enumerate token strings over {X,Y,Z} on the k~ support with I letters on
   the remaining q = k - k~ sites, evolve with the composed k-site class map
   (identities gain Z exactly per the paper's rule - machine-checked per
   string), histogram sizes as Fractions, assert equals Eq. (7) for every
   (k~, q) with k~ + q <= 7, k~ >= 1. [PROVED for the pinned range].
5. Random-Clifford baseline, re-derived from scratch: uniform non-identity
   image Pauli => pi_rc(m) = C(k,m) 3^m / (4^k - 1); w_rc = sum pi_rc/3^m
   = (2^k - 1)/(4^k - 1) = 1/(2^k + 1) (paper's value) - assert exact
   equality of the Fraction chain for k <= 16. [DERIVED].
6. (NUMERICAL, float only) Asymptote confirmation: emit
   r_k = w_ct^-1 / (2 x 1.8^k) for k = 1..16 and verify monotone -> 1 trend;
   and w_ct/w_rc growth ~ (2/1.8)^k/2. Float comparisons are labeled
   NUMERICAL and never promoted.
7. ([DERIVED]) Fig. 1b shape at exact level: odd-N_XY total mass == 1/2,
   even-part mass == 1/2 (exact Fractions), odd-part mode at
   floor/round(2k/3) for the checked k; random-Clifford pi_rc mode at
   floor(3k/4).

- Pass: all of 1-5 exact-pass and 6 shows the approaching trend.
- Falsification: any exact mismatch in checks 1-5, or r_k not approaching 1
  (trend clearly away across k = 10..16).

## Gate C -- bounded alternative-family search (scaffolded, pinned scope)

Scope pinned now (avoid unbounded search): evaluate, with the SAME exact
machinery (composed class maps + 3^k token enumeration, Fractions),
w_M(k) = E_{uniform token}[3^-size(S_M token)] for k <= 8 over the pinned
families: (a) uniform-letter meshes M_L = prod_{i<j} exp(i pi/4 L_i L_j),
L in {X, Y, Z} (Z-mesh = U_ct must reproduce Eq. (4) - consistency anchor);
(b) all products over perfect matchings of sites with arbitrary letters per
pair (k <= 8: 105 matchings x 3^4 letterings = 8505 composites);
(c) commuting deterministic layers CZ-mesh, CX-chain, and U_ct composed with
Metric: the paper-beating direction is a SMALLER shadow norm, i.e. a
LARGER Pauli weight: PASS if any pinned family attains
w_M > w_ct (equivalently w_M^-1 < w_ct^-1 = the paper's ||O||^2) at some
k <= 8 while reproducing the U_ct anchor, else NOT-FOUND-in-pinned-scope +
note. Explicit non-goal: the full space of
deterministic commuting Clifford families is not enumerated; a NOT-FOUND
verdict is scope-limited, not an obstruction proof.

## Cost budget

Gate A: ~65k matrix checks + 11520 label passes - seconds.
Gate B: 3^10 x (10x10 bit ops) + Fractions - well under 2 minutes.
Gate C: 8505+ composites x numpy-weight evaluation - minutes, bounded.
No RNG anywhere: gates A/B (and C) are seed-free exact computations; the
campaign stores integer/Fraction outputs only (floats only in tagged
NUMERICAL columns).

## Sources read first-hand (2026-08-29)

- npj QI 12, 86 (2026), DOI 10.1038/s41534-026-01227-w (= arXiv:2412.01850):
  full text + Methods + Table 2 landing page (published online 2026-04-04
  per nature.com metadata; arXiv v1 2024-11-28 per arXiv API).
- arXiv:2608.18935 abstract (2026-08-19) - follow-up campaign only.

## Signature

Gates and tolerances frozen 2026-08-29 before any campaign run. Revisions
must be appended below, never silently edited.

## Revision 1 (2026-08-29, after first partial gate-C run)

The Gate C pass/fail inequality as originally frozen was inverted (it
flagged maps with a LARGER shadow norm than U_ct, i.e. worse maps, as
"beats"). The intended scientific question - can any pinned family beat
the paper's 2 x 1.8^k - requires the opposite direction. Corrected above
to: PASS iff w_M > w_ct at some k <= 8 (smaller shadow norm). No other
gate is affected; gates A/B had already PASSED under the original text
before this revision, and their text is unchanged.
Reason: semantic check on the k <= 4 smoke run showed 39/39 family-(b)
composites flagged "beats" because a single disjoint pair rotation
scrambles, not contracts. Handedness of the comparison is a statement bug,
not a code bug; code and statement now agree on w_M > w_ct.


## Revision 2 (2026-08-30, appended before any k > 8 run)

Scope for the gate-C k-extension campaign (new folder under campaigns/ only;
Revisions 0/1 text above unchanged). Owner: shadows-ext.

### Questions and extension scope

Q1 (regression gate): recompute k = 4..8 with the extension code path and
assert exact Fraction equality against the frozen artifact
campaigns/2026-08-29T232134Z_4f1a9cab_e274dc65db35/gates_c.json for the
deterministic-layer maps (cz_mesh, cx_chain, cx_mesh_lex, uct_then_cz_mesh)
and the single-letter meshes (mesh_X, mesh_Y, mesh_Z). The extension's
k > 8 rows COUNT only if every regressed value matches; any mismatch aborts
the campaign before any k > 9 enumeration.

Q2 (extension): for each k in 9..12 compute w_M(k) and w_ct(k) with the
same exact machinery (GF(2) class maps, full 3^k token enumeration,
fractions.Fraction; no RNG, no floats in any verdict path) for the carrier
map M = cx_chain (gate-C deterministic-layer staircase of CNOTs) versus
U_ct = mesh_Z. Consistency assertion added now: for every completed k,
the enumerated w_ct(k) MUST equal Eq. (4) exactly (holds at k=11 in the
pre-run probe: 907505/1162261467 both ways), and mesh_Z == U_ct.

### Carrier-identity scope note (frozen before run, both facts named)

- Fact 1: at k = 4 the frozen artifact's improvement carrier is PURE
  cx_chain: w(cx_chain) = 361/6561 > w(mesh_Z = U_ct) = 353/6561, i.e.
  shadow-norm ratio ||cx_chain||^2 / ||U_ct||^2 = 353/361. The literal
  composite cx_chain after U_ct is a DIFFERENT map: at k = 4 its weight is
  1/81 < 353/6561 (norm LARGER), so it is NOT the improvement carrier and
  is excluded from the extension's per-k table. M in Q2 means cx_chain
  evaluated on the paper's normalized token ensemble.
- Fact 2: per the frozen artifact, cx_chain already loses to U_ct at
  k = 6, 7, 8 on weight (7937/531441 < 8177/531441; 12425/1594323 <
  13385/1594323; 175321/43046721 < 198593/43046721); the extension
  measures whether the k=4 edge returns at k = 9..12.

### Pinned pass/fail criterion (verdict per k in 9..12)

Direction pinned per Revision 1 semantics: EDGE-PERSISTS conveys a SMALLER
shadow norm, i.e. w_M(k) > w_ct(k) as exact Fractions (strictly). Verdicts:
- EDGE-PERSISTS iff w_M(k) > w_ct(k);
- EDGE-VANISHES iff w_M(k) <= w_ct(k);
- UNKNOWN-BUDGET iff k not reached before budget stop.
For every reported k, and also for the k <= 8 regression rows, the exact
per-k difference w_M(k) - w_ct(k) is reported so each verdict is
re-derivable without trust.

### Budget (pinned before run, single core, nice -n 10)

Total wall cap ~90 min. Measured pre-run (this machine): k=11 single-map
enumeration ~0.72 s/3^11 tokens (4.09 us/token); k=9..12 all four maps
together project to ~15 s total, so the budget is expected to be idle;
if any stage exceeds the cap, stop at the largest completed k and report
the remaining k as UNKNOWN-BUDGET (exact stopping statement required).

### Verdict propagation rule

The k <= 8 reproduction rows carry verdicts only as REGRESSION-OK/
REGRESSION-FAIL (they re-anchor the closed campaign; they are not new
edge claims). Only k = 9..12 rows carry EDGE-PERSISTS/EDGE-VANISHES/
UNKNOWN-BUDGET.
