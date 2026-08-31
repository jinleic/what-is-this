# delcap Pre-Statement — Open Item 1: (3,10) and the q=3 ladder via sparse orbit certificates

Committed BEFORE any computation of this campaign. Owner: DelcapNextGate. Date: 2026-08-30T17:50Z.

## Which open item is attacked

Open item 1 of `cs/delcap/README.md:572-588` (What remains open), in the order
the contract lists them:

    * `(q,n) = (3,10)` and the `q=3` ladder above `n=5` need either a sparse
      channel representation or a GPU; the dense exact-integer matrix is the wall.

Explicitly NOT attacked in this campaign: the Morozov-Duman Table III residual
(the contract itself pre-commits "it will not be guessed at"; no instrument here
— their tau/p intermediates or code — can establish its cause, rule 14); the
Lambda-search global-optimality gap at m=22,23 (an exhaustive certification over
2^24 subsets needs a certified subset-sum ranking this workstation cannot do in
the budget; the current rows remain valid bounds either way, and that is stated,
not "settled"); Pinto-Ribeiro n in {29,31} (same wall, one letter larger: 2^29
inputs); and the n -> infinity capacity (rule 14: nothing in the finite-n
instrument can establish an asymptotic claim, and rule 7/rule 9 forbid calling
a finite-n verification asymptotic — the new rows below are finite-n theorems
at the named n).

## Falsifiable outcome being sought

Claim shape. The q-ary deletion channel W: F_q^n x Sigma^{<=n} (W(y|x) =
A(x,y) d^(n-k) (1-d)^k, A = #subsequences) is equivariant under the group
named in the GROUP-OBJECT LINE below.

GROUP-OBJECT LINE (requested by Main, 2026-08-30, and confirmed independently
in exact arithmetic by Main on (3,2,1/3), (4,2,1/5), (3,3,2/7) with the
position-permutation control failing as expected): the group is
G = S_q x C_2 acting on SYMBOL VALUES (alphabet relabeling sigma plus word
reversal rho), applied simultaneously to input and output words. Deletion
acts on POSITIONS; value relabeling/reversal commute with position deletion,
so A(sigma x, sigma y) = A(x,y) and A(rho x, rho y) = A(x,y) — the channel is
G-equivariant. This is NOT the falsified gate-B S_n input-POSITION/type
reduction (position permutation breaks equivariance). Claim shape: with W
G-equivariant (machine-verified as an anchor, never assumed), a
capacity-attaining input distribution exists on the G-invariant simplex, and
both certificate functionals factor through G-orbits:

  (a) PRIMAL. For p constant on input orbits, I(p) = sum_i P_i A_i - sum_j s_j u_j log2(u_j),
      with P_i = orbit_i mass = |O_i| p_i, A_i = sum over output-orbit classes
      (j, c) of m(c alpha_jk / D) log2(c alpha_jk / D) (m = multiplicity of the
      (j,c) class in the orbit incidence), s_j |O|_j^out = output orbit size,
      u_j = sum over incidence classes of c p_i m_ij / Dden. EVERY quantity is
      an exact rational; I is one outward-rounded Arb pass over O(25k) log
      terms at (3,10) instead of 5.2e9 dense cells. Valid LOWER bound: it is
      the exact mutual information of one explicit rational p, for ANY p.
  (b) DUAL. On the FULL q^n x sum q^k matrix: max_x KL(W(.|x) || D') attained
      over orbit representatives WHEN D' is G-invariant; for G-invariant D',
      KL(W(.|sigma x)||D') = KL(W(.|x)||D') is an exact identity checked
      numerically as an anchor. The certificate remains max over the FULL
      simplex for ANY full-support D' — representatives suffice only for
      invariant D' and the invariant candidate family is a SUBSET of candidates,
      so the dual bound stays a valid upper bound under exactly the same
      Csiszar-Tusnady hypothesis as the frozen pipeline (rule 15: the certified
      domain of the reduction is stated; no domain is moved after seeing data).

Rows targeted (exact, rule 16 — fixed before any run):
  * (q,n,d) = (3,10,d) for d in {1/2, 1/5, 1/10, 1/20}  - the named wall row,
    each with certified [primal, dual] bracket of C_{3,10}(d) per symbol, plus
    exact-Arb recomputation of their LB1, LB+ (new: via orbit-incidence Delta_10),
    UB from arXiv:2607.19559 Table I (printed 3-dp values for (3,10) already
    transcribed in the frozen TNB campaign).
  * (3,6), (3,7), (3,8), (3,9) at d = 1/2 and at the paper grid d in
    {1/20, 1/10, 1/5} where memory allows, completing the q=3 ladder above n=5.
  * One cross-validation row the dense code CAN do: (3,6) and (3,7) at d=1/2
    certified BOTH by the frozen dense path (gate_c_dhalf.row) and by the new
    orbit path; the two certified intervals must overlap (both contain the same
    C), and a disjoint pair falsifies the orbit implementation and stops the
    campaign.
  * The +infinity dual no-claim branch and the uniform-D' fallback carry over
    unchanged: a snapped D' entry of zero where W(y|x) > 0 returns +infinity
    (no claim), never a silent skip.

Adjudication rule (fixed before the first run):
  * ROW LANDS: orbit-path certified interval overlaps the dense-path interval
    at every cross-validation row AND is contained inside [their LB+, their UB]
    at the comparison rows; widths reported as first-class numbers AS MEASURED.
  * ROW MISSES (interval NOT inside the published sandwich, or wider than the
    sandwich): reported as measured, labelled, escalated to Main via hub
    BEFORE any README write; no adjustment toward any published value.
  * ANCHOR FAILS (equivariance identity or orbit-size identity fails at any
    probed point, or dense-vs-orbit intervals are disjoint): the orbit
    reduction is FALSIFIED for this pipeline; campaign stops; the failure is
    recorded with the exact probed point and the two intervals; nothing is
    built on the reduction.

## Budget

  * Wall-clock: at most ~2 hours of compute for (3,10) rows; each orbit-path
    row is expected in minutes (sparse incidence ~1e6 entries, not 5e9).
  * Cross-validation rows (3,6),(3,7) d=1/2 on the dense path: minutes each.
  * mpmath locator NOT used (reuses the established deviation: float64 BA
    locating step feeding the exact-rational snap — validity is independent of
    the locator; only width is affected).
  * Hard stop: if the (3,10) orbit incidence build exceeds 30 minutes, the
    (3,10) rows are reported as not reached with the timing evidence and the
    ladder rows n=6..9 carry the result.

## Precision and box

  * Arb working precision 400 bits, outward rounding, all exact rationals fmpq —
    the frozen standard. Snap denominators: input 2^30, output 2^30 with +1
    bump and uniform fallback, exactly as frozen. log2 via mpmath at 80 dps for
    the sandwich closed forms (their LB1/LB+/UB), Arb .log()/log2 inside
    certificate functionals — the frozen standard. Certified interval widths at
    400 bits reported as first-class numbers for every row.
  * Box (rule 16): q = 3 fixed; n in {10} for the headline row plus {6,7,8,9}
    as the ladder; d in {1/2, 1/5, 1/10, 1/20} at n=10 and the paper-grid
    subset d in {1/20, 1/10, 1/5} on the ladder rows as memory allows. Nothing
    is added, removed or re-graded after the first result; any additional row
    would be a separate pre-registered statement reported alongside a miss.
  * NOT swept (rule 7 sentence is maintained in the README): q >= 4 at any
    n > 5 (orbit machinery exists but rows are not pre-registered here); any
    n > 10; Morozov-Duman m > 23; the asymptotic capacity C(BDC_d); the
    published Pinto-Ribeiro n in {29,31} rows.

## Anti-numerology clause

A search that halts when it hits a target value is numerology. Every value is
reported including misses; the (3,10) row is run regardless of what the ladder
rows show; the ladder rows are run regardless of what (3,10) shows. The
cross-validation rows are run FIRST and gate everything: if the orbit path and
the dense path disagree at (3,6) or (3,7) d=1/2, the campaign stops there.

Signed: DelcapNextGate. No orbit-path computation run yet as of commit time.
