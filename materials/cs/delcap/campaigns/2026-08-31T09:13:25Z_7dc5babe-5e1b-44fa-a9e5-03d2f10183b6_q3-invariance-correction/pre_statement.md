# delcap correction pre-statement — exact orbit-total-mass replay of the frozen q=3 ladder

Committed as the FIRST file in this correction campaign, before the correction runner exists and before any replay computation. Owner: DelcapNextQ. Date: 2026-08-31T09:13:25Z.

## Trigger and non-negotiable separation of roles

The q=4 campaign's v3 startup regression found that the frozen accepted OUTPUT reference distribution `D` at `(q,n,d)=(3,6,1/5)` was not invariant: its per-WORD largest-remainder snap split flat output orbit 51 (size 12) between exact integer counts 578028 and 578029. The split concerns the 1,093-entry output reference distribution used by the DUAL. It does not by itself assert a defect in the input distribution `p` used by the PRIMAL. The frozen v3 failure artifacts remain immutable in the parent q=4 campaign. No q=4 row may run until this separate correction campaign succeeds.

The previous blanket claim that equal invariant float masses necessarily receive equal largest-remainder integer counts is retracted: stable tie-breaking may split an orbit. A count of mixed orbits is diagnostic only. It is not a substitute for the exact invariance, full-alphabet, positivity, and normalization checks below.

## Fixed correction box and order

Replay exactly the 20 prior q=3 rows in this order: `n=6,7,8,9,10`, and within each n, `d=(1/2,1/5,1/10,1/20)`. Nothing is added after seeing a result. No q=4 row, q=3 n=11 row, other d, asymptotic claim, or TNB rounding-cause claim is in this correction campaign.

For each n, build the exact `G=S_3 x C_2` structure once, where `S_3` acts on SYMBOL VALUES and `C_2` reverses words. Deletion acts on positions. The convention is left-to-right tuples and

`A(x,y)=#{i_1<...<i_k: x[i_1]...x[i_k]=y}`,
`W(y|x)=A(x,y)d^(n-k)(1-d)^k`.

## Primal: exact invariant input p, with a direct full-alphabet equality check

The float64 BA locator supplies only a candidate vector of TOTAL INPUT-ORBIT masses. Largest-remainder snap that vector at scale `2^30` to nonnegative integers summing exactly `2^30`. For input orbit `O_j`, assign every word `x in O_j` the exact rational probability

`p_x = m_j / (2^30 |O_j|)`.

Zero input masses are allowed. Assert exact sum one and exact invariance under a transposition and 3-cycle of symbol values plus reversal before any representative compression. Induce the output marginal `D_p` exactly from this p; assert its exact sum one and output-orbit invariance.

Compute the primal mutual information by TWO routes at Arb 400-bit outward rounding:

1. the orbit-weighted representative route; and
2. a direct full-input-alphabet route over all `3^n` words, with the exact per-word p and exact induced `D_p`.

The two Arb enclosures must overlap. Every full-word KL enclosure must overlap its input-orbit representative's enclosure. Any disagreement is a correction failure and stops the campaign.

## Dual: snap OUTPUT-orbit total masses, never per-word ties

The float64 locator marginal is used only to form prospective TOTAL OUTPUT-ORBIT masses `M_t^float = |T_t| D_t^float`. The fixed BA dual candidate `ba_total_orbit_mass` is constructed by largest-remainder snapping the vector `(M_t^float)_t` at scale `2^30` so the integer orbit masses sum exactly `2^30`. If any orbit receives zero, add one to EVERY output-orbit mass and use the resulting exact total `S`; this is the frozen full-support bump in orbit coordinates. Set, exactly,

`D_y = mass(T_t) / (S |T_t|)` for y in T_t.

The second and only fallback candidate is the predeclared `uniform` candidate: total orbit mass `|T_t|`, hence `D_y=1/sum_{k=0}^n 3^k` for every output word. Candidate order is fixed as `(ba_total_orbit_mass, uniform)`. Evaluate both exactly; choose the finite candidate with the smaller Arb UPPER endpoint, breaking an exact upper-endpoint tie in that fixed order. A `+infinity` candidate is a no-claim candidate, never a finite bound.

The legacy per-WORD `ba_word` integer vector is reconstructed for audit and archived. Any split rejects it before representative compression and before it can enter the corrected bound. It is never selected, even if a rep-only objective looks favorable. There is no objective-dependent filtering.

For EACH admitted dual candidate, before selection:

1. archive every output-orbit total-mass numerator, the exact total denominator S, every orbit size, and the exact per-word formula `mass/(S|T|)`;
2. assert every numerator is positive, exact total mass is one, and no positive W term sees `D_y=0`;
3. expand exact masses to every output word and assert invariance under reversal, the symbol transposition `(0 1)`, and the 3-cycle `(0 1 2)`;
4. evaluate Arb KL at every input-orbit representative and directly at every one of the `3^n` input words; every full-word KL enclosure must overlap its representative enclosure and the two maxima must overlap.

Only a candidate passing all four checks may enter the upper bound.

## Exact sandwich and replacement adjudication

Recompute the TNB finite-n sandwich at Arb 400 bits from the exact orbit Phi route:

`LB+ = (1-d)log2 3 + H_Bin(n,1-d)/n - h2(d) + Delta_n(d)/n`,
`UB = (1-d)log2 3`.

For each row archive: frozen old interval and source path; legacy ba_word split audit; exact input orbit masses; full-vs-orbit primal checks; both admitted dual candidates and their full-vs-representative checks; chosen candidate; corrected interval; width sign; and exact published margins.

Fixed verdicts use Arb endpoints only:

* `CERT_LOWER_BEATS_LBplus` iff corrected primal/n lower endpoint is greater than LB+ upper endpoint;
* `CERT_UPPER_BEATS_UB` iff corrected dual/n upper endpoint is less than UB lower endpoint.

A legacy split is a recorded retraction trigger, not a failure of the correction instrument. The campaign continues so all 20 old rows are classified. A failure of any NEW exact sum, positivity, invariance, full-alphabet equality/overlap, width-sign, or sandwich-containment check stops immediately and is escalated before README text. A corrected interval outside `[LB+,UB]` also stops and escalates. No arithmetic is tuned toward an old or published value.

## Anchors, counterfactuals, and resources

Before row 1:

* assert the frozen `(3,6,1/2)` old interval and primal ball bit-for-bit;
* re-run frozen A1-A6 and the raw `(2,10,1/20)` `+infinity` no-claim plant (exactly one zero under raw `2^-40` snap);
* re-run the out-of-grid `(3,3,1/2)` split plant and require the direct full-alphabet dual to strictly exceed the inadmissible representative-only value;
* run the new total-orbit-mass construction at `(3,6,1/2)` and require exact primal/full and dual/full checks to pass.

One low-priority process only: `nice -n 10`; OpenBLAS/OMP/MKL/vecLib/NumExpr threads pinned to 1. Float locator: NumPy/SciPy float64, 4,000 BA iterations, label `COMPUTATIONAL-EVIDENCE`, affects width only. Every probability and every certificate is exact `fmpq` plus Arb 400-bit outward rounding, label `MACHINE-VERIFIED`. Per-row hard stop: 90 minutes; correction-stage hard stop: 4 hours. Resource timing is `COMPUTATIONAL-EVIDENCE` only.

Rows are appended with per-line SHA-256 and fsync. Resume skips only checksum-valid complete rows. Scripts are frozen as `.asrun`; attempts, tool versions, manifests, and final checksums are archived.

## Mandatory rule-7 scope sentence

The correction sweep is exactly q=3 with `(n,d)` in `{6,7,8,9,10} x {1/2,1/5,1/10,1/20}`, using `G=S_3 x C_2` on symbol VALUES and reversal, exact total input-orbit mass snap `2^30`, exact total output-orbit mass snap `2^30`, full finite channel support for every dual candidate, direct full input alphabets of size `3^n`, and Arb 400-bit outward rounding. NOT swept: q=4 or q>=4; q=3 n>=11; other d; per-word tie-breaking as an admitted bound; non-G-invariant correction candidates; TNB rounding cause; Morozov-Duman, Pinto-Ribeiro, and every n->infinity/asymptotic capacity claim. Every replacement is only a finite-n theorem.

The split hazard was discovered retrospectively in the parent q=4 startup audit. The correction box, candidate order, mass formulas, checks, verdict inequalities, precision, and stops above are prospective and frozen before the correction runner exists.

Signed: DelcapNextQ.
