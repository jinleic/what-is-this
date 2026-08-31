"""
gate_c_freeze.py — assemble the immutable Gate C freeze campaign.

Collects every raw output produced by the Gate C runners, emits the per-row
deliverable tables (CSV), the pre-registration deviation log, the rule-5
retraction record, tool versions, checksums and a manifest. Deterministic:
re-running it produces the same tables from the same inputs.

Source campaigns consumed (each immutable, left untouched):
  A) 2026-08-30T12:27:34Z_a09c3817…  gate_c_dhalf.py, 19 rows
  B) 2026-08-30T12:37:07Z_26256c1e…  gate_c_dhalf.py, 11 rows
  C) 2026-08-30T12:34:14Z_4cd10de4…  gate_c_md.py: LO-CVB, extension, GAVB, PR
"""
from __future__ import annotations
import glob, hashlib, json, os, shutil, sys, time, uuid

ROOT = '/Users/jinleic/jinleic-workspace/cs/delcap'
CAMP = os.path.join(ROOT, 'campaigns')
SRC = os.path.join(ROOT, 'src')

SOURCES = {
    'A_dhalf_run1': '2026-08-30T12:27:34Z_a09c3817-e4b6-4984-8707-a97bafaedea6_b3495d48bbc6',
    'B_dhalf_run2': '2026-08-30T12:37:07Z_26256c1e-b090-4012-add9-5b77b0dfe488_9bee55a5e137',
    'C_md_pr': '2026-08-30T12:34:14Z_4cd10de4-efae-43f4-b9a6-39cab1db46b5_8ecf82513989',
    'D_pr_fastpath': '2026-08-30T13:04:38Z_b22da63b-6285-45a2-9fca-49aa90f8f011_c92652d2f9a9',
}

DEVIATIONS = """# Pre-registration deviations — Gate C

`pre_statement.md` was committed BEFORE any Gate C computation (its own header
records this, and a byte-identical copy is frozen in this campaign). It is NOT
edited retroactively. Every difference between the committed grid and what was
actually run is listed here, with cause. Rule 7: the swept set below is the
exact claim; nothing outside it is claimed.

## (A) Tavakoli-Nguyen-Bose sandwich rows

Committed: q in {2,3}, n in {3,5,10}, d in {1/20, 1/10, 1/5} (18 points), plus
extension to q=4, to q=3 n in {7,8}, and to d=1/2.

RAN, exactly:
  d = 1/2 (the ticket's target): (q,n) = (2,2) (2,3) (2,4) (2,5) (2,6) (2,7)
      (2,8) (3,2) (3,3) (3,4) (3,5) (4,3) (4,4)                = 13 rows
  paper grid: (q,n) = (2,3) (2,5) (2,10) (3,3) (3,5), each at
      d in {1/20, 1/10, 1/5}                                    = 15 rows
  28 distinct rows. Two rows — (3,3,1/2) and (3,5,1/2) — were additionally
  re-run in a SEPARATE process as a reproducibility check and agreed to every
  printed digit (0.594406224780 / 0.594406260821 and 0.512042 both times).

NOT RUN, with cause:
  * (q,n) = (3,10) at d in {1/20,1/10,1/5} — 3 committed points dropped.
    Cause: the exact channel is 3^10 = 59049 inputs by sum_k 3^k = 88573
    outputs, a dense exact-integer matrix of 5.2e9 entries — beyond this
    workstation. The pre-statement already flagged this point as at-risk.
  * extension (q,n) = (3,7), (3,8) at d=1/2 — not run. Cause: 3^7 x 3280 and
    3^8 x 9841 exact-rational certificate loops exceed the wall-clock budget.
    The q=3 ladder stops at n=5; the q=2 ladder at n=8 (n=10 on the paper grid).

METHOD DEVIATION (material, stated plainly):
  The pre-statement fixes the locating step as "mpmath 150-dps BA". These rows
  were produced with a float64 numpy BA locating step instead.
  Cause: at q=3, n=5 the mpmath-150dps locator cost ~290 s per row and, under
  heavy deletion (d=1/20), its iterates underflowed and collapsed onto a
  near-degenerate p, giving a valid but very loose dual.
  WHY THIS CANNOT AFFECT VALIDITY: the locator only proposes candidates. The
  certified lower bound is the exact mutual information I(p*) of an exact
  rational p* = m/M — a true MI, hence <= C for ANY p*. The certified upper
  bound is max_x KL(W(.|x) || D') for an exact rational full-support D' —
  hence >= C for ANY such D'. Both are one-shot outward-rounded Arb
  evaluations of logs of exact rationals; there is no interval iteration. A
  worse locator costs interval WIDTH and nothing else. Achieved widths are
  2.9e-10 to 6.7e-6 bits/symbol, tighter than the mpmath path achieved.

## (B) Pinto-Ribeiro C_{n,k} rows

Committed: (n,k) with k <= n <= 12, listed explicitly up to n=8.
RAN: n = 6,7,8,9 for all 1 <= k <= n — 30 rows.
NOT RUN: n in {10,11,12}. Cause: wall-clock.
Pinto-Ribeiro publish C_{n,k} only at n in {29,31} (GPU), out of certified
reach here — stated in the pre-statement as a structural limit, not discovered
late. Our n=6..9 rows therefore have NO published analogue: they are new
certified values, not a reproduction.

BECAUSE there is no published value to check them against, each row is
produced TWICE by different machinery and the two are compared:
  route 1 (gate_c_pr.py): float64 BA locator -> exact rational snap (2^30 in,
     2^40 out) -> gate_c_dhalf.cert_primal / cert_dual, Arb at 400 bits
  route 2 (delcap_cert, frozen): mpmath 150-dps BA locator -> exact rational
     snap (2^24) -> certify_ab, Arb at 300 bits
The routes share only the exact integer channel matrix subseq_matrix, which is
itself anchored by the frozen gate-B table. Both intervals must contain the
same true C_{n,k}, so the test is that they OVERLAP; a disjoint pair would mean
at least one certificate is wrong. RESULT: 22 of 30 rows have both routes
(route 2 was stopped at n=9,k=1 for cost — an mpmath row at n=9 costs
~10-80 min against ~1 s for route 1); all 22 OVERLAP, 0 disjoint. The
remaining 8 rows (n=9, k=2..9) are route-1 only and are labelled as such in
TABLE_pr_cnk_rows.csv (verdict FASTPATH_ONLY_no_mpmath_row).

## (C) Morozov-Duman LO-CVB rows

Committed: m in {5, 10, 15, 23}, E(m,w) recomputed from scratch for m <= 10,
anchored on their eq (28) values E(5,2)=32, E(5,3)=52, E(5,4)=54.
RAN: m in {5, 22, 23}.
DEVIATION: m in {10,15} dropped, m=22 added. Cause: their Table III prints
columns for m = 5, 22 and 23 only, so m=10 and m=15 have no published value to
compare against while m=22 does. Choosing the comparable columns is a
deviation from the committed set and is recorded as one.
NOT SWEPT (rule 7): the Lambda search is EXHAUSTIVE over all 2^6 = 64 subsets
at m=5. At m=22 and m=23 all 2^23 and 2^24 subsets are enumerated and ranked
in float64, and only the top 24 per row are certified in Arb. Since EVERY
Lambda yields a valid converse bound, incomplete certification can only make
the reported bound larger (weaker), never invalid. The reported value is the
min over CERTIFIED candidates — NOT a proven global optimum over Lambda.
"""

MD_FINDING = """# Morozov-Duman Table III — 36 rows recomputed and certified

Cleared for the record by Main after independent re-derivation of the offset
band. This is a REPRODUCTION with a named printing convention and a quantified
precision difference. It is not a discrepancy, a disagreement or a tension, and
no published claim is contradicted.

## What was computed

All 36 rows of Morozov-Duman (arXiv:2504.20961) Table III — the layer-oriented
converse bound (LO-CVB) code rate, their eqs (24)-(30) — at m in {5, 22, 23},
n in {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, infinity}, delta = 1/5 and
eps = 1/5 as EXACT rationals, recomputed from the integer E(m,w) tables.

Every certified value is one outward-rounded Arb evaluation of log2 of an exact
rational at 400 bits. Interval widths: 0.0 at every printed digit, maximum ball
RADIUS over all 36 rows 2.84e-120 (reported as a first-class number rather than
as "zero"). Truncation order: none — the LO-CVB is a finite closed form in the
E(m,w) integers, so there is no series to truncate and the tail bound is
exactly 0. Label: MACHINE-VERIFIED.

## Agreement with the published table

In all 36 rows the printed value EXCEEDS the certified value:

    gap = printed - certified  in  [7.10008e-07, 1.03994e-05]
    spread of gap over 36 rows  =  9.689365e-06  <  1.0e-05 = one 5dp ulp

Solving printed = round(certified + Delta) at 5 decimals SIMULTANEOUSLY over
all 36 rows, under ceiling printing (round up — the conservative and correct
convention for printing an upper bound, since a rounded-up upper bound remains
a valid upper bound):

    Delta in (3.990e-07, 7.10008e-07]

This band is DETERMINED, NOT FITTED: under ceiling printing gap lies in
[Delta, Delta + ulp), so the MINIMUM gap caps Delta from above and the MAXIMUM
gap bounds it from below at max_gap - ulp. Nothing was tuned. The band is
non-empty for exactly one reason — the spread is less than one ulp — and
"spread < 1 ulp" is equivalent to the existence of a single-Delta ceiling fit.
Under round-half-up printing the corresponding band is [5.399e-06, 5.710e-06].

So the published values agree with the certified ones to within 7.1e-7
absolute in every row, and the printed digits are exactly what round-up at five
decimals produces from a value that close.

## Why this is a precision difference and not a formula difference

A missing or differing FORMULA term would be n-dependent or m-dependent. Such
a term would BREAK the single-Delta fit, because the same Delta must reproduce
rows spanning n = 1 to n = infinity and m = 5, 22, 23 simultaneously. The fit
holds across all 36 rows, which is positive evidence AGAINST a formula
difference. A uniform positive offset bounded by 7.1e-7 at rate values
0.55-0.82 is 0.9e-6 to 1.3e-6 relative — the size of a 6-to-7-significant-digit
intermediate. Mechanism label: INFERENCE. The convention identification
(round-up printing) is INFERENCE. Our own values are MACHINE-VERIFIED.

## Aggregate check on the E(m,w) transcription

The n -> infinity row equals log2(sum_w E(m,w) delta^{m-w}(1-delta)^w)/m, a
weighted sum over the WHOLE E(m,·) column, and it reproduces their printed
value at every m: 0.80272 (m=5), 0.73569 (m=22), 0.73414 (m=23). This is an
aggregate check on the typed Table I columns, not a per-cell one.

Per-cell verification of E(m,w):
  * m=5, all w: recomputed from scratch by brute force over all (x,y) pairs.
    E(5,·) = [1, 10, 32, 52, 54, 32]; their eq (28) prints E(5,2)=32,
    E(5,3)=52, E(5,4)=54 — reproduced exactly. MACHINE-VERIFIED.
  * m=20,21,22,23, columns w in {0,1,2,m}: re-derived independently and agree
    with the typed values. w=0: E=1. w=1: E=2m. w=m: E=2^m (their eq 28).
    w=2: closed form E(m,2) = 2*C(m,2) + 2*floor(m^2/4), derived and proved in
    gate_c_md.py and brute-force validated at m=4..12; it reproduces their
    typed 580, 640, 704, 770 exactly. MACHINE-VERIFIED for those cells.
  * all other cells of m=20..23: CITED-DEPENDENCY (typed from their Table I).

## Status

First certified enclosures of these 36 rows. No published claim is
contradicted; their converse stands. If anything the round-up printing
convention is the correct choice for a converse bound.
"""

RETRACTIONS = """# Retractions — Gate C (rule 5: recorded inline, never removed)

## R1 — "the LO-CVB rows disagree with the paper" — RETRACTED (cleared)

ORIGINAL CLAIM (mine, 2026-08-30, escalated to Main): 19 of 36 LO-CVB rows
"disagree" with Morozov-Duman's printed Table III, my certified value lying
below theirs by 5.0e-6 to 1.04e-5, and the cause is "last-digit rounding in a
printed table".

FALSIFIED IN TWO STAGES.

(a) The stated CAUSE was falsified by Main with exact decimal arithmetic. For
m=22,n=32 my certified 0.66389960 rounds to 0.66390 under BOTH round-half-up
and ceiling while the printed value is 0.66391; for m=5,n=1024 my certified
0.80277963 rounds to 0.80278 under both, printed 0.80279. The differences
1.040e-5 and 1.037e-5 exceed the 5dp ulp of 1.0e-5, and a printed-table
rounding artifact of my own UNSHIFTED value cannot exceed one ulp. Withdrawn.

(b) The 19-vs-17 PARTITION was then falsified by me. Sorting all 36 rows by
gap = printed - certified, my partition falls exactly at the sort boundary:
every row with gap < 5.04e-6 was labelled "reproduces", every row with
gap >= 5.04e-6 "strictly tighter". The split was a cut at the +-5e-6
comparison window I chose myself, and carries ZERO structural information. The
two classes are interleaved on every structural property: m in {5,22,23} on
both sides, full-support Lambda 11/19 vs 6/17, n from 1 to infinity on both
sides, |Lambda| from 5 to 24 on both sides. Correlations of gap are weak and
consistent with noise at n=36: n +0.40, m*n +0.28, is_full_set +0.20, m -0.08,
|Lambda| -0.03.

(c) A third hypothesis, Main's, was also falsified: that the gaps show a
SYSTEMATIC ADDITIVE finite-n term, inferred from two quoted gaps sitting
3.0e-8 apart. The gaps are a CONTINUUM — 35 distinct values among 36 rows at
1e-8 resolution — and those two are simply the top two of it (third:
1.0238e-5), a top-of-distribution proximity artifact.

(d) CORRECTION I VOLUNTEERED AGAINST MY OWN REPORT: the direction split is
36-0, not 19-0. In every one of the 36 rows the printed value exceeds the
certified value; minimum gap +7.1e-7.

WHAT REPLACED IT: a reproduction with a named printing convention. See
MD_TABLE_III_FINDING.md in this campaign. Delta band under ceiling printing
(3.990e-7, 7.10008e-7], determined by the min and max gaps rather than fitted,
non-empty precisely because the spread 9.689365e-6 is less than one 5dp ulp.
No published claim is contradicted.

## R2 — no other Gate C claim has been retracted.

The d=1/2 sandwich improvement, the paper-grid recomputation rows, the
Pinto-Ribeiro extension rows and the delta in {1/20,1/2,4/5} extension rows
stand as first reported.
"""

SOLVER_AUDIT = """# Solver / tolerance audit — Gate C

Prompted by a repo-wide warning from Main: a scipy HiGHS call elsewhere in the
repo returned a point violating its own bounds by 100% because the solver's
default primal feasibility tolerance (1e-7) EQUALLED the box radius being
enforced — an instrument that cannot resolve the quantity it is asked about.

## Finding: no solver anywhere in the Gate C certified path

grep for scipy / linprog / highs / optimize across gate_c_dhalf.py,
gate_c_md.py, gate_c_run.py and delcap_cert.py returns 0 hits. There is no
LP/QP/MILP in the certified path. Float appears in exactly two places, both
non-load-bearing: (a) a float64 numpy BA locating step, (b) float64 ranking of
Lambda subsets in the LO-CVB search. Both feed an exact-rational snap, after
which every certified number is an outward-rounded Arb evaluation of logs of
exact rationals at 400 bits.

Tolerance versus magnitude: certified per-symbol widths run 2.9e-10 to 6.7e-6.
No float tolerance is ever compared against anything at that scale. The
sandwich comparisons are Arb-endpoint inequalities (certified lower endpoint vs
the exact-Arb bound's upper endpoint) with margins of 0.02-0.13 bits/symbol —
four to six orders of magnitude above any interval width.

## A live invalid-certificate branch, found before it bit

The dual certificate is max_x KL(W(.|x) || D') for an exact rational D'. If
D'(y) = 0 while W(y|x) > 0 for some x, that KL is +infinity. Silently skipping
such a term would yield a too-SMALL "upper bound" — an invalid certificate that
looks perfectly clean.

This is not hypothetical. VERIFIED at q=2, n=10, d=1/20:
  * output alphabet 2047; BA output marginal minimum mass 9.765625e-14
  * snap resolution 2^-40 = 9.094947e-13, i.e. 9.3x COARSER than the mass it
    must represent
  * exactly 1 of 2047 entries snaps to zero, and that output IS reachable
    (W(y|x) > 0 for some x)
  * cert_dual on the raw snap returns +inf (no claim) — NOT a finite too-small
    number. The unsafe branch does not exist in the code.
  * with the +1 support bump: 8.529977690111 (valid, tight)
  * with uniform D': 9.672330910401 (valid, looser)
  * the row takes the min over candidates, so the tight valid one is used.

Same failure shape as the HiGHS trap — a resolution coarser than the quantity
being resolved — caught pre-emptively rather than after a retraction.
"""


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for blk in iter(lambda: f.read(1 << 20), b''):
            h.update(blk)
    return h.hexdigest()


def main():
    stamp = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    code_hash = hashlib.sha256(open(__file__, 'rb').read()).hexdigest()[:12]
    out_root = os.environ.get('GATEC_FREEZE_OUT', CAMP)   # dry-run override
    out = os.path.join(out_root, f'{stamp}_{uuid.uuid4()}_{code_hash}')
    os.makedirs(out, exist_ok=True)
    prov = {}

    for fn in ('gate_c_dhalf.py', 'gate_c_md.py', 'gate_c_pr.py',
               'gate_c_freeze.py', 'delcap_cert.py'):
        shutil.copy(os.path.join(SRC, fn), os.path.join(out, fn))
        prov[fn] = f'src/{fn} as run'
    shutil.copy(os.path.join(ROOT, 'pre_statement.md'),
                os.path.join(out, 'pre_statement.md'))
    prov['pre_statement.md'] = 'committed before any Gate C run; byte-identical'

    rows_all = []
    for tag, d in SOURCES.items():
        sd = os.path.join(CAMP, d)
        if not os.path.isdir(sd):
            print(f'WARNING missing source campaign {d}')
            continue
        for f in sorted(os.listdir(sd)):
            src = os.path.join(sd, f)
            if os.path.isfile(src) and f != 'checksums.sha256':
                shutil.copy(src, os.path.join(out, f'{tag}__{f}'))
                prov[f'{tag}__{f}'] = f'raw output, campaign {d}'
        rp = os.path.join(sd, 'rows.jsonl')
        if os.path.isfile(rp):
            rows_all += [json.loads(L) for L in open(rp)]

    seen, uniq = set(), []
    for r in rows_all:
        key = (r['q'], r['n'], r['d'])
        if key not in seen:
            seen.add(key)
            uniq.append(r)
    uniq.sort(key=lambda r: (r['q'], r['n'], r['d']))
    with open(os.path.join(out, 'TABLE_sandwich_rows.csv'), 'w') as f:
        f.write('q,n,d,their_LB1,their_LBplus,their_UB,cert_lo_per_symbol,'
                'cert_hi_per_symbol,interval_width_per_symbol,primal_ball_radius,'
                'dual_ball_radius,truncation_order,tail_bound,arb_prec_bits,'
                'snap_denom_in,snap_denom_out,nX,nY,verdict\n')
        for r in uniq:
            f.write(f"{r['q']},{r['n']},{r['d']},{r['lb1']:.10f},"
                    f"{r['lbplus']:.10f},{r['ub']:.10f},"
                    f"{r['cert_lo_per_symbol']:.12f},{r['cert_hi_per_symbol']:.12f},"
                    f"{r['cert_width_per_symbol']:.4e},{r['primal_rad']:.3e},"
                    f"{r['dual_rad']:.3e},none-exact-full-channel,0,400,"
                    f"{r['snap_in']},{r['snap_out']},{r['nX']},{r['nY']},"
                    f"{r['verdict']}\n")
    prov['TABLE_sandwich_rows.csv'] = f'deliverable table, {len(uniq)} distinct rows'

    md = os.path.join(CAMP, SOURCES['C_md_pr'])
    lp = os.path.join(md, 'locvb_rows.jsonl')
    if os.path.isfile(lp):
        lr = [json.loads(L) for L in open(lp)]
        with open(os.path.join(out, 'TABLE_locvb_rows.csv'), 'w') as f:
            f.write('m,n,delta,eps,their_printed,cert_lo,cert_hi,interval_width,'
                    'ball_radius,gap_printed_minus_cert,lambda_size,'
                    'lambda_is_full_set,lambda_search,truncation_order,tail_bound,'
                    'E_provenance\n')
            for r in sorted(lr, key=lambda r: (r['m'], 10 ** 9 if r['n'] == 'inf'
                                               else r['n'])):
                f.write(f"{r['m']},{r['n']},{r['delta']},{r['eps']},"
                        f"{r['printed']:.5f},{r['cert_lo']:.12f},{r['cert_hi']:.12f},"
                        f"{r['width']:.1e},{r['rad']:.3e},"
                        f"{r['printed'] - r['cert_hi']:.5e},{len(r['Lambda'])},"
                        f"{int(r['is_full_set'])},{r['lambda_search']},"
                        f"none-finite-closed-form,0,"
                        f"\"{r['E_provenance']}\"\n")
        prov['TABLE_locvb_rows.csv'] = 'deliverable table, 36 LO-CVB rows'

    ep = os.path.join(md, 'locvb_extension_rows.jsonl')
    if os.path.isfile(ep):
        er = [json.loads(L) for L in open(ep)]
        with open(os.path.join(out, 'TABLE_locvb_extension_rows.csv'), 'w') as f:
            f.write('m,n,delta,eps,cert_lo,cert_hi,interval_width,lambda_size,'
                    'lambda_is_full_set,published_analogue\n')
            for r in er:
                f.write(f"{r['m']},{r['n']},{r['delta']},{r['eps']},"
                        f"{r['cert_lo']:.12f},{r['cert_hi']:.12f},{r['width']:.1e},"
                        f"{len(r['Lambda'])},{int(r['is_full_set'])},"
                        f"none-MD-plot-only-Fig2\n")
        prov['TABLE_locvb_extension_rows.csv'] = f'{len(er)} new certified rows'

    prp = os.path.join(CAMP, SOURCES['D_pr_fastpath'], 'pr_rows_fastpath.jsonl')
    if os.path.isfile(prp):
        with open(os.path.join(out, 'TABLE_pr_cnk_rows.csv'), 'w') as f:
            f.write('n,k,cert_lo,cert_hi,interval_width,ball_radius_primal,'
                    'ball_radius_dual,route2_mpmath_lo,route2_mpmath_hi,'
                    'routes_overlap,published_analogue,truncation_order,'
                    'tail_bound,verdict\n')
            for L in open(prp):
                r = json.loads(L)
                r2lo = r.get('mpmath_route_lo')
                r2hi = r.get('mpmath_route_hi')
                f.write(f"{r['n']},{r['k']},{r['cert_lo']:.12f},"
                        f"{r['cert_hi']:.12f},{r['width']:.4e},"
                        f"{r['primal_rad']:.3e},{r['dual_rad']:.3e},"
                        f"{'' if r2lo is None else f'{r2lo:.10f}'},"
                        f"{'' if r2hi is None else f'{r2hi:.10f}'},"
                        f"{'' if 'routes_overlap' not in r else int(r['routes_overlap'])},"
                        f"none-PR-publish-only-n29-n31,none-exact-full-channel,0,"
                        f"{r['verdict']}\n")
        prov['TABLE_pr_cnk_rows.csv'] = ('deliverable table, PR extension rows, '
                                        'two routes where available')

    # ---- extra route-2 rows from a separate parallel mpmath run (n=9), kept as
    # raw evidence with their overlap verdict against the route-1 table
    extra = sorted(glob.glob('/tmp/pr9/r9_*.json'))
    if extra:
        fast = {}
        for L in open(prp):
            r = json.loads(L)
            fast[(r['n'], r['k'])] = (r['cert_lo'], r['cert_hi'])
        with open(os.path.join(out, 'pr_extra_mpmath_crosscheck.txt'), 'w') as f:
            f.write('Additional route-2 (mpmath 150-dps locator, delcap_cert '
                    'certify_capacity) rows at n=9 from a separate parallel\n'
                    'run, cross-checked against the route-1 fast-path table. '
                    'Both intervals must contain the same true C_{n,k},\n'
                    'so the test is overlap.\n\n')
            for p in extra:
                r = json.loads(open(p).read())
                key = (r['n'], r['k'])
                shutil.copy(p, os.path.join(out, f'E_pr9_mpmath__{os.path.basename(p)}'))
                prov[f'E_pr9_mpmath__{os.path.basename(p)}'] = 'route-2 raw row'
                if key in fast:
                    L1, H1 = fast[key]
                    ov = not (H1 < r['lo_f'] or L1 > r['hi_f'])
                    f.write(f'C_{r["n"]},{r["k"]}: route1 [{L1:.12f}, {H1:.12f}]'
                            f'  route2 [{r["lo_f"]:.12f}, {r["hi_f"]:.12f}]'
                            f'  overlap={ov}  route2_wall={r["wall_s"]}s\n')
        prov['pr_extra_mpmath_crosscheck.txt'] = 'extra two-route check at n=9'

    for fn, body in (('DEVIATIONS.md', DEVIATIONS),
                     ('RETRACTIONS.md', RETRACTIONS),
                     ('MD_TABLE_III_FINDING.md', MD_FINDING),
                     ('SOLVER_AUDIT.md', SOLVER_AUDIT)):
        with open(os.path.join(out, fn), 'w') as f:
            f.write(body)
        prov[fn] = 'record'

    # ---- ANCHOR: re-certify two already-frozen gate-B rows with the GATE-C
    # certificate code path (the one that produced every Gate C number), and
    # compare against the frozen intervals. This isolates the certificate
    # layer: the channel-builder layer is anchored separately by the
    # Tavakoli Examples 1-3 self-test inside gate_c_dhalf.py.
    sys.path.insert(0, SRC)
    import importlib.util
    import numpy as np
    sp0 = importlib.util.spec_from_file_location('gh0', os.path.join(SRC, 'gate_c_dhalf.py'))
    gh0 = importlib.util.module_from_spec(sp0)
    sp0.loader.exec_module(gh0)
    from delcap_cert import subseq_matrix
    FROZEN = {(3, 2): (1.4697819938, 1.4697820261, 1.469781993757231),
              (5, 3): (1.8715442105, 1.8715442453, None)}
    with open(os.path.join(out, 'anchor_check.txt'), 'w') as f:
        f.write('ANCHOR: frozen gate-B rows re-certified with the Gate C '
                'certificate code path (gate_c_dhalf.cert_primal / cert_dual)\n'
                'frozen source: campaigns/2026-08-30T01:43:09Z_eabfc720f2ef '
                'and campaigns/enclosure_table.csv (printed to 10 dp)\n\n')
        for (n, k), (flo, fhi, rc28_dual) in FROZEN.items():
            Ni = np.array(subseq_matrix(n, k)[0], dtype=object)
            Cnk = subseq_matrix(n, k)[1]
            Nf = np.array([[float(Ni[i, j]) / Cnk for j in range(Ni.shape[1])]
                           for i in range(Ni.shape[0])])
            pb, Db, _, _ = gh0.ba_float(Nf, iters=4000)
            mi = gh0.snap(pb, 1 << 30)
            mo = gh0.snap(Db, 1 << 40)
            if min(mo) == 0:
                mo = [v + 1 for v in mo]
            lo = gh0.cert_primal(Ni, Cnk, mi)
            hi = min([gh0.cert_dual(Ni, Cnk, mo),
                      gh0.cert_dual(Ni, Cnk, [1] * Ni.shape[1])],
                     key=lambda a: a.upper())
            L, H = float(lo.lower()), float(hi.upper())
            f.write(f'f({n},{k}):\n'
                    f'  frozen printed interval [{flo:.10f}, {fhi:.10f}]\n'
                    f'  gate-C code path        [{L:.13f}, {H:.13f}]\n'
                    f'  width {H - L:.4e}   contained in frozen interval: '
                    f'{L >= flo - 5e-11 and H <= fhi + 5e-11}\n'
                    f'  primal ball {str(lo)[:64]}\n'
                    f'  dual   ball {str(hi)[:64]}\n')
            if rc28_dual is not None:
                f.write(f'  independent frozen rc28 direct dual {rc28_dual:.15f}'
                        f'  |ours - rc28| = {abs(H - rc28_dual):.3e}\n')
            f.write('\n')
    prov['anchor_check.txt'] = 'anchor: frozen gate-B rows via gate-C cert code'

    # ---- SECOND ROUTE for the LB+ column: their closed form (Cor. 1) vs the
    # certified exact MI at uniform input. No shared machinery, so agreement
    # overdetermines the sandwich implementation, the channel builder and
    # cert_primal simultaneously.
    from flint import fmpq as _fmpq
    with open(os.path.join(out, 'second_route_check.txt'), 'w') as f:
        f.write('SECOND ROUTE, LB+ column: Tavakoli Cor. 1 closed form '
                '(h2, H_Bin, Delta_n from exact integer pattern counts)\n'
                'vs certified exact mutual information at UNIFORM input '
                '(full channel matrix, cert_primal). Their Cor. 1 asserts\n'
                'these coincide; a defect in either route breaks agreement.\n\n')
        for (q, n, dd) in [(2, 3, _fmpq(1, 2)), (2, 5, _fmpq(1, 2)),
                           (2, 3, _fmpq(1, 20)), (2, 5, _fmpq(1, 5)),
                           (3, 3, _fmpq(1, 2)), (3, 3, _fmpq(1, 10)),
                           (4, 3, _fmpq(1, 2))]:
            sw = gh0.sandwich(q, n, dd)
            Nc, Dd, _o = gh0.full_bdc_channel(q, n, dd)
            pu = gh0.cert_primal(Nc, Dd, [1] * Nc.shape[0])
            lbp = float(sw['lbplus'])
            lo2, hi2 = float((pu / n).lower()), float((pu / n).upper())
            f.write(f'q={q} n={n} d={dd}: route1 LB+ = {lbp:.15f}\n'
                    f'    route2 MI(uniform)/n = [{lo2:.15f}, {hi2:.15f}]  '
                    f'LB+ inside: {lo2 <= lbp <= hi2}  '
                    f'|diff| = {abs(lbp - lo2):.3e}\n')
    prov['second_route_check.txt'] = 'overdetermination check for the LB+ column'

    # ---- regenerate the support-guard evidence live, so the campaign carries
    # its own raw proof of the audit claim in SOLVER_AUDIT.md
    sys.path.insert(0, SRC)
    import importlib.util
    import numpy as np
    from flint import fmpq
    sp = importlib.util.spec_from_file_location('gh', os.path.join(SRC, 'gate_c_dhalf.py'))
    gh = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(gh)
    N, Dden, outs = gh.full_bdc_channel(2, 10, fmpq(1, 20))
    Nf = np.array([[float(N[i, j]) / Dden for j in range(N.shape[1])]
                   for i in range(N.shape[0])])
    p_ba, D_ba, _, _ = gh.ba_float(Nf, iters=1500)
    m_out = gh.snap(D_ba, 1 << 40)
    zeros = [j for j, v in enumerate(m_out) if v == 0]
    reachable = [j for j in zeros if Nf[:, j].max() > 0]
    raw = gh.cert_dual(N, Dden, m_out)
    bumped = gh.cert_dual(N, Dden, [v + 1 for v in m_out]) if zeros else raw
    unif = gh.cert_dual(N, Dden, [1] * N.shape[1])
    with open(os.path.join(out, 'support_guard_check.txt'), 'w') as f:
        f.write('D\'(y)=0 support-guard verification, q=2 n=10 d=1/20\n'
                f'output alphabet nY = {N.shape[1]}\n'
                f'min BA output mass = {float(D_ba.min()):.6e}\n'
                f'snap resolution 2^-40 = {2.0 ** -40:.6e}\n'
                f'ratio (resolution / mass) = {2.0 ** -40 / float(D_ba.min()):.3f}x coarser\n'
                f'entries snapped to zero = {len(zeros)}\n'
                f'of those REACHABLE (W(y|x)>0 for some x) = {len(reachable)}\n'
                f'cert_dual on raw snap  = {float(raw.upper())}  '
                f'(+inf means NO CLAIM, the safe branch)\n'
                f'cert_dual on +1 bumped = {float(bumped.upper())}\n'
                f'cert_dual on uniform   = {float(unif.upper())}\n'
                'the row takes min over candidates, so a valid tight dual is used\n')
    prov['support_guard_check.txt'] = 'live regeneration of the audit evidence'

    import flint, numpy, mpmath, platform
    with open(os.path.join(out, 'tool_versions.txt'), 'w') as f:
        f.write(f'python {sys.version}\nplatform {platform.platform()}\n'
                f'python-flint {flint.__version__}\nnumpy {numpy.__version__}\n'
                f'mpmath {mpmath.__version__}\narb working precision 400 bits\n')
    prov['tool_versions.txt'] = 'environment as run'

    man = dict(campaign=os.path.basename(out), stamp=stamp,
               freeze_script_sha256_12=code_hash, sources=SOURCES,
               n_sandwich_rows=len(uniq), owner='DelcapGateC',
               note='Gate C freeze. The 36 Morozov-Duman LO-CVB rows were '
                    'escalated to Main and CLEARED after independent '
                    're-derivation of the offset band; they are recorded as a '
                    'reproduction with a named round-up printing convention, '
                    'not as a discrepancy. See MD_TABLE_III_FINDING.md and '
                    'RETRACTIONS.md.',
               provenance=prov)
    with open(os.path.join(out, 'manifest.json'), 'w') as f:
        json.dump(man, f, indent=1)
    with open(os.path.join(out, 'checksums.sha256'), 'w') as f:
        for fn in sorted(os.listdir(out)):
            p = os.path.join(out, fn)
            if os.path.isfile(p) and fn != 'checksums.sha256':
                f.write(f'{sha256(p)}  {fn}\n')
    print(out)
    print(f'{len(uniq)} distinct sandwich rows; {len(prov)} provenance entries')


if __name__ == '__main__':
    main()
