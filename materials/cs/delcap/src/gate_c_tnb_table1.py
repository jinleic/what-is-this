"""
gate_c_tnb_table1.py — printing-convention test on Tavakoli-Nguyen-Bose
Table I (arXiv:2607.19559).

SCOPE (fixed by Main before running): transcribe the 18 printed entries, compare
each against the exact-Arb recomputation of THEIR OWN Cor.1/Thm.1 closed forms
already frozen in the Gate C campaign, report printed-versus-computed per row
WITH THE SIGN of the difference, and test whether a single-offset round-up fit
exists as it does for Morozov-Duman. This is a read plus a comparison: no new
certification, no new rows, no re-derivation. Certified values are READ from the
frozen campaign CSV.

WHY: two of the three primaries in this target print conservative round-UPS
rather than correctly-rounded values — Morozov-Duman (resolved this session: a
single offset <= 7.1e-7 plus 5-dp ceiling reproduces all 36 Table III rows) and
Pinto-Ribeiro (stated in their own text: "reported value = BA rate + tolerance",
an additive round-up). TNB is the third and unknown. So this tests whether
conservative round-up printing is a CONVENTION IN THIS LITERATURE.

Printed values transcribed first-hand from ar5iv TABLE I, "Capacity bounds
C_{q,n} for q in {2,3}", columns q | n | d | LB_1 (this paper) | LB_2 [12] |
KM [5] | LB+ (this paper) | C_{(q,n)} | UB [1], all to 3 decimals. Section V
states the ordering LB_1 <= LB_2 <= LB+ <= C_{q,n} <= UB.

NOTE on n-dependence: their LB_1 = (1-d)log2 q - h2(d) and UB = (1-d)log2 q are
n-INDEPENDENT, and their printed columns confirm this (identical down each
q,d group). So LB_1 and UB are compared on all 18 rows from values already
computed; LB+ carries Delta_n(d) and is compared on the 15 rows we computed.
"""
from __future__ import annotations
import os
import csv, hashlib, json, os, time, uuid
from decimal import Decimal, ROUND_CEILING, ROUND_HALF_UP, ROUND_FLOOR

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FROZEN = os.path.join(
    ROOT, 'campaigns',
    '2026-08-30T13:07:27Z_7988b619-e89e-4e77-8de7-2c2a37a429d9_23445e8dbeeb',
    'TABLE_sandwich_rows.csv')

# (q, n, d) -> (LB1, LB2, KM, LBplus, C_BA, UB); None where the paper prints "—"
PRINTED = {
    (2, 3, '1/20'): (0.664, 0.714, 0.730, 0.910, 0.910, 0.950),
    (2, 3, '1/10'): (0.431, 0.531, 0.569, 0.825, 0.827, 0.900),
    (2, 3, '1/5'): (0.078, 0.278, 0.372, 0.668, 0.676, 0.800),
    (2, 5, '1/20'): (0.664, 0.714, 0.730, 0.886, 0.887, 0.950),
    (2, 5, '1/10'): (0.431, 0.531, 0.569, 0.782, 0.786, 0.900),
    (2, 5, '1/5'): (0.078, 0.278, 0.372, 0.602, 0.613, 0.800),
    (2, 10, '1/20'): (0.664, 0.714, 0.730, 0.851, 0.852, 0.950),
    (2, 10, '1/10'): (0.431, 0.531, 0.569, 0.722, 0.728, 0.900),
    (2, 10, '1/5'): (0.078, 0.278, 0.372, 0.516, 0.531, 0.800),
    (3, 3, '1/20'): (1.219, 1.249, None, 1.453, 1.454, 1.506),
    (3, 3, '1/10'): (0.957, 1.016, None, 1.328, 1.329, 1.426),
    (3, 3, '1/5'): (0.546, 0.663, None, 1.095, 1.101, 1.268),
    (3, 5, '1/20'): (1.219, 1.249, None, 1.425, 1.426, 1.506),
    (3, 5, '1/10'): (0.957, 1.016, None, 1.276, 1.279, 1.426),
    (3, 5, '1/5'): (0.546, 0.663, None, 1.011, 1.020, 1.268),
    (3, 10, '1/20'): (1.219, 1.249, None, 1.386, 1.387, 1.506),
    (3, 10, '1/10'): (0.957, 1.016, None, 1.208, 1.212, 1.426),
    (3, 10, '1/5'): (0.546, 0.663, None, 0.908, 0.920, 1.268),
}
Q3 = Decimal('0.001')          # their printed precision: 3 decimals


def offset_band(pairs, mode):
    """Feasible additive offset Delta with printed = round_mode(computed+Delta)
    simultaneously over all (printed, computed) pairs. Returns (lo, hi, ok)."""
    lo, hi = Decimal('-1'), Decimal('1')
    for printed, computed in pairs:
        p = Decimal(str(printed))
        c = Decimal(repr(computed))
        if mode is ROUND_HALF_UP:
            a, b = p - Q3 / 2 - c, p + Q3 / 2 - c
        elif mode is ROUND_CEILING:
            a, b = p - Q3 - c, p - c
        else:                                  # truncation / floor
            a, b = p - c, p + Q3 - c
        lo = max(lo, a)
        hi = min(hi, b)
    return float(lo), float(hi), lo < hi


def main():
    stamp = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    code_hash = hashlib.sha256(open(__file__, 'rb').read()).hexdigest()[:12]
    camp = os.path.join(ROOT, 'campaigns', f'{stamp}_{uuid.uuid4()}_{code_hash}')
    os.makedirs(camp, exist_ok=True)

    cert = {(int(r['q']), int(r['n']), r['d']): r
            for r in csv.DictReader(open(FROZEN))}
    # LB1 and UB are n-independent: harvest one representative per (q,d)
    by_qd = {}
    for (q, n, d), r in cert.items():
        by_qd.setdefault((q, d), r)

    rows = []
    for key in sorted(PRINTED, key=lambda k: (k[0], k[1], k[2])):
        q, n, d = key
        lb1, lb2, km, lbp, cba, ub = PRINTED[key]
        rep = by_qd.get((q, d))
        c = cert.get(key)
        rec = dict(q=q, n=n, d=d, printed_LB1=lb1, printed_LBplus=lbp,
                   printed_C_BA=cba, printed_UB=ub)
        if rep is not None:
            my_lb1 = float(rep['their_LB1'])
            my_ub = float(rep['their_UB'])
            rec.update(my_LB1=my_lb1, my_UB=my_ub,
                       LB1_printed_minus_computed=lb1 - my_lb1,
                       UB_printed_minus_computed=ub - my_ub)
        if c is not None:
            my_lbp = float(c['their_LBplus'])
            rec.update(my_LBplus=my_lbp,
                       LBplus_printed_minus_computed=lbp - my_lbp,
                       cert_lo=float(c['cert_lo_per_symbol']),
                       cert_hi=float(c['cert_hi_per_symbol']),
                       cert_width=float(c['interval_width_per_symbol']))
            # separate observation: their BA column vs our certified enclosure
            rec['cert_lo_minus_printed_C_BA'] = rec['cert_lo'] - cba
            rec['printed_C_BA_below_certified_lower'] = rec['cert_lo'] > cba
        else:
            rec['LBplus_status'] = 'NOT_RECOMPUTED (q=3,n=10 dropped pre-registered)'
        rows.append(rec)

    # ---- convention test, per column
    report = {}
    for col, pkey, ckey in (('LB1', 'printed_LB1', 'my_LB1'),
                            ('LBplus', 'printed_LBplus', 'my_LBplus'),
                            ('UB', 'printed_UB', 'my_UB')):
        pairs = [(r[pkey], r[ckey]) for r in rows if ckey in r]
        diffs = [r[pkey] - r[ckey] for r in rows if ckey in r]
        # dedupe identical (printed, computed) pairs for the n-independent cols
        uniq = sorted(set(pairs))
        ent = dict(n_rows=len(pairs), n_distinct_pairs=len(uniq),
                   min_diff=min(diffs), max_diff=max(diffs),
                   n_printed_above=sum(1 for x in diffs if x > 0),
                   n_printed_below=sum(1 for x in diffs if x < 0),
                   n_exact=sum(1 for x in diffs if x == 0))
        for name, mode in (('round_half_up', ROUND_HALF_UP),
                           ('ceiling_roundup', ROUND_CEILING),
                           ('truncation', ROUND_FLOOR)):
            lo, hi, ok = offset_band(uniq, mode)
            ent[name] = dict(delta_lo=lo, delta_hi=hi, feasible=ok)
        report[col] = ent

    out = dict(stamp=stamp, code_sha256_12=code_hash,
               source='arXiv:2607.19559 TABLE I, transcribed first-hand',
               frozen_certified_table=FROZEN, printed_precision_dp=3,
               n_printed_rows=len(PRINTED), rows=rows, convention_test=report)
    with open(os.path.join(camp, 'tnb_table1_convention_test.json'), 'w') as f:
        json.dump(out, f, indent=1)
    with open(os.path.join(camp, 'tnb_table1_comparison.csv'), 'w') as f:
        f.write('q,n,d,printed_LB1,computed_LB1,LB1_printed_minus_computed,'
                'printed_LBplus,computed_LBplus,LBplus_printed_minus_computed,'
                'printed_UB,computed_UB,UB_printed_minus_computed,'
                'printed_C_BA,cert_lo,cert_hi,cert_width,'
                'cert_lo_minus_printed_C_BA\n')
        for r in rows:
            g = lambda k, fmt='%.9f': (fmt % r[k]) if k in r else ''
            f.write(f"{r['q']},{r['n']},{r['d']},{r['printed_LB1']},"
                    f"{g('my_LB1')},{g('LB1_printed_minus_computed','%+.6f')},"
                    f"{r['printed_LBplus']},{g('my_LBplus')},"
                    f"{g('LBplus_printed_minus_computed','%+.6f')},"
                    f"{r['printed_UB']},{g('my_UB')},"
                    f"{g('UB_printed_minus_computed','%+.6f')},"
                    f"{r['printed_C_BA']},{g('cert_lo')},{g('cert_hi')},"
                    f"{g('cert_width','%.3e')},"
                    f"{g('cert_lo_minus_printed_C_BA','%+.6f')}\n")

    print('=== PER-ROW: printed minus computed (their own formulas) ===')
    for r in rows:
        lbp = (f"{r['LBplus_printed_minus_computed']:+.6f}"
               if 'LBplus_printed_minus_computed' in r else '   (n/c)')
        print(f"  q{r['q']} n{r['n']:<2} d{r['d']:<5} "
              f"LB1 {r['printed_LB1']:.3f} vs {r.get('my_LB1', float('nan')):.6f} "
              f"({r.get('LB1_printed_minus_computed', float('nan')):+.6f})  "
              f"LB+ {r['printed_LBplus']:.3f} vs "
              f"{r.get('my_LBplus', float('nan')):.6f} ({lbp})  "
              f"UB {r['printed_UB']:.3f} vs {r.get('my_UB', float('nan')):.6f} "
              f"({r.get('UB_printed_minus_computed', float('nan')):+.6f})")
    print()
    print('=== CONVENTION TEST per column ===')
    for col, e in report.items():
        print(f"  {col}: rows={e['n_rows']} distinct={e['n_distinct_pairs']} "
              f"printed_above={e['n_printed_above']} "
              f"below={e['n_printed_below']} exact={e['n_exact']} "
              f"diff in [{e['min_diff']:+.6f}, {e['max_diff']:+.6f}]")
        for name in ('round_half_up', 'ceiling_roundup', 'truncation'):
            b = e[name]
            print(f"      {name:16s}: Delta in [{b['delta_lo']:+.3e}, "
                  f"{b['delta_hi']:+.3e}] feasible={b['feasible']}")
    print()
    below = [r for r in rows if r.get('printed_C_BA_below_certified_lower')]
    print(f'=== SEPARATE OBSERVATION (escalation, not part of the convention '
          f'test): printed C_(q,n) below our certified LOWER endpoint in '
          f'{len(below)} of {len([r for r in rows if "cert_lo" in r])} '
          f'compared rows ===')
    for r in below:
        print(f"  q{r['q']} n{r['n']:<2} d{r['d']:<5} printed C={r['printed_C_BA']:.3f} "
              f"certified [{r['cert_lo']:.9f}, {r['cert_hi']:.9f}] "
              f"cert_lo - printed = {r['cert_lo_minus_printed_C_BA']:+.6f}")
    print(camp)


if __name__ == '__main__':
    main()
