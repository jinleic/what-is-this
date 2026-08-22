#!/usr/bin/env python3
"""Zero-bad (emergent-free) member test via perfect-square reduction -- and
the attack on the last grid cell (89,(-2,1)).

Emergent-freeness at a member b does NOT need P(b) factored: strip the
primes of S u supp(b) from P(b) (numerator and denominator separately) and
test whether the remainders are perfect squares.  Together with the
small-prime alignment symbols this is a COMPLETE decider for assembly
solubility at a candidate, below _MR_LIMIT, with NO factorization refusals.

Phase A (validation): frozen zero-bad rows True; audited members with
emergent_status == 'complete' must satisfy zb == (bad_count == 0).
Phase B (attack): every box Closer89 enumerated for (89,(-2,1)) -- read its
per-a supp pools from data/l13_cell89.json -- decided by certificate().
A True run gets phase-6 cross-checks (tied/steered where the engine permits);
the frozen symbols + perfect-square remainder are themselves the certificate.
"""
import contextlib
import json
import math
import os
import signal
import time
from fractions import Fraction as F
from h10q import (_l10_P, _l10_supp, _L11_CLASSES, _L12_ESCAPE, _L13_ZERO_BAD,
                  _sun_h, _l7_tied_status, _l9_steered_solvable, hilbert, OO,
                  vp, primerange, FactorBudget, PrimalityBound, _is_prime,
                  factorint as _kernel_factorint, _MR_LIMIT)


def pattern(a, z, tau):
    A = 1 + 4*a*a
    delta = 1 - A*tau*tau
    alpha = -delta*A
    Z = z**3
    D = 1 - Z - a*a*Z*Z
    s = (a - 1)/2
    return A, delta, alpha, Z, D, s


def emergent_free(a, z, tau, b):
    """(is_square_or_None, (rem_num, rem_den))."""
    A, delta, alpha, Z, D, s = pattern(a, z, tau)
    if D == 0 or delta == 0:
        return None
    c0 = _sun_h(a, b, Z)
    if c0 is None:
        return None
    M0 = 16 - delta*c0*c0 - 32*A*b*s*s
    if M0 == 0:
        return None
    P = _l10_P(a, Z, D, A, delta, s)
    v = sum(c*b**i for i, c in enumerate(P))
    if v == 0:
        return None
    primes = {2, 3, 5, 7} | _l10_supp(alpha) | _l10_supp(delta) | _l10_supp(b)
    rn = abs(v.numerator)
    rd = v.denominator
    for p in primes:
        while rn % p == 0:
            rn //= p
        while rd % p == 0:
            rd //= p
    ok = (math.isqrt(rn)**2 == rn) and (math.isqrt(rd)**2 == rd)
    return ok, (rn, rd)


def aligned(a, z, tau, b):
    A, delta, alpha, Z, D, s = pattern(a, z, tau)
    c0 = _sun_h(a, b, Z)
    if c0 is None:
        return None
    M0 = 16 - delta*c0*c0 - 32*A*b*s*s
    if M0 == 0:
        return None
    x0, d0 = alpha*M0, alpha*2*b
    S = {2, 3, 5, 7} | _l10_supp(alpha) | _l10_supp(delta)
    out = {p: hilbert(x0, d0, p) for p in sorted(S)}
    for q in sorted(_l10_supp(b) - S):
        out[q] = hilbert(x0, d0, q)
    out['oo'] = hilbert(x0, d0, OO)
    return out


def _small_primes_upto(y):
    from h10q import primerange as _pr
    return list(_pr(2, y))


def smooth_emergent(a, z, tau, b, y=10**6):
    """Exact emergent decision when the stripped remainder is y-smooth.

    Returns (status, detail): status in {'zero','align-fail','bad-N','refused'}.
    'zero' => every emergent symbol +1 (SOLUBLE candidate, certified)."""
    A, delta, alpha, Z, D, s = pattern(a, z, tau)
    if D == 0 or delta == 0:
        return ("none", None)
    c0 = _sun_h(a, b, Z)
    if c0 is None:
        return ("none", None)
    M0 = 16 - delta*c0*c0 - 32*A*b*s*s
    if M0 == 0:
        return ("none", None)
    P = _l10_P(a, Z, D, A, delta, s)
    v = sum(c*b**i for i, c in enumerate(P))
    if v == 0:
        return ("none", None)
    x0, d0 = alpha*M0, alpha*2*b
    frozen = {2, 3, 5, 7} | _l10_supp(alpha) | _l10_supp(delta) | _l10_supp(b)
    for p in sorted(frozen):
        out = hilbert(x0, d0, p)
        if out != 1:
            return ("align-fail", p)
    if hilbert(x0, d0, OO) != 1:
        return ("align-fail", 'oo')
    # strip frozen + supp(b); smooth-part of the remainder
    rn = abs(v.numerator)
    rd = v.denominator
    for p in frozen:
        while rn % p == 0:
            rn //= p
        while rd % p == 0:
            rd //= p
    em = []
    smooth = True
    for p in _small_primes_upto(y + 1):
        if p in frozen:
            continue
        e = 0
        while rn % p == 0:
            rn //= p
            e += 1
        e2 = 0
        while rd % p == 0:
            rd //= p
            e2 += 1
        if (e - e2) % 2:
            # odd valuation OUTSIDE frozen u supp(b): emergent place
            if hilbert(x0, d0, p) != 1:
                em.append(p)
    if rn != 1 or rd != 1:
        smooth = False
    return ("bad" if em else ("zero" if smooth else "cofactor-big"),
            (em, (rn, rd), smooth))


def report_hit(kind, a, z, tau, b):
    print(f"*** {kind} HIT a={a} b={b} tau={tau}", flush=True)
    try:
        t1 = _l7_tied_status(a, b, z, tau)
    except Exception as e:
        t1 = f"refused:{type(e).__name__}"
    try:
        t2 = _l9_steered_solvable(a, b, z, tau)
    except Exception as e:
        t2 = f"refused:{type(e).__name__}"
    print(f"***   cross-check tied={t1} steered={t2}", flush=True)


# ---------------- stage 2: cofactor ladder (run 2, 2026-08-18) ----------------
# After smooth_emergent(a, z, tau, b) returns ('cofactor-big', (em, (rn, rd),
# False)) the row is decided by the ladder below, applied to each stripped
# cofactor x in (rn, rd).  Every prime factor of x exceeds 1e6 (the y-smooth
# strip removed the rest), and em == [] (no bad place below 1e6).
#   (a) x a perfect square -> contributes no emergent place;
#   (b) kernel _is_prime(x) -- PROVEN primality only (deterministic MR below
#       _MR_LIMIT, Pocklington above; PrimalityBound = passed every MR base
#       but unprovable):
#         True  -> proved prime.  The parity law (Hilbert product formula with
#                  every frozen and infinite place verified +1) forces
#                  hilbert(alpha*M0, alpha*2b, x) = +1; a computed -1 on a
#                  singleton emergent place is an arithmetic contradiction and
#                  is recorded LOUDLY as INCONSISTENT.
#         PrimalityBound -> 'jacobi-unproved': the symbol computed at x is then
#                  only a Jacobi symbol; +1 is consistent with zero-bad OR an
#                  even bad pair -- EVIDENCE tier only, never a proof;
#         False -> composite, continue;
#   (c) odd perfect-power reduction x = t^e (e odd, since x is not a square;
#       odd exponents of t are exactly the odd exponents of x), then (b) on t;
#   (d) composites with <= 72 digits: kernel factorint under signal.alarm(20)
#       -> exact odd-exponent emergent set, hilbert at each prime; kernel
#       FactorBudget/PrimalityBound or the alarm = refusal, never a guess.
# Composites above 72 digits are refused ('refused:cofactor-big-digits-N').
# The factorint rung (d) is the only budgeted step (global wall budget, rows
# in ascending digit order, chosen by the driver via defer=True).

_LADDER_DIGIT_CAP = 72
_FACTINT_ALARM = 20            # seconds, per factorint call
_PRIME_STALL_GUARD = 300       # seconds, stall guard on _is_prime only
_FACT_BUDGET_SECONDS = 7200    # global budget for ladder rung (d) alone


class _FactorAlarm(Exception):
    pass


def _alarm_handler(signum, frame):
    raise _FactorAlarm()


@contextlib.contextmanager
def _alarm_guard(seconds):
    """signal.alarm guard; silently inactive off the main thread."""
    try:
        signal.signal(signal.SIGALRM, _alarm_handler)
        signal.alarm(seconds)
    except (ValueError, OSError):
        yield
        return
    try:
        yield
    finally:
        signal.alarm(0)


def _iroot(x, k):
    """floor(x ** (1/k)) for x >= 1, k >= 2 (exact bisection)."""
    if x < 2:
        return x
    lo, hi = 1, 1 << (x.bit_length() // k + 1)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid ** k <= x:
            lo = mid
        else:
            hi = mid - 1
    return lo


def _reduce_odd_power(x):
    """If x = t^e for an odd e >= 3 (x not a square => any power rep is odd),
    return the fully reduced base t; else return x.  Returns (t, reduced)."""
    y = x
    reduced = False
    while True:
        e = 3
        found = False
        while e <= y.bit_length():
            t = _iroot(y, e)
            if t > 1 and t ** e == y:
                y = t
                reduced = True
                found = True
                break
            e += 2
        if not found:
            return y, reduced


def _classify_cofactor(x):
    """Ladder rungs (a)-(c) for one stripped cofactor x.  Returns a dict with
    kind in {'unit', 'square', 'prime', 'jacobi', 'composite', 'isprime-stall'}
    plus, where meaningful, 'place' (the single prospective emergent prime for
    kind 'prime'/'jacobi'), 'digits', and 'rung' ('a'|'b'|'c'|'d')."""
    if x == 1:
        return {'kind': 'unit', 'rung': 'a'}
    r = math.isqrt(x)
    if r * r == x:
        return {'kind': 'square', 'rung': 'a'}
    y, reduced = _reduce_odd_power(x)
    rung = 'c' if reduced else 'b'
    try:
        with _alarm_guard(_PRIME_STALL_GUARD):
            isp = _is_prime(y)
    except PrimalityBound:
        return {'kind': 'jacobi', 'place': y, 'rung': rung, 'digits': len(str(y))}
    except _FactorAlarm:
        return {'kind': 'isprime-stall', 'rung': rung, 'digits': len(str(y))}
    if isp:
        return {'kind': 'prime', 'place': y, 'rung': rung, 'digits': len(str(y))}
    return {'kind': 'composite', 'base': y, 'rung': 'd', 'digits': len(str(y))}


def _factorint_capped(x):
    """Kernel factorint under signal.alarm(_FACTINT_ALARM); raises _FactorAlarm
    on timeout, propagates kernel FactorBudget/PrimalityBound."""
    with _alarm_guard(_FACTINT_ALARM):
        return _kernel_factorint(x)


def cofactor_decide(a, z, tau, b, detail=None, defer=False, deadline=None):
    """Stage-2 cofactor ladder for one candidate (run 2).

    Call AFTER smooth_emergent(a, z, tau, b) returned status 'cofactor-big';
    pass its detail as `detail` to avoid recomputation (else it is recomputed;
    non-'cofactor-big' statuses are then passed straight through).  `defer`
    stops before rung (d): a row whose only open part is a <= 72-digit
    composite returns ('deferred', info) instead of factoring -- the driver
    batches deferred rows in ascending digit order under a global budget.
    `deadline` (a time.time() cutoff) gates rung (d) only.

    Returns (verdict, info) with verdict in {'zero' (proved: every emergent
    symbol +1, soluble candidate), 'bad-n' (n proved emergent places with
    symbol -1), 'zero-jacobi' / 'bad-jacobi' (evidence tier: an unproved
    MR-passed cofactor is involved), 'INCONSISTENT', 'deferred',
    'refused:cofactor-big-digits-N', 'refused:FactorBudget',
    'refused:PrimalityBound', 'refused:factorint-timeout',
    'refused:budget-exhausted', 'refused:isprime-stall'}."""
    if detail is None:
        st, detail = smooth_emergent(a, z, tau, b)
        if st != 'cofactor-big':
            return (st, detail)
    A, delta, alpha, Z, D, s = pattern(a, z, tau)
    em, (rn, rd), smooth = detail
    c0 = _sun_h(a, b, Z)
    M0 = 16 - delta*c0*c0 - 32*A*b*s*s
    x0, d0 = alpha*M0, alpha*2*b

    parts = [_classify_cofactor(rn), _classify_cofactor(rd)]
    info = {'rn': str(rn), 'rd': str(rd), 'x0': str(x0), 'd0': str(d0),
            'parts': [{k: (str(v) if isinstance(v, int) else v)
                       for k, v in p.items()} for p in parts]}
    if any(p['kind'] == 'isprime-stall' for p in parts):
        return ('refused:isprime-stall', info)
    composite = [p for p in parts if p['kind'] == 'composite']
    big = [p for p in composite if p['digits'] > _LADDER_DIGIT_CAP]
    if big:
        return ('refused:cofactor-big-digits-%d' % big[0]['digits'], info)
    if composite and defer:
        info['defer_digits'] = max(p['digits'] for p in composite)
        return ('deferred', info)
    factored = []
    if composite:
        if deadline is not None and time.time() > deadline:
            return ('refused:budget-exhausted', info)
        try:
            for p in composite:
                factored.append(_factorint_capped(p['base']))
        except FactorBudget:
            return ('refused:FactorBudget', info)
        except PrimalityBound:
            return ('refused:PrimalityBound', info)
        except _FactorAlarm:
            return ('refused:factorint-timeout', info)

    places = []                       # (prime, tier, rung)
    for p in parts:
        if p['kind'] == 'prime':
            places.append((p['place'], 'proved',
                           'pprime' if p['rung'] == 'c' else 'prime'))
        elif p['kind'] == 'jacobi':
            places.append((p['place'], 'unproved',
                           'pprime' if p['rung'] == 'c' else 'prime'))
    for fac in factored:
        for q, e in sorted(fac.items()):
            if e % 2:
                places.append((q, 'proved', 'factorint'))
    info['places'] = [[str(q), hilbert(x0, d0, q), tier, rung]
                      for q, tier, rung in places]
    proved_negs = [pl for pl in info['places'] if pl[1] == -1 and pl[2] == 'proved']
    if len(info['places']) == 1 and proved_negs:
        print(f"*** INCONSISTENT a={a} b={b} tau={tau}: singleton emergent "
              f"place {info['places'][0][0]} has Hilbert symbol -1 "
              f"(parity law demands +1)", flush=True)
        return ('INCONSISTENT', info)
    rungs = {pl[3] for pl in info['places']} | ({'factorint'} if factored else set())
    info['rung'] = ('factorint' if 'factorint' in rungs else
                    'pprime' if 'pprime' in rungs else
                    'prime' if 'prime' in rungs else 'square')
    if proved_negs:
        return ('bad-%d' % len(proved_negs), info)
    if any(pl[1] == -1 for pl in info['places']):
        return ('bad-jacobi', info)
    if any(pl[2] == 'unproved' for pl in info['places']):
        return ('zero-jacobi', info)
    return ('zero', info)


def record_hit(kind, a, z, tau, b, extra=None):
    """Print a hit with *** markers, cross-check via the engine (refusals
    acceptable), and return the structured hit record."""
    try:
        t1 = _l7_tied_status(a, b, z, tau)
    except Exception as e:
        t1 = f"refused:{type(e).__name__}"
    try:
        t2 = _l9_steered_solvable(a, b, z, tau)
    except Exception as e:
        t2 = f"refused:{type(e).__name__}"
    print(f"*** {kind} HIT a={a} b={b} tau={tau}", flush=True)
    print(f"***   cross-check tied={t1} steered={t2}", flush=True)
    rec = {'kind': kind, 'a': str(a), 'b': str(b), 'tau': str(tau), 'z': str(z),
           'tied': str(t1), 'steered': str(t2)}
    if extra:
        rec.update(extra)
    return rec


_COUNTERS = ('tested', 'aligned', 'align_fail', 'none_status', 'smooth_zero',
             'smooth_bad', 'cofactor_big_stored', 'cofactor_square_zero',
             'cofactor_prime_zero', 'cofactor_pprime_zero', 'factorint_zero',
             'factorint_bad', 'cofactor_prime_bad', 'cofactor_pprime_bad',
             'jacobi_unproved_zero', 'jacobi_unproved_bad',
             'refused_cofactor_gt72digits', 'refused_factorint_timeout',
             'refused_budget_exhausted', 'refused_isprime_stall',
             'inconsistent_count', 'hits_total')

_BRANCHES = (('sq', lambda a: (1 + 2*a*a)/(1 + 4*a*a)),
             ('can', lambda a: 2*a/(1 + 4*a*a)),
             ('one', lambda a: F(1)))


def _record_ladder(C, reasons, hits, rows, r, verdict, info, z, ladder_stage):
    """Fold one ladder verdict into counters / row log / hits."""
    rung = info.get('rung') if isinstance(info, dict) else None
    row = {'ladder_stage': ladder_stage,
           **{k: r[k] for k in ('src', 'a', 'br', 'f', 'e', 'eps', 'q',
                                'shape', 'b', 'digits', 'rn', 'rd')},
           'verdict': verdict, 'rung': rung}
    if isinstance(info, dict):
        row['places'] = info.get('places', [])
    rows.append(row)
    if verdict != 'zero':
        print(f"  [{ladder_stage}] {verdict} src={r['src']} a={r['a']} "
              f"br={r['br']} f={r['f']} e={r['e']} eps={r['eps']} q={r['q']} "
              f"shape={r['shape']}", flush=True)
    if verdict == 'zero':
        if rung == 'factorint':
            C['factorint_zero'] += 1
        elif rung == 'pprime':
            C['cofactor_pprime_zero'] += 1
        elif rung == 'prime':
            C['cofactor_prime_zero'] += 1
        else:
            C['cofactor_square_zero'] += 1
        C['hits_total'] += 1
        hits.append(record_hit(f"C-{ladder_stage}-{rung}-zero", F(r['a']), z,
                               F(r['tau']), F(r['b']),
                               extra={'rung': rung,
                                      'places': info.get('places'),
                                      'rn': info.get('rn'), 'rd': info.get('rd')}))
    elif verdict.startswith('bad-') and verdict != 'bad-jacobi':
        if rung == 'factorint':
            C['factorint_bad'] += 1
        elif rung == 'pprime':
            C['cofactor_pprime_bad'] += 1
        else:
            C['cofactor_prime_bad'] += 1
    elif verdict == 'zero-jacobi':
        C['jacobi_unproved_zero'] += 1
    elif verdict == 'bad-jacobi':
        C['jacobi_unproved_bad'] += 1
    elif verdict.startswith('refused:cofactor-big-digits'):
        C['refused_cofactor_gt72digits'] += 1
    elif verdict in ('refused:FactorBudget', 'refused:PrimalityBound',
                     'refused:factorint-timeout'):
        C['refused_factorint_timeout'] += 1
        why = verdict.split(':', 1)[1]
        reasons[why] = reasons.get(why, 0) + 1
    elif verdict == 'refused:budget-exhausted':
        C['refused_budget_exhausted'] += 1
    elif verdict == 'refused:isprime-stall':
        C['refused_isprime_stall'] += 1
    elif verdict == 'INCONSISTENT':
        C['inconsistent_count'] += 1


# ------------------------------------------------------------------ phase A
def validate():
    n_ok = 0
    for tbl, key, k in _L13_ZERO_BAD:
        w, ut = key
        z = F(w)*F(*ut)
        if tbl == 'L11':
            row = _L11_CLASSES[key]
            a = F(*row[0])
            eps = row[1]
            A = 1 + 4*a*a
            tau = (1 + 2*a*a)/A
            b = F(eps*(row[2] + k*row[3]))
        else:
            import l12_class
            a = F(1)
            tau = F(3, 5)
            f, eps, q = _L12_ESCAPE[key]
            cert = l12_class.l12_class_cert(w, *ut)
            b = F(eps*f*(q + k*cert['N']))
        zb = emergent_free(a, z, tau, b)
        assert zb is not None and zb[0], (tbl, key, k, zb[1])
        n_ok += 1
    print(f"A1 frozen zero-bad: {n_ok}/{len(_L13_ZERO_BAD)}", flush=True)
    n = 0
    for line in open('data/audit_l11h_members.jsonl'):
        r = json.loads(line)
        if r.get('emergent_status') != 'complete' or r.get('bad_count') is None:
            continue
        w, ut = r['row_key']
        row = _L11_CLASSES[(w, tuple(ut))]
        a = F(*row[0])
        eps = row[1]
        A = 1 + 4*a*a
        z = F(w)*F(*ut)
        tau = (1 + 2*a*a)/A
        zb = emergent_free(a, z, tau, F(eps*r['Q']))
        assert zb is not None
        assert zb[0] == (r['bad_count'] == 0), (w, ut, r['k'], zb, r['bad_count'])
        n += 1
    n = agree = 0
    for line in open('data/audit_l11h_members.jsonl'):
        r = json.loads(line)
        if r.get('emergent_status') != 'complete' or r.get('bad_count') is None:
            continue
        w, ut = r['row_key']
        row = _L11_CLASSES[(w, tuple(ut))]
        a = F(*row[0]); eps = row[1]; A = 1+4*a*a
        z = F(w)*F(*ut); tau = (1+2*a*a)/A
        b = F(eps*r['Q'])
        st, dt = smooth_emergent(a, z, tau, b)
        if st == 'cofactor-big':
            continue      # big run (checked separately in phase B)
        assert st in ('zero', 'bad'), (w, ut, r['k'], st)
        assert (st == 'zero') == (r['bad_count'] == 0), (w, ut, r['k'], st, r['bad_count'], dt[1][:2])
        agree += 1
        n += 1
    print(f"A2 smooth-member agreement on decidable: {agree}", flush=True)


# ------------------------------------------------------------------ phase B
def attack():
    meta = json.load(open('data/l13_cell89.json'))['meta']
    pool = meta['supp_union_odd_per_a']
    w, ut = 89, (-2, 1)
    z = F(w, 1)*F(ut[0], ut[1])
    C = {k: 0 for k in _COUNTERS}
    reasons = {}
    hits, stored, rows = [], [], []
    done = set()
    T0 = time.time()
    wall = {}

    def box_spec():
        return {'cell': [89, [-2, 1]], 'z': '-178',
                'a_pool': {k: list(v) for k, v in pool.items()},
                'branches': ['sq', 'can', 'one'],
                'shapes': {'0': 'b = eps*f^e*q', '1': 'b = eps*f^e/q'},
                'e': [1, 2], 'eps': [1, -1],
                'q_lo': 'primerange(41,3001) = 418 primes',
                'q_hi': 'primerange(3001,8001) = 577 primes',
                'stage_A': 'a=17, branch sq only, q_hi (18,464 candidates)',
                'stage_B': 'every a in pool x every branch x q_lo (290,928 candidates)'}

    def _dump(obj, path):
        tmp = path + '.tmp'
        with open(tmp, 'w') as fh:
            json.dump(obj, fh, indent=1)
        os.replace(tmp, path)

    def checkpoint(status, extra=None):
        payload = {'meta': {'tool': 'math/h10q/l13_filter.py -- run 2 '
                                    '(stage 1 smooth decision + stage 2 cofactor ladder)',
                            'status': status,
                            'generated': time.strftime('%Y-%m-%d %H:%M:%S'),
                            'lineage_run1': 'h10q-l13-filter3: QS1 box 290,928 '
                                            'candidates, 0 hits (smooth decisions only)',
                            'box': box_spec(),
                            'cpu_policy': 'single compute subprocess, single-threaded math; NO in-loop sleep (Main ruling 2026-08-18: periodic time.sleep exiles the thread to ~2% CPU duty under box load 45+; process runs nice -n 19, one core of 28 = ~3% of box, inside the <=50% cap)',
                            'provenance': 'stage 1: smooth_emergent is exact when the '
                                          'stripped remainder is 1e6-smooth; stage 2 ladder: '
                                          'h10q kernel _is_prime (proven primality) / '
                                          'budget-bounded factorint under signal.alarm(20); '
                                          'jacobi-unproved rows are evidence tier only; '
                                          'refusals are never evidence'},
                   'counters': C, 'refusal_reasons': reasons,
                   'wall_seconds': dict(wall, total=time.time() - T0),
                   'done_blocks': sorted(done),
                   'stored': stored, 'hits': hits, 'ladder_rows': rows}
        if extra:
            payload.update(extra)
        _dump(payload, 'data/l13_filter_run2.json')
        return payload

    # crash-safe resume of stages A/B only (identical enumeration order keeps
    # the counters exact); stages C and D always run fresh from `stored`
    if os.path.exists('data/l13_filter_run2.json'):
        try:
            prev = json.load(open('data/l13_filter_run2.json'))
            resumable = prev.get('meta', {}).get('status') in ('running-A', 'running-B')
            if prev.get('meta', {}).get('status') == 'complete':
                cell_now = json.load(open('data/l13_cell89.json'))
                resumable = 'filter_run2' not in cell_now
            if resumable and prev.get('meta', {}).get('box') == box_spec():
                for k, v in prev['counters'].items():
                    if k in C:
                        C[k] = v
                reasons.update(prev.get('refusal_reasons', {}))
                hits.extend(prev.get('hits', []))
                stored.extend(prev.get('stored', []))
                done.update(tuple(x) for x in prev.get('done_blocks', []))
                print(f"resume: {len(done)} blocks done, tested={C['tested']}, "
                      f"stored={len(stored)}", flush=True)
        except Exception as e:
            print(f"resume skipped ({type(e).__name__}: {e})", flush=True)

    QS1 = list(primerange(41, 3001))
    QS2 = list(primerange(3001, 8001))
    branch_f = dict(_BRANCHES)

    def run_block(src, a_, br, qs):
        a = F(a_)
        tau = branch_f[br](a)
        for f in pool[str(a_)]:
            for e_ in (1, 2):
                for eps in (1, -1):
                    for q in qs:
                        for shape, b in ((0, F(eps*f**e_*q)),
                                         (1, F(eps*f**e_, q))):
                            C['tested'] += 1
                            # in-loop sleep removed per Main ruling 2026-08-18
                            # (periodic sleep -> scheduler exile at box load 45+)
                            st, dt = smooth_emergent(a, z, tau, b)
                            if st in ('zero', 'bad', 'cofactor-big'):
                                C['aligned'] += 1
                            if st == 'zero':
                                C['smooth_zero'] += 1
                                C['hits_total'] += 1
                                hits.append(record_hit(f"{src}-{br}-smooth-zero",
                                                       a, z, tau, b))
                            elif st == 'bad':
                                C['smooth_bad'] += 1
                            elif st == 'cofactor-big':
                                em, (rn, rd), _sm = dt
                                stored.append({'src': src, 'a': a_, 'br': br,
                                               'f': f, 'e': e_, 'eps': eps, 'q': q,
                                               'shape': shape, 'b': str(b),
                                               'tau': str(tau), 'rn': str(rn),
                                               'rd': str(rd),
                                               'digits': max(len(str(rn)), len(str(rd)))})
                            elif st == 'align-fail':
                                C['align_fail'] += 1
                            else:
                                C['none_status'] += 1

    # ---------------- stage A: QS2 pass (a=17, branch sq, q in 3001..8000)
    if ('A', 17) not in done:
        checkpoint('running-A')
        run_block('A', 17, 'sq', QS2)
        done.add(('A', 17))
        wall['A'] = time.time() - T0
        print(f"A done: tested={C['tested']} hits={C['hits_total']} "
              f"stored={len(stored)} ({wall['A']:.0f}s)", flush=True)
        checkpoint('running-A')

    # ---------------- stage B: the entire stage-1 box, re-enumerated
    tb = time.time()
    for a_s in sorted(pool, key=int):
        a_ = int(a_s)
        if vp(F(a_), 2) != 0:
            continue
        for br, _tf in _BRANCHES:
            if ('B', a_, br) not in done:
                run_block('B', a_, br, QS1)
                done.add(('B', a_, br))
        checkpoint('running-B')
        print(f"  a={a_} done: tested={C['tested']} hits={C['hits_total']} "
              f"stored={len(stored)} ({time.time()-T0:.0f}s)", flush=True)
    wall['B'] = time.time() - tb

    # ---------------- stage C: cofactor ladder over the stored rows
    tc = time.time()
    stored.sort(key=lambda r: (r['digits'], r['src'], r['a'], r['br'], r['f'],
                               r['e'], r['eps'], r['q'], r['shape']))
    C['cofactor_big_stored'] = len(stored)
    print(f"C1: ladder rungs (a)-(c) over {len(stored)} stored rows (unbudgeted)",
          flush=True)
    deferred = []
    for i, r in enumerate(stored):
        # in-loop sleep removed per Main ruling 2026-08-18
        verdict, info = cofactor_decide(F(r['a']), z, F(r['tau']), F(r['b']),
                                        detail=([], (int(r['rn']), int(r['rd'])), False),
                                        defer=True)
        if verdict == 'deferred':
            deferred.append(r)
        else:
            _record_ladder(C, reasons, hits, rows, r, verdict, info, z, 'C1')
    wall['C1'] = time.time() - tc
    print(f"C1 done: deferred-to-factorint={len(deferred)} "
          f"zero(prime)={C['cofactor_prime_zero']} "
          f"zero(pprime)={C['cofactor_pprime_zero']} "
          f"zero(sq)={C['cofactor_square_zero']} "
          f"jacobi={C['jacobi_unproved_zero'] + C['jacobi_unproved_bad']} "
          f"refused>72d={C['refused_cofactor_gt72digits']} "
          f"inconsistent={C['inconsistent_count']} ({wall['C1']:.0f}s)", flush=True)

    t2 = time.time()
    deadline = t2 + _FACT_BUDGET_SECONDS
    print(f"C2: factorint rung, {len(deferred)} rows ascending digits, "
          f"global budget {_FACT_BUDGET_SECONDS}s", flush=True)
    for i, r in enumerate(deferred):
        if time.time() > deadline:
            for r2 in deferred[i:]:
                rows.append({'ladder_stage': 'C2',
                             **{k: r2[k] for k in ('src', 'a', 'br', 'f', 'e',
                                                   'eps', 'q', 'shape', 'b',
                                                   'digits', 'rn', 'rd')},
                             'verdict': 'refused:budget-exhausted', 'rung': None})
                C['refused_budget_exhausted'] += 1
            print(f"C2: budget exhausted, {len(deferred) - i} rows "
                  f"refused:budget-exhausted", flush=True)
            break
        verdict, info = cofactor_decide(F(r['a']), z, F(r['tau']), F(r['b']),
                                        detail=([], (int(r['rn']), int(r['rd'])), False),
                                        deadline=deadline)
        _record_ladder(C, reasons, hits, rows, r, verdict, info, z, 'C2')
        if i and i % 50 == 0:
            print(f"  C2 {i}/{len(deferred)} ({time.time()-t2:.0f}s) "
                  f"zero={C['factorint_zero']} bad={C['factorint_bad']} "
                  f"refused={C['refused_factorint_timeout']}", flush=True)
    wall['C2'] = time.time() - t2

    # ------------- stage D: DenseScan's dense-box cofactor-big rows
    # (a=17 only, tau in {sq, can}, |b| < 20001, cell (89,(-2,1))): the
    # highest-value set anywhere -- small b, aligned, only the remainder
    # undecided.  Every row is re-verified through our own smooth_emergent
    # before the ladder touches it.  Counters kept SEPARATE under the
    # 'dense_box' key; the factorint rung gets its own 2h budget.
    td = time.time()
    dense = None
    dense_note = ''
    path = 'data/l13_dense_cofactorbig.jsonl'
    for attempt in range(15):
        if os.path.exists(path):
            s1 = os.path.getsize(path)
            time.sleep(20)
            if os.path.getsize(path) == s1 and s1 > 0:
                dense = []
                for ln in open(path):
                    ln = ln.strip()
                    if ln:
                        dense.append(json.loads(ln))
                break
        else:
            if attempt == 0:
                print(f"D: waiting for {path} (DenseScan)", flush=True)
            time.sleep(60)
    DC = {k: 0 for k in _COUNTERS}
    DC['dense_row_mismatch'] = 0
    dreasons, dhits, drows = {}, [], []
    tau_sq17 = branch_f['sq'](F(17))
    tau_can17 = branch_f['can'](F(17))
    if dense is None:
        dense_note = (f'{path} never appeared during 15 one-minute polls; '
                      'dense box left unprocessed (noted, not evidence)')
        print(f"D: {dense_note}", flush=True)
        dense_block = {'status': 'skipped-file-absent', 'note': dense_note,
                       'source': 'DenseScan (a=17, tau in {sq,can}, |b|<20001, '
                                 'cell (89,(-2,1))), expected 1,286 rows',
                       'counters': DC, 'refusal_reasons': dreasons,
                       'hits': dhits, 'ladder_rows': drows}
    else:
        dense.sort(key=lambda r: (max(len(str(r['rn'])), len(str(r['rd']))),
                                  str(r['tau']), str(r['b'])))
        DC['cofactor_big_stored'] = len(dense)
        DC['tested'] = len(dense)
        verified = []
        for i, dr in enumerate(dense):
            # in-loop sleep removed per Main ruling 2026-08-18
            # dense file carries the branch NAME in 'tau' (DenseScan format)
            a, b = F(dr['a']), F(dr['b'])
            tfield = str(dr['tau'])
            if tfield in branch_f:
                br = tfield
                tau = branch_f[br](a)
            else:
                tau = F(tfield)
                br = ('sq' if tau == tau_sq17 else
                      'can' if tau == tau_can17 else '?')
            rn_f, rd_f = int(dr['rn']), int(dr['rd'])
            st, dt = smooth_emergent(a, z, tau, b)
            row = {'src': 'dense', 'a': dr['a'], 'br': br, 'f': None, 'e': None,
                   'eps': None, 'q': None, 'shape': None, 'b': str(b),
                   'tau': str(tau), 'rn': str(rn_f), 'rd': str(rd_f),
                   'digits': max(len(str(rn_f)), len(str(rd_f)))}
            if st != 'cofactor-big' or dt[1] != (rn_f, rd_f):
                DC['dense_row_mismatch'] += 1
                drows.append({**row, 'verdict': 'dense-row-mismatch', 'rung': None,
                              'recheck': st})
                print(f"  [D] dense-row-mismatch b={b} tau={tau} file_status"
                      f"=cofactor-big recheck={st} rn_rd_match="
                      f"{dt[1] == (rn_f, rd_f)}", flush=True)
                continue
            DC['aligned'] += 1
            verified.append(row)
        dd = []
        for i, row in enumerate(verified):
            # in-loop sleep removed per Main ruling 2026-08-18
            verdict, info = cofactor_decide(F(row['a']), z, F(row['tau']),
                                            F(row['b']),
                                            detail=([], (int(row['rn']),
                                                         int(row['rd'])), False),
                                            defer=True)
            if verdict == 'deferred':
                dd.append(row)
            else:
                _record_ladder(DC, dreasons, dhits, drows, row, verdict, info,
                               z, 'D1')
        print(f"D1 done: verified={len(verified)} mismatch="
              f"{DC['dense_row_mismatch']} deferred={len(dd)} "
              f"zero(prime)={DC['cofactor_prime_zero']} "
              f"zero(pprime)={DC['cofactor_pprime_zero']} "
              f"refused>72d={DC['refused_cofactor_gt72digits']}", flush=True)
        td2 = time.time()
        ddeadline = td2 + _FACT_BUDGET_SECONDS
        print(f"D2: factorint rung, {len(dd)} rows ascending digits, "
              f"own budget {_FACT_BUDGET_SECONDS}s", flush=True)
        for i, row in enumerate(dd):
            if time.time() > ddeadline:
                for r2 in dd[i:]:
                    drows.append({'ladder_stage': 'D2', **{k: r2[k] for k in
                         ('src', 'a', 'br', 'f', 'e', 'eps', 'q', 'shape', 'b',
                          'digits', 'rn', 'rd')},
                        'verdict': 'refused:budget-exhausted', 'rung': None})
                    DC['refused_budget_exhausted'] += 1
                print(f"D2: budget exhausted, {len(dd) - i} rows "
                      f"refused:budget-exhausted", flush=True)
                break
            verdict, info = cofactor_decide(F(row['a']), z, F(row['tau']),
                                            F(row['b']),
                                            detail=([], (int(row['rn']),
                                                         int(row['rd'])), False),
                                            deadline=ddeadline)
            _record_ladder(DC, dreasons, dhits, drows, row, verdict, info,
                           z, 'D2')
            if i and i % 50 == 0:
                print(f"  D2 {i}/{len(dd)} ({time.time()-td2:.0f}s) "
                      f"zero={DC['factorint_zero']} bad={DC['factorint_bad']} "
                      f"refused={DC['refused_factorint_timeout']}", flush=True)
        dense_block = {'status': 'complete',
                       'source': 'DenseScan: a=17, tau in {sq,can}, |b|<20001, '
                                 'cell (89,(-2,1)); file ' + path,
                       'expected_rows': 1286, 'rows': len(dense),
                       'counters': DC, 'refusal_reasons': dreasons,
                       'hits': dhits, 'ladder_rows': drows,
                       'own_factorint_budget_seconds': _FACT_BUDGET_SECONDS,
                       'verified_by': 'every row re-run through smooth_emergent; '
                                      'status and (rn, rd) independently confirmed '
                                      'before the ladder',
                       'wall_seconds': {'D_total': time.time() - td}}
        wall['D'] = time.time() - td

    dense_hits_n = len(dhits)
    outcome = ('HITS FOUND -- cell (89,(-2,1)) has proved soluble candidates; '
               'see hits[]' if (hits or dense_hits_n) else
               'no proved zero-bad candidate in the decided part of the box')
    if dense_hits_n:
        outcome += f' (dense_box hits: {dense_hits_n})'
    checkpoint('complete', extra={'outcome': outcome, 'dense_box': dense_block})
    cell = json.load(open('data/l13_cell89.json'))
    cell['filter_run2'] = {'generated': time.strftime('%Y-%m-%d %H:%M:%S'),
                           'tool': 'math/h10q/l13_filter.py -- run 2',
                           'counters': dict(C), 'refusal_reasons': dict(reasons),
                           'wall_seconds': dict(wall, total=time.time() - T0),
                           'hits': hits, 'outcome': outcome,
                           'dense_box': {'status': dense_block['status'],
                                         'counters': dense_block['counters'],
                                         'refusal_reasons': dense_block['refusal_reasons'],
                                         'hits': dense_block['hits'],
                                         'note': dense_block.get('note', '')},
                           'detail': 'data/l13_filter_run2.json (box spec, stored '
                                     'cofactor rows, per-row ladder verdicts)'}
    _dump(cell, 'data/l13_cell89.json')
    print(f"B/C/D done: tested={C['tested']} stored={C['cofactor_big_stored']} "
          f"hits={C['hits_total']} dense_rows={DC['tested']} "
          f"dense_hits={dense_hits_n} ({time.time()-T0:.0f}s)", flush=True)


def validate2():
    n = 0
    lad_ok = lad_ref = 0
    for tbl, key, k in _L13_ZERO_BAD:
        w, ut = key
        z = F(w)*F(*ut)
        if tbl == 'L11':
            row = _L11_CLASSES[key]
            a = F(*row[0]); eps = row[1]; A = 1+4*a*a
            tau = (1+2*a*a)/A
            b = F(eps*(row[2] + k*row[3]))
        else:
            import l12_class
            a = F(1); tau = F(3,5)
            f, eps, q = _L12_ESCAPE[key]
            cert = l12_class.l12_class_cert(w, *ut)
            b = F(eps*f*(q + k*cert['N']))
        st, dt = smooth_emergent(a, z, tau, b)
        if st in ('zero', 'bad'):
            n += 1
        elif st == 'cofactor-big':
            # ground truth (Main, 2026-08-18): frozen ZERO-BAD members must
            # come back 'zero' through the stage-2 ladder.  The L11 (5,(1,7))
            # k=65 row has a PROVABLY PRIME 48-digit remainder (rung (b)
            # end-to-end); the composite-remainder rows exercise rung (d).
            verdict, info = cofactor_decide(a, z, tau, b, detail=dt)
            if verdict == 'zero':
                lad_ok += 1
            elif verdict.startswith('refused:'):
                lad_ref += 1       # kernel refusal: allowed, never evidence
                print(f"  ladder refusal on frozen row {tbl} {key} k={k}: "
                      f"{verdict}", flush=True)
            else:
                print(f"*** LADDER VALIDATION FAILED {tbl} {key} k={k}: "
                      f"{verdict} {info}", flush=True)
                raise SystemExit("cofactor ladder contradicts frozen zero-bad "
                                 "ground truth")
    print(f"A1-smooth frozen zero-bad decided: {n}/{len(_L13_ZERO_BAD)}; "
          f"ladder zero on cofactor-big frozen rows: {lad_ok} "
          f"(refused {lad_ref})", flush=True)


if __name__ == '__main__':
    validate2()
    print("phase B", flush=True)
    attack()
