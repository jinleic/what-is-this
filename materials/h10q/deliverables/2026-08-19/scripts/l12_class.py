"""L12b VERIFIED CLASS ALIGNMENT for the 103 escape rows of h10q._L12_ESCAPE.

Setup.  a = 1 (A = 5), square branch tau = (1+2a^2)/A = 3/5: delta = -4/5,
alpha = -delta*A = 4, a PERFECT SQUARE; s_tie = (a-1)/2 = 0.  A row
(w, (u1, u2)) of h10q._L12_ESCAPE records the cell z = w*u1/u2 and a base
b0 = eps*f*q with f | numerator(z) (eps = +1 in all 103 rows).  Members of
the aligned class are b = eps*f*Q with Q prime, Q == q (mod N), where
N = lcm(8, prod_{p in S} p^k_p, 4A) and k_p = _l10_exponent(P, b0, p).

The minimal controlled set is
    S_min = {2,3,5,7} u supp(alpha) u supp(delta) u {f} = {2,3,5,7,f}
(supp(alpha) = {2}, supp(delta) = {2,5}); this is the _l10_class_cert frozen
set plus the primes of the FIXED factor f of b.  EXACT SCOPE of `ok=True`
(STEP (i) ONLY):

  * FROZEN SIDE: for every p in S_min, the exponent lemma with modulus
    p^k_p | N gives v_p(P(b)) = v_p(P(b0)) and P(b)/P(b0) a square in Q_p,
    so (x0, d0)_p is CONSTANT on the class and equal to its verified +1
    base value (symbols and valuations are re-verified on sampled members;
    b - b0 = eps*f*k*N is asserted 0 mod p^k_p).
  * WILD place Q: v_Q(d0) = 1, v_Q(x0) = -4 (v_Q(P(b)) = 0 because
    P(b) == -delta*Z^4 (mod Q) and Q is coprime to delta*Z), unit(x) == A
    mod squares, hence (x0,d0)_Q = (A|Q) = (A|q) = +1 along the class
    (5 | N fixes the residue mod 5).  Verified directly per member.
  * REAL place: d0 = 8b has sign(eps*f*Q) = +1 and, past the Cauchy bound
    |b| = f*|Q| > Q0 = 1 + max|c_i/c_lead|, sign(x0) is the asymptotic one
    sx = sign(alpha)*sign(lead)*eps^8 = +1, so (x0,d0)_R = +1 on the class
    (asserted f*Q > Q0 AND verified exactly on every sampled member).

NOT COVERED -- STEP (ii), OUTSIDE the class claim: the EMERGENT places,
odd primes l not in S_min u {Q} with v_l(P(b)) ODD.  There v_l(d0) = 0 and
(x0,d0)_l = (unit(d0)|l) may be -1.  They are PROBED (exact factorization
of P(b) where the factorization budget permits; exact valuation dichotomy
at every prime l < 10000 on sampled members) and REPORTED, never claimed.
FactorBudget / PrimalityBound refusals of the proven engine are logged and
counted; they are never evidence either way.  Consequently ok=True is NOT a
claim of global conic solubility (by reciprocity that would additionally
need the emergent -1 count to be even).

EXCLUDED MEMBERS (the class claim is void at these finitely many primes):
Q dividing  D_z * numerator(z) * a * delta,  i.e. with a = 1, delta = -4/5:
    excluded = supp(D_z) u supp(numerator(z)) u {2,5},
stored per row in the certificate under 'excluded'.  At such a Q the
uniform valuation facts (v_Q(P(b)) = 0, unit(P) == A mod squares) fail.
Note the base q itself lies in the excluded set for the 5 rows with 41 | D_z
((11,(7,3)), (19,(7,3)), (43,(1,7)), (43,(-2,7)), (47,(2,1))): the base
wild symbol there is a BASE-POINT fact verified directly, and class
transport to members is unaffected because Q == q (mod N) with q | N
false, so no member Q = q + kN (k >= 1) ever equals q.
"""

import math
from fractions import Fraction as F

from h10q import (OO, FactorBudget, PrimalityBound, _L12_ESCAPE, _is_prime,
                  _l10_P, _l10_exponent, _l10_supp, _sun_h, factorint,
                  hilbert, legendre, primerange, unit_part, vp)

# Branch constants: a = 1, A = 5, tau = 3/5, delta = -4/5, alpha = 4.
_A = F(5)
_DELTA = F(-4, 5)
_ALPHA = F(4)


def _eval_P(P, b):
    """Exact evaluation of the coefficient list P (lowest first) at b."""
    return sum(c * b ** i for i, c in enumerate(P) if c)


def l12_class_cert(w, u1, u2):
    """Aligned-class certificate for the escape-family member b = eps*f*Q.

    Mirroring h10q._l10_class_cert semantics: returns None if the row or its
    base point is degenerate for the class argument, else a dict:
        ok  -- every frozen symbol p in S_min, the wild symbol (A|q), and the
               asymptotic real symbol are all +1 (class-constancy theorems
               above; emergent places are step (ii), NOT covered);
        N   -- the alignment modulus lcm(8, prod p^k_p, 4A);
        S   -- the minimal controlled (frozen) set {2,3,5,7,f};
        ks  -- the per-prime exponent-lemma moduli k_p;
        syms -- the verified base/frozen + wild + asymptotic real symbols;
        Q0  -- the Cauchy bound for the asymptotic real sign;
        f, eps, q -- the escape data (b = eps*f*Q, Q == q mod N);
        excluded  -- the finite excluded-member set
               {Q : Q | D_z * numerator(z) * a * delta}.
    """
    row = _L12_ESCAPE.get((w, (u1, u2)))
    if row is None:
        return None
    f, eps, q = row
    a = F(1)
    A = 1 + 4 * a * a
    tau = (1 + 2 * a * a) / A
    delta = 1 - A * tau * tau
    if delta == 0:
        return None
    alpha = -delta * A
    z = F(w) * F(u1, u2)
    Z = z ** 3
    D = 1 - Z - a * a * Z * Z
    if a == 0 or z == 0 or D == 0:
        return None
    s = (a - 1) / 2
    P = _l10_P(a, Z, D, A, delta, s)
    b0 = F(eps * f * q)
    S = sorted({2, 3, 5, 7} | _l10_supp(alpha) | _l10_supp(delta) | {f})
    # guards (mirror _l10_class_cert where load-bearing for the class claims)
    if not _is_prime(q) or not _is_prime(f):
        return None                    # wild Legendre / frozen escape prime
    if q in S or vp(a, 2) != 0:
        return None
    if vp(z, f) <= 0:
        return None                    # the escape hypothesis: f | numerator(z)
    if any(vp(t, f) != 0 for t in (delta, A, a, D)):
        return None                    # f may share ONLY z with the data
    if any(vp(t, q) != 0 for t in (delta, A, a)):
        return None                    # wild class transport needs 5 | N, q | 2A
    # NOTE: q IS allowed to divide z or D_z (5 rows have 41 | D_z); the base
    # wild symbol is then still verified directly below, and no class member
    # Q = q + kN (k >= 1) can equal q since q | N is ruled out next.
    ks = {p: _l10_exponent(P, b0, p) for p in S}
    N = 8
    for p, k in ks.items():
        N = N * p ** k // math.gcd(N, p ** k)
    m = 4 * A.numerator * A.denominator
    N = N * m // math.gcd(N, m)
    if math.gcd(q, N) != 1:
        return None                    # the progression holds no other prime
    c0 = _sun_h(a, b0, Z)
    if c0 is None:
        return None
    M0 = 16 - delta * c0 * c0 - 32 * A * b0 * s * s
    if M0 == 0:
        return None
    Pb0 = _eval_P(P, b0)
    # kernel identity x = alpha*P(b)/(D^2 A^2 b^4), i.e. M0 = P(b)/(D^2 A^2 b^4):
    assert M0 == Pb0 / (b0 ** 4 * D * D * A * A), (w, u1, u2)
    x0, d0 = alpha * M0, alpha * 2 * b0
    syms = {p: hilbert(x0, d0, p) for p in S}
    nz = [i for i, c in enumerate(P) if c != 0]
    lead = P[nz[-1]]
    Q0 = 1 + max((abs(c / lead) for c in P[:nz[-1]]), default=F(0))
    sx = (1 if alpha > 0 else -1) * (1 if lead > 0 else -1) * (eps ** nz[-1])
    sd = 1 if alpha * eps > 0 else -1
    syms[OO] = -1 if (sx < 0 and sd < 0) else 1     # asymptotic past Q0
    syms['q1'] = legendre(A, q)                     # class value of (A|Q)
    excluded = sorted(_l10_supp(D) | _l10_supp(z.numerator)
                      | _l10_supp(delta) | _l10_supp(a))
    return {'ok': all(v == 1 for v in syms.values()), 'N': N, 'S': S,
            'ks': ks, 'syms': syms, 'Q0': Q0, 'f': f, 'eps': eps, 'q': q,
            'excluded': excluded}


# ---------------------------------------------------------------- sampling

_SPOT_PRIMES = [p for p in primerange(11, 10000)]


def _member_setup(P, D, f, eps, Q):
    """x0, d0 and P(b) at the member b = eps*f*Q (exact fractions)."""
    b = F(eps * f * Q)
    Pb = _eval_P(P, b)
    x0 = _ALPHA * Pb / (b ** 4 * D * D * _A * _A)
    d0 = _ALPHA * 2 * b
    return b, Pb, x0, d0


def _check_member(rec, cert, P, D, Pb0_vals, k, Q):
    """STEP (i) verification at one member b = eps*f*Q: frozen constancy,
    wild and real place.  Appends to rec['failures'] on any violation."""
    f, eps, q, N, S, Q0 = (cert['f'], cert['eps'], cert['q'], cert['N'],
                           cert['S'], cert['Q0'])
    b, Pb, x0, d0 = _member_setup(P, D, f, eps, Q)
    fail = rec['failures']

    def bad(desc):
        fail.append([k, str(Q), desc])

    # member eligibility (the class claim is void on the excluded set)
    if Q in cert['excluded']:
        rec['excluded_members'].append([k, Q])
        return False
    if math.gcd(Q, N) != 1:
        bad('gcd(Q, N) != 1')
        return False
    if Pb == 0:
        bad('P(b) = 0: x0 degenerate (impossible past Q0)')
        return False
    # modulus hypothesis of the exponent lemma, per frozen prime:
    b0 = F(eps * f * q)
    for p, kp in cert['ks'].items():
        if vp(b - b0, p) < kp:
            bad(f'modulus hypothesis fails at {p}: b !≡ b0 mod {p}^{kp}')
        v = vp(Pb, p)
        if v != Pb0_vals[p]:
            bad(f'class constancy fails at {p}: v_p(P(b))={v} != {Pb0_vals[p]}')
    # wild place Q: valuations, L9a, and the class symbol
    if vp(d0, Q) != 1 or vp(Pb, Q) != 0 or vp(D, Q) != 0 or vp(Pb, 5) != Pb0_vals[5]:
        bad('wild/side valuations off')
    if vp(x0, Q) != -4:
        bad(f'v_Q(x0) = {vp(x0, Q)} != -4')
    if hilbert(x0, d0, Q) != 1:
        bad('wild symbol (x0,d0)_Q != +1')
    if legendre(_A, Q) != 1 or legendre(_A, Q) != legendre(_A, q):
        bad('wild class symbol (A|Q) != (A|q) = +1')
    # exact support of d0 = 8*eps*f*Q is {2, f, Q}:
    sup_d0 = {p for p in (set(S) | {Q}) if vp(d0, p) != 0}
    if sup_d0 != {2, f, Q}:
        bad(f'supp(d0) = {sup_d0} != {{2, f, Q}}')
    # frozen symbols +1, equal to the recorded base symbols
    for p in S:
        h = hilbert(x0, d0, p)
        if h != 1 or h != cert['syms'][p]:
            bad(f'frozen symbol at {p}: {h}')
    # real place: the Cauchy bound |b| = f*|Q| > Q0 forces the asymptotic
    # sign; where it applies the exact symbol must equal cert's asymptotic
    # value.  The exact real symbol is checked unconditionally too.
    if F(f) * abs(Q) <= Q0:
        rec['cauchy_violations'] += 1
    elif hilbert(x0, d0, OO) != cert['syms'][OO]:
        bad(f'real symbol {hilbert(x0, d0, OO)} != asymptotic {cert["syms"][OO]}')
    if hilbert(x0, d0, OO) != 1:
        bad('real place != +1')
    rec['members'] += 1
    return True


def _probe_dichotomy(rec, cert, P, D, k, Q):
    """STEP (ii) probing at a sampled member: for every odd prime l not in
    S_min u {Q} with l < 10000 (plus every prime of the row's excluded set),
    the exact valuation dichotomy: hilbert is FORCED +1 by valuations iff
    v_l(P(b)) is even; odd v_l(P(b)) means emergent.  Then, where the
    factorization budget permits, the COMPLETE emergent set and the global
    reciprocity cross-check.  Refusals are counted, never evidence."""
    f, eps, q, S = cert['f'], cert['eps'], cert['q'], cert['S']
    Smin = set(S)
    b, Pb, x0, d0 = _member_setup(P, D, f, eps, Q)
    datum = [k, str(Q)]
    fail = rec['failures']
    # (1) exact dichotomy on the spot set (no factorization needed: the
    #     kernel's vp / hilbert / legendre are exact on any input)
    spot = [l for l in _SPOT_PRIMES if l not in Smin and l != Q]
    spot += [l for l in cert['excluded']
             if l not in Smin and l != Q and l > 7 and l not in spot]
    for l in spot:
        vpP, vx, vd = vp(Pb, l), vp(x0, l), vp(d0, l)
        if vd != 0:
            fail.append(datum + [f'v_{l}(d0) = {vd} != 0 off-controlled'])
        if vx % 2 != vpP % 2:
            fail.append(datum + [f'uniform-fact fails at {l}'])
        h = hilbert(x0, d0, l)
        if vx % 2 == 0:
            if h != 1:
                fail.append(datum + [f'forced place {l}: symbol {h}'])
        else:
            if h != legendre(unit_part(d0, l), l):
                fail.append(datum + [f'emergent formula mismatch at {l}'])
            rec['emergent'].append(datum + [l, vpP, vx, h, 'spot'])
            rec['emergent_count'] += 1
    # (2) complete enumeration of the emergent set, budget permitting
    num, den = abs(Pb.numerator), Pb.denominator
    for p in Smin | {Q}:
        while num % p == 0:
            num //= p
        while den % p == 0:
            den //= p
    try:
        fn = factorint(num)
        fd = factorint(den)
    except (FactorBudget, PrimalityBound):
        rec['refusals'] += 1
        return
    full_emergent = []
    complete = True
    for l in sorted(set(fn) | set(fd)):
        if l in Smin or l == Q:
            continue
        vpP = fn.get(l, 0) - fd.get(l, 0)
        vx = vp(x0, l)
        h = hilbert(x0, d0, l)
        if vpP % 2 == 0:
            if vx % 2 != 0 or h != 1:
                fail.append(datum + [f'forced place {l} violated (exact)'])
        else:
            if h != legendre(unit_part(d0, l), l):
                fail.append(datum + [f'emergent formula mismatch at {l} (exact)'])
            full_emergent.append((l, vpP, vx, h))
            rec['emergent'].append(datum + [l, vpP, vx, h, 'exact'])
    # reciprocity: with frozen, wild, real all +1 and non-emergent odd
    # places forced, the emergent -1 count must be even.
    if complete:
        rec['complete_facs'] += 1
        negativity = sum(1 for (_, _, _, h) in full_emergent if h == -1)
        if negativity % 2 != 0:
            fail.append(datum + [f'reciprocity: {negativity} emergent -1s'])


def verify_row(w, u1, u2, kmax=40, sample=12):
    """Run STEP (i) on every prime member Q = q + kN (k = 1..kmax) and the
    STEP (ii) probe on `sample` of them; returns the JSON-ready record."""
    key = (w, (u1, u2))
    cert = l12_class_cert(w, u1, u2)
    if cert is None:
        return {'row': str(key), 'cert': None,
                'failures': [[None, None, 'certificate refused']]}
    a = F(1)
    z = F(w) * F(u1, u2)
    Z = z ** 3
    D = 1 - Z - a * a * Z * Z
    P = _l10_P(a, Z, D, _A, _DELTA, F(0))
    b0 = F(cert['eps'] * cert['f'] * cert['q'])
    Pb0 = _eval_P(P, b0)
    Pb0_vals = {p: vp(Pb0, p) for p in cert['S']}
    rec = {'row': str(key), 'N': cert['N'], 'f': cert['f'], 'eps': cert['eps'],
           'q': cert['q'], 'S': cert['S'], 'ks': cert['ks'],
           'excluded': cert['excluded'], 'excluded_members': [],
           'ok_cert': cert['ok'], 'Q0': str(cert['Q0']),
           'members': 0, 'sampled': 0, 'failures': [], 'refusals': 0,
           'emergent': [], 'emergent_count': 0, 'complete_facs': 0,
           'cauchy_violations': 0, 'primality_refusals': 0}
    members = []
    for k in range(1, kmax + 1):
        Q = cert['q'] + k * cert['N']
        if math.gcd(Q, cert['N']) != 1:
            rec['failures'].append([k, str(Q), 'gcd(Q, N) != 1'])
            continue
        try:
            isp = _is_prime(Q)
        except (FactorBudget, PrimalityBound):
            rec['primality_refusals'] += 1
            rec['refusals'] += 1
            continue
        if not isp:
            continue
        if _check_member(rec, cert, P, D, Pb0_vals, k, Q):
            members.append((k, Q))
    for k, Q in members[:sample]:
        _probe_dichotomy(rec, cert, P, D, k, Q)
    rec['sampled'] = sum(1 for (k, _) in members[:sample]
                         if not any(fail[0] == k for fail in rec['failures']))
    return rec


def verify_all(out_path, procs=1, kmax=40, sample=12):
    """Certify all 103 rows and write the sample data file."""
    rows = sorted(_L12_ESCAPE)
    if procs > 1:
        import multiprocessing as mp
        with mp.Pool(procs) as pool:
            results = pool.starmap(verify_row,
                                   [(w, u1, u2, kmax, sample)
                                    for (w, (u1, u2)) in rows])
    else:
        results = [verify_row(w, u1, u2, kmax, sample) for (w, (u1, u2)) in rows]
    data = {}
    for (w, (u1, u2)), rec in zip(rows, results):
        data[str((w, (u1, u2)))] = rec
    with open(out_path, 'w') as fh:
        import json
        json.dump(data, fh, indent=1, sort_keys=True)
    return data


if __name__ == '__main__':
    import sys
    out = sys.argv[1] if len(sys.argv) > 1 else 'data/l12b_class_sample.json'
    procs = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    main_data = verify_all(out, procs=procs)
    n_ok = sum(1 for r in main_data.values() if r.get('ok_cert'))
    n_fail = sum(len(r['failures']) for r in main_data.values())
    print(f'certified {n_ok}/{len(main_data)} rows; failures: {n_fail}')
    sys.exit(0 if (n_ok == len(main_data) and n_fail == 0) else 1)
