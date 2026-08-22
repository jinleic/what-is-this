"""L17 matched-evaluation artifact regenerator (REPO-PRESERVED).

Rebuilds data/l17_matched.jsonl from:
  - data/l17_sieve.jsonl            (member rows, FullGridSieve)
  - data/l17_badroots_closures.jsonl (closure-authority roots, RateModel)

Rules (advisory-locked):
  * class small factor interval from bad_root_lists: known = certified root
    masses + resolved_mass prefixes (v2); residual = unresolved bounds;
    partial classes (excluded pairs) are labeled, their factor bias noted;
  * marginal expectation per hit = p * m_root (certified roots only);
  * union expectation per member = 1 - prod over hit roots (1 - p*m_root);
  * permutation null: per class-prime, observed odd-hit count held fixed,
    assignment randomized among eligible members; 2000 sims, seed 20260819.
Run: python3 l17_matched.py
"""
import json
import random
import statistics as st
from fractions import Fraction as Fr

SV = 'data/l17_sieve.jsonl'
BC = 'data/l17_badroots_closures.jsonl'
OUT = 'data/l17_matched.jsonl'
SEED = 20260819
SIMS = 2000

svcs = {(r['cell'][0], tuple(r['cell'][1])): r for r in
        (json.loads(l) for l in open(SV)) if r['type'] == 'class-summary'}
bcc = {(r['cell'][0], tuple(r['cell'][1])): r for r in
       (json.loads(l) for l in open(BC)) if r['type'] == 'class'}

rows_meta = [{
    'type': 'meta', 'layer': 'L17 matched evaluation', 'date': '2026-08-19',
    'author': 'lead (kernels: l17_sieve members + l17_badroots_closures)',
    'version': 'post-advisory-corrected-v2',
    'generator': 'math/h10q/l17_matched.py (repo-preserved)',
    'p_limit': 10000,
    'small_factor_rule': 'known_p = certified root masses + resolved_mass '
                         'prefixes (v2 artifact); window factor interval = '
                         '[prod(1-known_p-residual_p), prod(1-known_p)]; '
                         'f_p_table is NOT used (it blanks excluded entries); '
                         'p<=10000-window statistic; no p>10000 tail bound',
    'excluded_bound_semantics': 'per excluded root: known mass = resolved_mass prefix '
                         '(exact Fraction), residual = unfinished Taylor-tree remainder; '
                         'window factor interval = [prod(1-known-resid), prod(1-known)]; '
                         'max unresolved residual across the 67 pairs: 6.9e-15',
    'notes': 'Marginal-preserving permutation null; kind-aware matched '
             'expectation. Original-estimand +0.019*/3.6sigma = Haar(Z_p)-vs-'
             'k[0,119] window gap, characterized; NOT dependence.'},
    {'type': 'meta-methods',
     'matched_support': 'closure-authority progressions; identical member sample both sides',
     'roots': 'bad-signed root residues k mod p, certified + excluded-root resolved prefixes (exact Fractions), p<=10000',
     'odd_valuation_masses': 'simple roots 1/(p+1); step roots certified root-specific; excluded roots carry exact resolved_mass prefix + residual bound; conditional odd-prob per hit p*m_root',
     'union_expectation': 'sum over members of 1 - prod over hit roots (1 - p*m_root)',
     'permutation_null': f'per class-prime: observed odd-hit count fixed, assignment randomized among eligible members, {SIMS} sims, seed {SEED}',
     'original_estimand': 'class-weighted mean of obs_clean - recomputed p<=10000 class factor'}]

per_class = []
struct = []
for cell in sorted(svcs):
    sv = svcs[cell]
    bro = bcc[cell]
    mems = [m for m in sv['members'] if m['Q_prime_proved']]
    if not mems:
        continue
    roots = {}
    for e in bro['bad_root_lists']:
        p = e['p']
        if p > 10000:
            continue
        for x in e['roots']:
            if x['status'] == 'certified':
                roots.setdefault(p, []).append((x['r'], Fr(x['mass_num'], x['mass_den'])))
            elif x.get('resolved_mass_num') is not None:
                # excluded singular root: carry the exact resolved prefix;
                # residual interval tracked in the factor rows, negligible
                # for the point statistics below.
                roots.setdefault(p, []).append((x['r'], Fr(x['resolved_mass_num'],
                                                          x['resolved_mass_den'])))
    # Window factor from bad_root_lists (v2 artifact): known_p = certified
    # root masses + resolved_mass prefix on excluded (singular-branch) roots;
    # residual_p = sum of unresolved residuals. CLEAN window factor:
    # prod (1 - known_p)  is the UPPER bound; prod (1 - known_p - residual_p)
    # the LOWER bound.  All 67 excluded pairs: verified no certified mass
    # lost; resolved prefixes now persisted (RateModel v2 patch).
    # p<=10000 window statement only: no p>10000 tail bound.
    sf = sf_lo = 1.0
    for e in bro['bad_root_lists']:
        if e['p'] > 10000:
            continue
        kn = 0.0
        resid = 0.0
        for x in e['roots']:
            if x['status'] == 'certified':
                kn += x['mass_num']/x['mass_den']
            elif x.get('resolved_mass_num') is not None:
                kn += x['resolved_mass_num']/x['resolved_mass_den']
                resid += float(x.get('residual_mass_upper_bound_float', 0.0))
        sf *= 1 - kn
        sf_lo *= max(0.0, 1 - kn - resid)
    obs_hits = mexp_hits = mod_bad = obs_bad = 0
    per_p = []
    for p, rl in roots.items():
        rset = {r for r, _ in rl}
        elig = [i for i, m in enumerate(mems) if m['k'] % p in rset]
        nodd = sum(1 for m in mems if any(pp == p and sy == -1
                                          for pp, sy in m['em_small']))
        per_p.append((p, elig, nodd))
    obs_union = 0
    for m in mems:
        neg = {p for p, sy in m['em_small'] if sy == -1 and p <= 10000}
        obs_hits += len(neg); obs_bad += bool(neg); obs_union += bool(neg)
        cp = 1.0
        for p, rl in roots.items():
            for r, mass in rl:
                if m['k'] % p == r:
                    cond = mass * p
                    cp *= float(1 - cond)
                    mexp_hits += float(cond)
        mod_bad += 1 - cp
    per_class.append({'cell': list(cell), 'n_members': len(mems),
                      'observed_odd_hits': obs_hits,
                      'matched_marginal_expectation': round(mexp_hits, 3),
                      'observed_bad_members': obs_bad,
                      'union_expectation': round(mod_bad, 3),
                      'hit_conditioned_bad_residual': round((obs_bad - mod_bad)/len(mems), 5),
                      'obs_clean': round(1 - obs_bad/len(mems), 5),
                      'small_factor_recomputed_upper': round(sf, 6),
                      'small_factor_recomputed_lower': round(sf_lo, 6),
                      'small_factor_window': 10000,
                      'partial': bro['partial']})
    struct.append((cell, len(mems), per_p, obs_union, sf, sf_lo,
                   1 - obs_bad/len(mems)))

tot_members = sum(x['n_members'] for x in per_class)
tot_bad = sum(x['observed_bad_members'] for x in per_class)
tot_hits = sum(x['observed_odd_hits'] for x in per_class)
tot_mexp = sum(x['matched_marginal_expectation'] for x in per_class)
tot_union = sum(x['union_expectation'] for x in per_class)
res = [x['hit_conditioned_bad_residual'] for x in per_class]
def _z(v):
    m = st.mean(v); return round(m, 5), round(m/(st.stdev(v)/len(v)**0.5), 2)
# estimand e = obs_clean - factor is DECREASING in the factor:
# lower endpoint uses the factor's UPPER bound, upper endpoint the LOWER bound.
est_lo_mean, est_lo_z = _z([x['obs_clean'] - x['small_factor_recomputed_upper'] for x in per_class])
est_hi_mean, est_hi_z = _z([x['obs_clean'] - x['small_factor_recomputed_lower'] for x in per_class])

rng = random.Random(SEED)
sums = []
for _ in range(SIMS):
    t = 0
    for cell, n, per_p, _, _, _, _ in struct:
        bad = set()
        for p, elig, nodd in per_p:
            if nodd and elig:
                for i in rng.sample(elig, min(nodd, len(elig))):
                    bad.add(i)
        t += len(bad)
    sums.append(t)
mu, sd = st.mean(sums), st.stdev(sums)
zz = (tot_bad - mu)/sd
# Standard two-sided empirical tail p-value (integer support):
p_hi = (1 + sum(1 for s in sums if s >= tot_bad)) / (SIMS + 1)
p_lo = (1 + sum(1 for s in sums if s <= tot_bad)) / (SIMS + 1)
p_emp = min(2 * min(p_hi, p_lo), 1.0)

glob = {'type': 'global', 'members': tot_members,
        'observed_bad_members': tot_bad, 'observed_odd_hits': tot_hits,
        'matched_marginal_expectation': round(tot_mexp, 1),
        'marginal_ratio': round(tot_hits/tot_mexp, 5),
        'union_expectation_hit_conditioned': round(tot_union, 1),
        'union_residual': round(tot_bad - tot_union, 1),
        'per_class_hit_conditioned_residual_mean': round(st.mean(res), 5),
        'per_class_hit_conditioned_residual_sd': round(st.stdev(res), 4),
        'per_class_hit_conditioned_residual_max': max(res),
        'per_class_hit_conditioned_residual_min': min(res),
        'classes_beyond_0p1': sum(1 for x in res if abs(x) > 0.1),
        'original_estimand_interval': [est_lo_mean, est_hi_mean],
        'original_estimand_z_interval': [est_lo_z, est_hi_z],
        'original_estimand_note': 'truncated factor is an UPPER bound on clean mass; interval brackets the Haar-window estimand incl. excluded-pair residual bounds',
        'original_estimand_interpretation': 'Haar Z_p-uniform factor vs observed clean on k-window [0,119]; window is the matched support',
        'permutation_null_mean': round(mu, 1), 'permutation_null_sd': round(sd, 2),
        'permutation_z': round(zz, 3),
        'permutation_p_two_sided_empirical': round(p_emp, 4),
        'permutation_p_upper_tail': round(p_hi, 4),
        'permutation_p_lower_tail': round(p_lo, 4),
        'permutation_p_definition': '2*min(p_hi,p_lo) with p_hi=(1+#{T>=obs})/(SIMS+1), p_lo=(1+#{T<=obs})/(SIMS+1)',
        'permutation_sims': SIMS,
        'permutation_seed': SEED,
        'partial_classes': sum(1 for x in per_class if x['partial']),
        'partial_note': 'partial classes carry excluded (class,p) pairs with residual bounds; their truncated factors are UPPER bounds, intervalized by excluded_mass_upper_bound'}
with open(OUT, 'w') as out:
    for m in rows_meta:
        out.write(json.dumps(m) + '\n')
    out.write(json.dumps(glob) + '\n')
    for row in per_class:
        out.write(json.dumps({'type': 'class', **row}) + '\n')
print('wrote', OUT, ':', 3 + len(per_class), 'lines')
print(json.dumps({k: glob[k] for k in ('members', 'observed_bad_members',
      'marginal_ratio', 'union_residual', 'permutation_z',
      'original_estimand_interval', 'original_estimand_z_interval')}, indent=1))
