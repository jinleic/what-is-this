
import sys, os, json
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/mm3/src')
os.chdir('/Users/jinleic/jinleic-workspace/cs/mm3/src')
from gate_b_floor import prep, canon, factor_targets
from tensor_data import LEFT_SLP, RIGHT_SLP

def slp_env(slp, prefix):
    env = {}
    for name, terms in slp:
        acc = {}
        for sign, atom in terms:
            if atom.startswith(prefix):
                v=[0]*9; v[int(atom[1:])]=1; t=tuple(v)
            else: t=env[atom]
            for i,x in enumerate(t): acc[i]=acc.get(i,0)+sign*x
        env[name]=tuple(acc[i] for i in range(9))
    return env

Lenv = slp_env(LEFT_SLP, 'A'); Renv = slp_env(RIGHT_SLP, 'B')

def build_aux_cnf(classes, reps, T, aux_classes):
    """Aux-1 encoding: like floor CNF but value set includes aux classes, and classes
    with stale empty reps in `reps` get recomputed over the extended base."""
    n = len(classes)
    input_dirs = {canon(tuple(1 if j==i else 0 for j in range(9))) for i in range(9)}
    base = sorted(set(classes) | input_dirs | set(aux_classes))
    classes_ext = list(classes) + list(aux_classes)
    def reps_of(val):
        if val in input_dirs: return []
        found = []
        for a in base:
            for b in base:
                ok = False
                for sa in (1, -1):
                    for sb in (1, -1):
                        s = tuple(sa*a[i]+sb*b[i] for i in range(9))
                        if s == val or s == tuple(-x for x in val):
                            ok = True; break
                    if ok: break
                if ok: found.append((a, b))
        seen = set(); dd = []
        for pr in found:
            key = tuple(sorted(pr))
            if key not in seen: seen.add(key); dd.append(pr)
        return dd
    reps_cache = {val: reps_of(val) for val in classes_ext}
    # sanity: every target class must have reps (else the question is trivially no)
    clauses = []; ctr=[0]
    def newvar():
        ctr[0]+=1; return ctr[0]
    Y = {}
    for g in range(T):
        for i in range(len(classes_ext)):
            Y[(g,i)] = newvar()
    P = {}
    for g in range(T):
        for i, val in enumerate(classes_ext):
            for r in range(len(reps_cache[val])):
                P[(g,i,r)] = newvar()
    m = len(classes_ext)
    for g in range(T):
        for a in range(m):
            for b in range(a+1, m):
                clauses.append([-Y[(g,a)], -Y[(g,b)]])
    for c in range(n):
        clauses.append([Y[(g,c)] for g in range(T)])
    for g in range(T):
        for i, val in enumerate(classes_ext):
            rv = reps_cache[val]
            if not rv:
                if i < n:
                    clauses.append([-Y[(g,i)]])
                continue
            clauses.append([-Y[(g,i)]] + [P[(g,i,r)] for r in range(len(rv))])
            for r in range(len(rv)):
                clauses.append([-P[(g,i,r)], Y[(g,i)]])
                for opnd in rv[r]:
                    if opnd in input_dirs: continue
                    if opnd in classes_ext:
                        ii = classes_ext.index(opnd)
                        if g == 0:
                            clauses.append([-P[(g,i,r)]])
                        else:
                            clauses.append([-P[(g,i,r)]] + [Y[(h,ii)] for h in range(g)])
    return clauses, ctr[0], reps_cache

U_t, V_t, W_t = factor_targets()
from pysat.solvers import Cadical153
results = {}
for name, targets, envd, T, auxkey in (('U', U_t, Lenv, 13, 'u12'), ('V', V_t, Renv, 14, 'v9')):
    classes, reps = prep(targets)
    aux = [canon(envd[auxkey])]
    clauses, nv, rc = build_aux_cnf(classes, reps, T, aux)
    s = Cadical153()
    for cl in clauses: s.add_clause(cl)
    sat = s.solve()
    s.delete()
    print(name, 'T=', T, 'aux=', aux, '->', 'SAT' if sat else 'UNSAT', '(vars', nv, 'clauses', len(clauses), ')')
    results[name] = {'T': T, 'aux': [list(a) for a in aux], 'SAT': bool(sat), 'vars': nv, 'clauses': len(clauses)}
json.dump(results, open('/Users/jinleic/jinleic-workspace/cs/mm3/scratch/control_aux1.json','w'), indent=2)
