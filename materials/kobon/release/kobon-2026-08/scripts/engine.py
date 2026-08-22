"""Kobon triangle problem: general-configuration decision engine.

Decides, for n lines in the plane (parallels and multi-line intersection points
allowed), whether >= target pairwise interior-disjoint triangles ("Kobon
triangles", Clement-Bader convention: other lines MAY cross a selected
triangle) can be realized.

Architecture (soundness asymmetry):
  * The CNF is a RELAXATION: every clause is NECESSARY -- satisfied by the
    combinatorial data of every real line configuration after normalization
    (rotate so no vertical lines, sort indices by slope, parallel classes are
    index-contiguous).  Hence UNSAT at target T proves K(n) < T uncondition-
    ally (modulo solver trust; rerun with DRAT proof for certification).
  * A SAT model is only a CANDIDATE: it is realized by fixed-slope LP
    straightening and then verified with EXACT rational arithmetic; only the
    verified rational configuration constitutes a lower-bound proof.

Machine-verified lemmas backing the encoding (see session log / report):
  L1 (middle crossing): for slope-sorted a<b<c pairwise crossing, x_ac lies
      strictly between x_ab and x_bc.  11,810 triples, 0 failures.
  L2 (sigma 3-way): X(a;b,c) = X(b;a,c) = X(c;a,b)  (single orientation bit).
  L3 (sidedness): sigma=+  <=>  P_bc above a  <=>  P_ab above c
                            <=>  P_ac BELOW b.
  L4 (parallel linkage): u||w, line v crossing both meets the LOWER of {u,w}
      first iff slope(v) > slope(u).
  L5 (disjointness): rules R1-R4 (below) are exactly equivalent to exact
      interior-disjointness on 167,557 random triangle pairs in configurations
      with planted parallels and concurrencies (0 necessity violations,
      0 sufficiency gaps).

Validation ladder (this engine, SAT/UNSAT pairs):
  K(3)=1 K(4)=2 K(5)=5 K(6)=7 K(7)=11 K(8)=15 K(9)=21   -- all reproduced.
"""

from fractions import Fraction as Fr
from itertools import combinations, product
from math import lcm, comb
import random, time, os

from pysat.formula import CNF, IDPool
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver


# ============================= exact geometry ==============================

def gen_config(n, n_par_pairs=0, n_conc=0, coord=10**4, seed=0):
    """Random configuration y = m_i x + b_i with planted degeneracies."""
    rng = random.Random(seed)
    for _ in range(1000):
        ms = sorted(Fr(m) for m in rng.sample(range(-coord, coord), n))
        for _ in range(n_par_pairs):
            i = rng.randrange(n - 1); ms[i + 1] = ms[i]
        bs = [Fr(rng.randrange(-coord, coord)) for _ in range(n)]
        for _ in range(n_conc):
            i, j, k = rng.sample(range(n), 3)
            if ms[i] == ms[j] or ms[k] == ms[i] or ms[k] == ms[j]:
                continue
            x = (bs[j]-bs[i])/(ms[i]-ms[j]); bs[k] = ms[i]*x + bs[i] - ms[k]*x
        if any(ms[i] == ms[j] and bs[i] == bs[j]
               for i, j in combinations(range(n), 2)):
            continue
        return ms, bs
    raise RuntimeError("generation failed")


def crossings(n, ms, bs):
    Xp = {}
    for i, j in combinations(range(n), 2):
        if ms[i] == ms[j]:
            Xp[i, j] = None
        else:
            x = (bs[j]-bs[i])/(ms[i]-ms[j]); Xp[i, j] = (x, ms[i]*x + bs[i])
    return Xp


def cross(Xp, i, j): return Xp[min(i, j), max(i, j)]


def tri_ok(Xp, t):
    i, j, k = t
    A, B, C = cross(Xp, i, j), cross(Xp, i, k), cross(Xp, j, k)
    return A is not None and B is not None and C is not None and A != B


def tri_verts(Xp, t):
    i, j, k = t
    return (cross(Xp, i, j), cross(Xp, i, k), cross(Xp, j, k))


def face_multiplicity_counts(n, ms, bs, selected):
    """Return (line-through-triangle incidences, bounded arrangement faces).

    A line crosses a triangle's open interior exactly when the triangle has
    vertices on both strict sides of that line.  Vertex contacts are therefore
    handled without perturbation: a zero sign neither creates nor suppresses
    an incidence between the other two, oppositely signed vertices.
    """
    if len(ms) != n or len(bs) != n:
        raise ValueError("line count")
    ms = list(map(Fr, ms)); bs = list(map(Fr, bs))
    if len(set(zip(ms, bs))) != n:
        raise ValueError("duplicate lines")
    Xp = crossings(n, ms, bs)
    selected = tuple(tuple(t) for t in selected)
    if any(not tri_ok(Xp, t) for t in selected):
        raise ValueError("degenerate selected triangle")

    parallel_pairs = 0
    point_lines = {}
    for i, j in combinations(range(n), 2):
        p = cross(Xp, i, j)
        if p is None:
            parallel_pairs += 1
        else:
            point_lines.setdefault(p, set()).update((i, j))
    multiple_point_penalty = sum(comb(len(inc)-1, 2)
                                 for inc in point_lines.values())
    bounded_faces = (comb(n-1, 2) - parallel_pairs
                     - multiple_point_penalty)

    crossing_count = 0
    for t in selected:
        verts = tri_verts(Xp, t)
        for r in range(n):
            if r in t:
                continue
            signs = [y - ms[r]*x - bs[r] for x, y in verts]
            crossing_count += min(signs) < 0 < max(signs)
    return crossing_count, bounded_faces


def hom(p):
    x, y = p
    d = lcm(x.denominator, y.denominator)
    return (int(x*d), int(y*d), d)


def orient_h(P, Q, R):
    a1, b1, d1 = P; a2, b2, d2 = Q; a3, b3, d3 = R
    return (a2*d1-a1*d2)*(b3*d1-b1*d3) - (b2*d1-b1*d2)*(a3*d1-a1*d3)


def interiors_disjoint_h(H1, H2):
    """Separating-axis over Q: interiors of two nondegenerate closed triangles
    are disjoint iff some side-line of one weakly separates them."""
    for own, other in ((H1, H2), (H2, H1)):
        for A, B in ((0, 1), (0, 2), (1, 2)):
            so = [orient_h(own[A], own[B], p) for p in own]
            st = [orient_h(own[A], own[B], p) for p in other]
            if all(v >= 0 for v in so) and all(v <= 0 for v in st): return True
            if all(v <= 0 for v in so) and all(v >= 0 for v in st): return True
    return False


# ============================== CNF builder ================================

def build_model(n, target, rules=("R1", "R2", "R3", "R4")):
    """CNF relaxation: 'exists n-line configuration with >= target pairwise
    interior-disjoint triangles'.  Variables:
      P(i,j)   lines i,j cross (else parallel; classes index-contiguous)
      C(i,j,k) concurrent triple (common point)
      X(r,i,j) crossing of i strictly precedes crossing of j along line r
      U(i,j)   for parallel candidates: lower-index line lies above
      S(t)     triple t selected as a Kobon triangle
      SA/SB(p,q,r) vertex P_pq strictly above/below line r  (aux, one-sided)
      B(r,u,a,b)   point of u on r strictly between points of a and b (aux)
      IN(p,q,T)    vertex P_pq strictly inside triangle T (aux)
    """
    pool = IDPool(); cnf = CNF()
    sp = lambda *a: tuple(sorted(a))
    P  = lambda i, j: pool.id(("P",)+sp(i, j))
    Cv = lambda i, j, k: pool.id(("C",)+sp(i, j, k))
    X  = lambda r, i, j: pool.id(("X", r, i, j))
    U  = lambda i, j: pool.id(("U",)+sp(i, j))
    S  = lambda t: pool.id(("S",)+sp(*t))
    lines = range(n)
    trip = list(combinations(lines, 3))
    prs = list(combinations(lines, 2))

    # ---- A1: parallel classes contiguous + transitive;  A2: C -> crossings
    for i, j, k in trip:
        cnf.append([P(i, k), -P(i, j)]); cnf.append([P(i, k), -P(j, k)])
        cnf.append([P(i, j), P(j, k), -P(i, k)])
        for a, b in combinations((i, j, k), 2):
            cnf.append([-Cv(i, j, k), P(a, b)])
    # ---- A3: concurrency closure on 4-sets
    for q in combinations(lines, 4):
        ts = list(combinations(q, 3))
        for t1, t2 in combinations(ts, 2):
            for t3 in ts:
                if t3 not in (t1, t2):
                    cnf.append([-Cv(*t1), -Cv(*t2), Cv(*t3)])
    # ---- A4: totality with ties; A5: transitivity (incl. tie substitution)
    for r in lines:
        oth = [x for x in lines if x != r]
        for i, j in combinations(oth, 2):
            c = Cv(*sp(r, i, j))
            cnf.append([-P(*sp(r, i)), -P(*sp(r, j)), X(r, i, j), X(r, j, i), c])
            cnf.append([-X(r, i, j), -X(r, j, i)])
            cnf.append([-X(r, i, j), -c]); cnf.append([-X(r, j, i), -c])
            for w in (i, j):
                cnf.append([P(*sp(r, w)), -X(r, i, j)])
                cnf.append([P(*sp(r, w)), -X(r, j, i)])
        for i, j, k in combinations(oth, 3):
            for a, b, c2 in ((i,j,k),(i,k,j),(j,i,k),(j,k,i),(k,i,j),(k,j,i)):
                cnf.append([-X(r, a, b), -X(r, b, c2), X(r, a, c2)])
        for i, j in combinations(oth, 2):
            for k in oth:
                if k in (i, j): continue
                cnf.append([-Cv(*sp(r, i, j)), -X(r, j, k), X(r, i, k)])
                cnf.append([-Cv(*sp(r, i, j)), -X(r, k, j), X(r, k, i)])
    # ---- A6: sigma three-way equality (L2), guarded
    for a, b, c in trip:
        g = [-P(a, b), -P(a, c), -P(b, c), Cv(a, b, c)]
        for L1, L2 in ((X(a, b, c), X(b, a, c)), (X(b, a, c), X(c, a, b))):
            cnf.append(g+[-L1, L2]); cnf.append(g+[L1, -L2])
    # ---- A7: parallel-order linkage (L4) + U transitivity within classes
    for i, j in prs:
        for v in lines:
            if v in (i, j): continue
            pv = P(*sp(v, i))
            if v > j:
                cnf.append([P(i, j), -pv, -U(i, j), X(v, j, i)])
                cnf.append([P(i, j), -pv, U(i, j), X(v, i, j)])
            elif v < i:
                cnf.append([P(i, j), -pv, -U(i, j), X(v, i, j)])
                cnf.append([P(i, j), -pv, U(i, j), X(v, j, i)])
    for i, j, k in trip:
        g = [P(i, j), P(j, k)]
        cnf.append(g+[-U(i, j), -U(j, k), U(i, k)])
        cnf.append(g+[U(i, j), U(j, k), -U(i, k)])

    # ---- sidedness aux (L3 table; ties -> both literals false)
    def sigma_lits(p, q, r):
        a, b, c = sp(p, q, r)
        if r == a: return X(a, b, c), X(a, c, b)
        if r == b: return X(b, c, a), X(b, a, c)
        return X(c, a, b), X(c, b, a)
    SAv = lambda p, q, r: pool.id(("SA",)+sp(p, q)+(r,))
    SBv = lambda p, q, r: pool.id(("SB",)+sp(p, q)+(r,))
    side_done = set()
    def ensure_side(p, q, r):
        assert r not in (p, q) and p != q, ("ensure_side misuse", p, q, r)
        key = sp(p, q)+(r,)
        sa, sb = SAv(p, q, r), SBv(p, q, r)
        if key in side_done: return sa, sb
        side_done.add(key)
        la, lb = sigma_lits(p, q, r)
        g = [-P(*sp(p, q)), -P(*sp(p, r)), -P(*sp(q, r)), Cv(*sp(p, q, r))]
        cnf.append(g+[-la, sa]); cnf.append(g+[-lb, sb])
        for w in (p, q):        # w || r: P_pq on w; above r iff line w above r
            up = U(*sp(w, r))
            above = up if w < r else -up
            cnf.append([P(*sp(w, r)), -P(*sp(p, q)), -above, sa])
            cnf.append([P(*sp(w, r)), -P(*sp(p, q)), above, sb])
        # Reverse direction: SA/SB are exact whenever P_pq exists.  Exact
        # auxiliaries let disjointness clauses propagate back into the order
        # variables instead of leaving arbitrary one-sided truth values.
        cnf.append([-sa, P(*sp(p, q))]); cnf.append([-sb, P(*sp(p, q))])
        cnf.append([-sa, -sb])
        cnf.append([-P(*sp(p, q)), Cv(*sp(p, q, r)), sa, sb])
        cnf.append([-sa, -Cv(*sp(p, q, r))])
        cnf.append([-sb, -Cv(*sp(p, q, r))])
        cnf.append([-sa, -P(*sp(p, r)), -P(*sp(q, r)), la])
        cnf.append([-sb, -P(*sp(p, r)), -P(*sp(q, r)), lb])
        for w in (p, q):
            up = U(*sp(w, r))
            above = up if w < r else -up
            cnf.append([-sa, P(*sp(w, r)), above])
            cnf.append([-sb, P(*sp(w, r)), -above])
        return sa, sb

    BETv = lambda r, u, a, b: pool.id(("B", r, u)+sp(a, b))
    bet_done = set()
    def ensure_bet(r, u, a, b):
        assert u not in (a, b) and r not in (a, b) and r != u and a != b, \
            ("ensure_bet misuse", r, u, a, b)
        key = (r, u)+sp(a, b)
        bv = BETv(r, u, a, b)
        if key in bet_done: return bv
        bet_done.add(key)
        cnf.append([-X(r, a, u), -X(r, u, b), bv])
        cnf.append([-X(r, b, u), -X(r, u, a), bv])
        # B is exactly the disjunction a<u<b OR b<u<a.  Antisymmetry makes
        # the two compact reverse clauses exclude the crossed combinations.
        cnf.append([-bv, X(r, a, u), X(r, b, u)])
        cnf.append([-bv, X(r, u, b), X(r, u, a)])
        return bv

    # ---- selection base: nondegenerate triangle
    for t in trip:
        for a, b in combinations(t, 2):
            cnf.append([-S(t), P(a, b)])
        cnf.append([-S(t), -Cv(*t)])

    tri_with = {r: [t for t in trip if r in t] for r in lines}
    others = lambda t, r: tuple(x for x in t if x != r)

    # ---- FACE: selected triples are triangular arrangement regions
    if "FACE" in rules:
        for t in trip:
            for r in lines:
                if r in t:
                    continue
                above, below = [], []
                for p, q in combinations(t, 2):
                    sa, sb = ensure_side(p, q, r)
                    above.append(sa)
                    below.append(sb)
                for sa in above:
                    for sb in below:
                        cnf.append([-S(t), -sa, -sb])
            for side_line in t:
                a, b = others(t, side_line)
                for r in lines:
                    if r not in t:
                        between = ensure_bet(side_line, r, a, b)
                        cnf.append([-S(t), -between])

    # ---- R1: shared line, same side -> side-interval interiors disjoint
    if "R1" in rules:
        for r in lines:
            for T1, T2 in combinations(tri_with[r], 2):
                a, b = others(T1, r); c, d = others(T2, r)
                # T1 != T2 both contain r, so {a,b} != {c,d} always.  When the
                # triangles share a second line (say a == c), the endpoint-
                # coincidence disjunct below degenerates to the single literal
                # C(r,b,d) and no clause forbids equal base intervals on r; that
                # case is instead excluded by the R1 instance on the shared
                # line a, whose apexes P_rb = P_rd coincide (C(r,b,d)) and whose
                # base intervals properly overlap, firing a B literal there.
                sa1, sb1 = ensure_side(a, b, r)
                sa2, sb2 = ensure_side(c, d, r)
                dis = []
                for u, (e, f) in ((c, (a, b)), (d, (a, b)),
                                  (a, (c, d)), (b, (c, d))):
                    if u not in (e, f):
                        dis.append([ensure_bet(r, u, e, f)])
                if a != c and b != d:
                    dis.append([Cv(*sp(r, a, c)), Cv(*sp(r, b, d))])
                if a != d and b != c:
                    dis.append([Cv(*sp(r, a, d)), Cv(*sp(r, b, c))])
                for ss1, ss2 in ((sa1, sa2), (sb1, sb2)):
                    for D in dis:
                        cnf.append([-S(T1), -S(T2), -ss1, -ss2]
                                   + [-x for x in D])
    # ---- R2: no vertex strictly inside the other triangle
    if "R2" in rules:
        INdone = {}
        def ensure_in(p, q, T):
            key = sp(p, q)+sp(*T)
            if key in INdone: return INdone[key]
            iv = pool.id(("IN",)+key); INdone[key] = iv
            sides = [(ensure_side(p, q, r), ensure_side(*others(T, r), r))
                     for r in T]
            for pat in product((0, 1), repeat=3):
                cl = [iv]
                for (svv, sov), bit in zip(sides, pat):
                    cl.append(-svv[bit]); cl.append(-sov[bit])
                cnf.append(cl)
            return iv
        for T1, T2 in combinations(trip, 2):
            for A, B in ((T1, T2), (T2, T1)):
                for p, q in combinations(B, 2):
                    if p in A or q in A: continue   # vertex on A's line: never strictly inside
                    cnf.append([-S(A), -S(B), -ensure_in(p, q, A)])
    # ---- R3: no proper side x side transversal crossing
    if "R3" in rules:
        for T1, T2 in combinations(trip, 2):
            for r1 in T1:
                for r2 in T2:
                    if r1 == r2: continue
                    a, b = others(T1, r1); c, d = others(T2, r2)
                    if r2 in (a, b) or r1 in (c, d): continue  # Z is a vertex: not interior
                    b1 = ensure_bet(r1, r2, a, b)
                    b2 = ensure_bet(r2, r1, c, d)
                    cnf.append([-S(T1), -S(T2), -b1, -b2])
    # ---- R4: vertex on another's side (via concurrency) -> rest outside
    if "R4" in rules:
        for T1, T2 in combinations(trip, 2):
            for A, B in ((T1, T2), (T2, T1)):
                for r in A:
                    if r in B: continue
                    a, b = others(A, r)
                    saA, sbA = ensure_side(a, b, r)
                    for u, v in combinations(B, 2):
                        w = [x for x in B if x not in (u, v)][0]
                        # measure Z along r via a line outside {a,b}
                        um = u if u not in (a, b) else (v if v not in (a, b) else None)
                        if um is None: continue      # Z = shared vertex P_ab
                        cc = Cv(*sp(u, v, r))
                        bv = ensure_bet(r, um, a, b)
                        for p, q in ((u, w), (v, w)):
                            sav, sbv = ensure_side(p, q, r)
                            cnf.append([-S(A), -S(B), -cc, -bv, -sav, -saA])
                            cnf.append([-S(A), -S(B), -cc, -bv, -sbv, -sbA])

    cnf.extend(CardEnc.atleast(lits=[S(t) for t in trip], bound=target,
                               vpool=pool, encoding=EncType.seqcounter))
    return cnf, pool


def add_symbreak(cnf, pool):
    """WLOG sigma-flip (exactly ONE sigma bit may ever be fixed).

    Justification: the pi-rotation (x,y) -> (-x,-y), i.e. b_i -> -b_i, maps a
    configuration y = m_i x + b_i to y = m_i x - b_i.  Slopes, index order,
    parallel classes, concurrencies and interior-disjointness of selected
    triangles are all preserved, while every crossing x-coordinate negates, so
    every order literal X(r,i,j) swaps with X(r,j,i): ALL sigma bits flip
    simultaneously.  Hence WLOG sigma(0,1,2) = + whenever lines 0,1,2 cross
    pairwise and non-concurrently.  (The plain x-flip (x,y) -> (-x,y) does NOT
    justify this clause: it relabels i -> n-1-i and maps sigma(0,1,2) to
    sigma(n-3,n-2,n-1), acting as the identity on that bit for n = 3.  And
    because b -> -b flips all bits at once, fixing any SECOND sigma bit is not
    WLOG: for m = (-19,-3,10,18), b = (4,7,5,16) no image under the involution
    satisfies sigma(0,1,2) = sigma(0,1,3) = +.)"""
    sp = lambda *a: tuple(sorted(a))
    Pid = lambda i, j: pool.id(("P",)+sp(i, j))
    cnf.append([-Pid(0, 1), -Pid(0, 2), -Pid(1, 2),
                pool.id(("C", 0, 1, 2)), pool.id(("X", 0, 1, 2))])


def add_exact_selection(cnf, pool, n, target):
    """Require exactly ``target`` selected triangles.

    ``build_model`` already requires at least ``target``.  At most is WLOG:
    any larger pairwise interior-disjoint family contains a target-sized
    subfamily, and deselecting triangles only disables constraints.
    """
    selected = [pool.id(("S",)+t) for t in combinations(range(n), 3)]
    cnf.extend(CardEnc.atmost(lits=selected, bound=target, vpool=pool,
                              encoding=EncType.totalizer))


def add_triangle_crossing_indicators(cnf, pool, n):
    """Reify selected-triangle/open-interior crossing incidences.

    For a selected triangle T and line r outside T, TC(T,r) is true exactly
    when T has a vertex strictly above r and a vertex strictly below r.
    ``build_model`` has already made every required SA/SB literal exact.
    """
    crossing_lits = []
    for t in combinations(range(n), 3):
        selected = pool.obj2id.get(("S",)+t)
        if selected is None:
            raise ValueError("missing selected-triangle variable")
        for r in range(n):
            if r in t:
                continue
            above, below = [], []
            for p, q in combinations(t, 2):
                key = tuple(sorted((p, q)))+(r,)
                sa = pool.obj2id.get(("SA",)+key)
                sb = pool.obj2id.get(("SB",)+key)
                if sa is None or sb is None:
                    raise ValueError("missing exact sidedness variable")
                above.append(sa); below.append(sb)
            crossing = pool.id(("TC",)+t+(r,))
            crossing_lits.append(crossing)
            cnf.append([-crossing, selected])
            cnf.append([-crossing]+above)
            cnf.append([-crossing]+below)
            for sa in above:
                for sb in below:
                    cnf.append([-selected, -sa, -sb, crossing])
    return crossing_lits


def add_crossed_triangle_indicators(cnf, pool, n):
    """Flag every selected triangle crossed by at least one outside line.

    TX(T) is forced true by any strict straddling witness.  The converse is
    deliberately unnecessary: the face bound only needs a lower bound, and
    every geometric arrangement extends the CNF by setting TX exactly on its
    crossed selected triangles.
    """
    crossed_lits = []
    for t in combinations(range(n), 3):
        selected = pool.obj2id.get(("S",)+t)
        if selected is None:
            raise ValueError("missing selected-triangle variable")
        crossed = pool.id(("TX",)+t)
        crossed_lits.append(crossed)
        cnf.append([-crossed, selected])
        for r in range(n):
            if r in t:
                continue
            above, below = [], []
            for p, q in combinations(t, 2):
                key = tuple(sorted((p, q)))+(r,)
                sa = pool.obj2id.get(("SA",)+key)
                sb = pool.obj2id.get(("SB",)+key)
                if sa is None or sb is None:
                    raise ValueError("missing exact sidedness variable")
                above.append(sa)
                below.append(sb)
            for sa in above:
                for sb in below:
                    cnf.append([-selected, -sa, -sb, crossed])

    return crossed_lits


def add_no_concurrency_crossed_triangle_indicators(cnf, pool, n):
    """Compactly flag crossed selected triangles in a no-concurrency cube.

    Without triple concurrency, every vertex of a selected triangle lies
    strictly above or below an outside line.  If TX(T) is false, the three
    above-side literals must therefore agree for every outside line.  A mixed
    side pattern forces TX(T).  Callers must enforce no concurrency.
    """
    crossed_lits = []
    for t in combinations(range(n), 3):
        selected = pool.obj2id.get(("S",)+t)
        if selected is None:
            raise ValueError("missing selected-triangle variable")
        crossed = pool.id(("TX",)+t)
        crossed_lits.append(crossed)
        cnf.append([-crossed, selected])
        for r in range(n):
            if r in t:
                continue
            above = []
            for p, q in combinations(t, 2):
                key = tuple(sorted((p, q)))+(r,)
                side = pool.obj2id.get(("SA",)+key)
                if side is None:
                    raise ValueError("missing exact sidedness variable")
                above.append(side)
            pivot = above[0]
            for side in above[1:]:
                cnf.append([-selected, crossed, -pivot, side])
                cnf.append([-selected, crossed, pivot, -side])
    return crossed_lits


def add_face_bound(cnf, pool, n, target, crossing_lits=(), per_line=False):
    """Add the exact bounded-face penalty.

    For an essential arrangement of distinct affine lines,

      bounded_faces = C(n-1, 2) - Q - sum_p C(k_p-1, 2),

    where Q counts parallel pairs and k_p is the number of lines through a
    finite point p.  Each selected triangle consumes one bounded arrangement
    face.  Each supplied crossing literal certifies one further consumed face;
    callers may use either crossed-triangle flags or the stronger exact list
    of line-through-triangle incidences.  Interior-disjoint triangles consume
    disjoint face sets.  A positive target makes the arrangement essential.

    With ``per_line=True``, also add the universally necessary linewise
    capacity inequalities (the crossing-refined capability theorem, see
    report.md):

      selected sides on r (+ exact TC crossing incidences, when granted)
      <= n-2-q_r-sum_{p on r}(k_p-4).

    A k-line point is represented once: (a,b,c) are its three least-labelled
    lines, and u=k-3 later lines also pass through P_ab.  Its face penalty is
    encoded positively as 1 + 2u + C(u,2) = C(k-1,2).
    """
    if target < 1:
        raise ValueError("face bound requires a positive target")
    slack = comb(n-1, 2) - target
    if slack < 0:
        cnf.append([])
        return

    sp = lambda *a: tuple(sorted(a))
    Pid = lambda i, j: pool.id(("P",)+sp(i, j))
    Cid = lambda *a: pool.id(("C",)+sp(*a))
    penalty = list(crossing_lits)
    penalty.extend(-Pid(i, j) for i, j in combinations(range(n), 2))

    for a, b, c in combinations(range(n), 3):
        concurrent = Cid(a, b, c)
        earlier = [Cid(a, b, d) for d in range(c) if d not in (a, b)]
        first = pool.id(("ZP", a, b, c))
        cnf.append([-first, concurrent])
        for old in earlier:
            cnf.append([-first, -old])
        cnf.append([-concurrent]+earlier+[first])
        penalty.append(first)

        later = []
        for d in range(c+1, n):
            through = Cid(a, b, d)
            extra = pool.id(("ZU", a, b, c, d))
            cnf.append([-extra, first])
            cnf.append([-extra, through])
            cnf.append([-first, -through, extra])
            penalty.extend((extra, extra))
            later.append(extra)
        for d, e in combinations(later, 2):
            pair = pool.id(("ZV", a, b, c, d, e))
            cnf.append([-pair, d])
            cnf.append([-pair, e])
            cnf.append([-d, -e, pair])
            penalty.append(pair)

    cnf.extend(CardEnc.atmost(lits=penalty, bound=slack, vpool=pool,
                              encoding=EncType.seqcounter))
    if per_line:
        triples = list(combinations(range(n), 3))
        tc_lits = []
        if crossing_lits:
            tc_lits = [
                pool.id(("TC",)+t+(r,))
                for t in triples for r in range(n) if r not in t
            ]
            if set(crossing_lits) != set(tc_lits):
                raise ValueError(
                    "per_line requires the exact TC family from "
                    "add_triangle_crossing_indicators")
        base_count = comb(n-1, 2)
        for r in range(n):
            line_lits = [
                pool.id(("S",)+t) for t in triples if r in t
            ]
            if crossing_lits:
                # Couple open-interior crossing incidences into the
                # linewise slot injection: sigma_r + gamma_r <=
                # n-2-q_r-sum_{p on r}(k_p-4)  (the global sum of the
                # per-line shares is sum_p k_p(k_p-4)).
                line_lits.extend(
                    pool.id(("TC",)+t+(r,))
                    for t in triples if r not in t
                )
            line_lits.extend(
                -Pid(r, j) for j in range(n) if j != r
            )
            for a, b, c in triples:
                first = pool.id(("ZP", a, b, c))
                later = [
                    pool.id(("ZU", a, b, c, d))
                    for d in range(c+1, n)
                ]
                if r in (a, b, c):
                    line_lits.append(-first)
                    line_lits.extend(later)
                elif r > c:
                    incident = pool.id(("ZU", a, b, c, r))
                    for d, e in combinations(later, 2):
                        if incident in (d, e):
                            line_lits.append(
                                pool.id(("ZV", a, b, c, d, e)))
            bound = n-2 + base_count
            if bound < len(line_lits):
                cnf.extend(CardEnc.atmost(
                    lits=line_lits, bound=bound, vpool=pool,
                    encoding=EncType.seqcounter))


def add_no_concurrency_capacities(
        cnf, pool, n, parallel_lower_bounds, target=None):
    """Bound selected-side counts in a no-concurrency cube.

    If q_r lines are parallel to line r and all finite crossings on r are
    distinct, its selected triangle sides number at most n-2-q_r.  The supplied
    values may be lower bounds on q_r; weaker resulting capacities remain
    necessary.

    If ``target`` is supplied, the sum of the line counts is at least
    3*target, so the other lines' capacities also imply a lower bound on each
    individual line count.
    """
    if len(parallel_lower_bounds) != n:
        raise ValueError("one parallel lower bound is required per line")
    triples = list(combinations(range(n), 3))
    capacities = [n-2-q_min for q_min in parallel_lower_bounds]
    capacity_sum = sum(capacities)
    for r, cap in enumerate(capacities):
        if cap < 0:
            cnf.append([])
            continue
        selected = [pool.id(("S",)+t) for t in triples if r in t]
        cnf.extend(CardEnc.atmost(lits=selected, bound=cap, vpool=pool,
                                  encoding=EncType.seqcounter))
        if target is not None:
            lower = max(0, 3*target - (capacity_sum-cap))
            if lower > cap:
                cnf.append([])
            elif lower:
                cnf.extend(CardEnc.atleast(
                    lits=selected, bound=lower, vpool=pool,
                    encoding=EncType.seqcounter))


def add_no_concurrency_crossing_budget(
        cnf, pool, n, target, parallel_lower_bounds, crossing_lits=None,
        per_line=False):
    """Bound line-through-selected-triangle incidences without concurrency.

    Let q_r be the number of arrangement lines parallel to line r.  Selected
    sides inject into used elementary edges, while each crossing incidence
    contains a distinct unused elementary edge.  Therefore

        3*target + crossings <= sum_r (n-2-q_r).

    With ``per_line=True``, also encode the stronger local inequalities

        selected sides on r + crossings by r <= n-2-q_r.

    Lower bounds on q_r yield weaker necessary constraints.  This theorem is
    false when three or more arrangement lines may be concurrent, so callers
    must enforce no concurrency independently.
    """
    if len(parallel_lower_bounds) != n:
        raise ValueError("one parallel lower bound is required per line")
    slack = sum(n-2-q_min for q_min in parallel_lower_bounds) - 3*target
    if slack < 0:
        cnf.append([])
        return []
    supplied_crossings = crossing_lits is not None
    if supplied_crossings:
        crossing_lits = list(crossing_lits)
    else:
        crossing_lits = add_triangle_crossing_indicators(cnf, pool, n)
    if slack < len(crossing_lits):
        cnf.extend(CardEnc.atmost(
            lits=crossing_lits, bound=slack, vpool=pool,
            encoding=EncType.seqcounter))
    if per_line:
        if supplied_crossings:
            raise ValueError(
                "per-line budgets require exact generated crossing literals")
        triples = list(combinations(range(n), 3))
        for r, q_min in enumerate(parallel_lower_bounds):
            capacity = n-2-q_min
            if capacity < 0:
                cnf.append([])
                continue
            line_lits = [
                pool.id(("S",)+t) for t in triples if r in t
            ] + [
                pool.id(("TC",)+t+(r,)) for t in triples if r not in t
            ]
            if capacity < len(line_lits):
                cnf.extend(CardEnc.atmost(
                    lits=line_lits, bound=capacity, vpool=pool,
                    encoding=EncType.seqcounter))
    return crossing_lits


# ============================ decode & realize =============================

def extract(model, pool, n):
    tv = set(l for l in model if l > 0)
    sp = lambda *a: tuple(sorted(a))
    has = lambda key: pool.obj2id.get(key) in tv
    par = [(i, j) for i, j in combinations(range(n), 2)
           if not has(("P",)+sp(i, j))]
    conc = [t for t in combinations(range(n), 3) if has(("C",)+t)]
    sel = [t for t in combinations(range(n), 3) if has(("S",)+t)]
    orders = {}
    from functools import cmp_to_key
    for r in range(n):
        oth = [x for x in range(n) if x != r and (min(r, x), max(r, x)) not in
               set(par)]
        def cmpf(i, j):
            if has(("X", r, i, j)): return -1
            if has(("X", r, j, i)): return 1
            return 0
        orders[r] = sorted(oth, key=cmp_to_key(cmpf))
    return dict(par=par, conc=conc, sel=sel, orders=orders)


def straighten(n, info, tries=60, seed=0):
    """Fixed-slope LP realization + exact rationalization (see session log)."""
    import numpy as np
    from scipy.optimize import linprog
    rng = random.Random(seed)
    par = set(map(tuple, info["par"]))
    cls = list(range(n))
    def find(x):
        while cls[x] != x:
            cls[x] = cls[cls[x]]; x = cls[x]
        return x
    for i, j in par: cls[find(i)] = find(j)
    classes = {}
    for i in range(n): classes.setdefault(find(i), []).append(i)
    class_list = sorted(classes.values(), key=lambda c: min(c))
    conc = [tuple(t) for t in info["conc"]]
    conc_set = set(conc)

    for attempt in range(tries):
        deltas = [rng.randrange(3, 40) for _ in class_list]
        acc = -sum(deltas)//2
        mval = {}
        for cl, d in zip(class_list, deltas):
            acc += d
            for i in cl: mval[i] = Fr(acc)
        A_ub, b_ub, A_eq = [], [], []
        def xcoef(r, i):
            d = mval[r] - mval[i]
            row = [Fr(0)]*n; row[i] = Fr(1)/d; row[r] = Fr(-1)/d
            return row
        sub = lambda r1, r2: [a-b for a, b in zip(r1, r2)]
        for (a, b, c) in conc:
            A_eq.append(sub(xcoef(a, b), xcoef(a, c)))
        tied = lambda r, u, v: tuple(sorted((r, u, v))) in conc_set
        for r in range(n):
            sq = info["orders"][r]
            for u, v in zip(sq, sq[1:]):
                if tied(r, u, v):
                    A_eq.append(sub(xcoef(r, u), xcoef(r, v)))
                else:
                    A_ub.append(sub(xcoef(r, u), xcoef(r, v))); b_ub.append(-1)
        res = linprog([0]*n, A_ub=[[float(x) for x in r_] for r_ in A_ub],
                      b_ub=b_ub,
                      A_eq=([[float(x) for x in r_] for r_ in A_eq] or None),
                      b_eq=([0]*len(A_eq) if A_eq else None),
                      bounds=[(-10**6, 10**6)]*n, method="highs")
        if not res.success: continue
        bf = res.x
        rows = [row[:] for row in A_eq]
        piv_cols, rr = [], 0
        for col in range(n):
            p = next((k for k in range(rr, len(rows)) if rows[k][col] != 0),
                     None)
            if p is None: continue
            rows[rr], rows[p] = rows[p], rows[rr]
            pv = rows[rr][col]
            rows[rr] = [x/pv for x in rows[rr]]
            for k in range(len(rows)):
                if k != rr and rows[k][col] != 0:
                    f = rows[k][col]
                    rows[k] = [x - f*y for x, y in zip(rows[k], rows[rr])]
            piv_cols.append(col); rr += 1
            if rr == len(rows): break
        bq = [None]*n
        for c2 in range(n):
            if c2 not in piv_cols:
                bq[c2] = Fr(round(bf[c2]*1000), 1000)
        for k in range(len(piv_cols)-1, -1, -1):
            col = piv_cols[k]; row = rows[k]
            val = Fr(0)
            for c2 in range(col+1, n):
                if row[c2] != 0: val -= row[c2]*bq[c2]
            bq[col] = val
        ms = [mval[i] for i in range(n)]
        ok = True
        for row, rhs in zip(A_ub, b_ub):
            if not sum(r_*b_ for r_, b_ in zip(row, bq)) < 0:
                ok = False; break
        if ok: return ms, bq
    return None


def verify_selection(n, ms, bs, sel, minimum=0):
    """EXACT certificate check for a lower bound.

    Validates the line and selection containers before checking every selected
    triangle and every pair with homogeneous rational geometry.
    """
    if len(ms) != n or len(bs) != n:
        return False, "line count"
    lines = tuple(zip(map(Fr, ms), map(Fr, bs)))
    if len(set(lines)) != n:
        return False, "duplicate lines"
    selected = tuple(tuple(t) for t in sel)
    if len(selected) < minimum:
        return False, f"only {len(selected)} triangles"
    if len(set(selected)) != len(selected):
        return False, "duplicate triangle"
    if any(len(t) != 3 or tuple(sorted(t)) != t or
           any(not isinstance(i, int) or i < 0 or i >= n for i in t)
           for t in selected):
        return False, "invalid triangle index"
    Xp = crossings(n, list(map(Fr, ms)), list(map(Fr, bs)))
    for t in selected:
        if not tri_ok(Xp, t):
            return False, f"degenerate {t}"
    H = {t: tuple(hom(v) for v in tri_verts(Xp, t)) for t in selected}
    for t1, t2 in combinations(selected, 2):
        if not interiors_disjoint_h(H[t1], H[t2]):
            return False, f"overlap {t1} {t2}"
    return True, "ok"


N10_LOWER_BOUND_LINES = (
    (-39, 0), (-39, -9653), (-37, -7868), (-10, -589), (31, -2168),
    (33, -2158), (51, -2878), (53, 92), (58, 0), (100, 490),
)
N10_LOWER_BOUND_TRIANGLES = (
    (0, 2, 6), (0, 3, 5), (0, 3, 8), (0, 4, 5), (0, 4, 6),
    (0, 7, 8), (0, 7, 9), (1, 2, 3), (1, 4, 7), (1, 5, 7),
    (1, 5, 8), (1, 6, 8), (1, 6, 9), (2, 3, 7), (2, 4, 7),
    (2, 4, 8), (2, 5, 8), (2, 5, 9), (2, 6, 9), (3, 4, 6),
    (3, 5, 6), (3, 7, 9), (3, 8, 9), (4, 5, 9), (4, 8, 9),
)


def verify_n10_lower_bound():
    """Verify the explicit rational certificate K(10) >= 25."""
    ms, bs = zip(*N10_LOWER_BOUND_LINES)
    return verify_selection(10, list(map(Fr, ms)), list(map(Fr, bs)),
                            N10_LOWER_BOUND_TRIANGLES, minimum=25)


# K_gen(12) >= 38 certificate. Combinatorics: Kabanovitch's 12-line/38-triangle
# arrangement (Charade 6, June 1999; via Savchuk's LineOrder gallery table and
# SVG). Exact rational reconstruction 2026-08-15: slopes rationalized at
# denominator <= 100, then the two published triple points {0,3,5} and
# {0,8,10} imposed exactly by solving b5 and b10 from the concurrency
# equalities (linear in intercepts for fixed slopes).
N12_LOWER_BOUND_LINES = (
    (Fr(0), Fr(13201, 73)), (Fr(14, 59), Fr(10036, 67)),
    (Fr(49, 93), Fr(7519, 56)), (Fr(44, 43), Fr(-949, 94)),
    (Fr(21, 13), Fr(-8694, 43)), (Fr(169, 49), Fr(-6845645193, 14794472)),
    (Fr(357, 10), Fr(-701583, 95)), (Fr(-127, 32), Fr(70201, 79)),
    (Fr(-55, 32), Fr(4457, 8)), (Fr(-1), Fr(19303, 51)),
    (Fr(-31, 42), Fr(9623947, 28105)), (Fr(-7, 32), Fr(12252, 49)),
)
N12_LOWER_BOUND_TRIANGLES = (
    (0, 1, 2), (0, 1, 7), (0, 3, 7), (0, 3, 9), (0, 4, 10), (0, 4, 11),
    (0, 6, 9), (0, 6, 10), (1, 3, 8), (1, 3, 10), (1, 4, 11), (1, 5, 9),
    (1, 5, 10), (1, 6, 8), (1, 6, 11), (1, 7, 9), (2, 3, 4), (2, 3, 6),
    (2, 5, 6), (2, 5, 8), (2, 7, 8), (2, 7, 10), (2, 9, 10), (2, 9, 11),
    (3, 5, 7), (3, 5, 9), (3, 6, 11), (3, 8, 11), (4, 5, 7), (4, 6, 7),
    (4, 6, 9), (4, 8, 9), (4, 8, 10), (5, 8, 11), (5, 10, 11), (6, 8, 10),
    (7, 9, 11), (7, 10, 11),
)


def verify_n12_lower_bound():
    """Verify the explicit rational certificate K_gen(12) >= 38."""
    ms, bs = zip(*N12_LOWER_BOUND_LINES)
    return verify_selection(12, list(ms), list(bs),
                            N12_LOWER_BOUND_TRIANGLES, minimum=38)


# ========================== instance driver ================================

def decide(n, target, rules=("R1", "R2", "R3", "R4"), solver="cadical195",
           quiet=False):
    """In-process decision (no symmetry breaking): SAT model or None."""
    t0 = time.time()
    cnf, pool = build_model(n, target, rules)
    sat = Solver(name=solver, bootstrap_with=cnf.clauses)
    res = sat.solve()
    if not quiet:
        print(f"n={n} target={target}: {'SAT' if res else 'UNSAT'} "
              f"vars={pool.top} clauses={len(cnf.clauses)} "
              f"({time.time()-t0:.1f}s)")
    return res, (sat.get_model() if res else None), pool


def dihedral_orbits_3(n):
    """Orbit representatives of 3-subsets of Z_n under the dihedral group
    (shift i -> i+1 realized by plane rotation; reversal by the x-flip)."""
    unseen = set(combinations(range(n), 3)); reps = []
    while unseen:
        t = min(unseen)
        orb = {tuple(sorted(((-x if refl else x)+shift) % n for x in t))
               for shift in range(n) for refl in (False, True)}
        reps.append(min(orb)); unseen -= orb
    return reps


def dump_instances(n, target, out_dir, cubes=True):
    """Write the proof-relevant DIMACS family.

    The base is ``build_model`` plus three necessary accelerators: one WLOG
    sigma bit, exact selection, and the bounded-face degeneracy penalty.
    With ``cubes=True`` it writes this exhaustive family:

      * simple: all pairs cross and no concurrency;
      * par_noc: lines 0,1 parallel and no concurrency;
      * par_conc: lines 0,1 parallel and some concurrency;
      * c_i: all pairs cross and one concurrency from each D_n orbit.

    A parallel class is cyclically contiguous in slope order, so a plane
    rotation puts a parallel pair at {0,1}.  With no concurrency,
    3T <= n(n-2)-2Q; when this forces Q=1, ``par_noc`` fixes every other pair
    as crossing.  If all pairs cross, rotation/reflection moves any concurrent
    triple to a D_n representative.  The label-preserving b -> -b involution
    keeps ``add_symbreak`` sound in every cube.
    """
    os.makedirs(out_dir, exist_ok=True)
    cnf, pool = build_model(n, target)
    add_symbreak(cnf, pool)
    add_face_bound(cnf, pool, n, target)
    add_exact_selection(cnf, pool, n, target)

    sp2 = lambda *a: tuple(sorted(a))
    Pid = lambda i, j: pool.id(("P",)+sp2(i, j))
    Cid = lambda t: pool.id(("C",)+tuple(sorted(t)))
    all_pairs = list(combinations(range(n), 2))
    all_triples = list(combinations(range(n), 3))
    allP = [[Pid(i, j)] for i, j in all_pairs]
    noC = [[-Cid(t)] for t in all_triples]

    paths = {}
    prefix = os.path.join(out_dir, f"kobon_n{n}_t{target}_proof")
    base = prefix+".cnf"
    cnf.to_file(base)
    paths["base"] = base
    if not cubes:
        return paths

    c = cnf.copy()
    c.extend(allP)
    c.extend(noC)
    add_no_concurrency_capacities(c, pool, n, [0]*n)
    p = prefix+"_simple.cnf"
    c.to_file(p)
    paths["simple"] = p

    c = cnf.copy()
    c.append([-Pid(0, 1)])
    c.extend(noC)
    q_max = (n*(n-2)-3*target)//2
    if q_max < 1:
        c.append([])
    elif q_max == 1:
        c.extend([[Pid(i, j)] for i, j in all_pairs if (i, j) != (0, 1)])
    add_no_concurrency_capacities(c, pool, n, [1, 1]+[0]*(n-2))
    p = prefix+"_par_noc.cnf"
    c.to_file(p)
    paths["par_noc"] = p

    c = cnf.copy()
    c.append([-Pid(0, 1)])
    c.append([Cid(t) for t in all_triples])
    p = prefix+"_par_conc.cnf"
    c.to_file(p)
    paths["par_conc"] = p

    for idx, rep in enumerate(dihedral_orbits_3(n)):
        c = cnf.copy()
        c.extend(allP)
        c.append([Cid(rep)])
        p = prefix+f"_c{idx}.cnf"
        c.to_file(p)
        paths[f"c{idx}"] = p
    return paths
