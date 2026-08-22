"""Standalone replay of proofs/alternation_l5_char0.md.

Recomputes everything from the generator list alone (no experiment imports):
1. generator-monotone Pauli-string closure BFS at L=5 (exact ints);
2. set-law containment and count equality (reached == S_law = 262656);
3. signed alternating-form system over the 256 checkerboard pairs: full
   262656 x 512 equation grid, zero contradictions, one connected component;
4. sampled re-verification of the extracted form against the sign equations.
Run: PYTHONPATH=src python tests/test_alternation_l5.py
"""
import random

N = 10
J = (1 << N) - 1
D = 0b0101010101  # checkerboard on the 2x5 grid, odd weight, bit0 = 1
PAR = lambda a: bin(a).count("1") & 1


def main():
    idx = lambda r, c: r * 5 + c
    bonds = [(idx(r, c), idx(r, c + 1)) for r in range(2) for c in range(4)] + \
            [(idx(0, c), idx(1, c)) for c in range(5)]
    gens = [(1 << i, 0) for i in range(N)] + [(0, (1 << u) | (1 << v)) for u, v in bonds]
    assert len(gens) == 23

    def symp(a, b):
        return PAR(a[0] & b[1]) ^ PAR(a[1] & b[0])

    reached = set(gens)
    queue = list(reached)
    while queue:
        s = queue.pop()
        for g in gens:
            if symp(s, g):
                t = (s[0] ^ g[0], s[1] ^ g[1])
                if t not in reached:
                    reached.add(t)
                    queue.append(t)
    assert len(reached) == 262656, len(reached)
    print("ok: closure BFS reaches exactly 262656 strings")

    def in_law(s):
        x, z = s
        if PAR(z) != 0:
            return False
        if z == J:
            return True
        return PAR((J ^ z) & x) == 1 ^ PAR(z & D)

    assert all(in_law(s) for s in reached)
    law_count = sum(1 for z in range(1 << N) if PAR(z) == 0 and z != J) * (1 << (N - 1)) + (1 << N)
    assert law_count == 262656
    print("ok: every reached string satisfies the set law; counts match (set equality)")

    def crep(a):
        a &= 1023
        return a if (a >> 9) & 1 == 0 else a ^ J

    cid, classes, seen = {}, [], set()
    for a in range(1024):
        r = crep(a)
        if r not in seen:
            seen.add(r)
            cid[r] = len(classes)
            classes.append(r)
    assert len(classes) == 512
    picls = [cid[crep(r ^ D)] for r in classes]
    oriented = [(r & 1) == 0 for r in classes]
    assert all(oriented[picls[i]] != oriented[i] for i in range(512))
    orient = [i if oriented[i] else picls[i] for i in range(512)]
    t_bool = [0 if oriented[i] else 1 for i in range(512)]
    oidx = {}
    c = 0
    for i in range(512):
        if oriented[i]:
            oidx[i] = c
            c += 1

    parent = list(range(256))
    flp = [0] * 256

    def find(a):
        x = 0
        while parent[a] != a:
            x ^= flp[a]
            a = parent[a]
        return a, x

    def unite(a, b, w):
        ra, xa = find(a)
        rb, xb = find(b)
        if ra == rb:
            return (xa ^ xb) == w
        parent[ra] = rb
        flp[ra] = xa ^ xb ^ w
        return True

    edges = contradictions = 0
    for s in reached:
        x, z = s
        for Ai in range(512):
            a = classes[Ai]
            B = cid[crep(a ^ x)]
            w = PAR(z & a) ^ t_bool[B] ^ PAR(z & classes[picls[B]]) ^ t_bool[Ai] ^ 1
            edges += 1
            if not unite(oidx[orient[B]], oidx[orient[Ai]], w):
                contradictions += 1
                break
        if contradictions:
            break
    assert edges == 262656 * 512 == 134479872
    assert contradictions == 0
    print(f"ok: {edges} sign equations solved, 0 contradictions")

    roots = {find(i)[0] for i in range(256)}
    assert len(roots) == 1
    print("ok: alternating form unique up to scalar (1 component)")

    q = {i: find(oidx[i])[1] for i in oidx}
    sig = [(-1) ** q[orient[i]] for i in range(512)]
    random.seed(7)
    sample = random.sample(sorted(reached), 40)
    for s in sample:
        x, z = s
        for Ai in random.sample(range(512), 8):
            a = classes[Ai]
            B = cid[crep(a ^ x)]
            lhs = ((-1) ** PAR(z & a)) * (1 if oriented[B] else -1) * sig[orient[B]] \
                + ((-1) ** PAR(z & classes[picls[B]])) * (1 if oriented[Ai] else -1) * sig[orient[Ai]]
            assert lhs == 0
    print("ok: extracted form independently satisfies sampled equations")
    print("PASS: g_{2x5}(Q) = sp_512(Q) (+) sp_512(Q), dim 262656, exact char-0 replay")


if __name__ == "__main__":
    main()
