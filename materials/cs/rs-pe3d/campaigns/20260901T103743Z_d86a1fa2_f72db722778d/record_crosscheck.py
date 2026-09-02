import json
import itertools

P = NAX = 13


def subgroup(p, m):
    return sorted(a for a in range(1, p) if pow(a, m, p) == 1)


def kron_H(s):
    S = [subgroup(P, si) for si in s]
    N = s[0] * s[1] * s[2]
    Hs = [[[pow(x, k, P) for x in ax] for k in range(len(ax) - 1)] for ax in S]
    rows = []
    for r0 in Hs[0]:
        for r1 in Hs[1]:
            for r2 in Hs[2]:
                row = [0] * N
                for i0 in range(s[0]):
                    for i1 in range(s[1]):
                        for i2 in range(s[2]):
                            row[(i0 * s[1] + i1) * s[2] + i2] = (
                                r0[i0] * r1[i1] * r2[i2]) % P
                rows.append(row)
    return rows, S, N


def rank(rows, nc):
    M = [r[:] for r in rows]
    r = 0
    for c in range(nc):
        pr = -1
        for i in range(r, len(M)):
            if M[i][c] % P:
                pr = i
                break
        if pr < 0:
            continue
        M[r], M[pr] = M[pr], M[r]
        inv = pow(M[r][c], -1, P)
        M[r] = [(x * inv) % P for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % P:
                f = M[i][c]
                M[i] = [(M[i][j] - f * M[r][j]) % P for j in range(nc)]
        r += 1
        if r == len(M):
            break
    return r


def census(s):
    H, S, N = kron_H(s)
    d = min(s)
    tot = line = off = 0
    offs = []
    for Sq in itertools.combinations(range(N), d):
        sub = [[h[j] for j in Sq] for h in H]
        dd = len(Sq) - rank(sub, len(Sq))
        assert dd >= 0
        if dd > 0:
            tot += 1
            pts = [(j // (s[1] * s[2]), (j // s[2]) % s[1], j % s[2]) for j in Sq]
            ln = False
            for ax in (0, 1, 2):
                oth = [k for k in range(3) if k != ax]
                if all(p_[oth[0]] == pts[0][oth[0]] and p_[oth[1]] == pts[0][oth[1]]
                       for p_ in pts):
                    ln = True
            if ln:
                line += 1
            else:
                off += 1
                offs.append(list(Sq))
    return tot, line, off, offs


def verify_witnesses():
    raw = open("campaigns/20260901T103743Z_d86a1fa2_f72db722778d/witnesses.json").read()
    obj = json.loads(raw + "}")  # record defect RD-1: file missing one '}' closer
    e0 = obj["censuses"]["E0"]
    fw = e0["witness_primary"]["full_word"]
    s = (2, 2, 4)
    H, S, N = kron_H(s)
    for tag in ("sum0_1_ot_h_ot_delta1", "sum1_g_ot_1_ot_delta1", "sum"):
        w = fw[tag]
        assert all(sum(h[j] * w[j] for j in range(N)) % P == 0 for h in H), tag
    s0, s1, sm = fw["sum0_1_ot_h_ot_delta1"], fw["sum1_g_ot_1_ot_delta1"], fw["sum"]
    assert [(a + b) % P for a, b in zip(s0, s1)] == sm
    assert sorted(j for j in range(N) if sm[j]) == e0["witness_primary"]["pair_flat"] == [0, 12]
    # sum0 axis-1-constant, support in x2=1 plate; sum1 axis-0-constant
    for a0 in range(2):
        for a2 in range(4):
            assert len({s0[(a0 * 2 + a1) * 4 + a2] for a1 in range(2)}) == 1
    for a1 in range(2):
        for a2 in range(4):
            assert len({s1[(a0 * 2 + a1) * 4 + a2] for a0 in range(2)}) == 1
    t, l, o, offs = census((2, 2, 4))
    assert (t, l, o) == (e0["carriers_total"], e0["lines"], e0["offline"]) == (24, 16, 8)
    print("witness word in ker(H): OK; sum support", sorted(j for j in range(N) if sm[j]))
    print("E0 independent recount:", t, l, o)
    t1, l1, o1, _ = census((2, 3, 4))
    c1 = e0["C1"]
    assert (t1, l1, o1) == (c1["carriers_total"], c1["lines"], c1["offline"]) == (12, 12, 0)
    print("C1 independent recount:", t1, l1, o1)
    t2, l2, o2, _ = census((3, 3, 4))
    c2 = e0["C2"]
    assert (t2, l2, o2) == (c2["carriers_total"], c2["lines"], c2["offline"]) == (24, 24, 0)
    print("C2 independent recount:", t2, l2, o2)
    print("RECORD CROSS-CHECK PASS (witness words in kernel; all 3 censuses reproduce)")


if __name__ == "__main__":
    verify_witnesses()
