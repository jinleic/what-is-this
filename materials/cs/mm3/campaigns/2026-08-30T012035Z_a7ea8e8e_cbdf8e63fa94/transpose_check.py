"""Resolve the 57-vs-59 tension on MWS-59: transposition sanity check.

My printed-Table-2 output schedule counts 28 (9 v-gates + 19 C-assembly
additions). The Table-3 W block has d(W)=14 floor-impossible, so
C(Wfac) >= 15 and the transposition principle demands output >= 29.
These contradict, so something is wrong. This script builds the transposed
circuit from the output schedule CONSTRUCTIVELY and checks what it computes
and how many gates it needs.
"""
from collections import defaultdict

R = 23
D = 9

v_gates = [
    ("v0", [("M0", 1), ("M12", -1)]),
    ("v1", [("M16", 1), ("v0", 1)]),
    ("v2", [("M11", 1), ("M21", -1)]),
    ("v3", [("M14", 1), ("M13", -1)]),
    ("v4", [("M18", 1), ("v1", -1)]),
    ("v5", [("M2", 1), ("v4", 1)]),
    ("v6", [("M4", 1), ("M9", 1)]),
    ("v7", [("M15", 1), ("v3", 1)]),
    ("v8", [("M17", 1), ("v2", -1)]),
]
C_list = [
    [("M6", 1), ("M13", 1), ("M20", 1)],
    [("M0", 1), ("M1", 1), ("M5", 1)],
    [("M3", 1), ("M10", 1), ("v5", 1)],
    [("v6", 1), ("v8", 1)],
    [("M8", 1), ("v1", -1), ("v2", 1), ("v6", -1), ("v7", 1)],
    [("M9", 1), ("v4", -1), ("v7", -1)],
    [("M7", -1), ("M22", 1), ("v3", -1), ("v8", -1)],
    [("M7", 1), ("M19", 1), ("M21", 1), ("v0", 1)],
    [("M7", -1), ("v5", 1)],
]


def main():
    # forward value of every wire as an exact vector over the 23 products
    val = {f"M{r}": tuple(1 if j == r else 0 for j in range(R)) for r in range(R)}
    for name, terms in v_gates:
        acc = [0] * R
        for atom, sign in terms:
            src = val[atom]
            for j in range(R):
                acc[j] += sign * src[j]
        val[name] = tuple(acc)
    # output k's value:
    outs = []
    for k, terms in enumerate(C_list):
        acc = [0] * R
        for atom, sign in terms:
            src = val[atom]
            for j in range(R):
                acc[j] += sign * src[j]
        outs.append(tuple(acc))
    # out[k] must be the W row k (C-entry weights over products)
    # sanity: out[0] = M6+M13+M20 exactly.

    # ---- transposed circuit via reverse accumulation ----
    # tau: wire -> vector over the 9 outputs. Start: tau(C_k)=e_k.
    # Reverse the combined schedule (v-gates then final assemblies as gates).
    # Model final assembly of C_k as a chain of binary adds: building
    # C_k = t1 + t2 + ... incrementally, each step is a forward gate.
    # For transposition we treat each C_k's chain as gates too.
    gates = []  # (name, [(sign, atom)...2], ) binary only
    # expand C_k chains: sequential accumulation in printed order
    chain = {}
    for k, terms in enumerate(C_list):
        wires = [a for a, _ in terms]
        signs = [s for _, s in terms]
        cur = wires[0]
        cs = signs[0]
        for nxt, sn in zip(wires[1:], signs[1:]):
            gname = f"g{k}_{nxt}"
            gates.append((gname, cs, cur, sn, nxt))
    fwd = list(v_gates) + [(n, [(s1, a), (s2, b)]) for (n, s1, a, s2, b) in gates]

    # consumers: wire -> list of (gate name, its sign toward wire)
    consumers = defaultdict(list)
    for name, terms in fwd:
        for sign, atom in terms:
            consumers[atom].append((name, sign))
    # outputs: wire -> (C_k, sign)
    for k, (w, s) in chain.items():
        consumers[w].append((f"C{k}", s))

    # reverse pass
    tau = {f"C{k}": tuple(1 if j == k else 0 for j in range(D)) for k in range(D)}
    t_adds = []
    for name, terms in reversed(fwd):
        # tau(name) = sum over consumers (sign * tau(consumer)); assemble
        # with binary adds at cost (#consumers - 1) if #consumers > 0
        acc = None
        cs = []
        for cname, sign in consumers.get(name, []):
            t = tau.get(cname)
            if t is None:
                continue
            cs.append((sign, t))
        if not cs:
            tau[name] = None  # dead wire
            continue
        s0, t0v = cs[0]
        cur = tuple(s0 * x for x in t0v)
        for s, t in cs[1:]:
            cur = tuple(cur[j] + s * t[j] for j in range(D))
            t_adds.append(name)
        tau[name] = cur
    # tau(M_r) = Wfac row r
    Wfac = []
    for r in range(R):
        t = tau.get(f"M{r}")
        Wfac.append(list(t) if t is not None else [0] * D)
    # compare with the Table-3 W block (my W): W2 rows = C entries, so
    # Table-3 W = (my W2). Wfac rows r should equal column r of my W2.
    # verify Wfac == W2^T (as vectors over 9 entries)
    print("transposed gates:", len(t_adds))
    print("expected L' = L + m - n where L =", 9 + 19, " m=23 n=9 ->", 9 + 19 + 14, "(C-chain expansions included)")

    # Count of the transposed circuit per the standard model:
    # additions inside tau-gates = #t_adds. Plus the inputs of the transposed
    # circuit are the 9 C-entries (free). No other cost. So transposed circuit
    # has len(t_adds) additions.
    print("Wfac circuit additions =", len(t_adds), "(bound claim: >= 15)")
    # verify exactness: compare Wfac with the direct W2^T
    # rebuild W2 (C-entry weights over products) from the forward values:
    W2T = []
    for r in range(R):
        row = []
        for k in range(D):
            row.append(outs[k][r])
        W2T.append(row)
    same = all(tuple(Wfac[r]) == tuple(W2T[r]) for r in range(R))
    print("transposed circuit computes exactly the W factor map:", same)


if __name__ == "__main__":
    main()
