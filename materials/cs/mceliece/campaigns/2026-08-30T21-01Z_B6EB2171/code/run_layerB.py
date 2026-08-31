"""Gate B Layer B runner — adversarial + exhaustive G-space probes.

Registered probes (pre_statement + ADDENDUM 1, section 2-fix / A1.3):

  E1  EXHAUSTIVE census at (m, n, t) = (6, 64, 2): EVERY monic irreducible
      degree-2 polynomial over F_64 — exactly (64^2 - 64)/2 = 2016 of them —
      full census of the G-population at this cell.  A degeneracy rate over
      an EXHAUSTED population (census, not sample).  Support = full field.
      Y := RREF rows (deterministic; ADDENDUM A1.1 proves the verdict is
      G-, m-, t-dependent only — but the instance data is recorded per G
      regardless).

  E2  special-G families at registered cells:
      (i)  m=6, t=3: G = Z^3 + c for every c in F_64, keep the irreducibles
           (those c that are not cubes; gcd(3,63)=3 divides, so exactly the
           non-cubes: 42 of 64).
      (ii) m=6, t=3: G = Z^3 + Z + c for every c in F_64, keep irreducibles.
      (iii) m=7, t=2: G = Z^2 + Z + c for every c in F_128; irreducible iff
           Tr(c) = 1 (64 candidates).

  V1  invariance VERIFICATION probes (ADDENDUM A1.1 re-scope):
      (a) support-reordering: at (6,64,3), seed 1387's instance with its
          columns/support consistently permuted (3 distinct permutations).
      (b) Y-basis change: same instance, first two RREF rows swapped and
          with row0 -> row0 + row1 (still spans the SAME code -> per
          A1.1(c) the interpolant-span is unchanged -> verdict unchanged).
      Every probe must return the same verdict as the parent instance.

Checkpoint-after-every-record; DEGENERATE halts everything for escalation
(registered rule).  All arithmetic exact F_{2^m}.
"""
from __future__ import annotations
import json, os, sys, time, hashlib

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from gfield import GF
from fastfield import EField, vec_pis_irreducible
from instance import Instance
from engine import engine_verdict, cascade_verdict, brute_pair_values

OUT = os.path.dirname(HERE)
STATE = os.path.join(OUT, "layerB_state.json")
LOG = os.path.join(OUT, "layerB_log.jsonl")


def append_log(rec):
    with open(LOG, "a") as f:
        f.write(json.dumps(rec, default=str) + "\n")


def load_state():
    if os.path.exists(STATE):
        return json.load(open(STATE))
    return {"probes": {}, "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}


def save_state(st):
    tmp = STATE + ".tmp"
    json.dump(st, open(tmp, "w"), indent=1)
    os.replace(tmp, STATE)


G6 = GF(6)
E6 = EField(6)
G7 = GF(7)
E7 = EField(7)


def measure(inst, gf, ef, brute_confirm=True):
    ev = engine_verdict(inst.F, inst.Y, np.asarray(inst.lam, np.uint16),
                        inst.support, ef, gf)
    out = {"engine": {kk: ev[kk] for kk in ev if kk != "beta_at_points"}}
    if ev["verdict"] == "DEGENERATE" and brute_confirm:
        cv = cascade_verdict(gf, [gf.ptrim(list(f)) for f in inst.F])
        out["brute_recheck"] = cv
        out["verdict"] = "DEGENERATE" if cv["verdict"] == "DEGENERATE" else "INSTRUMENT_DISAGREEMENT"
        if out["verdict"] == "DEGENERATE":
            out["G"] = inst.G
            out["support"] = inst.support
            out["lam"] = [int(x) for x in inst.lam]
            out["F_degrees"] = [gf.pdeg(f) for f in inst.F]
    elif ev["verdict"] == "NONDEGENERATE":
        out["verdict"] = "NONDEGENERATE"
    else:
        out["verdict"] = ev["verdict"]
    return out


# ---------------------------------------------------------------- E1

def enum_irreducible_deg2_m6():
    """ALL monic irreducible degree-2 polys over F_64: Z^2 + aZ + b with
    discriminant condition in char 2: Z^2 + aZ + b irreducible iff a != 0
    and Tr(b/a^2) = 1.  Counted independently: (q^2-q)/2 = 2016."""
    q = 64
    # trace table over F_64: Tr(x) = sum_{i=0}^{5} x^{2^i}
    tr = np.zeros(q, dtype=np.uint8)
    for x in range(q):
        v, acc = x, 0
        for _ in range(6):
            acc ^= v
            v = E6.MUL[v, v]
        tr[x] = acc & 1
    out = []
    for a in range(1, q):
        a2 = int(E6.MUL[a, a])
        a2inv = int(E6.INV[a2])
        for b in range(q):
            if tr[int(E6.MUL[b, a2inv])]:
                out.append([b, a, 1])
    return out


def run_E1(st):
    if st["probes"].get("E1", {}).get("done"):
        return
    p = st["probes"].setdefault("E1", {"records": [], "next": 0, "done": False})
    Gs = enum_irreducible_deg2_m6()
    assert len(Gs) == 2016, len(Gs)
    t0 = time.time()
    n_deg = 0
    while p["next"] < len(Gs) and not st.get("halt_event"):
        G = Gs[p["next"]]
        try:
            inst = Instance(6, 64, 2, seed=777 + p["next"], forced_G=G)
            gok = all(inst.guards.get(x, False) for x in
                      ("alpha_identity", "gamma_k", "eps_gcd_const", "eps_maxdeg_is_D"))
            if not gok:
                rec = {"probe": "E1", "idx": p["next"], "G": G,
                       "verdict": "GUARD_FAIL",
                       "guard_detail": {kk: vv for kk, vv in inst.guards.items() if kk != "degs"}}
            else:
                mres = measure(inst, G6, E6)
                rec = {"probe": "E1", "idx": p["next"], "G": G,
                       "verdict": mres["verdict"], "engine_route": mres["engine"].get("route")}
                if mres["verdict"] == "DEGENERATE":
                    rec.update({kk: mres[kk] for kk in ("G", "support", "lam", "F_degrees")})
                    n_deg += 1
        except Exception as e:
            rec = {"probe": "E1", "idx": p["next"], "G": G,
                   "verdict": "BUILD_ERROR", "error": repr(e)[:200]}
        p["records"].append(rec)
        append_log(rec)
        p["next"] += 1
        save_state(st)
        if rec["verdict"] == "DEGENERATE":
            print(f"[E1] DEGENERATE at idx {p['next']-1} G={G} — HALT", flush=True)
            st["halt_event"] = rec
            save_state(st)
            return
        if p["next"] % 96 == 0:
            print(f"[E1] {p['next']}/2016 degenerate_so_far={n_deg} [{time.time()-t0:.0f}s]",
                  flush=True)
    p["done"] = p["next"] >= len(Gs)
    p["n_degenerate"] = n_deg
    p["elapsed_s"] = round(time.time() - t0, 1)
    save_state(st)
    print(f"[E1] DONE n={p['next']} degenerate={n_deg} [{p['elapsed_s']}s]", flush=True)


# ---------------------------------------------------------------- E2

def run_E2(st):
    if st["probes"].get("E2", {}).get("done"):
        return
    p = st["probes"].setdefault("E2", {"records": [], "next": 0, "done": False})
    fams = []
    # (i) Z^3 + c over F_64
    for c in range(64):
        fams.append(("E2i_Z3c", 6, 64, 3, [c, 0, 0, 1]))
    # (ii) Z^3 + Z + c over F_64
    for c in range(64):
        fams.append(("E2ii_Z3Zc", 6, 64, 3, [c, 1, 0, 1]))
    # (iii) Z^2 + Z + c over F_128
    for c in range(128):
        fams.append(("E2iii_Z2Zc_m7", 7, 128, 2, [c, 1, 1]))
    t0 = time.time()
    while p["next"] < len(fams) and not st.get("halt_event"):
        name, m, n, t, G = fams[p["next"]]
        gf, ef = (G6, E6) if m == 6 else (G7, E7)
        try:
            if not gf.pis_irreducible(G):
                rec = {"probe": name, "idx": p["next"], "G": G,
                       "verdict": "REDUCIBLE_SKIPPED"}
            else:
                inst = Instance(m, n, t, seed=555 + p["next"], forced_G=G)
                gok = all(inst.guards.get(x, False) for x in
                          ("alpha_identity", "gamma_k", "eps_gcd_const",
                           "eps_maxdeg_is_D"))
                if not gok:
                    rec = {"probe": name, "idx": p["next"], "G": G,
                           "verdict": "GUARD_FAIL",
                           "guard_detail": {kk: vv for kk, vv in inst.guards.items() if kk != "degs"}}
                else:
                    mres = measure(inst, gf, ef)
                    rec = {"probe": name, "idx": p["next"], "G": G,
                           "verdict": mres["verdict"],
                           "engine_route": mres["engine"].get("route")}
                    if mres["verdict"] == "DEGENERATE":
                        rec.update({kk: mres[kk] for kk in ("G", "support", "lam", "F_degrees")})
        except Exception as e:
            rec = {"probe": name, "idx": p["next"], "G": G,
                   "verdict": "BUILD_ERROR", "error": repr(e)[:200]}
        p["records"].append(rec)
        append_log(rec)
        p["next"] += 1
        save_state(st)
        if rec["verdict"] == "DEGENERATE":
            print(f"[{name}] DEGENERATE at G={G} — HALT", flush=True)
            st["halt_event"] = rec
            save_state(st)
            return
    p["done"] = p["next"] >= len(fams)
    from collections import Counter
    p["tally"] = dict(Counter(r["verdict"] for r in p["records"]))
    p["elapsed_s"] = round(time.time() - t0, 1)
    save_state(st)
    print(f"[E2] DONE tally={p['tally']} [{p['elapsed_s']}s]", flush=True)


# ---------------------------------------------------------------- V1

def run_V1(st):
    if st["probes"].get("V1", {}).get("done"):
        return
    p = st["probes"].setdefault("V1", {"records": [], "done": False})
    inst = Instance(6, 64, 3, 1387)
    parent = measure(inst, G6, E6)
    recs = [{"probe": "V1_parent", "verdict": parent["verdict"],
             "route": parent["engine"].get("route")}]
    rng = np.random.default_rng(4242)
    for k in range(3):
        perm = rng.permutation(64)
        inst2 = Instance.__new__(Instance)
        inst2.__dict__.update(inst.__dict__)
        inst2.support = [inst.support[i] for i in perm]
        inst2.lam = [inst.lam[i] for i in perm]
        inst2.Y = inst.Y[:, perm]
        m2 = measure(inst2, G6, E6)
        recs.append({"probe": f"V1_perm{k}", "verdict": m2["verdict"],
                     "same_as_parent": m2["verdict"] == parent["verdict"]})
    # basis change: row0 <-> row1, then row0 += row1 (span unchanged)
    for name, idxs in [("V1_swap01", [1, 0] + list(range(2, inst.k))),
                       ("V1_row0plus1", None)]:
        inst3 = Instance.__new__(Instance)
        inst3.__dict__.update(inst.__dict__)
        if idxs is not None:
            inst3.Y = inst.Y[idxs]
        else:
            Yn = inst.Y.copy()
            Yn[0] = Yn[0] ^ Yn[1]
            inst3.Y = Yn
        m3 = measure(inst3, G6, E6)
        recs.append({"probe": name, "verdict": m3["verdict"],
                     "same_as_parent": m3["verdict"] == parent["verdict"]})
    p["records"] = recs
    p["done"] = bool(all(r.get("same_as_parent", True) for r in recs))
    p["invariance_holds"] = p["done"]
    save_state(st)
    print(f"[V1] invariance_holds={p['invariance_holds']}", flush=True)


def main():
    st = load_state()
    for fn in (run_V1, run_E2, run_E1):
        if st.get("halt_event"):
            break
        fn(st)
    print(json.dumps({"status": "halt" if st.get("halt_event") else "complete",
                      "probes": {k: {kk: vv for kk, vv in v.items()
                                     if kk != "records"}
                                 for k, v in st["probes"].items()}},
                     indent=1, default=str))


if __name__ == "__main__":
    main()
