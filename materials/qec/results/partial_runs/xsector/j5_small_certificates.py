"""J.5 attempt — small-lattice certificates (3x3, 4x4, 3x5).

Machine-verified Route A lemmas + Route B rank-bound data:
  (A1) Zcen(Q) = ker_col H_X  =>  d_Z(Q) >= d_Z(P)     (Z-monotonicity)
  (A2) violation window: d_SX(P) <= w_dem < d_Z(P)     [definitional]
  (A3) cols of M are R-translates of col_0: col_i=0 <=> col_0=0 <=> M=0, hence
       single-row demotion requires M = 0               (sharpens J.4(iii))
  (A4) dim{(C,D): M=0} = 2*dim - rank H_Z; trivial family (B^T g, A^T g) demotes
       nothing; sigma = dim - rank H_Z = k_P when rank H_X = rank H_Z
  (B)  T_X vs rho_X tightness instances + dimension-twin non-forceability search
Also runs the sharp light-stabilizer probe on every parent with d_SX(P) < d_X(P)
(ell*m in {9,15,16} was NOT covered exhaustively by hunts 1-3).

GF(2)/numpy only, single thread, no SAT.
"""
from __future__ import annotations

import itertools
import json
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE.parent))

from qec_research.gf2.linalg import rank_np, nullspace_np  # noqa: E402
import j5core  # noqa: E402

OUT = HERE.parent / "j5_small_certificates.json"
MAX_PARENTS_PER_LATTICE = 300
ENUM_CAP = 20  # span enumerations max 2^cap


def span_enum(basis, fn):
    """Iterate the span of row-basis, calling fn(Vchunk)."""
    k = basis.shape[0]
    for start in range(1, 1 << k, 1 << 12):
        idx = np.arange(start, min(start + (1 << 12), 1 << k), dtype=np.int64)
        V = np.zeros((idx.size, basis.shape[1]), dtype=np.uint8)
        for j in range(k):
            V ^= (((idx >> j) & 1)[:, None] & basis[j]).astype(np.uint8)
        fn(V, idx)
        if start + (1 << 12) >= 1 << k:
            break


def min_weight_rowspace(R):
    """Min over the CODE rowspace excluding the zero vector (independent basis)."""
    Rr, pivots = nullspace_np(R), None
    # independent basis via rref rows
    from qec_research.gf2.linalg import rref_np
    Rr, _ = rref_np(R)
    best = [None]

    def fn(V, _):
        w = V.sum(axis=1)
        w = w[w > 0]
        if w.size:
            wm = int(w.min())
            if best[0] is None or wm < best[0]:
                best[0] = wm
    span_enum(Rr, fn)
    return best[0]


def parent_distance_TX(parent: j5core.ParentData):
    """Exact d_X(P), min-weight X-logicals count, and T_X(P). None if too big."""
    dim, N, lat = parent.dim, parent.N, parent.lat
    kerHZ = nullspace_np(parent.HZ)          # Xcen(P)
    dk = kerHZ.shape[0]
    if dk > ENUM_CAP:
        return None
    ns_SX = nullspace_np(parent.HX)          # x not in S_X(P) <=> x has odd dot here
    best = [None]

    def pass1(V, _):
        w = V.sum(axis=1)
        logical = ((V @ ns_SX.T) % 2).any(axis=1)
        if logical.any():
            wm = int(w[logical].min())
            if best[0] is None or wm < best[0]:
                best[0] = wm
    span_enum(kerHZ, pass1)
    if best[0] is None:
        return None
    mins = []

    def pass2(V, _):
        w = V.sum(axis=1)
        logical = ((V @ ns_SX.T) % 2).any(axis=1)
        sel = np.nonzero(logical & (w == best[0]))[0]
        mins.extend([V[i] for i in sel[:2048 - len(mins)]])
    if len(mins) < 2048:
        span_enum(kerHZ, pass2)
    orb = [parent.HX]
    for v in mins:
        orb.append(np.hstack([lat.sparse(v[:dim]), lat.sparse(v[dim:])]))
    MX = np.vstack(orb)
    TX = int(rank_np(MX) - parent.rHX)
    return {"d_X": int(best[0]), "T_X": TX, "min_logical_count": len(mins)}


def enumerate_parents(ell, m, rng, cap):
    grid = [(a, b) for a in range(ell) for b in range(m)]
    sets = [tuple(sorted(c)) for w in (1, 2, 3) for c in itertools.combinations(grid, w)]
    lat = j5core.Lattice(ell, m)
    order = rng.permutation(len(sets))
    out = []
    seen = set()
    for ai in order:
        At = sets[ai]
        for bt in sets:
            key = (At, bt)
            if key in seen:
                continue
            seen.add(key)
            p = j5core.ParentData(lat, At, bt)
            if p.kP < 2:
                continue
            out.append(p)
            if len(out) >= cap:
                return out
    return out


def analyse_instance(p, cd, dZ, ET):
    """Full demotion analysis of one valid (C,D). Returns record or None."""
    dim, N = p.dim, p.N
    C, D = j5core.split_CD(p, cd)
    sp = j5core.demotion_spaces(p, C, D)
    M = sp["M"]
    rec = {"valid": not bool((M ^ M.T).any())}
    if not rec["valid"]:
        return rec
    # (A3): columns of M are translates of col_0; hence no zero column unless M=0
    c0 = M[:, [0]]
    rec["circulant_cols_ok"] = all(
        np.array_equal((p.lat.X[k] @ c0) % 2, M[:, [k]]) for k in range(1, dim))
    rec["zero_col_with_nonzero_M"] = bool(
        M.any() and any(not M[:, i].any() for i in range(dim)))
    rec["M0"] = not bool(M.any())
    rec["delta_bar"] = sp["delta_bar"]
    rec["rho_X"] = int(rank_np(np.vstack([p.HZ, sp["CD"]])) - p.rHZ)
    kerM = sp["kerM"]
    rec["dim_kerM"] = int(kerM.shape[0])
    if rec["M0"] and kerM.shape[0] <= ET:
        # (A4): W = full rowspace HX  <=> rank(kerM @ HX) == rank HX
        rec["M0_W_full"] = int(rank_np((kerM @ p.HX) % 2)) == int(p.rHX)
    # exact demotion analysis
    dk = kerM.shape[0]
    if 0 < dk <= ET:
        probe = sp["probe_SZD"]
        CD = sp["CD"]
        HX, dimp = p.HX, p.dim
        best = [None]

        def fn(V, _):
            Lx = V
            Vx = (Lx @ HX) % 2
            Vz = (Lx @ CD) % 2
            w = Vx.sum(axis=1)
            dem = ((Vz @ probe.T) % 2).any(axis=1) & (w > 0)
            if dem.any():
                wm = int(w[dem].min())
                if best[0] is None or wm < best[0]:
                    best[0] = wm
        span_enum(kerM, fn)
        rec["w_dem"] = best[0]
    else:
        rec["w_dem"] = None
    return rec


def main():
    t0 = time.monotonic()
    rng = np.random.default_rng(0x5A5)
    out = {"schema": "j5-small-certificates-v1",
           "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "lattices": [], "sharp": {"probes": 0, "candidates": 0, "violations": 0},
           "violations": 0}
    for (ell, m) in [(3, 3), (3, 5), (4, 4)]:
        L = {"ell": ell, "m": m, "parents": 0, "instances": 0, "valid": 0,
             "checks": defaultdict(int), "sigma_stats": defaultdict(int),
             "window_parents": [], "tight_TX": [], "twins": []}
        twinmap = defaultdict(dict)
        parents = enumerate_parents(ell, m, rng, MAX_PARENTS_PER_LATTICE)
        L["parents"] = len(parents)
        for p in parents:
            if p.rHX == p.rHZ:
                L["checks"]["rank_eq_HX_HZ"] += 1
            L["checks"]["sigma_eq_kP"] += int(p.sigma == p.kP)
            L["sigma_stats"][f"sigma={p.sigma}"] += 1
            dd = parent_distance_TX(p)
            dSX = min_weight_rowspace(p.HX) if p.rHX <= ENUM_CAP else None
            dZ = dd["d_X"] if dd else None
            TX = dd["T_X"] if dd else None
            win = (dSX is not None and dZ is not None and dSX < dZ)
            # instance sampling: validity kernel combos + syzygy elements
            V = p.validity
            Vd = V.shape[0]
            samples = []
            if Vd <= 10:
                for s in range(1, 1 << Vd):
                    sel = np.array([(s >> j) & 1 for j in range(Vd)], np.uint8)
                    samples.append((sel @ V) % 2)
            else:
                for _ in range(200):
                    sel = rng.integers(0, 2, Vd).astype(np.uint8)
                    if sel.any():
                        samples.append((sel @ V) % 2)
            S = p.syzygy
            for r in range(S.shape[0]):
                samples.append(S[r])
            if S.shape[0]:
                for _ in range(40):
                    sel = rng.integers(0, 2, S.shape[0]).astype(np.uint8)
                    if sel.any():
                        samples.append((sel @ S) % 2)
            seen = set()
            margins = []
            for cd in samples:
                C0, D0 = j5core.split_CD(p, cd)
                key = (C0.tobytes(), D0.tobytes())
                if key in seen:
                    continue
                seen.add(key)
                rec = analyse_instance(p, cd, dZ, ENUM_CAP)
                L["instances"] += 1
                if not rec["valid"]:
                    continue
                L["valid"] += 1
                if rec["circulant_cols_ok"]:
                    L["checks"]["circulant_cols"] += 1
                if not rec["zero_col_with_nonzero_M"]:
                    L["checks"]["no_zero_col_unless_M0"] += 1
                if rec.get("M0_W_full"):
                    L["checks"]["M0_implies_W_full"] += 1
                wm = rec.get("w_dem")
                if TX is not None and rec["rho_X"] == TX and len(L["tight_TX"]) < 25:
                    L["tight_TX"].append(
                        {"A": p.A_terms, "B": p.B_terms, "T_X": TX, "rho_X": rec["rho_X"], "margin": (wm - dZ) if (wm is not None and dZ is not None) else None,
                         "w_dem": wm, "d_XP": dZ, "cd": cd.tolist()})
                if wm is not None and dZ is not None:
                    margins.append((wm, cd.copy(), rec))
                    if wm < dZ:
                        out["violations"] += 1
                        L.setdefault("VIOLATIONS", []).append(
                            {"A": p.A_terms, "B": p.B_terms, "w_dem": wm, "d_XP": dZ,
                             "cd": cd.tolist()})
                    # twin profile: all dimension data only
                    pkey = (int(p.kP), int(rec["delta_bar"]), int(rec["rho_X"]),
                            int(rec["dim_kerM"]), dSX, dZ, TX)
                    twinmap[pkey][wm - dZ] = (list(p.A_terms), list(p.B_terms),
                                              cd.tolist(), wm)
            if win:
                lights = j5core.light_stabilizers(p, dZ)
                probes, cands = j5core.sharp_probe(p, lights, dZ, rng)
                out["sharp"]["probes"] += len(probes)
                out["sharp"]["candidates"] += len(cands)
                vered = []
                for cnd in cands:
                    ver = j5core.verify_full(p, cnd["cd"], cnd["lam"], dZ)
                    if ver["verified"]:
                        vered.append({"cand": cnd, "verify": ver})
                if vered:
                    out["sharp"]["violations"] += len(vered)
                    out.setdefault("sharp_witnesses", []).append(
                        {"ell": ell, "m": m, "A": p.A_terms, "B": p.B_terms,
                         "d_XP": dZ, "witnesses": vered})
                L["window_parents"].append({
                    "A": p.A_terms, "B": p.B_terms, "k_P": int(p.kP),
                    "d_SX": dSX, "d_XP": dZ, "T_X": TX, "sigma": int(p.sigma),
                    "n_lights": len(lights), "probes": len(probes),
                    "cands": len(cands), "verified": len(vered),
                    "min_margin": min((w - dZ for (w, _, __) in margins), default=None)})
        for key, margin_d in twinmap.items():
            if len(margin_d) >= 2:
                items = list(margin_d.items())[:3]
                L["twins"].append({"profile": dict(zip(
                    ["k_P", "delta_bar", "rho_X", "dim_kerM", "d_SX", "d_X", "T_X"], key)),
                    "data": [{"margin": mg, "parent": v[0:2], "cd": v[2], "w_dem": v[3]}
                             for mg, v in items]})
        L["checks"] = dict(L["checks"])
        L["sigma_stats"] = dict(L["sigma_stats"])
        L["twins"] = L["twins"][:15]
        out["lattices"].append(L)
        print(f"lattice {ell}x{m}: parents={L['parents']} valid={L['valid']} "
              f"viol={out['violations']} sharpV={out['sharp']['violations']} "
              f"twins={len(L['twins'])} tight={len(L['tight_TX'])} "
              f"t={time.monotonic()-t0:.0f}s", flush=True)
        out["wall_seconds_partial"] = round(time.monotonic() - t0, 1)
        Path(str(OUT)).write_text(json.dumps(out, indent=1, default=str))
    out["wall_seconds"] = round(time.monotonic() - t0, 1)
    out.pop("wall_seconds_partial", None)
    Path(str(OUT)).write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps({"violations": out["violations"],
                      "sharp_violations": out["sharp"]["violations"],
                      "seconds": out["wall_seconds"]}, default=str))


if __name__ == "__main__":
    main()
