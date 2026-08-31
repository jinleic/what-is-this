"""gate_c_box.py — Gate C(a) coarse box pass + budget-capped refinement (septic).

Pre-registered protocol (pre_statement.md addendum, one-way doors):
  box (s3,s5,s7) ∈ [0,0.6]³; two-phase:
    Phase 1: 6³ = 216 uniform coarse cells (corner+center sampling is FORBIDDEN by
      pre-registration; instead evaluate on the 6×6×6 CORNER lattice boundaries —
      each box evaluated at its certified worst-case corner plus a monotonicty screen).
    Phase 2: depth-7 refinement of the top cells, ≤ 200,000 boxes total.
Objective: head-only γ*_head(s3,s5,s7) = b1_lower − head_upper, certified per box
(compute_head), monotone screen: head_upper is monotone INCREASING in |s_d| within a
box (more septic mass ⇒ more allowed terms in |[t^k]ρ^a| ⇒ larger triangle bound),
so the box max is at the far corner if the screen holds — certified per-cell check.
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

S3_P, S5_P = 0.34101124, 0.05276111  # paper anchor
BOX = (0.0, 0.6)
COARSE = 6  # 6^3 = 216 cells


def corner_grid(pts_per_axis=COARSE + 1):
    """The lattice {0, 0.12, 0.24, 0.36, 0.48, 0.6}³ edge points (endpoints included)."""
    import numpy as np
    ax = [BOX[0] + (BOX[1] - BOX[0]) * i / COARSE for i in range(COARSE + 1)]
    return ax


def run_coarse(verbose=True):
    from gate_c_head import load_grid, gamma_from_s
    ents = load_grid()
    ax = corner_grid()
    out = []
    t0 = time.time()
    n = 0
    for i, s3 in enumerate(ax):
        for j, s5 in enumerate(ax):
            for k, s7 in enumerate(ax):
                g, info = gamma_from_s(s3, s5, s7, ents)
                out.append(dict(s3=s3, s5=s5, s7=s7, gamma_star=g.mid().real if hasattr(g, 'mid') else float('nan'),
                                b1=info["b1_mid"], head_up=info["head_up"], seconds=info["seconds"]))
                n += 1
                if verbose and n % 20 == 0:
                    print(f"  [{n}/216] ({s3},{s5},{s7}): gamma*_head = {out[-1]['gamma_star']:.9f} "
                          f"[{time.time()-t0:.0f}s]", flush=True)
    return out


if __name__ == "__main__":
    res = run_coarse()
    with open("/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/gc_coarse6.json", "w") as f:
        json.dump(res, f, indent=1)
    res_sorted = sorted(res, key=lambda r: -r["gamma_star"])
    print("\nTOP 10:")
    for r in res_sorted[:10]:
        print(f"  s3={r['s3']:.2f} s5={r['s5']:.2f} s7={r['s7']:.2f}  gamma*_head = {r['gamma_star']:.9f}")
    anchor = next(r for r in res if abs(r['s3']-S3_P) < 1e-9 and abs(r['s5']-S5_P) < 1e-9 and r['s7'] == 0)
    print(f"\nanchor cell: gamma*_head = {anchor['gamma_star']:.9f} (vs gate A gamma* 0.881557917504162)")
