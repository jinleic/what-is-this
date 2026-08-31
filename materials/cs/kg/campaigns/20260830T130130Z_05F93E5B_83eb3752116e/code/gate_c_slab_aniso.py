import sys, json, time, collections
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/kg/src')
from flint import arb
from gate_c_bnb import box_bound, GAMMA
from gate_c_fast import load_grid_entries, build_Q_polys

ents = load_grid_entries(); Qs, a10, a01 = build_Q_polys(ents)

def split_axis(box):
    """Bisect the axis contributing most to the V-interval width (hi^2 - lo^2), since the
    b1 slack b1(V_min) - b1(V) over the box is what blocks certified pruning."""
    widths = [(hi*hi - lo*lo, i) for i, (lo, hi) in enumerate(box)]
    _, i = max(widths)
    lo, hi = box[i]
    mid = 0.5*(lo+hi)
    A = list(box); B = list(box)
    A[i] = (lo, mid); B[i] = (mid, hi)
    return tuple(A), tuple(B)

def run(root, ref, budget=200000, min_width=0.0):
    gref = GAMMA[ref]
    q = collections.deque([root])
    n = npr = 0
    pruned_vol = 0.0
    surv = []
    t0 = time.time()
    while q:
        if n >= budget:
            surv.extend(list(q)); break
        box = q.popleft()
        ub, info = box_bound(box, Qs, a10, a01)
        n += 1
        vol = 1.0
        for (lo,hi) in box: vol *= (hi-lo)
        if ub <= gref:
            npr += 1; pruned_vol += vol; continue
        # stop splitting if all axes are already tiny
        if max(hi-lo for (lo,hi) in box) <= min_width:
            surv.append(box); continue
        A,B = split_axis(box)
        q.append(A); q.append(B)
    return dict(n=n, npr=npr, pruned_vol=pruned_vol, surv=surv, secs=time.time()-t0)

R = ((0.325, 0.350), (0.0, 0.080), (0.0, 0.020))
volR = (0.35-0.325)*0.08*0.02
out = {}
for ref in ("gateA","qhead"):
    r = run(R, ref, budget=200000, min_width=1e-6)
    sv = sum((b[0][1]-b[0][0])*(b[1][1]-b[1][0])*(b[2][1]-b[2][0]) for b in r["surv"])
    s7max = max((b[2][1] for b in r["surv"]), default=0.0)
    s3h = (min((b[0][0] for b in r["surv"]), default=0), max((b[0][1] for b in r["surv"]), default=0))
    s5h = (min((b[1][0] for b in r["surv"]), default=0), max((b[1][1] for b in r["surv"]), default=0))
    print(f"anisotropic slab ref {ref}: evals {r['n']} pruned {r['npr']} "
          f"pruned_vol {r['pruned_vol']/volR*100:.4f}% of R  surv {len(r['surv'])} "
          f"surv_vol {sv:.3e} ({sv/volR*100:.4f}%)  s7max {s7max:.6f} "
          f"s3 [{s3h[0]:.6f},{s3h[1]:.6f}] s5 [{s5h[0]:.6f},{s5h[1]:.6f}] [{r['secs']:.0f}s]")
    out[ref] = dict(n_eval=r['n'], n_pruned=r['npr'], pruned_frac_of_R=r['pruned_vol']/volR,
                    n_surv=len(r['surv']), surv_vol=sv, surv_frac_of_R=sv/volR,
                    s7max=s7max, s3hull=list(s3h), s5hull=list(s5h), seconds=r['secs'])
json.dump(out, open('/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/gc_bnb_slab_aniso.json','w'), indent=1)
print("wrote gc_bnb_slab_aniso.json")
