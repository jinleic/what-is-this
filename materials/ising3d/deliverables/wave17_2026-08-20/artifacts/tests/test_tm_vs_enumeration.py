
import numpy as np, time
from ising.lattices import hyperrect
from ising.exact_enumeration import dos_bonds
from ising.transfer_matrix import box_broken_bond_poly, torus_broken_bond_poly

def enum_broken(shape, periodic):
    lat = hyperrect(shape, periodic)
    g = dos_bonds(lat)          # index = satisfied bonds
    nb = lat.n_bonds
    return [int(g[nb-q]) for q in range(nb+1)], nb

ok = 0
# --- tori ---
for shape in [(4,4),(5,4),(2,2,2),(2,2,3),(2,3,3),(3,3,2),(2,2,4),(3,2,4),(2,4,3)]:
    per = True
    e,nb = enum_broken(shape, per)
    t = torus_broken_bond_poly(shape)
    e_tr = e[:len(t)] + [0]*max(0,len(t)-len(e))
    assert list(t) == [x for x in e][:len(t)] or t == e[:len(t)], (shape, t[:8], e[:8])
    while len(e)>1 and e[-1]==0: e.pop()
    assert t == e, (shape, t[:10], e[:10], len(t), len(e))
    ok += 1
    print("torus", shape, "OK  nb=",nb, " first coeffs", t[:6])

# --- open boxes ---
for shape in [(3,3),(4,4),(3,3,2),(3,3,3),(4,3,2)]:
    per = False
    e,nb = enum_broken(shape, per)
    t = box_broken_bond_poly(shape)
    while len(e)>1 and e[-1]==0: e.pop()
    assert t == e, (shape, t[:10], e[:10])
    ok += 1
    print("open box", shape, "OK nb=",nb, "first", t[:6])

# --- mixed: periodic cross-section, open transfer direction (cylinder) ---
for shape,per in [((3,3,2),(True,True)),((3,2,3),(True,True)),((4,3),(True,))]:
    dim=len(shape)
    fullper = tuple(list(per)+[False])
    e,nb = enum_broken(shape, fullper)
    t = box_broken_bond_poly(shape, periodic=per)
    while len(e)>1 and e[-1]==0: e.pop()
    assert t == e, (shape, t[:10], e[:10])
    ok += 1
    print("cylinder", shape, per, "OK nb=", nb)
print("ALL", ok, "cross-checks passed")
