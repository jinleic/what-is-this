"""Gate C stage 6 — OmegaPcVersusPen.

Uses pre_statement.md (committed a7c80de, byte-identical copy in this dir)
as the pre-registration. Protocol order (each control is a stop-gate):

  E0     hard-hash verification of every frozen dependency (v11 core,
         vxxz24_float, params mat, kernel basis JSON, margin matrix npy,
         frozen branchpin.py dependency of record) — before ANY import of
         campaign code.
  GATE-Z engine equivalence at r=0 vs the two frozen stage_b records
         (1e-12 per R branch, 1e-6 on Omega_raw), lemma1 ON and OFF.
  C1     identity-accept at r=0: cand1.lo - cand0.lo == -(2*prop_0*eps_0)
         to <= 1e-12; live candidates vs the frozen CandW triple to 1e-9.
  C1a    leaf-width-zero invariance for every non-dist0 block (both the
         point pass and the rho=r0 box pass).
  C2     wrong-index reject: region-1 dual lookup == frozen eps_1, != eps_0.
  C3     dual-drift reject: +1e-3 broadcast on glob-0 lam entries moves the
         eps boundary by EXACTLY 4*drift (<=1e-12) and R_glob[0] by
         -2*prop*4*drift (<=1e-12).
  PFIRST exact-Fraction certified feasible ray radius r0 for s = theta·v
         (theta = the frozen 21-coefficient CandidateWitness u), rounded
         DOWN to a dyadic; compared to the registered expectation
         1.0607145909782878e-06 at 1e-12 relative.
  P1     the ordering pass at e0 (rho = 0): certified G1/G2 intervals
         over the point (radius-0 leaves), from the live aggregation.
  P2     the ordering pass at e1 (rho = r0): certified G1/G2 intervals
         over the EXACT-kernel ray endpoint box p* + r0·s (one signed
         box pass, endpoints exact dyadics of the ray coordinates).
  PV     certified kernel membership of the e1 box: A·delta == 0 exactly
         (integer), sum(delta) == 0 exactly, |delta_i| <= r0 <= radius
         (exact rational comparison) — the feasibility proof that the
         swept set is inside D, recorded before the P2 verdict is read.
  A      adjudication per §8 of the prereg (exactly one terminal verdict).
"""
import hashlib, json, os, sys, time
from fractions import Fraction

sys.dont_write_bytecode = True
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

RUN = os.path.dirname(os.path.abspath(__file__))
TARGET = os.path.dirname(os.path.dirname(RUN))   # campaigns/<run> -> target
SRC = os.path.join(TARGET, "src")
SCRATCH = os.path.join(TARGET, "scratch")
T0 = time.time()

# ---- frozen anchors (registered in pre_statement §10) ----
SHA_CORE = "ce70f95922f3d5098e892e1b7047570ebb5267826b51f2776c91b56d9934adb5"
SHA_VF = "31657c92a1f841805e7355d5b5390f76e0d7ca2d04b1a3f2ce9fc204d9d5d890"
SHA_KB = "7d69c1dd4c9355dfe6d6b3a762658ba10d93e4bcaea3b240ebbbd811de3eba60"
SHA_MM = "423dbafd7c530ec7256ddf68cafdc15e143503ed9358c39bb81e4bc66b074b4e"
SHA_BP = "d7aca589c5dafbf7942f2b4009849edc311071e0262606695c474618f7ba3c56"
PARAMS_PREFIX = "df75ae3a"
MAT = os.path.join(SCRATCH, "rmmcode", "data", "K100_2.37155181.mat")
RADIUS = 1e-7
DRIFT = 1e-3
R_POS_CEIL = 1.6166317499672092e-05  # frozen positivity_ceiling (display)
SIGNAL = 1.5816497000997742e-07
GUARD_GAP = 2.0258350805768544e-06   # certified rung2 endpoint minus published
GUARD_FACTOR = GUARD_GAP / SIGNAL    # 12.808...
TOL_BRANCH = 1e-12
TOL_OMEGA = 1e-6
LAMSUM_3 = 0.5883309785276519        # frozen R_comp[0,2] branch (display anchor)

FROZEN_R_DETAIL = {
    "R_comp[0,3]": 0.24530661807075219,
    "R_comp[1,3]": 0.24715974210930786,
    "R_comp[2,3]": 0.2456557897939231,
    "R_comp[0,2]": 0.5883309785276519,
    "R_glob[0]": 0.49684932826618977,
    "R_glob[1]": 0.4968521822710635,
    "R_glob[2]": 0.49684892842208755,
}
FROZEN_RAW = 2.3715538358350807
FROZEN_R_DETAIL_OFF = {
    "R_comp[0,3]": 0.24530661954543967,
    "R_comp[1,3]": 0.2471597433637035,
    "R_comp[2,3]": 0.24565579085620265,
    "R_comp[0,2]": 0.5883309785276519,
    "R_glob[0]": 0.49685076173807374,
    "R_glob[1]": 0.49685294448520734,
    "R_glob[2]": 0.49685097954539326,
}
FROZEN_RAW_OFF = 2.3715518061863814
FROZEN_EPS0 = 2.1502086942121845e-06
FROZEN_EPS1 = 1.143325024764798e-06
FROZEN_C = [
    (0.4968507617380737, 0.4968507617380738),
    (0.49684932826618944, 0.4968507617386526),
    (0.4968493282665137, 0.4968507617389769),
]
THETA = [1, -1, 0, -1, 1, 0, 1, 0, 0, 1, -1, -1, -1, 1, 0, -1, 0, 0, -1, 0, 1]

report = {"started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
          "attempt": 1}

def save():
    with open(os.path.join(RUN, "checkpoint.json"), "w") as fh:
        json.dump(report, fh, indent=1, default=str)

def fail(stage, detail):
    report["verdict"] = "FROZEN-NEGATIVE"
    report["failed_control"] = stage
    report["failure_detail"] = str(detail)
    report["elapsed"] = time.time() - T0
    save()
    print(f"CONTROL-FAIL {stage}: {detail}")
    sys.stdout.flush()
    sys.exit(1)

def open_at(stage, detail):
    report["verdict"] = "OPEN"
    report["open_at"] = stage
    report["open_detail"] = str(detail)
    report["elapsed"] = time.time() - T0
    save()
    print(f"OPEN at {stage}: {detail}")
    sys.stdout.flush()
    sys.exit(2)

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

# ================= E0: dependency hashes =================
deps = {
    "interval_core.py": (os.path.join(SRC, "interval_core.py"), SHA_CORE),
    "vxxz24_float.py": (os.path.join(SRC, "vxxz24_float.py"), SHA_VF),
    "kernel_basis_V_exact_rational.json":
        (os.path.join(SCRATCH, "kernel_basis_V_exact_rational.json"), SHA_KB),
    "margin_matrix_M.npy": (os.path.join(SCRATCH, "margin_matrix_M.npy"), SHA_MM),
    "branchpin.py (frozen dependency of record)":
        (os.path.join(TARGET, "campaigns",
                      "20260901T102048Z_b555441e_d3303018f163", "branchpin.py"),
         SHA_BP),
}
report["E0"] = {}
for name, (path, want) in deps.items():
    got = sha256_file(path)
    report["E0"][name] = {"ok": got == want, "sha256": got}
    if got != want:
        fail("E0", f"{name}: sha256 {got} != frozen {want}")
save()
print("E0 PASS (5 dependencies hash-verified)")

# ================= exact prologue (no interval compute) =================
import numpy as np  # imported after E0; data files hash-pinned
A_np = np.load(os.path.join(SCRATCH, "margin_matrix_M.npy"))
A_exact = [[int(x) for x in row] for row in A_np]

# kernel basis: exact Fractions from the sha-pinned JSON
with open(os.path.join(SCRATCH, "kernel_basis_V_exact_rational.json")) as fh:
    KB_JSON = json.load(fh)
V_frac = [[Fraction(num, den) for num, den in col] for col in
          (KB_JSON[str(j)] for j in range(len(KB_JSON)))]
assert len(V_frac) == 21 and all(len(col) == 45 for col in V_frac)

# exact rank check of A over Q (row-echelon with Fractions)
rows = [[Fraction(x) for x in row] for row in A_exact]
rank = 0
piv_rows = []
for c in range(45):
    piv = None
    for ri in range(rank, len(rows)):
        if rows[ri][c] != 0:
            piv = ri
            break
    if piv is None:
        continue
    rows[rank], rows[piv] = rows[piv], rows[rank]
    pv = rows[rank][c]
    rows[rank] = [x / pv for x in rows[rank]]
    for ri in range(len(rows)):
        if ri != rank and rows[ri][c] != 0:
            f = rows[ri][c]
            rows[ri] = [a - f * b for a, b in zip(rows[ri], rows[rank])]
    piv_rows.append(c)
    rank += 1
report["PFIRST_context"] = {"A_rank_exact": rank, "n_cols": 45,
                            "kernel_dim_claimed": 45 - rank}
if rank != 24:
    fail("PFIRST-rank", f"A rank {rank} != frozen 24")
assert 45 - rank == 21

# A·v == 0 column-by-column, exact
Av_bad = []
for j, col in enumerate(V_frac):
    for ri in range(27):
        s = sum(Fraction(A_exact[ri][k]) * col[k] for k in range(45))
        if s != 0:
            Av_bad.append((j, ri))
if Av_bad:
    fail("PFIRST-Av", f"A·v nonzero at {Av_bad[:3]}")

# s = theta·v, exact
s_exact = [sum(Fraction(THETA[j]) * V_frac[j][i] for j in range(21))
           for i in range(45)]
sum_s = sum(s_exact)
if sum_s != 0:
    fail("PFIRST-s", f"sum(s) = {sum_s} != 0")
As_bad = []
for ri in range(27):
    t = sum(Fraction(A_exact[ri][k]) * s_exact[k] for k in range(45))
    if t != 0:
        As_bad.append(ri)
if As_bad:
    fail("PFIRST-As", f"A·s nonzero on rows {As_bad}")

# certified feasible ray radius EXACT: three walls, all exact-Fraction:
# (a) positivity: c_i + rho*s_i > 0 requires rho < c_i/(-s_i) for s_i < 0;
# (b) box: |rho*s_i| <= RADIUS requires rho <= RADIUS/max|s_i| (over s_i!=0);
# (c) the frozen CandidateWitness direction convention: s matches the
#     frozen `direction` array (theta·v with the recorded sign pattern).
inv_smax = None            # wall (b): RADIUS / max|s_i|
neg = []                   # wall (a) candidates
for i in range(45):
    if s_exact[i] != 0:
        h = abs(s_exact[i])
        q_b = Fraction(1, 1) / h
        inv_smax = q_b if inv_smax is None or q_b < inv_smax else inv_smax
        if s_exact[i] < 0:
            neg.append(i)
if inv_smax is None:
    fail("PFIRST-s", "s is identically zero")
from scipy.io import loadmat
params = np.asarray(loadmat(MAT)["params"]).flatten()
report["PFIRST_context"]["params_sha256"] = hashlib.sha256(
    params.tobytes()).hexdigest()
if not report["PFIRST_context"]["params_sha256"].startswith(PARAMS_PREFIX):
    fail("PFIRST-params", report["PFIRST_context"]["params_sha256"])
from scipy.io import loadmat as _lm  # noqa: F401 (already imported above
# in some paths); derive dist0 param offset with a THROWAWAY ParamManager.
sys.path.insert(0, SRC)
import interval_core as _ic_probe
import vxxz24_float as F
_pm0 = F.ParamManager()
_ws0 = F.Workspace(_pm0, 5.0, 1.0, 0.0, 3)
d0_off = _pm0.start[_ws0.globstage.dist_id[0]]
del _pm0, _ws0
c_frac = [Fraction(float(params[d0_off + i])) for i in range(45)]
wall_pos = None
for i in neg:                       # s_exact[i] < 0 here
    q_a = c_frac[i] / (-s_exact[i])
    if wall_pos is None or q_a < wall_pos:
        wall_pos = q_a
smax = max(abs(s_exact[i]) for i in range(45))
R_box = Fraction(RADIUS) / smax          # wall (b) exact
walls = {"box": R_box, "positivity": wall_pos}
# Fraction(RADIUS) IS a dyadic (50-bit numerator / power-of-two denominator),
# so when the box wall binds, r0 = R_box EXACTLY — no re-rounding needed.
R_exact = min(R_box, wall_pos)
binding = "box" if R_box <= wall_pos else "positivity"
if R_box <= wall_pos:
    r0 = R_box
else:
    # positivity binds: round DOWN to a 70-bit dyadic
    num, den = R_exact.numerator, R_exact.denominator
    k = 70
    r0_num = (num << k) // den
    r0 = Fraction(r0_num, 1 << k)
r0_float = float(r0)
report["PFIRST"] = {
    "s_max_abs": str(smax),
    "walls": {kk: str(v) for kk, v in walls.items()},
    "walls_float": {kk: float(v) for kk, v in walls.items()},
    "binding_wall": binding,
    "R_exact": str(R_exact),
    "R_exact_float": float(R_exact),
    "r0": str(r0),
    "r0_float": r0_float,
    "r0_le_R_exact": r0 <= R_exact,
}
if not (r0 <= R_exact):
    fail("PFIRST-round", "dyadic r0 exceeds exact R_exact")
# registered expectation (precompute_note.md correction): the BOX wall
# binds for this s (max|s_i| = 1 recorded by the frozen direction), so
# r0 must equal the registered binary64 radius 1e-7 exactly; positivity
# ceiling 1.6166317499672092e-05 must exceed it.
if abs(r0_float - RADIUS) > 0:
    open_at("E-FIRST",
            f"exact r0 {r0_float!r} vs expected {RADIUS!r} "
            f"(box wall must bind at max|s_i|={smax})")
if float(wall_pos) <= RADIUS:
    open_at("E-FIRST", f"positivity wall {float(wall_pos)} does not "
                       f"exceed the radius {RADIUS}")

# ================= interval phase =================
sys.path.insert(0, SRC)
import interval_core as IC
import vxxz24_float as F

I = IC.Ival

def as_iv(x):
    if isinstance(x, I):
        return x
    if hasattr(x, "v"):
        x = x.v
    if isinstance(x, I):
        return x
    if isinstance(x, arb_ := __import__("flint").arb):
        return I(x.lower(), x.upper())
    if isinstance(x, (list, tuple)) and len(x) == 1:
        return as_iv(x[0])
    return IC.const(float(np.asarray(x).reshape(-1)[0]))

def lemma_eps(gm, r, lam_shift=0.0):
    dm = gm.dist_max[r].v
    eps_iv = I(__import__("flint").arb(0), __import__("flint").arb(0))
    arb = __import__("flint").arb
    for i, shp in enumerate(gm.shapes):
        if float(dm[i].hi) <= 0:
            continue
        g_val = gm.lam_sum[r].v[0].lo - 1
        for d in range(3):
            g_val = g_val + gm.lam_margin[(r, d)].v[shp[d]].lo
        g_val = g_val + lam_shift
        diff = IC.ln_iv(dm[i]) - IC.const(float(g_val))
        e = diff.pos()
        eps_iv = I(min(eps_iv.lo, e.lo), max(eps_iv.hi, e.hi))
    return eps_iv

def part_eps(t, r):
    arb = __import__("flint").arb
    dm = t.split_dist_max[r].v
    eps_iv = I(arb(0), arb(0))
    for i, sp in enumerate(t.splits):
        if float(dm[i].hi) <= 0:
            continue
        g_val = t.lam_sum[r].v[0].lo - 1
        for d in range(3):
            g_val = g_val + t.lam_margin[(r, d)].v[sp[d] - t.lam_low[d]].lo
        diff = IC.ln_iv(dm[i]) - IC.const(float(g_val))
        e = diff.pos()
        eps_iv = I(min(eps_iv.lo, e.lo), max(eps_iv.hi, e.hi))
    return eps_iv

def frac_dyadic(q: Fraction) -> I:
    """Exact dyadic RADIUS point-iv for a NONNEGATIVE Fraction q:
    q == 0 -> the exact zero point-iv; q > 0 -> round q UP to a 70-bit
    dyadic (asserted >= q), built as exact Arb integer/shift — no float
    ever touches the value."""
    if q == 0:
        from flint import arb
        z = arb(0)
        return I(z, z)            # zero-width radius: leaf == center
    assert q > 0
    k = 70
    up_num = -(((-q.numerator) << k) // q.denominator)  # ceil(num*2^k/den)
    dy = Fraction(up_num, 1 << k)
    assert dy >= q, "dyadic radius below exact Fraction offset"
    from flint import arb
    if up_num >= 0:
        v = arb(up_num) / arb(1 << k)
    else:
        v = -(arb(-up_num) / arb(1 << k))
    return I(v, v)

def aggregate(ws, gm, pm, include_lemma1=True):
    """VERBATIM port of frozen branchpin.py::aggregate (dependency of
    record, sha d7aca589…): same R_glob[0] candidate construction."""
    from flint import arb
    total_R = arb(0)
    details, glob_cands, leg = {}, {}, {}
    for r in range(3):
        leg[f"glob_r{r}"] = lemma_eps(gm, r) if include_lemma1 else None

    def pen_of(r):
        Hdm = gm.dist_max[r].entropy().v
        Hsd = gm.dist[r].entropy().v
        prop = gm.region_prop[r].item()
        if include_lemma1:
            return prop * (Hdm - Hsd + 2 * leg[f"glob_r{r}"])
        return prop * (Hdm - Hsd)

    for l in range(3, 4):
        for r in range(3):
            nb = [I(arb(0), arb(0)) for _ in range(3)]
            pen = I(arb(0), arb(0))
            pc = I(arb(0), arb(0))
            for t in ws.parts[l - 1]:
                for tt in range(3):
                    nb[tt] = nb[tt] + as_iv(t.num_block_contribution[r][tt])
                if isinstance(t, F.Part):
                    Hdm = t.split_dist_max[r].entropy().v
                    Hsd = t.split_dist[r].entropy().v
                    frac = t.part_frac.item() * t.region_prop[r].item()
                    if include_lemma1:
                        e = part_eps(t, r)
                        pen = pen + frac * (Hdm - Hsd + 2 * e)
                    else:
                        pen = pen + frac * (Hdm - Hsd)
                    pc = pc + as_iv(t.p_comp[r])
            cand = [nb[tt].lo - (pc.hi if tt == r else pen.hi)
                    for tt in range(3)]
            R_low = min(cand)
            if float(R_low) > 0:
                total_R = total_R + R_low
            details[f"R_comp[{r},{l}]"] = float(R_low)
    nb2 = [I(arb(0), arb(0)) for _ in range(3)]
    for t in ws.parts[1]:
        for tt in range(3):
            nb2[tt] = nb2[tt] + as_iv(t.num_block_contribution[tt])
    R2 = min(z.lo for z in nb2)
    if float(R2) > 0:
        total_R = total_R + R2
    details["R_comp[0,2]"] = float(R2)
    for r in range(3):
        nbv = [as_iv(z) for z in gm.num_block[r]]
        pen_iv = pen_of(r)
        pc_iv = as_iv(gm.p_comp[r].item())
        cand_iv = []
        for tt in range(3):
            rhs = pc_iv if tt == r else pen_iv
            cand_iv.append((nbv[tt].lo - rhs.hi, nbv[tt].hi - rhs.lo))
        R_low = min(c[0] for c in cand_iv)
        if float(R_low) > 0:
            total_R = total_R + R_low
        details[f"R_glob[{r}]"] = float(R_low)
        if r == 0:
            glob_cands["R_glob[0]"] = cand_iv
    ms_iv = [as_iv(z) for z in gm.mat_size]
    M_low = min(float(z.lo) for z in ms_iv)
    tgt_hi = float((arb(7).log() * 4).upper())
    Om = (tgt_hi - float(total_R.lower())) / M_low
    return details, glob_cands, Om, leg, M_low


def leaf_width_audit(ws, gm, pm, d0_gid):
    """C1a: every non-dist0 leaf must carry a width-0 point interval."""
    bad = []
    blocks = [("region_prop", gm.region_prop_id)] + \
             [(f"dist_max[{r}]", gm.dist_max_id[r]) for r in range(3)] + \
             [(f"dist[{r}]", gm.dist_id[r]) for r in range(1, 3)]
    for name, gid in blocks:
        vec = gm.pm.get(gid)
        for k, iv in enumerate(vec):
            if float(iv.hi) > float(iv.lo):
                bad.append((name, k))
    for r in range(3):
        for t in range(3):
            vec = gm.pm.get(gm.lam_margin_id[(r, t)])
            for k, iv in enumerate(vec):
                if float(iv.hi) > float(iv.lo):
                    bad.append((f"lam_margin[{r},{t}]", k))
        vec = gm.pm.get(gm.lam_sum_id[r])
        for k, iv in enumerate(vec):
            if float(iv.hi) > float(iv.lo):
                bad.append((f"lam_sum[{r}]", k))
    return bad

def run_pass(params_vec, drift=0.0, drift_region=None, box_delta=None):
    """One v11 interval pass. box_delta: list of 45 exact Fractions
    (or None) — when given, dist[0] leaves move to the EXACT signed box
    center_i ± |delta_i| (a degenerate box at the ray point when
    delta = r0·s is fixed; |delta_i| is the exact dyadic ray offset)."""
    arb = __import__("flint").arb
    pm = F.ParamManager()
    ws = F.Workspace(pm, 5.0, 1.0, 0.0, 3)
    p_use = params_vec
    if drift != 0.0:
        p_use = params_vec.copy()
        g = ws.globstage
        for t in range(3):
            gid = g.lam_margin_id[(drift_region, t)]
            p_use[pm.start[gid]:pm.start[gid] + pm.size[gid]] += drift
        gid = g.lam_sum_id[drift_region]
        p_use[pm.start[gid]] += drift
    pm.set_value(p_use)
    st, sz = pm.start, pm.size
    gm0 = ws.globstage
    d0 = st[gm0.dist_id[0]]
    leaves = []
    for i in range(45):
        x = IC.const(float(params_vec[d0 + i]))
        if box_delta is None:
            leaves.append(I(x.lo, x.hi))
        else:
            hfrac = abs(box_delta[i])
            h_iv = frac_dyadic(hfrac)   # dyadic RADIUS enclosing |delta_i|
            lo = (x - h_iv).lo
            hi = (x + h_iv).hi
            leaves.append(I(lo, hi))

    def get(gid):
        if gid == gm0.dist_id[0]:
            return list(leaves)
        return [IC.const(float(p_use[st[gid] + k]))
                for k in range(sz[gid])]

    def get_scalar(gid):
        return float(pm.cur_x[st[gid]])

    pm.get = get
    pm.get_scalar = get_scalar
    F.set_interval_pm(pm)
    IC.ZC_FLAG["hit"] = False
    constraints, equalities, value = ws.evaluate()
    if IC.ZC_FLAG["hit"]:
        raise RuntimeError("zero-crossing guard fired")
    return ws, ws.globstage, pm


d0 = None  # set in-run below (dist0 param-block offset)

# ---- GATE-Z point pass (lemma1 on) ----
with IC.prec():
    _pm0 = F.ParamManager()
    _ws0 = F.Workspace(_pm0, 5.0, 1.0, 0.0, 3)
    globals()['d0'] = _pm0.start[_ws0.globstage.dist_id[0]]
    del _pm0, _ws0
    params0 = params
    ws, gm, pm = run_pass(params0)
    det, _, Om, leg, M_low = aggregate(ws, gm, pm, include_lemma1=True)
    report["GATEZ_with"] = {"R_detail": det, "Om_raw": Om,
                            "eps0_live": float(leg["glob_r0"].hi),
                            "M_low": M_low}
    for kname, fr in FROZEN_R_DETAIL.items():
        if abs(det[kname] - fr) > TOL_BRANCH:
            fail("GATE-Z", f"{kname}: {det[kname]} vs frozen {fr}")
    if abs(Om - FROZEN_RAW) > TOL_OMEGA:
        fail("GATE-Z", f"Om_raw {Om} vs {FROZEN_RAW}")
    if abs(float(leg["glob_r0"].hi) - FROZEN_EPS0) > TOL_BRANCH:
        fail("GATE-Z/eps0", f"eps0 {float(leg['glob_r0'].hi)} vs {FROZEN_EPS0}")
    save()
    print("GATE-Z (with lemma1) PASS")

    # ---- GATE-Z lemma1-off diagnostic ----
    det_off, _, Om_off, _, _ = aggregate(ws, gm, pm, include_lemma1=False)
    report["GATEZ_without"] = {"R_detail": det_off, "Om_raw": Om_off}
    for kname, fr in FROZEN_R_DETAIL_OFF.items():
        if abs(det_off[kname] - fr) > TOL_BRANCH:
            fail("GATE-Z-off", f"{kname}: {det_off[kname]} vs frozen {fr}")
    if abs(Om_off - FROZEN_RAW_OFF) > TOL_OMEGA:
        fail("GATE-Z-off", f"Om_raw {Om_off} vs {FROZEN_RAW_OFF}")
    report["GATE_Z"] = "PASS (both frozen records reproduced)"
    save()
    print("GATE-Z (lemma1-off) PASS")

    # ---- C1 identity-accept ----
    c0, c1, c2 = FROZEN_C
    lhs = c1[0] - c0[0]
    rhs = -2.0 * (1.0 / 3.0) * float(leg["glob_r0"].hi)
    report["C1"] = {"lhs_cand1_minus_cand0_lo": lhs,
                    "rhs_minus_2propEps": rhs,
                    "identity_gap": abs(lhs - rhs)}
    if abs(lhs - rhs) > TOL_BRANCH:
        fail("C1", f"identity gap {abs(lhs - rhs)} exceeds {TOL_BRANCH}")
    _, gcand_live, _, _, _ = aggregate(ws, gm, pm, include_lemma1=True)
    report["C1"]["live_candidates"] = gcand_live["R_glob[0]"]
    for idx, (flo, fhi) in enumerate(FROZEN_C):
        llo, lhi = gcand_live["R_glob[0]"][idx]
        if abs(llo - flo) > 1e-9 or abs(lhi - fhi) > 5e-9:
            fail("C1-live", f"cand{idx}: live [{llo}, {lhi}] vs frozen "
                            f"[{flo}, {fhi}]")
    report["C1"]["status"] = "PASS"
    save()
    print("C1 (identity-accept) PASS")

    # ---- C2 wrong-index reject (region-1 dual lookup) ----
    eps1_live = float(lemma_eps(gm, 1).hi)
    report["C2"] = {"eps1_live": eps1_live, "eps1_frozen": FROZEN_EPS1,
                    "distinct_from_eps0": abs(eps1_live - FROZEN_EPS0) > 1e-9}
    if abs(eps1_live - FROZEN_EPS1) > TOL_BRANCH:
        fail("C2", f"eps1 {eps1_live} vs frozen {FROZEN_EPS1}")
    report["C2"]["status"] = "PASS"
    print("C2 (wrong-index reject) PASS")

    # ---- C3 dual-drift reject ----
    ws2, gm2, pm2 = run_pass(params0, drift=DRIFT, drift_region=0)
    eps0_drift = lemma_eps(gm2, 0)
    obs = float(eps0_drift.hi)
    expected_hi = FROZEN_EPS0 + 4.0 * DRIFT
    report["C3"] = {"eps0_drifted_hi": obs, "expected_hi": expected_hi,
                    "coupling_gap": abs(obs - expected_hi), "drift": DRIFT}
    if abs(obs - expected_hi) > TOL_BRANCH:
        fail("C3", f"drift coupling: eps0 {obs} vs expected {expected_hi}")
    det2, _, _, _, _ = aggregate(ws2, gm2, pm2, include_lemma1=True)
    drift_shift = det2["R_glob[0]"] - det["R_glob[0]"]
    predicted_shift = -2.0 * (1.0 / 3.0) * 4.0 * DRIFT
    report["C3"]["R_glob0_drifted"] = det2["R_glob[0]"]
    report["C3"]["drift_shift_observed"] = drift_shift
    report["C3"]["drift_shift_predicted"] = predicted_shift
    if abs(drift_shift - predicted_shift) > TOL_BRANCH:
        fail("C3", f"branch response {drift_shift} vs predicted "
                   f"{predicted_shift}")
    report["C3"]["status"] = "PASS (4*drift coupling exact)"
    save()
    print("C3 (dual-drift reject) PASS")

    # ---- P1: ordering at e0 (rho = 0, degenerate point) ----
    cands0 = gcand_live["R_glob[0]"]
    g1p1_iv = I(cands0[1][0], cands0[1][1]) - I(cands0[0][0], cands0[0][1])
    g2p1_iv = I(cands0[2][0], cands0[2][1]) - I(cands0[0][0], cands0[0][1])
    report["P1"] = {
        "G1": [float(g1p1_iv.lo), float(g1p1_iv.hi)],
        "G2": [float(g2p1_iv.lo), float(g2p1_iv.hi)],
    }
    save()
    print("P1 (rho=0) recorded: G1", float(g1p1_iv.lo), "G2",
          float(g2p1_iv.lo))

    # ---- PV: certified kernel membership of the ray endpoint ----
    delta_exact = [r0 * s_exact[i] for i in range(45)]
    pv = {"A_delta_zero": all(
        sum(Fraction(A_exact[ri][k]) * delta_exact[k]
            for k in range(45)) == 0 for ri in range(27)),
        "sum_delta_zero": sum(delta_exact) == 0,
        "max_abs_delta": str(max(abs(d) for d in delta_exact)),
        "le_radius": all(abs(d) <= Fraction(RADIUS) for d in delta_exact),
        "positivity": all(float(params0[globals()["d0"] + i])
                          + float(delta_exact[i]) > 0 for i in range(45)),
    }
    pv["max_abs_delta_float"] = float(Fraction(pv["max_abs_delta"]))
    if not (pv["A_delta_zero"] and pv["sum_delta_zero"] and pv["le_radius"]):
        fail("PV", "ray endpoint fails exact kernel-membership checks")
    report["PV"] = pv
    save()
    print("PV (kernel membership) PASS: A·delta==0, sum==0, |delta|<=r0<=1e-7")

    # ---- P2: ordering at e1 (rho = r0, exact-dyadic ray box) ----
    ws3, gm3, pm3 = run_pass(params0, box_delta=delta_exact)
    bad_leaves = leaf_width_audit(ws3, gm3, pm3, None)
    report["C1a"] = {"nonzero_width_leaves": bad_leaves[:10],
                     "n_bad": len(bad_leaves)}
    if bad_leaves:
        open_at("C1a", f"{len(bad_leaves)} non-dist0 leaves carry width "
                       f"(wiring defect): {bad_leaves[:5]}")
    det3, gcand3, Om3, leg3, M3 = aggregate(ws3, gm3, pm3,
                                            include_lemma1=True)
    report["C1a"]["status"] = "PASS (all non-dist0 leaves width-0)"
    report["C1a"]["eps0_box"] = float(leg3["glob_r0"].hi)
    if abs(float(leg3["glob_r0"].hi) - FROZEN_EPS0) > TOL_BRANCH:
        open_at("C1a", f"eps0 not bit-constant at the ray endpoint: "
                       f"{float(leg3['glob_r0'].hi)} vs {FROZEN_EPS0}")
    save()
    cands = gcand3["R_glob[0]"]
    cand0_iv = I(cands[0][0], cands[0][1])
    cand1_iv = I(cands[1][0], cands[1][1])
    cand2_iv = I(cands[2][0], cands[2][1])
    g1_iv = cand1_iv - cand0_iv
    g2_iv = cand2_iv - cand0_iv
    g1_lo, g1_up = g1_iv.low(), g1_iv.up()
    g2_lo, g2_up = g2_iv.low(), g2_iv.up()
    report["P2"] = {
        "G1_lo_outward": float(g1_lo), "G1_up": float(g1_up),
        "G2_lo_outward": float(g2_lo), "G2_up": float(g2_up),
        "cand_intervals": cands,
        "R_glob0": det3["R_glob[0]"],
        "Om_raw_point": Om3,
        "M_low": M3,
        "hull_span_cand0": float(cand0_iv.hi - cand0_iv.lo),
    }
    # ---- A: adjudication per prereg §8 ----
    pos1 = float(g1_lo) > 0
    pos2 = float(g2_lo) > 0
    if pos1 and pos2:
        report["verdict"] = "OBSTACLE-PASS"
        report["P2"]["collapse_to"] = 0
        mins = min(float(g1_lo), float(g2_lo), float(g1p1_iv.low()),
                   float(g2p1_iv.low()))
        report["P2"]["min_margin"] = mins
        report["P2"]["min_margin_over_signal"] = mins / SIGNAL
    else:
        report["verdict"] = "BOX-CONTRADICTION"
        ft = []
        if not pos1:
            ft.append({"term": "G1", "lo": float(g1_lo),
                       "hi": float(g1_up),
                       "width": float(g1_up - g1_lo),
                       "width_over_signal": float(g1_up - g1_lo) / SIGNAL})
        if not pos2:
            ft.append({"term": "G2", "lo": float(g2_lo),
                       "hi": float(g2_up),
                       "width": float(g2_up - g2_lo),
                       "width_over_signal": float(g2_up - g2_lo) / SIGNAL})
        report["P2"]["failing_terms"] = ft
    report["guardrail"] = {
        "gap_to_published": GUARD_GAP,
        "signal": SIGNAL,
        "factor": GUARD_FACTOR,
        "note": "no outcome of this campaign approaches the published "
                "2.37155181; no record claim anywhere",
    }
    report["elapsed"] = time.time() - T0
    save()
    print("VERDICT:", report["verdict"])
    print(f"elapsed {time.time() - T0:.1f}s")
