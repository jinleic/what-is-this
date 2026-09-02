"""Stage-(a) float64 transcription of the Vassilevska Williams-Xu-Xu-Zhou
2024 (SODA'24) released verifier `rmmcode` (osf.io/7wgh2), dual-mode so the
same code path runs under interval semantics in stage (b).

Transcribed files (scratch/rmmcode/src/, read first-hand 2026-08-30):
  evaluation/Workspace.m, GlobalStage.m, PartInfo.m, PartInfoLv2.m,
  PartInfoZero.m, FindOrCreatePart.m, FindPartByIdentifier.m;
  utils/{PrepareSplits,PrepareShapes,JointToMargin,MarginalDist,Rot3,Rot3c}.m;
  complete_split/{EncodeCSD,DecodeCSD,ConcatCSD}.m;
  autograd/GVar.m (value part only), verify/VerifyOmega.m, GetFeasibility.m.

Key structural deltas vs the Alman25 (rung-1) transcription
(src/alman25_float.py), per src/stage_b_rung2_notes.md and the MATLAB diff:
  * 3 regions (r = 1..3), no Dims/permutation table: region r's shared
    dimension IS r (GlobalStage: complete_split{r}{r}, dist_margins{r};
    PartInfo p_comp: shared dim = region).
  * ONE p_comp per region (Z-compatibility), not p_compY/p_compZ.
  * p_comp numerator: component-with-zero test is `min(ptr.shape) == 0`
    (PartInfo.m line 284) — there is no Y/Z special-case split.
  * Registration order (Workspace.Build): per r in 1..3: num_retain_comp
    groups for lv = 2..max_level, then num_retain_glob{r}; then
    single_mat_size; then omega (bounds 0..Inf; K registered omega-mode
    fixed = expinfo.K); then globstage.Build().
  * VXXZ24 extra linear constraints on the interval side's feasibility
    audit: global-stage region_prop Y/Z symmetry ([0,1,-1]) and (omega
    mode, K==1) X/Y symmetry ([1,-1,0]).
  * Schonhage line written explicitly: c_viol gets
    log(q+2) * 2^(max_level-1) - value.

Part naming: part  (VXXZ24 `part`) == term (Alman25 `term`).
Level-2 shapes {013,031,004} handled exactly as in PartInfoLv2.m (the
closed-form contribution structure is identical between the two releases).

The interval mode is exercised by stage_b_rung2.py; ParamManager.get is
monkey-patched there (IntervalTree pattern), exactly as in interval_core.py.
GVar is transcribed value-only: gradients are optimizer machinery, not part
of the verified constraint set.
"""

import numpy as np
from flint import arb
from scipy.io import loadmat

Q = 5.0
MAX_LEVEL = 3


def PrepareShapes(level):
    sc = 2 ** level
    return [(i, j, sc - i - j) for i in range(sc + 1) for j in range(sc - i + 1)]


def PrepareSplits(shape):
    """1-based 6-tuples (i,j,k,i2,j2,k2) like MATLAB rows (PrepareSplits.m)."""
    sc = sum(shape)
    sh = sc // 2
    out = []
    for i in range(min(sh, shape[0]) + 1):
        for j in range(min(sh - i, shape[1]) + 1):
            k = sh - i - j
            if k > shape[2]:
                continue
            out.append((i, j, k, shape[0] - i, shape[1] - j, shape[2] - k))
    return out


class JointToMargin:
    """JointToMargin.m: 3 sparse (n_col x (sum_col+1)) matrices."""

    def __init__(self, columns, sum_col):
        n = len(columns)
        self.mats = [np.zeros((n, sum_col + 1)) for _ in range(3)]
        for i, col in enumerate(columns):
            for t in range(3):
                self.mats[t][i, col[t]] = 1.0

    def apply(self, dist):
        # interval path: list of Ivals
        if isinstance(dist, list) and dist and dist[0].__class__.__name__ == "Ival":
            import interval_core as IC
            out = []
            for m in self.mats:
                v = []
                for j in range(m.shape[1]):
                    s = None
                    for i in np.nonzero(m[:, j])[0]:
                        s = dist[i] if s is None else s + dist[i]
                    if s is None:
                        s = IC.const(0.0)
                    v.append(s)
                out.append(v)
            return out
        return [dist @ m for m in self.mats]


def Rot3(x, n):
    """Rot3.m: rotate a 3-element vector LEFT by n."""
    n %= 3
    if n == 1:
        return (x[1], x[2], x[0])
    if n == 2:
        return (x[2], x[0], x[1])
    return x


def Rot3c(x, n):
    """Rot3c.m: rotate a list x RIGHT by n (VXXZ24 uses it on cell arrays of
    complete_split; identical semantics to Alman25's)."""
    n %= 3
    if n == 1:
        return (x[1], x[2], x[0])
    if n == 2:
        return (x[2], x[0], x[1])
    return x


def EncodeCSD(arr):
    """EncodeCSD.m; 0-based here (MATLAB returns 1-based = idx+1)."""
    idx = 0
    for a in arr:
        idx = idx * 3 + a
    return idx


def DecodeCSD(identity, power):
    """DecodeCSD.m inverse; identity 0-based here."""
    arr = [0] * power
    idx = identity
    for i in range(power - 1, -1, -1):
        arr[i] = idx % 3
        idx //= 3
    return arr


def kron_csd(lhs, rhs):
    """ConcatCSD.m value part."""
    if (isinstance(lhs, (list, np.ndarray)) and len(lhs) and
            getattr(lhs[0], '__class__', type).__name__ == 'Ival'):
        import interval_core as IC
        out = []
        for a in lhs:
            for b in rhs:
                out.append(a * IC.const(b) if not isinstance(b, IC.Ival) else a * b)
        return out
    return np.kron(lhs, rhs)


# ------------------------------------------------------------------ entropy
# GVar.Entropy value part: sum -x*log(x + (x<=0)) elementwise (GVar.m:201-205).

def entropy_vec(v):
    if isinstance(v, list) and v and v[0].__class__.__name__ == "Ival":
        import interval_core as IC
        R = IC.Ival(arb(0), arb(0))
        for wi in v:
            R = R + IC.ent_iv_point(wi)
        return R
    v = np.asarray(v, dtype=np.float64)
    shifted = np.where(v <= 0.0, v + 1.0, v)
    with np.errstate(divide="ignore", invalid="ignore"):
        logs = np.log(shifted)
    return float(np.sum(-v * logs))


def normalized_entropy(a, p):
    """GVar.NormalizedEntropy value part: sum -a_i ln(a_i/p) (GVar.m:213-224)."""
    a = np.asarray(a, dtype=np.float64)
    if p <= 0.0:
        return 0.0  # corner case at p==0: value 0 (all a_i==0)
    ratio = np.where(a <= 0.0, a / p + 1.0, a / p)
    with np.errstate(divide="ignore", invalid="ignore"):
        logs = np.log(ratio)
    return float(np.sum(-a * logs))


# -------------------------------------------------------------- param manager

class ParamManager:
    """ParamManager.m, value part only (no gradients, no SNOPT)."""

    def __init__(self):
        self.start = []      # 0-based start of each group
        self.size = []
        self.lb = []
        self.ub = []
        self.cur_x = np.zeros(0)
        # linear constraints (A x <= b, Aeq x = beq) recorded for the
        # stage-(b) feasibility audit; GetLinearConstraints drops empty rows.
        self.lin_A = []      # list of (dict gid->row_vector, b_row)
        self.lin_Aeq = []    # same

    def register(self, n, lb, ub):
        gid = len(self.start)
        self.start.append(sum(self.size) if self.size else 0)
        self.size.append(n)
        self.lb.extend([lb] * n if np.isscalar(lb) else lb)
        self.ub.extend([ub] * n if np.isscalar(ub) else ub)
        return gid

    @property
    def num_input(self):
        return sum(self.size)

    def set_value(self, x):
        x = np.asarray(x, dtype=np.float64)
        self.cur_x = np.clip(x, np.asarray(self.lb), np.asarray(self.ub))

    def get(self, gid):
        s = self.start[gid]
        n = self.size[gid]
        return self.cur_x[s:s + n].copy()

    def get_scalar(self, gid):
        return float(self.get(gid)[0])

    # -- linear constraint recording (AddLinearConstraintEq) ---------------
    def _lc_rows(self, entries, b):
        """AddLinearConstraint(Eq) row expansion. coeff is (n_group x r)."""
        nrow = 1 if np.isscalar(b) else len(np.atleast_1d(b))
        rows = [dict() for _ in range(nrow)]
        b = np.atleast_1d(b)
        for gid, coeff in entries:
            coeff = np.asarray(coeff)
            if coeff.ndim == 1:
                coeff = coeff.reshape(1, -1)
            if coeff.shape[1] == 1 and nrow > 1:
                coeff = np.tile(coeff, (1, nrow))
            for r in range(nrow):
                for c in range(self.size[gid]):
                    v = float(coeff[r, c]) if (r < coeff.shape[0] and
                                               c < coeff.shape[1]) else 0.0
                    if v != 0.0:
                        rows[r][self.start[gid] + c] = v
        return rows, b

    def add_lincon_eq(self, entries, b):
        rows, b = self._lc_rows(entries, b)
        self.lin_Aeq.extend(rows)
        self._lin_beq = b

    def add_lincon(self, entries, b):
        rows, b = self._lc_rows(entries, b)
        self.lin_A.extend(rows)
        self._lin_b = b

    def linear_violations(self):
        """max |Aeq x - beq| over rows (0 rows in A by construction here)."""
        out = 0.0
        for row, br in zip(self.lin_Aeq, np.atleast_1d(self._lin_beq)):
            s = sum(v * self.cur_x[c] for c, v in row.items()) - float(br)
            out = max(out, abs(s))
        return out


def zeros(n):
    # interval-mode coercion: when the active pm.get produces Ivals, yield a
    # point interval so downstream arithmetic stays in endpoint-pair land.
    if _PM_GLOBAL is not None:
        try:
            v0 = _PM_GLOBAL.get(0)[0]
            if v0.__class__.__name__ == "Ival":
                import interval_core as IC
                return G(IC.const(0.0))
        except Exception:
            pass
    return G(np.zeros(n))


_PM_GLOBAL = None


def _interval_mode_active():
    if _PM_GLOBAL is None:
        return False
    try:
        v0 = _PM_GLOBAL.get(0)[0]
        return v0.__class__.__name__ == "Ival"
    except Exception:
        return False


def set_interval_pm(pm):
    global _PM_GLOBAL
    _PM_GLOBAL = pm


class G:
    """Value holder emulating a GVar row vector (values only, dual-mode)."""
    __slots__ = ("v",)

    def __init__(self, v):
        if isinstance(v, arb):
            self.v = v
            return
        if v.__class__.__name__ == "Ival":
            self.v = v
            return
        if isinstance(v, list) and v and (isinstance(v[0], arb) or
                                          v[0].__class__.__name__ == "Ival"):
            self.v = list(v)
            return
        self.v = np.asarray(v, dtype=np.float64)

    @staticmethod
    def _comb(a, b, op):
        from interval_core import Ival, const
        is_iv = lambda z: z.__class__.__name__ == "Ival"
        is_arb = lambda z: isinstance(z, arb)
        as_list = lambda z: isinstance(z, list) or (
            type(z) is np.ndarray and z.dtype == object)
        if type(a) is np.ndarray and a.dtype != object:
            return op(a, b)
        if type(b) is np.ndarray and b.dtype != object:
            return op(a, b)
        if as_list(a) and as_list(b):
            return [op(x, y) for x, y in zip(a, b)]
        if as_list(a):
            bI = b if (is_iv(b) or is_arb(b)) else const(b)
            return [op(x, bI) for x in a]
        if as_list(b):
            aI = a if (is_iv(a) or is_arb(a)) else const(a)
            return [op(aI, y) for y in b]
        if is_iv(a) or is_iv(b):
            aI = a if is_iv(a) else const(a)
            bI = b if is_iv(b) else const(b)
            return op(aI, bI)
        if is_arb(a) or is_arb(b):
            aA = a if is_arb(a) else arb(float(a))
            bA = b if is_arb(b) else arb(float(b))
            return op(aA, bA)
        return op(a, b)

    def __add__(self, o):
        ov = o.v if isinstance(o, G) else o
        return G(self._comb(self.v, ov, lambda x, y: x + y))

    def __radd__(self, o):
        ov = o.v if isinstance(o, G) else o
        return G(self._comb(ov, self.v, lambda x, y: x + y))

    def __sub__(self, o):
        ov = o.v if isinstance(o, G) else o
        return G(self._comb(self.v, ov, lambda x, y: x - y))

    def __rsub__(self, o):
        ov = o.v if isinstance(o, G) else o
        return G(self._comb(ov, self.v, lambda x, y: x - y))

    def __mul__(self, o):
        ov = o.v if isinstance(o, G) else o
        return G(self._comb(self.v, ov, lambda x, y: x * y))

    def __rmul__(self, o):
        ov = o.v if isinstance(o, G) else o
        return G(self._comb(ov, self.v, lambda x, y: x * y))

    def __truediv__(self, o):
        ov = o.v if isinstance(o, G) else o
        return G(self._comb(self.v, ov, lambda x, y: x / y))

    def entropy(self):
        v = self.v
        if isinstance(v, list) and v and v[0].__class__.__name__ == "Ival":
            import interval_core as IC
            return G(IC.ent_vec(v))
        if isinstance(v, list):
            return G(ent_vec_iv(v))
        if isinstance(v, arb):
            return G(ent_iv(v))
        if v.__class__.__name__ == "Ival":
            import interval_core as IC
            return G(IC.ent_iv_point(v))
        return G(entropy_vec(self.v))

    def normalized_entropy(self, p):
        pv = p.v if isinstance(p, G) else p
        v = self.v
        if isinstance(v, list) and v and v[0].__class__.__name__ == "Ival":
            import interval_core as IC
            if pv.__class__.__name__ != "Ival":
                pv2 = IC.const(float(np.asarray(pv).reshape(-1)[0]))
            else:
                pv2 = pv
            return G(IC.nent_vec(v, pv2))
        if isinstance(v, list):
            return G(nent_vec_iv(v, pv))
        if isinstance(v, arb):
            return G(nent_iv(v, pv))
        if pv.__class__.__name__ == "Ival":
            # scalar interval x list-of-floats: wrap operands
            import interval_core as IC
            return G(IC.nent_vec([IC.const(z) for z in
                                  (self.v if isinstance(self.v, (list, np.ndarray))
                                   else list(np.atleast_1d(self.v)))], pv))
        return G(normalized_entropy(self.v, pv))

    def __matmul__(self, m):
        v = self.v
        if isinstance(v, list):
            out = []
            for j in range(len(m[0]) if len(m) else 0):
                s = None
                for i in range(len(m)):
                    if m[i][j] != 0:
                        s = v[i] if s is None else s + v[i]
                out.append(arb(0) if s is None else s)
            return G(out)
        return G(self.v @ m)

    def item(self):
        if self.v.__class__.__name__ == "Ival":
            return self.v
        if isinstance(self.v, arb):
            return self.v
        if isinstance(self.v, list):
            if len(self.v) == 1:
                inner = self.v[0]
                if inner.__class__.__name__ == "Ival":
                    return inner
                if isinstance(inner, arb):
                    return inner
                return float(inner)
            inner = self.v[0]
            return float(inner) if not isinstance(inner, arb) else inner
        return float(np.asarray(self.v).reshape(-1)[0])

    def __len__(self):
        return len(self.v)

    def __getitem__(self, k):
        v = self.v
        if isinstance(v, arb):
            return G(v)
        if isinstance(v, list):
            return G(v[k])
        return G(self.v[k])


def ent_iv(x: arb) -> arb:
    """Interval enclosure of -t*ln(t) for t in ball x (t >= 0)."""
    lo, hi = x.lower(), x.upper()
    if hi <= 0:
        return arb(0)
    if lo < 0:
        lo = arb(0)
    cands = []
    if lo > 0:
        cands.append(-lo * lo.log())
    if hi > 0:
        cands.append(-hi * hi.log())
    if not cands:
        return arb(0)
    lo_v = min(c.lower() for c in cands)
    hi_v = max(c.upper() for c in cands)
    return arb(lo_v, hi_v)


def ent_vec_iv(ws) -> arb:
    t = arb(0)
    for wi in ws:
        t = t + ent_iv(wi)
    return t


def nent_iv(x: arb, p: arb) -> arb:
    if x.upper() <= 0:
        return arb(0)
    return ent_iv(x) + x * p.log()


def nent_vec_iv(ws, p) -> arb:
    t = arb(0)
    for wi in ws:
        if wi.upper() <= 0:
            continue
        t = t + nent_iv(wi, p)
    return t


# ------------------------------------------------------------------- parts

class PartLv2:
    """PartInfoLv2.m — level 2 only."""

    def __init__(self, pm, q, level, part_id, shape, identifier):
        assert level == 2
        self.pm = pm
        self.q = q
        self.level = 2
        self.power = 2
        self.part_id = part_id
        self.shape = tuple(shape)
        self.identifier = tuple(identifier)
        mx, mn = max(shape), min(shape)
        if mx == 2 and mn == 0:
            stype, std = "022", (0, 2, 2)
        elif mx == 2:
            stype, std = "112", (1, 1, 2)
        elif mx == 3:
            if shape in [(0, 1, 3), (1, 3, 0), (3, 0, 1)]:
                stype, std = "013", (0, 1, 3)
            else:
                stype, std = "031", (0, 3, 1)
        else:
            stype, std = "004", (0, 0, 4)
        self.shape_type = stype
        rn = 0
        s = tuple(shape)
        while s != std:
            s = (s[2], s[0], s[1])
            rn += 1
        self.rotate_num = rn
        self.split_0_id = None
        if stype in ("112", "022"):
            self.split_0_id = pm.register(1, 0.0, 0.5)
        self.part_frac = None

    def evaluate_init(self):
        self.split_0 = (self.pm.get_scalar(self.split_0_id)
                        if self.split_0_id is not None else None)
        self.part_frac = zeros(1)

    def evaluate_pre(self):
        pass

    def evaluate_post(self):
        # PartInfoLv2.m EvaluatePost, verbatim structure
        import interval_core as IC
        if _interval_mode_active():
            cs1 = [IC.const(0.0) for _ in range(9)]
            cs2 = [IC.const(0.0) for _ in range(9)]
            cs3 = [IC.const(0.0) for _ in range(9)]
        else:
            cs1 = np.zeros(9); cs2 = np.zeros(9); cs3 = np.zeros(9)
        if self.shape_type == "022":
            self.num_block_contribution = (0.0, 0.0, 0.0)
            inner = entropy_vec(np.array([self.split_0, self.split_0,
                                          1 - 2 * self.split_0])) \
                + 2 * np.log(self.q) * (1 - 2 * self.split_0)
            self.mat_size_contribution = (0.0, 0.0, inner)
            cs1[EncodeCSD([0, 0])] = 1.0
            cs2[EncodeCSD([0, 2])] = self.split_0
            cs2[EncodeCSD([2, 0])] = self.split_0
            cs2[EncodeCSD([1, 1])] = 1 - 2 * self.split_0
            cs3 = cs2
        elif self.shape_type == "112":
            self.num_block_contribution = (np.log(2.0), np.log(2.0),
                                           entropy_vec(np.array([self.split_0,
                                                                 self.split_0,
                                                                 1 - 2 * self.split_0])))
            s0 = self.split_0
            self.mat_size_contribution = (2 * s0 * 0.0 + (1 - 2 * s0) * np.log(self.q),
                                          2 * s0 * np.log(self.q) + (1 - 2 * s0) * 0.0,
                                          2 * s0 * 0.0 + (1 - 2 * s0) * np.log(self.q))
            cs1[EncodeCSD([0, 1])] = 0.5
            cs1[EncodeCSD([1, 0])] = 0.5
            cs2 = cs1.copy()
            cs3[EncodeCSD([0, 2])] = s0
            cs3[EncodeCSD([2, 0])] = s0
            cs3[EncodeCSD([1, 1])] = 1 - 2 * s0
        elif self.shape_type in ("013", "031"):
            self.num_block_contribution = (0.0, 0.0, 0.0)
            inner = np.log(2.0) + np.log(self.q)
            self.mat_size_contribution = (0.0, 0.0, inner)
            cs1[EncodeCSD([0, 0])] = 1.0
            if self.shape_type == "013":
                cs2[EncodeCSD([0, 1])] = 0.5
                cs2[EncodeCSD([1, 0])] = 0.5
                cs3[EncodeCSD([1, 2])] = 0.5
                cs3[EncodeCSD([2, 1])] = 0.5
            else:
                cs2[EncodeCSD([1, 2])] = 0.5
                cs2[EncodeCSD([2, 1])] = 0.5
                cs3[EncodeCSD([0, 1])] = 0.5
                cs3[EncodeCSD([1, 0])] = 0.5
        else:  # 004
            self.num_block_contribution = (0.0, 0.0, 0.0)
            self.mat_size_contribution = (0.0, 0.0, 0.0)
            cs1[EncodeCSD([0, 0])] = 1.0
            cs2 = cs1.copy()
            cs3[EncodeCSD([2, 2])] = 1.0
        cs = [cs1, cs2, cs3]
        self.num_block_contribution = Rot3(self.num_block_contribution,
                                           self.rotate_num)
        self.mat_size_contribution = Rot3(self.mat_size_contribution,
                                          self.rotate_num)
        frac = self.part_frac.item()
        self.num_block_contribution = tuple(frac * x for x in
                                            self.num_block_contribution)
        self.mat_size_contribution = tuple(frac * x for x in
                                           self.mat_size_contribution)
        cs = Rot3c(cs, self.rotate_num)
        if _interval_mode_active():
            import interval_core as IC
            flat = [IC.const(v) for c in cs
                    for v in (c if isinstance(c, (list, np.ndarray)) else [c])]
            cs = [flat[i * 9:(i + 1) * 9] for i in range(3)]
        self.complete_split = [G(c) for c in cs]
        self.p_comp = 0.0  # PartInfoLv2.m: placeholder, always 0
        self.hash_penalty_term = 0.0


class PartZero:
    """PartInfoZero.m — level >= 3 with a zero coordinate."""

    def __init__(self, pm, q, level, part_id, shape, identifier):
        self.pm = pm
        self.q = q
        self.level = level
        self.power = 2 ** (level - 1)
        self.part_id = part_id
        self.shape = tuple(shape)
        self.identifier = tuple(identifier)
        if shape[0] == 0:
            self.zero_dim, self.nz1, self.nz2, self.base = 0, 1, 2, (0.0, 0.0, 1.0)
        elif shape[1] == 0:
            self.zero_dim, self.nz1, self.nz2, self.base = 1, 0, 2, (1.0, 0.0, 0.0)
        else:
            self.zero_dim, self.nz1, self.nz2, self.base = 2, 0, 1, (0.0, 1.0, 0.0)
        P = self.power
        self.complete_split_id = [np.full(3 ** P, -1, dtype=np.int64)
                                  for _ in range(3)]
        self.complete_split_id[self.zero_dim][0] = -2
        for idx in range(3 ** P):
            arr = DecodeCSD(idx, P)
            if sum(arr) != self.shape[self.nz1]:
                continue
            self.complete_split_id[self.nz1][idx] = pm.register(1, 0.0, 1.0)
            opp = [2 - a for a in arr]
            self.complete_split_id[self.nz2][EncodeCSD(opp)] = \
                self.complete_split_id[self.nz1][idx]
        self.part_frac = None
        # placeholders (PartInfoZero.m Build lines 95-98)
        self.num_block_contribution = [(0.0, 0.0, 0.0)] * 3
        self.hash_penalty_term = [0.0] * 3
        self.p_comp = [0.0] * 3

    def evaluate_init(self):
        P = self.power
        import interval_core as IC
        interval = _interval_mode_active()
        self.part_frac = zeros(1)
        vals0 = ([IC.const(0.0) for _ in range(3 ** P)] if interval
                 else np.zeros(3 ** P))
        self.complete_split = []
        for t in range(3):
            vals = list(vals0) if interval else np.array(vals0)
            for idx in range(3 ** P):
                gid = self.complete_split_id[t][idx]
                if gid == -1:
                    vals[idx] = IC.const(0.0) if interval else 0.0
                elif gid == -2:
                    vals[idx] = IC.const(1.0) if interval else 1.0
                else:
                    v = self.pm.get_scalar(gid)
                    vals[idx] = (IC.const(v) if interval else v)
            self.complete_split.append(G(vals))

    def evaluate_pre(self):
        pass

    def evaluate_post(self):
        # PartInfoZero.m EvaluatePost: inner_prod_size over compressed CSD
        csd = self.complete_split[self.nz1].v
        inner = entropy_vec(csd)
        for idx in range(len(csd)):
            if csd[idx] != 0.0:
                arr = DecodeCSD(idx, self.power)
                inner += csd[idx] * sum(1 for a in arr if a == 1) * np.log(self.q)
        base = self.part_frac * inner
        if base.__class__.__name__ == "Ival":
            self.mat_size_contribution = tuple(base.v * b for b in self.base)
        elif isinstance(base, G):
            if base.v.__class__.__name__ == "Ival":
                self.mat_size_contribution = tuple(base.v * b for b in self.base)
            else:
                self.mat_size_contribution = tuple(base.v * b for b in self.base)
        else:
            self.mat_size_contribution = tuple(base * b for b in self.base)
        # num_block_contribution / hash_penalty_term / p_comp stay 0


class Part:
    """PartInfo.m — level >= 3, strictly positive shape."""

    def __init__(self, pm, q, level, part_id, shape, identifier, parts):
        self.pm = pm
        self.q = q
        self.level = level
        self.power = 2 ** (level - 1)
        self.sum_col = 2 ** level
        self.sum_half = self.sum_col // 2
        self.part_id = part_id
        self.shape = tuple(shape)
        self.identifier = tuple(identifier)
        self.splits = PrepareSplits(shape)
        self.num_split = len(self.splits)
        self.j2m = JointToMargin(self.splits, self.sum_half)
        lows = [min(s[t] for s in self.splits) for t in range(3)]
        highs = [max(s[t] for s in self.splits) for t in range(3)]
        self.lam_low = lows
        self.lam_high = highs
        # Registration order (PartInfo.m RegisterVariablesAndLinearConstraints):
        # per r in 1..3: split_dist_id{r}, split_dist_max_id{r}; then
        # region_prop_id; then lam_margin{r,t} then lam_sum{r} per r.
        self.split_dist_id = [None] * 3
        self.split_dist_max_id = [None] * 3
        for r in range(3):
            self.split_dist_id[r] = pm.register(self.num_split, 0.0, 1.0)
            self.split_dist_max_id[r] = pm.register(self.num_split, 0.0, 1.0)
        self.region_prop_id = pm.register(3, 0.0, 1.0)
        self.lam_margin_id = {}
        self.lam_sum_id = [None] * 3
        for r in range(3):
            for t in range(3):
                w = highs[t] - lows[t] + 1
                self.lam_margin_id[(r, t)] = pm.register(w, -np.inf, np.inf)
            self.lam_sum_id[r] = pm.register(1, -np.inf, np.inf)
        # linear constraints (for the stage-(b) audit)
        for r in range(3):
            pm.add_lincon_eq([(self.split_dist_id[r], np.ones(self.num_split))], 1.0)
            for t in range(3):
                pm.add_lincon_eq([(self.split_dist_id[r], self.j2m.mats[t]),
                                  (self.split_dist_max_id[r], -self.j2m.mats[t])],
                                 np.zeros(self.sum_half + 1))
        pm.add_lincon_eq([(self.region_prop_id, np.ones(3))], 1.0)
        # children
        self.left_idx = [[None] * self.num_split for _ in range(3)]
        self.right_idx = [[None] * self.num_split for _ in range(3)]
        self.left_ptr = [[None] * self.num_split for _ in range(3)]
        self.right_ptr = [[None] * self.num_split for _ in range(3)]
        for r in range(3):
            ident = (self.part_id, r)  # MATLAB [part_id, r], 1-based part_id
            for i, sp in enumerate(self.splits):
                lid, lptr = find_or_create_part(level - 1, sp[0:3], ident,
                                                parts, pm, q)
                rid, rptr = find_or_create_part(level - 1, sp[3:6], ident,
                                                parts, pm, q)
                self.left_idx[r][i], self.left_ptr[r][i] = lid, lptr
                self.right_idx[r][i], self.right_ptr[r][i] = rid, rptr
        self.part_frac = None

    def evaluate_init(self):
        pm = self.pm
        self.region_prop = G(pm.get(self.region_prop_id))
        self.split_dist = [G(pm.get(self.split_dist_id[r])) for r in range(3)]
        self.split_dist_max = [G(pm.get(self.split_dist_max_id[r]))
                               for r in range(3)]
        self.lam_margin = {}
        self.lam_sum = [None] * 3
        for r in range(3):
            for t in range(3):
                self.lam_margin[(r, t)] = G(pm.get(self.lam_margin_id[(r, t)]))
            self.lam_sum[r] = G(pm.get(self.lam_sum_id[r]))
        self.part_frac = zeros(1)

    def evaluate_pre(self):
        # PartInfo.m EvaluatePre
        for i in range(self.num_split):
            for r in range(3):
                w = self.part_frac * self.split_dist[r][i] * self.region_prop[r]
                l = self.left_ptr[r][i]
                l.part_frac = l.part_frac + w
                rr = self.right_ptr[r][i]
                rr.part_frac = rr.part_frac + w

    def _pcomp_num(self, region):
        """PartInfo.m GetPcompNumerator: shared dim = region (0-based here)."""
        dimW = region
        sh = self.sum_half
        weighted = [zeros(3 ** (self.power // 2)) for _ in range(sh + 1)]
        prob_sum = np.zeros(sh + 1)
        used = np.zeros(sh + 1, dtype=bool)
        res = zeros(1)
        for i in range(self.num_split):
            for ptr in (self.left_ptr[region][i], self.right_ptr[region][i]):
                if ptr.shape[dimW] == 0:
                    continue
                if min(ptr.shape) == 0:
                    res = res + self.split_dist[region][i] \
                        * ptr.complete_split[dimW].entropy()
                else:
                    j = ptr.shape[dimW]
                    weighted[j] = weighted[j] + self.split_dist[region][i] \
                        * ptr.complete_split[dimW]
                    vv = self.split_dist[region][i]
                    if hasattr(vv, 'v'):
                        vv = vv.v
                    if vv.__class__.__name__ == 'Ival':
                        prob_sum[j] += float(vv.lo)
                    elif isinstance(vv, arb):
                        prob_sum[j] += float(vv)
                    else:
                        prob_sum[j] += float(vv)
                    used[j] = True
        for j in range(1, sh + 1):
            if not used[j]:
                continue
            res = res + weighted[j].normalized_entropy(prob_sum[j])
        return res

    def _pcomp_den(self, region):
        """PartInfo.m GetPcompDenominator: H(csr{r,r}) - H(marginal_r)."""
        res = self.complete_split_region[(region, region)].entropy()
        marg = self.j2m.apply(self.split_dist[region].v)
        res = res - G(marg[region]).entropy()
        return res

    def evaluate_post(self):
        self.hash_penalty_term = [None] * 3
        self.num_block_contribution = [None] * 3
        for r in range(3):
            hm = self.split_dist_max[r].entropy() - self.split_dist[r].entropy()
            self.hash_penalty_term[r] = hm * self.part_frac * self.region_prop[r]
            marg = self.j2m.apply(self.split_dist[r].v)
            ent = [G(m).entropy().item() for m in marg]
            self.num_block_contribution[r] = tuple(
                self.part_frac.item() * self.region_prop[r].item() * e
                for e in ent)
        self.complete_split_region = {}
        P = self.power
        for r in range(3):
            for t in range(3):
                acc = zeros(3 ** P)
                for i in range(self.num_split):
                    acc = acc + self.split_dist[r][i] * G(kron_csd(
                        self.left_ptr[r][i].complete_split[t].v,
                        self.right_ptr[r][i].complete_split[t].v))
                self.complete_split_region[(r, t)] = acc
        self.complete_split = [zeros(3 ** P) for _ in range(3)]
        for t in range(3):
            for r in range(3):
                self.complete_split[t] = self.complete_split[t] + \
                    self.complete_split_region[(r, t)] * self.region_prop[r]
        # PartInfo.m: mat_size_contribution should be [0,0,0] at level >= 3
        self.mat_size_contribution = (0.0, 0.0, 0.0)
        self.p_comp = [None] * 3
        for r in range(3):
            w = self.part_frac * self.region_prop[r]
            self.p_comp[r] = (self._pcomp_num(r) - self._pcomp_den(r)) * w

    def lagrange_constraints(self):
        """PartInfo.m GetLagrangeConstraints: exp(sum lam + lam_sum - 1) -
        split_dist_max(t) == 0 for each split left-shape."""
        out = []
        for r in range(3):
            for t, sp in enumerate(self.splits):
                total = self.lam_sum[r]
                for d in range(3):
                    g = self.lam_margin[(r, d)]
                    val = g.v[sp[d] - self.lam_low[d]]
                    if isinstance(val, G):
                        val = val.v
                    total = total + val
                tv = total.v if isinstance(total, G) else total
                if isinstance(tv, list):
                    tv = tv[0]
                if tv.__class__.__name__ == 'Ival':
                    from interval_core import ln_iv, const
                    dm = self.split_dist_max[r].v
                    dmv = dm[t]
                    if hasattr(dmv, 'v'):
                        dmv = dmv.v
                    # MATLAB form exp(tv - 1) - dmv == 0  <=>  log-form
                    # residual ln(dmv) - (tv - 1) == 0. (The +1 matters.)
                    e = ln_iv(dmv) - (tv - const(1.0))
                    out.append((float(e.lo), float(e.hi)))
                elif isinstance(tv, arb):
                    dm = self.split_dist_max[r].v
                    dmv = dm[t] if not isinstance(dm[t], G) else dm[t].v
                    import interval_core as IC
                    ce = IC.exp_iv(total) if False else tv.exp()
                    out.append(np.array(ce - dmv).reshape(-1)[0])
                else:
                    with np.errstate(over="ignore"):
                        constr = np.exp(tv - 1.0) - self.split_dist_max[r].v[t]
                    out.append(constr)
        return out


def find_or_create_part(level, shape, identifier, parts, pm, q):
    arr = parts[level - 1]
    for idx, t in enumerate(arr):
        if t.shape == tuple(shape) and t.identifier == tuple(identifier):
            return idx + 1, t
    if level == 2:
        obj = PartLv2(pm, q, level, len(arr) + 1, shape, identifier)
    elif min(shape) == 0:
        obj = PartZero(pm, q, level, len(arr) + 1, shape, identifier)
    else:
        obj = Part(pm, q, level, len(arr) + 1, shape, identifier, parts)
    arr.append(obj)
    return len(arr), obj


class GlobalStage:
    """GlobalStage.m — 3 regions, single p_comp."""

    def __init__(self, pm, q, max_level, parts):
        self.pm = pm
        self.q = q
        self.level = max_level
        self.power = 2 ** (max_level - 1)
        self.sum_col = 2 ** max_level
        self.shapes = PrepareShapes(max_level)
        self.n_shape = len(self.shapes)
        self.j2m = JointToMargin(self.shapes, self.sum_col)
        self.region_prop_id = pm.register(3, 0.0, 1.0)
        self.dist_id = [None] * 3
        self.dist_max_id = [None] * 3
        for r in range(3):
            self.dist_id[r] = pm.register(self.n_shape, 0.0, 1.0)
            self.dist_max_id[r] = pm.register(self.n_shape, 0.0, 1.0)
        self.lam_margin_id = {}
        self.lam_sum_id = [None] * 3
        for r in range(3):
            for t in range(3):
                self.lam_margin_id[(r, t)] = pm.register(self.sum_col + 1,
                                                         -np.inf, np.inf)
            self.lam_sum_id[r] = pm.register(1, -np.inf, np.inf)
        # extra linear constraints of VXXZ24 (GlobalStage.m:91-112)
        pm.add_lincon_eq([(self.region_prop_id, np.ones(3))], 1.0)
        for r in range(3):
            pm.add_lincon_eq([(self.dist_id[r], np.ones(self.n_shape))], 1.0)
            for t in range(3):
                pm.add_lincon_eq([(self.dist_id[r], self.j2m.mats[t]),
                                  (self.dist_max_id[r], -self.j2m.mats[t])],
                                 np.zeros(self.sum_col + 1))
        # Y and Z symmetric
        pm.add_lincon_eq([(self.region_prop_id, [0.0, 1.0, -1.0])], 0.0)
        # omega mode with K == 1: X, Y, Z all symmetric
        pm.add_lincon_eq([(self.region_prop_id, [1.0, -1.0, 0.0])], 0.0)
        # children per region
        self.part_ptr = [[None] * self.n_shape for _ in range(3)]
        for r in range(3):
            for i, shp in enumerate(self.shapes):
                arr = parts[self.level - 1]
                cur_id = len(arr) + 1
                shp_t = tuple(shp)
                if min(shp_t) == 0:
                    obj = PartZero(pm, q, self.level, cur_id, shp_t, (0, r))
                else:
                    obj = Part(pm, q, self.level, cur_id, shp_t, (0, r), parts)
                arr.append(obj)
                self.part_ptr[r][i] = obj

    def evaluate_init(self):
        pm = self.pm
        self.region_prop = G(pm.get(self.region_prop_id))
        self.dist = [G(pm.get(self.dist_id[r])) for r in range(3)]
        self.dist_max = [G(pm.get(self.dist_max_id[r])) for r in range(3)]
        self.lam_margin = {}
        for r in range(3):
            for t in range(3):
                self.lam_margin[(r, t)] = G(pm.get(self.lam_margin_id[(r, t)]))
        self.lam_sum = [G(pm.get(self.lam_sum_id[r])) for r in range(3)]

    def evaluate_pre(self):
        for r in range(3):
            for i in range(self.n_shape):
                self.part_ptr[r][i].part_frac = \
                    G(self.dist[r].v[i] * self.region_prop.v[r])

    def evaluate_post(self, parts):
        self.hash_penalty_term = [None] * 3
        self.num_block = [None] * 3
        for r in range(3):
            hm = (self.dist_max[r].entropy() - self.dist[r].entropy()) \
                * self.region_prop[r]
            self.hash_penalty_term[r] = hm
            marg = self.j2m.apply(self.dist[r].v)
            self.num_block[r] = tuple(
                G(m).entropy().item() * self.region_prop.v[r] for m in marg)
        # mat_size: sum over all parts of all levels >= 2
        self.mat_size = (0.0, 0.0, 0.0)
        for l in range(2, self.level + 1):
            for t in parts[l - 1]:
                ms = t.mat_size_contribution
                self.mat_size = tuple(a + b for a, b in zip(self.mat_size, ms))
        self.complete_split = [[zeros(3 ** self.power) for _ in range(3)]
                               for _ in range(3)]
        for r in range(3):
            for t in range(3):
                acc = zeros(3 ** self.power)
                for i in range(self.n_shape):
                    acc = acc + self.dist[r][i] * self.part_ptr[r][i].complete_split[t]
                self.complete_split[r][t] = acc
        self.p_comp = [None] * 3
        for r in range(3):
            self.p_comp[r] = (self._pcomp_num(r) - self._pcomp_den(r)) \
                * self.region_prop[r]

    def _pcomp_num(self, region):
        dimW = region
        weighted = [zeros(3 ** self.power) for _ in range(self.sum_col + 1)]
        prob_sum = np.zeros(self.sum_col + 1)
        used = np.zeros(self.sum_col + 1, dtype=bool)
        res = zeros(1)
        for i in range(self.n_shape):
            ptr = self.part_ptr[region][i]
            if ptr.shape[dimW] == 0:
                continue
            if min(ptr.shape) == 0:
                res = res + self.dist[region][i] \
                    * ptr.complete_split[dimW].entropy()
            else:
                k = ptr.shape[dimW]
                weighted[k] = weighted[k] + self.dist[region][i] \
                    * ptr.complete_split[dimW]
                vo = self.dist[region][i].v if hasattr(self.dist[region][i], 'v') \
                    else self.dist[region][i]
                if vo.__class__.__name__ == 'Ival':
                    prob_sum[k] += float(vo.lo)
                elif isinstance(vo, arb):
                    prob_sum[k] += float(vo)
                else:
                    prob_sum[k] += float(vo)
                used[k] = True
        for k in range(1, self.sum_col + 1):
            if not used[k]:
                continue
            res = res + weighted[k].normalized_entropy(prob_sum[k])
        return res

    def _pcomp_den(self, region):
        res = self.complete_split[region][region].entropy()
        marg = self.j2m.apply(self.dist[region].v)
        res = res - G(marg[region]).entropy()
        return res

    def lagrange_constraints(self):
        out = []
        for r in range(3):
            for t, shp in enumerate(self.shapes):
                total = self.lam_sum[r]
                for d in range(3):
                    total = total + self.lam_margin[(r, d)].v[shp[d]]
                tv = total.v if isinstance(total, G) else total
                if isinstance(tv, list):
                    tv = tv[0]
                dm_vec = self.dist_max[r].v
                dmv = dm_vec[t]
                if hasattr(dmv, 'v'):
                    dmv = dmv.v
                if tv.__class__.__name__ == 'Ival':
                    from interval_core import ln_iv, const
                    e = ln_iv(dmv) - (tv - const(1.0))
                    out.append((float(e.lo), float(e.hi)))
                elif isinstance(tv, arb):
                    ce = tv.exp()
                    out.append(np.array(ce - dmv).reshape(-1)[0])
                else:
                    with np.errstate(over="ignore"):
                        constr = np.exp(tv - 1.0) - dmv
                    out.append(constr)
        return out


class Workspace:
    """Workspace.m — 3-region version with explicit Schonhage line."""

    def __init__(self, pm, q, K, omega, max_level=MAX_LEVEL):
        self.pm = pm
        self.q = q
        self.K = K
        self.max_level = max_level
        self.parts = [[] for _ in range(max_level)]
        # VXXZ24 Workspace.Build registration order:
        # per r: num_retain_comp{r}{lv} for lv = 2..max_level; then
        # num_retain_glob{r}; then single_mat_size; then omega; then K.
        self.num_retain_comp_id = {}
        self.num_retain_glob_id = [None] * 3
        for r in range(3):
            for lv in range(2, max_level + 1):
                self.num_retain_comp_id[(r, lv)] = pm.register(1, 0.0, np.inf)
            self.num_retain_glob_id[r] = pm.register(1, 0.0, np.inf)
        self.single_mat_size_id = pm.register(1, 0.0, np.inf)
        self.omega_id = pm.register(1, 0.0, np.inf)
        self.K_id = pm.register(1, K, K)
        self.globstage = GlobalStage(pm, q, max_level, self.parts)

    def omega_value(self):
        return self.pm.get_scalar(self.omega_id)

    def evaluate(self):
        pm = self.pm
        L = self.max_level
        c_viol = []
        ceq_viol = []
        self.ret_glob = [G(pm.get(self.num_retain_glob_id[r])) for r in range(3)]
        self.ret_comp = {}
        for r in range(3):
            for lv in range(2, L + 1):
                self.ret_comp[(r, lv)] = G(pm.get(self.num_retain_comp_id[(r, lv)]))
        self.single_mat_size = G(pm.get(self.single_mat_size_id))
        self.omega = G(pm.get(self.omega_id))
        self.globstage.evaluate_init()
        for l in range(1, L + 1):
            for t in self.parts[l - 1]:
                t.evaluate_init()
        self.globstage.evaluate_pre()
        for l in range(L, 0, -1):
            for t in self.parts[l - 1]:
                t.evaluate_pre()
        for l in range(1, L + 1):
            for t in self.parts[l - 1]:
                t.evaluate_post()
        self.globstage.evaluate_post(self.parts)
        # Lagrange constraints (PartInfo only, per Workspace.m:158-168)
        for l in range(1, L + 1):
            for t in self.parts[l - 1]:
                if isinstance(t, Part):
                    ceq_viol.extend(t.lagrange_constraints())
        ceq_viol.extend(self.globstage.lagrange_constraints())
        # asymmetric hashing level >= 3 (Workspace.m:171-195)
        for l in range(3, L + 1):
            for r in range(3):
                num_block = [0.0, 0.0, 0.0]
                pen = zeros(1)
                pcomp = zeros(1)
                for t in self.parts[l - 1]:
                    nb = t.num_block_contribution[r]
                    num_block = [a + b for a, b in zip(num_block, nb)]
                    pen = pen + t.hash_penalty_term[r]
                    pcomp = pcomp + t.p_comp[r]
                for tt in range(3):
                    if tt == r:
                        c_viol.append((self.ret_comp[(r, l)] - num_block[tt]
                                       + pcomp).item())
                    else:
                        c_viol.append((self.ret_comp[(r, l)] - num_block[tt]
                                       + pen).item())
        # symmetric hashing at level 2 (Workspace.m:196-212)
        num_block = [0.0, 0.0, 0.0]
        for t in self.parts[1]:
            nb = t.num_block_contribution
            num_block = [a + b for a, b in zip(num_block, nb)]
        for tt in range(3):
            c_viol.append((self.ret_comp[(0, 2)] - num_block[tt]).item())
        for r in range(1, 3):
            ceq_viol.append(self.ret_comp[(r, 2)].item())
        # global hashing (Workspace.m:214-231)
        for r in range(3):
            num_block = self.globstage.num_block[r]
            pen = self.globstage.hash_penalty_term[r]
            pcomp = self.globstage.p_comp[r]
            for tt in range(3):
                if tt == r:
                    c_viol.append((self.ret_glob[r] - num_block[tt]
                                   + pcomp).item())
                else:
                    c_viol.append((self.ret_glob[r] - num_block[tt]
                                   + pen).item())
        # single matrix size (Workspace.m:233-239)
        ms = self.globstage.mat_size
        for tt in range(3):
            denom = 1.0
            if tt == 2:
                denom = self.K
            c_viol.append((self.single_mat_size - ms[tt] / denom).item())
        # value + Schonhage line (Workspace.m:241-252)
        value = (self.single_mat_size * self.omega).item()
        for r in range(3):
            value += self.ret_glob[r].item()
        for lv in range(2, L + 1):
            for r in range(3):
                value += self.ret_comp[(r, lv)].item()
        self.value = value
        from interval_core import Ival as _Iv
        target_iv = np.log(self.q + 2) * (2 ** (L - 1))
        if isinstance(value, _Iv):
            sv = _Iv(target_iv - value.hi, target_iv - value.lo)
        elif isinstance(value, tuple):
            sv = (target_iv - value[1], target_iv - value[0])
        else:
            sv = target_iv - float(value)
        c_viol.append(sv)

        def _scalarize(x):
            if isinstance(x, tuple):
                return x  # interval (lo, hi) pair — keep both endpoints
            if x.__class__.__name__ == 'Ival':
                return (float(x.lo), float(x.hi))
            if isinstance(x, list):
                z = x[0]
                if z.__class__.__name__ == 'Ival':
                    return (float(z.lo), float(z.hi))
            return float(np.asarray(x).reshape(-1)[0])
        c_pair = [_scalarize(x) for x in c_viol]
        ceq_pair = [_scalarize(x) for x in ceq_viol]
        c_viol_out = [p[0] if isinstance(p, tuple) else p for p in c_pair]
        ceq_out = [p for p in ceq_pair]
        return c_viol_out, ceq_out, value


def build_and_load(mat_path, q=Q, K=1.0, max_level=MAX_LEVEL):
    pm = ParamManager()
    ws = Workspace(pm, q, K, omega=0.0, max_level=max_level)
    data = loadmat(mat_path)
    params = np.asarray(data["params"]).flatten()
    if len(pm.lb) != len(params):
        raise SystemExit(f"param count mismatch: registered {len(pm.lb)} "
                         f"vs file {len(params)}")
    pm.set_value(params)
    return pm, ws, params


def max_violation(pm, ws):
    c, ceq, value = ws.evaluate()

    def _flo(z):
        if isinstance(z, tuple):
            return z[1]
        if z.__class__.__name__ == 'Ival':
            return float(z.hi)
        return float(z)
    cf = [_flo(z) for z in c if z is not None]
    cef = [abs(_flo(z)) for z in ceq if z is not None]
    cmax = max(cf) if cf else 0.0
    cemax = max(cef) if cef else 0.0
    linv = pm.linear_violations()
    return max(cmax, cemax, linv), cmax, cemax, value, c, ceq


if __name__ == "__main__":
    import sys
    mat = sys.argv[1] if len(sys.argv) > 1 else \
        "/Users/jinleic/jinleic-workspace/cs/omega/scratch/rmmcode/data/K100_2.37155181.mat"
    pm, ws, params = build_and_load(mat)
    print(f"registered params: {pm.num_input}  file params: {len(params)}")
    print(f"omega in file: {ws.omega_value():.10f}")
    v, cmax, cemax, value, c, ceq = max_violation(pm, ws)
    print(f"max c violation      : {cmax:.6e}")
    print(f"max |ceq| violation  : {cemax:.6e}")
    print(f"max linear violation : {pm.linear_violations():.6e}")
    print(f"GetFeasibility maxViol: {v:.6e}")
    print(f"value (Schonhage LHS): {float(value):.10f}")
    print(f"target 2^(L-1) log(q+2) = {np.log(7) * 4:.10f}")
