"""Stage-(a) float64 transcription of the Alman et al. 2025 (SODA'25) released
verification program (`code_matrix_mult.zip`, osf.io/mw5ak), used to certify
omega <= 2.371339 from data/W1.00_2.371339.mat (24855 float64 parameters).

Faithful transcription of:
  src/autograd/{GVar,ParamManager}.m   (param registration ORDER and bounds)
  src/evaluation/{Workspace,GlobalStage,TermInfo,TermInfoZero,TermInfoLv2,
                  FindOrCreateTerm,FindTermByIdentifier}.m
  src/complete_split/{EncodeCSD,DecodeCSD,ConcatCSD}.m
  src/utils/{PrepareShapes,PrepareSplits,JointToMargin,Dims,Rot3,Rot3c,
             MarginalDist}.m
  src/verify/{VerifyOmega,GetFeasibility}.m

Gradient machinery is dropped (verification needs values only). All logs are
MATLAB `log` = natural log (the released code is e-based throughout).
"""

import numpy as np
from flint import arb
from scipy.io import loadmat

Q = 5.0
MAX_LEVEL = 3

DIM_PERMS = [(1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1)]


def Dims(p):  # 1-based like MATLAB
    return DIM_PERMS[p - 1]


def PrepareShapes(level):
    sc = 2 ** level
    return [(i, j, sc - i - j) for i in range(sc + 1) for j in range(sc - i + 1)]


def PrepareSplits(shape):
    """1-based 6-tuples (i,j,k,i2,j2,k2) like MATLAB rows."""
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
    """column: list of shape tuples; returns list of 3 dense matrices
    (n_col x (sum_col+1))."""

    def __init__(self, columns, sum_col):
        n = len(columns)
        self.mats = [np.zeros((n, sum_col + 1)) for _ in range(3)]
        for i, col in enumerate(columns):
            for t in range(3):
                self.mats[t][i, col[t]] = 1.0

    def apply(self, dist):
        return [dist @ m for m in self.mats]  # 3 marginal vectors


def Rot3(x, n):
    n %= 3
    if n == 1:
        return (x[1], x[2], x[0])
    if n == 2:
        return (x[2], x[0], x[1])
    return x


def Rot3c(x, n):
    n %= 3
    if n == 1:
        return (x[1], x[2], x[0])
    if n == 2:
        return (x[2], x[0], x[1])
    return x


def EncodeCSD(arr):
    idx = 0
    for a in arr:
        idx = idx * 3 + a
    return idx  # 0-based; MATLAB returns idx+1


def DecodeCSD(identity, power):
    arr = [0] * power
    idx = identity  # 0-based
    for i in range(power - 1, -1, -1):
        arr[i] = idx % 3
        idx //= 3
    return arr


# ------------------------------------------------------------------ entropy

def entropy_vec(v):
    """GVar.Entropy value part: sum -x*log(x + (x<=0)) elementwise."""
    v = np.asarray(v, dtype=np.float64)
    shifted = np.where(v <= 0.0, v + 1.0, v)
    with np.errstate(divide="ignore", invalid="ignore"):
        logs = np.log(shifted)
    return float(np.sum(-v * logs))


def normalized_entropy(a, p):
    """GVar.NormalizedEntropy value: sum -a_i log(a_i / p)."""
    a = np.asarray(a, dtype=np.float64)
    if p <= 0.0:
        return 0.0
    ratio = np.where(a <= 0.0, a / p + 1.0, a / p)
    with np.errstate(divide="ignore", invalid="ignore"):
        logs = np.log(ratio)
    return float(np.sum(-a * logs))


def kron_csd(lhs, rhs):
    return np.kron(lhs, rhs)


# -------------------------------------------------------------- param manager

class ParamManager:
    def __init__(self):
        self.start = []      # 0-based start of each group
        self.size = []
        self.lb = []
        self.ub = []
        self.cur_x = np.zeros(0)

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


# ------------------------ interval (arb) primitive ops --------------------
# All run inside a caller-supplied flint ctx.workprec block. Outward rounding:
# every inequality consumes explicit endpoints; entropy intervals enclose
# -t ln t over the interval argument (t>0 assumed for support entries, exact
# zero handled separately).

def alman_interval_scalar(x):
    if isinstance(x, arb):
        return x
    a = np.asarray(x)
    return arb(float(a.reshape(-1)[0]))


def ent_iv(x: arb) -> arb:
    """Interval enclosure of -t*ln(t) for t in x (support >= 0)."""
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
    """-x ln(x/p) = -x ln x + x ln p for x interval."""
    if x.upper() <= 0:
        return arb(0)
    return ent_iv(x) + x * p.log()


def nent_vec_iv(ws, p: arb) -> arb:
    t = arb(0)
    for wi in ws:
        if wi.upper() <= 0:
            continue
        t = t + nent_iv(wi, p)
    return t


class G:
    __doc__ = "Value holder emulating a GVar row vector (values only)."
    __slots__ = ("v",)

    def __init__(self, v):
        # interval scalars/vectors pass through untouched
        if isinstance(v, arb):
            self.v = v
            return
        if isinstance(v, list) and v and isinstance(v[0], arb):
            self.v = list(v)
            return
        if type(v).__name__ in ("Vec", "Tag"):   # interval wrappers
            self.v = v
            return
        self.v = np.asarray(v, dtype=np.float64)

    # interval helpers -------------------------------------------------
    @staticmethod
    def _arb(v):
        return isinstance(v, arb)

    @staticmethod
    def _iv(x):
        """x is a scalar (python float / np scalar) — point interval."""
        if isinstance(x, np.ndarray):
            raise TypeError("array in scalar slot")
        return alman_interval_scalar(x)

    def _comb(self, a, b, op):
        """combine two operands (arb | float | list[arb] | np.ndarray) under op.
        np.ndarray operands (float mode) force float path."""
        if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
            return op(a, b)
        if isinstance(a, list) and isinstance(b, list):
            return [op(x, y) for x, y in zip(a, b)]
        if isinstance(a, list):
            b = self._iv(b) if not self._arb(b) else b
            return [op(x, b) for x in a]
        if isinstance(b, list):
            a = self._iv(a) if not self._arb(a) else a
            return [op(a, y) for y in b]
        if not self._arb(a):
            a = self._iv(a)
        if not self._arb(b):
            b = self._iv(b)
        return op(a, b)

    # arithmetic -------------------------------------------------------
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
        if isinstance(v, list):
            return G(ent_vec_iv(v))
        if isinstance(v, arb):
            return G(ent_iv(v))
        return G(entropy_vec(self.v))

    def normalized_entropy(self, p):
        pv = p.v if isinstance(p, G) else p
        v = self.v
        if isinstance(v, list):
            return G(nent_vec_iv(v, pv))
        if isinstance(v, arb):
            return G(nent_iv(v, pv))
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
        if isinstance(self.v, arb):
            return self.v
        if isinstance(self.v, list):
            assert len(self.v) == 1
            return self.v[0]
        return float(np.asarray(self.v).reshape(-1)[0])

    def __len__(self):
        return len(self.v)

    def __getitem__(self, k):
        v = self.v
        if isinstance(v, arb):
            return G(v)          # scalar already (rare)
        if isinstance(v, list):
            return G(v[k])       # index the interval vector
        return G(self.v[k])
        return G(self.v[k])

def zeros(n):
    return G(np.zeros(n))


# ------------------------------------------------------------------- terms

class TermLv2:
    """TermInfoLv2."""

    def __init__(self, pm, q, level, term_id, shape, identifier):
        assert level == 2
        self.pm = pm
        self.q = q
        self.level = 2
        self.power = 2
        self.term_id = term_id
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
        self.term_frac = None

    def evaluate_init(self):
        self.split_0 = (self.pm.get_scalar(self.split_0_id)
                        if self.split_0_id is not None else None)
        self.term_frac = zeros(1)

    def evaluate_pre(self):
        pass

    def evaluate_post(self):
        q = self.q
        s0 = self.split_0 if self.split_0 is not None else 0.0
        cs1 = np.zeros(9)
        cs2 = np.zeros(9)
        cs3 = np.zeros(9)
        if self.shape_type == "022":
            self.num_block_contrib_value = (0.0, 0.0, 0.0)
            self.mat_contrib_value = (0.0, 0.0,
                                      entropy_vec(np.array([s0, s0, 1 - 2 * s0]))
                                      + 2 * np.log(q) * (1 - 2 * s0))
            cs1[EncodeCSD([0, 0])] = 1.0
            cs2[EncodeCSD([0, 2])] = s0
            cs2[EncodeCSD([2, 0])] = s0
            cs2[EncodeCSD([1, 1])] = 1 - 2 * s0
            cs3 = cs2
        elif self.shape_type == "112":
            self.num_block_contrib_value = (np.log(2.0), np.log(2.0),
                                            entropy_vec(np.array([s0, s0, 1 - 2 * s0])))
            self.mat_contrib_value = (2 * s0 * 0.0 + (1 - 2 * s0) * np.log(q),
                                      2 * s0 * np.log(q) + (1 - 2 * s0) * 0.0,
                                      2 * s0 * 0.0 + (1 - 2 * s0) * np.log(q))
            cs1[EncodeCSD([0, 1])] = 0.5
            cs1[EncodeCSD([1, 0])] = 0.5
            cs2 = cs1.copy()
            cs3[EncodeCSD([0, 2])] = s0
            cs3[EncodeCSD([2, 0])] = s0
            cs3[EncodeCSD([1, 1])] = 1 - 2 * s0
        elif self.shape_type in ("013", "031"):
            self.num_block_contrib_value = (0.0, 0.0, 0.0)
            inner = np.log(2.0) + np.log(q)
            self.mat_contrib_value = (0.0, 0.0, inner)
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
            self.num_block_contrib_value = (0.0, 0.0, 0.0)
            self.mat_contrib_value = (0.0, 0.0, 0.0)
            cs1[EncodeCSD([0, 0])] = 1.0
            cs2 = cs1.copy()
            cs3[EncodeCSD([2, 2])] = 1.0
        cs = [cs1, cs2, cs3]
        self.num_block = Rot3(self.num_block_contrib_value, self.rotate_num)
        self.mat_size = Rot3(self.mat_contrib_value, self.rotate_num)
        frac = self.term_frac.item()
        self.num_block = tuple(frac * x for x in self.num_block)
        self.mat_size = tuple(frac * x for x in self.mat_size)
        cs = Rot3c(cs, self.rotate_num)
        self.complete_split = [G(c) for c in cs]


class TermZero:
    """TermInfoZero (level >= 3, shape has a zero coordinate)."""

    def __init__(self, pm, q, level, term_id, shape, identifier):
        self.pm = pm
        self.q = q
        self.level = level
        self.power = 2 ** (level - 1)
        self.term_id = term_id
        self.shape = tuple(shape)
        self.identifier = tuple(identifier)
        if shape[0] == 0:
            self.zero_dim, self.nz1, self.nz2, self.base = 0, 1, 2, (0.0, 0.0, 1.0)
        elif shape[1] == 0:
            self.zero_dim, self.nz1, self.nz2, self.base = 1, 0, 2, (1.0, 0.0, 0.0)
        else:
            self.zero_dim, self.nz1, self.nz2, self.base = 2, 0, 1, (0.0, 1.0, 0.0)
        P = self.power
        self.csd_id = [None] * 3
        for t in range(3):
            self.csd_id[t] = np.full(3 ** P, -1, dtype=np.int64)
        self.csd_id[self.zero_dim][0] = -2
        for idx in range(3 ** P):
            arr = DecodeCSD(idx, P)
            if sum(arr) != self.shape[self.nz1]:
                continue
            self.csd_id[self.nz1][idx] = pm.register(1, 0.0, 1.0)
            opp = [2 - a for a in arr]
            self.csd_id[self.nz2][EncodeCSD(opp)] = self.csd_id[self.nz1][idx]
        self.term_frac = None

    def evaluate_init(self):
        P = self.power
        self.term_frac = zeros(1)
        self.complete_split = []
        for t in range(3):
            vals = np.zeros(3 ** P)
            for idx in range(3 ** P):
                gid = self.csd_id[t][idx]
                if gid == -1:
                    vals[idx] = 0.0
                elif gid == -2:
                    vals[idx] = 1.0
                else:
                    vals[idx] = self.pm.get_scalar(gid)
            self.complete_split.append(G(vals))

    def evaluate_pre(self):
        pass

    def evaluate_post(self):
        q = self.q
        csd = self.complete_split[self.nz1].v
        inner = entropy_vec(csd)
        for idx in range(len(csd)):
            if csd[idx] != 0.0:
                arr = DecodeCSD(idx, self.power)
                inner += csd[idx] * sum(1 for a in arr if a == 1) * np.log(q)
        base = self.term_frac * inner
        self.mat_size = tuple(base.v * b for b in self.base)
        # MATLAB placeholder: num_block_contribution / hash_penalty_term /
        # p_compY / p_compZ are all-[0,0,0] / 0 for zero-shape terms.
        self.num_block_contribution = [(0.0, 0.0, 0.0)] * 6
        self.hash_penalty_term = [0.0] * 6
        self.p_compY = [0.0] * 6
        self.p_compZ = [0.0] * 6
        self.num_block = (0.0, 0.0, 0.0)


class Term:
    """TermInfo (level >= 3, strictly positive shape)."""
    def __init__(self, pm, q, level, term_id, shape, identifier, terms):
        self.pm = pm
        self.q = q
        self.level = level
        self.power = 2 ** (level - 1)
        self.sum_col = 2 ** level
        self.sum_half = self.sum_col // 2
        self.term_id = term_id
        self.shape = tuple(shape)
        self.identifier = tuple(identifier)
        self.splits = PrepareSplits(shape)
        self.num_split = len(self.splits)
        self.j2m = JointToMargin(self.splits, self.sum_half)
        lows = [min(s[t] for s in self.splits) for t in range(3)]
        highs = [max(s[t] for s in self.splits) for t in range(3)]
        self.lam_low = lows
        self.lam_high = highs
        # Registration order (exact MATLAB order):
        self.split_dist_id = [None] * 6
        self.split_dist_max_id = [None] * 6
        for r in range(6):
            self.split_dist_id[r] = pm.register(self.num_split, 0.0, 1.0)
            self.split_dist_max_id[r] = pm.register(self.num_split, 0.0, 1.0)
        self.region_prop_id = pm.register(6, 0.0, 1.0)
        self.lam_margin_id = {}
        self.lam_sum_id = [None] * 6
        for r in range(6):
            for t in range(3):
                w = highs[t] - lows[t] + 1
                self.lam_margin_id[(r, t)] = pm.register(w, -np.inf, np.inf)
            self.lam_sum_id[r] = pm.register(1, -np.inf, np.inf)
        # Build children (in exact order).
        self.left_idx = [[None] * self.num_split for _ in range(6)]
        self.right_idx = [[None] * self.num_split for _ in range(6)]
        self.left_ptr = [[None] * self.num_split for _ in range(6)]
        self.right_ptr = [[None] * self.num_split for _ in range(6)]
        for r in range(6):
            ident = (self.term_id, r)  # [term_id, r] with 1-based term_id
            for i, sp in enumerate(self.splits):
                lid, lptr = find_or_create_term(level - 1, sp[0:3], ident, terms, pm, q)
                rid, rptr = find_or_create_term(level - 1, sp[3:6], ident, terms, pm, q)
                self.left_idx[r][i], self.left_ptr[r][i] = lid, lptr
                self.right_idx[r][i], self.right_ptr[r][i] = rid, rptr
        self.term_frac = None

    def evaluate_init(self):
        pm = self.pm
        self.region_prop = G(pm.get(self.region_prop_id))
        self.split_dist = [G(pm.get(self.split_dist_id[r])) for r in range(6)]
        self.split_dist_max = [G(pm.get(self.split_dist_max_id[r])) for r in range(6)]
        self.lam_margin = {}
        self.lam_sum = [None] * 6
        for r in range(6):
            for t in range(3):
                g = self.lam_margin_id[(r, t)]
                lo = self.lam_low[t]
                hi = self.lam_high[t]
                self.lam_margin[(r, t)] = G(pm.get(g))
            self.lam_sum[r] = G(pm.get(self.lam_sum_id[r]))
        self.term_frac = zeros(1)

    def evaluate_pre(self):
        for i in range(self.num_split):
            for r in range(6):
                w = self.term_frac * self.split_dist[r][i] * self.region_prop[r]
                l = self.left_ptr[r][i]
                l.term_frac = l.term_frac + w
                rr = self.right_ptr[r][i]
                rr.term_frac = rr.term_frac + w

    def _pcomp_num(self, region, which):
        """Numerator part of the p_comp for compatibilities in dim W.
        which='Y': special when shape Z==0; which='Z': special when any 0.
        MATLAB dims are 1-based; shapes stored 0-based (indices shifted by 1)."""
        dimx, dimy, dimz = Dims(region + 1)
        dimW = (dimy - 1) if which == "Y" else (dimz - 1)
        dimS = (dimz - 1)
        sh = self.sum_half
        weighted = [zeros(3 ** (self.power // 2)) for _ in range(sh + 1)]
        prob_sum = np.zeros(sh + 1)
        used = np.zeros(sh + 1, dtype=bool)
        res = zeros(1)
        for i in range(self.num_split):
            for ptr in (self.left_ptr[region][i], self.right_ptr[region][i]):
                if ptr.shape[dimW] == 0:
                    continue
                if which == "Y":
                    special = ptr.shape[dimS] == 0
                else:
                    special = min(ptr.shape) == 0
                if special:
                    res = res + self.split_dist[region][i] \
                        * ptr.complete_split[dimW].entropy()
                else:
                    j = ptr.shape[dimW]
                    weighted[j] = weighted[j] + self.split_dist[region][i] \
                        * ptr.complete_split[dimW]
                    prob_sum[j] += self.split_dist[region][i].item()
                    used[j] = True
        for j in range(1, sh + 1):
            if not used[j]:
                continue
            res = res + weighted[j].normalized_entropy(prob_sum[j])
        return res

    def _pcomp_den(self, region, which):
        dimx, dimy, dimz = Dims(region + 1)
        dimW = (dimy - 1) if which == "Y" else (dimz - 1)
        res = self.complete_split_region[(region, dimW)].entropy()
        marg = self.j2m.apply(self.split_dist[region].v)
        res = res - G(marg[dimW]).entropy()
        return res

    def evaluate_post(self):
        # hash penalty
        self.hash_penalty_term = [None] * 6
        self.num_block_contribution = [None] * 6
        for r in range(6):
            hm = self.split_dist_max[r].entropy() - self.split_dist[r].entropy()
            self.hash_penalty_term[r] = hm * self.term_frac * self.region_prop[r]
            marg = self.j2m.apply(self.split_dist[r].v)
            ent = [G(m).entropy().item() for m in marg]
            self.num_block_contribution[r] = tuple(
                self.term_frac.item() * self.region_prop[r].item() * e for e in ent)
        # complete splits
        self.complete_split_region = {}
        P = self.power
        for r in range(6):
            for t in range(3):
                acc = zeros(3 ** P)
                for i in range(self.num_split):
                    acc = acc + self.split_dist[r][i] * G(kron_csd(
                        self.left_ptr[r][i].complete_split[t].v,
                        self.right_ptr[r][i].complete_split[t].v))
                self.complete_split_region[(r, t)] = acc
        self.complete_split = [zeros(3 ** P) for _ in range(3)]
        for t in range(3):
            for r in range(6):
                self.complete_split[t] = self.complete_split[t] + \
                    self.complete_split_region[(r, t)] * self.region_prop[r]
        # MATLAB: level>=3 strictly-positive terms do NOT contribute to
        # mat_size ("mat_size_contribution should be [0,0,0]"); only level-2
        # and zero-shape terms carry it. Keep the attribute for uniform use.
        self.mat_size = (0.0, 0.0, 0.0)
        # p_comp
        self.p_compY = [None] * 6
        self.p_compZ = [None] * 6
        for r in range(6):
            w = self.term_frac * self.region_prop[r]
            self.p_compY[r] = (self._pcomp_num(r, "Y") - self._pcomp_den(r, "Y")) * w
            self.p_compZ[r] = (self._pcomp_num(r, "Z") - self._pcomp_den(r, "Z")) * w

    def lagrange_constraints(self):
        out = []
        for r in range(6):
            for t, sp in enumerate(self.splits):
                total = self.lam_sum[r]
                for d in range(3):
                    g = self.lam_margin[(r, d)]
                    val = g.v[sp[d] - self.lam_low[d]]
                    if isinstance(val, G):
                        val = val.v
                    total = total + val
                # interval/flexible exp: exp(total - 1) - dist_max(t)
                tv = total.v if isinstance(total, G) else total
                if isinstance(tv, list):
                    tv = tv[0]
                if isinstance(tv, arb):
                    ce = tv.exp() if hasattr(tv, 'exp') else arb(tv).__exp__()
                    dm = self.split_dist_max[r].v
                    dmv = dm[t] if not isinstance(dm[t], G) else dm[t].v
                    out.append(np.array(ce - dmv).reshape(-1)[0])
                else:
                    constr = np.exp(tv - 1.0) - self.split_dist_max[r].v[t]
                    out.append(constr)
        return out


def find_or_create_term(level, shape, identifier, terms, pm, q):
    arr = terms[level - 1]
    for idx, t in enumerate(arr):
        if t.shape == tuple(shape) and t.identifier == tuple(identifier):
            return idx + 1, t
    if level == 2:
        obj = TermLv2(pm, q, level, len(arr) + 1, shape, identifier)
    elif min(shape) == 0:
        obj = TermZero(pm, q, level, len(arr) + 1, shape, identifier)
    else:
        obj = Term(pm, q, level, len(arr) + 1, shape, identifier, terms)
    arr.append(obj)
    return len(arr), obj


# -------------------------------------------------------------- global stage

class GlobalStage:
    def __init__(self, pm, q, max_level, terms):
        self.pm = pm
        self.q = q
        self.level = max_level
        self.power = 2 ** (max_level - 1)
        self.sum_col = 2 ** max_level
        self.shapes = PrepareShapes(max_level)
        self.n_shape = len(self.shapes)
        self.j2m = JointToMargin(self.shapes, self.sum_col)
        self.region_prop_id = pm.register(6, 0.0, 1.0)
        self.dist_id = [None] * 6
        self.dist_max_id = [None] * 6
        for r in range(6):
            self.dist_id[r] = pm.register(self.n_shape, 0.0, 1.0)
            self.dist_max_id[r] = pm.register(self.n_shape, 0.0, 1.0)
        self.lam_margin_id = {}
        self.lam_sum_id = [None] * 6
        for r in range(6):
            for t in range(3):
                self.lam_margin_id[(r, t)] = pm.register(self.sum_col + 1,
                                                         -np.inf, np.inf)
            self.lam_sum_id[r] = pm.register(1, -np.inf, np.inf)
        self.term_id = [[None] * self.n_shape for _ in range(6)]
        self.term_ptr = [[None] * self.n_shape for _ in range(6)]
        for r in range(6):
            for i, shp in enumerate(self.shapes):
                arr = terms[self.level - 1]
                cur_id = len(arr) + 1
                shp_t = tuple(shp)
                if min(shp_t) == 0:
                    obj = TermZero(pm, q, self.level, cur_id, shp_t, (0, r))
                else:
                    obj = Term(pm, q, self.level, cur_id, shp_t, (0, r), terms)
                arr.append(obj)
                self.term_id[r][i] = cur_id
                self.term_ptr[r][i] = obj

    def evaluate_init(self):
        pm = self.pm
        self.region_prop = G(pm.get(self.region_prop_id))
        self.dist = [G(pm.get(self.dist_id[r])) for r in range(6)]
        self.dist_max = [G(pm.get(self.dist_max_id[r])) for r in range(6)]
        self.lam_margin = {}
        for r in range(6):
            for t in range(3):
                self.lam_margin[(r, t)] = G(pm.get(self.lam_margin_id[(r, t)]))
        self.lam_sum = [G(pm.get(self.lam_sum_id[r])) for r in range(6)]

    def evaluate_pre(self):
        for r in range(6):
            for i in range(self.n_shape):
                self.term_ptr[r][i].term_frac = \
                    G(self.dist[r].v[i] * self.region_prop.v[r])

    def evaluate_post(self, terms):
        self.hash_penalty_term = [None] * 6
        self.num_block = [None] * 6
        for r in range(6):
            hm = (self.dist_max[r].entropy() - self.dist[r].entropy()) \
                * self.region_prop[r]
            self.hash_penalty_term[r] = hm
            marg = self.j2m.apply(self.dist[r].v)
            self.num_block[r] = tuple(
                G(m).entropy().item() * self.region_prop.v[r] for m in marg)
        self.mat_size = (0.0, 0.0, 0.0)
        for l in range(2, self.level + 1):
            for t in terms[l - 1]:
                ms = t.mat_size
                self.mat_size = tuple(a + b for a, b in zip(self.mat_size, ms))
        self.complete_split = [[zeros(3 ** self.power) for _ in range(3)]
                               for _ in range(6)]
        for r in range(6):
            for t in range(3):
                acc = zeros(3 ** self.power)
                for i in range(self.n_shape):
                    acc = acc + self.dist[r][i] * self.term_ptr[r][i].complete_split[t]
                self.complete_split[r][t] = acc
        self.p_compY = [None] * 6
        self.p_compZ = [None] * 6
        for r in range(6):
            self.p_compY[r] = (self._pcomp_num(r, "Y") - self._pcomp_den(r, "Y")) \
                * self.region_prop[r]
            self.p_compZ[r] = (self._pcomp_num(r, "Z") - self._pcomp_den(r, "Z")) \
                * self.region_prop[r]

    def _pcomp_num(self, region, which):
        dimx, dimy, dimz = Dims(region + 1)
        dimW = (dimy - 1) if which == "Y" else (dimz - 1)
        dimS = (dimz - 1)
        weighted = [zeros(3 ** self.power) for _ in range(self.sum_col + 1)]
        prob_sum = np.zeros(self.sum_col + 1)
        used = np.zeros(self.sum_col + 1, dtype=bool)
        res = zeros(1)
        for i in range(self.n_shape):
            ptr = self.term_ptr[region][i]
            if ptr.shape[dimW] == 0:
                continue
            if which == "Y":
                special = ptr.shape[dimS] == 0
            else:
                special = min(ptr.shape) == 0
            if special:
                res = res + self.dist[region][i] * ptr.complete_split[dimW].entropy()
            else:
                k = ptr.shape[dimW]
                weighted[k] = weighted[k] + self.dist[region][i] * ptr.complete_split[dimW]
                prob_sum[k] += self.dist[region][i].item()
                used[k] = True
        for k in range(1, self.sum_col + 1):
            if not used[k]:
                continue
            res = res + weighted[k].normalized_entropy(prob_sum[k])
        return res

    def _pcomp_den(self, region, which):
        dimx, dimy, dimz = Dims(region + 1)
        dimW = (dimy - 1) if which == "Y" else (dimz - 1)
        res = self.complete_split[region][dimW].entropy()
        marg = self.j2m.apply(self.dist[region].v)
        res = res - G(marg[dimW]).entropy()
        return res

    def lagrange_constraints(self):
        out = []
        for r in range(6):
            for t, shp in enumerate(self.shapes):
                total = self.lam_sum[r]
                for d in range(3):
                    total = total + self.lam_margin[(r, d)].v[shp[d]]
                tv = total.v if isinstance(total, G) else total
                if isinstance(tv, list):
                    tv = tv[0]
                if isinstance(tv, arb):
                    ce = tv.exp() if hasattr(tv, 'exp') else arb(tv).__exp__()
                    dm = self.dist_max[r].v
                    dmv = dm[t] if not isinstance(dm[t], G) else dm[t].v
                    out.append(np.array(ce - dmv).reshape(-1)[0])
                else:
                    constr = np.exp(tv - 1.0) - self.dist_max[r].v[t]
                    out.append(constr)
        return out


# ------------------------------------------------------------------ workspace

class Workspace:
    def __init__(self, pm, q, K, omega, max_level=MAX_LEVEL):
        self.pm = pm
        self.q = q
        self.K = K
        self.max_level = max_level
        self.terms = [[] for _ in range(max_level)]
        self.num_retain_comp_id = {}
        for r in range(6):
            for lv in range(2, max_level + 1):
                self.num_retain_comp_id[(r, lv)] = pm.register(1, 0.0, np.inf)
        self.num_retain_glob_id = [pm.register(1, 0.0, np.inf) for _ in range(6)]
        self.single_mat_size_id = pm.register(1, 0.0, np.inf)
        self.omega_id = pm.register(1, 0.0, np.inf)
        self.K_id = pm.register(1, K, K)
        self.globstage = GlobalStage(pm, q, max_level, self.terms)

    def omega_value(self):
        return self.pm.get_scalar(self.omega_id)

    def evaluate(self):
        pm = self.pm
        L = self.max_level
        c_viol = []     # values that must be <= 0
        ceq_viol = []   # values that must be = 0
        # init
        self.ret_glob = [G(pm.get(self.num_retain_glob_id[r])) for r in range(6)]
        self.ret_comp = {}
        for r in range(6):
            for lv in range(2, L + 1):
                self.ret_comp[(r, lv)] = G(pm.get(self.num_retain_comp_id[(r, lv)]))
        self.single_mat_size = G(pm.get(self.single_mat_size_id))
        self.omega = G(pm.get(self.omega_id))
        self.globstage.evaluate_init()
        for l in range(1, L + 1):
            for t in self.terms[l - 1]:
                t.evaluate_init() if hasattr(t, "evaluate_init") else None
        # pre
        self.globstage.evaluate_pre()
        for l in range(L, 0, -1):
            for t in self.terms[l - 1]:
                t.evaluate_pre()
        # post
        for l in range(1, L + 1):
            for t in self.terms[l - 1]:
                t.evaluate_post()
        self.globstage.evaluate_post(self.terms)
        # lagrange constraints (level>=3 TermInfo only)
        for l in range(1, L + 1):
            for t in self.terms[l - 1]:
                if isinstance(t, Term):
                    ceq_viol.extend(t.lagrange_constraints())
        ceq_viol.extend(self.globstage.lagrange_constraints())
        # asymmetric hashing level>=3
        for l in range(3, L + 1):
            for r in range(6):
                dimx, dimy, dimz = Dims(r + 1)
                num_block = [0.0, 0.0, 0.0]
                pen = zeros(1)
                pcy = zeros(1)
                pcz = zeros(1)
                for t in self.terms[l - 1]:
                    nb = t.num_block_contribution[r]
                    num_block = [a + b for a, b in zip(num_block, nb)]
                    pen = pen + t.hash_penalty_term[r]
                    pcy = pcy + t.p_compY[r]
                    pcz = pcz + t.p_compZ[r]
                for tt in range(3):
                    if tt == dimx - 1:
                        foo = pen
                        c_viol.append((self.ret_comp[(r, l)] - num_block[tt] + pen).item())
                    elif tt == dimy - 1:
                        foo = pcy
                        c_viol.append((self.ret_comp[(r, l)] - num_block[tt] + pcy).item())
                    else:
                        foo = pcz
                        c_viol.append((self.ret_comp[(r, l)] - num_block[tt] + pcz).item())
        # level-2 symmetric hashing
        num_block = [0.0, 0.0, 0.0]
        for t in self.terms[1]:
            nb = t.num_block
            num_block = [a + b for a, b in zip(num_block, nb)]
        for tt in range(3):
            c_viol.append((self.ret_comp[(0, 2)] - num_block[tt]).item())
        for r in range(1, 6):
            ceq_viol.append(self.ret_comp[(r, 2)].item())
        # global hashing
        for r in range(6):
            dimx, dimy, dimz = Dims(r + 1)
            num_block = self.globstage.num_block[r]
            pen = self.globstage.hash_penalty_term[r]
            pcy = self.globstage.p_compY[r]
            pcz = self.globstage.p_compZ[r]
            for tt in range(3):
                if tt == dimx - 1:
                    c_viol.append((self.ret_glob[r] - num_block[tt] + pen).item())
                elif tt == dimy - 1:
                    c_viol.append((self.ret_glob[r] - num_block[tt] + pcy).item())
                else:
                    c_viol.append((self.ret_glob[r] - num_block[tt] + pcz).item())
        # single matrix size
        ms = self.globstage.mat_size
        for tt in range(3):
            denom = 1.0
            if tt == 2:
                denom = self.K
            c_viol.append((self.single_mat_size - ms[tt] / denom).item())
        # final value
        value = (self.single_mat_size * self.omega).item()
        for r in range(6):
            value += self.ret_glob[r].item()
        for lv in range(2, L + 1):
            for r in range(6):
                value += self.ret_comp[(r, lv)].item()
        self.value = value
        c_viol = [float(np.asarray(x).reshape(-1)[0]) for x in c_viol]
        ceq_viol = [float(np.asarray(x).reshape(-1)[0]) for x in ceq_viol]
        c_viol.append(float(np.log(self.q + 2) * (2 ** (L - 1)) - value))
        return np.array(c_viol), np.array(ceq_viol), value


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
    with np.errstate(invalid="ignore"):
        cmax = np.nanmax(c) if len(c) else 0.0
        cemax = np.nanmax(np.abs(ceq)) if len(ceq) else 0.0
    return cmax, cemax, value, c, ceq


if __name__ == "__main__":
    import sys
    mat = sys.argv[1] if len(sys.argv) > 1 else \
        "../scratch/alman_code/data/W1.00_2.371339.mat"
    pm, ws, params = build_and_load(mat)
    print(f"registered params: {pm.num_input}  file params: {len(params)}")
    print(f"loaded omega in file: {ws.omega_value():.10f}")
    cmax, cemax, value, c, ceq = max_violation(pm, ws)
    print(f"max c violation      : {cmax:.6e}")
    print(f"max |ceq| violation  : {cemax:.6e}")
    print(f"value (Schonhage LHS): {float(value):.10f}")
    print(f"target 2^(L-1) log(q+2) = {np.log(7) * 4:.10f}")
    print(f"GetFeasibility maxViol (c, |ceq| combined): "
          f"{max(cmax, cemax):.6e}")
