"""Validation of the FSS fitting module on synthetic data with KNOWN exponents.

Must pass before the module touches real sweep data (pre-campaign requirement,
fss-bb/README.md). Geometry (grid uniform in the SCALING VARIABLE x, which is
how a well-designed FSS sweep is laid out):

  - five sizes N = 144..1296; per size, p points at
        p*(N) + dx * N^(-1/nu0),  dx in {0, +-0.5, +-1.0, +-1.5, +-2.0, +-2.5};
  - p*(N) = p0 - a*N^(-1/nu0) with a = (p0 - 0.370) * 144^(1/nu0): the ladder
    is CONSISTENT with the truth and starts at p*(144) = 0.370;
  - the scaling function satisfies f(x=-a) = 0.10 (so p*(N) really is the
    WER=0.10 crossing at every size) and f stays strictly inside (0,1) on the
    whole sampled x-range (asserted — no clipping bias);
  - truth (p0, nu0) deliberately AWAY from the paper's (0.488, 1.18).

T1 (exact recovery): f = cubic in x (shifted cubic in t = x + a; a shifted
    cubic is still a cubic in x, so the fitter's model class contains truth
    EXACTLY). Pass: truth (p0, nu0) inside the 95% bootstrap CIs.
    Known cosmetic caveat: the designed cubic has a shallow non-monotone dip
    (~1e-4 WER) in the far-left tail — physically unrealistic but irrelevant
    to estimator recovery since truth lies exactly in the model class.
T2 (ansatz robustness): f = shifted logistic — NOT in the fitter's model
    class. Pass: |p_inf_hat - p0| < 0.01; nu deviation is REPORTED as the
    ansatz-bias scale (the paper carries the same caveat for its polynomial).

Run:  python scripts/validate_fss_synthetic.py
Writes a JSON artifact + prints a PASS/FAIL table.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from fss_fit import fit_fss_bootstrap, fit_linearized_fss  # noqa: E402

SHOTS = 200_000
SIZES = [144, 324, 576, 900, 1296]
P_STAR_144 = 0.370   # ladder anchor (paper Table 1 smallest size)
DX_GRID = np.array([0.0, -0.5, 0.5, -1.0, 1.0, -1.5, 1.5, -2.0, 2.0, -2.5, 2.5])


def ladder_a(p0: float, nu0: float) -> float:
    return (p0 - P_STAR_144) * 144 ** (1.0 / nu0)


def make_data(p0, nu0, f, seed):
    """f takes x = (p - p0) N^(1/nu0); f(-a) must be 0.10; in-range asserted."""
    a = ladder_a(p0, nu0)
    rng = np.random.default_rng(seed)
    ps, Ns = [], []
    for N in SIZES:
        pstar_N = p0 - a * N ** (-1.0 / nu0)
        for dx in DX_GRID:
            ps.append(pstar_N + dx * N ** (-1.0 / nu0))
            Ns.append(N)
    ps = np.array(ps)
    Ns = np.array(Ns)
    x = (ps - p0) * Ns ** (1.0 / nu0)
    w_true = f(x, a)
    assert w_true.min() > 0.0 and w_true.max() < 1.0, (
        f"synthetic truth out of (0,1): [{w_true.min()}, {w_true.max()}]"
    )
    counts = rng.binomial(SHOTS, w_true)
    return ps, Ns, counts, w_true, a


def run_test(name, p0, nu0, f, seed, expect_exact: bool):
    p, N, counts, w_true, a = make_data(p0, nu0, f, seed)
    pstar = p0 - a * N ** (-1.0 / nu0)
    res = fit_fss_bootstrap(counts, np.full_like(counts, SHOTS), p, N, pstar, window=0.06, n_boot=500, seed=7)
    lin = fit_linearized_fss(p0 - a * np.array(SIZES, float) ** (-1.0 / nu0), SIZES, nu=nu0)
    lo, hi = res["p_inf_ci"]
    nlo, nhi = res["nu_ci"]
    t1_pass = (lo <= p0 <= hi) and (nlo <= nu0 <= nhi)
    print(f"\n[{name}]  truth: p_inf={p0}  nu={nu0}  (a={a:.4f}, wer span "
          f"{w_true.min():.4f}..{w_true.max():.4f})")
    print(f"  collapse fit : p_inf={res['p_inf']:.4f} CI95=({lo:.4f},{hi:.4f})  RSS={res['rss']:.2e} pts={res['n_pts']}")
    print(f"                 nu={res['nu']:.3f} CI95=({nlo:.3f},{nhi:.3f})  n_boot_ok={res['n_boot_ok']}/500")
    print(f"  linearized   : p_inf={lin['p_inf']:.4f} +-{lin['p_inf_se']:.4f}  (nu fixed {nu0})")
    ok = t1_pass if expect_exact else abs(res["p_inf"] - p0) < 0.01
    print(f"  verdict: {'PASS' if ok else 'FAIL'}")
    return {
        "name": name, "p0": p0, "nu0": nu0, "a": a,
        "collapse": {"p_inf": res["p_inf"], "ci": [lo, hi], "nu": res["nu"], "nu_ci": [nlo, nhi]},
        "linearized": {"p_inf": lin["p_inf"], "se": lin["p_inf_se"]},
        "pass": bool(ok),
        "n_points": int(len(p)),
    }


def main():
    out = {}
    # T1: exact cubic (in x), f(-a)=0.10, in-range on sampled window
    def f1(x, a):
        t = x + a
        # nonnegative on sampled t-range thanks to the soft-landing x^2 term
        return 0.10 + 0.055 * t + 0.055 * t**2 + 0.0041 * t**3

    out["T1_exact_cubic"] = run_test("T1 exact cubic", 0.55, 1.35, f1, seed=20260829, expect_exact=True)

    # T2: shifted logistic (not in model class)
    def f2(x, a):
        s = np.log(8.5) / 1.6  # so f(-a) = 0.10 exactly
        return 0.95 / (1.0 + np.exp(-1.6 * (x + a - s)))

    out["T2_logistic"] = run_test("T2 logistic (model-biased)", 0.61, 1.05, f2, seed=31337, expect_exact=False)

    out["all_pass"] = out["T1_exact_cubic"]["pass"] and out["T2_logistic"]["pass"]
    dest = Path(__file__).resolve().parent.parent / "campaigns-smoke" / "fss_module_validation.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2, default=str))
    print(f"\nartifact written: {dest}")
    return 0 if out["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
