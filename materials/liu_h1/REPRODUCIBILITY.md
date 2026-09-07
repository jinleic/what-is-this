# Liu Hypothesis 1 reproducibility

## Claim boundary

The continuum theorem is an ordinary **HUMAN-AUDITED** proof. The commands below reproduce:

1. **MACHINE-VERIFIED** exact algebraic identities;
2. **MACHINE-VERIFIED** finite Arb compression and finite-grid claims;
3. **COMPUTATIONAL-EVIDENCE** floating-point spectra and fits.

They do not formalize Taylor's theorem, logarithmic differentiation, convergence of the Gram series, Fubini, positive-kernel closure, or extension to finite signed Borel measures. Those steps are audited in [`AUDIT.md`](AUDIT.md). Hypothesis 2 is not proved by this H1 package; the separate current [H2 argument](../uc/H2_PAPER/REFEREE.md) has its own review obligations.

## Exact environment

The clean reproduction on 2026-08-26 used an isolated `uv` environment, CPython 3.14.3, and these exact packages:

| Component | Version |
|---|---|
| python-flint | 0.9.0 |
| mpmath | 1.3.0 |
| NumPy | 2.5.2 |
| SciPy | 1.18.0 |
| SymPy | 1.14.0 |
| platform | `macOS-26.5.2-arm64-arm-64bit-Mach-O` |

The machine-readable capture is [`verification/environment.json`](verification/environment.json). Every numerical thread pool was limited to one thread:

```sh
export PYTHONHASHSEED=0 PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
```

## Clean-room commands

Run from `math/`:

```bash
UV=(uv run --isolated --no-project --python 3.14 \
  --with python-flint==0.9.0 \
  --with mpmath==1.3.0 \
  --with numpy==2.5.2 \
  --with scipy==1.18.0 \
  --with sympy==1.14.0 \
  python -I -B)

"${UV[@]}" uc/liu_R_nsd.py
"${UV[@]}" uc/liu_kernel.py
"${UV[@]}" uc/liu_tail.py
"${UV[@]}" uc/liu_moment_cert.py
```

The isolated runner may reuse downloaded wheel archives, but it creates an environment independent of the repository `.venv`; Python `-I -B` disables user-site, environment-path, and bytecode-cache dependence.

Run the independent standard-library checker separately:

```sh
python3 -I -B LIU_H1/verification/independent_exact_checker.py
```

It imports no canonical `uc` module and no SymPy, FLINT, NumPy, SciPy, or mpmath code.

### Scaled-family extension — 2026-09-06

The current manuscript explicitly extends the Gram proof to every
`f_lambda(s)=lambda*s*(1-s)`, `0 <= lambda <= 1`. From `math/`:

```sh
nice -n 10 .venv/bin/python -I -B LIU_H1/verification/scaled_kernel_check.py
```

The new independent SymPy 1.14.0 calculation passed seven exact identities and
the lost-square negative control (`-1/219`), under a 15-second CPU / 30-second
wall cap. It imports no canonical proof code. This checks the new algebra;
the continuum conclusion is the manuscript's ordinary proof, not a formalization.
The unchanged standard-library checker above also passed its 16 identities.
Current source/PDF/checker pins are in the
[`2026-09-06 evidence record`](../../data/portfolio-supervisor-20260906-h1.json).
The earlier source and log hash table below is historical and does not pin the
now-extended manuscript.

## Observed clean results — historical 2026-08-26

All five commands exited `0`. Supervised wall times on this workstation were:

| Command | Wall time | Evidence log |
|---|---:|---|
| `uc/liu_R_nsd.py` | 52.6 s | [`liu_R_nsd.clean.log`](verification/logs/liu_R_nsd.clean.log) |
| `uc/liu_kernel.py` | 2.0 s | [`liu_kernel.clean.log`](verification/logs/liu_kernel.clean.log) |
| `uc/liu_tail.py` | 2.1 s | [`liu_tail.clean.log`](verification/logs/liu_tail.clean.log) |
| `uc/liu_moment_cert.py` | 1.3 s | [`liu_moment_cert.clean.log`](verification/logs/liu_moment_cert.clean.log) |
| independent exact checker | 0.03 s | [`independent_exact_checker.clean.log`](verification/logs/independent_exact_checker.clean.log) |

Times include isolated-runner startup with a warm local package cache and are not performance guarantees.

The canonical proof run ended with:

```text
HUMAN-AUDITED [integral Gram proof]: R <= 0 without projection.
HUMAN-AUDITED [restricted-form identity]: Liu's Hypothesis 1 is settled.
OPEN [scope]: Liu's Hypothesis 2 ... remains open.
MACHINE-VERIFIED [finite Arb]: degree-8, degree-12, and degree-16
   compression upper bounds are strictly negative corroborating checks.
```

The independent checker printed 16 `MACHINE-VERIFIED` identities and explicitly left the calculus and measure theory human-audited.

## Pinned source and log hashes

SHA-256 values for the reproduced inputs and outputs:

| Path | SHA-256 |
|---|---|
| `uc/liu_R_nsd.py` | `5e6197f18132cd8e139786bcbd51afd0c85b011a832ec71f895ebaea3024211c` |
| `uc/liu_kernel.py` | `de7e767b3a3f6c4f077efe04ff5e380fff4e19f63d1d57316c0ea6a6b2ff3529` |
| `uc/liu_tail.py` | `2eb61f9891feb080a82cd2457c05260f80592d2600b29f76b8e0faf7111218b4` |
| `uc/liu_moment_cert.py` | `f0f37d46d37a9d781b2e3fab12519d524968852357b7f0e629a32c61a79c4ce6` |
| `LIU_H1/verification/independent_exact_checker.py` | `1a7daa605359dea2629d0054a00f252fb85203324ea80e02ccbd1c6488491ab6` |
| `verification/logs/liu_R_nsd.clean.log` | `df6de77d661587cb654e8ef015b094bbbdbadb4784a22dc7d5682bcafb641d75` |
| `verification/logs/liu_kernel.clean.log` | `b16776b215475a75eb79c34d0fd3dc9e712968f0f6069769ee9e08105be36109` |
| `verification/logs/liu_tail.clean.log` | `759ac700858028328b2a0a2bb5b48632733db995f5ce44552b6f748923c3e7ee` |
| `verification/logs/liu_moment_cert.clean.log` | `f9669aee235aaeed2c8ba969aff06bf76dba58cb34204427ba3f26071627a8cf` |
| `verification/logs/independent_exact_checker.clean.log` | `c157a5316776292144f122813fac9b8dcc05bf61a1322069f6d9fa78824f874a` |
| `paper/main.tex` | `22a5b0483f00c6efebff548f268c84897a8d08c246c495bc023a95a67dcad49b` |
| `paper/refs.bib` | `674a263e3cc980bf3a0c46f40dd3e1b4987ef5056133c378edcac3ba590d2118` |
| `paper/main.pdf` | `180dc4773fd516eab11fd90f79c7893744b20c04bc7db73fb12618213edff907` |
| `verification/logs/codex_referee_2026-08-26.md` | `8b22810a5505882557a4c4e0f548d51b2a185177fa361e36768ee0f8512795e5` |
| `verification/logs/internal_referee_2026-08-26.md` | `de6d82a28f8522c8c9e600972636c54027ebfd7007d464827f18853622f93ebe` |
| `verification/logs/fable_referee_attempt.log` | `6af2b40439da28045b20609ac296755356bf504a11da2b29266871ef59c1cf7e` |
| `literature/search/arxiv_refresh_2026-08-26.json` | `48a0461b2dd1da44d5650da529b9f6510576e5da8576645f4dacea22888dfe69` |
| `literature/manifest.sha256` | `24d6ed4151d768cda12ac17acb8f15218da9b261ffabd4c14a8b8aab572f974e` |

A changed source or log must be re-run and re-hashed; these values are evidence for the exact bytes audited here.

## Manuscript reproduction

From `math/LIU_H1/paper/`:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

A successful build checks only manuscript syntax, references, and bibliography. It does not verify the theorem. The theorem-level literature search, saved primary sources, and novelty limits are in [`LITERATURE_ORIGINALITY.md`](LITERATURE_ORIGINALITY.md).

The repaired manuscript was rebuilt after two fresh referee passes.  The
independent Codex review and separate internal re-derivation found no blocker
or major defect; their four minor findings are resolved and preserved in the
two review logs above.  The Fable log records an availability failure before
any mathematical review and contributes no positive evidence.

## Acceptance checklist

- All exact checker assertions pass with exit code `0`.
- The canonical output uses `HUMAN-AUDITED`, `MACHINE-VERIFIED`, `NUMERICAL`, and `CONJECTURED` consistently; its final scope block separately marks Hypothesis 2 `OPEN`.
- Finite Arb results are reported only for their named finite polynomial subspaces or node sets.
- No finite spectrum, sampled ratio, or roundoff-scale eigenvalue is promoted to a continuum theorem.
- The manuscript states the finite-signed-measure theorem first and treats the `L^2` operator statement as a corollary.
- The H1 theorem and its scaled corollary alone imply no unconditional frequency improvement; do not cite them as a proof of Hypothesis 2.
