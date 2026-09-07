# Referee package — the union-closed constant, unconditionally, and improved

Two results are submitted for review. Both are machine-assisted; every non-elementary
step is a pinned artifact, and the whole chain replays from source in about ten minutes.

| | |
|---|---|
| Manuscript | [`main.tex`](main.tex), 25 pp. |
| Replay driver | [`../verification/replay_h2_chain.py`](../verification/replay_h2_chain.py) |
| Artifacts | eight JSON certificates under `uc/verification/results/` |
| Package generator | [`make_referee_package.py`](make_referee_package.py) |

The two top-level hashes are printed by the generator and listed in `MANIFEST.tsv`; they
are deliberately **not** hand-copied into this file, so that this file cannot go stale.
Run `./.venv/bin/python -I -B uc/H2_PAPER/make_referee_package.py` from `math/` to obtain
`dist/referee-package/MANIFEST.tsv` and `dist/referee-package.tar.gz`.

**Novelty gate update (2026-09-05):** the reported Costa preprint in the
user's research history has not yet been located and compared. See the
[targeted unresolved lead](../LITERATURE_ORIGINALITY.md#targeted-unresolved-prior-art-lead--2026-09-05).
The dated literature comparison in §2 is not a priority clearance. The package
remains available for proof review; neither outside refereeing nor submission
has occurred.

## 1. What is claimed

**Result 1 (unconditional $c'$).** Every union-closed family $\mathcal F$ of subsets of a
finite set, other than $\emptyset$ and $\{\emptyset\}$, has an element contained in at
least $c'|\mathcal F|$ of its members, where

$$c' = 1-m^\ast \in [0.382709087918735029930312902098972611626381433 \pm 4.95\times10^{-46}].$$

$m^\ast = p^\ast x^\ast$ with $x^\ast$ the unique root in $(0,1)$ of $z^4-2z^3+3z^2-1$ and
$p^\ast = h(x^\ast)/h(x^{\ast2})$ — Liu's constants, defined by his equations (87)–(90).

**Result 2 (a strictly better constant $c''$).** The same statement holds with

$$c'' = 1-m_{16/25} \in [0.38284565599065407172003132812201407735944023015042 \pm 5.34\times10^{-52}],$$
$$c''-c' \in [0.000136568071919041789718426023041 \pm 4.66\times10^{-34}] > 0.$$

$m_{16/25}$ is the analogous constant for the Example-5 protocol built from
$f(x)=\tfrac45x(1-x)$ instead of $f(x)=x(1-x)$: $x$ is the unique root in $(0,1)$ of
$2x^2+\tfrac{16}{25}x^2(1-x)^2-1$, $p = h(x)/h(x^2)$, $m = px$.

## 2. What is new relative to the literature

Liu, *Improving the Lower Bound for the Union-closed Sets Conjecture via Conditionally IID
Coupling*, arXiv:2306.08824v1 (2023), obtains $\approx0.38271$ but states it conditionally.
His abstract, verbatim:

> Under numerically verified hypotheses, the lower bound for the union-closed sets
> conjecture can be improved to approximately 0.38271, a number that can be defined as
> the solution to an analytic equation.

Those hypotheses are his Section V-A positive-semidefiniteness statement ("Hypothesis 1")
and his Section V-B global-minimiser structure ("Hypothesis 2"). Result 1 removes both:
the inequality behind Hypothesis 2 is proved for *every* pair of Borel probability measures
on $[0,1]$ and, more generally, for every conditionally i.i.d. coupling, with no mean
constraint — which is strong enough that Hypothesis 1 is never needed. Result 2 then
changes the protocol and exceeds the value itself.

Literature position as of 2026-09-03 (bounded search; no priority claim beyond it).
The published progression is Gilmer $0.01$ (arXiv:2211.09055), then
$(3-\sqrt5)/2\approx0.38197$ (arXiv:2211.11689, arXiv:2211.11731, arXiv:2211.11504), then
$\approx0.38234$ from the Sawin/Yu/Cambie convex combination — Liu's own attribution, "the
best constant obtainable through this approach is around 0.38234, as evaluated by Yu and
Cambie" (arXiv:2212.12500; Yu, *Entropy* 25 (2023) 767), certified locally as
$\psi+3\times10^{-4}=0.3823455$ — and then Liu's conditional $\approx0.38271$. That
ceiling is **0.38234, not 0.38284**: it is not a near-collision with $c''$, and a
search-engine synthesis conflating the two decimals was observed and is recorded as false.
The nearest published decimal to $c''$ is Liu's conditional $0.38271$.

The record is the arXiv API, not a search summary. The `all:"union-closed"` endpoint of
`uc/literature/search_protocol.json` step 1 was re-queried on 2026-09-03
(`submittedDate` descending) together with the pinned-record `id_list` query; the machine
record is `uc/literature/arxiv_refresh_2026-09-03.json` and the audit is
`uc/LITERATURE_ORIGINALITY.md`, which is the single authority for originality here. The
feed reports 104 records (103 at the prior 2026-08-26 refresh). The one new record,
arXiv:2608.25147v1 (Tian, *Frankl's Conjecture at Height Four and the Structure of
Height-Five Counterexamples*), proves Frankl for height $\le4$ (empty-set-free form) and
constrains height-five counterexamples; it proves no universal frequency constant, and its
full text contains no occurrence of *entropy*, *Gilmer*, *coupling*, or the string `0.38`.
arXiv:2306.08824 is still v1, so no revision has removed the conditionality, and the other
pinned records are unchanged. We found no published treatment of the scaled Example-5
family $f_\lambda(x)=\lambda x(1-x)$.

## 3. Replay in one command

From `math/`, with the repository's virtual environment:

```
nice -n 19 ./.venv/bin/python -I -B uc/verification/replay_h2_chain.py
```

Expected: exit status 0, final two lines

```
checks passed: 155  failed: 0
REPLAY_H2_CHAIN PASS
```

in roughly 600 s on one core (dominated by `psi-reduction`, ~300 s per run). The driver is
standard-library only. For each of the eight links it

* **[A]** requires each pinned SHA-256 (artifact file, internal digest, module source) to
  occur verbatim in `PROGRESS.md`, so the ledger and the code cannot drift apart silently;
* **[B]** recomputes those three hashes from disk, recomputes each artifact's internal
  digest in that module's own canonicalisation, and checks the cross-artifact pins;
* **[C]** re-runs the module twice from `math/` under the invoking interpreter with
  `-I -B` and requires both runs to reproduce the on-disk artifact byte for byte.

No run overwrites a pinned artifact: seven modules accept `--output` and write to a scratch
directory; the one fixed-path writer (`liu9_psi_reduction_audit.py`) is snapshotted before
each run and restored afterwards, including on interrupt.

**Environment.** Python 3.9.6; `python-flint==0.6.0` (Arb ball arithmetic), `sympy==1.14.0`,
`mpmath==1.3.0`, plus `numpy==2.0.2`/`scipy==1.13.1` used only by `liu9_binding.py`; see
`requirements-freeze.txt`. Byte-identical replay depends on those versions. A POSIX `nice`
and a venv at `math/.venv` are required, because `liu9_h2_phi_audit.py` shells out to
re-run the boundary module. Precisions: 400 bits for the twovar, boundary, phi-audit, lift,
mixture and scaled certificates; 480 for the reduction; 320 for its independent
re-derivation.

## 4. Independent-verification checklist

Executable without trusting any tooling in this repository.

1. **File hashes.** From the package root, `shasum -a 256 -c SHA256SUMS`. This covers all
   eight artifacts, all eleven modules, the driver, the manuscript and the ledger.
2. **Internal digests.** For each artifact, remove or blank the `report_sha256` field, then
   SHA-256 of `json.dumps(body, sort_keys=True, separators=(',',':'))` in the artifact's own
   scope. The scope per artifact is the `digest_scope` column of the driver's `CHAIN` table
   and is restated in `MANIFEST.tsv`: `omit` for `liu9-h2-reduction`, `empty` (set the field
   to `""`) for `liu9-psi-reduction`, `omit+nl` (append `"\n"`) for the other six.
   **Caveat, and it will cost you time if you miss it:** `liu9-h2-twovar.json`,
   `liu9-h2-boundary.json` and `liu9-h2-phi-audit.json` describe their own scope in prose as
   "with `report_sha256` omitted", without mentioning the trailing newline — but their pinned
   digests do include it. The three modules' canonicalisation is `omit+nl`; the prose is
   incomplete, the digests are correct. It is not fixed here because changing those strings
   would regenerate the artifacts and invalidate their pins throughout the chain.
3. **Ledger presence.** `grep -c <hash> PROGRESS.md` is at least 1 for each of the 24 pinned
   values. That is check [A], and it needs only `grep`.
4. **One module, without the driver.** For any `--output` module:
   `./.venv/bin/python -I -B uc/<module>.py --output /tmp/a.json`, twice, then `cmp` both
   against the pinned artifact. For `liu9_psi_reduction_audit.py`, copy the tree first: it
   writes its fixed path.
5. **Cross-artifact pins.** In the `dependencies` block of `liu9-h2-general-lift.json`,
   `liu9-h2-mixture-theorem.json` and `liu9-cprime-four-fifths-ab.json`, every entry records
   `match: true` and the hash it expects; eleven of those entries name a chain link and are
   compared against that link's pin by check [B]. The two that do not
   (`frontier_module`, `frontier_artifact`) are outside the chain and are asserted by the
   eighth module itself at run time, so [C] covers them.
6. **Cheap numeric spot-checks**, in any interval-arithmetic library or CAS.
   * $A(\tfrac12,\tfrac12)=R(\tfrac12,\tfrac12)\in[0.0067704710\pm2.79\times10^{-11}]>0$.
   * The interior zero: $A(x^\ast,x^\ast)$ and $\nabla A(x^\ast,x^\ast)$ enclose $0$ to
     radius $\le5\times10^{-68}$; at the scaled protocol, $A(x,x)\in[\pm1.19\times10^{-68}]$.
     This is why bisection alone cannot prove $A\ge0$ and the zeros are handled by exact
     algebra plus a Taylor bound.
   * $\kappa_{\mathrm{master}}$ minima, which show the covers are not over-tightened:
     $\kappa_{\mathrm{master}}(\tfrac{27}{50},\tfrac{27}{50}+10^{-6})\in[2.394193379778020415966765\pm3.03\times10^{-25}]$
     at $\kappa=1$, and
     $\kappa_{\mathrm{master}}(\tfrac{501}{1000},\tfrac{501001}{1000000})\in[2.275086765014144819939496\pm4.33\times10^{-25}]$
     at $\kappa=\tfrac{16}{25}$.
   * The $(1+\kappa ab)$ factor is load-bearing: $|\gamma'(\tfrac9{20})|\in[1.910518707203791212162174\pm3.80\times10^{-25}]<2$
     at $\kappa=1$, and $|\gamma'(\tfrac{54}{125})|\in[1.926050003455288550151689\pm2.56\times10^{-25}]<2$ at $\kappa=\tfrac{16}{25}$.
   * Endpoint bounds at $\kappa=\tfrac{16}{25}$: $\nu_0(\tfrac18)\in[2.9283853862304257001\pm3.12\times10^{-20}]$,
     $\nu_1(\tfrac78)\in[4.3347701765338265363\pm2.86\times10^{-20}]$, both below the exact
     $|\gamma'|$ there.
   * Polynomial identities, exactly, in any CAS: $x^2+\pi_\kappa(x,x)-1=2x^2+\kappa x^2(1-x)^2-1$,
     which is $x^4-2x^3+3x^2-1$ at $\kappa=1$;
     $(1+\kappa a^2)(1+\kappa b^2)-(1+\kappa ab)^2=\kappa(a-b)^2$; and both rescalings of
     Lemma 10.6.
7. **Falsification hooks.** The certificates are not vacuous: each module carries
   deliberate errors ("mutations") and refuses to certify while any is active. **54 across
   the chain** — 4 reduction, 4 psi-reduction, 10 twovar, 5 boundary, 6 phi-audit, 7
   general-lift, 6 mixture-theorem, 12 scaled (5 for $A$, 5 for $\Phi$, 2 for the mixture
   identity). Each is recorded in its artifact with the witness that refutes it. Two are
   worth reading first: a first version of the lift artifact declared 100-digit rational
   truncations of the constants to be "the constants", and at those rationals Liu's own
   minimiser has a *certified negative* gap $[-5.95\times10^{-102}\pm10^{-112}]$ — the
   theorem is a statement about the exact algebraic constants, not about nearby rationals;
   and the scaled N0 mutation claims a lower bound $[4.802552033417898148112837\pm1.3\times10^{-28}]$
   that exceeds the exact $|\gamma'(\tfrac18)|\in[3.249243440498020263844499\pm2.23\times10^{-25}]$,
   hence is certified false.

## 5. Limits, stated plainly

* **Neither result has been refereed outside this repository.** Three blank-context
  adversarial reviews are recorded (the $\Phi$ certificate; the lift/mixture artifacts and
  the Proposition-3 argument; the scaled certificate and the write-up), and they found no
  soundness defect in the final versions — but they are internal records, not refereeing.
* The results concern **the constant only** and say nothing about the conjectured value
  $\tfrac12$.
* **No optimality in $\lambda$ is claimed.** $\lambda=\tfrac45$ is one value where the
  endpoint constant exceeds $c'$ and the necessary face condition $c_0>0$ still holds. Each
  $\lambda$ needs its own root, its own $\beta$, and its own strata constants.
* $c'$ and $c''$ are each the **best constant their own protocol pair can give** (the
  inequality is an equality along a one-parameter family of laws), so any further
  improvement requires another protocol, not sharper analysis of these.
* Items labelled **computational-evidence** carry no proof weight and are used nowhere:
  the numerical minimum of $P_2$, the numerical minima of $\kappa_{\mathrm{master}}$, the
  random-pair samples, and the float discovery scans.
* Three open obligations are recorded and deliberately not attempted: the exact value of
  $m_{\rm full}(\beta)$ off $\beta^\ast$; the Liu-Lemma-8 / Sawin three-protocol route,
  which needs $\inf_{\mathcal C_2(\mu)}\mathbb E[g_2]$ controlled uniformly over all
  symmetric self-couplings; and $\lambda$-optimality.
* Superseded material inside the repository is **not** part of this package, and one item
  must not be cited as a live obligation: the entrywise target $\Psi_m\ge0$ on $[0,1]^4$ is
  refuted at $(0,\tfrac12,\tfrac12,0)$ and was replaced by the reduction actually used.

## 6. Package contents and provenance

`MANIFEST.tsv` lists every file with its role, size and SHA-256. The set is closed under
three relations, all re-checked by `make_referee_package.py --check` before anything is
written: every repo path the manuscript names via `\path{...}`; the eight `CHAIN` links of
the driver and the 24 hashes they pin; and the transitive local imports of those eight
modules, which reach exactly three further modules and one further artifact.

Three files carry the literature position rather than the mathematics:
`uc/LITERATURE_ORIGINALITY.md` (the audit, and the single authority for originality),
`uc/literature/search_protocol.json` (the endpoints queried and their coverage limits),
and `uc/literature/arxiv_refresh_2026-09-03.json` (the API record behind §2).

One disambiguation, because the repository contains a similarly named document:
`uc/REPRODUCIBILITY.md` is **not** about this chain. It records the clean-room replay of
Campaign I — a different target on the $\psi$-route, with its own 488-million-node trace
audit — and mentions none of the eight artifacts here. It is deliberately outside this
package; §3 above is the only replay recipe for these results.

`PROGRESS.md` is the repository ledger and is included in full. It is large and mostly
concerns unrelated work, but check [A] asserts the 24 pins against it, so shipping it is
what makes the anti-drift check runnable. The pins live in two blocks: the `H2-REFEREE-PACKAGE`
entry (21 values, the $c'$ chain) and the newest `H2-PAPER-CPRIME2` entry (the eighth link).
Where a hash was superseded, the ledger keeps the old value inside an explicit
"was `<hash>`" or "superseded" note rather than deleting it.

Nothing in this directory is a stored copy of a canonical file: the package is generated,
so a change anywhere upstream appears as a changed hash on the next run rather than as
silent drift.
