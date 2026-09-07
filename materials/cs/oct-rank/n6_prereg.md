# N6 pre-statement — gate `n6-qi-real13-descent`

Campaign: mint only after this file is committed, with
`python3 ../../scripts/campaign.py init --gate n6-qi-real13-descent --prereg n6_prereg.md --agent OctRankN5`
from `cs/oct-rank/`. Date: 2026-09-04. No N6 scientific computation may run
before init.

## 1. Question and why this is the shortest new frontier

Target `TF` is the real 3 x 8 x 8 tensor `blockdiag(tau,tau)`, rank-equivalent
to the `(L_1,L_i,L_j)` octonion triple under the frozen Route-F identities.
The standing exact interval is `rank_R(TF) in {13,14}`. This gate attacks only
the downward side: produce a certified real rank-13 decomposition.

This is a deterministic Galois-symmetric continuation from N3's exact Q(i)
rank-13 point to the real locus. It is not a rerun of N1: there are no random,
tau-merge, frozen-candidate, or TRF seeds and no least-squares seed sweep. It
is not N2: there is no Groebner/F4/elimination computation. Ordinary
Nullstellensatz infeasibility is structurally unavailable because N3 gives an
exact complex point. Degree-4/5 SOS is only RAM-unexcluded, not presently
runnable: the repository has no SDP generator, solver, rational recovery, or
exact replay for its 30,876-dimensional Gram matrix. The exact Q(i) point is
the only saved zero-residual initialization not already exhausted.

Frozen inputs, read-only and hash-pinned:

* N3 exact twelve-term witness
  `campaigns/20260904T034332Z_f5a61843_0333458ed434/witness_tf_c12_qi.json`,
  SHA-256
  `251e00577b43b400aa97d306bb923a98eb79ca7850ac67652d45288c2e570037`;
* N3 instrument (the registered T5 split specification), SHA-256
  `d529cce2938af06645529f624c58a877a9fc3c45f6bc1efc8affafa16cca94ef`;
* frozen real lower-13 replay, SHA-256
  `894b09096763bbdc88c38c9d315d546d4826b161849d06c4be8ba138afaaa590`;
* validated square-slice Krawczyk source, copied into the new run by value,
  source SHA-256
  `781e51ec1a72cfb8acf7e7c7b1e44a96cac011ec6f46dc61e88f5b7e01b87035`.

No existing saved real candidate is a witness. N1's best vector was not
serialized and recovering it would require the prohibited rerun.

## 2. Exact starting point and controls

Rebuild TF independently from the Hamilton table. Parse the N3 witness over
Q(i), require its hash, and substitute its twelve terms into all 192 entries.
Reconstruct N3's registered 13-term T5 point without search: split term zero
in the conjugated B-factor by fixed `W=e_3`, using the serialized S0 only to
map that fixed split vector back by `blockdiag(S0^-1,S0^-1)`. This is a
reconstruction of the frozen point, not a new S0 choice. Require exactly 13
nonzero terms and 192/192 exact substitution. Its coefficientwise conjugate
must also pass 192/192. Both must reject `TF[0,0,0]+=1` with the mismatch
recorded.

Before TF tracking, run a deterministic positive plant through the same
continuation code. Generate integer factor matrices of shapes 3x13, 8x13,
8x13 with NumPy PCG64 seed `20260904`, entries sampled uniformly from
`{-3,...,3}`; regenerate the entire fixed plant once with seed `20260905`
only if a factor column is zero or the initial complex Jacobian chart fails
(the two-seed list and first-valid rule are fixed now). Form its exact real
tensor, then disguise the point termwise over Q(i) by multiplying A-column s
by `g_s` and B-column s by `g_s^-1`, with `g_s=1+i` for even s and `g_s=2+i`
for odd s. The same two-sided continuation must meet at a real midpoint and
reproduce the plant. The exact real plant at its true solution must be
accepted by the copied square-slice Krawczyk core. Its one-coordinate-corrupt
target must fail substitution by that same witness. Any failed positive
control or false substitution accept is INVALID INSTRUMENT; TF is not run.

The frozen Krawczyk core is consumed as code, not as a theorem claim. Its
source bytes are copied and hash-checked before import. A candidate
certificate must be independently replayed by a separate `n6_replay.py`.

## 3. Fixed continuation

Pack factors in canonical order A (3x13), B (8x13), C (8x13), total 247
complex coordinates, with the 192 CP residuals in lexicographic `(p,b,c)`
order.

At the exact start `x0`, form the complex 192x247 Jacobian. Use SciPy QR with
column pivoting (SciPy 1.18.1, NumPy 2.5.2) only to propose the first 192
columns. Verify the proposed 192x192 complex chart is nonsingular exactly by
the 384x384 rational realification `[[Re,-Im],[Im,Re]]` using python-flint
0.9.0. If it is singular, stop INCONCLUSIVE-CHART; there is no pivot fallback.
The remaining 55 coordinates are free.

For each free coordinate prescribe the conjugation-symmetric straight path

`q(t)=(1-t)q0+t*conj(q0)`, `0<=t<=1`.

Track the 192 pivot coordinates from both exact endpoints to `t=1/2`: forward
from `(t,x)=(0,x0)` and backward from `(1,conj(x0))`. Each half has 64 fixed
base steps. At each target step use at most 8 complex Newton corrections on
the fixed square chart. If a step fails, bisect it deterministically to depth
at most 8; no other damping, seed, chart, or parameter change is allowed. A
Newton step is accepted only when max complex residual is at most `1e-12`.
Checkpoints are atomically persisted after every accepted base step.

A TF midpoint is admitted to certification only if both halves reach 1/2,
their max coordinate disagreement is at most `1e-10`, the averaged midpoint's
max imaginary coordinate is at most `1e-10`, its real CP relative residual is
at most `1e-12`, all factor-column norms are nonzero and at most `1e4`, and
the real 192x247 Jacobian has numerical rank 192 at threshold `1e-9`. These
are admission filters only and prove no rank statement.

## 4. Certification and replay

For an admitted midpoint, first convert every float64 coordinate to its exact
dyadic rational and test all 192 tensor equations exactly. If they vanish,
serialize the exact real 13-term witness. Otherwise pass the dyadic center to
the copied Route-F square-slice Krawczyk core with its fixed QRCP slice and
frozen radius ladder `1,10^-1,...,10^-12`. A rank claim requires strict
containment on all 192 solved coordinates and the core's outward-Arb check on
every strict margin. Serialize the center, selected/free columns, exact
bounds, and radius.

`n6_replay.py`, which does not import the continuation driver, must rebuild TF,
re-evaluate the exact center/residual/Jacobian and rerun all containment
margins from the serialized data. Only a clean replay is adjudicative.
Combined with the frozen lower-13 replay, either an exact real substitution or
a strict replayed Krawczyk containment yields `rank_R(TF)=13`
MACHINE-VERIFIED modulo the already frozen lower-13 chain labels.

Any path miss, chart failure, non-real midpoint, threshold miss, no
containment, replay failure, timeout, or crash is FAILURE TO CERTIFY. It is
not evidence for rank 14. “Absence of a 13-witness is NOT evidence for 14;
absence of an impossibility argument is NOT evidence for 13.” The interval
remains `{13,14}`, `18 <= R_R(T_O) <= 25` is untouched, and no published work
is refuted.

## 5. Runtime and lifecycle

Python 3.14.3, NumPy 2.5.2, SciPy 1.18.1, python-flint 0.9.0. Set
`PYTHONDONTWRITEBYTECODE=1`, `sys.dont_write_bytecode=True`, nice 10, and all
BLAS/thread variables to one before importing NumPy. Set equal soft/hard
`RLIMIT_CPU=1800` seconds. Enforce a 4 GiB RSS ceiling at every Newton/base
boundary from `ru_maxrss`, with an orderly persisted abort. One driver
process; the independent replay may be a second finite process only after an
admitted candidate. No random choice beyond the two fixed plant seeds, no
post-result adjustment, no TRF, no Groebner/F4, no SOS, and no writes to any
frozen run.

Commit this preregistration, init, copy it byte-identically, bind all source
hashes in provenance, create the campaign-local driver/core/replay, and launch
under process supervision as:

```
PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
nice -n 10 /Users/jinleic/jinleic-workspace/cs/.venv/bin/python \
  n6_qi_real13_descent.py
```

Freeze and close only after terminal results and any required replay. Preserve
all checkpoints, failed attempts, and defect logs.
