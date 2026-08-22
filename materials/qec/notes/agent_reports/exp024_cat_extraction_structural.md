NEGATIVE_STRUCTURAL — the specified unverified 4+4 two-ancilla scheme cannot recover Gross's seven-layer all-two-qubit round.

# EXP-024 structural result: two-ancilla cat extraction

## Scope and policy

This result concerns catalogue member `12_6_0193` with the exact
catalogue polynomials persisted in the JSON artifact.  Every weight-8
perturbed check uses an **unverified** two-qubit cat: H plus a cat-preparation
CNOT, four controlled-Pauli data couplings from each half, one merge CNOT, and
X measurement of the primary ancilla.  The 72 weight-6 checks remain
one-ancilla measurements.

The verdict metric counts every two-qubit gate layer, including cat
preparation and merge.  Single-qubit preparation/H and measurement are
bookkept separately.  This policy is not fault tolerant: one fault on a cat
ancilla can fan out to four data qubits.  No flag or cat-verification detector
is present.

## Exact and budgeted depth evidence

- Integrated all-two-qubit depth 7: **INFEASIBLE**,
  exact for the specified scheme.  The unrestricted model permits a separate
  arbitrary 4+4 partition for every check and non-translation-invariant data,
  preparation, and merge slots; it enforces physical-qubit disjointness,
  prep < data < merge, and the controlled-Pauli even-crossing criterion.
- Integrated depth 8: **UNKNOWN** after
  the preregistered 60 s
  budget.  This is budgeted inconclusive evidence, not an infeasibility proof.
- Achieved construction: total all-two-qubit depth
  **9**, a valid
  translation-invariant upper bound.
- Data-coupling-only model: depth **7**
  (OPTIMAL).  Excluding cat plumbing creates an
  apparent 8 -> 7 improvement, but that is a bookkeeping artifact and is not
  the verdict depth.
- Integrated optimum bracket for this scheme: **[8, 9]**.

The one-ancilla PBB depth 8 label is class-free exact: weight 8 is a lower
bound and the persisted schedule is a matching witness.  The Gross depth 7
schedule and the EXP-016 schedule populations are translation-invariant; our
local artifact certifies optimality within that class, while the unrestricted
depth-6 exclusion is external (ASC, arXiv:2603.21499).

## Measurement semantics

- Structural schedule valid: `True`.
- Physical-qubit collisions:
  0.
- Controlled-Pauli semantic parity defects:
  0.
- Noiseless sampling: **0**
  detector firings and **0**
  observable flips in 4096 shots.
- Cat and one-ancilla PBB observable definitions identical:
  `True`
  (SHA-256 `ad3e161c65dd4cf5a4733b5fc807e288c019a90ed9ad1fd257585ae4f5f0fd89`).

## Exact resource accounting per round

| design | data | ancillas | total qubits | data 2q | cat prep+merge | all 2q gates | data-only depth | all-2q depth | full operation depth | qubits x all-2q depth |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| two-ancilla PBB | 144 | 216 | 360 | 1008 | 144 | 1152 | 7 | 9 | 11 | 3240 |
| one-ancilla PBB | 144 | 144 | 288 | 1008 | 0 | 1008 | 8 | 8 | 10 | 2304 |
| Gross | 144 | 144 | 288 | 864 | 0 | 864 | 7 | 7 | 9 | 2016 |

## Verdict and limitations

**NEGATIVE_STRUCTURAL.**  This specified scheme cannot recover Gross's
seven-layer round: depth 7 is exactly infeasible and the achieved physical
construction has depth 9.  This is not a best-multi-ancilla-scheme claim, a
logical-error claim, or a fault-tolerance proof.  Flag-assisted and verified
cat variants remain open.

The 3-way BP+OSD benchmark is not reported here.  Its v1 feasibility failure
and frozen v2 protocol are recorded in
`results/partial_runs/exp024_benchmark_feasibility_amendment.json`; the
zero-byte canonical benchmark sentinel remains untouched.

## Reproduction

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
NUMEXPR_NUM_THREADS=1 PYTHONPATH=src \
.venv/bin/python experiments/exp024_cat_extraction.py --structural-only
```

Wall time: 70.180 s on
macOS-26.5.2-arm64-arm-64bit-Mach-O.  Depth-8 status is budget dependent;
no exhaustiveness claim is made for that search.
