# Compact finite-memory SAW crosswalk, k <= 12

Artifacts: `experiments/e229_finite_memory_saw14.py`,
`results/bounds/finite_memory_saw14.json`, and the independent verifier
`tests/test_finite_memory_saw14.py`.

## 1. Outcome and exact scope

**[THEOREM — compact crosswalk].** For each even memory
`k=4,6,8,10,12`, the deterministic packed first-use implementation constructs
the same ordered state and transition tables as the independently verified
`e227_finite_memory_saw.py` automaton. The match is content-equivalent and is
active state/transition digest, not merely state-count agreement.

**[COMPUTATION — resources].** The complete `k <= 12` in-gen bounding crosswalk
uses 5.180090 process CPU seconds and 179,814,400 peak RSS bytes in the stored
producer run, far below the 2 GiB gate. Per-memory counts are exactly

\[
\begin{array}{r|r|r}
k & \#\text{states} & \#\text{transitions} \\
\hline
4 & 3 & 7 \\
6 & 20 & 69 \\
8 & 205 & 805 \\
10 & 2722 & 11074 \\
12 & 41424 & 169975 .
\end{array}
\]

**[UNRESOLVED — the actual memory-14 frontier].** No memory-14 automaton,
spectral radius, Collatz vector, or critical endpoint is produced here. The
producer's `k=14` row is explicitly `not_launched` and its endpoint claim is
null. This note therefore does not change the certified critical interval.

## 2. Packed first-use normal form

Represent a translated suffix path by its step word, where each step is one of
six signed axes. The compact word stores a 4-bit length followed by 3 bits per
step. A step

\[
\pm e_a
 \]

is encoded as `2a` when its sign is positive and `2a+1` when negative.

**[LEMMA — uniqueness of first-use normal form].** Every orbit of finite
signed-axis words under coordinate permutations multiplied by independent
coordinate sign flips has exactly one first-use representative: at the first
occurrence of each source axis it is renamed to the first still-unused output
axis, and its sign is chosen so that the first output step in that axis is
positive.

**Proof.** Suppose two signed permutations give two normal representatives.
Before the first nonvanishing component appears on any axis, coordinatewise
comparisons are unaffected. At that occurrence, the normal-form rule forces the
same earliest free output axis and the same positive sign. Repeating the same
statement at each next first-used axis fixes the entire output word. `[]`

The verifier is clean-room: it rebuilds words with a distinct recursive
canonicalization and an explicit BFS obstacle-contact closability test, then
checks the same ordered e227 state and transition digests for every
`k=4,6,8,10,12`. No memory-14 launch or spectral claim enters either the
producer or the verifier.

## 3. Route boundary

The compact representation removes the two immediate blocks listed before the
wave: tuple-coordinate states and the unbounded closability cache. It is a
mechanical foundation for a later `k=14` run, not that run's certificate. Any
subsequent endpoint change must still construct `A_14`, verify the mass, prove
`componentwise rho(A_14)w<w` by integer Collatz rotation, and route the result
through `atanh` with directed rounding.
