# Nonstationary ladder certificates, failed transfers, and a finite L=9 bound

## 1. Scope and notation

Let \(\Lambda_L\) be the open two-leg ladder with \(L\) rungs, and let

\[
 A_L=\sum_{v\in\Lambda_L}X_v,
 \qquad B_L=\sum_{(u,v)\in E(\Lambda_L)}Z_uZ_v.
\]

A literal word in `A,B` is evaluated as a generator-left-nested commutator, with each
commutator divided by two in the \(Q_{(a\mid b)}=X^aZ^b\) Pauli convention.  Its length
is its bracket depth.  Write \(D_k(\Lambda_L)\) for the dimension of the span of words
of depth at most \(k\).

This note separates exact finite certificates from an all-length result.  In particular,
no result below asserts that the target \(D_{2L}(\Lambda_L)\geq2^L\) holds for every
\(L\).

## 2. Finite nine-rung certificate

**[COMPUTATION]** On \(\Lambda_3\), exact rational echelon reduction of every literal
word of depth at most six gives rank \(16\).  In discovery order, one exact basis is

\[
\begin{split}
& A,B,AB,AAB,BAB,ABAB,BBAB,AABAB,ABBAB,BABAB,BBBAB,\\
& AABBAB,ABABAB,ABBBAB,BABBAB,BBBBAB.
\end{split}
\]

For each of these roots, append six outer moves independently chosen from
\(\{AB,BB\}\), and evaluate on \(\Lambda_9\).  A depth-first enumeration (outer `AB`
before outer `BB`) examines 802 of the 1,024 resulting words and greedily selects 512
of them.  Every selected word has depth at most \(6+2\cdot6=18\).

Let \(V\) be the integer Pauli-coordinate matrix of these 512 selected operators.
The artifact uses the explicit integer map \(H\) from Pauli coordinates to 512
coordinates defined by

\[
 h(c)=\operatorname{mix64}(c),\qquad
 (Hq)_{h(c)\bmod512}\mathrel{+}=(-1)^{1-\operatorname{bit}_9(h(c))}q_c,
\]

where `mix64` is the recorded xor--shift--multiply map with fixed 64-bit constants.
The exact selected word list and the SHA-256 digest of the resulting integer
\(512\times512\) matrix \(HV\) are in
`results/algebra_growth/ladder_nonstat.json`.

The independently reconstructed ranks are

\[
 \operatorname{rank}_{2147483647}(HV)=512,
 \qquad
 \operatorname{rank}_{2147483629}(HV)=512.
\]

Therefore \(HV\) has a nonzero minor modulo either prime.  Since \(H\) is an integer
linear map, \(\operatorname{rank}_{\mathbb Q}V\geq
\operatorname{rank}_{\mathbb F_p}(HV)\).  Thus the finite certificate proves

\[
 \boxed{D_{18}(\Lambda_9)\geq512=2^9.}
\]

The two-prime calculation supplies independent lower-bound certificates; it is not
used to infer an exact rational rank.

## 3. Smaller growing-seed family

**[COMPUTATION]** Exact rational echelon at \(L=2\), depth at most four, selects the
six roots

\[
 A,B,AB,AAB,BAB,ABAB.
\]

Appending \(L-2\) outer moves from \(\{AB,BB\}\) gives the following modular rank
table.  The two recorded primes agree in every entry.

| \(L\) | number of words | maximum depth | rank modulo both primes | \(2^L\) |
|---:|---:|---:|---:|---:|
| 2 | 6 | 4 | 6 | 4 |
| 3 | 12 | 6 | 10 | 8 |
| 4 | 24 | 8 | 20 | 16 |
| 5 | 48 | 10 | 40 | 32 |
| 6 | 96 | 12 | 77 | 64 |
| 7 | 192 | 14 | 148 | 128 |
| 8 | 384 | 16 | 272 | 256 |

This is a finite nonstationary-seed certificate through \(L=8\), not an induction:
the same family was not shown to have an injective transition beyond the computed
range.

## 4. Exact failure of the first patched transfer at nine rungs

Let \(P_8\) be the known 248 retained `AB/BB` rows plus eight `BAB`-tail patches at
\(L=8\).  Its 256 operators have rank 256 modulo both recorded primes.

**[COMPUTATION]** The naive child family

\[
 \{ABw,BBw:w\in P_8\}
\]

contains 512 depth-at-most-18 words at \(L=9\), but has modular rank only
\(478\) modulo both primes.  Adding every unused seven-move `AB/BB` descendant of
`BAB` adds 112 words, for 624 total, and produces modular rank \(499\) modulo both
primes.  These ranks are lower-bound measurements only: they do **not** give a rational
upper bound, and they do not rule out another nonstationary patch or transfer.

## 5. A rigorous obstruction to the proposed restriction intertwiner

**[LEMMA]** For every \(L\geq2\), deleting the final rung is not an intertwiner for the
normalized \(B\)-adjoint.

Let \(E_{L+1}\) discard every Pauli term acting nontrivially on the added rung.  Let
\(v\) be the top site of the old last rung and \(u\) the top site of the new rung, and
put \(P=X_vZ_u\).  Then \(E_{L+1}(P)=0\).  The new horizontal bond \(Z_vZ_u\) has the
only contribution to \([B_{L+1},P]/2\) that becomes identity on the new rung, namely

\[
 E_{L+1}\!\left(\frac{[B_{L+1},P]}2\right)=-X_vZ_v\ne0.
\]

All other anticommuting bond contributions still carry a nonidentity Pauli on the new
rung and are erased.  Hence

\[
 E_{L+1}\,\operatorname{ad}_{B_{L+1}}(P)
 \ne \operatorname{ad}_{B_L}\,E_{L+1}(P)=0.
\]

For an explicit literal-word check, `ABABAB` evaluated at \(L=3\), then erased to
\(L=2\), differs coefficientwise from its direct \(L=2\) evaluation; the exact
coefficient-vector digest and all differing terms are recorded in the JSON artifact.
Consequently, this natural restriction map cannot justify the requested
block-triangular suffix induction by an adjoint-intertwining argument.  This lemma does
not exclude a different injective transition.

## 6. Longer-memory period-two candidates

**[FALSIFIED]** Six literal alternating one-move/three-move constructions were tested.
The `one_then_three` schedule starts from \(A,BA,ABA,BABA\), while
`three_then_one` starts from \(A,B,AB,AAB\).  Each then alternates a one-letter
move from \(\{A,B\}\) with a move from the displayed three-letter pair, in the
order named by its schedule.  At \(L=3\), exact rational ranks are strictly below
the required eight in all six cases, so each literal period-two family fails before any
all-length transfer question arises.

| schedule | three-letter pair | exact \(L=3\) rank over \(\mathbb Q\) | modular ranks for \(L=2,\ldots,8\), both primes |
|---|---|---:|---|
| one then three | `ABB`, `BAB` | 6 | 4, 6, 12, 24, 48, 96, 192 |
| one then three | `ABB`, `BBB` | 6 | 4, 6, 12, 24, 48, 96, 188 |
| one then three | `BAB`, `BBB` | 6 | 4, 6, 12, 24, 48, 96, 184 |
| three then one | `ABB`, `BAB` | 5 | 4, 5, 10, 20, 40, 80, 160 |
| three then one | `ABB`, `BBB` | 6 | 4, 6, 12, 24, 48, 96, 190 |
| three then one | `BAB`, `BBB` | 6 | 4, 6, 12, 24, 48, 96, 189 |

The full word counts, depths, both-prime ranks, and timings are machine-readable in the
artifact.

## 7. Conclusion

**[UNRESOLVED]** No integer unit-determinant transition matrix, explicit suffix map, or
all-\(L\) block-triangular minor was found.  The strongest new positive result is the
clean-room-verifiable finite certificate
\(D_{18}(\Lambda_9)\geq512\).  The strongest transfer information is negative but
narrow: the direct delete-a-rung intertwiner is impossible, the displayed patched
P8-to-P9 child construction has a modular deficit, and the six specified period-two
three-letter schedules already fail exactly at \(L=3\).

## 8. Reproduction

From the repository root:

```sh
timeout 7200 .venv/bin/python experiments/e91_ladder_nonstat.py
timeout 10800 .venv/bin/python tests/test_ladder_nonstat.py
```

The second command is a clean-room verifier: it does not import the producer and
rebuilds all Pauli vectors, characteristic-zero seed ranks, modular tables, the L9
integer hash minor, the patched transition ranks, and the restriction witness.
