# A strictly improved certified lower endpoint: K_c > 0.2138368062, by an exact finite-memory SAW automaton

Artifacts: `experiments/e227_finite_memory_saw.py` (producer),
`results/bounds/finite_memory_saw.json` (artifact),
`tests/test_finite_memory_saw.py` (clean-room verifier; imports no producer code).

Reproduce:

```sh
.venv/bin/python experiments/e227_finite_memory_saw.py --max-memory 12
.venv/bin/python tests/test_finite_memory_saw.py
```

**Status convention.** `[THEOREM]` proved here or in a cited repository note;
`[COMPUTATION]` exact finite calculation; `[EXTERNAL]` imported published
result with the exact location cited; `[UNRESOLVED]` deliberately undecided.

---

## 0. Result

**[THEOREM][COMPUTATION][EXTERNAL].** For the simple-cubic lattice,

```
mu  <  4747526/1000000  =  4.747526,
K_c >  atanh(1000000/4747526)
    >  0.213836806210869730430209429387937989133306728225224073121671,
```

so the certified critical interval improves from

```
[0.2122159753270231627267517174278577806020, 0.2527310098586630030260020266135701299926]
```

to

```
[0.2138368062108697304302094293879379891333, 0.2527310098586630030260020266135701299926]
```

a `4.00%` width reduction, and the first movement of either endpoint since the
wave-14 SAW sieve union.  The `mu` bound is marginally stronger than the
printed Poenitz-Tittmann value `4.7476` (their table rounds up to 4 decimals).

## 1. The three logical steps and their tags

1. **[THEOREM] (memory-k domination).** A *memory-k walk* revisits no vertex
   within `k` steps.  Every self-avoiding walk is a memory-k walk, so
   `c_n <= f^(k)(n)` for all `n`, and `mu <= rho_k := lim f^(k)(n)^(1/n)`.
2. **[THEOREM][EXTERNAL] (automaton).** Memory-k walks are counted exactly by
   a finite automaton whose states are the symmetry classes of the longest
   walk suffix still completable to a simple loop of length `<= k`
   (A. Poenitz and P. Tittmann, *Improved Upper Bounds for Self-Avoiding
   Walks in Z^d*, Electron. J. Combin. 7 (2000) #R21, DOI 10.37236/1499,
   Section 2; full text obtained and hashed in `sources/manifest.yaml`).
   `rho_k` is the Perron root of the transition matrix `A_k`.
   The kill test is exact in both directions: the retained suffix holds at
   most the last `k-1` vertices, and at a *first* violation the closing
   segment is a simple cycle, so its start is still retained.
3. **[THEOREM] (Collatz-Wielandt certificate).** If `w` is any strictly
   positive integer vector with `q * A_k w < p * w` componentwise, then
   `rho_k < p/q`.  This is one-sided and needs no eigenvalue computation:
   iterating gives `q^n (A_k^n w)_i < p^n w_i`, and column sums of `A_k^n`
   against `w >= 1` dominate `f^(k)(n)` up to the constant `max w / min w`.
   The stored `w` for `k = 12` has 41,424 entries; the minimum integer
   residual `p*w - q*A_k w` is `18,990,104 > 0`, replayed in exact integer
   arithmetic by producer and verifier on *independently rebuilt* transition
   tables.
4. **[THEOREM][EXTERNAL] (Ising transfer).** `chi(v) <= sum_n c_n v^n` for
   `0 < v < 1` (SAW domination of the two-point function), so `v_c >= 1/mu`
   and `K_c >= atanh(1/mu) > atanh(q/p)`.  This is the identical inference
   already used by `proofs/kc_bounds_improved.md` for the `c_36` route; only
   the `mu` upper bound is new.

The decimal `0.2138368...` is a directed floor of `atanh(q/p)` at 60 digits
(mpmath, 100 dps in the producer, re-derived at 120 and 200 dps in the
verifier with agreement below `1e-100`).

## 2. Float discipline

Floating point appears in exactly one place: proposing the positive vector
(power iteration on `A_k + I`; dead-end suffix states carry Perron weight
zero, so their components decay and are floored before integerization).  The
proposal is then rounded to integers and the strict inequality is *replayed
in exact integer arithmetic*; the float path can affect only whether a
certificate is found, never its validity.  The published-table controls are
consistency checks, not inputs.

## 3. Controls

* **Exact closed form.** The rebuilt `k = 4` automaton has 3 states and
  characteristic polynomial `x^3 - 4x^2 - 4x - 1`, exactly the reciprocal of
  the Poenitz-Tittmann denominator `1 - 4z - 4z^2 - z^3` at `d = 3`
  (their Section 2 displays the `d`-generic form).
* **Published table.** Float brackets at `k = 4, 6, 8, 10, 12` are
  `4.8645365, 4.8074109, 4.7779699, 4.7598376, 4.7475251`; each ceils to the
  printed Table 2 values `4.8646, 4.8075, 4.7780, 4.7599, 4.7476` ("the
  values shown are true upper bounds").
* **Direct enumeration.** For `k = 4, 6` and `n <= 7` the automaton counts
  equal a brute-force enumeration of memory-k walks; for `k = 12`, walks of
  length `<= 7` cannot contain a loop of length `<= 12` without being
  non-self-avoiding, and the counts equal the direct SAW counts
  `6, 30, 150, 726, 3534, 16926, 81390`.
* **Structure digests.** Producer and verifier build the automata with
  different closability searches (BFS with L1 pruning vs depth-first) and
  independently coded canonicalization; SHA-256 digests of the state and
  transition tables agree at every `k`.

## 4. Scope and honesty

* The bound is a *lower* endpoint improvement only; the upper endpoint
  `I_3/2` is untouched.
* `K_c = 0.221654626` was not used to select, tune, or stop anything; the
  memory cutoff `k = 12` was fixed by the resource budget (41,424 states,
  44.7 s CPU, 425 MB peak RSS) before the bound was known.
* `k = 14` is the next rung (Poenitz-Tittmann print `4.7387`, which would
  give `K_c > 0.2142`); projected ~15x states (~620k) and a growing
  closability cache. `[UNRESOLVED]` here; it is a resource frontier, not a
  mathematical one.  The published `mu < 4.7114` (MacDonald et al. 2000)
  remains excluded: its full text is still unobtained, so nothing here
  changes that audit verdict.
* The automaton-correctness step is tagged `[EXTERNAL]` on Poenitz-Tittmann
  Section 2 with the full text on hand; the in-repo controls above
  independently pin the construction at every memory used.

## 5. Claim ledger

| statement | status |
|---|---|
| `c_n <= f^(k)(n)` for all `n, k` | `[THEOREM]` (inclusion of walk classes) |
| automaton counts memory-k walks exactly | `[THEOREM][EXTERNAL]` P-T Sec. 2 + controls above |
| `q A w < p w` with positive integer `w` implies `rho < p/q` | `[THEOREM]` (Collatz-Wielandt, one-sided) |
| `mu < 4.747526` | `[THEOREM][COMPUTATION]` integer replay, both builds |
| `K_c > atanh(1000000/4747526)` | `[THEOREM][EXTERNAL]` SAW domination, as in `kc_bounds_improved.md` |
| 60-digit decimal floor | `[COMPUTATION]` directed, dual-precision agreement |
| `k = 14` extension | `[UNRESOLVED]` resource frontier |
