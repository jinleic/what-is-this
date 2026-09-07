# N7 pre-statement — gate `n7-split-exclusion-theorem`

## 1. Question and why this is the shortest new frontier

N6 froze `INCONCLUSIVE-CHART`: the fixed N3 T5 anchor (term-zero split,
`w = c_W = (S0^T) e_3`, i.e. row 3 of the serialized S0 per
`T5_SPLIT_AMENDMENT.md`) is an exact 13-term zero of the CP parametrization
whose full exact
complex Jacobian has rank 161/192 (realified 322/384, kernel 86) — no square
chart exists there. A disclosed pre-prereg probe (2026-09-04, exact
elimination; preserved below as a control) found rank 161-162 across 33
split candidates and never 192. This preregistered gate does not rerun that
finite population. It registers an analytic theorem that, given one exactly
computed premise, decides **every** single-term split witness at once — any
term index `s` and any split vector `w` over `Q(i)` — and either excludes
the entire single-term split mechanism permanently or records that the
mechanism remains open, with machine-verified evidence either way.

## 2. Theorem (field-independent, registered before compute)

**Statement.** Let `K` be any field, and fix the 13-slot CP parametrization
`Phi` on coordinate blocks `A in K^{m x 13}` (`m=3`), `B in K^{8 x 13}`,
`C in K^{8 x 13}` with residual map `(A,B,C) -> sum_s a_s (x) b_s (x) c_s`.
Fix twelve nonzero terms `t_1..t_12` and let `x12` be their packing (slot
13 zero). For any `s in {1..12}` and any `w in K^8`, let `x13(s,w)` be the
point obtained by replacing term `t_s = (a,b,c)` — where `a != 0`, `b != 0`
— by the two children `(a, b, c+w)` and `(-a, b, w)`, all other terms
unchanged. Then, with `T12 = image dPhi_{x12}`:

```
image dPhi_{x13(s,w)}  subseteq  T12  +  (K^m (x) b (x) w)  +  (a (x) K^8 (x) w),
```

and consequently, for `w != 0` (so that `a (x) b (x) w` is a nonzero
vector of the intersection), because the two added spaces intersect in at
least `span{a (x) b (x) w}`:

```
rank dPhi_{x13(s,w)}  <=  rank dPhi_{x12}  +  (m + 8 - 1)  =  rank dPhi_{x12} + 10.
```

For `w = 0` both added spaces are zero and the sharper bound
`rank dPhi_{x13(s,0)} <= rank dPhi_{x12}` holds directly (proof below).

**Proof (basis-inclusion argument).** The pair of children contributes to
the differential image exactly the vectors
```
da1 (x) b (x) (c+w) + a (x) db1 (x) (c+w) + a (x) b (x) dc1
(-da2) (x) b (x) w + a (x) (-db2) (x) w + a (x) b (x) (-dc2)
```

over arbitrary variations (child 2 carries factor `-a`, so its `a`-slot
variation contributes `-da2 (x) b (x) w` and its `c`-slot variation
`a (x) b (x) (-dc2)`). Expanding each `(c+w)` slot by bilinearity and
grouping:

```
= (da1 - da2) (x) b (x) w        in K^m (x) b (x) w
+ a (x) (db1 - db2) (x) w        in a (x) K^8 (x) w
+ (da1) (x) b (x) c + a (x) (db1) (x) c + a (x) b (x) (dc1 - dc2)
                                 in image of the ORIGINAL term s, hence in T12.
```

The ELEVEN unsplitted terms contribute their original variation images,
each contained in `T12` (identical terms at `x12`). Only bilinearity,
distributivity, and the additive inverse are used: the argument is valid
over every field, every characteristic, with no division.

**Dimension count, by cases on `w`.** If `w = 0`, both added spaces are
zero, so `image dPhi_{x13(s,0)} subseteq T12` and
`rank dPhi_{x13(s,0)} <= rank dPhi_{x12}` — stronger than the general
bound, with no subtraction needed. If `w != 0` (with `a != 0`, `b != 0` by
hypothesis), the added space `K^m (x) b (x) w` has dimension exactly `m`,
the added space `a (x) K^8 (x) w` has dimension exactly `8`, and their
intersection contains the nonzero vector `a (x) b (x) w` (take `da`
proportional to `a` and `db` proportional to `b`), so their sum adds at
most `m + 8 - 1 = 10` dimensions over `T12`. Substituting `m = 3` gives
`+10`. **End proof.**

**Decisive corollary.** If `rank dPhi_{x12} <= 181` then EVERY single-term
split witness — every `s`, every `w in Q(i)^8` — has rank at most `191 <
192`, so no square chart exists at any of them, and the entire single-term
split mechanism (including N6's fixed anchor) is permanently unchartable.
If `rank dPhi_{x12} >= 182`, the theorem does not decide the mechanism and
the run records that honestly as nondecisive.

## 3. Exact audit (registered population: the theorem premise and checks)

The driver, importing nothing from any N6 module:

1. **Bindings.** Hash-check the frozen N3 witness
   `251e00577b43b400aa97d306bb923a98eb79ca7850ac67652d45288c2e570037`, N3
   instrument `d529cce2938af06645529f624c58a877a9fc3c45f6bc1efc8affafa16cca94ef`,
   the N6 verdict `051bc4f0be98ded4afebfe9b157f78a0491809414a8dd3b23cdd4eeb2a30967f`
   (lifecycle context), this preregistration, and the N6
   `T5_SPLIT_AMENDMENT.md`
   `88378ef0dfb88921d97fa22a5d3a71e036ecda38cca4c8906ba01b433b08f361`.
2. **TF and witness gates.** Rebuild TF from the Hamilton table (24 nonzero
   entries); verify 192/192 exact substitution of the twelve-term witness,
   its coefficientwise conjugate, and singleton corruption rejects
   `TF[0,0,0] += 1` for both.
3. **Premise.** Pack the twelve terms (slot 13 zero), form the exact
   192x247 Jacobian over `Q(i)`, and compute its exact rank by
   self-contained Q(i) Gaussian elimination: `rank12`.
4. **Calibration (control).** Build the N6 anchor `x13(0, c_W)` with
   `c_W = (S0^T) e_3` exactly as in `T5_SPLIT_AMENDMENT.md`; require its
   exact rank to equal the N6-audited `161`, and require
   `161 <= rank12 + 10` (theorem consistency on the recorded anchor). Any
   mismatch is INVALID-INSTRUMENT.
5. **Basis-inclusion machine check.** At the same anchor, verify the
   theorem's containment EXACTLY: every one of the 247 columns of
   `dPhi_{x13}` must lie in `T12 + (K^3 (x) b_0 (x) w) + (a_0 (x) K^8 (x) w)`,
   decided by exact rational elimination against a serialized basis of that
   sum space. A single failure is INVALID-INSTRUMENT (it would falsify the
   registered proof).
6. **Terminal serialization.** `rank12`, the calibration ranks, the
   containment result, and the theorem verdict
   (`SPLIT-MECHANISM-EXCLUDED` iff `rank12 <= 181`, else
   `THEOREM-NONDECISIVE`) are written atomically to `n7_results.json` with
   `no_rank_inference: true`.

## 4. Preserved preprobe (motivation/control, not the registered population)

The 33-candidate exact-rank preprobe (term `s in {0,1,2}` crossed with
eleven exact rational `w`, all ranks 161-162, never 192) is preserved
verbatim as `n7_preprobe_control.json` with its generating code reference.
It motivated this gate; the registered population is the theorem's universal
quantifier, not those 33 points.

## 5. Runtime, lifecycle, and scope

Python 3.14.3; `PYTHONDONTWRITEBYTECODE=1`, `sys.dont_write_bytecode=True`,
nice 10, BLAS threads pinned before import, `RLIMIT_CPU=1800` soft/hard, and
a 4 GiB RSS ceiling with orderly persisted abort — all enforced from process
start. Launch uses the owner-directed scheduler wrapper (runtime-only):
`/usr/sbin/taskpolicy -a nice -n 10 <venv python> n7_split_exclusion.py`.
One finite process; no replay process (no certification is claimed); no
randomness; no writes to any frozen run. Commit this preregistration
path-scoped, init the run, copy byte-identically, bind source hashes in
provenance, and launch only after an independent proof-and-instrument
review returns PASS.

Either terminal state is final for this gate. No outcome makes a real-rank
claim: the `{13,14}` frontier, the frozen lower-13 chain, and
`18 <= R_R(T_O) <= 25` are untouched; an exclusion verdict is an instrument
impossibility theorem about single-term split witnesses, not about the rank
of TF.
