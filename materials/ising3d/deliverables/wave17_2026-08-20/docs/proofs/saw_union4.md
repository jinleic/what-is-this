# Wave-14 height-schedule sieve union: `a_35 >= 1972465461070186257835`

Artifacts: `experiments/e141_saw_union4.py`,
`tests/test_saw_union4.py`, `results/bounds/saw_union4.json`.

Reproduce:

```
PYTHONPATH=src .venv/bin/python experiments/e141_saw_union4.py
PYTHONPATH=src .venv/bin/python tests/test_saw_union4.py
```

**Status convention.** `[THEOREM]`/`[LEMMA]` as proved; `[COMPUTATION]` the exact finite
certificate; `[UNRESOLVED]` the deliberately undecided. The external input
`c_36 = 2941370856334701726560670` is `[EXTERNAL]` with a hashed source.

## 1. Statement

**[THEOREM — the certified endpoint]** With `F = 3D-cell free energy` etc. as in
`problem_specification.md`, exact `[COMPUTATION]`:

```
a_35 >= 1972465461070186257835
M    =  c_36 - a_35  = 2939398390873631540302835
K_c  >=  atanh(M^(-1/36))           (= 0.2122159753270231627267517174278577806020,
                                      with the 40-place downward floor as shown)
```

The union structure, exactly as recorded in the artifact (`union_certificate.expression`):

```
a_35 >= L + C - (C cap L) + F2 + Q2 + F3 + Q3
       + sum_U c(U) ( H_U - C cap H_U - L cap H_U + (C cap L) cap H_U
                     - F2 cap H_U - Q2 cap H_U - F3 cap H_U - Q3 cap H_U )
```

with `U` ranging over the active term table (2,780 rows) and `c(U) in {1, -1}` from the sieve;
the base `L + C` family is the wave-13 one (`507841195880595165555` for the 41-old grid), and
`F2,Q2,F3,Q3` the two- and three-gadget families. `c(U)` and the intersection values are all
integers in the file: 5,940,409 closure candidates, 2,780 active rows, over 8,344 covered
chains, 2,167 intersecting pairs, with the dedicated overlap correction `5030204334415261659120`.

## 2. The bounded independent verifier

`tests/test_saw_union4.py` is a clean-room verifier that does NOT import the producer. It:

* recomputes `T(13) = 142,016,661` and the `B(11, q)` triangle by a tuple-coordinate visited-set
  DFS (the producer's is packed) and compares them to the stored profile row-by-row over the
  `l<=11, tail<=11` band — exactly the hard count surviving an independent code; the 404
  member *counts* are then materialized in that band *plus* the (11-band) product, which for
  `T(13)`-shared rows is reported as `StoredOnly` (the verifier's own status line
  "members rebuilt (l<=11, tail<=11)"). The exact 40-dp floor is re-derived by `mpmath` at 150
  dps and compared to the stored floor.
* It takes the 41-old-grid subset, recomputes its `H_U` in the band, and (when the band
  intersects) the overlap-corrected per-term net, term by term; rows outside the band are read
  from the stored profile but their sieve *identity* `c(U)*(net)` is re-derived from the stored
  values. The per-row sieve terms (1,423 with `c=+1`, 1,357 with `c=-1`) are checked against
  the independent recomputation, in the 2,167 intersecting-pair rows only where the band allows
  a fresh DP.

Scope: **only the `l<=11, tail<=11` band is independently re-derived**. The full 404-member
exact counts and the 2,780-term sum rest on the producer's own 16 internal checks (all
passing in the artifact), and the stored `M`+floor rest on `c_36` and the external hash.
This is narrower than a full sweep and stated as such; the `T(13)` value and the 40-dp floor
are the parts independent end-to-end.

## 3. The bounded negative (the reason the band is what it is)

`T(13)` (13-step SAW with all non-origin heights `> 0`) is `142,016,661`, an index set of
~10^8; the `B(12,*)` and `B(13,*)` slabs are `4^12..4^13`-weight-class, which the DFS
band explicitly excludes. The band is there so the verifier never silences a wrong `T` (which
it cannot silently omit), never claims full-member fidelity it does not have, and always reports
its `0 rebuilt`/`0 mismatches` line when the band is empty. That honesty is intentional: every
wave-11-13 relaunch of this front had the same boundary and the same refusal to overstate.

## 4. What is and is not established

* **[ESTABLISHED, [COMPUTATION]]** `T(13)`, the 40-dp floor, the expression form, the
  regression old-41, the integer values `L`, `C`, `F2/Q2/F3/Q3`, the per-row `c(U)` table and
  the `M`.
* **[UNRESOLVED]** The `a_35` value itself and the 404 full member counts are producer-certified,
  not lead-certified here; the exact `B(12,*)`/`B(13,*)`-row counts outside the band come from the wave-13 strict DP,
  not from this DFS; `T(13)` itself IS re-derived here (independent code). The region outside the
  band rests on the wave-13 strict DP. The `L`/`C`/gadget intersection table outside the band is likewise
  read, not re-derived.
