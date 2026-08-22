# n=18 cube batch 2: validation and discovery launches

This ledger records the authoritative existing CNFs. **No cube build was run in
this batch**: `cube_build.py` writes these same five paths, so rebuilding would
overwrite the checked artifacts.

Generator contract checked from `scratch/kobon/n18/cube_build.py`:
`N=18`, `TARGET=94`, and cases `q0`, `q1`, `q2`, `q3paired`, `q3triad`.
The DIMACS headers below match the corresponding case files and target.

## CNF fingerprints

| case | file | DIMACS header | bytes | SHA-256 | batch-2 launch state |
|---|---|---:|---:|---|---|
| q0 | `n18_cube_q0_t94.cnf` | `p cnf 1708348 29932035` | 1,063,540,279 | `69b0380b0e79aa76f8a6863d352e3a13955d4da4e12b41900f4a707403605bbe` | `launched=false` here; pre-existing solver active |
| q1 | `n18_cube_q1_t94.cnf` | `p cnf 1682318 29879983` | 1,062,343,488 | `1cac82769b80b761f986b2276adb71946d976532326d61e6a8d08a2927a959cf` | `launched=false` here; pre-existing solvers active |
| q2 | `n18_cube_q2_t94.cnf` | `p cnf 1656280 29827915` | 1,061,146,329 | `d88e051c848a0840d376994c3e4017058de70ee9941c7872f067ef9f9d32e07b` | `launched=true` below |
| q3paired | `n18_cube_q3paired_t94.cnf` | `p cnf 1630234 29775831` | 1,059,948,802 | `4b957d5e4204cdc2627bf2e3716bca73d359e5d51aa59ad22b1469348c2a03c8` | `launched=false` (held; theorem-M Q=3 branch) |
| q3triad | `n18_cube_q3triad_t94.cnf` | `p cnf 1630228 29775819` | 1,059,948,532 | `1c48a99ce8d0380f9e58f2f051958f4816b23906abf7938d369d2059dd115382` | `launched=false` (held; theorem-M Q=3 branch) |

The q2 log independently confirms the parsed header and 1,061,146,329-byte
input. All five files were fingerprinted before the q2 launch.

## Process-state check

Before launching q2, the required CPU-filtered solver count was **12**, below
the `<14` ceiling. The equivalent check used was:

```sh
ps -axo %cpu,args | awk 'BEGIN{c=0} /kissat|cadical/ && $0 !~ /awk/ && $0 !~ /grep/ && $0 !~ /tee/ && $1 >= 1 {c++} END{print c}'
```

The n18 entries observed after launch were:

- q0: Kissat PID `97712`, active, `n18_cube_q0_t94.cnf`.
- q1: Kissat PID `96954`, active, `n18_cube_q1_t94.cnf`; Cadical wrapper PID
  `9055` (`run_cadical.py q1`) also active.
- q2: wrapper PID `10443` and Kissat PID `10488`, active,
  `n18_cube_q2_t94.cnf`.

## q2 discovery launch

Stable launch name: `n18-q2-discovery` (hub persistent process).

Exact command:

```sh
/Users/jinleic/jinleic-workspace/scratch/kobon/solve_sat_capture.sh \
  /Users/jinleic/jinleic-workspace/scratch/kobon/n18/n18_cube_q2_t94.cnf \
  /Users/jinleic/jinleic-workspace/scratch/kobon/n18/n18_cube_q2_t94.kissat0.log \
  0
```

This is discovery mode only: stdout is captured in
`n18_cube_q2_t94.kissat0.log`; no DRAT path or proof output was supplied.
The log banner and parsed DIMACS header are present.

No second launch was made. q3paired and q3triad remain explicitly queued as
`launched=false` and must not be started in this batch.
