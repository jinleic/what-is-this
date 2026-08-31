"""gate_b_ladder_report.py — summarize the rank 24..20 search rows into
the joint (residual, factor-norm) table Main requires. Read-only."""
import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
rows_out = []
for fn in sorted(glob.glob(os.path.join(HERE, "gate_b_rank*_results.json"))):
    d = json.load(open(fn))
    rows_out.append({k: d[k] for k in
                     ("rank", "starts_attempted", "best_rel",
                      "best_max_factor_norm", "best_factor_norms",
                      "best_max_column_product", "best_ratio", "seed",
                      "elapsed_s")})
print(json.dumps(rows_out, indent=1))
