"""Render conservative per-symbol endpoints from archived 400-bit Arb balls."""
import csv
import hashlib
import json
import math
from pathlib import Path

from flint import arb, ctx, fmpq

CAMPAIGN = Path(
    "/Users/jinleic/jinleic-workspace/cs/delcap/campaigns/"
    "2026-08-31T09:13:25Z_7dc5babe-5e1b-44fa-a9e5-03d2f10183b6_q3-invariance-correction"
)
ROWS = CAMPAIGN / "corrected_q3_rows.jsonl"
ROW_HASHES = CAMPAIGN / "corrected_q3_rows.line_checksums.json"
OUT = CAMPAIGN / "conservative_intervals.jsonl"
OUT_HASHES = CAMPAIGN / "conservative_intervals.line_checksums.json"
CSV_PATH = CAMPAIGN / "TABLE_conservative_intervals.csv"
REPORT = CAMPAIGN / "outward_rounding_report.json"
ctx.prec = 400


def key(row):
    return f"{row['q']}:{row['n']}:{row['d']}"


def exact_float_rational(value: float) -> fmpq:
    numerator, denominator = value.as_integer_ratio()
    return fmpq(numerator, denominator)


def outward_lower(endpoint):
    value = math.nextafter(float(endpoint), -math.inf)
    rational = exact_float_rational(value)
    assert arb(rational).upper() <= endpoint
    return value, rational


def outward_upper(endpoint):
    value = math.nextafter(float(endpoint), math.inf)
    rational = exact_float_rational(value)
    assert arb(rational).lower() >= endpoint
    return value, rational


source_ledger = json.loads(ROW_HASHES.read_text())
assert source_ledger["algorithm"] == "sha256"
source_lines = [raw for raw in ROWS.read_bytes().splitlines(keepends=True) if raw.strip()]
source_rows = []
for raw in source_lines:
    row = json.loads(raw)
    assert hashlib.sha256(raw).hexdigest() == source_ledger["rows"][key(row)]
    source_rows.append(row)
assert len(source_rows) == 20

rendered = []
for row in source_rows:
    n = row["n"]
    primal = arb(row["primal_orbit_ball"]) / n
    dual = arb(row["dual_full_alphabet_ball"]) / n
    assert primal.upper() <= dual.upper()
    lo, lo_q = outward_lower(primal.lower())
    hi, hi_q = outward_upper(dual.upper())
    width_endpoint = dual.upper() - primal.lower()
    width, width_q = outward_upper(width_endpoint)
    assert lo <= hi
    rendered.append(
        {
            "q": row["q"],
            "n": n,
            "d": row["d"],
            "primal_per_symbol_arb_ball": str(primal),
            "dual_per_symbol_arb_ball": str(dual),
            "certified_lower_outward_decimal": repr(lo),
            "certified_upper_outward_decimal": repr(hi),
            "certified_width_upper_outward_decimal": repr(width),
            "certified_lower_outward_binary_rational": str(lo_q),
            "certified_upper_outward_binary_rational": str(hi_q),
            "certified_width_upper_outward_binary_rational": str(width_q),
            "lower_containment_assertion": "outward_binary_rational <= Arb primal lower endpoint",
            "upper_containment_assertion": "outward_binary_rational >= Arb dual upper endpoint",
            "width_containment_assertion": "outward_binary_rational >= dual.upper-primal.lower",
            "original_row_float_fields": {
                "cert_lo_per_symbol": row["cert_lo_per_symbol"],
                "cert_hi_per_symbol": row["cert_hi_per_symbol"],
                "cert_width_per_symbol": row["cert_width_per_symbol"],
                "status": "REPORT_ONLY_BINARY64_NOT_CERTIFICATE_ENDPOINTS",
            },
            "lbplus_report_float": row["lbplus"],
            "ub_report_float": row["ub"],
            "lower_gain_report_float": row["lower_gain_over_lbplus"],
            "upper_gain_report_float": row["upper_gain_below_ub"],
            "verdict": row["verdict"],
            "legacy_mixed_output_orbits": row["legacy_output_reference_audit"]["mixed_output_orbits"],
            "dual_from": row["dual_from"],
            "endpoint_evidence_label": "MACHINE-VERIFIED",
            "published_margin_verdict_evidence_label": "MACHINE-VERIFIED",
            "report_float_evidence_label": "COMPUTATIONAL-EVIDENCE",
        }
    )

ledger = {"algorithm": "sha256", "rows": {}}
with OUT.open("wb") as handle:
    for row in rendered:
        raw = (json.dumps(row, sort_keys=True) + "\n").encode()
        handle.write(raw)
        ledger["rows"][key(row)] = hashlib.sha256(raw).hexdigest()
    handle.flush()
OUT_HASHES.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")

fields = [
    "n",
    "d",
    "certified_lower_outward_decimal",
    "certified_upper_outward_decimal",
    "certified_width_upper_outward_decimal",
    "lbplus_report_float",
    "ub_report_float",
    "lower_gain_report_float",
    "upper_gain_report_float",
    "legacy_mixed_output_orbits",
    "dual_from",
    "verdict",
]
with CSV_PATH.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    for row in sorted(rendered, key=lambda item: (item["n"], item["d"])):
        writer.writerow({field: row[field] for field in fields})

report = {
    "status": "PASS",
    "rows": len(rendered),
    "precision_bits": 400,
    "source_rows_sha256": hashlib.sha256(ROWS.read_bytes()).hexdigest(),
    "source_line_checksum_ledger_sha256": hashlib.sha256(ROW_HASHES.read_bytes()).hexdigest(),
    "output_rows_sha256": hashlib.sha256(OUT.read_bytes()).hexdigest(),
    "output_line_checksum_ledger_sha256": hashlib.sha256(OUT_HASHES.read_bytes()).hexdigest(),
    "csv_sha256": hashlib.sha256(CSV_PATH.read_bytes()).hexdigest(),
    "endpoint_rule": "parse archived block Arb balls at 400 bits; divide by n; nextafter one binary64 ULP outward; assert exact binary-rational containment",
    "original_float_fields": "REPORT_ONLY_BINARY64_NOT_CERTIFICATE_ENDPOINTS",
    "endpoint_evidence_label": "MACHINE-VERIFIED",
}
REPORT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
print(json.dumps(report, sort_keys=True))
