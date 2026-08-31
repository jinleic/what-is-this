# Gate A precision sweep (Main's corrected diagnostic, Delcap-style)

Working-precision response of the Gate-A margin chain:

| working bits | b1_rad | margin_rad | margin_low > 0 |
|---------------|--------|------------|----------------|
| 512 | 6.10e-121 | 9.97e-121 | YES |
| 256 | 4.85e-78  | 4.85e-78  | YES |
| 128 | 1.65e-39  | 1.65e-39  | YES |
| 96  | 7.08e-30  | 7.08e-30  | YES |

Monotone width growth with collapsing precision; margin low never crosses zero.
No float leak; interval path live end-to-end.  (At 256 bits the full-grid run
gateA_final.py reports margin width 2.35e-16 because that chain accumulates
31752 A-grid partials, each with ~4e-19 width; this sweep recomputes only
b1's two-cell closed form so its width is tighter; both are valid.)
