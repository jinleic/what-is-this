# CORRECTION — prereg-provenance hash drift in frozen run
# 20260901T132102Z_f17bb4cb_117823c32d87 (agent OmegaPcPen, 2026-09-01;
# owner Main audit finding, verbatim item 1)

## The defect

`campaigns/20260901T132102Z_f17bb4cb_117823c32d87/PREREG_PROVENANCE.json`
records:

- `source_sha256` = `443fc30235168daac6fa2f37c1807a79465e12ac842c6b8bf1ee
  1f7f02492e31`
- `init_prereg_sha256` = `443fc302…92e31` (same)
- `byte_identical` = `true`

That was true AT COPY TIME (the copy matched the init-minted manifest's
`prereg_sha256`, both `443fc302…92e31`). However, the in-dir
`pre_statement.md` was then AMENDED (Amendment 1, commit `2bc7628`, before
any verdict compute; Amendment 2, commit `7305623`, display-doc only), so
the file NOW hashes `03fae37c60dcc6a4c2480b00d215beb7f2de0607511d8bee71f3
7ccff2e36ec4` (pinned by the run's `sha256s.txt`). The provenance record
therefore reads as if the copy still matches init, which it no longer does.

## The corrections (exact)

- run id: `20260901T132102Z_f17bb4cb_117823c32d87`
- init hash (manifest prereg_sha256, claim-time bytes):
  `443fc30235168daac6fa2f37c1807a79465e12ac842c6b8bf1ee1f7f02492e31`
- amended in-dir hash (current bytes):
  `03fae37c60dcc6a4c2480b00d215beb7f2de0607511d8bee71f37ccff2e36ec4`
- Amendment 1 commit: `2bc7628`; Amendment 2 commit: `7305623`
- ordering: the amendment PRECEDED all verdict compute (attempt 4 onward);
  no verdict-bearing number depended on the pre-amendment bytes
- authoritative bytes for `pre_statement.md`: the frozen run's
  `sha256s.txt` entry pins `03fae37c…36ec4` — the amended file is what
  the campaign ran under; the frozen run dir itself is NOT modified (its
  checksums ledger would break).

## General rule for successors (binding on this target)

If a prereg is amended after `campaign.py init`, the provenance record
MUST carry BOTH hashes: the init-time hash (what the manifest minted)
AND the final amended hash (what the run actually executed under), with
the amendment commits in between, and the run's `sha256s.txt` is
authoritative for the final bytes. A provenance file that shows only the
init hash of an amended prereg is defective in exactly the way recorded
here.
