# Current Core First-Word Z2 Profile — 2026-11-02

## Handoff scope

This round starts from the saved Z2 v2 receipt and its per-source-unit
`feedback_Z2` rows 4 / 5 / 6. The six source-cell terms are:

`C0_ZZ W0 + 2 C0_Z W0_Z + C0 W0_ZZ + B1_ZZ W0_Z + 2 B1_Z W0_ZZ + B1 W0_ZZZ`.

From those terms, the run constructs a K0 conditional directed-interval
profile for eight W1_Z2 target cells and six fields. It uses `t = x / 2`;
`j < i` uses complete cells, while `j = i` uses the moving-upper-bound
triangular Green kernel. Axis continuation is exact. Physical `d x` and the
derivative scales are retained. Formal powers, the input receipt path, source
IDs, and inherited uncertainty labels are retained; the upstream receipt does
not provide a separate per-term error decomposition, and matrix-sector
subdivision is still opaque.

## Counts and interpretation

Across all eight target cells, the term counts are:

| Field | Directed-interval terms |
| --- | ---: |
| `F1` | 1608 |
| `Uz1` | 2124 |
| `K1` | 2124 |
| `P1` | 108 |
| `d_x_F1` | 1608 |
| `d_x_Uz1` | 2124 |

The Z2 `W0_Z3` quantity in the earlier record was a conservative per-source
bound, not a physical function value. Its error has not been removed.

A separate calculation covers the previous Z1 missing `B1[0] * W1_Z2`
contribution in row 5, columns 1 / 2 / 3. This document constructs only that
branch: it has not been merged with the old item, and same-source
correspondence has not been certified. For each of the `F1` and `d_x_F1`
fields there are 5856 directed-interval terms. These include all causal
internal sources from the target cells, but only the K0-selected finite words;
the old record neither summed nor cancelled them.

Saved family / log lineage does not establish a matrix-sector exact identity
or a common physical owner. `Kchi`, gauge, all-word tails, and old errors are
not handled. Therefore this profile must not be described as closed, and it
does not cancel the previous 24 unknowns.

Frozen v2 validation status: 35 synthetic checks passed; the actual profile
was called once and the seam was called once; term counts are unchanged; no
ancestor was invoked. Every item carries the classification label
`inherited:feedback_Z2_directed_enclosure`. The seam additionally carries
`inherited:raw_B1_directed_enclosure` and `uncertainty = true`. These are
classification labels, not newly proved error bounds. No missing-branch claim
is made here: only the branch described above was constructed, without merging
it into the old item or certifying same-source identity.

The v1 receipt remains only as a QA-rejected metadata archive and is not
re-entered. This status does not assert that the v1 numerical values were
wrong; it records that its metadata was unreliable because it was generated
from input with no actual error labels (`error_classes = []`) and reported
`seam uncertainty = false`.

## Receipt and hashes

- Input Z2 SHA: `c53ff32ad3ce5288279953e84bc7417b71477c82b00d94466d0d2880adbc7ac9`
- Input low SHA: `5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6`
- Frozen v2 [receipt](../experiments/root_st073/core_firstword_z2_profile_v2_receipt.json.gz)
- Frozen v2 receipt SHA: `bcbd81914ba8294a8da13790640d63ab12907fd41624fceff10f304a5ca12115`
- Frozen v2 source SHA: `1e07d0c7a72b6c12ee1cb0a9b06fd59a01dc5a7e8c7a11eae6deb84fa5049817`
- Frozen v2 source entry SHA: `66acf8217d4bcf9e605116b0404b18d8f2293765fb3dced2ae9ffcce83754cbd`
- Frozen v2 source result SHA: `001f7ecfee729cd339d5456528493e6c4dc82b5cd3c15056ce0d43bb1f66b7c2`
- Archived v1 receipt (QA-rejected metadata only): `experiments/root_st073/core_firstword_z2_profile_v1_receipt.json.gz`
- Archived v1 receipt SHA: `3ff8d17e8a5a930f1d8a0c1139be85bb61962137294dd412b3a9087cdb206642`

Both archives were byte-checked against their saved source, entry, result and
frozen executable files. Staged executable hashes match both manifests. The
frozen files retain trailing blank lines; the staged whitespace check disables
only blank-at-EOF checking.

## Next work and guardrail

Next, establish same-source provenance and correspondence to the old missing
unknowns, audit all errors, identify `Kchi` and the common physical owner,
avoid duplicate accounting, and then connect the result to core / collar.
The 212 ledger is unchanged. Do not execute ancestor computation.
