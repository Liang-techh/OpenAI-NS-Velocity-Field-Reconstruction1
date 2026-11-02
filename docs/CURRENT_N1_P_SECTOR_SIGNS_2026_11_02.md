# N1 P-sector signs and scope — 2026-11-02

Goal **ACTIVE / INCOMPLETE**. This handoff records a finite physical
`m15` P-compensation candidate; it is not five-moment closure, a complete
Kchi field, or a Navier–Stokes validation.

The fresh v3 consumer ran once. Its source manifest SHA-256 is
`76b2360eeac14972a4ae8eee956ac378332c862295ac0ed3b5fd86582370ff20`, its
exclusive entry SHA-256 is
`32f2d52fa31fdd260b39f6815816191d3cb6e2a53fff20fa51fa3f40bfc9ef12`, and its
result SHA-256 is
`bec7fb6e8124c47491cc4b187b6340bde28048b234ab995fda89c729b2c33176`. The
complete gzip receipt is retained at
`experiments/root_st073/n1_p_sector_signs_v3_receipt.json.gz` with
raw SHA-256
`fa210083418b9b7e72e5414dde6dcbf32144d7e4011c539db61d2b8ead8c2dc6`.

## Returned sector inventory

Each of the four finite branches has 36 Z0 sectors. The Z1 counts are:

| branch | Z0 | Z1 |
| --- | ---: | ---: |
| `G0_P_W0` | 36 | 144 |
| `G0_P_DeltaW0` | 36 | 180 |
| `D_P_W0` | 36 | 180 |
| `D_P_DeltaW0` | 36 | 216 |

There are 864 sectors in total. Every returned interval contains zero; none
has a strict positive or strict negative sign, and the error-only count is
zero. The records retain binary interval endpoints, original source IDs and
paths, both ordinary-Z1 derivative classes, and the already-applied physical
weights. No coefficient addition was performed, and no old term was added to
the 212-term ledger or its six signed inner/lower contributions. These are
physical `m15` P-compensation sectors, not `P1`.

The independent synthetic fixture passed 14/14 checks. Existing error records
were preserved without deletion. The common source labels and atom endpoints
remain provenance compatibility only; they do not prove one common function,
gauge owner, or a normalized-equation identity.

## Failed predecessors

The v1 attempt failed in its synthetic fixture before entering real data. Its
failure archive is retained as
`experiments/root_st073/n1_p_compensation_cell_sum_v1_failure.json.gz`, raw
SHA-256
`012e2d410a74ee16eb6eaf6859bac00884a29514c55a9330006ce09d54e421a1`.

The v2 attempt failed while saving an `ivmpf` value that was not JSON
serializable. The real leaf may have executed, but no result is accepted. Do
not repeat any v2 coefficient sum. Its failure archive is
`experiments/root_st073/n1_p_compensation_cell_sum_v2_failure.json.gz`, raw
SHA-256
`29e6215eb2137317585592acd7aa638cc84f212a1d7654331a247fe6bf172efa`.

## Next executable boundary

The next bounded step is a proof of same-source signed correlation and gauge
ownership for the retained endpoint sectors, followed by an E/P composition
that demonstrates non-duplicate counting. P-tail ownership, an independent
pulse match, and terminal cancellation remain unresolved. No error
cancellation may be claimed from interval overlap or from the sector signs.

