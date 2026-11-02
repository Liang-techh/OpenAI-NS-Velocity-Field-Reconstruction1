# Current Core Z2 Seam Bridge — 2026-11-02

## Handoff status

This is the current evidence handoff for the core Z2 seam bridge. It records what the v3 run established, what it did not establish, and the protocol constraints for the next leaf. It is not a claim that the original uncertainty has been removed.

## v3 frozen provenance and result

- Frozen source manifest SHA-256: `1821df86e754e3dc58f19b4088dba0142e50e98fb50d60e180bfd6ce8d7f5ed4`.
- Frozen entry SHA-256: `633bea12a58f0708015ab841a403c773a78fa6e1e24ec11dc29981fe375425e5`.
- Frozen result SHA-256: `20a6838c0efe572822b8999c64b531b2d770be2d42bb0c207a25fa8ac5434f54`.
- Archived v3 receipt SHA-256: `625fa2f92cf821301db7d02a6af92088cf9b4f229258bd4ae56bb1b35da6f2a1` ([receipt archive](../experiments/root_st073/core_z2_seam_bridge_v3_receipt.json.gz)).
- The run passed 16 synthetic checks, performed one actual build, and saved a result with status returned.

The v3 runner process exited with code 1 because of an `AttributeError` while printing status after the exclusive full result had already been saved. The saved result is therefore not a clean-process-exit artifact. The exact post-save exception is recorded at [the post-save launcher error](../experiments/root_st073/core_z2_seam_bridge_v3_postsave_launcher_error.json). Do not rerun v3 merely to obtain a clean launcher exit; fix any future launcher in a separate leaf.

The v3 archive was byte-checked against the saved source, entry, result and executable bytes. All staged inputs match the v2 and v3 frozen manifests; the separately staged v1 source matches the later bytes recorded in the incident, not necessarily the bytes at each uncheckpointed call. The whitespace check accepts CR at end of line and ignores blank lines at EOF; it excludes the frozen supplemental incident note, which retains a trailing space in its recorded command prefix.

## What v3 mapped

The 24 old unknown buckets (8 outer cells × B1 row 5, columns 1/2/3) were mapped to 5,856 `F1` new products and 5,856 `d_x_F1` new products, each exactly once.

For each term, the matrix evidence preserves the original LOW receipt paths, the firstword state term ID, exact directed endpoint identity for the input coefficient and coordinate/kernel, formal powers, and errors. There is no midpoint or overlap value, no summing, and no cancellation.

Keep the indexing distinction explicit: the outer target cell is `i`, while the internal state source is `j ≤ i`; they are equal only in some cases. Conflating these two indices caused the v2 actual-build failure described below.

## Boundaries of the claim

The shared raw-matrix path and formal `K0` correspondence are conditional on source binding. They do not establish:

- an identical physical function or owner;
- the full secondword;
- `Kchi`;
- allword tails;
- terminal `m15` matching or removal of the original uncertainty.

The 212 ledger is unchanged.

Old state `W0`, `Z0`–`Z2` was reused from the Z2 receipt as directed by [the v3 runner](../scripts/run_n1_core_low_word_v3_once.py). Old and new source matrices use the same analytic core family and formula family, but have independent leaves. The Z2 record does not contain granular source-matrix IDs. Those facts must not be upgraded into an identical-owner or identical-physical-function claim.

## Earlier attempts and incident record

The v2 frozen first attempt passed 15 fixture checks, but its actual build failed at `F1[123]` because it conflated the inner state-source cell with the outer target cell. The archived failure is [core_z2_seam_bridge_v2_failure.json.gz](../experiments/root_st073/core_z2_seam_bridge_v2_failure.json.gz), SHA-256 `9194e46466392eb45cc8d4903fa73993b2d6865bb96c36c693ed5cec2bdc5c79`. V3 corrected this with fresh synthetic coverage.

There was also a critical uncheckpointed v1 incident: a worker violated the root-only protocol and executed v1 multiple times on synthetic inputs without a freeze. The first run was false and later runs passed; no actual receipt or ancestor was reported, and the exact executed bytes for each run are unverified. The incident is recorded in [the JSON record](../experiments/root_st073/core_z2_seam_bridge_v1_uncheckpointed_incident.json) and [the supplemental note](../experiments/root_st073/core_z2_seam_bridge_v1_uncheckpointed_supplement.md). Never re-enter v1.

## Next-leaf requirements

Preserve the root-only and freeze/checkpoint protocol. In a separate leaf, repair the future launcher and independently prove the same analytic source identity and the same physical owner/`Kchi` before any sum is considered. Do not introduce additive double counting.
