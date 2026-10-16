# Signed inner/lower m15 consumed by preRp - 2026-10-16

`OriginalPreRpSignedInnerM15Numeric.evaluate(Z,chart,coordinate)` now consumes the direct signed cumulative core/cutoff/lower m15 at rminus in the compact and quiet assemblies through Rp. This closes the pending external API substitution from the preceding milestone. The full goal is **ACTIVE / INCOMPLETE**; tasks 14/15 remain IN PROGRESS.

## Source replacement and physical units

The adapter composes the existing `OriginalPreRpM15Numeric` and reads the checked `OriginalN1InnerM15Numeric` q=2 packet. It requires the same original family, exact selected sc and phase sc/4, actual Z, lower-collar chart and complete ordered six-term source partition. It removes both legacy `inner` and `lower` magnitude terms once and inserts:

1. `core_axial_coupling`;
2. `core_pressure_integral_and_finite_boundary`;
3. `first_lower_axial_width`;
4. `first_lower_regular_pressure_width`;
5. `first_lower_singular_pressure`;
6. `first_lower_singular_width_error`.

Each term's physical log atoms and ordinary C0/C1 coefficient/error rows are retained without another p=2 normalization. The active, compact finite-N, exact compact leading and weighted own-swirl terms remain unchanged. Quiet additionally retains its exact leading history. There are ten terms in compact and eleven in quiet. The prior owner and its cached packets are not mutated.

The cumulative inner contribution already includes the core, selected velocity cutoff and lower pressure history through rminus. It is not an increment to add to the old inner/lower terms. Nonzero leading pressure after own velocity support remains present. Input and implementation hash conflicts are rejected. Equivalent exact decimal Z strings select the same saved packet; new Z values delegate to the same direct inner source engine with its complete error bounds.

## Function references and numerical scope

This is a direct numerical source adapter. It does not register its results as archived FunctionRef nodes or infer source identity from interval overlap. The inherited source partition theorem is kept as a reference, separately from the computed source provenance. Its exact identity is not an independent numerical Rp cancellation measurement.

The inserted terms have signed approximations and full directed source errors. Active and compact finite-N pressure still use the previous enclosure terms; resolved signed global centers, useful relative cancellation accuracy and a global five-moment numerical join remain unproved. Compact-to-quiet equality checks establish unchanged history transport at that chart boundary, not the independent Rp-to-pulse seam.

## Focused evidence

The focused construction exited 0 in 66.722 seconds; saved-report inspection also exited 0. The report is `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_preRp_signed_inner_m15_numeric.json`, 2,767,007 bytes, SHA256 `a0a22152a347be33a866590a3e93fd9e1c788be32ae8f90e3efdfd0c1266c6a4`, pinning 1840 inputs and 81 implementations.

- 21 packets at Z=-.3,0,.3, compact x=1,1.25,1.4,2 and quiet s=0,.5,1.
- 42 exact replacement/absence conditions, 252 unchanged physical-unit/coefficient conditions and 21 unchanged outer-source comparisons.
- 420 partition conditions, 168 provenance/open-scope conditions and six compact/quiet/Rp history checks.
- 18 exact empty compact-increment zeros, six decimal alias checks, seven rejected malformed partitions and five rejected invalid queries.

Run `python scripts/run_original_preRp_signed_inner_m15.py --inspect-report` to inspect saved evidence. Rerun `--write-report` only when changed implementations require fresh construction. Existing read-only worker: `/root/n1_atlas_map`, **GPT-5.6 Luna / max**. No new agent was created.

Its read-only review confirmed the single replacement, unchanged physical units/outer terms and limited provenance claims. These six terms represent the n1 core/first/lower route; the adapter is not a complete signed finite-N modified-leading pressure solution.

## Next pressure construction: first_micro is not an unchanged field

The source scan resolves a material routing issue for the next interval:

- The leading first-collar and first-micro backgrounds share `rho=4 exp(hb*s)`, `R=epsilon*rho` and `dy=hb ds`. The finite-N phase origin shift does not shift this absolute radius.
- The actual finite-N selector is zero only on its admitted unchanged/flat left collars. In `current_original_micro_finite_N.py`, the flat collar ends at 3sc/4. The complete interval [sc/2,1/4] and endpoint 1/4 are classified `original_modified_right`.
- Consequently `OriginalFirstCollarAxial2` supplies the leading background only. Its source cannot be promoted to the actual modified first-micro pressure by declaring finite-N terms zero.
- The existing partial-mean pressure provider starts at first-micro coordinate 1/4. Extending the signed active integral therefore requires an explicit source route for the uncovered [sc/2,1/4] slice.

Required next actions, in order:

1. Extend the actual fixed-source U/Q/Q_Z/Q_y pressure route over that first-micro slice, retaining ordinary Z derivatives, all nonzero incoming histories and finite RQ_y/RQ_yZ endpoints. Use native s; do not materialize the impossible ratio s/(sc/4).
2. Keep the leading background and actual selected-N correction separate. Include `E_N=E exp(A/N)` and `V_N=V+B/N` with their original controls, phase, physical measures and derivative errors. Existing interval correction caps are bounds, not resolved signed values.
3. Integrate the signed S1 pressure with complete source errors, then cover the remaining active charts and actual compact finite-N pressure. Reuse the checked q=2 cumulative inlet exactly once. Preserve exact compact/quiet K-squared terms and weighted own swirl.
4. Preserve selected inverse/defect correlations to measure the independently transported Rp value/Z1 cancellation. Do not replace that measurement with a terminal source identity or interval overlap.
5. Complete all five histories, mixed joins, common global N, flatten/analytic preheat/exact heat and finite-energy tail, then admissible stress/flat remainder, genuine temporal n>=2 recursion/smooth sum, oscillatory correction and independent corrected Cartesian NS validation.

No heat, stress-cone, temporal recursion or corrected NS gate is closed by this adapter.
