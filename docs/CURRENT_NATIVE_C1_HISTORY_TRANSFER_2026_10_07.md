# Original inlet C1 memory and adjacent signed Duhamel transfer

> Successor: [CURRENT_NATIVE_TRUE_CHART_C1_TRANSFER_2026_10_07.md](CURRENT_NATIVE_TRUE_CHART_C1_TRANSFER_2026_10_07.md) ([45bcf343](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/45bcf343eb03513d137ab150db8170034580dc2d)) now supplies all17 true chart lengths, the factored C1 kernel adapter on five declared cells, and actual inlet-to3sc/4 whole-Z history transport. The active bridge and full inlet-to-O2/Rc path remain open.

Checked implementation: [df383c34](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/df383c34328c8ece9c81b2efd33daa6111cafb10). **The original five incoming histories and their Z rows, separate P0/P0_Z, common directed coordinate rebasing and adjacent C1 transfer now execute.** The serial operator covers the actual O2 radial coordinate[.12,.15] in three adjacent cells at Z=.5 and whole Z=[.49,.51]. This true width.03 is1,500,000 times the preceding microscopic cell width; it is expanded source coverage, not a percentage of the full reconstruction.

## Original inlet and pressure datum

The exact original inlet is r_minus=Ra*exp(hb*s_c/2), with phase coordinate s_c/2 strictly inside the first bridge. The existing NativeGenericSourcePackets.left_inlet path is reused. The actual radius binder uses selected_sc_multiple=1/2 and gives exact log(R/r_minus)=0 and fractional phase0 without adding a microscopic offset to a huge absolute radius log.

From the same unchanged packet, the implementation extracts packet.histories[j][0] and ordinary_axial_coefficient(row,1). The five normalized units are

```text
m = Mz/(R*Pstar)
h = Mtheta/(sqrt(2)*R^(3/2)*Pstar)
k = Mtheta_z/(sqrt(2)*R^(3/2)*Pstar²)
e = Mztheta/(R*Pstar²)
p = Mp/Pstar²
```

No additional width, radial or Pstar conversion is applied. The original incoming histories are retained even where nonzero. P0/Pstar² and P0_Z/Pstar² remain separate from p/p_Z. Original absolute-pressure function covers are P0+p and P0_Z+p_Z.

The initial *correction* histories at the actual inlet are exactly zero by the already checked strict flat-loop collar and Section11 initial condition. This is not a reset of the original background memory. The same collar proves q=A=B=0 without evaluating an unresolved active inverse. Both whole Z=[-1,1] and Z=[.49,.51] inlet queries execute. Across them,28 C0/Z original-history/datum/absolute-pressure rows are replayed.

## Common arithmetic coordinates

The old four log bases are width, Pstar², F0² and swirl-unit², with a separate radius log. Width differs between bridge/switch; swirl-unit and radius covers vary between queries; restore/reshape also use adapted amplitude bases. Matching exponent tuples alone is therefore insufficient to add cells from different packets.

CommonSourceCoordinates retains only the same accepted2logPstar as a shared symbolic base. Entire width/F0/swirl/radius log covers move into the formal offset. Signed coefficient intervals and all log intervals are copied explicitly to the common directed context. This handles the legacy bridge-width context while preserving full covers. No source factor is exponentiated, no endpoint/midpoint becomes a field value, and no radius cap is substituted for the actual radius.

The source family/axis datum must match. The2logPstar interval must be exactly the same accepted cover after context copying; differently rounded covers are rejected conservatively. Rebasing may widen dependencies and does not establish new joint source correlations. All C0 and Z rows are rebased individually to the same basis/ledger before addition.

## Actual serial C1 operator and original background

The native partition is [.12,.13,.14,.15] on original O2_slope, with dy=dcoordinate and three exact1/100 widths. Each cell uses the actual phase union at candidate N=1024, covering full periods rather than sampled phases, and the original five signed density C0/Z functions. Z endpoints, widths, radius phase and recovery rates are independent of Z.

For each original rate lambda=1,3/2,3/2,1,0:

```text
H_out = exp(-lambda*width)*H_in + I_cell
H_out_Z = exp(-lambda*width)*H_in_Z + I_cell_Z
```

Composing cells yields an affine operator whose incoming coefficient is exp(-lambda*.03), with a separately accumulated signed C0/Z contribution. The p channel has exact coefficient1; quiet cells preserve pressure memory. Positive decay factors remain formal even for unmaterializable huge widths. Both source signs and every original nonzero-V cross term are retained.

The right endpoint original background functions and separate P0/P0_Z are queried directly from the same unchanged source packet. The completion formula is

```text
own_history_at_right = original_background_history_at_right
                    + operator(actual_incoming_correction_at_local_left)
```

**The actual incoming correction at coordinate.12 is still missing.** The original inlet-to-.12 gap is not silently skipped or treated as quiet. No corrected total history is claimed from the affine operator alone. This is the next global assembly dependency.

The expanded cells have conservative whole-period signed covers that can span zero. Their existence does not prove a tight cancellation estimate or a net integral sign. Efficient signed oscillatory integration/cancellation remains required for quantitative Rc targets and cone margins.

```python
from lei_ren_part1_paper_compliant_current_native_C1_history_transfer import NativeC1HistoryTransfer

# c1_owner is the accepted same NativeDensityC1LocalIntegrals.
backend = NativeC1HistoryTransfer(c1_owner)
inlet = backend.inlet(Z=(-1, 1), N=1024)
P0_Z = inlet['functions']['P0_Z']
local = backend.serial(Z=('.49', '.51'),
    endpoints=['.12', '.13', '.14', '.15'], N=1024)
# Applying the local operator requires all five true incoming correction
# functions and Z rows at .12 in the same original source family.
```

Evidence: two actual inlet C1 queries/28 original rows; two native three-cell partitions/20 cumulative C0/Z affine increment rows; all six cells use actual full-period phase covers; two independent signed formal-rebase references, legacy-context entire-cover copy, tiny-scale/huge-quiet-decay retention, common-family/Pstar rejection, and30 independent analytic piecewise C0/Z comparisons with nonzero incoming memory and a quiet segment. Read-only review: **GPT-5.6 Luna / max**, no material mathematical blocker. 1047 working/index source hashes pass. Producer 21.844s; focused checker 25.890s on the already-live native source.

## Detailed next production tasks

- [x] **LEFT4c3-common-adjacent-cell-basis:** shared accepted Pstar² factor, directed complete-context copies, varying source/radius covers moved formally, signed C0/Z rows rebased before addition.
- [x] **LEFT4c3-local-serial-C1-Duhamel-transfer on[.12,.15]:** three true adjacent O2 cells, exact total width3/100, original rates, explicit incoming argument, preserved rate0 pressure memory.
- [x] **LEFT4c3-original-incoming-histories at actual r_minus:** all five original functions and Z rows on whole[-1,1] and[.49,.51]; original memory retained and separate checked correction initial zero.
- [x] **LEFT4c3-P0-and-P0_Z on original inlet/local right endpoint:** independent original datum/derivative, original absolute pressure C0/Z covers, no merger into p or double normalization.
- [x] **LEFT4c3-actual-chart-length-transfer backend/five declared cells (see successor); global route remains open:** use existing original radius offsets/Jacobians to construct true log-radius widths for each bridge, switch, macro, reshape, restoration, patch, O2 and O3 interval. Preserve microscopic width factors, long chart lengths and source-dependent endpoint definitions; quiet intervals still transport inherited correction memory.
- [ ] **LEFT4c3-supported-flat/quiet-segments:** derive exactly which regions have q=0 from the checked kappa/eta collar arguments. Carry all five incoming correction histories with the original decay; do not infer an uncomputed interval is quiet from a nearby point or a local flat slice.
- [ ] **LEFT4c3-full-Z-cutoff/angle-coverage:** resolve remaining bridge_first/O2_axial cutoff crossings and derivative coverage on fullZ, preserving tiny positive Delta/eta. Keep unresolved function domains explicit and original smooth cutoff behavior.
- [ ] **LEFT4c3-inlet-to-local-incoming-correction:** start from the checked actual inlet correction zero and integrate/transport the true signed original density functions through every intervening region to O2 coordinate.12. Return all five actual C0/Z incoming correction functions for the existing local affine operator; no skipped gap or invented initialization.
- [ ] **LEFT4c3-signed-oscillatory-cancellation:** exploit the original periodic phase/primitive identities to obtain whole-chart signed averages and integration-by-parts/cancellation bounds. Separate zero-mean leading terms from quadratic mean contributions; bound slow source variation and Z derivatives. Candidate full-period range hulls alone are too broad for tight terminal targets.
- [ ] **LEFT4c3-robust-global-common-coordinate-error:** quantify the dependency loss introduced by rebasing broad radius/amplitude log covers. Where bounds are too wide, use exact original factor/affine-radius correlations or certified subdivisions rather than choosing scalar caps or suppressing signed terms.
- [ ] **LEFT4c3-actual-cumulative-Rc-C1-functions:** propagate the actual incoming corrections through all remaining original charts to Rc, preserving source/phase/transfer/integration errors. Combine with original background histories and P0/P0_Z only after the correct endpoint argument is present.
- [ ] **LEFT4d-actual-Rc-targets:** build A,A_Z and divided(J-M)/mu from the same completed C1 histories and original pressure datum, keeping inverse-mu factors formal and source roles explicit.
- [ ] **LEFT4d-independent-repair-and-terminal-identities:** apply the reserved-band inverse, produce unique controls as functions of Z, and verify all five terminal identities. Continue changed radial velocity, absolute pressure and stress through downstream flatten/heat pieces.
- [ ] **LEFT4c2-higher-phase/primitive-jets + HIGH/CONT/OUTER/ENERGY/LEFT4e:** higher mixed recurrence/tails, same-function C4 seams, analytic pressure/exact heat/finite energy, derivative constants/common finite N and global admissible cone.
- [ ] **REC/WAVE/PHYS:** actual n-dependent coefficient recovery with independent repairs/smooth sum; oscillatory quadratic stress cancellation; full corrected Cartesian NS and multi-time scale/particle diagnostics.

Task completion requires exact code/result/check commits and actual function domains. Next work is the intervening actual chart transport and signed oscillation control needed to supply the local incoming correction, then actual Rc cumulative functions. Do not rerun unchanged local C1/q/phase stages. Every global completion gate remains false. Prior local derivative details are in the [signed-density/C1 handoff](CURRENT_NATIVE_DENSITY_C1_LOCAL_INTEGRALS_2026_10_07.md).

Scoped gate: `current_original_native_inlet_common_basis_and_adjacent_C1_Duhamel_operator_executed`. Five terminal moment identities, global repair/cone and coefficient recursion remain incomplete.
