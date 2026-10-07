# Actual inlet-to-phase2 C1 histories through both microscopic bridges

Checked implementation: [5e6f617c](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/5e6f617c296d50701fa74fedd79def6e1c05e91a). Actual five-moment correction and own-history C0/Z covers now extend from original inlet sc/2 through phase1 and the complete second bridge [1,2] to phase2. Executed domains are whole Z=[-1,1] and Z=[.49,.51], at candidate N=1024. The macro bridge and downstream route are still open. These are conservative function enclosures, not quantitative terminal closure or scale recursion.

## Implemented route and ownership

`NativeSecondBridgeC1Histories.second_bridge(Z, N=1024)` requires the checked original first-bridge owner. It recomputes that accepted route within the same source-family/basis/ledger, then uses its actual phase1 **correction** C0/Z functions as incoming data. It never passes phase1 own/background histories to the correction operator and never resets the correction to zero.

The second bridge uses its original native coordinate [1,2] and true nonzero formal log-radius width h_bridge. A fresh signed-source query covers the entire [1,2] x Z box; the earlier point fixture 1.831 is not used to define functions or bounds. All five signed nonlinear density and first-Z rows use the accepted original cutoff body/transition/flat and whole-period C1 primitive backend. The positive Duhamel mass includes the true width once; the decay multiplies the inherited phase1 correction once.

For each original rate lambda and five-moment row j:

```text
delta H_j(phase2,Z)   = exp(-lambda*h_bridge)*delta H_j(phase1,Z)   + I_j(Z)
delta H_j_Z(phase2,Z) = exp(-lambda*h_bridge)*delta H_j_Z(phase1,Z) + I_j_Z(Z)
own H_j(phase2,Z)     = original H_j(phase2,Z) + delta H_j(phase2,Z)
```

The width/endpoints are independent of Z. The rate-zero pressure row retains exact coefficient 1 on both inherited C0/Z functions. Original endpoint backgrounds and their ordinary Z rows are queried at bridge_second phase2. P0/P0_Z remain separate analytic data and are not double-added to a moment history.

## Uniform positive denominator and original seams

The original lower a>=h_bridge/(2*Kmax)>0 is uniform over bridge_second [1,2] and Z[-1,1]. `current_inner_relaxed_inputs.bounds()` certifies Ra<R<=110, all Z[-1,1], and derives the second bridge through convex interpolation of h_bridge*Dbar to 4/5. `current_generic_shear_source_bounds.positive_denominators()` binds this analytic inequality to the original source formulas. The downstream chart certificate has `source_function_positivity_not_inferred_from_saved_box=true`.

The new stage reuses that chart-specific a lower and still queries all varying signed C0/Z roots on the complete box. The previously checked exact original radius/Jacobian identities bind the phase1 seam and the phase2 endpoint; no absolute-radius subtraction or new selected source constant is introduced. Those identities alone do not claim new high-order velocity seam checks.

## Focused evidence and limitations

- **40** actual phase2 correction/own-history C0/Z rows replay identically across the two Z domains.
- **20** inherited nonzero phase1 correction C0/Z rows are checked against the committed first-bridge records. Every new row obeys the actual nonzero-incoming affine relation.
- Original rate-zero pressure memory, common basis/ledger, separate phase2 background/P0_Z, positive microscopic width, full coordinate domain and no global gate overclaim are checked.
- The accepted original loop/cutoff receipt supplies **144** scalar primitive/implicit derivative and **36** q/q_Z comparisons. These references were inherited, not rerun for this stage. The original true-width/nonzero-incoming affine theorem is also inherited.
- Read-only source and implementation review: **GPT-5.6 Luna / max**, no material blocker found in this scope.
- **1060** working/index source hashes match. Production 27.562s; focused checker 29.453s on the same live native source.

The derivative covers may be extremely wide. This stage establishes actual incoming data at phase2 and two gap-free microscopic bridge intervals. It does not install a persistent multi-cell operator cache, a tight inverse point evaluator, small moment errors, terminal repair controls, a global common N, stress cone, higher-order coefficient recursion or full corrected NS. Tight correlated/oscillatory integration is still required. Merely counting completed charts is not a percentage of the final analytic construction.

## Agent tasks and acceptance criteria

Mark a task complete only after committing its implementation, result and focused checker. Record exact coordinate/Z domains, inherited incoming data, remaining gaps and the source commit. Work from this report and the current source; older sections are historical. Avoid rerunning unchanged inlet, local C1 and first-bridge stages.

- [x] **LEFT4c3-true-chart geometry:** all 17 original positive lengths, 16 exact radius seams and true-width signed C1 adapter on five declared cells; see predecessor.
- [x] **LEFT4c3-actual-initial collar:** sc/2 -> 3sc/4, whole-Z original correction/background/P0_Z retained.
- [x] **LEFT4c3-active first bridge, whole-Z C1 covers:** sc/2 -> phase1 with genuine cutoff/phase crossings and actual inherited C0/Z histories. This item does not mean quantitative closure.
- [x] **LEFT4c3-second bridge conservative C1 covers (completed in this report):** reuse the original whole-period C1 cover only after proving its a-positive lower bound on bridge_second [1,2]. Integrate its true h_bridge width, start from the known phase1 correction C0/Z functions, and compare original adjacent radius seams. No reset or omitted interval. Produce actual phase2 own histories on whole Z and the declared interval Z.
- [ ] **LEFT4c3-macro bridge:** carry those phase2 rows through bridge_macro [0,1], using 4*logP+log(100/4)+1000-2*h_bridge. Keep microscopic terms and the original selected constants; handle all source/cutoff crossings on the full cell.
- [ ] **LEFT4c3-microswitch chain:** integrate switch_first [0,1], switch_second [1,2], switch_power [0,1] in order. Use the two h_switch lengths and log(110/100)-2*h_switch. Preserve inherited C0/Z pressure memory and original signed axial terms.
- [ ] **LEFT4c3-reshape and inner route:** integrate reshape [0,1] and inner_reference [0,1], including the intervals outside the already accepted [.12,.15] local cell. Query original endpoint backgrounds, carry actual corrections and report radius seams explicitly.
- [ ] **LEFT4c3-pressure restoration:** integrate axial_restore [0,1] and restore_buffer [-7,-6], preserving analytic P0/P0_Z as a separate datum. Verify ordinary-Z differentiation and the same source-family admission.
- [ ] **LEFT4c3-actual patch:** cover native R/Rpatch in [1,e] with exact original log-coordinate lengths; apply the Jacobian once. Integrate all intervals from inherited left correction rows and verify the exact patch/Rh radius seam.
- [ ] **LEFT4c3-Rh/O2 inlet:** integrate Rh_reference [-5,0] and O2_slope [0,1]. In particular supply actual incoming functions at O2 coordinate .12 before applying the existing [.12,.15] serial operator.
- [ ] **LEFT4c3-O2 axial and buffer:** cover O2_axial [0,1] with its exp(M)-1 true width and whole-Z cutoff crossings. Then cover O2_buffer [0,11], preserving the [0,9]/[9,11] admission distinction. Local cells cannot fill intervening gaps.
- [ ] **LEFT4c3-O3/right collar:** integrate the required O3_slope_mu and O3_power route to the actual Rc endpoint. Prove quiet pieces on their original support; carry inherited rate-zero pressure C0/Z memory even when local density is zero.
- [ ] **LEFT4c3-tight first-bridge bounds:** measure the widths of these covers and the repair-target budgets. Restore source correlations and use certified partitions or analytic phase averaging where broad denominator/derivative hulls dominate; do not cap field amplitudes or replace them with bound coordinates.
- [ ] **LEFT4c3-signed oscillatory integration:** derive phase primitive/zero-mean integration-by-parts bounds with quadratic means and first-Z slow variation separated. Bound truncation and boundary terms without enumerating astronomically many cycles.
- [ ] **LEFT4c3-global Rc histories:** complete every original inlet-to-Rc interval and combine correction rows with original endpoint backgrounds/P0_Z. Require function-domain C0/Z evidence and an explicit no-gap route ledger.
- [ ] **LEFT4d-terminal control and five identities:** derive A/A_Z and divided (J-M)/mu from completed histories using the reserved independent inverse and formal inverse-mu. Verify all five terminal identities as Z functions; sampled values are insufficient.
- [ ] **HIGH/common N/cone:** extend genuine higher jets and same-function C4 seams, then establish one finite N and stress-cone margins on all required regions. The current N=1024 first-bridge cover is not this admission.
- [ ] **OUTER/ENERGY/flat remainder:** recover compatible analytic preheat pressure, exact heat exterior, physical finite energy and flat-remainder decay from the accepted matched background.
- [ ] **REC/WAVE/PHYS:** implement n-dependent coefficient recovery, independent repairs and smooth summation; then quadratic stress-canceling oscillations, corrected Cartesian NS residuals and measured multi-time dynamics.


Next production starts at the known phase2 correction rows and integrates bridge_macro [0,1], then the switch chain, preserving the exact original route. Do not run accepted old stages as separate verification jobs; integrate them only when needed to build the same-source actual incoming functions. A reusable cached route/operator may improve this cost later.

Scoped gate: `current_original_inlet_to_phase2_second_bridge_C1_history_covers_executed`. Global inlet-to-Rc admission, terminal five moments, common finite N/cone, coefficient recursion and full NS remain false.
