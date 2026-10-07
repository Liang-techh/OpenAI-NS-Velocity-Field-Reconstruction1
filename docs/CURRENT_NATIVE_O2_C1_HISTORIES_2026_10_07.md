# Actual original complete O2 C1 histories

Checked implementation: [27c21b19](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/27c21b19b41aa7b5da9a5979fe0faf08b698dee9). Actual five-moment correction and own-history C0/Z covers now extend from the original sc/2 inlet through the inner/middle route and **the complete O2 slope, axial and buffer charts**, ending at the directly queried O3 transition coordinate0. Whole Z=[-1,1] and Z=[.49,.51] use candidate N=1024. O3-to-Rc and quantitative terminal matching remain incomplete.

## Original complete O2 route

| Original chart/cell | Native domain | True log-radius width |
|---|---|---|
| O2_slope initial | [0,.12] | .12 |
| O2_slope three adjacent local cells | [.12,.13], [.13,.14], [.14,.15] | .01 each |
| O2_slope remainder | [.15,1] | .85 |
| O2_axial | [0,1] | exp(Md)-1 |
| O2_buffer shared admitted part | selector[0,9] | 9 |
| O2_buffer remaining source part | selector[9,11] | 2 |

`NativeO2C1Histories.route(Z,N=1024)` carries the accepted actual O2 coordinate0 correction C0/Z functions through these eight cells. The actual correction at .12 supplies the three-cell [.12,.15] operator; its composite output is compared with the same serial route. The complete eight-cell operator is separately compared with its serial output. Widths are the original log-radius lengths, with the source Jacobian applied once. Every endpoint background and separate analytic P0/P0_Z is queried from the same original family. The terminal background is queried at O3_slope_mu coordinate0 and bound to the accepted buffer/transition source-function join.

The buffer selector maps to shared_offset=selector-11. The smaller relaxed admission shared_offset[-11,-2] covers selector[0,9]. Selector[9,11] is source-covered and integrated here, but no new cone admission is inferred. The original reference/slope, slope/axial, axial/buffer and buffer/transition source joins are bound separately from radius identities.

## Improved original periodic C1 bounds

The field, source roots, cutoff and phase inverse are unchanged. A new bound adapter keeps the signed parameter r=u/sqrt(1+u^2) and uses the exact Mobius identity

`D^2-(1-r^2)^2*sin(psi)^2=((1+r^2)*cos(psi)-2*r)^2`, where D=1-2*r*cos(psi)+r^2.

This gives |theta_u|<=2. Analytically split Fourier/Mobius bounds cover both |r|<=1/2 (including zero) and |r|>1/2, with |W1_u|<=16*pi and |(sW2)_u|<=1000*pi. Full fixed-angle parameter derivatives and the implicit fixed-fractional-phase psi_Z term are retained. The original identity A=a/2*(phi-psi/(2*pi)) and the correlated B angle term avoid large Poisson denominator powers. Parameter-sensitive C0 bounds intersect accepted universal caps; numerical caps bound functions and never define their values.

The original body/transition/flat cutoff union is retained. If the entire original source box proves q=q_Z=0, then A/B and all five density C0/Z increments are exactly zero. Inherited correction histories are still transported; rate-zero pressure retains coefficient1. Neither small q nor a selected sample triggers this branch.

## Focused evidence and limits

- **320** new correction/own C0/Z rows and **160** inherited rows from eight cells on two Z domains.
- **20** actual .12-.15 local and **20** complete O2 serial/composite comparisons; **20** directly queried O3 transition0 own C0/Z rows.
- **144** independent original-loop value/fixed-phase Z references and **36** cutoff q/q_Z references, with signed-r zero crossings, large positive/negative u and all cutoff branches. These fixtures diagnose bounds and do not define the native field.
- Exact-flat primitive and all five density C0/Z zero checks; original widths, source joins, P0/P0_Z, common family/bases/ledger and pressure memory retained.
- Read-only model metadata: **GPT-5.6 Luna / max**. Review PASS within first-Z/O2 scope; checker receipt is separate from producer output.
- **1073** working/index hashes pass. Production 72.407s; focused checker 75.203s on the same live native source.

These are conservative function-domain covers. The new bounds improve this O2 segment but do not retroactively tighten the inherited inner/middle histories. No claim of small terminal moment error, five terminal Z identities, global finite N, stress cone, higher jets, heat/energy, flat remainder, coefficient recursion or corrected NS is made.

## Agent tasks and acceptance criteria
Mark a task complete only after committing its implementation, result and focused checker. Record exact coordinate/Z domains, inherited incoming data, remaining gaps and the source commit. Work from this report and the current source; older sections are historical. Avoid rerunning unchanged inlet, local C1 and first-bridge stages.

- [x] **LEFT4c3-true-chart geometry:** all 17 original positive lengths, 16 exact radius seams and true-width signed C1 adapter on five declared cells; see predecessor.
- [x] **LEFT4c3-actual-initial collar:** sc/2 -> 3sc/4, whole-Z original correction/background/P0_Z retained.
- [x] **LEFT4c3-active first bridge, whole-Z C1 covers:** sc/2 -> phase1 with genuine cutoff/phase crossings and actual inherited C0/Z histories. This item does not mean quantitative closure.
- [x] **LEFT4c3-second bridge conservative C1 covers (completed in this report):** reuse the original whole-period C1 cover only after proving its a-positive lower bound on bridge_second [1,2]. Integrate its true h_bridge width, start from the known phase1 correction C0/Z functions, and compare original adjacent radius seams. No reset or omitted interval. Produce actual phase2 own histories on whole Z and the declared interval Z.
- [x] **LEFT4c3-macro bridge conservative C1 covers (completed here):** carry those phase2 rows through bridge_macro [0,1], using 4*logP+log(100/4)+1000-2*h_bridge. Keep microscopic terms and the original selected constants; handle all source/cutoff crossings on the full cell.
- [x] **LEFT4c3-microswitch chain conservative C1 covers (completed here):** integrate switch_first [0,1], switch_second [1,2], switch_power [0,1] in order. Use the two h_switch lengths and log(110/100)-2*h_switch. Preserve inherited C0/Z pressure memory and original signed axial terms.
- [x] **LEFT4c3-reshape and inner route conservative C1 covers (completed here):** integrate reshape [0,1] and inner_reference [0,1], including the intervals outside the already accepted [.12,.15] local cell. Query original endpoint backgrounds, carry actual corrections and report radius seams explicitly.
- [x] **LEFT4c3-pressure restoration conservative C1 covers (completed here):** integrate axial_restore [0,1] and restore_buffer [-7,-6], preserving analytic P0/P0_Z as a separate datum. Verify ordinary-Z differentiation and the same source-family admission.
- [x] **LEFT4c3-actual patch conservative C1 covers (completed here):** cover native R/Rpatch in [1,e] with exact original log-coordinate lengths; apply the Jacobian once. Integrate all intervals from inherited left correction rows and verify the exact patch/Rh radius seam.
- [x] **LEFT4c3-Rh reference to O2 inlet:** full Rh_reference [-5,0], actual O2_slope0 background and same-function source join; actual incoming correction C0/Z functions installed.
- [x] **LEFT4c3-O2 slope and local incoming conservative C1 covers (completed here):** integrate O2_slope [0,1] from these actual coordinate0 correction functions; explicitly supply coordinate .12, then apply the existing [.12,.15] operator and continue to1 without a gap.
- [x] **LEFT4c3-O2 axial and buffer conservative C1 covers (completed here):** cover O2_axial [0,1] with its exp(M)-1 true width and whole-Z cutoff crossings. Then cover O2_buffer [0,11], preserving the [0,9]/[9,11] admission distinction. Local cells cannot fill intervening gaps.
- [ ] **LEFT4c3-O3/right collar:** start from the actual O3_slope_mu coordinate0 correction functions. Integrate O3_slope_mu [0,1], then original power offsets [0,1] to r_plus and [1,2] to Rc=Rw*exp(2). Rc power phase is 2/Tw, not1. Preserve the original right collar offsets [.75,1.25] within the power chart. Prove quiet pieces on their original support; carry inherited rate-zero pressure C0/Z memory even when local density is zero.
- [x] **LEFT4c3-O2 parameter-uniform periodic bounds:** original signed-r Mobius/Fourier first-Z bounds installed; exact-flat primitive/density increments are zero without resetting memory. This is a bounds improvement, not a new field definition.
- [ ] **LEFT4c3-tight first-bridge and inherited bounds:** measure the widths of these covers and the repair-target budgets. Restore source correlations and use certified partitions or analytic phase averaging where broad denominator/derivative hulls dominate; do not cap field amplitudes or replace them with bound coordinates.
- [ ] **LEFT4c3-signed oscillatory integration:** derive phase primitive/zero-mean integration-by-parts bounds with quadratic means and first-Z slow variation separated. Bound truncation and boundary terms without enumerating astronomically many cycles.
- [ ] **LEFT4c3-global Rc histories:** complete every original inlet-to-Rc interval and combine correction rows with original endpoint backgrounds/P0_Z. Require function-domain C0/Z evidence and an explicit no-gap route ledger.
- [ ] **LEFT4d-terminal control and five identities:** derive A/A_Z and divided (J-M)/mu from completed histories using the reserved independent inverse and formal inverse-mu. Verify all five terminal identities as Z functions; sampled values are insufficient.
- [ ] **HIGH/common N/cone:** extend genuine higher jets and same-function C4 seams, then establish one finite N and stress-cone margins on all required regions. The current N=1024 first-bridge cover is not this admission.
- [ ] **OUTER/ENERGY/flat remainder:** recover compatible analytic preheat pressure, exact heat exterior, physical finite energy and flat-remainder decay from the accepted matched background.
- [ ] **REC/WAVE/PHYS:** implement n-dependent coefficient recovery, independent repairs and smooth summation; then quadratic stress-canceling oscillations, corrected Cartesian NS residuals and measured multi-time dynamics.

Next production starts with the actual O3_slope_mu coordinate0 corrections. Complete the transition and power offsets0->1->2 to Rc, with exact original flat support and pressure memory. The original right collar [.75,1.25] is inside the power chart. Rc->2Rc has length log2; the outer 2Rc->Rb admission remains open. Then tighten signed oscillatory/slow bounds and construct terminal repair controls; broad covers alone cannot prove closure.

Scoped gate: `current_original_inlet_to_O2_exit_C1_history_covers_and_parameter_bounds_executed`. Global inlet-to-Rc admission and all terminal/common-N/cone/recursion/full-NS gates remain false.
