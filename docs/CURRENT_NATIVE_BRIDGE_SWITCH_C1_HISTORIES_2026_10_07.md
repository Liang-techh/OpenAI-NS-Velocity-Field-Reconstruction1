# Actual inlet-to-R110 C1 histories through macro bridge and switch chain

> Successor: [CURRENT_NATIVE_MIDDLE_O2_INLET_C1_HISTORIES_2026_10_07.md](CURRENT_NATIVE_MIDDLE_O2_INLET_C1_HISTORIES_2026_10_07.md) ([caa69cab](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/caa69cabc652be58338ea8da12c7839d70e178af)) now carries actual R110 C0/Z correction functions through six complete middle/reference charts to the directly queried O2_slope0 inlet, retaining analytic patch exp(1), original backgrounds/P0_Z and seven same-function source joins. O2/O3-to-Rc and quantitative closure remain open.

Checked implementation: [bed4bb3b](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/bed4bb3b9b0ced79cfda2fca7860d02ee6ae9ffe). Actual five-moment correction and own-history C0/Z covers now extend without gaps from original sc/2 inlet to R110 through all six original inner charts. This stage adds complete bridge_macro, switch_first, switch_second and switch_power intervals to the accepted two microscopic bridges. Whole Z=[-1,1] and Z=[.49,.51] are executed at candidate N=1024. The route after R110 and quantitative terminal closure are incomplete.

## Original added route and true lengths

| Chart | Entire original native domain | True log-radius width |
|---|---|---|
| bridge_macro | fraction [0,1] | 4*logPstar+log(100/4)+1000-2*h_bridge |
| switch_first | phase [0,1] | h_switch |
| switch_second | phase [1,2] | h_switch |
| switch_power | fraction [0,1] | log(110/100)-2*h_switch |

`NativeBridgeSwitchC1Histories.route(Z,N=1024)` obtains the actual accepted phase2 correction C0/Z functions, then integrates these four full cells in order. Four exact original radius/periodic-phase seams are asserted. Original full-box signed roots, chart-uniform positive-a certificates and the accepted original cutoff/whole-period C1 backend are used on every cell. No active cutoff/phase interval is discarded.

The signed density rows already use ordinary y=log R. In particular, switch phase source derivatives have their microscopic inverse-width conversion from the existing adapter; macro and switch-power rows are ordinary y rows already. The true positive mass accounts for each width once. Macro and power widths retain their negative microscopic smoothing terms, which must not be dropped or replaced by local selector widths.

Each chart carries the five correction C0/Z functions as `exp(-lambda*width)*incoming + signed_integral`. The original endpoint background is then queried and added once for own histories. Analytic P0/P0_Z remains separate. Rate-zero pressure has exact unit incoming coefficient on every segment. All rows share the accepted original family, common formal basis and arithmetic ledger.

The same signed increments also form a reusable phase2-to-R110 affine C1 operator. Its final covers overlap the serially propagated covers; both enclose the same source functions. This overlap is a consistency check, not proof of small error or a tighter global bound.

## Evidence and limits

- **160** new correction/own-history C0/Z rows: four complete charts, two Z domains, 20 rows per chart/domain.
- **80** actual inherited C0/Z correction rows are retained across the four new chart transitions. The phase2 record is byte-identical to the accepted predecessor.
- **20** serial/composite final cover comparisons; four exact original adjacent radius/phase seams; true-width branch, positive-a, common-basis, P0 and pressure-memory guards.
- All eight added full-box native slow-q source queries returned enclosed. The general original scalar-loop, cutoff, true-width mass and nonzero-incoming affine theorems are inherited without separate old-stage reruns.
- Read-only mathematical review: **GPT-5.6 Luna / max**, no material enclosure or normalization error in this scope.
- **1064** working/index hashes pass. Production 41.765s; focused checker 43.782s using the same live original native source.

The source functions and actual correction decomposition are unchanged. Bounds may remain extremely wide; a box that contains zero does not prove a function vanishes. These results supply actual incoming functions at R110. They do not establish terminal five-moment identities, useful repair error sizes, original velocity/stress higher-derivative seam matching, a global finite N/cone, coefficient recursion or full corrected NS. Counting completed charts is not a percentage of those final requirements.

## Agent tasks and acceptance criteria

Mark a task complete only after committing its implementation, result and focused checker. Record exact coordinate/Z domains, inherited incoming data, remaining gaps and the source commit. Work from this report and the current source; older sections are historical. Avoid rerunning unchanged inlet, local C1 and first-bridge stages.

- [x] **LEFT4c3-true-chart geometry:** all 17 original positive lengths, 16 exact radius seams and true-width signed C1 adapter on five declared cells; see predecessor.
- [x] **LEFT4c3-actual-initial collar:** sc/2 -> 3sc/4, whole-Z original correction/background/P0_Z retained.
- [x] **LEFT4c3-active first bridge, whole-Z C1 covers:** sc/2 -> phase1 with genuine cutoff/phase crossings and actual inherited C0/Z histories. This item does not mean quantitative closure.
- [x] **LEFT4c3-second bridge conservative C1 covers (completed in this report):** reuse the original whole-period C1 cover only after proving its a-positive lower bound on bridge_second [1,2]. Integrate its true h_bridge width, start from the known phase1 correction C0/Z functions, and compare original adjacent radius seams. No reset or omitted interval. Produce actual phase2 own histories on whole Z and the declared interval Z.
- [x] **LEFT4c3-macro bridge conservative C1 covers (completed here):** carry those phase2 rows through bridge_macro [0,1], using 4*logP+log(100/4)+1000-2*h_bridge. Keep microscopic terms and the original selected constants; handle all source/cutoff crossings on the full cell.
- [x] **LEFT4c3-microswitch chain conservative C1 covers (completed here):** integrate switch_first [0,1], switch_second [1,2], switch_power [0,1] in order. Use the two h_switch lengths and log(110/100)-2*h_switch. Preserve inherited C0/Z pressure memory and original signed axial terms.
- [x] **LEFT4c3-reshape and inner route conservative C1 covers (see successor):** integrate reshape [0,1] and inner_reference [0,1], including the intervals outside the already accepted [.12,.15] local cell. Query original endpoint backgrounds, carry actual corrections and report radius seams explicitly.
- [x] **LEFT4c3-pressure restoration conservative C1 covers (see successor):** integrate axial_restore [0,1] and restore_buffer [-7,-6], preserving analytic P0/P0_Z as a separate datum. Verify ordinary-Z differentiation and the same source-family admission.
- [x] **LEFT4c3-actual patch conservative C1 covers (see successor):** cover native R/Rpatch in [1,e] with exact original log-coordinate lengths; apply the Jacobian once. Integrate all intervals from inherited left correction rows and verify the exact patch/Rh radius seam.
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



Next production starts from the known R110 correction C0/Z functions and covers reshape, inner_reference, axial_restore, restore_buffer, actual_patch and Rh_reference. The full actual patch has analytic endpoint exp(1) and exact true log-radius width 1; preserve that source expression rather than selecting a rounded endpoint or midpoint. Complete every intervening interval before claiming actual O2 incoming data.

Scoped gate: `current_original_inlet_to_R110_macro_bridge_switch_C1_history_covers_executed`. Global inlet-to-Rc admission, terminal five moments, common N/cone, recursion and full NS remain false.
