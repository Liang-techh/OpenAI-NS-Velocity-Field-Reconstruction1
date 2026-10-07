# Actual whole first-bridge C1 history covers

> Successor: [CURRENT_NATIVE_SECOND_BRIDGE_C1_HISTORIES_2026_10_07.md](CURRENT_NATIVE_SECOND_BRIDGE_C1_HISTORIES_2026_10_07.md) ([5e6f617c](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/5e6f617c296d50701fa74fedd79def6e1c05e91a)) now carries actual phase1 C0/Z correction functions through the full second bridge [1,2] to phase2, with original endpoint backgrounds and separate P0/P0_Z. Macro/downstream route and quantitative closure remain open.

Checked implementation: [c14d5cc7](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/c14d5cc75263869d79cb89c326ceab9d4f0ae9c4). The actual original inlet-to-phase1 route now has conservative five-moment correction and own-history C0/Z covers on whole Z=[-1,1] and Z=[.49,.51]. This completes the first bridge interval, including its genuine active cutoff crossings. The second bridge and remaining route to Rc are still missing. Wide covers do not prove small errors, terminal identities, a global common frequency or scale recursion.

## What is implemented

`NativeFirstBridgeC1Histories.first_bridge(Z, N=1024)` starts from the accepted original inlet sc/2 and carries its known correction initial condition through the flat collar to 3sc/4. It then integrates the entire remaining first-bridge coordinate interval [3sc/4,1]. Its true log-radius width is h_bridge*(1-3sc/4), retained as a nonzero formal factor. The original radius and phase definitions are unchanged.

The original q is bounded by treating its branches separately:

```text
kappa = a + b^2/a; Delta = kappa-2
q = sigma(1-Delta/eta)*sqrt((2eta-Delta)/(2a))  when Delta<eta
q = 0                                                when Delta>=eta
```

On active support, eta<=1/2 and kappa<2+eta imply 0<a<=3 and |b|<=3/2. Body Delta<=0 uses sigma=1; transition 0<Delta<eta uses the original sigma on [0,1] and its true first derivative; the flat branch has exactly q=q_Z=0. Full source a_Z, b_Z and Delta_Z rows remain present. The active restrictions bound functions only on their proved support and do not redefine field values.

The original periodic primitive and unique implicit inverse are covered for every fractional phase in [0,1], without selecting a phase or using samples as the native field. Original Fourier/Mobius identities give safe |W1|<=4pi and (1-r^2)*W2<=100pi, hence |A|<=159 and |B/E|<=100. Here E is the original Utheta/Pstar, so the B cover has the original normalization. Candidate N>=160 is required by this backend's |A|/N<1 bound; N=1024 is executed. This is a local cover condition, not an accepted global N.

Signed-r differentiation stays smooth through r=0. The positive original Poisson denominator is bounded without subtracting a rounded r from 1: 1-|r|>=1/[2(1+|u|)^2]. Fixed-phase first-Z derivatives include the original implicit angle derivative and every differentiated source factor. The cutoff is flat at both endpoints, so the branch union has no seam distribution terms. Higher jets are not claimed.

The existing signed-density first-Z kernels then retain all five nonlinear and cross terms. Original Duhamel masses integrate their whole-cell covers with the true width exactly once. The same correction C0/Z functions are carried from the known inlet; they are not reset at phase1. Original background moments and separate analytic P0/P0_Z are queried at the actual phase1 endpoint and added to obtain own histories. Every row uses the accepted common formal basis and arithmetic ledger.

These analytic covers are intentionally conservative and may be extremely wide. They establish a gap-free first-bridge function-domain enclosure and actual next-chart incoming data. They do not supply a tight phase point evaluator, useful terminal error sizes, oscillatory cancellation, admissible cone margins or a corrected NS solution.

## Evidence

- Two actual native queries: whole Z=[-1,1] and Z=[.49,.51]. **40** correction/own-history C0/Z rows replay identically.
- **144** independent comparisons with the original loop: A, B/Pstar and their implicit fixed-phase first-Z derivatives. Reference cases include signed r crossing zero and large positive/negative Poisson parameters.
- **36** independent original q/q_Z comparisons, including a single genuine body/transition/flat crossing. Tiny positive eta remains nonzero formally.
- The independent fixture uses exact rational interval arithmetic. Scalar references and finite differences are diagnostics; they never define native functions or certify higher jets.
- Read-only mathematical review: **GPT-5.6 Luna / max**, no material enclosure error found in this scope.
- **1056** working/index source hashes match. Production 21.375s; focused checker 23.860s on the same live native source.

## Agent tasks and acceptance criteria

Mark a task complete only after committing its implementation, result and focused checker. Record exact coordinate/Z domains, inherited incoming data, remaining gaps and the source commit. Work from this report and the current source; older sections are historical. Avoid rerunning unchanged inlet, local C1 and first-bridge stages.

- [x] **LEFT4c3-true-chart geometry:** all 17 original positive lengths, 16 exact radius seams and true-width signed C1 adapter on five declared cells; see predecessor.
- [x] **LEFT4c3-actual-initial collar:** sc/2 -> 3sc/4, whole-Z original correction/background/P0_Z retained.
- [x] **LEFT4c3-active first bridge, whole-Z C1 covers:** sc/2 -> phase1 with genuine cutoff/phase crossings and actual inherited C0/Z histories. This item does not mean quantitative closure.
- [x] **LEFT4c3-second bridge conservative C1 covers (see successor):** reuse the original whole-period C1 cover only after proving its a-positive lower bound on bridge_second [1,2]. Integrate its true h_bridge width, start from the known phase1 correction C0/Z functions, and compare original adjacent radius seams. No reset or omitted interval. Produce actual phase2 own histories on whole Z and the declared interval Z.
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

Scoped gate: `current_original_whole_first_bridge_cutoff_phase_C1_history_covers_executed`. Global inlet-to-Rc admission, terminal five moments, common finite N/cone, coefficient recursion and full NS remain false.
