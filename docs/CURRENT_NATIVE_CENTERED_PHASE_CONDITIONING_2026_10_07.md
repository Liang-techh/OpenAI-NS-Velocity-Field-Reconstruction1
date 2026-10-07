# Centered original phase and a*nu=v mixed derivative conditioning

Checked implementation: [4eb8b731](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/4eb8b73150a190e05603ff006076541e0c7b1d30). `NativeCenteredPhaseConditioning.route(Z)` rebuilds the same original inlet-to-Rc route on **24 true cells /17 charts** for Z=[-1,1] and Z=[.49,.51], for every integer N>=160. Centered phase residuals and the original relation a*nu=v remove a spurious large factor from A/B mixed inverse bounds. All five averaged target certificates improve in both domains. Maximum natural-log cap: approximately `1.106919531e+408906090034569677` -> `9.624298732e+408906090034569676`.

These are conservative logarithmic bounds, not actual error or field values. The direct certificate remains sharper at N=160 for every target. Actual source target evaluation, five controls, terminal closure and genuine scale recursion remain incomplete. The next production step connects exact N-dependent source functions to transport and repair; a sharper floor comparison is useful evidence, not a substitute for that interface.

## Same-source centered identities

Here v is the original loop parameter, distinct from the physical axial velocity V. Set nu=1+t0^2+2q^2, alpha=sqrt(a)*q, W=H-psi/2, J=t0*qP+q^2*W. Then

`T2=(nu-1)*psi+4J; S=nu*psi+4J=2pi*nu*phi; K=J/nu=(-b*q*P+a*q^2*W)/v`, with `v=a*nu`.

The original body has v=2+2eta and the transition has2<v<=2+2eta<=3. The alpha cutoff has no a denominator; native a/b/Delta/chi derivatives still enter its complete first/mixed chain. The derivative of v(kappa) is bounded through the original sigma derivatives; body, transition and exact flat branches remain the original smooth construction.

At the inverse, the FIXED-FREE-ANGLE residuals are

`L_i=-4nu*K_i`,

`L_yZ=-4*(nu_y*K_Z+nu_Z*K_y+nu*K_yZ)`.

Both first cross terms remain. Substituting the inverse phase equation before differentiating while holding the wrong variables fixed would lose one of them. The numerator -b*q*P+a*q^2*W and safe denominator v>=2 give complete K_i/K_yZ bounds. Weighted bq*P derivatives and aq^2*H derivatives retain q/chi correlations; no division by q or chi is used.

For F=1+t^2, psi_i=-4nu*K_i/F. In A_i and A_yZ, combine a*nu=v BEFORE bounding. This gives

`A_i=a_i*(phi-psi/(2pi))/2+v*K_i/(pi*F)`

and

`A_yZ=a_yZ*(phi-psi/(2pi))/2+(v*K_yZ+v_y*K_Z+v_Z*K_y)/(pi*F)-2v*t*(t_y*K_Z+t_Z*K_y)/(pi*F^2)+8v*nu*t*t_psi*K_y*K_Z/(pi*F^3)`.

All t slow derivatives above are fixed-free-angle derivatives. The original normalized curvature cap retains the final inverse term. No t bound such as |t|<=1 is assumed.

## Complete original B correction

Use M=-a*hatT1/(2pi)-b*phi, B=E*M/2. At fixed free angle, a*T1=-b*psi+2a*qP, so the first/mixed free terms are -b_i*(phi-psi/(2pi))-(a*qP)_i/pi. The full mixed inverse correction is

`2v*(1-t^2)*(t_y*K_Z+t_Z*K_y)/(pi*F^2)-8v*nu*(1-t^2)*t_psi*K_y*K_Z/(pi*F^3)+2t*(v*K_yZ+v_y*K_Z+v_Z*K_y)/(pi*F)`.

This identity comes from the complete original a*hatT1 product; it retains b sensitivity in t_y/t_Z. It does not replace the b-weighted terms of a different representation by a-weights. Independent symbolic checks verify this cancellation. The complete E*M product gives B_y/B_Z/B_yZ.

The original phase range gives |A|<=a/2<=3/2. Cauchy for T1/T2 gives |a*T1|/(2pi)<=sqrt(a*(v-a)), and b^2<=a*(v-a); hence |B|/E<=sqrt(a*(v-a))<=v/2<=3/2. Underlying E/V and all signed pressure/energy source terms are retained even on the flat modulation branch.

## Native route evidence and next source interface

True microscopic widths, one global fractional phase, original zero inlet, endpoint terms, all incoming histories, quiet pressure memory, P0/P0_Z and Rc A/A_Z/mu remain the same. Accepted source artifacts are immutable; only same-function bounds are intersected. The native cap table is for normalized N*r C1 at N=160, with directed full records in JSON.

| Z query | Target | Previous averaged log cap | Centered averaged log cap |
| --- | --- | --- | --- |
| whole_Z | M | 1.094044215e+408906090034569677 | 9.495545571e+408906090034569676 |
| whole_Z | D=(J-M)/mu | 1.094044215e+408906090034569677 | 9.495545571e+408906090034569676 |
| whole_Z | I | 1.086891261e+408906090034569677 | 9.424016038e+408906090034569676 |
| whole_Z | S | 1.094044215e+408906090034569677 | 9.495545571e+408906090034569676 |
| whole_Z | Cp | 1.106919531e+408906090034569677 | 9.624298732e+408906090034569676 |
| Z_interval | M | 1.094044215e+408906090034569677 | 9.495545571e+408906090034569676 |
| Z_interval | D=(J-M)/mu | 1.094044215e+408906090034569677 | 9.495545571e+408906090034569676 |
| Z_interval | I | 1.086891261e+408906090034569677 | 9.424016038e+408906090034569676 |
| Z_interval | S | 1.094044215e+408906090034569677 | 9.495545571e+408906090034569676 |
| Z_interval | Cp | 1.106919531e+408906090034569677 | 9.624298732e+408906090034569676 |

The first-bridge centered K_yZ log cap is `9.556345675e+408906090034569676`. Its remaining normalized inverse-curvature component is `1.108350121e+408906090034569677`. This is still the principal conservative loss. The actual SOURCE1/2/3 transport/target interface can proceed using the minimum of existing valid certificates and the exact conditional repair theorem; bounds must never become function coefficients.

## Focused evidence

- 960 C0/Z transport rows, 480 inherited rows, 40 exact quiet rows and 8 retained quiet pressure rows.
- 2576 signed source derivative terms, 112 weighted first-bridge terms and 40 target order rows.
- 6 independent symbolic alpha/A/M product identities, including complete M inverse sensitivity. 81 independent original-cutoff v-range/first/second derivative references.
- 80 independent original scalar-loop references cover signed-r, r crossing, transition and flat branches. Unchanged periodic kernel proofs/reference receipt are inherited with source hashes; no redundant kernel rerun was needed.
- Producer 69.547s; focused checker 77.968s on the warm original runtime; 1102 working/index dependency hashes. Read-only reviewer: **GPT-5.6 Luna / max**.

## Agent tasks and acceptance criteria

- [x] **COND6b3 centered phase residual bounds:** use K=J/nu=(-b*q*P+a*q^2*(H-psi/2))/v, alpha=sqrt(a)*q and v=a*nu on the original active support. Bind the body/transition v range2<=v<=3 and complete alpha/v/source mixed rules. Retain BOTH nu_y*K_Z and nu_Z*K_y in L_yZ, every inverse cross term and signed chi=p2/dstar. Whole fixed-free-angle K_i/K_yZ covers feed the same original inverse and five-target route.
- [x] **COND6b4 centered B bounds:** differentiate the ORIGINAL -a*hatT1/(2pi)-b*phi product. Cancel a*t0=-b in the fixed-free-angle product, then combine a*nu=v in its full inverse correction. Preserve b sensitivity inside t_y/t_Z and every inverse term. Independent symbolic identities verify the complete M_yZ correction; original scalar-loop references cover the resulting B_y/B_Z/B_yZ caps. These are function-domain bounds, not installed point controls.
- [x] **COND6b5 C0 normalization:** |A|<=a/2<=3/2; Cauchy bounds |B|/E<=sqrt(a*(v-a))<=v/2<=3/2. Preserve exact flat modulation zeros and nonzero underlying velocity. Intersect same-function bounds only; no source coefficient comes from a cap.
- [ ] **SOURCE1 exact N-dependent source transport:** connect the17 accepted original loop function graphs to the24 true-cell Duhamel route. Each cell must integrate the FULL signed delta_density_j(y,Z,N), not just its leading/remainder cap. Start with the actual zero correction at the original inlet; carry all preceding history, global fractional phase and original P0/P0_Z. Emit exact integral function nodes and complete first Z rules for all five Rc corrections. An interval coefficient or N^-2 cover must not be returned as a target function value.
- [ ] **SOURCE2 actual source/inverse evaluator:** add a source-derived point evaluator or a rigorously defined functional oracle for the original signed input DAG, phase inverse, A/B, deltaE/deltaV and densities. Preserve directed/exact source units and formal scales; do not take midpoints of saved covers. Test signed r, r=0, transition/flat branches, endpoints and nonzero E/V/b/p2. Keep bounds and function values as separate types; report precision limitations without changing the family.
- [ ] **SOURCE3 same-original target functions:** apply the current Rc transformation to the exact incoming functions: M, (J-M)/mu, I, S and Cp, using the same actual A/A_Z and mu and full signed divided-row numerator. Bind N as the same positive integer throughout the route. The existing conditional repair theorem applies to compatible C1 target functions; passing conservative covers alone does not construct them.
- [ ] **CONTROL1 actual five function controls:** solve B_exact(mu)*h+N*r+Q_exact(mu,h)/N=0 for the SOURCE3 functions with original five-bump/Gram functions. Use source-bound exact inverse/matrix functions, preserve incoming histories and distinguish h from -r. Give the actual convergent iteration/functional solution and its first Z derivative; do not replace it by a ball bound or sampled constant vector.
- [ ] **CONTROL2 finite-N and closure:** choose and justify one finite integer N compatible with actual repair/cone/higher-jet constraints, without materializing an astronomical integer unnecessarily. Prove all five terminal equalities and first Z derivatives throughout the required Z domain for the actual repaired field. Carry it through Rc->2Rc, matching/right collar and separate post-repair2Rc->Rb admission. Repair-only sufficient log conditions do not by themselves admit global N.
- [ ] **COND6c integrated variation:** where it improves the best retained certificate or practical precision, prove weighted G_y/G_yZ integrals/total variation on the original first bridge, keeping width, phase, endpoints and inherited memory. No sampled-only whole-domain claim.
- [ ] **COND4/5 partitions and seams:** derive no-gap source-defined body/transition/flat partitions and exact G/G_Z join identities before cancellation; keep all16 joins and the initial flat collar. Preserve competing certificates. Shared radii alone do not prove same-function endpoint cancellation.
- [ ] **COND6d retain best bounds:** compare direct and all averaged certificates at identical N over both Z domains. Current centered averaging remains worse at N=160 for every row. This comparison is diagnostic, not an independent requirement to postpone the actual SOURCE/CONTROL interface. Use whichever proved bound is tighter; separate cap improvement from measured actual error.
- [ ] **HIGH/OUTER/ENERGY:** actual higher slow/fast/physical jets and C4 seams, compatible pressure, whole-route stress/cone, exact heat exterior, physical tail energy and flat remainder. The C0/y/Z/yZ repair-bound stage does not admit these layers.
- [ ] **REC/WAVE/PHYS:** distinct n=1/n>=2 coefficient recovery with independent moment repairs, finite-order remainder/smooth sum, two actual oscillatory pulse families and mean lift, averaged stress cancellation, corrected Cartesian NS and measured scale recursion/material winding.


Scoped gate: `current_original_native_centered_phase_a_nu_correlation_and_five_target_bounds_executed`. Original centered mixed source bounds and five target covers are implemented. Actual target values/controls/terminal closure, point inverse/higher jets, one global finite N/cone, heat/energy, coefficient recursion, oscillatory stress cancellation and full corrected NS remain open.
