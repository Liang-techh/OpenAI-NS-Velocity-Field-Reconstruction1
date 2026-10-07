# Current original middle relaxed inputs: R110 through Rm

**Successor:** the current whole original implicit patch condition through Rh, covering supports and retained partial-history gaps, is implemented in [CURRENT_GENERIC_SHEAR_PATCH_RELAXED_2026_10_07.md](CURRENT_GENERIC_SHEAR_PATCH_RELAXED_2026_10_07.md). Whole source norm/scales, changed loop/transport/repair/N and true recursion remain open.

Checked implementation: [46b538f5](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/46b538f507e8a76e6a2f0aa386077f80ff293f9d). **LEFT4c1-gates-b-reshape, reference, restore and buffer are implemented on the current original analytic source.** Together with the existing whole open inner exit, the original relaxed-input gates now cover Ra<R<=Rm through the unchanged source functions. The active five-bump patch is the next original middle condition gap. Whole-upstream assembly, p1/p2 derivative norms/conservative scales, changed generic loop/history transport, new repair/common N, global strict tensor cone and coefficient recursion remain open. The full long-term goal stays active.

## Current source functions and executable domains

`CurrentReshapeRelaxedInputs().query(chart, coordinate, Z=(-1,1))` covers `reshape` and `inner_reference`, each selector phase in[0,1]. The physical domain is 110<=R<=Rz=exp(-8)*Rref, including the actual variable reshape on110..Rsh and the original power continuation Rsh..Rz. `CurrentRestoreRelaxedInputs().query(...)` covers `axial_restore` phase[0,1] and `restore_buffer` offset[-7,-6], Rz<=R<=Rm=exp(-6)*Rref. They return analytic source-domain certificates, not velocity point values. The active patch selector is rejected by both APIs.

The new modules consume the existing typed full signed generic inputs and current inner theorem. They bind actual R110 B/V/m/pressure inputs, exact source callables/ASTs, all five history ODEs and absolute pressure equation, same positive kernels, original T=400*A, exact Rsh six-function equality and135 mixed4 rows, the current correlated E=V110-4Z, unchanged restoration kernels and exact Rz/restore/buffer/Rm source interfaces. The pressure source, datum and KN/Kp are the same analytic pressure family. Cap intervals remain covers of the original source functions; none is selected as its defining value.

No historical finite-fixture cone is transferred. Parameter-family equality is supplemented by actual source/function identities and the exact primitive-to-stress equations. Neither producer/checker constructs an ancestor graph or materializes a radius, width, Pstar or microscopic source amplitude.

## Whole reshape/reference angular barrier

For the original log swirl source

    log(Utheta/Pstar)=B(Z)*(1-sigma(y/T))-log(1+Z^2)+y/10-logCstar-logPstar,
    T=400*A, ||B||C2<=2*A,

the actual a=.8+2B*sigma'(y/T)/T has |a-.8|<=.08; the conservative admitted cover is[.7,.9]. The actual V=V110 is y-constant, so b=0. The source A/T cancellation is done before enclosing; independent broad A and T boxes do not define this shear.

Let **Qangular** denote the paper angular inertial quantity, not the radial velocity recovery Q. Complete original angular inertial sectors give

    p1=Itheta/F=R*Qangular/L, L=1-delta*Z^2,
    D_y Qangular+(2-a/2)*Qangular=SQ,
    D_y p1+(1-a/2)*p1=R*SQ/L,
    SQ=-W*(1-a/2)-delta*(1-2ZV)/2-Hv*zeta,
    W=1-A_operator(Mz/R), Hv=(1-delta)*Z/2+(1-Z^2)*V.

The current exact mean is the positive Z-independent average `(110/R)*Mz110/110+(1-110/R)*V110`. The same current inlet velocity/mean C2 error budgets retain |W-(-3+4delta Z^2)|<=10epsilon0. Zeta convexly interpolates the current entry gradient and -2Z/(1+Z^2). With V110=4Z+e, |e|<epsilon0, the target Hv*zeta splits into a nonpositive leading term and `-2Z*(1-Z^2)*e/(1+Z^2)<=2epsilon0`. This proves the scalar max bound without losing the signed leading term.

The directed SQ lower bound is greater than7/5. At p1=3 the derivative is at least110*SQ-3>0. The same current R110 inlet has p1>3; hence the entire reshape/reference path has p1>3, H0=p1>3, kappa=a<1 and t0=0. Integration for one logarithmic unit gives Qangular>1/2 thereafter. The query API does not give that half-bound at the reshape inlet or on a box including it; every inner_reference source point is after the long phase.

## Full signed restoration/buffer condition

The actual restoring velocity is V=4Z+E*(1-sigma(t)), E=actualV110-4Z, with the original axial5 E source, flat cutoff and full inherited kernels. The original swirl stays a=.8, zeta=-2Z/(1+Z^2). After restoration V=4Z and b=0 exactly. The angular ODE remains valid for arbitrary V_y; Z-independent cutoff weights preserve the V/mean error budgets. The p1>3 and Qangular>1/2 crossing barriers therefore continue from Rz through Rm.

The full signed axial inertial source is retained through the paper N:

    Iz=sqrt(R/2)*N/L, p2=Iz/F=R*N/(L*Utheta),
    J=N/Utheta^2, b=2V_y/Utheta,
    bw=b*p2/p1=2J*V_y/Qangular,
    H0=p1*(a-bw)/a.

The producer independently replays all original raw axial pressure, full energy, mean and nonlinear meridional sectors into this N. It does not omit p2 or replace it by shear. The same original source satisfies the exact arbitrary-V_y axial ODE, pressure C1 bound Kp, inherited R110 N norm, reference-radius decay gate and KN budget. Thus |N|<=KN*Pstar^2 after Rz. Since Utheta/Pstar>=exp(-.8)/2, |J|<=20KN and |V_y|<=16epsilon0, one obtains |bw|<=1280KN*epsilon0. The original Pstar>=1 log gate and b*Pstar bound give kappa<1. The directed uniform H0-2 lower bound is positive, using p1>3 and a-|bw|>0.

The changed whole-family loop, repair and common finite N are not implied by this original source admission. Strict completed tensor registry counts remain unchanged; no new strict regions are admitted. These modules also do not supply all p1/p2 slow or axial derivative norms for generic loop scales.

## Focused evidence

Reshape/reference: 36 actual source conditions, 13 original AST bindings, 16 exact source/ODE/control identities, 60 coefficients from independently integrated moderate full physical primitives, 3 positive crossing margins, 4 source queries and 5 invalid-domain guards. The independent fixtures test units/formulas; the exact source equations and analytic budgets prove the whole source domain.

Restoration/buffer: 19 actual source conditions, 7 AST bindings, 10 full signed identities, 6 positive pressure/inherited-N norm gates, 5 positive relaxed margins, two original normalization budgets and six invalid-source/active-patch guards. Working/index audit matched 972 bound hashes. Read-only reviewer: GPT-5.6 Luna / max. Old scoped N>=68,533,403 remains historical and cannot select the new whole-loop N.

## Next production tasks

- [x] **LEFT4c1-gates-a:** all16 original typed source covers export full I/F, S/F and division-free signed relaxed invariants.
- [x] **LEFT4c1-gates-b-inner:** exact current source attachment and analytic Ra<R<=110 relaxed input; Ra zero endpoint separate.
- [x] **LEFT4c1-gates-b-reshape/reference:** actual variable B/T, full histories and positive mean source; current angular SQ/p1/Q barriers on110..Rz.
- [x] **LEFT4c1-gates-b-restore/buffer:** actual whole restoring E/pressure/N functions, signed bw and H0 condition onRz..Rm. Stops before active patch.
- [ ] **LEFT4c1-gates-b-patch-source:** bind the actual five-bump H,V functions, common pressure, full inherited moments and repair amplitudes. Retain D_y=xD_x, a=1-2xH_x/H and b=2xV_x/Utheta; do not replace active bumps by a quiet power field.
- [ ] **LEFT4c1-gates-b-patch-supports:** establish source conditions on each active support49..51,59..61,69..71, including complete signed I/F sectors and derivatives. Use actual same-family defect/correction bounds to derive a>0,H0>2 and any required signed D/Q branch; endpoint flatness is insufficient.
- [ ] **LEFT4c1-gates-b-patch-gaps:** carry all five original moment histories and pressure across quiet gaps and after supports. The velocity may be power there, but its retained moments differ from the unpatched reference continuation. Prove the source inequality for those actual histories.
- [ ] **LEFT4c1-gates-c-norms:** derive full source p1,p2,t0 and required slow/axial derivative bounds across actual domains. Keep radius and amplitude factors in exact logarithmic/factored form; do not invert denominator covers containing zero or materialize tiny a_min.
- [ ] **LEFT4c1-gates-c-scales:** combine lower a, interior H0 margin, strict right-edge excess and p/t0 norms into the generic conservative scales and an executable full current-source input bundle. Preserve original width correlation.
- [ ] **LEFT4c1-gates-d:** reserve the actual power interval for new terminal repair; derive shared loop support from strict inner collar/right edge with valid quiet strips. Do not taper a loop so its v>2 condition disappears.
- [ ] **LEFT4c1-live-use:** exercise selected arbitrary source queries when an existing checked owner graph is available; do not rebuild every ancestor merely for a cache-wrapper check. Core saved rho[.001,4] still does not admit the axis.
- [ ] **LEFT4c2-phase/jets:** install one N*log(R/r_minus) through original selectors; construct factored q/r/inversion/zero-mean primitives and current whole-cell increments through mixed4 with phase-held slow derivatives.
- [ ] **LEFT4c3-transport/tensor:** propagate all five own changed-history defects over exact formal widths, source seams and quiet gaps; supply explicit factor-basis rebase identities, unchanged P0 and physical Cartesian/time tensor exports.
- [ ] **LEFT4d-repair:** reconstruct the new five-bump functional Jacobian, inverse, nonlinear error/uniqueness bounds and all five terminal identities as functions of Z.
- [ ] **LEFT4d-new-N / e:** derive whole-family loop/repair perturbation bounds, select a new full finite N and close the whole modified stress cone. The old scoped N is not a substitute.
- [ ] **CONT / ENERGY:** complete changed source/function interfaces and finite energy over the actual space/time domain and exact heat exterior.
- [ ] **REC / WAVE / PHYS:** true n-dependent coefficient equations and independent repairs, smooth sum, mean/oscillatory averaged stress cancellation, corrected physical uvw/p/f, Cartesian residual and shrinking/elongation/winding diagnostics.

New true gates: `current_original_R110_Rz_relaxed_generic_input_certified` and `current_original_Rz_Rm_relaxed_generic_input_certified`. All whole-upstream/current-loop/changed-moments/new-repair/common-N/global-cone/actual-recursion/full-corrected-NS gates remain false.
