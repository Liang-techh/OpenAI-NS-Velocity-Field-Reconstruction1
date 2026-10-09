# Successor: actual Rm leading five-defect inverse and partial patch

Checked source [1d4903c6](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/1d4903c652476c44ac0d579b71c09e8971b42876); evidence/tasks [CURRENT_ORIGINAL_RM_DEFECT_PATCH_INVERSE_2026_10_09.md](CURRENT_ORIGINAL_RM_DEFECT_PATCH_INVERSE_2026_10_09.md). Actual Rm rows now feed canonical normalized leading repair Z0..5 and partial velocity/pressure/primitives at native0,.5. Full radial mixed4, whole-axis and finite-N Rc/all24/global N remain open. Historical text below is preserved.

---

# Actual Rsh reference / axial restoration functions through Rm

Checked source [0d365bc0](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/0d365bc065756ad5a9d6312d06bf2220b9191471). Full reconstruction **ACTIVE / INCOMPLETE**. The actual long-reshape terminal six histories now feed original reference transport, full-cutoff axial restoration and postrestore through the true Rm inlet. Actual normalized histories, Q, physical velocity/pressure/five primitives and mixed y/Z derivatives through total order4 are available at native Z=0,.5. All incoming tails and original analytic P0 remain. Finite-N density integrals, actual active patch, whole-axis source closure and actual n-recursion remain open.

Accepted predecessors: [actual long-reshape endpoint](CURRENT_ORIGINAL_LONG_RESHAPE_ENDPOINT_2026_10_09.md), [both switches / R110](CURRENT_ORIGINAL_SECOND_SWITCH_R110_2026_10_09.md), first switch, factored R100 interface, actual macro F/V and six moments, full finite terminal kernels, original reference/O2 point/integral/partial histories and general fixed-phase Z calculus. New functions reuse accepted caches and scalar integrals, without old ancestor producer/full-checker execution.

## True source geometry and Rsh inlet

The native factory rehydrates the actual long-reshape terminal shapes and V/P0 from its accepted report, in the same source basis and ledger. The original normalized P0 row must equal the canonical upstream exact record; overlap is insufficient. The same fixed T=400Abar tuple is required. E is the actual V_Rsh-4Z, with no fitted j/error cap substituted as the field value.

    Rsh=110exp(T); Rref=110exp(10(logCstar+logPstar))
    Rz=Rref exp(-8); Rrestore=Rref exp(-7); Rm=Rref exp(-6)
    gap=10(logCstar+logPstar)-T-8 >0; dy=dlogR.

The actual Rsh shapes are centered before transport:

    mean_error=mean-4Z; angular_error=theta-5/8
    mixed_error=theta_z-4Z*theta
    axial_square=axial-8Z*mean+16Z^2
    swirl_error=swirl-5/6; pressure_error=pressure-5.

No positivity cap chooses the centered row. Signed source intervals and derivative errors remain. These formulas and their inverse are original source identities, not data fits.

## Reference and exact Rz amplitude

For exact rational phase q in[0,1], g=q*gap. Each incoming reference multiplier exp(-rate*g) is a nonzero formal positive function, retaining sign and source coefficients. With E=V_Rsh-4Z,

    dmean=E+exp(-g)*(dmean_Rsh-E)
    dangular=exp(-1.6g)*dangular_Rsh
    dmixed=exp(-1.6g)*dmixed_Rsh+(5/8)E*(1-exp(-1.6g))
    daxial=E^2+exp(-g)*(daxial_Rsh-E^2)
    dswirl=exp(-1.2g)*dswirl_Rsh
    dpressure=exp(-.2g)*dpressure_Rsh.

V is unchanged throughout reference transport. Its exact log angular amplitude is

    log(Utheta/Pstar)=-log(1+Z^2)
       +(T/10-logCstar-logPstar)*(1-q)-.8q.

At q=1 the large common T/logC/logP terms cancel analytically before enclosure, giving -log(1+Z^2)-.8. Exact reference-offset physical factors collect powers of Pstar^2 before multiplication; no huge physical amplitude or radius is materialized. The phase0 shapes and V remain the exact incoming objects.

## Full original axial restoration and postrestore

For0<=t<=1, V=4Z+E*(1-sigma(t)) and log(Utheta/Pstar)=-log(1+Z^2)+(-8+t)/10. The complete finite kernels are

    J(k,j;t)=integral_0^t exp(-k*(t-s))*(1-sigma(s))^j ds
    (k,j)=(1,1),(1.6,1),(1,2).

Interior prefixes reuse original128 positive cells, with full uncertain endpoint-segment enclosures. At t=1 the same original raw scalar integrals from five_defect_admission are multiplied by exp(-k). The raw integral is integral_0^1 exp(k*s)*(1-sigma(s))^j ds; omitting that exp(-k) factor would change the source. The provider requires ordered finite strictly positive admitted boxes; zero, negative, sign-crossing and nonfinite boxes reject.

    dmean(t)=exp(-t)*dmean_Rz+E*J(1,1;t)
    dmixed(t)=exp(-1.6t)*dmixed_Rz+E*J(1.6,1;t)
    daxial(t)=exp(-t)*daxial_Rz+E^2*J(1,2;t).

Angular, swirl and pressure centered rows retain their original homogeneous rates1.6,1.2,.2. Postrestore offsets[-7,-6] multiply every actual restore-exit row by exp(-rate*(offset+7)); V=4Z exactly. This unpatched path stops at Rm and does not enter the unresolved active patch. Rz/restoration and restore/postrestore centered inlet objects are shared exactly.

The old patch receipt is used solely to bind the full original cutoff integrals and their source family. Its legacy patch-closure flags do not promote this new actual route to active-patch or finite-N closure.

## Physical mixed derivative source rows

The six own shape equations are differentiated through y order4 using the original E*alpha source, alpha=1-sigma and its ordinary y derivatives. Q uses actual V and mean with the original delta denominator. Only Q orders0..4 are used; total mixed order never exceeds4. All ordinary Z Taylor coefficients are converted to derivatives by n! only at the final grid.

Radial prefactors are differentiated before assembly. The binomial rates are .5 for Ur, .1 for Utheta,1.6 for Mtheta/Mtheta_z,1 for Mz and the axial contribution,1.2 for the swirl contribution, and .2 for pressure increment. Pstar^2*P0 appears only in y order0. Physical Mztheta retains its separate R*axial and -R*Utheta^2*swirl/2 terms. No positive exponential cap defines an amplitude. Original mixed-source Q/physical/primitive assignments are AST-bound and independent own-rate/prefactor checks validate the new factored assembly.

These are actual point-function/mixed-derivative packets for source frames0,.5. A uniform radial cell oracle and whole-axis source cells are still required for the complete finite-N five-density integrals and functional terminal closure. No Cartesian NS residual or admissible stress conclusion is inferred from these rows alone.

## Evidence and scope

Files: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_reference_restore_functions.py, .json, _check.py, _check.json. Producer 12.266s; checker 61.047s; terminal exit0. **352 staged exact dependency hashes PASS**. 55 symbolic reference ODE/inlet/centering/Rz geometry, prefactor Leibniz and ordinary-coefficient identities pass.

Two independent finite diagnostic source functions use signed derivative/incoming rows, true sigma, direct own-rate Volterra integrals and complete nonlinear V^2. 432 complete history, 2160 own-shape y-derivative, 1620 physical mixed4 and 12 full-cutoff kernel comparisons pass. Diagnostic finite geometry does not select native parameters. Independent numerical quadrature is reference evidence, not a rigorous interval certificate for the reference solver; production uses directed source/cell/admitted-integral enclosures.

Native checks preserve 14 exact source-object joins, 12 nonzero reference decay functions, 288 actual history Taylor records, 1080 factored physical mixed4 records and 4 P0 object memory checks. 19 invalid phase/offset/source/chart/precision/integral requests reject. Existing **GPT-5.6 Luna / max** read-only review confirmed transport/rate/normalization algebra; its provider-guard and mixed-source binding concerns were resolved by root. Root owns edits/compute/Git.

The new gate is actual conditional two-frame Rsh->Rm source functions. Broader R110 source-history closure, whole-Z oracle, active patch, finite-N/all24 integrals, five controls, global N, completed heat/stress/flat, actual n-recursion and corrected UVW remain false. Upstream core/micro cached inlets remain conditional. Full objective stays active.

## Detailed next production tasks

- [x] **ACTUAL Rsh -> Rz.** Actual centered incoming histories/E, signed nonzero long tails, original exact gap and Rz amplitude identity, derivative/physical units.
- [x] **ACTUAL Rz -> RESTORE EXIT -> Rm.** Full original sigma kernels and nonzero actual incoming histories, exact V restoration, P0/source/T objects, unpatched Rm limit and physical mixed4.
- [ ] **NEXT ACTUAL Rm TARGET / ACTIVE-PATCH INLET.** Translate the new centered Rm rows into the original normalized five-defect namespace with Am=exp(-.6)Pstar/(1+Z^2), exact Rm, mixed/axial units and separate P0. Preserve actual tails and compare defining source graphs, not legacy enclosure overlap. Bind the actual implicit patch target and outer feedback without selecting broad coefficient values.
- [ ] **NEXT SOURCE RADIAL CELLS / FINITE-N COMPILER.** Enclose the new reference/restoration/postrestore functions and mixed derivatives on whole radial cells, with exact source amplitude/radius units and true sigma prefix changes. Bind accepted original stress-direction/raw/inertial and general fixed-phase Z backends. Retain exact all-N phase, both N-dependent levels and denominator hypotheses. Integrate the missing Rsh->Rm five-density windows; the closed Rh_reference integrals already accepted are a different source window.
- [ ] **WHOLE-Z ACTUAL SOURCE PROVIDER.** Replace the two cached source frames with source-correlated axial cells from original analytic core/fixed-point functions, anchored pole primitive and fourteen pressure atoms. Preserve ordinary higher derivatives, source errors, sign/midplane splits and uniform operator bounds. Point Taylor rows cannot certify whole-axis functional identities.
- [ ] **CONTINUOUS CORE/MICRO/R100/RESHAPE PROVIDERS.** Recover core and both original microscopic smoothing charts as functions/cells with dy=hb*ds, true displacement and comparison/actual histories. Complete mixed derivatives and all0<=y<=T long-reshape kernels; cached exits/endpoints remain conditional.
- [ ] **UNIQUE ACTIVE PATCH.** Solve actual implicit repair coefficient functions with Jacobian, uniqueness, actual derivative/error bounds and outer feedback; feed raw_patch_rows and the original finite-N compiler. Keep admitted bounds distinct from solved function values.
- [ ] **ALL24 / Rc TARGETS.** Complete active[1,71/40] plus accepted terminal patch integrals, actual reference/O2/axial/buffer/O3 inlets, original Rc_E/Rc_E_Z and both N-dependent levels. Preserve incoming histories through all downstream paths.
- [ ] **FIVE TERMINAL FUNCTIONS / ONE GLOBAL N.** Close five moment identities as functions of Z, reconcile source/repair/cone/interface and centered complete-target N^-2/Picard/limit/tail inequalities, choose one compatible finite N. Sample closure and scaled incompatible N are insufficient.
- [ ] **MATCHED BACKGROUND / HEAT / ENERGY / STRESS / FLAT.** Complete annular/flatten/analytic-pressure joins, exact heat exterior, axis regularity, finite-energy tail and regional admissible stress/flat remainder.
- [ ] **REAL n-RECURSION / MEANS / TWO PULSE FAMILIES / CORRECTED UVW.** Implement actual n-dependent recovery/repairs and smooth summation with divergence-preserving cutoffs; averaged quadratic stress cancellation and independent corrected Cartesian NS. Export physical u/v/w and measure contraction, relative elongation, recursive scale relation and material winding separately.

Mark DONE with defining code, scoped report/receipt and a commit. Continue production after necessary changed-scope checks pass; retain the full objective until all required layers are achieved.
