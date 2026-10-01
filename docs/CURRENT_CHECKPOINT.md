## Pressure-to-finite-core propagation (2026-09-30)

All14 pressure approximation derivatives bounded uniformly through order24. Actual degree18 finite core at Z=.3 now has pressure-only coefficient and five core-moment error intervals. Uz error at R=4/Lambda <=2.253655468045584e-36. This excludes axis errors, radial truncation/RK continuation and full five-defect closure. See PRESSURE_HIGH_DERIVATIVE_CORE_PROPAGATION_2026_09_30.md.

## Uniform pressure approximation budget (2026-09-30)

All 14 stored-datum pressure stages now have an absolute error enclosure on |Z|<=0.8: weighted normalized C2 upper 3.079264485429091e-12. Flatten uses a conservative positive-mass bound; relative flatten error, original parameter errors, runtime evaluation roundoff, core/RK propagation and five-moment closure remain open. See POINTWISE_PREHEAT_ERROR_BUDGET_2026_09_30.md. Next: relative flatten coefficient error and pressure-to-core/moment propagation.

## Latest completed pressure error budget (2026-09-30)

All 14 pressure components now have a same-source, stored-parameter pointwise error budget at Z=0.3. Weighted normalized error upper: 5.6517545845415295e-13. This is not a uniform C2 bound or five-moment closure. See [the detailed receipt](POINTWISE_PREHEAT_ERROR_BUDGET_2026_09_30.md).

Next: uniform |Z|<=0.8 pressure approximation errors, then core/RK and five-defect propagation; retain the common analytic datum. Temporal recursion and oscillatory correction remain open.

# 2026-09-30: Both steep pressure transitions enclosed at high order

Shared order12 interval integration now handles the two post-flatten steep
transitions, including negative primitive coefficients and stored switch
origin/length offsets. Actual64-panel relative mass error upper bounds:
1.69391e-12 and1.09072e-12; both widths shrink32->64. Axial derivative
errors are zero in these beta0 preheat stages. First-two-stage regression
results remain unchanged. Variable-beta flatten and subsequent propagation
to the common core/five moment defects remain open.
See STEEP_PREHEAT_INTEGRALS_2026_09_30.md.

# 2026-09-30: High-order pressure integral error enclosure

Directed interval Taylor algebra and midpoint remainder integration now
replace low-order subdivision for the first two variable preheat stages.
Order12 with64 radial panels gives normalized mass error upper bounds
2.21850e-13 and2.97946e-41 against the retained192-point Gauss values.
Widths shrink32->64. Analytic flat-edge bounds preserve nonzero startup
contributions without evaluating singular derivative formulas at endpoints.
Independent Taylor and switch/primitive checks pass. Remaining variable
stages and propagation to the common source/five defects are still open.
See HIGH_ORDER_PREHEAT_INTEGRALS_2026_09_30.md.

# 2026-09-30: Convergent variable-stage integral/error intervals

Directed monotone switch and primitive grids now feed positive radial
interval subdivision. Actual first two transition widths shrink64->128
panels. Signed true-minus-finite Gauss mass intervals give absolute error
upper bounds0.00240550 and1.01380e-8; beta2 converts them to compact
axial derivative budgets. Six negative-slope stages have analytic integral
intervals and finite-mass error receipts. The first transition bound is
still too loose for the target five-defect hierarchy. Variable-beta
flatten errors, other transitions and source propagation remain open.
See VARIABLE_PREHEAT_INTERVAL_INTEGRALS_2026_09_30.md.

# 2026-09-30: Complete stored-parameter preheat norm bound

H=1 collar and exact exterior integral enclosures now use the shared stored
heat normalization; no pressure datum is reset. All eleven finite stages
plus reference and terminal pieces give a complete conservative normalized
C2 pressure norm upper bound47.82861181496495 on |Z|<=.8.
Independent resolved collar/exterior and normalization-shift checks pass.
JSON endpoints retain exact binary MP tuples as well as rounded displays.
This bounds the stored-parameter analytic pressure integral, not its Gauss
approximation error or the five centered defect norms. Core/RK propagation,
original parameter derivations and actual heat/energy gates remain open.
See COMPLETE_PREHEAT_ENCLOSURE_2026_09_30.md.

# 2026-09-30: Stored-schedule endpoint and radial-stage enclosures

Directed interval evaluation removes primitive quadrature from preheat
log-amplitude endpoint bounds. Exact switch symmetry gives J(1)=1/2.
Actual 11 finite pre-collar stages have positive mass/axial derivative
upper bounds; maximum endpoint log width2.39823e-205. Negative-slope
integral caps avoid huge pulse-length estimates. Resolved primitive and
independent positive-stage integration checks pass. Scope is relative to
stored Decimal parameters, not exact exp/log paper-input derivations.
Heat collar, tight quadrature remainders and core/RK errors remain open.
See SCHEDULE_ENDPOINT_ENCLOSURES_2026_09_30.md.

# 2026-09-30: Actual common preheat compact-interval envelope

Actual common schedule: 205 positive pressure atoms at radial Gauss order192,
260 digits, retained separately. Analytic axial factor bounds on |Z|<=.8
through derivative3 give normalized C2 envelope46.00696349. Radial
quadrature and scalar rounding errors remain unenclosed. This is a bound
for the declared finite atom representation, not full-source closure.
Positive radial-stage mass/variable-beta derivative bound APIs added;
their schedule endpoint inputs still need enclosure.
See PREHEAT_INTERVAL_BOUNDS_2026_09_30.md.

# 2026-09-30: Conditional interval C2 inverse majorant

Implemented factorial-weighted C2 Banach norm, exact common amplitude
inverse-square bound on compact |Z|<=a<1, and finite Picard tail estimate.
Resolved analytic family checks pass; original C1 smoke remains passing.
This does not enclose the actual 69-part source norm or quadrature constants.
See FIVE_BUMP_C2_MAJORANT_2026_09_30.md.

# 2026-09-30: Second-Z centered source and corrected field completed

All 69 labeled centered defect parts, flat kernels through derivative 13,
finite nonlinear inverse and 71-part physical moment reconstruction now
carry value/first/second axial derivatives. Actual Z=.3, x=1.25 completed;
P0/P0_ZZ preserved and Ur agrees with the shared stress recovery.
Maximum nominal second-Z inverse residual after ten updates: 5.85616e-209.
Resolved density second derivative discrepancy: 1.75147e-14; corrected
field scaled discrepancy: 9.64652e-32. These are finite construction checks.
Uniform axial bounds, source error enclosures, relative-flat/infinite inverse
closure, heat/energy matching and temporal recursion remain open.
See SECOND_CENTERED_DEFECTS_2026_09_30.md.

# 2026-09-30: Second axial source adapters and actual inner endpoint

Second-Z arithmetic, comparison, bridge, exponential continuation and R100--R110
switch adapters implemented. Resolved independent shifted-source/RK stencils
pass with max scaled error 1.513e-32, including five moments and Ur_Z.
Fresh actual degree18, Z-depth3 source completed at Z=.3; inner comparison
endpoint receipt saved. Actual collar/R100 propagation completed with second-Z fields/moments and
raw energies; an independent cache and post_collar_R100 receipt were saved.
Do not use the old first-Z cache as second-Z source data.
See SECOND_AXIAL_PROPAGATION_2026_09_30.md.
Actual R110 second-Z execution completed with inlet continuity checked.
Interval/source error
bounds and all previously open global matching/recursion gates remain open.

# 2026-09-30: Incremental analytic five-bump inverse installed

The callable same-source correction now supports incremental Picard inversion,
with individual updates and first-Z tangents preserved. Actual Z=.3, x=1.25:
10 nonlinear updates reduce maximum nominal centered terminal residual to
9.44313e-211. Same-source P0/P0_Z retained. This is the analytic local five-bump
inverse, not temporal n-recursion or uniform functional closure.
See FIVE_BUMP_INCREMENTAL_INVERSE_2026_09_30.md and inverse_check.json.
Actual uniform C1/C2 defects, infinite response/source error enclosures and
heat/energy matching remain open. Earlier cone samples concern degree-three
fields and are not transferred as a certificate of this new inverse field.

# 2026-09-30: Actual-source relaxed-cone sample and terminal tail

Actual same-source five-bump reconstruction at Z=.3 now has 38 interval samples.
All pass the kappa<=2 relaxed branch; paper claims only relaxed cone in this
region, not admissible cone. Full degree-three velocity yields nonzero degree
4--6 centered moment tails: row2 +8.93e-85, row4 -8.32e-71, row5 +5.56e-71.
These are scalar nominal input-composition diagnostics, not exact functional
closure or uniform-Z bounds. See FIVE_BUMP_CONE_AND_TERMINAL_TAIL_2026_09_30.md.
Uniform C1/C2 source bounds, response remainder, heat/energy matching and temporal
n-recursion remain open. Existing tiny source and response terms remain separate.

## 2026-09-30: actual common-source reference/bump join executed

ReferenceDefectBackground now reconstructs the terminal reference interval from the exact power moments plus all centered source defect parts, keeping the original P0. It includes reference_power in the authoritative moment part maps and does not infer separate raw energies from d4. Resolved independent radial-density and source-lineage checks pass (maximum relative error3.20e-30).

The actual warm-R100 run at Z=.3, x=1.25, P9/W2 completed: same P0/P0_Z, shared stress recovery, 69 source labels plus reference_power and five_bump_correction. A local ignored reference snapshot preserves the exact evaluated source for further diagnostics without repeated core/defect work. See REFERENCE_DEFECT_BACKGROUND_2026_09_30.md and five_moment_reference_background_check.json.

The corrected raw swirl convention is integral(u_theta^2/2)dR, matching the old common provider; Mztheta=raw_axial-raw_swirl. Joined part maps now include correction increments, including optional unavailable raw state handling. This is a finite pointwise join, NOT uniform five-moment closure, cone certification, full outer/heat-energy completion or temporal recursion.

## 2026-09-30: finite bump field and full quadratic moment recovery

The correction field adapter recovers u_theta, Uz, F, pressure/moments, raw axial/swirl quadratic changes and value-level Ur from the same first-Z data, preserves P0, enforces fixed Rm and updates corrected stress/shear inputs. The resolved join and independent density/first-Z radial checks pass. Ur_Z needs second-Z data and remains unavailable.

The actual serialized Z=.3 P9/W2 correction at x=1.25 retains 55 velocity monomials and 461 cumulative-moment monomials, including 406 nonzero terms above defect degree3. All quadratic terms through degree6 are retained, so finite response truncation is not mistaken for closure. The actual full baseline join remains open. See FIVE_BUMP_FIELD_2026_09_30.md.

Conditional C1 contraction/Catalan tail bounds are implemented for supplied uniform e and the finite fixed weights. Actual uniform e and analytic quadrature enclosures remain unverified. Five-moment functional closure, uniform C2/cone, energy/heat and genuine temporal recursion remain open.

## 2026-09-30: paper five-bump map and separated defect response

Implemented Section10.2 Eq.(10.8): two axial and three angular compact bumps, the fixed invertible matrix, explicit quadratic coupling, first-Z support, Jacobian and partial cumulative moment changes. Independent field-density/inverse/tangent replay at order96 has relative errors below4.41e-16; quadrature is unenclosed.

The formal five-defect inverse response now retains55 monomials through degree3. Independent coefficientwise residual~2.01e-87, scalar contraction coefficient comparison~1.67e-22 and first-Z check~1.48e-52. The actual serialized Z=.3 P9/W2 run retains55 monomial responses, separately nonzero d3/d5 contributions, and all69 original source labels in the linear response. Nonlinear source-stage products are not individually expanded. See FIVE_BUMP_MAP_2026_09_30.md and FIVE_BUMP_RESPONSE_2026_09_30.md.

These are finite map/response milestones, NOT actual functional five-moment closure or temporal recursion. Convergence/remainder control, uniform analytic input bounds, corrected-field installation, radial energy/heat and stress-cone requirements remain open.

## 2026-09-30: all five centered finite-component defect inputs computed

The common actual R110 source now feeds CenteredComponentDefects with pressure order 9, width order 2, and automatic first-Z tangents. The Z=.3 receipt preserves all five nonzero baseline defects and 69 separate source labels. d1~2.23437e-15, d2~5.68986e-16, d4~5.91501e-41; negative d3/d5 retain arbitrary-exponent values (log magnitudes ~-1.6e152/-2e151). Explicit P0 is unchanged. See CENTERED_COMPONENT_DEFECTS_2026_09_30.md and the centered_component_defects_check.json receipt.

This completes finite-component defect assembly, NOT five-moment closure. Uniform analytic smallness, higher axial derivatives, source/jet/quadrature bounds, actual five-bump repair, exact heat/energy/cone and temporal recursion remain open.

## 2026-09-30: flat reshape defects lifted to pressure9/width2 components

The scalar kernel now supports positive integer B derivatives with independent saddles. compose_flat_shape_defect retains the rectangular finite-ring value and first-Z tangent as separate derivative-order terms. Actual serialized common-source pressure9/width2 B/B_Z replay passes all three angular kernels: all30 atoms survive and the n11 (9,2) contribution remains separately nonzero; the tangent uses derivatives through12. Independent resolved P2/W1 coefficient integration discrepancies are3.73e-60 value and3.39e-61 tangent; scalar derivatives1..3 differ by3.83e-59. All source/finite-ring/window/quadrature errors remain unenclosed. Next merge these components with centered axial/inner sources into all five defect functions, then Section10 repair. Details: [COMPONENT_FLAT_DEFECT_2026_09_30.md](COMPONENT_FLAT_DEFECT_2026_09_30.md).

## 2026-09-30: nonzero flat reshape sources recovered in signed-log coordinates

New flat_shape_defect kernel integrates exp(-k ell) expm1(m B sigma(ell/T)) and its first B sensitivity using saddle-centered composite quadrature. Actual leading B/B_Z atoms from the committed R110 receipt retain nonzero angular, pressure and angular-energy defect sources at T=4e152. Their modes are around1e102, far beyond the old1842 cutoff. The independent leading saddle log formula differs by about3e-51 to6e-51. This is a scalar leading-atom source receipt, not the complete five functional defects; finite windows/quadrature/source errors are unenclosed. Details: [FLAT_SHAPE_DEFECT_2026_09_30.md](FLAT_SHAPE_DEFECT_2026_09_30.md). Next lift the kernel into pressure/width components and merge centered axial sources before Section10 repair.

## 2026-09-30: component axial restoration reaches Rh

PressureWidthAxialRestore implements Section9.38 from the common Rz entry through axial restoration and reference continuation to Rh (phases0..3), retaining prior moment/raw-quadratic parts and adding separate restoration/reference increments. Actual degree18 pressure9/width2 Z=.3 replay passes phases.5,1,3: entry field discrepancies are zero; Uz=4Z and Uz_Z=4 at phases1,3; axis P0 and prior parts remain unchanged. Independent resolved physical integral discrepancy5.99e-11, first-Z stencil4.16e-9 and receipt2.54e-9 remain unenclosed. The reference adapter now also exposes already-computed raw quadratic totals required by this chain. Next construct cancellation-resistant functional defects for Section10, not total-minus-target rounded zeros. Global matching, energy, cone, heat and temporal recursion remain open. Details: [COMPONENT_AXIAL_RESTORE_2026_09_30.md](COMPONENT_AXIAL_RESTORE_2026_09_30.md).

## 2026-09-30: component reference extension reaches axial restoration entry

PressureWidthReferenceExtension continues the common reshape endpoint to Rz=exp(-8)Rref with analytic five-moment/raw-quadratic increments and inherited first Z tangents. Original seeds, reshape and reference increments remain separate; axis P0 is unchanged. Actual degree18 pressure9/width2 Z=.3 replay reaches Rz and preserves nonzero Uz_Z pressure-tail atoms; recorded boundary field discrepancies are zero. Independent resolved integral error is 3.62e-100 and first-Z stencil discrepancy 1.02e-15. Global matching, terminal closure, energy, heat, cone and temporal recursion remain open. Details: [COMPONENT_REFERENCE_EXTENSION_2026_09_30.md](COMPONENT_REFERENCE_EXTENSION_2026_09_30.md).

## 2026-09-30: component long reshape reaches reference endpoint

PressureWidthLongReshape carries the common R110 state through Section9.30 with separate inherited five-moment/raw-quadratic seeds, increments and first Z tangents. Actual pressure9/width2 degree18 source at Z=.3 completes phases .5 and1 with T=4e152; phase-zero fields agree exactly at recorded precision, nonzero Uz_Z pressure-tail survives, and leading endpoint reference log/slope identities agree to arithmetic precision. Independent resolved T400 normalized quadrature differs by 3.13e-26 and first-Z stencil by 6.11e-17. Source/quadrature/jet errors remain unenclosed; global installation, reference extension, functional terminal moments, energy, cone and temporal recursion remain open. Details: [COMPONENT_LONG_RESHAPE_2026_09_30.md](COMPONENT_LONG_RESHAPE_2026_09_30.md).

## 2026-09-30: component shear switches reach R=110

PressureWidthExitSwitches now propagates the actual derivative-aware R100 state through both Section9.4 switches and the analytic a=4/5 constant-power segment. Actual pressure9/width2 degree18 replay completes: boundary fields/moments agree below 1e-200 scaled error, leading F110/F100=1.1^(-.4), Uz_Z pressure-tail survives, and raw axial-minus-swirl matches its moment. Independent resolved scalar physical-moment replays differ by about 9.61e-13 at hb=1e-6 and 1.02e-9 at hb=1e-5. Local Z=.3 source/reshape inputs fit their nominal bounds; uniform C2/cone, long reshape, functional terminal moments, energy and temporal recursion remain open. Details: [COMPONENT_SHEAR_SWITCHES_2026_09_30.md](COMPONENT_SHEAR_SWITCHES_2026_09_30.md). A local ignored R100 cache allows --resume of bounded switch work without recomputing the collar.

## 2026-09-30: positive-epsilon prescribed exit continued analytically to R=100

PressureWidthExitContinuation now integrates the prescribed field changes from the frozen auxiliary drivers as finite-ring exponential polynomials. It preserves epsilon-driven F/Uz and all five moment changes rather than freezing prescribed values. Actual pressure9/width2 degree18 replay completes R=1,10,100; first-width F/Uz increments and pressure-tail Uz_Z survive. Independent resolved scalar RK4 fields/five moments agree to 2.93e-13 scaled error. Independently differentiated analytic moment primitives give local finite-ring divergence/derivative identities at arithmetic precision. No full Cartesian/global certificate, source/collar error enclosure, switching/reshape, terminal moment closure, energy, cone or temporal recursion is claimed. Details: [ANALYTIC_POST_COLLAR_EXIT_2026_09_30.md](ANALYTIC_POST_COLLAR_EXIT_2026_09_30.md).

## 2026-09-30: automatic axial jets and moment-derived radial velocity at prescribed exit

AxialPressureWidthExitBridge differentiates the same core, auxiliary drivers and prescribed finite RK state. Actual degree18 pressure9/width2 replay at delta=1e-200 and h_b=exp(-100-1e154) retains nonzero Uz_Z/P_Z pressure-tail atoms and all five first-width moment-Z increments. Start-core fields and recovered Ur agree below 1e-200 scaled error. Local Cartesian chart and same-data stress evaluation are available. This is a local finite-jet construction: the ODE divergence identity is not independent Cartesian validation. Switching/reshape, functional terminal moments, finite radial energy, exact heat, cone and temporal recursion remain open. Details: [AXIAL_JOINT_EXIT_2026_09_30.md](AXIAL_JOINT_EXIT_2026_09_30.md).

## 2026-09-30: actual tiny-width atoms enter prescribed base exit ODE

PressureWidthExitBridge now integrates the original Section 9.25-9.26 base velocity/five-moment ODE in s=y/h_b, using the same bivariate auxiliary comparison. Actual degree-18 source replay with pressure order9/width order2 retains nonzero pressure-tail/axial atoms, nonzero width-driven log(F/Fa), and nonzero first-width increments of all five cumulative moments at s=2. Pressure orders3/9 overlap at serialized precision; this is not a remainder enclosure. Independent resolved-family scalar bridge replay differs by at most 4.1148e-18. Actual h_b=exp(-100-1e154) is retained separately rather than lost in each atom. Full prescribed Z tangents, switching/reshape, functional moment repair, cone and finite-energy closure remain open. Details: [JOINT_PRESSURE_WIDTH_EXIT_2026_09_30.md](JOINT_PRESSURE_WIDTH_EXIT_2026_09_30.md). Pressure/width Taylor order here is not temporal coefficient recursion.

## 2026-09-30: revised goal and pressure-component auxiliary exit

PROJECT_GOAL.md now records the user-updated seven-stage paper-faithful objective. ComponentExitComparison propagates complete-preheat core atoms through the original Section 9.23 auxiliary comparison, including five moments, pressure, Z jets, stress and frozen D/E drivers. Actual source replay retains nonzero pressure and axial atoms through five boundary/transition/frozen probes; zero-tail scalar replay differs by at most 4.50e-233. Independent resolved-scale fields/moments/driver replay differs by at most 5.74e-53. Pressure-parameter order 9 and ODE errors remain unenclosed; tiny actual exit-width increments can still be lost within each atom. This is neither the prescribed exit bridge nor annular completion. See [COMPONENT_AUXILIARY_EXIT_2026_09_30.md](COMPONENT_AUXILIARY_EXIT_2026_09_30.md). Next preserve width dependence and propagate the same data through the prescribed bridge/reshape and functional moment repair.

## 2026-09-30: complete preheat tail propagated into local nonlinear core

New component_pressure_core retains Pprefix + lambda*Ppost through the original nonlinear radial recurrence, without pressure-parameter truncation. Actual degree-18 source retains powers through 9; the nonzero post-Rv pressure tail and its first axial radial response survive even though nominal prefix+tail loses them. Local five moments and divergence-derived Ur are available; Lambda R=1,4 component-scaled divergence is at most 5.68e-259. Independent scalar/component recurrence replay differs by at most 2.57e-101. Component-valued signed angular solver also retains nonlinear r^2 and exp(-5e152) pressure atoms; nominal atom residuals do not certify the unresolved scalar/full-field row. Read [COMPONENT_PREHEAT_CORE_2026_09_30.md](COMPONENT_PREHEAT_CORE_2026_09_30.md). Next propagate component core and angular coefficients through exit/reshape/incoming corrections; the existing scalar joined field is still unmodified. Finite energy, cone and temporal recursion remain open.

## 2026-09-30: coherent continuous waiting rebuild removes oversized angular target

Opt-in build_joined_field(continuous_pressure=True, coherent_waiting=True) now computes waiting from the installed continuous angular primitive and retained collar deficit, then reconstructs schedule/pressure/core/inner dependents. Initial oversized target -1077.580145807 was traced in part to a missing (1-epsilon) conversion between preheat and heat normalization; it is not a certified physical defect. The corrected Z=.3 retained angular target is +2.31629269314e-837 (|r|/mu^29=4.56070569695e-45); nominal terminal pressure improves from 2.15130102859e-69 to 1.18757636672e-98. These are nominal/identity-based quantities, not enclosed terminal-moment closure. Signed coupled solve and analytic tangents are implemented; complete analytic preheat core components and tiny pressure-row representation remain unresolved. New report and next actions: [COHERENT_WAITING_COUPLED_ANGULAR_2026_09_30.md](COHERENT_WAITING_COUPLED_ANGULAR_2026_09_30.md). Finite energy, cone and recursion remain open.

## 2026-09-30: complete preheat datum derivation

Lei-Ren v2 equation (6.10) defines P0_pre by replacing H by 1 over the complete future angular profile, not merely cutting at Rv. The flatten interval retains a varying analytic Z exponent; all later preheat intervals are Z-independent. The implementation lane now extracts this datum as separate prefix/post-Rv atoms and arbitrary-center Taylor components. New adapter ContinuousPreheatPressure supplies the complete target and arbitrary-center stage Taylor jets. The actual-source check at Z=.3 retains nonzero post-Rv pressure and derivative with log absolute magnitude around -2.719157344816895e28; nominal total-minus-prefix loses that tail. The independent derivative discrepancy is 2.15e-60 relative. MP128/192 flatten change is 6.36e-39 relative (unenclosed). Receipt: experiments/root_st073/lei_ren_part1_paper_continuous_preheat_pressure_check.json. Do not install a rounded sum as a complete core jet. P1 derivation/adapter is complete; coupled bump solve and component-aware core reconstruction remain P2.

## 2026-09-30: two local papers revise the route; angular MP replay improved

Read Lei-Ren 2609.35406v2 (29 Sep 2026) and Duraiswami 2609.17642v1 from user-supplied local PDFs. Current route and prioritized tasks: [TWO_PAPER_ROUTE_2026_09_30.md](TWO_PAPER_ROUTE_2026_09_30.md). This supersedes earlier instructions to append raw post-Rv heat pressure to the axis datum: v2 Sections 7.4 and 12.4 require simultaneous angular-moment and analytic preheat-pressure restoration, with coherent dependent core reconstruction. The current pre-Rv pressure adapter is not yet proved equal to that complete analytic target.

Continuous angular moments/swirl now use coherent MP nodes; incoming inner offsets regenerate from the live reconstructed inner field rather than old seed receipts. Completed Z=.3 source replay has heat Ttheta/legacy approximately 8.0702e-69; flatten Ttheta/legacy remains approximately 1. P(infinity) stays approximately 2.1513e-69. Comparisons combine node and source changes; absolute stresses remain huge. This supersedes the older statement that angular stress is unchanged everywhere. Five functional moments, finite energy, admissible cone and scale recursion remain unresolved.

New angular tail-constant extractor retains actual anchor, bump, inner offsets and quadratic heat atoms. Independent exponential fixture passes at approximately 4.92e-91 relative error; tiny-delta denominator avoids cancellation. Exact heat remainder and interval bounds are not enclosed. No actual moment or field is overwritten.

## Rebuilt-source coupled replay completed

The public build_joined_field(continuous_pressure=True) path reproduces the nominal terminal pressure improvement and completes coupled actual five-moment stress evaluations at Rv+50 and Rtail+4. Relative total axial stress is 9.8244e-97 and 1.7797e-67 of the legacy probe, respectively. Angular stress remains unchanged at serialized precision. Comparisons include old binary64-96 versus new MP-192 pressure quadrature and the reconstructed inner source; neither candidate has a certified stress cone. Prioritize coherent angular/mixed-moment restoration next, retaining unresolved inner inputs and terminal radial energy obstruction.

## 2026-09-30: coherent pressure-datum core rebuild reduces nominal defect by 67 orders

MP pressure-stage orders 128 and 192 converge near the old P(infinity)=-0.01229506638282; the difference is 1.56e-15. Raising those orders cannot explain the 0.0123 mismatch. The core was constructed using a legacy pressure anchor before the outer continuous angular schedule was installed.

New opt-in path: build_joined_field(continuous_pressure=True), followed by ContinuousIncomingOuterField on that source. ContinuousAxisPressureJets integrates the same pre-Rv angular schedule with MP nodes and exact constant-stage atoms before computing the nonlinear core pressure Taylor coefficients. Exit, reshape, axial restoration, inner moment corrections and exterior axial inputs are reconstructed. Shared pressure queries select coherent MP quadrature automatically. No completed field's Z-dependent pressure is shifted.

Actual Z=.3 replay gives P(infinity)=2.15130102858891e-69, P(infinity)_Z=-2.36840480211623e-69 and pressure transport coefficient 3.44602898707912e-69. Relative to the old nominal pressure defect this is a factor 1.7796940893e-67. This is improvement, not exact closure: inner uncertainty entries 3 and 5 persist, the post-Rv axis pressure tail jet remains omitted, quadrature is unenclosed, and finite energy/stress cone/recursion remain unproved.

The stress evaluator/bundle now expose separate transport, linear, mixed/quadratic and pressure components plus the coupled angular Btheta target. The legacy two-point replay sums those pieces back to the original stresses; no moment is overwritten. New candidate receipt: experiments/root_st073/lei_ren_part1_paper_continuous_pressure_rebuild.{py,json,md}.

## 2026-09-30: infinite pressure target and node sensitivity measured

The pressure provider now exposes terminal_pressure_jet(Z), retaining actual inner anchor, preheat, compact bump and complete infinite heat contributions. The restoration target is extracted without changing the pressure datum or coefficients. New reproducible receipt: experiments/root_st073/lei_ren_part1_paper_continuous_pressure_target.{py,json,md}.

At Z=.3 nominal P(infinity)=-0.01208803828438 and P(infinity)_Z=+0.01330793205620. Replacing only pressure variable-stage binary64 Gauss nodes with same-order MP nodes changes these to -0.01229506811512 and +0.01353585480563 (about 1.7 percent shift). Totals before cancellation are about 4e12. Therefore materialized pressure mismatch is real in the current numerical evaluation, but its interpretation as an analytic target defect remains unproven. Resolve MP quadrature convergence and inherited inner pressure input before solving a late compact bump against this target. No pressure gauge reset, coefficient change, finite-energy or stress-cone certificate is claimed.

## 2026-09-30: actual five-moment bundle installed; pressure mismatch isolated

Shared profile pressure and partial quadratic APIs are installed. Nominal bounded partial-energy and pressure receipts are available. Integrated actual five-moment/pressure/velocity stress replay at Z=.3 (post-Rv flattening and heat) completes but yields enormous stress/shear ratios; coupled matching is not closed.

The post-support pressure combination H=2(1+delta) Z P-(1-Z^2) PZ is approximately -0.01936304114, with P=-0.01208803828 and PZ=0.01330793206 at both points. N_z equals the pressure contribution to I_z when axial velocity/jets vanish; the saved I_z and N_z agree to their serialized precision. This identifies pressure compatibility as an immediate actionable obstruction, while smaller moment terms remain unresolved. Do not remove it with a Z-dependent pressure gauge reset.

Artifacts: experiments/root_st073/lei_ren_part1_paper_continuous_{partial_axial_energy,pressure_moments,moment_bundle,bundle_stress_check,pressure_tail_defect}.{py,md,json} as applicable. Finite energy, admissible stress, full residual and scale recursion remain uncertified. The previously materialized nonzero radial-energy tail remains open.

## 2026-09-30: quantitative terminal radial-energy obstruction audit

Added continuous_radial_energy_tail.py/md/json for the current regenerated candidate. It retains Mz/Mz_Z and C=[(1-delta)Z Mz+(1-Z^2)Mz_Z]/(1-delta Z^2), giving physical ur=-nu*C/r. At fixed tau, radial kinetic-energy density per dZ dlogR is pi*nu^2/2*(dz/dZ)*C^2 with the explicit physical axial Jacobian.

At Z=.3, nu=1 and tau=exp(-4), materialized C remains nonzero and the pointwise logarithmic energy density is positive. Post-Rv and heat mass/transport agree nominally to about 1.35e-292; density to about 2.70e-292. These are arithmetic consistency results, not accuracy certificates for the residual. No Z interval is integrated, no terminal value is forced zero, and the exact functional coefficient/integral closure remains unproved.

This quantitative audit keeps the full finite-energy requirement visible while partial-energy and pressure worker checks continue. Next accept their separate-component results, run unified moment/stress diagnostics, and enclose actual integral inputs/terminal transport. The current numerical radial tail must be resolved before any global finite-energy or recursion claim.

## 2026-09-30: analytic radial velocity jets prepared for stress

Added velocity_radial_jets(profile,logR,Z) using the same installed angular slope, incoming cutoff derivative, live pulse product jet and owned end-bump derivative. It supplies Utheta_y and Uz_y in logR coordinates alongside the existing Z jets, without changing coefficients or terminal moments.

Independent shared-field point derivative checks passed in incoming, pulse, end-bump and heat regions; maximum local relative difference is 1.61e-16. These are nominal functional checks, not input uncertainty or stress-cone certification. The helper is available for the next same-candidate stress adapter after actual five-moment/P providers are accepted.

Partial axial energy and preheat pressure workers remain active. Their Python processes were confirmed live during this turn. No incomplete worker dependency or prepared local quadratic dispatch is published as completed work.

## 2026-09-30: shared pressure heat cumulative and infinite-tail atoms

Added pressure_increments(t,Z) and complete_pressure_heat_integral(Z) to the common heat-moment owner. They use the pressure weight Utheta^2/(2R), exponent -(1+delta), and the same retained heat polynomial/collar atoms as the installed angular field. Reference, tiny correction and correction Z remain separate. Infinite-tail analytic polynomial truncation is bounded; quadrature/arithmetic are not.

Shared-field finite pressure-integrand checks at t=.6 and 4 passed with maximum relative difference 9.50e-14. Infinite-tail Z difference passes nominally at working precision, and Rtail increments vanish. Actual preheat/inner pressure provider remains in progress; full pressure matching, five moments, finite energy and recursion are not claimed.

The partial-axial worker's Python PID44536 was confirmed live during this turn. Parent has prepared lazy quadratic dispatch locally but will publish it only with the accepted provider, keeping uploaded code free of unfinished dependencies.

## 2026-09-30: incoming swirl reference regenerated coherently

Replaced inherited incoming I_swirl by the same normalized ContinuousAngularMoments propagation used by actual angular cumulative moments. Reference energy is generated before actual inner offsets are applied once; energy target and live axial coefficients are then recomputed together. Reference swirl/Ep^2 remains Z independent, preserving the analytic tangent cancellation for this input.

Shared-field test passed: actual Rp swirl and Z primitive match regenerated reference plus measured inner offsets to about 1.17e-84 nominal relative difference. The change from inherited I_swirl is about 7.92e-19; this is not a global accuracy certificate. Updated linear-equation replay errors are about 6.87e-174 and 1.09e-173, and the live coefficient tangent evaluates. Quadrature/input uncertainty, terminal radial transport and complete five moments remain open.

Pressure cumulative implementation and partial axial-square integration are assigned to bounded Luna/max workers. Their unfinished files are not accepted or uploaded as completed work. Next integrate their verified results, complete full pressure/heat targets and stress diagnostics.

## 2026-09-30: cumulative heat angular moments and infinite swirl component

Extended actual seeded angular and swirl-square moments/Z jets beyond Rtail using the installed common heat polynomial. Collar Gauss atoms and analytic exterior exponential atoms retain reference, heat corrections and Z terms separately. Actual quadratic moments now dispatch through finite heat radii. Six forward-integrand checks at t=.6,4,30 pass with maximum relative difference 7.69e-14; Rtail increments are zero and corrections remain nonzero.

Added complete_swirl_heat_integral for Rtail..infinity Utheta^2 dR, retaining its reference/correction/Z atoms and an analytic integrated polynomial-truncation bound. Positive delta makes this swirl radial integral converge. This does not include physical Z weights or radial kinetic energy; total finite energy and recursive closure are still unproven. Quadrature/arithmetic and full heat-target uncertainty remain open.

Next regenerate exact heat-defect/angular/pressure targets from these shared atoms, regenerate inherited incoming swirl coherently, complete pressure moments/Z, finish partial axial energy and resolve terminal radial transport.

## 2026-09-30: continuous heat point jets installed

Installed the new common heat functional in ContinuousAngularSchedule and propagated separate heat deficit/log-amplitude/slope correction atoms to actual angular field adapters. Preserved nonzero heat Z derivatives beyond Decimal exponent storage using arbitrary-exponent MP values. Shared-field collar and exterior checks passed: retained-log radial derivative maximum relative difference 1.83e-12; Z derivative and actual field jet nominal agreement at working precision. The live nominal future-energy target is regenerated on construction; exact heat-defect integral targets still retain their explicitly declared Taylor approximations. Complete heat cumulative moments, input bounds, finite energy and recursion remain open.

Evidence: continuous_heat_install_check.py/md/json. Next integrate the retained heat atoms into full cumulative angular/quadratic/pressure moments, regenerate correction targets coherently, and complete partial axial energy.

## 2026-09-30: separate heat deficit/jet provider prepared

Added continuous_heat_kernel.py/md/json: one Gamma integral for H/H-prime/H-second, normalized positive deficit quadrature, separate logH and tiny derivative correction, plus a common small-x quadratic Taylor approximation with explicit analytic truncation bounds. A first run exposed absolute-tolerance loss in the tiny deficit integral; normalization by h fixed it. Functional derivative checks now pass with maximum relative difference 4.28e-19; exp(-1e6) argument retains nonzero deficit. Arithmetic/quadrature are not enclosed. Provider is not yet installed in the shared field or targets.

Next install using separate heat atoms, regenerate angular/pressure/energy targets consistently, and continue cumulative heat moments. Never interpret this standalone provider as full heat or finite-energy closure.

## 2026-09-30: mixed moments and complete exterior axial energy installed

Added actual inner-seeded mixed moments/Z jets and complete axial-square/Z jets, using the live continuous incoming, pulse and end atoms. Installed linear_axial_moments_jet and quadratic_moments_jet on ContinuousIncomingProfile. The quadratic provider currently supports Rv through Rtail. Actual inner offsets are applied once; materialized terminal residuals remain nonzero. Pulse/end supports are disjoint by xi = 13 + mu*t_v > 11 on the end band.

Fixed the tiny positive pulse startup primitive using a positive scaled endpoint integral. The startup checks retain a nonzero value at xi=1e-5; full source energy integration runs again. Mixed Rp matching is about 2.28e-81; local mixed integrands about 2.08e-23 and 1.19e-17. Exterior quadratic integrand checks are about 7.34e-15 and 1.88e-14; axial energy is constant after axial support. These are nominal consistency checks, not certified quadrature accuracy or NS momentum residual bounds.

Next: partial axial-square/Z primitives inside incoming/pulse/end; consistent heat H and derivative primitives; complete pressure/Z; regenerate inherited incoming swirl from the installed angular primitive and update targets coherently. Enclose inputs and coefficient tangents, resolve nonzero terminal radial transport, then finite energy, stress/remainder and scale recursion. No full closure claim.

## 2026-09-30: actual inner-seeded angular moments extended to Rtail

Added continuous_angular_moments.py/md/json and installed ContinuousIncomingProfile.angular_moments_jet. It transports actual Mtheta, full swirl-square integral and both Z derivatives from Rh through every preheat stage to Rtail. Actual inner raw swirl is half the square integral; the API uses twice that atom. Reference baseline, signed compact-bump contributions and actual inner offsets remain separate; offsets are reapplied once. Constant stages use exact exponential transfers and finite stages declared Gauss nodes on the installed continuous schedule. Absolute checkpoint guards avoid rejecting recorded endpoints after subtraction of enormous log radii.

Actual field checks passed: Rh seed nominal relative agreement about 3.21e-291 (not an input accuracy certificate); four forward moment integrands agree with actual velocity/velocity-Z in flattening with maximum relative difference about 6.83e-14. Separately retained tiny bump primitives agree with point correction integrands to about 1.11e-17. Samples include Rh, flattening, angular bump and Rtail. Quadrature and inherited inner uncertainty are still unenclosed. Complete five moments, heat continuation, finite energy, stress and recursion remain open.

Next: combine these angular rows with the same live axial row-2 and quadratic primitives plus actual inner seeds, producing mixed and z_theta moments/Z jets. Continue angular/energy primitives through heat using consistent H values and derivatives, then enclose inputs and propagate coefficient tangents. Preserve nonzero terminal radial residuals; do not infer finite energy from these local checks.

## 2026-09-30: continuous angular correction and pressure Z jet installed

Completed continuous_angular_correction.py with one continuous beta definition for point values, coefficient equations, partial pressure and full bump energy atoms. The actual field now includes amplitude*h_Z after flattening; tiny h/h_Z/h_t are retained separately. Coefficient Z derivatives use the implicit two-equation Jacobian, keyed by full MP Z strings. Parameter parsing preserves declared precision; atoms use 100 working digits, not a claim of 100-digit accuracy. Pressure exposes the backward bump contribution and its Z derivative separately from baseline pressure. Existing correction/schedule identities are retained and stale caches cleared. The energy target is regenerated from the current live tail nominal contribution instead of preserving the previous future term.

Actual shared-field run passed: five angular installation samples agree with the provider and the first bump retains nonzero Utheta_Z. Bump pressure Z derivative fourth-order difference relative error is 2.329932752163308432839506125565243891769e-20; Rp mass and Z matching remain about 6.92e-83. Pressure and energy providers confirm shared bump atoms. Standalone bump derivative, parity, support and endpoint checks passed. These are numerical consistency results, not tiny pressure-target cancellation, full pressure/moment closure, finite energy or recursive matching.

Next: continue full angular/pressure moments and their Z jets using these atoms; bound the inherited preheat baseline and heat H contributions; enclose actual inner/swirl/future energy inputs; propagate complete input intervals through coefficient tangents; resolve the actual nonzero terminal radial tail before any finite-energy or scale-recursion claim. Keep unrelated scale_reference edits untouched.

## 2026-09-30: continuous angular schedule installed in the shared field

Installed continuous_angular_schedule.py in the existing schedule object. Full pre-heat log A, switch slopes, Z flattening and public MP Z entry points now use continuous MP primitives. Interior primitive working precision is 100 digits, arithmetic 443; accuracy is not inferred from either. Angular field adapters consume the returned jets. Heat normalization and interpolation retain epsilon in MP while H remains inherited. Incoming log Ep is regenerated before solving. Five transition Z differences at 1e-4 and 1e-30 pass; heat interface log jump is about 4.97e-292; Rp mass and Z matching still pass around 6.92e-83. Shared runtime identities and nonzero terminal residual retention pass.

Next tasks: implement continuous relative angular correction and its coefficient Z jets; update angular cumulative and pressure moment providers together; enclose H and its energy contributions; bound actual inner offsets and inherited swirl/future energy; propagate full coefficient tangent intervals and resolve radial tail. No full five-moment, finite-energy, stress or scale-recursion certification.

## 2026-09-30: continuous incoming angular primitive and second-row bounds

Completed continuous_incoming_angular.py and continuous_incoming_angular_enclosure.py with exact J(1)=1/2 and full incoming mixed-row bounds. At 4096 panels mixed relative width is about 5.858e-4. regenerate_incoming now computes the nominal mixed factor with MP J instead of float-backed schedule._log_A; the input changes about -9.57e-17 relatively, at 100 working digits/order96 (not certified accuracy). Both reference matching rows now pass their quadrature uncertainty into the coefficient report, preserving inner offsets once. Updated field Rp mass and Z matching pass at about 6.92e-83; shared runtime identity checks pass and both terminal residuals remain nonzero.

Next tasks: install continuous J and consistent angular jets throughout the complete exterior schedule; update angular cumulative/pressure/energy atoms together; enclose measured inner offsets, inherited swirl and future heat contribution; propagate all Z jets through coefficient tangents. Incoming angular bounds are for the ideal source, not a certificate of the legacy complete angular velocity. Finite energy, stress and scale recursion remain open.

## 2026-09-30: reference incoming axial intervals installed in coefficient report

Completed continuous_incoming_enclosure.py with full-support bounds and analytic reference Z derivatives. At Md=.5 and Z=.3, 4096-panel relative widths are about 2.221e-4 for Iz and 2.156e-4 for Iuz2; the installed values are contained. The incoming provider exposes full_incoming_enclosure. The coefficient report propagates mass uncertainty into base1 and squared-velocity uncertainty with the correct negative sign into the target, preserving actual stored inner offsets once. Nominal amplitude and coefficients remain contained. Energy-target relative interval width from this component alone is about 9.49e-43; coefficient width remains pulse-dominated at about 2.49 percent.

Conditional gaps remain: mixed angular row2, inherited swirl energy, inner offsets, future heat energy, normalization-source uncertainty and their Z jets. Reference axial derivative intervals are available but not yet propagated through coefficient tangents. Next tasks: enclose the mixed angular primitive and its factor; enclose inner moment offsets and future energy; transport each source separately into base/target intervals and tangent equations; tighten pulse range bounds. No global mean, finite-energy, stress or scale-recursion certification.

## 2026-09-30: pulse energy uncertainty now propagated

Completed continuous_pulse_energy_enclosure.py, JSON and proof notes. Startup uses monotonic switch bounds to enclose its nested primitive; plateau is analytic; cutoff uses outward positive rectangles. Entire support is covered. At 4096 panels Kp relative width is about 4.10e-9 and the installed nominal value is contained. The coefficient interval solver now consumes this Kp range; amplitude relative width is about 2.05e-9, coefficients about 2.49 percent. Incoming rows and target remain fixed conditional inputs. No coefficients or terminal residuals were overwritten.

Next tasks: bound actual continuous incoming rows and squared-velocity energy; separate mixed angular uncertainty and target uncertainty; pass those intervals into the coefficient solve. Tighten correlated coefficient numerators where needed, then add incoming/target Z-derivative enclosures. Basis, pulse rows and Kp have bounds; full mean, finite energy, stress and scale recursion remain unproved.

## 2026-09-30: conditional coefficient interval solve

Completed continuous_solve_enclosure.py with exact endpoint JSON and notes. It propagates basis and full pulse bounds through the factored determinant and the weighted energy quadratic. The live runtime exposes coefficient_enclosure. Conditional unique positive root and nominal containment pass at 1024 and 4096 panels; coefficient relative widths narrow from about 10.36 to 2.49 percent. Incoming rows, target and Kp are still fixed stored dyadic parameters, with originating uncertainty explicitly unbounded. Fixed-coefficient residual intervals containing zero establish compatibility only; nonzero materialized terminal residuals remain.

Next tasks: (1) Enclose continuous incoming cutoff mass and squared-velocity integrals using a scaled endpoint transformation and explicit positive tails. (2) Enclose pulse energy Kp rather than using refinement differences as bounds. (3) Feed actual input intervals into solve_atoms and preserve each uncertainty source. (4) Tighten correlated inverse numerators if row dependency loss dominates. (5) Carry base/target Z derivative intervals into the tangent solve before certifying physical radial tail closure. Finite energy, all five moments, stress and scale recursion remain open.

## 2026-09-30: complete continuous pulse interval bounds

Added experiments/root_st073/lei_ren_part1_paper_continuous_pulse_enclosure.py and its JSON/notes. Both weighted rows have outward bounds for the full functional pulse support, including three positive omitted pieces. At 4096 panels the centered relative interval width is about 1.89 percent. The live shared runtime exposes pulse_enclosure and its provider passes independent containment. Basis and pulse integration are now bounded for declared decimal parameters; incoming/angular uncertainty, pulse energy and complete coefficient closure remain open. Terminal mass and Z residuals are still retained and nonzero; finite energy and scale recursion are not certified.

Next agent tasks: tighten the pulse interval width with analytic curvature or adaptive range panels; preserve exact dyadic endpoints and positive omitted support; propagate correlated basis and pulse bounds into the coefficient solve without subtracting nearly equal matrix products; enclose incoming rows and pulse energy before claiming full mean closure. Keep user changes in scale_reference files.

## Interval enclosures for continuous basis atoms - 2026-09-30

New continuous_basis_enclosure.py/md/json uses outward interval arithmetic
and an analytic midpoint remainder for the actual continuous bump definition.
It encloses I0, both weighted B rows/matrix atoms, and squared-bump energy Gram.
Exact dyadic endpoints are serialized; displayed decimal digits are not the
certificate. Real declared mu/ell are covered, inherited parameter errors not.

For 4096 panels: I0 width 1.216e-6, B1 width 5.706e-6; determinant width 1.477e-33.
The correlated analytic determinant identity (with positive sinh bounds)
proves this candidate's real basis matrix invertible despite tiny mu.
Nominal I0/matrix/determinant/Gram all lie inside independent intervals.
1024-to4096 panel refinement decreases the analytic remainder by 16.

SharedContinuousAxialRuntime now exposes optional basis_enclosure; its fixture
checks the actual owned basis against these bounds. Complete pulse quadrature,
incoming moments, coefficient errors and source target remain unenclosed.
The broad full-axial/global-energy certificates stay FALSE, and materialized
terminal residuals remain visible. This bounds basis integrals, not the whole
NS field or exact terminal cancellation.

Next: enclose pulse centered quadrature and incoming primitives; propagate
correlated atom/input errors through functional coefficients; establish both
terminal constraints without replacing materialized residuals. Full outer
moments/pressure, global energy, recursive matching and stress/remainder remain.

## Shared live continuous axial atoms installed - 2026-09-30

New continuous_axial_runtime.py/md/json owns ONE live basis/pulse, complete
matrix and pulse rows, numeric p, bump energy atoms and solved a/c. The point
component, algebra and complete primitive share objects by identity. Partial
integrals retain the same basis; complete end atoms use the stored matrix.

The actual ContinuousIncomingProfile now uses this owner for point values,
means and coefficient Z tangents. No coefficients or full integral atoms are
reconstructed from JSON inside this installed path. Kp retains its inherited
100-digit atom precision and open accuracy bounds; algebra uses 200 digits.

Full-field reruns: Rp mass/Z-mass matching 6.92e-83. The actual primitive and
independent full-atom replay (including Z derivative) now differ by zero at
reported precision. The SHARED materialized residual remains nonzero:
normalized mass balance 8.50e-202; Z balance 1.12e-201. Global finite energy
and recursive closure are NOT certified. No mass or derivative is reset.

Next tasks:
- [x] Share live complete atoms between solve, values and cumulative means.
- [x] Differentiate coefficients with the same live algebra and actual inputs.
- [ ] Enclose complete and partial integral errors, coefficient sensitivity and
      signed pulse omitted pieces without conflating them with exact identities.
- [ ] Represent and establish BOTH terminal constraints for the exact continuous
      functional coefficients, preserving materialized residual diagnostics.
- [ ] Finish angular jets, all outer moments/pressure and physical energy;
      recursive matching, stress/remainder and oscillatory corrections remain.

## Direct remaining end integrals installed - 2026-09-30

ContinuousAxialBump.tail uses reflection of the same even bump to integrate
remaining mass directly. Flat endpoint integrals factor out the endpoint
maximum before quadrature, avoiding absolute-error stopping on tiny values.
At precision100, s=.1499 retains log(tail)=-763.2991417 while full-minus-partial
rounds to zero. This is a numerical representation improvement, not closure.

ContinuousAxialCorrection now caches and RETAINS its evaluated terminal
balance rho. End-region cumulative rows use rho minus the direct remaining
end integral. The installed mean/Z-mean jet path uses the same representation
and retains rho_Z; no terminal mass or derivative is replaced with zero.

Actual full-field reruns: Rp mass and Z-mass matching remain 6.92e-83.
Installed normalized terminal mass balance is 1.11e-200; Z balance6.29e-201;
physical radial tail remains nonzero. Incoming scalar serialization now keeps
its declared100-digit quadrature precision, with separate source/algebra and
mixed-integral precision provenance; extra printed digits are not accuracy.

Next tasks:
- [x] Direct end-tail primitives and retained-residual cumulative representation.
- [x] Install representation in actual outer mean AND Z-mean jets.
- [x] Preserve incoming atom precision in regenerated receipts.
- [ ] Share complete continuous atoms with functional coefficient definitions,
      keeping evaluated coefficient/atom residuals and quadrature bounds apart.
- [ ] Establish both terminal constraints under that complete integral model.
- [ ] Finish angular jets, outer five moments and pressure, global energy,
      recursive matching, stress/remainder and oscillatory corrections.

## Actual terminal balance and energy-tail obstruction - 2026-09-30

New continuous_terminal_balance.py/md/json evaluates the INSTALLED continuous
incoming/exterior and actual input/coefficient Z jets, not a formal quadrature
graph. Installed terminal mass/sumabs balance is 1.73e-200; Z mass balance is
6.62e-201. Both remain nonzero. Independent full-atom replay differs slightly
from the cumulative primitive; that summation discrepancy is retained.

Recovered physical radial tail: ur=-nu*C(Z)/r, with
C=((1-delta)*Z*Mz+(1-Z^2)*Mz_Z)/(1-delta*Z^2).
The radial kinetic energy density per dZ dlogr is pi*nu^2*C^2*dz/dZ.
At Z=.3 the nominal C/Rh is nonzero, logabs about 1.1504127226124e28.
Thus increasing ordinary coefficient precision does not establish finite
GLOBAL energy: any surviving continuous 1/r tail makes that integral diverge.
Pulse omitted-piece bounds are transported separately and remain conditional
on nominal coefficients; other quadrature/input/coefficient errors are open.

Next: tie complete continuous integral atoms, coefficient constraints and
partial primitives through a cancellation-preserving definition; close BOTH
terminal mass and its Z derivative while preserving actual incoming seeds.
Do not assign mass zero, add an arbitrary radial mask, or use a favorable
alternate term replay as physical-field closure evidence. Angular primitive,
all five outer moments, stress/remainder and recursive matching remain open.

## Shared continuous incoming axial field installed - 2026-09-30

New continuous_incoming.py/md/json provides the same MP smooth cutoff for
Uz, Uz_y, Uz_Z and weighted cumulative mass/axial-energy integrals. It reuses
the pulse switch and evaluates the reflected cutoff to retain positive flat
tails. Incoming mean_Z is analytic; quadrature accuracy is not enclosed.

New continuous_incoming_outer.py/md/json installs that provider before Rp,
regenerates I_z, I_theta_z and I_uz2, reapplies measured inner offsets ONCE,
updates the prior-energy target and re-solves with retained continuous atoms.
Both incoming point values and mass/Z-mass now feed the joined field.
Rp mass and mass_Z relative matching: 6.92e-83 each. New linear row replay:
1.27e-173 / 1.24e-173; energy replay 1.63e-201 (numerical algebra only).

A float-rounded amplitude derivative produced a 2.26e-17 derivative jump.
Using the exact MP derivative of (1+Z^2)^-1 removes that discrepancy.
Zero-Z tangent factors are taken from the regenerated incoming receipt.

Remaining limits: the mixed integral retains the float-backed angular
primitive, swirl energy is inherited, future angular target derivative is
still finite-differenced, and integral/coefficient uncertainty is unenclosed.
No terminal mean is assigned zero. Global finite energy, full five moments,
NS stress/remainder and scale recursion remain incomplete.

Next actionable tasks:
- [x] Install continuous incoming axial values and matching mass/Z-mass.
- [x] Regenerate axial incoming rows/target and install updated coefficients.
- [ ] Replace retained angular primitive and cache/jet precision coherently.
- [ ] Resolve terminal mean AND its Z derivative with explicit integral bounds.
- [ ] Continue all five moments and pressure jets through outer/heat layers.
- [ ] Establish finite global energy and cross-scale matching before claiming
      recursive background closure; then stress/remainder and pulse corrections.

## Actual axial coefficient and cumulative mean Z derivatives - 2026-09-30

New continuous_axial_tangents.py/md/json differentiates the same continuous
linear rows and quadratic energy constraint. Independent declared-direction
fourth-order replay refines 2.19e-17 -> 1.37e-18; actual nominal input jets
replay derivative rows around 1e-201. Arithmetic replay is not a certificate.

New seeded_input_tangents.py/md/json transports actual corrected Rh momentsZ
and reference/source-factor derivatives into normalized incoming rows/target.
It does not recover tiny reference inputs by subtracting a large seed. Swirl
reference energy derivative cancels analytically after Ep normalization.
Unresolved inner entries [3,5], subtraction and quadrature uncertainty remain.
Future angular-energy derivative is explicitly float-backed fourth-order FD;
its relative two-step change is about 3e-12, not an error enclosure.

New continuous_seeded_outer_jets.py/md/json installs post-Rp Uz_Z and mean_Z
from those input/coefficient jets. Local radial divergence refines 6.24e-24
-> 3.90e-25 without neighboring full coefficient solves. Pulse integrals are
separate from cumulative means; each seed/amplitude is counted once.
Schedule/angular/tail/pressure identities and nonzero terminal rows remain.
Incoming before Rp retains a declared finite-difference fallback; angular
velocity Z jets remain incomplete. Signed pulse omission bounds are separate
from incoming/coefficient/quadrature uncertainty.

Next: analytic angular coefficient/heat target derivatives; shared continuous
incoming primitive and Z cache precision; actual mean and mean_Z closure;
outer five moments/jets and global energy; recursive matching, stress/remainder
and oscillatory correction. Global finite energy/scale recursion are not done.

## Continuous axial correction installed in joined exterior - 2026-09-30

New continuous_seeded_outer.py/md/json adapter installs continuous Uz AND
its cumulative seeded mean together. Original schedule/angular/tail/pressure
objects are preserved. Seven fixed-Z positions cover Rp, startup, plateau, cutoff,
first end bump, Rv and after-Rv. Post-Rv constant mass replay is 2.10e-443.
Actual Ur is recovered from the same mean. Independent radial divergence
relative cancellation refines 6.24e-24 -> 3.90e-25 as step halves.
These are local diagnostics; float-Z derivatives remain uncertified.

Startup weighted primitives now use endpoint scaling and integration by
parts, retaining the shared smooth sigma definition. Truncating a SUBTRACTED
sigma integral produces a negative correction: absolute omission bounds and
sign are propagated by the component and installed mean. Actual-mu startup
xi=.015 bound transport passes. Numerical loss of a positive primitive raises;
no nonzero terminal mean is replaced by zero. Quadrature and inherited input
uncertainty are not covered by omitted-piece bounds.

Next derive shared continuous incoming and Z derivatives (avoid float cache
loss), certify actual terminal mean/Z-mean closure, and finish exterior jets,
five moments and energy. The installed finite numerical tail remains nonzero;
full finite energy, recursive matching, stress/remainder and oscillatory
closure are still open. Do not equate tiny local divergence with NS closure.

## Continuous cutoff primitive and shared cumulative component - 2026-09-30

Completed cutoff partial integral branches for xi in (10,11): endpoint,
saddle and after-window evaluation, with positive omitted-piece bounds.
Four actual-candidate positions run; independent endpoint derivative replay
refines 3.19e-12 -> 7.98e-13 when step is halved. Startup (0,.02] remains
unsupported explicitly. Quadrature error is not enclosed by omitted bounds.

ContinuousAxialCorrection now uses the solved pulse/end coefficients for
point values and normalized cumulative rows, with the SAME continuous atoms.
Actual terminal nominal values remain nonzero; no mean/tail reset occurs.
This component is not yet installed in the global physical field.

Injection audit: preserve original schedule/angular/tail/pressure identities;
carry the entire seeded incoming receipt, not just linear base rows. Its
row_normalization is seeded but dimensionless_integrals.I_z remains unseeded.
Pre-Rp mass seed must be added once. Convert stored full-row pulse weight to
current-point weight by exp[lambda*(13-xi)/mu]. Replace velocity and primitive
providers together. Details in continuous_axial_solve.md.

Next complete startup weighted primitive and shared incoming/Z provenance;
install consistent axial values/means, then establish terminal mean AND its
Z derivative closure. Global finite energy and scale recursion remain open.

## Continuous pulse energy and re-solved axial correction - 2026-09-30

Completed: continuous_pulse_energy.py/md/json and continuous_axial_solve.py/md/json
under experiments/root_st073. Kp uses the same continuous pulse definition,
MP startup/cutoff nodes, and analytic plateau. Same-order 70/100 precision
refinement is 7.54e-72; order80/112 refinement is 2.76e-21 (not an enclosure).
Fixed an import-time decimal-boundary precision leak. Old float Kp differs
by 9.84e-15 relative. Working digits do not equal certified integral accuracy.

Re-solved actual seeded base rows/energy target with continuous bump, pulse
and energy atoms. Old coefficients have 6.68e-15 relative continuous-row
defects; new row arithmetic replay is about 1e-173. Energy algebra replay is
about 1e-201. Consistent solved end-bump value/derivative/weighted primitive
provider is available. These coefficients are not installed globally yet.
Continuous pulse partial integrals now cover the analytic plateau in addition
to the saddle window, retaining an explicit positive startup bound.

Next: complete startup/off-window cutoff and incoming primitives; derive
shared Z derivatives; install continuous coefficients in both point velocity
and cumulative means together. Establish actual mean/Z-mean closure before
any tail removal. Then exterior energy, outer jets/five moments, recursive
transitions, stress/remainder and oscillatory correction. Full goal stays open.

## Continuous integral providers and bounded core energy - 2026-09-30

New artifacts: lei_ren_part1_paper_continuous_axial_basis.py/md/json,
lei_ren_part1_paper_continuous_axial_pulse.py/md/json,
lei_ren_part1_paper_core_energy.py/md/json, under experiments/root_st073.
The continuous bump shares values, derivatives and full/partial primitives.
Its matrix differs from the old float-node canonical matrix by 2.21e-19.
The MP pulse shares pointwise values and full weighted rows; saddle-window
partial weighted primitives now include positive omitted-piece bounds and
an independent derivative diagnostic. Off-window partial evaluation raises
an explicit error and remains to be implemented. Pulse full-row log inputs
change by 6.68e-15 from the float-centered evaluator. No global installation.

Actual bounded core energy integration confirms declining kinetic energy
under the selected shrinking time scales, despite growing velocity amplitude.
This is core-only; the exterior nonzero radial tail remains the global energy
obstruction. Do not infer recursive scale closure or full NS residual closure.

Next tasks in order:
1. Complete off-window continuous pulse partial primitives with bounds.
2. Add continuous pulse energy and shared incoming primitive/Z derivatives.
3. Re-solve with the same continuous bump/pulse/incoming atoms used by velocity.
4. Establish actual terminal mean AND its Z derivative cancellation before
   modifying any 1/r tail; keep nonzero arithmetic tails visible.
5. Complete outer jets/five moments, annular/exterior energy and measured
   core-width profiles; then uniform stress/remainder and oscillatory layers.
The full goal remains active; local core energy is not global finite energy.

## Exact quadrature dependency graph and measured multi-time core geometry - 2026-09-30

Run lei_ren_part1_paper_axial_closure_graph.py for actual shared candidate
inputs. Canonical matrix and row atoms are now exported by the axial solver.
Exact rational Cramer expressions eliminate both quadrature row terms;
their formal Z derivatives also cancel for a fixed matrix and differentiable
shared atoms. This is a quadrature-model identity, not continuous closure.
Independent basis perturbations produce nonzero residual terms, not zero.
Rounded graph coefficient/energy replays remain separate. The physical
nonzero tail is unchanged; continuous integral/Z-derivative provenance and
the exact energy root still need work. Read closure_graph.md for the next
shared continuous/conservative basis requirements.

Run lei_ren_part1_paper_core_scale_geometry.py. It uses actual regular-core
physical callable roundtrips at logq=-alpha/delta, alpha=0,2,4,6. Measured
aspect gains are 1,e,e^2,e^3; swirl/axial amplitudes, local winding density
and axial vorticity component grow with the mapped scale laws. Log-aspect
errors are about 1e-244 and relative velocity roundtrip errors about 1e-247.
These extreme MP chart times are not ordinary simulation times. One core
point diagnoses geometry/coordinate scaling; it is not a measured vortex
core width, an integrated streamline, full vorticity, global energy, or a
completed recursive transition. Both Python/JSON/Markdown artifacts exist.

Next connect a shared continuous basis for values, partial/full primitives,
and Z jets; prove the actual mean and its derivative close before removing
the 1/r tail. Independently integrate energy over shrinking core/annuli and
exterior, measure core-width profiles across time, and finish outer jets and
five moments. Uniform stress/remainder and oscillatory correction remain
open. Do not treat formal row cancellation as a finite-energy certificate.

## Installed seeded axial field; mean conditioning exposed - 2026-09-30

Read lei_ren_part1_paper_seeded_outer_field.py/md/json and the regenerated
seeded_shared_candidate.json. New coefficients are now installed in both
axial velocity and the actual cumulative mean. Schedule/angular/pressure
identities are retained. Direct seeded transport avoids subtract/add loss.
Fixed Decimal stage-offset loss in axial_incoming._row, MP bump normalization,
and stored full pulse row reuse. Prior joined_outer.json is labelled historical.
At Z=.3 Rp relative mean jump=3.541e-260; Rv reported mean jump=0 at precision;
pulse independent divergence relative cancellation=3.899e-25. These are local
numerical diagnostics, not full momentum residuals or uniform certificates.

IMPORTANT: tail coefficient ratio is -2.590e-183, but log(|tail/Rh|) remains
1.1504e28. The numerical tail is still nonzero and the implemented radial
1/r energy obstruction remains. Increasing ordinary MP precision is not the
next solution: derive an exact closure identity for the same actual field,
with continuous integral/conservative-model provenance and uncertainty.
Do not set a nonzero terminal mean to zero or hide an arithmetic residual.
Then restore corrected outer jets, full five moments, energy/support tests,
actual multi-time scale diagnostics and stress/remainder. Oscillatory layers
and full residual <=1e-3 remain later requirements.

## Shared inner-through-heat callable and actual seeded exterior solve - 2026-09-30

Run lei_ren_part1_paper_joined_outer.py. JoinedOuterField keeps the SAME source
profile and the actual five-bump corrected Rh mass/pressure offsets. It now
provides a physical velocity callable from the regular core through heat.
The Z=.3 heat-chart roundtrip passed at working precision. This is a
provisional continuation: sampled nonzero radial transport still prevents
claiming finite energy. Float-Z derivative errors, missing corrected outer
velocity jets and full five-moment exterior stress remain explicit.

Run lei_ren_part1_paper_seeded_axial_inputs.py --shared to reproduce the actual
Rh-seeded Section 7.31/7.34 coefficient solve (seeded_shared_candidate.json).
Canonical linear row replay is 5.924e-416 and 1.945e-416; energy replay is
2.870e-444. These are finite-quadrature coefficient equations, NOT full NS
residuals or a continuous mean closure proof. New coefficients are not yet
installed in the joined velocity/primitive. The default command is a clearly
labelled historical fixture regression, not the actual shared candidate.

Next implement a seeded axial pulse/primitive provider together: before Rp,
retain old outer cumulative mass plus actual Rh offset; from Rp onward use
updated incoming m1 and the newly solved pulse/end-bump coefficients. Avoid
double-counting the inner offset in the joined radial transport. Recompute
Z jets from the same actual primitive, update incoming integral provenance,
measure terminal mean by replay rather than assign zero, and quantify
quadrature/precision refinements. Preserve shared swirl and axis pressure.
Then resolve finite energy and axial support before claiming temporal scale
recursion. Global norms, admissible stress/remainder and oscillatory layer
remain open.

## Actual inner seed transport and serialization precision - 2026-09-30

New lei_ren_part1_paper_seeded_axial_inputs.py maps explicit raw inner offsets
(z, theta_z, z_theta) into the two normalized Section 7.31 linear rows and
the Section 7.34 prior-energy input. Use the same schedule log_Rp and log_Ep.
No absent offset defaults to zero. This is preparation for a new outer
coefficient solve, not a completed repair or a replacement axial primitive.
The copied dimensionless_integrals remain historical reference inputs;
update the actual primitive consistently before using a seeded receipt in
CorrectedSourceProfile. A dimensional replay at logRp=5e152+17 passed.

from_signed_log now preserves arbitrary_exponent_value when available,
validates its sign, and retains the legacy log-only fallback. The stored
precision receipt shows relative error 2.24e-161 versus 2.23e-8 for a
log-only roundtrip. The existing actual-tail solve/replay at Z=.5 passes
(linear differences below 5.80e-133, energy difference 6.85e-160).

Next: complete the same-object inner/outer join, measure terminal transport,
then feed actual offsets into the pulse/end-bump coefficient solve and replay
the same primitive. Do not zero the residual mean. Finite energy, uniform
source constants, temporal scale recursion, and oscillatory residual remain
unestablished.

## Five-bump numerical correction and physical inner callable - 2026-09-30

Read experiments/root_st073/lei_ren_part1_paper_inner_corrected_field.py/md/json
and inner_moment_map.py/json. Equation (10.8) now solves all five representative
nonlinear moment equations with analytic coefficient Z derivatives. Actual
partial moments, pressure with unchanged P0, and continuity-derived Ur are
connected to the shared core-to-Rh candidate. Unknown defect entries3/5 use
explicit midpoint representatives and conditional value/Z intervals; they
are not claimed zero. Exact normalized bump mass identities are retained.
At Z=.3 eleven center/flank probes pass the relaxed cone. Higher quadrature
replay of the terminal representative equations gives normalized residuals
around 1e-33 or lower, distinct from the much smaller algebraic solve residual.
The local mapped-divergence flank replay gives q div(u) about 1.57e-24.
CorrectedInnerField.physical_field().velocity(x,y,z,t) now covers the actual
regular core through Rh and raises outside the constructed domain.
Next audit terminal fields/moments over Z, resolve or enclose input uncertainty,
then append the supplied corrected outer/heat profile with actual mean
identities. Tiny nonzero mass tails cannot imply finite energy. Uniform
source constants, inner admissible collar, global moment closure/energy,
temporal scale recursion and oscillatory-corrected residual remain open.
Detailed ordered tasks are in inner_corrected_field.md.

## Axial reference restored; actual defect inputs recorded - 2026-09-30

Read experiments/root_st073/lei_ren_part1_paper_axial_restore.py/md/json.
Section 9.38 transports actual moments and Z jets from Rz through Rh.
At Z=.3, V reaches 4Z=1.2, V_Z reaches 4, and V_y vanishes at the end.
Five relaxed-cone probes pass. Largest 32/64 moment quadrature difference
is 1.51e-32. Independent radial mapped-divergence replay at restoration
midpoint gives q div(u)=2.14e-24, relative cancellation 3.26e-25.
Centered defects d1,d2,d4 are about 2.234e-15,5.690e-16,5.915e-41.
Angular/pressure entries d3,d5 are unresolved by subtraction: nonzero
roundoff artifacts exceed their conditional physical bounds. Do not set
them to zero or use them as measured defects. The receipt preserves
conditional bounds and flags the missing uniform C1/e_star certificate.
Next implement the fixed Section10.8 five-bump matrix and nonlinear
moment correction, with stable signed-log/interval treatment of tiny
defects and actual partial bump moments. Detailed checklist is in the
axial_restore.md. Source requires relaxed cones in these later intervals;
strict admissibility is needed in the inner collar. Finite energy,
heat exterior matching and temporal scale recursion remain unfinished.

## Long reshape now reaches the axial restoration entry - 2026-09-30

Read experiments/root_st073/lei_ren_part1_paper_long_reshape.py/md/json.
Section 9.30 is connected to actual shared-candidate R=110 moments and
Z jets. Endpoint-normalized MP quadrature reaches Rsh and analytic
reference-power primitives continue to Rz=exp(-8)Rref. Five Z=.3 probes
pass the relaxed cone; full admissible stress remains unestablished.
Largest 16/32 normalized quadrature difference is 6.57e-25. The angular
reference log-value difference at Rsh/Rz is about 5.59e-263. Inherited
inner moments are preserved, and conditional value/Z tail bounds are
recorded separately from uncertified quadrature error.
Next implement axial restoration on [Rz,e Rz] with actual five moments
and Z jets, continue to Rh, compute actual Section 10 repair defects,
and implement moment repair. Complex A_Omega, mixed core A/K, pressure
Z-tail bounds, heat exterior, finite energy, admissible stress and temporal
scale recursion remain unfinished. Provisional source budgets are unchanged.

## New shared candidate reaches R=110 - 2026-09-30

Read experiments/root_st073/lei_ren_part1_paper_shared_candidate_1.md/json.
The reproducible Python script recomputes pressure and degree 18 core at
j=1e-14, Lambda=1e36, logCstar=5e151, logPstar=14, delta=1e-200.
All five necessary input gates and six sampled core signs pass. The four
new-parameter connection probes pass the relaxed cone; stricter admissible
cones are not established. Explicit phases resolve distinct switch shears
although the physical radius offsets round away. R=110 has a=.8,b=0.
Real-axis norm bounds are available; mixed core A, full K and complex
A_Omega bounds remain open. A=1e150 and logK_upper=1e152 are provisional.
Next implement long angular reshaping with actual moments/Z jets while
bounding these source constants, then axial restoration, moment repair
and heat exterior matching. Finite energy, temporal scale recursion and
full oscillatory-corrected residual remain unfinished. The detailed ordered
checklist is in shared_candidate_1.md; mark tasks complete only with receipts.

## Short source switches now reach R=110 - 2026-09-30

ExitSwitches implements both exact short switch intervals, actual moment
and Z-jet transport, and analytic a=.8,b=0 continuation to R=110. At Z=.3,
R=100 field/moment/jet/slope matching is 1.60e-160 relative. All six sampled
relaxed cones pass. The 16/32 switch endpoint difference is 7.03e-12; the
largest interior difference is 2.44e-9. The recorded b=0 cone margin gives
D(R=110) about 5.68e37, above 3. Read exit_switches.md/json under
experiments/root_st073. These remain development-fixture diagnostics.

Next priority is the shared parameter rebuild, not blindly extending
this fixture to Rsh. Its necessary parameter gate explicitly fails.
Follow connection_scale_gate.md for selecting j, controlling A/K,
choosing Cstar/Rref/delta and a common source collar, and recomputing
pressure/core. Then replay the reusable exit/switch modules, reshape
angular velocity, restore axial velocity and repair actual five moments.
Full matching, finite energy and temporal scale recursion remain open.

## Shared source parameters must be rebuilt - 2026-09-30

The current development fixture fails necessary Section9 input conditions:
Cstar/A, long angular-shape radius, axial-radius separation, j and short
collar width. Read connection_scale_gate.md/json under experiments/root_st073.
This is an explicit necessary-condition rejection, not just missing proof.
Real endpoint max(G) bounds A; it is not the complex A_Omega bound.

build_source_core now accepts shared j/logC/logPstar/delta and recomputes
outer schedule, pressure anchor, axis and core together. Its nondefault
wiring receipt passes local checks, but does not certify a replacement.
Do not reuse old fixture receipts after changing parameters. Next select
one compatible shared candidate with controlled A/K, then rerun pressure,
core, switches and long reshape. Existing local source modules remain
useful algorithm implementations; full matching and recursion remain open.

## Exit continued to R=100 - 2026-09-30

The frozen comparison now extends beyond the core and supplies exact
algebraic D/I_z drivers; independent formula replay through R=110 agrees
at about 1e-160 relative. ExitContinuation propagates actual endpoint
values and five moments to R=100, retaining positive prescribed shear
and conditional subprecision velocity/moment error bounds. Six samples
(two Z sections, R=1,10,100) pass relaxed cone tests. This does not
certify uniform Z jets, source constants or a completed outer match.

Read exit_continuation.md/json and connection_next.md under
experiments/root_st073. Next implement the source short shear switches
on 100..110. Before the long reshape, satisfy the shared source parameter
and radial-order gates; the existing lower-bound logCstar candidate is
not certified for those inequalities. Recompute shared pressure/core
when changing parameters. Temporal recursion and global energy remain open.

## Actual exit Z jets and three-component callable - 2026-09-30

ExitTangents now transports actual five-moment Z jets and recovers Ur
from the same Mz/Mz_Z. The shared stress evaluator accepts exact ODE
shears to preserve the extremely small terminal shear. All five starting
moments and their jets match the core. At Z=.3,y=.0025, independently
differenced relative divergence cancellation error is 3.95e-20; the
16/32-step terminal Ur relative difference is 6.05e-19. Driver Z-step
uncertainty is checked separately. These are finite local diagnostics.

The midpoint passes sampled admissible/relaxed cone tests; the endpoint
passes relaxed only. The core-join strict test has numerically zero
T dot S and is retained as false. The strict condition excludes the Ra boundary. Near-join interior margins
and uniform cone control remain unverified.
LocalCoreExitField.velocity(x,y,z,t) now covers core and initial exit;
physical roundtrips at two times pass. Full matching, finite energy and
temporal recursion remain open. See exit_field.md/json and
exit_callable.json under experiments/root_st073.

## Initial actual inner exit - 2026-09-30

The same source-scaled core now feeds Section 9.23 comparison and actual
Section 9.25 exit modules. F, Uz and five moments are jointly integrated
on 0<=log(R/Ra)<=.01, Ra=4/Lambda, Lambda=1e36. The actual-ODE receipt
compares 16/32/64 steps with fixed comparison discretization. At Z=.3,
the 32/64 log-amplitude difference is 3.05e-11 and the largest normalized
moment difference is 3.85e-11. This is an initial exit, not full matching
or temporal scale recursion. Actual Z-jet transport and stress/cone
are next. See experiments/root_st073/lei_ren_part1_paper_exit_bridge.md.

The future-pressure module supplies an uncorrected tail Z envelope.
Angular bump coefficient derivatives remain unresolved; corrected-tail
C1 bounds, source collar constants and global finite energy remain open.

## Source outer replacement and first correction inputs — 2026-09-29

2026-09-30: high-degree core continuation now reaches s=Lambda R=4.1
at sampled Z, with a same-polynomial velocity/pressure/five-moment/stress
adapter. AxisPressureJets supplies analytic dominant pressure jets from
the actual preflatten/tail decomposition; its normalized P/Pstar^2 and
physical units are explicit. Future tail values and heat bounds remain
separate; future-tail Z derivative bounds are NOT established.

The historical Lambda=1e18..1e24 nearby rounded-root point genuinely has
F_R>0 even though the exact root passes. It is preserved as an intentional
failure probe in `lei_ren_part1_paper_core_ra_experiment.json`. Choosing
Lambda=1e36 suppresses that failure in the tested narrow window
DeltaZ=c*|U1|/(Lambda*|H0'|), c=-16,-4,-1,0,1,4,16. The actual anchored
pressure adapter's 45 samples through s=4.1 retain F>0,F_R<0; its pressure
ratio U1^2/(Lambda*sigma0^2) is 2.55e-6. Prototype degree18/24 narrow-layer
differences are below displayed precision; at Z=.3 the degree18 nominal
core stress/shear angular ratio is at most 1.12e-28. These are finite
polynomial/sample diagnostics, not a uniform theorem, admissible global
cone or time-scale recursion. Next establish tail derivative/domain bounds,
check the core across the remaining Z domain, then perform Section 9/10
connection and moment repair using these actual core inputs.

2026-09-30: explicit schedule/correction/tail construction now supports a
shared source-scaled candidate with logRref=301.1814 (j=.02, Lambda=2500,
logCstar=2logLambda, logPstar=14). Its same-axis pressure is recomputed,
not inherited from the old independent core. The A_Omega=0 choice is only
a necessary lower-bound demonstration, not the full theorem parameter regime.
`lei_ren_part1_paper_axis_jets.py` implements the exact regular axis data
and Eq. (8.7). `lei_ren_part1_paper_core_recursion.py` implements the
nonlinear Eq. (8.2) radial Taylor recurrence; the actual-pressure third-degree
prefix is saved in `lei_ren_part1_paper_shared_pressure_core.json`.
At Z=.3 and R=1e-15, 5e-16, 2.5e-16 the local equation error exponents are
3.0000 angular and 3.00034 axial. The first axial slope agrees with the
independent Eq. (8.7) callback to 2.81e-13 relatively. This is a LOCAL
regular-core prefix, not time-scale recursion: extending to Ra=4/Lambda,
selecting Lambda from the actual pressure/contraction data, matching to Rh
and repairing the connection moments remain open.

2026-09-30: `SourcePulseMoments` now combines all five actual cumulative
moments from the reference through Rv, with actual source-pulse axial,
mixed and quadratic integrals and forward angular/energy propagation.
The bulk quadratic radial identity defect is 8.33e-15; the two linear
moment identities give 4.17e-14 in the bulk and 1.40e-11 at the end bumps.
The same-candidate bulk stress is saved in `lei_ren_part1_paper_pulse_moments.json`.
This does not close later flattening/heat moments or certify the cone.

A concrete core-scale obstruction is now explicit. Section 9 requires
Rref=110*(Cstar*Pstar)^10; the default demo logPstar=14/logRref=10 is
incompatible with the Section 8 core regime. For j=.02, Lambda>=2500
and even A_Omega>=0 imply logRref>=301.1814. The regular-core equations,
matching jets and Section 10 moment repair are extracted in
`experiments/root_st073/lei_ren_part1_paper_regular_core_plan.md`.
Do not reuse the old Lambda=10 core as the source candidate's core.

Newest implementation adds actual five cumulative moments through the initial
axial turnoff and arbitrary-exponent source stress from those moments and
the same candidate pressure. See `lei_ren_part1_paper_reference_moments.py`,
`lei_ren_part1_paper_mp_stress.py`, and `lei_ren_part1_paper_source_stress.py`.
The axial coefficient solve and end-bump primitive now share their effective
quadrature order and canonical full-support nodes. No nonzero exterior mean
is overwritten with zero: finite quadrature consistency is not exact
continuous closure and does not establish global finite energy. Later-stage
mixed/angular moments, regular axis core, stress cone and recursive PDE
coefficients remain open.
The new actual stress receipt at offsets -1, .5 and 2 with Z=.3 gives
maximum relative defect 1.06e-10 in I_theta_y+I_theta=N_theta and
I_z_y+I_z/2=N_z. This verifies those local moment equations, not an
admissibility cone or the full NS residual. Canonical axial linear replay
defects are about 3.1e-105 and 4.2e-106 for the inherited numerical inputs.

Latest unified candidate: `lei_ren_part1_paper_corrected_profile.py` now
provides shared angular/axial coefficients, a source axial primitive and
dependent radial velocity. It exposes Cartesian `velocity_from_tau` with
streamfunction axial localization and signed-log arbitrary-exponent output.
Independent bulk M_R=Uz error is 2.09e-15; radial Z-stencil refinement is
1.48e-14. The corresponding pressure adapter retains baseline and signed
bump correction separately. The baseline pressure now propagates at arbitrary
log radius; physical pressure shares the velocity chart and axial cutoff.
This temporary-reference candidate remains
unaccepted: regular axis core, global mean tails/energy, exact pressure
target, stress cone and recursive PDE corrections are unresolved.

Latest axial progress: actual incoming moments, row-normalized pulse
integrals and corrected angular energy tail now feed the source affine and
quadratic energy solve. a_p=1.0100502663 is within the source interval;
independent end-bump quadrature replays both normalized linear rows to
2.45e-19 relatively. This is algebra/quad evidence for the actual numerical
inputs, not PDE accuracy: incoming moment refinement is about 3.8e-7,
pulse refinement 3.38e-9 and angular-tail refinement 5.06e-13, with separate
heat bounds. Source profile assembly, common-pressure/stress recomputation,
coefficient derivatives, relaxed cone and new regular inner core remain open.

Latest: source angular coefficients now use arbitrary-exponent signed logs
and a scaled small-branch quadratic solve, with a callable multiplicative
angular bump representation. Actual bump quadrature replays its normalized
angular increment to 5.11e-15. This does NOT resolve cancellation at the much
smaller pressure-target scale, certify the approximate heat inputs, or remove
waiting-root uncertainty. The fixed axial pulse is also implemented:
Kp=.24504962020069448, source .24<Kp<.246, refinement 2.03e-15.
Actual axial RHS/tail energy, corrected pressure, derivatives/cone and the
new regular core remain required before accepting the source outer field.

Read `LEI_REN_SOURCE_OUTER.md`. The actual leading remainder now has an
axial-viscosity breakdown. At the sampled regular-core point its angular
and axial components are predominantly the source axial-viscosity terms;
diagnostic subtraction leaves the radial source rather than a flat field.
`lei_ren_part1_first_order_sources.py/.json` materializes Omega0 and the
known first-order forcing from the SAME core, with analytic radial polynomial
jets and independent physical comparisons at three scales. First-order
F1, Uz1 and P1 remain unsolved and must be coupled, extended and moment repaired.

The source outer schedule and exact Section 7 closure kernels are separate
from the old weak-anchor field. They require actual moment/pulse binding,
waiting-length closure and a new regular common-pressure core. No automatic
inheritance of source theorem constants or current-field diagnostics is allowed.

Source velocity-schedule checks now pass, including a 212-digit checkpoint
case preserving the 100-unit flattening and 3-unit heat connection after
an enormous logarithmic radius. Same-candidate normalized backward pressure
has a separate omitted-tail bound. The pre-heat waiting equation is also
numerically solved: declared example tau_wait=82.448458, independently
replayed normalized angular defect 2.68e-10. This waiting length is distinct
from physical time-to-blowup tau. The source axial closure also has stable
row scaling down to mu=1e-28 (condition about 4.32), preserving precise RHS
data and stage-local gamma offsets. Actual heat bumps, axial pulse moments,
all-five-moment binding and a new regular shared-pressure core remain open.
The source demonstration does not replace the callable field's h=.001 manifest.

Heat-replacement r_H and s_H now have candidate-derived logarithmic inputs
in `lei_ren_part1_paper_heat_defects.py/.json`. Their positive leading values
are retained with analytic Taylor remainder bounds instead of underflowing
to zero. Independent Gamma-expectation checks pass; pressure-collar integral
refinement is 4.24e-11. This is bounded heat-defect evaluation, not the exact
angular coefficient solve. Incoming Z-dependent r_pre is also retained by
the actual flattening difference ODE followed by logarithmic mu^(30(1-mu))
propagation. The next task is a coefficient representation that preserves
these tiny inputs, with waiting-root and quadrature uncertainties separated.

## Actual stress path and continuous-moment blocker — 2026-09-29

Continuous axial repair now re-solves the interpolated moment system at
each queried Z, with analytic derivatives of that same solve. The default
moment grid is 257 nodes; independent actual integrals give maximum mass
5.72e-14, mixed 8.36e-11 and quadratic 2.17e-11 residuals. The old direct
coefficient-interpolation receipt is retained separately. This removes one
finite terminal-moment obstruction, not the whole cone or recursion.
The validated replay seed avoids rebuilding all 257 radial integrations;
the physical field's default loader checks its pressure-receipt hash.

The updated physical-field receipt uses this same 257-node moment solve.
Its 12-point finest Cartesian divergence is 1.76e-7 (relative 1.70e-9),
axis parity errors are zero, and the finest-time full energy estimate is
0.6904 with 1.43% coarse/refined difference. These are finite diagnostics,
with imposed similarity exponents, not recursion or PDE acceptance.

After correcting the moment-derivative stencil, the actual collar/heat-tail
stress maximum is 4.00e-10, down from 1.972e-3 in the historical candidate.
This still exceeds quadrature-only support; derivative error must also be
accounted for before assigning a physical meaning to the remaining tail.
The transition remains inadmissible and does not supply the relaxed-cone
input required by the Section 11 shear loop. Tiny heat-edge stresses must
not be used as robust directional-cone evidence.

`lei_ren_part1_physical_remainder.py` now measures physical R_B, D T_B,
and E_B from the SAME callable velocity, pressure and actual moment stress.
It selects the explicit envelope B(z)^2, retains all velocity/pressure
cutoff derivatives in physical finite differences, and preserves radial
E_B=R_B because the completed tensor's radial divergence is zero.
The six fixed-sector sample receipt finds growing leading-order remainders,
not flatness: sampled core axial E_B is about 7.68, 61.32, 490.07 at
tau=1/8,1/32,1/128. Spatial/time step refinement is recorded; these are
neither global maxima nor volume-L2 norms. Implement compatible leading
transition and actual higher coefficients before claiming remainder decay.

The profile-derived stress evaluator is now implemented in
`experiments/root_st073/lei_ren_part1_profile_stress.py`. It evaluates the
source inertial terms and radial shear from one profile's actual five
cumulative moments, their Z derivatives, and its pressure. Independent
exact-heat and inward-collar checks validate this algebraic path; they do
not certify the assembled extended background.

The historical first assembled-profile stress diagnostic found only 4 of 27 sampled
points satisfying the cone and a heat-tail stress component of 1.972e-3,
compared with stress quadrature refinement of 1.34e-13. This receipt is a
failed candidate diagnostic, not a PDE pass. Direct interpolation of axial
repair coefficients left off-grid mixed/quadratic moment defects. The
continuous solve above addresses those defects. A later review also exposed
a reversed fourth-order moment-Z derivative stencil in the stress adapter;
its receipt is retained as a failed historical diagnostic. Recompute the
current stress with the corrected stencil and independent derivative checks
before attributing its heat-tail stress to a physical obstruction.

The extended swirl pressure is now integrated inward from the exact collar
pressure, avoiding cancellation against the larger axis datum. Its independent
radial pressure identity receipt has maximum relative error 3.28e-8.

After continuous moment closure, remeasure actual heat-tail stress and the
whole connection cone before implementing a higher-order remainder correction.
Do not replace measured terminal moments by theoretical heat targets inside
the diagnostic. Geometry fits and divergence passes do not establish scale
recursion, admissibility or the final full momentum gate.

## Extended exterior boundary data — 2026-09-29

Actual axial/mixed/quadratic repair is now bound to the rebuilt shared-pressure
swirl at five Z slices (-1,-.5,0,.5,1). Independent 128-node moment residual
is below 1.59e-13. A degenerate quadratic-root branch was repaired after
actual-profile integration exposed it. Next implement a smooth Z-dependent
repair with pinned branch orientation, between-slice moment checks and
dependent radial recovery. This is not yet the replacement physical field
or the whole admissible stress cone.

The actual extended swirl now has a rebuilt common-pressure finite core:
`lei_ren_part1_extended_pressure_core.json` converges in three iterations,
collocation defect 4.32e-10, independent 16-Z pressure holdout 7.05e-9.
Replay with `load_extended_profile()`, whose default requires passed pressure
holdout. The 65-node attempt failed off-grid at 0.00189012 and is retained
as `_coarse65.json`. Rebuilt-core angular closure remains 2.01e-11 over
33 Z points. Next bind actual axial/quadratic repair to this profile,
ensure smooth dependence on Z, and recover the radial component.

The actual heat collar pressure/shear/ODE stress implementation now passes
finite local cone and boundary-limit checks. The actual extended swirl in
`lei_ren_part1_extended_swirl.py` closes the angular target at 33 Z points
with independent integral defect 2.01e-11 and sampled positive F/negative
F_R. This is not the final stress cone over the whole outer connection.
Its actual pressure differs from the old core; the shared-pressure rebuild
is recorded above. The axial quadratic repair primitive is also
implemented, with synthetic independent moment residual below 4.20e-11.

Read `LEI_REN_EXTENDED_EXTERIOR.md`. The 33-point necessary angular radius
threshold is 999.6394381, with endpoint threshold 499.6662749. Both are
conditional on the saved core, not sufficient for closure.

Actual source collar moment targets are now transported inward by
`lei_ren_part1_exterior_targets.py`. Use Rb=2048, ell=.5, Ra=1242.17479109:
ell=.75 would put the inner boundary below the necessary angular threshold.
Actual collar targets at Ra pass the angular and endpoint necessary screens
at all 33 Z points. Independent heat-transfer identities agree within
8.53e-14; radial moment-density derivatives agree to relative 1.01e-8.

The quadratic target is about +2.46463, so zero axial/mixed moments alone
do not restore terminal data. The target-only receipt intentionally leaves
P0 unset; use the rebuilt profile's actual axis pressure when computing its
cumulative pressure target. Next restore the axial, mixed and quadratic
moments on this same extended swirl, then its full stress cone.
NS recursion, whole-background cone and higher-order flatness remain open.

## Live route after actual moment/field diagnostics — 2026-09-29

Read `LEI_REN_BACKGROUND_FIELD.md` and its retained JSON receipts. The
callable 3D kinematic seed now has repaired axial/mixed moments, dependent
radial velocity, physical axial localization, full radial-tail energy and
six-scale geometry diagnostics. Finest sampled Cartesian divergence is
9.34e-7 (relative 3.09e-9); full sampled energy is 0.002711–0.002890.
These results do not establish NS recursion or stress admissibility.

The angular budget rejects the short Rjoin=0.2 seed as the next admissible
source route under its current axis/pressure choices. Necessary angular and
pressure inequalities fail at all nine Z points. Next implement an extended
outer/collar family from Part I Sections 4–7, screen those inequalities,
recompute common pressure/core, and restore all five moments before cone and
stress/remainder claims. Do not continue unconstrained fits of the short join.
The earlier short-join tasks below are historical dependencies/diagnostics.

# Current checkpoint - constrained velocity-field integration

## Lei–Ren Part I integration — 2026-09-29

Read [LEI_REN_PART_I_INTEGRATION.md](LEI_REN_PART_I_INTEGRATION.md) before
continuing the background/stress route. It pins arXiv:2609.35406v1 and
introduces LR1-01 through LR1-13 without changing the final `1e-3` gate or
promoting a candidate. Separate background residual, admissible stress
divergence, flat remainder, and the full corrected-field residual.

Implemented Part I stress tensor completion (including its theta-theta
entry), projected divergence, radial-remainder preservation, core axis
slopes and fixed-sector bounds in `lei_ren_part1.py`. The independent
manufactured Cartesian divergence check has maximum error 1.01e-11.
This is an algebraic implementation check, not an existing-profile certificate.

The mapping delta=2h means the current h=0.005 lies outside the paper's
stated delta<1/200 range. Do not silently transfer its theorem to this
seed or silently change old parameters. LR1-03 now has a declared source-compatible parameter manifest (h=0.001).
LR1-04 has a finite-interval five-moment/pressure API and a local
PaperCoreReference receipt (pressure identity defect below 4e-17); actual shared outer
pressure and renormalized tails remain open. Follow the four acceptance
stages in PROJECT_GOAL.md; full 1e-3 acceptance follows reliable corrections. Part I does not complete oscillatory
cancellation or global smooth forcing. Keep ST006 as the retained numerical
baseline under its own protocol; ST073 metrics are not directly comparable.

## Actual exterior-pressure binding and core handoff — 2026-09-29

The joined swirl now uses the existing exact heat exterior outside R=0.2.
Its axis pressure is computed from the same swirl integral and actual heat
tail. Pressure matching at R=0.2 is within 7e-18; an independent radial
pressure derivative check has max error 1.76e-9. This independent smooth
blend is not the paper inward collar or full five-moment restoration.

The nonlinear finite core can now accept this pressure on its Chebyshev
grid. `lei_ren_part1_pressure_core.py/.json` records a degree-four
core/pressure iteration, including independent off-grid defects and saved
axis pressure values. `load_core()` rebuilds the local core and joined swirl
without repeating the iteration. Inspect the actual JSON convergence and
holdout flags before accepting a pressure handoff. The retained 257-node
receipt converges in seven iterations: collocation defect 5.47e-10 and
independent profile-pressure defect 2.14e-8 (not momentum residual). These are local numerical
compatibility results, not admissibility, finite global energy or recursion.

Next: compute and restore all five moment functions on this common profile,
including renormalized heat tails; construct axial/radial matching and
profile-derived stress. Then evaluate geometry and separate R_B, D T_B, E_B
across scales. Do not return to independent pressure fitting.

## Latest continuation checkpoint - 2026-09-27

Balanced-candidate improvements persist at all three tested scales. Physical
L2 at scales 1, 1/2, 1/4 is 2.6796083e6, 4.5122969e6, 7.5985389e6;
pulled-back L2 is 2.6796083e6, 2.6783825e6, 2.6771946e6. Against the
previous parent, L2 improves about 3.68% and maximum about 12.9% at each
scale. Against itself across scales, normalized residual remains nearly
constant: this is a better reference profile, not dynamical closure.
Sources: `scale_reference_refined_balanced_multiscale.json` and
`scale_reference_refined_balanced_rescaled_defect.json`.

The wave fit's quartic has a real stationary amplitude increment near
-0.404914, outside the earlier 1% trust bound. A separate larger-amplitude
experiment is now exploring wave amplitude [0.1, 1.1], holding background
and local corrections fixed. It must not be confused with the bounded
candidate or a uniform rescaling of the entire field; morphology, physical
constraints and independent residual checks will be required before use.

Actual-field replay of the corrected balanced candidate is COMPLETE:
44,400-point L2 2.67960826077e6 and maximum 2.41717297546e11, agreeing
with the analytic correction prediction. Both improve relative to the
five-percent parent (2.78208007209e6, 2.77457969824e11), and to the
original pressure-projected reference (2.93063976034e6, 2.70437112609e11).
The replay caches all 58 chunks for subsequent scale diagnostics.
Source: `scale_reference_refined_balanced_actual.json`.

The wave-amplitude alternative also improves on a newly generated 4,800
point support-split grid with angular shift 0.231: L2 ratio 0.977339,
maximum ratio 0.958547 versus the same parent evaluated on that grid.
This is fresh-grid evidence, not spatial convergence or a global PDE
certificate. The test spans the same full wave cylinder; no domain was
removed. Source: `scale_reference_wave_holdout_o4.json`.

The refined fitter's imaginary-sign bug is now fixed and both fits rerun.
Matrix-vs-basis and loader-vs-basis comparisons agree to 1.18e-15 and
5.00e-16 relative. Invalid reports are preserved with `_invalid_imag_sign`
suffixes. The peak-heavy direction still fails, but the corrected balanced
direction improves refined-grid L2 from 2.7820801e6 to 2.6796083e6 and
maximum from 2.77458e11 to approximately 2.41717e11. Additional velocity
L2 is 0.0021639143; cumulative correction L2 on the refined grid is
0.017183348 (7.941% of the recovered original reference norm on that grid).

Cross-grid transfer also improves both metrics: on the original 18,720
points, L2 changes 2.3750216e6 -> 2.3576868e6 (-0.730%), maximum
1.5730462e11 -> 1.5597136e11 (-0.848%). This grid was used for the parent
fit, so it is not a wholly unseen holdout. Actual-field replay of the new
candidate is in progress; no PDE or dynamical recursion acceptance yet.
Sources: `scale_reference_refined_balanced.json`,
`scale_reference_refined_balanced_transfer.json`.

Historical invalidation: review found an imaginary-sign error in the new refined fitter's
`_basis_delta` helper: it used +Im where the shared real-control basis uses
-Im. This invalidates that fitter's initial quadratic line-search numbers
and the first balanced follow-up's cumulative correction norm. The prior
refined-fit rejection paragraph below is historical, not authoritative
until the corrected rerun is recorded. The actual 44,400-point replay,
earlier assembly check and cross-grid helper use separate correct paths.
Do not interpret the erroneous 59% cumulative correction estimate.

The independently implemented original-wave amplitude fit contracts the
complex wave coefficients directly and is unaffected by that helper.
It gives refined-grid L2 2.713699e6 and maximum 2.684803e11, reductions
2.458% and 3.236% from the current corrected reference, with amplitude
increment -0.00637797 and an additional velocity norm 0.0021639143.
Pressure increments are added to the current corrected pressure once.
This remains a fitted diagnostic pending transfer and actual-field replay;
it is not an accepted dynamical recursion step.

The first refined-grid peak-weighted trust fit is complete and rejected.
All 13 exact-quadratic backtracking trials worsen physical-volume L2.
The smallest step (1/4096) gives L2 2.7824619e6 versus 2.7820801e6,
while peak decreases slightly to 2.7742807e11. The best sampled peak occurs
at 1/128 but raises L2 to 3.14802e6. `scale_reference_refined_fit.json`
therefore has `selected: null`. This rejects this weighted direction, not
the entire correction space. The next solves vary peak weights and test
an independent original-wave amplitude direction; both retain the joint
L2/maximum acceptance rule. Once a grid is used for fitting it is no longer
an independent holdout for that candidate.

A callable compact outer moment compensator is now available in
`scale_reference_outer_replay.load_compensated_reference()`, but is NOT
adopted. Its C-infinity axisymmetric swirl leaves the inner wave untouched.
Order-32 calibration and independent order-64 bump integration agree to
2.38e-8 relative; the reference moment still carries a 0.259% coarse/fine
quadrature discrepancy, so global neutrality is not certified.

The annular momentum replay rejects direct addition under the existing
scale generator, unchanged pressure and zero force: at order 20, outer L2
increases from 9.39788 to 2352.97 and maximum from 1.19186e5 to 2.03104e7.
Order 12 gives the same adverse direction but differs by about 10% in
corrected L2; these are not converged continuum norms. The correction's
self energy is 8.41143e-5; including the measured cross term, the total
energy change is 8.37790e-5 on this annulus. Moment compensation therefore
needs its own dynamical matching; it cannot be treated as a free repair.
Sources: `scale_reference_outer_moment.json`, `scale_reference_outer_replay.json`.

Pulled back to the reference coordinates with the isotropic acceleration
factor s^(3/2), the multiscale L2 values are 2.7820801e6, 2.7808352e6,
2.7796287e6; maxima are 2.7745797e11, 2.7738288e11, 2.7730838e11.
The scale-1/4 ratios are 0.999119 and 0.999461. These small sampled
decreases do not establish a uniform contraction sufficient for closure.
The anisotropic component normalization is bounded separately rather than
inferred exactly from aggregate scalar norms. See
`scale_reference_rescaled_defect.json`.

The corrected reference's refined-grid angular budget still assigns
50.164% of squared residual to the pressure-invariant mean azimuthal
component (L2 1.97045e6). Its convection is almost entirely fluctuation
momentum flux (1.96934e6), versus 1727.36 for convection of the angular
mean. This supports continuing velocity/wave-flux correction rather than
pressure-only fitting. See `scale_reference_refined_budget.json`; these
are sampled norms, not continuum lower bounds.

The assembly check is complete: rebuilding all 360 response columns gives
the stored 18,720-point fit L2 exactly, with maximum relative difference
1.94e-16. Eight refined-grid probes, including both peaks, agree with the
callable field jets; pressure-gradient relative difference is 2.54e-10.
No assembly mismatch was detected. The next optimization uses the refined
44,400-point jets and requires both sampled L2 and peak improvement.

The first cached multiscale replay gives (L2, maximum):
scale 1: (2.78208e6, 2.77458e11); scale 1/2: (4.68490e6, 7.84557e11);
scale 1/4: (7.88927e6, 2.21847e12). These use p_s=s^-1 p_0 and zero force.
Absolute growth alone does not quantify relative imbalance in a singular
scaling; pulled-back residual contraction must also be checked. None meets
the absolute target, and no dynamical recursion step is accepted.

Global axial angular momentum of this candidate is approximately
2.81266605e-8 (order 4 to 8 change 0.259%). The prescribed whole-field map
has M_z(s)=s^1.49 M_z(1), requiring nonzero net torque. An outer moment
compensation is being constructed as a necessary compatibility correction;
it does not establish local momentum balance or smooth forcing.
Sources: `scale_reference_assembly_check.json`,
`scale_reference_multiscale.json`, `scale_reference_angular_momentum.json`.

Corrected independent replay of the five-percent trust continuation is now
complete on 44,400 points, with pressure counted once: L2 2.78208007209e6
versus baseline 2.93063976034e6 (-5.069%), but maximum 2.77457969824e11
versus 2.70437112609e11 (+2.596%). Do not adopt this candidate as a joint
momentum improvement. The much better fit-grid numbers below do not carry
over directly. The completed assembly check supports a grid-dependent fit
effect rather than a detected implementation mismatch.
See `scale_reference_trust_nonlinear_actual.json`. Reusable per-chunk jets
are now saved during replay, and failure reporting preserves prior results.

The repeated trust-region continuation now reaches fit-grid L2 2.37502e6
and maximum 1.57305e11, about 18.87% and 40.68% below the common
pressure-projected scale-generator baseline, at cumulative 5% velocity
change. Source: `scale_reference_trust_nonlinear_fit.json`. Independent
momentum and physical gates remain pending. The solver now exposes explicit
warm-start, solver and output arguments, preserving earlier frozen reports.
The completed shape replay of this candidate retains radial contraction,
relative axial elongation, stronger winding and decreasing sampled energy
at prescribed scales 1, 1/2 and 1/4. This is kinematic evidence only.

The velocity-metric trust-region solve improves the one-percent update:
`scale_reference_trust_fit.json` gives fit-grid L2 2.69257e6 and maximum
1.76492e11, down about 8.02% and 33.44% from the pressure-projected baseline.
This is a reduced-rank constrained fit, not an independent PDE result.
The first independent replay contained a duplicated prior pressure
projection in its corrected case. That corrected result must not be used;
the uncorrected baseline case is separately reusable. The reference loader's
`base_reference` pressure was also fixed to contain only the prior pressure
projection, excluding the selected step's pressure increment. Numerical
shape observations using velocity alone are unaffected by this pressure fix.

Five bounded nonlinear reference updates reach fit-grid L2 2.77569e6 and
maximum 2.48638e11, reductions of about 5.18% and 6.24% against the original
pressure-projected scale-generator baseline. The cumulative velocity change
is 5%, versus 1% for the alternatives below. This candidate is in
`scale_reference_nonlinear_fit.json` and loads through the reference adapter;
independent momentum and physical gates remain unverified. These are
optimizer iterations, not accepted dynamical recursion steps.

A bounded original-wave amplitude plus pressure fit now improves the
scale-generator residual by 2.011% in L2 and 2.193% in maximum on the fit
grid, under the same 1% velocity-change budget. The mean is held fixed.
Use `scale_wave_reference.load_reference()` for this alternative; its
independent momentum and physical gates remain pending. Its fit-grid L2
2.86858e6 remains far above the requested 1e-3.

The first actual reference-velocity correction has been fitted with a 1%
velocity L2 trust bound. Its exact-quadratic fit-grid residual improves from
2.92745e6 to 2.89622e6 L2 and from 2.65177e11 to 2.61537e11 maximum.
Independent momentum and physical gates remain pending. The corrected
reference's prescribed scale family shows radial contraction, relative
axial elongation, stronger winding and decreasing sampled energy at scales
1, 1/2 and 1/4; these are kinematic morphology observations, not NS recursion.
See `scale_reference_velocity_step.json`, `scale_reference_candidate.py`
and `scale_reference_shape_replay.json` in `experiments/root_st073`.

Scale-compatible dynamics now has a direct diagnostic: replacing only the
fitted time derivative with the solenoidal scale generator increases inner
momentum L2 from 1.68732e6 to 2.92940e6 (+73.61%), maximum from 8.91218e10
to 2.66874e11. Pressure-only projection with the existing 90 local pressure
columns reduces the new L2 by just 0.0664%. Therefore the next correction
must address the reference velocity, not merely its independently fitted
time tangent. A cylindrical budget localizes 94.72% of squared residual
inside the two existing patches. About 51.18% is mode-0 azimuthal error,
whose true angular mean cannot be changed by a periodic pressure gradient.
See [dynamic matching](ST073_SCALE_DYNAMIC_MATCHING.md).
The prescribed scale map also admits a uniform energy bound tending to
zero as scale tends to zero; this analytic kinematic bound does not establish
NS dynamics or carry over automatically to new corrected reference fields.

New profile-aware continuation `localized_drift_1800_fit.json` reaches drift
1799.999991729 on the fitted interval, about 9.98% below the original parent.
The assembled refined-grid momentum maximum/L2 are 8.53686e10 / 1.67258e6;
all 81 sampled cones and endpoint shape gates pass. Independent 18,720-point
reference replay instead gives maximum/L2 8.91218e10 / 1.68732e6: L2 improves
over the first patch, but its maximum worsens by about 4.3%. Do not adopt it
as a joint momentum improvement. Completed endpoint replay also shows peak
regression: maximum/L2 9.01431e10 / 1.68797e6 (peak about 5.4% above the
first patch). Both rows are frozen in `localized_drift_1800_fit_actual_replay.json`.
Direct global drift is 1799.999991721 at
delta_k=1e-6, 1800.000663267 at 2e-6, and 1800.006038192 at 1e-5; the fitted
cap does not extend uniformly to these longer intervals.
The fit trades momentum quality for reduced profile drift; it is not an NS
step. The velocity scale map now has an implemented physical-time generator
in `scale_transport_generator.py`, checked against shrinking-time differences.
The next decisive calculation replaces the fitted time derivative with this
generator and measures the momentum defect with explicit pressure and zero
forcing. No recursive dynamical closure has been established.

The stricter `localized_parent_drift_fit.json` is now also independently
replayed: maximum/L2 8.21304e10 / 1.65408e6 at reference and 8.11335e10 /
1.65466e6 at the short endpoint. Its directly evaluated profile drift is
1999.619777783, restoring the original parent's level. It trades some
momentum reduction for this stricter profile cap; retain both candidates
below as documented alternatives. `global_two_patch_candidate.load(source)`
can load this stricter source explicitly. A velocity-only solenoidal scale
map is also available; see [scale transport](ST073_SOLENOIDAL_SCALE_TRANSPORT.md).
That map and its time derivative are kinematic, with no NS step accepted.

Latest preferred momentum/shape continuation: the drift-capped two-patch
candidate in `localized_two_patch_constrained.json`. Independent actual
18,720-point replay gives maximum/L2 7.49634e10 / 1.63977e6 at reference
and 7.40549e10 / 1.64031e6 at k0+1e-6, improving both metrics over the
one-patch baseline at both times. `global_two_patch_candidate.load()` is
the assembled global field interface. Direct profile drift is about
2122.104, matching the imposed cap but still above the unpatched parent's
1999.620. The stricter parent-level and 1800-drift fits above are now available
as alternatives. Full-support residual, critical-time behavior and scale
recursion remain unproved; zero recursion steps are accepted.

The preceding continuation baseline was the constrained localized patch
in `experiments/root_st073/localized_constrained_tangent.json`, with its
fixed enriched 264-control parent. Independent actual Cartesian momentum
replay on 18,720 points confirms simultaneous maximum/L2 improvement at
reference and k0+1e-6: corrected maxima 8.54481e10 / 8.54896e10 and volume
L2 1.81103e6 / 1.81146e6 (about 27.7% / 6.1% below the respective parent).
The physical time interval is only 1.6919e-10; no NS step or recursion is
accepted, and the residual target 1e-3 remains unmet by many orders.

`global_localized_candidate.load()` supplies the same inner candidate with
the global compact mean and disjoint exterior collar. Assembly probes
preserve inner/exterior values and vanish beyond the support union at both
tested times. This is a global diagnostic field, not a full-support residual
certificate or a critical-time extension. The older paragraph below about
missing axial localization describes the pre-extension base, not this field.

Prior acceleration fitting had a weighted-objective gradient bug, now fixed.
The corrected candidate improves independent L2 but increases the holdout
peak, so it is not adopted as a joint improvement. Prior unconstrained
acceleration also failed the aspect gate. Prioritize joint optimization of
the two compact spatial patches, retaining moments, positive cone margins,
endpoint shape constraints and a sampled peak cap. See
[current experiment evidence and next solve](ST073_LOCALIZED_CONSTRAINED_NEXT.md).

## Earlier ST073 research update - 2026-09-27
A globally localized diagnostic now combines the balanced compact wave,
inner acceleration correction and exterior collar tangent through one
callable field. Sampled assembly checks at reference and endpoint preserve
the constituent inner/outer values exactly and return finite axis values.
This does not establish whole-domain momentum or recursion. See
[unified global candidate](ST073_GLOBAL_ACCELERATION_CANDIDATE.md).

The inner acceleration diagnostic reduces sampled endpoint L2 by 4.78%; its
actual shape replay fails relative axial elongation, so constrained acceleration fitting is required before adoption. Independently, exterior collar
L2 falls 7.97% on a separate replay grid. These percentages must not be added.

Global-scope audit confirms a separate required gap: the current base is
registered only on abs(eta)<=0.5, and its radial heat exterior has no axial
localization. The candidate is not yet a globally defined finite-energy
field. A divergence-preserving axial exterior/pressure construction is
being investigated alongside local residual correction. See
[global exterior gap](ST073_GLOBAL_EXTERIOR_GAP.md).

Degree-3 mean enrichment with direct endpoint geometry now reduces
independent actual momentum L2 to 1.92894e6 and maximum to 1.18244e11
(about 27% and 30% below the prior 236-control candidate). All three
actual sampled forward geometric signs pass. Enlarged-mean independent
replay gives moment maximum 2.1084e-7 and 27/27 cones passing, minimum
margin 9.94824e-5. Higher-quadrature refitting now lowers same-grid L2 from 1.92798e6 to 1.89886e6, but increases the peak from 1.13931e11 to 1.15611e11. The peak-capped solve improves independent L2 to 1.89999e6 but worsens its maximum to 1.20729e11. It is not a joint improvement; peak constraints must cover both grids before a separate disjoint check. See [refined momentum](ST073_REFINED_MOMENTUM.md).
No NS step or scale recursion is accepted. See
[enriched mean endpoint](ST073_ENRICHED_MEAN_ENDPOINT.md).

Degree-3 mode-2 enrichment now combines with shape constraints: independent
actual-field momentum L2 is 2.63329e6 and max 1.68320e11, while all three
sampled forward geometric directions remain favorable over delta-k=1e-6.
Moments/cones pass the assembled solve. Mean harmonic now contributes
80.33% of squared residual and is the next enrichment target. This remains
an affine candidate, not a momentum-accepted trajectory or recursion. See
[enriched shape tangent](ST073_ENRICHED_SHAPE_TANGENT.md).

Larger shape-direction margin now gives all three requested sampled forward
changes over delta-k=1e-6: radial RMS decreases, aspect increases, and
weighted angular speed increases. Peak swirl still decreases. This is an
affine candidate increment over physical delta-t about 1.692e-10, not an
accepted NS step or scale recursion. Independent reference-time replay now
gives momentum L2 2.84751e6, max 2.10261e11, moment max 2.1176e-7 and all
27/27 cones passing. Endpoint momentum/matching remain unaccepted. See the updated
[shape-constrained result](ST073_SHAPE_CONSTRAINED_TANGENT.md).

Shape-constrained tangent now preserves independent moments and 27/27
cones while giving all three requested central derivative signs. However,
the forward delta-k=1e-6 interval still expands radius and lowers aspect:
curvature overwhelms the small direction margin. A stronger-margin
candidate is generated and under actual-field replay. No geometric step
or scale recursion is accepted. See
[shape-constrained tangent](ST073_SHAPE_CONSTRAINED_TANGENT.md).

Separately, spatial degree-3 mode-2 enrichment reduces independent-grid
L2 from 2.8472e6 to 2.6329e6; degree 4 is slightly worse on that grid.
This remains to be combined with shape constraints. See
[higher harmonic tangent](ST073_HIGHER_HARMONIC_TANGENT.md).

Target-direction check now rejects advancing the cone-compatible tangent
as successful recursion: fixed-cylinder radial RMS increases and axial/
radial aspect decreases, although weighted angular speed increases.
These are sampled affine-tangent trends, not an NS trajectory. Shape-rate
constraints must join the next tangent solve; residual and matching alone
are insufficient. See [wave shape direction](ST073_WAVE_SHAPE_DIRECTION.md).

Joint moment/cone correction now passes independent actual-field cone
replay at all 27/27 locations (81/81 inequalities, minimum margin 9.9759e-5).
Order-96 integral moment maximum is 2.1219e-7. Independent full momentum
L2 is 2.8472e6 and maximum 2.1019e11: restoring cones costs 6.98% in L2
against the moment-only candidate, retaining about 61.4% reduction versus
the old dense tangent baseline. No recursive step is accepted. The next
construction varies wave shape with both moments and cones in the inner
solve. See [joint moment/cone result](ST073_MOMENT_CONE_WAVE.md).

Moment-constrained wave correction now restores actual order-96 integrated
moments to max 2.54e-8, versus 583600.74 for the earlier candidate, while
independent full momentum L2 changes only from 2.6453e6 to 2.6614e6.
Cone replay improves to 12/27 locations passing but remains a failure.
Higher-grid wave growth is +997.0958. No recursive step is accepted.
See [moment-constrained wave result](ST073_MOMENT_CONSTRAINED_WAVE.md).

Momentum-aware wave-shape optimization now reduces independent actual
patch volume L2 from 7.3766e6 to 2.6453e6 and maximum from 5.6409e11 to
2.0098e11, about 64% for both. Five-node flux constraints remain matched.
The growth sign discrepancy is traced to coarse quadrature of nearly
cancelling production/dissipation terms. Higher z20/r11 quadrature gives
lambda +997.096 for the new wave; positive growth has this finite-grid
support, but the 1000 floor and continuum convergence are not verified.
Renewed mean compatibility remains required. See
[momentum-aware wave construction](ST073_MOMENTUM_AWARE_WAVE.md).

Denser harmonic fitting now improves the new order-13 independent physical
replay: frozen-wave patch L2 9.9155e6 becomes 7.3766e6, with maximum
5.6409e11. The earlier severe fitting instability is reduced, but the
absolute residual remains unacceptable. The next construction optimizes
wave shape against full momentum while retaining local flux/growth targets.
See the denser correction in [joint wave/mean status](ST073_JOINT_WAVE_MEAN_STATUS.md).

The old growing-wave joint 45-control mean fit has completed spatial replay:
annulus volume L2 worsens from 193257.9992 to 193600.2673, despite a small
patch improvement. Do not extend this candidate as a successful recursive
step. See [joint mean/wave spatial result](ST073_JOINT_WAVE_MEAN_STATUS.md).
Full harmonic tangent replay on the locked co-designed candidate now fails
the spatial holdout: volume L2 is 1.1175e8 versus frozen-wave 1.0036e7 and
mean-only 1.2146e5. The independent harmonic budget identifies the second
harmonic as 98.52% of the fitted squared residual. Denser spatial fitting
is next; no time step or scale recursion is accepted.

Joint polarization/growth fitting now matches the required radial flux at
five nodes with relative error 4.69e-11. Higher-quadrature growth lambda is
997.36, positive but below the imposed 1000 floor. Wave RMS is 4.347 times
mean RMS, so this remains a large nonlinear candidate requiring complete
residual correction. See [stress and growth co-design](ST073_STRESS_GROWTH_CODESIGN.md).
No accepted physical trajectory or scale recursion follows from this fit.

The actual physical mean-plus-wave replay is now available in
`mean_wave_replay.py/.json`. The old scalar-amplitude wave candidate fails
to improve the complete residual: independent patch volume L2 increases
from 120498.06 to 120955.20 after fitting oscillatory time/pressure controls.
This rejects a conclusion based on its small angular-mean improvement.
See [physical mean-wave replay](ST073_PHYSICAL_MEAN_WAVE_REPLAY.md).

The constrained short step retains a positive mode-1 instantaneous energy
rate on the fixed patch: split moderate/higher quadrature gives 13161 and
14501, versus initial 19027 and 20438. The 9.24% endpoint sensitivity is
not convergence evidence. The two unconstrained delta-k=0.01 paths instead
lose growth in all eight tested modes. See
[the evolved growth screen](ST073_EVOLVED_WAVE_GROWTH.md). Joint stress and
growth co-design is now the next construction; energy positivity alone is
not covariance matching or actual wave integration.

The first dynamically constrained state step is complete at delta k=0.0001.
Recomputed endpoint compatibility passes 27/27 sampled cones, with integrated
moment maximum 8.88285e-8; endpoint momentum maximum/L2 are 4.43477e9 and
188215.12. This is initial/endpoint local-tangent compatibility, not a
continuous-time or recursive certificate. Fixed-cylinder spin does not
increase. See `meridional_constrained_evolution.json` and the
[state evolution record](ST073_STATE_EVOLUTION.md). Next investigate actual
wave growth/coupling rather than extending mean relaxation alone.

State-dependent unconstrained evolution now improves on frozen slopes over
delta k=0.01: the independent endpoint maximum is 2.27427e9 versus
6.00640e10, and volume L2 is 84464.34 versus 1473368.61. This uses a
165-point fit and two explicit Euler steps, with the full quadratic
nonlinearity. It is not a matched or converged trajectory. A coarse
44-point fit gave misleadingly poor evolution because its axial sampling
aliased basis directions; the dense unregularized matrix is full rank.
See [state evolution and its limitations](ST073_STATE_EVOLUTION.md).

Important target check: on a fixed physical cylinder, that refreshed
delta-k=0.01 path broadens the enstrophy radial RMS from 0.0005595 to
0.0016775, lowers the axial/radial ratio from 0.9337 to 0.3904, and lowers
weighted angular speed from 2.125e6 to 0.968e6. The residual improvement is
not evidence of the requested vortex amplification. The diagnostic is
`vortex_state_observables.json`; it describes the sampled cylinder, not an
identified global core. Do not extend this unconstrained mean-relaxation
path as if it were successful scale recursion.

The constrained 44-direction replay is now complete. On the independent
176-point spatial holdout, momentum maximum decreases from 7.01009e9 to
4.50644e9 and physical-volume L2 from 237497.46 to 188297.75. All 27 sampled
stress cones pass at radial quadrature order 64; the four integrated moments
have maximum absolute value 1.47389e-5 at order 96. The selected solver point
is feasible and improves the objective, but SLSQP returned status 8, not
convergence. See `broad_meridional_constrained_replay.json`. This supersedes
the pending compatibility statement below, not the full PDE/recursion gates.

An opt-in grouped evaluator preserves the original fields and tested Cartesian
jets while reducing the cold 2000-point field benchmark from 5.2885 s to
0.1337 s (39.54 times). The full constrained replay took 121.88 s. This is
computational acceleration, not a reduction of the governing equations.
The next experiment evolves the velocity coefficients and recomputes their
slopes at each state; frozen initial slopes are already rejected below.

For that experiment, `integrated_state_moments.py` now assembles current-state
moment equations directly from fixed-radius conservation integrals. At a
nonzero velocity state, its four predicted moments agree with independent
actual-field integrated replay to 3.59e-8 (order 96). This removes the need
to use the less accurate order-24 pointwise-residual integral rows as the
evolution constraints. It does not establish continuum derivative accuracy.

The first full-vector meridional/pressure correction now reduces both metrics
on the same independent spatial holdout: maximum 7.01009e9 -> 4.38510e9,
physical-volume L2 237497.46 -> 186820.70. Its 44 directions preserve the
instantaneous broad-shear velocity and fit its time derivative plus pressure.
This is an unconstrained instantaneous experiment: the earlier moment/cone
passes do not transfer, and finite-difference divergence remains 0.154809.
The next solve must retain the new full-vector objective while reinstating
moment/cone compatibility. No scale recursion or PDE acceptance is established.
See [the meridional correction record](ST073_MERIDIONAL_MOMENTUM.md).

Off-time replay of that unconstrained affine-k candidate preserves its local
benefit through delta k=0.001, but direct extension to delta k=0.1 and 1
fails badly (momentum maxima 1.52546e12 and 3.82489e14). These are callable
evaluations, not time integration. The combined broad-amplitude slope crosses
zero at delta k=0.034086. Do not use frozen initial slopes as a recursive
trajectory; the next evolution needs state-dependent updates and full nonlinear
terms. Full records are `broad_meridional_time_audit.json` and
`broad_meridional_scale_audit.json`. A separate constrained 44-direction solve
is in progress; no compatibility result is claimed until its replay completes.

The next state-dependent evolution now has a factorized nonlinear residual
engine, `meridional_state_cache.py`. It supports 25 velocity values, their
25 time slopes and 19 pressure coefficients, retaining all quadratic
advection. Two nonzero-state five-point comparisons agree relatively to
2.25e-11, but absolute discrepancies reach 0.02481; this is not final-gate
accuracy or an integrated trajectory. See
[the state-cache record](ST073_MERIDIONAL_STATE_CACHE.md).

A new broad annular r^-2 swirl direction admits positive instantaneous
wave energy growth: at amplitude 44.8774 the mode-1 rate is 1.0217e4
on the refined split quadrature, with all 22 sampled centrifugal growth
conditions positive. This changes the earlier all-decaying energy result
by changing the mean, not by changing viscosity or numerical time steps.
Its added local swirl RMS is about 50 times the old mean velocity RMS.
It is not an accepted low-residual replacement, and modes 2 through 8
still decay. The mean compatibility result is recorded in
`broad_shear_dynamic_control.json`: adding da/dk=-227.669 fixes the wave
region's stress sign. Independent replay gives moment maximum 1.45e-5,
6/6 outer cone passes and 5/5 wave-region cone passes. The instantaneous
growth matrix is preserved, but full sampled momentum remains 7.46e9.
The earlier 18-control seed and its 0/5 local cone result are historical;
see [the original co-design record](ST073_SHEAR_CODESIGN.md).

Radial pressure primitives are now implemented. At the old centrifugal
peak an inner-datum repair lowers the norm from 7.48e9 to 1.55e9, but
creates exterior axial imbalance. A compact primitive restores exterior
pressure while retaining a 4.71e9 collar defect. The numerical radial
integral budget shows that pure swirl plus compact pressure cannot remove
this added radial imbalance. The next construction needs a divergence-free
meridional time correction and/or actual wave stresses, coupled to the
pressure source. See [dynamic and pressure matching](ST073_BROAD_DYNAMIC_MATCHING.md).

Full spatial potential evolution now includes the angular mean and modes
1 through 8, cutoff derivatives, nonlinear transport and viscosity.
An implicit weak-diffusion BDF solve fixes the explicit-step instability.
Three Cartesian time replays have momentum maxima 7.80e9, 3.79e9 and
1.91e9, still far worse than the background-only field. The oscillatory
energy proxy falls to 1.44% of its initial value; initial center covariance
is not retained. This is numerical integration progress, not successful
scale recursion or recursive amplification. See
[the full spatial evolution and decay audit](ST073_SPATIAL_FOURIER_EVOLUTION.md).
The initial-time energy screen now finds no positive growth direction
in any retained mode at support multipliers 1, 2, 4 or 8. At multiplier
8 the best combined amplitude rate is still -2.63e5, while the largest
strain-only rate is 1.42e4. See [the energy budget](ST073_WAVE_ENERGY_BUDGET.md).
The broad-shear seed above passes the initial energy-growth screen; the
remaining blocker is making that growth compatible with the full mean
and wave dynamics and the required two-component stress.
The paper requires growth followed by decay; decay at a pulse tail alone
is not a failure criterion.

A moving-normal principal amplitude/pressure inverse is implemented and
passes manufactured checks, but its current center-path spatial wave is
rejected. Full momentum on spatial holdouts grows from 5.91e9 for the
moving wave to 8.86e9-1.91e10 after its forced correction. The background
alone is about 4.50e5 on the same local grid. This trial does not replace
the feedback mean trajectory. See [the rejected wave experiment](ST073_MOVING_NORMAL_WAVE.md).

The new finite-dimensional spatial solve addresses these omitted operators
but does not provide the paper's supported inverse or recursive estimates.
Do not repeat center-only inverse sweeps or interpret a principal ODE check
as a full-field residual reduction.

State-dependent pressure/swirl-slope integration now completes the short
interval k in [11,11.001]. Three direct time replays have integrated moment
maxima 9.87e-4, 4.99e-4 and 9.57e-4, with 6/6 outer cone samples passing
each time. At k=11.0005, the matched fixed-slope comparator is 1.658,
so feedback reduces that moment drift by about 3,322 times. See
[the completed trajectory diagnostic](ST073_STATE_DEPENDENT_EVOLUTION.md).

This is sampled short-time compatibility, not scale recursion. The worst
moment replay is close to 1e-3 without a verified error margin, spatial
neighborhood holdouts have not been replayed along this trajectory, and
full sampled momentum remains 1.42e6. Full momentum max/volume-L2, finite
energy, forcing/domain requirements and recursive contraction remain
unestablished. No candidate is accepted.

The preceding spatial seed passes eight neighborhood holdouts and nine
inner nodes at k=11. Shared radial integration makes the 22-node benchmark
2.53 times faster. The code and all completed trajectory/comparator
artifacts are in experiments/root_st073.

The snapshot below is historical, not current ST073 acceptance.

Snapshot date: **2026-09-20**.

This integration-branch snapshot predates the 2026-09-22 `main` research
consolidation. On `main`, ST006 is the historical runnable API baseline,
while ST061-D/P are newer residual-oriented controls and ST063-G2R is the
latest documented geometry experiment. The paired ST061/ST063 samples
report full-vector maxima about `0.0225`â€“`0.0272` and volume L2 about
`0.033`â€“`0.034`; their seeds and candidate artifacts differ from the
ST006 protocol below. None meets both `1e-3` momentum gates. See
`main:docs/FINAL_RESEARCH_SNAPSHOT_2026-09-22.md` and
`main:docs/RESEARCH_STATUS.md` for the original study records and
candidate availability. Do not rank them against ST073 annulus screens
without a shared complete-field validation protocol.

Machine-readable authority: [`project_status.json`](../project_status.json).
Active integration branch at snapshot: `codex/cr001-constraints@9491103442e53190b7d90d6f0fc5f3165aacc697`.

This file is intentionally a **current** checkpoint. Historical optimization experiments, rejected basis studies, exact-reconstruction notes, and earlier route-local diagnostics remain available in Git history and artifacts, but they are not merge blockers for the current delivery chain unless they affect the frozen candidate being shipped.

## Final delivery target

The nine lanes converge on a directly usable time-dependent Cartesian velocity field

`velocity(x, y, z, t) -> [u, v, w]`

that can be called from Python, saved/reloaded, exported on reproducible grids for MATLAB/Python, and inspected with streamline/vorticity diagnostics. Its publicly observable geometry and time evolution should be made progressively closer to the public OpenAI velocity-field visualization.

This target is **not** a paper-exact reconstruction target and **not** a complete blow-up proof target.

The three project states are independent:

- `velocity_export_ready`
- `visualization_ready`
- `pde_validated`

CI success, rendering success, optimizer convergence, or visual resemblance does not promote the other states automatically.

## Repository-wide PDE baseline

The retained same-protocol numerical baseline is ST006 on `main`:

- candidate SHA-256: `6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3`
- held-out seed: `9172801`
- held-out Cartesian points: `4096`
- validation times: six
- finest spatial step: `0.005`
- momentum sampled max: `0.1082289305112118`
- momentum volume-L2: `0.10758432876230622`
- registered momentum target: `1e-3`
- `pde_validated=false`

Only directly comparable complete-residual protocols may claim improvement over this baseline. Scoped Kokuno terms, morphology diagnostics, sampled derivative audits, or render metrics are not directly comparable unless they reproduce the same residual contract.

## Canonical constrained delivery

The current canonical constrained delivery remains `eq45_supported_velocity_candidate_v1`.

- public evaluator: `openai_ns_reconstruction.eq45_supported_delivery:velocity`
- save/load: `Eq45SupportedDeliveryField.save_candidate(...)` / `Eq45SupportedDeliveryField.load_candidate(...)`
- grid evaluator: `default_field().grid`
- delivery capsule: `artifacts/constrained/eq45_supported_delivery_capsule.json`
- supported child SHA-256: `2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d`
- parent SHA-256: `48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7`
- registered box: `[-2,2]^3`
- registered time interval: `[0.25,0.75]`
- physical support connection: `r < 2`, `|z| < 2`

Current canonical truth state:

```text
velocity_export_ready = true
visualization_ready   = false
pde_validated         = false
```

`velocity_export_ready=true` means the canonical Eq45 field is callable/saveable/exportable. It does **not** mean visual correspondence, exact OpenAI-field identity, PDE acceptance, paper exactness, or blow-up.

## ST052-M visualization-candidate delivery state

The frozen ST052-M child is no longer awaiting materialization. The following delivery pieces are already live on constrained ancestry:

1. whole-child exact-source-backed runtime and save/load identity;
2. checksum-bound `33^3 x 5 times` sampled-grid export;
3. NPZ and MATLAB-v5 MAT outputs;
4. GNU Octave MAT-load/render smoke;
5. live delivery-identity reconciliation.

Key merged integration PRs:

- #665 â€” ST052 grid export materialization;
- #680 â€” GNU Octave consumer/load/render smoke;
- #706 â€” live ST052-M grid delivery identity reconciliation;
- #777 â€” machine-state synchronization removing stale rematerialization routing.

The #706 repaired exact head `62ed0677675cbba3ff409a6ba124bc157bc5b504` completed repository `tests` workflow run `35475995041` with conclusion `success` before merge.

ST052-M still has a narrower delivery gap than the canonical Eq45 field: the repository has not closed or explicitly accepted the **standalone-package parent-runtime/dependency identity seam** as sufficient for a standalone ST052 export claim. Native MATLAB scientific execution has also not been performed for the ST052 grid; the current evidence is GNU Octave consumption.

Therefore ST052-specific state remains:

```text
source_runtime_backed_callable_save_load_ready = true
exact_source_runtime_identity_closed           = true
grid_export_materialized                       = true
octave_mat_load_render_smoke_passed            = true
standalone_package_parent_runtime_ready         = false
velocity_export_ready                           = false
visualization_ready                             = false
visual_correspondence_verified                  = false
pde_validated                                   = false
```

Do **not** open another ST052 whole-child materialization, NPZ/MAT exporter, or Octave replay lane. Those seams are already closed.

## Current CI truth

After #777 merged, the live constrained head is

`9491103442e53190b7d90d6f0fc5f3165aacc697`.

Its repository `tests` push workflow is run **35507408812**. At this snapshot it is still **queued**, so the live merge head must not be described as newly CI-green yet.

Runner backlog is not scientific evidence. A queued workflow is neither PASS nor FAIL.

## Active sibling delivery assets

### Agent 9 â€” source-observable and delivery diagnostics

`main` already contains ST054 delivery plumbing for:

- continuous Python callable;
- Python streamline/vorticity rendering;
- deterministic NPZ/MAT handoff;
- VTK/ParaView handoff.

Current open Agent-9 work includes:

- #825: forward-port the official-public qualitative observable contract to current `main`;
- #834: renderer-independent cylindrical morphology fingerprint for ST054.

At this snapshot their current exact-head workflows are queued. These are useful sibling assets and reusable diagnostics, but they do not silently replace the constrained Eq45/ST052 candidate identity.

### Agent 7 â€” morphology/capacity lane

Open Agent-7 work studies outer-reservoir morphology, temporal curvature, representation equivalence, and a compact toroidal swirl preflight. These are candidate-capacity or morphology experiments. They should not be promoted into canonical delivery merely because a local capacity gate passes.

New basis growth is not the shortest constrained-delivery path while the frozen ST052 runtime/diagnostic chain is still incomplete.

### Kokuno Agents 1â€“5

The Kokuno lanes are advancing an executable **strict-inner** PA.10 contraction-center stack: velocity, time derivative, spatial derivatives, self-advection, oscillatory composition, cylindrical mean projections, and independent derivative/advection audits.

They remain inner-only and do not yet provide all of:

- outer/global leading join;
- corrected/global candidate velocity;
- matched pressure;
- preregistered restricted forcing;
- complete `velocity/pressure/forcing` candidate API;
- same-protocol complete NS residual;
- final divergence-L2 / momentum max/L2 acceptance.

These lanes are valuable upstream scientific work but are not blockers for exporting a clearly labeled constrained visualization candidate.

## Next integration task â€” shortest delivery chain

The next constrained integration task is **not** another basis experiment.

For the same frozen ST052-M candidate:

1. close **or explicitly accept** the standalone-package parent-runtime/dependency seam;
2. preserve one unified `velocity(x,y,z,t)` identity through save/load;
3. replay the already-live `33^3 x 5` NPZ/MAT export;
4. run fixed-seed/fixed-camera streamline and vorticity diagnostics;
5. run source-observable diagnostics without inventing public numerical targets;
6. produce one Python/MATLAB-facing integration report.

If that candidate becomes a stable visualization candidate, it may be exported and labeled as such even while `pde_validated=false`.

Only if this exact frozen candidate is selected for PDE work should the integration route rebuild compatible pressure/restricted forcing and run a fresh independent 4096-point momentum/divergence validation before any `pde_validated` promotion.

## Duplicate-work firewall

Do not duplicate:

- ST052 whole-child materialization;
- ST052 NPZ/MAT grid export;
- ST052 GNU Octave MAT consumer smoke;
- ST054 Python render / NPZ-MAT / VTK plumbing;
- Agent-9 public-observable contract work;
- Kokuno strict-inner derivative/advection subterms already owned by Agents 1â€“5.

New work should close a currently open seam in the shortest `[u,v,w]` delivery chain or provide a genuinely independent validator for a selected candidate.

## Truth boundary

The following implications are forbidden:

```text
CI green              != PDE-valid
optimizer converged   != PDE-valid
export works          != visual correspondence
visual resemblance    != PDE-valid
visual resemblance    != exact OpenAI field
sampled local audit   != whole-domain validation
PDE-valid             != blow-up proof
```

The project should prefer a reproducible, clearly labeled callable visualization candidate over blocking all export on unresolved historical proof goals, while keeping PDE and exact-source claims fail-closed.
