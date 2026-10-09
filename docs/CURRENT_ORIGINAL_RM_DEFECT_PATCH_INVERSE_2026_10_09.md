# Actual Rm leading five-defect inverse and partial patch

Checked source [1d4903c6](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/1d4903c652476c44ac0d579b71c09e8971b42876). Full reconstruction **ACTIVE / INCOMPLETE**. Actual Rm histories now define the original five defects, feed the anisotropically normalized leading implicit inverse, and provide partial bump velocities, pressure and five cumulative primitives through ordinary Z order5 (Ur through4). Native source frames0,.5 only. The genuine leading repair has strict inclusion, contraction and a common Jacobian for higher coefficients. Whole-axis source functions, full radial mixed4, finite-N Rc/all24/controls/global N and real n-recursion remain open.

Predecessor: [actual reference/restoration through Rm](CURRENT_ORIGINAL_REFERENCE_RESTORE_FUNCTIONS_2026_10_09.md). Accepted caches and original bump/cutoff receipts are reused. No old pressure/core/repair constructor, ancestor producer or full checker is rerun. Root owns edits, computation and Git; the existing **GPT-5.6 Luna / max** worker provided read-only normalization and physical-prefactor reviews.

## Actual five-defect source adapter

The live upstream owner supplies `postrestore((-6,1))`. Its exact offset, six centered names, context, formal basis/ledger and Z0..5 rows are required. Family, implicit source, original analytic pressure datum and dependency hashes remain the same accepted objects. The five rows at Rm are

    d0 = mean_error
    d1 = mixed_error
    d2 = angular_error
    d3 = axial_square*invAm2 - swirl_error/2
    d4 = pressure_error/2
    invAm2 = exp(1.2)*(1+Z^2)^2*canonical_Pstar^-2
    Am = exp(-.6)*Pstar/(1+Z^2)
    Rm = 110*exp(10(logCstar+logPstar)-6).

These correspond to the original Delta_z/Rm, centered mixed angular moment/(sqrt2*Rm^1.5*Am), angular defect/(sqrt2*Rm^1.5*Am), centered axial/swirl defect/(Rm*Am^2), and Delta_p/Am^2. Defining reference primitives recover these expressions independently. Original analytic P0 is separate and never enters d4. All actual reference/restoration tails remain signed, nonzero formal factors.

## Canonical anisotropic normalization before materialization

The native Pstar scale is astronomical. Materializing angular defects before dividing by their original tiny angular unit would erase their relative scale behind a directed exponential tail bound. The code instead divides each actual factored source row in its common basis, then materializes its normalized enclosure.

    t = exact fixed positive axial box coordinate from accepted original repair
    angularUnit = 10^10*t^2*canonical_Pstar^-2
    S = diag(t,t,angularUnit,angularUnit,angularUnit)
    h = S*u; normalized defects = S^-1*d.

t is an arithmetic box coordinate, never a physical field value. The angular unit uses the same canonical Pstar factor as invAm2, rather than an unrelated inverse-amplitude cap. Both units are Z-independent. The original matrix has commuting2+3 blocks, so S^-1*L*S=L. The exact rewritten quadratic map uses

    fg' = angularUnit*fg; gg' = gg
    ff' = angularUnit*ff; ff_over_x' = angularUnit*ff_over_x
    invAm2' = exp(1.2)*(1+Z^2)^2/10^10.

The ordinary angular weight cover may have lower bound0, because it encloses a positive factor vastly below the materialization threshold. This is only an outer coefficient enclosure for nonlinear weights. It is NEVER a divisor, a chosen angular unit, or a replacement source value. Division uses the strictly positive formal angular unit. Uniform strict self-map/contraction on the outer coefficient box includes the actual positive-source map. The checker verifies canonical formal positivity and this distinction.

## Same unique leading implicit coefficient family

The solved equation is the original leading map

    Lh+Q(h,Am^-2)+d=0.

The normalized unit box has strict inclusion and contraction. Native initial contraction upper bounds: Z=0, 6.17698344e-06; Z=.5, 9.65153663e-06. The first derivative and every ordinary Z coefficient2..5 use the same invertible Jacobian. Controls map back by S, including all derivatives with no extra scale-derivative terms. Midpoints choose a fixed matrix preconditioner only; data, weights and coefficients remain enclosing intervals. Zero containment is a diagnostic, while closure follows from strict inclusion and uniqueness.

These are conditional local/source-frame solution jets. They do not certify an entire axial function domain. The separate finite-N Rc map is `B*h+N*r(N,Z)+Q(h,h)/N=0`; its actual complete targets and24-cell source transport are still required. This leading inverse is not substituted for those targets or controls. Broad existing global flags remain false.

## Partial bump fields and five cumulative primitives

The callable patch uses exact rational x=R/Rm in[1,e], the original three centers5/4,3/2,7/4 and radius.025. It integrates original beta partial weights using directed closed cells, or exact zero/full-support identities. With h=(c1,c2,xi1,xi2,xi3), f=sum xi*beta and g=c1*beta_1+c2*beta_3. The original five partial changes retain axial, angular, axial-angular cross and quadratic terms. They are added to the actual incoming defects before recovering mean/theta/mixed/centered energy/pressure.

    m=4Z+d0(x)/x; theta=5/8*x^1.6+d2(x)
    mixed=4Z*theta+d1(x)
    s=d3(x)-5/12*x^1.2; pm=d4(x)+5/2*x^.2
    V=4Z+g; Utheta=Am*(x^.1+f)
    Q=(2Z*V-Z*m*(1-delta)-(1-Z^2)*m_Z)/(1-delta*Z^2).

At x>=71/40 all bumps are complete and the identical full-weight implicit map proves the terminal five defects vanish. The code retains the actual inlet and unrefined diagnostic; this is an exact identity of the unique solution, not an inlet reset or a zero-containment closure claim. This identity is restricted to the two source frames and the leading map.

Physical recovery uses FIXED Rm prefactors:

    Mz=R*m
    Mtheta=sqrt2*Rm^1.5*Am*theta
    Mtheta_z=sqrt2*Rm^1.5*Am*mixed
    centered_Mztheta=Rm*Am^2*s
    Mztheta=centered_Mztheta+8Z*Mz-16Z^2*R
    Mp=Am^2*pm; P=Pstar^2*P0+Mp
    Ur=sqrt(R/2)*Q.

The Z-dependent Am and Am^2 products are differentiated before output, and the centered energy is uncentered before reporting physical Mztheta. Multiplying theta by the current R^1.5 instead of fixed Rm^1.5 would count its explicit x evolution twice. Signed tiny correction sectors remain separately accessible even when adding them to a large baseline widens the total enclosure. Utheta/Uz first ordinary y derivatives are available. Full radial mixed4 and radial-cell integration are not yet available.

API: instantiate `OriginalRmDefectPatchInverse()` (accepted receipt required), then `evaluate('0',(5,4))`, `evaluate('.5',(3,2))`, or `owner(label).coefficients()`. The public factory validates provenance; the internal fixture operator is private. Unadmitted source frames, non-rational/out-of-domain coordinates, invalid cells/box units/noncommuting matrices reject.

## Evidence and scope

Files: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_defect_patch_inverse.py, .json, _check.py, _check.json. Producer 15.953s; checker 77.484s; terminal exit0. **356 staged exact dependency hashes PASS**.

66 symbolic block/map/Jacobian/amplitude/defect/physical-unit identities pass. Two independent finite fixtures manufacture a smooth implicit solution from independently integrated genuine beta weights and evaluate direct partial integrals and physical primitives/velocities/pressure. 60 coefficient Taylor, 210 partial integral and 770 physical Taylor comparisons pass. Reference quadrature/differentiation is numerical diagnostic evidence, not a rigorous reference certificate. Its roundoff comparison budget is explicitly separate from production directed enclosures. Z-independent radial weights are cached before high-order reference differentiation.

Native checks compare 60 exact control records and 60 exact normalized-defect records, preserve 16 inlet row joins and 4 P0 object checks, and validate 4 full-weight terminal identities. 15 invalid requests reject.

This closes the actual two-frame leading Rm target/inverse/partial-patch gap. Whole-Z functions, full radial derivative/cell providers, active finite-N patch and all24 integrals, one compatible global N, complete heat/stress/flat layers and actual n-dependent recursion/pulses/corrected NS remain open. No whole-axis, recursive-scale, admissible-stress or full corrected UVW percentage is inferred from this scoped gate.

## Detailed next production tasks

- [x] **ACTUAL Rm FIVE TARGETS.** Translate actual centered Rm histories into original normalized source rows, canonical Am/Rm/Pstar factors, separate P0 and signed tails. Bind live offset/name/basis/source rather than legacy enclosure overlap.
- [x] **NORMALIZED LEADING INVERSE Z0..5.** Perform common-basis anisotropic division before ordinary materialization; solve the original leading map with strict inclusion, contraction and the same Jacobian; map all controls/derivatives back. Keep finite-N targets distinct.
- [x] **PARTIAL LEADING PATCH / PHYSICAL PRIMITIVES.** Original beta f/g, directed partial weights, five changes, own physical moment prefactors, centered energy recovery, P0, Ur Z0..4 and Utheta/Uz first-y rows. Retain tiny correction sectors and conditional full-weight identities.
- [ ] **NEXT FULL RADIAL MIXED4 PATCH PROVIDER.** Recover beta ordinary radial derivatives through4 in directed form, including support edges and flat endpoint behavior. Differentiate partial primitive equations using the same coefficient functions, apply x*d/dx before physical prefactor dressing, recover Q through total mixed order4 and export raw_patch_rows. Do not claim full smooth patch from first-y rows alone.
- [ ] **NEXT PATCH / REFERENCE / RESTORATION RADIAL CELLS.** Extend exact rational point functions to ordered closed source cells. Preserve source logs, denominator positivity, signed remainder and true sigma/beta primitives. Connect original inertial/stress-direction/fixed-phase Z backends using the correct all-N phase and both N-dependent levels. Validate complete source integrals rather than fitting density caps.
- [ ] **NEXT FINITE-N Rsh->Rm FIVE DENSITIES.** Integrate genuine new reference, restoration and postrestore paths with incoming histories and analytic P0. Accepted Rh_reference integrals are a different source window. Add changed-scope directed receipts and source graphs.
- [ ] **WHOLE-Z ACTUAL CORE / MICRO SOURCE ORACLE.** Replace the two cached source frames by source-correlated axial cells from analytic core/fixed-point functions, anchored pole primitive and fourteen pressure atoms. Carry ordinary higher derivatives and errors, midplane/pole splits and uniform bounds. Sample jets cannot prove entire-axis closure.
- [ ] **CONTINUOUS MICRO / MACRO / SWITCH / LONG-RESHAPE PROVIDERS.** Complete both original micro smoothings as functions/cells with dy=hb*ds and their real inlets. Extend long reshape from endpoint kernels to all0<=y<=T using actual anchored B, full finite kernels, signed incoming tails and P0. Preserve exact current-mode powers through both switches.
- [ ] **FINITE-N UNIQUE Rc REPAIR / OUTER FEEDBACK.** Assemble true complete r(N,Z), then solve B*h+N*r+Q/N=0 with its own Jacobian, uniqueness, derivative/error bounds and feedback. Use the leading patch as background source only; never substitute leading h or broad caps for finite-N targets. Bind paired live owner, exact N and actual target integrator.
- [ ] **ALL24 SOURCE WINDOWS / TRUE Rc_E AND Rc_E_Z.** Complete active[1,71/40] and accepted terminal support, reference/O2/axial/buffer/O3 inlets with original raw/inertial terms and both N-dependent levels. Preserve source memory through each path, then compile the actual complete finite-N target.
- [ ] **FIVE TERMINAL FUNCTIONS / ONE COMPATIBLE GLOBAL N.** Close five identities as functions of Z and reconcile repair/cone/interface/centered-target N^-2/Picard/limit/tail inequalities for one finite N. Point leading closure does not satisfy this task.
- [ ] **MATCHED BACKGROUND / HEAT / ENERGY / STRESS / FLAT.** Complete annular, flatten and pressure joins, exact heat exterior, axis regularity, finite-energy tail and regional admissible stress/flat remainder. Background residual magnitude alone is not stress cancellation evidence.
- [ ] **REAL n-DEPENDENT RECURSION / MEANS / TWO PULSE FAMILIES.** Recover actual n=1 and n>=2 equations with independent repairs and divergence-preserving cutoffs/smooth summation; build averaged quadratic stress cancellation.
- [ ] **CORRECTED CARTESIAN NS / PHYSICAL DYNAMICS.** Export corrected u/v/w, independently validate the full forced residual and finite energy, and measure radial contraction, relative axial elongation, true scale recursion and material winding separately.

Mark DONE with defining code, scoped report/receipt and a commit. Continue production after necessary changed-scope checks; the full objective stays active.
