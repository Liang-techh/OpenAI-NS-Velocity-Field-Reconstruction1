# Actual Rm active source atlas and known five-density C0/Z contributions

Checked source [58a422f3](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/58a422f3a86e4d2e7f2df547367696385c0714ab). Full reconstruction **ACTIVE / INCOMPLETE**. Predecessor: [whole terminal local integrals](CURRENT_ORIGINAL_RM_TERMINAL_INTEGRALS_2026_10_09.md).

The actual finite-N density integration now reaches into the active angular-control region. At conditional Z=.5, all six active cells are enclosed; together with the accepted terminal operator this gives a complete **local integral bound on [1,e]**, for C0 and ordinary Z rows of all five densities at candidate N=257. At Z=0, twelve refined active cells are enclosed and four remain explicit unknown source integrals. Neither result supplies the real incoming finite-N correction at Rm, sharp five-terminal closure, a whole-Z field, a cone/global-N theorem or true n-recursion.

## Exact active atlas and same-owner integration

The three original gamma bumps have centers5/4,3/2,7/4 and radius1/40. The source support partition is

    [1,49/40,51/40,59/40,61/40,69/40,71/40].

It contains all six exact flat support edges. The actual mixed4 source preserves the flat value and derivative identities at those edges. Each nonempty cell is evaluated through the accepted original Rm radius phase, generic source quotients, q axial jets, Z-only original phase inverse and five signed densities. There is no old-owner/repair fallback, phase midpoint or sampled field value.

The accepted whole-patch quotient receipt proves Delta<0 and sigma=1 for both actual frames across [1,e]. A new q-cutoff branch adapter is unnecessary for this source. The remaining unknown cells concern the signed-u conditioning in the original phase inverse, not a failed q-cutoff theorem.

The new stateless integrator projects the accepted terminal `coordinate` and `integrate_source_cell` AST. Only the exact chart lower endpoint changes from71/40 to1; the integration body, rates, masses, widths, source-family/closed-cell checks, exact Rm*x/P0 bindings and Z-independent phase guards remain unchanged. The whole atlas must cover [1,71/40] and retain every original support edge. Refinement uses exact rational midpoints, preserves adjacent cells, and is bounded at declared depth3 in the native report.

For an active cell [a,b] and join q=71/40,

    w=log(b/a), suffix=log(q/b), dy=dx/x
    I_j=exp(-lambda_j*suffix)*integral_a^b exp(-lambda_j*log(b/x))*f_j(x,Z)*dx/x
    lambda_m=1, lambda_h=lambda_k=3/2, lambda_e=1, lambda_p=0.

Actual whole-source C0/Z ranges are multiplied by positive own-rate kernel masses. No extra R, Pstar or N normalization is introduced. E,V and every density already use the original common Pstar units. The exact same five formal bases, context, ledger and live P0 object remain attached; P0 is carried separately from the pressure-density integral. Pressure decay is exactly1 through the active join, terminal tail and Rm incoming operator.

Known active contributions are transported through the whole terminal suffix before adding the accepted terminal integral. If U_j denotes the sum of unresolved active integrals ending at71/40, the actual history still requires

    deltaH_j(Rh)=exp(-lambda_j)*deltaH_j(Rm)
                + known_j + exp(-lambda_j*log(e/(71/40)))*U_j.

The implementation records U_j as explicit unknown source terms. It never treats them as zero or includes them in a claimed complete integral. `transport_supplied_incoming` rejects an absent Rm correction and rejects any unresolved active integral, even if an incoming vector is supplied. Actual leading Rm, join and Rh histories are separate comparison/source memory packets; none is substituted for the finite-N correction inlet.

## Concrete unresolved cells at Z=0

Exact numerator/denominator pairs:

- [49, 40] to [197, 160]
- [203, 160] to [51, 40]
- [69, 40] to [277, 160]
- [283, 160] to [71, 40]

These cells touch support edges of the first and third bumps. Their original phase records report `requires_signed_source_refinement` for u=p2*q/dstar. This status means the interval does not certify either |u|<=1/4 or a strict sign with |u|>=1/8. The checker verifies that the u coefficient enclosure spans zero. **It does not establish an exact p2 root or an actual singularity.** Do not infer either from an interval sign.

The source supplies a valid closed-cell input record but no executable primitive/density enclosure on those boxes. Every finite closed partition still retains an endpoint-touching cell, so uniform bisection alone does not prove the required small-u or signed-source theorem. Focus on the original source carrier and flat edge estimates, preserving true nonzero Z rows. All broad/full-prefix/repair/whole-Z/global-N flags remain false.

The defining gamma_rows in current_original_Rm_patch_mixed4_cells.py and beta_jets in actual_patch_mixed_C4.py retain derivatives through order4. The accepted flat_pulse_derivatives.py endpoint bounds use w=1-r^2 and the common exp(-1/w) factor with original normalization and BETA_POLYNOMIALS. Preserve those correlated factors and radius40 derivative factors when bounding the edge source. An independent interval [0,upper] for each beta derivative loses this relationship. The existing partial_weight beta-power integrals are directed integral bounds, not point values or selected controls.

## Evidence and scope

Files: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_whole_density_integrals.py, .json.gz, _check.py, _check.json. Producer 137.047s; scoped checker 137.703s, terminal exit0; **1190 staged exact dependency hashes PASS**. Existing GPT-5.6 Luna/max performed read-only source review; root owns implementation/compute/Git.

The atlas has 18 resolved active source cells across both conditional frames, giving 180 weighted cell C0/Z rows. It reports 20 summed known active rows, 20 composed known whole-patch rows and 6 leading memory packets. All four unresolved source cells and their original status/geometry/P0 are retained. Only Z=.5 has a complete local whole-patch integral bound; Z=0 remains partial.

Independent signed affine-in-logx sources at two Z values and both signs give 240 comparisons against exact high-precision weighted integrals. They cover target-aware active/join weights, active-to-terminal composition and nonzero incoming C0/Z memory for all five rates. Invalid source-chart coordinates and partitions lacking edges/endpoints are rejected. These finite fixtures are diagnostic only. Live native outputs agree exactly with the scoped report and preserve all original source ownership. Bounds may be wide or span zero; completeness of a local integral enclosure does not prove a small moment defect.

## Detailed production tasks

- [x] **RM-A1 EXACT ACTIVE SOURCE ATLAS.** Install the exact support-edge partition, bounded rational refinement, actual phase/source cells and a gap-free ordered union of resolved plus explicit unknown cells. Preserve q, b/b_Z, true p2, P0 and signed densities. Keep current frame scope explicit.
- [x] **RM-A2 KNOWN ACTIVE / TERMINAL COMPOSITION.** Integrate every resolved cell with true dx/x, original own rates and suffix to71/40, then transport through the full terminal suffix. Retain unknown integral terms and Rm inlet memory. Z=.5 complete local [1,e] C0/Z bounds; Z=0 partial.
- [ ] **RM-A3 FIRST PRIORITY: FLAT SUPPORT EDGE SIGNED-u SOURCE CARRIER.** Resolve the four exact Z=0 cells above using the same p2/q/dstar source. Recover correlated C0 cancellation and original flat bump logarithmic bounds rather than an independent interval difference. Prove a true small-u central layer or strict signed branch on each region, with no midpoint/end-point choice as a field definition. Carry the actual p2_Z/q_Z rows across any small-u/zero regime. Use analytic exact thin-layer geometry if the original flat theorem supports it; do not spend repeated broad uniform bisections.
- [ ] **RM-A4 ALL-u ORIGINAL PRIMITIVE / Z COVER.** If a source cell genuinely covers small and both signed-u regimes, build a certified conditional union using the original small-r and Mobius formulas and their overlap. Establish correct positive denominators and directed factor lower bounds in each branch, then hull outputs in the same algebra. Do not convert a source-sign uncertainty into an asserted actual root, replace nonzero Z jets by zero, or transfer an old-owner margin.
- [ ] **RM-A5 COMPLETE Z=0 LOCAL INTEGRAL.** Enclose every remaining source cell and integrate its actual densities with each own-rate weight. Replace the explicit U_j terms only with proven source integrals. Keep complete-local flags false until the exact [1,e] atlas has no missing cells. Reconcile this with the already complete local Z=.5 frame.
- [ ] **RM-A6 SHARP INTEGRAL BOUNDS.** Measure widths in the original formal units. Refine actual source/radius/phase cells and prove any weighted phase-mean cancellation independently, keeping original slow terms and both N levels. No assumed zero oscillatory integral; no small-defect claim from a full-period range cover.
- [ ] **RM-A7 CUSTOM REFINEMENT TRANSPORT API.** Propagate the chosen max_depth through transport_supplied_incoming so a successful custom deeper atlas is reused instead of recomputing canonical depth3. The current scoped receipt covers canonical depth3; custom deeper full transport is not claimed by it.
- [ ] **RM-A8 REAL FINITE-N INCOMING AT Rm.** Build the upstream correction histories through inner/micro/macro/switch/reshape/reference/restoration, retaining every incoming tail and common P0. Supply genuine C0/Z correction bounds at Rm; the current leading memories and an ad hoc zero vector are not that source.
- [ ] **RM-A9 COMPLETE Rc CUMULATIVE FUNCTIONS / ALL24.** Combine complete active and terminal integrals with the real prefix and outer windows. Keep each own-rate suffix and both original levels; produce Rc_E/Rc_E_Z and all24 source functions under one family and explicit N.
- [ ] **RM-A10 GENUINE MIXED y/Z JETS.** Add original fixed-phase y derivatives and true fast N*A_phi/N*B_phi chains. Differentiate cutoffs on source-proven branches and attach each physical R/Pstar factor exactly once. Do not infer missing y rows from the Z-only adapter.
- [ ] **RM-A11 ACTUAL dstar/CONE AND FINITE-N REPAIR.** Prove correlated signed H0,D,J and relaxed-cone inequalities for the actual owner; solve B*h+N*r+Q/N=0 with the actual cumulative sources. Select one global frequency and prove terminal conditions as functions of Z, including poles/midplane. Selected eta/dstar remain parameters until these proofs close.
- [ ] **RM-A12 WHOLE-Z / PRESSURE / HEAT / ENERGY.** Install the functional whole-axis source controls, analytic preheat pressure and exact heat joins; verify finite energy and compatible tails without pressure-tail fitting.
- [ ] **RM-A13 STRESS/FLAT AND REAL n-RECURSION.** Construct admissible stress and flat remainder, then the actual n=1/n>=2 coefficient equations and independent moment repairs; use divergence-free truncation and controlled smooth summation. Coordinate scaling alone is insufficient.
- [ ] **RM-A14 PULSES / CORRECTED CARTESIAN NS / DYNAMICS.** Add mean and oscillatory corrections and averaged stress cancellation. Validate the corrected Cartesian equation, contraction, relative axial elongation, material winding and finite energy across physical times.

Mark DONE with defining source, scoped report/receipt and commit. Keep the full objective active. Proceed from the four actual unknown source cells; do not rerun ancestor producers/checkers or repair/pressure owners.
