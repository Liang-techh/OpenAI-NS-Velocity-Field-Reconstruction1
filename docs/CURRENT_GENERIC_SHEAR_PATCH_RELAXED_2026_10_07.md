# Current original implicit patch relaxed input through Rh

Checked implementation: [52d90755](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/52d90755923633eb19692eabf76a8925f0565b98). **LEFT4c1-gates-b-patch-source, patch-supports and patch-gaps are implemented for the current original analytic background.** The new whole Rm..Rh patch condition extends the already attached open inner/reshape/restoration conditions from Ra<R<=Rm to Ra<R<=Rh. This closes the original middle relaxed-input source gap; the complete generic loop input bundle still requires source derivative norms, conservative scales and support reservation. The full long-term goal remains active.

## Current source and domain

`CurrentPatchRelaxedInputs().query(x, Z=(-1,1))` returns an analytic source-domain certificate for x=R/Rm in[1,e], Rm=exp(-6)Rref and Rh=exp(-5)Rref. It covers all three active support groups [49/40,51/40], [59/40,61/40], [69/40,71/40], and the four quiet intervals between the inlet, those groups and Rh. It is not a point velocity evaluator.

The original functions are retained:

    H=x^.1+sum_k xi_k(Z)*beta_k(x), Utheta=Am*H,
    V=4Z+c1(Z)*beta_first(x)+c2(Z)*beta_third(x),
    D_y=x*D_x, a=1-2*x*H_x/H, b=2*x*V_x/Utheta.

All five partial primitive changes, inherited defects, full energy and analytic P0+Mp remain present. A quiet local bump does not zero previously accumulated partial histories or pressure. Only x>=71/40 uses the exact terminal identities of the same unique full-weight implicit map. The existing current mixed4 source theorem already attaches all support traces and the true Rh function join.

The module consumes current source/pressure/five-ODE and normalization receipts, unchanged original callables and ASTs, current exact transported defects/tails and the same implicit coefficient family through axial5. Fixed CA/CQ/CS are universal bump estimates; their historical source admission is not copied. Their delta, epsilon0, Pstar, entry mean/pressure, norm, radius and coefficient-smallness hypotheses are checked against the current source. The analytic summed C1 norm of the correlated true coefficient function is used, rather than the norm of arbitrary independent Taylor boxes. The original stored t is required to bound the exact current 2CAe; the original rounded CS*t cover must enclose recomputation.

## Complete signed condition

The same full primitive-to-stress equations give

    p1=R*Qangular/L, p2=R*N/(L*Utheta), L=1-delta*Z^2,
    J=N/Utheta^2, bw=b*p2/p1=2*J*V_y/Qangular,
    H0=p1*(a-bw)/a.

Qangular is the paper angular inertial quantity, not radial velocity. The fixed source estimate yields SQ>=1.8-delta/2-10epsilon0-CS*t>7/5. The same exact Q equation and current Rm inlet preserve Qangular>=1/2 through active supports and quiet gaps. The current radius proof Rm>=16 and 0<L<=1 imply p1>=8. The normalized swirl remains above1/8 and below1, with a in[.7,.9].

For the axial term, all signed pressure, energy, mean and meridional contributions are retained in the exact N equation. The same full unpatched reference-join pressure C1 bound plus the actual partial pressure change <=3t*Pstar^2 bounds the pressure over the whole patch. The complete N source budget is below2KN*Pstar^2 and the inherited Rm |N| is at most KN*Pstar^2. Positive scalar transport therefore gives |N|<=2KN*Pstar^2. Hence |J|<=128KN and |bw|<=512KN*CS*t. This deliberately derives a conservative full bound instead of transferring the earlier 256KN estimate.

Current directed bounds are CS*t<=1.947543898e-12, |bw|<=1.794856456e-05, kappa<1 and

    p1*(a-bw)-2a >=8*(.7-|bw|)-1.8>0,
    H0-2 >= 4.222062679.

These are sufficient for the original kappa<=2 relaxed branch. They do not construct a strict completed stress tensor, install a new generic loop or certify true n-dependent coefficient recursion. Whole p1/p2 derivative norms and conservative generic scales remain open; the strict tensor registry and historical scoped N>=68,533,403 are unchanged.

## Focused evidence

Passed: 33 current source conditions, 13 original AST bindings, 8 exact signed/partial-history identities, 12 positive current norm/smallness gates, 11 positive margins, 27 independent exact-rational signed H0 envelopes, all seven support/gap boxes and six invalid-domain guards. Working/index audit matched 976 bound hashes. No raw source ancestor graph was rebuilt, and no tiny width or enormous source radius/amplitude was materialized. Read-only review used GPT-5.6 Luna / max.

## Next production tasks

- [x] **LEFT4c1-gates-a:** typed original full I/F, S/F and signed division-free invariants on16 source charts.
- [x] **LEFT4c1-gates-b-inner/middle:** exact current source admission for open inner exit, long reshape/reference, restoration/buffer.
- [x] **LEFT4c1-gates-b-patch-source:** bind actual five-bump H/V, exact same current h(Z), transported defects, five ODEs and P0.
- [x] **LEFT4c1-gates-b-patch-supports:** original whole relaxed condition on all three active support groups, using full signed p2 and current C1 correction bounds.
- [x] **LEFT4c1-gates-b-patch-gaps:** retain actual partial histories/pressure in all four quiet intervals; terminal closure only after the final full support.
- [ ] **LEFT4c1-gates-c-norms-inventory:** list the exact slow and axial derivative orders consumed by the generic loop, inverse and zero-mean primitive. Map each to actual current packet rows and function-level analytic bounds, distinguishing row availability from whole quotient norm admission.
- [ ] **LEFT4c1-gates-c-norms-denominators:** supply positive source-function lower bounds for Utheta, a and any inverse denominator on the actual open loop domain. Keep loghb/logPstar/logR factored. A cache denominator box containing zero must not be inverted or replaced by a scalar cap.
- [ ] **LEFT4c1-gates-c-norms-p:** derive complete p1,p2 bounds and all required derivatives, including axial pressure/energy and meridional sectors, on every original source chart. Preserve signed cancellations before bounding.
- [ ] **LEFT4c1-gates-c-norms-t0:** derive b/a and its required slow/axial derivatives with actual width correlation through bridge, switches, reshape and patch. Do not materialize microscopic a_min or discard finite-width switches.
- [ ] **LEFT4c1-gates-c-scales:** combine actual a, H0 interior margin, strict right-edge excess and admitted p/t0 derivative norms into the conservative generic scales and an executable whole source input bundle. Keep whole_upstream_current_source_packet_assembled false until all hypotheses are supplied.
- [ ] **LEFT4c1-gates-d-support:** locate strict untouched inner collar and strict right edge, reserve an actual later power interval for the new five-moment repair, and choose common loop support/quiet strips. A tapered loop cannot silently lose its required v>2.
- [ ] **LEFT4c1-live-use:** execute selected source point queries with an already existing checked owner graph. Do not rebuild all ancestors merely to test a cache adapter. Core saved rho[.001,4] does not cover the axis.
- [ ] **LEFT4c2-phase:** use one N*log(R/r_minus) through actual chart selectors, with exact radius/width identities and phase held fixed during slow derivatives.
- [ ] **LEFT4c2-jets:** build factored q/r, inverses, zero-mean primitives and whole-cell increments through mixed4, retaining every pressure and meridional term.
- [ ] **LEFT4c3-transport:** propagate all five own changed-history defects over actual formal widths, source seams and quiet gaps. Add explicit rebase identities between factor bases and retain P0.
- [ ] **LEFT4c3-tensor:** export changed Cartesian/time velocity, pressure, divergence-form stress and remainder from those same histories, including full meridional decomposition and all interfaces.
- [ ] **LEFT4d-repair-map:** construct new terminal five-bump defect map, fixed/full functional Jacobian, directed inverse, nonlinear estimates and a uniqueness ball for the modified family.
- [ ] **LEFT4d-terminal-identities:** establish each new terminal moment identity as a function of Z and preserve pressure/source joins. Sample residual containment is insufficient.
- [ ] **LEFT4d-new-N / e:** derive whole modified loop plus repair perturbation budgets; select one new whole-family finite N and close its strict completed cone. Historical scoped N cannot select this N.
- [ ] **CONT / ENERGY:** finish changed function interfaces and finite energy over the actual physical space/time domain and exact heat exterior.
- [ ] **REC:** implement correct n=1 and n>=2 equations, common inner domain, independent order-specific moment repairs and a controlled smooth sum. Coordinate scaling of a single field is insufficient.
- [ ] **WAVE:** implement mean corrections/two-family pulses and prove the averaged quadratic stress cancels the constructed admissible stress.
- [ ] **PHYS:** deliver corrected uvw/p/f, independent Cartesian residual, and measured radial contraction, axial aspect ratio, swirl/vorticity growth and material winding.

New true gate: `current_original_Rm_Rh_relaxed_generic_input_certified`. Whole-upstream/current-loop/changed-moments/new-repair/common-N/global-cone/actual-recursion/full-corrected-NS gates remain false.
