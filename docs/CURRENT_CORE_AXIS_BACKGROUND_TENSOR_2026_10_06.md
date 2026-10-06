# Nonsingular Cartesian core axis T/E (2026-10-06)

The same admitted nonlinear core now supplies the full background stress and all six original NS remainder sectors on **rho[0,4], Z[-1,1]**, including the physical axis. The new Cartesian adapter cancels apparent inverse-radius terms before taking source bounds. It keeps nonzero transverse derivatives and the nonzero axial remainder on the axis.

The completed tensor inventory remains **33 regions, 32 adjacent traces, 10 internal traces**; the primitive atlas remains14/8. The axis is a coordinate extension of the existing core region. This is a spatial source-bound milestone at finite requested log tau and positive viscosity. Global cone/lift, independent terminal-time flatness, required-domain energy, resolved u/v/w/p and actual n-dependent recursion remain open.

## Implementation and use

- Operator: `experiments/root_st073/lei_ren_part1_paper_compliant_current_core_axis_operator.py`.
- Producer/data: `experiments/root_st073/lei_ren_part1_paper_compliant_current_core_axis_background_tensor.py` and deterministic complete `.json.gz`.
- Independent checker/receipt: matching `_check.py` and `_check.json`.
- Focused controller: `python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentcoreaxis`.
- Warm reuse: `CurrentCoreAxisBackgroundTensor(core_tensor=checked_current_core_background_tensor)`.
- `field.axis(Z, log_tau, viscosity)` computes the exact axis.
- `field.cartesian(Z, X, Y, log_tau, viscosity)` accepts the whole normalized Cartesian disk `X^2+Y^2<=8` and Z[-1,1]. Here X,Y are normalized source coordinates, not physical Cartesian coordinates.
- The original cylindrical `CurrentCoreBackgroundTensor.chart()` still requires rho>0. Route axis requests through the new adapter.

Reuse the checked core/first/bridge graph. Do not rerun historical producers for this focused stage. The output contains directed source enclosures and exact factor ledgers, not numerical point values obtained by choosing interval representatives.

## Coordinates and nonsingular remainder

R=epsilon*rho, rho=(X^2+Y^2)/2, X=x_phys/(sqrt(nu)*lambda*sqrt(epsilon)), Y likewise, Z=z_phys/(sqrt(nu)*lambda^(1-delta)). The implicit time map is lambda^2-lambda^(2delta)*z_phys^2/nu=tau.

Let L=1-delta*Z^2, d=1-Z^2 and A_gamma(g)=(gamma*Z*g+d*g_Z-2*Z*rho*g_rho)/L. Q and V=Uz come from the same core; f is Phi dressed once with the original relative F0 derivatives. The exact six source shapes are:

- Radial time: `(Q+rho*Q_rho+(1-delta)*Z*Q_Z/2)/L`, multiplied by `(X/2,Y/2)*sqrt(epsilon)`.
- Radial viscosity: `-(rho*Q_rhorho+2*Q_rho)`, multiplied by `(X,Y)/sqrt(epsilon)`.
- Radial nonlinear transport: `Q*(rho*Q_rho+Q/2)+V*(d*Q_Z-2*Z*Q-2*Z*rho*Q_rho)/L`, multiplied by `(X/2,Y/2)*sqrt(epsilon)`.
- Radial axial viscosity: `-A_(-3+delta)(A_(-2)(Q))`, multiplied by `(X/2,Y/2)*sqrt(epsilon)`.
- Theta axial viscosity: `-A_(-3)(A_(-2-delta)(f))`, multiplied by `(-Y,X)*sqrt(epsilon)*F0base`.
- Axial axial viscosity: `-A_(-2)(A_(-1-delta)(V))`, in E_z.

These source units still carry their original lambda powers. All remainder contributions and their mixed derivatives carry viscosity `nu^(1/2-Ntotal/2)`. Each transverse derivative contributes epsilon^(-1/2) and lambda^(-1); each axial derivative changes the lambda exponent by delta-1. The theta F0 base factor is retained and its relative axial dressing is differentiated exactly once.

The radial-viscosity cancellation is an exact function identity: the original inverse-R expression reduces to the smooth ordinary radial derivatives above. The original common pressure FTC, `C+rho*C_rho=Phi^2`, cancels pressure gradient against centrifugal force before any axis bound. Original P0, moments and all product/nonlinear tails are retained.

## Equality to the positive-radius original lift

The admission replays actual source programs, rather than comparing final interval boxes:

1. The original six raw sectors give36 true Euler/Z source-row identities through total order2. Raw modes, normalization and all110 Cartesian contribution units are source-bound.
2. Independent Cartesian chain rules give30 mixed2 identities against unchanged `core_step`.
3. Factored polar differentiation of X*g, Y*g and scalar g gives30 original cylindrical-to-Cartesian coefficient identities. Angular-unit derivatives are included.
4. Six normalized radial/Z rules are checked against AST-replayed original `physical_operators`.
5. The actual old and new `physical_source_row` call arguments, original `factor_logs`, and actual core packet are replayed for110 exact factor identities. `R=epsilon*rho`, Pstar, zero D/H, and the shared original epsilon/F0 source expressions are bound to production assignments.

The coefficient identity is `Cartesian coefficient=(2/rho)^(Nxy/2)*normalized polar coefficient`. Combining it with the old lift's `(2/R)^(Nxy/2)` factor leaves epsilon^(-Nxy/2). Cancellation occurs symbolically on rho>0; the resulting analytic Cartesian functions extend to rho0 without inverse coordinates in the production operator. Small-radius samples are not the axis limit proof.

Only the **total** original stress is zero, by the separately admitted integrated core equations and analytic source extension. The original eleven signed stress sectors remain in the checked positive-radius parent's ledger. No individual sector is set to zero. Cartesian stress mixed3 and divergence mixed2 are the exact zero total extension; all nonzero E contributions remain.

## Targeted evidence and scope

- Seven full views: whole axis, center, fresh axis, whole near-axis box, fresh near-axis/time/nu, fresh positive radius, and core exit.
- Each view includes120 full stress3,30 divergence2 and110 remainder2 contribution rows: **1,820 physical rows** total, plus45 original mixed4 Q/V/f source jets and their moment/nonlinear/product tails per view.
- Axis transverse vector values vanish by parity while nonzero first transverse derivatives remain; the axial remainder at the center is strictly nonzero.
- Independent checker clears symbolic proof and local adapter caches, rebuilds proofs/complete source outputs, checks original epsilon/F0/lambda/nu factors, and rejects invalid domains and foreign/unchecked owners.
- Focused producer/checker/controller and four changed Python compilation checks pass. Working/index source hashes and staged whitespace are checked before publication.
- Read-only source/math review: GPT-5.6 Luna / max. The final raw-row, angular pullback and physical-factor provenance gaps were resolved.
- Gates now true: `current_core_axis_Cartesian_full_tensor_remainder_mixed2_available`, `core_axis_tensor_remainder_limits_certified`, `current_core_full_background_tensor_available`.
- Global admissibility, stress lift, full NS, terminal-time remainder, required-domain energy, resolved points and temporal recursion stay false.

## Completed axis tasks

- [x] F57C4e-core-axis10 / axis1: same checked core/common/first source graph and family/datum/callable bindings.
- [x] axis2: all six original nonsingular remainder sectors; inverse-R cancellation before bounds.
- [x] axis3: unchanged infinite fixed-point Q/Phi/Uz mixed4 jets, F0 dressing and every source tail on the axis.
- [x] axis4: smooth Cartesian remainder with correct parity and nonzero transverse derivatives.
- [x] axis5: full total stress3, divergence2 and remainder2 extension, retaining the signed parent ledger.
- [x] axis6: true original source rows, polar/Cartesian derivatives and full original physical-factor equality on rho>0.
- [x] axis7 / core11/12-axis: independent complete data/receipt, focused controller and scoped full-core admission.

## Next construction tasks

- [ ] F57C-angular-full1a: identify the four unfinished angular internal seams and their exact native support coordinates. Reuse the checked CurrentAngularBackgroundStress and complete source graph. Map each seam to an existing source-function/flat-support proof; do not promote local difference bounds to a full tensor identity.
- [ ] F57C-angular-full1b: establish whole-source Phi/velocity, all five cumulative histories, full energy and absolute pressure equality across each seam with the actual original signed memory/amplitude. Prove the derivatives required for stress3 and remainder2 from the true programs.
- [ ] F57C-angular-full1c: run the unchanged full tensor, divergence, completion and remainder operators on those same source jets, then create canonical common completed-trace bounds for all components. Use exact function identity before forming triangle bounds.
- [ ] F57C-angular-full1d: independent four-seam producer/checker, fresh axial/time/viscosity requests, complete source hashes and focused controller. Update only the proved full tensor support count; keep the primitive atlas14/8 separate.
- [ ] F57C-global-cover1: one-family/source/datum dispatch of all regional tensor APIs and the core axis, with a complete chart/adjacent/internal seam table and no coordinate gap. Retain compact-time scope explicitly.
- [ ] F57C-cone1: original admissible stress-cone margins and lift on every required region, using full signed stress and exact amplitude factors. Zero core or exterior stress does not certify other regions.
- [ ] F57C-flat-time1: independent time-dependent remainder/derivative estimates and the claimed terminal flat decay, separate from these spatial source bounds.
- [ ] F57C-points1: convergent point representations of the shared nonlinear core and cumulative histories with truncation/error control; return physical u(x,y,z,t),v,w,p without interval midpoint substitution.
- [ ] F57C-energy1: state and integrate the required physical energy domain with original core and full heat tails/localization; do not relabel an unlocalized whole-space field finite-energy.
- [ ] F57D-recursion1: original n=1 coefficient recovery, independent per-order moment repair and pressure/history transport on the common core interval.
- [ ] F57D-recursion2: distinct n>=2 equations with their actual n-dependent operators and at least two nontrivial recovered orders. Leading-field coordinate scaling is not coefficient recursion.
- [ ] F57D-recursion3: finite-order remainder and the required smooth sum; cut off streamfunction/vector potential before curl to preserve divergence.
- [ ] F57E-oscillation1: both oscillatory pulse families and mean corrections; independently compute the averaged quadratic flux and cancel the admitted stress.
- [ ] F57E-corrected1: independent corrected Cartesian NS residual and full corrected L-infinity/L2 targets.
- [ ] F57F-dynamics1: quantitative radial contraction, relative axial slenderness, swirl/vorticity amplification and real material-line winding across times, with recursive coefficients measured separately.

Each task should produce implementation, full source evidence, an independent scoped receipt, completed task status and a commit. The long-term objective remains active.
