# Current O3 own repaired histories, pressure, radial source and functional exit

Implementation and source receipts: commit [90a587ac](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/90a587ac45de1a450af738c00b831543996c6d27).

## Constructed result and scope

The independent quiet-power two-axial/three-swirl repair is now installed as a callable source. This adds actual ordinary logR derivatives through order4, axial derivatives through order5, continuous partial bump primitives, all five repaired cumulative histories, same-axis pressure and radial recovery from the field's own repaired M. **At q=2, all five history increments and the velocity/pressure/radial increments equal zero as functions of Z[-1,1], through the retained derivative orders.**

The controls remain the same unique exact implicit vector defined by the actual modulated source defects. Returned boxes enclose those controls; they do not assign independent values for each query or replace the source by midpoints. Original incoming moments and the original pressure datum remain. The positive modulation kinetic history and nonnegative repair kinetic source remain, even when the signed total energy-moment defect cancels.

This is a scoped local field/moment/pressure/radial exit construction. The quiet-power callable domain ends at q=2; an open q>2 continuation and the full modified physical dispatcher are not yet installed. The full modified completed tensor, energy history, tensor/physical interface joins, downstream analytic heat compatibility, common finite cone N, O2 taper and modified closed O3 cone remain open. Original regional counts stay 15 strict nonzero whole regions plus a separate exact zero exterior, with17 original whole regions open. Inventory33/32/14 and primitive atlas14/8 are unchanged. Those original counts do not admit this modified field.

## Source API and dependencies

`CurrentO3RepairedHistories(repair=checked_independent_repair)` consumes the existing checked warm source graph. Call `.history(region,Z,coordinate)` for:

| Region | Coordinate | Construction |
| --- | --- | --- |
| O2_buffer | buffer offset[0,11] | Same actual modulated O2 field and own histories |
| O3_slope_mu | transition offset[0,1] | Same actual modulated O3 variable-slope field and own histories |
| quiet_O3_power | q[0,2] | Original power plus independent compact repair and own partial histories |

The focused controller stage is `currento3repairedhistories`. Files use `experiments/root_st073/lei_ren_part1_paper_compliant_current_O3_repaired_histories` with `_operator.py`, `.py`, `_check.py`, `.json`, `_check.json`, and `_views.json.gz`. The compressed views contain the whole repair band/full Z, fresh interiors of each support and the exact q=2 boundary.

Keep the parent chain:

`checked original current graph -> actual O3 direction/finite-N modulation -> own modulated histories -> independent source-correlated implicit repair -> own repaired histories`.

Reuse session64009 and its live checked graph on the current host when available. Session IDs are host-local; another host must not assume it can resume that process. Reconstruct only the necessary checked dependency chain if no warm graph exists. Do not run the cold all-stage controller merely to repeat a completed local stage.

## Fixed source and partial primitives

Let y=q-1=logx, x=R/R0, R0=Rw*exp(1), t=log(R/Rd)=2+y. The exact original amplitude is `A(Z)=Pstar*Ua*C(Z)*f2`, with C=1/(1+Z^2), f2=exp(-1-3mu/2). The three fixed logx centers are1/5,1/2,4/5 with radius1/40. The source is `g_i(x)=exp(-y)*beta(y-c_i)` using the same exact normalized raw flat bump. Axial controls use the first and last swirl supports, so both same-center products remain.

Ordinary logR rows use the product derivative

`d_y^n g_i=exp(-y)*sum_{k=0}^n binomial(n,k)*(-1)^(n-k)*beta^(k)(y-c_i)`.

Raw beta Taylor rows are converted to ordinary derivatives using k!, radius^(k+1) and the shared raw normalization. The exact flat-edge helper is reused; no inverse of an interval crossing a support edge is used.

Partial single weights integrate `beta(w)*exp(p*(c_i+w))dw` up to the actual endpoint. Partial product weights integrate `beta(w)^2*exp((p-1)*(c_i+w))dw`. Zero below support and the same checked H/H2 weights above full support are exact branches. The correlated axial weight is the signed partial integral of `(exp(-mu*logx)-1)/mu`, formed as `-logx*integral_0^1 exp(-mu*logx*r)dr` before enclosure. Its sign clamp requires positive mu, wholly positive logx support and a strictly negative full D weight. Partial values use128 directed range subdivisions per raw support; the same256-cell checked repair normalization/full weights are retained.

Write sa=sqrt(mu)/N, se=mu/N and h=(a0,a2,e0,e1,e2). For mass w, divided D, linear weights I/S/Cp, cross c, kinetic e and pressure p, the actual partial primitive increments are:

| Primitive | Partial increment |
| --- | --- |
| M | sa*(a0*w0+a2*w2) |
| J | sa*(a0*(w0+mu*D0)+a2*(w2+mu*D2))+sa*se*(a0*e0*c0+a2*e2*c2) |
| I | se*sum e_i*I_i |
| S | sa^2*(a0^2*e0_weight+a2^2*e2_weight)-se*sum e_i*S_i-se^2/2*sum e_i^2*energy_i |
| Cp | se*sum e_i*Cp_i+se^2/2*sum e_i^2*pressure_i |

The names e_i above are swirl controls; energy_i are distinct positive product weights. All overlap cross and both square terms are kept.

## Actual history units, pressure and radial recovery

If previous B denotes the actual modulated cumulative scalar history at t=2+y, the new scalar increments relative to the original source are:

`m=B_m+f2*M/x`, `h=B_h+f2*I/x^(3/2)`,

`k=B_k+f2^2*J/x^(3/2)`, `e=B_e+f2^2*S/x`, `p=B_p+f2^2*Cp`.

The exact five scalar factors in the normalized physical histories remain `Pstar*Ua*C`, `Ua*C`, `Pstar*Ua^2*C^2`, `Ua^2*C^2`, `Ua^2*C^2`, respectively for m,h,k,e,p. The original normalized m/k sector is Pstar^0 and the new increment is Pstar^1. Pressure is the original absolute axis datum plus original cumulative pressure plus the new own Cp increment.

The ordinary history derivatives come from their actual local density ODEs: m'=Vhat-m; h'=du-3h/2; k'=Unew*Vhat-3k/2; e'=Vhat^2-(2Uold*du+du^2)/2-e; p'=(2Uold*du+du^2)/2. The code retains each derivative row's axial Taylor source.

Raw physical Uz=Pstar*Vhat. The radial increment follows the own-moment divergence-free formula

`deltaQ/Pstar=[2Z*Vhat-(1-delta)*Z*delta_m-(1-Z^2)*d_Z(delta_m)]/(1-delta*Z^2)`.

The original nonzero incoming radial velocity is retained in Pstar^0; the new radial source is Pstar^1, with the original sqrt(R/2) factor and ordinary half-shift derivative rows. Never feed Vhat directly into a raw-Uz stress operator.

## Functional exact exit: source binding, not chosen defects

The actual independent repair has scaled defining equation `B*h+Q_N(h,h)=-d_N`, with d_N computed from the same field's `.scalars(2)`. Its correlated D row is the prior source-bound identity `actual_D_is_same_signed_history_difference_after_transport_and_scaling`. The new constructor explicitly requires that accepted identity and the unique actual-source implicit-root/no-midpoint certificates.

The proof binds the live AST chain `self.repair.histories -> field.scalars(2) -> correlated_scaled_defects -> self.repair.controls`, the parent's fixed support-end integral and actual post-support exp transports, and the exact full-support partial weights of this same map. At inlet:

`Bm2=f2*sa*dM`, `Bh2=f2*se*dI`, `Bk2=f2^2*sa*(dM+mu*dD)`, `Be2=f2^2*se*dS`, `Bp2=f2^2*se*dCp`.

After modulation support ends these transport by exp(-y), exp(-3y/2), exp(-3y/2), exp(-y),1. At full repair support, write r=B*h+Q_N(h,h)+d_N. The actual callable scalar totals factor as:

`m=f2*exp(-y)*sa*rM`, `h=f2*exp(-3y/2)*se*rI`,

`k=f2^2*exp(-3y/2)*sa*(rM+mu*rD)`, `e=f2^2*exp(-y)*se*rS`, `p=f2^2*se*rCp`.

The exact implicit vector has r=0. Defects are not defined as minus the map; they are the existing actual source defects. The q=2 field branch is therefore an algebraic rewrite of the same continuous field, with unreduced interval sums retained for diagnostics. A box crossing q=2 is not assigned zero throughout. Profiles are exactly flat before this exit, so the own ODEs imply vanishing derivative rows as well. Pressure and radial increments then vanish for all Z through the retained orders.

## Evidence and gates

Evidence: 94 source-bound symbolic identities/guards; actual full-Z exit and fresh support views; same implicit controls and full weights; before-support equality; positive modulation kinetic retained; exact-exit-crossing exclusion; producer/checker; checked constructor/API; focused controller and compilation. Working/index source hashes matched 840 dependency files. The same read-only GPT-5.6 Luna / max reviewer confirmed the actual-source connection gap is closed within this q=2 scope.

Four scoped gates are admitted by the new receipt: actual bump logR4/Z5 source, own repaired partial histories/pressure/radial, modified five-moment repair, and quiet-power functional exit to the original source. `current_O3_transition_finite_N_modified_profiles_installed` remains false because it denotes the complete modified physical dispatcher/common-frequency installation, beyond these local callable rows. Completed modified tensor/physical joins, modified analytic preheat compatibility, modified whole cones/global tensor/full NS remain false.

The demonstration N=10^12 still meets only the repair threshold27,303,666. It must not be promoted to the common stress-cone frequency.

## Executable tasks for the next agent

Mark a task complete only after its source operator, producer/receipt and scoped call path are connected. Record a commit and remaining scope. Preserve unrelated dirty files and keep the full long-term goal active.

- [x] **REPAIR5:** actual log-translated bump ordinary logR4 and original common-amplitude axial5 source rows with flat edges.
- [x] **REPAIR6:** continuous partial single/product/divided weights, exact same-map full branches and all overlap products.
- [x] **REPAIR7:** own five repaired histories with actual incoming data, homogeneous transport and correct R/Pstar factors.
- [x] **REPAIR8:** same-axis own pressure and radial recovery from repaired M, callable sectors/derivatives.
- [x] **REPAIR9a:** source-bound q=2 functional five-moment/velocity/pressure/radial exit through retained rows.
- [ ] **TENSOR1 — Raw stress sectors.** Own a new modified-pre-stress operator/class/checker; consume `CurrentO3RepairedHistories`, not original pre data. Assemble raw Uz=Pstar*Vhat, m/k Pstar^0+Pstar^1 and actual absolute pressure. Replay full paper stress (3.16)–(3.18) before collecting monomials. Theta retained k and u*m generate Pstar^1/Pstar^2; axial local/shear are Pstar^1, linear m has Pstar^0/Pstar^1, V*m has Pstar^1/Pstar^2, energy/pressure Pstar^2. Preserve signs and shifted ordinary stress rows through order3. Do not freeze local logarithmic amplitude derivatives or use Vhat as raw Uz.
- [ ] **TENSOR2 — Modified remainder.** Use radial velocity with both Pstar^0 and Pstar^1 sectors, theta and axial Pstar^1. Retain all Pstar^0/Pstar^1/Pstar^2 radial quadratic cross terms and local axial transport, plus each physical beta/radius exponent. Export signed ordinary remainder rows through order2 from the same velocity, including axial viscosities. Compare coefficient programs to the unchanged full remainder AST.
- [ ] **TENSOR3 — Own complete energy.** Rebuild the energy density and continuous partial/future energy integrals for the modified field, with every radial/axial square and cross term. Original interior energy cannot be copied after velocities change. Preserve the accepted original terminal energy constant only when justified by the exact unchanged future source; prove forward/backward integrals describe one modified function. Do not confuse repaired S moment closure with total kinetic-energy equality inside the band.
- [ ] **TENSOR4 — Complete diagonal/radial tensor.** Use own stress, remainder and energy in the original tensor completion/divergence formula. Keep all radial tensor components, scalar offsets and energy transport terms. Bind completed tensor to the actual modified source and full absolute pressure. Recheck only new operator identities and required source receipts.
- [ ] **TENSOR5/REPAIR9b — Open continuation and exit joins.** Install exact original power continuation for q>2 up to actual Tw, using the proven common q=2 trace/flat differences. Complete modified field dispatch on affected O2/O3/quiet-power charts and inherit the original remainder of the chain through source-function equality. Then prove completed tensor, physical spatial4/time1 exit traces and downstream analytic preheat/heat exterior compatibility. The q=2 local equality alone does not certify these.
- [ ] **TENSOR6 — All affected interfaces.** Prove O2 buffer entrance/taper boundaries, O2-to-O3 transition, slope-to-quiet-power and all repair support traces on the same source. Split boxes at chart/flat edges. Keep inherited original interface inventory distinct from newly certified modified interfaces. Use function identities or continuous bounds, not sampled endpoint matches.
- [ ] **TENSOR7 — Phase-aware derivative estimates.** Derive uniform bounds for NlogR modulation jets and repair profiles, retaining genuine N derivative growth. Identify shear terms, pressure/history terms and remainder scales that need correlated cancellation; generic pointwise O(1/N) amplitude alone is insufficient. Keep actual positive mu and original logarithmic slope.
- [ ] **REPAIR10/MOD7 — One finite common N.** Combine repair threshold with all derivative/shear, pressure, completion, interface, O2 taper and two-vector cone inequalities. Fix supports before selecting N and publish a sufficient finite integer. Do not use the repair-only threshold as cone admission.
- [ ] **REPAIR11a — Independent O2 taper.** Prove the affected taper throughout its closed relevant interval and all Z, including flat support edges, using full signed current stress/energy/pressure and correct source modes. No borrowed original-source cone gate.
- [ ] **REPAIR11b — Closed modified O3 cone.** Prove both strict signed two-vector inequalities on the whole modified transition and repair/power source with actual nonzero stress. Treat the original v=0 obstruction explicitly through the new source, retaining full theta and axial components, moment transport and shear. Only promote gates after continuous proof and common N.
- [ ] **GLOBAL1 — Remaining original/modified inner cones.** Complete O2/Rh, patch/restore/reshape/switch/bridge/core and any newly affected charts, preserving exact zero exterior separately. Construct one smooth global admissible completed tensor and its physical lift.
- [ ] **REC1 — True n=1 recovery.** Implement the first higher-order coefficient source/recovery equation with correct paper coefficients, common core interval and independent five-moment repair. Fix the known nonflat leading-origin remainder; coordinate scaling of leading field is insufficient.
- [ ] **REC2 — n>=2 recovery.** Implement n-dependent equations, per-order moment repairs and cutoff at streamfunction/vector-potential level. Verify source-function compatibility and exact divergence-free recovery at each order.
- [ ] **REC3 — Estimates and smooth sum.** Prove finite-order remainder bounds and a valid smooth sum; distinguish temporal flatness from compact spatial flatness. Record convergence domain, coefficient bounds and admissible forcing behavior.
- [ ] **WAVE1 — Two actual oscillatory families.** Construct both homogeneous pulse families with actual positive amplitudes and flat edge weights, matching the completed admissible stress. Bind covariance integrals and finite-frequency errors to their real velocity functions.
- [ ] **WAVE2 — Mean and quadratic cancellation.** Install signed divergence-preserving lifts and mean corrections; prove averaged quadratic momentum-flux cancels the background stress and control nonaveraged errors.
- [ ] **PHYS1 — Resolved u/v/w.** Export a common actual 3D/time velocity evaluator with similarity/Cartesian maps and pressure, retaining anisotropic scaling, axis regularity and physical tail. Remove unresolved parameter/callback substitutes before claiming resolved physical fields.
- [ ] **PHYS2 — Full corrected NS and energy.** Validate the actual corrected velocity, pressure and smooth forcing independently in Cartesian coordinates; report full residual Linfinity/L2 against the original1e-3 targets, finite total energy and tail behavior. Background stress decomposition is a preceding construction, not the final residual test.
- [ ] **DYNAMICS — Measured recursion and winding.** Measure radial contraction, axial/radial aspect ratio, swirl/axial/vorticity scaling and recursive relation across actual coefficient levels. Separate instantaneous streamlines from accumulated material-particle winding. Animation remains secondary.

Immediate owner: TENSOR1 and TENSOR2 on this checked modified source, then TENSOR3/TENSOR4 and REPAIR9b. The full long-term objective is active and incomplete.
