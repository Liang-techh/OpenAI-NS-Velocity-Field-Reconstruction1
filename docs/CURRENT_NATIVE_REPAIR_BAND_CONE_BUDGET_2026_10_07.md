# Original repair-band partial functions and quiet/power cone budgets

Successor: [CURRENT_NATIVE_Q_FLAT_CONE_BUDGET_2026_10_07.md](CURRENT_NATIVE_Q_FLAT_CONE_BUDGET_2026_10_07.md) completes the restricted original q-flat quantitative margins and merges them with the active/quiet/band C0 local budgets across 24 cells/35 branches. Common N, actual controls/tail, terminal closure, global exterior admission and genuine recursion remain open.

Checked source: [04b6bfda](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/04b6bfdafed4f087e6247605051c5b3b0b343a2d). Predecessor: [continuous modulation state errors](CURRENT_NATIVE_CONE_STATE_ERRORS_2026_10_07.md). Producer 3.110s, focused checker 47.531s. Working/index dependency audit: 1,144 matching hashes. Read-only mathematical reviewer: **gpt-5.6-luna / max**. Existing original owners and the accepted all-N route were reused; no source ancestor constructor was called.

The same-source construction now carries five exact finite-Picard partial history functions on the entire repair band Rc<=R<=2Rc and Z[-1,1], including their first ordinary-Z derivatives. It also supplies conditional control-ball bounds for those histories, full inertial/radial/pressure differences and normalized state errors. The unchanged quiet interval r_plus..Rc and the original baseline under the correction band now have quantitative G1..G4 margins from the full signed O3 power theorem. Sufficient local frequency bounds are combined with the previous source/repair/active-state bounds.

This is a completed local **function and budget layer**. It does not install a numerical original-source fixed point or establish actual terminal closure. The exact finite terminal histories equal the corresponding finite control residuals; they have not been set to zero. All-chart q-flat input margins, globally compatible N, actual controls/tail/field, outer joins and coefficient recursion remain open.

## Exact profiles and partial histories

The authoritative geometry is r_plus=Rw*exp(1), Rc=Rw*exp(2), 2Rc=Rw*exp(2+log2). The checked original reservation has Tw>2+log2, so these intervals remain inside the original O3 power chart. Let x=R/Rc, A=Ac_theta(Z)/Pstar, alpha=1/2+mu, and retain the same original strictly positive mu used by source and repair. The original profile is E=A*x^-alpha, V=0. Using the original compact bumps g_i and the existing finite Picard control functions,

```text
F=(e0*g0+e1*g1+e2*g2)/N
G=(a0*g0+a2*g2)/N
E_N=A*(x^-alpha+F), V_N=A*G
g_i(x)=beta((logx-ci)/ell)/(ell*J0*x)
ell=log2/40, ci=(log2/5,log2/2,4log2/5).
```

The signed original density increments are

```text
fm=dV, fh=dE, fk=(E+dE)*dV
fe=dV^2-E*dE-dE^2/2, fp=E*dE+dE^2/2.
D_j(x)=x^-rate_j*(D_j(Rc)+integral_1^x t^(rate_j-1)*f_j(t)dt)
rate=(1,3/2,3/2,1,0) for (m,h,k,e,p).
```

The implementation uses a distinct integration variable and endpoint, preserves incoming Rc functions from all 24 original cells and retains rate-zero pressure memory. It applies full C1 product rules to A and each control; no A_Z term or meridional cross term is discarded. Original P0 remains separate and cancels only in the stress difference. Defining functions use exact raw-bump integrals and original source graphs; ranges never become defining coefficient values.

For the existing finite control residual res=B*h+N*r+Q(h)/N, exact endpoint identities are

```text
Dm(2)=2^-1*A*res_M/N
Dh(2)=2^-3/2*A*res_I/N
De(2)=2^-1*A^2*res_S/N
Dp(2)=A^2*res_Cp/N
Dk(2)=2^-3/2*A^2*(res_M+mu*res_D)/N.
```

The same identities hold after ordinary-Z differentiation with the full A_Z quotient/product rules. In particular, the k row uses M+mu*D, not the divided D row alone. They turn terminal closure into a precise future fixed-point/tail test without claiming that a finite Picard iterate closes it.

## Conditional continuous band bounds

The existing same-source repair theorem supplies a common C1 ball radius rho_control for all five control value/Z rows. Within each axial or swirl family the three bump supports are disjoint. Since |beta|<=exp(-1) and |beta'|<=8/exp(2), the directed bounds are

```text
g_max=exp(-1)/(ell*J0)
g_y_max=8*exp(-2)/(ell^2*J0)+g_max
|F|,|G|,|F_Z|,|G_Z|<=rho_control*g_max/N
|F_y|,|G_y|<=rho_control*g_y_max/N.
```

These give profile and density C0/Z bounds with both actual A and A_Z. Each partial moment is bounded by its incoming all-N coefficient history plus its own full positive Duhamel mass on log2. The masses are log2 for pressure, (1-2^-1) for rate1 and (1-2^-3/2)/(3/2) for rate3/2. No history is reset before repair.

The previous full original recovery difference operator produces inertial, radial and pressure error polynomials. The denominator bound uses the actual positive A and power floor 2^-2/3, because 0<mu<1/6. The recorded positivity requirement ensures |F/x^-alpha|<=1/2. Consequently

```text
delta_p=(R/E*(delta_Ilinear+Pstar*delta_Iquadratic)-p*F/x^-alpha)/(1+F/x^-alpha)
delta_a=-2*(F_y+alpha*F)/(x^-alpha+F)
delta_b=2*G_y/(x^-alpha+F).
```

All resulting C0 state error majorants have only negative N powers. The budget applies to any same-source C1 control in the accepted ball, including the genuine fixed point once independently instantiated. It does not assign an enclosure endpoint to a control function and does not identify finite Picard functions with that fixed point.

## Original signed power margin and frequency combination

The accepted original whole O3 power theorem retains full energy, absolute pressure and axial direction correlations. It supplies D>=R*theta_floor and Q/D^2>=2-direction>0. With a=2+2mu, b=0 and kappa-2=2mu,

```text
G1=a>=2
G2=a*2mu>=4mu
G3=a*D>=2*D_floor
G4=a^3*Q>=8*D_floor^2*(2-direction).
```

For the quiet interval use R>=Rw*exp(1); for the band use R>=Rw*exp(2). The actual p1 and p2 covers bound H=2+max(3,|p1|,|p2|). The previously independently checked polynomial gradient estimate preserves every margin at least by half whenever the state error is below rho_cone=min(1,g/(512H^5)). The code consumes the full signed Q theorem; it never asserts p2=0 merely because original V=0.

Quiet-state errors reuse the accepted continuous power_to_Rc history bound: the source profile is unchanged, while incoming modulation histories and pressure remain. Band-state errors use the new conditional control-ball partial histories. Each negative-power polynomial gives a directed sufficient logN bound; the final recorded maximum includes source, repair, active-state, quiet-state, band-state and band positivity. It is a local collection of sufficient lower bounds, not a selected global frequency.

## Evidence and implementation

Focused checks passed 120 independent partial physical C0/Z quadrature comparisons and 40 finite terminal-residual comparisons, at both band endpoints and inside a bump, with nonconstant amplitude, nonzero incoming pressure and full Z product rules. Additional checks cover ten original density/C1 identities, two exact quotient/shear identities, 42 independent bump value/radial derivative comparisons and 16 signed normalized state comparisons. These scalar tests are explicitly manufactured operator references, not evaluated original point fields or fixed points. The live output is bound to the actual original source family, mu, amplitude, continuous route, own-rate masses and all five control-ball rows.

Implementation: `experiments/root_st073/lei_ren_part1_paper_compliant_current_native_repair_band_cone_budget.py`; result and independent checker use the same stem. The optional `current_native_band_function_evaluator.py` caches each expression by its free variables, preserving integration-variable binding and endpoint dependence; a fixed bump integral or control coefficient is reused across outer quadrature points. This improves nested exact-function evaluation without changing the explicit original source/integral oracle contract or installing a fixed-point certificate. Gate: `current_original_repair_band_partial_C1_functions_and_quiet_power_cone_budget_bound`.

The native producer completed with its existing live original owners. The focused checker subsequently uses the saved hash-bound native output and its exact 4,858-node original graph prefix; it does not rebuild live original owners. It independently recomputes directed local frequency/mass conditions at higher precision. The previous Python session 77153 has ended, so do not target it on continuation. Preserve and reuse current saved graph/range receipts; reconstruct actual live source services only when a production task needs fresh native queries, not to rerun accepted checks.

## Executable continuation tasks

Take one bounded production task, preserve the same original live/source definitions, publish its result and independent evidence, then check the box and name the commit. Do not rebuild expensive source ancestors or rerun accepted ancestor checks unless their source or dependencies change. Root owns construction/compute; the existing gpt-5.6-luna/max worker is read-only. Animation remains secondary to the velocity functions and recurrence.

- [x] **CONE3a-partial-functions:** exact repair-band density/C1 and five partial histories, with distinct integration variable and incoming memories.
- [x] **CONE3b-terminal-identities:** exact finite-Picard terminal history/residual relations, including joint k=M+muD and A_Z factors.
- [x] **CONE3c-control-ball-budget:** original bump/radial derivative caps; continuous partial C0/Z moment, full stress and state error polynomials.
- [x] **CONE2b-quiet / CONE3d-power-margin:** full signed original O3 power margins on quiet and repair intervals; combine local frequency conditions with source/repair/active bounds.
- [ ] **NEXT CONE2a-q-flat-domain:** bind original q=0 to kappa-2>=eta on each actual source chart/subset, including boundary values. Record empty subsets where kappa<=2; do not use an inner-exit collar as an entire-chart theorem.
- [ ] **CONE2a-O3-power:** apply the present complete signed power margin to the q-flat portion of power_to_r_plus, retaining its incoming history state errors. Combine with active coverage at the common boundary.
- [ ] **CONE2a-O3-transition:** consume original `O3_transition_full_Ttheta_over_F_log_lower` and `2-full_directional_expression_upper`; on 2mu*sigma>=eta derive G2>=2eta, G3>=2D_floor and G4>=8D_floor^2*(2-direction). Preserve the original endpoint/transition source, full Z and incoming pressure.
- [ ] **CONE2a-O2-axial:** on b^2/2>=eta consume whole-source D/theta and Q/theta^2 positive bounds and the actual canonical theta reserve. The b=0 edges/midplane are outside this subset; never incorrectly assign them strict G2.
- [ ] **CONE2a-inner-six:** derive quantitative D,Q and weighted G1..G4 lower bounds for the entire q-flat subset of bridge_first/second/macro and switch_first/second/power. Existing whole-chart signs and a selected exit collar are insufficient. Retain actual omega,Hbar,kappa and signed phase/Z correlations. If interval dependency prevents a bound, subdivide actual coordinates or prove source inequalities rather than substitute sign assertions for magnitudes.
- [ ] **CONE2c-complete-modulation-cover:** combine active and q-flat certificates across all 24 cells/35 original cutoff branches, preserving source-family identities at joins and a strictly positive quantitative tolerance on each applicable subset.
- [ ] **CONTROL1b-common-N:** enumerate every still-needed source, positivity, cone-region and join condition, keep higher-y terms in their correct N powers and prove compatibility. Produce a whole-frequency receipt only after coverage is complete; do not transfer old modified-source frequencies.
- [ ] **CONTROL1c-original-point-oracle:** implement the accepted original phase inverse and source point/integral oracle for one compatible N. It must evaluate the original graphs and actual coefficient functions rather than phase samples or upper-cover endpoints.
- [ ] **CONTROL1d-fixed-point/tail:** instantiate all five C1 controls at that N with the same target/mu/B/Q, preserve implicit Z derivatives and certify the convergent Picard tail. Report finite residual and tail separately. Replace conditional control-ball qualification only with genuine same-source solved-control evidence.
- [ ] **CONTROL2a-terminal-closure:** apply the exact endpoint identities and independently integrate all corrected five histories/Z rows with a certified tail, proving terminal functional closure throughout Z[-1,1]. A finite iteration or point residual does not close this task.
- [ ] **CONTROL2b-field-pressure-heat:** insert the solved band profile, recover radial velocity and absolute pressure with the unchanged analytic datum, establish compact support and function-level joins at Rc/2Rc through the original exterior.
- [ ] **HIGH/OUTER-global:** finite required higher mixed derivatives, changed source seams, outer annuli, heat collars/exterior, complete stress cone and physical finite-energy tail. Higher radial derivatives may grow with N; do not require their artificial simultaneous decay for C0 cone selection.
- [ ] **REC-coefficients:** implement the genuine n-dependent recovery equations and independent repair at each order on the common inner interval, then controlled truncation and flat summation. Coordinate rescaling of one profile is insufficient.
- [ ] **WAVE-cancellation:** construct the two oscillatory pulse families and mean corrections, measure averaged quadratic stress cancellation, and retain the flat forcing/remainder estimates.
- [ ] **PHYS-uvw/diagnostics:** provide corrected physical [u(x,y,z,t),v(x,y,z,t),w(x,y,z,t)], independent Cartesian NS residual, radial contraction/axial elongation/swirl-vorticity exponents and material winding across time. The full long-term goal remains active and incomplete.
