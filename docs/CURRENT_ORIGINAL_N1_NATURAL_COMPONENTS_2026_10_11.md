# Original n1 natural component enclosures and own-mean recovery

The axial, pressure and radial component conditioning blocker has advanced. At Z=371/1000, the core exit, actual first inlet and corrected selected inner exit each pass all **14 normalized and 14 physical value/Z1 width targets** for Uz1, K1, P1, d_x_Uz1, own mean1, Q1 and V1. Every enclosure includes the nonseed source errors and the full six-row feedback tail. The previously accepted four F1/d_x_F1 value/Z1 swirl rows remain available. These are selected-point component results; the full reconstruction goal remains **ACTIVE / INCOMPLETE**.

Read this checkpoint before [the previous full-feedback checkpoint](CURRENT_ORIGINAL_N1_FULL_WEIGHTED_FEEDBACK_2026_10_11.md). The historical 0/16 coupled-prefix result is superseded for selected component accuracy by the new conditioned enclosures. This does not certify a complete n1 coefficient function, its own five-moment repair, genuine n2 recursion, the smooth time sum or the full corrected NS residual.

## Genuine leading own mean and radial recovery

The new leading-mean adapter exposes the shared factors before interval subtraction:

```text
Uz = uaxis + epsilon*Psi
mean(Uz) = uaxis + epsilon*meanPsi
Q = Qaxis + epsilon*Qpsi
Qaxis = [(1+delta)*Z*uaxis - (1-Z^2)*uaxis_Z]/L
Qpsi = [2Z*Psi - (1-delta)Z*meanPsi - (1-Z^2)*meanPsi_Z]/L.
```

Here L=1-delta*Z^2 and epsilon=Lambda^-1. Pressure Pstar^2 remains a separate formal sector. The actual mean is its own integral history, using the original core inlet:

```text
rho*meanPsi_actual = 4*meanPsi_core(4) + integral_4^rho Psi_actual drho
drho = hb*rho*ds
meanPsi_r = (Psi_actual-meanPsi_actual)/rho
meanPsi_rr = (Psi_actual_r-2*meanPsi_r)/rho.
```

This does not replace the actual mean by a comparison field. The actual Psi and normalized Phi radial derivatives are reconstructed from the original prescribed-shear equations. Mixed Z grades 5/4/3 are available for Psi, meanPsi and Phi; Q retains grades 4/3/2. Core/actual inlet equality and the exact zero axis factor are checked. The original nonlinear correction, pressure sectors, positive hb and inverse hb in Psi_rr remain present. meanPsi_rr requires no inverse hb.

## Natural six-component normalization

The new normalization is

```text
W = Snew*What
Snew = diag(F0*Lambda, epsilon, epsilon, epsilon*Pstar^2,
            F0*Lambda, epsilon).
T = diag(1, epsilon, epsilon, epsilon*Pstar^2, 1, epsilon)
Snew = Sold*T.
```

T has no x or Z dependence, so each residual/B1 entry is multiplied by T_j/T_i, with no additional gauge. Akin and Achi are unchanged. The genuine normalized pressure and axial source rows become

```text
g4 = 2*x*N1p/Pstar^2
g6 = 2*N1z - epsilon*4Z*rho*N1p/L,
x=sqrt(rho).
```

Source bounds retain the genuine own-mean/Q epsilon factors. Core and actual collar bounds remain separate until integration: the core dx length is 2, while the collar has dx=hb*x*ds/2 and phase length at most 1/2. All inverse hb powers cancel before taking the source norm. The unchanged normalized g5 source is inherited from the accepted full-feedback adapter.

The original outer complex eta is preserved. Matrix bounds use eta/4; the initial M0 source bound uses eta/8; all-word iteration loses eta/16 and returns a value tail on eta/16. A separate eta/32 gap supplies ordinary Z1/Z2 bounds. The admitted all-word Volterra/Cauchy argument is reused with these newly derived norms, rather than copying the old shared absolute tail into smaller components.

| Uniform natural quantity | Conservative upper bound |
| --- | ---: |
| Initial M0 on eta/8 | 1.075649e73 |
| C0 on eta/4 | 1.300001e-319 |
| C1 on eta/4 | 7.000001e-320 |
| All-word ratio Q | 1.287827e-283 |
| Feedback value remainder on eta/16 | 7.948523e-209 |
| Feedback ordinary Z2 remainder on eta/32 | 7.104315e-137 |

These are upper bounds in the stated natural normalization, not physical NS residuals. The 1e-320 norm floor is used only after a directed proof from the original source logs. It changes an upper-bound rounding cover; it changes no source width, F0, pressure or Lambda parameter. Fmax appears only as the uniform complex amplitude norm proxy.

## Correlated signed source seeds and full error enclosures

The viscous pressure source is integrated exactly by FTC across the core and actual collar:

```text
N1p_viscous = Lambda*(2Q_r + rho*Q_rr)
           = 2Qpsi_r + rho*Qpsi_rr
P1_viscous/(epsilon*Pstar^2)
  = [Qpsi(rho)-Qpsi(0)+rho*Qpsi_r(rho)]/Pstar^2
Qpsi(0)=0.
```

The complete pressure also has an inertial source. Its integrated complex bound is included in the error; it is not set to zero. At the corrected exit the pressure value relative width is at most 1.557431e-28, and Z1 relative width at most 4.517835e-10. Pressure Z2 is supplied as an enclosure, but its seed interval can cross zero, so no relative Z2 width admission is claimed.

The dominant axial source is C(Z)=-2*T_-2*T_(-1-delta)uaxis. Its regular principal solution is

```text
Uz1/epsilon = C*rho/4
K1/epsilon = -C*rho/8
d_x_Uz1/epsilon = C*sqrt(rho)/2
mean1/epsilon = C*rho/8.
```

The nonaxis epsilon*Psi source and pressure-source contribution are integrated with complex bounds. Those errors and the full six-row feedback remainder are added to the signed seeds. No physical correction is fitted to improve the width result.

The exact temporal n=1 radial recovery is then applied to the own mean:

```text
mean1 = Uz1+K1
Q1 = [(1-delta)Z*Uz1-(1+delta)Z*K1-(1-Z^2)*mean1_Z]/L
V1 = rho*epsilon*Q1.
```

The physical Q1 scale has one epsilon, and V1 has two. The Q/V error includes a further mean derivative bound, so ordinary Z1 recovery uses the newly admitted higher Z error control. These scales are Z-independent; no F0 derivative gauge is added to axial/pressure/radial errors. At rho=0 all seven components and their Z derivatives are exactly zero by the regular axis data; zero values are not counted as relative-width successes.

| Selected query at Z=371/1000 | Normalized passes | Physical passes |
| --- | ---: | ---: |
| Core axis | Exact zero data | Exact zero data |
| Core exit rho=4 | 14/14 | 14/14 |
| Actual first inlet | 14/14 | 14/14 |
| Corrected selected inner exit | 14/14 | 14/14 |

Each nonaxis count comprises value and ordinary Z1 for seven components. The target is total interval width divided by its absolute lower magnitude <=0.001. A crossing-zero reference is a failed width claim, even if its absolute enclosure is valid. The four previously accepted swirl rows are separate and use their own admitted F0*Lambda transformation.

## Artifacts and acceptance

Two new quartets under `experiments/root_st073/` have stems

```text
lei_ren_part1_paper_compliant_current_original_n1_leading_mean_epsilon_factor
lei_ren_part1_paper_compliant_current_original_n1_natural_component_enclosure
```

Each contains a source `.py`, result `.json.gz`, independent `_check.py` and source-bound `_check.json`. The checks cover executable original own-mean bindings, exact FTC and recovery identities, extra scaling and source rows, independent signed polynomial fixtures, genuine actual residual execution, original tube losses, directed log covers, physical scales, zero axis data, selected total width budgets and rejected scope mutations. The existing read-only reviewer was reused: **GPT-5.6 Luna / max**. No new agent was spawned.

```python
mean_owner=CurrentOriginalN1LeadingMeanEpsilonFactor(accepted_axial_epsilon_owner)
owner=CurrentOriginalN1NaturalComponentEnclosure(mean_owner,accepted_full_feedback_owner)
v=owner.report(owner.evaluate('.371','1'))
assert v['current_original_n1_natural_component_feedback_and_selected_axial_pressure_enclosures_installed']
assert v['passing_normalized_width_rows']==14
assert v['passing_physical_width_rows']==14
assert not v['n1_full_inner_interval_solution_certified']
```

The two parents must be accepted and share the same original axial source owner. Ancestor hashes, actual coordinate dispatch, packet identity and result scopes remain bound.

## Next agent task list

Prioritize the function/equation contract next. Do not repeat old conditioning scans or alter accepted ancestor files/receipts. Create a new leaf adapter and receipt so existing source closures remain usable.

- [x] **CURRENT-RP-N1-LEADING-OWN-MEAN-PSI-FACTOR** Original own core/collar mean, epsilon/pressure factors and leading Q radial/Z grades delivered.
- [x] **CURRENT-RP-N1-AXIAL-PRESSURE-NATURAL-SCALES** Exact extra normalization, genuine source scales and newly conditioned matrix majorants delivered.
- [x] **CURRENT-RP-N1-COMPONENT-WEIGHTED-FULL-TAIL** All six couplings admitted in the natural normalization, through ordinary Z2 on the explicit smaller tubes.
- [x] **CURRENT-RP-N1-SELECTED-AXIAL-PRESSURE-SOURCE-ENCLOSURES** Exact correlated pressure FTC and axis-forced axial seeds, all nonseed source bounds and full feedback errors delivered at the three nonaxis queries. No pressure-Z2 relative-width claim.
- [x] **CURRENT-RP-N1-SELECTED-RADIAL-RECOVERY** Own mean1, Q1 and V1 recovered from the common axial seed plus full errors; all six value/Z1 budgets pass at each selected nonaxis query. This closes selected recovery only.
- [ ] **CURRENT-RP-N1-FUNCTION-OWNER** Publish source-bound coefficient functions on every required Z and radius of the original common inner interval. Define the convergent principal/residual series as the function contract, keeping the exact sources and the same original core/collar family. Provide interval query semantics and explicitly distinguish exact functions, partial centers and error covers.
- [ ] **CURRENT-RP-N1-SIX-EQUATION-CONTRACTS** Derive regular-axis uniqueness and termwise six-row equation identities from the convergent series. Bound differentiated series on explicit tubes. Supply the full physical W=Snew*What conversion, including F0 derivatives exactly once. Point width targets alone do not certify the equations.
- [ ] **CURRENT-RP-N1-MIXED-RADIAL-JETS** Recover ordinary x/rho derivatives using the true six equations and their differentiated forms. Include actual chiSpatial derivative/inverse hb, source derivatives and matrix derivatives. Supply the mixed orders required by stress and every join, with a genuine axis Taylor limit. Do not infer arbitrary radial derivatives from the three endpoint samples.
- [ ] **CURRENT-RP-N1-COMMON-Z-AND-INNER-COVERAGE** Extend usable function/derivative enclosures over the required real Z interval and original common inner radii. Handle symmetric, vanishing and crossing-zero coefficients by absolute bounds. Preserve the corrected Rin dispatch; do not use canonical phase 1/2 as the selected original exit.
- [ ] **CURRENT-RP-N1-OWN-CUMULATIVE-MOMENTS** Integrate the completed n1 function's own cumulative radial moments. Keep actual core and collar histories, signed products and mixed Z jets. Leading meanPsi exposure is not a substitute for these n1 repair moments.
- [ ] **CURRENT-RP-N1-OWN-FIVE-MOMENT-REPAIR** Solve n1's original moment repair from those own cumulative moments, preserving divergence through the common potential. Establish terminal functional identities over Z, original repair support, pressure compatibility, joins and finite-energy conditions. Do not copy the leading repair.
- [ ] **CURRENT-RP-N1-CHART-EXTENSION-AND-HEAT-COMPATIBILITY** Extend the repaired n1 field through the completed transition/flatten/heat charts with the prescribed shared data. Check needed derivative joins, support and exact exterior heat compatibility. Keep the independent original pressure datum.
- [ ] **CURRENT-RP-N2-GENUINE-FORCING** Build true n-dependent n2 forcing using the completed leading and repaired n1 functions, their products and mixed derivatives. Bind the actual temporal weights and n-dependent recovery equations. No re-labeled n1 forcing or coordinate rescaling.
- [ ] **CURRENT-RP-N2-SOLUTION-AND-OWN-REPAIR** Solve the genuine n2 system on the same common interval and perform its own cumulative moment repair. Publish coefficient/tail bounds and equation identities before claiming two-level recursion.
- [ ] **CURRENT-RP-FINITE-ORDER-TIME-REMAINDER** Assemble repaired coefficients with the correct temporal powers; compute a finite-order equation remainder and the original regional stress/flat-remainder decomposition. Keep background stress separate from the full corrected NS residual.
- [ ] **CURRENT-RP-SMOOTH-TIME-SUM** Implement the original smooth sum/cutoff schedule with bounds for all required time and spatial derivatives and joins. Prove convergence and terminal flatness in the stated domain.
- [ ] **CURRENT-RP-ADMISSIBLE-STRESS-REGIONS** Compute cone margins over the actual recursive field in core, exit, matching, pulse/end, flatten and heat regions. Report regional extrema and norms using the same original data.
- [ ] **CURRENT-RP-OSCILLATORY-CORRECTION** Add the original two-family oscillatory pulses and mean corrections when the background prerequisites hold. Measure averaged quadratic stress and its cancellation of admissible background stress; preserve actual supports and independent forcing.
- [ ] **CURRENT-RP-FULL-RESIDUAL-AND-PHYSICAL-DIAGNOSTICS** After correction, validate Cartesian residual norms against the full <0.001 target. Measure actual time-dependent radial/axial scales, aspect ratio, velocity/vorticity growth, material winding and energy. Selected n1 component enclosures alone are not evidence of a completed finite-time blow-up construction.
