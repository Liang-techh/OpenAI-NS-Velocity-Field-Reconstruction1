# Original n=1 full weighted feedback bound and selected coupled swirl

The complete residual matrix now has uniform complex bounds after source normalization and principal swirl resummation. A new all-word Volterra/Cauchy tail covers every remaining coupling among all six rows. That tail is consumed by the actual first-collar swirl enclosures. At Z=371/1000 and corrected Rin, all four normalized and all four physical swirl value/Z1 width targets pass with **all remaining feedback included**. This advances beyond the previous initial-integral-only result. The full reconstruction goal remains **ACTIVE / INCOMPLETE**.

Read this checkpoint before [the axial epsilon factor](CURRENT_ORIGINAL_N1_AXIAL_EPSILON_FACTOR_2026_10_11.md) and [the correlated initial source](CURRENT_ORIGINAL_N1_CORRELATED_INNER_SOURCE_2026_10_11.md). These new bounds do not certify axial/pressure/mean/Q relative accuracy, the full n1 function contracts, its own five-moment repair, or n2 recursion. The old wide 0/16 coupled-prefix result belongs to the historical route; it is not the new four-row coupled swirl result.

## Exact normalized residual matrix

The unknown normalization remains

```text
W=S What, S=diag(F0*Lambda,1,1,1,F0*Lambda,1), epsilon=Lambda^-1.
B0star=S^-1*B0*S+S^-1*B1*S_Z, B1star=S^-1*B1*S.
```

Known source derivatives of F=F0*Phi are cancelled symbolically before interval boxing. The unknown-weight derivative is added once. Akin has entries (1,5)=1, (2,6)=1 and (3,6)=-1; Achi has (5,1)=-2chi. Their supports are disjoint, so the principal infinity row norm is at most 2 on the original chi tube. Both are removed from the residual.

Put d=1-Z^2, L=1-delta*Z^2, x=sqrt(rho), H=Z*(1-delta)/2+d*Uz, E=Z*Uz-1/2, B=Q+(1-2Z*Uz)/L, Cphi=Phi+rho*Phi_r and m=d*Gprime/L. The exact residual includes

```text
B0res(4,1)=4*x*F0^2*Phi
B0res(5,1)=2*epsilon*(Q+(-2+delta)*E/L)+remaining_gauge
B0res(5,2)=epsilon^2*(2*Ta(Phi)-2*(-1+delta)*Z*Cphi/L)
             -2*epsilon*m*Phi
B0res(5,3)=-2*epsilon^2*(1+delta)*Z*Cphi/L
B0res(6,1)=-8*epsilon*F0^2*Z*rho*Phi/L.
```

The last row has one epsilon **after** normalization because the original executable B0(6,1) has two epsilon shifts: one inside its tuple value and one applied by the enclosing loop. The checker replays that executable source and all 16 nonzero normalized residual/B1 entries, rather than relying on a formula copied from a line fragment. The other exact entries, their actual normalized Phi callback, and the original Uz/Q inputs are installed in the new module.

Uniform matrix majorants are on eta/2. Original core/collar Phi, Uz and own mean bounds are retained, with Cauchy gaps eta/2 for their required Z derivatives and the exact radial recovery for Q. The (5,1) remaining gauge uses the admitted epsilon*Psi identity. F0-squared back-coupling uses the original uniform complex Fmax bound. Norm sectors explicitly label Fmax as a bound proxy; it is never selected as the source amplitude. Pressure, positive source width and Lambda sectors remain present.

## Uniform initial source bound and full tail

The new normalized g5 bound uses the already bound identity

```text
gstar5=-2*m^2*Phi+2*epsilon*cross-2*epsilon^2*Tb(Ta(Phi)).
```

m is bounded on eta/2 and its derivative on eta/4. Normalized Phi and the genuine axial/pressure forcing supply the other source rows. The core integration length is 2. Actual collar integration uses dx=hb*x*ds/2 with phase length at most 1/2; all inverse hb powers cancel before taking the norm. The L1 source norm sums rows 4, 5 and 6, because the principal evolution mixes rows. Multiplying by exp(2a) gives an initial source M0 on eta/4, where a=sqrt(4.1).

All scalar covers are directed. A sector may receive the small cover 1e-240 only after its original source log is proved no larger than log(1e-240). This changes an upper-bound rounding floor; it does not choose a larger physical width/amplitude or alter a source parameter.

For k further iterations, allocate k equal Cauchy gaps across loss ell=eta/8. Ordered real radial integration contributes a^k/k!. The principal Green factors telescope to exp(2a). Thus

```text
||V^k W0|| <= M0*exp(2a)*a^k*(C0+k*C1/ell)^k/k!
           <= M0*exp(2a)*Q^k
Q=e*a*(C0+C1/ell)<1.

tail after retained terms 0..N <= M0*exp(2a)*Q^(N+1)/(1-Q).
```

Every derivative word is covered; the old ceil(k/2) count is not reused for the resummed propagator. Nested Cauchy bounds apply to the entire inner function and therefore include derivatives of earlier Green factors and coefficients. The value tail lies on eta/8. A separate final gap eta/16 bounds ordinary Z1 on eta/16. The original outer eta is unchanged.

The new default retains only W0. All later feedback is enclosed, without assigning it the value zero:

| Uniform quantity | Conservative upper bound |
| --- | ---: |
| M0 on eta/4 | 3.554692e56 |
| C0 residual norm on eta/2 | 7.000001e-240 |
| C1 norm on eta/2 | 2.000001e-240 |
| Q | 1.839753e-204 |
| All feedback value remainder on eta/8 | 3.752493e-146 |
| All feedback ordinary Z1 remainder on eta/16 | 1.254275e-110 |

The numeric table rounds the directed interval outputs upward. These are absolute bounds in the stated six-component normalization. They are not full physical NS residuals, proof of blow-up or relative error bounds for every velocity/pressure component.

## Selected actual coupled swirl enclosures

The genuine signed core-plus-actual-collar source integrals remain the centers for F1 and d_x_F1. Their nonlinear, quadrature, kernel and homogeneous-transport errors are retained. The uniform feedback remainder is added to both rows, including its Z1 error. Physical errors use

```text
|delta W| <= |S|*|delta What|
|delta W_Z| <= |S|*|delta What_Z|+|S_Z|*|delta What|,
S_Z=-Lambda*Gprime*S.
```

The ordinary derivative gauge is applied once. The new budget rows say `remaining_B0_B1_feedback_included=true`. At the exact first inlet and corrected selected exit, four normalized and four physical F1/d_x_F1 value/Z1 rows pass the total-width/absolute-lower target 0.001. Only these selected swirl rows have this new component accuracy admission; no whole-Z relative-width claim is made.

Two quartets under `experiments/root_st073/` have stems

```text
lei_ren_part1_paper_compliant_current_original_n1_weighted_feedback_bound
lei_ren_part1_paper_compliant_current_original_n1_coupled_swirl_enclosure
```

Each contains a `.py`, `.json.gz`, `_check.py` and `_check.json`. The checks cover executable matrix normalization, principal supports, exact simplex/factorial premises, genuine actual Phi execution, original tube losses, every directed log cover, physical error transformation, actual coupled width rows, and packet/scope guards. The read-only reviewer is GPT-5.6 Luna / max; no new child was spawned.

```python
bound_owner=CurrentOriginalN1WeightedFeedbackBound(accepted_axial_epsilon_owner)
swirl_owner=CurrentOriginalN1CoupledSwirlEnclosure(bound_owner)
v=swirl_owner.report(swirl_owner.evaluate('.371','1'))
assert v['current_original_n1_selected_actual_swirl_enclosures_include_full_residual_feedback']
assert v['passing_normalized_coupled_width_rows']==4
assert v['passing_physical_coupled_width_rows']==4
assert not v['n1_full_inner_interval_solution_certified']
```

## Detailed next tasks

- [x] **CURRENT-RP-N1-WEIGHTED-RESIDUAL-MAJORANT** Exact normalized residual and B1 entries, genuine normalized Phi callback and uniform C0/C1 bounds delivered. Preserve this accepted source closure.
- [x] **CURRENT-RP-N1-RESUMMED-PROPAGATOR-FULL-TAIL** Uniform normalized source M0 and all-word tail delivered with the corrected source/value/Z1 tube schedule. The old nonadjacent-word count is superseded for this normalization.
- [x] **CURRENT-RP-N1-SELECTED-COUPLED-SWIRL-BOUND** The two actual swirl rows and ordinary Z1 include the full feedback tail and pass the selected 0.001 width target. This is a component enclosure, not the full solved n1 field.
- [ ] **CURRENT-RP-N1-AXIAL-PRESSURE-NATURAL-SCALES** Derive source scales for Uz1/K1/P1/d_x_Uz1 from their defining forcing before boxing. The current shared absolute feedback bound is adequate for selected swirl, but can swamp an exponentially smaller axial/pressure coefficient. Review candidate additional weights Uz1/K1/d_x_Uz1~epsilon and P1~epsilon*Pstar^2; derive their validity rather than assuming them. Retain hydro, pressure, swirl and own-mean contributions separately, including degenerate/crossing-zero cases.
- [ ] **CURRENT-RP-N1-LEADING-OWN-MEAN-PSI-FACTOR** Construct genuine leading mean-uaxis=epsilon*meanPsi from its own core/collar FTC history. Recover Q and Q_r/Q_rr from this shared field, exposing pressure/epsilon sectors before subtraction. The actual mean is not the comparison mean. Include inlet equality, radial moments and ordinary Z derivatives. This prevents the rough (Uz-mean)/rho bound from losing the small factor needed for pressure forcing.
- [ ] **CURRENT-RP-N1-AXIAL-PRESSURE-SIGNED-INTEGRALS** Compute correlated core plus actual collar signed g4/g6 integrals, their Z derivatives and full source/quadrature/kernel errors on the corrected interval. Keep the independent original pressure datum. Use original source scales for each component's reference norm, and handle a reference crossing zero as a width failure rather than asserting accuracy.
- [ ] **CURRENT-RP-N1-COMPONENT-WEIGHTED-FULL-TAIL** Conjugate the whole six-row system with any newly justified weights, rebinding each entry and source. Recompute uniform M0/C0/C1 on explicit tubes. A bound after multiplying by a component's scale is required before combining it with the signed component interval. Do not transplant the swirl normalized tail directly as an axial/pressure relative error claim.
- [ ] **CURRENT-RP-N1-RADIAL-RECOVERY-AND-MIXED-JETS** Combine the admitted coupled Uz1/K1 mean with the exact n1 Q1/V1 recovery. Supply the ordinary Z and radial derivatives needed by stress/matching and verify divergence through the common potential. Do not fit an independent radial coefficient.
- [ ] **CURRENT-RP-N1-FUNCTION-AND-EQUATION-CONTRACTS** Publish the actual common-inner coefficient as a source-bound function owner for all Z and all required real radii. Prove regular axis data, the full six equations, their mixed derivatives and common analytic extension from the convergent series. Tail norm admission and two selected endpoint widths alone do not close these contracts.
- [ ] **CURRENT-RP-N1-OWN-FIVE-MOMENT-REPAIR** Use the solved n1 field's own five cumulative radial moments to perform its divergence-preserving repair. Verify exact joins/supports/pressure/energy over the completed leading charts. No reuse of the leading repair as a substitute for n1's own history.
- [ ] **CURRENT-RP-N2-GENUINE-SOURCE-AND-OWN-REPAIR** Build the true n-dependent n2 forcing from the completed leading and repaired n1 fields, then solve/repair it on the same common interval. Merely adding n=2 to a label or reusing an n1 forcing packet does not count as scale recursion.
- [ ] **CURRENT-RP-TIME-SUM-AND-REMAINDER** Continue finite-order time remainder and the smooth sum with explicit derivative/convergence bounds through the chart interfaces. Preserve the original independent pressure and exact exterior heat field.
- [ ] **CURRENT-RP-FULL-PHYSICAL-DIAGNOSTICS** Measure core radial width, relative axial slenderness, velocity/vorticity growth and material winding on the actual time-dependent field. Then close the regional stress cone, flat remainder, oscillatory correction and its averaged quadratic stress cancellation, and test the full corrected NS residual. The current n1 bounds do not by themselves prove a finite-time singularity.
