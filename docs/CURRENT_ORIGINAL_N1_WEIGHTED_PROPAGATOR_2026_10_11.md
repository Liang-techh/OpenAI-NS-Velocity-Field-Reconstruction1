# Original n=1 weighted propagator and actual source integrals

The original n=1 regular system now has an executable F0*Lambda-normalized operator, an exact leading gauge potential in its radial propagator, and actual source integrals on the corrected core/first-collar interval. The four swirl kernel-truncation budgets retain their source sectors and are finite at the actual corrected exit. The full reconstruction goal remains **ACTIVE / INCOMPLETE**.

Read this checkpoint before [the local-dispatch/component-budget checkpoint](CURRENT_ORIGINAL_N1_LOCAL_DISPATCH_COMPONENT_BUDGET_2026_10_11.md). That earlier checkpoint's 0/16 full relative-width result remains authoritative for the old first coupled prefix. This new delivery is a differently preconditioned initial source integral; it is not a completed coupled n=1 solution and does not replace the old coupled prefix.

## Source normalization and exact gauge

Use the six original regular fields and the source weight

```text
W=(F1,Uz1,K1,P1,d_x_F1,d_x_Uz1)
S=diag(F0*Lambda,1,1,1,F0*Lambda,1)
W=S*What.
```

Only the swirl pair is normalized by F0*Lambda. Axial and pressure fields keep their separate zero/even-F0 sectors. F0, original hb and their astronomical scales remain formal source factors.

The conjugated operator is

```text
B0star = S^-1 B0 S + S^-1 B1 S_Z
B1star = S^-1 B1 S
gstar  = S^-1 g.
```

Ordinary Z jets of the transformed coefficients include the derivatives of S. `relative_amplitude` solves the exact Taylor recurrence `R'=-p*Lambda*Gprime*R`, `R(0)=1`, for the relative F0^p factor. Lambda powers remain algebraic throughout. Existing ordinary leading-source F derivatives are retained; they are not given a second gauge term. The extra term concerns the derivative of the unknown weight only.

The original source identity is

```text
Gprime = L*Haxis/(Haxis^2+sigma^2)
chi = Haxis^2/(Haxis^2+sigma^2)
epsilon*Lambda = 1
epsilon*(2*Hactual/L)*(-Lambda*Gprime)
    = -2chi - 2*(Hactual-Haxis)*Gprime/L.
```

The leading **-2chi** contribution is not a small epsilon perturbation. It is now included in the propagator together with the three original kinematic rows. The remaining matrix is explicitly `B0star-Akin-Achi`; the kinematic rows are removed from it to prevent counting them twice. The residual `(Hactual-Haxis)` still needs a tight original-source bound before a full weighted feedback-tail admission.

## Original analytic tube retained

For the original polynomial

```text
Haxis(Z)=-4Z^3-jZ^2+(4+(1-delta)/2)Z+j,
```

the original complex tube gives a directed bound b on |Im Haxis|. The accepted source satisfies `sigma^2-2b^2>0`. Consequently `Re(Haxis^2+sigma^2)>0` and

```text
|Haxis|^2 <= Re(Haxis^2+sigma^2)
|chi| <= 1.
```

The original common complex-Z radius is preserved. It is not replaced by an epsilon-sized tube or selected from a rounded amplitude enclosure. This bounds the potential propagator uniformly on that original tube without expanding F0.

## Resummed radial source solve

The normalized swirl block of the initial source solve is

```text
f_xx + 3*f_x/x + 2*chi*f = gstar5
f(0)=f_x(0)=0.
```

Let `K0(x,s)=s*(1-s^2/x^2)/2` and `K_j=K0` composed radially j additional times. The signed propagator is `sum_j>=0 (-2chi)^j K_j`. The implementation constructs the exact rational power/log kernel coefficients, retains indices 0..12, and bounds all later indices. It also computes the radial derivative kernel and first ordinary Z derivative directly through chi_Z.

For real `0<=s<=x<=a`, positivity of the bare kernels yields

```text
K_j <= a^(2j+1) / (2*4^j*j!*(j+1)!)
d_x K_j <= j*a^(2j) / (4^j*j!*(j+1)!), j>=1.
```

Including the factor 2^j from the potential gives the admitted factorial tail with r=a^2/2. Separate ratios cover the value, the j-weighted derivative and the j^2-weighted mixed derivative. The Z derivative tail uses the actual chi_Z jet, so it does not incur an unnecessary Cauchy loss.

The actual source delivery integrates all six gstar fields over 16 original core cells and four actual first-collar cells at Z=371/1000 and corrected Rin. The first-collar dx measure cancels inverse hb before enclosing; actual own source/moment data and the same pressure graph are retained. The normalized initial fields and their physical ordinary Z1 jets are returned. The core-axis initial trace is exactly zero. Outer endpoint intervals are rejected; internal integration cells remain interval enclosures of entire functions.

Kernel tails multiply the absolute source jets **per monomial sector**, then transform back with the proper weight and derivative gauge. No scalar zero-F0 tail is inserted into the swirl result.

At the corrected exit, the kernel-tail bound divided by each original source monomial and directed coefficient norm is:

| Normalized field | Z0 | Z1 |
|---|---:|---:|
| F1/(F0*Lambda) | 2.286e-17 | 4.572e-17 |
| d_x_F1/(F0*Lambda) | 3.590e-16 | 7.179e-16 |

These rounded values describe **potential-kernel truncation only**. The normalization is a source bound, not a selected field value or a measured relative error of the exact solution. The signed initial prefix still has a zero-crossing enclosure in these four budgets, so its relative field width is not certified. Source/whole-cell quadrature uncertainty and remaining B0star/B1star feedback are still unfinished.

## Files, public API and checks

The new quartet under `experiments/root_st073/` is

```text
lei_ren_part1_paper_compliant_current_original_n1_weighted_propagator.py
lei_ren_part1_paper_compliant_current_original_n1_weighted_propagator.json.gz
lei_ren_part1_paper_compliant_current_original_n1_weighted_propagator_check.py
lei_ren_part1_paper_compliant_current_original_n1_weighted_propagator_check.json.
```

```python
owner = CurrentOriginalN1WeightedPropagator(accepted_component_budget_owner)
view = owner.report(owner.evaluate('.371', '1', 'first'))
assert view['current_original_n1_weighted_potential_propagator_and_source_integrals_computed']
assert view['remaining_normalized_B0_B1_feedback_not_solved']
assert view['full_sector_preserving_coupled_tail_still_open']
assert not view['n1_full_inner_interval_solution_certified']
```

Acceptance checks include 234 exact rational moments derived independently from the radial ODE, 33 independent symbolic operator/forcing derivative rows, 24 high-precision Bessel Green-function consistency observations, the actual three-row Haxis*Gprime/L=chi identity, epsilon*Lambda=1, the 20 original source cells, axis zeros, removed kinematic feedback rows and ownership/mutation/endpoint guards. The Bessel observations are consistency checks; the rational identities and analytic kernel/tube bounds carry the tail admission.

Checked public construction/reporting returns the new local gate true. Global support, solved/full-inner, repair and n=2 gates stay false. The single reused read-only reviewer remains **GPT-5.6 Luna / max**; no new child was spawned.

## Ordered continuation tasks

- [x] **CURRENT-RP-N1-WEIGHTED-SOURCE-OPERATOR** Transform original B0/B1/g with the correct F0*Lambda block weights and genuine ordinary jets. Keep pressure and axial sectors separate and add the unknown gauge once.
- [x] **CURRENT-RP-N1-LEADING-CHI-PROPAGATOR** Put -2chi and Akin into the regular radial propagator, retain the original complex tube, and publish exact signed kernels plus sector-preserving value/Z1 truncation bounds.
- [x] **CURRENT-RP-N1-WEIGHTED-ACTUAL-SOURCE-INTEGRALS** Compute the new initial integral at the corrected exit and exact zero axis, including real actual collar data, original width cancellation, physical Z1 transformation and the four finite kernel-tail diagnostics.
- [ ] **CURRENT-RP-N1-ANALYTIC-DOMINANT-CORE-SOURCE** Extract the leading normalized theta forcing algebraically before boxing: `gstar5=-2*(d*Gprime/L)^2*Phi+lower Lambda sectors`. The original explicit model is `Phi_model=B(chi*rho/2)` and satisfies `(d_xx+3*d_x/x+2chi)Phi_model=0`; differentiating with respect to chi gives `Lchi(partial_chi Phi_model)=-2Phi_model`. Independently check these identities and compute the dominant core source integral as `(d*Gprime/L)^2*partial_chi Phi_model`, including the radial/Z derivatives. Enclose the original nonlinear Phi correction and lower forcing sectors separately through the same propagator. Keep the actual prescribed-shear collar integral separate; the analytic core extension is not the actual collar. This offers a correlated source integral without substituting the model for the genuine field.
- [ ] **CURRENT-RP-N1-WEIGHTED-RESIDUAL-MAJORANT** Recompute normalized forcing and residual-matrix majorants from original Phi/Uz/mean/pressure source bounds. Bound Hactual-Haxis as the original small correction, rather than subtracting independently boxed nearly equal velocities. Use explicit epsilon*Psi in the core and the actual integrated collar correction. Preserve width cancellation and F0^2 pressure/back-coupling sectors. Do not reuse the old unweighted Cp/C1 or scalar M0 after division by F0.
- [ ] **CURRENT-RP-N1-WEIGHTED-COUPLED-FEEDBACK** Integrate `(B0star-Akin-Achi)*What+B1star*d_Z What` with the new propagator and genuine initial Z derivatives. Compute the first actual residual-feedback term, retain derivative losses, and derive a compatible full weighted series tail. Transform errors back with the physical weight/gauge. Only that tail can replace the old unweighted swirl bound.
- [ ] **CURRENT-RP-N1-SIGNED-INTEGRAL-REFINEMENT** Attribute the zero-crossing weighted prefix to original source jets and kernel/whole-cell quadrature contributions. Preserve correlations in the exact power/log kernel and original source model. Consider extracting/integrating the dominant analytic Phi-model source and enclosing its genuine nonlinear correction separately; do not replace the original field by the model. Refine only material radial cells, report interval widths and distinguish them from kernel truncation.
- [ ] **CURRENT-RP-N1-AXIAL-PRESSURE-CORRELATION** Tighten the inherited 13.3%/20.4%/39.8% axial/pressure widths and mean/Q cancellation using original FTC/shared-source identities. Maintain pressure datum, incoming means and exact temporal n=1 radial recovery.
- [ ] **CURRENT-RP-N1-SUPPORT-EQUATIONS-WHOLE-Z** Complete the global same-source support dispatch, whole-Z/endpoints and derivative coverage, and independently assemble solved n=1 equations. Neither this initial integral nor an actual first feedback term alone proves the solved coefficient.
- [ ] **CURRENT-RP-N1-OWN-REPAIR-THEN-N2** After solved inner n=1, perform its own five-moment annular repair through a streamfunction/vector potential, retain joins/pressure/energy, and then build and solve actual n=2 sources with the true n-dependent recovery on the same interval.
- [ ] **CURRENT-RP-RECURSIVE-SUM-AND-PHYSICAL-DIAGNOSTICS** Continue finite-order/time remainder, smooth sum, core-width/aspect ratio, swirl/vorticity growth, material winding and energy measurements. Regional cone, flat remainder, oscillatory stress cancellation and full corrected NS validation remain required.

The next advance is the actual residual feedback and its compatible weighted majorant, alongside tighter signed source integration. This checkpoint does not establish completed scale recursion or full NS reconstruction.
