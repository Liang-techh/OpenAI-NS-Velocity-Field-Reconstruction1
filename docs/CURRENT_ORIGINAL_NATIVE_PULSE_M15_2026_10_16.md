# Original repaired weighted m15 through the selected pulse - 2026-10-16

`OriginalNativePulseM15Numeric.evaluate(Z,chart,coordinate)` evaluates the actual repaired weighted pressure moment m15 as factored ordinary C0/C1 function enclosures on the original entrance/main/gap/end pulse throughRv. It preserves the same original positive mu, selected native branch/N, true incoming histories, R*dR measure and finite R*Q_y endpoints. This extends the existing velocity and P1 work with the second pressure moment; it does not turn an m13/P1 packet into an m15 inlet. The full project remains **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS. Independent preRp numerical m15 joining, resolved signed centers/relative cancellation, global five-moment numerical admission, heat/stress/global-N, genuine temporal n>=2 and corrected NS remain open.

## Actual inlet from the repaired terminal theorem

The accepted complete `current_original_n1_pulse_numeric_bridge` encloses the actual full pulse increment in the Z-independent unit reference*A*Rv^2. It identifies the same true native source, original selected coefficients, all three charts, full incoming mean history, finite mu and physical endpoints. It is a pulse-only integral; inner/lower/preRp and the reserved repair have not been added to that component.

The separately accepted repaired terminal source proves the actual m15 value and its certified Z rows zero atRv. Its m15 handle is `terminal['m15']`, not the legacy P1 tail alias `terminal['m13']`. Since the compact n1 repair velocity/swirl is zero on this outer tail, the same pulse forcing remains. Therefore

```
actual repaired m15(Rp) = -complete original pulse increment.
```

This is the actual inlet function derived from the terminal condition. Its ordinary derivative rows come from the checked function enclosure; no interval endpoint is selected as the signed value. The extra minus is the terminal boundary orientation, not a second negation of the original weighted source. Full pulse bounds are used only for this complete boundary identity. Every main partial increment is computed separately and added to this inlet. An empty prefix equals the inlet exactly. A separate numerical integration of the repaired inner/preRp chain is still required for independent numerical seam admission and is explicitly not claimed here.

## Original weighted pressure identity

For L=1-delta*Z^2, d=1-Z^2, h=1-delta and the same fixed leading S used for P1,

```
m15_y = -R*P1_y,   P1_y=-R*S/2+Q_y+Q_yy
Delta m15 = S1/2-[R*Q_y]
```

The original S1 contains [R^2*Q]/L-I1_Q/L+hZ*I1_QZ/(2L), the finite quadratic boundaries, I1_Q_squared, I1_UQZ, I1_U_squared and I1_UUZ with their original coefficients. Here I1_f=integral R*f dR=integral R^2*f dt. It is not the unweighted I0 primitive and not an integral with R^2*dR measure.

## Main partial primitives and inherited mean

Let theta=Ep*exp(-b*t), b=1/2+mu, a0=1/2-mu, a1=3/2-mu, c=1-2mu, xi=mu*t and g=gp(xi). The original normalized Q recovery decomposes as

```
Q/theta = A+c0*exp(-a0*t)
A=bb*g+cc*k0
bb=2Z*ap/L
cc=-(hZ*ap+d*Dtheta(ap))/L
c0=-C1*(1+Z^2), Dtheta(c0)=0
k_r=integral_0^t exp(-r*v)*gp(xi-mu*v)dv
```

where Dtheta(f)=f_Z-2Z*f/(1+Z^2). The genuine generated I1 linear primitive is

```
I1_Q_generated/(R^2*theta) = bb*k_a1+cc*(k_a0-k_a1).
```

The kernel difference follows the actual two-rate convolution, a1-a0=1; it is not the previous I0 first-time kernel. Positive directed cells and nonzero omitted tails retain finite history. Dtheta(c0)=0 is a proved canonical inlet identity. Inherited linear terms cancel between the true S1 integral and its finite boundary, and inherited viscous terms cancel between the true RQ_y endpoints; they are not discarded by zeroing the mean.

The quadratic generated matrix integrates the actual (g,k_a0) products with rate c=1-2mu. The source identity gp(0)=0 and its pinned derivative bound give

```
k_a0 in g/a0 + [-89mu/a0^2,89mu/a0^2].
```

This is intersected with positivity and the monotone upper bound before the falling cutoff. Actual positive tail integrals remain. The matrix is a uniform partial-window enclosure, never a signed full-support cancellation. The inherited cross term is retained separately with coefficient c0*(A-bb*k_a0-cc*kv_a0)/2 and radial phase exp(a0*t). The pure inherited term is exactly

```
Delta m15_mean = -(Rp*Pstar*U0*C1)^2*t/4.
```

Its ordinary Z1 row is exactly zero from the canonical source identity. It is nonzero even when local axial velocity is inactive. The implementation retains the original negative logarithmic history rather than folding it into a cap or resetting it.

## End and inactive gap

The remaining end uses the same two original bumps and selected C5 coefficient functions. With H_j(s)=-integral_s^0 exp(a0*(v-s))*beta_j(v)dv, the full basis is beta1,beta2,H1,H2. The new positive-cell matrix covers the entire remaining interval, including inactive gaps. Thus it retains H_j*H_k contributions to Q^2 outside the local beta supports; the previous P1 beta-only quadratic helper cannot supply this weighted moment.

The a1-weighted beta mass is computed separately. The exact identity

```
integral_s^0 exp(a1*v)*H_j(v)dv
  =-W_a1_j(s)+exp(s)*W_a0_j(s)
```

provides both I1_Q and its ordinary Z derivative. Using the accepted actual repaired terminal value gives

```
m15(s)=-S1_remaining/2-R(s)*Q_y(s).
```

For the inactive gap, anchor at the actual end value at Rb=Rv*exp(-4), not atRv across active end bumps. With Q=C(Z)/R,

```
m15(R)=m15(Rb)-hZ*C_Z*(Rb-R)/(4L)+C^2*log(Rb/R)/4.
```

The projected Q coefficient includes its true Z derivative and the inherited means. Empty gap increments return the same end value before interval subtraction. Main xi=11 and gap xi=11 use the same actual remaining source. Beyond the last support all remaining weighted forcing is empty, giving the checked zero terminal value/rows without resetting earlier histories.

## Original units and numerical scope

Write F=Ep*exp(-(1/2+mu)*13/mu+logE). The end linear, quadratic and viscous units are Rv^2*F, Rv^2*F^2 and Rv*F. Their combined inverse-mu log coefficients are18.5,11 and5.5 respectively, with the original saddle/logmu/finite atoms preserved. The main-gap linear query reduces to (xi+5.5)/mu before boxing. Original logE is applied once and no hyper-exponential is materialized. The inlet reference*U0 factor is applied once in its fixed Z-independent p=2 unit; ordinary source derivatives are already in the coefficient rows.

This API reports the paper's weighted coefficient moment in its original R measure. It does not introduce an extra physical pressure time factor or claim a separately assembled physical-space moment diagnostic. Existing `OriginalNativePulsePressureVelocity` still supplies physical velocity/pressure through orders0/1. Signed point centers and useful cancellation of different exponential sectors remain unresolved.

## Focused evidence

The final focused run exited0 in210.985s. Saved-report inspection exited0 and verified1862 scientific input hashes and80 implementation hashes. The3,729,124-byte report SHA256 is `311c61f1bc38143498ba443acac266b7160001037ef4bc022e52d9faf52940f0`.

- 10 independent direct weighted-N1p value/ordinary-Z1 comparisons on main/end finite-mu cases. The reference integrates the original Q_y/Q_yy density with actual R weight, without the new S1 endpoint reduction or cell matrices. Split-support Gauss orders64/96 differ by at most1.319e-11 relative to1+absolute reference value; all reference values lie in the directed source intervals.
- 5 symbolic identities executing the actual new unit recipe and the accepted original p=2 recipe, including end linear/quadratic/viscous units and the combined gap-query inverse-mu atom;
- 2 positive off-support H_j*H_k tests, covering the history products that the former beta-only pressure helper cannot supply;
- 45 actual selected-parameter m15 packets,840 source/C1/scope conditions,12 terminal-derived inlet/chart seam conditions and18 nonzero inherited mean-log conditions;
- 36 source-certified terminal/empty-support zero conditions,9 weighted atom conditions and5 rejected invalid requests.

Use `python scripts/run_original_native_pulse_m15.py --inspect-report` to inspect the saved evidence without repeating the focused run. `--write-report` reruns only the new weighted-moment checks when implementation changes require new evidence.

Existing read-only worker `/root/n1_atlas_map`: GPT-5.6 Luna/max. It confirmed the actual terminal-derived inlet, two I1 rates, inherited cross/log terms, complete end history basis, backward/gap orientation and p=2 units. The independent reference initially omitted theta's radial-decay derivatives in Uy/Uyy; component comparison and an independent viscous endpoint identified the fixture error. The corrected reference restores these derivatives without changing the production formulas or relaxing their bounds. No new agent was spawned.

## Next construction work

The read-only source map for the next slice is now resolved. `current_original_n1_preRp_moment_functions.py:390-394` composes the signed active/compact/quiet functions; its compact source at295-330 retains the R*dR measure and finite RQ_y endpoints. Inner/lower/active directed components are separately available through `current_original_n1_full_five_numeric_bridge.py` and `current_original_n1_preRp_pressure_C1.py`. These enclosures do not yet provide a resolved signed numerical preRp center. The compact finite-N source exposes `perturbation_bounds.error_m15`, `endpoint_q_R_m15` and `signed_pressure_m15`; preserve these as errors/endpoints with their source ownership. The own-swirl term from `current_original_n1_reserved_own_fields.py:98-121` is -R*theta0*Utheta1*ds: consume `patch_histories['m15']` / `base_pressure_prefix['m15']` and the original p=2 unit. Neither the existing m13 swirl-prefix kernel nor its unit recipe can be reused unchanged.

1. Numerically integrate repaired preRp m15 from its actual inner/lower/preRp pressure components and compact repair response, including the exact compact/quiet K^2 terms and finite-N errors. Compare its actual source function with the terminal-derived inlet, not only overlap of unrelated interval boxes.
2. Assemble all five original n1 numerical histories at the required own-chart joins, with signed functions and required derivative/mixed bounds. Keep global five-moment numerical admission open until this evidence exists.
3. Continue the original flatten/heat field with compatible analytic pressure, then energy/radial tail, common global N, admissible stress and flat remainder.
4. Build independent temporal n>=2 recursion/smooth sum, oscillatory stress correction and independent corrected Cartesian NS validation.
