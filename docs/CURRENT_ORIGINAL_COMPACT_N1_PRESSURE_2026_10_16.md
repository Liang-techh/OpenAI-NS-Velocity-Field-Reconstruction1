# Original compact n1 pressure numerical assembly - 2026-10-16

`OriginalCompactN1PressureVelocity.evaluate(Z,x,log_tau="-2",azimuth=".25")` now exports the original compact physical velocity and pressure with temporal orders0/1 on x=R/Rc in[1,2]. The new pressure reader consumes the same checked signed source functions and their ordinary C0/C1 enclosures. It adds the complete incoming core/lower/active pressure, fixed leading compact N1p prefix and actual local swirl-pressure integral. The project remains **ACTIVE / INCOMPLETE**. Tasks14/15 remain IN PROGRESS. Signed point centers, useful relative pressure cancellation, long quiet continuation and genuine temporal n>=2 remain open.

## Actual pressure partition and source ownership

Let s=log(R/Rc). The authoritative source is `current_original_n1_reserved_own_fields.pressure_prefix_in_band` and its same-graph cumulative m13 function:

```
P1(s) = P1_inner + P1_lower + P1_active
        + integral_Rc^(Rc*exp(s)) N1p_leading dR
        + integral_0^s Utheta0(u)*Utheta1(u) du.
```

The inner and lower terms use the checked full-five numerical bridge in the original fixed m13 pressure unit. Active pressure uses the checked pre-Rp bridge's active piece only; compact, quiet and pulse pieces are not mistakenly added as incoming pressure. The core pressure datum, true lower source errors and all17 active-chart memories persist once. The fixed leading forcing is retained, rather than computed from the repaired n1 field. The original P0 datum remains once in the leading physical pressure.

The new output is a signed exponential sum with separate directed coefficient and log intervals, consistent with the existing physical velocity format. Magnitude caps enclose unknown signed functions; their midpoint is a numerical proxy and never defines the exact function. Hyper-exponential physical radii and amplitudes are not materialized.

## Nonzero leading pressure and finite-N partial errors

With the original k=K/Rc, the compact baseline has U0=0, q0=-k/x and exactly Z-independent k. The radial-free pressure identity gives

```
P1_compact_leading(s) = Rc*k^2/4*(1-exp(-s))
D_Z P1_compact_leading = 0.
```

The finite viscous endpoints q+R*q_R vanish for this baseline. At s=log2 the nonzero value is Rc*k^2/8, exactly the checked full-band baseline. This positive pressure is not discarded when the local axial velocity is zero.

The actual finite-N perturbation is separately enclosed. C0 uses `repair_band_pressure_source.perturbation_bounds.error_m13`; C1 uses `preRp_pressure_C1.source_data.original_target_and_control_C3_update.compact[*].normalized_first_Z_upper.m13`. These are termwise absolute endpoint and integral bounds uniform for any upper x in[1,2], including the lower endpoint. The boundary factors3/5 and positive integral-width majorants do not rely on cancellation at the full endpoint. A full-band signed value is not scaled into a partial-prefix value.

The original `preRp_C1_inputs.scale_ratios.m13` is log(Rc/fixed_pressure_unit), with shared Rc canceled before boxing. Both finite-N normalized derivative bounds are proved below-1000 before attaching the directed positive exp(-1000) cover. At x=1 the same-endpoint prefix is exactly zero. This finite-N cover is not differentiated as a defining source.

## Actual supported swirl-pressure primitive and C1 product

Use the same ell=log2/400, centers c_j=j*log2/40 and original positive normalizer Nhalf=integral beta(t)*exp(ell*t/2)dt. With alpha=1/2+mu,

```
shape_j(s) = beta((s-c_j)/ell)/(ell*exp(c_j/2)*Nhalf)
I_j(s) = integral_0^s exp(-alpha*u)*shape_j(u)du
       = exp(-(alpha+1/2)*c_j)/Nhalf
         * integral_-1^clip((s-c_j)/ell) beta(t)*exp(-alpha*ell*t)dt.
```

Only the three actual angular supports j=13,14,15 enter pressure. Their disjointness from the old leading repair supports makes Utheta0=Pstar*A_rc*exp(-alpha*s) exact on these integrands, while nonlocal histories remain. With the fixed original A_rc scale and normalized beta coefficients, the coefficient of each pressure term is A_rc*beta_hat*I_j, and its ordinary Z derivative is (A_rc_Z*beta_hat+A_rc*beta_hat_Z)*I_j. Source amplitude derivatives are retained; fixed physical units have exactly zero Z derivative. No extra sqrt(2R), inverse Lambda or physical dR factor is inserted into theta0*Utheta1 ds.

The original positive mu is retained. A directed degree3 midpoint Taylor integral plus the full-cell fourth-derivative/4! remainder encloses each weighted beta primitive, including its flat endpoints. The before-support prefix is exactly zero; after support its cumulative value persists. Bump histories are never reset.

## Physical pressure units and temporal weight

The fixed m13 unit is exactly Phi*A*Rv, Phi=10.99*sqrt(2*pi)*(exp5-exp3). Its source recipe is decomposed into separate logRc and scalar terms, with the common inverse-mu powers combined before interval boxing:

```
logUnit13 = logRc + logPhi + logP + logUp -15
            +11/(2mu) -3*(4mu)^(-2/3)
            -(184/3)*logmu +(2/3)*log4 -(1/2)*log6.
```

The swirl-pressure term keeps the existing Utheta1 log unit plus logP and the source-fixed A_rc log scale. The exact leading compact pressure has logRc+2*logk-log4 as its own unit. Separate units preserve terms that would disappear in a prematurely boxed sum. An independent symbolic lift proves equality with the original `pressure_unit_log(...,power=1)`.

Physical pressure is the existing temporal-order-zero pressure plus lambda^(-2)*P1. It receives no additional Lambda^-2 conversion. Ordinary and source-scaled times retain the same delta cancellation as the physical velocity. The public pressure is a factored enclosure, not a high-relative-accuracy resolved signed sum.

## Focused evidence

The focused run exits0 in101.451s and the saved-report inspection exits0. It exports33 physical packets at3 axial positions, before/in/after all three angular supports, partial angular bumps, Rc and2Rc, with ordinary and source-scaled times. Checks cover15 independent weighted-beta defining integrals, independent original m13-unit/baseline/viscous-endpoint algebra,120 source/history conditions,80 actual partial finite-N bounds,60 nonzero-baseline conditions,102 support/positive-integral conditions,36 incoming-pressure conditions,210 physical pressure order/unit conditions and24 scaled-time conditions, plus12 exact inlet-alias conditions. Maximum tested raw integral width is2.811897439742031e-8.

Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_compact_n1_pressure_numeric.json`; 30,381,627 bytes; SHA256 `f3400c58ab6267a96db7ab0787a953b9147a4eda33c1f79a57231dfabbc3ed3b`; 1825 scientific input pins and 73 implementation pins. Accepted scientific ancestor and implementation bytes remain unchanged.

The existing read-only `/root/n1_atlas_map`, GPT-5.6 Luna/max, confirmed the original partial-pressure constructor, nonzero baseline, termwise uniform finite-N partial bounds, lower endpoints, physical normalization and incoming-history ownership. Its implementation review found no substantive underbound or unit/ownership defect. An exact rational endpoint comparison also preserves the zero prefix for decimal aliases1.0 and1.000e0. No new worker was spawned.

## Next construction work and remaining limits

1. Continue the original nonzero-mu quiet background from2Rc toRp. Use stable analytic exponential/exprel own-rate formulas for length Tw-2-log2. Preserve the full leading m/h/k/e/p histories and P0, and the source-admitted zero n1 velocity tail. The short compact polynomial-moment series is not a long-interval oracle.
2. Continue P1 through the quiet band and original pulse toRv. Use the same original fixed N1p, complete source histories and pressure gauge. Velocity/mass closure at2Rc does not imply P1=0 there. A backward source representation may improve cancellation only with the actual full terminal theorem and full remaining forcing.
3. Resolve actual signed finite-N/repair coefficient values and pressure cancellation where viable. Current C0/C1 enclosures do not supply a useful combined relative pressure error, an independent signed terminal residual, Q1_Z, pressure radial derivatives or numerical C3 mixed norms.
4. Close global physical/high-order joins and common-N admission, axis/tail energy, compatible analytic preheat/exact heat exterior, admissible stress and independently measured flat remainder.
5. Build genuine independent n>=2 coefficient recursion and smooth sum, then the two-family oscillatory correction and independent corrected Cartesian NS validation. The new pressure layer does not close these gates.
