# Original O2 factored point inputs and true radius phase

Checked source: [d762ade6](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/d762ade6d02e9dc295ffe6ab7c570d5a79ac27ff). Predecessor: [original parameter frame and normalized pressure points](CURRENT_ORIGINAL_O2_SOURCE_PRESSURE_POINTS_2026_10_07.md). Both new checkers PASS. Git-index dependency audits cover 619 and 1049 hashes. Existing read-only reviewer: **GPT-5.6 Luna / max**; no new workers or ancestor reconstructions.

Two new production services now connect the original O2 slope chart to numerical source evaluation. Complete E,V,a,b,p1,p2 value/Z rows carry genuine source-derived finite coefficients, exact correlated factors and directed radial/pressure errors. A separate lightweight service evaluates the original spatial phase frac(N log(R/r_minus)) with directed periodic covers, including the positive microscopic origin offset. Native conditioned loop inversion, A/B evaluation, a signed density integral, all-chart oracle, controls, common N and actual coefficient recursion remain open.

## Complete O2 point inputs

`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_O2_factored_point_inputs.py/.json` and `_check.py/.json` implement `OriginalO2FactoredPointInputs(dps=50,cells=4096).evaluate(y=...,Z=...)` on exact y in[0,1], Z in[-1,1]. Output `inputs` contains value/ordinary-Z row pairs for E,V,a,b,p1,p2. The source coordinate is y=log(R/Rref), not log(R/110). Each row is a sum of finite coefficients times the same original

```text
R^r * Pstar^p * delta^d * L^ell, L=1-delta*Z².
```

The finite radial functions are the defining original J_sigma and three integrals M0,M1,M2:

```text
f=exp(y/10-(3/5)J_sigma(y))
H=(5/8+M0)*exp(-3y/2)
D=(5/12+M2/2)*exp(-y)
P=5/2+M1/2
E=f/(1+Z²), V=4Z/Pstar
a=4/5+(6/5)sigma(y), b=0.
```

Actual quadrature supplies finite point coefficients. Independent 4096-cell directed integrals enclose the true radial inputs. The checked pressure defining constant alpha and its independent directed enclosure are reused as a Z-independent source constant, not as a cover/cap point. Full original inertial formulas retain both axial Pstar^-1 and Pstar^+1 sectors, including the positive energy and absolute pressure terms. Ordinary Z rows include P0_ZZ and all delta/L derivative terms. The selected radius, Pstar, delta and L are never numerically expanded or reset.

Every term carries a finite coefficient absolute-error upper bound and, if pressure-dependent, a separate logarithmic late-source error upper bound. Consumers must propagate

```text
sum |R^r Pstar^p delta^d L^ell| *
    (finite_coefficient_error + exp(late_pressure_error_log_upper)).
```

An explicit zero late record means exactly zero; other tiny positive errors must remain positive. Coefficients are linear in P0/P0_Z/P0_ZZ, so the pressure sensitivity propagation is exact in those inputs. These errors apply to the full restored original quantities only after their correlated factors are applied; normalized errors alone do not bound the physical residual.

**Consumer rule:** the order-1 row is already the completed ordinary Z derivative. Its L^-2 factor belongs to that expression. Do not differentiate that restored factor a second time. `unscaled_scalar()` rejects unresolved native factors. The service mode deliberately differs from `original_source`, because the existing scalar transport evaluator cannot consume these factored rows directly.

Five point queries cover y=0,0.53,1 and signed/zero Z values. An independent five-history ODE supplies finite diagnostic-unit references for 40 full p1/p2 value/Z comparisons; maximum discrepancy is about9.38e-13. All76 finite coefficient error records and19 nonzero pressure-tail records are checked. Four exact factor-reconstruction identities and linear pressure sensitivities are checked. The ODE is a diagnostic reference; rigor comes from directed source intervals and the accepted pressure-tail theorem. Its finite R/Pstar/delta units never select native parameters.

## True radius phase with the original origin

`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_O2_radius_phase_points.py/.json` and `_check.py/.json` implement `OriginalO2RadiusPhasePoints(dps=80).evaluate(y=...,N=...)`. The source origin is

```text
Ra=4*epsilon_core
r_minus=Ra*exp(hb*s_c/2)
log(r_minus)=logRa+hb*s_c/2
hb=cstar*K^-100, same actual global physical norm sum
s_c=the explicitly selected original collar endpoint.
```

The exact phase argument on O2 is

```text
N * (log(110/4)+14*logPstar+10*logCstar+1000+y-hb*s_c/2).
```

The selected logCstar is an exact positive dyadic with positive integer exponent. For integer N, 10*N*logCstar and1000*N are integer periods. They are removed exactly by the existing binary modular recipe without allocating their enormous integer. This reduction occurs only inside frac; R, logRref and r_minus retain their original definitions. Common-frame logPstar and logRref recipes are checked explicitly.

The point approximation uses the actual remaining exp/log formula. Independent directed arithmetic encloses it. The original negative N*hb*s_c/2 term remains in the source; its checked log-width bound proves it lies below a positive representable error budget. Subtracting [0,budget] and then periodic projection gives a true phase cover, including wrapping across0. No saved radius/width cap is used as a point value and no microscopic term is declared zero.

The circular point-error contract is arithmetic_error + exp(origin_offset_error_log_upper). At y=0.53,N=7, the phase approximation is0.8834932861130901, arithmetic error upper about1.14e-63, and a proved origin-offset budget1e-100. Four queries use N=1,7,257,2^200+3. Independent 600-digit Decimal exp/ln references agree within the reported arithmetic errors; radius-origin algebra, the positive-offset bound and negative-offset wrap cover are checked.

**Scope:** queries require an explicit positive integer N with at most4096 bits and a successful directed comparison of the origin-offset bound to the chosen numerical budget. Precision adapts to the integer size. This is conditional candidate-N point evaluation, not an all-N phase theorem, globally admitted frequency, full numerical oracle or conditioned loop solution.

## Executable continuation

Treat this file as the current handoff. Completed boxes below are evidence-backed local production layers. Global gates remain false until the actual original functions and errors satisfy their complete domains. Carry forward previous original source hashes and accepted receipts; do not reconstruct large ancestors solely to repeat their checks.

- [x] **O2-FULL-FACTORED-POINT-INPUTS:** E,V,a,b,p1,p2 ordinary value/Z coefficient rows, correlated R/Pstar/delta/L restoration, directed radial and pressure budgets, exact reconstruction and independent ODE/unit references.
- [x] **O2-ORIGINAL-PHASE-ORIGIN:** actual r_minus/source-collar binding and exact integer-period removal.
- [x] **O2-CANDIDATE-PHASE-POINTS:** source-derived exp/log point phase and true directed periodic covers, with the positive microscopic origin term retained. Conditional explicit candidate queries only.
- [ ] **NEXT LOOP-SCALE-SOURCE-BINDING:** read the original `GenericLoopScales` selection and accepted whole-input/boundary receipts. Attach a,b,p1,p2,t0,margin and boundary-excess bounds from the same source family. Preserve the exact selected eta/d_star rules. A local point or a broad cap is insufficient to choose whole-chart/global loop scales.
- [ ] **LOOP-NATIVE-FACTOR-ARITHMETIC:** implement the arithmetic needed for p2/d_star and related terms from `FactoredPointRow`. Retain small complementary/Poisson widths in factored or logarithmic form. Do not evaluate the astronomical R, set delta=0, drop Pstar^-1 sectors or round a near-unit loop parameter to1.
- [ ] **LOOP-CONDITIONED-ANGLE-INVERSION:** adapt the original `GenericLoopPointZ` forward angle, implicit inverse and root brackets to the factored input and phase cover. Use one common source/scales/candidate N; return directed root and arithmetic errors. Preserve periodic branches at0/1.
- [ ] **LOOP-ORDINARY-Z-PRIMITIVES:** consume the completed E_Z,a_Z,b_Z,p2_Z rows exactly once. Compute original A,B_over_Pstar and their slow ordinary Z rows, including q_Z, inverse-angle Z and the defining primitive integrals. Check sign and physical Pstar units against the existing finite scalar backend without replacing native parameters.
- [ ] **LOOP-ERROR-STABILITY:** propagate finite radial/pressure, logarithmic late-source, phase, root and quadrature errors through denominators. Attach strictly positive denominator/root margins before reporting certified A/B values. Report unresolved conditions instead of using an interval endpoint as a point.
- [ ] **ONE-REAL-SIGNED-DENSITY-INTEGRAL:** evaluate one original n-dependent density from the accepted exact function graph using the same point inputs, selected scales, true radius phase and A/B functions. Integrate over its original coordinate/rate, with sign, Jacobian, point/root/quadrature/roundoff and tail errors. A source-density coefficient or saved enclosure is not an evaluated integral.
- [ ] **O2-SOURCE-ROLE-DISPATCH:** bind complete E/V/p1/p2 and A/B value/Z roles to the exact graph hashes. Expose a factored-aware transport wrapper; do not falsely label incomplete factored data `mode='original_source'` for the existing scalar `current_native_Rc_function_transport.evaluate()`.
- [ ] **CORRELATED-ABSOLUTE-PRESSURE:** evaluate cancellation-sensitive P0+cumulative-p sums from their common original integrals. The separate forward error budgets cannot establish exact terminal cancellation at an enormous radius. Keep original raw W and all fourteen pressure stages bound.
- [ ] **OTHER-CHART-POINT-RECIPES:** extend the same source-factor/error interface to remaining original bridge/switch/reshape/reference/O2-axial-buffer/O3 cells. Reuse exact coordinate/radius and derivative conversions. O2 slope is one chart; it does not finish the seventeen-chart source oracle.
- [ ] **ALL-CHART-INTEGRALS/TAILS:** attach original integrands, own-rate integrals, tail/stability/error bounds and same analytic datum over complete original domains. Keep repair-cell coverage and cone assumptions explicit.
- [ ] **COMMON-N/CONTROL-INSTALLATION:** after all required actual function/integral bounds are available, satisfy the original simultaneous finite-frequency constraints with one compatible common N; install all five new controls and functional terminal bounds. Candidate phase queries do not admit N.
- [ ] **CORRECTED-FIELD/JOINS:** recover corrected velocity and analytic pressure from actual controls; complete original core, matching, repair and exact heat-exterior joins with required smoothness and finite energy.
- [ ] **GLOBAL-STRESS/FLAT-REMAINDER:** establish the actual whole-domain stress cone and separate flat remainder with explicit scale-dependent bounds. Do not interpret the present normalized point errors as full momentum residuals.
- [ ] **ACTUAL-COEFFICIENT-RECURSION:** implement distinct n=1 and n>=2 recovery equations on the common core, independent moment repair, divergence-preserving truncation, remainder bounds and smooth summation. Coordinate-rescaled copies of one profile do not count as this recursion.
- [ ] **OSCILLATORY/PHYSICAL-FIELD:** build the original two-family pulses and mean corrections, evaluate averaged quadratic stress cancellation, then expose Cartesian u(x,y,z,t),v(x,y,z,t),w(x,y,z,t) and independent full residual/core-width/aspect-ratio/winding diagnostics.

Long-term goal remains active. No global blow-up, completed scale recursion or full corrected NS field is claimed by these services.
