# Original O2 signed densities and first genuine pressure own-rate integral

> Successor: [all five midplane own-rate integrals and incoming history](CURRENT_ORIGINAL_O2_FIVE_OWN_INTEGRALS_2026_10_07.md). The previous NEXT MIDPLANE-ALL-FIVE-OWN-RATES task is completed for exact Z0/C0; nonzero-Z and functional/global transport remain open.

Checked source: [2576ab3f](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/2576ab3ff5daccd00fc395966a6764b1cd8c6837). Predecessor: [original O2 phase-held Z primitives](CURRENT_ORIGINAL_O2_CONDITIONED_SLOW_Z_2026_10_07.md). Focused receipts PASS; direct Git-index audit covers 1070 dependency hashes. Existing read-only reviewer: **GPT-5.6 Luna / max**, no new workers or ancestor reconstructions.

Two production gaps are now closed locally. Actual original O2 source/error and true-radius inverse/primitive rows execute the accepted exact five signed density graph, including Z derivatives. A separate whole-cell source path evaluates one genuine pressure own-rate contribution over the complete O2 y in[0,1] window at exact Z=0 and diagnostic N=7. Its finest enclosure is strictly negative. Point samples are not used as proof of a continuous integral.

This does not close five-moment matching, all-chart own histories, terminal identities as functions of Z, common-N admission, corrected global stress, or coefficient recursion. The pressure contribution must be added to the incoming pressure history with the same P0.

## Executable interfaces

`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_O2_signed_densities.py/.json` and `_check.py/.json` provide:

```python
owner = OriginalO2SignedDensities()
point = owner.evaluate(y='.53', Z='.37', N=7, bits=80)
```

The output contains changed normalized E_N,V_N, their Z rows, both increments and their Z rows, all five signed normalized own-rate densities and their Z rows. Actual source point coefficients/errors, original selected scales, true radius phase and C0/phase-held Z primitives bind exact nodes of the accepted function graph. Density formulas are executed from that graph. Its source log-majorants are not selected as field values.

`BoundDensityGraph(..., jets=None)` explicitly supports value-only graph evaluation for whole-cell sources; its Z output method rejects missing Z bindings. This does not fabricate a zero derivative.

`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_O2_pressure_own_integral.py` and `_check.py/.json` provide:

```python
owner = OriginalO2PressureOwnIntegral()
contribution = owner.integrate(N=7, cells=8192, bits=24)
```

The producer's complete cell evidence is shipped as `lei_ren_part1_paper_compliant_current_original_O2_pressure_own_integral.json.gz`. Decompress with `gzip.decompress()` before JSON parsing. The checker creates this lossless artifact and binds its hash plus the uncompressed JSON hash. The raw 133.7 MB local JSON is ignored; the checked gzip is about 23.8 MB. No cell/source evidence is discarded by compression.

## Signed density graph and tiny changes

With E=Utheta/Pstar, V=Uz/Pstar, x=A/N, e=delta_E=E*expm1(x), and u=delta_V=B_over_Pstar/N, the accepted graph executes:

```text
delta_m = u
delta_h = e
delta_k = V*e + E*u + e*u
delta_e = 2*V*u + u² - E*e - e²/2
delta_p = E*e + e²/2.
```

It preserves nonzero original V, signed cross terms, and the quadratic increments. Z uses the completed ordinary source rows and the original phase-held primitive rows; total phase_Z=0. The E_Z product terms remain. Total-y and higher spatial derivatives are still a separate task: they require the fast A_phi/B_phi terms as well as the slow-y rows.

The normalized own rates are m=1, h=k=3/2, e=1, p=0. Their integration coordinate is y=log(R/Rref). The graph's normalized densities already contain the original radius normalization; applying another R or radial Jacobian is incorrect.

The numerical graph interpreter evaluates expm1(x)=x*integral_0^1 exp(t*x)dt. A 64-term series and a directed exp(1)/66! tail bound enclose the mean exponential. The original x remains factored. Original O2 a=.8+1.2*sigma in[.8,2] and A=a*(phi-psi/(2*pi))/2 give |A/N|<=1 for every candidate N>=1. A proven range intersection is not a selected point field. The positive tiny-q endpoint delta_E and pressure density do not become exact zero.

The density checker compares 180 components across ten independent finite-unit cases. Original scalar loop/slow-Z quadrature gives the reference modulation. Reference densities use original history_densities(changed)-history_densities(original); their Z references use symbolic Jacobians of those original histories. All are enclosed. These diagnostic units do not select native physical parameters.

## Genuine whole-cell pressure contribution

The exact midplane branch is bound to the full original inertial coefficient template and accepted pressure parity P0_Z(0)=0. Symbolic substitution proves original V=p2=b=t0=0, E=f and L=1 there, independently of rounded point outputs. It does not set p2_Z, B or delta_V to zero.

The source f=exp(y/10-3*J_sigma(y)/5), a=.8+1.2*sigma(y), and q²=(2+2*eta-a)/(2*a) retain the same original selected positive eta. The cutoff is exactly one because a<=2. At y=1 the exact source q²=eta/2 is positive. Conservative whole-cell q ranges may include zero but are never treated as an exact flat branch or a selected q point.

Directed monotone endpoint rectangles enclose J_sigma on ordered cells. Forward prefixes and backward suffixes intersect using exact J(1)=1/2. On each whole cell, J lies between its directed left lower and right upper endpoints. f is enclosed by the full y/J rectangle; endpoint f values alone would be unsafe because f is not globally monotone. Production performs no nested defining-mass quadratures at each node. The masses are not needed for this symmetry-bound pressure contribution; the general source/mass path remains open.

The actual radius phase at y=0, including the original positive microscopic origin offset, is propagated as phi(y)=frac(phi0+N*y). Wrapped cell phase ranges are split into closed period pieces using directed bounds. One N=7 is used throughout this integral. No approximate phase midpoint determines a crossing. The original monotone inverse encloses every source root consistent with each joint cell parameter/phase range.

At p2=0 the exact transformed antiderivative identity cancels equal linear angle terms before arithmetic:

```text
A = a*q²*sin(2*psi)/(4*pi*(1+2*q²))
B_over_Pstar = -a*E*q*sin(psi)/(2*pi).
```

This is the original primitive, with an exact symbolic equivalence check. Each phase piece feeds its primitive enclosure to the accepted signed density graph. The pressure density hull covers the entire cell; multiplying it by the exact cell width and summing encloses the genuine signed integral. Pressure rate0 has no Duhamel decay. P0 and incoming cumulative pressure are not reset.

The result is the normalized contribution delta_Mp/Pstar². The physical Pstar² factor remains formal. This number is neither a full momentum residual nor a complete pressure matching error.

The table rounds the endpoints; the receipt retains exact directed MPF bits.

| Directed source cells | Lower contribution | Upper contribution |
|---:|---:|---:|
| 64 | -8.42883e-4 | 8.14796e-4 |
| 256 | -2.33355e-4 | 2.09420e-4 |
| 2048 | -3.97432e-5 | 1.59003e-5 |
| 8192 | -1.88671e-5 | -4.95116e-6 |

The interval width falls from about1.658e-3 to1.392e-5. All refinement intervals overlap; the finest upper endpoint is negative. The four levels cover all seven true phase wrap cells. The fine 8192-cell computation took about146 seconds locally; the complete producer, including source attachment and other levels, took about201 seconds.

An independent defining-J/pressure ODE with scalar phase inversion gives approximately -1.1909127e-5, inside the directed native interval. Two numerical tolerances differ by about6.88e-13. This floating calculation uses the eta=0 limit only as an independent reference at the same fixed native phase, not as a selected native field. An explicit source derivative bound gives |integral_native-integral_eta0|<=2*eta; native eta is enclosed above by1e-520 in the reference bound. That bound applies to the normalized contribution. The certified result remains the actual positive-eta whole-cell enclosure. Additional checks enclose original defining J/f values and the accepted actual midplane density point.

The strict negative sign disproves an assumption that this O2 pressure contribution is exactly zero at N=7. It does not prove N=7 meets the global frequency constraints or determine the sign of the full transported pressure defect after preceding/following charts and repair.

## Next executable tasks

- [x] **O2-SIGNED-DENSITY-DISPATCH:** actual source, true common-candidate phase, C0/slow-Z rows and stable tiny expm1 execute all five signed density values/Z graph nodes. No cap/majorant point field.
- [x] **ONE-ACTUAL-OWN-RATE-INTEGRAL:** full original O2 pressure contribution at exact Z=0, N=7, y in[0,1], rate0, with whole-cell source/phase/inverse enclosures and preserved incoming/P0 memory.
- [x] **INTEGRAL-REFINEMENT/STRICT-SIGN:** 64/256/2048/8192 directed partitions contract and establish a negative normalized contribution; independent defining-ODE limit and source eta-error bound corroborate it.
- [x] **ORDERED-J-PRESSURE-PATH:** source J cells with forward/backward exact-endpoint intersections and full f rectangles, without repeated nested point quadratures. This applies to the exact midplane pressure route.
- [ ] **NEXT MIDPLANE-ALL-FIVE-OWN-RATES:** reuse the accepted compact whole-cell source/error cache as function ranges, not field points. Restore f,a,q, exact source symmetry and true phase unions in the same original basis. Rebuild conditioned primitives and execute all five density graph values; B and delta_V must remain. Enclose exp(-rate*(1-y))*density over each source cell before multiplying by dy; use rates1,3/2,3/2,1,0 exactly. Return five full-window contributions with incoming histories still separate. Independently integrate the five Duhamel ODEs and retain this checked pressure sign. Do not claim Z-functional closure from the midplane.
- [ ] **GENERAL-ORDERED-SOURCE-MASSES:** add genuine ordered J and M0/M1/M2 defining ODE/quadrature values and directed cell errors. Bind complete f,H,D,P to original p1/p2 templates, including both Pstar sectors, positive delta/L and all pressure-tail terms. Preserve the accepted independent pressure alpha. Avoid reconstructing ancestors or selecting coefficient/field cover midpoints.
- [ ] **NONMIDPLANE-WHOLE-CELL-DENSITIES:** extend actual joint source cells to fixed nonzero Z and then Z ranges. Carry native R=Rref*exp(y) in the common formal basis; preserve full p2/p2_Z and source-signed rho/s/inverse-h. Refine or use a valid regular branch across p2=0. All true phase unions must use the same N. Nonmidplane point enclosures alone cannot certify a continuous integral.
- [ ] **Z-DERIVATIVE-INTEGRALS/FUNCTIONAL-TERMINALS:** integrate completed density Z rows with common sources and own rates, including nonzero p2_Z at the midplane. Retain pressure and cumulative histories in a shared analytic representation for cancellation-sensitive exact terminal identities. A sampled axial grid or tiny pressure reference error is not a functional terminal identity.
- [ ] **ACTUAL-RADIAL/HIGHER-SPATIAL-JETS:** implement original slow-y primitive derivatives and fast phase derivatives, then total-y E_N/V_N and a_N/b_N. Supply all ordinary profile jets needed by recovery/stress through the prescribed orders. Fast y terms are order N^0 after the original1/N modulation; they cannot be omitted or inferred from chart selectors. Maintain exact E>0 and original positive denominator certificates.
- [ ] **OWN-HISTORIES/ALL-CHART-TRANSPORT:** carry nonzero incoming five defects through preceding/following charts, quiet gaps, buffer/flatten and exterior joins. Pressure rate0 preserves memory. Combine original inlet histories and the same P0 with changed contributions; do not restart histories at a convenient chart. Build the all-chart numerical source/integral oracle before declaring background matching.
- [ ] **COMMON-N/INDEPENDENT-FIVE-REPAIR:** solve simultaneous original whole-domain frequency, source/stability/tail and repair constraints for one compatible N. Current N=7,19,257 queries are diagnostics. Use correct n-dependent equations and independent five-moment repair; the strict pressure contribution is an input defect, not an automatically canceled quantity.
- [ ] **MATCHED-BACKGROUND/RECURSION/OSCILLATORY/UVW:** install functional five controls, analytic pressure/heat-exterior/finite-energy matching and global admissible stress/flat remainder. Then n=1 and n>=2 recovery, divergence-preserving truncation, flat summation, two-family oscillatory/mean corrections and stress cancellation. Finally deliver physical Cartesian uvw, full independent residual and core-width/aspect-ratio/material-winding diagnostics.

The full long-term objective remains active. Actual coefficient/scale recursion and full corrected NS reconstruction remain incomplete.
