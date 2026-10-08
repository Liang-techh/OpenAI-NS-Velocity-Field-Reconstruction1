# Original O2 true-radius phase-held Z primitives

Checked source: [db22ee20](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/db22ee20c298a71a29e40e456b18fe0df4210694). Predecessor: [original O2 conditioned C0 primitives](CURRENT_ORIGINAL_O2_CONDITIONED_PRIMITIVES_2026_10_07.md). Focused checker PASS: 48 independent original scalar-Z component comparisons, 16 inverse-coordinate chart comparisons, five native true-radius phase queries and a 1062-hash Git-index dependency audit. Existing read-only reviewer: **GPT-5.6 Luna / max**; no new workers or ancestor reconstructions.

The original O2 phase-held inverse-angle derivative, dimensionless A_Z and B_over_Pstar_Z now have executable source/error enclosures. The genuine p2_Z and E_Z rows are consumed once, with the same selected positive auxiliary scales and original radius/Pstar/delta/L factors. This completes the local derivative layer needed before evaluating signed changed-moment densities. It does not install an actual integral, the all-chart numerical oracle, five controls, a common admitted N, a corrected global field or coefficient recursion.

## Interface and source contract

`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_O2_conditioned_slow_Z.py/.json` and `_check.py/.json` implement:

```python
owner = OriginalO2ConditionedSlowZ()
result = owner.evaluate(y='.53', Z='.37', N=7, bits=80)
```

The output `actual_phase_held_Z_enclosures` binds the C0 inverse bracket and its true radius phase to `psi_Z`, `A_Z_slow`, and `B_Z_slow`. B_Z_slow is the derivative of B_over_Pstar because the source E already equals Utheta/Pstar. Source fields are never chosen from interval midpoints or caps. The outputs enclose functions; they do not select scalar physical field values.

The input y=log(R/Rref) and fast phase are held constant in this derivative. On the original O2 branch a_Z=b_Z=t0=q_Z=0, selected eta and d_star are global Z-independent constants, and the original radius phase has phase_Z=0. The point provider already supplies completed ordinary Z rows, including all L and pressure sectors; differentiating L again would double-count them.

The accepted C0 receipt, full same-family dependency closure, separate numerical precisions and selected original auxiliary log bits remain required. Accepted defining-quadrature approximate coefficients, exact saved bits and their directed finite/late-pressure errors are restored into the live factored point provider. Only identical precision/error contracts may reuse that cache. A new y/Z point still requires genuine source computation. Neither an ancestor build nor original defining quadratures are rerun for the five recorded queries.

## Derivative construction

Write u=p2*q/d_star, h=sqrt(1+u²), r=u/h and s=h^-2. The source u_Z=q*p2_Z/d_star stays factored. At fixed fast phase the implicit derivative is

```text
psi_Z = -T2_Z_fixed_angle / (1+t²)
A_Z = -a*psi_Z/(4*pi)
B_over_Pstar_Z = -a*[E_Z*T1 + E*(T1_Z_fixed_angle+t*psi_Z)]/(4*pi).
```

The strictly positive implicit denominator has exact lower bound 1. Its division retains the original factors. The E_Z product term is present.

At the exact midplane p2=0, the implementation uses the regular source identities

```text
T1 = 2*q*sin(psi)
T1_Z = q*u_Z*sin(2*psi)
T2_Z = 4*q²*u_Z*[sin(psi)+sin(3*psi)/3].
```

Thus p2=0 does not erase p2_Z or A_Z. Intervals containing zero use the regular Fourier branch, never a guessed sign or a 1/r expression. The 48-mode antiderivative derivative series has explicit geometric tail bounds for W1, W1_r and S2. Exact periodic and half-period O2 symmetry points return exact zero derivatives. At y=1 the exact a=2 identity and positive tiny q remain; the endpoint derivative is not replaced by zero.

For signed large u, the derivative of the transformed original antiderivatives uses both the original-angle and Mobius-angle coordinates. Positive rho=1/[h(h+abs(u))], s=h^-2 and inverse-h remain source factors. The original-angle branch forms stable source-signed denominators from rho² and half-angle squares. The Mobius branch retains fixed-original-angle chi_Z=2*(u_Z/h)*sin(chi), even though chi is the inverse's numerical coordinate. Direction and derivative formulas use correlated half-angle identities. No rounded subtraction 1-r² defines s, and no huge u or h is materialized.

Exact symbolic identities compare the differentiated transformed T1/T2 to the original formulas and confirm the Fourier derivative coefficients. The read-only Luna reviewer independently checked the regular series, its tails and the nonzero midplane limit.

## Evidence and limits

The five native queries are y=.53,Z=-.37,0,.37 at N=7; y=1,Z=.37 at N=257; and y=.53,Z=.37 at diagnostic N=19. The last genuine radius phase lies in its source-signed narrow peak and executes the Mobius-coordinate derivative. N=19 was chosen solely to exercise that branch; it is not a common frequency admitted by the global constraints.

Independent references restore the accepted complete source coefficients in diagnostic finite units (R,Pstar,delta)=(50,7,.1) and (70,13,.01). Original `GenericLoopPointZ` angle quadrature and implicit differentiation supply the reference values. All 48 comparisons of psi_Z, A_Z and B_Z lie in the new enclosures, across 16 coordinate-chart cases. Both signed charts, the exact nonzero-derivative midplane, the regular small-r endpoint and exact symmetry points are exercised. These finite units do not select native parameters.

The same separate precisions remain: 260 interval digits for auxiliary log bits and factored arithmetic, 50 defining coefficient digits, 90 directed coefficient digits, 80 radius phase digits with N guard digits, and 4096 radial error cells. This local check does not assert a complete derivative error bound over the entire O2 domain, momentum residual, five-moment terminal closure or full NS reconstruction.

## Next executable tasks

- [x] **O2-PHASE-HELD-Z-PRIMITIVES:** original implicit inverse and A_Z/B_over_Pstar_Z on actual source/error and true radius phase; complete p2_Z and E_Z factors retained.
- [x] **Z-CONDITIONED-DERIVATIVE-ARITHMETIC:** source-factored rho/s/inverse-h, both inverse coordinates, positive implicit division and original derivative identities.
- [x] **MIDPLANE/ENDPOINT-Z-BRANCHES:** regular p2=0 with p2_Z nonzero; positive tiny-q a(1)=2 endpoint; no sign choice on crossing intervals.
- [x] **ACCEPTED-POINT-CACHE-REUSE:** reuse only actual defining-quadrature coefficient bits and directed errors under the same accepted source/hash/precision contract. Five recorded queries do not rebuild ancestors or nested source quadratures.
- [ ] **NEXT O2-SIGNED-DENSITY-DISPATCH:** read the accepted exact n-dependent velocity/five-increment graph. Connect native E,V,A,B and completed phase-held Z rows at the same y,Z,N. Apply the actual radial phase derivative where required; slow-Z alone is not the full spatial jet. Preserve all five signed cross terms, full p1/p2 rows and denominator roots. Return five density enclosures with a common source family/basis/error ledger, and explicitly keep `actual_changed_five_moment_integral_evaluated=False` until integration exists. Compare finite diagnostic units against the existing exact graph, not a formula copied into the checker.
- [ ] **FAST-ACTUAL-RADIAL-POINT-PATH:** map the defining J_sigma and M0/M1/M2 equations and exact initial/endpoint data. Implement ordered-node ODE or equivalent defining quadrature reuse, with source values and directed errors. Avoid repeating nested quadratures at every y node. Approximate values must come from the defining equation and an error enclosure; do not select cover midpoints as physical values. Retain the independent original pressure alpha and cached accepted pressure-tail proof. Benchmark ordered nodes against accepted original coefficients and report runtime/error separately.
- [ ] **ONE-ACTUAL-OWN-RATE-INTEGRAL:** choose one genuine signed increment density and its original source window/rate. Substitute R=Rref*exp(y) in the common formal basis before applying the radius rate/Jacobian, so exact factors cancel before exponentiation. Use the same explicit N in every density node and actual frac(N*log(R/r_minus)) phase; handle period boundaries and narrow inverse peaks. Propagate coefficient, phase/origin, inverse, Fourier, integration and roundoff errors. Produce a signed integral enclosure tied to its own rate; one partial increment does not install all five controls or terminal identities.
- [ ] **DENSITY-AND-INTEGRAL-ERROR-REFINEMENT:** vary only numerical precision/cells/partition for the same source and candidate N; show useful contraction of the integral enclosure. Keep cancellations within their shared source representation. If a forward enclosure is too wide, improve correlation or arithmetic before increasing native physical parameters or claiming cancellation.
- [ ] **CORRELATED-PRESSURE/TERMINAL-Z:** preserve common pressure datum and cumulative radial histories for exact cancellation-sensitive terminal identities as functions of Z. Independent forward pressure error boxes at astronomical radius cannot certify exact terminal matching. Supply the original analytic identity/proof and a matching numerical representation.
- [ ] **REMAINING-CHARTS/WHOLE-ORACLE:** extend actual defining point/error recipes, completed spatial jets and densities to the remaining bridge, switch, reshape, reference, axial/buffer and O3 charts. Preserve chart joins, complete domains, tail and repair conditions. O2 sampled points do not close the all-chart original source oracle.
- [ ] **COMMON-N/CONTROL/CORRECTED-FIELD:** solve simultaneous original whole-domain frequency constraints with one N, install five controls and functional terminal bounds, recover corrected velocity/pressure and original heat-exterior/finite-energy matching. Admit global stress-cone and flat-remainder conditions separately.
- [ ] **ACTUAL-RECURSION/OSCILLATORY/UVW:** after background matching, implement n=1 and n>=2 recovery, independent moment repair, divergence-preserving truncation and flat summation. Then implement two-family oscillatory pulses/mean correction and averaged quadratic stress cancellation. Deliver physical Cartesian uvw, full independent residual and core-width/aspect-ratio/material-winding diagnostics.

The long-term goal remains active. Genuine coefficient/scale recursion and full corrected NS reconstruction are incomplete.
