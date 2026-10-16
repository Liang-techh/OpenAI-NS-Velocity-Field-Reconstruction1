Latest continuation: [Original O2-slope numerical/physical values](CURRENT_ORIGINAL_OUTER_SLOPE_2026_10_15.md) now carries this source through the next complete slope-transition chart with directed integration error. The accepted Rh/Rref scientific source and its limits below remain unchanged.

# Original selected-N Rh-to-Rref numerical and physical continuation - 2026-10-15

`OriginalOuterReferenceNumeric.evaluate(Z,offset)` and `OriginalOuterReferenceVelocity.evaluate(Z,offset,log_tau="-2",azimuth=".25")` now continue the same actual original selected-N construction across the complete reference annulus, offset in[-5,0]. Offset=-5 is Rh and offset=0 is Rref. The full reconstruction remains **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS.

This is the first of the six prescribed outer charts. It extends computational coverage beyond the active five-bump patch without substituting old candidate N257/N1024 data. It does not close the remaining outer charts, signed pressure matching, stress or temporal recursion.

## Same-source analytic background and geometry

The accepted original smooth Rh function identity gives, in generic S=Pstar units,

```
A = 1/(1+Z^2), V = 4Z/Pstar
E = A*exp(offset/10)
m = V, h = 5E/8, k = 5EV/8
e = V^2-5E^2/12, p = 5E^2/2
pressure = P0/Pstar^2+p
R = 110*exp(10logCstar+10logPstar+offset)
dy = d(offset)
```

The leading values and their TaylorZ2 rows follow this smooth-function identity, rather than an interval-overlap assertion. Original C6/C5 outer source and Rh-identity receipts are pinned. The independently normalized analytic P0 is retained and added once. Physical axial/swirl/mass values have Z2 enclosures and radial Q has Z1, recovered from the same modified V/m and full m_Z. No higher derivative row is padded.

The physical wrapper retains the same lambda/time map and Cartesian rotation as the accepted Rm/Rh wrapper. Radius geometry remains atomic with log110 included once. Ordinary log_tau and source-scaled S=-delta*log(tau)/2 both work; huge absolute physical exponentials are not materialized. Values remain sums of bounded coefficients times separate exponential factors, requiring native (Z,offset,azimuth,time) coordinates.

## Finite-N fields and five histories have different roles

Only the cumulative m/h/k/e/p errors follow Duhamel transport. They enter from the actual accepted Rh packet, already scaled by1/N. Their rates remain1,3/2,3/2,1,0. Each new signed driver is bounded by the original all-N magnitude archive divided by the same original DoubleDyadic N once. Pressure memory is exactly1. Five exact unit cells cover[-5,0]; overlapping Z covers cannot duplicate the physical measure or predecessor memory. Exact rational comparisons classify partial and near-boundary queries before rounding.

The local E_N=E*exp(A_loop/N), V_N=V+B_loop/N are algebraic functions at the current radius, not cumulative histories. At exact Rh the accepted local E/V enclosures are reused, preserving the shared endpoint. At later points the same-family local source bounds are evaluated there; no incoming E/V error is additively transported or double-counted. Magnitude cap endpoints are never chosen as signed source values. Leading point proxies retain full unknown-signed finite-N error enclosures; the actual fast phase and signed integral centers remain unresolved.

The actual n1 terminal primitive J and its Z derivatives continue unchanged. Local u1/f1=0 and the original pre-repair own-field identity imply J_y=0. The support check extends to Rref using log(Rw/Rref)=logPstar+1>0 and Rc=Rw*exp(2). Thus the local branch stays before reserved repair. Nonzero radial n1 velocity remains in the physical sum. This does not admit the still-open global high-order n1 field joins.

## Evidence and remaining pressure interface

The focused run exits0 in49.238s and saved-report inspection exits0. It exports23 actual physical packets over5 Z positions and ordinary/scaled time. Checks cover250 numerical Rh seam coefficients,40 signed physical sum seams,40 original-N/n1 history conditions,65 physical-measure/pressure-memory conditions,80 P0 checks,30 independently integrated own-density ODE values,21 scaled-time conditions,2 exact close-boundary queries and5 invalid-request rejections.

The independent ODE comparison uses the original higher-precision amplitude context to normalize exp(-logPstar) values. It propagates the accepted serialized inlet's directed uncertainty with its homogeneous kernel; a saved point proxy is not treated as exact. This avoids creating an artificial discrepancy by rounding the huge logPstar before exponentiation. It is a value-level consistency check, not certification of missing global norms.

Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_outer_reference_numeric.json`; 8,741,719 bytes; SHA256 `3e97916aab982bcf5b21690a3843787e630936c8a024806ccc643c0073c14db0`;1838 scientific input pins and54 implementation pins. Existing read-only GPT-5.6 Luna/max supplied the source/measure map and confirmed the local E/V versus cumulative-history distinction. Accepted ancestor source/report files are unchanged.

The generic temporal-order-zero pressure P0+p is available. The separate signed reduced-pressure integrals S0/S1 and boundary Q_y/RQ_y needed for m13/m15 have NOT been continued numerically. Neither these integrals nor full n1 pressure can be inferred from current Z rows. Exact heat/preheat pressure compatibility remains open.

## Next tasks for the original route

1. Continue Rref into O2_slope, native[0,1], physical offset6+x. Use E=A*exp(x/10-3*integral_0^x sigma/5), V=4Z/Pstar; carry all five actual incoming histories instead of replacing them with local equilibrium values.
2. Continue O2_axial using offset6+exp(40x), dy=40exp(40x)dx, then O2_buffer[0,11], O3_transition[0,1], O3_power[0,2] to Rc. Consume the original smooth sigma, original mu and exact power-flow source. Preserve every original incoming memory and the full physical measure.
3. Export same-family Rh signed S0/S1 and actual total Q_y/RQ_y numerical enclosures. Continue m13=-S0/2+(Q+Q_y-initial[Q+Q_y]) and m15=S1/2-(RQ_y-initial[RQ_y]) with their true physical radius weights. Keep P0 and full n1 pressure responsibilities separate.
4. Carry actual modified terminal defect intervals to reserved repair; complete high radial/fast-phase and higher-Z bounds and global field/pressure joins before claiming physical C2, common global N or five-moment closure.
5. Complete exact heat exterior, energy/radial tail, admissible stress and flat remainder, then independent temporal n>=2 coefficients, smooth summation, two-family oscillatory cancellation and independent fixed-forcing Cartesian NS validation. Current outer coverage is not temporal scale recursion.
