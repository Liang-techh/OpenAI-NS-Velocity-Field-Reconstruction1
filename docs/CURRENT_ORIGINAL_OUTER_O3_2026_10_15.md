# Original O3 transition/power numerical and physical continuation - 2026-10-15

`OriginalOuterO3Numeric.evaluate(Z,chart,x)` and `OriginalOuterO3Velocity.evaluate(Z,chart,x,log_tau="-2",azimuth=".25")` now continue the actual original selected-N construction through O3_transition x in[0,1] and O3_power x in[0,2], ending at Rc. The same incoming five moment histories, independent P0 and nonzero n1 primitive remain. All six original outer charts now have numerical/physical readers, with the limits below. The full reconstruction remains **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS.

## Original nonzero mu and normalized correction integrals

The original parameter is mu=exp(-4logPstar)/1000, with logPstar=exp40+11. It remains strictly positive in an unbounded-exponent representation. Tw=-60logmu retains the original Rc/repair reservation. The source's accepted mu admission remains separate from any new stress-cone claim.

During transition set J(t)=integral_0^t sigma(u)du, with the original flat sigmoid. During power, t=1+x and J=1/2+x by the exact transition symmetry. Define

```
U_c(t) = expm1(-c*mu*J(t))/mu, c=1 or2
C1(t) = integral_0^t exp(u)*U_1(u) du
C0(t) = integral_0^t U_2(u) du
Cm(t) = integral_0^t exp(-u)*U_2(u) du
```

These order-one normalized quantities preserve the effect of mu before any addition to a large field. They do not compute exp(-mu*J)-1 by subtraction. The reference field Fref below is the exact buffer formula analytically continued to b=11+t; it is an algebraic frame only, not a replacement mu=0 construction. With E0 the actual leading swirl at the buffer endpoint,

```
E = Eref + Eref*mu*U_1
h = href + E0*exp(-3t/2)*mu*C1
e = eref - E0^2*exp(-t)*mu*C0/2
p = pref + E0^2*mu*Cm/2
m = mref, k = kref, V = 0
pressure = P0/Pstar^2+p
```

The separate E/h corrections are negative, e positive and p negative at every nonzero sampled transition point and throughout the sampled power chart. Each is exported in `leading_small_sectors` through TaylorZ2. Larger total enclosures may still be too wide to resolve these tiny differences; exporting a signed sector does not certify relative accuracy of the full pressure or an admissible stress margin. Leading Q=-4*M_axial*exp(-(11+t)) remains Z-independent because the actual mass and V are unchanged by mu.

## Directed transition integration and original power flow

`OriginalO3Scalar` uses the six-state ODE

```
J'=sigma, U_c'=-c*sigma*(1+mu*U_c)
C1'=exp(t)*U_1, C0'=U_2, Cm'=exp(-t)*U_2.
```

A degree4 Taylor endpoint polynomial uses the ODE-derived fifth Taylor coefficient over a full-step value tube as remainder. The source sigmoid supplies jets0..4. Monotonicity and positive-integral comparisons provide additional sign-preserving bounds near the flat endpoint. The defining expm1 identity intersects the U enclosures. At transition endpoint, J=1/2 is imposed as a source symmetry identity. Point proxies remain distinct from directed enclosures.

The power chart is the original `quiet.power_flow`, rewritten in the same normalized variables. The transition C1/C0/Cm inlet is retained. Its local integrals expand

```
expm1(-c*mu*q)/mu = -c*q+c^2*mu*q^2/2-c^3*mu^2*q^3/6+r
0<=r<=c^4*mu^3*q^4/24, q=1/2+u.
```

The guard admits mu<=.001, c<=2 and q<=2.5, so c*mu*q<=.005; the one-sided remainder is valid. Positive weighted polynomial moments use100 exponential Taylor terms and an explicit full-interval remainder. All mu factors remain present. This implements the original finite-mu exprel formulas without subtractive divided differences or a resonance singularity.

## Selected N, endpoints and physical geometry

The transition source cells are [0,1/4], [1/4,1/2], [1/2,3/4], [3/4,1]; power cells are [0,1], [1,2]. Exact rational traversal merges overlapping axial covers and applies dx once. The selected original N is unchanged, incoming corrections already contain1/N, new drivers divide by that N once, and pressure memory is exactly1.

Power has exact same-source flat local A/B and local E/V errors equal to zero. This does not zero the incoming five errors: their homogeneous histories survive. The power origin uses that exact flat source identity rather than copying the transition's wider local error box. Their numerical/physical seams overlap; future correlation-sensitive adapters must retain this same-source endpoint identity rather than treating the two local error boxes as independent arbitrary choices.

Global y=log(R/110) is 10logCstar+10logPstar+exp40+11+x in transition and 10logCstar+10logPstar+exp40+12+x in power. All finite offsets are separate atoms; log110 is included once by physical assembly. Transition ends at Rw and power ends at Rc=Rw*exp2. The accepted all17 pre-Rc n1 source has local axial/swirl zero and constant meridional primitive J throughout both charts, including the closed Rc endpoint. No continuation beyond the Rc repair seam is admitted. The original time map, scaled-time comparison and Cartesian rotation remain.

## Focused evidence

The focused run exits0 in64.820s; saved-report inspection exits0. It exports58 original physical packets, with18 independent constant-sigmoid defining-integral values,9 independent power defining-integral values,16 independent DOP853 values for the original sigmoid,500 numerical/80 physical buffer and transition seam conditions,110 N/n1-history conditions,140 measure/pressure-memory checks,220 P0 conditions,255 nonzero-mu/signed-sector conditions,550 quiet-local-zero/retained-history conditions,4 exact close-boundary classifications,21 scaled-time conditions and7 rejected requests. Maximum tested normalized-profile interval width is1.561214281096204e-9. The independent float comparison is consistency evidence with tolerance5e-12, not the source of directed bounds; maximum point-proxy difference4.500690586188186e-10 is enclosed.

The existing read-only worker `/root/n1_atlas_map` reports GPT-5.6 Luna / max. It reviewed the original support admission, scalar remainder, normalized power formulas and field-sector signs without finding a material algebra or underbound error. Its endpoint-correlation caveat is recorded above. No new agent was spawned.

Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_outer_o3_numeric.json`;23,032,828 bytes; SHA256 `744ef25a25ef8c0e5888a56b0fad5cea10e3db844814fcd216894eb6a1c8d209`;1841 scientific input pins and65 implementation pins. Accepted ancestor scientific input bytes are unchanged.

## Next original construction work and limits

1. Feed the actual Rc leading fields and nonzero selected-N five corrections into the original compact Rc-to-2Rc reserved repair. Use `lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.py` and its original C3 repair source. Recover signed controls/terminal defects in their source-owned normalization; do not use magnitude-cap endpoints as signed values or re-solve a different leading branch.
2. Continue the repaired exit to Rp with the original finite-mu quiet flow and actual complete histories. Do not bypass compact repair by extending the current homogeneous pre-Rc n1 branch. Bind the pulse inlet to the resulting actual source.
3. Supply same-family signed S0/S1 and actual total Q_y/RQ_y for reduced m13/m15 pressure through Rc/repair, retaining physical R/R^2 measures. Current Z jets are not total radial/fast-phase derivatives. Independent P0 plus p does not certify cancellation-sensitive relative pressure accuracy.
4. Finish full n1/global high-order joins, physical C2/common N, compatible analytic preheat datum/exact heat exterior, physical energy/radial tail, admissible stress and flat remainder. No stress-cone margin has yet been established by these O3 field values.
5. Construct independent genuine temporal n>=2 coefficients, smooth sum and two-family oscillatory cancellation, then independently validate fixed-forcing corrected Cartesian NS. Spatial outer-chart completion is not temporal scale recursion or a full corrected NS solution.
