# Original O2-slope numerical and physical continuation - 2026-10-15

`OriginalOuterSlopeNumeric.evaluate(Z,x)` and `OriginalOuterSlopeVelocity.evaluate(Z,x,log_tau="-2",azimuth=".25")` now continue the same actual original selected-N field across the complete O2_slope chart, x in[0,1], R=Rref*exp(x). The full reconstruction remains **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS. Original parameters, source family, N, independent P0 and nonzero n1 primitive are retained.

## The original profile and all its leading memories

Set A=1/(1+Z^2), B=4Z/Pstar. The original sigma and exact original slope law reduce to four scalar functions:

```
F' = (1/10-3*sigma(x)/5)*F, F(0)=1
H' = F-3H/2,              H(0)=5/8
G' = F^2-G,               G(0)=5/6
C' = F^2,                 C(0)=5
E=A*F, V=m=B, h=A*H, k=A*B*H
e=B^2-A^2*G/2, p=A^2*C/2
pressure=P0/Pstar^2+p
```

These are the original own-rate equations, including the nonlocal h/k/e/p memories from Rref. They are not reinitialized to an equilibrium at each radius. Original symmetry sigma(1-x)=1-sigma(x) gives integral_0^1 sigma=1/2, so F(1)=exp(-1/5) exactly. Endpoint scalar proxies are F~0.8187307531, H~0.6321060845, G~0.8926802504, C~5.9589254140; the exported intervals, rather than extra printed digits, determine their accuracy.

## Directed integration rather than an unbounded ODE approximation

`OriginalSlopeScalar` uses512 fixed native cells, a degree4 Taylor endpoint polynomial and an interval fifth-coefficient remainder over the entire step. The positive state tube follows differential comparison with -1/2<=1/10-3sigma/5<=1/10. H/G/C bounds retain positivity and conservatively bound their forcing. Ordinary ODE recurrences derive the state fifth coefficients from the actual state tube and sigma Taylor0..4; no value-only tube is substituted for derivative bounds.

All recurrence coefficients include1/k!. The endpoint remainder is h^5 times the fifth Taylor coefficient, with no second factorial. Exact rational coordinates determine grid and partial steps. The exact endpoint F identity intersects the computed directed enclosure. Point states are midpoint-coefficient Taylor proxies; full error relative to each serialized point is exported separately.

The largest measured scalar directed interval width is7.02258532260907e-10. Eight checks enclose independent exact solutions for constant sigma0/1, including the resonant G equation. Twenty scalar values agree with an independent DOP853 solve (SciPy1.17.1), allowing5e-12 for that floating solver. Its maximum difference from the Taylor point proxy is2.834712464760969e-10. This floating comparison is a consistency check; the directed Taylor construction supplies the mathematical error bounds.

## Actual selected-N corrections and physical values

The numerical adapter transports the five incoming selected-N correction intervals from the accepted Rref packet, using dy=dx once, rates1,3/2,3/2,1,0, and the original all-N density caps. Incoming errors already include1/N; new drivers divide by the same original DoubleDyadic N once. Overlapping Z covers cannot double-count a native cell. Pressure memory remains exactly1, and P0 is added once.

E_N/V_N remain local algebraic source functions, separate from cumulative history transport. Their actual Rref local enclosures are reused at x=0; later points use only the same-family local source bounds at that radius. No candidate N257/N1024 data or magnitude cap endpoint is chosen as a source value. Directed leading integration error and unknown-signed finite-N errors both remain in field enclosures. Fast phase/signed finite-N centers are not resolved.

The shared physical map exports cylindrical and Cartesian leading+n1 sums with atomic y=log(R/110)=10logCstar+10logPstar+x. The factor110 and original time powers are included once. The unchanged n1 primitive still yields nonzero radial velocity; local n1 axial/swirl values remain zero before reserved repair. Its domain admission uses log(Rw/(e*Rref))=logPstar>0. Native coordinates and factored sums remain required; no huge physical exponential or direct Cartesian inverse is introduced.

## Focused evidence and limits

The focused run exits0 in50.790s and saved-report inspection exits0. It exports28 actual physical packets over5 Z positions and ordinary/scaled time, with250 numerical and40 physical Rref seam checks,55 N/n1/local-inlet checks,75 physical-measure/pressure-memory checks,100 P0 checks,21 scaled-time checks,2 close-boundary classifications and5 rejected requests. Existing read-only GPT-5.6 Luna/max independently confirmed the scalar source reduction, coefficient recurrence, value tube and fifth-coefficient remainder.

Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_outer_slope_numeric.json`;10,447,444 bytes; SHA256 `875e1c5723f301ef6da241168df30020596c719d67c4693143989581ff373273`;1838 scientific input pins and57 implementation pins. Accepted ancestor source and report bytes remain unchanged.

Public values include TaylorZ2 generic/axial/swirl/mass enclosures and radial Q through Z1. Total radial/fast-phase derivatives, complete physical C2 joins and full n1 pressure remain open. The separate signed reduced-pressure S0/S1 and Q_y/RQ_y interface for m13/m15 has not been supplied; generic P0+p does not replace it. No global energy, exact heat, stress/flat remainder or temporal n>=2 claim is made.

## Next original construction work

1. Feed the actual O2_slope endpoint F/H/G/C, finite-N five histories, P0 and J into O2_axial. Its native coordinate s in[0,1] uses physical log distance d=exp(40s)-1 and dy=40exp(40s)ds. The original beta=sigma(1-s) gradually removes axial velocity; do not integrate it with unit native measure.
2. Preserve the memory-sensitive axial moment convolutions and the resonant E^2 history. Local swirl has E=E1*exp(-d/2); h and pressure moments have exact homogeneous/forced recovery. Large d needs stable factored exponentials and controlled convolution tails, not zero resets.
3. Continue O2_buffer, O3_transition and original O3_power through Rc, then reserved moment repair. Retain the original mu and source power flow, all selected-N memories and the actual signed terminal defect enclosure.
4. Add same-family total Q_y and signed S0/S1 numerical pressure history at Rh/Rref and later interfaces, with correct R/R^2 weights. Complete full n1 pressure, high radial/Z bounds and global joins/common N before admitting physical C2 or five-moment closure.
5. Complete exact heat/analytic preheat compatibility, energy/radial tail, admissible stress/flat remainder, then genuinely independent temporal n>=2 coefficients, smooth summation, two-family oscillatory cancellation and fixed-forcing independent Cartesian NS validation. Spatial outer-chart continuation does not constitute temporal scale recursion.
