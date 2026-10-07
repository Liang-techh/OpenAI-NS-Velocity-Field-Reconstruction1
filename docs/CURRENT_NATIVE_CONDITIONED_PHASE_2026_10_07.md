# Native conditioned phase inverse and original A/B functions

> Successor: [CURRENT_NATIVE_CANDIDATE_DENSITIES_2026_10_07.md](CURRENT_NATIVE_CANDIDATE_DENSITIES_2026_10_07.md) ([23cc3d8a](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/23cc3d8aafe6259f3389c8e242be6a562b0ab307)) now executes native candidate velocities and all five signed C0 density kernels on these boxes. Spatial phase and actual integrals remain open.

Checked implementation: [50d5ccee](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/50d5ccee39750398e221347df47e227c4794c404). **The original Section11 phase inverse and A/B primitives now execute on five native source boxes, including two whose Poisson peaks cannot be resolved by materializing their width.** These are C0 function enclosures at explicit candidate periodic phases. Original spatial phase binding, slow derivatives, whole radial coverage, moment integration and actual coefficient recursion remain open.

## What was produced

The previous [correlated q backend](CURRENT_NATIVE_CORRELATED_SHEAR_Q_2026_10_07.md) supplies the same original a,b,p2,t0,E and q. The new backend computes u=p2*q/d_star with the same whole-family d_star, keeps rho=1/[h(h+abs(u))] and s=h^-2 positive and log-factored, and avoids subtracting rounded r^2 from1. It never exponentiates an enormous positive log or chooses a source midpoint.

For large fixed-sign u, the Mobius angle E_angle has positive derivative s/D. The original W2 expression is multiplied by its compensating s before evaluation:

T2=t0^2*psi+2*t0*q/(h*r)*(E_angle-psi)+q^2/r^2*((2-3*s)*E_angle+s*psi+2*r*sin(E_angle)).

The exact source phase is Phi=(psi+T2)/(2*pi*(1+t0^2+2*q^2)). The code uses both psi/(2*pi) and E_angle/(2*pi), reflects angles across the half period to keep the correct signed branch, and uses exact endpoint/half-period identities. Small signed r, including zero, uses the original Fourier antiderivatives with explicit geometric tail errors; the formulas do not divide by r. A source box whose u crosses the small/large sign regimes remains unresolved.

Each inverse has two independently directed monotone searches: one encloses the lower root and one encloses the upper root using actual source phase interval endpoints. The more informative angle chart is selected. A returned bracket encloses every source root in the box; its midpoint is never a field value.

Original A=a/2*(Phi-psi/(2*pi)) and B/Pstar=E/2*(-a*T1/(2*pi)-b*Phi) are evaluated by algebraically equivalent factored formulas. In particular B/Pstar=E*(t0*A-a*q*W1/(2*pi*h)). Keeping q and q^2 outside finite cancellation preserves the tiny nonzero buffer primitives. Periodic and half-period primitive values, and all flat-branch primitive values, are exactly zero.

## API and executed scope

```python
from lei_ren_part1_paper_compliant_current_native_conditioned_phase import NativeConditionedPhase

# q_owner is the accepted NativeCorrelatedShearQ using the same native seed.
backend = NativeConditionedPhase(q_owner)
source, loop = backend.query('O2_slope', Z=('.5', '.5'), coordinate='.1337')
result = loop.evaluate('.337')
values = loop.primitives(result['selected_inverse']['coordinate_interval'],
                         result['selected_inverse']['chart'])
```

| Original source box | C0 conditioning |
|---|---|
| inner_reference, Z=.5, coordinate .1337 | Negative u with enormous positive log; two-angle inverse |
| O2_slope, Z=.5, coordinate .1337 | Negative u with enormous positive log; two-angle inverse |
| O2_buffer, Z=.5, coordinate5.337 | Tiny q retained; small-r series with signed p2 cover |
| O3_slope_mu, Z=.5, offset .537 | Exact flat q, identity phase and zero primitives |
| O3_power, Z=.5, original phase cover of .537/Tw | Exact flat q, identity phase and zero primitives |

All five use phases0,.137,.337,.5,.663,.863,1:35 fresh queries, including12 nontrivial active inverse/primitive queries. The largest active phase-image width is below9.1e-13. These point-coordinate boxes do not certify the full Z/radial charts. The previous17 full-Z source q query boxes remain separately valid.

Evidence:16 independent original scalar phase/A/B comparisons across opposite signs, p2=0 and flat cases;8 unmaterializable signed-u inverse cases; positive factored rho/s; nonzero tiny-q primitives; rejection of an arbitrary sign for a wide signed-u box; five fresh native source and result replays. Working/index hash closure: 1018. Producer 14.5s; focused checker 18.4s with the already-live seed. Review worker routing: **GPT-5.6 Luna / max**, read-only.

New scoped gate: `current_native_conditioned_C0_phase_inverse_and_primitives_executed`. No spatial phase/N binding, slow derivative, changed density integration, moment repair, global finite N/cone, coefficient recursion or corrected NS admission follows. All global completion gates remain false.

## Detailed next production handoff

- [x] **LEFT4c2-conditioning-C0:** positive factored rho/s and small/large signed source conditioning.
- [x] **LEFT4c2-inverse-C0-boxes:** directed two-coordinate Phi inverse and original A/B C0 on the five boxes above.
- [x] **LEFT4c2-density-C0 on these five boxes (see successor; whole coverage remains open):** use these actual A/B enclosures to compute E_N=E*exp(A/N), V_N=V+(B/Pstar)/N for an explicitly declared positive integer candidate N. Preserve tiny expm1(A/N) in log coordinates. Compute all five signed increments using original V and E. Candidate phase/N must remain explicitly separate from global admission.
- [ ] **LEFT4c2-native-coverage:** subdivide original Z/coordinate boxes that cross u=0 or the small/large regimes; retain a complete union of sign branches and directed tails. Extend phase/A/B functions to all17 modification charts and inlet traces, not only the five samples.
- [ ] **LEFT4c2-q-jets:** execute original q_y/q_Z and necessary higher orders. Preserve the exact flat derivatives and conditional positive active gamma. The current backend admits q C0 only.
- [ ] **LEFT4c2-phase-slow-jets:** differentiate the original phase implicitly at fixed phi, using the exact positive Phi_psi. Keep the transformed coordinate correlations and original width/Pstar conversions; do not differentiate a selector or an inverse midpoint.
- [ ] **LEFT4c2-A-B-slow-jets:** original A_y/A_Z and B_y/B_Z, then needed orders for (11.18), stress and physical recovery. Reuse original signed inertial roots and absolute P0.
- [ ] **LEFT4c2-exact-spatial-phase:** evaluate original y=log(R/r_minus) while retaining hb*s and hb*s_c/2 separately from huge base logs. Bind one common integer N and its fractional phase consistently across seams. Do not floor an interval crossing an integer; split periodic cells.
- [ ] **LEFT4c2-source-seams:** execute same-function source traces, phase traces and derivatives at every bridge/micro/power/restore/patch/O2/O3 seam. Use original chart coordinates and actual lengths.
- [ ] **LEFT4c3-actual-integrals:** integrate all five signed density functions from the real left inlet toRc, with rates1,3/2,3/2,1,0. Preserve full meridional, pressure, energy and mixed terms, actual inlet histories and nondecaying pressure memory. Certify oscillation/interval errors.
- [ ] **LEFT4d-target-functions:** evaluate actual five Rc defects and C1 Z derivatives in the same units, including A_Z/A and divided (J-M)/mu. Conditional defect bounds alone are not target functions.
- [ ] **LEFT4d-controls:** solve actual five target equations with the reserved-band operator, bound unique controls and preserve all five terminal Z-function identities. Keep inverse-mu loss explicit.
- [ ] **LEFT4d-recovered-field:** recover modified radial velocity, pressure and stress from the actual modified histories and original absolute P0 on(Rc,2Rc); prove traces agree beyond2Rc.
- [ ] **HIGH / CONT / OUTER / ENERGY:** required derivative orders, seam continuation, exact heat exterior and finite energy for the changed family.
- [ ] **LEFT4e:** admit a single actual finite N satisfying every source/derivative/repair condition, then certify the full modified cone. Historical scoped N is not this N.
- [ ] **REC:** recover genuine n=1/n>=2 coefficients with each order's equations and independent moment repair, then finite remainder and smooth sum. Rescaling an existing field does not accomplish coefficient recursion.
- [ ] **WAVE / PHYS:** mean corrections and two-family oscillatory stress cancellation, flat remainder/forcing, complete physical uvw/p/f, independent Cartesian residuals and contraction/elongation/material winding diagnostics.

Pick the next bounded production item, execute it on the same native seed, commit actual code/results, and mark only its achieved scope complete. Avoid repeating unchanged proof stages. Preserve unrelated worktree edits.
