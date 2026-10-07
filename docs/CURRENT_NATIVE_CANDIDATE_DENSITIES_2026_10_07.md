# Native candidate velocities and five signed density kernels

> Successor: [CURRENT_NATIVE_SPATIAL_PHASE_2026_10_07.md](CURRENT_NATIVE_SPATIAL_PHASE_2026_10_07.md) ([4e302b5c](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/4e302b5c46aea4c5ee8e05e2d471d00fe76f9810)) binds original spatial phase to candidate velocities and five signed kernels on five declared source-coordinate requests. Whole coverage and actual cumulative integrals remain open.

Checked implementation: [23cc3d8a](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/23cc3d8aafe6259f3389c8e242be6a562b0ab307). **Actual source A/B now produce finite-N angular/axial velocity increments and all five signed density kernels on the five native boxes admitted by the [conditioned phase stage](CURRENT_NATIVE_CONDITIONED_PHASE_2026_10_07.md).** The original scalar formulas are no longer merely unused fixtures: the new kernel consumes actual native functions. Phase and N=1024 remain explicit candidate parameters. Original spatial phase, whole-chart coverage, slow derivatives and actual cumulative integrals remain open.

## Source, units and computation

E=Utheta/Pstar is the same positive original root. V=Uz/Pstar is recovered from the original packet's axial velocity row in exactly the existing source/radius bases and arithmetic ledger. No additional radius or Pstar multiplier is introduced. Original absolute P0 is retained.

The implementation computes deltaE=E*expm1(A/N), deltaV=(B/Pstar)/N, E_N=E+deltaE and V_N=V+deltaV. It preserves extremely tiny source increments by expm1(z)=z*integral_0^1 exp(t*z)dt. For the enclosing source exponent z in[-1,1], the integral has the directed positive range[exp(min(0,z_lower)),exp(max(0,z_upper))]. The original log-factored z stays outside that range; a tiny z is never substituted by zero. An unresolved or too-large source exponent requires refinement or a larger candidate N.

The five signed dimensionless Duhamel kernels are:

- m=deltaV.
- h=deltaE.
- k=V*deltaE+E*deltaV+deltaE*deltaV.
- e=2*V*deltaV+deltaV^2-E*deltaE-deltaE^2/2.
- p=E*deltaE+deltaE^2/2.

Each square uses the same signed source squared, with a nonnegative enclosure. All original V, energy, pressure and mixed cross terms are present. The corresponding logR recovery rates are1,3/2,3/2,1,0. These are normalized Duhamel kernels, not unweighted physical dR densities. The pressure rate0 means that quiet regions do not erase incoming pressure memory.

## Executed data and API

Same five original boxes: inner_reference/O2_slope at Z=.5 and coordinate.1337; O2_buffer at Z=.5 and coordinate5.337; O3_slope_mu at Z=.5 and offset.537; O3_power at Z=.5 and the original phase cover of.537/Tw. Same seven phases0,.137,.337,.5,.663,.863,1, candidate N=1024.

35 fresh native velocity/density queries produced175 signed kernel enclosures, including12 nontrivial active phase queries. Exact periodic/half-period values and the two flat O3 branches produce zero increments. The O2 buffer retains nonzero tiny angular/axial increments and all five tiny signed kernels at phase.137. These results do not certify every Z/radial coordinate or bind the candidate phase to a spatial phase.

```python
from lei_ren_part1_paper_compliant_current_native_candidate_densities import (
    NativeCandidateDensities, candidate_at_phase)

# phase_owner is the accepted NativeConditionedPhase on the same live seed.
backend = NativeCandidateDensities(phase_owner)
source, loop, V = backend.source('O2_buffer', Z=('.5', '.5'), coordinate='5.337')
candidate = candidate_at_phase(loop, V, phi='.137', N=1024)
deltaE = candidate['values']['deltaE']
five_kernels = candidate['densities']
```

Evidence: same-source native replay of35 queries and175 kernel enclosures;16 independent original finite-N velocity and physical moment-density comparisons with nonzero original axial velocity, opposite signed r, p2=0 and flat inputs; positive and negative unmaterializable expm1 scales retained; tiny five-kernel functions retained; invalid N rejected. Working/index hash closure: 1022. Producer 14.8s and checker 17.2s with the already-live seed. Review worker: read-only **GPT-5.6 Luna / max**.

New scoped gate: `current_native_candidate_velocity_and_five_signed_density_C0_functions_executed`. These are candidate C0 functions. No spatial phase, slow jets, integrated histories, repaired moment functions, common N, modified global cone, coefficient recursion or corrected NS field is admitted. All global completion gates remain false.

## Next production tasks

- [x] **LEFT4c2-density-C0-boxes:** actual native E_N/V_N and all five signed normalized density kernels at explicit candidate phases/N on the five admitted boxes.
- [x] **LEFT4c2-exact-radius-offset (17 source charts; see successor):** evaluate the unchanged original radius tree and r_minus=Ra*exp(hb*s_c/2). Keep hb*s and hb*s_c/2 separate from huge logRa. Bind y=log(R/r_minus) by cancellation of the same source base, not subtraction of independent conservative radius covers.
- [x] **LEFT4c2-spatial-phase on declared requests (see successor; whole coverage remains open):** bind phi=N*y for one explicitly declared candidate integer. Keep integer period separately, split cells that cross an integer, and retain tiny phase offsets near an endpoint. No interval midpoint/floor is a field value. A candidate N is not the theorem's global frequency admission.
- [ ] **LEFT4c2-whole-native-coverage:** extend conditioning, A/B and velocity/density functions across all17 source charts using original Z/coordinate subdivisions and sign unions. A source box crossing the small/large signed-u regimes must stay unresolved until covered.
- [ ] **LEFT4c2-q-slow-jets:** original q_y/q_Z and needed higher orders, with lazy exact flat derivatives, positive active gamma and retained Delta=0 tiny values.
- [ ] **LEFT4c2-phase-and-primitive-jets:** original implicit inverse derivatives at fixed phi, A_y/A_Z and B_y/B_Z, and needed higher orders. Keep hb/Pstar conversion exactly once and use original Phi_psi positivity.
- [ ] **LEFT4c2-density-jets:** derive original source C1 Z and required radial derivatives of all five signed kernels and modified E/V. Preserve every original term and the correlated square/expm1 factors.
- [ ] **LEFT4c2-seams-and-inlet:** same-function phase/density/derivative traces at every original seam and left inlet, with actual chart lengths. Exact zero local increments do not reset incoming histories.
- [ ] **LEFT4c3-integrated-functions:** use actual source phase and signed kernels to compute the five defect histories throughRc. Integrate oscillatory functions with certified quadrature/enclosure/oscillation error; propagate rates1,3/2,3/2,1,0 and retain absolute P0/pressure memory.
- [ ] **LEFT4d-targets:** actual five Rc target functions and C1 Z jets, actual A/A_Z normalization and divided(J-M)/mu. Bounds-only receipts are not target functions.
- [ ] **LEFT4d-controls-and-field:** execute the reserved-band inverse on those targets, bound unique controls and terminal Z identities, then recover the modified radial velocity/pressure/stress and their traces beyond2Rc.
- [ ] **HIGH / CONT / OUTER / ENERGY:** necessary higher derivatives, same-function seam continuation, heat-exterior compatibility and changed-family finite energy.
- [ ] **LEFT4e:** one actual common finite N covering every source/derivative/integration/repair condition, then the full modified cone. N=1024 above is not this admission.
- [ ] **REC / WAVE / PHYS:** actual n-dependent coefficient recovery and independent moment repairs/smooth sum; mean/two-family oscillatory stress cancellation; corrected uvw/p/f and independent Cartesian/dynamical evidence.

Continue with the next bounded production item and commit real function outputs. The [phase handoff](CURRENT_NATIVE_CONDITIONED_PHASE_2026_10_07.md) retains the longer dependency detail. Avoid repeated unchanged proof runs and preserve unrelated edits.
