# Original radius offsets and spatial candidate functions

Checked implementation: [4e302b5c](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/4e302b5c46aea4c5ee8e05e2d471d00fe76f9810). **The actual phase N*log(R/r_minus), with candidate N=1024, now feeds native modified angular/axial velocities and all five signed density kernels at five declared source-coordinate requests.** The previous free-phase stage has a spatial successor. The original17 radius charts have offset/Jacobian representations and bounded point phase enclosures. This does not admit whole-chart modified velocities, actual cumulative histories, repair, or a global common N.

## Executed source arithmetic

The same original left radius is r_minus=Ra*exp(hb_bridge*s_c/2). Ra=4*exp(-4*logPstar-1000), logRref=log110+10*(logCstar+logPstar). Matching source terms cancel before numerical arithmetic. Microscopic hb_bridge/hb_switch offsets remain positive log-factored functions and never enter a rounded huge absolute logR. The existing global radius operator supplies its original descriptor; its cap-based numeric radius cover is not consumed as a field value.

Original logCstar and T=400*Abar are source-selected singleton binary constants. The original selection assignments are AST-bound. Explicit exact rational coordinate coefficients permit exact modular arithmetic on their binary tuples, without allocating their huge integer cycles. Analytic logPstar/Tw remain intervals. An interval coordinate is never silently replaced by an exact literal or midpoint. For example, exact literal .1337 in the inner-reference chart has a bounded phase, whereas its finite-precision interval cover conservatively gives a full period because of the enormous selected T/C multipliers. These are different requested sets.

All17 affine offset identities,17 original native-coordinate Jacobians and16 adjacent radius/periodic-phase seam identities are checked symbolically. O2 axial uses actual_y=exp(Md*phase) and Jacobian Md*actual_y; buffer uses actual_y=exp(Md)+offset and Jacobian1. The O2/O3 seam uses original logPstar=exp(Md)+11. These seam checks do not certify modified velocity/stress derivative matching.

## Data and callable interface

17 original point-radius queries at Z=.5 retain bounded periodic cells. At the actual selected inlet s=s_c/2, the same radius cancellation gives exactly zero offset and phase. Positive microscopic bridge offsets at the two sampled bridge points remain nonzero formal sources.

At inner_reference/O2_slope coordinate.1337, O2_buffer offset5.337, O3_slope_mu offset.537, and O3_power original offset.537/Tw, the derived spatial phase feeds5 candidate velocity/density cells and25 signed kernel enclosures. Three cells have nonzero kernels; the two q-flat O3 cells have exact zero increments. The O3 explicit offset is divided by the same analytic Tw, so Tw*phase=.537 cancels as a source expression.

```python
from lei_ren_part1_paper_compliant_current_native_spatial_phase import (
    NativeSpatialPhase, spatial_candidate_functions)

# native_owner and density_owner share the already-live original seed.
binder = NativeSpatialPhase(native_owner)
result = spatial_candidate_functions(
    density_owner, binder, 'O2_buffer', ('.5', '.5'), '5.337', 1024)
geometry = result['actual_source_radius_and_phase']
cells = result['actual_spatial_candidate_cells']
```

Evidence: original same-source17-query replay, actual inlet,5 spatial candidate replays,25 signed kernels,17 Jacobian and16 seam identities, independent CRT for a selected binary exponent10^40, signed/negative-exponent modulus, wrap/full-period union and exact-vs-interval coordinate checks. 1023 working/index hashes pass. Producer 21.25s; focused checker 21.88s on the already-live source. Read-only routing: **GPT-5.6 Luna / max**.

## Next production tasks

- [x] **LEFT4c2-exact-radius-offset:** same-source affine log-radius offsets and coordinate Jacobians on17 original charts; microscopic positive widths retained.
- [x] **LEFT4c2-spatial-phase on declared requests:** actual candidate N*y periodic enclosures, endpoint unions, selected binary/source-expression cancellations and exact inlet phase.
- [x] **LEFT4c2-spatial-density on five requests:** same-source spatial phase feeds original modified E/V and all five signed kernels. Entire-chart coverage remains open.
- [ ] **LEFT4c3-local-signed-integral:** take an actual narrow active cell in O2_slope around.1337, query its entire coordinate interval and true spatial phase, and enclose each signed Duhamel contribution with the exact positive kernel mass. Report contribution functions, not a globally initialized history. Retain rates1,3/2,3/2,1,0 and absolute P0.
- [ ] **LEFT4c3-source-inlet-and-propagation:** attach original five incoming histories from the same source r_minus and original P0; transport the exact zero defect inlet through quiet collars without resetting pressure memory. Continue cell contributions across original chart lengths; split phase/sign/conditioning cells as needed.
- [ ] **LEFT4c2-whole-coverage:** cover every original Z/radial chart interval for q, conditioned inverse, A/B and modified velocities/densities. Exact point requests are not interval-wide coverage.
- [ ] **LEFT4c2-q-and-phase-slow-jets:** actual q_y/q_Z, implicit inverse/A/B derivatives, higher needed orders and same-function seam traces. Apply original radial/width/Pstar conversion exactly once.
- [ ] **LEFT4c3-C1-histories:** actual five cumulative defect functions and C1 Z derivatives through Rc, with certified signed integration error. Local C0 contributions cannot substitute for these functions.
- [ ] **LEFT4d-targets-and-controls:** compute actual Rc target functions/A/A_Z and divided(J-M)/mu, execute the reserved-band inverse, then prove five terminal Z-function identities and changed radial velocity/pressure/stress continuation beyond2Rc.
- [ ] **HIGH / CONT / OUTER / ENERGY / LEFT4e:** necessary higher derivatives, modified seam/heat-exterior/energy compatibility, one genuine common finite N and full modified stress cone.
- [ ] **REC / WAVE / PHYS:** n-dependent coefficient recovery with independent repairs and smooth sum; mean/two-family oscillatory quadratic stress cancellation; corrected Cartesian uvw/p/f and dynamical/full-residual evidence.

Scoped gate: `current_original_native_radius_offsets_and_candidate_spatial_phase_executed`. Every global completion gate remains false. Keep actual original source enclosures rather than midpoint field values, implement the next bounded function output, mark only its achieved scope complete and preserve unrelated edits. Longer dependency detail remains in the [candidate density handoff](CURRENT_NATIVE_CANDIDATE_DENSITIES_2026_10_07.md).
