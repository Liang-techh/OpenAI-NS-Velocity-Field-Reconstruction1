# Inner-collar swirl and pressure continuation

Added six compact pure-swirl rings, each with five nonzero quintic time coefficients. Retain the old60 coefficients and initialize the30 new ones at zero. The guard is1-standard_cutoff(r²/.02), exactly zero for r<=.1. All core probes have radius at most.086603, so this allows changes in the previously frozen annulus while preserving the probes and initial velocity. Coefficients stay in[-1,1]. Each raw mode is compensated against the existing outer angular-moment basis; forcing is unchanged.

The focused test verifies exact embedding of the prior candidate, a nonzero change in the formerly frozen region, unchanged core/initial values, numerical zero added moment at construction quadrature, divergence, and serialization. It passed with warnings treated as errors.

## Results

The point(r,z,t)=(.32,.8857142857,.75) previously had residual norm1.15115944 at h=.00125. The new inner_swirl candidate reduces it to.62614070. Full random-sampled maximum remains1.26462981, structural probes pass. A separate64x64 cylindrical grid locates a new final-time peak dominated by axial momentum: axial component1.15338 near(r,z)=(1.0154,.68794). Axisymmetric pure-swirl corrections leave the axial momentum equation unchanged when poloidal velocity, pressure and force are fixed.

Extended the existing bounded18-term pressure fitter to handle this wrapper and dense collocation, with an analytic Jacobian for a fourth-power residual objective. The resulting inner_swirl_pressure candidate lowers the full independent sampled maximum to1.24417262; sampled energy/core constraints remain satisfied. The training mean-square increased .07781193 -> .07892784 because this run optimizes fourth power, so it must not be described as a mean-square improvement. All PDE thresholds still fail .001.

## Quadrature limitation

At final time, added angular moment is approximately zero at construction order96 but -2.77825e-5 at144 and -2.81608e-5 at192. Therefore the stored mode compensation is not exact continuum torque closure. Energy estimates are .70897679, .70898116, .70898137 respectively. This refinement is recorded in inner_swirl_pressure/quadrature.json; it does not certify uniform time or spatial behavior. Higher-order mode normalization should accompany final numerical acceptance, without claiming that it resolves the much larger local PDE residual.

## Reproduction and next work

Run constrained_temporal_optimize.run(inner=True, fourth_power=True, dense_cylindrical=True, analytic_jacobian=True, coefficient_initial='artifacts/constrained/quintic_swirl/candidate.json', grid_resolution=48, output='artifacts/constrained/inner_swirl'). Then run constrained_outer_pressure.run(initial='artifacts/constrained/inner_swirl/candidate.json', output='artifacts/constrained/inner_swirl_pressure', dense_cylindrical=True, fourth_power=True). Validate with configs/constraints_inner_swirl.json and the standard validation CLI.

Next work: expand pressure's spatial representation or jointly vary poloidal velocity, targeting the axial/radial residual. Preserve the now-enabled inner swirl modes, fixed restricted force, dense sample coverage and all acceptance thresholds. Current pressure-only fitting is bounded but its global low-degree basis still leaves a large defect. Final blind audit remains pending.
