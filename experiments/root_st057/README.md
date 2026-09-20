# ST057 — full velocity redistribution and harmonic compatibility

Task #847. Start from the actual ST056-J2 raw field (SHA256 `31961e6c7ff58358954967308fdaa3dcda3d52dfe05c796e35d898203bf0375d`) and source PR #754 at `4df0221b428441b4d92afe2f994b024ee66b2146`.

## Implemented in the current bounded round

The new optimizer permits all 1728 velocity coefficients and 864 pressure coefficients in the SAME pre-existing compact axisymmetric family. Initial velocity is no longer frozen. Initial energy is still exactly normalized by the existing representation, with the normalization differentiated in the adjoint. Original viscosity 0.01, time [0.25,0.75], smooth compact velocity AND pressure, independently prescribed restricted force, support cylinder, parameter bounds and original 0.001 momentum gates are unchanged. Force coefficients stay fixed.

A separable tensor evaluator and analytic adjoint include all three nonlinear momentum components. Coupled Jacobian whitening includes velocity-pressure cross terms. Seven harmonic-gradient weak moments (degrees 2 through 8) are included at 13 Gauss time nodes plus endpoints. Exact polynomial Gram matrices normalize those moments; velocity integrals are numerical quadrature, not interval bounds. Additional pressure, core, shear, bias and vorticity-distribution restrictions are explicitly autonomous finite training safeguards.

Initial broad trial C completed 1200 iterations, but every candidate offered by its registered training selection failed at least one structural/no-worsening safeguard. Its output therefore falls back to the original mathematical parent; it is NOT a new improved candidate. Lower unconstrained training residuals are not being promoted.

H and K use stronger structure penalties and retain only the original isolated probe's raw velocities at 25 training times through a numerical nullspace. That small probe lock does NOT freeze the entire initial spatial field; the global normalization may still scale the probe. H/K differ in harmonic penalty weight. They are followed, if needed, by a separately registered TRAINING-only sequential half-space feasibility restoration in the coupled numerical metric. This is not exact constrained pressure elimination or a continuous minimax solver.

## Current evidence boundary

Eight focused derivative/normalization/operator/structure tests passed locally, including full-objective and independent Cartesian residual checks. Two attempted starts were rejected by the finite-difference gradient calibration at a hinge boundary; finer-step and off-boundary checks resolved the calibration issue without relaxing its error tolerance. Logs are retained. No complete PDE or full historical/Lean test claim.

The training checkpoints are not yet final scientific results at this commit. Fresh original-validator seeds `9205791` and `9205792` are reserved for after every selected coefficient vector is frozen. No holdout-driven coefficient tuning is permitted within this round. Full raw data, histories, executed source and final replays will be retained in the deliverable; this note alone does not substitute for them.

No main/default/publication/other-agent/schedule changes. `pde_validated=false`, `source_correspondence_verified=false`, `paper_exact=false`, `blowup_proved=false`.
