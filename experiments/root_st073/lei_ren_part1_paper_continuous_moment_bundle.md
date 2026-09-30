# Same-candidate five-moment and stress bundle

The shared profile exposes angular_moments_jet, linear_axial_moments_jet, quadratic_moments_jet and pressure_moments_jet. The quadratic API dispatches to the partial provider before Rv and the complete provider at/after Rv. Pressure is anchored once to the actual inner datum. The bundle assembles all five actual moments, their Z jets and analytic velocity jets without forcing terminal values to zero.

The integrated post-Rv probe was replayed after API wiring at Z=.3 at flattening (Rv+50 in log radius) and heat (Rtail+4). Both runs produced finite outputs but enormous total-stress/shear ratios. This is evidence of failed matching, not admissible stress or a PDE residual certificate. Incoming/pulse/end components have their own bounded receipts; the full bundle has not been independently replayed across those regions.

Partial axial component checks: Rv agreement 0, accepted radial/Z integrand maximum 4.70e-17, isolated end-square/Z checks 2.49e-22. Whole-sum end stencils report relative error 1 because the accumulated pulse dominates; those diagnostics remain recorded. Pressure checks include the actual inner anchor, separate signed bump and heat corrections; most huge-radius whole-sum differences are precision limited. Quadrature and arithmetic errors are not enclosed.

The pressure-tail audit identifies a nonzero combination 2(1+delta) Z P - (1-Z^2) PZ about -0.019363 at both probe points. After axial support, N_z is exactly the pressure contribution to I_z, and saved I_z/N_z agree at serialization precision. Subdominant moment contributions are not resolved by that comparison. Restore the actual inner-seeded cumulative pressure and angular targets together; changing a Z-dependent datum alone changes the momentum equation.

No claim of full moment closure, finite energy, stress cone, remainder decomposition, oscillatory residual or scale recursion is made.
