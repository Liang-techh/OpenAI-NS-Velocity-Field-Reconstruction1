# Local pressure expansion and pressure-independent obstruction

Added27 pressure coefficients in[-1,1]: nine Gaussian modes at (R²,Z²) centers{0,1.5,3}², width1.2, multiplied by the shrinking compact bump and tau^-1.01; each has three quadratic Bernstein time factors. Velocity, force, energy and core profiles are unchanged. Save/load retains the entire velocity candidate. The focused test checks unchanged velocity, compact support, reflection symmetry and serialization.

Dense fourth-power fitting reduced the independent sampled maximum only1.24417262 ->1.24361370. Structural probes pass; PDE fails .001. This is a small improvement, not a resolution of the axial/radial defect. Artifacts are local_pressure/. The pressure coefficient bounds were recorded before fitting in configs/constraints_local_pressure.json.

## Pressure-independent circulation

For fixed u and f, write A=u_t+(u dot grad)u-nu*laplacian(u)-f. On any closed spatial curve, integral grad(p) dot dl=0. Hence every smooth pressure satisfies sup_curve||A+grad(p)|| >= abs(integral A dot dl)/curve_length. This identity does not assume any finite pressure basis.

At t=.75, a meridional rectangle r in[.6,1.0],z in[.2,.6] has perimeter1.6. Numerically, its circulation is .57410265,.57231439,.57143677 for derivative steps .005,.0025,.00125 at fixed quadrature order96. At the finest step, orders32,64,96 agree to about6e-14. The resulting estimated lower bound is .35714798. Values are finite-difference/quadrature observations, not certified interval bounds; nonetheless this is strong evidence that pressure-only fitting cannot attain.001 for this fixed velocity/force.

The circulation operator has a manufactured polynomial test: u=(z,0,r²), whose meridional acceleration curl is2z and whose rectangular circulation is known exactly. That test passed. Reproduce with python -m openai_ns_reconstruction.constrained_pressure_circulation. Derivative and quadrature refinement are varied separately, and the data are saved in local_pressure/circulation.json.

## Required next representation change

Move to joint poloidal velocity and pressure optimization. Add compact divergence-free streamfunction corrections that vanish initially and at core probes, with explicit bounded coefficients and a well-resolved spatial guard. Use analytic parameter residual/energy derivatives where possible. Changing poloidal velocity may change convection and therefore the pressure-independent circulation, unlike pressure additions. Preserve the successful inner/outer swirl modes, fixed restricted force, core/energy constraints and final thresholds. Do not continue pressure-only basis growth as the main route.

Reproduce the pressure fit using constrained_outer_pressure.run(initial='artifacts/constrained/inner_swirl_pressure/candidate.json',output='artifacts/constrained/local_pressure',dense_cylindrical=True,fourth_power=True,local_pressure=True), then the independent validation CLI with configs/constraints_local_pressure.json.
