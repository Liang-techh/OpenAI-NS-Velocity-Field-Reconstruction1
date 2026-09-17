# Joint poloidal and pressure experiment

Added27 velocity coefficients from streamfunctions psi=r²*z*F(r²,z²)*B_k(2(t-.25)), with nine Gaussian centers and cubic Bernstein k=1,2,3. The induced Cartesian field is (-x*(F+2z²*F_v),-y*(F+2z²*F_v),2z*(F+r²*F_s)); these analytic streamfunction derivatives make the velocity divergence free. Optimize those27 coefficients jointly with the existing27 pressure coefficients, all bounded in[-1,1], with fixed restricted force and existing swirl parameters. Reuse residual/energy caches and analytic parameter Jacobians.

Two core-preservation variants were tested. The cutoff variant uses B(.01/r²), zero inside r<=.1. The smooth anchored variant uses q=(r²-.01*tau)²/(r²+.01*tau+.03)². On the moving core cylinder q=q_s=0, so the streamfunction and its first spatial derivatives vanish and velocity probes are unchanged. This does not assert zero velocity gradients or a zero correction throughout a core neighborhood. Initial velocity is retained by the time basis. Both variants are compact and axis-regular. Two parametrized tests verify core/initial preservation, divergence at interior fixtures and serialization.

## Results

The cutoff fit took27 calls and reached ftol. It reduced the existing circulation-bound estimate .35714798 ->.17468923, but worsened the full sampled maximum to1.31293603, with an axial-dominated peak near r=.1164 in the cutoff transition. Preserve it as a failed maximum-residual comparison.

The anchored fit took29 calls and reached ftol. Maximum at the standard h=.005 is1.00633693, about19.1% below the prior1.24361370. Energy range[.71778740,1.0], core drift .04949132 and core signs pass. The same loop estimate is.15954425 (numerical, not interval-certified), showing that poloidal changes address the defect that pressure alone could not.

## Divergence and refinement are separate from core/energy probes

At h=.005, the numerical divergence maximum is.00364145, above1e-5. This large discretization error was already present in the prior inner-swirl candidate; it must not be hidden by reporting core/energy success as all structural gates passing. For the fixed anchored candidate, refinement gives:

| Step | Momentum sampled max | Divergence sampled max |
| --- | ---: | ---: |
| .0025 | .97049888 | .0004612053 |
| .00125 | .95958807 | .0000499532 |
| .000625 | .95664924 | .0000043161 |

The finest sampled divergence maximum is below1e-5, while momentum remains far above.001. These are existing development samples, not a final blind audit or a uniform proof. No validation threshold or original configuration was relaxed. Standard validation remains failed; refinement is supplemental evidence.

## Reproduction and next work

Run python -m openai_ns_reconstruction.constrained_poloidal_optimize for the cutoff comparison. For the working reference call constrained_poloidal_optimize.run(output='artifacts/constrained/poloidal_anchor',anchor_core=True). Validate with configs/constraints_poloidal.json. Artifacts include the candidate, training, original validation, circulation and derivative refinement.

Next: retain the smooth core anchor and jointly free the successful swirl and poloidal coefficients with pressure, since changing poloidal transport changes the azimuthal equation. Compare full residual and dense grids, not circulation alone. Preserve initial/core constraints, bounded force, coefficient bounds and energy gates. Address derivative resolution of the narrow swirl guard during final convergence studies.
