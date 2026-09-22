# ST067 — source-mechanism audit, not a new optimized field

Task #1249. Manual user-authorized continuation; paused schedules were not restarted. Original ST066-B/W parameters, viscosity, fixed force and thresholds are unchanged. No candidate/default/visualization promotion.

## Source findings and conditional mathematical restriction

The OpenAI-hosted paper distinguishes an axisymmetric leading background from non-axisymmetric oscillatory corrections whose quadratic averages carry missing momentum flux. It uses height-dependent q(z,tau), an inner/annular/heat-exterior construction and an existential smooth compact force. Our entire ST064--66 candidate is axisymmetric, scales uniformly in height, and uses a pre-fixed two-parameter force. These are different architectures and different force-selection problems.

For a globally axisymmetric smooth exact solution with regular axis, bounded initial Gamma=r*u_theta, and integrable sup torque, the scalar swirl equation gives

`||Gamma(t)||_inf <= ||Gamma(t0)||_inf + integral ||r*f_theta||_inf dt`.

The particular nonvanishing infinite-scale law r~tau^(1/2),u_theta~tau^(-1/2-h),h>0 requires Gamma~tau^(-h) unbounded. It is incompatible with those whole-field-axisymmetry and bounded-torque assumptions. This is NOT a general NS regularity result or a finite-ladder impossibility proof. It does not apply unchanged to the source's non-axisymmetric full correction. The current finite B/W ladder did NOT exceed the numerical global torque budget; W at k6 has located |Gamma|max .492268 versus estimated budget .582584. Both are numerical, not interval-certified suprema.

## Actual frozen-field diagnostics

At k6, 96-order physical volume quadrature:

|Field|Radial residual L2|Angular residual L2|Axial residual L2|Full L2|Angular fraction of squared error|
|---|---:|---:|---:|---:|---:|
|ST066-B|32.24419255|62.44538186|39.55825520|80.64718981|59.9546%|
|ST066-W|47.51036255|77.98208200|51.97590702|105.07109295|55.0837%|

Axisymmetric pressure cannot change the angular residual. 40/64/96 quadrature levels are retained; no continuum error enclosure. At k6 the declared central observation cylinder alone still has residual L2 B18.59475328/W23.76966587, so annulus-only pulses cannot by themselves fix the interior. The observation cylinder is an autonomous diagnostic region, not the source core boundary.

Implemented Gamma identity including r*R_theta and independently checked Cartesian finite-difference ladders at k0,2.7,6. Fine relative torque differences <=4.19e-8. These are operator-consistency errors, NOT PDE errors.

## Actual implementation

`source_coordinates.py` implements q-z^2*q^(2h)=tau, eta=z/q^(1/2-h), X=r^2/(2q), with safeguarded inversion and first derivatives. This is a unit-viscosity coordinate building block, not a reconstructed velocity or a numerical solution to the source profile equations. For physical nu=.01 the separate viscosity transformation must still be used. At source z=.4, tau=.0078125 gives q=.16495499 rather than tau; q remains nonzero as tau tends to zero away from the origin.

The offline package also contains a radial-only inverse-stress diagnostic and a pointwise PSD minimum-trace completion. The inverse has nonzero outer radial traces: that gauge cannot satisfy zero stress at both boundaries. Axial transport/general stress freedom can change this conclusion. The PSD matrix is NOT a divergence-free wave field, and is never substituted for a freely chosen force or claimed to cancel the full NS residual.

## Other primary work inspected

- Duraiswami arXiv:2609.17642: directly relevant similarity-coordinate numerical exploration; mapped Chebyshev grids, pressure gauge, Newton/continuation and axis series. Porous-wall and moment-matching results are not a global compact full NS candidate. Raw author-code download attempts failed DNS, so no external solver execution or code-copy license claim.
- Chen--Hou arXiv:2305.05660: approximate-profile stability and rigorous error-control methodology, for inviscid/boundary problems, not our viscous whole-space benchmark.
- Hou arXiv:2405.10916: two-scale dynamic rescaling in generalized systems; do not transplant its dimension/viscosity modifications into the original fixed equation.
- Ozanski--Zajaczkowski arXiv:2405.16670v2: swirl equation (1.7) and maximum-principle role. Their conditional regularity theorem is not generalized here.
- OpenAI official NavierStokesAndEuler README describes Lean formalizations. Neither that build nor all166paper pages were independently checked this round.

Source links:
https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
https://github.com/openai/NavierStokesAndEuler
https://arxiv.org/abs/2609.17642
https://arxiv.org/abs/2305.05660
https://arxiv.org/abs/2405.10916
https://arxiv.org/abs/2405.16670

## Execution and delivery

21 focused tests passed2.88s. Independent clean-package tests21passed3.18s; Wk6 whole component/region dictionary and both models' three-resolution stress diagnostics reproduced exactly. All source input hashes verified. The full Gamma-extremum search was executed in the main work directory, not repeated in the clean smoke test. Two early foreground read-only audits reached tool time limits; completed rerun used unchanged field/numerical operations. Later source documentation and completed-output overwrite protection did not change mathematics.

Actual input archive SHA25645fad36ec3b4e946f94182fb24dabe4cc5df77dc7e6646a11e0d9d9a5b396153. The coordinate source and this record are on the branch. Complete diagnostic sources, minimal original evaluators, frozen B/W profiles, identities and numerical evidence are delivered in NS_ST067_Mechanism_Audit.zip, not claimed all uploaded here. Full original fitting histories remain in the previous complete ST066 archive.

No optimization, new candidate, native MATLAB, cloud CI, external numerical solver or Lean run. No residual decrease claimed. Original full1e-3 target remains UNMET.

Recommended next construction: solve a compatible inner profile in source coordinates; match its heat exterior and radial/axial moments; then attempt realizable non-axisymmetric stress-producing waves with independent three-dimensional residual checks. Keep the old fixed-force benchmark separate from any newly registered existential smooth-force construction.