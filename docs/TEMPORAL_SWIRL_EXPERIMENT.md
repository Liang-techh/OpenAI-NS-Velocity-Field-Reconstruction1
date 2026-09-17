# Temporal swirl correction: scope and experiment

The latest diagnostic identifies azimuthal momentum as the dominant residual channel. This experiment introduces 15 bounded coefficients while freezing the parent velocity parameters, pressure and restricted force.

Let V=(-y b,x b,0) be the existing smooth compact outer swirl basis, with b=0 in the core cylinder. For each q=(rÂ²/4)^i(zÂ²/4)^j, define W=V(q-m), where m=J(Vq)/J(V) and J integrates axial angular momentum. Each mode is axisymmetric pure swirl, so its divergence vanishes identically: the x and y derivatives cancel because b(q-m) depends only on rÂ² and z. If m is evaluated exactly, J(W)=0 by linearity. In code m is computed by quadrature, so moment cancellation is numerical evidence with quadrature error, not a certified identity for the stored constants.

The time multiplier is a cubic Bernstein polynomial B_k(2(t-.25)), k=1,2,3. Every multiplier vanishes at .25, so the correction preserves the initial field exactly. Spatial support and the core exclusion are inherited from V. The existing parent moment model remains unchanged because the additional zero-moment modes are outside that model; no quadratic-time extrapolation is applied to the new cubic correction.

Energy is not invariant under this addition and must be checked after fitting. Local momentum balance is not implied by zero total moment. The experiment uses all three momentum components in its fitting objective, so gains in swirl cannot hide a larger radial centripetal residual.

Run `python -m openai_ns_reconstruction.constrained_temporal_optimize`, then the independent validator with configs/constraints_temporal_swirl.json and the saved candidate/training paths. Training and validation use separate samples and differential operators; the validation set remains a development set and final blind acceptance remains pending.

## Results (development, not acceptance)

Uniform-time fit: 160 calls, training mean-square residual 0.10515165 -> 0.06328958. Independent maximum 1.95654830 (parent 1.94296483), worst at t=.75; energy [0.89923205,1.0], core drift .04949132, sampled structure passes. Cached quadratic residual agrees with the direct FD implementation within 3.83e-12.

Tail-adaptive comparison: selected 256 largest parent residuals at t=.748 from a new 8192-point training pool (seed 20261017), appended to 2048 baseline points. 176 calls, training mean-square 0.33565053 -> 0.09721559 on this different distribution. Independent maximum worsens to 2.44784120; energy [.78842774,1.00115946], structure samples pass. This is a preserved failed comparison, not the selected candidate.

Neither result meets .001. Retain continued_pressure as the best sampled-maximum reference. Temporal freedom improves average error, but these globally supported polynomial spatial corrections redistribute endpoint error and do not resolve the maximum. Next: identify the new endpoint hotspot and use localized spatial swirl modes or a direct azimuthal evolution solve, rather than expanding these coefficients' bounds. Eleven temporal/cache tests passed; worker also ran 17 temporal/outer/tensor tests.
