# Outer residual source and missing stress realizability

The saved short DAE trajectory can now be reconstructed by
`SavedOuterDAEField` without rebuilding its moment-map samples. The
reconstructed ODE solution matches every stored swirl trajectory node
within rtol=1e-9, atol=1e-11. It remains an experimental field restricted
to k in [13,13.1], not a globally accepted export.

At k=13.05, a radial/axial grid finds a full momentum peak 1.34590e8 at
eta=0 and bridge fraction y=0.4265, near the edge of the first outer bump.
The physical point is (0.001004093753,0,0). The axial component terms are:

| Term | Axial value |
| --- | ---: |
| Physical time derivative | 2.30275e6 |
| Advection | -2.81892e6 |
| Pressure gradient | 1.49087e4 |
| Negative viscous Laplacian | -1.34088e8 |
| Full axial residual | -1.34589e8 |

The decomposition sums back to the independently computed full momentum.
This sampled peak is dominated by axial viscous curvature, not an omitted
large pressure cancellation. It motivates redesigning the outer radial
shape before simply increasing optimization iterations. It does not prove
which new profile will satisfy the complete PDE constraints.

A split 64-point stress-primitive screen at outer-window fractions
0.5 and 0.75 and eta=-0.2,0,0.2 passes 0/6 cone nodes. Every node has
positive target_dot_N (1924 to 5488), contrary to the required negative
sign. Two nodes also have negative lambda_squared. Cone ratios of zero
at those nodes reflect the sqrt(max(lambda_squared,0)) convention and
are not passes. The earlier 9/9 inner-support result therefore cannot
justify placing the paper's pulse correction in these outer regions.

Next broaden/reshape the outer compact streamfunction modes to reduce
axial viscous curvature while restoring moments and retaining the inner
support gap. Evaluate the complete momentum cost, then repair the outer
mean/pressure stress direction before attempting a Section 7 pulse there.
The frozen pressure is only a restriction of these experiments, not of
the paper. Preserve time-interval constraints as profiles change.

This is a diagnostic grid and six sampled stress nodes, not a global
maximum or continuum impossibility proof. No full momentum/volume-L2,
finite-energy, or recursive contraction requirement is satisfied.

Reproduce: `python experiments/root_st073/midplane_outer_residual_source.py`.
Reusable experimental evaluator:
`from saved_outer_dae_field import SavedOuterDAEField` then
`field.fields(points, remaining_time)` within the saved interval.
