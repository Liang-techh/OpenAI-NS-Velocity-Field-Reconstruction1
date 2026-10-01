# C120-F13: full radial reference-annulus relaxed cone

The repaired reference-annulus field now has a directed cone certificate over the entire radial domain `x=R/Rm in [1,2]` and the accepted local axial family `Z in [0.49,0.51]`. All 40 gap-free radial cells certify the relaxed `kappa<=2` branch. This supersedes the five-point radial sampling limitation of C120-F12 for this annulus only. It does not establish strong admissibility or whole-axis coverage.

## Whole-cell bounds

The grid includes all three bump support boundaries and centres, and subdivides each resulting interval into four exact rational cells. A positive cumulative weight is zero before its support, is the frozen full-support weight after its support, and is bounded by `[0,upper(full)]` while crossing the support. Positivity makes this enclosure valid at every radial point, including partial supports. Signed linear rows negate the positive interval only after enclosure. Quadratic coupling and implicit C1 control intervals remain intact.

For the normalized bump `beta=exp(-1/q)/(r*N)`, `q=1-t^2`, `r=1/40`, the entire support has

`0 <= beta <= exp(-1)/(r*N_lower)`.

Since `q^-2*exp(-1/q) <= 4*exp(-2)` on `0<q<=1`, its derivative satisfies

`|beta_x| <= 8*exp(-2)/(r^2*N_lower)`.

Left-half support derivatives are nonnegative and right-half derivatives nonpositive; outside/support endpoints they vanish exactly. Interval negation preserves the directed lower endpoint independently of ambient scalar precision. These analytic bounds cover all radial points, rather than estimating maxima from dense samples.

The physical cell evaluator uses the same five source defects, implicit controls, moment scales, pressure `P0+Mp`, velocity reconstruction and stress formulas as the preceding point field. All radial powers and physical radius factors are interval evaluated over each full cell. Dependence may widen bounds, but no correlation is removed by midpoint substitution or fitted targets.

## Evidence and scope

- `lei_ren_part1_paper_interval_five_bump_cell_map.py/.json`: analytic primitive and bump envelopes, independent scalar integral/value/derivative containment fixtures, and point reduction.
- `lei_ren_part1_paper_interval_reference_annulus_cone_atlas.py/.json`: all 40 full cells, complete `[1,2]` coverage, normalized shear/stress and cone margins. No unresolved cell remains.
- `lei_ren_part1_paper_interval_reference_annulus_cone_atlas_check.py/.json`: 30 physical coefficients reduce exactly to the original adapter on a degenerate radial cell and are contained in a support-crossing full cell. This is a translation regression, not independent quadrature; independent physical quadrature remains in the C120-F12 fixture.

The cone statement is conditional on the accepted fixed construction parameters and the preceding local implicit family certificate. Original parameter remainders, C2/higher implicit controls, whole-axis matching, the callable `R110..Rm` connecting field, shear modulation to the required final strong stress, flatten/heat collar/exact exterior, finite energy of the complete field, flat remainder, genuine temporal recursion, oscillatory corrections and independent full Cartesian NS residual remain incomplete.

The next structural step is the actual `R110..Rm` long reshape and axial restoration field, with controlled cumulative moments and the same analytic pressure. An endpoint defect enclosure already exists; it must not be presented as the complete connecting field.
