# Coupled bridge bubble: cone gain versus momentum cost

The extended ST073 candidate's radial bridge has large momentum residual
and fails the [OpenAI paper's](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
Section-4 relaxed-cone analog. An azimuthal bubble alone changed the
swirl slope but left the meridional-shear term too large. This experiment
adds an axisymmetric poloidal streamfunction bubble inside the moving
radial bridge. Its radial shape and first four endpoint derivatives vanish;
the velocity is obtained from `u_r=-psi_z/r`, `u_z=psi_r/r`, including the
moving-coordinate and axial-cutoff derivatives. It is analytically
divergence-free and preserves the inner and exterior velocity traces.

The `curvature` shape also has `psi=psi_r=0` at `X=1`. At this point it
changes `U_X` without directly changing `U` or the cumulative `M=integral U
dX`, allowing a focused meridional-shear test. Combined with an outer-biased
swirl bubble of amplitude `1.1`, a poloidal curvature amplitude `-2.5`
makes the two sampled `X=1, eta=.2,.3` relaxed-cone inequalities positive
at `k=3,5.5,6`. At `k=5.5` their margins are about `0.689` and `0.683`.
This is a local diagnostic only: nearby bridge samples fail, including
`X=.5859375,.75,1.25` on both axial slices.

More decisively, fourth-order physical-coordinate momentum residuals at
eight bridge points rise from a baseline sampled maximum of `3946.7` to
`43927.3` for this coupled candidate. At `X=.5859375, eta=.2`, the
candidate's axial viscous term is about `41858`; only about `1032` comes
from axial second derivatives. The large new defect is mainly **radial
curvature of the poloidal correction**, so changing only the temporal
amplitude derivative cannot plausibly remove it without introducing an
extreme time scale. See `coupled_bridge_screen.json` and
`coupled_bridge_hotspot.json` in `experiments/root_st073`.

This candidate is rejected for the full momentum target. The useful result
is the coupled constraint: a bridge repair must control both the relaxed
stress cone and full physical momentum throughout an annulus, not at an
isolated point. The next mean-flow solve should use a wider or dynamically
evolved poloidal correction with explicit radial-viscosity control, keep
the five outgoing radial moments, and only then apply the paper's
nonaxisymmetric oscillatory stress mechanism. No current result proves a
continuous cone, complete residual maximum, volume L2, or smooth forcing
gate below `1e-3`.

The follow-up minimum-curvature and three-mode screens reduce the sampled
momentum penalty but still fail the full bridge requirements; see
[`ST073_BRIDGE_CURVATURE_SEARCH.md`](ST073_BRIDGE_CURVATURE_SEARCH.md).
