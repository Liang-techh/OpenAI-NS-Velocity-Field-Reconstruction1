# Remote meridional shape under both radial moments

This follows the radial-moment and stress-realization separation in the
[OpenAI Navier–Stokes construction](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).

The moment-matched remote `[4,16]` field now admits a second streamfunction
shape, `N(X)=B(y)(2y-1-alpha)`, alongside the original `B(X)`. Its axial
velocity contribution is proportional to `dN/dX`; because `N` vanishes at
both patch endpoints, this preserves the zero radial integral of the axial
correction. The radial velocity is differentiated from the full
streamfunction, preserving incompressibility. At each eta/time slice the
coefficient of `B` is re-solved from the quadratic kinetic moment. The
angular moment is unchanged because the new mode has no swirl component.

At `tau=.0084`, six meridional amplitudes `-5,-3,-1,0,1,3` are screened at
`eta=0,.2,.35` and `X=6,10,14`. The new shape is materially more effective
than the swirl-only null mode: at `X=10,eta=0`, the stress projection onto the
local cone normal changes from `+3.50` at amplitude zero to `-3.76` at
amplitude `-5`. However its cone ratio at that point is `4.91`, above the
strict limit `1`. At amplitude `-5`, the same projection is negative at
`eta=.2` but the ratio is `9.20`; at `eta=.35` the projection remains
positive. None of the sampled positions passes the entire local analog cone.
The direct moment defects over these slices remain below `9.88e-4`.
At `eta=.2`, amplitude `-5`, and `X=6,10,14`, sampled finite-difference
divergence has magnitude at most `3.56e-9`, while full momentum residual
norms are approximately `2100,7135,1664` respectively. The shape fixes an
orientation problem locally but not the PDE residual.

Increasing the negative amplitude is constrained by the kinetic moment:
the quadratic discriminant at amplitude `-8` is negative at the three
sampled eta slices, so this one-mode branch has no real moment-matching
coefficient there. This suggests searching a **multi-mode meridional
streamfunction**, including eta-dependent coefficients, rather than simply
increasing the magnitude of this one mode. The screen is not a normalized
Appendix-A cone calculation or a continuum exclusion, and no oscillatory
stress, finite-energy exterior, full momentum or force condition is met.

Run `python experiments/root_st073/remote_meridional_shape_screen.py`.
Results are in `experiments/root_st073/remote_meridional_shape_screen.json`,
marked `accepted: false`.
