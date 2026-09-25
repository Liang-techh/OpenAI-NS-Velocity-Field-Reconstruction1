# A wider radial cone near the midplane, limited by axial extent

`moment_bridge_cone_region_search.py` surveys 45 radial–axial nodes on
the three-knot moment-closed bridge at each of `k=11` and `k=19`.
Only three coordinates pass the sampled physical full-residual cone on
**both** scales: `(y,eta)=(0.05,-0.2)`, `(0.05,0.2)`, and `(0.35,0)`.
The new midplane point has cone ratios `0.01647` and `0.01597` on the
two scales, respectively. Its sampled complete momentum norms are still
`4.73e5` and `1.84e9`; pointwise cone feasibility is not momentum
improvement.

`moment_bridge_cone_support.py --midplane` refines the radial band at
`eta=0`. The passing samples run from `y=0.27` to `0.37` at `k=11`,
and `0.28` to `0.37` at `k=19`. The best-centered half-width *between
passing samples* is `2.07e-4` and `1.17e-5`. A pulse with time half-width
`0.04 tau` would need radial half-width `sqrt(nu*0.04 tau)`, namely
`3.13e-4` and `1.95e-5`, for its radial diffusion time to match the
pulse half-width. The gap is now only `1.51–1.68` in width, far smaller
than around the original `y=0.05` source.

The axial follow-up at `y=0.325` is less favorable. Passing samples at
both scales run only from `eta=-0.05` through `eta=0.025`; the next samples
at `-0.075` and `0.05` fail the local cone. The best-centered physical
axial half-width between passing samples is `6.11e-5` at `k=11` and
`3.93e-6` at `k=19`. The same diffusion-balanced target is **5.11** and
**4.97** times wider. Even the optimistic brackets between the first
failing samples leave a factor of about three in axial width.

All cone evaluations here use a Gauss-12 radial stress primitive and
finite-difference velocity jets. Passing nodes do not certify a connected
open region in space or time, and this physical cone analogue is not the
paper's normalized leading-profile condition. The midplane point is a
useful candidate for **background redesign** because its radial cone is
broader on two scales, but it cannot yet support a three-dimensional
pulse with the chosen `0.04 tau` duration and diffusion balance. The
next correction should widen the axial cone or jointly shorten the pulse
and control the resulting temporal cutoff, then revisit amplitude,
pressure, mean flow, and radial moments as one residual-improvement cycle.
