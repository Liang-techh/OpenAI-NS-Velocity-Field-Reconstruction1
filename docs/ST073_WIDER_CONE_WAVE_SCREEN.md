# Wider axial support does not close the scale-recursive wave residual

The three-knot mean repair opens the sampled positive axial cone point at
all integer scales `k=11..19`. This experiment tests the immediate
wave-level consequence: increase the exact-curl wave's axial half-width
by `1.4x` and its diffusion-scaled time half-width by `1.4^2`, using the
moment-repaired mean. The center stress target is recomputed from that
mean, and the two positive physical covariance weights are refitted.
The Kelvin wavevectors and polarizations remain those of the prior source.

At pulse center, full Cartesian momentum is evaluated on the same
21 radial-axial nodes and 16 angles for both widths, including nodes
near each cutoff. With wave amplitude multiplier `0.1`, the maximum
residual falls from `2.76e9` to `2.17e9` at `k=11` and from `1.02e13`
to `8.07e12` at `k=19`: reductions of `21.3%` and `20.9%`.
The widened-wave peak still grows by a factor `3715` across eight
halvings and remains thousands of times larger than the sampled mean
residual. The wider cone therefore reduces a cutoff cost but does not
produce residual cancellation or a recursive step.

The sampled cone pass has a narrow margin and is not a proof of support
at every radial point, time, or intermediate continuous scale. The next
implementation must solve for a spatially varying transported amplitude,
normal pressure, and mean/stress correction on the supported region,
then require direct full-residual reduction on adjacent scales and pulse
interior times. Merely widening this frozen wave is insufficient.

Reproduce with
`python experiments/root_st073/midplane_wider_cone_wave_screen.py`.
The JSON file records the common grid, refitted weights, widths, and
both residual budgets.
