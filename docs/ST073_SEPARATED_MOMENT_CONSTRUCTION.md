# Separated radial directions close the sampled outer moments

The previous 24-mode bridge had an eight-row physical outer-moment
Jacobian with sampled condition ratio about `3.0e6`, and fitting those
moments increased the momentum residual without closing them. Following
the separated-bump idea of the [OpenAI paper's Appendix A, Lemma A.1](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
`separated_moment_modes.py` adds two disjoint flat radial windows in the
bridge, `y in (.12,.38)` and `(.62,.88)`. Each has azimuthal and poloidal
streamfunction modes with even and odd axial factors. Two independent
coefficient groups are connected by a smooth, flat-ended function of
logarithmic time between `k=11` and `k=19`: 16 new amplitudes in all.
The poloidal correction is built as `u_r=-psi_z/r`, `u_z=psi_r/r`, so it
is spatially divergence-free by construction; the azimuthal correction
also has zero divergence under axisymmetry. The windows leave the earlier
cone hotspots around `y=.05` unchanged.

`separated_moment_fit.py` fits the eight physical tangential outer-moment
components at `k=11,19`, `eta=±.2`. Its normalized eight-by-16 moment
Jacobian has singular values from `155.69` to `0.2656`, a sampled
condition ratio about `586`. The unregularized moment solve reaches a
largest normalized quadrature residual of `1.8e-12`. Independent direct
finite-difference radial integration in `separated_moment_validate.py`
finds outer-stress norm ratios to the original cone-aware field below
`1.1e-9` at the four trained slices. At the untrained middle scale
`k=15`, the ratios remain about `4.2e-4`; they are small but not zero.

| Field | Maximum full momentum, `k=11` | `k=15` | `k=19` |
| --- | ---: | ---: | ---: |
| Before separated repair | `3.264e5` | `2.069e7` | `1.312e9` |
| Exact sampled outer-moment repair | `1.373e6` | `8.551e7` | `5.324e9` |

Thus this is a usable *moment-control direction*, not an accepted
velocity field: sampled full momentum worsens by roughly fourfold and
still scales upward with dyadic refinement. A balanced momentum/moment
fit either leaves `46%` of the normalized outer-moment defect (weight
`5`) or, when weighted more heavily, keeps about `0.23%` while still
worsening momentum. The local physical cone analogue at `y=.05` is
preserved because the new modes vanish in a neighborhood of it; this
does not give a continuous cone certificate on the new windows.

The repair covers two physical tangential integrals, not the paper's
five normalized leading-profile moments or its uniform parameter
families. The calculation uses finite quadrature and the adapter's
experimental extended `k <= 20` interval. It does not establish exact
continuum moment closure, the paper's recursive residual improvement,
finite total energy, or the required full-domain max and volume-L2
below `1e-3`. The next coupled step must retain these separated moment
directions while adding momentum-canceling wave/mean corrections across
each scale, rather than treating moment closure as the final correction.

## Constrained momentum recovery and a third radial window

`separated_moment_constrained.py` minimizes complete Cartesian momentum
over the original 16 directions while enforcing the eight sampled outer
moment equalities. An L4 residual objective reduces the two-scale peak
by about `2%–3%` compared with the moment-only coefficients and also
improves `k=13,17` holdouts; the remaining peak is still about four
times the pre-repair field.

`separated_moment_three_window.py` adds one independent flat bump on
`y in (.40,.60)`, leaving all three windows disjoint and the `y=.05`
cone neighborhood untouched. The resulting 24 amplitudes allow an
exact sampled moment solve followed by constrained L4 momentum
minimization. This substantially reduces the cost of moment closure:

| Field | `k=11` peak | `k=15` peak | `k=19` peak |
| --- | ---: | ---: | ---: |
| Pre-repair cone-aware | `3.264e5` | `2.069e7` | `1.312e9` |
| Two windows, exact sampled moments | `1.373e6` | `8.551e7` | `5.324e9` |
| Three windows, exact sampled moments + L4 momentum | `4.897e5` | `3.033e7` | `1.878e9` |

The three-window constrained optimizer reports a largest normalized
training-moment residual `5.35e-14`. Independent direct radial
integration gives outer-stress norms below `1.0e-9` of the original at
`k=11,19`. Time holdout momentum peaks are `3.862e6` at `k=13` and
`2.382e8` at `k=17`; both are about `64%` below the corresponding
two-window constrained peaks. The off-knot outer-stress ratios to the
original are about `0.45%` at `k=13`, `0.023%` at `k=15`, and `0.42%`
at `k=17` (`separated_moment_three_window_transfer.py`). Thus endpoint
closure does **not** imply exact closure throughout the scale interval.

The three-window absolute peak still grows with exponent about `1.488`
per halving from `k=11` to `k=19`, and is `43%–50%` larger than the
pre-repair candidate. The next scale-recursive step needs moment closure
at additional scale stages, a controlled interstage transfer, and
nonaxisymmetric wave/mean corrections to cancel the large remaining
momentum. None of these finite-grid numbers is a continuum acceptance
or a proof of the paper's recursive mechanism.
