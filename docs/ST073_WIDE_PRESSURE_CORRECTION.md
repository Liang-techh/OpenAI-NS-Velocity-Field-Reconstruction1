# Compact pressure correction on the widened ST073-V join

The widened order-14 ST073-V core can be joined kinematically to the pure
swirl heat exterior, but the join has late complete-momentum residuals of
several thousand. This experiment adds six axisymmetric pressure functions
supported strictly inside the ratio-4 radial bridge. Each uses
`256 y^4(1-y)^4` times one of two radial and three axial polynomials, so
the pressure and its first three radial endpoint derivatives remain matched.
The velocity and its solenoidal construction do not change.

At `k=5.5` the six pressure amplitudes were fitted to the full unforced
Cartesian momentum on 15 bridge points. Re-evaluating the resulting field
with the independent fourth-order finite-difference operator gives:

| Sample | Baseline maximum | Corrected maximum |
| --- | ---: | ---: |
| 15 training points | 7303.66 | 5942.90 |
| 6 disjoint spatial points at the same time | 7006.96 | 3902.26 |
| The same disjoint locations at `k=5.25` | 5414.92 | 3023.96 |

The analytic pressure-gradient prediction and complete finite-difference
re-evaluation agree within `3.0e-7` as absolute residual vectors on these
samples. This is a real local decrease on held-out points, but it is still
millions of times above the required `1e-3` pointwise threshold. Since an
axisymmetric pressure has no azimuthal gradient, angular momentum defects
remain: the corrected angular component is still above `1e3` on each set.
The fitted pressure also changes the radial pressure force without a matched
change to swirl/meridional flow. It is an exploratory transition profile,
not the paper's completed centrifugal-pressure or five-moment join.

The next step is to add endpoint-preserving solenoidal streamfunction and
swirl modes and solve their coupled momentum/pressure response over space
and several scales. Axial localization, finite global energy, nonaxisymmetric
stress correction, and full-domain acceptance remain open.

Run `python experiments/root_st073/wide_pressure_fit.py` from the repository
root. The coefficients and sampled results are in
`experiments/root_st073/wide_pressure_fit.json`. The candidate is callable
through `wide_pressure_fit.load_pressure_candidate().fields(points, tau)`.
