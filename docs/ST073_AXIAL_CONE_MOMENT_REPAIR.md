# Two-scale axial-cone change with sampled moment repair

This experiment changes the first radial window's odd-axial poloidal mean
coefficients at the `k=11` and `k=19` scale knots. Halving those two
coefficients opens the sampled cone point `y=0.325, eta=+0.05` at both
scales, but the uncorrected change destroys the outer moment conditions:
the largest absolute sampled moment rises from `2.64e-6` to `243`.

The repair holds those two coefficients fixed and uses the other 34
coefficients to restore the twelve sampled physical outer moments at
`k=11,15,19`. SLSQP converges in 19 iterations. The largest normalized
moment defect is `2.41e-15`; the largest absolute moment is `1.89e-11`.
At `y=0.325`, the sampled cone passes at `eta=-0.05,-0.0125,+0.025,+0.05`
for both `k=11` and `k=19`. The `eta=-0.075` point still fails. The
`eta=+0.05` cone ratios are `0.992` and `0.972`, so the positive-side
margin is narrow at the earlier scale.

This is a transferable *mean-geometry adjustment*, not a scale recursion.
The test uses five axial samples at two scales; it does not establish a
connected cone on a continuum of space, time, and scales. It also does not
reduce the full wave-plus-mean momentum residual. A next useful step is to
optimize the mean with a wider cone margin while maintaining the moments,
then directly evaluate a spatially supported exact-curl correction and
full residual on adjacent scales and interior pulse times. Recursive
progress requires residual contraction with each refinement, not merely
matching normalized shapes or sampled moments.

Reproduce with `python experiments/root_st073/midplane_axial_cone_knob.py`
and `python experiments/root_st073/midplane_axial_cone_moment_repair.py`.
The corresponding JSON files contain the coefficients, optimizer result,
and all sampled cone rows.

## Integer-scale holdout and middle-knot repair

Checking the first, two-coefficient repair on every integer scale from
`k=11` through `k=19` reveals that its `eta=+0.05` cone fails at
`k=12,13,14,15,16`. Its cone ratio rises from `0.992` at `k=11` to
`1.030` at `k=16` and passes again at `k=17`.

Halving the same first-window odd-poloidal coefficient at the middle
`k=15` knot as well, then reoptimizing the other 33 coefficients, restores
the twelve sampled moments with normalized maximum `1.50e-14`. The
`eta=+0.05` cone passes at each of the nine integer scales `k=11..19`.
Its weakest sampled margin is only `1 - 0.998181 = 0.001819` at `k=12`.
The maximum momentum norm over four local axial samples is `0.973–0.981`
times the original mean at most intermediate scales, but its absolute
value increases from `4.49e5` at `k=11` to `1.64e9` at `k=19`, a factor
of `3655`. This is a mean-only local screen; it does not show recursive
residual contraction. Continuous scale, radial, angular, and time support
also remain unverified.

Run the expanded experiment with
`python experiments/root_st073/midplane_axial_cone_moment_repair.py --all-knots`
and `python experiments/root_st073/midplane_axial_cone_scale_scan.py --all-knots`.
The default versions of those commands reproduce the earlier two-knot
candidate and its failing integer-scale holdout.
