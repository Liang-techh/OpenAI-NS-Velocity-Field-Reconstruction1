# Radial quadrature invalidates the earlier small-moment candidate

The candidate in `midplane_connected_cone_edge_repair.json` was fitted to
12-point Gauss quadrature on just two radial panels: the inner core and
the whole transition layer. That second panel does not adequately resolve
the separated compact bumps. Re-evaluating the identical field yields:

| Radial rule | k=11 max moment | k=15 max moment | k=19 max moment |
| --- | ---: | ---: | ---: |
| Original two panels, 12 points | `4.72e-8` | `2.83e-7` | `7.69e-6` |
| Original two panels, 48 points | `223.98` | `3543.39` | `55935.51` |
| Original two panels, 96 points | `212.09` | `3353.96` | `52951.33` |
| Split at bump edges/midpoints, 24 points/panel | `213.181790` | `3372.199231` | `53239.824212` |
| Same split, 48 points/panel | `213.181840` | `3372.200029` | `53239.835131` |

Thus the earlier small cached and directly replayed moments were both
artifacts of the same under-resolved quadrature. Direct replay alone did
not provide an independent numerical integration check. The candidate's
moment-closure claim is withdrawn; earlier reports using only the same
12-point rule also require a refinement check before acceptance.

The intermediate-scale screen still finds all four tested support corners
passing at `k=12,13,14,16,17,18`, but its original coarse moment values are
not reliable physical integrals. Independently increasing radial stress
quadrature from 12 to 96 points leaves the tested near-edge cone ratios
below one: `0.950 -> 0.9131` at k=11 and `0.950 -> 0.8938` at k=19.
This distinguishes the moment integration failure from the sampled
geometric improvement.

The moment cache now accepts explicit bridge-fraction panel breaks.
The corrected fitting route splits at every separated bump edge and
midpoint, fits with 24 points/panel, and replays with 48 points/panel.
It also builds the affine correction jets directly on a zero background,
avoiding repeated base evaluations and subtractive cancellation. Stress
constraints use 48-point integrals with a 96-point direct replay. These
are improved numerical rules, not certified error bounds; candidate
acceptance still depends on the independent replay and convergence.

Reproduce the audits:

```text
python experiments/root_st073/midplane_connected_moment_quadrature.py
python experiments/root_st073/midplane_connected_moment_quadrature.py --split
python experiments/root_st073/midplane_connected_cone_quadrature.py
python experiments/root_st073/midplane_connected_cone_interscale.py
```

The corrected constrained solve is
`python experiments/root_st073/midplane_connected_cone_repair.py --resolved`.
Its separate output must be assessed before replacing any mean field.

## Resolved joint solve outcome

The corrected SLSQP run reached its 300-iteration limit, rather than
establishing feasibility or infeasibility. Its normalized moment maximum
is 42.919704 (independent-order replay 42.919720); the direct absolute
moment maximum is 52266.8043. Each of k=11,15,19 passes only 8 of 9
sampled cone nodes, with maximum ratios 3.53104, 3.86884, 4.30800.
The candidate is rejected and does not establish scale recursion.

The next constructive diagnostic separates moment-only feasibility from
joint cone optimization, checks the resolved moment Jacobian and scaling,
and independently replays any resulting candidate. A failed bounded
optimizer is not evidence that the profile family is mathematically
infeasible. No field has been promoted on the basis of this run.

## Moment-only diagnosis

`midplane_resolved_feasibility.py` integrates the exact quadratic moment map
of the affine velocity jets and uses its analytic Jacobian. The initial
12-by-36 normalized Jacobian has row rank 12 and condition number 1.63068;
there is no initial linear rank deficiency. A tensor-layout error in the
first run was fixed before the reported numerical run; the corrected
polynomial was checked against direct cached moments and its Jacobian
against a directional finite difference.

The moment-only solve exhausted 180 evaluations. Its normalized maximum
falls only from 43.80267 to 41.81849. The remaining error is predominantly
axial: direct 96-point split replay gives maxima approximately 204.0601,
3215.4819, 50689.3895 at k=11,15,19. The corresponding 48-point replay
agrees in normalized maximum to about 1e-8, but the largest absolute
cross-order difference is 0.001753; this is not an absolute 1e-3 certificate.
Only 4/9, 6/9, 6/9 sampled cone nodes pass for this candidate.

The final Jacobian's smallest singular value falls to about 0.000580
(from 1.049 initially), suggesting the optimizer approaches a weak axial
response direction. It does not establish an unattainability theorem.
Next analyze the signed axial-moment quadratic forms and their extrema,
then decide whether the poloidal basis or axial/time profile needs to
change. Repeating the same joint cone optimization is not the next step.
The candidate remains rejected; scale recursion is unestablished.

Reproduce: `python experiments/root_st073/midplane_resolved_feasibility.py`.
