# Joint moment and cone correction

The fixed growing wave now has a mode-0 tangent correction satisfying all
four assembled integral moments and all 81 cone inequalities at 27 sampled
locations. Independent full-field momentum and moment/cone replays are
complete. No recursive step is accepted.

## Construction

`wave_mean_cone_projection.py` exports affine time/pressure control rows
and quadratic wave-force forms for the 27 saved locations. The wave uses
the full curl correction, and the cone calculation retains its mean force.
`wave_moment_cone_tangent.py` holds the wave and nonzero-harmonic tangent
controls fixed, eliminates four moment equalities, whitens the remaining
mode-0 least-squares objective, and solves a convex inequality-constrained
problem. A phase-I linear program establishes feasibility before SLSQP.
Only responsive cone rows receive an additional 1e-4 safety margin.

The assembled moment maximum is 6.8103e-9, the minimum cone margin is
9.9999994e-5, and all 27 locations / 81 inequalities pass. These numbers
alone are not independent physical-field evidence.

## Independent momentum result

The actual Cartesian field is evaluated at 18,720 points on the existing
order-13 spatial replay grid, unused by this fit, at one reference time.
Zero external forcing is used in this diagnostic.

| Candidate | Physical-volume momentum L2 | Momentum maximum |
|---|---:|---:|
| Moment-only correction | 2,661,414.30 | 2.0373372e11 |
| Joint moment/cone correction | 2,847,158.20 | 2.1019005e11 |

The L2 cost of restoring the sampled cones is about 6.98%. Relative to
the older dense tangent baseline (7,376,639.77), the L2 reduction remains
about 61.4%. Neither percentage is overall project completion.
The mean harmonic accounts for 68.72% of squared residual and mode 2
for 30.67%. These remain the dominant correction targets.

The saved higher-grid energy operators give growth lambda +997.0958.
They remain applicable because instantaneous velocity and geometry are
unchanged. This is finite-grid positive growth, not a verified 1000 floor
or an evolving shrinking vortex.

## Independent compatibility

The frozen separate actual-field replay gives an order-96 joint moment
maximum of 2.1218875e-7. All 27 locations and 81 inequalities pass the
order-64 cone replay, with minimum margin 9.9758657e-5 and 8,795 radial
quadrature points. The moment error is larger than the assembled value,
but remains small at this quadrature. Neither this sampled cone pass nor
the integral moment result is an assertion of pointwise PDE accuracy.

## Reproduction

After the existing wave momentum cache and moment candidate are generated:

```powershell
python experiments/root_st073/wave_mean_cone_projection.py
python experiments/root_st073/wave_moment_cone_tangent.py
python experiments/root_st073/wave_dynamics_replay.py --source experiments/root_st073/wave_moment_cone_tangent.json --output experiments/root_st073/wave_moment_cone_replay.json
python experiments/root_st073/wave_dynamics_mean_compatibility.py --candidate experiments/root_st073/wave_moment_cone_tangent.json --replay experiments/root_st073/wave_moment_cone_replay.json --output experiments/root_st073/wave_moment_cone_compatibility.json
python experiments/root_st073/cached_wave_growth.py experiments/root_st073/wave_moment_cone_tangent.json experiments/root_st073/wave_moment_cone_growth.json
```

Full momentum remains many orders of magnitude above 1e-3. Matching at
sampled locations and at a single time does not establish matching over
a continuous domain or time interval. The next construction varies wave
shape with the moment and cone constraints included in the inner solve;
it must retain nonlinear interaction residuals. Continuous time evolution,
scale recursion and simultaneous core contraction, relative elongation and
spin-up remain unestablished.

This follows the working interpretation of Sections 7–9 recorded in
`ST073_PAPER_ROUTE.md`: retain the complete curl-induced wave force and
nonlinear interactions, and revisit mean compatibility after changing a
wave. The sampled physical cone here is not the manuscript's normalized
leading-stress cone, and the finite QP is not its recursive construction.
Reference: [OpenAI-hosted manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).

## Bounded outer wave optimization

`wave_cone_codesign.py/.json` now implements a varying-wave outer solve
with all moment/cone constraints retained in the inner tangent fit. The
column-normalized harmonic elimination retains rank 144/144. The envelope
gradient agrees with independent inner re-solves to relative error 3.01e-7.
The bounded outer SLSQP run reaches its iteration limit without a feasible
improving step, so the selected normalized objective remains 1.0. This
adds a working constrained objective, not a new improved wave candidate.
Shape-rate constraints are handled separately in the fixed-wave construction.
