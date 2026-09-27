# Matching and shape-direction constrained tangent

The fixed-wave tangent solve now includes the three local geometric
directions requested by the project, in addition to four integral moments
and all 81 sampled cone inequalities. This is an instantaneous construction;
it does not establish an evolving NS field or scale recursion.

## Shape operator and solve

`wave_shape_tangent_rows.py` forms three affine rows on the 180 time/pressure
controls. The instantaneous velocity and vorticity are frozen on the same
fixed physical cylinder used by `wave_tangent_observables.py`. The rows
differentiate enstrophy-weighted radial RMS, axial/radial RMS ratio and
weighted angular speed. Pressure columns are exactly zero. Velocity-curl
responses use analytic Fourier jets. The baseline is inferred from the
earlier central time difference minus the selected analytic response; its
exact reconstruction is by construction, not independent validation.

`shape_constrained_wave_tangent.py` eliminates the four moment equalities
and solves a convex least-squares problem over all 180 controls with cone
and shape inequalities. It requires fractional contraction and elongation
rates of at least 0.01 per k and a fractional weighted-spin magnitude
increase of at least 0.01 per k. These are explicit small direction margins,
not derived recursive scaling exponents or a substitute for the PDE target.

Assembled results:

| Quantity | Value |
|---|---:|
| Training momentum volume L2 | 2,868,845.21 |
| Earlier cone-compatible training L2 | 2,868,781.90 |
| Integral moment maximum | 5.18e-9 |
| Cone passing locations | 27 / 27 |
| Minimum cone margin | 9.9999981e-5 |
| Radial RMS derivative per k | -1.1358189e-5 |
| Aspect derivative per k | +0.0020859563 |
| Weighted angular-speed derivative per k | +7.6915423e9 |

The training residual cost is about 0.0022%. Actual-field shape and
momentum/compatibility replays are separate reports; assembled signs must
not be reported as successful physical evolution before those results.

## Actual-field result and interval failure

The small-margin candidate's independent full-field replay gives momentum
L2 2,847,211.46 and maximum 2.1020837e11. Its order-96 moment maximum is
2.1012e-7 and all 27/27 cones pass (81/81 inequalities; minimum margin
9.9758682e-5).

All three central finite-difference shape directions pass on the actual
field. However, the forward interval from k0 to k0+1e-6 still expands the
radial RMS and decreases the aspect ratio. Curvature dominates the small
linear margins. This candidate is therefore not accepted even as a
successful forward geometric increment. `wave_tangent_observables.py`
now records one-sided rates and forward-interval flags in addition to the
central derivative, without rerunning expensive field evaluations when
three saved samples are already available.

A separate candidate `shape_direction_margin_tangent.json` increases the
explicit fractional direction floor to 10 per k. Its assembled rates are
[-0.01135819, +2.08595631, +7.7038227e9], while training L2 is 2,869,199.24
(about 0.0145% above the earlier cone-compatible candidate). All assembled
moments and cones still pass. Its actual interval and compatibility replays
are pending; the larger margin is an experimental response to observed
curvature, not a proven recursive exponent.

## Reproduction

```powershell
python experiments/root_st073/wave_shape_tangent_rows.py
python experiments/root_st073/shape_constrained_wave_tangent.py
python experiments/root_st073/wave_tangent_observables.py --candidate experiments/root_st073/shape_constrained_wave_tangent.json --output experiments/root_st073/shape_constrained_wave_observables.json
python experiments/root_st073/wave_dynamics_replay.py --source experiments/root_st073/shape_constrained_wave_tangent.json --output experiments/root_st073/shape_constrained_wave_replay.json
python experiments/root_st073/wave_dynamics_mean_compatibility.py --candidate experiments/root_st073/shape_constrained_wave_tangent.json --replay experiments/root_st073/shape_constrained_wave_replay.json --output experiments/root_st073/shape_constrained_wave_compatibility.json
```

Full momentum max/L2 below 1e-3, finite-time persistence, global finite
energy, complete inner/transition/exterior matching and recursive geometry
remain required. Sampled shape-direction inequalities alone prove none of
these. The geometry is a fixed-cylinder vorticity statistic, not an identified
global core boundary.
