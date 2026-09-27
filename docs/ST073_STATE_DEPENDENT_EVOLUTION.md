# State-dependent outer evolution diagnostic

The new `experiments/root_st073/outer_feedback_evolution.py` evolves nine
compact azimuthal velocity coefficients instead of holding their time slopes
fixed. At each midpoint stage it rebuilds the velocity, advection, shear,
pressure response and stress constraints at the actual scale and state.
Nine pressure coefficients and nine coefficient slopes are solved together.

The local control problem enforces four physical moment equations and cone
inequalities at the 22 spatial training nodes. A fixed-reference strictly
convex quadratic objective selects among feasible controls. This objective
keeps controls near the spatial seed; it does not minimize or certify full
momentum. The script stops and records failure if shear positivity or the
constrained control solve fails.

With remaining time tau and k = -log2(2 tau), physical time differentiation
is minus tau differentiation. Thus a coefficient slope a'(k) contributes
U a'(k)/(tau log(2)) to momentum. The swirl value and slope bases are kept
separate so that each new state contributes its changed nonlinear advection.
The compact axisymmetric swirl addition has zero divergence analytically.
This property alone does not establish the full field's energy bounds.

The first diagnostic uses explicit midpoint steps of 0.0005 over
k in [11, 11.001]. A cubic Hermite interpolant of the state, using solved
endpoint slopes, and a cubic interpolant of pressure give one callable
field. Interpolation is forbidden outside the saved interval. The field
is replayed at k = 11.00035, 11.0005, 11.00065; the finite-difference time
stencils stay inside the interval. These replays include independently
integrated moment identities, outer cones and full-vector sampled momentum.
Neither the interpolant nor a finite number of samples certifies an ODE
solution or a continuous interval bound.

The moment locations and six cone replay locations coincide spatially with
training locations; independence here is in time, quadrature and the moment
identity, not a new spatial support certificate. The momentum grid has
additional spatial holdouts. Replaying the eight spatial neighborhood
holdouts over time remains outstanding.

`outer_feedback_fixed_comparison.py` uses exactly the feedback run's initial
state and control but keeps its slopes fixed. It evaluates the same
order-96 integrated moments at k=11.0005, providing a matched comparator.

## Completed first trajectory

All five state solves completed without loss of constrained feasibility.
The reconstructed trajectory produced these independent time replays:

| k | Integrated moment max | Outer cone samples | Sampled full momentum max |
| --- | ---: | ---: | ---: |
| 11.00035 | 0.00098697034 | 6/6 | 1422537 |
| 11.00050 | 0.00049910112 | 6/6 | 1422759.2 |
| 11.00065 | 0.00095690791 | 6/6 | 1422981.5 |

At the shared midpoint k=11.0005, the matched fixed-slope comparator has
moment max 1.6580147; feedback has 0.00049910112, a
factor 3322.0 smaller. This is evidence that state-dependent controls
suppress the observed short-time moment drift. The worst replay is close
to 1e-3 and has no established numerical-error margin. These are moment
values, not the full momentum residual required by the goal.

The 18/18 cone sample passes and the three moment samples do not establish
continuous time/support validity. Full sampled momentum is still about
1.42e6, spatial-volume L2 and finite energy are not certified, and no
recursive contraction or nonaxisymmetric correction is established.

## Integration cost

`grouped_outer_cache.py` shares radial quadrature among equal physical-z
centers, preserving every target radius as a panel boundary. Its 22-node,
order-24 benchmark uses 2,062 points instead of 4,870 and takes 46.80 seconds
instead of 118.62 seconds (2.53 times faster in this run). Independent panel
layouts change quadrature rounding/truncation: the largest cone-output
difference is 1.8e-4, approximately 1.3e-6 relative in the affected ratio.
The benchmark is saved in `grouped_outer_cache_check.json`.
The cache accepts explicit additional radial cuts. The feedback caller
passes P_BREAKS, including the remote support at .93 to .99. These extra
cuts are beyond all current cone targets (y <= .755) and do not change
their quadrature; a direct node/weight equality check confirms this.

## Relation to the target

The [OpenAI paper, Section 7, Proposition 7.2 and equation (7.13)](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
constructs transverse amplitudes and pressure along pulse paths. Sections
8–9 couple mean corrections with residual improvement. The present
axisymmetric feedback is only a prerequisite experiment for a background
that remains compatible over time. It is not that pulse inverse, a
nonaxisymmetric realization, or a recursive contraction.

The authoritative run artifact is `outer_feedback_evolution.json`. A partial
artifact contains only completed stages; missing replay results mean the
run has not established a replay outcome. All reports retain
`accepted=false` and `scale_recursion_established=false`. The full goal
still requires full momentum max and spatial-volume L2 below 1e-3, finite
energy, prescribed forcing/domain, recursive improvement and deliverables.

Run from the repository root:

```powershell
python experiments/root_st073/outer_feedback_evolution.py
python experiments/root_st073/outer_feedback_fixed_comparison.py
```
