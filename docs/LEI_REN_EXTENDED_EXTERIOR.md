# Extended exterior and inward moment targets

This route replaces the incompatible short join as the next construction
candidate. The existing callable short-join field remains a geometry seed.
Source: Lei and Ren, [arXiv:2609.35406v1](https://arxiv.org/html/2609.35406v1),
terminal conditions (1.2), moments (2.21)--(2.22), and the inward heat collar.

## Radius selection

`lei_ren_part1_extended_radius_screen.py/.json` evaluates the saved finite
core on 33 points in `Z in [-1,1]`. It compares the renormalized angular
target with `F_core(.005,Z) Rb^2` and requires the heat endpoint to lie below
the core-exit value. Radii 64 and 256 fail; 1024 and 2048 pass these two
necessary conditions at all sampled points. Neither pressure compatibility
nor a sufficient monotone connection follows from this screen.

`lei_ren_part1_join_thresholds.py/.json` locates the necessary angular
threshold at 999.6394381 and endpoint threshold at 499.6662749, both worst
at Z=1 on the sampled grid.

The prototype uses `Rb=2048`, `ell=.5`, `h=.001`, `c_inf=.1`, and
`epsilon=delta/100=2e-5`. Its collar inner boundary is
`Ra=Rb exp(-ell)=1242.1747910914733`. The earlier ell=.75 choice would put
Ra=967.4067 below the necessary angular threshold; its actual inward target
also fails the core-exit monotone angular bound. A screen at Rb alone is
insufficient. With ell=.5, actual inward targets pass the angular and
endpoint necessary conditions at all 33 Z points, with minimum margins
42.40847573 and 8.41183917e-5 respectively. This remains a conditional
screen for the old core, which will change under recomputed pressure.

## Transfer the actual exterior data

`lei_ren_part1_exterior_targets.py` accepts the actual pure-swirl exterior
`F(R,Z)` between Ra and Rb and the exact heat exterior from Rb to infinity.
It transports terminal conditions inward with logarithmic radial quadrature:

- Angular target: heat renormalized target at Rb minus `integral_Ra^Rb 2 R F`.
- Axial and mixed targets: zero, since the entire supplied exterior has Uz=0.
- Quadratic target: minus the complete exterior integral of `-R F^2`.
- Pressure at Ra: heat pressure at Rb minus `integral_Ra^Rb F^2`.

The cumulative pressure moment is `P(Ra,Z)-P0(Z)`. It remains null unless
the same candidate's recomputed axis pressure P0 is supplied. The pressure
datum of the old short-join core must not be silently reused.

`lei_ren_part1_exterior_targets_checks.py/.json` evaluates the source collar
factor `1-epsilon exp(-4/y^2)`, `y=log(Rb/R)`, over 33 Z points. Independent
pure-heat transfer identities have maximum absolute discrepancy 8.53e-14;
32/64-node collar differences are below 2.85e-14. Finite-difference radial
derivatives agree with the local moment densities to relative 1.01e-8.

The JSON records the actual angular, quadratic and pressure targets at Ra.
These are dimensionless similarity-profile moments. The positive quadratic
target requires interior axial energy to balance the negative swirl term;
repairing only the zero axial and mixed moments is insufficient.

## Actual angular-matched candidate

The actual `HeatCollar` supplies pressure, shear and ODE-integrated stress,
including flat-factor-normalized values near the zero-stress heat boundary.
Its checks record pressure-identity relative error 1.33e-10, angular stress
ODE residual below 7.72e-15 and axial stress ODE residual below 4.87e-24.
Sampled normalized collar cone margins are positive; minimum
kappa-(2+delta/2) is 9.9679e-4. Edge limits are checked at y<=5e-5 and retain
their measured first-order approach. These are finite local collar checks.

`lei_ren_part1_extended_swirl.py` now constructs a real extended F from the
saved regular core, rather than merely screening radii. It preserves the
core through R=.005. A local ODE transports negative radial shear smoothly
to a weak power-law anchor by R=.01. A flat logarithmic blend connects that
anchor to the actual collar; its timing is solved from the inward angular
target separately at each Z. The blend starts where anchor and extended
heat target cross, so the difference term preserves negative radial shear.

`lei_ren_part1_extended_swirl_checks.py/.json` checks 33 Z points, using
256-node integration independent of the 128-node timing solve. Maximum
angular defect is 2.01e-11. Sampled F is positive and F_R negative in every
layer; independent radial derivatives agree to relative 5.48e-7.

The candidate supplies its actual full swirl-pressure integral. Its new
axis pressure differs materially from the old core datum: at Z=0 about
-0.09887287 versus -0.06934343, and at Z=1 about -2.31260e-5 versus
-0.03684420. These are a required core-rebuild handoff, not pressure closure.
The weak anchor does not itself have the final admissible stress cone.
Negative F_R alone must not be promoted to kappa>2 or a whole-background
stress claim. The actual angular profile is an intermediate connection
candidate whose pressure, axial moments and shear still require restoration.

`lei_ren_part1_axial_quadratic_repair.py` supplies a finite three-moment axial
repair: two linear constraints followed by the quadratic condition in a bump
nullspace. Synthetic nonconstant-swirl/nonzero-base cases use a nonzero
quadratic target, with independent 32-node residual below 4.20e-11. Those
checks establish the repair primitive, not an assembled background.

## Common-pressure rebuilt core

`lei_ren_part1_extended_pressure_core.py/.json` now rebuilds the finite
nonlinear core under the actual extended swirl's pressure integral and
resolves the angular blend after each rebuild. The 257-node iteration
converges in three steps with collocation pressure defect 4.32e-10.
Independent 16-point pressure holdout with refined radial quadrature has
maximum defect 7.05e-9. `load_extended_profile()` replays this candidate and
by default refuses a saved receipt that failed the pressure holdout.

The earlier 65-node attempt also converged at nodes but failed the independent
holdout at 0.00189012; it is retained as `_coarse65.json`. This is why node
convergence alone is insufficient. The rebuilt 257-node candidate retains
angular closure over 33 Z points at 2.01e-11, positive swirl and sampled
negative radial shear, with independent F_R error below 5.57e-7 in
`lei_ren_part1_extended_swirl_rebuilt_checks.json`.

## Next construction work

1. Use the actual collar at Ra as the outer jet and five-moment boundary data.
2. Retain the actual angular-matched, common-pressure candidate as the input
   to axial/quadratic repair; rebuild its pressure whenever swirl is modified.
3. Preserve the regular core and actual collar jets during repairs. Negative
   F_R of the intermediate weak anchor is not a final cone certificate.
4. Restore axial, mixed and quadratic moments on the same profile. Verify
   between interpolation nodes, including the axial endpoint behavior.
5. Compute actual profile-derived stress and its cone margins over the whole
   connection. Then implement higher-order coefficients and measure the
   physical background residual, stress divergence and remainder across scales.

Passing local collar or moment-transfer checks does not establish a global
admissible background, all-orders flatness, or NS scale-recursion closure.
