# Whole original inactive-gap physical tensor — 2026-10-04

The same complete inactive gap xi in [11,13-4mu], d=13-xi in [4mu,2],
Z in [-1,1] now has a completed symmetric physical stress tensor,
cylindrical stress mixed3, diagonal/divergence mixed2 and all three
physical remainder components mixed2. Cartesian tensor, divergence,
remainder and the decomposition -div(T)+E are exported as signed sectors.
This consumes the current five-moment/absolute-pressure companion.

## Full physical operator and units

Use the original implicit lambda relation and fixed positive viscosity nu:

`lambda^2-lambda^(2delta)*z^2/nu=tau=1-t`,
`R=r^2/(2nu*lambda^2)`,
`Z=z/(sqrt(nu)*lambda^(1-delta))`.

The original units remain
`ur=sqrt(nu)*lambda^(-1)*Ur`,
`utheta,uz=sqrt(nu)*lambda^(-1-delta)*(Utheta,Uz)`,
`p=nu*lambda^(-2-2delta)*P` and
`T=nu*lambda^(-2-delta)*Tprofile`.

The completed tensor has
`Trtheta=Ttheta`, `Trz=Tz`,
`Ttheta_theta=r*partial_z(Tz)`, symmetric partners, and other entries zero.
Its radial divergence is structurally zero. The same full meridional
momentum operator is consumed from the accepted pulse-end source, including
pressure/centrifugal cancellation, radial transport and radial viscosity;
its 28 identities are not rerun as a new upstream milestone.

## Nonzero radial history and the correct factor

The actual gap input has Uz=0, but nonzero Mz and radial velocity remain.
The radial history must use
`D1=D0*exp((.5-mu)*d/mu)`, with its exact reduced logarithm.
Using terminal D0 with the unchanged gap coefficient would lose the
backward history factor. Source velocity rows already contain all ordinary
derivatives; the frozen basepoint log is not differentiated again.

Direct replay of the current row algorithm proves Ur ordinary derivatives
are `(-.5)^j*Ur` in the gap, so its radial vector Laplacian is zero.
Conservative interval rows for this term may still enclose zero with width;
the symbolic identity and signed enclosures are recorded separately.
The actual physical remainder reduces to

`Er=partial_t(ur)+ur*partial_r(ur)-nu*partial_zz(ur)`,
`Etheta=-nu*partial_zz(utheta)`, `Ez=0`.

Time change, radial nonlinear transport and axial viscosity remain.
In particular, Er is not set to zero at the right gap/end interface.

## Same-source end join and focused evidence

The formal coupled right endpoint d=4mu maps exactly to end s=-4.
Actual original end R/B/H log recipes and the production physical radius route for pulse_gap/pulse_gap_end/pulse_end are replayed. The original full beta
moment histories acquire exp(4*lambda_i), matching gap D1/D2; both linear
and quadratic radial factors agree. Current functional
velocity4/moment4/stress3/pressure4 source joins compose through the same
nu, lambda, time, radius and pressure units to physical stress3,
diagonal/divergence2 and remainder2. This is a source-functional join,
not overlap of interval enclosures. Native functional moment/energy/pressure/angular/support and axial-order identities plus the canonical P0/FTC getter are consumed directly. Both actual pullback calls preserve log_tau/theta/viscosity arguments; the interface statement applies at identical finite log_tau, angle and fixed nu>0.

Focused checker PASS: 79 source/operator/interface identities,
427 current hashes, 792 finite signed physical rows and 24
structural zero rows. The current whole-gap report is independently recomputed.

An independent moderate fixture integrates both full original beta histories,
retains nonzero Mz/radial velocity and pressure memory, and differentiates
Cartesian time, convection, pressure, full Laplacian and completed tensor
divergence directly. At nu=.01/.7 it supplies 114 comparisons
through stress3 and diagonal/divergence/remainder2, with nonzero radial and
angular errors and zero axial error. Tolerance 1.0e-55; maximum positive
enclosure miss 3.14298702824123051704566361771e-79. These fixture parameters are not actual source
values and do not establish a cone, global flatness, or corrected NS accuracy.

Run only the bounded new stage:

`python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage pulsegapphysical`

The continuous whole-gap regional cone and end/tail composition are now
admitted in docs/PULSE_GAP_CONE_2026_10_04.md. Main pulse, exit, entrance, upstream finite-width bridge
feedback and completed global tensor admissibility remain. Global temporal
flatness/physical-volume/required-domain energy, true n-dependent recursion,
oscillatory stress correction and corrected Cartesian residual/dynamics
remain unfinished.
