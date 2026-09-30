# Source outer construction and first correction inputs

The source Section 6 candidate is a separate construction from the retained
weak-anchor physical field. It must not inherit that field's pressure/core
receipt or be accepted by changing candidate metadata.

`lei_ren_part1_paper_outer.py` implements the reference velocity
Utheta=Pstar/(1+Z^2)*(R/Rref)^(.1), Uz=4Z, then the source radial slope
schedule, axial turnoff, Z flattening and terminal heat connection.
Its reference F is singular at the axis; this temporary reference is not
the regular finite core. It is not yet a global callable physical solution.

The derived checkpoint includes log(Rv/Rp)=13/mu. Ordinary floating point
cannot preserve short transitions after that checkpoint for small mu.
Evaluate using stage-local logarithmic coordinates and Decimal checkpoint
positions. Physical radii must not be blindly exponentiated. The paper's
absolute constants are existential; a numerical input satisfying the
explicit scalar inequalities does not certify their required smallness.
No compressed numerical schedule may inherit the source cone theorem.

`lei_ren_part1_paper_outer_pressure.py` integrates the SAME candidate's
backward pressure in normalized P/Pstar^2 units. Beyond a pre-heat cutoff
y>=1, the slope schedule and monotone flattening give the omitted-tail bound
Utheta(y)^2/[2 Pstar^2 (1-epsilon)^2]. The finite-quadrature error remains
separate. Its reference-extension axis pressure is temporary: Section 7
pressure restoration and a regular shared-pressure core are still required.

`lei_ren_part1_paper_outer_closure.py` implements the exact two angular bump
equations (7.21), the axial affine system (7.31), and its quadratic energy
root (7.34). It also constructs the normalized source bump weights and
translated gamma matrix. These kernels require actual candidate moment
defects and pulse integrals. Synthetic algebra checks do not close the
source outer candidate. For very small mu, a divided-exponential row scaling
keeps the axial system condition near 4.32 instead of losing its determinant
in ordinary floating point. Tiny-mu right-hand sides must retain Decimal
precision; stage-local bump offsets replace collapsed absolute centers.

`lei_ren_part1_paper_waiting_length.py` integrates the actual normalized
incoming angular-moment ODE and solves the pre-heat scalar matching equation.
It retains the tiny heat deficit separately before exponential amplification.
For the declared demonstration input, waiting length is about 82.448458;
independent backward source-stage replay differs from the incoming ratio by
2.68e-10. `apply_waiting_root` builds the matched schedule. The exact heat
moment still requires angular bump correction, so this is pre-heat matching
only, not an all-five-moment or regular-core certificate.

`lei_ren_part1_paper_heat_defects.py` now evaluates the positive heat
replacement inputs r_H and s_H from (7.13),(7.17) in logarithmic form.
It retains the first heat Taylor term with an explicit analytic relative
remainder bound; it does not subtract two rounded copies of H=1. The
infinite angular tail is integrated analytically, including its 1/h factor.
The short collar is integrated from the actual sigma/flat cutoff. Independent
Gamma-expectation quadrature checks the Taylor inequality at representable
arguments. For the matched demonstration, both log defects are about
-2.7191573448e28: positive but impossible to materialize in binary64.
Pressure-collar quadrature refinement differs by 4.24e-11 in its normalized
integral. That quadrature uncertainty is distinct from the analytic Taylor
bound. The incoming Z-dependent r_pre is now evaluated by the actual
flattening difference ODE, then propagated by mu^(30(1-mu)) in logarithmic
form. Its initial difference is exactly zero because the Z amplitude
separates before flattening. Early common-moment contraction omissions are
inherited explicitly from the incoming receipt. The waiting-root numerical
error remains separate. Angular coefficients and the corrected profile still
need to be bound; these logs do not certify angular closure.

At Z=.5 the propagated preheat defect has log(r_pre) about -1925.4724;
the finite flattening difference is about -2.47e-17. It is integrated using
the positive variation-of-constants formula, so an absolute solver tolerance
cannot erase the difference. Combined ODE/quadrature refinement changes the
finite difference by about 2.66e-7 relatively. This is numerical evidence,
not an interval certificate for the coefficient functions or their derivatives.

## Angular coefficients and profile representation

`lei_ren_part1_paper_angular_correction.py/.json` now solves the small
quadratic branch for the candidate-derived numerical inputs using mpmath's
arbitrary exponent representation. It divides d_j by r before solving,
so tiny coefficients never become zero. The angular multiplicative bump
is callable in stage-local coordinates; `corrected_at_log_radius` retains
the base log amplitude and its tiny log1p correction as separate terms.
Both bump multipliers pass the positivity condition for the declared input.
Independent actual bump quadrature gives angular increment/r error 5.11e-15.

This is not an exact pressure-cancellation certificate. The weights retain
finite quadrature uncertainty, the heat inputs are bounded Taylor values,
and a finite precision subtraction cannot resolve the pressure target/r.
The receipt reports this unresolved scale explicitly. The numerical waiting
root's error also exceeds the tiny physical defect scales and must remain
separate; do not declare all moments closed from algebraic residuals.

## Actual axial pulse inputs still required

Section 7.5 fixes gp(xi)=[1-sigma(xi-10)]*integral_0^xi sigma(50v)dv
for xi>=0 and zero otherwise. Its support is [0,11]. The two end bumps
are centered at stage-local t=13/mu-3 and 13/mu-1. Kp is the actual
integral_0^13 exp(-2xi)*gp(xi)^2 dxi, with source bound .24<Kp<.246.
The stable affine API requires EACH source RHS row to be divided by
exp(13*lambda_i/mu), lambda_i=.5-i*mu. Do not feed unnormalized pulse
integrals into its row-scaled matrix. The energy target (7.34) includes
the actual future integral of the corrected angular profile from Rv to
infinity, so it cannot be substituted before the angular repair is bound.

`lei_ren_part1_paper_axial_pulse.py` now materializes this fixed pulse,
its derivative and actual Kp integral. Its checks give Kp=.24504962020069448
and quadrature refinement 2.03e-15, within the source range. Actual linear
RHS integrals and the corrected angular energy target are now materialized
by the following modules, retaining their uncertainties.

`lei_ren_part1_paper_axial_incoming.py/.json` integrates the SAME temporary
reference and source outer candidate through Rp, producing m1, m2 and
E_prior with signed logarithmic row normalizations. Actual axial/mixed
moment refinement is about 3.14e-7/3.78e-7 relatively at Z=.5; prior-energy
refinement is 1.37e-12. This reference remains singular in F at the axis
and is not the future regular core.

`lei_ren_part1_paper_axial_pulse_moments.py/.json` integrates both actual
row-normalized pulse terms around the flat-endpoint saddle. The radial
stage remains logarithmic; omitted positive pieces have separate bounds.
Source-input quadrature refinement is 3.38e-9 relatively. Independent
original-coordinate quadrature at mu=1e-6 differs by 1.21e-10 relatively.

`lei_ren_part1_paper_axial_energy_tail.py/.json` integrates the actual
future corrected angular energy through flattening, angular bumps,
drop/steep/restore/waiting and heat exterior. It retains heat-collar and
exterior deficit bounds separately instead of rounding their widths to zero.
Q at Z=0,.5 is about 533.0691,806.2893; maximum quadrature refinement is
5.06e-13 relatively. The Z=0 energy-target contribution is 6.51e-37.

`lei_ren_part1_paper_axial_correction.py/.json` binds these actual inputs
to the source affine and energy equations, with an analytic determinant
and arbitrary-exponent end coefficients. The selected a_p is approximately
1.0100502663, inside (.9,1.2). Independent translated-bump quadrature
replays both normalized linear rows to 2.45e-19 relatively. Algebraic
working-precision errors are reported separately and do not replace the
larger incoming/pulse/tail integration uncertainties. Pulse and end-bump
axial factors are callable in separate stage-local coordinates.

## Unified source candidate and axial primitive

`lei_ren_part1_paper_axial_primitive.py/.json` supplies the actual normalized
pulse primitive in startup, bulk, cutoff and post-pulse coordinates. Startup
and cutoff tail bounds are retained. The bulk particular solution satisfies
the source transport ODE; the post-pulse identity reuses the same actual
normalized pulse rows. Moderate-mu independent quadrature checks test the
boundary and centered branches. Very small startup values may remain below
the startup evaluator's tolerance; the unified profile rejects a purported
zero at a positive interior pulse coordinate.

`lei_ren_part1_paper_corrected_profile.py/.json` now assembles source angular
and axial coefficients in one candidate. It integrates the SAME axial
velocity for Mz/R and derives V/R from source (3.9), with a Z derivative
of that primitive. End-bump moments are independently reintegrated rather
than forcing the exterior mass to zero. The numerical exterior mass and
its derivative remain part of the open all-five-moment/energy acceptance.
In the sampled bulk pulse, the independent M_R=Uz check differs by 2.09e-15
relatively and the radial Z stencil changes by 1.48e-14 under refinement.

The callable `velocity_from_tau(x,y,z,tau)` returns Cartesian u,v,w as
signed-log dictionaries with arbitrary-exponent values. `cartesian_chart`
also gives the corresponding physical coordinates. The physical API applies
the existing smooth axial cutoff to the meridional streamfunction, including
its radial correction; chart APIs remain unlocalized. Inputs may be decimal
strings. The physical-coordinate/chart roundtrip is checked separately.
The current coefficient stencil rejects charts too near |Z|=1. This is a
computable candidate API, not a regular-core or finite-global-energy claim.

`velocity(x,y,z,t,T=1)` provides the time-coordinate wrapper;
`velocity_values_from_tau` returns numeric mpmath scalars rather than the
signed-log receipts. Run from the repository root, for example:

```python
import sys
sys.path.insert(0, "experiments/root_st073")
from lei_ren_part1_paper_corrected_profile import CorrectedSourceProfile
profile = CorrectedSourceProfile()
u, v, w = profile.velocity_values_from_tau(".1", ".2", ".01", ".1")
```

`lei_ren_part1_paper_corrected_pressure.py/.json` uses the same angular
correction object for actual full/truncated bump pressure integrals, keeping
the backward pressure change separate from its much larger baseline.
The sign is P_corrected-P_baseline=-Erel^2 Delta_p_bump. Exact pressure
target cancellation remains unresolved, including zero-target endpoints
when finite-precision bump integration leaves a nonzero residue.

The pressure adapter now propagates the baseline P/Utheta^2 backward at any
logarithmic radius. Constant slopes use analytic contractions; finite
transitions use actual source-profile quadrature. The H=1 heat nominal value
and propagated heat deficit bound remain separate. Its Rref comparison with
the existing pressure integral differs by 9.0e-17. The unified profile's
`pressure_at_log_radius` and `pressure_from_tau` use this same adapter;
physical pressure scales by nu*B(z)^2*q^(-1-delta). Velocity and pressure
share the same implicit physical chart. The sampled chart/physical-coordinate
roundtrips pass independently of the primitive check. These are coordinate
and same-profile checks, not a full momentum-residual or stress-cone claim.

Next recompute SAME-profile stress, check coefficient derivatives and the
relaxed cone, and replace the temporary reference with the regular shared
pressure core. Exact pressure-target cancellation and all-five-moment global
closure remain unproved; the retained physical field is not replaced yet.

The binding order is: solve waiting length from actual normalized angular
data; restore angular/pressure moments; integrate actual axial pulse and
solve its linear/quadratic constraints; recompute actual pressure and stress;
construct and connect the regular inner core; restore inner five moments;
then apply the Section 11 shear loop under its verified relaxed-cone inputs.

## Actual first-order core sources

`lei_ren_part1_first_order_sources.py` already materializes the source
Omega0 and the preceding angular/axial viscosity terms from the actual
retained common-pressure core. It uses (13.8),(13.16)-(13.19), including
analytic radial jets of the finite core polynomial. Physical-coordinate
measurements independently confirm these sources at three time scales.
The radial pressure equation is P1_R=2 F0 F1-Omega0/(2R); the known
Omega0 source cannot be used alone as a completed pressure correction.

The measured leading angular and axial remainders are predominantly axial
viscosity. Removing those terms only in the diagnostic leaves a radial
source that requires the coupled correction. No velocity or pressure has
been corrected by this subtraction. The next implementation must solve
F1, Uz1, P1 jointly with zero positive-order axis data, extend them, restore
their moments and demonstrate actual remainder-order improvement. Use the
existing full radial recurrence where suitable, with the SAME new axis
data; do not reuse the unrelated ST073 axis traces as a source certificate.
# Actual initial-stage stress and coherent end-bump replay

The new `SourcePulseMoments(profile)` callable extends all five actual
cumulative moments through Rv. `AngularCumulative` propagates normalized
angular and swirl-energy moments analytically on long power stages;
`PaperPulseCumulative` uses the actual row-specific pulse primitives and
end-bump integrals; `PulseEnergyCumulative` integrates the quadratic axial
moment. The pressure moment uses the same full nominal pressure difference
from its temporary reference-axis value. Incoming, pulse, heat and
quadrature uncertainties remain inherited. The bulk quadratic radial
identity defect is 8.33e-15, and actual bulk stress is now recorded.
Neither terminal moments nor their numerical defects are forced to zero.

The current demo reference radius is incompatible with the paper's regular
core scale relation. `lei_ren_part1_paper_core_scales.py` records the
necessary test; with j=.02 and Lambda=2500 even the lower bound A_Omega=0
requires logRref>=301.1814 when logPstar=14. The demo remains logRref=10.
Changing this relation must be followed by recomputing all actual moments
and shared pressure, constructing Eq. (8.2), joining its radial jets and
performing the Section 10 inner moment correction. The extracted equations
and requirements are in `lei_ren_part1_paper_regular_core_plan.md`.

`PaperReferenceMoments(profile)` in
`lei_ren_part1_paper_reference_moments.py` computes all five cumulative
moments of the same velocity candidate from the origin through the initial
axial turnoff. Its temporary power reference is integrable but is not a
smooth axis core. Queries outside the implemented domain fail explicitly.
`SourceStress(profile, moments)` in `lei_ren_part1_paper_source_stress.py`
uses these actual moments, their Z derivatives and the same pressure to
evaluate inertial stress, viscous shear and total stress at arbitrary MP
exponents. Its receipt checks the radial inertial equations as well as
derivative refinement; stress values are not full momentum residuals.
The saved stress receipt checks offsets -1, .5, 2 at Z=.3: the maximum
relative radial inertial-equation defect is 1.06e-10. These samples do not
cover the later pulse or heat region and do not certify the stress cone.

The axial coefficient solve now promotes bump quadrature to at least 192,
and its actual primitive uses that effective order with identical canonical
nodes on complete bump supports. `lei_ren_part1_paper_axial_tail_repair.py`
records the actual linear and energy replay. Any remaining exterior mean is
retained. An arbitrarily tiny nonzero mean cannot certify finite energy on
an unbounded radial domain. Exact continuous closure, all later-stage
moments, a regular core, the stress cone and scale recursion remain open.
