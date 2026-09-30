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
