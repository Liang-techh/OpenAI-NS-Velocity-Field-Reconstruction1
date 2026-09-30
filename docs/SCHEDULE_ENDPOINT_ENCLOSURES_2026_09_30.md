# Directed enclosures for stored-schedule preheat stages

The continuous preheat schedule endpoint evaluator now uses directed interval
arithmetic and exact outer branches of the flat switch primitive. This
removes endpoint primitive quadrature from the pre-collar mass-bound inputs.

For sigma(x)=exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)), extended by 0 and 1,
sigma(x)+sigma(1-x)=1. Thus J(1)=1/2 exactly, J(x)=0 for x<=0 and
J(x)=x-1/2 for x>=1. Named schedule endpoints do not require genuine partial
primitive integrals. Tiny interval overlap with 0 or 1 is enclosed with
monotonicity and the Lipschitz bound 0<=J'<=1. A general partial input is
conservatively bounded rather than evaluated with uncertified quadrature.

The normalized log-amplitude is evaluated using the same four primitives
and the same stored mu, delta, y_d, y_rel and Ts as the continuous schedule.
Stage bounds enforce the paper gate 0<mu<=1/60 and require 0<delta<1. Negative-slope stages
use an exponential integral cap

    mass <= density_factor * exp(2 ell(left)) / (-2 slope_upper).

This avoids enormous length-times-amplitude estimates on the very long
pulse stage. Transition stages use conservative slope envelopes; no exact
constant slope is assumed across a rounded switching boundary. The density
prefactor is 1/2 before and during flatten, and 1/8 after flatten. Bounds of
axial derivatives through order three use the positive-factor envelopes.

## Actual receipt

`lei_ren_part1_paper_schedule_endpoint_enclosures_check.py` writes enclosures
for 11 finite stages before heat_connection, plus named endpoint data.
The largest normalized log-amplitude endpoint interval width is approximately
2.39823e-205. All positive stage mass upper bounds remain separate,
including the very small flatten and later-stage masses. Zero-length waiting
is handled explicitly. The reference contribution on the inner power branch
is not part of these 11 finite-stage bounds.

## Scope and next step

The interval arithmetic encloses the continuous formulas **relative to the
stored Decimal schedule parameters**. It does not enclose the calculation
of those parameters from exp/log and the original paper constants. The
legacy floating quadrature value of J(1) is not used; this adapter targets
the shared continuous schedule with exact J(1)=1/2.

These are conservative positive-integral bounds, not tight Gauss quadrature
error estimates. The finite-atom approximation and the integral upper bound
must not be confused. Heat-collar/exterior normalization, core and RK source
propagation, full five-defect norms and endpoint strips are still open.
Next narrow variable-stage integrals using interval subdivision, preserve
derivative remainders of variable-beta stages, and add heat/exterior bounds
before assembling a full pressure-datum error enclosure.
