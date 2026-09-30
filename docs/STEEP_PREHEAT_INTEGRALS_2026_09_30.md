# High-order enclosures for both steep preheat transitions

The shared high-order integrator now covers the two post-flatten steep
transitions, with the same directed Taylor remainder machinery and the same
stored source normalization. The normalized angular pressure density after
flatten has prefactor1/8 and is independent of Z.

## Local formulas

At steep_transition_in, write x=y-y_rel. The normalized log-amplitude is

    ell(y)=ell(y_rel)+(-.5-mu)x-(1-mu)J(x).

At steep_transition_out, the stored left checkpoint y_q can differ from
y_rel+1+Ts through Decimal rounding. Retain

    shift=y_q-y_rel-1-Ts,
    ell(y)=ell(y_q)-1.5x+(1-delta/2)[J(x+shift)-J(shift)].

The local stage length is likewise retained as the interval of the stored
right-minus-left checkpoints. It is not replaced by an assumed unit length.
The rising slope correction is represented by a negative J coefficient.
Flat-edge exponential correction intervals handle both coefficient signs;
the positive-sign-only multiplier from the first transitions is insufficient
for the terminal slope restoration.

## Actual evidence

`lei_ren_part1_paper_steep_preheat_integrals_check.py` compares order12
enclosures at32 and64 radial panels with the existing192-point Gauss
representation. Both interval widths shrink. At64 panels, conservative
relative mass error upper bounds are approximately

- steep_transition_in: 1.69391e-12;
- steep_transition_out: 1.09072e-12.

The absolute masses and errors have enormous negative exponents and remain
separate atoms in the machine receipt. These are relative errors against the
positive integral lower endpoint, not a ratio to a rounded total pressure.
The axial derivative error budgets are exactly zero because beta=0 on both
stages in the analytic preheat construction. The first two transition
regression receipts also pass after this extension, with their earlier
error bounds unchanged.

## Remaining work

The variable-beta flatten integral and its axial derivative errors still need
a separate enclosed treatment. Pressure approximation error must then be
assembled from all labeled stage errors and propagated to the common core,
RK bridges, restoration and centered moment defects. These checks are
relative to stored parameters and do not close the exact paper parameter
derivations, uniform source bounds, five-moment inverse or temporal recursion.
