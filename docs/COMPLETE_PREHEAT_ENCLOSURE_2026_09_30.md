# Complete stored-parameter preheat pressure norm bound

The pressure datum now has a conservative directed interval norm bound
including the reference contribution, all finite pre-collar stages, the
H=1 heat connection, and the infinite H=1 exterior. This is the same shared
preheat pressure definition used by the existing reconstruction.

## Terminal formula and normalization

Let lambda=1+delta and t=y-y_tail. The stored source normalization gives

    S=exp(2(log_c_inf-logPstar)-log(2)-lambda logR_tail).

The normalized positive terminal density is S exp(-lambda t) K0(t)^2.
For 0<=t<=3, 1-epsilon<=K0<=1; beyond3, K0=1 exactly. Hence

    S(1-epsilon)^2 (1-exp(-3lambda))/lambda <= I_collar
        <= S (1-exp(-3lambda))/lambda,
    I_exterior=S exp(-3lambda)/lambda.

The interval arithmetic evaluates the actual stored log_c_inf. A separately
derived endpoint normalization is reported only as a residual diagnostic;
the source normalization is never reset. Both terminal pieces are independent
of Z after replacing H by1, so their axial derivative contributions vanish.
This replacement defines the analytic preheat datum; it does not certify the
actual heat-dependent velocity field or its energy.

## Actual result and independent check

`lei_ren_part1_paper_complete_preheat_enclosure_check.py` combines the
existing eleven stage upper bounds with the reference and terminal pieces.
On |Z|<=.8 the normalized absolute derivative upper bounds for orders0..3
are approximately 3.44587, 11.0268, 66.7120 and 537.225. Thus

    ||P0/Pstar^2||_C2 <= 47.8287,

where C2 means sup|f|+sup|f_Z|+sup|f_ZZ|/2. The machine receipt retains
the computed upper bound 47.82861181496495... and Pstar^2=exp(28).
Tiny terminal masses stay separately labeled rather than being recovered
by subtracting rounded totals.

`lei_ren_part1_paper_terminal_preheat_enclosure_fixture.py` independently
integrates a resolved collar and checks the exterior formula. It deliberately
shifts the stored log_c_inf by .01 and confirms that the retained scale
changes by exp(.02), demonstrating that normalization is preserved.
The checks pass. Independent quadrature is a regression comparison, not
the proof of the collar bounds.

The interval JSON encoder now includes exact binary MP endpoint tuples
(sign, mantissa, exponent, bitcount) as well as human-readable decimal strings.
The decimal strings are rounded display values; reconstruct certified
endpoints from the exact tuples, never from the display values. Scalar
upper bounds also retain their exact MP tuples.

## Remaining scope

This bounds the analytic integral defined by the stored schedule parameters.
The original exp/log parameter derivations are not enclosed. The bound is
on pressure magnitude and derivatives, not the approximation error of the
Gauss representation. It does not establish a small centered five-defect
norm, core/RK error propagation, uniform inverse closure, finite energy,
the actual heat exterior, temporal recursion or oscillatory correction.

Next derive narrow interval integrals for variable radial stages, compare
them to the retained Gauss values with explicit derivative error budgets,
and propagate those pressure errors into the common core and all five
centered moment defects. Relative-flat atoms must remain separate.
