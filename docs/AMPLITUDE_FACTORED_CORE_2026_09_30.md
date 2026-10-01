# Amplitude-factored core and radial positivity — 2026-09-30

The finite core now has an exact amplitude-factored recurrence F=F0(Z) A.
It removes the common exponentially small factor before interval evaluation,
without rescaling only one equation or dropping its nonlinear effects.

The new `lei_ren_part1_paper_amplitude_factored_core.py` applies
F_Z/F0=A_Z+ell A, where ell=F0'/F0. Axial swirl sources retain
S A_i A_j, with S=F0^2. Generated physical pressure rows also retain S;
P[0] is the same supplied analytic preheat datum. Axial derivative depth
is consumed exactly as in the original finite radial recursion.

An independent fixture reconstructs F0 A and compares every retained
F, Uz and P coefficient with the original unfactored recurrence at
Z=-.4,0,.3, radial degree4 and axial depth3. Maximum scaled discrepancy
is below1e-170. Independent pressure-integral checks and scalar-component
conversion also pass. This checks the coupled algebra, not omitted orders.

`lei_ren_part1_paper_factored_core_positivity.py` consumes the accepted
coherent pressure with its whole-axis integral error bounds and the actual
degree18 parameters. It converts each normalized radial polynomial on
x=R/R_a in [0,1] into directed Bernstein coefficients. Their convex hull
bounds every radius, rather than only sampled radii.

At Z=.3, F/F0 has a positive whole-radius lower bound approximately
0.282979986880543. Since the analytic F0 is positive, this establishes
positivity of this finite-model slice. The same slice also now supplies
normalized reciprocal-F axial derivatives through order3 over every radius,
with the inverse axis amplitude kept separate in logarithmic form.
These are axial derivative bounds at one Z, not mixed/global core norms.

The single interval Z=[-1,1] gives an inconclusive Bernstein lower bound
approximately -1.41e607. This is severe interval dependency in the rational
axis-gradient products; it is not an observed zero or a negative velocity.
The common amplitude dependency has been removed, but the H0 and gradient
dependencies remain. The global positivity flag is explicitly false.

The axis helper now accepts an optional axial subinterval while retaining
the same global analytic G bound. The accepted pressure rows are exposed
through one shared input helper, preserving accepted component-mass checks
and the error-receipt domain check for both core representations.

Next work must subdivide the entire axial domain with extra resolution
near the unique H0 root, or retain the rational products symbolically.
Use the factored recurrence and Bernstein bounds on every cell, combine
positive cell lower bounds into reciprocal C3 estimates, and then complete
the frozen-profile K contributions. Infinite radial remainder, functional
five-moment closure, full finite-width cone and true temporal recursion
remain open.
