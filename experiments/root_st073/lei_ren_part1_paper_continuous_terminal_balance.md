# Installed terminal balance and radial energy

Run `python experiments/root_st073/lei_ren_part1_paper_continuous_terminal_balance.py`
to evaluate the current continuous incoming/exterior adapter, including actual
inner seeds and coefficient Z tangents. This extends the earlier formal
quadrature graph to the installed continuous field. It retains every nominal
terminal mass, mass derivative, and signed pulse omission bound.

Beyond axial support, the source recovery formula gives

\[
C(Z)=\frac{(1-\delta)Z M_z+(1-Z^2)M_{z,Z}}{1-\delta Z^2},
\qquad u_r=-\frac{\nu C(Z)}r.
\]

At fixed time gap \(\tau=T-t\), the physical axial Jacobian is

\[
\frac{dz}{dZ}=\sqrt\nu\,\tau^{(1-\delta)/2}
\frac{1-\delta Z^2}{(1-Z^2)^{(3-\delta)/2}}.
\]

For kinetic energy defined by one half of the volume integral of squared
velocity, the radial component contributes

\[
\frac{dE_r}{dZ\,d\log r}=\pi\nu^2 C(Z)^2\frac{dz}{dZ}.
\]

A nonzero continuous coefficient therefore makes the integral over unbounded
radial logarithmic distance diverge. Small normalized residuals cannot prove
finite global energy. The receipt records actual mass and C relative to Rh
in signed logarithms, so it does not hide the enormous source scale factor.

The normalized balance decomposes incoming mass, pulse mass, and two end-bump
contributions, and separately decomposes their Z derivatives. It records both
independent term replay and the installed cumulative primitive. Their
difference is retained because finite-precision summation order changes a
tiny cancellation residual. The physical energy calculation uses the actual
installed mean, not a more favorable alternate replay.

The transported pulse bound covers omitted kernel pieces only, conditional
on nominal coefficients. It excludes quadrature, incoming-input, coefficient,
and subtraction errors. None of these numerical results certifies exact
terminal closure. The next correction must use a single continuous integral
definition for complete atoms, coefficient constraints, and partial primitives,
and preserve both the incoming seed and its Z derivative. Assigning terminal
mass to zero or imposing a radial cutoff would bypass that requirement.
