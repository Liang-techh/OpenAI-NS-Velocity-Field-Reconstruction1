# Solenoidal anisotropic scale transport

`experiments/root_st073/solenoidal_scale_transport.py` supplies a velocity-only
mapping between scales. It is not a time integrator or an accepted NS step.

For scale s in (0,1] and h=0.005, define

    S = diag(s^(1/2), s^(1/2), s^(1/2-h))
    u_P(x) = s^(1/2-h) S u_ref(S^-1 x) / det(S).

The Cartesian horizontal components scale by s^(-1/2), the axial component
by s^(-1/2-h). Divergence transforms as

    div u_P(x) = s^(1/2-h) / det(S) * div u_ref(S^-1 x).

To give the axisymmetric mean swirl the stronger s^(-1/2-h) amplification,
add (s^(-1/2-h)-s^(-1/2)) times the reference mean swirl at S^-1 x in the
azimuthal direction. This addition has zero divergence because its scalar
profile depends only on r and z. Fixed absolute quadrature angles ensure
that property even if the angular mean is only approximately integrated.
Five angles integrate the current cylindrical modes |m|<=2 exactly, up to
roundoff. Applying the stronger swirl factor to every nonaxisymmetric
azimuthal component would not generally preserve divergence.

For each fixed positive s, bounded linear coordinate/velocity factors
preserve finite energy and map compact support to S times its old support.
This does not establish a uniform bound as s tends to zero or a critical-time
extension. The map composes by multiplying scales; this is a kinematic
property, not dynamical recursive closure.

## Interface and evidence

From the experiment directory:

```python
from global_two_patch_candidate import load
from solenoidal_scale_transport import VelocityAdapter, SolenoidalScaleTransport
field, localized, snapshot, _, _ = load()
reference = VelocityAdapter(field, snapshot['inputs']['mean']['tau'])
next_scale = SolenoidalScaleTransport(reference, 0.5, localized.inner.h)
velocity = next_scale.velocity(points)  # Cartesian N by 3 array
```

The frozen report records the source hashes. Five-versus-twelve-angle mean
swirl differs by at most 3.91e-12 on the selected probes. The numerical
composition T_0.5(T_0.8 u) versus T_0.4 u has relative difference 4.32e-15.
Mapped probe values are finite and nonzero, and sampled outside-support
values are zero. Cartesian finite differences agree with the Piola
divergence prediction to maximum 1.26e-5 and 1.11e-5 at two step sizes;
this is a bounded sampling check, not a continuum divergence certificate.

The physical-time derivative is now implemented in
`experiments/root_st073/scale_transport_generator.py`. With t=-tau,
s=tau/tau0, reference Cartesian Jacobian J, and mean swirl m, it is

```
ut = (diag(0.5,0.5,0.5+h) u
      + J @ (x/2,y/2,(0.5-h)z) + h*m*e_theta) / tau0.
```

The frozen generator report checks nine axis/interior/transition/exterior
points: halving the Cartesian Jacobian step changes the generator relatively
by 3.24e-10; forward physical-time differences show first-order convergence
with error ratios approaching 2. This checks the derivative of the chosen
map, not the momentum equation.

No pressure transformation or force has been accepted with this mapping.
Anisotropic transport is not asserted to be a symmetry
of the fixed-viscosity NS equation. The remaining task is to construct
pressure and evolution corrections and measure the full dynamic residual
between successive mapped fields. No recursion step is accepted here.

## Uniform energy bound for the prescribed velocity map

The kinematic map has a stronger energy property than finiteness at each
fixed scale. Write the reference in cylindrical components, with m(r,z)
the angular mean swirl, and split

```
H = ur*e_r + (utheta-m)*e_theta
Z = m*e_theta + uz*e_z
T_s u(x) = s^(-1/2) H(S_s^-1 x)
           + s^(-1/2-h) Z(S_s^-1 x).
```

For a true angular mean, H and Z are orthogonal after integration over
the whole angle. With E_H=integral |H|^2 and E_Z=integral |Z|^2 over R^3,
the change of variables det(S_s)=s^(3/2-h) gives exactly

```
||T_s u||_L2^2 = s^(1/2-h)*E_H + s^(1/2-3h)*E_Z.
```

Thus for the current h=0.005 and 0<s<=1, this squared norm is bounded
by E_H+E_Z and tends to zero as s tends to zero. Kinetic energy is half
this expression. The reference is smooth and compactly supported, so both
constants are finite. The implemented five-angle projection agrees with
the true angular mean for the current resolved cylindrical modes |m|<=2.
For a more general smooth reference with a discrete mean, the cross term
need not vanish, but the bound
`||T_s u||^2 <= 2*s^(1/2-h)*E_H + 2*s^(1/2-3h)*E_Z` still applies whenever
these two reference norms are finite.

This is an analytic energy bound for the prescribed velocity family only.
It does not establish a smooth force, pressure balance, or critical-time
extension of an NS solution. Future dynamically fitted corrections must
satisfy their own uniform energy bounds; this estimate cannot be inherited
automatically by an independently optimized sequence of reference fields.
