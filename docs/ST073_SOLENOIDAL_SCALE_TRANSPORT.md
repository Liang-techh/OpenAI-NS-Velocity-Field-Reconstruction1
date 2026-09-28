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

No pressure transformation, force, or time derivative has been accepted
with this mapping. Anisotropic transport is not asserted to be a symmetry
of the fixed-viscosity NS equation. The remaining task is to construct
pressure and evolution corrections and measure the full dynamic residual
between successive mapped fields. No recursion step is accepted here.
