"""Reference-velocity correction for scale-compatible NS momentum.

All inputs are Cartesian values on the same points; derivatives may come
from analytic curl bases. Pressure and forcing are held fixed unless an
explicit pressure-gradient correction is supplied. No PDE acceptance here.
"""
import numpy as np

from scale_transport_generator import solenoidal_scale_generator


def velocity_correction_response(points, velocity, jacobian, delta_velocity,
                                 delta_jacobian, delta_laplacian, tau0, h,
                                 delta_mean_swirl, nu,
                                 delta_pressure_gradient=None):
    """Return linear response L(v) and exact quadratic remainder (v.grad)v.

    For R(u,p)=G(u)+Du*u-nu*lap(u)+grad(p), G is the linear
    solenoidal scale generator. Therefore
    R(u+a*v,p+a*q)-R(u,p)=a*L(v,q)+a**2*(Dv*v).
    This changes the reference velocity itself, unlike time-tangent patches.
    The caller must construct v divergence-free and preserve support/gates.
    """
    u = np.asarray(velocity, dtype=float)
    j = np.asarray(jacobian, dtype=float)
    v = np.asarray(delta_velocity, dtype=float)
    dj = np.asarray(delta_jacobian, dtype=float)
    lap = np.asarray(delta_laplacian, dtype=float)
    if u.shape != v.shape or lap.shape != v.shape:
        raise ValueError('velocity, correction, and Laplacian shapes must agree')
    if j.shape != dj.shape or j.shape != (len(u), 3, 3):
        raise ValueError('Jacobians must have matching (N,3,3) shapes')
    linear = solenoidal_scale_generator(
        points, v, dj, tau0, h, delta_mean_swirl)
    linear += np.einsum('nij,nj->ni', j, v)
    linear += np.einsum('nij,nj->ni', dj, u) - float(nu) * lap
    if delta_pressure_gradient is not None:
        gp = np.asarray(delta_pressure_gradient, dtype=float)
        if gp.shape != v.shape:
            raise ValueError('pressure gradient must match velocity shape')
        linear += gp
    quadratic = np.einsum('nij,nj->ni', dj, v)
    return linear, quadratic
