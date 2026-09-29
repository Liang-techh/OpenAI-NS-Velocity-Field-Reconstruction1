"""Explicit Part I interfaces, arXiv:2609.35406v1, Sections 2, 8 and 17.

These algebraic helpers do not construct the paper's profiles or certify
admissibility/flatness. Physical stress inputs must include cutoff derivatives.
"""
import numpy as np

SOURCE = 'https://arxiv.org/html/2609.35406v1'


def parameter_map(h):
    """Repo q=lambda^2, X=R, eta=Z after x/sqrt(nu); delta=2h."""
    delta = 2.0*float(h)
    return dict(h=float(h), delta=delta, A=(1+delta)/2, D=(1-delta)/2,
                in_part1_stated_delta_range=bool(0 < delta < 1/200),
                coordinates='q=lambda^2; X=R; eta=Z; source x=physical x/sqrt(nu)')


def completed_stress_tensor(points, angular_stress, axial_stress, axial_z_derivative):
    """Eq.17.24 in Cartesian coordinates; annular input r>0 required.

    The theta-theta entry r*d_z(Tz) cancels the unwanted radial divergence.
    This tensor is not the physical viscous Cauchy stress -pI+nu(Du+Du^T).
    """
    p = np.asarray(points, float)
    if p.ndim != 2 or p.shape[1] != 3 or not np.all(np.isfinite(p)):
        raise ValueError('points must be finite (N,3)')
    r = np.hypot(p[:, 0], p[:, 1])
    if np.any(r <= 0):
        raise ValueError('Annular interface requires r>0; axis extension is separate')
    a, b, bz = np.broadcast_arrays(angular_stress, axial_stress, axial_z_derivative, r)[:3]
    if not all(np.all(np.isfinite(v)) for v in (a, b, bz)):
        raise ValueError('Stress jets must be finite')
    er = np.column_stack((p[:, 0]/r, p[:, 1]/r, np.zeros(len(r))))
    et = np.column_stack((-p[:, 1]/r, p[:, 0]/r, np.zeros(len(r))))
    ez = np.zeros_like(er); ez[:, 2] = 1
    outer = lambda x, y: np.einsum('ni,nj->nij', x, y)
    return (a[:, None, None]*(outer(er, et)+outer(et, er))
            + b[:, None, None]*(outer(er, ez)+outer(ez, er))
            + (r*bz)[:, None, None]*outer(et, et))


def projected_divergence(r, angular_stress, axial_stress, angular_r_derivative, axial_r_derivative):
    """D T in the (r,theta,z) frame, Sections 2.4 and 17.2."""
    r, a, b, ar, br = np.broadcast_arrays(r, angular_stress, axial_stress,
                                         angular_r_derivative, axial_r_derivative)
    if np.any(r <= 0):
        raise ValueError('Annular interface requires r>0')
    return np.stack((np.zeros_like(r), ar+2*a/r, br+b/r), axis=-1)


def background_remainder(full_cylindrical_residual, stress_divergence):
    """E_B=R_B+D T; preserves the radial residual, never labels it zero."""
    residual, divergence = np.asarray(full_cylindrical_residual), np.asarray(stress_divergence)
    if residual.shape != divergence.shape or residual.shape[-1] != 3:
        raise ValueError('Both inputs must have matching (...,3) shapes')
    return residual+divergence


def linear_core_axis_slopes(Z, delta, F0, F0_Z, Uz0, Uz0_Z, P0, P0_Z):
    """Eq.8.7, exact axis slopes shared by the linear model and regular core.

    P0 is supplied by the same assembled outer profile, not freely fitted.
    The linear model alone is not a nonlinear stress-free solution.
    """
    Z = np.asarray(Z, float)
    if not (0 < delta < 1) or np.any(np.abs(Z) > 1) or np.any(np.asarray(F0) <= 0):
        raise ValueError('Require 0<delta<1, |Z|<=1 and F0>0')
    d, L = 1-Z**2, 1-delta*Z**2
    H0 = (1-delta)*Z/2+d*Uz0
    Fr = ((1+delta/2-Z*Uz0-d*Uz0_Z)*F0+H0*F0_Z)/(4*L)
    Uzr = ((1+delta)*(1-2*Z*Uz0)*Uz0/2+H0*Uz0_Z+d*P0_Z-2*(1+delta)*Z*P0)/(2*L)
    return Fr, Uzr


def sector_tau_bounds(lam, eta_limit):
    """Eq.17.28; excludes the endpoints where uniform tau-flatness fails."""
    if not 0 <= eta_limit < 1:
        raise ValueError('A fixed interior sector |Z|<=eta_limit<1 is required')
    lam = np.asarray(lam, float)
    if np.any(lam <= 0):
        raise ValueError('lambda must be positive')
    return (1-eta_limit**2)*lam**2, lam**2
