"""Explicit odd-poloidal seed for the source's central axial-stretch pattern.

Uses an existing Phi(0,1) mode, with unchanged swirl and physical support.
This is a qualitative design choice, not a numerical identification of the
OpenAI field or a pressure/forcing/PDE fit. It does not replace the default.
"""
from dataclasses import replace
import numpy as np
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate


def bipolar_seed(base: Eq45SupportedVelocityCandidate, phi01=1.):
    if not isinstance(base,Eq45SupportedVelocityCandidate):
        raise TypeError('base must be a supported Eq45 candidate')
    basis=base.parent.profile_basis
    if isinstance(phi01,(bool,np.bool_)) or not np.isscalar(phi01) or not np.isfinite(phi01):
        raise ValueError('phi01 must be finite')
    if not 0<float(phi01)<=basis.coefficient_limit:
        raise ValueError('phi01 must be positive and within existing coefficient limit')
    phi=[0.]*basis.mode_count
    phi[basis.mode_indices.index((0,1))]=float(phi01)
    return replace(base,parent=replace(base.parent,profile_basis=replace(
        basis,phi_coefficients=tuple(phi))))
