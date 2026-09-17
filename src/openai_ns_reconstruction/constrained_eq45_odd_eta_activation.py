"""Minimal bounded activation of the existing Eq. (4.5) odd-eta modes.

The frozen Eq45 seed already contains the tensor-product modes ``(0, 1)`` in
both ``Phi`` and ``F``, but initializes both coefficients to zero.  This helper
only changes those existing coefficients.  It does not grow the basis, fit a
visual target, alter forcing/pressure, or promote any scientific readiness
claim.

The returned object is an ordinary :class:`Eq45VelocityCandidate`, so its
existing public ``velocity/at_points/grid`` evaluator and JSON serialization
remain the only numerical delivery path.
"""
from __future__ import annotations

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_profile_basis import Eq45CompactProfileBasis


ODD_ETA_MODE = (0, 1)


def _bounded_optional(value, *, name: str, limit: float) -> float | None:
    if value is None:
        return None
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(f"{name} must be a finite scalar or None")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite scalar or None") from exc
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite")
    if abs(result) > limit:
        raise ValueError(f"{name} exceeds coefficient_limit={limit}")
    return result


def activate_eq45_odd_eta01(
    candidate: Eq45VelocityCandidate,
    *,
    phi01: float | None = None,
    swirl01: float | None = None,
) -> Eq45VelocityCandidate:
    """Return ``candidate`` with only existing ``(0,1)`` odd-eta values changed.

    ``None`` preserves the corresponding current coefficient.  At least one
    coefficient must be supplied explicitly.  Values remain inside the source
    basis' existing hard coefficient bound; all basis-shape/support parameters,
    every other coefficient, ``h`` and the delivery time interval are copied
    unchanged.

    This is a representation-capability operation, not a fit.  In particular,
    a nonzero ``phi01`` can remove the frozen seed's strict midplane-poloidal
    reflection constraint, while ``swirl01`` independently permits z-odd swirl.
    Whether either change improves correspondence to a public image is a
    separate measurement.
    """
    if not isinstance(candidate, Eq45VelocityCandidate):
        raise TypeError("candidate must be Eq45VelocityCandidate")
    if phi01 is None and swirl01 is None:
        raise ValueError("provide phi01 and/or swirl01 explicitly")

    basis = candidate.profile_basis
    try:
        index = basis.mode_indices.index(ODD_ETA_MODE)
    except ValueError as exc:
        raise ValueError("candidate profile basis does not contain mode (0, 1)") from exc

    phi_value = _bounded_optional(
        phi01, name="phi01", limit=float(basis.coefficient_limit)
    )
    swirl_value = _bounded_optional(
        swirl01, name="swirl01", limit=float(basis.coefficient_limit)
    )

    phi = list(basis.phi_coefficients)
    swirl = list(basis.swirl_coefficients)
    if phi_value is not None:
        phi[index] = phi_value
    if swirl_value is not None:
        swirl[index] = swirl_value

    updated_basis = Eq45CompactProfileBasis(
        radial_degree=basis.radial_degree,
        eta_degree=basis.eta_degree,
        phi_coefficients=tuple(phi),
        swirl_coefficients=tuple(swirl),
        x_cut=basis.x_cut,
        eta_cut=basis.eta_cut,
        cutoff_power=basis.cutoff_power,
        coefficient_limit=basis.coefficient_limit,
    )
    return Eq45VelocityCandidate(
        profile_basis=updated_basis,
        h=candidate.h,
        time_start=candidate.time_start,
        time_end=candidate.time_end,
    )
