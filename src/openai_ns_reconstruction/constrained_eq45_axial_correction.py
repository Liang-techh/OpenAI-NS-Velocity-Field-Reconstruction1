"""Minimal symmetry-preserving axial correction lift for the Eq. (4.5) candidate.

The frozen Eq45 candidate currently uses an autonomous compact Phi/F basis with
eta degree two.  Agent-7's candidate-specific local sensitivity screen identified
the even ``(0,4)`` Phi and swirl modes as two interpretable axial correction
channels.  This module turns only those two screened directions into a concrete
representation lift without changing the radial basis or activating every new
coefficient exposed by a blanket degree increase.

A zero correction is a function-preserving reparameterization: every existing
coefficient is copied by its ``(i,j)`` mode label, all newly exposed modes are
zero, and the returned public ``velocity(points,time)`` is therefore unchanged.
Nonzero correction values are autonomous bounded design choices.  Applying the
lift does not establish visual correspondence, physical support, PDE validity,
paper exactness, hidden-field recovery, singularity, or blow-up.
"""
from __future__ import annotations

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_profile_basis import Eq45CompactProfileBasis


TARGET_ETA_DEGREE = 4
TARGET_MODE = (0, 4)


def _finite_bounded(value: float, *, limit: float, name: str) -> float:
    value = float(value)
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite")
    if abs(value) > limit:
        raise ValueError(f"{name} exceeds profile coefficient_limit={limit}")
    return value


def lift_eq45_axial_eta4(
    candidate: Eq45VelocityCandidate,
    *,
    phi_eta4: float = 0.0,
    swirl_eta4: float = 0.0,
) -> Eq45VelocityCandidate:
    """Embed ``candidate`` into eta degree four and activate only ``(0,4)``.

    The source basis must not already contain eta-degree-four modes.  Existing
    Phi/F coefficients are copied by mode identity, not by flat-array position.
    Every newly introduced coefficient remains exactly zero except the two
    explicitly supplied ``Phi(0,4)`` and ``F(0,4)`` values.

    The returned object is an ordinary :class:`Eq45VelocityCandidate`, so its
    existing JSON save/load, canonical hash, public velocity API, nontriviality
    check, and truth-boundary metadata continue to apply unchanged.
    """
    if not isinstance(candidate, Eq45VelocityCandidate):
        raise TypeError("candidate must be Eq45VelocityCandidate")

    base = candidate.profile_basis
    if base.eta_degree >= TARGET_ETA_DEGREE:
        raise ValueError(
            "source profile already contains eta-degree-four modes; "
            "refuse to overwrite an existing correction channel"
        )

    phi_eta4 = _finite_bounded(
        phi_eta4, limit=base.coefficient_limit, name="phi_eta4"
    )
    swirl_eta4 = _finite_bounded(
        swirl_eta4, limit=base.coefficient_limit, name="swirl_eta4"
    )

    old_phi = dict(zip(base.mode_indices, base.phi_coefficients, strict=True))
    old_swirl = dict(zip(base.mode_indices, base.swirl_coefficients, strict=True))

    new_modes = tuple(
        (i, j)
        for i in range(base.radial_degree + 1)
        for j in range(TARGET_ETA_DEGREE + 1)
    )
    phi = [float(old_phi.get(mode, 0.0)) for mode in new_modes]
    swirl = [float(old_swirl.get(mode, 0.0)) for mode in new_modes]

    target_index = new_modes.index(TARGET_MODE)
    phi[target_index] = phi_eta4
    swirl[target_index] = swirl_eta4

    lifted_basis = Eq45CompactProfileBasis(
        radial_degree=base.radial_degree,
        eta_degree=TARGET_ETA_DEGREE,
        phi_coefficients=tuple(phi),
        swirl_coefficients=tuple(swirl),
        x_cut=base.x_cut,
        eta_cut=base.eta_cut,
        cutoff_power=base.cutoff_power,
        coefficient_limit=base.coefficient_limit,
    )
    return Eq45VelocityCandidate(
        profile_basis=lifted_basis,
        h=candidate.h,
        time_start=candidate.time_start,
        time_end=candidate.time_end,
    )
