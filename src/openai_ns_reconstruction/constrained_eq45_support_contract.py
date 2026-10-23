"""Physical-support compatibility audit for Eq. (4.5) similarity profiles.

The Eq. (4.5) lane parameterizes profiles in similarity coordinates ``(X, eta)``.
Compactness in those coordinates must not be silently promoted to the CR001
physical-space support claim.  This module checks that coordinate mapping only;
it does not evaluate a velocity field or establish PDE/visual correspondence.

For the repository-governed positive-exponent relation

    q - z**2 * q**(2*h) = 1 - t,
    eta = z / q**(1/2-h),

one has on the physical branch ``q > 0``

    q * (1 - eta**2) = 1 - t.

Therefore any profile support reaching ``|eta| = 1`` has unbounded physical
extent as the branch edge is approached.  A finite ``eta`` cutoff below one can
instead be mapped to conservative physical radial/axial extents over the
registered time window.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Any


@dataclass(frozen=True)
class Eq45SupportAudit:
    X_cutoff: float
    eta_cutoff: float
    h: float
    time_interval: tuple[float, float]
    physical_radial_support: float
    physical_axial_support: float
    finite_physical_extent: bool
    q_max: float | None
    radial_max: float | None
    axial_max: float | None
    eta_zero_radial_witness: float
    max_X_cutoff_at_eta_zero: float
    similarity_support_alone_satisfies_physical_support: bool
    reasons: tuple[str, ...]
    claim_scope: str = "eq45_similarity_to_physical_support_compatibility_only"
    velocity_export_ready: bool = False
    visual_correspondence_verified: bool = False
    pde_validated: bool = False
    paper_exact: bool = False

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _positive_finite(name: str, value: float) -> float:
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return result


def audit_eq45_similarity_support(
    *,
    X_cutoff: float,
    eta_cutoff: float,
    h: float,
    time_interval: tuple[float, float],
    physical_radial_support: float,
    physical_axial_support: float,
    require_physical_support: bool = False,
) -> Eq45SupportAudit:
    """Map a rectangular ``(X, eta)`` support to CR001 physical support.

    The returned extents are conservative over ``0 <= X <= X_cutoff`` and
    ``|eta| <= eta_cutoff``.  The current Eq. (4.5) regime requires
    ``0 <= h < 1/2`` so that the declared axial similarity exponent is positive.

    ``require_physical_support=True`` turns incompatibility into a fail-closed
    exception suitable for integration gates.  It does not alter any threshold.
    """
    X_cutoff = _positive_finite("X_cutoff", X_cutoff)
    eta_cutoff = _positive_finite("eta_cutoff", eta_cutoff)
    radial_support = _positive_finite("physical_radial_support", physical_radial_support)
    axial_support = _positive_finite("physical_axial_support", physical_axial_support)

    h = float(h)
    if not math.isfinite(h) or not (0.0 <= h < 0.5):
        raise ValueError("h must be finite and lie in [0, 0.5)")

    if len(time_interval) != 2:
        raise ValueError("time_interval must contain exactly two values")
    t_min, t_max = (float(time_interval[0]), float(time_interval[1]))
    if not (math.isfinite(t_min) and math.isfinite(t_max)):
        raise ValueError("time_interval must be finite")
    if not (0.0 <= t_min < t_max < 1.0):
        raise ValueError("time_interval must satisfy 0 <= t_min < t_max < 1")

    # At eta=0 the source relation reduces exactly to q=1-t.  This gives a
    # finite radial witness even when eta_cutoff reaches the |eta|=1 branch edge.
    q_eta_zero_max = 1.0 - t_min
    eta_zero_radial_witness = math.sqrt(2.0 * q_eta_zero_max * X_cutoff)
    max_X_cutoff_at_eta_zero = radial_support**2 / (2.0 * q_eta_zero_max)

    reasons: list[str] = []
    if eta_zero_radial_witness > radial_support:
        reasons.append("eta_zero_radial_extent_exceeds_registered_support")

    if eta_cutoff >= 1.0:
        finite_extent = False
        q_max = None
        radial_max = None
        axial_max = None
        reasons.append("eta_cutoff_reaches_or_crosses_unbounded_physical_branch_edge")
        compatible = False
    else:
        finite_extent = True
        # q=(1-t)/(1-eta^2), maximized at earliest t and largest |eta|.
        q_max = (1.0 - t_min) / (1.0 - eta_cutoff**2)
        radial_max = math.sqrt(2.0 * q_max * X_cutoff)
        axial_max = eta_cutoff * q_max ** (0.5 - h)
        if radial_max > radial_support:
            reasons.append("radial_extent_exceeds_registered_support")
        if axial_max > axial_support:
            reasons.append("axial_extent_exceeds_registered_support")
        compatible = radial_max <= radial_support and axial_max <= axial_support

    if compatible:
        reasons.append("similarity_support_maps_inside_registered_physical_support")

    report = Eq45SupportAudit(
        X_cutoff=X_cutoff,
        eta_cutoff=eta_cutoff,
        h=h,
        time_interval=(t_min, t_max),
        physical_radial_support=radial_support,
        physical_axial_support=axial_support,
        finite_physical_extent=finite_extent,
        q_max=q_max,
        radial_max=radial_max,
        axial_max=axial_max,
        eta_zero_radial_witness=eta_zero_radial_witness,
        max_X_cutoff_at_eta_zero=max_X_cutoff_at_eta_zero,
        similarity_support_alone_satisfies_physical_support=compatible,
        reasons=tuple(reasons),
    )

    if require_physical_support and not compatible:
        raise ValueError(
            "Eq. (4.5) similarity support does not by itself satisfy the "
            "registered physical support; add a governed physical taper/exterior "
            "connection or choose compatible similarity cutoffs"
        )
    return report
