"""Strict reflection-structure receipt for the frozen Eq. (4.5) seed.

This module certifies one candidate-specific structural property that is not a
Navier--Stokes validation result.  The frozen seed has exactly zero coefficients
for every odd power of eta in both Phi and F.  Together with the existing Eq.
(4.5) similarity relation and streamfunction coupling, that implies

    q(-z)=q(z), X(-z)=X(z), eta(-z)=-eta(z),
    Phi(-eta)=Phi(eta), F(-eta)=F(eta),
    U(-eta)=U(eta), v0(-eta)=-v0(eta),

and therefore cylindrical channel parities

    u_r(x,y,-z,t)     = -u_r(x,y,z,t),
    u_theta(x,y,-z,t) =  u_theta(x,y,z,t),
    w(x,y,-z,t)       =  w(x,y,z,t).

In particular the radial/poloidal velocity vanishes on the midplane z=0 for
this seed.  The coefficient parity is the strict evidence; public-velocity
mirror probes below are an implementation cross-check, not the proof itself.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate


SCHEMA = "eq45_seed_reflection_structure_v1"
FROZEN_CANDIDATE_SHA256 = (
    "48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7"
)
DEFAULT_CANDIDATE_PATH = Path("artifacts/constrained/eq45_velocity_candidate_seed.json")
PUBLIC_CROSSCHECK_TOLERANCE = 2e-11


def _odd_eta_modes(candidate: Eq45VelocityCandidate) -> tuple[dict[str, Any], ...]:
    basis = candidate.profile_basis
    rows: list[dict[str, Any]] = []
    for index, (i, j) in enumerate(basis.mode_indices):
        if j % 2 == 0:
            continue
        rows.append(
            {
                "mode": [int(i), int(j)],
                "phi": float(basis.phi_coefficients[index]),
                "F": float(basis.swirl_coefficients[index]),
            }
        )
    return tuple(rows)


def certify_even_eta_seed(candidate: Eq45VelocityCandidate) -> dict[str, Any]:
    """Certify exact even-in-eta Phi/F coefficient support for one candidate.

    This deliberately requires literal zero odd-eta coefficients.  A small but
    nonzero odd coefficient is a different representation and must not inherit
    the strict reflection certificate merely because a sampled parity defect is
    numerically small.
    """

    if not isinstance(candidate, Eq45VelocityCandidate):
        raise TypeError("candidate must be Eq45VelocityCandidate")
    odd_rows = _odd_eta_modes(candidate)
    nonzero = [row for row in odd_rows if row["phi"] != 0.0 or row["F"] != 0.0]
    if nonzero:
        raise ValueError(f"odd-eta profile coefficients must be exactly zero: {nonzero}")

    return {
        "odd_eta_modes": list(odd_rows),
        "odd_eta_coefficients_exactly_zero": True,
        "derived_profile_parity": {
            "Phi": "even_in_eta",
            "Phi_X": "even_in_eta",
            "Phi_eta": "odd_in_eta",
            "F": "even_in_eta",
            "U": "even_in_eta",
            "v0": "odd_in_eta",
        },
        "derived_physical_channel_parity": {
            "q": "even_in_z",
            "X": "even_in_z",
            "eta": "odd_in_z",
            "u_r": "odd_in_z",
            "u_theta": "even_in_z",
            "w": "even_in_z",
            "midplane_radial_velocity": "exactly_zero_at_z_equals_0",
        },
        "assumptions": [
            "positive single-valued q branch from q-z^2*q^(2h)=1-t",
            "Eq45 compact eta envelope depends on eta^2",
            "Eq45 streamfunction coupling U=Phi+X*Phi_X",
            "Eq45 streamfunction coupling formula for v0",
            "Cartesian Eq45 radial/swirl mixing",
        ],
    }


def _cylindrical_projection_numerators(points: np.ndarray, velocity: np.ndarray):
    x = points[:, 0]
    y = points[:, 1]
    u = velocity[:, 0]
    v = velocity[:, 1]
    radial = x * u + y * v
    swirl = -y * u + x * v
    return radial, swirl


def public_velocity_reflection_crosscheck(
    candidate: Eq45VelocityCandidate,
    *,
    tolerance: float = PUBLIC_CROSSCHECK_TOLERANCE,
) -> dict[str, float]:
    """Cross-check the strict parity consequence through public ``at_points``.

    The test uses cylindrical projection numerators so no division by radius is
    required.  It intentionally samples only the frozen public evaluator and
    does not read profile values or optimizer tensors.
    """

    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("tolerance must be positive and finite")

    points = np.array(
        [
            [0.31, 0.11, 0.06],
            [0.42, -0.17, 0.10],
            [-0.28, 0.24, 0.14],
            [0.15, 0.33, 0.18],
            [-0.37, -0.19, 0.22],
            [0.26, -0.29, 0.12],
        ],
        dtype=float,
    )
    times = np.array([0.25, 0.33, 0.41, 0.50, 0.62, 0.75], dtype=float)
    mirrored = points.copy()
    mirrored[:, 2] *= -1.0

    plus = candidate.at_points(points, times)
    minus = candidate.at_points(mirrored, times)
    radial_plus, swirl_plus = _cylindrical_projection_numerators(points, plus)
    radial_minus, swirl_minus = _cylindrical_projection_numerators(mirrored, minus)

    midplane = points.copy()
    midplane[:, 2] = 0.0
    mid_velocity = candidate.at_points(midplane, times)
    mid_radial, _ = _cylindrical_projection_numerators(midplane, mid_velocity)

    defects = {
        "radial_odd_max_abs": float(np.max(np.abs(radial_minus + radial_plus))),
        "swirl_even_max_abs": float(np.max(np.abs(swirl_minus - swirl_plus))),
        "axial_even_max_abs": float(np.max(np.abs(minus[:, 2] - plus[:, 2]))),
        "midplane_radial_max_abs": float(np.max(np.abs(mid_radial))),
    }
    if not all(np.isfinite(value) for value in defects.values()):
        raise RuntimeError("public reflection cross-check produced nonfinite defects")
    failed = {name: value for name, value in defects.items() if value > tolerance}
    if failed:
        raise RuntimeError(
            f"public velocity does not respect certified reflection structure: {failed}"
        )
    return defects


def audit_frozen_eq45_seed_reflection(
    candidate_path: str | Path = DEFAULT_CANDIDATE_PATH,
) -> dict[str, Any]:
    """Load the frozen artifact and return its strict reflection receipt."""

    candidate = Eq45VelocityCandidate.load_json(candidate_path)
    if candidate.sha256 != FROZEN_CANDIDATE_SHA256:
        raise ValueError(
            "frozen Eq45 candidate identity drifted: "
            f"expected {FROZEN_CANDIDATE_SHA256}, got {candidate.sha256}"
        )

    strict = certify_even_eta_seed(candidate)
    defects = public_velocity_reflection_crosscheck(candidate)
    truth = candidate.to_dict()["truth_boundary"]
    return {
        "schema": SCHEMA,
        "status": "strict_seed_reflection_structure_verified",
        "candidate_sha256": candidate.sha256,
        "strict_evidence": strict,
        "public_velocity_crosscheck": {
            "tolerance": PUBLIC_CROSSCHECK_TOLERANCE,
            **defects,
        },
        "direct_visualization_implication": (
            "for this frozen even-eta seed the poloidal radial channel reverses "
            "across z=0 and is exactly zero on the midplane, while swirl and "
            "axial channels are mirror-even"
        ),
        "representation_feedback": (
            "any target morphology requiring nonzero midplane radial flow needs "
            "an odd-eta profile channel; even-eta eta^4 growth alone cannot add it"
        ),
        "states": {
            "velocity_export_ready": bool(truth["velocity_export_ready"]),
            "visualization_ready": bool(truth["visualization_ready"]),
            "visual_correspondence_verified": bool(
                truth["visual_correspondence_verified"]
            ),
            "physical_support_validated": bool(truth["physical_support_validated"]),
            "pde_validated": bool(truth["pde_validated"]),
            "paper_exact": bool(truth["paper_exact"]),
            "openai_field_identified": bool(truth["openai_field_identified"]),
            "blowup_proved": bool(truth["blowup_proved"]),
        },
        "limitations": (
            "This is a strict structural property of the current autonomous seed "
            "under the existing Eq45 formulas. It is not physical-support, visual-"
            "correspondence, PDE, paper-exact, hidden-field, singularity or blow-up evidence."
        ),
    }
