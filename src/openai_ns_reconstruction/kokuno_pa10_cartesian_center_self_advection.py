"""Analytic self-advection for the executable PA.10 Cartesian center velocity.

Pinned public provenance:
KokunoYumeto/yang-mills-interacting-workbench@
143f6773feb424ad9ed3a8d116653200f20346b7,
navier-stokes/navier_stokes_workbench.tex, corrected 2026-09-09 reader.

Agent-1 #819 exposes the production analytic Jacobian

    J[i,j] = partial_j u_i

for the existing inner PA.10 contraction-center velocity.  This module adds the
next candidate-facing nonlinear building block

    (u . grad)u_i = sum_j u_j J[i,j]

by contracting that analytic Jacobian with the same public velocity.  No finite
difference is used in the production self-advection path.  Directional finite
differences appear only in focused/report cross-checks.

This remains the *inner PA.10 contraction center*.  It is not the corrected
fixed point, not an outer/global join, not a complete Navier--Stokes residual,
and not independent PDE validation.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa10_cartesian_center_spatial_derivatives import (
    KokunoPA10CartesianCenterSpatialDerivatives,
)
from .kokuno_pa10_cartesian_center_velocity import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)


SCHEMA = "kokuno-pa10-cartesian-center-self-advection-v1"

_SOURCE_FORMULAS = {
    "jacobian_convention": "J[component,axis]=partial_axis velocity_component",
    "self_advection": "(u.grad)u_i=sum_j u_j J[i,j]",
    "source_coordinate_scope": "inner PA.10 contraction-center domain only",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "inner_cartesian_center_velocity_inherited_from_agent1_803": True,
    "inner_cartesian_center_velocity_dt_inherited_from_agent1_811": True,
    "inner_cartesian_center_spatial_jacobian_inherited_from_agent1_819": True,
    "inner_cartesian_center_self_advection_executable": True,
    "inner_cartesian_center_self_advection_vectorized": True,
    "self_advection_uses_finite_difference_in_production": False,
    "source_center_is_final_corrected_fixed_point": False,
    "fixed_point_correction_materialized": False,
    "global_cartesian_spacetime_leading_velocity_materialized": False,
    "outer_join_localization_materialized": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "viscous_laplacian_materialized_here": False,
    "complete_candidate_api_ready": False,
    "complete_kokuno_composite_velocity": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _relative_max(actual: Any, reference: Any) -> float:
    actual_arr = np.asarray(actual, dtype=float)
    reference_arr = np.asarray(reference, dtype=float)
    scale = np.maximum(np.maximum(np.abs(actual_arr), np.abs(reference_arr)), 1.0)
    return float(np.max(np.abs(actual_arr - reference_arr) / scale))


def _directional_fd4(
    velocity_callable,
    point: np.ndarray,
    t: float,
    direction: np.ndarray,
    step: float,
) -> np.ndarray:
    """Verifier-only centered FD4 directional derivative at fixed time."""
    p = np.asarray(point, dtype=float)
    d = np.asarray(direction, dtype=float)
    h = float(step)
    fp2 = velocity_callable(*(p + 2.0 * h * d), t)
    fp1 = velocity_callable(*(p + h * d), t)
    fm1 = velocity_callable(*(p - h * d), t)
    fm2 = velocity_callable(*(p - 2.0 * h * d), t)
    return (-fp2 + 8.0 * fp1 - 8.0 * fm1 + fm2) / (12.0 * h)


@dataclass(frozen=True)
class KokunoPA10CartesianCenterSelfAdvection:
    """Expose analytic ``(u . grad)u`` for the current inner Cartesian center."""

    spatial: KokunoPA10CartesianCenterSpatialDerivatives = field(
        default_factory=KokunoPA10CartesianCenterSpatialDerivatives
    )

    def __post_init__(self) -> None:
        if not isinstance(self.spatial, KokunoPA10CartesianCenterSpatialDerivatives):
            raise TypeError("spatial must be KokunoPA10CartesianCenterSpatialDerivatives")

    @property
    def field(self):
        return self.spatial.field

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.spatial.velocity(x, y, z, t)

    def velocity_dt(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.spatial.velocity_dt(x, y, z, t)

    def velocity_jacobian(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.spatial.velocity_jacobian(x, y, z, t)

    def self_advection(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Return the analytic convective term ``(u . grad)u``.

        ``velocity_jacobian`` uses ``J[component,axis]``.  Contracting the
        final axis with the Cartesian velocity therefore gives exactly
        ``sum_axis u_axis * partial_axis u_component``.
        """
        velocity = self.velocity(x, y, z, t)
        jacobian = self.velocity_jacobian(x, y, z, t)
        out = np.einsum("...ij,...j->...i", jacobian, velocity)
        if np.any(~np.isfinite(out)):
            raise RuntimeError("analytic PA.10 center self-advection became nonfinite")
        return out

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "spatial_derivative_configuration": self.spatial.configuration(),
            "self_advection_realization": "analytic-jacobian-contraction-v1",
        }

    @property
    def field_sha256(self) -> str:
        return self.spatial.field_sha256

    @property
    def temporal_derivative_sha256(self) -> str:
        return self.spatial.temporal_derivative_sha256

    @property
    def spatial_derivative_sha256(self) -> str:
        return self.spatial.spatial_derivative_sha256

    @property
    def self_advection_sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.configuration()).encode("utf-8")).hexdigest()

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA10CartesianCenterSelfAdvection":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        if payload.get("schema") != SCHEMA:
            raise ValueError("configuration schema mismatch")
        if payload.get("self_advection_realization") != "analytic-jacobian-contraction-v1":
            raise ValueError("unsupported self-advection realization")
        return cls(
            spatial=KokunoPA10CartesianCenterSpatialDerivatives.from_configuration(
                payload["spatial_derivative_configuration"]
            )
        )

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoPA10CartesianCenterSelfAdvection":
        return cls.from_configuration(json.loads(Path(path).read_text()))


def build_report() -> dict[str, Any]:
    """Build a deterministic Agent-1 engineering receipt for this increment."""
    candidate = KokunoPA10CartesianCenterSelfAdvection()
    x0, x1 = candidate.spatial.source_X_interval
    span = x1 - x0
    if not span > 0.0:
        raise RuntimeError("source inner X interval is degenerate")

    X = x0 + span * np.array([0.18, 0.24, 0.30, 0.36, 0.21, 0.33])
    eta = np.array([-0.36, -0.22, -0.06, 0.11, 0.27, 0.39])
    t = np.array([0.44, 0.47, 0.50, 0.53, 0.56, 0.49])
    theta = np.array([0.17, 0.83, 1.41, 2.03, 2.67, 3.31])
    physical = candidate.field.cartesian_from_similarity(X, eta, t, theta)
    points = np.stack((physical["x"], physical["y"], physical["z"]), axis=-1)

    velocity = candidate.velocity(points[:, 0], points[:, 1], points[:, 2], t)
    jacobian = candidate.velocity_jacobian(
        points[:, 0], points[:, 1], points[:, 2], t
    )
    analytic = candidate.self_advection(points[:, 0], points[:, 1], points[:, 2], t)
    explicit = np.einsum("...ij,...j->...i", jacobian, velocity)

    coarse = []
    fine = []
    for point, time, vel in zip(points, t, velocity):
        speed = float(np.linalg.norm(vel))
        if not speed > 0.0:
            raise RuntimeError("self-advection verifier encountered zero velocity")
        direction = vel / speed
        coarse.append(
            speed * _directional_fd4(candidate.velocity, point, float(time), direction, 4.0e-4)
        )
        fine.append(
            speed * _directional_fd4(candidate.velocity, point, float(time), direction, 2.0e-4)
        )
    coarse_arr = np.asarray(coarse)
    fine_arr = np.asarray(fine)

    # Rotation covariance is checked independently on one off-axis source point.
    phi = 0.731
    c, s = np.cos(phi), np.sin(phi)
    p0 = points[2]
    p1 = np.array([c * p0[0] - s * p0[1], s * p0[0] + c * p0[1], p0[2]])
    a0 = candidate.self_advection(*p0, float(t[2]))
    a1 = candidate.self_advection(*p1, float(t[2]))
    a0_rot = np.array([c * a0[0] - s * a0[1], s * a0[0] + c * a0[1], a0[2]])

    axis = candidate.field.cartesian_from_similarity(
        np.zeros(3), np.array([-0.2, 0.0, 0.2]), np.full(3, 0.5), np.zeros(3)
    )
    axis_adv = candidate.self_advection(axis["x"], axis["y"], axis["z"], axis["t"])

    checks = {
        "production_matches_explicit_jacobian_contraction_abs_max": float(
            np.max(np.abs(analytic - explicit))
        ),
        "self_advection_vs_directional_fd4_fine_relative_max": _relative_max(
            analytic, fine_arr
        ),
        "directional_fd4_fine_vs_coarse_relative_max": _relative_max(
            fine_arr, coarse_arr
        ),
        "rotation_covariance_relative_max": _relative_max(a1, a0_rot),
        "axis_transverse_self_advection_abs_max": float(np.max(np.abs(axis_adv[:, :2]))),
        "self_advection_nontrivial_rms": float(np.sqrt(np.mean(analytic * analytic))),
        "all_probe_values_finite": bool(
            np.all(np.isfinite(velocity))
            and np.all(np.isfinite(jacobian))
            and np.all(np.isfinite(analytic))
            and np.all(np.isfinite(axis_adv))
        ),
    }

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "source": {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "corrected_release": CORRECTED_RELEASE,
            "corrected_release_date": CORRECTED_RELEASE_DATE,
        },
        "source_formulas": _SOURCE_FORMULAS,
        "configuration": candidate.configuration(),
        "field_sha256": candidate.field_sha256,
        "temporal_derivative_sha256": candidate.temporal_derivative_sha256,
        "spatial_derivative_sha256": candidate.spatial_derivative_sha256,
        "self_advection_sha256": candidate.self_advection_sha256,
        "machine_checks": checks,
        "truth_boundary": dict(_TRUTH_BOUNDARY),
        "scientific_gates": {
            "momentum_max_l2": 1.0e-3,
            "divergence_max_l2": 1.0e-5,
            "free_residual_defined_forcing_forbidden": True,
        },
    }
    payload["receipt_sha256"] = hashlib.sha256(
        _canonical_json(payload).encode("utf-8")
    ).hexdigest()
    return payload


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config-output", type=Path, required=True)
    args = parser.parse_args()

    candidate = KokunoPA10CartesianCenterSelfAdvection()
    report = build_report()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    candidate.save_configuration(args.config_output)
    print(json.dumps(report["machine_checks"], indent=2, sort_keys=True))


if __name__ == "__main__":
    _main()
