"""Controlled Cartesian Laplacian for the executable PA.10 center velocity.

Pinned public provenance:
KokunoYumeto/yang-mills-interacting-workbench@
143f6773feb424ad9ed3a8d116653200f20346b7,
navier-stokes/navier_stokes_workbench.tex, corrected 2026-09-09 reader.

Agent-1 #829 exposes the executable inner PA.10 contraction-center velocity,
its analytic time derivative, analytic first spatial derivatives, and analytic
self-advection.  The next candidate-facing viscous ingredient is

    Delta u = partial_xx u + partial_yy u + partial_zz u.

The corrected public reconstruction fixes the velocity/profile formulas, but
this repository does not yet expose the full second analytic profile chain
needed for a compact closed-form Cartesian Laplacian.  This module therefore
uses one frozen centered seven-point FD6 realization of the *public velocity*
as an explicitly repository-numerical derivative surface.  The production
step is fixed below and is not caller-tunable.  A different FD4 realization is
used only for engineering cross-checks.

This remains the inner PA.10 contraction center.  It does not invent an outer
continuation, pressure, forcing, corrected fixed point, complete momentum
residual, or independent PDE validation.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa10_cartesian_center_self_advection import (
    KokunoPA10CartesianCenterSelfAdvection,
)
from .kokuno_pa10_cartesian_center_velocity import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)


SCHEMA = "kokuno-pa10-cartesian-center-laplacian-v1"
PRODUCTION_SPATIAL_STEP = 1.0e-3

_SOURCE_FORMULAS = {
    "source_velocity": (
        "u1=(v0/(2q))x-q^(-A-1/2)Fy; "
        "u2=(v0/(2q))y+q^(-A-1/2)Fx; u3=q^(-A)U"
    ),
    "laplacian_definition": "Delta u=sum_axis partial_axis_axis u",
    "production_realization": (
        "repository numerical realization: centered seven-point FD6 pure "
        "second derivatives of the public inner velocity"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "inner_cartesian_center_velocity_inherited_from_agent1_803": True,
    "inner_cartesian_center_velocity_dt_inherited_from_agent1_811": True,
    "inner_cartesian_center_spatial_jacobian_inherited_from_agent1_819": True,
    "inner_cartesian_center_self_advection_inherited_from_agent1_829": True,
    "inner_cartesian_center_velocity_laplacian_executable": True,
    "inner_cartesian_center_velocity_laplacian_vectorized": True,
    "laplacian_uses_finite_difference_in_production": True,
    "laplacian_is_public_source_formula": False,
    "production_laplacian_step_caller_tunable": False,
    "outside_inner_domain_zero_extension_used": False,
    "source_center_is_final_corrected_fixed_point": False,
    "fixed_point_correction_materialized": False,
    "global_cartesian_spacetime_leading_velocity_materialized": False,
    "outer_join_localization_materialized": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
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


def _broadcast_xyz_t(
    x: Any, y: Any, z: Any, t: Any
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    arrays = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    if any(np.any(~np.isfinite(a)) for a in arrays):
        raise ValueError("x,y,z,t must contain only finite values")
    return arrays[0], arrays[1], arrays[2], arrays[3]


def _shifted_velocity(
    velocity_callable,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    axis: int,
    offset: float,
) -> np.ndarray:
    coords = [x, y, z]
    shifted = list(coords)
    shifted[axis] = shifted[axis] + float(offset)
    return np.asarray(velocity_callable(shifted[0], shifted[1], shifted[2], t), dtype=float)


def _fd6_laplacian(
    velocity_callable,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    step: float,
) -> np.ndarray:
    """Centered seven-point sixth-order Laplacian of a vector callable."""
    h = float(step)
    if not np.isfinite(h) or h <= 0.0:
        raise ValueError("step must be finite and positive")
    x_arr, y_arr, z_arr, t_arr = _broadcast_xyz_t(x, y, z, t)
    u0 = np.asarray(velocity_callable(x_arr, y_arr, z_arr, t_arr), dtype=float)
    lap = np.zeros_like(u0)
    for axis in range(3):
        fm3 = _shifted_velocity(velocity_callable, x_arr, y_arr, z_arr, t_arr, axis, -3.0 * h)
        fm2 = _shifted_velocity(velocity_callable, x_arr, y_arr, z_arr, t_arr, axis, -2.0 * h)
        fm1 = _shifted_velocity(velocity_callable, x_arr, y_arr, z_arr, t_arr, axis, -h)
        fp1 = _shifted_velocity(velocity_callable, x_arr, y_arr, z_arr, t_arr, axis, h)
        fp2 = _shifted_velocity(velocity_callable, x_arr, y_arr, z_arr, t_arr, axis, 2.0 * h)
        fp3 = _shifted_velocity(velocity_callable, x_arr, y_arr, z_arr, t_arr, axis, 3.0 * h)
        lap += (
            2.0 * fm3
            - 27.0 * fm2
            + 270.0 * fm1
            - 490.0 * u0
            + 270.0 * fp1
            - 27.0 * fp2
            + 2.0 * fp3
        ) / (180.0 * h * h)
    return lap


def _fd4_laplacian(
    velocity_callable,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    step: float,
) -> np.ndarray:
    """Verifier-only centered five-point fourth-order Laplacian."""
    h = float(step)
    if not np.isfinite(h) or h <= 0.0:
        raise ValueError("step must be finite and positive")
    x_arr, y_arr, z_arr, t_arr = _broadcast_xyz_t(x, y, z, t)
    u0 = np.asarray(velocity_callable(x_arr, y_arr, z_arr, t_arr), dtype=float)
    lap = np.zeros_like(u0)
    for axis in range(3):
        fm2 = _shifted_velocity(velocity_callable, x_arr, y_arr, z_arr, t_arr, axis, -2.0 * h)
        fm1 = _shifted_velocity(velocity_callable, x_arr, y_arr, z_arr, t_arr, axis, -h)
        fp1 = _shifted_velocity(velocity_callable, x_arr, y_arr, z_arr, t_arr, axis, h)
        fp2 = _shifted_velocity(velocity_callable, x_arr, y_arr, z_arr, t_arr, axis, 2.0 * h)
        lap += (-fp2 + 16.0 * fp1 - 30.0 * u0 + 16.0 * fm1 - fm2) / (
            12.0 * h * h
        )
    return lap


def _relative_max(actual: Any, reference: Any) -> float:
    a = np.asarray(actual, dtype=float)
    b = np.asarray(reference, dtype=float)
    scale = np.maximum(np.maximum(np.abs(a), np.abs(b)), 1.0)
    return float(np.max(np.abs(a - b) / scale))


@dataclass(frozen=True)
class KokunoPA10CartesianCenterLaplacian:
    """Expose a frozen numerical ``Delta u`` for the current inner center."""

    nonlinear: KokunoPA10CartesianCenterSelfAdvection = field(
        default_factory=KokunoPA10CartesianCenterSelfAdvection
    )

    def __post_init__(self) -> None:
        if not isinstance(self.nonlinear, KokunoPA10CartesianCenterSelfAdvection):
            raise TypeError("nonlinear must be KokunoPA10CartesianCenterSelfAdvection")

    @property
    def field(self):
        return self.nonlinear.field

    @property
    def spatial(self):
        return self.nonlinear.spatial

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.nonlinear.velocity(x, y, z, t)

    def velocity_dt(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.nonlinear.velocity_dt(x, y, z, t)

    def velocity_jacobian(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.nonlinear.velocity_jacobian(x, y, z, t)

    def self_advection(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.nonlinear.self_advection(x, y, z, t)

    def velocity_laplacian(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Return frozen-FD6 ``Delta u``; any stencil domain escape fails closed."""
        out = _fd6_laplacian(
            self.velocity,
            x,
            y,
            z,
            t,
            step=PRODUCTION_SPATIAL_STEP,
        )
        if np.any(~np.isfinite(out)):
            raise RuntimeError("PA.10 center velocity Laplacian became nonfinite")
        return out

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "self_advection_configuration": self.nonlinear.configuration(),
            "laplacian_realization": "centered-seven-point-fd6-public-velocity-v1",
            "production_spatial_step": PRODUCTION_SPATIAL_STEP,
        }

    @property
    def field_sha256(self) -> str:
        return self.nonlinear.field_sha256

    @property
    def temporal_derivative_sha256(self) -> str:
        return self.nonlinear.temporal_derivative_sha256

    @property
    def spatial_derivative_sha256(self) -> str:
        return self.nonlinear.spatial_derivative_sha256

    @property
    def self_advection_sha256(self) -> str:
        return self.nonlinear.self_advection_sha256

    @property
    def laplacian_sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.configuration()).encode("utf-8")).hexdigest()

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA10CartesianCenterLaplacian":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        if payload.get("schema") != SCHEMA:
            raise ValueError("configuration schema mismatch")
        if payload.get("laplacian_realization") != "centered-seven-point-fd6-public-velocity-v1":
            raise ValueError("unsupported Laplacian realization")
        if payload.get("production_spatial_step") != PRODUCTION_SPATIAL_STEP:
            raise ValueError("production Laplacian step drifted from the frozen value")
        return cls(
            nonlinear=KokunoPA10CartesianCenterSelfAdvection.from_configuration(
                payload["self_advection_configuration"]
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
    ) -> "KokunoPA10CartesianCenterLaplacian":
        return cls.from_configuration(json.loads(Path(path).read_text()))


def build_report() -> dict[str, Any]:
    """Build a deterministic engineering receipt without making a PDE claim."""
    candidate = KokunoPA10CartesianCenterLaplacian()
    x0, x1 = candidate.spatial.source_X_interval
    span = x1 - x0
    if not span > 0.0:
        raise RuntimeError("source inner X interval is degenerate")

    X = x0 + span * np.array([0.16, 0.20, 0.24, 0.28, 0.18, 0.26])
    eta = np.array([-0.30, -0.18, -0.05, 0.09, 0.21, 0.30])
    t = np.array([0.46, 0.48, 0.50, 0.52, 0.54, 0.49])
    theta = np.array([0.23, 0.79, 1.37, 1.97, 2.59, 3.17])
    physical = candidate.field.cartesian_from_similarity(X, eta, t, theta)

    production = candidate.velocity_laplacian(
        physical["x"], physical["y"], physical["z"], physical["t"]
    )
    fd4_coarse = _fd4_laplacian(
        candidate.velocity,
        physical["x"], physical["y"], physical["z"], physical["t"],
        step=2.0e-3,
    )
    fd4_fine = _fd4_laplacian(
        candidate.velocity,
        physical["x"], physical["y"], physical["z"], physical["t"],
        step=1.0e-3,
    )

    phi = 0.713
    c, s = np.cos(phi), np.sin(phi)
    p0 = np.array([physical["x"][2], physical["y"][2], physical["z"][2]])
    p1 = np.array([c * p0[0] - s * p0[1], s * p0[0] + c * p0[1], p0[2]])
    l0 = candidate.velocity_laplacian(*p0, float(t[2]))
    l1 = candidate.velocity_laplacian(*p1, float(t[2]))
    l0_rot = np.array([c * l0[0] - s * l0[1], s * l0[0] + c * l0[1], l0[2]])

    boundary = candidate.field.cartesian_from_similarity(x1, 0.0, 0.5, 0.0)
    boundary_fail_closed = False
    try:
        candidate.velocity_laplacian(
            boundary["x"], boundary["y"], boundary["z"], boundary["t"]
        )
    except ValueError:
        boundary_fail_closed = True

    checks = {
        "production_step": PRODUCTION_SPATIAL_STEP,
        "laplacian_vs_independent_fd4_fine_relative_max": _relative_max(
            production, fd4_fine
        ),
        "independent_fd4_fine_vs_coarse_relative_max": _relative_max(
            fd4_fine, fd4_coarse
        ),
        "rotation_covariance_relative_max": _relative_max(l1, l0_rot),
        "laplacian_nontrivial_rms": float(np.sqrt(np.mean(production * production))),
        "inner_boundary_stencil_fails_closed": boundary_fail_closed,
        "all_probe_values_finite": bool(
            np.all(np.isfinite(production))
            and np.all(np.isfinite(fd4_coarse))
            and np.all(np.isfinite(fd4_fine))
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
        "laplacian_sha256": candidate.laplacian_sha256,
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

    candidate = KokunoPA10CartesianCenterLaplacian()
    report = build_report()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    candidate.save_configuration(args.config_output)
    print(json.dumps(report["machine_checks"], indent=2, sort_keys=True))


if __name__ == "__main__":
    _main()
