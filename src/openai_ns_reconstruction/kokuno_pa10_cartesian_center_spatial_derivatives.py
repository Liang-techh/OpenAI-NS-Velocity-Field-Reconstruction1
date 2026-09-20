"""Analytic spatial derivatives for the executable PA.10 Cartesian center velocity.

Pinned public provenance:
KokunoYumeto/yang-mills-interacting-workbench@
143f6773feb424ad9ed3a8d116653200f20346b7,
navier-stokes/navier_stokes_workbench.tex, corrected 2026-09-09 reader.

Agent-1 #803 materializes the public inner contraction-center velocity

    tau=1-t,
    q-z^2 q^(2h)-tau=0,
    X=(x^2+y^2)/(2q),       eta=z/q^D,

    u1=(v0/(2q))x-q^(-A-1/2)Fy,
    u2=(v0/(2q))y+q^(-A-1/2)Fx,
    u3=q^(-A)U.

Agent-1 #811 supplies the analytic ``(v0)_eta`` required below.  At fixed
time the public coordinate identities imply

    q_x=q_y=0,
    q_z=2 z q^(2h)/L,
    X_x=x/q,  X_y=y/q,  X_z=-X q_z/q,
    eta_x=eta_y=0,
    eta_z=(1-eta^2) q^(-D)/L,

with L=1-2h eta^2.  This module differentiates the displayed Cartesian
velocity analytically, exposing the candidate-facing spatial Jacobian,
divergence and vorticity without finite differences in production.

Finite differences appear only in focused/report cross-checks.  This remains
the *inner PA.10 contraction center*.  It does not materialize the corrected
fixed point, an outer/global join, matched pressure, restricted forcing, or
an independent Navier--Stokes validation.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa10_cartesian_center_velocity import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
    KokunoPA10CartesianCenterVelocity,
)
from .kokuno_pa10_cartesian_center_velocity_dt import (
    KokunoPA10CartesianCenterVelocityTimeDerivative,
)


SCHEMA = "kokuno-pa10-cartesian-center-spatial-derivatives-v1"

_SOURCE_FORMULAS = {
    "coordinates": (
        "tau=1-t; q-z^2 q^(2h)-tau=0; X=(x^2+y^2)/(2q); "
        "eta=z/q^D; L=1-2h eta^2"
    ),
    "coordinate_spatial_derivatives": (
        "q_x=q_y=0; q_z=2 z q^(2h)/L; X_x=x/q; X_y=y/q; "
        "X_z=-X q_z/q; eta_x=eta_y=0; eta_z=(1-eta^2) q^(-D)/L"
    ),
    "cartesian_velocity": (
        "u1=(v0/(2q))x-q^(-A-1/2)Fy; "
        "u2=(v0/(2q))y+q^(-A-1/2)Fx; u3=q^(-A)U"
    ),
    "spatial_jacobian": (
        "differentiate the displayed Cartesian velocity using analytic "
        "F_X,F_eta,U_X,U_eta,(v0)_X,(v0)_eta and the coordinate derivatives"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "inner_cartesian_center_velocity_inherited_from_agent1_803": True,
    "inner_cartesian_center_velocity_dt_inherited_from_agent1_811": True,
    "source_coordinate_spatial_derivatives_analytic": True,
    "inner_cartesian_center_velocity_spatial_jacobian_executable": True,
    "inner_cartesian_center_velocity_spatial_jacobian_vectorized": True,
    "inner_cartesian_center_divergence_executable": True,
    "inner_cartesian_center_vorticity_executable": True,
    "spatial_derivatives_use_finite_difference_in_production": False,
    "source_center_is_final_corrected_fixed_point": False,
    "fixed_point_correction_materialized": False,
    "global_cartesian_spacetime_leading_velocity_materialized": False,
    "outer_join_localization_materialized": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "complete_candidate_api_ready": False,
    "complete_kokuno_composite_velocity": False,
    "unified_cartesian_velocity_export_ready": False,
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


@dataclass(frozen=True)
class KokunoPA10CartesianCenterSpatialDerivatives:
    """Add analytic spatial derivatives to the current inner Cartesian center."""

    derivative: KokunoPA10CartesianCenterVelocityTimeDerivative = field(
        default_factory=KokunoPA10CartesianCenterVelocityTimeDerivative
    )

    def __post_init__(self) -> None:
        if not isinstance(
            self.derivative, KokunoPA10CartesianCenterVelocityTimeDerivative
        ):
            raise TypeError(
                "derivative must be KokunoPA10CartesianCenterVelocityTimeDerivative"
            )

    @property
    def field(self) -> KokunoPA10CartesianCenterVelocity:
        return self.derivative.field

    @property
    def physical_profiles(self):
        return self.derivative.physical_profiles

    @property
    def source_X_interval(self) -> tuple[float, float]:
        return self.derivative.source_X_interval

    @property
    def A(self) -> float:
        return float(self.derivative.A)

    @property
    def D(self) -> float:
        return float(self.derivative.D)

    @property
    def h(self) -> float:
        return float(self.derivative.h)

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.field.velocity(x, y, z, t)

    def velocity_dt(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.derivative.velocity_dt(x, y, z, t)

    def coordinate_spatial_derivatives(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        """Return analytic spatial derivatives of the public similarity map."""
        coords = self.field.similarity_coordinates(x, y, z, t)
        q = coords["q"]
        eta = coords["eta"]
        L = 1.0 - 2.0 * self.h * eta * eta
        if np.any(L <= 0.0):
            raise RuntimeError("source coordinate L lost positivity")

        q_z = 2.0 * coords["z"] * np.power(q, 2.0 * self.h) / L
        X_x = coords["x"] / q
        X_y = coords["y"] / q
        X_z = -coords["X"] * q_z / q
        eta_z = (1.0 - eta * eta) * np.power(q, -self.D) / L
        zeros = np.zeros_like(q)
        return {
            **coords,
            "L": L,
            "q_x": zeros,
            "q_y": zeros,
            "q_z": q_z,
            "X_x": X_x,
            "X_y": X_y,
            "X_z": X_z,
            "eta_x": zeros,
            "eta_y": zeros,
            "eta_z": eta_z,
        }

    def values(self, x: Any, y: Any, z: Any, t: Any) -> dict[str, np.ndarray]:
        """Return the velocity plus analytic coefficient/spatial derivatives."""
        base = self.field.values(x, y, z, t)
        coord_d = self.coordinate_spatial_derivatives(x, y, z, t)
        X = base["X"]
        eta = base["eta"]
        x0, x1 = self.source_X_interval
        radial_tol = 64.0 * np.finfo(float).eps * max(1.0, abs(x1))
        if np.any(X < x0 - radial_tol) or np.any(X > x1 + radial_tol):
            raise ValueError("point lies outside the source inner X interval")
        X_eval = np.clip(X, x0, x1)

        profile_d = self.physical_profiles.derivatives(X_eval, eta)
        v0_eta = self.derivative.v0_eta(X_eval, eta)
        q = base["q"]
        q_z = coord_d["q_z"]

        v0_x = profile_d["v_0_X"] * coord_d["X_x"]
        v0_y = profile_d["v_0_X"] * coord_d["X_y"]
        v0_z = profile_d["v_0_X"] * coord_d["X_z"] + v0_eta * coord_d["eta_z"]
        F_x = profile_d["F_0_X"] * coord_d["X_x"]
        F_y = profile_d["F_0_X"] * coord_d["X_y"]
        F_z = (
            profile_d["F_0_X"] * coord_d["X_z"]
            + profile_d["F_0_eta"] * coord_d["eta_z"]
        )
        U_x = profile_d["U_0_X"] * coord_d["X_x"]
        U_y = profile_d["U_0_X"] * coord_d["X_y"]
        U_z = (
            profile_d["U_0_X"] * coord_d["X_z"]
            + profile_d["U_0_eta"] * coord_d["eta_z"]
        )

        radial = base["v_0"] / (2.0 * q)
        radial_x = v0_x / (2.0 * q)
        radial_y = v0_y / (2.0 * q)
        radial_z = v0_z / (2.0 * q) - base["v_0"] * q_z / (2.0 * q * q)

        swirl_scale = np.power(q, -self.A - 0.5)
        swirl = swirl_scale * base["F_0"]
        swirl_x = swirl_scale * F_x
        swirl_y = swirl_scale * F_y
        swirl_z = swirl_scale * (
            F_z + (-self.A - 0.5) * (q_z / q) * base["F_0"]
        )

        axial_scale = np.power(q, -self.A)
        axial_x = axial_scale * U_x
        axial_y = axial_scale * U_y
        axial_z = axial_scale * (U_z - self.A * (q_z / q) * base["U_0"])

        x_arr = base["x"]
        y_arr = base["y"]
        u_x = radial + x_arr * radial_x - y_arr * swirl_x
        u_y = x_arr * radial_y - swirl - y_arr * swirl_y
        u_z = x_arr * radial_z - y_arr * swirl_z
        v_x = y_arr * radial_x + swirl + x_arr * swirl_x
        v_y = radial + y_arr * radial_y + x_arr * swirl_y
        v_z = y_arr * radial_z + x_arr * swirl_z
        w_x = axial_x
        w_y = axial_y
        w_z = axial_z

        jacobian = np.stack(
            (
                np.stack((u_x, u_y, u_z), axis=-1),
                np.stack((v_x, v_y, v_z), axis=-1),
                np.stack((w_x, w_y, w_z), axis=-1),
            ),
            axis=-2,
        )
        divergence = u_x + v_y + w_z
        vorticity = np.stack((w_y - v_z, u_z - w_x, v_x - u_y), axis=-1)

        if np.any(~np.isfinite(jacobian)):
            raise RuntimeError("analytic spatial velocity Jacobian became nonfinite")
        if np.any(~np.isfinite(vorticity)):
            raise RuntimeError("analytic vorticity became nonfinite")

        return {
            **base,
            "L_coordinate": coord_d["L"],
            "q_z": q_z,
            "X_x": coord_d["X_x"],
            "X_y": coord_d["X_y"],
            "X_z": coord_d["X_z"],
            "eta_z": coord_d["eta_z"],
            "v_0_eta": v0_eta,
            "v_0_x": v0_x,
            "v_0_y": v0_y,
            "v_0_z": v0_z,
            "F_0_x": F_x,
            "F_0_y": F_y,
            "F_0_z": F_z,
            "U_0_x": U_x,
            "U_0_y": U_y,
            "U_0_z": U_z,
            "radial_coefficient": radial,
            "radial_coefficient_x": radial_x,
            "radial_coefficient_y": radial_y,
            "radial_coefficient_z": radial_z,
            "swirl_coefficient": swirl,
            "swirl_coefficient_x": swirl_x,
            "swirl_coefficient_y": swirl_y,
            "swirl_coefficient_z": swirl_z,
            "axial_coefficient_x": axial_x,
            "axial_coefficient_y": axial_y,
            "axial_coefficient_z": axial_z,
            "velocity_jacobian": jacobian,
            "divergence": divergence,
            "vorticity": vorticity,
        }

    def velocity_jacobian(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Return J[component,axis]=partial_axis velocity_component."""
        return self.values(x, y, z, t)["velocity_jacobian"]

    def divergence(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.values(x, y, z, t)["divergence"]

    def vorticity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.values(x, y, z, t)["vorticity"]

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "temporal_derivative_configuration": self.derivative.configuration(),
            "spatial_derivative_realization": "analytic-source-chain-rule-v1",
        }

    @property
    def field_sha256(self) -> str:
        return self.field.field_sha256

    @property
    def temporal_derivative_sha256(self) -> str:
        return self.derivative.derivative_sha256

    @property
    def spatial_derivative_sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.configuration()).encode("utf-8")).hexdigest()

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA10CartesianCenterSpatialDerivatives":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        if payload.get("schema") != SCHEMA:
            raise ValueError("configuration schema mismatch")
        if payload.get("spatial_derivative_realization") != "analytic-source-chain-rule-v1":
            raise ValueError("unsupported spatial derivative realization")
        return cls(
            derivative=KokunoPA10CartesianCenterVelocityTimeDerivative.from_configuration(
                payload["temporal_derivative_configuration"]
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
    ) -> "KokunoPA10CartesianCenterSpatialDerivatives":
        return cls.from_configuration(json.loads(Path(path).read_text()))

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    @staticmethod
    def _fd4_axis(
        field: KokunoPA10CartesianCenterVelocity,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        t: np.ndarray,
        axis: int,
        step: float,
    ) -> np.ndarray:
        coords = [np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float)]

        def shifted(multiplier: float) -> np.ndarray:
            shifted_coords = [value.copy() for value in coords]
            shifted_coords[axis] = shifted_coords[axis] + multiplier * step
            return field.velocity(shifted_coords[0], shifted_coords[1], shifted_coords[2], t)

        return (
            shifted(-2.0)
            - 8.0 * shifted(-1.0)
            + 8.0 * shifted(1.0)
            - shifted(2.0)
        ) / (12.0 * step)

    @classmethod
    def _fd4_jacobian(
        cls,
        field: KokunoPA10CartesianCenterVelocity,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        t: np.ndarray,
        step: float,
    ) -> np.ndarray:
        columns = [cls._fd4_axis(field, x, y, z, t, axis=axis, step=step) for axis in range(3)]
        return np.stack(columns, axis=-1)

    def report(self) -> dict[str, Any]:
        _, xmax = self.source_X_interval
        X = np.asarray([0.06, 0.14, 0.23, 0.31], dtype=float) * xmax
        eta = np.asarray([0.08, -0.17, 0.27, -0.36], dtype=float)
        t = np.asarray([0.34, 0.46, 0.54, 0.66], dtype=float)
        theta = np.asarray([0.25, 0.85, 1.45, 2.05], dtype=float)
        physical = self.field.cartesian_from_similarity(X, eta, t, theta)

        analytic = self.velocity_jacobian(physical["x"], physical["y"], physical["z"], physical["t"])
        fd_coarse = self._fd4_jacobian(self.field, physical["x"], physical["y"], physical["z"], physical["t"], 2.0e-4)
        fd_fine = self._fd4_jacobian(self.field, physical["x"], physical["y"], physical["z"], physical["t"], 1.0e-4)

        analytic_vorticity = self.vorticity(physical["x"], physical["y"], physical["z"], physical["t"])
        fd_vorticity = np.stack(
            (
                fd_fine[..., 2, 1] - fd_fine[..., 1, 2],
                fd_fine[..., 0, 2] - fd_fine[..., 2, 0],
                fd_fine[..., 1, 0] - fd_fine[..., 0, 1],
            ),
            axis=-1,
        )
        analytic_divergence = np.trace(analytic, axis1=-2, axis2=-1)
        jac_scale = max(1.0, float(np.max(np.abs(analytic))))
        divergence_relative = float(np.max(np.abs(analytic_divergence)) / jac_scale)

        coord = self.coordinate_spatial_derivatives(physical["x"], physical["y"], physical["z"], physical["t"])
        coord_step = 1.0e-5
        cplus = self.field.similarity_coordinates(physical["x"], physical["y"], physical["z"] + coord_step, physical["t"])
        cminus = self.field.similarity_coordinates(physical["x"], physical["y"], physical["z"] - coord_step, physical["t"])
        q_z_fd = (cplus["q"] - cminus["q"]) / (2.0 * coord_step)
        X_z_fd = (cplus["X"] - cminus["X"]) / (2.0 * coord_step)
        eta_z_fd = (cplus["eta"] - cminus["eta"]) / (2.0 * coord_step)

        Xr = 0.20 * xmax
        etar = 0.19
        tr = 0.51
        angle = 0.63
        p0 = self.field.cartesian_from_similarity(Xr, etar, tr, 0.0)
        p1 = self.field.cartesian_from_similarity(Xr, etar, tr, angle)
        j0 = self.velocity_jacobian(p0["x"], p0["y"], p0["z"], p0["t"])
        j1 = self.velocity_jacobian(p1["x"], p1["y"], p1["z"], p1["t"])
        c = math.cos(angle)
        s = math.sin(angle)
        rotation = np.asarray([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
        rotated = rotation @ j0 @ rotation.T

        axis = self.field.cartesian_from_similarity(
            np.zeros(3), np.asarray([-0.25, 0.0, 0.25]), np.asarray([0.36, 0.50, 0.64])
        )
        axis_jacobian = self.velocity_jacobian(axis["x"], axis["y"], axis["z"], axis["t"])

        payload: dict[str, Any] = {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "configuration": self.configuration(),
            "field_sha256": self.field_sha256,
            "temporal_derivative_sha256": self.temporal_derivative_sha256,
            "spatial_derivative_sha256": self.spatial_derivative_sha256,
            "machine_checks": {
                "velocity_jacobian_vs_fd4_fine_relative_max": _relative_max(analytic, fd_fine),
                "fd4_fine_vs_coarse_relative_max": _relative_max(fd_fine, fd_coarse),
                "vorticity_vs_fd4_curl_relative_max": _relative_max(analytic_vorticity, fd_vorticity),
                "analytic_divergence_sampled_abs_max": float(np.max(np.abs(analytic_divergence))),
                "analytic_divergence_relative_to_jacobian_scale": divergence_relative,
                "q_z_vs_centered_fd_relative_max": _relative_max(coord["q_z"], q_z_fd),
                "X_z_vs_centered_fd_relative_max": _relative_max(coord["X_z"], X_z_fd),
                "eta_z_vs_centered_fd_relative_max": _relative_max(coord["eta_z"], eta_z_fd),
                "rotation_covariance_relative_max": _relative_max(j1, rotated),
                "axis_jacobian_all_finite": bool(np.all(np.isfinite(axis_jacobian))),
                "vorticity_nontrivial_on_probe": bool(np.any(np.abs(analytic_vorticity) > 0.0)),
                "all_spatial_derivative_probe_values_finite": bool(
                    np.all(np.isfinite(analytic)) and np.all(np.isfinite(analytic_vorticity))
                ),
            },
            "scientific_gates": {
                "momentum_max_l2": 1.0e-3,
                "divergence_max_l2": 1.0e-5,
                "free_residual_defined_forcing_forbidden": True,
            },
            "truth_boundary": self.truth_boundary,
        }
        payload["receipt_sha256"] = hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()
        return payload

    def save_report(self, path: str | Path) -> dict[str, Any]:
        payload = self.report()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload


def _main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config-output", type=Path)
    args = parser.parse_args()
    spatial = KokunoPA10CartesianCenterSpatialDerivatives()
    payload = spatial.save_report(args.output)
    if args.config_output is not None:
        spatial.save_configuration(args.config_output)
    print("field_sha256=", payload["field_sha256"])
    print("temporal_derivative_sha256=", payload["temporal_derivative_sha256"])
    print("spatial_derivative_sha256=", payload["spatial_derivative_sha256"])
    print("receipt_sha256=", payload["receipt_sha256"])
    print("machine_checks=", payload["machine_checks"])


if __name__ == "__main__":
    _main()
